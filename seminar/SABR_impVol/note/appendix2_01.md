`````
.githooks/checks.py
"""Git がコミットまたは送信する内容を、差分に応じて検証する。"""

import json
import os

# Git が渡した引数と、リポジトリ内の検証コマンドをシェルなしで実行する。
import subprocess  # ruff: ignore[suspicious-subprocess-import]
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory
from types import MappingProxyType
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping

PYTHON_PROJECTS = (".githooks",)
CARGO_INPUTS = frozenset({"Cargo.toml", "Cargo.lock", "rust-toolchain.toml", "clippy.toml"})

type _Command = tuple[str | Path, ...]
type _FileKind = Literal["rust", "python", "toml", "json", "nix", "shell"]
type _FileCheck = tuple[str, _FileKind]
type _PythonCheck = tuple[str, tuple[str, ...]]

_SUFFIXES: Mapping[str, _FileKind] = MappingProxyType(
    {
        ".rs": "rust",
        ".py": "python",
        ".toml": "toml",
        ".json": "json",
        ".nix": "nix",
    }
)


@dataclass(frozen=True, slots=True)
class _CheckPlan:
    files: tuple[_FileCheck, ...]
    python: tuple[_PythonCheck, ...]
    crates: frozenset[str]
    workspace: bool

    def needs_push_checks(self) -> bool:
        """送信時にツリー全体が必要になるかを判定する。

        Returns:
            ツリー全体を使う検査が一つでもあれば真。

        """
        return bool(self.python or self.crates or self.workspace)


def _parse_changes(raw: bytes) -> dict[str, bool]:
    """削除や移動元を含む Git の差分を読む。

    Returns:
        パスから、変更後に通常ファイルとして読めるかどうかへの対応。

    """
    fields = raw.split(b"\0")
    return {
        os.fsdecode(path): header.split()[1] in {b"100644", b"100755"}
        for header, path in zip(fields[0:-1:2], fields[1::2], strict=True)
    }


def _file_kind(name: str) -> _FileKind | None:
    match name:
        case "flake.lock":
            return "json"
        case ".githooks/pre-commit" | ".githooks/pre-push":
            return "shell"
        case _:
            return _SUFFIXES.get(PurePosixPath(name).suffix)


def _crate_for_path(name: str) -> str | None:
    match PurePosixPath(name).parts:
        case ("crates", crate, _, *_):
            return crate
        case _:
            return None


def _python_files(directory: str, changed: Mapping[str, bool]) -> tuple[str, ...]:
    # 設定の変更は、変更していないファイルの検査にも影響する。
    if f"{directory}/pyproject.toml" in changed:
        return (".",)
    prefix = directory + "/"
    return tuple(
        path.removeprefix(prefix)
        for path, regular in changed.items()
        if regular and path.startswith(prefix) and path.endswith(".py")
    )


def _plan(changed: Mapping[str, bool]) -> _CheckPlan:
    return _CheckPlan(
        files=tuple(
            (path, kind)
            for path, regular in changed.items()
            if regular and (kind := _file_kind(path)) is not None
        ),
        python=tuple(
            (directory, files)
            for directory in PYTHON_PROJECTS
            if (files := _python_files(directory, changed))
        ),
        crates=frozenset(crate for path in changed if (crate := _crate_for_path(path)) is not None),
        workspace=bool(
            CARGO_INPUTS & changed.keys() or any(p.startswith(".cargo/") for p in changed)
        ),
    )


def _run(
    root: Path,
    command: _Command,
    *,
    env: Mapping[str, str] | None = None,
    quiet: bool = False,
) -> None:
    subprocess.run(  # ruff: ignore[subprocess-without-shell-equals-true]
        command,
        cwd=root,
        env=env,
        stdout=subprocess.DEVNULL if quiet else None,
        check=True,
    )


def _git(
    repo: Path,
    args: tuple[str, ...],
    *,
    env: Mapping[str, str] | None = None,
    data: bytes | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(  # ruff: ignore[subprocess-without-shell-equals-true]
        # 開発環境の PATH にある Git を使う。
        ("git", *args),  # ruff: ignore[start-process-with-partial-path]
        cwd=repo,
        env=env,
        input=data,
        stdout=subprocess.PIPE,
        check=check,
    )


def _read_toml(path: Path) -> Mapping[str, object]:
    with path.open("rb") as source:
        return tomllib.load(source)


def _check_config(root: Path, name: str, kind: Literal["toml", "json"]) -> None:
    try:
        if kind == "toml":
            _read_toml(root / name)
        else:
            json.loads((root / name).read_bytes())
    except ValueError as error:
        message = f"{name}: {error}"
        raise ValueError(message) from error


def _check_syntax(root: Path, name: str, kind: _FileKind) -> None:
    path = root / name
    match kind:
        case "python":
            compile(path.read_bytes(), name, "exec")
        case "toml" | "json":
            _check_config(root, name, kind)
        case "nix":
            _run(root, ("nix-instantiate", "--parse", path), quiet=True)
        case "shell":
            _run(root, ("sh", "-n", path))
        case "rust":
            pass  # Rust は全対象をまとめて検査する。


def _check_files(root: Path, files: tuple[_FileCheck, ...]) -> None:
    for name, kind in files:
        _check_syntax(root, name, kind)
    rust = tuple("./" + name for name, kind in files if kind == "rust")
    if rust:
        _run(
            root,
            ("rustfmt", "--check", "--edition", "2024", "--config", "skip_children=true", *rust),
        )


def _checkout(
    repo: Path, root: Path, files: tuple[_FileCheck, ...] | None, tip: str | None
) -> None:
    command = ("checkout-index", f"--prefix={root}/")
    env = None
    if tip is not None:
        # push の対象が別の参照でも、その commit の内容で検査する。
        env = {**os.environ, "GIT_INDEX_FILE": str(root.parent / "index")}
        _git(repo, ("read-tree", tip), env=env)
    if files is None:
        _git(repo, (*command, "--all"), env=env)
    else:
        # 構文チェックだけなら、対象のファイルだけ取り出す。
        paths = b"".join(os.fsencode(name) + b"\0" for name, _ in files)
        _git(repo, (*command, "-z", "--stdin"), env=env, data=paths)


def _python_commands(tools: object, files: tuple[str, ...]) -> Iterator[_Command]:
    if isinstance(tools, dict):
        if "ruff" in tools:
            yield ("ruff", "check", "--no-fix", "--", *files)
            yield ("ruff", "format", "--check", "--", *files)
        if "ty" in tools:
            yield ("ty", "check", "--", *files)


def _check_python(root: Path, project: _PythonCheck) -> None:
    directory, files = project
    config = root / directory / "pyproject.toml"
    if config.is_file():
        for command in _python_commands(_read_toml(config).get("tool"), files):
            _run(config.parent, command)


def _cargo_scope(root: Path, plan: _CheckPlan) -> tuple[str, ...]:
    if plan.workspace:
        return ("--workspace",)
    manifests = (root / "crates" / crate / "Cargo.toml" for crate in sorted(plan.crates))
    return tuple(arg for path in manifests if path.is_file() for arg in ("-p", _package_name(path)))


def _package_name(path: Path) -> str:
    package = _read_toml(path).get("package")
    name = package.get("name") if isinstance(package, dict) else None
    if not isinstance(name, str):
        message = f"{path}: package.name is missing or is not a string"
        raise TypeError(message)
    return name


def _cargo_commands(scope: tuple[str, ...]) -> tuple[_Command, ...]:
    flags = ("--locked", "--all-features", *scope)
    return (
        ("cargo", "clippy", *flags, "--all-targets", "--", "-D", "warnings"),
        ("cargo", "test", *flags),
        ("cargo", "doc", *flags, "--no-deps"),
    )


def _cargo_environment(target: Path, env: Mapping[str, str]) -> dict[str, str]:
    return {
        **env,
        "CARGO_TARGET_DIR": os.path.normpath(target),
        "RUSTDOCFLAGS": (env.get("RUSTDOCFLAGS", "") + " -D warnings").strip(),
    }


def _repository_target(repo: Path, env: Mapping[str, str]) -> Path:
    return repo / env.get("CARGO_TARGET_DIR", "target")


def _check_push(root: Path, plan: _CheckPlan, target: Path) -> None:
    for project in plan.python:
        _check_python(root, project)
    scope = _cargo_scope(root, plan)
    if scope:
        env = _cargo_environment(target, os.environ)
        for command in _cargo_commands(scope):
            _run(root, command, env=env)


def _clean_at_tip(repo: Path, tip: str) -> bool:
    return (
        _git(repo, ("rev-parse", "HEAD")).stdout.strip().decode() == tip
        and not _git(
            repo,
            (
                "status",
                "--porcelain=v1",
                "-z",
                "--untracked-files=normal",
                "--ignore-submodules=none",
            ),
        ).stdout
    )


def _validate(repo: Path, revisions: tuple[str, ...], tip: str | None = None) -> None:
    _run(repo, ("git", "diff", "--check", *revisions, "--"))
    raw = _git(repo, ("diff", "--raw", "--no-renames", "-z", *revisions, "--")).stdout
    plan = _plan(_parse_changes(raw))
    full_tip = tip if plan.needs_push_checks() else None
    if full_tip is not None and _clean_at_tip(repo, full_tip):
        # 通常の push は同じソースパスを使い、Cargo の増分ビルドを活かす。
        _check_files(repo, plan.files)
        _check_push(repo, plan, _repository_target(repo, os.environ))
    elif plan.files or full_tip is not None:
        with TemporaryDirectory(prefix="rhodium-hook-") as tmp:
            root = Path(tmp) / "src"
            root.mkdir()
            _checkout(repo, root, plan.files if full_tip is None else None, tip)
            _check_files(root, plan.files)
            if full_tip is not None:
                # 一時 checkout のテストは env!("CARGO_MANIFEST_DIR") などで一時パスを埋め込む。
                # リポジトリの target を共有すると、Cargo が新鮮と判断したその成果物を後の通常の
                # cargo test が再利用し、削除済みの一時パスを読んで失敗する。ビルド先を一時
                # ディレクトリの中に分け、checkout と一緒に消す。依存も作り直すので遅くなる。
                _check_push(root, plan, Path(tmp) / "target")


def _push_base(repo: Path, remote: str, tip: str, old: str) -> str | None:
    # 強制 push もリモートの旧ツリーと直接比較する。
    if old.strip("0"):
        return old
    history = _git(repo, ("rev-list", "--boundary", tip, "--not", "--remotes=" + remote)).stdout
    if not history:
        return None
    boundary = next(
        (row[1:].decode() for row in history.splitlines() if row.startswith(b"-")), None
    )
    if boundary is not None:
        return boundary
    return _git(repo, ("hash-object", "-t", "tree", "--stdin"), data=b"").stdout.strip().decode()


def _push_ranges(repo: Path, remote: str, updates: Iterable[str]) -> Iterator[tuple[str, str]]:
    for line in updates:
        _, local, _, old = line.split()
        if not local.strip("0"):
            continue
        # 注釈付きタグでも、タグの先にある commit のツリーを検査する。
        peeled = _git(repo, ("rev-parse", "-q", "--verify", local + "^{commit}"), check=False)
        if peeled.returncode == 0:
            tip = peeled.stdout.strip().decode()
            base = _push_base(repo, remote, tip, old)
            if base is not None:
                yield base, tip


def _main() -> None:
    repo = Path.cwd()
    match sys.argv[1:]:
        case ["pre-commit"]:
            _validate(repo, ("--cached",))
        case ["pre-push", remote]:
            for base, tip in dict.fromkeys(_push_ranges(repo, remote, sys.stdin)):
                _validate(repo, (base, tip), tip)
        case _:
            sys.exit("usage: checks.py pre-commit | pre-push <remote>")


if __name__ == "__main__":
    try:
        _main()
    except (OSError, ValueError, TypeError, SyntaxError, subprocess.CalledProcessError) as error:
        sys.exit(f"{sys.argv[1]}: {error}")
`````

`````
.githooks/pyproject.toml
[tool.ruff]
target-version = "py315"
line-length = 100

[tool.ruff.lint]
select = ["ALL"]
preview = true
ignore = [
    # フォーマッターが扱う空白、クォート、末尾のカンマ。
    "tab-indentation",
    "indentation-with-invalid-multiple",
    "indentation-with-invalid-multiple-comment",
    "over-indented",
    "incorrect-blank-line-before-class",
    "docstring-tab-indentation",
    "triple-single-quotes",
    "bad-quotes-inline-string",
    "bad-quotes-multiline-string",
    "bad-quotes-docstring",
    "avoidable-escaped-quote",
    "unnecessary-escaped-quote",
    "missing-trailing-comma",
    "prohibited-trailing-comma",
    # 複数行の docstring は 1 行目から書く。
    "multi-line-summary-second-line",
    # 著作権表示はリポジトリ単位で管理する。
    "missing-copyright-notice",
    # docstring の文末には日本語の句点を使う。
    "missing-trailing-period",
    "missing-terminal-punctuation",
]

[tool.ruff.lint.pycodestyle]
max-line-length = 120

[tool.ty.analysis]
respect-type-ignore-comments = false
strict-equality-semantics = true
strict-generic-narrowing = true

[tool.ty.environment]
python-version = "3.15"
python-platform = "all"

[tool.ty.rules]
all = "error"

[tool.ty.terminal]
error-on-warning = true
`````

`````
AGENTS.md
# rhodium 開発ガイド

rhodium は、金融商品を小さなコア IR で表し、価格・quote 感応度・XVA を計算する Rust ライブラリです。

## 設計書の読み方

初めて読むときは、[設計書の「読み方」](docs/design/README.md#読み方)から始める。[全体像](docs/design/01-overview.md)、[設計原則](docs/design/02-principles.md)、[アーキテクチャ](docs/design/03-architecture.md)を読んでから、担当する領域の章に進む。個々の作業では、次の文書を起点にする。

| 作業                   | 最初に確認すること                                                                                                                                                                                      |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 機能を実装する         | [ロードマップ](docs/design/12-roadmap.md)で対象の受け入れ基準（`M0-3` の形の ID）を特定し、関係する章の設計と「設計判断」を読む。                                                                       |
| 型や API を変更する    | [主要な型の地図](docs/design/03-architecture.md#主要な型の地図)で定義元を探し、[型の規約](docs/design/03-architecture.md#型の規約)と照らし合わせる。                                                      |
| 商品やモデルを追加する | [拡張の仕方](docs/design/03-architecture.md#拡張の仕方)で触る場所を決め、[拡張のチェックリスト](docs/design/11-quality.md#拡張のチェックリスト)で必要なテストと文書を確認する。                           |
| 設計を見直す           | 該当する文書の末尾の「設計判断」から元の判断と理由をたどる。旧設計からの経緯は [簡素化の計画](docs/SIMPLIFICATION_PLAN.md) にだけある。                                                                      |

## Milestone と issue で進める開発

- [Milestone 一覧](https://github.com/kakune/rhodium/milestones)から対象のマイルストーンを開き、Relationships の `Blocked by` に未完了の前提がない issue を選ぶ。並行して進められる作業は[ロードマップの依存図](docs/design/12-roadmap.md#全体像)に従う。
- Issue は独立して実装・検証できる成果ごとに作り、対象の milestone を設定する。本文には対象の crate やモジュール、実装する振る舞い、担当する受け入れ基準の ID を書き、設計書の節にリンクする。複数の issue で満たす基準は、最終的に合否を確認する issue を一つ決める。
- 前提があるときは後続 issue の `Blocked by` に直接の前提 issue を登録する。作業を分割したときや設計が変わったときは、関係する issue の担当基準、文書リンク、Relationships を更新する。
- 実装を完了する PR の本文には `Closes #番号`、担当する受け入れ基準、それを確かめるテストの名前と検証結果を書く。途中の PR では issue を参照し、残る作業が分かるようにする。
- マイルストーンの完了確認の issue に、受け入れ基準ごとの結果（テストの名前、実行結果、所要時間）を記録して閉じる。

## 実装と設計変更で守ること

- [設計原則](docs/design/02-principles.md)が衝突した場合は、番号の小さい原則を優先する。新しい型、trait、層、ファイル形式は、防ぐ誤りか得られる性能を具体的に示せるときだけ足す（P7）。型の名前と責務は設計書に従い、フィールドとメソッドの細部は rustdoc を正本とする。規約を変えるときは定義元の文書に判断と理由を記し、他の文書からはそこにリンクする。
- 計算は純粋関数として書き、入力と結果は不変に保つ。時計や乱数の seed、I/O、並列の実行器などの依存は引数で渡す。外部データは型付きの値を作る境界で一度だけ検証する。失敗しうる計算は `rhodium_core::Result<T>` で失敗を呼び出し元に返す。
- crate 間の依存は[crate の構成](docs/design/03-architecture.md#crate-の構成)の向きに限る。コア crate の外部依存は scioxide、serde、thiserror、chrono に限り（chrono はシステム時計を読む `clock` feature を使わない）、`#![forbid(unsafe_code)]` を適用する。ファイル I/O と CLI はファサードと `rhodium-cli` に置く。
- 営業日・日数計算には現地の `Date`、評価・fixing の公表・行使・担保決済には UTC の `Instant` を用いる。時刻と地域の解決は[日付とカレンダー](docs/design/04-foundation/01-dates-and-calendars.md#瞬間と現地時刻)に従い、評価に使うスナップショットは呼び出し側から受け取る。
- 商品の意味は[コア IR](docs/design/05-products/01-ir.md)に従う。標準商品は型付きの構造体として、独自の商品は DSL で表現する。組み込みのモデルやプライサーは閉じた enum に足し、利用者の拡張は trait を実装した値をジェネリックな API に渡す。
- 性能を著しく落とす簡素化をしない。[性能の規則](docs/design/03-architecture.md#性能の規則)の仕組み（準備の再利用、影響範囲だけの再校正、AAD の active set、パスのバッチ評価など）を保つ。
- 変換・抽出・集約には、`for` ループよりイテレータを優先する。列挙型の分岐には、`if` を重ねるより、各ケースを網羅する `match` を使う。単純な条件判定には `if` を、逐次処理の流れが明確になる場面では `for` を使う。
- doc に必要なコメントだけでなく、必要十分に、コードの意図や意味が伝わるコメントを書く。コメントは基本的に日本語だが、技術的な単語やコード中の変数名などをコメントで言及するときは英語でいい。テストやテストケースにおいても、何をテストしているか、その意図はなにか日本語コメントで記述する。
- 数学的な処理に関するところは、SciOxideに対応する実装がないか見て、ある場合積極的に使う。SciOxideへの開発やIssue起票などはしないが、必要そうであればユーザーに提案すること。
- Parse, don't validate を徹底する。外部入力や実行時に解決する情報は適切な境界で解析・検証し、不変条件が成立する型付きの値に変換する。内部 API は検証済みの型を受け取り、保証済みの条件を繰り返し検証しない。不変条件を破る値の構築・変更を防ぎ、未検証と検証済みのように中身の違う状態を型で区別する。実行時にしか判定できない制約や数値計算上の失敗は `Result` で扱う。
- 黙って誤った数字を出す演算は型で防ぐ。通貨、量の種類、FX の向きと受渡区分は[金額と量](docs/design/04-foundation/03-money-and-quantities.md)の型と許す演算に従い、ハンドルや添字の取り違えは[5 つのブランド lifetime](docs/design/03-architecture.md#ブランド)で防ぐ。`struct`、`enum`、newtype、phantom の型パラメータを適切に使うが、防ぐ誤りを示せない型の細分化（処理の段階ごとの typestate、次元の代数など）はしない。
- Rustの新機能 (1.100 の機能など) は有用な場合は積極的に使うこと。
- 様々なAPI や利用例は、使用するバージョンに合わせて Context7 で調べる。リポジトリ内の型・関数の定義や参照、trait の実装をたどるときは Serena を活用する。

## 検証

- 変更した振る舞いは、担当する受け入れ基準を確かめるテストで検証する。テストの名前は説明的な snake_case とし、オラクルと許容誤差は[品質](docs/design/11-quality.md)に従う。
- 数値の比較は[許容誤差](docs/design/11-quality.md#許容誤差)で行い、ビット単位の一致は求めない。並列の計算は、スレッド数や分割方法を変えても許容誤差の範囲で一致することを確かめる。
- 型で禁止する演算とブランドは、compile-fail テスト（trybuild）と正しい使い方のテストの対で確かめる。
- コミットの前に `cargo fmt --all --check`、`cargo clippy --locked --workspace --all-targets --all-features -- -D warnings`、`cargo test --locked --workspace --all-features`、`RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --all-features --no-deps` を通す。

## ブランチとプルリクエスト

- 最新の `main` から作業ブランチを切る。機能追加は `feat/...`、不具合修正は `fix/...`、振る舞いを変えない改修は `refactor/...`、文書変更は `docs/...` とし、`/` に続けて変更内容を短く書く。
- 作業と検証が済んだら、作業ブランチから `main` にプルリクエストを出す。本文には変更の目的と検証結果を記し、実行できなかった検証があれば理由を添える。
- コミットメッセージは `feat: ...` などのような形式の英語、PRのタイトルは英語、PRのbodyは日本語とする。
  - PR の body は、docs やコードへの参照に説明を委ねず、それ自体で変更の背景・目的・仕組み・影響を理解できるよう、自然な日本語で具体的かつ十分に詳しく記述すること。特に、どの変更箇所がどの機能追加・仕様変更に対応し、なぜ・どのように実装したのか、各テストが何を検証し何を保証するのかを明確に説明すること。可能な限りコミットごとに整理し、理解に必要なら PR の差分を超えた設計上の背景も補足すること。複雑な内容には数式・Mermaid 図・表・コード片・具体例などを適宜用い、単なる変更点やテスト名の列挙、過度な要約で済ませないこと。
- 意味の単位で頻繁にコミットする。コミットは適切なタイミングで自由に行ってよいので、**各コミットを小さく保つ**。 必要があればPRの前にコミットを分割して意味の単位に整理するなどする。
`````

`````
Cargo.toml
[workspace]
members = [
    "crates/rhodium-core",
    "crates/rhodium-products",
    "crates/rhodium-market",
    "crates/rhodium-pricing",
    "crates/rhodium-risk",
    "crates/rhodium-xva",
    "crates/rhodium",
    "crates/rhodium-cli",
    "crates/rhodium-testkit",
]
resolver = "3"

[workspace.package]
version = "0.1.0"
edition = "2024"
rust-version = "1.100.0"
license = "MIT OR Apache-2.0"
publish = false

[workspace.dependencies]
chrono = { version = "0.4.45", default-features = false, features = ["std"] }
criterion = { version = "0.8", default-features = false, features = ["cargo_bench_support"] }
csv = "1.4"
proptest = { version = "1.12", default-features = false, features = ["std"] }
scioxide-ad = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
scioxide-core = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false, features = ["std"] }
scioxide-interp = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
scioxide-linalg = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
scioxide-optimize = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
scioxide-rootfind = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
scioxide-sparse = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
rhodium = { path = "crates/rhodium" }
rhodium-core = { path = "crates/rhodium-core" }
rhodium-market = { path = "crates/rhodium-market" }
rhodium-pricing = { path = "crates/rhodium-pricing" }
rhodium-products = { path = "crates/rhodium-products" }
rhodium-risk = { path = "crates/rhodium-risk" }
rhodium-xva = { path = "crates/rhodium-xva" }
rhodium-testkit = { path = "crates/rhodium-testkit" }
serde = { version = "1", features = ["derive"] }
thiserror = "2"
trybuild = "1.0"

[workspace.lints.rust]
rust_2018_idioms = { level = "deny", priority = -1 }
missing_docs = "warn"
unfulfilled_lint_expectations = "deny"
unsafe_code = "forbid"
unused_must_use = "deny"

[workspace.lints.clippy]
all = { level = "deny", priority = -1 }
pedantic = { level = "warn", priority = -1 }
allow_attributes_without_reason = "deny"
cast_precision_loss = "allow"
dbg_macro = "deny"
disallowed_macros = "deny"
disallowed_methods = "deny"
expect_used = "deny"
many_single_char_names = "allow"
panic = "deny"
print_stderr = "deny"
print_stdout = "deny"
similar_names = "allow"
todo = "deny"
unimplemented = "deny"
unreachable = "deny"
unwrap_used = "deny"

[workspace.lints.rustdoc]
broken_intra_doc_links = "deny"
invalid_codeblock_attributes = "deny"
private_intra_doc_links = "deny"
`````

`````
README.md
# rhodium

rhodium は、金融商品を小さなコア IR で表し、価格、quote への感応度、XVA を、同じ market、同じ商品表現、同じ設定の上で計算する Rust ライブラリです。
数値計算の基盤には [scioxide](https://github.com/kakune/scioxide) を使います。

標準商品は型付きの Rust 構造体として、独自の商品は小さな DSL で書き、どちらも同じ IR に変換します。
曲線とプライサーは数値型 `R` に対するジェネリックなコードとして一度だけ書き、f64 で値を、AAD で微分を求めます。
quote への感応度は、取引の θ 勾配を集計してから校正の随伴で quote まで戻して求めます。
通貨の違う金額の和や、別の market の曲線ハンドルの使用のように黙って誤った数字を出す誤りは、型とブランド lifetime でコンパイル時に拒否します。

設計書は [docs/design/](docs/design/README.md) にあります。
旧設計からの簡素化と、milestone、issue、コードの移行の手順は [docs/SIMPLIFICATION_PLAN.md](docs/SIMPLIFICATION_PLAN.md) にあります。

## 構成

設計上の crate は次の 8 つです（[アーキテクチャ](docs/design/03-architecture.md)）。

| crate | 責務 |
|---|---|
| `rhodium-core` | 日付と時刻、カレンダー、日数計算、スケジュール、通貨と金額、量の型、参照データ、fixing、エラー、並列の分割 |
| `rhodium-products` | コア IR、DSL、標準商品、取引と担保契約 |
| `rhodium-market` | quote、市場の構造と θ、market view、曲線、ボラ面 |
| `rhodium-pricing` | 価格設定、校正、リニア評価、解析プライサー、モデル、モンテカルロと AMC |
| `rhodium-risk` | quote 感応度、ガンマ、シナリオ、theta と P&L 説明、付け替え、SIMM、FRTB-SA |
| `rhodium-xva` | エクスポージャー、担保、XVA の指標と感応度 |
| `rhodium` | ファサード、`Workspace`、`Session` |
| `rhodium-cli` | `rh` コマンド |

## 現状

マイルストーン M1（金利リニアの完成）の受け入れ基準 M1-1〜M1-8 を満たしました（[ロードマップ](docs/design/12-roadmap.md)、Issue #37）。
JPY（TONA の会合日のステップと年末のジャンプ、TIBOR 6M・3M のベーシス）、USD（SOFR）、USDJPY（FX の短期と MtM XCCY の長期による外貨担保の曲線）を一つの市場の定義で校正し、
FRA、金利先物、IRS、ベーシススワップ、OIS（lookback、observation shift、lockout、支払いの遅れ）、FX フォワードと FX スワップ、XCCY の PV と quote デルタを、CSA の 4 つの場合の割引で求めます。
所要時間は [M0 の計測](docs/benchmarks/m0.md) と [M1 の計測](docs/benchmarks/m1.md) にあります。
次のマイルストーンは M2（ボラと解析オプション）です。

ワークスペースは上の 8 crate と開発用の `rhodium-testkit` からなり、旧設計の crate は残っていません。

## サンプル

### M1 の金利リニア

最初の検証データセットの金利商品を評価するワークスペースです。[市場の定義と取引の説明](examples/m1-rates/README.md)

```sh
cargo run --locked -p rhodium-cli -- price examples/m1-rates --date 2026-11-24
cargo run --locked -p rhodium-cli -- risk examples/m1-rates --date 2026-11-24 --threads 4
```

### M0 の OIS

JPY TONA と USD SOFR の OIS を評価するサンプルワークスペースを同梱しています。

```sh
cargo run --locked -p rhodium-cli -- price examples/m0-ois --date 2026-11-24
cargo run --locked -p rhodium-cli -- risk examples/m0-ois --as-of 2026-11-24T17:00:00+09:00
```

この日付指定は価格設定の東京 17 時を使います。
コマンドは、取引別の PV または全取引の通貨別 quote デルタを CSV で標準出力に返します。[サンプルの入力と計算](examples/m0-ois/README.md)
`--threads N` を付けると取引を N 個のワーカーに分けて計算します（既定は 1 で逐次）。

## 開発環境

NixOS では、flake の開発シェルに Rust、Cargo、Clippy、rustfmt、Git、Python 3.15、Ruff、ty が揃います。
Rust のバージョンとコンポーネントは `rust-toolchain.toml`、Nix の入力は `flake.lock` で固定しています。
開発シェルでは `RUST_SRC_PATH` に同じツールチェーンの標準ライブラリのソースパスが設定されます。

direnv のシェル連携を有効にした環境では、リポジトリのルートで初回に `direnv allow` を実行します。
以降はディレクトリに入ると開発シェルの環境が読み込まれます。
`nix develop` でも同じ開発シェルに入れます。

## 開発時の検証

`rust-toolchain.toml` が指定するツールチェーンで、リポジトリのルートから実行します。

```sh
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets --all-features -- -D warnings
cargo test --locked --workspace --all-features
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --all-features --no-deps
```

数値の比較は [許容誤差](docs/design/11-quality.md#許容誤差) で行い、ビット単位の一致は求めません。

### Git フック

Git フックを有効にするには、クローンごとに `git config core.hooksPath .githooks` を実行します。
コミット時はステージした差分の空白と、変更ファイルの Rust 書式、Python・TOML・JSON・Nix・フックの構文を検査します。
push 時は送信先との差分に同じ検査を適用し、変更した crate の Clippy、テスト、rustdoc を実行します。
作業ツリーに送らない変更があると、送る commit を一時ディレクトリに取り出し、ビルド先もその中に分けて検査します（依存から作り直すので時間がかかります）。
ワークスペース共通の Cargo 設定を変えた場合は全 crate が対象です。
`.githooks/` の Python の変更には、その `pyproject.toml` に設定した Ruff と ty を適用します。
新規ブランチでは送信先の追跡 ref と共有する履歴を基準にし、共有履歴がなければ初回の push としてツリー全体を検査します。

## ライセンス

rhodium は、利用者が選択する [MIT License](LICENSE-MIT) または
[Apache License 2.0](LICENSE-APACHE) の条件で利用できます。
`````

`````
clippy.toml
msrv = "1.100.0"
allow-unwrap-in-tests = true
allow-expect-in-tests = true
allow-panic-in-tests = true

disallowed-macros = [
    { path = "std::env", reason = "Read environment values at the adapter boundary." },
    { path = "std::option_env", reason = "Read environment values at the adapter boundary." },
]

disallowed-methods = [
    { path = "std::time::Instant::now", reason = "Pass the time as an explicit input." },
    { path = "std::time::SystemTime::now", reason = "Pass the time as an explicit input." },
    { path = "chrono::Utc::now", reason = "Pass the time as an explicit input.", allow-invalid = true },
    { path = "chrono::Local::now", reason = "Pass the time as an explicit input.", allow-invalid = true },
    { path = "std::env::args", reason = "Read platform settings at the adapter boundary." },
    { path = "std::env::args_os", reason = "Read platform settings at the adapter boundary." },
    { path = "std::env::var", reason = "Read platform settings at the adapter boundary." },
    { path = "std::env::var_os", reason = "Read platform settings at the adapter boundary." },
    { path = "std::env::vars", reason = "Read platform settings at the adapter boundary." },
    { path = "std::env::vars_os", reason = "Read platform settings at the adapter boundary." },
    { path = "std::env::current_dir", reason = "Read platform settings at the adapter boundary." },
    { path = "std::env::current_exe", reason = "Read platform settings at the adapter boundary." },
    { path = "std::env::home_dir", reason = "Read platform settings at the adapter boundary." },
    { path = "std::env::temp_dir", reason = "Read platform settings at the adapter boundary." },
]
`````

`````
crates/rhodium-cli/Cargo.toml
[package]
name = "rhodium-cli"
description = "rhodium のコマンドラインインターフェース"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
publish.workspace = true

[lints]
workspace = true

[[bin]]
name = "rh"
path = "src/main.rs"

[dependencies]
clap = { version = "4.6", features = ["derive"] }
csv.workspace = true
rhodium.workspace = true

[dev-dependencies]
rhodium-testkit.workspace = true
tempfile = "3.27"
`````

`````
crates/rhodium-cli/src/main.rs
//! ワークスペースの PV と quote デルタを CSV で出力する。

use std::{
    error::Error,
    io::{self, Write},
    num::NonZeroUsize,
    path::PathBuf,
    process::ExitCode,
};

use clap::{Args, Parser, Subcommand};
use rhodium::{
    prelude::{ExecutionContext, ThreadExecutor},
    session::{GroupOutcome, Session},
    workspace::Workspace,
};

#[derive(Parser)]
#[command(
    version,
    about = "Price trades and compute quote delta from a rhodium workspace"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Write PV rows by trade and currency.
    Price(ValuationArgs),
    /// Write portfolio delta by quote and currency, in each quote's reporting unit.
    Risk(ValuationArgs),
}

#[derive(Args)]
struct ValuationArgs {
    /// Workspace directory.
    dir: PathBuf,
    #[command(flatten)]
    time: ValuationTime,
    /// Worker threads that split the trades; 1 runs sequentially.
    #[arg(long, value_name = "N", default_value = "1")]
    threads: NonZeroUsize,
}

#[derive(Args)]
#[group(required = true, multiple = false)]
struct ValuationTime {
    /// Instant with a UTC offset or a time zone in brackets.
    #[arg(long, value_name = "DATETIME")]
    as_of: Option<String>,
    /// Local date evaluated at the configured `valuation_close`.
    #[arg(long, value_name = "YYYY-MM-DD")]
    date: Option<String>,
}

// 終了コードと標準エラーに書くエラー。
type CommandError = (ExitCode, Box<dyn Error>);

// 入力ファイル、引数、評価時点の誤り。CSV を一行も出さずに 2 で終える。
fn input_error(error: impl Into<Box<dyn Error>>) -> CommandError {
    (ExitCode::from(2), error.into())
}

// 評価や出力の失敗。部分的に計算できた行を出した後の失敗と同じ 1 で終える。
fn execution_error(error: impl Into<Box<dyn Error>>) -> CommandError {
    (ExitCode::FAILURE, error.into())
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let mut stderr = io::stderr().lock();
    match run(&cli.command, io::stdout().lock(), &mut stderr) {
        Ok(code) => code,
        Err((code, error)) => {
            let _ = writeln!(stderr, "error: {error}");
            code
        }
    }
}

// 入力を読み、評価時点を解決し、校正と取引の診断を書いてから CSV を出す。
// 一部の校正や取引が失敗しても計算できた行は出し、戻り値の終了コードを 1 にする。
fn run(
    command: &Command,
    output: impl Write,
    diagnostics: &mut impl Write,
) -> Result<ExitCode, CommandError> {
    let args = match command {
        Command::Price(args) | Command::Risk(args) => args,
    };
    // 取引を分けるワーカーの数。既定の 1 は逐次で、スレッドを作らない。
    let executor = ThreadExecutor::new(args.threads);
    let execution = ExecutionContext::new(&executor);
    let workspace =
        Workspace::load_with_execution(&args.dir, execution).map_err(|error| match error {
            rhodium::prelude::Error::Execution(_) | rhodium::prelude::Error::Internal(_) => {
                execution_error(error)
            }
            _ => input_error(error),
        })?;
    let as_of = match (&args.time.as_of, &args.time.date) {
        (Some(at), None) => workspace.resolve_as_of(at),
        (None, Some(date)) => date.parse().and_then(|date| workspace.resolve_date(date)),
        _ => return Err(input_error("specify one of --as-of and --date")),
    }
    .map_err(input_error)?;
    // 校正グループ単位の失敗は GroupOutcome に残るため、Session::new 全体の失敗は
    // 評価時点が曲線のピラー以降にあるなど、入力どうしの不整合で市場を組み立てられない場合に当たる。
    // 計算できた行が一つも無いので、部分的な失敗（1）ではなく入力のエラー（2）として返す。
    // 並列の実行の失敗（Execution など）は入力の誤りではないので 1 にする。
    let session = Session::new(workspace, as_of, execution).map_err(|error| match error {
        rhodium::prelude::Error::InvalidInput(_)
        | rhodium::prelude::Error::InvalidData(_)
        | rhodium::prelude::Error::MissingData(_) => input_error(error),
        _ => execution_error(error),
    })?;

    let mut failed = false;
    for warning in session.calibrated_market().warnings() {
        writeln!(diagnostics, "warning: calibration: {warning:?}").map_err(execution_error)?;
    }
    for group in session.calibrated_market().groups() {
        let reason = match group.outcome() {
            GroupOutcome::Solved(_) => continue,
            GroupOutcome::Failed(error) => error.to_string(),
            GroupOutcome::Unavailable => "upstream calibration unavailable".into(),
        };
        let curves = group
            .blocks()
            .iter()
            .map(|block| block.curve().as_str())
            .collect::<Vec<_>>()
            .join(", ");
        writeln!(diagnostics, "error: calibration for {curves}: {reason}")
            .map_err(execution_error)?;
        failed = true;
    }
    for (id, warnings) in session.preparation_warnings() {
        for warning in warnings {
            writeln!(diagnostics, "warning: trade {id}: {warning:?}").map_err(execution_error)?;
        }
    }

    let mut writer = csv::WriterBuilder::new()
        .has_headers(false)
        .from_writer(output);
    failed |= match command {
        Command::Price(_) => price(&session, execution, &mut writer, diagnostics)?,
        Command::Risk(_) => risk(&session, execution, &mut writer, diagnostics)?,
    };
    writer.flush().map_err(execution_error)?;
    Ok(if failed {
        ExitCode::FAILURE
    } else {
        ExitCode::SUCCESS
    })
}

// 取引 ID 順、通貨順に PV を書く。失敗した取引は診断に回し、失敗の有無を返す。
fn price(
    session: &Session,
    execution: ExecutionContext<'_>,
    writer: &mut csv::Writer<impl Write>,
    diagnostics: &mut impl Write,
) -> Result<bool, CommandError> {
    // 計算に失敗しても終了コード 1 の出力は表の形を保つよう、ヘッダーを先に書く。
    writer
        .write_record(["trade_id", "currency", "pv"])
        .map_err(execution_error)?;
    let mut failed = false;
    for (trade, result) in session.pv(execution).map_err(execution_error)? {
        match result {
            Ok(pv) => {
                for money in pv.terms() {
                    writer
                        .serialize((&trade, money.ccy(), money.amount()))
                        .map_err(execution_error)?;
                }
            }
            Err(error) => {
                writeln!(diagnostics, "error: trade {trade}: {error}").map_err(execution_error)?;
                failed = true;
            }
        }
    }
    Ok(failed)
}

// 全取引の quote デルタを quote ID 順、通貨順に書く。集計から除いた取引は診断に回す。
fn risk(
    session: &Session,
    execution: ExecutionContext<'_>,
    writer: &mut csv::Writer<impl Write>,
    diagnostics: &mut impl Write,
) -> Result<bool, CommandError> {
    // 計算に失敗しても終了コード 1 の出力は表の形を保つよう、ヘッダーを先に書く。
    writer
        .write_record(["quote_id", "currency", "sensitivity"])
        .map_err(execution_error)?;
    let report = session.quote_delta(execution).map_err(execution_error)?;
    let mut rows: Vec<_> = report
        .by_ccy()
        .iter()
        .flat_map(|vector| {
            vector
                .entries()
                .iter()
                .map(|(quote, value)| (quote, vector.ccy(), value))
        })
        .collect();
    rows.sort_unstable_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
    for row in rows {
        writer.serialize(row).map_err(execution_error)?;
    }
    for (trade, error) in report.failures() {
        writeln!(diagnostics, "error: trade {trade}: {error}").map_err(execution_error)?;
    }
    Ok(!report.failures().is_empty())
}
`````

`````
crates/rhodium-cli/tests/workspace_cli.rs
//! 実ファイルのワークスペースをファサードと rh から評価する。

use std::{
    error::Error,
    fmt::Write as _,
    fs, io,
    num::NonZeroUsize,
    path::{Path, PathBuf},
    process::{Command, Output},
};

use rhodium::{
    prelude::{ExecutionContext, ThreadExecutor},
    session::Session,
    workspace::Workspace,
};
use rhodium_testkit::assert_close;
use tempfile::TempDir;

type TestResult<T = ()> = Result<T, Box<dyn Error>>;
type Rows = Vec<(String, String, f64)>;

#[expect(
    clippy::disallowed_macros,
    reason = "Cargo が渡す、試験データへのパスを取得する"
)]
fn sample() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m0-ois")
}

fn copy_dir(from: &Path, to: &Path) -> io::Result<()> {
    for entry in fs::read_dir(from)? {
        let entry = entry?;
        let target = to.join(entry.file_name());
        if entry.file_type()?.is_dir() {
            fs::create_dir(&target)?;
            copy_dir(&entry.path(), &target)?;
        } else {
            fs::copy(entry.path(), target)?;
        }
    }
    Ok(())
}

fn workspace_copy() -> io::Result<TempDir> {
    let dir = tempfile::tempdir()?;
    copy_dir(&sample(), dir.path())?;
    Ok(dir)
}

#[expect(
    clippy::disallowed_macros,
    reason = "Cargo がビルドした試験対象の実行ファイルを使う"
)]
fn run(command: &str, dir: &Path, args: &[&str]) -> io::Result<Output> {
    Command::new(env!("CARGO_BIN_EXE_rh"))
        .arg(command)
        .arg(dir)
        .args(args)
        .output()
}

fn rows(output: &Output, headers: &[&str]) -> TestResult<Rows> {
    let mut reader = csv::Reader::from_reader(output.stdout.as_slice());
    assert!(reader.headers()?.iter().eq(headers.iter().copied()));
    reader
        .records()
        .map(|row| {
            let row = row?;
            Ok((row[0].to_owned(), row[1].to_owned(), row[2].parse()?))
        })
        .collect()
}

fn diagnostic(dir: &Path) -> TestResult<String> {
    load_diagnostic(Workspace::load(dir))
}

fn load_diagnostic(workspace: rhodium::prelude::Result<Workspace>) -> TestResult<String> {
    match workspace {
        Err(error @ rhodium::prelude::Error::InvalidData(_)) => Ok(error.to_string()),
        Err(error) => Err(error.into()),
        Ok(_) => Err(io::Error::other("expected Workspace::load to fail").into()),
    }
}

fn ois(notional: f64, accrual: f64, fixed: f64, quote: f64) -> (f64, f64) {
    // 期首で評価する単一支払いの OIS は DF = 1 / (1 + q α)。
    // 固定受けの PV を q で微分し、金利 1 bp あたりに直す。
    let df = (1.0 + quote * accrual).recip();
    (
        notional * accrual * (fixed - quote) * df,
        -notional * accrual * (1.0 + fixed * accrual) * df.powi(2) * 1e-4,
    )
}

#[test]
fn price_and_risk_match_for_as_of_and_valuation_close() -> TestResult {
    // 取引と quote の入力順は逆にしておき、CSV が ID 順になることも確かめる。
    let expected = [
        (
            "JPY-SWAP",
            "JPY.TONA.OIS.1M",
            "JPY",
            ois(1e8, 30.0 / 365.0, 0.006, 0.0045),
            1e-6,
        ),
        (
            "USD-SWAP",
            "USD.SOFR.OIS.1M",
            "USD",
            ois(1e6, 30.0 / 360.0, 0.04, 0.035),
            1e-8,
        ),
    ];
    // 同じ瞬間（東京 17 時 = UTC 08 時）を、オフセット付き、地域付き、締め時刻付きの日付で指定する。
    let times: [&[&str]; 3] = [
        &["--as-of", "2026-11-24T17:00:00+09:00"],
        &["--as-of", "2026-11-24T17:00:00[Asia/Tokyo]"],
        &["--date", "2026-11-24"],
    ];
    for (command, headers, relative) in [
        ("price", ["trade_id", "currency", "pv"], 1e-12),
        ("risk", ["quote_id", "currency", "sensitivity"], 1e-10),
    ] {
        let outputs = times
            .iter()
            .map(|args| run(command, &sample(), args))
            .collect::<io::Result<Vec<_>>>()?;
        let tables = outputs
            .iter()
            .map(|output| {
                assert!(
                    output.status.success(),
                    "{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                assert_eq!(output.stderr, [] as [u8; 0]);
                rows(output, &headers)
            })
            .collect::<TestResult<Vec<_>>>()?;
        let [first, ..] = tables.as_slice() else {
            return Err(io::Error::other("no CLI output").into());
        };
        for table in &tables {
            assert_eq!(table.len(), expected.len());
            for (
                ((id, ccy, value), (first_id, first_ccy, first_value)),
                (trade, quote, oracle_ccy, (pv, delta), absolute),
            ) in table.iter().zip(first).zip(expected)
            {
                let (oracle_id, oracle) = if command == "price" {
                    (trade, pv)
                } else {
                    (quote, delta)
                };
                assert_eq!((id.as_str(), ccy.as_str()), (oracle_id, oracle_ccy));
                assert_eq!((id, ccy), (first_id, first_ccy));
                assert_close!(*value, *first_value, relative, absolute);
                assert_close!(*value, oracle, relative, absolute);
            }
        }
    }
    Ok(())
}

#[test]
fn cli_sums_quote_deltas_over_trades_of_a_currency() -> TestResult {
    // 同じ通貨の取引を足し、quote デルタが取引ごとの閉形式のデルタの和になることを確かめる。
    // 追加の JPY 取引は元本 5,000 万円、固定 0.50% の払いで、デルタの符号が受けと逆になる。
    let dir = workspace_copy()?;
    fs::write(
        dir.path().join("trades/more.toml"),
        r#"version = 1

[[trade]]
id = "JPY-SWAP-2"
margin_set = "JPY-CASH"
book = "JPY-SWAPS"
trade_date = 2026-11-24

[trade.product]
type = "swap"
convention = "JPY-TONA-OIS"
start = 2026-11-24
tenor = "1M"
notional = "50_000_000 JPY"
fixed_rate = "0.50%"
fixed_side = "pay"
"#,
    )?;
    let accrual = 30.0 / 365.0;
    let (receive_pv, receive_delta) = ois(1e8, accrual, 0.006, 0.0045);
    let (pay_pv, pay_delta) = ois(5e7, accrual, 0.005, 0.0045);
    let (usd_pv, usd_delta) = ois(1e6, 30.0 / 360.0, 0.04, 0.035);

    let price = run("price", dir.path(), &["--date", "2026-11-24"])?;
    assert!(
        price.status.success(),
        "{}",
        String::from_utf8_lossy(&price.stderr)
    );
    let pv = rows(&price, &["trade_id", "currency", "pv"])?;
    let expected_pv = [
        ("JPY-SWAP", "JPY", receive_pv, 1e-6),
        ("JPY-SWAP-2", "JPY", -pay_pv, 1e-6),
        ("USD-SWAP", "USD", usd_pv, 1e-8),
    ];
    assert_eq!(pv.len(), expected_pv.len());
    for ((id, ccy, value), (oracle_id, oracle_ccy, oracle, absolute)) in pv.iter().zip(expected_pv)
    {
        assert_eq!((id.as_str(), ccy.as_str()), (oracle_id, oracle_ccy));
        assert_close!(*value, oracle, 1e-12, absolute);
    }

    let risk = run("risk", dir.path(), &["--date", "2026-11-24"])?;
    assert!(
        risk.status.success(),
        "{}",
        String::from_utf8_lossy(&risk.stderr)
    );
    let delta = rows(&risk, &["quote_id", "currency", "sensitivity"])?;
    assert_eq!(delta.len(), 2);
    assert_eq!(
        (delta[0].0.as_str(), delta[0].1.as_str()),
        ("JPY.TONA.OIS.1M", "JPY")
    );
    // 打ち消し合う和なので、相対誤差の尺度に項の絶対値の和を使う。
    let gross = receive_delta.abs() + pay_delta.abs();
    assert_close!(delta[0].2, receive_delta - pay_delta, 1e-10, 1e-6, gross);
    assert_eq!(
        (delta[1].0.as_str(), delta[1].1.as_str()),
        ("USD.SOFR.OIS.1M", "USD")
    );
    assert_close!(delta[1].2, usd_delta, 1e-10, 1e-8);
    Ok(())
}

#[test]
fn cli_reports_workspace_load_errors_as_input_errors() -> TestResult {
    // 読み込みで見つかる誤り（未知のフィールド、解析できない率）は入力のエラー 2 で、CSV を一行も出さない。
    // 解析できない率は、どの取引のどのフィールドの何が誤りかを報告する（以前は SwapInput の
    // どの形にも一致しないとだけ表示し、原因が分からなかった）。
    for (from, to, expected) in [
        (
            "book = \"JPY-SWAPS\"",
            "book = \"JPY-SWAPS\"\ndesk = \"RATES\"",
            vec!["unknown field `desk`"],
        ),
        (
            "fixed_rate = \"4%\"",
            "fixed_rate = \"invalid_rate\"",
            vec![
                "trade[0]",
                "`fixed_rate`",
                "rate must end with % or bp",
                "product",
            ],
        ),
    ] {
        let dir = workspace_copy()?;
        let path = dir.path().join("trades/swaps.toml");
        let trades = fs::read_to_string(&path)?;
        assert!(trades.contains(from));
        fs::write(&path, trades.replacen(from, to, 1))?;
        for command in ["price", "risk"] {
            let output = run(command, dir.path(), &["--date", "2026-11-24"])?;
            let error = String::from_utf8_lossy(&output.stderr);
            assert_eq!(output.status.code(), Some(2), "{error}");
            assert!(
                error.starts_with("error:") && error.contains("swaps.toml"),
                "{error}"
            );
            assert!(expected.iter().all(|part| error.contains(part)), "{error}");
            assert_eq!(output.stdout, [] as [u8; 0]);
        }
    }
    Ok(())
}

#[test]
fn session_reuses_the_same_instant_for_price_and_risk() -> TestResult {
    // 同じセッションから両方を呼び出し、取引 ID と通貨別の結果を確認する。
    fn send_sync<T: Send + Sync>() {}
    send_sync::<Session>();
    let workspace = Workspace::load(sample())?;
    let as_of = workspace.resolve_date("2026-11-24".parse()?)?;
    assert_eq!(as_of.to_string(), "2026-11-24T08:00:00Z");
    let execution = ExecutionContext::sequential();
    let session = Session::new(workspace, as_of, execution)?;
    assert_eq!(session.as_of(), as_of);
    assert_eq!(session.calibrated_market().groups().len(), 2);
    let prices = session.pv(execution)?;
    assert_eq!(prices.len(), 2);
    assert!(prices.values().all(Result::is_ok));
    let risk = session.quote_delta(execution)?;
    assert!(risk.failures().is_empty());
    assert_eq!(risk.by_ccy().len(), 2);
    Ok(())
}

#[test]
fn workspace_reports_multiple_csv_errors_with_lines() -> TestResult {
    // レコード長の違い、重複 ID、非有限値のすべてを一度に報告する。
    let dir = workspace_copy()?;
    fs::write(
        dir.path().join("quotes/snapshot.csv"),
        "quote_id,value\nJPY.TONA.OIS.1M,45bp\nJPY.TONA.OIS.1M,0.006\nTOO-SHORT\nUSD.SOFR.OIS.1M,NaN\n",
    )?;
    let error = diagnostic(dir.path())?;
    assert!(error.contains("quotes/snapshot.csv"), "{error}");
    assert!(
        error.contains("line 3") && error.contains("duplicate ID"),
        "{error}"
    );
    assert!(error.contains("line: 4"), "{error}");
    assert!(
        error.contains("line 5") && error.contains("value: expected a finite decimal"),
        "{error}"
    );
    Ok(())
}

#[test]
fn workspace_reports_multiple_toml_rows_and_duplicate_trade_ids() -> TestResult {
    // 同じファイルの別々の不正な行と、ファイルをまたぐ取引 ID の重複を検出する。
    let dir = workspace_copy()?;
    let path = dir.path().join("trades/swaps.toml");
    let later = dir.path().join("trades/zz-duplicate.toml");
    let good = fs::read_to_string(&path)?;
    fs::write(&later, "]]")?;
    fs::write(
        &path,
        good.replace("book = \"USD-SWAPS\"", "bok = \"USD-SWAPS\"")
            .replace(
                "margin_set = \"JPY-CASH\"",
                "margin_set = \"NO-SUCH-MARGIN\"",
            ),
    )?;
    let error = diagnostic(dir.path())?;
    // 後の文書が壊れていても、先に読み込んだ文書の行エラーを元の順で返す。
    let executor = ThreadExecutor::new(NonZeroUsize::MIN.saturating_add(3));
    let parallel = || {
        load_diagnostic(Workspace::load_with_execution(
            dir.path(),
            ExecutionContext::new(&executor),
        ))
    };
    assert_eq!(parallel()?, error);
    assert!(error.contains("trades/swaps.toml:3:1: trade[0]"), "{error}");
    assert!(error.contains("unknown field `bok`"), "{error}");
    assert!(
        error.contains("trades/swaps.toml:18:1: trade[1]"),
        "{error}"
    );
    assert!(
        error.contains("margin_set: unknown margin set NO-SUCH-MARGIN"),
        "{error}"
    );
    fs::write(&path, &good)?;
    fs::write(later, good)?;
    let error = diagnostic(dir.path())?;
    assert_eq!(parallel()?, error);
    assert!(
        error.contains("zz-duplicate.toml") && error.contains("duplicate ID JPY-SWAP"),
        "{error}"
    );
    Ok(())
}

#[test]
fn workspace_rejects_invalid_version_and_fixing_header() -> TestResult {
    // 版を解釈できない文書と、日次指標に日中の列を書いたファイルを拒否する。
    let dir = workspace_copy()?;
    let path = dir.path().join("agreements.toml");
    let good = fs::read_to_string(&path)?;
    fs::write(&path, good.replace("version = 1", "version = 2"))?;
    let error = diagnostic(dir.path())?;
    assert!(
        error.contains("agreements.toml") && error.contains("unsupported version 2"),
        "{error}"
    );
    fs::write(&path, good)?;
    fs::write(
        dir.path().join("fixings/JPY-TONA.csv"),
        "at,value\n2026-11-20T08:00:00Z,0.0046\n",
    )?;
    let error = diagnostic(dir.path())?;
    assert!(
        error.contains("fixings/JPY-TONA.csv") && error.contains("expected CSV header date,value"),
        "{error}"
    );
    Ok(())
}

#[test]
fn workspace_reports_reference_errors_by_source_entry() -> TestResult {
    // 別々の指標が参照する未知カレンダーを、入力ファイルと指標 ID で報告する。
    let dir = workspace_copy()?;
    let path = dir.path().join("refdata/indices.toml");
    let indices = fs::read_to_string(&path)?
        .replace("calendar = \"JPTO\"", "calendar = \"UNKNOWN-JP\"")
        .replace("calendar = \"USNY\"", "calendar = \"UNKNOWN-US\"");
    fs::write(path, indices)?;
    let error = diagnostic(dir.path())?;
    assert!(
        error.contains("refdata/indices.toml: index[JPY-TONA]") && error.contains("UNKNOWN-JP"),
        "{error}"
    );
    assert!(
        error.contains("refdata/indices.toml: index[USD-SOFR]") && error.contains("UNKNOWN-US"),
        "{error}"
    );
    Ok(())
}

#[test]
fn cli_requires_exactly_one_valuation_time() -> TestResult {
    // clap の引数エラーで入力の読み込みや CSV の出力に進まないことを確認する。
    for command in ["price", "risk"] {
        for args in [
            &[][..],
            &["--date", "2026-11-24", "--as-of", "2026-11-24T08:00:00Z"][..],
        ] {
            let output = run(command, &sample(), args)?;
            let error = String::from_utf8_lossy(&output.stderr);
            assert_eq!(output.status.code(), Some(2), "{error}");
            assert!(
                error.contains("--date") && error.contains("--as-of"),
                "{error}"
            );
            assert_eq!(output.stdout, [] as [u8; 0]);
        }
    }
    Ok(())
}

#[test]
fn workspace_rejects_ambiguous_or_missing_local_times() -> TestResult {
    // 同じ地域の重複時刻は明示オフセットで区別し、欠落時刻は拒否する。
    let workspace = Workspace::load(sample())?;
    let instant = workspace.resolve_as_of("2026-11-24T17:00:00[Asia/Tokyo]")?;
    assert_eq!(
        instant,
        workspace.resolve_as_of("2026-11-24T17:00:00+09:00[Asia/Tokyo]")?
    );
    assert!(
        workspace
            .resolve_as_of("2026-11-24T17:00:00+08:00[Asia/Tokyo]")
            .is_err()
    );
    assert!(workspace.resolve_as_of("2026-11-24T17:00:00").is_err());
    assert!(
        workspace
            .resolve_as_of("2026-03-08T02:30:00[America/New_York]")
            .is_err()
    );
    assert!(
        workspace
            .resolve_as_of("2026-11-01T01:30:00[America/New_York]")
            .is_err()
    );
    let early = workspace.resolve_as_of("2026-11-01T01:30:00-04:00[America/New_York]")?;
    let late = workspace.resolve_as_of("2026-11-01T01:30:00-05:00[America/New_York]")?;
    assert_eq!(early.to_string(), "2026-11-01T05:30:00Z");
    assert_eq!(late.to_string(), "2026-11-01T06:30:00Z");
    // 参照データに無い地域は、他の地域の参照と同じく参照の欠けとして扱う。
    assert!(matches!(
        workspace.resolve_as_of("2026-11-24T17:00:00[Europe/Paris]"),
        Err(rhodium::prelude::Error::InvalidData(message)) if message.contains("Europe/Paris")
    ));
    Ok(())
}

#[test]
fn date_requires_an_unambiguous_configured_close() -> TestResult {
    // 日付指定には締め時刻が必要で、夏時間の重複と欠落は計算に進む前に拒否する。
    let dir = workspace_copy()?;
    let path = dir.path().join("config/pricing.toml");
    let config = fs::read_to_string(&path)?;
    let close = "valuation_close = { time = \"17:00\", zone = \"Asia/Tokyo\" }";
    for (replacement, date, diagnostic) in [
        ("", "2026-11-24", "requires valuation_close"),
        (
            "valuation_close = { time = \"01:30\", zone = \"America/New_York\" }",
            "2026-11-01",
            "ambiguous local time",
        ),
        (
            "valuation_close = { time = \"02:30\", zone = \"America/New_York\" }",
            "2026-03-08",
            "transition gap",
        ),
    ] {
        fs::write(&path, config.replace(close, replacement))?;
        let output = run("price", dir.path(), &["--date", date])?;
        let error = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(2), "{error}");
        assert!(error.contains(diagnostic), "{error}");
        assert_eq!(output.stdout, [] as [u8; 0]);
    }
    fs::write(path, config.replace(close, ""))?;
    let output = run("price", dir.path(), &["--as-of", "2026-11-24T08:00:00Z"])?;
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    Ok(())
}

#[test]
fn cli_reports_valuation_time_after_curve_pillars_as_input_error() -> TestResult {
    // サンプルの曲線のピラーは 2026-12-24T08:00Z だけなので、翌日の締め時刻では市場を組み立てられない。
    // 評価時点と市場定義の不整合は入力のエラー（2）で、部分的な失敗（1）と違い CSV を一行も出さない。
    for command in ["price", "risk"] {
        let output = run(command, &sample(), &["--date", "2026-12-25"])?;
        let error = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(2), "{error}");
        assert!(error.contains("pillar times after as-of"), "{error}");
        assert_eq!(output.stdout, [] as [u8; 0]);
    }
    Ok(())
}

#[test]
fn cli_keeps_independent_results_when_calibration_fails() -> TestResult {
    // TONA の quote を欠かせても、SOFR の PV とデルタは出力する。
    let dir = workspace_copy()?;
    fs::write(
        dir.path().join("quotes/snapshot.csv"),
        "quote_id,value\nUSD.SOFR.OIS.1M,0.035\n",
    )?;
    for (command, headers, id) in [
        ("price", ["trade_id", "currency", "pv"], "USD-SWAP"),
        (
            "risk",
            ["quote_id", "currency", "sensitivity"],
            "USD.SOFR.OIS.1M",
        ),
    ] {
        let output = run(command, dir.path(), &["--date", "2026-11-24"])?;
        let error = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{error}");
        assert!(
            error.contains("calibration for JPY-TONA") && error.contains("trade JPY-SWAP"),
            "{error}"
        );
        let rows = rows(&output, &headers)?;
        assert_eq!(rows.len(), 1);
        assert_eq!((rows[0].0.as_str(), rows[0].1.as_str()), (id, "USD"));
    }
    Ok(())
}

#[test]
fn cli_keeps_other_trades_when_a_published_fixing_is_missing() -> TestResult {
    // 11 月 20 日の TONA は 23 日の休日を挟み、24 日の締め時刻に公表済みになる。
    let dir = workspace_copy()?;
    let path = dir.path().join("trades/swaps.toml");
    let trades = fs::read_to_string(&path)?.replace(
        "convention = \"JPY-TONA-OIS\"\nstart = 2026-11-24\ntenor = \"1M\"",
        "convention = \"JPY-TONA-OIS\"\nstart = 2026-11-20\nmaturity = 2026-12-20",
    );
    fs::write(path, trades)?;
    let output = run("price", dir.path(), &["--date", "2026-11-24"])?;
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    fs::write(dir.path().join("fixings/JPY-TONA.csv"), "date,value\n")?;
    for (command, headers) in [
        ("price", ["trade_id", "currency", "pv"]),
        ("risk", ["quote_id", "currency", "sensitivity"]),
    ] {
        let output = run(command, dir.path(), &["--date", "2026-11-24"])?;
        let error = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{error}");
        assert!(
            error.contains("trade JPY-SWAP") && error.contains("missing published fixing JPY-TONA"),
            "{error}"
        );
        let rows = rows(&output, &headers)?;
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].1, "USD");
    }
    Ok(())
}

fn write_term_fallback_fixture(dir: &Path) -> TestResult {
    let indices = dir.join("refdata/indices.toml");
    let mut definitions = fs::read_to_string(&indices)?;
    definitions.push_str(r#"
[[index]]
id = "JPY-TERM"
calendar = "JPTO"
status = { ceased = { since = 2026-11-19 } }
publication = { business_days = 0, calendar = "JPTO", time = "17:00", zone = "Asia/Tokyo" }
observation_time = { time = "17:00", zone = "Asia/Tokyo" }
kind = { term_rate = { ccy = "JPY", day_count = "ACT365F", tenor = "2BD", value_offset = "0BD", bdc = "following", end_of_month = "ignore" } }

[[fallback]]
index = "JPY-TERM"
trigger = "cessation"
replacement = { index = "JPY-TONA", observation = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } }, annualization = "ACT365F", spread = "10bp" }
applies = [{ definitions = { version = "2021", traded_on_or_after = 2026-11-18 } }]
"#);
    fs::write(indices, definitions)?;

    fs::write(
        dir.join("trades/legacy.toml"),
        r#"
version = 1

[[trade]]
id = "JPY-LEGACY"
margin_set = "JPY-CASH"
book = "JPY-SWAPS"
trade_date = 2026-11-18
fallback = { definitions = "2021" }
fees = [{ at = 2026-11-24T08:00:00Z, amount = "10 JPY", side = "pay" }]

[trade.product]
type = "amended"
from = "2026-11-24T17:00:00[Asia/Tokyo]"
after = { scale = 0.5 }

[trade.product.before]
type = "swap"
legs = [
  { side = "receive", notional = { fixed = "1_000_000 JPY" }, day_count = "ACT365F", coupon = { term = { index = "JPY-TERM" } }, schedule = { explicit = { periods = [{ accrual = { start = 2026-11-19, end = 2026-11-24 }, payment = 2026-11-24, reset = 2026-11-19 }] } } },
  { side = "pay", notional = { fixed = "1_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "0bp" } }, schedule = { explicit = { periods = [{ accrual = { start = 2026-11-19, end = 2026-11-24 }, payment = 2026-11-24 }] } } },
]
"#,
    )?;
    Ok(())
}

#[test]
fn cli_prices_amended_term_coupon_with_fallback_and_isolates_missing_fixing() -> TestResult {
    // 月曜が休日なので 11 月 19 日の TONA は 1 日、20 日は 4 日分。支払いは評価時点と同じで DF は 1。
    let dir = workspace_copy()?;
    write_term_fallback_fixture(dir.path())?;
    let config = dir.path().join("config/pricing.toml");
    fs::write(
        &config,
        fs::read_to_string(&config)?
            .replace("version = 1", "version = 1\nsettlement = \"include\""),
    )?;
    let input = dir.path().join("fixings/JPY-TONA.csv");
    fs::write(&input, "date,value\n2026-11-19,1.2%\n2026-11-20,2.5%\n")?;
    let args = ["--as-of", "2026-11-24T08:00:00Z"];
    let output = run("price", dir.path(), &args)?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(
        stderr.contains("trade JPY-LEGACY: FallbackApplied"),
        "{stderr}"
    );
    let first = 0.012 / 365.0;
    let second = 0.025 * 4.0 / 365.0;
    let expected = 0.5 * 1e6 * (first + second + first * second + 0.001 * 5.0 / 365.0) - 10.0;
    let headers = ["trade_id", "currency", "pv"];
    let priced = rows(&output, &headers)?;
    let legacy = priced
        .iter()
        .find(|(id, _, _)| id == "JPY-LEGACY")
        .ok_or_else(|| io::Error::other("missing legacy trade"))?;
    assert_close!(legacy.2, expected, 1e-12, 1e-6);

    fs::write(&input, "date,value\n2026-11-19,1.2%\n")?;
    let output = run("price", dir.path(), &args)?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.contains("trade JPY-LEGACY") && stderr.contains("missing published fixing JPY-TONA"),
        "{stderr}"
    );
    let priced = rows(&output, &headers)?;
    assert_eq!(
        priced
            .iter()
            .map(|(id, _, _)| id.as_str())
            .collect::<Vec<_>>(),
        ["JPY-SWAP", "USD-SWAP"]
    );
    // 同じ瞬間を決済済みとして扱うと、欠損した fixing を読まずに価格計算が成功する。
    let settings = fs::read_to_string(&config)?
        .replace("settlement = \"include\"", "settlement = \"exclude\"");
    fs::write(&config, settings)?;
    let output = run("price", dir.path(), &args)?;
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let priced = rows(&output, &headers)?;
    assert_eq!(
        priced
            .iter()
            .map(|(id, _, _)| id.as_str())
            .collect::<Vec<_>>(),
        ["JPY-SWAP", "USD-SWAP"]
    );
    Ok(())
}

#[test]
fn cli_uses_only_the_replacement_payment_at_the_cutoff() -> TestResult {
    // 両方の預金が境界ちょうどに満期を迎え、新しい預金の元本だけを受け取る。
    let dir = workspace_copy()?;
    let config = dir.path().join("config/pricing.toml");
    fs::write(
        &config,
        fs::read_to_string(&config)?
            .replace("version = 1", "version = 1\nsettlement = \"include\""),
    )?;
    fs::write(
        dir.path().join("trades/replaced.toml"),
        r#"
version = 1

[[trade]]
id = "REPLACED"
margin_set = "JPY-CASH"
book = "JPY-SWAPS"
trade_date = 2026-11-20

[trade.product]
type = "amended"
from = 2026-11-24T08:00:00Z

[trade.product.before]
type = "deposit"
start = 2026-11-20
maturity = 2026-11-24
notional = "100 JPY"
rate = "0bp"
day_count = "ACT365F"

[trade.product.after.replace]
type = "deposit"
start = 2026-11-20
maturity = 2026-11-24
notional = "200 JPY"
rate = "0bp"
day_count = "ACT365F"
"#,
    )?;
    let output = run("price", dir.path(), &["--as-of", "2026-11-24T08:00:00Z"])?;
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let priced = rows(&output, &["trade_id", "currency", "pv"])?;
    let replaced = priced
        .iter()
        .find(|(id, _, _)| id == "REPLACED")
        .ok_or_else(|| io::Error::other("missing replaced trade"))?;
    assert_close!(replaced.2, 200.0, 1e-12, 1e-6);
    Ok(())
}

#[test]
fn cli_keeps_jpy_rows_when_a_usd_calibration_instrument_cannot_be_prepared() -> TestResult {
    // SOFR の慣行を 1 営業日の lookback にすると、11 月 24 日に始まる校正商品の最初の観測は
    // 23 日の SOFR で、評価の締め時刻には公表済みのはずの値がサンプルに無い。
    // 校正商品の準備の失敗は USD のグループの失敗として残し、独立した JPY の行を出力して
    // 終了コード 1 を返す（以前は校正全体のエラーで、終了コード 2 と出力ゼロになった）。
    let dir = workspace_copy()?;
    let path = dir.path().join("refdata/conventions.toml");
    let sofr = "observation = { days = { business_days = \"USNY\" }, aggregation";
    let conventions = fs::read_to_string(&path)?;
    assert!(conventions.contains(sofr));
    fs::write(
        &path,
        conventions.replace(
            sofr,
            "observation = { days = { business_days = \"USNY\" }, method = { lookback = 1 }, aggregation",
        ),
    )?;
    for (command, headers) in [
        ("price", ["trade_id", "currency", "pv"]),
        ("risk", ["quote_id", "currency", "sensitivity"]),
    ] {
        let output = run(command, dir.path(), &["--date", "2026-11-24"])?;
        let error = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{error}");
        // USD のグループの失敗には、準備できなかった校正商品の quote と欠けた fixing が付く。
        assert!(
            error.contains("calibration for USD-SOFR")
                && error.contains("quote USD.SOFR.OIS.1M")
                && error.contains("missing published fixing USD-SOFR"),
            "{error}"
        );
        let rows = rows(&output, &headers)?;
        assert!(!rows.is_empty(), "{command}: {error}");
        assert!(rows.iter().all(|row| row.1 == "JPY"), "{command}: {rows:?}");
    }
    Ok(())
}

#[test]
fn cli_warns_about_rate_quotes_written_without_a_unit() -> TestResult {
    // 3.5% を単位なしの 3.5 と書くと、値は 350% として使われる。rh はその quote の ID を付けた
    // 警告を標準エラーに出す。単位付きの 45bp と小数の 0.035 は警告しない。
    let dir = workspace_copy()?;
    fs::write(
        dir.path().join("quotes/snapshot.csv"),
        "quote_id,value\nUSD.SOFR.OIS.1M,3.5\nJPY.TONA.OIS.1M,45bp\n",
    )?;
    let output = run("price", dir.path(), &["--date", "2026-11-24"])?;
    let error = String::from_utf8_lossy(&output.stderr);
    assert!(
        error.contains("QuoteUnitSuspected") && error.contains("USD.SOFR.OIS.1M"),
        "{error}"
    );
    assert!(!error.contains("JPY.TONA.OIS.1M"), "{error}");
    // 元のサンプルでは警告が出ない。
    let output = run("price", &sample(), &["--date", "2026-11-24"])?;
    assert_eq!(String::from_utf8_lossy(&output.stderr), "");
    Ok(())
}

#[test]
fn cli_thread_count_keeps_price_and_risk_csv() -> TestResult {
    // サンプルに 1M の取引を 8 件足し、4 スレッドで取引を区間に分けても同じ行と値を出すことを確かめる。
    // 取引ごとの PV は分け方によらず、デルタは区間ごとの和の順だけが変わる。
    let dir = workspace_copy()?;
    let mut trades = String::from("version = 1\n");
    for (i, (ccy, convention, unit, side)) in [
        ("JPY", "JPY-TONA-OIS", 100_000_000.0, "pay"),
        ("USD", "USD-SOFR-OIS", 1_000_000.0, "receive"),
    ]
    .into_iter()
    .cycle()
    .take(8)
    .enumerate()
    {
        let scale = f64::from(u8::try_from(i + 1)?);
        write!(
            trades,
            r#"
[[trade]]
id = "EXTRA-{i}"
margin_set = "{ccy}-CASH"
book = "{ccy}-SWAPS"
trade_date = 2026-11-24

[trade.product]
type = "swap"
convention = "{convention}"
start = 2026-11-24
tenor = "1M"
notional = "{notional} {ccy}"
fixed_rate = "{rate}bp"
fixed_side = "{side}"
"#,
            notional = unit * scale,
            rate = 30 + 10 * i,
        )?;
    }
    fs::write(dir.path().join("trades/extra.toml"), trades)?;
    for (command, headers, relative) in [
        ("price", ["trade_id", "currency", "pv"], 1e-12),
        ("risk", ["quote_id", "currency", "sensitivity"], 1e-10),
    ] {
        let sequential = run(command, dir.path(), &["--date", "2026-11-24"])?;
        let parallel = run(
            command,
            dir.path(),
            &["--date", "2026-11-24", "--threads", "4"],
        )?;
        for output in [&sequential, &parallel] {
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
        }
        let expected = rows(&sequential, &headers)?;
        let actual = rows(&parallel, &headers)?;
        assert_eq!(expected.len(), if command == "price" { 10 } else { 2 });
        assert_eq!(actual.len(), expected.len());
        for ((id, ccy, value), (expected_id, expected_ccy, expected_value)) in
            actual.iter().zip(&expected)
        {
            assert_eq!((id, ccy), (expected_id, expected_ccy));
            let absolute = if ccy == "USD" { 1e-8 } else { 1e-6 };
            assert_close!(*value, *expected_value, relative, absolute);
        }
    }
    // ワーカーの数は 1 以上。0 は引数の誤りとして終了コード 2 になる。
    let zero = run(
        "price",
        dir.path(),
        &["--date", "2026-11-24", "--threads", "0"],
    )?;
    assert_eq!(zero.status.code(), Some(2));
    assert_eq!(zero.stdout, [] as [u8; 0]);
    Ok(())
}

#[test]
fn cli_prices_a_daily_settled_futures_trade_at_its_calibration_quote() -> TestResult {
    // 証券、清算指標、限月、LGM マーク、校正テンプレートを別々のファイルから読み、
    // 同じ限月を同じ quote で約定した取引を rh price で評価する。
    let dir = workspace_copy()?;
    let root = dir.path();
    let append = |name: &str, addition: &str| -> io::Result<()> {
        let path = root.join(name);
        let mut text = fs::read_to_string(&path)?;
        text.push_str(addition);
        fs::write(path, text)
    };
    fs::write(
        root.join("refdata/securities.toml"),
        r#"
version = 1
[[security]]
id = "USD-FUT-NOV"
kind = "listed_future"
currency = "USD"
unit = "contract"
valuation = "contract"
"#,
    )?;
    append(
        "refdata/indices.toml",
        r#"

[[index]]
id = "USD-TERM"
calendar = "USNY"
publication = { business_days = 0, calendar = "USNY", time = "10:00", zone = "America/New_York" }
observation_time = { time = "10:00", zone = "America/New_York" }
kind = { term_rate = { ccy = "USD", day_count = "ACT360", tenor = "1W", value_offset = "0D", bdc = "following", end_of_month = "ignore" } }

[[index]]
id = "USD-FUT-NOV"
calendar = "USNY"
publication = { business_days = 0, calendar = "USNY", time = "10:00", zone = "America/New_York" }
observation_time = { time = "10:00", zone = "America/New_York" }
kind = { rate_futures_settlement = { contract = "USD-FUT-NOV", underlying = "USD-TERM" } }
"#,
    )?;
    append(
        "refdata/conventions.toml",
        r#"

[[convention]]
id = "USD-FUTURES"
[convention.kind.rate_future]
multiplier = 1000.0
settlement_lag = "0D"
settlement_time = { time = "10:00", zone = "America/New_York" }
contracts = { "2026-11" = { settlement_index = "USD-FUT-NOV", final_clearing = 2026-11-30, rate = { term = { fixing = 2026-11-30 } } } }
"#,
    )?;
    let path = root.join("market/definition.toml");
    let market = fs::read_to_string(&path)?.replace(
        "quote = \"USD.SOFR.OIS.1M\"\nkind = \"rate\"\ntarget = \"par_rate\"\ntemplate = { ois = { convention = \"USD-SOFR-OIS\", tenor = \"1M\" } }",
        "quote = \"USD.FUT.NOV\"\nkind = \"futures_price\"\ntarget = \"price\"\ntemplate = { rate_future = { convention = \"USD-FUTURES\", delivery_month = \"2026-11\" } }",
    );
    fs::write(path, market)?;
    append(
        "market/definition.toml",
        r#"

[lgm.USD]
mean_reversion = ["USD.LGM.MR"]
volatility = ["USD.LGM.VOL"]
"#,
    )?;
    append("config/pricing.toml", "\nUSD-TERM = \"USD-SOFR\"\n")?;
    fs::write(
        root.join("quotes/snapshot.csv"),
        "quote_id,value\nJPY.TONA.OIS.1M,45bp\nUSD.FUT.NOV,0.96\nUSD.LGM.MR,0.0\nUSD.LGM.VOL,0.01\n",
    )?;
    fs::write(
        root.join("trades/futures.toml"),
        r#"
version = 1
[[trade]]
id = "USD-FUT-TRADE"
margin_set = "USD-CASH"
book = "USD-SWAPS"
trade_date = 2026-11-24
[trade.product]
type = "rate_future"
settlement_index = "USD-FUT-NOV"
contracts = 2
trade_price = 0.96
day_count = "ACT360"
"#,
    )?;
    let output = run("price", root, &["--as-of", "2026-11-24T08:00:00Z"])?;
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let price = rows(&output, &["trade_id", "currency", "pv"])?
        .into_iter()
        .find(|(id, _, _)| id == "USD-FUT-TRADE")
        .ok_or_else(|| io::Error::other("missing futures PV"))?;
    assert_eq!(price.1, "USD");
    assert_close!(price.2, 0.0_f64, 1e-12_f64, 1e-8_f64);
    Ok(())
}

#[expect(
    clippy::disallowed_macros,
    reason = "リポジトリのサンプルの場所をテストのソースの位置から決める"
)]
fn m1_sample() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m1-rates")
}

#[test]
fn cli_prices_and_risks_the_m1_rates_workspace() -> TestResult {
    // 最初の検証データセットのサンプルを rh price と rh risk で評価する。全 18 取引に PV の行があり、
    // 近似の警告は予想の測度に従う（TONA 担保の TIBOR の IRS は CrossCollateralProjection を出さず、
    // 25bp の CSA の IRS は出す）。ワーカー数を変えても CSV の値は許容誤差で一致する。
    let dir = m1_sample();
    let price = run("price", &dir, &["--date", "2026-11-24"])?;
    assert!(
        price.status.success(),
        "{}",
        String::from_utf8_lossy(&price.stderr)
    );
    let pvs = rows(&price, &["trade_id", "currency", "pv"])?;
    let trades: std::collections::BTreeSet<_> = pvs.iter().map(|(id, _, _)| id.as_str()).collect();
    assert_eq!(trades.len(), 18);
    let stderr = String::from_utf8(price.stderr)?;
    assert!(!stderr.contains("error"), "{stderr}");
    assert!(stderr.contains("trade JPY-IRS-5Y-CASH-25: CrossCollateralProjection"));
    assert!(!stderr.contains("trade JPY-IRS-5Y: CrossCollateralProjection"));
    assert!(stderr.contains("trade USD-OIS-JPY-CASH: CrossCollateralProjection"));

    let risk = run("risk", &dir, &["--date", "2026-11-24"])?;
    assert!(risk.status.success());
    let deltas = rows(&risk, &["quote_id", "currency", "sensitivity"])?;
    for quote in [
        "JPY.TONA.OIS.5Y",
        "JPY.TONA.TURN.2026",
        "JPY.TIBOR6M-TONA.BASIS.5Y",
        "JPY.TIBOR3M-TONA.BASIS.5Y",
        "USD.SOFR.OIS.2Y",
        "USDJPY.SPOT",
        "USDJPY.XCCY.2Y",
        "JPY.LGM.VOL",
    ] {
        assert!(
            deltas
                .iter()
                .any(|(id, _, value)| id == quote && value.abs() > 1e-6),
            "{quote}"
        );
    }

    for (command, values, headers) in [
        ("price", &pvs, ["trade_id", "currency", "pv"]),
        ("risk", &deltas, ["quote_id", "currency", "sensitivity"]),
    ] {
        let parallel = run(command, &dir, &["--date", "2026-11-24", "--threads", "2"])?;
        assert!(parallel.status.success());
        let parallel = rows(&parallel, &headers)?;
        assert_eq!(parallel.len(), values.len());
        for ((id, ccy, value), (parallel_id, parallel_ccy, parallel_value)) in
            values.iter().zip(&parallel)
        {
            assert_eq!((id, ccy), (parallel_id, parallel_ccy));
            assert_close!(*value, *parallel_value, 1e-10_f64, 1e-6_f64);
        }
    }
    Ok(())
}
`````

`````
crates/rhodium-core/Cargo.toml
[package]
name = "rhodium-core"
description = "日付、金額、参照データ、エラーなどの共通基盤"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
publish.workspace = true

[lints]
workspace = true

[dependencies]
chrono.workspace = true
scioxide-ad.workspace = true
scioxide-core.workspace = true
serde.workspace = true
thiserror.workspace = true

[dev-dependencies]
proptest.workspace = true
rhodium-testkit.workspace = true
`````

`````
crates/rhodium-core/src/brand.rs
//! ハンドルを作ったスコープを型で区別する。

use std::marker::PhantomData;

/// 別のスコープの値と揃えられない、不変な lifetime の印。
#[derive(Debug, Clone, Copy)]
pub struct Brand<'id>(PhantomData<fn(&'id ()) -> &'id ()>);

/// 新しい印を高階クロージャの中で使う。
pub fn with_brand<T>(f: impl for<'id> FnOnce(Brand<'id>) -> T) -> T {
    f(Brand(PhantomData))
}
`````

`````
crates/rhodium-core/src/calendar.rs
//! 休日集合、営業日調整、カレンダーを指定した日付の加算。

mod advance;
mod wire;

use std::{collections::BTreeMap, ops::RangeInclusive, sync::Arc};

use crate::{
    Error, Result, Warning,
    date::{Date, DateInterval, Weekday},
    id::CalendarId,
};

pub use advance::{
    AdvanceRule, BusinessDayConvention, EndOfMonthRule, adjust, advance, money_market_period,
};
pub(crate) use advance::{
    adjust_date, advance_date, is_last_business_day, shift_business, shift_calendar_date,
};

/// 基本カレンダーの週末に含む曜日。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct WeekendRule(u8);

impl WeekendRule {
    /// 指定の曜日を休日とする。
    pub fn new(days: impl IntoIterator<Item = Weekday>) -> Self {
        Self(
            days.into_iter()
                .fold(0, |mask, day| mask | weekday_bit(day)),
        )
    }

    /// 曜日が週末に含まれるか。
    #[must_use]
    pub const fn contains(self, day: Weekday) -> bool {
        self.0 & weekday_bit(day) != 0
    }
}

const fn weekday_bit(day: Weekday) -> u8 {
    1 << day.num_days_from_monday()
}

/// 基本カレンダーの休日表の範囲外で使う規則。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Extrapolation {
    /// 欠損として扱う。
    Strict,
    /// 週末の規則で判定し、外挿の警告を返す。
    WeekendsOnly,
}

/// 基本カレンダーの休日表と週末の規則。
#[derive(Debug, Clone)]
pub struct BaseCalendar {
    known: RangeInclusive<Date>,
    weekend: WeekendRule,
    holidays: Arc<[Date]>,
    extrapolation: Extrapolation,
}

impl BaseCalendar {
    /// 両端を含む既知範囲と、その中の休日の表を受け取る。
    ///
    /// # Errors
    ///
    /// 範囲が空か、休日の表に範囲外の日があれば [`Error::InvalidData`]。
    pub fn new(
        known: RangeInclusive<Date>,
        weekend: WeekendRule,
        mut holidays: Vec<Date>,
        extrapolation: Extrapolation,
    ) -> Result<Self> {
        if known.is_empty() || holidays.iter().any(|date| !known.contains(date)) {
            return Err(Error::InvalidData("invalid base calendar coverage".into()));
        }
        holidays.sort_unstable();
        holidays.dedup();
        Ok(Self {
            known,
            weekend,
            holidays: holidays.into(),
            extrapolation,
        })
    }

    fn is_business_day(
        &self,
        id: &CalendarId,
        date: Date,
        warnings: &mut Vec<Warning>,
    ) -> Result<bool> {
        let weekday_is_business = !self.weekend.contains(date.weekday());
        if self.known.contains(&date) {
            return Ok(weekday_is_business && self.holidays.binary_search(&date).is_err());
        }
        match self.extrapolation {
            Extrapolation::Strict => Err(Error::MissingData(format!(
                "calendar {id} does not cover {date}"
            ))),
            Extrapolation::WeekendsOnly => {
                warnings.push(Warning::CalendarExtrapolated {
                    calendar: id.clone(),
                    date,
                });
                Ok(weekday_is_business)
            }
        }
    }
}

/// 基本カレンダーから休日集合を作る式。
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CalendarSpec {
    /// 名前付きの基本カレンダー。
    Base(CalendarId),
    /// いずれかが休日なら休日。二つ以上の式を指定する。
    Union(Vec<Self>),
    /// すべてが休日のときだけ休日。二つ以上の式を指定する。
    Intersection(Vec<Self>),
    /// 休日の集合から日付を除き、追加する。
    Adjusted {
        /// 追加と除去を適用する前の式。
        inner: Box<Self>,
        /// 休日とする日付。
        add: Vec<Date>,
        /// 営業日とする日付。
        remove: Vec<Date>,
    },
}

/// 参照データから解決したカレンダー。既知範囲の営業日判定は O(1)。
#[derive(Debug, Clone)]
pub struct Calendar(Arc<Compiled>);

#[derive(Debug)]
struct Compiled {
    expression: ResolvedSpec,
    known: RangeInclusive<Date>,
    business: Vec<u64>,
}

#[derive(Debug)]
enum ResolvedSpec {
    Base(CalendarId, BaseCalendar),
    Union(Vec<Self>),
    Intersection(Vec<Self>),
    Adjusted {
        inner: Box<Self>,
        add: Vec<Date>,
        remove: Vec<Date>,
    },
}

impl Calendar {
    /// 式の参照を基本カレンダーで解決し、共通の既知範囲を事前計算する。
    ///
    /// # Errors
    ///
    /// 参照先、和・共通部分の要素数、休日の追加・除去が無効なら
    /// [`Error::InvalidData`]。
    pub fn compile(
        spec: &CalendarSpec,
        bases: &BTreeMap<CalendarId, BaseCalendar>,
    ) -> Result<Self> {
        // 既知範囲は、式が参照するすべての基本暦の既知範囲の共通部分。範囲内の営業日を
        // 1 日 1 ビットで前計算し、範囲外は外挿の警告や欠損を返すため照会のたびに式を評価する。
        // 基本暦の既知範囲が交わらなければ範囲は空になり、ビット集合も空になる。
        let mut known = Date::from_ymd(1, 1, 1)?..=Date::from_ymd(9999, 12, 31)?;
        let expression = ResolvedSpec::resolve(spec, bases, &mut known)?;
        let length = (*known.end()).days_since(*known.start()) + 1;
        let length = index(length.max(0))?;
        let mut business = vec![0_u64; length.div_ceil(64)];
        let mut warnings = Vec::new();
        for day_index in 0..length {
            let days = i32::try_from(day_index).map_err(|_| calendar_index_error())?;
            let date = known.start().checked_add_days(days)?;
            if expression.is_business_day(date, &mut warnings)? {
                let word = business
                    .get_mut(day_index / 64)
                    .ok_or_else(calendar_index_error)?;
                *word |= 1 << (day_index % 64);
            }
        }
        Ok(Self(Arc::new(Compiled {
            expression,
            known,
            business,
        })))
    }

    /// 営業日かどうかと、使用した外挿の警告を返す。
    ///
    /// # Errors
    ///
    /// 参照する基本暦の `Strict` 範囲外なら [`Error::MissingData`]。
    pub fn is_business_day(&self, date: Date) -> Result<(bool, Vec<Warning>)> {
        let mut warnings = Vec::new();
        let business = self.lookup(date, &mut warnings)?;
        Ok((business, unique_warnings(warnings)))
    }

    /// 半開区間 [start, end) の営業日数と外挿の警告を返す。
    ///
    /// # Errors
    ///
    /// 逆順なら [`Error::InvalidInput`]。
    /// 区間に参照する基本暦の `Strict` 範囲外の日があれば [`Error::MissingData`]。
    pub fn business_days_between(&self, start: Date, end: Date) -> Result<(u32, Vec<Warning>)> {
        let period = DateInterval::new(start, end)?;
        // 区間を start からの日数で、既知範囲より前 [0, known_start)、既知範囲内
        // [known_start, known_end)、既知範囲より後 [known_end, length) に分ける。
        // 既知範囲内はビットを数え、その外は一日ずつ式を評価して警告と欠損を集める。
        let length = period.end().days_since(period.start());
        let known_start = self.0.known.start().days_since(start).clamp(0, length);
        let known_end = (self.0.known.end().days_since(start) + 1).clamp(known_start, length);
        let mut business = 0;
        if known_start < known_end {
            let first = start
                .checked_add_days(known_start)?
                .days_since(*self.0.known.start());
            business = self
                .0
                .count_bits(index(first)?, index(first + known_end - known_start)?)?;
        }
        let mut warnings = Vec::new();
        for days in (0..known_start).chain(known_end..length) {
            let date = start.checked_add_days(days)?;
            business += u32::from(self.0.expression.is_business_day(date, &mut warnings)?);
        }
        Ok((business, unique_warnings(warnings)))
    }

    pub(crate) fn base_id(&self) -> Option<&CalendarId> {
        match &self.0.expression {
            ResolvedSpec::Base(id, _) => Some(id),
            _ => None,
        }
    }

    pub(crate) fn lookup(&self, date: Date, warnings: &mut Vec<Warning>) -> Result<bool> {
        if self.0.known.contains(&date) {
            let bit = index(date.days_since(*self.0.known.start()))?;
            let word = self
                .0
                .business
                .get(bit / 64)
                .ok_or_else(calendar_index_error)?;
            return Ok(word & (1 << (bit % 64)) != 0);
        }
        self.0.expression.is_business_day(date, warnings)
    }
}

impl Compiled {
    // ビット添字の半開区間 [first, after) で立っているビットを数える。呼び出し側は first < after を
    // 保証する。先頭の語は first_bit より下のビットを、末尾の語は last_bit より上のビットを落とす。
    fn count_bits(&self, first: usize, after: usize) -> Result<u32> {
        let first_word = first / 64;
        let last_word = (after - 1) / 64;
        let first_bit = first % 64;
        let last_bit = (after - 1) % 64;
        let words = self
            .business
            .get(first_word..=last_word)
            .ok_or_else(calendar_index_error)?;
        Ok(words
            .iter()
            .enumerate()
            .map(|(index, &word)| {
                let lower = if index == 0 {
                    u64::MAX << first_bit
                } else {
                    u64::MAX
                };
                let upper = if index == last_word - first_word {
                    u64::MAX >> (63 - last_bit)
                } else {
                    u64::MAX
                };
                (word & lower & upper).count_ones()
            })
            .sum())
    }
}

impl ResolvedSpec {
    fn resolve(
        spec: &CalendarSpec,
        bases: &BTreeMap<CalendarId, BaseCalendar>,
        known: &mut RangeInclusive<Date>,
    ) -> Result<Self> {
        match spec {
            CalendarSpec::Base(id) => {
                let base = bases
                    .get(id)
                    .ok_or_else(|| Error::InvalidData(format!("unknown calendar {id}")))?;
                *known = (*known.start()).max(*base.known.start())
                    ..=(*known.end()).min(*base.known.end());
                Ok(Self::Base(id.clone(), base.clone()))
            }
            CalendarSpec::Union(items) => Self::resolve_group(items, bases, known).map(Self::Union),
            CalendarSpec::Intersection(items) => {
                Self::resolve_group(items, bases, known).map(Self::Intersection)
            }
            CalendarSpec::Adjusted { inner, add, remove } => {
                let mut add = add.clone();
                let mut remove = remove.clone();
                add.sort_unstable();
                add.dedup();
                remove.sort_unstable();
                remove.dedup();
                if add.iter().any(|date| remove.binary_search(date).is_ok()) {
                    return Err(Error::InvalidData(
                        "calendar adds and removes the same holiday".into(),
                    ));
                }
                Ok(Self::Adjusted {
                    inner: Box::new(Self::resolve(inner, bases, known)?),
                    add,
                    remove,
                })
            }
        }
    }

    // 和と共通部分は二つ以上の式を要求し、各要素の参照で既知範囲を共通部分に絞る。
    fn resolve_group(
        items: &[CalendarSpec],
        bases: &BTreeMap<CalendarId, BaseCalendar>,
        known: &mut RangeInclusive<Date>,
    ) -> Result<Vec<Self>> {
        validate_group_len(items.len())?;
        items
            .iter()
            .map(|item| Self::resolve(item, bases, known))
            .collect()
    }

    fn is_business_day(&self, date: Date, warnings: &mut Vec<Warning>) -> Result<bool> {
        match self {
            Self::Base(id, base) => base.is_business_day(id, date, warnings),
            // 短絡させず、参照したすべての基本暦の欠損と外挿を伝える。
            Self::Union(items) => items.iter().try_fold(true, |result, item| {
                Ok(result & item.is_business_day(date, warnings)?)
            }),
            Self::Intersection(items) => items.iter().try_fold(false, |result, item| {
                Ok(result | item.is_business_day(date, warnings)?)
            }),
            Self::Adjusted { inner, add, remove } => {
                let business = inner.is_business_day(date, warnings)?;
                Ok(add.binary_search(&date).is_err()
                    && (remove.binary_search(&date).is_ok() || business))
            }
        }
    }
}

pub(crate) fn unique_warnings(mut warnings: Vec<Warning>) -> Vec<Warning> {
    warnings.sort_unstable();
    warnings.dedup();
    warnings
}

fn index(days: i32) -> Result<usize> {
    usize::try_from(days).map_err(|_| calendar_index_error())
}

fn calendar_index_error() -> Error {
    Error::Internal("calendar bit index outside the compiled range".into())
}

fn validate_group_len(length: usize) -> Result<()> {
    if length < 2 {
        return Err(Error::InvalidData(
            "calendar group needs at least two expressions".into(),
        ));
    }
    Ok(())
}
`````

`````
crates/rhodium-core/src/calendar/advance.rs
use chrono::Months;

use super::{Calendar, unique_warnings};
use crate::{
    Error, Result, Warning,
    date::{Date, DateInterval, Weekday, date_error},
    tenor::{DateOffset, Tenor, TenorUnit},
};

/// 営業日でない日を動かす規則。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BusinessDayConvention {
    /// 調整しない。
    Unadjusted,
    /// 次の営業日へ進める。
    Following,
    /// 次の営業日が翌月なら前の営業日に戻す。
    ModifiedFollowing,
    /// 前の営業日へ戻す。
    Preceding,
    /// 前の営業日が前月なら次の営業日に進める。
    ModifiedPreceding,
    /// 日曜か月曜なら先へ、それ以外なら前へ探す。
    Nearest,
}

/// 月または年を加算するときの月末規則。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EndOfMonthRule {
    /// 日を保ち、存在しない日だけ月末に詰める。
    Ignore,
    /// 起点が月末の暦日なら、結果も月末の暦日にする。
    LastCalendarDay,
    /// 起点が月末の営業日なら、結果も月末の営業日にする。
    LastBusinessDay,
}

/// 日付の加算に使う月末規則、営業日調整、解決済みのカレンダー。
#[derive(Debug, Clone, Copy)]
pub struct AdvanceRule<'a> {
    /// 月または年の加算時の月末規則。
    pub end_of_month: EndOfMonthRule,
    /// 暦日で加算した結果に適用する調整。
    pub bdc: BusinessDayConvention,
    /// 営業日の照会と調整に使うカレンダー。
    pub calendar: &'a Calendar,
}

/// 日付を営業日調整し、使用した外挿の警告を返す。
///
/// # Errors
///
/// 検索する日付が値域を超えるか、必要なカレンダーの範囲を外れた場合。
pub fn adjust(
    date: Date,
    bdc: BusinessDayConvention,
    calendar: &Calendar,
) -> Result<(Date, Vec<Warning>)> {
    let mut warnings = Vec::new();
    let adjusted = adjust_date(date, bdc, calendar, &mut warnings)?;
    Ok((adjusted, unique_warnings(warnings)))
}

/// 符号付きのオフセットで日付を進め、調整と外挿の警告を返す。
///
/// # Errors
///
/// 日付が値域を超えるか、規則に従う営業日を決められない場合。
pub fn advance(
    date: Date,
    offset: DateOffset,
    rule: &AdvanceRule<'_>,
) -> Result<(Date, Vec<Warning>)> {
    let mut warnings = Vec::new();
    let advanced = advance_date(
        date,
        i32::from(offset.amount()),
        offset.unit(),
        rule,
        &mut warnings,
    )?;
    Ok((advanced, unique_warnings(warnings)))
}

/// 短期市場テナーまたは期間の長さを、約定日・スポット日から半開区間にする。
///
/// # Errors
///
/// 日付かカレンダーの範囲を超えるか、テナーに従う区間を作れない場合。
pub fn money_market_period(
    trade_date: Date,
    tenor: Tenor,
    spot_lag: DateOffset,
    rule: &AdvanceRule<'_>,
) -> Result<(DateInterval, Vec<Warning>)> {
    let mut warnings = Vec::new();
    let shift =
        |date, amount, warnings: &mut Vec<_>| shift_business(date, amount, rule.calendar, warnings);
    let spot = |warnings: &mut Vec<_>| {
        advance_date(
            trade_date,
            i32::from(spot_lag.amount()),
            spot_lag.unit(),
            rule,
            warnings,
        )
    };
    let (start, end) = match tenor {
        Tenor::Overnight => (trade_date, shift(trade_date, 1, &mut warnings)?),
        Tenor::TomorrowNext => {
            let tomorrow = shift(trade_date, 1, &mut warnings)?;
            (tomorrow, shift(tomorrow, 1, &mut warnings)?)
        }
        Tenor::SpotNext => {
            let start = spot(&mut warnings)?;
            (start, shift(start, 1, &mut warnings)?)
        }
        Tenor::Length(length) => {
            let start = spot(&mut warnings)?;
            let end = advance_date(
                start,
                i32::from(length.amount()),
                length.unit(),
                rule,
                &mut warnings,
            )?;
            (start, end)
        }
    };
    Ok((DateInterval::new(start, end)?, unique_warnings(warnings)))
}

pub(crate) fn advance_date(
    date: Date,
    amount: i32,
    unit: TenorUnit,
    rule: &AdvanceRule<'_>,
    warnings: &mut Vec<Warning>,
) -> Result<Date> {
    let shifted = match unit {
        TenorUnit::BusinessDay => return shift_business(date, amount, rule.calendar, warnings),
        TenorUnit::Month | TenorUnit::Year
            if rule.end_of_month == EndOfMonthRule::LastBusinessDay
                && is_last_business_day(date, rule.calendar, warnings)? =>
        {
            let target = shift_calendar_date(date, amount, unit, false)?;
            search(
                target.month_end()?,
                -1,
                Some(target.month_start()?),
                rule.calendar,
                warnings,
            )?
            .ok_or_else(|| Error::InvalidInput("target month has no business day".into()))?
        }
        _ => shift_calendar_date(
            date,
            amount,
            unit,
            rule.end_of_month == EndOfMonthRule::LastCalendarDay,
        )?,
    };
    adjust_date(shifted, rule.bdc, rule.calendar, warnings)
}

pub(crate) fn shift_calendar_date(
    date: Date,
    amount: i32,
    unit: TenorUnit,
    preserve_month_end: bool,
) -> Result<Date> {
    let months = match unit {
        TenorUnit::Day => return date.checked_add_days(amount),
        TenorUnit::Week => {
            return date.checked_add_days(amount.checked_mul(7).ok_or_else(date_error)?);
        }
        TenorUnit::Month => amount,
        TenorUnit::Year => amount.checked_mul(12).ok_or_else(date_error)?,
        TenorUnit::BusinessDay => {
            return Err(Error::InvalidInput(
                "business day offsets require a calendar".into(),
            ));
        }
    };
    // chrono は同じ日が無ければ対象月の月末に詰める。
    let magnitude = Months::new(months.unsigned_abs());
    let target = if months.is_negative() {
        date.naive().checked_sub_months(magnitude)
    } else {
        date.naive().checked_add_months(magnitude)
    };
    let target = Date::from_naive(target.ok_or_else(date_error)?)?;
    if preserve_month_end && date == date.month_end()? {
        target.month_end()
    } else {
        Ok(target)
    }
}

pub(crate) fn is_last_business_day(
    date: Date,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<bool> {
    if !calendar.lookup(date, warnings)? {
        return Ok(false);
    }
    let last = date.month_end()?;
    if date == last {
        return Ok(true);
    }
    Ok(search(date.checked_add_days(1)?, 1, Some(last), calendar, warnings)?.is_none())
}

pub(crate) fn shift_business(
    mut date: Date,
    amount: i32,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<Date> {
    if amount == 0 {
        return search(date, 1, None, calendar, warnings)?.ok_or_else(date_error);
    }
    let direction = amount.signum();
    for _ in 0..amount.unsigned_abs() {
        date = search(
            date.checked_add_days(direction)?,
            direction,
            None,
            calendar,
            warnings,
        )?
        .ok_or_else(date_error)?;
    }
    Ok(date)
}

pub(crate) fn adjust_date(
    date: Date,
    bdc: BusinessDayConvention,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<Date> {
    use BusinessDayConvention::{
        Following, ModifiedFollowing, ModifiedPreceding, Nearest, Preceding, Unadjusted,
    };
    let (direction, boundary) = match bdc {
        Unadjusted => return Ok(date),
        Following => (1, None),
        ModifiedFollowing => (1, Some(date.month_end()?)),
        Preceding => (-1, None),
        ModifiedPreceding => (-1, Some(date.month_start()?)),
        Nearest => {
            let direction = if matches!(date.weekday(), Weekday::Sun | Weekday::Mon) {
                1
            } else {
                -1
            };
            (direction, None)
        }
    };
    if let Some(adjusted) = search(date, direction, boundary, calendar, warnings)? {
        return Ok(adjusted);
    }
    // 月内の元の向きに営業日がない場合、反対側に向き直す。
    search(
        date.checked_add_days(-direction)?,
        -direction,
        None,
        calendar,
        warnings,
    )?
    .ok_or_else(date_error)
}

fn search(
    mut date: Date,
    step: i32,
    boundary: Option<Date>,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<Option<Date>> {
    loop {
        if calendar.lookup(date, warnings)? {
            return Ok(Some(date));
        }
        if boundary == Some(date) {
            return Ok(None);
        }
        date = date.checked_add_days(step)?;
    }
}
`````

`````
crates/rhodium-core/src/calendar/wire.rs
use std::fmt;

use serde::{
    Deserialize, Deserializer, Serialize, Serializer,
    de::{Error as _, IgnoredAny, MapAccess, Visitor, value::MapAccessDeserializer},
};

use super::{CalendarSpec, validate_group_len};
use crate::date::Date;

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
enum CalendarOutput<'a> {
    Union(&'a [CalendarSpec]),
    Intersection(&'a [CalendarSpec]),
    Adjusted {
        inner: &'a CalendarSpec,
        add: &'a [Date],
        remove: &'a [Date],
    },
}

impl Serialize for CalendarSpec {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Base(id) => id.serialize(serializer),
            Self::Union(items) => CalendarOutput::Union(items).serialize(serializer),
            Self::Intersection(items) => CalendarOutput::Intersection(items).serialize(serializer),
            Self::Adjusted { inner, add, remove } => {
                CalendarOutput::Adjusted { inner, add, remove }.serialize(serializer)
            }
        }
    }
}

#[derive(Deserialize)]
#[serde(field_identifier, rename_all = "snake_case")]
enum CalendarKey {
    Union,
    Intersection,
    Adjusted,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Adjustment {
    inner: Box<CalendarSpec>,
    add: Vec<Date>,
    remove: Vec<Date>,
}

impl<'de> Deserialize<'de> for CalendarSpec {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(CalendarVisitor)
    }
}

struct CalendarVisitor;

impl<'de> Visitor<'de> for CalendarVisitor {
    type Value = CalendarSpec;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a calendar ID or a table with one calendar operation")
    }

    fn visit_str<E: serde::de::Error>(self, id: &str) -> Result<Self::Value, E> {
        id.parse().map(CalendarSpec::Base).map_err(E::custom)
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
        let key = map
            .next_key::<CalendarKey>()?
            .ok_or_else(|| A::Error::custom("calendar table needs one operation"))?;
        let spec = match key {
            CalendarKey::Union => CalendarSpec::Union(read_group(&mut map)?),
            CalendarKey::Intersection => CalendarSpec::Intersection(read_group(&mut map)?),
            CalendarKey::Adjusted => {
                let Adjustment { inner, add, remove } = map.next_value_seed(AdjustmentMap)?;
                CalendarSpec::Adjusted { inner, add, remove }
            }
        };
        if map.next_key::<IgnoredAny>()?.is_some() {
            return Err(A::Error::custom(
                "calendar table must contain exactly one operation",
            ));
        }
        Ok(spec)
    }
}

struct AdjustmentMap;

impl<'de> serde::de::DeserializeSeed<'de> for AdjustmentMap {
    type Value = Adjustment;

    fn deserialize<D: Deserializer<'de>>(self, deserializer: D) -> Result<Self::Value, D::Error> {
        deserializer.deserialize_map(self)
    }
}

impl<'de> Visitor<'de> for AdjustmentMap {
    type Value = Adjustment;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a table with inner, add, and remove")
    }

    fn visit_map<A: MapAccess<'de>>(self, map: A) -> Result<Self::Value, A::Error> {
        Adjustment::deserialize(MapAccessDeserializer::new(map))
    }
}

fn read_group<'de, A: MapAccess<'de>>(map: &mut A) -> Result<Vec<CalendarSpec>, A::Error> {
    let items: Vec<_> = map.next_value()?;
    validate_group_len(items.len()).map_err(A::Error::custom)?;
    Ok(items)
}
`````

`````
crates/rhodium-core/src/currency.rs
//! 通貨のコードと、受け渡す資産の種類。

use std::{fmt, str::FromStr};

use crate::{
    Error, Result,
    id::{CommodityId, SecurityId},
};

/// ASCII の英大文字 3 字から成る通貨コード。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Currency([u8; 3]);

impl Currency {
    /// 日本円。
    pub const JPY: Self = Self(*b"JPY");
    /// 米ドル。
    pub const USD: Self = Self(*b"USD");

    /// 通貨コードの文字列。
    #[must_use]
    pub fn as_str(&self) -> &str {
        // 構築時に英大文字だけを受け付ける。
        std::str::from_utf8(&self.0).unwrap_or_default()
    }
}

impl FromStr for Currency {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        let code: [u8; 3] = text.as_bytes().try_into().map_err(|_| currency_error())?;
        if !code.iter().all(u8::is_ascii_uppercase) {
            return Err(currency_error());
        }
        Ok(Self(code))
    }
}

impl fmt::Display for Currency {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

serde_text!(Currency);

fn currency_error() -> Error {
    Error::InvalidData("currency code must be three uppercase ASCII letters".into())
}

/// 受け渡しや保有の対象。
#[derive(
    Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "snake_case")]
pub enum Asset {
    /// 現金。
    Cash(Currency),
    /// 証券や株価指数の単位。
    Security(SecurityId),
    /// 銘柄と受渡場所の決まったコモディティ。
    Commodity(CommodityId),
}
`````

`````
crates/rhodium-core/src/date.rs
//! 営業日や利息の期間を表す現地の日付。

use std::{fmt, str::FromStr};

use chrono::{Datelike, Months, NaiveDate, TimeDelta};

use crate::{Error, Result};

pub use chrono::Weekday;

/// 先発グレゴリオ暦の現地の日付。年は 1〜9999。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Date(NaiveDate);

impl Date {
    /// 年月日から日付を作る。
    ///
    /// # Errors
    ///
    /// 年が値域外か、存在しない日付なら [`Error::InvalidInput`]。
    pub fn from_ymd(year: i32, month: u8, day: u8) -> Result<Self> {
        Self::from_naive(
            NaiveDate::from_ymd_opt(year, u32::from(month), u32::from(day))
                .ok_or_else(date_error)?,
        )
    }

    pub(crate) fn from_naive(date: NaiveDate) -> Result<Self> {
        if (1..=9999).contains(&date.year()) {
            Ok(Self::from_valid_naive(date))
        } else {
            Err(date_error())
        }
    }

    pub(crate) const fn from_valid_naive(date: NaiveDate) -> Self {
        Self(date)
    }

    pub(crate) const fn naive(self) -> NaiveDate {
        self.0
    }

    /// 年月日。
    #[must_use]
    #[allow(clippy::cast_possible_truncation, reason = "月は 12、日は 31 以下")]
    pub fn ymd(self) -> (i32, u8, u8) {
        (self.0.year(), self.0.month() as u8, self.0.day() as u8)
    }

    /// 現地の曜日。
    #[must_use]
    pub fn weekday(self) -> Weekday {
        self.0.weekday()
    }

    /// 符号付きの暦日数を足す。
    ///
    /// # Errors
    ///
    /// 結果が日付の値域を出ると [`Error::InvalidInput`]。
    pub fn checked_add_days(self, days: i32) -> Result<Self> {
        Self::from_naive(
            self.0
                .checked_add_signed(TimeDelta::days(i64::from(days)))
                .ok_or_else(date_error)?,
        )
    }

    /// 二つの日付の間の符号付き暦日数（self − earlier）。
    #[must_use]
    #[allow(
        clippy::cast_possible_truncation,
        reason = "表せる全期間でも i32 の暦日数に収まる"
    )]
    pub fn days_since(self, earlier: Self) -> i32 {
        self.0.signed_duration_since(earlier.0).num_days() as i32
    }

    pub(crate) fn month_start(self) -> Result<Self> {
        Ok(Self(self.0.with_day(1).ok_or_else(date_error)?))
    }

    pub(crate) fn month_end(self) -> Result<Self> {
        // 9999 年 12 月も chrono の翌月初から求められる。
        let last = self
            .month_start()?
            .0
            .checked_add_months(Months::new(1))
            .and_then(|next| next.pred_opt())
            .ok_or_else(date_error)?;
        Ok(Self(last))
    }
}

impl FromStr for Date {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        let bytes = text.as_bytes();
        if bytes.len() != 10
            || bytes.iter().enumerate().any(|(index, byte)| match index {
                4 | 7 => *byte != b'-',
                _ => !byte.is_ascii_digit(),
            })
        {
            return Err(date_error());
        }
        Self::from_naive(NaiveDate::parse_from_str(text, "%Y-%m-%d").map_err(|_| date_error())?)
    }
}

impl fmt::Display for Date {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

serde_text!(Date);

/// 日付の半開区間 [開始日, 終了日)。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(try_from = "DateIntervalInput")]
pub struct DateInterval {
    start: Date,
    end: Date,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct DateIntervalInput {
    start: Date,
    end: Date,
}

impl TryFrom<DateIntervalInput> for DateInterval {
    type Error = Error;

    fn try_from(input: DateIntervalInput) -> Result<Self> {
        Self::new(input.start, input.end)
    }
}

impl DateInterval {
    /// 開始日と終了日から区間を作る。同日は空の区間になる。
    ///
    /// # Errors
    ///
    /// 開始日が終了日より後なら [`Error::InvalidInput`]。
    pub fn new(start: Date, end: Date) -> Result<Self> {
        if start > end {
            return Err(Error::InvalidInput("date interval is reversed".into()));
        }
        Ok(Self { start, end })
    }

    /// 区間に含む開始日。
    #[must_use]
    pub const fn start(self) -> Date {
        self.start
    }

    /// 区間に含まない終了日。
    #[must_use]
    pub const fn end(self) -> Date {
        self.end
    }

    /// 空の区間か。
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.start == self.end
    }
}

pub(crate) fn date_error() -> Error {
    Error::InvalidInput("invalid date or date outside years 0001–9999".into())
}
`````

`````
crates/rhodium-core/src/day_count.rs
//! 利息の期間を日数基準に従って年数に換算する。

use std::{fmt, num::NonZeroU8, str::FromStr};

use chrono::{Datelike, NaiveDate};

use crate::{
    Error, Result, Warning,
    calendar::Calendar,
    date::{Date, DateInterval},
    id::CalendarId,
    quantity::YearFraction,
    tenor::Frequency,
};

/// 日数基準。BUS/252 では営業日を数える基本カレンダーを指定する。
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum DayCount {
    /// ACT/360。
    Act360,
    /// ACT/365 Fixed。
    Act365F,
    /// ACT/365 Leap。
    Act365L,
    /// ACT/364。
    Act364,
    /// ACT/ACT ISDA。
    ActActIsda,
    /// ACT/ACT ICMA。
    ActActIcma,
    /// ACT/ACT AFB。
    ActActAfb,
    /// NL/365。
    Nl365,
    /// 30/360 Bond Basis。
    Thirty360BondBasis,
    /// 30/360 US。
    Thirty360Us,
    /// 30E/360。
    ThirtyE360,
    /// 30E/360 ISDA。
    ThirtyE360Isda,
    /// 30E+/360。
    ThirtyEPlus360,
    /// BUS/252。
    Bus252 {
        /// 営業日を数える基本カレンダー。
        calendar: CalendarId,
    },
    /// 1/1。
    OneOne,
}

impl FromStr for DayCount {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        Ok(match text {
            "ACT360" => Self::Act360,
            "ACT365F" => Self::Act365F,
            "ACT365L" => Self::Act365L,
            "ACT364" => Self::Act364,
            "ACTACT_ISDA" => Self::ActActIsda,
            "ACTACT_ICMA" => Self::ActActIcma,
            "ACTACT_AFB" => Self::ActActAfb,
            "NL365" => Self::Nl365,
            "THIRTY_360_BOND" => Self::Thirty360BondBasis,
            "THIRTY_360_US" => Self::Thirty360Us,
            "THIRTY_E_360" => Self::ThirtyE360,
            "THIRTY_E_360_ISDA" => Self::ThirtyE360Isda,
            "THIRTY_E_PLUS_360" => Self::ThirtyEPlus360,
            "ONE_ONE" => Self::OneOne,
            _ => {
                let calendar = text
                    .strip_prefix("BUS252(")
                    .and_then(|id| id.strip_suffix(')'))
                    .ok_or_else(|| Error::InvalidData("invalid day count".into()))?;
                Self::Bus252 {
                    calendar: calendar.parse()?,
                }
            }
        })
    }
}

impl fmt::Display for DayCount {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Act360 => "ACT360",
            Self::Act365F => "ACT365F",
            Self::Act365L => "ACT365L",
            Self::Act364 => "ACT364",
            Self::ActActIsda => "ACTACT_ISDA",
            Self::ActActIcma => "ACTACT_ICMA",
            Self::ActActAfb => "ACTACT_AFB",
            Self::Nl365 => "NL365",
            Self::Thirty360BondBasis => "THIRTY_360_BOND",
            Self::Thirty360Us => "THIRTY_360_US",
            Self::ThirtyE360 => "THIRTY_E_360",
            Self::ThirtyE360Isda => "THIRTY_E_360_ISDA",
            Self::ThirtyEPlus360 => "THIRTY_E_PLUS_360",
            Self::Bus252 { calendar } => return write!(f, "BUS252({calendar})"),
            Self::OneOne => "ONE_ONE",
        })
    }
}

serde_text!(DayCount);

/// 利払いの期間と、それを覆う通常の期間および頻度。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PeriodFrame {
    period: DateInterval,
    regular: Vec<DateInterval>,
    frequency: Frequency,
}

impl PeriodFrame {
    /// 通常の期間の連続した列から利払いの枠を作る。
    ///
    /// 年あたりの回数 n を持つ頻度（ACT/ACT ICMA と ACT/365L が使えるもの）では、各通常期間が
    /// 頻度のおよそ一刻みであることも確かめる。日数 d の期間について、最も近い刻みの数
    /// round(d·n/365.25) が 1、つまり 0.5 ≤ d·n/365.25 < 1.5 であればよい。
    /// ロール（毎月第 3 水曜の 28 日の期間など）と営業日調整による数日のずれは通し、
    /// 二刻み以上や半刻み未満の参照期間は拒否する。半年払いに一年の参照期間を渡すと、
    /// 9 か月の期間の ICMA の係数が 3/4 ではなく 273/(2·365) ≈ 0.374 になるためである。
    /// 年あたりの回数を持たない頻度（週、日、5M など）の枠は日数計算が使えないので、長さを問わない。
    ///
    /// # Errors
    ///
    /// 通常の期間が空、連続していない、利払いの期間を覆っていない、または頻度の一刻みに
    /// 当たらない場合は [`Error::InvalidInput`]。
    pub fn new(
        period: DateInterval,
        regular: Vec<DateInterval>,
        frequency: Frequency,
    ) -> Result<Self> {
        let Some((first, last)) = regular.first().zip(regular.last()) else {
            return Err(Error::InvalidInput(
                "a period frame requires at least one regular period".into(),
            ));
        };
        if first.start() > period.start()
            || last.end() < period.end()
            || regular.iter().any(|item| item.is_empty())
            || regular
                .windows(2)
                .any(|pair| pair[0].end() != pair[1].start())
        {
            return Err(Error::InvalidInput(
                "regular periods must be contiguous and cover the accrual period".into(),
            ));
        }
        if let Some(per_year) = frequency.periods_per_year()
            && let Some(irregular) = regular.iter().find(|item| !is_one_step(**item, per_year))
        {
            return Err(Error::InvalidInput(format!(
                "regular period {irregular:?} is not one step of {frequency}"
            )));
        }
        Ok(Self {
            period,
            regular,
            frequency,
        })
    }

    /// 利払いの期間。
    #[must_use]
    pub const fn period(&self) -> DateInterval {
        self.period
    }

    /// 利払いの期間を覆う、通常の期間の列。
    #[must_use]
    pub fn regular(&self) -> &[DateInterval] {
        &self.regular
    }

    /// 利払いの頻度。
    #[must_use]
    pub const fn frequency(&self) -> Frequency {
        self.frequency
    }

    fn periods_per_year(&self) -> Result<NonZeroU8> {
        self.frequency.periods_per_year().ok_or_else(|| {
            Error::InvalidInput("day count requires a frequency with periods per year".into())
        })
    }
}

// 期間の長さに最も近い刻みの数が 1 か。平均の年の長さ 365.25 日で刻みの日数を見積もる。
fn is_one_step(interval: DateInterval, per_year: NonZeroU8) -> bool {
    let steps =
        f64::from(interval.end().days_since(interval.start())) * f64::from(per_year.get()) / 365.25;
    (0.5..1.5).contains(&steps)
}

/// 日程が月末の慣行に従っているか。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MonthEndSchedule {
    /// 月末にそろえる。
    Aligned,
    /// 月末にそろえない。
    #[default]
    NotAligned,
}

/// 日数基準が参照する利払いの事実。
#[derive(Debug, Clone, Copy, Default)]
pub struct AccrualContext<'a> {
    /// ACT/ACT ICMA と ACT/365L が使う期間と頻度。
    pub frame: Option<&'a PeriodFrame>,
    /// 30/360 US が使う月末の慣行。
    pub month_end: MonthEndSchedule,
    /// 30E/360 ISDA が使う満期日。
    pub maturity: Option<Date>,
    /// BUS/252 が使う、日数基準の ID と同じ基本カレンダー。
    pub calendar: Option<&'a Calendar>,
}

impl AccrualContext<'_> {
    fn covering_frame(&self, start: Date, end: Date) -> Result<&PeriodFrame> {
        self.frame
            .filter(|frame| frame.period.start() <= start && end <= frame.period.end())
            .ok_or_else(|| {
                Error::InvalidInput("day count requires a frame covering the accrual period".into())
            })
    }
}

/// 半開区間の日数係数と、計算で使用したカレンダー外挿の警告を返す。
///
/// # Errors
///
/// 区間が逆順か、非空の区間で日数基準に必要な文脈が揃わない場合は
/// [`Error::InvalidInput`]。BUS/252 の基本カレンダーが期間を覆えない場合は
/// [`Error::MissingData`]。
pub fn year_fraction(
    dc: &DayCount,
    start: Date,
    end: Date,
    ctx: &AccrualContext<'_>,
) -> Result<(YearFraction, Vec<Warning>)> {
    let period = DateInterval::new(start, end)?;
    let mut warnings = Vec::new();
    let value = if period.is_empty() {
        0.0
    } else {
        let days = f64::from(end.days_since(start));
        match dc {
            DayCount::Act360 => days / 360.0,
            DayCount::Act365F => days / 365.0,
            DayCount::Act365L => days / act_365l_denominator(ctx.covering_frame(start, end)?)?,
            DayCount::Act364 => days / 364.0,
            DayCount::ActActIsda => isda(start, end)?,
            DayCount::ActActIcma => icma(start, end, ctx.covering_frame(start, end)?)?,
            DayCount::ActActAfb => afb(start, end)?,
            DayCount::Nl365 => {
                let excluded: u32 = leap_dates(start, end)
                    .filter(|date| start.naive() < *date && *date <= end.naive())
                    .map(|_| 1)
                    .sum();
                (days - f64::from(excluded)) / 365.0
            }
            // d₁ は 31 を 30 にし、d₁ が 30 か 31 だったときだけ d₂ の 31 を 30 にする。
            DayCount::Thirty360BondBasis => {
                let first = start.ymd().2.min(30);
                let second = end.ymd().2;
                thirty(
                    start,
                    end,
                    first,
                    if first == 30 { second.min(30) } else { second },
                )
            }
            DayCount::Thirty360Us => thirty_us(start, end, ctx.month_end)?,
            DayCount::ThirtyE360 => thirty(start, end, start.ymd().2.min(30), end.ymd().2.min(30)),
            DayCount::ThirtyE360Isda => {
                let maturity = ctx
                    .maturity
                    .ok_or_else(|| Error::InvalidInput("30E/360 ISDA requires maturity".into()))?;
                thirty_e_isda(start, end, maturity)?
            }
            // 31 日を翌月 1 日に移すと、30 日基準の分子では +31 のままになる。
            DayCount::ThirtyEPlus360 => thirty(start, end, start.ymd().2.min(30), end.ymd().2),
            DayCount::Bus252 { calendar } => {
                let base = ctx
                    .calendar
                    .filter(|base| base.base_id() == Some(calendar))
                    .ok_or_else(|| {
                        Error::InvalidInput(format!("BUS/252 requires calendar {calendar}"))
                    })?;
                let (count, notices) = base.business_days_between(start, end)?;
                warnings = notices;
                f64::from(count) / 252.0
            }
            DayCount::OneOne => 1.0,
        }
    };
    Ok((YearFraction::new(value, dc.clone())?, warnings))
}

fn days_in_year(leap: bool) -> f64 {
    if leap { 366.0 } else { 365.0 }
}

fn leap_dates(start: Date, end: Date) -> impl Iterator<Item = NaiveDate> {
    (start.naive().year()..=end.naive().year())
        .filter_map(|year| NaiveDate::from_ymd_opt(year, 2, 29))
}

fn act_365l_denominator(frame: &PeriodFrame) -> Result<f64> {
    let start = frame.period.start();
    let end = frame.period.end();
    let leap = if frame.periods_per_year()?.get() == 1 {
        // 年払いは利払いの期間 (start, end] の閏日を使う。
        leap_dates(start, end).any(|date| start.naive() < date && date <= end.naive())
    } else {
        end.naive().leap_year()
    };
    Ok(days_in_year(leap))
}

fn isda(start: Date, end: Date) -> Result<f64> {
    let first_year = start.naive().year();
    let last_year = end.naive().year();
    let first_denominator = days_in_year(start.naive().leap_year());
    if first_year == last_year {
        return Ok(f64::from(end.days_since(start)) / first_denominator);
    }
    // 間にある完全な暦年の寄与はどの年も 1。
    let first_days = Date::from_ymd(first_year + 1, 1, 1)?.days_since(start);
    let last_days = end.days_since(Date::from_ymd(last_year, 1, 1)?);
    Ok(f64::from(first_days) / first_denominator
        + f64::from(last_year - first_year - 1)
        + f64::from(last_days) / days_in_year(end.naive().leap_year()))
}

// 各通常期間 R_k との重なりの日数を F × |R_k| で割って足す（ICMA Rule 251）。スタブでは
// 枠が格子を延ばした仮想の通常期間を持ち、重ならない通常期間は 0 日として寄与しない。
fn icma(start: Date, end: Date, frame: &PeriodFrame) -> Result<f64> {
    let frequency = f64::from(frame.periods_per_year()?.get());
    Ok(frame
        .regular
        .iter()
        .map(|period| {
            let overlap = end
                .min(period.end())
                .days_since(start.max(period.start()))
                .max(0);
            let regular_days = f64::from(period.end().days_since(period.start()));
            f64::from(overlap) / (frequency * regular_days)
        })
        .sum())
}

fn afb(start: Date, end: Date) -> Result<f64> {
    let mut cursor = end;
    let mut whole_years = 0;
    let first_year = start.naive().year();
    while cursor.naive().year() > first_year {
        let (year, month, day) = cursor.ymd();
        let previous = if month == 2 && day >= 28 {
            Date::from_ymd(year - 1, 2, 1)?.month_end()?
        } else {
            Date::from_ymd(year - 1, month, day)?
        };
        if previous < start {
            break;
        }
        cursor = previous;
        whole_years += 1;
    }
    // 余りの分母は利払いの区間 [start, cursor) に閏日があるかで決まる。
    let leap = leap_dates(start, cursor).any(|date| start.naive() <= date && date < cursor.naive());
    Ok(f64::from(whole_years) + f64::from(cursor.days_since(start)) / days_in_year(leap))
}

// 2 月末の置き換え（月末の日程だけ）、d₂ の 31（d₁ ≥ 30 のとき）、d₁ の 31 の順に適用する。
// d₁ の 2 月末を先に 30 にするので、2 月末に始まる期間でも d₂ の 31 を 30 にできる。
fn thirty_us(start: Date, end: Date, month_end: MonthEndSchedule) -> Result<f64> {
    let (_, start_month, mut first) = start.ymd();
    let (_, end_month, mut second) = end.ymd();
    if month_end == MonthEndSchedule::Aligned && start_month == 2 && start == start.month_end()? {
        first = 30;
        if end_month == 2 && end == end.month_end()? {
            second = 30;
        }
    }
    if first >= 30 {
        second = second.min(30);
    }
    Ok(thirty(start, end, first.min(30), second))
}

// 月末の日は 30 日に置き換える。ただし満期日が 2 月末なら、その終了日は元の日のまま数える。
fn thirty_e_isda(start: Date, end: Date, maturity: Date) -> Result<f64> {
    let first = if start == start.month_end()? {
        30
    } else {
        start.ymd().2
    };
    let (_, month, day) = end.ymd();
    let second = if end == end.month_end()? && !(end == maturity && month == 2) {
        30
    } else {
        day
    };
    Ok(thirty(start, end, first, second))
}

fn thirty(start: Date, end: Date, first: u8, second: u8) -> f64 {
    let (first_year, first_month, _) = start.ymd();
    let (last_year, last_month, _) = end.ymd();
    let numerator = 360 * (last_year - first_year)
        + 30 * (i32::from(last_month) - i32::from(first_month))
        + i32::from(second)
        - i32::from(first);
    f64::from(numerator) / 360.0
}
`````

`````
crates/rhodium-core/src/error.rs
//! ライブラリで共通に使うエラーと結果の型。

use std::{fmt::Display, num::NonZeroU32};

use crate::{
    date::Date,
    fixing::ObservationPoint,
    fx::CcyPair,
    id::{CalendarId, IndexId, QuoteId},
};

/// 計算や入力の読み込みに成功した値、またはその失敗理由。
pub type Result<T> = std::result::Result<T, Error>;

/// 成功した計算で使った近似や外挿の情報。
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Warning {
    /// 休日表の範囲を外れたため、週末の規則だけで判定した。
    CalendarExtrapolated {
        /// 外挿した基本カレンダーの ID。
        calendar: CalendarId,
        /// 外挿した日付。
        date: Date,
    },
    /// 確認書の支払日を休日のまま使った。
    ExplicitDateIsHoliday {
        /// 確認書の支払日。
        date: Date,
    },
    /// 公表前の値を、評価時点以降の残りの期間から予想する。
    UnpublishedFixingForecast {
        /// 予想する指標。
        index: IndexId,
        /// 予想の対象。
        point: ObservationPoint,
    },
    /// 公表前の観測に入力された推定値を使った。
    FixingEstimateUsed {
        /// 推定値を使った指標。
        index: IndexId,
        /// 推定値の対象。
        point: ObservationPoint,
    },
    /// 取引の規定に従って代替指標を使った。
    FallbackApplied {
        /// 元の指標。
        index: IndexId,
        /// 代替指標。
        replacement: IndexId,
    },
    /// 代表性を失った後も元の指標を使った。
    NonRepresentativeIndex {
        /// 代表性を失った指標。
        index: IndexId,
    },
    /// 観測期間の終了と支払日のずれに伴うタイミング調整を省いた。
    PaymentDelayIgnored,
    /// 異なる割引・予想曲線で、担保の違いによる予想の調整を省いた。
    CrossCollateralProjection,
    /// FRA の予想と割引の曲線差に伴う凸性を省いた。
    FraDiscounting,
    /// FX リセットと金利の共分散を省いた。
    MtmFxForwardProduct,
    /// 予想の日次因子を積にし、凸性の調整を省いた。
    DailyForwardProduct,
    /// 契約の丸めを予想した金額または率に適用した。
    RoundingAtForward,
    /// 金利かスプレッドとして読む quote の絶対値が 1（100%）を超える。
    ///
    /// 5% を単位なしの 5 と書いたような、単位の付け忘れの疑いがある。値はそのまま使う。
    QuoteUnitSuspected {
        /// 疑わしい quote。
        quote: QuoteId,
    },
    /// 直接の評価時点 FX と宣言した三角経路の値に差がある。
    FxCrossInconsistent {
        /// 直接のレートを使うペア。
        pair: CcyPair,
    },
    /// 最小二乗の校正で、残差が binding の許容幅を超えた（当てはめの許容幅の外）。
    ///
    /// 過剰決定の当てはめでは正しい解でも残差が残るので、解は使い、quote を報告する。
    CalibrationFitOutsideTolerance {
        /// 残差が許容幅を超えた quote。
        quote: QuoteId,
    },
    /// 最小二乗のグループの解の記録を Gauss–Newton（残差の 2 階の項を省く）で作った。
    ///
    /// このグループを通る quote の感応度は近似になる。厳密な記録は 2 階の項を含む。
    CalibrationGaussNewton {
        /// グループが読む quote。
        quote: QuoteId,
    },
}

/// DSL のソース上でエラーが起きた位置。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Span {
    /// ソースファイルの名前。
    pub file: String,
    /// 1 から数えた行番号。
    pub line: NonZeroU32,
    /// 1 から数えた列番号。
    pub column: NonZeroU32,
}

/// 入力、データ、計算の失敗を呼び出し元に伝えるエラー。
#[derive(Debug, Clone, thiserror::Error)]
pub enum Error {
    /// API に渡された入力が使い方の規則に反する。
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// 読み込んだデータが形式や整合性の規則に反する。
    #[error("invalid data: {0}")]
    InvalidData(String),
    /// 評価に必要な参照データ、fixing、quote などが揃っていない。
    #[error("missing data: {0}")]
    MissingData(String),
    /// 校正が収束しなかった、または解を求められなかった。
    #[error("calibration failed: {0}")]
    Calibration(String),
    /// 数値計算が定義域外の値、特異な行列などで失敗した。
    #[error("numerical failure: {0}")]
    Numerical(String),
    /// 指定した評価方式が商品や要求された計算を扱えない。
    #[error("unsupported: {0}")]
    Unsupported(String),
    /// DSL の構文、型、参照などに誤りがある。
    #[error("{message}")]
    Dsl {
        /// 誤りが起きたソース上の位置。
        span: Span,
        /// 誤りの説明。
        message: String,
        /// 修正方法を示せる場合の提案。
        help: Option<String>,
    },
    /// 並列の実行が止まったか、ワーカーを動かせなかった（取り消し、期限、スレッドや領域の確保の失敗）。
    #[error("execution failed: {0}")]
    Execution(String),
    /// ライブラリ内部で前提に反する状態が見つかった。
    #[error("internal error: {0}")]
    Internal(String),
}

impl Error {
    /// エラーの種類を保ったまま、対象の取引や処理などの文脈を説明の先頭に付ける。
    #[must_use]
    pub fn context(mut self, what: impl Display) -> Self {
        let message = match &mut self {
            Self::InvalidInput(message)
            | Self::InvalidData(message)
            | Self::MissingData(message)
            | Self::Calibration(message)
            | Self::Numerical(message)
            | Self::Unsupported(message)
            | Self::Dsl { message, .. }
            | Self::Execution(message)
            | Self::Internal(message) => message,
        };
        *message = format!("{what}: {message}");
        self
    }
}
`````

`````
crates/rhodium-core/src/fixing.rs
//! 指標の観測点。

use std::{collections::BTreeMap, fmt, str::FromStr, sync::Arc};

use crate::{
    Error, Result, Warning,
    date::{Date, DateInterval},
    id::{CreditKey, IndexId},
    quantity::RecoveryRate,
    refdata::RefData,
    time::Instant,
};

/// 月次指標を特定する参照月。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct YearMonth(Date);

impl YearMonth {
    /// 年と月から参照月を作る。
    ///
    /// # Errors
    ///
    /// 年月が日付の値域に収まらない場合。
    pub fn new(year: i32, month: u8) -> Result<Self> {
        Ok(Self(Date::from_ymd(year, month, 1)?))
    }

    /// 参照月の初日。
    #[must_use]
    pub const fn first_day(self) -> Date {
        self.0
    }
}

/// 指標の値を指定する観測日、瞬間、または参照月。
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "snake_case")]
pub enum ObservationPoint {
    /// 日次の fixing の日。
    Date(Date),
    /// 日中の観測の瞬間。
    Instant(Instant),
    /// 月次指標の参照月。
    Month(YearMonth),
}

impl FromStr for YearMonth {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        if text.len() != 7 {
            return Err(Error::InvalidData("reference month must be YYYY-MM".into()));
        }
        let date = format!("{text}-01").parse::<Date>()?;
        Ok(Self(date))
    }
}

impl fmt::Display for YearMonth {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (year, month, _) = self.0.ymd();
        write!(f, "{year:04}-{month:02}")
    }
}

serde_text!(YearMonth);

/// 一つの信用事由の発生日と種類。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CreditEventRecord {
    /// 信用事由が起きた現地の日付。
    pub event_date: Date,
    /// 信用契約で用いる事由の種類。
    pub kind: String,
}

/// 一つの指標と観測点について入力された値。
pub type FixingEntry = (IndexId, ObservationPoint, f64);

/// 評価時点に合わせて呼び出し側が選んだ、検証前の fixing のデータ。
#[derive(Debug, Clone, Default)]
pub struct FixingStoreParts {
    /// 公表された実績値。
    pub published: Vec<FixingEntry>,
    /// 公表前の観測の推定値。
    pub estimates: Vec<FixingEntry>,
    /// 発生した信用事由。
    pub credit_events: Vec<(CreditKey, CreditEventRecord)>,
    /// 信用事由の回収率。
    pub recoveries: Vec<(CreditKey, RecoveryRate)>,
}

type Fixings = BTreeMap<IndexId, BTreeMap<ObservationPoint, f64>>;

/// 公表値、推定値、信用事由と回収率の不変のスナップショット。
#[derive(Debug, Clone, Default)]
pub struct FixingStore {
    published: Arc<Fixings>,
    estimates: Arc<Fixings>,
    credit_events: Arc<BTreeMap<CreditKey, CreditEventRecord>>,
    recoveries: Arc<BTreeMap<CreditKey, RecoveryRate>>,
}

/// 一つの観測が対象とする期間、または日中の観測の瞬間。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ObservationPeriod {
    start: Instant,
    end: Instant,
}

impl ObservationPeriod {
    /// 観測期間の開始。
    #[must_use]
    pub const fn start(self) -> Instant {
        self.start
    }

    /// 観測期間の終了。
    #[must_use]
    pub const fn end(self) -> Instant {
        self.end
    }

    /// 開始と終了の瞬間で観測の期間を表す。
    ///
    /// # Errors
    ///
    /// 終了が開始より前なら `InvalidInput`。
    pub fn new(start: Instant, end: Instant) -> Result<Self> {
        if end < start {
            return Err(Error::InvalidInput(
                "observation period ends before it starts".into(),
            ));
        }
        Ok(Self { start, end })
    }

    /// 一つの瞬間を対象とする観測。
    #[must_use]
    pub const fn at(at: Instant) -> Self {
        Self { start: at, end: at }
    }
}

/// 評価時点、公表時刻、観測の期間で決まる値の取得方法。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FixingResolution {
    /// 公表済みの実績。
    Published(f64),
    /// 公表前の観測の推定値。
    Estimated(f64),
    /// 市場から観測期間の値を予想する。
    Forecast,
    /// 指定した瞬間から観測期間の終了までの市場予想を使う。
    ForecastFrom(Instant),
}

impl FixingStore {
    /// 参照データで指標と観測点を照合し、呼び出し側が選んだ値からスナップショットを作る。
    ///
    /// # Errors
    ///
    /// 未登録の指標、観測点の種類の不一致、非有限値、同じ表の重複は `InvalidData`。
    pub fn new(refdata: &RefData, parts: FixingStoreParts) -> Result<Self> {
        let FixingStoreParts {
            published,
            estimates,
            credit_events,
            recoveries,
        } = parts;
        Ok(Self {
            published: Arc::new(collect_fixings(refdata, published)?),
            estimates: Arc::new(collect_fixings(refdata, estimates)?),
            credit_events: Arc::new(collect_credit_events(refdata, credit_events)?),
            recoveries: Arc::new(collect_credit(refdata, recoveries)?),
        })
    }

    /// 指定した観測点の公表値を読む。
    #[must_use]
    pub fn published(&self, index: &IndexId, point: &ObservationPoint) -> Option<f64> {
        self.published.get(index)?.get(point).copied()
    }

    /// 指定した観測点の推定値を読む。公表判定を含む選択は resolve が行う。
    #[must_use]
    pub fn estimate(&self, index: &IndexId, point: &ObservationPoint) -> Option<f64> {
        self.estimates.get(index)?.get(point).copied()
    }

    /// 信用窓 (start, end] に起きた事由を読む。
    #[must_use]
    pub fn credit_event(
        &self,
        credit: &CreditKey,
        window: DateInterval,
    ) -> Option<&CreditEventRecord> {
        self.credit_events
            .get(credit)
            .filter(|event| window.start() < event.event_date && event.event_date <= window.end())
    }

    /// 信用事由の回収率を読む。
    #[must_use]
    pub fn recovery(&self, credit: &CreditKey) -> Option<RecoveryRate> {
        self.recoveries.get(credit).copied()
    }

    /// 推定値を追加または置換した新しいスナップショットを返す。
    ///
    /// # Errors
    ///
    /// 未登録の指標、観測点の不一致、非有限値、一度に入力した値の重複は `InvalidData`。
    pub fn with_estimates(
        &self,
        refdata: &RefData,
        estimates: impl IntoIterator<Item = FixingEntry>,
    ) -> Result<Self> {
        let mut updated = self.clone();
        merge_fixings(&mut updated.estimates, collect_fixings(refdata, estimates)?);
        Ok(updated)
    }

    /// シナリオの公表値と信用事由を追加または置換した、別のスナップショットを作る。
    ///
    /// # Errors
    ///
    /// 追加する値の参照先が存在しない、値や観測点が不適切、同じキーが入力に複数あれば `InvalidData`。
    pub fn with_assumed_published(
        &self,
        refdata: &RefData,
        fixings: impl IntoIterator<Item = FixingEntry>,
        events: impl IntoIterator<Item = (CreditKey, CreditEventRecord)>,
    ) -> Result<Self> {
        let fixings = collect_fixings(refdata, fixings)?;
        let events = collect_credit_events(refdata, events)?;
        let mut updated = self.clone();
        merge_fixings(&mut updated.published, fixings);
        if !events.is_empty() {
            Arc::make_mut(&mut updated.credit_events).extend(events);
        }
        Ok(updated)
    }

    /// 解決済みの公表の瞬間と観測期間から、評価で使う値の取得方法を決める。
    ///
    /// # Errors
    ///
    /// 公表後に実績が欠ける場合、または公表前で終了済みの観測の推定値が欠ける場合は `MissingData`。
    pub fn resolve(
        &self,
        index: &IndexId,
        point: ObservationPoint,
        as_of: Instant,
        publication: Instant,
        period: ObservationPeriod,
    ) -> Result<(FixingResolution, Option<Warning>)> {
        if as_of >= publication {
            return self
                .published(index, &point)
                .map(|value| (FixingResolution::Published(value), None))
                .ok_or_else(|| {
                    Error::MissingData(format!("missing published fixing {index} {point:?}"))
                });
        }
        if as_of <= period.start {
            return Ok((FixingResolution::Forecast, None));
        }
        if let Some(value) = self.estimate(index, &point) {
            return Ok((
                FixingResolution::Estimated(value),
                Some(Warning::FixingEstimateUsed {
                    index: index.clone(),
                    point,
                }),
            ));
        }
        if as_of < period.end {
            return Ok((
                FixingResolution::ForecastFrom(as_of),
                Some(Warning::UnpublishedFixingForecast {
                    index: index.clone(),
                    point,
                }),
            ));
        }
        Err(Error::MissingData(format!(
            "missing fixing estimate {index} {point:?}"
        )))
    }
}

fn collect_fixings(
    refdata: &RefData,
    rows: impl IntoIterator<Item = FixingEntry>,
) -> Result<Fixings> {
    let mut result = Fixings::new();
    for (index, point, value) in rows {
        refdata.validate_point(&index, point)?;
        if !value.is_finite() {
            return Err(Error::InvalidData(format!(
                "non-finite fixing {index} {point:?}"
            )));
        }
        if result
            .entry(index.clone())
            .or_default()
            .insert(point, value)
            .is_some()
        {
            return Err(Error::InvalidData(format!(
                "duplicate fixing {index} {point:?}"
            )));
        }
    }
    Ok(result)
}

fn merge_fixings(target: &mut Arc<Fixings>, added: Fixings) {
    if !added.is_empty() {
        let target = Arc::make_mut(target);
        for (index, points) in added {
            target.entry(index).or_default().extend(points);
        }
    }
}

fn collect_credit_events(
    refdata: &RefData,
    rows: impl IntoIterator<Item = (CreditKey, CreditEventRecord)>,
) -> Result<BTreeMap<CreditKey, CreditEventRecord>> {
    let rows = collect_credit(refdata, rows)?;
    if rows.values().any(|event| event.kind.trim().is_empty()) {
        return Err(Error::InvalidData("credit event kind is empty".into()));
    }
    Ok(rows)
}

fn collect_credit<T>(
    refdata: &RefData,
    rows: impl IntoIterator<Item = (CreditKey, T)>,
) -> Result<BTreeMap<CreditKey, T>> {
    let mut result = BTreeMap::new();
    for (key, value) in rows {
        refdata.currency(key.ccy)?;
        if result.insert(key, value).is_some() {
            return Err(Error::InvalidData("duplicate credit key".into()));
        }
    }
    Ok(result)
}
`````

`````
crates/rhodium-core/src/fx.rs
//! 通貨ペアの向きと受渡区分を持つ FX の量。

use std::{fmt, str::FromStr};

use scioxide_ad::{Mask, Real};

use crate::{
    Error, Result,
    currency::Currency,
    date::Date,
    money::{Money, checked_amount},
    time::Instant,
};

/// 1 base あたりの quote を表す、向き付きの通貨ペア。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CcyPair {
    base: Currency,
    quote: Currency,
}

impl CcyPair {
    /// 異なる二つの通貨を、base、quote の順に組み合わせる。
    ///
    /// # Errors
    ///
    /// 同じ通貨なら `InvalidInput`。
    pub fn new(base: Currency, quote: Currency) -> Result<Self> {
        if base == quote {
            return Err(Error::InvalidInput(format!("FX pair repeats {base}")));
        }
        Ok(Self { base, quote })
    }

    /// 1 単位を受け渡す側の通貨。
    #[must_use]
    pub const fn base(self) -> Currency {
        self.base
    }

    /// 価格を表す側の通貨。
    #[must_use]
    pub const fn quote(self) -> Currency {
        self.quote
    }

    /// base と quote を入れ替える。
    #[must_use]
    pub const fn inverse(self) -> Self {
        Self {
            base: self.quote,
            quote: self.base,
        }
    }
}

impl FromStr for CcyPair {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        let (base, quote) = text
            .split_at_checked(3)
            .ok_or_else(|| Error::InvalidData("invalid FX pair".into()))?;
        Self::new(base.parse()?, quote.parse()?)
            .map_err(|error| Error::InvalidData(error.to_string()))
    }
}

impl fmt::Display for CcyPair {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.base, self.quote)
    }
}

serde_text!(CcyPair);

/// 評価時点での受け渡し。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Today;

/// ペアのスポット日での受け渡し。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Spot(Date);

/// フォワード契約が受け渡される瞬間。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Forward(Instant);

/// 受渡区分と向きを持つ正のレート。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FxRate<D, R = f64> {
    pair: CcyPair,
    rate: R,
    delivery: D,
}

impl<D: Copy + Eq, R: Real> FxRate<D, R> {
    /// 通貨ペア。
    #[must_use]
    pub const fn pair(&self) -> CcyPair {
        self.pair
    }

    /// base 1 単位あたりの quote の量。
    #[must_use]
    pub const fn rate(&self) -> R {
        self.rate
    }

    /// 同じ受渡条件でレートの向きを反転する。
    ///
    /// # Errors
    ///
    /// 逆数が有限な正の数でなければ `Numerical`。
    pub fn invert(self) -> Result<Self> {
        let pair = self.pair.inverse();
        Self::checked(pair, R::one() / self.rate, self.delivery, numerical_rate)
    }

    /// 共通の中間通貨を持つ同じ受渡条件のレートをつなぐ。
    ///
    /// # Errors
    ///
    /// 通貨の向きか受渡が合わなければ `InvalidInput`、積が正の有限値でなければ `Numerical`。
    pub fn cross(self, next: &Self) -> Result<Self> {
        if self.pair.quote != next.pair.base || self.delivery != next.delivery {
            return Err(Error::InvalidInput(format!(
                "cannot cross FX rates {} and {} with these currencies and deliveries",
                self.pair, next.pair
            )));
        }
        let pair = CcyPair::new(self.pair.base, next.pair.quote)?;
        Self::checked(pair, self.rate * next.rate, self.delivery, numerical_rate)
    }

    fn checked(pair: CcyPair, rate: R, delivery: D, error: fn() -> Error) -> Result<Self> {
        if !rate.is_finite().all() || !rate.gt(R::zero()).all() {
            return Err(error());
        }
        Ok(Self {
            pair,
            rate,
            delivery,
        })
    }
}

impl<R: Real> FxRate<Today, R> {
    /// 評価時点での受け渡しに使うレートを作る。
    ///
    /// # Errors
    ///
    /// レートが正の有限値でなければ `InvalidInput`。
    pub fn today(pair: CcyPair, rate: R) -> Result<Self> {
        Self::checked(pair, rate, Today, invalid_rate)
    }

    /// base の金額を quote の金額に換算する。
    ///
    /// # Errors
    ///
    /// base の通貨でなければ `InvalidInput`、結果が非有限なら `Numerical`。
    pub fn convert(&self, money: Money<R>) -> Result<Money<R>> {
        if money.ccy() != self.pair.base {
            return Err(Error::InvalidInput(format!(
                "FX rate {} cannot convert {} from its base",
                self.pair,
                money.ccy()
            )));
        }
        checked_amount(self.pair.quote, money.amount() * self.rate)
    }
}

impl<R: Real> FxRate<Spot, R> {
    /// 指定したスポット日での受け渡しに使うレートを作る。
    ///
    /// # Errors
    ///
    /// レートが正の有限値でなければ `InvalidInput`。
    pub fn spot(pair: CcyPair, rate: R, spot_date: Date) -> Result<Self> {
        Self::checked(pair, rate, Spot(spot_date), invalid_rate)
    }

    /// スポット日。
    #[must_use]
    pub const fn spot_date(&self) -> Date {
        self.delivery.0
    }
}

impl<R: Real> FxRate<Forward, R> {
    /// 指定した受渡のフォワードを作る。
    ///
    /// # Errors
    /// レートが正の有限値でなければ `InvalidInput`。
    pub fn forward(pair: CcyPair, rate: R, at: Instant) -> Result<Self> {
        Self::checked(pair, rate, Forward(at), invalid_rate)
    }

    /// フォワードの受渡時刻。
    #[must_use]
    pub const fn settlement(&self) -> Instant {
        self.delivery.0
    }
}

/// スポット日からフォワード受渡日までのレート差。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FxForwardPoints {
    pair: CcyPair,
    points: f64,
    pip_exponent: u8,
    spot_date: Date,
    value_date: Date,
}

impl FxForwardPoints {
    /// ポイントの値と尺度、スポット日と受渡日を結び付ける。
    ///
    /// # Errors
    ///
    /// ポイントが非有限か、受渡日がスポット日より前なら `InvalidInput`。
    pub fn new(
        pair: CcyPair,
        points: f64,
        pip_exponent: u8,
        spot_date: Date,
        value_date: Date,
    ) -> Result<Self> {
        if !points.is_finite() || value_date < spot_date {
            return Err(Error::InvalidInput(
                "invalid forward points or value date".into(),
            ));
        }
        Ok(Self {
            pair,
            points,
            pip_exponent,
            spot_date,
            value_date,
        })
    }

    /// フォワード受渡日。
    #[must_use]
    pub const fn value_date(&self) -> Date {
        self.value_date
    }

    /// スポットに尺度を適用したポイントを加え、アウトライトを求める。
    ///
    /// # Errors
    ///
    /// ペアかスポット日が合わなければ `InvalidInput`、結果が正の有限値でなければ `Numerical`。
    pub fn outright(&self, spot: &FxRate<Spot>) -> Result<f64> {
        if self.pair != spot.pair || self.spot_date != spot.spot_date() {
            return Err(Error::InvalidInput(
                "forward points and spot must have the same pair and spot date".into(),
            ));
        }
        let rate = spot.rate + self.points * 10.0_f64.powi(-i32::from(self.pip_exponent));
        if !rate.is_finite() || rate <= 0.0 {
            return Err(numerical_rate());
        }
        Ok(rate)
    }
}

fn invalid_rate() -> Error {
    Error::InvalidInput("FX rate must be finite and positive".into())
}

fn numerical_rate() -> Error {
    Error::Numerical("calculated FX rate is not finite and positive".into())
}
`````

`````
crates/rhodium-core/src/id.rs
//! 参照データを引くための識別子。

use std::{fmt, str::FromStr, sync::Arc};

use crate::{Error, Result, currency::Currency};

fn valid_id(text: &str) -> bool {
    (1..=128).contains(&text.len())
        && text
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || b"-_.:/+".contains(&byte))
}

fn valid_entity(text: &str) -> bool {
    let Some(lei) = text.strip_prefix("LEI:") else {
        return true;
    };
    let digits = lei.get(18..).and_then(|suffix| suffix.parse::<u8>().ok());
    lei.len() == 20
        && lei
            .get(18..)
            .is_some_and(|suffix| suffix.bytes().all(|byte| byte.is_ascii_digit()))
        && digits.is_some_and(|check| (2..=98).contains(&check))
        && lei.bytes().try_fold(0_u16, |remainder, byte| {
            let (shift, digit) = match byte {
                b'0'..=b'9' => (10, byte - b'0'),
                b'A'..=b'Z' => (100, byte - b'A' + 10),
                _ => return None,
            };
            Some((remainder * shift + u16::from(digit)) % 97)
        }) == Some(1)
}

macro_rules! reference_ids {
    ($(#[$doc:meta] $name:ident, $label:literal $(, $validator:ident)?;)+) => {
        $(
            #[$doc]
            #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
            pub struct $name(Arc<str>);

            impl $name {
                /// ID の文字列。
                #[must_use]
                pub fn as_str(&self) -> &str {
                    &self.0
                }
            }

            impl FromStr for $name {
                type Err = Error;

                fn from_str(text: &str) -> Result<Self> {
                    if !(valid_id(text) $(&& $validator(text))?) {
                        return Err(Error::InvalidData(concat!("invalid ", $label, " ID").into()));
                    }
                    Ok(Self(text.into()))
                }
            }

            impl fmt::Display for $name {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    f.write_str(self.as_str())
                }
            }

            serde_text!($name);
        )+
    };
}

reference_ids! {
    /// 基本カレンダーの識別子。
    CalendarId, "calendar";
    /// 観測指標の識別子。
    IndexId, "index";
    /// 商品の慣行の識別子。
    ConventionId, "convention";
    /// 取引の識別子。
    TradeId, "trade";
    /// 取引が属するブックの識別子。
    BookId, "book";
    /// 取引が属するマージンセットの識別子。
    MarginSetId, "margin set";
    /// 信用主体の識別子。LEI: で始まる場合は検査桁を照合する。
    EntityId, "entity", valid_entity;
    /// 証券の識別子。
    SecurityId, "security";
    /// 受渡場所を含むコモディティの識別子。
    CommodityId, "commodity";
    /// 債務の優先順位の識別子。
    Seniority, "seniority";
    /// 信用契約の条項の識別子。
    DocClause, "documentation clause";
}

/// 2 個以上の区間を `.` で区切った quote の識別子。
///
/// 各区間には英大文字、数字、`_`、`+`、`-` を使える。値の意味は校正の binding が定める。
/// 市場の quote は `rhodium_market::quote` が扱うが、quote を指す警告（[`crate::Warning`]）が
/// 型付きの ID を持てるよう、型はこの crate で定義する。
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QuoteId(Arc<str>);

impl QuoteId {
    /// ID の文字列。
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }

    /// 文字列が quote ID の区間一つとして有効かを返す。
    #[must_use]
    pub fn is_valid_segment(segment: &str) -> bool {
        !segment.is_empty()
            && segment
                .bytes()
                .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b"_+-".contains(&b))
    }
}

impl FromStr for QuoteId {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        if !text.contains('.') || !text.split('.').all(Self::is_valid_segment) {
            return Err(Error::InvalidInput(format!("invalid quote ID: {text:?}")));
        }
        Ok(Self(text.into()))
    }
}

impl TryFrom<String> for QuoteId {
    type Error = Error;

    fn try_from(value: String) -> Result<Self> {
        value.parse()
    }
}

impl From<QuoteId> for String {
    fn from(id: QuoteId) -> Self {
        id.as_str().to_owned()
    }
}

impl fmt::Display for QuoteId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

serde_text!(QuoteId);

/// 信用主体と契約条件の組。信用事由と回収率の参照に使う。
#[derive(
    Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize,
)]
#[serde(deny_unknown_fields)]
pub struct CreditKey {
    /// 信用主体。
    pub entity: EntityId,
    /// 債務の優先順位。
    pub seniority: Seniority,
    /// 契約の条項。
    pub doc_clause: DocClause,
    /// 契約の通貨。
    pub ccy: Currency,
}
`````

`````
crates/rhodium-core/src/lib.rs
//! 日付、金額、参照データ、エラーなどの共通基盤。

#![forbid(unsafe_code)]

// serde の読み込みも文字列の検証済みコンストラクタを通す。
macro_rules! serde_text {
    ($type:ty) => {
        impl serde::Serialize for $type {
            fn serialize<S: serde::Serializer>(
                &self,
                serializer: S,
            ) -> std::result::Result<S::Ok, S::Error> {
                serializer.collect_str(self)
            }
        }

        impl<'de> serde::Deserialize<'de> for $type {
            fn deserialize<D: serde::Deserializer<'de>>(
                deserializer: D,
            ) -> std::result::Result<Self, D::Error> {
                <String as serde::Deserialize>::deserialize(deserializer)?
                    .parse()
                    .map_err(serde::de::Error::custom)
            }
        }
    };
}

pub mod brand;
pub mod calendar;
pub mod currency;
pub mod date;
pub mod day_count;
pub mod error;
pub mod fixing;
pub mod fx;
pub mod id;
pub mod money;
pub mod observation;
pub mod parallel;
pub mod quantity;
pub mod refdata;
pub mod rounding;
pub mod schedule;
pub mod tenor;
pub mod time;

pub use error::{Error, Result, Warning};
`````

`````
crates/rhodium-core/src/money.rs
//! 通貨タグ付きの金額と、通貨ごとの集計。

use std::collections::BTreeMap;

use scioxide_ad::{Mask, Real};

use crate::{
    Error, Result,
    currency::Currency,
    fx::{FxRate, Today},
};

/// 単一の通貨の金額。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Money<R = f64> {
    ccy: Currency,
    amount: R,
}

impl std::str::FromStr for Money {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        let mut parts = text.split_ascii_whitespace();
        let (Some(number), Some(currency), None) = (parts.next(), parts.next(), parts.next())
        else {
            return Err(Error::InvalidData(
                "money requires an amount and currency".into(),
            ));
        };
        let bytes = number.as_bytes();
        if bytes.iter().enumerate().any(|(i, byte)| {
            *byte == b'_'
                && !(i > 0
                    && bytes[i - 1].is_ascii_digit()
                    && bytes.get(i + 1).is_some_and(u8::is_ascii_digit))
        }) {
            return Err(Error::InvalidData(
                "money digit separators must be between digits".into(),
            ));
        }
        let number = number.replace('_', "");
        let unsigned = number.strip_prefix(['+', '-']).unwrap_or(&number);
        let (whole, fraction) = unsigned.split_once('.').unwrap_or((unsigned, ""));
        if whole.is_empty() && fraction.is_empty()
            || !whole
                .bytes()
                .chain(fraction.bytes())
                .all(|b| b.is_ascii_digit())
        {
            return Err(Error::InvalidData("money amount must be a decimal".into()));
        }
        let amount = number
            .parse::<f64>()
            .ok()
            .filter(|value| value.is_finite())
            .ok_or_else(|| Error::InvalidData("money amount must be finite".into()))?;
        Ok(Self::new(currency.parse()?, amount))
    }
}

impl<R: Real> Money<R> {
    /// 通貨を金額に結び付ける。
    #[must_use]
    pub const fn new(ccy: Currency, amount: R) -> Self {
        Self { ccy, amount }
    }

    /// 金額の通貨。
    #[must_use]
    pub const fn ccy(self) -> Currency {
        self.ccy
    }

    /// 符号付きの金額。
    #[must_use]
    pub const fn amount(self) -> R {
        self.amount
    }

    /// 通貨が同じ金額を足す。
    ///
    /// # Errors
    ///
    /// 通貨が違えば `InvalidInput`、結果が非有限なら `Numerical`。
    pub fn checked_add(self, other: Self) -> Result<Self> {
        if self.ccy != other.ccy {
            return Err(Error::InvalidInput(format!(
                "cannot add {} and {} amounts",
                self.ccy, other.ccy
            )));
        }
        checked_amount(self.ccy, self.amount + other.amount)
    }

    /// 同じ通貨の金額を倍率で拡大する。
    #[must_use]
    pub fn scale(self, factor: R) -> Self {
        Self::new(self.ccy, self.amount * factor)
    }
}

/// 通貨ごとの金額の和。値がゼロになっても感応度を保持するため項を残す。
#[derive(Debug, Clone)]
pub struct MoneyBag<R = f64>(BTreeMap<Currency, R>);

impl<R> Default for MoneyBag<R> {
    fn default() -> Self {
        Self(BTreeMap::new())
    }
}

impl<R: Real> MoneyBag<R> {
    /// 項を持たない集計を作る。
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// 金額を通貨ごとに足す。
    pub fn add(&mut self, money: Money<R>) {
        let entry = self.0.entry(money.ccy()).or_insert_with(R::zero);
        *entry += money.amount();
    }

    /// 通貨の項を読む。項が無ければゼロ。
    #[must_use]
    pub fn get(&self, currency: Currency) -> R {
        self.0.get(&currency).copied().unwrap_or_else(R::zero)
    }

    /// ゼロの項も含め、通貨順に集計済みの金額を返す。
    pub fn terms(&self) -> impl Iterator<Item = Money<R>> + '_ {
        self.0.iter().map(|(&ccy, &amount)| Money::new(ccy, amount))
    }

    /// 各通貨の金額を、指定した向きの直接レートで報告通貨に換算する。
    ///
    /// # Errors
    ///
    /// 直接レートの欠損は `MissingData`、重複は `InvalidInput`、非有限の結果は `Numerical`。
    pub fn convert(&self, to: Currency, rates: &[FxRate<Today, R>]) -> Result<Money<R>> {
        let total = self
            .terms()
            .try_fold(Money::new(to, R::zero()), |sum, money| {
                if money.ccy() == to {
                    return sum.checked_add(money);
                }
                let mut direct = rates
                    .iter()
                    .filter(|rate| rate.pair().base() == money.ccy() && rate.pair().quote() == to);
                let rate = direct.next().ok_or_else(|| {
                    Error::MissingData(format!("missing direct FX rate {}{to}", money.ccy()))
                })?;
                if direct.next().is_some() {
                    return Err(Error::InvalidInput(format!(
                        "duplicate direct FX rate {}{to}",
                        money.ccy()
                    )));
                }
                sum.checked_add(rate.convert(money)?)
            })?;
        Ok(total)
    }
}

impl PartialEq for MoneyBag<f64> {
    fn eq(&self, other: &Self) -> bool {
        // 値としての等値では、存在しない通貨をゼロとみなす。
        self.0.keys().chain(other.0.keys()).all(|&ccy| {
            self.get(ccy).partial_cmp(&other.get(ccy)) == Some(std::cmp::Ordering::Equal)
        })
    }
}

pub(crate) fn checked_amount<R: Real>(ccy: Currency, value: R) -> Result<Money<R>> {
    if !value.is_finite().all() {
        return Err(Error::Numerical(format!("non-finite {ccy} amount")));
    }
    Ok(Money::new(ccy, value))
}
`````

`````
crates/rhodium-core/src/observation.rs
//! 観測点と重みを日数基準とカレンダーから導く。

use std::{collections::BTreeMap, num::NonZeroU8};

use scioxide_ad::{Mask, Real};

use crate::{
    Error, Result, Warning,
    calendar::{BaseCalendar, Calendar, CalendarSpec, shift_business, unique_warnings},
    date::{Date, DateInterval},
    day_count::{AccrualContext, DayCount, year_fraction},
    fixing::{ObservationPoint, YearMonth},
    id::CalendarId,
    quantity::YearFraction,
    refdata::CalendarCache,
    rounding::Rounding,
    schedule::Period,
};

/// 一つの計算期間で観測点を選ぶ規則。
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum ObservationDays {
    /// 指標のカレンダーで営業日を選ぶ。
    BusinessDays(CalendarSpec),
    /// 各暦日を観測し、非営業日には直前の営業日の値を使う。
    AllDays {
        /// 指標の fixing カレンダー。
        fixing_calendar: CalendarSpec,
    },
    /// 計算期間と重なる参照月。
    Months,
    /// 複利指数の水準を観測する両端。
    Endpoints(EndpointDates),
}

/// 複利指数の観測期間。
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum EndpointDates {
    /// 計算期間の両端。
    Accrual,
    /// 両端を営業日数で前へずらす。
    Shifted {
        /// 複利指数のカレンダー。
        calendar: CalendarSpec,
        /// 前に戻る営業日数。
        lag: NonZeroU8,
    },
}

/// 計算期間から fixing 日と重みの区間を決める方法。
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, Default, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum ObservationMethod {
    /// 計算期間の観測点をそのまま使う。
    #[default]
    Plain,
    /// fixing 日だけを指定の観測日数だけ戻す。
    Lookback(NonZeroU8),
    /// 観測期間全体を指定の営業日数だけ戻す。
    ObservationShift(NonZeroU8),
    /// 最後の観測に締めの日の fixing を使う。値は締めの日の後の観測数。
    Lockout(NonZeroU8),
}

/// 観測値の集約に使う複利因子。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum CompoundingForm {
    /// 1 + 観測金利 × 重み。
    SimpleFactor,
    /// (1 + 観測金利) の重み乗。
    PowerFactor,
}

impl CompoundingForm {
    /// 丸めを適用した観測値と契約の重みから、一日の複利因子を求める。
    ///
    /// `spread` は日次因子に入れるスプレッドで、`SimpleFactor` の (r + spread)·w + 1 だけが使う。
    /// `PowerFactor` の (1 + r)^w はスプレッドを含まない。`PowerFactor` と `Inclusive` の組は
    /// [`ObservationRule::new`] が拒否するので、その規則から作った値は常に 0 になる。
    ///
    /// # Errors
    /// 因子が正の有限値にならなければ `Numerical`。
    pub fn factor<R: Real>(self, rate: R, weight: f64, spread: f64) -> Result<R> {
        let factor = match self {
            Self::SimpleFactor => (rate + spread) * weight + 1.0,
            Self::PowerFactor => {
                let base = rate + 1.0;
                if !base.gt(R::zero()).all() {
                    return Err(Error::Numerical("nonpositive compounding base".into()));
                }
                base.powf(R::constant(weight))
            }
        };
        if factor.is_finite().all() && factor.gt(R::zero()).all() {
            Ok(factor)
        } else {
            Err(Error::Numerical(
                "compounding factor must be finite and positive".into(),
            ))
        }
    }

    /// [`Self::factor`] − 1（因子の超過）を、1 に近い因子でも桁を落とさずに求める。
    ///
    /// 因子 1 + x を丸めてから 1 を引くと x の下位の桁が ulp(1) で失われるので、
    /// `SimpleFactor` は x = (r + spread)·w をそのまま返し（四則演算だけなので、どの環境でも同じ値）、
    /// `PowerFactor` は `exp_m1(w·ln_1p(r))` で求める。確定した因子を超過成長 Π − 1 に畳むときに使う。
    ///
    /// # Errors
    /// 因子が正の有限値にならなければ `Numerical`。
    pub fn factor_excess<R: Real>(self, rate: R, weight: f64, spread: f64) -> Result<R> {
        let excess = match self {
            Self::SimpleFactor => {
                let increment = (rate + spread) * weight;
                // 1 + x > 0 の条件は x > −1 と同じなので、1 を足さずに判定する。
                if !(increment.is_finite().all() && increment.gt(R::constant(-1.0)).all()) {
                    return Err(Error::Numerical(
                        "compounding factor must be finite and positive".into(),
                    ));
                }
                increment
            }
            Self::PowerFactor => (self.log_factor(rate, weight, spread)?).exp_m1(),
        };
        if excess.is_finite().all() {
            Ok(excess)
        } else {
            Err(Error::Numerical(
                "compounding factor must be finite and positive".into(),
            ))
        }
    }

    /// [`Self::factor`] の自然対数を、1 に近い因子でも桁を落とさずに求める。
    ///
    /// 因子 1 + x を丸めてから対数を取ると x の下位の桁が ulp(1) で失われるので、
    /// `SimpleFactor` は ln(1 + (r + spread)·w) を、`PowerFactor` は w·ln(1 + r) を
    /// `ln_1p` で直接求める。超過成長 Π 因子 − 1 を `exp_m1` で求めるときに使う。
    ///
    /// # Errors
    /// 因子が正の有限値にならなければ `Numerical`。
    pub fn log_factor<R: Real>(self, rate: R, weight: f64, spread: f64) -> Result<R> {
        // 1 + x > 0 の条件は x > −1 と同じなので、1 を足さずに判定する。
        let increment = match self {
            Self::SimpleFactor => (rate + spread) * weight,
            Self::PowerFactor => rate,
        };
        if !(increment.is_finite().all() && increment.gt(R::constant(-1.0)).all()) {
            return Err(Error::Numerical(match self {
                Self::SimpleFactor => "compounding factor must be finite and positive".into(),
                Self::PowerFactor => "nonpositive compounding base".into(),
            }));
        }
        let log = match self {
            Self::SimpleFactor => increment.ln_1p(),
            Self::PowerFactor => increment.ln_1p() * weight,
        };
        if log.is_finite().all() {
            Ok(log)
        } else {
            Err(Error::Numerical(
                "compounding factor must be finite and positive".into(),
            ))
        }
    }
}

/// 複利に対するスプレッドの扱い。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum SpreadTreatment {
    /// 年率換算後に足す。
    Exclusive,
    /// 各日の率に足す。
    Inclusive,
}

/// 平均を取るときの重み。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum AveragingWeights {
    /// 各観測が覆う区間の日数係数。
    CoverWeighted,
    /// 各観測に同じ重みを割り当てる。
    Equal,
}

/// 観測値の集約方法。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum Aggregation {
    /// 複利因子の積。
    Compound {
        /// 複利因子の式。
        form: CompoundingForm,
        /// スプレッドを掛ける段階。
        spread: SpreadTreatment,
    },
    /// 観測値の平均。
    Average(AveragingWeights),
    /// 複利指数の両端の比。
    IndexRatio,
}

/// 観測値ごとと集約後の年率に適用する丸め。
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, Default, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub struct ObservationRounding {
    /// 個別の観測値に適用する丸め。
    pub daily: Option<Rounding>,
    /// 集約後の年率に適用する丸め。
    pub aggregated: Option<Rounding>,
}

/// 観測日、方法、集約、丸めの検証済みの組。
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(try_from = "ObservationRuleInput")]
pub struct ObservationRule {
    days: ObservationDays,
    method: ObservationMethod,
    aggregation: Aggregation,
    rounding: ObservationRounding,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct ObservationRuleInput {
    days: ObservationDays,
    #[serde(default)]
