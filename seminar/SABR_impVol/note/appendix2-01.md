diff --git a/.githooks/checks.py b/.githooks/checks.py
index 70ac8ce7421e505f37354b8cfa428c80f31eca11..f364b0846c043bd13d4921cec45984553bddf26a 100644
--- a/.githooks/checks.py
+++ b/.githooks/checks.py
@@ -2,6 +2,7 @@
 
 import json
 import os
+import re
 
 # Git が渡した引数と、リポジトリ内の検証コマンドをシェルなしで実行する。
 import subprocess  # ruff: ignore[suspicious-subprocess-import]
@@ -18,6 +19,8 @@ if TYPE_CHECKING:
 
 PYTHON_PROJECTS = (".githooks",)
 CARGO_INPUTS = frozenset({"Cargo.toml", "Cargo.lock", "rust-toolchain.toml", "clippy.toml"})
+# examples/ のサンプルをテストとベンチマークで読む crate。CLI は逆依存の閉包で入る。
+EXAMPLE_READERS = ("rhodium",)
 
 type _Command = tuple[str | Path, ...]
 type _FileKind = Literal["rust", "python", "toml", "json", "nix", "shell"]
@@ -40,16 +43,249 @@ class _CheckPlan:
     files: tuple[_FileCheck, ...]
     python: tuple[_PythonCheck, ...]
     crates: frozenset[str]
+    examples: bool
     workspace: bool
 
-    def needs_push_checks(self) -> bool:
-        """送信時にツリー全体が必要になるかを判定する。
+    ui_pairs: bool = False
 
-        Returns:
-            ツリー全体を使う検査が一つでもあれば真。
 
-        """
-        return bool(self.python or self.crates or self.workspace)
+type _Mode = Literal["diff", "tree"]
+type _Scope = Literal["library", "all"]
+
+
+@dataclass(frozen=True, slots=True)
+class _Rule:
+    """型や lint で止められない実装の誤りを、Rust のソースの行の形で止める規則。
+
+    `diff` の規則は追加した行だけを照合するので、既存の違反を残したまま規則を入れられる。
+    担当の issue が既存の違反をなくしたら `tree` にし、ツリー全体を照合して後退を止める。
+    """
+
+    name: str
+    pattern: re.Pattern[str]
+    reason: str
+    # 規則の定義元の設計書の節。
+    doc: str
+    mode: _Mode = "diff"
+    # library はライブラリの src (テストのファイルと #[cfg(test)] の範囲を除く)、
+    # all は crates の下のすべての Rust のファイル。
+    scope: _Scope = "library"
+    # 空でなければ、この glob に当たるファイルだけを照合する。
+    paths: tuple[str, ...] = ()
+    # この glob に当たるファイルは照合しない。
+    allow: tuple[str, ...] = ()
+    # 当たっても違反にしない行の形。定義の行など。
+    exempt: re.Pattern[str] | None = None
+    # 違反にしない代わりに、当たった行とその前の window 行の中に要る文字列。
+    requires_nearby: str | None = None
+    window: int = 0
+    # コメントの中も照合するか。偽ならコードの部分だけを照合する。
+    comments: bool = False
+    # 関数の引数リスト内で照合し、該当する型の行を返す。
+    parameters: bool = False
+
+
+_ERRORS = "docs/design/04-foundation/04-reference-data.md#エラー"
+_TYPE_RULES = "docs/design/03-architecture.md#型の規約"
+_PERFORMANCE = "docs/design/03-architecture.md#性能の規則"
+_ONE_PLACE = "docs/design/02-principles.md#変更を一か所に収める規則"
+_CONST_LINE = re.compile(r"^\s*(?:pub(?:\([^)]*\))?\s+)?const\s")
+
+# 規則の追加と、diff から tree への切り替えの手順は docs/design/11-quality.md#ci にある。
+RULES = (
+    _Rule(
+        name="P-trybuild",
+        pattern=re.compile(r"trybuild::|\bTestCases\b"),
+        reason="compile-fail の例は compile_fail_pair! で _ok の例と対にして登録する",
+        doc="docs/design/11-quality.md#compile-fail-テスト",
+        mode="tree",
+        scope="all",
+        comments=True,
+    ),
+    _Rule(
+        name="P-ui-runner",
+        pattern=re.compile(r"\brun_compile_fail_pair\b"),
+        reason="型テストは compile_fail_pair! で対にして登録する",
+        doc="docs/design/11-quality.md#compile-fail-テスト",
+        mode="tree",
+        scope="all",
+        allow=("crates/rhodium-testkit/src/lib.rs", "crates/rhodium-testkit/src/ui.rs"),
+        comments=True,
+    ),
+    _Rule(
+        name="P-doctest-compile-fail",
+        pattern=re.compile(r"^\s*(?://[/!]|\*)\s*(?:```|~~~).*\bcompile_fail\b"),
+        reason="doctest の compile_fail は成功例と診断を対にできないので、testkit の対に移す",
+        doc="docs/design/11-quality.md#compile-fail-テスト",
+        mode="tree",
+        scope="all",
+        comments=True,
+    ),
+    _Rule(
+        name="P-internal",
+        pattern=re.compile(r"\bError::Internal\((?!_\))"),
+        reason="状態を型で表せないか確かめ、残すなら直前に「// 不変条件:」で保証の出どころを書く",
+        doc=_ERRORS,
+        requires_nearby="不変条件:",
+        window=3,
+    ),
+    _Rule(
+        name="P-path",
+        pattern=re.compile(r'#\[path\s*=\s*"\.\./'),
+        reason="他の crate やディレクトリのファイルを #[path] で取り込まず、共有の補助を使う",
+        doc="docs/design/11-quality.md#テストの補助",
+        scope="all",
+    ),
+    _Rule(
+        name="P-calendar-map",
+        pattern=re.compile(
+            r"BTreeMap<\s*(?:CalendarId\s*,\s*BaseCalendar|TzId\s*,\s*ZoneRules)\s*>"
+        ),
+        reason="カレンダーと地域は RefData と CalendarCache から解決し、生の表を受け渡さない",
+        doc="docs/design/04-foundation/02-schedules.md",
+        allow=(
+            "crates/rhodium-core/src/refdata.rs",
+            "crates/rhodium/src/workspace/refdata_files.rs",
+        ),
+    ),
+    _Rule(
+        name="P-daycount",
+        pattern=re.compile(r"/\s*(?:365|360|252)\.0\b|\b365\.25\b"),
+        reason="日数係数は DayCount と year_fraction で求め、日数の割り算を書かない",
+        doc="docs/design/04-foundation/01-dates-and-calendars.md#日数計算",
+        allow=("crates/rhodium-core/src/day_count.rs", "crates/rhodium-core/src/time.rs"),
+    ),
+    _Rule(
+        name="P-bus252",
+        pattern=re.compile(r"\bDayCount::Bus252\s*\{"),
+        reason="BUS/252 の暦は core が一度だけ解決し、下流で取り出さない",
+        doc="docs/design/04-foundation/01-dates-and-calendars.md#日数計算",
+        allow=("crates/rhodium-core/**",),
+    ),
+    _Rule(
+        name="P-let-else-kind",
+        pattern=re.compile(r"\blet\s+(?:IndexKind|Convention)::"),
+        reason="参照データの種類は種類つきの照会で一度だけ確かめ、下流で照合し直さない",
+        doc="docs/design/04-foundation/04-reference-data.md",
+        allow=("crates/rhodium-core/src/refdata.rs", "crates/rhodium-core/src/refdata/**"),
+    ),
+    _Rule(
+        name="P-zone",
+        pattern=re.compile(r"\.zone\(&"),
+        reason="現地時刻から瞬間への解決は core の一つの関数を通す",
+        doc="docs/design/04-foundation/01-dates-and-calendars.md#瞬間と現地時刻",
+        allow=("crates/rhodium-core/**",),
+    ),
+    _Rule(
+        name="P-tolerance",
+        pattern=re.compile(r"(?<![\w.])\d+(?:\.\d+)?e-\d+"),
+        reason="許容誤差は core::tolerance の定数を、アルゴリズムの定数は名前付きの const を使う",
+        doc="docs/design/11-quality.md#許容誤差",
+        allow=("crates/rhodium-core/src/tolerance.rs",),
+        exempt=_CONST_LINE,
+    ),
+    _Rule(
+        name="P-invalid",
+        pattern=re.compile(r"\bfn\s+invalid\s*\("),
+        reason="エラーの種類を名前に表さない補助を作らず、Error の関連関数を使う",
+        doc=_ERRORS,
+    ),
+    _Rule(
+        name="P-positive",
+        pattern=re.compile(r"<=\s*0\.0\b"),
+        reason="正の金額と倍率は newtype で受け、構築の後に検査し直さない",
+        doc="docs/design/04-foundation/03-money-and-quantities.md",
+        paths=("crates/rhodium-products/src/standard/**", "crates/rhodium-products/src/trade.rs"),
+    ),
+    _Rule(
+        name="P-to-product",
+        pattern=re.compile(r"\bpub\s+fn\s+to_product\b"),
+        reason="標準商品の入口は StandardProduct::to_product だけにする",
+        doc="docs/design/05-products/03-standard-products.md",
+        paths=("crates/rhodium-products/src/standard/**",),
+    ),
+    _Rule(
+        name="P-year-fraction-ir",
+        pattern=re.compile(r"\byear_fraction\s*\("),
+        reason="IR の検証は因果性と参照だけを確かめ、評価のための解析を持たない",
+        doc="docs/design/05-products/01-ir.md",
+        paths=("crates/rhodium-products/src/ir/**",),
+    ),
+    _Rule(
+        name="P-to-bits",
+        pattern=re.compile(r"\.to_bits\(\)\s*=="),
+        reason="IR の形をビットの比較で認識せず、作る側と同じ crate の分析を使う",
+        doc="docs/design/05-products/01-ir.md",
+    ),
+    _Rule(
+        name="P-caplet",
+        pattern=re.compile(r"(?i)caplet"),
+        reason="校正だけのための入口を取引の商品に足さず、標準商品とテンプレートで作る",
+        doc="docs/design/06-market/03-calibration.md",
+        paths=(
+            "crates/rhodium-products/src/standard/leg.rs",
+            "crates/rhodium-products/src/standard/rates.rs",
+        ),
+    ),
+    _Rule(
+        name="P-synthetic-id",
+        pattern=re.compile(r'format!\("(?:FX|LGM)/|"swaption\.|\.as_str\(\)\s*=='),
+        reason="ID の文字列を合成・比較して意味を持たせず、持ち主を enum で持つ",
+        doc=_TYPE_RULES,
+    ),
+    _Rule(
+        name="P-ptr-eq",
+        pattern=re.compile(r"\bArc::ptr_eq\b"),
+        reason="構造の同一性は一つの関数で確かめる",
+        doc="docs/design/06-market/02-curves.md",
+    ),
+    _Rule(
+        name="P-theta",
+        pattern=re.compile(r"\bMarket::from_(?:partial_)?theta\(|\.theta\(\)\["),
+        reason="校正のコールバックで θ 全体や市場を作り直さず、active set の市場を使う",
+        doc=_PERFORMANCE,
+        paths=("crates/rhodium-pricing/src/**",),
+    ),
+    _Rule(
+        name="P-cell-error",
+        pattern=re.compile(r"\bCell::new\(None\)"),
+        reason="コールバックのエラーの受け渡しは共通の部品を使う",
+        doc="docs/design/06-market/03-calibration.md",
+    ),
+    _Rule(
+        name="P-max-iter",
+        pattern=re.compile(r"with_max_iterations\(\s*(?:\d|Budget::\w+\.saturating_add\(\s*\d)"),
+        reason="反復の上限は名前付きの定数にする",
+        doc=_ONE_PLACE,
+    ),
+    _Rule(
+        name="P-validate-refs",
+        pattern=re.compile(r"\bvalidate_references\("),
+        reason="参照データとの照合は境界の一か所で行い、後段で照合し直さない",
+        doc="docs/design/06-market/03-calibration.md",
+        exempt=re.compile(r"\bfn\s+validate_references\b"),
+    ),
+    _Rule(
+        name="P-bool-arg",
+        parameters=True,
+        pattern=re.compile(r":\s*bool\b"),
+        reason="意味を持つ bool の引数を使わず、2 値の enum にする",
+        doc=_TYPE_RULES,
+    ),
+    _Rule(
+        name="P-resolved-requirement",
+        pattern=re.compile(r"\bResolvedRequirement\b"),
+        reason="呼び出し側が種類を知っている解決は、enum で返さず種類ごとの関数にする",
+        doc="docs/design/07-pricing/01-pricing.md",
+    ),
+    _Rule(
+        name="P-map-partitions",
+        pattern=re.compile(r"\bmap_partitions\s*\("),
+        reason="対象の分割は parallel の共通の関数を通す",
+        doc="docs/design/03-architecture.md#並列と再現性",
+        allow=("crates/rhodium-core/src/parallel.rs",),
+    ),
+)
 
 
 def _parse_changes(raw: bytes) -> dict[str, bool]:
@@ -109,6 +345,12 @@ def _plan(changed: Mapping[str, bool]) -> _CheckPlan:
             if (files := _python_files(directory, changed))
         ),
         crates=frozenset(crate for path in changed if (crate := _crate_for_path(path)) is not None),
+        ui_pairs=any(
+            PurePosixPath(path).full_match(glob)
+            for path in changed
+            for glob in ("crates/*/tests/compile_fail.rs", "crates/*/tests/ui/**")
+        ),
+        examples=any(path.startswith("examples/") for path in changed),
         workspace=bool(
             CARGO_INPUTS & changed.keys() or any(p.startswith(".cargo/") for p in changed)
         ),
@@ -135,8 +377,8 @@ def _git(
     repo: Path,
     args: tuple[str, ...],
     *,
-    env: Mapping[str, str] | None = None,
     data: bytes | None = None,
+    env: Mapping[str, str] | None = None,
     check: bool = True,
 ) -> subprocess.CompletedProcess[bytes]:
     return subprocess.run(  # ruff: ignore[subprocess-without-shell-equals-true]
@@ -192,54 +434,130 @@ def _check_files(root: Path, files: tuple[_FileCheck, ...]) -> None:
         )
 
 
-def _checkout(
-    repo: Path, root: Path, files: tuple[_FileCheck, ...] | None, tip: str | None
-) -> None:
-    command = ("checkout-index", f"--prefix={root}/")
-    env = None
-    if tip is not None:
-        # push の対象が別の参照でも、その commit の内容で検査する。
-        env = {**os.environ, "GIT_INDEX_FILE": str(root.parent / "index")}
-        _git(repo, ("read-tree", tip), env=env)
-    if files is None:
-        _git(repo, (*command, "--all"), env=env)
-    else:
-        # 構文チェックだけなら、対象のファイルだけ取り出す。
-        paths = b"".join(os.fsencode(name) + b"\0" for name, _ in files)
-        _git(repo, (*command, "-z", "--stdin"), env=env, data=paths)
-
-
-def _python_commands(tools: object, files: tuple[str, ...]) -> Iterator[_Command]:
+def _python_commands(tools: object, files: tuple[str, ...], *, tests: bool) -> Iterator[_Command]:
     if isinstance(tools, dict):
         if "ruff" in tools:
             yield ("ruff", "check", "--no-fix", "--", *files)
             yield ("ruff", "format", "--check", "--", *files)
         if "ty" in tools:
             yield ("ty", "check", "--", *files)
+    if tests:
+        # 検査の仕組みそのものの単体テスト。変更したファイルによらず全部を走らせる。
+        yield (sys.executable, "-I", "-m", "unittest", "discover", "-s", ".", "-p", "test_*.py")
 
 
 def _check_python(root: Path, project: _PythonCheck) -> None:
     directory, files = project
     config = root / directory / "pyproject.toml"
     if config.is_file():
-        for command in _python_commands(_read_toml(config).get("tool"), files):
+        tests = any(config.parent.glob("test_*.py"))
+        for command in _python_commands(_read_toml(config).get("tool"), files, tests=tests):
             _run(config.parent, command)
 
 
+@dataclass(frozen=True, slots=True)
+class _Workspace:
+    # crates/ の下のディレクトリ名から package 名への対応。
+    packages: Mapping[str, str]
+    # package 名から、それに依存する workspace の package 名の集合への対応。
+    dependents: Mapping[str, frozenset[str]]
+
+
+def _parse_workspace(metadata: object) -> _Workspace:
+    """`cargo metadata --no-deps` の出力から、workspace の逆依存の表を作る。
+
+    normal、dev、build の依存をすべて辺にする。dev の依存を含めるのは、下流の crate の
+    テストやベンチマークが上流のファイルや振る舞いに依存するため。
+
+    Returns:
+        ディレクトリ名と package 名の対応と、逆依存の表。
+
+    Raises:
+        TypeError: 出力の形が想定と違う場合。
+
+    """
+    if not isinstance(metadata, dict):
+        message = "cargo metadata: expected a JSON object"
+        raise TypeError(message)
+    members = metadata.get("workspace_members")
+    listed = metadata.get("packages")
+    if not isinstance(members, list) or not isinstance(listed, list):
+        message = "cargo metadata: missing workspace members or packages"
+        raise TypeError(message)
+    packages: dict[str, str] = {}
+    dependencies: dict[str, set[str]] = {}
+    for package in listed:
+        if not isinstance(package, dict) or package.get("id") not in members:
+            continue
+        name = package.get("name")
+        manifest = package.get("manifest_path")
+        if not isinstance(name, str) or not isinstance(manifest, str):
+            message = "cargo metadata: package without a name or manifest path"
+            raise TypeError(message)
+        used = package.get("dependencies")
+        if not isinstance(used, list):
+            message = f"cargo metadata: {name} has no dependency list"
+            raise TypeError(message)
+        packages[Path(manifest).parent.name] = name
+        dependencies[name] = {
+            target
+            for dependency in used
+            if isinstance(dependency, dict) and isinstance(target := dependency.get("name"), str)
+        }
+    names = frozenset(dependencies)
+    dependents = {
+        name: frozenset(user for user, used in dependencies.items() if name in used)
+        for name in names
+    }
+    return _Workspace(packages=MappingProxyType(packages), dependents=MappingProxyType(dependents))
+
+
+def _dependent_closure(workspace: _Workspace, seeds: Iterable[str]) -> frozenset[str]:
+    """起点の package と、それに直接または間接に依存する package の集合を返す。
+
+    Returns:
+        逆依存の閉包。
+
+    """
+    seen: set[str] = set()
+    pending = [name for name in seeds if name in workspace.dependents]
+    while pending:
+        name = pending.pop()
+        if name not in seen:
+            seen.add(name)
+            pending.extend(workspace.dependents[name])
+    return frozenset(seen)
+
+
+def _scope_arguments(workspace: _Workspace, plan: _CheckPlan) -> tuple[str, ...]:
+    """変更した crate の逆依存の閉包を、Cargo の package 指定にする。
+
+    Returns:
+        `--workspace` か `-p <package>` の並び。検査する crate が無ければ空。
+
+    """
+    seeds = [workspace.packages[crate] for crate in plan.crates if crate in workspace.packages]
+    if plan.examples:
+        seeds.extend(EXAMPLE_READERS)
+    packages = _dependent_closure(workspace, seeds)
+    if packages and packages == workspace.dependents.keys():
+        return ("--workspace",)
+    return tuple(arg for name in sorted(packages) for arg in ("-p", name))
+
+
 def _cargo_scope(root: Path, plan: _CheckPlan) -> tuple[str, ...]:
     if plan.workspace:
         return ("--workspace",)
-    manifests = (root / "crates" / crate / "Cargo.toml" for crate in sorted(plan.crates))
-    return tuple(arg for path in manifests if path.is_file() for arg in ("-p", _package_name(path)))
-
-
-def _package_name(path: Path) -> str:
-    package = _read_toml(path).get("package")
-    name = package.get("name") if isinstance(package, dict) else None
-    if not isinstance(name, str):
-        message = f"{path}: package.name is missing or is not a string"
-        raise TypeError(message)
-    return name
+    if not plan.crates and not plan.examples:
+        return ()
+    metadata = subprocess.run(
+        # 開発環境の PATH にある Cargo を使う。依存の解決は要らないので、ネットワークを使わない。
+        ("cargo", "metadata", "--format-version", "1", "--no-deps", "--offline"),  # ruff: ignore[start-process-with-partial-path]
+        cwd=root,
+        stdout=subprocess.PIPE,
+        check=True,
+    ).stdout
+    return _scope_arguments(_parse_workspace(json.loads(metadata)), plan)
 
 
 def _cargo_commands(scope: tuple[str, ...]) -> tuple[_Command, ...]:
@@ -251,65 +569,378 @@ def _cargo_commands(scope: tuple[str, ...]) -> tuple[_Command, ...]:
     )
 
 
-def _cargo_environment(target: Path, env: Mapping[str, str]) -> dict[str, str]:
-    return {
-        **env,
-        "CARGO_TARGET_DIR": os.path.normpath(target),
-        "RUSTDOCFLAGS": (env.get("RUSTDOCFLAGS", "") + " -D warnings").strip(),
-    }
-
-
-def _repository_target(repo: Path, env: Mapping[str, str]) -> Path:
-    return repo / env.get("CARGO_TARGET_DIR", "target")
-
-
-def _check_push(root: Path, plan: _CheckPlan, target: Path) -> None:
+def _check_push(root: Path, plan: _CheckPlan) -> None:
     for project in plan.python:
         _check_python(root, project)
     scope = _cargo_scope(root, plan)
+    env = {
+        **os.environ,
+        "RUSTDOCFLAGS": (os.environ.get("RUSTDOCFLAGS", "") + " -D warnings").strip(),
+    }
     if scope:
-        env = _cargo_environment(target, os.environ)
         for command in _cargo_commands(scope):
             _run(root, command, env=env)
+    if plan.ui_pairs and "--workspace" not in scope and "rhodium-testkit" not in scope:
+        _run(
+            root,
+            ("cargo", "test", "--locked", "-p", "rhodium-testkit", "--test", "compile_fail_pairs"),
+            env=env,
+        )
 
 
-def _clean_at_tip(repo: Path, tip: str) -> bool:
+@dataclass(frozen=True, slots=True)
+class _Violation:
+    rule: _Rule
+    path: str
+    line: int
+    text: str
+
+    def render(self) -> str:
+        """違反を、場所、規則、理由、定義元の一行と、該当する行にする。
+
+        Returns:
+            表示する文字列。
+
+        """
+        return (
+            f"{self.path}:{self.line}: [{self.rule.name}] {self.rule.reason} ({self.rule.doc})\n"
+            f"    {self.text.strip()}"
+        )
+
+
+def _in_scope(path: str, scope: _Scope) -> bool:
+    match PurePosixPath(path).parts:
+        case ("crates", _, "src", *rest) if path.endswith(".rs"):
+            # src/ir/tests.rs と src/ir/tests/ の下は単体テストのファイル。
+            return scope == "all" or not any(PurePosixPath(part).stem == "tests" for part in rest)
+        case ("crates", _, _, *_) if path.endswith(".rs"):
+            return scope == "all"
+        case _:
+            return False
+
+
+def _applies(rule: _Rule, path: str) -> bool:
+    location = PurePosixPath(path)
     return (
-        _git(repo, ("rev-parse", "HEAD")).stdout.strip().decode() == tip
-        and not _git(
-            repo,
-            (
-                "status",
-                "--porcelain=v1",
-                "-z",
-                "--untracked-files=normal",
-                "--ignore-submodules=none",
-            ),
-        ).stdout
+        _in_scope(path, rule.scope)
+        and (not rule.paths or any(location.full_match(glob) for glob in rule.paths))
+        and not any(location.full_match(glob) for glob in rule.allow)
     )
 
 
-def _validate(repo: Path, revisions: tuple[str, ...], tip: str | None = None) -> None:
-    _run(repo, ("git", "diff", "--check", *revisions, "--"))
+# リテラルの境界を先に認識し、内部の括弧とコメント記号を構文として数えない。
+_RUST_TOKEN = re.compile(
+    r"""//[^\n]*|/\*|(?:br|cr|r)(?P<hash>\#*)".*?"(?P=hash)|"(?:\\.|[^"\\])*"|'(?:\\(?:u\{[^}]*\}|x[\da-fA-F]{2}|.)|[^'\\\n])'""",
+    re.DOTALL,
+)
+_BLOCK_MARKER = re.compile(r"/\*|\*/")
+
+
+def _source_views(source: str) -> tuple[str, str]:
+    # コード照合用には文字列を残す。構文用にはリテラルも隠す。双方の行と列を保つ。
+    code: list[str] = []
+    syntax: list[str] = []
+    end = 0
+    while match := _RUST_TOKEN.search(source, end):
+        prefix = source[end : match.start()]
+        code.append(prefix)
+        syntax.append(prefix)
+        end = match.end()
+        if match[0] == "/*":
+            depth = 1
+            for marker in _BLOCK_MARKER.finditer(source, end):
+                depth += 1 if marker[0] == "/*" else -1
+                if depth == 0:
+                    end = marker.end()
+                    break
+            else:
+                end = len(source)
+        token = source[match.start() : end]
+        masked = re.sub(r"[^\n]+", lambda part: " " * len(part[0]), token)
+        syntax.append(masked)
+        code.append(masked if token.startswith(("//", "/*")) else token)
+    return "".join(code) + source[end:], "".join(syntax) + source[end:]
+
+
+_PARAM_TOKEN = re.compile(r"->|[<>()\[\]{};]|:\s*\bbool\b")
+
+
+def _bool_arguments(syntax: str, pattern: re.Pattern[str]) -> set[int]:
+    # 関数宣言の型引数と引数の括弧をたどる。関数ポインタなど入れ子の型は別の引数として数えない。
+    numbers: set[int] = set()
+    for function in re.finditer(r"\bfn\s+\w+\s*", syntax):
+        depth = 0
+        opened = False
+        for token in _PARAM_TOKEN.finditer(syntax, function.end()):
+            value = token[0]
+            if opened and depth == 1 and pattern.fullmatch(value):
+                numbers.add(syntax.count("\n", 0, token.end()) + 1)
+            if depth == 0 and value == "(":
+                opened = True
+            if value in {"<", "(", "[", "{"}:
+                depth += 1
+            elif value in {">", ")", "]", "}"}:
+                depth -= 1
+            if depth == 0 and (opened or value == ";"):
+                break
+    return numbers
+
+
+def _item_end(lines: list[str], start: int) -> int:
+    """`start` 行目から始まる item が終わる行の添字を返す。
+
+    波括弧の対応を数える。文字列とコメントの中の括弧は数えない。波括弧を持たない item は、
+    `;` で終わる行で終わる。
+
+    Returns:
+        item の最後の行の添字。
+
+    """
+    depth = 0
+    opened = False
+    for index in range(start, len(lines)):
+        code = lines[index]
+        opened |= "{" in code
+        depth += code.count("{") - code.count("}")
+        if (opened and depth <= 0) or (not opened and code.rstrip().endswith(";")):
+            return index
+    return len(lines) - 1
+
+
+def _test_lines(lines: list[str]) -> frozenset[int]:
+    """`#[cfg(test)]` を付けた item の行の添字を返す。
+
+    Returns:
+        テストのコードとして照合から外す行の添字の集合。
+
+    """
+    excluded: set[int] = set()
+    index = 0
+    while index < len(lines):
+        stripped = lines[index].strip()
+        if stripped.startswith("#[cfg(test)]"):
+            # 同じ行に item が続くときは、その行から item が始まる。
+            same_line = bool(stripped.removeprefix("#[cfg(test)]").strip())
+            item = index if same_line else index + 1
+            while (
+                not same_line
+                and item < len(lines)
+                and (not lines[item].strip() or lines[item].strip().startswith("#["))
+            ):
+                item += 1
+            end = _item_end(lines, min(item, len(lines) - 1))
+            excluded.update(range(index, end + 1))
+            index = end
+        index += 1
+    return frozenset(excluded)
+
+
+def _violations(
+    rules: Iterable[_Rule],
+    files: Mapping[str, str],
+    added: Mapping[str, frozenset[int]] | None,
+) -> list[_Violation]:
+    """規則に当たる行を集める。
+
+    Args:
+        rules: 照合する規則。
+        files: パスから、検査する版のファイルの内容への対応。
+        added: パスから、追加した行の番号 (1 始まり) への対応。None なら全規則をツリー
+            全体で照合する。監査に使う。
+
+    Returns:
+        違反の一覧。パスと行の順。
+
+    """
+    rules = tuple(rules)
+    found: list[_Violation] = []
+    for path in sorted(files):
+        applicable = [rule for rule in rules if _applies(rule, path)]
+        if not applicable:
+            continue
+        lines = files[path].split("\n")
+        code, syntax = _source_views(files[path])
+        code_lines = code.split("\n")
+        tests = _test_lines(syntax.split("\n"))
+        for rule in applicable:
+            if added is None or rule.mode == "tree":
+                numbers: Iterable[int] = range(1, len(lines) + 1)
+            else:
+                numbers = sorted(added.get(path, frozenset()))
+            targets = lines if rule.comments else code_lines
+            if rule.parameters:
+                numbers = _bool_arguments(syntax, rule.pattern).intersection(numbers)
+            else:
+                numbers = (
+                    n for n in numbers if n <= len(lines) and rule.pattern.search(targets[n - 1])
+                )
+            for number in numbers:
+                index = number - 1
+                if rule.scope == "library" and index in tests:
+                    continue
+                text = lines[index]
+                if rule.exempt is not None and rule.exempt.search(targets[index]):
+                    continue
+                window = lines[max(0, index - rule.window) : index + 1]
+                if rule.requires_nearby is not None and any(
+                    rule.requires_nearby in line for line in window
+                ):
+                    continue
+                found.append(_Violation(rule=rule, path=path, line=number, text=text))
+    return sorted(found, key=lambda violation: (violation.path, violation.line))
+
+
+def _added_lines(diff: str) -> dict[str, frozenset[int]]:
+    """`git diff -U0` の出力から、変更後のファイルで追加した行の番号を読む。
+
+    Returns:
+        パスから、追加した行の番号 (1 始まり) の集合への対応。
+
+    """
+    added: dict[str, set[int]] = {}
+    path: str | None = None
+    number = 0
+    for line in diff.splitlines():
+        if line.startswith("+++ "):
+            target = line[4:].removesuffix("\t")
+            path = target.removeprefix("b/") if target != "/dev/null" else None
+        elif line.startswith("@@ "):
+            # @@ -a,b +c,d @@ の c が、変更後のファイルの最初の行。
+            new = line.split(" ")[2]
+            number = int(new.removeprefix("+").split(",")[0])
+        elif line.startswith("+") and path is not None:
+            added.setdefault(path, set()).add(number)
+            number += 1
+    return {name: frozenset(numbers) for name, numbers in added.items()}
+
+
+def _index_sources(repo: Path) -> dict[str, str]:
+    # index の行は「mode blob stage\tpath」。
+    listing = _git(repo, ("ls-files", "-s", "-z", "--", "crates/")).stdout
+    entries = [entry.split(b"\t", 1) for entry in listing.split(b"\0") if entry]
+    blobs = [(header.split()[0], header.split()[1], path) for header, path in entries]
+    wanted = [
+        (blob, os.fsdecode(path))
+        for mode, blob, path in blobs
+        if mode in {b"100644", b"100755"} and path.endswith(b".rs")
+    ]
+    if not wanted:
+        return {}
+    output = _git(repo, ("cat-file", "--batch"), data=b"".join(blob + b"\n" for blob, _ in wanted))
+    files: dict[str, str] = {}
+    offset = 0
+    data = output.stdout
+    for _, path in wanted:
+        header_end = data.index(b"\n", offset)
+        size = int(data[offset:header_end].split()[2])
+        start = header_end + 1
+        files[path] = data[start : start + size].decode("utf-8", errors="replace")
+        offset = start + size + 1
+    return files
+
+
+def _check_patterns(
+    repo: Path,
+    revisions: tuple[str, ...],
+    files: Mapping[str, str],
+    paths: tuple[str, ...],
+    *,
+    env: Mapping[str, str] | None = None,
+) -> None:
+    diff = _git(
+        repo,
+        (
+            "-c",
+            "core.quotePath=false",
+            "--literal-pathspecs",
+            "diff",
+            "-U0",
+            "--inter-hunk-context=0",
+            "--no-color",
+            "--no-renames",
+            "--no-ext-diff",
+            "--src-prefix=a/",
+            "--dst-prefix=b/",
+            *revisions,
+            "--",
+            *paths,
+        ),
+        env=env,
+    ).stdout.decode("utf-8", errors="replace")
+    found = _violations(RULES, files, _added_lines(diff))
+    if found:
+        message = "forbidden patterns:\n" + "\n".join(violation.render() for violation in found)
+        raise ValueError(message)
+
+
+def _worktree_sources(repo: Path) -> dict[str, str]:
+    listing = _git(
+        repo, ("ls-files", "--cached", "--others", "--exclude-standard", "-z", "--", "crates/")
+    ).stdout
+    return {
+        name: (repo / name).read_text(encoding="utf-8", errors="replace")
+        for raw in listing.split(b"\0")
+        if raw.endswith(b".rs")
+        and (repo / (name := os.fsdecode(raw))).is_file()
+        and not (repo / name).is_symlink()
+    }
+
+
+def _audit(repo: Path) -> str:
+    """作業ツリーの全ファイルを、全規則のツリー全体で照合した結果をまとめる。
+
+    Returns:
+        規則ごとの件数と場所の一覧。
+
+    """
+    found = _violations(RULES, _worktree_sources(repo), None)
+    lines = [
+        f"{rule.name} ({rule.mode}): {sum(item.rule is rule for item in found)}" for rule in RULES
+    ]
+    lines.extend(violation.render() for violation in found)
+    return "\n".join(lines) + "\n"
+
+
+def _changes(repo: Path, revisions: tuple[str, ...]) -> dict[str, bool]:
     raw = _git(repo, ("diff", "--raw", "--no-renames", "-z", *revisions, "--")).stdout
-    plan = _plan(_parse_changes(raw))
-    full_tip = tip if plan.needs_push_checks() else None
-    if full_tip is not None and _clean_at_tip(repo, full_tip):
-        # 通常の push は同じソースパスを使い、Cargo の増分ビルドを活かす。
-        _check_files(repo, plan.files)
-        _check_push(repo, plan, _repository_target(repo, os.environ))
-    elif plan.files or full_tip is not None:
-        with TemporaryDirectory(prefix="rhodium-hook-") as tmp:
-            root = Path(tmp) / "src"
-            root.mkdir()
-            _checkout(repo, root, plan.files if full_tip is None else None, tip)
-            _check_files(root, plan.files)
-            if full_tip is not None:
-                # 一時 checkout のテストは env!("CARGO_MANIFEST_DIR") などで一時パスを埋め込む。
-                # リポジトリの target を共有すると、Cargo が新鮮と判断したその成果物を後の通常の
-                # cargo test が再利用し、削除済みの一時パスを読んで失敗する。ビルド先を一時
-                # ディレクトリの中に分け、checkout と一緒に消す。依存も作り直すので遅くなる。
-                _check_push(root, plan, Path(tmp) / "target")
+    return _parse_changes(raw)
+
+
+def _pre_commit(repo: Path) -> None:
+    revisions = ("--cached",)
+    _run(repo, ("git", "diff", "--check", *revisions, "--"))
+    _check_patterns(repo, revisions, _index_sources(repo), ("crates/",))
+    files = _plan(_changes(repo, revisions)).files
+    if files:
+        with TemporaryDirectory(prefix="rhodium-index-") as tmp:
+            paths = b"".join(os.fsencode(name) + b"\0" for name, _ in files)
+            _git(repo, ("checkout-index", f"--prefix={tmp}/", "-z", "--stdin"), data=paths)
+            _check_files(Path(tmp), files)
+
+
+def _pre_push(repo: Path, base: str, tip: str) -> None:
+    # 対象は送る差分から選び、検査する内容と行番号はともに作業ツリーから取る。
+    changed = _changes(repo, (base, tip))
+    if not changed:
+        return
+    paths = tuple(changed)
+    regular = {name: (repo / name).is_file() and not (repo / name).is_symlink() for name in paths}
+    # 別の ref にあるパスは、現在の index に無くても作業ツリーに存在しうる。
+    # 比較対象の名前だけを一時 index に載せる。ファイルのコピーとビルド先の変更は不要。
+    with TemporaryDirectory(prefix="rhodium-diff-") as tmp:
+        env = {**os.environ, "GIT_INDEX_FILE": str(Path(tmp) / "index")}
+        empty = _git(repo, ("hash-object", "--stdin"), data=b"").stdout.strip()
+        entries = b"".join(
+            b"100644 " + empty + b"\t" + os.fsencode(name) + b"\0"
+            for name, present in regular.items()
+            if present
+        )
+        _git(repo, ("update-index", "--info-only", "-z", "--index-info"), env=env, data=entries)
+        _run(repo, ("git", "--literal-pathspecs", "diff", "--check", base, "--", *paths), env=env)
+        _check_patterns(repo, (base,), _worktree_sources(repo), paths, env=env)
+    plan = _plan(regular)
+    _check_files(repo, plan.files)
+    _check_push(repo, plan)
 
 
 def _push_base(repo: Path, remote: str, tip: str, old: str) -> str | None:
@@ -345,12 +976,15 @@ def _main() -> None:
     repo = Path.cwd()
     match sys.argv[1:]:
         case ["pre-commit"]:
-            _validate(repo, ("--cached",))
+            _pre_commit(repo)
         case ["pre-push", remote]:
             for base, tip in dict.fromkeys(_push_ranges(repo, remote, sys.stdin)):
-                _validate(repo, (base, tip), tip)
+                _pre_push(repo, base, tip)
+        case ["audit"]:
+            # 既存の違反の残りを数える。失敗にはせず、担当の issue の進み具合を見るのに使う。
+            sys.stdout.write(_audit(repo))
         case _:
-            sys.exit("usage: checks.py pre-commit | pre-push <remote>")
+            sys.exit("usage: checks.py pre-commit | pre-push <remote> | audit")
 
 
 if __name__ == "__main__":
diff --git a/.githooks/pyproject.toml b/.githooks/pyproject.toml
index 894efa79745b4ed67ceada445e76df273ec274aa..03b115e0b3a66ee9a5da798e50168bf9c3ab2cd4 100644
--- a/.githooks/pyproject.toml
+++ b/.githooks/pyproject.toml
@@ -30,6 +30,18 @@ ignore = [
     "missing-terminal-punctuation",
 ]
 
+[tool.ruff.lint.per-file-ignores]
+# 単体テストは unittest で走らせ、判定は assert で書き、checks.py の非公開の関数を直接確かめる。
+"test_*.py" = [
+    # 標準の unittest で例外の種類とメッセージを検査する。
+    "pytest-unittest-raises-assertion",
+    "assert",
+    "private-member-access",
+    # unittest の test メソッドは self を使わなくてもクラスに置き、期待値は表にその場で書く。
+    "no-self-use",
+    "magic-value-comparison",
+]
+
 [tool.ruff.lint.pycodestyle]
 max-line-length = 120
 
diff --git a/.githooks/test_checks.py b/.githooks/test_checks.py
new file mode 100644
index 0000000000000000000000000000000000000000..e95803aa7090de434a99d10a82b1253919cdd00b
--- /dev/null
+++ b/.githooks/test_checks.py
@@ -0,0 +1,332 @@
+"""checks.py の検査範囲の計算と、禁止パターンの照合を確かめる。"""
+
+import dataclasses
+import re
+import unittest
+from pathlib import Path
+from tempfile import TemporaryDirectory
+from unittest.mock import patch
+
+import checks
+
+# 依存の向き: core <- products <- pricing <- rhodium <- cli。testkit はすべての dev 依存。
+_METADATA = {
+    "workspace_members": ["core", "testkit", "products", "pricing", "rhodium", "cli"],
+    "packages": [
+        {
+            "id": "core",
+            "name": "rhodium-core",
+            "manifest_path": "/w/crates/rhodium-core/Cargo.toml",
+            "dependencies": [{"name": "rhodium-testkit"}, {"name": "serde"}],
+        },
+        {
+            "id": "testkit",
+            "name": "rhodium-testkit",
+            "manifest_path": "/w/crates/rhodium-testkit/Cargo.toml",
+            "dependencies": [],
+        },
+        {
+            "id": "products",
+            "name": "rhodium-products",
+            "manifest_path": "/w/crates/rhodium-products/Cargo.toml",
+            "dependencies": [{"name": "rhodium-core"}],
+        },
+        {
+            "id": "pricing",
+            "name": "rhodium-pricing",
+            "manifest_path": "/w/crates/rhodium-pricing/Cargo.toml",
+            "dependencies": [{"name": "rhodium-products"}],
+        },
+        {
+            "id": "rhodium",
+            "name": "rhodium",
+            "manifest_path": "/w/crates/rhodium/Cargo.toml",
+            "dependencies": [{"name": "rhodium-pricing"}],
+        },
+        {
+            "id": "cli",
+            "name": "rhodium-cli",
+            "manifest_path": "/w/crates/rhodium-cli/Cargo.toml",
+            "dependencies": [{"name": "rhodium"}],
+        },
+        # workspace の外の package は無視する。
+        {
+            "id": "serde",
+            "name": "serde",
+            "manifest_path": "/registry/serde/Cargo.toml",
+            "dependencies": [],
+        },
+    ],
+}
+
+
+def _plan(crates: frozenset[str], *, examples: bool = False) -> checks._CheckPlan:
+    return checks._CheckPlan(files=(), python=(), crates=crates, examples=examples, workspace=False)
+
+
+class CargoScopeTest(unittest.TestCase):
+    """変更した crate の逆依存の閉包を検査範囲にすることを確かめる。"""
+
+    workspace = checks._parse_workspace(_METADATA)
+
+    def test_core_change_reaches_every_downstream_crate(self) -> None:
+        """Core を変えると、それに依存する crate をすべて検査する。testkit は core に依存しない。"""
+        scope = checks._scope_arguments(self.workspace, _plan(frozenset({"rhodium-core"})))
+        assert scope[1::2] == (
+            "rhodium",
+            "rhodium-cli",
+            "rhodium-core",
+            "rhodium-pricing",
+            "rhodium-products",
+        )
+
+    def test_dev_dependency_change_reaches_the_whole_workspace(self) -> None:
+        """開発用の依存の testkit を変えると全 crate が閉包に入るので、--workspace にする。"""
+        scope = checks._scope_arguments(self.workspace, _plan(frozenset({"rhodium-testkit"})))
+        assert scope == ("--workspace",)
+
+    def test_examples_change_checks_the_facade_and_the_cli(self) -> None:
+        """examples/ のサンプルは facade のテストが読み、CLI は facade に依存する。"""
+        scope = checks._scope_arguments(self.workspace, _plan(frozenset(), examples=True))
+        assert scope == ("-p", "rhodium", "-p", "rhodium-cli")
+
+    def test_ui_registration_includes_deletions(self) -> None:
+        """UI ケースや登録を取り除いたときも testkit の対応表を検査する。"""
+        for path in ["crates/c/tests/ui/a.json", "crates/c/tests/compile_fail.rs"]:
+            assert checks._plan({path: False}).ui_pairs
+        assert not checks._plan({"crates/c/tests/other.rs": True}).ui_pairs
+
+    def test_ui_registration_runs_once_outside_the_reverse_dependencies(self) -> None:
+        """通常の検査で testkit が範囲に入らない場合だけ、対応表を検査する。"""
+        plan = checks._plan({"crates/rhodium-core/tests/ui/case.json": True})
+        root = Path()
+        for scope in [("-p", "rhodium-core"), ("--workspace",)]:
+            with (
+                patch.object(checks, "_cargo_scope", return_value=scope),
+                patch.object(checks, "_run") as run,
+            ):
+                checks._check_push(root, plan)
+                commands = [call.args[1] for call in run.call_args_list]
+            expected = list(checks._cargo_commands(scope))
+            if scope != ("--workspace",):
+                expected.append(
+                    (
+                        "cargo",
+                        "test",
+                        "--locked",
+                        "-p",
+                        "rhodium-testkit",
+                        "--test",
+                        "compile_fail_pairs",
+                    )
+                )
+            assert commands == expected
+
+    def test_directories_without_a_package_are_ignored(self) -> None:
+        """Cargo.toml の無い crates/ の下の残骸は、検査範囲に入れない。"""
+        scope = checks._scope_arguments(self.workspace, _plan(frozenset({"rhodium-refdata"})))
+        assert scope == ()
+
+
+def _rule(**changes: object) -> checks._Rule:
+    # 照合の仕組みだけを確かめるため、単純な語の規則を作る。
+    base = checks._Rule(
+        name="P-test",
+        pattern=re.compile(r"\bforbidden\b"),
+        reason="テストの規則",
+        doc="docs/design/11-quality.md#ci",
+    )
+    return dataclasses.replace(base, **changes)
+
+
+_LIB = "crates/c/src/lib.rs"
+
+
+class HookInputsTest(unittest.TestCase):
+    """一時リポジトリで、実際の差分と読み出す版が一致することを確かめる。"""
+
+    def test_push_uses_worktree_and_commit_uses_index(self) -> None:
+        """部分ステージの index と作業ツリーを分け、文脈を結合する Git 設定でも行を特定する。"""
+        with TemporaryDirectory() as tmp:
+            repo = Path(tmp)
+            path = repo / _LIB
+            path.parent.mkdir(parents=True)
+            checks._git(repo, ("init", "-q"))
+            checks._git(repo, ("config", "diff.interHunkContext", "2"))
+            path.write_text("// base\n\nfn f() {}\n")
+            checks._git(repo, ("add", "."))
+            checks._git(
+                repo,
+                (
+                    # 開発者の global 設定の署名やフックに左右されずに、基準の commit を作る。
+                    "-c",
+                    "user.name=test",
+                    "-c",
+                    "user.email=test@example.org",
+                    "-c",
+                    "core.hooksPath=/dev/null",
+                    "-c",
+                    "commit.gpgsign=false",
+                    "commit",
+                    "-qm",
+                    "base",
+                ),
+            )
+            path.write_text("// changed\n\nfn f() { let _ = 1e-9; }\n")
+            checks._git(repo, ("add", "."))
+            tree = checks._git(repo, ("write-tree",)).stdout.strip().decode()
+            path.write_text("// changed\n\nfn f() { let _ = 1; }\n")
+            with self.assertRaisesRegex(ValueError, r"lib\.rs:3: \[P-tolerance\]"):
+                checks._pre_commit(repo)
+            with (
+                patch.object(checks, "_check_push") as push,
+                patch.object(checks, "_check_files") as files,
+            ):
+                checks._pre_push(repo, "HEAD", tree)
+                assert push.call_args.args[0] == repo
+                assert files.call_args.args == (repo, ((_LIB, "rust"),))
+            path.write_text("// changed\n\nfn f() { let _ = 1e-8; }\n")
+            with self.assertRaisesRegex(ValueError, r"lib\.rs:3: \[P-tolerance\]"):
+                checks._pre_push(repo, "HEAD", tree)
+
+            # 別ブランチで追加されたファイルが現在の index に無くても、作業ツリーの行を検査する。
+            new = path.with_name("new.rs")
+            new.write_text("// other branch\n")
+            checks._git(repo, ("add", str(new)))
+            other = checks._git(repo, ("write-tree",)).stdout.strip().decode()
+            checks._git(repo, ("update-index", "--force-remove", "--", str(new)))
+            path.write_text("// base\n\nfn f() {}\n")
+            new.write_text("fn g() { let _ = 1e-8; }\n")
+            with self.assertRaisesRegex(ValueError, r"new\.rs:1: \[P-tolerance\]"):
+                checks._pre_push(repo, "HEAD", other)
+
+
+class PatternTest(unittest.TestCase):
+    """禁止パターンの照合の範囲とモードを確かめる。"""
+
+    source = "fn a() {}\nlet forbidden = 1;\nlet forbidden = 2;\n"
+
+    def test_diff_mode_reports_only_added_lines(self) -> None:
+        """既存の違反 (2 行目) は許し、追加した 3 行目だけを報告する。"""
+        found = checks._violations([_rule()], {_LIB: self.source}, {_LIB: frozenset({3})})
+        assert [item.line for item in found] == [3]
+
+    def test_tree_mode_reports_the_whole_file(self) -> None:
+        """Tree の規則は、追加した行によらずファイル全体を照合する。"""
+        found = checks._violations([_rule(mode="tree")], {_LIB: self.source}, {})
+        assert [item.line for item in found] == [2, 3]
+
+    def test_audit_reports_every_rule_over_the_whole_tree(self) -> None:
+        """監査 (added が None) では、diff の規則もツリー全体で照合する。"""
+        found = checks._violations([_rule()], {_LIB: self.source}, None)
+        assert len(found) == 2
+
+    def test_nearby_marker_allows_the_line(self) -> None:
+        """直前の window 行に注記があれば違反にしない。窓の外の注記は数えない。"""
+        source = "// 不変条件: 作った場所で保証する\nlet x = 0;\nforbidden();\n\n\n\nforbidden();\n"
+        rule = _rule(requires_nearby="不変条件:", window=3)
+        found = checks._violations([rule], {_LIB: source}, None)
+        assert [item.line for item in found] == [7]
+
+    def test_library_scope_skips_test_code(self) -> None:
+        """ライブラリの規則は、テストのファイルと #[cfg(test)] の item を照合しない。"""
+        source = (
+            "fn a() {}\n"
+            "#[cfg(test)]\n"
+            "mod tests {\n"
+            '    const S: &str = "}";\n'
+            "    fn b() { forbidden(); }\n"
+            "}\n"
+            "fn c() { forbidden(); }\n"
+            "#[cfg(test)] extern crate self as c;\n"
+            "fn d() { forbidden(); }\n"
+        )
+        files = {
+            _LIB: source,
+            "crates/c/src/ir/tests.rs": "forbidden();\n",
+            "crates/c/src/ir/tests/standard.rs": "forbidden();\n",
+            "crates/c/tests/it.rs": "forbidden();\n",
+        }
+        found = checks._violations([_rule()], files, None)
+        assert [(item.path, item.line) for item in found] == [(_LIB, 7), (_LIB, 9)]
+
+    def test_rust_comments_and_multiline_arguments(self) -> None:
+        """コメント、リテラル、入れ子の型を区別し、bool の行だけを報告する。"""
+        source = (
+            'const S: &str = r#"hello " /* \n // world"#;\n'
+            "/* outer /* nested */ 1e-10\nfn ignored(x: bool) {} */\n"
+            "fn check<T: Fn()>(\n  x: bool,\n  f: fn(bool),\n  y: /* inline */ bool,\n) {}\n"
+        )
+        found = checks._violations(checks.RULES, {_LIB: source}, None)
+        assert [(item.rule.name, item.line) for item in found] == [
+            ("P-bool-arg", 6),
+            ("P-bool-arg", 8),
+        ]
+
+    def test_all_scope_includes_tests(self) -> None:
+        """All の規則は、テストのファイルも照合する。"""
+        files = {"crates/c/tests/it.rs": "forbidden\n"}
+        found = checks._violations([_rule(scope="all")], files, None)
+        assert len(found) == 1
+
+    def test_comments_paths_allow_and_exempt(self) -> None:
+        """コメントの中、許したファイル、対象外のファイル、exempt の形の行は違反にしない。"""
+        source = "// forbidden\nconst X: u8 = forbidden;\nlet y = forbidden; // ok\n"
+        rule = _rule(
+            exempt=re.compile(r"^\s*const\s"),
+            allow=("crates/allowed/**",),
+            paths=("crates/c/**", "crates/allowed/**"),
+        )
+        files = {
+            _LIB: source,
+            "crates/allowed/src/lib.rs": "forbidden\n",
+            "crates/other/src/lib.rs": "forbidden\n",
+        }
+        found = checks._violations([rule], files, None)
+        assert [(item.path, item.line) for item in found] == [(_LIB, 3)]
+
+    def test_added_lines_follow_hunk_headers(self) -> None:
+        """-U0 の hunk の見出しから、変更後のファイルの行番号を数える。削除したファイルは除く。"""
+        diff = (
+            "diff --git a/crates/c/src/lib.rs b/crates/c/src/lib.rs\n"
+            "--- a/crates/c/src/lib.rs\n"
+            "+++ b/crates/c/src/lib.rs\n"
+            "@@ -1,0 +2,2 @@\n"
+            "+one\n"
+            "+two\n"
+            "@@ -9 +10 @@\n"
+            "-old\n"
+            "+new\n"
+            "--- a/crates/c/src/gone.rs\n"
+            "+++ /dev/null\n"
+            "@@ -1 +0,0 @@\n"
+            "-gone\n"
+        )
+        assert checks._added_lines(diff) == {_LIB: frozenset({2, 3, 10})}
+
+    def test_repository_rules_have_unique_names_and_documents(self) -> None:
+        """規則の名前は表示と PR での参照に使うので一意で、定義元の文書を持つ。"""
+        names = [rule.name for rule in checks.RULES]
+        assert len(names) == len(set(names))
+        assert all(rule.doc.startswith("docs/design/") for rule in checks.RULES)
+
+    def test_internal_rule_ignores_patterns_in_match_arms(self) -> None:
+        """Error::Internal(_) の match の腕は構築ではないので照合しない。"""
+        rule = next(rule for rule in checks.RULES if rule.name == "P-internal")
+        source = 'Error::Internal(_) => 1,\nErr(Error::Internal("x".into()))\n'
+        found = checks._violations([rule], {_LIB: source}, None)
+        assert [item.line for item in found] == [2]
+
+    def test_doctest_rule_rejects_compile_fail_fences(self) -> None:
+        """文書の例の compile_fail の囲みを拒否し、compile_fail_pair! の呼び出しは許す。"""
+        rules = [rule for rule in checks.RULES if rule.name.startswith(("P-doctest", "P-trybuild"))]
+        files = {
+            _LIB: "/// ```compile_fail\n/// ```\n",
+            "crates/c/tests/compile_fail.rs": 'compile_fail_pair!("case");\n',
+        }
+        found = checks._violations(rules, files, None)
+        assert [(item.rule.name, item.line) for item in found] == [("P-doctest-compile-fail", 1)]
+
+
+if __name__ == "__main__":
+    unittest.main()
diff --git a/AGENTS.md b/AGENTS.md
index e643ade3751eb62e5b43f53d54122f3eff6470cb..8da52de6c06b03e9b27d2351698d73e5df8dc5b0 100644
--- a/AGENTS.md
+++ b/AGENTS.md
@@ -29,7 +29,7 @@ rhodium は、金融商品を小さなコア IR で表し、価格・quote 感
 - 営業日・日数計算には現地の `Date`、評価・fixing の公表・行使・担保決済には UTC の `Instant` を用いる。時刻と地域の解決は[日付とカレンダー](docs/design/04-foundation/01-dates-and-calendars.md#瞬間と現地時刻)に従い、評価に使うスナップショットは呼び出し側から受け取る。
 - 商品の意味は[コア IR](docs/design/05-products/01-ir.md)に従う。標準商品は型付きの構造体として、独自の商品は DSL で表現する。組み込みのモデルやプライサーは閉じた enum に足し、利用者の拡張は trait を実装した値をジェネリックな API に渡す。
 - 性能を著しく落とす簡素化をしない。[性能の規則](docs/design/03-architecture.md#性能の規則)の仕組み（準備の再利用、影響範囲だけの再校正、AAD の active set、パスのバッチ評価など）を保つ。
-- 変換・抽出・集約には、`for` ループよりイテレータを優先する。列挙型の分岐には、`if` を重ねるより、各ケースを網羅する `match` を使う。単純な条件判定には `if` を、逐次処理の流れが明確になる場面では `for` を使う。
+- 変換・抽出・集約には、`for` ループよりイテレータを優先する。単純な条件判定には `if` を、逐次処理の流れが明確になる場面では `for` を使う。列挙型の分岐は、下の「設計の自己点検」の S5 に従う。
 - doc に必要なコメントだけでなく、必要十分に、コードの意図や意味が伝わるコメントを書く。コメントは基本的に日本語だが、技術的な単語やコード中の変数名などをコメントで言及するときは英語でいい。テストやテストケースにおいても、何をテストしているか、その意図はなにか日本語コメントで記述する。
 - 数学的な処理に関するところは、SciOxideに対応する実装がないか見て、ある場合積極的に使う。SciOxideへの開発やIssue起票などはしないが、必要そうであればユーザーに提案すること。
 - Parse, don't validate を徹底する。外部入力や実行時に解決する情報は適切な境界で解析・検証し、不変条件が成立する型付きの値に変換する。内部 API は検証済みの型を受け取り、保証済みの条件を繰り返し検証しない。不変条件を破る値の構築・変更を防ぎ、未検証と検証済みのように中身の違う状態を型で区別する。実行時にしか判定できない制約や数値計算上の失敗は `Result` で扱う。
@@ -37,11 +37,31 @@ rhodium は、金融商品を小さなコア IR で表し、価格・quote 感
 - Rustの新機能 (1.100 の機能など) は有用な場合は積極的に使うこと。
 - 様々なAPI や利用例は、使用するバージョンに合わせて Context7 で調べる。リポジトリ内の型・関数の定義や参照、trait の実装をたどるときは Serena を活用する。
 
+## 設計の自己点検
+
+PR を出す前に、差分が次の項目に当たるかを確かめる。各項目は定義元の規則の索引で、判断は定義元に従う。項目の文言はここだけで管理し、規則の中身を変えるときは定義元の節を直す。
+
+- S1 同じ概念の型・関数・キャッシュを Serena で探し、横に別の入口を足さずにそれを使うか直す（[P7](docs/design/02-principles.md#p7-仕組みは理由があるときだけ足す)）。
+- S2 式、単位、分類、閾値を二か所目に書かず、定義元の関数か定数を呼ぶ（[変更を一か所に収める規則](docs/design/02-principles.md#変更を一か所に収める規則)）。
+- S3 検証済みの型を受け取った後で、同じ条件の検査や参照の引き直しをしない（[型の規約](docs/design/03-architecture.md#型の規約)）。
+- S4 `Error::Internal` を書く前に状態を型で表せないか確かめ、残すなら直前に `// 不変条件:` で保証の出どころを書く（[エラー](docs/design/04-foundation/04-reference-data.md#エラー)）。
+- S5 自作の enum を `_` や束縛で受けず、複数の variant にわたる分類は `matches!` や let-else でなく enum のメソッドの網羅した `match` に置く（[型の規約](docs/design/03-architecture.md#型の規約)）。
+- S6 ID の文字列を合成・分解・比較して意味を読まない（[型の規約](docs/design/03-architecture.md#型の規約)）。
+- S7 公開の関数や構築子を足すときは本番の呼び出し元を示し、テストも本番と同じ入口を通す（[P7](docs/design/02-principles.md#p7-仕組みは理由があるときだけ足す)）。
+- S8 取引・観測・残差の評価のループで、文字列の生成、`Arc` の複製、θ 全体の構築、カレンダーや地域の解決、準備や結合をしない（[性能の規則](docs/design/03-architecture.md#性能の規則)、[並列](docs/design/10-io-and-execution.md#並列)）。
+- S9 許容誤差と数値の閾値をリテラルで書かず、定義元のモジュールの定数を使う（[P2](docs/design/02-principles.md#p2-意味を一か所で定める)、[許容誤差](docs/design/11-quality.md#許容誤差)）。
+- S10 lint の例外は reason 付きの `#[expect]` を最小の範囲に付け、`clippy.toml` や検査の許可を足すときは PR に理由を書く（[lint](docs/design/03-architecture.md#lint)）。
+- S11 他の crate やディレクトリのファイルを `#[path]` で取り込まず、テストの共有の補助は[テストの補助](docs/design/11-quality.md#テストの補助)が定める場所に置く。
+- S12 テストはエラーを種類まで、集合と表を全体で比べ、fixture の書き換えを `if let` で黙って飛ばさない（[テストの種類](docs/design/11-quality.md#テストの種類)）。
+- S13 数秒を超えるテストは再校正を共有して短くするか、`#[ignore = "理由"]` を付ける（[CI](docs/design/11-quality.md#ci)）。
+- S14 型、モジュール、公開 API、ファイル形式、CLI を変えたら、地図と定義元の章を同じ PR で直す（[主要な型の地図](docs/design/03-architecture.md#主要な型の地図)）。
+- S15 `.githooks` の禁止パターンや lint に当たったら、検査を避ける書き換えでなく規則の意図どおりに直す（[CI](docs/design/11-quality.md#ci)）。
+
 ## 検証
 
 - 変更した振る舞いは、担当する受け入れ基準を確かめるテストで検証する。テストの名前は説明的な snake_case とし、オラクルと許容誤差は[品質](docs/design/11-quality.md)に従う。
 - 数値の比較は[許容誤差](docs/design/11-quality.md#許容誤差)で行い、ビット単位の一致は求めない。並列の計算は、スレッド数や分割方法を変えても許容誤差の範囲で一致することを確かめる。
-- 型で禁止する演算とブランドは、compile-fail テスト（trybuild）と正しい使い方のテストの対で確かめる。
+- 型で禁止する演算とブランドは、compile-fail テスト（testkit の `compile_fail_pair!`）と正しい使い方のテストの対で確かめる。
 - コミットの前に `cargo fmt --all --check`、`cargo clippy --locked --workspace --all-targets --all-features -- -D warnings`、`cargo test --locked --workspace --all-features`、`RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --all-features --no-deps` を通す。
 
 ## ブランチとプルリクエスト
@@ -50,4 +70,5 @@ rhodium は、金融商品を小さなコア IR で表し、価格・quote 感
 - 作業と検証が済んだら、作業ブランチから `main` にプルリクエストを出す。本文には変更の目的と検証結果を記し、実行できなかった検証があれば理由を添える。
 - コミットメッセージは `feat: ...` などのような形式の英語、PRのタイトルは英語、PRのbodyは日本語とする。
   - PR の body は、docs やコードへの参照に説明を委ねず、それ自体で変更の背景・目的・仕組み・影響を理解できるよう、自然な日本語で具体的かつ十分に詳しく記述すること。特に、どの変更箇所がどの機能追加・仕様変更に対応し、なぜ・どのように実装したのか、各テストが何を検証し何を保証するのかを明確に説明すること。可能な限りコミットごとに整理し、理解に必要なら PR の差分を超えた設計上の背景も補足すること。複雑な内容には数式・Mermaid 図・表・コード片・具体例などを適宜用い、単なる変更点やテスト名の列挙、過度な要約で済ませないこと。
+- 実装の PR の本文に「設計の自己点検」の節を置き、当たった項目ごとに結果（守った／例外）と根拠（file:line、テスト名、例外なら理由と追跡 issue）を表に書き、当たらない項目は「該当なし: S6、S11」のように番号だけを並べる。
 - 意味の単位で頻繁にコミットする。コミットは適切なタイミングで自由に行ってよいので、**各コミットを小さく保つ**。 必要があればPRの前にコミットを分割して意味の単位に整理するなどする。
diff --git a/Cargo.lock b/Cargo.lock
index a1e347284fe55c6d81e882aa28d25f70678cc5d8..08a7d72d3a48356c6b9a8acabec189164f3fa956 100644
--- a/Cargo.lock
+++ b/Cargo.lock
@@ -324,12 +324,6 @@ dependencies = [
  "rand_core",
 ]
 
-[[package]]
-name = "glob"
-version = "0.3.4"
-source = "registry+https://github.com/rust-lang/crates.io-index"
-checksum = "e4eba85ea1d0a966a983acd07deee566e67395d2d96b6fb39e62b5a833f1eb0b"
-
 [[package]]
 name = "half"
 version = "2.7.1"
@@ -600,8 +594,10 @@ dependencies = [
  "scioxide-ad",
  "scioxide-core",
  "scioxide-linalg",
+ "scioxide-optimize",
  "scioxide-rootfind",
  "scioxide-sparse",
+ "scioxide-special",
  "serde",
  "serde_json",
  "toml",
@@ -639,7 +635,9 @@ name = "rhodium-testkit"
 version = "0.1.0"
 dependencies = [
  "scioxide-core",
- "trybuild",
+ "serde",
+ "serde_json",
+ "tempfile",
 ]
 
 [[package]]
@@ -730,6 +728,22 @@ dependencies = [
  "scioxide-sparse",
 ]
 
+[[package]]
+name = "scioxide-optimize"
+version = "0.1.0"
+source = "git+https://github.com/kakune/scioxide?rev=a023b9d0bac667a6179ddbff2bf2d552f9d0a670#a023b9d0bac667a6179ddbff2bf2d552f9d0a670"
+dependencies = [
+ "scioxide-core",
+ "scioxide-differentiate",
+ "scioxide-interp",
+ "scioxide-linalg",
+ "scioxide-random",
+ "scioxide-sparse",
+ "scioxide-spatial",
+ "scioxide-special",
+ "scioxide-stats",
+]
+
 [[package]]
 name = "scioxide-poly"
 version = "0.1.0"
@@ -741,6 +755,18 @@ dependencies = [
  "scioxide-special",
 ]
 
+[[package]]
+name = "scioxide-random"
+version = "0.1.0"
+source = "git+https://github.com/kakune/scioxide?rev=a023b9d0bac667a6179ddbff2bf2d552f9d0a670#a023b9d0bac667a6179ddbff2bf2d552f9d0a670"
+dependencies = [
+ "scioxide-core",
+ "scioxide-differentiate",
+ "scioxide-linalg",
+ "scioxide-rootfind",
+ "scioxide-special",
+]
+
 [[package]]
 name = "scioxide-rootfind"
 version = "0.1.0"
@@ -753,6 +779,20 @@ dependencies = [
  "scioxide-special",
 ]
 
+[[package]]
+name = "scioxide-signal"
+version = "0.1.0"
+source = "git+https://github.com/kakune/scioxide?rev=a023b9d0bac667a6179ddbff2bf2d552f9d0a670#a023b9d0bac667a6179ddbff2bf2d552f9d0a670"
+dependencies = [
+ "scioxide-core",
+ "scioxide-fft",
+ "scioxide-interp",
+ "scioxide-linalg",
+ "scioxide-poly",
+ "scioxide-rootfind",
+ "scioxide-special",
+]
+
 [[package]]
 name = "scioxide-sparse"
 version = "0.1.0"
@@ -779,6 +819,19 @@ dependencies = [
  "scioxide-linalg",
 ]
 
+[[package]]
+name = "scioxide-stats"
+version = "0.1.0"
+source = "git+https://github.com/kakune/scioxide?rev=a023b9d0bac667a6179ddbff2bf2d552f9d0a670#a023b9d0bac667a6179ddbff2bf2d552f9d0a670"
+dependencies = [
+ "scioxide-core",
+ "scioxide-linalg",
+ "scioxide-random",
+ "scioxide-rootfind",
+ "scioxide-signal",
+ "scioxide-special",
+]
+
 [[package]]
 name = "serde"
 version = "1.0.229"
@@ -865,12 +918,6 @@ dependencies = [
  "unicode-ident",
 ]
 
-[[package]]
-name = "target-tuple"
-version = "1.0.2"
-source = "registry+https://github.com/rust-lang/crates.io-index"
-checksum = "876fef147edbcbddc8ac5cbbba92c7b86519e314e86638596c09673b2ed01e7f"
-
 [[package]]
 name = "tempfile"
 version = "3.27.0"
@@ -884,15 +931,6 @@ dependencies = [
  "windows-sys",
 ]
 
-[[package]]
-name = "termcolor"
-version = "1.4.1"
-source = "registry+https://github.com/rust-lang/crates.io-index"
-checksum = "06794f8f6c5c898b3275aebefa6b8a1cb24cd2c6c79397ab15774837a0bc5755"
-dependencies = [
- "winapi-util",
-]
-
 [[package]]
 name = "thiserror"
 version = "2.0.21"
@@ -962,21 +1000,6 @@ version = "1.2.0+spec-1.1.0"
 source = "registry+https://github.com/rust-lang/crates.io-index"
 checksum = "009dbb7700f11c55ef83aebbe4c0e9cf1f1ac0844bd8d6f00dc38f302026e13a"
 
-[[package]]
-name = "trybuild"
-version = "1.0.125"
-source = "registry+https://github.com/rust-lang/crates.io-index"
-checksum = "9838b494ee14d28f1eeea490e8f4f10248fa6dd371d4bc874b635de499408245"
-dependencies = [
- "glob",
- "serde",
- "serde_derive",
- "serde_json",
- "target-tuple",
- "termcolor",
- "toml",
-]
-
 [[package]]
 name = "unarray"
 version = "0.1.4"
diff --git a/Cargo.toml b/Cargo.toml
index d337b2674359aaee478fb182375008c70b0946d3..0592708867cf301ffdf2bd50f926411685af585c 100644
--- a/Cargo.toml
+++ b/Cargo.toml
@@ -28,8 +28,10 @@ scioxide-ad = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac66
 scioxide-core = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false, features = ["std"] }
 scioxide-interp = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
 scioxide-linalg = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
+scioxide-optimize = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
 scioxide-rootfind = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
 scioxide-sparse = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
+scioxide-special = { git = "https://github.com/kakune/scioxide", rev = "a023b9d0bac667a6179ddbff2bf2d552f9d0a670", default-features = false }
 rhodium = { path = "crates/rhodium" }
 rhodium-core = { path = "crates/rhodium-core" }
 rhodium-market = { path = "crates/rhodium-market" }
@@ -40,7 +42,6 @@ rhodium-xva = { path = "crates/rhodium-xva" }
 rhodium-testkit = { path = "crates/rhodium-testkit" }
 serde = { version = "1", features = ["derive"] }
 thiserror = "2"
-trybuild = "1.0"
 
 [workspace.lints.rust]
 rust_2018_idioms = { level = "deny", priority = -1 }
@@ -67,6 +68,7 @@ todo = "deny"
 unimplemented = "deny"
 unreachable = "deny"
 unwrap_used = "deny"
+wildcard_enum_match_arm = "deny"
 
 [workspace.lints.rustdoc]
 broken_intra_doc_links = "deny"
diff --git a/README.md b/README.md
index a0744f2ca67f8eedf16a24ed575c9966c4b32de7..9955a57930fdc3531ace827970f3d15ace7679d5 100644
--- a/README.md
+++ b/README.md
@@ -28,12 +28,26 @@ quote への感応度は、取引の θ 勾配を集計してから校正の随
 
 ## 現状
 
-マイルストーン M0（JPY と USD の OIS 曲線とスワップの PV と quote デルタのプロトタイプ）の受け入れ基準 M0-1〜M0-9 を満たしました（[ロードマップ](docs/design/12-roadmap.md)、Issue #16）。
-ワークスペースは上の 8 crate と開発用の `rhodium-testkit` からなり、旧設計の crate は残っていません。
-1 万件のスワップの PV と quote デルタの所要時間は [M0 の計測](docs/benchmarks/m0.md) にあります。
-次のマイルストーンは M1（金利リニアの完成）です。
+マイルストーン M1（金利リニアの完成）の受け入れ基準 M1-1〜M1-8 を満たしました（[ロードマップ](docs/design/12-roadmap.md)、Issue #37）。
+JPY（TONA の会合日のステップと年末のジャンプ、TIBOR 6M・3M のベーシス）、USD（SOFR）、USDJPY（FX の短期と MtM XCCY の長期による外貨担保の曲線）を一つの市場の定義で校正し、
+FRA、金利先物、IRS、ベーシススワップ、OIS（lookback、observation shift、lockout、支払いの遅れ）、FX フォワードと FX スワップ、XCCY の PV と quote デルタを、CSA の 4 つの場合の割引で求めます。
+所要時間は [M0 の計測](docs/benchmarks/m0.md) と [M1 の計測](docs/benchmarks/m1.md) にあります。
+次のマイルストーンは M2（ボラと解析オプション）です。
 
-## M0 のサンプル
+ワークスペースは上の 8 crate と開発用の `rhodium-testkit` からなり、旧設計の crate は残っていません。
+
+## サンプル
+
+### M1 の金利リニア
+
+最初の検証データセットの金利商品を評価するワークスペースです。[市場の定義と取引の説明](examples/m1-rates/README.md)
+
+```sh
+cargo run --locked -p rhodium-cli -- price examples/m1-rates --date 2026-11-24
+cargo run --locked -p rhodium-cli -- risk examples/m1-rates --date 2026-11-24 --threads 4
+```
+
+### M0 の OIS
 
 JPY TONA と USD SOFR の OIS を評価するサンプルワークスペースを同梱しています。
 
@@ -72,11 +86,12 @@ RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --all-features --no-de
 ### Git フック
 
 Git フックを有効にするには、クローンごとに `git config core.hooksPath .githooks` を実行します。
-コミット時はステージした差分の空白と、変更ファイルの Rust 書式、Python・TOML・JSON・Nix・フックの構文を検査します。
-push 時は送信先との差分に同じ検査を適用し、変更した crate の Clippy、テスト、rustdoc を実行します。
-作業ツリーに送らない変更があると、送る commit を一時ディレクトリに取り出し、ビルド先もその中に分けて検査します（依存から作り直すので時間がかかります）。
+コミット時はステージした差分の空白と、変更ファイルの Rust 書式、Python・TOML・JSON・Nix・フックの構文、追加した行の禁止パターン（[品質の CI](docs/design/11-quality.md#ci)）を検査します。
+push 時は送信先との差分に同じ検査を適用し、変更した crate とそれに依存する crate の Clippy、テスト、rustdoc を実行します。
+`python3 .githooks/checks.py audit` は、作業ツリーに残る禁止パターンの違反を規則ごとに数えます。
+pre-push の対象パスは送信する ref の差分で決めます。検査には未コミットの変更を含む現在の作業ツリーと、通常の Cargo ビルド先を使います。別ブランチを push するときも同じです。
 ワークスペース共通の Cargo 設定を変えた場合は全 crate が対象です。
-`.githooks/` の Python の変更には、その `pyproject.toml` に設定した Ruff と ty を適用します。
+`.githooks/` の Python の変更には、その `pyproject.toml` に設定した Ruff と ty と、`.githooks/test_*.py` の単体テストを適用します。
 新規ブランチでは送信先の追跡 ref と共有する履歴を基準にし、共有履歴がなければ初回の push としてツリー全体を検査します。
 
 ## ライセンス
diff --git a/crates/rhodium-cli/src/main.rs b/crates/rhodium-cli/src/main.rs
index 5c33ce2f157444ad2b668e5a4b212f98dbd18397..5718361e554009bb59b5f952b418cf718abfbd23 100644
--- a/crates/rhodium-cli/src/main.rs
+++ b/crates/rhodium-cli/src/main.rs
@@ -1,4 +1,4 @@
-//! ワークスペースの PV と quote デルタを CSV で出力する。
+//! ワークスペースの校正と評価の結果を CSV で出力する。
 
 use std::{
     error::Error,
@@ -10,15 +10,15 @@ use std::{
 
 use clap::{Args, Parser, Subcommand};
 use rhodium::{
-    prelude::{ExecutionContext, ThreadExecutor},
+    prelude::{ExecutionContext, ThreadExecutor, VolKind},
     session::{GroupOutcome, Session},
-    workspace::Workspace,
+    workspace::{CalibratedMarket, Workspace},
 };
 
 #[derive(Parser)]
 #[command(
     version,
-    about = "Price trades and compute quote delta from a rhodium workspace"
+    about = "Calibrate a rhodium market, price trades, or compute quote delta"
 )]
 struct Cli {
     #[command(subcommand)]
@@ -31,6 +31,8 @@ enum Command {
     Price(ValuationArgs),
     /// Write portfolio delta by quote and currency, in each quote's reporting unit.
     Risk(ValuationArgs),
+    /// Write the repricing residual of each calibrated volatility quote.
+    Calibrate(ValuationArgs),
 }
 
 #[derive(Args)]
@@ -82,13 +84,17 @@ fn main() -> ExitCode {
 
 // 入力を読み、評価時点を解決し、校正と取引の診断を書いてから CSV を出す。
 // 一部の校正や取引が失敗しても計算できた行は出し、戻り値の終了コードを 1 にする。
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#198 で各 variant を網羅する"
+)]
 fn run(
     command: &Command,
     output: impl Write,
     diagnostics: &mut impl Write,
 ) -> Result<ExitCode, CommandError> {
     let args = match command {
-        Command::Price(args) | Command::Risk(args) => args,
+        Command::Price(args) | Command::Risk(args) | Command::Calibrate(args) => args,
     };
     // 取引を分けるワーカーの数。既定の 1 は逐次で、スレッドを作らない。
     let executor = ThreadExecutor::new(args.threads);
@@ -106,22 +112,63 @@ fn run(
         _ => return Err(input_error("specify one of --as-of and --date")),
     }
     .map_err(input_error)?;
-    // 校正グループ単位の失敗は GroupOutcome に残るため、Session::new 全体の失敗は
-    // 評価時点が曲線のピラー以降にあるなど、入力どうしの不整合で市場を組み立てられない場合に当たる。
-    // 計算できた行が一つも無いので、部分的な失敗（1）ではなく入力のエラー（2）として返す。
-    // 並列の実行の失敗（Execution など）は入力の誤りではないので 1 にする。
-    let session = Session::new(workspace, as_of, execution).map_err(|error| match error {
+    let failed = match command {
+        Command::Calibrate(_) => {
+            let calibrated = workspace.calibrate(as_of).map_err(calibration_error)?;
+            let failed = calibration_diagnostics(&calibrated, diagnostics)?;
+            write_volatility_fit(&calibrated, output)?;
+            failed
+        }
+        Command::Price(_) | Command::Risk(_) => {
+            let session = Session::new(workspace, as_of, execution).map_err(calibration_error)?;
+            let mut failed = calibration_diagnostics(session.calibrated_market(), diagnostics)?;
+            for (id, warnings) in session.preparation_warnings() {
+                for warning in warnings {
+                    writeln!(diagnostics, "warning: trade {id}: {warning:?}")
+                        .map_err(execution_error)?;
+                }
+            }
+            let mut writer = csv::WriterBuilder::new()
+                .has_headers(false)
+                .from_writer(output);
+            failed |= if matches!(command, Command::Price(_)) {
+                price(&session, execution, &mut writer, diagnostics)?
+            } else {
+                risk(&session, execution, &mut writer, diagnostics)?
+            };
+            writer.flush().map_err(execution_error)?;
+            failed
+        }
+    };
+    Ok(if failed {
+        ExitCode::FAILURE
+    } else {
+        ExitCode::SUCCESS
+    })
+}
+
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#198 で各 variant を網羅する"
+)]
+fn calibration_error(error: rhodium::prelude::Error) -> CommandError {
+    match error {
         rhodium::prelude::Error::InvalidInput(_)
         | rhodium::prelude::Error::InvalidData(_)
         | rhodium::prelude::Error::MissingData(_) => input_error(error),
         _ => execution_error(error),
-    })?;
+    }
+}
 
+fn calibration_diagnostics(
+    calibrated: &CalibratedMarket,
+    diagnostics: &mut impl Write,
+) -> Result<bool, CommandError> {
     let mut failed = false;
-    for warning in session.calibrated_market().warnings() {
+    for warning in calibrated.warnings() {
         writeln!(diagnostics, "warning: calibration: {warning:?}").map_err(execution_error)?;
     }
-    for group in session.calibrated_market().groups() {
+    for group in calibrated.groups() {
         let reason = match group.outcome() {
             GroupOutcome::Solved(_) => continue,
             GroupOutcome::Failed(error) => error.to_string(),
@@ -137,25 +184,46 @@ fn run(
             .map_err(execution_error)?;
         failed = true;
     }
-    for (id, warnings) in session.preparation_warnings() {
-        for warning in warnings {
-            writeln!(diagnostics, "warning: trade {id}: {warning:?}").map_err(execution_error)?;
-        }
-    }
+    Ok(failed)
+}
 
+fn write_volatility_fit(
+    calibrated: &CalibratedMarket,
+    output: impl Write,
+) -> Result<(), CommandError> {
     let mut writer = csv::WriterBuilder::new()
         .has_headers(false)
         .from_writer(output);
-    failed |= match command {
-        Command::Price(_) => price(&session, execution, &mut writer, diagnostics)?,
-        Command::Risk(_) => risk(&session, execution, &mut writer, diagnostics)?,
-    };
-    writer.flush().map_err(execution_error)?;
-    Ok(if failed {
-        ExitCode::FAILURE
-    } else {
-        ExitCode::SUCCESS
-    })
+    writer
+        .write_record([
+            "surface",
+            "quote_id",
+            "vol_kind",
+            "target_vol",
+            "model_vol",
+            "residual",
+        ])
+        .map_err(execution_error)?;
+    let mut fits: Vec<_> = calibrated.volatility_fit().collect();
+    fits.sort_by_key(|fit| (fit.surface, fit.quote));
+    for fit in fits {
+        let kind = match fit.kind {
+            VolKind::Normal => "normal",
+            VolKind::Lognormal => "lognormal",
+            VolKind::ShiftedLognormal { .. } => "shifted_lognormal",
+        };
+        writer
+            .serialize((
+                fit.surface.as_str(),
+                fit.quote.as_str(),
+                kind,
+                fit.target,
+                fit.model,
+                fit.residual,
+            ))
+            .map_err(execution_error)?;
+    }
+    writer.flush().map_err(execution_error)
 }
 
 // 取引 ID 順、通貨順に PV を書く。失敗した取引は診断に回し、失敗の有無を返す。
diff --git a/crates/rhodium-cli/tests/workspace_cli.rs b/crates/rhodium-cli/tests/workspace_cli.rs
index 5b4db5950b882c1704f356865fe85a34d87a723c..ebe28a7fc7807a03fdcb6efec5bcc6a87e145fd3 100644
--- a/crates/rhodium-cli/tests/workspace_cli.rs
+++ b/crates/rhodium-cli/tests/workspace_cli.rs
@@ -1011,3 +1011,178 @@ day_count = "ACT360"
     assert_close!(price.2, 0.0_f64, 1e-12_f64, 1e-8_f64);
     Ok(())
 }
+
+#[expect(
+    clippy::disallowed_macros,
+    reason = "リポジトリのサンプルの場所をテストのソースの位置から決める"
+)]
+fn m1_sample() -> PathBuf {
+    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m1-rates")
+}
+
+#[test]
+fn cli_prices_and_risks_the_m1_rates_workspace() -> TestResult {
+    // 最初の検証データセットのサンプルを rh price と rh risk で評価する。全 18 取引に PV の行があり、
+    // 近似の警告は予想の測度に従う（TONA 担保の TIBOR の IRS は CrossCollateralProjection を出さず、
+    // 25bp の CSA の IRS は出す）。ワーカー数を変えても CSV の値は許容誤差で一致する。
+    let dir = m1_sample();
+    let price = run("price", &dir, &["--date", "2026-11-24"])?;
+    assert!(
+        price.status.success(),
+        "{}",
+        String::from_utf8_lossy(&price.stderr)
+    );
+    let pvs = rows(&price, &["trade_id", "currency", "pv"])?;
+    let trades: std::collections::BTreeSet<_> = pvs.iter().map(|(id, _, _)| id.as_str()).collect();
+    assert_eq!(trades.len(), 18);
+    let stderr = String::from_utf8(price.stderr)?;
+    assert!(!stderr.contains("error"), "{stderr}");
+    assert!(stderr.contains("trade JPY-IRS-5Y-CASH-25: CrossCollateralProjection"));
+    assert!(!stderr.contains("trade JPY-IRS-5Y: CrossCollateralProjection"));
+    assert!(stderr.contains("trade USD-OIS-JPY-CASH: CrossCollateralProjection"));
+
+    let risk = run("risk", &dir, &["--date", "2026-11-24"])?;
+    assert!(risk.status.success());
+    let deltas = rows(&risk, &["quote_id", "currency", "sensitivity"])?;
+    for quote in [
+        "JPY.TONA.OIS.5Y",
+        "JPY.TONA.TURN.2026",
+        "JPY.TIBOR6M-TONA.BASIS.5Y",
+        "JPY.TIBOR3M-TONA.BASIS.5Y",
+        "USD.SOFR.OIS.2Y",
+        "USDJPY.SPOT",
+        "USDJPY.XCCY.2Y",
+        "JPY.LGM.VOL",
+    ] {
+        assert!(
+            deltas
+                .iter()
+                .any(|(id, _, value)| id == quote && value.abs() > 1e-6),
+            "{quote}"
+        );
+    }
+
+    for (command, values, headers) in [
+        ("price", &pvs, ["trade_id", "currency", "pv"]),
+        ("risk", &deltas, ["quote_id", "currency", "sensitivity"]),
+    ] {
+        let parallel = run(command, &dir, &["--date", "2026-11-24", "--threads", "2"])?;
+        assert!(parallel.status.success());
+        let parallel = rows(&parallel, &headers)?;
+        assert_eq!(parallel.len(), values.len());
+        for ((id, ccy, value), (parallel_id, parallel_ccy, parallel_value)) in
+            values.iter().zip(&parallel)
+        {
+            assert_eq!((id, ccy), (parallel_id, parallel_ccy));
+            assert_close!(*value, *parallel_value, 1e-10_f64, 1e-6_f64);
+        }
+    }
+    Ok(())
+}
+
+#[test]
+fn calibrate_reads_rate_volatility_and_reports_partial_failures() -> TestResult {
+    // 同じワークスペースで、SABR の ATM とターム・RFR のキャップを校正する。
+    let dir = tempfile::tempdir()?;
+    copy_dir(&sample().with_file_name("m1-rates"), dir.path())?;
+    let path = dir.path().join("market/definition.toml");
+    let definitions = r#"
+
+[[volatilities]]
+id = "JPY-SABR"
+type = "swaption_cube"
+index = "JPY-TONA"
+convention = "JPY-TONA-OIS"
+collateral = "uncollateralized"
+vol_kind = "normal"
+expiry_time = { time = "11:00", zone = "Asia/Tokyo" }
+sabr = { beta = 0.5, shift = 0.01, fit = { atm_only = { rho = -0.15, nu = 0.4 } } }
+quotes = [{ quote = "JPY.SABR.1Y.2Y.ATM", expiry = "1Y", tenor = "2Y", strike = "atm" }]
+
+[[volatilities]]
+id = "JPY-TIBOR-CAP"
+type = "caplet"
+index = "JPY-TIBOR-6M"
+convention = "JPY-TIBOR6M-IRS"
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+vol_kind = "normal"
+quotes = [
+    { quote = "JPY.CAP.1Y.ATM", maturity = "1Y", strike = "atm" },
+    { quote = "JPY.CAP.2Y.ATM", maturity = "2Y", strike = "atm" },
+    { quote = "JPY.CAP.1Y.K1", maturity = "1Y", strike = { absolute = 0.01 } },
+    { quote = "JPY.CAP.2Y.K1", maturity = "2Y", strike = { absolute = 0.01 } },
+]
+
+[[volatilities]]
+id = "JPY-TONA-CAP"
+type = "caplet"
+index = "JPY-TONA"
+convention = "JPY-TONA-OIS"
+collateral = "uncollateralized"
+vol_kind = "normal"
+rfr_variance = "lyashenko_mercurio"
+quotes = [
+    { quote = "JPY.RFRCAP.1Y.ATM", maturity = "1Y", strike = "atm" },
+    { quote = "JPY.RFRCAP.2Y.ATM", maturity = "2Y", strike = "atm" },
+    { quote = "JPY.RFRCAP.1Y.K1", maturity = "1Y", strike = { absolute = 0.01 } },
+    { quote = "JPY.RFRCAP.2Y.K1", maturity = "2Y", strike = { absolute = 0.01 } },
+]
+"#;
+    let original = fs::read_to_string(&path)?;
+    fs::write(path, original + definitions)?;
+    let path = dir.path().join("quotes/snapshot.csv");
+    let mut quotes = fs::read_to_string(&path)?;
+    let good_quotes = [
+        "JPY.SABR.1Y.2Y.ATM",
+        "JPY.CAP.1Y.ATM",
+        "JPY.CAP.2Y.ATM",
+        "JPY.CAP.1Y.K1",
+        "JPY.CAP.2Y.K1",
+        "JPY.RFRCAP.1Y.ATM",
+        "JPY.RFRCAP.2Y.ATM",
+        "JPY.RFRCAP.1Y.K1",
+        "JPY.RFRCAP.2Y.K1",
+    ];
+    for quote in good_quotes {
+        writeln!(quotes, "{quote},60bp")?;
+    }
+    fs::write(&path, &quotes)?;
+    let output = run("calibrate", dir.path(), &["--date", "2026-11-24"])?;
+    assert!(
+        output.status.success(),
+        "{}",
+        String::from_utf8_lossy(&output.stderr)
+    );
+    let mut csv = csv::Reader::from_reader(output.stdout.as_slice());
+    assert!(csv.headers()?.iter().eq([
+        "surface",
+        "quote_id",
+        "vol_kind",
+        "target_vol",
+        "model_vol",
+        "residual"
+    ]));
+    let records = csv.records().collect::<Result<Vec<_>, _>>()?;
+    assert_eq!(records.len(), good_quotes.len());
+    assert!(records.is_sorted_by(|a, b| (&a[0], &a[1]) < (&b[0], &b[1])));
+    for row in &records {
+        assert_eq!(&row[2], "normal");
+        assert_close!(row[3].parse::<f64>()?, 0.006, 1e-12, 1e-14);
+        assert_close!(row[4].parse::<f64>()?, 0.006, 1e-12, 1e-14);
+        assert!(row[5].parse::<f64>()?.abs() <= 1e-14);
+    }
+
+    // タームキャップの quote を欠損させても、SABR と RFR キャップは出力される。
+    fs::write(path, quotes.replace("JPY.CAP.1Y.ATM,60bp\n", ""))?;
+    let output = run(
+        "calibrate",
+        dir.path(),
+        &["--as-of", "2026-11-24T08:00:00Z"],
+    )?;
+    assert_eq!(output.status.code(), Some(1));
+    let diagnostics = String::from_utf8(output.stderr)?;
+    assert!(diagnostics.contains("JPY-TIBOR-CAP") && diagnostics.contains("JPY.CAP.1Y.ATM"));
+    let mut csv = csv::Reader::from_reader(output.stdout.as_slice());
+    assert_eq!(csv.records().count(), 5);
+    Ok(())
+}
diff --git a/crates/rhodium-core/src/calendar.rs b/crates/rhodium-core/src/calendar.rs
index d66ea0e6afe75ce3015a53bf6e72bc1b718285f9..e5915e121fa0c0b0d33ae19748c50acea1707954 100644
--- a/crates/rhodium-core/src/calendar.rs
+++ b/crates/rhodium-core/src/calendar.rs
@@ -170,18 +170,19 @@ impl Calendar {
         // 基本暦の既知範囲が交わらなければ範囲は空になり、ビット集合も空になる。
         let mut known = Date::from_ymd(1, 1, 1)?..=Date::from_ymd(9999, 12, 31)?;
         let expression = ResolvedSpec::resolve(spec, bases, &mut known)?;
-        let length = (*known.end()).days_since(*known.start()) + 1;
-        let length = index(length.max(0))?;
-        let mut business = vec![0_u64; length.div_ceil(64)];
+        // 既知範囲の先頭からの日数を i32 で数え、ビットの位置は列挙の番号から取る。範囲が空なら
+        // 日数の範囲も空になる。語は 64 日ごとに足すので、立てるビットの語は常にある。
+        let mut business = Vec::new();
         let mut warnings = Vec::new();
-        for day_index in 0..length {
-            let days = i32::try_from(day_index).map_err(|_| calendar_index_error())?;
+        for (bit, days) in (0..=known.end().days_since(*known.start())).enumerate() {
             let date = known.start().checked_add_days(days)?;
-            if expression.is_business_day(date, &mut warnings)? {
-                let word = business
-                    .get_mut(day_index / 64)
-                    .ok_or_else(calendar_index_error)?;
-                *word |= 1 << (day_index % 64);
+            if bit % 64 == 0 {
+                business.push(0_u64);
+            }
+            if expression.is_business_day(date, &mut warnings)?
+                && let Some(word) = business.last_mut()
+            {
+                *word |= 1 << (bit % 64);
             }
         }
         Ok(Self(Arc::new(Compiled {
@@ -216,23 +217,32 @@ impl Calendar {
         let length = period.end().days_since(period.start());
         let known_start = self.0.known.start().days_since(start).clamp(0, length);
         let known_end = (self.0.known.end().days_since(start) + 1).clamp(known_start, length);
-        let mut business = 0;
-        if known_start < known_end {
-            let first = start
-                .checked_add_days(known_start)?
-                .days_since(*self.0.known.start());
-            business = self
-                .0
-                .count_bits(index(first)?, index(first + known_end - known_start)?)?;
-        }
+        let counted = if known_start < known_end {
+            self.0.count_known(
+                start.checked_add_days(known_start)?,
+                start.checked_add_days(known_end - 1)?,
+            )
+        } else {
+            Some(0)
+        };
+        // ビット集合は既知範囲を覆うように作るので、数えられないことは無い。数えられなければ、
+        // 既知範囲の日も式で評価する（同じ値になる）。
+        let (mut business, evaluated_until) = match counted {
+            Some(count) => (count, known_start),
+            None => (0, known_end),
+        };
         let mut warnings = Vec::new();
-        for days in (0..known_start).chain(known_end..length) {
+        for days in (0..evaluated_until).chain(known_end..length) {
             let date = start.checked_add_days(days)?;
             business += u32::from(self.0.expression.is_business_day(date, &mut warnings)?);
         }
         Ok((business, unique_warnings(warnings)))
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#198 で各 variant を網羅する"
+    )]
     pub(crate) fn base_id(&self) -> Option<&CalendarId> {
         match &self.0.expression {
             ResolvedSpec::Base(id, _) => Some(id),
@@ -241,48 +251,57 @@ impl Calendar {
     }
 
     pub(crate) fn lookup(&self, date: Date, warnings: &mut Vec<Warning>) -> Result<bool> {
-        if self.0.known.contains(&date) {
-            let bit = index(date.days_since(*self.0.known.start()))?;
-            let word = self
-                .0
-                .business
-                .get(bit / 64)
-                .ok_or_else(calendar_index_error)?;
-            return Ok(word & (1 << (bit % 64)) != 0);
+        // 既知範囲の外の日は、外挿の警告や欠損を返すために式を評価する。
+        match self.0.known_bit(date) {
+            Some(business) => Ok(business),
+            None => self.0.expression.is_business_day(date, warnings),
         }
-        self.0.expression.is_business_day(date, warnings)
     }
 }
 
 impl Compiled {
-    // ビット添字の半開区間 [first, after) で立っているビットを数える。呼び出し側は first < after を
-    // 保証する。先頭の語は first_bit より下のビットを、末尾の語は last_bit より上のビットを落とす。
-    fn count_bits(&self, first: usize, after: usize) -> Result<u32> {
-        let first_word = first / 64;
-        let last_word = (after - 1) / 64;
-        let first_bit = first % 64;
-        let last_bit = (after - 1) % 64;
-        let words = self
-            .business
-            .get(first_word..=last_word)
-            .ok_or_else(calendar_index_error)?;
-        Ok(words
-            .iter()
-            .enumerate()
-            .map(|(index, &word)| {
-                let lower = if index == 0 {
-                    u64::MAX << first_bit
-                } else {
-                    u64::MAX
-                };
-                let upper = if index == last_word - first_word {
-                    u64::MAX >> (63 - last_bit)
-                } else {
-                    u64::MAX
-                };
-                (word & lower & upper).count_ones()
-            })
-            .sum())
+    // 既知範囲の先頭からのビットの位置。範囲の外の日は None。
+    fn offset(&self, date: Date) -> Option<usize> {
+        if !self.known.contains(&date) {
+            return None;
+        }
+        usize::try_from(date.days_since(*self.known.start())).ok()
+    }
+
+    // 既知範囲の日の営業日のビット。範囲の外なら None で、呼び出し側が式を評価する。
+    fn known_bit(&self, date: Date) -> Option<bool> {
+        let bit = self.offset(date)?;
+        self.business
+            .get(bit / 64)
+            .map(|word| word & (1 << (bit % 64)) != 0)
+    }
+
+    // 既知範囲の日 first から last まで（両端を含む）の営業日を数える。first か last が範囲の外なら
+    // None。先頭の語は first のビットより下を、末尾の語は last のビットより上を落とす。
+    fn count_known(&self, first: Date, last: Date) -> Option<u32> {
+        let first = self.offset(first)?;
+        let last = self.offset(last)?;
+        let words = self.business.get(first / 64..=last / 64)?;
+        let last_index = words.len() - 1;
+        Some(
+            words
+                .iter()
+                .enumerate()
+                .map(|(index, &word)| {
+                    let lower = if index == 0 {
+                        u64::MAX << (first % 64)
+                    } else {
+                        u64::MAX
+                    };
+                    let upper = if index == last_index {
+                        u64::MAX >> (63 - last % 64)
+                    } else {
+                        u64::MAX
+                    };
+                    (word & lower & upper).count_ones()
+                })
+                .sum(),
+        )
     }
 }
 
@@ -364,14 +383,6 @@ pub(crate) fn unique_warnings(mut warnings: Vec<Warning>) -> Vec<Warning> {
     warnings
 }
 
-fn index(days: i32) -> Result<usize> {
-    usize::try_from(days).map_err(|_| calendar_index_error())
-}
-
-fn calendar_index_error() -> Error {
-    Error::Internal("calendar bit index outside the compiled range".into())
-}
-
 fn validate_group_len(length: usize) -> Result<()> {
     if length < 2 {
         return Err(Error::InvalidData(
diff --git a/crates/rhodium-core/src/calendar/advance.rs b/crates/rhodium-core/src/calendar/advance.rs
index cf98ce611972e196ba8a32d159d5daaca757af42..b9abace3f314f86d8065ad92b7e5ca2af2bcfee6 100644
--- a/crates/rhodium-core/src/calendar/advance.rs
+++ b/crates/rhodium-core/src/calendar/advance.rs
@@ -132,6 +132,10 @@ pub fn money_market_period(
     Ok((DateInterval::new(start, end)?, unique_warnings(warnings)))
 }
 
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#198 で各 variant を網羅する"
+)]
 pub(crate) fn advance_date(
     date: Date,
     amount: i32,
diff --git a/crates/rhodium-core/src/error.rs b/crates/rhodium-core/src/error.rs
index 4575b102b14c80a5f37299bde852af565023ac82..89888803d5e74a022babf253b52980775635add4 100644
--- a/crates/rhodium-core/src/error.rs
+++ b/crates/rhodium-core/src/error.rs
@@ -77,6 +77,25 @@ pub enum Warning {
         /// 直接のレートを使うペア。
         pair: CcyPair,
     },
+    /// ボラ校正のストライク格子でコール価格の凸性が失われた。
+    VolatilityButterfly {
+        /// ボラ面の ID。
+        surface: String,
+    },
+    /// 最小二乗の校正で、残差が binding の許容幅を超えた（当てはめの許容幅の外）。
+    ///
+    /// 過剰決定の当てはめでは正しい解でも残差が残るので、解は使い、quote を報告する。
+    CalibrationFitOutsideTolerance {
+        /// 残差が許容幅を超えた quote。
+        quote: QuoteId,
+    },
+    /// 最小二乗のグループの解の記録を Gauss–Newton（残差の 2 階の項を省く）で作った。
+    ///
+    /// このグループを通る quote の感応度は近似になる。厳密な記録は 2 階の項を含む。
+    CalibrationGaussNewton {
+        /// グループが読む quote。
+        quote: QuoteId,
+    },
 }
 
 /// DSL のソース上でエラーが起きた位置。
diff --git a/crates/rhodium-core/src/observation.rs b/crates/rhodium-core/src/observation.rs
index 54a4ebe3300e063faa160297aab4146f1edd4721..6cbd58c62d73df135eae1501d670c67dd34a5f46 100644
--- a/crates/rhodium-core/src/observation.rs
+++ b/crates/rhodium-core/src/observation.rs
@@ -431,6 +431,10 @@ pub fn fallback_observation_schedule(
     )
 }
 
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#197 で各 variant を網羅する"
+)]
 fn observation_schedule_using(
     accrual: DateInterval,
     rule: &ObservationRule,
@@ -529,6 +533,10 @@ fn shifted_interval(
     DateInterval::new(start, end)
 }
 
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#198 で各 variant を網羅する"
+)]
 fn business_observations(
     interval: DateInterval,
     method: ObservationMethod,
diff --git a/crates/rhodium-core/src/parallel.rs b/crates/rhodium-core/src/parallel.rs
index 9de8fdf32cda34894d4eb664c062e992a53dbbb7..9b9b4f2cd0f44432d9b1b484ecf181f3e7bb3e54 100644
--- a/crates/rhodium-core/src/parallel.rs
+++ b/crates/rhodium-core/src/parallel.rs
@@ -70,6 +70,8 @@ pub fn map_partitions<T: Send>(
     results
         .into_iter()
         .map(|result| {
+            // 不変条件: execute_partitioned が Ok を返したなら、実行器はすべての区間を一度ずつ
+            // 実行して出力の場所を埋めている（scioxide の実行器の約束）。
             result.ok_or_else(|| Error::Internal("executor left a partition without result".into()))
         })
         .collect()
@@ -95,6 +97,10 @@ pub fn check(execution: ExecutionContext<'_>) -> Result<()> {
 /// 実行の制御の失敗を対象ごとの失敗に変えると、止めた計算が一部の対象を除いた結果として
 /// 成功してしまうため、この分類で分ける。
 #[must_use]
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "外部の enum は variant が増えうるので、既定の腕で受ける"
+)]
 pub fn control_failure(error: &SciOxideError) -> Option<Error> {
     match error.view() {
         // 呼び出し側が渡した取り消し・期限・時計と、実行環境の資源の不足。
@@ -115,6 +121,7 @@ pub fn control_failure(error: &SciOxideError) -> Option<Error> {
 #[expect(clippy::needless_pass_by_value, reason = "map_err の変換関数に使う")]
 fn execution_error(error: SciOxideError) -> Error {
     // 実行の制御でなければ、コールバックの panic か分割の約束を破る実行器で、rhodium の誤り。
+    // 不変条件: 区間の処理は失敗を Result で返し panic しない。実行器は分割の約束を守る。
     control_failure(&error)
         .unwrap_or_else(|| Error::Internal(format!("parallel execution: {error}")))
 }
diff --git a/crates/rhodium-core/src/quantity.rs b/crates/rhodium-core/src/quantity.rs
index 248aae1c7002b9f4d7821fe0fe2ca84460ceae9d..3f5e09d7b9898170315543a92d65913f0f4d87e4 100644
--- a/crates/rhodium-core/src/quantity.rs
+++ b/crates/rhodium-core/src/quantity.rs
@@ -16,12 +16,30 @@ pub use volatility::{CallPut, LognormalVol, NormalVol, ShiftedLognormalVol, VolK
 use crate::{Error, Result, day_count::DayCount};
 
 /// 日数基準のタグを持つ日数係数。
-#[derive(Debug, Clone, PartialEq)]
+///
+/// ファイルでは `{ value = 0.0164, basis = "ACT365F" }` と書き、読み込みで有限性を確かめる。
+#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
+#[serde(try_from = "YearFractionInput")]
 pub struct YearFraction {
     value: f64,
     basis: DayCount,
 }
 
+#[derive(serde::Deserialize)]
+#[serde(deny_unknown_fields)]
+struct YearFractionInput {
+    value: f64,
+    basis: DayCount,
+}
+
+impl TryFrom<YearFractionInput> for YearFraction {
+    type Error = Error;
+
+    fn try_from(input: YearFractionInput) -> Result<Self> {
+        Self::new(input.value, input.basis)
+    }
+}
+
 impl YearFraction {
     /// 有限な日数係数に日数基準を指定する。
     ///
diff --git a/crates/rhodium-core/src/quantity/volatility.rs b/crates/rhodium-core/src/quantity/volatility.rs
index 0594dfd277d7caf5b734a1507de311731d6edf74..63ccc6eb0ce5b9dc1588f3476126b203aa21c635 100644
--- a/crates/rhodium-core/src/quantity/volatility.rs
+++ b/crates/rhodium-core/src/quantity/volatility.rs
@@ -106,10 +106,9 @@ impl ShiftedLognormalVol {
     /// シフトが異なれば `InvalidInput`。
     pub fn compare(self, other: Self) -> Result<Ordering> {
         self.check_shift(other)?;
-        // -0.0 と 0.0 は同じボラとして比較する。
-        self.value
-            .partial_cmp(&other.value)
-            .ok_or_else(|| Error::Internal("validated volatility is not comparable".into()))
+        // 構築時にボラを有限の非負に限り、-0.0 を 0.0 にそろえるので、全順序の比較が
+        // 数としての比較と一致し、-0.0 と 0.0 も同じボラになる。
+        Ok(self.value.total_cmp(&other.value))
     }
 
     /// 同じシフトを持つボラとの差。負の差も返す。
@@ -135,5 +134,6 @@ fn nonnegative(value: f64, name: &str) -> Result<f64> {
     if value < 0.0 {
         return Err(Error::InvalidInput(format!("{name} is negative")));
     }
-    Ok(value)
+    // 非負の値の絶対値は自身で、-0.0 だけが 0.0 になる。ボラの比較と表示で符号付きの零を区別しない。
+    Ok(value.abs())
 }
diff --git a/crates/rhodium-core/src/refdata.rs b/crates/rhodium-core/src/refdata.rs
index 6e51d7789b6f1883226841003da1765f55c25e19..dc1088651160557ffaa30c257ff62e244d632084 100644
--- a/crates/rhodium-core/src/refdata.rs
+++ b/crates/rhodium-core/src/refdata.rs
@@ -1,6 +1,7 @@
 //! 通貨、地域、カレンダー、指標と慣行の参照データ。
 
 mod fallback;
+mod fx_convention;
 mod index;
 mod rate_future;
 
@@ -31,6 +32,7 @@ use crate::{
 pub use fallback::{
     FallbackApplicability, FallbackDefinition, FallbackReplacement, FallbackTrigger,
 };
+pub use fx_convention::{FxConvention, FxSpotLag};
 pub use index::{
     DailyPublication, IndexDefinition, IndexKind, IndexStatus, InflationInterpolation, PriceKind,
     PublicationRule, Seasonality,
@@ -215,20 +217,16 @@ pub enum SwapCouponConvention {
     Compounded(OvernightLegConvention),
 }
 
-/// FX フォワードの通貨ペア、日程と共通の受渡時刻。
+/// FX フォワードのアウトライトの向きと満期の営業日規則。
+///
+/// スポットラグ、両通貨の休日の和の暦と受渡の時刻は、ペアの [`FxConvention`] から読む。
 #[derive(Debug, Clone, PartialEq, serde::Deserialize)]
 #[serde(deny_unknown_fields)]
 pub struct FxForwardConvention {
     /// レートの向き。
     pub pair: CcyPair,
-    /// 約定日からスポットまでのずれ。
-    pub spot_lag: DateOffset,
-    /// 受渡の暦。
-    pub calendar: CalendarSpec,
     /// 満期に適用する営業日規則。
     pub bdc: BusinessDayConvention,
-    /// 両通貨の受渡時刻。
-    pub settlement: LocalEventTime,
 }
 
 /// 預金の開始日・満期日と利息を決める慣行。
@@ -300,6 +298,8 @@ pub struct RefDataParts {
     pub fallbacks: BTreeMap<IndexId, FallbackDefinition>,
     /// 商品の慣行。
     pub conventions: BTreeMap<ConventionId, Convention>,
+    /// 通貨ペアの現物の受渡慣行。キーは慣行の `pair`（市場の標準の向き）。
+    pub fx_conventions: BTreeMap<CcyPair, FxConvention>,
     /// 証券の定義。
     pub securities: BTreeMap<SecurityId, SecurityDefinition>,
     /// コモディティの定義。
@@ -349,6 +349,8 @@ pub enum RefDataEntry {
     Index(IndexId),
     /// フォールバックの元の指標。
     Fallback(IndexId),
+    /// 通貨ペアの受渡慣行。
+    FxConvention(CcyPair),
 }
 
 impl RefData {
@@ -399,32 +401,11 @@ impl RefData {
                 errors.push((RefDataEntry::Currency(code), error));
             }
         }
-        for (id, security) in &parts.securities {
-            let result = (|| {
-                check_key(id, &security.id, "security")?;
-                lookup(&parts.currencies, &security.currency, "currency")?;
-                if security.kind == SecurityKind::Bond && security.issuer.is_none() {
-                    return Err(Error::InvalidData(format!("bond {id} requires an issuer")));
-                }
-                Ok(())
-            })();
-            if let Err(error) = result {
-                errors.push((RefDataEntry::Security(id.clone()), error));
-            }
-        }
-        for (id, commodity) in &parts.commodities {
-            let result = (|| {
-                check_key(id, &commodity.id, "commodity")?;
-                lookup(&parts.currencies, &commodity.currency, "currency")?;
-                if commodity.unit.trim().is_empty() {
-                    return Err(Error::InvalidData(format!(
-                        "commodity {id} requires a physical unit"
-                    )));
-                }
-                Ok(())
-            })();
-            if let Err(error) = result {
-                errors.push((RefDataEntry::Commodity(id.clone()), error));
+        validate_assets(&parts, &mut errors);
+        for (&pair, convention) in &parts.fx_conventions {
+            if let Err(error) = validate_fx_convention(pair, convention, &parts, &mut use_calendar)
+            {
+                errors.push((RefDataEntry::FxConvention(pair), error));
             }
         }
         let mut rate_futures = BTreeMap::new();
@@ -621,6 +602,24 @@ impl RefData {
         lookup(&self.0.parts.zones, id, "zone")
     }
 
+    /// 通貨ペアの受渡慣行を、向きによらずに返す。慣行の `pair` が市場の標準の向き。
+    ///
+    /// # Errors
+    ///
+    /// どちらの向きの慣行も無ければ `InvalidData`。
+    pub fn fx_convention(&self, pair: CcyPair) -> Result<&FxConvention> {
+        fx_convention(&self.0.parts, pair)
+    }
+
+    /// 通貨ペアの受渡日を数える、両通貨の決済カレンダーの和の式。向きによらず同じ式を返す。
+    ///
+    /// # Errors
+    ///
+    /// どちらかの通貨が未登録なら `InvalidData`。
+    pub fn fx_calendar_spec(&self, pair: CcyPair) -> Result<CalendarSpec> {
+        fx_calendar_spec(&self.0.parts, pair)
+    }
+
     /// 指標の観測点が公表される瞬間と、カレンダー外挿の警告を返す。
     ///
     /// # Errors
@@ -720,6 +719,78 @@ fn lookup<'a, K: Ord + std::fmt::Display, V>(
         .ok_or_else(|| Error::InvalidData(format!("unknown {kind} {id}")))
 }
 
+/// 証券とコモディティの定義の参照を照合し、失敗した項目を `errors` に足す。
+fn validate_assets(parts: &RefDataParts, errors: &mut Vec<(RefDataEntry, Error)>) {
+    for (id, security) in &parts.securities {
+        let result = (|| {
+            check_key(id, &security.id, "security")?;
+            lookup(&parts.currencies, &security.currency, "currency")?;
+            if security.kind == SecurityKind::Bond && security.issuer.is_none() {
+                return Err(Error::InvalidData(format!("bond {id} requires an issuer")));
+            }
+            Ok(())
+        })();
+        if let Err(error) = result {
+            errors.push((RefDataEntry::Security(id.clone()), error));
+        }
+    }
+    for (id, commodity) in &parts.commodities {
+        let result = (|| {
+            check_key(id, &commodity.id, "commodity")?;
+            lookup(&parts.currencies, &commodity.currency, "currency")?;
+            if commodity.unit.trim().is_empty() {
+                return Err(Error::InvalidData(format!(
+                    "commodity {id} requires a physical unit"
+                )));
+            }
+            Ok(())
+        })();
+        if let Err(error) = result {
+            errors.push((RefDataEntry::Commodity(id.clone()), error));
+        }
+    }
+}
+
+/// 表のキーと向き、逆向きの重複、受渡時刻の地域、両通貨の休日の和の暦を照合する。
+fn validate_fx_convention(
+    pair: CcyPair,
+    convention: &FxConvention,
+    parts: &RefDataParts,
+    use_calendar: &mut impl FnMut(&CalendarSpec) -> Result<Calendar>,
+) -> Result<()> {
+    check_key(&pair, &convention.pair, "FX convention")?;
+    // 同じペアの慣行を二つの向きで持つと、どちらを正本にするかが決まらない。
+    if parts.fx_conventions.contains_key(&pair.inverse()) {
+        return Err(Error::InvalidData(format!(
+            "FX conventions for both {pair} and {}",
+            pair.inverse()
+        )));
+    }
+    lookup(
+        &parts.zones,
+        &convention.settlement.zone,
+        "FX settlement zone",
+    )?;
+    use_calendar(&fx_calendar_spec(parts, pair)?)?;
+    Ok(())
+}
+
+fn fx_convention(parts: &RefDataParts, pair: CcyPair) -> Result<&FxConvention> {
+    parts
+        .fx_conventions
+        .get(&pair)
+        .or_else(|| parts.fx_conventions.get(&pair.inverse()))
+        .ok_or_else(|| Error::InvalidData(format!("unknown FX convention {pair}")))
+}
+
+fn fx_calendar_spec(parts: &RefDataParts, pair: CcyPair) -> Result<CalendarSpec> {
+    let [base, quote] = [pair.base(), pair.quote()].map(|ccy| {
+        lookup(&parts.currencies, &ccy, "currency")
+            .map(|definition| &definition.settlement_calendar)
+    });
+    Ok(fx_convention::pair_calendar(base?, quote?))
+}
+
 fn check_key<T: Eq + std::fmt::Display>(key: &T, id: &T, kind: &str) -> Result<()> {
     if key != id {
         return Err(Error::InvalidData(format!(
@@ -775,6 +846,8 @@ fn validate_index_kind(
         IndexKind::FxFixing { pair } => {
             currency(pair.base())?;
             currency(pair.quote())?;
+            // FX fixing の受渡日はペアのスポットラグで決まるので、慣行の無いペアは評価できない。
+            fx_convention(parts, *pair)?;
         }
         IndexKind::AssetPrice { asset, ccy, .. } => {
             currency(*ccy)?;
@@ -891,13 +964,16 @@ fn validate_convention(
                     return Err(Error::InvalidData(format!("{ccy} is not deliverable")));
                 }
             }
-            use_calendar(&fx.calendar)?;
-            lookup(&parts.zones, &fx.settlement.zone, "FX settlement zone")?;
+            fx_convention(parts, fx.pair)?;
         }
     }
     Ok(())
 }
 
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#197 で各 variant を網羅する"
+)]
 fn validate_swap_convention(
     convention: &Convention,
     swap: &SwapConvention,
diff --git a/crates/rhodium-core/src/refdata/fx_convention.rs b/crates/rhodium-core/src/refdata/fx_convention.rs
new file mode 100644
index 0000000000000000000000000000000000000000..f9e96d8bea94fa02634ff7dd8692a411c1c67743
--- /dev/null
+++ b/crates/rhodium-core/src/refdata/fx_convention.rs
@@ -0,0 +1,80 @@
+//! 通貨ペアの現物の受渡慣行。
+
+use crate::{
+    calendar::CalendarSpec,
+    fx::CcyPair,
+    id::CalendarId,
+    tenor::{DateOffset, TenorUnit},
+    time::LocalEventTime,
+};
+
+/// 現物 FX のスポットラグ。両通貨の決済カレンダーの和の営業日で数える。
+///
+/// 閉じた enum にして、スポットまでの短期の quote の連鎖（T+1 は ON だけ、T+2 は TN と ON）を
+/// 網羅して分岐し、市場に無い T+0 や T+5 の慣行を読み込みの境界で拒否する。
+#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
+#[serde(rename_all = "snake_case")]
+pub enum FxSpotLag {
+    /// 翌営業日がスポット（USDCAD など）。
+    OneBusinessDay,
+    /// 二営業日後がスポット（USDJPY、EURUSD など）。
+    TwoBusinessDays,
+}
+
+impl FxSpotLag {
+    /// スポットまでの営業日数。
+    #[must_use]
+    pub const fn business_days(self) -> i16 {
+        match self {
+            Self::OneBusinessDay => 1,
+            Self::TwoBusinessDays => 2,
+        }
+    }
+
+    /// 約定日からスポット日へ進めるずれ。
+    #[must_use]
+    pub const fn offset(self) -> DateOffset {
+        DateOffset::new(self.business_days(), TenorUnit::BusinessDay)
+    }
+}
+
+/// 通貨ペアの現物の受渡慣行。
+///
+/// 評価時点の換算レート X(0) の短期の連鎖、FX fixing の受渡日、`MtM` XCCY のリセットの
+/// 観測日、FX フォワードの校正商品の日程は、すべてこの慣行から決める（正本は一つ）。
+#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
+#[serde(deny_unknown_fields)]
+pub struct FxConvention {
+    /// 市場の標準の向き（USDJPY、EURUSD）。pip とポイントの quote はこの向きで読む。
+    pub pair: CcyPair,
+    /// スポットラグ。
+    pub spot_lag: FxSpotLag,
+    /// 各受渡日の受渡の瞬間。この瞬間を過ぎた当日の受渡はできない（同日受渡の締め）。
+    pub settlement: LocalEventTime,
+    /// 1 pip がレートの 10 の何乗分の一か（USDJPY は 2、EURUSD は 4）。
+    pub pip_exponent: u8,
+}
+
+impl FxConvention {
+    /// 1 pip のレートの大きさ 10^(−`pip_exponent`)。
+    #[must_use]
+    pub fn pip(&self) -> f64 {
+        10_f64.powi(-i32::from(self.pip_exponent))
+    }
+}
+
+/// 両通貨の決済カレンダーの和。ID の順に並べ、向きの違うペアで同じ式（同じキャッシュ）にする。
+pub(super) fn pair_calendar(base: &CalendarId, quote: &CalendarId) -> CalendarSpec {
+    if base == quote {
+        return CalendarSpec::Base(base.clone());
+    }
+    let (first, second) = if base < quote {
+        (base, quote)
+    } else {
+        (quote, base)
+    };
+    CalendarSpec::Union(vec![
+        CalendarSpec::Base(first.clone()),
+        CalendarSpec::Base(second.clone()),
+    ])
+}
diff --git a/crates/rhodium-core/src/refdata/index.rs b/crates/rhodium-core/src/refdata/index.rs
index 9d69162aa5be85ea7ebd9080dd0159e782d24757..fe204b5bc6dc515c9adfe00993dbf65da4b70b7e 100644
--- a/crates/rhodium-core/src/refdata/index.rs
+++ b/crates/rhodium-core/src/refdata/index.rs
@@ -306,6 +306,10 @@ impl ResolvedPublication {
 }
 
 impl IndexKind {
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#197 で各 variant を網羅する"
+    )]
     pub(super) const fn accepts_publication(&self, rule: &PublicationRule) -> bool {
         match self {
             Self::InflationIndex { .. } => matches!(rule, PublicationRule::Monthly { .. }),
diff --git a/crates/rhodium-core/src/refdata/rate_future.rs b/crates/rhodium-core/src/refdata/rate_future.rs
index 418b2cac69d1456e9e3faa9d49e34c7b74ca850d..4c11f45c73311c272115c361ebbe7991e967a69a 100644
--- a/crates/rhodium-core/src/refdata/rate_future.rs
+++ b/crates/rhodium-core/src/refdata/rate_future.rs
@@ -177,6 +177,10 @@ impl RefData {
     /// # Errors
     ///
     /// 清算指標が上場契約を参照しないか、原指標が金利でない場合。
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#197 で各 variant を網羅する"
+    )]
     pub fn rate_future_underlying(
         &self,
         index: &IndexId,
diff --git a/crates/rhodium-core/src/tenor.rs b/crates/rhodium-core/src/tenor.rs
index a7321121df90b1248dbac6c5c6db536bf5941052..ad668be0eed15cf0185e776eb4ba48a9da925fd4 100644
--- a/crates/rhodium-core/src/tenor.rs
+++ b/crates/rhodium-core/src/tenor.rs
@@ -54,6 +54,10 @@ pub struct TenorLength {
 impl TenorLength {
     /// 量と単位を正準化する。
     #[must_use]
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#198 で各 variant を網羅する"
+    )]
     pub const fn new(amount: u16, unit: TenorUnit) -> Self {
         let (amount, unit) = match unit {
             TenorUnit::BusinessDay => (amount, unit),
@@ -77,6 +81,10 @@ impl TenorLength {
         self.unit
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#198 で各 variant を網羅する"
+    )]
     fn expanded(self) -> (u32, TenorUnit) {
         let amount = u32::from(self.amount);
         match self.unit {
@@ -86,6 +94,10 @@ impl TenorLength {
         }
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#198 で各 variant を網羅する"
+    )]
     fn from_expanded(amount: u32, unit: TenorUnit) -> Result<Self> {
         let (amount, unit) = match unit {
             TenorUnit::Day if amount.is_multiple_of(7) => (amount / 7, TenorUnit::Week),
@@ -211,6 +223,10 @@ serde_text!(Tenor);
 impl FromStr for TenorLength {
     type Err = Error;
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#198 で各 variant を網羅する"
+    )]
     fn from_str(text: &str) -> Result<Self> {
         match text.parse()? {
             Tenor::Length(length) => Ok(length),
@@ -290,6 +306,10 @@ impl Frequency {
     /// # Errors
     ///
     /// ゼロまたは ON/TN/SN は [`Error::InvalidData`]。
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#198 で各 variant を網羅する"
+    )]
     pub fn every(step: Tenor) -> Result<Self> {
         match step {
             Tenor::Length(length) if length.amount != 0 => Ok(Self(Some(length))),
diff --git a/crates/rhodium-core/tests/calendar.rs b/crates/rhodium-core/tests/calendar.rs
index 652ff9559829721616e2899c166dcf1cd8a0bba5..9c858fa0e191f73695cbbe83e10fa8ee8a3915f9 100644
--- a/crates/rhodium-core/tests/calendar.rs
+++ b/crates/rhodium-core/tests/calendar.rs
@@ -285,6 +285,70 @@ fn business_day_count_respects_bit_boundaries_and_half_open_ranges() -> Result<(
     Ok(())
 }
 
+#[test]
+fn business_days_across_the_known_range_match_day_by_day_lookups() -> Result<()> {
+    // 既知範囲の内側はビット集合で、外側は式の評価（週末だけの外挿）で数える。区間が既知範囲の
+    // 前から後ろまでをまたいでも、営業日数と警告が一日ずつの照会の和と一致することを確かめる。
+    let cal = calendar(
+        "JPTO",
+        base(
+            "2027-01-04",
+            "2027-02-26",
+            &["2027-01-11", "2027-02-11"],
+            Extrapolation::WeekendsOnly,
+        )?,
+    )?;
+    let start: Date = "2026-12-28".parse()?;
+    let end: Date = "2027-03-08".parse()?;
+    let mut business = 0;
+    let mut warnings = Vec::new();
+    for days in 0..end.days_since(start) {
+        let (open, notices) = cal.is_business_day(start.checked_add_days(days)?)?;
+        business += u32::from(open);
+        warnings.extend(notices);
+    }
+    warnings.sort_unstable();
+    warnings.dedup();
+    assert_eq!(cal.business_days_between(start, end)?, (business, warnings));
+    // 既知範囲の最初と最後の日はビット集合から読み、外挿の警告を出さない。
+    assert_eq!(cal.is_business_day("2027-01-04".parse()?)?, (true, vec![]));
+    assert_eq!(cal.is_business_day("2027-02-26".parse()?)?, (true, vec![]));
+    // 範囲の外の平日は営業日で、外挿の警告を一つ出す。
+    let (open, notices) = cal.is_business_day("2027-03-01".parse()?)?;
+    assert!(open);
+    assert_eq!(notices.len(), 1);
+    Ok(())
+}
+
+#[test]
+fn calendar_with_disjoint_known_ranges_evaluates_every_day() -> Result<()> {
+    // 和の要素の既知範囲が交わらなければ共通の既知範囲は空で、ビット集合も空になる。その場合も
+    // 照会は式の評価で答え、範囲の外の日に Strict の基本暦は MissingData を返す。
+    let early: CalendarId = "EARLY".parse()?;
+    let late: CalendarId = "LATE".parse()?;
+    let bases = BTreeMap::from([
+        (
+            early.clone(),
+            base("2027-01-01", "2027-01-31", &[], Extrapolation::WeekendsOnly)?,
+        ),
+        (
+            late.clone(),
+            base("2027-03-01", "2027-03-31", &[], Extrapolation::Strict)?,
+        ),
+    ]);
+    let union = Calendar::compile(
+        &CalendarSpec::Union(vec![CalendarSpec::Base(early), CalendarSpec::Base(late)]),
+        &bases,
+    )?;
+    let (open, _) = union.is_business_day("2027-03-02".parse()?)?;
+    assert!(open);
+    assert!(matches!(
+        union.is_business_day("2027-02-02".parse()?),
+        Err(Error::MissingData(_))
+    ));
+    Ok(())
+}
+
 #[test]
 fn advance_respects_business_days_and_month_end() -> Result<()> {
     let cal = calendar(
diff --git a/crates/rhodium-core/tests/day_count.rs b/crates/rhodium-core/tests/day_count.rs
index 1565100c6f9708e28b186689b522b51f80b989bc..07e478e8b09f335811f29cd07a3fe340124446f7 100644
--- a/crates/rhodium-core/tests/day_count.rs
+++ b/crates/rhodium-core/tests/day_count.rs
@@ -173,6 +173,10 @@ fn all_day_counts_match_hand_calculations() -> Result<()> {
 }
 
 #[test]
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#196 で各 variant を網羅する"
+)]
 fn day_count_fractions_agree_with_published_oracles() -> Result<()> {
     // Strata の公開検算（固定コミット 910a144fc49532c1411ed95137556476d622c85c）の期間と値。
     // BUS/252 だけは QuantLib v1.40 の最初の例で、金曜から月曜までの区間で金曜だけを数える。
diff --git a/crates/rhodium-core/tests/money_and_quantities.rs b/crates/rhodium-core/tests/money_and_quantities.rs
index bfcde0a61117d3797707935d632c71a8bfc2e95e..fb81d7a29534f744bc6173efc4b2c15d6b11b08f 100644
--- a/crates/rhodium-core/tests/money_and_quantities.rs
+++ b/crates/rhodium-core/tests/money_and_quantities.rs
@@ -15,9 +15,9 @@ use rhodium_core::{
     id::{CalendarId, EntityId, IndexId},
     money::{Money, MoneyBag},
     quantity::{
-        Clean, Compounded, Continuous, Correlation, DiscountFactor, ParPrice, Probability, Rate,
-        RateFuturesPrice, RecoveryRate, RiskyAnnuity, RunningSpread, ShiftedLognormalVol, Simple,
-        Spread, SurvivalProbability, YearFraction, accrue,
+        Clean, Compounded, Continuous, Correlation, DiscountFactor, NormalVol, ParPrice,
+        Probability, Rate, RateFuturesPrice, RecoveryRate, RiskyAnnuity, RunningSpread,
+        ShiftedLognormalVol, Simple, Spread, SurvivalProbability, YearFraction, accrue,
     },
     rounding::{Rounding, RoundingMode},
 };
@@ -182,6 +182,17 @@ fn shift_mismatch() -> Result<()> {
     Ok(())
 }
 
+#[test]
+fn signed_zero_volatilities_compare_equal() -> Result<()> {
+    // 比較は全順序だが、構築時に -0.0 を 0.0 にそろえるので、符号付きの零は同じボラになる。
+    let negative = ShiftedLognormalVol::new(-0.0, 0.01)?;
+    let positive = ShiftedLognormalVol::new(0.0, 0.01)?;
+    assert_eq!(negative.compare(positive)?, Ordering::Equal);
+    assert!(negative.value().is_sign_positive());
+    assert!(NormalVol::new(-0.0)?.value().is_sign_positive());
+    Ok(())
+}
+
 #[test]
 fn fx_direction() -> Result<()> {
     // USDJPY は USD から JPY だけを換算する。逆向きの換算には反転したレートを渡す。
diff --git a/crates/rhodium-core/tests/refdata_and_fixing.rs b/crates/rhodium-core/tests/refdata_and_fixing.rs
index d8277dcdc3fb1fbf859181f424b667968d84b02f..bbbb830b70317d0d5bb36b162f03ae4a7c03ce16 100644
--- a/crates/rhodium-core/tests/refdata_and_fixing.rs
+++ b/crates/rhodium-core/tests/refdata_and_fixing.rs
@@ -12,12 +12,13 @@ use rhodium_core::{
         CreditEventRecord, FixingResolution, FixingStore, FixingStoreParts, ObservationPeriod,
         ObservationPoint, YearMonth,
     },
+    fx::CcyPair,
     id::{CalendarId, CommodityId, CreditKey, IndexId},
     quantity::RecoveryRate,
     refdata::{
-        CommodityDefinition, CurrencyDefinition, DailyPublication, Deliverability, IndexDefinition,
-        IndexKind, IndexStatus, InflationInterpolation, PriceKind, PublicationRule, RefData,
-        RefDataParts, Seasonality,
+        CommodityDefinition, CurrencyDefinition, DailyPublication, Deliverability, FxConvention,
+        FxSpotLag, IndexDefinition, IndexKind, IndexStatus, InflationInterpolation, PriceKind,
+        PublicationRule, RefData, RefDataEntry, RefDataParts, Seasonality,
     },
     rounding::{Rounding, RoundingMode},
     time::{Instant, LocalEventTime, UtcOffset, ZoneRules},
@@ -66,6 +67,8 @@ fn parts() -> Result<RefDataParts> {
             },
         );
     }
+    // USDJPY の FX fixing の受渡日は、ペアの慣行のスポットラグで決まる。
+    add_usd_jpy(&mut parts)?;
     let calendar = CalendarSpec::Base(calendar);
     let definitions = [
         (
@@ -122,6 +125,24 @@ fn parts() -> Result<RefDataParts> {
     Ok(parts)
 }
 
+/// USDJPY の受渡慣行（T+2、ニューヨーク 17 時、pip 0.01 円）を足す。
+fn add_usd_jpy(parts: &mut RefDataParts) -> Result<()> {
+    let usd_jpy: CcyPair = "USDJPY".parse()?;
+    parts.fx_conventions.insert(
+        usd_jpy,
+        FxConvention {
+            pair: usd_jpy,
+            spot_lag: FxSpotLag::TwoBusinessDays,
+            settlement: LocalEventTime {
+                time: "17:00:00".parse()?,
+                zone: "America/New_York".parse()?,
+            },
+            pip_exponent: 2,
+        },
+    );
+    Ok(())
+}
+
 fn daily(calendar: CalendarSpec, business_days: u16, time: &str) -> Result<PublicationRule> {
     Ok(PublicationRule::Daily(DailyPublication {
         business_days,
@@ -519,3 +540,139 @@ fn zone_lookup_returns_registered_rules_and_rejects_unknown_zones() -> Result<()
     ));
     Ok(())
 }
+
+/// USDJPY の受渡慣行に、休日の違う決済カレンダーの USD と JPY を足した参照データ。
+fn fx_parts() -> Result<RefDataParts> {
+    let mut parts = parts()?;
+    let tokyo: CalendarId = "JPTO".parse()?;
+    parts.calendars.insert(
+        tokyo.clone(),
+        BaseCalendar::new(
+            "2026-01-01".parse()?..="2026-12-31".parse()?,
+            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
+            vec!["2026-11-03".parse()?],
+            Extrapolation::Strict,
+        )?,
+    );
+    if let Some(jpy) = parts.currencies.get_mut(&Currency::JPY) {
+        jpy.settlement_calendar = tokyo;
+    }
+    Ok(parts)
+}
+
+#[test]
+fn fx_convention_lookup_accepts_either_direction() -> Result<()> {
+    // 慣行は市場の標準の向き（USDJPY）で一つだけ持ち、JPYUSD で引いても同じ慣行を返す。
+    // 受渡日を数える暦は両通貨の決済カレンダーの和で、向きによらず同じ式（同じキャッシュ）になる。
+    let refdata = RefData::new(fx_parts()?)?;
+    let usd_jpy: CcyPair = "USDJPY".parse()?;
+    let convention = refdata.fx_convention(usd_jpy)?;
+    assert!(ptr_eq(
+        convention,
+        refdata.fx_convention(usd_jpy.inverse())?
+    ));
+    assert_eq!(convention.pair, usd_jpy);
+    assert_eq!(convention.spot_lag.business_days(), 2);
+    assert_close!(convention.pip(), 0.01, 1e-15, 0.0);
+    let spec = refdata.fx_calendar_spec(usd_jpy)?;
+    assert_eq!(spec, refdata.fx_calendar_spec(usd_jpy.inverse())?);
+    assert_eq!(
+        spec,
+        CalendarSpec::Union(vec![
+            CalendarSpec::Base("JPTO".parse()?),
+            CalendarSpec::Base("USNY".parse()?),
+        ])
+    );
+    // 和の暦は、東京の祝日（11 月 3 日）とニューヨークの祝日（11 月 23 日）の両方で休む。
+    let calendar = refdata.calendar(&spec)?;
+    for holiday in ["2026-11-03", "2026-11-23"] {
+        assert!(!calendar.is_business_day(holiday.parse()?)?.0);
+    }
+    assert!(matches!(
+        refdata.fx_convention("EURUSD".parse()?),
+        Err(Error::InvalidData(message)) if message.contains("unknown FX convention EURUSD")
+    ));
+    Ok(())
+}
+
+fn ptr_eq<T>(lhs: &T, rhs: &T) -> bool {
+    std::ptr::eq(lhs, rhs)
+}
+
+#[test]
+fn refdata_rejects_fx_conventions_in_both_directions() -> Result<()> {
+    // 同じペアを両方の向きで持つと、pip とポイントをどちらの向きで読むかが決まらない。
+    let mut parts = fx_parts()?;
+    let usd_jpy: CcyPair = "USDJPY".parse()?;
+    let mut inverse = parts.fx_conventions[&usd_jpy].clone();
+    inverse.pair = usd_jpy.inverse();
+    parts.fx_conventions.insert(inverse.pair, inverse);
+    let Err(errors) = RefData::new_with_diagnostics(parts) else {
+        return Err(Error::Internal("both directions were accepted".into()));
+    };
+    assert!(errors.iter().all(|(entry, error)| matches!(
+        (entry, error),
+        (RefDataEntry::FxConvention(_), Error::InvalidData(message))
+            if message.contains("FX conventions for both")
+    )));
+    assert_eq!(errors.len(), 2);
+    Ok(())
+}
+
+#[test]
+fn refdata_requires_fx_convention_for_fx_fixing_index() -> Result<()> {
+    // FX fixing の受渡日はペアのスポットラグで決まるので、慣行の無いペアの fixing 指標は
+    // 参照データを作るときに拒否する。どちらの向きの慣行でもよい。
+    let mut parts = fx_parts()?;
+    parts.fx_conventions.clear();
+    let Err(errors) = RefData::new_with_diagnostics(parts) else {
+        return Err(Error::Internal(
+            "FX fixing without a convention was accepted".into(),
+        ));
+    };
+    assert!(errors.iter().any(|(entry, error)| matches!(
+        (entry, error),
+        (RefDataEntry::Index(id), Error::InvalidData(message))
+            if id.as_str() == "USDJPY-WMR" && message.contains("unknown FX convention USDJPY")
+    )));
+    let mut parts = fx_parts()?;
+    let usd_jpy: CcyPair = "USDJPY".parse()?;
+    let mut convention = parts
+        .fx_conventions
+        .remove(&usd_jpy)
+        .ok_or_else(|| Error::Internal("missing USDJPY convention".into()))?;
+    convention.pair = usd_jpy.inverse();
+    parts.fx_conventions.insert(convention.pair, convention);
+    RefData::new(parts)?;
+    Ok(())
+}
+
+#[test]
+fn refdata_rejects_fx_convention_key_mismatch_and_unknown_zone() -> Result<()> {
+    // 表のキーと慣行の向きが違う定義、受渡時刻の地域が未登録の定義を、慣行の項目として報告する。
+    let usd_jpy: CcyPair = "USDJPY".parse()?;
+    let mut parts = fx_parts()?;
+    let convention = parts.fx_conventions[&usd_jpy].clone();
+    parts.fx_conventions.insert("EURJPY".parse()?, convention);
+    let Err(errors) = RefData::new_with_diagnostics(parts) else {
+        return Err(Error::Internal("mismatched key was accepted".into()));
+    };
+    assert!(errors.iter().any(|(entry, _)| matches!(
+        entry,
+        RefDataEntry::FxConvention(pair) if pair.to_string() == "EURJPY"
+    )));
+
+    let mut parts = fx_parts()?;
+    if let Some(convention) = parts.fx_conventions.get_mut(&usd_jpy) {
+        convention.settlement.zone = "Asia/Tokyo".parse()?;
+    }
+    let Err(errors) = RefData::new_with_diagnostics(parts) else {
+        return Err(Error::Internal("unknown zone was accepted".into()));
+    };
+    assert!(errors.iter().any(|(entry, error)| matches!(
+        (entry, error),
+        (RefDataEntry::FxConvention(pair), Error::InvalidData(message))
+            if *pair == usd_jpy && message.contains("Asia/Tokyo")
+    )));
+    Ok(())
+}
diff --git a/crates/rhodium-core/tests/ui/clean_price_settlement.json b/crates/rhodium-core/tests/ui/clean_price_settlement.json
new file mode 100644
index 0000000000000000000000000000000000000000..db16badddbfe4b784b3fc3f849e81b5f0453903e
--- /dev/null
+++ b/crates/rhodium-core/tests/ui/clean_price_settlement.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0599",
+    "message": "no method named `settlement` found for struct `ParPrice<Clean>` in the current scope",
+    "spans": [
+      {
+        "file": "tests/ui/clean_price_settlement.rs",
+        "line": 11,
+        "column": 19
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-core/tests/ui/clean_price_settlement.stderr b/crates/rhodium-core/tests/ui/clean_price_settlement.stderr
deleted file mode 100644
index 2da01a9281ae120fb07d42433d7e5b36d4b6307a..0000000000000000000000000000000000000000
--- a/crates/rhodium-core/tests/ui/clean_price_settlement.stderr
+++ /dev/null
@@ -1,7 +0,0 @@
-error[E0599]: no method named `settlement` found for struct `ParPrice<Clean>` in the current scope
-  --> tests/ui/clean_price_settlement.rs:11:19
-   |
-11 |     let _ = clean.settlement(Money::new(Currency::JPY, 1000.0));
-   |                   ^^^^^^^^^^ method not found in `ParPrice<Clean>`
-   |
-   = note: the method was found for `ParPrice<Dirty>`
diff --git a/crates/rhodium-core/tests/ui/date_as_instant.json b/crates/rhodium-core/tests/ui/date_as_instant.json
new file mode 100644
index 0000000000000000000000000000000000000000..600dcc6b2faacaea6fbd26d99bafda0d4bcfc2a0
--- /dev/null
+++ b/crates/rhodium-core/tests/ui/date_as_instant.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0308",
+    "message": "mismatched types",
+    "spans": [
+      {
+        "file": "tests/ui/date_as_instant.rs",
+        "line": 7,
+        "column": 32
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-core/tests/ui/date_as_instant.stderr b/crates/rhodium-core/tests/ui/date_as_instant.stderr
deleted file mode 100644
index 7c37c8208236501088e7131d92525ecbaf660efd..0000000000000000000000000000000000000000
--- a/crates/rhodium-core/tests/ui/date_as_instant.stderr
+++ /dev/null
@@ -1,13 +0,0 @@
-error[E0308]: mismatched types
- --> tests/ui/date_as_instant.rs:7:32
-  |
-7 |     let _ = model_time(origin, at);
-  |             ----------         ^^ expected `Instant`, found `Date`
-  |             |
-  |             arguments to this function are incorrect
-  |
-note: function defined here
- --> src/time.rs
-  |
-  | pub fn model_time(origin: Instant, at: Instant) -> ModelTime {
-  |        ^^^^^^^^^^
diff --git a/crates/rhodium-core/tests/ui/date_plus_tenor.json b/crates/rhodium-core/tests/ui/date_plus_tenor.json
new file mode 100644
index 0000000000000000000000000000000000000000..77fcb0450cf6bbf2ef0c995372a58d8f7432ad82
--- /dev/null
+++ b/crates/rhodium-core/tests/ui/date_plus_tenor.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0369",
+    "message": "cannot add `Tenor` to `rhodium_core::date::Date`",
+    "spans": [
+      {
+        "file": "tests/ui/date_plus_tenor.rs",
+        "line": 4,
+        "column": 18
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-core/tests/ui/date_plus_tenor.stderr b/crates/rhodium-core/tests/ui/date_plus_tenor.stderr
deleted file mode 100644
index 8f30cfe902b95a4a9d408a263a12b2de5457fff2..0000000000000000000000000000000000000000
--- a/crates/rhodium-core/tests/ui/date_plus_tenor.stderr
+++ /dev/null
@@ -1,13 +0,0 @@
-error[E0369]: cannot add `Tenor` to `rhodium_core::date::Date`
- --> tests/ui/date_plus_tenor.rs:4:18
-  |
-4 |     let _ = date + tenor;
-  |             ---- ^ ----- Tenor
-  |             |
-  |             rhodium_core::date::Date
-  |
-note: `rhodium_core::date::Date` does not implement `Add<Tenor>`
- --> src/date.rs
-  |
-  | pub struct Date(NaiveDate);
-  | ^^^^^^^^^^^^^^^ `rhodium_core::date::Date` is defined in another crate
diff --git a/crates/rhodium-core/tests/ui/futures_price_as_rate.json b/crates/rhodium-core/tests/ui/futures_price_as_rate.json
new file mode 100644
index 0000000000000000000000000000000000000000..0cfaf07f36b35a73eeadfa43125d7a6cc85d843c
--- /dev/null
+++ b/crates/rhodium-core/tests/ui/futures_price_as_rate.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0308",
+    "message": "mismatched types",
+    "spans": [
+      {
+        "file": "tests/ui/futures_price_as_rate.rs",
+        "line": 11,
+        "column": 20
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-core/tests/ui/futures_price_as_rate.stderr b/crates/rhodium-core/tests/ui/futures_price_as_rate.stderr
deleted file mode 100644
index 95d82859798b3c87428d58e35a5e8c1c6296ed00..0000000000000000000000000000000000000000
--- a/crates/rhodium-core/tests/ui/futures_price_as_rate.stderr
+++ /dev/null
@@ -1,15 +0,0 @@
-error[E0308]: mismatched types
-  --> tests/ui/futures_price_as_rate.rs:11:20
-   |
-11 |     let _ = accrue(&price, &yf);
-   |             ------ ^^^^^^ expected `&Rate<_>`, found `&RateFuturesPrice`
-   |             |
-   |             arguments to this function are incorrect
-   |
-   = note: expected reference `&Rate<_>`
-              found reference `&RateFuturesPrice`
-note: function defined here
-  --> src/quantity/rates.rs
-   |
-   | pub fn accrue<K: RateKind>(rate: &Rate<K>, fraction: &YearFraction) -> Result<f64> {
-   |        ^^^^^^
diff --git a/crates/rhodium-core/tests/ui/fx_spot_as_today.json b/crates/rhodium-core/tests/ui/fx_spot_as_today.json
new file mode 100644
index 0000000000000000000000000000000000000000..08561b70e314268ae05a0149fb05c3aa2118bd4f
--- /dev/null
+++ b/crates/rhodium-core/tests/ui/fx_spot_as_today.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0599",
+    "message": "no method named `convert` found for struct `FxRate<Spot>` in the current scope",
+    "spans": [
+      {
+        "file": "tests/ui/fx_spot_as_today.rs",
+        "line": 6,
+        "column": 18
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-core/tests/ui/fx_spot_as_today.stderr b/crates/rhodium-core/tests/ui/fx_spot_as_today.stderr
deleted file mode 100644
index 0e81119b908527a074f42cda5f518862160f1d8d..0000000000000000000000000000000000000000
--- a/crates/rhodium-core/tests/ui/fx_spot_as_today.stderr
+++ /dev/null
@@ -1,11 +0,0 @@
-error[E0599]: no method named `convert` found for struct `FxRate<Spot>` in the current scope
- --> tests/ui/fx_spot_as_today.rs:6:18
-  |
-6 |     let _ = spot.convert(Money::new(Currency::USD, 2.0));
-  |                  ^^^^^^^
-  |
-help: there is a method `invert` with a similar name, but with different arguments
- --> src/fx.rs
-  |
-  |     pub fn invert(self) -> Result<Self> {
-  |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
diff --git a/crates/rhodium-core/tests/ui/normal_plus_lognormal.json b/crates/rhodium-core/tests/ui/normal_plus_lognormal.json
new file mode 100644
index 0000000000000000000000000000000000000000..8fa6121a5d135693f846530227414481f0d32642
--- /dev/null
+++ b/crates/rhodium-core/tests/ui/normal_plus_lognormal.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0308",
+    "message": "mismatched types",
+    "spans": [
+      {
+        "file": "tests/ui/normal_plus_lognormal.rs",
+        "line": 10,
+        "column": 22
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-core/tests/ui/normal_plus_lognormal.stderr b/crates/rhodium-core/tests/ui/normal_plus_lognormal.stderr
deleted file mode 100644
index 7cd432a339bf7568eb9bb74dbfd9c5328038272d..0000000000000000000000000000000000000000
--- a/crates/rhodium-core/tests/ui/normal_plus_lognormal.stderr
+++ /dev/null
@@ -1,7 +0,0 @@
-error[E0308]: mismatched types
-  --> tests/ui/normal_plus_lognormal.rs:10:22
-   |
-10 |     let _ = normal < lognormal;
-   |             ------   ^^^^^^^^^ expected `NormalVol`, found `LognormalVol`
-   |             |
-   |             expected because this is `NormalVol`
diff --git a/crates/rhodium-core/tests/ui/rate_kind_mismatch.json b/crates/rhodium-core/tests/ui/rate_kind_mismatch.json
new file mode 100644
index 0000000000000000000000000000000000000000..15bab94fc094bd7f78f2c94c1715b6afa11856a7
--- /dev/null
+++ b/crates/rhodium-core/tests/ui/rate_kind_mismatch.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0277",
+    "message": "cannot add `Rate<Continuous>` to `Rate<Simple>`",
+    "spans": [
+      {
+        "file": "tests/ui/rate_kind_mismatch.rs",
+        "line": 11,
+        "column": 20
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-core/tests/ui/rate_kind_mismatch.stderr b/crates/rhodium-core/tests/ui/rate_kind_mismatch.stderr
deleted file mode 100644
index d84270cafbb719b7e0d23f4330689462cad5e46e..0000000000000000000000000000000000000000
--- a/crates/rhodium-core/tests/ui/rate_kind_mismatch.stderr
+++ /dev/null
@@ -1,13 +0,0 @@
-error[E0277]: cannot add `Rate<Continuous>` to `Rate<Simple>`
-  --> tests/ui/rate_kind_mismatch.rs:11:20
-   |
-11 |     let _ = simple + continuous;
-   |                    ^ no implementation for `Rate<Simple> + Rate<Continuous>`
-   |
-help: the trait `Add<Rate<Continuous>>` is not implemented for `Rate<Simple>`
-      but trait `Add<Spread>` is implemented for it
-  --> src/quantity/rates.rs
-   |
-   | impl Add<Spread> for Rate<Simple> {
-   | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
-   = help: for that trait implementation, expected `Spread`, found `Rate<Continuous>`
diff --git a/crates/rhodium-core/tests/ui/rate_times_year_fraction.json b/crates/rhodium-core/tests/ui/rate_times_year_fraction.json
new file mode 100644
index 0000000000000000000000000000000000000000..35bf8e3fadae2a843a2a3a194fb14895f5098e12
--- /dev/null
+++ b/crates/rhodium-core/tests/ui/rate_times_year_fraction.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0369",
+    "message": "cannot multiply `Rate<Simple>` by `YearFraction`",
+    "spans": [
+      {
+        "file": "tests/ui/rate_times_year_fraction.rs",
+        "line": 11,
+        "column": 18
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-core/tests/ui/rate_times_year_fraction.stderr b/crates/rhodium-core/tests/ui/rate_times_year_fraction.stderr
deleted file mode 100644
index 1b4bc68c598b91b0898c707f2f07efa1fc6ab21e..0000000000000000000000000000000000000000
--- a/crates/rhodium-core/tests/ui/rate_times_year_fraction.stderr
+++ /dev/null
@@ -1,13 +0,0 @@
-error[E0369]: cannot multiply `Rate<Simple>` by `YearFraction`
-  --> tests/ui/rate_times_year_fraction.rs:11:18
-   |
-11 |     let _ = rate * yf;
-   |             ---- ^ -- YearFraction
-   |             |
-   |             Rate<Simple>
-   |
-note: `Rate<Simple>` does not implement `Mul<YearFraction>`
-  --> src/quantity/rates.rs
-   |
-   | pub struct Rate<K> {
-   | ^^^^^^^^^^^^^^^^^^ `Rate<Simple>` is defined in another crate
diff --git a/crates/rhodium-core/tests/ui/survival_as_discount.json b/crates/rhodium-core/tests/ui/survival_as_discount.json
new file mode 100644
index 0000000000000000000000000000000000000000..559ea00fcc356679664701e6a4efcbf55dcb121a
--- /dev/null
+++ b/crates/rhodium-core/tests/ui/survival_as_discount.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0308",
+    "message": "mismatched types",
+    "spans": [
+      {
+        "file": "tests/ui/survival_as_discount.rs",
+        "line": 10,
+        "column": 29
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-core/tests/ui/survival_as_discount.stderr b/crates/rhodium-core/tests/ui/survival_as_discount.stderr
deleted file mode 100644
index 0a6ab13d50e4bdaa55461daf930b863f6a76368e..0000000000000000000000000000000000000000
--- a/crates/rhodium-core/tests/ui/survival_as_discount.stderr
+++ /dev/null
@@ -1,15 +0,0 @@
-error[E0308]: mismatched types
-  --> tests/ui/survival_as_discount.rs:10:29
-   |
-10 |     let _ = start.growth_to(survival)?;
-   |                   --------- ^^^^^^^^ expected `DiscountFactor`, found `SurvivalProbability`
-   |                   |
-   |                   arguments to this method are incorrect
-   |
-   = note: expected struct `DiscountFactor`
-              found struct `SurvivalProbability`
-note: method defined here
-  --> src/quantity.rs
-   |
-   |     pub fn growth_to(self, end: Self) -> Result<R> {
-   |            ^^^^^^^^^
diff --git a/crates/rhodium-core/tests/ui/upfront_plus_running.json b/crates/rhodium-core/tests/ui/upfront_plus_running.json
new file mode 100644
index 0000000000000000000000000000000000000000..a90e5d42175929acc40844f14bfa0bfe6ceb9bae
--- /dev/null
+++ b/crates/rhodium-core/tests/ui/upfront_plus_running.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0369",
+    "message": "cannot add `RunningSpread` to `UpfrontFraction`",
+    "spans": [
+      {
+        "file": "tests/ui/upfront_plus_running.rs",
+        "line": 8,
+        "column": 41
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-core/tests/ui/upfront_plus_running.stderr b/crates/rhodium-core/tests/ui/upfront_plus_running.stderr
deleted file mode 100644
index 38fcfc491c6f9eb912ec7efb619a8d400109d58c..0000000000000000000000000000000000000000
--- a/crates/rhodium-core/tests/ui/upfront_plus_running.stderr
+++ /dev/null
@@ -1,23 +0,0 @@
-error[E0369]: cannot add `RunningSpread` to `UpfrontFraction`
- --> tests/ui/upfront_plus_running.rs:8:41
-  |
-8 |     let _ = UpfrontFraction::new(0.02)? + RunningSpread::new(0.01)?;
-  |             --------------------------- ^ ------------------------- RunningSpread
-  |             |
-  |             UpfrontFraction
-  |
-note: `UpfrontFraction` does not implement `Add<RunningSpread>`
- --> src/quantity/prices.rs
-  |
-  |               pub struct $name(f64);
-  |               ^^^^^^^^^^^^^^^^^^^^^^ `UpfrontFraction` is defined in another crate
-...
-  | / scalar_quantity! {
-  | |     /// CDS のランニング・スプレッド。
-  | |     RunningSpread, "running spread", |_| true;
-  | |     /// CDS の初期に支払う、元本に対する比率。
-... |
-  | |     Correlation, "correlation", |value| (-1.0..=1.0).contains(&value);
-  | | }
-  | |_- in this macro invocation
-  = note: this error originates in the macro `scalar_quantity` (in Nightly builds, run with -Z macro-backtrace for more info)
diff --git a/crates/rhodium-market/src/curve.rs b/crates/rhodium-market/src/curve.rs
index e6c658ed788775ed1952b03761429fcaed70a03c..a864f20de395224d810e8c35463f0032632da51b 100644
--- a/crates/rhodium-market/src/curve.rs
+++ b/crates/rhodium-market/src/curve.rs
@@ -11,6 +11,7 @@ use rhodium_core::{
     Error, Result,
     currency::Currency,
     fx::CcyPair,
+    quantity::Spread,
     time::{Instant, model_time},
 };
 use scioxide_ad::Real;
@@ -97,15 +98,22 @@ pub enum Extrapolation {
 }
 
 /// 割引因子または予想用の擬似割引因子の対数を作る部品。
-#[derive(Debug, Clone)]
+///
+/// ファイルの形式はこの型の serde の表現そのもので、瞬間は UTC の日時、`Flat` の率は
+/// 接尾辞付きの文字列（`"25bp"`、`"0.25%"`）で書く。外貨担保曲線（`FxImplied`）は市場の定義の
+/// FX から作るので、ファイルの曲線には書けない。
+#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
+#[serde(rename_all = "snake_case", deny_unknown_fields)]
 pub enum CurveExpr {
     /// 原点を固定し、ln DF を線形補間する節点。
     Nodes {
-        /// θ の部品名。
+        /// θ の部品名。省略すると `nodes`。
+        #[serde(default = "nodes_part")]
         part: String,
         /// UTC のピラー。
         pillars: Vec<Instant>,
         /// 最後のピラー以後の規則。
+        #[serde(default)]
         extrapolation: Extrapolation,
     },
     /// 同じ通貨の別の曲線の ln DF。
@@ -123,7 +131,7 @@ pub enum CurveExpr {
         body: Box<Self>,
     },
     /// 一定の瞬間フォワードで得る ln DF。
-    Flat(f64),
+    Flat(#[serde(deserialize_with = "scaled_rate")] f64),
     /// 格子の区間ごとに候補の平均フォワードの最大値を取り、区間内ではその値を一定とする。
     /// 同点の候補には感応度を等分し、格子の最終時刻より先の点は拒否する。
     MaxForward {
@@ -156,6 +164,7 @@ pub enum CurveExpr {
         dates: Vec<Instant>,
     },
     /// 同じ担保の反対通貨の割引曲線と FX の対数フォワード比から作る曲線。
+    #[serde(skip_deserializing)]
     FxImplied {
         /// θ の部品名。要素は 基準通貨と決済通貨の ln DF の差。
         part: String,
@@ -170,6 +179,22 @@ pub enum CurveExpr {
     },
 }
 
+fn nodes_part() -> String {
+    "nodes".into()
+}
+
+/// 率を接尾辞付きの文字列（`"25bp"`、`"0.25%"`）で読む。裸の数は単位の書き忘れとして拒否する。
+///
+/// 書式は取引や担保の率と同じ [`Spread`] の規則（符号つきの十進、空白なし）に揃える。
+fn scaled_rate<'de, D: serde::Deserializer<'de>>(
+    deserializer: D,
+) -> std::result::Result<f64, D::Error> {
+    let text = <String as serde::Deserialize>::deserialize(deserializer)?;
+    text.parse::<Spread>()
+        .map(Spread::value)
+        .map_err(|error| serde::de::Error::custom(format!("rate {text:?}: {error}")))
+}
+
 #[derive(Debug, Clone, Copy)]
 pub(crate) enum ParameterSpec<'a> {
     Nodes(&'a [Instant]),
@@ -317,6 +342,19 @@ impl CurveSpec {
         Ok(blocks)
     }
 
+    /// この曲線が定義する会合日の階段（`Steps`）のブロック。参照先のブロックは含まない。
+    pub fn step_blocks(&self) -> impl Iterator<Item = BlockId> {
+        let mut blocks = Vec::new();
+        // 構築時に部品名と日付を検査済み。
+        let _: Result<()> = self.expr.visit_blocks(&mut |part, spec| {
+            if matches!(spec, ParameterSpec::Steps(_)) {
+                blocks.push(BlockId::for_part(self.id.clone(), part));
+            }
+            Ok(())
+        });
+        blocks.into_iter()
+    }
+
     /// この曲線が定義する θ のブロックと要素数。参照先のブロックは含まない。
     pub fn parameter_blocks(&self) -> impl Iterator<Item = (BlockId, usize)> {
         let mut blocks = Vec::new();
@@ -328,6 +366,12 @@ impl CurveSpec {
         blocks.into_iter()
     }
 
+    /// 曲線を作る部品の式。
+    #[must_use]
+    pub const fn expr(&self) -> &CurveExpr {
+        &self.expr
+    }
+
     /// 曲線の ID。
     #[must_use]
     pub const fn id(&self) -> &CurveId {
@@ -348,6 +392,10 @@ impl CurveSpec {
 
     /// 単一の `Nodes` 曲線の UTC ピラー。合成曲線では空。
     #[must_use]
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#200 で各 variant を網羅する"
+    )]
     pub fn pillars(&self) -> &[Instant] {
         match &self.expr {
             CurveExpr::Nodes { part, pillars, .. } if part == "nodes" => pillars,
@@ -357,6 +405,10 @@ impl CurveSpec {
 
     /// 単一の `Nodes` 曲線の外挿規則。合成曲線では部品ごとに指定する。
     #[must_use]
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#200 で各 variant を網羅する"
+    )]
     pub const fn extrapolation(&self) -> Extrapolation {
         match &self.expr {
             CurveExpr::Nodes { extrapolation, .. } => *extrapolation,
diff --git a/crates/rhodium-market/src/lib.rs b/crates/rhodium-market/src/lib.rs
index a424577cca72533106413c99d44ecb1b10d258c0..6347df635d6a21d4fe56589b4b7ae0a4c2269929 100644
--- a/crates/rhodium-market/src/lib.rs
+++ b/crates/rhodium-market/src/lib.rs
@@ -8,3 +8,4 @@ pub mod market;
 pub mod quote;
 pub mod structure;
 pub mod view;
+pub mod vol;
diff --git a/crates/rhodium-market/src/structure.rs b/crates/rhodium-market/src/structure.rs
index 1c2726185172e5bcdbe9705e8cac914cf7f02c5b..110295888a4edcbb843c30f8ca319103b8f61c85 100644
--- a/crates/rhodium-market/src/structure.rs
+++ b/crates/rhodium-market/src/structure.rs
@@ -16,6 +16,7 @@ use rhodium_core::{
 
 use crate::curve::{CurveExpr, CurveId, CurveRole, CurveSpec, Extrapolation};
 use crate::fx::{self, FxSpec, ResolvedFx};
+use crate::vol::{VolId, VolSpec};
 
 /// 曲線のパラメータブロックの識別子。
 #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
@@ -152,6 +153,8 @@ pub(crate) struct Curve {
     pub(crate) spec: CurveSpec,
     pub(crate) expr: Arc<ResolvedExpr>,
     pub(crate) blocks: Vec<usize>,
+    /// 予想曲線の予想を定める担保の割引曲線（予想の測度）の位置。分からなければ `None`。
+    pub(crate) forecast_measure: Option<usize>,
 }
 
 #[derive(Debug)]
@@ -210,6 +213,7 @@ pub struct MarketStructure {
     pub(crate) grids: Vec<BlockGrid>,
     pub(crate) fx: BTreeMap<CcyPair, ResolvedFx>,
     layout: ParamLayout,
+    volatilities: BTreeMap<VolId, VolSpec>,
 }
 
 impl MarketStructure {
@@ -276,6 +280,7 @@ impl MarketStructure {
                     spec,
                     expr,
                     blocks: blocks.into_iter().collect(),
+                    forecast_measure: None,
                 }
             })
             .collect();
@@ -285,11 +290,120 @@ impl MarketStructure {
             grids,
             layout,
             fx: BTreeMap::new(),
+            volatilities: BTreeMap::new(),
         };
         structure.fx = fx::resolve(&structure, fx)?;
         Ok(structure)
     }
 
+    /// 市場にボラ面と対応する θ ブロックを追加する。
+    ///
+    /// # Errors
+    ///
+    /// ID が重複するか、格子を評価時点に結び付けられない場合。
+    pub fn with_volatilities(mut self, specs: impl IntoIterator<Item = VolSpec>) -> Result<Self> {
+        let mut blocks = BTreeMap::new();
+        for spec in specs {
+            let id = spec.id().clone();
+            if self.volatilities.contains_key(&id)
+                || self
+                    .curves
+                    .iter()
+                    .any(|curve| curve.spec.id().as_str() == id.as_str())
+            {
+                return Err(Error::InvalidInput(format!(
+                    "duplicate market object {}",
+                    id.as_str()
+                )));
+            }
+            blocks.extend(spec.blocks(self.as_of)?);
+            self.volatilities.insert(id, spec);
+        }
+        for (block, count) in blocks {
+            let end = self
+                .layout
+                .len
+                .checked_add(count)
+                .ok_or_else(|| Error::InvalidInput("theta layout is too large".into()))?;
+            if self
+                .layout
+                .block_positions
+                .insert(block.clone(), self.layout.ids.len())
+                .is_some()
+            {
+                return Err(Error::InvalidInput(format!(
+                    "duplicate volatility block {block:?}"
+                )));
+            }
+            self.layout.ids.push(block);
+            self.layout.ranges.push(self.layout.len..end);
+            self.layout.len = end;
+            self.grids.push(BlockGrid {
+                dates: Box::new([]),
+                times: Box::new([]),
+            });
+        }
+        Ok(self)
+    }
+
+    /// 指定したボラ面の構造。
+    #[must_use]
+    pub fn volatility(&self, id: &VolId) -> Option<&VolSpec> {
+        self.volatilities.get(id)
+    }
+
+    /// 予想曲線ごとに、その予想を定める担保の割引曲線（予想の測度）を記録する。
+    ///
+    /// 予想曲線は指標ごとに一本で、その予想は校正商品の担保の下で定まる。
+    /// 同じ担保で割り引く支払いでは予想と割引の担保が一致するので、曲線の ID が違っても
+    /// `CrossCollateralProjection` の近似は要らない。記録の無い予想曲線は、
+    /// 自身を割引にも使う場合だけを同じ担保とみなす。
+    ///
+    /// # Errors
+    ///
+    /// 曲線が無いか、予想曲線でないか、測度が同じ通貨の割引曲線でなければ `InvalidInput`。
+    pub fn with_forecast_measures(
+        mut self,
+        measures: impl IntoIterator<Item = (CurveId, CurveId)>,
+    ) -> Result<Self> {
+        for (projection, discount) in measures {
+            let position = self.role_position(&projection, CurveRole::Projection)?;
+            let measure = self.role_position(&discount, CurveRole::Discount)?;
+            let [projection_ccy, discount_ccy] =
+                [position, measure].map(|index| self.curves[index].spec.currency());
+            if projection_ccy != discount_ccy {
+                return Err(Error::InvalidInput(format!(
+                    "forecast measure {discount} of {projection} has currency {discount_ccy}, expected {projection_ccy}"
+                )));
+            }
+            self.curves[position].forecast_measure = Some(measure);
+        }
+        Ok(self)
+    }
+
+    /// 予想曲線の予想を定める担保の割引曲線。記録が無ければ `None`。
+    #[must_use]
+    pub fn forecast_measure(&self, projection: &CurveId) -> Option<&CurveId> {
+        let position = self.layout.position(projection)?;
+        self.curves[position]
+            .forecast_measure
+            .map(|measure| self.curves[measure].spec.id())
+    }
+
+    fn role_position(&self, id: &CurveId, role: CurveRole) -> Result<usize> {
+        let position = self
+            .layout
+            .position(id)
+            .ok_or_else(|| Error::InvalidInput(format!("unknown curve {id}")))?;
+        if self.curves[position].spec.has_role(role) {
+            Ok(position)
+        } else {
+            Err(Error::InvalidInput(format!(
+                "curve {id} has no {role:?} role"
+            )))
+        }
+    }
+
     /// モデル時間の原点となる評価の瞬間。
     #[must_use]
     pub const fn as_of(&self) -> Instant {
diff --git a/crates/rhodium-market/src/view.rs b/crates/rhodium-market/src/view.rs
index 99cc3c140520efdf1fcdbd5f734d25630995e802..31316ee9356afda412b73ad393205d40d933a281 100644
--- a/crates/rhodium-market/src/view.rs
+++ b/crates/rhodium-market/src/view.rs
@@ -56,6 +56,22 @@ impl<K> ResolvedCurve<K> {
     }
 }
 
+impl ResolvedCurve<Projection> {
+    /// この予想が `discount` と同じ担保の下で定まるか。
+    ///
+    /// 構造に予想の測度（[`MarketStructure::with_forecast_measures`]）があればそれと比べ、
+    /// 無ければ同じ曲線を予想と割引に使う場合だけを同じ担保とみなす。
+    /// 違う構造の曲線どうしは比べられないので、常に違う担保として扱う。
+    #[must_use]
+    pub fn is_forecast_under(&self, discount: &ResolvedCurve<Discount>) -> bool {
+        Arc::ptr_eq(&self.structure, &discount.structure)
+            && self.structure.curves[self.position]
+                .forecast_measure
+                .unwrap_or(self.position)
+                == discount.position
+    }
+}
+
 impl MarketStructure {
     /// 割引曲線の位置を解決する。
     ///
@@ -334,6 +350,55 @@ impl<R: Real> Market<R> {
 }
 
 impl<'v, R: Real> MarketView<'v, R> {
+    /// 面の ID をこの市場ビューのハンドルに結ぶ。
+    ///
+    /// # Errors
+    ///
+    /// 面が無ければ `MissingData`。
+    pub fn volatility(&self, id: &crate::vol::VolId) -> Result<crate::vol::VolRef<'v, '_>> {
+        let spec = self.market.structure().volatility(id).ok_or_else(|| {
+            Error::MissingData(format!("unknown volatility surface {}", id.as_str()))
+        })?;
+        Ok(crate::vol::VolRef {
+            spec,
+            _brand: self.brand,
+        })
+    }
+
+    /// 準備済みの点を同じ面と市場ビューに結ぶ。
+    ///
+    /// # Errors
+    ///
+    /// 点が別の面や構造に属するか、必要なブロックが利用できない場合。
+    pub fn bind_vol<'a>(
+        &'a self,
+        face: crate::vol::VolRef<'v, 'a>,
+        point: &'a crate::vol::VolPoint,
+    ) -> Result<crate::vol::VolAt<'v, 'a>> {
+        if !Arc::ptr_eq(point.structure(), self.market.structure())
+            || !point.belongs_to(face.spec.id())
+        {
+            return Err(Error::InvalidInput(
+                "volatility point belongs to a different surface or structure".into(),
+            ));
+        }
+        self.market.require_blocks(
+            point
+                .blocks()
+                .iter()
+                .filter_map(|block| self.market.structure().layout().block_position(block)),
+        )?;
+        Ok(crate::vol::VolAt { face, point })
+    }
+
+    /// 結合済みの点から、絶対ストライクでのボラを読む。
+    ///
+    /// # Errors
+    ///
+    /// ストライクやボラを式の定義域で表現できない場合。
+    pub fn vol_at(&self, at: crate::vol::VolAt<'v, '_>, strike: R) -> Result<R> {
+        crate::vol::read(&self.market, at, strike)
+    }
     /// 同じ構造で準備した LGM の θ ブロックを読む。
     ///
     /// # Errors
diff --git a/crates/rhodium-market/src/vol.rs b/crates/rhodium-market/src/vol.rs
new file mode 100644
index 0000000000000000000000000000000000000000..eb8363b2a2bfd2ea2c6b48e1ae5b32aeb13e75fa
--- /dev/null
+++ b/crates/rhodium-market/src/vol.rs
@@ -0,0 +1,600 @@
+//! 金利ボラ面の格子、準備済みの読み出しと SABR の展開。
+
+mod sabr;
+
+use std::{collections::BTreeSet, str::FromStr, sync::Arc};
+
+use rhodium_core::{
+    Error, Result,
+    brand::Brand,
+    id::IndexId,
+    quantity::{Rate, Simple, VolKind},
+    time::{Instant, model_time},
+};
+use scioxide_ad::{Mask, Real};
+
+use crate::{
+    curve::CurveId,
+    market::Market,
+    structure::{BlockId, MarketStructure},
+};
+
+pub use sabr::Sabr;
+
+/// ボラ面を参照する ID。
+#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Deserialize)]
+#[serde(try_from = "String")]
+pub struct VolId(CurveId);
+
+impl FromStr for VolId {
+    type Err = Error;
+
+    fn from_str(value: &str) -> Result<Self> {
+        value.parse().map(Self)
+    }
+}
+
+impl TryFrom<String> for VolId {
+    type Error = Error;
+
+    fn try_from(value: String) -> Result<Self> {
+        value.parse()
+    }
+}
+
+impl VolId {
+    /// ID の文字列。
+    #[must_use]
+    pub fn as_str(&self) -> &str {
+        self.0.as_str()
+    }
+
+    /// ボラの θ ブロック。キューブの部品はノードごとに区別する。
+    ///
+    /// # Errors
+    ///
+    /// 部品名が不正な場合。
+    pub fn block(&self, part: &str) -> Result<BlockId> {
+        BlockId::new(self.0.clone(), part)
+    }
+}
+
+/// スワップションの満期とテナーの間の補間。
+#[derive(Debug, Clone, Copy, Default, serde::Deserialize)]
+#[serde(rename_all = "snake_case")]
+pub enum CubeInterpolation {
+    /// 同じ ATM オフセットの総分散を補間する。
+    #[default]
+    TotalVariance,
+    /// α、ρ、ν とノードのフォワードを補間する。
+    Parameters,
+}
+
+/// SABR ノードで校正するパラメータの範囲。
+#[derive(Debug, Clone, Copy, serde::Deserialize)]
+#[serde(rename_all = "snake_case", deny_unknown_fields)]
+pub enum SabrFit {
+    /// α を ATM に当て、ρ と ν を定義に固定する。
+    AtmOnly {
+        /// 固定する ρ。
+        rho: f64,
+        /// 固定する ν。
+        nu: f64,
+    },
+    /// すべての行を重み付き最小二乗で当てる。
+    LeastSquares,
+    /// ATM を厳密に当て、ウィングを重み付き最小二乗で当てる。
+    ExactAtm,
+}
+
+impl SabrFit {
+    /// 固定する ρ と ν が Hagan の展開の定義域（|ρ| < 1、ν ≥ 0、有限）にあるか確かめる。
+    ///
+    /// 定義の読み込みと面の構築の両方がこの検査を使い、条件を一か所に保つ。
+    ///
+    /// # Errors
+    ///
+    /// `AtmOnly` の ρ か ν が定義域の外なら `InvalidInput`。
+    pub fn check(self) -> Result<Self> {
+        match self {
+            Self::AtmOnly { rho, nu }
+                if !(rho.is_finite() && rho.abs() < 1.0 && nu.is_finite() && nu >= 0.0) =>
+            {
+                Err(Error::InvalidInput(format!(
+                    "fixed SABR parameters outside the domain: rho {rho}, nu {nu}"
+                )))
+            }
+            Self::AtmOnly { .. } | Self::LeastSquares | Self::ExactAtm => Ok(self),
+        }
+    }
+
+    /// 校正する SABR 座標の個数。末尾には別にフォワードを置く。
+    #[must_use]
+    pub const fn dimension(self) -> usize {
+        match self {
+            Self::AtmOnly { .. } => 1,
+            Self::LeastSquares | Self::ExactAtm => 3,
+        }
+    }
+
+    /// 校正座標を α、ρ、ν に戻す。
+    #[must_use]
+    pub fn parameters<R: Real>(self, [log_alpha, rho, log_nu]: [R; 3]) -> [R; 3] {
+        let alpha = log_alpha.exp();
+        match self {
+            Self::AtmOnly { rho, nu } => [alpha, R::constant(rho), R::constant(nu)],
+            Self::LeastSquares | Self::ExactAtm => [alpha, rho.tanh(), log_nu.exp()],
+        }
+    }
+}
+
+/// 評価時点で満期とテナーを解決した金利ボラ面の構造。
+#[derive(Debug, Clone)]
+pub enum VolSpec {
+    /// 満期 × テナーの SABR キューブ。
+    SwaptionCube(SwaptionCube),
+    /// 満期 × ストライクのキャップレット面。
+    Caplet(CapletSurface),
+}
+
+/// 評価時点で解決した SABR キューブの座標。
+#[derive(Debug, Clone)]
+pub struct SwaptionCube {
+    id: VolId,
+    index: IndexId,
+    sabr: Sabr,
+    fit: SabrFit,
+    interpolation: CubeInterpolation,
+    expiries: Vec<Instant>,
+    tenors: Vec<f64>,
+}
+
+/// 一列分のキャップレット節点。節点の満期まで同じボラを使う。
+#[derive(Debug, Clone)]
+pub struct CapletColumn {
+    strike: Rate<Simple>,
+    ends: Vec<Instant>,
+}
+
+impl CapletColumn {
+    /// 絶対ストライクと昇順の節点満期を結ぶ。
+    ///
+    /// # Errors
+    ///
+    /// 満期が空か、昇順でない場合。
+    pub fn new(strike: Rate<Simple>, ends: Vec<Instant>) -> Result<Self> {
+        check_ends(&ends)?;
+        Ok(Self { strike, ends })
+    }
+}
+
+/// キャップレット面の構造。ATM のストライクは各区間の θ に置く。
+#[derive(Debug, Clone)]
+pub struct CapletSurface {
+    id: VolId,
+    index: IndexId,
+    kind: VolKind,
+    absolute: Vec<CapletColumn>,
+    atm: Vec<Instant>,
+}
+
+fn check_ends(ends: &[Instant]) -> Result<()> {
+    if ends.is_empty() || !ends.is_sorted_by(|a, b| a < b) {
+        return Err(Error::InvalidInput(
+            "volatility expiries must be nonempty and strictly increasing".into(),
+        ));
+    }
+    Ok(())
+}
+
+impl VolSpec {
+    /// SABR キューブを、昇順の満期とテナー（年）から作る。
+    ///
+    /// # Errors
+    ///
+    /// 軸または固定パラメータが不正な場合。
+    pub fn swaption_cube(
+        id: VolId,
+        index: IndexId,
+        sabr: Sabr,
+        fit: SabrFit,
+        interpolation: CubeInterpolation,
+        expiries: Vec<Instant>,
+        tenors: Vec<f64>,
+    ) -> Result<Self> {
+        check_ends(&expiries)?;
+        let fit = fit.check()?;
+        if tenors.is_empty()
+            || tenors.iter().any(|t| !t.is_finite() || *t <= 0.0)
+            || !tenors.is_sorted_by(|a, b| a < b)
+        {
+            return Err(Error::InvalidInput("invalid cube tenor axis".into()));
+        }
+        Ok(Self::SwaptionCube(SwaptionCube {
+            id,
+            index,
+            sabr,
+            fit,
+            interpolation,
+            expiries,
+            tenors,
+        }))
+    }
+
+    /// キャップレット面を、ストライクの昇順の列と任意の ATM 満期列から作る。
+    ///
+    /// # Errors
+    ///
+    /// 列が無いか、格子またはボラのシフトが不正な場合。
+    pub fn caplet(
+        id: VolId,
+        index: IndexId,
+        kind: VolKind,
+        absolute: Vec<CapletColumn>,
+        atm: Vec<Instant>,
+    ) -> Result<Self> {
+        if !atm.is_empty() {
+            check_ends(&atm)?;
+        }
+        if (absolute.is_empty() && atm.is_empty())
+            || !absolute.is_sorted_by(|a, b| a.strike.value() < b.strike.value())
+            || matches!(kind, VolKind::ShiftedLognormal { shift } if !shift.is_finite())
+        {
+            return Err(Error::InvalidInput(
+                "invalid caplet strike grid or volatility kind".into(),
+            ));
+        }
+        Ok(Self::Caplet(CapletSurface {
+            id,
+            index,
+            kind,
+            absolute,
+            atm,
+        }))
+    }
+
+    /// 面の ID。
+    #[must_use]
+    pub const fn id(&self) -> &VolId {
+        match self {
+            Self::SwaptionCube(spec) => &spec.id,
+            Self::Caplet(spec) => &spec.id,
+        }
+    }
+
+    /// 面が対象とする変動金利の指標。
+    #[must_use]
+    pub const fn index(&self) -> &IndexId {
+        match self {
+            Self::SwaptionCube(spec) => &spec.index,
+            Self::Caplet(spec) => &spec.index,
+        }
+    }
+
+    /// 面のボラの種類。
+    #[must_use]
+    pub const fn kind(&self) -> VolKind {
+        match self {
+            Self::SwaptionCube(spec) => spec.sabr.kind(),
+            Self::Caplet(spec) => spec.kind,
+        }
+    }
+
+    /// 満期、テナーの辞書順での SABR ノードのブロック。
+    ///
+    /// # Errors
+    ///
+    /// ブロック ID を作れない場合。
+    pub fn swaption_block(id: &VolId, expiry: usize, tenor: usize) -> Result<BlockId> {
+        id.block(&format!("swaption.{expiry}.{tenor}"))
+    }
+
+    /// キャップレット面全体のブロック。
+    ///
+    /// # Errors
+    ///
+    /// ブロック ID を作れない場合。
+    pub fn caplet_block(id: &VolId) -> Result<BlockId> {
+        id.block("caplet")
+    }
+
+    pub(crate) fn blocks(&self, as_of: Instant) -> Result<Vec<(BlockId, usize)>> {
+        match self {
+            Self::SwaptionCube(cube) => {
+                if cube.expiries[0] <= as_of {
+                    return Err(Error::InvalidInput(
+                        "swaption cube expiry must follow as-of".into(),
+                    ));
+                }
+                (0..cube.expiries.len())
+                    .flat_map(|e| (0..cube.tenors.len()).map(move |t| (e, t)))
+                    .map(|(e, t)| {
+                        Ok((
+                            Self::swaption_block(&cube.id, e, t)?,
+                            cube.fit.dimension() + 1,
+                        ))
+                    })
+                    .collect()
+            }
+            Self::Caplet(caplet) => {
+                let count = caplet.absolute.iter().map(|c| c.ends.len()).sum::<usize>()
+                    + 2 * caplet.atm.len();
+                Ok(vec![(Self::caplet_block(&caplet.id)?, count)])
+            }
+        }
+    }
+}
+
+/// 市場構造で補間位置を解決した読み出し点。
+#[derive(Debug, Clone)]
+pub struct VolPoint {
+    structure: Arc<MarketStructure>,
+    id: VolId,
+    blocks: Vec<BlockId>,
+    point: Point,
+}
+
+#[derive(Debug, Clone)]
+enum Point {
+    Cube {
+        time: f64,
+        nodes: Vec<(usize, f64, f64)>,
+    },
+    Caplet {
+        block: usize,
+        absolute: Vec<(f64, usize)>,
+        atm: Option<(usize, usize)>,
+    },
+}
+
+impl VolPoint {
+    /// この点で読み出す θ ブロック。
+    #[must_use]
+    pub fn blocks(&self) -> &[BlockId] {
+        &self.blocks
+    }
+
+    /// 点を読み出す市場構造。
+    #[must_use]
+    pub const fn structure(&self) -> &Arc<MarketStructure> {
+        &self.structure
+    }
+
+    pub(crate) fn belongs_to(&self, id: &VolId) -> bool {
+        self.id == *id
+    }
+}
+
+/// 市場ビューで解決した面のハンドル。
+#[derive(Debug, Clone, Copy)]
+pub struct VolRef<'v, 'a> {
+    pub(crate) spec: &'a VolSpec,
+    pub(crate) _brand: Brand<'v>,
+}
+
+impl VolRef<'_, '_> {
+    /// 面が対象とする変動金利の指標。
+    #[must_use]
+    pub const fn index(&self) -> &IndexId {
+        self.spec.index()
+    }
+
+    /// ボラの種類。
+    #[must_use]
+    pub const fn kind(self) -> VolKind {
+        self.spec.kind()
+    }
+}
+
+/// 市場ビュー内で面と点を結んだボラの読み出し。
+#[derive(Debug, Clone, Copy)]
+pub struct VolAt<'v, 'a> {
+    pub(crate) face: VolRef<'v, 'a>,
+    pub(crate) point: &'a VolPoint,
+}
+
+impl MarketStructure {
+    /// 満期と、キューブの場合はテナー（年）をボラ面の読み出し点にする。
+    ///
+    /// # Errors
+    ///
+    /// 面、満期、テナーの組が不正な場合。
+    pub fn prepare_vol_point(
+        self: &Arc<Self>,
+        id: &VolId,
+        expiry: Instant,
+        tenor: Option<f64>,
+    ) -> Result<VolPoint> {
+        let spec = self.volatility(id).ok_or_else(|| {
+            Error::InvalidInput(format!("unknown volatility surface {}", id.as_str()))
+        })?;
+        let block_position = |block: &BlockId| {
+            self.layout()
+                .block_position(block)
+                .ok_or_else(|| Error::Internal("volatility block missing from layout".into()))
+        };
+        let (point, blocks) = match (spec, tenor) {
+            (VolSpec::SwaptionCube(cube), Some(tenor))
+                if tenor.is_finite() && tenor > 0.0 && expiry > self.as_of() =>
+            {
+                let time = model_time(self.as_of(), expiry).years();
+                let times: Vec<_> = cube
+                    .expiries
+                    .iter()
+                    .map(|&end| model_time(self.as_of(), end).years())
+                    .collect();
+                let mut nodes = Vec::new();
+                let mut blocks = BTreeSet::new();
+                for (e, ew) in weights(&times, time) {
+                    for (t, tw) in weights(&cube.tenors, tenor) {
+                        let block = VolSpec::swaption_block(id, e, t)?;
+                        nodes.push((block_position(&block)?, times[e], ew * tw));
+                        blocks.insert(block);
+                    }
+                }
+                (
+                    Point::Cube {
+                        time: time.clamp(times[0], times[times.len() - 1]),
+                        nodes,
+                    },
+                    blocks.into_iter().collect(),
+                )
+            }
+            (VolSpec::Caplet(face), None) => {
+                let block = VolSpec::caplet_block(id)?;
+                let mut start = 0;
+                let absolute = face
+                    .absolute
+                    .iter()
+                    .map(|column| {
+                        let bucket = bucket(&column.ends, expiry);
+                        let element = start + bucket;
+                        start += column.ends.len();
+                        (column.strike.value(), element)
+                    })
+                    .collect();
+                let atm = (!face.atm.is_empty()).then(|| {
+                    let vol = start + bucket(&face.atm, expiry);
+                    (vol, vol + face.atm.len())
+                });
+                (
+                    Point::Caplet {
+                        block: block_position(&block)?,
+                        absolute,
+                        atm,
+                    },
+                    vec![block],
+                )
+            }
+            _ => {
+                return Err(Error::InvalidInput(format!(
+                    "invalid coordinates for volatility surface {}",
+                    id.as_str()
+                )));
+            }
+        };
+        Ok(VolPoint {
+            structure: Arc::clone(self),
+            id: id.clone(),
+            blocks,
+            point,
+        })
+    }
+}
+
+fn weights(axis: &[f64], value: f64) -> Vec<(usize, f64)> {
+    let upper = axis.partition_point(|&x| x < value);
+    if upper == 0 || upper == axis.len() || axis[upper] <= value {
+        return vec![(upper.min(axis.len() - 1), 1.0)];
+    }
+    let lower = upper - 1;
+    let weight = (value - axis[lower]) / (axis[upper] - axis[lower]);
+    vec![(lower, 1.0 - weight), (upper, weight)]
+}
+
+fn bucket(ends: &[Instant], expiry: Instant) -> usize {
+    ends.partition_point(|&end| end < expiry)
+        .min(ends.len() - 1)
+}
+
+pub(crate) fn read<R: Real>(market: &Market<R>, at: VolAt<'_, '_>, strike: R) -> Result<R> {
+    match (at.face.spec, &at.point.point) {
+        (VolSpec::SwaptionCube(cube), Point::Cube { time, nodes }) => {
+            let dimension = cube.fit.dimension();
+            let forward = nodes.iter().fold(R::zero(), |f, &(block, _, weight)| {
+                f + market.theta_block(block).at(dimension) * weight
+            });
+            match cube.interpolation {
+                CubeInterpolation::Parameters => {
+                    let mut parameters = [R::zero(); 3];
+                    for &(block, _, weight) in nodes {
+                        for (sum, value) in parameters
+                            .iter_mut()
+                            .zip(sabr_parameters(market, cube.fit, block))
+                        {
+                            *sum += value * weight;
+                        }
+                    }
+                    cube.sabr.volatility(forward, strike, *time, parameters)
+                }
+                CubeInterpolation::TotalVariance => {
+                    let offset = strike - forward;
+                    let mut variance = R::zero();
+                    for &(block, expiry, weight) in nodes {
+                        let node_forward = market.theta_block(block).at(dimension);
+                        let vol = cube.sabr.volatility(
+                            node_forward,
+                            node_forward + offset,
+                            expiry,
+                            sabr_parameters(market, cube.fit, block),
+                        )?;
+                        variance += vol.powi(2) * (expiry * weight);
+                    }
+                    Ok((variance / *time).sqrt())
+                }
+            }
+        }
+        (
+            VolSpec::Caplet(_),
+            Point::Caplet {
+                block,
+                absolute,
+                atm,
+            },
+        ) => {
+            if !strike.is_finite().all() {
+                return Err(Error::InvalidInput("nonfinite caplet strike".into()));
+            }
+            let theta = market.theta_block(*block);
+            let points = absolute
+                .iter()
+                .map(|&(k, vol)| (R::constant(k), theta.at(vol)))
+                .chain(atm.map(|(vol, k)| (theta.at(k), theta.at(vol))));
+            let value = strike_interpolation(points, strike);
+            if !(value.is_finite().all() && value.ge(R::zero()).all()) {
+                return Err(Error::Numerical("invalid caplet volatility".into()));
+            }
+            Ok(value)
+        }
+        _ => Err(Error::Internal("invalid bound volatility point".into())),
+    }
+}
+
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#201 で各 variant を網羅する"
+)]
+fn sabr_parameters<R: Real>(market: &Market<R>, fit: SabrFit, block: usize) -> [R; 3] {
+    let theta = market.theta_block(block);
+    match fit {
+        SabrFit::AtmOnly { .. } => fit.parameters([theta.at(0), R::zero(), R::zero()]),
+        _ => fit.parameters([theta.at(0), theta.at(1), theta.at(2)]),
+    }
+}
+
+// ATM のストライクは θ に依存するので、各レーンで近い左右の点を選ぶ。
+fn strike_interpolation<R: Real>(mut points: impl Iterator<Item = (R, R)>, strike: R) -> R {
+    let Some((first, vol)) = points.next() else {
+        return R::zero();
+    };
+    let mut lower = (first, vol);
+    let mut upper = (first, vol);
+    for (k, value) in points {
+        let nearer_lower = (k.le(strike) & (lower.0.gt(strike) | k.gt(lower.0)))
+            | (lower.0.gt(strike) & k.lt(lower.0));
+        let nearer_upper = (k.ge(strike) & (upper.0.lt(strike) | k.lt(upper.0)))
+            | (upper.0.lt(strike) & k.gt(upper.0));
+        lower = (
+            R::select(nearer_lower, k, lower.0),
+            R::select(nearer_lower, value, lower.1),
+        );
+        upper = (
+            R::select(nearer_upper, k, upper.0),
+            R::select(nearer_upper, value, upper.1),
+        );
+    }
+    let width = upper.0 - lower.0;
+    let weight = ((strike - lower.0) / R::select(width.gt(R::zero()), width, R::one()))
+        .clamp(R::zero(), R::one());
+    lower.1 + (upper.1 - lower.1) * weight
+}
diff --git a/crates/rhodium-market/src/vol/sabr.rs b/crates/rhodium-market/src/vol/sabr.rs
new file mode 100644
index 0000000000000000000000000000000000000000..b33788deecb9e8cce7af06ee544320c20087136e
--- /dev/null
+++ b/crates/rhodium-market/src/vol/sabr.rs
@@ -0,0 +1,129 @@
+//! Hagan の展開と、SABR の校正座標。
+
+use rhodium_core::{Error, Result, quantity::VolKind};
+use scioxide_ad::{Mask, Real};
+
+/// SABR のうち定義に固定する値。
+#[derive(Debug, Clone, Copy)]
+pub struct Sabr {
+    beta: f64,
+    shift: f64,
+    kind: VolKind,
+}
+
+impl Sabr {
+    /// Hagan の展開をボラの種類に対応させる。正規 SABR でもシフトを指定できる。
+    ///
+    /// # Errors
+    ///
+    /// β が [0, 1] の外、またはシフトの指定が種類と一致しなければ `InvalidInput`。
+    pub fn new(kind: VolKind, beta: f64, shift: f64) -> Result<Self> {
+        if !(0.0..=1.0).contains(&beta)
+            || !shift.is_finite()
+            || match kind {
+                VolKind::Normal => false,
+                VolKind::Lognormal => shift != 0.0,
+                VolKind::ShiftedLognormal { shift: quoted } => {
+                    (shift - quoted).abs() > 0.0 || !quoted.is_finite()
+                }
+            }
+        {
+            return Err(Error::InvalidInput("invalid SABR beta or shift".into()));
+        }
+        Ok(Self { beta, shift, kind })
+    }
+
+    /// ボラの種類。
+    #[must_use]
+    pub const fn kind(self) -> VolKind {
+        self.kind
+    }
+
+    /// 固定した β。
+    #[must_use]
+    pub const fn beta(self) -> f64 {
+        self.beta
+    }
+
+    /// 原資産とストライクに加えるシフト。
+    #[must_use]
+    pub const fn shift(self) -> f64 {
+        self.shift
+    }
+
+    /// α、ρ、ν に対する Hagan のインプライドボラ。
+    ///
+    /// # Errors
+    ///
+    /// パラメータ、時間、シフトした金利が定義域の外か、展開が正の有限なボラを返さない場合。
+    pub fn volatility<R: Real>(
+        self,
+        forward: R,
+        strike: R,
+        time: f64,
+        [alpha, rho, nu]: [R; 3],
+    ) -> Result<R> {
+        let f = forward + self.shift;
+        let k = strike + self.shift;
+        let beta = self.beta;
+        if !(time.is_finite()
+            && time >= 0.0
+            && f.gt(R::zero()).all()
+            && k.gt(R::zero()).all()
+            && alpha.gt(R::zero()).all()
+            && rho.abs().lt(R::one()).all()
+            && nu.ge(R::zero()).all())
+        {
+            return Err(Error::Numerical(
+                "SABR inputs outside the expansion domain".into(),
+            ));
+        }
+        let log_m = ((f - k) / k).ln_1p();
+        let log_fk = f.ln() + k.ln();
+        let one_minus_beta = 1.0 - beta;
+        let a = (log_fk * (one_minus_beta * 0.5)).exp();
+        let z = nu * a * log_m / alpha;
+        let small = z.abs().lt(R::constant(1e-5));
+        // sqrt(1 - 2ρz + z²) - 1 を有理化し、ATM 付近で log の桁落ちを避ける。
+        let root = ((z - rho).powi(2) + (R::one() - rho.powi(2))).sqrt();
+        let rational = (z - rho * 2.0) / (root + 1.0);
+        let safe_z = R::select(small, R::one(), z);
+        let log_arg = (rational + 1.0) * z / (R::one() - rho);
+        // 左ウィングで log1p の引数が -1 に丸まる場合は共役で同値な形を使う。
+        let left = z.lt(R::constant(-0.5));
+        let left_value = (R::select(left, root - z + rho, R::one()) / (rho + 1.0)).ln() * -1.0;
+        let direct = R::select(
+            left,
+            left_value,
+            R::select(small | left, R::one(), log_arg).ln_1p(),
+        );
+        let ratio = R::select(
+            small,
+            R::one() - rho * z * 0.5 + (R::constant(2.0) - rho.powi(2) * 3.0) * z.powi(2) / 12.0,
+            safe_z / R::select(small, R::one(), direct),
+        );
+        let alpha_a = alpha / a;
+        let polynomial = |square: R| R::one() + square / 24.0 + square.powi(2) / 1920.0;
+        let denominator = polynomial((log_m * one_minus_beta).powi(2));
+        let (leading, correction) = match self.kind {
+            VolKind::Normal => (
+                alpha * (log_fk * (beta * 0.5)).exp() * polynomial(log_m.powi(2)) / denominator,
+                beta * (beta - 2.0),
+            ),
+            VolKind::Lognormal | VolKind::ShiftedLognormal { .. } => {
+                (alpha_a / denominator, one_minus_beta * one_minus_beta)
+            }
+        };
+        let drift = alpha_a.powi(2) * (correction / 24.0)
+            + rho * nu * alpha_a * (beta / 4.0)
+            + (R::constant(2.0) - rho.powi(2) * 3.0) * nu.powi(2) / 24.0;
+        let vol = leading * ratio * (R::one() + drift * time);
+        if vol.is_finite().all() && vol.gt(R::zero()).all() {
+            Ok(vol)
+        } else {
+            Err(Error::Numerical(
+                "SABR expansion did not produce a positive volatility".into(),
+            ))
+        }
+    }
+}
diff --git a/crates/rhodium-market/tests/compile_fail.rs b/crates/rhodium-market/tests/compile_fail.rs
index 63d987e47a3ecdea72558cbc451d893d7ec5cb59..503a40fc37d1c07f293d965758fdd89c05fe4939 100644
--- a/crates/rhodium-market/tests/compile_fail.rs
+++ b/crates/rhodium-market/tests/compile_fail.rs
@@ -67,3 +67,9 @@ fn fx_ref_cannot_escape_view() {
     // FX ハンドルは結合先のビューの外へ持ち出せない。
     compile_fail_pair!("fx_ref_cannot_escape_view");
 }
+
+#[test]
+fn volatility_point_cross_view() {
+    // 同じ面と点でも、結合したボラの値はその市場ビューの内側でだけ読める。
+    compile_fail_pair!("volatility_point_cross_view");
+}
diff --git a/crates/rhodium-market/tests/ui/amount_cross_scope.json b/crates/rhodium-market/tests/ui/amount_cross_scope.json
new file mode 100644
index 0000000000000000000000000000000000000000..83e5f092f38968cd8f39fda057cc2366018fd367
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/amount_cross_scope.json
@@ -0,0 +1,29 @@
+[
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/amount_cross_scope.rs",
+        "line": 13,
+        "column": 34
+      }
+    ]
+  },
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/amount_cross_scope.rs",
+        "line": 13,
+        "column": 34
+      },
+      {
+        "file": "tests/ui/amount_cross_scope.rs",
+        "line": 13,
+        "column": 34
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/amount_cross_scope.stderr b/crates/rhodium-market/tests/ui/amount_cross_scope.stderr
deleted file mode 100644
index 354e63240632e761bd2da7a0937a6bff5b18146b..0000000000000000000000000000000000000000
--- a/crates/rhodium-market/tests/ui/amount_cross_scope.stderr
+++ /dev/null
@@ -1,28 +0,0 @@
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/amount_cross_scope.rs:13:34
-   |
-11 |             let amount = first.present(5.0);
-   |                 ------ `amount` declared here, outside of the closure body
-12 |             view.ccy_scope(Currency::JPY, curve, |second| {
-   |                                                   ------ `second` is only valid in the closure body
-13 |                 let _ = amount + second.present(1.0);
-   |                                  ^^^^^^^^^^^^^^^^^^^ `second` escapes the closure body here
-   |
-   = note: requirement occurs because of the type `CcyScope<'_, '_>`, which makes the generic argument `'_` invariant
-   = note: the struct `CcyScope<'c, 'v, R>` is invariant over the parameter `'c`
-   = help: see <https://doc.rust-lang.org/nomicon/subtyping.html> for more information about variance
-
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/amount_cross_scope.rs:13:34
-   |
-10 |         view.ccy_scope(Currency::JPY, curve, |first| {
-   |                                               -----
-   |                                               |
-   |                                               `first` is only valid in the closure body
-   |                                               has type `CcyScope<'1, '_>`
-...
-13 |                 let _ = amount + second.present(1.0);
-   |                                  ^^^^^^^^^^^^^^^^^^^
-   |                                  |
-   |                                  `first` escapes the closure body here
-   |                                  argument requires that `'1` must outlive `'static`
diff --git a/crates/rhodium-market/tests/ui/amount_plus_undiscounted.json b/crates/rhodium-market/tests/ui/amount_plus_undiscounted.json
new file mode 100644
index 0000000000000000000000000000000000000000..99f08c43cc4bfd2f00daa34a20e0324804a83da7
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/amount_plus_undiscounted.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0308",
+    "message": "mismatched types",
+    "spans": [
+      {
+        "file": "tests/ui/amount_plus_undiscounted.rs",
+        "line": 14,
+        "column": 30
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/amount_plus_undiscounted.stderr b/crates/rhodium-market/tests/ui/amount_plus_undiscounted.stderr
deleted file mode 100644
index fa9fc381fc2517a03723f342566b2a3d0d26755b..0000000000000000000000000000000000000000
--- a/crates/rhodium-market/tests/ui/amount_plus_undiscounted.stderr
+++ /dev/null
@@ -1,8 +0,0 @@
-error[E0308]: mismatched types
-  --> tests/ui/amount_plus_undiscounted.rs:14:30
-   |
-14 |             let _ = coupon + 100.0;
-   |                              ^^^^^ expected `Amount<'_>`, found floating-point number
-   |
-   = note: expected struct `Amount<'_>`
-                found type `{float}`
diff --git a/crates/rhodium-market/tests/ui/curve_ref_cross_view.json b/crates/rhodium-market/tests/ui/curve_ref_cross_view.json
new file mode 100644
index 0000000000000000000000000000000000000000..3856b4f707558408070b6238eac225376700ffe6
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/curve_ref_cross_view.json
@@ -0,0 +1,29 @@
+[
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/curve_ref_cross_view.rs",
+        "line": 10,
+        "column": 35
+      }
+    ]
+  },
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/curve_ref_cross_view.rs",
+        "line": 10,
+        "column": 35
+      },
+      {
+        "file": "tests/ui/curve_ref_cross_view.rs",
+        "line": 10,
+        "column": 35
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/curve_ref_cross_view.stderr b/crates/rhodium-market/tests/ui/curve_ref_cross_view.stderr
deleted file mode 100644
index 084162732d89b3197e3b87a8307b440fc7a1ba69..0000000000000000000000000000000000000000
--- a/crates/rhodium-market/tests/ui/curve_ref_cross_view.stderr
+++ /dev/null
@@ -1,28 +0,0 @@
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/curve_ref_cross_view.rs:10:35
-   |
- 9 |         let discount = first.discount_curve(&id)?;
-   |             -------- `discount` declared here, outside of the closure body
-10 |         market.with_view(|second| second.ccy_scope(Currency::JPY, discount, |_| Ok(())))
-   |                           ------  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ `second` escapes the closure body here
-   |                           |
-   |                           `second` is only valid in the closure body
-   |
-   = note: requirement occurs because of the type `CurveRef<'_, rhodium_market::view::Discount>`, which makes the generic argument `'_` invariant
-   = note: the struct `CurveRef<'v, K>` is invariant over the parameter `'v`
-   = help: see <https://doc.rust-lang.org/nomicon/subtyping.html> for more information about variance
-
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/curve_ref_cross_view.rs:10:35
-   |
- 8 |     market.with_view(|first| {
-   |                       -----
-   |                       |
-   |                       `first` is only valid in the closure body
-   |                       has type `MarketView<'1>`
- 9 |         let discount = first.discount_curve(&id)?;
-10 |         market.with_view(|second| second.ccy_scope(Currency::JPY, discount, |_| Ok(())))
-   |                                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
-   |                                   |
-   |                                   `first` escapes the closure body here
-   |                                   argument requires that `'1` must outlive `'static`
diff --git a/crates/rhodium-market/tests/ui/discount_point_cross_view.json b/crates/rhodium-market/tests/ui/discount_point_cross_view.json
new file mode 100644
index 0000000000000000000000000000000000000000..de21be2ec78b13eaf11646eebb9db268d9f3d1f6
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/discount_point_cross_view.json
@@ -0,0 +1,29 @@
+[
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/discount_point_cross_view.rs",
+        "line": 11,
+        "column": 43
+      }
+    ]
+  },
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/discount_point_cross_view.rs",
+        "line": 11,
+        "column": 43
+      },
+      {
+        "file": "tests/ui/discount_point_cross_view.rs",
+        "line": 11,
+        "column": 43
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/discount_point_cross_view.stderr b/crates/rhodium-market/tests/ui/discount_point_cross_view.stderr
deleted file mode 100644
index 2972bd638e8b53cb1fc3c30261f931294e5dd97b..0000000000000000000000000000000000000000
--- a/crates/rhodium-market/tests/ui/discount_point_cross_view.stderr
+++ /dev/null
@@ -1,28 +0,0 @@
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/discount_point_cross_view.rs:11:43
-   |
-10 |         let at = bound.discounts()[0].points()[0];
-   |             -- `at` declared here, outside of the closure body
-11 |         let _ = market.with_view(|second| second.df(at));
-   |                                   ------  ^^^^^^^^^^^^^ `second` escapes the closure body here
-   |                                   |
-   |                                   `second` is only valid in the closure body
-   |
-   = note: requirement occurs because of the type `CurveAt<'_, rhodium_market::view::Discount>`, which makes the generic argument `'_` invariant
-   = note: the struct `CurveAt<'v, K, R>` is invariant over the parameter `'v`
-   = help: see <https://doc.rust-lang.org/nomicon/subtyping.html> for more information about variance
-
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/discount_point_cross_view.rs:11:43
-   |
- 8 |     market.with_view(|first| {
-   |                       -----
-   |                       |
-   |                       `first` is only valid in the closure body
-   |                       has type `MarketView<'1>`
-...
-11 |         let _ = market.with_view(|second| second.df(at));
-   |                                           ^^^^^^^^^^^^^
-   |                                           |
-   |                                           `first` escapes the closure body here
-   |                                           argument requires that `'1` must outlive `'static`
diff --git a/crates/rhodium-market/tests/ui/forward_point_cross_view.json b/crates/rhodium-market/tests/ui/forward_point_cross_view.json
new file mode 100644
index 0000000000000000000000000000000000000000..f6a4b43148ac07c9413cbafdf59258676ab12984
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/forward_point_cross_view.json
@@ -0,0 +1,29 @@
+[
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/forward_point_cross_view.rs",
+        "line": 11,
+        "column": 43
+      }
+    ]
+  },
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/forward_point_cross_view.rs",
+        "line": 11,
+        "column": 43
+      },
+      {
+        "file": "tests/ui/forward_point_cross_view.rs",
+        "line": 11,
+        "column": 43
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/forward_point_cross_view.stderr b/crates/rhodium-market/tests/ui/forward_point_cross_view.stderr
deleted file mode 100644
index f66b1ba5aff9768cf64d63134e22a1db40ea52f1..0000000000000000000000000000000000000000
--- a/crates/rhodium-market/tests/ui/forward_point_cross_view.stderr
+++ /dev/null
@@ -1,28 +0,0 @@
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/forward_point_cross_view.rs:11:43
-   |
-10 |         let at = bound.projections()[0].points()[0];
-   |             -- `at` declared here, outside of the closure body
-11 |         let _ = market.with_view(|second| second.forward(at));
-   |                                   ------  ^^^^^^^^^^^^^^^^^^ `second` escapes the closure body here
-   |                                   |
-   |                                   `second` is only valid in the closure body
-   |
-   = note: requirement occurs because of the type `ForwardAt<'_>`, which makes the generic argument `'_` invariant
-   = note: the struct `ForwardAt<'v, R>` is invariant over the parameter `'v`
-   = help: see <https://doc.rust-lang.org/nomicon/subtyping.html> for more information about variance
-
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/forward_point_cross_view.rs:11:43
-   |
- 8 |     market.with_view(|first| {
-   |                       -----
-   |                       |
-   |                       `first` is only valid in the closure body
-   |                       has type `MarketView<'1>`
-...
-11 |         let _ = market.with_view(|second| second.forward(at));
-   |                                           ^^^^^^^^^^^^^^^^^^
-   |                                           |
-   |                                           `first` escapes the closure body here
-   |                                           argument requires that `'1` must outlive `'static`
diff --git a/crates/rhodium-market/tests/ui/fx_ref_cannot_escape_view.json b/crates/rhodium-market/tests/ui/fx_ref_cannot_escape_view.json
new file mode 100644
index 0000000000000000000000000000000000000000..afaef58c8c65062bccfc0a109f9977fdd16525c0
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/fx_ref_cannot_escape_view.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": null,
+    "message": "lifetime may not live long enough",
+    "spans": [
+      {
+        "file": "tests/ui/fx_ref_cannot_escape_view.rs",
+        "line": 5,
+        "column": 38
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/fx_ref_cannot_escape_view.stderr b/crates/rhodium-market/tests/ui/fx_ref_cannot_escape_view.stderr
deleted file mode 100644
index f21829f1d034e3a381ec7bdfd7bbaab176d22893..0000000000000000000000000000000000000000
--- a/crates/rhodium-market/tests/ui/fx_ref_cannot_escape_view.stderr
+++ /dev/null
@@ -1,12 +0,0 @@
-error: lifetime may not live long enough
- --> tests/ui/fx_ref_cannot_escape_view.rs:5:38
-  |
-5 |     let fx = market.with_view(|view| view.fx(pair, collateral))?;
-  |                                ----- ^^^^^^^^^^^^^^^^^^^^^^^^^ returning this value requires that `'1` must outlive `'2`
-  |                                |   |
-  |                                |   return type of closure is Result<FxRef<'2>, rhodium_core::Error>
-  |                                has type `MarketView<'1>`
-  |
-  = note: requirement occurs because of the type `FxRef<'_>`, which makes the generic argument `'_` invariant
-  = note: the struct `FxRef<'v, R>` is invariant over the parameter `'v`
-  = help: see <https://doc.rust-lang.org/nomicon/subtyping.html> for more information about variance
diff --git a/crates/rhodium-market/tests/ui/growth_point_cross_view.json b/crates/rhodium-market/tests/ui/growth_point_cross_view.json
new file mode 100644
index 0000000000000000000000000000000000000000..e22dd03ec635ecec54fa824bf69ac6a8ce185491
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/growth_point_cross_view.json
@@ -0,0 +1,29 @@
+[
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/growth_point_cross_view.rs",
+        "line": 16,
+        "column": 43
+      }
+    ]
+  },
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/growth_point_cross_view.rs",
+        "line": 16,
+        "column": 43
+      },
+      {
+        "file": "tests/ui/growth_point_cross_view.rs",
+        "line": 16,
+        "column": 43
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/growth_point_cross_view.stderr b/crates/rhodium-market/tests/ui/growth_point_cross_view.stderr
deleted file mode 100644
index 2e085a7c806c1228259bf961cce49532fc585951..0000000000000000000000000000000000000000
--- a/crates/rhodium-market/tests/ui/growth_point_cross_view.stderr
+++ /dev/null
@@ -1,28 +0,0 @@
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/growth_point_cross_view.rs:16:43
-   |
-15 |         let at = bound.projection_growth()[0].points()[0];
-   |             -- `at` declared here, outside of the closure body
-16 |         let _ = market.with_view(|second| second.growth(at));
-   |                                   ------  ^^^^^^^^^^^^^^^^^ `second` escapes the closure body here
-   |                                   |
-   |                                   `second` is only valid in the closure body
-   |
-   = note: requirement occurs because of the type `GrowthAt<'_>`, which makes the generic argument `'_` invariant
-   = note: the struct `GrowthAt<'v, R>` is invariant over the parameter `'v`
-   = help: see <https://doc.rust-lang.org/nomicon/subtyping.html> for more information about variance
-
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/growth_point_cross_view.rs:16:43
-   |
-13 |     market.with_view(|first| {
-   |                       -----
-   |                       |
-   |                       `first` is only valid in the closure body
-   |                       has type `MarketView<'1>`
-...
-16 |         let _ = market.with_view(|second| second.growth(at));
-   |                                           ^^^^^^^^^^^^^^^^^
-   |                                           |
-   |                                           `first` escapes the closure body here
-   |                                           argument requires that `'1` must outlive `'static`
diff --git a/crates/rhodium-market/tests/ui/pay_point_cross_scope.json b/crates/rhodium-market/tests/ui/pay_point_cross_scope.json
new file mode 100644
index 0000000000000000000000000000000000000000..c48e805c1aea62d0fccbb5f415a64567aaa91508
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/pay_point_cross_scope.json
@@ -0,0 +1,29 @@
+[
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/pay_point_cross_scope.rs",
+        "line": 13,
+        "column": 25
+      }
+    ]
+  },
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/pay_point_cross_scope.rs",
+        "line": 13,
+        "column": 25
+      },
+      {
+        "file": "tests/ui/pay_point_cross_scope.rs",
+        "line": 13,
+        "column": 25
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/pay_point_cross_scope.stderr b/crates/rhodium-market/tests/ui/pay_point_cross_scope.stderr
deleted file mode 100644
index 4d25992c2ecc5f9152de8873edcfc462078e2d22..0000000000000000000000000000000000000000
--- a/crates/rhodium-market/tests/ui/pay_point_cross_scope.stderr
+++ /dev/null
@@ -1,28 +0,0 @@
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/pay_point_cross_scope.rs:13:25
-   |
-11 |             let pay = first.at(&point)?;
-   |                 --- `pay` declared here, outside of the closure body
-12 |             view.ccy_scope(Currency::JPY, curve, |second| {
-   |                                                   ------ `second` is only valid in the closure body
-13 |                 let _ = second.pv(100.0, pay);
-   |                         ^^^^^^^^^^^^^^^^^^^^^ `second` escapes the closure body here
-   |
-   = note: requirement occurs because of the type `PayAt<'_>`, which makes the generic argument `'_` invariant
-   = note: the struct `PayAt<'c, R>` is invariant over the parameter `'c`
-   = help: see <https://doc.rust-lang.org/nomicon/subtyping.html> for more information about variance
-
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/pay_point_cross_scope.rs:13:25
-   |
-10 |         view.ccy_scope(Currency::JPY, curve, |first| {
-   |                                               -----
-   |                                               |
-   |                                               `first` is only valid in the closure body
-   |                                               has type `CcyScope<'1, '_>`
-...
-13 |                 let _ = second.pv(100.0, pay);
-   |                         ^^^^^^^^^^^^^^^^^^^^^
-   |                         |
-   |                         `first` escapes the closure body here
-   |                         argument requires that `'1` must outlive `'static`
