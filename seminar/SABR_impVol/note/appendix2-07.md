diff --git a/.githooks/checks.py b/.githooks/checks.py
index f364b0846c043bd13d4921cec45984553bddf26a..71fca10b57e712cd1416f23fe9d2f4150d1014a8 100644
--- a/.githooks/checks.py
+++ b/.githooks/checks.py
@@ -21,6 +21,8 @@ PYTHON_PROJECTS = (".githooks",)
 CARGO_INPUTS = frozenset({"Cargo.toml", "Cargo.lock", "rust-toolchain.toml", "clippy.toml"})
 # examples/ のサンプルをテストとベンチマークで読む crate。CLI は逆依存の閉包で入る。
 EXAMPLE_READERS = ("rhodium",)
+# テストからだけ使う開発専用の crate。src もテストのコードとして扱い、ライブラリの規則を当てない。
+DEV_CRATES = frozenset({"rhodium-testkit", "rhodium-fixtures"})
 
 type _Command = tuple[str | Path, ...]
 type _FileKind = Literal["rust", "python", "toml", "json", "nix", "shell"]
@@ -134,6 +136,7 @@ RULES = (
         pattern=re.compile(r'#\[path\s*=\s*"\.\./'),
         reason="他の crate やディレクトリのファイルを #[path] で取り込まず、共有の補助を使う",
         doc="docs/design/11-quality.md#テストの補助",
+        mode="tree",
         scope="all",
     ),
     _Rule(
@@ -610,9 +613,13 @@ class _Violation:
 
 def _in_scope(path: str, scope: _Scope) -> bool:
     match PurePosixPath(path).parts:
-        case ("crates", _, "src", *rest) if path.endswith(".rs"):
-            # src/ir/tests.rs と src/ir/tests/ の下は単体テストのファイル。
-            return scope == "all" or not any(PurePosixPath(part).stem == "tests" for part in rest)
+        case ("crates", crate, "src", *rest) if path.endswith(".rs"):
+            # src/ir/tests.rs と src/ir/tests/ の下は単体テストのファイル。開発専用の crate の src は
+            # テストの補助なので、テストのファイルと同じく all の規則だけを当てる。
+            return scope == "all" or (
+                crate not in DEV_CRATES
+                and not any(PurePosixPath(part).stem == "tests" for part in rest)
+            )
         case ("crates", _, _, *_) if path.endswith(".rs"):
             return scope == "all"
         case _:
diff --git a/.githooks/test_checks.py b/.githooks/test_checks.py
index e95803aa7090de434a99d10a82b1253919cdd00b..4bd2e52e9a3e7053962f631ea663279099d5c323 100644
--- a/.githooks/test_checks.py
+++ b/.githooks/test_checks.py
@@ -229,7 +229,7 @@ class PatternTest(unittest.TestCase):
         assert [item.line for item in found] == [7]
 
     def test_library_scope_skips_test_code(self) -> None:
-        """ライブラリの規則は、テストのファイルと #[cfg(test)] の item を照合しない。"""
+        """ライブラリの規則は、テストのファイル、#[cfg(test)] の item、開発専用の crate を照合しない。"""
         source = (
             "fn a() {}\n"
             "#[cfg(test)]\n"
@@ -246,9 +246,15 @@ class PatternTest(unittest.TestCase):
             "crates/c/src/ir/tests.rs": "forbidden();\n",
             "crates/c/src/ir/tests/standard.rs": "forbidden();\n",
             "crates/c/tests/it.rs": "forbidden();\n",
+            # 開発専用の crate の src はテストの補助なので、ライブラリの規則を当てない。
+            "crates/rhodium-fixtures/src/ois.rs": "forbidden();\n",
+            "crates/rhodium-testkit/src/lib.rs": "forbidden();\n",
         }
         found = checks._violations([_rule()], files, None)
         assert [(item.path, item.line) for item in found] == [(_LIB, 7), (_LIB, 9)]
+        # all の規則は開発専用の crate も照合する。
+        dev = {"crates/rhodium-fixtures/src/ois.rs": "forbidden();\n"}
+        assert len(checks._violations([_rule(scope="all")], dev, None)) == 1
 
     def test_rust_comments_and_multiline_arguments(self) -> None:
         """コメント、リテラル、入れ子の型を区別し、bool の行だけを報告する。"""
diff --git a/Cargo.lock b/Cargo.lock
index 08a7d72d3a48356c6b9a8acabec189164f3fa956..a32e38e5b1c9598fef83a92cb2a3e845bddaeef1 100644
--- a/Cargo.lock
+++ b/Cargo.lock
@@ -534,6 +534,7 @@ dependencies = [
  "criterion",
  "csv",
  "rhodium-core",
+ "rhodium-fixtures",
  "rhodium-market",
  "rhodium-pricing",
  "rhodium-products",
@@ -554,6 +555,7 @@ dependencies = [
  "clap",
  "csv",
  "rhodium",
+ "rhodium-fixtures",
  "rhodium-testkit",
  "tempfile",
 ]
@@ -571,6 +573,19 @@ dependencies = [
  "thiserror",
 ]
 
+[[package]]
+name = "rhodium-fixtures"
+version = "0.1.0"
+dependencies = [
+ "rhodium-core",
+ "rhodium-market",
+ "rhodium-pricing",
+ "rhodium-products",
+ "serde",
+ "serde_json",
+ "toml",
+]
+
 [[package]]
 name = "rhodium-market"
 version = "0.1.0"
@@ -588,6 +603,7 @@ name = "rhodium-pricing"
 version = "0.1.0"
 dependencies = [
  "rhodium-core",
+ "rhodium-fixtures",
  "rhodium-market",
  "rhodium-products",
  "rhodium-testkit",
@@ -619,6 +635,7 @@ version = "0.1.0"
 dependencies = [
  "proptest",
  "rhodium-core",
+ "rhodium-fixtures",
  "rhodium-market",
  "rhodium-pricing",
  "rhodium-products",
diff --git a/Cargo.toml b/Cargo.toml
index 0592708867cf301ffdf2bd50f926411685af585c..ac398b7d69919c4e7f0c827c199197a895886bf2 100644
--- a/Cargo.toml
+++ b/Cargo.toml
@@ -9,6 +9,7 @@ members = [
     "crates/rhodium",
     "crates/rhodium-cli",
     "crates/rhodium-testkit",
+    "crates/rhodium-fixtures",
 ]
 resolver = "3"
 
@@ -40,6 +41,7 @@ rhodium-products = { path = "crates/rhodium-products" }
 rhodium-risk = { path = "crates/rhodium-risk" }
 rhodium-xva = { path = "crates/rhodium-xva" }
 rhodium-testkit = { path = "crates/rhodium-testkit" }
+rhodium-fixtures = { path = "crates/rhodium-fixtures" }
 serde = { version = "1", features = ["derive"] }
 thiserror = "2"
 
@@ -53,6 +55,7 @@ unused_must_use = "deny"
 [workspace.lints.clippy]
 all = { level = "deny", priority = -1 }
 pedantic = { level = "warn", priority = -1 }
+allow_attributes = "deny"
 allow_attributes_without_reason = "deny"
 cast_precision_loss = "allow"
 dbg_macro = "deny"
diff --git a/README.md b/README.md
index 9955a57930fdc3531ace827970f3d15ace7679d5..f2fbfd9c9f12559fefa7d4550081564103944139 100644
--- a/README.md
+++ b/README.md
@@ -34,7 +34,7 @@ FRA、金利先物、IRS、ベーシススワップ、OIS（lookback、observati
 所要時間は [M0 の計測](docs/benchmarks/m0.md) と [M1 の計測](docs/benchmarks/m1.md) にあります。
 次のマイルストーンは M2（ボラと解析オプション）です。
 
-ワークスペースは上の 8 crate と開発用の `rhodium-testkit` からなり、旧設計の crate は残っていません。
+ワークスペースは上の 8 crate と、開発専用の `rhodium-testkit`（テストの数値の補助）と `rhodium-fixtures`（共有のテストの fixture）からなり、旧設計の crate は残っていません。
 
 ## サンプル
 
diff --git a/crates/rhodium-cli/Cargo.toml b/crates/rhodium-cli/Cargo.toml
index 22fcf937af17be25015a72a4bda1e8ea5e239c5f..36516d32b927b6559bbd1d89d8a0f80f8b63ca55 100644
--- a/crates/rhodium-cli/Cargo.toml
+++ b/crates/rhodium-cli/Cargo.toml
@@ -20,5 +20,6 @@ csv.workspace = true
 rhodium.workspace = true
 
 [dev-dependencies]
+rhodium-fixtures.workspace = true
 rhodium-testkit.workspace = true
 tempfile = "3.27"
diff --git a/crates/rhodium-cli/tests/workspace_cli.rs b/crates/rhodium-cli/tests/workspace_cli.rs
index ebe28a7fc7807a03fdcb6efec5bcc6a87e145fd3..4e5adb592bdae5069caf6ac5156d7ebe944f80c9 100644
--- a/crates/rhodium-cli/tests/workspace_cli.rs
+++ b/crates/rhodium-cli/tests/workspace_cli.rs
@@ -1,48 +1,27 @@
 //! 実ファイルのワークスペースをファサードと rh から評価する。
 
 use std::{
+    collections::{BTreeMap, BTreeSet},
     error::Error,
     fmt::Write as _,
     fs, io,
-    num::NonZeroUsize,
     path::{Path, PathBuf},
     process::{Command, Output},
 };
 
-use rhodium::{
-    prelude::{ExecutionContext, ThreadExecutor},
-    session::Session,
-    workspace::Workspace,
-};
-use rhodium_testkit::assert_close;
+use rhodium::{prelude::ExecutionContext, session::Session, workspace::Workspace};
+use rhodium_fixtures::files::{copy_dir, example};
+use rhodium_testkit::{PV_REL, SENSITIVITY_REL, assert_close, minor_unit_abs};
 use tempfile::TempDir;
 
 type TestResult<T = ()> = Result<T, Box<dyn Error>>;
 type Rows = Vec<(String, String, f64)>;
 
-#[expect(
-    clippy::disallowed_macros,
-    reason = "Cargo が渡す、試験データへのパスを取得する"
-)]
 fn sample() -> PathBuf {
-    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m0-ois")
+    example("m0-ois")
 }
 
-fn copy_dir(from: &Path, to: &Path) -> io::Result<()> {
-    for entry in fs::read_dir(from)? {
-        let entry = entry?;
-        let target = to.join(entry.file_name());
-        if entry.file_type()?.is_dir() {
-            fs::create_dir(&target)?;
-            copy_dir(&entry.path(), &target)?;
-        } else {
-            fs::copy(entry.path(), target)?;
-        }
-    }
-    Ok(())
-}
-
-fn workspace_copy() -> io::Result<TempDir> {
+fn workspace_copy() -> TestResult<TempDir> {
     let dir = tempfile::tempdir()?;
     copy_dir(&sample(), dir.path())?;
     Ok(dir)
@@ -72,18 +51,6 @@ fn rows(output: &Output, headers: &[&str]) -> TestResult<Rows> {
         .collect()
 }
 
-fn diagnostic(dir: &Path) -> TestResult<String> {
-    load_diagnostic(Workspace::load(dir))
-}
-
-fn load_diagnostic(workspace: rhodium::prelude::Result<Workspace>) -> TestResult<String> {
-    match workspace {
-        Err(error @ rhodium::prelude::Error::InvalidData(_)) => Ok(error.to_string()),
-        Err(error) => Err(error.into()),
-        Ok(_) => Err(io::Error::other("expected Workspace::load to fail").into()),
-    }
-}
-
 fn ois(notional: f64, accrual: f64, fixed: f64, quote: f64) -> (f64, f64) {
     // 期首で評価する単一支払いの OIS は DF = 1 / (1 + q α)。
     // 固定受けの PV を q で微分し、金利 1 bp あたりに直す。
@@ -120,8 +87,12 @@ fn price_and_risk_match_for_as_of_and_valuation_close() -> TestResult {
         &["--date", "2026-11-24"],
     ];
     for (command, headers, relative) in [
-        ("price", ["trade_id", "currency", "pv"], 1e-12),
-        ("risk", ["quote_id", "currency", "sensitivity"], 1e-10),
+        ("price", ["trade_id", "currency", "pv"], PV_REL),
+        (
+            "risk",
+            ["quote_id", "currency", "sensitivity"],
+            SENSITIVITY_REL,
+        ),
     ] {
         let outputs = times
             .iter()
@@ -278,142 +249,6 @@ fn cli_reports_workspace_load_errors_as_input_errors() -> TestResult {
     Ok(())
 }
 
-#[test]
-fn session_reuses_the_same_instant_for_price_and_risk() -> TestResult {
-    // 同じセッションから両方を呼び出し、取引 ID と通貨別の結果を確認する。
-    fn send_sync<T: Send + Sync>() {}
-    send_sync::<Session>();
-    let workspace = Workspace::load(sample())?;
-    let as_of = workspace.resolve_date("2026-11-24".parse()?)?;
-    assert_eq!(as_of.to_string(), "2026-11-24T08:00:00Z");
-    let execution = ExecutionContext::sequential();
-    let session = Session::new(workspace, as_of, execution)?;
-    assert_eq!(session.as_of(), as_of);
-    assert_eq!(session.calibrated_market().groups().len(), 2);
-    let prices = session.pv(execution)?;
-    assert_eq!(prices.len(), 2);
-    assert!(prices.values().all(Result::is_ok));
-    let risk = session.quote_delta(execution)?;
-    assert!(risk.failures().is_empty());
-    assert_eq!(risk.by_ccy().len(), 2);
-    Ok(())
-}
-
-#[test]
-fn workspace_reports_multiple_csv_errors_with_lines() -> TestResult {
-    // レコード長の違い、重複 ID、非有限値のすべてを一度に報告する。
-    let dir = workspace_copy()?;
-    fs::write(
-        dir.path().join("quotes/snapshot.csv"),
-        "quote_id,value\nJPY.TONA.OIS.1M,45bp\nJPY.TONA.OIS.1M,0.006\nTOO-SHORT\nUSD.SOFR.OIS.1M,NaN\n",
-    )?;
-    let error = diagnostic(dir.path())?;
-    assert!(error.contains("quotes/snapshot.csv"), "{error}");
-    assert!(
-        error.contains("line 3") && error.contains("duplicate ID"),
-        "{error}"
-    );
-    assert!(error.contains("line: 4"), "{error}");
-    assert!(
-        error.contains("line 5") && error.contains("value: expected a finite decimal"),
-        "{error}"
-    );
-    Ok(())
-}
-
-#[test]
-fn workspace_reports_multiple_toml_rows_and_duplicate_trade_ids() -> TestResult {
-    // 同じファイルの別々の不正な行と、ファイルをまたぐ取引 ID の重複を検出する。
-    let dir = workspace_copy()?;
-    let path = dir.path().join("trades/swaps.toml");
-    let later = dir.path().join("trades/zz-duplicate.toml");
-    let good = fs::read_to_string(&path)?;
-    fs::write(&later, "]]")?;
-    fs::write(
-        &path,
-        good.replace("book = \"USD-SWAPS\"", "bok = \"USD-SWAPS\"")
-            .replace(
-                "margin_set = \"JPY-CASH\"",
-                "margin_set = \"NO-SUCH-MARGIN\"",
-            ),
-    )?;
-    let error = diagnostic(dir.path())?;
-    // 後の文書が壊れていても、先に読み込んだ文書の行エラーを元の順で返す。
-    let executor = ThreadExecutor::new(NonZeroUsize::MIN.saturating_add(3));
-    let parallel = || {
-        load_diagnostic(Workspace::load_with_execution(
-            dir.path(),
-            ExecutionContext::new(&executor),
-        ))
-    };
-    assert_eq!(parallel()?, error);
-    assert!(error.contains("trades/swaps.toml:3:1: trade[0]"), "{error}");
-    assert!(error.contains("unknown field `bok`"), "{error}");
-    assert!(
-        error.contains("trades/swaps.toml:18:1: trade[1]"),
-        "{error}"
-    );
-    assert!(
-        error.contains("margin_set: unknown margin set NO-SUCH-MARGIN"),
-        "{error}"
-    );
-    fs::write(&path, &good)?;
-    fs::write(later, good)?;
-    let error = diagnostic(dir.path())?;
-    assert_eq!(parallel()?, error);
-    assert!(
-        error.contains("zz-duplicate.toml") && error.contains("duplicate ID JPY-SWAP"),
-        "{error}"
-    );
-    Ok(())
-}
-
-#[test]
-fn workspace_rejects_invalid_version_and_fixing_header() -> TestResult {
-    // 版を解釈できない文書と、日次指標に日中の列を書いたファイルを拒否する。
-    let dir = workspace_copy()?;
-    let path = dir.path().join("agreements.toml");
-    let good = fs::read_to_string(&path)?;
-    fs::write(&path, good.replace("version = 1", "version = 2"))?;
-    let error = diagnostic(dir.path())?;
-    assert!(
-        error.contains("agreements.toml") && error.contains("unsupported version 2"),
-        "{error}"
-    );
-    fs::write(&path, good)?;
-    fs::write(
-        dir.path().join("fixings/JPY-TONA.csv"),
-        "at,value\n2026-11-20T08:00:00Z,0.0046\n",
-    )?;
-    let error = diagnostic(dir.path())?;
-    assert!(
-        error.contains("fixings/JPY-TONA.csv") && error.contains("expected CSV header date,value"),
-        "{error}"
-    );
-    Ok(())
-}
-
-#[test]
-fn workspace_reports_reference_errors_by_source_entry() -> TestResult {
-    // 別々の指標が参照する未知カレンダーを、入力ファイルと指標 ID で報告する。
-    let dir = workspace_copy()?;
-    let path = dir.path().join("refdata/indices.toml");
-    let indices = fs::read_to_string(&path)?
-        .replace("calendar = \"JPTO\"", "calendar = \"UNKNOWN-JP\"")
-        .replace("calendar = \"USNY\"", "calendar = \"UNKNOWN-US\"");
-    fs::write(path, indices)?;
-    let error = diagnostic(dir.path())?;
-    assert!(
-        error.contains("refdata/indices.toml: index[JPY-TONA]") && error.contains("UNKNOWN-JP"),
-        "{error}"
-    );
-    assert!(
-        error.contains("refdata/indices.toml: index[USD-SOFR]") && error.contains("UNKNOWN-US"),
-        "{error}"
-    );
-    Ok(())
-}
-
 #[test]
 fn cli_requires_exactly_one_valuation_time() -> TestResult {
     // clap の引数エラーで入力の読み込みや CSV の出力に進まないことを確認する。
@@ -435,43 +270,6 @@ fn cli_requires_exactly_one_valuation_time() -> TestResult {
     Ok(())
 }
 
-#[test]
-fn workspace_rejects_ambiguous_or_missing_local_times() -> TestResult {
-    // 同じ地域の重複時刻は明示オフセットで区別し、欠落時刻は拒否する。
-    let workspace = Workspace::load(sample())?;
-    let instant = workspace.resolve_as_of("2026-11-24T17:00:00[Asia/Tokyo]")?;
-    assert_eq!(
-        instant,
-        workspace.resolve_as_of("2026-11-24T17:00:00+09:00[Asia/Tokyo]")?
-    );
-    assert!(
-        workspace
-            .resolve_as_of("2026-11-24T17:00:00+08:00[Asia/Tokyo]")
-            .is_err()
-    );
-    assert!(workspace.resolve_as_of("2026-11-24T17:00:00").is_err());
-    assert!(
-        workspace
-            .resolve_as_of("2026-03-08T02:30:00[America/New_York]")
-            .is_err()
-    );
-    assert!(
-        workspace
-            .resolve_as_of("2026-11-01T01:30:00[America/New_York]")
-            .is_err()
-    );
-    let early = workspace.resolve_as_of("2026-11-01T01:30:00-04:00[America/New_York]")?;
-    let late = workspace.resolve_as_of("2026-11-01T01:30:00-05:00[America/New_York]")?;
-    assert_eq!(early.to_string(), "2026-11-01T05:30:00Z");
-    assert_eq!(late.to_string(), "2026-11-01T06:30:00Z");
-    // 参照データに無い地域は、他の地域の参照と同じく参照の欠けとして扱う。
-    assert!(matches!(
-        workspace.resolve_as_of("2026-11-24T17:00:00[Europe/Paris]"),
-        Err(rhodium::prelude::Error::InvalidData(message)) if message.contains("Europe/Paris")
-    ));
-    Ok(())
-}
-
 #[test]
 fn date_requires_an_unambiguous_configured_close() -> TestResult {
     // 日付指定には締め時刻が必要で、夏時間の重複と欠落は計算に進む前に拒否する。
@@ -827,8 +625,7 @@ fn cli_warns_about_rate_quotes_written_without_a_unit() -> TestResult {
 
 #[test]
 fn cli_thread_count_keeps_price_and_risk_csv() -> TestResult {
-    // サンプルに 1M の取引を 8 件足し、4 スレッドで取引を区間に分けても同じ行と値を出すことを確かめる。
-    // 取引ごとの PV は分け方によらず、デルタは区間ごとの和の順だけが変わる。
+    // サンプルに 1M の取引を 8 件足し、--threads 4 で取引を区間に分けても同じ形の表を出すことを確かめる。
     let dir = workspace_copy()?;
     let mut trades = String::from("version = 1\n");
     for (i, (ccy, convention, unit, side)) in [
@@ -864,9 +661,9 @@ fixed_side = "{side}"
         )?;
     }
     fs::write(dir.path().join("trades/extra.toml"), trades)?;
-    for (command, headers, relative) in [
-        ("price", ["trade_id", "currency", "pv"], 1e-12),
-        ("risk", ["quote_id", "currency", "sensitivity"], 1e-10),
+    for (command, headers) in [
+        ("price", ["trade_id", "currency", "pv"]),
+        ("risk", ["quote_id", "currency", "sensitivity"]),
     ] {
         let sequential = run(command, dir.path(), &["--date", "2026-11-24"])?;
         let parallel = run(
@@ -881,17 +678,17 @@ fixed_side = "{side}"
                 String::from_utf8_lossy(&output.stderr)
             );
         }
-        let expected = rows(&sequential, &headers)?;
-        let actual = rows(&parallel, &headers)?;
+        // 値がスレッド数によらないことはファサードのテスト（parallel.rs）が確かめる。CLI は同じ行
+        // （ID と通貨の列の全体）を同じ順に書くことだけを見る。
+        let keys = |output| -> TestResult<Vec<(String, String)>> {
+            Ok(rows(output, &headers)?
+                .into_iter()
+                .map(|(id, ccy, _)| (id, ccy))
+                .collect())
+        };
+        let expected = keys(&sequential)?;
         assert_eq!(expected.len(), if command == "price" { 10 } else { 2 });
-        assert_eq!(actual.len(), expected.len());
-        for ((id, ccy, value), (expected_id, expected_ccy, expected_value)) in
-            actual.iter().zip(&expected)
-        {
-            assert_eq!((id, ccy), (expected_id, expected_ccy));
-            let absolute = if ccy == "USD" { 1e-8 } else { 1e-6 };
-            assert_close!(*value, *expected_value, relative, absolute);
-        }
+        assert_eq!(keys(&parallel)?, expected);
     }
     // ワーカーの数は 1 以上。0 は引数の誤りとして終了コード 2 になる。
     let zero = run(
@@ -1012,70 +809,97 @@ day_count = "ACT360"
     Ok(())
 }
 
-#[expect(
-    clippy::disallowed_macros,
-    reason = "リポジトリのサンプルの場所をテストのソースの位置から決める"
-)]
-fn m1_sample() -> PathBuf {
-    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m1-rates")
+/// 出力の数値の行を、品質の許容誤差で比べる表にする。CSV の 10 進の表現を往復しても値は変わらない。
+fn table(rows: &Rows) -> BTreeMap<(String, String), f64> {
+    rows.iter()
+        .map(|(id, ccy, value)| ((id.clone(), ccy.clone()), *value))
+        .collect()
 }
 
 #[test]
 fn cli_prices_and_risks_the_m1_rates_workspace() -> TestResult {
-    // 最初の検証データセットのサンプルを rh price と rh risk で評価する。全 18 取引に PV の行があり、
-    // 近似の警告は予想の測度に従う（TONA 担保の TIBOR の IRS は CrossCollateralProjection を出さず、
-    // 25bp の CSA の IRS は出す）。ワーカー数を変えても CSV の値は許容誤差で一致する。
-    let dir = m1_sample();
+    // 最初の検証データセットのサンプルを rh price と rh risk で評価し、出力の表と診断がファサードの
+    // 結果をそのまま書いたものであることを、行の集合と値の全体で比べる。取引ごとの近似の警告の意味
+    // （予想の測度など）とスレッド数によらないことは、ファサードのテスト（m1_rates.rs と parallel.rs）が
+    // 確かめるので、ここでは書き出しの形と中身の対応だけを見る。
+    let dir = example("m1-rates");
+    let workspace = Workspace::load(&dir)?;
+    let as_of = workspace.resolve_date("2026-11-24".parse()?)?;
+    let execution = ExecutionContext::sequential();
+    let session = Session::new(workspace, as_of, execution)?;
+
     let price = run("price", &dir, &["--date", "2026-11-24"])?;
     assert!(
         price.status.success(),
         "{}",
         String::from_utf8_lossy(&price.stderr)
     );
-    let pvs = rows(&price, &["trade_id", "currency", "pv"])?;
-    let trades: std::collections::BTreeSet<_> = pvs.iter().map(|(id, _, _)| id.as_str()).collect();
-    assert_eq!(trades.len(), 18);
+    let pvs = table(&rows(&price, &["trade_id", "currency", "pv"])?);
+    let expected = session
+        .pv(execution)?
+        .into_iter()
+        .map(|(trade, pv)| {
+            Ok(pv?
+                .terms()
+                .map(|money| ((trade.to_string(), money.ccy().to_string()), money.amount()))
+                .collect::<Vec<_>>())
+        })
+        .collect::<TestResult<Vec<_>>>()?
+        .into_iter()
+        .flatten()
+        .collect::<BTreeMap<_, _>>();
+    assert_eq!(
+        pvs.keys().collect::<Vec<_>>(),
+        expected.keys().collect::<Vec<_>>()
+    );
+    for (key, value) in &pvs {
+        assert_close!(*value, expected[key], PV_REL, minor_unit_abs(0); "{key:?}: ");
+    }
+    // 準備の警告は、取引ごとに "warning: trade <ID>: <警告>" の一行で、すべてを書く。
     let stderr = String::from_utf8(price.stderr)?;
     assert!(!stderr.contains("error"), "{stderr}");
-    assert!(stderr.contains("trade JPY-IRS-5Y-CASH-25: CrossCollateralProjection"));
-    assert!(!stderr.contains("trade JPY-IRS-5Y: CrossCollateralProjection"));
-    assert!(stderr.contains("trade USD-OIS-JPY-CASH: CrossCollateralProjection"));
+    let warned: BTreeSet<&str> = stderr
+        .lines()
+        .filter(|line| line.starts_with("warning: trade "))
+        .collect();
+    let expected_warnings: Vec<String> = session
+        .preparation_warnings()
+        .flat_map(|(id, warnings)| {
+            warnings
+                .iter()
+                .map(move |warning| format!("warning: trade {id}: {warning:?}"))
+        })
+        .collect();
+    assert_eq!(
+        warned,
+        expected_warnings.iter().map(String::as_str).collect()
+    );
+    assert!(warned.contains("warning: trade JPY-IRS-5Y-CASH-25: CrossCollateralProjection"));
 
     let risk = run("risk", &dir, &["--date", "2026-11-24"])?;
     assert!(risk.status.success());
-    let deltas = rows(&risk, &["quote_id", "currency", "sensitivity"])?;
-    for quote in [
-        "JPY.TONA.OIS.5Y",
-        "JPY.TONA.TURN.2026",
-        "JPY.TIBOR6M-TONA.BASIS.5Y",
-        "JPY.TIBOR3M-TONA.BASIS.5Y",
-        "USD.SOFR.OIS.2Y",
-        "USDJPY.SPOT",
-        "USDJPY.XCCY.2Y",
-        "JPY.LGM.VOL",
-    ] {
-        assert!(
-            deltas
+    let deltas = table(&rows(&risk, &["quote_id", "currency", "sensitivity"])?);
+    let report = session.quote_delta(execution)?;
+    assert!(report.failures().is_empty());
+    let expected = report
+        .by_ccy()
+        .iter()
+        .flat_map(|vector| {
+            vector
+                .entries()
                 .iter()
-                .any(|(id, _, value)| id == quote && value.abs() > 1e-6),
-            "{quote}"
-        );
-    }
-
-    for (command, values, headers) in [
-        ("price", &pvs, ["trade_id", "currency", "pv"]),
-        ("risk", &deltas, ["quote_id", "currency", "sensitivity"]),
-    ] {
-        let parallel = run(command, &dir, &["--date", "2026-11-24", "--threads", "2"])?;
-        assert!(parallel.status.success());
-        let parallel = rows(&parallel, &headers)?;
-        assert_eq!(parallel.len(), values.len());
-        for ((id, ccy, value), (parallel_id, parallel_ccy, parallel_value)) in
-            values.iter().zip(&parallel)
-        {
-            assert_eq!((id, ccy), (parallel_id, parallel_ccy));
-            assert_close!(*value, *parallel_value, 1e-10_f64, 1e-6_f64);
-        }
+                .map(|(quote, value)| ((quote.to_string(), vector.ccy().to_string()), *value))
+        })
+        .collect::<BTreeMap<_, _>>();
+    // 校正の全グループ（TONA、年末の跳び、TIBOR のベーシス、SOFR、FX、XCCY、LGM のマーク）の quote が
+    // 行に現れる。
+    assert!(expected.len() > 8);
+    assert_eq!(
+        deltas.keys().collect::<Vec<_>>(),
+        expected.keys().collect::<Vec<_>>()
+    );
+    for (key, value) in &deltas {
+        assert_close!(*value, expected[key], SENSITIVITY_REL, minor_unit_abs(0); "{key:?}: ");
     }
     Ok(())
 }
@@ -1084,7 +908,7 @@ fn cli_prices_and_risks_the_m1_rates_workspace() -> TestResult {
 fn calibrate_reads_rate_volatility_and_reports_partial_failures() -> TestResult {
     // 同じワークスペースで、SABR の ATM とターム・RFR のキャップを校正する。
     let dir = tempfile::tempdir()?;
-    copy_dir(&sample().with_file_name("m1-rates"), dir.path())?;
+    copy_dir(&example("m1-rates"), dir.path())?;
     let path = dir.path().join("market/definition.toml");
     let definitions = r#"
 
diff --git a/crates/rhodium-core/src/date.rs b/crates/rhodium-core/src/date.rs
index bb60bd8fb96d11e8886482d4d8d8b82c12df25e9..431969db5abb73ed041e87cc8ae3a73b3c6782af 100644
--- a/crates/rhodium-core/src/date.rs
+++ b/crates/rhodium-core/src/date.rs
@@ -43,7 +43,7 @@ impl Date {
 
     /// 年月日。
     #[must_use]
-    #[allow(clippy::cast_possible_truncation, reason = "月は 12、日は 31 以下")]
+    #[expect(clippy::cast_possible_truncation, reason = "月は 12、日は 31 以下")]
     pub fn ymd(self) -> (i32, u8, u8) {
         (self.0.year(), self.0.month() as u8, self.0.day() as u8)
     }
@@ -69,7 +69,7 @@ impl Date {
 
     /// 二つの日付の間の符号付き暦日数（self − earlier）。
     #[must_use]
-    #[allow(
+    #[expect(
         clippy::cast_possible_truncation,
         reason = "表せる全期間でも i32 の暦日数に収まる"
     )]
diff --git a/crates/rhodium-core/src/schedule/grid.rs b/crates/rhodium-core/src/schedule/grid.rs
index 5309833a97fab95666dab5ec4a0016fb3bae319c..c473170f763647e41662ca0c4f137695bae907ac 100644
--- a/crates/rhodium-core/src/schedule/grid.rs
+++ b/crates/rhodium-core/src/schedule/grid.rs
@@ -482,7 +482,7 @@ fn apply_roll(date: Date, roll: RollConvention, unit: TenorUnit) -> Result<Date>
     }
 }
 
-#[allow(clippy::cast_possible_wrap, reason = "曜日番号の範囲は 0〜6")]
+#[expect(clippy::cast_possible_wrap, reason = "曜日番号の範囲は 0〜6")]
 fn weekday_index(weekday: chrono::Weekday) -> i32 {
     weekday.num_days_from_monday() as i32
 }
diff --git a/crates/rhodium-core/tests/day_count.rs b/crates/rhodium-core/tests/day_count.rs
index 07e478e8b09f335811f29cd07a3fe340124446f7..40e6b7426b4005154a123fcacd14a544e9bc5189 100644
--- a/crates/rhodium-core/tests/day_count.rs
+++ b/crates/rhodium-core/tests/day_count.rs
@@ -11,6 +11,7 @@ use rhodium_core::{
     quantity::YearFraction,
     tenor::Frequency,
 };
+use rhodium_testkit::{MARKET_ABS, MARKET_REL, assert_close};
 
 fn period(start: &str, end: &str) -> Result<DateInterval> {
     DateInterval::new(start.parse()?, end.parse()?)
@@ -30,14 +31,6 @@ fn calendar(id: &str, extrapolation: Extrapolation) -> Result<Calendar> {
     )
 }
 
-fn assert_close(actual: f64, expected: f64, dc: &DayCount) {
-    let scale = actual.abs().max(expected.abs());
-    assert!(
-        (actual - expected).abs() <= (1e-12 * scale).max(1e-14),
-        "{dc}: {actual} != {expected}"
-    );
-}
-
 #[test]
 #[expect(
     clippy::too_many_lines,
@@ -155,7 +148,7 @@ fn all_day_counts_match_hand_calculations() -> Result<()> {
         };
         let (fraction, warnings) = year_fraction(&dc, accrual.start(), accrual.end(), &ctx)?;
         assert_eq!(fraction.basis(), &dc);
-        assert_close(fraction.value(), expected, &dc);
+        assert_close!(fraction.value(), expected, MARKET_REL, MARKET_ABS; "{}: ", dc);
         assert_eq!(warnings, []);
 
         // 空区間では暦、満期、頻度が渡されていなくても 0 で、タグを保つ。
@@ -165,18 +158,39 @@ fn all_day_counts_match_hand_calculations() -> Result<()> {
             accrual.start(),
             &AccrualContext::default(),
         )?;
-        assert_close(empty.value(), 0.0, &dc);
+        assert_close!(empty.value(), 0.0, MARKET_REL, MARKET_ABS; "{}: ", dc);
         assert_eq!(empty.basis(), &dc);
         assert_eq!(warnings, []);
     }
     Ok(())
 }
 
+/// 検算の表の行が読む参照枠。枠を読むのは ACT/ACT ICMA と ACT/365L だけ。
+fn oracle_frame<'a>(
+    dc: &DayCount,
+    icma: &'a PeriodFrame,
+    short_stub: &'a PeriodFrame,
+) -> Option<&'a PeriodFrame> {
+    match dc {
+        DayCount::ActActIcma => Some(icma),
+        DayCount::Act365L => Some(short_stub),
+        DayCount::Act360
+        | DayCount::Act365F
+        | DayCount::Act364
+        | DayCount::ActActIsda
+        | DayCount::ActActAfb
+        | DayCount::Nl365
+        | DayCount::Thirty360BondBasis
+        | DayCount::Thirty360Us
+        | DayCount::ThirtyE360
+        | DayCount::ThirtyE360Isda
+        | DayCount::ThirtyEPlus360
+        | DayCount::Bus252 { .. }
+        | DayCount::OneOne => None,
+    }
+}
+
 #[test]
-#[expect(
-    clippy::wildcard_enum_match_arm,
-    reason = "#196 で各 variant を網羅する"
-)]
 fn day_count_fractions_agree_with_published_oracles() -> Result<()> {
     // Strata の公開検算（固定コミット 910a144fc49532c1411ed95137556476d622c85c）の期間と値。
     // BUS/252 だけは QuantLib v1.40 の最初の例で、金曜から月曜までの区間で金曜だけを数える。
@@ -269,12 +283,7 @@ fn day_count_fractions_agree_with_published_oracles() -> Result<()> {
         "6M".parse()?,
     )?;
     for (dc, accrual, expected) in cases {
-        // 枠を読むのは ACT/ACT ICMA と ACT/365L だけ。
-        let frame = match dc {
-            DayCount::ActActIcma => Some(&icma_frame),
-            DayCount::Act365L => Some(&short_stub),
-            _ => None,
-        };
+        let frame = oracle_frame(&dc, &icma_frame, &short_stub);
         let ctx = AccrualContext {
             frame,
             month_end: MonthEndSchedule::Aligned,
@@ -282,7 +291,7 @@ fn day_count_fractions_agree_with_published_oracles() -> Result<()> {
             calendar: Some(&bus_calendar),
         };
         let (fraction, warnings) = year_fraction(&dc, accrual.start(), accrual.end(), &ctx)?;
-        assert_close(fraction.value(), expected, &dc);
+        assert_close!(fraction.value(), expected, MARKET_REL, MARKET_ABS; "{}: ", dc);
         assert_eq!(warnings, []);
     }
     Ok(())
@@ -360,7 +369,7 @@ fn documented_day_count_examples_hold() -> Result<()> {
         ),
     ] {
         let (fraction, warnings) = year_fraction(&dc, accrual.start(), accrual.end(), &ctx)?;
-        assert_close(fraction.value(), expected, &dc);
+        assert_close!(fraction.value(), expected, MARKET_REL, MARKET_ABS; "{}: ", dc);
         assert_eq!(warnings, []);
     }
     Ok(())
@@ -394,7 +403,7 @@ fn leap_day_boundaries_and_reference_periods_follow_each_basis() -> Result<()> {
             &AccrualContext::default(),
         )?
         .0;
-        assert_close(actual.value(), expected, &dc);
+        assert_close!(actual.value(), expected, MARKET_REL, MARKET_ABS; "{}: ", dc);
     }
 
     // 年払いは利払い区間が閏日を終端として含む場合に 366。年に複数回なら終了年が閏年かで決める。
@@ -432,12 +441,13 @@ fn leap_day_boundaries_and_reference_periods_follow_each_basis() -> Result<()> {
             ..AccrualContext::default()
         };
         let dc = DayCount::Act365L;
-        assert_close(
+        assert_close!(
             year_fraction(&dc, sampled.start(), sampled.end(), &ctx)?
                 .0
                 .value(),
             expected,
-            &dc,
+            MARKET_REL,
+            MARKET_ABS
         );
     }
 
@@ -457,7 +467,7 @@ fn leap_day_boundaries_and_reference_periods_follow_each_basis() -> Result<()> {
         ..AccrualContext::default()
     };
     let actual = year_fraction(&dc, accrued.start(), accrued.end(), &ctx)?.0;
-    assert_close(actual.value(), 60.0 / 364.0 + 43.0 / 364.0, &dc);
+    assert_close!(actual.value(), 60.0 / 364.0 + 43.0 / 364.0, MARKET_REL, MARKET_ABS; "{}: ", dc);
     Ok(())
 }
 
@@ -529,7 +539,7 @@ fn day_count_context_errors_and_calendar_warnings_propagate() -> Result<()> {
     ));
     // 月曜と水曜の 2 日を数える。火曜は休日、水曜だけが外挿される。
     let (actual, warnings) = year_fraction(&dc_bus, start, end, &ctx(&extra))?;
-    assert_close(actual.value(), 2.0 / 252.0, &dc_bus);
+    assert_close!(actual.value(), 2.0 / 252.0, MARKET_REL, MARKET_ABS; "{}: ", dc_bus);
     assert_eq!(
         warnings,
         [Warning::CalendarExtrapolated {
@@ -553,10 +563,11 @@ fn thirty_day_counts_apply_february_maturity_and_december_rollover() -> Result<(
     };
     for (ctx, expected) in [(aligned, 1.0), (AccrualContext::default(), 361.0 / 360.0)] {
         let dc = DayCount::Thirty360Us;
-        assert_close(
+        assert_close!(
             year_fraction(&dc, first, last, &ctx)?.0.value(),
             expected,
-            &dc,
+            MARKET_REL,
+            MARKET_ABS
         );
     }
     // 2 月の満期だけ元の日を使い、それ以外の月末は 30 日に置き換える。
@@ -566,10 +577,11 @@ fn thirty_day_counts_apply_february_maturity_and_december_rollover() -> Result<(
             maturity: Some(maturity),
             ..AccrualContext::default()
         };
-        assert_close(
+        assert_close!(
             year_fraction(&dc, first, last, &ctx)?.0.value(),
             expected,
-            &dc,
+            MARKET_REL,
+            MARKET_ABS
         );
     }
     // 12 月 31 日は翌年 1 月 1 日に置き換えるので 30 + 1 日。
@@ -581,6 +593,6 @@ fn thirty_day_counts_apply_february_maturity_and_december_rollover() -> Result<(
         &AccrualContext::default(),
     )?
     .0;
-    assert_close(actual.value(), 31.0 / 360.0, &dc);
+    assert_close!(actual.value(), 31.0 / 360.0, MARKET_REL, MARKET_ABS; "{}: ", dc);
     Ok(())
 }
diff --git a/crates/rhodium-core/tests/refdata_and_fixing.rs b/crates/rhodium-core/tests/refdata_and_fixing.rs
index bbbb830b70317d0d5bb36b162f03ae4a7c03ce16..1b882a3c7a59f29d1d12dad27bf97d960ae12b41 100644
--- a/crates/rhodium-core/tests/refdata_and_fixing.rs
+++ b/crates/rhodium-core/tests/refdata_and_fixing.rs
@@ -554,9 +554,12 @@ fn fx_parts() -> Result<RefDataParts> {
             Extrapolation::Strict,
         )?,
     );
-    if let Some(jpy) = parts.currencies.get_mut(&Currency::JPY) {
-        jpy.settlement_calendar = tokyo;
-    }
+    // fixture のキーが変わったときに書き換えを黙って飛ばさず、失敗させる。
+    parts
+        .currencies
+        .get_mut(&Currency::JPY)
+        .ok_or_else(|| Error::MissingData("fixture has no JPY currency".into()))?
+        .settlement_calendar = tokyo;
     Ok(parts)
 }
 
@@ -663,9 +666,12 @@ fn refdata_rejects_fx_convention_key_mismatch_and_unknown_zone() -> Result<()> {
     )));
 
     let mut parts = fx_parts()?;
-    if let Some(convention) = parts.fx_conventions.get_mut(&usd_jpy) {
-        convention.settlement.zone = "Asia/Tokyo".parse()?;
-    }
+    parts
+        .fx_conventions
+        .get_mut(&usd_jpy)
+        .ok_or_else(|| Error::MissingData("fixture has no USDJPY FX convention".into()))?
+        .settlement
+        .zone = "Asia/Tokyo".parse()?;
     let Err(errors) = RefData::new_with_diagnostics(parts) else {
         return Err(Error::Internal("unknown zone was accepted".into()));
     };
diff --git a/crates/rhodium-core/tests/schedules.rs b/crates/rhodium-core/tests/schedules/main.rs
similarity index 97%
rename from crates/rhodium-core/tests/schedules.rs
rename to crates/rhodium-core/tests/schedules/main.rs
index bfe89b390aee666fcf6f5f383029fbfb049cc086..36f10c7fa54aa19626bbb9e4c6bb98646b1e8a5e 100644
--- a/crates/rhodium-core/tests/schedules.rs
+++ b/crates/rhodium-core/tests/schedules/main.rs
@@ -1,6 +1,5 @@
 //! 計算期間と決済の瞬間を独立に求めた値と照合する。
 
-#[path = "schedules/observation.rs"]
 mod observation;
 
 use std::collections::BTreeMap;
@@ -26,6 +25,7 @@ use rhodium_core::{
     tenor::{Frequency, TenorLength, TenorUnit},
     time::{LocalEventTime, TzId, UtcOffset, ZoneRules},
 };
+use rhodium_testkit::{MARKET_ABS, MARKET_REL, assert_close};
 
 struct Dates {
     id: CalendarId,
@@ -150,14 +150,6 @@ fn generated(
     })
 }
 
-fn assert_close(actual: f64, expected: f64) {
-    let scale = actual.abs().max(expected.abs());
-    assert!(
-        (actual - expected).abs() <= (1e-12 * scale).max(1e-14),
-        "{actual} != {expected}"
-    );
-}
-
 proptest! {
     #![proptest_config(Config { cases: 256, rng_seed: RngSeed::Fixed(157), ..Config::default() })]
 
@@ -651,7 +643,12 @@ fn explicit_stub_reference_period() -> Result<()> {
                 .ok_or_else(|| Error::Internal("missing accrual context".into()))?;
             let (fraction, _) = year_fraction(&dc, stub.accrual.start(), stub.accrual.end(), &ctx)?;
             // 2 月 28 日から 3 月 31 日まで 31 日。半年の参照期間は 184 日。
-            assert_close(fraction.value(), 31.0 / (2.0 * 184.0));
+            assert_close!(
+                fraction.value(),
+                31.0 / (2.0 * 184.0),
+                MARKET_REL,
+                MARKET_ABS
+            );
         }
     }
     let regular = explicit("2028-02-29", "2028-08-31", "2028-08-31")?;
@@ -669,7 +666,7 @@ fn explicit_stub_reference_period() -> Result<()> {
         .ok_or_else(|| Error::Internal("missing accrual context".into()))?;
     let p = schedule.periods()[0].accrual();
     let (fraction, _) = year_fraction(&DayCount::Act365L, p.start(), p.end(), &context)?;
-    assert_close(fraction.value(), 184.0 / 366.0);
+    assert_close!(fraction.value(), 184.0 / 366.0, MARKET_REL, MARKET_ABS);
     Ok(())
 }
 
@@ -733,8 +730,13 @@ fn explicit_long_stub_carries_several_reference_periods() -> Result<()> {
         schedule.periods()[0].frame().map(PeriodFrame::regular),
         rule.periods()[0].frame().map(PeriodFrame::regular)
     );
-    assert_close(fractions(&schedule)?, 0.75);
-    assert_close(fractions(&rule)?, fractions(&schedule)?);
+    assert_close!(fractions(&schedule)?, 0.75, MARKET_REL, MARKET_ABS);
+    assert_close!(
+        fractions(&rule)?,
+        fractions(&schedule)?,
+        MARKET_REL,
+        MARKET_ABS
+    );
 
     // 二つの期間に分けた入力も書けるが、別の二つの計算期間になる（固定の係数の和は同じ 3/4）。
     let split = vec![
@@ -749,7 +751,7 @@ fn explicit_long_stub_carries_several_reference_periods() -> Result<()> {
     )?))?;
     assert_eq!(split.periods().len(), 2);
     assert_eq!(split.payments().len(), 1);
-    assert_close(fractions(&split)?, 0.75);
+    assert_close!(fractions(&split)?, 0.75, MARKET_REL, MARKET_ABS);
 
     // 区間を覆わない、連続しない、頻度の無い参照期間の列は拒否する。
     let invalid = [
@@ -831,12 +833,22 @@ fn explicit_reference_accepts_rolled_and_adjusted_regular_periods() -> Result<()
         Ok(dates.generate(&spec)?.0)
     };
     let explicit_rolled = explicit_schedule(rolled, "1M")?;
-    assert_close(fraction(&explicit_rolled, 0)?, fraction(&rule, 1)?);
-    assert_close(fraction(&explicit_rolled, 0)?, 1.0 / 12.0);
+    assert_close!(
+        fraction(&explicit_rolled, 0)?,
+        fraction(&rule, 1)?,
+        MARKET_REL,
+        MARKET_ABS
+    );
+    assert_close!(
+        fraction(&explicit_rolled, 0)?,
+        1.0 / 12.0,
+        MARKET_REL,
+        MARKET_ABS
+    );
     // 確認書が営業日調整後の日付で書いた半年の期間（01-04→07-01 の 178 日。暦の一刻みではない）も、
     // 通常期間として受理し、係数は 1/2 になる。
     let adjusted = explicit_schedule(interval("2027-01-04", "2027-07-01")?, "6M")?;
-    assert_close(fraction(&adjusted, 0)?, 0.5);
+    assert_close!(fraction(&adjusted, 0)?, 0.5, MARKET_REL, MARKET_ABS);
     Ok(())
 }
 
@@ -889,7 +901,7 @@ fn period_frame_rejects_reference_periods_that_are_not_one_step() -> Result<()>
         nine_months.end(),
         &context,
     )?;
-    assert_close(fraction.value(), 0.75);
+    assert_close!(fraction.value(), 0.75, MARKET_REL, MARKET_ABS);
     Ok(())
 }
 
diff --git a/crates/rhodium-core/tests/schedules/observation.rs b/crates/rhodium-core/tests/schedules/observation.rs
index 75b71ab00972866a869c28b00331c7517bd8a0ee..41e53e33abb2d298e21fb9b5ca7775f015a27f9e 100644
--- a/crates/rhodium-core/tests/schedules/observation.rs
+++ b/crates/rhodium-core/tests/schedules/observation.rs
@@ -17,7 +17,9 @@ use rhodium_core::{
     tenor::Frequency,
 };
 
-use super::{Dates, assert_close, generated, interval};
+use rhodium_testkit::{MARKET_ABS, MARKET_REL, assert_close};
+
+use super::{Dates, generated, interval};
 
 impl Dates {
     fn once(&self, start: &str, end: &str) -> Result<Schedule> {
@@ -63,7 +65,7 @@ fn check_daily(schedule: &ObservationSchedule, days: &[(&str, f64)]) -> Result<(
     assert_eq!(schedule.fixings().len(), days.len());
     for (&(point, weight), &(date, days)) in schedule.fixings().iter().zip(days) {
         assert_eq!(point, ObservationPoint::Date(date.parse()?));
-        assert_close(weight, days / 365.0);
+        assert_close!(weight, days / 365.0, MARKET_REL, MARKET_ABS);
     }
     Ok(())
 }
@@ -112,7 +114,12 @@ fn rfr_methods_move_fixing_dates_and_weights_independently() -> Result<()> {
     ] {
         let observations = observe(&dates, period, &compound(&dates, method)?)?;
         check_daily(&observations, &weights)?;
-        assert_close(observations.annualization().value(), 5.0 / 365.0);
+        assert_close!(
+            observations.annualization().value(),
+            5.0 / 365.0,
+            MARKET_REL,
+            MARKET_ABS
+        );
     }
     let weekend = dates.once("2027-01-08", "2027-01-11")?;
     let period = &weekend.periods()[0];
@@ -125,8 +132,13 @@ fn rfr_methods_move_fixing_dates_and_weights_independently() -> Result<()> {
         &AccrualContext::default(),
     )?;
     // 観測は木曜から金曜までの 1 日、利息計算は金曜から月曜までの 3 日。
-    assert_close(observations.annualization().value(), 1.0 / 365.0);
-    assert_close(coupon.value(), 3.0 / 365.0);
+    assert_close!(
+        observations.annualization().value(),
+        1.0 / 365.0,
+        MARKET_REL,
+        MARKET_ABS
+    );
+    assert_close!(coupon.value(), 3.0 / 365.0, MARKET_REL, MARKET_ABS);
     assert!(matches!(
         observe(&dates, period, &compound(&dates, Lockout(lag))?),
         Err(Error::InvalidData(_))
@@ -147,7 +159,12 @@ fn rfr_period_starting_on_holiday() -> Result<()> {
         interval("2027-01-11", "2027-01-12")?,
     ];
     assert_eq!(observations.coverage(), &coverage);
-    assert_close(observations.annualization().value(), 3.0 / 365.0);
+    assert_close!(
+        observations.annualization().value(),
+        3.0 / 365.0,
+        MARKET_REL,
+        MARKET_ABS
+    );
     let lag = NonZeroU8::MIN;
     let lookback = observe(
         &dates,
@@ -206,9 +223,24 @@ fn observation_weights_use_the_passed_snapshot_and_day_count() -> Result<()> {
         &updated.calendars,
     )?;
     // BUS/252 ではどちらの区間も営業日が 1 日ずつ。
-    assert_close(observations.fixings()[0].1, 1.0 / 252.0);
-    assert_close(observations.fixings()[1].1, 1.0 / 252.0);
-    assert_close(observations.annualization().value(), 2.0 / 252.0);
+    assert_close!(
+        observations.fixings()[0].1,
+        1.0 / 252.0,
+        MARKET_REL,
+        MARKET_ABS
+    );
+    assert_close!(
+        observations.fixings()[1].1,
+        1.0 / 252.0,
+        MARKET_REL,
+        MARKET_ABS
+    );
+    assert_close!(
+        observations.annualization().value(),
+        2.0 / 252.0,
+        MARKET_REL,
+        MARKET_ABS
+    );
     assert_eq!(observations.annualization().basis(), &basis);
     Ok(())
 }
@@ -248,7 +280,7 @@ fn calendar_days_months_and_index_endpoints_follow_the_rule() -> Result<()> {
             *point,
             ObservationPoint::Month(YearMonth::new(2027, month)?)
         );
-        assert_close(*weight, days / 365.0);
+        assert_close!(*weight, days / 365.0, MARKET_REL, MARKET_ABS);
     }
     assert_eq!(monthly.fixings().len(), 3);
     month_rule = ObservationRule::new(
@@ -259,7 +291,7 @@ fn calendar_days_months_and_index_endpoints_follow_the_rule() -> Result<()> {
     )?;
     let equal = observe(&dates, &months.periods()[0], &month_rule)?;
     for (_, weight) in equal.fixings() {
-        assert_close(*weight, 1.0);
+        assert_close!(*weight, 1.0, MARKET_REL, MARKET_ABS);
     }
     // 指数比はシフト後の木曜・金曜を使い、年率換算には 1 日分の係数を使う。
     let shifted = ObservationDays::Endpoints(EndpointDates::Shifted {
@@ -284,7 +316,12 @@ fn calendar_days_months_and_index_endpoints_follow_the_rule() -> Result<()> {
             ObservationPoint::Date("2027-01-08".parse()?),
         ]
     );
-    assert_close(endpoints.annualization().value(), 1.0 / 365.0);
+    assert_close!(
+        endpoints.annualization().value(),
+        1.0 / 365.0,
+        MARKET_REL,
+        MARKET_ABS
+    );
     Ok(())
 }
 
@@ -372,7 +409,7 @@ fn all_day_lockout_and_accrual_endpoints_follow_the_rule() -> Result<()> {
     .collect::<Result<Vec<_>>>()?;
     assert_eq!(fixings, expected);
     for &(_, weight) in observed.fixings() {
-        assert_close(weight, 1.0);
+        assert_close!(weight, 1.0, MARKET_REL, MARKET_ABS);
     }
     let covered = [
         "2027-01-04",
@@ -388,7 +425,12 @@ fn all_day_lockout_and_accrual_endpoints_follow_the_rule() -> Result<()> {
     })
     .collect::<Result<Vec<_>>>()?;
     assert_eq!(observed.coverage(), covered.as_slice());
-    assert_close(observed.annualization().value(), 5.0 / 365.0);
+    assert_close!(
+        observed.annualization().value(),
+        5.0 / 365.0,
+        MARKET_REL,
+        MARKET_ABS
+    );
 
     // 締めの日が残らない短い期間（2 行で Lockout(2)）は適用できない。
     let short = dates.once("2027-01-07", "2027-01-09")?;
@@ -424,9 +466,14 @@ fn all_day_lockout_and_accrual_endpoints_follow_the_rule() -> Result<()> {
         ]
     );
     for &(_, weight) in endpoints.fixings() {
-        assert_close(weight, 1.0);
+        assert_close!(weight, 1.0, MARKET_REL, MARKET_ABS);
     }
-    assert_close(endpoints.annualization().value(), 3.0 / 365.0);
+    assert_close!(
+        endpoints.annualization().value(),
+        3.0 / 365.0,
+        MARKET_REL,
+        MARKET_ABS
+    );
     Ok(())
 }
 
@@ -443,9 +490,11 @@ fn factor_excess_keeps_the_digits_with_arithmetic_only() -> Result<()> {
         CompoundingForm::SimpleFactor.factor_excess(0.01, 0.25, 0.002)?,
         (0.01 + 0.002) * 0.25
     );
-    assert_close(
+    assert_close!(
         CompoundingForm::PowerFactor.factor_excess(0.03, 0.25, 0.0)?,
         CompoundingForm::PowerFactor.factor(0.03, 0.25, 0.0)? - 1.0,
+        MARKET_REL,
+        MARKET_ABS
     );
     // 因子が正にならない値と非有限の値は、factor と同じく Numerical。
     for (form, rate, weight, spread) in [
@@ -479,11 +528,13 @@ fn log_factor_keeps_the_digits_that_one_plus_x_rounds_away() -> Result<()> {
 
     // スプレッドは SimpleFactor の (r + spread)·w に入り、PowerFactor は w·ln(1 + r) になる。
     let spread = CompoundingForm::SimpleFactor.log_factor(0.01, 0.5, 0.002)?;
-    assert_close(spread, (0.012_f64 * 0.5).ln_1p());
+    assert_close!(spread, (0.012_f64 * 0.5).ln_1p(), MARKET_REL, MARKET_ABS);
     let power = CompoundingForm::PowerFactor.log_factor(0.03, 0.25, 0.0)?;
-    assert_close(
+    assert_close!(
         power,
         CompoundingForm::PowerFactor.factor(0.03, 0.25, 0.0)?.ln(),
+        MARKET_REL,
+        MARKET_ABS
     );
 
     // 因子が正にならない値と非有限の値は、factor と同じく Numerical。
diff --git a/crates/rhodium-fixtures/Cargo.toml b/crates/rhodium-fixtures/Cargo.toml
new file mode 100644
index 0000000000000000000000000000000000000000..9781f7804b41eb79f69ec708b0f58c64e0e290d1
--- /dev/null
+++ b/crates/rhodium-fixtures/Cargo.toml
@@ -0,0 +1,20 @@
+[package]
+name = "rhodium-fixtures"
+description = "複数の crate のテストとベンチマークが共有する市場、参照データ、ワークスペースの fixture を持つ開発専用の crate"
+version.workspace = true
+edition.workspace = true
+rust-version.workspace = true
+license.workspace = true
+publish.workspace = true
+
+[lints]
+workspace = true
+
+[dependencies]
+rhodium-core.workspace = true
+rhodium-market.workspace = true
+rhodium-pricing.workspace = true
+rhodium-products.workspace = true
+serde.workspace = true
+serde_json = "1"
+toml = "1.2"
diff --git a/crates/rhodium-fixtures/quantlib/README.md b/crates/rhodium-fixtures/quantlib/README.md
new file mode 100644
index 0000000000000000000000000000000000000000..dc8bd9a611f54ba89bb66201f7e65755a72e006a
--- /dev/null
+++ b/crates/rhodium-fixtures/quantlib/README.md
@@ -0,0 +1,73 @@
+# QuantLib の OIS とターム指標の曲線
+
+QuantLib 1.42 の Python バインディングで生成した値。
+[QuantLib](https://github.com/lballabio/QuantLib/tree/v1.42) と [QuantLib-SWIG](https://github.com/lballabio/QuantLib-SWIG/tree/v1.42) は Modified BSD License で公開されている。
+生成データとスクリプトはこのリポジトリと同じライセンスで配布する。
+
+この 3 つの JSON は、複数の crate のテストが共有する市場の fixture（`rhodium-fixtures` の `ois`、`steps`、`multicurve`、`xccy`）が読む。
+同じディレクトリの `generate_ois.py` が生成する。
+スワップ PV（`swaps.json`）とボラ（`volatility.json`）は pricing のテストだけが読むので、`crates/rhodium-pricing/tests/fixtures/quantlib/` に置く。
+
+## TONA と SOFR の OIS 曲線
+
+`ois_curves.json` は、QuantLib 1.42 の `OISRateHelper` と
+`PiecewiseLogLinearDiscount` で 1M・6M・1Y・18M・2Y の quote をブートストラップした結果である。
+同じ曲線で割引と翌日物の予想を行う。
+約定日は 2025 年 1 月 7 日、スポットラグは 2 営業日、固定脚は年払いとする。
+両実装で満期から逆算して短い初回期間を作り、期間境界に Modified Following、支払いに Following を適用する。
+JPY は Japan カレンダーと ACT/365F で支払ラグは 0、USD は SOFR の fixing カレンダーと ACT/360 で支払ラグは 2 営業日とする。
+ピラーは QuantLib の `LastRelevantDate`（各商品の最終支払日）とする。
+評価、観測と支払いを UTC 00:00 にそろえ、曲線の時間軸は両実装で ACT/365F とする。
+
+```sh
+uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-fixtures/quantlib/generate_ois.py > crates/rhodium-fixtures/quantlib/ois_curves.json
+```
+
+Rust では各ピラーの DF と、評価日から最後のピラーまでの 5 つの区間に一つずつ置いた照合点（2025-01-24、2025-04-15、2025-08-12、2026-04-15、2026-10-15）の DF を、相対 1e-12、絶対 1e-14 の下限で照合する。
+両実装の停止判定と浮動小数点の演算順が異なるため、この許容誤差を使う。
+quote を固定金利にした校正商品を取引と同じ準備で評価し、PV が年金現価 × 2 × 残差の許容誤差（1e-14）以内で 0 になることも確かめる。
+
+同じ曲線で、QuantLib の `OvernightIndexedSwap`（日程は年払い、満期から逆算、期間境界に Modified Following、支払いに Following と各通貨の支払ラグ）を二つ評価し、PV と脚別の PV を `swaps` に置く。
+
+| 名前 | 開始 | テナー | 固定金利 | スプレッド | 固定脚 | 元本 |
+|---|---|---|---|---|---|---|
+| `spot` | 約定日（評価日）の 2 営業日後 | 18M | 18M の quote + 25bp | 0 | 受け | JPY 1 億、USD 100 万 |
+| `seasoned` | 2024-12-05 | 1Y | JPY 0.35%、USD 4.20% | 10bp | 払い | JPY 1 億、USD 100 万 |
+
+開始済みのスワップは、2024-12-05 から評価日の前日までの各営業日の fixing（`fixings`、JPY は 0.23% から、USD は 4.55% から 1bp ずつ 3 日周期で変える値）を確定値として読み、評価日の観測は曲線から予想する。
+rhodium の参照データでは fixing の公表を翌営業日の UTC 00:00 とするので、評価日の前日の値は評価の瞬間ちょうどに公表済みになり、評価日の観測期間は評価の瞬間に始まる。
+QuantLib はスプレッドを日次の複利の外で足す（rhodium の `Exclusive`）。生成スクリプトは、スプレッドの有無の変動脚の差が Σ 元本 × スプレッド × 計算期間 × DF に一致することを確かめてから出力する。
+Rust は QuantLib のノードの DF を θ に入れた市場で、同じ取引を `Swap::from_convention`（スポット起点は開始日を省いて慣行のスポットラグで求める）と同じ fixing で作り、PV、変動脚（固定金利 0 の PV）、固定脚（差）を照合する。
+許容誤差は相対 1e-12、絶対の下限は通貨の最小単位の 1e-6 倍とし、相対の尺度には両脚の絶対値の和と、変動脚の二つの項（N·P(s)/P(e) と N）の大きさとして元本の 2 倍を含める。
+
+## 会合日のステップと年末のターンを持つ TONA
+
+`ois_turn.json` は同じ TONA の quote を、年末のターンを含む `PiecewiseLogLinearDiscount` に校正した結果である。
+Japan カレンダーで期末の営業日から翌営業日までの日数係数を求め、50bp を掛けた値を ln DF の跳びにする。
+QuantLib のジャンプは指定日より後に効くため、QuantLib には発効の前日を指定する。
+rhodium は最初のピラーで会合のステップを切り替え、二つ目のピラーで長期の節点につなぐ。
+両実装のピラーと中間日の DF を相対 1e-12、絶対 1e-14 で比べ、JPY と USD の OIS quote を再評価する。
+
+```sh
+uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-fixtures/quantlib/generate_ois.py --turn > crates/rhodium-fixtures/quantlib/ois_turn.json
+```
+
+## TONA と TIBOR 6M・3M、SOFR と USD ターム指標の複数曲線
+
+`multicurve.json` は、QuantLib 1.42 の `OISRateHelper` で OIS 曲線を先に校正し、`OvernightIborBasisSwapRateHelper` で別のターム予想曲線を校正した結果である。
+OIS とベーシスには 6M、1Y、18M、2Y を使う。
+ベーシス商品は 2 営業日後に開始し、半年ごとにラグなしで支払う。
+翌日物脚に複利の外で spread を足して払い、ターム脚を受ける。
+ターム指標は JPY では QuantLib の TIBOR 6M と TIBOR 3M、USD では SOFR のカレンダーを使うテスト用の 6M 指標で、予想期間の端点は指標の fixing 規則で決める。
+両脚の支払いの頻度はターム指標のテナー（TIBOR 3M は四半期ごと）で、出力の `term_tenor` に記録する。
+配列の順は TIBOR 6M、USD のターム指標、TIBOR 3M（三つ目は M1 で足した）。
+OIS 曲線を両脚の割引と翌日物の予想に使う。
+
+JPY は割引とターム予想のピラーがそろっており、rhodium ではターム予想を TONA とベーシスの和で表す。
+共通のピラーの間ではどちらの実装も対数 DF が線形になる。
+USD の OIS は 2 営業日の支払ラグがあるため、ターム予想は独立した節点曲線で表す。
+Rust のテストはピラーと補間点の割引係数および予想の擬似 DF を相対 1e-12、絶対 1e-14 で比較し、元本 1 の校正商品を quote どおりに評価した PV を残差の許容誤差から求めた幅で照合する。
+
+```sh
+uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-fixtures/quantlib/generate_ois.py --multicurve > crates/rhodium-fixtures/quantlib/multicurve.json
+```
diff --git a/crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py b/crates/rhodium-fixtures/quantlib/generate_ois.py
similarity index 100%
rename from crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py
rename to crates/rhodium-fixtures/quantlib/generate_ois.py
diff --git a/crates/rhodium-pricing/tests/fixtures/quantlib/multicurve.json b/crates/rhodium-fixtures/quantlib/multicurve.json
similarity index 100%
rename from crates/rhodium-pricing/tests/fixtures/quantlib/multicurve.json
rename to crates/rhodium-fixtures/quantlib/multicurve.json
diff --git a/crates/rhodium-pricing/tests/fixtures/quantlib/ois_curves.json b/crates/rhodium-fixtures/quantlib/ois_curves.json
similarity index 100%
rename from crates/rhodium-pricing/tests/fixtures/quantlib/ois_curves.json
rename to crates/rhodium-fixtures/quantlib/ois_curves.json
diff --git a/crates/rhodium-pricing/tests/fixtures/quantlib/ois_turn.json b/crates/rhodium-fixtures/quantlib/ois_turn.json
similarity index 100%
rename from crates/rhodium-pricing/tests/fixtures/quantlib/ois_turn.json
rename to crates/rhodium-fixtures/quantlib/ois_turn.json
diff --git a/crates/rhodium/tests/support/desk.rs b/crates/rhodium-fixtures/src/desk.rs
similarity index 97%
rename from crates/rhodium/tests/support/desk.rs
rename to crates/rhodium-fixtures/src/desk.rs
index 0bef0c61225a03492fd15b38dd8bfae9b785a01c..b27e64cadfb9a93a2a4c2ddecd7b3808f9d9976e 100644
--- a/crates/rhodium/tests/support/desk.rs
+++ b/crates/rhodium-fixtures/src/desk.rs
@@ -5,8 +5,9 @@
 //! 取引は慣行どおりの固定対 RFR の OIS で、通貨、テナー、受払、スプレッド、
 //! スポット起点・開始済み・先日付を ID の順に均等に混ぜる。値はすべて決まった規則で作る。
 
-use std::{fmt::Write as _, fs, io, path::Path};
+use std::{fmt::Write as _, fs, path::Path};
 
+use crate::files::io_error;
 use rhodium_core::{
     Error, Result,
     date::{Date, Weekday},
@@ -111,11 +112,6 @@ pub fn write(dir: &Path, trades: usize) -> Result<()> {
     Ok(())
 }
 
-#[expect(clippy::needless_pass_by_value, reason = "map_err の変換関数に使う")]
-fn io_error(error: io::Error) -> Error {
-    Error::Internal(format!("desk workspace: {error}"))
-}
-
 fn format_error(error: std::fmt::Error) -> Error {
     Error::Internal(format!("desk workspace: {error}"))
 }
@@ -133,15 +129,11 @@ fn add_months(date: Date, months: i32) -> Result<Date> {
 }
 
 // 週末なら次の月曜に進める（週末だけのカレンダーの Following）。
-#[expect(
-    clippy::wildcard_enum_match_arm,
-    reason = "#196 で各 variant を網羅する"
-)]
 fn next_weekday(date: Date) -> Result<Date> {
     match date.weekday() {
         Weekday::Sat => date.checked_add_days(2),
         Weekday::Sun => date.checked_add_days(1),
-        _ => Ok(date),
+        Weekday::Mon | Weekday::Tue | Weekday::Wed | Weekday::Thu | Weekday::Fri => Ok(date),
     }
 }
 
diff --git a/crates/rhodium-fixtures/src/files.rs b/crates/rhodium-fixtures/src/files.rs
new file mode 100644
index 0000000000000000000000000000000000000000..aa0207aa30d5b7b40be0f884ea63a7df135a6d88
--- /dev/null
+++ b/crates/rhodium-fixtures/src/files.rs
@@ -0,0 +1,51 @@
+//! テストがリポジトリのサンプルのワークスペースを一時ディレクトリに複製し、書き換える補助。
+
+use std::{
+    fs,
+    io::{self, Write as _},
+    path::{Path, PathBuf},
+};
+
+use rhodium_core::{Error, Result};
+
+/// リポジトリの `examples/<name>` のディレクトリ（`m0-ois`、`m1-rates` など）。
+#[must_use]
+pub fn example(name: &str) -> PathBuf {
+    #[expect(
+        clippy::disallowed_macros,
+        reason = "リポジトリのサンプルの場所を fixture の crate のソースの位置から決める"
+    )]
+    let root = env!("CARGO_MANIFEST_DIR");
+    Path::new(root).join("../../examples").join(name)
+}
+
+/// テストの入出力の失敗を、テストが返す `rhodium_core::Error` にする。
+#[expect(clippy::needless_pass_by_value, reason = "map_err に渡す変換関数")]
+#[must_use]
+pub fn io_error(error: io::Error) -> Error {
+    Error::Internal(format!("test files: {error}"))
+}
+
+/// `from` のディレクトリの中身を、既存のディレクトリ `to` の下に再帰的に複製する。
+pub fn copy_dir(from: &Path, to: &Path) -> Result<()> {
+    for entry in fs::read_dir(from).map_err(io_error)? {
+        let entry = entry.map_err(io_error)?;
+        let target = to.join(entry.file_name());
+        if entry.file_type().map_err(io_error)?.is_dir() {
+            fs::create_dir(&target).map_err(io_error)?;
+            copy_dir(&entry.path(), &target)?;
+        } else {
+            fs::copy(entry.path(), target).map_err(io_error)?;
+        }
+    }
+    Ok(())
+}
+
+/// 既存のファイル `path` の末尾に `suffix` を書き足す。
+pub fn append(path: &Path, suffix: &str) -> Result<()> {
+    fs::OpenOptions::new()
+        .append(true)
+        .open(path)
+        .and_then(|mut file| file.write_all(suffix.as_bytes()))
+        .map_err(io_error)
+}
diff --git a/crates/rhodium-pricing/tests/fixtures/futures.rs b/crates/rhodium-fixtures/src/futures.rs
similarity index 80%
rename from crates/rhodium-pricing/tests/fixtures/futures.rs
rename to crates/rhodium-fixtures/src/futures.rs
index 8946183d481f121f62d88fdd9818dfd500383ad6..1f340ef4661fb14b69aa35e8b92eae17c472a53c 100644
--- a/crates/rhodium-pricing/tests/fixtures/futures.rs
+++ b/crates/rhodium-fixtures/src/futures.rs
@@ -41,19 +41,28 @@ use rhodium_pricing::{
 };
 use rhodium_products::standard::InstrumentTemplate;
 
-pub(crate) const START: &str = "2026-01-05";
-pub(crate) const AS_OF: &str = "2026-01-05T00:00:00Z";
-pub(crate) const FINAL: &str = "2026-01-07";
-pub(crate) const PILLAR: &str = "2026-05-08T16:00:00Z";
-pub(crate) const CLEARING: &str = "USD-FUT-JAN";
-pub(crate) const MULTIPLIER: f64 = 1_000.0;
-pub(crate) const RATE: f64 = 0.04;
+/// 先物の取引の約定日と評価日。
+pub const START: &str = "2026-01-05";
+/// 評価の瞬間（評価日の UTC 0 時）。
+pub const AS_OF: &str = "2026-01-05T00:00:00Z";
+/// 2026 年 1 月限の最終清算日。原資産のターム金利の fixing の日でもある。
+pub const FINAL: &str = "2026-01-07";
+/// USD 曲線の唯一のピラー。
+pub const PILLAR: &str = "2026-05-08T16:00:00Z";
+/// 限月の清算価格の指標と、上場先物の銘柄の ID。
+pub const CLEARING: &str = "USD-FUT-JAN";
+/// 先物の乗数。
+pub const MULTIPLIER: f64 = 1_000.0;
+/// [`market`] の曲線の連続複利のゼロレート。
+pub const RATE: f64 = 0.04;
 
-pub(crate) fn at(day: &str, time: &str) -> Result<Instant> {
+/// 日付と UTC の時刻（`HH:MM`）を瞬間にする。
+pub fn at(day: &str, time: &str) -> Result<Instant> {
     format!("{day}T{time}:00Z").parse()
 }
 
-pub(crate) fn refdata(delay: DateOffset, settlement_time: &str) -> Result<RefData> {
+/// 値洗いの受渡の遅れ `delay` と受渡の時刻 `settlement_time` の先物の慣行を持つ参照データ。
+pub fn refdata(delay: DateOffset, settlement_time: &str) -> Result<RefData> {
     RefData::new(reference_parts(delay, settlement_time)?)
 }
 
@@ -61,7 +70,8 @@ pub(crate) fn refdata(delay: DateOffset, settlement_time: &str) -> Result<RefDat
     clippy::too_many_lines,
     reason = "暦、原指標、清算と限月のテスト用参照データをここにまとめる"
 )]
-pub(crate) fn reference_parts(delay: DateOffset, settlement_time: &str) -> Result<RefDataParts> {
+/// [`refdata`] の参照先を解決する前の部品。テストが慣行を書き換えてから解決するときに使う。
+pub fn reference_parts(delay: DateOffset, settlement_time: &str) -> Result<RefDataParts> {
     let calendar: CalendarId = "USD".parse()?;
     let zone: TzId = "Etc/UTC".parse()?;
     let time = "16:00:00".parse()?;
@@ -174,7 +184,8 @@ pub(crate) fn reference_parts(delay: DateOffset, settlement_time: &str) -> Resul
     Ok(parts)
 }
 
-pub(crate) fn config(extra: &str) -> Result<PricingConfig> {
+/// USD の無担保の割引とターム指標の予想に USD 曲線を使う価格設定。`extra` は TOML の末尾に足す。
+pub fn config(extra: &str) -> Result<PricingConfig> {
     toml::from_str(&format!(
         r#"
 version = 1
@@ -191,7 +202,8 @@ TERM = "USD"
     .map_err(|error| Error::InvalidData(error.to_string()))
 }
 
-pub(crate) fn curve() -> Result<CurveSpec> {
+/// 一つのピラーを持つ USD の割引と予想の曲線。
+pub fn curve() -> Result<CurveSpec> {
     CurveSpec::new(
         "USD".parse()?,
         Currency::USD,
@@ -201,7 +213,8 @@ pub(crate) fn curve() -> Result<CurveSpec> {
     )
 }
 
-pub(crate) fn market(as_of: Instant, volatility: f64) -> Result<Market<f64>> {
+/// ゼロレート [`RATE`] の USD 曲線と、平均回帰 0、ボラ `volatility` の USD の LGM を持つ市場。
+pub fn market(as_of: Instant, volatility: f64) -> Result<Market<f64>> {
     let structure = Arc::new(MarketStructure::with_fx_and_lgm(
         as_of,
         [curve()?],
@@ -223,7 +236,7 @@ pub(crate) fn market(as_of: Instant, volatility: f64) -> Result<Market<f64>> {
 }
 
 /// 2026 年 1 月限の金利先物の校正テンプレート。
-pub(crate) fn futures_template() -> Result<InstrumentTemplate> {
+pub fn futures_template() -> Result<InstrumentTemplate> {
     Ok(InstrumentTemplate::RateFuture {
         convention: "USD-RATE-FUTURES".parse()?,
         delivery_month: "2026-01".parse()?,
@@ -232,7 +245,7 @@ pub(crate) fn futures_template() -> Result<InstrumentTemplate> {
 
 /// USD 曲線の一つのピラーを先物の価格で校正し、凸性の補正に使う LGM の平均回帰とボラを
 /// 直接マークにする市場の定義。
-pub(crate) fn futures_definition() -> Result<MarketDefinition> {
+pub fn futures_definition() -> Result<MarketDefinition> {
     let binding = InstrumentBinding::new(
         "USD.FUT.2026-01".parse()?,
         QuoteKind::FuturesPrice,
diff --git a/crates/rhodium-fixtures/src/lib.rs b/crates/rhodium-fixtures/src/lib.rs
new file mode 100644
index 0000000000000000000000000000000000000000..d5366b4f3efac64017ed98d68e5aba9887ea8fc8
--- /dev/null
+++ b/crates/rhodium-fixtures/src/lib.rs
@@ -0,0 +1,104 @@
+//! 複数の crate のテスト、ベンチマーク、例が共有する fixture を持つ開発専用の crate。
+//!
+//! 閉包と数値だけを受け取る比較の補助と許容誤差の定数は `rhodium-testkit` に置き、
+//! この crate は rhodium の型で書いた市場データ、参照データ、商品、ワークスペースの組み立てを持つ。
+//! 依存は core、products、market、pricing に限り、pricing、risk、ファサードの `[dev-dependencies]`
+//! から使う。risk やファサードの型を扱う補助は、それぞれの crate の `tests/` に置く。
+//!
+//! `quantlib/` の JSON は `QuantLib` で生成した外部の fixture で、出典と生成の手順は
+//! `quantlib/README.md` にある。テストの実行時に外部のライブラリは呼ばない。
+
+#![forbid(unsafe_code)]
+#![expect(
+    clippy::missing_errors_doc,
+    reason = "fixture の関数は、固定の ID と日付の解析や参照データの構築の失敗をそのまま返すだけで、\
+              関数ごとに書き分ける失敗の条件を持たない"
+)]
+
+pub mod desk;
+pub mod files;
+pub mod futures;
+pub mod multicurve;
+pub mod ois;
+pub mod steps;
+pub mod xccy;
+
+use rhodium_core::{Error, Result, date::Date};
+use rhodium_market::{
+    market::Market,
+    quote::{QuoteId, QuoteSet},
+};
+use rhodium_pricing::{
+    calibration::{CalibratedMarket, GroupOutcome},
+    context::Context,
+    prepare::{PrepareInput, Prepared, prepare},
+    resolve::CollateralKey,
+};
+use rhodium_products::standard::StandardProduct;
+
+/// 校正商品と取引の既定の担保。価格設定の無担保の割引の行で割り引く。
+pub const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
+
+/// `quotes` の `id` の値だけを `bump` 動かした quote の集合。中心差分の再校正の入力に使う。
+pub fn bumped(quotes: &QuoteSet, id: &QuoteId, bump: f64) -> Result<QuoteSet> {
+    QuoteSet::new(
+        quotes
+            .iter()
+            .map(|(key, value)| (key.clone(), if key == id { value + bump } else { value })),
+    )
+}
+
+/// 校正のすべてのグループが解けたことを確かめる。
+///
+/// 解けないグループがあれば、そのグループの記録を付けた `Calibration` のエラーを返す。
+/// 最小二乗のグループは停留条件で止まり、残差が残るので、残差の大きさは確かめない。
+pub fn solved(market: &CalibratedMarket) -> Result<()> {
+    market
+        .groups()
+        .iter()
+        .find(|group| !matches!(group.outcome(), GroupOutcome::Solved(_)))
+        .map_or(Ok(()), |group| {
+            Err(Error::Calibration(format!("{group:?}")))
+        })
+}
+
+/// [`solved`] に加えて、各グループの残差のノルムが `tolerance` 以下であることを確かめる。
+/// quote を再現する正方系だけの市場に使う。
+///
+/// # Panics
+///
+/// 解けたグループの残差のノルムが `tolerance` を超える場合。
+pub fn solved_within(market: &CalibratedMarket, tolerance: f64) -> Result<()> {
+    solved(market)?;
+    for group in market.groups() {
+        if let GroupOutcome::Solved(solution) = group.outcome() {
+            assert!(
+                solution.residual_norm() <= tolerance,
+                "{:?}: residual norm {} exceeds {tolerance}",
+                group.blocks(),
+                solution.residual_norm()
+            );
+        }
+    }
+    Ok(())
+}
+
+/// 標準商品を `trade_date` の約定として、`collateral` の割引と `market` の構造で準備する。
+pub fn prepare_standard(
+    product: &StandardProduct,
+    trade_date: Date,
+    collateral: &CollateralKey,
+    market: &Market<f64>,
+    ctx: &Context<'_>,
+) -> Result<Prepared> {
+    prepare(
+        PrepareInput::Standard {
+            product,
+            trade_date,
+        },
+        collateral,
+        &[],
+        market.structure(),
+        ctx,
+    )
+}
diff --git a/crates/rhodium-pricing/tests/fixtures/multicurve.rs b/crates/rhodium-fixtures/src/multicurve.rs
similarity index 89%
rename from crates/rhodium-pricing/tests/fixtures/multicurve.rs
rename to crates/rhodium-fixtures/src/multicurve.rs
index b53eb06eb050b95cb8916e25e86c2a3949c4bf5d..7830431694ed9291dcaa7e458c88646105b23cee 100644
--- a/crates/rhodium-pricing/tests/fixtures/multicurve.rs
+++ b/crates/rhodium-fixtures/src/multicurve.rs
@@ -23,15 +23,16 @@ use rhodium_pricing::{
     definition::{
         CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
     },
-    prepare::{PrepareInput, Prepared, prepare},
+    prepare::Prepared,
     resolve::CollateralKey,
 };
 use rhodium_products::standard::InstrumentTemplate;
 use serde::Deserialize;
 
-use super::ois::{Case, Inputs, Node, Sample, cases, fixture_quotes, midnight};
-
-const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
+use crate::{
+    COLLATERAL,
+    ois::{Case, Inputs, Node, Sample, cases, fixture_quotes, midnight},
+};
 
 #[derive(Deserialize)]
 struct QuantLibBasis {
@@ -58,15 +59,15 @@ pub struct MulticurveCase {
 }
 
 /// JPY の TONA/TIBOR 6M と USD の SOFR/テスト用ターム指標の入力を読む。
-pub(crate) fn multicurve_cases() -> Result<[MulticurveCase; 2]> {
+pub fn multicurve_cases() -> Result<[MulticurveCase; 2]> {
     let [jpy, usd, _] = all_multicurve_cases()?;
     Ok([jpy, usd])
 }
 
 /// TIBOR 6M、USD のターム指標、TIBOR 3M の三つのベーシスの入力を読む。
-pub(crate) fn all_multicurve_cases() -> Result<[MulticurveCase; 3]> {
+pub fn all_multicurve_cases() -> Result<[MulticurveCase; 3]> {
     let [jpy, usd, jpy_3m]: [QuantLibBasis; 3] =
-        serde_json::from_str(include_str!("quantlib/multicurve.json"))
+        serde_json::from_str(include_str!("../quantlib/multicurve.json"))
             .map_err(|error| Error::InvalidData(error.to_string()))?;
     let [tona, sofr] = cases()?;
     Ok([
@@ -93,12 +94,12 @@ impl QuantLibBasis {
 
 impl MulticurveCase {
     /// OIS にターム指標とベーシス商品の慣行を加える。
-    pub(crate) fn inputs(&self) -> Result<Inputs> {
+    pub fn inputs(&self) -> Result<Inputs> {
         self.inputs_with("")
     }
 
     /// [`Self::inputs`] の価格設定に、割引の行や近似の設定を加える。
-    pub(crate) fn inputs_with(&self, extra_config: &str) -> Result<Inputs> {
+    pub fn inputs_with(&self, extra_config: &str) -> Result<Inputs> {
         Inputs::with_forecasts(
             &[&self.overnight],
             self.overnight.as_of,
@@ -175,26 +176,29 @@ impl MulticurveCase {
     }
 
     /// 固定対ターム金利の IRS の慣行の ID。
-    pub(crate) fn irs_convention(&self) -> String {
+    #[must_use]
+    pub fn irs_convention(&self) -> String {
         format!("{}-IRS", self.term_index)
     }
 
     /// スポット起点の預金の慣行の ID。
-    pub(crate) fn deposit_convention(&self) -> String {
+    #[must_use]
+    pub fn deposit_convention(&self) -> String {
         format!("{}-DEPO", self.overnight.index)
     }
 
     /// ターム予想曲線の校正対象。
-    pub(crate) fn term_block(&self) -> Result<BlockId> {
+    pub fn term_block(&self) -> Result<BlockId> {
         BlockId::new(self.term_index.parse()?, "basis")
     }
 
     /// ベーシス商品の quote ID。
-    pub(crate) fn quote(&self, node: &Node) -> Result<QuoteId> {
+    pub fn quote(&self, node: &Node) -> Result<QuoteId> {
         format!("{}.BASIS.{}", self.term_index, node.tenor).parse()
     }
 
-    pub(crate) fn template(&self, node: &Node) -> Result<InstrumentTemplate> {
+    /// ピラーに対応するベーシススワップの校正テンプレート。
+    pub fn template(&self, node: &Node) -> Result<InstrumentTemplate> {
         Ok(InstrumentTemplate::BasisSwap {
             convention: self.term_index.parse()?,
             tenor: node.tenor,
@@ -202,12 +206,12 @@ impl MulticurveCase {
     }
 
     /// 共通の OIS 曲線と、それを上流に持つターム予想曲線を定義する。
-    pub(crate) fn definition(&self) -> Result<MarketDefinition> {
+    pub fn definition(&self) -> Result<MarketDefinition> {
         self.definition_with([])
     }
 
     /// [`Self::definition`] に、校正しない曲線（担保条件の違う割引曲線など）を加える。
-    pub(crate) fn definition_with(
+    pub fn definition_with(
         &self,
         extra: impl IntoIterator<Item = CurveSpec>,
     ) -> Result<MarketDefinition> {
@@ -215,7 +219,7 @@ impl MulticurveCase {
     }
 
     /// 共通の OIS 曲線に、ターム予想曲線のブロックを当てる任意の binding を組み合わせる。
-    pub(crate) fn definition_with_term_bindings(
+    pub fn definition_with_term_bindings(
         &self,
         term_bindings: impl IntoIterator<Item = InstrumentBinding>,
     ) -> Result<MarketDefinition> {
@@ -237,7 +241,7 @@ impl MulticurveCase {
     }
 
     /// ターム予想曲線。JPY は OIS 曲線とベーシスの和、USD は独立した節点の曲線。
-    pub(crate) fn term_curve(&self) -> Result<CurveSpec> {
+    pub fn term_curve(&self) -> Result<CurveSpec> {
         let term_nodes = CurveExpr::Nodes {
             part: "basis".into(),
             pillars: self
@@ -264,7 +268,7 @@ impl MulticurveCase {
 
     /// [`Self::definition_with`] で、`collateral` が返すベーシスの行（ノードの番号）だけ
     /// 校正商品の担保を変える。
-    pub(crate) fn definition_with_collateral(
+    pub fn definition_with_collateral(
         &self,
         extra: impl IntoIterator<Item = CurveSpec>,
         collateral: impl Fn(usize) -> Option<CollateralKey>,
@@ -322,7 +326,7 @@ impl MulticurveCase {
     }
 
     /// OIS とベーシスの入力 quote を作る。
-    pub(crate) fn quotes(&self) -> Result<QuoteSet> {
+    pub fn quotes(&self) -> Result<QuoteSet> {
         let term = self
             .nodes
             .iter()
@@ -335,7 +339,7 @@ impl MulticurveCase {
     }
 
     /// 指定のテナーと spread のベーシス取引を準備する。
-    pub(crate) fn prepare_basis(
+    pub fn prepare_basis(
         &self,
         node: &Node,
         spread: f64,
@@ -343,15 +347,6 @@ impl MulticurveCase {
         ctx: &Context<'_>,
     ) -> Result<Prepared> {
         let product = self.template(node)?.build(spread, ctx.refdata)?;
-        prepare(
-            PrepareInput::Standard {
-                product: &product,
-                trade_date: self.overnight.as_of,
-            },
-            &COLLATERAL,
-            &[],
-            market.structure(),
-            ctx,
-        )
+        crate::prepare_standard(&product, self.overnight.as_of, &COLLATERAL, market, ctx)
     }
 }
diff --git a/crates/rhodium-pricing/tests/fixtures/ois.rs b/crates/rhodium-fixtures/src/ois.rs
similarity index 93%
rename from crates/rhodium-pricing/tests/fixtures/ois.rs
rename to crates/rhodium-fixtures/src/ois.rs
index 864630f938b776d58e4e734bcef8b6b0524ea5b6..1829f7576d09bddddf29444dd0bf9c5e35438d40 100644
--- a/crates/rhodium-pricing/tests/fixtures/ois.rs
+++ b/crates/rhodium-fixtures/src/ois.rs
@@ -40,7 +40,7 @@ use rhodium_pricing::{
     definition::{
         CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
     },
-    prepare::{PrepareInput, Prepared, prepare},
+    prepare::Prepared,
     resolve::CollateralKey,
 };
 use rhodium_products::standard::{InstrumentTemplate, PayReceive};
@@ -132,19 +132,20 @@ pub struct Node {
 }
 
 /// TONA と SOFR のテスト入力を読み込む。
-pub(crate) fn cases() -> Result<[Case; 2]> {
-    serde_json::from_str(include_str!("quantlib/ois_curves.json"))
+pub fn cases() -> Result<[Case; 2]> {
+    serde_json::from_str(include_str!("../quantlib/ois_curves.json"))
         .map_err(|error| Error::InvalidData(error.to_string()))
 }
 
 /// 日付を UTC の 00:00 に解決する。
-pub(crate) fn midnight(date: Date) -> Result<Instant> {
+pub fn midnight(date: Date) -> Result<Instant> {
     format!("{date}T00:00:00Z").parse()
 }
 
 impl Case {
     /// 固定脚と観測に使う日数基準。
-    pub(crate) fn basis(&self) -> DayCount {
+    #[must_use]
+    pub fn basis(&self) -> DayCount {
         if self.currency == Currency::JPY {
             DayCount::Act365F
         } else {
@@ -153,12 +154,12 @@ impl Case {
     }
 
     /// 曲線の節点ブロック。
-    pub(crate) fn block(&self) -> Result<BlockId> {
+    pub fn block(&self) -> Result<BlockId> {
         Ok(BlockId::nodes(self.index.parse()?))
     }
 
     /// 割引と予想に使える曲線。
-    pub(crate) fn curve(&self) -> Result<CurveSpec> {
+    pub fn curve(&self) -> Result<CurveSpec> {
         CurveSpec::new(
             self.index.parse()?,
             self.currency,
@@ -172,7 +173,7 @@ impl Case {
     }
 
     /// ピラーに対応する慣行ベースの OIS。
-    pub(crate) fn template(&self, node: &Node) -> Result<InstrumentTemplate> {
+    pub fn template(&self, node: &Node) -> Result<InstrumentTemplate> {
         Ok(InstrumentTemplate::Ois {
             convention: self.index.parse()?,
             tenor: node.tenor,
@@ -180,12 +181,12 @@ impl Case {
     }
 
     /// ピラーの quote ID。
-    pub(crate) fn quote(&self, node: &Node) -> Result<QuoteId> {
+    pub fn quote(&self, node: &Node) -> Result<QuoteId> {
         format!("{}.OIS.{}", self.index, node.tenor).parse()
     }
 
     /// ピラーの校正条件。
-    pub(crate) fn binding(
+    pub fn binding(
         &self,
         node: &Node,
         collateral: &CollateralKey,
@@ -293,12 +294,12 @@ pub struct Inputs {
 
 impl Inputs {
     /// 同じ時刻に評価する曲線の入力をまとめる。
-    pub(crate) fn new(cases: &[&Case], as_of: Date, extra_discounting: &str) -> Result<Self> {
+    pub fn new(cases: &[&Case], as_of: Date, extra_discounting: &str) -> Result<Self> {
         Self::with_forecasts(cases, as_of, extra_discounting, |_, _| Ok(()))
     }
 
     /// OIS の入力に、別のテストで使う参照データと予想曲線の結び付けを加える。
-    pub(crate) fn with_forecasts(
+    pub fn with_forecasts(
         cases: &[&Case],
         as_of: Date,
         extra_discounting: &str,
@@ -392,7 +393,8 @@ impl Inputs {
     }
 
     /// 校正と準備に使うコンテキスト。
-    pub(crate) const fn context(&self) -> Context<'_> {
+    #[must_use]
+    pub const fn context(&self) -> Context<'_> {
         Context {
             as_of: self.as_of,
             refdata: &self.refdata,
@@ -407,7 +409,7 @@ fn format_error(error: std::fmt::Error) -> Error {
 }
 
 /// 指定した OIS のピラーと校正商品から市場定義を作る。
-pub(crate) fn definition(cases: &[(&Case, &CollateralKey)]) -> Result<MarketDefinition> {
+pub fn definition(cases: &[(&Case, &CollateralKey)]) -> Result<MarketDefinition> {
     let mut curves = Vec::new();
     let mut bindings = Vec::new();
     for &(case, collateral) in cases {
@@ -432,7 +434,7 @@ pub(crate) fn definition(cases: &[(&Case, &CollateralKey)]) -> Result<MarketDefi
 }
 
 /// ピラーに対応する入力 quote を返す。
-pub(crate) fn fixture_quotes(case: &Case) -> Result<Vec<(QuoteId, f64)>> {
+pub fn fixture_quotes(case: &Case) -> Result<Vec<(QuoteId, f64)>> {
     case.nodes
         .iter()
         .map(|node| Ok((case.quote(node)?, node.rate)))
@@ -440,7 +442,7 @@ pub(crate) fn fixture_quotes(case: &Case) -> Result<Vec<(QuoteId, f64)>> {
 }
 
 /// 指定のテナーと固定金利の OIS を準備する。
-pub(crate) fn prepare_ois(
+pub fn prepare_ois(
     case: &Case,
     node: &Node,
     rate: f64,
@@ -449,14 +451,5 @@ pub(crate) fn prepare_ois(
     ctx: &Context<'_>,
 ) -> Result<Prepared> {
     let product = case.template(node)?.build(rate, ctx.refdata)?;
-    prepare(
-        PrepareInput::Standard {
-            product: &product,
-            trade_date: case.as_of,
-        },
-        collateral,
-        &[],
-        market.structure(),
-        ctx,
-    )
+    crate::prepare_standard(&product, case.as_of, collateral, market, ctx)
 }
diff --git a/crates/rhodium-pricing/tests/fixtures/steps.rs b/crates/rhodium-fixtures/src/steps.rs
similarity index 81%
rename from crates/rhodium-pricing/tests/fixtures/steps.rs
rename to crates/rhodium-fixtures/src/steps.rs
index 461517f60a4d27f5a03f8010eebe070185e847a5..de398f2d6ba0d641e0aee5e9952950f1ce3ad56c 100644
--- a/crates/rhodium-pricing/tests/fixtures/steps.rs
+++ b/crates/rhodium-fixtures/src/steps.rs
@@ -20,7 +20,7 @@ use rhodium_pricing::{
 };
 use serde::Deserialize;
 
-use super::ois::{Case, Node, Sample, cases, fixture_quotes, midnight};
+use crate::ois::{Case, Node, Sample, cases, fixture_quotes, midnight};
 
 /// `QuantLib` で年末のターンを含めて校正した TONA の quote と DF。
 #[derive(Deserialize)]
@@ -49,8 +49,8 @@ pub const PENALTY_SCALE: f64 = 1e-4;
 const MEETINGS: [&str; 4] = ["2025-03-19", "2025-06-17", "2025-09-19", "2025-12-19"];
 
 /// 年末のターンの `QuantLib` の fixture を読む。
-pub(crate) fn quantlib_turn() -> Result<QuantLibTurn> {
-    serde_json::from_str(include_str!("quantlib/ois_turn.json"))
+fn quantlib_turn() -> Result<QuantLibTurn> {
+    serde_json::from_str(include_str!("../quantlib/ois_turn.json"))
         .map_err(|error| rhodium_core::Error::InvalidData(error.to_string()))
 }
 
@@ -64,7 +64,7 @@ pub struct PenalizedTona {
 
 impl PenalizedTona {
     /// TONA の慣行と、ターンの fixture の quote を読む。
-    pub(crate) fn new() -> Result<Self> {
+    pub fn new() -> Result<Self> {
         let [mut case, _] = cases()?;
         let turn = quantlib_turn()?;
         case.nodes.clone_from(&turn.nodes);
@@ -72,27 +72,27 @@ impl PenalizedTona {
     }
 
     /// 曲線の ID。
-    pub(crate) fn id(&self) -> Result<CurveId> {
+    pub fn id(&self) -> Result<CurveId> {
         self.case.index.parse()
     }
 
     /// 会合日の階段のブロック。
-    pub(crate) fn meetings(&self) -> Result<BlockId> {
+    pub fn meetings(&self) -> Result<BlockId> {
         BlockId::new(self.id()?, "meetings")
     }
 
     /// 長期の節点のブロック。
-    pub(crate) fn long(&self) -> Result<BlockId> {
+    pub fn long(&self) -> Result<BlockId> {
         BlockId::new(self.id()?, "long")
     }
 
     /// 年末の跳びのブロック。
-    pub(crate) fn turn_block(&self) -> Result<BlockId> {
+    pub fn turn_block(&self) -> Result<BlockId> {
         BlockId::new(self.id()?, "turn")
     }
 
     /// 階段、長期の節点、年末の跳びからなる TONA 曲線。
-    pub(crate) fn curve(&self) -> Result<CurveSpec> {
+    pub fn curve(&self) -> Result<CurveSpec> {
         let nodes = &self.turn.nodes;
         CurveSpec::from_expr(
             self.id()?,
@@ -126,7 +126,7 @@ impl PenalizedTona {
     }
 
     /// 罰則付きの市場の定義。`gauss_newton` なら最小二乗の記録を Gauss–Newton で作る。
-    pub(crate) fn definition(&self, gauss_newton: bool) -> Result<MarketDefinition> {
+    pub fn definition(&self, gauss_newton: bool) -> Result<MarketDefinition> {
         let collateral = CollateralKey::Uncollateralized;
         let mut bindings = self
             .turn
@@ -144,18 +144,7 @@ impl PenalizedTona {
                     .with_tolerance(OIS_TOLERANCE)
             })
             .collect::<Result<Vec<_>>>()?;
-        // 年末のスプレッドは直接マークで、変換 TurnSpreadToJump が s·τ_turn の跳びにする。
-        bindings.push(InstrumentBinding::new(
-            TURN_QUOTE.parse()?,
-            QuoteKind::Spread,
-            CalibrationTarget::DirectMark,
-            QuoteTransform::TurnSpreadToJump(YearFraction::new(
-                self.turn.fraction,
-                self.case.basis(),
-            )?),
-            collateral,
-            ElementSelector::Param(ParamId::new(self.turn_block()?, 0)),
-        )?);
+        bindings.push(self.turn_binding()?);
         let definition = MarketDefinition::new(
             self.case.currency.to_string().parse()?,
             "Etc/UTC".parse()?,
@@ -174,11 +163,26 @@ impl PenalizedTona {
         }
     }
 
-    /// 会合 1 回（1M のピラー）で区切った 2 区間の階段を 1M と 6M の OIS が、1Y で接続した長期の
-    /// 節点を 18M と 2Y の OIS が、年末の跳びを直接マークが決める正方系の市場の定義（R3）。
-    pub(crate) fn square_definition(&self) -> Result<MarketDefinition> {
+    /// 年末のスプレッドの直接マーク。変換 `TurnSpreadToJump` が `s·τ_turn` の跳びにする。
+    fn turn_binding(&self) -> Result<InstrumentBinding> {
+        InstrumentBinding::new(
+            TURN_QUOTE.parse()?,
+            QuoteKind::Spread,
+            CalibrationTarget::DirectMark,
+            QuoteTransform::TurnSpreadToJump(YearFraction::new(
+                self.turn.fraction,
+                self.case.basis(),
+            )?),
+            CollateralKey::Uncollateralized,
+            ElementSelector::Param(ParamId::new(self.turn_block()?, 0)),
+        )
+    }
+
+    /// 会合 1 回（1M のピラー）で区切った 2 区間の階段と、1Y で接続した長期の節点（18M、2Y）と、
+    /// 年末の跳びからなる TONA 曲線（R3）。
+    pub fn square_curve(&self) -> Result<CurveSpec> {
         let nodes = &self.turn.nodes;
-        let curve = CurveSpec::from_expr(
+        CurveSpec::from_expr(
             self.id()?,
             self.case.currency,
             [CurveRole::Discount, CurveRole::Projection],
@@ -203,9 +207,16 @@ impl PenalizedTona {
                     dates: vec![midnight(self.turn.effective)?],
                 },
             ]),
-        )?;
+        )
+    }
+
+    /// [`Self::square_curve`] の校正条件。階段を 1M と 6M の OIS が、長期の節点を 18M と 2Y の
+    /// OIS が、年末の跳びを直接マークが決める。
+    pub fn square_bindings(&self) -> Result<Vec<InstrumentBinding>> {
         let collateral = CollateralKey::Uncollateralized;
-        let mut bindings = nodes
+        let mut bindings = self
+            .turn
+            .nodes
             .iter()
             .enumerate()
             .map(|(position, node)| {
@@ -218,27 +229,22 @@ impl PenalizedTona {
                     .binding(node, &collateral, ElementSelector::Block(block))
             })
             .collect::<Result<Vec<_>>>()?;
-        bindings.push(InstrumentBinding::new(
-            TURN_QUOTE.parse()?,
-            QuoteKind::Spread,
-            CalibrationTarget::DirectMark,
-            QuoteTransform::TurnSpreadToJump(YearFraction::new(
-                self.turn.fraction,
-                self.case.basis(),
-            )?),
-            collateral,
-            ElementSelector::Param(ParamId::new(self.turn_block()?, 0)),
-        )?);
+        bindings.push(self.turn_binding()?);
+        Ok(bindings)
+    }
+
+    /// [`Self::square_curve`] と [`Self::square_bindings`] からなる正方系の市場の定義（R3）。
+    pub fn square_definition(&self) -> Result<MarketDefinition> {
         MarketDefinition::new(
             self.case.currency.to_string().parse()?,
             "Etc/UTC".parse()?,
-            [curve],
-            bindings,
+            [self.square_curve()?],
+            self.square_bindings()?,
         )
     }
 
     /// OIS とターンの quote。
-    pub(crate) fn quotes(&self) -> Result<QuoteSet> {
+    pub fn quotes(&self) -> Result<QuoteSet> {
         QuoteSet::new(
             fixture_quotes(&self.case)?
                 .into_iter()
diff --git a/crates/rhodium-pricing/tests/fixtures/xccy.rs b/crates/rhodium-fixtures/src/xccy.rs
similarity index 89%
rename from crates/rhodium-pricing/tests/fixtures/xccy.rs
rename to crates/rhodium-fixtures/src/xccy.rs
index 7178a552ab2db6f93ec1ecbf6db8db2c3eb5f382..6878a9b34fe32c60da3269a635d69e34197ee4a4 100644
--- a/crates/rhodium-pricing/tests/fixtures/xccy.rs
+++ b/crates/rhodium-fixtures/src/xccy.rs
@@ -38,14 +38,14 @@ use rhodium_products::{
     standard::InstrumentTemplate,
 };
 
-use super::ois::{Case, Inputs, cases, definition, fixture_quotes, midnight};
+use crate::ois::{Case, Inputs, cases, definition, fixture_quotes, midnight};
 
 /// JPY 現金担保の USD 割引曲線（外貨担保曲線）。
-pub const IMPLIED: &str = "USD-under-JPY";
+const IMPLIED: &str = "USD-under-JPY";
 /// `MtM` XCCY の慣行。
 pub const XCCY: &str = "USDJPY-XCCY";
 /// 元本のリセットに使う USDJPY の fixing。
-pub const FX_FIXING: &str = "USDJPY-WMR";
+const FX_FIXING: &str = "USDJPY-WMR";
 /// スポットのアウトライト。
 pub const SPOT: f64 = 150.0;
 /// 短期と日付付きの FX ポイント（pip、0.01 円）と、その受渡日。
@@ -70,7 +70,7 @@ pub struct XccyMarket {
 }
 
 /// JPY 現金担保（TONA 付利、スプレッド 0）。
-pub(crate) fn jpy_cash() -> Result<CollateralKey> {
+pub fn jpy_cash() -> Result<CollateralKey> {
     Ok(CollateralKey::Cash(CashCollateralTerms {
         ccy: Currency::JPY,
         remuneration: Remuneration {
@@ -87,22 +87,18 @@ fn write_error(error: std::fmt::Error) -> Error {
 
 impl XccyMarket {
     /// `QuantLib` の fixture の TONA と SOFR を読む。
-    pub(crate) fn new() -> Result<Self> {
+    pub fn new() -> Result<Self> {
         let [tona, sofr] = cases()?;
         Ok(Self { tona, sofr })
     }
 
     /// 参照データ（FX fixing と XCCY の慣行を含む）と、JPY 現金担保の割引の行を持つ入力。
-    pub(crate) fn inputs(&self, extra_config: &str) -> Result<Inputs> {
+    pub fn inputs(&self, extra_config: &str) -> Result<Inputs> {
         self.inputs_with_spot_lag(extra_config, FxSpotLag::TwoBusinessDays)
     }
 
     /// [`Self::inputs`] の USDJPY の慣行のスポットラグを `spot_lag` にした入力。
-    pub(crate) fn inputs_with_spot_lag(
-        &self,
-        extra_config: &str,
-        spot_lag: FxSpotLag,
-    ) -> Result<Inputs> {
+    pub fn inputs_with_spot_lag(&self, extra_config: &str, spot_lag: FxSpotLag) -> Result<Inputs> {
         let mut config = String::new();
         for (ccy, curve) in [("JPY", "JPY-TONA"), ("USD", IMPLIED)] {
             writeln!(
@@ -173,21 +169,26 @@ impl XccyMarket {
                         legs,
                     })),
                 );
-                if let Some(convention) = parts.fx_conventions.get_mut(&Self::pair()?) {
-                    convention.spot_lag = spot_lag;
-                }
+                // 慣行が無ければスポットラグを書き換えずに通さず、失敗させる。
+                parts
+                    .fx_conventions
+                    .get_mut(&Self::pair()?)
+                    .ok_or_else(|| {
+                        Error::MissingData("fixture has no USDJPY FX convention".into())
+                    })?
+                    .spot_lag = spot_lag;
                 Ok(())
             },
         )
     }
 
     /// 外貨担保曲線の ID。
-    pub(crate) fn implied() -> Result<CurveId> {
+    pub fn implied() -> Result<CurveId> {
         IMPLIED.parse()
     }
 
     /// USDJPY の FX ポイントと長期の XCCY。`long_end` を差し替えると長期の行を変えられる。
-    pub(crate) fn fx(&self, long_end: Vec<FxImpliedLongEnd>) -> Result<FxDefinition> {
+    pub fn fx(&self, long_end: Vec<FxImpliedLongEnd>) -> Result<FxDefinition> {
         let [on, tn, three, six] = POINTS;
         Ok(FxDefinition {
             pair: "USDJPY".parse()?,
@@ -226,7 +227,7 @@ impl XccyMarket {
     }
 
     /// 1Y と 18M の XCCY の長期の行。
-    pub(crate) fn long_end() -> Result<Vec<FxImpliedLongEnd>> {
+    pub fn long_end() -> Result<Vec<FxImpliedLongEnd>> {
         XCCY_QUOTES
             .iter()
             .map(|(quote, tenor, _, pillar)| {
@@ -242,7 +243,7 @@ impl XccyMarket {
     }
 
     /// TONA と SOFR の OIS 曲線と USDJPY の FX からなる市場の定義。
-    pub(crate) fn definition(&self) -> Result<MarketDefinition> {
+    pub fn definition(&self) -> Result<MarketDefinition> {
         let uncollateralized = CollateralKey::Uncollateralized;
         definition(&[
             (&self.tona, &uncollateralized),
@@ -252,7 +253,7 @@ impl XccyMarket {
     }
 
     /// OIS、FX と XCCY の quote。
-    pub(crate) fn quotes(&self) -> Result<QuoteSet> {
+    pub fn quotes(&self) -> Result<QuoteSet> {
         QuoteSet::new(
             fixture_quotes(&self.tona)?
                 .into_iter()
@@ -274,7 +275,7 @@ impl XccyMarket {
     }
 
     /// 長期の XCCY のテンプレート。
-    pub(crate) fn template(tenor: &str) -> Result<InstrumentTemplate> {
+    pub fn template(tenor: &str) -> Result<InstrumentTemplate> {
         Ok(InstrumentTemplate::XccyBasis {
             convention: XCCY.parse()?,
             tenor: tenor.parse::<TenorLength>()?,
@@ -282,12 +283,12 @@ impl XccyMarket {
     }
 
     /// 受渡日の 16:00 UTC（ペアの慣行の受渡時刻）。
-    pub(crate) fn settlement(date: &str) -> Result<Instant> {
+    pub fn settlement(date: &str) -> Result<Instant> {
         format!("{date}T16:00:00Z").parse()
     }
 
     /// 通貨ペア。
-    pub(crate) fn pair() -> Result<CcyPair> {
+    pub fn pair() -> Result<CcyPair> {
         "USDJPY".parse()
     }
 }
diff --git a/crates/rhodium-market/tests/curve_reads.rs b/crates/rhodium-market/tests/curve_reads.rs
index b39aebfba656914c08e335605b9d6c66e6f0e792..37302358c6b177c13e6884fd83bca2f25b42d3b1 100644
--- a/crates/rhodium-market/tests/curve_reads.rs
+++ b/crates/rhodium-market/tests/curve_reads.rs
@@ -15,7 +15,7 @@ use rhodium_market::{
     structure::{BlockId, MarketStructure},
     view::Slots,
 };
-use rhodium_testkit::assert_close;
+use rhodium_testkit::{MARKET_ABS, MARKET_REL, PV_REL, assert_close};
 use scioxide_ad::{Dual, Real, Tape};
 use scioxide_core::{ExecutionContext, InvalidInputReason, SciOxideError, argument};
 use scioxide_interp::{BSpline, ExtrapolationPolicy, Interpolant};
@@ -24,8 +24,6 @@ const ORIGIN: &str = "2026-01-01T00:00:00Z";
 const FIRST: &str = "2027-01-01T00:00:00Z";
 const SECOND: &str = "2028-01-01T00:00:00Z";
 const THIRD: &str = "2028-12-31T00:00:00Z";
-const REL: f64 = 1e-12;
-const ABS: f64 = 1e-14;
 
 fn structure(extrapolation: Extrapolation) -> Result<(Arc<MarketStructure>, CurveId)> {
     let id: CurveId = "JPY-TONA".parse()?;
@@ -117,7 +115,12 @@ fn max_forward_integrates_full_grid_intervals_and_rebinds() -> Result<()> {
         (SECOND, -0.13),
     ] {
         let slots = discount_slots(&structure, &maximum, date)?;
-        assert_close!(df(&market, &slots)?, expected_log.exp(), REL, ABS);
+        assert_close!(
+            df(&market, &slots)?,
+            expected_log.exp(),
+            MARKET_REL,
+            MARKET_ABS
+        );
     }
     assert!(matches!(
         structure.prepare_point(&maximum, "2028-01-01T00:00:00.000001Z".parse()?),
@@ -137,7 +140,12 @@ fn max_forward_integrates_full_grid_intervals_and_rebinds() -> Result<()> {
         Arc::clone(&structure),
         vec![-0.02, -0.02, -0.13, -0.025, -0.07, -0.08],
     )?;
-    assert_close!(df(&rebased, &slots)?, (-0.18_f64).exp(), REL, ABS);
+    assert_close!(
+        df(&rebased, &slots)?,
+        (-0.18_f64).exp(),
+        MARKET_REL,
+        MARKET_ABS
+    );
     let partial = Market::from_partial_theta(
         structure,
         vec![-0.02, -0.02, -0.13, f64::NAN, f64::NAN, f64::NAN],
@@ -197,8 +205,8 @@ fn scale_wraps_nested_max_forward_with_distinct_grids() -> Result<()> {
             assert_close!(
                 df(&market, &slots)?,
                 (-factor * forward_integral).exp(),
-                REL,
-                ABS
+                MARKET_REL,
+                MARKET_ABS
             );
         }
     }
@@ -293,18 +301,18 @@ fn max_forward_shares_ties_equally_in_aad_and_dual()
             &mut gradient,
             ExecutionContext::sequential(),
         )?;
-        assert_close!(value, expected, REL, ABS);
+        assert_close!(value, expected, MARKET_REL, MARKET_ABS);
         for (i, (derivative, share)) in gradient.into_iter().zip(shares).enumerate() {
             let expected_derivative = expected * share;
-            assert_close!(derivative, expected_derivative, REL, ABS);
+            assert_close!(derivative, expected_derivative, MARKET_REL, MARKET_ABS);
             let mut dual: Vec<_> = theta.iter().map(|&value| Dual::new(value, 0.0)).collect();
             dual[i] = Dual::variable(theta[i]);
             let market = Market::from_theta(Arc::clone(&structure), dual)?;
             assert_close!(
                 df(&market, &slots)?.derivative(),
                 expected_derivative,
-                REL,
-                ABS
+                MARKET_REL,
+                MARKET_ABS
             );
         }
     }
@@ -338,11 +346,16 @@ fn flat_forwards_matches_hand_calculated_discount_factors_and_forwards() -> Resu
         let bound = view.bind_slots(&slots)?;
         let expected = [1.0, 1.04_f64.sqrt(), 1.04, (1.04_f64 * 0.96).sqrt(), 0.96];
         for (&at, expected) in bound.discounts()[0].points().iter().zip(expected) {
-            assert_close!(view.df(at).value(), expected, REL, ABS);
+            assert_close!(view.df(at).value(), expected, MARKET_REL, MARKET_ABS);
         }
         // モデル時間が 1 年でも、予想の分母には指標の ACT/360 を使う。
         let simple = view.forward(bound.projections()[0].points()[0]);
-        assert_close!(simple, (1.04 / 0.96 - 1.0) * 360.0 / 365.0, REL, ABS);
+        assert_close!(
+            simple,
+            (1.04 / 0.96 - 1.0) * 360.0 / 365.0,
+            MARKET_REL,
+            MARKET_ABS
+        );
         Ok(())
     })
 }
@@ -411,7 +424,7 @@ fn flat_forwards_agree_with_scioxide_linear_spline_of_log_discount_factors()
             [model_time(as_of, at).years()],
             ExtrapolationPolicy::MethodDefined,
         )?;
-        assert_close!(df, log_df.exp(), REL, ABS);
+        assert_close!(df, log_df.exp(), MARKET_REL, MARKET_ABS);
     }
     Ok(())
 }
@@ -427,7 +440,7 @@ fn extrapolation_keeps_the_selected_forward_slope() -> Result<()> {
         let (structure, id) = structure(mode)?;
         let slots = discount_slots(&structure, &id, THIRD)?;
         let market = Market::from_theta(structure, vec![-0.02, -0.06])?;
-        assert_close!(df(&market, &slots)?, log_df.exp(), REL, ABS);
+        assert_close!(df(&market, &slots)?, log_df.exp(), MARKET_REL, MARKET_ABS);
     }
     let (structure, id) = structure(Extrapolation::Reject)?;
     assert!(matches!(
@@ -440,7 +453,12 @@ fn extrapolation_keeps_the_selected_forward_slope() -> Result<()> {
     ));
     let last = discount_slots(&structure, &id, SECOND)?;
     let market = Market::from_theta(structure, vec![-0.02, -0.06])?;
-    assert_close!(df(&market, &last)?, (-0.06_f64).exp(), REL, ABS);
+    assert_close!(
+        df(&market, &last)?,
+        (-0.06_f64).exp(),
+        MARKET_REL,
+        MARKET_ABS
+    );
     Ok(())
 }
 
@@ -450,7 +468,12 @@ fn instantaneous_extrapolation_uses_the_origin_for_one_pillar() -> Result<()> {
     let (structure, id, _) = two_curves()?;
     let slots = discount_slots(&structure, &id, SECOND)?;
     let market = Market::from_theta(structure, vec![0.98_f64.ln(), -0.03])?;
-    assert_close!(df(&market, &slots)?, 0.98_f64.powi(2), REL, ABS);
+    assert_close!(
+        df(&market, &slots)?,
+        0.98_f64.powi(2),
+        MARKET_REL,
+        MARKET_ABS
+    );
     Ok(())
 }
 
@@ -469,7 +492,7 @@ fn flat_last_segment_stays_flat_when_extrapolated_from_close_pillars() -> Result
     let slots = discount_slots(&structure, &id, THIRD)?;
     let log_df = 0.98_f64.ln();
     let market = Market::from_theta(structure, vec![log_df, log_df])?;
-    assert_close!(df(&market, &slots)?, 0.98, REL, ABS);
+    assert_close!(df(&market, &slots)?, 0.98, MARKET_REL, MARKET_ABS);
     Ok(())
 }
 
@@ -487,8 +510,8 @@ fn forward_retains_small_log_discount_differences() -> Result<()> {
         assert_close!(
             view.forward(bound.projections()[0].points()[0]),
             1e-11,
-            REL,
-            ABS
+            MARKET_REL,
+            MARKET_ABS
         );
         Ok(())
     })
@@ -508,8 +531,8 @@ fn growth_log_retains_small_log_discount_differences() -> Result<()> {
         let bound = view.bind_slots(&slots)?;
         let at = bound.projection_growth()[0].points()[0];
         assert_eq!(view.growth(at), 1.0);
-        assert_close!(view.log_growth(at), 1e-17, REL, 0.0);
-        assert_close!(view.log_growth(at).exp_m1(), 1e-17, REL, 0.0);
+        assert_close!(view.log_growth(at), 1e-17, MARKET_REL, 0.0);
+        assert_close!(view.log_growth(at).exp_m1(), 1e-17, MARKET_REL, 0.0);
         Ok(())
     })
 }
@@ -521,8 +544,8 @@ fn slots_rebind_on_the_same_structure_with_new_theta() -> Result<()> {
     let slots = discount_slots(&structure, &id, FIRST)?;
     let first = Market::from_theta(Arc::clone(&structure), vec![0.98_f64.ln(), -0.06])?;
     let second = Market::from_theta(structure, vec![0.97_f64.ln(), -0.06])?;
-    assert_close!(df(&first, &slots)?, 0.98, REL, ABS);
-    assert_close!(df(&second, &slots)?, 0.97, REL, ABS);
+    assert_close!(df(&first, &slots)?, 0.98, MARKET_REL, MARKET_ABS);
+    assert_close!(df(&second, &slots)?, 0.97, MARKET_REL, MARKET_ABS);
     Ok(())
 }
 
@@ -548,8 +571,13 @@ fn curve_reads_preserve_dual_sensitivities() -> Result<()> {
         let bound = view.bind_slots(&slots)?;
         let df = view.df(bound.discounts()[0].points()[0]).value();
         let forward = view.forward(bound.projections()[0].points()[0]);
-        assert_close!(df.derivative(), (0.98_f64 * 0.94).sqrt() / 2.0, REL, ABS);
-        assert_close!(forward.derivative(), 0.98 / 0.94, REL, ABS);
+        assert_close!(
+            df.derivative(),
+            (0.98_f64 * 0.94).sqrt() / 2.0,
+            MARKET_REL,
+            MARKET_ABS
+        );
+        assert_close!(forward.derivative(), 0.98 / 0.94, MARKET_REL, MARKET_ABS);
         Ok(())
     })
 }
@@ -585,9 +613,9 @@ fn currency_scope_preserves_the_aad_gradient() -> std::result::Result<(), Box<dy
         &mut gradient,
         ExecutionContext::sequential(),
     )?;
-    assert_close!(pv, 98.0, REL, 1e-6);
-    assert_close!(gradient[0], 98.0, REL, 1e-6);
-    assert_close!(gradient[1], 0.0, REL, 1e-6);
+    assert_close!(pv, 98.0, PV_REL, 1e-6);
+    assert_close!(gradient[0], 98.0, PV_REL, 1e-6);
+    assert_close!(gradient[1], 0.0, PV_REL, 1e-6);
     Ok(())
 }
 
@@ -710,8 +738,8 @@ fn curve_resolution_checks_role_and_market_availability() -> Result<()> {
         assert_close!(
             view.df(bound.discounts()[0].points()[0]).value(),
             0.97,
-            REL,
-            ABS
+            MARKET_REL,
+            MARKET_ABS
         );
         Ok(())
     })
@@ -737,7 +765,7 @@ fn currency_scope_discounts_and_tags_present_values() -> Result<()> {
                 scope.pv(100.0, pay) - scope.pv(2.0, pay).scale(3.0) + -scope.present(-5.0);
             let money = scope.into_money(amount);
             assert_eq!(money.ccy(), Currency::USD);
-            assert_close!(money.amount(), 95.24, REL, 1e-8);
+            assert_close!(money.amount(), 95.24, PV_REL, 1e-8);
             Ok(())
         })
     })
@@ -839,7 +867,7 @@ fn bound_pay_rejects_another_discount_curve_of_same_currency() -> Result<()> {
             let pv = scope
                 .into_money(scope.pv(100.0, scope.bind_pay(at)?))
                 .amount();
-            assert_close!(pv, 100.0 * (-0.02_f64).exp(), REL, 1e-6);
+            assert_close!(pv, 100.0 * (-0.02_f64).exp(), PV_REL, 1e-6);
             Ok(())
         })
     })
@@ -936,15 +964,15 @@ fn composed_curves_match_manual_log_discount_factors() -> Result<()> {
                 &discount_slots(&structure, &tona, &at.to_string())?
             )?,
             log_expected.exp(),
-            REL,
-            ABS
+            MARKET_REL,
+            MARKET_ABS
         );
     }
     assert_close!(
         df(&market, &discount_slots(&structure, &sofr, SECOND)?)?,
         (-0.08_f64).exp(),
-        REL,
-        ABS
+        MARKET_REL,
+        MARKET_ABS
     );
 
     let start = structure.prepare_point(&tibor, boundary)?;
@@ -962,8 +990,8 @@ fn composed_curves_match_manual_log_discount_factors() -> Result<()> {
         assert_close!(
             view.forward(bound.projections()[0].points()[0]),
             0.032_f64.exp_m1(),
-            REL,
-            ABS
+            MARKET_REL,
+            MARKET_ABS
         );
         Ok(())
     })?;
@@ -980,8 +1008,8 @@ fn composed_curves_match_manual_log_discount_factors() -> Result<()> {
             &discount_slots(&structure, &tona, &meeting.to_string())?
         )?,
         (-0.005_f64).exp(),
-        REL,
-        ABS
+        MARKET_REL,
+        MARKET_ABS
     );
     assert!(matches!(
         partial.with_view(|view| view.bind_slots(&slots).map(|_| ())),
diff --git a/crates/rhodium-pricing/Cargo.toml b/crates/rhodium-pricing/Cargo.toml
index ddfa6866ea830ba58bfe5d5b38310a963920ae8d..caa65f62f1fcef77af6ed9d8fb91e9c3d0fcb5a4 100644
--- a/crates/rhodium-pricing/Cargo.toml
+++ b/crates/rhodium-pricing/Cargo.toml
@@ -24,6 +24,7 @@ scioxide-special.workspace = true
 serde.workspace = true
 
 [dev-dependencies]
+rhodium-fixtures.workspace = true
 rhodium-testkit.workspace = true
 serde_json = "1"
 toml = "1.2"
diff --git a/crates/rhodium-pricing/src/calibration.rs b/crates/rhodium-pricing/src/calibration.rs
index ee0fdaed86ae4c91f3d962465c6cf9518c482132..2d45901dd28ce7726ee0bf476acb45275fbbee04 100644
--- a/crates/rhodium-pricing/src/calibration.rs
+++ b/crates/rhodium-pricing/src/calibration.rs
@@ -8,14 +8,6 @@ mod groups;
 mod solver;
 mod vol;
 
-#[cfg(test)]
-#[allow(
-    dead_code,
-    reason = "統合テスト用の QuantLib fixture のうち、再試行には TONA 曲線だけを使う"
-)]
-#[path = "../tests/fixtures/ois.rs"]
-mod ois_fixture;
-
 use std::{
     collections::{BTreeMap, BTreeSet},
     sync::Arc,
@@ -1102,29 +1094,206 @@ impl PreparedBinding {
 
 #[cfg(test)]
 mod tests {
-    use rhodium_market::quote::{QuoteOverride, QuoteOverrides};
-    use rhodium_testkit::assert_close;
-
-    use super::{
-        ois_fixture::{Inputs, cases, definition, fixture_quotes},
-        *,
+    use rhodium_core::{
+        calendar::{BaseCalendar, Extrapolation as CalendarExtrapolation, WeekendRule},
+        date::Weekday,
+        day_count::DayCount,
+        fixing::{FixingStore, FixingStoreParts},
+        id::{CalendarId, IndexId},
+        observation::{
+            Aggregation, CompoundingForm, ObservationDays, ObservationMethod, ObservationRounding,
+            ObservationRule, SpreadTreatment,
+        },
+        refdata::{
+            Convention, CurrencyDefinition, DailyPublication, Deliverability, FixedLegConvention,
+            IndexDefinition, IndexKind, IndexStatus, OisConvention, OvernightLegConvention,
+            PublicationRule, RefData, RefDataParts,
+        },
+        rounding::{Rounding, RoundingMode},
+        schedule::{PaymentGrouping, PaymentRule, PeriodAnchor},
+        tenor::{DateOffset, Frequency, TenorLength, TenorUnit},
+        time::{LocalEventTime, UtcOffset, ZoneRules},
     };
+    use rhodium_market::{
+        curve::{CurveSpec, Extrapolation},
+        quote::{QuoteOverride, QuoteOverrides},
+    };
+    use rhodium_products::standard::InstrumentTemplate;
+    use rhodium_testkit::{MARKET_ABS, MARKET_REL, assert_close};
+
+    use super::*;
+    use crate::{
+        config::PricingConfig,
+        definition::{
+            CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
+        },
+        resolve::CollateralKey,
+    };
+
+    /// 評価日。曲線のピラーと校正商品の約定日の起点。
+    const AS_OF: &str = "2025-01-07";
+    /// 校正商品の OIS のテナーとパーレート。慣行どおりの年払いで、ピラーは最終の支払日。
+    const NODES: [(&str, f64, &str); 2] =
+        [("1Y", 0.0025, "2026-01-09"), ("2Y", 0.0030, "2027-01-11")];
+
+    /// 週末だけの暦、UTC 0 時の公表と決済を持つ、JPY の TONA の OIS だけの参照データ。
+    ///
+    /// 単体テストは非公開のフィールドを書き換えるので、pricing に依存する共有の fixture
+    /// （`rhodium-fixtures`）を使えない。dev-dependency の循環で `rhodium_pricing` が二つになり、
+    /// 型が合わないためである。ここでは再試行の経路を通すのに要る最小の入力だけを組む。
+    fn tona_refdata() -> Result<RefData> {
+        let clock = LocalEventTime {
+            time: "00:00:00".parse()?,
+            zone: "Etc/UTC".parse()?,
+        };
+        let calendar_id: CalendarId = "JPY".parse()?;
+        let calendar = CalendarSpec::Base(calendar_id.clone());
+        let index: IndexId = "JPY-TONA".parse()?;
+        let mut parts = RefDataParts::default();
+        parts.zones.insert(
+            clock.zone.clone(),
+            ZoneRules::new(
+                "2024-01-01T00:00:00Z".parse()?..="2028-01-01T00:00:00Z".parse()?,
+                UtcOffset::from_seconds(0)?,
+                Vec::new(),
+            )?,
+        );
+        parts.calendars.insert(
+            calendar_id.clone(),
+            BaseCalendar::new(
+                "2024-01-01".parse()?..="2027-12-31".parse()?,
+                WeekendRule::new([Weekday::Sat, Weekday::Sun]),
+                Vec::new(),
+                CalendarExtrapolation::Strict,
+            )?,
+        );
+        parts.currencies.insert(
+            Currency::JPY,
+            CurrencyDefinition {
+                code: Currency::JPY,
+                minor_units: 0,
+                payment_rounding: Rounding::new(1.0, RoundingMode::HalfEven)?,
+                deliverability: Deliverability::Deliverable,
+                settlement_time: clock.clone(),
+                settlement_calendar: calendar_id,
+            },
+        );
+        parts.indices.insert(
+            index.clone(),
+            IndexDefinition {
+                id: index.clone(),
+                kind: IndexKind::OvernightRate {
+                    ccy: Currency::JPY,
+                    day_count: DayCount::Act365F,
+                },
+                calendar: calendar.clone(),
+                publication: PublicationRule::Daily(DailyPublication {
+                    business_days: 1,
+                    calendar: calendar.clone(),
+                    time: clock.time,
+                    zone: clock.zone.clone(),
+                }),
+                observation_time: Some(clock),
+                status: IndexStatus::Active,
+            },
+        );
+        let observation = ObservationRule::new(
+            ObservationDays::BusinessDays(calendar.clone()),
+            ObservationMethod::Plain,
+            Aggregation::Compound {
+                form: CompoundingForm::SimpleFactor,
+                spread: SpreadTreatment::Exclusive,
+            },
+            ObservationRounding::default(),
+        )?;
+        parts.conventions.insert(
+            "JPY-TONA".parse()?,
+            Convention::Ois(Box::new(OisConvention {
+                spot_lag: DateOffset::new(2, TenorUnit::BusinessDay),
+                fixed: FixedLegConvention {
+                    frequency: Frequency::ANNUAL,
+                    day_count: DayCount::Act365F,
+                    bdc: BusinessDayConvention::ModifiedFollowing,
+                    calendar: calendar.clone(),
+                },
+                float: OvernightLegConvention { index, observation },
+                payment: PaymentRule {
+                    anchor: PeriodAnchor::PeriodEnd,
+                    offset: DateOffset::new(0, TenorUnit::BusinessDay),
+                    calendar,
+                    bdc: BusinessDayConvention::Following,
+                    grouping: PaymentGrouping::default(),
+                    settlement: None,
+                },
+            })),
+        );
+        RefData::new(parts)
+    }
+
+    /// TONA の OIS の 2 本のピラーを、それぞれの OIS のパーレートで当てる市場の定義。
+    fn tona_definition() -> Result<MarketDefinition> {
+        let curve = CurveSpec::new(
+            "JPY-TONA".parse()?,
+            Currency::JPY,
+            [CurveRole::Discount, CurveRole::Projection],
+            NODES
+                .iter()
+                .map(|(_, _, pillar)| format!("{pillar}T00:00:00Z").parse())
+                .collect::<Result<Vec<_>>>()?,
+            Extrapolation::Reject,
+        )?;
+        let block = BlockId::nodes("JPY-TONA".parse()?);
+        let bindings = NODES
+            .iter()
+            .map(|(tenor, _, _)| {
+                InstrumentBinding::new(
+                    format!("JPY.TONA.OIS.{tenor}").parse()?,
+                    QuoteKind::Rate,
+                    CalibrationTarget::ParRate(InstrumentTemplate::Ois {
+                        convention: "JPY-TONA".parse()?,
+                        tenor: tenor.parse::<TenorLength>()?,
+                    }),
+                    QuoteTransform::Identity,
+                    CollateralKey::Uncollateralized,
+                    ElementSelector::Block(block.clone()),
+                )
+            })
+            .collect::<Result<Vec<_>>>()?;
+        MarketDefinition::new("JPY".parse()?, "Etc/UTC".parse()?, [curve], bindings)
+    }
 
     #[test]
     fn quote_update_retries_with_a_flat_seed_when_the_previous_theta_fails() -> Result<()> {
         // 前回の θ を、残差の最初の評価で DF を表現できない値に置き換える。
-        // 指定の OIS は一意で、quote の中央値を用いた初期値からは普通に解ける。
-        let [tona, _] = cases()?;
-        let collateral = crate::resolve::CollateralKey::Uncollateralized;
-        let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
-        let def = definition(&[(&tona, &collateral)])?;
-        let quotes = QuoteSet::new(fixture_quotes(&tona)?)?;
-        let ctx = inputs.context();
+        // OIS のピラーは一つずつ一意に決まるので更新は前の θ から warm start し、それに失敗すると
+        // quote の中央値を用いた初期値で解き直す。解き直した θ は全体の再校正と一致する。
+        let refdata = tona_refdata()?;
+        let fixings = FixingStore::new(&refdata, FixingStoreParts::default())?;
+        let config: PricingConfig = toml::from_str(
+            "version = 1\nreporting_currency = \"JPY\"\n\
+             [[discounting]]\ncurrency = \"JPY\"\ncollateral = \"uncollateralized\"\n\
+             curve = \"JPY-TONA\"\n[forecasting]\n\"JPY-TONA\" = \"JPY-TONA\"\n",
+        )
+        .map_err(|error| Error::InvalidData(error.to_string()))?;
+        let ctx = Context {
+            as_of: format!("{AS_OF}T00:00:00Z").parse()?,
+            refdata: &refdata,
+            fixings: &fixings,
+            config: &config,
+        };
+        let def = tona_definition()?;
+        let quotes = QuoteSet::new(
+            NODES
+                .iter()
+                .map(|(tenor, rate, _)| Ok((format!("JPY.TONA.OIS.{tenor}").parse()?, *rate)))
+                .collect::<Result<Vec<_>>>()?,
+        )?;
         let mut base = calibrate(&def, &quotes, &ctx)?;
         base.market = Market::from_theta(
             Arc::clone(base.market().structure()),
             vec![-2000.0; base.market().theta().len()],
         )?;
+        // 置き換えた θ では、校正商品の残差の評価そのものが失敗する。
         assert!(
             base.prepared
                 .rows
@@ -1137,8 +1306,8 @@ mod tests {
         );
         let overrides = QuoteOverrides::new(
             [QuoteOverride::Set {
-                id: tona.quote(&tona.nodes[0])?,
-                value: tona.nodes[0].rate + 0.0001,
+                id: format!("JPY.TONA.OIS.{}", NODES[0].0).parse()?,
+                value: NODES[0].1 + 0.0001,
             }],
             def.quote_kinds(),
         )?;
@@ -1149,8 +1318,9 @@ mod tests {
             updated.groups()[0].outcome(),
             GroupOutcome::Solved(_)
         ));
+        assert_eq!(updated.market().theta().len(), full.market().theta().len());
         for (&actual, &expected) in updated.market().theta().iter().zip(full.market().theta()) {
-            assert_close!(actual, expected, 1e-12, 1e-14);
+            assert_close!(actual, expected, MARKET_REL, MARKET_ABS);
         }
         Ok(())
     }
diff --git a/crates/rhodium-pricing/src/lib.rs b/crates/rhodium-pricing/src/lib.rs
index 55db1bfba1024f0b51ffa53bfaf70d843bc9b461..b6f23d3421833206fab0747df064359cff09bbe5 100644
--- a/crates/rhodium-pricing/src/lib.rs
+++ b/crates/rhodium-pricing/src/lib.rs
@@ -2,10 +2,6 @@
 
 #![forbid(unsafe_code)]
 
-// 校正の単体テストでも統合テストと同じ市場 fixture を使う。
-#[cfg(test)]
-extern crate self as rhodium_pricing;
-
 pub mod analytic;
 pub mod calibration;
 pub mod config;
diff --git a/crates/rhodium-pricing/tests/calibration_instruments.rs b/crates/rhodium-pricing/tests/calibration_instruments.rs
index 44dbce47249b110cc17f3b96696de6a5e1a3b10b..83d8fcaa1d5903bad66e53f224a44ceea119715f 100644
--- a/crates/rhodium-pricing/tests/calibration_instruments.rs
+++ b/crates/rhodium-pricing/tests/calibration_instruments.rs
@@ -5,50 +5,26 @@
 //! quote で約定した取引として評価すると、PV は校正の残差の許容誤差の内でゼロになる。
 
 use rhodium_core::{Error, Result, currency::Currency};
+use rhodium_fixtures::{COLLATERAL, solved_within};
+use rhodium_fixtures::{multicurve, ois};
 use rhodium_market::{
     market::Market,
     quote::{QuoteId, QuoteKind, QuoteSet},
 };
 use rhodium_pricing::{
-    calibration::{CalibratedMarket, GroupOutcome, calibrate},
+    calibration::{CalibratedMarket, calibrate},
     context::Context,
     definition::{CalibrationTarget, ElementSelector, InstrumentBinding, QuoteTransform},
     prepare::{PrepareInput, evaluate, prepare},
-    resolve::CollateralKey,
 };
 use rhodium_products::standard::InstrumentTemplate;
 use rhodium_testkit::assert_close;
 
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、校正商品の検査に使わない関数がある"
-)]
-#[path = "fixtures/ois.rs"]
-pub mod ois;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、校正商品の検査に使わない関数がある"
-)]
-#[path = "fixtures/multicurve.rs"]
-pub mod multicurve;
-
 use multicurve::{MulticurveCase, all_multicurve_cases};
 use ois::{Node, fixture_quotes};
 
-const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
 const TOLERANCE: f64 = 1e-14;
 
-fn solved(market: &CalibratedMarket) -> Result<()> {
-    for group in market.groups() {
-        let GroupOutcome::Solved(solution) = group.outcome() else {
-            return Err(Error::Calibration(format!("{group:?}")));
-        };
-        assert!(solution.residual_norm() <= TOLERANCE);
-    }
-    Ok(())
-}
-
 /// テンプレートの PV が quote K について affine なので、K = 0 と 1 の PV からパーと傾きを求める。
 fn par_and_slope(
     template: &InstrumentTemplate,
@@ -97,7 +73,7 @@ fn irs_quotes_calibrate_the_tibor_curve_and_reprice_at_par() -> Result<()> {
         let inputs = case.inputs()?;
         let ctx = inputs.context();
         let basis = calibrate(&case.definition()?, &case.quotes()?, &ctx)?;
-        solved(&basis)?;
+        solved_within(&basis, TOLERANCE)?;
         let term_block = ElementSelector::Block(case.term_block()?);
         let mut bindings = Vec::new();
         let mut quotes = fixture_quotes(&case.overnight)?;
@@ -121,7 +97,7 @@ fn irs_quotes_calibrate_the_tibor_curve_and_reprice_at_par() -> Result<()> {
             &quotes,
             &ctx,
         )?;
-        solved(&swaps)?;
+        solved_within(&swaps, TOLERANCE)?;
         assert_same_theta(&basis, &swaps);
         // 校正した IRS を quote の固定金利で約定した取引の PV は、残差の許容誤差 × 傾きの幅でゼロ。
         for node in &case.nodes {
@@ -199,7 +175,7 @@ fn deposit_template_calibrates_and_reprices_at_par() -> Result<()> {
             .chain([(deposit, par)]),
     )?;
     let deposits = calibrate(&definition, &quotes, &ctx)?;
-    solved(&deposits)?;
+    solved_within(&deposits, TOLERANCE)?;
     assert_same_theta(&original, &deposits);
     let (repriced, _) = par_and_slope(&template, &case, deposits.market(), &ctx)?;
     assert_close!(repriced, par, 0.0, TOLERANCE);
diff --git a/crates/rhodium-pricing/tests/fixtures/quantlib/README.md b/crates/rhodium-pricing/tests/fixtures/quantlib/README.md
index 5cf163a8ddf247d8d299d1319f527c8a100b4d9c..822fdff652c1a3c3940d39b5865bbcc45c775985 100644
--- a/crates/rhodium-pricing/tests/fixtures/quantlib/README.md
+++ b/crates/rhodium-pricing/tests/fixtures/quantlib/README.md
@@ -1,4 +1,6 @@
-# QuantLib のスワップ PV
+# QuantLib のスワップ PV とボラ
+
+TONA と SOFR の OIS 曲線、会合日のステップと年末のターン、ターム指標の複数曲線の fixture は、複数の crate のテストが共有するので `crates/rhodium-fixtures/quantlib/` にある（同じ場所の `README.md` に出典と生成の手順を書く）。
 
 QuantLib 1.42 の Python バインディングで、固定払い・翌日物複利受けのスワップを評価した。
 [QuantLib](https://github.com/lballabio/QuantLib/tree/v1.42) と [QuantLib-SWIG](https://github.com/lballabio/QuantLib-SWIG/tree/v1.42) は Modified BSD License で公開されている。
@@ -21,70 +23,6 @@ Rust のテストは保存した PV と脚別 PV を照合する。
 相対誤差は 1e-12、絶対誤差の下限は通貨の最小単位の 1e-6 倍とする。
 ネット PV の比較では両脚の絶対値の和も相対誤差の尺度に含める。
 
-## TONA と SOFR の OIS 曲線
-
-`ois_curves.json` は、QuantLib 1.42 の `OISRateHelper` と
-`PiecewiseLogLinearDiscount` で 1M・6M・1Y・18M・2Y の quote をブートストラップした結果である。
-同じ曲線で割引と翌日物の予想を行う。
-約定日は 2025 年 1 月 7 日、スポットラグは 2 営業日、固定脚は年払いとする。
-両実装で満期から逆算して短い初回期間を作り、期間境界に Modified Following、支払いに Following を適用する。
-JPY は Japan カレンダーと ACT/365F で支払ラグは 0、USD は SOFR の fixing カレンダーと ACT/360 で支払ラグは 2 営業日とする。
-ピラーは QuantLib の `LastRelevantDate`（各商品の最終支払日）とする。
-評価、観測と支払いを UTC 00:00 にそろえ、曲線の時間軸は両実装で ACT/365F とする。
-
-```sh
-uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py > crates/rhodium-pricing/tests/fixtures/quantlib/ois_curves.json
-```
-
-Rust では各ピラーの DF と、評価日から最後のピラーまでの 5 つの区間に一つずつ置いた照合点（2025-01-24、2025-04-15、2025-08-12、2026-04-15、2026-10-15）の DF を、相対 1e-12、絶対 1e-14 の下限で照合する。
-両実装の停止判定と浮動小数点の演算順が異なるため、この許容誤差を使う。
-quote を固定金利にした校正商品を取引と同じ準備で評価し、PV が年金現価 × 2 × 残差の許容誤差（1e-14）以内で 0 になることも確かめる。
-
-同じ曲線で、QuantLib の `OvernightIndexedSwap`（日程は年払い、満期から逆算、期間境界に Modified Following、支払いに Following と各通貨の支払ラグ）を二つ評価し、PV と脚別の PV を `swaps` に置く。
-
-| 名前 | 開始 | テナー | 固定金利 | スプレッド | 固定脚 | 元本 |
-|---|---|---|---|---|---|---|
-| `spot` | 約定日（評価日）の 2 営業日後 | 18M | 18M の quote + 25bp | 0 | 受け | JPY 1 億、USD 100 万 |
-| `seasoned` | 2024-12-05 | 1Y | JPY 0.35%、USD 4.20% | 10bp | 払い | JPY 1 億、USD 100 万 |
-
-開始済みのスワップは、2024-12-05 から評価日の前日までの各営業日の fixing（`fixings`、JPY は 0.23% から、USD は 4.55% から 1bp ずつ 3 日周期で変える値）を確定値として読み、評価日の観測は曲線から予想する。
-rhodium の参照データでは fixing の公表を翌営業日の UTC 00:00 とするので、評価日の前日の値は評価の瞬間ちょうどに公表済みになり、評価日の観測期間は評価の瞬間に始まる。
-QuantLib はスプレッドを日次の複利の外で足す（rhodium の `Exclusive`）。生成スクリプトは、スプレッドの有無の変動脚の差が Σ 元本 × スプレッド × 計算期間 × DF に一致することを確かめてから出力する。
-Rust は QuantLib のノードの DF を θ に入れた市場で、同じ取引を `Swap::from_convention`（スポット起点は開始日を省いて慣行のスポットラグで求める）と同じ fixing で作り、PV、変動脚（固定金利 0 の PV）、固定脚（差）を照合する。
-許容誤差は相対 1e-12、絶対の下限は通貨の最小単位の 1e-6 倍とし、相対の尺度には両脚の絶対値の和と、変動脚の二つの項（N·P(s)/P(e) と N）の大きさとして元本の 2 倍を含める。
-
-## 会合日のステップと年末のターンを持つ TONA
-
-`ois_turn.json` は同じ TONA の quote を、年末のターンを含む `PiecewiseLogLinearDiscount` に校正した結果である。
-Japan カレンダーで期末の営業日から翌営業日までの日数係数を求め、50bp を掛けた値を ln DF の跳びにする。
-QuantLib のジャンプは指定日より後に効くため、QuantLib には発効の前日を指定する。
-rhodium は最初のピラーで会合のステップを切り替え、二つ目のピラーで長期の節点につなぐ。
-両実装のピラーと中間日の DF を相対 1e-12、絶対 1e-14 で比べ、JPY と USD の OIS quote を再評価する。
-
-```sh
-uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py --turn > crates/rhodium-pricing/tests/fixtures/quantlib/ois_turn.json
-```
-
-## TONA と TIBOR 6M・3M、SOFR と USD ターム指標の複数曲線
-
-`multicurve.json` は、QuantLib 1.42 の `OISRateHelper` で OIS 曲線を先に校正し、`OvernightIborBasisSwapRateHelper` で別のターム予想曲線を校正した結果である。
-OIS とベーシスには 6M、1Y、18M、2Y を使う。
-ベーシス商品は 2 営業日後に開始し、半年ごとにラグなしで支払う。
-翌日物脚に複利の外で spread を足して払い、ターム脚を受ける。
-ターム指標は JPY では QuantLib の TIBOR 6M と TIBOR 3M、USD では SOFR のカレンダーを使うテスト用の 6M 指標で、予想期間の端点は指標の fixing 規則で決める。
-両脚の支払いの頻度はターム指標のテナー（TIBOR 3M は四半期ごと）で、出力の `term_tenor` に記録する。
-配列の順は TIBOR 6M、USD のターム指標、TIBOR 3M（三つ目は M1 で足した）。
-OIS 曲線を両脚の割引と翌日物の予想に使う。
-
-JPY は割引とターム予想のピラーがそろっており、rhodium ではターム予想を TONA とベーシスの和で表す。
-共通のピラーの間ではどちらの実装も対数 DF が線形になる。
-USD の OIS は 2 営業日の支払ラグがあるため、ターム予想は独立した節点曲線で表す。
-Rust のテストはピラーと補間点の割引係数および予想の擬似 DF を相対 1e-12、絶対 1e-14 で比較し、元本 1 の校正商品を quote どおりに評価した PV を残差の許容誤差から求めた幅で照合する。
-
-```sh
-uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py --multicurve > crates/rhodium-pricing/tests/fixtures/quantlib/multicurve.json
-```
-
 ## 金利ボラの SABR とキャップレット
 
 `volatility.json` は QuantLib 1.42 の `shiftedSabrVolatility` で作った正規・対数正規の Hagan ボラと、
diff --git a/crates/rhodium-pricing/tests/forecast_measure.rs b/crates/rhodium-pricing/tests/forecast_measure.rs
index bba9eac681ec1aba2f9f9fa8ad7c9e8e71915c07..bfb80d5c86b4b787e0cbe8f0a11b2ab0fc2d8d98 100644
--- a/crates/rhodium-pricing/tests/forecast_measure.rs
+++ b/crates/rhodium-pricing/tests/forecast_measure.rs
@@ -7,12 +7,14 @@
 use std::sync::Arc;
 
 use rhodium_core::{Error, Result, Warning, currency::Currency, money::Money, quantity::Spread};
+use rhodium_fixtures::solved;
+use rhodium_fixtures::{multicurve, ois};
 use rhodium_market::{
     curve::{CurveExpr, CurveId, CurveRole, CurveSpec},
     structure::MarketStructure,
 };
 use rhodium_pricing::{
-    calibration::{CalibratedMarket, GroupOutcome, calibrate},
+    calibration::{CalibratedMarket, calibrate},
     context::Context,
     prepare::{PrepareInput, Prepared, prepare},
     resolve::CollateralKey,
@@ -22,20 +24,6 @@ use rhodium_products::{
     standard::{Fra, Position, StandardProduct},
 };
 
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、担保の判定に使わない関数がある"
-)]
-#[path = "fixtures/ois.rs"]
-pub mod ois;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、担保の判定に使わない関数がある"
-)]
-#[path = "fixtures/multicurve.rs"]
-pub mod multicurve;
-
 use multicurve::{MulticurveCase, multicurve_cases};
 use ois::{Inputs, prepare_ois};
 
@@ -72,11 +60,7 @@ fn calibrated(
     ctx: &Context<'_>,
 ) -> Result<CalibratedMarket> {
     let result = calibrate(&case.definition_with(extra)?, &case.quotes()?, ctx)?;
-    for group in result.groups() {
-        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
-            return Err(Error::Calibration(format!("{group:?}")));
-        }
-    }
+    solved(&result)?;
     Ok(result)
 }
 
diff --git a/crates/rhodium-pricing/tests/fra.rs b/crates/rhodium-pricing/tests/fra.rs
new file mode 100644
index 0000000000000000000000000000000000000000..2b4d5d05b49b196c31729bf536c81c0cf4cc23e3
--- /dev/null
+++ b/crates/rhodium-pricing/tests/fra.rs
@@ -0,0 +1,215 @@
+//! FRA の期首の決済を、同じターム金利で割り戻した差額の手計算で照合する。FRA の恒等式を当てる
+//! 範囲（決済額の中の観測だけ）と、par quote での校正も確かめる。
+
+pub mod support;
+
+use rhodium_core::{
+    Result, Warning,
+    currency::{Asset, Currency},
+    fixing::{FixingStore, ObservationPoint},
+    id::IndexId,
+    money::Money,
+    refdata::{Convention, FraConvention, RefData},
+    time::{Instant, model_time},
+};
+use rhodium_fixtures::{COLLATERAL, futures};
+use rhodium_market::{
+    quote::{QuoteKind, QuoteSet},
+    structure::BlockId,
+};
+use rhodium_pricing::{
+    calibration::calibrate,
+    context::Context,
+    definition::{
+        CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
+    },
+    prepare::{PrepareInput, evaluate, prepare},
+    resolve::CollateralKey,
+};
+use rhodium_products::{
+    ir::{
+        Contract, Product,
+        analysis::{Linearity, linearity},
+    },
+    standard::{Fra, InstrumentTemplate, Position, StandardProduct},
+};
+use rhodium_testkit::{PV_REL, assert_close};
+
+use futures::{AS_OF, FINAL, RATE, START, at, config, curve, market, refdata, reference_parts};
+use support::{USD_ABS, price};
+
+#[test]
+fn fra_settlement_matches_start_discounted_cashflow() -> Result<()> {
+    // FRA は期末のクーポンではなく、同じターム金利で割り戻した差額を期首に受け渡す。
+    let as_of = AS_OF.parse()?;
+    let reference = refdata("0D".parse()?, "16:00:00")?;
+    let config = config("")?;
+    let fixings = FixingStore::default();
+    let context = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &fixings,
+    };
+    let market = market(as_of, 0.0)?;
+    let product = StandardProduct::Fra(Fra::new(
+        Money::new(Currency::USD, 100_000.0),
+        "TERM".parse()?,
+        FINAL.parse()?,
+        0.03,
+        Position::Long,
+    )?);
+    let tau = 7.0 / 365.0;
+    let forward = (RATE * tau).exp_m1() / tau;
+    let df = (-RATE * model_time(as_of, at(FINAL, "16:00")?).years()).exp();
+    let expected = df * 100_000.0 * tau * (forward - 0.03) / (1.0 + tau * forward);
+    assert_close!(
+        price(&product, &market, &context)?.get(Currency::USD),
+        expected,
+        PV_REL,
+        USD_ABS
+    );
+    Ok(())
+}
+
+/// 期首 `start` に受け渡す FRA の決済額 N·τ·(L − K)/(1 + τ'·L) と、同じ支払いに並べた
+/// 通常のタームクーポン N·τ·L₁ の独自契約。L は `fixing` の、L₁ は `coupon` の観測。
+fn fra_with_coupon(
+    reference: &RefData,
+    fixing: &str,
+    coupon: Option<&str>,
+    denominator: f64,
+    start: Instant,
+) -> Result<Product> {
+    let tau = 7.0 / 365.0;
+    let index: IndexId = "TERM".parse()?;
+    Product::new(
+        Contract::build(|builder| {
+            let floating = builder.fixing(&index, ObservationPoint::Date(fixing.parse()?));
+            let strike = builder.constant(0.03);
+            let difference = builder.sub(floating, strike);
+            let accrual = builder.constant(100_000.0 * tau);
+            let numerator = builder.mul([accrual, difference]);
+            let tau = builder.constant(denominator);
+            let interest = builder.mul([tau, floating]);
+            let one = builder.constant(1.0);
+            let discount = builder.add([one, interest]);
+            let mut amount = builder.div(numerator, discount);
+            if let Some(date) = coupon {
+                let rate = builder.fixing(&index, ObservationPoint::Date(date.parse()?));
+                let accrual = builder.constant(100_000.0 * 7.0 / 365.0);
+                let coupon = builder.mul([accrual, rate]);
+                amount = builder.add([amount, coupon]);
+            }
+            Ok(builder.pay(Asset::Cash(Currency::USD), amount, start))
+        })?,
+        reference,
+    )
+}
+
+#[test]
+fn fra_identity_applies_to_its_own_settlement_and_accrual_only() -> Result<()> {
+    // FRA の期首払いの扱い（支払いの遅れを期首から測る）は、その決済額の中の観測にだけ当てる。
+    // 01-14 の支払いに、01-14 に始まる FRA の決済と、01-07〜01-14 の通常のクーポンを並べる。
+    // クーポンは期末の 01-14 に払うので遅れは無く、支払い全体に FRA の扱いを広げると
+    // 期首の 01-07 から 5 営業日の遅れと誤って測り、近似を記録する（上限が小さければ拒否する）。
+    // また分母の τ' が指標の日数係数 τ と違う独自契約は 1/(1 + τ'L) が割引因子にならないので、
+    // 恒等式を当てずリニアでないと判定する。
+    let as_of = AS_OF.parse()?;
+    let reference = refdata("0D".parse()?, "16:00:00")?;
+    let config = config("")?;
+    let fixings = FixingStore::default();
+    let context = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &fixings,
+    };
+    let market = market(as_of, 0.0)?;
+    let start = at("2026-01-14", "16:00")?;
+    let tau = 7.0 / 365.0;
+    let mixed = fra_with_coupon(&reference, "2026-01-14", Some("2026-01-07"), tau, start)?;
+    let plan = prepare(
+        PrepareInput::Ir(&mixed),
+        &COLLATERAL,
+        &[],
+        market.structure(),
+        &context,
+    )?;
+    assert_eq!(plan.warnings(), &[] as &[Warning]);
+    let forward = (RATE * tau).exp_m1() / tau;
+    let df = (-RATE * model_time(as_of, start).years()).exp();
+    let expected = df * 100_000.0 * tau * ((forward - 0.03) / (1.0 + tau * forward) + forward);
+    assert_close!(
+        evaluate(&plan, &market)?.get(Currency::USD),
+        expected,
+        PV_REL,
+        USD_ABS
+    );
+
+    // 分母の日数係数を Act/360 で数えた独自契約は FRA の恒等式に当たらない。
+    let mismatched = fra_with_coupon(&reference, "2026-01-14", None, 7.0 / 360.0, start)?;
+    assert!(matches!(
+        linearity(&mismatched),
+        Linearity::NonLinear(reason) if reason.contains("division")
+    ));
+    // 同じ τ なら、小数で書いた τ の丸めの範囲で恒等式に当たる。
+    let decimal = fra_with_coupon(
+        &reference,
+        "2026-01-14",
+        None,
+        0.019_178_082_191_780_82,
+        start,
+    )?;
+    assert_eq!(linearity(&decimal), Linearity::Linear);
+    Ok(())
+}
+
+#[test]
+fn fra_quote_calibrates_the_trade_to_par() -> Result<()> {
+    // スポットから指定した 0D×1W FRA の par quote で曲線を解き、同一条件の取引を評価する。
+    let as_of = AS_OF.parse()?;
+    let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
+    parts.conventions.insert(
+        "USD-FRA".parse()?,
+        Convention::Fra(FraConvention {
+            index: "TERM".parse()?,
+            spot_lag: "2BD".parse()?,
+        }),
+    );
+    let reference = RefData::new(parts)?;
+    let config = config("")?;
+    let context = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &FixingStore::default(),
+    };
+    let template = InstrumentTemplate::Fra {
+        convention: "USD-FRA".parse()?,
+        start: "0D".parse()?,
+        end: "1W".parse()?,
+    };
+    let quote_id = "USD.FRA.1W".parse()?;
+    let binding = InstrumentBinding::new(
+        quote_id,
+        QuoteKind::Rate,
+        CalibrationTarget::ParRate(template.clone()),
+        QuoteTransform::Identity,
+        CollateralKey::Uncollateralized,
+        ElementSelector::Block(BlockId::nodes("USD".parse()?)),
+    )?;
+    let definition =
+        MarketDefinition::new("USD".parse()?, "Etc/UTC".parse()?, [curve()?], [binding])?;
+    let tau = 7.0 / 365.0;
+    let quote = (RATE * tau).exp_m1() / tau;
+    let calibrated = calibrate(
+        &definition,
+        &QuoteSet::new([("USD.FRA.1W".parse()?, quote)])?,
+        &context,
+    )?;
+    let trade = template.build_for_trade(quote, START.parse()?, &reference)?;
+    let pv = price(&trade, calibrated.market(), &context)?;
+    assert_close!(pv.get(Currency::USD), 0.0_f64, PV_REL, USD_ABS);
+    Ok(())
+}
diff --git a/crates/rhodium-pricing/tests/fx_curves.rs b/crates/rhodium-pricing/tests/fx_curves.rs
index f25bf0631703bd4a511b26661b46b7c9312f132f..46425823f3ff9eb225b5ee614a8ecb4c67c917c2 100644
--- a/crates/rhodium-pricing/tests/fx_curves.rs
+++ b/crates/rhodium-pricing/tests/fx_curves.rs
@@ -26,7 +26,7 @@ use rhodium_pricing::{
         DirectFx, FxDatedQuote, FxDefinition, FxImpliedCurve, FxQuotes, MarketDefinition,
     },
 };
-use rhodium_testkit::assert_close;
+use rhodium_testkit::{MARKET_ABS, MARKET_REL, assert_close};
 
 const BEFORE: &str = "2026-01-02T09:00:00Z";
 const JUST_BEFORE: &str = "2026-01-02T15:59:59Z";
@@ -37,8 +37,6 @@ const USDJPY_SPOT: &str = "2026-01-08T16:00:00Z";
 const USDCAD_SPOT: &str = "2026-01-06T16:00:00Z";
 const USDCAD_SN: &str = "2026-01-07T16:00:00Z";
 const LONG: &str = "2026-03-02T16:00:00Z";
-const REL: f64 = 1e-12;
-const ABS: f64 = 1e-14;
 
 fn clock() -> Result<LocalEventTime> {
     Ok(LocalEventTime {
@@ -288,42 +286,47 @@ fn fx_x0_short_dated_chain() -> Result<()> {
                 assert_close!(
                     jpy_market.forward(AT_CUTOFF.parse()?)?.rate(),
                     149.99,
-                    REL,
-                    ABS
+                    MARKET_REL,
+                    MARKET_ABS
                 );
                 assert_close!(
                     cad_market.forward(AT_CUTOFF.parse()?)?.rate(),
                     1.299,
-                    REL,
-                    ABS
+                    MARKET_REL,
+                    MARKET_ABS
                 );
                 (
                     hand_x0(now, AT_CUTOFF, USDJPY_TOMORROW, 149.99, 150.01)?,
                     hand_x0(now, AT_CUTOFF, USDCAD_SPOT, 1.299, 1.300)?,
                 )
             };
-            assert_close!(jpy_market.today().rate(), jpy_x0, REL, ABS);
-            assert_close!(jpy_other.today().rate(), jpy_x0, REL, ABS);
-            assert_close!(cad_market.today().rate(), cad_x0, REL, ABS);
+            assert_close!(jpy_market.today().rate(), jpy_x0, MARKET_REL, MARKET_ABS);
+            assert_close!(jpy_other.today().rate(), jpy_x0, MARKET_REL, MARKET_ABS);
+            assert_close!(cad_market.today().rate(), cad_x0, MARKET_REL, MARKET_ABS);
             for (at, expected) in [
                 (USDJPY_TOMORROW, 150.01),
                 (USDJPY_SPOT, 150.0),
                 ("2026-01-09T16:00:00Z", 150.03),
                 (LONG, 150.27),
             ] {
-                assert_close!(jpy_market.forward(at.parse()?)?.rate(), expected, REL, ABS);
+                assert_close!(
+                    jpy_market.forward(at.parse()?)?.rate(),
+                    expected,
+                    MARKET_REL,
+                    MARKET_ABS
+                );
             }
             assert_close!(
                 cad_market.forward(USDCAD_SPOT.parse()?)?.rate(),
                 1.300,
-                REL,
-                ABS
+                MARKET_REL,
+                MARKET_ABS
             );
             assert_close!(
                 cad_market.forward(USDCAD_SN.parse()?)?.rate(),
                 1.302,
-                REL,
-                ABS
+                MARKET_REL,
+                MARKET_ABS
             );
             let at = LONG.parse()?;
             let point = result
@@ -333,8 +336,8 @@ fn fx_x0_short_dated_chain() -> Result<()> {
             assert_close!(
                 jpy_other.bind_forward(&point)?.rate(),
                 jpy_x0 * (-0.032 * model_time(now, at).years()).exp(),
-                REL,
-                ABS
+                MARKET_REL,
+                MARKET_ABS
             );
             assert!(matches!(
                 jpy_market.bind_forward(&point),
@@ -396,7 +399,7 @@ fn fx_cross_direct_vs_triangulated() -> Result<()> {
         let usdjpy = view.fx_today("USDJPY".parse()?)?.rate();
         let usdcad = view.fx_today("USDCAD".parse()?)?.rate();
         let triangulated = view.fx_today(pair)?.rate();
-        assert_close!(triangulated, usdjpy / usdcad, REL, ABS);
+        assert_close!(triangulated, usdjpy / usdcad, MARKET_REL, MARKET_ABS);
         view.ccy_scope(
             "CAD".parse()?,
             view.discount_curve(&"CAD".parse()?)?,
@@ -409,8 +412,8 @@ fn fx_cross_direct_vs_triangulated() -> Result<()> {
                         assert_close!(
                             jpy_scope.into_money(value).amount(),
                             10.0 * triangulated,
-                            REL,
-                            ABS
+                            MARKET_REL,
+                            MARKET_ABS
                         );
                         Ok(())
                     },
@@ -443,7 +446,7 @@ fn fx_cross_direct_vs_triangulated() -> Result<()> {
     direct.market().with_view(|view| {
         let time = model_time(BEFORE.parse()?, "2026-01-05T16:00:00Z".parse()?).years();
         let actual = view.fx_today(pair)?.rate();
-        assert_close!(actual, 111.0 * (0.015 * time).exp(), REL, ABS);
+        assert_close!(actual, 111.0 * (0.015 * time).exp(), MARKET_REL, MARKET_ABS);
         assert!((actual - path).abs() > 1.0);
         Ok(())
     })
@@ -474,7 +477,7 @@ fn fx_quote_update_matches_full_recalibration() -> Result<()> {
         assert_eq!(ptr::eq(old.outcome(), new.outcome()), cad);
     }
     for (&actual, &expected) in updated.market().theta().iter().zip(fresh.market().theta()) {
-        assert_close!(actual, expected, REL, ABS);
+        assert_close!(actual, expected, MARKET_REL, MARKET_ABS);
     }
     Ok(())
 }
diff --git a/crates/rhodium-pricing/tests/fx_fixing.rs b/crates/rhodium-pricing/tests/fx_fixing.rs
new file mode 100644
index 0000000000000000000000000000000000000000..9db4e894d0cafae588a8b9349baa4a59b7e4a0fe
--- /dev/null
+++ b/crates/rhodium-pricing/tests/fx_fixing.rs
@@ -0,0 +1,468 @@
+//! FX の fixing、FX フォワードと FX スワップを、両通貨の受渡の割引と手計算のアウトライトで照合する。
+//! fixing の受渡日はペアのスポットラグで決まり、予想はペアの FX の担保区分の下で定まる。
+//! アウトライトの quote での校正（正方系と、価格の単位で残差を測る最小二乗）も確かめる。
+
+pub mod support;
+
+use rhodium_core::{
+    Error, Result, Warning,
+    currency::{Asset, Currency},
+    fixing::{FixingStore, ObservationPoint},
+    fx::{CcyPair, FxRate},
+    id::IndexId,
+    money::Money,
+    quantity::Spread,
+    refdata::FxSpotLag,
+    time::{Instant, model_time},
+};
+use rhodium_fixtures::{COLLATERAL, futures};
+use rhodium_market::{
+    curve::{CurveExpr, CurveRole, CurveSpec, Extrapolation},
+    fx::FxCollateral,
+    market::Market,
+    quote::{QuoteKind, QuoteSet},
+    structure::{BlockId, ParamId},
+};
+use rhodium_pricing::{
+    calibration::{GroupOutcome, calibrate},
+    config::PricingConfig,
+    context::Context,
+    definition::{
+        CalibrationTarget, DirectFx, ElementSelector, FxDefinition, FxQuotes, InstrumentBinding,
+        MarketDefinition, QuoteTransform,
+    },
+    prepare::{PrepareInput, evaluate, prepare},
+    resolve::CollateralKey,
+};
+use rhodium_products::{
+    agreement::{CashCollateralTerms, Remuneration},
+    ir::{Contract, Product},
+    standard::{FxForward, FxSwap, InstrumentTemplate, Position, StandardProduct},
+};
+use rhodium_testkit::{MARKET_ABS, MARKET_REL, PV_REL, assert_close};
+
+use futures::{AS_OF, PILLAR, RATE, START, at, curve};
+use support::{
+    JPY_ABS, USD_ABS, fx_config, fx_market, fx_market_with, fx_reference, fx_reference_with_lag,
+    price,
+};
+
+#[test]
+fn fx_fixing_spot_value_uses_the_pair_spot_lag() -> Result<()> {
+    // 受渡日の指定の無い FX fixing は、fixing 日からペアのスポットラグだけ進めた受渡日のレートを観測する。
+    // 01-12（月）の fixing の受渡日は、T+1 では 01-13 の祝日を飛ばして 01-14、T+2 では 01-15。
+    // 支払いは 01-20 で、その円の受渡を P_JPY で割り引く。F(v) = 150·exp(−(r_USD − r_JPY)·t_v)。
+    let as_of = AS_OF.parse()?;
+    let market = fx_market(as_of)?;
+    let config = fx_config("")?;
+    let fixings = FixingStore::default();
+    let pay = at("2026-01-20", "16:00")?;
+    let fx: IndexId = "FX-FIX".parse()?;
+    let discount = |rate: f64, at| (-rate * model_time(as_of, at).years()).exp();
+    let forward = |at| 150.0 * discount(RATE - 0.01, at);
+    let mut values = Vec::new();
+    for (spot_lag, value_date) in [
+        (FxSpotLag::OneBusinessDay, "2026-01-14"),
+        (FxSpotLag::TwoBusinessDays, "2026-01-15"),
+    ] {
+        let reference = fx_reference_with_lag(spot_lag)?;
+        let product = Product::new(
+            Contract::build(|builder| {
+                let rate = builder.fixing(&fx, ObservationPoint::Date("2026-01-12".parse()?));
+                Ok(builder.pay(Asset::Cash(Currency::JPY), rate, pay))
+            })?,
+            &reference,
+        )?;
+        let ctx = Context {
+            as_of,
+            refdata: &reference,
+            config: &config,
+            fixings: &fixings,
+        };
+        let plan = prepare(
+            PrepareInput::Ir(&product),
+            &COLLATERAL,
+            &[],
+            market.structure(),
+            &ctx,
+        )?;
+        // 受渡日と支払日が違うので、タイミング調整を省いた近似を記録する。警告の一覧全体を比べ、
+        // この近似だけを記録することを確かめる。
+        assert_eq!(
+            plan.warnings(),
+            [Warning::PaymentDelayIgnored],
+            "{spot_lag:?}"
+        );
+        let pv = evaluate(&plan, &market)?.get(Currency::JPY);
+        assert_close!(
+            pv,
+            forward(at(value_date, "16:00")?) * discount(0.01, pay),
+            PV_REL,
+            JPY_ABS
+        );
+        values.push(pv);
+    }
+    assert!((values[0] - values[1]).abs() > 1e-6);
+    Ok(())
+}
+
+#[test]
+fn fx_forward_and_swap_exchange_both_currencies() -> Result<()> {
+    // 両通貨の各受渡を割り引き、近い受渡と遠い受渡の向きが逆になることを確かめる。
+    let as_of = AS_OF.parse()?;
+    let market = fx_market(as_of)?;
+    let pair: CcyPair = "USDJPY".parse()?;
+    let reference = fx_reference(pair)?;
+    let config = fx_config("")?;
+    let context = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &FixingStore::default(),
+    };
+    let near = at("2026-01-07", "16:00")?;
+    let far = at("2026-01-21", "16:00")?;
+    let notional = Money::new(Currency::USD, 1000.0);
+    let first = FxRate::forward(pair, 151.0, near)?;
+    let second = FxRate::forward(pair, 152.0, far)?;
+    let forward = StandardProduct::FxForward(FxForward::new(notional, first, Position::Long)?);
+    let swap = StandardProduct::FxSwap(FxSwap::new(notional, first, second, Position::Long)?);
+    let df = |rate: f64, t| (-rate * model_time(as_of, t).years()).exp();
+
+    let pv = price(&forward, &market, &context)?;
+    assert_close!(
+        pv.get(Currency::USD),
+        1000.0 * df(RATE, near),
+        PV_REL,
+        USD_ABS
+    );
+    assert_close!(
+        pv.get(Currency::JPY),
+        -151_000.0 * df(0.01, near),
+        PV_REL,
+        JPY_ABS
+    );
+    let pv = price(&swap, &market, &context)?;
+    assert_close!(
+        pv.get(Currency::USD),
+        1000.0 * (df(RATE, near) - df(RATE, far)),
+        PV_REL,
+        USD_ABS
+    );
+    assert_close!(
+        pv.get(Currency::JPY),
+        -151_000.0 * df(0.01, near) + 152_000.0 * df(0.01, far),
+        PV_REL,
+        JPY_ABS
+    );
+
+    // 市場の二つのアウトライトで約定すると、各受渡の報告通貨建て PV がゼロになる。
+    let outright = |at| FxRate::forward(pair, 150.0 * df(RATE, at) / df(0.01, at), at);
+    let par = StandardProduct::FxSwap(FxSwap::new(
+        notional,
+        outright(near)?,
+        outright(far)?,
+        Position::Long,
+    )?);
+    let pv = price(&par, &market, &context)?;
+    assert_close!(
+        pv.get(Currency::JPY) + 150.0 * pv.get(Currency::USD),
+        0.0_f64,
+        PV_REL,
+        JPY_ABS
+    );
+    Ok(())
+}
+
+#[test]
+fn fx_outright_quote_calibrates_a_par_trade_in_both_currencies() -> Result<()> {
+    // outright の校正は base と quote の両方の PV を評価時点の FX で合算する。
+    let as_of = AS_OF.parse()?;
+    let pair: CcyPair = "USDJPY".parse()?;
+    let reference = fx_reference(pair)?;
+    let config = fx_config("")?;
+    let context = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &FixingStore::default(),
+    };
+    // 一つのアウトライトに許容幅を付けずに当てる正方系。
+    let definition = outright_definition(&[("USDJPY.FWD.1W", None)])?;
+    let quote = 150.0 * (-(RATE - 0.01) * 7.0 / 365.0).exp();
+    let calibrated = calibrate(
+        &definition,
+        &QuoteSet::new([
+            (
+                "USD.DF.1".parse()?,
+                -RATE * model_time(as_of, PILLAR.parse()?).years(),
+            ),
+            ("USDJPY.FWD.1W".parse()?, quote),
+            ("USDJPY.SPOT".parse()?, 150.0),
+        ])?,
+        &context,
+    )?;
+    let trade = outright_template()?.build_for_trade(quote, START.parse()?, &reference)?;
+    let market = calibrated.market();
+    let pv = price(&trade, market, &context)?;
+    let fx = market.with_view(|view| view.fx_today(pair))?;
+    assert_close!(
+        pv.convert(Currency::JPY, &[fx])?.amount(),
+        0.0_f64,
+        PV_REL,
+        JPY_ABS
+    );
+    Ok(())
+}
+
+/// 慣行 `USDJPY-FWD` の 1W の USDJPY アウトライトのテンプレート。
+fn outright_template() -> Result<InstrumentTemplate> {
+    Ok(InstrumentTemplate::FxForward {
+        convention: "USDJPY-FWD".parse()?,
+        tenor: "1W".parse()?,
+    })
+}
+
+/// USD のピラーを直接マークで、JPY のピラーを 1W の USDJPY アウトライトで当てる市場の定義。
+/// `outrights` は JPY の曲線に当てるアウトライトの quote の ID と許容幅 σ（無ければ正方系）。
+fn outright_definition(outrights: &[(&str, Option<f64>)]) -> Result<MarketDefinition> {
+    let pair: CcyPair = "USDJPY".parse()?;
+    let template = outright_template()?;
+    let discount = CurveSpec::new(
+        "JPY".parse()?,
+        Currency::JPY,
+        [CurveRole::Discount],
+        [PILLAR.parse()?],
+        Extrapolation::Reject,
+    )?;
+    let usd = InstrumentBinding::new(
+        "USD.DF.1".parse()?,
+        QuoteKind::Mark,
+        CalibrationTarget::DirectMark,
+        QuoteTransform::Identity,
+        CollateralKey::Uncollateralized,
+        ElementSelector::Param(ParamId::new(BlockId::nodes("USD".parse()?), 0)),
+    )?;
+    let jpy = outrights
+        .iter()
+        .map(|&(quote, sigma)| {
+            let binding = InstrumentBinding::new(
+                quote.parse()?,
+                QuoteKind::Price,
+                CalibrationTarget::Price(template.clone()),
+                QuoteTransform::Identity,
+                CollateralKey::Uncollateralized,
+                ElementSelector::Block(BlockId::nodes("JPY".parse()?)),
+            )?;
+            sigma.map_or(Ok(binding.clone()), |sigma| binding.with_tolerance(sigma))
+        })
+        .collect::<Result<Vec<_>>>()?;
+    MarketDefinition::new(
+        "USD".parse()?,
+        "Etc/UTC".parse()?,
+        [curve()?, discount],
+        std::iter::once(usd).chain(jpy),
+    )?
+    .with_fx([FxDefinition {
+        pair,
+        direct: Some(DirectFx {
+            market_collateral: "csa".parse()?,
+            quotes: FxQuotes {
+                spot: "USDJPY.SPOT".parse()?,
+                on: None,
+                tn: None,
+                sn: None,
+                forwards: vec![],
+            },
+            implied: None,
+        }),
+        collaterals: vec![FxCollateral {
+            id: "csa".parse()?,
+            base: "USD".parse()?,
+            quote: "JPY".parse()?,
+        }],
+        via: None,
+    }])
+}
+
+#[test]
+fn fx_outright_least_squares_measures_residuals_in_price_units() -> Result<()> {
+    // アウトライトの許容幅 σ は価格の単位なので、最小二乗の残差も価格の差 F − q で測る。
+    // 同じ受渡の二つのアウトライト q₁、q₂ = q₁ + 2 を σ₁ = 0.01、σ₂ = 1 で当てると、
+    // F* = (q₁/σ₁² + q₂/σ₂²)/(1/σ₁² + 1/σ₂²) = q₁ + 2/10001 で、残差は F* − q₁ ≈ 0.0002
+    // （σ₁ の 0.02 倍で内）と F* − q₂ ≈ −2（σ₂ の 2 倍で外）。外れた q₂ だけを警告する。
+    // 残差を ln(F/q) のまま 1/σ で重みを付けると、q₂ の外れは 0.013 と測られて警告されない。
+    // quote デルタ dθ/dq は、残差の quote 微分を −1 にした随伴が再校正の中心差分と一致する。
+    let as_of = AS_OF.parse()?;
+    let reference = fx_reference("USDJPY".parse()?)?;
+    let config = fx_config("")?;
+    let context = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &FixingStore::default(),
+    };
+    let definition = outright_definition(&[
+        ("USDJPY.FWD.1W.A", Some(0.01)),
+        ("USDJPY.FWD.1W.B", Some(1.0)),
+    ])?;
+    let q1 = 150.0 * (-(RATE - 0.01) * 7.0 / 365.0).exp();
+    let q2 = q1 + 2.0;
+    let quotes = |q1: f64| {
+        QuoteSet::new([
+            (
+                "USD.DF.1".parse()?,
+                -RATE * model_time(as_of, PILLAR.parse()?).years(),
+            ),
+            ("USDJPY.FWD.1W.A".parse()?, q1),
+            ("USDJPY.FWD.1W.B".parse()?, q2),
+            ("USDJPY.SPOT".parse()?, 150.0),
+        ])
+    };
+    let calibrated = calibrate(&definition, &quotes(q1)?, &context)?;
+    let fitted = (q1 * 1e4 + q2) / (1e4 + 1.0);
+    // F* で約定したフォワードは、報告通貨で合算した PV がゼロ。
+    let trade = outright_template()?.build_for_trade(fitted, START.parse()?, &reference)?;
+    let market = calibrated.market();
+    let fx = market.with_view(|view| view.fx_today("USDJPY".parse()?))?;
+    assert_close!(
+        price(&trade, market, &context)?
+            .convert(Currency::JPY, &[fx])?
+            .amount(),
+        0.0_f64,
+        PV_REL,
+        JPY_ABS
+    );
+    let jpy = BlockId::nodes("JPY".parse()?);
+    let group = calibrated
+        .groups()
+        .iter()
+        .position(|group| group.blocks().contains(&jpy))
+        .ok_or_else(|| Error::Internal("missing JPY group".into()))?;
+    let GroupOutcome::Solved(solution) = calibrated.groups()[group].outcome() else {
+        return Err(Error::Calibration("JPY group not solved".into()));
+    };
+    assert_close!(solution.residuals()[0], fitted - q1, 1e-9_f64, 1e-9_f64);
+    assert_close!(solution.residuals()[1], fitted - q2, 1e-9_f64, 1e-9_f64);
+    // 警告の一覧全体を比べ、外れた q₂ の警告だけがあり、q₁ の警告も他の警告も無いことを確かめる。
+    assert_eq!(
+        calibrated.warnings(),
+        [Warning::CalibrationFitOutsideTolerance {
+            quote: "USDJPY.FWD.1W.B".parse()?
+        }]
+    );
+
+    // dθ/dq₁ = −(∂G/∂θ)⁻¹·∂G/∂q₁（列は quote の ID の昇順で A が先）を、q₁ を動かした再校正と比べる。
+    let inverse = solution
+        .dg_dtheta()
+        .solve(&[1.0][..])
+        .map_err(|error| Error::Numerical(error.to_string()))?
+        .as_slice()[0];
+    let dg_dq1 = solution
+        .dg_dq()
+        .get(0, 0)
+        .map_err(|error| Error::Numerical(error.to_string()))?;
+    let position = calibrated
+        .market()
+        .structure()
+        .layout()
+        .range(&jpy)
+        .ok_or_else(|| Error::Internal("missing JPY block".into()))?
+        .start;
+    let bump = 1e-3;
+    let theta = |q1: f64| -> Result<f64> {
+        Ok(calibrate(&definition, &quotes(q1)?, &context)?
+            .market()
+            .theta()[position])
+    };
+    let central = (theta(q1 + bump)? - theta(q1 - bump)?) / (2.0 * bump);
+    assert_close!(-inverse * dg_dq1, central, 1e-6_f64, 1e-12_f64);
+    Ok(())
+}
+
+/// [`fx_market`] に、JPY を 25bp 上乗せした割引曲線 `JPY-25` を足した市場。
+fn fx_market_with_spread_curve(as_of: Instant) -> Result<Market<f64>> {
+    let spread = CurveSpec::from_expr(
+        "JPY-25".parse()?,
+        Currency::JPY,
+        [CurveRole::Discount],
+        CurveExpr::Sum(vec![
+            CurveExpr::Ref("JPY".parse()?),
+            CurveExpr::Flat(0.0025),
+        ]),
+    )?;
+    fx_market_with(as_of, [spread])
+}
+
+#[test]
+fn fx_fixing_under_another_collateral_records_cross_collateral_projection() -> Result<()> {
+    // FX の予想は、ペアの FX の担保区分（USD と JPY の曲線の組）の下で定まる。支払いを別の担保
+    // （JPY を 25bp 上乗せした曲線で割り引く CSA）で受け取ると、予想と割引の担保が違うので
+    // CrossCollateralProjection を記録し、禁止すると準備できない。
+    let as_of = AS_OF.parse()?;
+    let reference = fx_reference("USDJPY".parse()?)?;
+    let fixings = FixingStore::default();
+    let market = fx_market_with_spread_curve(as_of)?;
+    let product = Product::new(
+        Contract::build(|builder| {
+            let rate = builder.fixing(
+                &"FX-FIX".parse()?,
+                ObservationPoint::Date("2026-01-12".parse()?),
+            );
+            Ok(builder.pay(Asset::Cash(Currency::JPY), rate, at("2026-01-14", "16:00")?))
+        })?,
+        &reference,
+    )?;
+    let collateral = CollateralKey::Cash(CashCollateralTerms {
+        ccy: Currency::JPY,
+        remuneration: Remuneration {
+            index: "JPY-TERM".parse()?,
+            spread: Spread::new(0.0025)?,
+            floor: None,
+        },
+    });
+    let rules = r#"
+[[discounting]]
+currency = "JPY"
+collateral = { cash = { ccy = "JPY", index = "JPY-TERM", spread = "25bp" } }
+curve = "JPY-25"
+[[discounting]]
+currency = "USD"
+collateral = { cash = { ccy = "JPY", index = "JPY-TERM", spread = "25bp" } }
+curve = "USD"
+"#;
+    let prepare_under = |config: &PricingConfig, collateral: &CollateralKey| {
+        prepare(
+            PrepareInput::Ir(&product),
+            collateral,
+            &[],
+            market.structure(),
+            &Context {
+                as_of,
+                refdata: &reference,
+                config,
+                fixings: &fixings,
+            },
+        )
+    };
+    // 01-12 の fixing の受渡は T+2 の 01-14 の 16 時で、支払いと同じ瞬間なので、支払いの遅れの近似は
+    // 無い。警告の一覧全体を比べ、担保の違いの近似だけが増減することを確かめる。
+    let native = prepare_under(&fx_config("")?, &CollateralKey::Uncollateralized)?;
+    assert_eq!(native.warnings(), []);
+    let spread_config = fx_config(rules)?;
+    let foreign = prepare_under(&spread_config, &collateral)?;
+    assert_eq!(foreign.warnings(), [Warning::CrossCollateralProjection]);
+    // 予想は担保区分の曲線のフォワードのまま、割引だけが 25bp の曲線に変わる。
+    let ratio = evaluate(&foreign, &market)?.get(Currency::JPY)
+        / evaluate(&native, &market)?.get(Currency::JPY);
+    let pay = model_time(as_of, at("2026-01-14", "16:00")?).years();
+    assert_close!(ratio, (-0.0025 * pay).exp(), MARKET_REL, MARKET_ABS);
+    let strict = fx_config(&format!(
+        "{rules}[conventions]\ndisabled = ['cross_collateral_projection']"
+    ))?;
+    assert!(matches!(
+        prepare_under(&strict, &collateral),
+        Err(Error::Unsupported(_))
+    ));
+    Ok(())
+}
diff --git a/crates/rhodium-pricing/tests/least_squares_calibration.rs b/crates/rhodium-pricing/tests/least_squares_calibration.rs
index d51cc8c9883a36d5071f1adb4900e3e8598d3d6a..b8f75a84d27a42ff4bd8ffcfb523f27b4e6fd560 100644
--- a/crates/rhodium-pricing/tests/least_squares_calibration.rs
+++ b/crates/rhodium-pricing/tests/least_squares_calibration.rs
@@ -1,13 +1,15 @@
 //! 最小二乗の校正グループ: 手計算の検算、当てはめの許容幅の警告、warm start、定義の検査、
 //! 罰則付きの会合日の階段（R4）と Gauss–Newton の記録。
 
-use std::sync::Arc;
+use std::{collections::BTreeSet, sync::Arc};
 
 use rhodium_core::{Error, Result, Warning, currency::Currency};
+use rhodium_fixtures::COLLATERAL;
+use rhodium_fixtures::{ois, steps};
 use rhodium_market::{
     curve::{CurveRole, CurveSpec, Extrapolation},
     market::Market,
-    quote::{QuoteId, QuoteKind, QuoteSet},
+    quote::{QuoteKind, QuoteSet},
     structure::{BlockId, ParamId},
 };
 use rhodium_pricing::{
@@ -18,28 +20,12 @@ use rhodium_pricing::{
         StepPenalty,
     },
     prepare::evaluate,
-    resolve::CollateralKey,
 };
 use rhodium_testkit::assert_close;
 
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、最小二乗の検査に使わない関数がある"
-)]
-#[path = "fixtures/ois.rs"]
-pub mod ois;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、最小二乗の検査に使わない関数がある"
-)]
-#[path = "fixtures/steps.rs"]
-pub mod steps;
-
 use ois::{Case, Inputs, cases, midnight, prepare_ois};
 use steps::{OIS_TOLERANCE, PENALTY_SCALE, PENALTY_WEIGHT, PenalizedTona};
 
-const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
 const MARKS: [&str; 2] = ["MARK.1", "MARK.2"];
 
 /// 一つのピラーの曲線の ln DF を、許容幅 σ の二つの直接マークで当てる市場の定義。
@@ -89,14 +75,12 @@ fn marks(values: [f64; 2]) -> Result<QuoteSet> {
     )
 }
 
-#[expect(
-    clippy::wildcard_enum_match_arm,
-    reason = "#196 で各 variant を網羅する"
-)]
 fn solution(market: &CalibratedMarket, index: usize) -> Result<&GroupSolution> {
     match market.groups()[index].outcome() {
         GroupOutcome::Solved(solution) => Ok(solution),
-        other => Err(Error::Calibration(format!("{other:?}"))),
+        other @ (GroupOutcome::Failed(_) | GroupOutcome::Unavailable) => {
+            Err(Error::Calibration(format!("{other:?}")))
+        }
     }
 }
 
@@ -140,10 +124,6 @@ fn weighted_direct_marks_solve_the_hand_least_squares() -> Result<()> {
 }
 
 #[test]
-#[expect(
-    clippy::wildcard_enum_match_arm,
-    reason = "#196 で各 variant を網羅する"
-)]
 fn fit_outside_tolerance_warns_per_quote() -> Result<()> {
     // σ₁ = 1bp、σ₂ = 10bp で q₁ = −1%、q₂ = −2% を当てると、θ = (w₁²q₁ + w₂²q₂)/(w₁² + w₂²)
     // = −1.0099% で、残差は F₁ = −0.99bp（σ₁ の内）、F₂ = 99bp（σ₂ の外）。外れた quote だけを警告する。
@@ -153,22 +133,15 @@ fn fit_outside_tolerance_warns_per_quote() -> Result<()> {
     let ctx = inputs.context();
     let definition = direct_marks(&tona, [1e-4, 1e-3])?;
     let result = calibrate(&definition, &marks([-0.01, -0.02])?, &ctx)?;
-    let outside: Vec<_> = result
-        .warnings()
-        .iter()
-        .filter_map(|warning| match warning {
-            Warning::CalibrationFitOutsideTolerance { quote } => Some(quote.as_str()),
-            _ => None,
-        })
-        .collect();
-    assert_eq!(outside, ["MARK.2"]);
-    let consistent = result.update(&marks([-0.01, -0.01])?)?;
-    assert!(
-        !consistent
-            .warnings()
-            .iter()
-            .any(|warning| matches!(warning, Warning::CalibrationFitOutsideTolerance { .. }))
+    // 警告の一覧全体を比べ、外れた quote の警告だけがあり、他の警告も無いことを確かめる。
+    assert_eq!(
+        result.warnings(),
+        [Warning::CalibrationFitOutsideTolerance {
+            quote: "MARK.2".parse()?
+        }]
     );
+    let consistent = result.update(&marks([-0.01, -0.01])?)?;
+    assert_eq!(consistent.warnings(), []);
     Ok(())
 }
 
@@ -432,10 +405,6 @@ fn penalized_meeting_steps_fit_standard_tenor_ois() -> Result<()> {
 }
 
 #[test]
-#[expect(
-    clippy::wildcard_enum_match_arm,
-    reason = "#196 で各 variant を網羅する"
-)]
 fn gauss_newton_records_are_flagged_and_differ_from_exact() -> Result<()> {
     // Gauss–Newton の記録は残差の 2 階の項を省くので、残差の残る階段のグループでは ∂G/∂θ が
     // 厳密な記録と違う。解（θ）は同じで、近似の記録を解と警告の両方で報告する。
@@ -451,22 +420,26 @@ fn gauss_newton_records_are_flagged_and_differ_from_exact() -> Result<()> {
     }
     let gauss_newton = solution(&approximate, 1)?;
     assert_eq!(gauss_newton.hessian(), Some(HessianMode::GaussNewton));
-    let flagged: Vec<QuoteId> = approximate
-        .warnings()
-        .iter()
-        .filter_map(|warning| match warning {
-            Warning::CalibrationGaussNewton { quote } => Some(quote.clone()),
-            _ => None,
-        })
-        .collect();
-    let expected: Vec<QuoteId> = tona
+    // Gauss–Newton のグループのすべての OIS の quote に近似の警告が付く。厳密な記録の市場の警告を
+    // 基準にして、増えた警告がちょうどその集合であることを確かめる。
+    let gauss_newton_warnings = tona
         .case
         .nodes
         .iter()
-        .map(|node| tona.case.quote(node))
-        .collect::<Result<_>>()?;
-    assert_eq!(flagged.len(), expected.len());
-    assert!(expected.iter().all(|quote| flagged.contains(quote)));
+        .map(|node| {
+            Ok(Warning::CalibrationGaussNewton {
+                quote: tona.case.quote(node)?,
+            })
+        })
+        .collect::<Result<BTreeSet<_>>>()?;
+    let exact_warnings: BTreeSet<_> = exact.warnings().iter().cloned().collect();
+    let added: BTreeSet<_> = approximate
+        .warnings()
+        .iter()
+        .filter(|warning| !exact_warnings.contains(*warning))
+        .cloned()
+        .collect();
+    assert_eq!(added, gauss_newton_warnings);
     // 同じ右辺の逆作用素の値が、2 階の項の分だけ違う。
     let rhs = vec![1.0; 5];
     let [exact_inverse, approximate_inverse] = [solution(&exact, 1)?, gauss_newton].map(|s| {
diff --git a/crates/rhodium-pricing/tests/linear_pricing.rs b/crates/rhodium-pricing/tests/linear_pricing.rs
index be9ccd056c07acc63a8b65e2b58a2166aaffd651..12048ff92b406be7ddb6ef452a03182a2e6b959a 100644
--- a/crates/rhodium-pricing/tests/linear_pricing.rs
+++ b/crates/rhodium-pricing/tests/linear_pricing.rs
@@ -1,5 +1,7 @@
 //! 標準商品と検証済み IR を、本番の準備・結合・評価の入口で照合する。
 
+pub mod support;
+
 use std::{num::NonZeroU8, sync::Arc};
 
 use rhodium_core::{
@@ -31,8 +33,9 @@ use rhodium_core::{
         ScheduleSpec, SettlementContext, StartRule, StubLength, StubSpec, generate,
     },
     tenor::{DateOffset, Frequency, TenorUnit},
-    time::{Instant, LocalEventTime, UtcOffset, ZoneRules},
+    time::{Instant, LocalEventTime, UtcOffset, ZoneRules, model_time},
 };
+use rhodium_fixtures::futures;
 use rhodium_market::{
     curve::{CurveRole, CurveSpec, Extrapolation},
     market::Market,
@@ -49,13 +52,15 @@ use rhodium_products::{
     standard::{Coupon, Leg, NotionalExchange, PayReceive, PaymentRounding, StandardProduct, Swap},
     trade::{FallbackSelection, FallbackTerms, ProductSpec, Trade},
 };
-use rhodium_testkit::assert_close;
+use rhodium_testkit::{PV_REL, assert_close};
 use scioxide_ad::Tape;
 use scioxide_core::{ExecutionContext, InvalidInputReason, SciOxideError, argument};
 
+use support::{
+    JPY_ABS, USD_ABS, explicit_schedule, fx_config, fx_market, fx_reference, price, schedule, term,
+};
+
 const INDEX: &str = "JPY-TONA";
-const REL: f64 = 1e-12;
-const JPY_ABS: f64 = 1e-6;
 const NOTIONAL: f64 = 1_000_000.0;
 const DISABLE_DAILY_FORWARD_PRODUCT: &str = r#"disabled = ["daily_forward_product"]"#;
 
@@ -213,28 +218,6 @@ fn structure(
     Ok(Arc::new(MarketStructure::new(as_of, curves)?))
 }
 
-fn schedule(dates: &[Date], pays: &[Date]) -> Result<ScheduleSpec> {
-    let periods = dates
-        .windows(2)
-        .zip(pays)
-        .map(|(dates, &payment)| {
-            Ok(ExplicitPeriod {
-                accrual: DateInterval::new(dates[0], dates[1])?,
-                payment,
-                reset: None,
-                settlement: None,
-                references: Vec::new(),
-            })
-        })
-        .collect::<Result<_>>()?;
-    Ok(ScheduleSpec::Explicit(ExplicitPeriods::new(
-        periods,
-        ExplicitDatePolicy::KeepAsConfirmed,
-        PeriodCombination::Sum,
-        None,
-    )?))
-}
-
 fn rule(
     method: ObservationMethod,
     form: CompoundingForm,
@@ -391,7 +374,7 @@ fn compare_quantlib(currency: Currency) -> Result<()> {
     let scale = fixture.fixed_pv.abs() + fixture.floating_pv.abs();
     let market = Market::from_theta(structure.clone(), theta)?;
     let pv = evaluate(&prepared, &market)?.get(currency);
-    assert_close!(pv, fixture.pv, REL, abs.max(REL * scale));
+    assert_close!(pv, fixture.pv, PV_REL, abs.max(PV_REL * scale));
 
     // 固定金利が 0 の同じ契約で変動脚の PV を照合し、差から固定脚も確かめる。
     let float_only = swap(
@@ -405,8 +388,8 @@ fn compare_quantlib(currency: Currency) -> Result<()> {
     )?;
     let float_only = prepare_standard(&float_only, fixture.as_of, &structure, &context)?;
     let floating = evaluate(&float_only, &market)?.get(currency);
-    assert_close!(floating, fixture.floating_pv, REL, abs);
-    assert_close!(pv - floating, fixture.fixed_pv, REL, abs);
+    assert_close!(floating, fixture.floating_pv, PV_REL, abs);
+    assert_close!(pv - floating, fixture.fixed_pv, PV_REL, abs);
     assert_eq!(prepared.method(), Method::Linear);
     assert_eq!(
         prepared.blocks().len(),
@@ -570,12 +553,12 @@ fn fixed_prefix_and_intraday_forward_use_original_coupon_weights()
     let theta = -0.0018;
     let market = Market::from_theta(structure.clone(), vec![theta])?;
     let pv = evaluate(&prepared, &market)?.get(Currency::JPY);
-    assert_close!(pv, expected(theta), REL, JPY_ABS);
+    assert_close!(pv, expected(theta), PV_REL, JPY_ABS);
 
     // 同じ準備で θ を変えても、金曜の因子と月曜以降の比の両方が価格と感応度に反映される。
     let theta = -0.0042;
     let (pv, gradient) = jpy_pv_and_gradient(&structure, &prepared, theta)?;
-    assert_close!(pv, expected(theta), REL, JPY_ABS);
+    assert_close!(pv, expected(theta), PV_REL, JPY_ABS);
     // 金曜の因子と月曜の DF の積は q + (1 − q)·DF(月曜)。q は元の重みと残存期間の比。
     let derivative = NOTIONAL
         * (fixed_factor * (1.0 - full_to_remaining) * monday_time * (theta * monday_time).exp()
@@ -658,7 +641,7 @@ fn partially_fixed_excess_keeps_the_published_digits() -> Result<()> {
         // θ = 0 では割引因子が 1 なので、PV は丸めた受渡額そのもの。
         let pv =
             evaluate(&prepared, &Market::from_theta(structure, vec![0.0])?)?.get(Currency::JPY);
-        assert_close!(pv, 1002.0, REL, JPY_ABS);
+        assert_close!(pv, 1002.0, PV_REL, JPY_ABS);
     }
     Ok(())
 }
@@ -742,7 +725,7 @@ fn growth_and_excess_of_one_compounding_share_its_plan() -> Result<()> {
     let pv =
         evaluate(&prepared, &Market::from_theta(market_structure, vec![0.0])?)?.get(Currency::JPY);
     let excess = 0.005 / 365.0;
-    assert_close!(pv, n * excess + m * (1.0 + excess), REL, JPY_ABS);
+    assert_close!(pv, n * excess + m * (1.0 + excess), PV_REL, JPY_ABS);
     Ok(())
 }
 
@@ -833,7 +816,7 @@ fn explicit_long_stub_matches_generated_rfr_coupon() -> Result<()> {
         Ok(evaluate(&prepared, &market)?.get(Currency::JPY))
     };
     let generated = pv(&generated)?;
-    assert_close!(pv(&one_coupon)?, generated, REL, JPY_ABS);
+    assert_close!(pv(&one_coupon)?, generated, PV_REL, JPY_ABS);
     // 分けた二つの期間は別のクーポンで、差は交差項の数百円になる。
     let difference = generated - pv(&split)?;
     assert!(difference > 100.0, "{difference}");
@@ -936,7 +919,7 @@ fn intraday_forward_fraction_follows_dst_day_length() -> Result<()> {
         let expected = NOTIONAL * 3.0 / 360.0 * rate * discount;
         let market = Market::from_theta(structure, vec![theta])?;
         let pv = evaluate(&prepared, &market)?.get(Currency::USD);
-        assert_close!(pv, expected, REL, usd_abs);
+        assert_close!(pv, expected, PV_REL, usd_abs);
     }
     Ok(())
 }
@@ -980,7 +963,7 @@ fn fixing_publication_boundary_uses_estimate_before_and_published_at_or_after()
             &Market::from_theta(structure.clone(), vec![0.0])?,
         )?
         .get(Currency::JPY);
-        assert_close!(pv, NOTIONAL * rate / 365.0, REL, JPY_ABS);
+        assert_close!(pv, NOTIONAL * rate / 365.0, PV_REL, JPY_ABS);
         if rate > 0.02 {
             let missing = Context {
                 fixings: &only_estimate,
@@ -1057,7 +1040,7 @@ fn unpublished_term_fixing_uses_the_forward_over_the_index_deposit_period() -> R
     assert_close!(
         pv,
         NOTIONAL * (-0.04_f64).exp() * 0.02_f64.exp_m1(),
-        REL,
+        PV_REL,
         JPY_ABS
     );
     Ok(())
@@ -1162,7 +1145,7 @@ fn fallback_pv_compounds_a_published_day_with_a_projected_day() -> Result<()> {
     assert_close!(
         evaluate(&prepared, &market)?.get(ccy),
         expected,
-        REL,
+        PV_REL,
         JPY_ABS
     );
     Ok(())
@@ -1201,7 +1184,7 @@ fn settlement_boundary_and_binding_failures() -> Result<()> {
         )?;
         let market = Market::from_theta(market_structure.clone(), vec![-0.02])?;
         let pv = market.with_view(|view| prepared.bind(&view)?.evaluate())?;
-        assert_close!(pv.get(Currency::JPY), expected, REL, JPY_ABS);
+        assert_close!(pv.get(Currency::JPY), expected, PV_REL, JPY_ABS);
         assert_eq!(pv.terms().count(), currencies);
         assert_eq!(prepared.blocks(), []);
 
@@ -1260,13 +1243,13 @@ fn prepared_swap_reuses_structure_and_preserves_aad_gradient()
     assert_close!(
         evaluate(&prepared, &market)?.get(Currency::JPY),
         expected(-0.002),
-        REL,
+        PV_REL,
         JPY_ABS
     );
 
     let theta: f64 = -0.004;
     let (pv, gradient) = jpy_pv_and_gradient(&structure, &prepared, theta)?;
-    assert_close!(pv, expected(theta), REL, JPY_ABS);
+    assert_close!(pv, expected(theta), PV_REL, JPY_ABS);
     let derivative = NOTIONAL * (linear_df(as_of, end, theta, start) / 7.0 - factor * theta.exp());
     assert_close!(gradient, derivative, 1e-10, JPY_ABS);
     Ok(())
@@ -1312,7 +1295,7 @@ fn lookback_daily_product_preserves_weights_and_checks_policy() -> Result<()> {
     assert_close!(
         evaluate(&prepared, &market)?.get(Currency::JPY),
         expected,
-        REL,
+        PV_REL,
         JPY_ABS
     );
 
@@ -1409,7 +1392,7 @@ fn plain_holiday_start_and_lockout_preserve_daily_weights() -> Result<()> {
         assert_close!(
             evaluate(&prepared, &market)?.get(Currency::JPY),
             expected,
-            REL,
+            PV_REL,
             JPY_ABS
         );
 
@@ -1478,7 +1461,7 @@ fn observation_shift_telescopes_without_replacing_coupon_accrual() -> Result<()>
     let theta: f64 = -0.0018;
     let expected = NOTIONAL * 3.0 * (-theta / 6.0).exp_m1() * theta.exp();
     let pv = evaluate(&prepared, &Market::from_theta(structure, vec![theta])?)?.get(Currency::JPY);
-    assert_close!(pv, expected, REL, JPY_ABS);
+    assert_close!(pv, expected, PV_REL, JPY_ABS);
     Ok(())
 }
 
@@ -1526,7 +1509,7 @@ fn published_observation_shift_pays_three_yen() -> Result<()> {
     assert_close!(
         evaluate(&prepared, &market)?.get(Currency::JPY),
         3.0,
-        REL,
+        PV_REL,
         JPY_ABS
     );
     Ok(())
@@ -1616,7 +1599,7 @@ fn daily_rounding_inclusive_spread_and_power_use_daily_factors() -> Result<()> {
         assert_close!(
             evaluate(&prepared, &market)?.get(Currency::JPY),
             expected,
-            REL,
+            PV_REL,
             JPY_ABS
         );
 
@@ -1714,7 +1697,7 @@ fn different_forward_and_weight_bases_require_daily_factors() -> Result<()> {
         assert_close!(
             evaluate(&prepared, &market)?.get(Currency::JPY),
             expected,
-            REL,
+            PV_REL,
             JPY_ABS
         );
         // 基準が違えば NL/365 の場合も日次積が要るので、禁止すると準備で拒否する。
@@ -1782,7 +1765,7 @@ fn rounded_forward_compounds_before_aggregated_rounding() -> Result<()> {
         assert_close!(
             evaluate(&prepared, &market)?.get(Currency::JPY),
             expected,
-            REL,
+            PV_REL,
             JPY_ABS
         );
         let ctx = Context {
@@ -1851,6 +1834,131 @@ fn fixing_in_unresolved_choice_is_required_only_when_selected() -> Result<()> {
     record.option = 1;
     let prepared = prepare(&[record])?;
     let pv = evaluate(&prepared, &Market::from_theta(structure, vec![0.0])?)?.get(Currency::JPY);
-    assert_close!(pv, 100.0, REL, JPY_ABS);
+    assert_close!(pv, 100.0, PV_REL, JPY_ABS);
+    Ok(())
+}
+
+#[test]
+fn irs_and_basis_swap_match_coupon_cashflows() -> Result<()> {
+    // ターム対固定とターム対タームを同じ脚の構築 API で契約にし、スプレッドの符号を照合する。
+    let as_of = futures::AS_OF.parse()?;
+    let market = fx_market(as_of)?;
+    let reference = fx_reference("USDJPY".parse()?)?;
+    let config = fx_config("")?;
+    let context = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &FixingStore::default(),
+    };
+    let notional = Money::new(Currency::USD, 100_000.0);
+    let schedule = explicit_schedule()?;
+    let leg = |side, coupon| {
+        Leg::new(
+            side,
+            notional,
+            schedule.clone(),
+            DayCount::Act365F,
+            coupon,
+            NotionalExchange::None,
+            PaymentRounding::None,
+        )
+    };
+    let tau = 7.0 / 365.0;
+    let rate = (futures::RATE * tau).exp_m1() / tau;
+    let dfs = (-futures::RATE * model_time(as_of, futures::at("2026-01-14", "16:00")?).years())
+        .exp()
+        + (-futures::RATE * model_time(as_of, futures::at("2026-01-21", "16:00")?).years()).exp();
+    for (pay, expected_coupon) in [
+        (
+            Coupon::Fixed {
+                rate: Rate::new(0.03, DayCount::Act365F, Simple)?,
+            },
+            rate + 0.001 - 0.03,
+        ),
+        (term("TERM", 0.002)?, -0.001),
+    ] {
+        let swap = StandardProduct::Swap(Swap::from_legs([
+            leg(PayReceive::Pay, pay)?,
+            leg(PayReceive::Receive, term("TERM", 0.001)?)?,
+        ])?);
+        assert_close!(
+            price(&swap, &market, &context)?.get(Currency::USD),
+            notional.amount() * tau * dfs * expected_coupon,
+            PV_REL,
+            USD_ABS
+        );
+    }
+    Ok(())
+}
+
+#[test]
+fn term_coupon_payment_delay_is_recorded() -> Result<()> {
+    // ターム金利のクーポンを期末から 2 日遅れて払うと、自然な支払日（期末）とのずれのタイミング調整を
+    // 省く近似（PaymentDelayIgnored）を記録し、予想は期末の DF 比、割引は実際の支払日で求める。
+    // 近似を禁止すると準備できない。
+    let as_of = futures::AS_OF.parse()?;
+    let market = futures::market(as_of, 0.0)?;
+    let reference = futures::refdata("0D".parse()?, "16:00:00")?;
+    let fixings = FixingStore::default();
+    // 01-07〜01-14 の 1 期間を、期末から 2 日遅れた 01-16 に払う。
+    let schedule = schedule(
+        &["2026-01-07".parse()?, "2026-01-14".parse()?],
+        &["2026-01-16".parse()?],
+    )?;
+    let swap = StandardProduct::Swap(Swap::from_legs([
+        Leg::new(
+            PayReceive::Receive,
+            Money::new(Currency::USD, 100_000.0),
+            schedule.clone(),
+            DayCount::Act365F,
+            term("TERM", 0.0)?,
+            NotionalExchange::None,
+            PaymentRounding::None,
+        )?,
+        Leg::new(
+            PayReceive::Pay,
+            Money::new(Currency::USD, 100_000.0),
+            schedule,
+            DayCount::Act365F,
+            Coupon::Fixed {
+                rate: Rate::new(0.0, DayCount::Act365F, Simple)?,
+            },
+            NotionalExchange::None,
+            PaymentRounding::None,
+        )?,
+    ])?);
+    let prepare_with = |config: &PricingConfig| {
+        prepare_standard(
+            &swap,
+            futures::START.parse()?,
+            market.structure(),
+            &Context {
+                as_of,
+                refdata: &reference,
+                config,
+                fixings: &fixings,
+            },
+        )
+    };
+    let plan = prepare_with(&futures::config("")?)?;
+    // 警告の一覧全体を比べ、支払いの遅れの近似だけを記録することを確かめる。
+    assert_eq!(plan.warnings(), [Warning::PaymentDelayIgnored]);
+    let tau = 7.0 / 365.0;
+    let forward = (futures::RATE * tau).exp_m1() / tau;
+    let df =
+        (-futures::RATE * model_time(as_of, futures::at("2026-01-16", "16:00")?).years()).exp();
+    assert_close!(
+        evaluate(&plan, &market)?.get(Currency::USD),
+        100_000.0 * tau * forward * df,
+        PV_REL,
+        USD_ABS
+    );
+    assert!(matches!(
+        prepare_with(&futures::config(
+            "[conventions]\ndisabled = ['payment_delay_ignored']"
+        )?),
+        Err(Error::Unsupported(_))
+    ));
     Ok(())
 }
diff --git a/crates/rhodium-pricing/tests/m1_linear.rs b/crates/rhodium-pricing/tests/m1_linear.rs
deleted file mode 100644
index b9b675bc1888ec757e6907f4270152d6124a3d90..0000000000000000000000000000000000000000
--- a/crates/rhodium-pricing/tests/m1_linear.rs
+++ /dev/null
@@ -1,2168 +0,0 @@
-//! M1 の標準取引を独立したキャッシュフローと閉形式に照合する。
-
-use std::{num::NonZeroU32, sync::Arc};
-
-use rhodium_core::{
-    Error, Result, Warning,
-    calendar::{
-        BaseCalendar, BusinessDayConvention, Extrapolation as CalendarExtrapolation, WeekendRule,
-    },
-    currency::{Asset, Currency},
-    date::{Date, DateInterval, Weekday},
-    day_count::DayCount,
-    fixing::{FixingStore, FixingStoreParts, ObservationPoint},
-    fx::{CcyPair, FxRate},
-    id::{CalendarId, IndexId},
-    money::Money,
-    quantity::{Rate, RateFuturesPrice, Simple, Spread},
-    refdata::{
-        Convention, FraConvention, FxConvention, FxForwardConvention, FxSpotLag, IndexKind,
-        PublicationRule, RateFutureRate, RefData, RefDataParts,
-    },
-    schedule::{
-        ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination, ScheduleSpec,
-    },
-    time::{Instant, model_time},
-};
-use rhodium_market::{
-    curve::{CurveRole, CurveSpec, Extrapolation},
-    fx::{FxAnchor, FxCollateral, FxSpec},
-    market::Market,
-    quote::{QuoteKind, QuoteSet},
-    structure::{BlockId, MarketStructure},
-};
-use rhodium_pricing::{
-    calibration::{GroupOutcome, calibrate},
-    config::{AnalyticPricer, Method, PricingConfig},
-    context::Context,
-    definition::{
-        CalibrationTarget, DirectFx, ElementSelector, FxDefinition, FxQuotes, InstrumentBinding,
-        MarketDefinition, QuoteTransform,
-    },
-    prepare::{PrepareInput, Prepared, evaluate, prepare},
-    resolve::CollateralKey,
-};
-use rhodium_products::{
-    ir::{
-        Contract, Product,
-        analysis::{Linearity, Observation, linearity, observations},
-    },
-    standard::{
-        Coupon, Fra, FxForward, FxSwap, InstrumentTemplate, Leg, Notional, NotionalExchange,
-        PayReceive, PaymentRounding, Position, RateFuture, StandardProduct, Swap,
-    },
-};
-use rhodium_testkit::assert_close;
-use scioxide_ad::Tape;
-use scioxide_core::{ExecutionContext, InvalidInputReason, SciOxideError, argument};
-
-#[path = "fixtures/futures.rs"]
-mod futures;
-
-use futures::{
-    AS_OF, CLEARING, FINAL, MULTIPLIER, PILLAR, RATE, START, at, config, curve, futures_definition,
-    futures_template, market, refdata, reference_parts,
-};
-
-fn future(position: Position, price: f64) -> Result<StandardProduct> {
-    Ok(StandardProduct::RateFuture(RateFuture::new(
-        CLEARING.parse()?,
-        NonZeroU32::MIN,
-        RateFuturesPrice::new(price, DayCount::Act365F)?,
-        position,
-        None,
-    )))
-}
-
-fn futures_pv(
-    product: &StandardProduct,
-    context: &Context<'_>,
-    market: &Market<f64>,
-) -> Result<(f64, Vec<Warning>)> {
-    let trade_date = START.parse()?;
-    let analytic = prepare(
-        PrepareInput::Standard {
-            product,
-            trade_date,
-        },
-        &CollateralKey::Uncollateralized,
-        &[],
-        market.structure(),
-        context,
-    )?;
-    let contract = product.to_product(trade_date, context.refdata)?;
-    let ir = prepare(
-        PrepareInput::Ir(&contract),
-        &CollateralKey::Uncollateralized,
-        &[],
-        market.structure(),
-        context,
-    )?;
-    assert_eq!(
-        analytic.method(),
-        Method::Analytic(AnalyticPricer::RateFuture)
-    );
-    assert_eq!(ir.method(), Method::Linear);
-    let pv = evaluate(&analytic, market)?.get(Currency::USD);
-    assert_close!(
-        pv,
-        evaluate(&ir, market)?.get(Currency::USD),
-        1e-12_f64,
-        1e-8_f64
-    );
-    assert_eq!(analytic.warnings(), ir.warnings());
-    Ok((pv, analytic.warnings().to_vec()))
-}
-
-fn oracle_price(as_of: Instant, sessions: &[Instant], vol: f64) -> Result<f64> {
-    let ts: Vec<_> = std::iter::once(0.0)
-        .chain(sessions.iter().map(|&date| model_time(as_of, date).years()))
-        .collect();
-    let t1 = model_time(as_of, at(FINAL, "16:00")?).years();
-    let t2 = model_time(as_of, at("2026-01-14", "16:00")?).years();
-    let terminal = ts[ts.len() - 1];
-    // a=0 の H は t。共分散 σ² min(s,t) に各清算時点の係数を直接掛ける。
-    let covariance = ts
-        .windows(2)
-        .map(|w| (w[1] - w[0]) * w[0].min(t1))
-        .sum::<f64>()
-        - terminal * terminal.min(t1);
-    let tau = 7.0 / 365.0;
-    let log_mean = RATE * (t2 - t1) + vol * vol * (t2 - t1) * (t2 * t1 + covariance);
-    Ok(1.0 - log_mean.exp_m1() / tau)
-}
-
-#[test]
-fn futures_margin_ir_vs_analytic() -> Result<()> {
-    // 日次に共通の価格を 実装と別の式で手計算し、将来の二回目以後の値洗いが消えることを確かめる。
-    let as_of = AS_OF.parse()?;
-    let reference = refdata("0D".parse()?, "16:00:00")?;
-    let config = config("")?;
-    let fixings = FixingStore::default();
-    let context = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &fixings,
-    };
-    let margins = [
-        at(START, "16:00")?,
-        at("2026-01-06", "16:00")?,
-        at(FINAL, "16:00")?,
-    ];
-    let df = (-RATE * model_time(as_of, margins[0]).years()).exp();
-    for vol in [0.0, 0.5] {
-        let market = market(as_of, vol)?;
-        for (side, sign) in [(Position::Long, 1.0), (Position::Short, -1.0)] {
-            let (actual, warnings) = futures_pv(&future(side, 0.95)?, &context, &market)?;
-            let expected = sign * MULTIPLIER * df * (oracle_price(as_of, &margins, vol)? - 0.95);
-            assert_close!(actual, expected, 1e-12_f64, 1e-8_f64);
-            assert_eq!(warnings, []);
-        }
-    }
-    Ok(())
-}
-
-#[test]
-fn futures_settlement_delay() -> Result<()> {
-    // 同日内の遅れも翌日の受渡も許可制の近似になり、禁止した設定では評価を拒否する。
-    let as_of = AS_OF.parse()?;
-    let market = market(as_of, 0.0)?;
-    let config_ok = config("")?;
-    let config_disabled = config("[conventions]\ndisabled = ['payment_delay_ignored']")?;
-    for (lag, time) in [("0D", "18:00:00"), ("1BD", "08:00:00")] {
-        let reference = refdata(lag.parse()?, time)?;
-        let fixings = FixingStore::default();
-        let context = Context {
-            as_of,
-            refdata: &reference,
-            config: &config_ok,
-            fixings: &fixings,
-        };
-        let product = future(Position::Long, 0.95)?;
-        assert_eq!(
-            futures_pv(&product, &context, &market)?.1,
-            vec![Warning::PaymentDelayIgnored]
-        );
-        let disabled = Context {
-            config: &config_disabled,
-            ..context
-        };
-        assert!(matches!(
-            prepare(
-                PrepareInput::Standard {
-                    product: &product,
-                    trade_date: START.parse()?
-                },
-                &CollateralKey::Uncollateralized,
-                &[],
-                market.structure(),
-                &disabled,
-            ),
-            Err(Error::Unsupported(_))
-        ));
-    }
-    Ok(())
-}
-
-#[test]
-fn futures_vs_forward_settlement() -> Result<()> {
-    // 実装の先物価格が日次の値洗い口座 Q_M の式に一致し、最終日に一度だけ払う契約の
-    // 価格（最終受渡の割引債を数値基準とする式）とは凸性の分だけ異なることを確かめる。
-    // oracle 同士の比較だけでは、実装が満期一括の数値基準へ退行しても検出できないため、
-    // ライブラリの評価から先物価格を逆算して照合する。
-    let as_of = AS_OF.parse()?;
-    let reference = refdata("0D".parse()?, "16:00:00")?;
-    let config = config("")?;
-    let fixings = FixingStore::default();
-    let context = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &fixings,
-    };
-    let vol = 0.5;
-    let market = market(as_of, vol)?;
-    let first = at(START, "16:00")?;
-    let final_at = at(FINAL, "16:00")?;
-    let daily = [first, at("2026-01-06", "16:00")?, final_at];
-    // PV は次の値洗い M·P(s_1)·(F − K) だけなので、そこから実装の F を求める。
-    let trade_price = 0.95;
-    let df = (-RATE * model_time(as_of, first).years()).exp();
-    let (pv, _) = futures_pv(&future(Position::Long, trade_price)?, &context, &market)?;
-    let implied = pv / (MULTIPLIER * df) + trade_price;
-    assert_close!(
-        implied,
-        oracle_price(as_of, &daily, vol)?,
-        1e-12_f64,
-        1e-8_f64
-    );
-    let forward = oracle_price(as_of, &[final_at], vol)?;
-    assert!((implied - forward).abs() > 1e-6);
-    Ok(())
-}
-
-#[test]
-fn futures_publication_and_payment_boundaries_use_realized_margin() -> Result<()> {
-    // 確定済みで未払いの受渡は取引価格との差額。決済時点を除くと次回の差額だけ残る。
-    let as_of = at(START, "16:00")?;
-    let reference = refdata("0D".parse()?, "16:00:00")?;
-    let config = config("")?;
-    let product = future(Position::Long, 0.95)?;
-    let market = market(as_of, 0.0)?;
-    let missing = FixingStore::default();
-    let ctx = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &missing,
-    };
-    assert!(matches!(
-        prepare(
-            PrepareInput::Standard {
-                product: &product,
-                trade_date: START.parse()?
-            },
-            &CollateralKey::Uncollateralized,
-            &[],
-            market.structure(),
-            &ctx,
-        ),
-        Err(Error::MissingData(_))
-    ));
-    let fixings = FixingStore::new(
-        &reference,
-        FixingStoreParts {
-            published: vec![(
-                CLEARING.parse()?,
-                ObservationPoint::Date(START.parse()?),
-                0.954,
-            )],
-            ..FixingStoreParts::default()
-        },
-    )?;
-    let context = Context {
-        fixings: &fixings,
-        ..ctx
-    };
-    let (excluded, _) = futures_pv(&product, &context, &market)?;
-
-    let include: PricingConfig = toml::from_str(
-        r#"
-version = 1
-reporting_currency = "USD"
-settlement = "include"
-[[discounting]]
-currency = "USD"
-collateral = "uncollateralized"
-curve = "USD"
-[forecasting]
-TERM = "USD"
-"#,
-    )
-    .map_err(|error| Error::InvalidData(error.to_string()))?;
-    let included = futures_pv(
-        &product,
-        &Context {
-            config: &include,
-            ..context
-        },
-        &market,
-    )?
-    .0;
-    assert_close!(included - excluded, 4.0_f64, 1e-12_f64, 1e-8_f64);
-    Ok(())
-}
-
-#[test]
-fn fra_settlement_matches_start_discounted_cashflow() -> Result<()> {
-    // FRA は期末のクーポンではなく、同じターム金利で割り戻した差額を期首に受け渡す。
-    let as_of = AS_OF.parse()?;
-    let reference = refdata("0D".parse()?, "16:00:00")?;
-    let config = config("")?;
-    let fixings = FixingStore::default();
-    let context = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &fixings,
-    };
-    let market = market(as_of, 0.0)?;
-    let product = StandardProduct::Fra(Fra::new(
-        Money::new(Currency::USD, 100_000.0),
-        "TERM".parse()?,
-        FINAL.parse()?,
-        0.03,
-        Position::Long,
-    )?);
-    let plan = prepare(
-        PrepareInput::Standard {
-            product: &product,
-            trade_date: START.parse()?,
-        },
-        &CollateralKey::Uncollateralized,
-        &[],
-        market.structure(),
-        &context,
-    )?;
-    let tau = 7.0 / 365.0;
-    let forward = (RATE * tau).exp_m1() / tau;
-    let df = (-RATE * model_time(as_of, at(FINAL, "16:00")?).years()).exp();
-    let expected = df * 100_000.0 * tau * (forward - 0.03) / (1.0 + tau * forward);
-    assert_close!(
-        evaluate(&plan, &market)?.get(Currency::USD),
-        expected,
-        1e-12_f64,
-        1e-8_f64
-    );
-    Ok(())
-}
-
-/// 期首 `start` に受け渡す FRA の決済額 N·τ·(L − K)/(1 + τ'·L) と、同じ支払いに並べた
-/// 通常のタームクーポン N·τ·L₁ の独自契約。L は `fixing` の、L₁ は `coupon` の観測。
-fn fra_with_coupon(
-    reference: &RefData,
-    fixing: &str,
-    coupon: Option<&str>,
-    denominator: f64,
-    start: Instant,
-) -> Result<Product> {
-    let tau = 7.0 / 365.0;
-    let index: IndexId = "TERM".parse()?;
-    Product::new(
-        Contract::build(|builder| {
-            let floating = builder.fixing(&index, ObservationPoint::Date(fixing.parse()?));
-            let strike = builder.constant(0.03);
-            let difference = builder.sub(floating, strike);
-            let accrual = builder.constant(100_000.0 * tau);
-            let numerator = builder.mul([accrual, difference]);
-            let tau = builder.constant(denominator);
-            let interest = builder.mul([tau, floating]);
-            let one = builder.constant(1.0);
-            let discount = builder.add([one, interest]);
-            let mut amount = builder.div(numerator, discount);
-            if let Some(date) = coupon {
-                let rate = builder.fixing(&index, ObservationPoint::Date(date.parse()?));
-                let accrual = builder.constant(100_000.0 * 7.0 / 365.0);
-                let coupon = builder.mul([accrual, rate]);
-                amount = builder.add([amount, coupon]);
-            }
-            Ok(builder.pay(Asset::Cash(Currency::USD), amount, start))
-        })?,
-        reference,
-    )
-}
-
-#[test]
-fn fra_identity_applies_to_its_own_settlement_and_accrual_only() -> Result<()> {
-    // FRA の期首払いの扱い（支払いの遅れを期首から測る）は、その決済額の中の観測にだけ当てる。
-    // 01-14 の支払いに、01-14 に始まる FRA の決済と、01-07〜01-14 の通常のクーポンを並べる。
-    // クーポンは期末の 01-14 に払うので遅れは無く、支払い全体に FRA の扱いを広げると
-    // 期首の 01-07 から 5 営業日の遅れと誤って測り、近似を記録する（上限が小さければ拒否する）。
-    // また分母の τ' が指標の日数係数 τ と違う独自契約は 1/(1 + τ'L) が割引因子にならないので、
-    // 恒等式を当てずリニアでないと判定する。
-    let as_of = AS_OF.parse()?;
-    let reference = refdata("0D".parse()?, "16:00:00")?;
-    let config = config("")?;
-    let fixings = FixingStore::default();
-    let context = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &fixings,
-    };
-    let market = market(as_of, 0.0)?;
-    let start = at("2026-01-14", "16:00")?;
-    let tau = 7.0 / 365.0;
-    let mixed = fra_with_coupon(&reference, "2026-01-14", Some("2026-01-07"), tau, start)?;
-    let plan = prepare(
-        PrepareInput::Ir(&mixed),
-        &CollateralKey::Uncollateralized,
-        &[],
-        market.structure(),
-        &context,
-    )?;
-    assert_eq!(plan.warnings(), &[] as &[Warning]);
-    let forward = (RATE * tau).exp_m1() / tau;
-    let df = (-RATE * model_time(as_of, start).years()).exp();
-    let expected = df * 100_000.0 * tau * ((forward - 0.03) / (1.0 + tau * forward) + forward);
-    assert_close!(
-        evaluate(&plan, &market)?.get(Currency::USD),
-        expected,
-        1e-12_f64,
-        1e-8_f64
-    );
-
-    // 分母の日数係数を Act/360 で数えた独自契約は FRA の恒等式に当たらない。
-    let mismatched = fra_with_coupon(&reference, "2026-01-14", None, 7.0 / 360.0, start)?;
-    assert!(matches!(
-        linearity(&mismatched),
-        Linearity::NonLinear(reason) if reason.contains("division")
-    ));
-    // 同じ τ なら、小数で書いた τ の丸めの範囲で恒等式に当たる。
-    let decimal = fra_with_coupon(
-        &reference,
-        "2026-01-14",
-        None,
-        0.019_178_082_191_780_82,
-        start,
-    )?;
-    assert_eq!(linearity(&decimal), Linearity::Linear);
-    Ok(())
-}
-
-#[test]
-fn futures_quote_calibrates_the_trade_to_par() -> Result<()> {
-    // 先物のマークと独立した参照価格で校正し、同じ限月を約定した取引の PV がゼロであることを確かめる。
-    let as_of = AS_OF.parse()?;
-    let reference = refdata("0D".parse()?, "16:00:00")?;
-    let config = config("")?;
-    let context = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &FixingStore::default(),
-    };
-    let template = futures_template()?;
-    let definition = futures_definition()?;
-    let quote = oracle_price(
-        as_of,
-        &[
-            at(START, "16:00")?,
-            at("2026-01-06", "16:00")?,
-            at(FINAL, "16:00")?,
-        ],
-        0.5,
-    )?;
-    let market = calibrate(
-        &definition,
-        &QuoteSet::new([
-            ("USD.FUT.2026-01".parse()?, quote),
-            ("USD.LGM.MR".parse()?, 0.0),
-            ("USD.LGM.VOL".parse()?, 0.5),
-        ])?,
-        &context,
-    )?;
-    let product = template.build_for_valuation(quote, START.parse()?, as_of, &reference)?;
-    assert_close!(
-        futures_pv(&product, &context, market.market())?.0,
-        0.0_f64,
-        1e-12_f64,
-        1e-8_f64
-    );
-    Ok(())
-}
-
-fn fx_reference(pair: CcyPair) -> Result<RefData> {
-    RefData::new(fx_parts(pair)?)
-}
-
-fn fx_parts(pair: CcyPair) -> Result<RefDataParts> {
-    let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
-    let clock = parts.currencies[&Currency::USD].settlement_time.clone();
-    parts.conventions.insert(
-        "USDJPY-FWD".parse()?,
-        Convention::FxForward(Box::new(FxForwardConvention {
-            pair: "USDJPY".parse()?,
-            bdc: BusinessDayConvention::Following,
-        })),
-    );
-    // スポットラグ、受渡の暦（両通貨の休日の和）と受渡の時刻はペアの慣行から読む。
-    let usd_jpy: CcyPair = "USDJPY".parse()?;
-    parts.fx_conventions.insert(
-        usd_jpy,
-        FxConvention {
-            pair: usd_jpy,
-            spot_lag: FxSpotLag::TwoBusinessDays,
-            settlement: clock,
-            pip_exponent: 2,
-        },
-    );
-    let mut currency = parts.currencies[&Currency::USD].clone();
-    currency.code = Currency::JPY;
-    currency.minor_units = 0;
-    parts.currencies.insert(currency.code, currency);
-    let mut term = parts.indices[&"TERM".parse()?].clone();
-    term.id = "JPY-TERM".parse()?;
-    if let IndexKind::TermRate { ccy, .. } = &mut term.kind {
-        *ccy = Currency::JPY;
-    }
-    parts.indices.insert(term.id.clone(), term.clone());
-    term.id = "FX-FIX".parse()?;
-    term.kind = IndexKind::FxFixing { pair };
-    parts.indices.insert(term.id.clone(), term);
-    Ok(parts)
-}
-
-/// USDJPY の慣行のスポットラグを `spot_lag` にし、JPY の決済カレンダーに 2026-01-13（火）の
-/// 祝日を置いた参照データ。ペアの暦は USD と JPY の休日の和になる。
-fn fx_reference_with_lag(spot_lag: FxSpotLag) -> Result<RefData> {
-    let mut parts = fx_parts("USDJPY".parse()?)?;
-    let tokyo: CalendarId = "JPN".parse()?;
-    parts.calendars.insert(
-        tokyo.clone(),
-        BaseCalendar::new(
-            "2025-01-01".parse()?..="2026-12-31".parse()?,
-            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
-            vec!["2026-01-13".parse()?],
-            CalendarExtrapolation::Strict,
-        )?,
-    );
-    if let Some(jpy) = parts.currencies.get_mut(&Currency::JPY) {
-        jpy.settlement_calendar = tokyo;
-    }
-    if let Some(convention) = parts.fx_conventions.get_mut(&"USDJPY".parse()?) {
-        convention.spot_lag = spot_lag;
-    }
-    RefData::new(parts)
-}
-
-#[test]
-#[expect(
-    clippy::wildcard_enum_match_arm,
-    reason = "#196 で各 variant を網羅する"
-)]
-fn mtm_xccy_reset_fixing_dates_follow_the_pair_spot_lag() -> Result<()> {
-    // 元本のリセットは期首がスポット受渡になる日に観測する。期首は 2026-01-07（水）と 01-14（水）で、
-    // ペアの暦は USD と JPY の休日の和（JPY は 01-13 が祝日）。
-    // T+1: 01-07 の 1 営業日前は 01-06、01-14 の 1 営業日前は 01-13 を飛ばして 01-12。
-    // T+2: 01-07 の 2 営業日前は 01-05、01-14 の 2 営業日前は 01-12、01-09。
-    let swap = mtm_swap()?;
-    for (spot_lag, expected) in [
-        (FxSpotLag::OneBusinessDay, ["2026-01-06", "2026-01-12"]),
-        (FxSpotLag::TwoBusinessDays, ["2026-01-05", "2026-01-09"]),
-    ] {
-        let reference = fx_reference_with_lag(spot_lag)?;
-        let product = swap.to_product(START.parse()?, &reference)?;
-        let fx: IndexId = "FX-FIX".parse()?;
-        let mut dates: Vec<_> = observations(&product)
-            .into_iter()
-            .filter_map(|observation| match observation {
-                Observation::Fixing {
-                    index,
-                    point: ObservationPoint::Date(date),
-                } if *index == fx => Some(date),
-                _ => None,
-            })
-            .collect();
-        dates.sort_unstable();
-        let expected = expected
-            .iter()
-            .map(|date| date.parse())
-            .collect::<Result<Vec<Date>>>()?;
-        assert_eq!(dates, expected, "{spot_lag:?}");
-    }
-    Ok(())
-}
-
-#[test]
-fn fx_fixing_spot_value_uses_the_pair_spot_lag() -> Result<()> {
-    // 受渡日の指定の無い FX fixing は、fixing 日からペアのスポットラグだけ進めた受渡日のレートを観測する。
-    // 01-12（月）の fixing の受渡日は、T+1 では 01-13 の祝日を飛ばして 01-14、T+2 では 01-15。
-    // 支払いは 01-20 で、その円の受渡を P_JPY で割り引く。F(v) = 150·exp(−(r_USD − r_JPY)·t_v)。
-    let as_of = AS_OF.parse()?;
-    let market = fx_market(as_of)?;
-    let config = fx_config("")?;
-    let fixings = FixingStore::default();
-    let pay = at("2026-01-20", "16:00")?;
-    let fx: IndexId = "FX-FIX".parse()?;
-    let discount = |rate: f64, at| (-rate * model_time(as_of, at).years()).exp();
-    let forward = |at| 150.0 * discount(RATE - 0.01, at);
-    let mut values = Vec::new();
-    for (spot_lag, value_date) in [
-        (FxSpotLag::OneBusinessDay, "2026-01-14"),
-        (FxSpotLag::TwoBusinessDays, "2026-01-15"),
-    ] {
-        let reference = fx_reference_with_lag(spot_lag)?;
-        let product = Product::new(
-            Contract::build(|builder| {
-                let rate = builder.fixing(&fx, ObservationPoint::Date("2026-01-12".parse()?));
-                Ok(builder.pay(Asset::Cash(Currency::JPY), rate, pay))
-            })?,
-            &reference,
-        )?;
-        let ctx = Context {
-            as_of,
-            refdata: &reference,
-            config: &config,
-            fixings: &fixings,
-        };
-        let plan = prepare(
-            PrepareInput::Ir(&product),
-            &CollateralKey::Uncollateralized,
-            &[],
-            market.structure(),
-            &ctx,
-        )?;
-        // 受渡日と支払日が違うので、タイミング調整を省いた近似を記録する。
-        assert!(plan.warnings().contains(&Warning::PaymentDelayIgnored));
-        let pv = evaluate(&plan, &market)?.get(Currency::JPY);
-        assert_close!(
-            pv,
-            forward(at(value_date, "16:00")?) * discount(0.01, pay),
-            1e-12_f64,
-            1e-6_f64
-        );
-        values.push(pv);
-    }
-    assert!((values[0] - values[1]).abs() > 1e-6);
-    Ok(())
-}
-
-fn fx_market(as_of: Instant) -> Result<Market<f64>> {
-    let jpy = CurveSpec::new(
-        "JPY".parse()?,
-        Currency::JPY,
-        [CurveRole::Discount, CurveRole::Projection],
-        [PILLAR.parse()?],
-        Extrapolation::Reject,
-    )?;
-    let pair: CcyPair = "USDJPY".parse()?;
-    let structure = Arc::new(MarketStructure::with_fx(
-        as_of,
-        [curve()?, jpy],
-        [FxSpec {
-            pair,
-            anchor: FxAnchor::Direct {
-                spot: as_of,
-                short_starts: vec![],
-                market_collateral: "csa".parse()?,
-            },
-            collaterals: vec![FxCollateral {
-                id: "csa".parse()?,
-                base: "USD".parse()?,
-                quote: "JPY".parse()?,
-            }],
-            via: None,
-        }],
-    )?);
-    let t = model_time(as_of, PILLAR.parse()?).years();
-    let mut theta = vec![0.0; 3];
-    for (id, value) in [
-        (BlockId::nodes("USD".parse()?), -RATE * t),
-        (BlockId::nodes("JPY".parse()?), -0.01 * t),
-        (FxSpec::marks_block(pair)?, 150.0),
-    ] {
-        let range = structure
-            .layout()
-            .range(&id)
-            .ok_or_else(|| Error::Internal("missing FX test block".into()))?;
-        theta[range.start] = value;
-    }
-    Market::from_theta(structure, theta)
-}
-
-fn fx_config(conventions: &str) -> Result<PricingConfig> {
-    toml::from_str(&format!(
-        r#"
-version = 1
-reporting_currency = "USD"
-[[discounting]]
-currency = "USD"
-collateral = "uncollateralized"
-curve = "USD"
-[[discounting]]
-currency = "JPY"
-collateral = "uncollateralized"
-curve = "JPY"
-[forecasting]
-TERM = "USD"
-JPY-TERM = "JPY"
-{conventions}
-"#
-    ))
-    .map_err(|error| Error::InvalidData(error.to_string()))
-}
-
-fn price(
-    product: &StandardProduct,
-    market: &Market<f64>,
-    context: &Context<'_>,
-) -> Result<rhodium_core::money::MoneyBag> {
-    let plan = prepare(
-        PrepareInput::Standard {
-            product,
-            trade_date: START.parse()?,
-        },
-        &CollateralKey::Uncollateralized,
-        &[],
-        market.structure(),
-        context,
-    )?;
-    evaluate(&plan, market)
-}
-
-#[test]
-fn fx_forward_and_swap_exchange_both_currencies() -> Result<()> {
-    // 両通貨の各受渡を割り引き、近い受渡と遠い受渡の向きが逆になることを確かめる。
-    let as_of = AS_OF.parse()?;
-    let market = fx_market(as_of)?;
-    let pair: CcyPair = "USDJPY".parse()?;
-    let reference = fx_reference(pair)?;
-    let config = fx_config("")?;
-    let context = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &FixingStore::default(),
-    };
-    let near = at("2026-01-07", "16:00")?;
-    let far = at("2026-01-21", "16:00")?;
-    let notional = Money::new(Currency::USD, 1000.0);
-    let first = FxRate::forward(pair, 151.0, near)?;
-    let second = FxRate::forward(pair, 152.0, far)?;
-    let forward = StandardProduct::FxForward(FxForward::new(notional, first, Position::Long)?);
-    let swap = StandardProduct::FxSwap(FxSwap::new(notional, first, second, Position::Long)?);
-    let df = |rate: f64, t| (-rate * model_time(as_of, t).years()).exp();
-
-    let pv = price(&forward, &market, &context)?;
-    assert_close!(
-        pv.get(Currency::USD),
-        1000.0 * df(RATE, near),
-        1e-12_f64,
-        1e-8_f64
-    );
-    assert_close!(
-        pv.get(Currency::JPY),
-        -151_000.0 * df(0.01, near),
-        1e-12_f64,
-        1e-6_f64
-    );
-    let pv = price(&swap, &market, &context)?;
-    assert_close!(
-        pv.get(Currency::USD),
-        1000.0 * (df(RATE, near) - df(RATE, far)),
-        1e-12_f64,
-        1e-8_f64
-    );
-    assert_close!(
-        pv.get(Currency::JPY),
-        -151_000.0 * df(0.01, near) + 152_000.0 * df(0.01, far),
-        1e-12_f64,
-        1e-6_f64
-    );
-
-    // 市場の二つのアウトライトで約定すると、各受渡の報告通貨建て PV がゼロになる。
-    let outright = |at| FxRate::forward(pair, 150.0 * df(RATE, at) / df(0.01, at), at);
-    let par = StandardProduct::FxSwap(FxSwap::new(
-        notional,
-        outright(near)?,
-        outright(far)?,
-        Position::Long,
-    )?);
-    let pv = price(&par, &market, &context)?;
-    assert_close!(
-        pv.get(Currency::JPY) + 150.0 * pv.get(Currency::USD),
-        0.0_f64,
-        1e-12_f64,
-        1e-6_f64
-    );
-    Ok(())
-}
-
-fn explicit_schedule() -> Result<ScheduleSpec> {
-    let start = "2026-01-07".parse::<Date>()?;
-    let middle = "2026-01-14".parse::<Date>()?;
-    let end = "2026-01-21".parse::<Date>()?;
-    let periods = [(start, middle), (middle, end)]
-        .into_iter()
-        .map(|(start, end)| {
-            Ok(ExplicitPeriod {
-                accrual: DateInterval::new(start, end)?,
-                payment: end,
-                reset: None,
-                settlement: None,
-                references: vec![],
-            })
-        })
-        .collect::<Result<_>>()?;
-    Ok(ScheduleSpec::Explicit(ExplicitPeriods::new(
-        periods,
-        ExplicitDatePolicy::KeepAsConfirmed,
-        PeriodCombination::Sum,
-        None,
-    )?))
-}
-
-fn term(index: &str, spread: f64) -> Result<Coupon> {
-    Ok(Coupon::Term {
-        index: index.parse()?,
-        spread: Spread::new(spread)?,
-        gearing: 1.0,
-    })
-}
-
-#[test]
-fn irs_and_basis_swap_match_coupon_cashflows() -> Result<()> {
-    // ターム対固定とターム対タームを同じ脚の構築 API で契約にし、スプレッドの符号を照合する。
-    let as_of = AS_OF.parse()?;
-    let market = fx_market(as_of)?;
-    let reference = fx_reference("USDJPY".parse()?)?;
-    let config = fx_config("")?;
-    let context = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &FixingStore::default(),
-    };
-    let notional = Money::new(Currency::USD, 100_000.0);
-    let schedule = explicit_schedule()?;
-    let leg = |side, coupon| {
-        Leg::new(
-            side,
-            notional,
-            schedule.clone(),
-            DayCount::Act365F,
-            coupon,
-            NotionalExchange::None,
-            PaymentRounding::None,
-        )
-    };
-    let tau = 7.0 / 365.0;
-    let rate = (RATE * tau).exp_m1() / tau;
-    let dfs = (-RATE * model_time(as_of, at("2026-01-14", "16:00")?).years()).exp()
-        + (-RATE * model_time(as_of, at("2026-01-21", "16:00")?).years()).exp();
-    for (pay, expected_coupon) in [
-        (
-            Coupon::Fixed {
-                rate: Rate::new(0.03, DayCount::Act365F, Simple)?,
-            },
-            rate + 0.001 - 0.03,
-        ),
-        (term("TERM", 0.002)?, -0.001),
-    ] {
-        let swap = StandardProduct::Swap(Swap::from_legs([
-            leg(PayReceive::Pay, pay)?,
-            leg(PayReceive::Receive, term("TERM", 0.001)?)?,
-        ])?);
-        assert_close!(
-            price(&swap, &market, &context)?.get(Currency::USD),
-            notional.amount() * tau * dfs * expected_coupon,
-            1e-12_f64,
-            1e-8_f64
-        );
-    }
-    Ok(())
-}
-
-fn mtm_swap() -> Result<StandardProduct> {
-    let notional = Money::new(Currency::JPY, 150_000.0);
-    let schedule = explicit_schedule()?;
-    Ok(StandardProduct::Swap(Swap::from_legs([
-        Leg::new(
-            PayReceive::Receive,
-            notional,
-            schedule.clone(),
-            DayCount::Act365F,
-            term("JPY-TERM", 0.001)?,
-            NotionalExchange::InitialAndFinal,
-            PaymentRounding::None,
-        )?,
-        Leg::with_notional(
-            PayReceive::Pay,
-            Notional::FxReset {
-                base: notional,
-                fixing: "FX-FIX".parse()?,
-            },
-            schedule,
-            DayCount::Act365F,
-            term("TERM", 0.002)?,
-            NotionalExchange::InitialAndFinal,
-            PaymentRounding::None,
-        )?,
-    ])?))
-}
-
-#[test]
-fn mtm_xccy_fx_reset_matches_cashflows_in_both_directions() -> Result<()> {
-    // 同じ基準元本を正向きの fixing では割り、逆向きでは掛ける。期中の交換は元本の差額だけ。
-    let as_of = AS_OF.parse()?;
-    let market = fx_market(as_of)?;
-    let config = fx_config("")?;
-    let disabled = fx_config("[conventions]\ndisabled = ['mtm_fx_forward_product']")?;
-    let swap = mtm_swap()?;
-    let fixings = FixingStore::default();
-    let pay = ["2026-01-07", "2026-01-14", "2026-01-21"].map(|date| at(date, "16:00"));
-    let [start, middle, end] = pay;
-    let [start, middle, end] = [start?, middle?, end?];
-    let discount = |rate: f64, at| (-rate * model_time(as_of, at).years()).exp();
-    let loan = |at| 1000.0 * discount(0.01, at) / discount(RATE, at);
-    let [n0, n1] = [loan(start), loan(middle)];
-    let dollar_coupon = (RATE * (7.0 / 365.0)).exp_m1() + 0.002 * (7.0 / 365.0);
-    let yen_coupon = (0.01_f64 * (7.0 / 365.0)).exp_m1() + 0.001 * (7.0 / 365.0);
-    let usd = n0 * discount(RATE, start) + (n1 - n0 - n0 * dollar_coupon) * discount(RATE, middle)
-        - n1 * (1.0 + dollar_coupon) * discount(RATE, end);
-    let jpy = 150_000.0
-        * (-discount(0.01, start)
-            + yen_coupon * discount(0.01, middle)
-            + (1.0 + yen_coupon) * discount(0.01, end));
-
-    for pair in ["USDJPY", "JPYUSD"] {
-        let reference = fx_reference(pair.parse()?)?;
-        let ctx = Context {
-            as_of,
-            refdata: &reference,
-            config: &config,
-            fixings: &fixings,
-        };
-        let plan = prepare(
-            PrepareInput::Standard {
-                product: &swap,
-                trade_date: START.parse()?,
-            },
-            &CollateralKey::Uncollateralized,
-            &[],
-            market.structure(),
-            &ctx,
-        )?;
-        assert!(
-            plan.blocks()
-                .contains(&FxSpec::marks_block("USDJPY".parse()?)?)
-        );
-        assert_eq!(plan.warnings(), &[Warning::MtmFxForwardProduct]);
-        let actual = evaluate(&plan, &market)?;
-        assert_close!(actual.get(Currency::USD), usd, 1e-12_f64, 1e-8_f64);
-        assert_close!(actual.get(Currency::JPY), jpy, 1e-12_f64, 1e-6_f64);
-        assert!(matches!(
-            price(
-                &swap,
-                &market,
-                &Context {
-                    config: &disabled,
-                    ..ctx
-                }
-            ),
-            Err(Error::Unsupported(_))
-        ));
-    }
-    Ok(())
-}
-
-/// 固定クーポンの `MtM` の XCCY。JPY の脚は元本 150,000 円を固定し、USD の脚の元本を `FX-FIX` で
-/// リセットする。`reset` は 2 期目のリセットの明示日（無ければ規則どおり期首から遡る）。
-fn fixed_mtm_swap(dates: &[&str; 3], reset: Option<&str>) -> Result<StandardProduct> {
-    let schedule = |reset: Option<&str>| -> Result<ScheduleSpec> {
-        let periods = dates
-            .windows(2)
-            .enumerate()
-            .map(|(i, pair)| {
-                let (start, end) = (pair[0].parse::<Date>()?, pair[1].parse::<Date>()?);
-                Ok(ExplicitPeriod {
-                    accrual: DateInterval::new(start, end)?,
-                    payment: end,
-                    reset: if i == 1 {
-                        reset.map(str::parse).transpose()?
-                    } else {
-                        None
-                    },
-                    settlement: None,
-                    references: vec![],
-                })
-            })
-            .collect::<Result<_>>()?;
-        Ok(ScheduleSpec::Explicit(ExplicitPeriods::new(
-            periods,
-            ExplicitDatePolicy::KeepAsConfirmed,
-            PeriodCombination::Sum,
-            None,
-        )?))
-    };
-    let fixed = || -> Result<Coupon> {
-        Ok(Coupon::Fixed {
-            rate: Rate::new(0.03, DayCount::Act365F, Simple)?,
-        })
-    };
-    let notional = Money::new(Currency::JPY, 150_000.0);
-    Ok(StandardProduct::Swap(Swap::from_legs([
-        Leg::new(
-            PayReceive::Receive,
-            notional,
-            schedule(None)?,
-            DayCount::Act365F,
-            fixed()?,
-            NotionalExchange::InitialAndFinal,
-            PaymentRounding::None,
-        )?,
-        Leg::with_notional(
-            PayReceive::Pay,
-            Notional::FxReset {
-                base: notional,
-                fixing: "FX-FIX".parse()?,
-            },
-            schedule(reset)?,
-            DayCount::Act365F,
-            fixed()?,
-            NotionalExchange::InitialAndFinal,
-            PaymentRounding::None,
-        )?,
-    ])?))
-}
-
-#[test]
-fn mtm_xccy_reset_reads_the_fixing_at_the_pair_spot_delivery() -> Result<()> {
-    // 同じ指標・観測日の FX fixing の公表値は一つなので、元本リセットも通常の fixing と同じく
-    // ペアのスポットの受渡（ここでは UTC 06 時。東京 15 時に当たる）の F で予想する。
-    // 新しい元本を交換する瞬間（USD の決済時刻 UTC 16 時）とは別に扱う。
-    // 明示の reset 日 01-08 の fixing はスポットの受渡が 01-12 で、元本の交換 01-14 より早いので
-    // 支払いの遅れの近似を記録する。クーポンは固定にして、FX の予想だけを照合する。
-    let as_of = AS_OF.parse()?;
-    let market = fx_market(as_of)?;
-    let config = fx_config("")?;
-    let fixings = FixingStore::default();
-    let mut parts = fx_parts("USDJPY".parse()?)?;
-    if let Some(convention) = parts.fx_conventions.get_mut(&"USDJPY".parse()?) {
-        convention.settlement.time = "06:00:00".parse()?;
-    }
-    let reference = RefData::new(parts)?;
-    let ctx = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &fixings,
-    };
-    let discount = |rate: f64, at| (-rate * model_time(as_of, at).years()).exp();
-    // F(v) = 150·exp(−(r_USD − r_JPY)·t_v)。リセット後の USD の元本は 150,000 円 / F(v)。
-    let forward = |at| 150.0 * discount(RATE - 0.01, at);
-    let fx: IndexId = "FX-FIX".parse()?;
-
-    // 通常の fixing: 01-05 の fixing の受渡は T+2 の 01-07 06 時。
-    let value = at("2026-01-07", "06:00")?;
-    let pay = at("2026-01-07", "16:00")?;
-    let plain = Product::new(
-        Contract::build(|builder| {
-            let rate = builder.fixing(&fx, ObservationPoint::Date("2026-01-05".parse()?));
-            Ok(builder.pay(Asset::Cash(Currency::JPY), rate, pay))
-        })?,
-        &reference,
-    )?;
-    let plan = prepare(
-        PrepareInput::Ir(&plain),
-        &CollateralKey::Uncollateralized,
-        &[],
-        market.structure(),
-        &ctx,
-    )?;
-    assert_close!(
-        evaluate(&plan, &market)?.get(Currency::JPY),
-        forward(value) * discount(0.01, pay),
-        1e-12_f64,
-        1e-6_f64
-    );
-
-    let dates = ["2026-01-07", "2026-01-14", "2026-01-21"];
-    let [start, middle, end] = dates.map(|date| at(date, "16:00"));
-    let [start, middle, end] = [start?, middle?, end?];
-    let coupon = 0.03 * 7.0 / 365.0;
-    // 支払う USD の脚: 期首に元本 n0 を受け、期中に元本の差額 n1 − n0 を受けてクーポンを払い、
-    // 期末に n1 とクーポンを払う。
-    let usd = |n0: f64, n1: f64| {
-        n0 * discount(RATE, start) + (n1 - n0 - n0 * coupon) * discount(RATE, middle)
-            - n1 * (1.0 + coupon) * discount(RATE, end)
-    };
-    for (reset, second_value, warnings) in [
-        // 規則どおりの reset はスポットの受渡が期首の 01-14 と同じ日なので近似は無い。
-        (None, "2026-01-14", vec![]),
-        // 明示の reset 01-08（木）の T+2 の受渡は 01-12（月）。
-        (
-            Some("2026-01-08"),
-            "2026-01-12",
-            vec![Warning::PaymentDelayIgnored],
-        ),
-    ] {
-        let product = fixed_mtm_swap(&dates, reset)?;
-        let plan = prepare(
-            PrepareInput::Standard {
-                product: &product,
-                trade_date: START.parse()?,
-            },
-            &CollateralKey::Uncollateralized,
-            &[],
-            market.structure(),
-            &ctx,
-        )?;
-        assert_eq!(plan.warnings(), warnings.as_slice(), "{reset:?}");
-        let n0 = 150_000.0 / forward(value);
-        let n1 = 150_000.0 / forward(at(second_value, "06:00")?);
-        assert_close!(
-            evaluate(&plan, &market)?.get(Currency::USD),
-            usd(n0, n1),
-            1e-12_f64,
-            1e-8_f64
-        );
-    }
-    Ok(())
-}
-
-#[test]
-fn overnight_futures_compound_or_average_the_reference_period() -> Result<()> {
-    // 先物の清算指標を切り替えず、最終値の参照金利だけを複利・加重平均で組み立てる。
-    let as_of = AS_OF.parse()?;
-    let market = market(as_of, 0.0)?;
-    let config = config("")?;
-    let accrual = DateInterval::new("2026-01-05".parse()?, "2026-01-07".parse()?)?;
-    let tau = 1.0 / 365.0;
-    for (reference_rate, expected) in [
-        (
-            RateFutureRate::Compounded(accrual),
-            (RATE * 2.0 * tau).exp_m1() / (2.0 * tau),
-        ),
-        (
-            RateFutureRate::Averaged(accrual),
-            (RATE * tau).exp_m1() / tau,
-        ),
-    ] {
-        let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
-        if let Some(term) = parts.indices.get_mut(&"TERM".parse()?) {
-            term.kind = IndexKind::OvernightRate {
-                ccy: Currency::USD,
-                day_count: DayCount::Act365F,
-            };
-        }
-        if let Some(Convention::RateFuture(convention)) =
-            parts.conventions.get_mut(&"USD-RATE-FUTURES".parse()?)
-            && let Some(contract) = convention.contracts.values_mut().next()
-        {
-            contract.rate = reference_rate;
-        }
-        let reference = RefData::new(parts)?;
-        let ctx = Context {
-            as_of,
-            refdata: &reference,
-            config: &config,
-            fixings: &FixingStore::default(),
-        };
-        let df = (-RATE * model_time(as_of, at(START, "16:00")?).years()).exp();
-        let pv = futures_pv(&future(Position::Long, 0.95)?, &ctx, &market)?.0;
-        assert_close!(
-            pv,
-            (1.0 - expected - 0.95) * df * MULTIPLIER,
-            1e-12_f64,
-            1e-8_f64
-        );
-    }
-    Ok(())
-}
-
-#[test]
-fn overnight_futures_annualize_the_spanning_day_over_its_remaining_fraction() -> Result<()> {
-    // 評価時点 A が翌日物の観測の期間をまたぎ、その日の率が未公表のとき、[A, 終了) の forward の率を
-    // 残りの日数係数 τ_rem で年率化し、その率を元の 1 日の日数係数に当てる（参照データの章の
-    // fixing の規則）。成長を [A, 終了) に縮めたまま元の重みを掛けると、経過した分の利息が抜ける。
-    // 複利と加重平均のそれぞれを手計算と照合し、`futures_pv` で IR のリニア評価と解析の一致も確かめる。
-    let as_of = at(START, "20:00")?;
-    let market = market(as_of, 0.0)?;
-    let config = config("")?;
-    let accrual = DateInterval::new("2026-01-05".parse()?, "2026-01-07".parse()?)?;
-    let day = 1.0 / 365.0;
-    let first_end = at("2026-01-06", "16:00")?;
-    let second_end = at("2026-01-07", "16:00")?;
-    // 1 日目は A（20 時）から翌日の観測時刻（16 時）までの 20 時間だけが残る。
-    let remaining = 20.0 / 24.0 * day;
-    let first_growth = (RATE * model_time(as_of, first_end).years()).exp();
-    let second_growth = (RATE
-        * (model_time(as_of, second_end).years() - model_time(as_of, first_end).years()))
-    .exp();
-    let first_rate = (first_growth - 1.0) / remaining;
-    let second_rate = (second_growth - 1.0) / day;
-    for (reference_rate, expected) in [
-        (
-            RateFutureRate::Compounded(accrual),
-            ((1.0 + day * first_rate) * (1.0 + day * second_rate) - 1.0) / (2.0 * day),
-        ),
-        (
-            RateFutureRate::Averaged(accrual),
-            f64::midpoint(first_rate, second_rate),
-        ),
-    ] {
-        let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
-        if let Some(term) = parts.indices.get_mut(&"TERM".parse()?) {
-            term.kind = IndexKind::OvernightRate {
-                ccy: Currency::USD,
-                day_count: DayCount::Act365F,
-            };
-            // 翌日物は翌営業日の 16 時に公表する。A では 1 日目が未公表で、期間が A をまたぐ。
-            if let PublicationRule::Daily(rule) = &mut term.publication {
-                rule.business_days = 1;
-            }
-        }
-        if let Some(Convention::RateFuture(convention)) =
-            parts.conventions.get_mut(&"USD-RATE-FUTURES".parse()?)
-            && let Some(contract) = convention.contracts.values_mut().next()
-        {
-            contract.rate = reference_rate;
-        }
-        let reference = RefData::new(parts)?;
-        // 初日の清算価格は 16 時に公表済みなので、最新の決済価格 P_last として渡す。
-        let fixings = FixingStore::new(
-            &reference,
-            FixingStoreParts {
-                published: vec![(
-                    CLEARING.parse()?,
-                    ObservationPoint::Date(START.parse()?),
-                    0.954,
-                )],
-                ..FixingStoreParts::default()
-            },
-        )?;
-        let ctx = Context {
-            as_of,
-            refdata: &reference,
-            config: &config,
-            fixings: &fixings,
-        };
-        // 次の値洗いの受渡は翌日の 16 時。
-        let df = (-RATE * model_time(as_of, first_end).years()).exp();
-        let (pv, warnings) = futures_pv(&future(Position::Long, 0.95)?, &ctx, &market)?;
-        assert_close!(
-            pv,
-            (1.0 - expected - 0.954) * df * MULTIPLIER,
-            1e-12_f64,
-            1e-8_f64
-        );
-        assert!(
-            warnings
-                .iter()
-                .any(|warning| matches!(warning, Warning::UnpublishedFixingForecast { .. }))
-        );
-    }
-    Ok(())
-}
-
-fn check_gradient(
-    plan: &Prepared,
-    market: &Market<f64>,
-    block: &BlockId,
-    element: usize,
-    bump: f64,
-) -> Result<f64> {
-    let structure = market.structure();
-    let mut gradient = vec![0.0; market.theta().len()];
-    // 準備済みの商品をそのまま再利用し、AAD と中心差分の双方で同じマークを動かす。
-    Tape::default()
-        .gradient_into(
-            |theta| {
-                Market::from_theta(structure.clone(), theta.to_vec())
-                    .and_then(|market| evaluate(plan, &market))
-                    .map(|pv| pv.get(Currency::USD))
-                    .map_err(|error| {
-                        SciOxideError::invalid_input(
-                            argument!("market"),
-                            InvalidInputReason::Inconsistent,
-                            error.to_string(),
-                        )
-                    })
-            },
-            market.theta(),
-            &mut gradient,
-            ExecutionContext::sequential(),
-        )
-        .map_err(|error| Error::Numerical(error.to_string()))?;
-    let position = structure
-        .layout()
-        .range(block)
-        .ok_or_else(|| Error::Internal("missing AAD block".into()))?
-        .start
-        + element;
-    let evaluate_at = |delta: f64| -> Result<f64> {
-        let mut theta = market.theta().to_vec();
-        theta[position] += delta;
-        let market = Market::from_theta(structure.clone(), theta)?;
-        Ok(evaluate(plan, &market)?.get(Currency::USD))
-    };
-    let central = (evaluate_at(bump)? - evaluate_at(-bump)?) / (2.0 * bump);
-    assert_close!(gradient[position], central, 1e-7_f64, 5e-8_f64);
-    Ok(gradient[position])
-}
-
-#[test]
-fn futures_convexity_and_mtm_fx_marks_have_aad_gradients() -> Result<()> {
-    // LGM ボラと評価時点の FX マークをそれぞれ動かし、値洗いと元本リセットへ感応度が流れることを確かめる。
-    let as_of = AS_OF.parse()?;
-    let fixings = FixingStore::default();
-    let futures_market = market(as_of, 0.5)?;
-    let futures_reference = refdata("0D".parse()?, "16:00:00")?;
-    let futures_config = config("")?;
-    let futures_context = Context {
-        as_of,
-        refdata: &futures_reference,
-        config: &futures_config,
-        fixings: &fixings,
-    };
-    let product = future(Position::Long, 0.95)?;
-    let plan = prepare(
-        PrepareInput::Standard {
-            product: &product,
-            trade_date: START.parse()?,
-        },
-        &CollateralKey::Uncollateralized,
-        &[],
-        futures_market.structure(),
-        &futures_context,
-    )?;
-    assert!(
-        check_gradient(
-            &plan,
-            &futures_market,
-            &MarketStructure::lgm_block(Currency::USD)?,
-            1,
-            1e-3
-        )?
-        .abs()
-            > 1e-4
-    );
-
-    let fx = fx_market(as_of)?;
-    let reference = fx_reference("USDJPY".parse()?)?;
-    let fx_config = fx_config("")?;
-    let fx_context = Context {
-        as_of,
-        refdata: &reference,
-        config: &fx_config,
-        fixings: &fixings,
-    };
-    let swap = mtm_swap()?;
-    let plan = prepare(
-        PrepareInput::Standard {
-            product: &swap,
-            trade_date: START.parse()?,
-        },
-        &CollateralKey::Uncollateralized,
-        &[],
-        fx.structure(),
-        &fx_context,
-    )?;
-    // USD 脚にスプレッドを入れず同じ曲線で割り引くため、USD PV はほぼゼロ。支払いの連鎖を微分できることを照合する。
-    let _gradient = check_gradient(
-        &plan,
-        &fx,
-        &FxSpec::marks_block("USDJPY".parse()?)?,
-        0,
-        0.01,
-    )?;
-    Ok(())
-}
-
-#[test]
-fn fra_quote_calibrates_the_trade_to_par() -> Result<()> {
-    // スポットから指定した 0D×1W FRA の par quote で曲線を解き、同一条件の取引を評価する。
-    let as_of = AS_OF.parse()?;
-    let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
-    parts.conventions.insert(
-        "USD-FRA".parse()?,
-        Convention::Fra(FraConvention {
-            index: "TERM".parse()?,
-            spot_lag: "2BD".parse()?,
-        }),
-    );
-    let reference = RefData::new(parts)?;
-    let config = config("")?;
-    let context = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &FixingStore::default(),
-    };
-    let template = InstrumentTemplate::Fra {
-        convention: "USD-FRA".parse()?,
-        start: "0D".parse()?,
-        end: "1W".parse()?,
-    };
-    let quote_id = "USD.FRA.1W".parse()?;
-    let binding = InstrumentBinding::new(
-        quote_id,
-        QuoteKind::Rate,
-        CalibrationTarget::ParRate(template.clone()),
-        QuoteTransform::Identity,
-        CollateralKey::Uncollateralized,
-        ElementSelector::Block(BlockId::nodes("USD".parse()?)),
-    )?;
-    let definition =
-        MarketDefinition::new("USD".parse()?, "Etc/UTC".parse()?, [curve()?], [binding])?;
-    let tau = 7.0 / 365.0;
-    let quote = (RATE * tau).exp_m1() / tau;
-    let calibrated = calibrate(
-        &definition,
-        &QuoteSet::new([("USD.FRA.1W".parse()?, quote)])?,
-        &context,
-    )?;
-    let trade = template.build_for_trade(quote, START.parse()?, &reference)?;
-    let pv = price(&trade, calibrated.market(), &context)?;
-    assert_close!(pv.get(Currency::USD), 0.0_f64, 1e-12_f64, 1e-8_f64);
-    Ok(())
-}
-
-#[test]
-fn fx_outright_quote_calibrates_a_par_trade_in_both_currencies() -> Result<()> {
-    // outright の校正は base と quote の両方の PV を評価時点の FX で合算する。
-    let as_of = AS_OF.parse()?;
-    let pair: CcyPair = "USDJPY".parse()?;
-    let reference = fx_reference(pair)?;
-    let config = fx_config("")?;
-    let context = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &FixingStore::default(),
-    };
-    let template = InstrumentTemplate::FxForward {
-        convention: "USDJPY-FWD".parse()?,
-        tenor: "1W".parse()?,
-    };
-    let discount = CurveSpec::new(
-        "JPY".parse()?,
-        Currency::JPY,
-        [CurveRole::Discount],
-        [PILLAR.parse()?],
-        Extrapolation::Reject,
-    )?;
-    let usd = InstrumentBinding::new(
-        "USD.DF.1".parse()?,
-        QuoteKind::Mark,
-        CalibrationTarget::DirectMark,
-        QuoteTransform::Identity,
-        CollateralKey::Uncollateralized,
-        ElementSelector::Param(rhodium_market::structure::ParamId::new(
-            BlockId::nodes("USD".parse()?),
-            0,
-        )),
-    )?;
-    let jpy = InstrumentBinding::new(
-        "USDJPY.FWD.1W".parse()?,
-        QuoteKind::Price,
-        CalibrationTarget::Price(template.clone()),
-        QuoteTransform::Identity,
-        CollateralKey::Uncollateralized,
-        ElementSelector::Block(BlockId::nodes("JPY".parse()?)),
-    )?;
-    let definition = MarketDefinition::new(
-        "USD".parse()?,
-        "Etc/UTC".parse()?,
-        [curve()?, discount],
-        [usd, jpy],
-    )?
-    .with_fx([FxDefinition {
-        pair,
-        direct: Some(DirectFx {
-            market_collateral: "csa".parse()?,
-            quotes: FxQuotes {
-                spot: "USDJPY.SPOT".parse()?,
-                on: None,
-                tn: None,
-                sn: None,
-                forwards: vec![],
-            },
-            implied: None,
-        }),
-        collaterals: vec![FxCollateral {
-            id: "csa".parse()?,
-            base: "USD".parse()?,
-            quote: "JPY".parse()?,
-        }],
-        via: None,
-    }])?;
-    let quote = 150.0 * (-(RATE - 0.01) * 7.0 / 365.0).exp();
-    let calibrated = calibrate(
-        &definition,
-        &QuoteSet::new([
-            (
-                "USD.DF.1".parse()?,
-                -RATE * model_time(as_of, PILLAR.parse()?).years(),
-            ),
-            ("USDJPY.FWD.1W".parse()?, quote),
-            ("USDJPY.SPOT".parse()?, 150.0),
-        ])?,
-        &context,
-    )?;
-    let trade = template.build_for_trade(quote, START.parse()?, &reference)?;
-    let market = calibrated.market();
-    let pv = price(&trade, market, &context)?;
-    let fx = market.with_view(|view| view.fx_today(pair))?;
-    assert_close!(
-        pv.convert(Currency::JPY, &[fx])?.amount(),
-        0.0_f64,
-        1e-12_f64,
-        1e-6_f64
-    );
-    Ok(())
-}
-
-/// USD のピラーを直接マークで、JPY のピラーを 1W の USDJPY アウトライトで当てる市場の定義。
-/// `outrights` は JPY の曲線に当てるアウトライトの quote の ID と許容幅 σ（無ければ正方系）。
-fn outright_definition(outrights: &[(&str, Option<f64>)]) -> Result<MarketDefinition> {
-    let pair: CcyPair = "USDJPY".parse()?;
-    let template = InstrumentTemplate::FxForward {
-        convention: "USDJPY-FWD".parse()?,
-        tenor: "1W".parse()?,
-    };
-    let discount = CurveSpec::new(
-        "JPY".parse()?,
-        Currency::JPY,
-        [CurveRole::Discount],
-        [PILLAR.parse()?],
-        Extrapolation::Reject,
-    )?;
-    let usd = InstrumentBinding::new(
-        "USD.DF.1".parse()?,
-        QuoteKind::Mark,
-        CalibrationTarget::DirectMark,
-        QuoteTransform::Identity,
-        CollateralKey::Uncollateralized,
-        ElementSelector::Param(rhodium_market::structure::ParamId::new(
-            BlockId::nodes("USD".parse()?),
-            0,
-        )),
-    )?;
-    let jpy = outrights
-        .iter()
-        .map(|&(quote, sigma)| {
-            let binding = InstrumentBinding::new(
-                quote.parse()?,
-                QuoteKind::Price,
-                CalibrationTarget::Price(template.clone()),
-                QuoteTransform::Identity,
-                CollateralKey::Uncollateralized,
-                ElementSelector::Block(BlockId::nodes("JPY".parse()?)),
-            )?;
-            sigma.map_or(Ok(binding.clone()), |sigma| binding.with_tolerance(sigma))
-        })
-        .collect::<Result<Vec<_>>>()?;
-    MarketDefinition::new(
-        "USD".parse()?,
-        "Etc/UTC".parse()?,
-        [curve()?, discount],
-        std::iter::once(usd).chain(jpy),
-    )?
-    .with_fx([FxDefinition {
-        pair,
-        direct: Some(DirectFx {
-            market_collateral: "csa".parse()?,
-            quotes: FxQuotes {
-                spot: "USDJPY.SPOT".parse()?,
-                on: None,
-                tn: None,
-                sn: None,
-                forwards: vec![],
-            },
-            implied: None,
-        }),
-        collaterals: vec![FxCollateral {
-            id: "csa".parse()?,
-            base: "USD".parse()?,
-            quote: "JPY".parse()?,
-        }],
-        via: None,
-    }])
-}
-
-#[test]
-#[expect(
-    clippy::wildcard_enum_match_arm,
-    reason = "#196 で各 variant を網羅する"
-)]
-fn fx_outright_least_squares_measures_residuals_in_price_units() -> Result<()> {
-    // アウトライトの許容幅 σ は価格の単位なので、最小二乗の残差も価格の差 F − q で測る。
-    // 同じ受渡の二つのアウトライト q₁、q₂ = q₁ + 2 を σ₁ = 0.01、σ₂ = 1 で当てると、
-    // F* = (q₁/σ₁² + q₂/σ₂²)/(1/σ₁² + 1/σ₂²) = q₁ + 2/10001 で、残差は F* − q₁ ≈ 0.0002
-    // （σ₁ の 0.02 倍で内）と F* − q₂ ≈ −2（σ₂ の 2 倍で外）。外れた q₂ だけを警告する。
-    // 残差を ln(F/q) のまま 1/σ で重みを付けると、q₂ の外れは 0.013 と測られて警告されない。
-    // quote デルタ dθ/dq は、残差の quote 微分を −1 にした随伴が再校正の中心差分と一致する。
-    let as_of = AS_OF.parse()?;
-    let reference = fx_reference("USDJPY".parse()?)?;
-    let config = fx_config("")?;
-    let context = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &FixingStore::default(),
-    };
-    let definition = outright_definition(&[
-        ("USDJPY.FWD.1W.A", Some(0.01)),
-        ("USDJPY.FWD.1W.B", Some(1.0)),
-    ])?;
-    let q1 = 150.0 * (-(RATE - 0.01) * 7.0 / 365.0).exp();
-    let q2 = q1 + 2.0;
-    let quotes = |q1: f64| {
-        QuoteSet::new([
-            (
-                "USD.DF.1".parse()?,
-                -RATE * model_time(as_of, PILLAR.parse()?).years(),
-            ),
-            ("USDJPY.FWD.1W.A".parse()?, q1),
-            ("USDJPY.FWD.1W.B".parse()?, q2),
-            ("USDJPY.SPOT".parse()?, 150.0),
-        ])
-    };
-    let calibrated = calibrate(&definition, &quotes(q1)?, &context)?;
-    let fitted = (q1 * 1e4 + q2) / (1e4 + 1.0);
-    // F* で約定したフォワードは、報告通貨で合算した PV がゼロ。
-    let template = InstrumentTemplate::FxForward {
-        convention: "USDJPY-FWD".parse()?,
-        tenor: "1W".parse()?,
-    };
-    let trade = template.build_for_trade(fitted, START.parse()?, &reference)?;
-    let market = calibrated.market();
-    let fx = market.with_view(|view| view.fx_today("USDJPY".parse()?))?;
-    assert_close!(
-        price(&trade, market, &context)?
-            .convert(Currency::JPY, &[fx])?
-            .amount(),
-        0.0_f64,
-        1e-12_f64,
-        1e-6_f64
-    );
-    let jpy = BlockId::nodes("JPY".parse()?);
-    let group = calibrated
-        .groups()
-        .iter()
-        .position(|group| group.blocks().contains(&jpy))
-        .ok_or_else(|| Error::Internal("missing JPY group".into()))?;
-    let GroupOutcome::Solved(solution) = calibrated.groups()[group].outcome() else {
-        return Err(Error::Calibration("JPY group not solved".into()));
-    };
-    assert_close!(solution.residuals()[0], fitted - q1, 1e-9_f64, 1e-9_f64);
-    assert_close!(solution.residuals()[1], fitted - q2, 1e-9_f64, 1e-9_f64);
-    let outside: Vec<_> = calibrated
-        .warnings()
-        .iter()
-        .filter_map(|warning| match warning {
-            Warning::CalibrationFitOutsideTolerance { quote } => Some(quote.as_str()),
-            _ => None,
-        })
-        .collect();
-    assert_eq!(outside, ["USDJPY.FWD.1W.B"]);
-
-    // dθ/dq₁ = −(∂G/∂θ)⁻¹·∂G/∂q₁（列は quote の ID の昇順で A が先）を、q₁ を動かした再校正と比べる。
-    let inverse = solution
-        .dg_dtheta()
-        .solve(&[1.0][..])
-        .map_err(|error| Error::Numerical(error.to_string()))?
-        .as_slice()[0];
-    let dg_dq1 = solution
-        .dg_dq()
-        .get(0, 0)
-        .map_err(|error| Error::Numerical(error.to_string()))?;
-    let position = calibrated
-        .market()
-        .structure()
-        .layout()
-        .range(&jpy)
-        .ok_or_else(|| Error::Internal("missing JPY block".into()))?
-        .start;
-    let bump = 1e-3;
-    let theta = |q1: f64| -> Result<f64> {
-        Ok(calibrate(&definition, &quotes(q1)?, &context)?
-            .market()
-            .theta()[position])
-    };
-    let central = (theta(q1 + bump)? - theta(q1 - bump)?) / (2.0 * bump);
-    assert_close!(-inverse * dg_dq1, central, 1e-6_f64, 1e-12_f64);
-    Ok(())
-}
-
-#[test]
-fn fixed_notional_xccy_exchanges_both_currencies() -> Result<()> {
-    // 元本を固定した XCCY は、各通貨の脚が期首に元本を払い、期末に元本を受け取り、クーポンを
-    // それぞれの通貨の曲線で割り引く（FX のリセットは無い）。JPY の脚は 1% の固定、USD の脚は
-    // ターム金利の 1W（RATE の連続金利から求めた単利）にスプレッド 20bp。
-    let as_of = AS_OF.parse()?;
-    let market = fx_market(as_of)?;
-    let config = fx_config("")?;
-    let reference = fx_reference("USDJPY".parse()?)?;
-    let ctx = Context {
-        as_of,
-        refdata: &reference,
-        config: &config,
-        fixings: &FixingStore::default(),
-    };
-    let schedule = explicit_schedule()?;
-    let swap = StandardProduct::Swap(Swap::from_legs([
-        Leg::new(
-            PayReceive::Receive,
-            Money::new(Currency::JPY, 150_000.0),
-            schedule.clone(),
-            DayCount::Act365F,
-            Coupon::Fixed {
-                rate: Rate::new(0.01, DayCount::Act365F, Simple)?,
-            },
-            NotionalExchange::InitialAndFinal,
-            PaymentRounding::None,
-        )?,
-        Leg::new(
-            PayReceive::Pay,
-            Money::new(Currency::USD, 1_000.0),
-            schedule,
-            DayCount::Act365F,
-            term("TERM", 0.002)?,
-            NotionalExchange::InitialAndFinal,
-            PaymentRounding::None,
-        )?,
-    ])?);
-    let plan = prepare(
-        PrepareInput::Standard {
-            product: &swap,
-            trade_date: START.parse()?,
-        },
-        &CollateralKey::Uncollateralized,
-        &[],
-        market.structure(),
-        &ctx,
-    )?;
-    // 元本が固定なので、FX と金利の積の近似（MtmFxForwardProduct）は使わない。
-    assert!(!plan.warnings().contains(&Warning::MtmFxForwardProduct));
-    let pv = evaluate(&plan, &market)?;
-    let discount = |rate: f64, date: &str| -> Result<f64> {
-        Ok((-rate * model_time(as_of, at(date, "16:00")?).years()).exp())
-    };
-    let tau = 7.0 / 365.0;
-    let [start, middle, end] = ["2026-01-07", "2026-01-14", "2026-01-21"];
-    let jpy = 150_000.0
-        * (-discount(0.01, start)?
-            + 0.01 * tau * discount(0.01, middle)?
-            + (1.0 + 0.01 * tau) * discount(0.01, end)?);
-    let coupon = (RATE * tau).exp_m1() + 0.002 * tau;
-    let usd = 1_000.0
-        * (discount(RATE, start)?
-            - coupon * discount(RATE, middle)?
-            - (1.0 + coupon) * discount(RATE, end)?);
-    assert_close!(pv.get(Currency::JPY), jpy, 1e-12_f64, 1e-6_f64);
-    assert_close!(pv.get(Currency::USD), usd, 1e-12_f64, 1e-8_f64);
-    Ok(())
-}
-
-#[test]
-fn term_coupon_payment_delay_is_recorded() -> Result<()> {
-    // ターム金利のクーポンを期末から 2 日遅れて払うと、自然な支払日（期末）とのずれのタイミング調整を
-    // 省く近似（PaymentDelayIgnored）を記録し、予想は期末の DF 比、割引は実際の支払日で求める。
-    // 近似を禁止すると準備できない。
-    let as_of = AS_OF.parse()?;
-    let market = market(as_of, 0.0)?;
-    let reference = refdata("0D".parse()?, "16:00:00")?;
-    let fixings = FixingStore::default();
-    let schedule = ScheduleSpec::Explicit(ExplicitPeriods::new(
-        vec![ExplicitPeriod {
-            accrual: DateInterval::new("2026-01-07".parse()?, "2026-01-14".parse()?)?,
-            payment: "2026-01-16".parse()?,
-            reset: None,
-            settlement: None,
-            references: vec![],
-        }],
-        ExplicitDatePolicy::KeepAsConfirmed,
-        PeriodCombination::Sum,
-        None,
-    )?);
-    let swap = StandardProduct::Swap(Swap::from_legs([
-        Leg::new(
-            PayReceive::Receive,
-            Money::new(Currency::USD, 100_000.0),
-            schedule.clone(),
-            DayCount::Act365F,
-            term("TERM", 0.0)?,
-            NotionalExchange::None,
-            PaymentRounding::None,
-        )?,
-        Leg::new(
-            PayReceive::Pay,
-            Money::new(Currency::USD, 100_000.0),
-            schedule,
-            DayCount::Act365F,
-            Coupon::Fixed {
-                rate: Rate::new(0.0, DayCount::Act365F, Simple)?,
-            },
-            NotionalExchange::None,
-            PaymentRounding::None,
-        )?,
-    ])?);
-    let prepare_with = |config: &PricingConfig| {
-        prepare(
-            PrepareInput::Standard {
-                product: &swap,
-                trade_date: START.parse()?,
-            },
-            &CollateralKey::Uncollateralized,
-            &[],
-            market.structure(),
-            &Context {
-                as_of,
-                refdata: &reference,
-                config,
-                fixings: &fixings,
-            },
-        )
-    };
-    let plan = prepare_with(&config("")?)?;
-    assert!(plan.warnings().contains(&Warning::PaymentDelayIgnored));
-    let tau = 7.0 / 365.0;
-    let forward = (RATE * tau).exp_m1() / tau;
-    let df = (-RATE * model_time(as_of, at("2026-01-16", "16:00")?).years()).exp();
-    assert_close!(
-        evaluate(&plan, &market)?.get(Currency::USD),
-        100_000.0 * tau * forward * df,
-        1e-12_f64,
-        1e-8_f64
-    );
-    assert!(matches!(
-        prepare_with(&config(
-            "[conventions]\ndisabled = ['payment_delay_ignored']"
-        )?),
-        Err(Error::Unsupported(_))
-    ));
-    Ok(())
-}
-
-/// [`fx_market`] に、JPY を 25bp 上乗せした割引曲線 `JPY-25` を足した市場。
-fn fx_market_with_spread_curve(as_of: Instant, pair: CcyPair) -> Result<Market<f64>> {
-    let jpy = CurveSpec::new(
-        "JPY".parse()?,
-        Currency::JPY,
-        [CurveRole::Discount, CurveRole::Projection],
-        [PILLAR.parse()?],
-        Extrapolation::Reject,
-    )?;
-    let spread = CurveSpec::from_expr(
-        "JPY-25".parse()?,
-        Currency::JPY,
-        [CurveRole::Discount],
-        rhodium_market::curve::CurveExpr::Sum(vec![
-            rhodium_market::curve::CurveExpr::Ref("JPY".parse()?),
-            rhodium_market::curve::CurveExpr::Flat(0.0025),
-        ]),
-    )?;
-    let structure = Arc::new(MarketStructure::with_fx(
-        as_of,
-        [curve()?, jpy, spread],
-        [FxSpec {
-            pair,
-            anchor: FxAnchor::Direct {
-                spot: as_of,
-                short_starts: vec![],
-                market_collateral: "csa".parse()?,
-            },
-            collaterals: vec![FxCollateral {
-                id: "csa".parse()?,
-                base: "USD".parse()?,
-                quote: "JPY".parse()?,
-            }],
-            via: None,
-        }],
-    )?);
-    let t = model_time(as_of, PILLAR.parse()?).years();
-    let mut theta = vec![0.0; structure.layout().len()];
-    for (id, value) in [
-        (BlockId::nodes("USD".parse()?), -RATE * t),
-        (BlockId::nodes("JPY".parse()?), -0.01 * t),
-        (FxSpec::marks_block(pair)?, 150.0),
-    ] {
-        let range = structure
-            .layout()
-            .range(&id)
-            .ok_or_else(|| Error::Internal("missing FX test block".into()))?;
-        theta[range.start] = value;
-    }
-    let market = Market::from_theta(structure, theta)?;
-    Ok(market)
-}
-
-#[test]
-fn fx_fixing_under_another_collateral_records_cross_collateral_projection() -> Result<()> {
-    // FX の予想は、ペアの FX の担保区分（USD と JPY の曲線の組）の下で定まる。支払いを別の担保
-    // （JPY を 25bp 上乗せした曲線で割り引く CSA）で受け取ると、予想と割引の担保が違うので
-    // CrossCollateralProjection を記録し、禁止すると準備できない。
-    let as_of = AS_OF.parse()?;
-    let reference = fx_reference("USDJPY".parse()?)?;
-    let fixings = FixingStore::default();
-    let pair: CcyPair = "USDJPY".parse()?;
-    let market = fx_market_with_spread_curve(as_of, pair)?;
-    let product = Product::new(
-        Contract::build(|builder| {
-            let rate = builder.fixing(
-                &"FX-FIX".parse()?,
-                ObservationPoint::Date("2026-01-12".parse()?),
-            );
-            Ok(builder.pay(Asset::Cash(Currency::JPY), rate, at("2026-01-14", "16:00")?))
-        })?,
-        &reference,
-    )?;
-    let collateral = CollateralKey::Cash(rhodium_products::agreement::CashCollateralTerms {
-        ccy: Currency::JPY,
-        remuneration: rhodium_products::agreement::Remuneration {
-            index: "JPY-TERM".parse()?,
-            spread: Spread::new(0.0025)?,
-            floor: None,
-        },
-    });
-    let rules = r#"
-[[discounting]]
-currency = "JPY"
-collateral = { cash = { ccy = "JPY", index = "JPY-TERM", spread = "25bp" } }
-curve = "JPY-25"
-[[discounting]]
-currency = "USD"
-collateral = { cash = { ccy = "JPY", index = "JPY-TERM", spread = "25bp" } }
-curve = "USD"
-"#;
-    let prepare_under = |config: &PricingConfig, collateral: &CollateralKey| {
-        prepare(
-            PrepareInput::Ir(&product),
-            collateral,
-            &[],
-            market.structure(),
-            &Context {
-                as_of,
-                refdata: &reference,
-                config,
-                fixings: &fixings,
-            },
-        )
-    };
-    let native = prepare_under(&fx_config("")?, &CollateralKey::Uncollateralized)?;
-    assert!(
-        !native
-            .warnings()
-            .contains(&Warning::CrossCollateralProjection)
-    );
-    let spread_config = fx_config(rules)?;
-    let foreign = prepare_under(&spread_config, &collateral)?;
-    assert!(
-        foreign
-            .warnings()
-            .contains(&Warning::CrossCollateralProjection)
-    );
-    // 予想は担保区分の曲線のフォワードのまま、割引だけが 25bp の曲線に変わる。
-    let ratio = evaluate(&foreign, &market)?.get(Currency::JPY)
-        / evaluate(&native, &market)?.get(Currency::JPY);
-    let pay = model_time(as_of, at("2026-01-14", "16:00")?).years();
-    assert_close!(ratio, (-0.0025 * pay).exp(), 1e-12_f64, 1e-14_f64);
-    let strict = fx_config(&format!(
-        "{rules}[conventions]\ndisabled = ['cross_collateral_projection']"
-    ))?;
-    assert!(matches!(
-        prepare_under(&strict, &collateral),
-        Err(Error::Unsupported(_))
-    ));
-    Ok(())
-}
-
-#[test]
-fn futures_margin_index_cross_collateral_projection_is_recorded() -> Result<()> {
-    // 値洗いの口座の付利指標（margin_index）の予想が値洗いの割引と等しい率なら、値洗いの列の評価は
-    // 厳密（Q_M の測度と割引が一致）。付利指標を別の曲線で予想すると、その差を無視する近似として
-    // CrossCollateralProjection を記録する。
-    let as_of = AS_OF.parse()?;
-    let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
-    let mut margin = parts.indices[&"TERM".parse()?].clone();
-    margin.id = "MARGIN".parse()?;
-    margin.kind = IndexKind::OvernightRate {
-        ccy: Currency::USD,
-        day_count: DayCount::Act360,
-    };
-    parts.indices.insert(margin.id.clone(), margin);
-    if let Some(Convention::RateFuture(convention)) =
-        parts.conventions.get_mut(&"USD-RATE-FUTURES".parse()?)
-    {
-        convention.margin_index = Some("MARGIN".parse()?);
-    }
-    let reference = RefData::new(parts)?;
-    let alternative = CurveSpec::new(
-        "USD-ALT".parse()?,
-        Currency::USD,
-        [CurveRole::Projection],
-        [PILLAR.parse()?],
-        Extrapolation::Reject,
-    )?;
-    let structure = Arc::new(MarketStructure::with_fx_and_lgm(
-        as_of,
-        [curve()?, alternative],
-        [],
-        [(Currency::USD, vec![])],
-    )?);
-    let mut theta = vec![0.0; structure.layout().len()];
-    let t = model_time(as_of, PILLAR.parse()?).years();
-    for curve in ["USD", "USD-ALT"] {
-        let range = structure
-            .layout()
-            .range(&BlockId::nodes(curve.parse()?))
-            .ok_or_else(|| Error::Internal("missing curve block".into()))?;
-        theta[range.start] = -RATE * t;
-    }
-    let market = Market::from_theta(structure, theta)?;
-    let fixings = FixingStore::default();
-    let warnings = |margin_curve: &str| -> Result<Vec<Warning>> {
-        let config = config(&format!("MARGIN = \"{margin_curve}\""))?;
-        let ctx = Context {
-            as_of,
-            refdata: &reference,
-            config: &config,
-            fixings: &fixings,
-        };
-        Ok(futures_pv(&future(Position::Long, 0.95)?, &ctx, &market)?.1)
-    };
-    assert!(!warnings("USD")?.contains(&Warning::CrossCollateralProjection));
-    assert!(warnings("USD-ALT")?.contains(&Warning::CrossCollateralProjection));
-    Ok(())
-}
-
-/// USD の予想曲線を先物だけで、別の割引曲線 USD-OIS を直接マークで当てる市場の定義。
-fn futures_with_discount_definition() -> Result<MarketDefinition> {
-    let discount = CurveSpec::new(
-        "USD-OIS".parse()?,
-        Currency::USD,
-        [CurveRole::Discount],
-        [PILLAR.parse()?],
-        Extrapolation::Reject,
-    )?;
-    let mark = InstrumentBinding::new(
-        "USD-OIS.DF".parse()?,
-        QuoteKind::Mark,
-        CalibrationTarget::DirectMark,
-        QuoteTransform::Identity,
-        CollateralKey::Uncollateralized,
-        ElementSelector::Param(rhodium_market::structure::ParamId::new(
-            BlockId::nodes("USD-OIS".parse()?),
-            0,
-        )),
-    )?;
-    let future = InstrumentBinding::new(
-        "USD.FUT.2026-01".parse()?,
-        QuoteKind::FuturesPrice,
-        CalibrationTarget::Price(futures_template()?),
-        QuoteTransform::Identity,
-        CollateralKey::Uncollateralized,
-        ElementSelector::Block(BlockId::nodes("USD".parse()?)),
-    )?;
-    MarketDefinition::new(
-        "USD".parse()?,
-        "Etc/UTC".parse()?,
-        [curve()?, discount],
-        [future, mark],
-    )?
-    .with_lgm([(
-        Currency::USD,
-        rhodium_pricing::definition::LgmDefinition {
-            breaks: vec![],
-            mean_reversion: vec!["USD.LGM.MR".parse()?],
-            volatility: vec!["USD.LGM.VOL".parse()?],
-        },
-    )])
-}
-
-#[test]
-fn futures_identify_the_forecast_measure_only_through_their_margin_index() -> Result<()> {
-    // 先物の価格は値洗いの口座の測度で決まり、取引の担保の割引因子は打ち消し合うので、先物で当てた
-    // 予想曲線の予想の担保は、binding の担保からは分からない。USD の予想曲線を先物だけで当て、
-    // 無担保の行で別の割引曲線 USD-OIS に割り引くとき:
-    // - 付利指標の無い先物は根拠にならず測度は不明なので、ターム金利のクーポンは予想と割引の担保の
-    //   違いを決定論的なベーシスの仮定で扱い、CrossCollateralProjection を記録する。
-    // - 付利指標 MARGIN を USD-OIS で予想する先物は、口座の測度が USD-OIS の割引だと示すので、
-    //   同じ割引のクーポンは担保の調整を要しない。
-    let as_of = AS_OF.parse()?;
-    let definition = futures_with_discount_definition()?;
-    let config: PricingConfig = toml::from_str(
-        r#"
-version = 1
-reporting_currency = "USD"
-[[discounting]]
-currency = "USD"
-collateral = "uncollateralized"
-curve = "USD-OIS"
-[forecasting]
-TERM = "USD"
-MARGIN = "USD-OIS"
-"#,
-    )
-    .map_err(|error| Error::InvalidData(error.to_string()))?;
-    let quotes = QuoteSet::new([
-        ("USD.FUT.2026-01".parse()?, 0.96),
-        ("USD.LGM.MR".parse()?, 0.0),
-        ("USD.LGM.VOL".parse()?, 0.0),
-        (
-            "USD-OIS.DF".parse()?,
-            -RATE * model_time(as_of, PILLAR.parse()?).years(),
-        ),
-    ])?;
-    let fixings = FixingStore::default();
-    let warnings = |margin: bool| -> Result<Vec<Warning>> {
-        let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
-        let mut index = parts.indices[&"TERM".parse()?].clone();
-        index.id = "MARGIN".parse()?;
-        index.kind = IndexKind::OvernightRate {
-            ccy: Currency::USD,
-            day_count: DayCount::Act360,
-        };
-        parts.indices.insert(index.id.clone(), index);
-        if margin
-            && let Some(Convention::RateFuture(convention)) =
-                parts.conventions.get_mut(&"USD-RATE-FUTURES".parse()?)
-        {
-            convention.margin_index = Some("MARGIN".parse()?);
-        }
-        let reference = RefData::new(parts)?;
-        let ctx = Context {
-            as_of,
-            refdata: &reference,
-            config: &config,
-            fixings: &fixings,
-        };
-        let market = calibrate(&definition, &quotes, &ctx)?;
-        // 01-07 に観測する 1W のターム金利を期末の 01-14 に受け取るクーポン。
-        let term: IndexId = "TERM".parse()?;
-        let coupon = Product::new(
-            Contract::build(|builder| {
-                let rate = builder.fixing(&term, ObservationPoint::Date(FINAL.parse()?));
-                Ok(builder.pay(Asset::Cash(Currency::USD), rate, at("2026-01-14", "16:00")?))
-            })?,
-            &reference,
-        )?;
-        Ok(prepare(
-            PrepareInput::Ir(&coupon),
-            &CollateralKey::Uncollateralized,
-            &[],
-            market.market().structure(),
-            &ctx,
-        )?
-        .warnings()
-        .to_vec())
-    };
-    assert_eq!(warnings(false)?, [Warning::CrossCollateralProjection]);
-    assert_eq!(warnings(true)?, [] as [Warning; 0]);
-    Ok(())
-}
diff --git a/crates/rhodium-pricing/tests/ois_calibration.rs b/crates/rhodium-pricing/tests/ois_calibration.rs
index e43b1910acd5b816b10ae143f0631867294713ee..de077a9b858d5ac0ed8d585a06756785b71ae3fb 100644
--- a/crates/rhodium-pricing/tests/ois_calibration.rs
+++ b/crates/rhodium-pricing/tests/ois_calibration.rs
@@ -11,6 +11,7 @@ use rhodium_core::{
     rounding::{Rounding, RoundingMode},
     schedule::EndRule,
 };
+use rhodium_fixtures::{multicurve, ois, steps::PenalizedTona};
 use rhodium_market::{
     curve::CurveId,
     market::Market,
@@ -32,16 +33,6 @@ use rhodium_products::{
 use rhodium_testkit::assert_close;
 use scioxide_core::AdjointOperator;
 
-#[path = "fixtures/ois.rs"]
-pub mod ois;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、このテストで使わない関数がある"
-)]
-#[path = "fixtures/multicurve.rs"]
-pub mod multicurve;
-
 use ois::{Case, Inputs, cases, definition, fixture_quotes, midnight, prepare_ois};
 
 fn par_rates(
@@ -916,10 +907,6 @@ fn short_ois_calibrate_to_the_default_tolerance() -> Result<()> {
 }
 
 #[test]
-#[expect(
-    clippy::wildcard_enum_match_arm,
-    reason = "#196 で各 variant を網羅する"
-)]
 fn rate_quotes_above_one_warn_about_a_missing_unit() -> Result<()> {
     // 金利の binding が読む quote の絶対値が 1 を超えたら、校正の結果に quote ID 付きで警告する。
     // 5% を単位なしの 5 と書いたような誤りは、値だけでは 500% として校正されてしまうため。
@@ -947,15 +934,14 @@ fn rate_quotes_above_one_warn_about_a_missing_unit() -> Result<()> {
             (id, value)
         });
         let result = calibrate(&def, &QuoteSet::new(quotes)?, &ctx)?;
-        let warned: Vec<_> = result
-            .warnings()
+        // 警告の一覧全体を比べ、単位の疑いの警告だけがあることを確かめる。
+        let expected: Vec<_> = suspected
             .iter()
-            .filter_map(|warning| match warning {
-                Warning::QuoteUnitSuspected { quote } => Some(quote.clone()),
-                _ => None,
+            .map(|quote| Warning::QuoteUnitSuspected {
+                quote: quote.clone(),
             })
             .collect();
-        assert_eq!(warned, suspected, "{values:?}");
+        assert_eq!(result.warnings(), expected, "{values:?}");
     }
     Ok(())
 }
@@ -1107,91 +1093,16 @@ fn calibration_reports_instruments_beyond_the_curve_as_input_errors() -> Result<
     Ok(())
 }
 
-#[derive(serde::Deserialize)]
-struct QuantLibTurn {
-    effective: rhodium_core::date::Date,
-    spread: f64,
-    fraction: f64,
-    nodes: Vec<ois::Node>,
-    samples: Vec<ois::Sample>,
-}
-
-fn tona_with_steps_and_turn(
-    tona: &Case,
-    turn: &QuantLibTurn,
-) -> Result<(
-    rhodium_market::curve::CurveSpec,
-    Vec<rhodium_pricing::definition::InstrumentBinding>,
-)> {
-    use rhodium_core::quantity::YearFraction;
-    use rhodium_market::{
-        curve::{CurveExpr, CurveRole, CurveSpec, Extrapolation},
-        quote::QuoteKind,
-        structure::BlockId,
-    };
-    use rhodium_pricing::definition::{CalibrationTarget, InstrumentBinding, QuoteTransform};
-
-    let id: CurveId = tona.index.parse()?;
-    let curve = CurveSpec::from_expr(
-        id.clone(),
-        Currency::JPY,
-        [CurveRole::Discount, CurveRole::Projection],
-        CurveExpr::Sum(vec![
-            CurveExpr::Stitch {
-                short: Box::new(CurveExpr::Steps {
-                    part: "meetings".into(),
-                    meetings: vec![midnight(turn.nodes[0].date)?],
-                }),
-                long: Box::new(CurveExpr::Nodes {
-                    part: "long".into(),
-                    pillars: turn.nodes[2..]
-                        .iter()
-                        .map(|node| midnight(node.date))
-                        .collect::<Result<_>>()?,
-                    extrapolation: Extrapolation::Reject,
-                }),
-                at: midnight(turn.nodes[1].date)?,
-            },
-            CurveExpr::Jumps {
-                part: "turn".into(),
-                dates: vec![midnight(turn.effective)?],
-            },
-        ]),
-    )?;
-    let collateral = CollateralKey::Uncollateralized;
-    let mut bindings = turn
-        .nodes
-        .iter()
-        .enumerate()
-        .map(|(position, node)| {
-            let part = if position < 2 { "meetings" } else { "long" };
-            tona.binding(
-                node,
-                &collateral,
-                ElementSelector::Block(BlockId::new(id.clone(), part)?),
-            )
-        })
-        .collect::<Result<Vec<_>>>()?;
-    // 年末のスプレッドは直接マークで、変換 TurnSpreadToJump が s·τ_turn の跳びにする。
-    bindings.push(InstrumentBinding::new(
-        "JPY.TONA.TURN".parse()?,
-        QuoteKind::Spread,
-        CalibrationTarget::DirectMark,
-        QuoteTransform::TurnSpreadToJump(YearFraction::new(turn.fraction, tona.basis())?),
-        collateral,
-        ElementSelector::Param(ParamId::new(BlockId::new(id, "turn")?, 0)),
-    )?);
-    Ok((curve, bindings))
-}
-
 #[test]
 fn meeting_steps_and_turn_reprice_ois_quotes() -> Result<()> {
     // QuantLib の独立した OIS ブートストラップと、会合日を含む同じ区間の DF を照合する。
-    let [mut tona, sofr] = cases()?;
-    let turn: QuantLibTurn = serde_json::from_str(include_str!("fixtures/quantlib/ois_turn.json"))
-        .map_err(|error| Error::InvalidData(error.to_string()))?;
-    let (curve, mut bindings) = tona_with_steps_and_turn(&tona, &turn)?;
-    tona.nodes.clone_from(&turn.nodes);
+    // TONA は会合日の階段と年末の跳びを持つ正方系の曲線（fixture の `ois_turn.json` の quote）で、
+    // SOFR は通常の OIS 曲線で当てる。
+    let [_, sofr] = cases()?;
+    let penalized = PenalizedTona::new()?;
+    let (tona, turn) = (&penalized.case, &penalized.turn);
+    let curve = penalized.square_curve()?;
+    let mut bindings = penalized.square_bindings()?;
     let collateral = CollateralKey::Uncollateralized;
     let sofr_block = sofr.block()?;
     for node in &sofr.nodes {
@@ -1209,12 +1120,13 @@ fn meeting_steps_and_turn_reprice_ois_quotes() -> Result<()> {
         bindings,
     )?;
     let quotes = QuoteSet::new(
-        fixture_quotes(&tona)?
-            .into_iter()
-            .chain(fixture_quotes(&sofr)?)
-            .chain([("JPY.TONA.TURN".parse()?, turn.spread)]),
+        penalized
+            .quotes()?
+            .iter()
+            .map(|(id, value)| (id.clone(), value))
+            .chain(fixture_quotes(&sofr)?),
     )?;
-    let inputs = Inputs::new(&[&tona, &sofr], tona.as_of, "")?;
+    let inputs = Inputs::new(&[tona, &sofr], tona.as_of, "")?;
     let context = inputs.context();
     let result = calibrate(&definition, &quotes, &context)?;
     for group in result.groups() {
@@ -1243,7 +1155,7 @@ fn meeting_steps_and_turn_reprice_ois_quotes() -> Result<()> {
         }
         Ok(())
     })?;
-    for case in [&tona, &sofr] {
+    for case in [tona, &sofr] {
         for (actual, node) in par_rates(case, &collateral, market, &context)?
             .into_iter()
             .zip(&case.nodes)
diff --git a/crates/rhodium-pricing/tests/ois_futures.rs b/crates/rhodium-pricing/tests/ois_futures.rs
new file mode 100644
index 0000000000000000000000000000000000000000..361cd2105366af5d4e195444c4d6585836591c34
--- /dev/null
+++ b/crates/rhodium-pricing/tests/ois_futures.rs
@@ -0,0 +1,153 @@
+//! 翌日物の金利を参照する金利先物の最終値を、参照期間の複利と加重平均の手計算で照合する。
+//! 評価時点が参照期間の観測をまたぐときの、残りの日数係数による年率化も確かめる。
+
+pub mod support;
+
+use rhodium_core::{
+    Error, Result, Warning,
+    currency::Currency,
+    date::DateInterval,
+    day_count::DayCount,
+    fixing::{FixingStore, FixingStoreParts, ObservationPoint},
+    refdata::{IndexKind, PublicationRule, RateFutureRate, RefData},
+    time::model_time,
+};
+use rhodium_fixtures::futures;
+use rhodium_pricing::context::Context;
+use rhodium_products::standard::Position;
+use rhodium_testkit::{PV_REL, assert_close};
+
+use futures::{AS_OF, CLEARING, MULTIPLIER, RATE, START, at, config, market, reference_parts};
+use support::{USD_ABS, future, futures_pv, index_definition, january_contract};
+
+/// 先物の参照データの `TERM` を Act/365F の翌日物の指標にし、観測日から `publication_lag` 営業日後の
+/// 16 時に公表して、2026 年 1 月限の最終値の参照金利を `rate` にした参照データ。先物の清算指標は
+/// 切り替えない。
+fn overnight_reference(rate: RateFutureRate, publication_lag: u16) -> Result<RefData> {
+    let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
+    let term = index_definition(&mut parts, "TERM")?;
+    term.kind = IndexKind::OvernightRate {
+        ccy: Currency::USD,
+        day_count: DayCount::Act365F,
+    };
+    // 公表が日次でなくなったら、公表のラグを書き換えずに既定の場面で進まないようエラーにする。
+    let PublicationRule::Daily(rule) = &mut term.publication else {
+        return Err(Error::InvalidData(
+            "fixture index TERM is not published daily".into(),
+        ));
+    };
+    rule.business_days = publication_lag;
+    january_contract(&mut parts)?.rate = rate;
+    RefData::new(parts)
+}
+
+#[test]
+fn overnight_futures_compound_or_average_the_reference_period() -> Result<()> {
+    // 先物の清算指標を切り替えず、最終値の参照金利だけを複利・加重平均で組み立てる。
+    let as_of = AS_OF.parse()?;
+    let market = market(as_of, 0.0)?;
+    let config = config("")?;
+    let accrual = DateInterval::new("2026-01-05".parse()?, "2026-01-07".parse()?)?;
+    let tau = 1.0 / 365.0;
+    for (reference_rate, expected) in [
+        (
+            RateFutureRate::Compounded(accrual),
+            (RATE * 2.0 * tau).exp_m1() / (2.0 * tau),
+        ),
+        (
+            RateFutureRate::Averaged(accrual),
+            (RATE * tau).exp_m1() / tau,
+        ),
+    ] {
+        // 翌日物は fixture の既定どおり観測日の 16 時に公表するので、評価時点には何も公表されていない。
+        let reference = overnight_reference(reference_rate, 0)?;
+        let ctx = Context {
+            as_of,
+            refdata: &reference,
+            config: &config,
+            fixings: &FixingStore::default(),
+        };
+        let df = (-RATE * model_time(as_of, at(START, "16:00")?).years()).exp();
+        let pv = futures_pv(&future(Position::Long, 0.95)?, &ctx, &market)?.0;
+        assert_close!(
+            pv,
+            (1.0 - expected - 0.95) * df * MULTIPLIER,
+            PV_REL,
+            USD_ABS
+        );
+    }
+    Ok(())
+}
+
+#[test]
+fn overnight_futures_annualize_the_spanning_day_over_its_remaining_fraction() -> Result<()> {
+    // 評価時点 A が翌日物の観測の期間をまたぎ、その日の率が未公表のとき、[A, 終了) の forward の率を
+    // 残りの日数係数 τ_rem で年率化し、その率を元の 1 日の日数係数に当てる（参照データの章の
+    // fixing の規則）。成長を [A, 終了) に縮めたまま元の重みを掛けると、経過した分の利息が抜ける。
+    // 複利と加重平均のそれぞれを手計算と照合し、`futures_pv` で IR のリニア評価と解析の一致も確かめる。
+    let as_of = at(START, "20:00")?;
+    let market = market(as_of, 0.0)?;
+    let config = config("")?;
+    let accrual = DateInterval::new("2026-01-05".parse()?, "2026-01-07".parse()?)?;
+    let day = 1.0 / 365.0;
+    let first_end = at("2026-01-06", "16:00")?;
+    let second_end = at("2026-01-07", "16:00")?;
+    // 1 日目は A（20 時）から翌日の観測時刻（16 時）までの 20 時間だけが残る。
+    let remaining = 20.0 / 24.0 * day;
+    let first_growth = (RATE * model_time(as_of, first_end).years()).exp();
+    let second_growth = (RATE
+        * (model_time(as_of, second_end).years() - model_time(as_of, first_end).years()))
+    .exp();
+    let first_rate = (first_growth - 1.0) / remaining;
+    let second_rate = (second_growth - 1.0) / day;
+    for (reference_rate, expected) in [
+        (
+            RateFutureRate::Compounded(accrual),
+            ((1.0 + day * first_rate) * (1.0 + day * second_rate) - 1.0) / (2.0 * day),
+        ),
+        (
+            RateFutureRate::Averaged(accrual),
+            f64::midpoint(first_rate, second_rate),
+        ),
+    ] {
+        // 翌日物は翌営業日の 16 時に公表する。A では 1 日目が未公表で、期間が A をまたぐ。
+        let reference = overnight_reference(reference_rate, 1)?;
+        // 初日の清算価格は 16 時に公表済みなので、最新の決済価格 P_last として渡す。
+        let fixings = FixingStore::new(
+            &reference,
+            FixingStoreParts {
+                published: vec![(
+                    CLEARING.parse()?,
+                    ObservationPoint::Date(START.parse()?),
+                    0.954,
+                )],
+                ..FixingStoreParts::default()
+            },
+        )?;
+        let ctx = Context {
+            as_of,
+            refdata: &reference,
+            config: &config,
+            fixings: &fixings,
+        };
+        // 次の値洗いの受渡は翌日の 16 時。
+        let df = (-RATE * model_time(as_of, first_end).years()).exp();
+        let (pv, warnings) = futures_pv(&future(Position::Long, 0.95)?, &ctx, &market)?;
+        assert_close!(
+            pv,
+            (1.0 - expected - 0.954) * df * MULTIPLIER,
+            PV_REL,
+            USD_ABS
+        );
+        // 1 日目（01-05 の観測）は A で未公表なので、その率を予想した近似を記録する。2 日目は A の後に
+        // 始まる将来の観測なので記録しない。警告の一覧全体を比べ、他の近似が無いことも確かめる。
+        assert_eq!(
+            warnings,
+            [Warning::UnpublishedFixingForecast {
+                index: "TERM".parse()?,
+                point: ObservationPoint::Date(START.parse()?),
+            }]
+        );
+    }
+    Ok(())
+}
diff --git a/crates/rhodium-pricing/tests/rate_futures.rs b/crates/rhodium-pricing/tests/rate_futures.rs
new file mode 100644
index 0000000000000000000000000000000000000000..f45964366d7e2b056f25cc375e2eeb138a56f420
--- /dev/null
+++ b/crates/rhodium-pricing/tests/rate_futures.rs
@@ -0,0 +1,555 @@
+//! 金利先物の値洗いを、日次の値洗い口座の手計算、満期一括の契約との凸性の差、受渡の境界、
+//! 先物の quote での校正、値洗い口座の付利指標が定める予想の測度と AAD の感応度で照合する。
+
+pub mod support;
+
+use std::sync::Arc;
+
+use rhodium_core::{
+    Error, Result, Warning,
+    currency::{Asset, Currency},
+    day_count::DayCount,
+    fixing::{FixingStore, FixingStoreParts, ObservationPoint},
+    id::IndexId,
+    refdata::{IndexKind, RefData},
+    time::{Instant, model_time},
+};
+use rhodium_fixtures::{COLLATERAL, futures, prepare_standard};
+use rhodium_market::{
+    curve::{CurveRole, CurveSpec, Extrapolation},
+    fx::FxSpec,
+    market::Market,
+    quote::{QuoteKind, QuoteSet},
+    structure::{BlockId, MarketStructure, ParamId},
+};
+use rhodium_pricing::{
+    calibration::calibrate,
+    config::PricingConfig,
+    context::Context,
+    definition::{
+        CalibrationTarget, ElementSelector, InstrumentBinding, LgmDefinition, MarketDefinition,
+        QuoteTransform,
+    },
+    prepare::{PrepareInput, Prepared, evaluate, prepare},
+    resolve::CollateralKey,
+};
+use rhodium_products::{
+    ir::{Contract, Product},
+    standard::Position,
+};
+use rhodium_testkit::{PV_REL, SEQUENTIAL, assert_close};
+use scioxide_ad::Tape;
+use scioxide_core::{InvalidInputReason, SciOxideError, argument};
+
+use futures::{
+    AS_OF, CLEARING, FINAL, MULTIPLIER, PILLAR, RATE, START, at, config, curve, futures_definition,
+    futures_template, market, refdata, reference_parts,
+};
+use support::{
+    USD_ABS, future, futures_pv, fx_config, fx_market, fx_reference, mtm_swap,
+    rate_future_convention,
+};
+
+/// 清算の時点 `sessions` ごとに値洗いする 2026 年 1 月限の先物価格を、平均回帰 0、ボラ `vol` の
+/// LGM で実装と独立に手計算する。`sessions` を最終清算の一点にすると、満期に一度だけ払う契約の価格になる。
+fn oracle_price(as_of: Instant, sessions: &[Instant], vol: f64) -> Result<f64> {
+    let ts: Vec<_> = std::iter::once(0.0)
+        .chain(sessions.iter().map(|&date| model_time(as_of, date).years()))
+        .collect();
+    let t1 = model_time(as_of, at(FINAL, "16:00")?).years();
+    let t2 = model_time(as_of, at("2026-01-14", "16:00")?).years();
+    let terminal = ts[ts.len() - 1];
+    // a=0 の H は t。共分散 σ² min(s,t) に各清算時点の係数を直接掛ける。
+    let covariance = ts
+        .windows(2)
+        .map(|w| (w[1] - w[0]) * w[0].min(t1))
+        .sum::<f64>()
+        - terminal * terminal.min(t1);
+    let tau = 7.0 / 365.0;
+    let log_mean = RATE * (t2 - t1) + vol * vol * (t2 - t1) * (t2 * t1 + covariance);
+    Ok(1.0 - log_mean.exp_m1() / tau)
+}
+
+#[test]
+fn futures_margin_ir_vs_analytic() -> Result<()> {
+    // 日次に共通の価格を 実装と別の式で手計算し、将来の二回目以後の値洗いが消えることを確かめる。
+    let as_of = AS_OF.parse()?;
+    let reference = refdata("0D".parse()?, "16:00:00")?;
+    let config = config("")?;
+    let fixings = FixingStore::default();
+    let context = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &fixings,
+    };
+    let margins = [
+        at(START, "16:00")?,
+        at("2026-01-06", "16:00")?,
+        at(FINAL, "16:00")?,
+    ];
+    let df = (-RATE * model_time(as_of, margins[0]).years()).exp();
+    for vol in [0.0, 0.5] {
+        let market = market(as_of, vol)?;
+        for (side, sign) in [(Position::Long, 1.0), (Position::Short, -1.0)] {
+            let (actual, warnings) = futures_pv(&future(side, 0.95)?, &context, &market)?;
+            let expected = sign * MULTIPLIER * df * (oracle_price(as_of, &margins, vol)? - 0.95);
+            assert_close!(actual, expected, PV_REL, USD_ABS);
+            assert_eq!(warnings, []);
+        }
+    }
+    Ok(())
+}
+
+#[test]
+fn futures_settlement_delay() -> Result<()> {
+    // 同日内の遅れも翌日の受渡も許可制の近似になり、禁止した設定では評価を拒否する。
+    let as_of = AS_OF.parse()?;
+    let market = market(as_of, 0.0)?;
+    let config_ok = config("")?;
+    let config_disabled = config("[conventions]\ndisabled = ['payment_delay_ignored']")?;
+    for (lag, time) in [("0D", "18:00:00"), ("1BD", "08:00:00")] {
+        let reference = refdata(lag.parse()?, time)?;
+        let fixings = FixingStore::default();
+        let context = Context {
+            as_of,
+            refdata: &reference,
+            config: &config_ok,
+            fixings: &fixings,
+        };
+        let product = future(Position::Long, 0.95)?;
+        assert_eq!(
+            futures_pv(&product, &context, &market)?.1,
+            vec![Warning::PaymentDelayIgnored]
+        );
+        let disabled = Context {
+            config: &config_disabled,
+            ..context
+        };
+        assert!(matches!(
+            prepare_standard(&product, START.parse()?, &COLLATERAL, &market, &disabled),
+            Err(Error::Unsupported(_))
+        ));
+    }
+    Ok(())
+}
+
+#[test]
+fn futures_vs_forward_settlement() -> Result<()> {
+    // 実装の先物価格が日次の値洗い口座 Q_M の式に一致し、最終日に一度だけ払う契約の
+    // 価格（最終受渡の割引債を数値基準とする式）とは凸性の分だけ異なることを確かめる。
+    // oracle 同士の比較だけでは、実装が満期一括の数値基準へ退行しても検出できないため、
+    // ライブラリの評価から先物価格を逆算して照合する。
+    let as_of = AS_OF.parse()?;
+    let reference = refdata("0D".parse()?, "16:00:00")?;
+    let config = config("")?;
+    let fixings = FixingStore::default();
+    let context = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &fixings,
+    };
+    let vol = 0.5;
+    let market = market(as_of, vol)?;
+    let first = at(START, "16:00")?;
+    let final_at = at(FINAL, "16:00")?;
+    let daily = [first, at("2026-01-06", "16:00")?, final_at];
+    // PV は次の値洗い M·P(s_1)·(F − K) だけなので、そこから実装の F を求める。
+    let trade_price = 0.95;
+    let df = (-RATE * model_time(as_of, first).years()).exp();
+    let (pv, _) = futures_pv(&future(Position::Long, trade_price)?, &context, &market)?;
+    let implied = pv / (MULTIPLIER * df) + trade_price;
+    assert_close!(
+        implied,
+        oracle_price(as_of, &daily, vol)?,
+        1e-12_f64,
+        1e-8_f64
+    );
+    let forward = oracle_price(as_of, &[final_at], vol)?;
+    assert!((implied - forward).abs() > 1e-6);
+    Ok(())
+}
+
+#[test]
+fn futures_publication_and_payment_boundaries_use_realized_margin() -> Result<()> {
+    // 確定済みで未払いの受渡は取引価格との差額。決済時点を除くと次回の差額だけ残る。
+    let as_of = at(START, "16:00")?;
+    let reference = refdata("0D".parse()?, "16:00:00")?;
+    let config = config("")?;
+    let product = future(Position::Long, 0.95)?;
+    let market = market(as_of, 0.0)?;
+    let missing = FixingStore::default();
+    let ctx = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &missing,
+    };
+    assert!(matches!(
+        prepare_standard(&product, START.parse()?, &COLLATERAL, &market, &ctx),
+        Err(Error::MissingData(_))
+    ));
+    let fixings = FixingStore::new(
+        &reference,
+        FixingStoreParts {
+            published: vec![(
+                CLEARING.parse()?,
+                ObservationPoint::Date(START.parse()?),
+                0.954,
+            )],
+            ..FixingStoreParts::default()
+        },
+    )?;
+    let context = Context {
+        fixings: &fixings,
+        ..ctx
+    };
+    let (excluded, _) = futures_pv(&product, &context, &market)?;
+
+    let include: PricingConfig = toml::from_str(
+        r#"
+version = 1
+reporting_currency = "USD"
+settlement = "include"
+[[discounting]]
+currency = "USD"
+collateral = "uncollateralized"
+curve = "USD"
+[forecasting]
+TERM = "USD"
+"#,
+    )
+    .map_err(|error| Error::InvalidData(error.to_string()))?;
+    let included = futures_pv(
+        &product,
+        &Context {
+            config: &include,
+            ..context
+        },
+        &market,
+    )?
+    .0;
+    assert_close!(included - excluded, 4.0_f64, PV_REL, USD_ABS);
+    Ok(())
+}
+
+#[test]
+fn futures_quote_calibrates_the_trade_to_par() -> Result<()> {
+    // 先物のマークと独立した参照価格で校正し、同じ限月を約定した取引の PV がゼロであることを確かめる。
+    let as_of = AS_OF.parse()?;
+    let reference = refdata("0D".parse()?, "16:00:00")?;
+    let config = config("")?;
+    let context = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &FixingStore::default(),
+    };
+    let template = futures_template()?;
+    let definition = futures_definition()?;
+    let quote = oracle_price(
+        as_of,
+        &[
+            at(START, "16:00")?,
+            at("2026-01-06", "16:00")?,
+            at(FINAL, "16:00")?,
+        ],
+        0.5,
+    )?;
+    let market = calibrate(
+        &definition,
+        &QuoteSet::new([
+            ("USD.FUT.2026-01".parse()?, quote),
+            ("USD.LGM.MR".parse()?, 0.0),
+            ("USD.LGM.VOL".parse()?, 0.5),
+        ])?,
+        &context,
+    )?;
+    let product = template.build_for_valuation(quote, START.parse()?, as_of, &reference)?;
+    assert_close!(
+        futures_pv(&product, &context, market.market())?.0,
+        0.0_f64,
+        PV_REL,
+        USD_ABS
+    );
+    Ok(())
+}
+
+fn check_gradient(
+    plan: &Prepared,
+    market: &Market<f64>,
+    block: &BlockId,
+    element: usize,
+    bump: f64,
+) -> Result<f64> {
+    let structure = market.structure();
+    let mut gradient = vec![0.0; market.theta().len()];
+    // 準備済みの商品をそのまま再利用し、AAD と中心差分の双方で同じマークを動かす。
+    Tape::default()
+        .gradient_into(
+            |theta| {
+                Market::from_theta(structure.clone(), theta.to_vec())
+                    .and_then(|market| evaluate(plan, &market))
+                    .map(|pv| pv.get(Currency::USD))
+                    .map_err(|error| {
+                        SciOxideError::invalid_input(
+                            argument!("market"),
+                            InvalidInputReason::Inconsistent,
+                            error.to_string(),
+                        )
+                    })
+            },
+            market.theta(),
+            &mut gradient,
+            SEQUENTIAL,
+        )
+        .map_err(|error| Error::Numerical(error.to_string()))?;
+    let position = structure
+        .layout()
+        .range(block)
+        .ok_or_else(|| Error::Internal("missing AAD block".into()))?
+        .start
+        + element;
+    let evaluate_at = |delta: f64| -> Result<f64> {
+        let mut theta = market.theta().to_vec();
+        theta[position] += delta;
+        let market = Market::from_theta(structure.clone(), theta)?;
+        Ok(evaluate(plan, &market)?.get(Currency::USD))
+    };
+    let central = (evaluate_at(bump)? - evaluate_at(-bump)?) / (2.0 * bump);
+    assert_close!(gradient[position], central, 1e-7_f64, 5e-8_f64);
+    Ok(gradient[position])
+}
+
+#[test]
+fn futures_convexity_and_mtm_fx_marks_have_aad_gradients() -> Result<()> {
+    // LGM ボラと評価時点の FX マークをそれぞれ動かし、値洗いと元本リセットへ感応度が流れることを確かめる。
+    let as_of = AS_OF.parse()?;
+    let fixings = FixingStore::default();
+    let futures_market = market(as_of, 0.5)?;
+    let futures_reference = refdata("0D".parse()?, "16:00:00")?;
+    let futures_config = config("")?;
+    let futures_context = Context {
+        as_of,
+        refdata: &futures_reference,
+        config: &futures_config,
+        fixings: &fixings,
+    };
+    let product = future(Position::Long, 0.95)?;
+    let plan = prepare_standard(
+        &product,
+        START.parse()?,
+        &COLLATERAL,
+        &futures_market,
+        &futures_context,
+    )?;
+    assert!(
+        check_gradient(
+            &plan,
+            &futures_market,
+            &MarketStructure::lgm_block(Currency::USD)?,
+            1,
+            1e-3
+        )?
+        .abs()
+            > 1e-4
+    );
+
+    let fx = fx_market(as_of)?;
+    let reference = fx_reference("USDJPY".parse()?)?;
+    let fx_config = fx_config("")?;
+    let fx_context = Context {
+        as_of,
+        refdata: &reference,
+        config: &fx_config,
+        fixings: &fixings,
+    };
+    let swap = mtm_swap()?;
+    let plan = prepare_standard(&swap, START.parse()?, &COLLATERAL, &fx, &fx_context)?;
+    // USD 脚にスプレッドを入れず同じ曲線で割り引くため、USD PV はほぼゼロ。支払いの連鎖を微分できることを照合する。
+    let _gradient = check_gradient(
+        &plan,
+        &fx,
+        &FxSpec::marks_block("USDJPY".parse()?)?,
+        0,
+        0.01,
+    )?;
+    Ok(())
+}
+
+#[test]
+fn futures_margin_index_cross_collateral_projection_is_recorded() -> Result<()> {
+    // 値洗いの口座の付利指標（margin_index）の予想が値洗いの割引と等しい率なら、値洗いの列の評価は
+    // 厳密（Q_M の測度と割引が一致）。付利指標を別の曲線で予想すると、その差を無視する近似として
+    // CrossCollateralProjection を記録する。
+    let as_of = AS_OF.parse()?;
+    let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
+    let mut margin = parts.indices[&"TERM".parse()?].clone();
+    margin.id = "MARGIN".parse()?;
+    margin.kind = IndexKind::OvernightRate {
+        ccy: Currency::USD,
+        day_count: DayCount::Act360,
+    };
+    parts.indices.insert(margin.id.clone(), margin);
+    rate_future_convention(&mut parts)?.margin_index = Some("MARGIN".parse()?);
+    let reference = RefData::new(parts)?;
+    let alternative = CurveSpec::new(
+        "USD-ALT".parse()?,
+        Currency::USD,
+        [CurveRole::Projection],
+        [PILLAR.parse()?],
+        Extrapolation::Reject,
+    )?;
+    let structure = Arc::new(MarketStructure::with_fx_and_lgm(
+        as_of,
+        [curve()?, alternative],
+        [],
+        [(Currency::USD, vec![])],
+    )?);
+    let mut theta = vec![0.0; structure.layout().len()];
+    let t = model_time(as_of, PILLAR.parse()?).years();
+    for curve in ["USD", "USD-ALT"] {
+        let range = structure
+            .layout()
+            .range(&BlockId::nodes(curve.parse()?))
+            .ok_or_else(|| Error::Internal("missing curve block".into()))?;
+        theta[range.start] = -RATE * t;
+    }
+    let market = Market::from_theta(structure, theta)?;
+    let fixings = FixingStore::default();
+    let warnings = |margin_curve: &str| -> Result<Vec<Warning>> {
+        let config = config(&format!("MARGIN = \"{margin_curve}\""))?;
+        let ctx = Context {
+            as_of,
+            refdata: &reference,
+            config: &config,
+            fixings: &fixings,
+        };
+        Ok(futures_pv(&future(Position::Long, 0.95)?, &ctx, &market)?.1)
+    };
+    // 警告の一覧全体を比べ、付利指標の予想の曲線を変えたときに担保の違いの近似だけが増えることを確かめる。
+    assert_eq!(warnings("USD")?, []);
+    assert_eq!(warnings("USD-ALT")?, [Warning::CrossCollateralProjection]);
+    Ok(())
+}
+
+/// USD の予想曲線を先物だけで、別の割引曲線 USD-OIS を直接マークで当てる市場の定義。
+fn futures_with_discount_definition() -> Result<MarketDefinition> {
+    let discount = CurveSpec::new(
+        "USD-OIS".parse()?,
+        Currency::USD,
+        [CurveRole::Discount],
+        [PILLAR.parse()?],
+        Extrapolation::Reject,
+    )?;
+    let mark = InstrumentBinding::new(
+        "USD-OIS.DF".parse()?,
+        QuoteKind::Mark,
+        CalibrationTarget::DirectMark,
+        QuoteTransform::Identity,
+        CollateralKey::Uncollateralized,
+        ElementSelector::Param(ParamId::new(BlockId::nodes("USD-OIS".parse()?), 0)),
+    )?;
+    let future = InstrumentBinding::new(
+        "USD.FUT.2026-01".parse()?,
+        QuoteKind::FuturesPrice,
+        CalibrationTarget::Price(futures_template()?),
+        QuoteTransform::Identity,
+        CollateralKey::Uncollateralized,
+        ElementSelector::Block(BlockId::nodes("USD".parse()?)),
+    )?;
+    MarketDefinition::new(
+        "USD".parse()?,
+        "Etc/UTC".parse()?,
+        [curve()?, discount],
+        [future, mark],
+    )?
+    .with_lgm([(
+        Currency::USD,
+        LgmDefinition {
+            breaks: vec![],
+            mean_reversion: vec!["USD.LGM.MR".parse()?],
+            volatility: vec!["USD.LGM.VOL".parse()?],
+        },
+    )])
+}
+
+#[test]
+fn futures_identify_the_forecast_measure_only_through_their_margin_index() -> Result<()> {
+    // 先物の価格は値洗いの口座の測度で決まり、取引の担保の割引因子は打ち消し合うので、先物で当てた
+    // 予想曲線の予想の担保は、binding の担保からは分からない。USD の予想曲線を先物だけで当て、
+    // 無担保の行で別の割引曲線 USD-OIS に割り引くとき:
+    // - 付利指標の無い先物は根拠にならず測度は不明なので、ターム金利のクーポンは予想と割引の担保の
+    //   違いを決定論的なベーシスの仮定で扱い、CrossCollateralProjection を記録する。
+    // - 付利指標 MARGIN を USD-OIS で予想する先物は、口座の測度が USD-OIS の割引だと示すので、
+    //   同じ割引のクーポンは担保の調整を要しない。
+    let as_of = AS_OF.parse()?;
+    let definition = futures_with_discount_definition()?;
+    let config: PricingConfig = toml::from_str(
+        r#"
+version = 1
+reporting_currency = "USD"
+[[discounting]]
+currency = "USD"
+collateral = "uncollateralized"
+curve = "USD-OIS"
+[forecasting]
+TERM = "USD"
+MARGIN = "USD-OIS"
+"#,
+    )
+    .map_err(|error| Error::InvalidData(error.to_string()))?;
+    let quotes = QuoteSet::new([
+        ("USD.FUT.2026-01".parse()?, 0.96),
+        ("USD.LGM.MR".parse()?, 0.0),
+        ("USD.LGM.VOL".parse()?, 0.0),
+        (
+            "USD-OIS.DF".parse()?,
+            -RATE * model_time(as_of, PILLAR.parse()?).years(),
+        ),
+    ])?;
+    let fixings = FixingStore::default();
+    let warnings = |margin: bool| -> Result<Vec<Warning>> {
+        let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
+        let mut index = parts.indices[&"TERM".parse()?].clone();
+        index.id = "MARGIN".parse()?;
+        index.kind = IndexKind::OvernightRate {
+            ccy: Currency::USD,
+            day_count: DayCount::Act360,
+        };
+        parts.indices.insert(index.id.clone(), index);
+        if margin {
+            rate_future_convention(&mut parts)?.margin_index = Some("MARGIN".parse()?);
+        }
+        let reference = RefData::new(parts)?;
+        let ctx = Context {
+            as_of,
+            refdata: &reference,
+            config: &config,
+            fixings: &fixings,
+        };
+        let market = calibrate(&definition, &quotes, &ctx)?;
+        // 01-07 に観測する 1W のターム金利を期末の 01-14 に受け取るクーポン。
+        let term: IndexId = "TERM".parse()?;
+        let coupon = Product::new(
+            Contract::build(|builder| {
+                let rate = builder.fixing(&term, ObservationPoint::Date(FINAL.parse()?));
+                Ok(builder.pay(Asset::Cash(Currency::USD), rate, at("2026-01-14", "16:00")?))
+            })?,
+            &reference,
+        )?;
+        Ok(prepare(
+            PrepareInput::Ir(&coupon),
+            &COLLATERAL,
+            &[],
+            market.market().structure(),
+            &ctx,
+        )?
+        .warnings()
+        .to_vec())
+    };
+    assert_eq!(warnings(false)?, [Warning::CrossCollateralProjection]);
+    assert_eq!(warnings(true)?, [] as [Warning; 0]);
+    Ok(())
+}
diff --git a/crates/rhodium-pricing/tests/support/mod.rs b/crates/rhodium-pricing/tests/support/mod.rs
new file mode 100644
index 0000000000000000000000000000000000000000..68ec3bbff7e90c04749d950e1be423445a2b5a18
--- /dev/null
+++ b/crates/rhodium-pricing/tests/support/mod.rs
@@ -0,0 +1,446 @@
+//! pricing の統合テストが共有する補助。
+//!
+//! 各テストのファイルは `pub mod support;` で取り込む。補助を `pub` にしておけば、取り込んだファイルが
+//! 使わない補助で `dead_code` が出ない。複数の crate のテストが使う補助は `rhodium-fixtures` に、
+//! 数値の比較は `rhodium-testkit` に置く。
+//!
+//! ここの補助は [`rhodium_fixtures::futures`] の USD の市場と参照データを土台に、金利先物の取引と、
+//! JPY の曲線、USDJPY の FX、FX の fixing と XCCY の脚を足した市場と参照データを組み立てる。
+
+use std::{num::NonZeroU32, sync::Arc};
+
+use rhodium_core::{
+    Error, Result, Warning,
+    calendar::{
+        BaseCalendar, BusinessDayConvention, Extrapolation as CalendarExtrapolation, WeekendRule,
+    },
+    currency::Currency,
+    date::{Date, DateInterval, Weekday},
+    day_count::DayCount,
+    fx::CcyPair,
+    id::CalendarId,
+    money::{Money, MoneyBag},
+    quantity::{RateFuturesPrice, Spread},
+    refdata::{
+        Convention, FxConvention, FxForwardConvention, FxSpotLag, IndexDefinition, IndexKind,
+        RateFutureContract, RateFutureConvention, RefData, RefDataParts,
+    },
+    schedule::{
+        ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination, ScheduleSpec,
+    },
+    time::{Instant, model_time},
+};
+use rhodium_fixtures::{
+    COLLATERAL,
+    futures::{CLEARING, PILLAR, RATE, START, curve, reference_parts},
+    prepare_standard,
+};
+use rhodium_market::{
+    curve::{CurveRole, CurveSpec, Extrapolation},
+    fx::{FxAnchor, FxCollateral, FxSpec},
+    market::Market,
+    structure::{BlockId, MarketStructure},
+};
+use rhodium_pricing::{
+    config::{AnalyticPricer, Method, PricingConfig},
+    context::Context,
+    prepare::{PrepareInput, evaluate, prepare},
+};
+use rhodium_products::standard::{
+    Coupon, Leg, Notional, NotionalExchange, PayReceive, PaymentRounding, Position, RateFuture,
+    StandardProduct, Swap,
+};
+use rhodium_testkit::{PV_REL, assert_close, minor_unit_abs};
+
+/// USD の PV とキャッシュフローの絶対の下限 a。最小単位 0.01 ドルの 1e-6 倍。
+pub const USD_ABS: f64 = minor_unit_abs(2);
+/// JPY の PV とキャッシュフローの絶対の下限 a。最小単位 1 円の 1e-6 倍。
+pub const JPY_ABS: f64 = minor_unit_abs(0);
+
+/// 2026 年 1 月限の金利先物を、清算価格 `price` で 1 枚約定した取引。
+///
+/// # Errors
+///
+/// 清算指標の ID と先物価格の構築に失敗した場合。
+pub fn future(position: Position, price: f64) -> Result<StandardProduct> {
+    Ok(StandardProduct::RateFuture(RateFuture::new(
+        CLEARING.parse()?,
+        NonZeroU32::MIN,
+        RateFuturesPrice::new(price, DayCount::Act365F)?,
+        position,
+        None,
+    )))
+}
+
+/// 先物の取引を解析のプライサーと IR のリニア評価の両方で準備・評価し、USD の PV と警告を返す。
+///
+/// 二つの評価が同じ PV と同じ警告を出すことも確かめる。
+///
+/// # Errors
+///
+/// 取引の IR への変換、準備、評価のいずれかが失敗した場合。
+///
+/// # Panics
+///
+/// 二つの評価の方式、PV、警告のいずれかが期待と違う場合。
+pub fn futures_pv(
+    product: &StandardProduct,
+    context: &Context<'_>,
+    market: &Market<f64>,
+) -> Result<(f64, Vec<Warning>)> {
+    let trade_date = START.parse()?;
+    let analytic = prepare_standard(product, trade_date, &COLLATERAL, market, context)?;
+    let contract = product.to_product(trade_date, context.refdata)?;
+    let ir = prepare(
+        PrepareInput::Ir(&contract),
+        &COLLATERAL,
+        &[],
+        market.structure(),
+        context,
+    )?;
+    assert_eq!(
+        analytic.method(),
+        Method::Analytic(AnalyticPricer::RateFuture)
+    );
+    assert_eq!(ir.method(), Method::Linear);
+    let pv = evaluate(&analytic, market)?.get(Currency::USD);
+    assert_close!(
+        pv,
+        evaluate(&ir, market)?.get(Currency::USD),
+        PV_REL,
+        USD_ABS
+    );
+    assert_eq!(analytic.warnings(), ir.warnings());
+    Ok((pv, analytic.warnings().to_vec()))
+}
+
+/// [`fx_parts`] を解決した参照データ。
+///
+/// # Errors
+///
+/// 参照データの部品の組み立てか解決に失敗した場合。
+pub fn fx_reference(pair: CcyPair) -> Result<RefData> {
+    RefData::new(fx_parts(pair)?)
+}
+
+/// 先物の参照データに、JPY、USDJPY の慣行（T+2）とフォワードの慣行 `USDJPY-FWD`、JPY のターム指標
+/// `JPY-TERM` と、`pair` の向きで公表する FX の fixing `FX-FIX` を足した部品。
+///
+/// # Errors
+///
+/// 固定の ID の解析か、先物の参照データの部品の組み立てに失敗した場合。
+pub fn fx_parts(pair: CcyPair) -> Result<RefDataParts> {
+    let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
+    let clock = parts.currencies[&Currency::USD].settlement_time.clone();
+    parts.conventions.insert(
+        "USDJPY-FWD".parse()?,
+        Convention::FxForward(Box::new(FxForwardConvention {
+            pair: "USDJPY".parse()?,
+            bdc: BusinessDayConvention::Following,
+        })),
+    );
+    // スポットラグ、受渡の暦（両通貨の休日の和）と受渡の時刻はペアの慣行から読む。
+    let usd_jpy: CcyPair = "USDJPY".parse()?;
+    parts.fx_conventions.insert(
+        usd_jpy,
+        FxConvention {
+            pair: usd_jpy,
+            spot_lag: FxSpotLag::TwoBusinessDays,
+            settlement: clock,
+            pip_exponent: 2,
+        },
+    );
+    let mut currency = parts.currencies[&Currency::USD].clone();
+    currency.code = Currency::JPY;
+    currency.minor_units = 0;
+    parts.currencies.insert(currency.code, currency);
+    let mut term = index_definition(&mut parts, "TERM")?.clone();
+    term.id = "JPY-TERM".parse()?;
+    // TERM がターム金利でなくなったら、JPY の指標にせずに黙って進まないようエラーにする。
+    let IndexKind::TermRate { ccy, .. } = &mut term.kind else {
+        return Err(Error::InvalidData(
+            "fixture index TERM is not a term rate".into(),
+        ));
+    };
+    *ccy = Currency::JPY;
+    parts.indices.insert(term.id.clone(), term.clone());
+    term.id = "FX-FIX".parse()?;
+    term.kind = IndexKind::FxFixing { pair };
+    parts.indices.insert(term.id.clone(), term);
+    Ok(parts)
+}
+
+/// USDJPY の慣行のスポットラグを `spot_lag` にし、JPY の決済カレンダーに 2026-01-13（火）の
+/// 祝日を置いた参照データ。ペアの暦は USD と JPY の休日の和になる。
+///
+/// # Errors
+///
+/// 部品の組み立て、カレンダーの構築か参照データの解決に失敗した場合。
+pub fn fx_reference_with_lag(spot_lag: FxSpotLag) -> Result<RefData> {
+    let mut parts = fx_parts("USDJPY".parse()?)?;
+    let tokyo: CalendarId = "JPN".parse()?;
+    parts.calendars.insert(
+        tokyo.clone(),
+        BaseCalendar::new(
+            "2025-01-01".parse()?..="2026-12-31".parse()?,
+            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
+            vec!["2026-01-13".parse()?],
+            CalendarExtrapolation::Strict,
+        )?,
+    );
+    parts
+        .currencies
+        .get_mut(&Currency::JPY)
+        .ok_or_else(|| Error::MissingData("fixture has no JPY currency".into()))?
+        .settlement_calendar = tokyo;
+    usd_jpy_convention(&mut parts)?.spot_lag = spot_lag;
+    RefData::new(parts)
+}
+
+// 以下は、fixture の部品を書き換えるテストが書き換える先を借りる補助。`if let` で受けると、fixture の
+// ID や種類を変えたときに書き換えを黙って飛ばし、既定の場面のまま合格してしまう。見つからなければ
+// エラーを返し、テストを失敗させる。
+
+/// 部品の `id` の指標の定義。
+///
+/// # Errors
+///
+/// `id` が指標の ID として読めないか、部品にその指標が無い場合。
+pub fn index_definition<'a>(
+    parts: &'a mut RefDataParts,
+    id: &str,
+) -> Result<&'a mut IndexDefinition> {
+    parts
+        .indices
+        .get_mut(&id.parse()?)
+        .ok_or_else(|| Error::MissingData(format!("fixture has no index {id}")))
+}
+
+/// 部品の USDJPY の FX の慣行。
+///
+/// # Errors
+///
+/// 部品に USDJPY の慣行が無い場合。
+pub fn usd_jpy_convention(parts: &mut RefDataParts) -> Result<&mut FxConvention> {
+    parts
+        .fx_conventions
+        .get_mut(&"USDJPY".parse()?)
+        .ok_or_else(|| Error::MissingData("fixture has no USDJPY FX convention".into()))
+}
+
+/// 先物の参照データの部品の、金利先物の慣行 `USD-RATE-FUTURES`。
+///
+/// # Errors
+///
+/// 部品にその ID の慣行が無いか、金利先物の慣行でない場合。
+pub fn rate_future_convention(parts: &mut RefDataParts) -> Result<&mut RateFutureConvention> {
+    let Some(Convention::RateFuture(convention)) =
+        parts.conventions.get_mut(&"USD-RATE-FUTURES".parse()?)
+    else {
+        return Err(Error::MissingData(
+            "fixture has no USD-RATE-FUTURES rate future convention".into(),
+        ));
+    };
+    Ok(convention.as_mut())
+}
+
+/// [`rate_future_convention`] の 2026 年 1 月限。
+///
+/// # Errors
+///
+/// 部品に金利先物の慣行か、その 2026 年 1 月限が無い場合。
+pub fn january_contract(parts: &mut RefDataParts) -> Result<&mut RateFutureContract> {
+    rate_future_convention(parts)?
+        .contracts
+        .get_mut(&"2026-01".parse()?)
+        .ok_or_else(|| Error::MissingData("fixture has no 2026-01 rate future contract".into()))
+}
+
+/// ゼロレート [`RATE`] の USD 曲線、1% の JPY 曲線と、評価時点の USDJPY 150 のマークを持つ市場。
+/// FX の担保区分 `csa` は USD と JPY の曲線の組。
+///
+/// # Errors
+///
+/// 曲線や FX の仕様、市場の構造の構築に失敗した場合。
+pub fn fx_market(as_of: Instant) -> Result<Market<f64>> {
+    fx_market_with(as_of, [])
+}
+
+/// [`fx_market`] の曲線の後に、θ を持たない式の曲線 `extra` を足した市場。
+///
+/// # Errors
+///
+/// 曲線や FX の仕様、市場の構造の構築に失敗した場合。
+pub fn fx_market_with(
+    as_of: Instant,
+    extra: impl IntoIterator<Item = CurveSpec>,
+) -> Result<Market<f64>> {
+    let jpy = CurveSpec::new(
+        "JPY".parse()?,
+        Currency::JPY,
+        [CurveRole::Discount, CurveRole::Projection],
+        [PILLAR.parse()?],
+        Extrapolation::Reject,
+    )?;
+    let pair: CcyPair = "USDJPY".parse()?;
+    let structure = Arc::new(MarketStructure::with_fx(
+        as_of,
+        [curve()?, jpy].into_iter().chain(extra),
+        [FxSpec {
+            pair,
+            anchor: FxAnchor::Direct {
+                spot: as_of,
+                short_starts: vec![],
+                market_collateral: "csa".parse()?,
+            },
+            collaterals: vec![FxCollateral {
+                id: "csa".parse()?,
+                base: "USD".parse()?,
+                quote: "JPY".parse()?,
+            }],
+            via: None,
+        }],
+    )?);
+    let t = model_time(as_of, PILLAR.parse()?).years();
+    let mut theta = vec![0.0; structure.layout().len()];
+    for (id, value) in [
+        (BlockId::nodes("USD".parse()?), -RATE * t),
+        (BlockId::nodes("JPY".parse()?), -0.01 * t),
+        (FxSpec::marks_block(pair)?, 150.0),
+    ] {
+        let range = structure
+            .layout()
+            .range(&id)
+            .ok_or_else(|| Error::Internal("missing FX test block".into()))?;
+        theta[range.start] = value;
+    }
+    Market::from_theta(structure, theta)
+}
+
+/// USD と JPY の無担保の割引と、`TERM`、`JPY-TERM` の予想を [`fx_market`] の曲線で行う価格設定。
+/// `conventions` は TOML の末尾に足す。
+///
+/// # Errors
+///
+/// 価格設定の TOML が読めない場合。
+pub fn fx_config(conventions: &str) -> Result<PricingConfig> {
+    toml::from_str(&format!(
+        r#"
+version = 1
+reporting_currency = "USD"
+[[discounting]]
+currency = "USD"
+collateral = "uncollateralized"
+curve = "USD"
+[[discounting]]
+currency = "JPY"
+collateral = "uncollateralized"
+curve = "JPY"
+[forecasting]
+TERM = "USD"
+JPY-TERM = "JPY"
+{conventions}
+"#
+    ))
+    .map_err(|error| Error::InvalidData(error.to_string()))
+}
+
+/// 標準商品を先物の約定日に既定の担保で約定した取引として準備し、`market` で評価した通貨ごとの PV。
+///
+/// # Errors
+///
+/// 準備か評価に失敗した場合。
+pub fn price(
+    product: &StandardProduct,
+    market: &Market<f64>,
+    context: &Context<'_>,
+) -> Result<MoneyBag> {
+    let plan = prepare_standard(product, START.parse()?, &COLLATERAL, market, context)?;
+    evaluate(&plan, market)
+}
+
+/// `dates` の隣り合う日付を期間とし、各期間を `pays` の同じ位置の日に払う明示のスケジュール。
+///
+/// # Errors
+///
+/// 期間かスケジュールの構築に失敗した場合。
+pub fn schedule(dates: &[Date], pays: &[Date]) -> Result<ScheduleSpec> {
+    let periods = dates
+        .windows(2)
+        .zip(pays)
+        .map(|(dates, &payment)| {
+            Ok(ExplicitPeriod {
+                accrual: DateInterval::new(dates[0], dates[1])?,
+                payment,
+                reset: None,
+                settlement: None,
+                references: Vec::new(),
+            })
+        })
+        .collect::<Result<_>>()?;
+    Ok(ScheduleSpec::Explicit(ExplicitPeriods::new(
+        periods,
+        ExplicitDatePolicy::KeepAsConfirmed,
+        PeriodCombination::Sum,
+        None,
+    )?))
+}
+
+/// 2026-01-07 から 1 週ごとの 2 期間（01-07〜01-14、01-14〜01-21）を、各期末に払う明示のスケジュール。
+///
+/// # Errors
+///
+/// 日付の解析か、期間とスケジュールの構築に失敗した場合。
+pub fn explicit_schedule() -> Result<ScheduleSpec> {
+    let [start, middle, end] = ["2026-01-07", "2026-01-14", "2026-01-21"].map(str::parse::<Date>);
+    let [start, middle, end] = [start?, middle?, end?];
+    schedule(&[start, middle, end], &[middle, end])
+}
+
+/// 指標 `index` のターム金利にスプレッド `spread` を足す、ギアリング 1 のクーポン。
+///
+/// # Errors
+///
+/// 指標の ID の解析か、スプレッドの構築に失敗した場合。
+pub fn term(index: &str, spread: f64) -> Result<Coupon> {
+    Ok(Coupon::Term {
+        index: index.parse()?,
+        spread: Spread::new(spread)?,
+        gearing: 1.0,
+    })
+}
+
+/// [`explicit_schedule`] の `MtM` の XCCY。JPY の脚は元本 150,000 円を固定して `JPY-TERM` + 10bp を
+/// 受け、USD の脚は同じ基準元本を `FX-FIX` で期首ごとにリセットして `TERM` + 20bp を払う。
+/// 両脚とも期首と期末に元本を交換する。
+///
+/// # Errors
+///
+/// スケジュール、クーポン、脚かスワップの構築に失敗した場合。
+pub fn mtm_swap() -> Result<StandardProduct> {
+    let notional = Money::new(Currency::JPY, 150_000.0);
+    let schedule = explicit_schedule()?;
+    Ok(StandardProduct::Swap(Swap::from_legs([
+        Leg::new(
+            PayReceive::Receive,
+            notional,
+            schedule.clone(),
+            DayCount::Act365F,
+            term("JPY-TERM", 0.001)?,
+            NotionalExchange::InitialAndFinal,
+            PaymentRounding::None,
+        )?,
+        Leg::with_notional(
+            PayReceive::Pay,
+            Notional::FxReset {
+                base: notional,
+                fixing: "FX-FIX".parse()?,
+            },
+            schedule,
+            DayCount::Act365F,
+            term("TERM", 0.002)?,
+            NotionalExchange::InitialAndFinal,
+            PaymentRounding::None,
+        )?,
+    ])?))
+}
diff --git a/crates/rhodium-pricing/tests/volatility_calibration.rs b/crates/rhodium-pricing/tests/volatility_calibration.rs
index e31ab206c5bb852414d0b3f9067aa8a9791a693e..fa79e93382368e3598b8a35d576e4f707813de1f 100644
--- a/crates/rhodium-pricing/tests/volatility_calibration.rs
+++ b/crates/rhodium-pricing/tests/volatility_calibration.rs
@@ -1,6 +1,8 @@
 //! M2-1 と M2-2 の金利ボラ校正を独立の価格例と照合する。
 
 use rhodium_core::{Error, Result, Warning, quantity::VolKind, time::model_time};
+use rhodium_fixtures::{COLLATERAL, bumped};
+use rhodium_fixtures::{multicurve, ois};
 use rhodium_market::{
     curve::{CurveExpr, CurveRole, CurveSpec},
     quote::{QuoteId, QuoteSet},
@@ -16,24 +18,10 @@ use rhodium_pricing::{
             VolatilityDefinition, VolatilityMarket,
         },
     },
-    resolve::CollateralKey,
 };
 use rhodium_testkit::assert_close;
 use serde::Deserialize;
 
-#[allow(
-    dead_code,
-    reason = "共有 fixture の曲線監査用のヘルパーはここでは使わない"
-)]
-#[path = "fixtures/ois.rs"]
-mod ois;
-
-#[allow(dead_code, reason = "共有 fixture の慣行と参照データだけを使う")]
-#[path = "fixtures/multicurve.rs"]
-mod multicurve;
-
-const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
-
 #[derive(Deserialize)]
 struct Oracle {
     sabr: Vec<SabrCase>,
@@ -338,14 +326,6 @@ fn cube_groups(market: &CalibratedMarket) -> impl Iterator<Item = &CalibrationGr
     surface_groups(market, "JPY-CUBE")
 }
 
-fn bumped(quotes: &QuoteSet, id: &QuoteId, bump: f64) -> Result<QuoteSet> {
-    QuoteSet::new(
-        quotes
-            .iter()
-            .map(|(q, value)| (q.clone(), if q == id { value + bump } else { value })),
-    )
-}
-
 #[test]
 fn swaption_cube_reprices_market_quotes() -> Result<()> {
     // 先に ATM-only で各ノードの曲線由来の F を読む。次に既知の SABR で各ウィングを作る。
diff --git a/crates/rhodium-pricing/tests/xccy.rs b/crates/rhodium-pricing/tests/xccy.rs
new file mode 100644
index 0000000000000000000000000000000000000000..c266d6e6c06b1067bb9996c70774e02347268be9
--- /dev/null
+++ b/crates/rhodium-pricing/tests/xccy.rs
@@ -0,0 +1,346 @@
+//! 元本を FX でリセットする `MtM` の XCCY と元本を固定した XCCY を、両通貨の受渡の手計算で照合する。
+//! 元本のリセットの観測日はペアのスポットラグで遡り、その fixing はペアのスポットの受渡の FX で予想する。
+
+pub mod support;
+
+use std::collections::HashSet;
+
+use rhodium_core::{
+    Error, Result, Warning,
+    currency::{Asset, Currency},
+    date::{Date, DateInterval},
+    day_count::DayCount,
+    fixing::{FixingStore, ObservationPoint},
+    id::IndexId,
+    money::Money,
+    quantity::{Rate, Simple},
+    refdata::{FxSpotLag, RefData},
+    schedule::{
+        ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination, ScheduleSpec,
+    },
+    time::model_time,
+};
+use rhodium_fixtures::{COLLATERAL, futures, prepare_standard};
+use rhodium_market::fx::FxSpec;
+use rhodium_pricing::{
+    context::Context,
+    prepare::{PrepareInput, evaluate, prepare},
+};
+use rhodium_products::{
+    ir::{
+        Contract, Product,
+        analysis::{Observation, observations},
+    },
+    standard::{
+        Coupon, Leg, Notional, NotionalExchange, PayReceive, PaymentRounding, StandardProduct, Swap,
+    },
+};
+use rhodium_testkit::{PV_REL, assert_close};
+
+use futures::{AS_OF, RATE, START, at};
+use support::{
+    JPY_ABS, USD_ABS, explicit_schedule, fx_config, fx_market, fx_parts, fx_reference,
+    fx_reference_with_lag, mtm_swap, price, term, usd_jpy_convention,
+};
+
+#[test]
+fn mtm_xccy_reset_fixing_dates_follow_the_pair_spot_lag() -> Result<()> {
+    // 元本のリセットは期首がスポット受渡になる日に観測する。期首は 2026-01-07（水）と 01-14（水）で、
+    // ペアの暦は USD と JPY の休日の和（JPY は 01-13 が祝日）。
+    // T+1: 01-07 の 1 営業日前は 01-06、01-14 の 1 営業日前は 01-13 を飛ばして 01-12。
+    // T+2: 01-07 の 2 営業日前は 01-05、01-14 の 2 営業日前は 01-12、01-09。
+    // 両脚のターム金利のクーポン（TERM と JPY-TERM、値の開始のずれ 0D）は期首の 01-07 と 01-14 に
+    // 観測し、スポットラグでは動かない。観測の一覧全体を比べ、FX の観測がリセットの二つだけで、
+    // 別の日の FX や他の種類の観測が無いことも確かめる。
+    let swap = mtm_swap()?;
+    let [term, jpy_term, fx] = ["TERM", "JPY-TERM", "FX-FIX"].map(str::parse::<IndexId>);
+    let [term, jpy_term, fx] = [term?, jpy_term?, fx?];
+    for (spot_lag, resets) in [
+        (FxSpotLag::OneBusinessDay, ["2026-01-06", "2026-01-12"]),
+        (FxSpotLag::TwoBusinessDays, ["2026-01-05", "2026-01-09"]),
+    ] {
+        let reference = fx_reference_with_lag(spot_lag)?;
+        let product = swap.to_product(START.parse()?, &reference)?;
+        let actual: HashSet<_> = observations(&product).into_iter().collect();
+        let expected = [
+            (&term, "2026-01-07"),
+            (&term, "2026-01-14"),
+            (&jpy_term, "2026-01-07"),
+            (&jpy_term, "2026-01-14"),
+            (&fx, resets[0]),
+            (&fx, resets[1]),
+        ]
+        .into_iter()
+        .map(|(index, date)| {
+            Ok(Observation::Fixing {
+                index,
+                point: ObservationPoint::Date(date.parse()?),
+            })
+        })
+        .collect::<Result<HashSet<_>>>()?;
+        assert_eq!(actual, expected, "{spot_lag:?}");
+    }
+    Ok(())
+}
+
+#[test]
+fn mtm_xccy_fx_reset_matches_cashflows_in_both_directions() -> Result<()> {
+    // 同じ基準元本を正向きの fixing では割り、逆向きでは掛ける。期中の交換は元本の差額だけ。
+    let as_of = AS_OF.parse()?;
+    let market = fx_market(as_of)?;
+    let config = fx_config("")?;
+    let disabled = fx_config("[conventions]\ndisabled = ['mtm_fx_forward_product']")?;
+    let swap = mtm_swap()?;
+    let fixings = FixingStore::default();
+    let pay = ["2026-01-07", "2026-01-14", "2026-01-21"].map(|date| at(date, "16:00"));
+    let [start, middle, end] = pay;
+    let [start, middle, end] = [start?, middle?, end?];
+    let discount = |rate: f64, at| (-rate * model_time(as_of, at).years()).exp();
+    let loan = |at| 1000.0 * discount(0.01, at) / discount(RATE, at);
+    let [n0, n1] = [loan(start), loan(middle)];
+    let dollar_coupon = (RATE * (7.0 / 365.0)).exp_m1() + 0.002 * (7.0 / 365.0);
+    let yen_coupon = (0.01_f64 * (7.0 / 365.0)).exp_m1() + 0.001 * (7.0 / 365.0);
+    let usd = n0 * discount(RATE, start) + (n1 - n0 - n0 * dollar_coupon) * discount(RATE, middle)
+        - n1 * (1.0 + dollar_coupon) * discount(RATE, end);
+    let jpy = 150_000.0
+        * (-discount(0.01, start)
+            + yen_coupon * discount(0.01, middle)
+            + (1.0 + yen_coupon) * discount(0.01, end));
+
+    for pair in ["USDJPY", "JPYUSD"] {
+        let reference = fx_reference(pair.parse()?)?;
+        let ctx = Context {
+            as_of,
+            refdata: &reference,
+            config: &config,
+            fixings: &fixings,
+        };
+        let plan = prepare_standard(&swap, START.parse()?, &COLLATERAL, &market, &ctx)?;
+        assert!(
+            plan.blocks()
+                .contains(&FxSpec::marks_block("USDJPY".parse()?)?)
+        );
+        assert_eq!(plan.warnings(), &[Warning::MtmFxForwardProduct]);
+        let actual = evaluate(&plan, &market)?;
+        assert_close!(actual.get(Currency::USD), usd, PV_REL, USD_ABS);
+        assert_close!(actual.get(Currency::JPY), jpy, PV_REL, JPY_ABS);
+        assert!(matches!(
+            price(
+                &swap,
+                &market,
+                &Context {
+                    config: &disabled,
+                    ..ctx
+                }
+            ),
+            Err(Error::Unsupported(_))
+        ));
+    }
+    Ok(())
+}
+
+/// 固定クーポンの `MtM` の XCCY。JPY の脚は元本 150,000 円を固定し、USD の脚の元本を `FX-FIX` で
+/// リセットする。`reset` は 2 期目のリセットの明示日（無ければ規則どおり期首から遡る）。
+fn fixed_mtm_swap(dates: &[&str; 3], reset: Option<&str>) -> Result<StandardProduct> {
+    let schedule = |reset: Option<&str>| -> Result<ScheduleSpec> {
+        let periods = dates
+            .windows(2)
+            .enumerate()
+            .map(|(i, pair)| {
+                let (start, end) = (pair[0].parse::<Date>()?, pair[1].parse::<Date>()?);
+                Ok(ExplicitPeriod {
+                    accrual: DateInterval::new(start, end)?,
+                    payment: end,
+                    reset: if i == 1 {
+                        reset.map(str::parse).transpose()?
+                    } else {
+                        None
+                    },
+                    settlement: None,
+                    references: vec![],
+                })
+            })
+            .collect::<Result<_>>()?;
+        Ok(ScheduleSpec::Explicit(ExplicitPeriods::new(
+            periods,
+            ExplicitDatePolicy::KeepAsConfirmed,
+            PeriodCombination::Sum,
+            None,
+        )?))
+    };
+    let fixed = || -> Result<Coupon> {
+        Ok(Coupon::Fixed {
+            rate: Rate::new(0.03, DayCount::Act365F, Simple)?,
+        })
+    };
+    let notional = Money::new(Currency::JPY, 150_000.0);
+    Ok(StandardProduct::Swap(Swap::from_legs([
+        Leg::new(
+            PayReceive::Receive,
+            notional,
+            schedule(None)?,
+            DayCount::Act365F,
+            fixed()?,
+            NotionalExchange::InitialAndFinal,
+            PaymentRounding::None,
+        )?,
+        Leg::with_notional(
+            PayReceive::Pay,
+            Notional::FxReset {
+                base: notional,
+                fixing: "FX-FIX".parse()?,
+            },
+            schedule(reset)?,
+            DayCount::Act365F,
+            fixed()?,
+            NotionalExchange::InitialAndFinal,
+            PaymentRounding::None,
+        )?,
+    ])?))
+}
+
+#[test]
+fn mtm_xccy_reset_reads_the_fixing_at_the_pair_spot_delivery() -> Result<()> {
+    // 同じ指標・観測日の FX fixing の公表値は一つなので、元本リセットも通常の fixing と同じく
+    // ペアのスポットの受渡（ここでは UTC 06 時。東京 15 時に当たる）の F で予想する。
+    // 新しい元本を交換する瞬間（USD の決済時刻 UTC 16 時）とは別に扱う。
+    // 明示の reset 日 01-08 の fixing はスポットの受渡が 01-12 で、元本の交換 01-14 より早いので
+    // 支払いの遅れの近似を記録する。クーポンは固定にして、FX の予想だけを照合する。
+    let as_of = AS_OF.parse()?;
+    let market = fx_market(as_of)?;
+    let config = fx_config("")?;
+    let fixings = FixingStore::default();
+    let mut parts = fx_parts("USDJPY".parse()?)?;
+    usd_jpy_convention(&mut parts)?.settlement.time = "06:00:00".parse()?;
+    let reference = RefData::new(parts)?;
+    let ctx = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &fixings,
+    };
+    let discount = |rate: f64, at| (-rate * model_time(as_of, at).years()).exp();
+    // F(v) = 150·exp(−(r_USD − r_JPY)·t_v)。リセット後の USD の元本は 150,000 円 / F(v)。
+    let forward = |at| 150.0 * discount(RATE - 0.01, at);
+    let fx: IndexId = "FX-FIX".parse()?;
+
+    // 通常の fixing: 01-05 の fixing の受渡は T+2 の 01-07 06 時。
+    let value = at("2026-01-07", "06:00")?;
+    let pay = at("2026-01-07", "16:00")?;
+    let plain = Product::new(
+        Contract::build(|builder| {
+            let rate = builder.fixing(&fx, ObservationPoint::Date("2026-01-05".parse()?));
+            Ok(builder.pay(Asset::Cash(Currency::JPY), rate, pay))
+        })?,
+        &reference,
+    )?;
+    let plan = prepare(
+        PrepareInput::Ir(&plain),
+        &COLLATERAL,
+        &[],
+        market.structure(),
+        &ctx,
+    )?;
+    assert_close!(
+        evaluate(&plan, &market)?.get(Currency::JPY),
+        forward(value) * discount(0.01, pay),
+        PV_REL,
+        JPY_ABS
+    );
+
+    let dates = ["2026-01-07", "2026-01-14", "2026-01-21"];
+    let [start, middle, end] = dates.map(|date| at(date, "16:00"));
+    let [start, middle, end] = [start?, middle?, end?];
+    let coupon = 0.03 * 7.0 / 365.0;
+    // 支払う USD の脚: 期首に元本 n0 を受け、期中に元本の差額 n1 − n0 を受けてクーポンを払い、
+    // 期末に n1 とクーポンを払う。
+    let usd = |n0: f64, n1: f64| {
+        n0 * discount(RATE, start) + (n1 - n0 - n0 * coupon) * discount(RATE, middle)
+            - n1 * (1.0 + coupon) * discount(RATE, end)
+    };
+    for (reset, second_value, warnings) in [
+        // 規則どおりの reset はスポットの受渡が期首の 01-14 と同じ日なので近似は無い。
+        (None, "2026-01-14", vec![]),
+        // 明示の reset 01-08（木）の T+2 の受渡は 01-12（月）。
+        (
+            Some("2026-01-08"),
+            "2026-01-12",
+            vec![Warning::PaymentDelayIgnored],
+        ),
+    ] {
+        let product = fixed_mtm_swap(&dates, reset)?;
+        let plan = prepare_standard(&product, START.parse()?, &COLLATERAL, &market, &ctx)?;
+        assert_eq!(plan.warnings(), warnings.as_slice(), "{reset:?}");
+        let n0 = 150_000.0 / forward(value);
+        let n1 = 150_000.0 / forward(at(second_value, "06:00")?);
+        assert_close!(
+            evaluate(&plan, &market)?.get(Currency::USD),
+            usd(n0, n1),
+            PV_REL,
+            USD_ABS
+        );
+    }
+    Ok(())
+}
+
+#[test]
+fn fixed_notional_xccy_exchanges_both_currencies() -> Result<()> {
+    // 元本を固定した XCCY は、各通貨の脚が期首に元本を払い、期末に元本を受け取り、クーポンを
+    // それぞれの通貨の曲線で割り引く（FX のリセットは無い）。JPY の脚は 1% の固定、USD の脚は
+    // ターム金利の 1W（RATE の連続金利から求めた単利）にスプレッド 20bp。
+    let as_of = AS_OF.parse()?;
+    let market = fx_market(as_of)?;
+    let config = fx_config("")?;
+    let reference = fx_reference("USDJPY".parse()?)?;
+    let ctx = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &FixingStore::default(),
+    };
+    let schedule = explicit_schedule()?;
+    let swap = StandardProduct::Swap(Swap::from_legs([
+        Leg::new(
+            PayReceive::Receive,
+            Money::new(Currency::JPY, 150_000.0),
+            schedule.clone(),
+            DayCount::Act365F,
+            Coupon::Fixed {
+                rate: Rate::new(0.01, DayCount::Act365F, Simple)?,
+            },
+            NotionalExchange::InitialAndFinal,
+            PaymentRounding::None,
+        )?,
+        Leg::new(
+            PayReceive::Pay,
+            Money::new(Currency::USD, 1_000.0),
+            schedule,
+            DayCount::Act365F,
+            term("TERM", 0.002)?,
+            NotionalExchange::InitialAndFinal,
+            PaymentRounding::None,
+        )?,
+    ])?);
+    let plan = prepare_standard(&swap, START.parse()?, &COLLATERAL, &market, &ctx)?;
+    // 元本が固定なので、FX と金利の積の近似（MtmFxForwardProduct）は使わない。クーポンは期末に払い、
+    // 元本は各通貨の決済時刻に交換するので、他の近似も記録しない。
+    assert_eq!(plan.warnings(), []);
+    let pv = evaluate(&plan, &market)?;
+    let discount = |rate: f64, date: &str| -> Result<f64> {
+        Ok((-rate * model_time(as_of, at(date, "16:00")?).years()).exp())
+    };
+    let tau = 7.0 / 365.0;
+    let [start, middle, end] = ["2026-01-07", "2026-01-14", "2026-01-21"];
+    let jpy = 150_000.0
+        * (-discount(0.01, start)?
+            + 0.01 * tau * discount(0.01, middle)?
+            + (1.0 + 0.01 * tau) * discount(0.01, end)?);
+    let coupon = (RATE * tau).exp_m1() + 0.002 * tau;
+    let usd = 1_000.0
+        * (discount(RATE, start)?
+            - coupon * discount(RATE, middle)?
+            - (1.0 + coupon) * discount(RATE, end)?);
+    assert_close!(pv.get(Currency::JPY), jpy, PV_REL, JPY_ABS);
+    assert_close!(pv.get(Currency::USD), usd, PV_REL, USD_ABS);
+    Ok(())
+}
diff --git a/crates/rhodium-pricing/tests/xccy_curves.rs b/crates/rhodium-pricing/tests/xccy_curves.rs
index 9bec336156f2731b0444ad745dde1bff01597954..0a8cbbeed7948c4fbf33e60da7d9499e0d84c2b7 100644
--- a/crates/rhodium-pricing/tests/xccy_curves.rs
+++ b/crates/rhodium-pricing/tests/xccy_curves.rs
@@ -15,9 +15,11 @@ use rhodium_core::{
     },
     time::LocalEventTime,
 };
-use rhodium_market::quote::{QuoteId, QuoteSet};
+use rhodium_fixtures::{bumped, solved_within};
+use rhodium_fixtures::{ois, xccy};
+use rhodium_market::quote::QuoteId;
 use rhodium_pricing::{
-    calibration::{CalibratedMarket, GroupOutcome, calibrate},
+    calibration::{CalibratedMarket, calibrate},
     context::Context,
     definition::FxImpliedLongEnd,
     prepare::{PrepareInput, evaluate, prepare},
@@ -25,20 +27,6 @@ use rhodium_pricing::{
 use rhodium_products::standard::{FxForward, FxSwap, Position, StandardProduct};
 use rhodium_testkit::assert_close;
 
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、XCCY の校正に使わない関数がある"
-)]
-#[path = "fixtures/ois.rs"]
-pub mod ois;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、XCCY の校正に使わない関数がある"
-)]
-#[path = "fixtures/xccy.rs"]
-pub mod xccy;
-
 use ois::{Inputs, midnight};
 use xccy::{POINTS, SPOT, XCCY_QUOTES, XccyMarket, jpy_cash};
 
@@ -49,16 +37,6 @@ fn trade_date() -> Result<rhodium_core::date::Date> {
     "2025-01-07".parse()
 }
 
-fn solved(market: &CalibratedMarket) -> Result<()> {
-    for group in market.groups() {
-        let GroupOutcome::Solved(solution) = group.outcome() else {
-            return Err(Error::Calibration(format!("{group:?}")));
-        };
-        assert!(solution.residual_norm() <= TOLERANCE);
-    }
-    Ok(())
-}
-
 /// JPY 建ての PV。USD の受渡を評価時点の X(0) で換算する。
 fn jpy_pv(product: &StandardProduct, market: &CalibratedMarket, ctx: &Context<'_>) -> Result<f64> {
     let plan = prepare(
@@ -100,7 +78,7 @@ fn fx_implied_long_end_reproduces_fx_points_and_xccy_quotes() -> Result<()> {
     let inputs = fixture.inputs("")?;
     let ctx = inputs.context();
     let market = calibrate(&fixture.definition()?, &fixture.quotes()?, &ctx)?;
-    solved(&market)?;
+    solved_within(&market, TOLERANCE)?;
     let block = rhodium_market::structure::BlockId::new(XccyMarket::implied()?, "basis")?;
     let group = market
         .groups()
@@ -164,15 +142,6 @@ fn fx_implied_calibration_instruments_are_at_par_as_trades() -> Result<()> {
     Ok(())
 }
 
-fn bumped(quotes: &QuoteSet, id: &str, bump: f64) -> Result<QuoteSet> {
-    let id: QuoteId = id.parse()?;
-    QuoteSet::new(
-        quotes
-            .iter()
-            .map(|(key, value)| (key.clone(), if *key == id { value + bump } else { value })),
-    )
-}
-
 #[test]
 fn fx_implied_long_end_update_matches_full_recalibration() -> Result<()> {
     // XCCY の quote の変更は外貨担保曲線のグループだけを解き直し、TONA、SOFR と FX のマーク
@@ -198,10 +167,10 @@ fn fx_implied_long_end_update_matches_full_recalibration() -> Result<()> {
             ["USD-SOFR", "FX/USDJPY"].as_slice(),
         ),
     ] {
-        let quotes = bumped(&quotes, id, bump)?;
+        let quotes = bumped(&quotes, &id.parse()?, bump)?;
         let updated = base.update(&quotes)?;
         let full = calibrate(&definition, &quotes, &ctx)?;
-        solved(&updated)?;
+        solved_within(&updated, TOLERANCE)?;
         for (lhs, rhs) in updated.market().theta().iter().zip(full.market().theta()) {
             assert_close!(*lhs, *rhs, 1e-12, 1e-14);
         }
@@ -275,7 +244,7 @@ fn fx_implied_long_end_alone_survives_the_cutoff_of_a_t_plus_one_pair() -> Resul
             ..inputs.context()
         };
         let market = calibrate(&definition, &fixture.quotes()?, &ctx)?;
-        solved(&market)?;
+        solved_within(&market, TOLERANCE)?;
         for (_, tenor, quote, _) in XCCY_QUOTES {
             assert_close!(xccy_par(tenor, &market, &ctx)?, quote, 1e-12_f64, 1e-12_f64);
         }
diff --git a/crates/rhodium-products/src/ir/tests.rs b/crates/rhodium-products/src/ir/tests.rs
index c9abaee0eddb93e784e1d702fe3e4ff154c0b20c..69696c3f7607da29e3430c0fca0fa3f27466dc98 100644
--- a/crates/rhodium-products/src/ir/tests.rs
+++ b/crates/rhodium-products/src/ir/tests.rs
@@ -581,49 +581,6 @@ fn fixing_switches_at_publication_and_past_unpublished_requires_an_estimate() ->
     Ok(())
 }
 
-#[test]
-fn causality_uses_publication_instants_and_checks_each_reachable_branch() -> Result<()> {
-    let refdata = refdata()?;
-    let id: IndexId = "USDJPY-SPOT".parse()?;
-    let publication = "2027-01-08T01:00:00Z".parse::<Instant>()?;
-    let before = Instant::from_unix_micros(publication.unix_micros() - 1)?;
-    let at = "2027-01-12T03:00:00Z".parse()?;
-    for decision in [before, publication] {
-        let contract = Contract::build(|b| {
-            let spot = b.fixing(&id, ObservationPoint::Instant(publication));
-            let strike = b.constant(100.0);
-            let condition = b.gt(spot, strike);
-            let pay = b.pay(CASH, strike, at);
-            let zero = b.zero();
-            Ok(b.cond(condition, decision, pay, zero))
-        })?;
-        assert_eq!(
-            Product::new(contract, &refdata).is_ok(),
-            decision == publication
-        );
-    }
-    let invalid_amount = Contract::build(|b| {
-        let spot = b.fixing(&id, ObservationPoint::Instant(publication));
-        Ok(b.pay(CASH, spot, before))
-    })?;
-    assert!(matches!(
-        Product::new(invalid_amount, &refdata),
-        Err(Error::InvalidInput(_))
-    ));
-    for child_is_early in [true, false] {
-        let contract = Contract::build(|b| {
-            let spot = b.fixing(&id, ObservationPoint::Instant(publication));
-            let pay = b.pay(CASH, spot, if child_is_early { publication } else { at });
-            let zero = b.zero();
-            // 支払額の公表を判断より後にしても、支払時刻までに公表されればよい。
-            let decision = if child_is_early { at } else { before };
-            Ok(b.choose(Party::Us, "exercise", decision, [pay, zero]))
-        })?;
-        assert_eq!(Product::new(contract, &refdata).is_ok(), !child_is_early);
-    }
-    Ok(())
-}
-
 #[test]
 fn exercise_effectiveness_and_settled_payments_control_which_fixings_are_read() -> Result<()> {
     let refdata = refdata()?;
@@ -716,60 +673,6 @@ fn exercise_effectiveness_and_settled_payments_control_which_fixings_are_read()
     Ok(())
 }
 
-#[test]
-fn exercise_records_must_match_the_product_and_deadline() -> Result<()> {
-    let refdata = refdata()?;
-    let cutoff = "2027-01-08T01:00:00Z".parse()?;
-    let future = "2027-01-12T03:00:00Z".parse()?;
-    let product = Product::new(
-        Contract::build(|b| {
-            let value = b.constant(110.0);
-            let pay = b.pay(CASH, value, future);
-            let zero = b.zero();
-            Ok(b.choose(Party::Us, "exercise", cutoff, [pay, zero]))
-        })?,
-        &refdata,
-    )?;
-    let empty = FixingStore::default();
-    // キー、効力時刻、選択肢を商品に照合する。重複は選択肢が同じでも拒否する。
-    let record = ExerciseRecord {
-        label: "exercise".into(),
-        at: cutoff,
-        effective_at: cutoff,
-        option: 1,
-    };
-    let invalid = ExerciseRecord {
-        option: 2,
-        ..record.clone()
-    };
-    let late = ExerciseRecord {
-        effective_at: future,
-        ..record.clone()
-    };
-    let unknown = ExerciseRecord {
-        label: "other".into(),
-        ..record.clone()
-    };
-    for records in [
-        vec![invalid],
-        vec![late],
-        vec![unknown],
-        vec![record.clone(), record],
-    ] {
-        assert!(matches!(
-            reduce(
-                &product,
-                cutoff,
-                &empty,
-                &records,
-                SettlementBoundary::Exclude
-            ),
-            Err(Error::InvalidData(_))
-        ));
-    }
-    Ok(())
-}
-
 #[test]
 fn cond_select_and_logic_skip_unreached_fixings_and_arithmetic() -> Result<()> {
     let refdata = refdata()?;
@@ -1043,10 +946,6 @@ fn undecided_branches_defer_numerical_failures_of_published_aggregates() -> Resu
 }
 
 #[test]
-#[expect(
-    clippy::wildcard_enum_match_arm,
-    reason = "#196 で各 variant を網羅する"
-)]
 fn pending_condition_requires_its_published_fixings() -> Result<()> {
     // 到達した Cond の条件は判定時刻に必ず評価するので、判定が評価時点より後でも、
     // 公表済みの値が欠ければエラーにする。分かる項は畳み、まだ公表されない項だけを残す。
@@ -1086,23 +985,24 @@ fn pending_condition_requires_its_published_fixings() -> Result<()> {
         let ContractNode::Cond { cond, .. } = root(&reduced) else {
             return Err(Error::Internal("expected a pending condition".into()));
         };
-        match &reduced.contract.exprs[*cond as usize] {
+        let condition = &reduced.contract.exprs[*cond as usize];
+        if pending {
             // 150 > 100 は真なので、論理積はまだ公表されない比較だけになる。
-            ExprNode::Binary(nodes::Binary::Gt, [later, _]) => {
-                assert!(pending);
-                assert!(matches!(
-                    reduced.contract.exprs[*later as usize],
-                    ExprNode::Fixing {
-                        point: ObservationPoint::Instant(at),
-                        ..
-                    } if at == gate_at
-                ));
-            }
+            let ExprNode::Binary(nodes::Binary::Gt, [later, _]) = condition else {
+                return Err(Error::Internal(format!(
+                    "unexpected condition {condition:?}"
+                )));
+            };
+            assert!(matches!(
+                reduced.contract.exprs[*later as usize],
+                ExprNode::Fixing {
+                    point: ObservationPoint::Instant(at),
+                    ..
+                } if at == gate_at
+            ));
+        } else {
             // 50 > 100 は偽なので論理積は偽に決まる。判定時刻は先なので Cond は残す。
-            ExprNode::Bool(false) => assert!(!pending),
-            other => {
-                return Err(Error::Internal(format!("unexpected condition {other:?}")));
-            }
+            assert!(matches!(condition, ExprNode::Bool(false)), "{condition:?}");
         }
     }
     Ok(())
@@ -1354,32 +1254,6 @@ fn compound_factors_apply_daily_rounding_and_the_selected_spread_rule() -> Resul
     Ok(())
 }
 
-#[test]
-fn product_rejects_distinct_choices_with_the_same_key_and_empty_choices() -> Result<()> {
-    let refdata = refdata()?;
-    let gate = "2027-01-08T01:00:00Z".parse()?;
-    let at = "2027-01-11T01:00:00Z".parse()?;
-    let duplicate = Contract::build(|b| {
-        let x = b.constant(1.0);
-        let pay = b.pay(CASH, x, at);
-        let zero = b.zero();
-        let first = b.choose(Party::Us, "exercise", gate, [pay, zero]);
-        let second = b.choose(Party::Us, "exercise", gate, [zero, pay]);
-        Ok(b.all([first, second]))
-    })?;
-    // 同じ行使記録が異なる選択に一致する構成と、選べる枝のない構成を拒否する。
-    assert!(matches!(
-        Product::new(duplicate, &refdata),
-        Err(Error::InvalidInput(_))
-    ));
-    let empty = Contract::build(|b| Ok(b.choose(Party::Us, "exercise", gate, [])))?;
-    assert!(matches!(
-        Product::new(empty, &refdata),
-        Err(Error::InvalidInput(_))
-    ));
-    Ok(())
-}
-
 #[test]
 fn give_flips_signs_and_parties_by_parity() -> Result<()> {
     // 旧 IR の give の例を移す。give の偶奇は額の符号と選択の当事者だけに作用し、
@@ -1554,351 +1428,6 @@ fn node_budget_rejects_contracts_over_the_limit() -> Result<()> {
     Ok(())
 }
 
-#[test]
-fn observation_builders_reject_rules_that_did_not_make_the_schedule() -> Result<()> {
-    // 旧 IR の観測列の検査の例のうち、新しい builder に当たるものを移す。
-    // 複利は観測列を作った複利の規則とだけ、平均は平均の規則で作った観測列とだけ組める。
-    let refdata = refdata()?;
-    let id: IndexId = "JPY-TONA".parse()?;
-    let at = "2027-01-12T03:00:00Z".parse()?;
-    let compound = Aggregation::Compound {
-        form: CompoundingForm::SimpleFactor,
-        spread: SpreadTreatment::Exclusive,
-    };
-    let lookback = observation_rows(
-        &refdata,
-        ObservationMethod::Lookback(NonZeroU8::MIN),
-        compound,
-        ObservationRounding::default(),
-    )?;
-    let plain = observation_rows(
-        &refdata,
-        ObservationMethod::Plain,
-        compound,
-        ObservationRounding::default(),
-    )?;
-    let average = observation_rows(
-        &refdata,
-        ObservationMethod::Lookback(NonZeroU8::MIN),
-        Aggregation::Average(AveragingWeights::Equal),
-        ObservationRounding::default(),
-    )?;
-    for (observations, rule) in [(&lookback, plain.rule()), (&average, average.rule())] {
-        let result = Contract::build(|b| {
-            let gain = b.compounded(&id, observations, rule, Spread::new(0.0)?)?;
-            Ok(b.pay(CASH, gain, at))
-        });
-        assert!(matches!(result, Err(Error::InvalidInput(_))));
-    }
-    let result = Contract::build(|b| {
-        let mean = b.average(&id, &lookback)?;
-        Ok(b.pay(CASH, mean, at))
-    });
-    assert!(matches!(result, Err(Error::InvalidInput(_))));
-    Ok(())
-}
-
-mod lifecycle;
-mod standard;
-
-#[test]
-fn analysis_preserves_payment_occurrences_and_branch_paths() -> Result<()> {
-    // 同じ支払ノードが繰り返し現れても支払回数を保つ。式の観測だけを重複排除する。
-    use analysis::{Branch, Linearity, Observation, cashflows, linearity, observations};
-
-    let data = refdata()?;
-    let spot: IndexId = "USDJPY-SPOT".parse()?;
-    let observed = "2027-02-01T00:00:00Z".parse()?;
-    let decided = "2027-02-02T00:00:00Z".parse()?;
-    let paid = "2027-02-03T00:00:00Z".parse()?;
-    let point = ObservationPoint::Instant(observed);
-    let contract = Contract::build(|builder| {
-        let fixing = builder.fixing(&spot, point);
-        let one = builder.constant(1.0);
-        let pay = builder.pay(CASH, one, paid);
-        let conditional_amount = builder.pay(CASH, fixing, paid);
-        let choice = builder.choose(Party::Us, "exercise", decided, [conditional_amount, pay]);
-        let condition = builder.gt(fixing, one);
-        let contingent = builder.cond(condition, decided, choice, pay);
-        Ok(builder.all([pay, pay, contingent]))
-    })?;
-    let product = Product::new(contract, &data)?;
-    let flows = cashflows(&product);
-    assert_eq!(flows.len(), 5);
-    assert_eq!(flows[0].branches.len(), 0);
-    assert_eq!(flows[1].amount.index(), flows[0].amount.index());
-    assert!(matches!(
-        flows[2].branches.as_slice(),
-        [
-            Branch::Condition { when: true, .. },
-            Branch::Choice {
-                label: "exercise",
-                option: 0,
-                ..
-            }
-        ]
-    ));
-    assert!(matches!(
-        flows[4].branches.as_slice(),
-        [Branch::Condition { when: false, .. }]
-    ));
-    assert_eq!(
-        observations(&product),
-        [Observation::Fixing {
-            index: &spot,
-            point
-        }]
-    );
-    assert!(matches!(linearity(&product), Linearity::NonLinear(_)));
-    Ok(())
-}
-
-/// 予想する複利 g と定数 2 から支払額を組み、リニア方式への適用可否を判定する。
-fn linearity_of(
-    refdata: &RefData,
-    observations: &ObservationSchedule,
-    amount: impl for<'b> FnOnce(&mut Builder<'b>, ExprIx<'b>, ExprIx<'b>) -> ExprIx<'b>,
-) -> Result<analysis::Linearity> {
-    let id: IndexId = "JPY-TONA".parse()?;
-    let at = "2027-01-12T03:00:00Z".parse()?;
-    let product = Product::new(
-        Contract::build(|b| {
-            let growth = b.compounded(&id, observations, observations.rule(), Spread::new(0.0)?)?;
-            let constant = b.constant(2.0);
-            let amount = amount(b, growth, constant);
-            Ok(b.pay(CASH, amount, at))
-        })?,
-        refdata,
-    )?;
-    Ok(analysis::linearity(&product))
-}
-
-#[test]
-fn linearity_rejects_products_and_quotients_of_projected_growth() -> Result<()> {
-    // 予想する複利 g の定数倍、定数での割り算、丸め、符号の反転と差は、g について線形なので
-    // リニア評価で扱える。g 同士の積、g による割り算、g を含む最大値は線形でないため、
-    // 方式の選択がその理由とともに Unsupported を返せるよう NonLinear にする。
-    use analysis::Linearity;
-
-    let refdata = refdata()?;
-    let observations = observation_rows(
-        &refdata,
-        ObservationMethod::Lookback(NonZeroU8::MIN),
-        Aggregation::Compound {
-            form: CompoundingForm::SimpleFactor,
-            spread: SpreadTreatment::Exclusive,
-        },
-        ObservationRounding::default(),
-    )?;
-    let rounding = Rounding::new(0.0001, RoundingMode::HalfEven)?;
-
-    let scaled = linearity_of(&refdata, &observations, |b, g, c| {
-        let notional = b.constant(1_000.0);
-        let term = b.mul([g, notional]);
-        b.add([term, c])
-    })?;
-    assert_eq!(scaled, Linearity::Linear);
-    let divided = linearity_of(&refdata, &observations, |b, g, _| {
-        let denominator = b.constant(4.0);
-        b.div(g, denominator)
-    })?;
-    assert_eq!(divided, Linearity::Linear);
-    let rounded = linearity_of(&refdata, &observations, |b, g, _| b.round(g, rounding))?;
-    assert_eq!(rounded, Linearity::Linear);
-    let negated = linearity_of(&refdata, &observations, |b, g, c| {
-        let negative = b.neg(g);
-        b.sub(negative, c)
-    })?;
-    assert_eq!(negated, Linearity::Linear);
-
-    let squared = linearity_of(&refdata, &observations, |b, g, _| b.mul([g, g]))?;
-    assert_eq!(
-        squared,
-        Linearity::NonLinear("multiplication of two projected quantities")
-    );
-    let reciprocal = linearity_of(&refdata, &observations, |b, g, c| b.div(c, g))?;
-    assert_eq!(
-        reciprocal,
-        Linearity::NonLinear("division by a projected quantity")
-    );
-    let floored = linearity_of(&refdata, &observations, |b, g, c| b.max([g, c]))?;
-    assert_eq!(
-        floored,
-        Linearity::NonLinear("unsupported cashflow expression")
-    );
-    Ok(())
-}
-
-#[test]
-fn cond_rejects_a_branch_paying_before_its_decision() -> Result<()> {
-    // Cond は判断の瞬間より前に受け渡す枝を持てない（Choose と同じ因果性の規則）。
-    let refdata = refdata()?;
-    let id: IndexId = "USDJPY-SPOT".parse()?;
-    let publication = "2027-01-08T01:00:00Z".parse::<Instant>()?;
-    let before = Instant::from_unix_micros(publication.unix_micros() - 1)?;
-    for (paid, accepted) in [(before, false), (publication, true)] {
-        let contract = Contract::build(|b| {
-            let spot = b.fixing(&id, ObservationPoint::Instant(publication));
-            let strike = b.constant(100.0);
-            let condition = b.gt(spot, strike);
-            let pay = b.pay(CASH, strike, paid);
-            let zero = b.zero();
-            Ok(b.cond(condition, publication, pay, zero))
-        })?;
-        let product = Product::new(contract, &refdata);
-        if accepted {
-            // 判断と同じ瞬間の受渡は許す。
-            product?;
-        } else {
-            assert!(
-                matches!(&product, Err(Error::InvalidInput(message)) if message.contains("precedes decision")),
-                "{product:?}"
-            );
-        }
-    }
-    Ok(())
-}
-
-// Product::new が契約を拒否し、エラーの種類が期待どおりであることを確かめる。
-fn assert_rejected(
-    refdata: &RefData,
-    what: &str,
-    expected: fn(&Error) -> bool,
-    contract: Contract,
-) {
-    let result = Product::new(contract, refdata);
-    assert!(
-        matches!(&result, Err(error) if expected(error)),
-        "{what}: {result:?}"
-    );
-}
-
-fn is_input(error: &Error) -> bool {
-    matches!(error, Error::InvalidInput(_))
-}
-
-#[test]
-fn product_new_rejects_malformed_expressions() -> Result<()> {
-    // 非有限の定数と、項の無い min・max は builder で組めても検査で拒否する。
-    let refdata = refdata()?;
-    let at = "2027-01-12T03:00:00Z".parse()?;
-    assert_rejected(
-        &refdata,
-        "non-finite constant",
-        is_input,
-        Contract::build(|b| {
-            let nan = b.constant(f64::NAN);
-            Ok(b.pay(CASH, nan, at))
-        })?,
-    );
-    assert_rejected(
-        &refdata,
-        "empty min",
-        is_input,
-        Contract::build(|b| {
-            let min = b.min([]);
-            Ok(b.pay(CASH, min, at))
-        })?,
-    );
-    assert_rejected(
-        &refdata,
-        "empty max",
-        is_input,
-        Contract::build(|b| {
-            let max = b.max([]);
-            Ok(b.pay(CASH, max, at))
-        })?,
-    );
-    Ok(())
-}
-
-#[test]
-fn product_new_rejects_index_misuse() -> Result<()> {
-    // 指標の種類に合わない読み方を、参照データと照合して一度だけ拒否する。
-    let refdata = refdata()?;
-    let at = "2027-01-12T03:00:00Z".parse()?;
-    let event: IndexId = "JP-EVENT".parse()?;
-    let recovery: IndexId = "JP-RECOVERY".parse()?;
-    let day = ObservationPoint::Date("2027-01-08".parse()?);
-    let window = DateInterval::new("2027-01-08".parse()?, "2027-01-11".parse()?)?;
-    let empty = DateInterval::new("2027-01-08".parse()?, "2027-01-08".parse()?)?;
-    let rows = observation_rows(
-        &refdata,
-        ObservationMethod::Plain,
-        Aggregation::Compound {
-            form: CompoundingForm::SimpleFactor,
-            spread: SpreadTreatment::Exclusive,
-        },
-        ObservationRounding::default(),
-    )?;
-    // 信用事由の指標は数の fixing として読めず、信用窓が要る。
-    assert_rejected(
-        &refdata,
-        "credit index as a numeric fixing",
-        is_input,
-        Contract::build(|b| {
-            let value = b.fixing(&event, day);
-            Ok(b.pay(CASH, value, at))
-        })?,
-    );
-    // 複利は翌日物の指標だけに使える。FX の fixing は数でも複利にできない。
-    let compounded = |index: &IndexId| {
-        Contract::build(|b| {
-            let growth = b.compounded(index, &rows, rows.rule(), Spread::new(0.0)?)?;
-            Ok(b.pay(CASH, growth, at))
-        })
-    };
-    assert_rejected(
-        &refdata,
-        "compounding a non-overnight index",
-        is_input,
-        compounded(&"USDJPY-SPOT".parse()?)?,
-    );
-    Product::new(compounded(&"JPY-TONA".parse()?)?, &refdata)?;
-    // 事由の判定と回収率は、それぞれ対応する種類の指標と空でない窓で読む。
-    let defaulted = |index: &IndexId, window| {
-        Contract::build(|b| {
-            let defaulted = b.defaulted(index, window);
-            let one = b.constant(1.0);
-            let zero = b.constant(0.0);
-            let amount = b.select(defaulted, one, zero);
-            Ok(b.pay(CASH, amount, at))
-        })
-    };
-    assert_rejected(
-        &refdata,
-        "defaulted on a recovery index",
-        is_input,
-        defaulted(&recovery, window)?,
-    );
-    assert_rejected(
-        &refdata,
-        "empty credit window",
-        is_input,
-        defaulted(&event, empty)?,
-    );
-    assert_rejected(
-        &refdata,
-        "recovery on an event index",
-        is_input,
-        Contract::build(|b| {
-            let rate = b.recovery(&event, window);
-            Ok(b.pay(CASH, rate, at))
-        })?,
-    );
-    // 参照データに無い指標は、データの作り手が直す参照の欠け。
-    assert_rejected(
-        &refdata,
-        "unknown index",
-        |error| matches!(error, Error::InvalidData(_)),
-        Contract::build(|b| {
-            let value = b.fixing(&"EUR-ESTR".parse()?, day);
-            Ok(b.pay(CASH, value, at))
-        })?,
-    );
-    Ok(())
-}
-
 #[test]
 fn tona_fixing_switches_at_publication_and_forecasts_a_spanning_period() -> Result<()> {
     // 金曜 2027-01-08 の TONA。観測期間は金曜 09:00 から月曜 09:00（JST）まで、
diff --git a/crates/rhodium-products/tests/ir_validation.rs b/crates/rhodium-products/tests/ir_validation.rs
new file mode 100644
index 0000000000000000000000000000000000000000..ad05a004678e52ca0cb37708681e1f4929ef92b1
--- /dev/null
+++ b/crates/rhodium-products/tests/ir_validation.rs
@@ -0,0 +1,562 @@
+//! 契約の組み立てと `Product::new` の検査、検証済み IR の解析を公開 API で確かめる。
+
+pub mod support;
+
+use std::num::NonZeroU8;
+
+use rhodium_core::{
+    Error, Result,
+    calendar::CalendarSpec,
+    currency::Currency,
+    date::DateInterval,
+    day_count::{AccrualContext, DayCount},
+    fixing::{FixingStore, ObservationPoint},
+    id::IndexId,
+    observation::{
+        Aggregation, AveragingWeights, CompoundingForm, ObservationDays, ObservationMethod,
+        ObservationRounding, ObservationRule, ObservationSchedule, SpreadTreatment,
+        observation_schedule,
+    },
+    quantity::Spread,
+    refdata::RefData,
+    rounding::{Rounding, RoundingMode},
+    schedule::{
+        ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination, ScheduleSpec,
+        SettlementContext, generate,
+    },
+    time::Instant,
+};
+use rhodium_products::{
+    ir::{Builder, Contract, ExprIx, Party, Product, SettlementBoundary, analysis, reduce},
+    trade::ExerciseRecord,
+};
+use support::{CASH, refdata};
+
+// 2027-01-08 から 01-12 の一期間を JPTO の営業日で観測する列。複利か平均かと観測の規約を選べる。
+fn observation_rows(
+    refdata: &RefData,
+    method: ObservationMethod,
+    aggregation: Aggregation,
+    rounding: ObservationRounding,
+) -> Result<ObservationSchedule> {
+    let accrual = DateInterval::new("2027-01-08".parse()?, "2027-01-12".parse()?)?;
+    let spec = ScheduleSpec::Explicit(ExplicitPeriods::new(
+        vec![ExplicitPeriod {
+            accrual,
+            payment: accrual.end(),
+            reset: None,
+            settlement: None,
+            references: Vec::new(),
+        }],
+        ExplicitDatePolicy::KeepAsConfirmed,
+        PeriodCombination::Sum,
+        None,
+    )?);
+    let calendar = CalendarSpec::Base("JPTO".parse()?);
+    let (schedule, _) = generate(
+        &spec,
+        accrual.start(),
+        refdata.calendars(),
+        &SettlementContext {
+            zones: refdata.zones(),
+            default_time: Some(&refdata.currency(Currency::JPY)?.settlement_time),
+            calendar: Some(&calendar),
+        },
+    )?;
+    let rule = ObservationRule::new(
+        ObservationDays::BusinessDays(CalendarSpec::Base("JPTO".parse()?)),
+        method,
+        aggregation,
+        rounding,
+    )?;
+    observation_schedule(
+        &schedule.periods()[0],
+        &rule,
+        &DayCount::Act365F,
+        &AccrualContext::default(),
+        refdata.calendars(),
+    )
+    .map(|result| result.0)
+}
+
+#[test]
+fn causality_uses_publication_instants_and_checks_each_reachable_branch() -> Result<()> {
+    // 判断と支払額が読む fixing は、その瞬間までに公表されていなければならない。公表ちょうどの判断は
+    // 通り、1µs 前の判断と公表前の支払いは拒否する。選択の枝は到達しうるものをそれぞれ検査する。
+    // 拒否はどれも因果性の違反なので、エラーの種類（InvalidInput）まで確かめる。
+    let refdata = refdata()?;
+    let id: IndexId = "USDJPY-SPOT".parse()?;
+    let publication = "2027-01-08T01:00:00Z".parse::<Instant>()?;
+    let before = Instant::from_unix_micros(publication.unix_micros() - 1)?;
+    let at = "2027-01-12T03:00:00Z".parse()?;
+    for decision in [before, publication] {
+        let contract = Contract::build(|b| {
+            let spot = b.fixing(&id, ObservationPoint::Instant(publication));
+            let strike = b.constant(100.0);
+            let condition = b.gt(spot, strike);
+            let pay = b.pay(CASH, strike, at);
+            let zero = b.zero();
+            Ok(b.cond(condition, decision, pay, zero))
+        })?;
+        let product = Product::new(contract, &refdata);
+        if decision == publication {
+            product?;
+        } else {
+            assert!(
+                matches!(&product, Err(Error::InvalidInput(_))),
+                "{product:?}"
+            );
+        }
+    }
+    let invalid_amount = Contract::build(|b| {
+        let spot = b.fixing(&id, ObservationPoint::Instant(publication));
+        Ok(b.pay(CASH, spot, before))
+    })?;
+    assert!(matches!(
+        Product::new(invalid_amount, &refdata),
+        Err(Error::InvalidInput(_))
+    ));
+    for child_is_early in [true, false] {
+        let contract = Contract::build(|b| {
+            let spot = b.fixing(&id, ObservationPoint::Instant(publication));
+            let pay = b.pay(CASH, spot, if child_is_early { publication } else { at });
+            let zero = b.zero();
+            // 支払額の公表を判断より後にしても、支払時刻までに公表されればよい。
+            let decision = if child_is_early { at } else { before };
+            Ok(b.choose(Party::Us, "exercise", decision, [pay, zero]))
+        })?;
+        let product = Product::new(contract, &refdata);
+        if child_is_early {
+            // 判断より前に受け渡す枝を持つ選択は拒否する。
+            assert!(
+                matches!(&product, Err(Error::InvalidInput(_))),
+                "{product:?}"
+            );
+        } else {
+            product?;
+        }
+    }
+    Ok(())
+}
+
+#[test]
+fn exercise_records_must_match_the_product_and_deadline() -> Result<()> {
+    let refdata = refdata()?;
+    let cutoff = "2027-01-08T01:00:00Z".parse()?;
+    let future = "2027-01-12T03:00:00Z".parse()?;
+    let product = Product::new(
+        Contract::build(|b| {
+            let value = b.constant(110.0);
+            let pay = b.pay(CASH, value, future);
+            let zero = b.zero();
+            Ok(b.choose(Party::Us, "exercise", cutoff, [pay, zero]))
+        })?,
+        &refdata,
+    )?;
+    let empty = FixingStore::default();
+    // キー、効力時刻、選択肢を商品に照合する。重複は選択肢が同じでも拒否する。
+    let record = ExerciseRecord {
+        label: "exercise".into(),
+        at: cutoff,
+        effective_at: cutoff,
+        option: 1,
+    };
+    let invalid = ExerciseRecord {
+        option: 2,
+        ..record.clone()
+    };
+    let late = ExerciseRecord {
+        effective_at: future,
+        ..record.clone()
+    };
+    let unknown = ExerciseRecord {
+        label: "other".into(),
+        ..record.clone()
+    };
+    for records in [
+        vec![invalid],
+        vec![late],
+        vec![unknown],
+        vec![record.clone(), record],
+    ] {
+        assert!(matches!(
+            reduce(
+                &product,
+                cutoff,
+                &empty,
+                &records,
+                SettlementBoundary::Exclude
+            ),
+            Err(Error::InvalidData(_))
+        ));
+    }
+    Ok(())
+}
+
+#[test]
+fn product_rejects_distinct_choices_with_the_same_key_and_empty_choices() -> Result<()> {
+    let refdata = refdata()?;
+    let gate = "2027-01-08T01:00:00Z".parse()?;
+    let at = "2027-01-11T01:00:00Z".parse()?;
+    let duplicate = Contract::build(|b| {
+        let x = b.constant(1.0);
+        let pay = b.pay(CASH, x, at);
+        let zero = b.zero();
+        let first = b.choose(Party::Us, "exercise", gate, [pay, zero]);
+        let second = b.choose(Party::Us, "exercise", gate, [zero, pay]);
+        Ok(b.all([first, second]))
+    })?;
+    // 同じ行使記録が異なる選択に一致する構成と、選べる枝のない構成を拒否する。
+    assert!(matches!(
+        Product::new(duplicate, &refdata),
+        Err(Error::InvalidInput(_))
+    ));
+    let empty = Contract::build(|b| Ok(b.choose(Party::Us, "exercise", gate, [])))?;
+    assert!(matches!(
+        Product::new(empty, &refdata),
+        Err(Error::InvalidInput(_))
+    ));
+    Ok(())
+}
+
+#[test]
+fn observation_builders_reject_rules_that_did_not_make_the_schedule() -> Result<()> {
+    // 旧 IR の観測列の検査の例のうち、新しい builder に当たるものを移す。
+    // 複利は観測列を作った複利の規則とだけ、平均は平均の規則で作った観測列とだけ組める。
+    let refdata = refdata()?;
+    let id: IndexId = "JPY-TONA".parse()?;
+    let at = "2027-01-12T03:00:00Z".parse()?;
+    let compound = Aggregation::Compound {
+        form: CompoundingForm::SimpleFactor,
+        spread: SpreadTreatment::Exclusive,
+    };
+    let lookback = observation_rows(
+        &refdata,
+        ObservationMethod::Lookback(NonZeroU8::MIN),
+        compound,
+        ObservationRounding::default(),
+    )?;
+    let plain = observation_rows(
+        &refdata,
+        ObservationMethod::Plain,
+        compound,
+        ObservationRounding::default(),
+    )?;
+    let average = observation_rows(
+        &refdata,
+        ObservationMethod::Lookback(NonZeroU8::MIN),
+        Aggregation::Average(AveragingWeights::Equal),
+        ObservationRounding::default(),
+    )?;
+    for (observations, rule) in [(&lookback, plain.rule()), (&average, average.rule())] {
+        let result = Contract::build(|b| {
+            let gain = b.compounded(&id, observations, rule, Spread::new(0.0)?)?;
+            Ok(b.pay(CASH, gain, at))
+        });
+        assert!(matches!(result, Err(Error::InvalidInput(_))));
+    }
+    let result = Contract::build(|b| {
+        let mean = b.average(&id, &lookback)?;
+        Ok(b.pay(CASH, mean, at))
+    });
+    assert!(matches!(result, Err(Error::InvalidInput(_))));
+    Ok(())
+}
+
+#[test]
+fn analysis_preserves_payment_occurrences_and_branch_paths() -> Result<()> {
+    // 同じ支払ノードが繰り返し現れても支払回数を保つ。式の観測だけを重複排除する。
+    use analysis::{Branch, Linearity, Observation, cashflows, linearity, observations};
+
+    let data = refdata()?;
+    let spot: IndexId = "USDJPY-SPOT".parse()?;
+    let observed = "2027-02-01T00:00:00Z".parse()?;
+    let decided = "2027-02-02T00:00:00Z".parse()?;
+    let paid = "2027-02-03T00:00:00Z".parse()?;
+    let point = ObservationPoint::Instant(observed);
+    let contract = Contract::build(|builder| {
+        let fixing = builder.fixing(&spot, point);
+        let one = builder.constant(1.0);
+        let pay = builder.pay(CASH, one, paid);
+        let conditional_amount = builder.pay(CASH, fixing, paid);
+        let choice = builder.choose(Party::Us, "exercise", decided, [conditional_amount, pay]);
+        let condition = builder.gt(fixing, one);
+        let contingent = builder.cond(condition, decided, choice, pay);
+        Ok(builder.all([pay, pay, contingent]))
+    })?;
+    let product = Product::new(contract, &data)?;
+    let flows = cashflows(&product);
+    assert_eq!(flows.len(), 5);
+    assert_eq!(flows[0].branches.len(), 0);
+    assert_eq!(flows[1].amount.index(), flows[0].amount.index());
+    assert!(matches!(
+        flows[2].branches.as_slice(),
+        [
+            Branch::Condition { when: true, .. },
+            Branch::Choice {
+                label: "exercise",
+                option: 0,
+                ..
+            }
+        ]
+    ));
+    assert!(matches!(
+        flows[4].branches.as_slice(),
+        [Branch::Condition { when: false, .. }]
+    ));
+    assert_eq!(
+        observations(&product),
+        [Observation::Fixing {
+            index: &spot,
+            point
+        }]
+    );
+    assert!(matches!(linearity(&product), Linearity::NonLinear(_)));
+    Ok(())
+}
+
+/// 予想する複利 g と定数 2 から支払額を組み、リニア方式への適用可否を判定する。
+fn linearity_of(
+    refdata: &RefData,
+    observations: &ObservationSchedule,
+    amount: impl for<'b> FnOnce(&mut Builder<'b>, ExprIx<'b>, ExprIx<'b>) -> ExprIx<'b>,
+) -> Result<analysis::Linearity> {
+    let id: IndexId = "JPY-TONA".parse()?;
+    let at = "2027-01-12T03:00:00Z".parse()?;
+    let product = Product::new(
+        Contract::build(|b| {
+            let growth = b.compounded(&id, observations, observations.rule(), Spread::new(0.0)?)?;
+            let constant = b.constant(2.0);
+            let amount = amount(b, growth, constant);
+            Ok(b.pay(CASH, amount, at))
+        })?,
+        refdata,
+    )?;
+    Ok(analysis::linearity(&product))
+}
+
+#[test]
+fn linearity_rejects_products_and_quotients_of_projected_growth() -> Result<()> {
+    // 予想する複利 g の定数倍、定数での割り算、丸め、符号の反転と差は、g について線形なので
+    // リニア評価で扱える。g 同士の積、g による割り算、g を含む最大値は線形でないため、
+    // 方式の選択がその理由とともに Unsupported を返せるよう NonLinear にする。
+    use analysis::Linearity;
+
+    let refdata = refdata()?;
+    let observations = observation_rows(
+        &refdata,
+        ObservationMethod::Lookback(NonZeroU8::MIN),
+        Aggregation::Compound {
+            form: CompoundingForm::SimpleFactor,
+            spread: SpreadTreatment::Exclusive,
+        },
+        ObservationRounding::default(),
+    )?;
+    let rounding = Rounding::new(0.0001, RoundingMode::HalfEven)?;
+
+    let scaled = linearity_of(&refdata, &observations, |b, g, c| {
+        let notional = b.constant(1_000.0);
+        let term = b.mul([g, notional]);
+        b.add([term, c])
+    })?;
+    assert_eq!(scaled, Linearity::Linear);
+    let divided = linearity_of(&refdata, &observations, |b, g, _| {
+        let denominator = b.constant(4.0);
+        b.div(g, denominator)
+    })?;
+    assert_eq!(divided, Linearity::Linear);
+    let rounded = linearity_of(&refdata, &observations, |b, g, _| b.round(g, rounding))?;
+    assert_eq!(rounded, Linearity::Linear);
+    let negated = linearity_of(&refdata, &observations, |b, g, c| {
+        let negative = b.neg(g);
+        b.sub(negative, c)
+    })?;
+    assert_eq!(negated, Linearity::Linear);
+
+    let squared = linearity_of(&refdata, &observations, |b, g, _| b.mul([g, g]))?;
+    assert_eq!(
+        squared,
+        Linearity::NonLinear("multiplication of two projected quantities")
+    );
+    let reciprocal = linearity_of(&refdata, &observations, |b, g, c| b.div(c, g))?;
+    assert_eq!(
+        reciprocal,
+        Linearity::NonLinear("division by a projected quantity")
+    );
+    let floored = linearity_of(&refdata, &observations, |b, g, c| b.max([g, c]))?;
+    assert_eq!(
+        floored,
+        Linearity::NonLinear("unsupported cashflow expression")
+    );
+    Ok(())
+}
+
+#[test]
+fn cond_rejects_a_branch_paying_before_its_decision() -> Result<()> {
+    // Cond は判断の瞬間より前に受け渡す枝を持てない（Choose と同じ因果性の規則）。
+    let refdata = refdata()?;
+    let id: IndexId = "USDJPY-SPOT".parse()?;
+    let publication = "2027-01-08T01:00:00Z".parse::<Instant>()?;
+    let before = Instant::from_unix_micros(publication.unix_micros() - 1)?;
+    for (paid, accepted) in [(before, false), (publication, true)] {
+        let contract = Contract::build(|b| {
+            let spot = b.fixing(&id, ObservationPoint::Instant(publication));
+            let strike = b.constant(100.0);
+            let condition = b.gt(spot, strike);
+            let pay = b.pay(CASH, strike, paid);
+            let zero = b.zero();
+            Ok(b.cond(condition, publication, pay, zero))
+        })?;
+        let product = Product::new(contract, &refdata);
+        if accepted {
+            // 判断と同じ瞬間の受渡は許す。
+            product?;
+        } else {
+            assert!(
+                matches!(&product, Err(Error::InvalidInput(message)) if message.contains("precedes decision")),
+                "{product:?}"
+            );
+        }
+    }
+    Ok(())
+}
+
+// Product::new が契約を拒否し、エラーの種類が期待どおりであることを確かめる。
+fn assert_rejected(
+    refdata: &RefData,
+    what: &str,
+    expected: fn(&Error) -> bool,
+    contract: Contract,
+) {
+    let result = Product::new(contract, refdata);
+    assert!(
+        matches!(&result, Err(error) if expected(error)),
+        "{what}: {result:?}"
+    );
+}
+
+fn is_input(error: &Error) -> bool {
+    matches!(error, Error::InvalidInput(_))
+}
+
+#[test]
+fn product_new_rejects_malformed_expressions() -> Result<()> {
+    // 非有限の定数と、項の無い min・max は builder で組めても検査で拒否する。
+    let refdata = refdata()?;
+    let at = "2027-01-12T03:00:00Z".parse()?;
+    assert_rejected(
+        &refdata,
+        "non-finite constant",
+        is_input,
+        Contract::build(|b| {
+            let nan = b.constant(f64::NAN);
+            Ok(b.pay(CASH, nan, at))
+        })?,
+    );
+    assert_rejected(
+        &refdata,
+        "empty min",
+        is_input,
+        Contract::build(|b| {
+            let min = b.min([]);
+            Ok(b.pay(CASH, min, at))
+        })?,
+    );
+    assert_rejected(
+        &refdata,
+        "empty max",
+        is_input,
+        Contract::build(|b| {
+            let max = b.max([]);
+            Ok(b.pay(CASH, max, at))
+        })?,
+    );
+    Ok(())
+}
+
+#[test]
+fn product_new_rejects_index_misuse() -> Result<()> {
+    // 指標の種類に合わない読み方を、参照データと照合して一度だけ拒否する。
+    let refdata = refdata()?;
+    let at = "2027-01-12T03:00:00Z".parse()?;
+    let event: IndexId = "JP-EVENT".parse()?;
+    let recovery: IndexId = "JP-RECOVERY".parse()?;
+    let day = ObservationPoint::Date("2027-01-08".parse()?);
+    let window = DateInterval::new("2027-01-08".parse()?, "2027-01-11".parse()?)?;
+    let empty = DateInterval::new("2027-01-08".parse()?, "2027-01-08".parse()?)?;
+    let rows = observation_rows(
+        &refdata,
+        ObservationMethod::Plain,
+        Aggregation::Compound {
+            form: CompoundingForm::SimpleFactor,
+            spread: SpreadTreatment::Exclusive,
+        },
+        ObservationRounding::default(),
+    )?;
+    // 信用事由の指標は数の fixing として読めず、信用窓が要る。
+    assert_rejected(
+        &refdata,
+        "credit index as a numeric fixing",
+        is_input,
+        Contract::build(|b| {
+            let value = b.fixing(&event, day);
+            Ok(b.pay(CASH, value, at))
+        })?,
+    );
+    // 複利は翌日物の指標だけに使える。FX の fixing は数でも複利にできない。
+    let compounded = |index: &IndexId| {
+        Contract::build(|b| {
+            let growth = b.compounded(index, &rows, rows.rule(), Spread::new(0.0)?)?;
+            Ok(b.pay(CASH, growth, at))
+        })
+    };
+    assert_rejected(
+        &refdata,
+        "compounding a non-overnight index",
+        is_input,
+        compounded(&"USDJPY-SPOT".parse()?)?,
+    );
+    Product::new(compounded(&"JPY-TONA".parse()?)?, &refdata)?;
+    // 事由の判定と回収率は、それぞれ対応する種類の指標と空でない窓で読む。
+    let defaulted = |index: &IndexId, window| {
+        Contract::build(|b| {
+            let defaulted = b.defaulted(index, window);
+            let one = b.constant(1.0);
+            let zero = b.constant(0.0);
+            let amount = b.select(defaulted, one, zero);
+            Ok(b.pay(CASH, amount, at))
+        })
+    };
+    assert_rejected(
+        &refdata,
+        "defaulted on a recovery index",
+        is_input,
+        defaulted(&recovery, window)?,
+    );
+    assert_rejected(
+        &refdata,
+        "empty credit window",
+        is_input,
+        defaulted(&event, empty)?,
+    );
+    assert_rejected(
+        &refdata,
+        "recovery on an event index",
+        is_input,
+        Contract::build(|b| {
+            let rate = b.recovery(&event, window);
+            Ok(b.pay(CASH, rate, at))
+        })?,
+    );
+    // 参照データに無い指標は、データの作り手が直す参照の欠け。
+    assert_rejected(
+        &refdata,
+        "unknown index",
+        |error| matches!(error, Error::InvalidData(_)),
+        Contract::build(|b| {
+            let value = b.fixing(&"EUR-ESTR".parse()?, day);
+            Ok(b.pay(CASH, value, at))
+        })?,
+    );
+    Ok(())
+}
diff --git a/crates/rhodium-products/src/ir/tests/lifecycle.rs b/crates/rhodium-products/tests/lifecycle.rs
similarity index 86%
rename from crates/rhodium-products/src/ir/tests/lifecycle.rs
rename to crates/rhodium-products/tests/lifecycle.rs
index 45569717f9c44234bd2b967d83ca4a22037b72a3..727bfb7f777bae02164497dcc9b614211876bc97 100644
--- a/crates/rhodium-products/src/ir/tests/lifecycle.rs
+++ b/crates/rhodium-products/tests/lifecycle.rs
@@ -1,19 +1,41 @@
+//! 取引の解約、差し替え、行使とフォールバックを、取引から IR を作る公開の入口で確かめる。
+
+pub mod support;
+
 use std::collections::BTreeMap;
 
 use rhodium_core::{
-    calendar::{BusinessDayConvention, EndOfMonthRule},
+    Error, Result, Warning,
+    calendar::{BusinessDayConvention, CalendarSpec, EndOfMonthRule},
+    currency::Currency,
     date::Date,
+    day_count::DayCount,
+    fixing::{FixingStore, FixingStoreParts, ObservationPoint},
+    id::IndexId,
     money::Money,
-    refdata::{FallbackApplicability, FallbackDefinition, FallbackReplacement, FallbackTrigger},
+    observation::{
+        Aggregation, CompoundingForm, ObservationDays, ObservationMethod, ObservationRounding,
+        ObservationRule, SpreadTreatment,
+    },
+    refdata::{
+        FallbackApplicability, FallbackDefinition, FallbackReplacement, FallbackTrigger, IndexKind,
+        IndexStatus, PublicationRule, RefData,
+    },
+    time::Instant,
 };
-
-use super::*;
-use crate::{
-    ir::analysis::{Observation, cashflows, observations},
+use rhodium_products::{
+    ir::{
+        Contract, Party, Product, SettlementBoundary,
+        analysis::{Expression, Observation, cashflows, observations},
+        reduce,
+    },
     standard::PayReceive,
-    trade::{FallbackSelection, FallbackTerms, Fee, ProductSpec, Trade},
+    trade::{ExerciseRecord, FallbackSelection, FallbackTerms, Fee, ProductSpec, Trade},
 };
+use rhodium_testkit::assert_close;
+use support::{CASH, refdata, refdata_parts};
 
+// 契約だけを持ち、行使・手数料・フォールバックの指定がない取引。
 fn trade(contract: Contract) -> Result<Trade> {
     Ok(Trade {
         id: "T-001".parse()?,
@@ -28,20 +50,34 @@ fn trade(contract: Contract) -> Result<Trade> {
     })
 }
 
-#[expect(
-    clippy::wildcard_enum_match_arm,
-    reason = "#196 で各 variant を網羅する"
-)]
+/// 受渡を契約の順に並べ、どれも定数に畳まれていることを確かめて額を返す。
 fn amounts(product: &Product) -> Result<Vec<f64>> {
     cashflows(product)
         .iter()
-        .map(|payment| match payment.amount.expression() {
-            analysis::Expression::Constant(value) => Ok(value),
-            _ => Err(Error::Internal("expected a resolved amount".into())),
+        .map(|payment| {
+            let Expression::Constant(value) = payment.amount.expression() else {
+                return Err(Error::Internal("expected a resolved amount".into()));
+            };
+            Ok(value)
         })
         .collect()
 }
 
+/// 受渡がちょうど一つで分岐の外にあり、額が定数に畳まれていることを確かめて額を返す。
+fn const_amount(product: &Product) -> Result<f64> {
+    let flows = cashflows(product);
+    let [flow] = flows.as_slice() else {
+        return Err(Error::Internal("expected one payment".into()));
+    };
+    if !flow.branches.is_empty() {
+        return Err(Error::Internal("expected an unconditional payment".into()));
+    }
+    let Expression::Constant(value) = flow.amount.expression() else {
+        return Err(Error::Internal("expected a fixed payment".into()));
+    };
+    Ok(value)
+}
+
 #[test]
 fn partial_terminations_scale_final_amounts_at_the_cutoff_and_keep_fees() -> Result<()> {
     // 変動クーポンに上限をかけた後で解約する。元本を縮めてから上限をかけると 50 になる。
@@ -182,6 +218,8 @@ fn replacement_applies_pre_cutoff_exercise_before_splitting() -> Result<()> {
     Ok(())
 }
 
+// 当日公表のターム金利 JPY-TERM を足し、TONA の複利に 10bp を足す代替を定めた参照データ。
+// 代替は 2021 年版の定義（2027-01-05 以降の約定）、プロトコル p、個別の合意のどれかで適用する。
 fn fallback_data(status: IndexStatus) -> Result<(RefData, IndexId, IndexId)> {
     let mut parts = refdata_parts()?;
     let term: IndexId = "JPY-TERM".parse()?;
@@ -236,6 +274,7 @@ fn fallback_data(status: IndexStatus) -> Result<(RefData, IndexId, IndexId)> {
     Ok((RefData::new(parts)?, term, overnight))
 }
 
+// 指標の一つの観測をそのまま 2027-01-12 に払う取引。
 fn fixing_trade(index: &IndexId, date: Date) -> Result<Trade> {
     let payment = "2027-01-12T08:00:00Z".parse()?;
     trade(Contract::build(|b| {
diff --git a/crates/rhodium-products/src/ir/tests/standard.rs b/crates/rhodium-products/tests/standard_products.rs
similarity index 95%
rename from crates/rhodium-products/src/ir/tests/standard.rs
rename to crates/rhodium-products/tests/standard_products.rs
index 8463a1a57dddbb26c10260216af7b1e814ca876d..446ccf28a1d2ddc1fd11845158c49642a8b76749 100644
--- a/crates/rhodium-products/src/ir/tests/standard.rs
+++ b/crates/rhodium-products/tests/standard_products.rs
@@ -1,3 +1,5 @@
+//! 預金、OIS とスワップの標準商品を IR に下ろし、簡約後の受渡を手計算と照合する。
+
 use std::num::NonZeroU8;
 
 use rhodium_core::{
@@ -30,15 +32,20 @@ use rhodium_core::{
     tenor::Frequency,
     time::{Instant, LocalEventTime, UtcOffset, ZoneRules},
 };
+use rhodium_products::{
+    ir::{
+        Product, SettlementBoundary,
+        analysis::{Expression, cashflows},
+        reduce,
+    },
+    standard::{
+        Coupon, Deposit, DepositDates, InstrumentTemplate, Leg, NotionalExchange, PayReceive,
+        PaymentRounding, Position, ProductKind, StandardProduct, Swap,
+    },
+};
 use rhodium_testkit::assert_close;
 
-use crate::standard::{
-    Coupon, Deposit, DepositDates, InstrumentTemplate, Leg, NotionalExchange, PayReceive,
-    PaymentRounding, Position, ProductKind, StandardProduct, Swap,
-};
-
-use super::{ContractNode, ExprNode, Product, SettlementBoundary, reduce};
-
+// 受渡額を各通貨の最小単位へ切り捨てる参照データ。
 fn refdata() -> Result<RefData> {
     refdata_rounding(RoundingMode::TowardZero)
 }
@@ -212,11 +219,8 @@ fn fixings(refdata: &RefData, index: &str, rows: &[(&str, f64)]) -> Result<Fixin
 }
 
 /// `reduce` 後の受渡を (時刻, 確定した額) の列にする。未確定の額は `None`。
-/// `give` は構築時に額の符号反転になるので、契約は `All` と `Pay` だけからなる。
-#[expect(
-    clippy::wildcard_enum_match_arm,
-    reason = "#196 で各 variant を網羅する"
-)]
+/// `give` は構築時に額の符号反転になるので、標準商品の受渡は条件や選択の枝の中に現れない。
+/// 共有した支払ノードは出現ごとに数える（`cashflows` の規約）。
 fn reduced_payments(
     product: &Product,
     as_of: &str,
@@ -230,23 +234,20 @@ fn reduced_payments(
         &[],
         SettlementBoundary::Include,
     )?;
-    let mut work = vec![reduced.contract.root];
-    let mut payments = Vec::new();
-    while let Some(node) = work.pop() {
-        match &reduced.contract.contracts[node as usize] {
-            ContractNode::All(children) => work.extend(children),
-            ContractNode::Pay { asset, amount, at } => {
-                assert_eq!(asset, &Asset::Cash(ccy));
-                let value = match reduced.contract.exprs[*amount as usize] {
-                    ExprNode::Const(value) => Some(value.value()),
-                    _ => None,
-                };
-                payments.push((*at, value));
+    cashflows(&reduced)
+        .into_iter()
+        .map(|payment| {
+            assert_eq!(payment.asset, &Asset::Cash(ccy));
+            if !payment.branches.is_empty() {
+                return Err(Error::Internal("expected a linear cashflow".into()));
             }
-            _ => return Err(Error::Internal("expected a linear cashflow".into())),
-        }
-    }
-    Ok(payments)
+            // 額が定数に畳まれていれば確定、ほかの式なら未確定。
+            let Expression::Constant(value) = payment.amount.expression() else {
+                return Ok((payment.at, None));
+            };
+            Ok((payment.at, Some(value)))
+        })
+        .collect()
 }
 
 fn assert_payments(
@@ -312,6 +313,7 @@ fn deposit_short_tenors_cross_the_weekend_and_reverse_for_borrowing() -> Result<
     let refdata = refdata()?;
     let fixings = fixings(&refdata, "JPY-TONA", &[])?;
     // ON、TN、SN はそれぞれ当日、翌営業日、スポット日から金曜〜月曜の 3 日を数える。
+    // 借入（Short）なので、元本を受け取って元利を返す向きになる。
     for (trade, tenor) in [
         ("2027-01-08", "ON"),
         ("2027-01-07", "TN"),
@@ -360,6 +362,7 @@ fn confirmed_deposit_keeps_the_holiday_and_retains_its_warning_after_reduction()
     )?
     .to_product(start, &refdata)?;
     let as_of = "2027-01-08T00:00:00Z";
+    // 確定した日付は休日でも動かさず、その旨の警告を IR の構築時にも簡約の後にも返す。
     let expected = [Warning::ExplicitDateIsHoliday { date: start }];
     assert_eq!(product.warnings(), &expected);
     assert_eq!(
diff --git a/crates/rhodium-products/tests/support/mod.rs b/crates/rhodium-products/tests/support/mod.rs
new file mode 100644
index 0000000000000000000000000000000000000000..3483ce2eb0b17cb9a10520b11645b9b3de27ffb2
--- /dev/null
+++ b/crates/rhodium-products/tests/support/mod.rs
@@ -0,0 +1,178 @@
+//! 複数の統合テストが共有する参照データと定数。
+//!
+//! `src/ir/tests.rs` の単体テストも同じ参照データを持つ。単体テストは `tests/` の補助を
+//! `#[path]` なしで読めず、統合テストは `#[cfg(test)]` の補助を読めないので、両方に置く。
+
+use rhodium_core::{
+    Result,
+    calendar::{BaseCalendar, CalendarSpec, Extrapolation, WeekendRule},
+    currency::{Asset, Currency},
+    date::Weekday,
+    day_count::DayCount,
+    fx::CcyPair,
+    id::{CommodityId, CreditKey, IndexId},
+    refdata::{
+        CommodityDefinition, CurrencyDefinition, DailyPublication, Deliverability, FxConvention,
+        FxSpotLag, IndexDefinition, IndexKind, IndexStatus, PublicationRule, RefData, RefDataParts,
+    },
+    rounding::{Rounding, RoundingMode},
+    time::{LocalEventTime, UtcOffset, ZoneRules},
+};
+
+/// 円の現金。契約の受渡の既定の資産。
+pub const CASH: Asset = Asset::Cash(Currency::JPY);
+
+fn credit_key() -> Result<CreditKey> {
+    Ok(CreditKey {
+        entity: "EXAMPLE".parse()?,
+        seniority: "SENIOR".parse()?,
+        doc_clause: "CR".parse()?,
+        ccy: Currency::JPY,
+    })
+}
+
+/// [`refdata_parts`] を検証した参照データ。
+///
+/// # Errors
+///
+/// 参照データの組み立てか検証に失敗した場合。
+pub fn refdata() -> Result<RefData> {
+    RefData::new(refdata_parts()?)
+}
+
+/// USDJPY の受渡慣行（T+2、pip 0.01 円）を、通貨と同じ決済時刻で足す。
+fn add_usd_jpy(parts: &mut RefDataParts, clock: &LocalEventTime) -> Result<()> {
+    let usd_jpy: CcyPair = "USDJPY".parse()?;
+    parts.fx_conventions.insert(
+        usd_jpy,
+        FxConvention {
+            pair: usd_jpy,
+            spot_lag: FxSpotLag::TwoBusinessDays,
+            settlement: clock.clone(),
+            pip_exponent: 2,
+        },
+    );
+    Ok(())
+}
+
+/// 2027 年の JPTO の暦（土日休み）と UTC+9 の東京、円とドル、JPY-TONA・USDJPY-SPOT・
+/// 信用事由・回収率の指標、WTI を持つ参照データの部品。
+///
+/// フォールバックのテストは、この部品に指標と代替を足してから検証する。
+///
+/// # Errors
+///
+/// 識別子や日時の解析か、暦と地域の構築に失敗した場合。
+pub fn refdata_parts() -> Result<RefDataParts> {
+    let mut parts = RefDataParts::default();
+    let calendar = "JPTO".parse()?;
+    let clock = LocalEventTime {
+        time: "09:00:00".parse()?,
+        zone: "Asia/Tokyo".parse()?,
+    };
+    parts.calendars.insert(
+        calendar,
+        BaseCalendar::new(
+            "2027-01-01".parse()?..="2027-12-31".parse()?,
+            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
+            Vec::new(),
+            Extrapolation::Strict,
+        )?,
+    );
+    parts.zones.insert(
+        clock.zone.clone(),
+        ZoneRules::new(
+            "2027-01-01T00:00:00Z".parse()?..="2027-12-31T23:59:59Z".parse()?,
+            UtcOffset::from_seconds(9 * 3600)?,
+            Vec::new(),
+        )?,
+    );
+    for code in [Currency::JPY, Currency::USD] {
+        parts.currencies.insert(
+            code,
+            CurrencyDefinition {
+                code,
+                minor_units: u8::from(code == Currency::USD) * 2,
+                payment_rounding: Rounding::new(1.0, RoundingMode::HalfEven)?,
+                deliverability: Deliverability::Deliverable,
+                settlement_time: clock.clone(),
+                settlement_calendar: "JPTO".parse()?,
+            },
+        );
+    }
+    let calendar = CalendarSpec::Base("JPTO".parse()?);
+    let daily = |business_days| -> Result<_> {
+        Ok(PublicationRule::Daily(DailyPublication {
+            business_days,
+            calendar: calendar.clone(),
+            time: "10:00:00".parse()?,
+            zone: clock.zone.clone(),
+        }))
+    };
+    // USDJPY の FX fixing の受渡日は、ペアの慣行のスポットラグで決まる。
+    add_usd_jpy(&mut parts, &clock)?;
+    for (name, kind, publication, observation_time) in [
+        (
+            "JPY-TONA",
+            IndexKind::OvernightRate {
+                ccy: Currency::JPY,
+                day_count: DayCount::Act365F,
+            },
+            daily(1)?,
+            Some(clock.clone()),
+        ),
+        (
+            "USDJPY-SPOT",
+            IndexKind::FxFixing {
+                pair: "USDJPY".parse()?,
+            },
+            PublicationRule::Intraday { delay_micros: 0 },
+            None,
+        ),
+        (
+            "JP-EVENT",
+            IndexKind::CreditEvent {
+                credit: credit_key()?,
+            },
+            daily(0)?,
+            None,
+        ),
+        (
+            "JP-RECOVERY",
+            IndexKind::Recovery {
+                credit: credit_key()?,
+            },
+            daily(1)?,
+            None,
+        ),
+    ] {
+        let id: IndexId = name.parse()?;
+        parts.indices.insert(
+            id.clone(),
+            IndexDefinition {
+                id,
+                kind,
+                calendar: calendar.clone(),
+                publication,
+                observation_time,
+                status: IndexStatus::Active,
+            },
+        );
+    }
+    insert_commodity(&mut parts)?;
+    Ok(parts)
+}
+
+// 現物受渡の例が使うコモディティ。
+fn insert_commodity(parts: &mut RefDataParts) -> Result<()> {
+    let oil: CommodityId = "WTI-CUSHING".parse()?;
+    parts.commodities.insert(
+        oil.clone(),
+        CommodityDefinition {
+            id: oil,
+            unit: "bbl".into(),
+            currency: Currency::USD,
+        },
+    );
+    Ok(())
+}
diff --git a/crates/rhodium-risk/Cargo.toml b/crates/rhodium-risk/Cargo.toml
index a8b430aaaffcec4b08c51406cf0071c4bd4e9839..e8cc216309dc58f625bd0939c9dc908e44d5cbe9 100644
--- a/crates/rhodium-risk/Cargo.toml
+++ b/crates/rhodium-risk/Cargo.toml
@@ -18,6 +18,7 @@ scioxide-ad.workspace = true
 scioxide-core.workspace = true
 
 [dev-dependencies]
+rhodium-fixtures.workspace = true
 proptest.workspace = true
 rhodium-products.workspace = true
 rhodium-testkit.workspace = true
diff --git a/crates/rhodium-risk/tests/least_squares.rs b/crates/rhodium-risk/tests/least_squares.rs
index 7c790fd161903b240aa3491dd0aa8a168d388d62..d3762c060607a0156ebf0ead9c9ec76d4d414469 100644
--- a/crates/rhodium-risk/tests/least_squares.rs
+++ b/crates/rhodium-risk/tests/least_squares.rs
@@ -1,5 +1,7 @@
 //! 最小二乗の校正グループを通る quote デルタを、手計算と再校正の中心差分で照合する。
 
+pub mod support;
+
 use std::collections::BTreeMap;
 
 use rhodium_core::{
@@ -8,6 +10,8 @@ use rhodium_core::{
     date::Date,
     money::Money,
 };
+use rhodium_fixtures::{COLLATERAL, bumped, solved};
+use rhodium_fixtures::{ois, steps};
 use rhodium_market::{
     curve::{CurveId, CurveRole, CurveSpec, Extrapolation},
     quote::{QuoteId, QuoteKind, QuoteSet},
@@ -20,54 +24,19 @@ use rhodium_pricing::{
         CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
     },
     prepare::{PrepareInput, Prepared, evaluate, prepare},
-    resolve::CollateralKey,
 };
 use rhodium_products::{
     ir::{Contract, Product},
     standard::{Deposit, DepositDates, PaymentRounding, Position, StandardProduct},
 };
-use rhodium_risk::{delta::quote_delta, gradient::param_gradient, vector::RiskVector};
-use rhodium_testkit::{assert_close, check_aad_vs_bump};
-use scioxide_core::ExecutionContext;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、最小二乗の照合に使わない関数がある"
-)]
-#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
-pub mod ois;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、最小二乗の照合に使わない関数がある"
-)]
-#[path = "../../rhodium-pricing/tests/fixtures/steps.rs"]
-pub mod steps;
+use rhodium_risk::{delta::quote_delta, gradient::param_gradient};
+use rhodium_testkit::{SEQUENTIAL, assert_close, check_aad_vs_bump, entry_or_zero};
+use support::{BP, vector_for};
 
 use ois::{Inputs, Node, cases, midnight, prepare_ois};
 use steps::{PenalizedTona, TURN_QUOTE};
 
-const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
-const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
 const NOTIONAL: f64 = 1_000_000.0;
-const BP: f64 = 1e-4;
-
-fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
-    vector
-        .entries()
-        .binary_search_by(|(key, _)| key.cmp(id))
-        .ok()
-        .map_or(0.0, |index| vector.entries()[index].1)
-}
-
-fn jpy(report: &rhodium_risk::delta::DeltaReport) -> Result<&RiskVector<QuoteId>> {
-    report
-        .by_ccy()
-        .iter()
-        .find(|vector| vector.ccy() == Currency::JPY)
-        .ok_or_else(|| Error::Internal("no JPY quote deltas".into()))
-}
-
 #[test]
 fn lsq_adjoint_oracle() -> Result<()> {
     // 感応度の検証表の例: min (θ − q₁)² + 4(θ − q₂)² の解 θ = (q₁ + 4q₂)/5 で、
@@ -126,16 +95,41 @@ fn lsq_adjoint_oracle() -> Result<()> {
     let pv = evaluate(&trade, market.market())?.get(Currency::JPY);
     assert_close!(theta_gradient.1, pv, 1e-12, 1e-6);
     let report = quote_delta(std::slice::from_ref(&trade), &market, SEQUENTIAL)?;
-    let vector = jpy(&report)?;
+    let vector = vector_for(report.by_ccy(), Currency::JPY)?;
     // マークの感応度の単位は自然な単位 1。
-    assert_close!(sensitivity(vector, &marks[0]), pv / 5.0, 1e-10, 1e-6);
-    assert_close!(sensitivity(vector, &marks[1]), 4.0 * pv / 5.0, 1e-10, 1e-6);
+    assert_close!(
+        entry_or_zero(vector.entries(), &marks[0]),
+        pv / 5.0,
+        1e-10,
+        1e-6
+    );
+    assert_close!(
+        entry_or_zero(vector.entries(), &marks[1]),
+        4.0 * pv / 5.0,
+        1e-10,
+        1e-6
+    );
     Ok(())
 }
 
 /// バンプの再校正の許容誤差。最小二乗のグループは停留条件の残差を θ に戻した各行の動きで判定する。
 const BUMP_TOLERANCE: f64 = 1e-15;
 
+/// 残差のノルム `tolerance` で校正し、すべてのグループが解けた市場を返す。
+///
+/// 最小二乗のグループは停留条件で止まり、quote の残差が残るので、残差の大きさは確かめない
+/// （共有の `support::calibrated` は正方系の残差も確かめる）。
+fn calibrated(
+    definition: &MarketDefinition,
+    quotes: &QuoteSet,
+    ctx: &Context<'_>,
+    tolerance: f64,
+) -> Result<CalibratedMarket> {
+    let market = calibrate_with_residual_tolerance(definition, quotes, ctx, tolerance)?;
+    solved(&market)?;
+    Ok(market)
+}
+
 /// 校正の停止の誤差を一つの PV に伝えた幅 `Σ_j |∂V/∂θ_j|·|δθ_j|`（品質の「ソルバーの誤差」）。
 ///
 /// δθ = (∂G/∂θ)⁻¹·G は各グループの `theta_error`（列は θ の添字の昇順）。最小二乗の罰則だけで
@@ -183,29 +177,6 @@ fn solver_pv_error(trade: &Prepared, market: &CalibratedMarket) -> Result<f64> {
         .sum()
 }
 
-fn calibrated(
-    definition: &MarketDefinition,
-    quotes: &QuoteSet,
-    ctx: &Context<'_>,
-    tolerance: f64,
-) -> Result<CalibratedMarket> {
-    let market = calibrate_with_residual_tolerance(definition, quotes, ctx, tolerance)?;
-    for group in market.groups() {
-        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
-            return Err(Error::Calibration(format!("{group:?}")));
-        }
-    }
-    Ok(market)
-}
-
-fn bumped(quotes: &QuoteSet, id: &QuoteId, bump: f64) -> Result<QuoteSet> {
-    QuoteSet::new(
-        quotes
-            .iter()
-            .map(|(key, value)| (key.clone(), if key == id { value + bump } else { value })),
-    )
-}
-
 /// 階段の区間をまたぐ 9M と、接続点をまたぐ 18M の OIS（quote から外した固定金利）と、
 /// 年末をまたぐ預金。
 fn trades(
@@ -281,7 +252,7 @@ fn penalized_meeting_steps_quote_delta_matches_bump() -> Result<()> {
     let prepared = trades(&tona, &base, &ctx)?;
     for (index, trade) in prepared.iter().enumerate() {
         let report = quote_delta(std::slice::from_ref(trade), &base, SEQUENTIAL)?;
-        let vector = jpy(&report)?;
+        let vector = vector_for(report.by_ccy(), Currency::JPY)?;
         for id in &ids {
             // 中心差分が使う ±BP と ±BP/2 の再校正ごとに PV と、停止の誤差を θ 勾配で伝えた幅を
             // 求め、大きい方の幅を一つの PV の誤差にする（品質の c_s = 2 は照合の関数が掛ける）。
@@ -311,7 +282,7 @@ fn penalized_meeting_steps_quote_delta_matches_bump() -> Result<()> {
                         .map(|&(_, pv, _)| pv)
                         .ok_or_else(|| Error::Internal(format!("unexpected bump {h}")))
                 },
-                sensitivity(vector, id),
+                entry_or_zero(vector.entries(), id),
                 BP,
                 pv_error,
                 1e-6,
@@ -328,8 +299,11 @@ fn penalized_meeting_steps_quote_delta_matches_bump() -> Result<()> {
     let deviation = ids
         .iter()
         .map(|id| {
-            let lhs = sensitivity(jpy(&exact)?, id);
-            let rhs = sensitivity(jpy(&gauss_newton)?, id);
+            let lhs = entry_or_zero(vector_for(exact.by_ccy(), Currency::JPY)?.entries(), id);
+            let rhs = entry_or_zero(
+                vector_for(gauss_newton.by_ccy(), Currency::JPY)?.entries(),
+                id,
+            );
             Ok::<_, Error>((lhs - rhs).abs() / lhs.abs().max(1e-6))
         })
         .collect::<Result<Vec<_>>>()?
diff --git a/crates/rhodium-risk/tests/m1_quote_delta.rs b/crates/rhodium-risk/tests/m1_quote_delta.rs
index 3b76c7593d2a7bb845ac472742bc63fd2efd0b4b..d25364b675dd6f8d05fabc80d784a677cf6eb032 100644
--- a/crates/rhodium-risk/tests/m1_quote_delta.rs
+++ b/crates/rhodium-risk/tests/m1_quote_delta.rs
@@ -9,11 +9,9 @@
 //!
 //! 再校正した市場は quote とバンプの組ごとに一度だけ作り、同じ市場で全取引を評価する。
 
-use std::{
-    cell::RefCell,
-    collections::HashMap,
-    num::{NonZeroU8, NonZeroU32},
-};
+pub mod support;
+
+use std::num::{NonZeroU8, NonZeroU32};
 
 use rhodium_core::{
     Error, Result,
@@ -32,85 +30,30 @@ use rhodium_core::{
         ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination, ScheduleSpec,
     },
 };
+use rhodium_fixtures::COLLATERAL;
+use rhodium_fixtures::{futures, multicurve, ois, steps};
 use rhodium_market::quote::{QuoteId, QuoteKind, QuoteSet};
 use rhodium_pricing::{
-    calibration::{CalibratedMarket, GroupOutcome, calibrate_with_residual_tolerance},
+    calibration::CalibratedMarket,
     context::Context,
     definition::MarketDefinition,
     prepare::{PrepareInput, Prepared, evaluate, prepare},
-    resolve::CollateralKey,
 };
 use rhodium_products::standard::{
     Coupon, Deposit, DepositDates, Fra, Leg, NotionalExchange, PayReceive, PaymentRounding,
     Position, RateFuture, StandardProduct, Swap,
 };
-use rhodium_risk::{delta::quote_delta, vector::RiskVector};
-use rhodium_testkit::check_aad_vs_bump;
-use scioxide_core::ExecutionContext;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
-)]
-#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
-pub mod ois;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
-)]
-#[path = "../../rhodium-pricing/tests/fixtures/multicurve.rs"]
-pub mod multicurve;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
-)]
-#[path = "../../rhodium-pricing/tests/fixtures/steps.rs"]
-pub mod steps;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
-)]
-#[path = "../../rhodium-pricing/tests/fixtures/futures.rs"]
-pub mod futures;
+use rhodium_risk::delta::quote_delta;
+use rhodium_testkit::{SEQUENTIAL, check_aad_vs_bump, entry_or_zero};
+use support::{BP, BUMP_TOLERANCE, Recalibrations, calibrated};
 
 use multicurve::{MulticurveCase, all_multicurve_cases};
 use ois::{Inputs, Node, prepare_ois};
 use steps::PenalizedTona;
 
-const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
-const BP: f64 = 1e-4;
-const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
-/// バンプの再校正の残差の許容誤差（通常の 1e-14 の 1/100）。
-const BUMP_TOLERANCE: f64 = 1e-16;
-
 /// 市場から照合する取引を準備する関数。取引の名前と準備の組を返す。
 type Trades<'a> = dyn Fn(&CalibratedMarket) -> Result<Vec<(String, Prepared)>> + 'a;
 
-fn calibrated(
-    definition: &MarketDefinition,
-    quotes: &QuoteSet,
-    ctx: &Context<'_>,
-) -> Result<CalibratedMarket> {
-    let market = calibrate_with_residual_tolerance(definition, quotes, ctx, BUMP_TOLERANCE)?;
-    for group in market.groups() {
-        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
-            return Err(Error::Calibration(format!("{group:?}")));
-        }
-    }
-    Ok(market)
-}
-
-fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
-    vector
-        .entries()
-        .binary_search_by(|(key, _)| key.cmp(id))
-        .ok()
-        .map_or(0.0, |index| vector.entries()[index].1)
-}
-
 /// quote の報告単位（デルタの値の単位）と、バンプの幅。マークは報告単位 1 だが、1 だけ動かすと
 /// 非線形の項が大きいので、1bp の幅で照合し、AAD を幅あたりに換算する。
 #[expect(
@@ -138,23 +81,13 @@ fn check_all(
     floor: f64,
     trades: &Trades<'_>,
 ) -> Result<()> {
-    let base = calibrated(definition, quotes, ctx)?;
-    // 同じ (quote, 幅) の再校正は取引の間で共有する。
-    let cache: RefCell<HashMap<(QuoteId, u64), CalibratedMarket>> = RefCell::new(HashMap::new());
+    let base = calibrated(definition, quotes, ctx, BUMP_TOLERANCE)?;
+    // 同じ (quote, 幅) の再校正と、その市場での全取引の準備は、取引の間で共有する。
+    let recalibrations = Recalibrations::new(definition, quotes, ctx, BUMP_TOLERANCE, trades);
     let reprice = |id: &QuoteId, bump: f64, index: usize| -> Result<f64> {
-        let key = (id.clone(), bump.to_bits());
-        if !cache.borrow().contains_key(&key) {
-            let bumped =
-                QuoteSet::new(quotes.iter().map(|(key, value)| {
-                    (key.clone(), if key == id { value + bump } else { value })
-                }))?;
-            let market = calibrated(definition, &bumped, ctx)?;
-            cache.borrow_mut().insert(key.clone(), market);
-        }
-        let cache = cache.borrow();
-        let market = &cache[&key];
-        let (_, prepared) = &trades(market)?[index];
-        Ok(evaluate(prepared, market.market())?.get(ccy))
+        let entry = recalibrations.get(id, bump)?;
+        let (market, prepared) = &*entry;
+        Ok(evaluate(&prepared[index].1, market.market())?.get(ccy))
     };
     for (index, (name, prepared)) in trades(&base)?.iter().enumerate() {
         let report = quote_delta(std::slice::from_ref(prepared), &base, SEQUENTIAL)?;
@@ -175,7 +108,7 @@ fn check_all(
             let (report_unit, width) = units(&base, id)?;
             check_aad_vs_bump(
                 |h| reprice(id, h, index),
-                sensitivity(vector, id) * width / report_unit,
+                entry_or_zero(vector.entries(), id) * width / report_unit,
                 width,
                 BUMP_TOLERANCE * scale,
                 floor,
diff --git a/crates/rhodium-risk/tests/quote_delta.rs b/crates/rhodium-risk/tests/quote_delta.rs
index 0aaa049c09077fb5fde9ae0cb07889335d5e675a..600229893d2383577e8fded135814bcbc5918f56 100644
--- a/crates/rhodium-risk/tests/quote_delta.rs
+++ b/crates/rhodium-risk/tests/quote_delta.rs
@@ -1,5 +1,7 @@
 //! M0 の市場と商品で、quote デルタを全体再校正の中心差分と照合する。
 
+pub mod support;
+
 use std::{collections::BTreeSet, num::NonZeroUsize, sync::Arc};
 
 use proptest::{
@@ -20,6 +22,8 @@ use rhodium_core::{
         ScheduleSpec,
     },
 };
+use rhodium_fixtures::{COLLATERAL, bumped, prepare_standard};
+use rhodium_fixtures::{multicurve, ois};
 use rhodium_market::{
     curve::{CurveId, Extrapolation},
     fx::FxCollateral,
@@ -27,9 +31,9 @@ use rhodium_market::{
     quote::{QuoteId, QuoteSet},
 };
 use rhodium_pricing::{
-    calibration::{CalibratedMarket, GroupOutcome, calibrate, calibrate_with_residual_tolerance},
+    calibration::{CalibratedMarket, GroupOutcome, calibrate},
     context::Context,
-    definition::{DirectFx, FxDefinition, FxImpliedCurve, FxQuotes, MarketDefinition},
+    definition::{DirectFx, FxDefinition, FxImpliedCurve, FxQuotes},
     prepare::{PrepareInput, Prepared, evaluate, prepare},
     resolve::CollateralKey,
 };
@@ -44,33 +48,19 @@ use rhodium_products::{
 use rhodium_risk::{
     delta::{DeltaReport, quote_delta, quote_delta_in},
     gradient::param_gradient,
-    vector::RiskVector,
 };
-use rhodium_testkit::{TickClock, assert_close, check_aad_vs_bump};
+use rhodium_testkit::{
+    SEQUENTIAL, TickClock, assert_close, check_aad_vs_bump, entry_or_zero, minor_unit_abs,
+};
 use scioxide_core::{ClockBinding, Deadline};
-
-#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
-pub mod ois;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、このテストで使わない関数がある"
-)]
-#[path = "../../rhodium-pricing/tests/fixtures/multicurve.rs"]
-pub mod multicurve;
+use support::{BP, BUMP_TOLERANCE, Recalibrations, calibrated, vector_for};
 
 use ois::{Case, Inputs, cases, definition, fixture_quotes, prepare_ois};
 
-const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
-/// 一つの取引のデルタは逐次に求める。
-const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
-const BP: f64 = 1e-4;
 const NOTIONAL: f64 = 1_000_000.0;
 // バンプの再校正の残差の許容誤差。品質の照合手順どおり、通常の校正の許容誤差（1e-14）の 1/100。
 // 超過成長を対数から exp_m1 で求めるので、1M の OIS でも残差は約 1e-17 まで下がり、この値に届く。
 // この残差を δ_V に伝える。δ_V / h₂ は 1 bp あたりの感応度の約 2e-12 倍で、比べる相対誤差 1e-10 より小さい。
-const BUMP_TOLERANCE: f64 = 1e-16;
-
 /// 両通貨の現金担保の割引を、それぞれの OIS 曲線に結ぶ価格設定。
 const CASH_DISCOUNTING: &str = r#"
 [[discounting]]
@@ -84,25 +74,6 @@ collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
 curve = "USD-SOFR"
 "#;
 
-fn prepare_trade(
-    product: &StandardProduct,
-    date: Date,
-    collateral: &CollateralKey,
-    cm: &CalibratedMarket,
-    ctx: &Context<'_>,
-) -> Result<Prepared> {
-    prepare(
-        PrepareInput::Standard {
-            product,
-            trade_date: date,
-        },
-        collateral,
-        &[],
-        cm.market().structure(),
-        ctx,
-    )
-}
-
 fn cash(case: &Case) -> Result<CollateralKey> {
     Ok(CollateralKey::Cash(CashCollateralTerms {
         ccy: case.currency,
@@ -198,48 +169,8 @@ fn fixed_rfr_swap(case: &Case, ctx: &Context<'_>) -> Result<StandardProduct> {
     ])?))
 }
 
-fn calibrated(
-    def: &MarketDefinition,
-    quotes: &QuoteSet,
-    ctx: &Context<'_>,
-    tolerance: f64,
-) -> Result<CalibratedMarket> {
-    let cm = calibrate_with_residual_tolerance(def, quotes, ctx, tolerance)?;
-    for group in cm.groups() {
-        let GroupOutcome::Solved(solution) = group.outcome() else {
-            return Err(Error::Calibration(format!("risk test fixture: {group:?}")));
-        };
-        assert!(solution.residual_norm() <= tolerance);
-    }
-    Ok(cm)
-}
-
-fn bumped(quotes: &QuoteSet, id: &QuoteId, bump: f64) -> Result<QuoteSet> {
-    QuoteSet::new(
-        quotes
-            .iter()
-            .map(|(key, value)| (key.clone(), if key == id { value + bump } else { value })),
-    )
-}
-
-fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
-    vector
-        .entries()
-        .binary_search_by(|(key, _)| key.cmp(id))
-        .ok()
-        .map_or(0.0, |index| vector.entries()[index].1)
-}
-
-fn vector_for(vectors: &[RiskVector<QuoteId>], ccy: Currency) -> Result<&RiskVector<QuoteId>> {
-    vectors
-        .iter()
-        .find(|vector| vector.ccy() == ccy)
-        .ok_or_else(|| Error::Internal(format!("no quote deltas for {ccy}")))
-}
-
 fn absolute_tolerance(ctx: &Context<'_>, ccy: Currency) -> Result<f64> {
-    let digits = ctx.refdata.currency(ccy)?.minor_units;
-    Ok(1e-6 * 10.0_f64.powi(-i32::from(digits)))
+    Ok(minor_unit_abs(ctx.refdata.currency(ccy)?.minor_units))
 }
 
 /// バンプで照合する取引の名前、約定日、担保、商品。
@@ -303,52 +234,67 @@ fn quote_delta_vs_bump() -> Result<()> {
             .chain(fixture_quotes(&sofr)?),
     )?;
     let base = calibrate(&def, &quotes, &ctx)?;
-
-    for case in [&tona, &sofr] {
+    // 両通貨の取引を一つの並びにし、(quote, 幅) の組ごとの再校正の市場で全取引を一度だけ準備する。
+    // 取引ごと・quote ごとに解き直すと 2 通貨 × 5 取引 × 10 quote × 4 幅で約 400 回の全体の
+    // 再校正になるが、共有すれば 10 quote × 4 幅の 40 回で済む。
+    let trades = [&tona, &sofr]
+        .into_iter()
+        .map(|case| {
+            Ok(bump_trades(case, &ctx)?
+                .into_iter()
+                .map(move |trade| (case, trade)))
+        })
+        .collect::<Result<Vec<_>>>()?
+        .into_iter()
+        .flatten()
+        .collect::<Vec<_>>();
+    let prepare_all = |cm: &CalibratedMarket| -> Result<Vec<Prepared>> {
+        trades
+            .iter()
+            .map(|(_, (_, trade_date, collateral, product))| {
+                prepare_standard(product, *trade_date, collateral, cm.market(), &ctx)
+            })
+            .collect()
+    };
+    let recalibrations = Recalibrations::new(&def, &quotes, &ctx, BUMP_TOLERANCE, &prepare_all);
+    for (index, ((case, (name, ..)), plan)) in trades.iter().zip(prepare_all(&base)?).enumerate() {
         let absolute = absolute_tolerance(&ctx, case.currency)?;
         let block = case.block()?;
-        for (name, trade_date, collateral, product) in bump_trades(case, &ctx)? {
-            let prepare_one =
-                |cm: &CalibratedMarket| prepare_trade(&product, trade_date, &collateral, cm, &ctx);
-            let plan = prepare_one(&base)?;
-            let gradient = param_gradient(&plan, &base)?;
-            assert_close!(
-                gradient.pv().get(case.currency),
-                evaluate(&plan, base.market())?.get(case.currency),
-                1e-12,
-                absolute
-            );
-            let [active] = gradient.by_ccy() else {
-                return Err(Error::Internal("expected one currency in gradient".into()));
-            };
-            assert_ne!(active.entries(), []);
-            assert!(active.entries().iter().all(|(id, _)| id.block() == &block));
+        let gradient = param_gradient(&plan, &base)?;
+        assert_close!(
+            gradient.pv().get(case.currency),
+            evaluate(&plan, base.market())?.get(case.currency),
+            1e-12,
+            absolute
+        );
+        let [active] = gradient.by_ccy() else {
+            return Err(Error::Internal("expected one currency in gradient".into()));
+        };
+        assert_ne!(active.entries(), []);
+        assert!(active.entries().iter().all(|(id, _)| id.block() == &block));
 
-            let delta = quote_delta(&[plan], &base, SEQUENTIAL)?;
-            assert!(delta.failures().is_empty());
-            let vector = vector_for(delta.by_ccy(), case.currency)?;
-            assert_ne!(vector.entries(), []);
-            let scale: f64 = vector
-                .entries()
-                .iter()
-                .map(|(_, value)| value.abs() / BP)
-                .sum();
-            for (id, _) in quotes.iter() {
-                check_aad_vs_bump(
-                    |bump| {
-                        let cm =
-                            calibrated(&def, &bumped(&quotes, id, bump)?, &ctx, BUMP_TOLERANCE)?;
-                        Ok::<_, Error>(
-                            evaluate(&prepare_one(&cm)?, cm.market())?.get(case.currency),
-                        )
-                    },
-                    sensitivity(vector, id),
-                    BP,
-                    BUMP_TOLERANCE * scale,
-                    absolute,
-                    &format!("{name} {} {id}", case.currency),
-                )?;
-            }
+        let delta = quote_delta(&[plan], &base, SEQUENTIAL)?;
+        assert!(delta.failures().is_empty());
+        let vector = vector_for(delta.by_ccy(), case.currency)?;
+        assert_ne!(vector.entries(), []);
+        let scale: f64 = vector
+            .entries()
+            .iter()
+            .map(|(_, value)| value.abs() / BP)
+            .sum();
+        for (id, _) in quotes.iter() {
+            check_aad_vs_bump(
+                |bump| {
+                    let entry = recalibrations.get(id, bump)?;
+                    let (cm, prepared) = &*entry;
+                    Ok::<_, Error>(evaluate(&prepared[index], cm.market())?.get(case.currency))
+                },
+                entry_or_zero(vector.entries(), id),
+                BP,
+                BUMP_TOLERANCE * scale,
+                absolute,
+                &format!("{name} {} {id}", case.currency),
+            )?;
         }
     }
     Ok(())
@@ -395,7 +341,7 @@ fn rounded_payments_report_the_delta_of_the_unrounded_amount() -> Result<()> {
         ctx.refdata,
     )?;
     let pv = |product: &StandardProduct, cm: &CalibratedMarket| -> Result<f64> {
-        let plan = prepare_trade(product, tona.as_of, &COLLATERAL, cm, &ctx)?;
+        let plan = prepare_standard(product, tona.as_of, &COLLATERAL, cm.market(), &ctx)?;
         Ok(evaluate(&plan, cm.market())?.get(tona.currency))
     };
     let discount = |cm: &CalibratedMarket| -> Result<f64> {
@@ -412,7 +358,7 @@ fn rounded_payments_report_the_delta_of_the_unrounded_amount() -> Result<()> {
     // 丸めは実際に働いている（基準点で 1 円未満の差がある）。
     assert!(difference != 0.0 && difference.abs() < 1.0, "{difference}");
 
-    let plan = prepare_trade(&rounded, tona.as_of, &COLLATERAL, &base, &ctx)?;
+    let plan = prepare_standard(&rounded, tona.as_of, &COLLATERAL, base.market(), &ctx)?;
     let report = quote_delta(&[plan], &base, SEQUENTIAL)?;
     assert!(report.failures().is_empty());
     let vector = vector_for(report.by_ccy(), tona.currency)?;
@@ -429,14 +375,14 @@ fn rounded_payments_report_the_delta_of_the_unrounded_amount() -> Result<()> {
         pv(&rounded, &cm)
     };
     assert_eq!(rounded_pv(BP)?, rounded_pv(-BP)?);
-    assert!(sensitivity(vector, &six_months).abs() > absolute);
+    assert!(entry_or_zero(vector.entries(), &six_months).abs() > absolute);
     for (id, _) in quotes.iter() {
         check_aad_vs_bump(
             |bump| {
                 let cm = calibrated(&def, &bumped(&quotes, id, bump)?, &ctx, BUMP_TOLERANCE)?;
                 Ok::<_, Error>(pv(&unrounded, &cm)? + difference * discount(&cm)?)
             },
-            sensitivity(vector, id),
+            entry_or_zero(vector.entries(), id),
             BP,
             BUMP_TOLERANCE * scale,
             absolute,
@@ -493,8 +439,9 @@ fn short_tenor_quote_deltas_match_bumps() -> Result<()> {
                 convention_swap(&case, None, "2W", rates[2], PayReceive::Pay, 0.0005)?,
             ),
         ] {
-            let prepare_one =
-                |cm: &CalibratedMarket| prepare_trade(&product, case.as_of, &COLLATERAL, cm, &ctx);
+            let prepare_one = |cm: &CalibratedMarket| {
+                prepare_standard(&product, case.as_of, &COLLATERAL, cm.market(), &ctx)
+            };
             let delta = quote_delta(&[prepare_one(&base)?], &base, SEQUENTIAL)?;
             assert!(delta.failures().is_empty());
             let vector = vector_for(delta.by_ccy(), case.currency)?;
@@ -513,7 +460,7 @@ fn short_tenor_quote_deltas_match_bumps() -> Result<()> {
                             evaluate(&prepare_one(&cm)?, cm.market())?.get(case.currency),
                         )
                     },
-                    sensitivity(vector, id),
+                    entry_or_zero(vector.entries(), id),
                     BP,
                     BUMP_TOLERANCE * scale,
                     absolute,
@@ -559,7 +506,9 @@ fn quote_delta_stops_instead_of_dropping_trades_when_the_deadline_passes() -> Re
     )?;
     let trades = bump_trades(&tona, &ctx)?
         .into_iter()
-        .map(|(_, date, collateral, product)| prepare_trade(&product, date, &collateral, &cm, &ctx))
+        .map(|(_, date, collateral, product)| {
+            prepare_standard(&product, date, &collateral, cm.market(), &ctx)
+        })
         .collect::<Result<Vec<_>>>()?;
     let full = quote_delta(&trades, &cm, SEQUENTIAL)?;
     assert!(full.failures().is_empty());
@@ -655,7 +604,7 @@ fn check_linearity(
                 1 => convention_swap(case, None, "1Y", rate, side, 0.0)?,
                 _ => convention_swap(case, None, "2Y", rate, side, 0.001)?,
             };
-            prepare_trade(&product, case.as_of, &COLLATERAL, cm, ctx)
+            prepare_standard(&product, case.as_of, &COLLATERAL, cm.market(), ctx)
         })
         .collect::<Result<Vec<_>>>()?;
     // 合計は 3 ワーカーで区間ごとの和を足し、単独の値は逐次で求める。
@@ -686,9 +635,15 @@ fn check_linearity(
             .map(|(id, _)| id)
             .collect();
         for id in ids {
-            let sum: f64 = parts.iter().map(|part| sensitivity(part, id)).sum();
-            let gross: f64 = parts.iter().map(|part| sensitivity(part, id).abs()).sum();
-            assert_close!(sensitivity(total, id), sum, 1e-10, floor, gross);
+            let sum: f64 = parts
+                .iter()
+                .map(|part| entry_or_zero(part.entries(), id))
+                .sum();
+            let gross: f64 = parts
+                .iter()
+                .map(|part| entry_or_zero(part.entries(), id).abs())
+                .sum();
+            assert_close!(entry_or_zero(total.entries(), id), sum, 1e-10, floor, gross);
         }
     }
     // 通貨の違う感応度は足せない。
@@ -769,7 +724,7 @@ fn quote_delta_propagates_to_upstream_quotes() -> Result<()> {
                 let cm = calibrated(&def, &bumped(&quotes, &id, bump)?, &ctx, BUMP_TOLERANCE)?;
                 Ok::<_, Error>(evaluate(&prepare_one(&cm)?, cm.market())?.get(Currency::JPY))
             },
-            sensitivity(vector, &id),
+            entry_or_zero(vector.entries(), &id),
             BP,
             BUMP_TOLERANCE * scale,
             absolute_tolerance(&ctx, Currency::JPY)?,
@@ -810,14 +765,14 @@ fn tibor_basis_quote_delta_matches_full_recalibration() -> Result<()> {
                     calibrated(&definition, &bumped(&quotes, id, h)?, &ctx, BUMP_TOLERANCE)?;
                 reprice(market.market())
             },
-            sensitivity(vector, id),
+            entry_or_zero(vector.entries(), id),
             BP,
             BUMP_TOLERANCE * scale,
             floor,
             id.as_str(),
         )?;
     }
-    assert_ne!(sensitivity(vector, &jpy.quote(node)?), 0.0);
+    assert_ne!(entry_or_zero(vector.entries(), &jpy.quote(node)?), 0.0);
 
     // TONA を動かしたとき、基準のベーシス θ を固定すると下流の再校正を経由する寄与が消える。
     let basis = base
@@ -836,7 +791,8 @@ fn tibor_basis_quote_delta_matches_full_recalibration() -> Result<()> {
             reprice(&market)
         };
         let direct_delta = direct_pv(BP / 2.0)? - direct_pv(-BP / 2.0)?;
-        max_indirect = max_indirect.max((sensitivity(vector, &id) - direct_delta).abs());
+        max_indirect =
+            max_indirect.max((entry_or_zero(vector.entries(), &id) - direct_delta).abs());
     }
     // 比較の単位は元本 1。丸めや中心差分の打ち切りを超える間接寄与がある。
     assert!(
@@ -893,7 +849,7 @@ fn translation_fx() -> Result<()> {
     let usd_swap = convention_swap(&sofr, None, "1Y", 0.045, PayReceive::Receive, 0.0)?;
     let prepare = |market: &CalibratedMarket| {
         [&jpy_swap, &usd_swap]
-            .map(|swap| prepare_trade(swap, tona.as_of, &COLLATERAL, market, &ctx))
+            .map(|swap| prepare_standard(swap, tona.as_of, &COLLATERAL, market.market(), &ctx))
     };
     let front_jpy = tona.quote(&tona.nodes[0])?;
     let usd_rate = sofr.quote(&sofr.nodes[2])?;
@@ -962,15 +918,19 @@ fn translation_fx() -> Result<()> {
                 let foreign_pv = evaluate(foreign_trade, base.market())?.get(foreign);
                 let sign = if foreign == usd { 1.0 } else { -1.0 };
                 assert_close!(
-                    sensitivity(vector, &spot),
+                    entry_or_zero(vector.entries(), &spot),
                     sign * 0.01 * x0 * foreign_pv,
                     1e-10,
                     floor
                 );
-                let frozen_fx = sensitivity(vector_for(local.by_ccy(), reporting)?, &front_jpy)
-                    + x0 * sensitivity(vector_for(local.by_ccy(), foreign)?, &front_jpy);
+                let frozen_fx =
+                    entry_or_zero(vector_for(local.by_ccy(), reporting)?.entries(), &front_jpy)
+                        + x0 * entry_or_zero(
+                            vector_for(local.by_ccy(), foreign)?.entries(),
+                            &front_jpy,
+                        );
                 // スポット受渡日までの JPY の割引を経由する項が、現地デルタの換算に加わる。
-                assert!((sensitivity(vector, &front_jpy) - frozen_fx).abs() > floor);
+                assert!((entry_or_zero(vector.entries(), &front_jpy) - frozen_fx).abs() > floor);
             }
         }
     }
@@ -1007,8 +967,8 @@ fn check_translation_bumps(
     let parallel = vector_for(concurrent.by_ccy(), reporting)?;
     for (id, _) in base.quotes().iter() {
         assert_close!(
-            sensitivity(vector, id),
-            sensitivity(parallel, id),
+            entry_or_zero(vector.entries(), id),
+            entry_or_zero(parallel.entries(), id),
             1e-10,
             floor
         );
@@ -1026,7 +986,7 @@ fn check_translation_bumps(
         })
         .sum::<f64>();
     for &(id, unit) in targets {
-        let aad = sensitivity(vector, id);
+        let aad = entry_or_zero(vector.entries(), id);
         assert!(aad.abs() > floor, "{reporting} {id}: {aad:e}");
         check_aad_vs_bump(
             |h| reprice(id, h),
@@ -1075,8 +1035,8 @@ fn quote_delta_reporting_requires_fx_for_foreign_pv() -> Result<()> {
         let rhs = vector_for(same.by_ccy(), Currency::USD)?;
         for (id, _) in quotes.iter() {
             assert_close!(
-                sensitivity(lhs, id),
-                sensitivity(rhs, id),
+                entry_or_zero(lhs.entries(), id),
+                entry_or_zero(rhs.entries(), id),
                 1e-10,
                 absolute_tolerance(&ctx, Currency::USD)?
             );
@@ -1133,7 +1093,7 @@ fn settlement_at_valuation_instant_has_pv_without_delta() -> Result<()> {
     let def = definition(&[(&tona, &COLLATERAL)])?;
     let cm = calibrate(&def, &QuoteSet::new(fixture_quotes(&tona)?)?, &ctx)?;
     let product = deposit(&tona, "2024-12-16".parse()?, tona.as_of, Position::Long)?;
-    let settled = prepare_trade(&product, tona.as_of, &COLLATERAL, &cm, &ctx)?;
+    let settled = prepare_standard(&product, tona.as_of, &COLLATERAL, cm.market(), &ctx)?;
     let report = quote_delta(&[settled], &cm, SEQUENTIAL)?;
     assert!(report.by_ccy().is_empty() && report.failures().is_empty());
 
@@ -1153,7 +1113,7 @@ curve = "JPY-TONA"
         config: &config,
         ..ctx
     };
-    let due_today = prepare_trade(&product, tona.as_of, &COLLATERAL, &cm, &ctx)?;
+    let due_today = prepare_standard(&product, tona.as_of, &COLLATERAL, cm.market(), &ctx)?;
     assert_eq!(due_today.blocks(), []);
     let gradient = param_gradient(&due_today, &cm)?;
     assert!(gradient.pv().get(Currency::JPY) > 0.0);
diff --git a/crates/rhodium-risk/tests/support/mod.rs b/crates/rhodium-risk/tests/support/mod.rs
new file mode 100644
index 0000000000000000000000000000000000000000..a59bab76bbead0dd4c4eb1e36a9e86e2ea36ada8
--- /dev/null
+++ b/crates/rhodium-risk/tests/support/mod.rs
@@ -0,0 +1,119 @@
+//! risk の統合テストが共有する補助。
+//!
+//! 各テストのファイルは `pub mod support;` で取り込む。補助を `pub` にしておけば、取り込んだファイルが
+//! 使わない補助で `dead_code` が出ない。複数の crate のテストが使う補助は `rhodium-fixtures` に、
+//! 数値の比較は `rhodium-testkit` に置く。
+
+use std::{cell::RefCell, collections::HashMap, rc::Rc};
+
+use rhodium_core::{Error, Result, currency::Currency};
+use rhodium_fixtures::{bumped, solved_within};
+use rhodium_market::quote::{QuoteId, QuoteSet};
+use rhodium_pricing::{
+    calibration::{CalibratedMarket, calibrate_with_residual_tolerance},
+    context::Context,
+    definition::MarketDefinition,
+};
+use rhodium_risk::vector::RiskVector;
+
+/// 1bp。金利と spread の quote の感応度の単位で、中心差分のバンプの幅 h₁ にも使う。
+pub const BP: f64 = 1e-4;
+
+/// バンプの再校正の残差の許容誤差。品質の照合手順どおり、通常の校正の許容誤差（1e-14）の 1/100。
+///
+/// 超過成長を対数から `exp_m1` で求めるので、1M の OIS でも残差は約 1e-17 まで下がり、この値に届く。
+/// この残差を `δ_V` に伝える。`δ_V / h₂` は 1 bp あたりの感応度の約 2e-12 倍で、比べる相対誤差 1e-10 より
+/// 小さい。倍精度の分解能でこの値に届かない校正（XCCY のパー、最小二乗）は、届く最小の値を
+/// それぞれのテストで定める。
+pub const BUMP_TOLERANCE: f64 = 1e-16;
+
+/// 残差のノルム `tolerance` で校正し、すべてのグループが残差 `tolerance` 以内で解けた市場を返す。
+///
+/// 正方系だけの市場に使う。最小二乗のグループは停留条件で止まり残差が残るので、照合するテストが
+/// [`rhodium_fixtures::solved`] で解けたことだけを確かめる。
+///
+/// # Errors
+///
+/// 校正が失敗した場合と、解けないグループがある場合。
+pub fn calibrated(
+    definition: &MarketDefinition,
+    quotes: &QuoteSet,
+    ctx: &Context<'_>,
+    tolerance: f64,
+) -> Result<CalibratedMarket> {
+    let market = calibrate_with_residual_tolerance(definition, quotes, ctx, tolerance)?;
+    solved_within(&market, tolerance)?;
+    Ok(market)
+}
+
+/// 通貨ごとの quote デルタの中から `ccy` の報告通貨のベクトルを選ぶ。
+///
+/// # Errors
+///
+/// `ccy` のベクトルが無い場合。
+pub fn vector_for(vectors: &[RiskVector<QuoteId>], ccy: Currency) -> Result<&RiskVector<QuoteId>> {
+    vectors
+        .iter()
+        .find(|vector| vector.ccy() == ccy)
+        .ok_or_else(|| Error::Internal(format!("no quote deltas for {ccy}")))
+}
+
+/// 一つの (quote, 幅) の組で再校正した市場と、その市場で準備した取引。
+pub type Recalibrated<T> = Rc<(CalibratedMarket, T)>;
+
+/// 中心差分の再校正を、quote とバンプの幅の組ごとに一度だけ行い、その市場で準備した取引と共有する。
+///
+/// `check_aad_vs_bump` は一つの quote に 4 つの幅（±h₁、±h₂）で再校正を求める。取引ごと・quote ごとに
+/// 解き直すと、再校正の回数が取引の数に比例し、取引ごとに全取引を準備すると準備が取引の数の二乗に
+/// なる。同じ組の再校正と準備を一度だけ作り、全取引の照合で使い回す。
+pub struct Recalibrations<'a, T> {
+    definition: &'a MarketDefinition,
+    quotes: &'a QuoteSet,
+    ctx: &'a Context<'a>,
+    tolerance: f64,
+    prepare: &'a dyn Fn(&CalibratedMarket) -> Result<T>,
+    // 幅は f64 のビット列で区別する。照合の幅は固定の値なので、同じ幅は同じビット列になる。
+    cache: RefCell<HashMap<(QuoteId, u64), Recalibrated<T>>>,
+}
+
+impl<'a, T> Recalibrations<'a, T> {
+    /// 残差のノルム `tolerance` で [`calibrated`] し、`prepare` で取引を準備する再校正の入れ物。
+    pub fn new(
+        definition: &'a MarketDefinition,
+        quotes: &'a QuoteSet,
+        ctx: &'a Context<'a>,
+        tolerance: f64,
+        prepare: &'a dyn Fn(&CalibratedMarket) -> Result<T>,
+    ) -> Self {
+        Self {
+            definition,
+            quotes,
+            ctx,
+            tolerance,
+            prepare,
+            cache: RefCell::new(HashMap::new()),
+        }
+    }
+
+    /// `id` の quote を `bump` だけ動かして再校正した市場と、その市場で準備した取引。
+    ///
+    /// # Errors
+    ///
+    /// 再校正か取引の準備が失敗した場合。
+    pub fn get(&self, id: &QuoteId, bump: f64) -> Result<Recalibrated<T>> {
+        let key = (id.clone(), bump.to_bits());
+        if let Some(entry) = self.cache.borrow().get(&key) {
+            return Ok(Rc::clone(entry));
+        }
+        let market = calibrated(
+            self.definition,
+            &bumped(self.quotes, id, bump)?,
+            self.ctx,
+            self.tolerance,
+        )?;
+        let prepared = (self.prepare)(&market)?;
+        let entry = Rc::new((market, prepared));
+        self.cache.borrow_mut().insert(key, Rc::clone(&entry));
+        Ok(entry)
+    }
+}
diff --git a/crates/rhodium-risk/tests/xccy_quote_delta.rs b/crates/rhodium-risk/tests/xccy_quote_delta.rs
index c8064e62f1740427645de7142a004f505b1c37d6..6f633ac0ee5c059311735a984852e4e8915fc85d 100644
--- a/crates/rhodium-risk/tests/xccy_quote_delta.rs
+++ b/crates/rhodium-risk/tests/xccy_quote_delta.rs
@@ -1,71 +1,27 @@
 //! R6 の網羅表の行: 外貨担保曲線（FX ポイントと `MtM` XCCY）を通る `MtM` XCCY、FX フォワード、
 //! FX スワップの quote デルタ（Δq と Δx）を、quote ごとに全体を再校正した中心差分と照合する。
 
+pub mod support;
+
 use rhodium_core::{Error, Result, currency::Currency, date::Date, fx::FxRate, money::Money};
-use rhodium_market::quote::{QuoteId, QuoteKind, QuoteSet};
+use rhodium_fixtures::xccy;
+use rhodium_market::quote::{QuoteId, QuoteKind};
 use rhodium_pricing::{
-    calibration::{CalibratedMarket, GroupOutcome, calibrate_with_residual_tolerance},
+    calibration::CalibratedMarket,
     context::Context,
-    definition::MarketDefinition,
     prepare::{PrepareInput, Prepared, evaluate, prepare},
 };
 use rhodium_products::standard::{FxForward, FxSwap, Position, StandardProduct};
-use rhodium_risk::{delta::quote_delta_in, vector::RiskVector};
-use rhodium_testkit::check_aad_vs_bump;
-use scioxide_core::ExecutionContext;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、XCCY のデルタに使わない関数がある"
-)]
-#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
-pub mod ois;
-
-#[allow(
-    dead_code,
-    reason = "共有の fixture のうち、XCCY のデルタに使わない関数がある"
-)]
-#[path = "../../rhodium-pricing/tests/fixtures/xccy.rs"]
-pub mod xccy;
+use rhodium_risk::delta::quote_delta_in;
+use rhodium_testkit::{SEQUENTIAL, check_aad_vs_bump, entry_or_zero};
+use support::{BP, Recalibrations, calibrated};
 
 use xccy::{SPOT, XccyMarket, jpy_cash};
 
-const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
-const BP: f64 = 1e-4;
 /// バンプの再校正の残差の許容誤差。品質の照合手順は通常の 1e-14 の 1/100 を求めるが、XCCY のパーの
 /// 残差は倍精度の分解能で約 1.5e-16 までしか下がらない。届く最小の値で解き、`δ_V` に伝える。
 const BUMP_TOLERANCE: f64 = 1e-15;
 
-fn calibrated(
-    definition: &MarketDefinition,
-    quotes: &QuoteSet,
-    ctx: &Context<'_>,
-) -> Result<CalibratedMarket> {
-    let market = calibrate_with_residual_tolerance(definition, quotes, ctx, BUMP_TOLERANCE)?;
-    for group in market.groups() {
-        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
-            return Err(Error::Calibration(format!("{group:?}")));
-        }
-    }
-    Ok(market)
-}
-
-fn bumped(quotes: &QuoteSet, id: &QuoteId, bump: f64) -> Result<QuoteSet> {
-    QuoteSet::new(
-        quotes
-            .iter()
-            .map(|(key, value)| (key.clone(), if key == id { value + bump } else { value })),
-    )
-}
-
-fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
-    vector
-        .entries()
-        .binary_search_by(|(key, _)| key.cmp(id))
-        .ok()
-        .map_or(0.0, |index| vector.entries()[index].1)
-}
-
 /// 照合する取引と、その PV の絶対の下限（JPY）。
 fn products() -> Result<Vec<(&'static str, StandardProduct, f64)>> {
     let pair = XccyMarket::pair()?;
@@ -155,9 +111,21 @@ fn xccy_fx_forward_and_fx_swap_quote_deltas_match_bumps() -> Result<()> {
     let ctx = inputs.context();
     let definition = fixture.definition()?;
     let quotes = fixture.quotes()?;
-    let base = calibrated(&definition, &quotes, &ctx)?;
-    for (name, product, floor) in products()? {
-        let prepared = prepare_trade(&product, &base, &ctx)?;
+    let base = calibrated(&definition, &quotes, &ctx, BUMP_TOLERANCE)?;
+    // 同じ (quote, 幅) の再校正と、その市場での全取引の準備は、取引の間で共有する。
+    let products = products()?;
+    let prepare_all = |market: &CalibratedMarket| -> Result<Vec<Prepared>> {
+        products
+            .iter()
+            .map(|(_, product, _)| prepare_trade(product, market, &ctx))
+            .collect()
+    };
+    let recalibrations =
+        Recalibrations::new(&definition, &quotes, &ctx, BUMP_TOLERANCE, &prepare_all);
+    for (index, ((name, _, floor), prepared)) in
+        products.iter().zip(prepare_all(&base)?).enumerate()
+    {
+        let (name, floor) = (*name, *floor);
         let report = quote_delta_in(
             std::slice::from_ref(&prepared),
             &base,
@@ -178,10 +146,11 @@ fn xccy_fx_forward_and_fx_swap_quote_deltas_match_bumps() -> Result<()> {
             let unit = unit(&base, id)?;
             check_aad_vs_bump(
                 |h| {
-                    let market = calibrated(&definition, &bumped(&quotes, id, h)?, &ctx)?;
-                    jpy_pv(&prepare_trade(&product, &market, &ctx)?, &market)
+                    let entry = recalibrations.get(id, h)?;
+                    let (market, prepared) = &*entry;
+                    jpy_pv(&prepared[index], market)
                 },
-                sensitivity(vector, id),
+                entry_or_zero(vector.entries(), id),
                 unit,
                 BUMP_TOLERANCE * scale,
                 floor,
@@ -196,7 +165,7 @@ fn xccy_fx_forward_and_fx_swap_quote_deltas_match_bumps() -> Result<()> {
             "USDJPY.SPOT"
         };
         assert!(
-            sensitivity(vector, &expected.parse()?).abs() > floor,
+            entry_or_zero(vector.entries(), &expected.parse()?).abs() > floor,
             "{name} {expected}"
         );
     }
diff --git a/crates/rhodium-testkit/src/lib.rs b/crates/rhodium-testkit/src/lib.rs
index 4a5989eb1bdbd13d8fe47e91d1267fc2a112cbca..9aa5e6c62e391ce3ff4ba8faad5dd9c21792b133 100644
--- a/crates/rhodium-testkit/src/lib.rs
+++ b/crates/rhodium-testkit/src/lib.rs
@@ -37,17 +37,102 @@ macro_rules! compile_fail_pair {
     }};
 }
 
+/// PV、キャッシュフロー、エクスポージャーの点、分位点の相対の許容誤差 ε。
+///
+/// 許容誤差の定数は、品質の章の許容誤差の表（`docs/design/11-quality.md#許容誤差`）の ε と a を
+/// そのまま名前にしたもので、表を変えるときはここも同じ PR で直す。
+pub const PV_REL: f64 = 1e-12;
+/// 校正したパラメータと market の値（DF、forward、ボラ）の相対の許容誤差 ε。
+pub const MARKET_REL: f64 = 1e-12;
+/// 校正したパラメータと market の値の絶対の下限 a。
+pub const MARKET_ABS: f64 = 1e-14;
+/// 1 次の感応度（quote、パラメータ）と、XVA の値と感応度の相対の許容誤差 ε。
+pub const SENSITIVITY_REL: f64 = 1e-10;
+/// 陰的に解いた量（implied vol、パーレート）と回帰係数の相対の許容誤差 ε。
+pub const IMPLIED_REL: f64 = 1e-10;
+/// 陰的に解いた量と回帰係数の絶対の下限 a。
+pub const IMPLIED_ABS: f64 = 1e-12;
+/// 2 次の感応度の相対の許容誤差 ε。
+pub const SECOND_ORDER_REL: f64 = 1e-8;
+/// PV と感応度の絶対の下限 a を、通貨の最小単位の何倍にするか。
+const MINOR_UNIT_SHARE: f64 = 1e-6;
+
+/// PV と感応度の絶対の下限 a。通貨の小数の桁数 `minor_units`（JPY は 0、USD は 2）の最小単位の
+/// 1e-6 倍を返す。感応度では、感応度の単位（1bp など）あたりの値として使う。
+#[must_use]
+pub const fn minor_unit_abs(minor_units: u8) -> f64 {
+    // 定数の文脈でも使えるよう、桁数だけ 10 で割る。
+    let mut value = MINOR_UNIT_SHARE;
+    let mut digit = 0;
+    while digit < minor_units {
+        value /= 10.0;
+        digit += 1;
+    }
+    value
+}
+
+/// 決定的に一つのスレッドで実行する実行文脈。
+///
+/// テストの既定の実行文脈で、並列の結果と許容誤差の範囲で一致することを確かめる基準にもする。
+pub const SEQUENTIAL: scioxide_core::ExecutionContext<'static> =
+    scioxide_core::ExecutionContext::sequential();
+
+/// 疎なベクトルの `key` の値。キーの昇順に並んだ `entries` に無いキーは 0 とする。
+///
+/// quote デルタなどの疎なベクトルは、ちょうど 0 になった要素を持たないため、欠けたキーを 0 として
+/// 読む（品質の章の許容誤差）。
+#[must_use]
+pub fn entry_or_zero<K: Ord>(entries: &[(K, f64)], key: &K) -> f64 {
+    entries
+        .binary_search_by(|(entry, _)| entry.cmp(key))
+        .map_or(0.0, |index| entries[index].1)
+}
+
+/// 二つの疎なベクトルを、両側のキーの和集合で照合する。
+///
+/// 片側に無いキーは 0 として [`assert_close!`] と同じ式で比べるので、片側にだけある許容誤差を
+/// 超える項（余分な非ゼロ項や欠けた項）を検出する。`entries` はどちらもキーの昇順に並べる。
+///
+/// # Panics
+///
+/// いずれかのキーで値が許容誤差に収まらない場合。
+pub fn assert_sparse_close<K: Ord + std::fmt::Debug>(
+    actual: &[(K, f64)],
+    expected: &[(K, f64)],
+    relative: f64,
+    absolute: f64,
+) {
+    let keys: std::collections::BTreeSet<&K> =
+        actual.iter().chain(expected).map(|(key, _)| key).collect();
+    for key in keys {
+        let (left, right) = (entry_or_zero(actual, key), entry_or_zero(expected, key));
+        let tolerance = (relative * left.abs().max(right.abs())).max(absolute);
+        assert!(
+            left.is_finite() && right.is_finite() && (left - right).abs() <= tolerance,
+            "{key:?}: actual {left}, expected {right}, tolerance {tolerance}"
+        );
+    }
+}
+
 /// 数値を相対誤差と絶対誤差の大きい方で照合する。
 ///
 /// `|x − y| ≤ max(rel · M, abs)` を確かめる。M は `max(|x|, |y|)` で、5 番目の引数を渡すと
 /// それも M に含める。打ち消しの大きい和（ネットした PV やデルタ）は、和の項の絶対値の和を
 /// 5 番目に渡す（品質の許容誤差の M = max(|x|, |y|, Σ|項|)）。
+/// 最後に `; "{label}: "` のように書式を足すと、失敗の表示の先頭にその文脈を付ける。
 #[macro_export]
 macro_rules! assert_close {
     ($actual:expr, $expected:expr, $relative:expr, $absolute:expr $(,)?) => {
         $crate::assert_close!($actual, $expected, $relative, $absolute, 0.0)
     };
-    ($actual:expr, $expected:expr, $relative:expr, $absolute:expr, $scale:expr $(,)?) => {{
+    ($actual:expr, $expected:expr, $relative:expr, $absolute:expr; $($context:tt)+) => {
+        $crate::assert_close!($actual, $expected, $relative, $absolute, 0.0; $($context)+)
+    };
+    ($actual:expr, $expected:expr, $relative:expr, $absolute:expr, $scale:expr $(,)?) => {
+        $crate::assert_close!($actual, $expected, $relative, $absolute, $scale; "")
+    };
+    // `;` の後ろは失敗の表示の先頭に付ける文脈（表の行の名前など）で、format! と同じ書き方をする。
+    ($actual:expr, $expected:expr, $relative:expr, $absolute:expr, $scale:expr; $($context:tt)+) => {{
         let actual: f64 = $actual;
         let expected: f64 = $expected;
         let relative: f64 = $relative;
@@ -57,7 +142,8 @@ macro_rules! assert_close {
         let tolerance = (relative * magnitude).max(absolute);
         assert!(
             actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= tolerance,
-            "actual {actual}, expected {expected}, tolerance {tolerance}"
+            "{}actual {actual}, expected {expected}, tolerance {tolerance}",
+            format_args!($($context)+)
         );
     }};
 }
@@ -89,7 +175,7 @@ pub fn check_aad_vs_bump<E>(
     let fine = central(half)?;
     let truncation = (coarse - fine).abs() / 3.0;
     let solver = pv_error * (unit / half);
-    let rounding = (1e-10 * aad.abs().max(fine.abs())).max(absolute);
+    let rounding = (SENSITIVITY_REL * aad.abs().max(fine.abs())).max(absolute);
     let tolerance = 2.0 * truncation + solver + rounding;
     assert!(
         aad.is_finite() && fine.is_finite() && (aad - fine).abs() <= tolerance,
diff --git a/crates/rhodium-testkit/tests/assert_close.rs b/crates/rhodium-testkit/tests/assert_close.rs
index 06b90742f4f62c0ab3250b2cd972594974bc8c01..7b5d106bb1109216908995f44b8226b7489b480a 100644
--- a/crates/rhodium-testkit/tests/assert_close.rs
+++ b/crates/rhodium-testkit/tests/assert_close.rs
@@ -26,3 +26,11 @@ fn non_finite_values_never_match() {
     // NaN は許容誤差によらず一致しない。
     assert_close!(f64::NAN, f64::NAN, 1.0, 1.0);
 }
+
+#[test]
+#[should_panic(expected = "ACT/360: actual")]
+fn context_prefixes_the_failure_message() {
+    // 表の行ごとに照合するテストが、失敗した行を表示で特定できる。
+    let label = "ACT/360";
+    assert_close!(1.0 + 1e-11, 1.0, 1e-12, 0.0; "{label}: ");
+}
diff --git a/crates/rhodium/Cargo.toml b/crates/rhodium/Cargo.toml
index 6957fe7c312c497d7de5de601f40cd6ce92286cc..e03c0257805ef713cd403c86f5e15f41c5e24f1b 100644
--- a/crates/rhodium/Cargo.toml
+++ b/crates/rhodium/Cargo.toml
@@ -23,6 +23,7 @@ serde_json = "1"
 toml = "1.2"
 
 [dev-dependencies]
+rhodium-fixtures.workspace = true
 criterion.workspace = true
 rhodium-testkit.workspace = true
 scioxide-core.workspace = true
diff --git a/crates/rhodium/benches/active_market.rs b/crates/rhodium/benches/active_market.rs
index 0bb50a87b9b4477bf099bb77baa0e5dc053da447..fe2790bc462d38da9b930c5faa160d4fb09e8e52 100644
--- a/crates/rhodium/benches/active_market.rs
+++ b/crates/rhodium/benches/active_market.rs
@@ -1,13 +1,10 @@
 //! 固定した 1,000 件の OIS の quote デルタを、取引が読まない曲線の本数だけ変えて測る。
 
-#[allow(dead_code, reason = "校正用 fixture の定義関数は未使用曲線を足せない")]
-#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
-pub mod ois;
-
 use std::{hint::black_box, time::Duration};
 
 use criterion::{BenchmarkId, Criterion};
-use rhodium_core::{Result, parallel::ExecutionContext, time::Instant};
+use rhodium_core::{Result, time::Instant};
+use rhodium_fixtures::{COLLATERAL, ois};
 use rhodium_market::{
     curve::{CurveRole, CurveSpec, Extrapolation},
     quote::{QuoteId, QuoteKind, QuoteSet},
@@ -19,14 +16,11 @@ use rhodium_pricing::{
         CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
     },
     prepare::Prepared,
-    resolve::CollateralKey,
 };
 use rhodium_risk::delta::quote_delta;
+use rhodium_testkit::SEQUENTIAL;
 
 const TRADES: usize = 1_000;
-const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
-const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
-
 fn portfolio(extra_curves: usize) -> Result<(CalibratedMarket, Vec<Prepared>)> {
     let [case, _] = ois::cases()?;
     let inputs = ois::Inputs::new(&[&case], case.as_of, "")?;
diff --git a/crates/rhodium/benches/m0.rs b/crates/rhodium/benches/m0.rs
index 8ecdcbc02897876a9deb6f0439d219b7a47c61b2..01dc6aad9f1627a8937a16c38d214050ea72f9f0 100644
--- a/crates/rhodium/benches/m0.rs
+++ b/crates/rhodium/benches/m0.rs
@@ -4,13 +4,6 @@
 //! ディレクトリに書く合成のデスク（JPY と USD の 1M〜30Y の OIS 曲線と慣行どおりの OIS）。
 //! 結果は日付、コミット、マシンとともに `docs/benchmarks/` に残す。
 
-#[allow(
-    dead_code,
-    reason = "テストと共有する生成器のうち、計測で使わない定数がある"
-)]
-#[path = "../tests/support/desk.rs"]
-mod desk;
-
 use std::{hint::black_box, num::NonZeroUsize, path::Path};
 
 use criterion::{BatchSize, Criterion};
@@ -20,6 +13,7 @@ use rhodium::{
     workspace::Workspace,
 };
 use rhodium_core::time::Instant;
+use rhodium_fixtures::desk;
 use rhodium_market::quote::QuoteOverride;
 use tempfile::TempDir;
 
diff --git a/crates/rhodium/benches/m1.rs b/crates/rhodium/benches/m1.rs
index d6dd0c3e7c0a2563a6412ba5d46420f6cc96e202..8966959f8f3543c92ae4209563605c4b06026feb 100644
--- a/crates/rhodium/benches/m1.rs
+++ b/crates/rhodium/benches/m1.rs
@@ -6,12 +6,7 @@
 //! 結果は日付、コミット、マシンとともに `docs/benchmarks/m1.md` に残す。
 //! 計測の名前は `m1_rates/` で始め、M0 のデスクで M1-7 を測る `m0` の `m1/quote_update_10k` と区別する。
 
-use std::{
-    fs,
-    hint::black_box,
-    num::NonZeroUsize,
-    path::{Path, PathBuf},
-};
+use std::{fs, hint::black_box, num::NonZeroUsize, path::Path};
 
 use criterion::{BatchSize, Criterion};
 use rhodium::{
@@ -20,6 +15,7 @@ use rhodium::{
     workspace::Workspace,
 };
 use rhodium_core::time::Instant;
+use rhodium_fixtures::files::{copy_dir, example};
 use rhodium_market::quote::QuoteOverride;
 use tempfile::TempDir;
 
@@ -28,38 +24,13 @@ const DESK: usize = 10_000;
 /// 10 万件のデスク（性能の目安の「線形な取引 10 万件の PV と quote デルタ」）。
 const LARGE_DESK: usize = 100_000;
 
-#[expect(
-    clippy::disallowed_macros,
-    reason = "リポジトリのサンプルの場所をベンチのソースの位置から決める"
-)]
-fn sample() -> PathBuf {
-    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m1-rates")
-}
-
-#[expect(
-    clippy::expect_used,
-    reason = "計測の入力を作れなければ、計測を続けられない"
-)]
-fn copy_dir(from: &Path, to: &Path) {
-    for entry in fs::read_dir(from).expect("sample directory") {
-        let entry = entry.expect("sample entry");
-        let target = to.join(entry.file_name());
-        if entry.file_type().expect("entry type").is_dir() {
-            fs::create_dir(&target).expect("directory");
-            copy_dir(&entry.path(), &target);
-        } else {
-            fs::copy(entry.path(), target).expect("file");
-        }
-    }
-}
-
 /// サンプルの取引ファイルの `[[trade]]` ごとの塊（文書の `version` の行を除く）。
 #[expect(
     clippy::expect_used,
     reason = "計測の入力を作れなければ、計測を続けられない"
 )]
 fn trade_blocks() -> Vec<String> {
-    let mut paths: Vec<_> = fs::read_dir(sample().join("trades"))
+    let mut paths: Vec<_> = fs::read_dir(example("m1-rates").join("trades"))
         .expect("trades directory")
         .map(|entry| entry.expect("trade file").path())
         .collect();
@@ -83,7 +54,7 @@ fn trade_blocks() -> Vec<String> {
 )]
 fn desk(trades: usize) -> TempDir {
     let dir = tempfile::tempdir().expect("temporary directory");
-    copy_dir(&sample(), dir.path());
+    copy_dir(&example("m1-rates"), dir.path()).expect("copy the sample workspace");
     let trades_dir = dir.path().join("trades");
     fs::remove_dir_all(&trades_dir).expect("remove sample trades");
     fs::create_dir(&trades_dir).expect("trades directory");
diff --git a/crates/rhodium/examples/issue_180_workspace.rs b/crates/rhodium/examples/issue_180_workspace.rs
index 26bf460d690c2c44b00226ea375ace0a71932835..a66439b443f5c44973a3e9a327f58fd6375f021e 100644
--- a/crates/rhodium/examples/issue_180_workspace.rs
+++ b/crates/rhodium/examples/issue_180_workspace.rs
@@ -1,9 +1,5 @@
 //! Issue #180 のデスク生成、読み込み・準備・最大 RSS の計測、前後比較のための数値出力。
 
-#[allow(dead_code, reason = "計測では desk の締め日の定数を使わない")]
-#[path = "../tests/support/desk.rs"]
-mod desk;
-
 use std::{hint::black_box, num::NonZeroUsize, path::Path, time::Instant};
 
 use rhodium::{
@@ -11,6 +7,7 @@ use rhodium::{
     session::{GroupOutcome, Session},
     workspace::Workspace,
 };
+use rhodium_fixtures::desk;
 
 #[expect(
     clippy::disallowed_methods,
diff --git a/crates/rhodium/tests/calendar_io.rs b/crates/rhodium/tests/calendar_io.rs
index 25b0090efa0b00c0ff0ac12002f7c066658e6dab..b0e1438fb249a6f7386dfd173ae28f93ac270848 100644
--- a/crates/rhodium/tests/calendar_io.rs
+++ b/crates/rhodium/tests/calendar_io.rs
@@ -188,25 +188,40 @@ fn calendar_documents_reject_invalid_shapes_and_references() -> Result<()> {
             "{json}"
         );
     }
-    for invalid_json in [
-        r#"[1,"JPTO"]"#,
-        r#"{"version":1,"calendar":{"adjusted":["JPTO",[],[]]}}"#,
-        r#"{"version":1,"calendar":{"union":["JPTO","USNY"],"union":["USNY","JPTO"]}}"#,
-        r#"{"version":1,"calendar":{"union":["JPTO","USNY"],"intersection":["JPTO","USNY"]}}"#,
+    // 文書の形の誤りは、どれも理由を付けた InvalidData にする。配列の文書、表でない調整、
+    // 演算の重複と二つの演算、未知のフィールド。
+    for (invalid_json, reason) in [
+        (r#"[1,"JPTO"]"#, "calendar JSON document must be an object"),
+        (
+            r#"{"version":1,"calendar":{"adjusted":["JPTO",[],[]]}}"#,
+            "expected a table with inner, add, and remove",
+        ),
+        (
+            r#"{"version":1,"calendar":{"union":["JPTO","USNY"],"union":["USNY","JPTO"]}}"#,
+            "calendar table must contain exactly one operation",
+        ),
+        (
+            r#"{"version":1,"calendar":{"union":["JPTO","USNY"],"intersection":["JPTO","USNY"]}}"#,
+            "calendar table must contain exactly one operation",
+        ),
+        (
+            r#"{"version":1,"calendar":"JPTO","unknown":true}"#,
+            "unknown field `unknown`",
+        ),
     ] {
         assert!(
-            read_calendar_json(invalid_json, &refdata).is_err(),
+            matches!(
+                read_calendar_json(invalid_json, &refdata),
+                Err(Error::InvalidData(message)) if message.contains(reason)
+            ),
             "{invalid_json}"
         );
     }
-    assert!(read_calendar_toml("version = 2\ncalendar = \"JPTO\"", &refdata).is_err());
-    assert!(
-        read_calendar_json(
-            r#"{"version":1,"calendar":"JPTO","unknown":true}"#,
-            &refdata
-        )
-        .is_err()
-    );
+    // 読めない版の文書も同じく InvalidData にする。
+    assert!(matches!(
+        read_calendar_toml("version = 2\ncalendar = \"JPTO\"", &refdata),
+        Err(Error::InvalidData(message)) if message.contains("unsupported calendar document version 2")
+    ));
     Ok(())
 }
 
@@ -226,11 +241,15 @@ kind = { overnight_rate = { ccy = "JPY", day_count = "ACT365F" } }
     );
     let json = serde_json::to_string(&index)?;
     assert_eq!(serde_json::from_str::<IndexDefinition>(&json)?, index);
+    // 公表の規則の未知のフィールドは、どの形の公表にも当たらないので読めない。
+    let unknown = toml::from_str::<IndexDefinition>(
+        &text.replace("business_days = 1", "business_days = 1, unknown = 0"),
+    )
+    .map(drop)
+    .map_err(|error| error.to_string());
     assert!(
-        toml::from_str::<IndexDefinition>(
-            &text.replace("business_days = 1", "business_days = 1, unknown = 0")
-        )
-        .is_err()
+        matches!(&unknown, Err(message) if message.contains("did not match any variant of untagged enum PublicationRule")),
+        "{unknown:?}"
     );
     Ok(())
 }
diff --git a/crates/rhodium/tests/fx_conventions_io.rs b/crates/rhodium/tests/fx_conventions_io.rs
index c13bec6e7d78b8ed5cb44448461c91d66db28d37..cd5a7c3c141672c3bb99b2f4c04ef167bb0e0113 100644
--- a/crates/rhodium/tests/fx_conventions_io.rs
+++ b/crates/rhodium/tests/fx_conventions_io.rs
@@ -3,30 +3,14 @@
 //! スポットラグ、同日受渡の締めと pip の桁は参照データに一つだけ置き、市場の定義の `[[fx]]` は
 //! quote と担保だけを書く。古い形式で市場の定義に `settlement` を書いた入力は未知のフィールドとして拒否する。
 
-#[allow(
-    dead_code,
-    reason = "共有のデスクの生成器のうち、ここで使わない定数がある"
-)]
-#[path = "support/desk.rs"]
-mod desk;
-
-use std::{fs, io::Write as _, path::Path};
+use std::fs;
 
 use rhodium::{prelude::ExecutionContext, session::Session, workspace::Workspace};
 use rhodium_core::{Error, Result, fx::CcyPair};
-
-#[expect(clippy::needless_pass_by_value, reason = "map_err に渡す変換関数")]
-fn io_error(error: std::io::Error) -> Error {
-    Error::Internal(error.to_string())
-}
-
-fn append(path: &Path, suffix: &str) -> Result<()> {
-    fs::OpenOptions::new()
-        .append(true)
-        .open(path)
-        .and_then(|mut file| file.write_all(suffix.as_bytes()))
-        .map_err(io_error)
-}
+use rhodium_fixtures::{
+    desk,
+    files::{append, io_error},
+};
 
 /// 市場の定義の USDJPY のスポットだけの `[[fx]]`。受渡の慣行は書かない。
 const MARKET_FX: &str = r#"
diff --git a/crates/rhodium/tests/m1_rates.rs b/crates/rhodium/tests/m1_rates.rs
index fbe9551ab9dd7c933d8bd9034a274c1db9683b52..7c6560152e43961c6232c0231ff52b714d88d5a9 100644
--- a/crates/rhodium/tests/m1_rates.rs
+++ b/crates/rhodium/tests/m1_rates.rs
@@ -4,93 +4,45 @@
 //! USD の SOFR、USDJPY の外貨担保曲線（FX ポイントと `MtM` XCCY）を校正し、M1 のリニア商品と
 //! CSA の 4 つの場合の取引の PV、報告通貨の合計と quote デルタを求める。
 
+pub mod support;
+
 use std::{
-    collections::BTreeSet,
+    collections::{BTreeMap, BTreeSet},
     fs,
-    num::NonZeroUsize,
-    path::{Path, PathBuf},
+    path::Path,
     ptr,
 };
 
-use rhodium::{
-    prelude::{ExecutionContext, ThreadExecutor},
-    session::{GroupOutcome, PortfolioValuation, Session, Valuation},
-    workspace::Workspace,
-};
+use rhodium::session::Session;
 use rhodium_core::{Error, Result, Warning, currency::Currency, id::TradeId};
+use rhodium_fixtures::{
+    files::{copy_dir, example, io_error},
+    solved_within,
+};
 use rhodium_market::{
     quote::{QuoteId, QuoteKind, QuoteOverride},
+    structure::BlockId,
     view::Slots,
 };
 use rhodium_pricing::config::{AnalyticPricer, Method};
 use rhodium_risk::vector::RiskVector;
-use rhodium_testkit::{assert_close, check_aad_vs_bump};
+use rhodium_testkit::{
+    PV_REL, SEQUENTIAL, assert_close, check_aad_vs_bump, entry_or_zero, minor_unit_abs,
+};
+use support::{reporting_amounts, session_on, valuation};
 
-const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
-const REL: f64 = 1e-12;
-const JPY_ABS: f64 = 1e-6;
+const JPY_ABS: f64 = minor_unit_abs(0);
 const TONA_10Y: &str = "JPY.TONA.OIS.10Y";
-
-#[expect(
-    clippy::disallowed_macros,
-    reason = "リポジトリのサンプルの場所をテストのソースの位置から決める"
-)]
-fn sample() -> PathBuf {
-    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m1-rates")
-}
-
-#[expect(clippy::needless_pass_by_value, reason = "map_err に渡す変換関数")]
-fn io_error(error: std::io::Error) -> Error {
-    Error::Internal(error.to_string())
-}
-
-fn copy_dir(from: &Path, to: &Path) -> Result<()> {
-    for entry in fs::read_dir(from).map_err(io_error)? {
-        let entry = entry.map_err(io_error)?;
-        let target = to.join(entry.file_name());
-        if entry.file_type().map_err(io_error)?.is_dir() {
-            fs::create_dir(&target).map_err(io_error)?;
-            copy_dir(&entry.path(), &target)?;
-        } else {
-            fs::copy(entry.path(), target).map_err(io_error)?;
-        }
-    }
-    Ok(())
-}
-
-fn session(dir: &Path, execution: ExecutionContext<'_>) -> Result<Session> {
-    let workspace = Workspace::load(dir)?;
-    let as_of = workspace.resolve_date("2026-11-24".parse()?)?;
-    Session::new(workspace, as_of, execution)
-}
-
-fn valuation<'a>(report: &'a PortfolioValuation, id: &str) -> Result<&'a Valuation> {
-    report
-        .results
-        .iter()
-        .find(|(trade, _)| trade.as_str() == id)
-        .ok_or_else(|| Error::Internal(format!("missing result for {id}")))?
-        .1
-        .as_ref()
-        .map_err(Clone::clone)
-}
-
-fn has(value: &Valuation, warning: &Warning) -> bool {
-    value.warnings.contains(warning)
-}
+/// サンプルの評価日。価格設定の締め時刻（東京の 17 時）で評価の瞬間に解決する。
+const DATE: &str = "2026-11-24";
 
 #[test]
 fn m1_rates_workspace_calibrates_and_values_every_trade() -> Result<()> {
     // すべての校正グループが残差の許容誤差（1e-14）で解け、すべての取引が評価でき、合計から除く取引は無い。
-    let session = session(&sample(), SEQUENTIAL)?;
+    let session = session_on(&example("m1-rates"), DATE, SEQUENTIAL)?;
     let market = session.calibrated_market();
     assert_eq!(market.groups().len(), 9);
-    for group in market.groups() {
-        let GroupOutcome::Solved(solution) = group.outcome() else {
-            return Err(Error::Calibration(format!("{group:?}")));
-        };
-        assert!(solution.residual_norm() <= 1e-14, "{:?}", group.blocks());
-    }
+    solved_within(market, 1e-14)?;
     let report = session.value_portfolio(SEQUENTIAL)?;
     assert_eq!(report.excluded, [] as [TradeId; 0]);
     assert_eq!(report.results.len(), 18);
@@ -98,50 +50,84 @@ fn m1_rates_workspace_calibrates_and_values_every_trade() -> Result<()> {
     assert_eq!(total.ccy(), Currency::JPY);
     assert!(total.amount().is_finite());
 
-    // 先物は方式の表の既定で解析プライサー、他はリニア評価。
-    for (id, value) in &report.results {
-        let value = value.as_ref().map_err(Clone::clone)?;
-        let expected = if id.as_str() == "JPY-TONA-FUT" {
-            Method::Analytic(AnalyticPricer::RateFuture)
-        } else {
-            Method::Linear
-        };
-        assert_eq!(value.method, expected, "{id}");
-    }
-    // 予想の測度: TONA 担保の TIBOR の IRS は担保の違いの近似を使わない。25bp の CSA の IRS と、
-    // JPY 担保で割り引く SOFR の OIS と XCCY は使う。USD 担保と無担保の SOFR の OIS は使わない。
-    let cross = Warning::CrossCollateralProjection;
-    for (id, expected) in [
-        ("JPY-IRS-5Y", false),
-        ("JPY-BASIS-3M-5Y", false),
-        ("JPY-IRS-5Y-CASH-25", true),
-        ("USD-OIS-JPY-CASH", true),
-        ("USD-OIS-CTD", true),
-        ("USD-OIS-USD-CASH", false),
-        ("USD-OIS-UNCOLLATERALIZED", false),
-        ("USDJPY-MTM-XCCY-2Y", true),
-        ("USDJPY-FWD-1Y", false),
-    ] {
-        assert_eq!(has(valuation(&report, id)?, &cross), expected, "{id}");
-    }
-    assert!(has(
-        valuation(&report, "JPY-FRA-6X12")?,
-        &Warning::FraDiscounting
-    ));
-    for id in ["JPY-OIS-LOOKBACK", "JPY-OIS-LOCKOUT"] {
-        assert!(
-            has(valuation(&report, id)?, &Warning::DailyForwardProduct),
-            "{id}"
-        );
-    }
-    assert!(!has(
-        valuation(&report, "JPY-OIS-SHIFT")?,
-        &Warning::DailyForwardProduct
-    ));
-    assert!(has(
-        valuation(&report, "USDJPY-MTM-XCCY-2Y")?,
-        &Warning::MtmFxForwardProduct
-    ));
+    // 取引ごとの評価方式と近似の警告の全体を、期待の表と比べる（一部の取引や一部の警告だけを
+    // 確かめると、他の取引に増えた警告や消えた警告を見逃すため）。
+    // - 方式: 先物は方式の表の既定で解析プライサー、他はリニア評価。
+    // - 予想の測度（CrossCollateralProjection）: TONA 担保の TIBOR の IRS とベーシスは担保の違いの近似を
+    //   使わない。25bp の CSA の IRS と、JPY 担保で割り引く SOFR の OIS（CTD を含む）と二つの XCCY は使う。
+    //   USD 担保と無担保の SOFR の OIS と FX のフォワード・スワップは使わない。
+    // - 規約の近似: FRA は FraDiscounting、lookback と lockout は DailyForwardProduct（observation shift は
+    //   記録しない）、MtM XCCY は MtmFxForwardProduct。
+    // - PaymentDelayIgnored: 観測期間の終わりより後に払う脚（ターム金利のクーポン、2 営業日後払いの
+    //   OIS、先物の値洗い）。預金、FX、XCCY と、残りの受渡の日付が観測期間と揃う一部解約の IRS は記録しない。
+    let (linear, futures) = (Method::Linear, Method::Analytic(AnalyticPricer::RateFuture));
+    let (cross, delay) = (
+        Warning::CrossCollateralProjection,
+        Warning::PaymentDelayIgnored,
+    );
+    let expected = [
+        ("JPY-BASIS-3M-5Y", linear, vec![delay.clone()]),
+        ("JPY-DEPO-3M", linear, vec![]),
+        ("JPY-FRA-6X12", linear, vec![Warning::FraDiscounting]),
+        ("JPY-IRS-5Y", linear, vec![delay.clone()]),
+        (
+            "JPY-IRS-5Y-CASH-25",
+            linear,
+            vec![delay.clone(), cross.clone()],
+        ),
+        ("JPY-IRS-AMENDED", linear, vec![]),
+        (
+            "JPY-OIS-LOCKOUT",
+            linear,
+            vec![delay.clone(), Warning::DailyForwardProduct],
+        ),
+        (
+            "JPY-OIS-LOOKBACK",
+            linear,
+            vec![delay.clone(), Warning::DailyForwardProduct],
+        ),
+        ("JPY-OIS-SHIFT", linear, vec![delay.clone()]),
+        ("JPY-TONA-FUT", futures, vec![delay.clone()]),
+        ("USD-OIS-CTD", linear, vec![delay.clone(), cross.clone()]),
+        (
+            "USD-OIS-JPY-CASH",
+            linear,
+            vec![delay.clone(), cross.clone()],
+        ),
+        ("USD-OIS-UNCOLLATERALIZED", linear, vec![delay.clone()]),
+        ("USD-OIS-USD-CASH", linear, vec![delay]),
+        ("USDJPY-FWD-1Y", linear, vec![]),
+        ("USDJPY-FXSWAP-6M", linear, vec![]),
+        (
+            "USDJPY-MTM-XCCY-2Y",
+            linear,
+            vec![cross.clone(), Warning::MtmFxForwardProduct],
+        ),
+        ("USDJPY-XCCY-2Y", linear, vec![cross]),
+    ]
+    .into_iter()
+    .map(|(id, method, warnings)| {
+        Ok((
+            id.parse::<TradeId>()?,
+            (method, warnings.into_iter().collect::<BTreeSet<_>>()),
+        ))
+    })
+    .collect::<Result<BTreeMap<_, _>>>()?;
+    let actual = report
+        .results
+        .iter()
+        .map(|(id, value)| {
+            let value = value.as_ref().map_err(Clone::clone)?;
+            Ok((
+                id.clone(),
+                (
+                    value.method,
+                    value.warnings.iter().cloned().collect::<BTreeSet<_>>(),
+                ),
+            ))
+        })
+        .collect::<Result<BTreeMap<_, _>>>()?;
+    assert_eq!(actual, expected);
 
     // CSA の 4 つの場合: USD 担保と無担保（USD-SOFR で割り引く）は同じ PV、JPY 担保は外貨担保曲線で、
     // CTD は毎年の平均フォワードの大きい方で割り引くので、PV は JPY 担保以下。
@@ -149,7 +135,7 @@ fn m1_rates_workspace_calibrates_and_values_every_trade() -> Result<()> {
     assert_close!(
         usd("USD-OIS-USD-CASH")?,
         usd("USD-OIS-UNCOLLATERALIZED")?,
-        REL,
+        PV_REL,
         1e-8
     );
     assert!((usd("USD-OIS-JPY-CASH")? - usd("USD-OIS-USD-CASH")?).abs() > 1.0);
@@ -183,9 +169,9 @@ fn amended_trade_matches_the_scaled_original_and_its_fee() -> Result<()> {
     // 一部解約は効力の瞬間（評価時点）以後の受渡をすべて 0.4 倍にし、手数料 200 万円を別に受け取る。
     // 同じ IRS を変更せずに書いた写しの PV から、0.4 × PV と手数料の PV の和を手計算で求める。
     let dir = tempfile::tempdir().map_err(io_error)?;
-    copy_dir(&sample(), dir.path())?;
+    copy_dir(&example("m1-rates"), dir.path())?;
     let amended = valuation(
-        &session(dir.path(), SEQUENTIAL)?.value_portfolio(SEQUENTIAL)?,
+        &session_on(dir.path(), DATE, SEQUENTIAL)?.value_portfolio(SEQUENTIAL)?,
         "JPY-IRS-AMENDED",
     )?
     .pv
@@ -201,7 +187,7 @@ fn amended_trade_matches_the_scaled_original_and_its_fee() -> Result<()> {
         "type = \"amended\"\nfrom = 2026-11-24T08:00:00Z\nafter = { scale = 0.4 }\n\n[trade.product.before]\n",
         "",
     )?;
-    let original_session = session(dir.path(), SEQUENTIAL)?;
+    let original_session = session_on(dir.path(), DATE, SEQUENTIAL)?;
     let original = valuation(
         &original_session.value_portfolio(SEQUENTIAL)?,
         "JPY-IRS-AMENDED",
@@ -212,21 +198,13 @@ fn amended_trade_matches_the_scaled_original_and_its_fee() -> Result<()> {
     assert_close!(
         amended,
         0.4 * original + fee,
-        REL,
+        PV_REL,
         JPY_ABS,
         original.abs() + fee
     );
     Ok(())
 }
 
-fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
-    vector
-        .entries()
-        .iter()
-        .find_map(|(key, value)| (key == id).then_some(*value))
-        .unwrap_or(0.0)
-}
-
 fn reporting(session: &Session) -> Result<RiskVector<QuoteId>> {
     let report = session.quote_delta_reporting(SEQUENTIAL)?;
     assert!(report.failures().is_empty());
@@ -242,38 +220,48 @@ fn m1_rates_quote_update_re_solves_downstream_groups_only() -> Result<()> {
     // 外貨担保曲線だけを解き直し、他のグループ（年末の跳び、階段、SOFR、FX と LGM のマーク）と
     // 取引の準備を共有する。結果は更新後の quote から新規に作った Session と許容誤差で一致する。
     let dir = tempfile::tempdir().map_err(io_error)?;
-    copy_dir(&sample(), dir.path())?;
-    let base = session(dir.path(), SEQUENTIAL)?;
+    copy_dir(&example("m1-rates"), dir.path())?;
+    let base = session_on(dir.path(), DATE, SEQUENTIAL)?;
     let id: QuoteId = TONA_10Y.parse()?;
     let overrides = base.quote_overrides([QuoteOverride::Add {
         pattern: id.as_str().parse()?,
         shift: 1e-4,
     }])?;
     let updated = base.with_quotes(&overrides, SEQUENTIAL)?;
-    let resolved: BTreeSet<String> = base
-        .calibrated_market()
-        .groups()
+    // 解き直したグループは、結果の記録（outcome）が新しく作られる。グループの並びは同じ定義から
+    // 決まるので、両方の並びの長さを確かめてから位置で対応させる。
+    let [old_groups, new_groups] =
+        [&base, &updated].map(|session| session.calibrated_market().groups());
+    assert_eq!(old_groups.len(), new_groups.len());
+    let resolved: BTreeSet<BlockId> = old_groups
         .iter()
-        .zip(updated.calibrated_market().groups())
+        .zip(new_groups)
         .filter(|(old, new)| !ptr::eq(old.outcome(), new.outcome()))
-        .map(|(old, _)| format!("{} {}", old.blocks()[0].curve(), old.blocks()[0].part()))
+        .flat_map(|(old, _)| old.blocks().iter().cloned())
         .collect();
+    let expected = [
+        ("JPY-TIBOR-3M", "basis"),
+        ("JPY-TIBOR-6M", "basis"),
+        ("JPY-TONA", "long"),
+        ("USD-under-JPY", "basis"),
+    ]
+    .into_iter()
+    .map(|(curve, part)| BlockId::new(curve.parse()?, part))
+    .collect::<Result<BTreeSet<_>>>()?;
+    assert_eq!(resolved, expected);
+    // 取引の準備は共有する。取引の ID の列全体が同じで、各取引の警告の記録は同じ場所を指す。
+    let [before, after] = [&base, &updated].map(|session| {
+        session
+            .preparation_warnings()
+            .map(|(id, warnings)| (id.clone(), warnings))
+            .collect::<Vec<_>>()
+    });
     assert_eq!(
-        resolved,
-        [
-            "JPY-TIBOR-3M basis",
-            "JPY-TIBOR-6M basis",
-            "JPY-TONA long",
-            "USD-under-JPY basis"
-        ]
-        .map(String::from)
-        .into()
+        before.iter().map(|(id, _)| id).collect::<Vec<_>>(),
+        after.iter().map(|(id, _)| id).collect::<Vec<_>>()
     );
-    for ((before, _), (after, _)) in base
-        .preparation_warnings()
-        .zip(updated.preparation_warnings())
-    {
-        assert!(ptr::eq(before, after));
+    for ((_, before), (_, after)) in before.iter().zip(&after) {
+        assert!(ptr::eq(*before, *after));
     }
 
     let snapshot = dir.path().join("quotes/snapshot.csv");
@@ -283,29 +271,22 @@ fn m1_rates_quote_update_re_solves_downstream_groups_only() -> Result<()> {
         text.replace(&format!("{TONA_10Y},1.32%"), &format!("{TONA_10Y},1.33%")),
     )
     .map_err(io_error)?;
-    let fresh = session(dir.path(), SEQUENTIAL)?;
+    let fresh = session_on(dir.path(), DATE, SEQUENTIAL)?;
     // 更新は前の θ から、新規は固定の初期値から解くので、θ はソルバーの誤差を伝えた幅で一致する。
-    for (lhs, rhs) in updated
-        .calibrated_market()
-        .market()
-        .theta()
-        .iter()
-        .zip(fresh.calibrated_market().market().theta())
-    {
+    let [lhs, rhs] = [&updated, &fresh].map(|session| session.calibrated_market().market().theta());
+    assert_eq!(lhs.len(), rhs.len());
+    for (lhs, rhs) in lhs.iter().zip(rhs) {
         assert_close!(*lhs, *rhs, 1e-10, 1e-12);
     }
-    let [lhs, rhs] = [&updated, &fresh].map(|session| session.value_portfolio(SEQUENTIAL));
+    // 取引の ID の列全体と、各取引の報告通貨の PV を比べる。
+    let [lhs, rhs] = [&updated, &fresh].map(|session| reporting_amounts(session, SEQUENTIAL));
     let (lhs, rhs) = (lhs?, rhs?);
-    for ((id, left), (_, right)) in lhs.results.iter().zip(&rhs.results) {
-        let [left, right] = [left, right].map(|value| value.as_ref().map_err(Clone::clone));
-        assert_close!(
-            left?.reporting.amount(),
-            right?.reporting.amount(),
-            1e-10,
-            JPY_ABS,
-            1e9
-        );
-        let _ = id;
+    assert_eq!(
+        lhs.keys().collect::<Vec<_>>(),
+        rhs.keys().collect::<Vec<_>>()
+    );
+    for (id, value) in &lhs {
+        assert_close!(*value, rhs[id], 1e-10, JPY_ABS, 1e9; "{id}: ");
     }
     let [lhs, rhs] = [reporting(&updated)?, reporting(&fresh)?];
     let scale = rhs
@@ -315,8 +296,8 @@ fn m1_rates_quote_update_re_solves_downstream_groups_only() -> Result<()> {
         .sum::<f64>();
     for id in lhs.entries().iter().chain(rhs.entries()).map(|(id, _)| id) {
         assert_close!(
-            sensitivity(&lhs, id),
-            sensitivity(&rhs, id),
+            entry_or_zero(lhs.entries(), id),
+            entry_or_zero(rhs.entries(), id),
             1e-8,
             JPY_ABS,
             scale
@@ -325,38 +306,6 @@ fn m1_rates_quote_update_re_solves_downstream_groups_only() -> Result<()> {
     Ok(())
 }
 
-#[test]
-fn m1_rates_values_do_not_depend_on_thread_count() -> Result<()> {
-    // 取引の準備、PV、報告通貨のデルタを 3 つのワーカーに分けても、逐次と許容誤差で一致する。
-    let sequential = session(&sample(), SEQUENTIAL)?;
-    let executor = ThreadExecutor::new(NonZeroUsize::MIN.saturating_add(2));
-    let parallel_execution = ExecutionContext::new(&executor);
-    let parallel = session(&sample(), parallel_execution)?;
-    let [lhs, rhs] = [
-        sequential.value_portfolio(SEQUENTIAL)?,
-        parallel.value_portfolio(parallel_execution)?,
-    ];
-    for ((id, left), (right_id, right)) in lhs.results.iter().zip(&rhs.results) {
-        assert_eq!(id, right_id);
-        let [left, right] = [left, right].map(|value| value.as_ref().map_err(Clone::clone));
-        assert_close!(
-            left?.reporting.amount(),
-            right?.reporting.amount(),
-            REL,
-            JPY_ABS
-        );
-    }
-    let lhs = reporting(&sequential)?;
-    let report = parallel.quote_delta_reporting(parallel_execution)?;
-    let [rhs] = report.by_ccy() else {
-        return Err(Error::Internal("expected one reporting currency".into()));
-    };
-    for id in lhs.entries().iter().chain(rhs.entries()).map(|(id, _)| id) {
-        assert_close!(sensitivity(&lhs, id), sensitivity(rhs, id), 1e-10, JPY_ABS);
-    }
-    Ok(())
-}
-
 #[test]
 #[expect(
     clippy::wildcard_enum_match_arm,
@@ -366,7 +315,7 @@ fn m1_rates_portfolio_quote_delta_matches_bump_and_recalibrate() -> Result<()> {
     // 最初の検証データセットのポートフォリオ全体の報告通貨（JPY）の PV について、quote デルタ
     // （AAD、FX の換算の連鎖律、多段の校正の随伴）を、quote ごとに影響範囲を再校正した中心差分と照合する。
     // 再校正は Session::with_quotes（残差の許容誤差 1e-14）なので、その誤差を δ_V に伝える。
-    let base = session(&sample(), SEQUENTIAL)?;
+    let base = session_on(&example("m1-rates"), DATE, SEQUENTIAL)?;
     let delta = reporting(&base)?;
     let total = |session: &Session| -> Result<f64> {
         Ok(session.value_portfolio(SEQUENTIAL)?.total?.amount())
@@ -407,7 +356,7 @@ fn m1_rates_portfolio_quote_delta_matches_bump_and_recalibrate() -> Result<()> {
                 }])?;
                 total(&base.with_quotes(&overrides, SEQUENTIAL)?)
             },
-            sensitivity(&delta, id) * width / report_unit,
+            entry_or_zero(delta.entries(), id) * width / report_unit,
             width,
             1e-14 * scale,
             JPY_ABS,
diff --git a/crates/rhodium/tests/market_definition_io.rs b/crates/rhodium/tests/market_definition_io.rs
index 1daeab48410191cb3b4db3c672ed38b23a8f0662..2509d49fccfedb3d75a28ea9156cdaa4b7004a67 100644
--- a/crates/rhodium/tests/market_definition_io.rs
+++ b/crates/rhodium/tests/market_definition_io.rs
@@ -3,17 +3,16 @@
 //! ファイルの曲線は `CurveExpr` の serde の表現そのもので、瞬間は UTC の日時で明示する。
 //! 会合表や年末の規則を参照データから解決する書き方（M6）はまだ読まない。
 
-#[allow(
-    dead_code,
-    reason = "共有のデスクの生成器のうち、ここで使わない定数がある"
-)]
-#[path = "support/desk.rs"]
-mod desk;
+pub mod support;
 
-use std::{fmt::Write as _, fs, io::Write as _, path::Path};
+use std::{fmt::Write as _, fs};
 
-use rhodium::{prelude::ExecutionContext, session::Session, workspace::Workspace};
+use rhodium::workspace::Workspace;
 use rhodium_core::{Error, Result, time::Instant};
+use rhodium_fixtures::{
+    desk,
+    files::{append, io_error},
+};
 use rhodium_market::{
     curve::{CurveExpr, CurveId},
     market::Market,
@@ -22,26 +21,12 @@ use rhodium_market::{
 };
 use rhodium_pricing::calibration::{GroupOutcome, HessianMode};
 use rhodium_testkit::assert_close;
-
-const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
+use support::session;
 
 fn fmt_error(error: std::fmt::Error) -> Error {
     Error::Internal(error.to_string())
 }
 
-#[expect(clippy::needless_pass_by_value, reason = "map_err に渡す変換関数")]
-fn io_error(error: std::io::Error) -> Error {
-    Error::Internal(error.to_string())
-}
-
-fn append(path: &Path, suffix: &str) -> Result<()> {
-    fs::OpenOptions::new()
-        .append(true)
-        .open(path)
-        .and_then(|mut file| file.write_all(suffix.as_bytes()))
-        .map_err(io_error)
-}
-
 /// デスクの JPY の OIS の quote ID（1M、3M、6M、1Y、2Y … 30Y）とピラー（スポット + テナー + 7 日）。
 const SHORT: [&str; 4] = ["1M", "3M", "6M", "1Y"];
 const LONG: [(&str, &str); 11] = [
@@ -161,12 +146,6 @@ fn workspace(definition: &str) -> Result<tempfile::TempDir> {
     Ok(dir)
 }
 
-fn session(dir: &Path) -> Result<Session> {
-    let workspace = Workspace::load(dir)?;
-    let as_of = workspace.resolve_as_of(desk::AS_OF)?;
-    Session::new(workspace, as_of, SEQUENTIAL)
-}
-
 fn instants(dates: &[&str]) -> Result<Vec<Instant>> {
     dates
         .iter()
@@ -181,7 +160,7 @@ fn market_file_reads_tona_steps_and_turn() -> Result<()> {
     // 正方系の各グループは quote を残差の許容誤差で再現する。
     let meetings = ["2027-01-11", "2027-03-15", "2027-06-14"];
     let dir = workspace(&tona_definition(&meetings, "", "")?)?;
-    let session = session(dir.path())?;
+    let session = session(dir.path(), desk::AS_OF)?;
     let market = session.calibrated_market();
     let id: CurveId = "JPY-TONA".parse()?;
     let structure = market.market().structure();
@@ -278,7 +257,7 @@ max_forward = { parts = [{ flat = "0bp" }, { sum = [{ flat = "0bp" }, { scale =
         "",
         extra,
     )?)?;
-    let session = session(dir.path())?;
+    let session = session(dir.path(), desk::AS_OF)?;
     let market = session.calibrated_market().market();
     let at: Instant = "2027-12-03T00:00:00Z".parse()?;
     let t = rhodium_core::time::model_time(market.structure().as_of(), at).years();
@@ -331,7 +310,7 @@ scale = "1bp"
         ),
     ] {
         let dir = workspace(&tona_definition(&meetings, tolerance, &extra)?)?;
-        let session = session(dir.path())?;
+        let session = session(dir.path(), desk::AS_OF)?;
         let group = session
             .calibrated_market()
             .groups()
diff --git a/crates/rhodium/tests/parallel.rs b/crates/rhodium/tests/parallel.rs
index 8c455ad0dc940560cdb6329e33017518ff767450..be455d2b978190925e098870f74bb660e225ac0f 100644
--- a/crates/rhodium/tests/parallel.rs
+++ b/crates/rhodium/tests/parallel.rs
@@ -1,7 +1,6 @@
 //! スレッド数を変えても、取引ごとの PV と通貨別の quote デルタが許容誤差で一致することを確かめる（M0-9）。
 
-#[path = "support/desk.rs"]
-mod desk;
+pub mod support;
 
 use std::{
     collections::{BTreeMap, BTreeSet},
@@ -18,9 +17,14 @@ use rhodium::{
     workspace::Workspace,
 };
 use rhodium_core::{Error, Result, currency::Currency, id::TradeId};
+use rhodium_fixtures::{desk, files::example};
 use rhodium_market::quote::{QuoteId, QuoteKind, QuoteOverride, QuoteOverrides};
-use rhodium_testkit::{TickClock, assert_close};
+use rhodium_testkit::{
+    PV_REL, SENSITIVITY_REL, SEQUENTIAL, TickClock, assert_close, assert_sparse_close,
+    minor_unit_abs,
+};
 use scioxide_core::{ClockBinding, Deadline, ExecutionPlan, Executor, PartitionFn, SciOxideError};
+use support::{reporting_amounts, session_on};
 
 /// スレッドの実行器に処理を任せ、実行した計画の最大の区間の数を記録する。
 struct Recording {
@@ -150,6 +154,45 @@ fn pv_and_quote_delta_do_not_depend_on_thread_count() -> Result<()> {
     compare(36, &[2, 3, 4, 8])
 }
 
+#[test]
+fn m1_rates_values_do_not_depend_on_thread_count() -> Result<()> {
+    // M1 の検証データセット（先物、FRA、FX、XCCY、担保条件の違う OIS）でも、取引の準備、PV、報告通貨の
+    // デルタを 3 つのワーカーに分けた結果が、逐次と許容誤差で一致する。CLI のテストは同じ不変性を
+    // 出力の形だけで確かめるので、値の照合はここで行う。
+    let dir = example("m1-rates");
+    let executor = ThreadExecutor::new(NonZeroUsize::MIN.saturating_add(2));
+    let parallel = ExecutionContext::new(&executor);
+    let sequential_session = session_on(&dir, "2026-11-24", SEQUENTIAL)?;
+    let parallel_session = session_on(&dir, "2026-11-24", parallel)?;
+    // 取引の ID の列全体と、各取引の報告通貨の PV を比べる。
+    let lhs = reporting_amounts(&sequential_session, SEQUENTIAL)?;
+    let rhs = reporting_amounts(&parallel_session, parallel)?;
+    assert_eq!(lhs.len(), 18);
+    assert_eq!(
+        lhs.keys().collect::<Vec<_>>(),
+        rhs.keys().collect::<Vec<_>>()
+    );
+    for (id, value) in &lhs {
+        assert_close!(*value, rhs[id], PV_REL, minor_unit_abs(0); "{id}: ");
+    }
+    // 報告通貨のデルタは、欠けたキーを 0 として両方のキーの和集合で比べる。
+    let delta = |session: &Session, execution| -> Result<_> {
+        let report = session.quote_delta_reporting(execution)?;
+        assert!(report.failures().is_empty(), "{:?}", report.failures());
+        let [vector] = report.by_ccy() else {
+            return Err(Error::Internal("expected one reporting currency".into()));
+        };
+        Ok(vector.clone())
+    };
+    assert_sparse_close(
+        delta(&sequential_session, SEQUENTIAL)?.entries(),
+        delta(&parallel_session, parallel)?.entries(),
+        SENSITIVITY_REL,
+        minor_unit_abs(0),
+    );
+    Ok(())
+}
+
 #[test]
 #[ignore = "1 万件の取引を読み込んで準備する遅いテスト（M0-9 の確認で手で実行する）"]
 fn ten_thousand_swaps_do_not_depend_on_thread_count() -> Result<()> {
diff --git a/crates/rhodium/tests/portfolio.rs b/crates/rhodium/tests/portfolio.rs
index b500c79935e5fe75207142e5ef4cce71533f8d8c..b08b9642a0c11d33bce0113f5982d04f844c2b64 100644
--- a/crates/rhodium/tests/portfolio.rs
+++ b/crates/rhodium/tests/portfolio.rs
@@ -1,58 +1,27 @@
 //! ポートフォリオの集計と、準備に使う入力を変えた場合の評価を検証する。
 
-#[path = "support/desk.rs"]
-mod desk;
+pub mod support;
 
-use std::{fs, io::Write as _, num::NonZeroUsize, path::Path, ptr, sync::Arc};
+use std::{fs, num::NonZeroUsize, path::Path, ptr, sync::Arc};
 
 use rhodium::{
     prelude::{ExecutionContext, ThreadExecutor},
-    session::{PortfolioValuation, Session, Valuation},
-    workspace::Workspace,
+    session::{PortfolioValuation, Session},
 };
 use rhodium_core::{Error, Result, Warning, currency::Currency, id::TradeId, time::Instant};
+use rhodium_fixtures::{
+    desk,
+    files::{append, io_error},
+};
 use rhodium_market::{
     quote::{QuoteId, QuoteOverride},
     view::Slots,
 };
 use rhodium_pricing::config::Method;
-use rhodium_testkit::assert_close;
+use rhodium_testkit::{PV_REL, SEQUENTIAL, assert_close, minor_unit_abs};
+use support::{row, session, valuation};
 
-const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
-const REL: f64 = 1e-12;
-const JPY_ABS: f64 = 1e-6;
-
-#[expect(clippy::needless_pass_by_value, reason = "map_err に渡す変換関数")]
-fn io_error(error: std::io::Error) -> Error {
-    Error::Internal(error.to_string())
-}
-
-fn append(path: &Path, suffix: &str) -> Result<()> {
-    fs::OpenOptions::new()
-        .append(true)
-        .open(path)
-        .and_then(|mut file| file.write_all(suffix.as_bytes()))
-        .map_err(io_error)
-}
-
-fn session(root: &Path, as_of: &str) -> Result<Session> {
-    let workspace = Workspace::load(root)?;
-    let as_of = workspace.resolve_as_of(as_of)?;
-    Session::new(workspace, as_of, SEQUENTIAL)
-}
-
-fn row<'a>(report: &'a PortfolioValuation, id: &str) -> Result<&'a Result<Valuation>> {
-    report
-        .results
-        .iter()
-        .find(|(trade, _)| trade.as_str() == id)
-        .map(|(_, value)| value)
-        .ok_or_else(|| Error::Internal(format!("missing result for {id}")))
-}
-
-fn valuation<'a>(report: &'a PortfolioValuation, id: &str) -> Result<&'a Valuation> {
-    row(report, id)?.as_ref().map_err(Clone::clone)
-}
+const JPY_ABS: f64 = minor_unit_abs(0);
 
 fn excluded(report: &PortfolioValuation) -> Vec<&str> {
     report.excluded.iter().map(TradeId::as_str).collect()
@@ -133,16 +102,16 @@ fn portfolio_isolates_trade_failures_and_totals_successes() -> Result<()> {
         assert_eq!(value.method, Method::Linear);
         assert!(value.std_error.is_none());
         let expected = value.pv.get(Currency::JPY) + value.pv.get(Currency::USD) * today_fx;
-        assert_close!(value.reporting.amount(), expected, REL, JPY_ABS);
+        assert_close!(value.reporting.amount(), expected, PV_REL, JPY_ABS);
         let parallel_value = parallel_value.as_ref().map_err(Clone::clone)?;
-        assert_close!(parallel_value.reporting.amount(), expected, REL, JPY_ABS);
+        assert_close!(parallel_value.reporting.amount(), expected, PV_REL, JPY_ABS);
         sum += expected;
         scale += expected.abs();
     }
     for report in [report, parallel_report] {
         let total = report.total?;
         assert_eq!(total.ccy(), Currency::JPY);
-        assert_close!(total.amount(), sum, REL, JPY_ABS, scale);
+        assert_close!(total.amount(), sum, PV_REL, JPY_ABS, scale);
     }
 
     // USD の校正グループが失敗しても、JPY の独立した取引は評価できる。
@@ -164,7 +133,7 @@ fn portfolio_isolates_trade_failures_and_totals_successes() -> Result<()> {
     assert_eq!(excluded(&partial_report), ["T00001", "T00003", "T00004"]);
     let expected = valuation(&partial_report, "T00000")?.reporting.amount()
         + valuation(&partial_report, "T00002")?.reporting.amount();
-    assert_close!(partial_report.total?.amount(), expected, REL, JPY_ABS);
+    assert_close!(partial_report.total?.amount(), expected, PV_REL, JPY_ABS);
 
     // 両方の校正グループが失敗すると、全件が除外され合計は報告通貨のゼロになる。
     let remove_jpy = partial.quote_overrides([QuoteOverride::Remove {
@@ -177,7 +146,7 @@ fn portfolio_isolates_trade_failures_and_totals_successes() -> Result<()> {
     assert_eq!(all_excluded.results.len(), 5);
     let total = all_excluded.total?;
     assert_eq!(total.ccy(), Currency::JPY);
-    assert_close!(total.amount(), 0.0, REL, JPY_ABS);
+    assert_close!(total.amount(), 0.0, PV_REL, JPY_ABS);
 
     // quote を戻せば、同じ準備済みの取引が復旧する。欠けた fixing の取引は引き続き除外される。
     let restore = partial.quote_overrides([QuoteOverride::Set {
@@ -193,7 +162,7 @@ fn portfolio_isolates_trade_failures_and_totals_successes() -> Result<()> {
     }
     let report = restored.value_portfolio(SEQUENTIAL)?;
     assert_eq!(excluded(&report), ["T00004"]);
-    assert_close!(report.total?.amount(), sum, REL, JPY_ABS, scale);
+    assert_close!(report.total?.amount(), sum, PV_REL, JPY_ABS, scale);
     // 品質の「部分的な校正」の行: 補って再校正した市場の報告通貨のデルタも、全体を解き直した値と
     // 一致する（更新は前の θ から、全体は固定の初期値から解くので、ソルバーの誤差の幅で比べる）。
     assert_reporting_deltas_match(&restored, &complete)?;
@@ -291,7 +260,7 @@ fn portfolio_preserves_successful_trades_when_total_overflows() -> Result<()> {
     let report = session(root, desk::AS_OF)?.value_portfolio(SEQUENTIAL)?;
     assert_eq!(report.excluded, []);
     let expected = valuation(&report, "LONG-A")?.reporting.amount();
-    assert_close!(report.total?.amount(), expected, REL, JPY_ABS);
+    assert_close!(report.total?.amount(), expected, PV_REL, JPY_ABS);
     Ok(())
 }
 
@@ -390,12 +359,12 @@ fn assert_jpy_pv(
 ) -> Result<()> {
     let value = valuation(report, trade)?;
     assert_eq!(value.reporting.ccy(), Currency::JPY);
-    assert_close!(value.reporting.amount(), expected, REL, JPY_ABS);
+    assert_close!(value.reporting.amount(), expected, PV_REL, JPY_ABS);
     let pv = &session.pv(SEQUENTIAL)?[&trade.parse()?];
     assert_close!(
         pv.as_ref().map_err(Clone::clone)?.get(Currency::JPY),
         expected,
-        REL,
+        PV_REL,
         JPY_ABS
     );
     Ok(())
@@ -489,11 +458,11 @@ collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }"#;
     let terms_report = terms_changed.value_portfolio(SEQUENTIAL)?;
     let respread = 1_000_000.0 * discount(&terms_changed, "JPY-CSA", payment)?;
     assert_jpy_pv(&terms_changed, &terms_report, "CSA", respread)?;
-    assert_close!(respread, alternate, REL, JPY_ABS);
+    assert_close!(respread, alternate, PV_REL, JPY_ABS);
     assert_close!(
         corrected.value_portfolio(SEQUENTIAL)?.total?.amount(),
         corrected_report.total?.amount(),
-        REL,
+        PV_REL,
         JPY_ABS
     );
     Ok(())
diff --git a/crates/rhodium/tests/pricing_methods.rs b/crates/rhodium/tests/pricing_methods.rs
index 7caaba68e732502bc4868f9ba75114feb387dd93..c944c5b9bda056a89e6f769aae2203daae09c1eb 100644
--- a/crates/rhodium/tests/pricing_methods.rs
+++ b/crates/rhodium/tests/pricing_methods.rs
@@ -5,71 +5,30 @@
 //! 同じ式になる（`futures_margin_ir_vs_analytic`）ので、手数料を含めて PV が一致する。
 //! 変更した取引（一部解約）は標準商品の条件と一致しないので種類を持たず、表に当たらない。
 
-use std::{
-    fs,
-    io::{self, Write as _},
-    path::Path,
-};
+pub mod support;
 
-use rhodium::{
-    prelude::ExecutionContext,
-    session::{PortfolioValuation, Session, Valuation},
-    workspace::Workspace,
-};
-use rhodium_core::{Error, Result, currency::Currency};
+use std::{fs, path::Path};
+
+use rhodium::session::PortfolioValuation;
+use rhodium_core::{Result, currency::Currency};
+use rhodium_fixtures::files::{append, copy_dir, example, io_error};
 use rhodium_pricing::config::{AnalyticPricer, Method};
-use rhodium_testkit::assert_close;
+use rhodium_testkit::{PV_REL, SEQUENTIAL, assert_close, minor_unit_abs};
+use support::{session, valuation};
 
 const AS_OF: &str = "2026-11-24T08:00:00Z";
-const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
-const REL: f64 = 1e-12;
 /// USD の最小単位（セント）の 1e-6 倍。
-const USD_ABS: f64 = 1e-8;
-
-#[expect(clippy::needless_pass_by_value, reason = "map_err に渡す変換関数")]
-fn io_error(error: io::Error) -> Error {
-    Error::Internal(error.to_string())
-}
-
-fn copy_dir(from: &Path, to: &Path) -> io::Result<()> {
-    for entry in fs::read_dir(from)? {
-        let entry = entry?;
-        let target = to.join(entry.file_name());
-        if entry.file_type()?.is_dir() {
-            fs::create_dir(&target)?;
-            copy_dir(&entry.path(), &target)?;
-        } else {
-            fs::copy(entry.path(), target)?;
-        }
-    }
-    Ok(())
-}
-
-fn append(path: &Path, suffix: &str) -> Result<()> {
-    fs::OpenOptions::new()
-        .append(true)
-        .open(path)
-        .and_then(|mut file| file.write_all(suffix.as_bytes()))
-        .map_err(io_error)
-}
+const USD_ABS: f64 = minor_unit_abs(2);
 
 /// M0 のサンプルに、USD の金利先物の証券、清算指標、限月、校正の quote と LGM のマーク、
 /// 取引（手数料付き、手数料なし、一部解約）を足したワークスペースを書く。
 fn futures_workspace(root: &Path, extra_pricing: &str) -> Result<()> {
-    copy_dir(&sample(), root).map_err(io_error)?;
+    copy_dir(&example("m0-ois"), root)?;
     add_futures_refdata(root)?;
     add_futures_market(root, extra_pricing)?;
     write_futures_trades(root)
 }
 
-#[expect(
-    clippy::disallowed_macros,
-    reason = "リポジトリのサンプルの場所をテストのソースの位置から決める"
-)]
-fn sample() -> std::path::PathBuf {
-    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m0-ois")
-}
-
 /// 先物の証券、原指標と清算指標、限月の慣行を参照データに足す。
 fn add_futures_refdata(root: &Path) -> Result<()> {
     fs::write(
@@ -200,20 +159,7 @@ after = {{ scale = 0.5 }}
 fn portfolio(extra_pricing: &str) -> Result<PortfolioValuation> {
     let dir = tempfile::tempdir().map_err(io_error)?;
     futures_workspace(dir.path(), extra_pricing)?;
-    let workspace = Workspace::load(dir.path())?;
-    let as_of = workspace.resolve_as_of(AS_OF)?;
-    Session::new(workspace, as_of, SEQUENTIAL)?.value_portfolio(SEQUENTIAL)
-}
-
-fn valuation<'a>(report: &'a PortfolioValuation, id: &str) -> Result<&'a Valuation> {
-    report
-        .results
-        .iter()
-        .find(|(trade, _)| trade.as_str() == id)
-        .ok_or_else(|| Error::Internal(format!("missing result for {id}")))?
-        .1
-        .as_ref()
-        .map_err(Clone::clone)
+    session(dir.path(), AS_OF)?.value_portfolio(SEQUENTIAL)
 }
 
 fn usd(report: &PortfolioValuation, id: &str) -> Result<f64> {
@@ -247,7 +193,7 @@ fn workspace_methods_table_selects_linear_for_rate_futures() -> Result<()> {
         assert_eq!(valuation(&linear, id)?.method, Method::Linear, "{id}");
         let pv = usd(&linear, id)?;
         assert!(pv.abs() > 1.0, "{id}: {pv}");
-        assert_close!(pv, usd(&analytic, id)?, REL, USD_ABS);
+        assert_close!(pv, usd(&analytic, id)?, PV_REL, USD_ABS);
     }
     Ok(())
 }
@@ -262,7 +208,7 @@ fn amended_rate_future_trade_uses_the_linear_method() -> Result<()> {
     assert_close!(
         half.pv.get(Currency::USD),
         0.5 * usd(&report, "FUT-PLAIN")?,
-        REL,
+        PV_REL,
         USD_ABS
     );
     Ok(())
diff --git a/crates/rhodium/tests/support/mod.rs b/crates/rhodium/tests/support/mod.rs
new file mode 100644
index 0000000000000000000000000000000000000000..932026ea03afa629494e605ed551ed662aa51f8b
--- /dev/null
+++ b/crates/rhodium/tests/support/mod.rs
@@ -0,0 +1,78 @@
+//! ファサードの統合テストが共有する補助。
+//!
+//! 各テストのファイルは `pub mod support;` で取り込む。補助を `pub` にしておけば、取り込んだファイルが
+//! 使わない補助で `dead_code` が出ない。サンプルの複製などファイルの補助は `rhodium_fixtures::files` に、
+//! 数値の比較は `rhodium-testkit` に置く。
+
+use std::{collections::BTreeMap, path::Path};
+
+use rhodium::{
+    prelude::ExecutionContext,
+    session::{PortfolioValuation, Session, Valuation},
+    workspace::Workspace,
+};
+use rhodium_core::{Error, Result, id::TradeId};
+use rhodium_testkit::SEQUENTIAL;
+
+/// `dir` のワークスペースを読み、`as_of`（オフセットか地域名つきの日時）に逐次に評価するセッション。
+///
+/// # Errors
+///
+/// ワークスペースを読めないか、評価の瞬間を解決できないか、セッションを作れない場合。
+pub fn session(dir: &Path, as_of: &str) -> Result<Session> {
+    let workspace = Workspace::load(dir)?;
+    let as_of = workspace.resolve_as_of(as_of)?;
+    Session::new(workspace, as_of, SEQUENTIAL)
+}
+
+/// `dir` のワークスペースを読み、日付 `date` を価格設定の締め時刻で解決した瞬間に、`execution` で
+/// 準備するセッション。
+///
+/// # Errors
+///
+/// ワークスペースを読めないか、日付を締め時刻で解決できないか、セッションを作れない場合。
+pub fn session_on(dir: &Path, date: &str, execution: ExecutionContext<'_>) -> Result<Session> {
+    let workspace = Workspace::load(dir)?;
+    let as_of = workspace.resolve_date(date.parse()?)?;
+    Session::new(workspace, as_of, execution)
+}
+
+/// ポートフォリオの評価の中の、取引 `id` の行（成功した評価か、その取引の失敗）。
+///
+/// # Errors
+///
+/// `id` が取引の ID として解析できないか、結果に `id` の行が無い場合。
+pub fn row<'a>(report: &'a PortfolioValuation, id: &str) -> Result<&'a Result<Valuation>> {
+    let id: TradeId = id.parse()?;
+    report
+        .results
+        .iter()
+        .find_map(|(trade, result)| (*trade == id).then_some(result))
+        .ok_or_else(|| Error::Internal(format!("missing result for {id}")))
+}
+
+/// 取引 `id` の成功した評価。
+///
+/// # Errors
+///
+/// [`row`] が失敗した場合と、その取引の評価が失敗した場合（取引のエラーを返す）。
+pub fn valuation<'a>(report: &'a PortfolioValuation, id: &str) -> Result<&'a Valuation> {
+    row(report, id)?.as_ref().map_err(Clone::clone)
+}
+
+/// ポートフォリオを `execution` で評価し、取引ごとの報告通貨の PV を取引の ID の順に返す。
+///
+/// # Errors
+///
+/// 評価が失敗した場合と、評価できない取引がある場合（その取引のエラーを返す）。
+pub fn reporting_amounts(
+    session: &Session,
+    execution: ExecutionContext<'_>,
+) -> Result<BTreeMap<TradeId, f64>> {
+    session
+        .value_portfolio(execution)?
+        .results
+        .into_iter()
+        .map(|(id, value)| Ok((id, value?.reporting.amount())))
+        .collect()
+}
diff --git a/crates/rhodium/tests/workspace_diagnostics.rs b/crates/rhodium/tests/workspace_diagnostics.rs
new file mode 100644
index 0000000000000000000000000000000000000000..ebec00054a915cff1fe6511e07ae88613ebf53bd
--- /dev/null
+++ b/crates/rhodium/tests/workspace_diagnostics.rs
@@ -0,0 +1,220 @@
+//! ワークスペースの読み込みの診断と、評価の瞬間の解決をファサードの公開 API で確かめる。
+//!
+//! どれも `rh` を起動せず `Workspace` と `Session` だけを使うので、CLI のテストから移した。
+//! 入力ファイルの不正は、ファイルと行（または参照元の項目）を付けた `InvalidData` の診断として
+//! 一度にすべて報告する。
+
+use std::{fs, num::NonZeroUsize, path::Path};
+
+use rhodium::{
+    prelude::{ExecutionContext, ThreadExecutor},
+    session::Session,
+    workspace::Workspace,
+};
+use rhodium_core::Error;
+use rhodium_fixtures::files::{copy_dir, example};
+use tempfile::TempDir;
+
+type TestResult<T = ()> = Result<T, Box<dyn std::error::Error>>;
+
+fn workspace_copy() -> TestResult<TempDir> {
+    let dir = tempfile::tempdir()?;
+    copy_dir(&example("m0-ois"), dir.path())?;
+    Ok(dir)
+}
+
+/// `dir` のワークスペースの読み込みが返す `InvalidData` の診断の文。
+fn diagnostic(dir: &Path) -> TestResult<String> {
+    load_diagnostic(Workspace::load(dir))
+}
+
+/// 読み込みの結果が `InvalidData` の診断であることを確かめ、その文を返す。他の種類のエラーや成功は失敗にする。
+fn load_diagnostic(workspace: rhodium_core::Result<Workspace>) -> TestResult<String> {
+    match workspace {
+        Err(error @ Error::InvalidData(_)) => Ok(error.to_string()),
+        Err(error) => Err(error.into()),
+        Ok(_) => Err("expected Workspace::load to fail".into()),
+    }
+}
+
+#[test]
+fn session_reuses_the_same_instant_for_price_and_risk() -> TestResult {
+    // 同じセッションから両方を呼び出し、取引 ID と通貨別の結果を確認する。
+    fn send_sync<T: Send + Sync>() {}
+    send_sync::<Session>();
+    let workspace = Workspace::load(example("m0-ois"))?;
+    let as_of = workspace.resolve_date("2026-11-24".parse()?)?;
+    assert_eq!(as_of.to_string(), "2026-11-24T08:00:00Z");
+    let execution = ExecutionContext::sequential();
+    let session = Session::new(workspace, as_of, execution)?;
+    assert_eq!(session.as_of(), as_of);
+    assert_eq!(session.calibrated_market().groups().len(), 2);
+    let prices = session.pv(execution)?;
+    assert_eq!(prices.len(), 2);
+    assert!(prices.values().all(Result::is_ok));
+    let risk = session.quote_delta(execution)?;
+    assert!(risk.failures().is_empty());
+    assert_eq!(risk.by_ccy().len(), 2);
+    Ok(())
+}
+
+#[test]
+fn workspace_reports_multiple_csv_errors_with_lines() -> TestResult {
+    // レコード長の違い、重複 ID、非有限値のすべてを一度に報告する。
+    let dir = workspace_copy()?;
+    fs::write(
+        dir.path().join("quotes/snapshot.csv"),
+        "quote_id,value\nJPY.TONA.OIS.1M,45bp\nJPY.TONA.OIS.1M,0.006\nTOO-SHORT\nUSD.SOFR.OIS.1M,NaN\n",
+    )?;
+    let error = diagnostic(dir.path())?;
+    assert!(error.contains("quotes/snapshot.csv"), "{error}");
+    assert!(
+        error.contains("line 3") && error.contains("duplicate ID"),
+        "{error}"
+    );
+    assert!(error.contains("line: 4"), "{error}");
+    assert!(
+        error.contains("line 5") && error.contains("value: expected a finite decimal"),
+        "{error}"
+    );
+    Ok(())
+}
+
+#[test]
+fn workspace_reports_multiple_toml_rows_and_duplicate_trade_ids() -> TestResult {
+    // 同じファイルの別々の不正な行と、ファイルをまたぐ取引 ID の重複を検出する。
+    let dir = workspace_copy()?;
+    let path = dir.path().join("trades/swaps.toml");
+    let later = dir.path().join("trades/zz-duplicate.toml");
+    let good = fs::read_to_string(&path)?;
+    fs::write(&later, "]]")?;
+    fs::write(
+        &path,
+        good.replace("book = \"USD-SWAPS\"", "bok = \"USD-SWAPS\"")
+            .replace(
+                "margin_set = \"JPY-CASH\"",
+                "margin_set = \"NO-SUCH-MARGIN\"",
+            ),
+    )?;
+    let error = diagnostic(dir.path())?;
+    // 後の文書が壊れていても、先に読み込んだ文書の行エラーを元の順で返す。
+    let executor = ThreadExecutor::new(NonZeroUsize::MIN.saturating_add(3));
+    let parallel = || {
+        load_diagnostic(Workspace::load_with_execution(
+            dir.path(),
+            ExecutionContext::new(&executor),
+        ))
+    };
+    assert_eq!(parallel()?, error);
+    assert!(error.contains("trades/swaps.toml:3:1: trade[0]"), "{error}");
+    assert!(error.contains("unknown field `bok`"), "{error}");
+    assert!(
+        error.contains("trades/swaps.toml:18:1: trade[1]"),
+        "{error}"
+    );
+    assert!(
+        error.contains("margin_set: unknown margin set NO-SUCH-MARGIN"),
+        "{error}"
+    );
+    fs::write(&path, &good)?;
+    fs::write(later, good)?;
+    let error = diagnostic(dir.path())?;
+    assert_eq!(parallel()?, error);
+    assert!(
+        error.contains("zz-duplicate.toml") && error.contains("duplicate ID JPY-SWAP"),
+        "{error}"
+    );
+    Ok(())
+}
+
+#[test]
+fn workspace_rejects_invalid_version_and_fixing_header() -> TestResult {
+    // 版を解釈できない文書と、日次指標に日中の列を書いたファイルを拒否する。
+    let dir = workspace_copy()?;
+    let path = dir.path().join("agreements.toml");
+    let good = fs::read_to_string(&path)?;
+    fs::write(&path, good.replace("version = 1", "version = 2"))?;
+    let error = diagnostic(dir.path())?;
+    assert!(
+        error.contains("agreements.toml") && error.contains("unsupported version 2"),
+        "{error}"
+    );
+    fs::write(&path, good)?;
+    fs::write(
+        dir.path().join("fixings/JPY-TONA.csv"),
+        "at,value\n2026-11-20T08:00:00Z,0.0046\n",
+    )?;
+    let error = diagnostic(dir.path())?;
+    assert!(
+        error.contains("fixings/JPY-TONA.csv") && error.contains("expected CSV header date,value"),
+        "{error}"
+    );
+    Ok(())
+}
+
+#[test]
+fn workspace_reports_reference_errors_by_source_entry() -> TestResult {
+    // 別々の指標が参照する未知カレンダーを、入力ファイルと指標 ID で報告する。
+    let dir = workspace_copy()?;
+    let path = dir.path().join("refdata/indices.toml");
+    let indices = fs::read_to_string(&path)?
+        .replace("calendar = \"JPTO\"", "calendar = \"UNKNOWN-JP\"")
+        .replace("calendar = \"USNY\"", "calendar = \"UNKNOWN-US\"");
+    fs::write(path, indices)?;
+    let error = diagnostic(dir.path())?;
+    assert!(
+        error.contains("refdata/indices.toml: index[JPY-TONA]") && error.contains("UNKNOWN-JP"),
+        "{error}"
+    );
+    assert!(
+        error.contains("refdata/indices.toml: index[USD-SOFR]") && error.contains("UNKNOWN-US"),
+        "{error}"
+    );
+    Ok(())
+}
+
+#[test]
+fn workspace_rejects_ambiguous_or_missing_local_times() -> TestResult {
+    // 同じ地域の重複時刻は明示オフセットで区別し、欠落時刻は拒否する。
+    let workspace = Workspace::load(example("m0-ois"))?;
+    let instant = workspace.resolve_as_of("2026-11-24T17:00:00[Asia/Tokyo]")?;
+    assert_eq!(
+        instant,
+        workspace.resolve_as_of("2026-11-24T17:00:00+09:00[Asia/Tokyo]")?
+    );
+    // 不正な時刻は、どれも入力の誤りとして理由を付けて拒否する。地域の規則と違うオフセット、
+    // オフセットも地域も無い日時、春の欠落した時刻、秋の重複した時刻をオフセットなしで書いたもの。
+    for (text, reason) in [
+        (
+            "2026-11-24T17:00:00+08:00[Asia/Tokyo]",
+            "UTC offset does not match zone rules",
+        ),
+        ("2026-11-24T17:00:00", "invalid date/time"),
+        (
+            "2026-03-08T02:30:00[America/New_York]",
+            "local time is in a zone transition gap",
+        ),
+        (
+            "2026-11-01T01:30:00[America/New_York]",
+            "ambiguous local time requires an explicit UTC offset",
+        ),
+    ] {
+        assert!(
+            matches!(
+                workspace.resolve_as_of(text),
+                Err(Error::InvalidInput(message)) if message.contains(reason)
+            ),
+            "{text}"
+        );
+    }
+    let early = workspace.resolve_as_of("2026-11-01T01:30:00-04:00[America/New_York]")?;
+    let late = workspace.resolve_as_of("2026-11-01T01:30:00-05:00[America/New_York]")?;
+    assert_eq!(early.to_string(), "2026-11-01T05:30:00Z");
+    assert_eq!(late.to_string(), "2026-11-01T06:30:00Z");
+    // 参照データに無い地域は、他の地域の参照と同じく参照の欠けとして扱う。
+    assert!(matches!(
+        workspace.resolve_as_of("2026-11-24T17:00:00[Europe/Paris]"),
+        Err(Error::InvalidData(message)) if message.contains("Europe/Paris")
+    ));
+    Ok(())
+}
diff --git a/docs/benchmarks/issue-180.md b/docs/benchmarks/issue-180.md
index 54811159f686862bafd64a1aea00a5603169548e..e870d13d350db23c6d3d8073226547ae8516cba2 100644
--- a/docs/benchmarks/issue-180.md
+++ b/docs/benchmarks/issue-180.md
@@ -35,7 +35,7 @@ Criterion 0.8.2 で 10 標本を採り、中央の推定値と 95% 信頼区間
 
 ## ワークスペースの読み込みと準備
 
-[M0 と同じ生成器](../../crates/rhodium/tests/support/desk.rs) で JPY と USD のスワップを作りました。
+[M0 と同じ生成器](../../crates/rhodium-fixtures/src/desk.rs) で JPY と USD のスワップを作りました。
 入力の生成を先に済ませ、読み込みのみと、読み込み後に校正・準備まで行う計測を別々のプロセスで 1 回ずつ実行しました。
 時間はプロセス内の各工程で測り、最大 RSS は GNU time の `%M` で取得しました。
 変更前の読み込みは逐次で、ワーカー数は準備に適用されます。
diff --git a/docs/benchmarks/m0.md b/docs/benchmarks/m0.md
index 8f4d5fca19ed36dcbd0cd68ea097aa5c7ec4278f..e26f7897c0b54c86b387cc65d5f082225a277dba 100644
--- a/docs/benchmarks/m0.md
+++ b/docs/benchmarks/m0.md
@@ -20,7 +20,7 @@
 
 ## ワークロード
 
-入力は `crates/rhodium/tests/support/desk.rs` が一時ディレクトリに書く合成のデスクで、M0-9 のテストと同じ生成器です。
+入力は `crates/rhodium-fixtures/src/desk.rs` が一時ディレクトリに書く合成のデスクで、M0-9 のテストと同じ生成器です。
 
 - 曲線: JPY-TONA と USD-SOFR の OIS 曲線。各 15 ピラー（1M、3M、6M、1Y〜5Y、7Y、10Y、12Y、15Y、20Y、25Y、30Y）を、同じテナーの OIS の quote（計 30 本）で校正します。
 - 参照データ: 週末だけのカレンダー（2024〜2060 年）、東京と夏時間のあるニューヨークの地域、約 13 か月の TONA と SOFR の fixing。OIS の慣行はスポット 2 営業日、固定脚は年払い、SOFR の支払いは 2 営業日後です。
diff --git a/docs/design/03-architecture.md b/docs/design/03-architecture.md
index be4be108e23454cea928f2f4ddca7f8443c631ae..9f5b165849ddbb135c0c5582b58817e1361cea74 100644
--- a/docs/design/03-architecture.md
+++ b/docs/design/03-architecture.md
@@ -1,6 +1,6 @@
 # アーキテクチャ
 
-> rhodium は、ライブラリの 8 つの crate と開発用の `rhodium-testkit` からなる Cargo ワークスペースです。
+> rhodium は、ライブラリの 8 つの crate と、開発専用の `rhodium-testkit` と `rhodium-fixtures` からなる Cargo ワークスペースです。
 > 商品（`rhodium-products`）と市場（`rhodium-market`）は互いを知らず、評価（`rhodium-pricing`）が両方を結びます。
 > この文書は crate とモジュール、主要な型の定義元、型の規約とブランド、数値型と AD、エラー、依存と lint、並列、性能の規則、拡張の仕方を定めます。
 
@@ -16,6 +16,8 @@
 | `rhodium-xva` | 時間グリッド、パス上の価値、担保と close-out、IM、XVA の指標、XVA の感応度、配分、incremental XVA | risk まで |
 | `rhodium` | ファサード: `prelude`、`Workspace`（ファイルの読み込み）、`Session`（高水準 API） | すべて、toml、csv、serde_json |
 | `rhodium-cli` | `rh` コマンド | rhodium |
+| `rhodium-testkit`（開発専用） | 閉包と数値だけを受け取るテストの補助（数値の照合、許容誤差の定数、AAD とバンプの照合、型テストの対、決定的な時計） | scioxide-core、serde、serde_json、tempfile |
+| `rhodium-fixtures`（開発専用） | 複数の crate のテスト、ベンチマーク、例が共有する市場・参照データ・商品・ワークスペースの fixture と QuantLib の値 | core、products、market、pricing、serde、serde_json、toml |
 
 ```mermaid
 flowchart BT
@@ -37,6 +39,10 @@ flowchart BT
 - 商品は market を知らず、market は商品を知りません。
   校正は「校正商品を取引と同じ標準商品とプライサーで評価する」ため、市場の定義（構造と校正の指定）と校正は `rhodium-pricing` に置きます。
 - ファイルの読み込み（TOML、CSV、JSON）はファサードと CLI だけが行い、コア crate は I/O を持ちません。
+- 開発専用の二つの crate は `publish = false` で、ライブラリの crate からは `[dev-dependencies]` でだけ使います（図には描きません）。
+  `rhodium-fixtures` は pricing に依存するので、pricing の単体テスト（`src` の `#[cfg(test)]`）からは使いません。
+  dev-dependency の循環で `rhodium_pricing` が二つになり、型が合わないためです。
+  役割の分け方と共有の置き場所は [品質のテストの補助](11-quality.md#テストの補助) にあります。
 
 ## モジュールの地図
 
@@ -223,8 +229,8 @@ impl<R: Real> Market<R> {
   toml、csv、serde_json はファサードと CLI だけが使い、引数の解析の clap は CLI だけが使います。
   開発用には proptest（性質テスト）、criterion（ベンチマーク、`default-features = false`）、tempfile（ファイルのワークスペースのテストと、型テストのケースごとの出力先）を使い、テストの外部 fixture の読み込みには serde_json と toml も使います。
   testkit の型テストのランナーは、rustc と Cargo の JSON の出力と、期待する診断の `case.json` の解析に serde_json を使います。
-- 開発専用の crate `rhodium-testkit` が、テストの補助を一か所に持ちます（一覧と、まだ無い補助の担当は [品質のテストの補助](11-quality.md#テストの補助)）。
-  ライブラリの crate が testkit に依存するときは dev-dependency だけにします（テストを持たない間の rhodium-xva は依存しません）。
+- 開発専用の crate `rhodium-testkit` が数値と閉包だけを受け取るテストの補助を、`rhodium-fixtures` が rhodium の型で書いた共有の fixture を持ちます（一覧と、まだ無い補助の担当は [品質のテストの補助](11-quality.md#テストの補助)）。
+  ライブラリの crate がこれらに依存するときは dev-dependency だけにします（テストを持たない間の rhodium-xva は依存しません）。
 - feature は加算的にし、feature の有無で計算結果を変えません。
 - 公開する値の型は `Send + Sync` です（`Var<'t, _>` とブランド付きの値だけが例外で、スコープの中に閉じます）。
 
@@ -238,7 +244,7 @@ impl<R: Real> Market<R> {
   この lint は `matches!`、let-else、タプルの `match` を見ないので、複数の variant にわたる分類はそれらで書かず、enum のメソッドの網羅した `match` に置きます。
 - `clippy.toml` の disallowed は `std::time::SystemTime::now`、`std::time::Instant::now`、`chrono::Utc::now`、`chrono::Local::now`、`std::env` の読み出しだけにします（[P1](02-principles.md#p1-純粋関数と明示的な入力)）。
 - 型や lint で止められない実装の誤りの形（ID の文字列の合成、許容誤差のリテラル、下流での参照データの照合し直しなど）は、`.githooks` の禁止パターンの表で止めます（規則の書き方と diff・tree のモードは[品質の CI](11-quality.md#ci)）。
-- lint の例外は reason 付きの `#[expect]` にし、効く範囲を最小の item に絞ります。`#[allow]` は、例外が要らなくなっても `unfulfilled_lint_expectations` で気づけないので避けます。
+- lint の例外は reason 付きの `#[expect]` にし、効く範囲を最小の item に絞ります。`#[allow]` は、例外が要らなくなっても `unfulfilled_lint_expectations` で気づけないので使わず、`clippy::allow_attributes` を deny にして書けなくします。
 - `missing_docs` は公開 API に warn とし、doc のコメントと意図を説明するコメントは日本語で書きます。
 
 ## 並列と再現性
diff --git a/docs/design/11-quality.md b/docs/design/11-quality.md
index 2c936c10eb342084308b8191ad9f66573bc6bfe1..82a4ae6cd3fe87f918c82160502dcb25a16ea134 100644
--- a/docs/design/11-quality.md
+++ b/docs/design/11-quality.md
@@ -14,7 +14,7 @@
 | 性質テスト | 法則（スケジュールの不変条件、カレンダーの合成、`MoneyBag` の和、随伴の線形性） | proptest |
 | compile-fail テスト | 型で禁止した演算がコンパイルできず、正しい使い方はコンパイルできる | testkit の `_ok` と失敗の対 |
 | 一致のテスト | 同じ商品の二つの評価（IR のリニア評価と解析プライサー、MC と閉形式）が一致する | [オラクル](#オラクル) |
-| 外部の fixture | QuantLib、ORE などの値と一致する | 事前に生成した値を `tests/fixtures/` に置く |
+| 外部の fixture | QuantLib、ORE などの値と一致する | 事前に生成した値を、一つの crate のテストだけが読むなら `tests/fixtures/` に、複数の crate のテストが読むなら `rhodium-fixtures` に置く |
 | AAD とバンプ | AAD の感応度が、再校正込みのバンプの差分に一致する | `check_aad_vs_bump` |
 | MC の推定量 | MC の値が閉形式か別の推定量と標準誤差の範囲で一致する | `assert_mc!` |
 | golden | 結果が以前の値から変わっていない | `tests/golden/` の CSV |
@@ -23,6 +23,16 @@
   テストを ID で登録する仕組みは持ちません。
   受け入れ基準の ID（`M0-3` など）は [ロードマップ](12-roadmap.md) と issue・PR の本文に書き、どのテストで確かめたかをそこに記します。
 - テストは、何を確かめ、なぜその値になるかを日本語のコメントで書き、乱数は固定の seed を使います。
+- テストは確かめる振る舞いの単位でファイルに分け、マイルストーンの単位でまとめません（例: pricing の `rate_futures.rs`、`fra.rs`、`fx_fixing.rs`）。
+  同じシナリオは一か所で確かめます。
+  スレッド数によらないことはファサードの `parallel.rs` が値で確かめ、CLI のテストは `--threads` で同じ形の表を出すことだけを確かめます。
+- 公開 API だけを使うテストは crate の `tests/` の統合テストに置き、`src` の `#[cfg(test)]` には非公開の構造（IR のノードの共有、校正の内部の状態など）を見るテストだけを置きます。
+  統合テストは利用者と同じ入口を通るので、内部の書き換えでテストを直す必要が無く、公開 API の不足にも気づけます。
+- 判定は、テストが落ちるべきときに落ちる強さで書きます。
+  - エラーは `is_err()` でなく種類まで確かめ（`matches!(result, Err(Error::InvalidData(message)) if message.contains(..))`）、診断の文は理由の部分を照合します。
+  - 集合と表は全体で比べます。取引ごとの結果は両側の取引 ID の列全体と各取引の値を比べ（`zip` だけでは片側の末尾の余分な結果を見逃す）、警告は一部の取引の有無でなく、全取引の警告の集合を期待の表と比べます。
+  - fixture の書き換えは `get_mut(..).expect(..)` などで、キーが無ければ失敗させます。`if let` で黙って飛ばすと、fixture のキーの名前が変わったときに既定のシナリオのまま通ります。
+  - 判定を厳しくしたときは、期待値を一時的に変えて落ちることを確かめ、その手順を PR に書きます。
 - golden の値を更新するときは、変わった理由をコミットのメッセージに書きます。
 - 外部の fixture は、出典（ライブラリと版、公表資料）、生成のスクリプト、生成の設定、ライセンスを、fixture のディレクトリの `README.md` に書きます。
   テストの実行時に外部のライブラリを呼びません。
@@ -30,16 +40,48 @@
 
 ### テストの補助
 
-テストの補助は開発専用の crate `rhodium-testkit`（dev-dependency だけで使い、公開しない）に置き、閉包と数値だけを受け取るため、どの crate のテストからも使えます。
+テストの補助は、開発専用の二つの crate（`publish = false`、dev-dependency だけで使う）と、各 crate の `tests/support/` に置きます。
+置き場所は、補助が受け取る型と、使う crate の範囲で決めます。
+
+| 置き場所 | 置くもの | 使う側 |
+|---|---|---|
+| `rhodium-testkit` | 閉包と数値だけを受け取る補助と、許容誤差の定数。金融の crate に依存しない | どの crate のテストからも |
+| `rhodium-fixtures` | 複数の crate のテストが使い、core・products・market・pricing の型だけで書ける市場データ、参照データ、商品、ワークスペースの組み立て、QuantLib の値（`quantlib/`、出典は `quantlib/README.md`） | pricing、risk、ファサード、CLI の統合テスト、ベンチマーク、例 |
+| 各 crate の `tests/support/mod.rs` | その crate の統合テストだけが使う補助と、risk やファサードの型を扱う補助 | 同じ crate の統合テスト |
+| テストのファイルの中 | 一つのファイルだけが使う補助 | そのファイル |
+
+- 他の crate やディレクトリのファイルを `#[path]` で取り込みません（禁止パターン `P-path`）。
+  取り込みは Cargo の依存グラフに現れず、[CI](#ci) の逆依存の検査で利用側を追えないうえ、使わない関数ごとに `dead_code` の例外が要るためです。
+- `tests/support/mod.rs` は各テストのファイルで `pub mod support;` と宣言し、補助を `pub` にします。
+  テストの crate の `pub` の item は使われなくても `dead_code` にならないので、ファイルごとに使う補助が違っても例外が要りません。
+- `rhodium-fixtures` は pricing に依存するので、pricing の単体テストからは使いません（dev-dependency の循環で `rhodium_pricing` の型が二つになる）。
+  非公開の状態を書き換える単体テストは、要る最小の入力をそのテストの中で組みます。
+- 公開の関数にするのは、テストから実際に使うものだけです。
+
+`rhodium-testkit` の補助:
 
 | 補助 | 役割 |
 |---|---|
-| `assert_close!(x, y, rel, abs)`、`assert_close!(x, y, rel, abs, scale)` | `abs(x − y) ≤ max(rel·M, abs)` を確かめる。M は `max(abs(x), abs(y))` で、5 引数の形は和の項の絶対値の和などの `scale` も M に含める（[許容誤差](#許容誤差)） |
+| `assert_close!(x, y, rel, abs)`、`assert_close!(x, y, rel, abs, scale)` | `abs(x − y) ≤ max(rel·M, abs)` を確かめる。M は `max(abs(x), abs(y))` で、5 引数の形は和の項の絶対値の和などの `scale` も M に含める（[許容誤差](#許容誤差)）。末尾に `; "{label}: "` を足すと、失敗の表示の先頭に表の行などの文脈を付ける |
+| `PV_REL`、`MARKET_REL`、`MARKET_ABS`、`SENSITIVITY_REL`、`IMPLIED_REL`、`IMPLIED_ABS`、`SECOND_ORDER_REL`、`minor_unit_abs(digits)` | [許容誤差](#許容誤差)の表の ε と a。`minor_unit_abs` は通貨の最小単位の 1e-6 倍（JPY は 1e-6、USD は 1e-8）。表を変えるときは定数も同じ PR で直す |
+| `entry_or_zero(entries, key)`、`assert_sparse_close(actual, expected, rel, abs)` | 疎なベクトル（`RiskVector::entries()` など）を、欠けたキーを 0 として読み、両側のキーの和集合で比べる |
+| `SEQUENTIAL` | 逐次の実行文脈。並列の結果と比べる基準にも使う |
 | `assert_mc!(estimate, se, reference, k, a)` | `abs(estimate − reference) ≤ k·se + a` を確かめる。MC とともに M4 で足す |
 | `check_aad_vs_bump(reprice, aad, unit, pv_error, absolute, label)` | quote の増分（自然な単位）から再校正込みの価格を返す閉包、報告単位あたりの AAD の感応度、報告単位の幅 h₁、再校正の PV の誤差 δ_V、絶対の下限を受け取り、下の手順で照合する |
 | `compile_fail_pair!("case")` | `tests/ui/case_ok.rs` を実行し、`tests/ui/case.rs` のコンパイル診断を `case.json` と照合する |
 | `TickClock` | 読むたびに 1 ns 進む決定的な時計。`TickClock::after(n)` を期限にすると n 回目の読み出しで止まるので、区間の開始後の取消と期限の伝わり方をスレッドの速さによらず確かめられる |
 
+`rhodium-fixtures` のモジュール:
+
+| モジュール | 内容 |
+|---|---|
+| `ois`、`multicurve`、`steps` | QuantLib で校正した TONA と SOFR の OIS 曲線、ターム指標の複数曲線、会合日の階段と年末の跳びの TONA（R1〜R5） |
+| `xccy` | FX のポイントと MtM XCCY で校正する JPY 担保の USD 曲線（R6） |
+| `futures` | 金利先物の清算指標と限月、LGM のマークを持つ USD の市場 |
+| `desk` | M0-9 とベンチマークの合成のデスクのワークスペース |
+| `files` | リポジトリの `examples/` のサンプルの場所、複製、追記 |
+| 直下 | `COLLATERAL`、`bumped`、`solved`、`solved_within`、`prepare_standard` |
+
 ## オラクル
 
 正しさの根拠には、評価のコードと独立したものを使います。
@@ -283,7 +325,7 @@ M = max(|x|, |y|, Σ|項|)
   型テストの登録や `tests/ui/` の追加・変更・削除があれば、testkit の登録検査も実行します。
 - 型や lint で止められない実装の誤りの形は、`.githooks/checks.py` の禁止パターンの表（`RULES`）で止めます。
   各行は、名前（`P-…`）、Rust の行に当てる正規表現、対象と例外のパスの glob、守る理由、規則の定義元の設計書の節を持ちます。
-  ライブラリの規則はテストのファイル（`tests/`、`src/**/tests.rs`、`src/**/tests/`）と `#[cfg(test)]` の item を照合せず、コメントの中も照合しません。
+  ライブラリの規則はテストのファイル（`tests/`、`src/**/tests.rs`、`src/**/tests/`）、`#[cfg(test)]` の item と、開発専用の crate（`rhodium-testkit`、`rhodium-fixtures`）を照合せず、コメントの中も照合しません。
   `requires_nearby` を持つ行は、当たった行とその直前の数行に決まった注記（`P-internal` の `// 不変条件:`）があれば違反にしません。
   行のモードは二つです。
   `diff` の行は、commit と push で追加した行だけを照合するので、既存の違反を残したまま規則を入れられます。
@@ -294,7 +336,9 @@ M = max(|x|, |y|, Σ|項|)
   規則を足すときは、まず `audit` で既存の違反を数え、`diff` で入れる理由か、`tree` で入れられることを PR に書きます。
   違反に当たったときは、検査を避ける書き換えではなく、規則の定義元の設計書に従って直します。
 - マイルストーンの完了確認では、ワークスペース全体で `cargo fmt --all --check`、`cargo clippy --locked --workspace --all-targets --all-features -- -D warnings`、`cargo test --locked --workspace --all-features`、`RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --all-features --no-deps` を実行し、結果を完了確認の issue に記録します。
-- 大きな MC、外部 fixture の全件、XVA の規模のテストは `#[ignore = "理由"]` を付けた遅いテストにし、`cargo test --workspace --all-features -- --ignored` で定期的に実行します。
+- 一つのテストは数秒で終えます。
+  中心差分の照合は、同じ (quote, 幅) の再校正とその市場での取引の準備を、取引の間で共有します（risk のテストの `support::Recalibrations`）。
+  共有しても数秒を超える大きな MC、外部 fixture の全件、XVA の規模のテストは `#[ignore = "理由"]` を付けた遅いテストにし、`cargo test --workspace --all-features -- --ignored` で定期的に実行します。
 - ベンチマークはテストと分け、[入出力と実行](10-io-and-execution.md) の手順で測ります。
 
 ## 拡張のチェックリスト
@@ -338,3 +382,12 @@ M = max(|x|, |y|, Σ|項|)
 - 判断: 型とブランドで禁止した演算は、`_ok` の例と対にした compile-fail テストで確かめます。
 - 理由: 型の禁止は、便利な演算子を足すと黙って外れますが、テストがあれば気づけます。
 - 採らない案: 型の禁止を文書だけで示す案は、退行を検出できないため採りません。
+
+### D5 共有の fixture は開発専用の crate に置く
+
+- 判断: 複数の crate のテストが使う市場や商品の fixture は、開発専用の crate `rhodium-fixtures` に置き、利用側の dev-dependency にします。
+  数値と閉包だけの補助は `rhodium-testkit` に残し、testkit は金融の crate に依存させません。
+- 理由: Cargo の依存として宣言すれば、逆依存の検査が利用側のテストを追え、使わない関数の `dead_code` の例外も要りません。
+  testkit を金融の crate に依存させないことで、core や market のテストからも testkit を使えます。
+- 採らない案: pricing の `tests/fixtures/*.rs` を他の crate が `#[path]` で取り込む案は、依存グラフに現れず、取り込むファイルごとに例外が要るため採りません。
+  fixture を testkit に置く案は、testkit が pricing に依存し、core と market のテストが testkit を使えなくなるため採りません。
