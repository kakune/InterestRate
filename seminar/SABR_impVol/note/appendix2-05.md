diff --git a/crates/rhodium/src/workspace.rs b/crates/rhodium/src/workspace.rs
index a1dc744a3ef568d1683807087411e23c968a01e3..2dfba2174b540131196eee4e15377d37a2469ca3 100644
--- a/crates/rhodium/src/workspace.rs
+++ b/crates/rhodium/src/workspace.rs
@@ -19,9 +19,12 @@ use rhodium_core::{
 };
 use rhodium_market::quote::QuoteSet;
 use rhodium_pricing::{
-    config::PricingConfig, definition::MarketDefinition, resolve::CollateralKey,
+    calibration::calibrate, config::PricingConfig, context::Context, definition::MarketDefinition,
+    resolve::CollateralKey,
 };
-use rhodium_products::ir::Product;
+use rhodium_products::trade::ResolvedTrade;
+
+pub use rhodium_pricing::calibration::CalibratedMarket;
 
 use self::{input::Source, trade_files::LoadedTrade};
 
@@ -85,6 +88,24 @@ impl Workspace {
         })
     }
 
+    /// 指定の瞬間の市場を校正する。
+    ///
+    /// # Errors
+    ///
+    /// 評価時点と入力の組から市場構造または校正条件を組み立てられない場合。
+    pub fn calibrate(&self, as_of: Instant) -> Result<CalibratedMarket> {
+        calibrate(&self.inputs.market, self.quotes(), &self.context(as_of))
+    }
+
+    pub(crate) fn context(&self, as_of: Instant) -> Context<'_> {
+        Context {
+            as_of,
+            refdata: &self.inputs.refdata,
+            fixings: &self.inputs.fixings,
+            config: &self.inputs.config,
+        }
+    }
+
     /// このワークスペースの quote スナップショット。
     #[must_use]
     pub const fn quotes(&self) -> &QuoteSet {
@@ -127,14 +148,19 @@ impl Workspace {
     /// 取引のブックを返す。指定の取引が無ければ `None`。
     #[must_use]
     pub fn book(&self, trade: &TradeId) -> Option<&BookId> {
-        self.inputs.trades.get(trade).map(|trade| &trade.book)
+        self.inputs
+            .trades
+            .get(trade)
+            .map(|loaded| &loaded.resolved.trade().book)
     }
 
-    pub(crate) fn trades(&self) -> impl Iterator<Item = (&TradeId, &Product, &CollateralKey)> {
+    pub(crate) fn trades(
+        &self,
+    ) -> impl Iterator<Item = (&TradeId, &ResolvedTrade, &CollateralKey)> {
         self.inputs
             .trades
             .iter()
-            .map(|(id, trade)| (id, &trade.product, &trade.collateral))
+            .map(|(id, loaded)| (id, &loaded.resolved, &loaded.collateral))
     }
 }
 
diff --git a/crates/rhodium/src/workspace/csv_files.rs b/crates/rhodium/src/workspace/csv_files.rs
index 6d30e0a3cf0dd34c3f4bb3f67c2d5fd927324981..8af50540b85c403a02fc2b12bc083850e87fd946 100644
--- a/crates/rhodium/src/workspace/csv_files.rs
+++ b/crates/rhodium/src/workspace/csv_files.rs
@@ -124,7 +124,7 @@ fn rows(
     finish(errors)
 }
 
-fn parse_value(text: &str) -> Result<f64> {
+pub(super) fn parse_value(text: &str) -> Result<f64> {
     // 接尾辞付きの値も十進で尺度を変えてから一度だけ f64 に丸める。
     let scaled;
     let decimal = if let Some(number) = text.strip_suffix("bp") {
diff --git a/crates/rhodium/src/workspace/input.rs b/crates/rhodium/src/workspace/input.rs
index c9cb02a09af3cca1819fe39ea93865fce532b46a..528eb3a40ee09ea84f5e8d141bc6e49456769a81 100644
--- a/crates/rhodium/src/workspace/input.rs
+++ b/crates/rhodium/src/workspace/input.rs
@@ -27,7 +27,7 @@ impl Source {
         toml::from_str(&self.text).map_err(|error| self.invalid(error))
     }
 
-    // 文書の版を照合する。M0 が解釈できるのは版 1 だけ。
+    // 文書の版を照合する。1.0 までの文書の版は 1 だけ（設計書「ファイルの形式」）。
     pub(super) fn version(&self, version: Spanned<u32>) -> Result<()> {
         let start = version.span().start;
         let number = version.into_inner();
@@ -103,6 +103,10 @@ pub(super) fn invalid(path: &Path, error: impl Display) -> Error {
 }
 
 // 集めたエラーを改行区切りの一つの InvalidData にする。InvalidData の接頭辞は外側で一度だけ付ける。
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#198 で各 variant を網羅する"
+)]
 pub(super) fn finish(errors: Vec<Error>) -> Result<()> {
     if errors.is_empty() {
         return Ok(());
diff --git a/crates/rhodium/src/workspace/market_file.rs b/crates/rhodium/src/workspace/market_file.rs
index 30f50615396d0b21461f368720bb0dc607155341..fb4bbab223c193a8c11769abb4ec9beee6c04ac9 100644
--- a/crates/rhodium/src/workspace/market_file.rs
+++ b/crates/rhodium/src/workspace/market_file.rs
@@ -1,7 +1,7 @@
 use std::{collections::BTreeMap, path::PathBuf};
 
 use rhodium_core::{
-    Result,
+    Error, Result,
     calendar::CalendarSpec,
     currency::Currency,
     id::CalendarId,
@@ -9,14 +9,14 @@ use rhodium_core::{
     time::{Instant, TzId},
 };
 use rhodium_market::{
-    curve::{CurveId, CurveRole, CurveSpec, Extrapolation},
+    curve::{CurveExpr, CurveId, CurveRole, CurveSpec, Extrapolation},
     quote::{QuoteId, QuoteKind},
     structure::{BlockId, ParamId},
 };
 use rhodium_pricing::{
     definition::{
         CalibrationTarget, ElementSelector, FxDefinition, InstrumentBinding, LgmDefinition,
-        MarketDefinition, QuoteTransform,
+        MarketDefinition, QuoteTransform, StepPenalty, vol::VolatilityDefinition,
     },
     resolve::CollateralKey,
 };
@@ -24,7 +24,10 @@ use rhodium_products::standard::InstrumentTemplate;
 use serde::Deserialize;
 use toml::{Spanned, Value};
 
-use super::input::{Source, finish};
+use super::{
+    csv_files,
+    input::{Source, finish},
+};
 
 // 曲線と binding は行ごとの診断を返せるよう、位置付きの値として受けてから個別に復号する。
 #[derive(Deserialize)]
@@ -38,17 +41,27 @@ struct MarketDocument {
     #[serde(default)]
     fx: Vec<Spanned<Value>>,
     #[serde(default)]
+    volatilities: Vec<Spanned<Value>>,
+    #[serde(default)]
     lgm: BTreeMap<String, Spanned<Value>>,
+    #[serde(default)]
+    step_penalties: Vec<Spanned<Value>>,
+    #[serde(default)]
+    gauss_newton: Vec<Spanned<Value>>,
 }
 
+// 曲線は明示したピラーの節点（`pillars`）か、部品の式（`expr`、`CurveExpr` の serde の表現）で書く。
 #[derive(Deserialize)]
 #[serde(deny_unknown_fields)]
 struct CurveInput {
     currency: Currency,
     roles: Vec<CurveRole>,
-    pillars: Vec<Instant>,
     #[serde(default)]
-    extrapolation: Extrapolation,
+    pillars: Option<Vec<Instant>>,
+    #[serde(default)]
+    extrapolation: Option<Extrapolation>,
+    #[serde(default)]
+    expr: Option<CurveExpr>,
 }
 
 #[derive(Deserialize)]
@@ -57,11 +70,16 @@ struct InstrumentInput {
     quote: QuoteId,
     kind: QuoteKind,
     target: Target,
-    template: InstrumentTemplate,
+    #[serde(default)]
+    template: Option<InstrumentTemplate>,
     #[serde(default)]
     transform: QuoteTransform,
-    collateral: CollateralKey,
+    #[serde(default)]
+    collateral: Option<CollateralKey>,
     determines: Determines,
+    /// 最小二乗の当てはめの許容幅（quote の単位。`"0.5bp"` などの接尾辞付きか十進）。
+    #[serde(default)]
+    tolerance: Option<Value>,
 }
 
 #[derive(Deserialize)]
@@ -69,41 +87,150 @@ struct InstrumentInput {
 enum Target {
     ParRate,
     Price,
+    DirectMark,
 }
 
 #[derive(Deserialize)]
 #[serde(deny_unknown_fields)]
 struct Determines {
     curve: CurveId,
+    #[serde(default)]
+    part: Option<String>,
     element: Option<u32>,
 }
 
+#[derive(Deserialize)]
+#[serde(deny_unknown_fields)]
+struct BlockInput {
+    curve: CurveId,
+    part: String,
+}
+
+#[derive(Deserialize)]
+#[serde(deny_unknown_fields)]
+struct PenaltyInput {
+    curve: CurveId,
+    part: String,
+    weight: f64,
+    /// 参照の幅 s（`"1bp"` などの接尾辞付きか十進）。
+    scale: Value,
+}
+
+/// 接尾辞付き（`"0.5bp"`、`"1%"`）の文字列か十進の数を読む。
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "外部の enum は variant が増えうるので、既定の腕で受ける"
+)]
+fn scaled(value: &Value, field: &str) -> Result<f64> {
+    match value {
+        Value::String(text) => csv_files::parse_value(text),
+        Value::Float(number) => Ok(*number),
+        Value::Integer(number) => i32::try_from(*number)
+            .map(f64::from)
+            .map_err(|_| Error::InvalidData(format!("{field}: number out of range"))),
+        _ => Err(Error::InvalidData(format!(
+            "{field}: expected a number or a string with bp or %"
+        ))),
+    }
+    .map_err(|error| error.context(field))
+}
+
+/// 曲線の行を節点か部品の式の曲線にする。二つの書き方は排他で、外挿は節点の曲線か式の `nodes` に書く。
+fn curve_spec(id: &str, row: CurveInput) -> Result<CurveSpec> {
+    match (row.pillars, row.extrapolation, row.expr) {
+        (Some(pillars), extrapolation, None) => CurveSpec::new(
+            id.parse()?,
+            row.currency,
+            row.roles,
+            pillars,
+            extrapolation.unwrap_or_default(),
+        ),
+        (None, None, Some(expr)) => {
+            CurveSpec::from_expr(id.parse()?, row.currency, row.roles, expr)
+        }
+        (Some(_), _, Some(_)) => Err(Error::InvalidData("pillars and expr are exclusive".into())),
+        (None, Some(_), Some(_)) => Err(Error::InvalidData(
+            "extrapolation belongs to the nodes inside expr".into(),
+        )),
+        (None, _, None) => Err(Error::InvalidData("curve requires pillars or expr".into())),
+    }
+}
+
+/// binding の行を校正の条件にする。商品の目標はテンプレートと担保を、直接マークはどちらも取らない。
+fn binding(row: InstrumentInput) -> Result<InstrumentBinding> {
+    let block = BlockId::new(
+        row.determines.curve,
+        row.determines.part.as_deref().unwrap_or("nodes"),
+    )?;
+    let determines = match row.determines.element {
+        Some(element) => ElementSelector::Param(ParamId::new(block, element)),
+        None => ElementSelector::Block(block),
+    };
+    let (target, collateral) = match (row.target, row.template, row.collateral) {
+        (Target::ParRate, Some(template), Some(collateral)) => {
+            (CalibrationTarget::ParRate(template), collateral)
+        }
+        (Target::Price, Some(template), Some(collateral)) => {
+            (CalibrationTarget::Price(template), collateral)
+        }
+        (Target::DirectMark, None, None) => (
+            CalibrationTarget::DirectMark,
+            CollateralKey::Uncollateralized,
+        ),
+        (Target::DirectMark, _, _) => {
+            return Err(Error::InvalidData(
+                "direct_mark takes neither a template nor a collateral".into(),
+            ));
+        }
+        (Target::ParRate | Target::Price, _, _) => {
+            return Err(Error::InvalidData(
+                "par_rate and price require a template and a collateral".into(),
+            ));
+        }
+    };
+    let binding = InstrumentBinding::new(
+        row.quote,
+        row.kind,
+        target,
+        row.transform,
+        collateral,
+        determines,
+    )?;
+    match row.tolerance {
+        Some(tolerance) => binding.with_tolerance(scaled(&tolerance, "tolerance")?),
+        None => Ok(binding),
+    }
+}
+
+fn validate_valuation_references(
+    refdata: &RefData,
+    source: &Source,
+    document: &MarketDocument,
+    errors: &mut Vec<Error>,
+) {
+    if let Err(error) = refdata.calendar(&CalendarSpec::Base(document.valuation_calendar.clone())) {
+        errors.push(source.invalid(error.context("valuation_calendar")));
+    }
+    if let Err(error) = refdata.zone(&document.valuation_zone) {
+        errors.push(source.invalid(error.context("valuation_zone")));
+    }
+}
+
 // 市場定義を読み、参照データにある通貨、カレンダー、地域、慣行、担保の指標と照合する。
 // 曲線と binding の整合（未知の曲線や未割当の要素）は MarketDefinition::new が確かめる。
 pub(super) fn load(path: PathBuf, refdata: &RefData) -> Result<MarketDefinition> {
     let source = Source::read(path)?;
     let document: MarketDocument = source.decode()?;
-    source.version(document.version)?;
+    source.version(document.version.clone())?;
     let mut errors = Vec::new();
-    if let Err(error) = refdata.calendar(&CalendarSpec::Base(document.valuation_calendar.clone())) {
-        errors.push(source.invalid(error.context("valuation_calendar")));
-    }
-    if let Err(error) = refdata.zone(&document.valuation_zone) {
-        errors.push(source.invalid(error.context("valuation_zone")));
-    }
+    validate_valuation_references(refdata, &source, &document, &mut errors);
 
     let mut curves = Vec::new();
     if let Err(error) = source.entries("curves", document.curves, |id, row: CurveInput| {
         refdata
             .currency(row.currency)
             .map_err(|error| error.context("currency"))?;
-        curves.push(CurveSpec::new(
-            id.parse()?,
-            row.currency,
-            row.roles,
-            row.pillars,
-            row.extrapolation,
-        )?);
+        curves.push(curve_spec(id, row)?);
         Ok(())
     }) {
         errors.push(error);
@@ -114,39 +241,63 @@ pub(super) fn load(path: PathBuf, refdata: &RefData) -> Result<MarketDefinition>
         "instruments",
         document.instruments,
         |row: InstrumentInput| {
-            row.template.validate(refdata)?;
-            row.collateral
-                .validate_references(refdata)
-                .map_err(|error| error.context("collateral"))?;
-            let block = BlockId::nodes(row.determines.curve);
-            let determines = match row.determines.element {
-                Some(element) => ElementSelector::Param(ParamId::new(block, element)),
-                None => ElementSelector::Block(block),
-            };
-            let target = match row.target {
-                Target::ParRate => CalibrationTarget::ParRate(row.template),
-                Target::Price => CalibrationTarget::Price(row.template),
-            };
-            bindings.push(InstrumentBinding::new(
-                row.quote,
-                row.kind,
-                target,
-                row.transform,
-                row.collateral,
-                determines,
+            if let Some(template) = &row.template {
+                template.validate(refdata)?;
+            }
+            if let Some(collateral) = &row.collateral {
+                collateral
+                    .validate_references(refdata)
+                    .map_err(|error| error.context("collateral"))?;
+            }
+            bindings.push(binding(row)?);
+            Ok(())
+        },
+    ) {
+        errors.push(error);
+    }
+    let mut penalties = Vec::new();
+    if let Err(error) = source.rows(
+        "step_penalties",
+        document.step_penalties,
+        |row: PenaltyInput| {
+            penalties.push(StepPenalty::new(
+                BlockId::new(row.curve, &row.part)?,
+                row.weight,
+                scaled(&row.scale, "scale")?,
             )?);
             Ok(())
         },
     ) {
         errors.push(error);
     }
+    let mut gauss_newton = Vec::new();
+    if let Err(error) = source.rows("gauss_newton", document.gauss_newton, |row: BlockInput| {
+        gauss_newton.push(BlockId::new(row.curve, &row.part)?);
+        Ok(())
+    }) {
+        errors.push(error);
+    }
     let mut fx = Vec::new();
     if let Err(error) = source.rows("fx", document.fx, |row: FxDefinition| {
+        // スポットラグと短期の quote の連鎖、quote の向きは参照データのペアの慣行と照らす。
+        row.validate_references(refdata)?;
         fx.push(row);
         Ok(())
     }) {
         errors.push(error);
     }
+    let mut volatilities = Vec::new();
+    if let Err(error) = source.rows(
+        "volatilities",
+        document.volatilities,
+        |row: VolatilityDefinition| {
+            row.validate_references(refdata)?;
+            volatilities.push(row);
+            Ok(())
+        },
+    ) {
+        errors.push(error);
+    }
     let mut lgm = Vec::new();
     if let Err(error) = source.entries("lgm", document.lgm, |key, row: LgmDefinition| {
         let ccy = key.parse()?;
@@ -165,5 +316,8 @@ pub(super) fn load(path: PathBuf, refdata: &RefData) -> Result<MarketDefinition>
     )
     .and_then(|definition| definition.with_fx(fx))
     .and_then(|definition| definition.with_lgm(lgm))
+    .and_then(|definition| definition.with_volatilities(volatilities))
+    .and_then(|definition| definition.with_step_penalties(penalties))
+    .and_then(|definition| definition.with_gauss_newton(gauss_newton))
     .map_err(|error| source.invalid(error))
 }
diff --git a/crates/rhodium/src/workspace/refdata_files.rs b/crates/rhodium/src/workspace/refdata_files.rs
index 8efa567a791ef3a4d5147ea34166f7659c815412..c80e901ce077de585aa1fceaba14631eceb00f80 100644
--- a/crates/rhodium/src/workspace/refdata_files.rs
+++ b/crates/rhodium/src/workspace/refdata_files.rs
@@ -6,8 +6,8 @@ use rhodium_core::{
     date::{Date, Weekday},
     id::{CalendarId, ConventionId},
     refdata::{
-        Convention, CurrencyDefinition, FallbackDefinition, IndexDefinition, RefData, RefDataEntry,
-        RefDataParts, SecurityDefinition,
+        Convention, CurrencyDefinition, FallbackDefinition, FxConvention, IndexDefinition, RefData,
+        RefDataEntry, RefDataParts, SecurityDefinition,
     },
     time::{Instant, TzId, UtcOffset, ZoneRules},
 };
@@ -34,9 +34,19 @@ macro_rules! document {
 document!(Currencies, currency);
 document!(Calendars, calendar);
 document!(Zones, zone);
-document!(Conventions, convention);
 document!(Securities, security);
 
+// 商品の慣行と、通貨ペアの現物の受渡慣行（`[[fx]]`）を同じファイルに書く。
+#[derive(Deserialize)]
+#[serde(deny_unknown_fields)]
+struct Conventions {
+    version: Spanned<u32>,
+    #[serde(default)]
+    convention: Vec<Spanned<Value>>,
+    #[serde(default)]
+    fx: Vec<Spanned<Value>>,
+}
+
 #[derive(Deserialize)]
 #[serde(deny_unknown_fields)]
 struct Indices {
@@ -126,11 +136,15 @@ pub(super) fn load(dir: &Path) -> Result<RefData> {
     let Conventions {
         version,
         convention,
+        fx,
     } = source.decode()?;
     source.version(version)?;
     source.rows("convention", convention, |row: ConventionInput| {
         add_unique(&mut parts.conventions, row.id, row.kind)
     })?;
+    source.rows("fx", fx, |row: FxConvention| {
+        add_unique(&mut parts.fx_conventions, row.pair, row)
+    })?;
 
     let securities = dir.join("securities.toml");
     if securities
@@ -158,6 +172,7 @@ pub(super) fn load(dir: &Path) -> Result<RefData> {
                     RefDataEntry::Fallback(id) => ("indices.toml", format!("fallback[{id}]")),
                     RefDataEntry::Security(id) => ("securities.toml", format!("security[{id}]")),
                     RefDataEntry::Commodity(id) => ("", format!("commodity[{id}]")),
+                    RefDataEntry::FxConvention(pair) => ("conventions.toml", format!("fx[{pair}]")),
                 };
                 format!("{}: {field}: {error}", dir.join(file).display())
             })
diff --git a/crates/rhodium/src/workspace/trade_files.rs b/crates/rhodium/src/workspace/trade_files.rs
index 16bb8b67edf5fbab6af77f882c775a61b1e5d31b..3ffc1ce475ff45cea4a7da700e2d4c1270fb9835 100644
--- a/crates/rhodium/src/workspace/trade_files.rs
+++ b/crates/rhodium/src/workspace/trade_files.rs
@@ -12,10 +12,10 @@ use rhodium_core::{
 };
 use rhodium_pricing::resolve::CollateralKey;
 use rhodium_products::{
-    ir::Product,
     standard::PayReceive,
     trade::{
-        Amendment, ExerciseRecord, FallbackTerms, Fee, ProductSpec, Trade, TradeProductResolver,
+        Amendment, ExerciseRecord, FallbackTerms, Fee, ProductSpec, ResolvedTrade, Trade,
+        TradeProductResolver,
     },
 };
 use serde::Deserialize;
@@ -27,10 +27,10 @@ use super::{
 };
 
 // 読み込み時に IR へ変換し、マージンセットから担保条件を引いた取引。
+// 取引の条件は、方式の表の照合と解析プライサーの入力（標準商品、約定日、手数料）に使う。
 #[derive(Debug)]
 pub(super) struct LoadedTrade {
-    pub(super) book: BookId,
-    pub(super) product: Product,
+    pub(super) resolved: ResolvedTrade,
     pub(super) collateral: CollateralKey,
 }
 
@@ -263,14 +263,13 @@ fn resolve_trade(
     let trade = row
         .resolve(refdata)
         .map_err(|error| error.context(format!("trade {id}")))?;
-    let product = resolver
-        .resolve(&trade)
+    let resolved = resolver
+        .resolve(trade)
         .map_err(|error| error.context(format!("product of {id}")))?;
     Ok((
-        trade.id,
+        resolved.trade().id.clone(),
         LoadedTrade {
-            book: trade.book,
-            product,
+            resolved,
             collateral: collateral.clone(),
         },
     ))
diff --git a/crates/rhodium/tests/fx_conventions_io.rs b/crates/rhodium/tests/fx_conventions_io.rs
new file mode 100644
index 0000000000000000000000000000000000000000..c13bec6e7d78b8ed5cb44448461c91d66db28d37
--- /dev/null
+++ b/crates/rhodium/tests/fx_conventions_io.rs
@@ -0,0 +1,129 @@
+//! 通貨ペアの受渡慣行を参照データの `refdata/conventions.toml` の `[[fx]]` から読む。
+//!
+//! スポットラグ、同日受渡の締めと pip の桁は参照データに一つだけ置き、市場の定義の `[[fx]]` は
+//! quote と担保だけを書く。古い形式で市場の定義に `settlement` を書いた入力は未知のフィールドとして拒否する。
+
+#[allow(
+    dead_code,
+    reason = "共有のデスクの生成器のうち、ここで使わない定数がある"
+)]
+#[path = "support/desk.rs"]
+mod desk;
+
+use std::{fs, io::Write as _, path::Path};
+
+use rhodium::{prelude::ExecutionContext, session::Session, workspace::Workspace};
+use rhodium_core::{Error, Result, fx::CcyPair};
+
+#[expect(clippy::needless_pass_by_value, reason = "map_err に渡す変換関数")]
+fn io_error(error: std::io::Error) -> Error {
+    Error::Internal(error.to_string())
+}
+
+fn append(path: &Path, suffix: &str) -> Result<()> {
+    fs::OpenOptions::new()
+        .append(true)
+        .open(path)
+        .and_then(|mut file| file.write_all(suffix.as_bytes()))
+        .map_err(io_error)
+}
+
+/// 市場の定義の USDJPY のスポットだけの `[[fx]]`。受渡の慣行は書かない。
+const MARKET_FX: &str = r#"
+[[fx]]
+pair = "USDJPY"
+collaterals = [{ id = "market", base = "USD-SOFR", quote = "JPY-TONA" }]
+direct = { market_collateral = "market", quotes = { spot = "USDJPY.SPOT" } }
+"#;
+
+/// デスクに USDJPY のスポットを足し、参照データのスポットラグを `spot_lag` にした X(0)。
+fn x0(spot_lag: &str) -> Result<f64> {
+    let dir = tempfile::tempdir().map_err(io_error)?;
+    let root = dir.path();
+    desk::write(root, 0)?;
+    append(&root.join("market/definition.toml"), MARKET_FX)?;
+    append(&root.join("quotes/snapshot.csv"), "USDJPY.SPOT,150\n")?;
+    let conventions = root.join("refdata/conventions.toml");
+    let text = fs::read_to_string(&conventions).map_err(io_error)?;
+    fs::write(
+        &conventions,
+        text.replace(
+            "spot_lag = \"two_business_days\"",
+            &format!("spot_lag = \"{spot_lag}\""),
+        ),
+    )
+    .map_err(io_error)?;
+    let workspace = Workspace::load(root)?;
+    let as_of = workspace.resolve_as_of(desk::AS_OF)?;
+    let session = Session::new(workspace, as_of, ExecutionContext::sequential())?;
+    let pair: CcyPair = "USDJPY".parse()?;
+    session
+        .calibrated_market()
+        .market()
+        .with_view(|view| Ok(view.fx_today(pair)?.rate()))
+}
+
+#[test]
+fn workspace_reads_fx_conventions_from_refdata() -> Result<()> {
+    // デスクの参照データは USDJPY を T+2、東京 17 時の締め、pip 0.01 円とする。
+    // 市場の定義はスポットの quote と担保だけを書き、X(0) = X(S)·P_JPY(t_S)/P_USD(t_S) は
+    // スポットの受渡から評価時点まで両通貨の OIS（JPY 0.45%、USD 3.5% 付近）で割り戻した値になる。
+    // P_JPY/P_USD は 1 日あたり約 (3.5% − 0.45%)/365 ≈ 8.4e-5 だけ 1 を上回るので、X(0) はスポットより
+    // わずかに大きく、参照データのスポットラグだけを T+1 に変えると割り戻す日数が減って小さくなる。
+    let t_plus_two = x0("two_business_days")?;
+    let t_plus_one = x0("one_business_day")?;
+    for rate in [t_plus_two, t_plus_one] {
+        assert!(rate > 150.0 && rate < 150.1, "{rate}");
+    }
+    assert!(t_plus_two > t_plus_one, "{t_plus_two} {t_plus_one}");
+    Ok(())
+}
+
+#[test]
+fn workspace_rejects_settlement_in_market_fx() -> Result<()> {
+    // 受渡の慣行の正本は参照データだけにする。市場の定義に書いた `settlement` は、
+    // 行の位置付きの未知のフィールドとして拒否する。
+    let dir = tempfile::tempdir().map_err(io_error)?;
+    let root = dir.path();
+    desk::write(root, 0)?;
+    append(
+        &root.join("market/definition.toml"),
+        &MARKET_FX.replace(
+            "quotes = { spot = \"USDJPY.SPOT\" } }",
+            "quotes = { spot = \"USDJPY.SPOT\" }, settlement = { spot_lag = \"two_business_days\", settlement = { time = \"17:00:00\", zone = \"Asia/Tokyo\" }, pip_exponent = 2 } }",
+        ),
+    )?;
+    let Err(Error::InvalidData(message)) = Workspace::load(root) else {
+        return Err(Error::Internal(
+            "settlement in the market FX was accepted".into(),
+        ));
+    };
+    assert!(message.contains("definition.toml"), "{message}");
+    assert!(message.contains("fx[0]"), "{message}");
+    assert!(message.contains("unknown field `settlement`"), "{message}");
+
+    // 参照データの慣行の誤り（未登録の地域）は conventions.toml の fx[USDJPY] として報告する。
+    let dir = tempfile::tempdir().map_err(io_error)?;
+    let root = dir.path();
+    desk::write(root, 0)?;
+    let conventions = root.join("refdata/conventions.toml");
+    let text = fs::read_to_string(&conventions).map_err(io_error)?;
+    fs::write(
+        &conventions,
+        text.replace(
+            "settlement = { time = \"17:00:00\", zone = \"Asia/Tokyo\" }",
+            "settlement = { time = \"17:00:00\", zone = \"Europe/London\" }",
+        ),
+    )
+    .map_err(io_error)?;
+    let Err(Error::InvalidData(message)) = Workspace::load(root) else {
+        return Err(Error::Internal(
+            "unknown FX settlement zone was accepted".into(),
+        ));
+    };
+    assert!(
+        message.contains("conventions.toml: fx[USDJPY]"),
+        "{message}"
+    );
+    Ok(())
+}
diff --git a/crates/rhodium/tests/m1_rates.rs b/crates/rhodium/tests/m1_rates.rs
new file mode 100644
index 0000000000000000000000000000000000000000..fbe9551ab9dd7c933d8bd9034a274c1db9683b52
--- /dev/null
+++ b/crates/rhodium/tests/m1_rates.rs
@@ -0,0 +1,418 @@
+//! M1 の最初の検証データセット（`examples/m1-rates`）を、ファサードの `Session` で校正・評価する。
+//!
+//! 一つの市場の定義で JPY の TONA（会合日の階段と年末の跳び）、TIBOR 6M・3M のベーシス、
+//! USD の SOFR、USDJPY の外貨担保曲線（FX ポイントと `MtM` XCCY）を校正し、M1 のリニア商品と
+//! CSA の 4 つの場合の取引の PV、報告通貨の合計と quote デルタを求める。
+
+use std::{
+    collections::BTreeSet,
+    fs,
+    num::NonZeroUsize,
+    path::{Path, PathBuf},
+    ptr,
+};
+
+use rhodium::{
+    prelude::{ExecutionContext, ThreadExecutor},
+    session::{GroupOutcome, PortfolioValuation, Session, Valuation},
+    workspace::Workspace,
+};
+use rhodium_core::{Error, Result, Warning, currency::Currency, id::TradeId};
+use rhodium_market::{
+    quote::{QuoteId, QuoteKind, QuoteOverride},
+    view::Slots,
+};
+use rhodium_pricing::config::{AnalyticPricer, Method};
+use rhodium_risk::vector::RiskVector;
+use rhodium_testkit::{assert_close, check_aad_vs_bump};
+
+const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
+const REL: f64 = 1e-12;
+const JPY_ABS: f64 = 1e-6;
+const TONA_10Y: &str = "JPY.TONA.OIS.10Y";
+
+#[expect(
+    clippy::disallowed_macros,
+    reason = "リポジトリのサンプルの場所をテストのソースの位置から決める"
+)]
+fn sample() -> PathBuf {
+    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m1-rates")
+}
+
+#[expect(clippy::needless_pass_by_value, reason = "map_err に渡す変換関数")]
+fn io_error(error: std::io::Error) -> Error {
+    Error::Internal(error.to_string())
+}
+
+fn copy_dir(from: &Path, to: &Path) -> Result<()> {
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
+fn session(dir: &Path, execution: ExecutionContext<'_>) -> Result<Session> {
+    let workspace = Workspace::load(dir)?;
+    let as_of = workspace.resolve_date("2026-11-24".parse()?)?;
+    Session::new(workspace, as_of, execution)
+}
+
+fn valuation<'a>(report: &'a PortfolioValuation, id: &str) -> Result<&'a Valuation> {
+    report
+        .results
+        .iter()
+        .find(|(trade, _)| trade.as_str() == id)
+        .ok_or_else(|| Error::Internal(format!("missing result for {id}")))?
+        .1
+        .as_ref()
+        .map_err(Clone::clone)
+}
+
+fn has(value: &Valuation, warning: &Warning) -> bool {
+    value.warnings.contains(warning)
+}
+
+#[test]
+fn m1_rates_workspace_calibrates_and_values_every_trade() -> Result<()> {
+    // すべての校正グループが残差の許容誤差（1e-14）で解け、すべての取引が評価でき、合計から除く取引は無い。
+    let session = session(&sample(), SEQUENTIAL)?;
+    let market = session.calibrated_market();
+    assert_eq!(market.groups().len(), 9);
+    for group in market.groups() {
+        let GroupOutcome::Solved(solution) = group.outcome() else {
+            return Err(Error::Calibration(format!("{group:?}")));
+        };
+        assert!(solution.residual_norm() <= 1e-14, "{:?}", group.blocks());
+    }
+    let report = session.value_portfolio(SEQUENTIAL)?;
+    assert_eq!(report.excluded, [] as [TradeId; 0]);
+    assert_eq!(report.results.len(), 18);
+    let total = report.total.clone()?;
+    assert_eq!(total.ccy(), Currency::JPY);
+    assert!(total.amount().is_finite());
+
+    // 先物は方式の表の既定で解析プライサー、他はリニア評価。
+    for (id, value) in &report.results {
+        let value = value.as_ref().map_err(Clone::clone)?;
+        let expected = if id.as_str() == "JPY-TONA-FUT" {
+            Method::Analytic(AnalyticPricer::RateFuture)
+        } else {
+            Method::Linear
+        };
+        assert_eq!(value.method, expected, "{id}");
+    }
+    // 予想の測度: TONA 担保の TIBOR の IRS は担保の違いの近似を使わない。25bp の CSA の IRS と、
+    // JPY 担保で割り引く SOFR の OIS と XCCY は使う。USD 担保と無担保の SOFR の OIS は使わない。
+    let cross = Warning::CrossCollateralProjection;
+    for (id, expected) in [
+        ("JPY-IRS-5Y", false),
+        ("JPY-BASIS-3M-5Y", false),
+        ("JPY-IRS-5Y-CASH-25", true),
+        ("USD-OIS-JPY-CASH", true),
+        ("USD-OIS-CTD", true),
+        ("USD-OIS-USD-CASH", false),
+        ("USD-OIS-UNCOLLATERALIZED", false),
+        ("USDJPY-MTM-XCCY-2Y", true),
+        ("USDJPY-FWD-1Y", false),
+    ] {
+        assert_eq!(has(valuation(&report, id)?, &cross), expected, "{id}");
+    }
+    assert!(has(
+        valuation(&report, "JPY-FRA-6X12")?,
+        &Warning::FraDiscounting
+    ));
+    for id in ["JPY-OIS-LOOKBACK", "JPY-OIS-LOCKOUT"] {
+        assert!(
+            has(valuation(&report, id)?, &Warning::DailyForwardProduct),
+            "{id}"
+        );
+    }
+    assert!(!has(
+        valuation(&report, "JPY-OIS-SHIFT")?,
+        &Warning::DailyForwardProduct
+    ));
+    assert!(has(
+        valuation(&report, "USDJPY-MTM-XCCY-2Y")?,
+        &Warning::MtmFxForwardProduct
+    ));
+
+    // CSA の 4 つの場合: USD 担保と無担保（USD-SOFR で割り引く）は同じ PV、JPY 担保は外貨担保曲線で、
+    // CTD は毎年の平均フォワードの大きい方で割り引くので、PV は JPY 担保以下。
+    let usd = |id| valuation(&report, id).map(|value| value.pv.get(Currency::USD));
+    assert_close!(
+        usd("USD-OIS-USD-CASH")?,
+        usd("USD-OIS-UNCOLLATERALIZED")?,
+        REL,
+        1e-8
+    );
+    assert!((usd("USD-OIS-JPY-CASH")? - usd("USD-OIS-USD-CASH")?).abs() > 1.0);
+    assert!(usd("USD-OIS-CTD")? <= usd("USD-OIS-JPY-CASH")? + 1e-8);
+    Ok(())
+}
+
+/// サンプルの写しの取引ファイルを書き換える。
+fn replace_in(path: &Path, from: &str, to: &str) -> Result<()> {
+    let text = fs::read_to_string(path).map_err(io_error)?;
+    assert!(text.contains(from), "{from}");
+    fs::write(path, text.replace(from, to)).map_err(io_error)
+}
+
+/// 一点の割引因子を読む。
+fn discount(session: &Session, curve: &str, at: &str) -> Result<f64> {
+    let market = session.calibrated_market().market();
+    let structure = market.structure();
+    let curve = structure.resolve_discount(&curve.parse()?)?;
+    let point = curve.prepare_point(at.parse()?)?;
+    let slots = Slots::from_resolved(std::sync::Arc::clone(structure), [(curve, vec![point])], [])?;
+    market.with_view(|view| {
+        Ok(view
+            .df(view.bind_slots(&slots)?.discounts()[0].points()[0])
+            .value())
+    })
+}
+
+#[test]
+fn amended_trade_matches_the_scaled_original_and_its_fee() -> Result<()> {
+    // 一部解約は効力の瞬間（評価時点）以後の受渡をすべて 0.4 倍にし、手数料 200 万円を別に受け取る。
+    // 同じ IRS を変更せずに書いた写しの PV から、0.4 × PV と手数料の PV の和を手計算で求める。
+    let dir = tempfile::tempdir().map_err(io_error)?;
+    copy_dir(&sample(), dir.path())?;
+    let amended = valuation(
+        &session(dir.path(), SEQUENTIAL)?.value_portfolio(SEQUENTIAL)?,
+        "JPY-IRS-AMENDED",
+    )?
+    .pv
+    .get(Currency::JPY);
+    let rates = dir.path().join("trades/rates.toml");
+    replace_in(
+        &rates,
+        "fees = [{ amount = \"2_000_000 JPY\", side = \"receive\", at = 2026-11-26T08:00:00Z }]\n",
+        "",
+    )?;
+    replace_in(
+        &rates,
+        "type = \"amended\"\nfrom = 2026-11-24T08:00:00Z\nafter = { scale = 0.4 }\n\n[trade.product.before]\n",
+        "",
+    )?;
+    let original_session = session(dir.path(), SEQUENTIAL)?;
+    let original = valuation(
+        &original_session.value_portfolio(SEQUENTIAL)?,
+        "JPY-IRS-AMENDED",
+    )?
+    .pv
+    .get(Currency::JPY);
+    let fee = 2_000_000.0 * discount(&original_session, "JPY-TONA", "2026-11-26T08:00:00Z")?;
+    assert_close!(
+        amended,
+        0.4 * original + fee,
+        REL,
+        JPY_ABS,
+        original.abs() + fee
+    );
+    Ok(())
+}
+
+fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
+    vector
+        .entries()
+        .iter()
+        .find_map(|(key, value)| (key == id).then_some(*value))
+        .unwrap_or(0.0)
+}
+
+fn reporting(session: &Session) -> Result<RiskVector<QuoteId>> {
+    let report = session.quote_delta_reporting(SEQUENTIAL)?;
+    assert!(report.failures().is_empty());
+    let [vector] = report.by_ccy() else {
+        return Err(Error::Internal("expected one reporting currency".into()));
+    };
+    Ok(vector.clone())
+}
+
+#[test]
+fn m1_rates_quote_update_re_solves_downstream_groups_only() -> Result<()> {
+    // M1-7: TONA の 10Y の quote を 1bp 上げると、TONA の長期と、その下流の TIBOR 6M・3M と
+    // 外貨担保曲線だけを解き直し、他のグループ（年末の跳び、階段、SOFR、FX と LGM のマーク）と
+    // 取引の準備を共有する。結果は更新後の quote から新規に作った Session と許容誤差で一致する。
+    let dir = tempfile::tempdir().map_err(io_error)?;
+    copy_dir(&sample(), dir.path())?;
+    let base = session(dir.path(), SEQUENTIAL)?;
+    let id: QuoteId = TONA_10Y.parse()?;
+    let overrides = base.quote_overrides([QuoteOverride::Add {
+        pattern: id.as_str().parse()?,
+        shift: 1e-4,
+    }])?;
+    let updated = base.with_quotes(&overrides, SEQUENTIAL)?;
+    let resolved: BTreeSet<String> = base
+        .calibrated_market()
+        .groups()
+        .iter()
+        .zip(updated.calibrated_market().groups())
+        .filter(|(old, new)| !ptr::eq(old.outcome(), new.outcome()))
+        .map(|(old, _)| format!("{} {}", old.blocks()[0].curve(), old.blocks()[0].part()))
+        .collect();
+    assert_eq!(
+        resolved,
+        [
+            "JPY-TIBOR-3M basis",
+            "JPY-TIBOR-6M basis",
+            "JPY-TONA long",
+            "USD-under-JPY basis"
+        ]
+        .map(String::from)
+        .into()
+    );
+    for ((before, _), (after, _)) in base
+        .preparation_warnings()
+        .zip(updated.preparation_warnings())
+    {
+        assert!(ptr::eq(before, after));
+    }
+
+    let snapshot = dir.path().join("quotes/snapshot.csv");
+    let text = fs::read_to_string(&snapshot).map_err(io_error)?;
+    fs::write(
+        &snapshot,
+        text.replace(&format!("{TONA_10Y},1.32%"), &format!("{TONA_10Y},1.33%")),
+    )
+    .map_err(io_error)?;
+    let fresh = session(dir.path(), SEQUENTIAL)?;
+    // 更新は前の θ から、新規は固定の初期値から解くので、θ はソルバーの誤差を伝えた幅で一致する。
+    for (lhs, rhs) in updated
+        .calibrated_market()
+        .market()
+        .theta()
+        .iter()
+        .zip(fresh.calibrated_market().market().theta())
+    {
+        assert_close!(*lhs, *rhs, 1e-10, 1e-12);
+    }
+    let [lhs, rhs] = [&updated, &fresh].map(|session| session.value_portfolio(SEQUENTIAL));
+    let (lhs, rhs) = (lhs?, rhs?);
+    for ((id, left), (_, right)) in lhs.results.iter().zip(&rhs.results) {
+        let [left, right] = [left, right].map(|value| value.as_ref().map_err(Clone::clone));
+        assert_close!(
+            left?.reporting.amount(),
+            right?.reporting.amount(),
+            1e-10,
+            JPY_ABS,
+            1e9
+        );
+        let _ = id;
+    }
+    let [lhs, rhs] = [reporting(&updated)?, reporting(&fresh)?];
+    let scale = rhs
+        .entries()
+        .iter()
+        .map(|(_, value)| value.abs())
+        .sum::<f64>();
+    for id in lhs.entries().iter().chain(rhs.entries()).map(|(id, _)| id) {
+        assert_close!(
+            sensitivity(&lhs, id),
+            sensitivity(&rhs, id),
+            1e-8,
+            JPY_ABS,
+            scale
+        );
+    }
+    Ok(())
+}
+
+#[test]
+fn m1_rates_values_do_not_depend_on_thread_count() -> Result<()> {
+    // 取引の準備、PV、報告通貨のデルタを 3 つのワーカーに分けても、逐次と許容誤差で一致する。
+    let sequential = session(&sample(), SEQUENTIAL)?;
+    let executor = ThreadExecutor::new(NonZeroUsize::MIN.saturating_add(2));
+    let parallel_execution = ExecutionContext::new(&executor);
+    let parallel = session(&sample(), parallel_execution)?;
+    let [lhs, rhs] = [
+        sequential.value_portfolio(SEQUENTIAL)?,
+        parallel.value_portfolio(parallel_execution)?,
+    ];
+    for ((id, left), (right_id, right)) in lhs.results.iter().zip(&rhs.results) {
+        assert_eq!(id, right_id);
+        let [left, right] = [left, right].map(|value| value.as_ref().map_err(Clone::clone));
+        assert_close!(
+            left?.reporting.amount(),
+            right?.reporting.amount(),
+            REL,
+            JPY_ABS
+        );
+    }
+    let lhs = reporting(&sequential)?;
+    let report = parallel.quote_delta_reporting(parallel_execution)?;
+    let [rhs] = report.by_ccy() else {
+        return Err(Error::Internal("expected one reporting currency".into()));
+    };
+    for id in lhs.entries().iter().chain(rhs.entries()).map(|(id, _)| id) {
+        assert_close!(sensitivity(&lhs, id), sensitivity(rhs, id), 1e-10, JPY_ABS);
+    }
+    Ok(())
+}
+
+#[test]
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#203 で各 variant を網羅する"
+)]
+fn m1_rates_portfolio_quote_delta_matches_bump_and_recalibrate() -> Result<()> {
+    // 最初の検証データセットのポートフォリオ全体の報告通貨（JPY）の PV について、quote デルタ
+    // （AAD、FX の換算の連鎖律、多段の校正の随伴）を、quote ごとに影響範囲を再校正した中心差分と照合する。
+    // 再校正は Session::with_quotes（残差の許容誤差 1e-14）なので、その誤差を δ_V に伝える。
+    let base = session(&sample(), SEQUENTIAL)?;
+    let delta = reporting(&base)?;
+    let total = |session: &Session| -> Result<f64> {
+        Ok(session.value_portfolio(SEQUENTIAL)?.total?.amount())
+    };
+    let market = base.calibrated_market();
+    let unit = |id: &QuoteId| -> Result<f64> {
+        let quote = market
+            .quotes()
+            .get(id)
+            .ok_or_else(|| Error::Internal(format!("missing quote {id}")))?;
+        Ok(
+            match market
+                .quote_kind(id)
+                .ok_or_else(|| Error::Internal(format!("unknown kind of {id}")))?
+            {
+                QuoteKind::FxSpot => 0.01 * quote,
+                QuoteKind::FxForwardPoints | QuoteKind::Mark => 1.0,
+                _ => 1e-4,
+            },
+        )
+    };
+    let scale = delta
+        .entries()
+        .iter()
+        .map(|(id, value)| Ok(value.abs() / unit(id)?))
+        .sum::<Result<f64>>()?;
+    for (id, _) in market.quotes().iter() {
+        let (report_unit, width) = match market.quote_kind(id) {
+            // マークは報告単位 1 だが、1 だけ動かすと凸性の非線形が大きいので 1bp の幅で照合する。
+            Some(QuoteKind::Mark) => (1.0, 1e-4),
+            _ => (unit(id)?, unit(id)?),
+        };
+        check_aad_vs_bump(
+            |h| {
+                let overrides = base.quote_overrides([QuoteOverride::Add {
+                    pattern: id.as_str().parse()?,
+                    shift: h,
+                }])?;
+                total(&base.with_quotes(&overrides, SEQUENTIAL)?)
+            },
+            sensitivity(&delta, id) * width / report_unit,
+            width,
+            1e-14 * scale,
+            JPY_ABS,
+            id.as_str(),
+        )?;
+    }
+    Ok(())
+}
diff --git a/crates/rhodium/tests/market_definition_io.rs b/crates/rhodium/tests/market_definition_io.rs
new file mode 100644
index 0000000000000000000000000000000000000000..1daeab48410191cb3b4db3c672ed38b23a8f0662
--- /dev/null
+++ b/crates/rhodium/tests/market_definition_io.rs
@@ -0,0 +1,435 @@
+//! 市場の定義 `market/definition.toml` の曲線の式、部品を指す当てはめ、直接マーク、最小二乗の指定。
+//!
+//! ファイルの曲線は `CurveExpr` の serde の表現そのもので、瞬間は UTC の日時で明示する。
+//! 会合表や年末の規則を参照データから解決する書き方（M6）はまだ読まない。
+
+#[allow(
+    dead_code,
+    reason = "共有のデスクの生成器のうち、ここで使わない定数がある"
+)]
+#[path = "support/desk.rs"]
+mod desk;
+
+use std::{fmt::Write as _, fs, io::Write as _, path::Path};
+
+use rhodium::{prelude::ExecutionContext, session::Session, workspace::Workspace};
+use rhodium_core::{Error, Result, time::Instant};
+use rhodium_market::{
+    curve::{CurveExpr, CurveId},
+    market::Market,
+    structure::BlockId,
+    view::Slots,
+};
+use rhodium_pricing::calibration::{GroupOutcome, HessianMode};
+use rhodium_testkit::assert_close;
+
+const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
+
+fn fmt_error(error: std::fmt::Error) -> Error {
+    Error::Internal(error.to_string())
+}
+
+#[expect(clippy::needless_pass_by_value, reason = "map_err に渡す変換関数")]
+fn io_error(error: std::io::Error) -> Error {
+    Error::Internal(error.to_string())
+}
+
+fn append(path: &Path, suffix: &str) -> Result<()> {
+    fs::OpenOptions::new()
+        .append(true)
+        .open(path)
+        .and_then(|mut file| file.write_all(suffix.as_bytes()))
+        .map_err(io_error)
+}
+
+/// デスクの JPY の OIS の quote ID（1M、3M、6M、1Y、2Y … 30Y）とピラー（スポット + テナー + 7 日）。
+const SHORT: [&str; 4] = ["1M", "3M", "6M", "1Y"];
+const LONG: [(&str, &str); 11] = [
+    ("2Y", "2028-12-03"),
+    ("3Y", "2029-12-03"),
+    ("4Y", "2030-12-03"),
+    ("5Y", "2031-12-03"),
+    ("7Y", "2033-12-03"),
+    ("10Y", "2036-12-03"),
+    ("12Y", "2038-12-03"),
+    ("15Y", "2041-12-03"),
+    ("20Y", "2046-12-03"),
+    ("25Y", "2051-12-03"),
+    ("30Y", "2056-12-03"),
+];
+/// 1Y のピラー（接続点）。
+const STITCH: &str = "2027-12-03T00:00:00Z";
+/// 年末の跳びの発効（2027 年の最初の営業日）と、その受渡期間の日数係数 4/365。
+const TURN: &str = "2027-01-04T00:00:00Z";
+const TURN_FRACTION: f64 = 4.0 / 365.0;
+
+/// JPY-TONA を、会合日の階段と長期の節点の接続と年末の跳びの和で書いた市場の定義。
+/// `meetings` の会合で階段を区切り、`short_rows` に短期の当てはめの行の追加の指定を書く。
+fn tona_definition(meetings: &[&str], short_extra: &str, extra: &str) -> Result<String> {
+    let mut text = String::from(
+        r#"version = 1
+valuation_calendar = "JPTO"
+valuation_zone = "Asia/Tokyo"
+
+[curves.JPY-TONA]
+currency = "JPY"
+roles = ["discount", "projection"]
+
+[[curves.JPY-TONA.expr.sum]]
+"#,
+    );
+    writeln!(text, "stitch.at = {STITCH}").map_err(fmt_error)?;
+    writeln!(
+        text,
+        "stitch.short.steps = {{ part = \"meetings\", meetings = [{}] }}",
+        meetings
+            .iter()
+            .map(|date| format!("{date}T00:00:00Z"))
+            .collect::<Vec<_>>()
+            .join(", ")
+    )
+    .map_err(fmt_error)?;
+    writeln!(
+        text,
+        "stitch.long.nodes = {{ part = \"long\", pillars = [{}] }}",
+        LONG.iter()
+            .map(|(_, date)| format!("{date}T00:00:00Z"))
+            .collect::<Vec<_>>()
+            .join(", ")
+    )
+    .map_err(fmt_error)?;
+    writeln!(
+        text,
+        "\n[[curves.JPY-TONA.expr.sum]]\njumps = {{ part = \"turns\", dates = [{TURN}] }}"
+    )
+    .map_err(fmt_error)?;
+    text.push_str(
+        r#"
+[curves.USD-SOFR]
+currency = "USD"
+roles = ["discount", "projection"]
+pillars = [2026-12-26T00:00:00Z, 2027-03-05T00:00:00Z, 2027-06-02T00:00:00Z, 2027-12-03T00:00:00Z, 2028-12-03T00:00:00Z]
+"#,
+    );
+    let ois = |tenor: &str, part: &str, extra: &str| {
+        format!(
+            "\n[[instruments]]\nquote = \"JPY.TONA.OIS.{tenor}\"\nkind = \"rate\"\ntarget = \"par_rate\"\n\
+             template = {{ ois = {{ convention = \"JPY-TONA-OIS\", tenor = \"{tenor}\" }} }}\n\
+             collateral = {{ cash = {{ ccy = \"JPY\", index = \"JPY-TONA\", spread = \"0bp\" }} }}\n\
+             determines = {{ curve = \"JPY-TONA\", part = \"{part}\" }}\n{extra}"
+        )
+    };
+    for tenor in SHORT {
+        text.push_str(&ois(tenor, "meetings", short_extra));
+    }
+    for (tenor, _) in LONG {
+        text.push_str(&ois(tenor, "long", ""));
+    }
+    for (element, tenor) in ["1M", "3M", "6M", "1Y", "2Y"].iter().enumerate() {
+        writeln!(text,
+            "\n[[instruments]]\nquote = \"USD.SOFR.OIS.{tenor}\"\nkind = \"rate\"\ntarget = \"par_rate\"\n\
+             template = {{ ois = {{ convention = \"USD-SOFR-OIS\", tenor = \"{tenor}\" }} }}\n\
+             collateral = {{ cash = {{ ccy = \"USD\", index = \"USD-SOFR\", spread = \"0bp\" }} }}\n\
+             determines = {{ curve = \"USD-SOFR\", element = {element} }}"
+        ).map_err(fmt_error)?;
+    }
+    write!(
+        text,
+        r#"
+[[instruments]]
+quote = "JPY.TONA.TURN"
+kind = "spread"
+target = "direct_mark"
+transform = {{ turn_spread_to_jump = {{ value = {TURN_FRACTION}, basis = "ACT365F" }} }}
+determines = {{ curve = "JPY-TONA", part = "turns", element = 0 }}
+"#
+    )
+    .map_err(fmt_error)?;
+    text.push_str(extra);
+    Ok(text)
+}
+
+/// デスクのワークスペースの市場の定義を差し替え、年末のスプレッドの quote を足す。
+fn workspace(definition: &str) -> Result<tempfile::TempDir> {
+    let dir = tempfile::tempdir().map_err(io_error)?;
+    desk::write(dir.path(), 0)?;
+    fs::write(dir.path().join("market/definition.toml"), definition).map_err(io_error)?;
+    append(
+        &dir.path().join("quotes/snapshot.csv"),
+        "JPY.TONA.TURN,25bp\n",
+    )?;
+    Ok(dir)
+}
+
+fn session(dir: &Path) -> Result<Session> {
+    let workspace = Workspace::load(dir)?;
+    let as_of = workspace.resolve_as_of(desk::AS_OF)?;
+    Session::new(workspace, as_of, SEQUENTIAL)
+}
+
+fn instants(dates: &[&str]) -> Result<Vec<Instant>> {
+    dates
+        .iter()
+        .map(|date| format!("{date}T00:00:00Z").parse())
+        .collect()
+}
+
+#[test]
+fn market_file_reads_tona_steps_and_turn() -> Result<()> {
+    // 3 回の会合で区切った 4 区間の階段を 1M、3M、6M、1Y の OIS が、長期の 11 の節点を 2Y〜30Y の OIS が
+    // 決め、年末のスプレッド 25bp が跳び 25bp × 4/365 を直接決める。読んだ式は Rust で組んだ式と一致し、
+    // 正方系の各グループは quote を残差の許容誤差で再現する。
+    let meetings = ["2027-01-11", "2027-03-15", "2027-06-14"];
+    let dir = workspace(&tona_definition(&meetings, "", "")?)?;
+    let session = session(dir.path())?;
+    let market = session.calibrated_market();
+    let id: CurveId = "JPY-TONA".parse()?;
+    let structure = market.market().structure();
+    let expected = CurveExpr::Sum(vec![
+        CurveExpr::Stitch {
+            short: Box::new(CurveExpr::Steps {
+                part: "meetings".into(),
+                meetings: instants(&meetings)?,
+            }),
+            long: Box::new(CurveExpr::Nodes {
+                part: "long".into(),
+                pillars: instants(&LONG.map(|(_, date)| date))?,
+                extrapolation: rhodium_market::curve::Extrapolation::default(),
+            }),
+            at: STITCH.parse()?,
+        },
+        CurveExpr::Jumps {
+            part: "turns".into(),
+            dates: vec![TURN.parse()?],
+        },
+    ]);
+    assert_eq!(
+        structure
+            .curve(&id)
+            .ok_or_else(|| Error::Internal("missing TONA".into()))?
+            .expr(),
+        &expected
+    );
+    for group in market.groups() {
+        let GroupOutcome::Solved(solution) = group.outcome() else {
+            return Err(Error::Calibration(format!("{group:?}")));
+        };
+        assert!(solution.residual_norm() <= 1e-14);
+        assert_eq!(solution.hessian(), None);
+    }
+    let turn = structure
+        .layout()
+        .range(&BlockId::new(id, "turns")?)
+        .ok_or_else(|| Error::Internal("missing turn block".into()))?;
+    assert_close!(
+        market.market().theta()[turn.start],
+        0.0025 * TURN_FRACTION,
+        1e-12,
+        1e-14
+    );
+    Ok(())
+}
+
+/// 一点の割引因子を読む。
+fn df(market: &Market<f64>, curve: &str, at: Instant) -> Result<f64> {
+    let structure = market.structure();
+    let curve = structure.resolve_discount(&curve.parse()?)?;
+    let point = curve.prepare_point(at)?;
+    let slots = Slots::from_resolved(std::sync::Arc::clone(structure), [(curve, vec![point])], [])?;
+    market.with_view(|view| {
+        Ok(view
+            .df(view.bind_slots(&slots)?.discounts()[0].points()[0])
+            .value())
+    })
+}
+
+#[test]
+fn market_file_reads_spread_ctd_and_floor_curves() -> Result<()> {
+    // 担保条件ごとの割引曲線を式で書く。TONA + 25bp の付利の曲線、SOFR と SOFR + 10bp の CTD、
+    // 0% のフロア付きの TONA の付利の曲線。どれも θ を持たず、校正した曲線から手計算で求まる。
+    let extra = r#"
+[curves.JPY-TONA-25]
+currency = "JPY"
+roles = ["discount"]
+expr.sum = [{ ref = "JPY-TONA" }, { flat = "25bp" }]
+
+[curves.USD-SOFR-10]
+currency = "USD"
+roles = ["discount"]
+expr.sum = [{ ref = "USD-SOFR" }, { flat = "0.1%" }]
+
+[curves.USD-CTD]
+currency = "USD"
+roles = ["discount"]
+expr.max_forward = { parts = [{ ref = "USD-SOFR" }, { ref = "USD-SOFR-10" }], grid = [2027-06-02T00:00:00Z, 2027-12-03T00:00:00Z] }
+
+[curves.JPY-TONA-FLOOR]
+currency = "JPY"
+roles = ["discount"]
+
+[[curves.JPY-TONA-FLOOR.expr.sum]]
+ref = "JPY-TONA"
+
+[[curves.JPY-TONA-FLOOR.expr.sum]]
+max_forward = { parts = [{ flat = "0bp" }, { sum = [{ flat = "0bp" }, { scale = { factor = -1.0, body = { forward_ref = "JPY-TONA" } } }] }], grid = [2027-06-02T00:00:00Z, 2027-12-03T00:00:00Z] }
+"#;
+    let dir = workspace(&tona_definition(
+        &["2027-01-11", "2027-03-15", "2027-06-14"],
+        "",
+        extra,
+    )?)?;
+    let session = session(dir.path())?;
+    let market = session.calibrated_market().market();
+    let at: Instant = "2027-12-03T00:00:00Z".parse()?;
+    let t = rhodium_core::time::model_time(market.structure().as_of(), at).years();
+    let tona = df(market, "JPY-TONA", at)?;
+    let sofr = df(market, "USD-SOFR", at)?;
+    // 一定のスプレッドの曲線は DF に exp(−s·t) を掛ける。
+    assert_close!(
+        df(market, "JPY-TONA-25", at)?,
+        tona * (-0.0025 * t).exp(),
+        1e-12,
+        1e-14
+    );
+    assert_close!(
+        df(market, "USD-SOFR-10", at)?,
+        sofr * (-0.001 * t).exp(),
+        1e-12,
+        1e-14
+    );
+    // 各区間で SOFR + 10bp のフォワードが大きいので、CTD は格子の点で SOFR + 10bp と一致する。
+    assert_close!(
+        df(market, "USD-CTD", at)?,
+        sofr * (-0.001 * t).exp(),
+        1e-12,
+        1e-14
+    );
+    // TONA のフォワードは正なので、0% のフロアは効かず、補正 max(0, 0 − f) は 0。
+    assert_close!(df(market, "JPY-TONA-FLOOR", at)?, tona, 1e-12, 1e-14);
+    Ok(())
+}
+
+#[test]
+fn market_file_reads_least_squares_options() -> Result<()> {
+    // 4 回の会合（5 区間）を 4 本の OIS で当てるには未知数が足りないので、短期の行に許容幅を付け、
+    // 階段に罰則を掛けた最小二乗にする。`gauss_newton` を書くと記録を Gauss–Newton で作る。
+    let meetings = ["2027-01-11", "2027-03-15", "2027-06-14", "2027-09-13"];
+    let penalty = r#"
+[[step_penalties]]
+curve = "JPY-TONA"
+part = "meetings"
+weight = 1.0
+scale = "1bp"
+"#;
+    let tolerance = "tolerance = \"0.5bp\"\n";
+    let meetings_block = BlockId::new("JPY-TONA".parse()?, "meetings")?;
+    for (extra, mode) in [
+        (penalty.to_owned(), HessianMode::Exact),
+        (
+            format!("{penalty}\n[[gauss_newton]]\ncurve = \"JPY-TONA\"\npart = \"meetings\"\n"),
+            HessianMode::GaussNewton,
+        ),
+    ] {
+        let dir = workspace(&tona_definition(&meetings, tolerance, &extra)?)?;
+        let session = session(dir.path())?;
+        let group = session
+            .calibrated_market()
+            .groups()
+            .iter()
+            .find(|group| group.blocks() == [meetings_block.clone()])
+            .ok_or_else(|| Error::Internal("missing meetings group".into()))?;
+        let GroupOutcome::Solved(solution) = group.outcome() else {
+            return Err(Error::Calibration(format!("{group:?}")));
+        };
+        assert_eq!(solution.hessian(), Some(mode));
+    }
+    Ok(())
+}
+
+#[test]
+fn market_file_reports_invalid_curve_expressions_with_positions() -> Result<()> {
+    // 書き方の誤りは、ファイルと行（曲線の ID か instruments の番号）と理由を付けた入力のエラーにする。
+    let meetings = ["2027-01-11", "2027-03-15", "2027-06-14"];
+    let base = tona_definition(&meetings, "", "")?;
+    let cases: Vec<(String, &str)> = vec![
+        (
+            base.replacen(
+                "roles = [\"discount\", \"projection\"]\n\n[[curves.JPY-TONA.expr.sum]]",
+                "roles = [\"discount\", \"projection\"]\npillars = [2027-12-03T00:00:00Z]\n\n[[curves.JPY-TONA.expr.sum]]",
+                1,
+            ),
+            "pillars and expr are exclusive",
+        ),
+        (
+            base.replace("stitch.long.nodes", "stitch.long.nodez"),
+            "unknown variant",
+        ),
+        (
+            base.replace(
+                "jumps = { part = \"turns\"",
+                "fx_implied = { part = \"turns\"",
+            ),
+            "unknown variant",
+        ),
+        (
+            base.replace("2027-03-15T00:00:00Z", "2027-01-01T00:00:00Z"),
+            "curves[JPY-TONA]",
+        ),
+        (
+            format!(
+                "{base}\n[curves.JPY-FLAT]\ncurrency = \"JPY\"\nroles = [\"discount\"]\nexpr.flat = 0.0025\n"
+            ),
+            "curves[JPY-FLAT]",
+        ),
+        // 曲線の率は取引の率と同じ書式で読み、単位の前の空白を受け付けない。
+        (
+            format!(
+                "{base}\n[curves.JPY-FLAT]\ncurrency = \"JPY\"\nroles = [\"discount\"]\nexpr.flat = \"25 bp\"\n"
+            ),
+            "rate must be a scaled decimal",
+        ),
+        (
+            format!(
+                "{base}\n[curves.JPY-GRID]\ncurrency = \"JPY\"\nroles = [\"discount\"]\nexpr.max_forward = {{ parts = [{{ ref = \"JPY-TONA\" }}], grid = [2027-06-02] }}\n"
+            ),
+            "invalid TOML date or instant literal in grid",
+        ),
+        (
+            base.replace(
+                "target = \"direct_mark\"\n",
+                "target = \"direct_mark\"\ntemplate = { ois = { convention = \"JPY-TONA-OIS\", tenor = \"1M\" } }\n",
+            ),
+            "direct_mark takes neither",
+        ),
+        (
+            base.replacen(
+                "template = { ois = { convention = \"JPY-TONA-OIS\", tenor = \"1M\" } }\n",
+                "",
+                1,
+            ),
+            "require a template",
+        ),
+        (
+            base.replace(
+                &format!("{{ value = {TURN_FRACTION}, basis = \"ACT365F\" }}"),
+                "0.0109",
+            ),
+            "instruments[",
+        ),
+        (
+            base.replacen("part = \"long\" }", "part = \"longer\" }", 1),
+            "unknown curve",
+        ),
+    ];
+    for (definition, expected) in cases {
+        let dir = workspace(&definition)?;
+        let Err(Error::InvalidData(message)) = Workspace::load(dir.path()) else {
+            return Err(Error::Internal(format!(
+                "accepted a definition with {expected}"
+            )));
+        };
+        assert!(message.contains("definition.toml"), "{message}");
+        assert!(message.contains(expected), "{expected}: {message}");
+    }
+    Ok(())
+}
diff --git a/crates/rhodium/tests/parallel.rs b/crates/rhodium/tests/parallel.rs
index aad527e43e3323cf9ba4c4c9cfb741a9939713fa..8c455ad0dc940560cdb6329e33017518ff767450 100644
--- a/crates/rhodium/tests/parallel.rs
+++ b/crates/rhodium/tests/parallel.rs
@@ -277,12 +277,14 @@ fn session_quote_delta_uses_configured_reporting_currency() -> Result<()> {
     .map_err(io_error)?;
     let definition = dir.path().join("market/definition.toml");
     let mut market = fs::read_to_string(&definition).map_err(io_error)?;
-    market.push_str(r#"
+    market.push_str(
+        r#"
 [[fx]]
 pair = "USDJPY"
 collaterals = [{ id = "market", base = "USD-SOFR", quote = "JPY-TONA" }]
-direct = { market_collateral = "market", quotes = { spot = "USDJPY.SPOT" }, settlement = { spot_lag = "two_business_days", settlement = { time = "17:00:00", zone = "Asia/Tokyo" }, pip_exponent = 2 } }
-"#);
+direct = { market_collateral = "market", quotes = { spot = "USDJPY.SPOT" } }
+"#,
+    );
     fs::write(definition, market).map_err(io_error)?;
     let snapshot = dir.path().join("quotes/snapshot.csv");
     let mut quotes = fs::read_to_string(&snapshot).map_err(io_error)?;
diff --git a/crates/rhodium/tests/portfolio.rs b/crates/rhodium/tests/portfolio.rs
index 07ba2fa4909afe5f6489c97942d01ed6756bb2ae..b500c79935e5fe75207142e5ef4cce71533f8d8c 100644
--- a/crates/rhodium/tests/portfolio.rs
+++ b/crates/rhodium/tests/portfolio.rs
@@ -78,7 +78,7 @@ fn add_usdjpy(root: &Path) -> Result<()> {
 [[fx]]
 pair = "USDJPY"
 collaterals = [{ id = "market", base = "USD-SOFR", quote = "JPY-TONA" }]
-direct = { market_collateral = "market", quotes = { spot = "USDJPY.SPOT" }, settlement = { spot_lag = "two_business_days", settlement = { time = "17:00:00", zone = "Asia/Tokyo" }, pip_exponent = 2 } }
+direct = { market_collateral = "market", quotes = { spot = "USDJPY.SPOT" } }
 "#,
     )?;
     append(&root.join("quotes/snapshot.csv"), "USDJPY.SPOT,150\n")
@@ -194,6 +194,55 @@ fn portfolio_isolates_trade_failures_and_totals_successes() -> Result<()> {
     let report = restored.value_portfolio(SEQUENTIAL)?;
     assert_eq!(excluded(&report), ["T00004"]);
     assert_close!(report.total?.amount(), sum, REL, JPY_ABS, scale);
+    // 品質の「部分的な校正」の行: 補って再校正した市場の報告通貨のデルタも、全体を解き直した値と
+    // 一致する（更新は前の θ から、全体は固定の初期値から解くので、ソルバーの誤差の幅で比べる）。
+    assert_reporting_deltas_match(&restored, &complete)?;
+    Ok(())
+}
+
+/// 二つのセッションの報告通貨の quote デルタが、ソルバーの誤差の幅で一致することを確かめる。
+fn assert_reporting_deltas_match(restored: &Session, complete: &Session) -> Result<()> {
+    let restored_delta = restored.quote_delta_reporting(SEQUENTIAL)?;
+    let complete_delta = complete.quote_delta_reporting(SEQUENTIAL)?;
+    assert!(
+        restored_delta
+            .failures()
+            .keys()
+            .eq(complete_delta.failures().keys())
+    );
+    let [restored_vector] = restored_delta.by_ccy() else {
+        return Err(Error::Internal("expected one reporting currency".into()));
+    };
+    let [complete_vector] = complete_delta.by_ccy() else {
+        return Err(Error::Internal("expected one reporting currency".into()));
+    };
+    let ids: std::collections::BTreeSet<_> = restored_vector
+        .entries()
+        .iter()
+        .chain(complete_vector.entries())
+        .map(|(id, _)| id)
+        .collect();
+    let value = |vector: &rhodium_risk::vector::RiskVector<QuoteId>, id: &QuoteId| {
+        vector
+            .entries()
+            .iter()
+            .find_map(|(key, delta)| (key == id).then_some(*delta))
+            .unwrap_or(0.0)
+    };
+    let delta_scale = complete_vector
+        .entries()
+        .iter()
+        .map(|(_, delta)| delta.abs())
+        .sum::<f64>();
+    for id in ids {
+        assert_close!(
+            value(restored_vector, id),
+            value(complete_vector, id),
+            1e-10,
+            1e-6,
+            delta_scale
+        );
+    }
     Ok(())
 }
 
@@ -421,6 +470,26 @@ fn session_reprepares_for_intraday_instant_fixings_and_collateral() -> Result<()
     assert!((standard - alternate).abs() > 1_000.0);
     assert_jpy_pv(&corrected, &corrected_report, "CSA", standard)?;
     assert_jpy_pv(&collateral_changed, &collateral_report, "CSA", alternate)?;
+
+    // マージンセットは変えず、その CSA の付利スプレッドだけを 0bp から 25bp に変えても、
+    // 同じ取引の担保の鍵が変わり、準備し直して 25bp の割引曲線で評価する。
+    fs::write(&trades, &input).map_err(io_error)?;
+    let agreements = root.join("agreements.toml");
+    let terms = fs::read_to_string(&agreements).map_err(io_error)?;
+    let cash = r#"id = "JPY-CASH"
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }"#;
+    assert!(terms.contains(cash));
+    fs::write(
+        &agreements,
+        terms.replace(cash, &cash.replace("0bp", "25bp")),
+    )
+    .map_err(io_error)?;
+    let terms_changed = session(root, intraday_payment)?;
+    assert_new_preparation(&corrected, &terms_changed);
+    let terms_report = terms_changed.value_portfolio(SEQUENTIAL)?;
+    let respread = 1_000_000.0 * discount(&terms_changed, "JPY-CSA", payment)?;
+    assert_jpy_pv(&terms_changed, &terms_report, "CSA", respread)?;
+    assert_close!(respread, alternate, REL, JPY_ABS);
     assert_close!(
         corrected.value_portfolio(SEQUENTIAL)?.total?.amount(),
         corrected_report.total?.amount(),
diff --git a/crates/rhodium/tests/pricing_methods.rs b/crates/rhodium/tests/pricing_methods.rs
new file mode 100644
index 0000000000000000000000000000000000000000..7caaba68e732502bc4868f9ba75114feb387dd93
--- /dev/null
+++ b/crates/rhodium/tests/pricing_methods.rs
@@ -0,0 +1,269 @@
+//! ワークスペースの取引の評価方式を、価格設定の方式の表（`[[methods]]`）で選ぶ。
+//!
+//! 方式の表は取引の標準商品の種類で引く（設計書「方式」）。金利先物の既定は解析プライサーで、
+//! 表に `linear` を書けば値洗いの列の IR をリニア評価する。二つの評価は値洗いの口座の測度で
+//! 同じ式になる（`futures_margin_ir_vs_analytic`）ので、手数料を含めて PV が一致する。
+//! 変更した取引（一部解約）は標準商品の条件と一致しないので種類を持たず、表に当たらない。
+
+use std::{
+    fs,
+    io::{self, Write as _},
+    path::Path,
+};
+
+use rhodium::{
+    prelude::ExecutionContext,
+    session::{PortfolioValuation, Session, Valuation},
+    workspace::Workspace,
+};
+use rhodium_core::{Error, Result, currency::Currency};
+use rhodium_pricing::config::{AnalyticPricer, Method};
+use rhodium_testkit::assert_close;
+
+const AS_OF: &str = "2026-11-24T08:00:00Z";
+const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
+const REL: f64 = 1e-12;
+/// USD の最小単位（セント）の 1e-6 倍。
+const USD_ABS: f64 = 1e-8;
+
+#[expect(clippy::needless_pass_by_value, reason = "map_err に渡す変換関数")]
+fn io_error(error: io::Error) -> Error {
+    Error::Internal(error.to_string())
+}
+
+fn copy_dir(from: &Path, to: &Path) -> io::Result<()> {
+    for entry in fs::read_dir(from)? {
+        let entry = entry?;
+        let target = to.join(entry.file_name());
+        if entry.file_type()?.is_dir() {
+            fs::create_dir(&target)?;
+            copy_dir(&entry.path(), &target)?;
+        } else {
+            fs::copy(entry.path(), target)?;
+        }
+    }
+    Ok(())
+}
+
+fn append(path: &Path, suffix: &str) -> Result<()> {
+    fs::OpenOptions::new()
+        .append(true)
+        .open(path)
+        .and_then(|mut file| file.write_all(suffix.as_bytes()))
+        .map_err(io_error)
+}
+
+/// M0 のサンプルに、USD の金利先物の証券、清算指標、限月、校正の quote と LGM のマーク、
+/// 取引（手数料付き、手数料なし、一部解約）を足したワークスペースを書く。
+fn futures_workspace(root: &Path, extra_pricing: &str) -> Result<()> {
+    copy_dir(&sample(), root).map_err(io_error)?;
+    add_futures_refdata(root)?;
+    add_futures_market(root, extra_pricing)?;
+    write_futures_trades(root)
+}
+
+#[expect(
+    clippy::disallowed_macros,
+    reason = "リポジトリのサンプルの場所をテストのソースの位置から決める"
+)]
+fn sample() -> std::path::PathBuf {
+    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m0-ois")
+}
+
+/// 先物の証券、原指標と清算指標、限月の慣行を参照データに足す。
+fn add_futures_refdata(root: &Path) -> Result<()> {
+    fs::write(
+        root.join("refdata/securities.toml"),
+        r#"
+version = 1
+[[security]]
+id = "USD-FUT-NOV"
+kind = "listed_future"
+currency = "USD"
+unit = "contract"
+valuation = "contract"
+"#,
+    )
+    .map_err(io_error)?;
+    append(
+        &root.join("refdata/indices.toml"),
+        r#"
+
+[[index]]
+id = "USD-TERM"
+calendar = "USNY"
+publication = { business_days = 0, calendar = "USNY", time = "10:00", zone = "America/New_York" }
+observation_time = { time = "10:00", zone = "America/New_York" }
+kind = { term_rate = { ccy = "USD", day_count = "ACT360", tenor = "1W", value_offset = "0D", bdc = "following", end_of_month = "ignore" } }
+
+[[index]]
+id = "USD-FUT-NOV"
+calendar = "USNY"
+publication = { business_days = 0, calendar = "USNY", time = "10:00", zone = "America/New_York" }
+observation_time = { time = "10:00", zone = "America/New_York" }
+kind = { rate_futures_settlement = { contract = "USD-FUT-NOV", underlying = "USD-TERM" } }
+"#,
+    )?;
+    append(
+        &root.join("refdata/conventions.toml"),
+        r#"
+
+[[convention]]
+id = "USD-FUTURES"
+[convention.kind.rate_future]
+multiplier = 1000.0
+settlement_lag = "0D"
+settlement_time = { time = "10:00", zone = "America/New_York" }
+contracts = { "2026-11" = { settlement_index = "USD-FUT-NOV", final_clearing = 2026-11-30, rate = { term = { fixing = 2026-11-30 } } } }
+"#,
+    )
+}
+
+/// USD の曲線を先物の価格で校正し、LGM のマークと価格設定を足す。
+fn add_futures_market(root: &Path, extra_pricing: &str) -> Result<()> {
+    let path = root.join("market/definition.toml");
+    let market = fs::read_to_string(&path).map_err(io_error)?.replace(
+        "quote = \"USD.SOFR.OIS.1M\"\nkind = \"rate\"\ntarget = \"par_rate\"\ntemplate = { ois = { convention = \"USD-SOFR-OIS\", tenor = \"1M\" } }",
+        "quote = \"USD.FUT.NOV\"\nkind = \"futures_price\"\ntarget = \"price\"\ntemplate = { rate_future = { convention = \"USD-FUTURES\", delivery_month = \"2026-11\" } }",
+    );
+    fs::write(&path, market).map_err(io_error)?;
+    append(
+        &path,
+        "\n[lgm.USD]\nmean_reversion = [\"USD.LGM.MR\"]\nvolatility = [\"USD.LGM.VOL\"]\n",
+    )?;
+    // 先物は USD 建てなので、FX を使わずに集計できるよう報告通貨を USD にする（JPY の取引は換算の
+    // ペアが無く除外されるが、ここでは先物の行だけを読む）。
+    let pricing = root.join("config/pricing.toml");
+    let text = fs::read_to_string(&pricing).map_err(io_error)?.replace(
+        "reporting_currency = \"JPY\"",
+        "reporting_currency = \"USD\"",
+    );
+    fs::write(
+        &pricing,
+        format!("{text}\nUSD-TERM = \"USD-SOFR\"\n{extra_pricing}"),
+    )
+    .map_err(io_error)?;
+    // 凸性の補正が効くよう LGM のボラを 1% にする。
+    fs::write(
+        root.join("quotes/snapshot.csv"),
+        "quote_id,value\nJPY.TONA.OIS.1M,45bp\nUSD.FUT.NOV,0.96\nUSD.LGM.MR,0.0\nUSD.LGM.VOL,0.01\n",
+    )
+    .map_err(io_error)
+}
+
+/// 手数料付き、手数料なし、一部解約の三つの先物の取引を書く。
+fn write_futures_trades(root: &Path) -> Result<()> {
+    let future = r#"type = "rate_future"
+settlement_index = "USD-FUT-NOV"
+contracts = 2
+trade_price = 0.955
+day_count = "ACT360""#;
+    fs::write(
+        root.join("trades/futures.toml"),
+        format!(
+            r#"
+version = 1
+[[trade]]
+id = "FUT-FEE"
+margin_set = "USD-CASH"
+book = "USD-FUTURES"
+trade_date = 2026-11-24
+fees = [{{ amount = "100 USD", side = "pay", at = 2026-11-27T15:00:00Z }}]
+[trade.product]
+{future}
+
+[[trade]]
+id = "FUT-PLAIN"
+margin_set = "USD-CASH"
+book = "USD-FUTURES"
+trade_date = 2026-11-24
+[trade.product]
+{future}
+
+[[trade]]
+id = "FUT-HALF"
+margin_set = "USD-CASH"
+book = "USD-FUTURES"
+trade_date = 2026-11-24
+[trade.product]
+type = "amended"
+from = {AS_OF}
+after = {{ scale = 0.5 }}
+[trade.product.before]
+{future}
+"#
+        ),
+    )
+    .map_err(io_error)
+}
+
+fn portfolio(extra_pricing: &str) -> Result<PortfolioValuation> {
+    let dir = tempfile::tempdir().map_err(io_error)?;
+    futures_workspace(dir.path(), extra_pricing)?;
+    let workspace = Workspace::load(dir.path())?;
+    let as_of = workspace.resolve_as_of(AS_OF)?;
+    Session::new(workspace, as_of, SEQUENTIAL)?.value_portfolio(SEQUENTIAL)
+}
+
+fn valuation<'a>(report: &'a PortfolioValuation, id: &str) -> Result<&'a Valuation> {
+    report
+        .results
+        .iter()
+        .find(|(trade, _)| trade.as_str() == id)
+        .ok_or_else(|| Error::Internal(format!("missing result for {id}")))?
+        .1
+        .as_ref()
+        .map_err(Clone::clone)
+}
+
+fn usd(report: &PortfolioValuation, id: &str) -> Result<f64> {
+    Ok(valuation(report, id)?.pv.get(Currency::USD))
+}
+
+#[test]
+fn workspace_rate_future_is_priced_analytically_by_default() -> Result<()> {
+    // 方式の表が空なら、標準商品の金利先物は解析プライサーで評価する。手数料は解析の入力にも入り、
+    // 手数料の有無の差は手数料の受渡（100 USD 払い）の PV だけになる。
+    let report = portfolio("")?;
+    for id in ["FUT-FEE", "FUT-PLAIN"] {
+        assert_eq!(
+            valuation(&report, id)?.method,
+            Method::Analytic(AnalyticPricer::RateFuture),
+            "{id}"
+        );
+    }
+    let fee = usd(&report, "FUT-FEE")? - usd(&report, "FUT-PLAIN")?;
+    // 割引因子は 1 に近いので、手数料の PV は −100 USD に近く、それより絶対値が小さい。
+    assert!(fee < -99.9 && fee > -100.0, "{fee}");
+    Ok(())
+}
+
+#[test]
+fn workspace_methods_table_selects_linear_for_rate_futures() -> Result<()> {
+    // 表に linear を書くと同じ取引を値洗いの列の IR でリニア評価し、手数料を含めて解析の PV と一致する。
+    let analytic = portfolio("")?;
+    let linear = portfolio("\n[[methods]]\nproduct = \"rate_future\"\nmethod = \"linear\"\n")?;
+    for id in ["FUT-FEE", "FUT-PLAIN"] {
+        assert_eq!(valuation(&linear, id)?.method, Method::Linear, "{id}");
+        let pv = usd(&linear, id)?;
+        assert!(pv.abs() > 1.0, "{id}: {pv}");
+        assert_close!(pv, usd(&analytic, id)?, REL, USD_ABS);
+    }
+    Ok(())
+}
+
+#[test]
+fn amended_rate_future_trade_uses_the_linear_method() -> Result<()> {
+    // 一部解約した先物は種類を持たないので方式の表に当たらず、変更後の IR をリニア評価する。
+    // 効力の瞬間（評価時点）以後の値洗いをすべて半分にするので、PV は元の先物の半分。
+    let report = portfolio("")?;
+    let half = valuation(&report, "FUT-HALF")?;
+    assert_eq!(half.method, Method::Linear);
+    assert_close!(
+        half.pv.get(Currency::USD),
+        0.5 * usd(&report, "FUT-PLAIN")?,
+        REL,
+        USD_ABS
+    );
+    Ok(())
+}
diff --git a/crates/rhodium/tests/standard_product_io.rs b/crates/rhodium/tests/standard_product_io.rs
index 0bf1f377c896be6a4e88239a241967417ab66249..087324bd7759c361aae6f9814936f99ddc22d46b 100644
--- a/crates/rhodium/tests/standard_product_io.rs
+++ b/crates/rhodium/tests/standard_product_io.rs
@@ -238,3 +238,82 @@ far = { rate = 151.0, settlement = 2026-01-07T16:00:00Z }
     ));
     Ok(())
 }
+
+#[test]
+fn m1_products_reject_invalid_terms() {
+    // M1 の商品の TOML の不正な入力は、読み込みの境界で InvalidData にする（参照データを要る照合は
+    // IR を作るときに行い、ここでは型付きの値が作れない条件を確かめる）。
+    let fra = r#"
+type = "fra"
+notional = "100000 USD"
+index = "USD-TERM"
+fixing = 2026-01-06
+fixed_rate = "4%"
+"#;
+    let future = r#"
+type = "rate_future"
+settlement_index = "USD-RATE-FUT-JAN"
+contracts = 2
+trade_price = 0.96
+day_count = "ACT360"
+"#;
+    let forward = r#"
+type = "fx_forward"
+pair = "USDJPY"
+notional = "1000 USD"
+rate = 150.0
+settlement = 2026-01-07T16:00:00Z
+"#;
+    let mtm = r#"
+type = "swap"
+
+[[legs]]
+side = "receive"
+notional = { fixed = "150_000 JPY" }
+coupon = { fixed = { rate = "1%" } }
+day_count = "ACT365F"
+schedule = { explicit = { periods = [{ accrual = { start = 2026-01-07, end = 2026-04-07 }, payment = 2026-04-07 }] } }
+exchange = "initial_and_final"
+
+[[legs]]
+side = "pay"
+notional = { fx_reset = { base = "150_000 JPY", fixing = "USDJPY-WMR" } }
+coupon = { fixed = { rate = "3%" } }
+day_count = "ACT360"
+schedule = { explicit = { periods = [{ accrual = { start = 2026-01-07, end = 2026-04-07 }, payment = 2026-04-07 }] } }
+exchange = "initial_and_final"
+"#;
+    for invalid in [
+        // FRA の固定金利は尺度付き、元本は正。
+        fra.replace("\"4%\"", "0.04"),
+        fra.replace("100000 USD", "-100000 USD"),
+        fra.replace("fixing = 2026-01-06", "fixing = 2026-01-06T00:00:00Z"),
+        // 先物の枚数は正の整数、価格は有限、日数基準は既知。
+        future.replace("contracts = 2", "contracts = 0"),
+        future.replace("contracts = 2", "contracts = -2"),
+        future.replace("trade_price = 0.96", "trade_price = \"96\""),
+        future.replace("ACT360", "ACT999"),
+        // FX フォワードの元本はペアの base の通貨、レートは正、受渡は UTC の瞬間。
+        forward.replace("1000 USD", "1000 JPY"),
+        forward.replace("rate = 150.0", "rate = -150.0"),
+        forward.replace("2026-01-07T16:00:00Z", "2026-01-07"),
+        // FX リセットの脚は元本の交換を要する。リセットの脚を二本持つ XCCY は作れない。
+        mtm.replacen("exchange = \"initial_and_final\"\n", "", 2),
+        mtm.replace(
+            "notional = { fixed = \"150_000 JPY\" }",
+            "notional = { fx_reset = { base = \"1_000 USD\", fixing = \"USDJPY-WMR\" } }",
+        ),
+    ] {
+        assert!(
+            matches!(
+                read_standard_product_toml(&invalid),
+                Err(Error::InvalidData(_))
+            ),
+            "{invalid}"
+        );
+    }
+    // 両方の形は正しい入力として読める（不正な入力が別の理由で落ちていないことの対照）。
+    for valid in [fra, future, forward, mtm] {
+        assert!(read_standard_product_toml(valid).is_ok(), "{valid}");
+    }
+}
diff --git a/crates/rhodium/tests/support/desk.rs b/crates/rhodium/tests/support/desk.rs
index 4034856e46d0833aba7cd350fe8281a9d04b883a..0bef0c61225a03492fd15b38dd8bfae9b785a01c 100644
--- a/crates/rhodium/tests/support/desk.rs
+++ b/crates/rhodium/tests/support/desk.rs
@@ -133,6 +133,10 @@ fn add_months(date: Date, months: i32) -> Result<Date> {
 }
 
 // 週末なら次の月曜に進める（週末だけのカレンダーの Following）。
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#196 で各 variant を網羅する"
+)]
 fn next_weekday(date: Date) -> Result<Date> {
     match date.weekday() {
         Weekday::Sat => date.checked_add_days(2),
@@ -384,6 +388,13 @@ payment = { anchor = "period_end", offset = "2BD", calendar = "USNY", bdc = "fol
 [convention.kind.ois.float]
 index = "USD-SOFR"
 observation = { days = { business_days = "USNY" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } }
+
+# USDJPY の受渡慣行。FX のマークを足したテストが評価時点の換算レートに使う。
+[[fx]]
+pair = "USDJPY"
+spot_lag = "two_business_days"
+settlement = { time = "17:00:00", zone = "Asia/Tokyo" }
+pip_exponent = 2
 "#;
 
 const AGREEMENTS: &str = r#"version = 1
diff --git a/docs/benchmarks/README.md b/docs/benchmarks/README.md
index 0f04ca4bac8fbc3dae4ebe7e3e309d264c8b2c70..b9d10b44b072216d8d3b6ab826c1db6fa7687d86 100644
--- a/docs/benchmarks/README.md
+++ b/docs/benchmarks/README.md
@@ -6,4 +6,4 @@
 | マイルストーン | 記録 | 計測のコマンド |
 |---|---|---|
 | M0 | [m0.md](m0.md) | `cargo bench --locked -p rhodium --bench m0` |
-| M1 | [#180 の計測](issue-180.md) | `cargo bench --locked -p rhodium --bench active_market`、`issue_180_workspace` |
+| M1 | [m1.md](m1.md)、[M1-7 の計測](m1-7.md)、[#180 の計測](issue-180.md) | `cargo bench --locked -p rhodium --bench m1`、`--bench m0 -- m1/quote_update_10k`、`--bench active_market`、`issue_180_workspace` |
diff --git a/docs/benchmarks/m1.md b/docs/benchmarks/m1.md
new file mode 100644
index 0000000000000000000000000000000000000000..f365aa82dc0e5aba741879d6fb8fd36988e36d62
--- /dev/null
+++ b/docs/benchmarks/m1.md
@@ -0,0 +1,99 @@
+# M1 の計測
+
+> M1 の性能の目安の記録です（[ロードマップ](../design/12-roadmap.md#m1-金利リニアの完成)、Issue #37）。
+> 最初の検証データセットの市場（[`examples/m1-rates`](../../examples/m1-rates/README.md)）の校正と、その取引を繰り返した 1 万件と 10 万件のデスクの読み込み、準備、PV、quote デルタ、quote 1 本の更新を測りました。
+> どの項目も目安を大きく下回り、AAD の quote デルタは PV の 1.4 倍で済みました。
+
+## 環境
+
+| 項目 | 内容 |
+|---|---|
+| 日付 | 2026-11-26 |
+| コミット | b0a367d（`bench: measure the M1 validation dataset`）。計測した HEAD は、その後の文書だけのコミットを含む f6842eb |
+| マシン | Intel Core i7-12700F（12 コア: P コア 8・E コア 4、20 スレッド、最大 4.9 GHz）、メモリ 31 GiB |
+| OS | NixOS 26.05、Linux 6.18.46 |
+| ツールチェーン | rustc 1.100.0、`cargo bench` の既定のプロファイル（release、opt-level 3） |
+| 計測 | criterion 0.8.2、各 10 標本、`cargo bench --locked -p rhodium --bench m1` |
+
+設計書の目安は 32 コアのマシンを「1 台」とします（[性能の目安](../design/10-io-and-execution.md#性能の目安)）。
+このマシンは 12 コア 20 スレッドなので、「1 台」の値には 20 スレッドの値を、「1 コア」の値には 1 スレッドの値を使います。
+
+## ワークロード
+
+- 市場: `examples/m1-rates` の市場の定義そのものです。
+  JPY-TONA（3 回の会合の 4 区間の階段、1 年で接続する 2Y〜10Y の長期の節点、2026 年末の跳び）、TIBOR 6M と 3M（TONA へのベーシス、各 4 ピラー）、USD-SOFR（7 ピラー）、USDJPY の外貨担保曲線（TN・3M・6M のポイントと 1Y・2Y・5Y の MtM XCCY）、FX と LGM のマークの 9 グループで、quote は 35 本です。
+- デスク: サンプルの 18 件の取引（預金、FRA、金利先物、IRS、ベーシス、lookback・observation shift・lockout の OIS、一部解約の IRS、FX フォワード、FX スワップ、MtM XCCY、固定元本の XCCY、CSA の 4 つの場合の SOFR の OIS と 25bp の CSA の IRS）を、ID だけを変えて 1 万件と 10 万件まで繰り返したものです。
+  同じ条件の繰り返しなので、商品の多様さはサンプルの 18 種類に限られます。
+- 評価時点: 2026-11-24 17:00 JST。
+
+| 計測の名前 | 測る処理 | スレッド |
+|---|---|---|
+| `m1_rates/calibration` | 取引の無いデスクの `Session::new`。9 グループの校正商品の準備と校正 | 1 |
+| `m1_rates/workspace_load_{10k,100k}` | `Workspace::load`。ファイルの解析と取引の IR への変換 | 1（読み込みは逐次） |
+| `m1_rates/session_new_{10k,100k}` | `Session::new`。校正と取引の準備 | 20 |
+| `m1_rates/pv_{10k,100k}` | 準備済みのセッションの `value_portfolio`（取引ごとの PV と報告通貨 JPY への換算と合計） | 20 |
+| `m1_rates/quote_delta_{10k,100k}` | 準備済みのセッションの `quote_delta_reporting`（θ 勾配の AAD、報告通貨への換算の連鎖律、校正の随伴） | 20 |
+| `m1_rates/quote_update_{10k,100k}` | `JPY.TONA.OIS.10Y` を 1.32% から 1.33% にする `with_quotes`（TONA の長期、TIBOR 6M・3M、外貨担保曲線の 4 グループの再校正）と、その後の `value_portfolio` と `quote_delta_reporting` | 20 |
+
+## 結果
+
+値は criterion の推定値（10 標本の中央の推定）で、括弧は 95% の区間です。
+
+| 計測 | 1 万件 | 10 万件 |
+|---|---|---|
+| `m1_rates/calibration`（1 スレッド、取引なし） | 30.6 ms（30.3 ms〜30.9 ms） | — |
+| `m1_rates/workspace_load`（1 スレッド） | 930 ms（929 ms〜931 ms） | 9.37 s（9.31 s〜9.43 s） |
+| `m1_rates/session_new`（20 スレッド） | 429 ms（387 ms〜469 ms） | 4.53 s（4.22 s〜4.75 s） |
+| `m1_rates/pv`（20 スレッド） | 23.0 ms（23.0 ms〜23.1 ms） | 233 ms（233 ms〜234 ms） |
+| `m1_rates/quote_delta`（20 スレッド） | 31.7 ms（30.2 ms〜33.1 ms） | 327 ms（315 ms〜338 ms） |
+| `m1_rates/quote_update`（20 スレッド） | 63.2 ms（62.7 ms〜63.5 ms） | 576 ms（567 ms〜589 ms） |
+
+- 計測の全体（criterion の反復と準備を含む）は約 7 分でした。
+- PV と quote デルタは取引の数にほぼ比例し（1 万件から 10 万件で 10.1 倍と 10.3 倍）、quote デルタは PV の 1.38 倍（1 万件）と 1.40 倍（10 万件）です。
+  デルタは取引ごとの θ 勾配の AAD（読む θ ブロックだけを変数にする active set）と、報告通貨への換算の連鎖律、校正の随伴の一回の逆解きからなり、随伴の費用は取引の数によりません。
+- quote の更新は、4 グループ（TONA の長期、TIBOR 6M・3M、外貨担保曲線）の再校正と PV・デルタの和で、PV とデルタの和（1 万件で 54.7 ms）との差 8.5 ms が再校正の時間の目安です。
+- 端から端までの時間（読み込み、校正、準備、評価、CSV の出力）を、同じ 10 万件のデスクで release の `rh` に GNU time を当てて測りました。
+  `--threads` は `Workspace::load_with_execution` の取引のデコードと IR への変換も分けます。
+
+  | コマンド | 1 スレッド | `--threads 20` |
+  |---|---|---|
+  | `rh price` | 25.1 s | 8.86 s（最大の常駐メモリ 9.4 GiB） |
+  | `rh risk` | 26.8 s | 9.00 s（最大の常駐メモリ 9.4 GiB） |
+
+  `rh price` の CSV は 1 スレッドと 20 スレッドでバイト単位で同じでした（12 万 2,222 行。FX と XCCY の取引は 2 通貨の行を持ちます）。
+
+## 目安との比較
+
+| 目安（[性能の目安](../design/10-io-and-execution.md#性能の目安)） | M1 の値 | 評価 |
+|---|---|---|
+| JPY と USD の全曲線の校正（最初の検証データセットの市場の定義）、100 ms 未満（1 コア） | 30.6 ms（9 グループ、35 quote） | 目安内 |
+| quote 1 本の変更による再校正と、取引 1 万件のデスクの PV とデルタの再計算、1 秒未満（1 台） | 63.2 ms | 目安内 |
+| 線形な取引 10 万件の PV と quote デルタ、10 秒未満（1 台） | 準備済みのセッションで PV 233 ms とデルタ 327 ms の和 560 ms。読み込みと準備を含む端から端までは `rh risk --threads 20` で 9.00 s | 目安内。準備済みのセッションの PV とデルタは目安の 6%。端から端までの時間とメモリ（9.4 GiB）の大半は、取引のデコードと IR への変換（逐次で 9.4 s）と準備が占める。これらの縮小は M6 の目安の全項目の計測で見直す |
+| AAD によるリスクの費用、PV の 5 倍以内（リニア評価と解析プライサー） | 1.38 倍（1 万件）、1.40 倍（10 万件） | 目安内 |
+
+## M0 の読み込みと校正への影響
+
+M1 では、ワークスペースの取引に方式の表を当てるため、読み込みの結果に IR と一緒に取引（`Trade`）を持たせ、校正に予想の測度の記録と最小二乗の分岐を足しました。
+これが M0 のワークロードを遅くしていないことを、同じ日に同じマシンで、分岐元の main（81adcee）と計測した HEAD（f6842eb）を交互に 2 回ずつ測って確かめました（`cargo bench --locked -p rhodium --bench m0 -- 'm0/(calibration|workspace_load_10k)$'`）。
+
+| 計測 | main（1 回目、2 回目） | HEAD（1 回目、2 回目） |
+|---|---|---|
+| `m0/calibration` | 65.4 ms、64.6 ms | 64.2 ms、63.6 ms |
+| `m0/workspace_load_10k` | 2.195 s、2.187 s | 2.150 s、2.150 s |
+
+差は 2% 以内で、どちらも HEAD の方が短く、遅くなっていません。
+M0 の記録（[m0.md](m0.md)）の読み込み 10.27 s・校正 87.1 ms より短いのは、記録の後に main に入った改善（読み込みは [#180 の計測](issue-180.md) の観測の列の共有など）によります。
+
+## レビュー対応後の確認
+
+PR #192 のレビューへの対応（期間をまたぐ先物の 1 日、FX fixing の受渡、最小二乗の単位と採用条件、FRA の恒等式の式ごとの判定、予想の測度の根拠など）は、校正、準備、評価の経路に小さく触れます。
+同じマシンで d2317a1 の校正と 1 万件のデスクを測り直し（`cargo bench --locked -p rhodium --bench m1 -- 'm1_rates/(calibration|pv_10k|quote_delta_10k|quote_update_10k)$'`）、上の記録と変わらないことを確かめました。
+
+| 計測 | 上の記録 | d2317a1 |
+|---|---|---|
+| `m1_rates/calibration`（1 スレッド） | 30.6 ms | 30.8 ms（30.6 ms〜30.8 ms） |
+| `m1_rates/pv_10k`（20 スレッド） | 23.0 ms | 23.2 ms（22.9 ms〜23.5 ms） |
+| `m1_rates/quote_delta_10k`（20 スレッド） | 31.7 ms | 30.9 ms（30.7 ms〜31.1 ms） |
+| `m1_rates/quote_update_10k`（20 スレッド） | 63.2 ms | 63.5 ms（63.2 ms〜63.9 ms） |
+
+差は 3% 以内で、区間が重なるか接しています。
diff --git a/docs/design/01-overview.md b/docs/design/01-overview.md
index 41fb763094526be54918e874a37d3086772222af..74963b8a411f084f7f00a4ac6c4d59e2429379f9 100644
--- a/docs/design/01-overview.md
+++ b/docs/design/01-overview.md
@@ -25,6 +25,7 @@ rhodium は商品表現（コア IR）、market、価格設定を三つの計算
 - 計算: PV とキャッシュフロー、quote への 1 次と 2 次の感応度、ヘッジへの付け替え、シナリオ、theta と P&L 説明、ゴールシーク、SIMM、FRTB-SA、SA-CCR、XVA（CVA、DVA、FCA/FBA、ColVA、MVA、KVA）と XVA の感応度、配分、incremental XVA。
 - 評価方式: リニア評価、解析プライサー、モンテカルロ、AMC（回帰によるモンテカルロ）。
 - 最初の検証データセット: JPY（TONA、TIBOR）と USD（SOFR）の曲線、USDJPY の XCCY と FX、その金利商品と XVA。
+  金利商品までの部分は [M1 のサンプル](../../examples/m1-rates/README.md) にあります。
   通貨と資産は参照データと market の定義を足すだけで増やします。
 - どの順で実装するかは [ロードマップ](12-roadmap.md) が定めます。
 
@@ -119,10 +120,10 @@ flowchart LR
 
 固定金利と TONA の複利を交換する、開始済みのスワップの PV と、quote へのデルタを求めます。
 
-1. 市場の定義: `market.toml` は JPY-TONA 曲線（ln DF のノードと補間）と、TONA の OIS quote（`JPY.TONA.OIS.1Y` から `JPY.TONA.OIS.30Y`）を OIS の校正商品に結ぶ指定を持ちます。
+1. 市場の定義: `market/definition.toml` は JPY-TONA 曲線（ln DF のノードと補間）と、TONA の OIS quote（`JPY.TONA.OIS.1Y` から `JPY.TONA.OIS.30Y`）を OIS の校正商品に結ぶ指定を持ちます。
 2. 校正: `calibrate` が OIS のパーレートと quote の差を残差として、Newton 法で θ を解きます。
    解の記録として、残差の θ と quote に関するヤコビアンの分解を `CalibratedMarket` に残します。
-3. 商品: `trades.toml` の取引は `StandardProduct::Swap`（固定脚と TONA 複利の変動脚）です。
+3. 商品: `trades/*.toml` の取引は `StandardProduct::Swap`（固定脚と TONA 複利の変動脚）です。
    `to_product` が参照データのカレンダーと慣行から日程を作り、`Pay` と `Compounded` からなる IR にします。
 4. 準備: `prepare` は評価時点に公表済みの TONA の fixing を定数に置き換え、支払い済みのクーポンを除き、IR が線形なのでリニア評価の計画にします。
    公表前の観測は [参照データ](04-foundation/04-reference-data.md#fixingstore) の規則で予想か推定値にします。
diff --git a/docs/design/03-architecture.md b/docs/design/03-architecture.md
index 29bace1fb48c73c9491d39e9977395f98252dc9e..dc49026981e1e6305c2445ba7d009d3382dce98f 100644
--- a/docs/design/03-architecture.md
+++ b/docs/design/03-architecture.md
@@ -41,13 +41,16 @@ flowchart BT
 | crate | モジュール |
 |---|---|
 | core | `date`、`time`、`tenor`、`calendar`、`day_count`、`schedule`、`observation`、`currency`、`money`、`quantity`、`fx`、`rounding`、`id`、`refdata`、`fixing`、`brand`、`error`、`parallel` |
-| products | `ir`（`nodes`、`builder`、`validate`、`reduce`、`analysis`）、`dsl`（`syntax`、`check`、`lower`）、`standard`（`rates`、`fx`、`credit`、`equity`、`inflation`、`commodity`、`leg`、`template`、`wire`）、`trade`、`agreement` |
-| market | `quote`、`structure`、`market`、`view`、`curve`、`credit`、`inflation`、`asset`、`vol`、`interp` |
-| pricing | `config`、`context`、`resolve`、`definition`、`calibration`、`prepare`、`linear`、`par`、`analytic`、`model`、`mc`、`amc`、`valuation` |
+| products | `ir`（`nodes`、`builder`、`validate`、`reduce`、`analysis`、`amendment`、`fallback`）、`dsl`（`syntax`、`check`、`lower`）、`standard`（`rates`、`fra`、`future`、`fx`、`credit`、`equity`、`inflation`、`commodity`、`leg`、`template`、`wire`）、`trade`、`agreement` |
+| market | `quote`、`structure`、`market`、`view`、`curve`、`fx`、`credit`、`inflation`、`asset`、`vol`、`interp` |
+| pricing | `config`、`context`、`resolve`、`definition`、`calibration`、`prepare`、`linear`（`rfr`、`fx`、`future`）、`par`、`analytic`、`model`、`mc`、`amc`、`valuation` |
 | risk | `vector`、`gradient`、`adjoint`、`delta`、`gamma`、`bump`、`scenario`、`theta`、`explain`、`goal_seek`、`rebase`、`simm`、`frtb`、`saccr` |
 | xva | `grid`、`exposure`、`collateral`、`closeout`、`im`、`metrics`、`sensitivity`、`allocation`、`incremental` |
 | rhodium | `prelude`、`workspace`、`session`、`io` |
 
+表は 1.0 までの配置で、まだ無いモジュールは担当するマイルストーンで足します。
+M1 の唯一の解析プライサー（金利先物）は、値洗いの列をリニア評価と共有するため `linear::future` に置き、M2 のスワップション・キャップ・FX オプションの解析プライサーとともに `analytic` へ分けます。
+
 ## 主要な型の地図
 
 型は、それを定義できる最も下位の crate に置きます。
@@ -66,7 +69,7 @@ flowchart BT
 | `CcyPair`、`FxRate<D>`、`FxForwardPoints`、`Rounding` | core::fx、core::rounding | [金額と量](04-foundation/03-money-and-quantities.md) |
 | `VolKind`、`CallPut` | core::quantity | [金額と量](04-foundation/03-money-and-quantities.md) |
 | `IndexId`、`CalendarId`、`ConventionId`、`EntityId`、`SecurityId`、`CommodityId` | core::id | [参照データ](04-foundation/04-reference-data.md) |
-| `RefData`、`IndexDefinition`、`IndexKind`、`PublicationRule`、`Convention`、`CreditKey` | core::refdata、core::id | [参照データ](04-foundation/04-reference-data.md) |
+| `RefData`、`IndexDefinition`、`IndexKind`、`PublicationRule`、`Convention`、`FxConvention`、`FxSpotLag`、`CreditKey` | core::refdata、core::id | [参照データ](04-foundation/04-reference-data.md) |
 | `FixingStore`、`ObservationPoint` | core::fixing | [参照データ](04-foundation/04-reference-data.md#fixingstore) |
 | `Brand<'id>` | core::brand | この文書 |
 | `Error`、`Result<T>`、`Warning` | core::error | [参照データ](04-foundation/04-reference-data.md#エラー) |
@@ -75,28 +78,28 @@ flowchart BT
 | `ScriptSource`、`ScriptLibrary`、DSL の構文木 | products::dsl | [DSL](05-products/02-dsl.md) |
 | `StandardProduct`、`ProductKind`、各商品の構造体、`Leg`、`InstrumentTemplate` | products::standard | [標準商品](05-products/03-standard-products.md) |
 | `TradeId`、`BookId`、`MarginSetId` | core::id | [取引](05-products/04-trades.md) |
-| `Trade`、`ProductSpec`、`Amendment`、`ExerciseRecord`、`Fee`、`FallbackTerms`、`FallbackSelection`、`TradeProductResolver`、`Portfolio` | products::trade | [取引](05-products/04-trades.md) |
+| `Trade`、`ProductSpec`、`Amendment`、`ExerciseRecord`、`Fee`、`FallbackTerms`、`FallbackSelection`、`TradeProductResolver`、`ResolvedTrade`、`Portfolio` | products::trade | [取引](05-products/04-trades.md) |
 | `Counterparty`、`NettingSet`、`MarginSet`、`CsaTerms`、`CashCollateralTerms`、`ImTerms` | products::agreement | [取引](05-products/04-trades.md) |
-| `QuoteId`、`QuotePattern`、`QuoteSet`、`QuoteKind`、`QuoteOverrides` | market::quote（`QuoteId` は quote を指す `Warning` が持てるよう core::id で定義し、market::quote から再公開） | [quote](06-market/01-quotes.md) |
+| `QuoteId`、`QuotePattern`、`QuoteSet`、`QuoteKind`、`QuoteOverride`、`QuoteOverrides` | market::quote（`QuoteId` は quote を指す `Warning` が持てるよう core::id で定義し、market::quote から再公開） | [quote](06-market/01-quotes.md) |
 | `MarketStructure`、`CurveId`、`VolId`、`ModelId`、`CurveSpec`、`CurveExpr`、`CurvePoint`、`ParamId`、`BlockId`、`ParamLayout` | market::structure、market::curve | [曲線と market](06-market/02-curves.md) |
 | `Market<R>`、`MarketShifts`、`MarketView<'v, R>`、`ResolvedCurve<K>`、`CurveRef<'v, K>`、`ScalarId`、`ScalarRef<'v>`、`Slots`、`BoundSlots<'v, R>`、`CurveAt<'v, K, R>`、`ForwardAt<'v, R>`、`GrowthAt<'v, R>`、`CcyScope<'c, 'v, R>`、`PayAt<'c, R>`、`Amount<'c, R>` | market::market、market::view | [曲線と market](06-market/02-curves.md) |
 | `FxSpec`、`FxAnchor`、`FxCollateral`、`FxRef<'v>`、`FxForwardPoint`、`FxForwardInterval` | market::fx | [曲線と market](06-market/02-curves.md#fx-のフォワード) |
-| `VolSpec`、`VolRef<'v>`、`VolPoint`、`VolAt<'v>`、`SmileModel`、`SmileDynamics` | market::vol | [ボラティリティ](06-market/04-volatility.md) |
-| `PricingConfig`、`FundingCurves`、`Method`、`Context`、`Requirement`、`VolKey`、`CollateralKey` | pricing::config、pricing::context、pricing::resolve | [価格計算](07-pricing/01-pricing.md) |
-| `MarketDefinition`、`InstrumentBinding`、`CalibrationTarget`、`QuoteTransform`、`FxDefinition` | pricing::definition | [校正](06-market/03-calibration.md) |
-| `CalibratedMarket`、`CalibrationGroup`、`GroupOutcome`、`GroupSolution` | pricing::calibration | [校正](06-market/03-calibration.md) |
-| `Prepared`、`AnalyticPricer`、`Valuation`、`PortfolioValuation` | pricing::prepare、pricing::valuation | [価格計算](07-pricing/01-pricing.md) |
-| `CashflowPlan` | pricing::linear | [リニアと解析](07-pricing/02-linear-and-analytic.md) |
+| `VolSpec`、`VolRef<'v, 'a>`、`VolPoint`、`VolAt<'v, 'a>`、`SmileModel`、`SmileDynamics` | market::vol | [ボラティリティ](06-market/04-volatility.md) |
+| `PricingConfig`、`FundingCurves`、`Method`、`AnalyticPricer`、`LinearConvention`、`Context`、`Requirement`、`VolKey`、`CollateralKey` | pricing::config、pricing::context、pricing::resolve | [価格計算](07-pricing/01-pricing.md) |
+| `MarketDefinition`、`InstrumentBinding`、`CalibrationTarget`、`QuoteTransform`、`ElementSelector`、`StepPenalty`、`FxDefinition`、`FxImpliedCurve`、`FxImpliedLongEnd` | pricing::definition | [校正](06-market/03-calibration.md) |
+| `CalibratedMarket`、`CalibrationGroup`、`GroupOutcome`、`GroupSolution`、`HessianMode` | pricing::calibration | [校正](06-market/03-calibration.md) |
+| `Prepared`、`PrepareInput`、`Valuation`、`PortfolioValuation` | pricing::prepare、pricing::valuation | [価格計算](07-pricing/01-pricing.md) |
+| `CashflowPlan`（crate の内部。公開するのは `Prepared` だけ） | pricing::linear | [リニアと解析](07-pricing/02-linear-and-analytic.md) |
 | `PathModel`、`ModelSpec`、`Lgm`、`CrossCurrencyLgm` | pricing::model | [モデル](07-pricing/03-models.md) |
 | `McSettings`、`Kernel`、`PathBatch`、`AmcSettings` | pricing::mc、pricing::amc | [モンテカルロ](07-pricing/04-monte-carlo.md) |
-| `RiskVector<K>`、`ParamGradient`、`HedgeId`、`DeltaOptions`、`GammaMatrix<K>` | risk::vector、risk::delta、risk::gamma | [感応度](08-risk/01-sensitivities.md) |
+| `RiskVector<K>`、`ParamGradient`、`DeltaReport`、`HedgeId`、`DeltaOptions`、`GammaMatrix<K>` | risk::vector、risk::delta、risk::gamma | [感応度](08-risk/01-sensitivities.md) |
 | `Scenario`、`QuoteShift`、`MarketShift`、`BlockSelector`、`AssumedEvent`、`RollMode`、`ExerciseAssumption` | risk::scenario、risk::theta | [シナリオと P&L](08-risk/02-scenarios-and-pnl.md) |
 | `HedgeSet`、`HedgeEntry`、`RebaseMetrics`、`Rebased`、`SimmKey`、`CrifRow` | risk::rebase、risk::simm | [付け替えと規制](08-risk/03-hedging-and-regulatory.md) |
 | `TimeGrid`、`GridSpec`、`GridView<'g>`、`GridIx<'g>` | xva::grid | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md) |
 | `MporModel`、`CollateralBalance`、`PendingTransfer` | xva::collateral、xva::closeout | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md) |
 | `XvaMetric`、`XvaSettings`、`ImModel`、`XvaResult`、`XvaReport` | xva::metrics、xva::im | [XVA の指標](09-xva/02-metrics.md) |
 | `XvaRun`、`IncrementalXva` | xva::incremental | [XVA のリスク](09-xva/03-xva-risk.md) |
-| `Workspace`、`Session`、`Selection` | rhodium::workspace、rhodium::session | [入出力と実行](10-io-and-execution.md) |
+| `Workspace`、`Session`、`QuoteDeltaReport`、`Selection` | rhodium::workspace、rhodium::session | [入出力と実行](10-io-and-execution.md) |
 
 ## 型の規約
 
@@ -130,7 +133,7 @@ impl<R: Real> Market<R> {
 
 | lifetime | スコープの開き方 | 守ること | 文書 |
 |---|---|---|---|
-| `'v` | `Market::with_view` | `CurveRef<'v, K>`、`FxRef<'v>`、`VolRef<'v>` と、`bind_slots` が曲線・面と点を組にした `CurveAt<'v, K, R>`、`ForwardAt<'v, R>`、`GrowthAt<'v, R>`、`VolAt<'v>` はその view の中でだけ使え、読み出しは失敗しない | [曲線と market](06-market/02-curves.md) |
+| `'v` | `Market::with_view` | `CurveRef<'v, K>`、`FxRef<'v>`、`VolRef<'v, 'a>` と、`bind_slots` が曲線と点を組にした `CurveAt<'v, K, R>`、`ForwardAt<'v, R>`、`GrowthAt<'v, R>`、`bind_vol` が面と準備済みの点を組にした `VolAt<'v, 'a>` はその view の中でだけ使える。曲線の読み出しは失敗せず、ボラの読み出しはストライクが展開の定義域を外れるとエラーを返す（[ボラティリティ](06-market/04-volatility.md#読み出し)） | [曲線と market](06-market/02-curves.md) |
 | `'c` | `MarketView::ccy_scope` | `Amount<'c, R>` は同じ通貨スコープの値とだけ足せ、`PayAt<'c, R>`（スコープの割引曲線に結んだ支払いの点）はそのスコープの `pv` にだけ渡せる | [金額と量](04-foundation/03-money-and-quantities.md) |
 | `'b` | `Contract::build` | `ExprIx<'b>`、`ContractIx<'b>` は同じ builder の中でだけ組み合わせられる | [コア IR](05-products/01-ir.md) |
 | `'g` | `TimeGrid::with_view` | `GridIx<'g>` で読むパスの値は同じグリッドの添字でだけ読める | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md) |
@@ -193,7 +196,7 @@ impl<R: Real> Market<R> {
 - ポートフォリオの計算は取引ごと（XVA はネッティングセットごと）の `Result` を返し、一件の失敗で全体を止めません。
   集計は除いた対象の一覧を必ず伴います。
 - ライブラリのコードは利用者の入力で panic しません。
-  内部の不変条件の破れ（バグ）だけは `debug_assert!` で検出します。
+  内部の不変条件の破れ（バグ）は、まず型で表せないかを確かめ、表せないものだけを注記付きの `Error::Internal` で報告します（[参照データのエラー](04-foundation/04-reference-data.md#エラー)の三段の方針）。
 
 ## 依存とビルド
 
@@ -204,7 +207,8 @@ impl<R: Real> Market<R> {
   scioxide-core は `ThreadExecutor` のために `std` feature で使います。
   `std` は scioxide の初等関数の実装も切り替えるので、crate ごとに付けず `[workspace.dependencies]` で一度だけ指定し、どのビルドでも同じ実装にします。
   toml、csv、serde_json はファサードと CLI だけが使い、引数の解析の clap は CLI だけが使います。
-  開発用には trybuild（compile-fail テスト）、proptest（性質テスト）、criterion（ベンチマーク、`default-features = false`）、tempfile（ファイルのワークスペースのテスト）を使い、テストの外部 fixture の読み込みには serde_json と toml も使います。
+  開発用には proptest（性質テスト）、criterion（ベンチマーク、`default-features = false`）、tempfile（ファイルのワークスペースのテストと、型テストのケースごとの出力先）を使い、テストの外部 fixture の読み込みには serde_json と toml も使います。
+  testkit の型テストのランナーは、rustc と Cargo の JSON の出力と、期待する診断の `case.json` の解析に serde_json を使います。
 - 開発専用の crate `rhodium-testkit` が、テストの補助（`assert_close!`、`assert_mc!`、`check_aad_vs_bump`、`compile_fail_pair!`、決定的な時計 `TickClock`）を一か所に持ちます（[品質](11-quality.md)）。
   ライブラリの 8 crate はこれに dev-dependency としてだけ依存します。
 - feature は加算的にし、feature の有無で計算結果を変えません。
@@ -214,7 +218,13 @@ impl<R: Real> Market<R> {
 
 - `unsafe_code = "forbid"`、`clippy::all` を deny、`clippy::pedantic` を warn とし、数値コードに合わない項目（`similar_names`、`many_single_char_names`、`cast_precision_loss` など）は workspace で allow にします。
 - ライブラリのコードで `unwrap_used`、`expect_used`、`panic` を deny にし、テストでは許します。
+- `clippy::wildcard_enum_match_arm` を deny にし、enum を `_` や束縛の腕で受けません（variant を足したときに既定の腕へ黙って落ちるのを防ぐため）。
+  外部の enum に既定の腕が要る箇所と、担当の issue で直す既存の箇所だけは、外側の fn に理由付きの `#[expect]` を付けます。
+  既存の箇所の理由には担当の issue の番号を書くので、その issue が分岐を網羅すると `unfulfilled_lint_expectations` が expect の外し忘れを止めます。
+  この lint は `matches!`、let-else、タプルの `match` を見ないので、複数の variant にわたる分類はそれらで書かず、enum のメソッドの網羅した `match` に置きます。
 - `clippy.toml` の disallowed は `std::time::SystemTime::now`、`std::time::Instant::now`、`chrono::Utc::now`、`chrono::Local::now`、`std::env` の読み出しだけにします（[P1](02-principles.md#p1-純粋関数と明示的な入力)）。
+- 型や lint で止められない実装の誤りの形（ID の文字列の合成、許容誤差のリテラル、下流での参照データの照合し直しなど）は、`.githooks` の禁止パターンの表で止めます（規則の書き方と diff・tree のモードは[品質の CI](11-quality.md#ci)）。
+- lint の例外は reason 付きの `#[expect]` にし、効く範囲を最小の item に絞ります。`#[allow]` は、例外が要らなくなっても `unfulfilled_lint_expectations` で気づけないので避けます。
 - `missing_docs` は公開 API に warn とし、doc のコメントと意図を説明するコメントは日本語で書きます。
 
 ## 並列と再現性
diff --git a/docs/design/04-foundation/02-schedules.md b/docs/design/04-foundation/02-schedules.md
index 1663357724aece8fe39713077960dd91789da241..75bfde54f0620b537967f01e2701f8d094ca9e68 100644
--- a/docs/design/04-foundation/02-schedules.md
+++ b/docs/design/04-foundation/02-schedules.md
@@ -10,7 +10,7 @@
 |---|---|---|---|
 | 計算期間 | 利息を数える期間と日数係数 | `ScheduleSpec` | 半年ごとの固定脚の期間 |
 | 支払い | 一つ以上の計算期間の額を払う日 | `PaymentRule` | 期間末の 2 営業日後、四半期の期間を年に一度まとめて払う |
-| リセット | 期間の値を一度に決める fixing の日 | `ResetRule` | TIBOR の期首の 2 営業日前、MtM XCCY の元本の FX のリセット |
+| リセット | 期間の値を一度に決める fixing の日 | `ResetRule` | TIBOR の期首の 2 営業日前、MtM XCCY の元本の FX のリセット（ペアのスポットラグだけ期首から遡った日） |
 | 観測 | 期間の中で集約する観測の列と重み | `ObservationRule` | TONA の日次の複利、原油の月中の平均 |
 
 どの指標を観測するかは商品が決め、日程は日付と重みだけを決めます。
@@ -180,7 +180,12 @@ pub struct ResetRule {
     pub calendar: CalendarSpec,             // fixing のカレンダー
     pub bdc: BusinessDayConvention,
 }
-pub fn reset_dates(schedule: &Schedule, rule: &ResetRule, refdata: &RefData) -> Result<Vec<Date>>;
+// 期間ごとのリセット日と、カレンダーの外挿などの警告を返す。カレンダーは参照データの基本カレンダーの表から引く。
+pub fn reset_dates(
+    schedule: &Schedule,
+    rule: &ResetRule,
+    calendars: &BTreeMap<CalendarId, BaseCalendar>,
+) -> Result<(Vec<Date>, Vec<Warning>)>;
 ```
 
 | 慣行 | 規則 |
@@ -210,11 +215,11 @@ pub enum ObservationDays {
     Months,                                             // 期間に含まれる参照月（月次の指標）
     Endpoints(EndpointDates),                           // 公表された複利の指数の両端（Accrual か Shifted）
 }
-pub enum ObservationMethod {        // ISDA 2021 の方法に対応する。n は観測日の数（1 以上）
+pub enum ObservationMethod {        // ISDA 2021 の方法に対応する。n は観測日の数（型で 1 以上）
     Plain,
-    Lookback(u8),                   // 各日の値に n 観測日前の値を使う。重みは計算期間の日
-    ObservationShift(u8),           // 観測期間を n 観測日前へずらし、値も重みも観測期間で数える
-    Lockout(u8),                    // 最後の n 観測日は締めの日の値を使う。重みは計算期間の日。確認書の lockout の営業日数 k は Lockout(k − 1)
+    Lookback(NonZeroU8),            // 各日の値に n 観測日前の値を使う。重みは計算期間の日
+    ObservationShift(NonZeroU8),    // 観測期間を n 観測日前へずらし、値も重みも観測期間で数える
+    Lockout(NonZeroU8),             // 最後の n 観測日は締めの日の値を使う。重みは計算期間の日。確認書の lockout の営業日数 k は Lockout(k − 1)
 }
 pub enum Aggregation {
     Compound { form: CompoundingForm, spread: SpreadTreatment },   // SimpleFactor か PowerFactor、Exclusive か Inclusive
diff --git a/docs/design/04-foundation/03-money-and-quantities.md b/docs/design/04-foundation/03-money-and-quantities.md
index 00a57f4050ecb91bc9f291df4735cf8d8c4bcc75..3d427b594bdd635da82c05e61cc2fef2ec58039e 100644
--- a/docs/design/04-foundation/03-money-and-quantities.md
+++ b/docs/design/04-foundation/03-money-and-quantities.md
@@ -167,7 +167,7 @@ pip の指数は通貨ペアの慣行のデータです（USDJPY は 2、EURUSD
 | θ 勾配と quote 感応度の和 | `RiskVector<K>` のキー型（compile-fail） | `quote_delta_with_param_gradient` |
 | 尺度の要るフィールドに尺度の無い数を書く（`spread = 25`） | 境界の解析（`Err`） | `pricing_config_rejects_invalid_documents`、`standard_product_toml_rejects_ambiguous_dates_and_malformed_terms` |
 
-- compile-fail の行は trybuild のテストで、同じ名前の `_ok` の正しい使い方のテストと対にします。
+- compile-fail の行は testkit の `compile_fail_pair!` で、同じ名前の `_ok` の正しい使い方のテストと対にします。
   `_ok` も通すことで、誤用の側の失敗が無関係なコンパイルエラーでないことを確かめます（[品質](../11-quality.md)）。
 - 実行時のタグの行は、エラーの種類と内容を確かめる単体テストです。
 - 表は混同しやすく黙って誤るものに絞ります。
diff --git a/docs/design/04-foundation/04-reference-data.md b/docs/design/04-foundation/04-reference-data.md
index 966463a81a10978b1c3e60c5dc121471a24356d1..27338356f883b9754009c0dcdb157bf8859b3b21 100644
--- a/docs/design/04-foundation/04-reference-data.md
+++ b/docs/design/04-foundation/04-reference-data.md
@@ -108,7 +108,7 @@ pub enum IndexKind {
 | `OvernightRate` | 翌日物の率 | 指標の予想曲線。期間の複利は条件付きで DF の比に畳む（[スケジュール](02-schedules.md#df-の比に畳める条件)） |
 | `CompoundedRateIndex` | 翌日物の複利の指数の水準 | 元の翌日物の指標の予想曲線から水準の比を導く |
 | `SwapRate` | 慣行のスワップのパーレート | 慣行の変動の指標の予想曲線と、慣行の担保条件の割引曲線 |
-| `FxFixing` | FX のスポットのレート | 担保条件ごとの FX フォワード |
+| `FxFixing` | FX のスポットのレート。受渡日は fixing 日からペアの慣行のスポットラグだけ進めた日（[FX の慣行](#fx-の慣行)） | 担保条件ごとの FX フォワード |
 | `AssetPrice` | 資産の単価 | 資産のフォワード（配当、借株、保管費用の曲線） |
 | `EquityIndexLevel` | 指数の水準 | 指数単位の証券のフォワードをポイントの価値で割る |
 | `FuturesSettlement` | 先物の清算価格（単価） | 最終の清算日を含むどの観測日も現在の先物価格（コモディティは限月の先物の曲線の値、金融先物は原資産のフォワード） |
@@ -159,16 +159,13 @@ kind = { overnight_rate = { ccy = "JPY", day_count = "ACT365F" } }
 pub enum Convention {
     Deposit(DepositConvention),
     Fra(FraConvention),
-    Ois(OisConvention),
-    Irs(IrsConvention),             // 固定対タームレート
-    BasisSwap(BasisSwapConvention),
-    XccySwap(XccyConvention),
-    FxForward(FxForwardConvention),
-    RateFuture(RateFutureConvention),
-    CapFloor(CapFloorConvention),
-    Swaption(SwaptionConvention),
-    Cds(CdsConvention),
-    // 資産クラスを足すときに variant を足す
+    Ois(Box<OisConvention>),
+    Irs(Box<SwapConvention>),       // 固定対タームレート
+    BasisSwap(Box<SwapConvention>), // 同一通貨の変動対変動
+    XccySwap(Box<SwapConvention>),  // 通貨ごとの脚、MtM なら脚の fx_reset
+    FxForward(Box<FxForwardConvention>),
+    RateFuture(Box<RateFutureConvention>),
+    // M2 以降: CapFloor、Swaption、Cds など、資産クラスを足すときに variant を足す
 }
 pub struct OisConvention {
     pub spot_lag: DateOffset,
@@ -176,16 +173,49 @@ pub struct OisConvention {
     pub float: OvernightLegConvention,      // 指標、観測の規則
     pub payment: PaymentRule,
 }
+pub struct SwapConvention {
+    pub spot_lag: DateOffset,
+    pub legs: [SwapLegConvention; 2],       // 脚ごとの通貨、期間、支払い、クーポン（fixed、term、compounded）、FX リセット
+}
 ```
 
+IRS、ベーシススワップ、XCCY は脚の並びだけが違うので、一つの `SwapConvention` を共有し、variant が脚の組の検査（IRS は固定とターム金利、ベーシスは同一通貨の二つの変動、XCCY は二通貨）を決めます。
+`refdata/conventions.toml` では `[[convention]]` の `id` と `kind` の表に書きます。
+
 ```toml
-[convention."JPY-TONA-OIS".ois]
-spot_lag = "+2BD"
+[[convention]]
+id = "JPY-TONA-OIS"
+
+[convention.kind.ois]
+spot_lag = "2BD"
 fixed = { frequency = "1Y", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
-float = { index = "JPY-TONA", observation = { days = { business_days = "JPTO" }, method = "plain", aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } }
-payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" }
+payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
+
+[convention.kind.ois.float]
+index = "JPY-TONA"
+observation = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } }
+
+[[convention]]
+id = "JPY-TIBOR6M-IRS"
+
+[convention.kind.irs]
+spot_lag = "2BD"
+
+[[convention.kind.irs.legs]]
+currency = "JPY"
+periods = { frequency = "6M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
+payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
+coupon = "fixed"
+
+[[convention.kind.irs.legs]]
+currency = "JPY"
+periods = { frequency = "6M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
+payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
+coupon = { term = { index = "JPY-TIBOR-6M" } }
 ```
 
+全種類の例は [M1 のサンプル](../../../examples/m1-rates/refdata/conventions.toml) にあります。
+
 - 取引は慣行の ID を参照し、必要なフィールドだけを取引の条件で上書きします（上書きが優先します）。
 - 慣行の中身を変えると、それを参照する取引の解釈が変わります。
   既存の取引の解釈を変えずに慣行を改めるときは、新しい ID の慣行を足します。
@@ -193,19 +223,39 @@ payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "fo
 
 ### FX の慣行
 
+通貨ペアの現物の受渡慣行は、ペアごとに一つの `FxConvention` だけが持ちます。
+評価時点の換算レート X(0) の短期の連鎖（[FX のフォワード](../06-market/02-curves.md#fx-のフォワード)）、FX fixing の受渡日、MtM XCCY の元本リセットの観測日、FX フォワードの校正商品の日程は、すべてこの慣行から決めます。
+
 ```rust
 pub struct FxConvention {
-    pub pair: CcyPair,                      // 市場の標準の向き（EURUSD、USDJPY）
-    pub spot: SpotRule,                     // 営業日の数、数えるカレンダー、有効な日のカレンダー
-    pub pip_exponent: u8,                   // USDJPY は 2、EURUSD は 4
-    pub smile: Vec<(Tenor, FxSmileConvention)>,  // 満期の区分ごとの ATM、デルタ、バタフライの定義
-    pub triangulation: FxTriangulation,     // Direct か Via(USD)
+    pub pair: CcyPair,              // 市場の標準の向き（EURUSD、USDJPY）。pip とポイントはこの向きで読む
+    pub spot_lag: FxSpotLag,        // OneBusinessDay（USDCAD）か TwoBusinessDays（USDJPY）
+    pub settlement: LocalEventTime, // 各受渡日の受渡の瞬間。過ぎると当日の受渡はできない（同日受渡の締め）
+    pub pip_exponent: u8,           // USDJPY は 2、EURUSD は 4
 }
 ```
 
-- スポット日は `count_on` のカレンダーで営業日を数え、`valid_on` の営業日に調整した日です（USDJPY は T+2）。
-- スマイルの慣行の意味（USDJPY はプレミアムを含めたデルタの delta-neutral の ATM で、1Y 以下はスポットのデルタ、1Y を超えるとフォワードのデルタ）と、quote を校正の目標にする規則は [ボラティリティ](../06-market/04-volatility.md) が定めます。
-- quote の無い通貨ペアの X(0) は、`Via` の通貨を通る X(0) どうしの積で作ります。
+- 営業日は両通貨の決済カレンダーの和で数えます（`RefData::fx_calendar_spec`）。
+  スポット日は約定日から `spot_lag` の営業日だけ進めた日、FX fixing の受渡日は fixing 日から同じだけ進めた日、MtM XCCY のリセットの観測日は期首から同じだけ遡った日です（明示日程の `reset` があればその日）。
+- `RefData::fx_convention` は USDJPY と JPYUSD のどちらで引いても同じ慣行を返します。
+  同じペアの慣行を両方の向きで持つこと、表のキーと `pair` の向きの食い違い、未登録の受渡の地域は、参照データを作るときの `InvalidData` です。
+- FX fixing の指標（`IndexKind::FxFixing`）は、そのペアの慣行を要します。
+  慣行が無ければ、指標の定義の参照の欠けとして `InvalidData` です。
+- 市場の定義の `[[fx]]` は quote と担保の区分だけを書き、スポットラグ、締め、pip の桁を持ちません。
+  宣言された短期の quote が慣行のスポットラグで連鎖するか（T+1 は ON だけ、T+2 は TN と ON）と、quote の向きが慣行の向きかは、校正とワークスペースの読み込みで照合します（`FxDefinition::validate_references`）。
+- FX フォワードの慣行（`Convention::FxForward`）は、向きと満期の営業日規則だけを持ちます。
+- スマイルの慣行（USDJPY はプレミアムを含めたデルタの delta-neutral の ATM で、1Y 以下はスポットのデルタ、1Y を超えるとフォワードのデルタ）は、FX のボラ面を足す M2（ATM）と M7（スマイル）でこの慣行に足します（[ボラティリティ](../06-market/04-volatility.md)）。
+- quote の無い通貨ペアの X(0) は、市場の定義の `via` の通貨を通る X(0) どうしの積で作ります（三角換算の経路は市場の定義が宣言します）。
+
+ワークスペースでは `refdata/conventions.toml` に書きます。
+
+```toml
+[[fx]]
+pair = "USDJPY"
+spot_lag = "two_business_days"
+settlement = { time = "15:00:00", zone = "Asia/Tokyo" }
+pip_exponent = 2
+```
 
 ### CDS の指数
 
@@ -347,21 +397,52 @@ impl Error { pub fn context(self, what: impl Display) -> Error; }   // 「取引
 | `Execution` | 並列の実行器が取り消しや期限で止まった、ワーカーのスレッドや領域を確保できなかった | 呼び出し側（実行器の設定と実行環境） |
 | `Internal` | scioxide が `InvalidInput` を返した、並列のコールバックが panic した | rhodium の開発者 |
 
+内部の不変条件（同じ crate が作った値の前提）の扱いは、次の順で決めます。
+
+1. 型で表します。
+   呼び出し側が種類を知っているなら、enum で返して照合し直さず、種類ごとの関数か型を返します。
+   検証済みの前提は、検証した型に持たせます。
+2. 型で表せない前提（実行器が区間の結果を返すという約束、グラフのトポロジカル順、lock の poison など）だけは、破れたときに `Error::Internal` を返します。
+   直前の行に `// 不変条件: <保証の出どころ>` を書き、どこがその前提を保証するかを示します。
+   注記の無い `Internal` は `.githooks` の禁止パターン `P-internal` が拒否します（[品質の CI](../11-quality.md#ci)）。
+3. 同じモジュールが作った表への添字は許し、表を作った場所を `// 不変条件:` で示します。
+
+- 不変条件の破れを `debug_assert!` だけで検出することはしません。
+  release のビルドで誤った値のまま計算が続くためです。
+- `Internal` を制御の流れ（`.ok()?` で捨てる失敗など）に使いません。
+
 - 警告はエラーと別の値 `Warning` で返し、計算を止めません。
 
 ```rust
 pub enum Warning {
+    // 参照データと fixing
     CalendarExtrapolated { calendar: CalendarId, date: Date },
     ExplicitDateIsHoliday { date: Date },
     UnpublishedFixingForecast { index: IndexId, point: ObservationPoint },
     FixingEstimateUsed { index: IndexId, point: ObservationPoint },
     FallbackApplied { index: IndexId, replacement: IndexId },
     NonRepresentativeIndex { index: IndexId },
-    CalibrationFitOutsideTolerance { block: BlockId, quote: QuoteId, residual: f64 },
+    // リニア評価の名前付きの近似（価格設定の conventions で禁止できる）
+    PaymentDelayIgnored,
+    CrossCollateralProjection,
+    FraDiscounting,
+    MtmFxForwardProduct,
+    DailyForwardProduct,
+    RoundingAtForward,
+    // quote と校正
     QuoteUnitSuspected { quote: QuoteId },
+    FxCrossInconsistent { pair: CcyPair },
+    CalibrationFitOutsideTolerance { quote: QuoteId },
+    CalibrationGaussNewton { quote: QuoteId },
 }
 ```
 
+- 近似の 6 つの variant の意味と使う場面は [リニアと解析の「近似の規約」](../07-pricing/02-linear-and-analytic.md#近似の規約)が定めます。
+  取引の準備（`Prepared::warnings`）と校正商品の準備（`CalibratedMarket::warnings`）が記録します。
+- `CalibrationFitOutsideTolerance` は最小二乗のグループで、当てはめの残差が binding の許容幅 σ を超えた quote ごとに出します。
+  `CalibrationGaussNewton` は、解の記録を Gauss–Newton で作ったグループの quote ごとに出し、その quote の感応度が近似であることを示します（[校正](../06-market/03-calibration.md)）。
+- 個々のフィールドの正本は rustdoc です。
+
 - ポートフォリオの計算は取引ごと（XVA はネッティングセットごと）の `Result` を返し、一件の失敗で全体を止めません。
   集計は除いた取引の一覧を必ず伴い、成功した値だけを黙って足しません。
   校正グループの失敗は、そのグループと下流のグループを読む取引だけを失敗にします。
@@ -404,3 +485,19 @@ DSL の誤りの表示（位置と提案）は [DSL](../05-products/02-dsl.md) 
 - 判断: すべての crate が `rhodium_core::Error` を使い、variant は誰が直すかで分けます。
 - 理由: エラーの型の変換が要らず、利用者は variant を見るだけで対処が分かります。
 - 採らない案: crate ごとのエラー型や番号付きのエラーコードの表は、変換と番号の管理が要るため採りません。
+
+### D5 FX の慣行は参照データに一つだけ置く
+
+- 判断: 通貨ペアのスポットラグ、同日受渡の締め、pip の桁を、参照データの `FxConvention` だけに置き、X(0)、FX fixing、MtM XCCY のリセット、FX フォワードの校正商品がそれを読みます。
+- 理由: 商品は市場の定義を読めないので、慣行を市場の定義に置くと商品の側に別の値（固定の T+2 など）が要り、T+1 のペアで受渡日を黙って誤ります。
+  正本が一つなら、同じペアの日付が評価の経路によって食い違いません。
+- 採らない案: 市場の定義の `[[fx]]` に受渡の慣行を書く案は、商品が読めず二重の正本になるため採りません。
+  FX fixing の指標ごとにスポットラグを書く案も、同じペアの慣行を指標の数だけ重ねるため採りません。
+
+### D6 内部の不変条件は型で表し、表せないものだけ注記付きの `Internal` にする
+
+- 判断: 同じ crate が保証した前提は型に持たせ、型で表せない前提だけを、`// 不変条件:` の注記を付けた `Error::Internal` で報告します。
+- 理由: 検証済みの値を実行時に照合し直すと、ありえない状態の分岐と `Internal` が増え、variant を足したときに黙って既定の分岐に落ちます。
+  型で表せない前提も、release で止まり、どこが保証するかを読めるようにします。
+- 採らない案: `debug_assert!` だけで検出する案は、release で誤った値のまま計算が続くため採りません。
+  添字の panic に任せる案は、利用者の入力で panic しない規則と区別できないため、同じモジュールが作った表に限ります。
diff --git a/docs/design/05-products/01-ir.md b/docs/design/05-products/01-ir.md
index 95ee7542deefe2549bc76a652206763e578e3615..ce8f71104562d006b9c883abbdad2a20280c1c16 100644
--- a/docs/design/05-products/01-ir.md
+++ b/docs/design/05-products/01-ir.md
@@ -191,7 +191,8 @@ V = Σ_{Pay(a, x, t)} E^Q[ D(t) · x ]   （a ごとに、通貨なら割引、
 
 `Linear` は、`Choose` を含まず、条件が確定済みか信用の恒等式で扱える事由で、すべての額がリニア評価の知っている項（定数、`Fixing`、FX の `Fixing` の逆数、`Compounded`、`Average`、生存で重み付けた額、回収率の項）の線形結合であることを意味します。
 項の照合の規則は [リニアと解析](../07-pricing/02-linear-and-analytic.md) が定めます。
-M0 の `linearity` は、簡約済みの通貨受渡に対して、定数と翌日物の複利を線形結合し、契約の丸めを適用できるかを判定します。
+`linearity` は、簡約済みの通貨受渡に対して、受渡額が定数、ターム金利と翌日物の複利、FX の fixing（とその逆数）、金利先物の清算価格の線形結合で、契約の丸めを適用できるかを判定します。
+予想どうしの積は FX と金利の一つずつ（MtM XCCY の元本リセット、規約 `MtmFxForwardProduct`）だけ、予想での割り算は FX の逆数と FRA の決済（N·(1 − (1 + τK)·P_p(e)/P_p(s))）だけを許します。
 曲線の解決や近似の可否は、価格計算の準備で確かめます。
 
 ## 評価時点での簡約
diff --git a/docs/design/05-products/03-standard-products.md b/docs/design/05-products/03-standard-products.md
index 05e81e2113dbd334cfda3b8f2d498bf32ac403b1..bd3b65651083598622e567dd8f2f2a8be203eb6a 100644
--- a/docs/design/05-products/03-standard-products.md
+++ b/docs/design/05-products/03-standard-products.md
@@ -53,27 +53,27 @@ impl StandardProduct {
 pub enum PayReceive { Pay, Receive }
 pub enum Position { Long, Short }   // Short は Long の IR の give
 
-pub struct Leg {
-    pub side: PayReceive,
-    pub notional: Notional,
-    pub schedule: ScheduleSpec,          // core::schedule
-    pub day_count: DayCount,
-    pub coupon: Coupon,
-    pub exchange: NotionalExchange,
+pub struct Leg {                         // フィールドは非公開。Leg::new か serde の境界で検証して作る
+    side: PayReceive,
+    notional: Notional,
+    schedule: ScheduleSpec,              // core::schedule
+    day_count: DayCount,
+    coupon: Coupon,
+    exchange: NotionalExchange,
+    payment_rounding: PaymentRounding,   // 受渡額の丸め（None か通貨の規則）
 }
 
 pub enum Notional {
     Fixed(Money),
-    Amortizing(Vec<(Date, Money)>),      // 各日以降の元本
     FxReset { base: Money, fixing: IndexId }, // MtM の XCCY: 各期の元本 = 基準元本を期首の FX fixing で脚の通貨に換えた額
+    // 未実装: Amortizing(Vec<(Date, Money)>)（各日以降の元本）
 }
 
 pub enum Coupon {
     Fixed { rate: Rate<Simple> },
     Term { index: IndexId, spread: Spread, gearing: f64 },
     Compounded { index: IndexId, spread: Spread, rule: ObservationRule },
-    Averaged { index: IndexId, spread: Spread, rule: ObservationRule },
-    Cms { index: IndexId, spread: Spread, gearing: f64 },
+    // 未実装: Averaged { index, spread, rule }（RFR 平均）、Cms { index, spread, gearing }（M2）
 }
 
 pub enum NotionalExchange { None, Final, InitialAndFinal }
@@ -87,9 +87,11 @@ pub enum NotionalExchange { None, Final, InitialAndFinal }
 | RFR 平均の脚 | `N·τ_i·(Average + s)` | builder は平均の重みを全期間で正規化して IR に渡す |
 | CMS の脚 | `N·τ_i·(g·Fixing(SwapRate, f_i) + s)` | 評価は凸性の補正が要る（解析） |
 | 元本交換 | 期首に −N、期末に +N（`NotionalExchange` による） | XCCY |
-| FX リセット | 各期の元本を、基準元本の通貨から脚の通貨へ fixing のペアの向きに従って換えた額（基準元本の通貨がペアの base なら掛け、quote なら割る）にし、元本の差額を期首に交換する。fixing は fixing 指標の暦で期首の 2 営業日前に観測し（スポット受渡が期首になる日）、明示日程の `reset` があればその日を使う（ターム金利の fixing 日と共通） | MtM の XCCY。USDJPY で JPY の元本を固定し USD の脚をリセットするなら N_USD = N_JPY / Fixing(USDJPY, 期首の 2 営業日前) |
+| FX リセット | 各期の元本を、基準元本の通貨から脚の通貨へ fixing のペアの向きに従って換えた額（基準元本の通貨がペアの base なら掛け、quote なら割る）にし、元本の差額を期首に交換する。fixing はスポット受渡が期首になる日（ペアの慣行のスポットラグを両通貨の休日の和で期首から遡った日。[FX の慣行](../04-foundation/04-reference-data.md#fx-の慣行)）に観測し、明示日程の `reset` があればその日を使う（ターム金利の fixing 日と共通）。新しい元本は期首の脚の通貨の決済時刻に交換し、fixing の予想はペアのスポットの受渡で読む | MtM の XCCY。USDJPY（T+2）で JPY の元本を固定し USD の脚をリセットするなら N_USD = N_JPY / Fixing(USDJPY, 期首の 2 営業日前) |
 
 `gearing` は無次元の係数です。
+CMS の脚は CMS スワップとともに M2 で足します。
+償却元本と RFR 平均の脚は M1 の範囲（最初の検証データセットの商品）に無く、それを要する商品を足すときに variant と表の行を実装します（IR は既存の `Pay`、`Average` で表せます）。
 スプレッドが 0 のときは `+ 0` の項を作らず、リニア評価が標準の形として認識できる IR にします。
 
 ### M0 の預金とスワップ
@@ -168,7 +170,7 @@ Pay(c, N·(P_k − P_{k−1}), s_k)     k = 0, 1, …, K
   値洗いの口座の測度 Q_M と三つの成立条件の下で、リニア評価では最初の未公表の値洗いだけが残り、PV = N·P_d(s_next)·(F(A) − P_last) です（[金利先物](../07-pricing/02-linear-and-analytic.md#金利先物)）。
   P_last は A までに公表された最新の清算価格で、無ければ取引価格です。
 - 解析プライサーは同じ式を使い、金利先物の F(A) を凸性の補正で求めます（[金利先物](../07-pricing/02-linear-and-analytic.md#金利先物)）。
-  モデルでは清算価格が E^{Q_M}_{d_k}[最終決済値] なので、MC で IR を評価しても同じ値になります（テスト `futures_margin_ir_vs_analytic`）。
+  モデルでは清算価格が E^{Q_M}_{d_k}[最終決済値] なので、MC で IR を評価しても同じ値になります（M1 はリニア評価との一致 `futures_margin_ir_vs_analytic`、M4 は MC との一致 `futures_margin_mc_vs_analytic`）。
 
 取引の `settlement_index` は限月固有の日次指標を指定します。`contracts` は正の整数で、`trade_price` は内部の価格単位（例として 96 表示なら `0.96`）です。`day_count` は原指標の率の日数基準に合わせます。約定日の清算後に建てた取引は `first_clearing` で次の清算日を指定できます。値洗いの暦と公表規則は指標、乗数と受渡時刻は金利先物の慣行から取ります。
 
diff --git a/docs/design/05-products/04-trades.md b/docs/design/05-products/04-trades.md
index 52fc1ba484467b639f63e2bcbaf4555d3d60bf6a..0ff0a07ed2642f4fdec5ba3bc170c63d414b343f 100644
--- a/docs/design/05-products/04-trades.md
+++ b/docs/design/05-products/04-trades.md
@@ -76,7 +76,8 @@ impl Trade {
   `to_product` は、`by_index` の鍵がターム金利でない場合、`bilateral` の指標が無い場合、`Custom` の代替が参照データと合わない場合に `Error::InvalidData` を返します。
 - Rust の builder で組み立てた `Contract` は `ProductSpec::Contract` で渡せます。
   `Trade::to_product(&RefData)` が返す `Product` は取引の行使記録を持ち、`prepare` や `reduce` には空の追加記録を渡せます。
-- `TradeProductResolver` は一つの `RefData` で多くの取引を `Product` にします。
+- `TradeProductResolver` は一つの `RefData` で多くの取引を `Product` に下ろし、取引と組にした `ResolvedTrade` を返します。
+  `ResolvedTrade` はフィールドを公開せず、`TradeProductResolver::resolve` と `ResolvedTrade::new` だけが作るので、取引と別の取引の IR を組めません（`PrepareInput::Trade` の入力。[prepare と evaluate](../07-pricing/01-pricing.md#prepare-と-evaluate)）。
   同じ指標・観測点・重み・対象期間を持つ RFR の観測の列と、その公表時刻と観測期間の解決を取引の間で共有します。
   検証と失敗は `Trade::to_product` と同じで、複数のワーカーから同時に呼べます。
   ワークスペースの読み込みはこれを使います（[入出力と実行のメモリ](../10-io-and-execution.md#メモリ)）。
diff --git a/docs/design/06-market/01-quotes.md b/docs/design/06-market/01-quotes.md
index 2a1cfd4bb810f375d12a3d2d58a5f5f0171d1717..8d4ce2b4c2c0596919b90f32f27477ad5f7981f8 100644
--- a/docs/design/06-market/01-quotes.md
+++ b/docs/design/06-market/01-quotes.md
@@ -86,7 +86,7 @@ impl QuoteSet {
 ```
 
 - 金利や価格は decimal で持ちます。bp、%、32 分の 1 のような尺度は、読み込みの境界で一度だけ decimal に換えます（[金額と量](../04-foundation/03-money-and-quantities.md)）。
-- FX フォワードポイントは pip 数で持ちます。校正はペアの受渡慣行から pip 幅を読み、レートの差へ換算します。
+- FX フォワードポイントは pip 数で持ちます。校正は参照データのペアの受渡慣行から pip 幅を読み、レートの差へ換算します。
 - bid と ask だけが与えられた quote は、境界で mid = (bid + ask)/2 にします。
   bid > ask は入力のエラーです。
 - quote が欠けていても `QuoteSet` は作れます。
@@ -230,7 +230,7 @@ pub enum QuoteKind {
 - 判断: `QuoteSet` の値は decimal で、bp や % は読み込みのときに換えます。
   FX フォワードポイントだけは pip 数のまま持ち、校正が binding の通貨ペアの pip 指数で一度だけレートの差に換えます。
 - 理由: 内部で尺度を持ち回ると、尺度の取り違えが計算のどこでも起きえます。
-  pip の桁は通貨ペアの慣行のデータで、quote の意味と同じく市場の定義が与えます（D1）。
+  pip の桁は通貨ペアの慣行のデータで、参照データの `FxConvention` が与えます（[FX の慣行](../04-foundation/04-reference-data.md#fx-の慣行)）。
   市場の定義を知らない読み込みの境界では換えられないため、換算は FX の定義を解決する校正の一か所に置きます。
 - 採らない案: quote ごとに尺度を付けて持つ案は、すべての読み出しで尺度を確かめる必要があるため採りません。
 
diff --git a/docs/design/06-market/02-curves.md b/docs/design/06-market/02-curves.md
index 04e2ea230d01e272c43f39d8ae419e960651bd45..8d8384d9a3b5daa28df46bad127f10c1e17f8855 100644
--- a/docs/design/06-market/02-curves.md
+++ b/docs/design/06-market/02-curves.md
@@ -188,9 +188,11 @@ MaxForward:  f̄_{k,i} = −(L_k(t_{i+1}) − L_k(t_i)) / (t_{i+1} − t_i),   L
   付利のフロアを `MaxForward` で表すのは本質的価値だけの決定論的な近似で、フロアのオプション価値は含みません。
   付利を割引に入れる扱いでは、割引が数えた本質的価値とパス上の付利の差（時間価値）を ColVA が数えます（[XVA の指標](../09-xva/02-metrics.md#指標の式)）。
 - `Steps` の発効日は参照データの会合日の表が決めます。
-  会合日の無い定型のテナーの OIS だけで階段を当てる場合は、罰則付きの最小二乗で解きます（[校正](03-calibration.md#残差の形)）。
+  会合日の無い定型のテナーの OIS だけで階段を当てる場合は、罰則付きの最小二乗で解きます（`StepPenalty`。[校正](03-calibration.md#残差の形)）。
 - 部品の式は構造を作るときに展開し、`Ref` の循環と役割の矛盾を拒否します。
-- M1 の Rust API では `CurveSpec::from_expr` に曲線の式と名前付き部品を渡します。
+- M1 の Rust API では `CurveSpec::from_expr` に曲線の式と名前付き部品を渡し、ワークスペースの `market/definition.toml` では同じ式を serde の表現のまま `expr` に書きます（[校正](03-calibration.md#ワークスペースの市場定義)）。
+  M1 で実装した部品の形は、`Nodes { part, pillars, extrapolation }`（補間は ln DF の線形、`FlatForwards`）、`Stitch { short, long, at }`（`blend` なし）、`Steps { part, meetings }`（発効の瞬間の列）、`Jumps { part, dates }`、`Scale { factor, body }`（有限な定数の倍率）、`FxImplied { part, pair, base, pillars, extrapolation }` です。
+  表の `interp`、`blend`、θ の要素の倍率の `Scale`、会合表と年末の規則からの日付の解決は、それを使う商品や曲線とともに足します（非線形の補間とプロキシの β は M2 以降、ファイルの日付の規則は M6-4）。
   `Scale` の倍率には有限な定数を指定します。
   会合の発効時刻、年末の跳びの時刻、接続点、ピラーを、参照データに従ってあらかじめ解決した UTC の瞬間で指定します。
   接続は境界で値を連続にし、長期側の傾きへ切り替えます。
@@ -210,9 +212,10 @@ MaxForward:  f̄_{k,i} = −(L_k(t_{i+1}) − L_k(t_i)) / (t_{i+1} − t_i),   L
 
 - 原点は値 0 の固定の節点です。
   `ZeroRate` は原点と最初のピラー t₁ の間を z(t) = z(t₁) とし、L(t) = −z(t₁)·t で原点で 0 にします。
-- M0 の補間には `FlatForwards` を用います。
+- M1 までの補間は `FlatForwards` だけで、表の他の方式は非線形の補間を要する曲線とともに M2 以降で足します。
 - 最後のピラーより後は `FlatInstantaneousForward`（既定）、`FlatAverageForward`、`Reject` のどれかで外挿します。
 - ピラーは `FromBindings`（そのピラーを受け持つ校正商品がこの曲線を読む最後の日付）か明示の日付・テナーで置きます。
+  M1 で実装したのは明示の UTC の瞬間だけで、`FromBindings` とテナーからの解決はファイル形式の見直し（M6-4）で足します（外貨担保の曲線の短期のピラーは FX の受渡の区間から決めます）。
   支払いの遅れのある OIS では、予想の曲線のピラーは観測期間の終わり、割引の曲線のピラーは支払日です。
 - 非局所的な補間は、ある満期の quote のリスクを遠い満期へ漏らすので、リスクの報告に使う曲線には局所的な方式を勧めます。
 - 評価の経路の補間は、AD の型を通すために rhodium が `R: Real` のジェネリックなコードとして持ちます。
@@ -249,6 +252,7 @@ MaxForward:  f̄_{k,i} = −(L_k(t_{i+1}) − L_k(t_i)) / (t_{i+1} − t_i),   L
 - 付利のスプレッドとフロア（割引か ColVA）と資金調達のコスト（割引か FCA/FBA）は、価格設定がどちらか一方でだけ数えます。
 - 指標の予想の曲線は、担保条件によらず指標ごとに一本です（担保ごとの予想は、観測できる校正商品が無く校正できないため）。
   別の担保条件の取引に同じ予想を使うのは、担保の違いによる予想の調整を無視する決定論的なベーシスの仮定です。
+  予想の担保は曲線を当てる校正商品の担保で、校正がその割引曲線を予想の測度として市場構造に記録します（[担保の違いの判定](../07-pricing/02-linear-and-analytic.md#担保の違いの判定)）。
 
 ## FX のフォワード
 
@@ -264,19 +268,20 @@ n_c(t) = ln( F_c(t) / X(0) ) = L_b|c(t) − L_q|c(t)
 ```
 
 - S はスポット受渡日で、pts は価格の単位に直したポイントです。
-- 短期の quote（ON、TN、あれば SN）は、ペアの慣行（両通貨の休日を合わせたカレンダー）で実際の開始と終了の受渡日に解決してから連鎖させます。
+- 短期の quote（ON、TN、あれば SN）は、参照データのペアの慣行（[FX の慣行](../04-foundation/04-reference-data.md#fx-の慣行)。スポットラグ、締め、pip の桁、両通貨の休日を合わせたカレンダー）で実際の開始と終了の受渡日に解決してから連鎖させます。
   T+2 のペア（USDJPY）ではスポット日から TN、ON の順にたどり、T+1 のペア（USDCAD）ではスポット日が翌営業日なので ON だけをたどります。
   たどるのは、評価時点 A 以後に受渡できる最も早い受渡日 d_e までです。
   当日の受渡の締めを過ぎていれば当日は受渡できないので、ON は使わず、d_e は翌営業日です。
 - 受渡日の瞬間 t_e と t_S は [日付を瞬間にする規則](../04-foundation/01-dates-and-calendars.md#日付を瞬間にする規則) の受け渡しの行で決め、評価時点 A まで市場の担保の曲線で割り戻します。
   当日の受渡が午後なら、午前の評価時点からその瞬間までの日中の割引も入ります（テスト `fx_x0_short_dated_chain`）。
 - スポット X(S) と短期のポイントは direct mark の θ で、X(0) は導いた値です。
-- FX の Rust API は `MarketDefinition::with_fx` で受渡慣行と quote を受け取ります。
+- FX の Rust API は `MarketDefinition::with_fx` で quote と担保の区分を受け取り、受渡慣行は校正のときに参照データから読みます。
+  短期の quote が慣行のスポットラグで連鎖するかと quote の向きは、校正とワークスペースの読み込みで照合します。
   ON と TN はスポットに向かう区間の差、SN とそれ以後のポイントはスポットを始点とする区間の差です。
   受渡時刻ちょうどは当日の受渡を含めます。
   スポットから評価時点までの曲線が必要で、使用する短期区間が無ければスポットの受渡を起点にします。
   宣言された有効な区間の quote が欠けた場合は、該当する校正グループが失敗します。
-  FX ポイントと `FxImplied` の曲線の組（ポイントがあれば曲線が要り、曲線には締め後も残る TN、SN か日付付きの quote が要る）は、締めの前後で成否が変わらないよう、宣言された quote で市場の定義を作るときに確かめます。
+  FX ポイントと `FxImplied` の曲線の組（ポイントがあれば曲線が要り、曲線には締め後も残る TN、SN、日付付きの quote か長期の XCCY が要る）は、締めの前後で成否が変わらないよう、宣言された quote で市場の定義を作るときに確かめます。
 - X(0) は通貨ペアごとに一つで、担保によりません。
   評価時点ちょうどの受渡には、担保が効く期間が無いためです。
   通貨スコープの換算、三角換算、報告通貨への換算は、すべてこの X(0) を使います。
@@ -298,6 +303,11 @@ n_c(t) = ln( F_c(t) / X(0) ) = L_b|c(t) − L_q|c(t)
   直接の quote と三角換算の値の食い違いは、quote の検査の警告として報告し、換算には直接の quote を使います（[quote](01-quotes.md)、テスト `fx_cross_direct_vs_triangulated`）。
 - 外貨担保の割引曲線（例: `USD/csa:JPY-TONA`）は `FxImplied` で、n(t) の節点を θ にし、L_USD|JPY = L_JPY + n（実装は向きを構造を作るときに一度だけ決めます）。
   短期は FX スワップ、長期は MtM XCCY に校正し、接続点をまたぐ商品があれば両方を一つのグループで解きます。
+  Rust API とファイルでは、`FxImpliedCurve` の `long_end` にピラーごとの XCCY（quote、`XccySwap` の慣行、テナー、担保）を書きます。
+  慣行の 2 通貨はペアの 2 通貨と一致しなければならず、参照データに照らす検査（`FxDefinition::validate_references`）で確かめます（テスト `fx_implied_long_end_rejects_an_xccy_convention_of_another_pair`）。
+  長期の行は締めによらずピラーを一つずつ決めるので、T+1 のペアの spot、ON と長期の XCCY だけの定義も、締めの前後で校正できます（テスト `fx_implied_long_end_alone_survives_the_cutoff_of_a_t_plus_one_pair`）。
+  ピラーは FX ポイントの受渡のピラーの後に昇順で続き、XCCY の quote は外貨担保曲線の `basis` ブロックを当てる `ParRate(XccyBasis)` の binding になります（ブロックの要素の番号は締めの前後で ON が抜けるとずれるので、一般の binding から要素を指しません）。
+  FX ポイントと XCCY の行は一つの正方系のグループになり、TONA、SOFR と FX のマークが上流です（テスト `fx_implied_long_end_reproduces_fx_points_and_xccy_quotes`）。
   FX フォワード、XCCY ベーシス、両通貨の割引の無裁定の関係は構造の上でいつも成り立ちます。
 - 担保の無い取引の FX フォワードは、価格設定の `fx_collateral` が指定する担保で作ります（既定はペアの市場の担保 m）。
   無担保の曲線を両通貨に使うと、市場の FX フォワードを再現しないためです。
diff --git a/docs/design/06-market/03-calibration.md b/docs/design/06-market/03-calibration.md
index db1483063ac598d5bc5b4f0a4af06a2b76b05e32..d1be303bd584bb228ccae584e51c160b5a10ddf7 100644
--- a/docs/design/06-market/03-calibration.md
+++ b/docs/design/06-market/03-calibration.md
@@ -42,31 +42,102 @@ pub enum CalibrationTarget {
 ```
 
 - 評価日は `valuation_zone` で評価時点を現地の日付にした日で、校正商品の取引日はそれを `valuation_calendar` の営業日に調整した日です（[日付を瞬間にする規則](../04-foundation/01-dates-and-calendars.md#日付を瞬間にする規則)）。
-- 目標とテンプレートの組が許されるか（`ImpliedVol` はオプションのテンプレートだけ、`ParRate` は quote の値を一つ持つテンプレートだけ）と、テンプレートの PV が quote の値について affine か（`Shape`）は、`MarketDefinition` を作るときに一度だけ確かめて記録します。
+- 目標とテンプレートの組が許されるか（`ImpliedVol` はオプションのテンプレートだけ、`ParRate` は quote の値を一つ持つテンプレートだけ）は、`MarketDefinition` を作るときに一度だけ確かめます。
+  テンプレートと慣行の組（観測金利の丸めや quote のスプレッドの扱い）は参照データを要するので、校正の入口で一度だけ確かめます（`InstrumentTemplate::validate`）。
 - FX スワップや MtM XCCY のように二つの通貨を受け渡す校正商品は、支払い通貨ごとに `collateral` の下の割引曲線を価格設定の `discounting` の表で解決します（[価格計算](../07-pricing/01-pricing.md#要求の解決)）。
 - `InstrumentTemplate` は、慣行（`ConventionId`）と座標（テナー、満期、ストライク、限月）から標準商品を作る閉じた enum です（[標準商品](../05-products/03-standard-products.md)）。
   テンプレートの引数は慣行、座標、リテラルだけで、quote の値は `ParRate` の評価のときにだけ入ります。
 - パーは [パーレートの解](../07-pricing/02-linear-and-analytic.md#パーレートの解) で求めます。
   affine な商品（固定金利、`Exclusive` のスプレッド）は K = 0 と K = 1 の二回の評価からパーを正確に得て、この式は `R: Real` のまま評価できます。
-  affine でない商品（ZC インフレスワップの固定脚 N·((1 + K)^n − 1)、`Inclusive` のスプレッド）は f64 のニュートン法で解き、θ についての微分は `implicit_solution_into` の陰関数定理で 2 次以上まで戻します。
+  M1 の校正テンプレート（預金、FRA、金利先物、OIS、IRS、ベーシス、FX フォワード、MtM XCCY）はすべて affine で、quote のスプレッドは `Exclusive` に限ります（[標準商品](../05-products/03-standard-products.md#校正商品のテンプレート)）。
+  affine でない商品（ZC インフレスワップの固定脚 N·((1 + K)^n − 1)、`Inclusive` のスプレッド）は f64 のニュートン法で解き、θ についての微分は `implicit_solution_into` の陰関数定理で 2 次以上まで戻します（`par::solve_nonlinear`、テスト `par_nonlinear_lsq_quote_delta`）。
+  テンプレートの PV が quote について affine かの区別（`Shape`）と、それを校正のパーに使う経路は、最初の affine でないテンプレート（M9 の ZC インフレスワップ）とともに足します。
 - ボラの逆算は f64 で解き、θ についての微分は陰関数定理で戻します（[アーキテクチャ](../03-architecture.md#数値型と自動微分)）。
 
-### M0 の市場定義
+### ワークスペースの市場定義
 
 ワークスペースの `market/definition.toml` では、評価地域と営業日カレンダーを指定します。
-`[curves.<id>]` には `currency`、`roles`、UTC の瞬間を昇順に並べた `pillars` と、任意の `extrapolation` を指定します。
+`[curves.<id>]` には `currency`、`roles` と、曲線の形を次のどちらかで書きます。
 
-`[[instruments]]` の各行は一つの OIS パーレートを結び付けます。
-`quote` はスナップショットの ID、`kind` は `"rate"`、`target` は `"par_rate"`、`template` は `{ ois = { convention = "JPY-TONA-OIS", tenor = "1M" } }` の形です。
-`collateral` は価格設定と共通の担保条件、`determines` は `{ curve = "JPY-TONA" }` または `{ curve = "JPY-TONA", element = 0 }` を取り、`element` は 0 始まりのピラー番号です。
-`transform` の既定は `"identity"` です。
-JPY と USD の完全な入力は[サンプルワークスペース](../../../examples/m0-ois/market/definition.toml)を参照してください。
+- `pillars`（UTC の瞬間の昇順）と任意の `extrapolation`: 一つの `nodes` ブロックの節点の曲線（M0 の形）。
+- `expr`: 部品の式 `CurveExpr` の serde の表現そのもの（[曲線と market](02-curves.md#log-curve-の代数)）。
+  `nodes`（`part`、`pillars`、`extrapolation`）、`ref`、`forward_ref`、`sum`、`scale`（`factor`、`body`）、`flat`、`max_forward`（`parts`、`grid`）、`stitch`（`short`、`long`、`at`）、`steps`（`part`、`meetings`）、`jumps`（`part`、`dates`）を入れ子に書けます。
+  瞬間は UTC の日時で明示し、`flat` の率は `"25bp"` のように接尾辞付きの文字列で書きます（裸の数は単位の書き忘れとして拒否します）。
+  外貨担保曲線（`FxImplied`）は `[[fx]]` の `direct.implied` が作るので、`[curves]` には書けません。
 
-M1 の `instruments` には FRA、IRS、ベーシス、MtM XCCY、金利先物、FX フォワードのテンプレートも書けます。金利先物の行は `target = "price"` と `kind = "futures_price"`、アウトライトの FX は `target = "price"` と `kind = "price"` を使います。FX のスポットとポイントは `[[fx]]` で定義できます。
+`[[instruments]]` の各行は一つの quote を一つの目標に結び付けます。
+
+- `quote` はスナップショットの ID、`kind` は quote の種類、`target` は `"par_rate"`、`"price"`、`"direct_mark"` のどれかです。
+- `par_rate` と `price` は `template`（例: `{ ois = { convention = "JPY-TONA-OIS", tenor = "1M" } }`）と、価格設定と共通の担保条件 `collateral` を取ります。
+  `direct_mark` はどちらも取りません。
+- `determines` は `{ curve = "JPY-TONA" }`（`nodes` ブロック）、`{ curve = "JPY-TONA", part = "long" }`（名前の付いた部品のブロック）、`element`（0 始まりの要素番号）を組み合わせます。
+- `transform` の既定は `"identity"` で、年末の跳びは `{ turn_spread_to_jump = { value = 0.0109589, basis = "ACT365F" } }`（期末の受渡期間の日数係数を明示）と書きます。
+- `tolerance`（`"0.5bp"` などの quote の単位）を書いた行だけのグループは最小二乗で解きます（[残差の形](#残差の形)）。
+  会合日の階段の罰則は `[[step_penalties]]`（`curve`、`part`、`weight`、`scale`）、Gauss–Newton の記録は `[[gauss_newton]]`（`curve`、`part`）に書きます。
+
+JPY と USD の完全な入力は[M0 のサンプル](../../../examples/m0-ois/market/definition.toml)と[M1 のサンプル](../../../examples/m1-rates/market/definition.toml)を参照してください。
+
+M1 の `instruments` には FRA、IRS、ベーシス、MtM XCCY、金利先物、FX フォワードのテンプレートも書けます。金利先物の行は `target = "price"` と `kind = "futures_price"`、アウトライトの FX は `target = "price"` と `kind = "price"` を使います。FX のスポットとポイントは `[[fx]]` で定義でき、外貨担保曲線の長期を決める MtM XCCY は `[[fx]]` の `direct.implied.long_end` に書きます（[FX のフォワード](02-curves.md#fx-のフォワード)）。
 
 `[lgm.USD]` の `mean_reversion` と `volatility` は、区間順に並べた `mark` quote の ID です。それぞれの個数は `breaks` に並べた UTC の区切りの数に 1 を足した値です。区切りがなければ、一組のマークを全期間で使います。
 
-### 曲線の部品とテナーの展開を使う例
+### ワークスペースの金利ボラ
+
+`[[volatilities]]` に面の ID、指標、原資産の IRS または OIS の慣行と担保を指定します。
+ボラの種類は `normal`、`lognormal`、`shifted_lognormal` です。シフト付き対数正規では面の `shift` を指定します。
+個々の入力 ID と座標は面の `quotes` に置き、値は既存の `quotes/snapshot.csv` に置きます。
+ストライクの数値と当てはめの許容幅は十進表記です。以下の例は `examples/m1-rates` の参照データと曲線を使います。
+
+```toml
+[[volatilities]]
+id = "JPY-SABR"
+type = "swaption_cube"
+index = "JPY-TONA"
+convention = "JPY-TONA-OIS"
+collateral = "uncollateralized"
+vol_kind = "normal"
+expiry_time = { time = "11:00", zone = "Asia/Tokyo" }
+sabr = { beta = 0.5, shift = 0.01, fit = "exact_atm" }
+quotes = [
+    { quote = "JPY.SABR.1Y.2Y.ATM", expiry = "1Y", tenor = "2Y", strike = "atm" },
+    { quote = "JPY.SABR.1Y.2Y.LOW", expiry = "1Y", tenor = "2Y", strike = { atm_offset = -0.0025 }, atm_quote = "JPY.SABR.1Y.2Y.ATM", tolerance = 0.0001 },
+    { quote = "JPY.SABR.1Y.2Y.HIGH", expiry = "1Y", tenor = "2Y", strike = { atm_offset = 0.0025 }, atm_quote = "JPY.SABR.1Y.2Y.ATM", tolerance = 0.0001 },
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
+    { quote = "JPY.CAP.1Y.ATM", maturity = "1Y", strike = "atm" },
+    { quote = "JPY.CAP.2Y.ATM", maturity = "2Y", strike = "atm" },
+    { quote = "JPY.CAP.1Y.K1", maturity = "1Y", strike = { absolute = 0.01 } },
+    { quote = "JPY.CAP.2Y.K1", maturity = "2Y", strike = { absolute = 0.01 } },
+]
+```
+
+- 満期テナーは市場の約定日から進め、評価カレンダーの Modified Following で営業日調整します。
+  スワップのテナー軸は、同じ約定日から調整済み終点までの ACT/365F です。原資産の開始日には満期日の約定に慣行のスポットラグを適用します。
+- `fit = "least_squares"` では各 quote に、`fit = "exact_atm"` ではウィングに `tolerance` を書きます。
+  ATM だけなら `fit = { atm_only = { rho = -0.15, nu = 0.4 } }` と書きます。
+  `atm_quote` を指定した行では二つの quote を加えて絶対ボラにし、指定した両方の quote から校正を更新できます。
+  `sabr.shift` は正規 SABR 用です。補間の既定は `total_variance` で、`sabr.interpolation = "parameters"` も選べます。
+- 後決め RFR には `vol_kind = "normal"` と `rfr_variance = "lyashenko_mercurio"` を指定します。
+  ターム金利のキャップには IRS の慣行を使い、タームの fixing 時刻でボラの満期を決めます。
+
+たとえば ATM ボラに `60bp`、キャップの絶対ストライクのボラに `60bp`、SABR の LOW・HIGH にそれぞれ `2bp` と `1bp` を入力できます。
+`rh calibrate <workspace> --date 2026-11-24` は、面 ID・quote ID 順に
+`surface,quote_id,vol_kind,target_vol,model_vol,residual` を標準出力に書きます。
+`--as-of` で瞬間を指定することもできます。
+差分の行でも目標と再現値は絶対ボラです。
+失敗したグループがあれば、成功した行を出力し、診断を標準エラーに書いて終了コード 1 を返します。
+市場を組み立てられない入力エラーの終了コードは 2 です。
+
+### 曲線の部品を使う例
 
 ```toml
 version = 1
@@ -74,39 +145,51 @@ valuation_calendar = "JPTO"
 valuation_zone = "Asia/Tokyo"
 
 [curves.JPY-TONA]
-roles = ["discount", "projection"]
 currency = "JPY"
-expr.stitch = { at = "1Y", short = { steps = { meetings = "BOJ", count = 4 } }, long = { nodes = { interp = "monotone_convex" } } }
-jumps = { dates = "turn_of_year", part = "turns" }
+roles = ["discount", "projection"]
+
+[[curves.JPY-TONA.expr.sum]]
+stitch.at = 2027-12-03T00:00:00Z
+stitch.short.steps = { part = "meetings", meetings = [2026-12-19T00:00:00Z, 2027-01-23T00:00:00Z, 2027-03-19T00:00:00Z] }
+stitch.long.nodes = { part = "long", pillars = [2028-12-04T00:00:00Z, 2031-12-03T00:00:00Z, 2036-12-03T00:00:00Z] }
+
+[[curves.JPY-TONA.expr.sum]]
+jumps = { part = "turns", dates = [2027-01-04T00:00:00Z] }
 
 [curves.JPY-TIBOR-6M]
-roles = ["projection"]
 currency = "JPY"
-expr.sum = [{ ref = "JPY-TONA" }, { nodes = { interp = "flat_forwards", part = "basis" } }]
+roles = ["projection"]
+expr.sum = [{ ref = "JPY-TONA" }, { nodes = { part = "basis", pillars = [2028-12-04T00:00:00Z, 2031-12-03T00:00:00Z] } }]
 
 [[instruments]]
-quote = "JPY.TONA.OIS.{tenor}"
-tenors = ["1Y", "2Y", "3Y", "5Y", "7Y", "10Y", "15Y", "20Y", "30Y"]
-template = { ois = { convention = "JPY-TONA-OIS" } }
-target = "par_rate"
+quote = "JPY.TONA.OIS.10Y"
 kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "JPY-TONA-OIS", tenor = "10Y" } }
 collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
 determines = { curve = "JPY-TONA", part = "long" }
 
 [[instruments]]
-quote = "JPY.TIBOR6M-TONA.BASIS.{tenor}"
-tenors = ["2Y", "5Y", "10Y", "20Y", "30Y"]
-template = { basis_swap = { convention = "JPY-TIBOR6M-TONA-BASIS" } }
-target = "par_rate"
+quote = "JPY.TONA.TURN"
 kind = "spread"
+target = "direct_mark"
+transform = { turn_spread_to_jump = { value = 0.0109589, basis = "ACT365F" } }
+determines = { curve = "JPY-TONA", part = "turns", element = 0 }
+
+[[instruments]]
+quote = "JPY.TIBOR6M-TONA.BASIS.5Y"
+kind = "spread"
+target = "par_rate"
+template = { basis_swap = { convention = "JPY-TIBOR6M-TONA-BASIS", tenor = "5Y" } }
 collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
 determines = { curve = "JPY-TIBOR-6M", part = "basis" }
 ```
 
 - ファイルでは `target` と `template` を別のキーで書き、読み込みのときに一つの `CalibrationTarget` に組みます。
-
-- `quote` の `{tenor}` は `tenors` の値から ID を作る雛形で、ID の文字列を解析して座標を読むのではありません。
 - 年末の跳びのマーク、会合日の短期の階段を決める短期の OIS なども、同じ `instruments` の行で書きます。
+- 行ごとの誤り（式の形、目標とテンプレートの組、未知の部品）は、ファイルと行（`curves[ID]`、`instruments[i]`）を付けた `InvalidData` です（テスト `market_file_reports_invalid_curve_expressions_with_positions`）。
+- 会合表（`meetings = "BOJ"`）、年末などの日付の規則（`turn_of_year`）、接続点のテナー（`at = "1Y"`）を参照データから解決する書き方、`tenors` と `{tenor}` の雛形による行の展開、非線形の補間（`interp`）は、ファイル形式の見直し（M6-4）で足します。
+  M1 は解決済みの UTC の瞬間を書きます。
 
 ## 校正グループ
 
@@ -118,8 +201,12 @@ determines = { curve = "JPY-TIBOR-6M", part = "basis" }
 4. θ の各要素を決めるグループを一つに定めます。
    グループ内の複数の残差が同じ要素を対象にできます。
    未割当の要素、同じ要素への二重の直接指定、直接指定と当てはめの競合は定義を解決するときに拒否します。
+   ただし許容幅を持つ直接マークどうしは同じ要素に重ねられ、最小二乗で一つの値に決まります（下の検算）。
 
-- 正方系で解くグループでは、未知数と残差の数が一致しなければなりません。
+- グループの形は binding の許容幅で決まります（`InstrumentBinding::with_tolerance`）。
+  許容幅を持たない binding だけのグループは正方系で、未知数と残差の数が一致しなければなりません。
+  すべての binding が許容幅を持つグループは最小二乗で、残差と罰則の行の数が未知数以上でなければなりません。
+  許容幅を持つ binding と持たない binding が一つのグループに混ざる定義は `InvalidInput` です。
 - 曲線、ボラ、モデルのパラメータ、ハザードは同じグラフに入ります。
   quote から XVA のモデルまでが一つの DAG なので、XVA の感応度も quote まで戻せます（[XVA のリスク](../09-xva/03-xva-risk.md)）。
 - 同じ深さの独立なグループは並列に解けます。
@@ -179,12 +266,17 @@ F_i(θ, q) = m_i(θ) − q̃_i(q),    ∂F/∂q = −D,   D = ∂q̃/∂q（疎
 
 - 残差ごとに、読み出す全 quote と当てはめの対象を記録します。
   例えば ATM と差分の二つを使う残差は、どちらの変更でも再校正の対象です。
-- 最小二乗の重みは W = diag(1/σ_i²)（σ_i は binding の許容幅）で、θ と q に依存しません。
+- 最小二乗の重みは W = diag(1/σ̃_i²) で、θ と q に依存しません。
+  σ_i は binding の quote の単位の許容幅で、σ̃_i = σ_i·D_i が残差の単位の幅です（D_i は変換 q̃ の係数）。
+  最小二乗の残差は quote の単位で測り、q̃ を q の定数倍に保ちます。
+  FX アウトライトの正方系は対数の残差 ln F(θ) − ln q で解きますが、最小二乗では価格の差 F(θ) − q = q·expm1(ln F − ln q) を使います（D_i = 1）。
+  対数のまま 1/σ で重みを付けると、許容幅の単位が q 倍ずれます。q = 150、F = 150.10、σ = 0.01 の外れ（10σ）は 0.067σ と測られて警告されず、他の行との相対の重みも q² だけ狂います（テスト `fx_outright_least_squares_measures_residuals_in_price_units`）。
 - 最小二乗の ∂G/∂θ は J_Fᵀ W J_F + Σ_i (WF)_i ∇²F_i + ∇²P で、既定は 2 階の項を含む厳密な形です（このときだけ 1 次の quote 感応度が厳密です）。
-  要素の多い当てはめでは Gauss–Newton（残差の 2 階項を省き、二次罰則の ∇²P は残す）を選べ、選んだことを `GroupSolution` に記録します。
-- 会合日の階段を定型のテナーの OIS だけで当てる場合は、罰則 P(θ) = (λ/2) Σ_i ((f_{i+1} − f_i)/s)² を加えた最小二乗にします（s は参照の幅、λ > 0）。
-  罰則だけで決まる水準も曲線に使いますが、市場が織り込む政策金利として報告しません。
-- 検算: min (θ − q₁)² + 4(θ − q₂)² の解は θ = (q₁ + 4q₂)/5 で、quote への勾配は (1/5, 4/5) です。
+  要素の多い当てはめでは Gauss–Newton（残差の 2 階項を省き、二次罰則の ∇²P は残す）を選べ（`MarketDefinition::with_gauss_newton`）、選んだことを `GroupSolution::hessian` と警告 `Warning::CalibrationGaussNewton` に記録します。
+- 会合日の階段を定型のテナーの OIS だけで当てる場合は、罰則 P(θ) = (λ/2) Σ_i ((f_{i+1} − f_i)/s)² を加えた最小二乗にします（s は参照の幅、λ > 0。`StepPenalty` と `MarketDefinition::with_step_penalties`）。
+  罰則は残差の行 √λ·(f_{i+1} − f_i)/s として目的 ½‖r̃‖² に加えます。
+  罰則だけで決まる水準も曲線に使いますが、市場が織り込む政策金利として報告しません（テスト `penalized_meeting_steps_fit_standard_tenor_ois`）。
+- 検算: min (θ − q₁)² + 4(θ − q₂)² の解は θ = (q₁ + 4q₂)/5 で、quote への勾配は (1/5, 4/5) です（σ₁ = 1、σ₂ = 1/2 の二つの直接マーク。テスト `weighted_direct_marks_solve_the_hand_least_squares`、`lsq_adjoint_oracle`）。
 
 ## 解法
 
@@ -197,9 +289,17 @@ F_i(θ, q) = m_i(θ) − q̃_i(q),    ∂F/∂q = −D,   D = ∂q̃/∂q（疎
 | `Direct` | direct mark | 変換そのもの |
 | `Explicit` | 上流だけで陽に決まるブロック | 陽な写像を `R` で評価する |
 
-- 解法の希望は市場の定義に書け、書かなければグループの形（正方か、局所的か、陽か）から選びます。
+- M1 は、正方系を `Newton`（`TrustRegionSystem`）で、最小二乗を `LeastSquares`（scioxide-optimize の Levenberg–Marquardt）で、すべての残差が直接マークの正方系を `Direct` で解きます。
+  形はグループの binding の許容幅から決まり（[残差の形](#残差の形)）、解法の希望を市場の定義に書く仕組み、`Bootstrap`、`IterativeBootstrap`、`Explicit` は、局所的な補間や陽な写像のブロックを足すとき（M2 以降）に足します。
 - 残差のヤコビアン J_F は AD で求めてソルバーに渡します。
   LM は最後のヤコビアンを返さないので、解の記録は解いた点で線形化し直して作ります。
+- 最小二乗の停留は、記録の ∂G/∂θ で求めた補正 δθ = (∂G/∂θ)⁻¹·G が各行を自然な単位（正方系の解法と同じ変換後の残差、罰則はフォワードの差）で動かす量 ‖Ĵ δθ‖ を校正の許容誤差と比べて判定します（正方系では残差のノルムに一致します）。
+  FX アウトライトの行は価格の単位の動きを 1/q 倍して対数の単位に戻します。150 前後の価格の丸め（1 ulp が約 3e-14）は、絶対値の許容誤差より大きいからです。
+  LM は残差の残る問題で一次収束にとどまるので、LM の後に厳密なヘシアンのニュートン補正で停留条件を磨きます。
+  ∂G/∂θ が正定値でない（極小でない、識別できない向きがある）解は `Failed` です。
+  Gauss–Newton を指定したグループの記録 J̃ᵀJ̃ は常に半正定値で、識別できない向きしか検出できません。
+  そこで磨きと記録は J̃ᵀJ̃ のまま、停留を確かめた点で厳密なヘシアンを一度だけ求めて正定値を判定し、解を採用する条件を近似した記録から分けます（反復ごとの Hessian–vector product は省いたままです）。
+- 直接マークの入力検査（有限で、FX スポットの種類なら正）は、正方系か最小二乗かを分ける前に行います（テスト `least_squares_direct_marks_check_quotes_before_fitting`）。
 - 制約は解く座標の変換で表します（ハザード h = exp(u)、SABR の α = exp(a)、ρ = tanh(b)、ν = exp(c)）。
   θ には自然なパラメータを入れ、記録は θ の座標の ∂G/∂θ で残すので、陰関数定理の結果は座標の取り方によりません。
   解いた後に、θ の座標での停留（尺度を調整した ‖G‖ が許容値の内）と、変換の飽和（|ρ| が 1 に張り付く、α が 0 に張り付く）を検査します。
@@ -208,11 +308,13 @@ F_i(θ, q) = m_i(θ) − q̃_i(q),    ∂F/∂q = −D,   D = ∂q̃/∂q（疎
 ### 初期値と warm start
 
 - 解が一意なグループ（単調な bootstrap、強凸な最小二乗）は、前の `CalibratedMarket` の θ を初期値にします（warm start）。
+  M1 の強凸な最小二乗は、すべての行が直接マークのグループです。
   評価時点が変わったとき（同じ日の中を含む）は、ピラーの日付が同じブロックは前の値を写し、日付が変わったブロックは前の曲線を新しいピラーで評価した値を初期値にします。
+  この評価時点をまたぐ warm start（`calibrate` の `warm` 引数）は、時間を進める theta と P&L 説明（M3）で足し、M1 は同じ評価時点の `update` だけで前の θ を使います。
   warm start から収束しなければ、quote からの簡易な推定値で一度だけ解き直します。
 - 解の一意性が保証されないグループ（SABR、非線形の補間を含む Newton）は、固定の初期値の集合から解き、目的関数の値が最小の候補を選びます（warm start は使いません）。
 - warm start の有無で結果が違うのはソルバーの許容誤差の範囲だけです（[品質](../11-quality.md#許容誤差)）。
-  監査の再計算のために、warm start を使わない指定も選べます。
+  監査の再計算には、`update` ではなく `calibrate` で全体を解き直せば、warm start を使わない結果が得られます。
 - 現行の `update` は、各残差が順に異なるピラーを決め、そのピラーについて厳密に単調な OIS グループで前回の θ を使います。
   同一時刻に一回支払う両脚では割引因子が消え、複数回支払う場合は同一曲線で割り引いた変動脚が端点の差に縮約でき、始点が既知で quote が非負なら単調性を確認できます。
   循環グループや、この条件を満たさないグループには固定の初期値を使います。
@@ -232,19 +334,20 @@ pub struct CalibrationGroup {
 pub enum GroupOutcome { Solved(Box<GroupSolution>), Failed(Error), Unavailable }
 
 /// 解いたグループの、随伴に使う線形化。
-/// 行は row_quotes、列は CalibrationGroup の blocks、quotes、upstream の順。
+/// 行は正方系では row_quotes、最小二乗では blocks の θ。列は CalibrationGroup の blocks、quotes、upstream の順。
 pub struct GroupSolution {
     dg_dtheta: LuDecomposition<f64>,  // ∂G_k/∂θ_k の分解（scioxide-linalg）
     dg_dq: CsrMatrix<f64>,            // ∂G_k/∂q_k（scioxide-sparse）
     dg_dup: CsrMatrix<f64>,           // ∂G_k/∂θ_up
-    hessian: HessianMode,             // 最小二乗で厳密か Gauss–Newton か
+    hessian: Option<HessianMode>,     // 最小二乗で厳密か Gauss–Newton か。正方系は None
     residual_norm: f64,
     iterations: usize,
+    residuals: Vec<f64>,              // 解での残差 F（重みを掛けない）。行は row_quotes
+    theta_error: Vec<f64>,            // 停止の誤差 δθ = (∂G/∂θ)⁻¹·G
 }
 
-/// 市場の定義と quote から校正する。warm は前の校正（初期値と、一意なグループの warm start に使う）。
-pub fn calibrate(def: &MarketDefinition, quotes: &QuoteSet, ctx: &Context<'_>,
-                 warm: Option<&CalibratedMarket>) -> Result<CalibratedMarket>;
+/// 市場の定義と quote から、固定の初期値で校正する。
+pub fn calibrate(def: &MarketDefinition, quotes: &QuoteSet, ctx: &Context<'_>) -> Result<CalibratedMarket>;
 
 pub struct CalibratedMarket {
     market: Market<f64>,
@@ -258,9 +361,14 @@ pub struct CalibratedMarket {
 - 依存と quote の対応は `GroupSolution` ではなく `CalibrationGroup` が持ちます。
   `Failed` と `Unavailable` のグループにも依存と quote の対応を残す必要があるからです（[失敗](#失敗)）。
   `GroupSolution` は `Solved` のグループだけが持つ線形化です。
-- 現行の `calibrate` は評価時点と市場定義から構造を組み立て、固定の初期値から解きます。
+- `calibrate` は評価時点と市場定義から構造を組み立て、固定の初期値から解きます。
   バンプで照合するときは `calibrate_with_residual_tolerance` で残差の絶対許容誤差を指定できます。
-  最小二乗のグループを足すまで、`GroupSolution` は `hessian` を持ちません。
+- 最小二乗の G は停留条件 J̃ᵀ r̃（r̃ は重みを掛けた残差と罰則の行）で、行は自グループの θ です。
+  厳密な ∂G/∂θ と ∂G/∂θ_up は、目的 Φ = ½‖r̃‖² の Hessian–vector product を自グループの n 方向だけ求めて作ります（∇²Φ の上の n 行が ∂G/∂θ、下の行が対称性から ∂G/∂θ_up）。
+  ∂G/∂q = −J̃ᵀ W D は J̃ が q によらないので 1 階で厳密です。
+  形が違っても記録の型は同じなので、随伴（[感応度](../08-risk/01-sensitivities.md#3-校正の随伴)）は正方系と同じ式で流れます。
+- 最小二乗のグループで残差 |F_i| が σ̃_i を超えた quote は警告 `Warning::CalibrationFitOutsideTolerance` にします。
+  警告は解から求めるので、`update` で再利用したグループでも毎回求め直します。
 
 陰関数定理により、解の感応度は次のとおりです。
 
@@ -284,7 +392,7 @@ impl CalibratedMarket {
 ```
 
 - 更新には同じ評価時点、参照データ、fixing、価格設定と校正条件を使います。
-  これらが変わった場合は、`calibrate` で市場を校正し直します（前の校正を `warm` に渡せます）。
+  これらが変わった場合は、`calibrate` で市場を校正し直します（前の校正を初期値に渡す `warm` 引数は M3 で足します。[初期値と warm start](#初期値と-warm-start)）。
 - 変わった quote（`QuoteSet::with_overrides` が返す集合か、前の `QuoteSet` との差）を読むグループと、DAG でその下流のグループだけを解き直します。
   他のグループの θ の値を引き継ぎ、`GroupSolution` を共有します。
   準備済みの校正商品と依存順も共有し、quote の値に由来する警告を新しいスナップショットで計算します。
@@ -335,7 +443,7 @@ impl CalibratedMarket {
 | 確かめること | 対象のグループ | 条件 | 満たさないとき |
 |---|---|---|---|
 | quote の再現 | 完全に決まるグループ（正方系、bootstrap、direct mark、陽な写像） | すべての残差がソルバーの許容誤差の内 | 実装の誤り |
-| 停留と当てはめ | 最小二乗のグループ | 尺度を調整した ‖G‖ が許容値の内（停留）。残差 F_i を binding の許容幅 σ_i と比べて記録する（当てはめ） | 停留しなければ `Failed`。残差が σ_i を超えた quote は警告 `Warning::CalibrationFitOutsideTolerance` |
+| 停留と当てはめ | 最小二乗のグループ | 補正 δθ = (∂G/∂θ)⁻¹·G が各行を動かす量 ‖Ĵ δθ‖ が許容値の内（停留）。残差 F_i を binding の許容幅 σ̃_i と比べて記録する（当てはめ） | 停留しなければ `Failed`。残差が σ̃_i を超えた quote は警告 `Warning::CalibrationFitOutsideTolerance` |
 | 校正と取引の一致 | すべてのグループ | 同じ θ で校正商品を取引として評価した値が、残差のモデルの量 m_i(θ) と一致する | 校正と取引の評価の経路が分かれた誤り |
 
 - 最小二乗のグループは、正しい解でも quote の残差が残るので、quote の再現を求めません。
diff --git a/docs/design/06-market/04-volatility.md b/docs/design/06-market/04-volatility.md
index 5d10c824736cea55de6e301c3696b3a19914da34..b320eefd0a349c3d057764415b53585d26d0f66f 100644
--- a/docs/design/06-market/04-volatility.md
+++ b/docs/design/06-market/04-volatility.md
@@ -27,7 +27,7 @@ pub enum Expiry { Tenor(Tenor), Date(Date), Contract(ContractMonth) }
 ```
 
 - 境界のボラの値（quote、面の節点、結果）の型は `NormalVol`、`LognormalVol`、`ShiftedLognormalVol`（[金額と量](../04-foundation/03-money-and-quantities.md)）で、種類の違うボラを比べたり足したりする誤りは型で起きません。
-  評価の内側の読み出しは `R` を返し、種類の対応は下の [読み出し](#読み出し) のとおり `prepare` が一度だけ確かめます。
+  評価の内側の読み出しは `R` で値を計算し、種類の対応は下の [読み出し](#読み出し) のとおり `prepare` が一度だけ確かめます。
   シフトの違うシフト付き対数正規ボラは比べられません。
 - 面の軸は、満期 × ストライク（キャップレット、FX、株式、コモディティ）か、満期 × テナー × ストライク（スワップション）です。
 - RR と BF は二つのデルタのストライクの組の量で、単独の座標を持ちません（下の FX のスマイル）。
@@ -64,6 +64,8 @@ pub enum SabrFit { AtmOnly, LeastSquares, ExactAtm }
 
 - キューブのノードは（満期, テナー）の組で、各ノードが自分の SABR のブロック（α、ρ、ν）を持つ独立した校正グループです。
   ノードは互いに独立なので並列に解けます。
+  校正グループの θ には、同じ曲線で求めたノードのフォワードも派生値として記録します。
+  面の読み出しが節点の原資産を評価し直さず、曲線の感応度はこの派生値を通して戻ります。
 - 各ノードの残差は、そのノードのフォワードスワップレートと年金現価を読むので、割引と予想の曲線に依存します（曲線のグループの下流です）。
 - `LeastSquares` は各ストライクのボラを目標にした 3 変数の LM です。
   ストライクの ATM オフセットとボラの quote の表現は独立で、絶対ボラなら `Identity`、ATM ボラとの差なら `AtmPlusOffset` を使います。
@@ -79,7 +81,12 @@ pub enum SabrFit { AtmOnly, LeastSquares, ExactAtm }
 1. 絶対のストライクの列ごとに、満期の短い順に、区間ごとに一定のキャップレットのボラを、そのストライクのキャップの価格（フラットボラで評価した価格）を再現するように逐次に解きます。
 2. ATM の列は満期ごとにストライクが違うので、手順 1 の列をストライクの方向に補間した面の上で、区間の ATM の水準だけを未知数にして解きます。
    ATM のストライクは、構造の曲線から作るモデルのフォワード F(θ) です。
-3. 解いたボラを満期 × ストライクの節点として θ に入れます。
+3. 解いたボラと ATM のストライクを θ に入れます。絶対ストライク列は列ごとの満期区間で一定とし、
+   同じキャップレット満期のストライク方向はボラを線形補間します。格子の外は端の値を使います。
+
+金利キャップの ATM は、対象クーポンの割引済み年金現価で重み付けしたフォワードです。
+キャップごとのクーポン列は、長いキャップの先頭と一致する必要があります。
+同じ ATM ストライクに絶対ストライクの節点があると、ATM の水準を独立に決められません。
 
 - キャップの最初の期を含めるかは指標の種類が決めます。
   ターム金利（TIBOR）は最初の期の fixing が取引時に確定するので除き、後決めの RFR（TONA、SOFR）は観測期間の間は不確実なので含めます。
@@ -139,23 +146,25 @@ FX の ATM、デルタ、バタフライの慣行は、通貨ペアの慣行（
 
 ```rust
 impl<'v, R: Real> MarketView<'v, R> {
-    /// 面の満期と（テナーの）断面に結んだ点で、ストライク strike（面のストライクの座標の絶対値）のボラを読む。失敗しない。
-    pub fn vol_at(&self, at: VolAt<'v>, strike: R) -> R;
+    /// 面と結んだ点で、絶対ストライクのボラを読む。
+    /// ストライクが展開の定義域から外れるとエラーになる。
+    pub fn vol_at(&self, at: VolAt<'v, '_>, strike: R) -> Result<R>;
     /// ボラとストライク・満期の微分（局所ボラの Dupire の式に使う）。
-    pub fn smile_at(&self, at: VolAt<'v>, strike: R) -> SmileValue<R>;
+    pub fn smile_at(&self, at: VolAt<'v, '_>, strike: R) -> SmileValue<R>;
 }
 
 pub struct SmileValue<R> { vol: R, d_vol_dk: R, d2_vol_dk2: R, d_total_variance_dt: R }
 ```
 
 - `VolPoint` は満期と（スワップションの）テナーの断面と補間の重みを計算済みの点で、`prepare` が面ごとに作って `Prepared` に保存します。
-  `VolAt<'v>` は面と点を `bind_slots` が組にして検査した値で、読み出しはこの組だけを受け取ります（[曲線と market](02-curves.md#market-view)）。
+  `VolAt<'v, 'a>` は面のハンドルと点を `MarketView::bind_vol` が組にして検査した値で、読み出しはこの組だけを受け取ります（[曲線と market](02-curves.md#market-view)）。
+  点を `Slots` に持たせて `bind_slots` で他の点とまとめて結ぶ形は、ボラを読むプライサーを足すときに加えます。
   別の面の点で読む誤りは binding のエラーになり、組を作らずに読むコードはコンパイルできません。
 - ストライクは面に依存しない数なので読み出しの引数にし、フォワードに依存するストライクも評価のときに渡せます。
-- `VolRef<'v>` と `VolAt<'v>` はブランド `'v` の view の中でだけ使えます。
+- `VolRef<'v, 'a>` と `VolAt<'v, 'a>` はブランド `'v` の view の中でだけ使えます。
 - ボラの種類は面の定義の `VolKind` が持ち、プライサーは種類を見て Black、Bachelier、シフト付き Black を選びます。
   種類とプライサーの基本式の対応は、`prepare` が要求を面に解決するときに一度だけ確かめ、対応しない組（正規ボラしか使えない式に対数正規の面など）は `Error::InvalidInput` です（テスト `vol_kind_formula_mismatch`）。
-  評価の内側の読み出しは `R` を返し、種類を型では区別しません。
+  評価の内側の読み出しでは `R` で値を計算し、種類を型では区別しません。
 - ストライクがフォワードに依存する座標（ATM、ATM からの差、マネーネス）では、評価のときのフォワード（θ の関数）からストライクを作って読み出しに渡すので、曲線の動きがボラに入ります。
 
 ## 裁定の検査
diff --git a/docs/design/07-pricing/01-pricing.md b/docs/design/07-pricing/01-pricing.md
index 11552560292fdf6f8cc469ef7ecc6e7dc96d4a3a..6f12604eccc5b454d0f05b203fa1b2d1be4b6b6d 100644
--- a/docs/design/07-pricing/01-pricing.md
+++ b/docs/design/07-pricing/01-pricing.md
@@ -67,7 +67,7 @@ pub struct Context<'a> {
 | `volatility` | `BTreeMap<VolKey, VolId>` | ボラの要求（対象と商品の区分）からボラ面への表 | 空（ボラが必要な要求の規則は必須） |
 | `credit` | `BTreeMap<CreditKey, CurveId>` | 信用参照から生存曲線（と回収率）への表 | 空（信用が必要な要求の規則は必須） |
 | `assets` | `BTreeMap<Asset, CurveId>` | 資産からフォワードの曲線への表 | 空（資産のフォワードが必要な要求の規則は必須） |
-| `fx_collateral` | `FxCollateralRule` | 担保の無い FX フォワードを作る担保（評価時点の X(0) は担保によらない） | 市場の定義が通貨ペアごとに宣言した市場の担保 |
+| `fx_collateral` | `BTreeMap<CcyPair, CurveId>` | 担保の無い FX フォワードを作る FX の担保区分 ID（ペアはどちらの向きで書いてもよい。評価時点の X(0) は担保によらない） | 市場の定義が通貨ペアごとに宣言した市場の担保 |
 | `funding` | `FundingTreatment` | 資金調達の費用を割引に入れるか、FCA/FBA として XVA で数えるか（[XVA の指標](../09-xva/02-metrics.md)）。`InDiscounting` は `funding_curves` の借入と運用が同じ曲線の通貨だけに使える | `InXva` |
 | `funding_curves` | `BTreeMap<Currency, FundingCurves>` | 通貨ごとの自社の借入と運用の曲線。`InDiscounting` の資金調達の割引と、XVA の FCA と FBA が同じ値を読む | 空（`InDiscounting` の割引か FCA/FBA に使う通貨は必須） |
 | `remuneration` | `RemunerationTreatment` | 担保の付利を割引に入れるか、ColVA として数えるか（[XVA の指標](../09-xva/02-metrics.md)） | `InDiscounting` |
@@ -94,9 +94,9 @@ pub struct FundingCurves { pub borrow: CurveId, pub lend: CurveId }
   上から順に照合し、最初に一致した一件だけを適用します。
   複数の上書きを合成する規則は持ちません。
 
-### M0 で読み込める設定
+### M1 で読み込める設定
 
-M0 では以下の形式で割引、予想、評価時点の受渡、リニア方式、リニア評価の近似を指定します。
+M1 では以下の形式で割引、予想、評価時点の受渡、評価方式、リニア評価の近似を指定します（M0 の形に CTD と付利のフロアの担保条件、金利先物の方式、FX の担保区分を足した形）。
 無担保の割引に自社の資金調達費用を入れるときは、`funding = "in_discounting"` と `funding_curves` に対応する通貨の借入・運用の曲線を指定します。
 現金担保の照合には通貨、付利指標、スプレッドを使います。
 
@@ -118,6 +118,10 @@ curve = "JPY-TONA"
 product = "swap"
 method = "linear"
 
+[[methods]]
+product = "rate_future"
+method = { analytic = "rate_future" }
+
 [conventions]
 disabled = []
 max_payment_delay_business_days = 10
@@ -229,6 +233,10 @@ pub enum CollateralKey {
 }
 ```
 
+- M1 で実装した variant は `Discounting` と `Forecast` です。
+  FX のフォワードは要求を作らず、受け渡しの担保から FX の担保区分を直接解きます（担保付きは、その担保の下の両通貨の割引曲線の組と一致する区分。無担保は `fx_collateral`、無ければ市場の担保。一致する区分が無ければ市場の担保（区分が一つならその区分）で予想して `CrossCollateralProjection` を記録します）。
+  金利先物の解析プライサーは LGM のマーク（`[lgm.<通貨>]`）のブロックを直接読み、`Model` の要求は M4 のモデルとともに足します。
+  `Volatility` は M2、`Credit` は M5、`AssetForward` は M8 で足します。
 - 要求は `prepare` が一度だけ作ります。
   信用の要求は主体・優先順位・条項・通貨まで、担保の要求は付利条件や CTD の候補まで区別して解決します。
   割引は受け渡しごとの支払い通貨と、その受け渡しのマージンセットの担保から作るため、ノベーションなどで担保が受け渡しごとに違っても扱えます。
@@ -239,7 +247,8 @@ pub enum CollateralKey {
   ここで構造の `Arc` の同一性と、必要な全ブロックがこの市場で利用可能なことを検査します。
   成功した場合だけ、役割別のハンドルと、曲線・面と点を組にして検査した `CurveAt`、`ForwardAt`、`VolAt` を作り、評価ループに渡します（[曲線と market](../06-market/02-curves.md#market-view)）。
   評価ループの読み出しはこの組だけを受け取るため、曲線と点を組み替えられません。
-- 担保の CSA が違えば割引曲線は違いますが、予想曲線は指標ごとに一つで、担保の違いによる予想の調整はリニア評価の近似 `CrossCollateralProjection` として扱います（[リニアと解析](02-linear-and-analytic.md)）。
+- 担保の CSA が違えば割引曲線は違いますが、予想曲線は指標ごとに一つで、担保の違いによる予想の調整はリニア評価の近似 `CrossCollateralProjection` として扱います。
+  担保が違うかどうかは、予想曲線の予想の測度（校正商品の担保の割引曲線）と支払いの割引曲線で判定します（[担保の違いの判定](02-linear-and-analytic.md#担保の違いの判定)）。
 - 規則の無い要求は `Error::MissingData` です。
   黙って既定の曲線を使うことはしません。
 
@@ -265,6 +274,10 @@ pub enum Method {
 3. どちらも無く、簡約した IR のすべての項と条件が [リニア評価の恒等式](02-linear-and-analytic.md#恒等式) に該当するなら `Linear` を使います。
 4. それ以外はエラー（`Error::Unsupported`）にし、候補の方式を挙げます。
 
+- 商品の種類を持つのは型付きの標準商品の取引だけです。
+  builder の契約と、変更した取引（一部解約や置き換えの `ProductSpec::Amended`）は種類を持たず、`methods` の表に当たらないので、3 と 4 で決まります（例: 一部解約した金利先物は、変更後の IR をリニア評価します）。
+- 解析プライサーは型付きの標準商品を読み（[リニアと解析の D3](02-linear-and-analytic.md#d3-解析プライサーは型付きの標準商品を読む)）、取引の約定日と未決済の手数料も入力にします。
+  標準商品を持たない入力に解析方式を選ぶと `Error::Unsupported` です。
 - 選んだ方式が商品に適用できないとき（`Linear` に `Choose` を含む IR、`Analytic` に対応しない商品、`MonteCarlo` に `Choose` を含む IR、モデルが商品の観測を扱えない）は、`prepare` が `Error::Unsupported` を返します。
   別の方式へ黙って切り替えることはしません。
 - 方式の適用可能性はこの規則と各方式の前提だけで決め、候補の順位付けや検証状況の表は持ちません。
@@ -273,10 +286,20 @@ pub enum Method {
 
 ## prepare と evaluate
 
-M0 の `prepare` は `PrepareInput::Standard { product, trade_date }` または `PrepareInput::Ir(&Product)` と、担保の鍵、行使記録、市場構造、`Context` を受け取ります。
-IR を直接渡す場合は、簡約前の `Product` と、構築時と同じ参照データを使います。
+`prepare` は商品の入力 `PrepareInput` と、担保の鍵、行使記録、市場構造、`Context` を受け取ります。
+入力は三つで、方式の表を引く商品の種類もここから決まります。
+
+| `PrepareInput` | 内容 | 商品の種類 |
+|---|---|---|
+| `Standard { product, trade_date }` | 標準商品と約定日。準備のときに IR へ下ろす | `product.kind()` |
+| `Ir(&Product)` | 簡約前の検証済み IR。構築時と同じ参照データを使う | なし |
+| `Trade(&ResolvedTrade)` | 取引と、同じ参照データで下ろした IR の組（`TradeProductResolver::resolve` か `ResolvedTrade::new` だけが作る）。IR を作り直さずに使い、解析プライサーには取引の標準商品、約定日、手数料を渡す | `Trade::product_kind()` |
+
+`Session` は、読み込みで下ろした IR と取引を `PrepareInput::Trade` で渡します（ワークスペースの取引にも方式の表が当たります）。
+解析は取引を、リニア評価は IR を読むので、二つが別の取引から来ると方式を切り替えただけで別の取引を評価します。
+`ResolvedTrade` はフィールドを公開せず、取引を下ろす関数だけが組を作ることで、この取り違えを型で防ぎます（compile-fail テスト `resolved_trade_from_parts`）。
 `prepare_with_margin_set` は、担保の鍵の代わりに `MarginSet` を受け取り、VM から鍵を一度導いて同じ準備を行います。
-`prepare` と `evaluate` は取引一件の関数で、M0 の `Session` は取引の列を `rhodium_core::parallel::map_partitions` で実行器のワーカーに区間で分けて呼びます（[並列と再現性](../03-architecture.md#並列と再現性)）。
+`prepare` と `evaluate` は取引一件の関数で、`Session` は取引の列を `rhodium_core::parallel::map_partitions` で実行器のワーカーに区間で分けて呼びます（[並列と再現性](../03-architecture.md#並列と再現性)）。
 取引ごとの準備と PV は他の取引によらないので、結果は分け方によりません。
 取引と所属先を受け取る高水準の準備は、次の形に拡張します。
 
@@ -313,6 +336,9 @@ pub fn prepare(trade: &Trade, margin_set: &MarginSet, scripts: &ScriptLibrary,
 pub fn evaluate<R: Real>(p: &Prepared, market: &Market<R>) -> Result<MoneyBag<R>>;
 ```
 
+- M1 の `Prepared` は取引 ID を持たず（`Session` が取引 ID と組にします）、計画は `CashflowPlan` だけです。
+  金利先物の解析プライサーも、値洗いの列を同じ cashflow の計画で表し、最初の未決済の値洗いの予想を凸性の補正で置き換えます。
+  `Plan` の enum は M2 の解析プライサー（`AnalyticInput`）と M4 の `Kernel` で入れます。
 - `prepare` には取引と所属先、解析済みのスクリプトを渡します。
   直接呼ぶ場合も所属 ID と構造の評価時点が入力と一致することを確かめます。
   商品と手数料を同じ IR にまとめ、簡約、方式の選択、要求の解決、計画の作成、読む θ ブロックの収集を行います。
diff --git a/docs/design/07-pricing/02-linear-and-analytic.md b/docs/design/07-pricing/02-linear-and-analytic.md
index b46881669ba971aca44ef45432f1fdb71250cde3..ebef6049c714b0dd54837f08635940ca4341e66d 100644
--- a/docs/design/07-pricing/02-linear-and-analytic.md
+++ b/docs/design/07-pricing/02-linear-and-analytic.md
@@ -27,10 +27,10 @@
 | 恒等式 | 対象の項 | 評価式 | 厳密になる条件 | マイルストーン |
 |---|---|---|---|---|
 | `Fixed` | 定数 a の資産を時刻 t に受け渡す | a·P_d(t)（資産なら a·F_A(t)·P_d(t)） | なし | M0 |
-| `TermForward` | 定数 × ターム金利の観測 | N·τ·F_p(s, e)·P_d(t_pay) | 支払日 = e、予想と割引が同じ担保 | M1 |
+| `TermForward` | 定数 × ターム金利の観測 | N·τ·F_p(s, e)·P_d(t_pay) | 支払日 = e、予想と割引が同じ担保（[担保の違いの判定](#担保の違いの判定)） | M1 |
 | `RfrTelescoping` | 定数 × RFR の複利 | 係数 × Π_fix·Π_A(θ)·P_p(s_u)/P_p(e_obs)·P_d(t_pay) | 下の 5 条件、P_p = P_d、支払日が e_obs | M0 |
 | `RfrDailyProduct` | 定数 × 端点に畳めない成長率 | 日次丸め、SimpleFactor・PowerFactor、Inclusive を反映した成長率に定数と P_d(t_pay) を掛ける | 規約 `DailyForwardProduct` | M0 |
-| `FraSettlement` | FRA の決済額（期間の始まりに支払う） | N·(1 − (1 + τK)·P_p(e)/P_p(s))·P_d(s) | P_p = P_d。違えば規約 `FraDiscounting` | M1 |
+| `FraSettlement` | FRA の決済額（期間の始まりに支払う） | N·(1 − (1 + τK)·P_p(e)/P_p(s))·P_d(s) | P_p = P_d。違えば規約 `FraDiscounting`。予想と割引の担保が違えば `CrossCollateralProjection` も | M1 |
 | `FxForward` | 定数 × FX の観測かその逆数 | 支払通貨 c を相手とする向き X_{o→c} に正規化し、N·F_{o→c}(t → v)·P_c(t_pay) | 支払日 = v、同じ担保。違えば `PaymentDelayIgnored`、`CrossCollateralProjection` | M1 |
 | `FxRateProduct` | 定数 × FX の観測 × 金利の観測（MtM XCCY の元本リセット） | 二つのフォワードの積 | 規約 `MtmFxForwardProduct` | M1 |
 | `AssetForward` | 定数 × 資産価格の観測 | N·F_A(t)·P_d(t_pay) | 支払日 = 受渡日、同じ担保 | M8、M10 |
@@ -41,9 +41,15 @@
 
 - 公表済みの fixing と推定値は簡約で定数になっているため、恒等式は確定値と予想値を区別しません（[公表時刻の規則](../04-foundation/04-reference-data.md#fixingstore)）。
   評価時点をまたぐ未公表の観測は、曲線から予想する側です。
+- FRA の決済額 τ(L − K)/(1 + τL) は、分母と分子の L が同じ観測で、分母の τ が L の預金期間を指標の日数係数で測った長さと一致するとき（相対 1e-12 以内）だけ `FraSettlement` に当てます。
+  τ が違う独自契約では 1/(1 + τL) が期首から期末への割引因子にならないので、リニアでないと判定します。
+  期首払いとしての支払いの遅れの照合は、その決済額の中の観測にだけ当て、同じ支払いに並ぶ通常のタームクーポンは期末払いとして照合します（テスト `fra_identity_applies_to_its_own_settlement_and_accrual_only`）。
 - FX の観測の逆数 1/X_{c→o} は、逆向きのペアの観測 X_{o→c} に正規化します（例: NDF の N·(1 − K/X)）。
   支払通貨 c を相手とする向きの観測の c 払いの期待値は、c 建てのフォワード F_{o→c} = 1/F_{c→o} で、凸性の調整は要りません。
   支払日が受渡日 v と違うときは `PaymentDelayIgnored`、予想と割引の担保が違うときは `CrossCollateralProjection` を記録します。
+- FX の観測の受渡 v は、利用箇所によらず fixing 日からペアのスポットラグだけ進めた日の、ペアの慣行の受渡時刻です。
+  同じ指標・観測日の公表値は一つなので、MtM XCCY の元本リセットも通常の fixing と同じ F(t → v) で予想します。
+  新しい元本を交換する瞬間（期首の脚の通貨の決済時刻）は v と別に扱い、日付が違うとき（明示の `reset` 日）は `PaymentDelayIgnored` を記録します（テスト `mtm_xccy_reset_reads_the_fixing_at_the_pair_spot_delivery`）。
 
 ### RFR の telescoping
 
@@ -98,6 +104,28 @@ RFR の複利 Π_i(1 + w_i·r_i) を端点の割引因子の比 P_p(s)/P_p(e) 
 | `RfrCapletLmVariance` | RFR キャップレットの観測期間内のボラの形（線形に減ると仮定） | キャップの quote がこの規約で値付けされている |
 | `RoundingAtForward` | 丸めの非線形性（丸めた値の期待値と期待値を丸めた値の差） | 取引条件が丸めを宣言している |
 
+### 担保の違いの判定
+
+`CrossCollateralProjection` は、予想の担保と支払いの担保が違うときだけ記録します。
+予想曲線は指標ごとに一本で（[担保による割引](../06-market/02-curves.md#担保による割引)）、その予想は曲線を当てる校正商品の担保の下で定まります。
+この担保の割引曲線を、予想曲線の**予想の測度**と呼びます。
+
+- 校正は、予想曲線のブロックを当てる商品の binding のうち予想の担保を識別できるものがすべて、同じ割引曲線 D を指すとき、D を予想の測度として市場構造に記録します（`MarketStructure::with_forecast_measures`）。
+  先物以外の商品は、binding の担保で割り引く PV をパーにするので、価格設定の `discounting` の表で（曲線の通貨, binding の担保）から引いた割引曲線を指します。
+  金利先物の価格は値洗いの口座の測度で決まり、取引の担保の割引因子は打ち消し合うので（[金利先物](#金利先物)）、binding の担保は予想の担保を識別しません。
+  値洗いの口座の付利指標（`margin_index`）を持つ先物だけが、その指標の予想曲線を口座の測度の割引曲線として指し、持たない先物は根拠から除きます（テスト `futures_identify_the_forecast_measure_only_through_their_margin_index`）。
+  根拠が食い違う、割引の規則が無い、根拠が一つも無い（直接マークだけ、付利指標の無い先物だけ）予想曲線は測度が分からず、予想曲線そのものを割引にも使う場合だけを同じ担保とみなします。
+- 例: TONA 担保のベーシスで校正した TIBOR 6M の測度は JPY-TONA です。
+  TONA 担保の TIBOR のスワップは予想曲線と割引曲線の ID が違っても厳密で、25bp の付利の CSA の取引（TONA + 25bp で割り引く）と、JPY 担保の USD の支払いで割り引く SOFR の脚だけが規約を使います。
+- 恒等式ごとに比べるもの:
+
+| 恒等式 | 比べるもの | 違えば |
+|---|---|---|
+| `TermForward`、`RfrTelescoping` | 予想の測度と支払いの割引曲線 | `CrossCollateralProjection` |
+| `FraSettlement` | 予想曲線と割引曲線（凸性）、予想の測度と割引曲線（担保） | `FraDiscounting`、`CrossCollateralProjection` |
+| `FxForward` | 支払いの両通貨の割引曲線と、ペアの FX の担保区分の割引の組 | `CrossCollateralProjection` |
+| 金利先物の値洗い | 値洗いの口座の指標の予想と割引曲線（二つの率が等しいこと。[金利先物](#金利先物)） | `CrossCollateralProjection` |
+
 - 閾値は `LinearConventions` の値で、結合モデル（LGM など）との差の測定で見直します。
   M0 では、ずれを支払通貨の決済カレンダーで数え、使用した近似を `Prepared::warnings` で返します。
   使った恒等式と規約は cashflow の表の観測の列に記録します。
@@ -255,6 +283,9 @@ PV = N·P_d(s_next)·(F(A) − P_last)
 - 後決め複利先物の最終決済値は 1 − R、R = (Π − 1)/τ_a で、Π は観測期間の 1 日物の率の複利です。
   モデルの 1 日物の率は 1 日の割引債から作るので（[モデル](03-models.md)）、先物価格は 1 − (E^{Q_M}[Π_i 1/P(d_i, d_{i+1})] − 1)/τ_a です。
   連続時間の短期金利の積分ではなく IR の観測と同じ離散の積を使うので、解析式と MC は同じ確率変数の期待値です。
+- 評価時点 A が未公表の 1 日物の観測の期間をまたぐときは、[A, 終了) の forward の率を残りの日数係数 τ_rem で年率化し、その率を元の日の日数係数に当てます（[fixing の規則](../04-foundation/04-reference-data.md)、RFR の telescoping と同じ規則）。
+  後決め複利ではその日の因子が 1 + τ·(G(A, 終了) − 1)/τ_rem、加重平均では重みが τ/τ_rem 倍になります。
+  G(A, 終了) は評価時点で決まるので、凸性の項に入りません（テスト `overnight_futures_annualize_the_spanning_day_over_its_remaining_fraction`）。
 - 値洗いを連続とみなす式は Q_M の式の極限で、外部実装との照合の目安にだけ使います。
   日次との差は 1 日の分散の次数で、照合の許容誤差にこの差の見積もりを含めます（テスト `futures_continuous_marking_gap`）。
 
@@ -316,7 +347,8 @@ par   = Prot / (N·RPV01)
 ### D1 リニア評価は恒等式の照合で行い、端点の比は条件を満たすときだけ使う
 
 - 判断: 線形な IR は項ごとに恒等式を照合し、RFR の telescoping は 5 条件を満たすときだけ使います。
-- 理由: lookback、lockout、支払いの遅れの下で端点の比を使うと値が誤り、校正と取引の両方が同じ近道を使うとその誤りが見えなくなります。
+- 理由: lookback と lockout の下で端点の比を使うと値が誤り、校正と取引の両方が同じ近道を使うとその誤りが見えなくなります。
+  支払いの遅れは 5 条件に含めず、端点の比のまま支払日の割引で評価し、タイミングの調整を省いたことを名前付きの近似 `PaymentDelayIgnored` として記録します（遅れの上限を超えれば `Unsupported`）。
 - 採らない案: 商品ごとに専用のスワップの公式を書く案は、商品の意味が IR と二か所に分かれるため採りません。
 
 ### D2 近似は名前付きの規約として許可制にする
@@ -337,3 +369,11 @@ par   = Prot / (N·RPV01)
 - 判断: 反復を tape に記録せず、収束した点の陰関数定理で 1 命令として記録します。
 - 理由: tape が反復の回数によらず小さく、微分が打ち切りに左右されません。
 - 採らない案: 反復をそのまま AAD で微分する案は採りません。
+
+### D5 担保の違いは予想の測度で判定する
+
+- 判断: `CrossCollateralProjection` の要否を、予想曲線の予想の測度（曲線を当てる校正商品の担保の割引曲線）と支払いの割引曲線の比較で決めます。
+  測度は校正が binding と価格設定から導き、市場構造に記録します。
+- 理由: 予想曲線と割引曲線の ID を比べると、同じ担保の TIBOR のスワップやベーシスまで近似として報告し、規約を禁止すると厳密に評価できる取引と校正商品まで評価できなくなります。
+- 採らない案: 予想曲線と割引曲線の ID の比較は、複数曲線のすべての取引で誤って警告するため採りません。
+  曲線ごとに測度を市場の定義に書かせる案は、binding の担保と食い違いうる入力が増えるため採りません。
diff --git a/docs/design/08-risk/01-sensitivities.md b/docs/design/08-risk/01-sensitivities.md
index 21dd22b6b160ffdb6a6ff44550edd3bf6ff345f4..097bc86f3a467c2bc57937fac647baf55ffeb2f4 100644
--- a/docs/design/08-risk/01-sensitivities.md
+++ b/docs/design/08-risk/01-sensitivities.md
@@ -59,7 +59,7 @@ flowchart LR
 ### 1 θ 勾配
 
 - `Prepared` が知っている、取引が読む θ ブロックだけを `Var<'t, f64>` の変数にし、残りの θ は定数のまま評価します（active set、[アーキテクチャ](../03-architecture.md#数値型と自動微分)）。
-  M0 は残りの θ も取引ごとに定数として写して市場を作り直します（M0 の制限と後続の変更は [曲線と market](../06-market/02-curves.md#構造と-θ)）。
+  `Market<f64>::with_active_theta` は校正済みの構造、θ、利用可否を共有し、passive なブロックを取引ごとに写しません（[曲線と market](../06-market/02-curves.md#構造と-θ)）。
 - reverse sweep で通貨ごとの勾配 g_c = ∂V_c/∂θ を得て、通貨ごとの PV と一緒に `ParamGradient` にします。
   勾配はブロックの単位で疎であり、TONA と TIBOR だけを読む JPY スワップは二本の曲線のブロックにだけ値を持ちます。
 - tape はワーカーごとに一つを使い回し、取引ごとに並列に計算します。
@@ -124,7 +124,10 @@ for k in reverse_topological(groups):
 - 上流へ送る項により、TIBOR 曲線の感応度の一部が TONA の OIS quote の感応度になります。
 - 同じ quote を複数のグループが読むことがあるため（ATM との差の quote が ATM を参照する場合）、quote への随伴は加算で累積します。
 - 最小二乗のグループを通した 1 次感応度が厳密になるのは、∂G/∂θ と ∂G/∂θ_up を 2 階の項まで含めて記録したときです。
-  Gauss–Newton の記録では近似になり、結果にその旨を記録します。
+  校正は目的 Φ = ½‖r̃‖² の Hessian–vector product を自グループの θ の方向だけ求め、∂G/∂θ と（∇²Φ の対称性から）∂G/∂θ_up を同時に記録します（[解の記録](../06-market/03-calibration.md#解の記録)）。
+  記録の型は正方系と同じ（∂G/∂θ の LU 分解と疎な ∂G/∂q、∂G/∂θ_up）なので、随伴の式は変わりません。
+  Gauss–Newton の記録では近似になり、`GroupSolution::hessian` と校正の警告 `Warning::CalibrationGaussNewton`（グループが読む quote ごと）にその旨を記録します。
+  残差の残る当てはめでは、Gauss–Newton の随伴のデルタは再校正のバンプと一致しません（テスト `penalized_meeting_steps_quote_delta_matches_bump` の対照）。
 - シナリオの上書き（`MarketShift::Override`）は校正グループではなく、シフトと同じく読み出しに重ねる定数です（[シナリオと P&L](02-scenarios-and-pnl.md#市場のシフト)）。
   上書きした読み出しへの価格の勾配は 0 として 2 段目の Ḡ に入れず、校正の随伴はそのまま流します。
   そのため、上書きしたブロックを上流に持つ下流のグループは、再校正しない限り元の quote への依存を保ちます。
diff --git a/docs/design/10-io-and-execution.md b/docs/design/10-io-and-execution.md
index 5daa05f7969f64eabd94336a35a15abaa3e750c3..8b370b43475cf467d42fdc0fe03d78318d93dcfb 100644
--- a/docs/design/10-io-and-execution.md
+++ b/docs/design/10-io-and-execution.md
@@ -20,9 +20,9 @@ workspace/
     zones.toml             業務で使う地域の UTC オフセットと遷移表
     holidays/JPTO.csv      カレンダーごとの休日の列
     indices.toml           指標の定義（種類、通貨、テナー、日数計算、カレンダー、fixing のラグ）
-    conventions.toml       商品の慣行（OIS、IRS、FX スポットなど）
+    conventions.toml       商品の慣行（OIS、IRS など）と通貨ペアの受渡慣行（[[fx]]）
   market/
-    definition.toml        市場の定義（曲線とボラ面の構造、校正商品の指定、校正グループ）
+    definition.toml        市場の定義（曲線とボラ面の構造、校正商品の指定、最小二乗の重みと罰則）
   quotes/
     snapshot.csv           呼び出し側が選んだ評価時点の quote 一式
   fixings/
@@ -49,26 +49,29 @@ workspace/
 - rhodium はワークスペースの過去の版を管理しません。
   訂正した値を使うときは、呼び出し側がファイルを書き換えるか、訂正後のスナップショットを渡します（[全体像](01-overview.md#d2-データの履歴は呼び出し側が持つ)）。
 
-### M0 のワークスペース
+### M1 のワークスペース
 
-[OIS のサンプル](../../examples/m0-ois/) は M0 が読む TOML と CSV の一式です。
+[OIS のサンプル](../../examples/m0-ois/) は M0 の、[金利リニアのサンプル](../../examples/m1-rates/) は M1 の `Workspace::load` が読む TOML と CSV の一式です。
 各 TOML 文書は `version = 1` を持ちます。
-M0 の `Workspace::load` は TOML だけを読み、同じ型の JSON の代替（[共通の規則](#共通の規則)）は M6 のファイル形式の見直しまでに足します。
+`Workspace::load` は TOML だけを読み、同じ型の JSON の代替（[共通の規則](#共通の規則)）は M6 のファイル形式の見直しまでに足します。
+`collateral/`、`products/`、`config/xva.toml`、`regulatory/` は、それぞれを使うマイルストーン（M5、M4、M5、M8）で読みます。
 
-| ファイル | M0 で読む項目 |
+| ファイル | M1 で読む項目 |
 |---|---|
 | `refdata/currencies.toml` | `[[currency]]` の通貨定義 |
 | `refdata/calendars.toml` | `[[calendar]]` の `id`、既知範囲の両端 `from` と `through`、`weekend` の曜日、任意の `extrapolation`。既定は `strict` |
 | `refdata/holidays/<id>.csv` | 各基本カレンダーの休日。休日が無い場合も `date` ヘッダーを置く |
 | `refdata/zones.toml` | `[[zone]]` の `id`、UTC の適用範囲 `from` と `through`、`initial_offset_seconds`、任意の `transitions`（`at` と `offset_seconds`） |
-| `refdata/indices.toml`、`refdata/conventions.toml` | `[[index]]` の指標定義、`[[convention]]` の `id` と `kind`（`ois` または `deposit`） |
-| `market/definition.toml`、`quotes/snapshot.csv` | 明示した UTC ピラーと OIS の binding、`quote_id,value` のスナップショット。形式は[校正](06-market/03-calibration.md#m0-の市場定義)で定める |
+| `refdata/indices.toml` | `[[index]]` の指標定義（翌日物、ターム金利、FX の fixing、金利先物の清算価格）と `[[fallback]]` |
+| `refdata/conventions.toml` | `[[convention]]` の `id` と `kind`（`deposit`、`ois`、`irs`、`basis_swap`、`xccy_swap`、`fra`、`fx_forward`、`rate_future`）と、通貨ペアの受渡慣行 `[[fx]]`（[FX の慣行](04-foundation/04-reference-data.md#fx-の慣行)） |
+| `refdata/securities.toml` | 任意。`[[security]]` の証券の定義（金利先物の限月の上場契約など） |
+| `market/definition.toml`、`quotes/snapshot.csv` | 曲線（明示した UTC ピラーか部品の式 `expr`）、校正の binding、`[[fx]]`、`[lgm.*]`、最小二乗の指定と、`quote_id,value` のスナップショット。形式は[校正](06-market/03-calibration.md#ワークスペースの市場定義)で定める |
 | `fixings/<id>.csv` | 指標ごとの公表済み fixing。ヘッダーは指標の公表規則に従う |
 | `agreements.toml`、`trades/*.toml` | `[[margin_set]]` の `id` と `collateral`。各 `[[trade]]` は下の例と同じ形で標準商品、約定日、ブック、マージンセットを指定する |
-| `config/pricing.toml` | [M0 の価格設定](07-pricing/01-pricing.md#m0-で読み込める設定) |
+| `config/pricing.toml` | [M1 の価格設定](07-pricing/01-pricing.md#m1-で読み込める設定) |
 
 マージンセットには、現金担保の条件、CTD の候補（`cheapest_to_deliver`）か `"uncollateralized"` を指定します。
-書式は価格設定の `discounting` の `collateral` と同じで、現金担保には任意の `floor` を書けます（[M0 の価格設定](07-pricing/01-pricing.md#m0-で読み込める設定)）。
+書式は価格設定の `discounting` の `collateral` と同じで、現金担保には任意の `floor` を書けます（[M1 の価格設定](07-pricing/01-pricing.md#m1-で読み込める設定)）。
 ワークスペースの読み込み時に、標準商品の約定条件と参照先を解決して IR を作ります。
 取引の準備では、その IR を評価時点で利用できる fixing と決済条件に基づいて簡約します。
 
@@ -253,6 +256,15 @@ let bumped = session.with_quotes(&overrides, execution)?;
 
 ### 選択や再校正を行う API
 
+この節は 1.0 までに足す API を含みます。
+M1 までに実装したのは `with_quotes` と上の現行の API で、残りは次のマイルストーンで足します。
+
+| API | マイルストーン |
+|---|---|
+| `Selection`、`DeltaOptions`、`cashflows`、`gamma`、`scenarios`、`theta`、`explain`、`with_as_of` | M3（リスク計算ツール） |
+| `xva`、`incremental_xva`、`with_trades` | M5（XVA）、取引前の評価の高速化は M6 |
+| `with_config`、公開 API の見直し | M6 |
+
 ```rust
 let ws = Workspace::load("workspace")?;
 let as_of = ws.resolve_as_of("2026-11-20T17:00:00+09:00")?;
@@ -267,7 +279,7 @@ let xva = session.xva(&Selection::netting_set("NS-ACME-1"), execution)?;  // Xva
 
 | メソッド | 返すもの | 文書 |
 |---|---|---|
-| `pv` | 取引ごとの `Valuation`（通貨ごとの PV、報告通貨の PV、方式）と合計、除いた取引の一覧 | [価格計算](07-pricing/01-pricing.md) |
+| `pv` | 取引ごとの `Valuation`（通貨ごとの PV、報告通貨の PV、方式）と合計、除いた取引の一覧。M1 までは選択を取らず、受渡通貨別の `pv` と報告通貨の `value_portfolio` に分かれる | [価格計算](07-pricing/01-pricing.md) |
 | `cashflows` | 取引ごとのキャッシュフロー表（受渡時刻、通貨、期待額、DF、PV、確定か予想か推定か） | [価格計算](07-pricing/01-pricing.md#cashflow-の表) |
 | `quote_delta` | `DeltaOptions` の集計単位ごとの `RiskVector<QuoteId>` | [感応度](08-risk/01-sensitivities.md) |
 | `gamma` | 指定した方向かブロックの 2 次の感応度 | [感応度](08-risk/01-sensitivities.md) |
@@ -306,12 +318,14 @@ let xva = session.xva(&Selection::netting_set("NS-ACME-1"), execution)?;  // Xva
 
 `rh` はワークスペースを読み、結果を `output/` か標準出力に書く薄いコマンドです。
 
-### M0 のコマンド
+### M1 までのコマンド
 
 - `rh price <dir>` は `trade_id,currency,pv` を取引 ID、通貨の順に標準出力へ書きます。
 - `rh risk <dir>` は `quote_id,currency,sensitivity` を quote ID、通貨の順に標準出力へ書きます。感応度は全取引の通貨別の合計で、[quote の種類が定める変化幅](06-market/01-quotes.md#感応度の単位)あたりの額です。
 - どちらも `--as-of` または `--date` を指定します。地域を使う場合は `--as-of '2026-11-24T17:00:00[Asia/Tokyo]'` の形で指定します。`--date 2026-11-24` は価格設定の `valuation_close` で瞬間に解決します。
 - `--threads N`（既定 1）は読み込みから評価まで取引を N 個のワーカーに分けます。1 は逐次で、0 は引数の誤り（終了コード 2）です。取引ごとの PV は N によらず同じで、デルタは通貨別の和の順だけが変わるので [許容誤差](11-quality.md#許容誤差) の範囲で一致します。
+- 準備と校正の警告（近似の規約、校正の当てはめの許容誤差の超過など）は `warning: trade <ID>: <警告>` と `warning: calibration: <警告>` の形で標準エラーに書き、終了コードには影響しません。
+- `--format`、`--by` と `check`、`xva` のコマンドは下の拡張で足します（M3 以降。CLI の形は M6 の公開 API の見直しで確定します）。
 - 診断は標準エラーに書きます。校正や取引の評価が部分的に失敗したときは、計算できた行を出力して終了コード `1` を返します。一つの通貨の校正商品の準備の失敗（fixing の欠損など）も、その曲線のグループの失敗として扱い、独立した通貨の行を出力します。入力のエラーは `2`、すべて成功した場合は `0` です。評価時点が曲線の最初のピラー以降にあるなど、入力どうしの不整合で市場を組み立てられない場合も入力のエラーとして扱い、CSV を出力しません。
 
 ### コマンドの拡張
@@ -344,7 +358,7 @@ let xva = session.xva(&Selection::netting_set("NS-ACME-1"), execution)?;  // Xva
 - 並列は一段だけにします。
   対象の数（取引、ネッティングセット）がワーカー数の 2 倍以上なら対象で分け、少なければパスのバッチで分けます。
   パスを持たない計算（取引の準備、リニア評価、AAD の θ 勾配）は常に対象で分け、区間の数は min(ワーカー数, 対象の数) です（[並列と再現性](03-architecture.md#並列と再現性)）。
-  M0 の `Session` は、取引の準備、PV、quote デルタの θ 勾配を、取引 ID の順の連続した区間に分けます。
+  `Session` は、取引の準備、PV、quote デルタの θ 勾配を、取引 ID の順の連続した区間に分けます。
 - 取消と期限（実行文脈の `CancellationToken` と `Deadline`）は、区間の開始前に加えて取引の境界と AD の記録・逆伝播の中で確かめます。
   `Session::new` は校正の前にも確かめ、取引が 0 件の準備・PV・quote デルタも並列の入口で確かめます（取消した実行が空の結果で成功しないように）。
   止まれば、取引の失敗として結果に残さず、呼び出し全体を `Execution` にします（[並列と再現性](03-architecture.md#並列と再現性)）。
@@ -352,13 +366,14 @@ let xva = session.xva(&Selection::netting_set("NS-ACME-1"), execution)?;  // Xva
   上位の単位の指標は、メンバーのパス上の量を同じバッチで足してから非線形の関数を当てるためです（[XVA の指標](09-xva/02-metrics.md#集計の単位)）。
   入れ子の並列は、ワーカーの取り合いとメモリの増加を招くため使いません。
 - ワーカーは、作業領域 `Scratch`（評価用のバッファ）と AD の tape を一つずつ持ち、対象をまたいで使い回します。
-  M0 は tape をワーカーごとに一つ持ちます。
+  M1 までは tape をワーカーごとに一つ持ちます。
   リニア評価の小さなバッファは取引ごとに確保し、1 万件の PV は数十 ms で済むため（[M0 の計測](../benchmarks/m0.md)）、`Scratch` はホットループのある MC のカーネル（M4）とともに入れます。
 - 取引をまたいで共有する値（指標の ID の `Arc<str>`、市場の構造の `Arc`）を、日次の観測の行のような細かい単位で複製しません。
   共有の参照数の更新はワーカーの間で競合し、並列の伸びを止めます（M0 の準備で 20 スレッドの時間が半分になった例は [M0 の計測](../benchmarks/m0.md)）。
 - 縮約の順序は固定しません。
   スレッド数を変えた結果の一致は [許容誤差](11-quality.md#許容誤差) で判定します。
 - 校正はグループの依存の順に解き、依存し合わないグループは並列に解けます。
+  M1 までの校正は逐次です（最初の検証データセットの全 9 グループの校正が 1 コアで目安の 100 ms を大きく下回るため、[M1 の計測](../benchmarks/m1.md)）。グループの並列はボラ面の校正が重くなる M2 以降で測って入れます。
 
 ## メモリ
 
diff --git a/docs/design/11-quality.md b/docs/design/11-quality.md
index 435975216feb37ea47f372a432cca6d9d850d586..4cbfe65b10bb409243dfc961b94596802fce7f04 100644
--- a/docs/design/11-quality.md
+++ b/docs/design/11-quality.md
@@ -12,7 +12,7 @@
 | 単体テスト | 関数の振る舞いと境界の検証（不正な入力のエラー） | 各モジュールの `#[cfg(test)]` か crate の `tests/` |
 | 手計算の例 | 式と規約（日数計算、スケジュール、CSA の額、close-out）が手で計算した値に一致する | 例の導出をテストのコメントに書く |
 | 性質テスト | 法則（スケジュールの不変条件、カレンダーの合成、`MoneyBag` の和、随伴の線形性） | proptest |
-| compile-fail テスト | 型で禁止した演算がコンパイルできず、正しい使い方はコンパイルできる | trybuild の `_ok` と失敗の対 |
+| compile-fail テスト | 型で禁止した演算がコンパイルできず、正しい使い方はコンパイルできる | testkit の `_ok` と失敗の対 |
 | 一致のテスト | 同じ商品の二つの評価（IR のリニア評価と解析プライサー、MC と閉形式）が一致する | [オラクル](#オラクル) |
 | 外部の fixture | QuantLib、ORE などの値と一致する | 事前に生成した値を `tests/fixtures/` に置く |
 | AAD とバンプ | AAD の感応度が、再校正込みのバンプの差分に一致する | `check_aad_vs_bump` |
@@ -37,7 +37,7 @@
 | `assert_close!(x, y, rel, abs)`、`assert_close!(x, y, rel, abs, scale)` | `abs(x − y) ≤ max(rel·M, abs)` を確かめる。M は `max(abs(x), abs(y))` で、5 引数の形は和の項の絶対値の和などの `scale` も M に含める（[許容誤差](#許容誤差)） |
 | `assert_mc!(estimate, se, reference, k, a)` | `abs(estimate − reference) ≤ k·se + a` を確かめる |
 | `check_aad_vs_bump(reprice, aad, unit, pv_error, absolute, label)` | quote の増分（自然な単位）から再校正込みの価格を返す閉包、報告単位あたりの AAD の感応度、報告単位の幅 h₁、再校正の PV の誤差 δ_V、絶対の下限を受け取り、下の手順で照合する |
-| `compile_fail_pair!("case")` | `tests/ui/case_ok.rs` が通り、`tests/ui/case.rs` がコンパイルに失敗することを trybuild で確かめる |
+| `compile_fail_pair!("case")` | `tests/ui/case_ok.rs` を実行し、`tests/ui/case.rs` のコンパイル診断を `case.json` と照合する |
 | `TickClock` | 読むたびに 1 ns 進む決定的な時計。`TickClock::after(n)` を期限にすると n 回目の読み出しで止まるので、区間の開始後の取消と期限の伝わり方をスレッドの速さによらず確かめられる |
 
 ## オラクル
@@ -75,14 +75,19 @@
 | 地域の時刻（M0） | Europe/Berlin の春の欠落は入力エラー。秋の重複は明示オフセットにより 1 時間違う二つの瞬間になる | [日付とカレンダー](04-foundation/01-dates-and-calendars.md#瞬間と現地時刻) |
 | RFR の観測シフト（M1） | 観測期間 1 日・計算期間 3 日、元本 3,650 円、10% の支払いは 3 円 | [スケジュール](04-foundation/02-schedules.md#参照の式) |
 | 非営業日に始まる RFR の期間（M1） | 土曜に始まり火曜に終わる期間の因子は (1 + r_金·2/365)(1 + r_月·1/365) で、土日に金曜の値が掛かり全日を覆う（テスト `rfr_period_starting_on_holiday`） | [スケジュール](04-foundation/02-schedules.md) |
-| 金利先物の値洗い（M1） | 値洗いの列の IR のリニア評価は次の未決済の値洗いだけを残し、PV = N·P(s_next)·(F − P_last) が解析の価格と一致する | [リニアと解析](07-pricing/02-linear-and-analytic.md#金利先物) |
-| 先物の値洗いの測度（M4） | 確率金利で先物価格と割引が相関する LGM で、値洗いの列の IR の MC が Q_M の解析の価格と標準誤差の内で一致する（`futures_margin_ir_vs_analytic`）。最終日に一度だけ払う契約は凸性の分だけ違い（`futures_vs_forward_settlement`）、連続の値洗いの式との差は見積もった日次の差の内に収まる（`futures_continuous_marking_gap`） | [リニアと解析](07-pricing/02-linear-and-analytic.md#金利先物) |
+| 金利先物の値洗い（M1） | 値洗いの列の IR のリニア評価は次の未決済の値洗いだけを残し、PV = N·P(s_next)·(F − P_last) が解析の価格と一致する（`futures_margin_ir_vs_analytic`）。清算の瞬間と値洗いの受渡がずれると `PaymentDelayIgnored` を記録し、禁止すると評価しない（`futures_settlement_delay`） | [リニアと解析](07-pricing/02-linear-and-analytic.md#金利先物) |
+| 期間をまたぐ先物の 1 日（M1） | 評価時点 20 時が未公表の翌日物の期間（16 時〜翌 16 時）をまたぐと、[A, 終了) の成長を残りの 20/24 日で年率化し、元の 1 日に当てる。複利と加重平均の価格が手計算と一致し、IR と解析も一致する（`overnight_futures_annualize_the_spanning_day_over_its_remaining_fraction`） | [リニアと解析](07-pricing/02-linear-and-analytic.md#金利先物) |
+| FRA の恒等式（M1） | 一つの支払いに FRA の決済と通常のタームクーポンを並べても、期首払いの扱いは FRA の観測にだけ当たり、近似を記録しない。分母の τ を Act/360 で数えた独自契約はリニアでない（`fra_identity_applies_to_its_own_settlement_and_accrual_only`） | [リニアと解析](07-pricing/02-linear-and-analytic.md) |
+| 先物の値洗いの測度（M4） | 確率金利で先物価格と割引が相関する LGM で、値洗いの列の IR の MC が Q_M の解析の価格と標準誤差の内で一致する（`futures_margin_mc_vs_analytic`）。最終日に一度だけ払う契約は凸性の分だけ違い（`futures_vs_forward_settlement`）、連続の値洗いの式との差は見積もった日次の差の内に収まる（`futures_continuous_marking_gap`） | [リニアと解析](07-pricing/02-linear-and-analytic.md#金利先物) |
 | 明示の期間のスタブ（M1） | 半年払いで 1 か月の短いスタブは、`reference` に通常の半年の期間を書けば ACT/ACT ICMA の日数係数がスタブの日数 / (2 × 参照期間の日数) になり、省略すると leg の構築が `InvalidData` になる（`explicit_stub_reference_period`） | [スケジュール](04-foundation/02-schedules.md#明示の期間) |
 | 曲線と点の組（M0） | 予想曲線の点を割引曲線と組にする binding は `InvalidInput`、組を作らずに読むコードと別のスコープの `PayAt` はコンパイルできない（`bind_rejects_point_of_other_curve`、`read_with_unbound_point`、`pay_point_cross_scope`） | [曲線と market](06-market/02-curves.md#market-view) |
 | 資金調達曲線の一か所の指定（M1） | 資金調達の `InDiscounting` で、無担保の割引の行が `funding_curves` の `borrow` と違えば `prepare` が `InvalidInput`（`funding_curves_single_source`） | [価格計算](07-pricing/01-pricing.md#pricingconfig) |
 | FX の短期の連鎖（M1） | T+2 のペアはスポットから TN、ON を、T+1 のペアは ON だけを実際の受渡日で連鎖させ、当日の受渡の締めの後は翌営業日の受渡から割り戻す。スポット 1.300、ON 0.001 の T+1 のペアの当日の受渡のレートは 1.299（`fx_x0_short_dated_chain`） | [曲線と market](06-market/02-curves.md#fx-のフォワード) |
-| 先物の受渡の遅れ（M4） | 清算価格を決める瞬間と値洗いの受渡が一致すれば IR・解析・MC が一致し、翌朝の受渡では差が遅れの近似の内に収まる（`futures_settlement_delay`） | [リニアと解析](07-pricing/02-linear-and-analytic.md#金利先物) |
+| 先物の受渡の遅れ（M4） | 清算価格を決める瞬間と値洗いの受渡が一致すれば IR・解析・MC が一致し、翌朝の受渡では差が遅れの近似の内に収まる（`futures_settlement_delay_mc`） | [リニアと解析](07-pricing/02-linear-and-analytic.md#金利先物) |
 | 非線形のパーの微分（M1） | PV = p² − θ の θ = 1 で d²p/dθ² = −1/4（`par_second_derivative`）。残差が 0 でない最小二乗の校正で、非線形のパーの quote デルタが再最適化の有限差分と一致する（`par_nonlinear_lsq_quote_delta`） | [リニアと解析](07-pricing/02-linear-and-analytic.md#パーレートの解) |
+| 最小二乗の単位と採用（M1） | 同じ受渡の二つのアウトライト q₁、q₁ + 2 を σ = 0.01、1 で当てると F* = q₁ + 2/10001 で、外れた q₂ だけを警告し、quote デルタが再校正の差分と一致する（`fx_outright_least_squares_measures_residuals_in_price_units`）。σ を付けた FX スポットの直接マーク −1 は解く前に `InvalidInput`（`least_squares_direct_marks_check_quotes_before_fitting`）。Gauss–Newton の解も厳密なヘシアンが正定値のときだけ採用する（`acceptance_requires_a_positive_definite_hessian`） | [校正](06-market/03-calibration.md#残差の形) |
+| 外貨担保曲線の長期（M1） | T+1 のペアの spot、ON と長期の XCCY だけの定義が当日の締めの前後で校正でき、XCCY の quote を再現する（`fx_implied_long_end_alone_survives_the_cutoff_of_a_t_plus_one_pair`）。EUR と JPY の XCCY の慣行を USDJPY の長期に書くと `InvalidData`（`fx_implied_long_end_rejects_an_xccy_convention_of_another_pair`） | [曲線と market](06-market/02-curves.md#fx-のフォワード) |
+| FX の慣行の正本（M1） | MtM XCCY のリセットの観測日と FX fixing の受渡日は参照データのペアのスポットラグで決まり、JPY の祝日を含む暦で T+1 と T+2 の日付と PV が手計算と一致する。USDCAD に TN を足した定義は T+1 の慣行では連鎖せず、慣行だけを T+2 にすると校正できる（`mtm_xccy_reset_fixing_dates_follow_the_pair_spot_lag`、`fx_fixing_spot_value_uses_the_pair_spot_lag`、`fx_x0_chain_follows_refdata_spot_lag`）。元本リセットも通常の fixing と同じくペアのスポットの受渡（東京 15 時）の F で予想し、新しい元本の交換（USD の決済時刻）とは分ける。明示の reset 日でスポットの受渡が交換日より早ければ `PaymentDelayIgnored` を記録する（`mtm_xccy_reset_reads_the_fixing_at_the_pair_spot_delivery`） | [参照データ](04-foundation/04-reference-data.md#fx-の慣行) |
 | X(0) の正本（M1） | ON・TN の quote があり、市場の担保のフォワードが `FxImplied` で ON・TN・スポットに校正されない定義は `InvalidData`（`fx_x0_on_tn_consistency`）。直接のクロスの quote と三角換算が食い違えば警告を出し、直接の quote を使う（`fx_cross_direct_vs_triangulated`） | [曲線と market](06-market/02-curves.md#fx-のフォワード) |
 | ボラの種類と式（M2） | 正規ボラの面を Black だけの式に、対数正規の面を Bachelier だけの式に解決すると `prepare` が `InvalidInput`（`vol_kind_formula_mismatch`） | [ボラティリティ](06-market/04-volatility.md#読み出し) |
 | 相対尺度のガンマ（M3） | V(q) = q のスポットガンマは 0。V(q) = q² なら相対 1% の二軸で 2·(0.01q)² | [感応度](08-risk/01-sensitivities.md#ガンマ) |
@@ -90,8 +95,12 @@
 | smile dynamics の付け直し（M2） | 市場の残差が 0 でない過剰決定の SABR で、`StickyStrike` と `StickyMoneyness` の B が基準の面のモデルのボラへ付け直した有限差分と一致する（`smile_dynamics_overdetermined_refit`） | [感応度](08-risk/01-sensitivities.md#smile-dynamics) |
 | Philox の時刻追加（M4） | 同一 seed で 10 時の一点を求めた場合と 15 時も求めた場合の、10 時のモデル状態が許容誤差内で一致 | [モンテカルロ](07-pricing/04-monte-carlo.md#乱数) |
 | 日中の VM（M5） | コール 10 時と受渡 15 時の直前・ちょうど・直後で未決済と保有を照合し、同時刻の別周期のコールと cutoff が決済後の保有を読む | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md#イベントの順序) |
-| 部分的な校正（M1） | USD の quote を除くと USD を読む取引だけが失敗する。補って再校正すれば保存済みの準備を再び使え、全体を解き直した価格・デルタと一致する | [校正](06-market/03-calibration.md#失敗) |
-| 同じ通貨で担保条件が違う取引（M1） | 付利スプレッドが 0bp と 25bp の CSA を対応する別々の割引曲線に解決する。必要な規則が欠ければその取引をエラーにする | [価格計算](07-pricing/01-pricing.md#要求の解決) |
+| 部分的な校正（M1） | USD の quote を除くと USD を読む取引だけが失敗する。補って再校正すれば保存済みの準備を再び使え、全体を解き直した価格・デルタと一致する（`portfolio_isolates_trade_failures_and_totals_successes`） | [校正](06-market/03-calibration.md#失敗) |
+| 同じ通貨で担保条件が違う取引（M1） | 付利スプレッドが 0bp と 25bp の CSA を対応する別々の割引曲線に解決する。必要な規則が欠ければその取引をエラーにする。取引のマージンセットの付け替えでも、同じマージンセットの CSA の条件の変更でも準備し直し、新しい割引曲線の値になる（`session_reprepares_for_intraday_instant_fixings_and_collateral`） | [価格計算](07-pricing/01-pricing.md#要求の解決) |
+| 予想の測度（M1） | TONA 担保のベーシスで校正した TIBOR 6M の測度は JPY-TONA で、TONA 担保の TIBOR の取引は `CrossCollateralProjection` を記録せず、規約を禁止しても評価できる。25bp の CSA の取引と JPY 担保の SOFR の脚は記録する。担保の違う binding が混ざれば測度は記録しない（`tibor_swap_under_calibration_collateral_needs_no_cross_collateral_convention`、`tibor_swap_under_a_spread_csa_records_cross_collateral_projection`、`sofr_leg_under_jpy_collateral_records_cross_collateral_projection`、`forecast_measure_is_unknown_when_bindings_disagree`）。先物だけで当てた曲線は binding の担保で測度を決めず、別の割引曲線で割り引くクーポンは記録する。付利指標を割引曲線で予想する先物なら記録しない（`futures_identify_the_forecast_measure_only_through_their_margin_index`） | [リニアと解析](07-pricing/02-linear-and-analytic.md#担保の違いの判定) |
+| 固定元本の XCCY（M1） | 各通貨の脚が期首に元本を払い期末に受け取り、クーポンと元本をそれぞれの通貨の曲線で割り引いた和に一致し、`MtmFxForwardProduct` を記録しない（`fixed_notional_xccy_exchanges_both_currencies`） | [標準商品](05-products/03-standard-products.md) |
+| 近似の規約の記録（M1） | 期末より後に払うターム金利のクーポンは `PaymentDelayIgnored` を記録して支払日で割り引き（`term_coupon_payment_delay_is_recorded`）、別の担保で受け取る FX の fixing は担保区分の予想のまま割引だけが変わって `CrossCollateralProjection` を記録し（`fx_fixing_under_another_collateral_records_cross_collateral_projection`）、前の二つは規約を禁止すると `Unsupported`。値洗い口座の付利指標を割引と同じ曲線で予想する先物は記録せず、別の曲線で予想すると記録する（`futures_margin_index_cross_collateral_projection_is_recorded`） | [リニアと解析](07-pricing/02-linear-and-analytic.md#近似の規約) |
+| 最初の検証データセット（M1） | `examples/m1-rates` の 9 グループが残差 1e-14 で解け、18 取引のすべてを評価し、近似の警告が予想の測度に従い、一部解約の PV が元の IRS の 0.4 倍と手数料の PV の和に一致する（`m1_rates_workspace_calibrates_and_values_every_trade`、`amended_trade_matches_the_scaled_original_and_its_fee`） | [M1 のサンプル](../../examples/m1-rates/README.md) |
 | lookback の週末 | 木曜の金利が金曜から月曜までの 3 日に掛かる係数 1 + r·3/365 | [スケジュール](04-foundation/02-schedules.md) |
 | 回帰の随伴 | A = (1, 2)ᵀ、y = (1, 0)ᵀ で β = 1/5、ȳ = (1/5, 2/5)、Ā = (3/25, −4/25) | [XVA のリスク](09-xva/03-xva-risk.md) |
 | CSA の額 | 相手の閾値 50、自分の独立担保 30、V = 40 なら受け渡しは無い | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md) |
@@ -140,7 +149,9 @@ M = max(|x|, |y|, Σ|項|)
 
 - τ_i は quote i の校正の許容誤差です。
   両方の値が独立に解いた校正から来るため、最悪の差は 2 倍になり、c_s = 2 とします。
-- 最小二乗のグループは正しい解でも quote の残差が残るため、quote の許容誤差ではなく、停留条件の残差 G と校正の 2 階の行列 H から θ の誤差 δθ = H⁻¹ G を求めて伝えます。
+- 最小二乗のグループは正しい解でも quote の残差が残るため、quote の許容誤差ではなく、停留条件の残差 G と校正の 2 階の行列 H から θ の誤差 δθ = H⁻¹ G を求めて伝えます（`GroupSolution::theta_error`）。
+  バンプの照合でも同じ式を使い、中心差分が使う各再校正の `theta_error` と、その市場での PV の θ 勾配から Σ_j |∂V/∂θ_j|·|δθ_j| を求め、大きい方を一つの PV の誤差にします（テスト `penalized_meeting_steps_quote_delta_matches_bump`）。
+  校正の停止の判定 ‖Ĵ δθ‖ ≤ τ は各行の動きを抑えますが、罰則だけで決まる向きには対応する市場の quote が無いので、「許容誤差 × quote の感応度の和」ではこの誤差を覆えません。
 - 校正そのものの検証は、完全に決まるグループの quote の再現、最小二乗のグループの停留条件と当てはめの許容幅、同じ θ で校正商品を取引として評価した値とモデルの quote の一致の三つに分けます（[校正](06-market/03-calibration.md)）。
 - バンプの差分商では、二つの PV のソルバーの幅 δ_V から c_s · δ_V / (2h) を使います。
 
@@ -209,11 +220,11 @@ M = max(|x|, |y|, Σ|項|)
 | 記号 | 組み方 | マイルストーン |
 |---|---|---|
 | R1 | 単一の OIS 曲線（正方系の Newton） | M0 |
-| R2 | OIS の割引とターム指標の予想の 2 本の曲線 | M1 |
-| R3 | 会合日のステップと年末のジャンプを持つ OIS | M1 |
-| R4 | 罰則付き最小二乗の過剰決定の曲線 | M1 |
-| R5 | 基準の曲線へのスプレッドの曲線 | M1 |
-| R6 | FX から導く曲線と XCCY のベーシス | M1 |
+| R2 | OIS の割引とターム指標の予想の 2 本の曲線 | M1（`term_products_quote_deltas_match_bumps_on_term_curves`、USD の独立したターム曲線を含む） |
+| R3 | 会合日のステップと年末のジャンプを持つ OIS | M1（`meeting_steps_and_turn_quote_deltas_match_bumps`） |
+| R4 | 罰則付き最小二乗の過剰決定の曲線 | M1（`penalized_meeting_steps_quote_delta_matches_bump`） |
+| R5 | 基準の曲線へのスプレッドの曲線 | M1（`term_products_quote_deltas_match_bumps_on_term_curves`、`tibor_basis_quote_delta_matches_full_recalibration`） |
+| R6 | FX から導く曲線と XCCY のベーシス | M1（`xccy_fx_forward_and_fx_swap_quote_deltas_match_bumps`） |
 | R7 | スワップションのキューブ（ATM とオフセット、SABR） | M2 |
 | R8 | キャップレットのストリッピング | M2 |
 | R9 | LGM と FX のモデルの校正バスケット | M4 |
@@ -232,8 +243,8 @@ M = max(|x|, |y|, Σ|項|)
 | 商品の族 | 組み方 | 量 | マイルストーン |
 |---|---|---|---|
 | 預金、固定対 RFR スワップ（TONA、SOFR の OIS を含む） | R1 | Δq | M0 |
-| OIS、IRS（TIBOR）、ベーシス、FRA、金利先物 | R1〜R5 | Δq（Γ は M3） | M1 |
-| XCCY（MtM のリセット付き）、FX フォワードと FX スワップ | R6 | Δq、Δx | M1 |
+| OIS、IRS（TIBOR）、ベーシス、FRA、金利先物 | R1〜R5 | Δq（Γ は M3） | M1（`term_products_quote_deltas_match_bumps_on_term_curves`、`rate_future_quote_deltas_match_bumps`（LGM のマークを含む）、`meeting_steps_and_turn_quote_deltas_match_bumps`（lookback、observation shift、lockout、支払いの遅れ、年末をまたぐ預金）。最初の検証データセットの全取引の報告通貨のデルタは `m1_rates_portfolio_quote_delta_matches_bump_and_recalibrate`） |
+| XCCY（MtM のリセット付き）、FX フォワードと FX スワップ | R6 | Δq、Δx | M1（`xccy_fx_forward_and_fx_swap_quote_deltas_match_bumps`） |
 | キャップとフロア（ターム金利と後決めの RFR） | R2、R8 | Δq、ν | M2 |
 | スワップション、CMS、CMS スプレッドオプション | R2、R7 | Δq、ν（Γ は M3） | M2 |
 | FX バニラ（ATM） | R6 | Δx、ν | M2 |
@@ -254,17 +265,34 @@ M = max(|x|, |y|, Σ|項|)
 - [金額と量](04-foundation/03-money-and-quantities.md) の混同表の行のうち、型で守る行と、[5 つのブランド](03-architecture.md#ブランド) には、compile-fail テストを書きます。
 - compile-fail テストは必ず、同じ形で正しく書いた `_ok` の例と対にします。
   `_ok` が通ることで、失敗の例が無関係なエラーで落ちていないことを確かめます。
-- コンパイラの版で trybuild の `.stderr` が変わったときは、`.stderr` だけを更新します。
+- 対の登録（`tests/ui` の三つのファイルと `compile_fail_pair!`）は testkit のテストが照合し、登録を通らない書き方（doctest の `compile_fail`、ランナーの直接呼び出し、trybuild の使用）は [CI](#ci) の禁止パターン `P-doctest-compile-fail`、`P-ui-runner`、`P-trybuild` が拒否します。
+- testkit は正常例の終了ステータスを確認し、失敗例では rustc のエラーを全件照合します。`case.json` は各エラーの `code`（コードがなければ `null`）、`message`、すべての主診断位置 `spans`（crate 相対の `file`、1 始まりの `line` と `column`）を配列に記録します。ツールチェーン更新時は実際の診断を確認して期待値を更新します。
 - 実行時のタグで守る行（通貨のタグ、日数基準のタグ）は、エラーを返すことを確かめる単体テストにします。
 
 ## CI
 
 - 検査はリポジトリの `.githooks` のフック（`git config core.hooksPath .githooks`）で行い、外部の CI サービスは使いません。
-  pre-commit は、コミットする Rust のファイルの rustfmt と、TOML・JSON・Python・Nix の構文を検査します。
-  pre-push は、送る commit を一時ディレクトリに取り出し、変更した crate（`Cargo.toml`、`Cargo.lock`、ツールチェーン、`clippy.toml` を変えたときはワークスペース全体）に `cargo clippy --locked --all-features --all-targets -- -D warnings`、`cargo test --locked --all-features`、`RUSTDOCFLAGS="-D warnings" cargo doc --locked --all-features --no-deps` を実行します。
-  作業ツリーが送る commit と同じなら、リポジトリの `target` で増分ビルドします。
-  一時ディレクトリに取り出したときは、ビルド先もその一時ディレクトリの中に分けます。
-  テストは `env!("CARGO_MANIFEST_DIR")` でソースのパスを埋め込むので、一時の取り出しの成果物をリポジトリの `target` に置くと、後の通常の `cargo test` がそれを再利用し、削除済みの一時パスを読んで失敗するためです。
+  pre-commit は、コミットする Rust のファイルの rustfmt と、TOML・JSON・Python・Nix の構文と、下の禁止パターンを検査します。
+  pre-push は、送信する差分から対象を決め、現在の作業ツリーで検査範囲の crate に `cargo clippy --locked --all-features --all-targets -- -D warnings`、`cargo test --locked --all-features`、`RUSTDOCFLAGS="-D warnings" cargo doc --locked --all-features --no-deps` を実行します。
+  検査範囲は、変更した crate と、それに直接・間接に依存する crate（`cargo metadata` の normal・dev・build の依存の逆向きの閉包）です。
+  下流の crate のテストは上流の振る舞いやファイルに依存するので、変更した crate だけを検査すると下流の失敗を見逃すためです。
+  `examples/` を変えたときは、それを読む `rhodium` とその逆依存（CLI）を検査します。
+  `Cargo.toml`、`Cargo.lock`、ツールチェーン、`clippy.toml` を変えたときと、閉包がすべての crate になったときはワークスペース全体を検査します。
+  `.githooks` の Python を変えたときは、Ruff、ty と、`.githooks/test_*.py` の単体テストを実行します。
+  未コミットの変更や別ブランチの push でも、空白と禁止パターンの差分は送信先の基準版と作業ツリーから取り、ビルドには通常の Cargo キャッシュを使います。
+  型テストの登録や `tests/ui/` の追加・変更・削除があれば、testkit の登録検査も実行します。
+- 型や lint で止められない実装の誤りの形は、`.githooks/checks.py` の禁止パターンの表（`RULES`）で止めます。
+  各行は、名前（`P-…`）、Rust の行に当てる正規表現、対象と例外のパスの glob、守る理由、規則の定義元の設計書の節を持ちます。
+  ライブラリの規則はテストのファイル（`tests/`、`src/**/tests.rs`、`src/**/tests/`）と `#[cfg(test)]` の item を照合せず、コメントの中も照合しません。
+  `requires_nearby` を持つ行は、当たった行とその直前の数行に決まった注記（`P-internal` の `// 不変条件:`）があれば違反にしません。
+  行のモードは二つです。
+  `diff` の行は、commit と push で追加した行だけを照合するので、既存の違反を残したまま規則を入れられます。
+  `tree` の行は、送る版のツリー全体を照合し、既存の違反も許しません。
+  規則は、既存の違反があっても `diff` で先に入れ、担当の issue が既存の違反をすべてなくした PR で `tree` に切り替えます。
+  切り替えた後は同じ誤りを書き戻せません。
+  `python3 .githooks/checks.py audit` は、作業ツリーの全ファイルを全規則で照合し、規則ごとの違反の件数と場所を表示します（失敗にはしません）。
+  規則を足すときは、まず `audit` で既存の違反を数え、`diff` で入れる理由か、`tree` で入れられることを PR に書きます。
+  違反に当たったときは、検査を避ける書き換えではなく、規則の定義元の設計書に従って直します。
 - マイルストーンの完了確認では、ワークスペース全体で `cargo fmt --all --check`、`cargo clippy --locked --workspace --all-targets --all-features -- -D warnings`、`cargo test --locked --workspace --all-features`、`RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --all-features --no-deps` を実行し、結果を完了確認の issue に記録します。
 - 大きな MC、外部 fixture の全件、XVA の規模のテストは `#[ignore = "理由"]` を付けた遅いテストにし、`cargo test --workspace --all-features -- --ignored` で定期的に実行します。
 - ベンチマークはテストと分け、[入出力と実行](10-io-and-execution.md) の手順で測ります。
diff --git a/docs/design/12-roadmap.md b/docs/design/12-roadmap.md
index 2bbcee32989046838867575844a2edc352f1cb83..ac78029a209ca00908b44d3ca498a345788cb34d 100644
--- a/docs/design/12-roadmap.md
+++ b/docs/design/12-roadmap.md
@@ -35,7 +35,9 @@ M3 と M4 は M2 の後に並行して進められます。
   例外は 1.0 の判定の M6-4（公開 API とファイル形式の見直し）と M6-5（性能の目安の測定）で、この二つは手続きを含みます。
 - 各基準のテストは [品質](11-quality.md) のオラクルと許容誤差に従います。
   テストの名前は説明的な snake_case とし、基準の ID はテストのコードではなく、この文書、issue、PR に書きます。
-- 「quote を再現する」は、完全に決まる校正グループではソルバーの許容誤差での再現を、最小二乗のグループでは停留条件と当てはめの許容幅を満たすことを指します（[校正の検証](06-market/03-calibration.md#校正の検証)）。
+- 「quote を再現する」は、完全に決まる校正グループではソルバーの許容誤差での再現を指します。
+  最小二乗のグループでは、停留条件を満たし、残差が許容幅を超えた quote を警告 `CalibrationFitOutsideTolerance` で報告することを指します（許容幅の超過は警告付きの成功で、失敗ではありません。[校正の検証](06-market/03-calibration.md#校正の検証)）。
+- 「校正商品を取引として評価するとパーになる」も同じく、完全に決まるグループでは取引のパーが quote に一致することを、最小二乗のグループでは取引のパーが quote と残差の和 q + F に一致することを指します（テスト `penalized_meeting_steps_fit_standard_tenor_ois`、`fx_outright_least_squares_measures_residuals_in_price_units`）。
 - 性能は各マイルストーンで所要時間を記録し、目安（[入出力と実行](10-io-and-execution.md#性能の目安)）と比べます。
   目安は判定の門にしませんが、所要時間やメモリが目安を大きく超えるときは原因を issue に記録します。
 
@@ -165,7 +167,7 @@ XVA を MVA と KVA まで広げ、取引前の評価を速くし、1.0 とし
 | M6-1 | 動的 SIMM の時点 0 の値が SIMM と一致し、MVA が手計算の例（V = 2qz）と一致する |
 | M6-2 | SA-CCR の EAD が BCBS の例と一致し、BA-CVA と KVA が手計算の例と一致する |
 | M6-3 | 候補取引の incremental XVA と XVA 込みのゴールシークが動き、所要時間を記録する |
-| M6-4 | ファサードの公開 API とファイル形式を見直し、1.0 以降の互換性の方針（SemVer、ファイルの `version` と移行）を文書にする |
+| M6-4 | ファサードの公開 API とファイル形式を見直し、1.0 以降の互換性の方針（SemVer、ファイルの `version` と移行）を文書にする。M1 までに後へ送った形式（会合表・年末の規則・テナーによる日付の解決、`tenors` の雛形の展開、JSON の代替、quote の TOML と bid/ask）を決める |
 | M6-5 | 性能の目安の全項目を測り、報告を残す |
 
 ## 1.0 の後のマイルストーン
diff --git a/docs/design/13-glossary.md b/docs/design/13-glossary.md
index 179d24c20a18e6a97bee2bf638966e0dc130d8db..6d9b6c9cf496356672a653706a23151941010d52 100644
--- a/docs/design/13-glossary.md
+++ b/docs/design/13-glossary.md
@@ -92,6 +92,7 @@
 | FX の換算の担保 | `fx_collateral` | 担保の無い FX フォワードに使う担保。報告通貨への換算は通貨ペアごとに一つの X(0) を使い、担保によらない | [07-1](07-pricing/01-pricing.md) |
 | 要求 | requirement | 取引が必要とする市場の対象（この通貨と担保の割引、この指標の予想） | [07-1](07-pricing/01-pricing.md) |
 | 要求の解決 | resolution | 要求を価格設定の表で曲線などに結ぶこと | [07-1](07-pricing/01-pricing.md) |
+| 予想の測度 | forecast measure | 予想曲線の予想を定める担保の割引曲線。曲線を当てる校正商品の担保から校正が導く | [07-2](07-pricing/02-linear-and-analytic.md#担保の違いの判定) |
 | 評価方式 | method | リニア、解析、モンテカルロ、AMC のどれで評価するか | [07-1](07-pricing/01-pricing.md) |
 | 準備 | Prepared | θ に依存しない評価の準備。quote やシナリオの間で使い回す | [07-1](07-pricing/01-pricing.md) |
 | リニア評価 | linear evaluation | 線形な IR を DF と forward の積和で評価すること | [07-2](07-pricing/02-linear-and-analytic.md) |
diff --git a/examples/m1-rates/README.md b/examples/m1-rates/README.md
new file mode 100644
index 0000000000000000000000000000000000000000..434ff3326625a1d5b856fa1a5110b61ca9f0c41c
--- /dev/null
+++ b/examples/m1-rates/README.md
@@ -0,0 +1,42 @@
+# M1 の金利リニアのワークスペース（最初の検証データセット）
+
+JPY（TONA、TIBOR 6M と 3M）と USD（SOFR）の曲線、USDJPY の FX と外貨担保曲線を一つの市場の定義で校正し、M1 のリニア商品と CSA の 4 つの場合の取引を評価します。
+評価時点は 2026 年 11 月 24 日の東京 17 時（UTC 08 時）です。
+quote、fixing、休日、取引の値は検証のためのサンプルで、2028 年以後の休日は週末だけです。
+
+リポジトリのルートから実行します。
+
+```sh
+cargo run --locked -p rhodium-cli -- price examples/m1-rates --date 2026-11-24
+cargo run --locked -p rhodium-cli -- risk examples/m1-rates --date 2026-11-24 --threads 4
+```
+
+## 市場の定義
+
+| 曲線 | 形 | 校正 |
+|---|---|---|
+| `JPY-TONA` | 1 年までの 3 回の会合で区切った 4 区間の階段（`steps`）を 1 年の点で長期の節点（`nodes`）に接続し、2026 年末の跳び（`jumps`）を足す | 階段は 1M、3M、6M、1Y の OIS、長期は 2Y〜10Y の OIS、跳びは年末のスプレッド 50bp × 5/365 の直接マーク |
+| `JPY-TIBOR-6M`、`JPY-TIBOR-3M` | `JPY-TONA` と `basis` の節点の和 | TONA 担保の TIBOR 対 TONA のベーシススワップ（1Y、2Y、5Y、10Y） |
+| `USD-SOFR` | 明示した 7 つのピラーの節点 | USD 担保の SOFR の OIS（支払いは 2 営業日後） |
+| `USD-under-JPY` | `[[fx]]` の外貨担保曲線（`FxImplied`） | 短期は USDJPY の TN と 3M・6M のポイント、長期は MtM XCCY のベーシス（1Y、2Y、5Y） |
+| `JPY-TONA-25` | `JPY-TONA` と一定の 25bp の和 | 校正しない（25bp の付利の CSA の割引） |
+| `USD-CTD` | `USD-SOFR` と `USD-under-JPY` の毎年の平均フォワードの大きい方（`max_forward`） | 校正しない（JPY と USD の現金担保の CTD の割引） |
+
+- USDJPY の受渡慣行（T+2、東京 15 時、pip は 0.01 円）は `refdata/conventions.toml` の `[[fx]]` にあります。
+  評価の 17 時は当日の受渡の締めを過ぎているので、ON のポイントは連鎖に使わず、X(0) は翌受渡（TN）からスポット、ポイントの順に連鎖させて割り戻します。
+- 校正のグループは、年末の跳び、TONA の階段、TONA の長期、SOFR、FX のマーク、TIBOR 6M・3M、外貨担保曲線（FX ポイントと XCCY の一つのグループ）、LGM のマークです。
+  TONA の 10Y の quote を変えると、TONA の長期と、その下流の TIBOR 6M・3M と外貨担保曲線だけを解き直します。
+
+## 取引
+
+| ファイル | 取引 | 確かめること |
+|---|---|---|
+| `trades/rates.toml` | 3M の預金、TIBOR 6M の 6×12 FRA、2027 年 3 月限の TONA 3M 先物 10 枚、5Y の固定対 TIBOR 6M の IRS、5Y の TIBOR 3M 対 TONA のベーシス、lookback・observation shift・lockout の 2Y OIS（支払いは 2 営業日後）、開始済みの IRS を評価時点から 4 割に一部解約し手数料を受け取る取引 | 先物は解析プライサー（LGM の凸性の補正）で評価する。FRA は `FraDiscounting`、lookback と lockout は `DailyForwardProduct` を記録する。一部解約の PV は元の IRS の 0.4 倍と手数料の PV の和 |
+| `trades/fx.toml` | 1Y の FX フォワード、スポット対 6M の FX スワップ | JPY 現金担保の下で USD を外貨担保曲線で割り引く |
+| `trades/xccy.toml` | 2Y の MtM XCCY（USD の元本を USDJPY の東京 15 時の fixing でリセット）と固定元本の XCCY | MtM は `MtmFxForwardProduct`、どちらも SOFR の予想と JPY 担保の割引の担保が違うので `CrossCollateralProjection` を記録する |
+| `trades/csa.toml` | 同じ 2Y の SOFR の OIS を JPY 現金担保（外貨担保）、USD 現金担保（自国通貨）、CTD、無担保で、5Y の TIBOR 6M の IRS を 25bp の付利の CSA で | CSA の 4 つの場合の割引。TONA 担保の TIBOR の IRS は `CrossCollateralProjection` を記録せず、25bp の CSA の IRS と JPY 担保の SOFR の OIS だけが記録する |
+
+## 後のマイルストーンに送ったもの
+
+会合表（`meetings = "BOJ"`）や年末の規則を参照データから解決する書き方、`tenors` による行の展開は、ファイル形式の見直し（M6）で足します。
+このワークスペースは、解決済みの UTC の瞬間と日数係数を直接書いています（[校正の設計書](../../docs/design/06-market/03-calibration.md#ワークスペースの市場定義)）。
diff --git a/examples/m1-rates/agreements.toml b/examples/m1-rates/agreements.toml
new file mode 100644
index 0000000000000000000000000000000000000000..e537b53d3be54cfee784f2354bcf91890998b86e
--- /dev/null
+++ b/examples/m1-rates/agreements.toml
@@ -0,0 +1,21 @@
+version = 1
+
+[[margin_set]]
+id = "JPY-CASH"
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+
+[[margin_set]]
+id = "JPY-CASH-25"
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "25bp" } }
+
+[[margin_set]]
+id = "USD-CASH"
+collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
+
+[[margin_set]]
+id = "CTD-JPY-USD"
+collateral = { cheapest_to_deliver = [{ ccy = "JPY", index = "JPY-TONA", spread = "0bp" }, { ccy = "USD", index = "USD-SOFR", spread = "0bp" }] }
+
+[[margin_set]]
+id = "UNCOLLATERALIZED"
+collateral = "uncollateralized"
diff --git a/examples/m1-rates/config/pricing.toml b/examples/m1-rates/config/pricing.toml
new file mode 100644
index 0000000000000000000000000000000000000000..3825fbd4b5d3f73749aa692a40939577f12d9338
--- /dev/null
+++ b/examples/m1-rates/config/pricing.toml
@@ -0,0 +1,44 @@
+version = 1
+reporting_currency = "JPY"
+valuation_close = { time = "17:00", zone = "Asia/Tokyo" }
+
+[[discounting]]
+currency = "JPY"
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+curve = "JPY-TONA"
+
+[[discounting]]
+currency = "USD"
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+curve = "USD-under-JPY"
+
+[[discounting]]
+currency = "JPY"
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "25bp" } }
+curve = "JPY-TONA-25"
+
+[[discounting]]
+currency = "USD"
+collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
+curve = "USD-SOFR"
+
+[[discounting]]
+currency = "USD"
+collateral = { cheapest_to_deliver = [{ ccy = "JPY", index = "JPY-TONA", spread = "0bp" }, { ccy = "USD", index = "USD-SOFR", spread = "0bp" }] }
+curve = "USD-CTD"
+
+[[discounting]]
+currency = "JPY"
+collateral = "uncollateralized"
+curve = "JPY-TONA"
+
+[[discounting]]
+currency = "USD"
+collateral = "uncollateralized"
+curve = "USD-SOFR"
+
+[forecasting]
+JPY-TONA = "JPY-TONA"
+USD-SOFR = "USD-SOFR"
+JPY-TIBOR-6M = "JPY-TIBOR-6M"
+JPY-TIBOR-3M = "JPY-TIBOR-3M"
diff --git a/examples/m1-rates/fixings/JPY-TIBOR-3M.csv b/examples/m1-rates/fixings/JPY-TIBOR-3M.csv
new file mode 100644
index 0000000000000000000000000000000000000000..883a1934a052f1f2b327608622e9b055235cc2d7
--- /dev/null
+++ b/examples/m1-rates/fixings/JPY-TIBOR-3M.csv
@@ -0,0 +1,2 @@
+date,value
+2026-11-24,0.68%
diff --git a/examples/m1-rates/fixings/JPY-TIBOR-6M.csv b/examples/m1-rates/fixings/JPY-TIBOR-6M.csv
new file mode 100644
index 0000000000000000000000000000000000000000..af76af36f4bb7aa1eba971939223bb1b0758ff5c
--- /dev/null
+++ b/examples/m1-rates/fixings/JPY-TIBOR-6M.csv
@@ -0,0 +1,3 @@
+date,value
+2026-05-28,0.70%
+2026-11-24,0.80%
diff --git a/examples/m1-rates/fixings/JPY-TONA.csv b/examples/m1-rates/fixings/JPY-TONA.csv
new file mode 100644
index 0000000000000000000000000000000000000000..702e2aeb879f36c174b65cb6c3fbbbac19b2297c
--- /dev/null
+++ b/examples/m1-rates/fixings/JPY-TONA.csv
@@ -0,0 +1,146 @@
+date,value
+2026-05-01,0.477%
+2026-05-04,0.477%
+2026-05-05,0.477%
+2026-05-06,0.477%
+2026-05-07,0.477%
+2026-05-08,0.477%
+2026-05-11,0.477%
+2026-05-12,0.477%
+2026-05-13,0.477%
+2026-05-14,0.477%
+2026-05-15,0.477%
+2026-05-18,0.477%
+2026-05-19,0.477%
+2026-05-20,0.477%
+2026-05-21,0.477%
+2026-05-22,0.477%
+2026-05-25,0.477%
+2026-05-26,0.477%
+2026-05-27,0.477%
+2026-05-28,0.477%
+2026-05-29,0.477%
+2026-06-01,0.477%
+2026-06-02,0.477%
+2026-06-03,0.477%
+2026-06-04,0.477%
+2026-06-05,0.477%
+2026-06-08,0.477%
+2026-06-09,0.477%
+2026-06-10,0.477%
+2026-06-11,0.477%
+2026-06-12,0.477%
+2026-06-15,0.477%
+2026-06-16,0.477%
+2026-06-17,0.477%
+2026-06-18,0.477%
+2026-06-19,0.477%
+2026-06-22,0.477%
+2026-06-23,0.477%
+2026-06-24,0.477%
+2026-06-25,0.477%
+2026-06-26,0.477%
+2026-06-29,0.477%
+2026-06-30,0.477%
+2026-07-01,0.477%
+2026-07-02,0.477%
+2026-07-03,0.477%
+2026-07-06,0.477%
+2026-07-07,0.477%
+2026-07-08,0.477%
+2026-07-09,0.477%
+2026-07-10,0.477%
+2026-07-13,0.477%
+2026-07-14,0.477%
+2026-07-15,0.477%
+2026-07-16,0.477%
+2026-07-17,0.477%
+2026-07-20,0.477%
+2026-07-21,0.477%
+2026-07-22,0.477%
+2026-07-23,0.477%
+2026-07-24,0.477%
+2026-07-27,0.477%
+2026-07-28,0.477%
+2026-07-29,0.477%
+2026-07-30,0.477%
+2026-07-31,0.477%
+2026-08-03,0.477%
+2026-08-04,0.477%
+2026-08-05,0.477%
+2026-08-06,0.477%
+2026-08-07,0.477%
+2026-08-10,0.477%
+2026-08-11,0.477%
+2026-08-12,0.477%
+2026-08-13,0.477%
+2026-08-14,0.477%
+2026-08-17,0.477%
+2026-08-18,0.477%
+2026-08-19,0.477%
+2026-08-20,0.477%
+2026-08-21,0.477%
+2026-08-24,0.477%
+2026-08-25,0.477%
+2026-08-26,0.477%
+2026-08-27,0.477%
+2026-08-28,0.477%
+2026-08-31,0.477%
+2026-09-01,0.477%
+2026-09-02,0.477%
+2026-09-03,0.477%
+2026-09-04,0.477%
+2026-09-07,0.477%
+2026-09-08,0.477%
+2026-09-09,0.477%
+2026-09-10,0.477%
+2026-09-11,0.477%
+2026-09-14,0.477%
+2026-09-15,0.477%
+2026-09-16,0.477%
+2026-09-17,0.477%
+2026-09-18,0.477%
+2026-09-21,0.477%
+2026-09-22,0.477%
+2026-09-23,0.477%
+2026-09-24,0.477%
+2026-09-25,0.477%
+2026-09-28,0.477%
+2026-09-29,0.477%
+2026-09-30,0.477%
+2026-10-01,0.477%
+2026-10-02,0.477%
+2026-10-05,0.477%
+2026-10-06,0.477%
+2026-10-07,0.477%
+2026-10-08,0.477%
+2026-10-09,0.477%
+2026-10-12,0.477%
+2026-10-13,0.477%
+2026-10-14,0.477%
+2026-10-15,0.477%
+2026-10-16,0.477%
+2026-10-19,0.477%
+2026-10-20,0.477%
+2026-10-21,0.477%
+2026-10-22,0.477%
+2026-10-23,0.477%
+2026-10-26,0.477%
+2026-10-27,0.477%
+2026-10-28,0.477%
+2026-10-29,0.477%
+2026-10-30,0.477%
+2026-11-02,0.477%
+2026-11-04,0.477%
+2026-11-05,0.477%
+2026-11-06,0.477%
+2026-11-09,0.477%
+2026-11-10,0.477%
+2026-11-11,0.477%
+2026-11-12,0.477%
+2026-11-13,0.477%
+2026-11-16,0.477%
+2026-11-17,0.477%
+2026-11-18,0.477%
+2026-11-19,0.477%
+2026-11-20,0.477%
diff --git a/examples/m1-rates/fixings/JPY-TONA3M-FUT-2027-03.csv b/examples/m1-rates/fixings/JPY-TONA3M-FUT-2027-03.csv
new file mode 100644
index 0000000000000000000000000000000000000000..da8ea1980fef11c360961c1f2bafe8b89fb51a65
--- /dev/null
+++ b/examples/m1-rates/fixings/JPY-TONA3M-FUT-2027-03.csv
@@ -0,0 +1,2 @@
+date,value
+2026-11-24,0.9935
diff --git a/examples/m1-rates/fixings/USD-SOFR.csv b/examples/m1-rates/fixings/USD-SOFR.csv
new file mode 100644
index 0000000000000000000000000000000000000000..bd9c3b15cf0e1b21b1f876536ce1091bb713f72a
--- /dev/null
+++ b/examples/m1-rates/fixings/USD-SOFR.csv
@@ -0,0 +1,146 @@
+date,value
+2026-05-01,3.75%
+2026-05-04,3.75%
+2026-05-05,3.75%
+2026-05-06,3.75%
+2026-05-07,3.75%
+2026-05-08,3.75%
+2026-05-11,3.75%
+2026-05-12,3.75%
+2026-05-13,3.75%
+2026-05-14,3.75%
+2026-05-15,3.75%
+2026-05-18,3.75%
+2026-05-19,3.75%
+2026-05-20,3.75%
+2026-05-21,3.75%
+2026-05-22,3.75%
+2026-05-25,3.75%
+2026-05-26,3.75%
+2026-05-27,3.75%
+2026-05-28,3.75%
+2026-05-29,3.75%
+2026-06-01,3.75%
+2026-06-02,3.75%
+2026-06-03,3.75%
+2026-06-04,3.75%
+2026-06-05,3.75%
+2026-06-08,3.75%
+2026-06-09,3.75%
+2026-06-10,3.75%
+2026-06-11,3.75%
+2026-06-12,3.75%
+2026-06-15,3.75%
+2026-06-16,3.75%
+2026-06-17,3.75%
+2026-06-18,3.75%
+2026-06-19,3.75%
+2026-06-22,3.75%
+2026-06-23,3.75%
+2026-06-24,3.75%
+2026-06-25,3.75%
+2026-06-26,3.75%
+2026-06-29,3.75%
+2026-06-30,3.75%
+2026-07-01,3.75%
+2026-07-02,3.75%
+2026-07-03,3.75%
+2026-07-06,3.75%
+2026-07-07,3.75%
+2026-07-08,3.75%
+2026-07-09,3.75%
+2026-07-10,3.75%
+2026-07-13,3.75%
+2026-07-14,3.75%
+2026-07-15,3.75%
+2026-07-16,3.75%
+2026-07-17,3.75%
+2026-07-20,3.75%
+2026-07-21,3.75%
+2026-07-22,3.75%
+2026-07-23,3.75%
+2026-07-24,3.75%
+2026-07-27,3.75%
+2026-07-28,3.75%
+2026-07-29,3.75%
+2026-07-30,3.75%
+2026-07-31,3.75%
+2026-08-03,3.75%
+2026-08-04,3.75%
+2026-08-05,3.75%
+2026-08-06,3.75%
+2026-08-07,3.75%
+2026-08-10,3.75%
+2026-08-11,3.75%
+2026-08-12,3.75%
+2026-08-13,3.75%
+2026-08-14,3.75%
+2026-08-17,3.75%
+2026-08-18,3.75%
+2026-08-19,3.75%
+2026-08-20,3.75%
+2026-08-21,3.75%
+2026-08-24,3.75%
+2026-08-25,3.75%
+2026-08-26,3.75%
+2026-08-27,3.75%
+2026-08-28,3.75%
+2026-08-31,3.75%
+2026-09-01,3.75%
+2026-09-02,3.75%
+2026-09-03,3.75%
+2026-09-04,3.75%
+2026-09-07,3.75%
+2026-09-08,3.75%
+2026-09-09,3.75%
+2026-09-10,3.75%
+2026-09-11,3.75%
+2026-09-14,3.75%
+2026-09-15,3.75%
+2026-09-16,3.75%
+2026-09-17,3.75%
+2026-09-18,3.75%
+2026-09-21,3.75%
+2026-09-22,3.75%
+2026-09-23,3.75%
+2026-09-24,3.75%
+2026-09-25,3.75%
+2026-09-28,3.75%
+2026-09-29,3.75%
+2026-09-30,3.75%
+2026-10-01,3.75%
+2026-10-02,3.75%
+2026-10-05,3.75%
+2026-10-06,3.75%
+2026-10-07,3.75%
+2026-10-08,3.75%
+2026-10-09,3.75%
+2026-10-12,3.75%
+2026-10-13,3.75%
+2026-10-14,3.75%
+2026-10-15,3.75%
+2026-10-16,3.75%
+2026-10-19,3.75%
+2026-10-20,3.75%
+2026-10-21,3.75%
+2026-10-22,3.75%
+2026-10-23,3.75%
+2026-10-26,3.75%
+2026-10-27,3.75%
+2026-10-28,3.75%
+2026-10-29,3.75%
+2026-10-30,3.75%
+2026-11-02,3.75%
+2026-11-03,3.75%
+2026-11-04,3.75%
+2026-11-05,3.75%
+2026-11-06,3.75%
+2026-11-09,3.75%
+2026-11-10,3.75%
+2026-11-12,3.75%
+2026-11-13,3.75%
+2026-11-16,3.75%
+2026-11-17,3.75%
+2026-11-18,3.75%
+2026-11-19,3.75%
+2026-11-20,3.75%
diff --git a/examples/m1-rates/fixings/USDJPY-TKY-1500.csv b/examples/m1-rates/fixings/USDJPY-TKY-1500.csv
new file mode 100644
index 0000000000000000000000000000000000000000..d2fa9a9944c5241923b7f4e7c15c925e50e3d85c
--- /dev/null
+++ b/examples/m1-rates/fixings/USDJPY-TKY-1500.csv
@@ -0,0 +1,2 @@
+date,value
+2026-11-24,152.35
diff --git a/examples/m1-rates/market/definition.toml b/examples/m1-rates/market/definition.toml
new file mode 100644
index 0000000000000000000000000000000000000000..faaa1ba9aeab2eed536690eedd06cb5abfbf91a2
--- /dev/null
+++ b/examples/m1-rates/market/definition.toml
@@ -0,0 +1,267 @@
+version = 1
+valuation_calendar = "JPTO"
+valuation_zone = "Asia/Tokyo"
+
+# TONA: 1 年までの 3 回の会合で区切った 4 区間の階段、1 年で接続する長期の節点、2026 年末の跳び。
+[curves.JPY-TONA]
+currency = "JPY"
+roles = ["discount", "projection"]
+
+[[curves.JPY-TONA.expr.sum]]
+stitch.at = 2027-12-03T00:00:00Z
+stitch.short.steps = { part = "meetings", meetings = [2026-12-18T00:00:00Z, 2027-01-22T00:00:00Z, 2027-03-18T00:00:00Z] }
+stitch.long.nodes = { part = "long", pillars = [2028-12-04T00:00:00Z, 2029-12-03T00:00:00Z, 2031-12-03T00:00:00Z, 2033-12-05T00:00:00Z, 2036-12-03T00:00:00Z] }
+
+[[curves.JPY-TONA.expr.sum]]
+jumps = { part = "turns", dates = [2027-01-03T15:00:00Z] }
+
+# TIBOR 6M: TONA の曲線に、TONA 担保のベーシスで決める節点を足す。
+[curves.JPY-TIBOR-6M]
+currency = "JPY"
+roles = ["projection"]
+expr.sum = [{ ref = "JPY-TONA" }, { nodes = { part = "basis", pillars = [2027-12-03T00:00:00Z, 2028-12-04T00:00:00Z, 2031-12-03T00:00:00Z, 2036-12-03T00:00:00Z] } }]
+
+# TIBOR 3M: TONA の曲線に、TONA 担保のベーシスで決める節点を足す。
+[curves.JPY-TIBOR-3M]
+currency = "JPY"
+roles = ["projection"]
+expr.sum = [{ ref = "JPY-TONA" }, { nodes = { part = "basis", pillars = [2027-12-03T00:00:00Z, 2028-12-04T00:00:00Z, 2031-12-03T00:00:00Z, 2036-12-03T00:00:00Z] } }]
+
+[curves.USD-SOFR]
+currency = "USD"
+roles = ["discount", "projection"]
+pillars = [2027-01-10T00:00:00Z, 2027-03-13T00:00:00Z, 2027-06-10T00:00:00Z, 2027-12-11T00:00:00Z, 2028-12-11T00:00:00Z, 2031-12-11T00:00:00Z, 2036-12-11T00:00:00Z]
+
+# 25bp の付利の JPY 現金担保で割り引く曲線。
+[curves.JPY-TONA-25]
+currency = "JPY"
+roles = ["discount"]
+expr.sum = [{ ref = "JPY-TONA" }, { flat = "25bp" }]
+
+# USD の支払いの、JPY と USD の現金担保の CTD（各年の平均フォワードの大きい方）。
+[curves.USD-CTD]
+currency = "USD"
+roles = ["discount"]
+expr.max_forward = { parts = [{ ref = "USD-SOFR" }, { ref = "USD-under-JPY" }], grid = [2027-11-27T00:00:00Z, 2028-11-27T00:00:00Z, 2029-11-27T00:00:00Z, 2030-11-27T00:00:00Z, 2031-11-27T00:00:00Z, 2032-11-27T00:00:00Z, 2033-11-27T00:00:00Z, 2034-11-27T00:00:00Z, 2035-11-27T00:00:00Z, 2036-11-27T00:00:00Z] }
+
+[[instruments]]
+quote = "JPY.TONA.OIS.1M"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "JPY-TONA-OIS", tenor = "1M" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TONA", part = "meetings" }
+
+[[instruments]]
+quote = "JPY.TONA.OIS.3M"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "JPY-TONA-OIS", tenor = "3M" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TONA", part = "meetings" }
+
+[[instruments]]
+quote = "JPY.TONA.OIS.6M"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "JPY-TONA-OIS", tenor = "6M" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TONA", part = "meetings" }
+
+[[instruments]]
+quote = "JPY.TONA.OIS.1Y"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "JPY-TONA-OIS", tenor = "1Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TONA", part = "meetings" }
+
+[[instruments]]
+quote = "JPY.TONA.OIS.2Y"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "JPY-TONA-OIS", tenor = "2Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TONA", part = "long" }
+
+[[instruments]]
+quote = "JPY.TONA.OIS.3Y"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "JPY-TONA-OIS", tenor = "3Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TONA", part = "long" }
+
+[[instruments]]
+quote = "JPY.TONA.OIS.5Y"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "JPY-TONA-OIS", tenor = "5Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TONA", part = "long" }
+
+[[instruments]]
+quote = "JPY.TONA.OIS.7Y"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "JPY-TONA-OIS", tenor = "7Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TONA", part = "long" }
+
+[[instruments]]
+quote = "JPY.TONA.OIS.10Y"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "JPY-TONA-OIS", tenor = "10Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TONA", part = "long" }
+
+[[instruments]]
+quote = "JPY.TONA.TURN.2026"
+kind = "spread"
+target = "direct_mark"
+determines = { curve = "JPY-TONA", part = "turns", element = 0 }
+transform = { turn_spread_to_jump = { value = 0.0136986301369863, basis = "ACT365F" } }
+
+[[instruments]]
+quote = "JPY.TIBOR6M-TONA.BASIS.1Y"
+kind = "spread"
+target = "par_rate"
+template = { basis_swap = { convention = "JPY-TIBOR6M-TONA-BASIS", tenor = "1Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TIBOR-6M", part = "basis" }
+
+[[instruments]]
+quote = "JPY.TIBOR6M-TONA.BASIS.2Y"
+kind = "spread"
+target = "par_rate"
+template = { basis_swap = { convention = "JPY-TIBOR6M-TONA-BASIS", tenor = "2Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TIBOR-6M", part = "basis" }
+
+[[instruments]]
+quote = "JPY.TIBOR6M-TONA.BASIS.5Y"
+kind = "spread"
+target = "par_rate"
+template = { basis_swap = { convention = "JPY-TIBOR6M-TONA-BASIS", tenor = "5Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TIBOR-6M", part = "basis" }
+
+[[instruments]]
+quote = "JPY.TIBOR6M-TONA.BASIS.10Y"
+kind = "spread"
+target = "par_rate"
+template = { basis_swap = { convention = "JPY-TIBOR6M-TONA-BASIS", tenor = "10Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TIBOR-6M", part = "basis" }
+
+[[instruments]]
+quote = "JPY.TIBOR3M-TONA.BASIS.1Y"
+kind = "spread"
+target = "par_rate"
+template = { basis_swap = { convention = "JPY-TIBOR3M-TONA-BASIS", tenor = "1Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TIBOR-3M", part = "basis" }
+
+[[instruments]]
+quote = "JPY.TIBOR3M-TONA.BASIS.2Y"
+kind = "spread"
+target = "par_rate"
+template = { basis_swap = { convention = "JPY-TIBOR3M-TONA-BASIS", tenor = "2Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TIBOR-3M", part = "basis" }
+
+[[instruments]]
+quote = "JPY.TIBOR3M-TONA.BASIS.5Y"
+kind = "spread"
+target = "par_rate"
+template = { basis_swap = { convention = "JPY-TIBOR3M-TONA-BASIS", tenor = "5Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TIBOR-3M", part = "basis" }
+
+[[instruments]]
+quote = "JPY.TIBOR3M-TONA.BASIS.10Y"
+kind = "spread"
+target = "par_rate"
+template = { basis_swap = { convention = "JPY-TIBOR3M-TONA-BASIS", tenor = "10Y" } }
+collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
+determines = { curve = "JPY-TIBOR-3M", part = "basis" }
+
+[[instruments]]
+quote = "USD.SOFR.OIS.1M"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "USD-SOFR-OIS", tenor = "1M" } }
+collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
+determines = { curve = "USD-SOFR", element = 0 }
+
+[[instruments]]
+quote = "USD.SOFR.OIS.3M"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "USD-SOFR-OIS", tenor = "3M" } }
+collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
+determines = { curve = "USD-SOFR", element = 1 }
+
+[[instruments]]
+quote = "USD.SOFR.OIS.6M"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "USD-SOFR-OIS", tenor = "6M" } }
+collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
+determines = { curve = "USD-SOFR", element = 2 }
+
+[[instruments]]
+quote = "USD.SOFR.OIS.1Y"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "USD-SOFR-OIS", tenor = "1Y" } }
+collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
+determines = { curve = "USD-SOFR", element = 3 }
+
+[[instruments]]
+quote = "USD.SOFR.OIS.2Y"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "USD-SOFR-OIS", tenor = "2Y" } }
+collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
+determines = { curve = "USD-SOFR", element = 4 }
+
+[[instruments]]
+quote = "USD.SOFR.OIS.5Y"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "USD-SOFR-OIS", tenor = "5Y" } }
+collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
+determines = { curve = "USD-SOFR", element = 5 }
+
+[[instruments]]
+quote = "USD.SOFR.OIS.10Y"
+kind = "rate"
+target = "par_rate"
+template = { ois = { convention = "USD-SOFR-OIS", tenor = "10Y" } }
+collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
+determines = { curve = "USD-SOFR", element = 6 }
+
+# USDJPY: スポット、TN、3M と 6M のポイントと、外貨担保曲線の長期の MtM XCCY。
+# 評価の 17 時は当日の受渡の締め（15 時）を過ぎているので、ON は連鎖に使わない。
+[[fx]]
+pair = "USDJPY"
+collaterals = [{ id = "market", base = "USD-under-JPY", quote = "JPY-TONA" }]
+
+[fx.direct]
+market_collateral = "market"
+quotes = { spot = "USDJPY.SPOT", on = "USDJPY.ON", tn = "USDJPY.TN", forwards = [{ settlement = 2027-02-26, quote = "USDJPY.3M" }, { settlement = 2027-05-27, quote = "USDJPY.6M" }] }
+
+[fx.direct.implied]
+id = "USD-under-JPY"
+long_end = [
+    { pillar = 2027-12-11T00:00:00Z, quote = "USDJPY.XCCY.1Y", convention = "USDJPY-MTM-XCCY", tenor = "1Y", collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } } },
+    { pillar = 2028-12-11T00:00:00Z, quote = "USDJPY.XCCY.2Y", convention = "USDJPY-MTM-XCCY", tenor = "2Y", collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } } },
+    { pillar = 2031-12-11T00:00:00Z, quote = "USDJPY.XCCY.5Y", convention = "USDJPY-MTM-XCCY", tenor = "5Y", collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } } },
+]
+
+# 金利先物の凸性の補正に使う JPY の LGM のマーク。
+[lgm.JPY]
+mean_reversion = ["JPY.LGM.MR"]
+volatility = ["JPY.LGM.VOL"]
diff --git a/examples/m1-rates/quotes/snapshot.csv b/examples/m1-rates/quotes/snapshot.csv
new file mode 100644
index 0000000000000000000000000000000000000000..b5303d6185dcab50cf64259dce673c20a52c99f5
--- /dev/null
+++ b/examples/m1-rates/quotes/snapshot.csv
@@ -0,0 +1,36 @@
+quote_id,value
+JPY.TONA.OIS.1M,0.48%
+JPY.TONA.OIS.3M,0.52%
+JPY.TONA.OIS.6M,0.58%
+JPY.TONA.OIS.1Y,0.65%
+JPY.TONA.OIS.2Y,0.78%
+JPY.TONA.OIS.3Y,0.88%
+JPY.TONA.OIS.5Y,1.02%
+JPY.TONA.OIS.7Y,1.15%
+JPY.TONA.OIS.10Y,1.32%
+JPY.TONA.TURN.2026,50bp
+JPY.TIBOR6M-TONA.BASIS.1Y,18bp
+JPY.TIBOR6M-TONA.BASIS.2Y,20bp
+JPY.TIBOR6M-TONA.BASIS.5Y,24bp
+JPY.TIBOR6M-TONA.BASIS.10Y,27bp
+JPY.TIBOR3M-TONA.BASIS.1Y,12bp
+JPY.TIBOR3M-TONA.BASIS.2Y,14bp
+JPY.TIBOR3M-TONA.BASIS.5Y,17bp
+JPY.TIBOR3M-TONA.BASIS.10Y,19bp
+USD.SOFR.OIS.1M,3.70%
+USD.SOFR.OIS.3M,3.62%
+USD.SOFR.OIS.6M,3.50%
+USD.SOFR.OIS.1Y,3.35%
+USD.SOFR.OIS.2Y,3.25%
+USD.SOFR.OIS.5Y,3.35%
+USD.SOFR.OIS.10Y,3.55%
+USDJPY.SPOT,152.40
+USDJPY.ON,-2.5
+USDJPY.TN,-2.4
+USDJPY.3M,-112
+USDJPY.6M,-218
+USDJPY.XCCY.1Y,-35bp
+USDJPY.XCCY.2Y,-40bp
+USDJPY.XCCY.5Y,-45bp
+JPY.LGM.MR,0.03
+JPY.LGM.VOL,0.006
diff --git a/examples/m1-rates/refdata/calendars.toml b/examples/m1-rates/refdata/calendars.toml
new file mode 100644
index 0000000000000000000000000000000000000000..c43442d3056bc5c19ebec73d7ba049100898214e
--- /dev/null
+++ b/examples/m1-rates/refdata/calendars.toml
@@ -0,0 +1,14 @@
+version = 1
+
+# 2026 年 11 月から 2027 年 12 月までは主な祝日を、それ以後は週末だけを休日にするサンプルの暦。
+[[calendar]]
+id = "JPTO"
+from = 2025-01-01
+through = 2037-12-31
+weekend = ["Sat", "Sun"]
+
+[[calendar]]
+id = "USNY"
+from = 2025-01-01
+through = 2037-12-31
+weekend = ["Sat", "Sun"]
diff --git a/examples/m1-rates/refdata/conventions.toml b/examples/m1-rates/refdata/conventions.toml
new file mode 100644
index 0000000000000000000000000000000000000000..4cfdefccad485647f1d197e9fdbee304da51d257
--- /dev/null
+++ b/examples/m1-rates/refdata/conventions.toml
@@ -0,0 +1,122 @@
+version = 1
+
+[[convention]]
+id = "JPY-TONA-OIS"
+
+[convention.kind.ois]
+spot_lag = "2BD"
+fixed = { frequency = "1Y", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
+payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
+
+[convention.kind.ois.float]
+index = "JPY-TONA"
+observation = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } }
+
+[[convention]]
+id = "USD-SOFR-OIS"
+
+[convention.kind.ois]
+spot_lag = "2BD"
+fixed = { frequency = "1Y", day_count = "ACT360", bdc = "modified_following", calendar = "USNY" }
+payment = { anchor = "period_end", offset = "2BD", calendar = "USNY", bdc = "following" }
+
+[convention.kind.ois.float]
+index = "USD-SOFR"
+observation = { days = { business_days = "USNY" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } }
+
+[[convention]]
+id = "JPY-TIBOR6M-TONA-BASIS"
+
+[convention.kind.basis_swap]
+spot_lag = "2BD"
+
+[[convention.kind.basis_swap.legs]]
+currency = "JPY"
+periods = { frequency = "6M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
+payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
+coupon = { compounded = { index = "JPY-TONA", observation = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }
+
+[[convention.kind.basis_swap.legs]]
+currency = "JPY"
+periods = { frequency = "6M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
+payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
+coupon = { term = { index = "JPY-TIBOR-6M" } }
+
+[[convention]]
+id = "JPY-TIBOR3M-TONA-BASIS"
+
+[convention.kind.basis_swap]
+spot_lag = "2BD"
+
+[[convention.kind.basis_swap.legs]]
+currency = "JPY"
+periods = { frequency = "3M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
+payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
+coupon = { compounded = { index = "JPY-TONA", observation = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }
+
+[[convention.kind.basis_swap.legs]]
+currency = "JPY"
+periods = { frequency = "3M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
+payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
+coupon = { term = { index = "JPY-TIBOR-3M" } }
+
+[[convention]]
+id = "JPY-TIBOR6M-IRS"
+
+[convention.kind.irs]
+spot_lag = "2BD"
+
+[[convention.kind.irs.legs]]
+currency = "JPY"
+periods = { frequency = "6M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
+payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
+coupon = "fixed"
+
+[[convention.kind.irs.legs]]
+currency = "JPY"
+periods = { frequency = "6M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
+payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
+coupon = { term = { index = "JPY-TIBOR-6M" } }
+
+[[convention]]
+id = "USDJPY-MTM-XCCY"
+
+[convention.kind.xccy_swap]
+spot_lag = "2BD"
+
+[[convention.kind.xccy_swap.legs]]
+currency = "JPY"
+periods = { frequency = "3M", day_count = "ACT365F", bdc = "modified_following", calendar = { union = ["JPTO", "USNY"] } }
+payment = { anchor = "period_end", offset = "0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" }
+coupon = { compounded = { index = "JPY-TONA", observation = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }
+
+[[convention.kind.xccy_swap.legs]]
+currency = "USD"
+periods = { frequency = "3M", day_count = "ACT360", bdc = "modified_following", calendar = { union = ["JPTO", "USNY"] } }
+payment = { anchor = "period_end", offset = "0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" }
+coupon = { compounded = { index = "USD-SOFR", observation = { days = { business_days = "USNY" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }
+fx_reset = "USDJPY-TKY-1500"
+
+[[convention]]
+id = "JPY-TIBOR6M-FRA"
+kind.fra = { index = "JPY-TIBOR-6M", spot_lag = "2BD" }
+
+[[convention]]
+id = "USDJPY-FWD"
+kind.fx_forward = { pair = "USDJPY", bdc = "modified_following" }
+
+[[convention]]
+id = "JPY-TONA3M-FUT"
+
+[convention.kind.rate_future]
+multiplier = 25_000_000.0
+settlement_lag = "1BD"
+settlement_time = { time = "10:00", zone = "Asia/Tokyo" }
+contracts = { "2027-03" = { settlement_index = "JPY-TONA3M-FUT-2027-03", final_clearing = 2027-06-15, rate = { compounded = { start = 2027-03-17, end = 2027-06-16 } } } }
+
+# USDJPY は T+2、東京 15 時の受渡、pip は 0.01 円。評価の 17 時には当日の受渡の締めを過ぎている。
+[[fx]]
+pair = "USDJPY"
+spot_lag = "two_business_days"
+settlement = { time = "15:00", zone = "Asia/Tokyo" }
+pip_exponent = 2
diff --git a/examples/m1-rates/refdata/currencies.toml b/examples/m1-rates/refdata/currencies.toml
new file mode 100644
index 0000000000000000000000000000000000000000..1ea9c184c156568b02958fd9c9ee9dba714df87d
--- /dev/null
+++ b/examples/m1-rates/refdata/currencies.toml
@@ -0,0 +1,17 @@
+version = 1
+
+[[currency]]
+code = "JPY"
+minor_units = 0
+payment_rounding = { increment = 1.0, mode = "half_even" }
+deliverability = "deliverable"
+settlement_time = { time = "17:00", zone = "Asia/Tokyo" }
+settlement_calendar = "JPTO"
+
+[[currency]]
+code = "USD"
+minor_units = 2
+payment_rounding = { increment = 0.01, mode = "half_even" }
+deliverability = "deliverable"
+settlement_time = { time = "03:00", zone = "America/New_York" }
+settlement_calendar = "USNY"
diff --git a/examples/m1-rates/refdata/holidays/JPTO.csv b/examples/m1-rates/refdata/holidays/JPTO.csv
new file mode 100644
index 0000000000000000000000000000000000000000..6342db3e3a3f9fc92aab6dabcd5ed17175ad3add
--- /dev/null
+++ b/examples/m1-rates/refdata/holidays/JPTO.csv
@@ -0,0 +1,21 @@
+date
+2026-11-03
+2026-11-23
+2026-12-31
+2027-01-01
+2027-01-11
+2027-02-11
+2027-02-23
+2027-03-22
+2027-04-29
+2027-05-03
+2027-05-04
+2027-05-05
+2027-07-19
+2027-08-11
+2027-09-20
+2027-09-23
+2027-10-11
+2027-11-03
+2027-11-23
+2027-12-31
diff --git a/examples/m1-rates/refdata/holidays/USNY.csv b/examples/m1-rates/refdata/holidays/USNY.csv
new file mode 100644
index 0000000000000000000000000000000000000000..3c2ff4ece1a0aca64ab1ef3377412b875b18704e
--- /dev/null
+++ b/examples/m1-rates/refdata/holidays/USNY.csv
@@ -0,0 +1,15 @@
+date
+2026-11-11
+2026-11-26
+2026-12-25
+2027-01-01
+2027-01-18
+2027-02-15
+2027-05-31
+2027-06-18
+2027-07-05
+2027-09-06
+2027-10-11
+2027-11-11
+2027-11-25
+2027-12-24
diff --git a/examples/m1-rates/refdata/indices.toml b/examples/m1-rates/refdata/indices.toml
new file mode 100644
index 0000000000000000000000000000000000000000..8378736f98a2532234fec5a1513b40f9142ced3f
--- /dev/null
+++ b/examples/m1-rates/refdata/indices.toml
@@ -0,0 +1,45 @@
+version = 1
+
+[[index]]
+id = "JPY-TONA"
+calendar = "JPTO"
+publication = { business_days = 1, calendar = "JPTO", time = "17:00", zone = "Asia/Tokyo" }
+observation_time = { time = "17:00", zone = "Asia/Tokyo" }
+kind = { overnight_rate = { ccy = "JPY", day_count = "ACT365F" } }
+
+[[index]]
+id = "USD-SOFR"
+calendar = "USNY"
+publication = { business_days = 1, calendar = "USNY", time = "03:00", zone = "America/New_York" }
+observation_time = { time = "03:00", zone = "America/New_York" }
+kind = { overnight_rate = { ccy = "USD", day_count = "ACT360" } }
+
+[[index]]
+id = "JPY-TIBOR-6M"
+calendar = "JPTO"
+publication = { business_days = 0, calendar = "JPTO", time = "11:00", zone = "Asia/Tokyo" }
+observation_time = { time = "11:00", zone = "Asia/Tokyo" }
+kind = { term_rate = { ccy = "JPY", day_count = "ACT365F", tenor = "6M", value_offset = "2BD", bdc = "modified_following", end_of_month = "ignore" } }
+
+[[index]]
+id = "JPY-TIBOR-3M"
+calendar = "JPTO"
+publication = { business_days = 0, calendar = "JPTO", time = "11:00", zone = "Asia/Tokyo" }
+observation_time = { time = "11:00", zone = "Asia/Tokyo" }
+kind = { term_rate = { ccy = "JPY", day_count = "ACT365F", tenor = "3M", value_offset = "2BD", bdc = "modified_following", end_of_month = "ignore" } }
+
+# 東京 15 時の USDJPY の公表値。MtM XCCY の元本のリセットに使う。
+[[index]]
+id = "USDJPY-TKY-1500"
+calendar = { union = ["JPTO", "USNY"] }
+publication = { business_days = 0, calendar = { union = ["JPTO", "USNY"] }, time = "15:00", zone = "Asia/Tokyo" }
+observation_time = { time = "15:00", zone = "Asia/Tokyo" }
+kind = { fx_fixing = { pair = "USDJPY" } }
+
+# 2027 年 3 月限の 3 か月 TONA 先物の日次清算価格。
+[[index]]
+id = "JPY-TONA3M-FUT-2027-03"
+calendar = "JPTO"
+publication = { business_days = 0, calendar = "JPTO", time = "15:30", zone = "Asia/Tokyo" }
+observation_time = { time = "15:30", zone = "Asia/Tokyo" }
+kind = { rate_futures_settlement = { contract = "JPY-TONA3M-2027-03", underlying = "JPY-TONA" } }
diff --git a/examples/m1-rates/refdata/securities.toml b/examples/m1-rates/refdata/securities.toml
new file mode 100644
index 0000000000000000000000000000000000000000..0b980752835c5d66922054895e4bf4d7bab274b5
--- /dev/null
+++ b/examples/m1-rates/refdata/securities.toml
@@ -0,0 +1,8 @@
+version = 1
+
+[[security]]
+id = "JPY-TONA3M-2027-03"
+kind = "listed_future"
+currency = "JPY"
+unit = "contract"
+valuation = "contract"
diff --git a/examples/m1-rates/refdata/zones.toml b/examples/m1-rates/refdata/zones.toml
new file mode 100644
index 0000000000000000000000000000000000000000..4f23aebcbd737b2094e987140dabe3f5f9c61680
--- /dev/null
+++ b/examples/m1-rates/refdata/zones.toml
@@ -0,0 +1,41 @@
+version = 1
+
+[[zone]]
+id = "Asia/Tokyo"
+from = 2025-01-01T00:00:00Z
+through = 2038-01-01T00:00:00Z
+initial_offset_seconds = 32400
+
+[[zone]]
+id = "America/New_York"
+from = 2025-01-01T00:00:00Z
+through = 2038-01-01T00:00:00Z
+initial_offset_seconds = -18000
+transitions = [
+    { at = 2025-03-09T07:00:00Z, offset_seconds = -14400 },
+    { at = 2025-11-02T06:00:00Z, offset_seconds = -18000 },
+    { at = 2026-03-08T07:00:00Z, offset_seconds = -14400 },
+    { at = 2026-11-01T06:00:00Z, offset_seconds = -18000 },
+    { at = 2027-03-14T07:00:00Z, offset_seconds = -14400 },
+    { at = 2027-11-07T06:00:00Z, offset_seconds = -18000 },
+    { at = 2028-03-12T07:00:00Z, offset_seconds = -14400 },
+    { at = 2028-11-05T06:00:00Z, offset_seconds = -18000 },
+    { at = 2029-03-11T07:00:00Z, offset_seconds = -14400 },
+    { at = 2029-11-04T06:00:00Z, offset_seconds = -18000 },
+    { at = 2030-03-10T07:00:00Z, offset_seconds = -14400 },
+    { at = 2030-11-03T06:00:00Z, offset_seconds = -18000 },
+    { at = 2031-03-09T07:00:00Z, offset_seconds = -14400 },
+    { at = 2031-11-02T06:00:00Z, offset_seconds = -18000 },
+    { at = 2032-03-14T07:00:00Z, offset_seconds = -14400 },
+    { at = 2032-11-07T06:00:00Z, offset_seconds = -18000 },
+    { at = 2033-03-13T07:00:00Z, offset_seconds = -14400 },
+    { at = 2033-11-06T06:00:00Z, offset_seconds = -18000 },
+    { at = 2034-03-12T07:00:00Z, offset_seconds = -14400 },
+    { at = 2034-11-05T06:00:00Z, offset_seconds = -18000 },
+    { at = 2035-03-11T07:00:00Z, offset_seconds = -14400 },
+    { at = 2035-11-04T06:00:00Z, offset_seconds = -18000 },
+    { at = 2036-03-09T07:00:00Z, offset_seconds = -14400 },
+    { at = 2036-11-02T06:00:00Z, offset_seconds = -18000 },
+    { at = 2037-03-08T07:00:00Z, offset_seconds = -14400 },
+    { at = 2037-11-01T06:00:00Z, offset_seconds = -18000 },
+]
diff --git a/examples/m1-rates/trades/csa.toml b/examples/m1-rates/trades/csa.toml
new file mode 100644
index 0000000000000000000000000000000000000000..7eeaa2773f5f8fdda1438107d0fbfc11fa8262bc
--- /dev/null
+++ b/examples/m1-rates/trades/csa.toml
@@ -0,0 +1,70 @@
+version = 1
+
+[[trade]]
+id = "USD-OIS-JPY-CASH"
+margin_set = "JPY-CASH"
+book = "CSA"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+convention = "USD-SOFR-OIS"
+tenor = "2Y"
+notional = "10_000_000 USD"
+fixed_rate = "3.40%"
+fixed_side = "receive"
+
+[[trade]]
+id = "USD-OIS-USD-CASH"
+margin_set = "USD-CASH"
+book = "CSA"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+convention = "USD-SOFR-OIS"
+tenor = "2Y"
+notional = "10_000_000 USD"
+fixed_rate = "3.40%"
+fixed_side = "receive"
+
+[[trade]]
+id = "USD-OIS-CTD"
+margin_set = "CTD-JPY-USD"
+book = "CSA"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+convention = "USD-SOFR-OIS"
+tenor = "2Y"
+notional = "10_000_000 USD"
+fixed_rate = "3.40%"
+fixed_side = "receive"
+
+[[trade]]
+id = "USD-OIS-UNCOLLATERALIZED"
+margin_set = "UNCOLLATERALIZED"
+book = "CSA"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+convention = "USD-SOFR-OIS"
+tenor = "2Y"
+notional = "10_000_000 USD"
+fixed_rate = "3.40%"
+fixed_side = "receive"
+
+[[trade]]
+id = "JPY-IRS-5Y-CASH-25"
+margin_set = "JPY-CASH-25"
+book = "CSA"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+legs = [
+  { side = "pay", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "1.10%" } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
+  { side = "receive", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { term = { index = "JPY-TIBOR-6M" } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
+]
diff --git a/examples/m1-rates/trades/fx.toml b/examples/m1-rates/trades/fx.toml
new file mode 100644
index 0000000000000000000000000000000000000000..ae0644c541d44c4842cc2a89c2d57c51d9e351a5
--- /dev/null
+++ b/examples/m1-rates/trades/fx.toml
@@ -0,0 +1,27 @@
+version = 1
+
+[[trade]]
+id = "USDJPY-FWD-1Y"
+margin_set = "JPY-CASH"
+book = "FX"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "fx_forward"
+pair = "USDJPY"
+notional = "10_000_000 USD"
+rate = 148.0
+settlement = 2027-11-29T06:00:00Z
+
+[[trade]]
+id = "USDJPY-FXSWAP-6M"
+margin_set = "JPY-CASH"
+book = "FX"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "fx_swap"
+pair = "USDJPY"
+notional = "10_000_000 USD"
+near = { rate = 152.40, settlement = 2026-11-27T06:00:00Z }
+far = { rate = 150.20, settlement = 2027-05-27T06:00:00Z }
diff --git a/examples/m1-rates/trades/rates.toml b/examples/m1-rates/trades/rates.toml
new file mode 100644
index 0000000000000000000000000000000000000000..16634bcd11f9f67f2549d2b3c1113ac7a9b6e1a7
--- /dev/null
+++ b/examples/m1-rates/trades/rates.toml
@@ -0,0 +1,125 @@
+version = 1
+
+[[trade]]
+id = "JPY-DEPO-3M"
+margin_set = "JPY-CASH"
+book = "JPY-MM"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "deposit"
+start = 2026-11-26
+maturity = 2027-02-26
+notional = "1_000_000_000 JPY"
+rate = "0.50%"
+day_count = "ACT365F"
+
+[[trade]]
+id = "JPY-FRA-6X12"
+margin_set = "JPY-CASH"
+book = "JPY-MM"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "fra"
+notional = "1_000_000_000 JPY"
+index = "JPY-TIBOR-6M"
+fixing = 2027-05-24
+fixed_rate = "0.75%"
+
+[[trade]]
+id = "JPY-TONA-FUT"
+margin_set = "JPY-CASH"
+book = "JPY-FUTURES"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "rate_future"
+settlement_index = "JPY-TONA3M-FUT-2027-03"
+contracts = 10
+trade_price = 0.9935
+day_count = "ACT365F"
+
+[[trade]]
+id = "JPY-IRS-5Y"
+margin_set = "JPY-CASH"
+book = "JPY-SWAPS"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+legs = [
+  { side = "pay", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "1.10%" } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
+  { side = "receive", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { term = { index = "JPY-TIBOR-6M" } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
+]
+
+[[trade]]
+id = "JPY-BASIS-3M-5Y"
+margin_set = "JPY-CASH"
+book = "JPY-SWAPS"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+legs = [
+  { side = "pay", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", spread = "20bp", rule = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "3M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
+  { side = "receive", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { term = { index = "JPY-TIBOR-3M" } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "3M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
+]
+
+[[trade]]
+id = "JPY-OIS-LOOKBACK"
+margin_set = "JPY-CASH"
+book = "JPY-SWAPS"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+legs = [
+  { side = "receive", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "0.70%" } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
+  { side = "pay", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", rule = { days = { business_days = "JPTO" }, method = { lookback = 2 }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
+]
+
+[[trade]]
+id = "JPY-OIS-SHIFT"
+margin_set = "JPY-CASH"
+book = "JPY-SWAPS"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+legs = [
+  { side = "receive", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "0.70%" } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
+  { side = "pay", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", rule = { days = { business_days = "JPTO" }, method = { observation_shift = 2 }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
+]
+
+[[trade]]
+id = "JPY-OIS-LOCKOUT"
+margin_set = "JPY-CASH"
+book = "JPY-SWAPS"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+legs = [
+  { side = "receive", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "0.70%" } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
+  { side = "pay", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", rule = { days = { business_days = "JPTO" }, method = { lockout = 1 }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
+]
+
+[[trade]]
+id = "JPY-IRS-AMENDED"
+margin_set = "JPY-CASH"
+book = "JPY-SWAPS"
+trade_date = 2026-05-28
+fees = [{ amount = "2_000_000 JPY", side = "receive", at = 2026-11-26T08:00:00Z }]
+
+[trade.product]
+type = "amended"
+from = 2026-11-24T08:00:00Z
+after = { scale = 0.4 }
+
+[trade.product.before]
+type = "swap"
+legs = [
+  { side = "pay", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "0.90%" } }, schedule = { start = 2026-06-01, end = { after_start = "3Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
+  { side = "receive", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { term = { index = "JPY-TIBOR-6M" } }, schedule = { start = 2026-06-01, end = { after_start = "3Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
+]
diff --git a/examples/m1-rates/trades/xccy.toml b/examples/m1-rates/trades/xccy.toml
new file mode 100644
index 0000000000000000000000000000000000000000..3e47ae80b4800b9a6f447600359940a6b3988055
--- /dev/null
+++ b/examples/m1-rates/trades/xccy.toml
@@ -0,0 +1,27 @@
+version = 1
+
+[[trade]]
+id = "USDJPY-MTM-XCCY-2Y"
+margin_set = "JPY-CASH"
+book = "XCCY"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+legs = [
+  { side = "pay", notional = { fixed = "1_524_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", spread = "-38bp", rule = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-27, end = { after_start = "2Y" }, step = "3M", stub = "no_stub", calendar = { union = ["JPTO", "USNY"] }, bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" } }, exchange = "initial_and_final" },
+  { side = "receive", notional = { fx_reset = { base = "1_524_000_000 JPY", fixing = "USDJPY-TKY-1500" } }, day_count = "ACT360", coupon = { compounded = { index = "USD-SOFR", rule = { days = { business_days = "USNY" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-27, end = { after_start = "2Y" }, step = "3M", stub = "no_stub", calendar = { union = ["JPTO", "USNY"] }, bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" } }, exchange = "initial_and_final" },
+]
+
+[[trade]]
+id = "USDJPY-XCCY-2Y"
+margin_set = "JPY-CASH"
+book = "XCCY"
+trade_date = 2026-11-24
+
+[trade.product]
+type = "swap"
+legs = [
+  { side = "pay", notional = { fixed = "1_524_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", spread = "-40bp", rule = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-27, end = { after_start = "2Y" }, step = "3M", stub = "no_stub", calendar = { union = ["JPTO", "USNY"] }, bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" } }, exchange = "initial_and_final" },
+  { side = "receive", notional = { fixed = "10_000_000 USD" }, day_count = "ACT360", coupon = { compounded = { index = "USD-SOFR", rule = { days = { business_days = "USNY" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-27, end = { after_start = "2Y" }, step = "3M", stub = "no_stub", calendar = { union = ["JPTO", "USNY"] }, bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" } }, exchange = "initial_and_final" },
+]
