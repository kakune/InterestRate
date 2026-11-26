diff --git a/crates/rhodium-market/tests/ui/projection_curve_as_discount.json b/crates/rhodium-market/tests/ui/projection_curve_as_discount.json
new file mode 100644
index 0000000000000000000000000000000000000000..f19f3758bcfbb737fdb5c082c0f683fbca1e6819
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/projection_curve_as_discount.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0308",
+    "message": "mismatched types",
+    "spans": [
+      {
+        "file": "tests/ui/projection_curve_as_discount.rs",
+        "line": 10,
+        "column": 39
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/projection_curve_as_discount.stderr b/crates/rhodium-market/tests/ui/projection_curve_as_discount.stderr
deleted file mode 100644
index 62cc11b9d0c7ec74b1c1b51368c4925b11a0914e..0000000000000000000000000000000000000000
--- a/crates/rhodium-market/tests/ui/projection_curve_as_discount.stderr
+++ /dev/null
@@ -1,15 +0,0 @@
-error[E0308]: mismatched types
-  --> tests/ui/projection_curve_as_discount.rs:10:39
-   |
-10 |         view.ccy_scope(Currency::JPY, projection, |_| Ok(()))
-   |              ---------                ^^^^^^^^^^ expected `CurveRef<'_, Discount>`, found `CurveRef<'_, Projection>`
-   |              |
-   |              arguments to this method are incorrect
-   |
-   = note: expected struct `CurveRef<'_, rhodium_market::view::Discount>`
-              found struct `CurveRef<'_, rhodium_market::view::Projection>`
-note: method defined here
-  --> src/view.rs
-   |
-   |     pub fn ccy_scope<T>(
-   |            ^^^^^^^^^
diff --git a/crates/rhodium-market/tests/ui/projection_position_as_discount_slot.json b/crates/rhodium-market/tests/ui/projection_position_as_discount_slot.json
new file mode 100644
index 0000000000000000000000000000000000000000..753d07eadba3694568e939f0733de93c2d46791f
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/projection_position_as_discount_slot.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0271",
+    "message": "type mismatch resolving `<[(ResolvedCurve<Projection>, Vec<_>); 1] as IntoIterator>::Item == (ResolvedCurve<Discount>, Vec<CurvePoint>)`",
+    "spans": [
+      {
+        "file": "tests/ui/projection_position_as_discount_slot.rs",
+        "line": 11,
+        "column": 45
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/projection_position_as_discount_slot.stderr b/crates/rhodium-market/tests/ui/projection_position_as_discount_slot.stderr
deleted file mode 100644
index 001dfa0d64585b8b29108dd507e33ebc341d3019..0000000000000000000000000000000000000000
--- a/crates/rhodium-market/tests/ui/projection_position_as_discount_slot.stderr
+++ /dev/null
@@ -1,19 +0,0 @@
-error[E0271]: type mismatch resolving `<[(ResolvedCurve<Projection>, Vec<_>); 1] as IntoIterator>::Item == (ResolvedCurve<Discount>, Vec<CurvePoint>)`
-  --> tests/ui/projection_position_as_discount_slot.rs:11:45
-   |
-11 |     Slots::from_resolved(structure.clone(), [(projection, Vec::new())], [])?;
-   |     --------------------                    ^^^^^^^^^^^^^^^^^^^^^^^^^^ expected `(ResolvedCurve<Discount>, Vec<_>)`, found `(ResolvedCurve<Projection>, Vec<_>)`
-   |     |
-   |     required by a bound introduced by this call
-   |
-   = note: expected tuple `(ResolvedCurve<rhodium_market::view::Discount>, Vec<CurvePoint>)`
-              found tuple `(ResolvedCurve<rhodium_market::view::Projection>, Vec<_>)`
-note: required by a bound in `Slots::from_resolved`
-  --> src/view.rs
-   |
-   | impl Slots {
-...
-   |     pub fn from_resolved(
-   |         structure: Arc<MarketStructure>,
-   |         discounts: impl IntoIterator<Item = (ResolvedCurve<Discount>, Vec<CurvePoint>)>,
-   |                                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ required by this bound in `Slots::from_resolved`
diff --git a/crates/rhodium-market/tests/ui/read_with_unbound_point.json b/crates/rhodium-market/tests/ui/read_with_unbound_point.json
new file mode 100644
index 0000000000000000000000000000000000000000..b2acabc13c05bbeb53dbdfe15a14e6fdd464adaf
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/read_with_unbound_point.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0308",
+    "message": "mismatched types",
+    "spans": [
+      {
+        "file": "tests/ui/read_with_unbound_point.rs",
+        "line": 9,
+        "column": 25
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/read_with_unbound_point.stderr b/crates/rhodium-market/tests/ui/read_with_unbound_point.stderr
deleted file mode 100644
index e2d3385f1df67be2acb415cd404b69841c0dbac5..0000000000000000000000000000000000000000
--- a/crates/rhodium-market/tests/ui/read_with_unbound_point.stderr
+++ /dev/null
@@ -1,15 +0,0 @@
-error[E0308]: mismatched types
- --> tests/ui/read_with_unbound_point.rs:9:25
-  |
-9 |         let _ = view.df(point);
-  |                      -- ^^^^^ expected `CurveAt<'_, Discount>`, found `CurvePoint`
-  |                      |
-  |                      arguments to this method are incorrect
-  |
-  = note: expected struct `CurveAt<'_, rhodium_market::view::Discount>`
-             found struct `CurvePoint`
-note: method defined here
- --> src/view.rs
-  |
-  |     pub const fn df(&self, at: CurveAt<'v, Discount, R>) -> DiscountFactor<R> {
-  |                  ^^
diff --git a/crates/rhodium-market/tests/ui/volatility_point_cross_view.json b/crates/rhodium-market/tests/ui/volatility_point_cross_view.json
new file mode 100644
index 0000000000000000000000000000000000000000..77be63f0a37f3eaca3a751bfcee280e07dd3505d
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/volatility_point_cross_view.json
@@ -0,0 +1,29 @@
+[
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/volatility_point_cross_view.rs",
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
+        "file": "tests/ui/volatility_point_cross_view.rs",
+        "line": 10,
+        "column": 35
+      },
+      {
+        "file": "tests/ui/volatility_point_cross_view.rs",
+        "line": 10,
+        "column": 35
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-market/tests/ui/volatility_point_cross_view.rs b/crates/rhodium-market/tests/ui/volatility_point_cross_view.rs
new file mode 100644
index 0000000000000000000000000000000000000000..9458aad3a208abef7d02afd0a52825478458e17b
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/volatility_point_cross_view.rs
@@ -0,0 +1,16 @@
+use rhodium_core::Result;
+use rhodium_market::{
+    market::Market,
+    vol::{VolId, VolPoint},
+};
+
+fn volatility(market: &Market, id: &VolId, point: &VolPoint) -> Result<f64> {
+    market.with_view(|first| {
+        let at = first.bind_vol(first.volatility(id)?, point)?;
+        market.with_view(|second| second.vol_at(at, 0.01))
+    })
+}
+
+fn main() {
+    let _ = volatility;
+}
diff --git a/crates/rhodium-market/tests/ui/volatility_point_cross_view_ok.rs b/crates/rhodium-market/tests/ui/volatility_point_cross_view_ok.rs
new file mode 100644
index 0000000000000000000000000000000000000000..5dd9727c335504fd0fe2b37c1c3ce95034779556
--- /dev/null
+++ b/crates/rhodium-market/tests/ui/volatility_point_cross_view_ok.rs
@@ -0,0 +1,16 @@
+use rhodium_core::Result;
+use rhodium_market::{
+    market::Market,
+    vol::{VolId, VolPoint},
+};
+
+fn volatility(market: &Market, id: &VolId, point: &VolPoint) -> Result<f64> {
+    market.with_view(|view| {
+        let at = view.bind_vol(view.volatility(id)?, point)?;
+        view.vol_at(at, 0.01)
+    })
+}
+
+fn main() {
+    let _ = volatility;
+}
diff --git a/crates/rhodium-pricing/Cargo.toml b/crates/rhodium-pricing/Cargo.toml
index ff79b52a12ec1267d1666f48d1bdb4b0a8e1c36e..ddfa6866ea830ba58bfe5d5b38310a963920ae8d 100644
--- a/crates/rhodium-pricing/Cargo.toml
+++ b/crates/rhodium-pricing/Cargo.toml
@@ -17,8 +17,10 @@ rhodium-market.workspace = true
 scioxide-ad.workspace = true
 scioxide-core.workspace = true
 scioxide-linalg.workspace = true
+scioxide-optimize.workspace = true
 scioxide-rootfind.workspace = true
 scioxide-sparse.workspace = true
+scioxide-special.workspace = true
 serde.workspace = true
 
 [dev-dependencies]
diff --git a/crates/rhodium-pricing/src/analytic.rs b/crates/rhodium-pricing/src/analytic.rs
new file mode 100644
index 0000000000000000000000000000000000000000..c0d16d23580a8bce36234220f435ff6ab590b116
--- /dev/null
+++ b/crates/rhodium-pricing/src/analytic.rs
@@ -0,0 +1,86 @@
+//! 金利オプションの校正と解析価格に共通する計算。
+
+use rhodium_core::{
+    Error, Result,
+    quantity::VolKind,
+    time::{Instant, model_time},
+};
+use scioxide_ad::{Mask, Real};
+
+/// Lyashenko–Mercurio の後決め RFR の実効分散時間（年）。
+///
+/// 観測開始前は通常の時間に観測期間の 1/3 を加え、期間内は残り時間の
+/// 三乗に比例する。観測終了後は 0 を返す。
+///
+/// # Errors
+///
+/// 観測期間が正でない場合。
+pub fn rfr_caplet_variance_time(as_of: Instant, start: Instant, end: Instant) -> Result<f64> {
+    let length = model_time(start, end).years();
+    if length <= 0.0 {
+        return Err(Error::InvalidInput(
+            "RFR observation interval must be positive".into(),
+        ));
+    }
+    if as_of <= start {
+        Ok(model_time(as_of, start).years() + length / 3.0)
+    } else {
+        let remaining = model_time(as_of, end).years().max(0.0);
+        Ok(remaining * (remaining / length).powi(2) / 3.0)
+    }
+}
+
+// 同じ F と K の二つの価格の差では本源的価値が消える。時間価値を直接求めると
+// ディープ ITM で、二つの大きな価格から小さい残差を引き出す必要が無い。
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#202 で各 variant を網羅する"
+)]
+pub(crate) fn option_time_value<R: Real>(
+    kind: VolKind,
+    forward: R,
+    strike: R,
+    stddev: R,
+) -> Result<R> {
+    let call = forward.le(strike);
+    let sign = R::select(call, R::one(), -R::one());
+    let active = stddev.gt(R::zero());
+    let safe = R::select(active, stddev, R::one());
+    let diff = forward - strike;
+    let value = match kind {
+        VolKind::Normal => {
+            let d = sign * diff / safe;
+            sign * diff * cdf(d)? + stddev * density(d)
+        }
+        VolKind::Lognormal | VolKind::ShiftedLognormal { .. } => {
+            let shift = match kind {
+                VolKind::ShiftedLognormal { shift } => shift,
+                _ => 0.0,
+            };
+            let f = forward + shift;
+            let k = strike + shift;
+            if !(f.gt(R::zero()) & k.gt(R::zero())).all() {
+                return Err(Error::Numerical(
+                    "Black requires a positive shifted forward and strike".into(),
+                ));
+            }
+            let d = ((f - k) / k).ln_1p() / safe;
+            let d1 = d + safe * 0.5;
+            let d2 = d - safe * 0.5;
+            sign * (f * cdf(sign * d1)? - k * cdf(sign * d2)?)
+        }
+    };
+    Ok(R::select(active, value, R::zero()))
+}
+
+fn density<R: Real>(x: R) -> R {
+    (x.powi(2) * -0.5).exp() * (1.0 / (2.0 * std::f64::consts::PI).sqrt())
+}
+
+fn cdf<R: Real>(x: R) -> Result<R> {
+    x.try_map_scalar(|v| {
+        let phi = density(v);
+        Ok([scioxide_special::normal_cdf(v)?, phi, -v * phi])
+    })
+    .map_err(|error| Error::Numerical(format!("normal CDF: {error}")))
+}
diff --git a/crates/rhodium-pricing/src/calibration.rs b/crates/rhodium-pricing/src/calibration.rs
index c0857ee0d7d3abcb9809537281879f3fd88ecc36..78d97fb9480370dfc2c3dc69fbe15a7dba338e57 100644
--- a/crates/rhodium-pricing/src/calibration.rs
+++ b/crates/rhodium-pricing/src/calibration.rs
@@ -2,6 +2,7 @@
 
 mod groups;
 mod solver;
+mod vol;
 
 #[cfg(test)]
 #[allow(
@@ -24,6 +25,7 @@ use rhodium_core::{
     fx::CcyPair,
 };
 use rhodium_market::{
+    curve::CurveRole,
     fx::{FxForwardInterval, FxSpec, FxTodayPlan},
     market::Market,
     quote::{QuoteId, QuoteKind, QuoteSet},
@@ -36,15 +38,29 @@ use scioxide_sparse::CsrMatrix;
 
 use crate::{
     context::Context,
-    definition::{ElementSelector, InstrumentBinding, MarketDefinition},
+    definition::{ElementSelector, InstrumentBinding, MarketDefinition, ResolvedFxDefinition},
     linear::OisWarmStart,
     par,
     prepare::{PrepareInput, Prepared, evaluate, prepare},
 };
 
-use self::groups::GroupPlan;
+use self::groups::{GroupForm, GroupPlan};
 
-/// 正方校正グループの解の線形化。行は [`CalibrationGroup::row_quotes`] の順。
+/// 最小二乗のグループの解の記録の 2 階の項の扱い。
+#[derive(Debug, Clone, Copy, PartialEq, Eq)]
+pub enum HessianMode {
+    /// `∂G/∂θ = J_Fᵀ W J_F + Σ_i (WF)_i ∇²F_i + ∇²P`（既定）。1 次の quote 感応度が厳密になる。
+    Exact,
+    /// 残差の 2 階の項を省く `J_Fᵀ W J_F + ∇²P`。quote 感応度は近似になる。
+    GaussNewton,
+}
+
+/// 校正グループの解の線形化 `G(θ_k, q, θ_up) = 0`。
+///
+/// 正方系では G は残差 F で、quote の行は [`CalibrationGroup::row_quotes`] の順。
+/// ボラ校正は、末尾に曲線に依存するフォワードか ATM ストライクの方程式を持つ。
+/// 最小二乗では G は停留条件 J̃ᵀ r̃（r̃ は重み付けした残差と罰則の行）で、行は
+/// [`CalibrationGroup::blocks`] の θ の順。どちらの形でも随伴は同じ式で流れる。
 #[derive(Debug)]
 pub struct GroupSolution {
     dg_dtheta: LuDecomposition<f64>,
@@ -52,6 +68,10 @@ pub struct GroupSolution {
     dg_dup: CsrMatrix<f64>,
     residual_norm: f64,
     iterations: usize,
+    hessian: Option<HessianMode>,
+    residuals: Vec<f64>,
+    theta_error: Vec<f64>,
+    vol_butterfly_violation: bool,
 }
 
 impl GroupSolution {
@@ -89,6 +109,26 @@ impl GroupSolution {
     pub const fn iterations(&self) -> usize {
         self.iterations
     }
+
+    /// 最小二乗のグループの解の記録の 2 階の項の扱い。正方系は `None`。
+    #[must_use]
+    pub const fn hessian(&self) -> Option<HessianMode> {
+        self.hessian
+    }
+
+    /// 解での残差 `F_i = m_i(θ) − q̃_i`（重みを掛けない）。行は [`CalibrationGroup::row_quotes`] の順。
+    #[must_use]
+    pub fn residuals(&self) -> &[f64] {
+        &self.residuals
+    }
+
+    /// 数値の停止の誤差の見積もり δθ = (∂G/∂θ)⁻¹·G。列は [`CalibrationGroup::blocks`] の θ の順。
+    ///
+    /// ソルバーの誤差を感応度で伝えて比較の幅にするのに使う（品質の「ソルバーの誤差」）。
+    #[must_use]
+    pub fn theta_error(&self) -> &[f64] {
+        &self.theta_error
+    }
 }
 
 /// グループの校正の結果。
@@ -152,6 +192,23 @@ pub struct CalibratedMarket {
     prepared: Arc<PreparedCalibration>,
 }
 
+/// 校正された一行の市場ボラとその残差。
+#[derive(Debug, Clone, Copy)]
+pub struct VolatilityFit<'a> {
+    /// 面の ID。
+    pub surface: &'a rhodium_market::vol::VolId,
+    /// 入力 quote の ID。
+    pub quote: &'a QuoteId,
+    /// 校正に使うボラの種類。
+    pub kind: rhodium_core::quantity::VolKind,
+    /// 差分を絶対ボラに直した目標。
+    pub target: f64,
+    /// 市場が再現したボラ。
+    pub model: f64,
+    /// モデルのボラから目標を引いた値。
+    pub residual: f64,
+}
+
 impl CalibratedMarket {
     /// 利用できるブロックを持つ市場。失敗したグループのブロックは利用できない。
     #[must_use]
@@ -177,6 +234,19 @@ impl CalibratedMarket {
         &self.groups
     }
 
+    /// 解けた金利ボラの各行。失敗したグループは `groups()` から確認できる。
+    pub fn volatility_fit(&self) -> impl Iterator<Item = VolatilityFit<'_>> {
+        self.groups
+            .iter()
+            .filter_map(|group| match (&group.plan.form, group.outcome()) {
+                (GroupForm::Volatility(volatility), GroupOutcome::Solved(solution)) => {
+                    Some(volatility.fit(&self.quotes, solution))
+                }
+                _ => None,
+            })
+            .flatten()
+    }
+
     /// 約定日の調整、商品の準備、入力 quote の検査で発生した警告。
     #[must_use]
     pub fn warnings(&self) -> &[Warning] {
@@ -265,7 +335,7 @@ pub fn calibrate_with_residual_tolerance(
     let calendar = ctx
         .refdata
         .calendar(&CalendarSpec::Base(def.valuation_calendar.clone()))?;
-    let (trade_date, warnings) = adjust(date, BusinessDayConvention::Following, &calendar)?;
+    let (trade_date, mut warnings) = adjust(date, BusinessDayConvention::Following, &calendar)?;
 
     // テンプレートに合う慣行や観測金利の丸めは市場の定義の問題なので、
     // 準備の前に全 binding を確かめ、校正全体のエラーにする。
@@ -281,27 +351,14 @@ pub fn calibrate_with_residual_tolerance(
         .iter()
         .map(|spec| spec.resolve(ctx.as_of, ctx.refdata))
         .collect::<Result<Vec<_>>>()?;
-    let curves = def.curves.iter().cloned().chain(
-        fx.iter()
-            .filter_map(|pair| pair.fit.as_ref().map(|fit| fit.curve.clone())),
-    );
-    let structure = Arc::new(MarketStructure::with_fx_and_lgm(
-        ctx.as_of,
-        curves,
-        fx.iter().map(|pair| pair.spec.clone()),
-        def.lgm
-            .iter()
-            .map(|(&ccy, model)| (ccy, model.breaks.clone())),
-    )?);
-    // 校正商品の準備の失敗（公表済みの fixing の欠損、設定が許さない近似など）は、その binding が
-    // 当てはめるグループの失敗として残し、独立した曲線（別の通貨など）の校正を続ける。
+    let (structure, curve_blocks, volatility) =
+        build_market_structure(def, &fx, trade_date, &calendar, ctx)?;
     let mut rows: Vec<_> = def
         .bindings
         .iter()
         .map(|binding| BindingRow::new(binding, trade_date, &structure, ctx))
         .collect::<Result<_>>()?;
     add_lgm_mark_rows(def, &structure, ctx, &mut rows)?;
-    let mut warnings = warnings;
     for pair in fx {
         warnings.extend(pair.warnings);
         if !pair.marks.is_empty() {
@@ -335,18 +392,33 @@ pub fn calibrate_with_residual_tolerance(
                         pip: fit.pip,
                     },
                     warm_start: None,
+                    tolerance: None,
                 })));
             }
         }
     }
-    let warnings = warnings
+    let mut warnings: BTreeSet<_> = warnings
         .into_iter()
         .chain(rows.iter().flat_map(|row| row.warnings().cloned()))
         .collect();
-    let plans = groups::plan(structure.layout().blocks(), &rows, structure.layout())?
-        .into_iter()
-        .map(Arc::new)
-        .collect();
+    let least_squares = groups::LeastSquaresSpec {
+        penalties: &def.penalties,
+        gauss_newton: &def.gauss_newton,
+    };
+    let mut plans = groups::plan(
+        &structure.layout().blocks()[..curve_blocks],
+        &rows,
+        structure.layout(),
+        &least_squares,
+    )?
+    .into_iter()
+    .map(Arc::new)
+    .collect::<Vec<_>>();
+    for surface in volatility {
+        let (groups, notices) = surface.prepare(&structure, ctx)?;
+        plans.extend(groups.into_iter().map(Arc::new));
+        warnings.extend(notices);
+    }
     Arc::new(PreparedCalibration {
         structure,
         rows,
@@ -358,6 +430,67 @@ pub fn calibrate_with_residual_tolerance(
     .run(quotes, None)
 }
 
+type PreparedStructure = (Arc<MarketStructure>, usize, Vec<vol::ResolvedSurface>);
+
+fn build_market_structure(
+    def: &MarketDefinition,
+    fx: &[ResolvedFxDefinition],
+    trade_date: Date,
+    calendar: &rhodium_core::calendar::Calendar,
+    ctx: &Context<'_>,
+) -> Result<PreparedStructure> {
+    let volatility = def
+        .volatilities
+        .iter()
+        .map(|definition| {
+            vol::ResolvedSurface::resolve(definition, trade_date, calendar, ctx)
+                .map_err(|error| error.context(definition.id().as_str()))
+        })
+        .collect::<Result<Vec<_>>>()?;
+    let structure = build_structure(def, fx, ctx)?;
+    let curve_blocks = structure.layout().blocks().len();
+    let structure = Arc::new(
+        structure.with_volatilities(volatility.iter().map(|surface| surface.spec.clone()))?,
+    );
+    Ok((structure, curve_blocks, volatility))
+}
+
+/// 市場定義、評価時点で解決した FX の定義から、予想の測度を記録した市場構造を作る。
+fn build_structure(
+    def: &MarketDefinition,
+    fx: &[ResolvedFxDefinition],
+    ctx: &Context<'_>,
+) -> Result<MarketStructure> {
+    let curves = def.curves.iter().cloned().chain(
+        fx.iter()
+            .filter_map(|pair| pair.fit.as_ref().map(|fit| fit.curve.clone())),
+    );
+    let structure = MarketStructure::with_fx_and_lgm(
+        ctx.as_of,
+        curves,
+        fx.iter().map(|pair| pair.spec.clone()),
+        def.lgm
+            .iter()
+            .map(|(&ccy, model)| (ccy, model.breaks.clone())),
+    )?;
+    // 予想と割引の担保の違いは、予想曲線を当てる校正商品の担保で判定する。校正商品の準備も
+    // 同じ判定を使うので、binding の準備より前に記録する。割引の規則が構造に無い曲線や、
+    // 通貨の違う曲線を指す測度は、その曲線を読む取引の準備が報告するので、ここでは記録しない。
+    let measures: Vec<_> = def
+        .forecast_measures(ctx.config, ctx.refdata)
+        .into_iter()
+        .filter(|(projection, discount)| {
+            structure.curve(discount).is_some_and(|spec| {
+                spec.has_role(CurveRole::Discount)
+                    && structure
+                        .curve(projection)
+                        .is_some_and(|projection| projection.currency() == spec.currency())
+            })
+        })
+        .collect();
+    structure.with_forecast_measures(measures)
+}
+
 fn add_lgm_mark_rows(
     def: &MarketDefinition,
     structure: &MarketStructure,
@@ -432,12 +565,15 @@ impl PreparedCalibration {
             });
         }
         let market = Market::from_partial_theta(Arc::clone(&self.structure), theta, available)?;
+        // 最小二乗の当てはめの警告は解に由来するので、再利用したグループを含めて毎回求める。
+        let fit_warnings: Vec<_> = groups.iter().flat_map(least_squares_warnings).collect();
         let warnings = self
             .warnings
             .iter()
             .cloned()
             .chain(unit_warnings(&self.quote_kinds, quotes))
             .chain(market.fx_cross_warnings())
+            .chain(fit_warnings)
             .collect::<BTreeSet<_>>()
             .into_iter()
             .collect();
@@ -458,6 +594,11 @@ impl PreparedCalibration {
         quotes: &QuoteSet,
         previous_theta: Option<&[f64]>,
     ) -> GroupOutcome {
+        if let GroupForm::Volatility(volatility) = &plan.form
+            && let Some(error) = volatility.failure()
+        {
+            return GroupOutcome::Failed(error.clone());
+        }
         if let Some(error) = plan.rows.iter().find_map(|&row| self.rows[row].failure()) {
             return GroupOutcome::Failed(error.clone());
         }
@@ -484,6 +625,40 @@ impl PreparedCalibration {
     }
 }
 
+/// 解いた最小二乗のグループで、残差が許容幅 σ̃ を超えた quote と、Gauss–Newton で記録した quote。
+fn least_squares_warnings(group: &CalibrationGroup) -> Vec<Warning> {
+    if let (GroupForm::Volatility(volatility), GroupOutcome::Solved(solution)) =
+        (&group.plan.form, group.outcome())
+    {
+        return volatility.warnings(solution);
+    }
+    let (GroupForm::LeastSquares(plan), GroupOutcome::Solved(solution)) =
+        (&group.plan.form, group.outcome.as_ref())
+    else {
+        return Vec::new();
+    };
+    let outside = group
+        .plan
+        .row_quotes
+        .iter()
+        .zip(&plan.weights)
+        .zip(solution.residuals())
+        // 重みは 1/σ̃ なので、|F|·w > 1 が |F| > σ̃ と同じ。
+        .filter(|&((_, weight), residual)| (residual * weight).abs() > 1.0)
+        .map(|((quote, _), _)| Warning::CalibrationFitOutsideTolerance {
+            quote: quote.clone(),
+        });
+    let approximate = group
+        .plan
+        .row_quotes
+        .iter()
+        .filter(|_| plan.hessian == HessianMode::GaussNewton)
+        .map(|quote| Warning::CalibrationGaussNewton {
+            quote: quote.clone(),
+        });
+    outside.chain(approximate).collect()
+}
+
 /// 金利かスプレッドとして読む quote のうち、絶対値が 1 を超えるものを警告する。
 ///
 /// quote の種類は binding が決めるので、`QuoteSet` ではなく校正の入口で検査する。
@@ -602,6 +777,18 @@ impl BindingRow {
     }
 }
 
+/// 残差 m(θ) − q̃(q) を測る単位。
+///
+/// 正方系の解は残差の尺度によらないので、解きやすい変換後の単位（FX アウトライトは対数）を使う。
+/// 最小二乗の許容幅 σ は quote の単位なので、重み 1/σ̃、目的関数、警告を quote の単位に揃える。
+#[derive(Clone, Copy, Debug, PartialEq, Eq)]
+pub(crate) enum ResidualUnit {
+    /// 変換後の単位。
+    Transformed,
+    /// quote の単位。FX アウトライトの残差は価格の差 F − q になる。
+    Quote,
+}
+
 /// 残差を同じ市場で繰り返し評価するための準備。
 #[derive(Debug)]
 struct PreparedBinding {
@@ -609,6 +796,8 @@ struct PreparedBinding {
     target: BlockId,
     model: PreparedTarget,
     warm_start: Option<OisWarmStart>,
+    /// 最小二乗の許容幅 σ（quote の単位）。正方系の行は `None`。
+    tolerance: Option<f64>,
 }
 
 #[derive(Debug)]
@@ -652,9 +841,18 @@ impl PreparedBinding {
         ctx: &Context<'_>,
     ) -> Result<Self> {
         Self::prepare(binding, trade_date, structure, ctx)
+            .map(|prepared| Self {
+                tolerance: binding.tolerance(),
+                ..prepared
+            })
             .map_err(|error| error.context(format_args!("quote {}", binding.quote)))
     }
 
+    /// 直接マーク（θ − D·q）の行か。全行が直接マークの最小二乗は θ について強凸な二次式。
+    const fn is_direct(&self) -> bool {
+        matches!(self.model, PreparedTarget::Direct { .. })
+    }
+
     fn prepare(
         binding: &InstrumentBinding,
         trade_date: Date,
@@ -716,6 +914,7 @@ impl PreparedBinding {
                 quote,
                 target,
                 warm_start: None,
+                tolerance: None,
                 model: PreparedTarget::FxOutright {
                     ccy,
                     unit_strike: at_rate(1.0)?,
@@ -737,6 +936,7 @@ impl PreparedBinding {
             quote,
             target,
             warm_start,
+            tolerance: None,
             model: PreparedTarget::ParRate {
                 ccy,
                 zero,
@@ -755,6 +955,7 @@ impl PreparedBinding {
             quote,
             target,
             warm_start: None,
+            tolerance: None,
             model: PreparedTarget::Direct {
                 index,
                 scale,
@@ -816,14 +1017,25 @@ impl PreparedBinding {
         }
     }
 
-    fn quote_derivative(&self, quote: f64) -> f64 {
-        if matches!(self.model, PreparedTarget::FxOutright { .. }) {
-            quote.recip()
-        } else {
-            self.quote_scale()
+    /// quote の単位の残差の変化を、変換後の単位の残差の変化に直す係数（FX アウトライトは 1/q）。
+    fn transformed_per_quote(&self, quote: f64) -> f64 {
+        self.quote_derivative(quote, ResidualUnit::Transformed)
+            / self.quote_derivative(quote, ResidualUnit::Quote)
+    }
+
+    /// 変換 q̃ の quote についての微分 D。
+    fn quote_derivative(&self, quote: f64, unit: ResidualUnit) -> f64 {
+        match (&self.model, unit) {
+            (PreparedTarget::FxOutright { .. }, ResidualUnit::Transformed) => quote.recip(),
+            (PreparedTarget::FxOutright { .. }, ResidualUnit::Quote) => 1.0,
+            _ => self.quote_scale(),
         }
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#201 で各 variant を網羅する"
+    )]
     fn seed_rate(&self, quote: f64) -> f64 {
         match self.model {
             PreparedTarget::FxPoints { .. }
@@ -836,7 +1048,7 @@ impl PreparedBinding {
     }
 
     /// 残差 m(θ) − q̃(q)。市場から読む値に数値型を通す。
-    fn residual<R: Real>(&self, market: &Market<R>, quote: f64) -> Result<R> {
+    fn residual<R: Real>(&self, market: &Market<R>, quote: f64, unit: ResidualUnit) -> Result<R> {
         let model = match &self.model {
             PreparedTarget::ParRate {
                 ccy, zero, one, fx, ..
@@ -870,7 +1082,12 @@ impl PreparedBinding {
                 let pv = evaluate(unit_strike, market)?;
                 let log_model = pv.get(fx.other).ln() + fx.today.rate(market)?.rate().ln()
                     - (-pv.get(*ccy)).ln();
-                return Ok(log_model - R::constant(quote.ln()));
+                let log = log_model - R::constant(quote.ln());
+                return Ok(match unit {
+                    ResidualUnit::Transformed => log,
+                    // F − q = q·(exp(ln F − ln q) − 1)。差を対数のまま expm1 に通し、桁を落とさない。
+                    ResidualUnit::Quote => log.exp_m1() * quote,
+                });
             }
             PreparedTarget::Direct { index, .. } => market.theta()[*index],
             PreparedTarget::FxPoints { points, .. } => points.points(market)?,
@@ -909,7 +1126,10 @@ mod tests {
                 .rows
                 .iter()
                 .filter_map(BindingRow::prepared)
-                .any(|row| { row.residual(base.market(), 0.0).is_err() })
+                .any(|row| {
+                    row.residual(base.market(), 0.0, ResidualUnit::Transformed)
+                        .is_err()
+                })
         );
         let overrides = QuoteOverrides::new(
             [QuoteOverride::Set {
diff --git a/crates/rhodium-pricing/src/calibration/groups.rs b/crates/rhodium-pricing/src/calibration/groups.rs
index eee70f9742656d9eb828c6fb38b4cf291ce33f4a..8d13762b89270ee66888b2ce41766bea729487b1 100644
--- a/crates/rhodium-pricing/src/calibration/groups.rs
+++ b/crates/rhodium-pricing/src/calibration/groups.rs
@@ -10,8 +10,41 @@ use rhodium_market::{
 use scioxide_sparse::{CsrMatrix, GraphDirection, GraphPlan};
 
 use super::BindingRow;
+use crate::definition::StepPenalty;
 
-/// 一つの強連結成分を正方系として解くための、ブロック・θ の添字・残差の対応。
+/// グループの方程式の形。
+#[derive(Debug)]
+pub(super) enum GroupForm {
+    /// 残差の数が未知数と等しい正方系 F = 0。
+    Square,
+    /// 許容幅で重み付けした残差と罰則の二乗和を最小にする（停留条件 G = 0）。
+    LeastSquares(LeastSquaresPlan),
+    /// 金利ボラ。専用の解法は同じグループの記録を返す。
+    Volatility(Box<super::vol::VolGroup>),
+}
+
+/// 最小二乗のグループの重み、罰則の行、解の記録の作り方。
+#[derive(Debug)]
+pub(super) struct LeastSquaresPlan {
+    /// 残差の行ごとの重み 1/σ̃。σ̃ は残差の単位の許容幅。
+    pub(super) weights: Vec<f64>,
+    /// 罰則の行 √λ·(θ[upper] − θ[lower])/s。添字は自グループの未知数の中の位置。
+    pub(super) penalties: Vec<PenaltyRow>,
+    /// 解の記録の 2 階の項の扱い。
+    pub(super) hessian: super::HessianMode,
+    /// すべての行が直接マークで、目的が θ について強凸な二次式か（解が一意で warm start できる）。
+    pub(super) convex: bool,
+}
+
+/// 隣り合う会合の区間のフォワードの差の罰則の一行。
+#[derive(Debug, Clone, Copy)]
+pub(super) struct PenaltyRow {
+    pub(super) lower: usize,
+    pub(super) upper: usize,
+    pub(super) factor: f64,
+}
+
+/// 一つの強連結成分を解くための、ブロック・θ の添字・残差の対応。
 ///
 /// 各添字列の順は、そのまま [`super::GroupSolution`] の行列の行と列の順になる。
 #[derive(Debug)]
@@ -30,12 +63,22 @@ pub(super) struct GroupPlan {
     pub(super) row_quotes: Vec<QuoteId>,
     /// `rows` が読む quote。ID の昇順で重複がなく、`dg_dq` の列の順。
     pub(super) quotes: Vec<QuoteId>,
+    /// 方程式の形。
+    pub(super) form: GroupForm,
+}
+
+/// 最小二乗の指定（罰則と Gauss–Newton のブロック）。
+pub(super) struct LeastSquaresSpec<'a> {
+    pub(super) penalties: &'a BTreeMap<BlockId, StepPenalty>,
+    pub(super) gauss_newton: &'a BTreeSet<BlockId>,
 }
 
 /// 校正商品の足跡から依存グラフを作り、上流から解く順にグループを並べる。
 ///
 /// 辺は、残差が読むブロックから当てはめ先のブロックへ張る（上流から下流）。
-/// 強連結成分を一つのグループとし、未知数と残差の数が一致しない成分は拒否する。
+/// 強連結成分を一つのグループとする。許容幅を持たない binding の成分は正方系で、
+/// 未知数と残差の数が一致しなければ拒否する。すべての binding が許容幅を持つ成分は最小二乗で、
+/// 残差と罰則の行の数が未知数以上でなければ拒否する。
 ///
 /// 準備に失敗した binding は足跡が分からないので辺を持たない。その binding を含む成分は
 /// 校正で `Failed` になり、下流（推移的）は `Unavailable` になるので、これらの成分では
@@ -45,6 +88,7 @@ pub(super) fn plan(
     blocks: &[BlockId],
     rows: &[BindingRow],
     layout: &ParamLayout,
+    least_squares: &LeastSquaresSpec<'_>,
 ) -> Result<Vec<GroupPlan>> {
     // 頂点が無いと依存グラフの行列を作れないので、曲線が無い定義は先に返す。
     if blocks.is_empty() {
@@ -112,20 +156,12 @@ pub(super) fn plan(
                 .collect();
             let indices = parameters(&blocks)?;
             let row_indices = std::mem::take(&mut rows_by_component[component]);
-            // 正方系として解くので、未知数と残差の数が一致しなければならない。
-            // 失敗した binding の影響を受ける成分は解かないので、数の食い違いを定義の誤りとしない。
-            if !tainted[component] && row_indices.len() != indices.len() {
-                return Err(Error::InvalidInput(format!(
-                    "calibration group [{}] has {} parameters and {} residuals",
-                    blocks
-                        .iter()
-                        .map(|block| block.curve().as_str())
-                        .collect::<Vec<_>>()
-                        .join(", "),
-                    indices.len(),
-                    row_indices.len()
-                )));
-            }
+            // 失敗した binding の影響を受ける成分は解かないので、形の検査をしない。
+            let form = if tainted[component] {
+                GroupForm::Square
+            } else {
+                group_form(&blocks, &indices, &row_indices, rows, layout, least_squares)?
+            };
             let row_quotes: Vec<_> = row_indices
                 .iter()
                 .map(|&index| rows[index].quote().clone())
@@ -144,11 +180,126 @@ pub(super) fn plan(
                 quotes,
                 blocks,
                 upstream,
+                form,
             })
         })
         .collect()
 }
 
+/// 成分の binding の許容幅と罰則から方程式の形を決め、行の数を検査する。
+fn group_form(
+    blocks: &[BlockId],
+    parameters: &[usize],
+    row_indices: &[usize],
+    rows: &[BindingRow],
+    layout: &ParamLayout,
+    spec: &LeastSquaresSpec<'_>,
+) -> Result<GroupForm> {
+    let names = || {
+        blocks
+            .iter()
+            .map(|block| block.curve().as_str())
+            .collect::<Vec<_>>()
+            .join(", ")
+    };
+    let prepared = row_indices
+        .iter()
+        .map(|&index| {
+            rows[index].prepared().ok_or_else(|| {
+                Error::Internal("untainted calibration group has a failed binding".into())
+            })
+        })
+        .collect::<Result<Vec<_>>>()?;
+    // 罰則の行は、ブロックの隣り合う区間ごとに一つ。添字は自グループの未知数の中の位置。
+    let penalties = blocks
+        .iter()
+        .filter_map(|block| spec.penalties.get(block).map(|penalty| (block, penalty)))
+        .map(|(block, penalty)| {
+            let range = layout
+                .range(block)
+                .ok_or_else(|| Error::Internal("penalized block has no range".into()))?;
+            let offset = parameters
+                .iter()
+                .position(|&index| index == range.start)
+                .ok_or_else(|| Error::Internal("penalized block outside its group".into()))?;
+            Ok((0..range.len() - 1).map(move |element| PenaltyRow {
+                lower: offset + element,
+                upper: offset + element + 1,
+                factor: penalty.factor(),
+            }))
+        })
+        .collect::<Result<Vec<_>>>()?
+        .into_iter()
+        .flatten()
+        .collect::<Vec<_>>();
+    let gauss_newton = blocks.iter().any(|block| spec.gauss_newton.contains(block));
+    let tolerances: Vec<_> = prepared.iter().map(|row| row.tolerance).collect();
+    match (
+        tolerances.iter().all(Option::is_some),
+        tolerances.iter().all(Option::is_none),
+    ) {
+        // 許容幅を持たない binding だけの成分は正方系。未知数と残差の数が一致しなければならない。
+        (_, true) => {
+            if !penalties.is_empty() || gauss_newton {
+                return Err(Error::InvalidInput(format!(
+                    "calibration group [{}] has a step penalty or Gauss-Newton option but no binding tolerances",
+                    names()
+                )));
+            }
+            if prepared.len() != parameters.len() {
+                return Err(Error::InvalidInput(format!(
+                    "calibration group [{}] has {} parameters and {} residuals",
+                    names(),
+                    parameters.len(),
+                    prepared.len()
+                )));
+            }
+            Ok(GroupForm::Square)
+        }
+        (true, false) => {
+            if prepared.len() + penalties.len() < parameters.len() {
+                return Err(Error::InvalidInput(format!(
+                    "least-squares calibration group [{}] has {} parameters but {} residual and penalty rows",
+                    names(),
+                    parameters.len(),
+                    prepared.len() + penalties.len()
+                )));
+            }
+            // 残差の単位の許容幅 σ̃ = σ·D（D は quote の変換の係数）の逆数を重みにする。
+            // 最小二乗の残差は quote の単位で測るので（`ResidualUnit::Quote`）、D は quote によらない。
+            let weights = prepared
+                .iter()
+                .zip(&tolerances)
+                .map(|(row, tolerance)| {
+                    tolerance
+                        .map(|sigma| (sigma * row.quote_scale()).recip())
+                        .ok_or_else(|| Error::Internal("missing binding tolerance".into()))
+                })
+                .collect::<Result<Vec<_>>>()?;
+            Ok(GroupForm::LeastSquares(LeastSquaresPlan {
+                weights,
+                penalties,
+                hessian: if gauss_newton {
+                    super::HessianMode::GaussNewton
+                } else {
+                    super::HessianMode::Exact
+                },
+                convex: prepared.iter().all(|row| row.is_direct()),
+            }))
+        }
+        (false, false) => Err(Error::InvalidInput(format!(
+            "calibration group [{}] mixes bindings with and without tolerances: {}",
+            names(),
+            prepared
+                .iter()
+                .filter(|row| row.tolerance.is_none())
+                .map(|row| row.quote.as_str())
+                .collect::<Vec<_>>()
+                .join(", ")
+        ))),
+    }
+}
+
 /// `len` 個のブロックと有向辺から強連結成分を求め、上流から解く順に並べる。
 ///
 /// 返り値の第一要素は、並べた順の各成分に属するブロックの添字（昇順）。
@@ -198,6 +349,7 @@ fn ordered_components(
         }
     }
     if result.len() != graph.count() {
+        // 不変条件: 成分のグラフは強連結成分を縮めた DAG なので、トポロジカル順はすべての成分を並べる。
         return Err(Error::Internal("cyclic calibration component graph".into()));
     }
     Ok((result, ordered_labels))
diff --git a/crates/rhodium-pricing/src/calibration/solver.rs b/crates/rhodium-pricing/src/calibration/solver.rs
index 0e5c79dda090e068355550a0582883d4797cf48e..2a3392d7c9110c6ba5273c90cb67eec2aee10d52 100644
--- a/crates/rhodium-pricing/src/calibration/solver.rs
+++ b/crates/rhodium-pricing/src/calibration/solver.rs
@@ -1,4 +1,8 @@
-//! 正方系のソルバーと、解におけるヤコビアンの記録。
+//! 校正グループのソルバー（正方系と最小二乗）と、解における線形化の記録。
+
+mod least_squares;
+
+pub(super) use least_squares::require_positive_definite;
 
 use std::{cell::Cell, collections::BTreeMap, sync::Arc};
 
@@ -18,7 +22,10 @@ use scioxide_rootfind::{
 };
 use scioxide_sparse::CsrMatrix;
 
-use super::{BindingRow, GroupSolution, PreparedBinding, groups::GroupPlan};
+use super::{
+    BindingRow, GroupSolution, PreparedBinding, ResidualUnit,
+    groups::{GroupForm, GroupPlan},
+};
 
 /// 一つのグループの正方系を解き、自グループの θ と解の記録を返す。
 ///
@@ -46,6 +53,9 @@ pub(super) fn solve(
             missing.join(", ")
         )));
     }
+    if let GroupForm::Volatility(volatility) = &plan.form {
+        return volatility.solve(plan, structure, theta, quotes, tolerance.absolute());
+    }
     let targets = plan
         .rows
         .iter()
@@ -72,11 +82,16 @@ pub(super) fn solve(
         targets,
         tolerance,
     };
+    // 直接マークの入力検査（有限で、FX スポットなら正）は解き方によらず先に行う。最小二乗の行でも
+    // 非正のスポットを平均に混ぜたり、解いた後に別の場所で失敗させたりしない。
     let direct = problem
         .targets
         .iter()
         .map(|&(row, quote)| row.direct_value(quote))
         .collect::<Result<Vec<_>>>()?;
+    if let GroupForm::LeastSquares(least_squares) = &plan.form {
+        return problem.solve_least_squares(least_squares, previous_theta);
+    }
     if direct.iter().all(Option::is_some) {
         let mut point = vec![0.0; plan.parameters.len()];
         for (index, value) in direct.into_iter().flatten() {
@@ -215,8 +230,13 @@ impl Problem<'_> {
         Ok(seed)
     }
 
-    /// θ の添字 `active` に `values` を入れた市場で、残差 F = par(θ) − q を行の順に返す。
-    fn residuals<R: Real>(&self, active: &[usize], values: &[R]) -> Result<Vec<R>> {
+    /// θ の添字 `active` に `values` を入れた市場で、残差 F = par(θ) − q を行の順に `unit` の単位で返す。
+    fn residuals<R: Real>(
+        &self,
+        active: &[usize],
+        values: &[R],
+        unit: ResidualUnit,
+    ) -> Result<Vec<R>> {
         // 未知数と必要な上流だけが AD の入力。他の値は定数として市場を組み立てる。
         let mut theta: Vec<_> = self.theta.iter().copied().map(R::constant).collect();
         for (&index, &value) in active.iter().zip(values) {
@@ -227,7 +247,7 @@ impl Problem<'_> {
         self.targets
             .iter()
             .map(|&(row, quote)| {
-                row.residual(&market, quote)
+                row.residual(&market, quote, unit)
                     .map_err(|error| error.context(format_args!("quote {}", row.quote)))
             })
             .collect()
@@ -249,7 +269,10 @@ impl Problem<'_> {
         let active = &self.plan.parameters;
         let callback_error = Cell::new(None);
         let function = ResidualFn::new(n, |x: &[f64], out: &mut [f64]| {
-            let residuals = forward_error(self.residuals(active, x), &callback_error)?;
+            let residuals = forward_error(
+                self.residuals(active, x, ResidualUnit::Transformed),
+                &callback_error,
+            )?;
             out.copy_from_slice(&residuals);
             Ok(())
         })
@@ -257,7 +280,7 @@ impl Problem<'_> {
             let matrix = jacobian(
                 |values, residuals| {
                     residuals.extend(forward_error(
-                        self.residuals(active, values),
+                        self.residuals(active, values, ResidualUnit::Transformed),
                         &callback_error,
                     )?);
                     Ok(())
@@ -324,7 +347,7 @@ impl Problem<'_> {
         let jac = jacobian(
             |values, residuals| {
                 residuals.extend(forward_error(
-                    self.residuals(&active, values),
+                    self.residuals(&active, values, ResidualUnit::Transformed),
                     &callback_error,
                 )?);
                 Ok(())
@@ -366,6 +389,12 @@ impl Problem<'_> {
         }
         let dg_dup = upstream_jacobian(n, active.len(), jac.as_slice())?;
         let dg_dq = self.quote_jacobian()?;
+        // 正方系の停止の誤差は、残差を一回のニュートン補正で θ に戻した δθ = J⁻¹F。
+        let theta_error = dg_dtheta
+            .solve(residual)
+            .map_err(calibration_error)?
+            .as_slice()
+            .to_vec();
         point.truncate(n);
         Ok((
             point,
@@ -375,6 +404,10 @@ impl Problem<'_> {
                 dg_dup,
                 residual_norm,
                 iterations,
+                hessian: None,
+                vol_butterfly_violation: false,
+                residuals: residual.to_vec(),
+                theta_error,
             },
         ))
     }
@@ -403,7 +436,7 @@ impl Problem<'_> {
             &self
                 .targets
                 .iter()
-                .map(|&(row, quote)| -row.quote_derivative(quote))
+                .map(|&(row, quote)| -row.quote_derivative(quote, ResidualUnit::Transformed))
                 .collect::<Vec<_>>(),
         )
         .map_err(calibration_error)
@@ -424,7 +457,11 @@ impl Problem<'_> {
 /// `columns` 列の行優先のヤコビアンから、先頭 n 列を除いた上流の列を疎行列にする。
 ///
 /// 残差が読まない上流の節点の列は厳密にゼロなので、格納しない。
-fn upstream_jacobian(n: usize, columns: usize, jacobian: &[f64]) -> Result<CsrMatrix<f64>> {
+pub(super) fn upstream_jacobian(
+    n: usize,
+    columns: usize,
+    jacobian: &[f64],
+) -> Result<CsrMatrix<f64>> {
     let mut rows = Vec::new();
     let mut cols = Vec::new();
     let mut values = Vec::new();
@@ -443,13 +480,13 @@ fn upstream_jacobian(n: usize, columns: usize, jacobian: &[f64]) -> Result<CsrMa
 }
 
 #[expect(clippy::needless_pass_by_value, reason = "map_err の変換関数に使う")]
-fn calibration_error(error: SciOxideError) -> Error {
+pub(super) fn calibration_error(error: SciOxideError) -> Error {
     Error::Calibration(error.to_string())
 }
 
 /// scioxide のコールバックは [`SciOxideError`] しか返せないので、
 /// 元の [`Error`] を `failure` に退避し、ソルバーが失敗したときに取り出せるようにする。
-fn forward_error<T>(
+pub(super) fn forward_error<T>(
     result: Result<T>,
     failure: &Cell<Option<Error>>,
 ) -> std::result::Result<T, SciOxideError> {
diff --git a/crates/rhodium-pricing/src/calibration/solver/least_squares.rs b/crates/rhodium-pricing/src/calibration/solver/least_squares.rs
new file mode 100644
index 0000000000000000000000000000000000000000..ca2d57bb52cae695fb27f0fd496091b445beaf6a
--- /dev/null
+++ b/crates/rhodium-pricing/src/calibration/solver/least_squares.rs
@@ -0,0 +1,484 @@
+//! 過剰決定と罰則付きの当てはめを最小二乗で解き、停留条件の線形化を記録する。
+//!
+//! 目的は `Φ(θ) = ½‖r̃(θ)‖²` で、`r̃` は重みを掛けた残差 `w_i·F_i`（`w_i = 1/σ̃_i`）と
+//! 罰則の行 `√λ·(f_{i+1} − f_i)/s` を並べたもの。停留条件 `G = ∇Φ = J̃ᵀ r̃ = J_Fᵀ W F + ∇P` が
+//! グループの方程式で、解の記録は G の線形化（`∂G/∂θ`、`∂G/∂q`、`∂G/∂θ_up`）にする。
+
+use std::cell::Cell;
+
+use rhodium_core::{Error, Result};
+use scioxide_ad::{Real, hvp, jacobian};
+use scioxide_core::Budget;
+use scioxide_linalg::{
+    CholeskyDecomposition, CholeskyPivoting, HermitianInput, HermitianProjectionPolicy,
+    LuDecomposition, Matrix, RankPolicy,
+};
+use scioxide_optimize::{LeastSquares, LeastSquaresFn, LeastSquaresOptions};
+use scioxide_sparse::CsrMatrix;
+
+use super::{Problem, calibration_error, forward_error};
+use crate::calibration::{GroupSolution, HessianMode, ResidualUnit, groups::LeastSquaresPlan};
+
+/// LM の後に、厳密なヘシアンのニュートン補正で停留条件を磨く回数の上限。
+///
+/// 残差が残る最小二乗では LM の収束は一次にとどまり、校正の許容誤差の手前で止まることがある。
+/// 停留点の近くのニュートン補正は二次で収束するので、数回で倍精度の分解能に届く。
+const POLISH_STEPS: usize = 8;
+
+/// 一つの点での最小二乗のグループの線形化。
+struct Linearization {
+    solution: GroupSolution,
+    /// 停留条件の誤差 δθ = (∂G/∂θ)⁻¹·G が各行を動かす量の、行の自然な単位でのノルム ‖Ĵ δθ‖。
+    movement: f64,
+    /// 目的 Φ = ½‖r̃‖²。初期値の候補から解を選ぶのに使う。
+    cost: f64,
+}
+
+/// 行の二乗和の半分 `½ Σ r_i²`。
+fn objective<R: Real>(rows: Vec<R>) -> R {
+    rows.into_iter()
+        .fold(R::constant(0.0), |sum, row| sum + row * row)
+        * R::constant(0.5)
+}
+
+impl Problem<'_> {
+    /// 最小二乗のグループを解き、停留点での線形化を返す。
+    ///
+    /// 目的が強凸（直接マークだけ）なら前の θ から解き、収束しなければ quote からの初期値で一度だけ
+    /// 解き直す。それ以外は固定の初期値の集合から解き、目的の値が最小の解を選ぶ。
+    pub(super) fn solve_least_squares(
+        &self,
+        plan: &LeastSquaresPlan,
+        previous: Option<&[f64]>,
+    ) -> Result<(Vec<f64>, GroupSolution)> {
+        let mut failures = Vec::new();
+        let flat = self.flat_seed()?;
+        let seeds = if let Some(theta) = previous.filter(|_| plan.convex) {
+            let warm: Vec<_> = self
+                .plan
+                .parameters
+                .iter()
+                .map(|&index| theta[index])
+                .collect();
+            match self.least_squares_from(plan, &warm) {
+                Ok((point, linear)) => return Ok((point, linear.solution)),
+                Err(error) => failures.push(error.to_string()),
+            }
+            vec![flat]
+        } else {
+            let mut seeds = vec![vec![0.0; self.plan.parameters.len()]];
+            if flat.iter().any(|value| *value != 0.0) {
+                seeds.push(flat);
+            }
+            seeds
+        };
+        let mut best = None::<(Vec<f64>, Linearization)>;
+        for seed in seeds {
+            match self.least_squares_from(plan, &seed) {
+                Ok((point, linear))
+                    if best
+                        .as_ref()
+                        .is_none_or(|(_, previous)| linear.cost < previous.cost) =>
+                {
+                    best = Some((point, linear));
+                }
+                Ok(_) => {}
+                Err(error) => failures.push(error.to_string()),
+            }
+        }
+        best.map(|(point, linear)| (point, linear.solution))
+            .ok_or_else(|| {
+                Error::Calibration(format!(
+                    "least-squares group [{}]: {}",
+                    self.plan
+                        .blocks
+                        .iter()
+                        .map(|block| block.curve().as_str())
+                        .collect::<Vec<_>>()
+                        .join(", "),
+                    failures.join(" | ")
+                ))
+            })
+    }
+
+    /// 重みを掛けた残差 `w_i·F_i` と罰則の行を、θ の添字 `active` に `values` を入れて求める。
+    ///
+    /// 罰則の行の添字は自グループの未知数の中の位置で、`active` の先頭が自グループの未知数。
+    fn weighted_rows<R: Real>(
+        &self,
+        plan: &LeastSquaresPlan,
+        active: &[usize],
+        values: &[R],
+    ) -> Result<Vec<R>> {
+        // 許容幅 σ は quote の単位なので、残差も quote の単位で測る（FX アウトライトは F − q）。
+        let residuals = self.residuals(active, values, ResidualUnit::Quote)?;
+        let data = residuals
+            .into_iter()
+            .zip(&plan.weights)
+            .map(|(residual, &weight)| residual * R::constant(weight));
+        let penalties = plan
+            .penalties
+            .iter()
+            .map(|row| (values[row.upper] - values[row.lower]) * R::constant(row.factor));
+        Ok(data.chain(penalties).collect())
+    }
+
+    /// 一つの初期値から Levenberg–Marquardt で解き、厳密なヘシアンで停留条件まで磨く。
+    fn least_squares_from(
+        &self,
+        plan: &LeastSquaresPlan,
+        seed: &[f64],
+    ) -> Result<(Vec<f64>, Linearization)> {
+        let own = &self.plan.parameters;
+        let rows = plan.weights.len() + plan.penalties.len();
+        let callback_error = Cell::new(None);
+        let function = LeastSquaresFn::new(seed.len(), rows, |x: &[f64], out: &mut [f64]| {
+            let values = forward_error(self.weighted_rows(plan, own, x), &callback_error)?;
+            out.copy_from_slice(&values);
+            Ok(())
+        })
+        .map_err(calibration_error)?
+        .with_jacobian(|x: &[f64], out: &mut Matrix<f64>| {
+            let matrix = jacobian(
+                |values, rows| {
+                    rows.extend(forward_error(
+                        self.weighted_rows(plan, own, values),
+                        &callback_error,
+                    )?);
+                    Ok(())
+                },
+                x,
+            )?;
+            out.as_mut_slice().copy_from_slice(matrix.as_slice());
+            Ok(())
+        });
+        let solver = LeastSquares::new(
+            LeastSquaresOptions::default().with_max_iterations(Budget::ONE.saturating_add(199)),
+        )
+        .map_err(calibration_error)?;
+        let result = solver.solve(&function, seed).map_err(|error| {
+            callback_error
+                .take()
+                .unwrap_or_else(|| calibration_error(error))
+        })?;
+        // 収束は LM の終了理由ではなく、記録の ∂G/∂θ による停留条件で判定する。
+        let mut point = result.x().to_vec();
+        for step in 0..=POLISH_STEPS {
+            let linear = self.linearize_least_squares(plan, &point, result.iterations() + step)?;
+            if linear.movement <= self.tolerance.absolute() {
+                // Gauss–Newton の記録は常に半正定値なので、極小かどうかは判定できない。解を採用する
+                // 条件は近似した記録と分け、停留点で一度だけ厳密なヘシアンの正定値を確かめる。
+                if plan.hessian == HessianMode::GaussNewton {
+                    self.require_exact_minimum(plan, &point)?;
+                }
+                return Ok((point, linear));
+            }
+            for (value, error) in point.iter_mut().zip(linear.solution.theta_error()) {
+                *value -= error;
+            }
+        }
+        Err(Error::Calibration(format!(
+            "least-squares stationarity not reached after LM {:?}",
+            result.termination()
+        )))
+    }
+
+    /// 点 `local` で停留条件 G を線形化し、解の記録と停止の誤差を作る。
+    ///
+    /// 列の順は自グループの n 個のあとに上流の要素が続く。厳密な記録は ∇²Φ の自グループの n 列を
+    /// Hessian–vector product で求め、その上の n 行が `∂G/∂θ`、下の行が（対称性から）`∂G/∂θ_up` の転置。
+    /// `∇²Φ = J̃ᵀJ̃ + Σ_i r̃_i ∇²r̃_i` なので、これは `J_Fᵀ W J_F + Σ_i (WF)_i ∇²F_i + ∇²P` と
+    /// 上流との交差項 `J_Fᵀ W ∂F/∂θ_up + Σ_i (WF)_i ∂²F_i/∂θ_k∂θ_up` を同時に与える。
+    fn linearize_least_squares(
+        &self,
+        plan: &LeastSquaresPlan,
+        local: &[f64],
+        iterations: usize,
+    ) -> Result<Linearization> {
+        let n = self.plan.parameters.len();
+        let mut active = self.plan.parameters.clone();
+        active.extend_from_slice(&self.plan.upstream_parameters);
+        let columns = active.len();
+        let mut point = local.to_vec();
+        point.extend(
+            self.plan
+                .upstream_parameters
+                .iter()
+                .map(|&index| self.theta[index]),
+        );
+        let callback_error = Cell::new(None);
+        let lift = |error| {
+            callback_error
+                .take()
+                .unwrap_or_else(|| calibration_error(error))
+        };
+        let jac = jacobian(
+            |values, rows| {
+                rows.extend(forward_error(
+                    self.weighted_rows(plan, &active, values),
+                    &callback_error,
+                )?);
+                Ok(())
+            },
+            &point,
+        )
+        .map_err(lift)?;
+        let weighted = self.weighted_rows::<f64>(plan, &active, &point)?;
+        let entry = |row: usize, column: usize| jac.as_slice()[row * columns + column];
+        let gradient: Vec<f64> = (0..n)
+            .map(|j| (0..weighted.len()).map(|i| entry(i, j) * weighted[i]).sum())
+            .collect();
+        let (hessian, cross) = self.second_derivatives(
+            plan,
+            plan.hessian,
+            &active,
+            &point,
+            jac.as_slice(),
+            &weighted,
+        )?;
+        let matrix = Matrix::from_vec(n, n, hessian).map_err(calibration_error)?;
+        // 停留点が極小でない（鞍点や極大）か、当てはめで識別できない向きがあれば解とみなさない。
+        // 分解は正定値の判定だけに使い、随伴には同じ行列の LU 分解を残す（正方系と同じ記録の型）。
+        // Gauss–Newton の J̃ᵀJ̃ では識別できない向きだけを検出し、極小の判定は採用の前に厳密に行う。
+        require_positive_definite(&matrix)?;
+        let dg_dtheta = LuDecomposition::factor(matrix.view()).map_err(calibration_error)?;
+        let rank = dg_dtheta
+            .analyze_rank(RankPolicy::DimensionScaledEpsilon)
+            .map_err(calibration_error)?
+            .rank();
+        if rank < n {
+            return Err(Error::Calibration(format!(
+                "least-squares Hessian rank {rank}/{n}"
+            )));
+        }
+        let theta_error = dg_dtheta
+            .solve(gradient.as_slice())
+            .map_err(calibration_error)?
+            .as_slice()
+            .to_vec();
+        let data = plan.weights.len();
+        // 停止の判定は正方系と同じ変換後の単位で測る（FX アウトライトは対数）。価格の単位のままでは
+        // 大きな値の丸め（150 の 1 ulp は約 3e-14）が校正の許容誤差を超え、停留を判定できない。
+        let units: Vec<f64> = self
+            .targets
+            .iter()
+            .map(|&(row, quote)| row.transformed_per_quote(quote))
+            .collect();
+        let movement = movement(plan, &units, jac.as_slice(), columns, &theta_error);
+        let residuals: Vec<f64> = weighted[..data]
+            .iter()
+            .zip(&plan.weights)
+            .map(|(row, weight)| row / weight)
+            .collect();
+        let residual_norm = residuals
+            .iter()
+            .map(|value| value * value)
+            .sum::<f64>()
+            .sqrt();
+        let cost = objective(weighted);
+        let solution = GroupSolution {
+            dg_dtheta,
+            dg_dq: self.least_squares_quote_jacobian(plan, n, columns, jac.as_slice())?,
+            dg_dup: dense_to_csr(n, columns - n, &cross)?,
+            residual_norm,
+            iterations,
+            hessian: Some(plan.hessian),
+            vol_butterfly_violation: false,
+            residuals,
+            theta_error,
+        };
+        Ok(Linearization {
+            solution,
+            movement,
+            cost,
+        })
+    }
+
+    /// Gauss–Newton で解いた停留点 `local` が目的関数の極小であることを、厳密なヘシアン
+    /// ∇²Φ の自グループの n 列を一度だけ求めて確かめる。
+    fn require_exact_minimum(&self, plan: &LeastSquaresPlan, local: &[f64]) -> Result<()> {
+        let n = self.plan.parameters.len();
+        // 上流は定数のまま、自グループの列だけを Hessian–vector product で回す。
+        let (hessian, _) = self.second_derivatives(
+            plan,
+            HessianMode::Exact,
+            &self.plan.parameters,
+            local,
+            &[],
+            &[],
+        )?;
+        require_positive_definite(&Matrix::from_vec(n, n, hessian).map_err(calibration_error)?)
+    }
+
+    /// 停留条件 G の自グループと上流の θ についての微分を、行優先の n × n と n × u の配列で返す。
+    ///
+    /// `hessian[i·n + j] = ∂G_i/∂θ_j`、`cross[j·u + k] = ∂G_j/∂θ_up,k`。`mode` が Gauss–Newton なら
+    /// 残差の 2 階の項を省き、`jacobian` と `weighted` から J̃ᵀJ̃ を作る（厳密なら両者を使わない）。
+    fn second_derivatives(
+        &self,
+        plan: &LeastSquaresPlan,
+        mode: HessianMode,
+        active: &[usize],
+        point: &[f64],
+        jacobian: &[f64],
+        weighted: &[f64],
+    ) -> Result<(Vec<f64>, Vec<f64>)> {
+        let n = self.plan.parameters.len();
+        let columns = active.len();
+        let upstream = columns - n;
+        let entry = |row: usize, column: usize| jacobian[row * columns + column];
+        let callback_error = Cell::new(None);
+        let mut hessian = vec![0.0; n * n];
+        let mut cross = vec![0.0; n * upstream];
+        for j in 0..n {
+            let column: Vec<f64> = match mode {
+                HessianMode::Exact => {
+                    let mut direction = vec![0.0; columns];
+                    direction[j] = 1.0;
+                    hvp(
+                        |values| {
+                            Ok(objective(forward_error(
+                                self.weighted_rows(plan, active, values),
+                                &callback_error,
+                            )?))
+                        },
+                        point,
+                        &direction,
+                    )
+                    .map_err(|error| {
+                        callback_error
+                            .take()
+                            .unwrap_or_else(|| calibration_error(error))
+                    })?
+                    .as_slice()
+                    .to_vec()
+                }
+                // 残差の 2 階の項を省く。罰則の行は θ について線形なので ∇²P は厳密に残る。
+                HessianMode::GaussNewton => (0..columns)
+                    .map(|k| (0..weighted.len()).map(|i| entry(i, k) * entry(i, j)).sum())
+                    .collect(),
+            };
+            for (i, &value) in column[..n].iter().enumerate() {
+                hessian[i * n + j] = value;
+            }
+            cross[j * upstream..(j + 1) * upstream].copy_from_slice(&column[n..]);
+        }
+        Ok((hessian, cross))
+    }
+
+    /// `∂G/∂q = −J_Fᵀ W D`（行は自グループの θ、列は `plan.quotes`）。
+    ///
+    /// `r̃_i = w_i·(m_i(θ) − D_i·q_i)` なので `∂r̃_i/∂q_i = −w_i·D_i` で、`J̃` は q に依存しない。
+    fn least_squares_quote_jacobian(
+        &self,
+        plan: &LeastSquaresPlan,
+        n: usize,
+        columns: usize,
+        jacobian: &[f64],
+    ) -> Result<CsrMatrix<f64>> {
+        let mut rows = Vec::new();
+        let mut cols = Vec::new();
+        let mut values = Vec::new();
+        for (i, (&(row, quote), weight)) in self.targets.iter().zip(&plan.weights).enumerate() {
+            let column = self.plan.quotes.binary_search(&row.quote).map_err(|_| {
+                Error::Internal(format!(
+                    "calibration group does not record quote {}",
+                    row.quote
+                ))
+            })?;
+            let scale = -weight * row.quote_derivative(quote, ResidualUnit::Quote);
+            for j in 0..n {
+                let value = jacobian[i * columns + j] * scale;
+                if value != 0.0 {
+                    rows.push(j);
+                    cols.push(column);
+                    values.push(value);
+                }
+            }
+        }
+        // 同じ quote を読む行（同じ要素への重ねた直接マークなど）の寄与は座標ごとに足される。
+        CsrMatrix::from_triplets(n, self.plan.quotes.len(), &rows, &cols, &values)
+            .map_err(calibration_error)
+    }
+}
+
+/// 対称な `matrix` が正定値でなければ、校正の失敗にする。
+pub(in crate::calibration) fn require_positive_definite(matrix: &Matrix<f64>) -> Result<()> {
+    CholeskyDecomposition::factor(
+        matrix.view(),
+        HermitianInput::Projected(HermitianProjectionPolicy::Always),
+        CholeskyPivoting::None,
+    )
+    .map(drop)
+    .map_err(|error| {
+        Error::Calibration(format!(
+            "least-squares Hessian is not positive definite: {error}"
+        ))
+    })
+}
+
+/// 停留条件の誤差 δθ が各行を動かす量のノルム ‖Ĵ δθ‖。
+///
+/// 行を自然な単位（正方系の解法と同じ変換後の残差、罰則はフォワードの差）に戻して測る。
+/// `units` は quote の単位の残差を変換後の単位に直す係数。正方系では J δθ = F なので、
+/// この量は残差のノルムに一致し、校正の許容誤差と同じ尺度で比べられる。
+fn movement(
+    plan: &LeastSquaresPlan,
+    units: &[f64],
+    jacobian: &[f64],
+    columns: usize,
+    theta_error: &[f64],
+) -> f64 {
+    plan.weights
+        .iter()
+        .zip(units)
+        .map(|(weight, unit)| weight / unit)
+        .chain(plan.penalties.iter().map(|row| row.factor))
+        .enumerate()
+        .map(|(i, scale)| {
+            let moved: f64 = theta_error
+                .iter()
+                .enumerate()
+                .map(|(j, error)| jacobian[i * columns + j] * error)
+                .sum();
+            (moved / scale).powi(2)
+        })
+        .sum::<f64>()
+        .sqrt()
+}
+
+/// 行優先の密行列のうちゼロでない要素を疎行列にする。
+fn dense_to_csr(rows: usize, columns: usize, values: &[f64]) -> Result<CsrMatrix<f64>> {
+    let (mut row_indices, mut column_indices, mut entries) = (Vec::new(), Vec::new(), Vec::new());
+    for (index, &value) in values
+        .iter()
+        .enumerate()
+        .filter(|(_, value)| **value != 0.0)
+    {
+        row_indices.push(index / columns.max(1));
+        column_indices.push(index % columns.max(1));
+        entries.push(value);
+    }
+    CsrMatrix::from_triplets(rows, columns, &row_indices, &column_indices, &entries)
+        .map_err(calibration_error)
+}
+
+#[cfg(test)]
+mod tests {
+    use super::*;
+
+    #[test]
+    fn acceptance_requires_a_positive_definite_hessian() -> Result<()> {
+        // 解を採用する条件は厳密なヘシアンの正定値で、Gauss–Newton の記録 J̃ᵀJ̃ ではない。
+        // ∇²Φ = J̃ᵀJ̃ + Σ r̃ ∇²r̃ が鞍点（固有値 3 と −1）なら、J̃ᵀJ̃ が正定値でも採用しない。
+        let gauss_newton =
+            Matrix::from_vec(2, 2, vec![2.0, 0.0, 0.0, 1.0]).map_err(calibration_error)?;
+        let saddle = Matrix::from_vec(2, 2, vec![1.0, 2.0, 2.0, 1.0]).map_err(calibration_error)?;
+        require_positive_definite(&gauss_newton)?;
+        assert!(matches!(
+            require_positive_definite(&saddle),
+            Err(Error::Calibration(message)) if message.contains("not positive definite")
+        ));
+        Ok(())
+    }
+}
diff --git a/crates/rhodium-pricing/src/calibration/vol.rs b/crates/rhodium-pricing/src/calibration/vol.rs
new file mode 100644
index 0000000000000000000000000000000000000000..a67eb95f6cdb55876acb332c6f630688d29f6635
--- /dev/null
+++ b/crates/rhodium-pricing/src/calibration/vol.rs
@@ -0,0 +1,412 @@
+//! 金利ボラ面を既存の校正グループと解の記録に接続する。
+
+mod caplet;
+mod swaption;
+
+use std::{cell::Cell, collections::BTreeSet, sync::Arc};
+
+use rhodium_core::{
+    Error, Result, Warning,
+    date::Date,
+    quantity::{LognormalVol, NormalVol, ShiftedLognormalVol, VolKind},
+};
+use rhodium_market::{
+    market::Market,
+    quote::{QuoteId, QuoteSet},
+    structure::{BlockId, MarketStructure, ParamLayout},
+    vol::{VolId, VolSpec},
+};
+use scioxide_ad::{Real, jacobian};
+use scioxide_linalg::{LuDecomposition, Matrix, RankPolicy};
+use scioxide_sparse::CsrMatrix;
+
+use crate::{
+    context::Context,
+    definition::{QuoteTransform, vol::VolatilityDefinition},
+};
+
+use super::{
+    GroupSolution, HessianMode,
+    groups::{GroupForm, GroupPlan},
+    solver::{calibration_error, forward_error, upstream_jacobian},
+};
+
+pub(super) struct ResolvedSurface {
+    pub(super) spec: VolSpec,
+    groups: Vec<UnboundGroup>,
+    pub(super) warnings: Vec<Warning>,
+}
+
+#[derive(Debug)]
+pub(super) struct VolGroup {
+    pub(super) surface: VolId,
+    pub(super) kind: VolKind,
+    quotes: Vec<VolQuote>,
+    model: Result<Model>,
+}
+
+struct UnboundGroup {
+    block: BlockId,
+    quotes: Vec<VolQuote>,
+    model: UnboundModel,
+}
+
+enum UnboundModel {
+    Swaption(swaption::UnboundSwaption),
+    Caplet(caplet::UnboundCaplet),
+}
+
+#[derive(Debug)]
+enum Model {
+    Swaption(Box<swaption::SwaptionNode>),
+    Caplet(Box<caplet::CapletGrid>),
+}
+
+#[derive(Debug)]
+struct VolQuote {
+    id: QuoteId,
+    transform: QuoteTransform,
+    tolerance: Option<f64>,
+}
+
+impl VolQuote {
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#201 で各 variant を網羅する"
+    )]
+    fn atm(&self) -> Option<&QuoteId> {
+        match &self.transform {
+            QuoteTransform::AtmPlusOffset(quote) => Some(quote),
+            _ => None,
+        }
+    }
+
+    fn target(&self, quotes: &QuoteSet) -> Option<f64> {
+        let value = quotes.get(&self.id)?;
+        match self.atm() {
+            Some(atm) => Some(value + quotes.get(atm)?),
+            None => Some(value),
+        }
+    }
+
+    fn columns<'a>(&'a self, columns: &'a [QuoteId]) -> impl Iterator<Item = Result<usize>> + 'a {
+        std::iter::once(&self.id).chain(self.atm()).map(|id| {
+            columns
+                .binary_search(id)
+                .map_err(|_| Error::Internal(format!("missing volatility quote dependency {id}")))
+        })
+    }
+}
+
+impl ResolvedSurface {
+    pub(super) fn resolve(
+        definition: &VolatilityDefinition,
+        trade_date: Date,
+        calendar: &rhodium_core::calendar::Calendar,
+        ctx: &Context<'_>,
+    ) -> Result<Self> {
+        definition.validate_references(ctx.refdata)?;
+        match definition {
+            VolatilityDefinition::SwaptionCube {
+                surface,
+                expiry_time,
+                sabr,
+                quotes,
+            } => swaption::resolve(
+                surface,
+                expiry_time,
+                *sabr,
+                quotes,
+                trade_date,
+                calendar,
+                ctx,
+            ),
+            VolatilityDefinition::Caplet {
+                surface, quotes, ..
+            } => caplet::resolve(surface, quotes, trade_date, ctx),
+        }
+    }
+
+    pub(super) fn prepare(
+        self,
+        structure: &Arc<MarketStructure>,
+        ctx: &Context<'_>,
+    ) -> Result<(Vec<GroupPlan>, Vec<Warning>)> {
+        let mut warnings = self.warnings;
+        let groups = self
+            .groups
+            .into_iter()
+            .map(|group| {
+                let result = match group.model {
+                    UnboundModel::Swaption(node) => node
+                        .prepare(structure, ctx)
+                        .map(|(node, warnings)| (Model::Swaption(Box::new(node)), warnings)),
+                    UnboundModel::Caplet(grid) => grid
+                        .prepare(structure, ctx)
+                        .map(|(grid, warnings)| (Model::Caplet(Box::new(grid)), warnings)),
+                };
+                let (model, upstream) = match result {
+                    Ok((model, notices)) => {
+                        warnings.extend(notices);
+                        let blocks = match &model {
+                            Model::Swaption(node) => node.blocks(),
+                            Model::Caplet(grid) => grid.blocks(),
+                        };
+                        (Ok(model), blocks.into_iter().collect::<Vec<_>>())
+                    }
+                    Err(
+                        error @ (Error::MissingData(_)
+                        | Error::Unsupported(_)
+                        | Error::Numerical(_)
+                        | Error::Calibration(_)),
+                    ) => (Err(error.context(self.spec.id().as_str())), vec![]),
+                    Err(error) => return Err(error.context(self.spec.id().as_str())),
+                };
+                let layout = structure.layout();
+                let parameters = indices(layout, std::slice::from_ref(&group.block))?;
+                let upstream_parameters = indices(layout, &upstream)?;
+                let row_quotes = group.quotes.iter().map(|q| q.id.clone()).collect();
+                let quotes = group
+                    .quotes
+                    .iter()
+                    .flat_map(|q| std::iter::once(q.id.clone()).chain(q.atm().cloned()))
+                    .collect::<BTreeSet<_>>()
+                    .into_iter()
+                    .collect();
+                Ok(GroupPlan {
+                    blocks: vec![group.block],
+                    upstream,
+                    parameters,
+                    upstream_parameters,
+                    rows: vec![],
+                    quotes,
+                    row_quotes,
+                    form: GroupForm::Volatility(Box::new(VolGroup {
+                        surface: self.spec.id().clone(),
+                        kind: self.spec.kind(),
+                        quotes: group.quotes,
+                        model,
+                    })),
+                })
+            })
+            .collect::<Result<Vec<_>>>()?;
+        Ok((groups, warnings))
+    }
+}
+
+fn indices(layout: &ParamLayout, blocks: &[BlockId]) -> Result<Vec<usize>> {
+    let mut indices = Vec::new();
+    for block in blocks {
+        indices.extend(layout.range(block).ok_or_else(|| {
+            Error::Internal(format!("missing volatility footprint block {block:?}"))
+        })?);
+    }
+    Ok(indices)
+}
+
+impl VolGroup {
+    pub(super) fn failure(&self) -> Option<&Error> {
+        self.model.as_ref().err()
+    }
+
+    pub(super) fn fit<'a>(
+        &'a self,
+        quotes: &'a QuoteSet,
+        solution: &'a GroupSolution,
+    ) -> impl Iterator<Item = super::VolatilityFit<'a>> {
+        self.quotes
+            .iter()
+            .zip(solution.residuals())
+            .filter_map(move |(row, &residual)| {
+                let target = row.target(quotes)?;
+                Some(super::VolatilityFit {
+                    surface: &self.surface,
+                    quote: &row.id,
+                    kind: self.kind,
+                    target,
+                    model: target + residual,
+                    residual,
+                })
+            })
+    }
+
+    pub(super) fn warnings(&self, solution: &GroupSolution) -> Vec<Warning> {
+        let mut warnings: Vec<_> = self
+            .quotes
+            .iter()
+            .zip(solution.residuals())
+            .filter(|(row, residual)| row.tolerance.is_some_and(|width| residual.abs() > width))
+            .map(|(row, _)| Warning::CalibrationFitOutsideTolerance {
+                quote: row.id.clone(),
+            })
+            .collect();
+        if solution.vol_butterfly_violation {
+            warnings.push(Warning::VolatilityButterfly {
+                surface: self.surface.as_str().to_owned(),
+            });
+        }
+        warnings
+    }
+
+    pub(super) fn solve(
+        &self,
+        plan: &GroupPlan,
+        structure: &Arc<MarketStructure>,
+        theta: &[f64],
+        quotes: &QuoteSet,
+        tolerance: f64,
+    ) -> Result<(Vec<f64>, GroupSolution)> {
+        let targets = self
+            .quotes
+            .iter()
+            .map(|row| {
+                let target = row.target(quotes).ok_or_else(|| {
+                    Error::MissingData(format!("missing volatility quote {}", row.id))
+                })?;
+                let parsed = match self.kind {
+                    VolKind::Normal => NormalVol::new(target).map(NormalVol::value),
+                    VolKind::Lognormal => LognormalVol::new(target).map(LognormalVol::value),
+                    VolKind::ShiftedLognormal { shift } => {
+                        ShiftedLognormalVol::new(target, shift).map(ShiftedLognormalVol::value)
+                    }
+                };
+                parsed.map_err(|error| {
+                    Error::InvalidData(format!("volatility quote {}: {error}", row.id))
+                })
+            })
+            .collect::<Result<Vec<_>>>()?;
+        let problem = VolProblem {
+            plan,
+            structure,
+            theta,
+            targets,
+            rows: &self.quotes,
+            tolerance,
+        };
+        match self.model.as_ref().map_err(Clone::clone)? {
+            Model::Swaption(node) => node.solve(&problem),
+            Model::Caplet(grid) => grid.solve(&problem),
+        }
+        .map_err(|error| error.context(self.surface.as_str()))
+    }
+}
+
+struct VolProblem<'a> {
+    plan: &'a GroupPlan,
+    structure: &'a Arc<MarketStructure>,
+    theta: &'a [f64],
+    targets: Vec<f64>,
+    rows: &'a [VolQuote],
+    tolerance: f64,
+}
+
+impl VolProblem<'_> {
+    fn market<R: Real>(&self, local: &[R], upstream: Option<&[R]>) -> Result<Market<R>> {
+        let mut theta: Vec<_> = self.theta.iter().copied().map(R::constant).collect();
+        for (&index, &value) in self.plan.parameters.iter().zip(local) {
+            theta[index] = value;
+        }
+        if let Some(upstream) = upstream {
+            for (&index, &value) in self.plan.upstream_parameters.iter().zip(upstream) {
+                theta[index] = value;
+            }
+        }
+        Market::from_partial_theta(
+            Arc::clone(self.structure),
+            theta,
+            self.plan.blocks.iter().chain(&self.plan.upstream).cloned(),
+        )
+    }
+
+    fn point(&self, local: &[f64]) -> Vec<f64> {
+        local
+            .iter()
+            .copied()
+            .chain(self.plan.upstream_parameters.iter().map(|&i| self.theta[i]))
+            .collect()
+    }
+
+    fn linearize(&self, linear: VolLinearization<'_>) -> Result<GroupSolution> {
+        let n = self.plan.parameters.len();
+        let columns = linear.theta.cols();
+        let own = linear
+            .theta
+            .as_slice()
+            .chunks_exact(columns)
+            .flat_map(|row| row[..n].iter().copied())
+            .collect();
+        let own = Matrix::from_vec(n, n, own).map_err(calibration_error)?;
+        let lu = LuDecomposition::factor(own.view()).map_err(calibration_error)?;
+        let rank = lu
+            .analyze_rank(RankPolicy::DimensionScaledEpsilon)
+            .map_err(calibration_error)?
+            .rank();
+        if rank != n {
+            return Err(Error::Calibration(format!(
+                "volatility calibration rank {rank}/{n}"
+            )));
+        }
+        let theta_error = lu
+            .solve(linear.equations.as_slice())
+            .map_err(calibration_error)?
+            .as_slice()
+            .to_vec();
+        let dg_dup = upstream_jacobian(n, columns, linear.theta.as_slice())?;
+        let mut rows = Vec::new();
+        let mut cols = Vec::new();
+        let mut coefficients = Vec::new();
+        for (i, coefficients_by_row) in linear.quote.iter().enumerate() {
+            for (row, &value) in self.rows.iter().zip(coefficients_by_row) {
+                if value != 0.0 {
+                    for column in row.columns(&self.plan.quotes) {
+                        rows.push(i);
+                        cols.push(column?);
+                        coefficients.push(value);
+                    }
+                }
+            }
+        }
+        let dg_dq =
+            CsrMatrix::from_triplets(n, self.plan.quotes.len(), &rows, &cols, &coefficients)
+                .map_err(calibration_error)?;
+        Ok(GroupSolution {
+            dg_dtheta: lu,
+            dg_dq,
+            dg_dup,
+            residual_norm: linear
+                .residuals
+                .iter()
+                .fold(0.0_f64, |norm, &r| norm.hypot(r)),
+            iterations: linear.iterations,
+            hessian: linear.hessian,
+            residuals: linear.residuals,
+            theta_error,
+            vol_butterfly_violation: false,
+        })
+    }
+}
+
+struct VolLinearization<'a> {
+    theta: &'a Matrix<f64>,
+    quote: Vec<Vec<f64>>,
+    equations: Vec<f64>,
+    residuals: Vec<f64>,
+    iterations: usize,
+    hessian: Option<HessianMode>,
+}
+
+// SciOxide のクロージャを抜けたところで元のドメインエラーを復元する。
+fn differentiate<F>(point: &[f64], f: F) -> Result<Matrix<f64>>
+where
+    F: for<'t> Fn(&[scioxide_ad::Var<'t, f64>]) -> Result<Vec<scioxide_ad::Var<'t, f64>>>,
+{
+    let error = Cell::new(None);
+    jacobian(
+        |x, out| {
+            out.extend(forward_error(f(x), &error)?);
+            Ok(())
+        },
+        point,
+    )
+    .map_err(|other| error.take().unwrap_or_else(|| calibration_error(other)))
+}
diff --git a/crates/rhodium-pricing/src/calibration/vol/caplet.rs b/crates/rhodium-pricing/src/calibration/vol/caplet.rs
new file mode 100644
index 0000000000000000000000000000000000000000..9eb895674b1ec37d7cf1a67d195092448a49b03e
--- /dev/null
+++ b/crates/rhodium-pricing/src/calibration/vol/caplet.rs
@@ -0,0 +1,686 @@
+//! キャップのフラットボラを再現するキャップレットの区分一定ボラ。
+
+use std::{
+    cell::Cell,
+    collections::{BTreeMap, BTreeSet},
+    sync::Arc,
+};
+
+use rhodium_core::{
+    Error, Result, Warning,
+    currency::Currency,
+    date::Date,
+    quantity::{Rate, Simple, VolKind},
+    refdata::{IndexKind, RefData},
+    tenor::Tenor,
+    time::{Instant, model_time},
+};
+use rhodium_market::{
+    market::Market,
+    structure::{BlockId, MarketStructure},
+    view::Slots,
+    vol::{CapletColumn, VolId, VolPoint, VolSpec},
+};
+use rhodium_products::standard::{
+    StandardProduct,
+    leg::{CapletCoupon, CapletObservation},
+};
+use scioxide_ad::{Dual, Real};
+use scioxide_core::{Interval, ResidualTolerance, StepTolerance};
+use scioxide_rootfind::{Brent, ScalarFn, ScalarRootOptions, ScalarRootSolver};
+
+use crate::{
+    analytic::{option_time_value, rfr_caplet_variance_time},
+    context::Context,
+    definition::QuoteTransform,
+    definition::vol::{CapVolQuote, RateStrike, VolatilityMarket},
+    prepare::{PrepareInput, Prepared, evaluate, prepare},
+    resolve::{Requirement, ResolvedRequirement, resolve as resolve_requirement},
+};
+
+use super::{
+    GroupSolution, ResolvedSurface, UnboundGroup, UnboundModel, VolLinearization, VolProblem,
+    VolQuote, calibration_error, differentiate, forward_error,
+};
+
+pub(super) struct UnboundCaplet {
+    surface: VolatilityMarket,
+    coupons: Vec<CapletCoupon>,
+    timing: Vec<Timing>,
+    rows: Vec<CapRow>,
+    atm_count: usize,
+    ccy: Currency,
+}
+
+#[derive(Debug)]
+pub(super) struct CapletGrid {
+    id: VolId,
+    kind: VolKind,
+    ccy: Currency,
+    coupons: Vec<Prepared>,
+    points: Vec<VolPoint>,
+    timing: Vec<Timing>,
+    rows: Vec<CapRow>,
+    absolute_strikes: Vec<f64>,
+    discounts: Slots,
+    fractions: Vec<f64>,
+    atm_count: usize,
+}
+
+#[derive(Debug, Clone, Copy)]
+struct Timing {
+    expiry: Instant,
+    variance_sqrt: f64,
+}
+
+#[derive(Debug)]
+struct CapRow {
+    strike: RateStrike,
+    start: usize,
+    end: usize,
+    // この列の各キャップレットが使う過去の θ の位置。
+    previous: Vec<usize>,
+    atm_strike: Option<usize>,
+}
+
+pub(super) fn resolve(
+    surface: &VolatilityMarket,
+    inputs: &[CapVolQuote],
+    trade_date: Date,
+    ctx: &Context<'_>,
+) -> Result<ResolvedSurface> {
+    let kind = surface.kind()?;
+    let definition = ctx.refdata.index(&surface.index)?;
+    let (IndexKind::TermRate { ccy, .. } | IndexKind::OvernightRate { ccy, .. }) = definition.kind
+    else {
+        return Err(Error::InvalidData(
+            "cap requires a term or overnight index".into(),
+        ));
+    };
+    let mut warnings = Vec::new();
+    let CapSchedules { coupons, counts } =
+        schedules(surface, inputs, trade_date, ctx, &mut warnings)?;
+    let timing = timings(ctx.as_of, &coupons, surface, ctx.refdata)?;
+    let mut rows: Vec<_> = inputs
+        .iter()
+        .map(|input| (input, counts[&Tenor::Length(input.maturity)]))
+        .collect();
+    rows.sort_by(|(a, end_a), (b, end_b)| {
+        strike_order(a.strike, b.strike).then_with(|| end_a.cmp(end_b))
+    });
+    let atm_count = rows
+        .iter()
+        .filter(|(q, _)| matches!(q.strike, RateStrike::Atm))
+        .count();
+    let mut model_rows = Vec::new();
+    let mut quotes = Vec::new();
+    let mut columns = Vec::new();
+    let mut atm = Vec::new();
+    for group in rows.chunk_by(|(a, _), (b, _)| strike_order(a.strike, b.strike).is_eq()) {
+        let strike = group[0].0.strike;
+        let mut previous = Vec::new();
+        let mut expiries = Vec::with_capacity(group.len());
+        for &(row, end) in group {
+            let start = previous.len();
+            if end <= start
+                || !timing[start..end]
+                    .iter()
+                    .any(|point| point.variance_sqrt > 0.0)
+            {
+                return Err(Error::InvalidData(format!(
+                    "missing independent caplet interval for {}",
+                    row.quote
+                )));
+            }
+            let atm_strike =
+                matches!(strike, RateStrike::Atm).then_some(model_rows.len() + atm_count);
+            model_rows.push(CapRow {
+                strike,
+                start,
+                end,
+                previous: previous.clone(),
+                atm_strike,
+            });
+            quotes.push(VolQuote {
+                id: row.quote.clone(),
+                transform: QuoteTransform::Identity,
+                tolerance: None,
+            });
+            previous.resize(end, model_rows.len() - 1);
+            expiries.push(timing[end - 1].expiry);
+        }
+        match strike {
+            RateStrike::Absolute(value) => columns.push(CapletColumn::new(
+                Rate::new(value, rhodium_core::day_count::DayCount::Act365F, Simple)?,
+                expiries,
+            )?),
+            RateStrike::Atm => atm = expiries,
+            RateStrike::AtmOffset(_) => {
+                return Err(Error::Internal(
+                    "cap ATM offset reached resolved grid".into(),
+                ));
+            }
+        }
+    }
+    let spec = VolSpec::caplet(
+        surface.id.clone(),
+        surface.index.clone(),
+        kind,
+        columns,
+        atm,
+    )?;
+    let group = UnboundGroup {
+        block: VolSpec::caplet_block(&surface.id)?,
+        quotes,
+        model: UnboundModel::Caplet(UnboundCaplet {
+            surface: surface.clone(),
+            coupons,
+            timing,
+            rows: model_rows,
+            atm_count,
+            ccy,
+        }),
+    };
+    Ok(ResolvedSurface {
+        spec,
+        groups: vec![group],
+        warnings,
+    })
+}
+
+// 最長のキャップと共通のクーポンを再利用できるか、満期ごとに日程を照合する。
+struct CapSchedules {
+    coupons: Vec<CapletCoupon>,
+    counts: BTreeMap<Tenor, usize>,
+}
+
+fn schedules(
+    surface: &VolatilityMarket,
+    inputs: &[CapVolQuote],
+    trade_date: Date,
+    ctx: &Context<'_>,
+    warnings: &mut Vec<Warning>,
+) -> Result<CapSchedules> {
+    let mut caps = BTreeMap::<_, Vec<_>>::new();
+    for input in inputs {
+        if let std::collections::btree_map::Entry::Vacant(entry) =
+            caps.entry(Tenor::Length(input.maturity))
+        {
+            let template = surface.template(input.maturity, ctx.refdata)?;
+            let StandardProduct::Swap(swap) =
+                template.build_for_trade(0.0, trade_date, ctx.refdata)?
+            else {
+                return Err(Error::Internal("cap expected a floating swap leg".into()));
+            };
+            let (coupons, notices) =
+                swap.caplet_coupons(&surface.index, trade_date, ctx.refdata)?;
+            warnings.extend(notices);
+            if coupons.is_empty() {
+                return Err(Error::InvalidData(format!(
+                    "empty cap after fixing convention at {}",
+                    input.quote
+                )));
+            }
+            entry.insert(coupons);
+        }
+    }
+    let longest = caps
+        .iter()
+        .max_by_key(|(_, coupons)| coupons.len())
+        .map(|(&tenor, _)| tenor)
+        .ok_or_else(|| Error::InvalidData("empty cap grid".into()))?;
+    let coupons = caps
+        .remove(&longest)
+        .ok_or_else(|| Error::Internal("missing longest cap".into()))?;
+    let mut counts = BTreeMap::from([(longest, coupons.len())]);
+    for (tenor, cap) in caps {
+        if !cap.iter().zip(&coupons).all(|(a, b)| {
+            a.pay == b.pay && a.observation == b.observation && a.fraction == b.fraction
+        }) {
+            return Err(Error::InvalidData(format!(
+                "cap maturity {tenor} does not share the preceding caplet schedule"
+            )));
+        }
+        counts.insert(tenor, cap.len());
+    }
+    Ok(CapSchedules { coupons, counts })
+}
+
+fn strike_order(a: RateStrike, b: RateStrike) -> std::cmp::Ordering {
+    match (a, b) {
+        (RateStrike::Absolute(a), RateStrike::Absolute(b)) => (a + 0.0).total_cmp(&(b + 0.0)),
+        (RateStrike::Absolute(_), _) => std::cmp::Ordering::Less,
+        (_, RateStrike::Absolute(_)) => std::cmp::Ordering::Greater,
+        _ => std::cmp::Ordering::Equal,
+    }
+}
+
+fn timings(
+    as_of: Instant,
+    coupons: &[CapletCoupon],
+    surface: &VolatilityMarket,
+    refdata: &RefData,
+) -> Result<Vec<Timing>> {
+    let event = refdata
+        .index(&surface.index)?
+        .observation_time
+        .as_ref()
+        .ok_or_else(|| {
+            Error::InvalidData(format!("missing observation time for {}", surface.index))
+        })?;
+    let zone = refdata.zone(&event.zone)?;
+    let clock = |date| zone.resolve_convention(event.on_date(date));
+    let timing = coupons
+        .iter()
+        .map(|coupon| {
+            // ターム金利の時計は fixing 自体の時刻。観測と公表が別の日でもここでは延ばさない。
+            let (expiry, variance_time) = match coupon.observation {
+                CapletObservation::Term(date) => {
+                    let expiry = clock(date)?;
+                    (expiry, model_time(as_of, expiry).years().max(0.0))
+                }
+                CapletObservation::Overnight(dates) => {
+                    let expiry = clock(dates.end())?;
+                    (
+                        expiry,
+                        rfr_caplet_variance_time(as_of, clock(dates.start())?, expiry)?,
+                    )
+                }
+            };
+            Ok(Timing {
+                expiry,
+                variance_sqrt: variance_time.sqrt(),
+            })
+        })
+        .collect::<Result<Vec<_>>>()?;
+    if !timing.is_sorted_by(|a, b| a.expiry < b.expiry) {
+        return Err(Error::InvalidData(
+            "caplet expiries must be strictly increasing".into(),
+        ));
+    }
+    Ok(timing)
+}
+
+impl UnboundCaplet {
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#201 で各 variant を網羅する"
+    )]
+    pub(super) fn prepare(
+        self,
+        structure: &Arc<MarketStructure>,
+        ctx: &Context<'_>,
+    ) -> Result<(CapletGrid, Vec<Warning>)> {
+        let ResolvedRequirement::Discounting(curve) = resolve_requirement(
+            &Requirement::Discounting {
+                ccy: self.ccy,
+                collateral: self.surface.collateral.clone(),
+            },
+            ctx.config,
+            structure,
+            ctx.refdata,
+        )?
+        else {
+            return Err(Error::Internal("expected cap discount curve".into()));
+        };
+        let mut warnings = Vec::new();
+        let mut coupons = Vec::new();
+        let mut points = Vec::new();
+        let mut discount_points = Vec::new();
+        let mut fractions = Vec::new();
+        let id = self.surface.id.clone();
+        for (coupon, timing) in self.coupons.into_iter().zip(&self.timing) {
+            if coupon.pay <= ctx.as_of {
+                return Err(Error::InvalidInput("cap payment must follow as-of".into()));
+            }
+            let prepared = prepare(
+                PrepareInput::Ir(&coupon.product),
+                &self.surface.collateral,
+                &[],
+                structure,
+                ctx,
+            )?;
+            warnings.extend_from_slice(prepared.warnings());
+            coupons.push(prepared);
+            points.push(structure.prepare_vol_point(&id, timing.expiry, None)?);
+            discount_points.push(curve.prepare_point(coupon.pay)?);
+            fractions.push(coupon.fraction.value());
+        }
+        let mut absolute_strikes: Vec<_> = self
+            .rows
+            .iter()
+            .filter_map(|row| match row.strike {
+                RateStrike::Absolute(value) => Some(value),
+                _ => None,
+            })
+            .collect();
+        absolute_strikes.dedup_by(|a, b| (*a + 0.0).total_cmp(&(*b + 0.0)).is_eq());
+        let discounts =
+            Slots::from_resolved(Arc::clone(structure), [(curve, discount_points)], [])?;
+        Ok((
+            CapletGrid {
+                kind: self.surface.kind()?,
+                id,
+                ccy: self.ccy,
+                coupons,
+                points,
+                timing: self.timing,
+                rows: self.rows,
+                absolute_strikes,
+                discounts,
+                fractions,
+                atm_count: self.atm_count,
+            },
+            warnings,
+        ))
+    }
+}
+
+struct CouponValue<R> {
+    forward: R,
+    annuity: R,
+    sqrt_time: f64,
+}
+
+impl CapletGrid {
+    pub(super) fn blocks(&self) -> BTreeSet<BlockId> {
+        self.coupons
+            .iter()
+            .flat_map(Prepared::blocks)
+            .cloned()
+            .collect()
+    }
+
+    fn values<R: Real>(&self, market: &Market<R>) -> Result<Vec<CouponValue<R>>> {
+        market.with_view(|view| {
+            let slots = view.bind_slots(&self.discounts)?;
+            let points = slots.discounts()[0].points();
+            self.coupons
+                .iter()
+                .zip(&self.fractions)
+                .zip(points)
+                .zip(&self.timing)
+                .map(|(((coupon, &fraction), &discount), timing)| {
+                    let annuity = view.df(discount).value() * fraction;
+                    let forward = evaluate(coupon, market)?.get(self.ccy) / annuity;
+                    Ok(CouponValue {
+                        forward,
+                        annuity,
+                        sqrt_time: timing.variance_sqrt,
+                    })
+                })
+                .collect()
+        })
+    }
+
+    fn atm<R: Real>(values: &[CouponValue<R>]) -> R {
+        let (pv, annuity) = values
+            .iter()
+            .fold((R::zero(), R::zero()), |(pv, annuity), item| {
+                (pv + item.forward * item.annuity, annuity + item.annuity)
+            });
+        pv / annuity
+    }
+
+    fn price<R: Real>(
+        &self,
+        values: &[CouponValue<R>],
+        k: R,
+        vols: impl IntoIterator<Item = R>,
+    ) -> Result<R> {
+        values
+            .iter()
+            .zip(vols)
+            .try_fold(R::zero(), |sum, (item, vol)| {
+                Ok(sum
+                    + item.annuity
+                        * option_time_value(self.kind, item.forward, k, vol * item.sqrt_time)?)
+            })
+    }
+
+    fn smile<R: Real>(&self, market: &Market<R>, k: R, count: usize) -> Result<Vec<R>> {
+        market.with_view(|view| {
+            let face = view.volatility(&self.id)?;
+            self.points[..count]
+                .iter()
+                .map(|point| view.vol_at(view.bind_vol(face, point)?, k))
+                .collect()
+        })
+    }
+
+    pub(super) fn solve(&self, problem: &VolProblem<'_>) -> Result<(Vec<f64>, GroupSolution)> {
+        let mut point = vec![0.0; self.rows.len() + self.atm_count];
+        let values = self.values(&problem.market(&point, None)?)?;
+        let mut iterations = 0;
+        for (i, row) in self.rows.iter().enumerate() {
+            let strike = row
+                .atm_strike
+                .map_or_else(|| row.strike.value(0.0), |_| Self::atm(&values[..row.end]));
+            if let Some(index) = row.atm_strike {
+                if self
+                    .absolute_strikes
+                    .iter()
+                    .any(|&k| (k - strike).abs() <= 1e-14)
+                {
+                    return Err(Error::Calibration(format!(
+                        "cap ATM strike coincides with an absolute node for {}",
+                        problem.rows[i].id
+                    )));
+                }
+                point[index] = strike;
+            }
+            let fixed_vols = if row.atm_strike.is_some() {
+                self.smile(&problem.market(&point, None)?, strike, row.start)?
+            } else {
+                row.previous.iter().map(|&index| point[index]).collect()
+            };
+            let fixed_price = self.price(&values[..row.start], strike, fixed_vols)?;
+            let target = self.price(
+                &values[..row.end],
+                strike,
+                std::iter::repeat(problem.targets[i]),
+            )?;
+            // 早い区間の価格を固定したまま、増加する最後の区間の価格を合わせる。
+            let price_at =
+                |vol| self.price(&values[row.start..row.end], strike, std::iter::repeat(vol));
+            let (vol, evaluations) =
+                positive_vol(&price_at, target - fixed_price, problem.targets[i])
+                    .map_err(|error| error.context(format_args!("quote {}", problem.rows[i].id)))?;
+            point[i] = vol;
+            iterations += evaluations;
+        }
+        self.finish(problem, point, &values, iterations)
+    }
+
+    fn finish(
+        &self,
+        problem: &VolProblem<'_>,
+        point: Vec<f64>,
+        values: &[CouponValue<f64>],
+        iterations: usize,
+    ) -> Result<(Vec<f64>, GroupSolution)> {
+        let market = problem.market(&point, None)?;
+        let n = point.len();
+        let mut quote = vec![vec![0.0; self.rows.len()]; n];
+        let mut residuals = Vec::new();
+        let mut vega = Vec::new();
+        for (i, row) in self.rows.iter().enumerate() {
+            let k = row
+                .atm_strike
+                .map_or_else(|| row.strike.value(0.0), |index| point[index]);
+            let price = self.price(&values[..row.end], k, self.smile(&market, k, row.end)?)?;
+            let flat_price = |sigma| self.price(&values[..row.end], k, std::iter::repeat(sigma));
+            let (model, _) = positive_vol(&flat_price, price, problem.targets[i])?;
+            residuals.push(model - problem.targets[i]);
+            quote[i][i] = -1.0;
+
+            let dual = values[..row.end]
+                .iter()
+                .map(|item| CouponValue {
+                    forward: Dual::lift(item.forward),
+                    annuity: Dual::lift(item.annuity),
+                    sqrt_time: item.sqrt_time,
+                })
+                .collect::<Vec<_>>();
+            let slope = self
+                .price(
+                    &dual,
+                    Dual::lift(k),
+                    std::iter::repeat(Dual::variable(model)),
+                )?
+                .derivative();
+            if !(slope.is_finite() && slope > 0.0) {
+                return Err(Error::Calibration(format!(
+                    "unidentified cap flat volatility for {}",
+                    problem.rows[i].id
+                )));
+            }
+            vega.push((model, slope));
+        }
+        let norm = residuals.iter().fold(0.0_f64, |norm, &r| norm.hypot(r));
+        if norm > problem.tolerance {
+            return Err(Error::Calibration(format!(
+                "stripped cap quote residual {norm:e} exceeds tolerance"
+            )));
+        }
+        let jac = differentiate(&problem.point(&point), |x| {
+            let (local, upstream) = x.split_at(n);
+            let market = problem.market(local, Some(upstream))?;
+            let values = self.values(&market)?;
+            let mut equations = Vec::with_capacity(n);
+            let mut anchors = Vec::with_capacity(self.atm_count);
+            for (row, &(model, slope)) in self.rows.iter().zip(&vega) {
+                let k = row.atm_strike.map_or_else(
+                    || row.strike.value(scioxide_ad::Var::zero()),
+                    |index| local[index],
+                );
+                let value = self.price(&values[..row.end], k, self.smile(&market, k, row.end)?)?;
+                let flat = self.price(
+                    &values[..row.end],
+                    k,
+                    std::iter::repeat(scioxide_ad::Var::constant(model)),
+                )?;
+                // PV の差をフラットボラのベガで割るのが逆算ボラの陰関数微分。
+                equations.push((value - flat) / slope);
+                if let Some(index) = row.atm_strike {
+                    anchors.push(local[index] - Self::atm(&values[..row.end]));
+                }
+            }
+            equations.extend(anchors);
+            Ok(equations)
+        })?;
+        let equations = residuals
+            .iter()
+            .copied()
+            .chain(std::iter::repeat_n(0.0, self.atm_count))
+            .collect();
+        let mut solution = problem.linearize(VolLinearization {
+            theta: &jac,
+            quote,
+            equations,
+            residuals,
+            iterations,
+            hessian: None,
+        })?;
+        solution.vol_butterfly_violation = self.butterfly(&market, &point, values)?;
+        Ok((point, solution))
+    }
+
+    fn butterfly(
+        &self,
+        market: &Market<f64>,
+        theta: &[f64],
+        values: &[CouponValue<f64>],
+    ) -> Result<bool> {
+        for (i, value) in values.iter().enumerate() {
+            if value.sqrt_time == 0.0 {
+                continue;
+            }
+            let atm = self
+                .rows
+                .iter()
+                .find(|row| row.atm_strike.is_some() && row.end > i)
+                .or_else(|| self.rows.iter().rev().find(|row| row.atm_strike.is_some()))
+                .and_then(|row| row.atm_strike)
+                .map(|k| theta[k]);
+            let points = market.with_view(|view| {
+                let face = view.volatility(&self.id)?;
+                let point = view.bind_vol(face, &self.points[i])?;
+                self.absolute_strikes
+                    .iter()
+                    .copied()
+                    .chain(atm)
+                    .map(|k| Ok((k, view.vol_at(point, k)? * value.sqrt_time)))
+                    .collect::<Result<Vec<_>>>()
+            })?;
+            if butterfly(self.kind, value.forward, &points)? {
+                return Ok(true);
+            }
+        }
+        Ok(false)
+    }
+}
+
+// 不可能な価格とベガの飽和を、無制限な外挿にせず校正の失敗にする。
+fn positive_vol(f: &impl Fn(f64) -> Result<f64>, price: f64, seed: f64) -> Result<(f64, usize)> {
+    if price < 0.0 {
+        return Err(Error::Calibration(
+            "cap price is below the value of already stripped coupons".into(),
+        ));
+    }
+    if price == 0.0 {
+        return Ok((0.0, 1));
+    }
+    let mut upper = (seed * 2.0).max(1e-8);
+    for i in 1..=64 {
+        if f(upper)? >= price {
+            let calls = Cell::new(i);
+            let residual = |sigma| {
+                calls.set(calls.get() + 1);
+                Ok(f(sigma)? - price)
+            };
+            let sigma = root(&residual, 0.0, upper)?;
+            return Ok((sigma, calls.get()));
+        }
+        upper *= 2.0;
+    }
+    Err(Error::Calibration(
+        "cap price lies outside the option premium range".into(),
+    ))
+}
+
+pub(super) fn root(f: &impl Fn(f64) -> Result<f64>, lo: f64, hi: f64) -> Result<f64> {
+    let error = Cell::new(None);
+    let function = ScalarFn::new(|x| forward_error(f(x), &error));
+    let options = ScalarRootOptions::default()
+        .with_step_tolerance(StepTolerance::new(1e-17, 1e-15).map_err(calibration_error)?)
+        .with_residual_tolerance(ResidualTolerance::new(0.0, 0.0).map_err(calibration_error)?);
+    Brent
+        .solve(
+            &function,
+            Interval::new(lo, hi).map_err(calibration_error)?,
+            &options,
+        )
+        .map(|r| r.point())
+        .map_err(|other| error.take().unwrap_or_else(|| calibration_error(other)))
+}
+
+pub(super) fn butterfly(kind: VolKind, forward: f64, points: &[(f64, f64)]) -> Result<bool> {
+    let mut calls = points
+        .iter()
+        .map(|&(k, sd)| {
+            Ok((
+                k,
+                option_time_value(kind, forward, k, sd)? + (forward - k).max(0.0),
+            ))
+        })
+        .collect::<Result<Vec<_>>>()?;
+    calls.sort_by(|a, b| a.0.total_cmp(&b.0));
+    calls.dedup_by(|a, b| (a.0 + 0.0).total_cmp(&(b.0 + 0.0)).is_eq());
+    Ok(calls.windows(3).any(|p| {
+        let width = p[2].0 - p[0].0;
+        let weight = (p[1].0 - p[0].0) / width;
+        let chord = p[0].1 + (p[2].1 - p[0].1) * weight;
+        p[1].1 - chord > 1e-14_f64.max(1e-12 * p[1].1.abs())
+    }))
+}
diff --git a/crates/rhodium-pricing/src/calibration/vol/swaption.rs b/crates/rhodium-pricing/src/calibration/vol/swaption.rs
new file mode 100644
index 0000000000000000000000000000000000000000..c2bada1b994c7a99de84f2656559af0220c0f20f
--- /dev/null
+++ b/crates/rhodium-pricing/src/calibration/vol/swaption.rs
@@ -0,0 +1,672 @@
+//! スワップションの SABR ノードの準備、最小二乗と ATM 制約。
+
+use std::{
+    cell::Cell,
+    collections::{BTreeMap, BTreeSet},
+    sync::Arc,
+};
+
+use rhodium_core::{
+    Error, Result, Warning,
+    calendar::{AdvanceRule, BusinessDayConvention, Calendar, EndOfMonthRule, money_market_period},
+    currency::Currency,
+    date::Date,
+    quantity::VolKind,
+    tenor::{DateOffset, Tenor, TenorLength, TenorUnit},
+    time::{LocalEventTime, model_time},
+};
+use rhodium_market::{
+    market::Market,
+    structure::{BlockId, MarketStructure},
+    vol::{Sabr, SabrFit, VolSpec},
+};
+use scioxide_ad::{Dual, Mask, Real};
+use scioxide_core::Budget;
+use scioxide_linalg::Matrix;
+use scioxide_optimize::{LeastSquares, LeastSquaresFn, LeastSquaresOptions};
+
+use crate::{
+    calibration::HessianMode,
+    context::Context,
+    definition::QuoteTransform,
+    definition::vol::{
+        RateStrike, SabrCalibration, SwaptionVolQuote, VolatilityDefinition, VolatilityMarket,
+    },
+    par,
+    prepare::{PrepareInput, Prepared, evaluate, prepare},
+};
+
+use super::{
+    GroupSolution, ResolvedSurface, UnboundGroup, UnboundModel, VolLinearization, VolProblem,
+    VolQuote, calibration_error, differentiate, forward_error,
+};
+
+pub(super) struct UnboundSwaption {
+    surface: VolatilityMarket,
+    time: f64,
+    exercise_date: Date,
+    tenor: TenorLength,
+    sabr: Sabr,
+    fit: SabrFit,
+    strikes: Vec<RateStrike>,
+}
+
+#[derive(Debug)]
+pub(super) struct SwaptionNode {
+    zero: Prepared,
+    one: Prepared,
+    ccy: Currency,
+    time: f64,
+    sabr: Sabr,
+    fit: SabrFit,
+    strikes: Vec<RateStrike>,
+    atm: Option<usize>,
+}
+
+struct FitEquations<R> {
+    equations: Vec<R>,
+    quote_jacobian: Vec<Vec<R>>,
+    residuals: Vec<R>,
+}
+
+pub(super) fn resolve(
+    surface: &VolatilityMarket,
+    expiry_time: &LocalEventTime,
+    calibration: SabrCalibration,
+    quotes: &[SwaptionVolQuote],
+    trade_date: Date,
+    calendar: &Calendar,
+    ctx: &Context<'_>,
+) -> Result<ResolvedSurface> {
+    let zone = ctx.refdata.zone(&expiry_time.zone)?;
+    let sabr = VolatilityDefinition::sabr(surface, calibration)?;
+    let mut warnings = Vec::new();
+    let mut expiries = BTreeMap::new();
+    let mut tenors = BTreeMap::new();
+    let mut quotes_by_node = BTreeMap::<_, Vec<_>>::new();
+    let maturity = |tenor| {
+        money_market_period(
+            trade_date,
+            tenor,
+            DateOffset::new(0, TenorUnit::BusinessDay),
+            &AdvanceRule {
+                end_of_month: EndOfMonthRule::Ignore,
+                bdc: BusinessDayConvention::ModifiedFollowing,
+                calendar,
+            },
+        )
+    };
+    for quote in quotes {
+        let tenor = Tenor::Length(quote.tenor);
+        let expiry = Tenor::Length(quote.expiry);
+        if let std::collections::btree_map::Entry::Vacant(entry) = expiries.entry(expiry) {
+            let (date, notices) = maturity(expiry)?;
+            warnings.extend(notices);
+            entry.insert((
+                zone.resolve_convention(expiry_time.on_date(date.end()))?,
+                date.end(),
+            ));
+        }
+        if let std::collections::btree_map::Entry::Vacant(entry) = tenors.entry(tenor) {
+            let (dates, notices) = maturity(tenor)?;
+            warnings.extend(notices);
+            entry.insert(f64::from(dates.end().days_since(dates.start())) / 365.0);
+        }
+        quotes_by_node
+            .entry((expiry, tenor))
+            .or_default()
+            .push(quote);
+    }
+    let mut expiry_axis: Vec<_> = expiries.into_iter().collect();
+    expiry_axis.sort_by_key(|(_, (at, _))| *at);
+    let mut tenor_axis: Vec<_> = tenors.into_iter().collect();
+    tenor_axis.sort_by(|(_, a), (_, b)| a.total_cmp(b));
+    let spec = VolSpec::swaption_cube(
+        surface.id.clone(),
+        surface.index.clone(),
+        sabr,
+        calibration.fit,
+        calibration.interpolation,
+        expiry_axis.iter().map(|(_, (at, _))| *at).collect(),
+        tenor_axis.iter().map(|(_, years)| *years).collect(),
+    )?;
+    let mut groups = Vec::with_capacity(quotes_by_node.len());
+    for (e, (expiry, (at, date))) in expiry_axis.into_iter().enumerate() {
+        for (t, &(tenor, _)) in tenor_axis.iter().enumerate() {
+            let quotes = quotes_by_node.remove(&(expiry, tenor)).ok_or_else(|| {
+                Error::InvalidData(format!("missing SABR cube node {expiry} {tenor}"))
+            })?;
+            validate_node(calibration.fit, &quotes)?;
+            let Tenor::Length(length) = tenor else {
+                return Err(Error::Internal("expected swap length".into()));
+            };
+            groups.push(UnboundGroup {
+                block: VolSpec::swaption_block(&surface.id, e, t)?,
+                model: UnboundModel::Swaption(UnboundSwaption {
+                    surface: surface.clone(),
+                    time: model_time(ctx.as_of, at).years(),
+                    exercise_date: date,
+                    tenor: length,
+                    sabr,
+                    fit: calibration.fit,
+                    strikes: quotes.iter().map(|row| row.strike).collect(),
+                }),
+                quotes: quotes
+                    .iter()
+                    .map(|row| VolQuote {
+                        id: row.quote.clone(),
+                        transform: row
+                            .atm_quote
+                            .clone()
+                            .map_or(QuoteTransform::Identity, QuoteTransform::AtmPlusOffset),
+                        tolerance: row.tolerance,
+                    })
+                    .collect(),
+            });
+        }
+    }
+    Ok(ResolvedSurface {
+        spec,
+        groups,
+        warnings,
+    })
+}
+
+fn validate_node(fit: SabrFit, rows: &[&SwaptionVolQuote]) -> Result<()> {
+    let mut coordinates = BTreeSet::new();
+    let atm = rows
+        .iter()
+        .find(|row| matches!(row.strike, RateStrike::Atm));
+    let valid = match fit {
+        SabrFit::AtmOnly { .. } => rows.len() == 1 && atm.is_some(),
+        SabrFit::LeastSquares => rows.len() >= 3,
+        SabrFit::ExactAtm => rows.len() >= 3 && atm.is_some(),
+    };
+    if !valid {
+        return Err(Error::InvalidData(format!(
+            "not enough SABR quotes for {fit:?}"
+        )));
+    }
+    for row in rows {
+        let coord = match row.strike {
+            RateStrike::Atm => (0, 0),
+            RateStrike::AtmOffset(x) => (1, x.to_bits()),
+            RateStrike::Absolute(x) => (2, if x == 0.0 { 0 } else { x.to_bits() }),
+        };
+        if !coordinates.insert(coord)
+            || row.atm_quote.as_ref().is_some_and(|id| {
+                atm.is_none_or(|base| base.quote != *id || base.quote == row.quote)
+            })
+        {
+            return Err(Error::InvalidData(format!(
+                "duplicate strike or invalid ATM quote reference for {}",
+                row.quote
+            )));
+        }
+    }
+    Ok(())
+}
+
+impl UnboundSwaption {
+    pub(super) fn prepare(
+        self,
+        structure: &Arc<MarketStructure>,
+        ctx: &Context<'_>,
+    ) -> Result<(SwaptionNode, Vec<Warning>)> {
+        let template = self.surface.template(self.tenor, ctx.refdata)?;
+        let ccy = template.currency(ctx.refdata)?;
+        let at = |rate| {
+            let product = template.build_for_trade(rate, self.exercise_date, ctx.refdata)?;
+            prepare(
+                PrepareInput::Standard {
+                    product: &product,
+                    trade_date: self.exercise_date,
+                },
+                &self.surface.collateral,
+                &[],
+                structure,
+                ctx,
+            )
+        };
+        let zero = at(0.0)?;
+        let one = at(1.0)?;
+        let warnings = zero
+            .warnings()
+            .iter()
+            .chain(one.warnings())
+            .cloned()
+            .collect();
+        let atm = self
+            .strikes
+            .iter()
+            .position(|strike| matches!(strike, RateStrike::Atm));
+        Ok((
+            SwaptionNode {
+                zero,
+                one,
+                ccy,
+                time: self.time,
+                sabr: self.sabr,
+                fit: self.fit,
+                strikes: self.strikes,
+                atm,
+            },
+            warnings,
+        ))
+    }
+}
+
+impl SwaptionNode {
+    pub(super) fn blocks(&self) -> BTreeSet<BlockId> {
+        self.zero
+            .blocks()
+            .iter()
+            .chain(self.one.blocks())
+            .cloned()
+            .collect()
+    }
+
+    fn forward<R: Real>(&self, market: &Market<R>) -> Result<R> {
+        par::affine(
+            evaluate(&self.zero, market)?.get(self.ccy),
+            evaluate(&self.one, market)?.get(self.ccy),
+        )
+    }
+
+    fn vols<R: Real>(&self, local: &[R]) -> Result<Vec<R>> {
+        let forward = local[self.fit.dimension()];
+        let coords = if self.fit.dimension() == 1 {
+            [local[0], R::zero(), R::zero()]
+        } else {
+            [local[0], local[1], local[2]]
+        };
+        let params = self.fit.parameters(coords);
+        self.strikes
+            .iter()
+            .map(|strike| {
+                self.sabr
+                    .volatility(forward, strike.value(forward), self.time, params)
+            })
+            .collect()
+    }
+
+    fn atm_root(&self, forward: f64, target: f64, rho: f64, nu: f64) -> Result<f64> {
+        let f = forward + self.sabr.shift();
+        let beta = self.sabr.beta();
+        let normal = matches!(self.sabr.kind(), VolKind::Normal);
+        let scaled_target = if normal { target / f } else { target };
+        let a = (if normal {
+            beta * (beta - 2.0)
+        } else {
+            (1.0 - beta).powi(2)
+        }) * self.time
+            / 24.0;
+        let b = rho * beta * nu * self.time / 4.0;
+        let c = 1.0 + (2.0 - 3.0 * rho * rho) * nu * nu * self.time / 24.0;
+        let root = positive_root(a, b, c, scaled_target)?;
+        Ok((root * f.powf(1.0 - beta)).ln())
+    }
+
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#201 で各 variant を網羅する"
+    )]
+    fn project(&self, forward: f64, target: f64, x: &[f64]) -> Result<Vec<f64>> {
+        match self.fit {
+            SabrFit::ExactAtm => Ok(vec![
+                self.atm_root(forward, target, x[0].tanh(), x[1].exp())?,
+                x[0],
+                x[1],
+                forward,
+            ]),
+            _ => Ok(x.iter().copied().chain([forward]).collect()),
+        }
+    }
+
+    pub(super) fn solve(&self, problem: &VolProblem<'_>) -> Result<(Vec<f64>, GroupSolution)> {
+        let forward = self.forward(&problem.market(&[] as &[f64], None)?)?;
+        if !(forward + self.sabr.shift() > 0.0 && problem.targets.iter().all(|&x| x > 0.0)) {
+            return Err(Error::Calibration(
+                "SABR requires positive target volatilities and shifted forward".into(),
+            ));
+        }
+        let mut coordinates: Vec<_> = self.strikes.iter().map(|k| k.value(forward)).collect();
+        coordinates.sort_by(f64::total_cmp);
+        if coordinates.windows(2).any(|pair| pair[0] >= pair[1]) {
+            return Err(Error::Calibration("coincident SABR strikes".into()));
+        }
+        let target = self.atm.map_or_else(
+            || {
+                let mut sorted = problem.targets.clone();
+                sorted.sort_by(f64::total_cmp);
+                sorted[sorted.len() / 2]
+            },
+            |i| problem.targets[i],
+        );
+        if let SabrFit::AtmOnly { rho, nu } = self.fit {
+            let point = vec![self.atm_root(forward, target, rho, nu)?, forward];
+            return self.finish(problem, point, &[1.0], 0);
+        }
+        // 正規化して LM の勾配停止判定を入力のボラ単位に左右させない。
+        let weights: Vec<_> = problem
+            .rows
+            .iter()
+            .map(|row| row.tolerance.map_or(0.0, f64::recip))
+            .collect();
+        let mut best = None::<(Vec<f64>, GroupSolution)>;
+        let mut failure = None;
+        for (rho, nu) in [(0.0_f64, 0.30_f64), (-0.40, 0.65), (0.40, 0.65)] {
+            let seed = [rho.atanh(), nu.ln()];
+            let result = (|| {
+                let start = if matches!(self.fit, SabrFit::ExactAtm) {
+                    seed.to_vec()
+                } else {
+                    let alpha = self.atm_root(forward, target, rho, nu).unwrap_or_else(|_| {
+                        let f = forward + self.sabr.shift();
+                        let exponent = if matches!(self.sabr.kind(), VolKind::Normal) {
+                            -self.sabr.beta()
+                        } else {
+                            1.0 - self.sabr.beta()
+                        };
+                        (target * f.powf(exponent)).ln()
+                    });
+                    vec![alpha, seed[0], seed[1]]
+                };
+                let (point, iterations) =
+                    self.fit_from(problem, &weights, forward, target, &start)?;
+                self.finish(problem, point, &weights, iterations)
+            })();
+            match result {
+                Ok((point, solution)) => {
+                    let cost = |solution: &GroupSolution| {
+                        solution
+                            .residuals()
+                            .iter()
+                            .zip(&weights)
+                            .map(|(r, w)| (r * w).powi(2))
+                            .sum::<f64>()
+                    };
+                    if best
+                        .as_ref()
+                        .is_none_or(|(_, best)| cost(&solution) < cost(best))
+                    {
+                        best = Some((point, solution));
+                    }
+                }
+                Err(error) => failure = Some(error),
+            }
+        }
+        best.ok_or_else(|| {
+            failure.unwrap_or_else(|| Error::Calibration("SABR did not converge".into()))
+        })
+    }
+
+    fn fit_from(
+        &self,
+        problem: &VolProblem<'_>,
+        weights: &[f64],
+        forward: f64,
+        target: f64,
+        start: &[f64],
+    ) -> Result<(Vec<f64>, usize)> {
+        let active: Vec<_> = weights
+            .iter()
+            .enumerate()
+            .filter_map(|(i, &w)| (w > 0.0).then_some(i))
+            .collect();
+        let error = Cell::new(None);
+        let function =
+            LeastSquaresFn::new(start.len(), active.len(), |x: &[f64], out: &mut [f64]| {
+                let local = forward_error(self.project(forward, target, x), &error)?;
+                let vols = forward_error(self.vols(&local), &error)?;
+                for (out, &i) in out.iter_mut().zip(&active) {
+                    *out = (vols[i] - problem.targets[i]) * weights[i];
+                }
+                Ok(())
+            })
+            .map_err(calibration_error)?
+            .with_jacobian(|x: &[f64], out: &mut Matrix<f64>| {
+                let local = forward_error(self.project(forward, target, x), &error)?;
+                let partials = forward_error(self.partials(&local), &error)?;
+                for (row, &i) in active.iter().enumerate() {
+                    let reduced = self.reduced(&partials, i);
+                    for (j, d) in reduced.into_iter().enumerate() {
+                        out.as_mut_slice()[row * start.len() + j] = d * weights[i];
+                    }
+                }
+                Ok(())
+            });
+        let solver = LeastSquares::new(
+            LeastSquaresOptions::default().with_max_iterations(Budget::ONE.saturating_add(199)),
+        )
+        .map_err(calibration_error)?;
+        let result = solver
+            .solve(&function, start)
+            .map_err(|other| error.take().unwrap_or_else(|| calibration_error(other)))?;
+        Ok((
+            self.project(forward, target, result.x())?,
+            result.iterations(),
+        ))
+    }
+
+    // 列を Dual で作る。外側の Var から停留条件を微分すると厳密な二次項も入る。
+    fn partials<R: Real>(&self, local: &[R]) -> Result<Vec<Vec<R>>> {
+        (0..self.fit.dimension())
+            .map(|j| {
+                let dual: Vec<_> = local
+                    .iter()
+                    .enumerate()
+                    .map(|(i, &x)| Dual::new(x, if i == j { R::one() } else { R::zero() }))
+                    .collect();
+                Ok(self
+                    .vols(&dual)?
+                    .into_iter()
+                    .map(Dual::derivative)
+                    .collect())
+            })
+            .collect()
+    }
+
+    fn reduced<R: Real>(&self, partials: &[Vec<R>], row: usize) -> Vec<R> {
+        if let (SabrFit::ExactAtm, Some(atm)) = (self.fit, self.atm) {
+            let ratio = partials[0][row] / partials[0][atm];
+            partials[1..]
+                .iter()
+                .map(|col| col[row] - ratio * col[atm])
+                .collect()
+        } else {
+            partials.iter().map(|col| col[row]).collect()
+        }
+    }
+
+    fn equations<R: Real>(
+        &self,
+        local: &[R],
+        weights: &[f64],
+        targets: &[f64],
+        forward: R,
+    ) -> Result<FitEquations<R>> {
+        let residuals: Vec<_> = self
+            .vols(local)?
+            .into_iter()
+            .zip(targets)
+            .map(|(m, &q)| m - q)
+            .collect();
+        let dimension = self.fit.dimension();
+        let mut d = vec![vec![R::zero(); targets.len()]; dimension + 1];
+        match self.fit {
+            SabrFit::AtmOnly { .. } => d[0][0] = -R::one(),
+            SabrFit::LeastSquares | SabrFit::ExactAtm => {
+                let partials = self.partials(local)?;
+                let offset = usize::from(matches!(self.fit, SabrFit::ExactAtm));
+                if offset != 0 {
+                    let atm = self
+                        .atm
+                        .ok_or_else(|| Error::Internal("ATM constraint lacks a quote".into()))?;
+                    if !partials[0][atm].gt(R::zero()).all() {
+                        return Err(Error::Calibration("singular ATM constraint".into()));
+                    }
+                    d[0][atm] = -R::one();
+                }
+                for (i, &weight) in weights.iter().enumerate().filter(|(_, w)| **w > 0.0) {
+                    for (j, value) in self.reduced(&partials, i).into_iter().enumerate() {
+                        d[j + offset][i] = -value * (weight * weight);
+                    }
+                }
+            }
+        }
+        let mut equations: Vec<_> = d[..dimension]
+            .iter()
+            .map(|row| {
+                row.iter()
+                    .zip(&residuals)
+                    .fold(R::zero(), |g, (&d, &r)| g - d * r)
+            })
+            .collect();
+        equations.push(local[dimension] - forward);
+        Ok(FitEquations {
+            equations,
+            quote_jacobian: d,
+            residuals,
+        })
+    }
+
+    fn require_minimum(&self, jac: &Matrix<f64>) -> Result<()> {
+        let exact = usize::from(matches!(self.fit, SabrFit::ExactAtm));
+        if matches!(self.fit, SabrFit::AtmOnly { .. }) {
+            return Ok(());
+        }
+        let dimension = self.fit.dimension();
+        let row = |i| &jac.as_slice()[i * jac.cols()..(i + 1) * jac.cols()];
+        // ExactAtm のヘシアンは、停留条件の微分を ATM 制約の接空間へ制限する。
+        let mut entries = Vec::new();
+        for i in exact..dimension {
+            for j in exact..dimension {
+                let reduced = if exact == 0 {
+                    row(i)[j]
+                } else {
+                    row(i)[j] - row(i)[0] * row(0)[j] / row(0)[0]
+                };
+                entries.push(reduced);
+            }
+        }
+        let size = dimension - exact;
+        let matrix = Matrix::from_vec(size, size, entries).map_err(calibration_error)?;
+        crate::calibration::solver::require_positive_definite(&matrix)
+    }
+
+    fn finish(
+        &self,
+        problem: &VolProblem<'_>,
+        mut point: Vec<f64>,
+        weights: &[f64],
+        iterations: usize,
+    ) -> Result<(Vec<f64>, GroupSolution)> {
+        let n = point.len();
+        let forward = self.forward(&problem.market(&[], None)?)?;
+        for step in 0..=8 {
+            let FitEquations {
+                equations,
+                quote_jacobian: quote,
+                residuals,
+            } = self.equations(&point, weights, &problem.targets, forward)?;
+            let inputs = problem.point(&point);
+            let jac = differentiate(&inputs, |x| {
+                let (local, upstream) = x.split_at(n);
+                let forward = self.forward(&problem.market(local, Some(upstream))?)?;
+                self.equations(local, weights, &problem.targets, forward)
+                    .map(|x| x.equations)
+            })?;
+            let mut solution = problem.linearize(VolLinearization {
+                theta: &jac,
+                quote,
+                equations,
+                residuals,
+                iterations: iterations + step,
+                hessian: (!matches!(self.fit, SabrFit::AtmOnly { .. }))
+                    .then_some(HessianMode::Exact),
+            })?;
+            let vols = self.vols(&point)?;
+            let movement = self
+                .vols(
+                    &point
+                        .iter()
+                        .zip(solution.theta_error())
+                        .map(|(x, d)| x - d)
+                        .collect::<Vec<_>>(),
+                )?
+                .iter()
+                .zip(&vols)
+                .map(|(a, b)| (a - b).powi(2))
+                .sum::<f64>()
+                .sqrt();
+            if movement <= problem.tolerance
+                && (!matches!(self.fit, SabrFit::AtmOnly { .. } | SabrFit::ExactAtm)
+                    || solution.residuals()[self.atm.unwrap_or(0)].abs() <= problem.tolerance)
+            {
+                self.require_minimum(&jac)?;
+                solution.vol_butterfly_violation = super::caplet::butterfly(
+                    self.sabr.kind(),
+                    forward,
+                    &self
+                        .strikes
+                        .iter()
+                        .zip(vols)
+                        .map(|(strike, vol)| (strike.value(forward), vol * self.time.sqrt()))
+                        .collect::<Vec<_>>(),
+                )?;
+                return Ok((point, solution));
+            }
+            for (x, d) in point.iter_mut().zip(solution.theta_error()) {
+                *x -= d;
+            }
+        }
+        Err(Error::Calibration(
+            "SABR stationarity tolerance not reached".into(),
+        ))
+    }
+}
+
+// ATM は x = α / F^(1-β) の三次式。極値で区切り、最初の増加する根を選ぶ。
+fn positive_root(a: f64, b: f64, c: f64, target: f64) -> Result<f64> {
+    let f = |x: f64| ((a * x + b) * x + c) * x - target;
+    let b = 2.0 * b;
+    let a = 3.0 * a;
+    let mut edges = Vec::new();
+    if a == 0.0 {
+        if b != 0.0 {
+            edges.push(-c / b);
+        }
+    } else {
+        let disc = b * b - 4.0 * a * c;
+        if disc >= 0.0 {
+            let q = -0.5 * (b + disc.sqrt().copysign(b));
+            edges.push(q / a);
+            if q != 0.0 {
+                edges.push(c / q);
+            }
+        }
+    }
+    edges.retain(|&x| x.is_finite() && x > 0.0);
+    edges.sort_by(f64::total_cmp);
+    let mut low: f64 = 0.0;
+    for high in edges.iter().copied().chain(std::iter::once(
+        low.max(edges.last().copied().unwrap_or(0.0)).max(target) * 2.0,
+    )) {
+        if f(low) < 0.0 && f(high) > 0.0 {
+            return super::caplet::root(&|x| Ok(f(x)), low, high);
+        }
+        low = high;
+    }
+    // 最後の極値より先が増加する場合だけ、根を挟むまで区間を倍にする。
+    for _ in 0..64 {
+        let high = low * 2.0;
+        if f(low) < 0.0 && f(high) > 0.0 {
+            return super::caplet::root(&|x| Ok(f(x)), low, high);
+        }
+        low = high;
+    }
+    Err(Error::Calibration(
+        "SABR ATM has no identified positive root".into(),
+    ))
+}
diff --git a/crates/rhodium-pricing/src/config.rs b/crates/rhodium-pricing/src/config.rs
index b5edc482a2685afda604386af6310da0127dd7a7..5c9b2b0283258493666d98f82e617357d453f1ad 100644
--- a/crates/rhodium-pricing/src/config.rs
+++ b/crates/rhodium-pricing/src/config.rs
@@ -10,7 +10,10 @@ use rhodium_core::{
     Error, Result, Warning, currency::Currency, id::IndexId, refdata::RefData, time::LocalEventTime,
 };
 use rhodium_market::curve::CurveId;
-use rhodium_products::{ir::SettlementBoundary, standard::ProductKind};
+use rhodium_products::{
+    ir::{SettlementBoundary, analysis::Linearity},
+    standard::ProductKind,
+};
 use serde::{Deserialize, Deserializer, de};
 
 use crate::{definition::MarketDefinition, resolve::CollateralKey};
@@ -66,7 +69,7 @@ pub enum AnalyticPricer {
     RateFuture,
 }
 
-/// M0 のリニア評価で使える近似。
+/// リニア評価で使える名前付きの近似（設計書の「近似の規約」）。
 #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deserialize)]
 #[serde(rename_all = "snake_case")]
 pub enum LinearConvention {
@@ -366,8 +369,8 @@ impl PricingConfig {
 
     /// 方式の明示行を優先し、指定がなければ線形な IR をリニア方式で評価する。
     ///
-    /// `product` は標準商品の種類。呼び出し側は簡約済みの IR から、リニア方式が
-    /// 適用できるかを判定して渡す。
+    /// `product` は標準商品の種類。`linearity` は簡約済みの IR にリニア方式を
+    /// 適用できるかで、呼び出し側が判定して渡す。
     ///
     /// # Errors
     /// 選んだ方式を適用できない場合や、適用できる方式が指定されていない場合は
@@ -375,14 +378,14 @@ impl PricingConfig {
     pub fn select_method(
         &self,
         product: Option<ProductKind>,
-        linear_supported: bool,
+        linearity: Linearity,
     ) -> Result<Method> {
         let explicit = self
             .0
             .methods
             .iter()
             .find(|rule| Some(rule.product) == product);
-        match (explicit.map(|rule| rule.method), product, linear_supported) {
+        match (explicit.map(|rule| rule.method), product, linearity) {
             (
                 Some(Method::Analytic(AnalyticPricer::RateFuture)) | None,
                 Some(ProductKind::RateFuture),
@@ -391,11 +394,11 @@ impl PricingConfig {
             (Some(Method::Analytic(_)), _, _) => Err(Error::Unsupported(
                 "analytic pricer does not support this product".into(),
             )),
-            (Some(Method::Linear) | None, _, true) => Ok(Method::Linear),
-            (Some(Method::Linear), _, false) => Err(Error::Unsupported(
+            (Some(Method::Linear) | None, _, Linearity::Linear) => Ok(Method::Linear),
+            (Some(Method::Linear), _, Linearity::NonLinear(_)) => Err(Error::Unsupported(
                 "configured linear method cannot evaluate this product".into(),
             )),
-            (None, _, false) => Err(Error::Unsupported(
+            (None, _, Linearity::NonLinear(_)) => Err(Error::Unsupported(
                 "the default linear method cannot evaluate this product".into(),
             )),
         }
diff --git a/crates/rhodium-pricing/src/definition.rs b/crates/rhodium-pricing/src/definition.rs
index 63da3993efd962b7b924131c41f076681588865a..20d3b9c93b5fad68422650a69630afa1273741ba 100644
--- a/crates/rhodium-pricing/src/definition.rs
+++ b/crates/rhodium-pricing/src/definition.rs
@@ -1,27 +1,28 @@
 //! 市場の構造と校正 quote の定義。
 
 mod fx;
-pub use fx::{
-    DirectFx, FxDatedQuote, FxDefinition, FxImpliedCurve, FxQuotes, FxSettlement, FxSpotLag,
-};
+pub mod vol;
+pub(crate) use fx::ResolvedFxDefinition;
+pub use fx::{DirectFx, FxDatedQuote, FxDefinition, FxImpliedCurve, FxImpliedLongEnd, FxQuotes};
 
-use std::collections::BTreeMap;
+use std::collections::{BTreeMap, BTreeSet};
 
 use rhodium_core::{
     Error, Result,
     currency::Currency,
     id::CalendarId,
     quantity::YearFraction,
+    refdata::{Convention, RefData},
     time::{Instant, TzId},
 };
 use rhodium_market::{
-    curve::{CurveId, CurveSpec},
+    curve::{CurveId, CurveRole, CurveSpec},
     quote::{QuoteId, QuoteKind},
     structure::{BlockId, ParamId},
 };
 use rhodium_products::standard::InstrumentTemplate;
 
-use crate::resolve::CollateralKey;
+use crate::{config::PricingConfig, resolve::CollateralKey};
 
 /// 市場から測る校正目標。
 #[derive(Debug, Clone)]
@@ -45,24 +46,26 @@ impl CalibrationTarget {
 
 /// 入力 quote を目標の単位に変換する規則 q̃。
 #[derive(Debug, Clone, Default, PartialEq, serde::Deserialize)]
-#[serde(rename_all = "snake_case")]
+#[serde(rename_all = "snake_case", deny_unknown_fields)]
 pub enum QuoteTransform {
     /// 入力値をそのまま使う。
     #[default]
     Identity,
+    /// ボラの差分に、同じノードの ATM quote を足す。
+    AtmPlusOffset(QuoteId),
     /// 期末のスプレッド s に受渡期間の日数係数 `τ_turn` を掛け、ln DF の跳び `s·τ_turn` にする。
     ///
-    /// `τ_turn` は期末の営業日から翌営業日までを参照データで解決した値で、呼び出し側が渡す。
-    /// ファイルの binding からは読まない。
-    #[serde(skip)]
+    /// `τ_turn` は期末の営業日から翌営業日までの日数係数で、M1 では呼び出し側（ファイルでは
+    /// `{ turn_spread_to_jump = { value = 0.0164, basis = "ACT365F" } }`）が明示する。
+    /// 参照データの期末の規則から解決する書き方は、ファイル形式の見直し（M6）で足す。
     TurnSpreadToJump(YearFraction),
 }
 
 impl QuoteTransform {
-    /// 入力 quote q に掛けて目標の値 q̃ = D·q にする係数 D。どの変換も quote について線形。
+    /// 変換後の値に対する、この行が指定する quote の係数。すべての変換は quote について線形。
     pub(crate) fn scale(&self) -> f64 {
         match self {
-            Self::Identity => 1.0,
+            Self::Identity | Self::AtmPlusOffset(_) => 1.0,
             Self::TurnSpreadToJump(fraction) => fraction.value(),
         }
     }
@@ -95,6 +98,7 @@ pub struct InstrumentBinding {
     pub(crate) determines: ElementSelector,
     kind: QuoteKind,
     transform: QuoteTransform,
+    tolerance: Option<f64>,
 }
 
 impl InstrumentBinding {
@@ -176,6 +180,7 @@ impl InstrumentBinding {
                 determines,
                 kind,
                 transform,
+                tolerance: None,
             })
         } else {
             Err(Error::InvalidInput(format!(
@@ -201,6 +206,74 @@ impl InstrumentBinding {
     pub const fn transform(&self) -> &QuoteTransform {
         &self.transform
     }
+
+    /// 最小二乗の当てはめの許容幅 σ を quote の単位で与える。
+    ///
+    /// 許容幅を持つ binding だけからなるグループは、重み W = diag(1/σ̃²) の最小二乗で解く
+    /// （σ̃ = σ·D は残差の単位の幅で、D は quote の変換の係数）。許容幅を持たない binding の
+    /// グループは正方系で解き、未知数と残差の数の一致を求める。
+    ///
+    /// # Errors
+    ///
+    /// 許容幅が有限な正の数でなければ `InvalidInput`。
+    pub fn with_tolerance(mut self, tolerance: f64) -> Result<Self> {
+        if !(tolerance.is_finite() && tolerance > 0.0) {
+            return Err(Error::InvalidInput(format!(
+                "calibration tolerance of {} must be finite and positive",
+                self.quote
+            )));
+        }
+        self.tolerance = Some(tolerance);
+        Ok(self)
+    }
+
+    /// 最小二乗の当てはめの許容幅 σ（quote の単位）。正方系の binding は `None`。
+    #[must_use]
+    pub const fn tolerance(&self) -> Option<f64> {
+        self.tolerance
+    }
+}
+
+/// 会合日の階段（`Steps` の部品）の隣り合う区間のフォワードの差に掛ける罰則。
+///
+/// 罰則 `P(θ) = (λ/2) Σ_i ((f_{i+1} − f_i)/s)²` を最小二乗の目的に足す（λ は重み、s は参照の幅）。
+/// 定型のテナーの OIS だけで階段を当てると区間ごとの水準が決まらないので、
+/// 罰則が隣り合う区間を近づけて解を一意にする。罰則だけで決まる水準は政策金利の織り込みとして報告しない。
+#[derive(Debug, Clone, PartialEq)]
+pub struct StepPenalty {
+    pub(crate) block: BlockId,
+    pub(crate) weight: f64,
+    pub(crate) scale: f64,
+}
+
+impl StepPenalty {
+    /// 罰則を掛けるブロック、重み λ と参照の幅 s を指定する。
+    ///
+    /// # Errors
+    ///
+    /// λ か s が有限な正の数でなければ `InvalidInput`。
+    pub fn new(block: BlockId, weight: f64, scale: f64) -> Result<Self> {
+        if !([weight, scale]
+            .iter()
+            .all(|value| value.is_finite() && *value > 0.0))
+        {
+            return Err(Error::InvalidInput(format!(
+                "step penalty of {} {} needs a finite positive weight and scale",
+                block.curve(),
+                block.part()
+            )));
+        }
+        Ok(Self {
+            block,
+            weight,
+            scale,
+        })
+    }
+
+    /// 罰則の行 `√λ·(f_{i+1} − f_i)/s` の係数 `√λ/s`。
+    pub(crate) fn factor(&self) -> f64 {
+        self.weight.sqrt() / self.scale
+    }
 }
 
 /// 曲線の構造と校正条件。参照データと評価時点は校正時に受け取る。
@@ -211,10 +284,119 @@ pub struct MarketDefinition {
     pub(crate) curves: Vec<CurveSpec>,
     pub(crate) bindings: Vec<InstrumentBinding>,
     pub(crate) fx: Vec<FxDefinition>,
+    pub(crate) volatilities: Vec<vol::VolatilityDefinition>,
     pub(crate) lgm: BTreeMap<Currency, LgmDefinition>,
+    pub(crate) penalties: BTreeMap<BlockId, StepPenalty>,
+    pub(crate) gauss_newton: BTreeSet<BlockId>,
     quote_kinds: BTreeMap<QuoteId, QuoteKind>,
 }
 
+impl MarketDefinition {
+    /// 予想曲線ごとに、その予想を定める担保の割引曲線（予想の測度）を求める。
+    ///
+    /// 予想曲線のブロックを当てる商品の binding のうち、予想の担保を識別できるものがすべて
+    /// 同じ割引曲線を指すとき、その曲線が測度になる（[`Self::measure_evidence`]）。
+    /// 根拠の食い違い、割引の規則が無い根拠、根拠が一つも無い（直接マークだけ、付利指標の無い
+    /// 先物だけ）予想曲線は測度が分からないので返さない。
+    pub(crate) fn forecast_measures(
+        &self,
+        config: &PricingConfig,
+        refdata: &RefData,
+    ) -> Vec<(CurveId, CurveId)> {
+        self.curves
+            .iter()
+            .filter(|curve| curve.has_role(CurveRole::Projection))
+            .filter_map(|curve| {
+                let mut measures = self
+                    .bindings
+                    .iter()
+                    .filter(|binding| binding.determines.block().curve() == curve.id())
+                    .filter_map(|binding| {
+                        match Self::measure_evidence(binding, curve.currency(), config, refdata) {
+                            MeasureEvidence::Silent => None,
+                            MeasureEvidence::Discount(curve) => Some(Some(curve)),
+                            MeasureEvidence::Unknown => Some(None),
+                        }
+                    });
+                let first = measures.next()??;
+                measures
+                    .all(|measure| measure == Some(first))
+                    .then(|| (curve.id().clone(), first.clone()))
+            })
+            .collect()
+    }
+
+    /// 一つの binding が予想の測度について与える根拠。
+    ///
+    /// - 直接マークは担保を持たないので根拠にならない。
+    /// - 先物の価格は値洗いの口座の測度で決まり、取引の担保の割引因子は打ち消し合う
+    ///   （設計書の値洗いの列の評価）。同じ先物を
+    ///   どの担保で校正しても曲線は変わらないので、binding の担保は根拠にならない。値洗いの口座の
+    ///   付利指標（`margin_index`）があれば、その予想曲線が口座の測度の割引曲線を示す。
+    /// - それ以外の商品は、binding の担保で割り引く PV をパーにするので、その割引曲線が根拠になる。
+    fn measure_evidence<'c>(
+        binding: &InstrumentBinding,
+        currency: Currency,
+        config: &'c PricingConfig,
+        refdata: &RefData,
+    ) -> MeasureEvidence<'c> {
+        let found = |curve: Option<&'c CurveId>| {
+            curve.map_or(MeasureEvidence::Unknown, MeasureEvidence::Discount)
+        };
+        match binding.target.template() {
+            None => MeasureEvidence::Silent,
+            Some(InstrumentTemplate::RateFuture { convention, .. }) => {
+                match refdata.convention(convention) {
+                    Ok(Convention::RateFuture(future)) => future
+                        .margin_index
+                        .as_ref()
+                        .map_or(MeasureEvidence::Silent, |index| {
+                            found(config.forecast_curve(index))
+                        }),
+                    _ => MeasureEvidence::Unknown,
+                }
+            }
+            Some(_) => found(config.discount_curve(currency, &binding.collateral)),
+        }
+    }
+}
+
+/// 一つの binding が予想曲線の予想の測度について与える根拠。
+enum MeasureEvidence<'c> {
+    /// 予想の担保を識別しない（直接マーク、付利指標の無い先物）。
+    Silent,
+    /// この割引曲線の担保の下の予想だと示す。
+    Discount(&'c CurveId),
+    /// 担保を識別するが、割引曲線を決められない（割引や予想の規則が無い）。測度を不明にする。
+    Unknown,
+}
+
+/// 曲線のブロックのすべての要素が、当てはめか直接マークに割り当てられているかを確かめる。
+fn check_assigned(assignment: &BTreeMap<BlockId, Vec<Option<Assigned>>>) -> Result<()> {
+    if let Some((block, element)) = assignment.iter().find_map(|(block, flags)| {
+        flags
+            .iter()
+            .position(Option::is_none)
+            .map(|element| (block, element))
+    }) {
+        return Err(Error::InvalidInput(format!(
+            "unassigned element {element} of curve {} part {}",
+            block.curve(),
+            block.part()
+        )));
+    }
+    Ok(())
+}
+
+/// θ の要素への割り当て。直接マークと商品の当てはめの競合を検出する。
+#[derive(Debug, Clone, Copy)]
+enum Assigned {
+    /// 商品の残差が当てはめる。
+    Fit,
+    /// 直接マークが値を決める。`fitted` は最小二乗の許容幅を持つか。
+    Direct { fitted: bool },
+}
+
 /// 通貨別の区間一定 LGM に渡す市場マーク。
 #[derive(Debug, Clone, serde::Deserialize)]
 #[serde(deny_unknown_fields)]
@@ -277,7 +459,15 @@ impl MarketDefinition {
                     block.curve()
                 ))
             })?;
-            let direct = matches!(binding.target, CalibrationTarget::DirectMark);
+            // 直接マークは要素の値を決める。許容幅を持つ直接マークどうしだけは同じ要素を重ねて当て、
+            // 最小二乗で一つの値に決める（手計算の検算 min (θ − q₁)² + 4(θ − q₂)²）。
+            let current = if matches!(binding.target, CalibrationTarget::DirectMark) {
+                Assigned::Direct {
+                    fitted: binding.tolerance.is_some(),
+                }
+            } else {
+                Assigned::Fit
+            };
             let selected = match &binding.determines {
                 ElementSelector::Block(_) => assigned.as_mut_slice(),
                 ElementSelector::Param(param) => {
@@ -296,38 +486,104 @@ impl MarketDefinition {
                 }
             };
             for status in selected {
-                if *status == Some(true) || (direct && status.is_some()) {
-                    return Err(Error::InvalidInput(format!(
-                        "conflicting direct mark for {} {}",
-                        block.curve(),
-                        block.part()
-                    )));
+                match (*status, current) {
+                    (None, _) => *status = Some(current),
+                    (Some(Assigned::Fit), Assigned::Fit)
+                    | (
+                        Some(Assigned::Direct { fitted: true }),
+                        Assigned::Direct { fitted: true },
+                    ) => {}
+                    (Some(Assigned::Fit | Assigned::Direct { .. }), _) => {
+                        return Err(Error::InvalidInput(format!(
+                            "conflicting direct mark for {} {}",
+                            block.curve(),
+                            block.part()
+                        )));
+                    }
                 }
-                *status = Some(direct);
             }
         }
-        if let Some((block, element)) = assignment.iter().find_map(|(block, flags)| {
-            flags
-                .iter()
-                .position(Option::is_none)
-                .map(|element| (block, element))
-        }) {
-            return Err(Error::InvalidInput(format!(
-                "unassigned element {element} of curve {}",
-                block.curve()
-            )));
-        }
+        check_assigned(&assignment)?;
         Ok(Self {
             valuation_calendar,
             valuation_zone,
             curves,
             bindings,
             fx: Vec::new(),
+            volatilities: Vec::new(),
             lgm: BTreeMap::new(),
+            penalties: BTreeMap::new(),
+            gauss_newton: BTreeSet::new(),
             quote_kinds,
         })
     }
 
+    /// 会合日の階段のブロックに罰則を掛け、そのブロックを最小二乗で当てる。
+    ///
+    /// # Errors
+    ///
+    /// ブロックが市場の定義の曲線の `Steps` の部品でないか、区間が二つ未満で罰則の行が無いか、
+    /// 同じブロックに罰則を二度掛けた場合は `InvalidInput`。
+    pub fn with_step_penalties(
+        mut self,
+        penalties: impl IntoIterator<Item = StepPenalty>,
+    ) -> Result<Self> {
+        for penalty in penalties {
+            let steps = self
+                .curves
+                .iter()
+                .find(|curve| curve.id() == penalty.block.curve())
+                .and_then(|curve| {
+                    curve
+                        .parameter_blocks()
+                        .find(|(block, _)| *block == penalty.block)
+                        .filter(|_| curve.step_blocks().any(|block| block == penalty.block))
+                });
+            if steps.is_none_or(|(_, count)| count < 2) {
+                return Err(Error::InvalidInput(format!(
+                    "step penalty needs a steps part with two or more intervals: {} {}",
+                    penalty.block.curve(),
+                    penalty.block.part()
+                )));
+            }
+            if self
+                .penalties
+                .insert(penalty.block.clone(), penalty)
+                .is_some()
+            {
+                return Err(Error::InvalidInput("duplicate step penalty".into()));
+            }
+        }
+        Ok(self)
+    }
+
+    /// 最小二乗のグループの解の記録を Gauss–Newton（残差の 2 階の項を省く）で作るブロック。
+    ///
+    /// 既定は 2 階の項を含む厳密な記録で、このときだけ 1 次の quote 感応度が厳密になる。
+    /// 要素の多い当てはめで記録の費用を下げるために選べ、選んだことを解の記録と警告に残す。
+    ///
+    /// # Errors
+    ///
+    /// ブロックが市場の定義の曲線に無ければ `InvalidInput`。
+    pub fn with_gauss_newton(mut self, blocks: impl IntoIterator<Item = BlockId>) -> Result<Self> {
+        for block in blocks {
+            if !self
+                .curves
+                .iter()
+                .flat_map(CurveSpec::parameter_blocks)
+                .any(|(known, _)| known == block)
+            {
+                return Err(Error::InvalidInput(format!(
+                    "Gauss-Newton block {} {} is not in the market definition",
+                    block.curve(),
+                    block.part()
+                )));
+            }
+            self.gauss_newton.insert(block);
+        }
+        Ok(self)
+    }
+
     /// LGM の区間と quote ID を市場定義に追加する。
     ///
     /// # Errors
@@ -374,6 +630,14 @@ impl MarketDefinition {
     pub fn with_fx(mut self, definitions: impl IntoIterator<Item = FxDefinition>) -> Result<Self> {
         for definition in definitions {
             definition.validate()?;
+            // 外貨担保曲線の長期の XCCY は、その曲線の basis ブロックを当てる binding になる。
+            if let Some(implied) = definition
+                .direct
+                .as_ref()
+                .and_then(|direct| direct.implied.as_ref())
+            {
+                self.bindings.extend(implied.long_end_bindings()?);
+            }
             for (quote, kind) in definition.quote_kinds() {
                 if self
                     .quote_kinds
diff --git a/crates/rhodium-pricing/src/definition/fx.rs b/crates/rhodium-pricing/src/definition/fx.rs
index 90da2ebc692b101c2da2b59711bad199c9f1ade2..3558c7445f4571e9e93ead83f3b8282e80d036c7 100644
--- a/crates/rhodium-pricing/src/definition/fx.rs
+++ b/crates/rhodium-pricing/src/definition/fx.rs
@@ -2,14 +2,13 @@
 
 use rhodium_core::{
     Error, Result, Warning,
-    calendar::{
-        AdvanceRule, BusinessDayConvention, Calendar, CalendarSpec, EndOfMonthRule, adjust, advance,
-    },
+    calendar::{AdvanceRule, BusinessDayConvention, Calendar, EndOfMonthRule, adjust, advance},
     currency::Currency,
     date::Date,
     fx::CcyPair,
-    refdata::RefData,
-    tenor::{DateOffset, TenorUnit},
+    id::ConventionId,
+    refdata::{FxSpotLag, RefData},
+    tenor::{DateOffset, TenorLength, TenorUnit},
     time::{Instant, LocalEventTime, ZoneRules},
 };
 use rhodium_market::{
@@ -18,28 +17,10 @@ use rhodium_market::{
     quote::{QuoteId, QuoteKind},
     structure::BlockId,
 };
+use rhodium_products::standard::InstrumentTemplate;
 
-/// 現物 FX のスポットラグ。
-#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
-#[serde(rename_all = "snake_case")]
-pub enum FxSpotLag {
-    /// 翌営業日がスポット。
-    OneBusinessDay,
-    /// 二営業日後がスポット。
-    TwoBusinessDays,
-}
-
-/// 通貨の休日の和を使う現物 FX の受渡慣行。
-#[derive(Debug, Clone, serde::Deserialize)]
-#[serde(deny_unknown_fields)]
-pub struct FxSettlement {
-    /// スポットラグ。
-    pub spot_lag: FxSpotLag,
-    /// 受渡時刻。この瞬間を過ぎた当日の受渡は連鎖から外す。
-    pub settlement: LocalEventTime,
-    /// 1 pip がレートの 10 の負の何乗か。
-    pub pip_exponent: u8,
-}
+use super::{CalibrationTarget, ElementSelector, InstrumentBinding, QuoteTransform};
+use crate::resolve::CollateralKey;
 
 /// スポットから受渡日までの FX ポイント。
 #[derive(Debug, Clone, serde::Deserialize)]
@@ -68,23 +49,81 @@ pub struct FxQuotes {
     pub forwards: Vec<FxDatedQuote>,
 }
 
-/// FX ポイントから作る外貨担保割引曲線。
+/// FX ポイントと `MtM` XCCY から作る外貨担保割引曲線。
+///
+/// 短期は FX スワップのポイント（スポットと短期の quote、日付付きのポイント）、長期は `MtM` XCCY の
+/// ベーシスに校正する。どちらの行も同じ `basis` ブロックの節点を決め、一つのグループで解く。
 #[derive(Debug, Clone, serde::Deserialize)]
 #[serde(deny_unknown_fields)]
 pub struct FxImpliedCurve {
     /// 市場の担保の割引の組にある、今回作る側の曲線 ID。
     pub id: CurveId,
-    /// 最後の FX ポイントの受渡以後の補間規則。
+    /// 最後のピラー以後の補間規則。
     #[serde(default)]
     pub extrapolation: Extrapolation,
+    /// FX ポイントより先の、`MtM` XCCY のベーシスで決めるピラー。ピラーの昇順。
+    #[serde(default)]
+    pub long_end: Vec<FxImpliedLongEnd>,
+}
+
+/// 外貨担保曲線の長期の一つのピラーと、それを決める `MtM` XCCY のベーシスの quote。
+///
+/// XCCY の quote は外貨担保曲線の `basis` ブロックを当てる。ブロックの要素の番号は当日受渡の締めの
+/// 前後で（ON の区間が抜けて）ずれるので、一般の binding から要素を指さず、ここで宣言する。
+#[derive(Debug, Clone, serde::Deserialize)]
+#[serde(deny_unknown_fields)]
+pub struct FxImpliedLongEnd {
+    /// ピラーの UTC の瞬間。最後の FX ポイントの受渡より後。
+    pub pillar: Instant,
+    /// ベーシスの quote（スプレッド）。
+    pub quote: QuoteId,
+    /// `XccySwap` の慣行。払い脚にベーシスを付け、もう一方の脚の元本を FX でリセットする。
+    pub convention: ConventionId,
+    /// スポットからの満期。
+    pub tenor: TenorLength,
+    /// 校正商品の担保。両通貨の割引は価格設定の割引の表で解決する。
+    pub collateral: CollateralKey,
+}
+
+impl FxImpliedCurve {
+    /// FX ポイントと XCCY の行が当てる θ ブロック。
+    ///
+    /// # Errors
+    ///
+    /// 曲線 ID からブロック ID を作れない場合。
+    pub fn block(&self) -> Result<BlockId> {
+        BlockId::new(self.id.clone(), "basis")
+    }
+
+    /// 長期の XCCY の行を、外貨担保曲線のブロックを当てる校正の binding にする。
+    pub(crate) fn long_end_bindings(&self) -> Result<Vec<InstrumentBinding>> {
+        let block = self.block()?;
+        self.long_end
+            .iter()
+            .map(|row| {
+                InstrumentBinding::new(
+                    row.quote.clone(),
+                    QuoteKind::Spread,
+                    CalibrationTarget::ParRate(InstrumentTemplate::XccyBasis {
+                        convention: row.convention.clone(),
+                        tenor: row.tenor,
+                    }),
+                    QuoteTransform::Identity,
+                    row.collateral.clone(),
+                    ElementSelector::Block(block.clone()),
+                )
+            })
+            .collect()
+    }
 }
 
 /// 一つの直接ペアの起点と、その市場担保の校正。
+///
+/// スポットラグ、同日受渡の締めと pip の桁は、参照データのペアの
+/// [`FxConvention`](rhodium_core::refdata::FxConvention) から読む。
 #[derive(Debug, Clone, serde::Deserialize)]
 #[serde(deny_unknown_fields)]
 pub struct DirectFx {
-    /// 受渡慣行。
-    pub settlement: FxSettlement,
     /// ペアの市場の担保区分。
     pub market_collateral: CurveId,
     /// 市場の FX quote。
@@ -135,13 +174,28 @@ pub(crate) struct FxQuoteInterval {
 }
 
 impl FxDefinition {
-    /// 宣言された quote だけで決まる定義の整合性を、評価時点によらずに確かめる。
-    pub(crate) fn validate(&self) -> Result<()> {
+    /// 参照データのペアの慣行と、宣言された短期の quote が合うかを確かめる。
+    ///
+    /// quote は慣行の標準の向きで読むので、逆向きのペアの定義は拒否する。
+    /// T+1 のペアは ON だけ、T+2 のペアは TN（と ON）でスポットへ連鎖する。
+    ///
+    /// # Errors
+    ///
+    /// ペアの慣行が無いか、向きが慣行と違うか、短期の quote がスポットへ連鎖しなければ
+    /// `InvalidData`。
+    pub fn validate_references(&self, refdata: &RefData) -> Result<()> {
         let Some(direct) = &self.direct else {
             return Ok(());
         };
+        let convention = refdata.fx_convention(self.pair)?;
+        if convention.pair != self.pair {
+            return Err(Error::InvalidData(format!(
+                "market FX quotes for {} must use the convention direction {}",
+                self.pair, convention.pair
+            )));
+        }
         let quotes = &direct.quotes;
-        let chain = match direct.settlement.spot_lag {
+        let chain = match convention.spot_lag {
             FxSpotLag::OneBusinessDay => quotes.tn.is_none(),
             FxSpotLag::TwoBusinessDays => quotes.on.is_none() || quotes.tn.is_some(),
         };
@@ -151,13 +205,64 @@ impl FxDefinition {
                 self.pair
             )));
         }
+        // 長期の行は一般の binding になるとペアを失うので、ここで XCCY の慣行の 2 通貨を
+        // ペアに照らす。別のペアの慣行を書くと、外貨担保曲線を別の通貨の XCCY で当ててしまう。
+        let expected = [self.pair.base(), self.pair.quote()];
+        for row in direct.implied.iter().flat_map(|implied| &implied.long_end) {
+            let template = InstrumentTemplate::XccyBasis {
+                convention: row.convention.clone(),
+                tenor: row.tenor,
+            };
+            let legs = [
+                Some(template.currency(refdata)?),
+                template.other_currency(refdata)?,
+            ];
+            let matches = legs
+                .iter()
+                .all(|ccy| ccy.is_some_and(|ccy| expected.contains(&ccy)))
+                && legs[0] != legs[1];
+            if !matches {
+                return Err(Error::InvalidData(format!(
+                    "long-end XCCY convention {} of {} does not exchange {} and {}",
+                    row.convention,
+                    row.quote,
+                    self.pair.base(),
+                    self.pair.quote()
+                )));
+            }
+        }
+        Ok(())
+    }
+
+    /// 宣言された quote だけで決まる定義の整合性を、評価時点と参照データによらずに確かめる。
+    pub(crate) fn validate(&self) -> Result<()> {
+        let Some(direct) = &self.direct else {
+            return Ok(());
+        };
+        if let Some(implied) = &direct.implied
+            && !implied
+                .long_end
+                .windows(2)
+                .all(|pair| pair[0].pillar < pair[1].pillar)
+        {
+            return Err(Error::InvalidData(format!(
+                "long-end pillars of {} must be strictly increasing",
+                implied.id
+            )));
+        }
+        let quotes = &direct.quotes;
         // ON は当日受渡の締めを過ぎると連鎖から外れる。締めの前後で定義の成否が変わらないよう、
         // 外貨担保曲線の要否は使う区間ではなく宣言された quote で決める。
         // TN の区間は翌営業日に始まり、SN と日付付きの区間はスポットに始まるので、締め後も残る。
-        let after_cutoff =
+        // 長期の XCCY の行も締めによらずピラーを一つずつ決めるので、締め後に残る校正の対象に数える。
+        let surviving_points =
             quotes.tn.is_some() || quotes.sn.is_some() || !quotes.forwards.is_empty();
-        let points = quotes.on.is_some() || after_cutoff;
-        match (&direct.implied, points, after_cutoff) {
+        let points = quotes.on.is_some() || surviving_points;
+        let long_end = direct
+            .implied
+            .as_ref()
+            .is_some_and(|implied| !implied.long_end.is_empty());
+        match (&direct.implied, points, surviving_points || long_end) {
             // ポイントを市場の担保のフォワードで再現する自由度が無い（設計の fx_x0_on_tn_consistency）。
             (None, true, _) => Err(Error::InvalidData(format!(
                 "market FX points for {} require an FX-implied curve",
@@ -165,7 +270,7 @@ impl FxDefinition {
             ))),
             // T+1 のペアの ON だけでは、締め後に外貨担保曲線の校正区間が一つも残らない。
             (Some(_), _, false) => Err(Error::InvalidData(format!(
-                "FX-implied curve for {} needs a TN, SN or dated quote that remains after the same-day cutoff",
+                "FX-implied curve for {} needs a TN, SN, dated or long-end quote that remains after the same-day cutoff",
                 self.pair
             ))),
             _ => Ok(()),
@@ -175,13 +280,20 @@ impl FxDefinition {
     pub(crate) fn quote_kinds(&self) -> impl Iterator<Item = (&QuoteId, QuoteKind)> {
         self.direct.iter().flat_map(|direct| {
             let quotes = &direct.quotes;
-            std::iter::once((&quotes.spot, QuoteKind::FxSpot)).chain(
-                [&quotes.on, &quotes.tn, &quotes.sn]
-                    .into_iter()
-                    .filter_map(Option::as_ref)
-                    .chain(quotes.forwards.iter().map(|row| &row.quote))
-                    .map(|quote| (quote, QuoteKind::FxForwardPoints)),
-            )
+            let long_end = direct
+                .implied
+                .iter()
+                .flat_map(|implied| &implied.long_end)
+                .map(|row| (&row.quote, QuoteKind::Spread));
+            std::iter::once((&quotes.spot, QuoteKind::FxSpot))
+                .chain(
+                    [&quotes.on, &quotes.tn, &quotes.sn]
+                        .into_iter()
+                        .filter_map(Option::as_ref)
+                        .chain(quotes.forwards.iter().map(|row| &row.quote))
+                        .map(|quote| (quote, QuoteKind::FxForwardPoints)),
+                )
+                .chain(long_end)
         })
     }
 
@@ -204,11 +316,12 @@ impl FxDefinition {
                 warnings: Vec::new(),
             });
         };
+        self.validate_references(refdata)?;
         let mut warnings = Vec::new();
         let pair = self.pair;
-        let settlement = &direct.settlement;
-        let pip = 10_f64.powi(-i32::from(settlement.pip_exponent));
-        let calendar = FxCalendar::new(pair, &settlement.settlement, refdata)?;
+        let convention = refdata.fx_convention(pair)?;
+        let pip = convention.pip();
+        let calendar = FxCalendar::new(pair, &convention.settlement, refdata)?;
         let local_date = calendar.zone.to_local(as_of)?.date();
         let (trade, adjust_warnings) = adjust(
             local_date,
@@ -219,7 +332,7 @@ impl FxDefinition {
         let trade_at = calendar.instant(trade)?;
         let tomorrow = calendar.next(trade, &mut warnings)?;
         let tomorrow_at = calendar.instant(tomorrow)?;
-        let spot_date = match settlement.spot_lag {
+        let spot_date = match convention.spot_lag {
             FxSpotLag::OneBusinessDay => tomorrow,
             FxSpotLag::TwoBusinessDays => calendar.next(tomorrow, &mut warnings)?,
         };
@@ -297,9 +410,15 @@ impl FxDefinition {
         pip: f64,
     ) -> Result<Option<FxFit>> {
         let pair = self.pair;
-        match (&direct.implied, intervals.is_empty()) {
-            (None, true) => Ok(None),
-            (Some(implied), false) => {
+        // 締め後に FX ポイントの区間が残らなくても、長期の XCCY がピラーを決める。
+        let targets = !intervals.is_empty()
+            || direct
+                .implied
+                .as_ref()
+                .is_some_and(|implied| !implied.long_end.is_empty());
+        match (&direct.implied, targets) {
+            (None, false) => Ok(None),
+            (Some(implied), true) => {
                 let collateral = self
                     .collaterals
                     .iter()
@@ -317,7 +436,17 @@ impl FxDefinition {
                         )));
                     }
                 };
-                let target = BlockId::new(implied.id.clone(), "basis")?;
+                let target = implied.block()?;
+                // 短期は FX ポイントの受渡、長期は XCCY のピラー。ピラーの順に一つずつの行が決める。
+                let last = intervals.iter().map(|interval| interval.end).max();
+                if let (Some(last), Some(first)) = (last, implied.long_end.first())
+                    && first.pillar <= last
+                {
+                    return Err(Error::InvalidData(format!(
+                        "long-end pillar of {} must follow the last FX point settlement",
+                        implied.id
+                    )));
+                }
                 let curve = CurveSpec::from_expr(
                     implied.id.clone(),
                     ccy,
@@ -326,7 +455,11 @@ impl FxDefinition {
                         part: target.part().into(),
                         pair,
                         base: opposite.clone(),
-                        pillars: intervals.iter().map(|interval| interval.end).collect(),
+                        pillars: intervals
+                            .iter()
+                            .map(|interval| interval.end)
+                            .chain(implied.long_end.iter().map(|row| row.pillar))
+                            .collect(),
                         extrapolation: implied.extrapolation,
                     },
                 )?;
@@ -339,10 +472,10 @@ impl FxDefinition {
                 }))
             }
             // `validate` が宣言された quote で拒否するので、次の二つは定義の検査を経ない場合の防御。
-            (None, false) => Err(Error::InvalidData(format!(
+            (None, true) => Err(Error::InvalidData(format!(
                 "market FX points for {pair} require an FX-implied curve"
             ))),
-            (Some(_), true) => Err(Error::InvalidData(format!(
+            (Some(_), false) => Err(Error::InvalidData(format!(
                 "FX-implied curve for {pair} has no eligible calibration interval"
             ))),
         }
@@ -357,15 +490,8 @@ struct FxCalendar<'a> {
 
 impl<'a> FxCalendar<'a> {
     fn new(pair: CcyPair, settlement: &'a LocalEventTime, refdata: &'a RefData) -> Result<Self> {
-        let base = refdata.currency(pair.base())?.settlement_calendar.clone();
-        let quote = refdata.currency(pair.quote())?.settlement_calendar.clone();
-        let spec = if base == quote {
-            CalendarSpec::Base(base)
-        } else {
-            CalendarSpec::Union(vec![CalendarSpec::Base(base), CalendarSpec::Base(quote)])
-        };
         Ok(Self {
-            joint: refdata.calendar(&spec)?,
+            joint: refdata.calendar(&refdata.fx_calendar_spec(pair)?)?,
             settlement,
             zone: refdata.zone(&settlement.zone)?,
         })
diff --git a/crates/rhodium-pricing/src/definition/vol.rs b/crates/rhodium-pricing/src/definition/vol.rs
new file mode 100644
index 0000000000000000000000000000000000000000..35980236454af0a45377ea617e217edcad8682e1
--- /dev/null
+++ b/crates/rhodium-pricing/src/definition/vol.rs
@@ -0,0 +1,419 @@
+//! 金利ボラ面を市場 quote と慣行に結ぶ指定。
+
+use std::collections::{BTreeMap, BTreeSet};
+
+use rhodium_core::{
+    Error, Result,
+    id::{ConventionId, IndexId},
+    quantity::VolKind,
+    refdata::{Convention, RefData, SwapCouponConvention},
+    tenor::TenorLength,
+    time::LocalEventTime,
+};
+use rhodium_market::{
+    quote::{QuoteId, QuoteKind},
+    vol::{CubeInterpolation, Sabr, SabrFit, VolId},
+};
+
+use crate::resolve::CollateralKey;
+
+/// ボラの市場での測り方。シフト付きの場合のシフトは面に指定する。
+#[derive(Debug, Clone, Copy, serde::Deserialize)]
+#[serde(rename_all = "snake_case")]
+pub enum RateVolKind {
+    /// Bachelier の絶対ボラ。
+    Normal,
+    /// Black の相対ボラ。
+    Lognormal,
+    /// フォワードとストライクをシフトした Black の相対ボラ。
+    ShiftedLognormal,
+}
+
+impl RateVolKind {
+    fn kind(self, shift: Option<f64>) -> Result<VolKind> {
+        match (self, shift) {
+            (Self::Normal, None) => Ok(VolKind::Normal),
+            (Self::Lognormal, None) => Ok(VolKind::Lognormal),
+            (Self::ShiftedLognormal, Some(shift)) if shift.is_finite() => {
+                Ok(VolKind::ShiftedLognormal { shift })
+            }
+            _ => Err(Error::InvalidData(
+                "specify a finite shift only for shifted_lognormal".into(),
+            )),
+        }
+    }
+}
+
+/// ストライクの座標。数値は金利の十進表記。
+#[derive(Debug, Clone, Copy, serde::Deserialize)]
+#[serde(rename_all = "snake_case")]
+pub enum RateStrike {
+    /// フォワードの位置。
+    Atm,
+    /// フォワードからの金利差。
+    AtmOffset(f64),
+    /// 金利で直接指定する行使価格。
+    Absolute(f64),
+}
+
+impl RateStrike {
+    pub(crate) fn value<R: scioxide_ad::Real>(self, forward: R) -> R {
+        match self {
+            Self::Atm => forward,
+            Self::AtmOffset(offset) => forward + offset,
+            Self::Absolute(strike) => R::constant(strike),
+        }
+    }
+
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#201 で各 variant を網羅する"
+    )]
+    fn check(self, cap: bool) -> Result<()> {
+        match self {
+            Self::Atm => Ok(()),
+            Self::AtmOffset(value) if !cap && value.is_finite() && value != 0.0 => Ok(()),
+            Self::Absolute(value) if value.is_finite() => Ok(()),
+            _ => Err(Error::InvalidData(
+                "invalid or unsupported volatility strike".into(),
+            )),
+        }
+    }
+}
+
+/// スワップションの一つの市場 quote。差分の基準は同じノードの ATM。
+#[derive(Debug, Clone, serde::Deserialize)]
+#[serde(deny_unknown_fields)]
+pub struct SwaptionVolQuote {
+    /// 入力 quote。
+    pub quote: QuoteId,
+    /// 取引日からオプション満期まで。
+    pub expiry: TenorLength,
+    /// 満期で取引する原資産スワップのテナー。
+    pub tenor: TenorLength,
+    /// 行使価格の座標。
+    pub strike: RateStrike,
+    /// ボラの quote 自体が差分のとき参照する ATM の入力 quote。
+    pub atm_quote: Option<QuoteId>,
+    /// 最小二乗行の当てはめ幅。同じボラの単位で指定する。
+    pub tolerance: Option<f64>,
+}
+
+/// キャップのフラットボラの quote。
+#[derive(Debug, Clone, serde::Deserialize)]
+#[serde(deny_unknown_fields)]
+pub struct CapVolQuote {
+    /// 入力 quote。
+    pub quote: QuoteId,
+    /// スポットのキャップの満期。
+    pub maturity: TenorLength,
+    /// ATM または絶対ストライク。
+    pub strike: RateStrike,
+}
+
+/// 後決めキャップの市場 quote の分散規約。
+#[derive(Debug, Clone, Copy, serde::Deserialize)]
+#[serde(rename_all = "snake_case")]
+pub enum RfrVariance {
+    /// 観測期間でボラが線形に減る。
+    LyashenkoMercurio,
+}
+
+/// ボラ面と校正用クーポンの共通の定義。
+#[derive(Debug, Clone, serde::Deserialize)]
+pub struct VolatilityMarket {
+    /// 面の ID。
+    pub id: VolId,
+    /// 変動金利の指標。
+    pub index: IndexId,
+    /// 参照する IRS または OIS の慣行。
+    pub convention: ConventionId,
+    /// 割引を決める担保。
+    pub collateral: CollateralKey,
+    /// 市場 quote のボラの種類。
+    pub vol_kind: RateVolKind,
+    /// シフト付き対数正規 quote のシフト。
+    pub shift: Option<f64>,
+}
+
+impl VolatilityMarket {
+    pub(crate) fn kind(&self) -> Result<VolKind> {
+        self.vol_kind.kind(self.shift)
+    }
+
+    fn quote_kind(&self, difference: bool) -> QuoteKind {
+        match (self.vol_kind, difference) {
+            (RateVolKind::Normal, false) => QuoteKind::NormalVol,
+            (RateVolKind::Normal, true) => QuoteKind::NormalVolDifference,
+            (RateVolKind::Lognormal, false) => QuoteKind::LognormalVol,
+            (RateVolKind::ShiftedLognormal, false) => QuoteKind::ShiftedLognormalVol,
+            (_, true) => QuoteKind::LognormalVolDifference,
+        }
+    }
+
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#197 で各 variant を網羅する"
+    )]
+    pub(crate) fn template(
+        &self,
+        tenor: TenorLength,
+        refdata: &RefData,
+    ) -> Result<rhodium_products::standard::InstrumentTemplate> {
+        use rhodium_products::standard::InstrumentTemplate;
+        let invalid = || {
+            Error::InvalidData(format!(
+                "{} does not identify an IRS/OIS on {}",
+                self.convention, self.index
+            ))
+        };
+        let convention = self.convention.clone();
+        let template = match refdata.convention(&convention)? {
+            Convention::Ois(ois) if ois.float.index == self.index => {
+                InstrumentTemplate::Ois { convention, tenor }
+            }
+            Convention::Irs(swap) => {
+                let has_fixed = swap
+                    .legs
+                    .iter()
+                    .any(|leg| matches!(leg.coupon, SwapCouponConvention::Fixed));
+                let has_index = swap.legs.iter().any(|leg| match &leg.coupon {
+                    SwapCouponConvention::Term { index } => *index == self.index,
+                    SwapCouponConvention::Compounded(rule) => rule.index == self.index,
+                    SwapCouponConvention::Fixed => false,
+                });
+                if !has_fixed || !has_index {
+                    return Err(invalid());
+                }
+                InstrumentTemplate::Swap { convention, tenor }
+            }
+            _ => return Err(invalid()),
+        };
+        Ok(template)
+    }
+}
+
+/// Hagan の展開と最小二乗の指定。
+#[derive(Debug, Clone, Copy, serde::Deserialize)]
+#[serde(deny_unknown_fields)]
+pub struct SabrCalibration {
+    /// 固定する弾性 β。
+    pub beta: f64,
+    /// 校正するパラメータの範囲。
+    pub fit: SabrFit,
+    /// 正規 SABR のシフト。対数正規のシフトは面で指定する。
+    #[serde(default)]
+    pub shift: f64,
+    /// キューブの補間。
+    #[serde(default)]
+    pub interpolation: CubeInterpolation,
+}
+
+/// 約定日から格子と校正商品を解決する金利ボラ面の定義。
+#[derive(Debug, Clone, serde::Deserialize)]
+#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
+pub enum VolatilityDefinition {
+    /// スワップションの SABR キューブ。
+    SwaptionCube {
+        /// 面と原資産の指定。
+        #[serde(flatten)]
+        surface: VolatilityMarket,
+        /// 各満期の現地カット。
+        expiry_time: LocalEventTime,
+        /// Hagan の展開と当てはめの指定。
+        sabr: SabrCalibration,
+        /// 満期、テナーとストライクで指定した quote。
+        quotes: Vec<SwaptionVolQuote>,
+    },
+    /// キャップのフラットボラから作るキャップレット面。
+    Caplet {
+        /// 面と原資産の指定。
+        #[serde(flatten)]
+        surface: VolatilityMarket,
+        /// 後決め RFR の市場 quote が使う実効分散の定義。
+        rfr_variance: Option<RfrVariance>,
+        /// キャップの quote。
+        quotes: Vec<CapVolQuote>,
+    },
+}
+
+impl VolatilityDefinition {
+    /// 面の ID。
+    #[must_use]
+    pub fn id(&self) -> &VolId {
+        &self.surface().id
+    }
+
+    pub(crate) const fn surface(&self) -> &VolatilityMarket {
+        match self {
+            Self::SwaptionCube { surface, .. } | Self::Caplet { surface, .. } => surface,
+        }
+    }
+
+    pub(crate) fn sabr(surface: &VolatilityMarket, spec: SabrCalibration) -> Result<Sabr> {
+        let kind = surface.kind()?;
+        let shift = match kind {
+            VolKind::Normal => spec.shift,
+            VolKind::Lognormal | VolKind::ShiftedLognormal { .. } if spec.shift != 0.0 => {
+                return Err(Error::InvalidData(
+                    "use surface shift for shifted lognormal SABR".into(),
+                ));
+            }
+            VolKind::Lognormal => 0.0,
+            VolKind::ShiftedLognormal { shift } => shift,
+        };
+        Sabr::new(kind, spec.beta, shift)
+    }
+
+    /// 参照データと慣行・指標の組を照合する。
+    ///
+    /// # Errors
+    ///
+    /// 未登録の参照先または互換性のない指標・分散規約・ボラ種別。
+    pub fn validate_references(&self, refdata: &RefData) -> Result<()> {
+        use rhodium_core::refdata::IndexKind;
+        let surface = self.surface();
+        surface.collateral.validate_references(refdata)?;
+        let index = refdata.index(&surface.index)?;
+        let tenor =
+            match self {
+                Self::SwaptionCube {
+                    expiry_time,
+                    quotes,
+                    ..
+                } => {
+                    refdata.zone(&expiry_time.zone)?;
+                    quotes.first().map(|q| q.tenor)
+                }
+                Self::Caplet {
+                    rfr_variance,
+                    quotes,
+                    ..
+                } => {
+                    match (&index.kind, rfr_variance, surface.kind()?) {
+                        (
+                            IndexKind::OvernightRate { .. },
+                            Some(RfrVariance::LyashenkoMercurio),
+                            VolKind::Normal,
+                        )
+                        | (IndexKind::TermRate { .. }, None, _) => {}
+                        _ => return Err(Error::InvalidData(
+                            "cap index requires compatible volatility and RFR variance conventions"
+                                .into(),
+                        )),
+                    }
+                    quotes.first().map(|q| q.maturity)
+                }
+            }
+            .ok_or_else(|| Error::InvalidData("empty volatility quote grid".into()))?;
+        surface.template(tenor, refdata)?.validate(refdata)
+    }
+
+    pub(crate) fn quote_kinds(&self) -> Result<BTreeMap<QuoteId, QuoteKind>> {
+        let mut kinds = BTreeMap::new();
+        let surface = self.surface();
+        surface.kind()?;
+        let mut add = |quote: &QuoteId, difference| {
+            if kinds
+                .insert(quote.clone(), surface.quote_kind(difference))
+                .is_some()
+            {
+                return Err(Error::InvalidData(format!(
+                    "duplicate volatility quote {quote}"
+                )));
+            }
+            Ok(())
+        };
+        match self {
+            Self::SwaptionCube { sabr, quotes, .. } => {
+                Self::sabr(surface, *sabr)?;
+                // 固定の ρ・ν も β と同じく読み込みの境界で検証し、校正まで持ち越さない。
+                sabr.fit.check()?;
+                for row in quotes {
+                    row.strike.check(false)?;
+                    if row.expiry.amount() == 0 || row.tenor.amount() == 0 {
+                        return Err(Error::InvalidData(format!(
+                            "zero swaption tenor for {}",
+                            row.quote
+                        )));
+                    }
+                    let weighted = match sabr.fit {
+                        SabrFit::AtmOnly { .. } => false,
+                        SabrFit::LeastSquares => true,
+                        SabrFit::ExactAtm => !matches!(row.strike, RateStrike::Atm),
+                    };
+                    if weighted != row.tolerance.is_some()
+                        || row.tolerance.is_some_and(|x| !(x.is_finite() && x > 0.0))
+                    {
+                        return Err(Error::InvalidData(format!(
+                            "invalid fit tolerance for {}",
+                            row.quote
+                        )));
+                    }
+                    add(&row.quote, row.atm_quote.is_some())?;
+                }
+            }
+            Self::Caplet { quotes, .. } => {
+                for row in quotes {
+                    row.strike.check(true)?;
+                    if row.maturity.amount() == 0 {
+                        return Err(Error::InvalidData(format!(
+                            "zero cap maturity for {}",
+                            row.quote
+                        )));
+                    }
+                    add(&row.quote, false)?;
+                }
+            }
+        }
+        if kinds.is_empty() {
+            return Err(Error::InvalidData("empty volatility quote grid".into()));
+        }
+        Ok(kinds)
+    }
+}
+
+impl super::MarketDefinition {
+    /// 既存の市場定義に金利ボラ面とその quote binding を追加する。
+    ///
+    /// # Errors
+    ///
+    /// 面の ID が重複するか、格子や quote の種類が不正な場合。
+    pub fn with_volatilities(
+        mut self,
+        definitions: impl IntoIterator<Item = VolatilityDefinition>,
+    ) -> Result<Self> {
+        let mut ids: BTreeSet<_> = self
+            .volatilities
+            .iter()
+            .map(|def| def.id().clone())
+            .collect();
+        for def in definitions {
+            if !ids.insert(def.id().clone())
+                || self
+                    .curves
+                    .iter()
+                    .any(|curve| curve.id().as_str() == def.id().as_str())
+            {
+                return Err(Error::InvalidData(format!(
+                    "duplicate market object {}",
+                    def.id().as_str()
+                )));
+            }
+            for (quote, kind) in def.quote_kinds()? {
+                if self
+                    .quote_kinds
+                    .insert(quote.clone(), kind)
+                    .is_some_and(|previous| previous != kind)
+                {
+                    return Err(Error::InvalidData(format!(
+                        "inconsistent kinds for quote {quote}"
+                    )));
+                }
+            }
+            self.volatilities.push(def);
+        }
+        self.volatilities.sort_by(|a, b| a.id().cmp(b.id()));
+        Ok(self)
+    }
+}
diff --git a/crates/rhodium-pricing/src/lib.rs b/crates/rhodium-pricing/src/lib.rs
index 5b4cda5574b5aa6c73fc915accf51df22e2c99bf..55db1bfba1024f0b51ffa53bfaf70d843bc9b461 100644
--- a/crates/rhodium-pricing/src/lib.rs
+++ b/crates/rhodium-pricing/src/lib.rs
@@ -6,6 +6,7 @@
 #[cfg(test)]
 extern crate self as rhodium_pricing;
 
+pub mod analytic;
 pub mod calibration;
 pub mod config;
 pub mod context;
