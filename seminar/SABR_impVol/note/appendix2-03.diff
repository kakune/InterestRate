diff --git a/crates/rhodium-pricing/src/linear.rs b/crates/rhodium-pricing/src/linear.rs
index 10304f2afd8831b27b5a90166dfd3c63b3e64917..96b47124728c10b375bd162c1f473a93bfb07b67 100644
--- a/crates/rhodium-pricing/src/linear.rs
+++ b/crates/rhodium-pricing/src/linear.rs
@@ -5,7 +5,7 @@ mod fx;
 mod rfr;
 
 use std::{
-    collections::{BTreeMap, BTreeSet, HashMap, btree_map::Entry, hash_map},
+    collections::{BTreeMap, BTreeSet, HashMap, HashSet, btree_map::Entry, hash_map},
     sync::Arc,
 };
 
@@ -29,7 +29,10 @@ use rhodium_products::ir::{
     Product,
     analysis::{Cashflow, ExprRef, Expression, FxFixing, cashflows},
 };
-use rhodium_products::standard::RateFutureMargin;
+use rhodium_products::{
+    standard::{PayReceive, RateFutureMargin},
+    trade::Fee,
+};
 use scioxide_ad::{Mask, Real};
 
 use crate::{
@@ -369,6 +372,7 @@ pub(super) fn prepare(
 
 pub(super) fn prepare_rate_future(
     margin: &RateFutureMargin<'_>,
+    fees: &[Fee],
     warnings: Vec<Warning>,
     collateral: &CollateralKey,
     structure: &Arc<MarketStructure>,
@@ -377,6 +381,23 @@ pub(super) fn prepare_rate_future(
     use rhodium_core::fixing::{FixingResolution, ObservationPeriod};
 
     let mut builder = PlanBuilder::new(collateral, structure, ctx, warnings);
+    // 取引の手数料は値洗いと独立の受渡で、IR の評価と同じく決済の境界より後のものだけを数える。
+    for fee in fees
+        .iter()
+        .filter(|fee| ctx.config.settlement().includes(ctx.as_of, fee.at))
+    {
+        let signed = match fee.side {
+            PayReceive::Pay => -fee.amount.amount(),
+            PayReceive::Receive => fee.amount.amount(),
+        };
+        let discount = builder.discount(fee.amount.ccy())?;
+        let point = builder.payment_point(discount, fee.at)?;
+        builder.discounts[discount].group.payments.push(Payment {
+            point,
+            amount: vec![Instruction::Constant(signed)],
+            compounds: Vec::new(),
+        });
+    }
     let first = margin
         .clearings
         .partition_point(|row| !ctx.config.settlement().includes(ctx.as_of, row.settlement));
@@ -533,7 +554,11 @@ impl<'a> PlanBuilder<'a> {
         pay: Instant,
     ) -> Result<usize> {
         let value = fx::value_time(fixing, self.ctx, &mut self.warnings)?;
-        if fixing.value_at.is_none() && value != pay {
+        // 予想はスポットの受渡で読むので、その受渡と実際に受け渡す瞬間の遅れを照合する。
+        // 通常の fixing は支払いの瞬間、元本リセットは新しい元本を交換する瞬間と比べる。リセットの
+        // 元本に掛かる期末のクーポンの遅れは、金利の予想との積の規約 MtmFxForwardProduct が扱う。
+        let exchanged = fixing.value_at.unwrap_or(pay);
+        if value != exchanged {
             let terms = self.ctx.refdata.currency(currency)?;
             let date = self
                 .ctx
@@ -541,7 +566,8 @@ impl<'a> PlanBuilder<'a> {
                 .zone(&terms.settlement_time.zone)?
                 .to_local(value)?
                 .date();
-            let warning = rfr::payment_delay(self.ctx, currency, date, pay, &mut self.warnings)?;
+            let warning =
+                rfr::payment_delay(self.ctx, currency, date, exchanged, &mut self.warnings)?;
             self.warnings.extend(warning);
         }
         match self
@@ -615,6 +641,10 @@ impl<'a> PlanBuilder<'a> {
         clippy::too_many_lines,
         reason = "式の種類ごとの命令を一箇所で組み立てる"
     )]
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#202 で各 variant を網羅する"
+    )]
     fn amount(
         &mut self,
         product: &Product,
@@ -666,7 +696,7 @@ impl<'a> PlanBuilder<'a> {
         let mut instructions = Vec::with_capacity(reached.len());
         let mut dynamic = Vec::with_capacity(reached.len());
         let mut compounds = CompoundPlans::default();
-        let is_fra = reached.values().any(|expr| expr.fra_settlement().is_some());
+        let term_uses = term_uses(root);
         for expr in reached.into_values() {
             if let Some(compound) = excess_growth(expr) {
                 // 複利から 1 を引く式は、成長率を丸めてから引くと短い期間の超過の桁を落とすので、
@@ -728,11 +758,23 @@ impl<'a> PlanBuilder<'a> {
                     let slot = self.forecast(fixing.index)?;
                     let forecast = &mut self.forecasts[slot];
                     let discount = &self.discounts[discount].curve;
-                    let point = if is_fra {
-                        forecast.fra_point(fixing, discount, at, self.ctx, &mut self.warnings)?
-                    } else {
-                        forecast.term_point(fixing, discount, at, self.ctx, &mut self.warnings)?
-                    };
+                    // 支払いの遅れは使い方ごとに照合する。FRA の決済の中では期首払い、それ以外は
+                    // 期末のクーポンとして扱う。両方で使う観測は両方を照合し、予想は同じ点を共有する。
+                    let uses = term_uses.get(&expr.index()).copied().unwrap_or_default();
+                    let fra = uses
+                        .fra
+                        .then(|| {
+                            forecast.fra_point(fixing, discount, at, self.ctx, &mut self.warnings)
+                        })
+                        .transpose()?;
+                    let plain = (uses.plain || fra.is_none())
+                        .then(|| {
+                            forecast.term_point(fixing, discount, at, self.ctx, &mut self.warnings)
+                        })
+                        .transpose()?;
+                    let point = plain
+                        .or(fra)
+                        .ok_or_else(|| Error::Internal("unclassified term fixing".into()))?;
                     Instruction::Forward { slot, point }
                 }
                 (Expression::FxFixing(fixing), []) => {
@@ -854,3 +896,38 @@ impl<'a> PlanBuilder<'a> {
         })
     }
 }
+
+/// 支払いの式の中で、一つのターム金利の観測を使う式の種類。
+#[derive(Clone, Copy, Debug, Default)]
+struct TermUse {
+    /// 認識した FRA の決済額の外（期末のクーポンなど）で使う。
+    plain: bool,
+    /// 認識した FRA の決済額（期首払い）の中で使う。
+    fra: bool,
+}
+
+/// 根からの経路ごとに、各ターム金利の観測が認識した FRA の決済額の中にあるかを調べる。
+///
+/// 一つの支払いに FRA の決済と通常のタームクーポンが並んでも、FRA の期首払いの扱いは
+/// その FRA の観測にだけ当てる。共有した観測は、両方の経路で使えば両方の印を持つ。
+fn term_uses(root: ExprRef<'_>) -> HashMap<usize, TermUse> {
+    let mut uses: HashMap<usize, TermUse> = HashMap::new();
+    let mut seen = HashSet::new();
+    let mut work = vec![(root, false)];
+    while let Some((expr, inside)) = work.pop() {
+        if !seen.insert((expr.index(), inside)) {
+            continue;
+        }
+        if matches!(expr.expression(), Expression::TermFixing(_)) {
+            let entry = uses.entry(expr.index()).or_default();
+            if inside {
+                entry.fra = true;
+            } else {
+                entry.plain = true;
+            }
+        }
+        let inside = inside || expr.fra_settlement().is_some();
+        work.extend(expr.children().map(|child| (child, inside)));
+    }
+    uses
+}
diff --git a/crates/rhodium-pricing/src/linear/future.rs b/crates/rhodium-pricing/src/linear/future.rs
index dcbbec08b1391e96e0ca5b8d542d1ab1364bbaa9..8735a0f1218d09820014dd61da51401c40985966 100644
--- a/crates/rhodium-pricing/src/linear/future.rs
+++ b/crates/rhodium-pricing/src/linear/future.rs
@@ -58,6 +58,19 @@ struct RateLeg {
     start: usize,
     end: usize,
     weight: f64,
+    span: Span,
+}
+
+/// 予想する 1 本の率が覆う観測の期間。
+#[derive(Clone, Copy, Debug)]
+enum Span {
+    /// 観測の期間の全体の成長 G(開始, 終了) から率を作る。
+    Whole,
+    /// 評価時点 A が期間をまたぐ。[A, 終了) の成長を残りの日数係数 `τ_rem` で年率化し、
+    /// その率を元の期間の日数係数 τ に当てる（参照データの章の fixing の規則）。
+    /// `scale` は `τ/τ_rem` で、元の期間の増分は `τ·(G − 1)/τ_rem = scale·(G − 1)`。
+    /// 成長は評価時点で決まるので、fixing の瞬間を評価時点に置く。
+    FromAsOf { scale: f64 },
 }
 
 #[derive(Debug, Clone, Copy)]
@@ -72,6 +85,7 @@ struct ForecastLeg {
     start: Instant,
     end: Instant,
     weight: f64,
+    span: Span,
 }
 
 impl FuturePrice {
@@ -102,8 +116,13 @@ impl FuturePrice {
                 let (period, notices) = ctx.refdata.term_rate_dates(underlying, fixing)?;
                 warnings.extend(notices);
                 let tau = fraction(basis, period, ctx, warnings)?;
-                let fixed = read_rate(underlying, fixing, tau.recip(), &mut legs, ctx, warnings)?
-                    .unwrap_or(0.0);
+                let read = RateRead {
+                    weight: tau.recip(),
+                    accrual: tau,
+                    basis,
+                };
+                let fixed =
+                    read_rate(underlying, fixing, read, &mut legs, ctx, warnings)?.unwrap_or(0.0);
                 RatePlan::Sum { fixed }
             }
             RateFutureRate::Compounded(period) | RateFutureRate::Averaged(period) => {
@@ -256,6 +275,7 @@ impl ProjectedPrice {
                     start: position(leg.start)?,
                     end: position(leg.end)?,
                     weight: leg.weight,
+                    span: leg.span,
                 })
             })
             .collect::<Result<_>>()?;
@@ -332,12 +352,18 @@ impl ProjectedPrice {
         };
         let rate =
             match self.rate {
+                // 期間をまたぐ脚は fixing が評価時点（ζ = 0、傾き 0）なので成長は決定論的で、
+                // 期待値は G そのもの。増分 G − 1 を scale 倍して元の期間の増分にする。
                 RatePlan::Sum { fixed } => self.legs.iter().copied().enumerate().fold(
                     R::constant(fixed),
                     |rate, (i, leg)| {
                         let (log, loading) = log_expected(leg, i);
+                        let scale = match leg.span {
+                            Span::Whole => 1.0,
+                            Span::FromAsOf { scale } => scale,
+                        };
                         rate + (log + loading * loading * zeta[leg.fixing] * 0.5).exp_m1()
-                            * leg.weight
+                            * (leg.weight * scale)
                     },
                 ),
                 RatePlan::Compound { log_fixed, tau } => {
@@ -345,8 +371,17 @@ impl ProjectedPrice {
                     let mut log = R::constant(log_fixed);
                     for (i, leg) in self.legs.iter().copied().enumerate() {
                         let (component, loading) = log_expected(leg, i);
-                        log += component;
-                        weights[leg.fixing] += loading;
+                        match leg.span {
+                            Span::Whole => {
+                                log += component;
+                                weights[leg.fixing] += loading;
+                            }
+                            // 元の期間の複利因子 1 + τ·r = 1 + scale·(G − 1)。対数正規ではないが、
+                            // 決定論的なので共分散の項に加えない。
+                            Span::FromAsOf { scale } => {
+                                log += (component.exp_m1() * scale).ln_1p();
+                            }
+                        }
                     }
                     let mut suffix = R::zero();
                     for i in (1..n).rev() {
@@ -401,12 +436,16 @@ fn overnight_plan(
                     .into(),
             ));
         }
-        let weight = if compound {
-            1.0
-        } else {
-            covered / (actual * tau)
+        let read = RateRead {
+            weight: if compound {
+                1.0
+            } else {
+                covered / (actual * tau)
+            },
+            accrual: actual,
+            basis,
         };
-        if let Some(value) = read_rate(underlying, start, weight, legs, ctx, warnings)? {
+        if let Some(value) = read_rate(underlying, start, read, legs, ctx, warnings)? {
             if compound {
                 let simple = actual * value;
                 if simple <= -1.0 {
@@ -439,6 +478,10 @@ fn integration_factor<R: Real>(value: R) -> R {
     R::select(small, series, -(-value).exp_m1() / safe)
 }
 
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#197 で各 variant を網羅する"
+)]
 fn fraction(
     basis: &DayCount,
     interval: DateInterval,
@@ -466,10 +509,21 @@ fn fraction(
     Ok(fraction.value())
 }
 
+/// 1 本の率の読み方。
+#[derive(Clone, Copy)]
+struct RateRead<'a> {
+    /// 率の増分に掛ける重み。
+    weight: f64,
+    /// 率を当てる観測の期間の全体の日数係数 τ。
+    accrual: f64,
+    /// 日数係数の基準。
+    basis: &'a DayCount,
+}
+
 fn read_rate(
     index: &IndexId,
     date: Date,
-    weight: f64,
+    read: RateRead<'_>,
     legs: &mut Vec<ForecastLeg>,
     ctx: &Context<'_>,
     warnings: &mut BTreeSet<Warning>,
@@ -494,16 +548,41 @@ fn read_rate(
     warnings.extend(warning);
     match resolution {
         FixingResolution::Published(rate) | FixingResolution::Estimated(rate) => Ok(Some(rate)),
-        FixingResolution::Forecast | FixingResolution::ForecastFrom(_) => {
-            let start = match resolution {
-                FixingResolution::ForecastFrom(at) => at,
-                _ => period.start(),
-            };
+        FixingResolution::Forecast => {
             legs.push(ForecastLeg {
                 fixing: fixing.max(ctx.as_of),
-                start,
+                start: period.start(),
                 end: period.end(),
-                weight,
+                weight: read.weight,
+                span: Span::Whole,
+            });
+            Ok(None)
+        }
+        // 評価時点が期間をまたぐ。[A, 終了) の forward の率を期間全体の率とみなす（RFR の
+        // `insert_forward` と同じ規則）。残りの日数係数は観測時刻の間を時間で按分して求める。
+        FixingResolution::ForecastFrom(at) => {
+            let zone = ctx.refdata.zone(&event.zone)?;
+            let end = zone.to_local(period.end())?.date();
+            let remaining = super::rfr::remaining_fraction(at, end, zone, event, |from, to| {
+                if from >= to {
+                    Ok(0.0)
+                } else {
+                    fraction(read.basis, DateInterval::new(from, to)?, ctx, warnings)
+                }
+            })?;
+            if remaining <= 0.0 {
+                return Err(Error::Numerical(
+                    "futures forecast requires a positive remaining fraction".into(),
+                ));
+            }
+            legs.push(ForecastLeg {
+                fixing: at,
+                start: at,
+                end: period.end(),
+                weight: read.weight,
+                span: Span::FromAsOf {
+                    scale: read.accrual / remaining,
+                },
             });
             Ok(None)
         }
diff --git a/crates/rhodium-pricing/src/linear/fx.rs b/crates/rhodium-pricing/src/linear/fx.rs
index c9628cd2aa6b2745e2548fe86f8789750eb60625..0f6e34fd7a4fa23141e38c85f1e8331b5226c578 100644
--- a/crates/rhodium-pricing/src/linear/fx.rs
+++ b/crates/rhodium-pricing/src/linear/fx.rs
@@ -6,7 +6,6 @@ use rhodium_core::{
     Error, Result, Warning,
     calendar::{AdvanceRule, BusinessDayConvention, EndOfMonthRule, advance},
     fixing::ObservationPoint,
-    tenor::{DateOffset, TenorUnit},
     time::Instant,
 };
 use rhodium_market::{fx::FxForwardPlan, structure::MarketStructure};
@@ -48,27 +47,20 @@ pub(super) fn prepare(
     structure.fx_forward_plan(pair, id, value)
 }
 
+/// FX fixing の予想を読む受渡の瞬間。fixing 日からペアのスポットラグを両通貨の休日の和で
+/// 進めた日の、ペアの慣行の受渡時刻とする。
+///
+/// 同じ指標・観測日の公表値は一つなので、元本リセットの fixing もこの瞬間で読む。新しい元本を
+/// 受け渡す瞬間（[`FxFixing::value_at`]）は脚の通貨の決済時刻で、予想の受渡とは別に扱う。
 pub(super) fn value_time(
     fixing: FxFixing<'_>,
     ctx: &Context<'_>,
     warnings: &mut BTreeSet<Warning>,
-) -> Result<Instant> {
-    match fixing.value_at {
-        Some(at) => Ok(at),
-        None => spot_value(fixing, ctx, warnings),
-    }
-}
-
-fn spot_value(
-    fixing: FxFixing<'_>,
-    ctx: &Context<'_>,
-    warnings: &mut BTreeSet<Warning>,
 ) -> Result<Instant> {
     let pair = fixing.pair;
-    let definition = ctx.refdata.index(fixing.index)?;
-    let calendar = ctx.refdata.calendar(&definition.calendar)?;
-    let quote = ctx.refdata.currency(pair.quote())?;
-    let event = &quote.settlement_time;
+    let convention = ctx.refdata.fx_convention(pair)?;
+    let calendar = ctx.refdata.calendar(&ctx.refdata.fx_calendar_spec(pair)?)?;
+    let event = &convention.settlement;
     let zone = ctx.refdata.zone(&event.zone)?;
     let date = match fixing.point {
         ObservationPoint::Date(date) => date,
@@ -77,7 +69,7 @@ fn spot_value(
     };
     let (value_date, notices) = advance(
         date,
-        DateOffset::new(2, TenorUnit::BusinessDay),
+        convention.spot_lag.offset(),
         &AdvanceRule {
             calendar: &calendar,
             bdc: BusinessDayConvention::Following,
diff --git a/crates/rhodium-pricing/src/linear/rfr.rs b/crates/rhodium-pricing/src/linear/rfr.rs
index 94b6d51227aed699307291743489fcc930b7456a..4958cf9c2913994241dd34c02045dad517be5020 100644
--- a/crates/rhodium-pricing/src/linear/rfr.rs
+++ b/crates/rhodium-pricing/src/linear/rfr.rs
@@ -121,6 +121,15 @@ impl CompoundPlan {
     }
 }
 
+/// ターム金利の観測を受け渡す時点。
+#[derive(Debug, Clone, Copy)]
+enum TermPayment {
+    /// 観測期間の終わりに払うクーポン（`TermForward`）。
+    PeriodEnd,
+    /// 期首に割り引いて払う FRA の決済（`FraSettlement`）。
+    FraStart,
+}
+
 pub(super) struct Forecast<'a> {
     pub(super) curve: ResolvedCurve<Projection>,
     pub(super) rates: Vec<(CurvePoint, CurvePoint, YearFraction)>,
@@ -143,6 +152,10 @@ struct ProjectionRow {
 }
 
 impl<'a> Forecast<'a> {
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#197 で各 variant を網羅する"
+    )]
     pub(super) fn new(
         index: &IndexId,
         curve: ResolvedCurve<Projection>,
@@ -191,14 +204,7 @@ impl<'a> Forecast<'a> {
         ctx: &Context<'_>,
         warnings: &mut BTreeSet<Warning>,
     ) -> Result<usize> {
-        self.term_point_with_convention(
-            fixing,
-            discount,
-            pay,
-            ctx,
-            warnings,
-            LinearConvention::CrossCollateralProjection,
-        )
+        self.term_point_paid(fixing, discount, pay, ctx, warnings, TermPayment::PeriodEnd)
     }
 
     pub(super) fn fra_point(
@@ -209,24 +215,17 @@ impl<'a> Forecast<'a> {
         ctx: &Context<'_>,
         warnings: &mut BTreeSet<Warning>,
     ) -> Result<usize> {
-        self.term_point_with_convention(
-            fixing,
-            discount,
-            pay,
-            ctx,
-            warnings,
-            LinearConvention::FraDiscounting,
-        )
+        self.term_point_paid(fixing, discount, pay, ctx, warnings, TermPayment::FraStart)
     }
 
-    fn term_point_with_convention(
+    fn term_point_paid(
         &mut self,
         fixing: TermFixing<'_>,
         discount: &ResolvedCurve<Discount>,
         pay: Instant,
         ctx: &Context<'_>,
         warnings: &mut BTreeSet<Warning>,
-        convention: LinearConvention,
+        payment: TermPayment,
     ) -> Result<usize> {
         if fixing.publication <= ctx.as_of || fixing.period.end() <= ctx.as_of {
             return Err(Error::MissingData(format!(
@@ -239,13 +238,21 @@ impl<'a> Forecast<'a> {
                 "term projection in another payment currency".into(),
             ));
         }
-        if self.curve.spec().id() != discount.spec().id() {
-            warnings.insert(ctx.config.conventions().require(convention)?);
+        let conventions = ctx.config.conventions();
+        // 予想の担保（予想の測度）と支払いの担保が違うときだけ、担保の違いによる予想の調整を無視する。
+        if !self.curve.is_forecast_under(discount) {
+            warnings.insert(conventions.require(LinearConvention::CrossCollateralProjection)?);
         }
-        let natural = if convention == LinearConvention::FraDiscounting {
-            fixing.period.start()
-        } else {
-            fixing.period.end()
+        let natural = match payment {
+            TermPayment::PeriodEnd => fixing.period.end(),
+            // FRA の決済 N·(1 − (1 + τK)·P_p(e)/P_p(s))·P_d(s) は P_p = P_d のときだけ厳密で、
+            // 同じ担保でも予想曲線と割引曲線が違えば凸性を無視する。
+            TermPayment::FraStart => {
+                if self.curve.spec().id() != discount.spec().id() {
+                    warnings.insert(conventions.require(LinearConvention::FraDiscounting)?);
+                }
+                fixing.period.start()
+            }
         };
         let date = self.zone.to_local(natural)?.date();
         let delay = payment_delay(ctx, discount.spec().currency(), date, pay, warnings)?;
@@ -291,7 +298,8 @@ impl<'a> Forecast<'a> {
             return Ok(plan);
         };
         let conventions = ctx.config.conventions();
-        if self.curve.spec().id() != discount.spec().id() {
+        // 予想の担保（予想の測度）と支払いの担保が違うときだけ、担保の違いによる予想の調整を無視する。
+        if !self.curve.is_forecast_under(discount) {
             warnings.insert(conventions.require(LinearConvention::CrossCollateralProjection)?);
         }
         if rounding.is_some() {
@@ -406,27 +414,10 @@ impl<'a> Forecast<'a> {
         end: Date,
         warnings: &mut BTreeSet<Warning>,
     ) -> Result<YearFraction> {
-        let day = self.zone.to_local(at)?.date();
-        let cutoff = |date| {
-            self.zone
-                .resolve_convention(self.observation_time.on_date(date))
-        };
-        let start = if at < cutoff(day)? {
-            day.checked_add_days(-1)?
-        } else {
-            day
-        };
-        let next = start.checked_add_days(1)?;
-        let from = cutoff(start)?.unix_micros();
-        let to = cutoff(next)?.unix_micros();
-        let remaining = (to - at.unix_micros()) as f64 / (to - from) as f64;
-        // 観測時刻から翌日の観測時刻までの長さを使う。夏時間の開始・終了日は 24 時間と限らない。
-        let first = self.fraction(start, end, warnings)?.value();
-        let second = self.fraction(next, end, warnings)?.value();
-        YearFraction::new(
-            second + remaining * (first - second),
-            self.day_count.clone(),
-        )
+        let value = remaining_fraction(at, end, self.zone, self.observation_time, |from, to| {
+            Ok(self.fraction(from, to, warnings)?.value())
+        })?;
+        YearFraction::new(value, self.day_count.clone())
     }
 
     fn rate_point(
@@ -536,6 +527,35 @@ fn numerical(error: impl std::fmt::Display) -> Error {
     Error::Numerical(error.to_string())
 }
 
+/// 評価時点 `at` から観測の期間の終了日 `end` までの残りの日数係数。
+///
+/// `at` を含む暦日の観測時刻から翌日の観測時刻までを時間で按分し、その前後の暦日から `end` までの
+/// 日数係数を補間する。観測時刻から翌日の観測時刻までの長さを使うので、夏時間の開始・終了日
+/// （24 時間と限らない日）でも按分の比は 0 から 1 に収まる。`fraction(start, end)` は日付の間の
+/// 日数係数で、`start == end` のときは 0 を返すこと。
+pub(super) fn remaining_fraction(
+    at: Instant,
+    end: Date,
+    zone: &ZoneRules,
+    observation_time: &LocalEventTime,
+    mut fraction: impl FnMut(Date, Date) -> Result<f64>,
+) -> Result<f64> {
+    let day = zone.to_local(at)?.date();
+    let cutoff = |date| zone.resolve_convention(observation_time.on_date(date));
+    let start = if at < cutoff(day)? {
+        day.checked_add_days(-1)?
+    } else {
+        day
+    };
+    let next = start.checked_add_days(1)?;
+    let from = cutoff(start)?.unix_micros();
+    let to = cutoff(next)?.unix_micros();
+    let remaining = (to - at.unix_micros()) as f64 / (to - from) as f64;
+    let first = fraction(start, end)?;
+    let second = fraction(next, end)?;
+    Ok(second + remaining * (first - second))
+}
+
 pub(super) fn payment_delay(
     ctx: &Context<'_>,
     currency: rhodium_core::currency::Currency,
diff --git a/crates/rhodium-pricing/src/prepare.rs b/crates/rhodium-pricing/src/prepare.rs
index 5149949c0a8f3c14fc0f690870e088bb37d35cc3..7a37aec5b73d15bc08eacbe1ae684503de720f9e 100644
--- a/crates/rhodium-pricing/src/prepare.rs
+++ b/crates/rhodium-pricing/src/prepare.rs
@@ -15,19 +15,19 @@ use rhodium_products::{
         analysis::{Linearity, linearity},
         reduce,
     },
-    standard::StandardProduct,
-    trade::ExerciseRecord,
+    standard::{ProductKind, StandardProduct},
+    trade::{ExerciseRecord, ProductSpec, ResolvedTrade},
 };
 use scioxide_ad::Real;
 
 use crate::{
-    config::Method,
+    config::{AnalyticPricer, Method},
     context::Context,
     linear::{self, CashflowPlan},
     resolve::CollateralKey,
 };
 
-/// M0 のリニア評価に渡す商品。
+/// 準備する商品。方式の表は商品の種類（標準商品の `ProductKind`）で引く。
 #[derive(Debug, Clone, Copy)]
 pub enum PrepareInput<'a> {
     /// 標準商品と、日程の解決に使う約定日。
@@ -37,8 +37,25 @@ pub enum PrepareInput<'a> {
         /// 日程を生成する起点の約定日。
         trade_date: Date,
     },
-    /// 簡約前の検証済み IR。構築時と同じ参照データを評価にも使う。
+    /// 簡約前の検証済み IR。構築時と同じ参照データを評価にも使う。種類を持たない。
     Ir(&'a Product),
+    /// 取引と、同じ参照データでその取引を下ろした検証済み IR の組。
+    ///
+    /// 下ろした IR を作り直さずに使い、方式の表は取引の商品の種類で引く。解析プライサーは
+    /// 取引の標準商品、約定日と未決済の手数料を読む。組は [`ResolvedTrade`] だけが作るので、
+    /// 無関係な取引と IR を渡せない。
+    Trade(&'a ResolvedTrade),
+}
+
+impl PrepareInput<'_> {
+    /// 方式の表が照合する商品の種類。IR だけを渡した場合と、変更した取引は持たない。
+    fn kind(&self) -> Option<ProductKind> {
+        match self {
+            Self::Standard { product, .. } => Some(product.kind()),
+            Self::Ir(_) => None,
+            Self::Trade(resolved) => resolved.trade().product_kind(),
+        }
+    }
 }
 
 /// 市場のパラメータ θ に依存しない、評価に必要な計画。
@@ -133,36 +150,25 @@ pub fn prepare(
             "market structure and context have different valuation instants".into(),
         ));
     }
-    if let PrepareInput::Standard {
-        product: StandardProduct::RateFuture(future),
-        trade_date,
-    } = input
-        && let Method::Analytic(_) = ctx.config.select_method(
-            Some(rhodium_products::standard::ProductKind::RateFuture),
-            true,
-        )?
+    let kind = input.kind();
+    // 方式の表で解析プライサーを選んだ標準商品は、IR を簡約せずに型付きの商品から準備する。
+    // 表の照合はリニア評価の可否によらないので、ここでは可否を問わずに引く。
+    if kind.is_some()
+        && let Method::Analytic(pricer) = ctx.config.select_method(kind, Linearity::Linear)?
     {
-        let (margin, warnings) = future.margin(trade_date, ctx.refdata)?;
-        let plan = linear::prepare_rate_future(&margin, warnings, collateral, structure, ctx)?;
-        return Ok(Prepared {
-            as_of: ctx.as_of,
-            method: Method::Analytic(crate::config::AnalyticPricer::RateFuture),
-            slots: plan.slots,
-            plan: plan.plan,
-            blocks: plan.blocks,
-            warnings: plan.warnings,
-        });
+        return prepare_analytic(pricer, input, collateral, structure, ctx);
     }
     let lowered;
-    let (product, kind) = match input {
+    let product = match input {
         PrepareInput::Standard {
             product,
             trade_date,
         } => {
             lowered = product.to_product(trade_date, ctx.refdata)?;
-            (&lowered, Some(product.kind()))
+            &lowered
         }
-        PrepareInput::Ir(product) => (product, None),
+        PrepareInput::Ir(product) => product,
+        PrepareInput::Trade(resolved) => resolved.product(),
     };
     let simplified = reduce(
         product,
@@ -174,7 +180,7 @@ pub fn prepare(
     let applicable = linearity(&simplified);
     let method = ctx
         .config
-        .select_method(kind, applicable == Linearity::Linear)
+        .select_method(kind, applicable)
         .map_err(|error| match applicable {
             Linearity::NonLinear(reason) => error.context(reason),
             Linearity::Linear => error,
@@ -198,6 +204,58 @@ pub fn prepare(
     })
 }
 
+/// 解析プライサーの入力を型付きの標準商品から作る。プライサーを網羅して分岐し、
+/// 足したプライサーの入力を決めないままリニアの計画に落とさない。
+///
+/// 取引の入力では、約定日と未決済の手数料も読む（IR の評価と同じ受渡を数える）。
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#202 で各 variant を網羅する"
+)]
+fn prepare_analytic(
+    pricer: AnalyticPricer,
+    input: PrepareInput<'_>,
+    collateral: &CollateralKey,
+    structure: &Arc<MarketStructure>,
+    ctx: &Context<'_>,
+) -> Result<Prepared> {
+    match pricer {
+        AnalyticPricer::RateFuture => {
+            let (future, trade_date, fees) = match input {
+                PrepareInput::Standard {
+                    product: StandardProduct::RateFuture(future),
+                    trade_date,
+                } => (future, trade_date, &[][..]),
+                PrepareInput::Trade(resolved) => match &resolved.trade().product {
+                    ProductSpec::Standard(StandardProduct::RateFuture(future)) => {
+                        let trade = resolved.trade();
+                        (future, trade.trade_date, trade.fees.as_slice())
+                    }
+                    _ => return Err(rate_future_input()),
+                },
+                PrepareInput::Standard { .. } | PrepareInput::Ir(_) => {
+                    return Err(rate_future_input());
+                }
+            };
+            let (margin, warnings) = future.margin(trade_date, ctx.refdata)?;
+            let plan =
+                linear::prepare_rate_future(&margin, fees, warnings, collateral, structure, ctx)?;
+            Ok(Prepared {
+                as_of: ctx.as_of,
+                method: Method::Analytic(pricer),
+                slots: plan.slots,
+                plan: plan.plan,
+                blocks: plan.blocks,
+                warnings: plan.warnings,
+            })
+        }
+    }
+}
+
+fn rate_future_input() -> Error {
+    Error::Unsupported("the rate-future analytic pricer requires a standard rate future".into())
+}
+
 /// マージンセットの VM 条件から担保の鍵を導いて、価格評価を準備する。
 ///
 /// 支払いごとの割引規則は、支払い通貨と導出した担保の鍵で照合する。
diff --git a/crates/rhodium-pricing/src/resolve.rs b/crates/rhodium-pricing/src/resolve.rs
index d974a5113b0a87db703fc86f069d6db4431f0ed4..981ce8f97c39783fbed7be941447664804973e4c 100644
--- a/crates/rhodium-pricing/src/resolve.rs
+++ b/crates/rhodium-pricing/src/resolve.rs
@@ -208,16 +208,9 @@ fn deserialize_cash<'de, D: Deserializer<'de>>(
         .map_err(de::Error::custom)
 }
 
+/// 担保の付利の率を、取引の率と同じ [`Spread`] の書式（`"25bp"`、`"0.25%"`）で読む。
 fn parse_collateral_rate(text: &str) -> Result<f64> {
-    let (number, scale) = text
-        .strip_suffix("bp")
-        .map(|number| (number, 4))
-        .or_else(|| text.strip_suffix('%').map(|number| (number, 2)))
-        .ok_or_else(|| Error::InvalidData("cash collateral rate requires bp or %".into()))?;
-    // 十進のまま尺度を変え、単位の違う等価な値を同じ f64 に丸める。
-    format!("{number}e-{scale}")
-        .parse()
-        .map_err(|_| Error::InvalidData(format!("invalid cash collateral rate {text}")))
+    text.parse::<Spread>().map(Spread::value)
 }
 
 /// 商品の評価に必要な市場の要求。
@@ -253,6 +246,10 @@ pub enum ResolvedRequirement {
 /// `MissingData`。曲線の役割や通貨、無担保の資金調達の指定が矛盾すれば
 /// `InvalidInput`。指標定義が無ければ `InvalidData`、予想が未対応の指標なら
 /// `Unsupported`。
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#197 で各 variant を網羅する"
+)]
 pub fn resolve(
     requirement: &Requirement,
     config: &PricingConfig,
diff --git a/crates/rhodium-pricing/tests/calibration_instruments.rs b/crates/rhodium-pricing/tests/calibration_instruments.rs
new file mode 100644
index 0000000000000000000000000000000000000000..44dbce47249b110cc17f3b96696de6a5e1a3b10b
--- /dev/null
+++ b/crates/rhodium-pricing/tests/calibration_instruments.rs
@@ -0,0 +1,207 @@
+//! IRS と預金のテンプレートを校正商品に使い、校正と取引の評価が一致することを確かめる。
+//!
+//! ベーシスで校正した市場で IRS（預金）のパーを求め、それを quote にして曲線を校正し直すと、
+//! 同じ θ に戻る（解が一意なので、元の曲線が新しい方程式も満たす）。校正した IRS と預金を
+//! quote で約定した取引として評価すると、PV は校正の残差の許容誤差の内でゼロになる。
+
+use rhodium_core::{Error, Result, currency::Currency};
+use rhodium_market::{
+    market::Market,
+    quote::{QuoteId, QuoteKind, QuoteSet},
+};
+use rhodium_pricing::{
+    calibration::{CalibratedMarket, GroupOutcome, calibrate},
+    context::Context,
+    definition::{CalibrationTarget, ElementSelector, InstrumentBinding, QuoteTransform},
+    prepare::{PrepareInput, evaluate, prepare},
+    resolve::CollateralKey,
+};
+use rhodium_products::standard::InstrumentTemplate;
+use rhodium_testkit::assert_close;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、校正商品の検査に使わない関数がある"
+)]
+#[path = "fixtures/ois.rs"]
+pub mod ois;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、校正商品の検査に使わない関数がある"
+)]
+#[path = "fixtures/multicurve.rs"]
+pub mod multicurve;
+
+use multicurve::{MulticurveCase, all_multicurve_cases};
+use ois::{Node, fixture_quotes};
+
+const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
+const TOLERANCE: f64 = 1e-14;
+
+fn solved(market: &CalibratedMarket) -> Result<()> {
+    for group in market.groups() {
+        let GroupOutcome::Solved(solution) = group.outcome() else {
+            return Err(Error::Calibration(format!("{group:?}")));
+        };
+        assert!(solution.residual_norm() <= TOLERANCE);
+    }
+    Ok(())
+}
+
+/// テンプレートの PV が quote K について affine なので、K = 0 と 1 の PV からパーと傾きを求める。
+fn par_and_slope(
+    template: &InstrumentTemplate,
+    case: &MulticurveCase,
+    market: &Market<f64>,
+    ctx: &Context<'_>,
+) -> Result<(f64, f64)> {
+    let pv = |rate| {
+        let product = template.build(rate, ctx.refdata)?;
+        let plan = prepare(
+            PrepareInput::Standard {
+                product: &product,
+                trade_date: case.overnight.as_of,
+            },
+            &COLLATERAL,
+            &[],
+            market.structure(),
+            ctx,
+        )?;
+        Ok::<_, Error>(evaluate(&plan, market)?.get(case.overnight.currency))
+    };
+    let zero = pv(0.0)?;
+    let slope = pv(1.0)? - zero;
+    Ok((-zero / slope, slope))
+}
+
+fn irs(case: &MulticurveCase, node: &Node) -> Result<InstrumentTemplate> {
+    Ok(InstrumentTemplate::Swap {
+        convention: case.irs_convention().parse()?,
+        tenor: node.tenor,
+    })
+}
+
+fn assert_same_theta(lhs: &CalibratedMarket, rhs: &CalibratedMarket) {
+    // 二つの校正はそれぞれ残差 1e-14 で止まるので、θ はその誤差を伝えた幅で一致する。
+    for (lhs, rhs) in lhs.market().theta().iter().zip(rhs.market().theta()) {
+        assert_close!(*lhs, *rhs, 1e-10, 1e-12);
+    }
+}
+
+#[test]
+fn irs_quotes_calibrate_the_tibor_curve_and_reprice_at_par() -> Result<()> {
+    // JPY の TIBOR 6M と 3M（TONA へのベーシス）と USD のターム指標（独立の節点）で、ベーシスの代わりに
+    // 固定対ターム金利の IRS のパーレートでターム予想曲線を当てる。
+    for case in all_multicurve_cases()? {
+        let inputs = case.inputs()?;
+        let ctx = inputs.context();
+        let basis = calibrate(&case.definition()?, &case.quotes()?, &ctx)?;
+        solved(&basis)?;
+        let term_block = ElementSelector::Block(case.term_block()?);
+        let mut bindings = Vec::new();
+        let mut quotes = fixture_quotes(&case.overnight)?;
+        for node in &case.nodes {
+            let template = irs(&case, node)?;
+            let (par, _) = par_and_slope(&template, &case, basis.market(), &ctx)?;
+            let quote: QuoteId = format!("{}.IRS.{}", case.term_index, node.tenor).parse()?;
+            bindings.push(InstrumentBinding::new(
+                quote.clone(),
+                QuoteKind::Rate,
+                CalibrationTarget::ParRate(template),
+                QuoteTransform::Identity,
+                COLLATERAL,
+                term_block.clone(),
+            )?);
+            quotes.push((quote, par));
+        }
+        let quotes = QuoteSet::new(quotes)?;
+        let swaps = calibrate(
+            &case.definition_with_term_bindings(bindings)?,
+            &quotes,
+            &ctx,
+        )?;
+        solved(&swaps)?;
+        assert_same_theta(&basis, &swaps);
+        // 校正した IRS を quote の固定金利で約定した取引の PV は、残差の許容誤差 × 傾きの幅でゼロ。
+        for node in &case.nodes {
+            let template = irs(&case, node)?;
+            let id = format!("{}.IRS.{}", case.term_index, node.tenor).parse()?;
+            let quote = quotes
+                .get(&id)
+                .ok_or_else(|| Error::Internal("missing IRS quote".into()))?;
+            let (par, slope) = par_and_slope(&template, &case, swaps.market(), &ctx)?;
+            assert_close!(par, quote, 0.0, TOLERANCE);
+            assert_close!(
+                slope * (quote - par),
+                0.0,
+                0.0,
+                2.0 * slope.abs() * TOLERANCE
+            );
+        }
+    }
+    Ok(())
+}
+
+#[test]
+fn deposit_template_calibrates_and_reprices_at_par() -> Result<()> {
+    // TONA の最初のピラー（6M）を OIS の代わりにスポット起点の 6M の預金で当てる。
+    // 預金の満期は OIS の最終の支払日と同じなので、預金のパーを quote にすると同じ曲線に戻る。
+    let [case, _, _] = all_multicurve_cases()?;
+    let inputs = case.inputs()?;
+    let ctx = inputs.context();
+    let original = calibrate(&case.definition()?, &case.quotes()?, &ctx)?;
+    let first = &case.overnight.nodes[0];
+    let template = InstrumentTemplate::Deposit {
+        convention: case.deposit_convention().parse()?,
+        tenor: first.tenor.to_string().parse()?,
+        currency: Currency::JPY,
+    };
+    let (par, _) = par_and_slope(&template, &case, original.market(), &ctx)?;
+    let deposit: QuoteId = "JPY.DEPO.6M".parse()?;
+    let overnight_block = ElementSelector::Block(case.overnight.block()?);
+    let term_block = ElementSelector::Block(case.term_block()?);
+    let bindings = std::iter::once(InstrumentBinding::new(
+        deposit.clone(),
+        QuoteKind::Rate,
+        CalibrationTarget::ParRate(template.clone()),
+        QuoteTransform::Identity,
+        COLLATERAL,
+        overnight_block.clone(),
+    ))
+    .chain(case.overnight.nodes[1..].iter().map(|node| {
+        case.overnight
+            .binding(node, &COLLATERAL, overnight_block.clone())
+    }))
+    .chain(case.nodes.iter().map(|node| {
+        InstrumentBinding::new(
+            case.quote(node)?,
+            QuoteKind::Spread,
+            CalibrationTarget::ParRate(case.template(node)?),
+            QuoteTransform::Identity,
+            COLLATERAL,
+            term_block.clone(),
+        )
+    }))
+    .collect::<Result<Vec<_>>>()?;
+    let definition = rhodium_pricing::definition::MarketDefinition::new(
+        "JPY".parse()?,
+        "Etc/UTC".parse()?,
+        [case.term_curve()?, case.overnight.curve()?],
+        bindings,
+    )?;
+    let first_quote = case.overnight.quote(first)?;
+    let quotes = QuoteSet::new(
+        case.quotes()?
+            .iter()
+            .filter(|(id, _)| **id != first_quote)
+            .map(|(id, value)| (id.clone(), value))
+            .chain([(deposit, par)]),
+    )?;
+    let deposits = calibrate(&definition, &quotes, &ctx)?;
+    solved(&deposits)?;
+    assert_same_theta(&original, &deposits);
+    let (repriced, _) = par_and_slope(&template, &case, deposits.market(), &ctx)?;
+    assert_close!(repriced, par, 0.0, TOLERANCE);
+    Ok(())
+}
diff --git a/crates/rhodium-pricing/tests/config_and_resolve.rs b/crates/rhodium-pricing/tests/config_and_resolve.rs
index 1fd31f0c4fc0328b36154b031e524880a70f9a52..77bd98d964229d53d1f8becc007b127f461503f9 100644
--- a/crates/rhodium-pricing/tests/config_and_resolve.rs
+++ b/crates/rhodium-pricing/tests/config_and_resolve.rs
@@ -9,12 +9,13 @@ use rhodium_core::{
     date::Weekday,
     day_count::DayCount,
     fixing::FixingStore,
+    fx::CcyPair,
     id::{CalendarId, IndexId, SecurityId},
     quantity::{Rate, Simple, Spread, YearFraction},
     refdata::{
-        AssetValuation, CurrencyDefinition, DailyPublication, Deliverability, IndexDefinition,
-        IndexKind, IndexStatus, PublicationRule, RefData, RefDataParts, SecurityDefinition,
-        SecurityKind, SecurityUnit,
+        AssetValuation, CurrencyDefinition, DailyPublication, Deliverability, FxConvention,
+        FxSpotLag, IndexDefinition, IndexKind, IndexStatus, PublicationRule, RefData, RefDataParts,
+        SecurityDefinition, SecurityKind, SecurityUnit,
     },
     rounding::{Rounding, RoundingMode},
     time::{Instant, LocalEventTime, UtcOffset, ZoneRules},
@@ -36,7 +37,7 @@ use rhodium_products::{
         CashCollateralTerms, CollateralAsset, CsaTerms, Eligibility, EligibleCollateral, MarginSet,
         PerParty, Remuneration,
     },
-    ir::{Contract, Product, SettlementBoundary},
+    ir::{Contract, Product, SettlementBoundary, analysis::Linearity},
     standard::ProductKind,
 };
 use rhodium_testkit::assert_close;
@@ -139,6 +140,17 @@ fn fixture_refdata() -> Result<RefData> {
             },
         );
     }
+    // USDJPY の FX fixing の受渡日は、ペアの慣行のスポットラグで決まる。
+    let usd_jpy: CcyPair = "USDJPY".parse()?;
+    parts.fx_conventions.insert(
+        usd_jpy,
+        FxConvention {
+            pair: usd_jpy,
+            spot_lag: FxSpotLag::TwoBusinessDays,
+            settlement: parts.currencies[&Currency::JPY].settlement_time.clone(),
+            pip_exponent: 2,
+        },
+    );
     let fx: IndexId = "USDJPY-WMR".parse()?;
     parts.indices.insert(
         fx.clone(),
@@ -326,19 +338,22 @@ fn method_selection_checks_explicit_and_default_linear_applicability() -> Result
     // スワップには明示行があり、預金や標準商品に分類されない IR では既定を使う。
     let config = parse(CONFIG)?;
     assert_eq!(
-        config.select_method(Some(ProductKind::Swap), true)?,
+        config.select_method(Some(ProductKind::Swap), Linearity::Linear)?,
         Method::Linear
     );
     assert_eq!(
-        config.select_method(Some(ProductKind::Deposit), true)?,
+        config.select_method(Some(ProductKind::Deposit), Linearity::Linear)?,
         Method::Linear
     );
-    assert_eq!(config.select_method(None, true)?, Method::Linear);
-    assert!(
-        matches!(config.select_method(Some(ProductKind::Swap), false), Err(Error::Unsupported(error)) if error.contains("configured"))
+    assert_eq!(
+        config.select_method(None, Linearity::Linear)?,
+        Method::Linear
     );
     assert!(
-        matches!(config.select_method(Some(ProductKind::Deposit), false), Err(Error::Unsupported(error)) if error.contains("default"))
+        matches!(config.select_method(Some(ProductKind::Swap), Linearity::NonLinear("choice")), Err(Error::Unsupported(error)) if error.contains("configured"))
+    );
+    assert!(
+        matches!(config.select_method(Some(ProductKind::Deposit), Linearity::NonLinear("choice")), Err(Error::Unsupported(error)) if error.contains("default"))
     );
     Ok(())
 }
diff --git a/crates/rhodium-pricing/tests/fixtures/futures.rs b/crates/rhodium-pricing/tests/fixtures/futures.rs
new file mode 100644
index 0000000000000000000000000000000000000000..8946183d481f121f62d88fdd9818dfd500383ad6
--- /dev/null
+++ b/crates/rhodium-pricing/tests/fixtures/futures.rs
@@ -0,0 +1,252 @@
+//! 一つの USD 曲線、金利先物の清算指標と限月、LGM のマークからなる M1 の先物と FX の試験用の市場。
+//!
+//! 週末だけの暦（2025〜2026 年）、UTC の 16 時の公表と受渡、1W の USD ターム指標と、その指標を
+//! 原資産とする 2026 年 1 月限の先物（乗数 1,000、最終清算 2026-01-07）を持つ。
+
+use std::{collections::BTreeMap, sync::Arc};
+
+use rhodium_core::{
+    Error, Result,
+    calendar::{
+        BaseCalendar, BusinessDayConvention, CalendarSpec, EndOfMonthRule,
+        Extrapolation as CalendarExtrapolation, WeekendRule,
+    },
+    currency::Currency,
+    date::Weekday,
+    day_count::DayCount,
+    id::{CalendarId, IndexId, SecurityId},
+    refdata::{
+        AssetValuation, Convention, CurrencyDefinition, DailyPublication, Deliverability,
+        IndexDefinition, IndexKind, IndexStatus, PublicationRule, RateFutureContract,
+        RateFutureConvention, RateFutureRate, RefData, RefDataParts, SecurityDefinition,
+        SecurityKind, SecurityUnit,
+    },
+    rounding::{Rounding, RoundingMode},
+    tenor::DateOffset,
+    time::{Instant, LocalEventTime, TzId, UtcOffset, ZoneRules, model_time},
+};
+use rhodium_market::{
+    curve::{CurveRole, CurveSpec, Extrapolation},
+    market::Market,
+    quote::QuoteKind,
+    structure::{BlockId, MarketStructure},
+};
+use rhodium_pricing::{
+    config::PricingConfig,
+    definition::{
+        CalibrationTarget, ElementSelector, InstrumentBinding, LgmDefinition, MarketDefinition,
+        QuoteTransform,
+    },
+    resolve::CollateralKey,
+};
+use rhodium_products::standard::InstrumentTemplate;
+
+pub(crate) const START: &str = "2026-01-05";
+pub(crate) const AS_OF: &str = "2026-01-05T00:00:00Z";
+pub(crate) const FINAL: &str = "2026-01-07";
+pub(crate) const PILLAR: &str = "2026-05-08T16:00:00Z";
+pub(crate) const CLEARING: &str = "USD-FUT-JAN";
+pub(crate) const MULTIPLIER: f64 = 1_000.0;
+pub(crate) const RATE: f64 = 0.04;
+
+pub(crate) fn at(day: &str, time: &str) -> Result<Instant> {
+    format!("{day}T{time}:00Z").parse()
+}
+
+pub(crate) fn refdata(delay: DateOffset, settlement_time: &str) -> Result<RefData> {
+    RefData::new(reference_parts(delay, settlement_time)?)
+}
+
+#[expect(
+    clippy::too_many_lines,
+    reason = "暦、原指標、清算と限月のテスト用参照データをここにまとめる"
+)]
+pub(crate) fn reference_parts(delay: DateOffset, settlement_time: &str) -> Result<RefDataParts> {
+    let calendar: CalendarId = "USD".parse()?;
+    let zone: TzId = "Etc/UTC".parse()?;
+    let time = "16:00:00".parse()?;
+    let clock = LocalEventTime {
+        time,
+        zone: zone.clone(),
+    };
+    let mut parts = RefDataParts::default();
+    parts.zones.insert(
+        zone.clone(),
+        ZoneRules::new(
+            "2025-01-01T00:00:00Z".parse()?..="2027-01-01T00:00:00Z".parse()?,
+            UtcOffset::from_seconds(0)?,
+            vec![],
+        )?,
+    );
+    parts.calendars.insert(
+        calendar.clone(),
+        BaseCalendar::new(
+            "2025-01-01".parse()?..="2026-12-31".parse()?,
+            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
+            vec![],
+            CalendarExtrapolation::Strict,
+        )?,
+    );
+    parts.currencies.insert(
+        Currency::USD,
+        CurrencyDefinition {
+            code: Currency::USD,
+            minor_units: 2,
+            payment_rounding: Rounding::new(0.01, RoundingMode::HalfEven)?,
+            deliverability: Deliverability::Deliverable,
+            settlement_time: clock.clone(),
+            settlement_calendar: calendar.clone(),
+        },
+    );
+    let security: SecurityId = CLEARING.parse()?;
+    parts.securities.insert(
+        security.clone(),
+        SecurityDefinition {
+            id: security.clone(),
+            kind: SecurityKind::ListedFuture,
+            currency: Currency::USD,
+            unit: SecurityUnit::Contract,
+            issuer: None,
+            valuation: AssetValuation::Forward,
+        },
+    );
+    let calendar = CalendarSpec::Base(calendar);
+    for (name, kind) in [
+        (
+            "TERM",
+            IndexKind::TermRate {
+                ccy: Currency::USD,
+                day_count: DayCount::Act365F,
+                tenor: "1W".parse()?,
+                value_offset: "0D".parse()?,
+                bdc: BusinessDayConvention::Following,
+                end_of_month: EndOfMonthRule::Ignore,
+            },
+        ),
+        (
+            CLEARING,
+            IndexKind::RateFuturesSettlement {
+                contract: security,
+                underlying: "TERM".parse()?,
+            },
+        ),
+    ] {
+        let id: IndexId = name.parse()?;
+        parts.indices.insert(
+            id.clone(),
+            IndexDefinition {
+                id,
+                kind,
+                calendar: calendar.clone(),
+                status: IndexStatus::Active,
+                observation_time: Some(clock.clone()),
+                publication: PublicationRule::Daily(DailyPublication {
+                    business_days: 0,
+                    calendar: calendar.clone(),
+                    time,
+                    zone: zone.clone(),
+                }),
+            },
+        );
+    }
+    parts.conventions.insert(
+        "USD-RATE-FUTURES".parse()?,
+        Convention::RateFuture(Box::new(RateFutureConvention {
+            multiplier: MULTIPLIER,
+            settlement_lag: delay,
+            settlement_time: LocalEventTime {
+                time: settlement_time.parse()?,
+                zone,
+            },
+            margin_index: None,
+            contracts: BTreeMap::from([(
+                "2026-01".parse()?,
+                RateFutureContract {
+                    settlement_index: CLEARING.parse()?,
+                    final_clearing: FINAL.parse()?,
+                    rate: RateFutureRate::Term {
+                        fixing: FINAL.parse()?,
+                    },
+                },
+            )]),
+        })),
+    );
+    Ok(parts)
+}
+
+pub(crate) fn config(extra: &str) -> Result<PricingConfig> {
+    toml::from_str(&format!(
+        r#"
+version = 1
+reporting_currency = "USD"
+[[discounting]]
+currency = "USD"
+collateral = "uncollateralized"
+curve = "USD"
+[forecasting]
+TERM = "USD"
+{extra}
+"#
+    ))
+    .map_err(|error| Error::InvalidData(error.to_string()))
+}
+
+pub(crate) fn curve() -> Result<CurveSpec> {
+    CurveSpec::new(
+        "USD".parse()?,
+        Currency::USD,
+        [CurveRole::Discount, CurveRole::Projection],
+        [PILLAR.parse()?],
+        Extrapolation::Reject,
+    )
+}
+
+pub(crate) fn market(as_of: Instant, volatility: f64) -> Result<Market<f64>> {
+    let structure = Arc::new(MarketStructure::with_fx_and_lgm(
+        as_of,
+        [curve()?],
+        [],
+        [(Currency::USD, vec![])],
+    )?);
+    let curve = structure
+        .layout()
+        .range(&BlockId::nodes("USD".parse()?))
+        .ok_or_else(|| Error::Internal("missing curve block".into()))?;
+    let model = structure
+        .layout()
+        .range(&MarketStructure::lgm_block(Currency::USD)?)
+        .ok_or_else(|| Error::Internal("missing model block".into()))?;
+    let mut theta = vec![0.0; 3];
+    theta[curve.start] = -RATE * model_time(as_of, PILLAR.parse()?).years();
+    theta[model.start + 1] = volatility;
+    Market::from_theta(structure, theta)
+}
+
+/// 2026 年 1 月限の金利先物の校正テンプレート。
+pub(crate) fn futures_template() -> Result<InstrumentTemplate> {
+    Ok(InstrumentTemplate::RateFuture {
+        convention: "USD-RATE-FUTURES".parse()?,
+        delivery_month: "2026-01".parse()?,
+    })
+}
+
+/// USD 曲線の一つのピラーを先物の価格で校正し、凸性の補正に使う LGM の平均回帰とボラを
+/// 直接マークにする市場の定義。
+pub(crate) fn futures_definition() -> Result<MarketDefinition> {
+    let binding = InstrumentBinding::new(
+        "USD.FUT.2026-01".parse()?,
+        QuoteKind::FuturesPrice,
+        CalibrationTarget::Price(futures_template()?),
+        QuoteTransform::Identity,
+        CollateralKey::Uncollateralized,
+        ElementSelector::Block(BlockId::nodes("USD".parse()?)),
+    )?;
+    MarketDefinition::new("USD".parse()?, "Etc/UTC".parse()?, [curve()?], [binding])?.with_lgm([(
+        Currency::USD,
+        LgmDefinition {
+            breaks: vec![],
+            mean_reversion: vec!["USD.LGM.MR".parse()?],
+            volatility: vec!["USD.LGM.VOL".parse()?],
+        },
+    )])
+}
diff --git a/crates/rhodium-pricing/tests/fixtures/multicurve.rs b/crates/rhodium-pricing/tests/fixtures/multicurve.rs
index 1124095ccb17c7a0755f0a45a519ae9997f1a1fa..b53eb06eb050b95cb8916e25e86c2a3949c4bf5d 100644
--- a/crates/rhodium-pricing/tests/fixtures/multicurve.rs
+++ b/crates/rhodium-pricing/tests/fixtures/multicurve.rs
@@ -6,7 +6,10 @@ use rhodium_core::{
     Error, Result,
     calendar::{BusinessDayConvention, EndOfMonthRule},
     currency::Currency,
-    refdata::{Convention, IndexKind, SwapConvention, SwapCouponConvention, SwapLegConvention},
+    refdata::{
+        Convention, DepositConvention, IndexKind, SwapConvention, SwapCouponConvention,
+        SwapLegConvention,
+    },
     tenor::{DateOffset, TenorUnit},
 };
 use rhodium_market::{
@@ -35,6 +38,7 @@ struct QuantLibBasis {
     overnight_nodes: Vec<Node>,
     overnight_samples: Vec<Sample>,
     term_index: String,
+    term_tenor: String,
     nodes: Vec<Node>,
     samples: Vec<Sample>,
 }
@@ -45,19 +49,31 @@ pub struct MulticurveCase {
     pub overnight: Case,
     /// ターム指標とその予想曲線の ID。
     pub term_index: String,
+    /// ターム指標のテナー（ベーシスの両脚の支払いの頻度）。
+    pub term_tenor: String,
     /// ベーシス商品のテナーと spread、ターム予想曲線のピラー。
     pub nodes: Vec<Node>,
     /// ターム予想曲線のピラー以外の比較点。
     pub samples: Vec<Sample>,
 }
 
-/// JPY の TONA/TIBOR と USD の SOFR/テスト用ターム指標の入力を読む。
+/// JPY の TONA/TIBOR 6M と USD の SOFR/テスト用ターム指標の入力を読む。
 pub(crate) fn multicurve_cases() -> Result<[MulticurveCase; 2]> {
-    let [jpy, usd]: [QuantLibBasis; 2] =
+    let [jpy, usd, _] = all_multicurve_cases()?;
+    Ok([jpy, usd])
+}
+
+/// TIBOR 6M、USD のターム指標、TIBOR 3M の三つのベーシスの入力を読む。
+pub(crate) fn all_multicurve_cases() -> Result<[MulticurveCase; 3]> {
+    let [jpy, usd, jpy_3m]: [QuantLibBasis; 3] =
         serde_json::from_str(include_str!("quantlib/multicurve.json"))
             .map_err(|error| Error::InvalidData(error.to_string()))?;
     let [tona, sofr] = cases()?;
-    Ok([jpy.combine(tona), usd.combine(sofr)])
+    Ok([
+        jpy.combine(tona.clone()),
+        usd.combine(sofr),
+        jpy_3m.combine(tona),
+    ])
 }
 
 impl QuantLibBasis {
@@ -68,6 +84,7 @@ impl QuantLibBasis {
         MulticurveCase {
             overnight,
             term_index: self.term_index,
+            term_tenor: self.term_tenor,
             nodes: self.nodes,
             samples: self.samples,
         }
@@ -77,10 +94,15 @@ impl QuantLibBasis {
 impl MulticurveCase {
     /// OIS にターム指標とベーシス商品の慣行を加える。
     pub(crate) fn inputs(&self) -> Result<Inputs> {
+        self.inputs_with("")
+    }
+
+    /// [`Self::inputs`] の価格設定に、割引の行や近似の設定を加える。
+    pub(crate) fn inputs_with(&self, extra_config: &str) -> Result<Inputs> {
         Inputs::with_forecasts(
             &[&self.overnight],
             self.overnight.as_of,
-            "",
+            extra_config,
             |parts, forecasting| {
                 let overnight_id = self.overnight.index.parse()?;
                 let term_id = self.term_index.parse()?;
@@ -92,7 +114,7 @@ impl MulticurveCase {
                 term.kind = IndexKind::TermRate {
                     ccy: self.overnight.currency,
                     day_count: self.overnight.basis(),
-                    tenor: "6M".parse()?,
+                    tenor: self.term_tenor.parse()?,
                     value_offset: DateOffset::new(2, TenorUnit::BusinessDay),
                     bdc: BusinessDayConvention::ModifiedFollowing,
                     end_of_month: EndOfMonthRule::Ignore,
@@ -104,21 +126,44 @@ impl MulticurveCase {
                     coupon: SwapCouponConvention::Compounded(ois.float.clone()),
                     fx_reset: None,
                 };
-                // QuantLib のヘルパーと同じく、両脚を半年ごとにラグなしで払う。
-                pay.periods.frequency = "6M".parse()?;
+                // QuantLib のヘルパーと同じく、両脚をターム指標のテナーごとにラグなしで払う。
+                pay.periods.frequency = self.term_tenor.parse()?;
                 pay.payment.offset = DateOffset::new(0, TenorUnit::BusinessDay);
                 let mut receive = pay.clone();
                 receive.coupon = SwapCouponConvention::Term {
                     index: term.id.clone(),
                 };
+                // 固定対ターム金利の IRS（固定を払い、ターム金利を受ける。両脚はターム指標のテナーごと）。
+                let mut fixed = pay.clone();
+                fixed.coupon = SwapCouponConvention::Fixed;
+                let spot_lag = ois.spot_lag;
+                // スポット起点の預金（TONA の暦と日数基準）。
+                let deposit = DepositConvention {
+                    spot_lag,
+                    day_count: self.overnight.basis(),
+                    bdc: BusinessDayConvention::ModifiedFollowing,
+                    end_of_month: EndOfMonthRule::Ignore,
+                    calendar: ois.fixed.calendar.clone(),
+                };
                 parts.indices.insert(term.id.clone(), term);
                 parts.conventions.insert(
                     self.term_index.parse()?,
                     Convention::BasisSwap(Box::new(SwapConvention {
-                        spot_lag: ois.spot_lag,
-                        legs: [pay, receive],
+                        spot_lag,
+                        legs: [pay, receive.clone()],
                     })),
                 );
+                parts.conventions.insert(
+                    self.irs_convention().parse()?,
+                    Convention::Irs(Box::new(SwapConvention {
+                        spot_lag,
+                        legs: [fixed, receive],
+                    })),
+                );
+                parts.conventions.insert(
+                    self.deposit_convention().parse()?,
+                    Convention::Deposit(deposit),
+                );
                 writeln!(
                     forecasting,
                     "\"{}\" = \"{}\"",
@@ -129,6 +174,16 @@ impl MulticurveCase {
         )
     }
 
+    /// 固定対ターム金利の IRS の慣行の ID。
+    pub(crate) fn irs_convention(&self) -> String {
+        format!("{}-IRS", self.term_index)
+    }
+
+    /// スポット起点の預金の慣行の ID。
+    pub(crate) fn deposit_convention(&self) -> String {
+        format!("{}-DEPO", self.overnight.index)
+    }
+
     /// ターム予想曲線の校正対象。
     pub(crate) fn term_block(&self) -> Result<BlockId> {
         BlockId::new(self.term_index.parse()?, "basis")
@@ -139,7 +194,7 @@ impl MulticurveCase {
         format!("{}.BASIS.{}", self.term_index, node.tenor).parse()
     }
 
-    fn template(&self, node: &Node) -> Result<InstrumentTemplate> {
+    pub(crate) fn template(&self, node: &Node) -> Result<InstrumentTemplate> {
         Ok(InstrumentTemplate::BasisSwap {
             convention: self.term_index.parse()?,
             tenor: node.tenor,
@@ -148,6 +203,72 @@ impl MulticurveCase {
 
     /// 共通の OIS 曲線と、それを上流に持つターム予想曲線を定義する。
     pub(crate) fn definition(&self) -> Result<MarketDefinition> {
+        self.definition_with([])
+    }
+
+    /// [`Self::definition`] に、校正しない曲線（担保条件の違う割引曲線など）を加える。
+    pub(crate) fn definition_with(
+        &self,
+        extra: impl IntoIterator<Item = CurveSpec>,
+    ) -> Result<MarketDefinition> {
+        self.definition_with_collateral(extra, |_| None)
+    }
+
+    /// 共通の OIS 曲線に、ターム予想曲線のブロックを当てる任意の binding を組み合わせる。
+    pub(crate) fn definition_with_term_bindings(
+        &self,
+        term_bindings: impl IntoIterator<Item = InstrumentBinding>,
+    ) -> Result<MarketDefinition> {
+        let overnight_block = ElementSelector::Block(self.overnight.block()?);
+        let bindings = term_bindings
+            .into_iter()
+            .map(Ok)
+            .chain(self.overnight.nodes.iter().map(|node| {
+                self.overnight
+                    .binding(node, &COLLATERAL, overnight_block.clone())
+            }))
+            .collect::<Result<Vec<_>>>()?;
+        MarketDefinition::new(
+            self.overnight.currency.to_string().parse()?,
+            "Etc/UTC".parse()?,
+            [self.term_curve()?, self.overnight.curve()?],
+            bindings,
+        )
+    }
+
+    /// ターム予想曲線。JPY は OIS 曲線とベーシスの和、USD は独立した節点の曲線。
+    pub(crate) fn term_curve(&self) -> Result<CurveSpec> {
+        let term_nodes = CurveExpr::Nodes {
+            part: "basis".into(),
+            pillars: self
+                .nodes
+                .iter()
+                .map(|node| midnight(node.date))
+                .collect::<Result<_>>()?,
+            extrapolation: Extrapolation::Reject,
+        };
+        let term_curve = match self.overnight.currency {
+            Currency::JPY => CurveExpr::Sum(vec![
+                CurveExpr::Ref(self.overnight.index.parse()?),
+                term_nodes,
+            ]),
+            _ => term_nodes,
+        };
+        CurveSpec::from_expr(
+            self.term_index.parse()?,
+            self.overnight.currency,
+            [CurveRole::Projection],
+            term_curve,
+        )
+    }
+
+    /// [`Self::definition_with`] で、`collateral` が返すベーシスの行（ノードの番号）だけ
+    /// 校正商品の担保を変える。
+    pub(crate) fn definition_with_collateral(
+        &self,
+        extra: impl IntoIterator<Item = CurveSpec>,
+        collateral: impl Fn(usize) -> Option<CollateralKey>,
+    ) -> Result<MarketDefinition> {
         let overnight = self.overnight.curve()?;
         let term_nodes = CurveExpr::Nodes {
             part: "basis".into(),
@@ -176,13 +297,14 @@ impl MulticurveCase {
         let bindings = self
             .nodes
             .iter()
-            .map(|node| {
+            .enumerate()
+            .map(|(index, node)| {
                 InstrumentBinding::new(
                     self.quote(node)?,
                     QuoteKind::Spread,
                     CalibrationTarget::ParRate(self.template(node)?),
                     QuoteTransform::Identity,
-                    COLLATERAL,
+                    collateral(index).unwrap_or(COLLATERAL),
                     term_block.clone(),
                 )
             })
@@ -194,7 +316,7 @@ impl MulticurveCase {
         MarketDefinition::new(
             self.overnight.currency.to_string().parse()?,
             "Etc/UTC".parse()?,
-            [term, overnight],
+            [term, overnight].into_iter().chain(extra),
             bindings,
         )
     }
diff --git a/crates/rhodium-pricing/tests/fixtures/ois.rs b/crates/rhodium-pricing/tests/fixtures/ois.rs
index a94eee260e90c4b024b1d59744cf0a3dec6cbb29..864630f938b776d58e4e734bcef8b6b0524ea5b6 100644
--- a/crates/rhodium-pricing/tests/fixtures/ois.rs
+++ b/crates/rhodium-pricing/tests/fixtures/ois.rs
@@ -12,6 +12,7 @@ use rhodium_core::{
     date::{Date, Weekday},
     day_count::DayCount,
     fixing::{FixingStore, FixingStoreParts, ObservationPoint},
+    fx::CcyPair,
     id::{CalendarId, IndexId},
     observation::{
         Aggregation, CompoundingForm, ObservationDays, ObservationMethod, ObservationRounding,
@@ -19,8 +20,8 @@ use rhodium_core::{
     },
     refdata::{
         Convention, CurrencyDefinition, DailyPublication, Deliverability, FixedLegConvention,
-        IndexDefinition, IndexKind, IndexStatus, OisConvention, OvernightLegConvention,
-        PublicationRule, RefData, RefDataParts,
+        FxConvention, FxSpotLag, IndexDefinition, IndexKind, IndexStatus, OisConvention,
+        OvernightLegConvention, PublicationRule, RefData, RefDataParts,
     },
     rounding::{Rounding, RoundingMode},
     schedule::{PaymentGrouping, PaymentRule, PeriodAnchor},
@@ -332,6 +333,26 @@ impl Inputs {
             writeln!(forecasting, "\"{}\" = \"{}\"", case.index, case.index)
                 .map_err(format_error)?;
         }
+        // 両通貨の曲線があれば USDJPY の受渡慣行（T+2、16:00 UTC、pip は 0.01 円）を置き、
+        // 評価時点の換算レートと FX のテストが参照データから読めるようにする。
+        if [Currency::JPY, Currency::USD]
+            .iter()
+            .all(|ccy| parts.currencies.contains_key(ccy))
+        {
+            let usd_jpy: CcyPair = "USDJPY".parse()?;
+            parts.fx_conventions.insert(
+                usd_jpy,
+                FxConvention {
+                    pair: usd_jpy,
+                    spot_lag: FxSpotLag::TwoBusinessDays,
+                    settlement: LocalEventTime {
+                        time: "16:00:00".parse()?,
+                        zone: "Etc/UTC".parse()?,
+                    },
+                    pip_exponent: 2,
+                },
+            );
+        }
         add_forecasts(&mut parts, &mut forecasting)?;
         let mut config = String::from("version = 1\nreporting_currency = \"JPY\"\n");
         for (ccy, curve) in discounting {
diff --git a/crates/rhodium-pricing/tests/fixtures/quantlib/README.md b/crates/rhodium-pricing/tests/fixtures/quantlib/README.md
index 1d417ee6a5a48af2e773c3b45166af6b0646dfd2..5cf163a8ddf247d8d299d1319f527c8a100b4d9c 100644
--- a/crates/rhodium-pricing/tests/fixtures/quantlib/README.md
+++ b/crates/rhodium-pricing/tests/fixtures/quantlib/README.md
@@ -65,13 +65,15 @@ rhodium は最初のピラーで会合のステップを切り替え、二つ目
 uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py --turn > crates/rhodium-pricing/tests/fixtures/quantlib/ois_turn.json
 ```
 
-## TONA と TIBOR 6M、SOFR と USD ターム指標の複数曲線
+## TONA と TIBOR 6M・3M、SOFR と USD ターム指標の複数曲線
 
 `multicurve.json` は、QuantLib 1.42 の `OISRateHelper` で OIS 曲線を先に校正し、`OvernightIborBasisSwapRateHelper` で別のターム予想曲線を校正した結果である。
 OIS とベーシスには 6M、1Y、18M、2Y を使う。
 ベーシス商品は 2 営業日後に開始し、半年ごとにラグなしで支払う。
 翌日物脚に複利の外で spread を足して払い、ターム脚を受ける。
-ターム指標は JPY では QuantLib の TIBOR 6M、USD では SOFR のカレンダーを使うテスト用の 6M 指標で、予想期間の端点は指標の fixing 規則で決める。
+ターム指標は JPY では QuantLib の TIBOR 6M と TIBOR 3M、USD では SOFR のカレンダーを使うテスト用の 6M 指標で、予想期間の端点は指標の fixing 規則で決める。
+両脚の支払いの頻度はターム指標のテナー（TIBOR 3M は四半期ごと）で、出力の `term_tenor` に記録する。
+配列の順は TIBOR 6M、USD のターム指標、TIBOR 3M（三つ目は M1 で足した）。
 OIS 曲線を両脚の割引と翌日物の予想に使う。
 
 JPY は割引とターム予想のピラーがそろっており、rhodium ではターム予想を TONA とベーシスの和で表す。
@@ -82,3 +84,20 @@ Rust のテストはピラーと補間点の割引係数および予想の擬似
 ```sh
 uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py --multicurve > crates/rhodium-pricing/tests/fixtures/quantlib/multicurve.json
 ```
+
+## 金利ボラの SABR とキャップレット
+
+`volatility.json` は QuantLib 1.42 の `shiftedSabrVolatility` で作った正規・対数正規の Hagan ボラと、
+`bachelierBlackFormula` で値付けした合成キャップのフラットボラを持つ。
+SABR はシフトの有無、ATM、ATM の両側の近接点とウィングを照合する。
+キャップは Japan カレンダー、2025-01-07 評価、2 営業日スポット、ACT/365F、予想ゼロ、割引係数 1 とする。
+半年払いのタームは最初の期を除き、年払いの後決め RFR は最初の期を含める。
+RFR の分散時間は QuantLib へ渡す前に観測区間から計算し、各キャップに既知の異なる区間ボラを使って価格を出す。
+フラットボラは、その価格を再現する値を二分法で逆算する。
+
+```sh
+uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-pricing/tests/fixtures/quantlib/generate_volatility.py > crates/rhodium-pricing/tests/fixtures/quantlib/volatility.json
+```
+
+Rust ではボラと元本 1 のキャップ価格を相対 1e-12 と絶対 1e-14 の幅で照合する。
+許容誤差は両実装の演算順と逆算の停止判定の差を含む。
diff --git a/crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py b/crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py
index 9ed141f5fe7bb5f868e4afed8aed32c4a53bd627..c89e4f4417b6fe66cc3fb8052e9f8ed9bc97632e 100644
--- a/crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py
+++ b/crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py
@@ -149,7 +149,8 @@ def case(currency, name, index, payment_lag, rates, seasoned, turn=None, tenors=
     }, curve
 
 
-def multicurve_case(currency, overnight_id, index, term_id, term, lag, rates, seasoned, spreads):
+def multicurve_case(currency, overnight_id, index, term_id, term, lag, rates, seasoned, spreads,
+                    term_tenor="6M"):
     """OIS で割引と翌日物の予想を先に解き、ベーシスでターム指標の予想を解く。"""
     tenors = TENORS[1:]
     overnight, discount = case(currency, overnight_id, index, lag, rates, seasoned, tenors=tenors)
@@ -175,6 +176,7 @@ def multicurve_case(currency, overnight_id, index, term_id, term, lag, rates, se
         "overnight_nodes": overnight["nodes"],
         "overnight_samples": overnight["samples"],
         "term_index": term_id,
+        "term_tenor": term_tenor,
         "nodes": [
             {"tenor": tenor, "rate": spread, "date": day.ISO(), "df": df}
             for tenor, spread, (day, df) in zip(tenors, spreads, nodes, strict=True)
@@ -216,6 +218,10 @@ elif "--multicurve" in sys.argv:
         multicurve_case("USD", "USD-SOFR", sofr, "USD-TERM-6M", usd_term, 2,
                         [0.042, 0.041, 0.040, 0.039], (0.0420, 0.0010, 0.0455),
                         [-0.0007, -0.0002, 0.0004, 0.0008]),
+        # TIBOR 3M も TONA へのベーシスで、両脚を四半期ごとに払う。
+        multicurve_case("JPY", "JPY-TONA", tona, "JPY-TIBOR-3M", ql.Tibor(ql.Period("3M")), 0,
+                        tona_rates[1:], tona_seasoned, [0.0008, 0.0015, 0.0011, 0.0019],
+                        term_tenor="3M"),
     ], indent=2))
 else:
     print(json.dumps([
diff --git a/crates/rhodium-pricing/tests/fixtures/quantlib/generate_volatility.py b/crates/rhodium-pricing/tests/fixtures/quantlib/generate_volatility.py
new file mode 100644
index 0000000000000000000000000000000000000000..5d402c7e99860b2a50c5f5d7a6ad13cbb09670fb
--- /dev/null
+++ b/crates/rhodium-pricing/tests/fixtures/quantlib/generate_volatility.py
@@ -0,0 +1,75 @@
+"""QuantLib 1.42 で SABR と合成キャップの校正用ボラを作る。"""
+
+import json
+
+import QuantLib as ql
+
+
+def sabr_cases():
+    """ATM、ATM 近傍、ウィングを Hagan の 3 種類のボラで計算する。"""
+    for kind, forward, shift, formula, alpha in [
+        ("normal", 0.034, 0.0, ql.Normal, 0.028),
+        ("lognormal", 0.034, 0.0, ql.ShiftedLognormal, 0.035),
+        ("shifted_lognormal", -0.004, 0.024, ql.ShiftedLognormal, 0.027),
+        ("shifted_normal", -0.004, 0.024, ql.Normal, 0.027),
+    ]:
+        beta, nu, rho, time = 0.63, 0.52, -0.37, 2.4
+        for offset in [-0.013, -1e-10, 0.0, 1e-10, 0.009]:
+            strike = forward + offset
+            yield dict(
+                kind=kind, forward=forward, strike=strike, shift=shift, time=time,
+                beta=beta, alpha=alpha, nu=nu, rho=rho,
+                vol=ql.shiftedSabrVolatility(strike, forward, time, alpha, beta, nu, rho, shift, formula),
+            )
+
+
+def cap_cases():
+    """割引係数 1、フォワード 0 のキャップを QuantLib の Bachelier で値付けする。"""
+    cal = ql.Japan()
+    as_of = ql.Date(7, 1, 2025)
+    spot = cal.advance(as_of, 2, ql.Days)
+    for label, frequency, maturities in [
+        ("term", 6, [12, 18]),
+        ("overnight", 12, [12, 24]),
+    ]:
+        ends = [cal.advance(spot, i, ql.Months, ql.ModifiedFollowing)
+                for i in range(frequency, maturities[-1] + 1, frequency)]
+        schedule = list(zip([spot] + ends[:-1], ends))
+        if label == "term":
+            schedule = schedule[1:]
+        for strike in [0.002, 0.0]:
+            nodes = []
+            for i, (start, end) in enumerate(schedule):
+                if label == "term":
+                    expiry = cal.advance(start, -2, ql.Days)
+                    time = (expiry - as_of) / 365
+                else:
+                    expiry = end
+                    time = (start - as_of + (end - start) / 3) / 365
+                annuity = (end - start) / 365
+                vol = ([0.006, 0.011] if strike else [0.008, 0.009])[i]
+                nodes.append(dict(expiry=expiry.ISO(), annuity=annuity, time=time, vol=vol))
+
+            def cap_value(nodes, flat=None):
+                return sum(
+                    ql.bachelierBlackFormula(ql.Option.Call, strike, 0.0,
+                                             (n["vol"] if flat is None else flat) * n["time"]**0.5,
+                                             n["annuity"])
+                    for n in nodes
+                )
+
+            quotes = []
+            for i, maturity in enumerate(maturities, 1):
+                pv = cap_value(nodes[:i])
+                lo, hi = 0.0, 0.1
+                for _ in range(64):
+                    mid = (lo + hi) / 2
+                    if cap_value(nodes[:i], mid) < pv:
+                        lo = mid
+                    else:
+                        hi = mid
+                quotes.append(dict(maturity=f"{maturity}M", vol=(lo + hi) / 2, pv=pv))
+            yield dict(kind=label, strike=strike, nodes=nodes, quotes=quotes)
+
+
+print(json.dumps(dict(sabr=list(sabr_cases()), caps=list(cap_cases())), indent=2))
diff --git a/crates/rhodium-pricing/tests/fixtures/quantlib/multicurve.json b/crates/rhodium-pricing/tests/fixtures/quantlib/multicurve.json
index 1138a1c573c015859b917cdfeca806cb6c2f4176..a2ee64799d9b36ac9c538781a330cfefee6da652 100644
--- a/crates/rhodium-pricing/tests/fixtures/quantlib/multicurve.json
+++ b/crates/rhodium-pricing/tests/fixtures/quantlib/multicurve.json
@@ -49,6 +49,7 @@
       }
     ],
     "term_index": "JPY-TIBOR-6M",
+    "term_tenor": "6M",
     "nodes": [
       {
         "tenor": "6M",
@@ -148,6 +149,7 @@
       }
     ],
     "term_index": "USD-TERM-6M",
+    "term_tenor": "6M",
     "nodes": [
       {
         "tenor": "6M",
@@ -196,5 +198,105 @@
         "df": 0.9318918355727929
       }
     ]
+  },
+  {
+    "overnight_nodes": [
+      {
+        "tenor": "6M",
+        "rate": 0.002,
+        "date": "2025-07-09",
+        "df": 0.9989982592720291
+      },
+      {
+        "tenor": "1Y",
+        "rate": 0.004,
+        "date": "2026-01-09",
+        "df": 0.9960050264807817
+      },
+      {
+        "tenor": "18M",
+        "rate": 0.005,
+        "date": "2026-07-09",
+        "df": 0.992549331641807
+      },
+      {
+        "tenor": "2Y",
+        "rate": 0.007,
+        "date": "2027-01-12",
+        "df": 0.9860578741851718
+      }
+    ],
+    "overnight_samples": [
+      {
+        "date": "2025-01-24",
+        "df": 0.9999068998111723
+      },
+      {
+        "date": "2025-04-15",
+        "df": 0.9994634238175699
+      },
+      {
+        "date": "2025-08-12",
+        "df": 0.9984444851883086
+      },
+      {
+        "date": "2026-04-15",
+        "df": 0.994170676035945
+      },
+      {
+        "date": "2026-10-15",
+        "df": 0.9891420796400195
+      }
+    ],
+    "term_index": "JPY-TIBOR-3M",
+    "term_tenor": "3M",
+    "nodes": [
+      {
+        "tenor": "6M",
+        "rate": 0.0008,
+        "date": "2025-07-09",
+        "df": 0.9985978836829034
+      },
+      {
+        "tenor": "1Y",
+        "rate": 0.0015,
+        "date": "2026-01-09",
+        "df": 0.9945091015080203
+      },
+      {
+        "tenor": "18M",
+        "rate": 0.0011,
+        "date": "2026-07-09",
+        "df": 0.9909165069568865
+      },
+      {
+        "tenor": "2Y",
+        "rate": 0.0019,
+        "date": "2027-01-12",
+        "df": 0.9822993956366051
+      }
+    ],
+    "samples": [
+      {
+        "date": "2025-01-24",
+        "df": 0.9998696658525201
+      },
+      {
+        "date": "2025-04-15",
+        "df": 0.9992488952271131
+      },
+      {
+        "date": "2025-08-12",
+        "df": 0.997841083751653
+      },
+      {
+        "date": "2026-04-15",
+        "df": 0.9926020178702298
+      },
+      {
+        "date": "2026-10-15",
+        "df": 0.9863912022439585
+      }
+    ]
   }
 ]
diff --git a/crates/rhodium-pricing/tests/fixtures/quantlib/volatility.json b/crates/rhodium-pricing/tests/fixtures/quantlib/volatility.json
new file mode 100644
index 0000000000000000000000000000000000000000..e1256a2a4fc8dde1476a4223cc314f7dc83e13e4
--- /dev/null
+++ b/crates/rhodium-pricing/tests/fixtures/quantlib/volatility.json
@@ -0,0 +1,366 @@
+{
+  "sabr": [
+    {
+      "kind": "normal",
+      "forward": 0.034,
+      "strike": 0.021000000000000005,
+      "shift": 0.0,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.028,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.004971786612410266
+    },
+    {
+      "kind": "normal",
+      "forward": 0.034,
+      "strike": 0.0339999999,
+      "shift": 0.0,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.028,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.0034430546895146516
+    },
+    {
+      "kind": "normal",
+      "forward": 0.034,
+      "strike": 0.034,
+      "shift": 0.0,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.028,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.003443054682763394
+    },
+    {
+      "kind": "normal",
+      "forward": 0.034,
+      "strike": 0.034000000100000004,
+      "shift": 0.0,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.028,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.0034430546760121357
+    },
+    {
+      "kind": "normal",
+      "forward": 0.034,
+      "strike": 0.043000000000000003,
+      "shift": 0.0,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.028,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.0038045960184543317
+    },
+    {
+      "kind": "lognormal",
+      "forward": 0.034,
+      "strike": 0.021000000000000005,
+      "shift": 0.0,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.035,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.21077457381858106
+    },
+    {
+      "kind": "lognormal",
+      "forward": 0.034,
+      "strike": 0.0339999999,
+      "shift": 0.0,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.035,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.12649142258690763
+    },
+    {
+      "kind": "lognormal",
+      "forward": 0.034,
+      "strike": 0.034,
+      "shift": 0.0,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.035,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.12649142222600424
+    },
+    {
+      "kind": "lognormal",
+      "forward": 0.034,
+      "strike": 0.034000000100000004,
+      "shift": 0.0,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.035,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.12649142186510084
+    },
+    {
+      "kind": "lognormal",
+      "forward": 0.034,
+      "strike": 0.043000000000000003,
+      "shift": 0.0,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.035,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.11906366107410718
+    },
+    {
+      "kind": "shifted_lognormal",
+      "forward": -0.004,
+      "strike": -0.017,
+      "shift": 0.024,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.027,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.3037451524439348
+    },
+    {
+      "kind": "shifted_lognormal",
+      "forward": -0.004,
+      "strike": -0.0040000001,
+      "shift": 0.024,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.027,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.11880679923942589
+    },
+    {
+      "kind": "shifted_lognormal",
+      "forward": -0.004,
+      "strike": -0.004,
+      "shift": 0.024,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.027,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.11880679863263671
+    },
+    {
+      "kind": "shifted_lognormal",
+      "forward": -0.004,
+      "strike": -0.0039999999000000005,
+      "shift": 0.024,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.027,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.11880679802584754
+    },
+    {
+      "kind": "shifted_lognormal",
+      "forward": -0.004,
+      "strike": 0.004999999999999999,
+      "shift": 0.024,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.027,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.12443715598354288
+    },
+    {
+      "kind": "shifted_normal",
+      "forward": -0.004,
+      "strike": -0.017,
+      "shift": 0.024,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.027,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.0037541994217355173
+    },
+    {
+      "kind": "shifted_normal",
+      "forward": -0.004,
+      "strike": -0.0040000001,
+      "shift": 0.024,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.027,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.0023731092187857534
+    },
+    {
+      "kind": "shifted_normal",
+      "forward": -0.004,
+      "strike": -0.004,
+      "shift": 0.024,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.027,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.002373109212603822
+    },
+    {
+      "kind": "shifted_normal",
+      "forward": -0.004,
+      "strike": -0.0039999999000000005,
+      "shift": 0.024,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.027,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.0023731092064218917
+    },
+    {
+      "kind": "shifted_normal",
+      "forward": -0.004,
+      "strike": 0.004999999999999999,
+      "shift": 0.024,
+      "time": 2.4,
+      "beta": 0.63,
+      "alpha": 0.027,
+      "nu": 0.52,
+      "rho": -0.37,
+      "vol": 0.003010768048874448
+    }
+  ],
+  "caps": [
+    {
+      "kind": "term",
+      "strike": 0.002,
+      "nodes": [
+        {
+          "expiry": "2025-07-07",
+          "annuity": 0.5041095890410959,
+          "time": 0.4958904109589041,
+          "vol": 0.006
+        },
+        {
+          "expiry": "2026-01-07",
+          "annuity": 0.4958904109589041,
+          "time": 1.0,
+          "vol": 0.011
+        }
+      ],
+      "quotes": [
+        {
+          "maturity": "12M",
+          "vol": 0.006,
+          "pv": 0.00043907486933762453
+        },
+        {
+          "maturity": "18M",
+          "vol": 0.008979871180954953,
+          "pv": 0.0021552032520583365
+        }
+      ]
+    },
+    {
+      "kind": "term",
+      "strike": 0.0,
+      "nodes": [
+        {
+          "expiry": "2025-07-07",
+          "annuity": 0.5041095890410959,
+          "time": 0.4958904109589041,
+          "vol": 0.008
+        },
+        {
+          "expiry": "2026-01-07",
+          "annuity": 0.4958904109589041,
+          "time": 1.0,
+          "vol": 0.009
+        }
+      ],
+      "quotes": [
+        {
+          "maturity": "12M",
+          "vol": 0.008,
+          "pv": 0.0011329685815496468
+        },
+        {
+          "maturity": "18M",
+          "vol": 0.008582795886875414,
+          "pv": 0.002913453443943986
+        }
+      ]
+    },
+    {
+      "kind": "overnight",
+      "strike": 0.002,
+      "nodes": [
+        {
+          "expiry": "2026-01-09",
+          "annuity": 1.0,
+          "time": 0.33881278538812787,
+          "vol": 0.006
+        },
+        {
+          "expiry": "2027-01-12",
+          "annuity": 1.0082191780821919,
+          "time": 1.3415525114155251,
+          "vol": 0.011
+        }
+      ],
+      "quotes": [
+        {
+          "maturity": "12M",
+          "vol": 0.006,
+          "pv": 0.0006157037922225922
+        },
+        {
+          "maturity": "24M",
+          "vol": 0.009434634591593306,
+          "pv": 0.004795112935279358
+        }
+      ]
+    },
+    {
+      "kind": "overnight",
+      "strike": 0.0,
+      "nodes": [
+        {
+          "expiry": "2026-01-09",
+          "annuity": 1.0,
+          "time": 0.33881278538812787,
+          "vol": 0.008
+        },
+        {
+          "expiry": "2027-01-12",
+          "annuity": 1.0082191780821919,
+          "time": 1.3415525114155251,
+          "vol": 0.009
+        }
+      ],
+      "quotes": [
+        {
+          "maturity": "12M",
+          "vol": 0.008,
+          "pv": 0.0018577186797102136
+        },
+        {
+          "maturity": "24M",
+          "vol": 0.00866735652251448,
+          "pv": 0.006050588394704123
+        }
+      ]
+    }
+  ]
+}
diff --git a/crates/rhodium-pricing/tests/fixtures/steps.rs b/crates/rhodium-pricing/tests/fixtures/steps.rs
new file mode 100644
index 0000000000000000000000000000000000000000..461517f60a4d27f5a03f8010eebe070185e847a5
--- /dev/null
+++ b/crates/rhodium-pricing/tests/fixtures/steps.rs
@@ -0,0 +1,248 @@
+//! 定型テナーの OIS だけで会合日の階段を当てる、罰則付きの最小二乗の TONA 曲線（R4）。
+//!
+//! 曲線は Sum[Stitch{Steps(1 年までの 4 回の会合), Nodes(18M, 2Y), 1Y で接続}, Jumps(年末)]。
+//! 1M、6M、1Y の OIS が 5 区間の階段を、18M と 2Y の OIS が長期の節点を当てる。
+//! 階段は quote より未知数が多いので、隣り合う区間のフォワードの差の罰則で一意にする。
+//! 罰則を強くすると OIS の quote を再現しきれず、残差の残る過剰決定の当てはめになる。
+
+use rhodium_core::{Result, date::Date, quantity::YearFraction};
+use rhodium_market::{
+    curve::{CurveExpr, CurveId, CurveRole, CurveSpec, Extrapolation},
+    quote::{QuoteKind, QuoteSet},
+    structure::{BlockId, ParamId},
+};
+use rhodium_pricing::{
+    definition::{
+        CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
+        StepPenalty,
+    },
+    resolve::CollateralKey,
+};
+use serde::Deserialize;
+
+use super::ois::{Case, Node, Sample, cases, fixture_quotes, midnight};
+
+/// `QuantLib` で年末のターンを含めて校正した TONA の quote と DF。
+#[derive(Deserialize)]
+pub struct QuantLibTurn {
+    /// 年末の跳びが発効する日。
+    pub effective: Date,
+    /// 年末のスプレッド。
+    pub spread: f64,
+    /// 年末の受渡期間の日数係数。
+    pub fraction: f64,
+    /// OIS の quote とピラーの DF。
+    pub nodes: Vec<Node>,
+    /// ピラー以外の照合点の DF。
+    pub samples: Vec<Sample>,
+}
+
+/// 年末のターンの quote の ID。
+pub const TURN_QUOTE: &str = "JPY.TONA.TURN";
+/// OIS の当てはめの許容幅 σ（0.5bp）。
+pub const OIS_TOLERANCE: f64 = 0.5e-4;
+/// 罰則の重み λ。
+pub const PENALTY_WEIGHT: f64 = 1.0;
+/// 罰則の参照の幅 s（1bp）。隣の区間との 1bp の差が、0.5bp の残差の 1/4 の重さになる。
+pub const PENALTY_SCALE: f64 = 1e-4;
+/// 1 年までの会合の発効日。
+const MEETINGS: [&str; 4] = ["2025-03-19", "2025-06-17", "2025-09-19", "2025-12-19"];
+
+/// 年末のターンの `QuantLib` の fixture を読む。
+pub(crate) fn quantlib_turn() -> Result<QuantLibTurn> {
+    serde_json::from_str(include_str!("quantlib/ois_turn.json"))
+        .map_err(|error| rhodium_core::Error::InvalidData(error.to_string()))
+}
+
+/// 罰則付きの TONA の入力。`case.nodes` は `QuantLib` のターンの fixture の 5 テナー。
+pub struct PenalizedTona {
+    /// TONA の OIS の慣行と quote。
+    pub case: Case,
+    /// 年末のターン。
+    pub turn: QuantLibTurn,
+}
+
+impl PenalizedTona {
+    /// TONA の慣行と、ターンの fixture の quote を読む。
+    pub(crate) fn new() -> Result<Self> {
+        let [mut case, _] = cases()?;
+        let turn = quantlib_turn()?;
+        case.nodes.clone_from(&turn.nodes);
+        Ok(Self { case, turn })
+    }
+
+    /// 曲線の ID。
+    pub(crate) fn id(&self) -> Result<CurveId> {
+        self.case.index.parse()
+    }
+
+    /// 会合日の階段のブロック。
+    pub(crate) fn meetings(&self) -> Result<BlockId> {
+        BlockId::new(self.id()?, "meetings")
+    }
+
+    /// 長期の節点のブロック。
+    pub(crate) fn long(&self) -> Result<BlockId> {
+        BlockId::new(self.id()?, "long")
+    }
+
+    /// 年末の跳びのブロック。
+    pub(crate) fn turn_block(&self) -> Result<BlockId> {
+        BlockId::new(self.id()?, "turn")
+    }
+
+    /// 階段、長期の節点、年末の跳びからなる TONA 曲線。
+    pub(crate) fn curve(&self) -> Result<CurveSpec> {
+        let nodes = &self.turn.nodes;
+        CurveSpec::from_expr(
+            self.id()?,
+            self.case.currency,
+            [CurveRole::Discount, CurveRole::Projection],
+            CurveExpr::Sum(vec![
+                CurveExpr::Stitch {
+                    short: Box::new(CurveExpr::Steps {
+                        part: "meetings".into(),
+                        meetings: MEETINGS
+                            .iter()
+                            .map(|date| midnight(date.parse()?))
+                            .collect::<Result<_>>()?,
+                    }),
+                    long: Box::new(CurveExpr::Nodes {
+                        part: "long".into(),
+                        pillars: nodes[3..]
+                            .iter()
+                            .map(|node| midnight(node.date))
+                            .collect::<Result<_>>()?,
+                        extrapolation: Extrapolation::Reject,
+                    }),
+                    at: midnight(nodes[2].date)?,
+                },
+                CurveExpr::Jumps {
+                    part: "turn".into(),
+                    dates: vec![midnight(self.turn.effective)?],
+                },
+            ]),
+        )
+    }
+
+    /// 罰則付きの市場の定義。`gauss_newton` なら最小二乗の記録を Gauss–Newton で作る。
+    pub(crate) fn definition(&self, gauss_newton: bool) -> Result<MarketDefinition> {
+        let collateral = CollateralKey::Uncollateralized;
+        let mut bindings = self
+            .turn
+            .nodes
+            .iter()
+            .enumerate()
+            .map(|(position, node)| {
+                let block = if position < 3 {
+                    self.meetings()?
+                } else {
+                    self.long()?
+                };
+                self.case
+                    .binding(node, &collateral, ElementSelector::Block(block))?
+                    .with_tolerance(OIS_TOLERANCE)
+            })
+            .collect::<Result<Vec<_>>>()?;
+        // 年末のスプレッドは直接マークで、変換 TurnSpreadToJump が s·τ_turn の跳びにする。
+        bindings.push(InstrumentBinding::new(
+            TURN_QUOTE.parse()?,
+            QuoteKind::Spread,
+            CalibrationTarget::DirectMark,
+            QuoteTransform::TurnSpreadToJump(YearFraction::new(
+                self.turn.fraction,
+                self.case.basis(),
+            )?),
+            collateral,
+            ElementSelector::Param(ParamId::new(self.turn_block()?, 0)),
+        )?);
+        let definition = MarketDefinition::new(
+            self.case.currency.to_string().parse()?,
+            "Etc/UTC".parse()?,
+            [self.curve()?],
+            bindings,
+        )?
+        .with_step_penalties([StepPenalty::new(
+            self.meetings()?,
+            PENALTY_WEIGHT,
+            PENALTY_SCALE,
+        )?])?;
+        if gauss_newton {
+            definition.with_gauss_newton([self.meetings()?, self.long()?])
+        } else {
+            Ok(definition)
+        }
+    }
+
+    /// 会合 1 回（1M のピラー）で区切った 2 区間の階段を 1M と 6M の OIS が、1Y で接続した長期の
+    /// 節点を 18M と 2Y の OIS が、年末の跳びを直接マークが決める正方系の市場の定義（R3）。
+    pub(crate) fn square_definition(&self) -> Result<MarketDefinition> {
+        let nodes = &self.turn.nodes;
+        let curve = CurveSpec::from_expr(
+            self.id()?,
+            self.case.currency,
+            [CurveRole::Discount, CurveRole::Projection],
+            CurveExpr::Sum(vec![
+                CurveExpr::Stitch {
+                    short: Box::new(CurveExpr::Steps {
+                        part: "meetings".into(),
+                        meetings: vec![midnight(nodes[0].date)?],
+                    }),
+                    long: Box::new(CurveExpr::Nodes {
+                        part: "long".into(),
+                        pillars: nodes[2..]
+                            .iter()
+                            .map(|node| midnight(node.date))
+                            .collect::<Result<_>>()?,
+                        extrapolation: Extrapolation::Reject,
+                    }),
+                    at: midnight(nodes[1].date)?,
+                },
+                CurveExpr::Jumps {
+                    part: "turn".into(),
+                    dates: vec![midnight(self.turn.effective)?],
+                },
+            ]),
+        )?;
+        let collateral = CollateralKey::Uncollateralized;
+        let mut bindings = nodes
+            .iter()
+            .enumerate()
+            .map(|(position, node)| {
+                let block = if position < 2 {
+                    self.meetings()?
+                } else {
+                    self.long()?
+                };
+                self.case
+                    .binding(node, &collateral, ElementSelector::Block(block))
+            })
+            .collect::<Result<Vec<_>>>()?;
+        bindings.push(InstrumentBinding::new(
+            TURN_QUOTE.parse()?,
+            QuoteKind::Spread,
+            CalibrationTarget::DirectMark,
+            QuoteTransform::TurnSpreadToJump(YearFraction::new(
+                self.turn.fraction,
+                self.case.basis(),
+            )?),
+            collateral,
+            ElementSelector::Param(ParamId::new(self.turn_block()?, 0)),
+        )?);
+        MarketDefinition::new(
+            self.case.currency.to_string().parse()?,
+            "Etc/UTC".parse()?,
+            [curve],
+            bindings,
+        )
+    }
+
+    /// OIS とターンの quote。
+    pub(crate) fn quotes(&self) -> Result<QuoteSet> {
+        QuoteSet::new(
+            fixture_quotes(&self.case)?
+                .into_iter()
+                .chain([(TURN_QUOTE.parse()?, self.turn.spread)]),
+        )
+    }
+}
diff --git a/crates/rhodium-pricing/tests/fixtures/xccy.rs b/crates/rhodium-pricing/tests/fixtures/xccy.rs
new file mode 100644
index 0000000000000000000000000000000000000000..7178a552ab2db6f93ec1ecbf6db8db2c3eb5f382
--- /dev/null
+++ b/crates/rhodium-pricing/tests/fixtures/xccy.rs
@@ -0,0 +1,293 @@
+//! TONA と SOFR の OIS、USDJPY の FX ポイントと `MtM` XCCY で校正する、JPY 現金担保の USD 曲線（R6）。
+//!
+//! 外貨担保曲線 `USD-under-JPY` は `FxImplied` で、短期は FX スワップのポイント（ON、TN と
+//! 3M・6M の日付付きのポイント）、長期は TONA 対 SOFR の `MtM` XCCY のベーシス（1Y、18M）が決める。
+//! XCCY は JPY の TONA 脚にベーシスを付けて払い、USD の SOFR 脚の元本を USDJPY の fixing で
+//! 期首ごとにリセットする。JPY 現金担保の下で JPY は JPY-TONA、USD は外貨担保曲線で割り引く。
+
+use std::fmt::Write as _;
+
+use rhodium_core::{
+    Error, Result,
+    calendar::CalendarSpec,
+    currency::Currency,
+    fx::CcyPair,
+    id::IndexId,
+    quantity::Spread,
+    refdata::{
+        Convention, DailyPublication, FxSpotLag, IndexDefinition, IndexKind, IndexStatus,
+        PublicationRule, SwapConvention, SwapCouponConvention, SwapLegConvention,
+    },
+    tenor::TenorLength,
+    time::{Instant, LocalEventTime},
+};
+use rhodium_market::{
+    curve::{CurveId, Extrapolation},
+    fx::FxCollateral,
+    quote::QuoteSet,
+};
+use rhodium_pricing::{
+    definition::{
+        DirectFx, FxDatedQuote, FxDefinition, FxImpliedCurve, FxImpliedLongEnd, FxQuotes,
+        MarketDefinition,
+    },
+    resolve::CollateralKey,
+};
+use rhodium_products::{
+    agreement::{CashCollateralTerms, Remuneration},
+    standard::InstrumentTemplate,
+};
+
+use super::ois::{Case, Inputs, cases, definition, fixture_quotes, midnight};
+
+/// JPY 現金担保の USD 割引曲線（外貨担保曲線）。
+pub const IMPLIED: &str = "USD-under-JPY";
+/// `MtM` XCCY の慣行。
+pub const XCCY: &str = "USDJPY-XCCY";
+/// 元本のリセットに使う USDJPY の fixing。
+pub const FX_FIXING: &str = "USDJPY-WMR";
+/// スポットのアウトライト。
+pub const SPOT: f64 = 150.0;
+/// 短期と日付付きの FX ポイント（pip、0.01 円）と、その受渡日。
+pub const POINTS: [(&str, Option<&str>, f64); 4] = [
+    ("USDJPY.ON", None, -2.0),
+    ("USDJPY.TN", None, -2.5),
+    ("USDJPY.3M", Some("2025-04-09"), -95.0),
+    ("USDJPY.6M", Some("2025-07-09"), -190.0),
+];
+/// 長期の XCCY のテナー、quote（TONA 脚のベーシス）とピラー（最終の支払日）。
+pub const XCCY_QUOTES: [(&str, &str, f64, &str); 2] = [
+    ("USDJPY.XCCY.1Y", "1Y", -0.0035, "2026-01-13"),
+    ("USDJPY.XCCY.18M", "18M", -0.0040, "2026-07-13"),
+];
+
+/// TONA と SOFR の慣行と quote。
+pub struct XccyMarket {
+    /// JPY の OIS。
+    pub tona: Case,
+    /// USD の OIS。
+    pub sofr: Case,
+}
+
+/// JPY 現金担保（TONA 付利、スプレッド 0）。
+pub(crate) fn jpy_cash() -> Result<CollateralKey> {
+    Ok(CollateralKey::Cash(CashCollateralTerms {
+        ccy: Currency::JPY,
+        remuneration: Remuneration {
+            index: "JPY-TONA".parse()?,
+            spread: Spread::new(0.0)?,
+            floor: None,
+        },
+    }))
+}
+
+fn write_error(error: std::fmt::Error) -> Error {
+    Error::Internal(error.to_string())
+}
+
+impl XccyMarket {
+    /// `QuantLib` の fixture の TONA と SOFR を読む。
+    pub(crate) fn new() -> Result<Self> {
+        let [tona, sofr] = cases()?;
+        Ok(Self { tona, sofr })
+    }
+
+    /// 参照データ（FX fixing と XCCY の慣行を含む）と、JPY 現金担保の割引の行を持つ入力。
+    pub(crate) fn inputs(&self, extra_config: &str) -> Result<Inputs> {
+        self.inputs_with_spot_lag(extra_config, FxSpotLag::TwoBusinessDays)
+    }
+
+    /// [`Self::inputs`] の USDJPY の慣行のスポットラグを `spot_lag` にした入力。
+    pub(crate) fn inputs_with_spot_lag(
+        &self,
+        extra_config: &str,
+        spot_lag: FxSpotLag,
+    ) -> Result<Inputs> {
+        let mut config = String::new();
+        for (ccy, curve) in [("JPY", "JPY-TONA"), ("USD", IMPLIED)] {
+            writeln!(
+                config,
+                "[[discounting]]\ncurrency = \"{ccy}\"\ncollateral = {{ cash = {{ ccy = \"JPY\", index = \"JPY-TONA\", spread = \"0bp\" }} }}\ncurve = \"{curve}\""
+            )
+            .map_err(write_error)?;
+        }
+        config.push_str(extra_config);
+        Inputs::with_forecasts(
+            &[&self.tona, &self.sofr],
+            self.tona.as_of,
+            &config,
+            |parts, _| {
+                let index: IndexId = FX_FIXING.parse()?;
+                let calendar = CalendarSpec::Union(vec![
+                    CalendarSpec::Base("JPY".parse()?),
+                    CalendarSpec::Base("USD".parse()?),
+                ]);
+                let clock = LocalEventTime {
+                    time: "16:00:00".parse()?,
+                    zone: "Etc/UTC".parse()?,
+                };
+                parts.indices.insert(
+                    index.clone(),
+                    IndexDefinition {
+                        id: index.clone(),
+                        kind: IndexKind::FxFixing {
+                            pair: "USDJPY".parse()?,
+                        },
+                        calendar: calendar.clone(),
+                        publication: PublicationRule::Daily(DailyPublication {
+                            business_days: 0,
+                            calendar,
+                            time: clock.time,
+                            zone: clock.zone.clone(),
+                        }),
+                        observation_time: Some(clock),
+                        status: IndexStatus::Active,
+                    },
+                );
+                let leg = |case: &Case, fx_reset: Option<IndexId>| -> Result<SwapLegConvention> {
+                    let Some(Convention::Ois(ois)) = parts.conventions.get(&case.index.parse()?)
+                    else {
+                        return Err(Error::Internal("missing OIS test convention".into()));
+                    };
+                    let mut periods = ois.fixed.clone();
+                    periods.frequency = "3M".parse()?;
+                    Ok(SwapLegConvention {
+                        currency: case.currency,
+                        periods,
+                        payment: ois.payment.clone(),
+                        coupon: SwapCouponConvention::Compounded(ois.float.clone()),
+                        fx_reset,
+                    })
+                };
+                // 払い脚は JPY の TONA（ベーシスを付ける）、受け脚は USD の SOFR（元本をリセットする）。
+                let legs = [leg(&self.tona, None)?, leg(&self.sofr, Some(index))?];
+                let Some(Convention::Ois(ois)) = parts.conventions.get(&self.tona.index.parse()?)
+                else {
+                    return Err(Error::Internal("missing OIS test convention".into()));
+                };
+                let swap_lag = ois.spot_lag;
+                parts.conventions.insert(
+                    XCCY.parse()?,
+                    Convention::XccySwap(Box::new(SwapConvention {
+                        spot_lag: swap_lag,
+                        legs,
+                    })),
+                );
+                if let Some(convention) = parts.fx_conventions.get_mut(&Self::pair()?) {
+                    convention.spot_lag = spot_lag;
+                }
+                Ok(())
+            },
+        )
+    }
+
+    /// 外貨担保曲線の ID。
+    pub(crate) fn implied() -> Result<CurveId> {
+        IMPLIED.parse()
+    }
+
+    /// USDJPY の FX ポイントと長期の XCCY。`long_end` を差し替えると長期の行を変えられる。
+    pub(crate) fn fx(&self, long_end: Vec<FxImpliedLongEnd>) -> Result<FxDefinition> {
+        let [on, tn, three, six] = POINTS;
+        Ok(FxDefinition {
+            pair: "USDJPY".parse()?,
+            direct: Some(DirectFx {
+                market_collateral: "market".parse()?,
+                quotes: FxQuotes {
+                    spot: "USDJPY.SPOT".parse()?,
+                    on: Some(on.0.parse()?),
+                    tn: Some(tn.0.parse()?),
+                    sn: None,
+                    forwards: [three, six]
+                        .iter()
+                        .map(|(quote, date, _)| {
+                            Ok(FxDatedQuote {
+                                settlement: date
+                                    .ok_or_else(|| Error::Internal("dated FX point".into()))?
+                                    .parse()?,
+                                quote: quote.parse()?,
+                            })
+                        })
+                        .collect::<Result<_>>()?,
+                },
+                implied: Some(FxImpliedCurve {
+                    id: Self::implied()?,
+                    extrapolation: Extrapolation::default(),
+                    long_end,
+                }),
+            }),
+            collaterals: vec![FxCollateral {
+                id: "market".parse()?,
+                base: IMPLIED.parse()?,
+                quote: self.tona.index.parse()?,
+            }],
+            via: None,
+        })
+    }
+
+    /// 1Y と 18M の XCCY の長期の行。
+    pub(crate) fn long_end() -> Result<Vec<FxImpliedLongEnd>> {
+        XCCY_QUOTES
+            .iter()
+            .map(|(quote, tenor, _, pillar)| {
+                Ok(FxImpliedLongEnd {
+                    pillar: midnight(pillar.parse()?)?,
+                    quote: quote.parse()?,
+                    convention: XCCY.parse()?,
+                    tenor: tenor.parse()?,
+                    collateral: jpy_cash()?,
+                })
+            })
+            .collect()
+    }
+
+    /// TONA と SOFR の OIS 曲線と USDJPY の FX からなる市場の定義。
+    pub(crate) fn definition(&self) -> Result<MarketDefinition> {
+        let uncollateralized = CollateralKey::Uncollateralized;
+        definition(&[
+            (&self.tona, &uncollateralized),
+            (&self.sofr, &uncollateralized),
+        ])?
+        .with_fx([self.fx(Self::long_end()?)?])
+    }
+
+    /// OIS、FX と XCCY の quote。
+    pub(crate) fn quotes(&self) -> Result<QuoteSet> {
+        QuoteSet::new(
+            fixture_quotes(&self.tona)?
+                .into_iter()
+                .chain(fixture_quotes(&self.sofr)?)
+                .chain([("USDJPY.SPOT".parse()?, SPOT)])
+                .chain(
+                    POINTS
+                        .iter()
+                        .map(|(quote, _, value)| Ok::<_, Error>((quote.parse()?, *value)))
+                        .collect::<Result<Vec<_>>>()?,
+                )
+                .chain(
+                    XCCY_QUOTES
+                        .iter()
+                        .map(|(quote, _, value, _)| Ok::<_, Error>((quote.parse()?, *value)))
+                        .collect::<Result<Vec<_>>>()?,
+                ),
+        )
+    }
+
+    /// 長期の XCCY のテンプレート。
+    pub(crate) fn template(tenor: &str) -> Result<InstrumentTemplate> {
+        Ok(InstrumentTemplate::XccyBasis {
+            convention: XCCY.parse()?,
+            tenor: tenor.parse::<TenorLength>()?,
+        })
+    }
+
+    /// 受渡日の 16:00 UTC（ペアの慣行の受渡時刻）。
+    pub(crate) fn settlement(date: &str) -> Result<Instant> {
+        format!("{date}T16:00:00Z").parse()
+    }
+
+    /// 通貨ペア。
+    pub(crate) fn pair() -> Result<CcyPair> {
+        "USDJPY".parse()
+    }
+}
diff --git a/crates/rhodium-pricing/tests/forecast_measure.rs b/crates/rhodium-pricing/tests/forecast_measure.rs
new file mode 100644
index 0000000000000000000000000000000000000000..bba9eac681ec1aba2f9f9fa8ad7c9e8e71915c07
--- /dev/null
+++ b/crates/rhodium-pricing/tests/forecast_measure.rs
@@ -0,0 +1,319 @@
+//! 予想と割引の担保の違い（`CrossCollateralProjection`）を、予想曲線の予想の測度で判定する。
+//!
+//! 予想曲線は指標ごとに一本で、その予想は曲線を当てる校正商品の担保の下で定まる
+//! （設計書「担保による割引」）。同じ担保で割り引く取引は、予想曲線と割引曲線の ID が違っても
+//! 予想の担保の調整を要しない。担保の違う取引だけが決定論的なベーシスの仮定に頼る。
+
+use std::sync::Arc;
+
+use rhodium_core::{Error, Result, Warning, currency::Currency, money::Money, quantity::Spread};
+use rhodium_market::{
+    curve::{CurveExpr, CurveId, CurveRole, CurveSpec},
+    structure::MarketStructure,
+};
+use rhodium_pricing::{
+    calibration::{CalibratedMarket, GroupOutcome, calibrate},
+    context::Context,
+    prepare::{PrepareInput, Prepared, prepare},
+    resolve::CollateralKey,
+};
+use rhodium_products::{
+    agreement::{CashCollateralTerms, Remuneration},
+    standard::{Fra, Position, StandardProduct},
+};
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、担保の判定に使わない関数がある"
+)]
+#[path = "fixtures/ois.rs"]
+pub mod ois;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、担保の判定に使わない関数がある"
+)]
+#[path = "fixtures/multicurve.rs"]
+pub mod multicurve;
+
+use multicurve::{MulticurveCase, multicurve_cases};
+use ois::{Inputs, prepare_ois};
+
+/// 校正商品と同じ担保（無担保の行で TONA に割り引く）。
+const CALIBRATION: CollateralKey = CollateralKey::Uncollateralized;
+
+/// 25bp の付利スプレッドの JPY 現金担保で割り引く曲線（TONA + 25bp）。
+const SPREAD_CURVE: &str = "JPY-TONA-25";
+
+fn cash(ccy: Currency, index: &str, spread: f64) -> Result<CollateralKey> {
+    Ok(CollateralKey::Cash(CashCollateralTerms {
+        ccy,
+        remuneration: Remuneration {
+            index: index.parse()?,
+            spread: Spread::new(spread)?,
+            floor: None,
+        },
+    }))
+}
+
+/// 担保条件の違う割引曲線 `L(t) = L_base(t) − s·t`。校正しないので θ を持たない。
+fn spread_curve(id: &str, ccy: Currency, base: &str, spread: f64) -> Result<CurveSpec> {
+    CurveSpec::from_expr(
+        id.parse()?,
+        ccy,
+        [CurveRole::Discount],
+        CurveExpr::Sum(vec![CurveExpr::Ref(base.parse()?), CurveExpr::Flat(spread)]),
+    )
+}
+
+fn calibrated(
+    case: &MulticurveCase,
+    extra: Vec<CurveSpec>,
+    ctx: &Context<'_>,
+) -> Result<CalibratedMarket> {
+    let result = calibrate(&case.definition_with(extra)?, &case.quotes()?, ctx)?;
+    for group in result.groups() {
+        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
+            return Err(Error::Calibration(format!("{group:?}")));
+        }
+    }
+    Ok(result)
+}
+
+fn basis_trade(
+    case: &MulticurveCase,
+    collateral: &CollateralKey,
+    market: &CalibratedMarket,
+    ctx: &Context<'_>,
+) -> Result<Prepared> {
+    let node = &case.nodes[1];
+    let product = case.template(node)?.build(node.rate + 0.001, ctx.refdata)?;
+    prepare(
+        PrepareInput::Standard {
+            product: &product,
+            trade_date: case.overnight.as_of,
+        },
+        collateral,
+        &[],
+        market.market().structure(),
+        ctx,
+    )
+}
+
+fn has(warnings: &[Warning], expected: &Warning) -> bool {
+    warnings.contains(expected)
+}
+
+#[test]
+fn tibor_swap_under_calibration_collateral_needs_no_cross_collateral_convention() -> Result<()> {
+    // TIBOR 6M は TONA 担保のベーシスで校正するので、予想の測度は JPY-TONA。
+    // 同じ担保で割り引く TIBOR の取引と校正商品は、予想曲線と割引曲線の ID が違っても
+    // 担保の調整を要さず、警告を出さない。
+    let [jpy, _] = multicurve_cases()?;
+    let inputs = jpy.inputs()?;
+    let ctx = inputs.context();
+    let market = calibrated(&jpy, vec![], &ctx)?;
+    let term: CurveId = jpy.term_index.parse()?;
+    let overnight: CurveId = jpy.overnight.index.parse()?;
+    assert_eq!(
+        market.market().structure().forecast_measure(&term),
+        Some(&overnight)
+    );
+    assert!(!has(market.warnings(), &Warning::CrossCollateralProjection));
+    let trade = basis_trade(&jpy, &CALIBRATION, &market, &ctx)?;
+    assert!(!has(trade.warnings(), &Warning::CrossCollateralProjection));
+
+    // 近似を禁止しても、厳密に評価できる TIBOR のグループは解け、取引も準備できる。
+    let strict = multicurve_inputs(
+        &jpy,
+        "[conventions]\ndisabled = [\"cross_collateral_projection\"]\n",
+    )?;
+    let strict_ctx = strict.context();
+    let strict_market = calibrated(&jpy, vec![], &strict_ctx)?;
+    basis_trade(&jpy, &CALIBRATION, &strict_market, &strict_ctx)?;
+    Ok(())
+}
+
+#[test]
+fn tibor_swap_under_a_spread_csa_records_cross_collateral_projection() -> Result<()> {
+    // 25bp の付利の CSA は TONA + 25bp の曲線で割り引くので、TONA 担保で定まる TIBOR の予想と
+    // 担保が違う。この取引だけが CrossCollateralProjection を記録し、禁止すると準備できない。
+    let [jpy, _] = multicurve_cases()?;
+    let rule = format!(
+        "[[discounting]]\ncurrency = \"JPY\"\ncollateral = {{ cash = {{ ccy = \"JPY\", index = \"{}\", spread = \"25bp\" }} }}\ncurve = \"{SPREAD_CURVE}\"\n",
+        jpy.overnight.index
+    );
+    let extra = || spread_curve(SPREAD_CURVE, Currency::JPY, &jpy.overnight.index, 0.0025);
+    let inputs = multicurve_inputs(&jpy, &rule)?;
+    let ctx = inputs.context();
+    let market = calibrated(&jpy, vec![extra()?], &ctx)?;
+    let collateral = cash(Currency::JPY, &jpy.overnight.index, 0.0025)?;
+    let trade = basis_trade(&jpy, &collateral, &market, &ctx)?;
+    assert!(has(trade.warnings(), &Warning::CrossCollateralProjection));
+    // 校正商品の担保の取引は、同じ市場でも警告しない。
+    let native = basis_trade(&jpy, &CALIBRATION, &market, &ctx)?;
+    assert!(!has(native.warnings(), &Warning::CrossCollateralProjection));
+
+    let strict = multicurve_inputs(
+        &jpy,
+        &format!("{rule}[conventions]\ndisabled = [\"cross_collateral_projection\"]\n"),
+    )?;
+    let strict_ctx = strict.context();
+    let strict_market = calibrated(&jpy, vec![extra()?], &strict_ctx)?;
+    assert!(matches!(
+        basis_trade(&jpy, &collateral, &strict_market, &strict_ctx),
+        Err(Error::Unsupported(_))
+    ));
+    Ok(())
+}
+
+#[test]
+fn sofr_leg_under_jpy_collateral_records_cross_collateral_projection() -> Result<()> {
+    // SOFR の予想は USD 担保の OIS で定まる。JPY 現金担保の USD の支払いは別の USD 曲線で
+    // 割り引くので、SOFR の複利の予想と割引の担保が違う。
+    let [_, usd] = multicurve_cases()?;
+    let csa = "USD-JPYCSA";
+    let rule = format!(
+        "[[discounting]]\ncurrency = \"USD\"\ncollateral = {{ cash = {{ ccy = \"JPY\", index = \"JPY-TONA\", spread = \"0bp\" }} }}\ncurve = \"{csa}\"\n"
+    );
+    let inputs = multicurve_inputs(&usd, &rule)?;
+    let ctx = inputs.context();
+    let market = calibrated(
+        &usd,
+        vec![spread_curve(
+            csa,
+            Currency::USD,
+            &usd.overnight.index,
+            0.001,
+        )?],
+        &ctx,
+    )?;
+    let node = &usd.overnight.nodes[2];
+    let foreign = cash(Currency::JPY, "JPY-TONA", 0.0)?;
+    let trade = prepare_ois(
+        &usd.overnight,
+        node,
+        node.rate,
+        &foreign,
+        market.market(),
+        &ctx,
+    )?;
+    assert!(has(trade.warnings(), &Warning::CrossCollateralProjection));
+    let native = prepare_ois(
+        &usd.overnight,
+        node,
+        node.rate,
+        &CALIBRATION,
+        market.market(),
+        &ctx,
+    )?;
+    assert!(!has(native.warnings(), &Warning::CrossCollateralProjection));
+    Ok(())
+}
+
+#[test]
+fn fra_under_calibration_collateral_records_only_fra_discounting() -> Result<()> {
+    // FRA の決済は期首に P_d(s) で割り引くので、予想曲線と割引曲線が違えば同じ担保でも
+    // FraDiscounting の凸性を無視する。担保は同じなので CrossCollateralProjection は記録しない。
+    let [jpy, _] = multicurve_cases()?;
+    let inputs = jpy.inputs()?;
+    let ctx = inputs.context();
+    let market = calibrated(&jpy, vec![], &ctx)?;
+    let product = StandardProduct::Fra(Fra::new(
+        Money::new(Currency::JPY, 100_000_000.0),
+        jpy.term_index.parse()?,
+        jpy.nodes[1].date,
+        0.01,
+        Position::Long,
+    )?);
+    let trade = prepare(
+        PrepareInput::Standard {
+            product: &product,
+            trade_date: jpy.overnight.as_of,
+        },
+        &CALIBRATION,
+        &[],
+        market.market().structure(),
+        &ctx,
+    )?;
+    assert!(has(trade.warnings(), &Warning::FraDiscounting));
+    assert!(!has(trade.warnings(), &Warning::CrossCollateralProjection));
+    Ok(())
+}
+
+#[test]
+fn forecast_measure_is_unknown_when_bindings_disagree() -> Result<()> {
+    // TIBOR のベーシスの一つを 25bp の CSA で校正すると、予想を定める担保が一つに決まらない。
+    // 測度を記録せず、予想曲線と割引曲線が同じかどうかの比較に戻るので、TONA で割り引く
+    // 取引も CrossCollateralProjection を記録する。
+    let [jpy, _] = multicurve_cases()?;
+    let rule = format!(
+        "[[discounting]]\ncurrency = \"JPY\"\ncollateral = {{ cash = {{ ccy = \"JPY\", index = \"{}\", spread = \"25bp\" }} }}\ncurve = \"{SPREAD_CURVE}\"\n",
+        jpy.overnight.index
+    );
+    let inputs = multicurve_inputs(&jpy, &rule)?;
+    let ctx = inputs.context();
+    let collateral = cash(Currency::JPY, &jpy.overnight.index, 0.0025)?;
+    let definition = jpy.definition_with_collateral(
+        [spread_curve(
+            SPREAD_CURVE,
+            Currency::JPY,
+            &jpy.overnight.index,
+            0.0025,
+        )?],
+        |index| (index == 0).then(|| collateral.clone()),
+    )?;
+    let market = calibrate(&definition, &jpy.quotes()?, &ctx)?;
+    let term: CurveId = jpy.term_index.parse()?;
+    assert_eq!(market.market().structure().forecast_measure(&term), None);
+    let trade = basis_trade(&jpy, &CALIBRATION, &market, &ctx)?;
+    assert!(has(trade.warnings(), &Warning::CrossCollateralProjection));
+    Ok(())
+}
+
+#[test]
+fn market_structure_rejects_inconsistent_forecast_measures() -> Result<()> {
+    // 測度は同じ通貨の割引曲線でなければならず、予想の役割を持たない曲線には付けられない。
+    let [jpy, usd] = multicurve_cases()?;
+    let as_of = jpy.overnight.as_of;
+    let structure = || {
+        MarketStructure::new(
+            ois::midnight(as_of)?,
+            [
+                jpy.overnight.curve()?,
+                usd.overnight.curve()?,
+                spread_curve(SPREAD_CURVE, Currency::JPY, &jpy.overnight.index, 0.0025)?,
+            ],
+        )
+    };
+    let tona: CurveId = jpy.overnight.index.parse()?;
+    let sofr: CurveId = usd.overnight.index.parse()?;
+    let spread: CurveId = SPREAD_CURVE.parse()?;
+    for measures in [
+        vec![(tona.clone(), sofr.clone())],
+        vec![(spread.clone(), tona.clone())],
+        vec![("UNKNOWN".parse()?, tona.clone())],
+    ] {
+        assert!(matches!(
+            structure()?.with_forecast_measures(measures),
+            Err(Error::InvalidInput(_))
+        ));
+    }
+    let recorded = Arc::new(structure()?.with_forecast_measures([(tona.clone(), spread.clone())])?);
+    assert_eq!(recorded.forecast_measure(&tona), Some(&spread));
+    // 測度を付けた予想は測度の曲線と同じ担保で、予想曲線自身とは違う担保になる。
+    let projection = recorded.resolve_projection(&tona)?;
+    assert!(projection.is_forecast_under(&recorded.resolve_discount(&spread)?));
+    assert!(!projection.is_forecast_under(&recorded.resolve_discount(&tona)?));
+    // 記録の無い構造では、同じ曲線だけを同じ担保とみなす。
+    let plain = Arc::new(structure()?);
+    let projection = plain.resolve_projection(&tona)?;
+    assert!(projection.is_forecast_under(&plain.resolve_discount(&tona)?));
+    assert!(!projection.is_forecast_under(&plain.resolve_discount(&spread)?));
+    Ok(())
+}
+
+/// multicurve の入力に、割引の行や近似の設定を加える。
+fn multicurve_inputs(case: &MulticurveCase, extra: &str) -> Result<Inputs> {
+    case.inputs_with(extra)
+}
diff --git a/crates/rhodium-pricing/tests/fx_curves.rs b/crates/rhodium-pricing/tests/fx_curves.rs
index 73a3883ae0ebfcebede3e1f39b867dfed07f8855..f25bf0631703bd4a511b26661b46b7c9312f132f 100644
--- a/crates/rhodium-pricing/tests/fx_curves.rs
+++ b/crates/rhodium-pricing/tests/fx_curves.rs
@@ -9,7 +9,7 @@ use rhodium_core::{
     date::Weekday,
     fixing::{FixingStore, FixingStoreParts},
     fx::CcyPair,
-    refdata::{CurrencyDefinition, Deliverability, RefData, RefDataParts},
+    refdata::{CurrencyDefinition, Deliverability, FxConvention, FxSpotLag, RefData, RefDataParts},
     rounding::{Rounding, RoundingMode},
     time::{Instant, LocalEventTime, UtcOffset, ZoneRules, model_time},
 };
@@ -23,8 +23,7 @@ use rhodium_pricing::{
     config::PricingConfig,
     context::Context,
     definition::{
-        DirectFx, FxDatedQuote, FxDefinition, FxImpliedCurve, FxQuotes, FxSettlement, FxSpotLag,
-        MarketDefinition,
+        DirectFx, FxDatedQuote, FxDefinition, FxImpliedCurve, FxQuotes, MarketDefinition,
     },
 };
 use rhodium_testkit::assert_close;
@@ -54,8 +53,19 @@ struct Inputs {
     config: PricingConfig,
 }
 
+/// 既定のペアの慣行（USDJPY は T+2、USDCAD と CADJPY は T+1）。
+const CONVENTIONS: [(&str, FxSpotLag, u8); 3] = [
+    ("USDJPY", FxSpotLag::TwoBusinessDays, 2),
+    ("USDCAD", FxSpotLag::OneBusinessDay, 4),
+    ("CADJPY", FxSpotLag::OneBusinessDay, 2),
+];
+
 impl Inputs {
     fn new() -> Result<Self> {
+        Self::with_conventions(&CONVENTIONS)
+    }
+
+    fn with_conventions(conventions: &[(&str, FxSpotLag, u8)]) -> Result<Self> {
         let clock = clock()?;
         let mut parts = RefDataParts::default();
         parts.zones.insert(
@@ -95,6 +105,19 @@ impl Inputs {
                 },
             );
         }
+        // スポットラグ、同日受渡の締めと pip の桁は参照データのペアの慣行が決める。
+        for &(pair, spot_lag, pip_exponent) in conventions {
+            let pair: CcyPair = pair.parse()?;
+            parts.fx_conventions.insert(
+                pair,
+                FxConvention {
+                    pair,
+                    spot_lag,
+                    settlement: clock.clone(),
+                    pip_exponent,
+                },
+            );
+        }
         let refdata = RefData::new(parts)?;
         let fixings = FixingStore::new(&refdata, FixingStoreParts::default())?;
         let config = toml::from_str("version = 1\nreporting_currency = 'USD'")
@@ -137,11 +160,6 @@ fn pairs() -> Result<[FxDefinition; 2]> {
     let usd_jpy = FxDefinition {
         pair: "USDJPY".parse()?,
         direct: Some(DirectFx {
-            settlement: FxSettlement {
-                spot_lag: FxSpotLag::TwoBusinessDays,
-                settlement: clock()?,
-                pip_exponent: 2,
-            },
             market_collateral: "market".parse()?,
             quotes: FxQuotes {
                 spot: "USDJPY.SPOT".parse()?,
@@ -156,6 +174,7 @@ fn pairs() -> Result<[FxDefinition; 2]> {
             implied: Some(FxImpliedCurve {
                 id: "USD-under-JPY".parse()?,
                 extrapolation: Extrapolation::default(),
+                long_end: vec![],
             }),
         }),
         collaterals: vec![
@@ -167,11 +186,6 @@ fn pairs() -> Result<[FxDefinition; 2]> {
     let usd_cad = FxDefinition {
         pair: "USDCAD".parse()?,
         direct: Some(DirectFx {
-            settlement: FxSettlement {
-                spot_lag: FxSpotLag::OneBusinessDay,
-                settlement: clock()?,
-                pip_exponent: 4,
-            },
             market_collateral: "market".parse()?,
             quotes: FxQuotes {
                 spot: "USDCAD.SPOT".parse()?,
@@ -183,6 +197,7 @@ fn pairs() -> Result<[FxDefinition; 2]> {
             implied: Some(FxImpliedCurve {
                 id: "CAD-under-USD".parse()?,
                 extrapolation: Extrapolation::default(),
+                long_end: vec![],
             }),
         }),
         collaterals: vec![collateral("market", "USD", "CAD-under-USD")?],
@@ -405,11 +420,6 @@ fn fx_cross_direct_vs_triangulated() -> Result<()> {
         Ok(triangulated)
     })?;
     cross.direct = Some(DirectFx {
-        settlement: FxSettlement {
-            spot_lag: FxSpotLag::OneBusinessDay,
-            settlement: clock()?,
-            pip_exponent: 2,
-        },
         market_collateral: "market".parse()?,
         quotes: FxQuotes {
             spot: "CADJPY.SPOT".parse()?,
@@ -504,3 +514,40 @@ fn fx_missing_active_short_quote_only_blocks_its_pair() -> Result<()> {
     })?;
     assert_solved(&inputs.calibrate(&def, &quotes, AFTER)?)
 }
+
+#[test]
+fn fx_x0_chain_follows_refdata_spot_lag() -> Result<()> {
+    // スポットまでの短期の連鎖は参照データのペアの慣行で決まる。USDCAD に TN を足すと、
+    // T+1 の慣行ではスポットの後の区間になり連鎖にならないので、校正の入口で InvalidData になる。
+    // 市場の定義を変えずに参照データの慣行だけを T+2 にすると、同じ定義で校正できる。
+    let [_, mut cad] = pairs()?;
+    let direct = cad
+        .direct
+        .as_mut()
+        .ok_or_else(|| Error::Internal("test needs direct FX".into()))?;
+    direct.quotes.tn = Some("USDCAD.TN".parse()?);
+    let definition = definition([cad])?;
+    let quotes = quotes(&[("USDCAD.TN", 5.0)])?;
+    assert!(matches!(
+        Inputs::new()?.calibrate(&definition, &quotes, BEFORE),
+        Err(Error::InvalidData(message)) if message.contains("chain to spot")
+    ));
+    let t_plus_two = Inputs::with_conventions(&[("USDCAD", FxSpotLag::TwoBusinessDays, 4)])?;
+    let market = t_plus_two.calibrate(&definition, &quotes, BEFORE)?;
+    assert_solved(&market)?;
+    Ok(())
+}
+
+#[test]
+fn fx_definition_rejects_the_inverse_market_direction() -> Result<()> {
+    // pip とポイントは慣行の標準の向き（USDJPY）で読むので、逆向きの JPYUSD で書いた市場の
+    // quote は、ペアの慣行があっても校正の入口で拒否する。
+    let [mut jpy, _] = pairs()?;
+    jpy.pair = "JPYUSD".parse()?;
+    let definition = definition([jpy])?;
+    assert!(matches!(
+        Inputs::new()?.calibrate(&definition, &quotes(&[])?, BEFORE),
+        Err(Error::InvalidData(message)) if message.contains("convention direction USDJPY")
+    ));
+    Ok(())
+}
diff --git a/crates/rhodium-pricing/tests/least_squares_calibration.rs b/crates/rhodium-pricing/tests/least_squares_calibration.rs
new file mode 100644
index 0000000000000000000000000000000000000000..d51cc8c9883a36d5071f1adb4900e3e8598d3d6a
--- /dev/null
+++ b/crates/rhodium-pricing/tests/least_squares_calibration.rs
@@ -0,0 +1,485 @@
+//! 最小二乗の校正グループ: 手計算の検算、当てはめの許容幅の警告、warm start、定義の検査、
+//! 罰則付きの会合日の階段（R4）と Gauss–Newton の記録。
+
+use std::sync::Arc;
+
+use rhodium_core::{Error, Result, Warning, currency::Currency};
+use rhodium_market::{
+    curve::{CurveRole, CurveSpec, Extrapolation},
+    market::Market,
+    quote::{QuoteId, QuoteKind, QuoteSet},
+    structure::{BlockId, ParamId},
+};
+use rhodium_pricing::{
+    calibration::{CalibratedMarket, GroupOutcome, GroupSolution, HessianMode, calibrate},
+    context::Context,
+    definition::{
+        CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
+        StepPenalty,
+    },
+    prepare::evaluate,
+    resolve::CollateralKey,
+};
+use rhodium_testkit::assert_close;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、最小二乗の検査に使わない関数がある"
+)]
+#[path = "fixtures/ois.rs"]
+pub mod ois;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、最小二乗の検査に使わない関数がある"
+)]
+#[path = "fixtures/steps.rs"]
+pub mod steps;
+
+use ois::{Case, Inputs, cases, midnight, prepare_ois};
+use steps::{OIS_TOLERANCE, PENALTY_SCALE, PENALTY_WEIGHT, PenalizedTona};
+
+const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
+const MARKS: [&str; 2] = ["MARK.1", "MARK.2"];
+
+/// 一つのピラーの曲線の ln DF を、許容幅 σ の二つの直接マークで当てる市場の定義。
+fn direct_marks(tona: &Case, sigmas: [f64; 2]) -> Result<MarketDefinition> {
+    direct_marks_of_kind(tona, sigmas, QuoteKind::Mark)
+}
+
+/// [`direct_marks`] の quote の種類を `kind` にした定義。
+fn direct_marks_of_kind(
+    tona: &Case,
+    sigmas: [f64; 2],
+    kind: QuoteKind,
+) -> Result<MarketDefinition> {
+    let id: rhodium_market::curve::CurveId = tona.index.parse()?;
+    let curve = CurveSpec::new(
+        tona.index.parse()?,
+        Currency::JPY,
+        [CurveRole::Discount, CurveRole::Projection],
+        [midnight(tona.nodes[2].date)?],
+        Extrapolation::Reject,
+    )?;
+    let bindings = MARKS
+        .iter()
+        .zip(sigmas)
+        .map(|(quote, sigma)| {
+            InstrumentBinding::new(
+                quote.parse()?,
+                kind,
+                CalibrationTarget::DirectMark,
+                QuoteTransform::Identity,
+                COLLATERAL,
+                ElementSelector::Param(ParamId::new(BlockId::nodes(id.clone()), 0)),
+            )?
+            .with_tolerance(sigma)
+        })
+        .collect::<Result<Vec<_>>>()?;
+    MarketDefinition::new("JPY".parse()?, "Etc/UTC".parse()?, [curve], bindings)
+}
+
+fn marks(values: [f64; 2]) -> Result<QuoteSet> {
+    QuoteSet::new(
+        MARKS
+            .iter()
+            .zip(values)
+            .map(|(quote, value)| Ok((quote.parse()?, value)))
+            .collect::<Result<Vec<_>>>()?,
+    )
+}
+
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#196 で各 variant を網羅する"
+)]
+fn solution(market: &CalibratedMarket, index: usize) -> Result<&GroupSolution> {
+    match market.groups()[index].outcome() {
+        GroupOutcome::Solved(solution) => Ok(solution),
+        other => Err(Error::Calibration(format!("{other:?}"))),
+    }
+}
+
+fn linear(error: impl std::fmt::Display) -> Error {
+    Error::Numerical(error.to_string())
+}
+
+#[test]
+fn weighted_direct_marks_solve_the_hand_least_squares() -> Result<()> {
+    // 品質の手計算の例: min (θ − q₁)² + 4(θ − q₂)²（σ₁ = 1、σ₂ = 1/2）の解は θ = (q₁ + 4q₂)/5。
+    // 停留条件 G = (θ − q₁) + 4(θ − q₂) の微分は ∂G/∂θ = 5、∂G/∂q = (−1, −4) なので、
+    // quote への勾配は dθ/dq = −(∂G/∂θ)⁻¹·∂G/∂q = (1/5, 4/5)。
+    let [tona, _] = cases()?;
+    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
+    let ctx = inputs.context();
+    let [q1, q2] = [-0.01, -0.02];
+    let result = calibrate(&direct_marks(&tona, [1.0, 0.5])?, &marks([q1, q2])?, &ctx)?;
+    assert_eq!(result.groups().len(), 1);
+    let solution = solution(&result, 0)?;
+    assert_eq!(solution.hessian(), Some(HessianMode::Exact));
+    let theta = result.market().theta()[0];
+    assert_close!(theta, (q1 + 4.0 * q2) / 5.0, 1e-12, 1e-14);
+    let inverse = solution
+        .dg_dtheta()
+        .solve(&[1.0][..])
+        .map_err(linear)?
+        .as_slice()[0];
+    assert_close!(inverse, 0.2, 1e-10, 1e-12);
+    // 列は quote の ID の昇順（MARK.1、MARK.2）。
+    let dg_dq = [0, 1].map(|column| solution.dg_dq().get(0, column).map_err(linear));
+    let [dg_dq1, dg_dq2] = [dg_dq[0].clone()?, dg_dq[1].clone()?];
+    assert_close!(dg_dq1, -1.0, 1e-10, 1e-12);
+    assert_close!(dg_dq2, -4.0, 1e-10, 1e-12);
+    assert_close!(-inverse * dg_dq1, 0.2, 1e-10, 1e-12);
+    assert_close!(-inverse * dg_dq2, 0.8, 1e-10, 1e-12);
+    // 残差は重みを掛けない θ − q_i で、停止の誤差は停留点ではほぼ 0。
+    assert_close!(solution.residuals()[0], theta - q1, 1e-12, 1e-14);
+    assert_close!(solution.residuals()[1], theta - q2, 1e-12, 1e-14);
+    assert!(solution.theta_error()[0].abs() <= 1e-15);
+    Ok(())
+}
+
+#[test]
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#196 で各 variant を網羅する"
+)]
+fn fit_outside_tolerance_warns_per_quote() -> Result<()> {
+    // σ₁ = 1bp、σ₂ = 10bp で q₁ = −1%、q₂ = −2% を当てると、θ = (w₁²q₁ + w₂²q₂)/(w₁² + w₂²)
+    // = −1.0099% で、残差は F₁ = −0.99bp（σ₁ の内）、F₂ = 99bp（σ₂ の外）。外れた quote だけを警告する。
+    // 整合した quote に更新すると、再利用したグループの警告も解に合わせて消える。
+    let [tona, _] = cases()?;
+    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
+    let ctx = inputs.context();
+    let definition = direct_marks(&tona, [1e-4, 1e-3])?;
+    let result = calibrate(&definition, &marks([-0.01, -0.02])?, &ctx)?;
+    let outside: Vec<_> = result
+        .warnings()
+        .iter()
+        .filter_map(|warning| match warning {
+            Warning::CalibrationFitOutsideTolerance { quote } => Some(quote.as_str()),
+            _ => None,
+        })
+        .collect();
+    assert_eq!(outside, ["MARK.2"]);
+    let consistent = result.update(&marks([-0.01, -0.01])?)?;
+    assert!(
+        !consistent
+            .warnings()
+            .iter()
+            .any(|warning| matches!(warning, Warning::CalibrationFitOutsideTolerance { .. }))
+    );
+    Ok(())
+}
+
+#[test]
+fn least_squares_direct_marks_check_quotes_before_fitting() -> Result<()> {
+    // 直接マークの入力検査（FX スポットの種類なら正）は、正方系でも最小二乗でも解く前に行う。
+    // 許容幅を付けた二つのマーク −1 と 150 を平均すると θ は正になり、非正の quote が検査を
+    // 通り抜けて曲線を黙って動かす。検査は quote の種類だけで決まるので、曲線の要素に当てて確かめる。
+    let [tona, _] = cases()?;
+    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
+    let ctx = inputs.context();
+    let definition = direct_marks_of_kind(&tona, [1.0, 1.0], QuoteKind::FxSpot)?;
+    let fitted = calibrate(&definition, &marks([140.0, 150.0])?, &ctx)?;
+    assert_close!(fitted.market().theta()[0], 145.0, 1e-12, 1e-12);
+    let rejected = calibrate(&definition, &marks([-1.0, 150.0])?, &ctx)?;
+    assert!(
+        matches!(
+            rejected.groups()[0].outcome(),
+            GroupOutcome::Failed(Error::InvalidInput(message)) if message.contains("MARK.1")
+        ),
+        "{:?}",
+        rejected.groups()[0].outcome()
+    );
+    Ok(())
+}
+
+#[test]
+fn least_squares_update_warm_starts_a_strongly_convex_group() -> Result<()> {
+    // 直接マークだけの最小二乗は θ の強凸な二次式で解が一意なので、前の θ から解き直す。
+    // 更新の結果は、新しい quote から固定の初期値で解いた全体の校正と一致する。
+    let [tona, _] = cases()?;
+    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
+    let ctx = inputs.context();
+    let definition = direct_marks(&tona, [1.0, 0.5])?;
+    let base = calibrate(&definition, &marks([-0.01, -0.02])?, &ctx)?;
+    let quotes = marks([-0.01, -0.0195])?;
+    let updated = base.update(&quotes)?;
+    let full = calibrate(&definition, &quotes, &ctx)?;
+    solution(&updated, 0)?;
+    assert_close!(
+        updated.market().theta()[0],
+        full.market().theta()[0],
+        1e-12,
+        1e-14
+    );
+    assert_close!(
+        updated.market().theta()[0],
+        (-0.01 + 4.0 * -0.0195) / 5.0,
+        1e-12,
+        1e-14
+    );
+    Ok(())
+}
+
+#[test]
+fn least_squares_definitions_reject_inconsistent_options() -> Result<()> {
+    // 許容幅・罰則・Gauss–Newton の指定は、定義を作る境界と校正の入口で検査する。
+    let [tona, _] = cases()?;
+    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
+    let ctx = inputs.context();
+    let binding = || {
+        tona.binding(
+            &tona.nodes[0],
+            &COLLATERAL,
+            ElementSelector::Block(tona.block()?),
+        )
+    };
+    for sigma in [0.0, -1e-4, f64::NAN, f64::INFINITY] {
+        assert!(matches!(
+            binding()?.with_tolerance(sigma),
+            Err(Error::InvalidInput(_))
+        ));
+    }
+    for (weight, scale) in [(0.0, 1e-4), (1.0, f64::INFINITY), (-1.0, 1e-4)] {
+        assert!(matches!(
+            StepPenalty::new(tona.block()?, weight, scale),
+            Err(Error::InvalidInput(_))
+        ));
+    }
+    let square = ois::definition(&[(&tona, &COLLATERAL)])?;
+    // 罰則は Steps の部品にだけ掛けられる。
+    assert!(matches!(
+        square
+            .clone()
+            .with_step_penalties([StepPenalty::new(tona.block()?, 1.0, 1e-4)?]),
+        Err(Error::InvalidInput(message)) if message.contains("steps part")
+    ));
+    assert!(matches!(
+        square
+            .clone()
+            .with_gauss_newton([BlockId::new(tona.index.parse()?, "missing")?]),
+        Err(Error::InvalidInput(_))
+    ));
+    let quotes = QuoteSet::new(ois::fixture_quotes(&tona)?)?;
+    // 正方系（許容幅の無い binding）に Gauss–Newton を指定するのは定義の誤り。
+    let gauss_newton = square.clone().with_gauss_newton([tona.block()?])?;
+    assert!(matches!(
+        calibrate(&gauss_newton, &quotes, &ctx),
+        Err(Error::InvalidInput(message)) if message.contains("no binding tolerances")
+    ));
+    // 許容幅を持つ binding と持たない binding が同じグループにあると、形が一つに決まらない。
+    let block = tona.block()?;
+    let mixed = tona
+        .nodes
+        .iter()
+        .enumerate()
+        .map(|(index, node)| {
+            let element = u32::try_from(index).map_err(linear)?;
+            let binding = tona.binding(
+                node,
+                &COLLATERAL,
+                ElementSelector::Param(ParamId::new(block.clone(), element)),
+            )?;
+            if index == 0 {
+                binding.with_tolerance(1e-4)
+            } else {
+                Ok(binding)
+            }
+        })
+        .collect::<Result<Vec<_>>>()?;
+    let mixed = MarketDefinition::new("JPY".parse()?, "Etc/UTC".parse()?, [tona.curve()?], mixed)?;
+    assert!(matches!(
+        calibrate(&mixed, &quotes, &ctx),
+        Err(Error::InvalidInput(message)) if message.contains("mixes bindings")
+    ));
+    // 最小二乗でも残差と罰則の行が未知数より少なければ解が決まらない（5 要素を 4 本で当てる）。
+    let short = tona.nodes[..4]
+        .iter()
+        .map(|node| {
+            tona.binding(node, &COLLATERAL, ElementSelector::Block(block.clone()))?
+                .with_tolerance(1e-4)
+        })
+        .collect::<Result<Vec<_>>>()?;
+    let short = MarketDefinition::new("JPY".parse()?, "Etc/UTC".parse()?, [tona.curve()?], short)?;
+    assert!(matches!(
+        calibrate(&short, &quotes, &ctx),
+        Err(Error::InvalidInput(message)) if message.contains("least-squares")
+    ));
+    Ok(())
+}
+
+fn penalized(gauss_newton: bool) -> Result<(PenalizedTona, Inputs, CalibratedMarket)> {
+    let tona = PenalizedTona::new()?;
+    let inputs = Inputs::new(&[&tona.case], tona.case.as_of, "")?;
+    let result = calibrate(
+        &tona.definition(gauss_newton)?,
+        &tona.quotes()?,
+        &inputs.context(),
+    )?;
+    Ok((tona, inputs, result))
+}
+
+/// 曲線の θ を置き換えた市場で、OIS のパーレートを求める。
+fn par_rates(case: &Case, market: &Market<f64>, ctx: &Context<'_>) -> Result<Vec<f64>> {
+    case.nodes
+        .iter()
+        .map(|node| {
+            let value = |rate| {
+                let plan = prepare_ois(case, node, rate, &COLLATERAL, market, ctx)?;
+                Ok(evaluate(&plan, market)?.get(case.currency))
+            };
+            let zero = value(0.0)?;
+            let one: f64 = value(1.0)?;
+            Ok(-zero / (one - zero))
+        })
+        .collect()
+}
+
+#[test]
+fn penalized_meeting_steps_fit_standard_tenor_ois() -> Result<()> {
+    // R4: 1 年までの 4 回の会合（5 区間）を 1M、6M、1Y の OIS だけで当てる。未知数より quote が
+    // 少ないので、隣り合う区間のフォワードの差の罰則 (λ/2)Σ((f_{i+1} − f_i)/s)² で一意にする。
+    // 罰則を強くすると quote を再現しきれず、残差の残る当てはめになる。
+    let (tona, inputs, result) = penalized(false)?;
+    let ctx = inputs.context();
+    let blocks: Vec<_> = result
+        .groups()
+        .iter()
+        .map(|group| group.blocks().to_vec())
+        .collect();
+    // 年末の跳び（直接マーク）を上流に、階段、長期の節点の順に解く。
+    assert_eq!(
+        blocks,
+        [
+            vec![tona.turn_block()?],
+            vec![tona.meetings()?],
+            vec![tona.long()?]
+        ]
+    );
+    let [turn, meetings, long] = [0, 1, 2].map(|index| solution(&result, index));
+    let [turn, meetings, long] = [turn?, meetings?, long?];
+    assert_eq!(turn.hessian(), None);
+    assert_eq!(meetings.hessian(), Some(HessianMode::Exact));
+    assert_eq!(long.hessian(), Some(HessianMode::Exact));
+    assert!(
+        meetings
+            .residuals()
+            .iter()
+            .any(|residual| residual.abs() > 1e-6),
+        "{:?}",
+        meetings.residuals()
+    );
+
+    // 校正と取引の一致: OIS を取引として評価したパーは、quote に解の残差を足した値。
+    let market = result.market();
+    let pars = par_rates(&tona.case, market, &ctx)?;
+    let residuals: Vec<f64> = meetings
+        .residuals()
+        .iter()
+        .chain(long.residuals())
+        .copied()
+        .collect();
+    for ((par, node), residual) in pars.iter().zip(&tona.case.nodes).zip(&residuals) {
+        assert_close!(*par, node.rate + residual, 1e-10, 1e-12);
+    }
+
+    // 停留条件を独立に検算する。∂F_i/∂f_j を θ の中心差分で求め、
+    // G_j = Σ_i (F_i/σ²)·∂F_i/∂f_j + (λ/s²)·((f_j − f_{j−1}) − (f_{j+1} − f_j)) が 0 になる。
+    let range = market
+        .structure()
+        .layout()
+        .range(&tona.meetings()?)
+        .ok_or_else(|| Error::Internal("missing meetings block".into()))?;
+    let forwards = &market.theta()[range.clone()];
+    let h = 1e-7;
+    let penalty = PENALTY_WEIGHT / (PENALTY_SCALE * PENALTY_SCALE);
+    let last = forwards.len() - 1;
+    for (j, index) in range.enumerate() {
+        let bumped = |shift: f64| {
+            let mut theta = market.theta().to_vec();
+            theta[index] += shift;
+            par_rates(
+                &tona.case,
+                &Market::from_theta(Arc::clone(market.structure()), theta)?,
+                &ctx,
+            )
+        };
+        let [up, down] = [bumped(h)?, bumped(-h)?];
+        let data: Vec<f64> = (0..3)
+            .map(|i| residuals[i] / (OIS_TOLERANCE * OIS_TOLERANCE) * (up[i] - down[i]) / (2.0 * h))
+            .collect();
+        let left = if j > 0 {
+            forwards[j] - forwards[j - 1]
+        } else {
+            0.0
+        };
+        let right = if j < last {
+            forwards[j + 1] - forwards[j]
+        } else {
+            0.0
+        };
+        let gradient = data.iter().sum::<f64>() + penalty * (left - right);
+        let scale =
+            data.iter().map(|term| term.abs()).sum::<f64>() + penalty * (left.abs() + right.abs());
+        assert!(
+            gradient.abs() <= 1e-6 * scale,
+            "G[{j}] = {gradient:e} of {scale:e}"
+        );
+    }
+    Ok(())
+}
+
+#[test]
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#196 で各 variant を網羅する"
+)]
+fn gauss_newton_records_are_flagged_and_differ_from_exact() -> Result<()> {
+    // Gauss–Newton の記録は残差の 2 階の項を省くので、残差の残る階段のグループでは ∂G/∂θ が
+    // 厳密な記録と違う。解（θ）は同じで、近似の記録を解と警告の両方で報告する。
+    let (_, _, exact) = penalized(false)?;
+    let (tona, _, approximate) = penalized(true)?;
+    for (lhs, rhs) in exact
+        .market()
+        .theta()
+        .iter()
+        .zip(approximate.market().theta())
+    {
+        assert_close!(*lhs, *rhs, 1e-10, 1e-12);
+    }
+    let gauss_newton = solution(&approximate, 1)?;
+    assert_eq!(gauss_newton.hessian(), Some(HessianMode::GaussNewton));
+    let flagged: Vec<QuoteId> = approximate
+        .warnings()
+        .iter()
+        .filter_map(|warning| match warning {
+            Warning::CalibrationGaussNewton { quote } => Some(quote.clone()),
+            _ => None,
+        })
+        .collect();
+    let expected: Vec<QuoteId> = tona
+        .case
+        .nodes
+        .iter()
+        .map(|node| tona.case.quote(node))
+        .collect::<Result<_>>()?;
+    assert_eq!(flagged.len(), expected.len());
+    assert!(expected.iter().all(|quote| flagged.contains(quote)));
+    // 同じ右辺の逆作用素の値が、2 階の項の分だけ違う。
+    let rhs = vec![1.0; 5];
+    let [exact_inverse, approximate_inverse] = [solution(&exact, 1)?, gauss_newton].map(|s| {
+        s.dg_dtheta()
+            .solve(rhs.as_slice())
+            .map(|vector| vector.as_slice().to_vec())
+            .map_err(linear)
+    });
+    let difference = exact_inverse?
+        .iter()
+        .zip(approximate_inverse?)
+        .map(|(lhs, rhs)| (lhs - rhs).abs() / lhs.abs().max(1e-300))
+        .fold(0.0_f64, f64::max);
+    assert!(difference > 1e-6, "{difference:e}");
+    Ok(())
+}
diff --git a/crates/rhodium-pricing/tests/m1_linear.rs b/crates/rhodium-pricing/tests/m1_linear.rs
index ea7f098233f76039a3c33c6b4f3545f42c8e2b99..b9b675bc1888ec757e6907f4270152d6124a3d90 100644
--- a/crates/rhodium-pricing/tests/m1_linear.rs
+++ b/crates/rhodium-pricing/tests/m1_linear.rs
@@ -1,33 +1,28 @@
 //! M1 の標準取引を独立したキャッシュフローと閉形式に照合する。
 
-use std::{collections::BTreeMap, num::NonZeroU32, sync::Arc};
+use std::{num::NonZeroU32, sync::Arc};
 
 use rhodium_core::{
     Error, Result, Warning,
     calendar::{
-        BaseCalendar, BusinessDayConvention, CalendarSpec, EndOfMonthRule,
-        Extrapolation as CalendarExtrapolation, WeekendRule,
+        BaseCalendar, BusinessDayConvention, Extrapolation as CalendarExtrapolation, WeekendRule,
     },
-    currency::Currency,
+    currency::{Asset, Currency},
     date::{Date, DateInterval, Weekday},
     day_count::DayCount,
     fixing::{FixingStore, FixingStoreParts, ObservationPoint},
     fx::{CcyPair, FxRate},
-    id::{CalendarId, IndexId, SecurityId},
+    id::{CalendarId, IndexId},
     money::Money,
     quantity::{Rate, RateFuturesPrice, Simple, Spread},
     refdata::{
-        AssetValuation, Convention, CurrencyDefinition, DailyPublication, Deliverability,
-        FraConvention, FxForwardConvention, IndexDefinition, IndexKind, IndexStatus,
-        PublicationRule, RateFutureContract, RateFutureConvention, RateFutureRate, RefData,
-        RefDataParts, SecurityDefinition, SecurityKind, SecurityUnit,
+        Convention, FraConvention, FxConvention, FxForwardConvention, FxSpotLag, IndexKind,
+        PublicationRule, RateFutureRate, RefData, RefDataParts,
     },
-    rounding::{Rounding, RoundingMode},
     schedule::{
         ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination, ScheduleSpec,
     },
-    tenor::DateOffset,
-    time::{Instant, LocalEventTime, TzId, UtcOffset, ZoneRules, model_time},
+    time::{Instant, model_time},
 };
 use rhodium_market::{
     curve::{CurveRole, CurveSpec, Extrapolation},
@@ -37,204 +32,37 @@ use rhodium_market::{
     structure::{BlockId, MarketStructure},
 };
 use rhodium_pricing::{
-    calibration::calibrate,
+    calibration::{GroupOutcome, calibrate},
     config::{AnalyticPricer, Method, PricingConfig},
     context::Context,
     definition::{
-        CalibrationTarget, DirectFx, ElementSelector, FxDefinition, FxQuotes, FxSettlement,
-        FxSpotLag, InstrumentBinding, LgmDefinition, MarketDefinition, QuoteTransform,
+        CalibrationTarget, DirectFx, ElementSelector, FxDefinition, FxQuotes, InstrumentBinding,
+        MarketDefinition, QuoteTransform,
     },
     prepare::{PrepareInput, Prepared, evaluate, prepare},
     resolve::CollateralKey,
 };
-use rhodium_products::standard::{
-    Coupon, Fra, FxForward, FxSwap, InstrumentTemplate, Leg, Notional, NotionalExchange,
-    PayReceive, PaymentRounding, Position, RateFuture, StandardProduct, Swap,
+use rhodium_products::{
+    ir::{
+        Contract, Product,
+        analysis::{Linearity, Observation, linearity, observations},
+    },
+    standard::{
+        Coupon, Fra, FxForward, FxSwap, InstrumentTemplate, Leg, Notional, NotionalExchange,
+        PayReceive, PaymentRounding, Position, RateFuture, StandardProduct, Swap,
+    },
 };
 use rhodium_testkit::assert_close;
 use scioxide_ad::Tape;
 use scioxide_core::{ExecutionContext, InvalidInputReason, SciOxideError, argument};
 
-const START: &str = "2026-01-05";
-const AS_OF: &str = "2026-01-05T00:00:00Z";
-const FINAL: &str = "2026-01-07";
-const PILLAR: &str = "2026-05-08T16:00:00Z";
-const CLEARING: &str = "USD-FUT-JAN";
-const MULTIPLIER: f64 = 1_000.0;
-const RATE: f64 = 0.04;
+#[path = "fixtures/futures.rs"]
+mod futures;
 
-fn at(day: &str, time: &str) -> Result<Instant> {
-    format!("{day}T{time}:00Z").parse()
-}
-
-fn refdata(delay: DateOffset, settlement_time: &str) -> Result<RefData> {
-    RefData::new(reference_parts(delay, settlement_time)?)
-}
-
-#[expect(
-    clippy::too_many_lines,
-    reason = "暦、原指標、清算と限月のテスト用参照データをここにまとめる"
-)]
-fn reference_parts(delay: DateOffset, settlement_time: &str) -> Result<RefDataParts> {
-    let calendar: CalendarId = "USD".parse()?;
-    let zone: TzId = "Etc/UTC".parse()?;
-    let time = "16:00:00".parse()?;
-    let clock = LocalEventTime {
-        time,
-        zone: zone.clone(),
-    };
-    let mut parts = RefDataParts::default();
-    parts.zones.insert(
-        zone.clone(),
-        ZoneRules::new(
-            "2025-01-01T00:00:00Z".parse()?..="2027-01-01T00:00:00Z".parse()?,
-            UtcOffset::from_seconds(0)?,
-            vec![],
-        )?,
-    );
-    parts.calendars.insert(
-        calendar.clone(),
-        BaseCalendar::new(
-            "2025-01-01".parse()?..="2026-12-31".parse()?,
-            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
-            vec![],
-            CalendarExtrapolation::Strict,
-        )?,
-    );
-    parts.currencies.insert(
-        Currency::USD,
-        CurrencyDefinition {
-            code: Currency::USD,
-            minor_units: 2,
-            payment_rounding: Rounding::new(0.01, RoundingMode::HalfEven)?,
-            deliverability: Deliverability::Deliverable,
-            settlement_time: clock.clone(),
-            settlement_calendar: calendar.clone(),
-        },
-    );
-    let security: SecurityId = CLEARING.parse()?;
-    parts.securities.insert(
-        security.clone(),
-        SecurityDefinition {
-            id: security.clone(),
-            kind: SecurityKind::ListedFuture,
-            currency: Currency::USD,
-            unit: SecurityUnit::Contract,
-            issuer: None,
-            valuation: AssetValuation::Forward,
-        },
-    );
-    let calendar = CalendarSpec::Base(calendar);
-    for (name, kind) in [
-        (
-            "TERM",
-            IndexKind::TermRate {
-                ccy: Currency::USD,
-                day_count: DayCount::Act365F,
-                tenor: "1W".parse()?,
-                value_offset: "0D".parse()?,
-                bdc: BusinessDayConvention::Following,
-                end_of_month: EndOfMonthRule::Ignore,
-            },
-        ),
-        (
-            CLEARING,
-            IndexKind::RateFuturesSettlement {
-                contract: security,
-                underlying: "TERM".parse()?,
-            },
-        ),
-    ] {
-        let id: IndexId = name.parse()?;
-        parts.indices.insert(
-            id.clone(),
-            IndexDefinition {
-                id,
-                kind,
-                calendar: calendar.clone(),
-                status: IndexStatus::Active,
-                observation_time: Some(clock.clone()),
-                publication: PublicationRule::Daily(DailyPublication {
-                    business_days: 0,
-                    calendar: calendar.clone(),
-                    time,
-                    zone: zone.clone(),
-                }),
-            },
-        );
-    }
-    parts.conventions.insert(
-        "USD-RATE-FUTURES".parse()?,
-        Convention::RateFuture(Box::new(RateFutureConvention {
-            multiplier: MULTIPLIER,
-            settlement_lag: delay,
-            settlement_time: LocalEventTime {
-                time: settlement_time.parse()?,
-                zone,
-            },
-            margin_index: None,
-            contracts: BTreeMap::from([(
-                "2026-01".parse()?,
-                RateFutureContract {
-                    settlement_index: CLEARING.parse()?,
-                    final_clearing: FINAL.parse()?,
-                    rate: RateFutureRate::Term {
-                        fixing: FINAL.parse()?,
-                    },
-                },
-            )]),
-        })),
-    );
-    Ok(parts)
-}
-
-fn config(extra: &str) -> Result<PricingConfig> {
-    toml::from_str(&format!(
-        r#"
-version = 1
-reporting_currency = "USD"
-[[discounting]]
-currency = "USD"
-collateral = "uncollateralized"
-curve = "USD"
-[forecasting]
-TERM = "USD"
-{extra}
-"#
-    ))
-    .map_err(|error| Error::InvalidData(error.to_string()))
-}
-
-fn curve() -> Result<CurveSpec> {
-    CurveSpec::new(
-        "USD".parse()?,
-        Currency::USD,
-        [CurveRole::Discount, CurveRole::Projection],
-        [PILLAR.parse()?],
-        Extrapolation::Reject,
-    )
-}
-
-fn market(as_of: Instant, volatility: f64) -> Result<Market<f64>> {
-    let structure = Arc::new(MarketStructure::with_fx_and_lgm(
-        as_of,
-        [curve()?],
-        [],
-        [(Currency::USD, vec![])],
-    )?);
-    let curve = structure
-        .layout()
-        .range(&BlockId::nodes("USD".parse()?))
-        .ok_or_else(|| Error::Internal("missing curve block".into()))?;
-    let model = structure
-        .layout()
-        .range(&MarketStructure::lgm_block(Currency::USD)?)
-        .ok_or_else(|| Error::Internal("missing model block".into()))?;
-    let mut theta = vec![0.0; 3];
-    theta[curve.start] = -RATE * model_time(as_of, PILLAR.parse()?).years();
-    theta[model.start + 1] = volatility;
-    Market::from_theta(structure, theta)
-}
+use futures::{
+    AS_OF, CLEARING, FINAL, MULTIPLIER, PILLAR, RATE, START, at, config, curve, futures_definition,
+    futures_template, market, refdata, reference_parts,
+};
 
 fn future(position: Position, price: f64) -> Result<StandardProduct> {
     Ok(StandardProduct::RateFuture(RateFuture::new(
@@ -530,6 +358,99 @@ fn fra_settlement_matches_start_discounted_cashflow() -> Result<()> {
     Ok(())
 }
 
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
+        &CollateralKey::Uncollateralized,
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
+        1e-12_f64,
+        1e-8_f64
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
 #[test]
 fn futures_quote_calibrates_the_trade_to_par() -> Result<()> {
     // 先物のマークと独立した参照価格で校正し、同じ限月を約定した取引の PV がゼロであることを確かめる。
@@ -542,28 +463,8 @@ fn futures_quote_calibrates_the_trade_to_par() -> Result<()> {
         config: &config,
         fixings: &FixingStore::default(),
     };
-    let template = InstrumentTemplate::RateFuture {
-        convention: "USD-RATE-FUTURES".parse()?,
-        delivery_month: "2026-01".parse()?,
-    };
-    let binding = InstrumentBinding::new(
-        "USD.FUT.2026-01".parse()?,
-        QuoteKind::FuturesPrice,
-        CalibrationTarget::Price(template.clone()),
-        QuoteTransform::Identity,
-        CollateralKey::Uncollateralized,
-        ElementSelector::Block(BlockId::nodes("USD".parse()?)),
-    )?;
-    let definition =
-        MarketDefinition::new("USD".parse()?, "Etc/UTC".parse()?, [curve()?], [binding])?
-            .with_lgm([(
-                Currency::USD,
-                LgmDefinition {
-                    breaks: vec![],
-                    mean_reversion: vec!["USD.LGM.MR".parse()?],
-                    volatility: vec!["USD.LGM.VOL".parse()?],
-                },
-            )])?;
+    let template = futures_template()?;
+    let definition = futures_definition()?;
     let quote = oracle_price(
         as_of,
         &[
@@ -593,18 +494,30 @@ fn futures_quote_calibrates_the_trade_to_par() -> Result<()> {
 }
 
 fn fx_reference(pair: CcyPair) -> Result<RefData> {
+    RefData::new(fx_parts(pair)?)
+}
+
+fn fx_parts(pair: CcyPair) -> Result<RefDataParts> {
     let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
     let clock = parts.currencies[&Currency::USD].settlement_time.clone();
     parts.conventions.insert(
         "USDJPY-FWD".parse()?,
         Convention::FxForward(Box::new(FxForwardConvention {
             pair: "USDJPY".parse()?,
-            spot_lag: "2BD".parse()?,
-            calendar: CalendarSpec::Base("USD".parse()?),
             bdc: BusinessDayConvention::Following,
-            settlement: clock,
         })),
     );
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
     let mut currency = parts.currencies[&Currency::USD].clone();
     currency.code = Currency::JPY;
     currency.minor_units = 0;
@@ -618,9 +531,124 @@ fn fx_reference(pair: CcyPair) -> Result<RefData> {
     term.id = "FX-FIX".parse()?;
     term.kind = IndexKind::FxFixing { pair };
     parts.indices.insert(term.id.clone(), term);
+    Ok(parts)
+}
+
+/// USDJPY の慣行のスポットラグを `spot_lag` にし、JPY の決済カレンダーに 2026-01-13（火）の
+/// 祝日を置いた参照データ。ペアの暦は USD と JPY の休日の和になる。
+fn fx_reference_with_lag(spot_lag: FxSpotLag) -> Result<RefData> {
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
+    if let Some(jpy) = parts.currencies.get_mut(&Currency::JPY) {
+        jpy.settlement_calendar = tokyo;
+    }
+    if let Some(convention) = parts.fx_conventions.get_mut(&"USDJPY".parse()?) {
+        convention.spot_lag = spot_lag;
+    }
     RefData::new(parts)
 }
 
+#[test]
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#196 で各 variant を網羅する"
+)]
+fn mtm_xccy_reset_fixing_dates_follow_the_pair_spot_lag() -> Result<()> {
+    // 元本のリセットは期首がスポット受渡になる日に観測する。期首は 2026-01-07（水）と 01-14（水）で、
+    // ペアの暦は USD と JPY の休日の和（JPY は 01-13 が祝日）。
+    // T+1: 01-07 の 1 営業日前は 01-06、01-14 の 1 営業日前は 01-13 を飛ばして 01-12。
+    // T+2: 01-07 の 2 営業日前は 01-05、01-14 の 2 営業日前は 01-12、01-09。
+    let swap = mtm_swap()?;
+    for (spot_lag, expected) in [
+        (FxSpotLag::OneBusinessDay, ["2026-01-06", "2026-01-12"]),
+        (FxSpotLag::TwoBusinessDays, ["2026-01-05", "2026-01-09"]),
+    ] {
+        let reference = fx_reference_with_lag(spot_lag)?;
+        let product = swap.to_product(START.parse()?, &reference)?;
+        let fx: IndexId = "FX-FIX".parse()?;
+        let mut dates: Vec<_> = observations(&product)
+            .into_iter()
+            .filter_map(|observation| match observation {
+                Observation::Fixing {
+                    index,
+                    point: ObservationPoint::Date(date),
+                } if *index == fx => Some(date),
+                _ => None,
+            })
+            .collect();
+        dates.sort_unstable();
+        let expected = expected
+            .iter()
+            .map(|date| date.parse())
+            .collect::<Result<Vec<Date>>>()?;
+        assert_eq!(dates, expected, "{spot_lag:?}");
+    }
+    Ok(())
+}
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
+            &CollateralKey::Uncollateralized,
+            &[],
+            market.structure(),
+            &ctx,
+        )?;
+        // 受渡日と支払日が違うので、タイミング調整を省いた近似を記録する。
+        assert!(plan.warnings().contains(&Warning::PaymentDelayIgnored));
+        let pv = evaluate(&plan, &market)?.get(Currency::JPY);
+        assert_close!(
+            pv,
+            forward(at(value_date, "16:00")?) * discount(0.01, pay),
+            1e-12_f64,
+            1e-6_f64
+        );
+        values.push(pv);
+    }
+    assert!((values[0] - values[1]).abs() > 1e-6);
+    Ok(())
+}
+
 fn fx_market(as_of: Instant) -> Result<Market<f64>> {
     let jpy = CurveSpec::new(
         "JPY".parse()?,
@@ -950,6 +978,161 @@ fn mtm_xccy_fx_reset_matches_cashflows_in_both_directions() -> Result<()> {
     Ok(())
 }
 
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
+    if let Some(convention) = parts.fx_conventions.get_mut(&"USDJPY".parse()?) {
+        convention.settlement.time = "06:00:00".parse()?;
+    }
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
+        &CollateralKey::Uncollateralized,
+        &[],
+        market.structure(),
+        &ctx,
+    )?;
+    assert_close!(
+        evaluate(&plan, &market)?.get(Currency::JPY),
+        forward(value) * discount(0.01, pay),
+        1e-12_f64,
+        1e-6_f64
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
+        let plan = prepare(
+            PrepareInput::Standard {
+                product: &product,
+                trade_date: START.parse()?,
+            },
+            &CollateralKey::Uncollateralized,
+            &[],
+            market.structure(),
+            &ctx,
+        )?;
+        assert_eq!(plan.warnings(), warnings.as_slice(), "{reset:?}");
+        let n0 = 150_000.0 / forward(value);
+        let n1 = 150_000.0 / forward(at(second_value, "06:00")?);
+        assert_close!(
+            evaluate(&plan, &market)?.get(Currency::USD),
+            usd(n0, n1),
+            1e-12_f64,
+            1e-8_f64
+        );
+    }
+    Ok(())
+}
+
 #[test]
 fn overnight_futures_compound_or_average_the_reference_period() -> Result<()> {
     // 先物の清算指標を切り替えず、最終値の参照金利だけを複利・加重平均で組み立てる。
@@ -1000,6 +1183,91 @@ fn overnight_futures_compound_or_average_the_reference_period() -> Result<()> {
     Ok(())
 }
 
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
+        let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
+        if let Some(term) = parts.indices.get_mut(&"TERM".parse()?) {
+            term.kind = IndexKind::OvernightRate {
+                ccy: Currency::USD,
+                day_count: DayCount::Act365F,
+            };
+            // 翌日物は翌営業日の 16 時に公表する。A では 1 日目が未公表で、期間が A をまたぐ。
+            if let PublicationRule::Daily(rule) = &mut term.publication {
+                rule.business_days = 1;
+            }
+        }
+        if let Some(Convention::RateFuture(convention)) =
+            parts.conventions.get_mut(&"USD-RATE-FUTURES".parse()?)
+            && let Some(contract) = convention.contracts.values_mut().next()
+        {
+            contract.rate = reference_rate;
+        }
+        let reference = RefData::new(parts)?;
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
+            1e-12_f64,
+            1e-8_f64
+        );
+        assert!(
+            warnings
+                .iter()
+                .any(|warning| matches!(warning, Warning::UnpublishedFixingForecast { .. }))
+        );
+    }
+    Ok(())
+}
+
 fn check_gradient(
     plan: &Prepared,
     market: &Market<f64>,
@@ -1215,11 +1483,6 @@ fn fx_outright_quote_calibrates_a_par_trade_in_both_currencies() -> Result<()> {
     .with_fx([FxDefinition {
         pair,
         direct: Some(DirectFx {
-            settlement: FxSettlement {
-                spot_lag: FxSpotLag::TwoBusinessDays,
-                settlement: reference.currency(Currency::USD)?.settlement_time.clone(),
-                pip_exponent: 2,
-            },
             market_collateral: "csa".parse()?,
             quotes: FxQuotes {
                 spot: "USDJPY.SPOT".parse()?,
@@ -1262,3 +1525,644 @@ fn fx_outright_quote_calibrates_a_par_trade_in_both_currencies() -> Result<()> {
     );
     Ok(())
 }
+
+/// USD のピラーを直接マークで、JPY のピラーを 1W の USDJPY アウトライトで当てる市場の定義。
+/// `outrights` は JPY の曲線に当てるアウトライトの quote の ID と許容幅 σ（無ければ正方系）。
+fn outright_definition(outrights: &[(&str, Option<f64>)]) -> Result<MarketDefinition> {
+    let pair: CcyPair = "USDJPY".parse()?;
+    let template = InstrumentTemplate::FxForward {
+        convention: "USDJPY-FWD".parse()?,
+        tenor: "1W".parse()?,
+    };
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
+        ElementSelector::Param(rhodium_market::structure::ParamId::new(
+            BlockId::nodes("USD".parse()?),
+            0,
+        )),
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
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#196 で各 variant を網羅する"
+)]
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
+    let template = InstrumentTemplate::FxForward {
+        convention: "USDJPY-FWD".parse()?,
+        tenor: "1W".parse()?,
+    };
+    let trade = template.build_for_trade(fitted, START.parse()?, &reference)?;
+    let market = calibrated.market();
+    let fx = market.with_view(|view| view.fx_today("USDJPY".parse()?))?;
+    assert_close!(
+        price(&trade, market, &context)?
+            .convert(Currency::JPY, &[fx])?
+            .amount(),
+        0.0_f64,
+        1e-12_f64,
+        1e-6_f64
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
+    let outside: Vec<_> = calibrated
+        .warnings()
+        .iter()
+        .filter_map(|warning| match warning {
+            Warning::CalibrationFitOutsideTolerance { quote } => Some(quote.as_str()),
+            _ => None,
+        })
+        .collect();
+    assert_eq!(outside, ["USDJPY.FWD.1W.B"]);
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
+    let plan = prepare(
+        PrepareInput::Standard {
+            product: &swap,
+            trade_date: START.parse()?,
+        },
+        &CollateralKey::Uncollateralized,
+        &[],
+        market.structure(),
+        &ctx,
+    )?;
+    // 元本が固定なので、FX と金利の積の近似（MtmFxForwardProduct）は使わない。
+    assert!(!plan.warnings().contains(&Warning::MtmFxForwardProduct));
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
+    assert_close!(pv.get(Currency::JPY), jpy, 1e-12_f64, 1e-6_f64);
+    assert_close!(pv.get(Currency::USD), usd, 1e-12_f64, 1e-8_f64);
+    Ok(())
+}
+
+#[test]
+fn term_coupon_payment_delay_is_recorded() -> Result<()> {
+    // ターム金利のクーポンを期末から 2 日遅れて払うと、自然な支払日（期末）とのずれのタイミング調整を
+    // 省く近似（PaymentDelayIgnored）を記録し、予想は期末の DF 比、割引は実際の支払日で求める。
+    // 近似を禁止すると準備できない。
+    let as_of = AS_OF.parse()?;
+    let market = market(as_of, 0.0)?;
+    let reference = refdata("0D".parse()?, "16:00:00")?;
+    let fixings = FixingStore::default();
+    let schedule = ScheduleSpec::Explicit(ExplicitPeriods::new(
+        vec![ExplicitPeriod {
+            accrual: DateInterval::new("2026-01-07".parse()?, "2026-01-14".parse()?)?,
+            payment: "2026-01-16".parse()?,
+            reset: None,
+            settlement: None,
+            references: vec![],
+        }],
+        ExplicitDatePolicy::KeepAsConfirmed,
+        PeriodCombination::Sum,
+        None,
+    )?);
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
+        prepare(
+            PrepareInput::Standard {
+                product: &swap,
+                trade_date: START.parse()?,
+            },
+            &CollateralKey::Uncollateralized,
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
+    let plan = prepare_with(&config("")?)?;
+    assert!(plan.warnings().contains(&Warning::PaymentDelayIgnored));
+    let tau = 7.0 / 365.0;
+    let forward = (RATE * tau).exp_m1() / tau;
+    let df = (-RATE * model_time(as_of, at("2026-01-16", "16:00")?).years()).exp();
+    assert_close!(
+        evaluate(&plan, &market)?.get(Currency::USD),
+        100_000.0 * tau * forward * df,
+        1e-12_f64,
+        1e-8_f64
+    );
+    assert!(matches!(
+        prepare_with(&config(
+            "[conventions]\ndisabled = ['payment_delay_ignored']"
+        )?),
+        Err(Error::Unsupported(_))
+    ));
+    Ok(())
+}
+
+/// [`fx_market`] に、JPY を 25bp 上乗せした割引曲線 `JPY-25` を足した市場。
+fn fx_market_with_spread_curve(as_of: Instant, pair: CcyPair) -> Result<Market<f64>> {
+    let jpy = CurveSpec::new(
+        "JPY".parse()?,
+        Currency::JPY,
+        [CurveRole::Discount, CurveRole::Projection],
+        [PILLAR.parse()?],
+        Extrapolation::Reject,
+    )?;
+    let spread = CurveSpec::from_expr(
+        "JPY-25".parse()?,
+        Currency::JPY,
+        [CurveRole::Discount],
+        rhodium_market::curve::CurveExpr::Sum(vec![
+            rhodium_market::curve::CurveExpr::Ref("JPY".parse()?),
+            rhodium_market::curve::CurveExpr::Flat(0.0025),
+        ]),
+    )?;
+    let structure = Arc::new(MarketStructure::with_fx(
+        as_of,
+        [curve()?, jpy, spread],
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
+    let market = Market::from_theta(structure, theta)?;
+    Ok(market)
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
+    let pair: CcyPair = "USDJPY".parse()?;
+    let market = fx_market_with_spread_curve(as_of, pair)?;
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
+    let collateral = CollateralKey::Cash(rhodium_products::agreement::CashCollateralTerms {
+        ccy: Currency::JPY,
+        remuneration: rhodium_products::agreement::Remuneration {
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
+    let native = prepare_under(&fx_config("")?, &CollateralKey::Uncollateralized)?;
+    assert!(
+        !native
+            .warnings()
+            .contains(&Warning::CrossCollateralProjection)
+    );
+    let spread_config = fx_config(rules)?;
+    let foreign = prepare_under(&spread_config, &collateral)?;
+    assert!(
+        foreign
+            .warnings()
+            .contains(&Warning::CrossCollateralProjection)
+    );
+    // 予想は担保区分の曲線のフォワードのまま、割引だけが 25bp の曲線に変わる。
+    let ratio = evaluate(&foreign, &market)?.get(Currency::JPY)
+        / evaluate(&native, &market)?.get(Currency::JPY);
+    let pay = model_time(as_of, at("2026-01-14", "16:00")?).years();
+    assert_close!(ratio, (-0.0025 * pay).exp(), 1e-12_f64, 1e-14_f64);
+    let strict = fx_config(&format!(
+        "{rules}[conventions]\ndisabled = ['cross_collateral_projection']"
+    ))?;
+    assert!(matches!(
+        prepare_under(&strict, &collateral),
+        Err(Error::Unsupported(_))
+    ));
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
+    if let Some(Convention::RateFuture(convention)) =
+        parts.conventions.get_mut(&"USD-RATE-FUTURES".parse()?)
+    {
+        convention.margin_index = Some("MARGIN".parse()?);
+    }
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
+    assert!(!warnings("USD")?.contains(&Warning::CrossCollateralProjection));
+    assert!(warnings("USD-ALT")?.contains(&Warning::CrossCollateralProjection));
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
+        ElementSelector::Param(rhodium_market::structure::ParamId::new(
+            BlockId::nodes("USD-OIS".parse()?),
+            0,
+        )),
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
+        rhodium_pricing::definition::LgmDefinition {
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
+        if margin
+            && let Some(Convention::RateFuture(convention)) =
+                parts.conventions.get_mut(&"USD-RATE-FUTURES".parse()?)
+        {
+            convention.margin_index = Some("MARGIN".parse()?);
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
+            &CollateralKey::Uncollateralized,
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
