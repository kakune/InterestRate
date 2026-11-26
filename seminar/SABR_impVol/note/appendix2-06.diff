diff --git a/crates/rhodium-pricing/src/calibration.rs b/crates/rhodium-pricing/src/calibration.rs
index 78d97fb9480370dfc2c3dc69fbe15a7dba338e57..ee0fdaed86ae4c91f3d962465c6cf9518c482132 100644
--- a/crates/rhodium-pricing/src/calibration.rs
+++ b/crates/rhodium-pricing/src/calibration.rs
@@ -1,4 +1,8 @@
-//! OIS のパーレートから曲線を解き、校正の線形化を保存する。
+//! 市場の定義と quote から θ を解き、解での線形化を記録する。
+//!
+//! 曲線の部品、担保別の FX、金利ボラ面、LGM のマークを、依存の順に並べた校正グループごとに、
+//! 正方系か（許容幅を持つ binding のグループは）最小二乗で解く。
+//! 記録した線形化は、quote 感応度の随伴と、変わった quote の影響範囲だけの再校正に使う。
 
 mod groups;
 mod solver;
diff --git a/crates/rhodium-products/Cargo.toml b/crates/rhodium-products/Cargo.toml
index 252bd4b14438fae0b84b95761f3552ccdf1eca63..0d4f6530cf74cd9261605ab5e2aaabcf9468d0f9 100644
--- a/crates/rhodium-products/Cargo.toml
+++ b/crates/rhodium-products/Cargo.toml
@@ -1,6 +1,6 @@
 [package]
 name = "rhodium-products"
-description = "コア IR、DSL、標準商品と取引の定義"
+description = "コア IR、標準商品と校正商品のテンプレート、取引と担保条件の定義（DSL は M4 で足す）"
 version.workspace = true
 edition.workspace = true
 rust-version.workspace = true
diff --git a/crates/rhodium-products/src/lib.rs b/crates/rhodium-products/src/lib.rs
index 10bbe0244458702b829038e567b385f8d94e7d83..5351f5aff70b9170d3cf4053e123133c7c025de7 100644
--- a/crates/rhodium-products/src/lib.rs
+++ b/crates/rhodium-products/src/lib.rs
@@ -1,4 +1,8 @@
-//! コア IR、DSL、標準商品と取引の定義。
+//! コア IR、標準商品と校正商品のテンプレート、取引と担保条件の定義。
+//!
+//! IR は契約と式の builder、因果性の検査、評価時点での簡約と解析を持つ。
+//! 標準商品は M1 の金利と FX のリニア商品で、金利オプションは M2、ほかの資産クラスは各マイルストーンで足す。
+//! 独自の商品を書く DSL は M4、取引先とネッティングセットは M5 で足す。
 
 #![forbid(unsafe_code)]
 
diff --git a/crates/rhodium-products/src/standard/leg.rs b/crates/rhodium-products/src/standard/leg.rs
index dfb1829277937da050a3e411a89deff0dadf0537..787f8109498cb69b955b905ff99d0e98fab85cbf 100644
--- a/crates/rhodium-products/src/standard/leg.rs
+++ b/crates/rhodium-products/src/standard/leg.rs
@@ -1,4 +1,7 @@
-//! 固定元本の固定脚と RFR 複利脚を組み立てる。
+//! 固定元本か FX リセットの元本と、固定・ターム金利・RFR 複利のクーポンを持つ脚を IR に組み立てる。
+//!
+//! キャップレットの校正が読む変動クーポンの抽出も持つ。
+//! 校正の経路はキャップの標準商品とともに #42 で整理する。
 
 use rhodium_core::{
     Error, Result, Warning,
diff --git a/crates/rhodium-risk/Cargo.toml b/crates/rhodium-risk/Cargo.toml
index 5c36e7b504c6ea300160e12ca97b5c5c4061fee2..a8b430aaaffcec4b08c51406cf0071c4bd4e9839 100644
--- a/crates/rhodium-risk/Cargo.toml
+++ b/crates/rhodium-risk/Cargo.toml
@@ -1,6 +1,6 @@
 [package]
 name = "rhodium-risk"
-description = "価格の感応度、シナリオと規制計算"
+description = "市場パラメータの勾配と quote 感応度（シナリオと規制計算は M3 以降で足す）"
 version.workspace = true
 edition.workspace = true
 rust-version.workspace = true
diff --git a/crates/rhodium-risk/src/lib.rs b/crates/rhodium-risk/src/lib.rs
index 4c3501a22541a927f91cdd9ca9961257f7bbef4d..df63c8ff1e65d7652d58eeed99a77a9d6cd57d65 100644
--- a/crates/rhodium-risk/src/lib.rs
+++ b/crates/rhodium-risk/src/lib.rs
@@ -1,4 +1,7 @@
-//! 価格の感応度、シナリオと規制計算。
+//! 市場パラメータ θ に対する AAD の勾配と、校正の随伴による quote 感応度。
+//!
+//! quote デルタは受渡通貨別か報告通貨で集計する。ガンマ、シナリオ、theta と P&L 説明、
+//! ゴールシーク、ヘッジへの付け替え、SIMM、FRTB-SA は M3、SA-CCR は M6 で足す。
 
 #![forbid(unsafe_code)]
 
diff --git a/crates/rhodium-xva/Cargo.toml b/crates/rhodium-xva/Cargo.toml
index 2ddef30b7502f211360e8d497725b023ab06bba8..fe0a736ac7d2274e0ce8443272d08591d885160d 100644
--- a/crates/rhodium-xva/Cargo.toml
+++ b/crates/rhodium-xva/Cargo.toml
@@ -1,6 +1,6 @@
 [package]
 name = "rhodium-xva"
-description = "エクスポージャー、担保、XVA の指標と感応度"
+description = "エクスポージャー、担保、XVA の指標と感応度（M5 と M6 で足す。現在は空）"
 version.workspace = true
 edition.workspace = true
 rust-version.workspace = true
diff --git a/crates/rhodium-xva/src/lib.rs b/crates/rhodium-xva/src/lib.rs
index 5407cc5aac8e48e9e8810854381399b77bcfb23f..75e5548d3d93d455e4b16975bdc72e7c3773c518 100644
--- a/crates/rhodium-xva/src/lib.rs
+++ b/crates/rhodium-xva/src/lib.rs
@@ -1,4 +1,7 @@
-//! エクスポージャー、担保、XVA の指標と感応度。
+//! エクスポージャー、担保、XVA の指標と感応度を置く crate。
+//!
+//! 中身は M5（エクスポージャー、担保、CVA などの XVA、XVA の感応度）と M6（IM、MVA、KVA）で足す。
+//! 現在は空で、risk の上に置く crate の依存の向きだけを保つ。
 
 #![forbid(unsafe_code)]
 
diff --git a/docs/design/03-architecture.md b/docs/design/03-architecture.md
index dc49026981e1e6305c2445ba7d009d3382dce98f..be4be108e23454cea928f2f4ddca7f8443c631ae 100644
--- a/docs/design/03-architecture.md
+++ b/docs/design/03-architecture.md
@@ -30,6 +30,8 @@ flowchart BT
   cli[rhodium-cli] --> facade
 ```
 
+- 表の責務は 1.0 までの範囲です。
+  まだ無いモジュールと型は、下の[モジュールの地図](#モジュールの地図)と[主要な型の地図](#主要な型の地図)に担当のマイルストーンか issue を書き、各 crate の doc（`lib.rs`）は現在の範囲と、残りが入るマイルストーンを書きます。
 - crate の依存は各 `Cargo.toml` に宣言し、追加時に表の向きを確認します。
   Cargo は依存の循環を拒否します。
 - 商品は market を知らず、market は商品を知りません。
@@ -38,68 +40,80 @@ flowchart BT
 
 ## モジュールの地図
 
-| crate | モジュール |
-|---|---|
-| core | `date`、`time`、`tenor`、`calendar`、`day_count`、`schedule`、`observation`、`currency`、`money`、`quantity`、`fx`、`rounding`、`id`、`refdata`、`fixing`、`brand`、`error`、`parallel` |
-| products | `ir`（`nodes`、`builder`、`validate`、`reduce`、`analysis`、`amendment`、`fallback`）、`dsl`（`syntax`、`check`、`lower`）、`standard`（`rates`、`fra`、`future`、`fx`、`credit`、`equity`、`inflation`、`commodity`、`leg`、`template`、`wire`）、`trade`、`agreement` |
-| market | `quote`、`structure`、`market`、`view`、`curve`、`fx`、`credit`、`inflation`、`asset`、`vol`、`interp` |
-| pricing | `config`、`context`、`resolve`、`definition`、`calibration`、`prepare`、`linear`（`rfr`、`fx`、`future`）、`par`、`analytic`、`model`、`mc`、`amc`、`valuation` |
-| risk | `vector`、`gradient`、`adjoint`、`delta`、`gamma`、`bump`、`scenario`、`theta`、`explain`、`goal_seek`、`rebase`、`simm`、`frtb`、`saccr` |
-| xva | `grid`、`exposure`、`collateral`、`closeout`、`im`、`metrics`、`sensitivity`、`allocation`、`incremental` |
-| rhodium | `prelude`、`workspace`、`session`、`io` |
+| crate | モジュール | まだ無いモジュール（担当） |
+|---|---|---|
+| core | `date`、`time`、`tenor`、`calendar`、`day_count`、`schedule`、`observation`、`currency`、`money`、`quantity`、`fx`、`rounding`、`id`、`refdata`、`fixing`、`brand`、`error`、`parallel` | なし |
+| products | `ir`（`nodes`、`builder`、`validate`、`reduce`、`analysis`、`amendment`、`fallback`）、`dsl`（`syntax`、`check`、`lower`）、`standard`（`rates`、`fra`、`future`、`fx`、`credit`、`equity`、`inflation`、`commodity`、`leg`、`template`、`wire`）、`trade`、`agreement` | `dsl`（M4）、`standard::credit`（M5）、`standard::equity`（M8）、`standard::inflation`（M9）、`standard::commodity`（M10） |
+| market | `quote`、`structure`、`market`、`view`、`curve`、`fx`、`credit`、`inflation`、`asset`、`vol`（`sabr`）、`interp` | `credit`（M5）、`inflation`（M9）、`asset`（M8 と M10）、`interp`（M7 以降。[非線形の補間](06-market/02-curves.md#補間と外挿)とともに足す） |
+| pricing | `config`、`context`、`resolve`、`definition`（`fx`、`vol`）、`calibration`（`groups`、`solver`、`vol`）、`prepare`、`linear`（`rfr`、`fx`、`future`）、`par`、`analytic`、`model`、`mc`、`amc`、`valuation` | `model`、`mc`、`amc`（M4） |
+| risk | `vector`、`gradient`、`adjoint`、`delta`、`gamma`、`bump`、`scenario`、`theta`、`explain`、`goal_seek`、`rebase`、`simm`、`frtb`、`saccr` | `gamma`、`bump`、`scenario`、`theta`、`explain`、`goal_seek`、`rebase`、`simm`、`frtb`（M3）、`saccr`（M6） |
+| xva | `grid`、`exposure`、`collateral`、`closeout`、`im`、`metrics`、`sensitivity`、`allocation`、`incremental` | すべて（`im` は M6、ほかは M5）。現在の crate は空で、依存の向きだけを保つ |
+| rhodium | `prelude`、`workspace`、`session`、`io` | なし |
 
-表は 1.0 までの配置で、まだ無いモジュールは担当するマイルストーンで足します。
-M1 の唯一の解析プライサー（金利先物）は、値洗いの列をリニア評価と共有するため `linear::future` に置き、M2 のスワップション・キャップ・FX オプションの解析プライサーとともに `analytic` へ分けます。
+表は 1.0 までの配置です。
+まだ無いモジュールは、右の列の担当で足します。
+M1 の唯一の解析プライサー（金利先物）は、値洗いの列をリニア評価と共有するため `linear::future` に置いています。
+`analytic` は現在、金利オプションの校正と解析価格が共有する計算（オプションの時間価値の式と、後決め RFR の分散時間 `rfr_caplet_variance_time`）だけを持ちます。
+M2 のスワップション・キャップ・FX オプションの解析プライサーは [#41](https://github.com/kakune/rhodium/issues/41) と [#42](https://github.com/kakune/rhodium/issues/42) で `analytic` に足し、金利先物の解析プライサーもそのときに `analytic` へ分けます。
 
 ## 主要な型の地図
 
 型は、それを定義できる最も下位の crate に置きます。
 表の「文書」が型の意味と不変条件を定めます。
+まだ無い型は「まだ無い型（担当）」の列に、担当のマイルストーンか issue を必ず書きます。
+型を足す、移す、名前を変える PR は、同じ PR でこの表の行を直します（[設計の自己点検](../../AGENTS.md#設計の自己点検)の S14）。
 
-| 型 | 場所 | 文書 |
-|---|---|---|
-| `Date`、`Weekday`、`Tenor`、`Frequency` | core::date、core::tenor | [日付とカレンダー](04-foundation/01-dates-and-calendars.md) |
-| `Instant`、`LocalTime`、`LocalDateTime`、`LocalEventTime`、`TzId`、`ZoneRules`、`ModelTime` | core::time | [日付とカレンダー](04-foundation/01-dates-and-calendars.md#瞬間と現地時刻) |
-| `Calendar`、`BusinessDayConvention` | core::calendar | [日付とカレンダー](04-foundation/01-dates-and-calendars.md) |
-| `DayCount`、`PeriodFrame` | core::day_count | [日数計算](04-foundation/01-dates-and-calendars.md#日数計算) |
-| `ScheduleSpec`、`Schedule`、`Period`、`ExplicitPeriods`、`StubSpec`、`RollConvention`、`PaymentRule`、`ResetRule` | core::schedule | [スケジュール](04-foundation/02-schedules.md) |
-| `ObservationRule`、`ObservationMethod`、`Aggregation` | core::observation | [スケジュール](04-foundation/02-schedules.md) |
-| `Currency`、`Asset`、`Money<R>`、`MoneyBag<R>` | core::currency、core::money | [金額と量](04-foundation/03-money-and-quantities.md) |
-| `Rate<K>`、`Spread`、`YearFraction`、`DiscountFactor<R>`、`SurvivalProbability<R>`、vol と価格の型 | core::quantity | [金額と量](04-foundation/03-money-and-quantities.md) |
-| `CcyPair`、`FxRate<D>`、`FxForwardPoints`、`Rounding` | core::fx、core::rounding | [金額と量](04-foundation/03-money-and-quantities.md) |
-| `VolKind`、`CallPut` | core::quantity | [金額と量](04-foundation/03-money-and-quantities.md) |
-| `IndexId`、`CalendarId`、`ConventionId`、`EntityId`、`SecurityId`、`CommodityId` | core::id | [参照データ](04-foundation/04-reference-data.md) |
-| `RefData`、`IndexDefinition`、`IndexKind`、`PublicationRule`、`Convention`、`FxConvention`、`FxSpotLag`、`CreditKey` | core::refdata、core::id | [参照データ](04-foundation/04-reference-data.md) |
-| `FixingStore`、`ObservationPoint` | core::fixing | [参照データ](04-foundation/04-reference-data.md#fixingstore) |
-| `Brand<'id>` | core::brand | この文書 |
-| `Error`、`Result<T>`、`Warning` | core::error | [参照データ](04-foundation/04-reference-data.md#エラー) |
-| `Contract`、`Product`、`Builder<'b>`、`ContractIx<'b>`、`ExprIx<'b>`、`Party` | products::ir | [コア IR](05-products/01-ir.md) |
-| `SettlementBoundary` | products::ir | [コア IR](05-products/01-ir.md#評価時点での簡約) |
-| `ScriptSource`、`ScriptLibrary`、DSL の構文木 | products::dsl | [DSL](05-products/02-dsl.md) |
-| `StandardProduct`、`ProductKind`、各商品の構造体、`Leg`、`InstrumentTemplate` | products::standard | [標準商品](05-products/03-standard-products.md) |
-| `TradeId`、`BookId`、`MarginSetId` | core::id | [取引](05-products/04-trades.md) |
-| `Trade`、`ProductSpec`、`Amendment`、`ExerciseRecord`、`Fee`、`FallbackTerms`、`FallbackSelection`、`TradeProductResolver`、`ResolvedTrade`、`Portfolio` | products::trade | [取引](05-products/04-trades.md) |
-| `Counterparty`、`NettingSet`、`MarginSet`、`CsaTerms`、`CashCollateralTerms`、`ImTerms` | products::agreement | [取引](05-products/04-trades.md) |
-| `QuoteId`、`QuotePattern`、`QuoteSet`、`QuoteKind`、`QuoteOverride`、`QuoteOverrides` | market::quote（`QuoteId` は quote を指す `Warning` が持てるよう core::id で定義し、market::quote から再公開） | [quote](06-market/01-quotes.md) |
-| `MarketStructure`、`CurveId`、`VolId`、`ModelId`、`CurveSpec`、`CurveExpr`、`CurvePoint`、`ParamId`、`BlockId`、`ParamLayout` | market::structure、market::curve | [曲線と market](06-market/02-curves.md) |
-| `Market<R>`、`MarketShifts`、`MarketView<'v, R>`、`ResolvedCurve<K>`、`CurveRef<'v, K>`、`ScalarId`、`ScalarRef<'v>`、`Slots`、`BoundSlots<'v, R>`、`CurveAt<'v, K, R>`、`ForwardAt<'v, R>`、`GrowthAt<'v, R>`、`CcyScope<'c, 'v, R>`、`PayAt<'c, R>`、`Amount<'c, R>` | market::market、market::view | [曲線と market](06-market/02-curves.md) |
-| `FxSpec`、`FxAnchor`、`FxCollateral`、`FxRef<'v>`、`FxForwardPoint`、`FxForwardInterval` | market::fx | [曲線と market](06-market/02-curves.md#fx-のフォワード) |
-| `VolSpec`、`VolRef<'v, 'a>`、`VolPoint`、`VolAt<'v, 'a>`、`SmileModel`、`SmileDynamics` | market::vol | [ボラティリティ](06-market/04-volatility.md) |
-| `PricingConfig`、`FundingCurves`、`Method`、`AnalyticPricer`、`LinearConvention`、`Context`、`Requirement`、`VolKey`、`CollateralKey` | pricing::config、pricing::context、pricing::resolve | [価格計算](07-pricing/01-pricing.md) |
-| `MarketDefinition`、`InstrumentBinding`、`CalibrationTarget`、`QuoteTransform`、`ElementSelector`、`StepPenalty`、`FxDefinition`、`FxImpliedCurve`、`FxImpliedLongEnd` | pricing::definition | [校正](06-market/03-calibration.md) |
-| `CalibratedMarket`、`CalibrationGroup`、`GroupOutcome`、`GroupSolution`、`HessianMode` | pricing::calibration | [校正](06-market/03-calibration.md) |
-| `Prepared`、`PrepareInput`、`Valuation`、`PortfolioValuation` | pricing::prepare、pricing::valuation | [価格計算](07-pricing/01-pricing.md) |
-| `CashflowPlan`（crate の内部。公開するのは `Prepared` だけ） | pricing::linear | [リニアと解析](07-pricing/02-linear-and-analytic.md) |
-| `PathModel`、`ModelSpec`、`Lgm`、`CrossCurrencyLgm` | pricing::model | [モデル](07-pricing/03-models.md) |
-| `McSettings`、`Kernel`、`PathBatch`、`AmcSettings` | pricing::mc、pricing::amc | [モンテカルロ](07-pricing/04-monte-carlo.md) |
-| `RiskVector<K>`、`ParamGradient`、`DeltaReport`、`HedgeId`、`DeltaOptions`、`GammaMatrix<K>` | risk::vector、risk::delta、risk::gamma | [感応度](08-risk/01-sensitivities.md) |
-| `Scenario`、`QuoteShift`、`MarketShift`、`BlockSelector`、`AssumedEvent`、`RollMode`、`ExerciseAssumption` | risk::scenario、risk::theta | [シナリオと P&L](08-risk/02-scenarios-and-pnl.md) |
-| `HedgeSet`、`HedgeEntry`、`RebaseMetrics`、`Rebased`、`SimmKey`、`CrifRow` | risk::rebase、risk::simm | [付け替えと規制](08-risk/03-hedging-and-regulatory.md) |
-| `TimeGrid`、`GridSpec`、`GridView<'g>`、`GridIx<'g>` | xva::grid | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md) |
-| `MporModel`、`CollateralBalance`、`PendingTransfer` | xva::collateral、xva::closeout | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md) |
-| `XvaMetric`、`XvaSettings`、`ImModel`、`XvaResult`、`XvaReport` | xva::metrics、xva::im | [XVA の指標](09-xva/02-metrics.md) |
-| `XvaRun`、`IncrementalXva` | xva::incremental | [XVA のリスク](09-xva/03-xva-risk.md) |
-| `Workspace`、`Session`、`QuoteDeltaReport`、`Selection` | rhodium::workspace、rhodium::session | [入出力と実行](10-io-and-execution.md) |
+| 型 | 場所 | 文書 | まだ無い型（担当） |
+|---|---|---|---|
+| `Date`、`Weekday`、`Tenor`、`Frequency` | core::date、core::tenor | [日付とカレンダー](04-foundation/01-dates-and-calendars.md) | なし |
+| `Instant`、`LocalTime`、`LocalDateTime`、`LocalEventTime`、`TzId`、`ZoneRules`、`ModelTime` | core::time | [日付とカレンダー](04-foundation/01-dates-and-calendars.md#瞬間と現地時刻) | なし |
+| `Calendar`、`BusinessDayConvention` | core::calendar | [日付とカレンダー](04-foundation/01-dates-and-calendars.md) | なし |
+| `DayCount`、`PeriodFrame` | core::day_count | [日数計算](04-foundation/01-dates-and-calendars.md#日数計算) | なし |
+| `ScheduleSpec`、`Schedule`、`Period`、`ExplicitPeriods`、`StubSpec`、`RollConvention`、`PaymentRule`、`ResetRule` | core::schedule | [スケジュール](04-foundation/02-schedules.md) | なし |
+| `ObservationRule`、`ObservationMethod`、`Aggregation` | core::observation | [スケジュール](04-foundation/02-schedules.md) | なし |
+| `Currency`、`Asset`、`Money<R>`、`MoneyBag<R>` | core::currency、core::money | [金額と量](04-foundation/03-money-and-quantities.md) | なし |
+| `Rate<K>`、`Spread`、`YearFraction`、`DiscountFactor<R>`、`SurvivalProbability<R>`、vol と価格の型 | core::quantity | [金額と量](04-foundation/03-money-and-quantities.md) | なし |
+| `CcyPair`、`FxRate<D>`、`FxForwardPoints`、`Rounding` | core::fx、core::rounding | [金額と量](04-foundation/03-money-and-quantities.md) | なし |
+| `VolKind`、`CallPut` | core::quantity | [金額と量](04-foundation/03-money-and-quantities.md) | なし |
+| `IndexId`、`CalendarId`、`ConventionId`、`EntityId`、`SecurityId`、`CommodityId` | core::id | [参照データ](04-foundation/04-reference-data.md) | なし |
+| `RefData`、`IndexDefinition`、`IndexKind`、`PublicationRule`、`Convention`、`FxConvention`、`FxSpotLag`、`CreditKey` | core::refdata、core::id | [参照データ](04-foundation/04-reference-data.md) | なし |
+| `FixingStore`、`ObservationPoint` | core::fixing | [参照データ](04-foundation/04-reference-data.md#fixingstore) | なし |
+| `Brand<'id>` | core::brand | この文書 | なし |
+| `Error`、`Result<T>`、`Warning` | core::error | [参照データ](04-foundation/04-reference-data.md#エラー) | なし |
+| `Contract`、`Product`、`Builder<'b>`、`ContractIx<'b>`、`ExprIx<'b>`、`Party` | products::ir | [コア IR](05-products/01-ir.md) | なし |
+| `SettlementBoundary` | products::ir | [コア IR](05-products/01-ir.md#評価時点での簡約) | なし |
+| `ScriptSource`、`ScriptLibrary`、DSL の構文木 | products::dsl | [DSL](05-products/02-dsl.md) | すべて（M4） |
+| `StandardProduct`、`ProductKind`、各商品の構造体、`Leg`、`InstrumentTemplate` | products::standard | [標準商品](05-products/03-standard-products.md) | M2 以降の商品の構造体と variant（[標準商品の表](05-products/03-standard-products.md)の各行のマイルストーン） |
+| `TradeId`、`BookId`、`MarginSetId` | core::id | [取引](05-products/04-trades.md) | なし |
+| `Trade`、`ProductSpec`、`Amendment`、`ExerciseRecord`、`Fee`、`FallbackTerms`、`FallbackSelection`、`TradeProductResolver`、`ResolvedTrade`、`Portfolio` | products::trade | [取引](05-products/04-trades.md) | `Portfolio`（M5） |
+| `Counterparty`、`NettingSet`、`MarginSet`、`CsaTerms`、`CashCollateralTerms`、`ImTerms` | products::agreement | [取引](05-products/04-trades.md) | `Counterparty`、`NettingSet`、`ImTerms`（M5） |
+| `QuoteId`、`QuotePattern`、`QuoteSet`、`QuoteKind`、`QuoteOverride`、`QuoteOverrides` | market::quote（`QuoteId` は quote を指す `Warning` が持てるよう core::id で定義し、market::quote から再公開） | [quote](06-market/01-quotes.md) | なし |
+| `MarketStructure`、`CurveId`、`ModelId`、`CurveSpec`、`CurveExpr`、`CurvePoint`、`ParamId`、`BlockId`、`ParamLayout` | market::structure、market::curve | [曲線と market](06-market/02-curves.md) | `ModelId`（M4） |
+| `Market<R>`、`MarketShifts`、`MarketView<'v, R>`、`ResolvedCurve<K>`、`CurveRef<'v, K>`、`ScalarId`、`ScalarRef<'v>`、`Slots`、`BoundSlots<'v, R>`、`CurveAt<'v, K, R>`、`ForwardAt<'v, R>`、`GrowthAt<'v, R>`、`CcyScope<'c, 'v, R>`、`PayAt<'c, R>`、`Amount<'c, R>` | market::market、market::view | [曲線と market](06-market/02-curves.md) | `MarketShifts`（M3）、`ScalarId`、`ScalarRef<'v>`（スカラーを読む最初の資産クラス。M5 の回収率） |
+| `FxSpec`、`FxAnchor`、`FxCollateral`、`FxRef<'v, R>`、`FxForwardPoint`、`FxForwardInterval` | market::fx | [曲線と market](06-market/02-curves.md#fx-のフォワード) | なし |
+| `VolId`、`VolSpec`、`SwaptionCube`、`CapletSurface`、`CapletColumn`、`Sabr`、`SabrFit`、`CubeInterpolation`、`VolPoint`、`VolRef<'v, 'a>`、`VolAt<'v, 'a>`、`SmileModel`、`SmileDynamics` | market::vol（`VolId` は [#200](https://github.com/kakune/rhodium/issues/200) で core::id へ移す） | [ボラティリティ](06-market/04-volatility.md) | `SmileDynamics`（[#43](https://github.com/kakune/rhodium/issues/43)）、`SmileModel`（FX のスマイルとともに M7） |
+| `PricingConfig`、`FundingCurves`、`Method`、`AnalyticPricer`、`LinearConvention`、`Context`、`Requirement`、`VolKey`、`CollateralKey` | pricing::config、pricing::context、pricing::resolve | [価格計算](07-pricing/01-pricing.md) | `VolKey`（[#42](https://github.com/kakune/rhodium/issues/42)） |
+| `MarketDefinition`、`InstrumentBinding`、`CalibrationTarget`、`QuoteTransform`、`ElementSelector`、`StepPenalty`、`FxDefinition`、`FxImpliedCurve`、`FxImpliedLongEnd`、`VolatilityDefinition`、`RateStrike`、`RfrVariance`、`LgmDefinition` | pricing::definition（ボラの指定は pricing::definition::vol） | [校正](06-market/03-calibration.md) | なし |
+| `CalibratedMarket`、`CalibrationGroup`、`GroupOutcome`、`GroupSolution`、`HessianMode`、`VolatilityFit<'a>` | pricing::calibration | [校正](06-market/03-calibration.md) | なし |
+| `Prepared`、`PrepareInput`、`Valuation`、`PortfolioValuation` | pricing::prepare、pricing::valuation | [価格計算](07-pricing/01-pricing.md) | なし |
+| `CashflowPlan`（crate の内部。公開するのは `Prepared` だけ） | pricing::linear | [リニアと解析](07-pricing/02-linear-and-analytic.md) | なし |
+| `PathModel`、`ModelSpec`、`Lgm`、`CrossCurrencyLgm` | pricing::model | [モデル](07-pricing/03-models.md) | すべて（M4） |
+| `McSettings`、`Kernel`、`PathBatch`、`AmcSettings` | pricing::mc、pricing::amc | [モンテカルロ](07-pricing/04-monte-carlo.md) | すべて（M4） |
+| `RiskVector<K>`、`ParamGradient`、`DeltaReport`、`HedgeId`、`DeltaOptions`、`GammaMatrix<K>` | risk::vector、risk::gradient、risk::delta、risk::gamma | [感応度](08-risk/01-sensitivities.md) | `HedgeId`、`DeltaOptions`、`GammaMatrix<K>`（M3） |
+| `Scenario`、`QuoteShift`、`MarketShift`、`BlockSelector`、`AssumedEvent`、`RollMode`、`ExerciseAssumption` | risk::scenario、risk::theta | [シナリオと P&L](08-risk/02-scenarios-and-pnl.md) | すべて（M3） |
+| `HedgeSet`、`HedgeEntry`、`RebaseMetrics`、`Rebased`、`SimmKey`、`CrifRow` | risk::rebase、risk::simm | [付け替えと規制](08-risk/03-hedging-and-regulatory.md) | すべて（M3） |
+| `TimeGrid`、`GridSpec`、`GridView<'g>`、`GridIx<'g>` | xva::grid | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md) | すべて（M5） |
+| `MporModel`、`CollateralBalance`、`PendingTransfer` | xva::collateral、xva::closeout | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md) | すべて（M5） |
+| `XvaMetric`、`XvaSettings`、`ImModel`、`XvaResult`、`XvaReport` | xva::metrics、xva::im | [XVA の指標](09-xva/02-metrics.md) | すべて（`ImModel` は M6、ほかは M5） |
+| `XvaRun`、`IncrementalXva` | xva::incremental | [XVA のリスク](09-xva/03-xva-risk.md) | すべて（M5） |
+| `Workspace`、`Session`、`QuoteDeltaReport`、`Selection` | rhodium::workspace、rhodium::session | [入出力と実行](10-io-and-execution.md) | `Selection`（M3） |
+| pricing の型の再公開: `CalibratedMarket`、`GroupOutcome`、`Valuation`、`PortfolioValuation` | rhodium::workspace（`CalibratedMarket`）、rhodium::session（ほかの三つ） | [入出力と実行](10-io-and-execution.md#session) | なし |
+
+- M1 の LGM は、パスのモデルではなく、金利先物のコンベクシティ調整が読む通貨別のマークです。
+  θ のブロック（`MarketStructure::lgm_block`）と読み出し（`MarketView::lgm_parameters`）が market にあり、マークの quote をそのブロックに結ぶ指定が pricing の `LgmDefinition` です（[金利先物](07-pricing/02-linear-and-analytic.md#金利先物)）。
+  パスのモデルとしての `Lgm` と `CrossCurrencyLgm` は、M4 で pricing::model に足します（[モデル](07-pricing/03-models.md)）。
+- ファサードは、利用者が rhodium-pricing に直接依存せずに結果を扱えるよう、上の pricing の型を再公開します。
+  定義元は pricing のままで、ファサードは別の型を作りません。
 
 ## 型の規約
 
@@ -133,7 +147,7 @@ impl<R: Real> Market<R> {
 
 | lifetime | スコープの開き方 | 守ること | 文書 |
 |---|---|---|---|
-| `'v` | `Market::with_view` | `CurveRef<'v, K>`、`FxRef<'v>`、`VolRef<'v, 'a>` と、`bind_slots` が曲線と点を組にした `CurveAt<'v, K, R>`、`ForwardAt<'v, R>`、`GrowthAt<'v, R>`、`bind_vol` が面と準備済みの点を組にした `VolAt<'v, 'a>` はその view の中でだけ使える。曲線の読み出しは失敗せず、ボラの読み出しはストライクが展開の定義域を外れるとエラーを返す（[ボラティリティ](06-market/04-volatility.md#読み出し)） | [曲線と market](06-market/02-curves.md) |
+| `'v` | `Market::with_view` | `CurveRef<'v, K>`、`FxRef<'v, R>`、`VolRef<'v, 'a>` と、`bind_slots` が曲線と点を組にした `CurveAt<'v, K, R>`、`ForwardAt<'v, R>`、`GrowthAt<'v, R>`、`bind_vol` が面と準備済みの点を組にした `VolAt<'v, 'a>` はその view の中でだけ使える。曲線の読み出しは失敗せず、ボラの読み出しはストライクが展開の定義域を外れるとエラーを返す（[ボラティリティ](06-market/04-volatility.md#読み出し)） | [曲線と market](06-market/02-curves.md) |
 | `'c` | `MarketView::ccy_scope` | `Amount<'c, R>` は同じ通貨スコープの値とだけ足せ、`PayAt<'c, R>`（スコープの割引曲線に結んだ支払いの点）はそのスコープの `pv` にだけ渡せる | [金額と量](04-foundation/03-money-and-quantities.md) |
 | `'b` | `Contract::build` | `ExprIx<'b>`、`ContractIx<'b>` は同じ builder の中でだけ組み合わせられる | [コア IR](05-products/01-ir.md) |
 | `'g` | `TimeGrid::with_view` | `GridIx<'g>` で読むパスの値は同じグリッドの添字でだけ読める | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md) |
@@ -179,7 +193,7 @@ impl<R: Real> Market<R> {
 | 根の探索と方程式系 | scioxide-rootfind（`Brent`、`Toms748`、`SafeguardedNewton`、`NewtonSystem`、`TrustRegionSystem`） |
 | 最小二乗 | scioxide-optimize（`least_squares` の Levenberg–Marquardt） |
 | 線形代数 | scioxide-linalg（LU、QR、Cholesky） |
-| 補間 | 曲線とボラ面の評価の経路は、AD の数値型を通すため rhodium が `R: Real` で書く（monotone convex を含む）。scioxide-interp（`CubicSpline`、`PchipInterpolator`）は f64 の初期値と照合に使う |
+| 補間 | 曲線とボラ面の評価の経路は、AD の数値型を通すため rhodium が `R: Real` で書く。現在の曲線の補間は `FlatForwards` だけで、monotone convex などの非線形の補間は M7 以降に足す（[補間と外挿](06-market/02-curves.md#補間と外挿)）。scioxide-interp（`CubicSpline`、`PchipInterpolator`）は f64 の初期値と照合に使う |
 | 特殊関数 | scioxide-special（`normal_cdf`、`normal_quantile`、`owens_t`） |
 | 乱数と準乱数 | scioxide-random（`Philox4x64`、`QmcPlan`、`BrownianBridgePlan`） |
 | 回帰 | scioxide-stats（`ordinary_least_squares`）か scioxide-linalg の QR |
@@ -209,8 +223,8 @@ impl<R: Real> Market<R> {
   toml、csv、serde_json はファサードと CLI だけが使い、引数の解析の clap は CLI だけが使います。
   開発用には proptest（性質テスト）、criterion（ベンチマーク、`default-features = false`）、tempfile（ファイルのワークスペースのテストと、型テストのケースごとの出力先）を使い、テストの外部 fixture の読み込みには serde_json と toml も使います。
   testkit の型テストのランナーは、rustc と Cargo の JSON の出力と、期待する診断の `case.json` の解析に serde_json を使います。
-- 開発専用の crate `rhodium-testkit` が、テストの補助（`assert_close!`、`assert_mc!`、`check_aad_vs_bump`、`compile_fail_pair!`、決定的な時計 `TickClock`）を一か所に持ちます（[品質](11-quality.md)）。
-  ライブラリの 8 crate はこれに dev-dependency としてだけ依存します。
+- 開発専用の crate `rhodium-testkit` が、テストの補助を一か所に持ちます（一覧と、まだ無い補助の担当は [品質のテストの補助](11-quality.md#テストの補助)）。
+  ライブラリの crate が testkit に依存するときは dev-dependency だけにします（テストを持たない間の rhodium-xva は依存しません）。
 - feature は加算的にし、feature の有無で計算結果を変えません。
 - 公開する値の型は `Send + Sync` です（`Var<'t, _>` とブランド付きの値だけが例外で、スコープの中に閉じます）。
 
@@ -253,7 +267,7 @@ impl<R: Real> Market<R> {
 
 | 仕組み | 形 | 文書 |
 |---|---|---|
-| 準備の再利用 | `prepare` が θ に依存しない `Prepared`（要求を構造の添字に解決、cashflow の計画、線形な曲線の読み出しの補間の重み、MC 用の `Kernel`、読む θ ブロック）を作り、quote、シナリオ、バンプの間で使い回す。`prepared.bind(&view)` で構造の一致とブロックの利用可否を一度だけ確かめる | [価格計算](07-pricing/01-pricing.md) |
+| 準備の再利用 | `prepare` が θ に依存しない `Prepared`（要求を構造の添字に解決、cashflow の計画、線形な曲線の読み出しの補間の重み、MC 用の `Kernel`（M4 で足す）、読む θ ブロック）を作り、quote、シナリオ、バンプの間で使い回す。`prepared.bind(&view)` で構造の一致とブロックの利用可否を一度だけ確かめる | [価格計算](07-pricing/01-pricing.md) |
 | 影響範囲だけの再校正 | `CalibratedMarket::update` が、変わった quote を読む校正グループとその下流だけを、前の θ から warm start で解き直す | [校正](06-market/03-calibration.md) |
 | AAD の active set | 取引が読む θ ブロックだけを tape の変数にし、勾配を疎に持つ | [感応度](08-risk/01-sensitivities.md) |
 | 集計してから随伴 | θ 勾配を集計単位ごとに足してから、保存した分解で随伴を解く | [感応度](08-risk/01-sensitivities.md) |
@@ -262,7 +276,7 @@ impl<R: Real> Market<R> {
 | 照会時刻の追加で変わらないパス | Philox の鍵を日次・日内の橋のノードに付け、時刻を足しても既存のパスを変えない。pre-deal の候補や日付の違うネッティングセットを再シミュレーションなしに評価する | [モンテカルロ](07-pricing/04-monte-carlo.md) |
 | ネッティングセット単位の評価 | バッチを外側、取引を内側のループにし、シミュレーションした状態をキャッシュに載せたまま全取引を評価する | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md) |
 | 閉形式のパス上の評価 | LGM の下のスワップと FX フォワードは条件付きの閉形式で、それ以外は AMC で評価する | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md) |
-| 作業領域の再利用 | ワーカーごとの作業領域 `Scratch`（評価用のバッファ）と tape を使い回す | この文書 |
+| 作業領域の再利用 | ワーカーごとに tape を使い回す。作業領域 `Scratch`（評価用のバッファ）は、ホットループのある MC のカーネルとともに M4 で足す（[並列](10-io-and-execution.md#並列)） | この文書 |
 
 性能の目安と測り方は [入出力と実行](10-io-and-execution.md#性能の目安) にあります。
 
diff --git a/docs/design/06-market/01-quotes.md b/docs/design/06-market/01-quotes.md
index 8d4ce2b4c2c0596919b90f32f27477ad5f7981f8..fb0094ed09cf3ea4244f517335553071c9fa55a2 100644
--- a/docs/design/06-market/01-quotes.md
+++ b/docs/design/06-market/01-quotes.md
@@ -43,7 +43,8 @@ impl FromStr for QuoteId { /* 区間の文字と個数を検査する。違反
 | `EQ.NKY.EQVOL.1Y.K38000` | 日経平均の 1 年、ストライク 38000 のボラ |
 | `MARK.JPY-LGM.MEANREV` | JPY の LGM の平均回帰のマーク |
 
-ストライクの区間は次の字句で書き、市場の定義の binding がそれを座標として解釈します（[ボラティリティ](04-volatility.md)）。
+ストライクの区間は慣行として次の字句で書きます。
+rhodium はこの字句を解釈せず、ストライクの座標は市場の定義の quote の行の `strike` に書きます（[ボラティリティ](04-volatility.md#ボラの種類と座標)）。
 
 | 字句 | 座標 |
 |---|---|
diff --git a/docs/design/06-market/02-curves.md b/docs/design/06-market/02-curves.md
index 8d8384d9a3b5daa28df46bad127f10c1e17f8855..0aac5068e4574898c7901569e4a52832855293fb 100644
--- a/docs/design/06-market/02-curves.md
+++ b/docs/design/06-market/02-curves.md
@@ -8,11 +8,14 @@
 ## 構造と θ
 
 ```rust
-/// 検証済みの市場の構造。曲線、FX、ボラ、スカラー、モデルのパラメータの定義と、θ の配置を持つ。quote に依存しない。
+/// 検証済みの市場の構造。曲線、FX、金利ボラ面の定義と、θ の配置を持つ。quote に依存しない。
 pub struct MarketStructure {
-    objects: Vec<MarketObjectSpec>,   // CurveSpec、FxSpec、VolSpec、ScalarSpec、ModelParamSpec の enum
-    layout: ParamLayout,              // ブロックごとの θ の範囲と要素の制約
     as_of: Instant,
+    curves: Vec<Curve>,                      // 式を展開して θ のブロックの位置に結んだ曲線（予想の測度を含む）
+    grids: Vec<BlockGrid>,                   // ブロックごとの時刻の格子（ピラー、会合、LGM の区切り）
+    fx: BTreeMap<CcyPair, ResolvedFx>,       // 担保別の FX のフォワードと、X(0) のマークの位置
+    layout: ParamLayout,                     // ブロックごとの θ の範囲
+    volatilities: BTreeMap<VolId, VolSpec>,  // 金利ボラ面（スワップションのキューブとキャップレット）
 }
 
 /// θ の 1 要素の名前。感応度の行（θ 空間）とヤコビアンの添字に使う。
@@ -22,19 +25,32 @@ pub struct ParamId { block: BlockId, element: u32 }
 /// 構造と θ。Market<f64> は校正か明示的な θ からだけ作れる。
 pub struct Market<R = f64> {
     structure: Arc<MarketStructure>,
-    theta: Vec<R>,
-    available: BlockSet,          // 今回の校正で利用できるブロック
-    shifts: Arc<MarketShifts>,    // 読み出しに重ねる市場シフト（既定は空）
+    state: Arc<MarketState<R>>,
+}
+
+// θ の持ち方（crate の内部）。読み出しはどちらの形かを意識しない。
+enum MarketState<R> {
+    /// θ の全体と、配置の位置ごとの利用可否。校正やバンプの市場。
+    Dense { theta: Vec<R>, available: Vec<bool> },
+    /// AAD の市場。校正済みの passive な市場を共有し、取引が読むブロックだけを tape の変数で持つ。
+    Active { passive: Market<f64>, blocks: Vec<(usize, Box<[R]>)>, /* 全要素の写しは呼ばれたときだけ作る */ },
 }
 ```
 
-- `shifts` はシナリオの市場シフト（ゼロレートのピラーのシフト、ボラのシフト、ブロックの上書き）を読み出しに重ねる補正です（[シナリオと P&L](../08-risk/02-scenarios-and-pnl.md)）。
+- 構造の部品は、曲線、担保別の FX、金利ボラ面と、金利先物のコンベクシティ調整が読む通貨別の LGM のマーク（`MarketStructure::with_fx_and_lgm`）です。
+  スカラー（回収率など）は M5、パスのモデルのパラメータは M4 で足します（[主要な型の地図](../03-architecture.md#主要な型の地図)）。
+- 市場シフト `MarketShifts`（ゼロレートのピラーのシフト、ボラのシフト、ブロックの上書き）は、M3 のシナリオとともに `Market` の読み出しに重ねる補正として足します（[シナリオと P&L](../08-risk/02-scenarios-and-pnl.md)）。
   構造の `Arc` を共有したまま重ねるので、`bind_slots` を通り、`Prepared` を作り直さずに使えます。
   AAD の `Market<Var>` も同じ補正を持ちます。
 
-- オブジェクトの ID は種類ごとの newtype です（曲線は `CurveId`、ボラ面は `VolId`、モデルのパラメータは `ModelId`）。
-- θ のブロックは「オブジェクトの ID と部品の名前」（例: `JPY-TONA` の `long`）で識別し、要素の制約（`Free`、`Positive`、`NonNegative`、`UnitInterval`、`Correlation`）を持ちます。
-- M0 のイールドカーブは明示した UTC のピラーと、ln DF を並べた `nodes` ブロックを持ちます。ブロックは曲線 ID の辞書順に置きます。
+- オブジェクトの ID は種類ごとの newtype です（曲線は `CurveId`、ボラ面は `VolId`。モデルのパラメータの `ModelId` は M4 で足します）。
+- θ のブロックは「オブジェクトの ID と部品の名前」（例: `JPY-TONA` の `long`）で識別します。
+  現在の `BlockId` は持ち主を曲線の ID の形で持ち、ボラ面のブロックは面の ID、FX のマークは `FX/<ペア>`、LGM は `LGM/<通貨>` をその位置に置きます。
+  持ち主を種類ごとの enum（曲線、ボラ面、FX のマーク、LGM）にする変更は [#200](https://github.com/kakune/rhodium/issues/200) で行います。
+- 要素の制約（`Free`、`Positive`、`NonNegative`、`UnitInterval`、`Correlation`）は、[#200](https://github.com/kakune/rhodium/issues/200) で θ の配置に持たせます。
+  現在の配置は制約を持たず、FX のスポットが正であることなどは校正の側で個別に確かめています。
+- M0 のイールドカーブは明示した UTC のピラーと、ln DF を並べた `nodes` ブロックを持ちます。
+  曲線、FX のマーク、LGM のブロックは `BlockId` の辞書順に置き、ボラ面のブロックはその後ろに `MarketStructure::with_volatilities` で `BlockId` の順に足します（構造を一度に組み立てる形は [#200](https://github.com/kakune/rhodium/issues/200)）。
   割引と予想を兼ねる OIS 曲線では、二つの役割が同じブロックを読みます。
 - θ の各要素を解くグループは一つです。
   一つのグループは複数の商品の残差に同時に適合できます（[校正](03-calibration.md)）。
@@ -49,7 +65,7 @@ pub struct Market<R = f64> {
   勾配は読んだブロックだけの疎な `ParamGradient` です（[感応度](../08-risk/01-sensitivities.md)）。
   `Market<f64>::with_active_theta` は校正済みの構造、θ とブロックの利用可否を共有し、指定したブロックだけを tape の変数に置き換えます。
   ブロックの読み出しは元の市場の利用可否に従い、passive なブロックの値は必要になったときに定数に変換します。
-  取引と報告通貨への FX 換算の AAD に同じ仕組みを使います（[計測](../benchmarks/issue-180.md)）。
+  取引と報告通貨への FX 換算の AAD に同じ仕組みを使います（[計測](../../benchmarks/issue-180.md)）。
 
 ## market view
 
@@ -63,9 +79,11 @@ impl<'v, R: Real> MarketView<'v, R> {
     pub fn discount_curve(&self, id: &CurveId) -> Result<CurveRef<'v, Discount>>;
     pub fn projection_curve(&self, id: &CurveId) -> Result<CurveRef<'v, Projection>>;
     pub fn survival_curve(&self, id: &CurveId) -> Result<CurveRef<'v, Survival>>;
-    pub fn fx(&self, pair: CcyPair, collateral: &CurveId) -> Result<FxRef<'v>>;
-    pub fn vol(&self, id: &VolId) -> Result<VolRef<'v>>;
-    /// FX スポットや回収率のようなスカラーのマーク。
+    pub fn fx(&self, pair: CcyPair, collateral: &CurveId) -> Result<FxRef<'v, R>>;
+    pub fn volatility(&self, id: &VolId) -> Result<VolRef<'v, '_>>;
+    /// 面と準備済みの点を組にする（[ボラティリティ](04-volatility.md#読み出し)）。
+    pub fn bind_vol<'a>(&'a self, face: VolRef<'v, 'a>, point: &'a VolPoint) -> Result<VolAt<'v, 'a>>;
+    /// 回収率のようなスカラーのマーク（M5）。
     pub fn scalar(&self, id: &ScalarId) -> Result<ScalarRef<'v>>;
     /// Slots は market の型。構造が検証した位置と役割を、この view のハンドルと、曲線・面と点の組にする。
     pub fn bind_slots(&self, slots: &Slots) -> Result<BoundSlots<'v, R>>;
@@ -82,7 +100,7 @@ impl<'v, R: Real> MarketView<'v, R> {
 
 - ID からハンドルへの解決は一度だけで、曲線が無い、役割が違う、校正グループが失敗した、のどれもここでエラーになります。
 - M0 の view は割引と予想の曲線（`discount_curve`、`projection_curve`、`df`、`forward`、`growth`）だけを持ちます。
-  生存曲線（`survival_curve`、`survival`）は M5、FX（`fx`、`convert_from`）は M1、ボラ（`vol`）は M2、スカラー（`scalar`）はそれを読む資産クラスのマイルストーンで足します。
+  生存曲線（`survival_curve`、`survival`）は M5、FX（`fx`、`convert_from`）は M1、ボラ（`volatility`、`bind_vol`、`vol_at`）は M2、スカラー（`scalar`）はそれを読む最初の資産クラス（M5 の回収率）で足します。
 - `Slots` は `Arc<MarketStructure>` と検証済みの位置・役割を保持します。
   `bind_slots` は同じ構造から作った市場であることと現在のブロックの利用可否を一度だけ検査します。
 - `MarketStructure::resolve_discount` と `resolve_projection` は、役割を検査した構造の位置 `ResolvedCurve<K>` を返します。
@@ -91,7 +109,8 @@ impl<'v, R: Real> MarketView<'v, R> {
 - `CurveRef<'v, K>` は役割 `K` とブランド `'v` を持ち、同じ view の中でだけ使えます（[アーキテクチャ](../03-architecture.md#ブランド)）。
 - `CurvePoint` は年数、区間の位置、補間の重みを計算済みの点で、`prepare` が作って `Prepared` に保存するので、ループの中で区間を探しません。
   区間や重みは曲線固有なので、計画は曲線の slot と点を対で保存します。
-- `bind_slots` は slot の曲線と点の組を一度だけ検査し、`CurveAt<'v, K, R>`（曲線と一つの点）、`ForwardAt<'v, R>`（予想曲線と期間の両端の点と τ）、`GrowthAt<'v, R>`（予想曲線と期間の両端の点）、`VolAt<'v>`（面と点、[ボラティリティ](04-volatility.md#読み出し)）を返します。
+- `bind_slots` は slot の曲線と点の組を一度だけ検査し、`CurveAt<'v, K, R>`（曲線と一つの点）、`ForwardAt<'v, R>`（予想曲線と期間の両端の点と τ）、`GrowthAt<'v, R>`（予想曲線と期間の両端の点）を返します。
+  ボラの面と点の組 `VolAt<'v, 'a>` は `bind_vol` が作ります（[ボラティリティ](04-volatility.md#読み出し)）。
   `Slots` は割引の点、単利の予想の期間、DF 比の期間を曲線の ID ごとに持ち、DF 比の期間は `Slots::with_projection_growth` で足します。
   読み出しはこれらの組だけを受け取るので、曲線と点を別々に渡して組み替えるコードはコンパイルできず、別の曲線や構造の点を使う組は binding のエラーになります。
   θ について線形な補間では、点は疎な重みを持ちます。曲線の結合では、その重みで ln DF を求め、割引因子や単利フォワードを計算して検査します。
@@ -192,7 +211,7 @@ MaxForward:  f̄_{k,i} = −(L_k(t_{i+1}) − L_k(t_i)) / (t_{i+1} − t_i),   L
 - 部品の式は構造を作るときに展開し、`Ref` の循環と役割の矛盾を拒否します。
 - M1 の Rust API では `CurveSpec::from_expr` に曲線の式と名前付き部品を渡し、ワークスペースの `market/definition.toml` では同じ式を serde の表現のまま `expr` に書きます（[校正](03-calibration.md#ワークスペースの市場定義)）。
   M1 で実装した部品の形は、`Nodes { part, pillars, extrapolation }`（補間は ln DF の線形、`FlatForwards`）、`Stitch { short, long, at }`（`blend` なし）、`Steps { part, meetings }`（発効の瞬間の列）、`Jumps { part, dates }`、`Scale { factor, body }`（有限な定数の倍率）、`FxImplied { part, pair, base, pillars, extrapolation }` です。
-  表の `interp`、`blend`、θ の要素の倍率の `Scale`、会合表と年末の規則からの日付の解決は、それを使う商品や曲線とともに足します（非線形の補間とプロキシの β は M2 以降、ファイルの日付の規則は M6-4）。
+  表の `interp`、`blend`、θ の要素の倍率の `Scale`、会合表と年末の規則からの日付の解決は、それを使う商品や曲線とともに足します（非線形の補間は M7 以降、プロキシの β はそれを使う曲線とともに、ファイルの日付の規則は M6-4）。
   `Scale` の倍率には有限な定数を指定します。
   会合の発効時刻、年末の跳びの時刻、接続点、ピラーを、参照データに従ってあらかじめ解決した UTC の瞬間で指定します。
   接続は境界で値を連続にし、長期側の傾きへ切り替えます。
@@ -212,7 +231,8 @@ MaxForward:  f̄_{k,i} = −(L_k(t_{i+1}) − L_k(t_i)) / (t_{i+1} − t_i),   L
 
 - 原点は値 0 の固定の節点です。
   `ZeroRate` は原点と最初のピラー t₁ の間を z(t) = z(t₁) とし、L(t) = −z(t₁)·t で原点で 0 にします。
-- M1 までの補間は `FlatForwards` だけで、表の他の方式は非線形の補間を要する曲線とともに M2 以降で足します。
+- 現在の補間は `FlatForwards` だけです（M2 の金利ボラも曲線の補間を増やしません）。
+  表の他の方式は、非線形の補間を要する曲線とともに M7 以降で足し、そのときに market の `interp` モジュールを作ります。
 - 最後のピラーより後は `FlatInstantaneousForward`（既定）、`FlatAverageForward`、`Reject` のどれかで外挿します。
 - ピラーは `FromBindings`（そのピラーを受け持つ校正商品がこの曲線を読む最後の日付）か明示の日付・テナーで置きます。
   M1 で実装したのは明示の UTC の瞬間だけで、`FromBindings` とテナーからの解決はファイル形式の見直し（M6-4）で足します（外貨担保の曲線の短期のピラーは FX の受渡の区間から決めます）。
@@ -293,7 +313,7 @@ n_c(t) = ln( F_c(t) / X(0) ) = L_b|c(t) − L_q|c(t)
   短期の quote から X(0) を決めるときは、m のフォワードの両通貨の曲線の一方が、短期の受渡日とスポットの点を校正の目標に含む `FxImplied` の曲線である必要があります。
   そうでない市場の定義（両通貨とも翌日物の曲線だけで m のフォワードが決まり、短期の quote と独立に X(S) を再現できない）は、作るときに `Error::InvalidData` です（テスト `fx_x0_on_tn_consistency`）。
   条件を満たす定義では、食い違いは校正の残差として解消され、値を黙って選び直すことはありません。
-- `FxRef<'v>` は担保 c のフォワード `forward(t)` と、共通の X(0) の `today()`（`FxRate<Today, R>`）を返します。
+- `FxRef<'v, R>` は担保 c のフォワード `forward(t)` と、共通の X(0) の `today()`（`FxRate<Today, R>`）を返します。
   スポット受渡のレートを評価時点の換算に使う誤りは、受渡区分の型で起きません（[金額と量](../04-foundation/03-money-and-quantities.md)）。
 - フォワードポイントは F(T) − X(S) をペアの pip の桁で表した値です。
   スポットより前の短期の quote は隣り合う受渡日の差で、上の連鎖のとおり X(0) を決めます。
@@ -392,7 +412,7 @@ X(t) = ( V − Σ_{0 < e_i ≤ t} D_i · exp(−L(p_i)) ) · exp(L(t)),     L =
 
 ### D6 読み出しは曲線と点の組だけを受け取る
 
-- 判断: `bind_slots` と `CcyScope::at` が、曲線・面と準備した点の組を一度だけ検査して `CurveAt`、`ForwardAt`、`VolAt`、`PayAt` を作り、読み出しと割引はこの組だけを受け取ります。
+- 判断: `bind_slots`、`bind_vol`、`CcyScope::at` が、曲線・面と準備した点の組を一度だけ検査して `CurveAt`、`ForwardAt`、`VolAt`、`PayAt` を作り、読み出しと割引はこの組だけを受け取ります。
 - 理由: 点の区間と重みは曲線ごとに違うので、別の曲線の点で読むと黙って誤った値が出ます。
   有限な θ でも指数や割り算の結果を表現できない場合があります。`bind_slots` と `CcyScope::at` で計算と数値の検査を済ませ、結果を保持することで、評価ループの読み出しは失敗しません。
 - 採らない案: 読み出しのたびに点の曲線を照合する案は、評価の内側に検査と `Result` が入るため採りません。
diff --git a/docs/design/06-market/03-calibration.md b/docs/design/06-market/03-calibration.md
index d1be303bd584bb228ccae584e51c160b5a10ddf7..058120fde6a47953fc80d8862128112a6d277604 100644
--- a/docs/design/06-market/03-calibration.md
+++ b/docs/design/06-market/03-calibration.md
@@ -7,16 +7,22 @@
 
 ## 市場の定義
 
-`MarketDefinition` は `rhodium-pricing::definition` の型で、構造（曲線、FX、ボラ、スカラー、モデルのパラメータ）と、quote を θ に結ぶ instrument binding の列からなります。
+`MarketDefinition` は `rhodium-pricing::definition` の型で、構造（曲線、FX、金利ボラ面、LGM のマーク）と、quote を θ に結ぶ instrument binding の列からなります。
+スカラーは M5、パスのモデルのパラメータは M4 で足します（[主要な型の地図](../03-architecture.md#主要な型の地図)）。
 構造の部分は `MarketStructure` になり（[曲線と market](02-curves.md)）、binding の部分は校正だけが読みます。
 市場の定義は商品を知る必要があるので、商品を知らない `rhodium-market` ではなく評価の crate に置きます。
 
 ```rust
 pub struct MarketDefinition {
-    valuation_calendar: CalendarId,       // 評価日を営業日に調整するカレンダー
-    valuation_zone: TzId,                 // 評価時点を評価日にする地域
-    objects: Vec<MarketObjectSpec>,       // 校正時に評価時点と参照データで構造を具体化する
-    bindings: Vec<InstrumentBinding>,
+    valuation_calendar: CalendarId,               // 評価日を営業日に調整するカレンダー
+    valuation_zone: TzId,                         // 評価時点を評価日にする地域
+    curves: Vec<CurveSpec>,                       // 曲線の式と θ の部品
+    bindings: Vec<InstrumentBinding>,             // 曲線の部品を当てる binding
+    fx: Vec<FxDefinition>,                        // ペアごとの市場入力（スポットとフォワード）、担保別の割引、三角経路
+    volatilities: Vec<VolatilityDefinition>,      // 金利ボラ面。校正時に評価時点と参照データで格子を解決する
+    lgm: BTreeMap<Currency, LgmDefinition>,       // 通貨別の LGM のマーク
+    penalties: BTreeMap<BlockId, StepPenalty>,    // 罰則付きの最小二乗の指定
+    gauss_newton: BTreeSet<BlockId>,              // 解の記録を Gauss–Newton で作るブロック
 }
 
 /// quote を θ の要素に結ぶ指定。校正商品を使う目標は、商品のテンプレートを目標の中に持つ。
@@ -125,17 +131,12 @@ quotes = [
 - `fit = "least_squares"` では各 quote に、`fit = "exact_atm"` ではウィングに `tolerance` を書きます。
   ATM だけなら `fit = { atm_only = { rho = -0.15, nu = 0.4 } }` と書きます。
   `atm_quote` を指定した行では二つの quote を加えて絶対ボラにし、指定した両方の quote から校正を更新できます。
-  `sabr.shift` は正規 SABR 用です。補間の既定は `total_variance` で、`sabr.interpolation = "parameters"` も選べます。
-- 後決め RFR には `vol_kind = "normal"` と `rfr_variance = "lyashenko_mercurio"` を指定します。
+  `sabr.shift` は正規 SABR 用です（シフトの意味は[ボラティリティの D6](04-volatility.md#d6-sabr-のシフトの意味をボラの種類で分ける)）。補間の既定は `total_variance` で、`sabr.interpolation = "parameters"` も選べます。
+- 後決め RFR には `vol_kind = "normal"` と `rfr_variance = "lyashenko_mercurio"` を指定します（[後決め RFR の分散規約](04-volatility.md#後決め-rfr-の分散規約)）。
   ターム金利のキャップには IRS の慣行を使い、タームの fixing 時刻でボラの満期を決めます。
 
 たとえば ATM ボラに `60bp`、キャップの絶対ストライクのボラに `60bp`、SABR の LOW・HIGH にそれぞれ `2bp` と `1bp` を入力できます。
-`rh calibrate <workspace> --date 2026-11-24` は、面 ID・quote ID 順に
-`surface,quote_id,vol_kind,target_vol,model_vol,residual` を標準出力に書きます。
-`--as-of` で瞬間を指定することもできます。
-差分の行でも目標と再現値は絶対ボラです。
-失敗したグループがあれば、成功した行を出力し、診断を標準エラーに書いて終了コード 1 を返します。
-市場を組み立てられない入力エラーの終了コードは 2 です。
+校正したボラ面が各 quote を再現するかは `rh calibrate` で確かめられます（出力の形と終了コードは [CLI の現行のコマンド](../10-io-and-execution.md#現行のコマンド)）。
 
 ### 曲線の部品を使う例
 
@@ -209,7 +210,7 @@ determines = { curve = "JPY-TIBOR-6M", part = "basis" }
   許容幅を持つ binding と持たない binding が一つのグループに混ざる定義は `InvalidInput` です。
 - 曲線、ボラ、モデルのパラメータ、ハザードは同じグラフに入ります。
   quote から XVA のモデルまでが一つの DAG なので、XVA の感応度も quote まで戻せます（[XVA のリスク](../09-xva/03-xva-risk.md)）。
-- 同じ深さの独立なグループは並列に解けます。
+- 同じ深さの独立なグループは並列に解ける構造ですが、現状は依存の順に逐次に解きます（[並列](../10-io-and-execution.md#並列)）。
 
 ```mermaid
 flowchart LR
diff --git a/docs/design/06-market/04-volatility.md b/docs/design/06-market/04-volatility.md
index b320eefd0a349c3d057764415b53585d26d0f66f..37ae7fe418e180b3ddf4344a4c12ff4adff5d2c9 100644
--- a/docs/design/06-market/04-volatility.md
+++ b/docs/design/06-market/04-volatility.md
@@ -1,78 +1,153 @@
 # ボラティリティ
 
 > ボラ面は quote 空間（ATM、ATM からの差、RR、BF、絶対ストライクのボラ）に校正し、SABR や SVI のパラメータはその中間の表現として θ に持ちます。
-> 面の種類はスワップションのキューブ、キャップレット、FX のスマイル、株式とコモディティの面の閉じた enum で、評価コードは `VolRef<'v>` から読みます。
+> 面の種類は閉じた enum `VolSpec` で、現在はスワップションのキューブとキャップレットの二つです。FX、株式、コモディティの面は各マイルストーンで足します。
+> 評価コードは、面のハンドル `VolRef<'v, 'a>` と準備した点を組にした `VolAt<'v, 'a>` から読みます。
 > ベガは quote への感応度として報告し、フォワードが動いたときにボラの何を固定するか（smile dynamics）はリスクの方針として扱います。
 > この文書はボラの種類と座標、面の構造、各面の校正、読み出し、裁定の検査、ベガ、smile dynamics を定めます。
 
 ## ボラの種類と座標
 
-ボラの種類 `VolKind` とオプションの向き `CallPut` は [金額と量](../04-foundation/03-money-and-quantities.md#量の型) が定めます。
+ボラの種類 `VolKind`（`Normal`、`Lognormal`、`ShiftedLognormal { shift }`）とオプションの向き `CallPut` は、商品と市場の両方が使うため rhodium-core の型です（[金額と量](../04-foundation/03-money-and-quantities.md#量の型)）。
+
+金利ボラの quote の座標は、市場の定義の quote の行が持ちます（`rhodium-pricing::definition::vol`）。
 
 ```rust
-// ボラの種類 `VolKind`（Normal、Lognormal、ShiftedLognormal { shift }）と `CallPut` は、
-// 商品と市場の両方が使うため rhodium-core の型（金額と量の文書）。
+/// ストライクの座標。数値は金利の十進表記。
+pub enum RateStrike {
+    Atm,              // フォワードの位置
+    AtmOffset(f64),   // フォワードからの金利差。0 は書かない（Atm を使う）
+    Absolute(f64),    // 金利で直接指定する行使価格
+}
 
-/// ストライクの座標。quote の ID のストライクの字句を binding がこれに解釈する（[quote](01-quotes.md#quoteid)）。
+/// スワップションの一つの quote。満期とテナーは市場の約定日から数える。
+pub struct SwaptionVolQuote {
+    quote: QuoteId,
+    expiry: TenorLength,           // 約定日からオプション満期まで
+    tenor: TenorLength,            // 満期で取引する原資産スワップのテナー
+    strike: RateStrike,
+    atm_quote: Option<QuoteId>,    // ボラの quote が ATM との差のとき、基準の ATM の quote
+    tolerance: Option<f64>,        // 最小二乗の行の当てはめ幅（ボラの単位）
+}
+
+/// キャップのフラットボラの quote。
+pub struct CapVolQuote { quote: QuoteId, maturity: TenorLength, strike: RateStrike /* Atm か Absolute */ }
+```
+
+- ストライクは quote の ID の字句から読まず（[quote の ID](01-quotes.md#quoteid) は意味を読まない文字列です）、quote の行の `strike` に書きます。
+  `Atm` と `AtmOffset` のストライクはフォワード F(θ) の関数で、評価のときに F から作ります。`Absolute` は定数です。
+  キャップの quote は `Atm` と `Absolute` だけを受け付けます。
+- 満期とテナーの `TenorLength` は、校正の入口が評価時点と参照データで瞬間と年数に解決します（約定日、営業日調整、テナーの年数の規則は [ワークスペースの市場定義](03-calibration.md#ワークスペースの市場定義)）。
+  スワップションの満期は面の `expiry_time`、キャップレットの満期は指標の観測の時刻で瞬間にします（[日付を瞬間にする規則](../04-foundation/01-dates-and-calendars.md#日付を瞬間にする規則)、下の [キャップレットのストリッピング](#キャップレットのストリッピング)）。
+  時間の座標は評価時点からのモデル時間 `ModelTime` です。
+- テナーの軸は現在は裸の `f64` の年数で、読み出しの点を作る側も同じ規約で年数を作る必要があります。
+  テナーの座標を型にする変更は [#201](https://github.com/kakune/rhodium/issues/201) の 8.7 で行います。
+- 境界のボラの値（quote、面の節点、結果）の型は `NormalVol`、`LognormalVol`、`ShiftedLognormalVol`（[金額と量](../04-foundation/03-money-and-quantities.md)）で、種類の違うボラを比べたり足したりする誤りは型で起きません。
+  評価の内側の読み出しは `R` で値を計算し、種類の対応は下の [読み出し](#読み出し) のとおり `prepare` が一度だけ確かめます。
+  シフトの違うシフト付き対数正規ボラは比べられません。
+- 面の軸は、満期 × ストライク（キャップレット、FX、株式、コモディティ）か、満期 × テナー × ストライク（スワップション）です。
+
+### 将来の座標
+
+FX、株式、コモディティの面は、金利の三つの座標のほかに、マネーネスとデルタのストライクと、日付や限月の満期を要します。
+
+```rust
+/// 一般のストライクの座標。RateStrike の三つを含み、マネーネスとデルタを加える。
 pub enum StrikeCoordinate {
     Atm,
-    AtmOffset(f64),                         // ATM からの差（decimal）
+    AtmOffset(f64),
     Absolute(f64),
     Moneyness(f64),                         // フォワードに対する比
     Delta { delta: f64, side: CallPut },    // 25DC、25DP
 }
 
-/// 満期の座標。
+/// 一般の満期の座標。
 pub enum Expiry { Tenor(Tenor), Date(Date), Contract(ContractMonth) }
 ```
 
-- 境界のボラの値（quote、面の節点、結果）の型は `NormalVol`、`LognormalVol`、`ShiftedLognormalVol`（[金額と量](../04-foundation/03-money-and-quantities.md)）で、種類の違うボラを比べたり足したりする誤りは型で起きません。
-  評価の内側の読み出しは `R` で値を計算し、種類の対応は下の [読み出し](#読み出し) のとおり `prepare` が一度だけ確かめます。
-  シフトの違うシフト付き対数正規ボラは比べられません。
-- 面の軸は、満期 × ストライク（キャップレット、FX、株式、コモディティ）か、満期 × テナー × ストライク（スワップション）です。
-- RR と BF は二つのデルタのストライクの組の量で、単独の座標を持ちません（下の FX のスマイル）。
-- 満期の日付は、面の慣行の満期の時刻（FX のカット、取引所の満了時刻）で瞬間にし（[日付を瞬間にする規則](../04-foundation/01-dates-and-calendars.md#日付を瞬間にする規則)）、時間の座標は評価時点からのモデル時間 `ModelTime` です。
+- この二つは、最初に要る機能（付け替えの汎関数 `MarketFunctional::Volatility` の M3 か、FX のスマイルの M7 の早い方）とともに足します（[付け替えと規制](../08-risk/03-hedging-and-regulatory.md)）。
+  足すときは `RateStrike` をこの型に置き換え、同じ概念の型を二つ持ちません。
+  限月の満期はコモディティの面（M10）とともに足します。
+- RR と BF は二つのデルタのストライクの組の量で、単独の座標を持ちません（下の [FX のスマイル](#fx-のスマイル)）。
 
 ## 面の構造
 
+面の構造は `rhodium-market::vol` の型で、校正の入口が評価時点で満期とテナーを解決してから作り、`MarketStructure` が面の ID ごとに持ちます（[曲線と market](02-curves.md#構造と-θ)）。
+
 ```rust
-/// ボラ面の定義。MarketStructure の object の一つ。
+/// 評価時点で満期とテナーを解決した金利ボラ面の構造。
 pub enum VolSpec {
-    SwaptionCube { index: IndexId, kind: VolKind, expiries: Vec<Tenor>, tenors: Vec<Tenor>, smile: SmileModel },
-    Caplet { index: IndexId, kind: VolKind, expiries: Vec<Tenor>, strikes: StrikeGrid },
-    FxSmile { pair: CcyPair, expiries: Vec<Tenor>, smile: SmileModel },
-    Asset { asset: Asset, expiries: Vec<Expiry>, smile: SmileModel },
+    SwaptionCube(SwaptionCube),   // 満期 × テナーの SABR キューブ
+    Caplet(CapletSurface),        // 満期 × ストライクのキャップレット面
 }
 
-/// 満期（とテナー）のノードごとのスマイル。
-pub enum SmileModel {
-    /// SABR。β とシフトは定義に固定し、α、ρ、ν を θ に持つ。
-    Sabr { beta: f64, shift: f64, expansion: SabrExpansion /* HaganLognormal | HaganNormal */, fit: SabrFit },
-    /// raw SVI。a、b、ρ、m、σ を θ に持つ。
-    Svi,
-    /// ストライクの節点のボラを直接 θ に持ち、ストライクの方向に補間する。
-    StrikeNodes { interp: StrikeInterp },
+pub struct SwaptionCube {
+    id: VolId,
+    index: IndexId,
+    sabr: Sabr,                        // β、シフト、ボラの種類
+    fit: SabrFit,                      // ノードで校正するパラメータの範囲
+    interpolation: CubeInterpolation,  // ノードの間の補間
+    expiries: Vec<Instant>,            // 解決済みの満期（狭義の昇順、評価時点より後）
+    tenors: Vec<f64>,                  // 原資産のテナー（ACT/365F の年数、狭義の昇順）
 }
 
-pub enum SabrFit { AtmOnly, LeastSquares, ExactAtm }
+pub struct CapletSurface {
+    id: VolId,
+    index: IndexId,
+    kind: VolKind,
+    absolute: Vec<CapletColumn>,       // 絶対ストライクの昇順の列
+    atm: Vec<Instant>,                 // ATM の区間の満期（空なら ATM の列を持たない）
+}
+
+/// 一列分のキャップレットの節点。節点の満期まで同じボラを使う。
+pub struct CapletColumn { strike: Rate<Simple>, ends: Vec<Instant> }
+
+/// SABR のうち定義に固定する値。α、ρ、ν はノードの θ に持つ。
+pub struct Sabr { beta: f64, shift: f64, kind: VolKind }
+
+/// SABR のノードで校正するパラメータの範囲。
+pub enum SabrFit {
+    AtmOnly { rho: f64, nu: f64 },     // α だけを ATM に当て、ρ と ν は定義に固定する
+    LeastSquares,                      // α、ρ、ν をすべての行の重み付き最小二乗で当てる
+    ExactAtm,                          // ATM を厳密に当て、ウィングを重み付き最小二乗で当てる
+}
+
+/// スワップションの満期とテナーの間の補間。
+pub enum CubeInterpolation {
+    TotalVariance,                     // 既定。同じ ATM オフセットの総分散 σ²T を補間する
+    Parameters,                        // α、ρ、ν とノードのフォワードを補間する
+}
 ```
 
-- 満期の間は、同じ座標（同じ対数マネーネス、同じ ATM からの差）の総分散 σ²T を線形に補間するか、スマイルのパラメータを補間するかを定義で選びます。
-- ノードの間で長方形でない格子（満期ごとにストライクの節点が違う）を許します。
+- θ のブロックは面ごとに次の形です。
+
+| 面 | ブロック | 要素の並び |
+|---|---|---|
+| スワップションのキューブ | （満期, テナー）のノードごとに一つ | `AtmOnly` は [ln α, F]、`LeastSquares` と `ExactAtm` は [ln α, b, ln ν, F]（ρ = tanh b）。F はノードのフォワードスワップレート |
+| キャップレット | 面全体で一つ | 絶対ストライクの列ごとの節点のボラを列の順に並べ、その後に ATM の区間のボラ、ATM のストライクを並べる |
+
+- ブロックの並びは、読む側（market の `vol`）と書く側（pricing の校正）が同じ規約で持っています。
+  並びを一か所で定める変更は [#201](https://github.com/kakune/rhodium/issues/201) の 8.7 で行います。
+- キャップレットの列は満期ごとに違う節点を持てるので、長方形でない格子を表せます。
+- キューブの満期とテナーの間は `CubeInterpolation` で補間し、キャップレットは満期の方向に区間ごとに一定、ストライクの方向に線形です（下の [読み出し](#読み出し)）。
+- FX の ATM ボラの面は FX バニラとともに [#42](https://github.com/kakune/rhodium/issues/42) で、FX のスマイルは M7、株式とコモディティの面は M8 と M10 で、`VolSpec` の variant として足します。
+  FX と株式の面のスマイルの選び方（β = 1 の SABR、raw SVI、ストライクの節点）は、それらの面とともに `SmileModel` として定めます（下の各節）。
 
 ## スワップションのキューブと SABR
 
-- キューブのノードは（満期, テナー）の組で、各ノードが自分の SABR のブロック（α、ρ、ν）を持つ独立した校正グループです。
-  ノードは互いに独立なので並列に解けます。
-  校正グループの θ には、同じ曲線で求めたノードのフォワードも派生値として記録します。
+- キューブのノードは（満期, テナー）の組で、各ノードが自分の SABR のブロックを持つ独立した校正グループです。
+  ノードは互いに独立なので並列に解ける構造ですが、現状は逐次に解きます（[並列](../10-io-and-execution.md#並列)）。
+  ブロックには、同じ曲線で求めたノードのフォワード F も派生値として記録します。
   面の読み出しが節点の原資産を評価し直さず、曲線の感応度はこの派生値を通して戻ります。
 - 各ノードの残差は、そのノードのフォワードスワップレートと年金現価を読むので、割引と予想の曲線に依存します（曲線のグループの下流です）。
 - `LeastSquares` は各ストライクのボラを目標にした 3 変数の LM です。
-  ストライクの ATM オフセットとボラの quote の表現は独立で、絶対ボラなら `Identity`、ATM ボラとの差なら `AtmPlusOffset` を使います。
+  ストライクの ATM オフセットとボラの quote の表現は独立で、絶対ボラなら `Identity`、ATM ボラとの差なら `AtmPlusOffset` を使います（quote の行の `atm_quote`）。
   α = exp(a)、ρ = tanh(b)、ν = exp(c) の座標で解きます（[校正](03-calibration.md#解法)）。
   `AtmOnly` は ATM だけに α の 1 変数の求根、`ExactAtm` は ATM を厳密に満たし残りを最小二乗で当てる制約付きの形です。
+- SABR のボラは Hagan の展開で、正規ボラの面は正規の展開、対数正規とシフト付き対数正規の面は対数正規の展開を使います（`Sabr::volatility`）。
+  シフト d はフォワードとストライクの両方に加え、F + d > 0 と K + d > 0 を要します。シフトの意味はボラの種類ごとに違います（[D6](#d6-sabr-のシフトの意味をボラの種類で分ける)）。
 - 正規、対数正規、シフト付き対数正規のキューブは同じ構造で、ボラの種類と quote の値の量だけが違います。
-- CMS の複製のように広いストライクの範囲を読む用途では、SABR の低いストライクの密度が負になりうるので、複製の範囲で密度を検査し、負なら警告します（[リニアと解析](../07-pricing/02-linear-and-analytic.md)）。
+- CMS の複製のように広いストライクの範囲を読む用途では、SABR の低いストライクの密度が負になりうるので、複製の範囲で密度を検査し、負なら警告します（[#42](https://github.com/kakune/rhodium/issues/42)、[リニアと解析](../07-pricing/02-linear-and-analytic.md#cmstsr-の静的複製)）。
 
 ## キャップレットのストリッピング
 
@@ -90,17 +165,36 @@ pub enum SabrFit { AtmOnly, LeastSquares, ExactAtm }
 
 - キャップの最初の期を含めるかは指標の種類が決めます。
   ターム金利（TIBOR）は最初の期の fixing が取引時に確定するので除き、後決めの RFR（TONA、SOFR）は観測期間の間は不確実なので含めます。
-- 後決めの RFR のキャップレットは、観測期間の中で分散が減ることを実効的な分散の時間で扱います（Lyashenko–Mercurio）。
+- キャップレットの満期は、ターム金利では fixing の時刻、後決めの RFR では観測期間の終わりの観測の時刻です。
+  ターム金利の分散の時間は評価時点から満期までのモデル時間で、後決めの RFR は下の規約の実効分散時間です。
+- キャップのベガは、ストリッピングを通してキャップのフラットボラの quote に戻ります。
+- 校正したキャップレットの面でキャップの取引を評価するのは [#41](https://github.com/kakune/rhodium/issues/41)（後決め RFR のキャップレット）と [#42](https://github.com/kakune/rhodium/issues/42)（キャップとフロアの標準商品とターム金利）で、校正と同じ規約と式を使います。
+
+### 後決め RFR の分散規約
+
+後決めの RFR のキャップの quote は、観測期間の中で分散が減ることを実効的な分散の時間で扱う規約（Lyashenko–Mercurio）で値付けされています。
+この規約は、キャップの quote を解釈する面の属性として、市場の定義のキャップレットの面に `rfr_variance = "lyashenko_mercurio"`（`RfrVariance::LyashenkoMercurio`）と書きます（[D5](#d5-後決め-rfr-の分散規約を面の属性にする)）。
+
+観測期間 [T_s, T_e] のキャップレットの、評価時点 t での実効分散時間 v は次の式です（`analytic::rfr_caplet_variance_time`）。
+観測期間の中でボラが (T_e − u)/(T_e − T_s) の割合で線形に減ると仮定し、Bachelier の分散を σ²·v とします。
 
 ```text
-観測期間の前（t ≤ T_s）:  v = (T_s − t) + (T_e − T_s)/3
-観測期間の中（T_s < t < T_e）: v = (T_e − t)³ / (3 (T_e − T_s)²)
+観測期間の前（t ≤ T_s）:        v = (T_s − t) + (T_e − T_s)/3
+観測期間の中（T_s < t < T_e）:  v = (T_e − t)³ / (3 (T_e − T_s)²)
+観測期間の後（t ≥ T_e）:        v = 0
 ```
 
-- キャップのベガは、ストリッピングを通してキャップのフラットボラの quote に戻ります。
+- 時間は評価時点からのモデル時間で、観測期間の端は指標の観測の時刻で瞬間にします。
+  観測期間が正の長さを持たなければ `Error::InvalidInput` です。
+- 後決めの RFR の指標のキャップレットの面は、この規約と正規ボラ（`vol_kind = "normal"`）を要します。
+  ターム金利の指標の面はこの規約を書かず、書けば市場の定義の検査で `InvalidData` です。
+- 校正（[#38](https://github.com/kakune/rhodium/issues/38) のキャップの再価格付け）と取引の評価（[#41](https://github.com/kakune/rhodium/issues/41)）は同じ関数で v を求めるので、校正したボラは取引の評価でも同じ意味を持ちます（テスト `rfr_caplet_variance_matches_hand_calculation`、`caplet_stripping_reprices_term_and_rfr_caps`）。
 
 ## FX のスマイル
 
+M2 の FX の面は ATM ボラだけで、FX バニラ（ATM）とともに [#42](https://github.com/kakune/rhodium/issues/42) で足します。
+この節のスマイル（RR、BF とデルタの慣行）は M7 で足します。
+
 FX の ATM、デルタ、バタフライの慣行は、通貨ペアの慣行（参照データの `FxConvention` の満期の区分ごとの値）が持ちます（[参照データ](../04-foundation/04-reference-data.md)）。
 
 | 慣行 | 値 | 意味 |
@@ -133,6 +227,8 @@ FX の ATM、デルタ、バタフライの慣行は、通貨ペアの慣行（
 
 ## 株式とコモディティの面
 
+株式の面は M8、コモディティの面は M10 で足します。
+
 - 株式の面は対数正規で、満期は日付かテナー、ストライクはフォワードに対するマネーネスか絶対のストライクです。
   フォワードは配当と借株を含む anchor 付きの曲線（[曲線と market](02-curves.md#株式)）なので、現金配当の権利落ち日をまたぐ満期でもマネーネスの座標は連続です。
 - 各満期のスマイルは raw SVI（w(k) = a + b(ρ(k − m) + √((k − m)² + σ²))、k は対数マネーネス、w は総分散）か、β = 1 の SABR で、満期の間は同じ対数マネーネスの総分散を線形に補間します。
@@ -145,41 +241,68 @@ FX の ATM、デルタ、バタフライの慣行は、通貨ペアの慣行（
 ## 読み出し
 
 ```rust
+impl MarketStructure {
+    /// 満期と、キューブならテナー（ACT/365F の年数）を、面の補間の位置と重みを解決した読み出し点にする。
+    pub fn prepare_vol_point(self: &Arc<Self>, id: &VolId, expiry: Instant, tenor: Option<f64>) -> Result<VolPoint>;
+}
+
 impl<'v, R: Real> MarketView<'v, R> {
+    /// ID から面のハンドルを得る。
+    pub fn volatility(&self, id: &VolId) -> Result<VolRef<'v, '_>>;
+    /// 面と点の組を検査する。別の面や構造の点、利用できないブロックはエラー。
+    pub fn bind_vol<'a>(&'a self, face: VolRef<'v, 'a>, point: &'a VolPoint) -> Result<VolAt<'v, 'a>>;
     /// 面と結んだ点で、絶対ストライクのボラを読む。
     /// ストライクが展開の定義域から外れるとエラーになる。
     pub fn vol_at(&self, at: VolAt<'v, '_>, strike: R) -> Result<R>;
-    /// ボラとストライク・満期の微分（局所ボラの Dupire の式に使う）。
-    pub fn smile_at(&self, at: VolAt<'v, '_>, strike: R) -> SmileValue<R>;
+}
+```
+
+- `VolPoint` は満期と（スワップションの）テナーの断面と補間の重みを計算済みの点です。
+  キューブの点は、満期のモデル時間とテナーの年数の両方向の線形の重みで最大 4 つのノードを持ち、軸の外は端のノードを使います。
+  キャップレットの点は、各列と ATM の列で、満期を含む区間（最後の節点より後は最後の区間）の位置を持ちます。
+- 現在 `VolPoint` を作るのは、キャップレットの校正（キャップの再価格付け）だけです。スワップションのノードの校正は、点を作らずにノードの SABR の式を直接評価します。
+  `prepare` が取引の要求を面に解決して点を作り、`Prepared` に保存する経路は、解析プライサーとともに [#42](https://github.com/kakune/rhodium/issues/42) で足します（[#41](https://github.com/kakune/rhodium/issues/41) の後決め RFR のキャップレットも同じ経路を使います）。
+  `VolAt<'v, 'a>` は面のハンドルと点を `MarketView::bind_vol` が組にして検査した値で、読み出しはこの組だけを受け取ります（[曲線と market](02-curves.md#market-view)）。
+  点を `Slots` に持たせて `bind_slots` で他の点とまとめて結ぶ形は、準備のときに位置を解決した点を結ぶ変更（[#200](https://github.com/kakune/rhodium/issues/200)）とともに入れます。
+  別の面の点で読む誤りは binding のエラーになり、組を作らずに読むコードはコンパイルできません。
+- ストライクは面に依存しない数なので読み出しの引数にし、フォワードに依存するストライクも評価のときに渡せます。
+- キューブの読み出しは `CubeInterpolation` に従います。
+  `TotalVariance` は、点のフォワード F（ノードのフォワードの重み付きの和）からのオフセット K − F を各ノードで同じに保ち、ノードの SABR のボラ σ_i(F_i, F_i + K − F) の総分散 σ_i²·T_i を重みで足して、点の満期 T で割ります。
+  `Parameters` は、α、ρ、ν とフォワードを重みで足してから、点の満期で SABR のボラを評価します。
+- キャップレットの読み出しは、各列の（ストライク, ボラ）と ATM の（θ のストライク, θ のボラ）の点から、ストライクの左右で最も近い二点を線形に補間し、外側は端の値にします。
+  ATM のストライクは θ に依存するので、二点の選択は数値型のレーンごとに `R::select` で行います。
+- `VolRef<'v, 'a>` と `VolAt<'v, 'a>` はブランド `'v` の view の中でだけ使えます。
+- ボラの種類は面の定義の `VolKind` が持ち、プライサーは種類を見て Black、Bachelier、シフト付き Black を選びます。
+  種類とプライサーの基本式の対応は、`prepare` が要求を面に解決するときに一度だけ確かめ、対応しない組（正規ボラしか使えない式に対数正規の面など）は `Error::InvalidInput` です（テスト `vol_kind_formula_mismatch`。[#42](https://github.com/kakune/rhodium/issues/42) で足す）。
+  評価の内側の読み出しでは `R` で値を計算し、種類を型では区別しません。
+- ストライクがフォワードに依存する座標（ATM、ATM からの差、マネーネス）では、評価のときのフォワード（θ の関数）からストライクを作って読み出しに渡すので、曲線の動きがボラに入ります。
+- ボラとストライク・満期の微分を返す読み出しは、CMS の複製の密度の検査とともに [#42](https://github.com/kakune/rhodium/issues/42) で足します。
+
+```rust
+impl<'v, R: Real> MarketView<'v, R> {
+    /// ボラとストライク・満期の微分（複製の密度の検査と、局所ボラの Dupire の式に使う）。
+    pub fn smile_at(&self, at: VolAt<'v, '_>, strike: R) -> Result<SmileValue<R>>;
 }
 
 pub struct SmileValue<R> { vol: R, d_vol_dk: R, d2_vol_dk2: R, d_total_variance_dt: R }
 ```
 
-- `VolPoint` は満期と（スワップションの）テナーの断面と補間の重みを計算済みの点で、`prepare` が面ごとに作って `Prepared` に保存します。
-  `VolAt<'v, 'a>` は面のハンドルと点を `MarketView::bind_vol` が組にして検査した値で、読み出しはこの組だけを受け取ります（[曲線と market](02-curves.md#market-view)）。
-  点を `Slots` に持たせて `bind_slots` で他の点とまとめて結ぶ形は、ボラを読むプライサーを足すときに加えます。
-  別の面の点で読む誤りは binding のエラーになり、組を作らずに読むコードはコンパイルできません。
-- ストライクは面に依存しない数なので読み出しの引数にし、フォワードに依存するストライクも評価のときに渡せます。
-- `VolRef<'v, 'a>` と `VolAt<'v, 'a>` はブランド `'v` の view の中でだけ使えます。
-- ボラの種類は面の定義の `VolKind` が持ち、プライサーは種類を見て Black、Bachelier、シフト付き Black を選びます。
-  種類とプライサーの基本式の対応は、`prepare` が要求を面に解決するときに一度だけ確かめ、対応しない組（正規ボラしか使えない式に対数正規の面など）は `Error::InvalidInput` です（テスト `vol_kind_formula_mismatch`）。
-  評価の内側の読み出しでは `R` で値を計算し、種類を型では区別しません。
-- ストライクがフォワードに依存する座標（ATM、ATM からの差、マネーネス）では、評価のときのフォワード（θ の関数）からストライクを作って読み出しに渡すので、曲線の動きがボラに入ります。
-
 ## 裁定の検査
 
-| 検査 | 対象 | 意味 |
-|---|---|---|
-| バタフライ | すべての面の各ノード | ストライクの格子でリスク中立の密度 ∂²C/∂K² が負にならない。SVI では g(k) ≥ 0 |
-| カレンダー | 原資産が満期によらない面（FX、株式） | 同じマネーネスの総分散 σ²T が満期について減らない |
-| 満期の方向の単調性 | 原資産が満期ごとに違う面（スワップション、キャップレット、コモディティ） | 無裁定の条件ではない経験的な検査 |
+| 検査 | 対象 | 意味 | 実装 |
+|---|---|---|---|
+| バタフライ | すべての面の各ノード | ストライクの格子でリスク中立の密度 ∂²C/∂K² が負にならない。SVI では g(k) ≥ 0 | スワップションのノードとキャップレットの面は、校正のストライクの格子でコールの価格の凸性を確かめ、`Warning::VolatilityButterfly` を返す。SVI は M8 |
+| カレンダー | 原資産が満期によらない面（FX、株式） | 同じマネーネスの総分散 σ²T が満期について減らない | M7（FX のスマイル） |
+| 満期の方向の単調性 | 原資産が満期ごとに違う面（スワップション、キャップレット、コモディティ） | 無裁定の条件ではない経験的な検査 | 満期の方向の検査の仕組みとともに M7 |
 
 - 違反は `Warning` にし、校正を止めません。
   面を使うかどうかは市場の定義の判断です。
 
 ## ベガ
 
+ベガを quote の単位で出し、バンプと再校正の値と照合するのは [#44](https://github.com/kakune/rhodium/issues/44) です。
+校正はボラのグループの線形化も記録しているので、ボラを読む取引（#41、#42）が入れば、ベガは quote デルタと同じ随伴で求まります。
+
 - ベガはボラの quote への感応度として報告します。
   SABR や SVI のパラメータへの感応度は、校正の随伴で ATM と各ストライクの quote に戻ります（[感応度](../08-risk/01-sensitivities.md)）。
 - キャップのベガはキャップのフラットボラの quote に、FX のベガは ATM、RR、BF の quote に戻ります。
@@ -187,6 +310,9 @@ pub struct SmileValue<R> { vol: R, d_vol_dk: R, d2_vol_dk2: R, d_total_variance_
 
 ## smile dynamics
 
+この節の方針 `SmileDynamics` は [#43](https://github.com/kakune/rhodium/issues/43) で足します。
+それまでの quote デルタは、既定の `StickyQuote`（校正の随伴そのもの）だけです。
+
 市場は一つだけ校正し、フォワードが動いたときにボラの何を固定するかは、quote デルタの方針 `SmileDynamics` として扱います（[感応度](../08-risk/01-sensitivities.md)）。
 quote と時間を動かすシナリオはボラの quote を固定して再校正するので、既定の `StickyQuote` に従います（[シナリオと P&L](../08-risk/02-scenarios-and-pnl.md)）。
 市場のシフトは再校正しないので、曲線をシフトしてもボラの θ（面のパラメータ）はそのままです。
@@ -210,6 +336,7 @@ quote のストライクの座標によらず、同じ面のすべての quote 
 | デルタで quote する面（FX） | デルタ | デルタを固定して解き直したストライク |
 
 - この表が `StickyMoneyness` の唯一の定義で、quote デルタの随伴とバンプによる照合は同じ規則を使います（[感応度](../08-risk/01-sensitivities.md#smile-dynamics)）。
+- シフト付き対数正規の面のシフト d は面の `VolKind::ShiftedLognormal { shift }` で、正規 SABR のモデル用のシフトはこの表に入りません（[D6](#d6-sabr-のシフトの意味をボラの種類で分ける)）。
 
 - 既定の `StickyQuote` は校正の随伴をそのまま使います（曲線の quote への感応度は、曲線の動きでフォワードが動き、ボラの quote を固定したまま面が付き直る効果を含みます）。
 - quote デルタでは、各方針を、上流が動いたときのボラのパラメータの応答 B = ∂θ_v/∂θ_up として求めます（[感応度](../08-risk/01-sensitivities.md#smile-dynamics)）。
@@ -244,8 +371,42 @@ quote のストライクの座標によらず、同じ面のすべての quote 
 - 採らない案: 方針ごとの残差のヤコビアンを J_up にだけ差し替える案は、方針がボラへのヤコビアンも変えるとき、保存した分解との組が誤った応答を返すため採りません。
   方針ごとの残差を校正の随伴に組み込む案は、随伴の形が方針ごとに変わり複雑になるため採りません。
 
-### D4 面の種類を閉じた enum にする
+### D4 面の種類を閉じた enum にし、スマイルのモデルは面の構造体が持つ
 
-- 判断: 面の種類を `VolSpec` の variant、スマイルを `SmileModel` の variant にします。
+- 判断: 面の種類を `VolSpec` の variant にします。
+  スマイルのモデルは面の構造体のフィールドにし、金利のキューブは `Sabr` と `SabrFit` を直接持ちます。
+  モデルを選べる面（FX と株式）は、その面とともに `SmileModel` の enum を足します。
 - 理由: 面の種類は少なく、`match` で網羅でき、種類を足す場所が一つに決まります。
+  金利のキューブのスマイルは SABR だけなので、選べない選択肢の enum を先に置くと、使わない variant の検査と分岐が増えます。
 - 採らない案: 原資産、種類、軸の三つの型引数で面の型を作る案は、型が複雑になる割に防ぐ誤りが少ないため採りません（種類とプライサーの式の対応は `prepare` が一度だけ確かめます）。
+  すべての面に共通の `SmileModel` を M2 で置く案は、キューブが使わない SVI とストライクの節点の variant を持つことになるため採りません。
+
+### D5 後決め RFR の分散規約を面の属性にする
+
+- 判断: 後決めの RFR のキャップレットの実効分散時間（Lyashenko–Mercurio）は、キャップの quote を解釈する面の属性 `RfrVariance::LyashenkoMercurio` として市場の定義に書きます。
+  リニア評価の近似の規約（`LinearConventions`）には入れません。
+  校正のキャップの再価格付けと、取引のキャップレットの評価（[#41](https://github.com/kakune/rhodium/issues/41)）は、同じ規約と同じ関数 `analytic::rfr_caplet_variance_time` を使います。
+- 理由: 実効分散時間は、市場がキャップの quote を値付けするときの約束で、面のボラはこの約束の下でだけ意味を持ちます。
+  面の属性にすると、校正したボラと取引の評価が同じ約束を使うことが構造で保証され、規約を選ぶ場所が市場の定義の一か所になります。
+  ターム金利の面に書けば検査で拒否し、後決めの RFR の面に書き忘れても拒否するので、規約の無い RFR の面は作れません。
+- 採らない案: 近似の規約 `RfrCapletLmVariance` として価格設定で許可・禁止する案は採りません。
+  禁止すると、同じ面のボラを別の分散で読むことになり、校正と取引の評価が食い違った数字を黙って出します。
+  また、この規約は厳密なモデルの近似ではなく quote の約束なので、禁止したときに切り替える「厳密な評価」がありません。
+
+### D6 SABR のシフトの意味をボラの種類で分ける
+
+- 判断: SABR のシフトは、ボラの種類ごとに次の意味にします。
+
+| ボラの種類 | シフトの出どころ | 意味 |
+|---|---|---|
+| `ShiftedLognormal { shift }` | 面の `vol_kind = "shifted_lognormal"` と `shift` | quote のシフト。SABR の展開も同じ値を使う |
+| `Lognormal` | なし | 0 |
+| `Normal` | キューブの `sabr.shift`（既定 0） | 正規 SABR のモデル用のシフト。Hagan の展開が要する F + d > 0 と K + d > 0 を負の金利でも満たすための量で、quote（正規ボラ）はシフトを持たない |
+
+- 入力の検査: 対数正規とシフト付き対数正規の面に 0 でない `sabr.shift` を書くと `InvalidData` にし、二つの入力が同じ量を指す食い違いを読み込みで拒否します。
+- 理由: シフト付き対数正規ボラの quote はシフトと組で初めて意味を持つので、シフトは面の種類 `VolKind` が一度だけ持つべきです。
+  正規ボラの quote はシフトを持たないため、正規 SABR のシフトは quote と独立なモデルの入力です。
+- 現状と後続: 現在の `Sabr` はシフト付き対数正規のシフトを `kind` と `shift` の二つに持ち、`Sabr::new` で一致を検査しています。
+  この二重の保持は [#201](https://github.com/kakune/rhodium/issues/201) の 8.7 で解消し、`Sabr` はシフト付き対数正規のシフトを `VolKind` から読み、正規 SABR のシフトだけを自分で持ちます。
+  入力の形式（面の `shift` と `sabr.shift`）と上の表の対応は変えません。
+- 採らない案: SABR に一つの `shift` を常に持たせ、種類に関係なく展開に使う案は、シフト付き対数正規の quote のシフトと SABR のシフトが別々に書けてしまい、食い違うと同じ面の quote と読み出しが違うシフトの Black で解釈されるため採りません。
diff --git a/docs/design/07-pricing/02-linear-and-analytic.md b/docs/design/07-pricing/02-linear-and-analytic.md
index ebef6049c714b0dd54837f08635940ca4341e66d..11f8a6e8806a94293c881ea706afba980491f05c 100644
--- a/docs/design/07-pricing/02-linear-and-analytic.md
+++ b/docs/design/07-pricing/02-linear-and-analytic.md
@@ -101,9 +101,12 @@ RFR の複利 Π_i(1 + w_i·r_i) を端点の割引因子の比 P_p(s)/P_p(e) 
 | `MtmFxForwardProduct` | MtM XCCY の元本リセットの FX と金利の共分散 | なし |
 | `DeterministicCarry` | 金利と配当の確率的な変動による資産価格の比の調整 | なし |
 | `IrrAnnuityAtForward` | IRR 差金決済のアニュイティの凸性 | なし |
-| `RfrCapletLmVariance` | RFR キャップレットの観測期間内のボラの形（線形に減ると仮定） | キャップの quote がこの規約で値付けされている |
 | `RoundingAtForward` | 丸めの非線形性（丸めた値の期待値と期待値を丸めた値の差） | 取引条件が丸めを宣言している |
 
+- 実装の `LinearConvention` にまだ無い規約は、`IrrAnnuityAtForward`（IRR 差金決済のスワップションとともに [#42](https://github.com/kakune/rhodium/issues/42)）と `DeterministicCarry`（資産のフォワードとともに M8）です。
+- 後決め RFR のキャップレットの観測期間内の分散（Lyashenko–Mercurio）は、この表の近似の規約ではありません。
+  キャップの quote を解釈する面の規約 `RfrVariance::LyashenkoMercurio` が決め、価格設定では許可・禁止しません（[ボラティリティの D5](../06-market/04-volatility.md#d5-後決め-rfr-の分散規約を面の属性にする)）。
+
 ### 担保の違いの判定
 
 `CrossCollateralProjection` は、予想の担保と支払いの担保が違うときだけ記録します。
@@ -163,6 +166,7 @@ R への復元:              scioxide-ad の implicit_solution_into(p*, PV)   p*
 IR の形を認識する仕組みは持ちません。
 
 ```rust
+// 現在の variant は RateFuture だけで、ほかは下の表のマイルストーンで足す。
 pub enum AnalyticPricer {
     RateFuture,
     Swaption,
@@ -178,11 +182,11 @@ pub enum AnalyticPricer {
 | プライサー | 対象の標準商品 | 式 | マイルストーン |
 |---|---|---|---|
 | `RateFuture` | 短期金利先物 | 先物価格（フォワードに LGM の閉形式のコンベクシティ調整）と次の値洗い | M1 |
-| `Swaption` | ヨーロピアンスワップション（現物、IRR 差金、担保付き現金） | N·A(0)·π(F, K, σ, T_e) | M2 |
-| `CapFloor` | キャップとフロア（ターム金利、後決め RFR） | キャップレットごとの基本式 | M2 |
-| `CmsTsr` | CMS クーポン、CMS キャップとフロア | TSR の静的複製 | M2 |
-| `CmsSpreadCopula` | CMS スプレッドオプション | 周辺分布のガウスコピュラ | M2 |
-| `FxVanilla` | FX ヨーロピアン | Garman–Kohlhagen | M2（ATM）、M7（スマイル） |
+| `Swaption` | ヨーロピアンスワップション（現物、IRR 差金、担保付き現金） | N·A(0)·π(F, K, σ, T_e) | M2（[#42](https://github.com/kakune/rhodium/issues/42)） |
+| `CapFloor` | キャップとフロア（ターム金利、後決め RFR） | キャップレットごとの基本式 | M2（後決め RFR のキャップレットは [#41](https://github.com/kakune/rhodium/issues/41)、標準商品とターム金利は [#42](https://github.com/kakune/rhodium/issues/42)） |
+| `CmsTsr` | CMS クーポン、CMS キャップとフロア | TSR の静的複製 | M2（[#42](https://github.com/kakune/rhodium/issues/42)） |
+| `CmsSpreadCopula` | CMS スプレッドオプション | 周辺分布のガウスコピュラ | M2（[#42](https://github.com/kakune/rhodium/issues/42)） |
+| `FxVanilla` | FX ヨーロピアン | Garman–Kohlhagen | M2（ATM、[#42](https://github.com/kakune/rhodium/issues/42)）、M7（スマイル） |
 | `AssetVanilla` | 株式、株価指数、コモディティのヨーロピアン | 資産フォワードの Black（コモディティの正規ボラは Bachelier） | M8、M10 |
 | `InflationYoy` | YoY スワップ、YoY キャップとフロア | 比のフォワードと凸性の調整 | M9 |
 
@@ -212,12 +216,9 @@ Black:         π = ω(F·Φ(ω·d₁) − K·Φ(ω·d₂))                  d
 
 - ターム金利のキャップレット: P_d(t_pay)·τ·N·π(F_p(s, e), K, σ(T_fix, K), T_fix)。
   支払日が e と違えば `PaymentDelayIgnored` を使います。
-- 後決め RFR のキャップレット（Lyashenko–Mercurio）: 観測期間 [T_s, T_e] の中でボラが (T_e − u)/(T_e − T_s) の割合で線形に減ると仮定し、Bachelier の分散を σ²·v とします。
-
-```text
-t ≤ T_s:         v = (T_s − t) + (T_e − T_s)/3
-T_s < t < T_e:   v = (T_e − t)³ / (3·(T_e − T_s)²)
-```
+- 後決め RFR のキャップレット: 面の規約 `RfrVariance::LyashenkoMercurio` の実効分散時間 v で、Bachelier の分散を σ²·v とします。
+  v の式と規約の意味は [後決め RFR の分散規約](../06-market/04-volatility.md#後決め-rfr-の分散規約) が定め、校正のキャップの再価格付けと取引の評価（[#41](https://github.com/kakune/rhodium/issues/41)）は同じ関数 `analytic::rfr_caplet_variance_time` を使います。
+  σ は面の観測期間の終わりの満期のボラで、面が正規ボラであることは市場の定義の検査が保証します。
 
 - 観測期間に入ったキャップレットのフォワードは、確定した日の積と残りの期間の telescoping の比の積です。
 
diff --git a/docs/design/10-io-and-execution.md b/docs/design/10-io-and-execution.md
index 8b370b43475cf467d42fdc0fe03d78318d93dcfb..c59cebe0182d3ac91a7ce9065e2340841db5b3ed 100644
--- a/docs/design/10-io-and-execution.md
+++ b/docs/design/10-io-and-execution.md
@@ -318,12 +318,16 @@ let xva = session.xva(&Selection::netting_set("NS-ACME-1"), execution)?;  // Xva
 
 `rh` はワークスペースを読み、結果を `output/` か標準出力に書く薄いコマンドです。
 
-### M1 までのコマンド
+### 現行のコマンド
 
-- `rh price <dir>` は `trade_id,currency,pv` を取引 ID、通貨の順に標準出力へ書きます。
-- `rh risk <dir>` は `quote_id,currency,sensitivity` を quote ID、通貨の順に標準出力へ書きます。感応度は全取引の通貨別の合計で、[quote の種類が定める変化幅](06-market/01-quotes.md#感応度の単位)あたりの額です。
-- どちらも `--as-of` または `--date` を指定します。地域を使う場合は `--as-of '2026-11-24T17:00:00[Asia/Tokyo]'` の形で指定します。`--date 2026-11-24` は価格設定の `valuation_close` で瞬間に解決します。
-- `--threads N`（既定 1）は読み込みから評価まで取引を N 個のワーカーに分けます。1 は逐次で、0 は引数の誤り（終了コード 2）です。取引ごとの PV は N によらず同じで、デルタは通貨別の和の順だけが変わるので [許容誤差](11-quality.md#許容誤差) の範囲で一致します。
+| コマンド | 標準出力に書く CSV | 並び | 入ったマイルストーン |
+|---|---|---|---|
+| `rh price <dir>` | `trade_id,currency,pv`。取引ごと、通貨ごとの PV | 取引 ID、通貨の順 | M0 |
+| `rh risk <dir>` | `quote_id,currency,sensitivity`。全取引の通貨別の合計の quote デルタで、[quote の種類が定める変化幅](06-market/01-quotes.md#感応度の単位)あたりの額 | quote ID、通貨の順 | M0 |
+| `rh calibrate <dir>` | `surface,quote_id,vol_kind,target_vol,model_vol,residual`。ボラ面の quote ごとの目標のボラ、校正した面で再価格付けしたボラ、残差（再現値 − 目標）。差分の quote の行でも目標と再現値は絶対ボラ（[ボラの校正の書き方](06-market/03-calibration.md#ワークスペースの市場定義)） | 面 ID、quote ID の順 | M2 |
+
+- どのコマンドも `--as-of` または `--date` を指定します。地域を使う場合は `--as-of '2026-11-24T17:00:00[Asia/Tokyo]'` の形で指定します。`--date 2026-11-24` は価格設定の `valuation_close` で瞬間に解決します。
+- `--threads N`（既定 1）は読み込みから評価まで取引を N 個のワーカーに分けます。1 は逐次で、0 は引数の誤り（終了コード 2）です。`rh calibrate` は取引を評価しないので、N は読み込みにだけ効きます（校正は逐次です。[並列](#並列)）。取引ごとの PV は N によらず同じで、デルタは通貨別の和の順だけが変わるので [許容誤差](11-quality.md#許容誤差) の範囲で一致します。
 - 準備と校正の警告（近似の規約、校正の当てはめの許容誤差の超過など）は `warning: trade <ID>: <警告>` と `warning: calibration: <警告>` の形で標準エラーに書き、終了コードには影響しません。
 - `--format`、`--by` と `check`、`xva` のコマンドは下の拡張で足します（M3 以降。CLI の形は M6 の公開 API の見直しで確定します）。
 - 診断は標準エラーに書きます。校正や取引の評価が部分的に失敗したときは、計算できた行を出力して終了コード `1` を返します。一つの通貨の校正商品の準備の失敗（fixing の欠損など）も、その曲線のグループの失敗として扱い、独立した通貨の行を出力します。入力のエラーは `2`、すべて成功した場合は `0` です。評価時点が曲線の最初のピラー以降にあるなど、入力どうしの不整合で市場を組み立てられない場合も入力のエラーとして扱い、CSV を出力しません。
@@ -372,8 +376,11 @@ let xva = session.xva(&Selection::netting_set("NS-ACME-1"), execution)?;  // Xva
   共有の参照数の更新はワーカーの間で競合し、並列の伸びを止めます（M0 の準備で 20 スレッドの時間が半分になった例は [M0 の計測](../benchmarks/m0.md)）。
 - 縮約の順序は固定しません。
   スレッド数を変えた結果の一致は [許容誤差](11-quality.md#許容誤差) で判定します。
-- 校正はグループの依存の順に解き、依存し合わないグループは並列に解けます。
-  M1 までの校正は逐次です（最初の検証データセットの全 9 グループの校正が 1 コアで目安の 100 ms を大きく下回るため、[M1 の計測](../benchmarks/m1.md)）。グループの並列はボラ面の校正が重くなる M2 以降で測って入れます。
+- 校正はグループを依存の順に解きます。
+  依存し合わないグループ（同じ深さのグループや、スワップションのキューブの各ノード）は並列に解ける構造ですが、現状の校正は逐次です。
+  `calibrate` と `CalibratedMarket::update` は実行器（`ExecutionContext`）を受けず、グループを依存の順に一つずつ解きます。
+  M1 の最初の検証データセットの全 9 グループの校正は、1 コアで目安の 100 ms を大きく下回ります（[M1 の計測](../benchmarks/m1.md)）。
+  グループの並列化は将来の性能の項目です。ボラ面の校正の所要時間を測り、必要なら実行器を引数で受ける形で入れます。
 
 ## メモリ
 
diff --git a/docs/design/11-quality.md b/docs/design/11-quality.md
index 4cbfe65b10bb409243dfc961b94596802fce7f04..2c936c10eb342084308b8191ad9f66573bc6bfe1 100644
--- a/docs/design/11-quality.md
+++ b/docs/design/11-quality.md
@@ -35,7 +35,7 @@
 | 補助 | 役割 |
 |---|---|
 | `assert_close!(x, y, rel, abs)`、`assert_close!(x, y, rel, abs, scale)` | `abs(x − y) ≤ max(rel·M, abs)` を確かめる。M は `max(abs(x), abs(y))` で、5 引数の形は和の項の絶対値の和などの `scale` も M に含める（[許容誤差](#許容誤差)） |
-| `assert_mc!(estimate, se, reference, k, a)` | `abs(estimate − reference) ≤ k·se + a` を確かめる |
+| `assert_mc!(estimate, se, reference, k, a)` | `abs(estimate − reference) ≤ k·se + a` を確かめる。MC とともに M4 で足す |
 | `check_aad_vs_bump(reprice, aad, unit, pv_error, absolute, label)` | quote の増分（自然な単位）から再校正込みの価格を返す閉包、報告単位あたりの AAD の感応度、報告単位の幅 h₁、再校正の PV の誤差 δ_V、絶対の下限を受け取り、下の手順で照合する |
 | `compile_fail_pair!("case")` | `tests/ui/case_ok.rs` を実行し、`tests/ui/case.rs` のコンパイル診断を `case.json` と照合する |
 | `TickClock` | 読むたびに 1 ns 進む決定的な時計。`TickClock::after(n)` を期限にすると n 回目の読み出しで止まるので、区間の開始後の取消と期限の伝わり方をスレッドの速さによらず確かめられる |
@@ -89,10 +89,10 @@
 | 外貨担保曲線の長期（M1） | T+1 のペアの spot、ON と長期の XCCY だけの定義が当日の締めの前後で校正でき、XCCY の quote を再現する（`fx_implied_long_end_alone_survives_the_cutoff_of_a_t_plus_one_pair`）。EUR と JPY の XCCY の慣行を USDJPY の長期に書くと `InvalidData`（`fx_implied_long_end_rejects_an_xccy_convention_of_another_pair`） | [曲線と market](06-market/02-curves.md#fx-のフォワード) |
 | FX の慣行の正本（M1） | MtM XCCY のリセットの観測日と FX fixing の受渡日は参照データのペアのスポットラグで決まり、JPY の祝日を含む暦で T+1 と T+2 の日付と PV が手計算と一致する。USDCAD に TN を足した定義は T+1 の慣行では連鎖せず、慣行だけを T+2 にすると校正できる（`mtm_xccy_reset_fixing_dates_follow_the_pair_spot_lag`、`fx_fixing_spot_value_uses_the_pair_spot_lag`、`fx_x0_chain_follows_refdata_spot_lag`）。元本リセットも通常の fixing と同じくペアのスポットの受渡（東京 15 時）の F で予想し、新しい元本の交換（USD の決済時刻）とは分ける。明示の reset 日でスポットの受渡が交換日より早ければ `PaymentDelayIgnored` を記録する（`mtm_xccy_reset_reads_the_fixing_at_the_pair_spot_delivery`） | [参照データ](04-foundation/04-reference-data.md#fx-の慣行) |
 | X(0) の正本（M1） | ON・TN の quote があり、市場の担保のフォワードが `FxImplied` で ON・TN・スポットに校正されない定義は `InvalidData`（`fx_x0_on_tn_consistency`）。直接のクロスの quote と三角換算が食い違えば警告を出し、直接の quote を使う（`fx_cross_direct_vs_triangulated`） | [曲線と market](06-market/02-curves.md#fx-のフォワード) |
-| ボラの種類と式（M2） | 正規ボラの面を Black だけの式に、対数正規の面を Bachelier だけの式に解決すると `prepare` が `InvalidInput`（`vol_kind_formula_mismatch`） | [ボラティリティ](06-market/04-volatility.md#読み出し) |
+| ボラの種類と式（M2、[#42](https://github.com/kakune/rhodium/issues/42)） | 正規ボラの面を Black だけの式に、対数正規の面を Bachelier だけの式に解決すると `prepare` が `InvalidInput`（`vol_kind_formula_mismatch`） | [ボラティリティ](06-market/04-volatility.md#読み出し) |
 | 相対尺度のガンマ（M3） | V(q) = q のスポットガンマは 0。V(q) = q² なら相対 1% の二軸で 2·(0.01q)² | [感応度](08-risk/01-sensitivities.md#ガンマ) |
 | 上書きと下流の随伴（M3） | θ₁ = q、θ₂ = θ₁、V = θ₂ で θ₁ を上書きしても dV/dq = 1 | [感応度](08-risk/01-sensitivities.md#3-校正の随伴) |
-| smile dynamics の付け直し（M2） | 市場の残差が 0 でない過剰決定の SABR で、`StickyStrike` と `StickyMoneyness` の B が基準の面のモデルのボラへ付け直した有限差分と一致する（`smile_dynamics_overdetermined_refit`） | [感応度](08-risk/01-sensitivities.md#smile-dynamics) |
+| smile dynamics の付け直し（M2、[#43](https://github.com/kakune/rhodium/issues/43)） | 市場の残差が 0 でない過剰決定の SABR で、`StickyStrike` と `StickyMoneyness` の B が基準の面のモデルのボラへ付け直した有限差分と一致する（`smile_dynamics_overdetermined_refit`） | [感応度](08-risk/01-sensitivities.md#smile-dynamics) |
 | Philox の時刻追加（M4） | 同一 seed で 10 時の一点を求めた場合と 15 時も求めた場合の、10 時のモデル状態が許容誤差内で一致 | [モンテカルロ](07-pricing/04-monte-carlo.md#乱数) |
 | 日中の VM（M5） | コール 10 時と受渡 15 時の直前・ちょうど・直後で未決済と保有を照合し、同時刻の別周期のコールと cutoff が決済後の保有を読む | [エクスポージャーと担保](09-xva/01-exposure-and-collateral.md#イベントの順序) |
 | 部分的な校正（M1） | USD の quote を除くと USD を読む取引だけが失敗する。補って再校正すれば保存済みの準備を再び使え、全体を解き直した価格・デルタと一致する（`portfolio_isolates_trade_failures_and_totals_successes`） | [校正](06-market/03-calibration.md#失敗) |
