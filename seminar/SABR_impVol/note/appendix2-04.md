diff --git a/crates/rhodium-pricing/tests/ois_calibration.rs b/crates/rhodium-pricing/tests/ois_calibration.rs
index 884469ab5d34c2707441004ef7fd93260147cbc1..e43b1910acd5b816b10ae143f0631867294713ee 100644
--- a/crates/rhodium-pricing/tests/ois_calibration.rs
+++ b/crates/rhodium-pricing/tests/ois_calibration.rs
@@ -35,6 +35,10 @@ use scioxide_core::AdjointOperator;
 #[path = "fixtures/ois.rs"]
 pub mod ois;
 
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、このテストで使わない関数がある"
+)]
 #[path = "fixtures/multicurve.rs"]
 pub mod multicurve;
 
@@ -134,7 +138,8 @@ fn sofr_ois_calibration_matches_quantlib() -> Result<()> {
 #[test]
 fn multicurve_calibration_matches_quantlib() -> Result<()> {
     // OIS を先に解き、OIS で割り引くベーシス商品でターム指標の予想を解く。
-    for case in multicurve::multicurve_cases()? {
+    // JPY の TIBOR 6M と 3M、USD のテスト用ターム指標（6M）の三つを QuantLib と照合する。
+    for case in multicurve::all_multicurve_cases()? {
         let inputs = case.inputs()?;
         let ctx = inputs.context();
         let result = calibrate(&case.definition()?, &case.quotes()?, &ctx)?;
@@ -911,6 +916,10 @@ fn short_ois_calibrate_to_the_default_tolerance() -> Result<()> {
 }
 
 #[test]
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#196 で各 variant を網羅する"
+)]
 fn rate_quotes_above_one_warn_about_a_missing_unit() -> Result<()> {
     // 金利の binding が読む quote の絶対値が 1 を超えたら、校正の結果に quote ID 付きで警告する。
     // 5% を単位なしの 5 と書いたような誤りは、値だけでは 500% として校正されてしまうため。
diff --git a/crates/rhodium-pricing/tests/par.rs b/crates/rhodium-pricing/tests/par.rs
index e85c0a930ec52f9faa12171b8e885dfa1d8f951d..fa4ceb80ab907f6cfa764db1a8c2152aac8df578 100644
--- a/crates/rhodium-pricing/tests/par.rs
+++ b/crates/rhodium-pricing/tests/par.rs
@@ -63,3 +63,122 @@ fn par_rejects_unbracketed_or_singular_roots() -> Result<()> {
     ));
     Ok(())
 }
+
+/// quote の値 p について affine でない PV: (1 + p)^n − a·exp(b·θ)。ZC インフレスワップの
+/// 固定脚 N·((1 + K)^n − 1) と同じ形で、パーは閉形式 p* = (a·e^{bθ})^{1/n} − 1 を持つ。
+struct ZeroCoupon {
+    years: f64,
+    scale: f64,
+    slope: f64,
+}
+
+impl Residual for ZeroCoupon {
+    fn residual<R: Real<f64, Mask = bool>>(
+        &self,
+        rate: &[R],
+        parameters: &[R],
+        out: &mut [R],
+    ) -> std::result::Result<(), SciOxideError> {
+        out[0] = (rate[0] + R::constant(1.0)).powf(R::constant(self.years))
+            - (parameters[0] * R::constant(self.slope)).exp() * R::constant(self.scale);
+        Ok(())
+    }
+}
+
+impl ZeroCoupon {
+    /// 閉形式のパーと、θ についての 1 階・2 階の微分。
+    fn exact(&self, theta: f64) -> [f64; 3] {
+        let growth = (self.scale * (self.slope * theta).exp()).powf(self.years.recip());
+        let rate = self.slope / self.years;
+        [growth - 1.0, rate * growth, rate * rate * growth]
+    }
+}
+
+/// 二つの非線形のパーを重み 1 で当てる最小二乗の停留点 θ*(q) を、閉形式のパーのニュートン法で解く。
+fn reoptimize(instruments: &[ZeroCoupon; 2], quotes: [f64; 2]) -> f64 {
+    let mut theta = 0.0;
+    for _ in 0..100 {
+        let (gradient, hessian) = instruments.iter().zip(quotes).fold(
+            (0.0, 0.0),
+            |(gradient, hessian), (instrument, quote)| {
+                let [m, first, second] = instrument.exact(theta);
+                (
+                    gradient + (m - quote) * first,
+                    hessian + first * first + (m - quote) * second,
+                )
+            },
+        );
+        let step = gradient / hessian;
+        theta -= step;
+        if step.abs() < 1e-17 {
+            break;
+        }
+    }
+    theta
+}
+
+#[test]
+fn par_nonlinear_lsq_quote_delta() -> Result<()> {
+    // 残差の残る最小二乗 min ½Σ(m_i(θ) − q_i)² のモデルの quote m_i が非線形のパー。停留条件
+    // G = Σ(m_i − q_i)·m_i′ の θ 微分 H = Σ(m_i′² + (m_i − q_i)·m_i″) はパーの 2 階の微分を含むので、
+    // quote デルタ dθ/dq_j = m_j′/H を正しく求めるには m_i″ まで正しい必要がある。
+    // 陰関数定理を入れ子の AD の次数まで当てる solve_nonlinear は再最適化の有限差分と一致し、
+    // 1 回のニュートン補正で R の値を作る近道（2 階の微分を誤る）は一致しない。
+    let instruments = [
+        ZeroCoupon {
+            years: 2.0,
+            scale: 1.02,
+            slope: 1.0,
+        },
+        ZeroCoupon {
+            years: 5.0,
+            scale: 1.1,
+            slope: 3.0,
+        },
+    ];
+    let quotes = [0.03, 0.0];
+    let theta = reoptimize(&instruments, quotes);
+    let seeded = Dual::new(Dual::new(theta, 1.0), Dual::lift(1.0));
+    let mut implicit = Vec::new();
+    let mut shortcut = Vec::new();
+    for (instrument, quote) in instruments.iter().zip(quotes) {
+        let [m, first, second] = instrument.exact(theta);
+        // 停留点でも残差は 0 にならない（当てはめの残差が 2 階の項に効く）。
+        assert!((m - quote).abs() > 1e-3);
+        let par = solve_nonlinear(instrument, &[seeded], initial(-0.5, 1.0)?)?;
+        assert_close!(par.value().value(), m, 1e-12, 1e-14);
+        assert_close!(par.derivative().value(), first, 1e-10, 1e-12);
+        assert_close!(par.derivative().derivative(), second, 1e-8, 1e-12);
+        implicit.push([m, first, par.derivative().derivative()]);
+        // 近道: p = p₀ − PV(p₀, θ)/∂PV/∂p(p₀, θ)。1 階の微分は正しいが、2 階の微分は n 倍になる。
+        let mut pv = [Dual::lift(Dual::lift(0.0))];
+        instrument
+            .residual(&[Dual::lift(Dual::lift(m))], &[seeded], &mut pv)
+            .map_err(|error| Error::Numerical(error.to_string()))?;
+        let slope = instrument.years * (1.0 + m).powf(instrument.years - 1.0);
+        let approximate = Dual::lift(Dual::lift(m)) - pv[0] * Dual::lift(Dual::lift(slope.recip()));
+        assert_close!(approximate.derivative().value(), first, 1e-10, 1e-12);
+        shortcut.push([m, first, approximate.derivative().derivative()]);
+    }
+    let delta = |parts: &[[f64; 3]], j: usize| {
+        let hessian: f64 = parts
+            .iter()
+            .zip(quotes)
+            .map(|(&[m, first, second], quote)| first * first + (m - quote) * second)
+            .sum();
+        parts[j][1] / hessian
+    };
+    let h = 1e-5;
+    for j in 0..2 {
+        let bumped = |shift: f64| {
+            let mut quotes = quotes;
+            quotes[j] += shift;
+            reoptimize(&instruments, quotes)
+        };
+        let finite_difference = (bumped(h) - bumped(-h)) / (2.0 * h);
+        assert_close!(delta(&implicit, j), finite_difference, 1e-7, 1e-12);
+        let error = (delta(&shortcut, j) - finite_difference).abs() / finite_difference.abs();
+        assert!(error > 1e-4, "shortcut quote delta error {error:e}");
+    }
+    Ok(())
+}
diff --git a/crates/rhodium-pricing/tests/volatility_calibration.rs b/crates/rhodium-pricing/tests/volatility_calibration.rs
new file mode 100644
index 0000000000000000000000000000000000000000..e31ab206c5bb852414d0b3f9067aa8a9791a693e
--- /dev/null
+++ b/crates/rhodium-pricing/tests/volatility_calibration.rs
@@ -0,0 +1,626 @@
+//! M2-1 と M2-2 の金利ボラ校正を独立の価格例と照合する。
+
+use rhodium_core::{Error, Result, Warning, quantity::VolKind, time::model_time};
+use rhodium_market::{
+    curve::{CurveExpr, CurveRole, CurveSpec},
+    quote::{QuoteId, QuoteSet},
+    vol::{CubeInterpolation, Sabr, SabrFit, VolSpec},
+};
+use rhodium_pricing::{
+    analytic::rfr_caplet_variance_time,
+    calibration::{CalibratedMarket, CalibrationGroup, GroupOutcome, calibrate},
+    definition::{
+        MarketDefinition,
+        vol::{
+            CapVolQuote, RateStrike, RateVolKind, RfrVariance, SabrCalibration, SwaptionVolQuote,
+            VolatilityDefinition, VolatilityMarket,
+        },
+    },
+    resolve::CollateralKey,
+};
+use rhodium_testkit::assert_close;
+use serde::Deserialize;
+
+#[allow(
+    dead_code,
+    reason = "共有 fixture の曲線監査用のヘルパーはここでは使わない"
+)]
+#[path = "fixtures/ois.rs"]
+mod ois;
+
+#[allow(dead_code, reason = "共有 fixture の慣行と参照データだけを使う")]
+#[path = "fixtures/multicurve.rs"]
+mod multicurve;
+
+const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
+
+#[derive(Deserialize)]
+struct Oracle {
+    sabr: Vec<SabrCase>,
+    caps: Vec<CapCase>,
+}
+
+#[derive(Deserialize)]
+struct SabrCase {
+    kind: String,
+    forward: f64,
+    strike: f64,
+    shift: f64,
+    time: f64,
+    beta: f64,
+    alpha: f64,
+    nu: f64,
+    rho: f64,
+    vol: f64,
+}
+
+#[derive(Deserialize)]
+struct CapCase {
+    kind: String,
+    strike: f64,
+    nodes: Vec<CapNode>,
+    quotes: Vec<CapQuote>,
+}
+
+#[derive(Deserialize)]
+struct CapNode {
+    expiry: rhodium_core::date::Date,
+    annuity: f64,
+    time: f64,
+    vol: f64,
+}
+
+#[derive(Deserialize)]
+struct CapQuote {
+    maturity: rhodium_core::tenor::TenorLength,
+    vol: f64,
+    pv: f64,
+}
+
+fn fixture() -> Result<Oracle> {
+    serde_json::from_str(include_str!("fixtures/quantlib/volatility.json"))
+        .map_err(|error| Error::InvalidData(error.to_string()))
+}
+
+#[test]
+fn hagan_sabr_matches_quantlib() -> Result<()> {
+    // シフトなし・ありの Hagan 式を ATM の両側と離れたストライクで直接比較する。
+    for case in fixture()?.sabr {
+        let kind = match case.kind.as_str() {
+            "normal" | "shifted_normal" => VolKind::Normal,
+            "lognormal" => VolKind::Lognormal,
+            _ => VolKind::ShiftedLognormal { shift: case.shift },
+        };
+        let vol = Sabr::new(kind, case.beta, case.shift)?.volatility(
+            case.forward,
+            case.strike,
+            case.time,
+            [case.alpha, case.rho, case.nu],
+        )?;
+        assert_close!(vol, case.vol, 1e-12, 1e-14);
+    }
+    Ok(())
+}
+
+#[test]
+fn rfr_caplet_variance_matches_hand_calculation() -> Result<()> {
+    // 30 日の期間。開始の 5 日前、開始時、10 日経過時、終了時、終了後を比較する。
+    let start = "2025-01-06T00:00:00Z".parse()?;
+    let end = "2025-02-05T00:00:00Z".parse()?;
+    for (as_of, expected_days) in [
+        ("2025-01-01T00:00:00Z", 15.0),
+        ("2025-01-06T00:00:00Z", 10.0),
+        (
+            "2025-01-16T00:00:00Z",
+            20.0 * 20.0 * 20.0 / (3.0 * 30.0 * 30.0),
+        ),
+        ("2025-02-05T00:00:00Z", 0.0),
+        ("2025-02-06T00:00:00Z", 0.0),
+    ] {
+        assert_close!(
+            rfr_caplet_variance_time(as_of.parse()?, start, end)?,
+            expected_days / 365.0,
+            1e-12,
+            1e-14
+        );
+    }
+    assert!(rfr_caplet_variance_time(start, end, start).is_err());
+    Ok(())
+}
+
+fn surface(id: &str, index: &str, convention: &str) -> Result<VolatilityMarket> {
+    Ok(VolatilityMarket {
+        id: id.parse()?,
+        index: index.parse()?,
+        convention: convention.parse()?,
+        collateral: COLLATERAL,
+        vol_kind: RateVolKind::Normal,
+        shift: None,
+    })
+}
+
+#[test]
+fn caplet_stripping_reprices_term_and_rfr_caps() -> Result<()> {
+    // 予想金利ゼロ、割引係数 1 の例で、異なる区間ボラをフラットボラから復元する。
+    let [case, _] = multicurve::multicurve_cases()?;
+    let inputs = case.inputs()?;
+    let curves = [&case.overnight.index, &case.term_index]
+        .into_iter()
+        .map(|id| {
+            CurveSpec::from_expr(
+                id.parse()?,
+                case.overnight.currency,
+                [CurveRole::Discount, CurveRole::Projection],
+                CurveExpr::Flat(0.0),
+            )
+        })
+        .collect::<Result<Vec<_>>>()?;
+    let cases = fixture()?.caps;
+    let mut marks = Vec::new();
+    let mut surfaces = Vec::new();
+    for (kind, index, convention) in [
+        ("term", case.term_index.clone(), case.irs_convention()),
+        (
+            "overnight",
+            case.overnight.index.clone(),
+            case.overnight.index.clone(),
+        ),
+    ] {
+        let mut quotes = Vec::new();
+        for cap in cases.iter().filter(|cap| cap.kind == kind) {
+            for (i, row) in cap.quotes.iter().enumerate() {
+                let coordinate = if cap.strike == 0.0 { "ATM" } else { "K" };
+                let id: QuoteId =
+                    format!("CAP.{}.{i}.{coordinate}", kind.to_ascii_uppercase()).parse()?;
+                marks.push((id.clone(), row.vol));
+                let strike = if cap.strike == 0.0 {
+                    RateStrike::Atm
+                } else {
+                    RateStrike::Absolute(cap.strike)
+                };
+                quotes.push(CapVolQuote {
+                    quote: id,
+                    maturity: row.maturity,
+                    strike,
+                });
+            }
+        }
+        surfaces.push(VolatilityDefinition::Caplet {
+            surface: surface(kind, &index, &convention)?,
+            rfr_variance: (kind == "overnight").then_some(RfrVariance::LyashenkoMercurio),
+            quotes,
+        });
+    }
+    let definition = MarketDefinition::new("JPY".parse()?, "Etc/UTC".parse()?, curves, [])?
+        .with_volatilities(surfaces)?;
+    let calibrated = calibrate(&definition, &QuoteSet::new(marks)?, &inputs.context())?;
+    check_solved(&calibrated)?;
+    let market = calibrated.market();
+    for cap in &cases {
+        let id = cap.kind.parse()?;
+        let mut pv = 0.0;
+        for (node, cap_quote) in cap.nodes.iter().zip(&cap.quotes) {
+            let point =
+                market
+                    .structure()
+                    .prepare_vol_point(&id, ois::midnight(node.expiry)?, None)?;
+            let vol = market.with_view(|view| {
+                view.vol_at(view.bind_vol(view.volatility(&id)?, &point)?, cap.strike)
+            })?;
+            assert_close!(vol, node.vol, 1e-12, 1e-14);
+            // 独立の Bachelier 式で再現した価格を QuantLib のキャップ価格と比べる。
+            let sigma = vol * node.time.sqrt();
+            let d = -cap.strike / sigma;
+            let phi = (-d * d / 2.0).exp() / (2.0 * std::f64::consts::PI).sqrt();
+            let cdf = scioxide_special::normal_cdf(d)
+                .map_err(|error| Error::Numerical(error.to_string()))?;
+            pv += node.annuity * (sigma * phi - cap.strike * cdf);
+            assert_close!(pv, cap_quote.pv, 1e-12, 1e-14);
+            market.with_view(|view| {
+                let other = if cap.kind == "term" {
+                    "overnight"
+                } else {
+                    "term"
+                };
+                assert!(matches!(
+                    view.bind_vol(view.volatility(&other.parse()?)?, &point),
+                    Err(Error::InvalidInput(_))
+                ));
+                Ok::<_, Error>(())
+            })?;
+        }
+    }
+    // 保存したストリッピングの線形化（逆算ボラの陰関数微分）を再校正の差分と照合する。
+    // 先頭のキャップの quote は、前の区間の価格を固定して解く 2 本目の区間ボラにも効くので、
+    // 逐次の依存も確かめられる。絶対ストライク列と ATM 列、ターム金利と後決め RFR を選ぶ。
+    for (surface, id) in [
+        ("term", "CAP.TERM.0.K"),
+        ("term", "CAP.TERM.0.ATM"),
+        ("overnight", "CAP.OVERNIGHT.0.K"),
+        ("overnight", "CAP.OVERNIGHT.0.ATM"),
+    ] {
+        check_quote_sensitivity(&calibrated, surface, &id.parse()?)?;
+    }
+    check_impossible_cap(&calibrated)?;
+    Ok(())
+}
+
+fn check_impossible_cap(calibrated: &CalibratedMarket) -> Result<()> {
+    // 2 本目のキャップが先にストリップ済みの期より安くなる quote は解けない。
+    // この暗黙の裁定が一面にあっても、独立の後決め RFR の面には影響しない。
+    let id: QuoteId = "CAP.TERM.1.K".parse()?;
+    let old = calibrated
+        .quotes()
+        .get(&id)
+        .ok_or_else(|| Error::Internal("missing test cap quote".into()))?;
+    let changed = calibrated.update(&bumped(calibrated.quotes(), &id, 1e-6 - old)?)?;
+    for (base, next) in calibrated.groups().iter().zip(changed.groups()) {
+        if next.blocks()[0].curve().as_str() == "term" {
+            assert!(matches!(
+                next.outcome(),
+                GroupOutcome::Failed(Error::Calibration(_))
+            ));
+        } else {
+            assert!(std::ptr::eq(base.outcome(), next.outcome()));
+        }
+    }
+    Ok(())
+}
+
+fn check_solved(calibrated: &CalibratedMarket) -> Result<()> {
+    for group in calibrated.groups() {
+        if let GroupOutcome::Failed(error) = group.outcome() {
+            return Err(error.clone());
+        }
+        assert!(
+            matches!(group.outcome(), GroupOutcome::Solved(_)),
+            "{group:?}"
+        );
+    }
+    Ok(())
+}
+
+fn sabr_definition(case: &ois::Case, fit: SabrFit) -> Result<VolatilityDefinition> {
+    let mut quotes = Vec::new();
+    for expiry in ["3M", "6M"] {
+        for tenor in ["6M", "1Y"] {
+            let atm: QuoteId = format!("SABR.{expiry}.{tenor}.ATM").parse()?;
+            quotes.push(SwaptionVolQuote {
+                quote: atm.clone(),
+                expiry: expiry.parse()?,
+                tenor: tenor.parse()?,
+                strike: RateStrike::Atm,
+                atm_quote: None,
+                tolerance: matches!(fit, SabrFit::LeastSquares).then_some(1e-5),
+            });
+            if !matches!(fit, SabrFit::AtmOnly { .. }) {
+                for (j, offset) in [-0.002, -0.001, 0.001, 0.003].into_iter().enumerate() {
+                    quotes.push(SwaptionVolQuote {
+                        quote: format!("SABR.{expiry}.{tenor}.W{j}").parse()?,
+                        expiry: expiry.parse()?,
+                        tenor: tenor.parse()?,
+                        strike: RateStrike::AtmOffset(offset),
+                        atm_quote: Some(atm.clone()),
+                        tolerance: Some(1e-5),
+                    });
+                }
+            }
+        }
+    }
+    Ok(VolatilityDefinition::SwaptionCube {
+        surface: surface("JPY-CUBE", &case.index, &case.index)?,
+        sabr: SabrCalibration {
+            beta: 0.5,
+            shift: 0.01,
+            fit,
+            interpolation: CubeInterpolation::TotalVariance,
+        },
+        expiry_time: rhodium_core::time::LocalEventTime {
+            time: "00:00:00".parse()?,
+            zone: "Etc/UTC".parse()?,
+        },
+        quotes,
+    })
+}
+
+// ボラのブロックは面の ID を曲線の位置に持つので、面ごとのグループをその ID で選ぶ。
+fn surface_groups<'a>(
+    market: &'a CalibratedMarket,
+    surface: &'a str,
+) -> impl Iterator<Item = &'a CalibrationGroup> {
+    market
+        .groups()
+        .iter()
+        .filter(move |g| g.blocks()[0].curve().as_str() == surface)
+}
+
+fn cube_groups(market: &CalibratedMarket) -> impl Iterator<Item = &CalibrationGroup> {
+    surface_groups(market, "JPY-CUBE")
+}
+
+fn bumped(quotes: &QuoteSet, id: &QuoteId, bump: f64) -> Result<QuoteSet> {
+    QuoteSet::new(
+        quotes
+            .iter()
+            .map(|(q, value)| (q.clone(), if q == id { value + bump } else { value })),
+    )
+}
+
+#[test]
+fn swaption_cube_reprices_market_quotes() -> Result<()> {
+    // 先に ATM-only で各ノードの曲線由来の F を読む。次に既知の SABR で各ウィングを作る。
+    let [case, _] = ois::cases()?;
+    let inputs = ois::Inputs::new(&[&case], case.as_of, "")?;
+    let ctx = inputs.context();
+    let base = ois::definition(&[(&case, &COLLATERAL)])?;
+    let rates = ois::fixture_quotes(&case)?;
+    let atm = sabr_definition(
+        &case,
+        SabrFit::AtmOnly {
+            rho: -0.25,
+            nu: 0.4,
+        },
+    )?;
+    let VolatilityDefinition::SwaptionCube { quotes, .. } = &atm else {
+        return Err(Error::Internal("expected cube".into()));
+    };
+    let marks = QuoteSet::new(
+        rates
+            .iter()
+            .cloned()
+            .chain(quotes.iter().map(|q| (q.quote.clone(), 0.002))),
+    )?;
+    let first = calibrate(&base.clone().with_volatilities([atm])?, &marks, &ctx)?;
+    check_solved(&first)?;
+    for fit in first.volatility_fit() {
+        assert_close!(fit.model, 0.002, 1e-12, 1e-14);
+    }
+    let sabr = Sabr::new(VolKind::Normal, 0.5, 0.01)?;
+    let upstream_quote = rates[2].0.clone();
+    let mut quotes = rates;
+    for (e, (expiry, date)) in [("3M", "2025-04-07"), ("6M", "2025-07-07")]
+        .into_iter()
+        .enumerate()
+    {
+        let time = model_time(ctx.as_of, ois::midnight(date.parse()?)?).years();
+        for (t, tenor) in ["6M", "1Y"].into_iter().enumerate() {
+            let block = VolSpec::swaption_block(&"JPY-CUBE".parse()?, e, t)?;
+            let range = first
+                .market()
+                .structure()
+                .layout()
+                .range(&block)
+                .ok_or_else(|| Error::Internal("missing cube block".into()))?;
+            let f = first.market().theta()[range.start + 1];
+            let params = [
+                0.019 + 0.002 * f64::from(u32::try_from(e + t).unwrap_or_default()),
+                -0.25,
+                0.4,
+            ];
+            let atm = sabr.volatility(f, f, time, params)?;
+            quotes.push((format!("SABR.{expiry}.{tenor}.ATM").parse()?, atm));
+            for (j, offset) in [-0.002, -0.001, 0.001, 0.003].into_iter().enumerate() {
+                let difference = sabr.volatility(f, f + offset, time, params)? - atm;
+                quotes.push((format!("SABR.{expiry}.{tenor}.W{j}").parse()?, difference));
+            }
+        }
+    }
+    let quotes = QuoteSet::new(quotes)?;
+    for fit in [SabrFit::LeastSquares, SabrFit::ExactAtm] {
+        let definition = base
+            .clone()
+            .with_volatilities([sabr_definition(&case, fit)?])?;
+        let result = calibrate(&definition, &quotes, &ctx)?;
+        check_solved(&result)?;
+        assert_eq!(cube_groups(&result).count(), 4);
+        for row in result.volatility_fit() {
+            assert_close!(row.model, row.target, 1e-12, 1e-14);
+        }
+        // ATM quote は絶対 ATM と各差分行に効く。他のノードの解は再利用される。
+        check_cube_surface_reads(result.market())?;
+        check_quote_sensitivity(&result, "JPY-CUBE", &"SABR.3M.6M.ATM".parse()?)?;
+        check_cube_curve_sensitivity(&result, &quotes, &upstream_quote)?;
+        let outlier = bumped(&quotes, &"SABR.3M.6M.W0".parse()?, 0.0004)?;
+        let perturbed = result.update(&outlier)?;
+        check_solved(&perturbed)?;
+        assert!(
+            perturbed
+                .warnings()
+                .iter()
+                .any(|w| matches!(w, Warning::CalibrationFitOutsideTolerance { .. }))
+        );
+        if matches!(fit, SabrFit::ExactAtm) {
+            for row in perturbed
+                .volatility_fit()
+                .filter(|row| row.quote.as_str().rsplit('.').next() == Some("ATM"))
+            {
+                assert_close!(row.model, row.target, 1e-12, 1e-14);
+            }
+        }
+    }
+    Ok(())
+}
+
+#[test]
+fn fixed_sabr_parameters_are_checked_when_defined() -> Result<()> {
+    // ATM だけを当てるときに固定する ρ・ν は、校正を待たずに定義の読み込みで拒否する。
+    // 定義域の端（|ρ| = 1、負の ν、非有限）を拒み、内側の値は受け付ける。
+    let [case, _] = ois::cases()?;
+    let base = ois::definition(&[(&case, &COLLATERAL)])?;
+    for (rho, nu) in [
+        (1.0, 0.4),
+        (-1.5, 0.4),
+        (0.2, -0.1),
+        (f64::NAN, 0.4),
+        (0.2, f64::INFINITY),
+    ] {
+        let definition = sabr_definition(&case, SabrFit::AtmOnly { rho, nu })?;
+        assert!(
+            matches!(
+                base.clone().with_volatilities([definition]),
+                Err(Error::InvalidInput(_))
+            ),
+            "rho {rho}, nu {nu}"
+        );
+    }
+    let definition = sabr_definition(
+        &case,
+        SabrFit::AtmOnly {
+            rho: -0.99,
+            nu: 0.0,
+        },
+    )?;
+    base.with_volatilities([definition])?;
+    Ok(())
+}
+
+fn check_cube_surface_reads(market: &rhodium_market::market::Market) -> Result<()> {
+    // 同じ ATM オフセットで、実際に面を読む経路が節点の値と総分散を再現する。
+    let structure = market.structure();
+    let id = "JPY-CUBE".parse()?;
+    let times =
+        ["2025-04-07", "2025-05-07", "2025-07-07"].map(|date| date.parse().and_then(ois::midnight));
+    let [first, middle, last] = times;
+    let (first, middle, last) = (first?, middle?, last?);
+    let tenor = model_time(structure.as_of(), last).years();
+    let offsets = [-0.002, 0.0, 0.003];
+    let tenor_index = 0;
+    let t0 = model_time(structure.as_of(), first).years();
+    let t1 = model_time(structure.as_of(), last).years();
+    let t = model_time(structure.as_of(), middle).years();
+    let weight = (t - t0) / (t1 - t0);
+    let sabr = Sabr::new(VolKind::Normal, 0.5, 0.01)?;
+    market.with_view(|view| {
+        let mut forward = [0.0; 2];
+        let mut node_vols = [[0.0; 3]; 2];
+        let face = view.volatility(&id)?;
+        for (e, at) in [first, last].into_iter().enumerate() {
+            let block = VolSpec::swaption_block(&id, e, tenor_index)?;
+            let range = structure
+                .layout()
+                .range(&block)
+                .ok_or_else(|| Error::Internal("missing cube block".into()))?;
+            let f = market.theta()[range.start + 3];
+            forward[e] = f;
+            let point = structure.prepare_vol_point(&id, at, Some(tenor))?;
+            let bound = view.bind_vol(face, &point)?;
+            let alpha = if e == 0 { 0.019 } else { 0.021 };
+            for (j, offset) in offsets.into_iter().enumerate() {
+                let want = sabr.volatility(
+                    f,
+                    f + offset,
+                    model_time(structure.as_of(), at).years(),
+                    [alpha, -0.25, 0.4],
+                )?;
+                let vol = view.vol_at(bound, f + offset)?;
+                assert_close!(vol, want, 1e-12, 1e-14);
+                node_vols[e][j] = vol;
+            }
+        }
+        let point = structure.prepare_vol_point(&id, middle, Some(tenor))?;
+        let bound = view.bind_vol(face, &point)?;
+        let f = forward[0] * (1.0 - weight) + forward[1] * weight;
+        for (j, offset) in offsets.into_iter().enumerate() {
+            let variance = node_vols[0][j].powi(2) * t0 * (1.0 - weight)
+                + node_vols[1][j].powi(2) * t1 * weight;
+            assert_close!(
+                view.vol_at(bound, f + offset)?,
+                (variance / t).sqrt(),
+                1e-12,
+                1e-14
+            );
+        }
+        Ok(())
+    })
+}
+
+// 保存した解の線形化から作る θ の接ベクトル −G_θ⁻¹·G_q を、quote を上下に動かして
+// 再校正した θ の中心差分と比べる。quote を読まないグループは解き直さず、前の解を再利用する。
+fn check_quote_sensitivity(market: &CalibratedMarket, surface: &str, id: &QuoteId) -> Result<()> {
+    let quotes = market.quotes();
+    let up = market.update(&bumped(quotes, id, 1e-7)?)?;
+    let down = market.update(&bumped(quotes, id, -1e-7)?)?;
+    check_solved(&up)?;
+    check_solved(&down)?;
+    let mut read = false;
+    for (base, changed) in surface_groups(market, surface).zip(surface_groups(&up, surface)) {
+        if let Ok(column) = base.quotes().binary_search(id) {
+            read = true;
+            let GroupOutcome::Solved(solution) = base.outcome() else {
+                return Err(Error::Internal("expected solved volatility group".into()));
+            };
+            let range = market
+                .market()
+                .structure()
+                .layout()
+                .range(&base.blocks()[0])
+                .ok_or_else(|| Error::Internal("missing volatility parameters".into()))?;
+            let rhs = (0..range.len())
+                .map(|row| solution.dg_dq().get(row, column).map(|value| -value))
+                .collect::<std::result::Result<Vec<_>, _>>()
+                .map_err(|error| Error::Numerical(error.to_string()))?;
+            let tangent = solution
+                .dg_dtheta()
+                .solve(rhs.as_slice())
+                .map_err(|error| Error::Numerical(error.to_string()))?;
+            for (i, &expected) in range.zip(tangent.as_slice()) {
+                let numerical = (up.market().theta()[i] - down.market().theta()[i]) / 2e-7;
+                assert_close!(numerical, expected, 5e-6, 5e-7);
+            }
+        } else {
+            assert!(std::ptr::eq(base.outcome(), changed.outcome()));
+        }
+    }
+    assert!(read, "{surface} のどのグループも {id} を読まない");
+    Ok(())
+}
+
+// 上流曲線の変化はフォワードの派生値にも伝わる。保存したヤコビアンと再校正を比較する。
+fn check_cube_curve_sensitivity(
+    market: &CalibratedMarket,
+    quotes: &QuoteSet,
+    id: &QuoteId,
+) -> Result<()> {
+    let bump = 2e-8;
+    let up = market.update(&bumped(quotes, id, bump)?)?;
+    let down = market.update(&bumped(quotes, id, -bump)?)?;
+    check_solved(&up)?;
+    check_solved(&down)?;
+    let layout = market.market().structure().layout();
+    let numerical = |i| (up.market().theta()[i] - down.market().theta()[i]) / (2.0 * bump);
+    let mut moved = false;
+    for (base, changed) in cube_groups(market).zip(cube_groups(&up)) {
+        assert!(!std::ptr::eq(base.outcome(), changed.outcome()));
+        let GroupOutcome::Solved(solution) = base.outcome() else {
+            return Err(Error::Internal("expected solved SABR".into()));
+        };
+        let mut rates = Vec::new();
+        for block in base.upstream() {
+            let range = layout
+                .range(block)
+                .ok_or_else(|| Error::Internal("missing upstream curve parameters".into()))?;
+            rates.extend(range.map(numerical));
+        }
+        let range = layout
+            .range(&base.blocks()[0])
+            .ok_or_else(|| Error::Internal("missing SABR parameters".into()))?;
+        let rhs = (0..range.len())
+            .map(|row| {
+                rates.iter().enumerate().try_fold(0.0, |sum, (col, rate)| {
+                    Ok(sum - solution.dg_dup().get(row, col)? * rate)
+                })
+            })
+            .collect::<std::result::Result<Vec<_>, scioxide_core::SciOxideError>>()
+            .map_err(|error| Error::Numerical(error.to_string()))?;
+        let tangent = solution
+            .dg_dtheta()
+            .solve(rhs.as_slice())
+            .map_err(|error| Error::Numerical(error.to_string()))?;
+        moved |= tangent.as_slice().last().is_some_and(|x| x.abs() > 1e-10);
+        for (i, &expected) in range.zip(tangent.as_slice()) {
+            assert_close!(numerical(i), expected, 5e-6, 5e-7);
+        }
+    }
+    assert!(moved, "少なくとも一つのノードでフォワードが動く");
+    Ok(())
+}
diff --git a/crates/rhodium-pricing/tests/xccy_curves.rs b/crates/rhodium-pricing/tests/xccy_curves.rs
new file mode 100644
index 0000000000000000000000000000000000000000..9bec336156f2731b0444ad745dde1bff01597954
--- /dev/null
+++ b/crates/rhodium-pricing/tests/xccy_curves.rs
@@ -0,0 +1,379 @@
+//! 外貨担保曲線の長期を `MtM` XCCY のベーシスで校正し、短期の FX ポイントと一つのグループで解く（R6）。
+
+use std::ptr;
+
+use rhodium_core::{
+    Error, Result,
+    calendar::CalendarSpec,
+    currency::Currency,
+    fx::{CcyPair, FxRate},
+    id::IndexId,
+    money::Money,
+    refdata::{
+        Convention, DailyPublication, FxConvention, FxSpotLag, IndexDefinition, IndexKind,
+        IndexStatus, PublicationRule, SwapConvention, SwapCouponConvention, SwapLegConvention,
+    },
+    time::LocalEventTime,
+};
+use rhodium_market::quote::{QuoteId, QuoteSet};
+use rhodium_pricing::{
+    calibration::{CalibratedMarket, GroupOutcome, calibrate},
+    context::Context,
+    definition::FxImpliedLongEnd,
+    prepare::{PrepareInput, evaluate, prepare},
+};
+use rhodium_products::standard::{FxForward, FxSwap, Position, StandardProduct};
+use rhodium_testkit::assert_close;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、XCCY の校正に使わない関数がある"
+)]
+#[path = "fixtures/ois.rs"]
+pub mod ois;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、XCCY の校正に使わない関数がある"
+)]
+#[path = "fixtures/xccy.rs"]
+pub mod xccy;
+
+use ois::{Inputs, midnight};
+use xccy::{POINTS, SPOT, XCCY_QUOTES, XccyMarket, jpy_cash};
+
+const TOLERANCE: f64 = 1e-14;
+
+/// 評価日（2025-01-07、評価時点は UTC の 0 時）。校正商品の約定日と同じ。
+fn trade_date() -> Result<rhodium_core::date::Date> {
+    "2025-01-07".parse()
+}
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
+/// JPY 建ての PV。USD の受渡を評価時点の X(0) で換算する。
+fn jpy_pv(product: &StandardProduct, market: &CalibratedMarket, ctx: &Context<'_>) -> Result<f64> {
+    let plan = prepare(
+        PrepareInput::Standard {
+            product,
+            trade_date: trade_date()?,
+        },
+        &jpy_cash()?,
+        &[],
+        market.market().structure(),
+        ctx,
+    )?;
+    let pv = evaluate(&plan, market.market())?;
+    let x0 = market
+        .market()
+        .structure()
+        .fx_today_plan(XccyMarket::pair()?)?
+        .rate(market.market())?
+        .rate();
+    Ok(pv.get(Currency::JPY) + pv.get(Currency::USD) * x0)
+}
+
+/// `MtM` XCCY のテンプレートのパー（JPY の TONA 脚のベーシス）。PV は K について affine。
+fn xccy_par(tenor: &str, market: &CalibratedMarket, ctx: &Context<'_>) -> Result<f64> {
+    let template = XccyMarket::template(tenor)?;
+    let pv = |spread| {
+        let product = template.build_for_trade(spread, trade_date()?, ctx.refdata)?;
+        jpy_pv(&product, market, ctx)
+    };
+    let zero = pv(0.0)?;
+    Ok(-zero / (pv(1.0)? - zero))
+}
+
+#[test]
+fn fx_implied_long_end_reproduces_fx_points_and_xccy_quotes() -> Result<()> {
+    // FX ポイント（ON、TN、3M、6M）と XCCY（1Y、18M）が、外貨担保曲線の basis ブロックの
+    // 6 つのピラーを一つのグループで決める。TONA と SOFR の OIS と USDJPY のスポットが上流。
+    let fixture = XccyMarket::new()?;
+    let inputs = fixture.inputs("")?;
+    let ctx = inputs.context();
+    let market = calibrate(&fixture.definition()?, &fixture.quotes()?, &ctx)?;
+    solved(&market)?;
+    let block = rhodium_market::structure::BlockId::new(XccyMarket::implied()?, "basis")?;
+    let group = market
+        .groups()
+        .iter()
+        .find(|group| group.blocks() == [block.clone()])
+        .ok_or_else(|| Error::Internal("missing FX-implied group".into()))?;
+    let rows: Vec<&str> = group.row_quotes().iter().map(QuoteId::as_str).collect();
+    for quote in POINTS
+        .iter()
+        .map(|(quote, _, _)| *quote)
+        .chain(XCCY_QUOTES.iter().map(|(quote, _, _, _)| *quote))
+    {
+        assert!(rows.contains(&quote), "{quote} not in {rows:?}");
+    }
+    assert_eq!(rows.len(), 6);
+    // XCCY を同じ評価で建て直したパーは quote（残差の許容誤差 1e-14 の範囲）。
+    for (_, tenor, quote, _) in XCCY_QUOTES {
+        assert_close!(xccy_par(tenor, &market, &ctx)?, quote, 0.0, 1e-12);
+    }
+    Ok(())
+}
+
+#[test]
+fn fx_implied_calibration_instruments_are_at_par_as_trades() -> Result<()> {
+    // 校正商品を取引として評価すると PV はゼロ（校正の残差を PV に伝えた幅の内）。
+    // FX フォワードは 3M のアウトライト X(S) + pts·pip、FX スワップはスポットの近い受渡と 6M の遠い受渡、
+    // `MtM` XCCY は quote のベーシスで約定する。元本は 100 万 USD（XCCY は JPY 1 円あたり）。
+    let fixture = XccyMarket::new()?;
+    let inputs = fixture.inputs("")?;
+    let ctx = inputs.context();
+    let market = calibrate(&fixture.definition()?, &fixture.quotes()?, &ctx)?;
+    let pair = XccyMarket::pair()?;
+    let notional = Money::new(Currency::USD, 1_000_000.0);
+    let pip = 0.01;
+    let [_, _, three, six] = POINTS;
+    let outright = |points: f64| SPOT + points * pip;
+    let settle = |date: Option<&str>| {
+        XccyMarket::settlement(date.ok_or_else(|| Error::Internal("dated point".into()))?)
+    };
+    let forward = StandardProduct::FxForward(FxForward::new(
+        notional,
+        FxRate::forward(pair, outright(three.2), settle(three.1)?)?,
+        Position::Long,
+    )?);
+    let swap = StandardProduct::FxSwap(FxSwap::new(
+        notional,
+        FxRate::forward(pair, SPOT, XccyMarket::settlement("2025-01-09")?)?,
+        FxRate::forward(pair, outright(six.2), settle(six.1)?)?,
+        Position::Long,
+    )?);
+    // 1 pip（0.01 円）の誤差は 100 万 USD で 1 万円なので、校正の許容誤差 1e-14 の伝わる幅は 1e-6 円の程度。
+    for product in [&forward, &swap] {
+        assert_close!(jpy_pv(product, &market, &ctx)?, 0.0, 0.0, 1e-4);
+    }
+    for (_, tenor, quote, _) in XCCY_QUOTES {
+        let product =
+            XccyMarket::template(tenor)?.build_for_trade(quote, trade_date()?, ctx.refdata)?;
+        // 元本は JPY 1 円。ベーシスの PV の傾き（年金現価、約 1 円/1.0）に許容誤差を掛けた幅。
+        assert_close!(jpy_pv(&product, &market, &ctx)?, 0.0, 0.0, 1e-12);
+    }
+    Ok(())
+}
+
+fn bumped(quotes: &QuoteSet, id: &str, bump: f64) -> Result<QuoteSet> {
+    let id: QuoteId = id.parse()?;
+    QuoteSet::new(
+        quotes
+            .iter()
+            .map(|(key, value)| (key.clone(), if *key == id { value + bump } else { value })),
+    )
+}
+
+#[test]
+fn fx_implied_long_end_update_matches_full_recalibration() -> Result<()> {
+    // XCCY の quote の変更は外貨担保曲線のグループだけを解き直し、TONA、SOFR と FX のマーク
+    // （スポットと短期のポイントの直接マーク）の解を共有する。TONA の quote の変更は TONA と、その下流の
+    // 外貨担保曲線を解き直し、SOFR と FX のマークを共有する。
+    // どちらも全体を解き直した校正と θ が一致する。
+    let fixture = XccyMarket::new()?;
+    let inputs = fixture.inputs("")?;
+    let ctx = inputs.context();
+    let definition = fixture.definition()?;
+    let quotes = fixture.quotes()?;
+    let base = calibrate(&definition, &quotes, &ctx)?;
+    let tona_quote = fixture.tona.quote(&fixture.tona.nodes[2])?;
+    for (id, bump, shared) in [
+        (
+            "USDJPY.XCCY.18M",
+            0.0002,
+            ["JPY-TONA", "USD-SOFR", "FX/USDJPY"].as_slice(),
+        ),
+        (
+            tona_quote.as_str(),
+            0.0001,
+            ["USD-SOFR", "FX/USDJPY"].as_slice(),
+        ),
+    ] {
+        let quotes = bumped(&quotes, id, bump)?;
+        let updated = base.update(&quotes)?;
+        let full = calibrate(&definition, &quotes, &ctx)?;
+        solved(&updated)?;
+        for (lhs, rhs) in updated.market().theta().iter().zip(full.market().theta()) {
+            assert_close!(*lhs, *rhs, 1e-12, 1e-14);
+        }
+        for (old, new) in base.groups().iter().zip(updated.groups()) {
+            let curve = old.blocks()[0].curve().as_str();
+            assert_eq!(
+                ptr::eq(old.outcome(), new.outcome()),
+                shared.contains(&curve),
+                "{id}: {curve}"
+            );
+        }
+    }
+    Ok(())
+}
+
+#[test]
+fn fx_implied_long_end_rejects_pillars_before_fx_points() -> Result<()> {
+    // 長期のピラーは昇順で、最後の FX ポイントの受渡より後でなければならない。
+    let fixture = XccyMarket::new()?;
+    let inputs = fixture.inputs("")?;
+    let ctx = inputs.context();
+    let mut long_end = XccyMarket::long_end()?;
+    long_end.reverse();
+    let uncollateralized = rhodium_pricing::resolve::CollateralKey::Uncollateralized;
+    let base = || {
+        ois::definition(&[
+            (&fixture.tona, &uncollateralized),
+            (&fixture.sofr, &uncollateralized),
+        ])
+    };
+    assert!(matches!(
+        base()?.with_fx([fixture.fx(long_end)?]),
+        Err(Error::InvalidData(message)) if message.contains("strictly increasing")
+    ));
+    let early = vec![FxImpliedLongEnd {
+        pillar: midnight("2025-06-09".parse()?)?,
+        ..XccyMarket::long_end()?.remove(0)
+    }];
+    let definition = base()?.with_fx([fixture.fx(early)?])?;
+    assert!(matches!(
+        calibrate(&definition, &fixture.quotes()?, &ctx),
+        Err(Error::InvalidData(message)) if message.contains("follow the last FX point")
+    ));
+    Ok(())
+}
+
+#[test]
+fn fx_implied_long_end_alone_survives_the_cutoff_of_a_t_plus_one_pair() -> Result<()> {
+    // T+1 の USDJPY を spot、ON と長期の XCCY だけで当てる。締めの前は ON の区間と XCCY のピラー、
+    // 当日受渡の締めの後は ON が連鎖から外れても XCCY の方程式とピラーが残るので、どちらも校正でき、
+    // 外貨担保曲線で XCCY の quote を再現する。TN・SN・日付付きのポイントが無いことは拒否の理由にならない。
+    let fixture = XccyMarket::new()?;
+    let inputs = fixture.inputs_with_spot_lag("", FxSpotLag::OneBusinessDay)?;
+    let mut fx = fixture.fx(XccyMarket::long_end()?)?;
+    let direct = fx
+        .direct
+        .as_mut()
+        .ok_or_else(|| Error::Internal("test needs direct FX".into()))?;
+    direct.quotes.tn = None;
+    direct.quotes.forwards.clear();
+    let uncollateralized = rhodium_pricing::resolve::CollateralKey::Uncollateralized;
+    let definition = ois::definition(&[
+        (&fixture.tona, &uncollateralized),
+        (&fixture.sofr, &uncollateralized),
+    ])?
+    .with_fx([fx])?;
+    // 受渡の時刻は UTC の 16 時。評価日の 0 時は締めの前、16 時 1 秒は締めの後。
+    for as_of in [inputs.context().as_of, "2025-01-07T16:00:01Z".parse()?] {
+        let ctx = Context {
+            as_of,
+            ..inputs.context()
+        };
+        let market = calibrate(&definition, &fixture.quotes()?, &ctx)?;
+        solved(&market)?;
+        for (_, tenor, quote, _) in XCCY_QUOTES {
+            assert_close!(xccy_par(tenor, &market, &ctx)?, quote, 1e-12_f64, 1e-12_f64);
+        }
+    }
+    Ok(())
+}
+
+#[test]
+fn fx_implied_long_end_rejects_an_xccy_convention_of_another_pair() -> Result<()> {
+    // 長期の行は一般の XCCY の binding になるとペアを失うので、参照データに照らす検査で慣行の
+    // 2 通貨をペアと比べる。EUR と JPY を交換する慣行を USDJPY の外貨担保曲線に書くと、
+    // 別の通貨の XCCY で USD の曲線を当ててしまうので拒否する。USDJPY の慣行なら通る。
+    let fixture = XccyMarket::new()?;
+    let mut euro = fixture.sofr.clone();
+    euro.currency = "EUR".parse()?;
+    euro.index = "EUR-ESTR".into();
+    let other = "EURJPY-XCCY";
+    let inputs = Inputs::with_forecasts(
+        &[&fixture.tona, &fixture.sofr, &euro],
+        fixture.tona.as_of,
+        "",
+        |parts, _| {
+            // XCCY の慣行は元本をリセットする脚を一つ持つので、EURJPY の慣行と fixing 指標も置く。
+            let pair: CcyPair = "EURJPY".parse()?;
+            let clock = LocalEventTime {
+                time: "16:00:00".parse()?,
+                zone: "Etc/UTC".parse()?,
+            };
+            parts.fx_conventions.insert(
+                pair,
+                FxConvention {
+                    pair,
+                    spot_lag: FxSpotLag::TwoBusinessDays,
+                    settlement: clock.clone(),
+                    pip_exponent: 2,
+                },
+            );
+            let fixing: IndexId = "EURJPY-FIX".parse()?;
+            let calendar = CalendarSpec::Union(vec![
+                CalendarSpec::Base("JPY".parse()?),
+                CalendarSpec::Base("EUR".parse()?),
+            ]);
+            parts.indices.insert(
+                fixing.clone(),
+                IndexDefinition {
+                    id: fixing.clone(),
+                    kind: IndexKind::FxFixing { pair },
+                    calendar: calendar.clone(),
+                    publication: PublicationRule::Daily(DailyPublication {
+                        business_days: 0,
+                        calendar,
+                        time: clock.time,
+                        zone: clock.zone.clone(),
+                    }),
+                    observation_time: Some(clock),
+                    status: IndexStatus::Active,
+                },
+            );
+            let leg = |case: &ois::Case, fx_reset| -> Result<SwapLegConvention> {
+                let Some(Convention::Ois(ois)) = parts.conventions.get(&case.index.parse()?) else {
+                    return Err(Error::Internal("missing OIS test convention".into()));
+                };
+                let mut periods = ois.fixed.clone();
+                periods.frequency = "3M".parse()?;
+                Ok(SwapLegConvention {
+                    currency: case.currency,
+                    periods,
+                    payment: ois.payment.clone(),
+                    coupon: SwapCouponConvention::Compounded(ois.float.clone()),
+                    fx_reset,
+                })
+            };
+            let legs = [leg(&fixture.tona, None)?, leg(&euro, Some(fixing))?];
+            parts.conventions.insert(
+                other.parse()?,
+                Convention::XccySwap(Box::new(SwapConvention {
+                    spot_lag: "2BD".parse()?,
+                    legs,
+                })),
+            );
+            Ok(())
+        },
+    )?;
+    let misplaced = XccyMarket::long_end()?
+        .into_iter()
+        .map(|row| -> Result<_> {
+            Ok(FxImpliedLongEnd {
+                convention: other.parse()?,
+                ..row
+            })
+        })
+        .collect::<Result<Vec<_>>>()?;
+    assert!(matches!(
+        fixture.fx(misplaced)?.validate_references(inputs.context().refdata),
+        Err(Error::InvalidData(message)) if message.contains("does not exchange USD and JPY")
+    ));
+    fixture
+        .fx(XccyMarket::long_end()?)?
+        .validate_references(fixture.inputs("")?.context().refdata)?;
+    Ok(())
+}
diff --git a/crates/rhodium-products/src/ir.rs b/crates/rhodium-products/src/ir.rs
index a3768283c37f4e00f7d73b4545afbff03d17a21b..1ec9728e78efdb71f3d15d19f97a5ef4a50648b1 100644
--- a/crates/rhodium-products/src/ir.rs
+++ b/crates/rhodium-products/src/ir.rs
@@ -189,7 +189,10 @@ struct ResolvedCompounding {
 
 #[derive(Debug, Clone, Copy)]
 enum FixingKind {
-    TermRate,
+    /// ターム金利。`accrual` は指標の日数係数で測った預金期間の長さ（FRA の決済額の分母を照合する）。
+    TermRate {
+        accrual: f64,
+    },
     Fx(rhodium_core::fx::CcyPair),
     RateFuture,
     Other,
diff --git a/crates/rhodium-products/src/ir/amendment.rs b/crates/rhodium-products/src/ir/amendment.rs
index a4da5f4eee591ddbcbfa39ffe7261b573ff7d5d8..3d769b59b791b64792fe6c705e5e782cd3cd62f0 100644
--- a/crates/rhodium-products/src/ir/amendment.rs
+++ b/crates/rhodium-products/src/ir/amendment.rs
@@ -157,6 +157,10 @@ impl Contract {
     /// `side` の側の受渡だけを残す。`Before` は `from` より前、`After` は `from` 以後の受渡を残す。
     ///
     /// `exercises` を渡すと、`from` より前の判断で `from` 以後の受渡に関わるものを記録の選択肢に畳む。
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
     fn retain(self, from: Instant, side: Side, exercises: Option<&Exercises<'_>>) -> Result<Self> {
         // 各ノードの下に `from` 以後の受渡があるか。arena は子が親より前に並ぶので、前から一度で求まる。
         let reachable = reachable_contracts(&self.contracts, self.root);
diff --git a/crates/rhodium-products/src/ir/analysis.rs b/crates/rhodium-products/src/ir/analysis.rs
index 778fb2e8f7d45c11c4b5086c02326c2cdf5a6b03..0dab8cbd75ae0c2e94c8676d03fcfe39f2eafbf7 100644
--- a/crates/rhodium-products/src/ir/analysis.rs
+++ b/crates/rhodium-products/src/ir/analysis.rs
@@ -20,6 +20,12 @@ use super::{
     reachable_contracts,
 };
 
+/// FRA の決済額の分母の τ と、指標の預金期間の日数係数の一致を判定する相対の幅。
+///
+/// 標準の FRA は同じ計算で τ を作るのでビット単位で一致する。DSL で小数として書いた τ の
+/// 丸め（17 桁の 10 進表記）を受け入れ、日数の数え方の違い（1 日で 1e-3 程度）は区別する。
+const FRA_ACCRUAL_TOLERANCE: f64 = 1e-12;
+
 /// 根から到達する一つの受渡。
 #[derive(Debug)]
 pub struct Cashflow<'a> {
@@ -161,6 +167,10 @@ pub enum Observation<'a> {
 
 /// 到達しうる式が読む観測を、同じ観測を一度ずつ含む一覧にする。
 #[must_use]
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#199 で各 variant を網羅する"
+)]
 pub fn observations(product: &Product) -> Vec<Observation<'_>> {
     let contract = &product.contract;
     let reached = reachable_contracts(&contract.contracts, contract.root);
@@ -225,6 +235,10 @@ impl<'a> ExprRef<'a> {
 
     /// 式の演算と、商品内の子への参照を返す。
     #[must_use]
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
     pub fn expression(self) -> Expression<'a> {
         match self.node() {
             ExprNode::Const(value) => Expression::Constant(value.value()),
@@ -241,12 +255,13 @@ impl<'a> ExprRef<'a> {
             } => {
                 let reference = &self.product.references.fixings[index][point];
                 match reference.kind {
-                    FixingKind::TermRate => Expression::TermFixing(TermFixing {
+                    FixingKind::TermRate { accrual } => Expression::TermFixing(TermFixing {
                         index,
                         point: *point,
                         period: reference.period,
                         publication: reference.publication,
                         forecast_from: self.product.forecast_from(index, *point),
+                        accrual,
                     }),
                     FixingKind::Fx(pair) => Expression::FxFixing(FxFixing {
                         index,
@@ -292,6 +307,10 @@ impl<'a> ExprRef<'a> {
     }
 
     /// 分母が 1 + τL で、分子が同じ L と固定率の差の定数倍なら、L と τ を返す。
+    ///
+    /// 分母と分子の L は同じ観測（同じノード）で、τ は L の預金期間を指標の日数係数で測った長さと
+    /// 一致しなければならない。そのときだけ 1/(1 + τL) が期首から期末までの割引因子になり、
+    /// 期末のクーポンを期首に割り戻す恒等式が成り立つ。違う τ の独自契約は恒等式に当てない。
     #[must_use]
     pub fn fra_settlement(self) -> Option<(TermFixing<'a>, f64)> {
         let [numerator, denominator] = self.binary_children(Expression::Divide)?;
@@ -318,7 +337,8 @@ impl<'a> ExprRef<'a> {
             _ => return None,
         };
         let (tau, rate) = factors(interest)?;
-        (tau > 0.0 && rate.index() == observed.index()).then_some((fixing, tau))
+        let accrual = (tau - fixing.accrual).abs() <= FRA_ACCRUAL_TOLERANCE * fixing.accrual;
+        (tau > 0.0 && accrual && rate.index() == observed.index()).then_some((fixing, tau))
     }
 
     fn binary_children(self, operator: Expression<'_>) -> Option<[Self; 2]> {
@@ -374,7 +394,10 @@ pub struct FxFixing<'a> {
     pub point: ObservationPoint,
     /// スポットを観測する瞬間。
     pub observed_at: Instant,
-    /// 値洗いする元本をこの fixing で換算する受渡の瞬間。
+    /// 元本リセットで、この fixing で換算した新しい元本を交換する瞬間（脚の通貨の決済時刻）。
+    ///
+    /// fixing の予想を読む受渡ではない。同じ指標・観測日の公表値は一つなので、予想は通常の
+    /// fixing と同じくペアのスポットの受渡で読み、この瞬間とのずれは支払いの遅れとして扱う。
     pub value_at: Option<Instant>,
     /// 観測のペアの向き。
     pub pair: CcyPair,
@@ -404,6 +427,8 @@ pub struct TermFixing<'a> {
     pub publication: Instant,
     /// 期間の途中から予想する場合の開始の瞬間。
     pub forecast_from: Option<Instant>,
+    /// 指標の日数係数で測った預金期間の長さ。
+    pub accrual: f64,
 }
 
 /// 簡約済み IR に残る翌日物の複利。
@@ -477,17 +502,25 @@ pub struct RfrObservation {
     pub publication: Instant,
 }
 
-/// M0 のリニア評価への適用可否。
+/// リニア評価への適用可否。
 #[derive(Debug, Clone, Copy, PartialEq, Eq)]
 pub enum Linearity {
-    /// M0 のリニア評価で扱える式。
+    /// リニア評価で扱える式。
     Linear,
-    /// M0 のリニア評価を適用できない理由。
+    /// リニア評価を適用できない理由。
     NonLinear(&'static str),
 }
 
-/// 簡約後の契約に M0 のリニア方式を適用できるかを判定する。
+/// 簡約後の契約にリニア方式を適用できるかを判定する。
+///
+/// 受渡額が定数、ターム金利と翌日物の複利、FX の fixing（その逆数を含む）、金利先物の清算価格の
+/// 線形結合で、予想の積は FX と金利の一つずつ（MtM XCCY の元本リセット）、予想での割り算は FX の
+/// 逆数と FRA の決済だけのときにリニアとする。
 #[must_use]
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#199 で各 variant を網羅する"
+)]
 pub fn linearity(product: &Product) -> Linearity {
     // 受渡額が含む予想の種類を追い、FX と金利の一つずつの積だけを近似の対象にする。
     const RATE: u8 = 1;
diff --git a/crates/rhodium-products/src/ir/builder.rs b/crates/rhodium-products/src/ir/builder.rs
index b55760a9db624eee7b6c76762c8fc3a98eb1d0a6..c5ea05b7cd3a60f38eb6e33d79903d43353aee5d 100644
--- a/crates/rhodium-products/src/ir/builder.rs
+++ b/crates/rhodium-products/src/ir/builder.rs
@@ -149,6 +149,8 @@ impl<'b> Builder<'b> {
     }
 
     /// 元本換算に使う FX fixing と、新しい元本の受渡の瞬間を指定する。
+    ///
+    /// `value_at` は元本を交換する瞬間で、fixing の予想はペアのスポットの受渡で読む。
     pub fn fx_reset(
         &mut self,
         index: &IndexId,
diff --git a/crates/rhodium-products/src/ir/nodes.rs b/crates/rhodium-products/src/ir/nodes.rs
index b7f4eb948753d51d4b01bcdeab65d86a460d6e1e..ec8dd242a4900bf3d1540e86cc2b47a2545f339d 100644
--- a/crates/rhodium-products/src/ir/nodes.rs
+++ b/crates/rhodium-products/src/ir/nodes.rs
@@ -55,6 +55,10 @@ impl ContractNode {
         }
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
     pub(super) fn visit_children(&self, mut visit: impl FnMut(Ix)) {
         match self {
             Self::Cond { yes, no, .. } => {
@@ -270,6 +274,10 @@ impl ExprNode {
         }
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
     pub(super) fn children(&self) -> &[Ix] {
         match self {
             Self::Unary(_, xs) | Self::WithWarnings(_, xs) => xs,
@@ -280,6 +288,10 @@ impl ExprNode {
         }
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
     pub(super) fn try_map_children(&mut self, mut map: impl FnMut(Ix) -> Result<Ix>) -> Result<()> {
         let children: &mut [Ix] = match self {
             Self::Unary(_, xs) | Self::WithWarnings(_, xs) => xs,
diff --git a/crates/rhodium-products/src/ir/reduce.rs b/crates/rhodium-products/src/ir/reduce.rs
index f24e358c9619538671e0b3d8a1fe50f399a9bdd3..0fec72b74554de60442de07257e98d5a8f4537a5 100644
--- a/crates/rhodium-products/src/ir/reduce.rs
+++ b/crates/rhodium-products/src/ir/reduce.rs
@@ -181,6 +181,7 @@ impl Copies {
     }
 
     fn mapped(&self, reach: Reach, ix: Ix) -> Result<Ix> {
+        // 不変条件: 作業の stack は子を親より先に簡約するので、親を組むときに子の結果がある。
         self.get(reach, ix)
             .ok_or_else(|| Error::Internal(format!("unvisited IR node {ix}")))
     }
@@ -270,6 +271,10 @@ impl<'a> Reducer<'a> {
         self.contract_ids.mapped(Reach::Required, root)
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
     fn visit_contract(
         &mut self,
         node: &ContractNode,
@@ -303,6 +308,8 @@ impl<'a> Reducer<'a> {
                 // 判定が評価時点より後でも Required のまま、欠損をエラーにする。
                 let cond = self.expr(*cond, reach)?;
                 if *at <= self.as_of {
+                    // 不変条件: 判定時刻を過ぎた Required の条件の fixing は、欠ければ self.expr が
+                    // エラーを返すので、ここでは真偽の定数に畳まれている。
                     let ExprNode::Bool(value) = self.exprs.nodes[cond as usize] else {
                         return Err(Error::Internal(format!(
                             "condition at {at} remains unresolved"
@@ -470,6 +477,10 @@ impl<'a> Reducer<'a> {
         self.expr_ids.mapped(root_reach, root)
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
     fn visit_expr(
         &mut self,
         node: &'a ExprNode,
@@ -577,6 +588,10 @@ impl<'a> Reducer<'a> {
         Ok(())
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
     fn logic(
         &mut self,
         node: &ExprNode,
@@ -630,6 +645,10 @@ impl<'a> Reducer<'a> {
         Ok(())
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
     fn fold_expr(&mut self, node: ExprNode, reach: Reach) -> Result<Ix> {
         let constant = |ix| match self.exprs.nodes[ix as usize] {
             ExprNode::Const(value) => Some(value.value()),
@@ -843,6 +862,10 @@ impl<'a> Reducer<'a> {
         })))
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
     fn credit(
         &mut self,
         node: &ExprNode,
@@ -887,6 +910,10 @@ fn known_constant(reach: Reach, value: f64) -> Result<ExprNode> {
     })
 }
 
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#199 で各 variant を網羅する"
+)]
 fn select_children(node: &ExprNode) -> Result<[Ix; 3]> {
     match node {
         ExprNode::Select(children) => Ok(*children),
diff --git a/crates/rhodium-products/src/ir/tests.rs b/crates/rhodium-products/src/ir/tests.rs
index 9f38fc367470415e679d33f8ab0eb991d3d69ae8..c9abaee0eddb93e784e1d702fe3e4ff154c0b20c 100644
--- a/crates/rhodium-products/src/ir/tests.rs
+++ b/crates/rhodium-products/src/ir/tests.rs
@@ -45,6 +45,21 @@ fn refdata() -> Result<RefData> {
     RefData::new(refdata_parts()?)
 }
 
+/// USDJPY の受渡慣行（T+2、pip 0.01 円）を、通貨と同じ決済時刻で足す。
+fn add_usd_jpy(parts: &mut RefDataParts, clock: &LocalEventTime) -> Result<()> {
+    let usd_jpy: rhodium_core::fx::CcyPair = "USDJPY".parse()?;
+    parts.fx_conventions.insert(
+        usd_jpy,
+        rhodium_core::refdata::FxConvention {
+            pair: usd_jpy,
+            spot_lag: rhodium_core::refdata::FxSpotLag::TwoBusinessDays,
+            settlement: clock.clone(),
+            pip_exponent: 2,
+        },
+    );
+    Ok(())
+}
+
 fn refdata_parts() -> Result<RefDataParts> {
     let mut parts = RefDataParts::default();
     let calendar = "JPTO".parse()?;
@@ -91,6 +106,8 @@ fn refdata_parts() -> Result<RefDataParts> {
             zone: clock.zone.clone(),
         }))
     };
+    // USDJPY の FX fixing の受渡日は、ペアの慣行のスポットラグで決まる。
+    add_usd_jpy(&mut parts, &clock)?;
     for (name, kind, publication, observation_time) in [
         (
             "JPY-TONA",
@@ -1026,6 +1043,10 @@ fn undecided_branches_defer_numerical_failures_of_published_aggregates() -> Resu
 }
 
 #[test]
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#196 で各 variant を網羅する"
+)]
 fn pending_condition_requires_its_published_fixings() -> Result<()> {
     // 到達した Cond の条件は判定時刻に必ず評価するので、判定が評価時点より後でも、
     // 公表済みの値が欠ければエラーにする。分かる項は畳み、まだ公表されない項だけを残す。
@@ -1632,7 +1653,7 @@ fn analysis_preserves_payment_occurrences_and_branch_paths() -> Result<()> {
     Ok(())
 }
 
-/// 予想する複利 g と定数 2 から支払額を組み、M0 のリニア方式への適用可否を判定する。
+/// 予想する複利 g と定数 2 から支払額を組み、リニア方式への適用可否を判定する。
 fn linearity_of(
     refdata: &RefData,
     observations: &ObservationSchedule,
diff --git a/crates/rhodium-products/src/ir/tests/lifecycle.rs b/crates/rhodium-products/src/ir/tests/lifecycle.rs
index 23f524325941290edf4c1c19b8329ffa44aef633..45569717f9c44234bd2b967d83ca4a22037b72a3 100644
--- a/crates/rhodium-products/src/ir/tests/lifecycle.rs
+++ b/crates/rhodium-products/src/ir/tests/lifecycle.rs
@@ -28,6 +28,10 @@ fn trade(contract: Contract) -> Result<Trade> {
     })
 }
 
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#196 で各 variant を網羅する"
+)]
 fn amounts(product: &Product) -> Result<Vec<f64>> {
     cashflows(product)
         .iter()
@@ -398,3 +402,61 @@ fn fallback_terms_use_explicit_precedence_and_warn_on_the_legacy_rate() -> Resul
         .insert(term, FallbackSelection::IndexDefault);
     requires_replacement(&trade)
 }
+
+#[test]
+fn trade_fallback_terms_reject_invalid_references() -> Result<()> {
+    // 取引のフォールバックの指定は IR を作る前に参照データと照らす。ターム金利でない指標への個別の指定、
+    // 未知の指標との個別の合意、参照データと合わない固有の代替は、どれも InvalidData。
+    let (data, term, overnight) = fallback_data(IndexStatus::Active)?;
+    let mut trade = fixing_trade(&term, "2027-01-08".parse()?)?;
+    let custom = |index: &str| -> Result<FallbackSelection> {
+        let replacement = data
+            .fallback(&term)
+            .ok_or_else(|| Error::Internal("missing default fallback".into()))?
+            .replacement
+            .clone();
+        Ok(FallbackSelection::Custom {
+            trigger: FallbackTrigger::Cessation,
+            replacement: FallbackReplacement {
+                index: index.parse()?,
+                ..replacement
+            },
+        })
+    };
+    let reject = |trade: &Trade, expected: &str| {
+        assert!(
+            matches!(trade.to_product(&data), Err(Error::InvalidData(message)) if message.contains(expected)),
+            "{expected}"
+        );
+    };
+    // 個別の指定はターム金利だけに書ける（翌日物の TONA は対象外）。
+    trade
+        .fallback
+        .by_index
+        .insert(overnight.clone(), FallbackSelection::Disabled);
+    reject(&trade, "is not a term rate");
+    trade.fallback.by_index.clear();
+    // 個別に合意した指標は参照データに無ければならない。
+    trade.fallback.bilateral.insert("JPY-UNKNOWN".parse()?);
+    reject(&trade, "JPY-UNKNOWN");
+    trade.fallback.bilateral.clear();
+    // 取引固有の代替の指標も参照データの翌日物でなければならない。
+    trade
+        .fallback
+        .by_index
+        .insert(term.clone(), custom("JPY-UNKNOWN")?);
+    reject(&trade, "JPY-UNKNOWN");
+    trade
+        .fallback
+        .by_index
+        .insert(term.clone(), custom(term.as_str())?);
+    assert!(matches!(
+        trade.to_product(&data),
+        Err(Error::InvalidData(_))
+    ));
+    // 正しい固有の代替は受け付ける。
+    let valid = custom(overnight.as_str())?;
+    trade.fallback.by_index.insert(term.clone(), valid);
+    trade.to_product(&data)?;
+    Ok(())
+}
diff --git a/crates/rhodium-products/src/ir/tests/standard.rs b/crates/rhodium-products/src/ir/tests/standard.rs
index cff5ab07f04c7fa5e0d80f7363eeb6003451636f..8463a1a57dddbb26c10260216af7b1e814ca876d 100644
--- a/crates/rhodium-products/src/ir/tests/standard.rs
+++ b/crates/rhodium-products/src/ir/tests/standard.rs
@@ -213,6 +213,10 @@ fn fixings(refdata: &RefData, index: &str, rows: &[(&str, f64)]) -> Result<Fixin
 
 /// `reduce` 後の受渡を (時刻, 確定した額) の列にする。未確定の額は `None`。
 /// `give` は構築時に額の符号反転になるので、契約は `All` と `Pay` だけからなる。
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#196 で各 variant を網羅する"
+)]
 fn reduced_payments(
     product: &Product,
     as_of: &str,
diff --git a/crates/rhodium-products/src/ir/validate.rs b/crates/rhodium-products/src/ir/validate.rs
index d16aa792a94ed543874ad4292df2d46ae5cfadff..7b80f09bbc0e3d8c20078bcf1aedeb92e311c426 100644
--- a/crates/rhodium-products/src/ir/validate.rs
+++ b/crates/rhodium-products/src/ir/validate.rs
@@ -5,8 +5,10 @@ use std::{
 
 use rhodium_core::{
     Error, Result, Warning,
+    calendar::CalendarSpec,
     currency::Asset,
     date::DateInterval,
+    day_count::{AccrualContext, DayCount, year_fraction},
     fixing::ObservationPoint,
     id::IndexId,
     refdata::{IndexDefinition, IndexKind, IndexStatus, RefData},
@@ -36,6 +38,7 @@ impl ObservationCache {
             let mut entries = self
                 .0
                 .lock()
+                // 不変条件: lock を持つ間は表の参照を出し入れするだけで panic しないので、poison しない。
                 .map_err(|_| Error::Internal("RFR observation cache lock poisoned".into()))?;
             match entries.entry((definition.id.clone(), Arc::clone(rows))) {
                 Entry::Occupied(entry) => {
@@ -105,6 +108,10 @@ pub(super) fn validate(
     Ok((references, warnings))
 }
 
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#199 で各 variant を網羅する"
+)]
 fn known_at(
     node: &mut ExprNode,
     known: &[Option<Instant>],
@@ -241,6 +248,10 @@ fn resolve_compounding(
     }))
 }
 
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#197 で各 variant を網羅する"
+)]
 fn fixing(
     refdata: &RefData,
     definition: &IndexDefinition,
@@ -250,6 +261,33 @@ fn fixing(
     let (publication, mut warnings) = refdata.publication(id, point)?;
     let (period, notices) = refdata.observation_period(id, point)?;
     warnings.extend(notices);
+    let kind = match (&definition.kind, point) {
+        (IndexKind::TermRate { day_count, .. }, ObservationPoint::Date(date)) => {
+            // 標準の FRA と同じく、預金期間を指標の日数係数で測る（`Fra::build_contract`）。
+            let (deposit, notices) = refdata.term_rate_dates(id, date)?;
+            warnings.extend(notices);
+            let calendar = match day_count {
+                DayCount::Bus252 { calendar } => {
+                    Some(refdata.calendar(&CalendarSpec::Base(calendar.clone()))?)
+                }
+                _ => None,
+            };
+            let context = AccrualContext {
+                calendar: calendar.as_ref(),
+                maturity: Some(deposit.end()),
+                ..AccrualContext::default()
+            };
+            let (accrual, notices) =
+                year_fraction(day_count, deposit.start(), deposit.end(), &context)?;
+            warnings.extend(notices);
+            super::FixingKind::TermRate {
+                accrual: accrual.value(),
+            }
+        }
+        (IndexKind::FxFixing { pair }, _) => super::FixingKind::Fx(*pair),
+        (IndexKind::RateFuturesSettlement { .. }, _) => super::FixingKind::RateFuture,
+        _ => super::FixingKind::Other,
+    };
     Ok(ResolvedFixing {
         publication,
         period,
@@ -258,12 +296,7 @@ fn fixing(
             (IndexStatus::Ceased { since }, ObservationPoint::Date(date)) if date >= since),
         non_representative: matches!((definition.status, point),
             (IndexStatus::NonRepresentative { since }, ObservationPoint::Date(date)) if date >= since),
-        kind: match (&definition.kind, point) {
-            (IndexKind::TermRate { .. }, ObservationPoint::Date(_)) => super::FixingKind::TermRate,
-            (IndexKind::FxFixing { pair }, _) => super::FixingKind::Fx(*pair),
-            (IndexKind::RateFuturesSettlement { .. }, _) => super::FixingKind::RateFuture,
-            _ => super::FixingKind::Other,
-        },
+        kind,
     })
 }
 
diff --git a/crates/rhodium-products/src/standard/leg.rs b/crates/rhodium-products/src/standard/leg.rs
index 75ff080565058a49ed68ea7358027373569e86ed..dfb1829277937da050a3e411a89deff0dadf0537 100644
--- a/crates/rhodium-products/src/standard/leg.rs
+++ b/crates/rhodium-products/src/standard/leg.rs
@@ -4,24 +4,25 @@ use rhodium_core::{
     Error, Result, Warning,
     calendar::{BusinessDayConvention, Calendar, CalendarSpec},
     currency::{Asset, Currency},
-    date::Date,
-    day_count::{DayCount, year_fraction},
+    date::{Date, DateInterval},
+    day_count::{AccrualContext, DayCount, year_fraction},
     fixing::ObservationPoint,
     id::IndexId,
     money::Money,
     observation::{
         Aggregation, ObservationRule, SpreadTreatment, observation_schedule_with_calendars,
     },
-    quantity::{Rate, Simple, Spread, accrue},
+    quantity::{Rate, Simple, Spread, YearFraction, accrue},
     refdata::{CalendarCache, IndexKind, RefData},
     schedule::{
-        PeriodAnchor, PeriodCombination, ResetRule, Schedule, ScheduleSpec, SettlementContext,
-        generate_with_calendars, reset_dates,
+        Period, PeriodAnchor, PeriodCombination, ResetRule, Schedule, ScheduleSpec,
+        SettlementContext, generate_with_calendars, reset_dates,
     },
     tenor::{DateOffset, TenorUnit},
+    time::Instant,
 };
 
-use crate::ir::{Builder, ContractIx, ExprIx};
+use crate::ir::{Builder, Contract, ContractIx, ExprIx, Product};
 
 use super::{PaymentRounding, invalid, settlement, validate_notional};
 
@@ -224,24 +225,7 @@ impl Leg {
     ) -> Result<ContractIx<'b>> {
         let currency = refdata.currency(self.notional.currency(refdata)?)?;
         let asset = Asset::Cash(currency.code);
-        let settlement_calendar = CalendarSpec::Base(currency.settlement_calendar.clone());
-        let (schedule, notices) = generate_with_calendars(
-            &self.schedule,
-            trade_date,
-            &SettlementContext {
-                zones: refdata.zones(),
-                default_time: Some(&currency.settlement_time),
-                calendar: Some(&settlement_calendar),
-            },
-            calendars,
-        )?;
-        warnings.extend(notices);
-        if schedule.combine() != PeriodCombination::Sum {
-            return Err(Error::Unsupported(
-                "leg payments with compounding across accrual periods".into(),
-            ));
-        }
-        schedule.validate_day_count(&self.day_count)?;
+        let schedule = self.resolve_schedule(trade_date, currency, refdata, calendars, warnings)?;
         let coupon = self.resolve_coupon(refdata, &schedule, currency.code, warnings)?;
         let calendar = day_count_calendar(&self.day_count, refdata)?;
         let notionals = self.notionals(builder, &schedule, currency, refdata, warnings)?;
@@ -304,6 +288,104 @@ impl Leg {
         })
     }
 
+    fn resolve_schedule(
+        &self,
+        trade_date: Date,
+        currency: &rhodium_core::refdata::CurrencyDefinition,
+        refdata: &RefData,
+        calendars: &mut CalendarCache<'_>,
+        warnings: &mut Vec<Warning>,
+    ) -> Result<Schedule> {
+        let settlement_calendar = CalendarSpec::Base(currency.settlement_calendar.clone());
+        let (schedule, notices) = generate_with_calendars(
+            &self.schedule,
+            trade_date,
+            &SettlementContext {
+                zones: refdata.zones(),
+                default_time: Some(&currency.settlement_time),
+                calendar: Some(&settlement_calendar),
+            },
+            calendars,
+        )?;
+        warnings.extend(notices);
+        if schedule.combine() != PeriodCombination::Sum {
+            return Err(Error::Unsupported(
+                "leg payments with compounding across accrual periods".into(),
+            ));
+        }
+        schedule.validate_day_count(&self.day_count)?;
+        Ok(schedule)
+    }
+
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
+    pub(super) fn caplet_coupons(
+        &self,
+        index: &IndexId,
+        trade_date: Date,
+        refdata: &RefData,
+    ) -> Result<(Vec<CapletCoupon>, Vec<Warning>)> {
+        let currency = refdata.currency(self.notional.currency(refdata)?)?;
+        let mut calendars = refdata.calendar_cache();
+        let mut warnings = Vec::new();
+        let schedule =
+            self.resolve_schedule(trade_date, currency, refdata, &mut calendars, &mut warnings)?;
+        let coupon = self.resolve_coupon(refdata, &schedule, currency.code, &mut warnings)?;
+        let first = match &coupon {
+            ResolvedCoupon::Term {
+                index: observed, ..
+            } if *observed == index => 1,
+            ResolvedCoupon::Compounded {
+                index: observed, ..
+            } if *observed == index => 0,
+            _ => return Err(invalid("cap convention does not match its floating index")),
+        };
+        let calendar = day_count_calendar(&self.day_count, refdata)?;
+        let coupons = schedule
+            .payments()
+            .iter()
+            .skip(first)
+            .map(|payment| {
+                let periods = payment.periods();
+                if periods.len() != 1 {
+                    return Err(invalid("cap calibration requires one period per payment"));
+                }
+                let i = periods.start;
+                let period = &schedule.periods()[i];
+                let context = schedule
+                    .accrual_context(i, calendar.as_ref())
+                    .ok_or_else(|| Error::Internal("missing schedule accrual context".into()))?;
+                let fraction = self.coupon_fraction(period, &context, &mut warnings)?;
+                let mut observation = None;
+                let contract = Contract::build(|builder| {
+                    let (amount, dates) = Self::coupon(
+                        builder,
+                        &coupon,
+                        i,
+                        period,
+                        &context,
+                        &fraction,
+                        &mut calendars,
+                        &mut warnings,
+                    )?;
+                    observation = dates;
+                    Ok(builder.pay(Asset::Cash(currency.code), amount, payment.at()))
+                })?;
+                let observation = observation
+                    .ok_or_else(|| Error::Internal("caplet requires a floating coupon".into()))?;
+                Ok(CapletCoupon {
+                    product: Product::new(contract, refdata)?,
+                    pay: payment.at(),
+                    fraction,
+                    observation,
+                })
+            })
+            .collect::<Result<Vec<_>>>()?;
+        Ok((coupons, warnings))
+    }
+
     fn notionals<'b>(
         &self,
         builder: &mut Builder<'b>,
@@ -322,15 +404,15 @@ impl Leg {
                 let IndexKind::FxFixing { pair } = definition.kind else {
                     return Err(invalid(format!("{fixing} is not an FX fixing")));
                 };
-                // スポット受渡（2 営業日後）が期首になる日に FX を観測し、換算した元本を期首に受け渡す。
-                // スポットラグは参照データに持たないので、T+1 のペアなど別の慣行の fixing 日は
-                // 明示日程の `reset` で指定する（この日はターム金利の fixing と共通になる）。
+                // 期首がスポット受渡になる日（ペアの慣行のスポットラグを両通貨の休日の和で遡った日）に
+                // FX を観測し、換算した元本を期首に受け渡す。明示日程の `reset` はこの規則より優先する。
+                let spot_lag = refdata.fx_convention(pair)?.spot_lag;
                 let (dates, notices) = reset_dates(
                     schedule,
                     &ResetRule {
                         anchor: PeriodAnchor::PeriodStart,
-                        offset: DateOffset::new(-2, TenorUnit::BusinessDay),
-                        calendar: definition.calendar.clone(),
+                        offset: DateOffset::new(-spot_lag.business_days(), TenorUnit::BusinessDay),
+                        calendar: refdata.fx_calendar_spec(pair)?,
                         bdc: BusinessDayConvention::Preceding,
                     },
                     refdata.calendars(),
@@ -452,73 +534,135 @@ impl Leg {
                 let context = schedule
                     .accrual_context(i, calendar)
                     .ok_or_else(|| Error::Internal("missing schedule accrual context".into()))?;
-                let interval = period.accrual();
-                let (fraction, notices) =
-                    year_fraction(&self.day_count, interval.start(), interval.end(), &context)?;
-                warnings.extend(notices);
-                let interest = match coupon {
-                    ResolvedCoupon::Fixed(rate) => builder.constant(accrue(rate, &fraction)?),
-                    ResolvedCoupon::Term {
-                        index,
-                        dates,
-                        gearing,
-                        spread,
-                    } => {
-                        let mut rate = builder.fixing(index, ObservationPoint::Date(dates[i]));
-                        if gearing.to_bits() != 1.0_f64.to_bits() {
-                            let coefficient = builder.constant(*gearing);
-                            rate = builder.mul([coefficient, rate]);
-                        }
-                        if spread.value() != 0.0 {
-                            let spread = builder.constant(spread.value());
-                            rate = builder.add([rate, spread]);
-                        }
-                        let coupon_fraction = builder.constant(fraction.value());
-                        builder.mul([rate, coupon_fraction])
-                    }
-                    ResolvedCoupon::Compounded {
-                        index,
-                        spread,
-                        rule,
-                        basis,
-                    } => {
-                        let (observations, notices) = observation_schedule_with_calendars(
-                            period, rule, basis, &context, calendars,
-                        )?;
-                        warnings.extend(notices);
-                        let annualization = observations.annualization().value();
-                        if annualization <= 0.0 {
-                            return Err(invalid(
-                                "RFR coupon requires a positive annualization fraction",
-                            ));
-                        }
-                        let growth = builder.compounded(index, &observations, rule, *spread)?;
-                        let one = builder.constant(1.0);
-                        let excess = builder.sub(growth, one);
-                        let denominator = builder.constant(annualization);
-                        let mut rate = builder.div(excess, denominator);
-                        if matches!(
-                            rule.aggregation(),
-                            Aggregation::Compound {
-                                spread: SpreadTreatment::Exclusive,
-                                ..
-                            }
-                        ) && spread.value() != 0.0
-                        {
-                            let margin = builder.constant(spread.value());
-                            rate = builder.add([rate, margin]);
-                        }
-                        if let Some(rounding) = rule.rounding().aggregated {
-                            rate = builder.round(rate, rounding);
-                        }
-                        let coupon_fraction = builder.constant(fraction.value());
-                        builder.mul([rate, coupon_fraction])
-                    }
-                };
+                let fraction = self.coupon_fraction(period, &context, warnings)?;
+                let (interest, _) = Self::coupon(
+                    builder, coupon, i, period, &context, &fraction, calendars, warnings,
+                )?;
                 Ok(builder.mul([notionals[i], interest]))
             })
             .collect()
     }
+
+    fn coupon_fraction(
+        &self,
+        period: &Period,
+        context: &AccrualContext<'_>,
+        warnings: &mut Vec<Warning>,
+    ) -> Result<YearFraction> {
+        let interval = period.accrual();
+        let (fraction, notices) =
+            year_fraction(&self.day_count, interval.start(), interval.end(), context)?;
+        warnings.extend(notices);
+        Ok(fraction)
+    }
+
+    #[expect(
+        clippy::too_many_arguments,
+        reason = "標準脚とキャップレットで同じクーポン式を使う"
+    )]
+    fn coupon<'b>(
+        builder: &mut Builder<'b>,
+        coupon: &ResolvedCoupon<'_>,
+        i: usize,
+        period: &Period,
+        context: &AccrualContext<'_>,
+        fraction: &YearFraction,
+        calendars: &mut CalendarCache<'_>,
+        warnings: &mut Vec<Warning>,
+    ) -> Result<(ExprIx<'b>, Option<CapletObservation>)> {
+        Ok(match coupon {
+            ResolvedCoupon::Fixed(rate) => (builder.constant(accrue(rate, fraction)?), None),
+            ResolvedCoupon::Term {
+                index,
+                dates,
+                gearing,
+                spread,
+            } => {
+                let mut rate = builder.fixing(index, ObservationPoint::Date(dates[i]));
+                if gearing.to_bits() != 1.0_f64.to_bits() {
+                    let coefficient = builder.constant(*gearing);
+                    rate = builder.mul([coefficient, rate]);
+                }
+                if spread.value() != 0.0 {
+                    let spread = builder.constant(spread.value());
+                    rate = builder.add([rate, spread]);
+                }
+                let coupon_fraction = builder.constant(fraction.value());
+                (
+                    builder.mul([rate, coupon_fraction]),
+                    Some(CapletObservation::Term(dates[i])),
+                )
+            }
+            ResolvedCoupon::Compounded {
+                index,
+                spread,
+                rule,
+                basis,
+            } => {
+                let (observations, notices) =
+                    observation_schedule_with_calendars(period, rule, basis, context, calendars)?;
+                warnings.extend(notices);
+                let annualization = observations.annualization().value();
+                if annualization <= 0.0 {
+                    return Err(invalid(
+                        "RFR coupon requires a positive annualization fraction",
+                    ));
+                }
+                let growth = builder.compounded(index, &observations, rule, *spread)?;
+                let one = builder.constant(1.0);
+                let excess = builder.sub(growth, one);
+                let denominator = builder.constant(annualization);
+                let mut rate = builder.div(excess, denominator);
+                if matches!(
+                    rule.aggregation(),
+                    Aggregation::Compound {
+                        spread: SpreadTreatment::Exclusive,
+                        ..
+                    }
+                ) && spread.value() != 0.0
+                {
+                    let margin = builder.constant(spread.value());
+                    rate = builder.add([rate, margin]);
+                }
+                if let Some(rounding) = rule.rounding().aggregated {
+                    rate = builder.round(rate, rounding);
+                }
+                let coverage = observations.coverage();
+                let dates = coverage
+                    .first()
+                    .zip(coverage.last())
+                    .ok_or_else(|| Error::InvalidData("empty RFR observation period".into()))?;
+                let observation = CapletObservation::Overnight(DateInterval::new(
+                    dates.0.start(),
+                    dates.1.end(),
+                )?);
+                let coupon_fraction = builder.constant(fraction.value());
+                (builder.mul([rate, coupon_fraction]), Some(observation))
+            }
+        })
+    }
+}
+
+/// キャップレットに対応する、元本 1 の変動クーポンとその時間軸。
+#[derive(Debug)]
+pub struct CapletCoupon {
+    /// 既存のリニア評価で変動クーポンを評価するための IR。
+    pub product: Product,
+    /// 割引に使う支払の瞬間。
+    pub pay: Instant,
+    /// クーポンの計算期間の日数係数。
+    pub fraction: YearFraction,
+    /// 固定日または翌日物の観測期間。
+    pub observation: CapletObservation,
+}
+
+/// キャップレットの不確実性を測る時間の境界。
+#[derive(Debug, Clone, Copy, PartialEq, Eq)]
+pub enum CapletObservation {
+    /// ターム金利の fixing 日。
+    Term(Date),
+    /// 翌日物を集約する観測期間。
+    Overnight(DateInterval),
 }
 
 enum ResolvedCoupon<'a> {
@@ -549,6 +693,10 @@ fn validate_explicit_basis(schedule: &ScheduleSpec, basis: &DayCount) -> Result<
     Ok(())
 }
 
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#197 で各 variant を網羅する"
+)]
 pub(super) fn day_count_calendar(basis: &DayCount, refdata: &RefData) -> Result<Option<Calendar>> {
     match basis {
         DayCount::Bus252 { calendar } => refdata
diff --git a/crates/rhodium-products/src/standard/rates.rs b/crates/rhodium-products/src/standard/rates.rs
index 1ecd743ce189663d5f2d2980cf46ee170aea7e51..953822fe503334288b5eb96a83736a81f88b9cff 100644
--- a/crates/rhodium-products/src/standard/rates.rs
+++ b/crates/rhodium-products/src/standard/rates.rs
@@ -115,6 +115,10 @@ impl Deposit {
             .map_err(|error| error.context("deposit"))
     }
 
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#197 で各 variant を網羅する"
+    )]
     pub(super) fn build_contract(
         &self,
         trade_date: Date,
@@ -338,6 +342,49 @@ impl Swap {
         Ok(Self(SwapTerms::Legs(Box::new(legs))))
     }
 
+    /// 標準スワップの変動脚から、元本 1 のキャップ校正用クーポンを作る。
+    ///
+    /// # Errors
+    ///
+    /// 変動脚を一意に選べないか、指標や日程がキャップの条件に合わない場合。
+    pub fn caplet_coupons(
+        &self,
+        index: &rhodium_core::id::IndexId,
+        trade_date: Date,
+        refdata: &RefData,
+    ) -> Result<(Vec<super::leg::CapletCoupon>, Vec<Warning>)> {
+        let resolved;
+        let legs = match &self.0 {
+            SwapTerms::Convention(ois) => {
+                resolved = ois.legs(refdata)?;
+                &resolved
+            }
+            SwapTerms::Legs(legs) => legs.as_ref(),
+        };
+        let floating = match legs {
+            [
+                Leg {
+                    coupon: Coupon::Fixed { .. },
+                    ..
+                },
+                second,
+            ] => second,
+            [
+                first,
+                Leg {
+                    coupon: Coupon::Fixed { .. },
+                    ..
+                },
+            ] => first,
+            _ => {
+                return Err(invalid(
+                    "cap calibration requires one fixed and one floating leg",
+                ));
+            }
+        };
+        floating.caplet_coupons(index, trade_date, refdata)
+    }
+
     /// 各脚の日程を解決して検証済み IR にする。
     ///
     /// # Errors
diff --git a/crates/rhodium-products/src/standard/template.rs b/crates/rhodium-products/src/standard/template.rs
index e9e7a629c29a222ad197616ac88e58cb794d4bfd..87b4347f851e4d072e158026eb4056a470172bff 100644
--- a/crates/rhodium-products/src/standard/template.rs
+++ b/crates/rhodium-products/src/standard/template.rs
@@ -266,11 +266,13 @@ impl InstrumentTemplate {
                 )?))
             }
             (Self::FxForward { tenor, .. }, Convention::FxForward(fx)) => {
-                let calendar = refdata.calendar(&fx.calendar)?;
-                let (period, _) = fx_dates(trade_date()?, *tenor, fx, &calendar)?;
+                // スポットラグ、両通貨の休日の和と受渡の時刻はペアの慣行から読む。
+                let spot = refdata.fx_convention(fx.pair)?;
+                let calendar = refdata.calendar(&refdata.fx_calendar_spec(fx.pair)?)?;
+                let (period, _) = fx_dates(trade_date()?, *tenor, spot.spot_lag, fx, &calendar)?;
                 let settlement = refdata
-                    .zone(&fx.settlement.zone)?
-                    .resolve_convention(fx.settlement.on_date(period.end()))?;
+                    .zone(&spot.settlement.zone)?
+                    .resolve_convention(spot.settlement.on_date(period.end()))?;
                 Ok(StandardProduct::FxForward(FxForward::new(
                     Money::new(fx.pair.base(), 1.0),
                     FxRate::forward(fx.pair, value, settlement)?,
@@ -288,6 +290,10 @@ impl InstrumentTemplate {
     /// # Errors
     ///
     /// テンプレートに合う慣行がないか、慣行に校正できない丸めやスプレッドの規約がある場合。
+    #[expect(
+        clippy::wildcard_enum_match_arm,
+        reason = "#199 で各 variant を網羅する"
+    )]
     pub fn currency(&self, refdata: &RefData) -> Result<Currency> {
         let id = self.convention();
         let convention = refdata.convention(id).map_err(|error| {
@@ -377,13 +383,14 @@ fn tenor_offset(length: TenorLength) -> Result<DateOffset> {
 fn fx_dates(
     trade_date: Date,
     tenor: Tenor,
+    spot_lag: rhodium_core::refdata::FxSpotLag,
     fx: &rhodium_core::refdata::FxForwardConvention,
     calendar: &rhodium_core::calendar::Calendar,
 ) -> Result<(rhodium_core::date::DateInterval, Vec<rhodium_core::Warning>)> {
     money_market_period(
         trade_date,
         tenor,
-        fx.spot_lag,
+        spot_lag.offset(),
         &AdvanceRule {
             calendar,
             bdc: fx.bdc,
diff --git a/crates/rhodium-products/src/trade.rs b/crates/rhodium-products/src/trade.rs
index 056bc8a2fda38a58b02ff223b634ff7255a08006..026e002b9d0e96864ce19f7c594c37813d130001 100644
--- a/crates/rhodium-products/src/trade.rs
+++ b/crates/rhodium-products/src/trade.rs
@@ -95,17 +95,66 @@ impl<'r> TradeProductResolver<'r> {
         }
     }
 
-    /// 取引を検証済み IR にする。複数のワーカーから呼び出せる。
+    /// 取引を検証済み IR に下ろし、取引と組にする。複数のワーカーから呼び出せる。
     ///
     /// # Errors
     ///
     /// [`Trade::to_product`] と同じ失敗。
-    pub fn resolve(&self, trade: &Trade) -> Result<Product> {
-        trade.to_product_using(self.refdata, Some(&self.observations))
+    pub fn resolve(&self, trade: Trade) -> Result<ResolvedTrade> {
+        let product = trade.to_product_using(self.refdata, Some(&self.observations))?;
+        Ok(ResolvedTrade { trade, product })
+    }
+}
+
+/// 取引と、同じ参照データでその取引を下ろした検証済み IR の組。
+///
+/// 評価の方式によって読む側が変わる（解析プライサーは取引の標準商品・約定日・手数料を、
+/// リニア評価は IR を読む）ので、無関係な取引と IR を組むと方式を切り替えただけで別の取引を
+/// 評価してしまう。フィールドを公開せず、[`TradeProductResolver::resolve`] と
+/// [`ResolvedTrade::new`] だけが作ることで、二つが同じ取引から来ることを保証する。
+#[derive(Debug, Clone)]
+pub struct ResolvedTrade {
+    trade: Trade,
+    product: Product,
+}
+
+impl ResolvedTrade {
+    /// 取引を `refdata` で検証済み IR に下ろし、取引と組にする。
+    ///
+    /// # Errors
+    ///
+    /// [`Trade::to_product`] と同じ失敗。
+    pub fn new(trade: Trade, refdata: &RefData) -> Result<Self> {
+        let product = trade.to_product(refdata)?;
+        Ok(Self { trade, product })
+    }
+
+    /// 取引の条件。
+    #[must_use]
+    pub const fn trade(&self) -> &Trade {
+        &self.trade
+    }
+
+    /// 取引を下ろした検証済み IR。
+    #[must_use]
+    pub const fn product(&self) -> &Product {
+        &self.product
     }
 }
 
 impl Trade {
+    /// 価格設定の方式の表が照合する商品の種類。
+    ///
+    /// 型付きの標準商品だけが種類を持つ。builder の契約と変更した商品（一部解約、置き換え）は
+    /// 標準商品の条件と一致しないので種類を持たず、方式の表の行に当たらない。
+    #[must_use]
+    pub fn product_kind(&self) -> Option<crate::standard::ProductKind> {
+        match &self.product {
+            ProductSpec::Standard(product) => Some(product.kind()),
+            ProductSpec::Contract(_) | ProductSpec::Amended { .. } => None,
+        }
+    }
+
     /// 商品と手数料を検証済み IR にし、評価に使う行使記録も格納する。
     ///
     /// # Errors
diff --git a/crates/rhodium-products/tests/compile_fail.rs b/crates/rhodium-products/tests/compile_fail.rs
index 6b951ce71dcfb8c0d6e97c8dcd502bfbceb7069c..6e9b93c6e807858e131c133952849fc27200bb29 100644
--- a/crates/rhodium-products/tests/compile_fail.rs
+++ b/crates/rhodium-products/tests/compile_fail.rs
@@ -19,3 +19,9 @@ fn boolean_expression_cannot_be_paid() {
     // 真偽値の式は支払額に使えない。
     compile_fail_pair!("boolean_payment");
 }
+
+#[test]
+fn resolved_trade_pairs_only_a_trade_with_its_own_product() {
+    // 取引と下ろした IR の組は、取引を下ろす関数だけが作る。
+    compile_fail_pair!("resolved_trade_from_parts");
+}
diff --git a/crates/rhodium-products/tests/ui/boolean_payment.json b/crates/rhodium-products/tests/ui/boolean_payment.json
new file mode 100644
index 0000000000000000000000000000000000000000..4f10e207db4008891bc2418774141c8299631a11
--- /dev/null
+++ b/crates/rhodium-products/tests/ui/boolean_payment.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0308",
+    "message": "mismatched types",
+    "spans": [
+      {
+        "file": "tests/ui/boolean_payment.rs",
+        "line": 11,
+        "column": 52
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-products/tests/ui/boolean_payment.stderr b/crates/rhodium-products/tests/ui/boolean_payment.stderr
deleted file mode 100644
index fe2dcc986d481cbdd521308b69dae589f380ddc1..0000000000000000000000000000000000000000
--- a/crates/rhodium-products/tests/ui/boolean_payment.stderr
+++ /dev/null
@@ -1,15 +0,0 @@
-error[E0308]: mismatched types
-  --> tests/ui/boolean_payment.rs:11:52
-   |
-11 |         Ok(builder.pay(Asset::Cash(Currency::JPY), amount, at))
-   |                    ---                             ^^^^^^ expected `ExprIx<'_>`, found `ExprIx<'_, Bool>`
-   |                    |
-   |                    arguments to this method are incorrect
-   |
-   = note: expected struct `ExprIx<'_, Num>`
-              found struct `ExprIx<'_, rhodium_products::ir::Bool>`
-note: method defined here
-  --> src/ir/builder.rs
-   |
-   |     pub fn pay(&mut self, asset: Asset, amount: ExprIx<'b>, at: Instant) -> ContractIx<'b> {
-   |            ^^^
diff --git a/crates/rhodium-products/tests/ui/cross_builder_contract.json b/crates/rhodium-products/tests/ui/cross_builder_contract.json
new file mode 100644
index 0000000000000000000000000000000000000000..1188123b605aa3381d3c49b79200c66cee0d595b
--- /dev/null
+++ b/crates/rhodium-products/tests/ui/cross_builder_contract.json
@@ -0,0 +1,29 @@
+[
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/cross_builder_contract.rs",
+        "line": 7,
+        "column": 44
+      }
+    ]
+  },
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/cross_builder_contract.rs",
+        "line": 7,
+        "column": 44
+      },
+      {
+        "file": "tests/ui/cross_builder_contract.rs",
+        "line": 7,
+        "column": 44
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-products/tests/ui/cross_builder_contract.stderr b/crates/rhodium-products/tests/ui/cross_builder_contract.stderr
deleted file mode 100644
index 66716f225538f766ce7ec42fce67ad45ef14bb78..0000000000000000000000000000000000000000
--- a/crates/rhodium-products/tests/ui/cross_builder_contract.stderr
+++ /dev/null
@@ -1,28 +0,0 @@
-error[E0521]: borrowed data escapes outside of closure
- --> tests/ui/cross_builder_contract.rs:7:44
-  |
-6 |         let child = outer.zero();
-  |             ----- `child` declared here, outside of the closure body
-7 |         let _ = Contract::build(|inner| Ok(inner.all([child])))?;
-  |                                  -----     ^^^^^^^^^^^^^^^^^^ `inner` escapes the closure body here
-  |                                  |
-  |                                  `inner` is only valid in the closure body
-  |
-  = note: requirement occurs because of the type `ContractIx<'_>`, which makes the generic argument `'_` invariant
-  = note: the struct `ContractIx<'b>` is invariant over the parameter `'b`
-  = help: see <https://doc.rust-lang.org/nomicon/subtyping.html> for more information about variance
-
-error[E0521]: borrowed data escapes outside of closure
- --> tests/ui/cross_builder_contract.rs:7:44
-  |
-5 |     let _ = Contract::build(|outer| {
-  |                              -----
-  |                              |
-  |                              `outer` is only valid in the closure body
-  |                              has type `&mut rhodium_products::ir::Builder<'1>`
-6 |         let child = outer.zero();
-7 |         let _ = Contract::build(|inner| Ok(inner.all([child])))?;
-  |                                            ^^^^^^^^^^^^^^^^^^
-  |                                            |
-  |                                            `outer` escapes the closure body here
-  |                                            argument requires that `'1` must outlive `'static`
diff --git a/crates/rhodium-products/tests/ui/cross_builder_expression.json b/crates/rhodium-products/tests/ui/cross_builder_expression.json
new file mode 100644
index 0000000000000000000000000000000000000000..d8b7c7cf9c85f5cfafe07656c92636299e9f73cb
--- /dev/null
+++ b/crates/rhodium-products/tests/ui/cross_builder_expression.json
@@ -0,0 +1,29 @@
+[
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/cross_builder_expression.rs",
+        "line": 12,
+        "column": 25
+      }
+    ]
+  },
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/cross_builder_expression.rs",
+        "line": 12,
+        "column": 25
+      },
+      {
+        "file": "tests/ui/cross_builder_expression.rs",
+        "line": 12,
+        "column": 25
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-products/tests/ui/cross_builder_expression.stderr b/crates/rhodium-products/tests/ui/cross_builder_expression.stderr
deleted file mode 100644
index 144498c4021545f42629b80758ec87eafa232a40..0000000000000000000000000000000000000000
--- a/crates/rhodium-products/tests/ui/cross_builder_expression.stderr
+++ /dev/null
@@ -1,28 +0,0 @@
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/cross_builder_expression.rs:12:25
-   |
-10 |         let amount = outer.constant(42.0);
-   |             ------ `amount` declared here, outside of the closure body
-11 |         let _ = Contract::build(|inner| {
-   |                                  ----- `inner` is only valid in the closure body
-12 |             let local = inner.constant(1.0);
-   |                         ^^^^^^^^^^^^^^^^^^^ `inner` escapes the closure body here
-   |
-   = note: requirement occurs because of the type `ExprIx<'_>`, which makes the generic argument `'_` invariant
-   = note: the struct `ExprIx<'b, T>` is invariant over the parameter `'b`
-   = help: see <https://doc.rust-lang.org/nomicon/subtyping.html> for more information about variance
-
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/cross_builder_expression.rs:12:25
-   |
- 9 |     let _ = Contract::build(|outer| {
-   |                              -----
-   |                              |
-   |                              `outer` is only valid in the closure body
-   |                              has type `&mut rhodium_products::ir::Builder<'1>`
-...
-12 |             let local = inner.constant(1.0);
-   |                         ^^^^^^^^^^^^^^^^^^^
-   |                         |
-   |                         `outer` escapes the closure body here
-   |                         argument requires that `'1` must outlive `'static`
diff --git a/crates/rhodium-products/tests/ui/resolved_trade_from_parts.json b/crates/rhodium-products/tests/ui/resolved_trade_from_parts.json
new file mode 100644
index 0000000000000000000000000000000000000000..f2fc48a7d5c9fb1100a27f4ef08841425f566105
--- /dev/null
+++ b/crates/rhodium-products/tests/ui/resolved_trade_from_parts.json
@@ -0,0 +1,18 @@
+[
+  {
+    "code": "E0451",
+    "message": "fields `trade` and `product` of struct `ResolvedTrade` are private",
+    "spans": [
+      {
+        "file": "tests/ui/resolved_trade_from_parts.rs",
+        "line": 9,
+        "column": 21
+      },
+      {
+        "file": "tests/ui/resolved_trade_from_parts.rs",
+        "line": 9,
+        "column": 28
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-products/tests/ui/resolved_trade_from_parts.rs b/crates/rhodium-products/tests/ui/resolved_trade_from_parts.rs
new file mode 100644
index 0000000000000000000000000000000000000000..6c0d3ca6307ccd849d67c7dc148bafdf0cd71fd2
--- /dev/null
+++ b/crates/rhodium-products/tests/ui/resolved_trade_from_parts.rs
@@ -0,0 +1,14 @@
+use rhodium_products::{
+    ir::Product,
+    trade::{ResolvedTrade, Trade},
+};
+
+// 取引と IR を別々に渡して組を作ることはできない。解析は取引を、リニアは IR を読むので、
+// 無関係な二つを組むと方式を切り替えただけで別の取引を評価してしまう。
+fn pair(trade: Trade, product: Product) -> ResolvedTrade {
+    ResolvedTrade { trade, product }
+}
+
+fn main() {
+    let _ = pair;
+}
diff --git a/crates/rhodium-products/tests/ui/resolved_trade_from_parts_ok.rs b/crates/rhodium-products/tests/ui/resolved_trade_from_parts_ok.rs
new file mode 100644
index 0000000000000000000000000000000000000000..396815e789edb0c9ba663defc9780871a6358a04
--- /dev/null
+++ b/crates/rhodium-products/tests/ui/resolved_trade_from_parts_ok.rs
@@ -0,0 +1,30 @@
+use std::collections::BTreeMap;
+
+use rhodium_core::{
+    Result,
+    refdata::{RefData, RefDataParts},
+};
+use rhodium_products::{
+    ir::Contract,
+    trade::{FallbackTerms, ProductSpec, ResolvedTrade, Trade},
+};
+
+fn main() -> Result<()> {
+    // 組は取引を下ろして作るので、IR は必ずその取引から来る。
+    let refdata = RefData::new(RefDataParts::default())?;
+    let trade = Trade {
+        id: "T1".parse()?,
+        product: ProductSpec::Contract(Contract::build(|builder| Ok(builder.zero()))?),
+        margin_set: "MS".parse()?,
+        book: "BOOK".parse()?,
+        trade_date: "2027-01-04".parse()?,
+        attributes: BTreeMap::new(),
+        exercises: Vec::new(),
+        fees: Vec::new(),
+        fallback: FallbackTerms::default(),
+    };
+    let resolved = ResolvedTrade::new(trade, &refdata)?;
+    assert_eq!(resolved.trade().id.as_str(), "T1");
+    let _ = resolved.product();
+    Ok(())
+}
diff --git a/crates/rhodium-risk/tests/least_squares.rs b/crates/rhodium-risk/tests/least_squares.rs
new file mode 100644
index 0000000000000000000000000000000000000000..7c790fd161903b240aa3491dd0aa8a168d388d62
--- /dev/null
+++ b/crates/rhodium-risk/tests/least_squares.rs
@@ -0,0 +1,340 @@
+//! 最小二乗の校正グループを通る quote デルタを、手計算と再校正の中心差分で照合する。
+
+use std::collections::BTreeMap;
+
+use rhodium_core::{
+    Error, Result,
+    currency::{Asset, Currency},
+    date::Date,
+    money::Money,
+};
+use rhodium_market::{
+    curve::{CurveId, CurveRole, CurveSpec, Extrapolation},
+    quote::{QuoteId, QuoteKind, QuoteSet},
+    structure::{BlockId, ParamId},
+};
+use rhodium_pricing::{
+    calibration::{CalibratedMarket, GroupOutcome, calibrate, calibrate_with_residual_tolerance},
+    context::Context,
+    definition::{
+        CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
+    },
+    prepare::{PrepareInput, Prepared, evaluate, prepare},
+    resolve::CollateralKey,
+};
+use rhodium_products::{
+    ir::{Contract, Product},
+    standard::{Deposit, DepositDates, PaymentRounding, Position, StandardProduct},
+};
+use rhodium_risk::{delta::quote_delta, gradient::param_gradient, vector::RiskVector};
+use rhodium_testkit::{assert_close, check_aad_vs_bump};
+use scioxide_core::ExecutionContext;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、最小二乗の照合に使わない関数がある"
+)]
+#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
+pub mod ois;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、最小二乗の照合に使わない関数がある"
+)]
+#[path = "../../rhodium-pricing/tests/fixtures/steps.rs"]
+pub mod steps;
+
+use ois::{Inputs, Node, cases, midnight, prepare_ois};
+use steps::{PenalizedTona, TURN_QUOTE};
+
+const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
+const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
+const NOTIONAL: f64 = 1_000_000.0;
+const BP: f64 = 1e-4;
+
+fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
+    vector
+        .entries()
+        .binary_search_by(|(key, _)| key.cmp(id))
+        .ok()
+        .map_or(0.0, |index| vector.entries()[index].1)
+}
+
+fn jpy(report: &rhodium_risk::delta::DeltaReport) -> Result<&RiskVector<QuoteId>> {
+    report
+        .by_ccy()
+        .iter()
+        .find(|vector| vector.ccy() == Currency::JPY)
+        .ok_or_else(|| Error::Internal("no JPY quote deltas".into()))
+}
+
+#[test]
+fn lsq_adjoint_oracle() -> Result<()> {
+    // 感応度の検証表の例: min (θ − q₁)² + 4(θ − q₂)² の解 θ = (q₁ + 4q₂)/5 で、
+    // quote への勾配 (1/5, 4/5) を校正の随伴が再現する。θ はピラーの ln DF なので、
+    // ピラーで 100 万円を受け取る取引の V = N·exp(θ) の θ 勾配 V を、二つのマークに 1:4 で配る。
+    let [tona, _] = cases()?;
+    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
+    let ctx = inputs.context();
+    let id: CurveId = tona.index.parse()?;
+    let pillar = midnight(tona.nodes[2].date)?;
+    let curve = CurveSpec::new(
+        id.clone(),
+        Currency::JPY,
+        [CurveRole::Discount, CurveRole::Projection],
+        [pillar],
+        Extrapolation::Reject,
+    )?;
+    let marks: [QuoteId; 2] = ["MARK.1".parse()?, "MARK.2".parse()?];
+    let bindings = marks
+        .iter()
+        .zip([1.0, 0.5])
+        .map(|(quote, sigma)| {
+            InstrumentBinding::new(
+                quote.clone(),
+                QuoteKind::Mark,
+                CalibrationTarget::DirectMark,
+                QuoteTransform::Identity,
+                COLLATERAL,
+                ElementSelector::Param(ParamId::new(BlockId::nodes(id.clone()), 0)),
+            )?
+            .with_tolerance(sigma)
+        })
+        .collect::<Result<Vec<_>>>()?;
+    let definition = MarketDefinition::new("JPY".parse()?, "Etc/UTC".parse()?, [curve], bindings)?;
+    let quotes = QuoteSet::new([(marks[0].clone(), -0.01), (marks[1].clone(), -0.02)])?;
+    let market = calibrate(&definition, &quotes, &ctx)?;
+    let product = Product::new(
+        Contract::build(|builder| {
+            let amount = builder.constant(NOTIONAL);
+            Ok(builder.pay(Asset::Cash(Currency::JPY), amount, pillar))
+        })?,
+        ctx.refdata,
+    )?;
+    let trade = prepare(
+        PrepareInput::Ir(&product),
+        &COLLATERAL,
+        &[],
+        market.market().structure(),
+        &ctx,
+    )?;
+    let gradient = param_gradient(&trade, &market)?;
+    let [theta_gradient] = gradient.by_ccy()[0].entries() else {
+        return Err(Error::Internal("expected one theta gradient".into()));
+    };
+    // V = N·exp(θ) なので ∂V/∂θ = V。
+    let pv = evaluate(&trade, market.market())?.get(Currency::JPY);
+    assert_close!(theta_gradient.1, pv, 1e-12, 1e-6);
+    let report = quote_delta(std::slice::from_ref(&trade), &market, SEQUENTIAL)?;
+    let vector = jpy(&report)?;
+    // マークの感応度の単位は自然な単位 1。
+    assert_close!(sensitivity(vector, &marks[0]), pv / 5.0, 1e-10, 1e-6);
+    assert_close!(sensitivity(vector, &marks[1]), 4.0 * pv / 5.0, 1e-10, 1e-6);
+    Ok(())
+}
+
+/// バンプの再校正の許容誤差。最小二乗のグループは停留条件の残差を θ に戻した各行の動きで判定する。
+const BUMP_TOLERANCE: f64 = 1e-15;
+
+/// 校正の停止の誤差を一つの PV に伝えた幅 `Σ_j |∂V/∂θ_j|·|δθ_j|`（品質の「ソルバーの誤差」）。
+///
+/// δθ = (∂G/∂θ)⁻¹·G は各グループの `theta_error`（列は θ の添字の昇順）。最小二乗の罰則だけで
+/// 決まる向きには対応する市場の quote が無いので、quote の感応度の和では誤差を覆えない。
+fn solver_pv_error(trade: &Prepared, market: &CalibratedMarket) -> Result<f64> {
+    let gradient = param_gradient(trade, market)?;
+    let slopes: BTreeMap<&ParamId, f64> = gradient
+        .by_ccy()
+        .iter()
+        .flat_map(|vector| {
+            vector
+                .entries()
+                .iter()
+                .map(|(param, slope)| (param, *slope))
+        })
+        .collect();
+    let layout = market.market().structure().layout();
+    market
+        .groups()
+        .iter()
+        .map(|group| {
+            let GroupOutcome::Solved(solution) = group.outcome() else {
+                return Err(Error::Calibration(format!("{group:?}")));
+            };
+            let mut params = Vec::new();
+            for block in group.blocks() {
+                let range = layout
+                    .range(block)
+                    .ok_or_else(|| Error::Internal(format!("missing block {block:?}")))?;
+                for (element, index) in range.enumerate() {
+                    let element = u32::try_from(element)
+                        .map_err(|error| Error::Internal(error.to_string()))?;
+                    params.push((index, ParamId::new(block.clone(), element)));
+                }
+            }
+            params.sort_by_key(|&(index, _)| index);
+            Ok(params
+                .iter()
+                .zip(solution.theta_error())
+                .map(|((_, param), error)| {
+                    slopes.get(param).map_or(0.0, |slope| slope.abs()) * error.abs()
+                })
+                .sum::<f64>())
+        })
+        .sum()
+}
+
+fn calibrated(
+    definition: &MarketDefinition,
+    quotes: &QuoteSet,
+    ctx: &Context<'_>,
+    tolerance: f64,
+) -> Result<CalibratedMarket> {
+    let market = calibrate_with_residual_tolerance(definition, quotes, ctx, tolerance)?;
+    for group in market.groups() {
+        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
+            return Err(Error::Calibration(format!("{group:?}")));
+        }
+    }
+    Ok(market)
+}
+
+fn bumped(quotes: &QuoteSet, id: &QuoteId, bump: f64) -> Result<QuoteSet> {
+    QuoteSet::new(
+        quotes
+            .iter()
+            .map(|(key, value)| (key.clone(), if key == id { value + bump } else { value })),
+    )
+}
+
+/// 階段の区間をまたぐ 9M と、接続点をまたぐ 18M の OIS（quote から外した固定金利）と、
+/// 年末をまたぐ預金。
+fn trades(
+    tona: &PenalizedTona,
+    market: &CalibratedMarket,
+    ctx: &Context<'_>,
+) -> Result<Vec<Prepared>> {
+    let node = |tenor: &str| -> Result<Node> {
+        Ok(Node {
+            tenor: tenor.parse()?,
+            rate: 0.0,
+            date: tona.case.as_of,
+            df: 1.0,
+        })
+    };
+    let mut prepared = [("9M", 0.003), ("18M", 0.0065)]
+        .into_iter()
+        .map(|(tenor, rate)| {
+            prepare_ois(
+                &tona.case,
+                &node(tenor)?,
+                rate,
+                &COLLATERAL,
+                market.market(),
+                ctx,
+            )
+        })
+        .collect::<Result<Vec<_>>>()?;
+    let start: Date = "2025-07-09".parse()?;
+    let deposit = StandardProduct::Deposit(Deposit::new(
+        Money::new(Currency::JPY, NOTIONAL),
+        0.004,
+        DepositDates::Confirmed {
+            start,
+            maturity: "2026-04-09".parse()?,
+        },
+        None,
+        Some(tona.case.basis()),
+        Position::Long,
+        PaymentRounding::None,
+    )?);
+    prepared.push(prepare(
+        PrepareInput::Standard {
+            product: &deposit,
+            trade_date: tona.case.as_of,
+        },
+        &COLLATERAL,
+        &[],
+        market.market().structure(),
+        ctx,
+    )?);
+    Ok(prepared)
+}
+
+#[test]
+fn penalized_meeting_steps_quote_delta_matches_bump() -> Result<()> {
+    // R4 の網羅表の行: 罰則付きの最小二乗で当てた会合日の階段を通る quote デルタ（AAD と
+    // 2 階の項を含む校正の随伴）が、quote ごとに全体を再校正した中心差分と一致する。
+    // 残差が残る当てはめなので、随伴が残差の 2 階の項を落とすと一致しない（Gauss–Newton の記録）。
+    let tona = PenalizedTona::new()?;
+    let inputs = Inputs::new(&[&tona.case], tona.case.as_of, "")?;
+    let ctx = inputs.context();
+    let definition = tona.definition(false)?;
+    let quotes = tona.quotes()?;
+    let base = calibrated(&definition, &quotes, &ctx, BUMP_TOLERANCE)?;
+    let ids: Vec<QuoteId> = tona
+        .case
+        .nodes
+        .iter()
+        .map(|node| tona.case.quote(node))
+        .chain([TURN_QUOTE.parse()])
+        .collect::<Result<_>>()?;
+    let prepared = trades(&tona, &base, &ctx)?;
+    for (index, trade) in prepared.iter().enumerate() {
+        let report = quote_delta(std::slice::from_ref(trade), &base, SEQUENTIAL)?;
+        let vector = jpy(&report)?;
+        for id in &ids {
+            // 中心差分が使う ±BP と ±BP/2 の再校正ごとに PV と、停止の誤差を θ 勾配で伝えた幅を
+            // 求め、大きい方の幅を一つの PV の誤差にする（品質の c_s = 2 は照合の関数が掛ける）。
+            let repriced = [BP, -BP, BP / 2.0, -BP / 2.0]
+                .into_iter()
+                .map(|h| {
+                    let market =
+                        calibrated(&definition, &bumped(&quotes, id, h)?, &ctx, BUMP_TOLERANCE)?;
+                    let trade = &trades(&tona, &market, &ctx)?[index];
+                    Ok((
+                        h,
+                        evaluate(trade, market.market())?.get(Currency::JPY),
+                        solver_pv_error(trade, &market)?,
+                    ))
+                })
+                .collect::<Result<Vec<_>>>()?;
+            let pv_error = repriced
+                .iter()
+                .map(|&(_, _, error)| error)
+                .fold(0.0, f64::max);
+            check_aad_vs_bump(
+                |h| {
+                    repriced
+                        .iter()
+                        // 照合の関数は同じ式（BP と BP/2）で幅を作るので、ビットで一致する。
+                        .find(|&&(bump, _, _)| bump.to_bits() == h.to_bits())
+                        .map(|&(_, pv, _)| pv)
+                        .ok_or_else(|| Error::Internal(format!("unexpected bump {h}")))
+                },
+                sensitivity(vector, id),
+                BP,
+                pv_error,
+                1e-6,
+                &format!("trade {index} {id}"),
+            )?;
+        }
+    }
+    // 対照: Gauss–Newton の記録の随伴は残差の 2 階の項を落とすので、同じ取引のデルタが厳密な値から
+    // 比べる許容誤差（相対 1e-10）よりはるかに大きくずれる。上の照合がこの差を検出できることを示す。
+    let approximate = calibrated(&tona.definition(true)?, &quotes, &ctx, BUMP_TOLERANCE)?;
+    let trade = &trades(&tona, &approximate, &ctx)?[0];
+    let exact = quote_delta(std::slice::from_ref(&prepared[0]), &base, SEQUENTIAL)?;
+    let gauss_newton = quote_delta(std::slice::from_ref(trade), &approximate, SEQUENTIAL)?;
+    let deviation = ids
+        .iter()
+        .map(|id| {
+            let lhs = sensitivity(jpy(&exact)?, id);
+            let rhs = sensitivity(jpy(&gauss_newton)?, id);
+            Ok::<_, Error>((lhs - rhs).abs() / lhs.abs().max(1e-6))
+        })
+        .collect::<Result<Vec<_>>>()?
+        .into_iter()
+        .fold(0.0_f64, f64::max);
+    assert!(deviation > 1e-6, "{deviation:e}");
+    Ok(())
+}
diff --git a/crates/rhodium-risk/tests/m1_quote_delta.rs b/crates/rhodium-risk/tests/m1_quote_delta.rs
new file mode 100644
index 0000000000000000000000000000000000000000..3b76c7593d2a7bb845ac472742bc63fd2efd0b4b
--- /dev/null
+++ b/crates/rhodium-risk/tests/m1_quote_delta.rs
@@ -0,0 +1,434 @@
+//! AAD とバンプの網羅表の M1 の行: 商品の族と曲線の組み方ごとに、quote デルタ（AAD と校正の随伴）を
+//! quote ごとに全体を再校正した中心差分と照合する。
+//!
+//! - R2・R5: OIS の割引とターム指標の予想（TONA へのベーシスの TIBOR 6M・3M、USD の独立したターム曲線）
+//!   の上の IRS、FRA、ベーシススワップ。
+//! - R1 と LGM のマーク: 金利先物（解析プライサー、凸性の補正）の先物価格と LGM のマークへの感応度。
+//! - R3: 会合日の階段と年末の跳びを持つ TONA の上の OIS（lookback、observation shift、lockout、
+//!   支払いの遅れ）と、年末をまたぐ預金。
+//!
+//! 再校正した市場は quote とバンプの組ごとに一度だけ作り、同じ市場で全取引を評価する。
+
+use std::{
+    cell::RefCell,
+    collections::HashMap,
+    num::{NonZeroU8, NonZeroU32},
+};
+
+use rhodium_core::{
+    Error, Result,
+    calendar::CalendarSpec,
+    currency::Currency,
+    date::{Date, DateInterval},
+    day_count::DayCount,
+    fixing::FixingStore,
+    money::Money,
+    observation::{
+        Aggregation, CompoundingForm, ObservationDays, ObservationMethod, ObservationRounding,
+        ObservationRule, SpreadTreatment,
+    },
+    quantity::{Rate, RateFuturesPrice, Simple, Spread},
+    schedule::{
+        ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination, ScheduleSpec,
+    },
+};
+use rhodium_market::quote::{QuoteId, QuoteKind, QuoteSet};
+use rhodium_pricing::{
+    calibration::{CalibratedMarket, GroupOutcome, calibrate_with_residual_tolerance},
+    context::Context,
+    definition::MarketDefinition,
+    prepare::{PrepareInput, Prepared, evaluate, prepare},
+    resolve::CollateralKey,
+};
+use rhodium_products::standard::{
+    Coupon, Deposit, DepositDates, Fra, Leg, NotionalExchange, PayReceive, PaymentRounding,
+    Position, RateFuture, StandardProduct, Swap,
+};
+use rhodium_risk::{delta::quote_delta, vector::RiskVector};
+use rhodium_testkit::check_aad_vs_bump;
+use scioxide_core::ExecutionContext;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
+)]
+#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
+pub mod ois;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
+)]
+#[path = "../../rhodium-pricing/tests/fixtures/multicurve.rs"]
+pub mod multicurve;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
+)]
+#[path = "../../rhodium-pricing/tests/fixtures/steps.rs"]
+pub mod steps;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
+)]
+#[path = "../../rhodium-pricing/tests/fixtures/futures.rs"]
+pub mod futures;
+
+use multicurve::{MulticurveCase, all_multicurve_cases};
+use ois::{Inputs, Node, prepare_ois};
+use steps::PenalizedTona;
+
+const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
+const BP: f64 = 1e-4;
+const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
+/// バンプの再校正の残差の許容誤差（通常の 1e-14 の 1/100）。
+const BUMP_TOLERANCE: f64 = 1e-16;
+
+/// 市場から照合する取引を準備する関数。取引の名前と準備の組を返す。
+type Trades<'a> = dyn Fn(&CalibratedMarket) -> Result<Vec<(String, Prepared)>> + 'a;
+
+fn calibrated(
+    definition: &MarketDefinition,
+    quotes: &QuoteSet,
+    ctx: &Context<'_>,
+) -> Result<CalibratedMarket> {
+    let market = calibrate_with_residual_tolerance(definition, quotes, ctx, BUMP_TOLERANCE)?;
+    for group in market.groups() {
+        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
+            return Err(Error::Calibration(format!("{group:?}")));
+        }
+    }
+    Ok(market)
+}
+
+fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
+    vector
+        .entries()
+        .binary_search_by(|(key, _)| key.cmp(id))
+        .ok()
+        .map_or(0.0, |index| vector.entries()[index].1)
+}
+
+/// quote の報告単位（デルタの値の単位）と、バンプの幅。マークは報告単位 1 だが、1 だけ動かすと
+/// 非線形の項が大きいので、1bp の幅で照合し、AAD を幅あたりに換算する。
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#203 で各 variant を網羅する"
+)]
+fn units(market: &CalibratedMarket, id: &QuoteId) -> Result<(f64, f64)> {
+    Ok(
+        match market
+            .quote_kind(id)
+            .ok_or_else(|| Error::Internal(format!("unknown kind of {id}")))?
+        {
+            QuoteKind::Mark => (1.0, BP),
+            _ => (BP, BP),
+        },
+    )
+}
+
+/// 全取引と全 quote で、AAD の quote デルタを再校正の中心差分と照合する。
+fn check_all(
+    definition: &MarketDefinition,
+    quotes: &QuoteSet,
+    ctx: &Context<'_>,
+    ccy: Currency,
+    floor: f64,
+    trades: &Trades<'_>,
+) -> Result<()> {
+    let base = calibrated(definition, quotes, ctx)?;
+    // 同じ (quote, 幅) の再校正は取引の間で共有する。
+    let cache: RefCell<HashMap<(QuoteId, u64), CalibratedMarket>> = RefCell::new(HashMap::new());
+    let reprice = |id: &QuoteId, bump: f64, index: usize| -> Result<f64> {
+        let key = (id.clone(), bump.to_bits());
+        if !cache.borrow().contains_key(&key) {
+            let bumped =
+                QuoteSet::new(quotes.iter().map(|(key, value)| {
+                    (key.clone(), if key == id { value + bump } else { value })
+                }))?;
+            let market = calibrated(definition, &bumped, ctx)?;
+            cache.borrow_mut().insert(key.clone(), market);
+        }
+        let cache = cache.borrow();
+        let market = &cache[&key];
+        let (_, prepared) = &trades(market)?[index];
+        Ok(evaluate(prepared, market.market())?.get(ccy))
+    };
+    for (index, (name, prepared)) in trades(&base)?.iter().enumerate() {
+        let report = quote_delta(std::slice::from_ref(prepared), &base, SEQUENTIAL)?;
+        assert!(report.failures().is_empty(), "{name}");
+        let vector = report
+            .by_ccy()
+            .iter()
+            .find(|vector| vector.ccy() == ccy)
+            .ok_or_else(|| Error::Internal(format!("{name}: no {ccy} deltas")))?;
+        // 再校正の残差の許容誤差を、quote ごとの感応度（自然な単位）の和で PV の誤差に伝える。
+        let scale = vector
+            .entries()
+            .iter()
+            .map(|(id, delta)| Ok(delta.abs() / units(&base, id)?.0))
+            .sum::<Result<f64>>()?;
+        assert!(scale > 0.0, "{name} has no quote deltas");
+        for (id, _) in quotes.iter() {
+            let (report_unit, width) = units(&base, id)?;
+            check_aad_vs_bump(
+                |h| reprice(id, h, index),
+                sensitivity(vector, id) * width / report_unit,
+                width,
+                BUMP_TOLERANCE * scale,
+                floor,
+                &format!("{name} {id}"),
+            )?;
+        }
+    }
+    Ok(())
+}
+
+fn standard(
+    product: &StandardProduct,
+    trade_date: Date,
+    market: &CalibratedMarket,
+    ctx: &Context<'_>,
+) -> Result<Prepared> {
+    prepare(
+        PrepareInput::Standard {
+            product,
+            trade_date,
+        },
+        &COLLATERAL,
+        &[],
+        market.market().structure(),
+        ctx,
+    )
+}
+
+/// 固定 1% の IRS、ターム指標の FRA、quote から 10bp 外したベーシススワップ。元本は 1。
+fn term_products(
+    case: &MulticurveCase,
+    ctx: &Context<'_>,
+) -> Result<Vec<(String, StandardProduct)>> {
+    let node = &case.nodes[2];
+    let irs = rhodium_products::standard::InstrumentTemplate::Swap {
+        convention: case.irs_convention().parse()?,
+        tenor: node.tenor,
+    }
+    .build(0.01, ctx.refdata)?;
+    let fra = StandardProduct::Fra(Fra::new(
+        Money::new(case.overnight.currency, 1.0),
+        case.term_index.parse()?,
+        "2025-07-07".parse()?,
+        0.01,
+        Position::Long,
+    )?);
+    let basis = case.template(node)?.build(node.rate + 0.001, ctx.refdata)?;
+    Ok(vec![
+        (format!("{} IRS", case.term_index), irs),
+        (format!("{} FRA", case.term_index), fra),
+        (format!("{} basis", case.term_index), basis),
+    ])
+}
+
+#[test]
+fn term_products_quote_deltas_match_bumps_on_term_curves() -> Result<()> {
+    // R2・R5: TIBOR 6M と 3M は TONA の曲線とベーシスの和、USD のターム曲線は独立の節点。
+    // どれも OIS の曲線で割り引き、ターム曲線のリスクは随伴でベーシスと OIS の quote に戻る。
+    for case in all_multicurve_cases()? {
+        let inputs = case.inputs()?;
+        let ctx = inputs.context();
+        let products = term_products(&case, &ctx)?;
+        let trade_date = case.overnight.as_of;
+        let trades = |market: &CalibratedMarket| {
+            products
+                .iter()
+                .map(|(name, product)| {
+                    Ok((name.clone(), standard(product, trade_date, market, &ctx)?))
+                })
+                .collect()
+        };
+        check_all(
+            &case.definition()?,
+            &case.quotes()?,
+            &ctx,
+            case.overnight.currency,
+            1e-14,
+            &trades,
+        )?;
+    }
+    Ok(())
+}
+
+#[test]
+fn rate_future_quote_deltas_match_bumps() -> Result<()> {
+    // 金利先物は解析プライサー（LGM の凸性の補正）で評価する。曲線を先物の価格で校正し、
+    // 凸性に効く LGM の平均回帰とボラは直接マーク。先物価格への感応度は曲線の随伴で、
+    // マークへの感応度は凸性の補正で、取引から quote へ戻る。
+    let as_of = futures::AS_OF.parse()?;
+    let reference = futures::refdata("0D".parse()?, "16:00:00")?;
+    let config = futures::config("")?;
+    let fixings = FixingStore::default();
+    let ctx = Context {
+        as_of,
+        refdata: &reference,
+        config: &config,
+        fixings: &fixings,
+    };
+    let quotes = QuoteSet::new([
+        ("USD.FUT.2026-01".parse()?, 0.96),
+        ("USD.LGM.MR".parse()?, 0.03),
+        ("USD.LGM.VOL".parse()?, 0.01),
+    ])?;
+    let future = StandardProduct::RateFuture(RateFuture::new(
+        futures::CLEARING.parse()?,
+        NonZeroU32::MIN,
+        RateFuturesPrice::new(0.955, DayCount::Act365F)?,
+        Position::Long,
+        None,
+    ));
+    let trade_date = futures::START.parse()?;
+    let trades = |market: &CalibratedMarket| {
+        let prepared = standard(&future, trade_date, market, &ctx)?;
+        Ok(vec![("rate future".to_owned(), prepared)])
+    };
+    check_all(
+        &futures::futures_definition()?,
+        &quotes,
+        &ctx,
+        Currency::USD,
+        1e-8,
+        &trades,
+    )
+}
+
+/// 期間 2025-06-10〜2026-03-10（年末の跳びをまたぐ）の固定 1% 払い・TONA 複利受けの OIS で、
+/// 観測の規約 `method` と、期間の終わりから 2 営業日後の支払い（支払いの遅れ）を持つ。
+fn rfr_swap(case: &ois::Case, method: ObservationMethod) -> Result<StandardProduct> {
+    let calendar = CalendarSpec::Base(case.currency.to_string().parse()?);
+    let rule = ObservationRule::new(
+        ObservationDays::BusinessDays(calendar),
+        method,
+        Aggregation::Compound {
+            form: CompoundingForm::SimpleFactor,
+            spread: SpreadTreatment::Exclusive,
+        },
+        ObservationRounding::default(),
+    )?;
+    let schedule = ScheduleSpec::Explicit(ExplicitPeriods::new(
+        vec![ExplicitPeriod {
+            accrual: DateInterval::new("2025-06-10".parse()?, "2026-03-10".parse()?)?,
+            payment: "2026-03-12".parse()?,
+            reset: None,
+            settlement: None,
+            references: Vec::new(),
+        }],
+        ExplicitDatePolicy::KeepAsConfirmed,
+        PeriodCombination::Sum,
+        None,
+    )?);
+    let notional = Money::new(case.currency, 1_000_000.0);
+    let leg = |side, coupon| {
+        Leg::new(
+            side,
+            notional,
+            schedule.clone(),
+            case.basis(),
+            coupon,
+            NotionalExchange::None,
+            PaymentRounding::None,
+        )
+    };
+    Ok(StandardProduct::Swap(Swap::from_legs([
+        leg(
+            PayReceive::Pay,
+            Coupon::Fixed {
+                rate: Rate::new(0.01, case.basis(), Simple)?,
+            },
+        )?,
+        leg(
+            PayReceive::Receive,
+            Coupon::Compounded {
+                index: case.index.parse()?,
+                rule,
+                spread: Spread::new(0.0)?,
+            },
+        )?,
+    ])?))
+}
+
+#[test]
+fn meeting_steps_and_turn_quote_deltas_match_bumps() -> Result<()> {
+    // R3: 会合日の階段と年末の跳びの TONA。OIS の観測の規約 4 種（Plain は慣行の OIS）と支払いの遅れ、
+    // 年末をまたぐ預金のデルタが、年末のスプレッドを含む全 quote の再校正と一致する。
+    // lookback と lockout は日次積（DailyForwardProduct）で、observation shift は端点の比で評価する。
+    let tona = PenalizedTona::new()?;
+    let inputs = Inputs::new(&[&tona.case], tona.case.as_of, "")?;
+    let ctx = inputs.context();
+    let node = |tenor: &str| -> Result<Node> {
+        Ok(Node {
+            tenor: tenor.parse()?,
+            rate: 0.0,
+            date: tona.case.as_of,
+            df: 1.0,
+        })
+    };
+    let lag = |days: u8| NonZeroU8::new(days).ok_or_else(|| Error::Internal("zero lag".into()));
+    let swaps = [
+        ("lookback", ObservationMethod::Lookback(lag(2)?)),
+        (
+            "observation shift",
+            ObservationMethod::ObservationShift(lag(2)?),
+        ),
+        ("lockout", ObservationMethod::Lockout(lag(1)?)),
+    ]
+    .into_iter()
+    .map(|(name, method)| Ok((name.to_owned(), rfr_swap(&tona.case, method)?)))
+    .collect::<Result<Vec<_>>>()?;
+    let deposit = StandardProduct::Deposit(Deposit::new(
+        Money::new(Currency::JPY, 1_000_000.0),
+        0.004,
+        DepositDates::Confirmed {
+            start: "2025-07-09".parse()?,
+            maturity: "2026-04-09".parse()?,
+        },
+        None,
+        Some(tona.case.basis()),
+        Position::Long,
+        PaymentRounding::None,
+    )?);
+    let trade_date = tona.case.as_of;
+    let trades = |market: &CalibratedMarket| {
+        let mut prepared = [("9M", 0.003), ("18M", 0.0065)]
+            .into_iter()
+            .map(|(tenor, rate)| {
+                Ok((
+                    format!("OIS {tenor}"),
+                    prepare_ois(
+                        &tona.case,
+                        &node(tenor)?,
+                        rate,
+                        &COLLATERAL,
+                        market.market(),
+                        &ctx,
+                    )?,
+                ))
+            })
+            .collect::<Result<Vec<_>>>()?;
+        for (name, product) in &swaps {
+            prepared.push((name.clone(), standard(product, trade_date, market, &ctx)?));
+        }
+        prepared.push((
+            "deposit".to_owned(),
+            standard(&deposit, trade_date, market, &ctx)?,
+        ));
+        Ok(prepared)
+    };
+    check_all(
+        &tona.square_definition()?,
+        &tona.quotes()?,
+        &ctx,
+        Currency::JPY,
+        1e-6,
+        &trades,
+    )
+}
diff --git a/crates/rhodium-risk/tests/quote_delta.rs b/crates/rhodium-risk/tests/quote_delta.rs
index d1b4a20dab281288c8ba46605f655ab730801718..0aaa049c09077fb5fde9ae0cb07889335d5e675a 100644
--- a/crates/rhodium-risk/tests/quote_delta.rs
+++ b/crates/rhodium-risk/tests/quote_delta.rs
@@ -19,7 +19,6 @@ use rhodium_core::{
         EndRule, ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination,
         ScheduleSpec,
     },
-    time::LocalEventTime,
 };
 use rhodium_market::{
     curve::{CurveId, Extrapolation},
@@ -30,9 +29,7 @@ use rhodium_market::{
 use rhodium_pricing::{
     calibration::{CalibratedMarket, GroupOutcome, calibrate, calibrate_with_residual_tolerance},
     context::Context,
-    definition::{
-        DirectFx, FxDefinition, FxImpliedCurve, FxQuotes, FxSettlement, FxSpotLag, MarketDefinition,
-    },
+    definition::{DirectFx, FxDefinition, FxImpliedCurve, FxQuotes, MarketDefinition},
     prepare::{PrepareInput, Prepared, evaluate, prepare},
     resolve::CollateralKey,
 };
@@ -55,6 +52,10 @@ use scioxide_core::{ClockBinding, Deadline};
 #[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
 pub mod ois;
 
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、このテストで使わない関数がある"
+)]
 #[path = "../../rhodium-pricing/tests/fixtures/multicurve.rs"]
 pub mod multicurve;
 
@@ -854,14 +855,6 @@ fn usd_jpy(tona: &Case, sofr: &Case, implied: Option<FxImpliedCurve>) -> Result<
     Ok(FxDefinition {
         pair: "USDJPY".parse()?,
         direct: Some(DirectFx {
-            settlement: FxSettlement {
-                spot_lag: FxSpotLag::TwoBusinessDays,
-                settlement: LocalEventTime {
-                    time: "16:00:00".parse()?,
-                    zone: "Etc/UTC".parse()?,
-                },
-                pip_exponent: 2,
-            },
             market_collateral: collateral.clone(),
             quotes: FxQuotes {
                 spot: "USDJPY.SPOT".parse()?,
@@ -910,6 +903,7 @@ fn translation_fx() -> Result<()> {
         Some(FxImpliedCurve {
             id: "USD-under-JPY".parse()?,
             extrapolation: Extrapolation::default(),
+            long_end: vec![],
         }),
     ] {
         let short = implied.is_some();
diff --git a/crates/rhodium-risk/tests/ui/quote_delta_with_param_gradient.json b/crates/rhodium-risk/tests/ui/quote_delta_with_param_gradient.json
new file mode 100644
index 0000000000000000000000000000000000000000..465de77f078ce87211b8d816bf9505a7bee64b0f
--- /dev/null
+++ b/crates/rhodium-risk/tests/ui/quote_delta_with_param_gradient.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0308",
+    "message": "mismatched types",
+    "spans": [
+      {
+        "file": "tests/ui/quote_delta_with_param_gradient.rs",
+        "line": 6,
+        "column": 31
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-risk/tests/ui/quote_delta_with_param_gradient.stderr b/crates/rhodium-risk/tests/ui/quote_delta_with_param_gradient.stderr
deleted file mode 100644
index 4009b9436742553f2ed0d5d32dcc38d1d5a9194c..0000000000000000000000000000000000000000
--- a/crates/rhodium-risk/tests/ui/quote_delta_with_param_gradient.stderr
+++ /dev/null
@@ -1,15 +0,0 @@
-error[E0308]: mismatched types
- --> tests/ui/quote_delta_with_param_gradient.rs:6:31
-  |
-6 |     let _ = delta.checked_add(gradient);
-  |                   ----------- ^^^^^^^^ expected `&RiskVector<QuoteId>`, found `&RiskVector<ParamId>`
-  |                   |
-  |                   arguments to this method are incorrect
-  |
-  = note: expected reference `&RiskVector<QuoteId>`
-             found reference `&RiskVector<ParamId>`
-note: method defined here
- --> src/vector.rs
-  |
-  |     pub fn checked_add(&self, other: &Self) -> Result<Self> {
-  |            ^^^^^^^^^^^
diff --git a/crates/rhodium-risk/tests/ui/tape_escape.json b/crates/rhodium-risk/tests/ui/tape_escape.json
new file mode 100644
index 0000000000000000000000000000000000000000..4c5c4488d24c5465760eb9c1e2c85b67bd141b9f
--- /dev/null
+++ b/crates/rhodium-risk/tests/ui/tape_escape.json
@@ -0,0 +1,13 @@
+[
+  {
+    "code": "E0521",
+    "message": "borrowed data escapes outside of closure",
+    "spans": [
+      {
+        "file": "tests/ui/tape_escape.rs",
+        "line": 12,
+        "column": 13
+      }
+    ]
+  }
+]
diff --git a/crates/rhodium-risk/tests/ui/tape_escape.stderr b/crates/rhodium-risk/tests/ui/tape_escape.stderr
deleted file mode 100644
index bad1b8485cee4f989c0ec0e2ef00bda21fe03678..0000000000000000000000000000000000000000
--- a/crates/rhodium-risk/tests/ui/tape_escape.stderr
+++ /dev/null
@@ -1,15 +0,0 @@
-error[E0521]: borrowed data escapes outside of closure
-  --> tests/ui/tape_escape.rs:12:13
-   |
- 8 |     let escaped = Cell::new(None::<Var<'_, f64>>);
-   |         ------- `escaped` declared here, outside of the closure body
- 9 |     let _ = tape.linearize(
-10 |         |inputs, outputs| {
-   |          ------ `inputs` is only valid in the closure body
-11 |             // 記録中の変数を外へ保存すると、tape を使い回した後にも読めてしまう。
-12 |             escaped.set(inputs.first().copied());
-   |             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ `inputs` escapes the closure body here
-   |
-   = note: requirement occurs because of the type `Cell<Option<Var<'_, f64>>>`, which makes the generic argument `Option<Var<'_, f64>>` invariant
-   = note: the struct `Cell<T>` is invariant over the parameter `T`
-   = help: see <https://doc.rust-lang.org/nomicon/subtyping.html> for more information about variance
diff --git a/crates/rhodium-risk/tests/xccy_quote_delta.rs b/crates/rhodium-risk/tests/xccy_quote_delta.rs
new file mode 100644
index 0000000000000000000000000000000000000000..c8064e62f1740427645de7142a004f505b1c37d6
--- /dev/null
+++ b/crates/rhodium-risk/tests/xccy_quote_delta.rs
@@ -0,0 +1,204 @@
+//! R6 の網羅表の行: 外貨担保曲線（FX ポイントと `MtM` XCCY）を通る `MtM` XCCY、FX フォワード、
+//! FX スワップの quote デルタ（Δq と Δx）を、quote ごとに全体を再校正した中心差分と照合する。
+
+use rhodium_core::{Error, Result, currency::Currency, date::Date, fx::FxRate, money::Money};
+use rhodium_market::quote::{QuoteId, QuoteKind, QuoteSet};
+use rhodium_pricing::{
+    calibration::{CalibratedMarket, GroupOutcome, calibrate_with_residual_tolerance},
+    context::Context,
+    definition::MarketDefinition,
+    prepare::{PrepareInput, Prepared, evaluate, prepare},
+};
+use rhodium_products::standard::{FxForward, FxSwap, Position, StandardProduct};
+use rhodium_risk::{delta::quote_delta_in, vector::RiskVector};
+use rhodium_testkit::check_aad_vs_bump;
+use scioxide_core::ExecutionContext;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、XCCY のデルタに使わない関数がある"
+)]
+#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
+pub mod ois;
+
+#[allow(
+    dead_code,
+    reason = "共有の fixture のうち、XCCY のデルタに使わない関数がある"
+)]
+#[path = "../../rhodium-pricing/tests/fixtures/xccy.rs"]
+pub mod xccy;
+
+use xccy::{SPOT, XccyMarket, jpy_cash};
+
+const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
+const BP: f64 = 1e-4;
+/// バンプの再校正の残差の許容誤差。品質の照合手順は通常の 1e-14 の 1/100 を求めるが、XCCY のパーの
+/// 残差は倍精度の分解能で約 1.5e-16 までしか下がらない。届く最小の値で解き、`δ_V` に伝える。
+const BUMP_TOLERANCE: f64 = 1e-15;
+
+fn calibrated(
+    definition: &MarketDefinition,
+    quotes: &QuoteSet,
+    ctx: &Context<'_>,
+) -> Result<CalibratedMarket> {
+    let market = calibrate_with_residual_tolerance(definition, quotes, ctx, BUMP_TOLERANCE)?;
+    for group in market.groups() {
+        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
+            return Err(Error::Calibration(format!("{group:?}")));
+        }
+    }
+    Ok(market)
+}
+
+fn bumped(quotes: &QuoteSet, id: &QuoteId, bump: f64) -> Result<QuoteSet> {
+    QuoteSet::new(
+        quotes
+            .iter()
+            .map(|(key, value)| (key.clone(), if key == id { value + bump } else { value })),
+    )
+}
+
+fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
+    vector
+        .entries()
+        .binary_search_by(|(key, _)| key.cmp(id))
+        .ok()
+        .map_or(0.0, |index| vector.entries()[index].1)
+}
+
+/// 照合する取引と、その PV の絶対の下限（JPY）。
+fn products() -> Result<Vec<(&'static str, StandardProduct, f64)>> {
+    let pair = XccyMarket::pair()?;
+    let notional = Money::new(Currency::USD, 1_000_000.0);
+    // quote から 10bp 外した 18M の `MtM` XCCY。テンプレートの元本は JPY 1 円。
+    let xccy = XccyMarket::template("18M")?.build_for_trade(
+        -0.0030,
+        trade_date()?,
+        &XccyMarket::new()?.inputs("")?.context().refdata.clone(),
+    )?;
+    let forward = StandardProduct::FxForward(FxForward::new(
+        notional,
+        FxRate::forward(pair, 147.0, XccyMarket::settlement("2025-10-09")?)?,
+        Position::Long,
+    )?);
+    let swap = StandardProduct::FxSwap(FxSwap::new(
+        notional,
+        FxRate::forward(pair, SPOT, XccyMarket::settlement("2025-01-09")?)?,
+        FxRate::forward(pair, 145.0, XccyMarket::settlement("2026-01-09")?)?,
+        Position::Long,
+    )?);
+    Ok(vec![
+        ("MtM XCCY", xccy, 1e-14),
+        ("FX forward", forward, 1e-6),
+        ("FX swap", swap, 1e-6),
+    ])
+}
+
+fn trade_date() -> Result<Date> {
+    "2025-01-07".parse()
+}
+
+fn prepare_trade(
+    product: &StandardProduct,
+    market: &CalibratedMarket,
+    ctx: &Context<'_>,
+) -> Result<Prepared> {
+    prepare(
+        PrepareInput::Standard {
+            product,
+            trade_date: trade_date()?,
+        },
+        &jpy_cash()?,
+        &[],
+        market.market().structure(),
+        ctx,
+    )
+}
+
+/// 評価時点の X(0) で USD の受渡を JPY に換算した PV。
+fn jpy_pv(prepared: &Prepared, market: &CalibratedMarket) -> Result<f64> {
+    let pv = evaluate(prepared, market.market())?;
+    let x0 = market
+        .market()
+        .structure()
+        .fx_today_plan(XccyMarket::pair()?)?
+        .rate(market.market())?
+        .rate();
+    Ok(pv.get(Currency::JPY) + pv.get(Currency::USD) * x0)
+}
+
+/// quote の種類の感応度の単位（金利とスプレッドは 1bp、FX スポットは相対 1%、FX ポイントは 1 pip）。
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "#203 で各 variant を網羅する"
+)]
+fn unit(market: &CalibratedMarket, id: &QuoteId) -> Result<f64> {
+    Ok(
+        match market
+            .quote_kind(id)
+            .ok_or_else(|| Error::Internal(format!("unknown kind of {id}")))?
+        {
+            QuoteKind::FxSpot => 0.01 * SPOT,
+            QuoteKind::FxForwardPoints => 1.0,
+            _ => BP,
+        },
+    )
+}
+
+#[test]
+fn xccy_fx_forward_and_fx_swap_quote_deltas_match_bumps() -> Result<()> {
+    // JPY 現金担保の下で、USD の受渡は外貨担保曲線、JPY の受渡は TONA で割り引き、報告通貨 JPY に
+    // 評価時点の X(0) で換算する。Δq は TONA、SOFR の OIS と XCCY のベーシスへ、Δx は USDJPY の
+    // スポット（相対 1%）と ON、TN、日付付きのポイント（1 pip）へ、校正の随伴で戻る。
+    let fixture = XccyMarket::new()?;
+    let inputs = fixture.inputs("")?;
+    let ctx = inputs.context();
+    let definition = fixture.definition()?;
+    let quotes = fixture.quotes()?;
+    let base = calibrated(&definition, &quotes, &ctx)?;
+    for (name, product, floor) in products()? {
+        let prepared = prepare_trade(&product, &base, &ctx)?;
+        let report = quote_delta_in(
+            std::slice::from_ref(&prepared),
+            &base,
+            Currency::JPY,
+            SEQUENTIAL,
+        )?;
+        assert!(report.failures().is_empty());
+        let [vector] = report.by_ccy() else {
+            return Err(Error::Internal("expected one reporting currency".into()));
+        };
+        // 再校正の残差の許容誤差を、quote ごとの感応度（自然な単位）の和で PV の誤差に伝える。
+        let scale = vector
+            .entries()
+            .iter()
+            .map(|(id, delta)| Ok(delta.abs() / unit(&base, id)?))
+            .sum::<Result<f64>>()?;
+        for (id, _) in quotes.iter() {
+            let unit = unit(&base, id)?;
+            check_aad_vs_bump(
+                |h| {
+                    let market = calibrated(&definition, &bumped(&quotes, id, h)?, &ctx)?;
+                    jpy_pv(&prepare_trade(&product, &market, &ctx)?, &market)
+                },
+                sensitivity(vector, id),
+                unit,
+                BUMP_TOLERANCE * scale,
+                floor,
+                &format!("{name} {id}"),
+            )?;
+        }
+        // 照合が空振りでないことを確かめる。`MtM` XCCY はベーシスに、FX の取引はスポットに感応する
+        // （`MtM` XCCY は元本のリセットで FX の露出がほぼ打ち消されるので、スポットは確かめない）。
+        let expected = if name == "MtM XCCY" {
+            "USDJPY.XCCY.18M"
+        } else {
+            "USDJPY.SPOT"
+        };
+        assert!(
+            sensitivity(vector, &expected.parse()?).abs() > floor,
+            "{name} {expected}"
+        );
+    }
+    Ok(())
+}
diff --git a/crates/rhodium-testkit/Cargo.toml b/crates/rhodium-testkit/Cargo.toml
index cd3bd821e20de8aa7ee39f5bc43fbb105f43240e..e88fe27a26218d8f4e0e4dd3ac75a4eb30f79351 100644
--- a/crates/rhodium-testkit/Cargo.toml
+++ b/crates/rhodium-testkit/Cargo.toml
@@ -12,4 +12,6 @@ workspace = true
 
 [dependencies]
 scioxide-core.workspace = true
-trybuild.workspace = true
+serde.workspace = true
+serde_json = "1"
+tempfile = "3.27"
diff --git a/crates/rhodium-testkit/src/lib.rs b/crates/rhodium-testkit/src/lib.rs
index 23d236a0cc4ebaf5462c8dcf62cbf582c85af064..4a5989eb1bdbd13d8fe47e91d1267fc2a112cbca 100644
--- a/crates/rhodium-testkit/src/lib.rs
+++ b/crates/rhodium-testkit/src/lib.rs
@@ -1,5 +1,10 @@
 //! テストの補助と独立なオラクルを持つ開発専用の crate。
 
+mod ui;
+
+#[doc(hidden)]
+pub use ui::run_compile_fail_pair;
+
 /// 呼び出し元の `tests/ui/` にある `<case>_ok.rs` を実行し、`<case>.rs` の診断を照合する。
 ///
 /// ケース名は ASCII 小文字・数字・アンダースコアからなる空でない文字列リテラルで指定する。
@@ -24,18 +29,14 @@ macro_rules! compile_fail_pair {
             }
             $case
         };
-        $crate::run_compile_fail_pair(CASE);
+        #[expect(clippy::disallowed_macros, reason = "型テストの fixture は呼び出し元の crate に属する")]
+        const ROOT: &str = env!("CARGO_MANIFEST_DIR");
+        if let Err(error) = $crate::run_compile_fail_pair(ROOT, CASE) {
+            panic!("UI case {CASE}: {error}");
+        }
     }};
 }
 
-#[doc(hidden)]
-pub fn run_compile_fail_pair(case: &str) {
-    // trybuild はテストを実行する crate のマニフェストを起点にパスを解決する。
-    let cases = trybuild::TestCases::new();
-    cases.pass(format!("tests/ui/{case}_ok.rs"));
-    cases.compile_fail(format!("tests/ui/{case}.rs"));
-}
-
 /// 数値を相対誤差と絶対誤差の大きい方で照合する。
 ///
 /// `|x − y| ≤ max(rel · M, abs)` を確かめる。M は `max(|x|, |y|)` で、5 番目の引数を渡すと
diff --git a/crates/rhodium-testkit/src/ui.rs b/crates/rhodium-testkit/src/ui.rs
new file mode 100644
index 0000000000000000000000000000000000000000..e195f3ad20109aae2d07968d951935f1f2ee2774
--- /dev/null
+++ b/crates/rhodium-testkit/src/ui.rs
@@ -0,0 +1,345 @@
+//! 型テストのコンパイル、正常例の実行、失敗例の診断照合。
+
+use std::{
+    collections::{BTreeMap, BTreeSet},
+    env::consts::{DLL_EXTENSION, EXE_SUFFIX},
+    error::Error,
+    fs,
+    path::{Path, PathBuf},
+    process::{Command, Output},
+    sync::OnceLock,
+};
+
+use serde::{Deserialize, Serialize};
+use serde_json::Value;
+
+type Result<T> = std::result::Result<T, Box<dyn Error>>;
+
+fn checked(command: &mut Command) -> Result<Output> {
+    let output = command.output()?;
+    if !output.status.success() {
+        return Err(format!(
+            "{command:?} exited with {}\n{}\n{}",
+            output.status,
+            String::from_utf8_lossy(&output.stdout),
+            String::from_utf8_lossy(&output.stderr)
+        )
+        .into());
+    }
+    Ok(output)
+}
+
+#[derive(Deserialize)]
+struct Node {
+    id: String,
+    deps: Vec<Dependency>,
+}
+
+#[derive(Deserialize)]
+struct Dependency {
+    name: String,
+    pkg: String,
+    dep_kinds: Vec<DependencyKind>,
+}
+
+#[derive(Deserialize)]
+struct DependencyKind {
+    kind: Option<String>,
+}
+
+#[derive(Deserialize)]
+struct Artifact {
+    package_id: String,
+    manifest_path: PathBuf,
+    filenames: Vec<PathBuf>,
+    target: Target,
+}
+
+#[derive(Deserialize)]
+struct Target {
+    name: String,
+    kind: Vec<String>,
+    edition: String,
+}
+
+struct Compiler {
+    args: Vec<String>,
+}
+
+// Cargo がこのテスト用に解決した package ID と成果物を結ぶ。同じ target に残る別の版は読まない。
+impl Compiler {
+    #[expect(
+        clippy::disallowed_macros,
+        reason = "開発用の型テストを起動した Cargo を使う"
+    )]
+    fn prepare(root: &Path) -> Result<Self> {
+        let cargo = env!("CARGO");
+        let build = checked(Command::new(cargo).current_dir(root).args([
+            "build",
+            "--locked",
+            "--profile",
+            "test",
+            "--test",
+            "compile_fail",
+            "--all-features",
+            "--message-format=json-render-diagnostics",
+        ]))?;
+        let mut artifacts = Vec::new();
+        let mut search = BTreeSet::new();
+        for line in build
+            .stdout
+            .split(|byte| *byte == b'\n')
+            .filter(|row| !row.is_empty())
+        {
+            let msg: Value = serde_json::from_slice(line)?;
+            if msg["reason"] == "compiler-artifact" {
+                artifacts.push(serde_json::from_value::<Artifact>(msg)?);
+            } else if msg["reason"] == "build-script-executed" {
+                for path in serde_json::from_value::<Vec<String>>(msg["linked_paths"].clone())? {
+                    search.insert(path);
+                }
+            }
+        }
+        let owner = artifacts
+            .iter()
+            .find(|a| {
+                a.manifest_path == root.join("Cargo.toml")
+                    && a.target.kind.iter().any(|k| k == "lib")
+            })
+            .ok_or("Cargo did not build the library under test")?;
+        // 実行する正常例と同じホストの依存だけを解決する。
+        let host = checked(
+            Command::new("rustc")
+                .current_dir(root)
+                .args(["--print", "host-tuple"]),
+        )?;
+        let host = std::str::from_utf8(&host.stdout)?.trim();
+        let metadata = checked(Command::new(cargo).current_dir(root).args([
+            "metadata",
+            "--filter-platform",
+            host,
+            "--locked",
+            "--all-features",
+            "--format-version",
+            "1",
+        ]))?;
+        let mut meta: Value = serde_json::from_slice(&metadata.stdout)?;
+        let nodes: Vec<Node> = serde_json::from_value(meta["resolve"]["nodes"].take())?;
+        let node = nodes
+            .into_iter()
+            .find(|node| node.id == owner.package_id)
+            .ok_or("package missing from Cargo dependency graph")?;
+        let mut direct: BTreeMap<_, _> = node
+            .deps
+            .into_iter()
+            .filter(|dep| {
+                dep.dep_kinds
+                    .iter()
+                    .any(|k| k.kind.as_deref() != Some("build"))
+            })
+            .map(|dep| (dep.pkg, dep.name))
+            .collect();
+        direct.insert(owner.package_id.clone(), owner.target.name.clone());
+        let mut args = vec![format!("--edition={}", owner.target.edition)];
+        for artifact in &artifacts {
+            for file in &artifact.filenames {
+                let extension = file.extension().and_then(|ext| ext.to_str());
+                if extension == Some("rlib") || extension == Some(DLL_EXTENSION) {
+                    let dir = file.parent().ok_or("artifact missing parent directory")?;
+                    search.insert(format!("dependency={}", dir.display()));
+                    if let Some(name) = direct.get(&artifact.package_id) {
+                        args.extend(["--extern".into(), format!("{name}={}", file.display())]);
+                    }
+                }
+            }
+        }
+        for path in search {
+            args.extend(["-L".into(), path]);
+        }
+        Ok(Self { args })
+    }
+
+    fn command(&self, root: &Path, source: &Path) -> Command {
+        let mut cmd = Command::new("rustc");
+        cmd.current_dir(root).args(&self.args).arg(source);
+        cmd
+    }
+}
+
+#[derive(Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord)]
+#[serde(deny_unknown_fields)]
+struct Expected {
+    code: Option<String>,
+    message: String,
+    spans: Vec<Location>,
+}
+
+#[derive(Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord)]
+#[serde(deny_unknown_fields)]
+struct Location {
+    file: String,
+    line: usize,
+    column: usize,
+}
+
+#[derive(Deserialize)]
+struct Diagnostic {
+    level: String,
+    message: String,
+    code: Option<Code>,
+    spans: Vec<Span>,
+}
+
+#[derive(Deserialize)]
+struct Code {
+    code: String,
+}
+
+#[derive(Deserialize)]
+struct Span {
+    is_primary: bool,
+    file_name: String,
+    line_start: usize,
+    column_start: usize,
+}
+
+fn diagnostics(output: &Output) -> Result<Vec<Expected>> {
+    let mut errors = Vec::new();
+    for row in output
+        .stderr
+        .split(|byte| *byte == b'\n')
+        .filter(|row| !row.is_empty())
+    {
+        let d: Diagnostic = serde_json::from_slice(row)?;
+        if d.level != "error" || (d.spans.is_empty() && d.message.starts_with("aborting due to ")) {
+            continue;
+        }
+        errors.push(Expected {
+            code: d.code.map(|c| c.code),
+            message: d.message,
+            spans: d
+                .spans
+                .into_iter()
+                .filter(|s| s.is_primary)
+                .map(|s| Location {
+                    file: s.file_name.replace('\\', "/"),
+                    line: s.line_start,
+                    column: s.column_start,
+                })
+                .collect(),
+        });
+    }
+    errors.sort();
+    Ok(errors)
+}
+
+/// 呼び出し元の統合テストごとに Cargo の準備を共有し、ケースごとに独立した出力先を使う。
+#[doc(hidden)]
+pub fn run_compile_fail_pair(manifest_dir: &str, case: &str) -> Result<()> {
+    static COMPILER: OnceLock<std::result::Result<Compiler, String>> = OnceLock::new();
+    let root = Path::new(manifest_dir);
+    let compiler = COMPILER
+        .get_or_init(|| Compiler::prepare(root).map_err(|e| e.to_string()))
+        .as_ref()
+        .map_err(Clone::clone)?;
+    check_pair(compiler, root, case)
+}
+
+fn check_pair(compiler: &Compiler, root: &Path, case: &str) -> Result<()> {
+    let tmp = tempfile::tempdir()?;
+    let executable = tmp.path().join(format!("positive{EXE_SUFFIX}"));
+    let ui = Path::new("tests/ui");
+    checked(
+        compiler
+            .command(root, &ui.join(format!("{case}_ok.rs")))
+            .arg("-o")
+            .arg(&executable),
+    )?;
+    checked(Command::new(&executable).current_dir(root))?;
+
+    let source = ui.join(format!("{case}.rs"));
+    let output = compiler
+        .command(root, &source)
+        .args([
+            "--error-format=json",
+            "--diagnostic-width=300",
+            "--emit=metadata",
+            "--out-dir",
+        ])
+        .arg(tmp.path())
+        .output()?;
+    // rustc はコンパイルエラーで 1、ICE で 101 を返す。1 以外では診断を照合できないので、
+    // 原因を追えるよう rustc の出力をそのまま添える。
+    if output.status.code() != Some(1) {
+        return Err(format!(
+            "expected a rustc error for {}, got {}\n{}",
+            source.display(),
+            output.status,
+            String::from_utf8_lossy(&output.stderr)
+        )
+        .into());
+    }
+    let actual = diagnostics(&output)?;
+    let mut expected: Vec<Expected> =
+        serde_json::from_slice(&fs::read(root.join(source.with_extension("json")))?)?;
+    expected.sort();
+    if actual.is_empty() || actual != expected {
+        return Err(format!(
+            "diagnostics differ for {}\nexpected: {}\nactual: {}",
+            root.join(source).display(),
+            serde_json::to_string(&expected)?,
+            serde_json::to_string(&actual)?
+        )
+        .into());
+    }
+    Ok(())
+}
+
+#[cfg(test)]
+mod tests {
+    use super::*;
+
+    // 同じ失敗例で診断の個数・位置、失敗例のコンパイルの成否、正常例の終了コードを変え、
+    // それぞれが失敗になることを確認する。
+    #[test]
+    fn pair_requires_a_successful_run_and_all_expected_error_locations() -> Result<()> {
+        let tmp = tempfile::tempdir()?;
+        let root = tmp.path();
+        let ui = root.join("tests/ui");
+        fs::create_dir_all(&ui)?;
+        let ok = ui.join("sample_ok.rs");
+        let fail = ui.join("sample.rs");
+        fs::write(&ok, "fn main() {}")?;
+        fs::write(&fail, "compile_error!(\"first\");\nfn main() {}")?;
+        let expected = ui.join("sample.json");
+        let diagnostic = serde_json::json!({
+            "code": null, "message": "first",
+            "spans": [{"file": "tests/ui/sample.rs", "line": 1, "column": 1}]
+        });
+        fs::write(&expected, serde_json::to_vec(&[&diagnostic])?)?;
+        let compiler = Compiler { args: vec![] };
+        check_pair(&compiler, root, "sample")?;
+
+        fs::write(&expected, serde_json::to_vec(&[&diagnostic, &diagnostic])?)?;
+        let mismatch = || {
+            check_pair(&compiler, root, "sample")
+                .unwrap_err()
+                .to_string()
+        };
+        assert!(mismatch().contains("diagnostics differ"));
+        fs::write(&expected, serde_json::to_vec(&[&diagnostic])?)?;
+        fs::write(&fail, "\ncompile_error!(\"first\");\nfn main() {}")?;
+        assert!(mismatch().contains("diagnostics differ"));
+        fs::write(
+            &fail,
+            "compile_error!(\"first\");\ncompile_error!(\"extra\");\nfn main() {}",
+        )?;
+        assert!(mismatch().contains("diagnostics differ"));
+        // 失敗例がコンパイルに通ったら、診断を比べる前に失敗にする。
+        fs::write(&fail, "fn main() {}")?;
+        assert!(mismatch().contains("expected a rustc error"));
+        fs::write(&ok, "fn main() { std::process::exit(17); }")?;
+        assert!(mismatch().contains("17"));
+        Ok(())
+    }
+}
diff --git a/crates/rhodium-testkit/tests/compile_fail_pairs.rs b/crates/rhodium-testkit/tests/compile_fail_pairs.rs
index 17c90575aef7b364948c720aaee68b4f297504ca..51be07ef7aad9207e7367e1cca788b783690a73b 100644
--- a/crates/rhodium-testkit/tests/compile_fail_pairs.rs
+++ b/crates/rhodium-testkit/tests/compile_fail_pairs.rs
@@ -6,9 +6,11 @@ use std::{
     path::{Path, PathBuf},
 };
 
-// 各ケースが成功例と stderr を登録していること、登録されていないケースがないことを確認する。
-// また、compile-fail の例を対の登録の外に置く書き方（doctest の compile_fail と trybuild の直接呼び出し）を拒否する。
+// 各ケースが成功例と診断 JSON を登録していること、登録されていないケースがないことを確認する。
 // 各 crate のケースの実行と診断の照合は compile_fail_pair! に任せる。
+// compile-fail の例を対の登録の外に置く書き方（doctest の compile_fail、ランナーの直接呼び出し、
+// trybuild の使用）は、.githooks/checks.py の禁止パターン（P-doctest-compile-fail、P-ui-runner、
+// P-trybuild）が拒否する。
 #[test]
 #[expect(
     clippy::disallowed_macros,
@@ -21,9 +23,6 @@ fn each_compile_failure_has_an_executed_positive_case() -> io::Result<()> {
             continue;
         }
         check_ui(&directory)?;
-        for root in ["src", "tests", "examples", "benches", "build.rs"] {
-            check_sources(&directory.join(root))?;
-        }
     }
     Ok(())
 }
@@ -49,7 +48,7 @@ fn check_ui(directory: &Path) -> io::Result<()> {
             expected.extend([
                 ui.join(format!("{name}.rs")),
                 ui.join(format!("{name}_ok.rs")),
-                ui.join(format!("{name}.stderr")),
+                ui.join(format!("{name}.json")),
             ]);
         }
         assert!(
@@ -66,7 +65,7 @@ fn check_ui(directory: &Path) -> io::Result<()> {
                 path.is_file()
                     && matches!(
                         path.extension().and_then(|ext| ext.to_str()),
-                        Some("rs" | "stderr")
+                        Some("rs" | "json" | "stderr")
                     )
             })
             .collect()
@@ -83,52 +82,6 @@ fn check_ui(directory: &Path) -> io::Result<()> {
     Ok(())
 }
 
-// 対の登録を通らない compile-fail の例を拒否する。doctest の compile_fail は成功例と診断を
-// 対にできず、trybuild の直接呼び出しは tests/ui の外のファイルを検査から外せるため。
-fn check_sources(path: &Path) -> io::Result<()> {
-    if path.is_dir() {
-        for child in children(path)? {
-            check_sources(&child)?;
-        }
-    } else if path.is_file()
-        && path.extension().is_some_and(|ext| ext == "rs")
-        // この検査自身は探す文字列をソースに含むので対象から外す。
-        && !path.ends_with("rhodium-testkit/tests/compile_fail_pairs.rs")
-    {
-        // trybuild の直接呼び出しはマクロの実装ファイルだけに許可する。
-        let allows_trybuild = path.ends_with("rhodium-testkit/src/lib.rs");
-        let source = fs::read_to_string(path)?;
-        for (number, line) in (1_u64..).zip(source.lines()) {
-            let line = line.trim();
-            assert!(
-                allows_trybuild
-                    || !(line.contains("trybuild::")
-                        || line.contains("TestCases")
-                        || line.contains("run_compile_fail_pair")),
-                "{}:{number}: register the case with compile_fail_pair! instead of calling trybuild",
-                path.display(),
-            );
-            let body = line
-                .strip_prefix("///")
-                .or_else(|| line.strip_prefix("//!"))
-                .or_else(|| line.strip_prefix('*'))
-                .map(str::trim_start);
-            let fence = body.and_then(|body| {
-                body.strip_prefix("```")
-                    .or_else(|| body.strip_prefix("~~~"))
-            });
-            assert!(
-                !fence.is_some_and(|header| header
-                    .split([',', ' ', '\t'])
-                    .any(|flag| flag == "compile_fail")),
-                "{}:{number}: move the compile-fail example to a paired trybuild case",
-                path.display(),
-            );
-        }
-    }
-    Ok(())
-}
-
 fn children(directory: &Path) -> io::Result<BTreeSet<PathBuf>> {
     fs::read_dir(directory)?
         .map(|entry| entry.map(|entry| entry.path()))
diff --git a/crates/rhodium/Cargo.toml b/crates/rhodium/Cargo.toml
index 84271ca59f465010b6bebe2c88f03c640f027af5..6957fe7c312c497d7de5de601f40cd6ce92286cc 100644
--- a/crates/rhodium/Cargo.toml
+++ b/crates/rhodium/Cargo.toml
@@ -35,3 +35,7 @@ harness = false
 [[bench]]
 name = "active_market"
 harness = false
+
+[[bench]]
+name = "m1"
+harness = false
diff --git a/crates/rhodium/benches/m1.rs b/crates/rhodium/benches/m1.rs
new file mode 100644
index 0000000000000000000000000000000000000000..d6dd0c3e7c0a2563a6412ba5d46420f6cc96e202
--- /dev/null
+++ b/crates/rhodium/benches/m1.rs
@@ -0,0 +1,254 @@
+//! M1 の所要時間: 最初の検証データセットの校正と、M1 の商品を混ぜたデスクの準備、PV、quote デルタ、
+//! quote の更新。
+//!
+//! `cargo bench -p rhodium --bench m1` で実行する。入力は `examples/m1-rates` で、デスクはその
+//! 18 件の取引を ID を変えて繰り返したもの（同じ商品の繰り返しなので、観測の列の共有が最大に効く）。
+//! 結果は日付、コミット、マシンとともに `docs/benchmarks/m1.md` に残す。
+//! 計測の名前は `m1_rates/` で始め、M0 のデスクで M1-7 を測る `m0` の `m1/quote_update_10k` と区別する。
+
+use std::{
+    fs,
+    hint::black_box,
+    num::NonZeroUsize,
+    path::{Path, PathBuf},
+};
+
+use criterion::{BatchSize, Criterion};
+use rhodium::{
+    prelude::{ExecutionContext, ThreadExecutor},
+    session::{GroupOutcome, Session},
+    workspace::Workspace,
+};
+use rhodium_core::time::Instant;
+use rhodium_market::quote::QuoteOverride;
+use tempfile::TempDir;
+
+/// 1 万件のデスク。
+const DESK: usize = 10_000;
+/// 10 万件のデスク（性能の目安の「線形な取引 10 万件の PV と quote デルタ」）。
+const LARGE_DESK: usize = 100_000;
+
+#[expect(
+    clippy::disallowed_macros,
+    reason = "リポジトリのサンプルの場所をベンチのソースの位置から決める"
+)]
+fn sample() -> PathBuf {
+    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m1-rates")
+}
+
+#[expect(
+    clippy::expect_used,
+    reason = "計測の入力を作れなければ、計測を続けられない"
+)]
+fn copy_dir(from: &Path, to: &Path) {
+    for entry in fs::read_dir(from).expect("sample directory") {
+        let entry = entry.expect("sample entry");
+        let target = to.join(entry.file_name());
+        if entry.file_type().expect("entry type").is_dir() {
+            fs::create_dir(&target).expect("directory");
+            copy_dir(&entry.path(), &target);
+        } else {
+            fs::copy(entry.path(), target).expect("file");
+        }
+    }
+}
+
+/// サンプルの取引ファイルの `[[trade]]` ごとの塊（文書の `version` の行を除く）。
+#[expect(
+    clippy::expect_used,
+    reason = "計測の入力を作れなければ、計測を続けられない"
+)]
+fn trade_blocks() -> Vec<String> {
+    let mut paths: Vec<_> = fs::read_dir(sample().join("trades"))
+        .expect("trades directory")
+        .map(|entry| entry.expect("trade file").path())
+        .collect();
+    paths.sort();
+    paths
+        .iter()
+        .flat_map(|path| {
+            let text = fs::read_to_string(path).expect("trade file");
+            text.split("\n[[trade]]\n")
+                .skip(1)
+                .map(|block| format!("\n[[trade]]\n{block}"))
+                .collect::<Vec<_>>()
+        })
+        .collect()
+}
+
+/// サンプルの市場と参照データに、取引を `trades` 件まで繰り返したデスクを置く。
+#[expect(
+    clippy::expect_used,
+    reason = "計測の入力を作れなければ、計測を続けられない"
+)]
+fn desk(trades: usize) -> TempDir {
+    let dir = tempfile::tempdir().expect("temporary directory");
+    copy_dir(&sample(), dir.path());
+    let trades_dir = dir.path().join("trades");
+    fs::remove_dir_all(&trades_dir).expect("remove sample trades");
+    fs::create_dir(&trades_dir).expect("trades directory");
+    let blocks = trade_blocks();
+    let mut text = String::from("version = 1\n");
+    for (index, block) in blocks.iter().cycle().take(trades).enumerate() {
+        // ID の行だけを変え、条件は同じ取引を繰り返す。
+        let renamed = block.replacen("\nid = \"", &format!("\nid = \"D{index:06}-"), 1);
+        text.push_str(&renamed);
+    }
+    fs::write(trades_dir.join("desk.toml"), text).expect("desk trades");
+    dir
+}
+
+#[expect(
+    clippy::expect_used,
+    reason = "計測の入力を読めなければ、計測を続けられない"
+)]
+fn load(dir: &Path) -> (Workspace, Instant) {
+    let workspace = Workspace::load(dir).expect("m1 workspace loads");
+    let as_of = workspace
+        .resolve_date("2026-11-24".parse().expect("valuation date"))
+        .expect("valuation close");
+    (workspace, as_of)
+}
+
+fn threads() -> NonZeroUsize {
+    std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN)
+}
+
+/// 計測の前に、校正の全グループが解け、全取引の PV が成功し、デルタが取引を除かないことを確かめる。
+///
+/// 受渡通貨別の PV とデルタに加えて、計測する報告通貨の集計（合計と、報告通貨のデルタ）も
+/// 確かめる。換算に失敗した取引は合計から黙って除かれるので、除外が無いことを計測の外で見る。
+#[expect(
+    clippy::expect_used,
+    reason = "計測の前提が崩れていれば、計測を続けられない"
+)]
+fn validate(session: &Session, execution: ExecutionContext<'_>, trades: usize) {
+    assert!(
+        session
+            .calibrated_market()
+            .groups()
+            .iter()
+            .all(|group| matches!(group.outcome(), GroupOutcome::Solved(_)))
+    );
+    let pv = session.pv(execution).expect("pv");
+    assert_eq!(pv.len(), trades);
+    assert!(pv.values().all(Result::is_ok));
+    let portfolio = session.value_portfolio(execution).expect("portfolio");
+    assert!(portfolio.excluded.is_empty(), "{:?}", portfolio.excluded);
+    assert!(portfolio.total.is_ok(), "{:?}", portfolio.total);
+    if trades > 0 {
+        let delta = session.quote_delta(execution).expect("quote delta");
+        assert!(delta.failures().is_empty(), "{:?}", delta.failures());
+        let reporting = session
+            .quote_delta_reporting(execution)
+            .expect("reporting quote delta");
+        assert!(
+            reporting.failures().is_empty(),
+            "{:?}",
+            reporting.failures()
+        );
+    }
+}
+
+#[expect(
+    clippy::expect_used,
+    reason = "計測の入力から市場を作れなければ、計測を続けられない"
+)]
+fn calibration(c: &mut Criterion) {
+    // 取引の無いサンプルで Session::new を 1 スレッドで測ると、最初の検証データセットの全曲線の
+    // 校正（9 グループ、TONA の階段と跳び、TIBOR 6M・3M、SOFR、FX ポイントと XCCY）だけになる。
+    let dir = desk(0);
+    let (workspace, as_of) = load(dir.path());
+    let session = Session::new(workspace, as_of, ExecutionContext::sequential()).expect("session");
+    validate(&session, ExecutionContext::sequential(), 0);
+    c.bench_function("m1_rates/calibration", |b| {
+        b.iter_batched(
+            || load(dir.path()),
+            |(workspace, as_of)| {
+                Session::new(workspace, as_of, ExecutionContext::sequential()).expect("session")
+            },
+            BatchSize::PerIteration,
+        );
+    });
+}
+
+#[expect(
+    clippy::expect_used,
+    reason = "計測の入力から市場を作れなければ、計測を続けられない"
+)]
+fn desk_benchmarks(c: &mut Criterion, trades: usize, name: &str) {
+    // 読み込み、準備（校正込み）、PV、quote デルタ、TONA の 10Y の変更からの再計算を、使える全スレッドで測る。
+    let dir = desk(trades);
+    let executor = ThreadExecutor::new(threads());
+    let execution = ExecutionContext::new(&executor);
+    c.bench_function(&format!("m1_rates/workspace_load_{name}"), |b| {
+        b.iter(|| black_box(load(dir.path())));
+    });
+    c.bench_function(&format!("m1_rates/session_new_{name}"), |b| {
+        b.iter_batched(
+            || load(dir.path()),
+            |(workspace, as_of)| Session::new(workspace, as_of, execution).expect("session"),
+            BatchSize::PerIteration,
+        );
+    });
+    let (workspace, as_of) = load(dir.path());
+    let session = Session::new(workspace, as_of, execution).expect("session");
+    validate(&session, execution, trades);
+    c.bench_function(&format!("m1_rates/pv_{name}"), |b| {
+        b.iter(|| black_box(session.value_portfolio(execution).expect("portfolio")));
+    });
+    c.bench_function(&format!("m1_rates/quote_delta_{name}"), |b| {
+        b.iter(|| {
+            black_box(
+                session
+                    .quote_delta_reporting(execution)
+                    .expect("quote delta"),
+            )
+        });
+    });
+    let overrides = session
+        .quote_overrides([QuoteOverride::Set {
+            id: "JPY.TONA.OIS.10Y".parse().expect("quote ID"),
+            value: 0.0133,
+        }])
+        .expect("quote override");
+    validate(
+        &session.with_quotes(&overrides, execution).expect("update"),
+        execution,
+        trades,
+    );
+    c.bench_function(&format!("m1_rates/quote_update_{name}"), |b| {
+        b.iter(|| {
+            let bumped = session
+                .with_quotes(black_box(&overrides), execution)
+                .expect("quote update");
+            black_box(bumped.value_portfolio(execution).expect("portfolio"));
+            black_box(
+                bumped
+                    .quote_delta_reporting(execution)
+                    .expect("quote delta"),
+            );
+        });
+    });
+}
+
+fn desk_10k(c: &mut Criterion) {
+    desk_benchmarks(c, DESK, "10k");
+}
+
+fn desk_100k(c: &mut Criterion) {
+    desk_benchmarks(c, LARGE_DESK, "100k");
+}
+
+// criterion の生成する関数を公開しないよう、非公開のモジュールに置く。
+mod harness {
+    use criterion::Criterion;
+
+    criterion::criterion_group! {
+        name = benches;
+        config = Criterion::default().sample_size(10);
+        targets = super::calibration, super::desk_10k, super::desk_100k
+    }
+}
+
+criterion::criterion_main!(harness::benches);
diff --git a/crates/rhodium/src/io.rs b/crates/rhodium/src/io.rs
index ac93b9c940c19c5c9101f9600e99f62da92ec000..5439fb04f58d263f52bcb68eaece93939cf04fcf 100644
--- a/crates/rhodium/src/io.rs
+++ b/crates/rhodium/src/io.rs
@@ -73,6 +73,10 @@ pub(crate) fn normalize_workspace_datetimes(value: &mut toml::Value) -> Result<(
 }
 
 // ドメイン型の共通表現は文字列。日付と瞬間のフィールドだけでリテラルを変換する。
+#[expect(
+    clippy::wildcard_enum_match_arm,
+    reason = "外部の enum は variant が増えうるので、既定の腕で受ける"
+)]
 fn normalize_datetimes(value: &mut toml::Value, field: &str, workspace: bool) -> Result<()> {
     match value {
         toml::Value::Datetime(datetime) => {
@@ -104,7 +108,16 @@ fn normalize_datetimes(value: &mut toml::Value, field: &str, workspace: bool) ->
                 || workspace
                     && matches!(
                         field,
-                        "from" | "through" | "at" | "effective_at" | "pillars" | "breaks"
+                        "from"
+                            | "through"
+                            | "at"
+                            | "effective_at"
+                            | "pillars"
+                            | "pillar"
+                            | "breaks"
+                            | "grid"
+                            | "meetings"
+                            | "dates"
                     );
             let text = match (datetime.date, datetime.time, datetime.offset) {
                 (Some(date), None, None) if date_field => date.to_string(),
diff --git a/crates/rhodium/src/lib.rs b/crates/rhodium/src/lib.rs
index 3e753e74d7880d9d2bffd2a97681d9ab5230dc58..9fd9f260402bc8d1ce14b32a41cf667c1c16b64e 100644
--- a/crates/rhodium/src/lib.rs
+++ b/crates/rhodium/src/lib.rs
@@ -18,5 +18,6 @@ pub mod prelude {
     pub use rhodium_core::{
         Error, Result,
         parallel::{ExecutionContext, ThreadExecutor},
+        quantity::VolKind,
     };
 }
diff --git a/crates/rhodium/src/session.rs b/crates/rhodium/src/session.rs
index 41396e28cb1f07bd5dd8cd283df9a25108d1e521..e7f4e23f764aeba9f4bf99eef42905c8ac56ab9a 100644
--- a/crates/rhodium/src/session.rs
+++ b/crates/rhodium/src/session.rs
@@ -11,8 +11,7 @@ use rhodium_core::{
 };
 use rhodium_market::quote::{QuoteId, QuoteOverride, QuoteOverrides};
 use rhodium_pricing::{
-    calibration::{CalibratedMarket, calibrate},
-    context::Context,
+    calibration::CalibratedMarket,
     prepare::{PrepareInput, Prepared, evaluate, prepare},
     valuation,
 };
@@ -85,19 +84,20 @@ impl Session {
         as_of: Instant,
         execution: ExecutionContext<'_>,
     ) -> Result<Self> {
-        let context = context(&workspace, as_of);
+        let context = workspace.context(as_of);
         // 校正は実行器を使わないので、取り消した実行で校正を始めないよう先に確かめる。
         parallel::check(execution)?;
-        let calibrated = calibrate(&workspace.inputs.market, workspace.quotes(), &context)?;
+        let calibrated = workspace.calibrate(as_of)?;
         let structure = calibrated.market().structure();
         let trades: Vec<_> = workspace.trades().collect();
         let plans = map_partitions(execution, trades.len(), |range, child| {
             trades[range]
                 .iter()
-                .map(|&(_, product, collateral)| {
+                .map(|&(_, resolved, collateral)| {
                     parallel::check(child)?;
+                    // 下ろした IR を使い回し、方式の表は取引の商品の種類で引く。
                     Ok(prepare(
-                        PrepareInput::Ir(product),
+                        PrepareInput::Trade(resolved),
                         collateral,
                         &[],
                         structure,
@@ -262,7 +262,7 @@ impl Session {
         valuation::value_portfolio(
             &prepared,
             &self.calibrated,
-            &context(&self.workspace, self.as_of()),
+            &self.workspace.context(self.as_of()),
             execution,
         )
     }
@@ -315,12 +315,3 @@ impl Session {
         Ok(QuoteDeltaReport { delta, failures })
     }
 }
-
-fn context(workspace: &Workspace, as_of: Instant) -> Context<'_> {
-    Context {
-        as_of,
-        refdata: &workspace.inputs.refdata,
-        fixings: &workspace.inputs.fixings,
-        config: &workspace.inputs.config,
-    }
-}
