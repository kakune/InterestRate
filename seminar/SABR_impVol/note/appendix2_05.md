    )?;
    assert!(matches!(df(&partial, &slots), Err(Error::MissingData(_))));
    Ok(())
}

#[test]
fn scale_wraps_nested_max_forward_with_distinct_grids() -> Result<()> {
    // 内側の最大フォワードは 4%、6%。外側の最初の半年は 4.5% を採り、
    // 残りの 1 年半は内側の後半だけを平均する。負の倍率は外側の選択の後に掛ける。
    let half = "2026-07-02T12:00:00Z".parse()?;
    let first = FIRST.parse()?;
    let second = SECOND.parse()?;
    let nested = CurveExpr::MaxForward {
        parts: vec![
            CurveExpr::MaxForward {
                parts: vec![
                    CurveExpr::Flat(0.04),
                    CurveExpr::Stitch {
                        short: Box::new(CurveExpr::Flat(0.03)),
                        long: Box::new(CurveExpr::Flat(0.06)),
                        at: first,
                    },
                ],
                grid: vec![first, second],
            },
            CurveExpr::Flat(0.045),
        ],
        grid: vec![half, second],
    };
    let cases = [("NEGATIVE", -2.0), ("ZERO", 0.0)];
    let curves = cases.map(|(id, factor)| {
        CurveSpec::from_expr(
            id.parse()?,
            Currency::JPY,
            [CurveRole::Discount],
            CurveExpr::Scale {
                factor,
                body: Box::new(nested.clone()),
            },
        )
    });
    let [negative, zero] = curves;
    let structure = Arc::new(MarketStructure::new(ORIGIN.parse()?, [negative?, zero?])?);
    let market = Market::from_theta(Arc::clone(&structure), vec![])?;
    let first_half = 0.045 * 0.5;
    let last_three_halves = 0.04 * 0.5 + 0.06;
    for (date, forward_integral) in [
        (FIRST, first_half + last_three_halves / 3.0),
        (SECOND, first_half + last_three_halves),
    ] {
        for (id, factor) in cases {
            let slots = discount_slots(&structure, &id.parse()?, date)?;
            assert_close!(
                df(&market, &slots)?,
                (-factor * forward_integral).exp(),
                REL,
                ABS
            );
        }
    }
    Ok(())
}

#[test]
fn max_forward_reports_numerical_failure_in_a_candidate() -> Result<()> {
    // 二つ目の変化は +∞ になり、有限な候補が勝つように見えても計算は失敗する。
    let id: CurveId = "MAX".parse()?;
    let spec = CurveSpec::from_expr(
        id.clone(),
        Currency::JPY,
        [CurveRole::Discount],
        CurveExpr::MaxForward {
            parts: vec![
                CurveExpr::Flat(0.02),
                CurveExpr::Scale {
                    factor: f64::MAX,
                    body: Box::new(CurveExpr::Flat(-2.0)),
                },
            ],
            grid: vec![FIRST.parse()?],
        },
    )?;
    let structure = Arc::new(MarketStructure::new(ORIGIN.parse()?, [spec])?);
    let slots = discount_slots(&structure, &id, FIRST)?;
    let market = Market::<f64>::from_theta(structure, vec![])?;
    assert!(matches!(df(&market, &slots), Err(Error::Numerical(_))));
    Ok(())
}

#[test]
fn max_forward_shares_ties_equally_in_aad_and_dual()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // 候補が 2 個のときは scioxide の min、3 個以上のときは同点の平均を走査する分岐を通る。
    // どちらでも同点の m 個に 1/m ずつ配り、選ばれない候補の微分はゼロになることを確かめる。
    // θ は曲線 ID の辞書順に並ぶので、候補の順序を変えても配分の添字は変わらない。
    let first = FIRST.parse()?;
    let maximum: CurveId = "MAX".parse()?;
    let third = 1.0 / 3.0;
    let four = [-0.04_f64, -0.04, -0.04, -0.02];
    let cases: [(&[&str], &[f64], &[f64]); 3] = [
        // A と B が同点で 1/2 ずつ。
        (&["B", "A"], &[-0.04, -0.04], &[0.5, 0.5]),
        // A、B、C が同点で 1/3 ずつ、D は選ばれない。先頭が勝者でない順序も試す。
        (&["A", "B", "C", "D"], &four, &[third, third, third, 0.0]),
        (&["D", "C", "A", "B"], &four, &[third, third, third, 0.0]),
    ];
    for (order, theta, shares) in cases {
        let expected = theta[0].exp();
        let mut specs = order
            .iter()
            .map(|id| {
                CurveSpec::new(
                    id.parse()?,
                    Currency::JPY,
                    [CurveRole::Discount],
                    [first],
                    Extrapolation::Reject,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        specs.push(CurveSpec::from_expr(
            maximum.clone(),
            Currency::JPY,
            [CurveRole::Discount],
            CurveExpr::MaxForward {
                parts: order
                    .iter()
                    .map(|id| id.parse().map(CurveExpr::Ref))
                    .collect::<Result<_>>()?,
                grid: vec![first],
            },
        )?);
        let structure = Arc::new(MarketStructure::new(ORIGIN.parse()?, specs)?);
        let slots = discount_slots(&structure, &maximum, FIRST)?;
        let mut gradient = vec![0.0; theta.len()];
        let value = Tape::default().gradient_into(
            |active| {
                Market::from_theta(Arc::clone(&structure), active.to_vec())
                    .and_then(|market| df(&market, &slots))
                    .map_err(|error| {
                        SciOxideError::invalid_input(
                            argument!("market"),
                            InvalidInputReason::Inconsistent,
                            error.to_string(),
                        )
                    })
            },
            theta,
            &mut gradient,
            ExecutionContext::sequential(),
        )?;
        assert_close!(value, expected, REL, ABS);
        for (i, (derivative, share)) in gradient.into_iter().zip(shares).enumerate() {
            let expected_derivative = expected * share;
            assert_close!(derivative, expected_derivative, REL, ABS);
            let mut dual: Vec<_> = theta.iter().map(|&value| Dual::new(value, 0.0)).collect();
            dual[i] = Dual::variable(theta[i]);
            let market = Market::from_theta(Arc::clone(&structure), dual)?;
            assert_close!(
                df(&market, &slots)?.derivative(),
                expected_derivative,
                REL,
                ABS
            );
        }
    }
    Ok(())
}

#[test]
fn flat_forwards_matches_hand_calculated_discount_factors_and_forwards() -> Result<()> {
    // t = 0.5 の DF は 1 と 1.04 の幾何平均で、t = 1.5 では 1.04 と 0.96 の幾何平均。
    let (structure, id) = structure(Extrapolation::default())?;
    let dates = [
        ORIGIN,
        "2026-07-02T12:00:00Z",
        FIRST,
        "2027-07-02T12:00:00Z",
        SECOND,
    ];
    let points = dates
        .into_iter()
        .map(|at| structure.prepare_point(&id, at.parse()?))
        .collect::<Result<Vec<_>>>()?;
    let tau = YearFraction::new(365.0 / 360.0, DayCount::Act360)?;
    let forward = (points[2].clone(), points[4].clone(), tau);
    let slots = Slots::new(
        Arc::clone(&structure),
        [(id.clone(), points)],
        [(id, vec![forward])],
    )?;
    let market = Market::from_theta(structure, vec![1.04_f64.ln(), 0.96_f64.ln()])?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        let expected = [1.0, 1.04_f64.sqrt(), 1.04, (1.04_f64 * 0.96).sqrt(), 0.96];
        for (&at, expected) in bound.discounts()[0].points().iter().zip(expected) {
            assert_close!(view.df(at).value(), expected, REL, ABS);
        }
        // モデル時間が 1 年でも、予想の分母には指標の ACT/360 を使う。
        let simple = view.forward(bound.projections()[0].points()[0]);
        assert_close!(simple, (1.04 / 0.96 - 1.0) * 360.0 / 365.0, REL, ABS);
        Ok(())
    })
}

#[test]
fn flat_forwards_agree_with_scioxide_linear_spline_of_log_discount_factors()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // FlatForwards の L(t) は (0, 0) と (t_i, θ_i) を結ぶ折れ線で、FlatInstantaneousForward は
    // 最終区間をそのまま延ばす。scioxide-interp の線形 B-spline を MethodDefined
    // （最外の区間の延長）で評価した値と同じになるので、準備した重みの独立な検算に使う。
    let as_of: Instant = ORIGIN.parse()?;
    let id: CurveId = "USD-SOFR".parse()?;
    // 間隔の違うピラー（90 日、1 年、閏年を含む 3 年）で、区間ごとに異なる重みを作る。
    let pillars = ["2026-04-01T00:00:00Z", FIRST, "2029-01-01T00:00:00Z"]
        .into_iter()
        .map(str::parse)
        .collect::<Result<Vec<Instant>>>()?;
    let spec = CurveSpec::new(
        id.clone(),
        Currency::USD,
        [CurveRole::Discount],
        pillars,
        Extrapolation::FlatInstantaneousForward,
    )?;
    let structure = Arc::new(MarketStructure::new(as_of, [spec])?);
    // 原点と最初のピラーの間、ピラー上、ピラーの間、最後のピラーより後の点を読む。
    let queries = [
        "2026-02-15T00:00:00Z",
        "2026-04-01T00:00:00Z",
        "2026-09-01T00:00:00Z",
        FIRST,
        "2028-03-01T06:00:00Z",
        "2029-01-01T00:00:00Z",
        "2031-07-01T00:00:00Z",
    ]
    .into_iter()
    .map(str::parse)
    .collect::<Result<Vec<Instant>>>()?;
    let points = queries
        .iter()
        .map(|&at| structure.prepare_point(&id, at))
        .collect::<Result<Vec<_>>>()?;
    let slots = Slots::new(Arc::clone(&structure), [(id.clone(), points)], [])?;

    let theta = vec![-0.002, -0.011, -0.047];
    let times = structure
        .pillar_times(&id)
        .ok_or("pillar times of the curve")?;
    let nodes: Vec<f64> = std::iter::once(0.0)
        .chain(times.iter().map(|time| time.years()))
        .collect();
    let values: Vec<f64> = std::iter::once(0.0).chain(theta.iter().copied()).collect();
    let spline = BSpline::linear(nodes, values)?;

    let market = Market::from_theta(structure, theta)?;
    let dfs = market.with_view(|view| -> Result<Vec<f64>> {
        let bound = view.bind_slots(&slots)?;
        Ok(bound.discounts()[0]
            .points()
            .iter()
            .map(|&at| view.df(at).value())
            .collect())
    })?;
    for (&at, df) in queries.iter().zip(dfs) {
        let log_df = spline.evaluate(
            [model_time(as_of, at).years()],
            ExtrapolationPolicy::MethodDefined,
        )?;
        assert_close!(df, log_df.exp(), REL, ABS);
    }
    Ok(())
}

#[test]
fn extrapolation_keeps_the_selected_forward_slope() -> Result<()> {
    // L(1) = -0.02、L(2) = -0.06 なら、L(3) は最終区間の延長で -0.10、
    // 原点からの平均の延長で -0.09 になる。
    for (mode, log_df) in [
        (Extrapolation::FlatInstantaneousForward, -0.10_f64),
        (Extrapolation::FlatAverageForward, -0.09_f64),
    ] {
        let (structure, id) = structure(mode)?;
        let slots = discount_slots(&structure, &id, THIRD)?;
        let market = Market::from_theta(structure, vec![-0.02, -0.06])?;
        assert_close!(df(&market, &slots)?, log_df.exp(), REL, ABS);
    }
    let (structure, id) = structure(Extrapolation::Reject)?;
    assert!(matches!(
        structure.prepare_point(&id, THIRD.parse()?),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        structure.prepare_point(&id, "2025-12-31T23:59:59.999999Z".parse()?),
        Err(Error::InvalidInput(_))
    ));
    let last = discount_slots(&structure, &id, SECOND)?;
    let market = Market::from_theta(structure, vec![-0.02, -0.06])?;
    assert_close!(df(&market, &last)?, (-0.06_f64).exp(), REL, ABS);
    Ok(())
}

#[test]
fn instantaneous_extrapolation_uses_the_origin_for_one_pillar() -> Result<()> {
    // 1 年に唯一の節点があるとき、最終区間の始点は原点。2 年では L(2) = 2 L(1)。
    let (structure, id, _) = two_curves()?;
    let slots = discount_slots(&structure, &id, SECOND)?;
    let market = Market::from_theta(structure, vec![0.98_f64.ln(), -0.03])?;
    assert_close!(df(&market, &slots)?, 0.98_f64.powi(2), REL, ABS);
    Ok(())
}

#[test]
fn flat_last_segment_stays_flat_when_extrapolated_from_close_pillars() -> Result<()> {
    // 1 µs しか離れていない節点に同じ ln DF があれば、最終区間の傾きはゼロ。
    let id: CurveId = "JPY-TONA".parse()?;
    let spec = CurveSpec::new(
        id.clone(),
        Currency::JPY,
        [CurveRole::Discount],
        [FIRST.parse()?, "2027-01-01T00:00:00.000001Z".parse()?],
        Extrapolation::default(),
    )?;
    let structure = Arc::new(MarketStructure::new(ORIGIN.parse()?, [spec])?);
    let slots = discount_slots(&structure, &id, THIRD)?;
    let log_df = 0.98_f64.ln();
    let market = Market::from_theta(structure, vec![log_df, log_df])?;
    assert_close!(df(&market, &slots)?, 0.98, REL, ABS);
    Ok(())
}

#[test]
fn forward_retains_small_log_discount_differences() -> Result<()> {
    // δ = 1e-17 では exp(δ) の丸め結果は 1 だが、τ = 1e-6 の単利は約 1e-11。
    let (structure, id) = structure(Extrapolation::default())?;
    let start = structure.prepare_point(&id, ORIGIN.parse()?)?;
    let end = structure.prepare_point(&id, FIRST.parse()?)?;
    let tau = YearFraction::new(1e-6, DayCount::Act365F)?;
    let slots = Slots::new(Arc::clone(&structure), [], [(id, vec![(start, end, tau)])])?;
    let market = Market::from_theta(structure, vec![-1e-17, -0.06])?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        assert_close!(
            view.forward(bound.projections()[0].points()[0]),
            1e-11,
            REL,
            ABS
        );
        Ok(())
    })
}

#[test]
fn growth_log_retains_small_log_discount_differences() -> Result<()> {
    // δ = 1e-17 の DF 比は丸めると 1 になり、1 を引いた超過は 0 になる。
    // 対数のまま読めば δ が残り、exp_m1 で超過 1e-17 を倍精度の精度で得られる。
    let (structure, id) = structure(Extrapolation::default())?;
    let start = structure.prepare_point(&id, ORIGIN.parse()?)?;
    let end = structure.prepare_point(&id, FIRST.parse()?)?;
    let slots = Slots::from_resolved(Arc::clone(&structure), [], [])?
        .with_projection_growth([(structure.resolve_projection(&id)?, vec![(start, end)])])?;
    let market = Market::from_theta(structure, vec![-1e-17, -0.06])?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        let at = bound.projection_growth()[0].points()[0];
        assert_eq!(view.growth(at), 1.0);
        assert_close!(view.log_growth(at), 1e-17, REL, 0.0);
        assert_close!(view.log_growth(at).exp_m1(), 1e-17, REL, 0.0);
        Ok(())
    })
}

#[test]
fn slots_rebind_on_the_same_structure_with_new_theta() -> Result<()> {
    // 同じ点とスロットから、毎回その市場の θ の値を読む。
    let (structure, id) = structure(Extrapolation::default())?;
    let slots = discount_slots(&structure, &id, FIRST)?;
    let first = Market::from_theta(Arc::clone(&structure), vec![0.98_f64.ln(), -0.06])?;
    let second = Market::from_theta(structure, vec![0.97_f64.ln(), -0.06])?;
    assert_close!(df(&first, &slots)?, 0.98, REL, ABS);
    assert_close!(df(&second, &slots)?, 0.97, REL, ABS);
    Ok(())
}

#[test]
fn curve_reads_preserve_dual_sensitivities() -> Result<()> {
    // 中点 L = (θ₁ + θ₂)/2 の DF の θ₁ 微分は DF/2。
    // F(1, 2) の θ₁ 微分は DF(1)/DF(2) ÷ τ。
    let (structure, id) = structure(Extrapolation::default())?;
    let middle = structure.prepare_point(&id, "2027-07-02T12:00:00Z".parse()?)?;
    let start = structure.prepare_point(&id, FIRST.parse()?)?;
    let end = structure.prepare_point(&id, SECOND.parse()?)?;
    let tau = YearFraction::new(1.0, DayCount::Act365F)?;
    let slots = Slots::new(
        Arc::clone(&structure),
        [(id.clone(), vec![middle])],
        [(id, vec![(start, end, tau)])],
    )?;
    let market = Market::from_theta(
        structure,
        vec![Dual::new(0.98_f64.ln(), 1.0), Dual::new(0.94_f64.ln(), 0.0)],
    )?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        let df = view.df(bound.discounts()[0].points()[0]).value();
        let forward = view.forward(bound.projections()[0].points()[0]);
        assert_close!(df.derivative(), (0.98_f64 * 0.94).sqrt() / 2.0, REL, ABS);
        assert_close!(forward.derivative(), 0.98 / 0.94, REL, ABS);
        Ok(())
    })
}

#[test]
fn currency_scope_preserves_the_aad_gradient() -> std::result::Result<(), Box<dyn std::error::Error>>
{
    // 1 年の節点で支払う 100 円の PV と θ₁ 微分はともに 100·exp(θ₁)。
    // 支払点が読まない 2 年の節点の微分はゼロ。
    let (structure, id) = structure(Extrapolation::default())?;
    let point = structure.prepare_point(&id, FIRST.parse()?)?;
    let mut gradient = [0.0; 2];
    let pv = Tape::default().gradient_into(
        |theta| {
            let result = (|| {
                let market = Market::from_theta(Arc::clone(&structure), theta.to_vec())?;
                market.with_view(|view| {
                    view.ccy_scope(Currency::JPY, view.discount_curve(&id)?, |scope| {
                        let value = scope.pv(scioxide_ad::Real::constant(100.0), scope.at(&point)?);
                        Ok(scope.into_money(value).amount())
                    })
                })
            })();
            result.map_err(|error: Error| {
                SciOxideError::invalid_input(
                    argument!("market"),
                    InvalidInputReason::Inconsistent,
                    error.to_string(),
                )
            })
        },
        &[0.98_f64.ln(), -0.06],
        &mut gradient,
        ExecutionContext::sequential(),
    )?;
    assert_close!(pv, 98.0, REL, 1e-6);
    assert_close!(gradient[0], 98.0, REL, 1e-6);
    assert_close!(gradient[1], 0.0, REL, 1e-6);
    Ok(())
}

#[test]
fn bind_rejects_point_of_other_curve() -> Result<()> {
    // 節点の日程と通貨が同じでも、予想曲線の点を割引に渡すと拒否する。
    let (structure, discount, projection) = two_curves()?;
    let start = structure.prepare_point(&projection, ORIGIN.parse()?)?;
    let end = structure.prepare_point(&projection, FIRST.parse()?)?;
    let other = structure.prepare_point(&discount, FIRST.parse()?)?;
    let invalid_discount = Slots::from_resolved(
        Arc::clone(&structure),
        [(structure.resolve_discount(&discount)?, vec![end.clone()])],
        [],
    )?;
    let invalid_growth = Slots::from_resolved(Arc::clone(&structure), [], [])?
        .with_projection_growth([(
            structure.resolve_projection(&projection)?,
            vec![(start.clone(), other.clone())],
        )])?;
    let invalid_forward = Slots::new(
        Arc::clone(&structure),
        [],
        [(
            projection,
            vec![(start, other, YearFraction::new(1.0, DayCount::Act360)?)],
        )],
    )?;
    let market = Market::from_theta(structure, vec![-0.02, -0.03])?;
    market.with_view(|view| {
        assert!(matches!(
            view.bind_slots(&invalid_discount),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            view.bind_slots(&invalid_forward),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            view.bind_slots(&invalid_growth),
            Err(Error::InvalidInput(_))
        ));
        view.ccy_scope(Currency::USD, view.discount_curve(&discount)?, |scope| {
            assert!(matches!(scope.at(&end), Err(Error::InvalidInput(_))));
            Ok(())
        })
    })
}

#[test]
fn binding_checks_structure_identity_even_for_identical_definitions() -> Result<()> {
    // 曲線 ID や定義が一致しても、別の構造で作った点やスロットは使えない。
    let mode = Extrapolation::default();
    let (current, id) = structure(mode)?;
    let (other, _) = structure(mode)?;
    let foreign_curve = other.resolve_discount(&id)?;
    assert!(matches!(
        Slots::from_resolved(Arc::clone(&current), [(foreign_curve, vec![])], []),
        Err(Error::InvalidInput(_))
    ));
    let foreign_slots = discount_slots(&other, &id, FIRST)?;
    let foreign_point = other.prepare_point(&id, FIRST.parse()?)?;
    let mixed = Slots::new(
        Arc::clone(&current),
        [(id.clone(), vec![foreign_point.clone()])],
        [],
    )?;
    let market = Market::from_theta(current, vec![-0.02, -0.06])?;
    market.with_view(|view| {
        for slots in [&foreign_slots, &mixed] {
            assert!(matches!(
                view.bind_slots(slots),
                Err(Error::InvalidInput(_))
            ));
        }
        view.ccy_scope(Currency::JPY, view.discount_curve(&id)?, |scope| {
            assert!(matches!(
                scope.at(&foreign_point),
                Err(Error::InvalidInput(_))
            ));
            Ok(())
        })
    })
}

#[test]
fn curve_resolution_checks_role_and_market_availability() -> Result<()> {
    // 利用できない予想ブロックがあっても、利用可能な割引スロットは読める。
    let (structure, discount, projection) = two_curves()?;
    let discounts = discount_slots(&structure, &discount, FIRST)?;
    let projections = Slots::new(Arc::clone(&structure), [], [(projection.clone(), vec![])])?;
    assert!(matches!(
        Slots::new(Arc::clone(&structure), [(projection.clone(), vec![])], []),
        Err(Error::InvalidInput(_))
    ));
    let partial = Market::from_partial_theta(
        structure,
        vec![0.97_f64.ln(), f64::NAN],
        [BlockId::nodes(discount.clone())],
    )?;
    let absent = "USD-OTHER".parse()?;
    partial.with_view(|view| {
        assert!(matches!(
            view.discount_curve(&absent),
            Err(Error::MissingData(_))
        ));
        assert!(matches!(
            view.projection_curve(&discount),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            view.projection_curve(&projection),
            Err(Error::MissingData(_))
        ));
        assert!(matches!(
            view.bind_slots(&projections),
            Err(Error::MissingData(_))
        ));
        let bound = view.bind_slots(&discounts)?;
        assert_close!(
            view.df(bound.discounts()[0].points()[0]).value(),
            0.97,
            REL,
            ABS
        );
        Ok(())
    })
}

#[test]
fn currency_scope_discounts_and_tags_present_values() -> Result<()> {
    // 各支払いを割り引き、計算済みの現在価値と足し引きして通貨付きの金額を得る。
    let (structure, id, _) = two_curves()?;
    let point = structure.prepare_point(&id, FIRST.parse()?)?;
    let slots = Slots::new(Arc::clone(&structure), [(id, vec![])], [])?;
    let market = Market::from_theta(structure, vec![0.96_f64.ln(), -0.03])?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        let curve = bound.discounts()[0].curve();
        assert!(matches!(
            view.ccy_scope(Currency::JPY, curve, |_| Ok(())),
            Err(Error::InvalidInput(_))
        ));
        view.ccy_scope(Currency::USD, curve, |scope| {
            let pay = scope.at(&point)?;
            let amount =
                scope.pv(100.0, pay) - scope.pv(2.0, pay).scale(3.0) + -scope.present(-5.0);
            let money = scope.into_money(amount);
            assert_eq!(money.ccy(), Currency::USD);
            assert_close!(money.amount(), 95.24, REL, 1e-8);
            Ok(())
        })
    })
}

#[test]
fn bind_rejects_nonpositive_accrual_and_invalid_forward_periods() -> Result<()> {
    // 同時刻と逆順の期間、ゼロと負の日数係数を予想の境界で拒否する。
    let (structure, id) = structure(Extrapolation::default())?;
    let start = structure.prepare_point(&id, ORIGIN.parse()?)?;
    let end = structure.prepare_point(&id, FIRST.parse()?)?;
    let market = Market::from_theta(Arc::clone(&structure), vec![0.0, 0.0])?;
    for (start, end, tau) in [
        (&start, &start, 1.0),
        (&end, &start, 1.0),
        (&start, &end, 0.0),
        (&start, &end, -1.0),
    ] {
        let points = (
            start.clone(),
            end.clone(),
            YearFraction::new(tau, DayCount::Act365F)?,
        );
        let slots = Slots::new(Arc::clone(&structure), [], [(id.clone(), vec![points])])?;
        market.with_view(|view| {
            assert!(matches!(
                view.bind_slots(&slots),
                Err(Error::InvalidInput(_))
            ));
        });
    }
    Ok(())
}

#[test]
fn binding_reports_unrepresentable_discount_factors_and_forwards() -> Result<()> {
    // θ が有限でも exp のオーバーフローやゼロへのアンダーフローはあり得る。
    let (structure, id) = structure(Extrapolation::default())?;
    let discount = discount_slots(&structure, &id, FIRST)?;
    let point = structure.prepare_point(&id, FIRST.parse::<Instant>()?)?;
    for log_df in [1000.0, -1000.0] {
        let market = Market::from_theta(Arc::clone(&structure), vec![log_df, -0.06])?;
        market.with_view(|view| {
            assert!(matches!(
                view.bind_slots(&discount),
                Err(Error::Numerical(_))
            ));
            view.ccy_scope(Currency::JPY, view.discount_curve(&id)?, |scope| {
                assert!(matches!(scope.at(&point), Err(Error::Numerical(_))));
                Ok(())
            })
        })?;
    }
    let start = structure.prepare_point(&id, ORIGIN.parse()?)?;
    let tau = YearFraction::new(1.0, DayCount::Act365F)?;
    let projection = Slots::new(
        Arc::clone(&structure),
        [],
        [(id, vec![(start, point, tau)])],
    )?;
    let market = Market::from_theta(structure, vec![-1000.0, -0.06])?;
    market.with_view(|view| {
        assert!(matches!(
            view.bind_slots(&projection),
            Err(Error::Numerical(_))
        ));
    });
    Ok(())
}

#[test]
fn bound_pay_rejects_another_discount_curve_of_same_currency() -> Result<()> {
    // 型と通貨が同じ割引点についても、結合した曲線がスコープと一致する必要がある。
    let ids = ["first".parse::<CurveId>()?, "second".parse::<CurveId>()?];
    let at = FIRST.parse()?;
    let specs = ids
        .iter()
        .map(|id| {
            CurveSpec::new(
                id.clone(),
                Currency::JPY,
                [CurveRole::Discount],
                [at],
                Extrapolation::default(),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let structure = Arc::new(MarketStructure::new(ORIGIN.parse()?, specs)?);
    let slots = discount_slots(&structure, &ids[0], FIRST)?;
    let market = Market::from_theta(structure, vec![-0.02, -0.03])?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        let at = bound.discounts()[0].points()[0];
        view.ccy_scope(Currency::JPY, view.discount_curve(&ids[1])?, |scope| {
            assert!(matches!(scope.bind_pay(at), Err(Error::InvalidInput(_))));
            Ok(())
        })?;
        view.ccy_scope(Currency::JPY, bound.discounts()[0].curve(), |scope| {
            let pv = scope
                .into_money(scope.pv(100.0, scope.bind_pay(at)?))
                .amount();
            assert_close!(pv, 100.0 * (-0.02_f64).exp(), REL, 1e-6);
            Ok(())
        })
    })
}

const MEETING: &str = "2026-07-02T12:00:00Z";

fn composed_curve_structure() -> Result<(Arc<MarketStructure>, [CurveId; 3])> {
    use rhodium_market::curve::CurveExpr;

    let as_of = ORIGIN.parse()?;
    let meeting = MEETING.parse()?;
    let boundary = FIRST.parse()?;
    let later = SECOND.parse()?;
    let tona: CurveId = "JPY-TONA".parse()?;
    let tibor: CurveId = "JPY-TIBOR-6M".parse()?;
    let sofr: CurveId = "USD-SOFR".parse()?;
    let structure = Arc::new(MarketStructure::new(
        as_of,
        [
            CurveSpec::from_expr(
                tona.clone(),
                Currency::JPY,
                [CurveRole::Discount, CurveRole::Projection],
                CurveExpr::Sum(vec![
                    CurveExpr::Stitch {
                        short: Box::new(CurveExpr::Steps {
                            part: "meetings".into(),
                            meetings: vec![meeting],
                        }),
                        long: Box::new(CurveExpr::Nodes {
                            part: "long".into(),
                            pillars: vec![later],
                            extrapolation: Extrapolation::default(),
                        }),
                        at: boundary,
                    },
                    CurveExpr::Jumps {
                        part: "turn".into(),
                        dates: vec![boundary],
                    },
                ]),
            )?,
            CurveSpec::from_expr(
                tibor.clone(),
                Currency::JPY,
                [CurveRole::Projection],
                CurveExpr::Sum(vec![
                    CurveExpr::Ref(tona.clone()),
                    CurveExpr::Nodes {
                        part: "basis".into(),
                        pillars: vec![later],
                        extrapolation: Extrapolation::default(),
                    },
                ]),
            )?,
            CurveSpec::from_expr(
                sofr.clone(),
                Currency::USD,
                [CurveRole::Discount],
                CurveExpr::Flat(0.04),
            )?,
        ],
    )?);
    Ok((structure, [tona, tibor, sofr]))
}

#[test]
fn composed_curves_match_manual_log_discount_factors() -> Result<()> {
    // 会合の前後で傾きが変わり、ターンは発効時刻から効く。長期部分は接続点の値を引く。
    let (structure, [tona, tibor, sofr]) = composed_curve_structure()?;
    let as_of = structure.as_of();
    let meeting = MEETING.parse()?;
    let boundary = FIRST.parse()?;
    let later = SECOND.parse()?;
    let market = Market::from_theta(
        Arc::clone(&structure),
        vec![-0.004, -0.06, 0.01, 0.02, 0.001],
    )?;
    let before: Instant = "2026-12-31T23:59:59Z".parse()?;
    let cases = [
        (meeting, -0.005),
        (
            before,
            -0.005 - 0.02 * (model_time(as_of, before).years() - 0.5),
        ),
        (boundary, -0.016),
        (later, -0.046),
    ];
    for (at, log_expected) in cases {
        assert_close!(
            df(
                &market,
                &discount_slots(&structure, &tona, &at.to_string())?
            )?,
            log_expected.exp(),
            REL,
            ABS
        );
    }
    assert_close!(
        df(&market, &discount_slots(&structure, &sofr, SECOND)?)?,
        (-0.08_f64).exp(),
        REL,
        ABS
    );

    let start = structure.prepare_point(&tibor, boundary)?;
    let end = structure.prepare_point(&tibor, later)?;
    let slots = Slots::new(
        Arc::clone(&structure),
        [],
        [(
            tibor,
            vec![(start, end, YearFraction::new(1.0, DayCount::Act365F)?)],
        )],
    )?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        assert_close!(
            view.forward(bound.projections()[0].points()[0]),
            0.032_f64.exp_m1(),
            REL,
            ABS
        );
        Ok(())
    })?;

    // 参照先の長期ブロックが欠けても、短期の準備済みの点には必要ない。
    let partial = Market::from_partial_theta(
        Arc::clone(&structure),
        market.theta().to_vec(),
        [BlockId::new(tona.clone(), "meetings")?],
    )?;
    assert_close!(
        df(
            &partial,
            &discount_slots(&structure, &tona, &meeting.to_string())?
        )?,
        (-0.005_f64).exp(),
        REL,
        ABS
    );
    assert!(matches!(
        partial.with_view(|view| view.bind_slots(&slots).map(|_| ())),
        Err(Error::MissingData(_))
    ));
    Ok(())
}
`````

`````
crates/rhodium-market/tests/quotes.rs
//! quote の識別子とスナップショットの入力境界を検証する。

use rhodium_core::{Error, Result};
use rhodium_market::quote::{
    QuoteId, QuoteKind, QuoteOverride, QuoteOverrides, QuotePattern, QuoteSet,
};
use serde::{Deserialize, de::value::StrDeserializer};

#[test]
fn quote_id_parses_and_deserializes_with_the_same_grammar() -> Result<()> {
    // + と - は区間の一部になり、表記の大小や区切りは正規化しない。
    let valid = "JPY.TONA-SOFR.ATM+50BP";
    let parsed: QuoteId = valid.parse()?;
    let deserialize =
        |text| QuoteId::deserialize(StrDeserializer::<serde::de::value::Error>::new(text));
    assert_eq!(parsed.to_string(), valid);
    assert_eq!(deserialize(valid).ok(), Some(parsed));
    for invalid in [
        "",
        "JPY",
        "jpy.TONA",
        ".TONA",
        "JPY.",
        "JPY..TONA",
        "JPY.TON A",
        "JPY.円",
    ] {
        assert!(matches!(
            invalid.parse::<QuoteId>(),
            Err(Error::InvalidInput(_))
        ));
        assert!(deserialize(invalid).is_err());
    }
    Ok(())
}

#[test]
fn quote_set_keeps_finite_decimal_values_in_id_order() -> Result<()> {
    // 符号の制約は校正の binding が決める。スナップショットは負の quote を保持できる。
    let jpy: QuoteId = "JPY.TONA.OIS.1Y".parse()?;
    let usd: QuoteId = "USD.SOFR.OIS.1Y".parse()?;
    let quotes = QuoteSet::new([(usd.clone(), 0.0375), (jpy.clone(), -0.0005)])?;
    assert_eq!(quotes.get(&jpy), Some(-0.0005));
    assert_eq!(quotes.get(&usd), Some(0.0375));
    assert_eq!(quotes.get(&"JPY.TONA.OIS.2Y".parse()?), None);
    assert_eq!(
        quotes.iter().map(|(id, _)| id).collect::<Vec<_>>(),
        [&jpy, &usd]
    );
    Ok(())
}

#[test]
fn quote_set_rejects_duplicate_and_non_finite_entries() -> Result<()> {
    // 欠落や上書きとは区別して、入力列の不正をその場で返す。
    let id: QuoteId = "USD.SOFR.OIS.1Y".parse()?;
    assert!(matches!(
        QuoteSet::new([(id.clone(), 0.04), (id.clone(), 0.04)]),
        Err(Error::InvalidInput(_))
    ));
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            QuoteSet::new([(id.clone(), value)]),
            Err(Error::InvalidInput(_))
        ));
    }
    Ok(())
}

#[test]
fn quote_patterns_match_whole_segments() -> Result<()> {
    // `**` は追加区間がない場合にも一致し、`*` はちょうど一区間を消費する。
    let one_year = "JPY.TONA.OIS.1Y".parse()?;
    let two_year = "JPY.TONA.OIS.2Y".parse()?;
    let prefix = "JPY.TONA.OIS".parse()?;
    let short = "JPY.TONA.OIS.1".parse()?;
    let alternatives: QuotePattern = "JPY.*.OIS.{1Y,2Y}".parse()?;
    assert!(alternatives.matches(&one_year));
    assert!(alternatives.matches(&two_year));
    assert!(!alternatives.matches(&prefix));
    assert!(!alternatives.matches(&short));
    let tail: QuotePattern = "JPY.TONA.OIS.**".parse()?;
    assert!(tail.matches(&one_year));
    assert!(tail.matches(&prefix));
    let exact: QuotePattern = "JPY.TONA.OIS.1Y".parse()?;
    assert!(exact.matches(&one_year));
    assert!(!exact.matches(&two_year));
    for invalid in [
        "",
        "JPY.**.OIS",
        "JPY..OIS",
        "JPY.TON*",
        "JPY.{A,}",
        "JPY.{}",
        "JPY.{A,b}",
    ] {
        assert!(matches!(
            invalid.parse::<QuotePattern>(),
            Err(Error::InvalidInput(_))
        ));
    }
    Ok(())
}

#[test]
fn quote_overrides_apply_in_order_and_report_net_changes() -> Result<()> {
    // 追加した値は後続のパターンに入り、元に戻した値や追加後に削除した ID は差分に入らない。
    let one: QuoteId = "JPY.TONA.OIS.1Y".parse()?;
    let two: QuoteId = "JPY.TONA.OIS.2Y".parse()?;
    let three: QuoteId = "JPY.TONA.OIS.3Y".parse()?;
    let four: QuoteId = "JPY.TONA.OIS.4Y".parse()?;
    let price: QuoteId = "JPY.ASSET.PRICE".parse()?;
    let kinds = [
        (one.clone(), QuoteKind::Rate),
        (two.clone(), QuoteKind::Rate),
        (three.clone(), QuoteKind::Rate),
        (four.clone(), QuoteKind::Rate),
        (price.clone(), QuoteKind::Price),
    ]
    .into();
    let base = QuoteSet::new([
        (one.clone(), 0.125),
        (two.clone(), 0.25),
        (price.clone(), 100.0),
    ])?;
    let overrides = QuoteOverrides::new(
        [
            QuoteOverride::Set {
                id: three.clone(),
                value: 0.125,
            },
            QuoteOverride::Add {
                pattern: "JPY.TONA.OIS.*".parse()?,
                shift: -0.25,
            },
            QuoteOverride::Set {
                id: one.clone(),
                value: 0.125,
            },
            QuoteOverride::Remove { id: two.clone() },
            QuoteOverride::Set {
                id: four.clone(),
                value: 0.25,
            },
            QuoteOverride::Remove { id: four.clone() },
            QuoteOverride::Scale {
                pattern: "JPY.ASSET.PRICE".parse()?,
                factor: 1.5,
            },
        ],
        &kinds,
    )?;
    let (updated, changed) = base.with_overrides(&overrides)?;
    assert_eq!(changed, [two.clone(), three.clone(), price.clone()].into());
    assert_eq!(updated.get(&one), Some(0.125));
    assert_eq!(updated.get(&two), None);
    assert_eq!(updated.get(&three), Some(-0.125));
    assert_eq!(updated.get(&price), Some(150.0));
    assert_eq!(base.get(&two), Some(0.25));
    assert_eq!(base.get(&three), None);
    Ok(())
}

#[test]
fn quote_overrides_reject_invalid_targets_and_transforms() -> Result<()> {
    // 有限な引数でも、適用時のオーバーフローは入力のエラーになる。
    let rate: QuoteId = "JPY.TONA.OIS.1Y".parse()?;
    let price: QuoteId = "JPY.ASSET.PRICE".parse()?;
    let unknown: QuoteId = "JPY.UNKNOWN".parse()?;
    let kinds = [
        (rate.clone(), QuoteKind::Rate),
        (price.clone(), QuoteKind::Price),
    ]
    .into();
    let base = QuoteSet::new([(price.clone(), f64::MAX)])?;
    for operation in [
        QuoteOverride::Remove {
            id: unknown.clone(),
        },
        QuoteOverride::Set {
            id: rate.clone(),
            value: f64::NAN,
        },
        QuoteOverride::Add {
            pattern: "USD.**".parse()?,
            shift: 0.0,
        },
        QuoteOverride::Add {
            pattern: "JPY.TONA.**".parse()?,
            shift: f64::INFINITY,
        },
        QuoteOverride::Scale {
            pattern: "JPY.*.**".parse()?,
            factor: 1.1,
        },
        QuoteOverride::Scale {
            pattern: "JPY.ASSET.*".parse()?,
            factor: 0.0,
        },
    ] {
        assert!(matches!(
            QuoteOverrides::new([operation], &kinds),
            Err(Error::InvalidInput(_))
        ));
    }
    for operation in [
        QuoteOverride::Add {
            pattern: "JPY.TONA.**".parse()?,
            shift: 0.0,
        },
        QuoteOverride::Add {
            pattern: "JPY.ASSET.*".parse()?,
            shift: f64::MAX,
        },
        QuoteOverride::Scale {
            pattern: "JPY.ASSET.*".parse()?,
            factor: 2.0,
        },
    ] {
        let overrides = QuoteOverrides::new([operation], &kinds)?;
        assert!(matches!(
            base.with_overrides(&overrides),
            Err(Error::InvalidInput(_))
        ));
    }
    let broad_add = QuoteOverrides::new(
        [QuoteOverride::Add {
            pattern: "JPY.**".parse()?,
            shift: 0.0,
        }],
        &kinds,
    )?;
    assert!(matches!(
        QuoteSet::new([(unknown, 0.0)])?.with_overrides(&broad_add),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}
`````

`````
crates/rhodium-market/tests/structure.rs
//! 曲線の ID と定義、評価時点を原点とする節点、平坦な θ の境界を検証する。

use std::sync::Arc;

use rhodium_core::{Error, Result, currency::Currency, time::Instant};
use rhodium_market::{
    curve::{CurveExpr, CurveId, CurveRole, CurveSpec, Extrapolation},
    market::Market,
    structure::{BlockId, MarketStructure, ParamId},
};
use serde::{Deserialize, de::value::StrDeserializer};

fn curve(name: &str, ccy: Currency, dates: &[&str]) -> Result<CurveSpec> {
    CurveSpec::new(
        name.parse()?,
        ccy,
        [CurveRole::Discount, CurveRole::Projection],
        dates
            .iter()
            .map(|date| date.parse())
            .collect::<Result<Vec<Instant>>>()?,
        Extrapolation::default(),
    )
}

#[test]
fn curve_id_parses_and_deserializes_with_the_same_grammar() -> Result<()> {
    // 担保付きの曲線名のような `/` と `:` を含む ID も受け付け、表記はそのまま保つ。
    // 文字列の解析と serde は同じ検査を通るので、片方だけが通る ID はない。
    let deserialize =
        |text| CurveId::deserialize(StrDeserializer::<serde::de::value::Error>::new(text));
    let longest = "A".repeat(128);
    for valid in [
        "JPY-TONA",
        "USD/csa:JPY-TONA",
        "usd.sofr_1+x",
        longest.as_str(),
    ] {
        let parsed: CurveId = valid.parse()?;
        assert_eq!(parsed.to_string(), valid);
        assert_eq!(deserialize(valid).ok(), Some(parsed));
    }
    // 空、129 バイト、空白、許さない記号、非 ASCII を拒否する。
    let too_long = "A".repeat(129);
    for invalid in ["", too_long.as_str(), "JPY TONA", "JPY-TONA#1", "JPY-円"] {
        assert!(matches!(
            invalid.parse::<CurveId>(),
            Err(Error::InvalidInput(_))
        ));
        assert!(deserialize(invalid).is_err());
    }
    Ok(())
}

#[test]
fn curve_spec_requires_at_least_one_role() -> Result<()> {
    // 役割の無い曲線はどのハンドルからも読めないので、定義の時点で拒否する。
    let spec = CurveSpec::new(
        "JPY-TONA".parse()?,
        Currency::JPY,
        [],
        ["2027-01-01T00:00:00Z".parse()?],
        Extrapolation::default(),
    );
    assert!(matches!(spec, Err(Error::InvalidInput(_))));
    Ok(())
}

#[test]
fn theta_layout_orders_curves_and_excludes_the_origin() -> Result<()> {
    // 定義を逆順に渡しても、JPY の 2 点、USD の 1 点の順に並ぶ。
    let as_of = "2026-01-01T00:00:00Z".parse()?;
    let usd = curve("USD-SOFR", Currency::USD, &["2027-01-01T00:00:00Z"])?;
    let jpy = curve(
        "JPY-TONA",
        Currency::JPY,
        &["2027-01-01T00:00:00Z", "2028-01-01T00:00:00Z"],
    )?;
    let jpy_block = BlockId::nodes(jpy.id().clone());
    let usd_block = BlockId::nodes(usd.id().clone());
    let structure = MarketStructure::new(as_of, [usd, jpy])?;
    let layout = structure.layout();
    assert_eq!(layout.len(), 3);
    assert_eq!(layout.range(&jpy_block), Some(0..2));
    assert_eq!(layout.range(&usd_block), Some(2..3));
    assert_eq!(layout.index(&ParamId::new(jpy_block.clone(), 1)), Some(1));
    assert_eq!(layout.index(&ParamId::new(usd_block.clone(), 0)), Some(2));
    assert_eq!(layout.index(&ParamId::new(jpy_block, 2)), None);
    assert_eq!(layout.index(&ParamId::new(usd_block, u32::MAX)), None);
    Ok(())
}

#[test]
fn market_structure_rejects_invalid_pillars_and_duplicate_ids() -> Result<()> {
    // 原点は暗黙の節点なので、評価時点やそれ以前を陽に置けない。
    let origin = "2026-01-01T01:00:00Z";
    let as_of = origin.parse()?;
    let early = "2026-01-01T00:59:59.999999Z";
    let first = "2026-01-01T01:00:00.000001Z";
    let second = "2027-01-01T00:00:00Z";
    for invalid in [&[][..], &[first, first], &[second, first]] {
        assert!(matches!(
            curve("JPY-TONA", Currency::JPY, invalid),
            Err(Error::InvalidInput(_))
        ));
    }
    for invalid in [origin, early] {
        let spec = curve("JPY-TONA", Currency::JPY, &[invalid, first])?;
        assert!(matches!(
            MarketStructure::new(as_of, [spec]),
            Err(Error::InvalidInput(_))
        ));
    }
    // UTC の瞬間として異なっても、分解能を失うほど遠い節点どうしは区間を作れない。
    let close = curve(
        "JPY-TONA",
        Currency::JPY,
        &["9999-01-01T00:00:00Z", "9999-01-01T00:00:00.000001Z"],
    )?;
    assert!(matches!(
        MarketStructure::new(as_of, [close]),
        Err(Error::InvalidInput(_))
    ));
    let one = curve("JPY-TONA", Currency::JPY, &[first])?;
    let two = curve("JPY-TONA", Currency::JPY, &[second])?;
    assert!(matches!(
        MarketStructure::new(as_of, [one, two]),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn max_forward_and_scale_validate_curve_definitions() -> Result<()> {
    // 時刻の順序と候補の有無は曲線の定義時に、評価時点との関係は構造の構築時に検査する。
    let as_of = "2026-01-01T00:00:00Z".parse()?;
    let first = "2027-01-01T00:00:00Z".parse()?;
    let second = "2028-01-01T00:00:00Z".parse()?;
    let max = |parts, grid| {
        CurveSpec::from_expr(
            "MAX".parse()?,
            Currency::JPY,
            [CurveRole::Discount],
            CurveExpr::MaxForward { parts, grid },
        )
    };
    let candidates = || vec![CurveExpr::Flat(0.02)];
    for spec in [
        max(vec![], vec![first]),
        max(candidates(), vec![]),
        max(candidates(), vec![first, first]),
        max(candidates(), vec![second, first]),
    ] {
        assert!(matches!(spec, Err(Error::InvalidInput(_))));
    }
    let close = [
        "9999-01-01T00:00:00Z".parse()?,
        "9999-01-01T00:00:00.000001Z".parse()?,
    ];
    for grid in [
        vec!["2025-12-31T23:59:59Z".parse()?, first],
        vec![as_of, first],
        close.to_vec(),
    ] {
        let spec = max(candidates(), grid)?;
        assert!(matches!(
            MarketStructure::new(as_of, [spec]),
            Err(Error::InvalidInput(_))
        ));
    }
    let scale = |factor, body| CurveExpr::Scale {
        factor,
        body: Box::new(body),
    };
    assert!(matches!(
        max(
            vec![scale(f64::INFINITY, CurveExpr::Flat(0.02))],
            vec![first]
        ),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        max(vec![scale(0.5, CurveExpr::Flat(f64::NAN))], vec![first]),
        Err(Error::InvalidInput(_))
    ));

    // 候補の内側で付けた部品名も、曲線全体で一意である必要がある。
    let nodes = CurveExpr::Nodes {
        part: "nodes".into(),
        pillars: vec![first],
        extrapolation: Extrapolation::Reject,
    };
    assert!(matches!(
        max(vec![nodes.clone(), scale(0.5, nodes)], vec![first]),
        Err(Error::InvalidInput(message)) if message.contains("duplicate curve part")
    ));
    Ok(())
}

#[test]
fn nonlinear_curves_check_references_and_projection_roles() -> Result<()> {
    // 候補の内側の参照も通常の循環・欠損チェックを通す。
    let as_of = "2026-01-01T00:00:00Z".parse()?;
    let first = "2027-01-01T00:00:00Z".parse()?;
    let composite = |body| {
        CurveSpec::from_expr(
            "USD-MAX".parse()?,
            Currency::USD,
            [CurveRole::Discount],
            CurveExpr::MaxForward {
                parts: vec![CurveExpr::Scale {
                    factor: -1.0,
                    body: Box::new(body),
                }],
                grid: vec![first],
            },
        )
    };
    for (source, diagnostic) in [("MISSING", "unknown"), ("USD-MAX", "cyclic")] {
        let spec = composite(CurveExpr::Ref(source.parse()?))?;
        assert!(matches!(
            MarketStructure::new(as_of, [spec]),
            Err(Error::InvalidData(message)) if message.contains(diagnostic)
        ));
    }

    // 通常の参照は外貨を拒否し、付利のフォワード参照には予想の役割が要る。
    let id: CurveId = "JPY-TONA".parse()?;
    let jpy = |role| CurveSpec::from_expr(id.clone(), Currency::JPY, [role], CurveExpr::Flat(0.01));
    for (reference, role, diagnostic) in [
        (
            CurveExpr::Ref(id.clone()),
            CurveRole::Discount,
            "invalid curve reference",
        ),
        (
            CurveExpr::ForwardRef(id.clone()),
            CurveRole::Discount,
            "requires a projection curve",
        ),
    ] {
        assert!(matches!(
            MarketStructure::new(as_of, [composite(reference)?, jpy(role)?]),
            Err(Error::InvalidData(message)) if message.contains(diagnostic)
        ));
    }
    let spec = composite(CurveExpr::ForwardRef(id.clone()))?;
    assert!(MarketStructure::new(as_of, [spec, jpy(CurveRole::Projection)?]).is_ok());
    Ok(())
}

#[test]
fn market_checks_theta_length_and_only_the_available_values() -> Result<()> {
    // 利用不可のブロックにも添字を残し、利用できるブロックの配置を一定にする。
    let as_of = "2026-01-01T00:00:00Z".parse()?;
    let jpy = curve("JPY-TONA", Currency::JPY, &["2027-01-01T00:00:00Z"])?;
    let usd = curve("USD-SOFR", Currency::USD, &["2027-01-01T00:00:00Z"])?;
    let jpy_block = BlockId::nodes(jpy.id().clone());
    let usd_block = BlockId::nodes(usd.id().clone());
    let structure = Arc::new(MarketStructure::new(as_of, [jpy, usd])?);
    assert!(matches!(
        Market::from_theta(Arc::clone(&structure), vec![0.0]),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        Market::from_theta(Arc::clone(&structure), vec![0.0, f64::INFINITY]),
        Err(Error::InvalidInput(_))
    ));
    let market = Market::from_partial_theta(
        Arc::clone(&structure),
        vec![0.0, f64::NAN],
        [jpy_block.clone()],
    )?;
    assert!(Arc::ptr_eq(market.structure(), &structure));
    assert_eq!(market.theta().len(), 2);
    assert!(market.block_is_available(&jpy_block));
    assert!(!market.block_is_available(&usd_block));
    let unknown = BlockId::nodes("EUR-ESTR".parse()?);
    assert!(matches!(
        Market::from_partial_theta(structure, vec![0.0, 0.0], [unknown]),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}
`````

`````
crates/rhodium-market/tests/ui/amount_cross_scope.rs
mod common;

use common::fixture;
use rhodium_core::{Result, currency::Currency};

fn main() -> Result<()> {
    let (market, id, _, _) = fixture()?;
    market.with_view(|view| {
        let curve = view.discount_curve(&id)?;
        view.ccy_scope(Currency::JPY, curve, |first| {
            let amount = first.present(5.0);
            view.ccy_scope(Currency::JPY, curve, |second| {
                let _ = amount + second.present(1.0);
                Ok(())
            })
        })
    })
}
`````

`````
crates/rhodium-market/tests/ui/amount_cross_scope_ok.rs
mod common;

use common::fixture;
use rhodium_core::{Result, currency::Currency};

fn main() -> Result<()> {
    let (market, id, _, _) = fixture()?;
    market.with_view(|view| {
        let curve = view.discount_curve(&id)?;
        view.ccy_scope(Currency::JPY, curve, |scope| {
            let _ = scope.present(5.0) + scope.present(1.0);
            Ok(())
        })
    })
}
`````

`````
crates/rhodium-market/tests/ui/amount_plus_undiscounted.rs
mod common;

use common::fixture;
use rhodium_core::{Result, currency::Currency};

fn main() -> Result<()> {
    let (market, id, pay, _) = fixture()?;
    market.with_view(|view| {
        let curve = view.discount_curve(&id)?;
        view.ccy_scope(Currency::JPY, curve, |scope| {
            let at = scope.at(&pay)?;
            let coupon = scope.pv(5.0, at);
            // 将来の未割引の額 100.0 を、そのまま現在価値に足すことはできない。
            let _ = coupon + 100.0;
            Ok(())
        })
    })
}
`````

`````
crates/rhodium-market/tests/ui/amount_plus_undiscounted_ok.rs
mod common;

use common::fixture;
use rhodium_core::{Result, currency::Currency};

fn main() -> Result<()> {
    let (market, id, pay, _) = fixture()?;
    market.with_view(|view| {
        let curve = view.discount_curve(&id)?;
        view.ccy_scope(Currency::JPY, curve, |scope| {
            let at = scope.at(&pay)?;
            let coupon = scope.pv(5.0, at);
            // 将来の額は支払いの点で割り引いてから足す。
            let _ = coupon + scope.pv(100.0, at);
            Ok(())
        })
    })
}
`````

`````
crates/rhodium-market/tests/ui/common/mod.rs
use std::sync::Arc;

use rhodium_core::{Result, currency::Currency, day_count::DayCount, quantity::YearFraction};
use rhodium_market::{
    curve::{CurveId, CurvePoint, CurveRole, CurveSpec, Extrapolation},
    market::Market,
    structure::MarketStructure,
    view::Slots,
};

// 同じ曲線を割引・予想に使い、準備済みの点とスロットを全ての型テストから参照する。
pub fn fixture() -> Result<(Market, CurveId, CurvePoint, Slots)> {
    let id: CurveId = "JPY-TONA".parse()?;
    let start = "2026-01-01T00:00:00Z".parse()?;
    let end = "2027-01-01T00:00:00Z".parse()?;
    let curve = CurveSpec::new(
        id.clone(),
        Currency::JPY,
        [CurveRole::Discount, CurveRole::Projection],
        [end],
        Extrapolation::default(),
    )?;
    let structure = Arc::new(MarketStructure::new(start, [curve])?);
    let pay = structure.prepare_point(&id, end)?;
    let first = structure.prepare_point(&id, start)?;
    let tau = YearFraction::new(1.0, DayCount::Act365F)?;
    let slots = Slots::new(
        Arc::clone(&structure),
        [(id.clone(), vec![pay.clone()])],
        [(id.clone(), vec![(first, pay.clone(), tau)])],
    )?;
    Ok((Market::from_theta(structure, vec![-0.02])?, id, pay, slots))
}
`````

`````
crates/rhodium-market/tests/ui/curve_ref_cross_view.rs
mod common;

use common::fixture;
use rhodium_core::{Result, currency::Currency};

fn main() -> Result<()> {
    let (market, id, _, _) = fixture()?;
    market.with_view(|first| {
        let discount = first.discount_curve(&id)?;
        market.with_view(|second| second.ccy_scope(Currency::JPY, discount, |_| Ok(())))
    })
}
`````

`````
crates/rhodium-market/tests/ui/curve_ref_cross_view_ok.rs
mod common;

use common::fixture;
use rhodium_core::{Result, currency::Currency};

fn main() -> Result<()> {
    let (market, id, _, _) = fixture()?;
    market.with_view(|first| {
        first.ccy_scope(Currency::JPY, first.discount_curve(&id)?, |_| Ok(()))?;
        market.with_view(|second| {
            second.ccy_scope(Currency::JPY, second.discount_curve(&id)?, |_| Ok(()))
        })
    })
}
`````

`````
crates/rhodium-market/tests/ui/discount_point_cross_view.rs
mod common;

use common::fixture;
use rhodium_core::Result;

fn main() -> Result<()> {
    let (market, _, _, slots) = fixture()?;
    market.with_view(|first| {
        let bound = first.bind_slots(&slots)?;
        let at = bound.discounts()[0].points()[0];
        let _ = market.with_view(|second| second.df(at));
        Ok(())
    })
}
`````

`````
crates/rhodium-market/tests/ui/discount_point_cross_view_ok.rs
mod common;

use common::fixture;
use rhodium_core::Result;

fn main() -> Result<()> {
    let (market, _, _, slots) = fixture()?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        let _ = view.df(bound.discounts()[0].points()[0]);
        Ok(())
    })
}
`````

`````
crates/rhodium-market/tests/ui/forward_point_cross_view.rs
mod common;

use common::fixture;
use rhodium_core::Result;

fn main() -> Result<()> {
    let (market, _, _, slots) = fixture()?;
    market.with_view(|first| {
        let bound = first.bind_slots(&slots)?;
        let at = bound.projections()[0].points()[0];
        let _ = market.with_view(|second| second.forward(at));
        Ok(())
    })
}
`````

`````
crates/rhodium-market/tests/ui/forward_point_cross_view_ok.rs
mod common;

use common::fixture;
use rhodium_core::Result;

fn main() -> Result<()> {
    let (market, _, _, slots) = fixture()?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        let _ = view.forward(bound.projections()[0].points()[0]);
        Ok(())
    })
}
`````

`````
crates/rhodium-market/tests/ui/fx_ref_cannot_escape_view.rs
use rhodium_core::{Result, fx::CcyPair, time::Instant};
use rhodium_market::{curve::CurveId, market::Market};

fn forward(market: &Market, pair: CcyPair, collateral: &CurveId, at: Instant) -> Result<f64> {
    let fx = market.with_view(|view| view.fx(pair, collateral))?;
    Ok(fx.forward(at)?.rate())
}

fn main() {}
`````

`````
crates/rhodium-market/tests/ui/fx_ref_cannot_escape_view_ok.rs
use rhodium_core::{Result, fx::CcyPair, time::Instant};
use rhodium_market::{curve::CurveId, market::Market};

fn forward(market: &Market, pair: CcyPair, collateral: &CurveId, at: Instant) -> Result<f64> {
    market.with_view(|view| Ok(view.fx(pair, collateral)?.forward(at)?.rate()))
}

fn main() {}
`````

`````
crates/rhodium-market/tests/ui/growth_point_cross_view.rs
mod common;

use common::fixture;
use rhodium_core::Result;
use rhodium_market::view::Slots;

fn main() -> Result<()> {
    let (market, id, end, _) = fixture()?;
    let structure = market.structure();
    let start = structure.prepare_point(&id, structure.as_of())?;
    let slots = Slots::from_resolved(structure.clone(), [], [])?
        .with_projection_growth([(structure.resolve_projection(&id)?, vec![(start, end)])])?;
    market.with_view(|first| {
        let bound = first.bind_slots(&slots)?;
        let at = bound.projection_growth()[0].points()[0];
        let _ = market.with_view(|second| second.growth(at));
        Ok(())
    })
}
`````

`````
crates/rhodium-market/tests/ui/growth_point_cross_view_ok.rs
mod common;

use common::fixture;
use rhodium_core::Result;
use rhodium_market::view::Slots;

fn main() -> Result<()> {
    let (market, id, end, _) = fixture()?;
    let structure = market.structure();
    let start = structure.prepare_point(&id, structure.as_of())?;
    let slots = Slots::from_resolved(structure.clone(), [], [])?
        .with_projection_growth([(structure.resolve_projection(&id)?, vec![(start, end)])])?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        let value = view.growth(bound.projection_growth()[0].points()[0]);
        assert!((value - 0.02_f64.exp()).abs() <= 1e-12);
        Ok(())
    })
}
`````

`````
crates/rhodium-market/tests/ui/pay_point_cross_scope.rs
mod common;

use common::fixture;
use rhodium_core::{Result, currency::Currency};

fn main() -> Result<()> {
    let (market, id, point, _) = fixture()?;
    market.with_view(|view| {
        let curve = view.discount_curve(&id)?;
        view.ccy_scope(Currency::JPY, curve, |first| {
            let pay = first.at(&point)?;
            view.ccy_scope(Currency::JPY, curve, |second| {
                let _ = second.pv(100.0, pay);
                Ok(())
            })
        })
    })
}
`````

`````
crates/rhodium-market/tests/ui/pay_point_cross_scope_ok.rs
mod common;

use common::fixture;
use rhodium_core::{Result, currency::Currency};

fn main() -> Result<()> {
    let (market, id, point, _) = fixture()?;
    market.with_view(|view| {
        let curve = view.discount_curve(&id)?;
        view.ccy_scope(Currency::JPY, curve, |scope| {
            let pay = scope.at(&point)?;
            let _ = scope.pv(100.0, pay);
            Ok(())
        })
    })
}
`````

`````
crates/rhodium-market/tests/ui/projection_curve_as_discount.rs
mod common;

use common::fixture;
use rhodium_core::{Result, currency::Currency};

fn main() -> Result<()> {
    let (market, id, _, _) = fixture()?;
    market.with_view(|view| {
        let projection = view.projection_curve(&id)?;
        view.ccy_scope(Currency::JPY, projection, |_| Ok(()))
    })
}
`````

`````
crates/rhodium-market/tests/ui/projection_curve_as_discount_ok.rs
mod common;

use common::fixture;
use rhodium_core::{Result, currency::Currency};

fn main() -> Result<()> {
    let (market, id, _, _) = fixture()?;
    market.with_view(|view| {
        let discount = view.discount_curve(&id)?;
        view.ccy_scope(Currency::JPY, discount, |_| Ok(()))
    })
}
`````

`````
crates/rhodium-market/tests/ui/projection_position_as_discount_slot.rs
mod common;

use common::fixture;
use rhodium_core::Result;
use rhodium_market::view::Slots;

fn main() -> Result<()> {
    let (market, id, _, _) = fixture()?;
    let structure = market.structure();
    let projection = structure.resolve_projection(&id)?;
    Slots::from_resolved(structure.clone(), [(projection, Vec::new())], [])?;
    Ok(())
}
`````

`````
crates/rhodium-market/tests/ui/projection_position_as_discount_slot_ok.rs
mod common;

use common::fixture;
use rhodium_core::Result;
use rhodium_market::view::Slots;

fn main() -> Result<()> {
    let (market, id, _, _) = fixture()?;
    let structure = market.structure();
    let discount = structure.resolve_discount(&id)?;
    let projection = structure.resolve_projection(&id)?;
    let slots = Slots::from_resolved(
        structure.clone(),
        [(discount, Vec::new())],
        [(projection, Vec::new())],
    )?;
    market.with_view(|view| {
        view.bind_slots(&slots)?;
        Ok(())
    })
}
`````

`````
crates/rhodium-market/tests/ui/read_with_unbound_point.rs
mod common;

use common::fixture;
use rhodium_core::Result;

fn main() -> Result<()> {
    let (market, _, point, _) = fixture()?;
    market.with_view(|view| {
        let _ = view.df(point);
    });
    Ok(())
}
`````

`````
crates/rhodium-market/tests/ui/read_with_unbound_point_ok.rs
mod common;

use common::fixture;
use rhodium_core::Result;

fn main() -> Result<()> {
    let (market, _, _, slots) = fixture()?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        let _ = view.df(bound.discounts()[0].points()[0]);
        Ok(())
    })
}
`````

`````
crates/rhodium-pricing/Cargo.toml
[package]
name = "rhodium-pricing"
description = "市場の校正、評価の準備と価格計算"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
publish.workspace = true

[lints]
workspace = true

[dependencies]
rhodium-core.workspace = true
rhodium-products.workspace = true
rhodium-market.workspace = true
scioxide-ad.workspace = true
scioxide-core.workspace = true
scioxide-linalg.workspace = true
scioxide-optimize.workspace = true
scioxide-rootfind.workspace = true
scioxide-sparse.workspace = true
serde.workspace = true

[dev-dependencies]
rhodium-testkit.workspace = true
serde_json = "1"
toml = "1.2"
`````

`````
crates/rhodium-pricing/src/calibration.rs
//! OIS のパーレートから曲線を解き、校正の線形化を保存する。

mod groups;
mod solver;

#[cfg(test)]
#[allow(
    dead_code,
    reason = "統合テスト用の QuantLib fixture のうち、再試行には TONA 曲線だけを使う"
)]
#[path = "../tests/fixtures/ois.rs"]
mod ois_fixture;

use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use rhodium_core::{
    Error, Result, Warning,
    calendar::{BusinessDayConvention, CalendarSpec, adjust},
    currency::Currency,
    date::Date,
    fx::CcyPair,
};
use rhodium_market::{
    curve::CurveRole,
    fx::{FxForwardInterval, FxSpec, FxTodayPlan},
    market::Market,
    quote::{QuoteId, QuoteKind, QuoteSet},
    structure::{BlockId, MarketStructure},
};
use scioxide_ad::Real;
use scioxide_core::ResidualTolerance;
use scioxide_linalg::LuDecomposition;
use scioxide_sparse::CsrMatrix;

use crate::{
    context::Context,
    definition::{ElementSelector, InstrumentBinding, MarketDefinition, ResolvedFxDefinition},
    linear::OisWarmStart,
    par,
    prepare::{PrepareInput, Prepared, evaluate, prepare},
};

use self::groups::{GroupForm, GroupPlan};

/// 最小二乗のグループの解の記録の 2 階の項の扱い。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HessianMode {
    /// `∂G/∂θ = J_Fᵀ W J_F + Σ_i (WF)_i ∇²F_i + ∇²P`（既定）。1 次の quote 感応度が厳密になる。
    Exact,
    /// 残差の 2 階の項を省く `J_Fᵀ W J_F + ∇²P`。quote 感応度は近似になる。
    GaussNewton,
}

/// 校正グループの解の線形化 `G(θ_k, q, θ_up) = 0`。
///
/// 正方系では G は残差 F そのもので、行は [`CalibrationGroup::row_quotes`] の順。
/// 最小二乗では G は停留条件 J̃ᵀ r̃（r̃ は重み付けした残差と罰則の行）で、行は
/// [`CalibrationGroup::blocks`] の θ の順。どちらの形でも随伴は同じ式で流れる。
#[derive(Debug)]
pub struct GroupSolution {
    dg_dtheta: LuDecomposition<f64>,
    dg_dq: CsrMatrix<f64>,
    dg_dup: CsrMatrix<f64>,
    residual_norm: f64,
    iterations: usize,
    hessian: Option<HessianMode>,
    residuals: Vec<f64>,
    theta_error: Vec<f64>,
}

impl GroupSolution {
    /// 自グループの θ に対するヤコビアンの LU 分解。
    ///
    /// 列は [`CalibrationGroup::blocks`] のブロック順、各ブロック内は節点順。
    /// 転置系には逆作用素の随伴を使える。
    #[must_use]
    pub const fn dg_dtheta(&self) -> &LuDecomposition<f64> {
        &self.dg_dtheta
    }

    /// 入力 quote に対する残差の微分。列は [`CalibrationGroup::quotes`] の順。
    #[must_use]
    pub const fn dg_dq(&self) -> &CsrMatrix<f64> {
        &self.dg_dq
    }

    /// 上流の θ に対する残差の微分。
    ///
    /// 列は [`CalibrationGroup::upstream`] のブロック順、各ブロック内は節点順。
    #[must_use]
    pub const fn dg_dup(&self) -> &CsrMatrix<f64> {
        &self.dg_dup
    }

    /// 解での残差のユークリッドノルム。
    #[must_use]
    pub const fn residual_norm(&self) -> f64 {
        self.residual_norm
    }

    /// 採用した解を得るまでのソルバーの反復回数。
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }

    /// 最小二乗のグループの解の記録の 2 階の項の扱い。正方系は `None`。
    #[must_use]
    pub const fn hessian(&self) -> Option<HessianMode> {
        self.hessian
    }

    /// 解での残差 `F_i = m_i(θ) − q̃_i`（重みを掛けない）。行は [`CalibrationGroup::row_quotes`] の順。
    #[must_use]
    pub fn residuals(&self) -> &[f64] {
        &self.residuals
    }

    /// 数値の停止の誤差の見積もり δθ = (∂G/∂θ)⁻¹·G。列は [`CalibrationGroup::blocks`] の θ の順。
    ///
    /// ソルバーの誤差を感応度で伝えて比較の幅にするのに使う（品質の「ソルバーの誤差」）。
    #[must_use]
    pub fn theta_error(&self) -> &[f64] {
        &self.theta_error
    }
}

/// グループの校正の結果。
#[derive(Debug)]
pub enum GroupOutcome {
    /// 曲線とその線形化が利用できる。
    Solved(Box<GroupSolution>),
    /// 必要な quote が無いか、校正商品の準備（fixing の欠損、許されない近似など）、
    /// 残差の計算または求根に失敗した。
    Failed(Error),
    /// 必要な上流ブロックが利用できない。
    Unavailable,
}

/// 依存関係と入力を記録した一つの校正グループ。
#[derive(Debug, Clone)]
pub struct CalibrationGroup {
    plan: Arc<GroupPlan>,
    outcome: Arc<GroupOutcome>,
}

impl CalibrationGroup {
    /// このグループが解くブロック。市場の θ 配置順。
    #[must_use]
    pub fn blocks(&self) -> &[BlockId] {
        &self.plan.blocks
    }

    /// 残差が読む上流ブロック。市場の θ 配置順。
    #[must_use]
    pub fn upstream(&self) -> &[BlockId] {
        &self.plan.upstream
    }

    /// 残差が読む quote。ID の昇順で、同じ ID は一度だけ現れる。
    #[must_use]
    pub fn quotes(&self) -> &[QuoteId] {
        &self.plan.quotes
    }

    /// 一つの残差につき一つの quote ID。グループ内の binding の定義順。
    #[must_use]
    pub fn row_quotes(&self) -> &[QuoteId] {
        &self.plan.row_quotes
    }

    /// 校正の結果。
    #[must_use]
    pub fn outcome(&self) -> &GroupOutcome {
        &self.outcome
    }
}

/// 校正済みの市場、その入力 quote、グループごとの結果。
#[derive(Debug, Clone)]
pub struct CalibratedMarket {
    market: Market<f64>,
    quotes: QuoteSet,
    groups: Vec<CalibrationGroup>,
    warnings: Vec<Warning>,
    prepared: Arc<PreparedCalibration>,
}

impl CalibratedMarket {
    /// 利用できるブロックを持つ市場。失敗したグループのブロックは利用できない。
    #[must_use]
    pub const fn market(&self) -> &Market<f64> {
        &self.market
    }

    /// 校正に渡された入力 quote。
    #[must_use]
    pub const fn quotes(&self) -> &QuoteSet {
        &self.quotes
    }

    /// 市場定義が校正入力に割り当てた quote の種類。
    #[must_use]
    pub fn quote_kind(&self, id: &QuoteId) -> Option<QuoteKind> {
        self.prepared.quote_kinds.get(id).copied()
    }

    /// 依存順で並んだグループ。独立のグループは曲線 ID の昇順で選ぶ。
    #[must_use]
    pub fn groups(&self) -> &[CalibrationGroup] {
        &self.groups
    }

    /// 約定日の調整、商品の準備、入力 quote の検査で発生した警告。
    #[must_use]
    pub fn warnings(&self) -> &[Warning] {
        &self.warnings
    }

    /// 同じ評価時点と校正条件で、変更した quote を読むグループとその下流を解き直す。
    ///
    /// 独立したグループの結果と市場構造を共有する。解が一意と分かるグループは前の θ を使う。
    /// 必要な quote の欠損や収束の失敗はグループごとの結果に記録する。
    ///
    /// # Errors
    ///
    /// 市場のパラメータを組み立てられない場合。
    pub fn update(&self, quotes: &QuoteSet) -> Result<Self> {
        let changed = self.quotes.changed_ids(quotes);
        if !changed
            .iter()
            .any(|id| self.prepared.quote_kinds.contains_key(id))
        {
            return Ok(Self {
                quotes: quotes.clone(),
                ..self.clone()
            });
        }
        self.prepared.run(quotes, Some((self, &changed)))
    }
}

/// 同じ評価時点・参照データ・fixing・設定で繰り返し使える校正の準備。
#[derive(Debug)]
struct PreparedCalibration {
    structure: Arc<MarketStructure>,
    rows: Vec<BindingRow>,
    plans: Vec<Arc<GroupPlan>>,
    quote_kinds: BTreeMap<QuoteId, QuoteKind>,
    warnings: BTreeSet<Warning>,
    tolerance: ResidualTolerance<f64>,
}

/// 校正グループの正方系を依存順に解き、利用できるブロックと解の記録を返す。
///
/// グループ単位の失敗は [`GroupOutcome`] に記録する。校正商品の準備の失敗のうち、実行時のデータや
/// 設定によるもの（`MissingData`: 公表済みの fixing や規則の欠損、`Unsupported`: 設定が許さない近似、
/// `Numerical`、`Calibration`）も、その binding が当てはめるグループの失敗にする。
///
/// # Errors
///
/// 市場構造を組み立てられない場合。校正商品の準備が入力どうしの不整合（`InvalidInput`、
/// `InvalidData`: 曲線のピラーを越える校正商品、通貨の食い違いなど）で失敗した場合。
/// 失敗した binding の影響を受けないグループで、正方系として解くのに必要な残差の数がそろわない場合。
/// 校正商品のテンプレートと慣行が適合しないか、観測金利を丸める場合は、quote の ID を付けた `InvalidData`
/// （[`rhodium_products::standard::InstrumentTemplate::validate`]）。
/// 実行の制御（`Execution`）と rhodium の誤り（`Internal`）はそのまま返す。
pub fn calibrate(
    def: &MarketDefinition,
    quotes: &QuoteSet,
    ctx: &Context<'_>,
) -> Result<CalibratedMarket> {
    calibrate_with_residual_tolerance(def, quotes, ctx, 1e-14)
}

/// 各グループの正方系を、指定した絶対許容誤差で校正する。
///
/// 許容誤差はグループの残差のユークリッドノルムに適用する。
/// 残差の各行は、その行の目標の単位で測る（OIS はパーレート、FX ポイントはペアのレートの差、
/// 直接マークは θ の要素）。
/// 中心差分との照合では、通常の [`calibrate`] より厳しい許容誤差を指定できる。
///
/// # Errors
///
/// [`calibrate`] と同じ。許容誤差が有限な非負値でなければ `InvalidInput`。
pub fn calibrate_with_residual_tolerance(
    def: &MarketDefinition,
    quotes: &QuoteSet,
    ctx: &Context<'_>,
    absolute_tolerance: f64,
) -> Result<CalibratedMarket> {
    let tolerance = ResidualTolerance::new(absolute_tolerance, 0.0)
        .map_err(|error| Error::InvalidInput(format!("calibration residual tolerance: {error}")))?;
    let zone = ctx
        .refdata
        .zone(&def.valuation_zone)
        .map_err(|error| error.context("valuation zone"))?;
    let date = zone.to_local(ctx.as_of)?.date();
    let calendar = ctx
        .refdata
        .calendar(&CalendarSpec::Base(def.valuation_calendar.clone()))?;
    let (trade_date, warnings) = adjust(date, BusinessDayConvention::Following, &calendar)?;

    // テンプレートに合う慣行や観測金利の丸めは市場の定義の問題なので、
    // 準備の前に全 binding を確かめ、校正全体のエラーにする。
    for binding in &def.bindings {
        if let Some(template) = binding.target.template() {
            template
                .validate(ctx.refdata)
                .map_err(|error| error.context(format_args!("quote {}", binding.quote)))?;
        }
    }
    let fx = def
        .fx
        .iter()
        .map(|spec| spec.resolve(ctx.as_of, ctx.refdata))
        .collect::<Result<Vec<_>>>()?;
    let structure = build_structure(def, &fx, ctx)?;
    // 校正商品の準備の失敗（公表済みの fixing の欠損、設定が許さない近似など）は、その binding が
    // 当てはめるグループの失敗として残し、独立した曲線（別の通貨など）の校正を続ける。
    let mut rows: Vec<_> = def
        .bindings
        .iter()
        .map(|binding| BindingRow::new(binding, trade_date, &structure, ctx))
        .collect::<Result<_>>()?;
    add_lgm_mark_rows(def, &structure, ctx, &mut rows)?;
    let mut warnings = warnings;
    for pair in fx {
        warnings.extend(pair.warnings);
        if !pair.marks.is_empty() {
            let target = FxSpec::marks_block(pair.spec.pair)?;
            let range = structure
                .layout()
                .range(&target)
                .ok_or_else(|| Error::Internal("missing FX mark block".into()))?;
            rows.extend(range.zip(pair.marks).map(|(index, mark)| {
                BindingRow::Prepared(Box::new(PreparedBinding::direct(
                    mark.quote,
                    target.clone(),
                    index,
                    mark.scale,
                    mark.positive,
                )))
            }));
        }
        if let Some(fit) = pair.fit {
            for row in fit.intervals {
                let points = structure
                    .prepare_fx_interval(pair.spec.pair, &fit.collateral, row.start, row.end)
                    .map_err(|error| {
                        Error::InvalidData(format!("FX quote {}: {error}", row.quote))
                    })?;
                rows.push(BindingRow::Prepared(Box::new(PreparedBinding {
                    quote: row.quote,
                    target: fit.target.clone(),
                    model: PreparedTarget::FxPoints {
                        points: Box::new(points),
                        pip: fit.pip,
                    },
                    warm_start: None,
                    tolerance: None,
                })));
            }
        }
    }
    let warnings = warnings
        .into_iter()
        .chain(rows.iter().flat_map(|row| row.warnings().cloned()))
        .collect();
    let least_squares = groups::LeastSquaresSpec {
        penalties: &def.penalties,
        gauss_newton: &def.gauss_newton,
    };
    let plans = groups::plan(
        structure.layout().blocks(),
        &rows,
        structure.layout(),
        &least_squares,
    )?
    .into_iter()
    .map(Arc::new)
    .collect();
    Arc::new(PreparedCalibration {
        structure,
        rows,
        plans,
        quote_kinds: def.quote_kinds().clone(),
        warnings,
        tolerance,
    })
    .run(quotes, None)
}

/// 市場定義、評価時点で解決した FX の定義から、予想の測度を記録した市場構造を作る。
fn build_structure(
    def: &MarketDefinition,
    fx: &[ResolvedFxDefinition],
    ctx: &Context<'_>,
) -> Result<Arc<MarketStructure>> {
    let curves = def.curves.iter().cloned().chain(
        fx.iter()
            .filter_map(|pair| pair.fit.as_ref().map(|fit| fit.curve.clone())),
    );
    let structure = MarketStructure::with_fx_and_lgm(
        ctx.as_of,
        curves,
        fx.iter().map(|pair| pair.spec.clone()),
        def.lgm
            .iter()
            .map(|(&ccy, model)| (ccy, model.breaks.clone())),
    )?;
    // 予想と割引の担保の違いは、予想曲線を当てる校正商品の担保で判定する。校正商品の準備も
    // 同じ判定を使うので、binding の準備より前に記録する。割引の規則が構造に無い曲線や、
    // 通貨の違う曲線を指す測度は、その曲線を読む取引の準備が報告するので、ここでは記録しない。
    let measures: Vec<_> = def
        .forecast_measures(ctx.config)
        .into_iter()
        .filter(|(projection, discount)| {
            structure.curve(discount).is_some_and(|spec| {
                spec.has_role(CurveRole::Discount)
                    && structure
                        .curve(projection)
                        .is_some_and(|projection| projection.currency() == spec.currency())
            })
        })
        .collect();
    Ok(Arc::new(structure.with_forecast_measures(measures)?))
}

fn add_lgm_mark_rows(
    def: &MarketDefinition,
    structure: &MarketStructure,
    ctx: &Context<'_>,
    rows: &mut Vec<BindingRow>,
) -> Result<()> {
    for (&ccy, model) in &def.lgm {
        ctx.refdata.currency(ccy)?;
        let target = MarketStructure::lgm_block(ccy)?;
        let range = structure
            .layout()
            .range(&target)
            .ok_or_else(|| Error::Internal("missing LGM marks block".into()))?;
        rows.extend(
            range
                .zip(model.mean_reversion.iter().chain(&model.volatility))
                .map(|(index, quote)| {
                    BindingRow::Prepared(Box::new(PreparedBinding::direct(
                        quote.clone(),
                        target.clone(),
                        index,
                        1.0,
                        false,
                    )))
                }),
        );
    }
    Ok(())
}

impl PreparedCalibration {
    fn run(
        self: &Arc<Self>,
        quotes: &QuoteSet,
        previous: Option<(&CalibratedMarket, &BTreeSet<QuoteId>)>,
    ) -> Result<CalibratedMarket> {
        // 利用できる θ だけを書き込む。解けなかったグループの古い値を下流が読まないようにする。
        let mut theta = vec![0.0; self.structure.layout().len()];
        let mut available = BTreeSet::new();
        let mut changed_blocks = BTreeSet::new();
        let mut groups = Vec::with_capacity(self.plans.len());
        for (index, plan) in self.plans.iter().enumerate() {
            let outcome = match previous {
                Some((old, changed))
                    if !plan.quotes.iter().any(|id| changed.contains(id))
                        && !plan.upstream.iter().any(|id| changed_blocks.contains(id)) =>
                {
                    let outcome = &old.groups[index].outcome;
                    if matches!(outcome.as_ref(), GroupOutcome::Solved(_)) {
                        for &index in &plan.parameters {
                            theta[index] = old.market.theta()[index];
                        }
                    }
                    Arc::clone(outcome)
                }
                _ => {
                    // 成否にかかわらず下流を更新する。欠損の補完で利用可能になった場合も同じ順序を使う。
                    changed_blocks.extend(plan.blocks.iter());
                    let previous_theta = previous.and_then(|(old, _)| {
                        matches!(old.groups[index].outcome(), GroupOutcome::Solved(_))
                            .then_some(old.market.theta())
                    });
                    Arc::new(self.solve_group(plan, &mut theta, &available, quotes, previous_theta))
                }
            };
            if matches!(outcome.as_ref(), GroupOutcome::Solved(_)) {
                available.extend(plan.blocks.iter().cloned());
            }
            groups.push(CalibrationGroup {
                plan: Arc::clone(plan),
                outcome,
            });
        }
        let market = Market::from_partial_theta(Arc::clone(&self.structure), theta, available)?;
        // 最小二乗の当てはめの警告は解に由来するので、再利用したグループを含めて毎回求める。
        let fit_warnings: Vec<_> = groups.iter().flat_map(least_squares_warnings).collect();
        let warnings = self
            .warnings
            .iter()
            .cloned()
            .chain(unit_warnings(&self.quote_kinds, quotes))
            .chain(market.fx_cross_warnings())
            .chain(fit_warnings)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        Ok(CalibratedMarket {
            market,
            quotes: quotes.clone(),
            groups,
            warnings,
            prepared: Arc::clone(self),
        })
    }

    fn solve_group(
        &self,
        plan: &GroupPlan,
        theta: &mut [f64],
        available: &BTreeSet<BlockId>,
        quotes: &QuoteSet,
        previous_theta: Option<&[f64]>,
    ) -> GroupOutcome {
        if let Some(error) = plan.rows.iter().find_map(|&row| self.rows[row].failure()) {
            return GroupOutcome::Failed(error.clone());
        }
        if !plan.upstream.iter().all(|block| available.contains(block)) {
            return GroupOutcome::Unavailable;
        }
        match solver::solve(
            plan,
            &self.rows,
            &self.structure,
            theta,
            quotes,
            self.tolerance,
            previous_theta,
        ) {
            Ok((point, solution)) => {
                for (&index, value) in plan.parameters.iter().zip(point) {
                    theta[index] = value;
                }
                GroupOutcome::Solved(Box::new(solution))
            }
            Err(error) => GroupOutcome::Failed(error),
        }
    }
}

/// 解いた最小二乗のグループで、残差が許容幅 σ̃ を超えた quote と、Gauss–Newton で記録した quote。
fn least_squares_warnings(group: &CalibrationGroup) -> Vec<Warning> {
    let (GroupForm::LeastSquares(plan), GroupOutcome::Solved(solution)) =
        (&group.plan.form, group.outcome.as_ref())
    else {
        return Vec::new();
    };
    let outside = group
        .plan
        .row_quotes
        .iter()
        .zip(&plan.weights)
        .zip(solution.residuals())
        // 重みは 1/σ̃ なので、|F|·w > 1 が |F| > σ̃ と同じ。
        .filter(|&((_, weight), residual)| (residual * weight).abs() > 1.0)
        .map(|((quote, _), _)| Warning::CalibrationFitOutsideTolerance {
            quote: quote.clone(),
        });
    let approximate = group
        .plan
        .row_quotes
        .iter()
        .filter(|_| plan.hessian == HessianMode::GaussNewton)
        .map(|quote| Warning::CalibrationGaussNewton {
            quote: quote.clone(),
        });
    outside.chain(approximate).collect()
}

/// 金利かスプレッドとして読む quote のうち、絶対値が 1 を超えるものを警告する。
///
/// quote の種類は binding が決めるので、`QuoteSet` ではなく校正の入口で検査する。
/// 表記（`5`、`"5%"`）は読み込みで失われるので、正規化した値の絶対値で判定する。
/// 100% を超える金利やスプレッドは、単位を付けた入力でも疑わしい。
fn unit_warnings<'a>(
    kinds: &'a BTreeMap<QuoteId, QuoteKind>,
    quotes: &'a QuoteSet,
) -> impl Iterator<Item = Warning> + 'a {
    kinds
        .iter()
        .filter(|&(_, kind)| {
            matches!(
                kind,
                QuoteKind::Rate | QuoteKind::Spread | QuoteKind::CreditSpread
            )
        })
        .filter(|&(id, _)| quotes.get(id).is_some_and(|value| value.abs() > 1.0))
        .map(|(id, _)| Warning::QuoteUnitSuspected { quote: id.clone() })
}

/// 一つの binding の校正商品の準備の結果。
///
/// 実行時のデータや設定のために準備できなかった binding は、その当てはめ先のグループを
/// `GroupOutcome::Failed` にする（入力どうしの不整合は校正全体のエラー）。
/// 準備できなかった商品の足跡（読む曲線）は分からないので、依存グラフの辺を持たない。
#[derive(Debug)]
enum BindingRow {
    /// 二つの準備を持つので失敗の場合より大きく、箱に入れる（binding の数だけの小さな割り当て）。
    Prepared(Box<PreparedBinding>),
    Failed {
        quote: QuoteId,
        target: BlockId,
        error: Error,
    },
}

impl BindingRow {
    fn new(
        binding: &InstrumentBinding,
        trade_date: Date,
        structure: &Arc<MarketStructure>,
        ctx: &Context<'_>,
    ) -> Result<Self> {
        match PreparedBinding::new(binding, trade_date, structure, ctx) {
            Ok(prepared) => Ok(Self::Prepared(Box::new(prepared))),
            // 実行時のデータや設定のために、この校正商品だけを準備できない（公表済みの fixing や
            // 割引・予想の規則の欠損、設定が許さない近似、数値の失敗）。他の曲線の校正は続ける。
            Err(
                error @ (Error::MissingData(_)
                | Error::Unsupported(_)
                | Error::Numerical(_)
                | Error::Calibration(_)),
            ) => Ok(Self::Failed {
                quote: binding.quote.clone(),
                target: binding.determines.block().clone(),
                error,
            }),
            // 入力どうしの不整合（曲線のピラーを越える校正商品、通貨の食い違い、不正な商品の条件）は
            // 市場の定義の誤りなので、従来どおり校正全体のエラーにする。実行の制御と rhodium の誤りも
            // 呼び出し全体を止める。
            Err(
                error @ (Error::InvalidInput(_)
                | Error::InvalidData(_)
                | Error::Dsl { .. }
                | Error::Execution(_)
                | Error::Internal(_)),
            ) => Err(error),
        }
    }

    /// 比べる値を読む quote。
    fn quote(&self) -> &QuoteId {
        match self {
            Self::Prepared(prepared) => &prepared.quote,
            Self::Failed { quote, .. } => quote,
        }
    }

    /// 残差が当てはめるブロック。
    fn target(&self) -> &BlockId {
        match self {
            Self::Prepared(prepared) => &prepared.target,
            Self::Failed { target, .. } => target,
        }
    }

    /// 残差が読む θ ブロック。準備に失敗した binding は空。
    fn footprint(&self) -> impl Iterator<Item = &BlockId> {
        self.prepared()
            .into_iter()
            .flat_map(PreparedBinding::footprint)
    }

    /// 準備した校正商品。
    const fn prepared(&self) -> Option<&PreparedBinding> {
        match self {
            Self::Prepared(prepared) => Some(prepared),
            Self::Failed { .. } => None,
        }
    }

    /// 準備の失敗。
    const fn failure(&self) -> Option<&Error> {
        match self {
            Self::Prepared(_) => None,
            Self::Failed { error, .. } => Some(error),
        }
    }

    /// 準備で出た警告。
    fn warnings(&self) -> impl Iterator<Item = &Warning> {
        self.prepared()
            .into_iter()
            .flat_map(PreparedBinding::warnings)
    }
}

/// 残差を同じ市場で繰り返し評価するための準備。
#[derive(Debug)]
struct PreparedBinding {
    quote: QuoteId,
    target: BlockId,
    model: PreparedTarget,
    warm_start: Option<OisWarmStart>,
    /// 最小二乗の許容幅 σ（quote の単位）。正方系の行は `None`。
    tolerance: Option<f64>,
}

#[derive(Debug)]
enum PreparedTarget {
    /// affine な OIS のパーを PV(0) と PV(1) から得る。
    ParRate {
        ccy: Currency,
        zero: Box<Prepared>,
        one: Box<Prepared>,
        price_quote: bool,
        fx: Option<Box<ParFx>>,
    },
    FxOutright {
        ccy: Currency,
        unit_strike: Box<Prepared>,
        fx: Box<ParFx>,
    },
    Direct {
        index: usize,
        scale: f64,
        positive: bool,
    },
    FxPoints {
        points: Box<FxForwardInterval>,
        pip: f64,
    },
}

#[derive(Debug)]
struct ParFx {
    other: Currency,
    today: FxTodayPlan,
    blocks: Vec<BlockId>,
}

impl PreparedBinding {
    fn new(
        binding: &InstrumentBinding,
        trade_date: Date,
        structure: &Arc<MarketStructure>,
        ctx: &Context<'_>,
    ) -> Result<Self> {
        Self::prepare(binding, trade_date, structure, ctx)
            .map(|prepared| Self {
                tolerance: binding.tolerance(),
                ..prepared
            })
            .map_err(|error| error.context(format_args!("quote {}", binding.quote)))
    }

    /// 直接マーク（θ − D·q）の行か。全行が直接マークの最小二乗は θ について強凸な二次式。
    const fn is_direct(&self) -> bool {
        matches!(self.model, PreparedTarget::Direct { .. })
    }

    fn prepare(
        binding: &InstrumentBinding,
        trade_date: Date,
        structure: &Arc<MarketStructure>,
        ctx: &Context<'_>,
    ) -> Result<Self> {
        let quote = binding.quote.clone();
        let target = binding.determines.block().clone();
        let Some(template) = binding.target.template() else {
            let ElementSelector::Param(param) = &binding.determines else {
                return Err(Error::Internal("direct mark has no parameter".into()));
            };
            let index = structure.layout().index(param).ok_or_else(|| {
                Error::InvalidData(format!("unknown direct mark parameter {param:?}"))
            })?;
            // 直接マークの残差は θ − D·q。D は quote の変換が決める（ターンなら τ_turn）。
            return Ok(Self::direct(
                quote,
                target,
                index,
                binding.transform().scale(),
                binding.kind() == QuoteKind::FxSpot,
            ));
        };
        let at_rate = |rate| {
            let product = template.build_for_valuation(rate, trade_date, ctx.as_of, ctx.refdata)?;
            prepare(
                PrepareInput::Standard {
                    product: &product,
                    trade_date,
                },
                &binding.collateral,
                &[],
                structure,
                ctx,
            )
            .map(Box::new)
        };
        let ccy = template.currency(ctx.refdata)?;
        let fx = template
            .other_currency(ctx.refdata)?
            .map(|other| {
                let today = structure.fx_today_plan(CcyPair::new(other, ccy)?)?;
                let blocks = today.blocks().into_iter().collect();
                Ok::<_, Error>(Box::new(ParFx {
                    other,
                    today,
                    blocks,
                }))
            })
            .transpose()?;
        if matches!(
            template,
            rhodium_products::standard::InstrumentTemplate::FxForward { .. }
        ) {
            let fx =
                fx.ok_or_else(|| Error::Internal("FX outright missing base currency".into()))?;
            return Ok(Self {
                quote,
                target,
                warm_start: None,
                tolerance: None,
                model: PreparedTarget::FxOutright {
                    ccy,
                    unit_strike: at_rate(1.0)?,
                    fx,
                },
            });
        }
        let zero = at_rate(0.0)?;
        let one = at_rate(1.0)?;
        let warm_start = if matches!(
            template,
            rhodium_products::standard::InstrumentTemplate::Ois { .. }
        ) {
            one.ois_warm_start(&target)
        } else {
            None
        };
        Ok(Self {
            quote,
            target,
            warm_start,
            tolerance: None,
            model: PreparedTarget::ParRate {
                ccy,
                zero,
                one,
                price_quote: matches!(
                    binding.target,
                    crate::definition::CalibrationTarget::Price(_)
                ),
                fx,
            },
        })
    }

    fn direct(quote: QuoteId, target: BlockId, index: usize, scale: f64, positive: bool) -> Self {
        Self {
            quote,
            target,
            warm_start: None,
            tolerance: None,
            model: PreparedTarget::Direct {
                index,
                scale,
                positive,
            },
        }
    }

    /// 残差が読む θ ブロック。直接指定は自分の要素だけを読む。
    fn footprint(&self) -> impl Iterator<Item = &BlockId> {
        let parts: [&[BlockId]; 3] = match &self.model {
            PreparedTarget::ParRate { zero, one, fx, .. } => [
                zero.blocks(),
                one.blocks(),
                fx.as_ref().map_or(&[], |plan| plan.blocks.as_slice()),
            ],
            PreparedTarget::FxOutright {
                unit_strike, fx, ..
            } => [unit_strike.blocks(), &fx.blocks, &[]],
            PreparedTarget::Direct { .. } => [&[], &[], &[]],
            PreparedTarget::FxPoints { points, .. } => [points.blocks(), &[], &[]],
        };
        parts.into_iter().flatten()
    }

    fn warnings(&self) -> impl Iterator<Item = &Warning> {
        let parts: [&[Warning]; 2] = match &self.model {
            PreparedTarget::ParRate { zero, one, .. } => [zero.warnings(), one.warnings()],
            PreparedTarget::FxOutright { unit_strike, .. } => [unit_strike.warnings(), &[]],
            PreparedTarget::Direct { .. } | PreparedTarget::FxPoints { .. } => [&[], &[]],
        };
        parts.into_iter().flatten()
    }

    fn direct_value(&self, quote: f64) -> Result<Option<(usize, f64)>> {
        let PreparedTarget::Direct {
            index,
            scale,
            positive,
        } = self.model
        else {
            return Ok(None);
        };
        let value = quote * scale;
        if !value.is_finite() || (positive && value <= 0.0) {
            return Err(Error::InvalidInput(format!(
                "invalid direct mark {}",
                self.quote
            )));
        }
        Ok(Some((index, value)))
    }

    fn quote_scale(&self) -> f64 {
        match self.model {
            PreparedTarget::Direct { scale, .. } => scale,
            PreparedTarget::ParRate { .. } | PreparedTarget::FxOutright { .. } => 1.0,
            PreparedTarget::FxPoints { pip, .. } => pip,
        }
    }

    fn quote_derivative(&self, quote: f64) -> f64 {
        if matches!(self.model, PreparedTarget::FxOutright { .. }) {
            quote.recip()
        } else {
            self.quote_scale()
        }
    }

    fn seed_rate(&self, quote: f64) -> f64 {
        match self.model {
            PreparedTarget::FxPoints { .. }
            | PreparedTarget::FxOutright { .. }
            | PreparedTarget::ParRate {
                price_quote: true, ..
            } => 0.0,
            _ => quote,
        }
    }

    /// 残差 m(θ) − q̃(q)。市場から読む値に数値型を通す。
    fn residual<R: Real>(&self, market: &Market<R>, quote: f64) -> Result<R> {
        let model = match &self.model {
            PreparedTarget::ParRate {
                ccy, zero, one, fx, ..
            } => {
                let zero = evaluate(zero, market)?;
                let one = evaluate(one, market)?;
                let (zero, one) = if let Some(fx) = fx {
                    let rate = fx.today.rate(market)?.rate();
                    (
                        zero.get(*ccy) + zero.get(fx.other) * rate,
                        one.get(*ccy) + one.get(fx.other) * rate,
                    )
                } else {
                    (zero.get(*ccy), one.get(*ccy))
                };
                par::affine(zero, one)?
            }
            PreparedTarget::FxOutright {
                ccy,
                unit_strike,
                fx,
            } => {
                if quote <= 0.0 {
                    return Err(Error::InvalidData(format!(
                        "nonpositive FX outright for {}",
                        self.quote
                    )));
                }
                // 両通貨の受渡を価格評価で割り引いてから対数で照合する。
                // 大きなアウトライトでも、二つの合算 PV の差で傾きや精度を失わない。
                let pv = evaluate(unit_strike, market)?;
                let log_model = pv.get(fx.other).ln() + fx.today.rate(market)?.rate().ln()
                    - (-pv.get(*ccy)).ln();
                return Ok(log_model - R::constant(quote.ln()));
            }
            PreparedTarget::Direct { index, .. } => market.theta()[*index],
            PreparedTarget::FxPoints { points, .. } => points.points(market)?,
        };
        Ok(model - R::constant(self.quote_scale() * quote))
    }
}

#[cfg(test)]
mod tests {
    use rhodium_market::quote::{QuoteOverride, QuoteOverrides};
    use rhodium_testkit::assert_close;

    use super::{
        ois_fixture::{Inputs, cases, definition, fixture_quotes},
        *,
    };

    #[test]
    fn quote_update_retries_with_a_flat_seed_when_the_previous_theta_fails() -> Result<()> {
        // 前回の θ を、残差の最初の評価で DF を表現できない値に置き換える。
        // 指定の OIS は一意で、quote の中央値を用いた初期値からは普通に解ける。
        let [tona, _] = cases()?;
        let collateral = crate::resolve::CollateralKey::Uncollateralized;
        let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
        let def = definition(&[(&tona, &collateral)])?;
        let quotes = QuoteSet::new(fixture_quotes(&tona)?)?;
        let ctx = inputs.context();
        let mut base = calibrate(&def, &quotes, &ctx)?;
        base.market = Market::from_theta(
            Arc::clone(base.market().structure()),
            vec![-2000.0; base.market().theta().len()],
        )?;
        assert!(
            base.prepared
                .rows
                .iter()
                .filter_map(BindingRow::prepared)
                .any(|row| { row.residual(base.market(), 0.0).is_err() })
        );
        let overrides = QuoteOverrides::new(
            [QuoteOverride::Set {
                id: tona.quote(&tona.nodes[0])?,
                value: tona.nodes[0].rate + 0.0001,
            }],
            def.quote_kinds(),
        )?;
        let (bumped, _) = quotes.with_overrides(&overrides)?;
        let updated = base.update(&bumped)?;
        let full = calibrate(&def, &bumped, &ctx)?;
        assert!(matches!(
            updated.groups()[0].outcome(),
            GroupOutcome::Solved(_)
        ));
        for (&actual, &expected) in updated.market().theta().iter().zip(full.market().theta()) {
            assert_close!(actual, expected, 1e-12, 1e-14);
        }
        Ok(())
    }
}
`````

`````
crates/rhodium-pricing/src/calibration/groups.rs
//! 準備した評価の足跡から、曲線ブロックの依存順を求める。

use std::collections::{BTreeMap, BTreeSet};

use rhodium_core::{Error, Result};
use rhodium_market::{
    quote::QuoteId,
    structure::{BlockId, ParamLayout},
};
use scioxide_sparse::{CsrMatrix, GraphDirection, GraphPlan};

use super::BindingRow;
use crate::definition::StepPenalty;

/// グループの方程式の形。
#[derive(Debug)]
pub(super) enum GroupForm {
    /// 残差の数が未知数と等しい正方系 F = 0。
    Square,
    /// 許容幅で重み付けした残差と罰則の二乗和を最小にする（停留条件 G = 0）。
    LeastSquares(LeastSquaresPlan),
}

/// 最小二乗のグループの重み、罰則の行、解の記録の作り方。
#[derive(Debug)]
pub(super) struct LeastSquaresPlan {
    /// 残差の行ごとの重み 1/σ̃。σ̃ は残差の単位の許容幅。
    pub(super) weights: Vec<f64>,
    /// 罰則の行 √λ·(θ[upper] − θ[lower])/s。添字は自グループの未知数の中の位置。
    pub(super) penalties: Vec<PenaltyRow>,
    /// 解の記録の 2 階の項の扱い。
    pub(super) hessian: super::HessianMode,
    /// すべての行が直接マークで、目的が θ について強凸な二次式か（解が一意で warm start できる）。
    pub(super) convex: bool,
}

/// 隣り合う会合の区間のフォワードの差の罰則の一行。
#[derive(Debug, Clone, Copy)]
pub(super) struct PenaltyRow {
    pub(super) lower: usize,
    pub(super) upper: usize,
    pub(super) factor: f64,
}

/// 一つの強連結成分を解くための、ブロック・θ の添字・残差の対応。
///
/// 各添字列の順は、そのまま [`super::GroupSolution`] の行列の行と列の順になる。
#[derive(Debug)]
pub(super) struct GroupPlan {
    /// このグループが解くブロック。曲線 ID の昇順で、θ の配置順と一致する。
    pub(super) blocks: Vec<BlockId>,
    /// 自グループの残差が読む、他の成分のブロック。曲線 ID の昇順。
    pub(super) upstream: Vec<BlockId>,
    /// `blocks` の θ の添字。ブロック順・節点順で、`dg_dtheta` の列の順。
    pub(super) parameters: Vec<usize>,
    /// `upstream` の θ の添字。ブロック順・節点順で、`dg_dup` の列の順。
    pub(super) upstream_parameters: Vec<usize>,
    /// 当てはめ先がこのグループに属する binding の添字。定義順で、残差の行の順。
    pub(super) rows: Vec<usize>,
    /// 残差の行ごとの quote。binding の定義順。
    pub(super) row_quotes: Vec<QuoteId>,
    /// `rows` が読む quote。ID の昇順で重複がなく、`dg_dq` の列の順。
    pub(super) quotes: Vec<QuoteId>,
    /// 方程式の形。
    pub(super) form: GroupForm,
}

/// 最小二乗の指定（罰則と Gauss–Newton のブロック）。
pub(super) struct LeastSquaresSpec<'a> {
    pub(super) penalties: &'a BTreeMap<BlockId, StepPenalty>,
    pub(super) gauss_newton: &'a BTreeSet<BlockId>,
}

/// 校正商品の足跡から依存グラフを作り、上流から解く順にグループを並べる。
///
/// 辺は、残差が読むブロックから当てはめ先のブロックへ張る（上流から下流）。
/// 強連結成分を一つのグループとする。許容幅を持たない binding の成分は正方系で、
/// 未知数と残差の数が一致しなければ拒否する。すべての binding が許容幅を持つ成分は最小二乗で、
/// 残差と罰則の行の数が未知数以上でなければ拒否する。
///
/// 準備に失敗した binding は足跡が分からないので辺を持たない。その binding を含む成分は
/// 校正で `Failed` になり、下流（推移的）は `Unavailable` になるので、これらの成分では
/// 未知数と残差の数を比べない。足跡が欠けて本来一つの成分が分かれても、分かれた側は
/// 失敗した成分の下流になり、解かれない。
pub(super) fn plan(
    blocks: &[BlockId],
    rows: &[BindingRow],
    layout: &ParamLayout,
    least_squares: &LeastSquaresSpec<'_>,
) -> Result<Vec<GroupPlan>> {
    // 頂点が無いと依存グラフの行列を作れないので、曲線が無い定義は先に返す。
    if blocks.is_empty() {
        return Ok(Vec::new());
    }
    // 市場構造は定義の曲線だけから作るので、当てはめ先も足跡もこの表に必ずある。
    let block_positions: BTreeMap<_, _> = blocks.iter().zip(0..).collect();
    let mut edges = BTreeSet::new();
    for row in rows {
        let target = block_positions[row.target()];
        for dependency in row.footprint() {
            let source = block_positions[dependency];
            if source != target {
                edges.insert((source, target));
            }
        }
    }
    let (members, labels) = ordered_components(blocks.len(), &edges)?;
    // 残差は当てはめ先のブロックが属する成分に入る。binding の定義順を保つ。
    let mut rows_by_component = vec![Vec::new(); members.len()];
    for (index, row) in rows.iter().enumerate() {
        rows_by_component[labels[block_positions[row.target()]]].push(index);
    }
    // 準備に失敗した binding を含むか、その下流にある成分。成分は上流から並ぶので、順に伝える。
    let mut tainted: Vec<_> = rows_by_component
        .iter()
        .map(|component| component.iter().any(|&row| rows[row].failure().is_some()))
        .collect();
    for component in 0..members.len() {
        if edges.iter().any(|&(source, target)| {
            labels[target] == component && labels[source] != component && tainted[labels[source]]
        }) {
            tainted[component] = true;
        }
    }
    // ブロックの列を、θ の配置での添字の列に展開する。
    let parameters = |blocks: &[BlockId]| {
        blocks.iter().try_fold(Vec::new(), |mut indices, block| {
            let range = layout.range(block).ok_or_else(|| {
                Error::Internal(format!(
                    "calibration block {} has no parameter range",
                    block.curve()
                ))
            })?;
            indices.extend(range);
            Ok(indices)
        })
    };
    members
        .into_iter()
        .enumerate()
        .map(|(component, members)| {
            // 足跡のうち他の成分に属するブロックが上流。トポロジカル順により先に解かれている。
            let upstream: Vec<_> = rows_by_component[component]
                .iter()
                .flat_map(|&row| rows[row].footprint())
                .filter(|&block| labels[block_positions[block]] != component)
                .cloned()
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            let blocks: Vec<_> = members
                .into_iter()
                .map(|index| blocks[index].clone())
                .collect();
            let indices = parameters(&blocks)?;
            let row_indices = std::mem::take(&mut rows_by_component[component]);
            // 失敗した binding の影響を受ける成分は解かないので、形の検査をしない。
            let form = if tainted[component] {
                GroupForm::Square
            } else {
                group_form(&blocks, &indices, &row_indices, rows, layout, least_squares)?
            };
            let row_quotes: Vec<_> = row_indices
                .iter()
                .map(|&index| rows[index].quote().clone())
                .collect();
            let quotes = row_quotes
                .iter()
                .cloned()
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            Ok(GroupPlan {
                upstream_parameters: parameters(&upstream)?,
                parameters: indices,
                rows: row_indices,
                row_quotes,
                quotes,
                blocks,
                upstream,
                form,
            })
        })
        .collect()
}

/// 成分の binding の許容幅と罰則から方程式の形を決め、行の数を検査する。
fn group_form(
    blocks: &[BlockId],
    parameters: &[usize],
    row_indices: &[usize],
    rows: &[BindingRow],
    layout: &ParamLayout,
    spec: &LeastSquaresSpec<'_>,
) -> Result<GroupForm> {
    let names = || {
        blocks
            .iter()
            .map(|block| block.curve().as_str())
            .collect::<Vec<_>>()
            .join(", ")
    };
    let prepared = row_indices
        .iter()
        .map(|&index| {
            rows[index].prepared().ok_or_else(|| {
                Error::Internal("untainted calibration group has a failed binding".into())
            })
        })
        .collect::<Result<Vec<_>>>()?;
    // 罰則の行は、ブロックの隣り合う区間ごとに一つ。添字は自グループの未知数の中の位置。
    let penalties = blocks
        .iter()
        .filter_map(|block| spec.penalties.get(block).map(|penalty| (block, penalty)))
        .map(|(block, penalty)| {
            let range = layout
                .range(block)
                .ok_or_else(|| Error::Internal("penalized block has no range".into()))?;
            let offset = parameters
                .iter()
                .position(|&index| index == range.start)
                .ok_or_else(|| Error::Internal("penalized block outside its group".into()))?;
            Ok((0..range.len() - 1).map(move |element| PenaltyRow {
                lower: offset + element,
                upper: offset + element + 1,
                factor: penalty.factor(),
            }))
        })
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    let gauss_newton = blocks.iter().any(|block| spec.gauss_newton.contains(block));
    let tolerances: Vec<_> = prepared.iter().map(|row| row.tolerance).collect();
    match (
        tolerances.iter().all(Option::is_some),
        tolerances.iter().all(Option::is_none),
    ) {
        // 許容幅を持たない binding だけの成分は正方系。未知数と残差の数が一致しなければならない。
        (_, true) => {
            if !penalties.is_empty() || gauss_newton {
                return Err(Error::InvalidInput(format!(
                    "calibration group [{}] has a step penalty or Gauss-Newton option but no binding tolerances",
                    names()
                )));
            }
            if prepared.len() != parameters.len() {
                return Err(Error::InvalidInput(format!(
                    "calibration group [{}] has {} parameters and {} residuals",
                    names(),
                    parameters.len(),
                    prepared.len()
                )));
            }
            Ok(GroupForm::Square)
        }
        (true, false) => {
            if prepared.len() + penalties.len() < parameters.len() {
                return Err(Error::InvalidInput(format!(
                    "least-squares calibration group [{}] has {} parameters but {} residual and penalty rows",
                    names(),
                    parameters.len(),
                    prepared.len() + penalties.len()
                )));
            }
            // 残差の単位の許容幅 σ̃ = σ·D（D は quote の変換の係数）の逆数を重みにする。
            let weights = prepared
                .iter()
                .zip(&tolerances)
                .map(|(row, tolerance)| {
                    tolerance
                        .map(|sigma| (sigma * row.quote_scale()).recip())
                        .ok_or_else(|| Error::Internal("missing binding tolerance".into()))
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(GroupForm::LeastSquares(LeastSquaresPlan {
                weights,
                penalties,
                hessian: if gauss_newton {
                    super::HessianMode::GaussNewton
                } else {
                    super::HessianMode::Exact
                },
                convex: prepared.iter().all(|row| row.is_direct()),
            }))
        }
        (false, false) => Err(Error::InvalidInput(format!(
            "calibration group [{}] mixes bindings with and without tolerances: {}",
            names(),
            prepared
                .iter()
                .filter(|row| row.tolerance.is_none())
                .map(|row| row.quote.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        ))),
    }
}

/// `len` 個のブロックと有向辺から強連結成分を求め、上流から解く順に並べる。
///
/// 返り値の第一要素は、並べた順の各成分に属するブロックの添字（昇順）。
/// 第二要素は、各ブロックが属する成分の並べた後の番号。
fn ordered_components(
    len: usize,
    edges: &BTreeSet<(usize, usize)>,
) -> Result<(Vec<Vec<usize>>, Vec<usize>)> {
    let (sources, targets): (Vec<_>, Vec<_>) = edges.iter().copied().unzip();
    let graph = CsrMatrix::from_triplets(len, len, &sources, &targets, &vec![1.0; edges.len()])
        .and_then(|matrix| GraphPlan::new(matrix.view(), GraphDirection::Directed))
        .and_then(|plan| plan.strong_components())
        .map_err(|error| Error::Numerical(format!("calibration dependency graph: {error}")))?;
    let labels = graph.labels();
    let mut members = vec![Vec::new(); graph.count()];
    for (block, &label) in labels.iter().enumerate() {
        members[label].push(block);
    }

    // 成分を一つの頂点に縮めた DAG を走査する。独立な成分は先頭ブロックの順で選ぶ。
    let mut downstream = vec![BTreeSet::new(); graph.count()];
    let mut incoming = vec![0; graph.count()];
    for &(source, target) in edges {
        let (source, target) = (labels[source], labels[target]);
        if source != target && downstream[source].insert(target) {
            incoming[target] += 1;
        }
    }
    let mut ready: BTreeSet<_> = incoming
        .iter()
        .enumerate()
        .filter(|&(_, &count)| count == 0)
        .map(|(label, _)| (members[label][0], label))
        .collect();
    let mut result = Vec::with_capacity(graph.count());
    let mut ordered_labels = vec![0; len];
    while let Some((_, label)) = ready.pop_first() {
        for &block in &members[label] {
            ordered_labels[block] = result.len();
        }
        result.push(std::mem::take(&mut members[label]));
        for &target in &downstream[label] {
            incoming[target] -= 1;
            if incoming[target] == 0 {
                ready.insert((members[target][0], target));
            }
        }
    }
    if result.len() != graph.count() {
        return Err(Error::Internal("cyclic calibration component graph".into()));
    }
    Ok((result, ordered_labels))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cyclic_blocks_form_one_group_in_dependency_order() -> Result<()> {
        // 0 と 2 を同時に解く。どちらのブロックからも読む 1 はその後になる。
        let edges = BTreeSet::from([(0, 2), (2, 0), (0, 1), (2, 1)]);
        let (components, labels) = ordered_components(4, &edges)?;
        assert_eq!(components, [vec![0, 2], vec![1], vec![3]]);
        assert_eq!(labels, [0, 1, 0, 2]);
        Ok(())
    }
}
`````

`````
crates/rhodium-pricing/src/calibration/solver.rs
//! 校正グループのソルバー（正方系と最小二乗）と、解における線形化の記録。

mod least_squares;

use std::{cell::Cell, collections::BTreeMap, sync::Arc};

use rhodium_core::{Error, Result};
use rhodium_market::{
    market::Market,
    quote::{QuoteId, QuoteSet},
    structure::MarketStructure,
};
use scioxide_ad::{Real, jacobian};
use scioxide_core::{
    InvalidInputReason, ResidualTolerance, SciOxideError, StepTolerance, argument,
};
use scioxide_linalg::{LuDecomposition, Matrix, RankPolicy};
use scioxide_rootfind::{
    ResidualFn, SystemRootResult, SystemRootSolver, TrustRegionSystem, TrustRegionSystemOptions,
};
use scioxide_sparse::CsrMatrix;

use super::{
    BindingRow, GroupSolution, PreparedBinding,
    groups::{GroupForm, GroupPlan},
};

/// 一つのグループの正方系を解き、自グループの θ と解の記録を返す。
///
/// `theta` は上流のグループまで解いた θ で、自グループの要素は使わない。
/// 必要な quote が欠けていれば解かずに [`Error::MissingData`] を返す。
/// 一意な系には前の θ を使う。それ以外は固定の初期値から解き、残差ノルムが最小の解を採る。
pub(super) fn solve(
    plan: &GroupPlan,
    rows: &[BindingRow],
    structure: &Arc<MarketStructure>,
    theta: &[f64],
    quotes: &QuoteSet,
    tolerance: ResidualTolerance<f64>,
    previous_theta: Option<&[f64]>,
) -> Result<(Vec<f64>, GroupSolution)> {
    let missing: Vec<_> = plan
        .quotes
        .iter()
        .filter(|id| quotes.get(id).is_none())
        .map(QuoteId::as_str)
        .collect();
    if !missing.is_empty() {
        return Err(Error::MissingData(format!(
            "missing calibration quotes: {}",
            missing.join(", ")
        )));
    }
    let targets = plan
        .rows
        .iter()
        .map(|&index| {
            // 準備に失敗した binding を含むグループは解かずに Failed にするので、ここには来ない。
            let row = rows[index].prepared().ok_or_else(|| {
                Error::Internal(format!(
                    "calibration instrument {} was not prepared",
                    rows[index].quote()
                ))
            })?;
            quotes
                .get(&row.quote)
                .map(|quote| (row, quote))
                .ok_or_else(|| {
                    Error::Internal(format!("calibration quote {} has no value", row.quote))
                })
        })
        .collect::<Result<_>>()?;
    let problem = Problem {
        plan,
        structure,
        theta,
        targets,
        tolerance,
    };
    if let GroupForm::LeastSquares(least_squares) = &plan.form {
        return problem.solve_least_squares(least_squares, previous_theta);
    }
    let direct = problem
        .targets
        .iter()
        .map(|&(row, quote)| row.direct_value(quote))
        .collect::<Result<Vec<_>>>()?;
    if direct.iter().all(Option::is_some) {
        let mut point = vec![0.0; plan.parameters.len()];
        for (index, value) in direct.into_iter().flatten() {
            let position = plan
                .parameters
                .binary_search(&index)
                .map_err(|_| Error::Internal("direct parameter outside its group".into()))?;
            point[position] = value;
        }
        return problem.linearize(&point, 0.0, 0, &vec![0.0; point.len()]);
    }
    let mut failures = Vec::new();
    let seeds = if let Some(warm) = previous_theta.and_then(|point| problem.warm_seed(point)) {
        match problem.solve_from(&warm) {
            Ok(solution) => return Ok(solution),
            Err(error) => failures.push(error.to_string()),
        }
        // 前の θ で収束しなければ、quote から作った初期値で一度だけ試す。
        vec![problem.flat_seed()?]
    } else {
        // 一意性が分からないときは二つの固定初期値で過去の校正結果から独立に解く。
        let flat = problem.flat_seed()?;
        let mut seeds = vec![vec![0.0; plan.parameters.len()]];
        if flat.iter().any(|value| *value != 0.0) {
            seeds.push(flat);
        }
        seeds
    };
    let mut best = None::<(Vec<f64>, GroupSolution)>;
    for seed in seeds {
        match problem.solve_from(&seed) {
            Ok((point, solution))
                if best.as_ref().is_none_or(|(_, previous)| {
                    solution.residual_norm < previous.residual_norm
                }) =>
            {
                best = Some((point, solution));
            }
            Ok(_) => {}
            Err(error) => failures.push(error.to_string()),
        }
    }
    best.ok_or_else(|| {
        Error::Calibration(format!(
            "group [{}]: {}",
            plan.blocks
                .iter()
                .map(|block| block.curve().as_str())
                .collect::<Vec<_>>()
                .join(", "),
            failures.join(" | ")
        ))
    })
}

/// 一つのグループの残差を、自グループ（と線形化のときは上流）の θ の関数として評価する。
struct Problem<'a> {
    plan: &'a GroupPlan,
    structure: &'a Arc<MarketStructure>,
    /// 上流まで解いた θ。未知数として渡さない要素は定数として読む。
    theta: &'a [f64],
    /// 残差の行の順に、校正商品とその quote の値。
    targets: Vec<(&'a PreparedBinding, f64)>,
    /// 解とみなす残差のユークリッドノルムの上限。各行の残差は、その行の目標の単位
    /// （パーレート、FX のレートの差、θ の要素）で測る。
    tolerance: ResidualTolerance<f64>,
}

impl Problem<'_> {
    /// 最新のピラーが各行で単調に効き、一つずつ一意に決まるときに前の θ を取り出す。
    fn warm_seed(&self, previous: &[f64]) -> Option<Vec<f64>> {
        if self.plan.blocks.len() != 1 {
            return None;
        }
        let mut pillars = self
            .targets
            .iter()
            .map(|&(row, quote)| {
                row.warm_start
                    .filter(|shape| !shape.requires_nonnegative_quote || quote >= 0.0)
                    .map(|shape| shape.pillar)
            })
            .collect::<Option<Vec<_>>>()?;
        pillars.sort_unstable();
        pillars
            .into_iter()
            .eq(0..self.plan.parameters.len())
            .then(|| {
                self.plan
                    .parameters
                    .iter()
                    .map(|&index| previous[index])
                    .collect()
            })
    }

    /// 各ブロックを、そのブロックに当てはめる quote の中央値 r で一定の金利にした初期値。
    ///
    /// ピラーの時刻 t での ln DF を −r·t とする。
    fn flat_seed(&self) -> Result<Vec<f64>> {
        let mut rates = BTreeMap::new();
        for &(row, quote) in &self.targets {
            rates
                .entry(&row.target)
                .or_insert_with(Vec::new)
                .push(row.seed_rate(quote));
        }
        let mut seed = Vec::with_capacity(self.plan.parameters.len());
        for block in &self.plan.blocks {
            let rates = rates.get_mut(block).ok_or_else(|| {
                Error::Internal(format!("no quotes determine curve {}", block.curve()))
            })?;
            rates.sort_by(f64::total_cmp);
            let middle = rates.len() / 2;
            let median = if rates.len().is_multiple_of(2) {
                f64::midpoint(rates[middle - 1], rates[middle])
            } else {
                rates[middle]
            };
            let times = self.structure.block_times(block).ok_or_else(|| {
                Error::Internal(format!("no pillar times for curve {}", block.curve()))
            })?;
            let count = self
                .structure
                .layout()
                .range(block)
                .ok_or_else(|| Error::Internal(format!("no range for block {block:?}")))?
                .len();
            if times.len() < count {
                // 会合の前後の区間ごとに、瞬間フォワードそのものを持つ。
                seed.extend(std::iter::repeat_n(median, count));
            } else {
                seed.extend(times.iter().map(|time| -median * time.years()));
            }
        }
        Ok(seed)
    }

    /// θ の添字 `active` に `values` を入れた市場で、残差 F = par(θ) − q を行の順に返す。
    fn residuals<R: Real>(&self, active: &[usize], values: &[R]) -> Result<Vec<R>> {
        // 未知数と必要な上流だけが AD の入力。他の値は定数として市場を組み立てる。
        let mut theta: Vec<_> = self.theta.iter().copied().map(R::constant).collect();
        for (&index, &value) in active.iter().zip(values) {
            theta[index] = value;
        }
        let available = self.plan.blocks.iter().chain(&self.plan.upstream).cloned();
        let market = Market::from_partial_theta(Arc::clone(self.structure), theta, available)?;
        self.targets
            .iter()
            .map(|&(row, quote)| {
                row.residual(&market, quote)
                    .map_err(|error| error.context(format_args!("quote {}", row.quote)))
            })
            .collect()
    }

    /// 一つの初期値から trust-region Newton で自グループの θ を解き、解の記録を作る。
    ///
    /// ヤコビアンは近似せず、残差を AD で微分して渡す。
    fn solve_from(&self, seed: &[f64]) -> Result<(Vec<f64>, GroupSolution)> {
        let n = seed.len();
        // ln DF のステップが小さくても、短期金利の残差は許容値に達していないことがある。
        let step = StepTolerance::new(0.0, 0.0).map_err(calibration_error)?;
        let solver = TrustRegionSystem::new(
            TrustRegionSystemOptions::default()
                .with_residual_tolerance(self.tolerance)
                .with_step_tolerance(step),
        )
        .map_err(calibration_error)?;
        let active = &self.plan.parameters;
        let callback_error = Cell::new(None);
        let function = ResidualFn::new(n, |x: &[f64], out: &mut [f64]| {
            let residuals = forward_error(self.residuals(active, x), &callback_error)?;
            out.copy_from_slice(&residuals);
            Ok(())
        })
        .with_jacobian(|x: &[f64], out: &mut [f64]| {
            let matrix = jacobian(
                |values, residuals| {
                    residuals.extend(forward_error(
                        self.residuals(active, values),
                        &callback_error,
                    )?);
                    Ok(())
                },
                x,
            )?;
            out.copy_from_slice(matrix.as_slice());
            Ok(())
        });
        // 残差の評価で失敗したときは、ソルバーの汎用のエラーより元のエラーを返す。
        let root = solver.solve(&function, seed).map_err(|error| {
            callback_error
                .take()
                .unwrap_or_else(|| calibration_error(error))
        })?;
        self.record(&root)
    }

    /// 収束を確かめ、解の点で線形化した [`GroupSolution`] と自グループの θ を返す。
    ///
    /// 残差ノルムが許容誤差を超えるか、∂G/∂θ が数値的にランク落ちしていれば、
    /// 解を識別できないので失敗とする。
    fn record(&self, root: &SystemRootResult<f64>) -> Result<(Vec<f64>, GroupSolution)> {
        let residual_norm = root.residual_norm().unwrap_or(f64::INFINITY);
        if !root.is_success()
            || !residual_norm.is_finite()
            || residual_norm > self.tolerance.absolute()
        {
            return Err(Error::Calibration(format!(
                "Newton {:?}, residual norm {residual_norm:e}{}",
                root.termination(),
                self.largest_residual(root.residual())
            )));
        }
        self.linearize(
            root.point(),
            residual_norm,
            root.iterations(),
            root.residual(),
        )
    }

    fn linearize(
        &self,
        local: &[f64],
        residual_norm: f64,
        iterations: usize,
        residual: &[f64],
    ) -> Result<(Vec<f64>, GroupSolution)> {
        let n = self.plan.parameters.len();
        // 列の順は、自グループの n 個のあとに上流の要素が続く。
        let mut active = self.plan.parameters.clone();
        active.extend_from_slice(&self.plan.upstream_parameters);
        let mut point = local.to_vec();
        point.extend(
            self.plan
                .upstream_parameters
                .iter()
                .map(|&index| self.theta[index]),
        );

        // ソルバーの内部近似に頼らず、解で自グループと上流の列を同時に線形化する。
        let callback_error = Cell::new(None);
        let jac = jacobian(
            |values, residuals| {
                residuals.extend(forward_error(
                    self.residuals(&active, values),
                    &callback_error,
                )?);
                Ok(())
            },
            &point,
        )
        .map_err(|error| {
            callback_error
                .take()
                .unwrap_or_else(|| calibration_error(error))
        })?;

        // 行優先のヤコビアンから、自グループの列だけを取り出した正方行列を分解する。
        let own = jac
            .as_slice()
            .chunks_exact(active.len())
            .flat_map(|row| row[..n].iter().copied())
            .collect();
        let own = Matrix::from_vec(n, n, own).map_err(calibration_error)?;
        let dg_dtheta = LuDecomposition::factor(own.view()).map_err(|error| {
            Error::Calibration(format!(
                "residual Jacobian LU: {error}{}",
                self.largest_residual(residual)
            ))
        })?;
        let rank = dg_dtheta
            .analyze_rank(RankPolicy::DimensionScaledEpsilon)
            .map_err(calibration_error)?
            .rank();
        if rank < n {
            let rcond = dg_dtheta.condition_estimate().map_or_else(
                |error| error.to_string(),
                |estimate| format!("{estimate:?}"),
            );
            return Err(Error::Calibration(format!(
                "residual Jacobian rank {rank}/{n}, reciprocal condition {rcond}{}",
                self.largest_residual(residual)
            )));
        }
        let dg_dup = upstream_jacobian(n, active.len(), jac.as_slice())?;
        let dg_dq = self.quote_jacobian()?;
        // 正方系の停止の誤差は、残差を一回のニュートン補正で θ に戻した δθ = J⁻¹F。
        let theta_error = dg_dtheta
            .solve(residual)
            .map_err(calibration_error)?
            .as_slice()
            .to_vec();
        point.truncate(n);
        Ok((
            point,
            GroupSolution {
                dg_dtheta,
                dg_dq,
                dg_dup,
                residual_norm,
                iterations,
                hessian: None,
                residuals: residual.to_vec(),
                theta_error,
            },
        ))
    }

    /// ∂G/∂q を作る。行は残差の順、列は `plan.quotes` の順。
    fn quote_jacobian(&self) -> Result<CsrMatrix<f64>> {
        let columns = self
            .targets
            .iter()
            .map(|&(row, _)| {
                self.plan.quotes.binary_search(&row.quote).map_err(|_| {
                    Error::Internal(format!(
                        "calibration group does not record quote {}",
                        row.quote
                    ))
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let n = columns.len();
        // 各残差が一つの入力 quote を読み、変換した値を引く。
        CsrMatrix::from_triplets(
            n,
            self.plan.quotes.len(),
            &(0..n).collect::<Vec<_>>(),
            &columns,
            &self
                .targets
                .iter()
                .map(|&(row, quote)| -row.quote_derivative(quote))
                .collect::<Vec<_>>(),
        )
        .map_err(calibration_error)
    }

    /// 失敗の報告に添える、絶対値が最大の残差とその quote。
    fn largest_residual(&self, residuals: &[f64]) -> String {
        self.targets
            .iter()
            .zip(residuals)
            .max_by(|(_, lhs), (_, rhs)| lhs.abs().total_cmp(&rhs.abs()))
            .map_or_else(String::new, |((row, _), residual)| {
                format!(", largest residual {}={residual:e}", row.quote)
            })
    }
}

/// `columns` 列の行優先のヤコビアンから、先頭 n 列を除いた上流の列を疎行列にする。
///
/// 残差が読まない上流の節点の列は厳密にゼロなので、格納しない。
fn upstream_jacobian(n: usize, columns: usize, jacobian: &[f64]) -> Result<CsrMatrix<f64>> {
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut values = Vec::new();
    for (row, data) in jacobian.chunks_exact(columns).enumerate() {
        for (col, &value) in data[n..]
            .iter()
            .enumerate()
            .filter(|&(_, &value)| value != 0.0)
        {
            rows.push(row);
            cols.push(col);
            values.push(value);
        }
    }
    CsrMatrix::from_triplets(n, columns - n, &rows, &cols, &values).map_err(calibration_error)
}

#[expect(clippy::needless_pass_by_value, reason = "map_err の変換関数に使う")]
fn calibration_error(error: SciOxideError) -> Error {
    Error::Calibration(error.to_string())
}

/// scioxide のコールバックは [`SciOxideError`] しか返せないので、
/// 元の [`Error`] を `failure` に退避し、ソルバーが失敗したときに取り出せるようにする。
fn forward_error<T>(
    result: Result<T>,
    failure: &Cell<Option<Error>>,
) -> std::result::Result<T, SciOxideError> {
    result.map_err(|error| {
        failure.set(Some(error));
        SciOxideError::invalid_input(
            argument!("theta"),
            InvalidInputReason::OutOfRange,
            "calibration residual evaluation failed",
        )
    })
}
`````

`````
crates/rhodium-pricing/src/calibration/solver/least_squares.rs
//! 過剰決定と罰則付きの当てはめを最小二乗で解き、停留条件の線形化を記録する。
//!
//! 目的は `Φ(θ) = ½‖r̃(θ)‖²` で、`r̃` は重みを掛けた残差 `w_i·F_i`（`w_i = 1/σ̃_i`）と
//! 罰則の行 `√λ·(f_{i+1} − f_i)/s` を並べたもの。停留条件 `G = ∇Φ = J̃ᵀ r̃ = J_Fᵀ W F + ∇P` が
//! グループの方程式で、解の記録は G の線形化（`∂G/∂θ`、`∂G/∂q`、`∂G/∂θ_up`）にする。

use std::cell::Cell;

use rhodium_core::{Error, Result};
use scioxide_ad::{Real, hvp, jacobian};
use scioxide_core::Budget;
use scioxide_linalg::{
    CholeskyDecomposition, CholeskyPivoting, HermitianInput, HermitianProjectionPolicy,
    LuDecomposition, Matrix, RankPolicy,
};
use scioxide_optimize::{LeastSquares, LeastSquaresFn, LeastSquaresOptions};
use scioxide_sparse::CsrMatrix;

use super::{Problem, calibration_error, forward_error};
use crate::calibration::{GroupSolution, HessianMode, groups::LeastSquaresPlan};

/// LM の後に、厳密なヘシアンのニュートン補正で停留条件を磨く回数の上限。
///
/// 残差が残る最小二乗では LM の収束は一次にとどまり、校正の許容誤差の手前で止まることがある。
/// 停留点の近くのニュートン補正は二次で収束するので、数回で倍精度の分解能に届く。
const POLISH_STEPS: usize = 8;

/// 一つの点での最小二乗のグループの線形化。
struct Linearization {
    solution: GroupSolution,
    /// 停留条件の誤差 δθ = (∂G/∂θ)⁻¹·G が各行を動かす量の、行の自然な単位でのノルム ‖Ĵ δθ‖。
    movement: f64,
    /// 目的 Φ = ½‖r̃‖²。初期値の候補から解を選ぶのに使う。
    cost: f64,
}

/// 行の二乗和の半分 `½ Σ r_i²`。
fn objective<R: Real>(rows: Vec<R>) -> R {
    rows.into_iter()
        .fold(R::constant(0.0), |sum, row| sum + row * row)
        * R::constant(0.5)
}

impl Problem<'_> {
    /// 最小二乗のグループを解き、停留点での線形化を返す。
    ///
    /// 目的が強凸（直接マークだけ）なら前の θ から解き、収束しなければ quote からの初期値で一度だけ
    /// 解き直す。それ以外は固定の初期値の集合から解き、目的の値が最小の解を選ぶ。
    pub(super) fn solve_least_squares(
        &self,
        plan: &LeastSquaresPlan,
        previous: Option<&[f64]>,
    ) -> Result<(Vec<f64>, GroupSolution)> {
        let mut failures = Vec::new();
        let flat = self.flat_seed()?;
        let seeds = if let Some(theta) = previous.filter(|_| plan.convex) {
            let warm: Vec<_> = self
                .plan
                .parameters
                .iter()
                .map(|&index| theta[index])
                .collect();
            match self.least_squares_from(plan, &warm) {
                Ok((point, linear)) => return Ok((point, linear.solution)),
                Err(error) => failures.push(error.to_string()),
            }
            vec![flat]
        } else {
            let mut seeds = vec![vec![0.0; self.plan.parameters.len()]];
            if flat.iter().any(|value| *value != 0.0) {
                seeds.push(flat);
            }
            seeds
        };
        let mut best = None::<(Vec<f64>, Linearization)>;
        for seed in seeds {
            match self.least_squares_from(plan, &seed) {
                Ok((point, linear))
                    if best
                        .as_ref()
                        .is_none_or(|(_, previous)| linear.cost < previous.cost) =>
                {
                    best = Some((point, linear));
                }
                Ok(_) => {}
                Err(error) => failures.push(error.to_string()),
            }
        }
        best.map(|(point, linear)| (point, linear.solution))
            .ok_or_else(|| {
                Error::Calibration(format!(
                    "least-squares group [{}]: {}",
                    self.plan
                        .blocks
                        .iter()
                        .map(|block| block.curve().as_str())
                        .collect::<Vec<_>>()
                        .join(", "),
                    failures.join(" | ")
                ))
            })
    }

    /// 重みを掛けた残差 `w_i·F_i` と罰則の行を、θ の添字 `active` に `values` を入れて求める。
    ///
    /// 罰則の行の添字は自グループの未知数の中の位置で、`active` の先頭が自グループの未知数。
    fn weighted_rows<R: Real>(
        &self,
        plan: &LeastSquaresPlan,
        active: &[usize],
        values: &[R],
    ) -> Result<Vec<R>> {
        let residuals = self.residuals(active, values)?;
        let data = residuals
            .into_iter()
            .zip(&plan.weights)
            .map(|(residual, &weight)| residual * R::constant(weight));
        let penalties = plan
            .penalties
            .iter()
            .map(|row| (values[row.upper] - values[row.lower]) * R::constant(row.factor));
        Ok(data.chain(penalties).collect())
    }

    /// 一つの初期値から Levenberg–Marquardt で解き、厳密なヘシアンで停留条件まで磨く。
    fn least_squares_from(
        &self,
        plan: &LeastSquaresPlan,
        seed: &[f64],
    ) -> Result<(Vec<f64>, Linearization)> {
        let own = &self.plan.parameters;
        let rows = plan.weights.len() + plan.penalties.len();
        let callback_error = Cell::new(None);
        let function = LeastSquaresFn::new(seed.len(), rows, |x: &[f64], out: &mut [f64]| {
            let values = forward_error(self.weighted_rows(plan, own, x), &callback_error)?;
            out.copy_from_slice(&values);
            Ok(())
        })
        .map_err(calibration_error)?
        .with_jacobian(|x: &[f64], out: &mut Matrix<f64>| {
            let matrix = jacobian(
                |values, rows| {
                    rows.extend(forward_error(
                        self.weighted_rows(plan, own, values),
                        &callback_error,
                    )?);
                    Ok(())
                },
                x,
            )?;
            out.as_mut_slice().copy_from_slice(matrix.as_slice());
            Ok(())
        });
        let solver = LeastSquares::new(
            LeastSquaresOptions::default().with_max_iterations(Budget::ONE.saturating_add(199)),
        )
        .map_err(calibration_error)?;
        let result = solver.solve(&function, seed).map_err(|error| {
            callback_error
                .take()
                .unwrap_or_else(|| calibration_error(error))
        })?;
        // 収束は LM の終了理由ではなく、厳密なヘシアンの停留条件で判定する。
        let mut point = result.x().to_vec();
        for step in 0..=POLISH_STEPS {
            let linear = self.linearize_least_squares(plan, &point, result.iterations() + step)?;
            if linear.movement <= self.tolerance.absolute() {
                return Ok((point, linear));
            }
            for (value, error) in point.iter_mut().zip(linear.solution.theta_error()) {
                *value -= error;
            }
        }
        Err(Error::Calibration(format!(
            "least-squares stationarity not reached after LM {:?}",
            result.termination()
        )))
    }

    /// 点 `local` で停留条件 G を線形化し、解の記録と停止の誤差を作る。
    ///
    /// 列の順は自グループの n 個のあとに上流の要素が続く。厳密な記録は ∇²Φ の自グループの n 列を
    /// Hessian–vector product で求め、その上の n 行が `∂G/∂θ`、下の行が（対称性から）`∂G/∂θ_up` の転置。
    /// `∇²Φ = J̃ᵀJ̃ + Σ_i r̃_i ∇²r̃_i` なので、これは `J_Fᵀ W J_F + Σ_i (WF)_i ∇²F_i + ∇²P` と
    /// 上流との交差項 `J_Fᵀ W ∂F/∂θ_up + Σ_i (WF)_i ∂²F_i/∂θ_k∂θ_up` を同時に与える。
    fn linearize_least_squares(
        &self,
        plan: &LeastSquaresPlan,
        local: &[f64],
        iterations: usize,
    ) -> Result<Linearization> {
        let n = self.plan.parameters.len();
        let mut active = self.plan.parameters.clone();
        active.extend_from_slice(&self.plan.upstream_parameters);
        let columns = active.len();
        let mut point = local.to_vec();
        point.extend(
            self.plan
                .upstream_parameters
                .iter()
                .map(|&index| self.theta[index]),
        );
        let callback_error = Cell::new(None);
        let lift = |error| {
            callback_error
                .take()
                .unwrap_or_else(|| calibration_error(error))
        };
        let jac = jacobian(
            |values, rows| {
                rows.extend(forward_error(
                    self.weighted_rows(plan, &active, values),
                    &callback_error,
                )?);
                Ok(())
            },
            &point,
        )
        .map_err(lift)?;
        let weighted = self.weighted_rows::<f64>(plan, &active, &point)?;
        let entry = |row: usize, column: usize| jac.as_slice()[row * columns + column];
        let gradient: Vec<f64> = (0..n)
            .map(|j| (0..weighted.len()).map(|i| entry(i, j) * weighted[i]).sum())
            .collect();
        let (hessian, cross) =
            self.second_derivatives(plan, &active, &point, jac.as_slice(), &weighted)?;
        let matrix = Matrix::from_vec(n, n, hessian).map_err(calibration_error)?;
        // 停留点が極小でない（鞍点や極大）か、当てはめで識別できない向きがあれば解とみなさない。
        // 分解は正定値の判定だけに使い、随伴には同じ行列の LU 分解を残す（正方系と同じ記録の型）。
        let _definite = CholeskyDecomposition::factor(
            matrix.view(),
            HermitianInput::Projected(HermitianProjectionPolicy::Always),
            CholeskyPivoting::None,
        )
        .map_err(|error| {
            Error::Calibration(format!(
                "least-squares Hessian is not positive definite: {error}"
            ))
        })?;
        let dg_dtheta = LuDecomposition::factor(matrix.view()).map_err(calibration_error)?;
        let rank = dg_dtheta
            .analyze_rank(RankPolicy::DimensionScaledEpsilon)
            .map_err(calibration_error)?
            .rank();
        if rank < n {
            return Err(Error::Calibration(format!(
                "least-squares Hessian rank {rank}/{n}"
            )));
        }
        let theta_error = dg_dtheta
            .solve(gradient.as_slice())
            .map_err(calibration_error)?
            .as_slice()
            .to_vec();
        let data = plan.weights.len();
        let movement = movement(plan, jac.as_slice(), columns, &theta_error);
        let residuals: Vec<f64> = weighted[..data]
            .iter()
            .zip(&plan.weights)
            .map(|(row, weight)| row / weight)
            .collect();
        let residual_norm = residuals
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        let cost = objective(weighted);
        let solution = GroupSolution {
            dg_dtheta,
            dg_dq: self.least_squares_quote_jacobian(plan, n, columns, jac.as_slice())?,
            dg_dup: dense_to_csr(n, columns - n, &cross)?,
            residual_norm,
            iterations,
            hessian: Some(plan.hessian),
            residuals,
            theta_error,
        };
        Ok(Linearization {
            solution,
            movement,
            cost,
        })
    }

    /// 停留条件 G の自グループと上流の θ についての微分を、行優先の n × n と n × u の配列で返す。
    ///
    /// `hessian[i·n + j] = ∂G_i/∂θ_j`、`cross[j·u + k] = ∂G_j/∂θ_up,k`。
    fn second_derivatives(
        &self,
        plan: &LeastSquaresPlan,
        active: &[usize],
        point: &[f64],
        jacobian: &[f64],
        weighted: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        let n = self.plan.parameters.len();
        let columns = active.len();
        let upstream = columns - n;
        let entry = |row: usize, column: usize| jacobian[row * columns + column];
        let callback_error = Cell::new(None);
        let mut hessian = vec![0.0; n * n];
        let mut cross = vec![0.0; n * upstream];
        for j in 0..n {
            let column: Vec<f64> = match plan.hessian {
                HessianMode::Exact => {
                    let mut direction = vec![0.0; columns];
                    direction[j] = 1.0;
                    hvp(
                        |values| {
                            Ok(objective(forward_error(
                                self.weighted_rows(plan, active, values),
                                &callback_error,
                            )?))
                        },
                        point,
                        &direction,
                    )
                    .map_err(|error| {
                        callback_error
                            .take()
                            .unwrap_or_else(|| calibration_error(error))
                    })?
                    .as_slice()
                    .to_vec()
                }
                // 残差の 2 階の項を省く。罰則の行は θ について線形なので ∇²P は厳密に残る。
                HessianMode::GaussNewton => (0..columns)
                    .map(|k| (0..weighted.len()).map(|i| entry(i, k) * entry(i, j)).sum())
                    .collect(),
            };
            for (i, &value) in column[..n].iter().enumerate() {
                hessian[i * n + j] = value;
            }
            cross[j * upstream..(j + 1) * upstream].copy_from_slice(&column[n..]);
        }
        Ok((hessian, cross))
    }

    /// `∂G/∂q = −J_Fᵀ W D`（行は自グループの θ、列は `plan.quotes`）。
    ///
    /// `r̃_i = w_i·(m_i(θ) − D_i·q_i)` なので `∂r̃_i/∂q_i = −w_i·D_i` で、`J̃` は q に依存しない。
    fn least_squares_quote_jacobian(
        &self,
        plan: &LeastSquaresPlan,
        n: usize,
        columns: usize,
        jacobian: &[f64],
    ) -> Result<CsrMatrix<f64>> {
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut values = Vec::new();
        for (i, (&(row, quote), weight)) in self.targets.iter().zip(&plan.weights).enumerate() {
            let column = self.plan.quotes.binary_search(&row.quote).map_err(|_| {
                Error::Internal(format!(
                    "calibration group does not record quote {}",
                    row.quote
                ))
            })?;
            let scale = -weight * row.quote_derivative(quote);
            for j in 0..n {
                let value = jacobian[i * columns + j] * scale;
                if value != 0.0 {
                    rows.push(j);
                    cols.push(column);
                    values.push(value);
                }
            }
        }
        // 同じ quote を読む行（同じ要素への重ねた直接マークなど）の寄与は座標ごとに足される。
        CsrMatrix::from_triplets(n, self.plan.quotes.len(), &rows, &cols, &values)
            .map_err(calibration_error)
    }
}

/// 停留条件の誤差 δθ が各行を動かす量のノルム ‖Ĵ δθ‖。
///
/// 行を自然な単位（quote の残差、罰則はフォワードの差）に戻して測る。正方系では J δθ = F なので、
/// この量は残差のノルムに一致し、校正の許容誤差と同じ尺度で比べられる。
fn movement(plan: &LeastSquaresPlan, jacobian: &[f64], columns: usize, theta_error: &[f64]) -> f64 {
    plan.weights
        .iter()
        .copied()
        .chain(plan.penalties.iter().map(|row| row.factor))
        .enumerate()
        .map(|(i, scale)| {
            let moved: f64 = theta_error
                .iter()
                .enumerate()
                .map(|(j, error)| jacobian[i * columns + j] * error)
                .sum();
            (moved / scale).powi(2)
        })
        .sum::<f64>()
        .sqrt()
}

/// 行優先の密行列のうちゼロでない要素を疎行列にする。
fn dense_to_csr(rows: usize, columns: usize, values: &[f64]) -> Result<CsrMatrix<f64>> {
    let (mut row_indices, mut column_indices, mut entries) = (Vec::new(), Vec::new(), Vec::new());
    for (index, &value) in values
        .iter()
        .enumerate()
        .filter(|(_, value)| **value != 0.0)
    {
        row_indices.push(index / columns.max(1));
        column_indices.push(index % columns.max(1));
        entries.push(value);
    }
    CsrMatrix::from_triplets(rows, columns, &row_indices, &column_indices, &entries)
        .map_err(calibration_error)
}
`````

`````
crates/rhodium-pricing/src/config.rs
//! 価格設定と方式の選択。

use std::{
    collections::{BTreeMap, BTreeSet, btree_map::Entry},
    fmt,
    marker::PhantomData,
};

use rhodium_core::{
    Error, Result, Warning, currency::Currency, id::IndexId, refdata::RefData, time::LocalEventTime,
};
use rhodium_market::curve::CurveId;
use rhodium_products::{
    ir::{SettlementBoundary, analysis::Linearity},
    standard::ProductKind,
};
use serde::{Deserialize, Deserializer, de};

use crate::{definition::MarketDefinition, resolve::CollateralKey};

/// 価格設定の割引表の一行。
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiscountRule {
    /// 支払いの通貨。
    pub currency: Currency,
    /// 支払いに適用する担保条件。
    pub collateral: CollateralKey,
    /// 指定した通貨と担保条件で使う割引曲線。
    pub curve: CurveId,
}

/// 資金調達費用を数える場所。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FundingTreatment {
    /// FCA と FBA で数える。
    #[default]
    InXva,
    /// 無担保の支払いの割引に入れる。
    InDiscounting,
}

/// 通貨ごとの自社の資金調達曲線。
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FundingCurves {
    /// 自社の借入に使う曲線。
    pub borrow: CurveId,
    /// 自社の運用に使う曲線。
    pub lend: CurveId,
}

/// 商品の評価方式。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Method {
    /// cashflow の期待値と割引の積和で評価する。
    Linear,
    /// 型付き標準商品から解析評価する。
    Analytic(AnalyticPricer),
}

/// 組み込みの解析プライサー。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AnalyticPricer {
    /// 値洗い口座の測度で金利先物を評価する。
    RateFuture,
}

/// リニア評価で使える名前付きの近似（設計書の「近似の規約」）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LinearConvention {
    /// 支払日のずれに伴うタイミング調整を省く。
    PaymentDelayIgnored,
    /// 割引と予想が別曲線のとき、担保の違いによる予想の調整を省く。
    CrossCollateralProjection,
    /// 日次因子の積で成長率を予想する。
    DailyForwardProduct,
    /// 契約の丸めを予想値に適用する。
    RoundingAtForward,
    /// FRA の予想曲線と割引曲線の差に伴う凸性を省く。
    FraDiscounting,
    /// FX リセットと金利の共分散を省く。
    MtmFxForwardProduct,
}

/// リニア評価で許す近似と営業日数の上限。
#[derive(Debug, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct LinearConventions {
    disabled: BTreeSet<LinearConvention>,
    max_payment_delay_business_days: u32,
    max_lookback_lockout_business_days: u8,
}

impl Default for LinearConventions {
    fn default() -> Self {
        Self {
            disabled: BTreeSet::new(),
            max_payment_delay_business_days: 10,
            max_lookback_lockout_business_days: 5,
        }
    }
}

impl LinearConventions {
    /// 設定で指定の近似を許しているか。
    #[must_use]
    pub fn allows(&self, convention: LinearConvention) -> bool {
        !self.disabled.contains(&convention)
    }

    /// 観測期間の終了と支払日の間に許す営業日数。
    #[must_use]
    pub const fn max_payment_delay_business_days(&self) -> u32 {
        self.max_payment_delay_business_days
    }

    /// 日次積における lookback または lockout に許す営業日数。
    #[must_use]
    pub const fn max_lookback_lockout_business_days(&self) -> u8 {
        self.max_lookback_lockout_business_days
    }

    pub(crate) fn require(&self, convention: LinearConvention) -> Result<Warning> {
        if !self.allows(convention) {
            return Err(Error::Unsupported(format!(
                "linear convention {convention:?} is disabled"
            )));
        }
        Ok(match convention {
            LinearConvention::PaymentDelayIgnored => Warning::PaymentDelayIgnored,
            LinearConvention::CrossCollateralProjection => Warning::CrossCollateralProjection,
            LinearConvention::DailyForwardProduct => Warning::DailyForwardProduct,
            LinearConvention::RoundingAtForward => Warning::RoundingAtForward,
            LinearConvention::FraDiscounting => Warning::FraDiscounting,
            LinearConvention::MtmFxForwardProduct => Warning::MtmFxForwardProduct,
        })
    }
}

/// 商品の種類を評価方式に対応付ける一行。
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MethodRule {
    /// 対象の標準商品の種類。
    pub product: ProductKind,
    /// この商品に適用する評価方式。
    pub method: Method,
}

/// Serde で読み込む価格設定。文書の版は 1 で、報告通貨は必須。
///
/// 割引表は支払い通貨と担保条件の組ごとに一行、方式の表は商品の種類ごとに一行。
/// 市場の定義に依存する曲線の照合は [`crate::resolve::resolve`] が行う。
#[derive(Debug, Deserialize)]
#[serde(try_from = "PricingConfigDocument")]
pub struct PricingConfig(PricingConfigDocument);

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct PricingConfigDocument {
    version: u32,
    reporting_currency: Currency,
    #[serde(default)]
    valuation_close: Option<LocalEventTime>,
    #[serde(default)]
    discounting: Vec<DiscountRule>,
    #[serde(default, deserialize_with = "unique_rules")]
    forecasting: BTreeMap<IndexId, CurveId>,
    #[serde(default, deserialize_with = "unique_rules")]
    fx_collateral: BTreeMap<rhodium_core::fx::CcyPair, CurveId>,
    #[serde(default)]
    methods: Vec<MethodRule>,
    #[serde(default)]
    conventions: LinearConventions,
    #[serde(default)]
    settlement: SettlementBoundary,
    #[serde(default)]
    funding: FundingTreatment,
    #[serde(default, deserialize_with = "unique_rules")]
    funding_curves: BTreeMap<Currency, FundingCurves>,
}

fn unique_rules<'de, D, K, V>(deserializer: D) -> std::result::Result<BTreeMap<K, V>, D::Error>
where
    D: Deserializer<'de>,
    K: Deserialize<'de> + Ord + fmt::Display,
    V: Deserialize<'de>,
{
    struct Rules<K, V>(PhantomData<(K, V)>);

    impl<'de, K: Deserialize<'de> + Ord + fmt::Display, V: Deserialize<'de>> de::Visitor<'de>
        for Rules<K, V>
    {
        type Value = BTreeMap<K, V>;

        fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            formatter.write_str("a rule table with unique keys")
        }

        fn visit_map<A: de::MapAccess<'de>>(
            self,
            mut map: A,
        ) -> std::result::Result<Self::Value, A::Error> {
            let mut rules = BTreeMap::new();
            while let Some(key) = map.next_key()? {
                match rules.entry(key) {
                    Entry::Vacant(entry) => {
                        entry.insert(map.next_value()?);
                    }
                    Entry::Occupied(entry) => {
                        return Err(de::Error::custom(format!(
                            "duplicate rule for {}",
                            entry.key()
                        )));
                    }
                }
            }
            Ok(rules)
        }
    }

    // JSON の同名キーも、後の行で黙って上書きできないようにする。
    deserializer.deserialize_map(Rules(PhantomData))
}

impl TryFrom<PricingConfigDocument> for PricingConfig {
    type Error = Error;

    fn try_from(document: PricingConfigDocument) -> Result<Self> {
        if document.version != 1 {
            return Err(Error::InvalidData(format!(
                "unsupported pricing document version {}",
                document.version
            )));
        }
        // 割引規則は支払い通貨と担保の鍵を比べて重複を探す。率の異なる単位表記と
        // CTD の候補順は、担保の鍵を読み込むときにそろえてある。
        for (position, rule) in document.discounting.iter().enumerate() {
            if document.discounting[..position]
                .iter()
                .any(|other| rule.currency == other.currency && rule.collateral == other.collateral)
            {
                return Err(Error::InvalidData(format!(
                    "duplicate discounting rule for {} under {:?}",
                    rule.currency, rule.collateral
                )));
            }
        }
        for (position, rule) in document.methods.iter().enumerate() {
            if document.methods[..position]
                .iter()
                .any(|other| rule.product == other.product)
            {
                return Err(Error::InvalidData(format!(
                    "duplicate method rule for {:?}",
                    rule.product
                )));
            }
        }
        Ok(Self(document))
    }
}

impl PricingConfig {
    /// 外部の設定が参照する通貨、地域、指標、曲線の ID を照合する。
    ///
    /// # Errors
    ///
    /// 指定された ID が参照データまたは市場定義に無ければ `InvalidData`。
    pub fn validate_references(&self, refdata: &RefData, market: &MarketDefinition) -> Result<()> {
        let require_curve = |id: &CurveId| {
            if market.contains_curve(id) {
                Ok(())
            } else {
                Err(Error::InvalidData(format!("unknown curve {id}")))
            }
        };
        let mut errors = Vec::new();
        let mut record = |field: &str, result: Result<()>| {
            if let Err(error) = result {
                errors.push(format!("{field}: {error}"));
            }
        };
        record(
            "reporting_currency",
            refdata.currency(self.reporting_currency()).map(|_| ()),
        );
        if let Some(close) = self.valuation_close() {
            record(
                "valuation_close.zone",
                refdata.zone(&close.zone).map(|_| ()),
            );
        }
        for (i, rule) in self.0.discounting.iter().enumerate() {
            record(
                &format!("discounting[{i}]"),
                (|| {
                    refdata.currency(rule.currency)?;
                    rule.collateral.validate_references(refdata)?;
                    require_curve(&rule.curve)
                })(),
            );
        }
        for (index, curve) in &self.0.forecasting {
            record(
                &format!("forecasting.{index}"),
                (|| {
                    refdata.index(index)?;
                    require_curve(curve)
                })(),
            );
        }
        for (&ccy, curves) in &self.0.funding_curves {
            record(
                &format!("funding_curves.{ccy}"),
                (|| {
                    refdata.currency(ccy)?;
                    require_curve(&curves.borrow)?;
                    require_curve(&curves.lend)
                })(),
            );
        }
        if !errors.is_empty() {
            return Err(Error::InvalidData(errors.join("\n")));
        }
        Ok(())
    }

    /// PV を集計する報告通貨。
    #[must_use]
    pub const fn reporting_currency(&self) -> Currency {
        self.0.reporting_currency
    }

    /// 評価の日付を評価の瞬間に解決する締め時刻と地域。
    #[must_use]
    pub const fn valuation_close(&self) -> Option<&LocalEventTime> {
        self.0.valuation_close.as_ref()
    }

    /// 評価時点ちょうどの受渡の扱い。
    #[must_use]
    pub const fn settlement(&self) -> SettlementBoundary {
        self.0.settlement
    }

    /// リニア評価で許す近似と上限。
    #[must_use]
    pub const fn conventions(&self) -> &LinearConventions {
        &self.0.conventions
    }

    /// 資金調達費用を数える場所。
    #[must_use]
    pub const fn funding(&self) -> FundingTreatment {
        self.0.funding
    }

    /// 割引と XVA が共通に参照する通貨ごとの資金調達曲線。
    #[must_use]
    pub const fn funding_curves(&self) -> &BTreeMap<Currency, FundingCurves> {
        &self.0.funding_curves
    }

    /// 方式の明示行を優先し、指定がなければ線形な IR をリニア方式で評価する。
    ///
    /// `product` は標準商品の種類。`linearity` は簡約済みの IR にリニア方式を
    /// 適用できるかで、呼び出し側が判定して渡す。
    ///
    /// # Errors
    /// 選んだ方式を適用できない場合や、適用できる方式が指定されていない場合は
    /// `Unsupported`。
    pub fn select_method(
        &self,
        product: Option<ProductKind>,
        linearity: Linearity,
    ) -> Result<Method> {
        let explicit = self
            .0
            .methods
            .iter()
            .find(|rule| Some(rule.product) == product);
        match (explicit.map(|rule| rule.method), product, linearity) {
            (
                Some(Method::Analytic(AnalyticPricer::RateFuture)) | None,
                Some(ProductKind::RateFuture),
                _,
            ) => Ok(Method::Analytic(AnalyticPricer::RateFuture)),
            (Some(Method::Analytic(_)), _, _) => Err(Error::Unsupported(
                "analytic pricer does not support this product".into(),
            )),
            (Some(Method::Linear) | None, _, Linearity::Linear) => Ok(Method::Linear),
            (Some(Method::Linear), _, Linearity::NonLinear(_)) => Err(Error::Unsupported(
                "configured linear method cannot evaluate this product".into(),
            )),
            (None, _, Linearity::NonLinear(_)) => Err(Error::Unsupported(
                "the default linear method cannot evaluate this product".into(),
            )),
        }
    }

    pub(crate) fn discount_curve(
        &self,
        ccy: Currency,
        collateral: &CollateralKey,
    ) -> Option<&CurveId> {
        self.0.discounting.iter().find_map(|rule| {
            (rule.currency == ccy && rule.collateral == *collateral).then_some(&rule.curve)
        })
    }

    pub(crate) fn forecast_curve(&self, index: &IndexId) -> Option<&CurveId> {
        self.0.forecasting.get(index)
    }

    pub(crate) fn fx_collateral(&self, pair: rhodium_core::fx::CcyPair) -> Option<&CurveId> {
        self.0
            .fx_collateral
            .get(&pair)
            .or_else(|| self.0.fx_collateral.get(&pair.inverse()))
    }
}
`````

`````
crates/rhodium-pricing/src/context.rs
//! 評価時点と評価に使う不変の入力。

use rhodium_core::{fixing::FixingStore, refdata::RefData, time::Instant};

use crate::config::PricingConfig;

/// 評価に使うスナップショットと設定。
#[derive(Debug)]
pub struct Context<'a> {
    /// 評価の瞬間。
    pub as_of: Instant,
    /// 評価に使う参照データ。
    pub refdata: &'a RefData,
    /// 評価時点で利用できる fixing。
    pub fixings: &'a FixingStore,
    /// 方式、曲線および近似の設定。
    pub config: &'a PricingConfig,
}
`````

`````
crates/rhodium-pricing/src/definition.rs
//! 市場の構造と校正 quote の定義。

mod fx;
pub(crate) use fx::ResolvedFxDefinition;
pub use fx::{DirectFx, FxDatedQuote, FxDefinition, FxImpliedCurve, FxImpliedLongEnd, FxQuotes};

use std::collections::{BTreeMap, BTreeSet};

use rhodium_core::{
    Error, Result,
    currency::Currency,
    id::CalendarId,
    quantity::YearFraction,
    time::{Instant, TzId},
};
use rhodium_market::{
    curve::{CurveId, CurveRole, CurveSpec},
    quote::{QuoteId, QuoteKind},
    structure::{BlockId, ParamId},
};
use rhodium_products::standard::InstrumentTemplate;

use crate::{config::PricingConfig, resolve::CollateralKey};

/// 市場から測る校正目標。
#[derive(Debug, Clone)]
pub enum CalibrationTarget {
    /// 取引の PV がゼロになる固定金利。
    ParRate(InstrumentTemplate),
    /// 先物の価格または FX のアウトライト。
    Price(InstrumentTemplate),
    /// quote を変換 q̃ してそのまま θ の一要素にする（マーク、年末のターン）。商品を持たない。
    DirectMark,
}

impl CalibrationTarget {
    pub(crate) const fn template(&self) -> Option<&InstrumentTemplate> {
        match self {
            Self::ParRate(template) | Self::Price(template) => Some(template),
            Self::DirectMark => None,
        }
    }
}

/// 入力 quote を目標の単位に変換する規則 q̃。
#[derive(Debug, Clone, Default, PartialEq, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum QuoteTransform {
    /// 入力値をそのまま使う。
    #[default]
    Identity,
    /// 期末のスプレッド s に受渡期間の日数係数 `τ_turn` を掛け、ln DF の跳び `s·τ_turn` にする。
    ///
    /// `τ_turn` は期末の営業日から翌営業日までの日数係数で、M1 では呼び出し側（ファイルでは
    /// `{ turn_spread_to_jump = { value = 0.0164, basis = "ACT365F" } }`）が明示する。
    /// 参照データの期末の規則から解決する書き方は、ファイル形式の見直し（M6）で足す。
    TurnSpreadToJump(YearFraction),
}

impl QuoteTransform {
    /// 入力 quote q に掛けて目標の値 q̃ = D·q にする係数 D。どの変換も quote について線形。
    pub(crate) fn scale(&self) -> f64 {
        match self {
            Self::Identity => 1.0,
            Self::TurnSpreadToJump(fraction) => fraction.value(),
        }
    }
}

/// 残差が当てはめる市場パラメータの範囲。
#[derive(Debug, Clone)]
pub enum ElementSelector {
    /// ブロックの全要素。
    Block(BlockId),
    /// ブロック内の一要素。
    Param(ParamId),
}

impl ElementSelector {
    pub(crate) const fn block(&self) -> &BlockId {
        match self {
            Self::Block(block) => block,
            Self::Param(param) => param.block(),
        }
    }
}

/// 一つの quote を市場のパラメータに結ぶ校正条件。
#[derive(Debug, Clone)]
pub struct InstrumentBinding {
    pub(crate) quote: QuoteId,
    pub(crate) target: CalibrationTarget,
    pub(crate) collateral: CollateralKey,
    pub(crate) determines: ElementSelector,
    kind: QuoteKind,
    transform: QuoteTransform,
    tolerance: Option<f64>,
}

impl InstrumentBinding {
    /// 校正目標、quote の変換と当てはめの対象を指定する。
    ///
    /// OIS のパーレートは `Rate` の quote を変換せずにブロックか要素へ当てはめる。
    /// 直接マークは θ の一要素だけを決め、年末のターンは `Spread` の quote を
    /// [`QuoteTransform::TurnSpreadToJump`] で跳びの大きさに変える。
    ///
    /// # Errors
    ///
    /// 目標と quote の種類・変換・当てはめの対象が合わないか、OIS のテナーがゼロか、
    /// ターンの日数係数が正でない場合。
    pub fn new(
        quote: QuoteId,
        kind: QuoteKind,
        target: CalibrationTarget,
        transform: QuoteTransform,
        collateral: CollateralKey,
        determines: ElementSelector,
    ) -> Result<Self> {
        let valid = match (&target, kind, &transform, &determines) {
            (
                CalibrationTarget::ParRate(InstrumentTemplate::Ois { tenor, .. }),
                QuoteKind::Rate,
                QuoteTransform::Identity,
                _,
            ) => tenor.amount() > 0,
            (
                CalibrationTarget::ParRate(
                    InstrumentTemplate::Deposit { .. }
                    | InstrumentTemplate::Fra { .. }
                    | InstrumentTemplate::Swap { .. },
                ),
                QuoteKind::Rate,
                QuoteTransform::Identity,
                _,
            )
            | (
                CalibrationTarget::ParRate(
                    InstrumentTemplate::BasisSwap { .. } | InstrumentTemplate::XccyBasis { .. },
                ),
                QuoteKind::Spread,
                QuoteTransform::Identity,
                _,
            )
            | (
                CalibrationTarget::Price(InstrumentTemplate::RateFuture { .. }),
                QuoteKind::FuturesPrice,
                QuoteTransform::Identity,
                _,
            )
            | (
                CalibrationTarget::Price(InstrumentTemplate::FxForward { .. }),
                QuoteKind::Price,
                QuoteTransform::Identity,
                _,
            )
            | (
                CalibrationTarget::DirectMark,
                _,
                QuoteTransform::Identity,
                ElementSelector::Param(_),
            ) => true,
            // 跳びは期間の長さを掛けたスプレッドなので、日数係数が正でなければ意味を持たない。
            (
                CalibrationTarget::DirectMark,
                QuoteKind::Spread,
                QuoteTransform::TurnSpreadToJump(fraction),
                ElementSelector::Param(_),
            ) => fraction.value() > 0.0,
            _ => false,
        };
        if valid {
            Ok(Self {
                quote,
                target,
                collateral,
                determines,
                kind,
                transform,
                tolerance: None,
            })
        } else {
            Err(Error::InvalidInput(format!(
                "invalid calibration target, quote kind or parameter selection for {quote}"
            )))
        }
    }

    /// 入力 quote の ID。
    #[must_use]
    pub const fn quote(&self) -> &QuoteId {
        &self.quote
    }

    /// 入力 quote の種類。
    #[must_use]
    pub const fn kind(&self) -> QuoteKind {
        self.kind
    }

    /// 入力 quote の変換。
    #[must_use]
    pub const fn transform(&self) -> &QuoteTransform {
        &self.transform
    }

    /// 最小二乗の当てはめの許容幅 σ を quote の単位で与える。
    ///
    /// 許容幅を持つ binding だけからなるグループは、重み W = diag(1/σ̃²) の最小二乗で解く
    /// （σ̃ = σ·D は残差の単位の幅で、D は quote の変換の係数）。許容幅を持たない binding の
    /// グループは正方系で解き、未知数と残差の数の一致を求める。
    ///
    /// # Errors
    ///
    /// 許容幅が有限な正の数でなければ `InvalidInput`。
    pub fn with_tolerance(mut self, tolerance: f64) -> Result<Self> {
        if !(tolerance.is_finite() && tolerance > 0.0) {
            return Err(Error::InvalidInput(format!(
                "calibration tolerance of {} must be finite and positive",
                self.quote
            )));
        }
        self.tolerance = Some(tolerance);
        Ok(self)
    }

    /// 最小二乗の当てはめの許容幅 σ（quote の単位）。正方系の binding は `None`。
    #[must_use]
    pub const fn tolerance(&self) -> Option<f64> {
        self.tolerance
    }
}

/// 会合日の階段（`Steps` の部品）の隣り合う区間のフォワードの差に掛ける罰則。
///
/// 罰則 `P(θ) = (λ/2) Σ_i ((f_{i+1} − f_i)/s)²` を最小二乗の目的に足す（λ は重み、s は参照の幅）。
/// 定型のテナーの OIS だけで階段を当てると区間ごとの水準が決まらないので、
/// 罰則が隣り合う区間を近づけて解を一意にする。罰則だけで決まる水準は政策金利の織り込みとして報告しない。
#[derive(Debug, Clone, PartialEq)]
pub struct StepPenalty {
    pub(crate) block: BlockId,
    pub(crate) weight: f64,
    pub(crate) scale: f64,
}

impl StepPenalty {
    /// 罰則を掛けるブロック、重み λ と参照の幅 s を指定する。
    ///
    /// # Errors
    ///
    /// λ か s が有限な正の数でなければ `InvalidInput`。
    pub fn new(block: BlockId, weight: f64, scale: f64) -> Result<Self> {
        if !([weight, scale]
            .iter()
            .all(|value| value.is_finite() && *value > 0.0))
        {
            return Err(Error::InvalidInput(format!(
                "step penalty of {} {} needs a finite positive weight and scale",
                block.curve(),
                block.part()
            )));
        }
        Ok(Self {
            block,
            weight,
            scale,
        })
    }

    /// 罰則の行 `√λ·(f_{i+1} − f_i)/s` の係数 `√λ/s`。
    pub(crate) fn factor(&self) -> f64 {
        self.weight.sqrt() / self.scale
    }
}

/// 曲線の構造と校正条件。参照データと評価時点は校正時に受け取る。
#[derive(Debug, Clone)]
pub struct MarketDefinition {
    pub(crate) valuation_calendar: CalendarId,
    pub(crate) valuation_zone: TzId,
    pub(crate) curves: Vec<CurveSpec>,
    pub(crate) bindings: Vec<InstrumentBinding>,
    pub(crate) fx: Vec<FxDefinition>,
    pub(crate) lgm: BTreeMap<Currency, LgmDefinition>,
    pub(crate) penalties: BTreeMap<BlockId, StepPenalty>,
    pub(crate) gauss_newton: BTreeSet<BlockId>,
    quote_kinds: BTreeMap<QuoteId, QuoteKind>,
}

impl MarketDefinition {
    /// 予想曲線ごとに、その予想を定める担保の割引曲線（予想の測度）を求める。
    ///
    /// 予想曲線のブロックを当てる商品の binding がすべて、価格設定の割引の表で
    /// 曲線の通貨と binding の担保から同じ割引曲線に解決されるとき、その曲線が測度になる。
    /// 担保の違う binding が混ざる、割引の規則が無い、商品で当てる binding が無い
    /// （直接マークだけ）予想曲線は測度が分からないので返さない。
    pub(crate) fn forecast_measures(&self, config: &PricingConfig) -> Vec<(CurveId, CurveId)> {
        self.curves
            .iter()
            .filter(|curve| curve.has_role(CurveRole::Projection))
            .filter_map(|curve| {
                let mut measures = self
                    .bindings
                    .iter()
                    .filter(|binding| {
                        binding.target.template().is_some()
                            && binding.determines.block().curve() == curve.id()
                    })
                    .map(|binding| config.discount_curve(curve.currency(), &binding.collateral));
                let first = measures.next()??;
                measures
                    .all(|measure| measure == Some(first))
                    .then(|| (curve.id().clone(), first.clone()))
            })
            .collect()
    }
}

/// 曲線のブロックのすべての要素が、当てはめか直接マークに割り当てられているかを確かめる。
fn check_assigned(assignment: &BTreeMap<BlockId, Vec<Option<Assigned>>>) -> Result<()> {
    if let Some((block, element)) = assignment.iter().find_map(|(block, flags)| {
        flags
            .iter()
            .position(Option::is_none)
            .map(|element| (block, element))
    }) {
        return Err(Error::InvalidInput(format!(
            "unassigned element {element} of curve {} part {}",
            block.curve(),
            block.part()
        )));
    }
    Ok(())
}

/// θ の要素への割り当て。直接マークと商品の当てはめの競合を検出する。
#[derive(Debug, Clone, Copy)]
enum Assigned {
    /// 商品の残差が当てはめる。
    Fit,
    /// 直接マークが値を決める。`fitted` は最小二乗の許容幅を持つか。
    Direct { fitted: bool },
}

/// 通貨別の区間一定 LGM に渡す市場マーク。
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LgmDefinition {
    /// 区間を区切る瞬間。最後のマークは末尾の区間で使う。
    #[serde(default)]
    pub breaks: Vec<Instant>,
    /// 区間順の平均回帰のマーク。
    pub mean_reversion: Vec<QuoteId>,
    /// 区間順の LGM ボラのマーク。
    pub volatility: Vec<QuoteId>,
}

impl MarketDefinition {
    /// 明示した UTC ピラーを持つ曲線と、その全節点を決める binding を受け取る。
    ///
    /// # Errors
    ///
    /// 曲線 ID が重複するか、binding が未知の要素を指すか、未割当の要素がある場合。
