            Some(monthly),
        )?))
    };
    let one_coupon = explicit(vec![period(
        "2024-12-02",
        "2025-01-15",
        vec![regular, virtual_regular],
    )?])?;
    let split = explicit(vec![
        period("2024-12-02", "2025-01-02", Vec::new())?,
        period("2025-01-02", "2025-01-15", vec![virtual_regular])?,
    ])?;
    let pv = |schedule: &ScheduleSpec| -> Result<f64> {
        let product = swap(
            Money::new(Currency::JPY, 100_000_000.0),
            schedule,
            &DayCount::Act365F,
            INDEX,
            0.0,
            0.0,
            plain()?,
        )?;
        let prepared = prepare_standard(&product, "2024-11-28".parse()?, &structure, &context)?;
        Ok(evaluate(&prepared, &market)?.get(Currency::JPY))
    };
    let generated = pv(&generated)?;
    assert_close!(pv(&one_coupon)?, generated, REL, JPY_ABS);
    // 分けた二つの期間は別のクーポンで、差は交差項の数百円になる。
    let difference = generated - pv(&split)?;
    assert!(difference > 100.0, "{difference}");
    Ok(())
}

#[test]
fn intraday_forward_fraction_follows_dst_day_length() -> Result<()> {
    // 評価時点をまたぐ金曜の SOFR は、残りの時間を前後の観測時刻の UTC の差で割って日数係数を補間する。
    // 夏時間の切り替え日は 23 時間か 25 時間なので、24 時間を仮定すると残りの割合が 12/24 や
    // 13/24 になり、PV が数 % ずれて許容誤差を外れる。
    let edt = UtcOffset::from_seconds(-4 * 3600)?;
    let est = UtcOffset::from_seconds(-5 * 3600)?;
    let new_york = ZoneRules::new(
        "2024-10-01T00:00:00Z".parse()?..="2025-05-02T00:00:00Z".parse()?,
        edt,
        vec![
            ("2024-11-03T06:00:00Z".parse()?, est),
            ("2025-03-09T07:00:00Z".parse()?, edt),
        ],
    )?;
    let index = "USD-SOFR";
    let data = refdata_in_zone(
        Currency::USD,
        index,
        DayCount::Act360,
        Vec::new(),
        "America/New_York",
        new_york,
    )?;
    // またぐ要素の因子 Π_A は厳密な式なので、日次積の近似を禁止しても準備できる。
    let config = config(
        Currency::USD,
        index,
        true,
        "exclude",
        DISABLE_DAILY_FORWARD_PRODUCT,
    )?;
    // USD の最小単位 0.01 の 1e-6 倍。
    let usd_abs = 1e-8;
    let theta = -0.004;
    for (trade, start, end, as_of, period_end, pillar, remaining, day_length) in [
        // 2025-03-09（日）は 00:00 EST = 05:00Z から翌日 00:00 EDT = 04:00Z までの 23 時間。
        // 12:00 EDT の評価時点から月曜の観測時刻まで残り 12 時間。
        (
            "2025-03-06",
            "2025-03-07",
            "2025-03-10",
            "2025-03-09T16:00:00Z",
            "2025-03-10T04:00:00Z",
            "2025-03-24T00:00:00Z",
            12.0,
            23.0,
        ),
        // 2024-11-03（日）は 00:00 EDT = 04:00Z から翌日 00:00 EST = 05:00Z までの 25 時間。
        // 11:00 EST の評価時点から月曜の観測時刻まで残り 13 時間。
        (
            "2024-10-31",
            "2024-11-01",
            "2024-11-04",
            "2024-11-03T16:00:00Z",
            "2024-11-04T05:00:00Z",
            "2024-11-18T00:00:00Z",
            13.0,
            25.0,
        ),
    ] {
        let as_of: Instant = as_of.parse()?;
        let pillar: Instant = pillar.parse()?;
        let structure = structure(Currency::USD, as_of, &[pillar], true)?;
        let ctx = Context {
            as_of,
            refdata: &data,
            fixings: &FixingStore::default(),
            config: &config,
        };
        let dates = [start.parse()?, end.parse()?];
        let product = swap(
            Money::new(Currency::USD, NOTIONAL),
            &schedule(&dates, &dates[1..])?,
            &DayCount::Act360,
            index,
            0.0,
            0.0,
            plain()?,
        )?;
        let prepared = prepare_standard(&product, trade.parse()?, &structure, &ctx)?;
        assert!(
            prepared
                .warnings()
                .iter()
                .any(|warning| matches!(warning, Warning::UnpublishedFixingForecast { .. }))
        );

        // 金曜から月曜の元の重み 3/360 に、評価時点から月曜の観測時刻までの forward を掛ける。
        // 単一曲線で支払いは月曜の決済時刻なので、変動脚の PV は N·(3/360)·r·P(e)。
        let discount = linear_df(as_of, pillar, theta, period_end.parse()?);
        let tau = remaining / day_length / 360.0;
        let rate = (discount.recip() - 1.0) / tau;
        let expected = NOTIONAL * 3.0 / 360.0 * rate * discount;
        let market = Market::from_theta(structure, vec![theta])?;
        let pv = evaluate(&prepared, &market)?.get(Currency::USD);
        assert_close!(pv, expected, REL, usd_abs);
    }
    Ok(())
}

#[test]
fn fixing_publication_boundary_uses_estimate_before_and_published_at_or_after() -> Result<()> {
    // 公表直前は期間の途中なので推定値を使う。ちょうどの時点から実績値を要求する。
    let data = jpy()?;
    let store = fixings(&data, &[("2025-01-09", 0.07)], &[("2025-01-09", 0.02)])?;
    let only_estimate = fixings(&data, &[], &[("2025-01-09", 0.02)])?;
    let config = config(Currency::JPY, INDEX, true, "exclude", "")?;
    let product = jpy_swap(
        &["2025-01-09", "2025-01-10"],
        &["2025-01-13"],
        0.0,
        0.0,
        plain()?,
    )?;
    let date = "2025-01-08".parse()?;
    for (instant, rate) in [
        ("2025-01-09T23:59:59.999999Z", 0.02),
        ("2025-01-10T00:00:00Z", 0.07),
        ("2025-01-10T00:00:00.000001Z", 0.07),
    ] {
        let as_of = instant.parse()?;
        let structure = structure(
            Currency::JPY,
            as_of,
            &["2025-02-10T00:00:00Z".parse()?],
            true,
        )?;
        let context = Context {
            as_of,
            refdata: &data,
            fixings: &store,
            config: &config,
        };
        let prepared = prepare_standard(&product, date, &structure, &context)?;
        let pv = evaluate(
            &prepared,
            &Market::from_theta(structure.clone(), vec![0.0])?,
        )?
        .get(Currency::JPY);
        assert_close!(pv, NOTIONAL * rate / 365.0, REL, JPY_ABS);
        if rate > 0.02 {
            let missing = Context {
                fixings: &only_estimate,
                ..context
            };
            assert!(matches!(
                prepare_standard(&product, date, &structure, &missing),
                Err(Error::MissingData(_))
            ));
        }
    }
    Ok(())
}

#[test]
fn unpublished_term_fixing_uses_the_forward_over_the_index_deposit_period() -> Result<()> {
    // 元のターム金利が有効な期間の手計算。指標の預金開始は fixing の翌営業日。
    let index: IndexId = "JPY-TERM".parse()?;
    let ccy = Currency::JPY;
    let mut parts = refdata_parts(
        ccy,
        "JPY-TERM",
        DayCount::Act365F,
        Vec::new(),
        "Etc/UTC",
        ZoneRules::new(
            "2025-01-01T00:00:00Z".parse()?..="2025-04-30T23:59:59Z".parse()?,
            UtcOffset::from_seconds(0)?,
            Vec::new(),
        )?,
    )?;
    let definition = parts
        .indices
        .get_mut(&index)
        .ok_or_else(|| Error::Internal("test index".into()))?;
    definition.kind = IndexKind::TermRate {
        ccy,
        day_count: DayCount::Act365F,
        tenor: "2BD".parse()?,
        value_offset: "1BD".parse()?,
        bdc: BusinessDayConvention::Following,
        end_of_month: EndOfMonthRule::Ignore,
    };
    let refdata = RefData::new(parts)?;
    let date = ObservationPoint::Date("2025-01-07".parse()?);
    let as_of = "2025-01-06T00:00:00Z".parse()?;
    let pay = "2025-01-10T00:00:00Z".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let rate = b.fixing(&index, date);
            let notional_tau = b.constant(NOTIONAL * 2.0 / 365.0);
            let amount = b.mul([notional_tau, rate]);
            Ok(b.pay(Asset::Cash(ccy), amount, pay))
        })?,
        &refdata,
    )?;
    let structure = structure(ccy, as_of, &[pay], true)?;
    let ctx = Context {
        as_of,
        refdata: &refdata,
        fixings: &FixingStore::default(),
        config: &config(ccy, index.as_str(), true, "exclude", "")?,
    };
    let prepared = prepare(
        PrepareInput::Ir(&product),
        &CollateralKey::Uncollateralized,
        &[],
        &structure,
        &ctx,
    )?;
    let market = Market::from_theta(structure, vec![-0.04])?;
    let pv = evaluate(&prepared, &market)?.get(ccy);
    // 期間は水曜から金曜。金曜の ln DF が -0.04、水曜は評価と金曜の中間で -0.02。
    assert_close!(
        pv,
        NOTIONAL * (-0.04_f64).exp() * 0.02_f64.exp_m1(),
        REL,
        JPY_ABS
    );
    Ok(())
}

#[test]
fn fallback_pv_compounds_a_published_day_with_a_projected_day() -> Result<()> {
    // 最初の TONA は実績、後半 1 日は予想。評価時点の曲線で後半の成長と受渡の割引を求める。
    let ccy = Currency::JPY;
    let term: IndexId = "JPY-TERM".parse()?;
    let overnight: IndexId = INDEX.parse()?;
    let reset = "2025-01-08".parse()?;
    let basis = DayCount::Act365F;
    let mut parts = refdata_parts(
        ccy,
        INDEX,
        basis.clone(),
        Vec::new(),
        "Etc/UTC",
        ZoneRules::new(
            "2025-01-01T00:00:00Z".parse()?..="2025-04-30T23:59:59Z".parse()?,
            UtcOffset::from_seconds(0)?,
            Vec::new(),
        )?,
    )?;
    let mut definition = parts.indices[&overnight].clone();
    definition.id = term.clone();
    definition.status = IndexStatus::Ceased { since: reset };
    definition.kind = IndexKind::TermRate {
        ccy,
        day_count: basis.clone(),
        tenor: "2BD".parse()?,
        value_offset: "0BD".parse()?,
        bdc: BusinessDayConvention::Following,
        end_of_month: EndOfMonthRule::Ignore,
    };
    parts.indices.insert(term.clone(), definition);
    parts.fallbacks.insert(
        term.clone(),
        FallbackDefinition {
            index: term.clone(),
            trigger: FallbackTrigger::Cessation,
            replacement: FallbackReplacement {
                index: overnight.clone(),
                observation: plain()?,
                annualization: basis,
                spread: "10bp".parse()?,
            },
            applies: Vec::new(),
        },
    );
    let refdata = RefData::new(parts)?;
    let as_of = "2025-01-09T00:00:00Z".parse()?;
    let pay = "2025-01-10T00:00:00Z".parse()?;
    let product = Contract::build(|b| {
        let rate = b.fixing(&term, ObservationPoint::Date(reset));
        let c = b.constant(NOTIONAL * 2.0 / 365.0);
        let amount = b.mul([c, rate]);
        Ok(b.pay(Asset::Cash(ccy), amount, pay))
    })?;
    let trade = Trade {
        id: "fallback".parse()?,
        margin_set: "margin".parse()?,
        book: "book".parse()?,
        trade_date: "2025-01-07".parse()?,
        product: ProductSpec::Contract(product),
        attributes: std::collections::BTreeMap::new(),
        exercises: Vec::new(),
        fees: Vec::new(),
        fallback: FallbackTerms {
            by_index: [(term, FallbackSelection::IndexDefault)].into(),
            ..FallbackTerms::default()
        },
    };
    let product = trade.to_product(&refdata)?;
    let structure = structure(ccy, as_of, &[pay], true)?;
    let fixed = FixingStore::new(
        &refdata,
        FixingStoreParts {
            published: vec![(overnight, ObservationPoint::Date(reset), 0.073)],
            ..Default::default()
        },
    )?;
    let ctx = Context {
        as_of,
        refdata: &refdata,
        fixings: &fixed,
        config: &config(ccy, INDEX, true, "exclude", "")?,
    };
    let prepared = prepare(
        PrepareInput::Ir(&product),
        &CollateralKey::Uncollateralized,
        &[],
        &structure,
        &ctx,
    )?;
    let market = Market::from_theta(structure, vec![-0.0002])?;
    let first = 0.073 / 365.0;
    let second = 0.0002_f64.exp_m1();
    let expected =
        NOTIONAL * (-0.0002_f64).exp() * (first + second + first * second + 0.001 * 2.0 / 365.0);
    assert_close!(
        evaluate(&prepared, &market)?.get(ccy),
        expected,
        REL,
        JPY_ABS
    );
    Ok(())
}

#[test]
fn settlement_boundary_and_binding_failures() -> Result<()> {
    // 評価時点と同時刻の受け渡しの DF は 1。別の構造と未校正の必須ブロックを区別する。
    let as_of = "2025-01-07T00:00:00Z".parse()?;
    let pillar = "2025-02-10T00:00:00Z".parse()?;
    let market_structure = structure(Currency::JPY, as_of, &[pillar], true)?;
    let other = structure(Currency::JPY, as_of, &[pillar], true)?;
    let data = jpy()?;
    let product = Product::new(
        Contract::build(|builder| {
            let amount = builder.constant(100.0);
            Ok(builder.pay(Asset::Cash(Currency::JPY), amount, as_of))
        })?,
        &data,
    )?;
    let store = FixingStore::default();
    for (setting, expected, currencies) in [("include", 100.0, 1), ("exclude", 0.0, 0)] {
        let config = config(Currency::JPY, INDEX, true, setting, "")?;
        let ctx = Context {
            as_of,
            config: &config,
            refdata: &data,
            fixings: &store,
        };
        let prepared = prepare(
            PrepareInput::Ir(&product),
            &CollateralKey::Uncollateralized,
            &[],
            &market_structure,
            &ctx,
        )?;
        let market = Market::from_theta(market_structure.clone(), vec![-0.02])?;
        let pv = market.with_view(|view| prepared.bind(&view)?.evaluate())?;
        assert_close!(pv.get(Currency::JPY), expected, REL, JPY_ABS);
        assert_eq!(pv.terms().count(), currencies);
        assert_eq!(prepared.blocks(), []);

        assert!(matches!(
            evaluate(&prepared, &Market::from_theta(other.clone(), vec![-0.02])?),
            Err(Error::InvalidInput(_))
        ));
        let partial = Market::from_partial_theta(market_structure.clone(), vec![f64::NAN], [])?;
        if currencies == 1 {
            assert!(matches!(
                evaluate(&prepared, &partial),
                Err(Error::MissingData(_))
            ));
        } else {
            assert_eq!(evaluate(&prepared, &partial)?.terms().count(), 0);
        }
    }
    Ok(())
}

#[test]
fn prepared_swap_reuses_structure_and_preserves_aad_gradient()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // 単一曲線では変動脚の PV は N(P(開始) - P(終了))。割引と予想の感応度が両方流れる。
    let as_of = "2025-01-07T00:00:00Z".parse()?;
    let start = "2025-01-08T00:00:00Z".parse()?;
    let end = "2025-01-14T00:00:00Z".parse()?;
    let structure = structure(Currency::JPY, as_of, &[end], true)?;
    let data = jpy()?;
    let config = config(
        Currency::JPY,
        INDEX,
        true,
        "exclude",
        DISABLE_DAILY_FORWARD_PRODUCT,
    )?;
    let context = Context {
        as_of,
        refdata: &data,
        fixings: &FixingStore::default(),
        config: &config,
    };
    let product = jpy_swap(
        &["2025-01-08", "2025-01-14"],
        &["2025-01-14"],
        0.02,
        0.0,
        plain()?,
    )?;
    let prepared = prepare_standard(&product, "2025-01-07".parse()?, &structure, &context)?;
    assert_eq!(prepared.blocks(), &[BlockId::nodes("discount".parse()?)]);
    let factor = 1.0 + 0.02 * 6.0 / 365.0;
    let expected =
        |theta| NOTIONAL * (linear_df(as_of, end, theta, start) - factor * f64::exp(theta));
    let market = Market::from_theta(structure.clone(), vec![-0.002])?;
    assert_close!(
        evaluate(&prepared, &market)?.get(Currency::JPY),
        expected(-0.002),
        REL,
        JPY_ABS
    );

    let theta: f64 = -0.004;
    let (pv, gradient) = jpy_pv_and_gradient(&structure, &prepared, theta)?;
    assert_close!(pv, expected(theta), REL, JPY_ABS);
    let derivative = NOTIONAL * (linear_df(as_of, end, theta, start) / 7.0 - factor * theta.exp());
    assert_close!(gradient, derivative, 1e-10, JPY_ABS);
    Ok(())
}

#[test]
fn lookback_daily_product_preserves_weights_and_checks_policy() -> Result<()> {
    // 1 営業日の lookback では、木曜の率を金曜から月曜の 3 日に、金曜の率を月曜の 1 日に掛ける。
    let as_of = "2025-01-07T00:00:00Z".parse()?;
    let pillars = [
        "2025-01-10T00:00:00Z".parse()?,
        "2025-01-14T00:00:00Z".parse()?,
    ];
    let structure = structure(Currency::JPY, as_of, &pillars, true)?;
    let date = "2025-01-07".parse()?;
    let data = jpy()?;
    let lookback = rule(
        ObservationMethod::Lookback(NonZeroU8::MIN),
        CompoundingForm::SimpleFactor,
        SpreadTreatment::Exclusive,
        ObservationRounding::default(),
    )?;
    let product = jpy_swap(
        &["2025-01-10", "2025-01-14"],
        &["2025-01-14"],
        0.0,
        0.0,
        lookback,
    )?;
    let defaults = config(Currency::JPY, INDEX, true, "exclude", "")?;
    let ctx = Context {
        as_of,
        refdata: &data,
        fixings: &FixingStore::default(),
        config: &defaults,
    };
    let prepared = prepare_standard(&product, date, &structure, &ctx)?;
    assert!(prepared.warnings().contains(&Warning::DailyForwardProduct));
    let market = Market::from_theta(structure.clone(), vec![-0.0003, -0.0015])?;
    let thursday = (0.0001_f64).exp_m1();
    let friday = (0.0009_f64).exp_m1() / 3.0;
    let expected = NOTIONAL * ((1.0 + 3.0 * thursday) * (1.0 + friday) - 1.0) * (-0.0015_f64).exp();
    assert_close!(
        evaluate(&prepared, &market)?.get(Currency::JPY),
        expected,
        REL,
        JPY_ABS
    );

    // 近似の禁止、lookback の上限、支払日のずれの上限をそれぞれ準備時に適用する。
    for settings in [
        DISABLE_DAILY_FORWARD_PRODUCT,
        "max_lookback_lockout_business_days = 0",
        "max_payment_delay_business_days = 0",
    ] {
        let policy = config(Currency::JPY, INDEX, true, "exclude", settings)?;
        let ctx = Context {
            config: &policy,
            ..ctx
        };
        assert!(matches!(
            prepare_standard(&product, date, &structure, &ctx),
            Err(Error::Unsupported(_))
        ));
    }
    Ok(())
}

#[test]
fn plain_holiday_start_and_lockout_preserve_daily_weights() -> Result<()> {
    // 非営業日に始まる Plain と Lockout は金利の期間と重みの期間がずれるため、端点の比に畳めない。
    // 日次因子の積から求めた PV が手計算と一致し、日次積の禁止と lookback・lockout の営業日数の上限が
    // 準備のときに効くことを確かめる。
    let date = "2025-01-07".parse()?;
    let as_of = at(date)?;
    let pillars = [
        "2025-01-10T00:00:00Z".parse()?,
        "2025-01-13T00:00:00Z".parse()?,
        "2025-01-14T00:00:00Z".parse()?,
    ];
    let structure = structure(Currency::JPY, as_of, &pillars, true)?;
    let data = jpy()?;
    let defaults = config(Currency::JPY, INDEX, true, "exclude", "")?;
    let disabled = config(
        Currency::JPY,
        INDEX,
        true,
        "exclude",
        DISABLE_DAILY_FORWARD_PRODUCT,
    )?;
    let zero_limit = config(
        Currency::JPY,
        INDEX,
        true,
        "exclude",
        "max_lookback_lockout_business_days = 0",
    )?;
    let context = Context {
        as_of,
        refdata: &data,
        fixings: &FixingStore::default(),
        config: &defaults,
    };
    // 金曜・月曜・火曜の ln DF。評価時点からの線形補間で、日ごとの DF 比は
    // 木曜 exp(0.0001)、金曜から月曜 exp(0.0009)、月曜 exp(0.0012) になる。
    let pay_ln_df: f64 = -0.0024;
    let market = Market::from_theta(structure.clone(), vec![-0.0003, -0.0012, pay_ln_df])?;
    let thursday_excess = 0.0001_f64.exp_m1();
    let friday_daily_excess = 0.0009_f64.exp_m1() / 3.0;
    let monday_excess = 0.0012_f64.exp_m1();
    // 4 番目の要素は、営業日数の上限 0 で準備を拒否するか。上限は lookback と lockout にだけ掛かる。
    for (method, start, growth, exceeds_limit) in [
        // 土曜始まりでは、金曜の 3 日物の率を土日 2 日の重みで使う。
        (
            ObservationMethod::Plain,
            "2025-01-11",
            (1.0 + 2.0 * friday_daily_excess) * (1.0 + monday_excess),
            false,
        ),
        // 木曜始まりの Lockout(1) では、金曜の率を金曜からの 3 日分と月曜の 1 日分に使う。
        (
            ObservationMethod::Lockout(NonZeroU8::MIN),
            "2025-01-09",
            (1.0 + thursday_excess)
                * (1.0 + 3.0 * friday_daily_excess)
                * (1.0 + friday_daily_excess),
            true,
        ),
    ] {
        let rule = rule(
            method,
            CompoundingForm::SimpleFactor,
            SpreadTreatment::Exclusive,
            ObservationRounding::default(),
        )?;
        let product = jpy_swap(&[start, "2025-01-14"], &["2025-01-14"], 0.0, 0.0, rule)?;
        let prepared = prepare_standard(&product, date, &structure, &context)?;
        assert!(prepared.warnings().contains(&Warning::DailyForwardProduct));
        let expected = NOTIONAL * (growth - 1.0) * pay_ln_df.exp();
        assert_close!(
            evaluate(&prepared, &market)?.get(Currency::JPY),
            expected,
            REL,
            JPY_ABS
        );

        // 端点の比に畳めないので、日次積を禁止すると準備で拒否する。
        let disabled_context = Context {
            config: &disabled,
            ..context
        };
        assert!(matches!(
            prepare_standard(&product, date, &structure, &disabled_context),
            Err(Error::Unsupported(message)) if message.contains("DailyForwardProduct")
        ));
        // 上限 0 では Lockout(1) を拒否し、休日始まりの Plain は上限に関係なく日次積で準備できる。
        let limited_context = Context {
            config: &zero_limit,
            ..context
        };
        let limited = prepare_standard(&product, date, &structure, &limited_context);
        if exceeds_limit {
            assert!(matches!(
                limited,
                Err(Error::Unsupported(message)) if message.contains("lockout exceeds")
            ));
        } else {
            assert!(limited?.warnings().contains(&Warning::DailyForwardProduct));
        }
    }
    Ok(())
}

#[test]
fn observation_shift_telescopes_without_replacing_coupon_accrual() -> Result<()> {
    // 1 営業日のシフトで金曜から月曜の計算期間が木曜から金曜の観測期間になる。
    // 元の期間に払うため、端点比で得た 1 日分の成長率を 3 日分に換算する。
    let as_of = "2025-01-07T00:00:00Z".parse()?;
    let end = "2025-01-13T00:00:00Z".parse()?;
    let structure = structure(Currency::JPY, as_of, &[end], true)?;
    let data = jpy()?;
    let shift = rule(
        ObservationMethod::ObservationShift(NonZeroU8::MIN),
        CompoundingForm::SimpleFactor,
        SpreadTreatment::Exclusive,
        ObservationRounding::default(),
    )?;
    let product = jpy_swap(
        &["2025-01-10", "2025-01-13"],
        &["2025-01-13"],
        0.0,
        0.0,
        shift,
    )?;
    let config = config(
        Currency::JPY,
        INDEX,
        true,
        "exclude",
        DISABLE_DAILY_FORWARD_PRODUCT,
    )?;
    let ctx = Context {
        as_of,
        refdata: &data,
        config: &config,
        fixings: &FixingStore::default(),
    };
    let prepared = prepare_standard(&product, "2025-01-07".parse()?, &structure, &ctx)?;
    let theta: f64 = -0.0018;
    let expected = NOTIONAL * 3.0 * (-theta / 6.0).exp_m1() * theta.exp();
    let pv = evaluate(&prepared, &Market::from_theta(structure, vec![theta])?)?.get(Currency::JPY);
    assert_close!(pv, expected, REL, JPY_ABS);
    Ok(())
}

#[test]
fn published_observation_shift_pays_three_yen() -> Result<()> {
    // 木曜の fixing 10% は金曜の評価時点で公表済み。木曜から金曜の 1 日で年率化し、
    // 金曜から月曜の 3 日に利息を払う。
    let start = "2025-01-10".parse()?;
    let end = "2025-01-13".parse()?;
    let as_of = at(start)?;
    let structure = structure(Currency::JPY, as_of, &[at(end)?], true)?;
    let data = jpy()?;
    let fixings = fixings(&data, &[("2025-01-09", 0.10)], &[])?;
    let config = config(
        Currency::JPY,
        INDEX,
        true,
        "exclude",
        DISABLE_DAILY_FORWARD_PRODUCT,
    )?;
    let context = Context {
        as_of,
        refdata: &data,
        fixings: &fixings,
        config: &config,
    };
    let shift = rule(
        ObservationMethod::ObservationShift(NonZeroU8::MIN),
        CompoundingForm::SimpleFactor,
        SpreadTreatment::Exclusive,
        ObservationRounding::default(),
    )?;
    let product = swap(
        Money::new(Currency::JPY, 3650.0),
        &schedule(&[start, end], &[end])?,
        &DayCount::Act365F,
        INDEX,
        0.0,
        0.0,
        shift,
    )?;
    let prepared = prepare_standard(&product, "2025-01-07".parse()?, &structure, &context)?;
    // θ = 0 では割引因子が 1 なので、PV が支払額 3,650 × 10% × 3/365 に等しい。
    let market = Market::from_theta(structure, vec![0.0])?;
    assert_close!(
        evaluate(&prepared, &market)?.get(Currency::JPY),
        3.0,
        REL,
        JPY_ABS
    );
    Ok(())
}

#[test]
fn daily_rounding_inclusive_spread_and_power_use_daily_factors() -> Result<()> {
    // 日次丸め、Inclusive のスプレッド、PowerFactor では各因子が 1 + w·r の形にならず、
    // telescoping の条件 1 を満たさない。木曜と金曜に異なる金利を与えて、丸めの位置、スプレッドを
    // 足す位置、重みを指数にする形がそれぞれ日次因子の積に一致し、日次積を禁止すると拒否することを確かめる。
    let date = "2025-01-07".parse()?;
    let as_of = at(date)?;
    let pillars = [
        "2025-01-10T00:00:00Z".parse()?,
        "2025-01-13T00:00:00Z".parse()?,
    ];
    let structure = structure(Currency::JPY, as_of, &pillars, true)?;
    let data = jpy()?;
    let defaults = config(Currency::JPY, INDEX, true, "exclude", "")?;
    let disabled = config(
        Currency::JPY,
        INDEX,
        true,
        "exclude",
        DISABLE_DAILY_FORWARD_PRODUCT,
    )?;
    let context = Context {
        as_of,
        refdata: &data,
        fixings: &FixingStore::default(),
        config: &defaults,
    };
    // 木曜の 1 日物は 3.4%、金曜の 3 日物は 4.7%。それぞれの金利から市場の端点を置く。
    let friday_ln_df = -3.0 * (0.034_f64 / 365.0).ln_1p();
    let monday_ln_df = friday_ln_df - (0.047_f64 * 3.0 / 365.0).ln_1p();
    let market = Market::from_theta(structure.clone(), vec![friday_ln_df, monday_ln_df])?;
    let rounded = Some(Rounding::new(0.01, RoundingMode::HalfEven)?);
    for (form, spread, daily, margin, growth) in [
        // 日次丸めでは 3.4% → 3%、4.7% → 5% の後に因子を作る。
        (
            CompoundingForm::SimpleFactor,
            SpreadTreatment::Exclusive,
            rounded,
            0.0,
            (1.0 + 0.03 / 365.0) * (1.0 + 0.05 * 3.0 / 365.0),
        ),
        // Inclusive ではスプレッド 1% を足してから各因子を作る。
        (
            CompoundingForm::SimpleFactor,
            SpreadTreatment::Inclusive,
            None,
            0.01,
            (1.0 + 0.044 / 365.0) * (1.0 + 0.057 * 3.0 / 365.0),
        ),
        // PowerFactor では各観測の重みが指数になる。
        (
            CompoundingForm::PowerFactor,
            SpreadTreatment::Exclusive,
            None,
            0.0,
            1.034_f64.powf(1.0 / 365.0) * 1.047_f64.powf(3.0 / 365.0),
        ),
    ] {
        let daily_rule = rule(
            ObservationMethod::Plain,
            form,
            spread,
            ObservationRounding {
                daily,
                aggregated: None,
            },
        )?;
        let product = jpy_swap(
            &["2025-01-09", "2025-01-13"],
            &["2025-01-13"],
            0.0,
            margin,
            daily_rule,
        )?;
        let prepared = prepare_standard(&product, date, &structure, &context)?;
        assert!(prepared.warnings().contains(&Warning::DailyForwardProduct));
        assert_eq!(
            prepared.warnings().contains(&Warning::RoundingAtForward),
            daily.is_some()
        );
        let expected = NOTIONAL * (growth - 1.0) * monday_ln_df.exp();
        assert_close!(
            evaluate(&prepared, &market)?.get(Currency::JPY),
            expected,
            REL,
            JPY_ABS
        );

        // どの規約も端点の比に畳めないので、日次積を禁止すると準備で拒否する。
        let disabled_context = Context {
            config: &disabled,
            ..context
        };
        assert!(matches!(
            prepare_standard(&product, date, &structure, &disabled_context),
            Err(Error::Unsupported(message)) if message.contains("DailyForwardProduct")
        ));
    }
    Ok(())
}

#[test]
fn different_forward_and_weight_bases_require_daily_factors() -> Result<()> {
    // telescoping の条件 5: 重みの日数基準が予想の日数基準（指標の ACT/365F）と違うときは、
    // 係数の値が同じでも端点の比に畳まず、契約の基準の重みで日次因子の積を評価することを確かめる。
    let data = jpy()?;
    let date = "2025-01-07".parse()?;
    let as_of = at(date)?;
    let end = "2025-01-10".parse()?;
    let spec = schedule(&["2025-01-08".parse()?, end], &[end])?;
    let currency = data.currency(Currency::JPY)?;
    let (schedule, _) = generate(
        &spec,
        date,
        data.calendars(),
        &SettlementContext {
            zones: data.zones(),
            default_time: Some(&currency.settlement_time),
            calendar: Some(&calendar()?),
        },
    )?;
    let structure = structure(Currency::JPY, as_of, &[at(end)?], true)?;
    let defaults = config(Currency::JPY, INDEX, true, "exclude", "")?;
    let disabled = config(
        Currency::JPY,
        INDEX,
        true,
        "exclude",
        DISABLE_DAILY_FORWARD_PRODUCT,
    )?;
    let context = Context {
        as_of,
        refdata: &data,
        fixings: &FixingStore::default(),
        config: &defaults,
    };
    let theta = -0.0018_f64;
    let market = Market::from_theta(structure.clone(), vec![theta])?;
    let rule = plain()?;
    let index = INDEX.parse()?;
    // 標準脚は指標の日数基準を使うため、重みの基準だけを変えたクーポンを IR で組む。
    // 1 月の NL/365 は ACT/365F と係数が同じなので、基準自体の違いでも日次積になることを確かめる。
    for (basis, denominator) in [(DayCount::Act360, 360.0), (DayCount::Nl365, 365.0)] {
        let (observations, _) = observation_schedule(
            &schedule.periods()[0],
            &rule,
            &basis,
            &AccrualContext::default(),
            data.calendars(),
        )?;
        let product = Product::new(
            Contract::build(|builder| {
                let growth = builder.compounded(&index, &observations, &rule, Spread::new(0.0)?)?;
                let one = builder.constant(1.0);
                let excess = builder.sub(growth, one);
                let notional = builder.constant(NOTIONAL);
                let amount = builder.mul([notional, excess]);
                Ok(builder.pay(
                    Asset::Cash(Currency::JPY),
                    amount,
                    schedule.payments()[0].at(),
                ))
            })?,
            &data,
        )?;
        let prepare_ir = |context| {
            prepare(
                PrepareInput::Ir(&product),
                &CollateralKey::Uncollateralized,
                &[],
                &structure,
                context,
            )
        };
        let prepared = prepare_ir(&context)?;
        assert!(prepared.warnings().contains(&Warning::DailyForwardProduct));
        // 連続する 1 日物の forward はどちらも 365·(exp(0.0006) − 1)。重みには契約の分母を使う。
        let factor = 1.0 + 365.0 * (-theta / 3.0).exp_m1() / denominator;
        let expected = NOTIONAL * (factor.powi(2) - 1.0) * theta.exp();
        assert_close!(
            evaluate(&prepared, &market)?.get(Currency::JPY),
            expected,
            REL,
            JPY_ABS
        );
        // 基準が違えば NL/365 の場合も日次積が要るので、禁止すると準備で拒否する。
        assert!(matches!(
            prepare_ir(&Context { config: &disabled, ..context }),
            Err(Error::Unsupported(message)) if message.contains("DailyForwardProduct")
        ));
    }
    Ok(())
}

#[test]
fn rounded_forward_compounds_before_aggregated_rounding() -> Result<()> {
    // 予想 3.4% を日次で 3% に丸めた後で因子にする。
    // Inclusive の単利因子なら年率 4%、PowerFactor ならスプレッド加算後に 3.96% に丸まる。
    let as_of = "2025-01-07T00:00:00Z".parse()?;
    let end = "2025-01-09T00:00:00Z".parse()?;
    let structure = structure(Currency::JPY, as_of, &[end], true)?;
    let theta = -2.0 * (0.034_f64 / 365.0).ln_1p();
    let data = jpy()?;
    let date = "2025-01-07".parse()?;
    let defaults = config(Currency::JPY, INDEX, true, "exclude", "")?;
    let disabled = config(
        Currency::JPY,
        INDEX,
        true,
        "exclude",
        r#"disabled = ["rounding_at_forward"]"#,
    )?;
    let ctx = Context {
        as_of,
        refdata: &data,
        fixings: &FixingStore::default(),
        config: &defaults,
    };
    let market = Market::from_theta(structure.clone(), vec![theta])?;
    for (form, spread, expected_rate) in [
        (
            CompoundingForm::SimpleFactor,
            SpreadTreatment::Inclusive,
            0.04,
        ),
        (
            CompoundingForm::PowerFactor,
            SpreadTreatment::Exclusive,
            0.0396,
        ),
    ] {
        let rounding = ObservationRounding {
            daily: Some(Rounding::new(0.01, RoundingMode::HalfEven)?),
            aggregated: Some(Rounding::new(0.0001, RoundingMode::HalfEven)?),
        };
        let rule = rule(ObservationMethod::Plain, form, spread, rounding)?;
        let product = jpy_swap(
            &["2025-01-08", "2025-01-09"],
            &["2025-01-09"],
            0.0,
            0.01,
            rule,
        )?;
        let prepared = prepare_standard(&product, date, &structure, &ctx)?;
        assert!(prepared.warnings().contains(&Warning::RoundingAtForward));
        assert!(prepared.warnings().contains(&Warning::DailyForwardProduct));
        let expected = NOTIONAL * expected_rate / 365.0 * theta.exp();
        assert_close!(
            evaluate(&prepared, &market)?.get(Currency::JPY),
            expected,
            REL,
            JPY_ABS
        );
        let ctx = Context {
            config: &disabled,
            ..ctx
        };
        assert!(matches!(
            prepare_standard(&product, date, &structure, &ctx),
            Err(Error::Unsupported(_))
        ));
    }
    Ok(())
}

#[test]
fn fixing_in_unresolved_choice_is_required_only_when_selected() -> Result<()> {
    // 未行使の枝に残った公表済みの欠損を予想しない。行使が決まれば選んだ枝だけ評価する。
    use rhodium_products::{ir::Party, trade::ExerciseRecord};

    let as_of = "2025-01-10T00:00:00Z".parse()?;
    let exercise_at = "2025-01-13T00:00:00Z".parse()?;
    let pay_at = "2025-01-14T00:00:00Z".parse()?;
    let observed = ObservationPoint::Date("2025-01-08".parse()?);
    let index: IndexId = INDEX.parse()?;
    let data = jpy()?;
    let product = Product::new(
        Contract::build(|builder| {
            let asset = Asset::Cash(Currency::JPY);
            let fixing = builder.fixing(&index, observed);
            let projected = builder.pay(asset.clone(), fixing, pay_at);
            let constant = builder.constant(100.0);
            let fixed = builder.pay(asset, constant, pay_at);
            Ok(builder.choose(Party::Us, "switch", exercise_at, [projected, fixed]))
        })?,
        &data,
    )?;

    let config = config(Currency::JPY, INDEX, true, "exclude", "")?;
    let structure = structure(Currency::JPY, as_of, &[pay_at], true)?;
    let ctx = Context {
        as_of,
        config: &config,
        refdata: &data,
        fixings: &FixingStore::default(),
    };
    let prepare = |records: &[ExerciseRecord]| {
        prepare(
            PrepareInput::Ir(&product),
            &CollateralKey::Uncollateralized,
            records,
            &structure,
            &ctx,
        )
    };
    assert!(matches!(prepare(&[]), Err(Error::Unsupported(_))));
    let mut record = ExerciseRecord {
        label: "switch".into(),
        at: exercise_at,
        effective_at: as_of,
        option: 0,
    };
    assert!(matches!(
        prepare(std::slice::from_ref(&record)),
        Err(Error::MissingData(_))
    ));
    record.option = 1;
    let prepared = prepare(&[record])?;
    let pv = evaluate(&prepared, &Market::from_theta(structure, vec![0.0])?)?.get(Currency::JPY);
    assert_close!(pv, 100.0, REL, JPY_ABS);
    Ok(())
}
`````

`````
crates/rhodium-pricing/tests/m1_linear.rs
//! M1 の標準取引を独立したキャッシュフローと閉形式に照合する。

use std::{num::NonZeroU32, sync::Arc};

use rhodium_core::{
    Error, Result, Warning,
    calendar::{
        BaseCalendar, BusinessDayConvention, Extrapolation as CalendarExtrapolation, WeekendRule,
    },
    currency::{Asset, Currency},
    date::{Date, DateInterval, Weekday},
    day_count::DayCount,
    fixing::{FixingStore, FixingStoreParts, ObservationPoint},
    fx::{CcyPair, FxRate},
    id::{CalendarId, IndexId},
    money::Money,
    quantity::{Rate, RateFuturesPrice, Simple, Spread},
    refdata::{
        Convention, FraConvention, FxConvention, FxForwardConvention, FxSpotLag, IndexKind,
        PublicationRule, RateFutureRate, RefData, RefDataParts,
    },
    schedule::{
        ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination, ScheduleSpec,
    },
    time::{Instant, model_time},
};
use rhodium_market::{
    curve::{CurveRole, CurveSpec, Extrapolation},
    fx::{FxAnchor, FxCollateral, FxSpec},
    market::Market,
    quote::{QuoteKind, QuoteSet},
    structure::{BlockId, MarketStructure},
};
use rhodium_pricing::{
    calibration::calibrate,
    config::{AnalyticPricer, Method, PricingConfig},
    context::Context,
    definition::{
        CalibrationTarget, DirectFx, ElementSelector, FxDefinition, FxQuotes, InstrumentBinding,
        MarketDefinition, QuoteTransform,
    },
    prepare::{PrepareInput, Prepared, evaluate, prepare},
    resolve::CollateralKey,
};
use rhodium_products::{
    ir::{
        Contract, Product,
        analysis::{Observation, observations},
    },
    standard::{
        Coupon, Fra, FxForward, FxSwap, InstrumentTemplate, Leg, Notional, NotionalExchange,
        PayReceive, PaymentRounding, Position, RateFuture, StandardProduct, Swap,
    },
};
use rhodium_testkit::assert_close;
use scioxide_ad::Tape;
use scioxide_core::{ExecutionContext, InvalidInputReason, SciOxideError, argument};

#[path = "fixtures/futures.rs"]
mod futures;

use futures::{
    AS_OF, CLEARING, FINAL, MULTIPLIER, PILLAR, RATE, START, at, config, curve, futures_definition,
    futures_template, market, refdata, reference_parts,
};

fn future(position: Position, price: f64) -> Result<StandardProduct> {
    Ok(StandardProduct::RateFuture(RateFuture::new(
        CLEARING.parse()?,
        NonZeroU32::MIN,
        RateFuturesPrice::new(price, DayCount::Act365F)?,
        position,
        None,
    )))
}

fn futures_pv(
    product: &StandardProduct,
    context: &Context<'_>,
    market: &Market<f64>,
) -> Result<(f64, Vec<Warning>)> {
    let trade_date = START.parse()?;
    let analytic = prepare(
        PrepareInput::Standard {
            product,
            trade_date,
        },
        &CollateralKey::Uncollateralized,
        &[],
        market.structure(),
        context,
    )?;
    let contract = product.to_product(trade_date, context.refdata)?;
    let ir = prepare(
        PrepareInput::Ir(&contract),
        &CollateralKey::Uncollateralized,
        &[],
        market.structure(),
        context,
    )?;
    assert_eq!(
        analytic.method(),
        Method::Analytic(AnalyticPricer::RateFuture)
    );
    assert_eq!(ir.method(), Method::Linear);
    let pv = evaluate(&analytic, market)?.get(Currency::USD);
    assert_close!(
        pv,
        evaluate(&ir, market)?.get(Currency::USD),
        1e-12_f64,
        1e-8_f64
    );
    assert_eq!(analytic.warnings(), ir.warnings());
    Ok((pv, analytic.warnings().to_vec()))
}

fn oracle_price(as_of: Instant, sessions: &[Instant], vol: f64) -> Result<f64> {
    let ts: Vec<_> = std::iter::once(0.0)
        .chain(sessions.iter().map(|&date| model_time(as_of, date).years()))
        .collect();
    let t1 = model_time(as_of, at(FINAL, "16:00")?).years();
    let t2 = model_time(as_of, at("2026-01-14", "16:00")?).years();
    let terminal = ts[ts.len() - 1];
    // a=0 の H は t。共分散 σ² min(s,t) に各清算時点の係数を直接掛ける。
    let covariance = ts
        .windows(2)
        .map(|w| (w[1] - w[0]) * w[0].min(t1))
        .sum::<f64>()
        - terminal * terminal.min(t1);
    let tau = 7.0 / 365.0;
    let log_mean = RATE * (t2 - t1) + vol * vol * (t2 - t1) * (t2 * t1 + covariance);
    Ok(1.0 - log_mean.exp_m1() / tau)
}

#[test]
fn futures_margin_ir_vs_analytic() -> Result<()> {
    // 日次に共通の価格を 実装と別の式で手計算し、将来の二回目以後の値洗いが消えることを確かめる。
    let as_of = AS_OF.parse()?;
    let reference = refdata("0D".parse()?, "16:00:00")?;
    let config = config("")?;
    let fixings = FixingStore::default();
    let context = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &fixings,
    };
    let margins = [
        at(START, "16:00")?,
        at("2026-01-06", "16:00")?,
        at(FINAL, "16:00")?,
    ];
    let df = (-RATE * model_time(as_of, margins[0]).years()).exp();
    for vol in [0.0, 0.5] {
        let market = market(as_of, vol)?;
        for (side, sign) in [(Position::Long, 1.0), (Position::Short, -1.0)] {
            let (actual, warnings) = futures_pv(&future(side, 0.95)?, &context, &market)?;
            let expected = sign * MULTIPLIER * df * (oracle_price(as_of, &margins, vol)? - 0.95);
            assert_close!(actual, expected, 1e-12_f64, 1e-8_f64);
            assert_eq!(warnings, []);
        }
    }
    Ok(())
}

#[test]
fn futures_settlement_delay() -> Result<()> {
    // 同日内の遅れも翌日の受渡も許可制の近似になり、禁止した設定では評価を拒否する。
    let as_of = AS_OF.parse()?;
    let market = market(as_of, 0.0)?;
    let config_ok = config("")?;
    let config_disabled = config("[conventions]\ndisabled = ['payment_delay_ignored']")?;
    for (lag, time) in [("0D", "18:00:00"), ("1BD", "08:00:00")] {
        let reference = refdata(lag.parse()?, time)?;
        let fixings = FixingStore::default();
        let context = Context {
            as_of,
            refdata: &reference,
            config: &config_ok,
            fixings: &fixings,
        };
        let product = future(Position::Long, 0.95)?;
        assert_eq!(
            futures_pv(&product, &context, &market)?.1,
            vec![Warning::PaymentDelayIgnored]
        );
        let disabled = Context {
            config: &config_disabled,
            ..context
        };
        assert!(matches!(
            prepare(
                PrepareInput::Standard {
                    product: &product,
                    trade_date: START.parse()?
                },
                &CollateralKey::Uncollateralized,
                &[],
                market.structure(),
                &disabled,
            ),
            Err(Error::Unsupported(_))
        ));
    }
    Ok(())
}

#[test]
fn futures_vs_forward_settlement() -> Result<()> {
    // 実装の先物価格が日次の値洗い口座 Q_M の式に一致し、最終日に一度だけ払う契約の
    // 価格（最終受渡の割引債を数値基準とする式）とは凸性の分だけ異なることを確かめる。
    // oracle 同士の比較だけでは、実装が満期一括の数値基準へ退行しても検出できないため、
    // ライブラリの評価から先物価格を逆算して照合する。
    let as_of = AS_OF.parse()?;
    let reference = refdata("0D".parse()?, "16:00:00")?;
    let config = config("")?;
    let fixings = FixingStore::default();
    let context = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &fixings,
    };
    let vol = 0.5;
    let market = market(as_of, vol)?;
    let first = at(START, "16:00")?;
    let final_at = at(FINAL, "16:00")?;
    let daily = [first, at("2026-01-06", "16:00")?, final_at];
    // PV は次の値洗い M·P(s_1)·(F − K) だけなので、そこから実装の F を求める。
    let trade_price = 0.95;
    let df = (-RATE * model_time(as_of, first).years()).exp();
    let (pv, _) = futures_pv(&future(Position::Long, trade_price)?, &context, &market)?;
    let implied = pv / (MULTIPLIER * df) + trade_price;
    assert_close!(
        implied,
        oracle_price(as_of, &daily, vol)?,
        1e-12_f64,
        1e-8_f64
    );
    let forward = oracle_price(as_of, &[final_at], vol)?;
    assert!((implied - forward).abs() > 1e-6);
    Ok(())
}

#[test]
fn futures_publication_and_payment_boundaries_use_realized_margin() -> Result<()> {
    // 確定済みで未払いの受渡は取引価格との差額。決済時点を除くと次回の差額だけ残る。
    let as_of = at(START, "16:00")?;
    let reference = refdata("0D".parse()?, "16:00:00")?;
    let config = config("")?;
    let product = future(Position::Long, 0.95)?;
    let market = market(as_of, 0.0)?;
    let missing = FixingStore::default();
    let ctx = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &missing,
    };
    assert!(matches!(
        prepare(
            PrepareInput::Standard {
                product: &product,
                trade_date: START.parse()?
            },
            &CollateralKey::Uncollateralized,
            &[],
            market.structure(),
            &ctx,
        ),
        Err(Error::MissingData(_))
    ));
    let fixings = FixingStore::new(
        &reference,
        FixingStoreParts {
            published: vec![(
                CLEARING.parse()?,
                ObservationPoint::Date(START.parse()?),
                0.954,
            )],
            ..FixingStoreParts::default()
        },
    )?;
    let context = Context {
        fixings: &fixings,
        ..ctx
    };
    let (excluded, _) = futures_pv(&product, &context, &market)?;

    let include: PricingConfig = toml::from_str(
        r#"
version = 1
reporting_currency = "USD"
settlement = "include"
[[discounting]]
currency = "USD"
collateral = "uncollateralized"
curve = "USD"
[forecasting]
TERM = "USD"
"#,
    )
    .map_err(|error| Error::InvalidData(error.to_string()))?;
    let included = futures_pv(
        &product,
        &Context {
            config: &include,
            ..context
        },
        &market,
    )?
    .0;
    assert_close!(included - excluded, 4.0_f64, 1e-12_f64, 1e-8_f64);
    Ok(())
}

#[test]
fn fra_settlement_matches_start_discounted_cashflow() -> Result<()> {
    // FRA は期末のクーポンではなく、同じターム金利で割り戻した差額を期首に受け渡す。
    let as_of = AS_OF.parse()?;
    let reference = refdata("0D".parse()?, "16:00:00")?;
    let config = config("")?;
    let fixings = FixingStore::default();
    let context = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &fixings,
    };
    let market = market(as_of, 0.0)?;
    let product = StandardProduct::Fra(Fra::new(
        Money::new(Currency::USD, 100_000.0),
        "TERM".parse()?,
        FINAL.parse()?,
        0.03,
        Position::Long,
    )?);
    let plan = prepare(
        PrepareInput::Standard {
            product: &product,
            trade_date: START.parse()?,
        },
        &CollateralKey::Uncollateralized,
        &[],
        market.structure(),
        &context,
    )?;
    let tau = 7.0 / 365.0;
    let forward = (RATE * tau).exp_m1() / tau;
    let df = (-RATE * model_time(as_of, at(FINAL, "16:00")?).years()).exp();
    let expected = df * 100_000.0 * tau * (forward - 0.03) / (1.0 + tau * forward);
    assert_close!(
        evaluate(&plan, &market)?.get(Currency::USD),
        expected,
        1e-12_f64,
        1e-8_f64
    );
    Ok(())
}

#[test]
fn futures_quote_calibrates_the_trade_to_par() -> Result<()> {
    // 先物のマークと独立した参照価格で校正し、同じ限月を約定した取引の PV がゼロであることを確かめる。
    let as_of = AS_OF.parse()?;
    let reference = refdata("0D".parse()?, "16:00:00")?;
    let config = config("")?;
    let context = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &FixingStore::default(),
    };
    let template = futures_template()?;
    let definition = futures_definition()?;
    let quote = oracle_price(
        as_of,
        &[
            at(START, "16:00")?,
            at("2026-01-06", "16:00")?,
            at(FINAL, "16:00")?,
        ],
        0.5,
    )?;
    let market = calibrate(
        &definition,
        &QuoteSet::new([
            ("USD.FUT.2026-01".parse()?, quote),
            ("USD.LGM.MR".parse()?, 0.0),
            ("USD.LGM.VOL".parse()?, 0.5),
        ])?,
        &context,
    )?;
    let product = template.build_for_valuation(quote, START.parse()?, as_of, &reference)?;
    assert_close!(
        futures_pv(&product, &context, market.market())?.0,
        0.0_f64,
        1e-12_f64,
        1e-8_f64
    );
    Ok(())
}

fn fx_reference(pair: CcyPair) -> Result<RefData> {
    RefData::new(fx_parts(pair)?)
}

fn fx_parts(pair: CcyPair) -> Result<RefDataParts> {
    let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
    let clock = parts.currencies[&Currency::USD].settlement_time.clone();
    parts.conventions.insert(
        "USDJPY-FWD".parse()?,
        Convention::FxForward(Box::new(FxForwardConvention {
            pair: "USDJPY".parse()?,
            bdc: BusinessDayConvention::Following,
        })),
    );
    // スポットラグ、受渡の暦（両通貨の休日の和）と受渡の時刻はペアの慣行から読む。
    let usd_jpy: CcyPair = "USDJPY".parse()?;
    parts.fx_conventions.insert(
        usd_jpy,
        FxConvention {
            pair: usd_jpy,
            spot_lag: FxSpotLag::TwoBusinessDays,
            settlement: clock,
            pip_exponent: 2,
        },
    );
    let mut currency = parts.currencies[&Currency::USD].clone();
    currency.code = Currency::JPY;
    currency.minor_units = 0;
    parts.currencies.insert(currency.code, currency);
    let mut term = parts.indices[&"TERM".parse()?].clone();
    term.id = "JPY-TERM".parse()?;
    if let IndexKind::TermRate { ccy, .. } = &mut term.kind {
        *ccy = Currency::JPY;
    }
    parts.indices.insert(term.id.clone(), term.clone());
    term.id = "FX-FIX".parse()?;
    term.kind = IndexKind::FxFixing { pair };
    parts.indices.insert(term.id.clone(), term);
    Ok(parts)
}

/// USDJPY の慣行のスポットラグを `spot_lag` にし、JPY の決済カレンダーに 2026-01-13（火）の
/// 祝日を置いた参照データ。ペアの暦は USD と JPY の休日の和になる。
fn fx_reference_with_lag(spot_lag: FxSpotLag) -> Result<RefData> {
    let mut parts = fx_parts("USDJPY".parse()?)?;
    let tokyo: CalendarId = "JPN".parse()?;
    parts.calendars.insert(
        tokyo.clone(),
        BaseCalendar::new(
            "2025-01-01".parse()?..="2026-12-31".parse()?,
            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
            vec!["2026-01-13".parse()?],
            CalendarExtrapolation::Strict,
        )?,
    );
    if let Some(jpy) = parts.currencies.get_mut(&Currency::JPY) {
        jpy.settlement_calendar = tokyo;
    }
    if let Some(convention) = parts.fx_conventions.get_mut(&"USDJPY".parse()?) {
        convention.spot_lag = spot_lag;
    }
    RefData::new(parts)
}

#[test]
fn mtm_xccy_reset_fixing_dates_follow_the_pair_spot_lag() -> Result<()> {
    // 元本のリセットは期首がスポット受渡になる日に観測する。期首は 2026-01-07（水）と 01-14（水）で、
    // ペアの暦は USD と JPY の休日の和（JPY は 01-13 が祝日）。
    // T+1: 01-07 の 1 営業日前は 01-06、01-14 の 1 営業日前は 01-13 を飛ばして 01-12。
    // T+2: 01-07 の 2 営業日前は 01-05、01-14 の 2 営業日前は 01-12、01-09。
    let swap = mtm_swap()?;
    for (spot_lag, expected) in [
        (FxSpotLag::OneBusinessDay, ["2026-01-06", "2026-01-12"]),
        (FxSpotLag::TwoBusinessDays, ["2026-01-05", "2026-01-09"]),
    ] {
        let reference = fx_reference_with_lag(spot_lag)?;
        let product = swap.to_product(START.parse()?, &reference)?;
        let fx: IndexId = "FX-FIX".parse()?;
        let mut dates: Vec<_> = observations(&product)
            .into_iter()
            .filter_map(|observation| match observation {
                Observation::Fixing {
                    index,
                    point: ObservationPoint::Date(date),
                } if *index == fx => Some(date),
                _ => None,
            })
            .collect();
        dates.sort_unstable();
        let expected = expected
            .iter()
            .map(|date| date.parse())
            .collect::<Result<Vec<Date>>>()?;
        assert_eq!(dates, expected, "{spot_lag:?}");
    }
    Ok(())
}

#[test]
fn fx_fixing_spot_value_uses_the_pair_spot_lag() -> Result<()> {
    // 受渡日の指定の無い FX fixing は、fixing 日からペアのスポットラグだけ進めた受渡日のレートを観測する。
    // 01-12（月）の fixing の受渡日は、T+1 では 01-13 の祝日を飛ばして 01-14、T+2 では 01-15。
    // 支払いは 01-20 で、その円の受渡を P_JPY で割り引く。F(v) = 150·exp(−(r_USD − r_JPY)·t_v)。
    let as_of = AS_OF.parse()?;
    let market = fx_market(as_of)?;
    let config = fx_config("")?;
    let fixings = FixingStore::default();
    let pay = at("2026-01-20", "16:00")?;
    let fx: IndexId = "FX-FIX".parse()?;
    let discount = |rate: f64, at| (-rate * model_time(as_of, at).years()).exp();
    let forward = |at| 150.0 * discount(RATE - 0.01, at);
    let mut values = Vec::new();
    for (spot_lag, value_date) in [
        (FxSpotLag::OneBusinessDay, "2026-01-14"),
        (FxSpotLag::TwoBusinessDays, "2026-01-15"),
    ] {
        let reference = fx_reference_with_lag(spot_lag)?;
        let product = Product::new(
            Contract::build(|builder| {
                let rate = builder.fixing(&fx, ObservationPoint::Date("2026-01-12".parse()?));
                Ok(builder.pay(Asset::Cash(Currency::JPY), rate, pay))
            })?,
            &reference,
        )?;
        let ctx = Context {
            as_of,
            refdata: &reference,
            config: &config,
            fixings: &fixings,
        };
        let plan = prepare(
            PrepareInput::Ir(&product),
            &CollateralKey::Uncollateralized,
            &[],
            market.structure(),
            &ctx,
        )?;
        // 受渡日と支払日が違うので、タイミング調整を省いた近似を記録する。
        assert!(plan.warnings().contains(&Warning::PaymentDelayIgnored));
        let pv = evaluate(&plan, &market)?.get(Currency::JPY);
        assert_close!(
            pv,
            forward(at(value_date, "16:00")?) * discount(0.01, pay),
            1e-12_f64,
            1e-6_f64
        );
        values.push(pv);
    }
    assert!((values[0] - values[1]).abs() > 1e-6);
    Ok(())
}

fn fx_market(as_of: Instant) -> Result<Market<f64>> {
    let jpy = CurveSpec::new(
        "JPY".parse()?,
        Currency::JPY,
        [CurveRole::Discount, CurveRole::Projection],
        [PILLAR.parse()?],
        Extrapolation::Reject,
    )?;
    let pair: CcyPair = "USDJPY".parse()?;
    let structure = Arc::new(MarketStructure::with_fx(
        as_of,
        [curve()?, jpy],
        [FxSpec {
            pair,
            anchor: FxAnchor::Direct {
                spot: as_of,
                short_starts: vec![],
                market_collateral: "csa".parse()?,
            },
            collaterals: vec![FxCollateral {
                id: "csa".parse()?,
                base: "USD".parse()?,
                quote: "JPY".parse()?,
            }],
            via: None,
        }],
    )?);
    let t = model_time(as_of, PILLAR.parse()?).years();
    let mut theta = vec![0.0; 3];
    for (id, value) in [
        (BlockId::nodes("USD".parse()?), -RATE * t),
        (BlockId::nodes("JPY".parse()?), -0.01 * t),
        (FxSpec::marks_block(pair)?, 150.0),
    ] {
        let range = structure
            .layout()
            .range(&id)
            .ok_or_else(|| Error::Internal("missing FX test block".into()))?;
        theta[range.start] = value;
    }
    Market::from_theta(structure, theta)
}

fn fx_config(conventions: &str) -> Result<PricingConfig> {
    toml::from_str(&format!(
        r#"
version = 1
reporting_currency = "USD"
[[discounting]]
currency = "USD"
collateral = "uncollateralized"
curve = "USD"
[[discounting]]
currency = "JPY"
collateral = "uncollateralized"
curve = "JPY"
[forecasting]
TERM = "USD"
JPY-TERM = "JPY"
{conventions}
"#
    ))
    .map_err(|error| Error::InvalidData(error.to_string()))
}

fn price(
    product: &StandardProduct,
    market: &Market<f64>,
    context: &Context<'_>,
) -> Result<rhodium_core::money::MoneyBag> {
    let plan = prepare(
        PrepareInput::Standard {
            product,
            trade_date: START.parse()?,
        },
        &CollateralKey::Uncollateralized,
        &[],
        market.structure(),
        context,
    )?;
    evaluate(&plan, market)
}

#[test]
fn fx_forward_and_swap_exchange_both_currencies() -> Result<()> {
    // 両通貨の各受渡を割り引き、近い受渡と遠い受渡の向きが逆になることを確かめる。
    let as_of = AS_OF.parse()?;
    let market = fx_market(as_of)?;
    let pair: CcyPair = "USDJPY".parse()?;
    let reference = fx_reference(pair)?;
    let config = fx_config("")?;
    let context = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &FixingStore::default(),
    };
    let near = at("2026-01-07", "16:00")?;
    let far = at("2026-01-21", "16:00")?;
    let notional = Money::new(Currency::USD, 1000.0);
    let first = FxRate::forward(pair, 151.0, near)?;
    let second = FxRate::forward(pair, 152.0, far)?;
    let forward = StandardProduct::FxForward(FxForward::new(notional, first, Position::Long)?);
    let swap = StandardProduct::FxSwap(FxSwap::new(notional, first, second, Position::Long)?);
    let df = |rate: f64, t| (-rate * model_time(as_of, t).years()).exp();

    let pv = price(&forward, &market, &context)?;
    assert_close!(
        pv.get(Currency::USD),
        1000.0 * df(RATE, near),
        1e-12_f64,
        1e-8_f64
    );
    assert_close!(
        pv.get(Currency::JPY),
        -151_000.0 * df(0.01, near),
        1e-12_f64,
        1e-6_f64
    );
    let pv = price(&swap, &market, &context)?;
    assert_close!(
        pv.get(Currency::USD),
        1000.0 * (df(RATE, near) - df(RATE, far)),
        1e-12_f64,
        1e-8_f64
    );
    assert_close!(
        pv.get(Currency::JPY),
        -151_000.0 * df(0.01, near) + 152_000.0 * df(0.01, far),
        1e-12_f64,
        1e-6_f64
    );

    // 市場の二つのアウトライトで約定すると、各受渡の報告通貨建て PV がゼロになる。
    let outright = |at| FxRate::forward(pair, 150.0 * df(RATE, at) / df(0.01, at), at);
    let par = StandardProduct::FxSwap(FxSwap::new(
        notional,
        outright(near)?,
        outright(far)?,
        Position::Long,
    )?);
    let pv = price(&par, &market, &context)?;
    assert_close!(
        pv.get(Currency::JPY) + 150.0 * pv.get(Currency::USD),
        0.0_f64,
        1e-12_f64,
        1e-6_f64
    );
    Ok(())
}

fn explicit_schedule() -> Result<ScheduleSpec> {
    let start = "2026-01-07".parse::<Date>()?;
    let middle = "2026-01-14".parse::<Date>()?;
    let end = "2026-01-21".parse::<Date>()?;
    let periods = [(start, middle), (middle, end)]
        .into_iter()
        .map(|(start, end)| {
            Ok(ExplicitPeriod {
                accrual: DateInterval::new(start, end)?,
                payment: end,
                reset: None,
                settlement: None,
                references: vec![],
            })
        })
        .collect::<Result<_>>()?;
    Ok(ScheduleSpec::Explicit(ExplicitPeriods::new(
        periods,
        ExplicitDatePolicy::KeepAsConfirmed,
        PeriodCombination::Sum,
        None,
    )?))
}

fn term(index: &str, spread: f64) -> Result<Coupon> {
    Ok(Coupon::Term {
        index: index.parse()?,
        spread: Spread::new(spread)?,
        gearing: 1.0,
    })
}

#[test]
fn irs_and_basis_swap_match_coupon_cashflows() -> Result<()> {
    // ターム対固定とターム対タームを同じ脚の構築 API で契約にし、スプレッドの符号を照合する。
    let as_of = AS_OF.parse()?;
    let market = fx_market(as_of)?;
    let reference = fx_reference("USDJPY".parse()?)?;
    let config = fx_config("")?;
    let context = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &FixingStore::default(),
    };
    let notional = Money::new(Currency::USD, 100_000.0);
    let schedule = explicit_schedule()?;
    let leg = |side, coupon| {
        Leg::new(
            side,
            notional,
            schedule.clone(),
            DayCount::Act365F,
            coupon,
            NotionalExchange::None,
            PaymentRounding::None,
        )
    };
    let tau = 7.0 / 365.0;
    let rate = (RATE * tau).exp_m1() / tau;
    let dfs = (-RATE * model_time(as_of, at("2026-01-14", "16:00")?).years()).exp()
        + (-RATE * model_time(as_of, at("2026-01-21", "16:00")?).years()).exp();
    for (pay, expected_coupon) in [
        (
            Coupon::Fixed {
                rate: Rate::new(0.03, DayCount::Act365F, Simple)?,
            },
            rate + 0.001 - 0.03,
        ),
        (term("TERM", 0.002)?, -0.001),
    ] {
        let swap = StandardProduct::Swap(Swap::from_legs([
            leg(PayReceive::Pay, pay)?,
            leg(PayReceive::Receive, term("TERM", 0.001)?)?,
        ])?);
        assert_close!(
            price(&swap, &market, &context)?.get(Currency::USD),
            notional.amount() * tau * dfs * expected_coupon,
            1e-12_f64,
            1e-8_f64
        );
    }
    Ok(())
}

fn mtm_swap() -> Result<StandardProduct> {
    let notional = Money::new(Currency::JPY, 150_000.0);
    let schedule = explicit_schedule()?;
    Ok(StandardProduct::Swap(Swap::from_legs([
        Leg::new(
            PayReceive::Receive,
            notional,
            schedule.clone(),
            DayCount::Act365F,
            term("JPY-TERM", 0.001)?,
            NotionalExchange::InitialAndFinal,
            PaymentRounding::None,
        )?,
        Leg::with_notional(
            PayReceive::Pay,
            Notional::FxReset {
                base: notional,
                fixing: "FX-FIX".parse()?,
            },
            schedule,
            DayCount::Act365F,
            term("TERM", 0.002)?,
            NotionalExchange::InitialAndFinal,
            PaymentRounding::None,
        )?,
    ])?))
}

#[test]
fn mtm_xccy_fx_reset_matches_cashflows_in_both_directions() -> Result<()> {
    // 同じ基準元本を正向きの fixing では割り、逆向きでは掛ける。期中の交換は元本の差額だけ。
    let as_of = AS_OF.parse()?;
    let market = fx_market(as_of)?;
    let config = fx_config("")?;
    let disabled = fx_config("[conventions]\ndisabled = ['mtm_fx_forward_product']")?;
    let swap = mtm_swap()?;
    let fixings = FixingStore::default();
    let pay = ["2026-01-07", "2026-01-14", "2026-01-21"].map(|date| at(date, "16:00"));
    let [start, middle, end] = pay;
    let [start, middle, end] = [start?, middle?, end?];
    let discount = |rate: f64, at| (-rate * model_time(as_of, at).years()).exp();
    let loan = |at| 1000.0 * discount(0.01, at) / discount(RATE, at);
    let [n0, n1] = [loan(start), loan(middle)];
    let dollar_coupon = (RATE * (7.0 / 365.0)).exp_m1() + 0.002 * (7.0 / 365.0);
    let yen_coupon = (0.01_f64 * (7.0 / 365.0)).exp_m1() + 0.001 * (7.0 / 365.0);
    let usd = n0 * discount(RATE, start) + (n1 - n0 - n0 * dollar_coupon) * discount(RATE, middle)
        - n1 * (1.0 + dollar_coupon) * discount(RATE, end);
    let jpy = 150_000.0
        * (-discount(0.01, start)
            + yen_coupon * discount(0.01, middle)
            + (1.0 + yen_coupon) * discount(0.01, end));

    for pair in ["USDJPY", "JPYUSD"] {
        let reference = fx_reference(pair.parse()?)?;
        let ctx = Context {
            as_of,
            refdata: &reference,
            config: &config,
            fixings: &fixings,
        };
        let plan = prepare(
            PrepareInput::Standard {
                product: &swap,
                trade_date: START.parse()?,
            },
            &CollateralKey::Uncollateralized,
            &[],
            market.structure(),
            &ctx,
        )?;
        assert!(
            plan.blocks()
                .contains(&FxSpec::marks_block("USDJPY".parse()?)?)
        );
        assert_eq!(plan.warnings(), &[Warning::MtmFxForwardProduct]);
        let actual = evaluate(&plan, &market)?;
        assert_close!(actual.get(Currency::USD), usd, 1e-12_f64, 1e-8_f64);
        assert_close!(actual.get(Currency::JPY), jpy, 1e-12_f64, 1e-6_f64);
        assert!(matches!(
            price(
                &swap,
                &market,
                &Context {
                    config: &disabled,
                    ..ctx
                }
            ),
            Err(Error::Unsupported(_))
        ));
    }
    Ok(())
}

/// 固定クーポンの `MtM` の XCCY。JPY の脚は元本 150,000 円を固定し、USD の脚の元本を `FX-FIX` で
/// リセットする。`reset` は 2 期目のリセットの明示日（無ければ規則どおり期首から遡る）。
fn fixed_mtm_swap(dates: &[&str; 3], reset: Option<&str>) -> Result<StandardProduct> {
    let schedule = |reset: Option<&str>| -> Result<ScheduleSpec> {
        let periods = dates
            .windows(2)
            .enumerate()
            .map(|(i, pair)| {
                let (start, end) = (pair[0].parse::<Date>()?, pair[1].parse::<Date>()?);
                Ok(ExplicitPeriod {
                    accrual: DateInterval::new(start, end)?,
                    payment: end,
                    reset: if i == 1 {
                        reset.map(str::parse).transpose()?
                    } else {
                        None
                    },
                    settlement: None,
                    references: vec![],
                })
            })
            .collect::<Result<_>>()?;
        Ok(ScheduleSpec::Explicit(ExplicitPeriods::new(
            periods,
            ExplicitDatePolicy::KeepAsConfirmed,
            PeriodCombination::Sum,
            None,
        )?))
    };
    let fixed = || -> Result<Coupon> {
        Ok(Coupon::Fixed {
            rate: Rate::new(0.03, DayCount::Act365F, Simple)?,
        })
    };
    let notional = Money::new(Currency::JPY, 150_000.0);
    Ok(StandardProduct::Swap(Swap::from_legs([
        Leg::new(
            PayReceive::Receive,
            notional,
            schedule(None)?,
            DayCount::Act365F,
            fixed()?,
            NotionalExchange::InitialAndFinal,
            PaymentRounding::None,
        )?,
        Leg::with_notional(
            PayReceive::Pay,
            Notional::FxReset {
                base: notional,
                fixing: "FX-FIX".parse()?,
            },
            schedule(reset)?,
            DayCount::Act365F,
            fixed()?,
            NotionalExchange::InitialAndFinal,
            PaymentRounding::None,
        )?,
    ])?))
}

#[test]
fn mtm_xccy_reset_reads_the_fixing_at_the_pair_spot_delivery() -> Result<()> {
    // 同じ指標・観測日の FX fixing の公表値は一つなので、元本リセットも通常の fixing と同じく
    // ペアのスポットの受渡（ここでは UTC 06 時。東京 15 時に当たる）の F で予想する。
    // 新しい元本を交換する瞬間（USD の決済時刻 UTC 16 時）とは別に扱う。
    // 明示の reset 日 01-08 の fixing はスポットの受渡が 01-12 で、元本の交換 01-14 より早いので
    // 支払いの遅れの近似を記録する。クーポンは固定にして、FX の予想だけを照合する。
    let as_of = AS_OF.parse()?;
    let market = fx_market(as_of)?;
    let config = fx_config("")?;
    let fixings = FixingStore::default();
    let mut parts = fx_parts("USDJPY".parse()?)?;
    if let Some(convention) = parts.fx_conventions.get_mut(&"USDJPY".parse()?) {
        convention.settlement.time = "06:00:00".parse()?;
    }
    let reference = RefData::new(parts)?;
    let ctx = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &fixings,
    };
    let discount = |rate: f64, at| (-rate * model_time(as_of, at).years()).exp();
    // F(v) = 150·exp(−(r_USD − r_JPY)·t_v)。リセット後の USD の元本は 150,000 円 / F(v)。
    let forward = |at| 150.0 * discount(RATE - 0.01, at);
    let fx: IndexId = "FX-FIX".parse()?;

    // 通常の fixing: 01-05 の fixing の受渡は T+2 の 01-07 06 時。
    let value = at("2026-01-07", "06:00")?;
    let pay = at("2026-01-07", "16:00")?;
    let plain = Product::new(
        Contract::build(|builder| {
            let rate = builder.fixing(&fx, ObservationPoint::Date("2026-01-05".parse()?));
            Ok(builder.pay(Asset::Cash(Currency::JPY), rate, pay))
        })?,
        &reference,
    )?;
    let plan = prepare(
        PrepareInput::Ir(&plain),
        &CollateralKey::Uncollateralized,
        &[],
        market.structure(),
        &ctx,
    )?;
    assert_close!(
        evaluate(&plan, &market)?.get(Currency::JPY),
        forward(value) * discount(0.01, pay),
        1e-12_f64,
        1e-6_f64
    );

    let dates = ["2026-01-07", "2026-01-14", "2026-01-21"];
    let [start, middle, end] = dates.map(|date| at(date, "16:00"));
    let [start, middle, end] = [start?, middle?, end?];
    let coupon = 0.03 * 7.0 / 365.0;
    // 支払う USD の脚: 期首に元本 n0 を受け、期中に元本の差額 n1 − n0 を受けてクーポンを払い、
    // 期末に n1 とクーポンを払う。
    let usd = |n0: f64, n1: f64| {
        n0 * discount(RATE, start) + (n1 - n0 - n0 * coupon) * discount(RATE, middle)
            - n1 * (1.0 + coupon) * discount(RATE, end)
    };
    for (reset, second_value, warnings) in [
        // 規則どおりの reset はスポットの受渡が期首の 01-14 と同じ日なので近似は無い。
        (None, "2026-01-14", vec![]),
        // 明示の reset 01-08（木）の T+2 の受渡は 01-12（月）。
        (
            Some("2026-01-08"),
            "2026-01-12",
            vec![Warning::PaymentDelayIgnored],
        ),
    ] {
        let product = fixed_mtm_swap(&dates, reset)?;
        let plan = prepare(
            PrepareInput::Standard {
                product: &product,
                trade_date: START.parse()?,
            },
            &CollateralKey::Uncollateralized,
            &[],
            market.structure(),
            &ctx,
        )?;
        assert_eq!(plan.warnings(), warnings.as_slice(), "{reset:?}");
        let n0 = 150_000.0 / forward(value);
        let n1 = 150_000.0 / forward(at(second_value, "06:00")?);
        assert_close!(
            evaluate(&plan, &market)?.get(Currency::USD),
            usd(n0, n1),
            1e-12_f64,
            1e-8_f64
        );
    }
    Ok(())
}

#[test]
fn overnight_futures_compound_or_average_the_reference_period() -> Result<()> {
    // 先物の清算指標を切り替えず、最終値の参照金利だけを複利・加重平均で組み立てる。
    let as_of = AS_OF.parse()?;
    let market = market(as_of, 0.0)?;
    let config = config("")?;
    let accrual = DateInterval::new("2026-01-05".parse()?, "2026-01-07".parse()?)?;
    let tau = 1.0 / 365.0;
    for (reference_rate, expected) in [
        (
            RateFutureRate::Compounded(accrual),
            (RATE * 2.0 * tau).exp_m1() / (2.0 * tau),
        ),
        (
            RateFutureRate::Averaged(accrual),
            (RATE * tau).exp_m1() / tau,
        ),
    ] {
        let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
        if let Some(term) = parts.indices.get_mut(&"TERM".parse()?) {
            term.kind = IndexKind::OvernightRate {
                ccy: Currency::USD,
                day_count: DayCount::Act365F,
            };
        }
        if let Some(Convention::RateFuture(convention)) =
            parts.conventions.get_mut(&"USD-RATE-FUTURES".parse()?)
            && let Some(contract) = convention.contracts.values_mut().next()
        {
            contract.rate = reference_rate;
        }
        let reference = RefData::new(parts)?;
        let ctx = Context {
            as_of,
            refdata: &reference,
            config: &config,
            fixings: &FixingStore::default(),
        };
        let df = (-RATE * model_time(as_of, at(START, "16:00")?).years()).exp();
        let pv = futures_pv(&future(Position::Long, 0.95)?, &ctx, &market)?.0;
        assert_close!(
            pv,
            (1.0 - expected - 0.95) * df * MULTIPLIER,
            1e-12_f64,
            1e-8_f64
        );
    }
    Ok(())
}

#[test]
fn overnight_futures_annualize_the_spanning_day_over_its_remaining_fraction() -> Result<()> {
    // 評価時点 A が翌日物の観測の期間をまたぎ、その日の率が未公表のとき、[A, 終了) の forward の率を
    // 残りの日数係数 τ_rem で年率化し、その率を元の 1 日の日数係数に当てる（参照データの章の
    // fixing の規則）。成長を [A, 終了) に縮めたまま元の重みを掛けると、経過した分の利息が抜ける。
    // 複利と加重平均のそれぞれを手計算と照合し、`futures_pv` で IR のリニア評価と解析の一致も確かめる。
    let as_of = at(START, "20:00")?;
    let market = market(as_of, 0.0)?;
    let config = config("")?;
    let accrual = DateInterval::new("2026-01-05".parse()?, "2026-01-07".parse()?)?;
    let day = 1.0 / 365.0;
    let first_end = at("2026-01-06", "16:00")?;
    let second_end = at("2026-01-07", "16:00")?;
    // 1 日目は A（20 時）から翌日の観測時刻（16 時）までの 20 時間だけが残る。
    let remaining = 20.0 / 24.0 * day;
    let first_growth = (RATE * model_time(as_of, first_end).years()).exp();
    let second_growth = (RATE
        * (model_time(as_of, second_end).years() - model_time(as_of, first_end).years()))
    .exp();
    let first_rate = (first_growth - 1.0) / remaining;
    let second_rate = (second_growth - 1.0) / day;
    for (reference_rate, expected) in [
        (
            RateFutureRate::Compounded(accrual),
            ((1.0 + day * first_rate) * (1.0 + day * second_rate) - 1.0) / (2.0 * day),
        ),
        (
            RateFutureRate::Averaged(accrual),
            f64::midpoint(first_rate, second_rate),
        ),
    ] {
        let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
        if let Some(term) = parts.indices.get_mut(&"TERM".parse()?) {
            term.kind = IndexKind::OvernightRate {
                ccy: Currency::USD,
                day_count: DayCount::Act365F,
            };
            // 翌日物は翌営業日の 16 時に公表する。A では 1 日目が未公表で、期間が A をまたぐ。
            if let PublicationRule::Daily(rule) = &mut term.publication {
                rule.business_days = 1;
            }
        }
        if let Some(Convention::RateFuture(convention)) =
            parts.conventions.get_mut(&"USD-RATE-FUTURES".parse()?)
            && let Some(contract) = convention.contracts.values_mut().next()
        {
            contract.rate = reference_rate;
        }
        let reference = RefData::new(parts)?;
        // 初日の清算価格は 16 時に公表済みなので、最新の決済価格 P_last として渡す。
        let fixings = FixingStore::new(
            &reference,
            FixingStoreParts {
                published: vec![(
                    CLEARING.parse()?,
                    ObservationPoint::Date(START.parse()?),
                    0.954,
                )],
                ..FixingStoreParts::default()
            },
        )?;
        let ctx = Context {
            as_of,
            refdata: &reference,
            config: &config,
            fixings: &fixings,
        };
        // 次の値洗いの受渡は翌日の 16 時。
        let df = (-RATE * model_time(as_of, first_end).years()).exp();
        let (pv, warnings) = futures_pv(&future(Position::Long, 0.95)?, &ctx, &market)?;
        assert_close!(
            pv,
            (1.0 - expected - 0.954) * df * MULTIPLIER,
            1e-12_f64,
            1e-8_f64
        );
        assert!(
            warnings
                .iter()
                .any(|warning| matches!(warning, Warning::UnpublishedFixingForecast { .. }))
        );
    }
    Ok(())
}

fn check_gradient(
    plan: &Prepared,
    market: &Market<f64>,
    block: &BlockId,
    element: usize,
    bump: f64,
) -> Result<f64> {
    let structure = market.structure();
    let mut gradient = vec![0.0; market.theta().len()];
    // 準備済みの商品をそのまま再利用し、AAD と中心差分の双方で同じマークを動かす。
    Tape::default()
        .gradient_into(
            |theta| {
                Market::from_theta(structure.clone(), theta.to_vec())
                    .and_then(|market| evaluate(plan, &market))
                    .map(|pv| pv.get(Currency::USD))
                    .map_err(|error| {
                        SciOxideError::invalid_input(
                            argument!("market"),
                            InvalidInputReason::Inconsistent,
                            error.to_string(),
                        )
                    })
            },
            market.theta(),
            &mut gradient,
            ExecutionContext::sequential(),
        )
        .map_err(|error| Error::Numerical(error.to_string()))?;
    let position = structure
        .layout()
        .range(block)
        .ok_or_else(|| Error::Internal("missing AAD block".into()))?
        .start
        + element;
    let evaluate_at = |delta: f64| -> Result<f64> {
        let mut theta = market.theta().to_vec();
        theta[position] += delta;
        let market = Market::from_theta(structure.clone(), theta)?;
        Ok(evaluate(plan, &market)?.get(Currency::USD))
    };
    let central = (evaluate_at(bump)? - evaluate_at(-bump)?) / (2.0 * bump);
    assert_close!(gradient[position], central, 1e-7_f64, 5e-8_f64);
    Ok(gradient[position])
}

#[test]
fn futures_convexity_and_mtm_fx_marks_have_aad_gradients() -> Result<()> {
    // LGM ボラと評価時点の FX マークをそれぞれ動かし、値洗いと元本リセットへ感応度が流れることを確かめる。
    let as_of = AS_OF.parse()?;
    let fixings = FixingStore::default();
    let futures_market = market(as_of, 0.5)?;
    let futures_reference = refdata("0D".parse()?, "16:00:00")?;
    let futures_config = config("")?;
    let futures_context = Context {
        as_of,
        refdata: &futures_reference,
        config: &futures_config,
        fixings: &fixings,
    };
    let product = future(Position::Long, 0.95)?;
    let plan = prepare(
        PrepareInput::Standard {
            product: &product,
            trade_date: START.parse()?,
        },
        &CollateralKey::Uncollateralized,
        &[],
        futures_market.structure(),
        &futures_context,
    )?;
    assert!(
        check_gradient(
            &plan,
            &futures_market,
            &MarketStructure::lgm_block(Currency::USD)?,
            1,
            1e-3
        )?
        .abs()
            > 1e-4
    );

    let fx = fx_market(as_of)?;
    let reference = fx_reference("USDJPY".parse()?)?;
    let fx_config = fx_config("")?;
    let fx_context = Context {
        as_of,
        refdata: &reference,
        config: &fx_config,
        fixings: &fixings,
    };
    let swap = mtm_swap()?;
    let plan = prepare(
        PrepareInput::Standard {
            product: &swap,
            trade_date: START.parse()?,
        },
        &CollateralKey::Uncollateralized,
        &[],
        fx.structure(),
        &fx_context,
    )?;
    // USD 脚にスプレッドを入れず同じ曲線で割り引くため、USD PV はほぼゼロ。支払いの連鎖を微分できることを照合する。
    let _gradient = check_gradient(
        &plan,
        &fx,
        &FxSpec::marks_block("USDJPY".parse()?)?,
        0,
        0.01,
    )?;
    Ok(())
}

#[test]
fn fra_quote_calibrates_the_trade_to_par() -> Result<()> {
    // スポットから指定した 0D×1W FRA の par quote で曲線を解き、同一条件の取引を評価する。
    let as_of = AS_OF.parse()?;
    let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
    parts.conventions.insert(
        "USD-FRA".parse()?,
        Convention::Fra(FraConvention {
            index: "TERM".parse()?,
            spot_lag: "2BD".parse()?,
        }),
    );
    let reference = RefData::new(parts)?;
    let config = config("")?;
    let context = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &FixingStore::default(),
    };
    let template = InstrumentTemplate::Fra {
        convention: "USD-FRA".parse()?,
        start: "0D".parse()?,
        end: "1W".parse()?,
    };
    let quote_id = "USD.FRA.1W".parse()?;
    let binding = InstrumentBinding::new(
        quote_id,
        QuoteKind::Rate,
        CalibrationTarget::ParRate(template.clone()),
        QuoteTransform::Identity,
        CollateralKey::Uncollateralized,
        ElementSelector::Block(BlockId::nodes("USD".parse()?)),
    )?;
    let definition =
        MarketDefinition::new("USD".parse()?, "Etc/UTC".parse()?, [curve()?], [binding])?;
    let tau = 7.0 / 365.0;
    let quote = (RATE * tau).exp_m1() / tau;
    let calibrated = calibrate(
        &definition,
        &QuoteSet::new([("USD.FRA.1W".parse()?, quote)])?,
        &context,
    )?;
    let trade = template.build_for_trade(quote, START.parse()?, &reference)?;
    let pv = price(&trade, calibrated.market(), &context)?;
    assert_close!(pv.get(Currency::USD), 0.0_f64, 1e-12_f64, 1e-8_f64);
    Ok(())
}

#[test]
fn fx_outright_quote_calibrates_a_par_trade_in_both_currencies() -> Result<()> {
    // outright の校正は base と quote の両方の PV を評価時点の FX で合算する。
    let as_of = AS_OF.parse()?;
    let pair: CcyPair = "USDJPY".parse()?;
    let reference = fx_reference(pair)?;
    let config = fx_config("")?;
    let context = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &FixingStore::default(),
    };
    let template = InstrumentTemplate::FxForward {
        convention: "USDJPY-FWD".parse()?,
        tenor: "1W".parse()?,
    };
    let discount = CurveSpec::new(
        "JPY".parse()?,
        Currency::JPY,
        [CurveRole::Discount],
        [PILLAR.parse()?],
        Extrapolation::Reject,
    )?;
    let usd = InstrumentBinding::new(
        "USD.DF.1".parse()?,
        QuoteKind::Mark,
        CalibrationTarget::DirectMark,
        QuoteTransform::Identity,
        CollateralKey::Uncollateralized,
        ElementSelector::Param(rhodium_market::structure::ParamId::new(
            BlockId::nodes("USD".parse()?),
            0,
        )),
    )?;
    let jpy = InstrumentBinding::new(
        "USDJPY.FWD.1W".parse()?,
        QuoteKind::Price,
        CalibrationTarget::Price(template.clone()),
        QuoteTransform::Identity,
        CollateralKey::Uncollateralized,
        ElementSelector::Block(BlockId::nodes("JPY".parse()?)),
    )?;
    let definition = MarketDefinition::new(
        "USD".parse()?,
        "Etc/UTC".parse()?,
        [curve()?, discount],
        [usd, jpy],
    )?
    .with_fx([FxDefinition {
        pair,
        direct: Some(DirectFx {
            market_collateral: "csa".parse()?,
            quotes: FxQuotes {
                spot: "USDJPY.SPOT".parse()?,
                on: None,
                tn: None,
                sn: None,
                forwards: vec![],
            },
            implied: None,
        }),
        collaterals: vec![FxCollateral {
            id: "csa".parse()?,
            base: "USD".parse()?,
            quote: "JPY".parse()?,
        }],
        via: None,
    }])?;
    let quote = 150.0 * (-(RATE - 0.01) * 7.0 / 365.0).exp();
    let calibrated = calibrate(
        &definition,
        &QuoteSet::new([
            (
                "USD.DF.1".parse()?,
                -RATE * model_time(as_of, PILLAR.parse()?).years(),
            ),
            ("USDJPY.FWD.1W".parse()?, quote),
            ("USDJPY.SPOT".parse()?, 150.0),
        ])?,
        &context,
    )?;
    let trade = template.build_for_trade(quote, START.parse()?, &reference)?;
    let market = calibrated.market();
    let pv = price(&trade, market, &context)?;
    let fx = market.with_view(|view| view.fx_today(pair))?;
    assert_close!(
        pv.convert(Currency::JPY, &[fx])?.amount(),
        0.0_f64,
        1e-12_f64,
        1e-6_f64
    );
    Ok(())
}

#[test]
fn fixed_notional_xccy_exchanges_both_currencies() -> Result<()> {
    // 元本を固定した XCCY は、各通貨の脚が期首に元本を払い、期末に元本を受け取り、クーポンを
    // それぞれの通貨の曲線で割り引く（FX のリセットは無い）。JPY の脚は 1% の固定、USD の脚は
    // ターム金利の 1W（RATE の連続金利から求めた単利）にスプレッド 20bp。
    let as_of = AS_OF.parse()?;
    let market = fx_market(as_of)?;
    let config = fx_config("")?;
    let reference = fx_reference("USDJPY".parse()?)?;
    let ctx = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &FixingStore::default(),
    };
    let schedule = explicit_schedule()?;
    let swap = StandardProduct::Swap(Swap::from_legs([
        Leg::new(
            PayReceive::Receive,
            Money::new(Currency::JPY, 150_000.0),
            schedule.clone(),
            DayCount::Act365F,
            Coupon::Fixed {
                rate: Rate::new(0.01, DayCount::Act365F, Simple)?,
            },
            NotionalExchange::InitialAndFinal,
            PaymentRounding::None,
        )?,
        Leg::new(
            PayReceive::Pay,
            Money::new(Currency::USD, 1_000.0),
            schedule,
            DayCount::Act365F,
            term("TERM", 0.002)?,
            NotionalExchange::InitialAndFinal,
            PaymentRounding::None,
        )?,
    ])?);
    let plan = prepare(
        PrepareInput::Standard {
            product: &swap,
            trade_date: START.parse()?,
        },
        &CollateralKey::Uncollateralized,
        &[],
        market.structure(),
        &ctx,
    )?;
    // 元本が固定なので、FX と金利の積の近似（MtmFxForwardProduct）は使わない。
    assert!(!plan.warnings().contains(&Warning::MtmFxForwardProduct));
    let pv = evaluate(&plan, &market)?;
    let discount = |rate: f64, date: &str| -> Result<f64> {
        Ok((-rate * model_time(as_of, at(date, "16:00")?).years()).exp())
    };
    let tau = 7.0 / 365.0;
    let [start, middle, end] = ["2026-01-07", "2026-01-14", "2026-01-21"];
    let jpy = 150_000.0
        * (-discount(0.01, start)?
            + 0.01 * tau * discount(0.01, middle)?
            + (1.0 + 0.01 * tau) * discount(0.01, end)?);
    let coupon = (RATE * tau).exp_m1() + 0.002 * tau;
    let usd = 1_000.0
        * (discount(RATE, start)?
            - coupon * discount(RATE, middle)?
            - (1.0 + coupon) * discount(RATE, end)?);
    assert_close!(pv.get(Currency::JPY), jpy, 1e-12_f64, 1e-6_f64);
    assert_close!(pv.get(Currency::USD), usd, 1e-12_f64, 1e-8_f64);
    Ok(())
}

#[test]
fn term_coupon_payment_delay_is_recorded() -> Result<()> {
    // ターム金利のクーポンを期末から 2 日遅れて払うと、自然な支払日（期末）とのずれのタイミング調整を
    // 省く近似（PaymentDelayIgnored）を記録し、予想は期末の DF 比、割引は実際の支払日で求める。
    // 近似を禁止すると準備できない。
    let as_of = AS_OF.parse()?;
    let market = market(as_of, 0.0)?;
    let reference = refdata("0D".parse()?, "16:00:00")?;
    let fixings = FixingStore::default();
    let schedule = ScheduleSpec::Explicit(ExplicitPeriods::new(
        vec![ExplicitPeriod {
            accrual: DateInterval::new("2026-01-07".parse()?, "2026-01-14".parse()?)?,
            payment: "2026-01-16".parse()?,
            reset: None,
            settlement: None,
            references: vec![],
        }],
        ExplicitDatePolicy::KeepAsConfirmed,
        PeriodCombination::Sum,
        None,
    )?);
    let swap = StandardProduct::Swap(Swap::from_legs([
        Leg::new(
            PayReceive::Receive,
            Money::new(Currency::USD, 100_000.0),
            schedule.clone(),
            DayCount::Act365F,
            term("TERM", 0.0)?,
            NotionalExchange::None,
            PaymentRounding::None,
        )?,
        Leg::new(
            PayReceive::Pay,
            Money::new(Currency::USD, 100_000.0),
            schedule,
            DayCount::Act365F,
            Coupon::Fixed {
                rate: Rate::new(0.0, DayCount::Act365F, Simple)?,
            },
            NotionalExchange::None,
            PaymentRounding::None,
        )?,
    ])?);
    let prepare_with = |config: &PricingConfig| {
        prepare(
            PrepareInput::Standard {
                product: &swap,
                trade_date: START.parse()?,
            },
            &CollateralKey::Uncollateralized,
            &[],
            market.structure(),
            &Context {
                as_of,
                refdata: &reference,
                config,
                fixings: &fixings,
            },
        )
    };
    let plan = prepare_with(&config("")?)?;
    assert!(plan.warnings().contains(&Warning::PaymentDelayIgnored));
    let tau = 7.0 / 365.0;
    let forward = (RATE * tau).exp_m1() / tau;
    let df = (-RATE * model_time(as_of, at("2026-01-16", "16:00")?).years()).exp();
    assert_close!(
        evaluate(&plan, &market)?.get(Currency::USD),
        100_000.0 * tau * forward * df,
        1e-12_f64,
        1e-8_f64
    );
    assert!(matches!(
        prepare_with(&config(
            "[conventions]\ndisabled = ['payment_delay_ignored']"
        )?),
        Err(Error::Unsupported(_))
    ));
    Ok(())
}

/// [`fx_market`] に、JPY を 25bp 上乗せした割引曲線 `JPY-25` を足した市場。
fn fx_market_with_spread_curve(as_of: Instant, pair: CcyPair) -> Result<Market<f64>> {
    let jpy = CurveSpec::new(
        "JPY".parse()?,
        Currency::JPY,
        [CurveRole::Discount, CurveRole::Projection],
        [PILLAR.parse()?],
        Extrapolation::Reject,
    )?;
    let spread = CurveSpec::from_expr(
        "JPY-25".parse()?,
        Currency::JPY,
        [CurveRole::Discount],
        rhodium_market::curve::CurveExpr::Sum(vec![
            rhodium_market::curve::CurveExpr::Ref("JPY".parse()?),
            rhodium_market::curve::CurveExpr::Flat(0.0025),
        ]),
    )?;
    let structure = Arc::new(MarketStructure::with_fx(
        as_of,
        [curve()?, jpy, spread],
        [FxSpec {
            pair,
            anchor: FxAnchor::Direct {
                spot: as_of,
                short_starts: vec![],
                market_collateral: "csa".parse()?,
            },
            collaterals: vec![FxCollateral {
                id: "csa".parse()?,
                base: "USD".parse()?,
                quote: "JPY".parse()?,
            }],
            via: None,
        }],
    )?);
    let t = model_time(as_of, PILLAR.parse()?).years();
    let mut theta = vec![0.0; structure.layout().len()];
    for (id, value) in [
        (BlockId::nodes("USD".parse()?), -RATE * t),
        (BlockId::nodes("JPY".parse()?), -0.01 * t),
        (FxSpec::marks_block(pair)?, 150.0),
    ] {
        let range = structure
            .layout()
            .range(&id)
            .ok_or_else(|| Error::Internal("missing FX test block".into()))?;
        theta[range.start] = value;
    }
    let market = Market::from_theta(structure, theta)?;
    Ok(market)
}

#[test]
fn fx_fixing_under_another_collateral_records_cross_collateral_projection() -> Result<()> {
    // FX の予想は、ペアの FX の担保区分（USD と JPY の曲線の組）の下で定まる。支払いを別の担保
    // （JPY を 25bp 上乗せした曲線で割り引く CSA）で受け取ると、予想と割引の担保が違うので
    // CrossCollateralProjection を記録し、禁止すると準備できない。
    let as_of = AS_OF.parse()?;
    let reference = fx_reference("USDJPY".parse()?)?;
    let fixings = FixingStore::default();
    let pair: CcyPair = "USDJPY".parse()?;
    let market = fx_market_with_spread_curve(as_of, pair)?;
    let product = Product::new(
        Contract::build(|builder| {
            let rate = builder.fixing(
                &"FX-FIX".parse()?,
                ObservationPoint::Date("2026-01-12".parse()?),
            );
            Ok(builder.pay(Asset::Cash(Currency::JPY), rate, at("2026-01-14", "16:00")?))
        })?,
        &reference,
    )?;
    let collateral = CollateralKey::Cash(rhodium_products::agreement::CashCollateralTerms {
        ccy: Currency::JPY,
        remuneration: rhodium_products::agreement::Remuneration {
            index: "JPY-TERM".parse()?,
            spread: Spread::new(0.0025)?,
            floor: None,
        },
    });
    let rules = r#"
[[discounting]]
currency = "JPY"
collateral = { cash = { ccy = "JPY", index = "JPY-TERM", spread = "25bp" } }
curve = "JPY-25"
[[discounting]]
currency = "USD"
collateral = { cash = { ccy = "JPY", index = "JPY-TERM", spread = "25bp" } }
curve = "USD"
"#;
    let prepare_under = |config: &PricingConfig, collateral: &CollateralKey| {
        prepare(
            PrepareInput::Ir(&product),
            collateral,
            &[],
            market.structure(),
            &Context {
                as_of,
                refdata: &reference,
                config,
                fixings: &fixings,
            },
        )
    };
    let native = prepare_under(&fx_config("")?, &CollateralKey::Uncollateralized)?;
    assert!(
        !native
            .warnings()
            .contains(&Warning::CrossCollateralProjection)
    );
    let spread_config = fx_config(rules)?;
    let foreign = prepare_under(&spread_config, &collateral)?;
    assert!(
        foreign
            .warnings()
            .contains(&Warning::CrossCollateralProjection)
    );
    // 予想は担保区分の曲線のフォワードのまま、割引だけが 25bp の曲線に変わる。
    let ratio = evaluate(&foreign, &market)?.get(Currency::JPY)
        / evaluate(&native, &market)?.get(Currency::JPY);
    let pay = model_time(as_of, at("2026-01-14", "16:00")?).years();
    assert_close!(ratio, (-0.0025 * pay).exp(), 1e-12_f64, 1e-14_f64);
    let strict = fx_config(&format!(
        "{rules}[conventions]\ndisabled = ['cross_collateral_projection']"
    ))?;
    assert!(matches!(
        prepare_under(&strict, &collateral),
        Err(Error::Unsupported(_))
    ));
    Ok(())
}

#[test]
fn futures_margin_index_cross_collateral_projection_is_recorded() -> Result<()> {
    // 値洗いの口座の付利指標（margin_index）の予想が値洗いの割引と等しい率なら、値洗いの列の評価は
    // 厳密（Q_M の測度と割引が一致）。付利指標を別の曲線で予想すると、その差を無視する近似として
    // CrossCollateralProjection を記録する。
    let as_of = AS_OF.parse()?;
    let mut parts = reference_parts("0D".parse()?, "16:00:00")?;
    let mut margin = parts.indices[&"TERM".parse()?].clone();
    margin.id = "MARGIN".parse()?;
    margin.kind = IndexKind::OvernightRate {
        ccy: Currency::USD,
        day_count: DayCount::Act360,
    };
    parts.indices.insert(margin.id.clone(), margin);
    if let Some(Convention::RateFuture(convention)) =
        parts.conventions.get_mut(&"USD-RATE-FUTURES".parse()?)
    {
        convention.margin_index = Some("MARGIN".parse()?);
    }
    let reference = RefData::new(parts)?;
    let alternative = CurveSpec::new(
        "USD-ALT".parse()?,
        Currency::USD,
        [CurveRole::Projection],
        [PILLAR.parse()?],
        Extrapolation::Reject,
    )?;
    let structure = Arc::new(MarketStructure::with_fx_and_lgm(
        as_of,
        [curve()?, alternative],
        [],
        [(Currency::USD, vec![])],
    )?);
    let mut theta = vec![0.0; structure.layout().len()];
    let t = model_time(as_of, PILLAR.parse()?).years();
    for curve in ["USD", "USD-ALT"] {
        let range = structure
            .layout()
            .range(&BlockId::nodes(curve.parse()?))
            .ok_or_else(|| Error::Internal("missing curve block".into()))?;
        theta[range.start] = -RATE * t;
    }
    let market = Market::from_theta(structure, theta)?;
    let fixings = FixingStore::default();
    let warnings = |margin_curve: &str| -> Result<Vec<Warning>> {
        let config = config(&format!("MARGIN = \"{margin_curve}\""))?;
        let ctx = Context {
            as_of,
            refdata: &reference,
            config: &config,
            fixings: &fixings,
        };
        Ok(futures_pv(&future(Position::Long, 0.95)?, &ctx, &market)?.1)
    };
    assert!(!warnings("USD")?.contains(&Warning::CrossCollateralProjection));
    assert!(warnings("USD-ALT")?.contains(&Warning::CrossCollateralProjection));
    Ok(())
}
`````

`````
crates/rhodium-pricing/tests/ois_calibration.rs
//! 標準 OIS の校正を `QuantLib` と照合し、上流への線形化を検証する。

use std::{ptr, sync::Arc};

use rhodium_core::{
    Error, Result, Warning,
    currency::Currency,
    money::Money,
    observation::{ObservationMethod, ObservationRounding},
    quantity::Spread,
    rounding::{Rounding, RoundingMode},
    schedule::EndRule,
};
use rhodium_market::{
    curve::CurveId,
    market::Market,
    quote::{QuoteId, QuoteOverride, QuoteOverrides, QuoteSet},
    structure::{MarketStructure, ParamId},
    view::Slots,
};
use rhodium_pricing::{
    calibration::{CalibratedMarket, GroupOutcome, calibrate, calibrate_with_residual_tolerance},
    context::Context,
    definition::{ElementSelector, MarketDefinition},
    prepare::{PrepareInput, evaluate, prepare},
    resolve::CollateralKey,
};
use rhodium_products::{
    agreement::{CashCollateralTerms, Remuneration},
    standard::{PaymentRounding, StandardProduct, Swap},
};
use rhodium_testkit::assert_close;
use scioxide_core::AdjointOperator;

#[path = "fixtures/ois.rs"]
pub mod ois;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、このテストで使わない関数がある"
)]
#[path = "fixtures/multicurve.rs"]
pub mod multicurve;

use ois::{Case, Inputs, cases, definition, fixture_quotes, midnight, prepare_ois};

fn par_rates(
    case: &Case,
    collateral: &CollateralKey,
    market: &Market<f64>,
    ctx: &Context<'_>,
) -> Result<Vec<f64>> {
    case.nodes
        .iter()
        .map(|node| {
            let value = |rate| {
                let plan = prepare_ois(case, node, rate, collateral, market, ctx)?;
                Ok(evaluate(&plan, market)?.get(case.currency))
            };
            let zero = value(0.0)?;
            let one: f64 = value(1.0)?;
            Ok(-zero / (one - zero))
        })
        .collect()
}

fn assert_quantlib(case: &Case) -> Result<()> {
    let collateral = CollateralKey::Uncollateralized;
    let inputs = Inputs::new(&[case], case.as_of, "")?;
    let ctx = inputs.context();
    let result = calibrate(
        &definition(&[(case, &collateral)])?,
        &QuoteSet::new(fixture_quotes(case)?)?,
        &ctx,
    )?;
    let group = &result.groups()[0];
    let GroupOutcome::Solved(solution) = group.outcome() else {
        return Err(Error::Calibration(format!("{group:?}")));
    };
    assert!(solution.residual_norm() <= 1e-14);
    let market = result.market();
    let id: CurveId = case.index.parse()?;
    let expected: Vec<_> = case
        .nodes
        .iter()
        .map(|node| (node.date, node.df))
        .chain(case.samples.iter().map(|sample| (sample.date, sample.df)))
        .collect();
    let points = expected
        .iter()
        .map(|&(date, _)| market.structure().prepare_point(&id, midnight(date)?))
        .collect::<Result<Vec<_>>>()?;
    let slots = Slots::new(Arc::clone(market.structure()), [(id, points)], [])?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        for (&point, &(_, df)) in bound.discounts()[0].points().iter().zip(&expected) {
            assert_close!(view.df(point).value(), df, 1e-12, 1e-14);
        }
        Ok(())
    })?;
    // 校正商品を取引と同じ準備で評価し、quote の固定金利で PV がほぼ 0 になることを確かめる。
    // 元本 1 の PV は年金現価 A を使って A·(q − par) なので、品質の「ソルバーの誤差」の
    // c_s·|∂PV/∂q|·τ = 2·A·τ（τ は残差の許容誤差 1e-14）で比べる。
    for (node, rate) in case
        .nodes
        .iter()
        .zip(par_rates(case, &collateral, market, &ctx)?)
    {
        assert_close!(rate, node.rate, 0.0, 1e-14);
        let value = |fixed_rate| -> Result<f64> {
            let pricing = prepare_ois(case, node, fixed_rate, &collateral, market, &ctx)?;
            Ok(evaluate(&pricing, market)?.get(case.currency))
        };
        let annuity = (value(1.0)? - value(0.0)?).abs();
        assert!(annuity > 0.0);
        let pv = value(node.rate)?;
        assert!(
            pv.abs() <= 2.0 * annuity * 1e-14,
            "{} PV at the quote {pv:e}, annuity {annuity}",
            node.tenor
        );
    }
    Ok(())
}

#[test]
fn tona_ois_calibration_matches_quantlib() -> Result<()> {
    // 負金利も含む TONA の quote を再現し、ログ線形補間の DF を独立した解と照合する。
    assert_quantlib(&cases()?[0])
}

#[test]
fn sofr_ois_calibration_matches_quantlib() -> Result<()> {
    // SOFR のピラーは満期から 2 営業日後の支払日。割引にも予想にも同じ曲線を使う。
    assert_quantlib(&cases()?[1])
}

#[test]
fn multicurve_calibration_matches_quantlib() -> Result<()> {
    // OIS を先に解き、OIS で割り引くベーシス商品でターム指標の予想を解く。
    // JPY の TIBOR 6M と 3M、USD のテスト用ターム指標（6M）の三つを QuantLib と照合する。
    for case in multicurve::all_multicurve_cases()? {
        let inputs = case.inputs()?;
        let ctx = inputs.context();
        let result = calibrate(&case.definition()?, &case.quotes()?, &ctx)?;
        let [ois, basis] = result.groups() else {
            return Err(Error::Internal("expected OIS and basis groups".into()));
        };
        assert_eq!(ois.blocks(), &[case.overnight.block()?]);
        assert_eq!(basis.blocks(), &[case.term_block()?]);
        assert_eq!(basis.upstream(), ois.blocks());
        let tolerance = 1e-14;
        for group in [ois, basis] {
            let GroupOutcome::Solved(solution) = group.outcome() else {
                return Err(Error::Calibration(format!("{group:?}")));
            };
            assert!(solution.residual_norm() <= tolerance);
        }

        let market = result.market();
        let collateral = CollateralKey::Uncollateralized;
        let ccy = case.overnight.currency;
        let ois_quote =
            |node, rate| prepare_ois(&case.overnight, node, rate, &collateral, market, &ctx);
        let basis_quote = |node, rate| case.prepare_basis(node, rate, market, &ctx);
        for (nodes, prepare_quote) in [
            (&case.overnight.nodes, &ois_quote as &dyn Fn(_, _) -> _),
            (&case.nodes, &basis_quote),
        ] {
            for node in nodes {
                let pv =
                    |rate| Ok::<_, Error>(evaluate(&prepare_quote(node, rate)?, market)?.get(ccy));
                let zero = pv(0.0)?;
                let slope = pv(1.0)? - zero;
                assert!(slope.abs() > 0.0);
                assert_close!(-zero / slope, node.rate, 0.0, tolerance);
                // 残差の許容誤差をスプレッドまたは固定金利 1 あたりの PV で伝える。
                assert_close!(pv(node.rate)?, 0.0, 0.0, 2.0 * slope.abs() * tolerance);
            }
        }

        let structure = market.structure();
        let discount: CurveId = case.overnight.index.parse()?;
        let term: CurveId = case.term_index.parse()?;
        let expected_discount = case
            .overnight
            .nodes
            .iter()
            .map(|n| (n.date, n.df))
            .chain(case.overnight.samples.iter().map(|s| (s.date, s.df)))
            .collect::<Vec<_>>();
        let expected_term = case
            .nodes
            .iter()
            .map(|n| (n.date, n.df))
            .chain(case.samples.iter().map(|s| (s.date, s.df)))
            .collect::<Vec<_>>();
        let discount_points = expected_discount
            .iter()
            .map(|&(date, _)| structure.prepare_point(&discount, midnight(date)?))
            .collect::<Result<Vec<_>>>()?;
        // 予想曲線の擬似 DF は、評価時点からの growth の逆数で読む。
        let origin = structure.prepare_point(&term, ctx.as_of)?;
        let term_points = expected_term
            .iter()
            .map(|&(date, _)| {
                Ok((
                    origin.clone(),
                    structure.prepare_point(&term, midnight(date)?)?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        let slots = Slots::new(Arc::clone(structure), [(discount, discount_points)], [])?
            .with_projection_growth([(structure.resolve_projection(&term)?, term_points)])?;
        market.with_view(|view| {
            let bound = view.bind_slots(&slots)?;
            for (&point, &(_, df)) in bound.discounts()[0].points().iter().zip(&expected_discount) {
                assert_close!(view.df(point).value(), df, 1e-12, 1e-14);
            }
            for (&point, &(_, df)) in bound.projection_growth()[0]
                .points()
                .iter()
                .zip(&expected_term)
            {
                assert_close!(view.growth(point).recip(), df, 1e-12, 1e-14);
            }
            Ok(())
        })?;
    }
    Ok(())
}

fn with_theta(result: &CalibratedMarket, theta: Vec<f64>) -> Result<Market<f64>> {
    Market::from_theta(Arc::clone(result.market().structure()), theta)
}

fn bumped_rates(
    result: &CalibratedMarket,
    index: usize,
    case: &Case,
    ctx: &Context<'_>,
) -> Result<Vec<f64>> {
    let h = 1e-6;
    let rates = |bump| {
        let mut theta = result.market().theta().to_vec();
        theta[index] += bump;
        par_rates(
            case,
            &CollateralKey::Uncollateralized,
            &with_theta(result, theta)?,
            ctx,
        )
    };
    let up = rates(h)?;
    let down = rates(-h)?;
    Ok(up
        .into_iter()
        .zip(down)
        .map(|(up, down)| (up - down) / (2.0 * h))
        .collect())
}

#[test]
fn upstream_jacobian_and_transposed_lu_match_bumps() -> Result<()> {
    // 曲線 ID の順とは逆に、割引曲線を先に解いて下流の予想曲線に伝える。
    let [tona, _] = cases()?;
    let mut forecast = tona.clone();
    forecast.index = "JPY-FORECAST".into();
    for node in &mut forecast.nodes {
        node.rate += 0.002;
    }
    let uncollateralized = CollateralKey::Uncollateralized;
    let inputs = Inputs::new(&[&tona, &forecast], tona.as_of, "")?;
    let ctx = inputs.context();
    let def = definition(&[(&forecast, &uncollateralized), (&tona, &uncollateralized)])?;
    let result = calibrate(
        &def,
        &QuoteSet::new(
            fixture_quotes(&tona)?
                .into_iter()
                .chain(fixture_quotes(&forecast)?),
        )?,
        &ctx,
    )?;
    let [upstream, downstream] = result.groups() else {
        return Err(Error::Internal("expected two groups".into()));
    };
    assert_eq!(upstream.blocks(), &[tona.block()?]);
    assert_eq!(downstream.blocks(), &[forecast.block()?]);
    assert_eq!(downstream.upstream(), upstream.blocks());
    let GroupOutcome::Solved(solution) = downstream.outcome() else {
        return Err(Error::Calibration(format!("{downstream:?}")));
    };
    let layout = result.market().structure().layout();
    let upstream_indices = layout
        .range(&tona.block()?)
        .ok_or_else(|| Error::Internal("no upstream block".into()))?;
    let own_indices = layout
        .range(&forecast.block()?)
        .ok_or_else(|| Error::Internal("no solved block".into()))?;
    let n = own_indices.len();
    let rhs = vec![1.0; n];
    let mut lambda = vec![0.0; n];
    solution
        .dg_dtheta()
        .inverse_operator()
        .apply_adjoint_into(&rhs, &mut lambda, &mut vec![0.0; n])
        .map_err(|error| Error::Numerical(error.to_string()))?;
    for (col, index) in own_indices.enumerate() {
        let derivatives = bumped_rates(&result, index, &forecast, &ctx)?;
        let got: f64 = derivatives
            .iter()
            .zip(&lambda)
            .map(|(derivative, weight)| derivative * weight)
            .sum();
        assert_close!(got, rhs[col], 1e-8, 1e-9);
    }
    assert!(solution.dg_dup().nnz() > 0);
    for (col, index) in upstream_indices.enumerate() {
        for (row, expected) in bumped_rates(&result, index, &forecast, &ctx)?
            .into_iter()
            .enumerate()
        {
            let actual = solution
                .dg_dup()
                .get(row, col)
                .map_err(|error| Error::Numerical(error.to_string()))?;
            assert_close!(actual, expected, 1e-8, 1e-10);
        }
    }
    assert_eq!(solution.dg_dq().nnz(), n);
    for (row, col, value) in solution.dg_dq().iter() {
        assert_eq!(&downstream.row_quotes()[row], &downstream.quotes()[col]);
        assert_close!(value, -1.0, 0.0, 1e-15);
    }
    Ok(())
}

#[test]
fn cyclic_discounts_and_forecasts_calibrate_together() -> Result<()> {
    // 二つの OIS が互いの曲線で割り引くときは、全節点を含む一つの正方系になる。
    let [tona, _] = cases()?;
    let mut second = tona.clone();
    second.index = "JPY-SECOND".into();
    for node in &mut second.nodes {
        node.rate += 0.003;
    }
    let cash = CollateralKey::Cash(CashCollateralTerms {
        ccy: Currency::JPY,
        remuneration: Remuneration {
            index: tona.index.parse()?,
            spread: Spread::new(0.0)?,
            floor: None,
        },
    });
    let uncollateralized = CollateralKey::Uncollateralized;
    let inputs = Inputs::new(
        &[&tona, &second],
        tona.as_of,
        r#"
[[discounting]]
currency = "JPY"
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
curve = "JPY-SECOND"
"#,
    )?;
    let ctx = inputs.context();
    let def = definition(&[(&tona, &cash), (&second, &uncollateralized)])?;
    let result = calibrate(
        &def,
        &QuoteSet::new(
            fixture_quotes(&tona)?
                .into_iter()
                .chain(fixture_quotes(&second)?),
        )?,
        &ctx,
    )?;
    let [group] = result.groups() else {
        return Err(Error::Internal("expected one cyclic group".into()));
    };
    assert_eq!(group.blocks(), &[second.block()?, tona.block()?]);
    assert_eq!(group.upstream(), []);
    let GroupOutcome::Solved(solution) = group.outcome() else {
        return Err(Error::Calibration(format!("{group:?}")));
    };
    assert_eq!(
        solution.dg_dtheta().order(),
        tona.nodes.len() + second.nodes.len()
    );
    assert!(solution.residual_norm() <= 1e-14);

    // 前の θ のままで許容誤差を満たす変化でも、循環した系は全体校正と同じ固定初期値を使う。
    let overrides = QuoteOverrides::new(
        [QuoteOverride::Set {
            id: tona.quote(&tona.nodes[0])?,
            value: tona.nodes[0].rate + 2e-16,
        }],
        def.quote_kinds(),
    )?;
    let (bumped, _) = result.quotes().with_overrides(&overrides)?;
    let updated = result.update(&bumped)?;
    let full = calibrate(&def, &bumped, &ctx)?;
    let (GroupOutcome::Solved(updated_solution), GroupOutcome::Solved(full_solution)) =
        (updated.groups()[0].outcome(), full.groups()[0].outcome())
    else {
        return Err(Error::Calibration(
            "cyclic group did not solve after bump".into(),
        ));
    };
    assert!(full_solution.iterations() > 0);
    assert_eq!(updated_solution.iterations(), full_solution.iterations());
    assert_market_parameters_match(&updated, &full);
    Ok(())
}

#[test]
fn missing_quotes_leave_only_dependents_unavailable() -> Result<()> {
    // TONA の quote が揃わない場合も、依存しない SOFR は通常どおり利用できる。
    let [tona, sofr] = cases()?;
    let mut forecast = tona.clone();
    forecast.index = "JPY-FORECAST".into();
    let collateral = CollateralKey::Uncollateralized;
    let inputs = Inputs::new(&[&tona, &forecast, &sofr], tona.as_of, "")?;
    let def = definition(&[
        (&forecast, &collateral),
        (&sofr, &collateral),
        (&tona, &collateral),
    ])?;
    let result = calibrate(
        &def,
        &QuoteSet::new(
            fixture_quotes(&forecast)?
                .into_iter()
                .chain(fixture_quotes(&sofr)?),
        )?,
        &inputs.context(),
    )?;
    let [failed, blocked, independent] = result.groups() else {
        return Err(Error::Internal("expected three groups".into()));
    };
    let GroupOutcome::Failed(Error::MissingData(message)) = failed.outcome() else {
        return Err(Error::Calibration(format!("{failed:?}")));
    };
    assert!(
        failed
            .quotes()
            .iter()
            .all(|id| message.contains(id.as_str()))
    );
    assert!(matches!(blocked.outcome(), GroupOutcome::Unavailable));
    assert!(matches!(independent.outcome(), GroupOutcome::Solved(_)));
    assert!(!result.market().block_is_available(&tona.block()?));
    assert!(!result.market().block_is_available(&forecast.block()?));
    assert!(result.market().block_is_available(&sofr.block()?));
    Ok(())
}

#[test]
fn quote_update_recalibrates_transitive_dependents_and_reuses_other_solutions() -> Result<()> {
    // TONA → FORECAST → TAIL と、独立した SOFR を構成する。TAIL は変更した quote を直接読まない。
    let [tona, sofr] = cases()?;
    let mut forecast = tona.clone();
    forecast.index = "JPY-FORECAST".into();
    for node in &mut forecast.nodes {
        node.rate += 0.002;
    }
    let mut tail = forecast.clone();
    tail.index = "JPY-TAIL".into();
    let cash = CollateralKey::Cash(CashCollateralTerms {
        ccy: Currency::JPY,
        remuneration: Remuneration {
            index: tona.index.parse()?,
            spread: Spread::new(0.0)?,
            floor: None,
        },
    });
    let uncollateralized = CollateralKey::Uncollateralized;
    let inputs = Inputs::new(
        &[&tona, &forecast, &tail, &sofr],
        tona.as_of,
        r#"
[[discounting]]
currency = "JPY"
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
curve = "JPY-FORECAST"
"#,
    )?;
    let ctx = inputs.context();
    let def = definition(&[
        (&tail, &cash),
        (&forecast, &uncollateralized),
        (&tona, &uncollateralized),
        (&sofr, &uncollateralized),
    ])?;
    let entries = [&tona, &forecast, &tail, &sofr]
        .into_iter()
        .map(fixture_quotes)
        .collect::<Result<Vec<_>>>()?;
    let quotes = QuoteSet::new(entries.into_iter().flatten())?;
    let base = calibrate(&def, &quotes, &ctx)?;
    let [upstream, middle, downstream, independent] = base.groups() else {
        return Err(Error::Internal("expected four groups".into()));
    };
    assert_eq!(upstream.blocks(), &[tona.block()?]);
    assert_eq!(middle.upstream(), upstream.blocks());
    assert_eq!(downstream.upstream(), middle.blocks());
    assert_eq!(independent.blocks(), &[sofr.block()?]);

    let bumped_id = upstream.quotes()[0].clone();
    let overrides = QuoteOverrides::new(
        [QuoteOverride::Set {
            id: bumped_id.clone(),
            value: 0.009,
        }],
        def.quote_kinds(),
    )?;
    let (bumped, _) = quotes.with_overrides(&overrides)?;
    let updated = base.update(&bumped)?;
    let full = calibrate(&def, &bumped, &ctx)?;
    assert!(Arc::ptr_eq(
        base.market().structure(),
        updated.market().structure()
    ));
    for (index, (before, after)) in base.groups().iter().zip(updated.groups()).enumerate() {
        assert!(matches!(before.outcome(), GroupOutcome::Solved(_)));
        assert!(matches!(after.outcome(), GroupOutcome::Solved(_)));
        assert_eq!(ptr::eq(before.outcome(), after.outcome()), index == 3);
    }
    assert_market_parameters_match(&updated, &full);

    // 削除・復元で 2 段の下流まで利用可否が更新され、独立した解は同じ分解を持ち続ける。
    let remove = QuoteOverrides::new([QuoteOverride::Remove { id: bumped_id }], def.quote_kinds())?;
    let (missing, _) = bumped.with_overrides(&remove)?;
    let partial = updated.update(&missing)?;
    assert!(matches!(
        partial.groups()[0].outcome(),
        GroupOutcome::Failed(Error::MissingData(_))
    ));
    assert!(matches!(
        partial.groups()[1].outcome(),
        GroupOutcome::Unavailable
    ));
    assert!(matches!(
        partial.groups()[2].outcome(),
        GroupOutcome::Unavailable
    ));
    assert!(ptr::eq(
        partial.groups()[3].outcome(),
        independent.outcome()
    ));
    for block in [tona.block()?, forecast.block()?, tail.block()?] {
        assert!(!partial.market().block_is_available(&block));
    }
    let restored = partial.update(&bumped)?;
    assert_market_parameters_match(&restored, &full);
    assert!(ptr::eq(
        restored.groups()[3].outcome(),
        independent.outcome()
    ));
    Ok(())
}

fn assert_market_parameters_match(actual: &CalibratedMarket, expected: &CalibratedMarket) {
    assert_eq!(
        actual.market().theta().len(),
        expected.market().theta().len()
    );
    for (&actual, &expected) in actual
        .market()
        .theta()
        .iter()
        .zip(expected.market().theta())
    {
        assert_close!(actual, expected, 1e-12, 1e-14);
    }
}

#[test]
fn quote_update_reuses_the_market_for_unused_quotes_and_refreshes_unit_warnings() -> Result<()> {
    // 未使用の quote を受け取っても、構造と θ と解の記録を共有する。
    let [tona, _] = cases()?;
    let collateral = CollateralKey::Uncollateralized;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let def = definition(&[(&tona, &collateral)])?;
    let quotes = QuoteSet::new(fixture_quotes(&tona)?)?;
    let base = calibrate(&def, &quotes, &inputs.context())?;
    let unused: QuoteId = "JPY.UNUSED".parse()?;
    let mixed = QuoteSet::new(
        quotes
            .iter()
            .map(|(id, value)| (id.clone(), value))
            .chain([(unused.clone(), 4.0)]),
    )?;
    let unchanged = base.update(&mixed)?;
    assert_eq!(unchanged.quotes().get(&unused), Some(4.0));
    assert!(ptr::eq(base.market().theta(), unchanged.market().theta()));
    assert!(ptr::eq(
        base.groups()[0].outcome(),
        unchanged.groups()[0].outcome()
    ));
    assert_eq!(base.warnings(), unchanged.warnings());

    let rate = base.groups()[0].quotes()[0].clone();
    let warning = Warning::QuoteUnitSuspected {
        quote: rate.clone(),
    };
    let suspicious = QuoteOverrides::new(
        [QuoteOverride::Set {
            id: rate,
            value: 1.01,
        }],
        def.quote_kinds(),
    )?;
    let (suspicious, _) = quotes.with_overrides(&suspicious)?;
    let updated = unchanged.update(&suspicious)?;
    assert!(updated.warnings().contains(&warning));
    assert_eq!(
        updated.warnings(),
        calibrate(&def, &suspicious, &inputs.context())?.warnings()
    );
    assert_eq!(updated.update(&quotes)?.warnings(), base.warnings());
    Ok(())
}

#[test]
fn quote_update_starts_a_monotone_ois_group_from_the_previous_theta() -> Result<()> {
    // 残差の許容誤差の内側にある小さな変化なら、前の θ は初期点で収束する。
    // 固定初期値では反復を要する 2 本以上のピラーを使い、短期と複数クーポンをともに含める。
    let [tona, _] = cases()?;
    let collateral = CollateralKey::Uncollateralized;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let ctx = inputs.context();
    let def = definition(&[(&tona, &collateral)])?;
    let quotes = QuoteSet::new(fixture_quotes(&tona)?)?;
    let base = calibrate(&def, &quotes, &ctx)?;
    let id = tona.quote(&tona.nodes[0])?;
    let overrides = QuoteOverrides::new(
        [QuoteOverride::Set {
            id,
            value: tona.nodes[0].rate + 2e-16,
        }],
        def.quote_kinds(),
    )?;
    let (bumped, _) = quotes.with_overrides(&overrides)?;
    let updated = base.update(&bumped)?;
    let full = calibrate(&def, &bumped, &ctx)?;
    let GroupOutcome::Solved(warm) = updated.groups()[0].outcome() else {
        return Err(Error::Calibration(format!(
            "{:?}",
            updated.groups()[0].outcome()
        )));
    };
    let GroupOutcome::Solved(cold) = full.groups()[0].outcome() else {
        return Err(Error::Calibration(format!(
            "{:?}",
            full.groups()[0].outcome()
        )));
    };
    assert_eq!(warm.iterations(), 0);
    assert!(cold.iterations() > 0);
    assert_market_parameters_match(&updated, &full);
    Ok(())
}

#[test]
fn singular_residual_jacobian_is_not_a_solved_curve() -> Result<()> {
    // すべての行を同じ 1M OIS にすると残差がゼロでも残りの節点を識別できない。
    let [tona, sofr] = cases()?;
    let collateral = CollateralKey::Uncollateralized;
    let inputs = Inputs::new(&[&tona, &sofr], tona.as_of, "")?;
    let mut repeated = tona.clone();
    for node in &mut repeated.nodes {
        node.tenor = tona.nodes[0].tenor;
        node.rate = 0.0;
    }
    let quotes = QuoteSet::new(
        [(repeated.quote(&repeated.nodes[0])?, 0.0)]
            .into_iter()
            .chain(fixture_quotes(&sofr)?),
    )?;
    let result = calibrate(
        &definition(&[(&repeated, &collateral), (&sofr, &collateral)])?,
        &quotes,
        &inputs.context(),
    )?;
    assert!(matches!(
        result.groups()[0].outcome(),
        GroupOutcome::Failed(Error::Calibration(_))
    ));
    assert!(matches!(
        result.groups()[1].outcome(),
        GroupOutcome::Solved(_)
    ));
    assert!(!result.market().block_is_available(&tona.block()?));
    Ok(())
}

#[test]
fn definition_rejects_unassigned_and_out_of_range_nodes() -> Result<()> {
    // 指定先が存在しない要素や、どの残差も当てはめない要素を受理しない。
    let [tona, _] = cases()?;
    let collateral = CollateralKey::Uncollateralized;
    for element in [0, 99] {
        let binding = tona.binding(
            &tona.nodes[0],
            &collateral,
            ElementSelector::Param(ParamId::new(tona.block()?, element)),
        )?;
        let result = MarketDefinition::new(
            "JPY".parse()?,
            "Etc/UTC".parse()?,
            [tona.curve()?],
            [binding],
        );
        assert!(matches!(result, Err(Error::InvalidInput(_))));
    }
    Ok(())
}

#[test]
fn calibration_rejects_a_non_square_group() -> Result<()> {
    // 全節点を指定していても、残差が一つでは五つの節点を解けない。
    let [tona, _] = cases()?;
    let collateral = CollateralKey::Uncollateralized;
    let node = &tona.nodes[0];
    let binding = tona.binding(node, &collateral, ElementSelector::Block(tona.block()?))?;
    let def = MarketDefinition::new(
        "JPY".parse()?,
        "Etc/UTC".parse()?,
        [tona.curve()?],
        [binding],
    )?;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let quotes = QuoteSet::new([(tona.quote(node)?, node.rate)])?;
    assert!(matches!(
        calibrate(&def, &quotes, &inputs.context()),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn calibration_rejects_an_invalid_residual_tolerance() -> Result<()> {
    // 残差の許容誤差は入口で一度だけ解析し、有限な非負値でなければ求根を始めない。
    let [tona, _] = cases()?;
    let collateral = CollateralKey::Uncollateralized;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let def = definition(&[(&tona, &collateral)])?;
    let quotes = QuoteSet::new(fixture_quotes(&tona)?)?;
    for tolerance in [-1e-14, f64::NAN, f64::INFINITY] {
        assert!(
            matches!(
                calibrate_with_residual_tolerance(&def, &quotes, &inputs.context(), tolerance),
                Err(Error::InvalidInput(_))
            ),
            "tolerance {tolerance}"
        );
    }
    Ok(())
}

#[test]
fn calibration_uses_the_local_business_trade_date() -> Result<()> {
    // UTC 金曜の夕方は東京の土曜。月曜も休場なので火曜約定の OIS を評価する。
    let [mut tona, _] = cases()?;
    tona.as_of = "2025-01-14".parse()?;
    tona.nodes.truncate(1);
    tona.nodes[0].date = "2025-03-03".parse()?;
    let mut inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    inputs.as_of = "2025-01-10T18:00:00Z".parse()?;
    let ctx = inputs.context();
    let known_theta = -0.0005;
    let market = Market::from_theta(
        Arc::new(MarketStructure::new(ctx.as_of, [tona.curve()?])?),
        vec![known_theta],
    )?;
    let collateral = CollateralKey::Uncollateralized;
    let node = &tona.nodes[0];
    let [par] = par_rates(&tona, &collateral, &market, &ctx)?[..] else {
        return Err(Error::Internal("expected one quote".into()));
    };
    let binding = tona.binding(node, &collateral, ElementSelector::Block(tona.block()?))?;
    let def = MarketDefinition::new(
        "JPY".parse()?,
        "Asia/Tokyo".parse()?,
        [tona.curve()?],
        [binding],
    )?;
    let result = calibrate(&def, &QuoteSet::new([(tona.quote(node)?, par)])?, &ctx)?;
    assert!(matches!(
        result.groups()[0].outcome(),
        GroupOutcome::Solved(_)
    ));
    assert_close!(result.market().theta()[0], known_theta, 1e-12, 1e-14);
    Ok(())
}

// QuantLib のノードの割引係数を θ にそのまま入れ、曲線の校正の差を除いて、
// 慣行からの日程（スポットラグ、満期からの後ろ向きの生成、支払いの遅れ）、確定した fixing、
// 複利の外のスプレッドを含むスワップの評価だけを QuantLib と照合する。
fn assert_quantlib_swaps(case: &Case) -> Result<()> {
    let collateral = CollateralKey::Uncollateralized;
    let inputs = Inputs::new(&[case], case.as_of, "")?;
    let ctx = inputs.context();
    let structure = Arc::new(MarketStructure::new(ctx.as_of, [case.curve()?])?);
    let theta = case.nodes.iter().map(|node| node.df.ln()).collect();
    let market = Market::from_theta(structure, theta)?;
    let minor_units = ctx.refdata.currency(case.currency)?.minor_units;
    let absolute = 1e-6 * 10.0_f64.powi(-i32::from(minor_units));
    assert_eq!(case.swaps.len(), 2);
    for swap in &case.swaps {
        // スポット起点の取引は開始日を省き、約定日（評価日）から慣行の 2 営業日で求めさせる。
        let (start, trade_date) = match swap.name.as_str() {
            "spot" => (None, case.as_of),
            _ => (Some(swap.start), swap.start),
        };
        let price = |fixed_rate| -> Result<f64> {
            let product = StandardProduct::Swap(Swap::from_convention(
                case.index.parse()?,
                start,
                EndRule::AfterStart(swap.tenor),
                Money::new(case.currency, swap.notional),
                fixed_rate,
                swap.fixed_side,
                Spread::new(swap.spread)?,
                PaymentRounding::None,
            )?);
            let prepared = prepare(
                PrepareInput::Standard {
                    product: &product,
                    trade_date,
                },
                &collateral,
                &[],
                market.structure(),
                &ctx,
            )?;
            Ok(evaluate(&prepared, &market)?.get(case.currency))
        };
        let pv = price(swap.fixed_rate)?;
        // 固定金利 0 の同じ取引は変動脚だけの PV で、差が固定脚の PV。
        let floating = price(0.0)?;
        // 変動脚は元本ほどの大きさの二つの項（N·P(s)/P(e) と N）の差なので、
        // 品質の許容誤差の和の項の絶対値に、両脚の PV の絶対値と元本の 2 倍を入れる。
        let scale = swap.fixed_pv.abs() + swap.floating_pv.abs() + 2.0 * swap.notional;
        assert_close!(pv, swap.pv, 1e-12, absolute, scale);
        assert_close!(floating, swap.floating_pv, 1e-12, absolute, scale);
        assert_close!(pv - floating, swap.fixed_pv, 1e-12, absolute, scale);
    }
    Ok(())
}

#[test]
fn tona_convention_ois_pv_matches_quantlib() -> Result<()> {
    // TONA のスポット起点の off-market OIS と、fixing とスプレッドを持つ開始済みの OIS。
    assert_quantlib_swaps(&cases()?[0])
}

#[test]
fn sofr_convention_ois_pv_matches_quantlib() -> Result<()> {
    // SOFR は支払いを観測期間の末尾の 2 営業日後に置く。
    assert_quantlib_swaps(&cases()?[1])
}

#[test]
fn short_ois_calibrate_to_the_default_tolerance() -> Result<()> {
    // 1D・1W・1M の TONA の OIS だけの曲線を、既定の残差の許容誤差 1e-14 で校正する。
    // 1D の期間（τ = 1/365）の DF 比を丸めてから 1 を引くと、パーは約 8e-14 刻みになり
    // 1e-14 に届かない。超過成長を対数から exp_m1 で求めれば、quote が 0.5%、0 付近、
    // 負のどれでも、パーが quote を 1e-14 以内で再現する（短期の桁落ちの回帰テスト）。
    // 評価日 2025-01-07（火）のスポットは 01-09 で、ピラーは各 OIS の最終支払日。
    let [tona, _] = cases()?;
    let collateral = CollateralKey::Uncollateralized;
    let node = |tenor: &str, rate: f64, date: &str| -> Result<ois::Node> {
        Ok(ois::Node {
            tenor: tenor.parse()?,
            rate,
            date: date.parse()?,
            df: f64::NAN,
        })
    };
    for rates in [
        [0.005, 0.0051, 0.0052],
        [1e-7, 0.0, -1e-7],
        [-0.0005, -0.0004, -0.0003],
    ] {
        // 参照データ（暦、支払ラグ）は QuantLib のケースのまま、ピラーだけを短期にする。
        let mut case = tona.clone();
        case.nodes = vec![
            node("1D", rates[0], "2025-01-10")?,
            node("1W", rates[1], "2025-01-16")?,
            node("1M", rates[2], "2025-02-10")?,
        ];
        case.samples.clear();
        let inputs = Inputs::new(&[&case], case.as_of, "")?;
        let ctx = inputs.context();
        let result = calibrate(
            &definition(&[(&case, &collateral)])?,
            &QuoteSet::new(fixture_quotes(&case)?)?,
            &ctx,
        )?;
        let group = &result.groups()[0];
        let GroupOutcome::Solved(solution) = group.outcome() else {
            return Err(Error::Calibration(format!("{rates:?}: {group:?}")));
        };
        assert!(solution.residual_norm() <= 1e-14);
        for (node, rate) in
            case.nodes
                .iter()
                .zip(par_rates(&case, &collateral, result.market(), &ctx)?)
        {
            assert_close!(rate, node.rate, 0.0, 1e-14);
        }
    }
    Ok(())
}

#[test]
fn rate_quotes_above_one_warn_about_a_missing_unit() -> Result<()> {
    // 金利の binding が読む quote の絶対値が 1 を超えたら、校正の結果に quote ID 付きで警告する。
    // 5% を単位なしの 5 と書いたような誤りは、値だけでは 500% として校正されてしまうため。
    // 0.05（5%）のような通常の値と、絶対値が 1 ちょうどの値は警告しない。
    let [tona, _] = cases()?;
    let collateral = CollateralKey::Uncollateralized;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let ctx = inputs.context();
    let def = definition(&[(&tona, &collateral)])?;
    let first = tona.quote(&tona.nodes[0])?;
    let second = tona.quote(&tona.nodes[1])?;
    for (values, suspected) in [
        ([5.0, 0.05], vec![first.clone()]),
        ([0.05, -1.5], vec![second.clone()]),
        ([0.05, 1.0], vec![]),
    ] {
        let quotes = fixture_quotes(&tona)?.into_iter().map(|(id, value)| {
            let value = if id == first {
                values[0]
            } else if id == second {
                values[1]
            } else {
                value
            };
            (id, value)
        });
        let result = calibrate(&def, &QuoteSet::new(quotes)?, &ctx)?;
        let warned: Vec<_> = result
            .warnings()
            .iter()
            .filter_map(|warning| match warning {
                Warning::QuoteUnitSuspected { quote } => Some(quote.clone()),
                _ => None,
            })
            .collect();
        assert_eq!(warned, suspected, "{values:?}");
    }
    Ok(())
}

#[test]
fn calibration_rejects_conventions_that_round_observed_rates() -> Result<()> {
    // TONA の集約金利を 10 bp 刻みに丸める慣行では、同じ丸め値を返す曲線が一意に決まらない
    // （45 bp の quote は残差 5 bp で解けず、50 bp は解けてデルタを返し、50.5 bp では再び解けない）。
    // 日次の丸めも集約後の丸めも、校正の入口で quote 付きの InvalidData として拒否する。
    // 対になる正常例として、丸めの無い同じ慣行は校正できる。
    let [tona, _] = cases()?;
    let collateral = CollateralKey::Uncollateralized;
    let ten_bp = Rounding::new(0.001, RoundingMode::HalfEven)?;
    for rounding in [
        ObservationRounding {
            daily: Some(ten_bp),
            aggregated: None,
        },
        ObservationRounding {
            daily: None,
            aggregated: Some(ten_bp),
        },
    ] {
        let mut rounded = tona.clone();
        rounded.observation_rounding = rounding;
        let inputs = Inputs::new(&[&rounded], rounded.as_of, "")?;
        let def = definition(&[(&rounded, &collateral)])?;
        let quote = rounded.quote(&rounded.nodes[0])?;
        let result = calibrate(
            &def,
            &QuoteSet::new(fixture_quotes(&rounded)?)?,
            &inputs.context(),
        );
        assert!(
            matches!(&result, Err(Error::InvalidData(message))
                if message.contains("rounds observed rates") && message.contains(quote.as_str())),
            "{result:?}"
        );
    }
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let result = calibrate(
        &definition(&[(&tona, &collateral)])?,
        &QuoteSet::new(fixture_quotes(&tona)?)?,
        &inputs.context(),
    )?;
    assert!(
        result
            .groups()
            .iter()
            .all(|group| matches!(group.outcome(), GroupOutcome::Solved(_)))
    );
    Ok(())
}

#[test]
fn calibration_keeps_other_groups_when_an_instrument_cannot_be_prepared() -> Result<()> {
    // 校正商品の準備の失敗は、その binding が当てはめる曲線のグループの失敗として残し、
    // 独立した通貨の曲線は校正を続ける（以前は校正全体のエラーで、JPY の結果まで消えた）。
    // 一つ目は SOFR を 5 営業日の lookback にして公表済みの fixing を欠かした場合で、スポット起点の
    // OIS も評価日より前の観測を読むので MissingData になる。二つ目は支払いの遅れの近似
    // （SOFR は 2 営業日後払い）を設定で禁止した場合で、Unsupported になる。
    let [tona, sofr] = cases()?;
    let collateral = CollateralKey::Uncollateralized;
    let mut lookback = sofr.clone();
    lookback.observation_method = ObservationMethod::Lookback(
        std::num::NonZeroU8::new(5).ok_or_else(|| Error::Internal("nonzero lookback".into()))?,
    );
    lookback.fixings.clear();
    let scenarios = [
        (lookback, "", "missing fixing"),
        (
            sofr.clone(),
            "[conventions]\ndisabled = [\"payment_delay_ignored\"]\n",
            "payment_delay_ignored",
        ),
    ];
    for (usd, conventions, reason) in scenarios {
        let inputs = Inputs::new(&[&tona, &usd], tona.as_of, conventions)?;
        let def = definition(&[(&tona, &collateral), (&usd, &collateral)])?;
        let quotes = QuoteSet::new(
            fixture_quotes(&tona)?
                .into_iter()
                .chain(fixture_quotes(&usd)?),
        )?;
        let result = calibrate(&def, &quotes, &inputs.context())?;
        let outcome = |case: &Case| -> Result<&GroupOutcome> {
            let block = case.block()?;
            result
                .groups()
                .iter()
                .find(|group| group.blocks().contains(&block))
                .map(rhodium_pricing::calibration::CalibrationGroup::outcome)
                .ok_or_else(|| Error::Internal(format!("no group for {}", case.index)))
        };
        assert!(
            matches!(outcome(&tona)?, GroupOutcome::Solved(_)),
            "{reason}"
        );
        let GroupOutcome::Failed(error) = outcome(&usd)? else {
            return Err(Error::Calibration(format!(
                "{reason}: {:?}",
                outcome(&usd)?
            )));
        };
        // 失敗には、準備できなかった校正商品の quote の ID が付く。
        let message = error.to_string();
        assert!(message.contains("USD-SOFR.OIS."), "{reason}: {message}");
        assert!(
            matches!(
                (reason, error),
                ("missing fixing", Error::MissingData(_))
                    | ("payment_delay_ignored", Error::Unsupported(_))
            ),
            "{reason}: {error:?}"
        );
        assert!(result.market().block_is_available(&tona.block()?));
        assert!(!result.market().block_is_available(&usd.block()?));
    }
    Ok(())
}

#[test]
fn calibration_reports_instruments_beyond_the_curve_as_input_errors() -> Result<()> {
    // 準備の失敗のうち、実行時のデータや設定によるもの（fixing の欠損、許されない近似）だけを
    // グループの失敗にし、入力どうしの不整合は従来どおり校正全体のエラーにする。
    // TONA の最後のピラーを 2Y の校正商品の満期より前（2026-10-09）に置くと、外挿を拒否する曲線の先を
    // 読むので、独立した SOFR が解けても、quote 付きの InvalidInput を返す（CLI では入力のエラー 2）。
    let [tona, sofr] = cases()?;
    let collateral = CollateralKey::Uncollateralized;
    let mut beyond = tona.clone();
    let last = beyond
        .nodes
        .last_mut()
        .ok_or_else(|| Error::Internal("empty TONA fixture".into()))?;
    last.date = "2026-10-09".parse()?;
    let inputs = Inputs::new(&[&beyond, &sofr], beyond.as_of, "")?;
    let def = definition(&[(&beyond, &collateral), (&sofr, &collateral)])?;
    let quotes = QuoteSet::new(
        fixture_quotes(&beyond)?
            .into_iter()
            .chain(fixture_quotes(&sofr)?),
    )?;
    let result = calibrate(&def, &quotes, &inputs.context());
    assert!(
        matches!(&result, Err(Error::InvalidInput(message))
            if message.contains("JPY-TONA.OIS.2Y") && message.contains("beyond last pillar")),
        "{result:?}"
    );
    Ok(())
}

#[derive(serde::Deserialize)]
struct QuantLibTurn {
    effective: rhodium_core::date::Date,
    spread: f64,
    fraction: f64,
    nodes: Vec<ois::Node>,
    samples: Vec<ois::Sample>,
}

fn tona_with_steps_and_turn(
    tona: &Case,
    turn: &QuantLibTurn,
) -> Result<(
    rhodium_market::curve::CurveSpec,
    Vec<rhodium_pricing::definition::InstrumentBinding>,
)> {
    use rhodium_core::quantity::YearFraction;
    use rhodium_market::{
        curve::{CurveExpr, CurveRole, CurveSpec, Extrapolation},
        quote::QuoteKind,
        structure::BlockId,
    };
    use rhodium_pricing::definition::{CalibrationTarget, InstrumentBinding, QuoteTransform};

    let id: CurveId = tona.index.parse()?;
    let curve = CurveSpec::from_expr(
        id.clone(),
        Currency::JPY,
        [CurveRole::Discount, CurveRole::Projection],
        CurveExpr::Sum(vec![
            CurveExpr::Stitch {
                short: Box::new(CurveExpr::Steps {
                    part: "meetings".into(),
                    meetings: vec![midnight(turn.nodes[0].date)?],
                }),
                long: Box::new(CurveExpr::Nodes {
                    part: "long".into(),
                    pillars: turn.nodes[2..]
                        .iter()
                        .map(|node| midnight(node.date))
                        .collect::<Result<_>>()?,
                    extrapolation: Extrapolation::Reject,
                }),
                at: midnight(turn.nodes[1].date)?,
            },
            CurveExpr::Jumps {
                part: "turn".into(),
                dates: vec![midnight(turn.effective)?],
            },
        ]),
    )?;
    let collateral = CollateralKey::Uncollateralized;
    let mut bindings = turn
        .nodes
        .iter()
        .enumerate()
        .map(|(position, node)| {
            let part = if position < 2 { "meetings" } else { "long" };
            tona.binding(
                node,
                &collateral,
                ElementSelector::Block(BlockId::new(id.clone(), part)?),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    // 年末のスプレッドは直接マークで、変換 TurnSpreadToJump が s·τ_turn の跳びにする。
    bindings.push(InstrumentBinding::new(
        "JPY.TONA.TURN".parse()?,
        QuoteKind::Spread,
        CalibrationTarget::DirectMark,
        QuoteTransform::TurnSpreadToJump(YearFraction::new(turn.fraction, tona.basis())?),
        collateral,
        ElementSelector::Param(ParamId::new(BlockId::new(id, "turn")?, 0)),
    )?);
    Ok((curve, bindings))
}

#[test]
fn meeting_steps_and_turn_reprice_ois_quotes() -> Result<()> {
    // QuantLib の独立した OIS ブートストラップと、会合日を含む同じ区間の DF を照合する。
    let [mut tona, sofr] = cases()?;
    let turn: QuantLibTurn = serde_json::from_str(include_str!("fixtures/quantlib/ois_turn.json"))
        .map_err(|error| Error::InvalidData(error.to_string()))?;
    let (curve, mut bindings) = tona_with_steps_and_turn(&tona, &turn)?;
    tona.nodes.clone_from(&turn.nodes);
    let collateral = CollateralKey::Uncollateralized;
    let sofr_block = sofr.block()?;
    for node in &sofr.nodes {
        bindings.push(sofr.binding(
            node,
            &collateral,
            ElementSelector::Block(sofr_block.clone()),
        )?);
    }
    let id = curve.id().clone();
    let definition = MarketDefinition::new(
        "JPY".parse()?,
        "Etc/UTC".parse()?,
        [curve, sofr.curve()?],
        bindings,
    )?;
    let quotes = QuoteSet::new(
        fixture_quotes(&tona)?
            .into_iter()
            .chain(fixture_quotes(&sofr)?)
            .chain([("JPY.TONA.TURN".parse()?, turn.spread)]),
    )?;
    let inputs = Inputs::new(&[&tona, &sofr], tona.as_of, "")?;
    let context = inputs.context();
    let result = calibrate(&definition, &quotes, &context)?;
    for group in result.groups() {
        let GroupOutcome::Solved(solution) = group.outcome() else {
            return Err(Error::Calibration(format!("{group:?}")));
        };
        assert!(solution.residual_norm() <= 1e-14);
    }
    let market = result.market();
    let expected = turn
        .nodes
        .iter()
        .map(|node| (node.date, node.df))
        .chain(turn.samples.iter().map(|sample| (sample.date, sample.df)))
        .collect::<Vec<_>>();
    let structure = market.structure();
    let points = expected
        .iter()
        .map(|&(date, _)| structure.prepare_point(&id, midnight(date)?))
        .collect::<Result<Vec<_>>>()?;
    let slots = Slots::new(Arc::clone(structure), [(id, points)], [])?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        for (&point, &(_, expected)) in bound.discounts()[0].points().iter().zip(&expected) {
            assert_close!(view.df(point).value(), expected, 1e-12, 1e-14);
        }
        Ok(())
    })?;
    for case in [&tona, &sofr] {
        for (actual, node) in par_rates(case, &collateral, market, &context)?
            .into_iter()
            .zip(&case.nodes)
        {
            assert_close!(actual, node.rate, 1e-10, 1e-12);
        }
    }
    Ok(())
}
`````

`````
crates/rhodium-pricing/tests/par.rs
//! 非線形のパーを AD の入れ子で復元する。

use rhodium_core::{Error, Result};
use rhodium_pricing::par::solve_nonlinear;
use rhodium_testkit::assert_close;
use scioxide_ad::{Dual, Real, Residual};
use scioxide_core::{Interval, SciOxideError};
use scioxide_rootfind::SafeguardedInitialGuess;

struct Par;

impl Residual for Par {
    fn residual<R: Real<f64, Mask = bool>>(
        &self,
        rate: &[R],
        parameters: &[R],
        out: &mut [R],
    ) -> std::result::Result<(), SciOxideError> {
        out[0] = rate[0] * rate[0] - parameters[0];
        Ok(())
    }
}

fn initial(lower: f64, upper: f64) -> Result<SafeguardedInitialGuess<f64>> {
    Interval::new(lower, upper)
        .and_then(|interval| SafeguardedInitialGuess::new(interval, f64::midpoint(lower, upper)))
        .map_err(|error| Error::InvalidInput(error.to_string()))
}

#[test]
fn par_second_derivative() -> Result<()> {
    // p² − θ = 0 の正の根は θ¹ᐟ²。θ = 1 で二階微分は −1/4。
    let theta = Dual::new(Dual::new(1.0, 1.0), Dual::lift(1.0));
    let p = solve_nonlinear(&Par, &[theta], initial(0.1, 3.0)?)?;
    assert_close!(p.value().value(), 1.0, 1e-12, 1e-14);
    assert_close!(p.derivative().value(), 0.5, 1e-10, 1e-12);
    assert_close!(p.derivative().derivative(), -0.25, 1e-8, 1e-12);
    Ok(())
}

#[test]
fn par_third_derivative() -> Result<()> {
    // 同じ正の根の三階微分は 3/8。R への復元で三階の項まで伝わる。
    let theta = Dual::new(
        Dual::new(Dual::new(1.0, 1.0), Dual::lift(1.0)),
        Dual::lift(Dual::lift(1.0)),
    );
    let p = solve_nonlinear(&Par, &[theta], initial(0.1, 3.0)?)?;
    assert_close!(p.derivative().derivative().derivative(), 0.375, 1e-8, 1e-12);
    Ok(())
}

#[test]
fn par_rejects_unbracketed_or_singular_roots() -> Result<()> {
    // 根が区間外にある場合と、ゼロ残差でも傾きがゼロの場合を拒否する。
    assert!(matches!(
        solve_nonlinear(&Par, &[1.0], initial(2.0, 3.0)?),
        Err(Error::Calibration(_))
    ));
    assert!(matches!(
        solve_nonlinear(&Par, &[0.0], initial(0.0, 2.0)?),
        Err(Error::Calibration(_))
    ));
    Ok(())
}

/// quote の値 p について affine でない PV: (1 + p)^n − a·exp(b·θ)。ZC インフレスワップの
/// 固定脚 N·((1 + K)^n − 1) と同じ形で、パーは閉形式 p* = (a·e^{bθ})^{1/n} − 1 を持つ。
struct ZeroCoupon {
    years: f64,
    scale: f64,
    slope: f64,
}

impl Residual for ZeroCoupon {
    fn residual<R: Real<f64, Mask = bool>>(
        &self,
        rate: &[R],
        parameters: &[R],
        out: &mut [R],
    ) -> std::result::Result<(), SciOxideError> {
        out[0] = (rate[0] + R::constant(1.0)).powf(R::constant(self.years))
            - (parameters[0] * R::constant(self.slope)).exp() * R::constant(self.scale);
        Ok(())
    }
}

impl ZeroCoupon {
    /// 閉形式のパーと、θ についての 1 階・2 階の微分。
    fn exact(&self, theta: f64) -> [f64; 3] {
        let growth = (self.scale * (self.slope * theta).exp()).powf(self.years.recip());
        let rate = self.slope / self.years;
        [growth - 1.0, rate * growth, rate * rate * growth]
    }
}

/// 二つの非線形のパーを重み 1 で当てる最小二乗の停留点 θ*(q) を、閉形式のパーのニュートン法で解く。
fn reoptimize(instruments: &[ZeroCoupon; 2], quotes: [f64; 2]) -> f64 {
    let mut theta = 0.0;
    for _ in 0..100 {
        let (gradient, hessian) = instruments.iter().zip(quotes).fold(
            (0.0, 0.0),
            |(gradient, hessian), (instrument, quote)| {
                let [m, first, second] = instrument.exact(theta);
                (
                    gradient + (m - quote) * first,
                    hessian + first * first + (m - quote) * second,
                )
            },
        );
        let step = gradient / hessian;
        theta -= step;
        if step.abs() < 1e-17 {
            break;
        }
    }
    theta
}

#[test]
fn par_nonlinear_lsq_quote_delta() -> Result<()> {
    // 残差の残る最小二乗 min ½Σ(m_i(θ) − q_i)² のモデルの quote m_i が非線形のパー。停留条件
    // G = Σ(m_i − q_i)·m_i′ の θ 微分 H = Σ(m_i′² + (m_i − q_i)·m_i″) はパーの 2 階の微分を含むので、
    // quote デルタ dθ/dq_j = m_j′/H を正しく求めるには m_i″ まで正しい必要がある。
    // 陰関数定理を入れ子の AD の次数まで当てる solve_nonlinear は再最適化の有限差分と一致し、
    // 1 回のニュートン補正で R の値を作る近道（2 階の微分を誤る）は一致しない。
    let instruments = [
        ZeroCoupon {
            years: 2.0,
            scale: 1.02,
            slope: 1.0,
        },
        ZeroCoupon {
            years: 5.0,
            scale: 1.1,
            slope: 3.0,
        },
    ];
    let quotes = [0.03, 0.0];
    let theta = reoptimize(&instruments, quotes);
    let seeded = Dual::new(Dual::new(theta, 1.0), Dual::lift(1.0));
    let mut implicit = Vec::new();
    let mut shortcut = Vec::new();
    for (instrument, quote) in instruments.iter().zip(quotes) {
        let [m, first, second] = instrument.exact(theta);
        // 停留点でも残差は 0 にならない（当てはめの残差が 2 階の項に効く）。
        assert!((m - quote).abs() > 1e-3);
        let par = solve_nonlinear(instrument, &[seeded], initial(-0.5, 1.0)?)?;
        assert_close!(par.value().value(), m, 1e-12, 1e-14);
        assert_close!(par.derivative().value(), first, 1e-10, 1e-12);
        assert_close!(par.derivative().derivative(), second, 1e-8, 1e-12);
        implicit.push([m, first, par.derivative().derivative()]);
        // 近道: p = p₀ − PV(p₀, θ)/∂PV/∂p(p₀, θ)。1 階の微分は正しいが、2 階の微分は n 倍になる。
        let mut pv = [Dual::lift(Dual::lift(0.0))];
        instrument
            .residual(&[Dual::lift(Dual::lift(m))], &[seeded], &mut pv)
            .map_err(|error| Error::Numerical(error.to_string()))?;
        let slope = instrument.years * (1.0 + m).powf(instrument.years - 1.0);
        let approximate = Dual::lift(Dual::lift(m)) - pv[0] * Dual::lift(Dual::lift(slope.recip()));
        assert_close!(approximate.derivative().value(), first, 1e-10, 1e-12);
        shortcut.push([m, first, approximate.derivative().derivative()]);
    }
    let delta = |parts: &[[f64; 3]], j: usize| {
        let hessian: f64 = parts
            .iter()
            .zip(quotes)
            .map(|(&[m, first, second], quote)| first * first + (m - quote) * second)
            .sum();
        parts[j][1] / hessian
    };
    let h = 1e-5;
    for j in 0..2 {
        let bumped = |shift: f64| {
            let mut quotes = quotes;
            quotes[j] += shift;
            reoptimize(&instruments, quotes)
        };
        let finite_difference = (bumped(h) - bumped(-h)) / (2.0 * h);
        assert_close!(delta(&implicit, j), finite_difference, 1e-7, 1e-12);
        let error = (delta(&shortcut, j) - finite_difference).abs() / finite_difference.abs();
        assert!(error > 1e-4, "shortcut quote delta error {error:e}");
    }
    Ok(())
}
`````

`````
crates/rhodium-pricing/tests/xccy_curves.rs
//! 外貨担保曲線の長期を `MtM` XCCY のベーシスで校正し、短期の FX ポイントと一つのグループで解く（R6）。

use std::ptr;

use rhodium_core::{Error, Result, currency::Currency, fx::FxRate, money::Money};
use rhodium_market::quote::{QuoteId, QuoteSet};
use rhodium_pricing::{
    calibration::{CalibratedMarket, GroupOutcome, calibrate},
    context::Context,
    definition::FxImpliedLongEnd,
    prepare::{PrepareInput, evaluate, prepare},
};
use rhodium_products::standard::{FxForward, FxSwap, Position, StandardProduct};
use rhodium_testkit::assert_close;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、XCCY の校正に使わない関数がある"
)]
#[path = "fixtures/ois.rs"]
pub mod ois;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、XCCY の校正に使わない関数がある"
)]
#[path = "fixtures/xccy.rs"]
pub mod xccy;

use ois::midnight;
use xccy::{POINTS, SPOT, XCCY_QUOTES, XccyMarket, jpy_cash};

const TOLERANCE: f64 = 1e-14;

/// 評価日（2025-01-07、評価時点は UTC の 0 時）。校正商品の約定日と同じ。
fn trade_date() -> Result<rhodium_core::date::Date> {
    "2025-01-07".parse()
}

fn solved(market: &CalibratedMarket) -> Result<()> {
    for group in market.groups() {
        let GroupOutcome::Solved(solution) = group.outcome() else {
            return Err(Error::Calibration(format!("{group:?}")));
        };
        assert!(solution.residual_norm() <= TOLERANCE);
    }
    Ok(())
}

/// JPY 建ての PV。USD の受渡を評価時点の X(0) で換算する。
fn jpy_pv(product: &StandardProduct, market: &CalibratedMarket, ctx: &Context<'_>) -> Result<f64> {
    let plan = prepare(
        PrepareInput::Standard {
            product,
            trade_date: trade_date()?,
        },
        &jpy_cash()?,
        &[],
        market.market().structure(),
        ctx,
    )?;
    let pv = evaluate(&plan, market.market())?;
    let x0 = market
        .market()
        .structure()
        .fx_today_plan(XccyMarket::pair()?)?
        .rate(market.market())?
        .rate();
    Ok(pv.get(Currency::JPY) + pv.get(Currency::USD) * x0)
}

/// `MtM` XCCY のテンプレートのパー（JPY の TONA 脚のベーシス）。PV は K について affine。
fn xccy_par(tenor: &str, market: &CalibratedMarket, ctx: &Context<'_>) -> Result<f64> {
    let template = XccyMarket::template(tenor)?;
    let pv = |spread| {
        let product = template.build_for_trade(spread, trade_date()?, ctx.refdata)?;
        jpy_pv(&product, market, ctx)
    };
    let zero = pv(0.0)?;
    Ok(-zero / (pv(1.0)? - zero))
}

#[test]
fn fx_implied_long_end_reproduces_fx_points_and_xccy_quotes() -> Result<()> {
    // FX ポイント（ON、TN、3M、6M）と XCCY（1Y、18M）が、外貨担保曲線の basis ブロックの
    // 6 つのピラーを一つのグループで決める。TONA と SOFR の OIS と USDJPY のスポットが上流。
    let fixture = XccyMarket::new()?;
    let inputs = fixture.inputs("")?;
    let ctx = inputs.context();
    let market = calibrate(&fixture.definition()?, &fixture.quotes()?, &ctx)?;
    solved(&market)?;
    let block = rhodium_market::structure::BlockId::new(XccyMarket::implied()?, "basis")?;
    let group = market
        .groups()
        .iter()
        .find(|group| group.blocks() == [block.clone()])
        .ok_or_else(|| Error::Internal("missing FX-implied group".into()))?;
    let rows: Vec<&str> = group.row_quotes().iter().map(QuoteId::as_str).collect();
    for quote in POINTS
        .iter()
        .map(|(quote, _, _)| *quote)
        .chain(XCCY_QUOTES.iter().map(|(quote, _, _, _)| *quote))
    {
        assert!(rows.contains(&quote), "{quote} not in {rows:?}");
    }
    assert_eq!(rows.len(), 6);
    // XCCY を同じ評価で建て直したパーは quote（残差の許容誤差 1e-14 の範囲）。
    for (_, tenor, quote, _) in XCCY_QUOTES {
        assert_close!(xccy_par(tenor, &market, &ctx)?, quote, 0.0, 1e-12);
    }
    Ok(())
}

#[test]
fn fx_implied_calibration_instruments_are_at_par_as_trades() -> Result<()> {
    // 校正商品を取引として評価すると PV はゼロ（校正の残差を PV に伝えた幅の内）。
    // FX フォワードは 3M のアウトライト X(S) + pts·pip、FX スワップはスポットの近い受渡と 6M の遠い受渡、
    // `MtM` XCCY は quote のベーシスで約定する。元本は 100 万 USD（XCCY は JPY 1 円あたり）。
    let fixture = XccyMarket::new()?;
    let inputs = fixture.inputs("")?;
    let ctx = inputs.context();
    let market = calibrate(&fixture.definition()?, &fixture.quotes()?, &ctx)?;
    let pair = XccyMarket::pair()?;
    let notional = Money::new(Currency::USD, 1_000_000.0);
    let pip = 0.01;
    let [_, _, three, six] = POINTS;
    let outright = |points: f64| SPOT + points * pip;
    let settle = |date: Option<&str>| {
        XccyMarket::settlement(date.ok_or_else(|| Error::Internal("dated point".into()))?)
    };
    let forward = StandardProduct::FxForward(FxForward::new(
        notional,
        FxRate::forward(pair, outright(three.2), settle(three.1)?)?,
        Position::Long,
    )?);
    let swap = StandardProduct::FxSwap(FxSwap::new(
        notional,
        FxRate::forward(pair, SPOT, XccyMarket::settlement("2025-01-09")?)?,
        FxRate::forward(pair, outright(six.2), settle(six.1)?)?,
        Position::Long,
    )?);
    // 1 pip（0.01 円）の誤差は 100 万 USD で 1 万円なので、校正の許容誤差 1e-14 の伝わる幅は 1e-6 円の程度。
    for product in [&forward, &swap] {
        assert_close!(jpy_pv(product, &market, &ctx)?, 0.0, 0.0, 1e-4);
    }
    for (_, tenor, quote, _) in XCCY_QUOTES {
        let product =
            XccyMarket::template(tenor)?.build_for_trade(quote, trade_date()?, ctx.refdata)?;
        // 元本は JPY 1 円。ベーシスの PV の傾き（年金現価、約 1 円/1.0）に許容誤差を掛けた幅。
        assert_close!(jpy_pv(&product, &market, &ctx)?, 0.0, 0.0, 1e-12);
    }
    Ok(())
}

fn bumped(quotes: &QuoteSet, id: &str, bump: f64) -> Result<QuoteSet> {
    let id: QuoteId = id.parse()?;
    QuoteSet::new(
        quotes
            .iter()
            .map(|(key, value)| (key.clone(), if *key == id { value + bump } else { value })),
    )
}

#[test]
fn fx_implied_long_end_update_matches_full_recalibration() -> Result<()> {
    // XCCY の quote の変更は外貨担保曲線のグループだけを解き直し、TONA、SOFR と FX のマーク
    // （スポットと短期のポイントの直接マーク）の解を共有する。TONA の quote の変更は TONA と、その下流の
    // 外貨担保曲線を解き直し、SOFR と FX のマークを共有する。
    // どちらも全体を解き直した校正と θ が一致する。
    let fixture = XccyMarket::new()?;
    let inputs = fixture.inputs("")?;
    let ctx = inputs.context();
    let definition = fixture.definition()?;
    let quotes = fixture.quotes()?;
    let base = calibrate(&definition, &quotes, &ctx)?;
    let tona_quote = fixture.tona.quote(&fixture.tona.nodes[2])?;
    for (id, bump, shared) in [
        (
            "USDJPY.XCCY.18M",
            0.0002,
            ["JPY-TONA", "USD-SOFR", "FX/USDJPY"].as_slice(),
        ),
        (
            tona_quote.as_str(),
            0.0001,
            ["USD-SOFR", "FX/USDJPY"].as_slice(),
        ),
    ] {
        let quotes = bumped(&quotes, id, bump)?;
        let updated = base.update(&quotes)?;
        let full = calibrate(&definition, &quotes, &ctx)?;
        solved(&updated)?;
        for (lhs, rhs) in updated.market().theta().iter().zip(full.market().theta()) {
            assert_close!(*lhs, *rhs, 1e-12, 1e-14);
        }
        for (old, new) in base.groups().iter().zip(updated.groups()) {
            let curve = old.blocks()[0].curve().as_str();
            assert_eq!(
                ptr::eq(old.outcome(), new.outcome()),
                shared.contains(&curve),
                "{id}: {curve}"
            );
        }
    }
    Ok(())
}

#[test]
fn fx_implied_long_end_rejects_pillars_before_fx_points() -> Result<()> {
    // 長期のピラーは昇順で、最後の FX ポイントの受渡より後でなければならない。
    let fixture = XccyMarket::new()?;
    let inputs = fixture.inputs("")?;
    let ctx = inputs.context();
    let mut long_end = XccyMarket::long_end()?;
    long_end.reverse();
    let uncollateralized = rhodium_pricing::resolve::CollateralKey::Uncollateralized;
    let base = || {
        ois::definition(&[
            (&fixture.tona, &uncollateralized),
            (&fixture.sofr, &uncollateralized),
        ])
    };
    assert!(matches!(
        base()?.with_fx([fixture.fx(long_end)?]),
        Err(Error::InvalidData(message)) if message.contains("strictly increasing")
    ));
    let early = vec![FxImpliedLongEnd {
        pillar: midnight("2025-06-09".parse()?)?,
        ..XccyMarket::long_end()?.remove(0)
    }];
    let definition = base()?.with_fx([fixture.fx(early)?])?;
    assert!(matches!(
        calibrate(&definition, &fixture.quotes()?, &ctx),
        Err(Error::InvalidData(message)) if message.contains("follow the last FX point")
    ));
    Ok(())
}
`````

`````
crates/rhodium-products/Cargo.toml
[package]
name = "rhodium-products"
description = "コア IR、DSL、標準商品と取引の定義"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
publish.workspace = true

[lints]
workspace = true

[dependencies]
rhodium-core.workspace = true
scioxide-core.workspace = true
serde.workspace = true

[dev-dependencies]
rhodium-testkit.workspace = true
`````

`````
crates/rhodium-products/src/agreement.rs
//! 取引の担保条件。

use rhodium_core::{
    Error, Result,
    calendar::CalendarSpec,
    currency::Currency,
    day_count::DayCount,
    id::{IndexId, MarginSetId, SecurityId},
    quantity::{Rate, Simple, Spread},
    refdata::RefData,
};

/// 取引に適用するマージンセットの VM 条件。
#[derive(Debug, Clone, PartialEq)]
pub struct MarginSet {
    /// マージンセットの ID。
    pub id: MarginSetId,
    /// 変動証拠金の CSA。
    pub vm: Option<CsaTerms>,
}

/// CSA の適格担保条件。
#[derive(Debug, Clone, PartialEq)]
pub struct CsaTerms {
    eligible: Vec<EligibleCollateral>,
}

impl CsaTerms {
    /// 適格担保とその参照先を照合し、CSA を構築する。
    ///
    /// # Errors
    /// 適格担保が空、または参照先が存在しなければ `InvalidData`。
    pub fn new(eligible: Vec<EligibleCollateral>, refdata: &RefData) -> Result<Self> {
        if eligible.is_empty() {
            return Err(Error::InvalidData(
                "a VM CSA requires eligible collateral".into(),
            ));
        }
        for collateral in &eligible {
            match collateral.asset() {
                CollateralAsset::Cash(terms) => terms.validate_references(refdata)?,
                CollateralAsset::Security(id) => {
                    refdata.security(id)?;
                }
            }
        }
        Ok(Self { eligible })
    }

    /// CSA が指定する適格担保。
    #[must_use]
    pub fn eligible(&self) -> &[EligibleCollateral] {
        &self.eligible
    }
}

/// 自社と取引先の条件。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PerParty<T> {
    /// 自社の条件。
    pub own: T,
    /// 取引先の条件。
    pub cpty: T,
}

/// 当事者が担保を差し入れられるか。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Eligibility {
    /// 差し入れられる。
    Eligible,
    /// 差し入れられない。
    Ineligible,
}

/// CSA が指定する担保資産。
#[derive(Debug, Clone, PartialEq)]
pub enum CollateralAsset {
    /// 付利条件付きの現金。
    Cash(CashCollateralTerms),
    /// 証券。
    Security(SecurityId),
}

/// 担保資産、ヘアカット、差し入れられる当事者の組。
#[derive(Debug, Clone, PartialEq)]
pub struct EligibleCollateral {
    asset: CollateralAsset,
    haircut: f64,
    posters: PerParty<Eligibility>,
}

impl EligibleCollateral {
    /// ヘアカットを検証し、適格担保の条件を構築する。
    ///
    /// # Errors
    /// ヘアカットが 0 以上 1 未満でなければ `InvalidData`。
    pub fn new(
        asset: CollateralAsset,
        haircut: f64,
        posters: PerParty<Eligibility>,
    ) -> Result<Self> {
        if !(0.0..1.0).contains(&haircut) {
            return Err(Error::InvalidData(
                "collateral haircut must be in [0, 1)".into(),
            ));
        }
        Ok(Self {
            asset,
            haircut,
            posters,
        })
    }

    /// 差し入れる資産。
    #[must_use]
    pub const fn asset(&self) -> &CollateralAsset {
        &self.asset
    }

    /// 担保の評価に適用するヘアカット。
    #[must_use]
    pub const fn haircut(&self) -> f64 {
        self.haircut
    }

    /// 当事者ごとの差し入れの適格性。
    #[must_use]
    pub const fn posters(&self) -> PerParty<Eligibility> {
        self.posters
    }
}

/// 価格設定の割引の照合に使う現金担保条件。
#[derive(Debug, Clone, PartialEq)]
pub struct CashCollateralTerms {
    /// 担保の通貨。
    pub ccy: Currency,
    /// 現金担保に適用する付利条件。
    pub remuneration: Remuneration,
}

impl CashCollateralTerms {
    /// 通貨、付利指標、フロアに必要なカレンダーを参照データに照合する。
    ///
    /// # Errors
    /// 参照先が存在しなければ `InvalidData`。
    pub fn validate_references(&self, refdata: &RefData) -> Result<()> {
        refdata.currency(self.ccy)?;
        refdata.index(&self.remuneration.index)?;
        if let Some(DayCount::Bus252 { calendar }) =
            self.remuneration.floor.as_ref().map(Rate::basis)
        {
            refdata.calendar(&CalendarSpec::Base(calendar.clone()))?;
        }
        Ok(())
    }
}

/// 現金担保の付利。指標にスプレッドを加え、指定があればフロアを適用する。
#[derive(Debug, Clone, PartialEq)]
pub struct Remuneration {
    /// 現金担保に付利する指標。
    pub index: IndexId,
    /// 付利指標に加えるスプレッド。
    pub spread: Spread,
    /// 付利率の下限。
    pub floor: Option<Rate<Simple>>,
}
`````

`````
crates/rhodium-products/src/ir.rs
//! 契約の組み立て、因果性の検査、評価時点での簡約。

mod amendment;
pub mod analysis;
mod builder;
mod fallback;
mod nodes;
mod reduce;
#[cfg(test)]
mod tests;
mod validate;

use std::{collections::HashMap, marker::PhantomData, sync::Arc};

use rhodium_core::{
    Result, Warning,
    brand::{Brand, with_brand},
    date::DateInterval,
    fixing::{ObservationPeriod, ObservationPoint},
    id::{CreditKey, IndexId},
    refdata::RefData,
    time::Instant,
};

pub(crate) use amendment::Exercises;
pub use builder::Builder;
pub use reduce::reduce;
pub(crate) use validate::ObservationCache;

use nodes::{ContractNode, ExprNode, Ix, MAX_NODES, RfrRows, node_limit};

/// 数の式の種類。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Num {}

/// 真偽値の式の種類。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bool {}

/// 同じ builder に属する式への参照。
#[derive(Debug, Clone, Copy)]
pub struct ExprIx<'b, T = Num> {
    ix: Ix,
    _brand: Brand<'b>,
    sort: PhantomData<T>,
}

/// 同じ builder に属する契約への参照。
#[derive(Debug, Clone, Copy)]
pub struct ContractIx<'b> {
    ix: Ix,
    _brand: Brand<'b>,
}

/// 選択を行う当事者。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Party {
    /// ブックの側。
    Us,
    /// 取引先。
    Counterparty,
}

impl Party {
    const fn other(self) -> Self {
        match self {
            Self::Us => Self::Counterparty,
            Self::Counterparty => Self::Us,
        }
    }
}

/// builder で組み立てた、因果性の検査前の契約。
#[derive(Debug, Clone)]
pub struct Contract {
    contracts: Vec<ContractNode>,
    exprs: Vec<ExprNode>,
    root: Ix,
}

impl Contract {
    /// 新しい builder を開き、その中の契約を根にする。
    ///
    /// # Errors
    ///
    /// クロージャのエラーを返す。ノード数が上限を超える場合は `InvalidInput`。
    pub fn build<F>(f: F) -> Result<Self>
    where
        F: for<'b> FnOnce(&mut Builder<'b>) -> Result<ContractIx<'b>>,
    {
        with_brand(|brand| {
            let mut builder = Builder::new(brand);
            let root = f(&mut builder)?;
            builder.finish(root)
        })
    }

    fn reachable_exprs(&self, contracts: &[bool]) -> Vec<bool> {
        let mut reached = vec![false; self.exprs.len()];
        for (node, needed) in self.contracts.iter().zip(contracts) {
            if *needed {
                match node {
                    ContractNode::Pay { amount, .. } => reached[*amount as usize] = true,
                    ContractNode::Cond { cond, .. } => reached[*cond as usize] = true,
                    ContractNode::All(_) | ContractNode::Choose { .. } => {}
                }
            }
        }
        for (pos, node) in self.exprs.iter().enumerate().rev() {
            if reached[pos] {
                for &ix in node.children() {
                    reached[ix as usize] = true;
                }
            }
        }
        reached
    }

    // 簡約の途中で使った条件や定数のうち、結果の根から参照されるものだけを残す。
    fn compact(self) -> Result<Self> {
        let contracts_needed = reachable_contracts(&self.contracts, self.root);
        let exprs_needed = self.reachable_exprs(&contracts_needed);
        let mut expr_ids = vec![0; exprs_needed.len()];
        let mut exprs = Vec::new();
        for (pos, (mut node, needed)) in self.exprs.into_iter().zip(exprs_needed).enumerate() {
            if needed {
                node.try_map_children(|ix| Ok(expr_ids[ix as usize]))?;
                expr_ids[pos] = Ix::try_from(exprs.len()).map_err(|_| node_limit())?;
                exprs.push(node);
            }
        }
        let mut contract_ids = vec![0; contracts_needed.len()];
        let mut contracts = Vec::new();
        for (pos, (mut node, needed)) in
            self.contracts.into_iter().zip(contracts_needed).enumerate()
        {
            if needed {
                match &mut node {
                    ContractNode::Pay { amount, .. } => *amount = expr_ids[*amount as usize],
                    ContractNode::Cond { cond, yes, no, .. } => {
                        *cond = expr_ids[*cond as usize];
                        *yes = contract_ids[*yes as usize];
                        *no = contract_ids[*no as usize];
                    }
                    ContractNode::All(children)
                    | ContractNode::Choose {
                        options: children, ..
                    } => {
                        for ix in children {
                            *ix = contract_ids[*ix as usize];
                        }
                    }
                }
                contract_ids[pos] = Ix::try_from(contracts.len()).map_err(|_| node_limit())?;
                contracts.push(node);
            }
        }
        if contracts.len() + exprs.len() > MAX_NODES {
            return Err(node_limit());
        }
        Ok(Self {
            contracts,
            exprs,
            root: contract_ids[self.root as usize],
        })
    }
}

#[derive(Debug)]
struct ResolvedFixing {
    publication: Instant,
    period: ObservationPeriod,
    warnings: Vec<Warning>,
    ceased: bool,
    non_representative: bool,
    kind: FixingKind,
}

// 一つの RfrRows の公表時刻と観測期間の解決。同じ列を読む取引と簡約後の商品で共有する。
#[derive(Debug)]
struct ResolvedCompounding {
    // RfrRows::entries と同じ順の各行の解決。CompoundingRows の位置でそのまま引く。
    observations: Box<[ResolvedFixing]>,
    // 列の最後の公表時刻。契約の因果性の検査に使う。
    latest: Option<Instant>,
    // 各行の警告を重複なく集めたもの。
    warnings: Vec<Warning>,
}

#[derive(Debug, Clone, Copy)]
enum FixingKind {
    TermRate,
    Fx(rhodium_core::fx::CcyPair),
    RateFuture,
    Other,
}

#[derive(Debug)]
struct ResolvedCredit {
    key: CreditKey,
    publication: Instant,
    warnings: Vec<Warning>,
}

type ByIndex<K, V> = HashMap<IndexId, HashMap<K, V>>;
type ByChoice<V> = HashMap<String, HashMap<Instant, V>>;

#[derive(Debug, Default)]
struct References {
    fixings: ByIndex<ObservationPoint, ResolvedFixing>,
    compounded: ByIndex<Arc<RfrRows>, Arc<ResolvedCompounding>>,
    credit: ByIndex<DateInterval, ResolvedCredit>,
    choices: ByChoice<usize>,
}

/// 因果性と参照データとの整合を検査した契約。
#[derive(Debug, Clone)]
pub struct Product {
    contract: Contract,
    references: Arc<References>,
    forecasts: ByIndex<ObservationPoint, Instant>,
    construction_warnings: Arc<[Warning]>,
    warnings: Vec<Warning>,
    trade_exercises: Arc<[crate::trade::ExerciseRecord]>,
}

impl Product {
    /// 契約を検査し、観測期間と公表時刻を一度だけ解決する。
    ///
    /// # Errors
    ///
    /// 因果性、選択または式の構造に違反すれば `InvalidInput`。
    /// 未登録または種類の合わない参照データがあれば `InvalidData`。
    pub fn new(contract: Contract, refdata: &RefData) -> Result<Self> {
        Self::with_construction_warnings(contract, refdata, Vec::new())
    }

    pub(crate) fn with_construction_warnings(
        contract: Contract,
        refdata: &RefData,
        construction_warnings: Vec<Warning>,
    ) -> Result<Self> {
        Self::validated(contract, refdata, construction_warnings, None)
    }

    fn validated(
        mut contract: Contract,
        refdata: &RefData,
        mut construction_warnings: Vec<Warning>,
        cache: Option<&ObservationCache>,
    ) -> Result<Self> {
        let (references, mut warnings) = validate::validate(&mut contract, refdata, cache)?;
        construction_warnings.sort_unstable();
        construction_warnings.dedup();
        if !construction_warnings.is_empty() {
            warnings.extend_from_slice(&construction_warnings);
            warnings.sort_unstable();
            warnings.dedup();
        }
        Ok(Self {
            contract,
            references: Arc::new(references),
            forecasts: HashMap::new(),
            construction_warnings: construction_warnings.into(),
            warnings,
            trade_exercises: Arc::default(),
        })
    }

    pub(crate) fn with_fallback(
        contract: Contract,
        refdata: &RefData,
        warnings: Vec<Warning>,
        terms: &crate::trade::FallbackTerms,
        trade_date: rhodium_core::date::Date,
        cache: Option<&ObservationCache>,
    ) -> Result<Self> {
        let contract = fallback::lower(contract, refdata, terms, trade_date)?;
        Self::validated(contract, refdata, warnings, cache)
    }

    pub(crate) fn with_trade_exercises(mut self, records: &[crate::trade::ExerciseRecord]) -> Self {
        self.trade_exercises = records
            .iter()
            .filter(|record| {
                self.references
                    .choices
                    .get(&record.label)
                    .is_some_and(|times| times.contains_key(&record.at))
            })
            .cloned()
            .collect();
        self
    }
