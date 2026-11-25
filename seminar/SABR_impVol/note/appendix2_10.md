    );
    Ok((RefData::new(parts)?, term, overnight))
}

fn fixing_trade(index: &IndexId, date: Date) -> Result<Trade> {
    let payment = "2027-01-12T08:00:00Z".parse()?;
    trade(Contract::build(|b| {
        let rate = b.fixing(index, ObservationPoint::Date(date));
        Ok(b.pay(CASH, rate, payment))
    })?)
}

#[test]
fn fallback_switches_to_replacement_fixings_and_their_publications() -> Result<()> {
    // 金曜の観測は 3 日、月曜は 1 日分。月曜だけ予想に残る場合と、両方確定した場合を比べる。
    let since = "2027-01-08".parse()?;
    let (data, term, overnight) = fallback_data(IndexStatus::Ceased { since })?;
    let date = ObservationPoint::Date(since);
    let monday = ObservationPoint::Date("2027-01-11".parse()?);
    let previous = "2027-01-07".parse()?;
    let rows = vec![
        (overnight.clone(), date, 0.03),
        (term.clone(), ObservationPoint::Date(previous), 0.02),
    ];
    let partial = FixingStore::new(
        &data,
        FixingStoreParts {
            published: rows.clone(),
            ..Default::default()
        },
    )?;
    let full = FixingStore::new(
        &data,
        FixingStoreParts {
            published: rows
                .into_iter()
                .chain([(overnight.clone(), monday, 0.06)])
                .collect(),
            ..Default::default()
        },
    )?;
    let mut trade = fixing_trade(&term, since)?;
    trade.fallback.protocols.insert("p".into());
    let product = trade.to_product(&data)?;
    let reduced = reduce(
        &product,
        "2027-01-11T01:00:00Z".parse()?,
        &partial,
        &[],
        SettlementBoundary::Exclude,
    )?;
    assert_eq!(
        observations(&reduced),
        [Observation::Fixing {
            index: &overnight,
            point: monday
        }]
    );
    assert!(reduced.warnings().contains(&Warning::FallbackApplied {
        index: term.clone(),
        replacement: overnight
    }));

    let as_of = "2027-01-12T01:00:00Z".parse()?;
    assert!(
        matches!(reduce(&product, as_of, &partial, &[], SettlementBoundary::Exclude), Err(Error::MissingData(message)) if message.contains("JPY-TONA"))
    );
    let resolved = reduce(&product, as_of, &full, &[], SettlementBoundary::Exclude)?;
    let first = 0.03 * 3.0 / 365.0;
    let second = 0.06 / 365.0;
    let expected = (first + second + first * second) * 365.0 / 4.0 + 0.001;
    assert_close!(const_amount(&resolved)?, expected, 1e-12, 1e-6);

    // 取引固有のスプレッドを優先する。適用する定義が無い場合は停止した指標を予想しない。
    let mut replacement = data
        .fallback(&term)
        .ok_or_else(|| Error::Internal("test fallback".into()))?
        .replacement
        .clone();
    replacement.spread = "20bp".parse()?;
    trade.fallback.by_index.insert(
        term.clone(),
        FallbackSelection::Custom {
            trigger: FallbackTrigger::Cessation,
            replacement,
        },
    );
    let adjusted = reduce(
        &trade.to_product(&data)?,
        as_of,
        &full,
        &[],
        SettlementBoundary::Exclude,
    )?;
    assert_close!(const_amount(&adjusted)?, expected + 0.001, 1e-12, 1e-6);
    trade.fallback.by_index.clear();
    trade.fallback.protocols.clear();
    assert!(
        matches!(reduce(&trade.to_product(&data)?, as_of, &full, &[], SettlementBoundary::Exclude), Err(Error::MissingData(message)) if message.contains("no fallback for ceased index"))
    );

    trade.product = fixing_trade(&term, previous)?.product;
    let before = reduce(
        &trade.to_product(&data)?,
        as_of,
        &full,
        &[],
        SettlementBoundary::Exclude,
    )?;
    assert_close!(const_amount(&before)?, 0.02, 1e-12, 1e-6);
    assert_eq!(before.warnings(), []);
    Ok(())
}

#[test]
fn fallback_terms_use_explicit_precedence_and_warn_on_the_legacy_rate() -> Result<()> {
    // 適用可能な取引でも明示の無効化を優先する。既定を指定すれば版の条件を上書きできる。
    let since = "2027-01-08".parse()?;
    let (data, term, _) = fallback_data(IndexStatus::NonRepresentative { since })?;
    let as_of = "2027-01-11T01:00:00Z".parse()?;
    let mut trade = fixing_trade(&term, since)?;
    trade.fallback.definitions = Some("2021".into());
    let legacy = FixingStore::new(
        &data,
        FixingStoreParts {
            published: vec![(term.clone(), ObservationPoint::Date(since), 0.017)],
            ..Default::default()
        },
    )?;
    let check_legacy = |trade: &Trade| -> Result<()> {
        let legacy = reduce(
            &trade.to_product(&data)?,
            as_of,
            &legacy,
            &[],
            SettlementBoundary::Exclude,
        )?;
        assert_close!(const_amount(&legacy)?, 0.017, 1e-12, 1e-6);
        assert_eq!(
            legacy.warnings(),
            [Warning::NonRepresentativeIndex {
                index: term.clone()
            }]
        );
        Ok(())
    };
    check_legacy(&trade)?;
    trade.trade_date = "2027-01-05".parse()?;
    let requires_replacement = |trade: &Trade| -> Result<()> {
        assert!(
            matches!(reduce(&trade.to_product(&data)?, as_of, &legacy, &[], SettlementBoundary::Exclude), Err(Error::MissingData(message)) if message.contains("JPY-TONA"))
        );
        Ok(())
    };
    requires_replacement(&trade)?;
    trade.fallback.bilateral.insert(term.clone());
    trade.fallback.definitions = None;
    requires_replacement(&trade)?;
    trade
        .fallback
        .by_index
        .insert(term.clone(), FallbackSelection::Disabled);
    check_legacy(&trade)?;
    trade.fallback.bilateral.clear();
    trade
        .fallback
        .by_index
        .insert(term, FallbackSelection::IndexDefault);
    requires_replacement(&trade)
}

#[test]
fn trade_fallback_terms_reject_invalid_references() -> Result<()> {
    // 取引のフォールバックの指定は IR を作る前に参照データと照らす。ターム金利でない指標への個別の指定、
    // 未知の指標との個別の合意、参照データと合わない固有の代替は、どれも InvalidData。
    let (data, term, overnight) = fallback_data(IndexStatus::Active)?;
    let mut trade = fixing_trade(&term, "2027-01-08".parse()?)?;
    let custom = |index: &str| -> Result<FallbackSelection> {
        let replacement = data
            .fallback(&term)
            .ok_or_else(|| Error::Internal("missing default fallback".into()))?
            .replacement
            .clone();
        Ok(FallbackSelection::Custom {
            trigger: FallbackTrigger::Cessation,
            replacement: FallbackReplacement {
                index: index.parse()?,
                ..replacement
            },
        })
    };
    let reject = |trade: &Trade, expected: &str| {
        assert!(
            matches!(trade.to_product(&data), Err(Error::InvalidData(message)) if message.contains(expected)),
            "{expected}"
        );
    };
    // 個別の指定はターム金利だけに書ける（翌日物の TONA は対象外）。
    trade
        .fallback
        .by_index
        .insert(overnight.clone(), FallbackSelection::Disabled);
    reject(&trade, "is not a term rate");
    trade.fallback.by_index.clear();
    // 個別に合意した指標は参照データに無ければならない。
    trade.fallback.bilateral.insert("JPY-UNKNOWN".parse()?);
    reject(&trade, "JPY-UNKNOWN");
    trade.fallback.bilateral.clear();
    // 取引固有の代替の指標も参照データの翌日物でなければならない。
    trade
        .fallback
        .by_index
        .insert(term.clone(), custom("JPY-UNKNOWN")?);
    reject(&trade, "JPY-UNKNOWN");
    trade
        .fallback
        .by_index
        .insert(term.clone(), custom(term.as_str())?);
    assert!(matches!(
        trade.to_product(&data),
        Err(Error::InvalidData(_))
    ));
    // 正しい固有の代替は受け付ける。
    let valid = custom(overnight.as_str())?;
    trade.fallback.by_index.insert(term.clone(), valid);
    trade.to_product(&data)?;
    Ok(())
}
`````

`````
crates/rhodium-products/src/ir/tests/standard.rs
use std::num::NonZeroU8;

use rhodium_core::{
    Error, Result, Warning,
    calendar::{
        BaseCalendar, BusinessDayConvention, CalendarSpec, EndOfMonthRule, Extrapolation,
        WeekendRule,
    },
    currency::{Asset, Currency},
    date::{Date, DateInterval, Weekday},
    day_count::DayCount,
    fixing::{FixingStore, FixingStoreParts, ObservationPoint},
    id::{CalendarId, IndexId},
    money::Money,
    observation::{
        Aggregation, CompoundingForm, ObservationDays, ObservationMethod, ObservationRounding,
        ObservationRule, SpreadTreatment,
    },
    quantity::{Rate, Simple, Spread},
    refdata::{
        Convention, CurrencyDefinition, DailyPublication, Deliverability, DepositConvention,
        FixedLegConvention, IndexDefinition, IndexKind, IndexStatus, OisConvention,
        OvernightLegConvention, PublicationRule, RefData, RefDataParts,
    },
    rounding::{Rounding, RoundingMode},
    schedule::{
        EndRule, ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PaymentGrouping, PaymentRule,
        PeriodAnchor, PeriodCombination, ScheduleSpec,
    },
    tenor::Frequency,
    time::{Instant, LocalEventTime, UtcOffset, ZoneRules},
};
use rhodium_testkit::assert_close;

use crate::standard::{
    Coupon, Deposit, DepositDates, InstrumentTemplate, Leg, NotionalExchange, PayReceive,
    PaymentRounding, Position, ProductKind, StandardProduct, Swap,
};

use super::{ContractNode, ExprNode, Product, SettlementBoundary, reduce};

fn refdata() -> Result<RefData> {
    refdata_rounding(RoundingMode::TowardZero)
}

/// 受渡額の丸め方だけを指定した参照データ。増分は各通貨の最小単位。
fn refdata_rounding(mode: RoundingMode) -> Result<RefData> {
    let mut parts = RefDataParts::default();
    for (ccy, id, zone, offset, basis, name, holidays, minor_units) in [
        (
            Currency::JPY,
            "JPTO",
            "Asia/Tokyo",
            9 * 3600,
            DayCount::Act365F,
            "JPY-TONA",
            Vec::new(),
            0,
        ),
        (
            Currency::USD,
            "USNY",
            "America/New_York",
            -4 * 3600,
            DayCount::Act360,
            "USD-SOFR",
            vec!["2027-07-05".parse()?],
            2,
        ),
    ] {
        let id: CalendarId = id.parse()?;
        let clock = LocalEventTime {
            time: "15:00:00".parse()?,
            zone: zone.parse()?,
        };
        // 先頭のスタブは満期から 1 刻み遡った仮想の通常期間を参照枠にするため、
        // 2027 年初めに始まる OIS でも 2026 年の営業日判定が要る。
        parts.calendars.insert(
            id.clone(),
            BaseCalendar::new(
                "2026-01-01".parse()?..="2028-12-31".parse()?,
                WeekendRule::new([Weekday::Sat, Weekday::Sun]),
                holidays,
                Extrapolation::Strict,
            )?,
        );
        parts.zones.insert(
            clock.zone.clone(),
            ZoneRules::new(
                "2027-01-01T00:00:00Z".parse()?..="2029-01-01T00:00:00Z".parse()?,
                UtcOffset::from_seconds(offset)?,
                Vec::new(),
            )?,
        );
        let payment_rounding = Rounding::new(0.1_f64.powi(i32::from(minor_units)), mode)?;
        parts.currencies.insert(
            ccy,
            CurrencyDefinition {
                code: ccy,
                minor_units,
                payment_rounding,
                deliverability: Deliverability::Deliverable,
                settlement_time: clock.clone(),
                settlement_calendar: id.clone(),
            },
        );
        let calendar = CalendarSpec::Base(id);
        let index: IndexId = name.parse()?;
        parts.indices.insert(
            index.clone(),
            IndexDefinition {
                id: index.clone(),
                kind: IndexKind::OvernightRate {
                    ccy,
                    day_count: basis.clone(),
                },
                calendar: calendar.clone(),
                publication: PublicationRule::Daily(DailyPublication {
                    business_days: 1,
                    calendar: calendar.clone(),
                    time: "10:00:00".parse()?,
                    zone: clock.zone.clone(),
                }),
                observation_time: Some(clock),
                status: IndexStatus::Active,
            },
        );
        insert_conventions(&mut parts, ccy, index, basis, calendar)?;
    }
    RefData::new(parts)
}

fn insert_conventions(
    parts: &mut RefDataParts,
    ccy: Currency,
    index: IndexId,
    basis: DayCount,
    calendar: CalendarSpec,
) -> Result<()> {
    parts.conventions.insert(
        format!("{ccy}-DEP").parse()?,
        Convention::Deposit(DepositConvention {
            spot_lag: "+2BD".parse()?,
            day_count: basis.clone(),
            bdc: BusinessDayConvention::ModifiedFollowing,
            end_of_month: EndOfMonthRule::LastBusinessDay,
            calendar: calendar.clone(),
        }),
    );
    parts.conventions.insert(
        format!("{index}-OIS").parse()?,
        Convention::Ois(Box::new(OisConvention {
            spot_lag: "+2BD".parse()?,
            fixed: FixedLegConvention {
                frequency: Frequency::ANNUAL,
                day_count: basis,
                calendar: calendar.clone(),
                bdc: BusinessDayConvention::ModifiedFollowing,
            },
            float: OvernightLegConvention {
                index,
                observation: rule(
                    calendar.clone(),
                    ObservationMethod::Plain,
                    SpreadTreatment::Exclusive,
                    ObservationRounding::default(),
                )?,
            },
            payment: PaymentRule {
                anchor: PeriodAnchor::PeriodEnd,
                offset: "+2BD".parse()?,
                calendar,
                bdc: BusinessDayConvention::Following,
                grouping: PaymentGrouping::EachPeriod,
                settlement: None,
            },
        })),
    );
    Ok(())
}

fn rule(
    calendar: CalendarSpec,
    method: ObservationMethod,
    spread: SpreadTreatment,
    rounding: ObservationRounding,
) -> Result<ObservationRule> {
    ObservationRule::new(
        ObservationDays::BusinessDays(calendar),
        method,
        Aggregation::Compound {
            form: CompoundingForm::SimpleFactor,
            spread,
        },
        rounding,
    )
}

fn fixings(refdata: &RefData, index: &str, rows: &[(&str, f64)]) -> Result<FixingStore> {
    let index: IndexId = index.parse()?;
    let published = rows
        .iter()
        .map(|&(date, value)| Ok((index.clone(), ObservationPoint::Date(date.parse()?), value)))
        .collect::<Result<_>>()?;
    FixingStore::new(
        refdata,
        FixingStoreParts {
            published,
            ..FixingStoreParts::default()
        },
    )
}

/// `reduce` 後の受渡を (時刻, 確定した額) の列にする。未確定の額は `None`。
/// `give` は構築時に額の符号反転になるので、契約は `All` と `Pay` だけからなる。
fn reduced_payments(
    product: &Product,
    as_of: &str,
    fixings: &FixingStore,
    ccy: Currency,
) -> Result<Vec<(Instant, Option<f64>)>> {
    let reduced = reduce(
        product,
        as_of.parse()?,
        fixings,
        &[],
        SettlementBoundary::Include,
    )?;
    let mut work = vec![reduced.contract.root];
    let mut payments = Vec::new();
    while let Some(node) = work.pop() {
        match &reduced.contract.contracts[node as usize] {
            ContractNode::All(children) => work.extend(children),
            ContractNode::Pay { asset, amount, at } => {
                assert_eq!(asset, &Asset::Cash(ccy));
                let value = match reduced.contract.exprs[*amount as usize] {
                    ExprNode::Const(value) => Some(value.value()),
                    _ => None,
                };
                payments.push((*at, value));
            }
            _ => return Err(Error::Internal("expected a linear cashflow".into())),
        }
    }
    Ok(payments)
}

fn assert_payments(
    product: &Product,
    as_of: &str,
    fixings: &FixingStore,
    ccy: Currency,
    expected: &[(&str, f64)],
) -> Result<()> {
    let mut actual = reduced_payments(product, as_of, fixings, ccy)?
        .into_iter()
        .map(|(at, value)| {
            value
                .map(|value| (at, value))
                .ok_or_else(|| Error::Internal("expected a resolved cashflow amount".into()))
        })
        .collect::<Result<Vec<_>>>()?;
    actual.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.total_cmp(&b.1)));
    let mut expected = expected
        .iter()
        .map(|&(at, amount)| -> Result<(Instant, f64)> { Ok((at.parse()?, amount)) })
        .collect::<Result<Vec<_>>>()?;
    expected.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.total_cmp(&b.1)));
    assert_eq!(actual.len(), expected.len());
    let absolute = if ccy == Currency::USD { 1e-8 } else { 1e-6 };
    for ((actual_date, actual), (expected_date, expected)) in actual.into_iter().zip(expected) {
        assert_eq!(actual_date, expected_date);
        assert_close!(actual, expected, 1e-12, absolute);
    }
    Ok(())
}

#[test]
fn deposit_uses_spot_and_month_end_to_compute_principal_and_interest() -> Result<()> {
    let refdata = refdata()?;
    let deposit = StandardProduct::Deposit(Deposit::new(
        Money::new(Currency::JPY, 3650.0),
        0.10,
        DepositDates::Market("1M".parse()?),
        Some("JPY-DEP".parse()?),
        None,
        Position::Long,
        PaymentRounding::None,
    )?);
    assert_eq!(deposit.kind(), ProductKind::Deposit);
    let product = deposit.to_product("2027-02-24".parse()?, &refdata)?;
    let fixings = fixings(&refdata, "JPY-TONA", &[])?;
    // 水曜の 2 営業日後は 2 月の最終営業日。1 か月後は 3 月 31 日で利息は 33 日分。
    assert_payments(
        &product,
        "2027-02-24T00:00:00Z",
        &fixings,
        Currency::JPY,
        &[
            ("2027-02-26T06:00:00Z", -3650.0),
            ("2027-03-31T06:00:00Z", 3683.0),
        ],
    )
}

#[test]
fn deposit_short_tenors_cross_the_weekend_and_reverse_for_borrowing() -> Result<()> {
    let refdata = refdata()?;
    let fixings = fixings(&refdata, "JPY-TONA", &[])?;
    // ON、TN、SN はそれぞれ当日、翌営業日、スポット日から金曜〜月曜の 3 日を数える。
    for (trade, tenor) in [
        ("2027-01-08", "ON"),
        ("2027-01-07", "TN"),
        ("2027-01-06", "SN"),
    ] {
        let product = Deposit::new(
            Money::new(Currency::JPY, 3650.0),
            0.10,
            DepositDates::Market(tenor.parse()?),
            Some("JPY-DEP".parse()?),
            None,
            Position::Short,
            PaymentRounding::None,
        )?
        .to_product(trade.parse()?, &refdata)?;
        assert_payments(
            &product,
            "2027-01-06T00:00:00Z",
            &fixings,
            Currency::JPY,
            &[
                ("2027-01-08T06:00:00Z", 3650.0),
                ("2027-01-11T06:00:00Z", -3653.0),
            ],
        )?;
    }
    Ok(())
}

#[test]
fn confirmed_deposit_keeps_the_holiday_and_retains_its_warning_after_reduction() -> Result<()> {
    let refdata = refdata()?;
    let start = "2027-01-09".parse()?;
    let fixings = fixings(&refdata, "JPY-TONA", &[])?;
    let product = Deposit::new(
        Money::new(Currency::JPY, 3650.0),
        0.10,
        DepositDates::Confirmed {
            start,
            maturity: "2027-01-12".parse()?,
        },
        None,
        Some(DayCount::Act365F),
        Position::Long,
        PaymentRounding::None,
    )?
    .to_product(start, &refdata)?;
    let as_of = "2027-01-08T00:00:00Z";
    let expected = [Warning::ExplicitDateIsHoliday { date: start }];
    assert_eq!(product.warnings(), &expected);
    assert_eq!(
        reduce(
            &product,
            as_of.parse()?,
            &fixings,
            &[],
            SettlementBoundary::Include
        )?
        .warnings(),
        &expected
    );
    // 土曜に元本を出し、火曜までの 3 日に対する利息を払う。
    assert_payments(
        &product,
        as_of,
        &fixings,
        Currency::JPY,
        &[
            ("2027-01-09T06:00:00Z", -3650.0),
            ("2027-01-12T06:00:00Z", 3653.0),
        ],
    )
}

#[test]
fn tona_ois_matches_weekend_compounding_and_payment_lag() -> Result<()> {
    let refdata = refdata()?;
    let standard = StandardProduct::Swap(Swap::from_convention(
        "JPY-TONA-OIS".parse()?,
        None,
        EndRule::Date("2027-01-12".parse()?),
        Money::new(Currency::JPY, 3650.0),
        0.05,
        PayReceive::Receive,
        Spread::new(0.0)?,
        PaymentRounding::None,
    )?);
    assert_eq!(standard.kind(), ProductKind::Swap);
    let product = standard.to_product("2027-01-06".parse()?, &refdata)?;
    let fixings = fixings(
        &refdata,
        "JPY-TONA",
        &[("2027-01-08", 0.10), ("2027-01-11", 0.20)],
    )?;
    // 2 営業日後の金曜に始まり、金曜の fixing は 3 日分、月曜は 1 日分。両脚は満期の 2 営業日後払い。
    let compounded = 3650.0 * ((1.0 + 0.10 * 3.0 / 365.0) * (1.0 + 0.20 / 365.0) - 1.0);
    assert_payments(
        &product,
        "2027-01-12T02:00:00Z",
        &fixings,
        Currency::JPY,
        &[
            ("2027-01-14T06:00:00Z", 2.0),
            ("2027-01-14T06:00:00Z", -compounded),
        ],
    )
}

#[test]
fn tona_ois_puts_the_odd_period_first_like_backward_generation() -> Result<()> {
    let refdata = refdata()?;
    let product = Swap::from_convention(
        "JPY-TONA-OIS".parse()?,
        None,
        EndRule::AfterStart("18M".parse()?),
        Money::new(Currency::JPY, 3650.0),
        0.05,
        PayReceive::Receive,
        Spread::new(0.0)?,
        PaymentRounding::None,
    )?
    .to_product("2027-01-06".parse()?, &refdata)?;
    let fixings = fixings(&refdata, "JPY-TONA", &[])?;
    // スポット 2027-01-08（金）から 18M 後の 2028-07-08（土）は修正翌営業日で 07-10（月）。
    // 満期から 1 年遡った 2027-07-08（木）が境界で、先頭が 181 日、次が閏日を含む 368 日。
    // 開始日から進める生成なら境界は 2028-01-10 になり、支払日も額も変わる。
    // 固定脚は 3650 × 5% × 日数 / 365 を受け取り、変動脚は未確定の複利を同じ日に払う。
    let mut actual = reduced_payments(&product, "2027-01-06T00:00:00Z", &fixings, Currency::JPY)?;
    actual.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1.is_some().cmp(&b.1.is_some()))
    });
    let first: Instant = "2027-07-12T06:00:00Z".parse()?;
    let second: Instant = "2028-07-12T06:00:00Z".parse()?;
    assert_eq!(actual.len(), 4);
    for ((at, amount), (expected_at, expected)) in actual.into_iter().zip([
        (first, None),
        (first, Some(90.5)),
        (second, None),
        (second, Some(184.0)),
    ]) {
        assert_eq!(at, expected_at);
        match (amount, expected) {
            (Some(amount), Some(expected)) => assert_close!(amount, expected, 1e-12, 1e-9),
            (None, None) => {}
            _ => panic!("unexpected cashflow at {at:?}: {amount:?}"),
        }
    }
    Ok(())
}

fn period(start: &str, end: &str, payment: &str) -> Result<ExplicitPeriod> {
    Ok(ExplicitPeriod {
        accrual: DateInterval::new(start.parse()?, end.parse()?)?,
        payment: payment.parse()?,
        reset: None,
        settlement: None,
        references: Vec::new(),
    })
}

fn explicit(periods: Vec<ExplicitPeriod>, combine: PeriodCombination) -> Result<ScheduleSpec> {
    Ok(ScheduleSpec::Explicit(ExplicitPeriods::new(
        periods,
        ExplicitDatePolicy::KeepAsConfirmed,
        combine,
        None,
    )?))
}

fn fixed_leg(
    notional: Money,
    schedule: ScheduleSpec,
    basis: DayCount,
    rounding: PaymentRounding,
) -> Result<Leg> {
    Leg::new(
        PayReceive::Pay,
        notional,
        schedule,
        basis.clone(),
        Coupon::Fixed {
            rate: Rate::new(0.15, basis, Simple)?,
        },
        NotionalExchange::None,
        rounding,
    )
}

fn floating_leg(
    notional: Money,
    schedule: ScheduleSpec,
    index: &str,
    rule: ObservationRule,
    spread: f64,
) -> Result<Leg> {
    Leg::new(
        PayReceive::Receive,
        notional,
        schedule,
        DayCount::Act365F,
        Coupon::Compounded {
            index: index.parse()?,
            rule,
            spread: Spread::new(spread)?,
        },
        NotionalExchange::None,
        PaymentRounding::None,
    )
}

#[test]
fn sofr_swap_uses_shifted_annualization_and_rounds_the_annual_rate() -> Result<()> {
    let refdata = refdata()?;
    let schedule = explicit(
        vec![period("2027-07-02", "2027-07-06", "2027-07-07")?],
        PeriodCombination::Sum,
    )?;
    let notional = Money::new(Currency::USD, 3650.0);
    let fixed = fixed_leg(
        notional,
        schedule.clone(),
        DayCount::Act365F,
        PaymentRounding::None,
    )?;
    let rounding = Rounding::new(0.01, RoundingMode::HalfEven)?;
    let floating = floating_leg(
        notional,
        schedule,
        "USD-SOFR",
        rule(
            CalendarSpec::Base("USNY".parse()?),
            ObservationMethod::ObservationShift(NonZeroU8::MIN),
            SpreadTreatment::Exclusive,
            ObservationRounding {
                daily: Some(rounding),
                aggregated: Some(rounding),
            },
        )?,
        0.017,
    )?;
    let product =
        Swap::from_legs([fixed, floating])?.to_product("2027-07-01".parse()?, &refdata)?;
    let fixings = fixings(&refdata, "USD-SOFR", &[("2027-07-01", 0.104)])?;
    // 7/5 が休日なので観測は 7/1〜7/2 の 1 日。0.104 → 0.10、スプレッド込み 0.117 → 0.12。
    // 受渡額は指標の ACT/360 で年率化した後、脚の ACT/365F で 4 日分を数える。
    assert_payments(
        &product,
        "2027-07-02T15:00:00Z",
        &fixings,
        Currency::USD,
        &[
            ("2027-07-07T19:00:00Z", -6.0),
            ("2027-07-07T19:00:00Z", 4.8),
        ],
    )
}

#[test]
fn rfr_spread_is_applied_inside_or_outside_each_daily_factor() -> Result<()> {
    let refdata = refdata()?;
    let schedule = explicit(
        vec![period("2027-01-08", "2027-01-12", "2027-01-14")?],
        PeriodCombination::Sum,
    )?;
    let notional = Money::new(Currency::JPY, 3650.0);
    let fixings = fixings(
        &refdata,
        "JPY-TONA",
        &[("2027-01-08", 0.10), ("2027-01-11", 0.20)],
    )?;
    let spread = 0.05;
    // 金曜と月曜の日次因子の積に、スプレッドを加える位置だけを変えて比べる。
    let exclusive = 3650.0 * ((1.0 + 0.10 * 3.0 / 365.0) * (1.0 + 0.20 / 365.0) - 1.0) + 2.0;
    let inclusive = 3650.0 * ((1.0 + 0.15 * 3.0 / 365.0) * (1.0 + 0.25 / 365.0) - 1.0);
    for (treatment, expected) in [
        (SpreadTreatment::Exclusive, exclusive),
        (SpreadTreatment::Inclusive, inclusive),
    ] {
        let fixed = fixed_leg(
            notional,
            schedule.clone(),
            DayCount::Act365F,
            PaymentRounding::None,
        )?;
        let floating = floating_leg(
            notional,
            schedule.clone(),
            "JPY-TONA",
            rule(
                CalendarSpec::Base("JPTO".parse()?),
                ObservationMethod::Plain,
                treatment,
                ObservationRounding::default(),
            )?,
            spread,
        )?;
        let product =
            Swap::from_legs([fixed, floating])?.to_product("2027-01-07".parse()?, &refdata)?;
        assert_payments(
            &product,
            "2027-01-12T02:00:00Z",
            &fixings,
            Currency::JPY,
            &[
                ("2027-01-14T06:00:00Z", -6.0),
                ("2027-01-14T06:00:00Z", expected),
            ],
        )?;
    }
    Ok(())
}

#[test]
fn grouped_coupons_are_summed_before_currency_rounding() -> Result<()> {
    let refdata = refdata()?;
    let periods = vec![
        period("2027-01-07", "2027-01-08", "2027-01-11")?,
        period("2027-01-08", "2027-01-09", "2027-01-11")?,
    ];
    let schedule = explicit(periods, PeriodCombination::Sum)?;
    let notional = Money::new(Currency::JPY, 1460.0);
    let fixed = fixed_leg(
        notional,
        schedule.clone(),
        DayCount::Act365F,
        PaymentRounding::CurrencyConvention,
    )?;
    let floating = floating_leg(
        notional,
        schedule,
        "JPY-TONA",
        rule(
            CalendarSpec::Base("JPTO".parse()?),
            ObservationMethod::Plain,
            SpreadTreatment::Exclusive,
            ObservationRounding::default(),
        )?,
        0.0,
    )?;
    let product =
        Swap::from_legs([fixed, floating])?.to_product("2027-01-06".parse()?, &refdata)?;
    let fixings = fixings(
        &refdata,
        "JPY-TONA",
        &[("2027-01-07", 0.10), ("2027-01-08", 0.10)],
    )?;
    // 固定脚の 0.6 円を 2 回合算し、1.2 円を切り捨てる。RFR 脚も各計算期間を単純合算する。
    assert_payments(
        &product,
        "2027-01-11T02:00:00Z",
        &fixings,
        Currency::JPY,
        &[
            ("2027-01-11T06:00:00Z", -1.0),
            ("2027-01-11T06:00:00Z", 0.8),
        ],
    )
}

#[test]
fn published_tona_interest_keeps_the_half_yen_at_the_rounding_boundary() -> Result<()> {
    // 元本 73,109,500 円、ACT/365F の 1 日、公表済みの 0.5% の利息はちょうど 1,001.5 円で、
    // 円単位の HalfEven では偶数の 1,002 円になる。成長率 1 + 0.005/365 を丸めてから 1 を引くと
    // 超過の相対誤差が約 7e-12 になり 1,001.4999999932 円、つまり 1,001 円になってしまう
    // （確定分の超過成長を積に丸めずに求めることの回帰テスト）。確定因子の超過は e ← e + x + e·x の
    // 四則演算だけで畳むので、1 日分の超過は r·w そのもので、ちょうど半分の値も環境（libm）によらない。
    let refdata = refdata_rounding(RoundingMode::HalfEven)?;
    let schedule = explicit(
        vec![period("2027-01-07", "2027-01-08", "2027-01-08")?],
        PeriodCombination::Sum,
    )?;
    let notional = Money::new(Currency::JPY, 73_109_500.0);
    let fixings = fixings(&refdata, "JPY-TONA", &[("2027-01-07", 0.005)])?;
    for (rounding, expected) in [
        (PaymentRounding::None, 1001.5),
        (PaymentRounding::CurrencyConvention, 1002.0),
    ] {
        let fixed = fixed_leg(notional, schedule.clone(), DayCount::Act365F, rounding)?;
        let floating = Leg::new(
            PayReceive::Receive,
            notional,
            schedule.clone(),
            DayCount::Act365F,
            Coupon::Compounded {
                index: "JPY-TONA".parse()?,
                rule: rule(
                    CalendarSpec::Base("JPTO".parse()?),
                    ObservationMethod::Plain,
                    SpreadTreatment::Exclusive,
                    ObservationRounding::default(),
                )?,
                spread: Spread::new(0.0)?,
            },
            NotionalExchange::None,
            rounding,
        )?;
        let product =
            Swap::from_legs([fixed, floating])?.to_product("2027-01-06".parse()?, &refdata)?;
        // 支払いの瞬間ちょうどに評価し、Include で受渡を残す。固定脚は 0.15 の 1 日分 30,045 円。
        let mut payments =
            reduced_payments(&product, "2027-01-08T06:00:00Z", &fixings, Currency::JPY)?
                .into_iter()
                .map(|(_, value)| {
                    value.ok_or_else(|| Error::Internal("expected a fixed amount".into()))
                })
                .collect::<Result<Vec<_>>>()?;
        payments.sort_by(f64::total_cmp);
        let [paid, received] = payments.as_slice() else {
            return Err(Error::Internal("expected two payments".into()));
        };
        assert_close!(*paid, -30_045.0, 1e-12, 1e-6);
        // 丸める前の額は 1,001.5（丸めの境界）を保つ。
        assert_close!(*received, expected, 1e-15, 0.0);
    }
    Ok(())
}

#[test]
fn swap_rejects_missing_accrual_context_invalid_currency_and_grouping() -> Result<()> {
    let refdata = refdata()?;
    let trade: Date = "2027-01-07".parse()?;
    let periods = vec![period("2027-01-08", "2027-01-12", "2027-01-14")?];
    let schedule = explicit(periods.clone(), PeriodCombination::Sum)?;
    let notional = Money::new(Currency::JPY, 3650.0);
    let floating = |ccy, index: &str| {
        floating_leg(
            Money::new(ccy, 3650.0),
            schedule.clone(),
            index,
            rule(
                CalendarSpec::Base("JPTO".parse()?),
                ObservationMethod::Plain,
                SpreadTreatment::Exclusive,
                ObservationRounding::default(),
            )?,
            0.0,
        )
    };

    // 明示日程には ICMA と ACT/365L の枠がなく、30/360 US の月末慣行も分からない。
    assert!(matches!(
        fixed_leg(
            notional,
            schedule.clone(),
            DayCount::Thirty360Us,
            PaymentRounding::None
        ),
        Err(Error::InvalidData(_))
    ));
    for basis in [DayCount::ActActIcma, DayCount::Act365L] {
        let fixed = fixed_leg(notional, schedule.clone(), basis, PaymentRounding::None)?;
        assert!(
            matches!(Swap::from_legs([fixed, floating(Currency::JPY, "JPY-TONA")?])?.to_product(trade, &refdata),
            Err(Error::InvalidData(message)) if message.contains("reference period"))
        );
    }

    let fixed = fixed_leg(
        notional,
        schedule.clone(),
        DayCount::Act365F,
        PaymentRounding::None,
    )?;
    assert!(matches!(
        Swap::from_legs([fixed.clone(), floating(Currency::USD, "USD-SOFR")?]),
        Err(Error::InvalidData(_))
    ));
    assert!(
        matches!(Swap::from_legs([fixed, floating(Currency::JPY, "USD-SOFR")?])?.to_product(trade, &refdata),
        Err(Error::InvalidData(message)) if message.contains("index USD-SOFR currency"))
    );

    // 期間同士を複利にする指定を単純合算として受け取ると、金額が変わるので拒否する。
    let schedule = explicit(periods, PeriodCombination::StraightCompounding)?;
    let fixed = fixed_leg(notional, schedule, DayCount::Act365F, PaymentRounding::None)?;
    assert!(matches!(
        Swap::from_legs([fixed, floating(Currency::JPY, "JPY-TONA")?])?.to_product(trade, &refdata),
        Err(Error::Unsupported(_))
    ));
    Ok(())
}

#[test]
fn notional_exchanges_use_period_boundaries_instead_of_coupon_paydays() -> Result<()> {
    let refdata = refdata()?;
    let notional = Money::new(Currency::JPY, 3650.0);
    let fixed = Leg::new(
        PayReceive::Pay,
        notional,
        explicit(
            vec![period("2027-01-11", "2027-01-12", "2027-01-14")?],
            PeriodCombination::Sum,
        )?,
        DayCount::Act365F,
        Coupon::Fixed {
            rate: Rate::new(0.1, DayCount::Act365F, Simple)?,
        },
        NotionalExchange::InitialAndFinal,
        PaymentRounding::None,
    )?;
    let floating = floating_leg(
        notional,
        explicit(
            vec![period("2027-01-07", "2027-01-08", "2027-01-13")?],
            PeriodCombination::Sum,
        )?,
        "JPY-TONA",
        rule(
            CalendarSpec::Base("JPTO".parse()?),
            ObservationMethod::Plain,
            SpreadTreatment::Exclusive,
            ObservationRounding::default(),
        )?,
        0.0,
    )?;
    let product =
        Swap::from_legs([fixed, floating])?.to_product("2027-01-06".parse()?, &refdata)?;
    let fixings = fixings(&refdata, "JPY-TONA", &[("2027-01-07", 0.1)])?;
    // 固定金利を払う脚は期首に元本を受け取り、期末に返す。二つのクーポンは別の日に支払う。
    assert_payments(
        &product,
        "2027-01-08T02:00:00Z",
        &fixings,
        Currency::JPY,
        &[
            ("2027-01-11T06:00:00Z", 3650.0),
            ("2027-01-12T06:00:00Z", -3650.0),
            ("2027-01-13T06:00:00Z", 1.0),
            ("2027-01-14T06:00:00Z", -1.0),
        ],
    )
}

#[test]
fn instrument_template_builds_a_unit_payer_ois_and_rejects_other_conventions() -> Result<()> {
    // 校正のテンプレートは慣行どおりの OIS を、元本 1・固定払い・スプレッドなしで作る。
    // パーは固定金利について affine な PV から求めるので、quote の値は固定金利として入る。
    let refdata = refdata()?;
    let tenor: rhodium_core::tenor::TenorLength = "1Y".parse()?;
    let template = InstrumentTemplate::Ois {
        convention: "JPY-TONA-OIS".parse()?,
        tenor,
    };
    assert_eq!(template.currency(&refdata)?, Currency::JPY);
    template.validate(&refdata)?;
    let expected = StandardProduct::Swap(Swap::from_convention(
        "JPY-TONA-OIS".parse()?,
        None,
        EndRule::AfterStart(tenor),
        Money::new(Currency::JPY, 1.0),
        0.0125,
        PayReceive::Pay,
        Spread::new(0.0)?,
        PaymentRounding::None,
    )?);
    assert_eq!(template.build(0.0125, &refdata)?, expected);

    // 預金の慣行は OIS の校正商品にできない。データの誤りとして報告する。
    let deposit = InstrumentTemplate::Ois {
        convention: "JPY-DEP".parse()?,
        tenor,
    };
    assert!(matches!(
        deposit.currency(&refdata),
        Err(Error::InvalidData(_))
    ));
    assert!(matches!(
        deposit.validate(&refdata),
        Err(Error::InvalidData(_))
    ));
    assert!(matches!(
        deposit.build(0.0125, &refdata),
        Err(Error::InvalidData(_))
    ));
    Ok(())
}
`````

`````
crates/rhodium-products/src/ir/validate.rs
use std::{
    collections::{BTreeSet, HashMap, hash_map::Entry},
    sync::{Arc, Mutex, OnceLock},
};

use rhodium_core::{
    Error, Result, Warning,
    currency::Asset,
    date::DateInterval,
    fixing::ObservationPoint,
    id::IndexId,
    refdata::{IndexDefinition, IndexKind, IndexStatus, RefData},
    time::Instant,
};

use super::{
    Contract, References, ResolvedCompounding, ResolvedCredit, ResolvedFixing,
    nodes::{ContractNode, ExprNode, MAX_NODES, Nary, RfrRows, node_limit},
    reachable_contracts,
};

type CachedCompounding = Arc<OnceLock<Result<Arc<ResolvedCompounding>>>>;

// 同じ参照データで取引を解決する間だけ保持する。公表時刻の計算は表のロックを外して行う。
#[derive(Debug, Default)]
pub(crate) struct ObservationCache(Mutex<HashMap<(IndexId, Arc<RfrRows>), CachedCompounding>>);

impl ObservationCache {
    fn resolve(
        &self,
        rows: &mut Arc<RfrRows>,
        definition: &IndexDefinition,
        refdata: &RefData,
    ) -> Result<Arc<ResolvedCompounding>> {
        let entry = {
            let mut entries = self
                .0
                .lock()
                .map_err(|_| Error::Internal("RFR observation cache lock poisoned".into()))?;
            match entries.entry((definition.id.clone(), Arc::clone(rows))) {
                Entry::Occupied(entry) => {
                    *rows = Arc::clone(&entry.key().1);
                    Arc::clone(entry.get())
                }
                Entry::Vacant(entry) => Arc::clone(entry.insert(Arc::default())),
            }
        };
        entry
            .get_or_init(|| resolve_compounding(rows, definition, refdata))
            .clone()
    }
}

pub(super) fn validate(
    contract: &mut Contract,
    refdata: &RefData,
    cache: Option<&ObservationCache>,
) -> Result<(References, Vec<Warning>)> {
    if contract.contracts.len() + contract.exprs.len() > MAX_NODES {
        return Err(node_limit());
    }
    let needed = reachable_contracts(&contract.contracts, contract.root);
    let exprs_needed = contract.reachable_exprs(&needed);

    let mut references = References::default();
    let mut known = vec![None; exprs_needed.len()];
    for (pos, (node, needed)) in contract.exprs.iter_mut().zip(exprs_needed).enumerate() {
        if needed {
            known[pos] = known_at(node, &known, &mut references, refdata, cache)
                .map_err(|error| error.context(format_args!("expression {pos}")))?;
        }
    }

    let mut earliest = vec![None; needed.len()];
    for (pos, (node, needed)) in contract.contracts.iter().zip(needed).enumerate() {
        if needed {
            earliest[pos] = check_contract(node, &known, &earliest, &mut references, refdata)
                .map_err(|error| error.context(format_args!("contract {pos}")))?;
        }
    }

    let warnings = references
        .fixings
        .values()
        .flat_map(|items| items.values())
        .flat_map(|item| &item.warnings)
        .chain(
            references
                .compounded
                .values()
                .flat_map(|items| items.values())
                .flat_map(|item| &item.warnings),
        )
        .chain(
            references
                .credit
                .values()
                .flat_map(|items| items.values())
                .flat_map(|item| &item.warnings),
        )
        .cloned()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    Ok((references, warnings))
}

fn known_at(
    node: &mut ExprNode,
    known: &[Option<Instant>],
    references: &mut References,
    refdata: &RefData,
    cache: Option<&ObservationCache>,
) -> Result<Option<Instant>> {
    let defaulted = matches!(node, ExprNode::Defaulted { .. });
    match node {
        ExprNode::Const(value) if !value.value().is_finite() => {
            Err(invalid("constant must be finite"))
        }
        ExprNode::Nary(Nary::Min | Nary::Max, xs) if xs.is_empty() => {
            Err(invalid("min and max require at least one operand"))
        }
        ExprNode::Fixing {
            index,
            point,
            fx_reset_at,
        } => {
            let definition = refdata.index(index)?;
            check_numeric_index(&definition.kind, index)?;
            let publication = resolve_fixing(references, refdata, definition, *point)?;
            if let Some(at) = fx_reset_at
                && (!matches!(definition.kind, IndexKind::FxFixing { .. }) || *at < publication)
            {
                return Err(invalid(
                    "FX reset requires an FX fixing known by its value time",
                ));
            }
            Ok(Some(publication))
        }
        ExprNode::Compounded { index, rows, .. } => {
            let definition = refdata.index(index)?;
            check_numeric_index(&definition.kind, index)?;
            if !matches!(definition.kind, IndexKind::OvernightRate { .. }) {
                return Err(invalid(format!(
                    "compounding requires an overnight rate, got {index}"
                )));
            }
            let entries = references.compounded.entry(index.clone()).or_default();
            if let Some((source, resolved)) = entries.get_key_value(&rows.source) {
                rows.source = Arc::clone(source);
                return Ok(resolved.latest);
            }
            let resolved = match cache {
                Some(cache) => cache.resolve(&mut rows.source, definition, refdata)?,
                None => resolve_compounding(&rows.source, definition, refdata)?,
            };
            let latest = resolved.latest;
            entries.insert(Arc::clone(&rows.source), resolved);
            Ok(latest)
        }
        ExprNode::Average { index, rows, .. } => {
            let definition = refdata.index(index)?;
            check_numeric_index(&definition.kind, index)?;
            rows.iter().try_fold(None, |latest, &(point, _)| {
                Ok(latest.max(Some(resolve_fixing(
                    references, refdata, definition, point,
                )?)))
            })
        }
        ExprNode::Defaulted { index, window } | ExprNode::Recovery { index, window } => {
            let definition = refdata.index(index)?;
            let ((true, IndexKind::CreditEvent { credit: key })
            | (false, IndexKind::Recovery { credit: key })) = (defaulted, &definition.kind)
            else {
                return Err(invalid(format!("wrong credit index kind: {index}")));
            };
            if window.is_empty() {
                return Err(invalid("credit event window is empty"));
            }
            Ok(Some(resolve_credit(
                references, refdata, index, *window, key,
            )?))
        }
        _ => Ok(node
            .children()
            .iter()
            .filter_map(|&ix| known[ix as usize])
            .max()),
    }
}

fn check_numeric_index(kind: &IndexKind, index: &IndexId) -> Result<()> {
    if matches!(
        kind,
        IndexKind::CreditEvent { .. } | IndexKind::Recovery { .. }
    ) {
        return Err(invalid(format!(
            "credit index {index} requires a credit window"
        )));
    }
    Ok(())
}

fn resolve_fixing(
    references: &mut References,
    refdata: &RefData,
    definition: &IndexDefinition,
    point: ObservationPoint,
) -> Result<Instant> {
    let id = &definition.id;
    let observations = references.fixings.entry(id.clone()).or_default();
    match observations.entry(point) {
        Entry::Occupied(entry) => Ok(entry.get().publication),
        Entry::Vacant(entry) => Ok(entry
            .insert(fixing(refdata, definition, point)?)
            .publication),
    }
}

fn resolve_compounding(
    rows: &RfrRows,
    definition: &IndexDefinition,
    refdata: &RefData,
) -> Result<Arc<ResolvedCompounding>> {
    let observations = rows
        .entries
        .iter()
        .map(|&((point, _), _)| fixing(refdata, definition, point))
        .collect::<Result<Box<[_]>>>()?;
    let latest = observations.iter().map(|item| item.publication).max();
    let warnings = observations
        .iter()
        .flat_map(|item| item.warnings.iter().cloned())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    Ok(Arc::new(ResolvedCompounding {
        observations,
        latest,
        warnings,
    }))
}

fn fixing(
    refdata: &RefData,
    definition: &IndexDefinition,
    point: ObservationPoint,
) -> Result<ResolvedFixing> {
    let id = &definition.id;
    let (publication, mut warnings) = refdata.publication(id, point)?;
    let (period, notices) = refdata.observation_period(id, point)?;
    warnings.extend(notices);
    Ok(ResolvedFixing {
        publication,
        period,
        warnings,
        ceased: matches!((definition.status, point),
            (IndexStatus::Ceased { since }, ObservationPoint::Date(date)) if date >= since),
        non_representative: matches!((definition.status, point),
            (IndexStatus::NonRepresentative { since }, ObservationPoint::Date(date)) if date >= since),
        kind: match (&definition.kind, point) {
            (IndexKind::TermRate { .. }, ObservationPoint::Date(_)) => super::FixingKind::TermRate,
            (IndexKind::FxFixing { pair }, _) => super::FixingKind::Fx(*pair),
            (IndexKind::RateFuturesSettlement { .. }, _) => super::FixingKind::RateFuture,
            _ => super::FixingKind::Other,
        },
    })
}

fn resolve_credit(
    references: &mut References,
    refdata: &RefData,
    id: &IndexId,
    window: DateInterval,
    key: &rhodium_core::id::CreditKey,
) -> Result<Instant> {
    match references
        .credit
        .entry(id.clone())
        .or_default()
        .entry(window)
    {
        Entry::Occupied(entry) => Ok(entry.get().publication),
        Entry::Vacant(entry) => {
            let (publication, warnings) = refdata.credit_publication(id, window)?;
            Ok(entry
                .insert(ResolvedCredit {
                    key: key.clone(),
                    publication,
                    warnings,
                })
                .publication)
        }
    }
}

fn check_contract(
    node: &ContractNode,
    known: &[Option<Instant>],
    earliest: &[Option<Instant>],
    references: &mut References,
    refdata: &RefData,
) -> Result<Option<Instant>> {
    let at = match node {
        ContractNode::Pay { asset, amount, at } => {
            match asset {
                Asset::Cash(ccy) => {
                    refdata.currency(*ccy)?;
                }
                Asset::Security(id) => {
                    refdata.security(id)?;
                }
                Asset::Commodity(id) => {
                    refdata.commodity(id)?;
                }
            }
            check_known(known[*amount as usize], *at, "payment amount")?;
            *at
        }
        ContractNode::All(children) => {
            return Ok(children
                .iter()
                .filter_map(|&ix| earliest[ix as usize])
                .min());
        }
        ContractNode::Cond { cond, at, .. } => {
            check_known(known[*cond as usize], *at, "condition")?;
            *at
        }
        ContractNode::Choose {
            label, at, options, ..
        } => {
            if options.is_empty() {
                return Err(invalid(format!("choice {label:?} at {at} has no options")));
            }
            match references
                .choices
                .entry(label.clone())
                .or_default()
                .entry(*at)
            {
                Entry::Occupied(_) => {
                    return Err(invalid(format!("duplicate choice {label:?} at {at}")));
                }
                Entry::Vacant(entry) => {
                    entry.insert(options.len());
                }
            }
            *at
        }
    };
    let mut first = None;
    node.visit_children(|ix| {
        if let Some(instant) = earliest[ix as usize] {
            first = Some(first.map_or(instant, |current: Instant| current.min(instant)));
        }
    });
    if let Some(first) = first.filter(|instant| *instant < at) {
        return Err(invalid(format!(
            "branch event at {first} precedes decision at {at}"
        )));
    }
    Ok(Some(at))
}

fn check_known(known: Option<Instant>, at: Instant, what: &str) -> Result<()> {
    if let Some(known) = known.filter(|known| *known > at) {
        return Err(invalid(format!("{what} is known at {known}, after {at}")));
    }
    Ok(())
}

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidInput(message.into())
}
`````

`````
crates/rhodium-products/src/lib.rs
//! コア IR、DSL、標準商品と取引の定義。

#![forbid(unsafe_code)]

pub mod agreement;
pub mod ir;
pub mod standard;
pub mod trade;
`````

`````
crates/rhodium-products/src/standard.rs
//! 確認書の条件を持つ標準商品と、検証済み IR への変換。

mod fra;
mod future;
mod fx;
pub mod leg;
mod rates;
mod template;
mod wire;

use rhodium_core::{
    Error, Result, Warning,
    date::Date,
    money::Money,
    refdata::{CurrencyDefinition, RefData},
    time::Instant,
};

use crate::ir::{Builder, Contract, ExprIx, Product};

pub use fra::Fra;
pub use future::{RateFuture, RateFutureMargin};
pub use fx::{FxForward, FxSwap};
pub use leg::{Coupon, Leg, Notional, NotionalExchange, PayReceive};
pub use rates::{Deposit, DepositDates, Swap};
pub use template::InstrumentTemplate;

/// 確認書の条件を持つ標準商品。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum StandardProduct {
    /// 預金。
    Deposit(Deposit),
    /// 期首決済の先渡金利契約。
    Fra(Fra),
    /// 日次値洗いの上場金利先物。
    RateFuture(RateFuture),
    /// OIS、IRS、ベーシススワップと通貨スワップ。
    Swap(Swap),
    /// 現物受渡の FX フォワード。
    FxForward(FxForward),
    /// 近い日と遠い日の二通貨交換。
    FxSwap(FxSwap),
}

impl StandardProduct {
    /// 価格設定が照合する商品の種類。
    #[must_use]
    pub const fn kind(&self) -> ProductKind {
        match self {
            Self::Deposit(_) => ProductKind::Deposit,
            Self::Fra(_) => ProductKind::Fra,
            Self::RateFuture(_) => ProductKind::RateFuture,
            Self::Swap(_) => ProductKind::Swap,
            Self::FxForward(_) => ProductKind::FxForward,
            Self::FxSwap(_) => ProductKind::FxSwap,
        }
    }

    /// 約定日と参照データから日程と観測を解決し、契約を IR にする。
    ///
    /// # Errors
    ///
    /// 参照データ、日程、受渡や観測の条件が不正か、計算に失敗した場合。
    pub fn to_product(&self, trade_date: Date, refdata: &RefData) -> Result<Product> {
        let (contract, warnings) = self.build_contract(trade_date, refdata)?;
        Product::with_construction_warnings(contract, refdata, warnings)
            .map_err(|error| error.context(self.label()))
    }

    pub(crate) fn build_contract(
        &self,
        trade_date: Date,
        refdata: &RefData,
    ) -> Result<(Contract, Vec<Warning>)> {
        let result = match self {
            Self::Deposit(deposit) => deposit.build_contract(trade_date, refdata),
            Self::Fra(fra) => fra.build_contract(refdata),
            Self::RateFuture(future) => future.build_contract(trade_date, refdata),
            Self::Swap(swap) => swap.build_contract(trade_date, refdata),
            Self::FxForward(forward) => forward.build_contract(refdata),
            Self::FxSwap(swap) => swap.build_contract(refdata),
        };
        result.map_err(|error| error.context(self.label()))
    }

    const fn label(&self) -> &'static str {
        match self {
            Self::Deposit(_) => "deposit",
            Self::Fra(_) => "fra",
            Self::RateFuture(_) => "rate_future",
            Self::Swap(_) => "swap",
            Self::FxForward(_) => "fx_forward",
            Self::FxSwap(_) => "fx_swap",
        }
    }
}

/// 標準商品を確認書どおりの向きで持つか、反対向きで持つか。
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Position {
    /// 預金では期首に元本を払い、満期に元利金を受け取る。
    #[default]
    Long,
    /// Long と反対方向に受け渡す。
    Short,
}

/// 受渡額に適用する丸め。
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PaymentRounding {
    /// 受渡額を丸めない。
    #[default]
    None,
    /// 支払通貨に定義された規則で丸める。
    CurrencyConvention,
}

impl PaymentRounding {
    fn apply<'b>(
        self,
        builder: &mut Builder<'b>,
        amount: ExprIx<'b>,
        currency: &CurrencyDefinition,
    ) -> ExprIx<'b> {
        match self {
            Self::None => amount,
            Self::CurrencyConvention => builder.round(amount, currency.payment_rounding),
        }
    }
}

fn validate_notional(notional: Money) -> Result<()> {
    if !notional.amount().is_finite() || notional.amount() <= 0.0 {
        return Err(invalid("notional must be positive and finite"));
    }
    Ok(())
}

fn settlement(date: Date, ccy: &CurrencyDefinition, refdata: &RefData) -> Result<Instant> {
    // RefData の構築が通貨の決済時刻の地域を検査済みなので、引けないことはない。
    let event = &ccy.settlement_time;
    refdata
        .zone(&event.zone)?
        .resolve_convention(event.on_date(date))
}

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidData(message.into())
}

/// 価格設定の方式の表が照合する標準商品の種類。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ProductKind {
    /// 預金。
    Deposit,
    /// 期首決済の先渡金利契約。
    Fra,
    /// 日次値洗いの上場金利先物。
    RateFuture,
    /// OIS を含む金利スワップ。
    Swap,
    /// 現物受渡の FX フォワード。
    FxForward,
    /// FX スワップ。
    FxSwap,
}
`````

`````
crates/rhodium-products/src/standard/fra.rs
//! 期首決済の先渡金利契約。

use rhodium_core::{
    Result, Warning,
    currency::Asset,
    date::Date,
    day_count::{AccrualContext, year_fraction},
    fixing::ObservationPoint,
    id::IndexId,
    money::Money,
    refdata::{IndexKind, RefData},
};

use crate::ir::{Contract, Product};

use super::{Position, invalid, leg::day_count_calendar, settlement, validate_notional};

/// ターム金利の fixing を用い、対象期間の開始日に差金決済する FRA。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(try_from = "super::wire::FraInput", into = "super::wire::FraInput")]
pub struct Fra {
    pub(super) notional: Money,
    pub(super) index: IndexId,
    pub(super) fixing: Date,
    pub(super) fixed_rate: f64,
    pub(super) position: Position,
}

impl Fra {
    /// 金利を小数で指定する。買いはターム金利受け・固定金利払い。
    ///
    /// # Errors
    ///
    /// 元本が正の有限値でないか、固定金利が有限でない場合。
    pub fn new(
        notional: Money,
        index: IndexId,
        fixing: Date,
        fixed_rate: f64,
        position: Position,
    ) -> Result<Self> {
        validate_notional(notional)?;
        if !fixed_rate.is_finite() {
            return Err(invalid("FRA fixed rate must be finite"));
        }
        Ok(Self {
            notional,
            index,
            fixing,
            fixed_rate,
            position,
        })
    }

    /// 指標の預金期間の期首に一回の差金決済を組み立てる。
    ///
    /// # Errors
    ///
    /// 指標と元本の通貨が合わないか、日程、日数係数、既知時刻を解決できない場合。
    pub fn to_product(&self, refdata: &RefData) -> Result<Product> {
        let (contract, warnings) = self.build_contract(refdata)?;
        Product::with_construction_warnings(contract, refdata, warnings)
    }

    pub(super) fn build_contract(&self, refdata: &RefData) -> Result<(Contract, Vec<Warning>)> {
        let IndexKind::TermRate { ccy, day_count, .. } = &refdata.index(&self.index)?.kind else {
            return Err(invalid(format!("{} is not a term rate", self.index)));
        };
        if *ccy != self.notional.ccy() {
            return Err(invalid("FRA index and notional currencies differ"));
        }
        let (period, mut warnings) = refdata.term_rate_dates(&self.index, self.fixing)?;
        let calendar = day_count_calendar(day_count, refdata)?;
        let context = AccrualContext {
            calendar: calendar.as_ref(),
            maturity: Some(period.end()),
            ..AccrualContext::default()
        };
        let (fraction, notices) = year_fraction(day_count, period.start(), period.end(), &context)?;
        warnings.extend(notices);
        let tau = fraction.value();
        if tau <= 0.0 {
            return Err(invalid("FRA requires a positive accrual fraction"));
        }
        let at = settlement(period.start(), refdata.currency(*ccy)?, refdata)?;
        let contract = Contract::build(|builder| {
            let floating = builder.fixing(&self.index, ObservationPoint::Date(self.fixing));
            let strike = builder.constant(self.fixed_rate);
            let difference = builder.sub(floating, strike);
            let notional_accrual = builder.constant(self.notional.amount() * tau);
            let numerator = builder.mul([notional_accrual, difference]);
            let tau = builder.constant(tau);
            let interest = builder.mul([tau, floating]);
            let one = builder.constant(1.0);
            let denominator = builder.add([one, interest]);
            let amount = builder.div(numerator, denominator);
            let pay = builder.pay(Asset::Cash(*ccy), amount, at);
            Ok(match self.position {
                Position::Long => pay,
                Position::Short => builder.give(pay),
            })
        })?;
        Ok((contract, warnings))
    }
}
`````

`````
crates/rhodium-products/src/standard/future.rs
//! 金利先物の日次値洗い。

use std::num::NonZeroU32;

use rhodium_core::{
    Result, Warning,
    currency::{Asset, Currency},
    date::Date,
    fixing::ObservationPoint,
    id::IndexId,
    quantity::RateFuturesPrice,
    refdata::{RateFutureClearing, RefData},
};

use crate::ir::{Contract, Product};

use super::{Position, invalid};

/// 確認書の値洗い条件を持つ上場金利先物。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(
    try_from = "super::wire::RateFutureInput",
    into = "super::wire::RateFutureInput"
)]
pub struct RateFuture {
    pub(super) settlement_index: IndexId,
    pub(super) contracts: NonZeroU32,
    pub(super) trade_price: RateFuturesPrice,
    pub(super) position: Position,
    pub(super) first_clearing: Option<Date>,
}

/// 確認書と清算慣行から解決した値洗い。
#[derive(Debug)]
pub struct RateFutureMargin<'a> {
    /// 限月固有の清算価格の指標。
    pub index: &'a IndexId,
    /// 値洗いの通貨。
    pub currency: Currency,
    /// 枚数と乗数の積。売り持ちでは負。
    pub signed_multiplier: f64,
    /// 初回値洗いの比較対象となる取引価格。
    pub trade_price: f64,
    /// 最初から最終までの清算日程。
    pub clearings: Vec<RateFutureClearing>,
}

impl RateFuture {
    /// 限月固有の清算価格指標、枚数と取引価格で組み立てる。
    ///
    /// 初回清算日を省略すると約定日以降の最初の清算日から値洗いを行う。
    #[must_use]
    pub const fn new(
        settlement_index: IndexId,
        contracts: NonZeroU32,
        trade_price: RateFuturesPrice,
        position: Position,
        first_clearing: Option<Date>,
    ) -> Self {
        Self {
            settlement_index,
            contracts,
            trade_price,
            position,
            first_clearing,
        }
    }

    /// 日次の清算と受渡を解決する。
    ///
    /// # Errors
    ///
    /// 限月、原指標、日程、取引価格の日数基準または乗数が不正な場合。
    pub fn margin(
        &self,
        trade_date: Date,
        refdata: &RefData,
    ) -> Result<(RateFutureMargin<'_>, Vec<Warning>)> {
        let (convention, _) = refdata.rate_future_contract(&self.settlement_index)?;
        let (_, currency, day_count) = refdata.rate_future_underlying(&self.settlement_index)?;
        if self.trade_price.basis() != day_count {
            return Err(invalid(
                "futures price day count differs from the underlying",
            ));
        }
        let start = self.first_clearing.unwrap_or(trade_date);
        if start < trade_date {
            return Err(invalid("first futures clearing precedes trade date"));
        }
        let (clearings, warnings) = refdata.rate_future_schedule(&self.settlement_index, start)?;
        if self.first_clearing.is_some() && clearings[0].date != start {
            return Err(invalid(
                "explicit first futures clearing is not a clearing day",
            ));
        }
        let quantity = f64::from(self.contracts.get()) * convention.multiplier;
        if !quantity.is_finite() {
            return Err(invalid(
                "futures multiplier times contracts is outside the finite range",
            ));
        }
        let signed_multiplier = match self.position {
            Position::Long => quantity,
            Position::Short => -quantity,
        };
        Ok((
            RateFutureMargin {
                index: &self.settlement_index,
                currency,
                signed_multiplier,
                trade_price: self.trade_price.value(),
                clearings,
            },
            warnings,
        ))
    }

    /// 取引価格から最終清算までの日次の差額を IR にする。
    ///
    /// # Errors
    ///
    /// 限月、原指標、日程、受渡額または IR が不正な場合。
    pub fn to_product(&self, trade_date: Date, refdata: &RefData) -> Result<Product> {
        let (contract, warnings) = self.build_contract(trade_date, refdata)?;
        Product::with_construction_warnings(contract, refdata, warnings)
    }

    pub(super) fn build_contract(
        &self,
        trade_date: Date,
        refdata: &RefData,
    ) -> Result<(Contract, Vec<Warning>)> {
        let (margin, warnings) = self.margin(trade_date, refdata)?;
        let contract = Contract::build(|builder| {
            let mut previous = builder.constant(margin.trade_price);
            let multiplier = builder.constant(margin.signed_multiplier.abs());
            let mut payments = Vec::with_capacity(margin.clearings.len());
            for row in &margin.clearings {
                let price = builder.fixing(margin.index, ObservationPoint::Date(row.date));
                let change = builder.sub(price, previous);
                let amount = builder.mul([multiplier, change]);
                payments.push(builder.pay(Asset::Cash(margin.currency), amount, row.settlement));
                previous = price;
            }
            let long = builder.all(payments);
            Ok(match self.position {
                Position::Long => long,
                Position::Short => builder.give(long),
            })
        })?;
        Ok((contract, warnings))
    }
}
`````

`````
crates/rhodium-products/src/standard/fx.rs
//! 二通貨を現物で受け渡す FX 商品。

use rhodium_core::{
    Result, Warning,
    currency::Asset,
    fx::{Forward, FxRate},
    money::Money,
    refdata::{Deliverability, RefData},
};

use crate::ir::{Builder, Contract, ContractIx, Product};

use super::{Position, invalid, validate_notional};

/// 同じ瞬間に二通貨を約定レートで交換する FX フォワード。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(
    try_from = "super::wire::FxForwardInput",
    into = "super::wire::FxForwardInput"
)]
pub struct FxForward {
    pub(super) notional: Money,
    pub(super) rate: FxRate<Forward>,
    pub(super) position: Position,
}

impl FxForward {
    /// Long はペアの base 通貨を受け取り、quote 通貨を支払う。
    ///
    /// # Errors
    ///
    /// 元本が正の有限値でない、ペアの base と違う、または換算額が非有限の場合。
    pub fn new(notional: Money, rate: FxRate<Forward>, position: Position) -> Result<Self> {
        validate_notional(notional)?;
        if notional.ccy() != rate.pair().base() {
            return Err(invalid("FX notional must be in the base currency"));
        }
        if !(notional.amount() * rate.rate()).is_finite() {
            return Err(invalid("FX quote amount must be finite"));
        }
        Ok(Self {
            notional,
            rate,
            position,
        })
    }

    /// 受渡日の二通貨の支払いを組み立てる。
    ///
    /// # Errors
    ///
    /// 通貨の参照データ、現物受渡区分または IR が不正な場合。
    pub fn to_product(&self, refdata: &RefData) -> Result<Product> {
        let (contract, warnings) = self.build_contract(refdata)?;
        Product::with_construction_warnings(contract, refdata, warnings)
    }

    pub(super) fn build_contract(&self, refdata: &RefData) -> Result<(Contract, Vec<Warning>)> {
        self.validate_currencies(refdata)?;
        Ok((
            Contract::build(|builder| Ok(self.exchange(builder)))?,
            Vec::new(),
        ))
    }

    fn validate_currencies(&self, refdata: &RefData) -> Result<()> {
        for ccy in [self.rate.pair().base(), self.rate.pair().quote()] {
            if refdata.currency(ccy)?.deliverability != Deliverability::Deliverable {
                return Err(invalid(format!("{ccy} is not deliverable")));
            }
        }
        Ok(())
    }

    fn exchange<'b>(&self, builder: &mut Builder<'b>) -> ContractIx<'b> {
        let pair = self.rate.pair();
        let at = self.rate.settlement();
        let base = builder.constant(self.notional.amount());
        let quote = builder.constant(-self.notional.amount() * self.rate.rate());
        let base = builder.pay(Asset::Cash(pair.base()), base, at);
        let quote = builder.pay(Asset::Cash(pair.quote()), quote, at);
        let contract = builder.all([base, quote]);
        match self.position {
            Position::Long => contract,
            Position::Short => builder.give(contract),
        }
    }
}

/// 近い日と遠い日で反対方向に二通貨を交換する FX スワップ。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(
    try_from = "super::wire::FxSwapInput",
    into = "super::wire::FxSwapInput"
)]
pub struct FxSwap {
    pub(super) near: FxForward,
    pub(super) far: FxForward,
}

impl FxSwap {
    /// 元本と向きは近い受渡について指定する。
    ///
    /// # Errors
    ///
    /// レートのペアが異なる、受渡が昇順でない、または元本や換算額が不正な場合。
    pub fn new(
        notional: Money,
        near: FxRate<Forward>,
        far: FxRate<Forward>,
        position: Position,
    ) -> Result<Self> {
        if near.pair() != far.pair() || near.settlement() >= far.settlement() {
            return Err(invalid(
                "FX swap requires a single pair and increasing settlement instants",
            ));
        }
        let other_side = match position {
            Position::Long => Position::Short,
            Position::Short => Position::Long,
        };
        Ok(Self {
            near: FxForward::new(notional, near, position)?,
            far: FxForward::new(notional, far, other_side)?,
        })
    }

    /// 近い受渡と遠い受渡をまとめて組み立てる。
    ///
    /// # Errors
    ///
    /// 通貨の参照データ、現物受渡区分または IR が不正な場合。
    pub fn to_product(&self, refdata: &RefData) -> Result<Product> {
        let (contract, warnings) = self.build_contract(refdata)?;
        Product::with_construction_warnings(contract, refdata, warnings)
    }

    pub(super) fn build_contract(&self, refdata: &RefData) -> Result<(Contract, Vec<Warning>)> {
        self.near.validate_currencies(refdata)?;
        let contract = Contract::build(|builder| {
            let near = self.near.exchange(builder);
            let far = self.far.exchange(builder);
            Ok(builder.all([near, far]))
        })?;
        Ok((contract, Vec::new()))
    }
}
`````

`````
crates/rhodium-products/src/standard/leg.rs
//! 固定元本の固定脚と RFR 複利脚を組み立てる。

use rhodium_core::{
    Error, Result, Warning,
    calendar::{BusinessDayConvention, Calendar, CalendarSpec},
    currency::{Asset, Currency},
    date::Date,
    day_count::{DayCount, year_fraction},
    fixing::ObservationPoint,
    id::IndexId,
    money::Money,
    observation::{
        Aggregation, ObservationRule, SpreadTreatment, observation_schedule_with_calendars,
    },
    quantity::{Rate, Simple, Spread, accrue},
    refdata::{CalendarCache, IndexKind, RefData},
    schedule::{
        PeriodAnchor, PeriodCombination, ResetRule, Schedule, ScheduleSpec, SettlementContext,
        generate_with_calendars, reset_dates,
    },
    tenor::{DateOffset, TenorUnit},
};

use crate::ir::{Builder, ContractIx, ExprIx};

use super::{PaymentRounding, invalid, settlement, validate_notional};

/// 脚の受け払い。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PayReceive {
    /// クーポンを支払う。
    Pay,
    /// クーポンを受け取る。
    Receive,
}

impl PayReceive {
    pub(super) const fn opposite(self) -> Self {
        match self {
            Self::Pay => Self::Receive,
            Self::Receive => Self::Pay,
        }
    }
}

/// 脚の元本の受渡。
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NotionalExchange {
    /// クーポンだけを受け渡す。
    #[default]
    None,
    /// 満期に元本を受け取る（支払脚では支払う）。
    Final,
    /// 開始時と満期に反対方向に元本を受け渡す。
    InitialAndFinal,
}

/// 脚の元本。
#[derive(Debug, Clone, PartialEq)]
pub enum Notional {
    /// 全期間で同じ元本。
    Fixed(Money),
    /// 基準元本を各期の FX fixing で相手通貨に換算する。
    FxReset {
        /// 換算元の正の元本。
        base: Money,
        /// 基準通貨と支払通貨のペアを指定する指標。
        fixing: IndexId,
    },
}

impl Notional {
    pub(super) fn currency(&self, refdata: &RefData) -> Result<Currency> {
        match self {
            Self::Fixed(amount) => Ok(amount.ccy()),
            Self::FxReset { base, fixing } => {
                let IndexKind::FxFixing { pair } = refdata.index(fixing)?.kind else {
                    return Err(invalid(format!("{fixing} is not an FX fixing")));
                };
                if base.ccy() == pair.base() {
                    Ok(pair.quote())
                } else if base.ccy() == pair.quote() {
                    Ok(pair.base())
                } else {
                    Err(invalid(format!(
                        "FX fixing {fixing} does not contain the base notional currency"
                    )))
                }
            }
        }
    }
}

/// 脚のクーポン条件。
#[derive(Debug, Clone, PartialEq)]
pub enum Coupon {
    /// 単利の固定クーポン。
    Fixed {
        /// 日数基準付きの固定の年率。
        rate: Rate<Simple>,
    },
    /// リセット日に読むターム金利のクーポン。
    Term {
        /// ターム金利の指標。
        index: IndexId,
        /// クーポンのスプレッド。
        spread: Spread,
        /// 指標の率に掛ける係数。
        gearing: f64,
    },
    /// 翌日物金利の複利クーポン。
    Compounded {
        /// 翌日物金利の指標。
        index: IndexId,
        /// 観測・集約・丸めの規約。
        rule: ObservationRule,
        /// クーポンに適用するスプレッド。
        spread: Spread,
    },
}

/// 固定元本または FX リセット元本でクーポンを受け渡す脚。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(try_from = "super::wire::LegInput", into = "super::wire::LegInput")]
pub struct Leg {
    pub(super) side: PayReceive,
    pub(super) notional: Notional,
    pub(super) schedule: ScheduleSpec,
    pub(super) day_count: DayCount,
    pub(super) coupon: Coupon,
    pub(super) exchange: NotionalExchange,
    pub(super) payment_rounding: PaymentRounding,
}

impl Leg {
    /// 確認書の条件から脚を組み立てる。
    ///
    /// # Errors
    ///
    /// 元本が正の有限値でない、率と脚の日数基準が異なる、RFR の規約が
    /// 複利用でない、または明示日程に 30/360 US を指定した場合は
    /// [`Error::InvalidData`]。
    pub fn new(
        side: PayReceive,
        notional: Money,
        schedule: ScheduleSpec,
        day_count: DayCount,
        coupon: Coupon,
        exchange: NotionalExchange,
        payment_rounding: PaymentRounding,
    ) -> Result<Self> {
        Self::with_notional(
            side,
            Notional::Fixed(notional),
            schedule,
            day_count,
            coupon,
            exchange,
            payment_rounding,
        )
    }

    /// 固定または FX リセットの元本を指定する。
    ///
    /// # Errors
    ///
    /// 元本、脚の規約または日程が不正な場合。FX リセットは開始時・途中・満期に元本を受け渡す。
    pub fn with_notional(
        side: PayReceive,
        notional: Notional,
        schedule: ScheduleSpec,
        day_count: DayCount,
        coupon: Coupon,
        exchange: NotionalExchange,
        payment_rounding: PaymentRounding,
    ) -> Result<Self> {
        validate_notional(match &notional {
            Notional::Fixed(amount) | Notional::FxReset { base: amount, .. } => *amount,
        })?;
        if matches!(notional, Notional::FxReset { .. })
            && exchange != NotionalExchange::InitialAndFinal
        {
            return Err(invalid(
                "FX reset requires initial and final notional exchanges",
            ));
        }
        match &coupon {
            Coupon::Fixed { rate } if rate.basis() != &day_count => {
                return Err(invalid("fixed rate and leg day counts differ"));
            }
            Coupon::Term { gearing, .. } if !gearing.is_finite() => {
                return Err(invalid("term coupon gearing must be finite"));
            }
            Coupon::Compounded { rule, .. }
                if !matches!(rule.aggregation(), Aggregation::Compound { .. }) =>
            {
                return Err(invalid(
                    "RFR compounding requires a compounding observation rule",
                ));
            }
            Coupon::Fixed { .. } | Coupon::Compounded { .. } | Coupon::Term { .. } => {}
        }
        validate_explicit_basis(&schedule, &day_count)?;
        Ok(Self {
            side,
            notional,
            schedule,
            day_count,
            coupon,
            exchange,
            payment_rounding,
        })
    }

    pub(super) fn build<'b>(
        &self,
        builder: &mut Builder<'b>,
        trade_date: Date,
        refdata: &RefData,
        calendars: &mut CalendarCache<'_>,
        warnings: &mut Vec<Warning>,
    ) -> Result<ContractIx<'b>> {
        let currency = refdata.currency(self.notional.currency(refdata)?)?;
        let asset = Asset::Cash(currency.code);
        let settlement_calendar = CalendarSpec::Base(currency.settlement_calendar.clone());
        let (schedule, notices) = generate_with_calendars(
            &self.schedule,
            trade_date,
            &SettlementContext {
                zones: refdata.zones(),
                default_time: Some(&currency.settlement_time),
                calendar: Some(&settlement_calendar),
            },
            calendars,
        )?;
        warnings.extend(notices);
        if schedule.combine() != PeriodCombination::Sum {
            return Err(Error::Unsupported(
                "leg payments with compounding across accrual periods".into(),
            ));
        }
        schedule.validate_day_count(&self.day_count)?;
        let coupon = self.resolve_coupon(refdata, &schedule, currency.code, warnings)?;
        let calendar = day_count_calendar(&self.day_count, refdata)?;
        let notionals = self.notionals(builder, &schedule, currency, refdata, warnings)?;
        let amounts = self.coupons(
            builder,
            &schedule,
            &coupon,
            calendar.as_ref(),
            calendars,
            warnings,
            &notionals,
        )?;
        let mut payments = Vec::with_capacity(schedule.payments().len() + 2);

        // 期ごとの額は先に合算し、確認書が指定した場合だけ受渡額を丸める。
        for payment in schedule.payments() {
            let amounts = &amounts[payment.periods()];
            let total = match amounts {
                [amount] => *amount,
                _ => builder.add(amounts.iter().copied()),
            };
            let total = self.payment_rounding.apply(builder, total, currency);
            payments.push(builder.pay(asset.clone(), total, payment.at()));
        }
        match self.exchange {
            NotionalExchange::None => {}
            NotionalExchange::Final | NotionalExchange::InitialAndFinal => {
                // 計算期間の始点・終点に交換する。クーポンの支払ラグは使わない。
                let periods = schedule.periods();
                if self.exchange == NotionalExchange::InitialAndFinal {
                    let start = periods[0].adjusted().start();
                    let initial = self.payment_rounding.apply(builder, notionals[0], currency);
                    let outgoing = builder.neg(initial);
                    payments.push(builder.pay(
                        asset.clone(),
                        outgoing,
                        settlement(start, currency, refdata)?,
                    ));
                }
                if matches!(self.notional, Notional::FxReset { .. }) {
                    // リセット済みの元本と新しい元本との差額だけを交換する。
                    for (index, period) in periods.iter().enumerate().skip(1) {
                        let difference = builder.sub(notionals[index - 1], notionals[index]);
                        let difference = self.payment_rounding.apply(builder, difference, currency);
                        let at = settlement(period.adjusted().start(), currency, refdata)?;
                        payments.push(builder.pay(asset.clone(), difference, at));
                    }
                }
                let end = periods[periods.len() - 1].adjusted().end();
                let notional =
                    self.payment_rounding
                        .apply(builder, notionals[periods.len() - 1], currency);
                payments.push(builder.pay(asset, notional, settlement(end, currency, refdata)?));
            }
        }
        let contract = builder.all(payments);
        Ok(match self.side {
            PayReceive::Pay => builder.give(contract),
            PayReceive::Receive => contract,
        })
    }

    fn notionals<'b>(
        &self,
        builder: &mut Builder<'b>,
        schedule: &Schedule,
        currency: &rhodium_core::refdata::CurrencyDefinition,
        refdata: &RefData,
        warnings: &mut Vec<Warning>,
    ) -> Result<Vec<ExprIx<'b>>> {
        match &self.notional {
            Notional::Fixed(amount) => Ok(vec![
                builder.constant(amount.amount());
                schedule.periods().len()
            ]),
            Notional::FxReset { base, fixing } => {
                let definition = refdata.index(fixing)?;
                let IndexKind::FxFixing { pair } = definition.kind else {
                    return Err(invalid(format!("{fixing} is not an FX fixing")));
                };
                // 期首がスポット受渡になる日（ペアの慣行のスポットラグを両通貨の休日の和で遡った日）に
                // FX を観測し、換算した元本を期首に受け渡す。明示日程の `reset` はこの規則より優先する。
                let spot_lag = refdata.fx_convention(pair)?.spot_lag;
                let (dates, notices) = reset_dates(
                    schedule,
                    &ResetRule {
                        anchor: PeriodAnchor::PeriodStart,
                        offset: DateOffset::new(-spot_lag.business_days(), TenorUnit::BusinessDay),
                        calendar: refdata.fx_calendar_spec(pair)?,
                        bdc: BusinessDayConvention::Preceding,
                    },
                    refdata.calendars(),
                )?;
                warnings.extend(notices);
                let amount = builder.constant(base.amount());
                dates
                    .into_iter()
                    .zip(schedule.periods())
                    .map(|(date, period)| {
                        let at = settlement(period.adjusted().start(), currency, refdata)?;
                        let observed = builder.fx_reset(fixing, ObservationPoint::Date(date), at);
                        Ok(if base.ccy() == pair.base() {
                            builder.mul([amount, observed])
                        } else {
                            builder.div(amount, observed)
                        })
                    })
                    .collect()
            }
        }
    }

    fn resolve_coupon<'a>(
        &'a self,
        refdata: &'a RefData,
        schedule: &Schedule,
        currency: Currency,
        warnings: &mut Vec<Warning>,
    ) -> Result<ResolvedCoupon<'a>> {
        Ok(match &self.coupon {
            Coupon::Fixed { rate } => ResolvedCoupon::Fixed(rate),
            Coupon::Term {
                index,
                spread,
                gearing,
            } => {
                let definition = refdata.index(index)?;
                let IndexKind::TermRate {
                    ccy, value_offset, ..
                } = definition.kind
                else {
                    return Err(invalid(format!("index {index} is not a term rate")));
                };
                if ccy != currency {
                    return Err(invalid(format!(
                        "index {index} currency differs from leg currency"
                    )));
                }
                let offset = value_offset
                    .amount()
                    .checked_neg()
                    .ok_or_else(|| invalid("term index value offset cannot be reversed"))?;
                let (dates, notices) = reset_dates(
                    schedule,
                    &ResetRule {
                        anchor: PeriodAnchor::PeriodStart,
                        offset: DateOffset::new(offset, value_offset.unit()),
                        calendar: definition.calendar.clone(),
                        bdc: BusinessDayConvention::Preceding,
                    },
                    refdata.calendars(),
                )?;
                warnings.extend(notices);
                ResolvedCoupon::Term {
                    index,
                    dates,
                    spread: *spread,
                    gearing: *gearing,
                }
            }
            Coupon::Compounded {
                index,
                spread,
                rule,
            } => {
                let IndexKind::OvernightRate { ccy, day_count } = &refdata.index(index)?.kind
                else {
                    return Err(invalid(format!("index {index} is not an overnight rate")));
                };
                if *ccy != currency {
                    return Err(invalid(format!(
                        "index {index} currency differs from leg currency"
                    )));
                }
                validate_explicit_basis(&self.schedule, day_count)?;
                if day_count != &self.day_count {
                    schedule.validate_day_count(day_count)?;
                }
                ResolvedCoupon::Compounded {
                    index,
                    spread: *spread,
                    rule,
                    basis: day_count,
                }
            }
        })
    }

    #[expect(
        clippy::too_many_arguments,
        reason = "日程と観測を解決済みの条件で評価する"
    )]
    fn coupons<'b>(
        &self,
        builder: &mut Builder<'b>,
        schedule: &Schedule,
        coupon: &ResolvedCoupon<'_>,
        calendar: Option<&Calendar>,
        calendars: &mut CalendarCache<'_>,
        warnings: &mut Vec<Warning>,
        notionals: &[ExprIx<'b>],
    ) -> Result<Vec<ExprIx<'b>>> {
        schedule
            .periods()
            .iter()
            .enumerate()
            .map(|(i, period)| {
                let context = schedule
                    .accrual_context(i, calendar)
                    .ok_or_else(|| Error::Internal("missing schedule accrual context".into()))?;
                let interval = period.accrual();
                let (fraction, notices) =
                    year_fraction(&self.day_count, interval.start(), interval.end(), &context)?;
                warnings.extend(notices);
                let interest = match coupon {
                    ResolvedCoupon::Fixed(rate) => builder.constant(accrue(rate, &fraction)?),
                    ResolvedCoupon::Term {
                        index,
                        dates,
                        gearing,
                        spread,
                    } => {
                        let mut rate = builder.fixing(index, ObservationPoint::Date(dates[i]));
                        if gearing.to_bits() != 1.0_f64.to_bits() {
                            let coefficient = builder.constant(*gearing);
                            rate = builder.mul([coefficient, rate]);
                        }
                        if spread.value() != 0.0 {
                            let spread = builder.constant(spread.value());
                            rate = builder.add([rate, spread]);
                        }
                        let coupon_fraction = builder.constant(fraction.value());
                        builder.mul([rate, coupon_fraction])
                    }
                    ResolvedCoupon::Compounded {
                        index,
                        spread,
                        rule,
                        basis,
                    } => {
                        let (observations, notices) = observation_schedule_with_calendars(
                            period, rule, basis, &context, calendars,
                        )?;
                        warnings.extend(notices);
                        let annualization = observations.annualization().value();
                        if annualization <= 0.0 {
                            return Err(invalid(
                                "RFR coupon requires a positive annualization fraction",
                            ));
                        }
                        let growth = builder.compounded(index, &observations, rule, *spread)?;
                        let one = builder.constant(1.0);
                        let excess = builder.sub(growth, one);
                        let denominator = builder.constant(annualization);
                        let mut rate = builder.div(excess, denominator);
                        if matches!(
                            rule.aggregation(),
                            Aggregation::Compound {
                                spread: SpreadTreatment::Exclusive,
                                ..
                            }
                        ) && spread.value() != 0.0
                        {
                            let margin = builder.constant(spread.value());
                            rate = builder.add([rate, margin]);
                        }
                        if let Some(rounding) = rule.rounding().aggregated {
                            rate = builder.round(rate, rounding);
                        }
                        let coupon_fraction = builder.constant(fraction.value());
                        builder.mul([rate, coupon_fraction])
                    }
                };
                Ok(builder.mul([notionals[i], interest]))
            })
            .collect()
    }
}

enum ResolvedCoupon<'a> {
    Fixed(&'a Rate<Simple>),
    Term {
        index: &'a IndexId,
        dates: Vec<Date>,
        spread: Spread,
        gearing: f64,
    },
    Compounded {
        index: &'a IndexId,
        spread: Spread,
        rule: &'a ObservationRule,
        basis: &'a DayCount,
    },
}

fn validate_explicit_basis(schedule: &ScheduleSpec, basis: &DayCount) -> Result<()> {
    if matches!(
        (schedule, basis),
        (ScheduleSpec::Explicit(_), DayCount::Thirty360Us)
    ) {
        return Err(invalid(
            "30/360 US requires a schedule with a known month-end convention",
        ));
    }
    Ok(())
}

pub(super) fn day_count_calendar(basis: &DayCount, refdata: &RefData) -> Result<Option<Calendar>> {
    match basis {
        DayCount::Bus252 { calendar } => refdata
            .calendar(&CalendarSpec::Base(calendar.clone()))
            .map(Some),
        _ => Ok(None),
    }
}
`````

`````
crates/rhodium-products/src/standard/rates.rs
use rhodium_core::{
    Error, Result, Warning,
    calendar::{AdvanceRule, CalendarSpec, EndOfMonthRule, money_market_period},
    currency::Asset,
    date::{Date, DateInterval},
    day_count::{AccrualContext, DayCount, MonthEndSchedule, year_fraction},
    id::ConventionId,
    money::Money,
    quantity::{Rate, Simple, Spread, accrue},
    refdata::{Convention, IndexKind, RefData},
    schedule::{
        AccrualDates, EndRule, GeneratedPeriods, ScheduleSpec, StartRule, StubLength, StubSpec,
    },
    tenor::{Tenor, TenorUnit},
};

use crate::ir::{Contract, Product};

use super::{
    PaymentRounding, Position, invalid,
    leg::{Coupon, Leg, Notional, NotionalExchange, PayReceive, day_count_calendar},
    settlement, validate_notional,
};

/// 預金の期間の指定。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DepositDates {
    /// 確認済みの開始日と満期日。
    Confirmed {
        /// 元本を受け渡す日。
        start: Date,
        /// 元利金を受け渡す日。
        maturity: Date,
    },
    /// 約定日と預金の慣行から求めるテナー。
    Market(Tenor),
}

/// 単利の預金。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(
    try_from = "super::wire::DepositInput",
    into = "super::wire::DepositInput"
)]
pub struct Deposit {
    pub(super) notional: Money,
    pub(super) rate: f64,
    pub(super) dates: DepositDates,
    pub(super) convention: Option<ConventionId>,
    pub(super) day_count: Option<DayCount>,
    pub(super) position: Position,
    pub(super) payment_rounding: PaymentRounding,
}

impl Deposit {
    /// 年率を小数で受け取り、確認書の条件を組み立てる。
    ///
    /// 日数基準の指定があれば慣行より優先する。テナー指定には慣行が必要。
    ///
    /// # Errors
    ///
    /// 元本・年率・期間が不正か、必要な慣行または日数基準が指定されていない場合。
    pub fn new(
        notional: Money,
        rate: f64,
        dates: DepositDates,
        convention: Option<ConventionId>,
        day_count: Option<DayCount>,
        position: Position,
        payment_rounding: PaymentRounding,
    ) -> Result<Self> {
        validate_notional(notional)?;
        if !rate.is_finite() {
            return Err(invalid("deposit rate must be finite"));
        }
        match dates {
            DepositDates::Confirmed { start, maturity } if start >= maturity => {
                return Err(invalid("deposit start must precede maturity"));
            }
            DepositDates::Market(_) if convention.is_none() => {
                return Err(invalid("a deposit tenor requires a convention"));
            }
            DepositDates::Market(Tenor::Length(length)) if length.amount() == 0 => {
                return Err(invalid("deposit tenor must be positive"));
            }
            DepositDates::Confirmed { .. } | DepositDates::Market(_) => {}
        }
        if convention.is_none() && day_count.is_none() {
            return Err(invalid("deposit requires a day count or convention"));
        }
        if let Some(basis) = &day_count {
            validate_deposit_basis(dates, basis)?;
        }
        Ok(Self {
            notional,
            rate,
            dates,
            convention,
            day_count,
            position,
            payment_rounding,
        })
    }

    /// 約定日と参照データから元本と元利金の受渡を組み立てる。
    ///
    /// # Errors
    ///
    /// 慣行、日数基準、日程や決済を解決できないか、金額が値域を超えた場合。
    pub fn to_product(&self, trade_date: Date, refdata: &RefData) -> Result<Product> {
        self.build_contract(trade_date, refdata)
            .and_then(|(contract, warnings)| {
                Product::with_construction_warnings(contract, refdata, warnings)
            })
            .map_err(|error| error.context("deposit"))
    }

    pub(super) fn build_contract(
        &self,
        trade_date: Date,
        refdata: &RefData,
    ) -> Result<(Contract, Vec<Warning>)> {
        let convention = self
            .convention
            .as_ref()
            .map(|id| match refdata.convention(id)? {
                Convention::Deposit(convention) => Ok(convention),
                _ => Err(invalid(format!("{id} is not a deposit convention"))),
            })
            .transpose()?;
        let basis = if let Some(basis) = &self.day_count {
            basis
        } else {
            let basis = &convention
                .ok_or_else(|| Error::Internal("missing deposit day count".into()))?
                .day_count;
            validate_deposit_basis(self.dates, basis)?;
            basis
        };

        let currency = refdata.currency(self.notional.ccy())?;
        let mut context = AccrualContext::default();
        let (period, mut warnings) = match self.dates {
            DepositDates::Confirmed { start, maturity } => {
                let calendar =
                    refdata.calendar(&CalendarSpec::Base(currency.settlement_calendar.clone()))?;
                let mut warnings = Vec::new();
                for date in [start, maturity] {
                    match calendar.is_business_day(date) {
                        Ok((open, notices)) => {
                            warnings.extend(notices);
                            if !open {
                                warnings.push(Warning::ExplicitDateIsHoliday { date });
                            }
                        }
                        // 確認書の受渡日は休日表の収録範囲が変わっても確定している。
                        Err(Error::MissingData(_)) => {}
                        Err(error) => return Err(error),
                    }
                }
                (DateInterval::new(start, maturity)?, warnings)
            }
            DepositDates::Market(tenor) => {
                let convention = convention
                    .ok_or_else(|| Error::Internal("missing deposit convention".into()))?;
                let calendar = refdata.calendar(&convention.calendar)?;
                if convention.end_of_month != EndOfMonthRule::Ignore
                    && matches!(tenor, Tenor::Length(length) if matches!(length.unit(), TenorUnit::Month | TenorUnit::Year))
                {
                    context.month_end = MonthEndSchedule::Aligned;
                }
                money_market_period(
                    trade_date,
                    tenor,
                    convention.spot_lag,
                    &AdvanceRule {
                        calendar: &calendar,
                        bdc: convention.bdc,
                        end_of_month: convention.end_of_month,
                    },
                )?
            }
        };
        if period.is_empty() {
            return Err(invalid("deposit accrual period is empty after adjustment"));
        }
        let calendar = day_count_calendar(basis, refdata)?;
        context.calendar = calendar.as_ref();
        context.maturity = Some(period.end());
        let (fraction, notices) = year_fraction(basis, period.start(), period.end(), &context)?;
        warnings.extend(notices);
        let rate = Rate::new(self.rate, basis.clone(), Simple)?;
        let principal = self.notional.amount();
        let redemption = principal * (1.0 + accrue(&rate, &fraction)?);
        if !redemption.is_finite() {
            return Err(Error::Numerical(
                "deposit redemption is outside the finite range".into(),
            ));
        }
        let start = settlement(period.start(), currency, refdata)?;
        let end = settlement(period.end(), currency, refdata)?;
        let asset = Asset::Cash(currency.code);
        let contract = Contract::build(|builder| {
            let initial = builder.constant(principal);
            let initial = self.payment_rounding.apply(builder, initial, currency);
            let initial = builder.neg(initial);
            let final_amount = builder.constant(redemption);
            let final_amount = self.payment_rounding.apply(builder, final_amount, currency);
            let start = builder.pay(asset.clone(), initial, start);
            let end = builder.pay(asset, final_amount, end);
            let contract = builder.all([start, end]);
            Ok(match self.position {
                Position::Long => contract,
                Position::Short => builder.give(contract),
            })
        })?;
        Ok((contract, warnings))
    }
}

fn validate_deposit_basis(dates: DepositDates, basis: &DayCount) -> Result<()> {
    match (dates, basis) {
        (_, DayCount::ActActIcma | DayCount::Act365L) => Err(invalid(
            "deposit has no regular reference period or frequency",
        )),
        (DepositDates::Confirmed { .. }, DayCount::Thirty360Us) => {
            Err(invalid("30/360 US requires a known month-end convention"))
        }
        _ => Ok(()),
    }
}

/// 固定・変動の脚を組み合わせる金利スワップまたは通貨スワップ。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(try_from = "super::wire::SwapInput", into = "super::wire::SwapInput")]
pub struct Swap(pub(super) SwapTerms);

#[derive(Debug, Clone, PartialEq)]
pub(super) enum SwapTerms {
    Convention(Ois),
    Legs(Box<[Leg; 2]>),
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct Ois {
    pub(super) convention: ConventionId,
    pub(super) start: Option<Date>,
    pub(super) end: EndRule,
    pub(super) notional: Money,
    pub(super) fixed_rate: f64,
    pub(super) fixed_side: PayReceive,
    pub(super) spread: Spread,
    pub(super) payment_rounding: PaymentRounding,
}

impl Swap {
    /// OIS 慣行から二つの脚を導く。開始日を省略すると約定日からスポット日を求める。
    /// 指定した開始日と満期日は、営業日調整前の日程境界として扱う。
    /// 日程は満期から遡って作り、端数期間は先頭の短い期間にする。
    ///
    /// # Errors
    ///
    /// 元本・固定金利・期間の指定が不正な場合。
    #[expect(
        clippy::too_many_arguments,
        reason = "確認書が独立に指定する OIS の条件"
    )]
    pub fn from_convention(
        convention: ConventionId,
        start: Option<Date>,
        end: EndRule,
        notional: Money,
        fixed_rate: f64,
        fixed_side: PayReceive,
        spread: Spread,
        payment_rounding: PaymentRounding,
    ) -> Result<Self> {
        validate_notional(notional)?;
        if !fixed_rate.is_finite() {
            return Err(invalid("swap fixed rate must be finite"));
        }
        match (start, end) {
            (Some(start), EndRule::Date(end)) if start >= end => {
                return Err(invalid("swap start must precede maturity"));
            }
            (_, EndRule::AfterStart(length)) if length.amount() == 0 => {
                return Err(invalid("swap tenor must be positive"));
            }
            _ => {}
        }
        Ok(Self(SwapTerms::Convention(Ois {
            convention,
            start,
            end,
            notional,
            fixed_rate,
            fixed_side,
            spread,
            payment_rounding,
        })))
    }

    /// 確認書で条件を指定した二つの脚を組にする。
    ///
    /// # Errors
    ///
    /// 二つの脚の向きが同じか、脚または元本交換の組み合わせが不正な場合。
    pub fn from_legs(legs: [Leg; 2]) -> Result<Self> {
        let [first, second] = &legs;
        if first.side == second.side {
            return Err(invalid("swap legs must have opposite sides"));
        }
        match (&first.notional, &second.notional) {
            (Notional::Fixed(a), Notional::Fixed(b)) if a.ccy() == b.ccy() => {
                if matches!(
                    (&first.coupon, &second.coupon),
                    (Coupon::Fixed { .. }, Coupon::Fixed { .. })
                ) {
                    return Err(invalid("single-currency swap requires a floating leg"));
                }
            }
            (Notional::Fixed(_), Notional::Fixed(_)) => {
                validate_cross_currency_exchange(first, second)?;
            }
            (Notional::Fixed(fixed), Notional::FxReset { base, .. })
            | (Notional::FxReset { base, .. }, Notional::Fixed(fixed)) => {
                if fixed.ccy() != base.ccy() {
                    return Err(invalid(
                        "FX reset reference currency differs from the fixed notional currency",
                    ));
                }
                validate_cross_currency_exchange(first, second)?;
            }
            (Notional::FxReset { .. }, Notional::FxReset { .. }) => {
                return Err(invalid("MtM swap requires one fixed notional leg"));
            }
        }
        Ok(Self(SwapTerms::Legs(Box::new(legs))))
    }

    /// 各脚の日程を解決して検証済み IR にする。
    ///
    /// # Errors
    ///
    /// 慣行、指標、日程や受渡の条件が不正か、計算に失敗した場合。
    pub fn to_product(&self, trade_date: Date, refdata: &RefData) -> Result<Product> {
        self.build_contract(trade_date, refdata)
            .and_then(|(contract, warnings)| {
                Product::with_construction_warnings(contract, refdata, warnings)
            })
            .map_err(|error| error.context("swap"))
    }

    pub(super) fn build_contract(
        &self,
        trade_date: Date,
        refdata: &RefData,
    ) -> Result<(Contract, Vec<Warning>)> {
        let resolved;
        let legs = match &self.0 {
            SwapTerms::Legs(legs) => legs.as_ref(),
            SwapTerms::Convention(ois) => {
                resolved = ois.legs(refdata)?;
                &resolved
            }
        };
        // 慣行外の暦もここで一度だけ解決し、両脚と全観測期間で共有する。
        let mut calendars = refdata.calendar_cache();
        let mut warnings = Vec::new();
        let contract = Contract::build(|builder| {
            let first = legs[0]
                .build(builder, trade_date, refdata, &mut calendars, &mut warnings)
                .map_err(|error| error.context("first leg"))?;
            let second = legs[1]
                .build(builder, trade_date, refdata, &mut calendars, &mut warnings)
                .map_err(|error| error.context("second leg"))?;
            Ok(builder.all([first, second]))
        })?;
        Ok((contract, warnings))
    }
}

fn validate_cross_currency_exchange(first: &Leg, second: &Leg) -> Result<()> {
    if [first, second]
        .iter()
        .all(|leg| leg.exchange == NotionalExchange::InitialAndFinal)
    {
        Ok(())
    } else {
        Err(invalid(
            "cross-currency swap requires initial and final exchanges on both legs",
        ))
    }
}

impl Ois {
    fn legs(&self, refdata: &RefData) -> Result<[Leg; 2]> {
        let Convention::Ois(ois) = refdata.convention(&self.convention)? else {
            return Err(invalid(format!(
                "{} is not an OIS convention",
                self.convention
            )));
        };
        let IndexKind::OvernightRate { day_count, .. } = &refdata.index(&ois.float.index)?.kind
        else {
            return Err(Error::Internal(
                "OIS convention has no overnight index".into(),
            ));
        };
        let start = self
            .start
            .map_or(StartRule::FromTradeDate(ois.spot_lag), StartRule::Date);
        // 標準の OIS は満期から遡って日程を作り、端数期間を先頭の短い期間にする。
        // QuantLib の MakeOIS（DateGeneration::Backward）とも同じで、18M などの
        // 年の整数倍でないテナーでも、quote と同じ期間の分け方になる。
        let schedule = ScheduleSpec::Generated {
            periods: GeneratedPeriods::new(
                start,
                self.end,
                ois.fixed.frequency,
                None,
                StubSpec::Initial(StubLength::Short),
                ois.fixed.calendar.clone(),
                ois.fixed.bdc,
                AccrualDates::Adjusted,
            )?,
            payment: ois.payment.clone(),
        };
        let fixed_basis = &ois.fixed.day_count;
        let fixed = Leg::new(
            self.fixed_side,
            self.notional,
            schedule.clone(),
            fixed_basis.clone(),
            Coupon::Fixed {
                rate: Rate::new(self.fixed_rate, fixed_basis.clone(), Simple)?,
            },
            NotionalExchange::None,
            self.payment_rounding,
        )?;
        let floating = Leg::new(
            self.fixed_side.opposite(),
            self.notional,
            schedule,
            day_count.clone(),
            Coupon::Compounded {
                index: ois.float.index.clone(),
                rule: ois.float.observation.clone(),
                spread: self.spread,
            },
            NotionalExchange::None,
            self.payment_rounding,
        )?;
        Ok([fixed, floating])
    }
}
`````

`````
crates/rhodium-products/src/standard/template.rs
//! 慣行と座標から校正商品の標準商品を作るテンプレート。

use std::num::NonZeroU32;

use rhodium_core::{
    Error, Result,
    calendar::{AdvanceRule, EndOfMonthRule, advance, money_market_period},
    currency::Currency,
    date::Date,
    fixing::YearMonth,
    fx::FxRate,
    id::ConventionId,
    money::Money,
    observation::{Aggregation, SpreadTreatment},
    quantity::{Rate, RateFuturesPrice, Simple, Spread},
    refdata::{Convention, IndexKind, RefData, SwapConvention, SwapCouponConvention},
    schedule::{
        AccrualDates, EndRule, GeneratedPeriods, ScheduleSpec, StartRule, StubLength, StubSpec,
    },
    tenor::{DateOffset, Tenor, TenorLength},
    time::Instant,
};

use super::{
    Coupon, Deposit, DepositDates, Fra, FxForward, Leg, Notional, NotionalExchange, PayReceive,
    PaymentRounding, Position, RateFuture, StandardProduct, Swap, invalid,
};

/// 慣行とテナーまたは限月から校正商品を作る指定。
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum InstrumentTemplate {
    /// スポット起点の預金。
    Deposit {
        /// 商品条件の参照先。
        convention: ConventionId,
        /// 満期を決めるテナー。
        tenor: Tenor,
        /// 取引の通貨。
        currency: Currency,
    },
    /// スポットから期首と期末を数える FRA。
    Fra {
        /// 商品条件の参照先。
        convention: ConventionId,
        /// スポットから期首までの長さ。
        start: TenorLength,
        /// スポットから期末までの長さ。
        end: TenorLength,
    },
    /// 限月の金利先物。
    RateFuture {
        /// 商品条件の参照先。
        convention: ConventionId,
        /// 対象の限月。
        delivery_month: YearMonth,
    },
    /// スポット起点の OIS。
    Ois {
        /// 商品条件の参照先。
        convention: ConventionId,
        /// 満期を決めるテナー。
        tenor: TenorLength,
    },
    /// 固定対ターム金利スワップ。
    Swap {
        /// 商品条件の参照先。
        convention: ConventionId,
        /// 満期を決めるテナー。
        tenor: TenorLength,
    },
    /// 同一通貨のベーシススワップ。
    BasisSwap {
        /// 商品条件の参照先。
        convention: ConventionId,
        /// 満期を決めるテナー。
        tenor: TenorLength,
    },
    /// アウトライトで条件を指定する FX フォワード。
    FxForward {
        /// 商品条件の参照先。
        convention: ConventionId,
        /// 満期を決めるテナー。
        tenor: Tenor,
    },
    /// XCCY の払い脚にベーシスを適用するスワップ。
    XccyBasis {
        /// 商品条件の参照先。
        convention: ConventionId,
        /// 満期を決めるテナー。
        tenor: TenorLength,
    },
}

impl InstrumentTemplate {
    /// 慣行の識別子。
    #[must_use]
    pub const fn convention(&self) -> &ConventionId {
        match self {
            Self::Deposit { convention, .. }
            | Self::Fra { convention, .. }
            | Self::RateFuture { convention, .. }
            | Self::Ois { convention, .. }
            | Self::Swap { convention, .. }
            | Self::BasisSwap { convention, .. }
            | Self::FxForward { convention, .. }
            | Self::XccyBasis { convention, .. } => convention,
        }
    }

    /// 約定日を使わないテンプレートから、受渡額の丸めなしの標準商品を作る。
    ///
    /// # Errors
    ///
    /// 慣行や商品の条件が不正か、FRA・FX の受渡を求めるための約定日が必要な場合。
    pub fn build(&self, value: f64, refdata: &RefData) -> Result<StandardProduct> {
        self.build_with_dates(value, None, None, refdata)
    }

    /// 約定日を指定して、校正値を組み入れた確認書の条件を作る。
    ///
    /// # Errors
    ///
    /// 慣行や商品の条件が不正か、日程の解決に失敗した場合。
    pub fn build_for_trade(
        &self,
        value: f64,
        trade_date: Date,
        refdata: &RefData,
    ) -> Result<StandardProduct> {
        self.build_with_dates(value, Some(trade_date), None, refdata)
    }

    /// 校正時点に新たに建てる標準商品を作る。
    ///
    /// 先物の最初の清算は、評価時点より後に価格を決める最初の回とする。
    ///
    /// # Errors
    ///
    /// 慣行や商品の条件、清算日程を解決できない場合。
    pub fn build_for_valuation(
        &self,
        value: f64,
        trade_date: Date,
        as_of: Instant,
        refdata: &RefData,
    ) -> Result<StandardProduct> {
        self.build_with_dates(value, Some(trade_date), Some(as_of), refdata)
    }

    #[expect(
        clippy::too_many_lines,
        reason = "各テンプレートから標準商品への対応を一か所に示す"
    )]
    fn build_with_dates(
        &self,
        value: f64,
        trade_date: Option<Date>,
        as_of: Option<Instant>,
        refdata: &RefData,
    ) -> Result<StandardProduct> {
        let ccy = self.currency(refdata)?;
        let unit = Money::new(ccy, 1.0);
        let id = self.convention();
        let convention = refdata.convention(id)?;
        let trade_date = || trade_date.ok_or_else(|| invalid("template requires a trade date"));
        match (self, convention) {
            (Self::Deposit { tenor, .. }, Convention::Deposit(_)) => {
                Ok(StandardProduct::Deposit(Deposit::new(
                    unit,
                    value,
                    DepositDates::Market(*tenor),
                    Some(id.clone()),
                    None,
                    Position::Long,
                    PaymentRounding::None,
                )?))
            }
            (Self::Ois { tenor, .. }, Convention::Ois(_)) => {
                Ok(StandardProduct::Swap(Swap::from_convention(
                    id.clone(),
                    None,
                    EndRule::AfterStart(*tenor),
                    unit,
                    value,
                    PayReceive::Pay,
                    Spread::new(0.0)?,
                    PaymentRounding::None,
                )?))
            }
            (Self::Fra { start, end, .. }, Convention::Fra(fra)) => {
                let definition = refdata.index(&fra.index)?;
                let IndexKind::TermRate {
                    value_offset,
                    bdc,
                    end_of_month,
                    ..
                } = definition.kind
                else {
                    return Err(Error::Internal("FRA has no term index".into()));
                };
                let calendar = refdata.calendar(&definition.calendar)?;
                let roll = AdvanceRule {
                    bdc,
                    end_of_month,
                    calendar: &calendar,
                };
                let spot = advance(trade_date()?, fra.spot_lag, &roll)?.0;
                let first = advance(spot, tenor_offset(*start)?, &roll)?.0;
                let last = advance(spot, tenor_offset(*end)?, &roll)?.0;
                let reverse = value_offset
                    .amount()
                    .checked_neg()
                    .ok_or_else(|| invalid("FRA fixing offset cannot be reversed"))?;
                let fixing =
                    advance(first, DateOffset::new(reverse, value_offset.unit()), &roll)?.0;
                let (period, _) = refdata.term_rate_dates(&fra.index, fixing)?;
                if first >= last || period.start() != first || period.end() != last {
                    return Err(invalid(
                        "FRA start and end do not match the fixing's term period",
                    ));
                }
                Ok(StandardProduct::Fra(Fra::new(
                    unit,
                    fra.index.clone(),
                    fixing,
                    value,
                    Position::Long,
                )?))
            }
            (Self::RateFuture { delivery_month, .. }, Convention::RateFuture(future)) => {
                let index = &future.contracts[delivery_month].settlement_index;
                let (_, _, basis) = refdata.rate_future_underlying(index)?;
                let first = if let Some(as_of) = as_of {
                    let (rows, _) = refdata.rate_future_schedule(index, trade_date()?)?;
                    Some(
                        rows.iter()
                            .find(|row| row.observed_at > as_of)
                            .ok_or_else(|| {
                                Error::MissingData(format!("no future clearing for {index}"))
                            })?
                            .date,
                    )
                } else {
                    None
                };
                Ok(StandardProduct::RateFuture(RateFuture::new(
                    index.clone(),
                    NonZeroU32::MIN,
                    RateFuturesPrice::new(value, basis.clone())?,
                    Position::Long,
                    first,
                )))
            }
            (Self::Swap { tenor, .. }, Convention::Irs(swap)) => Ok(StandardProduct::Swap(
                make_swap(swap, *tenor, value, QuotedCoupon::Fixed, refdata)?,
            )),
            (Self::BasisSwap { tenor, .. }, Convention::BasisSwap(swap))
            | (Self::XccyBasis { tenor, .. }, Convention::XccySwap(swap)) => {
                Ok(StandardProduct::Swap(make_swap(
                    swap,
                    *tenor,
                    value,
                    QuotedCoupon::PaySpread,
                    refdata,
                )?))
            }
            (Self::FxForward { tenor, .. }, Convention::FxForward(fx)) => {
                // スポットラグ、両通貨の休日の和と受渡の時刻はペアの慣行から読む。
                let spot = refdata.fx_convention(fx.pair)?;
                let calendar = refdata.calendar(&refdata.fx_calendar_spec(fx.pair)?)?;
                let (period, _) = fx_dates(trade_date()?, *tenor, spot.spot_lag, fx, &calendar)?;
                let settlement = refdata
                    .zone(&spot.settlement.zone)?
                    .resolve_convention(spot.settlement.on_date(period.end()))?;
                Ok(StandardProduct::FxForward(FxForward::new(
                    Money::new(fx.pair.base(), 1.0),
                    FxRate::forward(fx.pair, value, settlement)?,
                    Position::Long,
                )?))
            }
            _ => Err(Error::Internal(
                "validated template has another convention kind".into(),
            )),
        }
    }

    /// 校正商品の PV を読む通貨。FX はペアの quote、XCCY は払い脚の通貨。
    ///
    /// # Errors
    ///
    /// テンプレートに合う慣行がないか、慣行に校正できない丸めやスプレッドの規約がある場合。
    pub fn currency(&self, refdata: &RefData) -> Result<Currency> {
        let id = self.convention();
        let convention = refdata.convention(id).map_err(|error| {
            if matches!(self, Self::Ois { .. }) {
                error.context("template.ois.convention")
            } else {
                error.context("template.convention")
            }
        })?;
        match (self, convention) {
            (Self::Ois { .. }, Convention::Ois(ois)) => {
                check_observation(&ois.float, false, id)?;
                match refdata.index(&ois.float.index)?.kind {
                    IndexKind::OvernightRate { ccy, .. } => Ok(ccy),
                    _ => Err(Error::Internal(format!("{id} has no overnight index"))),
                }
            }
            (Self::Deposit { currency, .. }, Convention::Deposit(_)) => {
                refdata.currency(*currency)?;
                Ok(*currency)
            }
            (Self::Fra { .. }, Convention::Fra(fra)) => match refdata.index(&fra.index)?.kind {
                IndexKind::TermRate { ccy, .. } => Ok(ccy),
                _ => Err(Error::Internal(format!("{id} has no term index"))),
            },
            (Self::RateFuture { delivery_month, .. }, Convention::RateFuture(future)) => {
                let contract = future.contracts.get(delivery_month).ok_or_else(|| {
                    invalid(format!("unknown delivery month {delivery_month} for {id}"))
                })?;
                Ok(refdata
                    .rate_future_underlying(&contract.settlement_index)?
                    .1)
            }
            (Self::Swap { .. }, Convention::Irs(swap))
            | (Self::BasisSwap { .. }, Convention::BasisSwap(swap))
            | (Self::XccyBasis { .. }, Convention::XccySwap(swap)) => {
                for (i, leg) in swap.legs.iter().enumerate() {
                    if let SwapCouponConvention::Compounded(overnight) = &leg.coupon {
                        check_observation(
                            overnight,
                            i == 0 && !matches!(self, Self::Swap { .. }),
                            id,
                        )?;
                    }
                }
                Ok(swap.legs[0].currency)
            }
            (Self::FxForward { .. }, Convention::FxForward(fx)) => Ok(fx.pair.quote()),
            (Self::Ois { .. }, _) => Err(invalid(format!(
                "template.ois.convention: {id} is not an OIS convention"
            ))),
            _ => Err(invalid(format!(
                "template.convention: {id} is not a convention for this template"
            ))),
        }
    }

    /// テンプレートを参照データと照らして検証する。
    ///
    /// # Errors
    ///
    /// テナーが不正、慣行が適合しない、または観測に丸めがある場合。
    pub fn validate(&self, refdata: &RefData) -> Result<()> {
        self.currency(refdata).map(drop)
    }

    /// 二通貨のテンプレートで PV 換算が必要なもう一つの通貨を返す。
    ///
    /// # Errors
    ///
    /// テンプレートの慣行が参照データにない場合。
    pub fn other_currency(&self, refdata: &RefData) -> Result<Option<Currency>> {
        Ok(match (self, refdata.convention(self.convention())?) {
            (Self::FxForward { .. }, Convention::FxForward(fx)) => Some(fx.pair.base()),
            (Self::XccyBasis { .. }, Convention::XccySwap(swap)) => Some(swap.legs[1].currency),
            _ => None,
        })
    }
}

fn tenor_offset(length: TenorLength) -> Result<DateOffset> {
    let amount = i16::try_from(length.amount())
        .map_err(|_| invalid("template tenor exceeds date offset range"))?;
    Ok(DateOffset::new(amount, length.unit()))
}

fn fx_dates(
    trade_date: Date,
    tenor: Tenor,
    spot_lag: rhodium_core::refdata::FxSpotLag,
    fx: &rhodium_core::refdata::FxForwardConvention,
    calendar: &rhodium_core::calendar::Calendar,
) -> Result<(rhodium_core::date::DateInterval, Vec<rhodium_core::Warning>)> {
    money_market_period(
        trade_date,
        tenor,
        spot_lag.offset(),
        &AdvanceRule {
            calendar,
            bdc: fx.bdc,
            end_of_month: EndOfMonthRule::Ignore,
        },
    )
}

fn check_observation(
    overnight: &rhodium_core::refdata::OvernightLegConvention,
    quoted_spread: bool,
    id: &ConventionId,
) -> Result<()> {
    let rounding = overnight.observation.rounding();
    if rounding.daily.is_some() || rounding.aggregated.is_some() {
        return Err(invalid(format!(
            "{id} rounds observed rates and cannot define a calibration instrument"
        )));
    }
    if quoted_spread
        && !matches!(
            overnight.observation.aggregation(),
            Aggregation::Compound {
                spread: SpreadTreatment::Exclusive,
                ..
            }
        )
    {
        return Err(invalid(format!("{id} requires an exclusive quoted spread")));
    }
    Ok(())
}

#[derive(Clone, Copy)]
enum QuotedCoupon {
    Fixed,
    PaySpread,
}

fn make_swap(
    convention: &SwapConvention,
    tenor: TenorLength,
    quote: f64,
    quoted_coupon: QuotedCoupon,
    refdata: &RefData,
) -> Result<Swap> {
    let [first, second] = &convention.legs;
    let exchange = if first.currency == second.currency {
        NotionalExchange::None
    } else {
        NotionalExchange::InitialAndFinal
    };
    let leg = |i: usize, side: PayReceive| -> Result<Leg> {
        let spec = &convention.legs[i];
        let other = &convention.legs[1 - i];
        let basis = &spec.periods.day_count;
        let spread = if matches!(quoted_coupon, QuotedCoupon::PaySpread) && i == 0 {
            Spread::new(quote)?
        } else {
            Spread::new(0.0)?
        };
        let coupon = match &spec.coupon {
            SwapCouponConvention::Fixed => Coupon::Fixed {
                rate: Rate::new(quote, basis.clone(), Simple)?,
            },
            SwapCouponConvention::Term { index } => Coupon::Term {
                index: index.clone(),
                spread,
                gearing: 1.0,
            },
            SwapCouponConvention::Compounded(overnight) => Coupon::Compounded {
                index: overnight.index.clone(),
                spread,
                rule: overnight.observation.clone(),
            },
        };
        let notional = if let Some(fixing) = &spec.fx_reset {
            Notional::FxReset {
                base: Money::new(other.currency, 1.0),
                fixing: fixing.clone(),
            }
        } else {
            Notional::Fixed(Money::new(spec.currency, 1.0))
        };
        let schedule = ScheduleSpec::Generated {
            periods: GeneratedPeriods::new(
                StartRule::FromTradeDate(convention.spot_lag),
                EndRule::AfterStart(tenor),
                spec.periods.frequency,
                None,
                StubSpec::Initial(StubLength::Short),
                spec.periods.calendar.clone(),
                spec.periods.bdc,
                AccrualDates::Adjusted,
            )?,
            payment: spec.payment.clone(),
        };
        refdata.currency(spec.currency)?;
        Leg::with_notional(
            side,
            notional,
            schedule,
            basis.clone(),
            coupon,
            exchange,
            PaymentRounding::None,
        )
    };
    Swap::from_legs([leg(0, PayReceive::Pay)?, leg(1, PayReceive::Receive)?])
}
`````

`````
crates/rhodium-products/src/standard/wire.rs
//! 商品の serde 表現を構築 API に接続する。

use std::{fmt, num::NonZeroU32, str::FromStr};

use rhodium_core::{
    Error, Result,
    date::Date,
    day_count::DayCount,
    fx::{CcyPair, FxRate},
    id::{ConventionId, IndexId},
    money::Money,
    observation::ObservationRule,
    quantity::{Rate, RateFuturesPrice, Simple, Spread},
    schedule::{EndRule, ScheduleSpec},
    tenor::{Tenor, TenorLength},
    time::Instant,
};
use serde::{
    Deserialize, Deserializer, Serialize, Serializer,
    de::{Error as _, MapAccess, Visitor},
};

use super::{
    Coupon, Deposit, DepositDates, Fra, FxForward, FxSwap, Leg, Notional, NotionalExchange,
    PayReceive, PaymentRounding, Position, RateFuture, Swap, invalid, rates::SwapTerms,
};

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct RateFutureInput {
    settlement_index: IndexId,
    contracts: NonZeroU32,
    trade_price: f64,
    day_count: DayCount,
    first_clearing: Option<Date>,
    #[serde(default)]
    position: Position,
}

impl TryFrom<RateFutureInput> for RateFuture {
    type Error = Error;

    fn try_from(input: RateFutureInput) -> Result<Self> {
        Ok(Self::new(
            input.settlement_index,
            input.contracts,
            RateFuturesPrice::new(input.trade_price, input.day_count)?,
            input.position,
            input.first_clearing,
        ))
    }
}

impl From<RateFuture> for RateFutureInput {
    fn from(future: RateFuture) -> Self {
        Self {
            settlement_index: future.settlement_index,
            contracts: future.contracts,
            trade_price: future.trade_price.value(),
            day_count: future.trade_price.basis().clone(),
            first_clearing: future.first_clearing,
            position: future.position,
        }
    }
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct FraInput {
    notional: MoneyText,
    index: IndexId,
    fixing: Date,
    fixed_rate: ScaledRate,
    #[serde(default)]
    position: Position,
}

impl TryFrom<FraInput> for Fra {
    type Error = Error;

    fn try_from(input: FraInput) -> Result<Self> {
        Self::new(
            input.notional.0,
            input.index,
            input.fixing,
            input.fixed_rate.0,
            input.position,
        )
    }
}

impl From<Fra> for FraInput {
    fn from(fra: Fra) -> Self {
        Self {
            notional: MoneyText(fra.notional),
            index: fra.index,
            fixing: fra.fixing,
            fixed_rate: ScaledRate(fra.fixed_rate),
            position: fra.position,
        }
    }
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct FxForwardInput {
    pair: CcyPair,
    notional: MoneyText,
    rate: f64,
    settlement: Instant,
    #[serde(default)]
    position: Position,
}

impl TryFrom<FxForwardInput> for FxForward {
    type Error = Error;

    fn try_from(input: FxForwardInput) -> Result<Self> {
        Self::new(
            input.notional.0,
            FxRate::forward(input.pair, input.rate, input.settlement)?,
            input.position,
        )
    }
}

impl From<FxForward> for FxForwardInput {
    fn from(fx: FxForward) -> Self {
        Self {
            pair: fx.rate.pair(),
            notional: MoneyText(fx.notional),
            rate: fx.rate.rate(),
            settlement: fx.rate.settlement(),
            position: fx.position,
        }
    }
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct FxExchangeInput {
    rate: f64,
    settlement: Instant,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct FxSwapInput {
    pair: CcyPair,
    notional: MoneyText,
    near: FxExchangeInput,
    far: FxExchangeInput,
    #[serde(default)]
    position: Position,
}

impl TryFrom<FxSwapInput> for FxSwap {
    type Error = Error;

    fn try_from(input: FxSwapInput) -> Result<Self> {
        Self::new(
            input.notional.0,
            FxRate::forward(input.pair, input.near.rate, input.near.settlement)?,
            FxRate::forward(input.pair, input.far.rate, input.far.settlement)?,
            input.position,
        )
    }
}

impl From<FxSwap> for FxSwapInput {
    fn from(fx: FxSwap) -> Self {
        Self {
            pair: fx.near.rate.pair(),
            notional: MoneyText(fx.near.notional),
            near: FxExchangeInput {
                rate: fx.near.rate.rate(),
                settlement: fx.near.rate.settlement(),
            },
            far: FxExchangeInput {
                rate: fx.far.rate.rate(),
                settlement: fx.far.rate.settlement(),
            },
            position: fx.near.position,
        }
    }
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct DepositInput {
    notional: MoneyText,
    rate: ScaledRate,
    convention: Option<ConventionId>,
    day_count: Option<DayCount>,
    start: Option<Date>,
    maturity: Option<Date>,
    tenor: Option<Tenor>,
    #[serde(default)]
    position: Position,
    #[serde(default)]
    payment_rounding: PaymentRounding,
}

impl TryFrom<DepositInput> for Deposit {
    type Error = Error;

    fn try_from(input: DepositInput) -> Result<Self> {
        let dates = match (input.start, input.maturity, input.tenor) {
            (Some(start), Some(maturity), None) => DepositDates::Confirmed { start, maturity },
            (None, None, Some(tenor)) => DepositDates::Market(tenor),
            _ => return Err(invalid("deposit requires start and maturity, or a tenor")),
        };
        Self::new(
            input.notional.0,
            input.rate.0,
            dates,
            input.convention,
            input.day_count,
            input.position,
            input.payment_rounding,
        )
    }
}

impl From<Deposit> for DepositInput {
    fn from(deposit: Deposit) -> Self {
        let (start, maturity, tenor) = match deposit.dates {
            DepositDates::Confirmed { start, maturity } => (Some(start), Some(maturity), None),
            DepositDates::Market(tenor) => (None, None, Some(tenor)),
        };
        Self {
            notional: MoneyText(deposit.notional),
            rate: ScaledRate(deposit.rate),
            convention: deposit.convention,
            day_count: deposit.day_count,
            start,
            maturity,
            tenor,
            position: deposit.position,
            payment_rounding: deposit.payment_rounding,
        }
    }
}

/// スワップの入力。`legs` があれば二つの脚で、無ければ慣行の短縮形（OIS）で書いたものとして読む。
///
/// 形ごとの untagged の enum にすると、どの形にも一致しなかったことしか報告できず、
/// `fixed_rate = "invalid_rate"` のようなフィールドの誤りが分からない。そこで一つの struct で
/// フィールドを読み、値の誤りにはフィールド名を付け、形は読み終えてから `legs` の有無で決める。
#[derive(Default, Serialize)]
pub(super) struct SwapInput {
    #[serde(skip_serializing_if = "Option::is_none")]
    legs: Option<Box<[Leg; 2]>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    convention: Option<ConventionId>,
    #[serde(skip_serializing_if = "Option::is_none")]
    start: Option<Date>,
    #[serde(skip_serializing_if = "Option::is_none")]
    maturity: Option<Date>,
    #[serde(skip_serializing_if = "Option::is_none")]
    tenor: Option<TenorLength>,
    #[serde(skip_serializing_if = "Option::is_none")]
    notional: Option<MoneyText>,
    #[serde(skip_serializing_if = "Option::is_none")]
    fixed_rate: Option<ScaledRate>,
    #[serde(skip_serializing_if = "Option::is_none")]
    fixed_side: Option<PayReceive>,
    #[serde(skip_serializing_if = "Option::is_none")]
    spread: Option<ScaledRate>,
    #[serde(skip_serializing_if = "Option::is_none")]
    payment_rounding: Option<PaymentRounding>,
}

const SWAP_FIELDS: &[&str] = &[
    "legs",
    "convention",
    "start",
    "maturity",
    "tenor",
    "notional",
    "fixed_rate",
    "fixed_side",
    "spread",
    "payment_rounding",
];

impl<'de> Deserialize<'de> for SwapInput {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> std::result::Result<Self, D::Error> {
        deserializer.deserialize_map(SwapVisitor)
    }
}

struct SwapVisitor;

impl<'de> Visitor<'de> for SwapVisitor {
    type Value = SwapInput;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("swap legs or OIS convention terms")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> std::result::Result<SwapInput, A::Error> {
        let mut input = SwapInput::default();
        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "legs" => field(&mut map, &key, &mut input.legs)?,
                "convention" => field(&mut map, &key, &mut input.convention)?,
                "start" => field(&mut map, &key, &mut input.start)?,
                "maturity" => field(&mut map, &key, &mut input.maturity)?,
                "tenor" => field(&mut map, &key, &mut input.tenor)?,
                "notional" => field(&mut map, &key, &mut input.notional)?,
                "fixed_rate" => field(&mut map, &key, &mut input.fixed_rate)?,
                "fixed_side" => field(&mut map, &key, &mut input.fixed_side)?,
                "spread" => field(&mut map, &key, &mut input.spread)?,
                "payment_rounding" => field(&mut map, &key, &mut input.payment_rounding)?,
                other => return Err(A::Error::unknown_field(other, SWAP_FIELDS)),
            }
        }
        Ok(input)
    }
}

/// 一つのフィールドの値を読み、誤りにフィールド名を付ける。
///
/// 商品は `type` で選ぶ内部タグの enum の中にあり、serde は値を一度バッファしてから読むので、
/// 入力形式（TOML など）のキーの経路は商品の位置までしか届かない。
/// 他の入力の型の `Option` のフィールドと同じく、JSON の `null` は書かなかったものとして扱う
/// （以前の版は省略した `start` などを `null` として書いた）。
fn field<'de, A: MapAccess<'de>, T: Deserialize<'de>>(
    map: &mut A,
    key: &str,
    slot: &mut Option<T>,
) -> std::result::Result<(), A::Error> {
    if slot.is_some() {
        return Err(A::Error::custom(format_args!("duplicate field `{key}`")));
    }
    *slot = map
        .next_value::<Option<T>>()
        .map_err(|error| A::Error::custom(format_args!("`{key}`: {error}")))?;
    Ok(())
}

impl TryFrom<SwapInput> for Swap {
    type Error = Error;

    fn try_from(input: SwapInput) -> Result<Self> {
        let SwapInput {
            legs,
            convention,
            start,
            maturity,
            tenor,
            notional,
            fixed_rate,
            fixed_side,
            spread,
            payment_rounding,
        } = input;
        if let Some(legs) = legs {
            // 脚で書いたスワップは、短縮形のフィールドを持たない（脚の条件と食い違いうるため）。
            let terms = [
                ("convention", convention.is_some()),
                ("start", start.is_some()),
                ("maturity", maturity.is_some()),
                ("tenor", tenor.is_some()),
                ("notional", notional.is_some()),
                ("fixed_rate", fixed_rate.is_some()),
                ("fixed_side", fixed_side.is_some()),
                ("spread", spread.is_some()),
                ("payment_rounding", payment_rounding.is_some()),
            ];
            if let Some((name, _)) = terms.iter().find(|(_, present)| *present) {
                return Err(invalid(format!(
                    "swap with legs does not take the convention field `{name}`"
                )));
            }
            return Self::from_legs(*legs);
        }
        let required = |name: &str| {
            invalid(format!(
                "OIS requires `{name}` (or write the swap with `legs`)"
            ))
        };
        let convention = convention.ok_or_else(|| required("convention"))?;
        let end = match (maturity, tenor) {
            (Some(date), None) => EndRule::Date(date),
            (None, Some(tenor)) => EndRule::AfterStart(tenor),
            _ => return Err(invalid("OIS requires exactly one of maturity and tenor")),
        };
        Self::from_convention(
            convention,
            start,
            end,
            notional.ok_or_else(|| required("notional"))?.0,
            fixed_rate.ok_or_else(|| required("fixed_rate"))?.0,
            fixed_side.ok_or_else(|| required("fixed_side"))?,
            Spread::new(spread.map_or(0.0, |spread| spread.0))?,
            payment_rounding.unwrap_or_default(),
        )
    }
}

impl From<Swap> for SwapInput {
    fn from(swap: Swap) -> Self {
        match swap.0 {
            SwapTerms::Legs(legs) => Self {
                legs: Some(legs),
                ..Self::default()
            },
            SwapTerms::Convention(ois) => {
                let (maturity, tenor) = match ois.end {
                    EndRule::Date(date) => (Some(date), None),
                    EndRule::AfterStart(tenor) => (None, Some(tenor)),
                };
                Self {
                    legs: None,
                    convention: Some(ois.convention),
                    start: ois.start,
                    maturity,
                    tenor,
                    notional: Some(MoneyText(ois.notional)),
                    fixed_rate: Some(ScaledRate(ois.fixed_rate)),
                    fixed_side: Some(ois.fixed_side),
                    spread: Some(ScaledRate(ois.spread.value())),
                    payment_rounding: Some(ois.payment_rounding),
                }
            }
        }
    }
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct LegInput {
    side: PayReceive,
    notional: NotionalInput,
    schedule: ScheduleSpec,
    day_count: DayCount,
    coupon: CouponInput,
    #[serde(default)]
    exchange: NotionalExchange,
    #[serde(default)]
    payment_rounding: PaymentRounding,
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
enum NotionalInput {
    Fixed(MoneyText),
    FxReset { base: MoneyText, fixing: IndexId },
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
enum CouponInput {
    Fixed {
        rate: ScaledRate,
    },
    Term {
        index: IndexId,
        #[serde(default)]
        spread: ScaledRate,
        #[serde(default = "unit_gearing")]
        gearing: f64,
    },
    Compounded {
        index: IndexId,
        rule: ObservationRule,
        #[serde(default)]
        spread: ScaledRate,
    },
}

const fn unit_gearing() -> f64 {
    1.0
}

impl TryFrom<LegInput> for Leg {
    type Error = Error;

    fn try_from(input: LegInput) -> Result<Self> {
        let coupon = match input.coupon {
            CouponInput::Fixed { rate } => Coupon::Fixed {
                rate: Rate::new(rate.0, input.day_count.clone(), Simple)?,
            },
            CouponInput::Term {
                index,
                spread,
                gearing,
            } => Coupon::Term {
                index,
                spread: Spread::new(spread.0)?,
                gearing,
            },
            CouponInput::Compounded {
                index,
                rule,
                spread,
            } => Coupon::Compounded {
                index,
                rule,
                spread: Spread::new(spread.0)?,
            },
        };
        Self::with_notional(
            input.side,
            match input.notional {
                NotionalInput::Fixed(money) => Notional::Fixed(money.0),
                NotionalInput::FxReset { base, fixing } => Notional::FxReset {
                    base: base.0,
                    fixing,
                },
            },
            input.schedule,
            input.day_count,
            coupon,
            input.exchange,
            input.payment_rounding,
        )
    }
}

impl From<Leg> for LegInput {
    fn from(leg: Leg) -> Self {
        let coupon = match leg.coupon {
            Coupon::Fixed { rate } => CouponInput::Fixed {
                rate: ScaledRate(rate.value()),
            },
            Coupon::Term {
                index,
                spread,
                gearing,
            } => CouponInput::Term {
                index,
                spread: ScaledRate(spread.value()),
                gearing,
            },
            Coupon::Compounded {
                index,
                rule,
                spread,
            } => CouponInput::Compounded {
                index,
                rule,
                spread: ScaledRate(spread.value()),
            },
        };
        Self {
            side: leg.side,
            notional: match leg.notional {
                Notional::Fixed(amount) => NotionalInput::Fixed(MoneyText(amount)),
                Notional::FxReset { base, fixing } => NotionalInput::FxReset {
                    base: MoneyText(base),
                    fixing,
                },
            },
            schedule: leg.schedule,
            day_count: leg.day_count,
            coupon,
            exchange: leg.exchange,
            payment_rounding: leg.payment_rounding,
        }
    }
}

struct MoneyText(Money);

impl FromStr for MoneyText {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        text.parse()
            .map(Self)
            .map_err(|error: Error| error.context("notional"))
    }
}

impl fmt::Display for MoneyText {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.0.amount(), self.0.ccy())
    }
}

#[derive(Default)]
struct ScaledRate(f64);

impl FromStr for ScaledRate {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        Ok(Self(text.parse::<Spread>()?.value()))
    }
}

impl fmt::Display for ScaledRate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // f64 の乗算を使わず、最短の十進表記から小数点を二桁動かす。
        let value = self.0.to_string();
        let (sign, value) = value
            .strip_prefix('-')
            .map_or(("", value.as_str()), |s| ("-", s));
        let (whole, fraction) = value.split_once('.').unwrap_or((value, ""));
        let (prefix, rest) = fraction.split_at(fraction.len().min(2));
        let whole = format!("{whole}{prefix}{}", &"00"[prefix.len()..]);
        let whole = whole.trim_start_matches('0');
        let whole = if whole.is_empty() { "0" } else { whole };
        if rest.is_empty() {
            write!(f, "{sign}{whole}%")
        } else {
            write!(f, "{sign}{whole}.{rest}%")
        }
    }
}

macro_rules! serde_text {
    ($type:ty) => {
        impl Serialize for $type {
            fn serialize<S: Serializer>(
                &self,
                serializer: S,
            ) -> std::result::Result<S::Ok, S::Error> {
                serializer.collect_str(self)
            }
        }

        impl<'de> Deserialize<'de> for $type {
            fn deserialize<D: Deserializer<'de>>(
                deserializer: D,
            ) -> std::result::Result<Self, D::Error> {
                String::deserialize(deserializer)?
                    .parse()
                    .map_err(D::Error::custom)
            }
        }
    };
}

serde_text!(MoneyText);
serde_text!(ScaledRate);
`````

`````
crates/rhodium-products/src/trade.rs
//! 取引のライフサイクルに関わる値。

use std::collections::{BTreeMap, BTreeSet};

use rhodium_core::{
    Error, Result, Warning,
    currency::Asset,
    date::Date,
    id::{BookId, IndexId, MarginSetId, TradeId},
    money::Money,
    refdata::{FallbackApplicability, FallbackReplacement, FallbackTrigger, RefData},
    time::Instant,
};

use crate::{
    ir::{Contract, Exercises, ObservationCache, Product},
    standard::{PayReceive, StandardProduct},
};

/// 商品の条件と取引に記録された行使、手数料、フォールバック。
#[derive(Debug, Clone)]
pub struct Trade {
    /// 取引の識別子。
    pub id: TradeId,
    /// 商品の条件。
    pub product: ProductSpec,
    /// 所属するマージンセット。
    pub margin_set: MarginSetId,
    /// 集計に使うブック。
    pub book: BookId,
    /// 確認書の約定日。
    pub trade_date: Date,
    /// 集計や価格設定の照合に使う属性。
    pub attributes: BTreeMap<String, String>,
    /// 取引に記録された行使。
    pub exercises: Vec<ExerciseRecord>,
    /// 商品の受渡に加える手数料。
    pub fees: Vec<Fee>,
    /// 取引上のフォールバックの指定と適用条件。
    pub fallback: FallbackTerms,
}

/// 取引の条件を IR に変換する入力。
#[derive(Debug, Clone)]
pub enum ProductSpec {
    /// 型付きの標準商品。
    Standard(StandardProduct),
    /// Rust の builder で組み立てた契約。
    Contract(Contract),
    /// 指定の瞬間以後に効力を持つ商品の変更。
    Amended {
        /// 変更前の条件。
        before: Box<ProductSpec>,
        /// 変更に含む最初の瞬間。
        from: Instant,
        /// 変更の内容。
        after: Amendment,
    },
}

/// 効力の瞬間以後の商品への変更。
#[derive(Debug, Clone)]
pub enum Amendment {
    /// 受渡の最終額に掛ける係数。0 以上 1 未満。
    Scale(f64),
    /// 新しい商品で条件を置き換える。
    Replace(Box<ProductSpec>),
}

/// 効力の瞬間に独立して支払う手数料。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fee {
    /// 正の手数料金額。
    pub amount: Money,
    /// 手数料の受け払い。
    pub side: PayReceive,
    /// 決済の瞬間。
    pub at: Instant,
}

/// 同じ参照データで複数の取引を解決し、RFR の観測列と公表時刻を共有する。
#[derive(Debug)]
pub struct TradeProductResolver<'r> {
    refdata: &'r RefData,
    observations: ObservationCache,
}

impl<'r> TradeProductResolver<'r> {
    /// 読み込みに使う参照データに共有の範囲を固定する。
    #[must_use]
    pub fn new(refdata: &'r RefData) -> Self {
        Self {
            refdata,
            observations: ObservationCache::default(),
        }
    }

    /// 取引を検証済み IR にする。複数のワーカーから呼び出せる。
    ///
    /// # Errors
    ///
    /// [`Trade::to_product`] と同じ失敗。
    pub fn resolve(&self, trade: &Trade) -> Result<Product> {
        trade.to_product_using(self.refdata, Some(&self.observations))
    }
}

impl Trade {
    /// 価格設定の方式の表が照合する商品の種類。
    ///
    /// 型付きの標準商品だけが種類を持つ。builder の契約と変更した商品（一部解約、置き換え）は
    /// 標準商品の条件と一致しないので種類を持たず、方式の表の行に当たらない。
    #[must_use]
    pub fn product_kind(&self) -> Option<crate::standard::ProductKind> {
        match &self.product {
            ProductSpec::Standard(product) => Some(product.kind()),
            ProductSpec::Contract(_) | ProductSpec::Amended { .. } => None,
        }
    }

    /// 商品と手数料を検証済み IR にし、評価に使う行使記録も格納する。
    ///
    /// # Errors
    ///
    /// 参照、契約の因果性、変更の条件や行使記録が不正ならエラー。
    /// 分割前に確定する必要がある選択に記録が無ければ `Unsupported`。
    pub fn to_product(&self, refdata: &RefData) -> Result<Product> {
        self.to_product_using(refdata, None)
    }

    fn to_product_using(
        &self,
        refdata: &RefData,
        cache: Option<&ObservationCache>,
    ) -> Result<Product> {
        self.fallback.validate(refdata)?;
        let mut exercises = Exercises::new(&self.exercises)?;
        let mut warnings = Vec::new();
        let mut contract =
            self.product
                .lower(self.trade_date, refdata, &mut exercises, &mut warnings)?;
        exercises.finish()?;
        if !self.fees.is_empty() {
            let fees = Contract::build(|builder| {
                let payments = self
                    .fees
                    .iter()
                    .map(|fee| {
                        let amount = fee.amount.amount();
                        if !amount.is_finite() || amount <= 0.0 {
                            return Err(Error::InvalidData(
                                "fee must be positive and finite".into(),
                            ));
                        }
                        let signed = match fee.side {
                            PayReceive::Pay => -amount,
                            PayReceive::Receive => amount,
                        };
                        let signed = builder.constant(signed);
                        Ok(builder.pay(Asset::Cash(fee.amount.ccy()), signed, fee.at))
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(builder.all(payments))
            })?;
            contract = contract.combine(fees)?;
        }
        Ok(Product::with_fallback(
            contract,
            refdata,
            warnings,
            &self.fallback,
            self.trade_date,
            cache,
        )?
        .with_trade_exercises(&self.exercises))
    }

    /// 行使記録を加えた取引を返す。
    #[must_use]
    pub fn exercised(&self, record: ExerciseRecord) -> Self {
        let mut next = self.clone();
        next.exercises.push(record);
        next
    }

    /// 現在の持分の一部を解約した取引を返す。
    ///
    /// # Errors
    ///
    /// 解約する割合が 0 より大きく 1 未満でなければ `InvalidInput`。
    pub fn partially_terminated(
        &self,
        from: Instant,
        fraction: f64,
        fee: Option<Fee>,
    ) -> Result<Self> {
        let remaining = 1.0 - fraction;
        if !(fraction > 0.0 && remaining > 0.0 && remaining < 1.0) {
            return Err(Error::InvalidInput(
                "partial termination requires a representable fraction in (0, 1)".into(),
            ));
        }
        Ok(self.changed(from, Amendment::Scale(remaining), fee))
    }

    /// 全部を解約した取引を返す。
    #[must_use]
    pub fn terminated(&self, from: Instant, fee: Option<Fee>) -> Self {
        self.changed(from, Amendment::Scale(0.0), fee)
    }

    /// 商品を置き換えた取引を返す。
    #[must_use]
    pub fn amended(&self, from: Instant, product: ProductSpec) -> Self {
        self.changed(from, Amendment::Replace(Box::new(product)), None)
    }

    fn changed(&self, from: Instant, after: Amendment, fee: Option<Fee>) -> Self {
        let mut next = self.clone();
        next.product = ProductSpec::Amended {
            before: Box::new(next.product),
            from,
            after,
        };
        next.fees.extend(fee);
        next
    }
}

impl ProductSpec {
    fn lower(
        &self,
        trade_date: Date,
        refdata: &RefData,
        exercises: &mut Exercises<'_>,
        warnings: &mut Vec<Warning>,
    ) -> Result<Contract> {
        let contract = match self {
            Self::Standard(product) => {
                let (contract, notices) = product.build_contract(trade_date, refdata)?;
                warnings.extend(notices);
                contract
            }
            Self::Contract(contract) => contract.clone(),
            Self::Amended {
                before,
                from,
                after,
            } => {
                let before = before.lower(trade_date, refdata, exercises, warnings)?;
                return match after {
                    Amendment::Scale(factor) => before.scale_from(*from, *factor),
                    Amendment::Replace(product) => {
                        let new = product.lower(trade_date, refdata, exercises, warnings)?;
                        before.replace_from(new, *from, exercises)
                    }
                };
            }
        };
        exercises.observe(&contract)?;
        Ok(contract)
    }
}

/// 一つの選択について、取引の行使記録が効力を持ったことを表す。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExerciseRecord {
    /// 商品の `Choose` のラベル。
    pub label: String,
    /// 商品に定めた判断の締切。
    pub at: Instant,
    /// 記録が効力を持った瞬間。
    pub effective_at: Instant,
    /// 選択肢のゼロ始まりの位置。
    pub option: usize,
}

/// 確認書と当事者の合意によるフォールバックの適用条件。
#[derive(Debug, Clone, Default, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FallbackTerms {
    /// 確認書が指定する定義書の版。
    pub definitions: Option<String>,
    /// 両当事者が批准したプロトコル。
    #[serde(default)]
    pub protocols: BTreeSet<String>,
    /// 個別に既定のフォールバックに合意した指標。
    #[serde(default)]
    pub bilateral: BTreeSet<IndexId>,
    /// 指標ごとの取引固有の指定。
    #[serde(default)]
    pub by_index: BTreeMap<IndexId, FallbackSelection>,
}

/// 既定の適用条件に優先する、取引のフォールバック指定。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum FallbackSelection {
    /// 指標の既定のフォールバックを契約に適用する。
    IndexDefault,
    /// フォールバックを適用しない。
    Disabled,
    /// 取引に固有の代替を適用する。
    Custom {
        /// この取引で適用する開始事由。
        trigger: FallbackTrigger,
        /// この取引に固有の代替の規則。
        replacement: FallbackReplacement,
    },
}

impl FallbackTerms {
    fn validate(&self, refdata: &RefData) -> Result<()> {
        for id in &self.bilateral {
            refdata.index(id)?;
        }
        for (id, selection) in &self.by_index {
            if !matches!(
                refdata.index(id)?.kind,
                rhodium_core::refdata::IndexKind::TermRate { .. }
            ) {
                return Err(Error::InvalidData(format!(
                    "fallback index {id} is not a term rate"
                )));
            }
            if let FallbackSelection::Custom { replacement, .. } = selection {
                refdata.validate_fallback_replacement(id, replacement)?;
            }
        }
        Ok(())
    }

    pub(crate) fn select<'a>(
        &'a self,
        id: &IndexId,
        trade_date: Date,
        refdata: &'a RefData,
    ) -> Option<(FallbackTrigger, &'a FallbackReplacement)> {
        match self.by_index.get(id) {
            Some(FallbackSelection::Disabled) => None,
            Some(FallbackSelection::Custom {
                trigger,
                replacement,
            }) => Some((*trigger, replacement)),
            Some(FallbackSelection::IndexDefault) => {
                let definition = refdata.fallback(id)?;
                Some((definition.trigger, &definition.replacement))
            }
            None => {
                let definition = refdata.fallback(id)?;
                definition
                    .applies
                    .iter()
                    .any(|applies| match applies {
                        FallbackApplicability::Definitions {
                            version,
                            traded_on_or_after,
                        } => {
                            self.definitions.as_ref() == Some(version)
                                && trade_date >= *traded_on_or_after
                        }
                        FallbackApplicability::Protocol { protocol } => {
                            self.protocols.contains(protocol)
                        }
                        FallbackApplicability::Bilateral => self.bilateral.contains(id),
                    })
                    .then_some((definition.trigger, &definition.replacement))
            }
        }
    }
}
`````

`````
crates/rhodium-products/tests/compile_fail.rs
//! コア IR の添字を型で区別する。

use rhodium_testkit::compile_fail_pair;

#[test]
fn contract_index_cannot_cross_builder_boundaries() {
    // 別の builder で作った契約は子契約にできない。
    compile_fail_pair!("cross_builder_contract");
}

#[test]
fn expression_index_cannot_cross_builder_boundaries() {
    // 別の builder の式を同じ算術式に混ぜられない。
    compile_fail_pair!("cross_builder_expression");
}

#[test]
fn boolean_expression_cannot_be_paid() {
    // 真偽値の式は支払額に使えない。
    compile_fail_pair!("boolean_payment");
}
`````

`````
crates/rhodium-products/tests/settlement_boundary.rs
//! 評価時点ちょうどの受渡を含めるかを確かめる。

use rhodium_core::{Result, time::Instant};
use rhodium_products::ir::SettlementBoundary::{self, Exclude, Include};

#[test]
fn settlement_boundary_preserves_legacy_cutoff_meaning() -> Result<()> {
    // 旧 Inclusive は cutoff と等しい受渡を落とす。旧 Exclusive は等しい受渡を残す。
    let as_of: Instant = "2026-11-20T00:00:00Z".parse()?;
    let before = Instant::from_unix_micros(as_of.unix_micros() - 1)?;
    let after = Instant::from_unix_micros(as_of.unix_micros() + 1)?;
    assert_eq!(
        [before, as_of, after].map(|at| Exclude.includes(as_of, at)),
        [false, false, true]
    );
    assert_eq!(
        [before, as_of, after].map(|at| Include.includes(as_of, at)),
        [false, true, true]
    );
    assert_eq!(SettlementBoundary::default(), Exclude);
    Ok(())
}
`````

`````
crates/rhodium-products/tests/ui/boolean_payment.rs
use rhodium_core::{
    Result,
    currency::{Asset, Currency},
};
use rhodium_products::ir::Contract;

fn main() -> Result<()> {
    let at = "2027-01-04T00:00:00Z".parse()?;
    let _ = Contract::build(|builder| {
        let amount = builder.boolean(false);
        Ok(builder.pay(Asset::Cash(Currency::JPY), amount, at))
    })?;
    Ok(())
}
`````

`````
crates/rhodium-products/tests/ui/boolean_payment_ok.rs
use rhodium_core::{
    Result,
    currency::{Asset, Currency},
};
use rhodium_products::ir::Contract;

fn main() -> Result<()> {
    let at = "2027-01-04T00:00:00Z".parse()?;
    let _ = Contract::build(|builder| {
        let amount = builder.constant(0.0);
        Ok(builder.pay(Asset::Cash(Currency::JPY), amount, at))
    })?;
    Ok(())
}
`````

`````
crates/rhodium-products/tests/ui/cross_builder_contract.rs
use rhodium_core::Result;
use rhodium_products::ir::Contract;

fn main() -> Result<()> {
    let _ = Contract::build(|outer| {
        let child = outer.zero();
        let _ = Contract::build(|inner| Ok(inner.all([child])))?;
        Ok(child)
    })?;
    Ok(())
}
`````

`````
crates/rhodium-products/tests/ui/cross_builder_contract_ok.rs
use rhodium_core::Result;
use rhodium_products::ir::Contract;

fn main() -> Result<()> {
    let _ = Contract::build(|outer| {
        let child = outer.zero();
        let _ = Contract::build(|inner| {
            let local = inner.zero();
            Ok(inner.all([local]))
        })?;
        Ok(child)
    })?;
    Ok(())
}
`````

`````
crates/rhodium-products/tests/ui/cross_builder_expression.rs
use rhodium_core::{
    Result,
    currency::{Asset, Currency},
};
use rhodium_products::ir::Contract;

fn main() -> Result<()> {
    let at = "2027-01-04T00:00:00Z".parse()?;
    let _ = Contract::build(|outer| {
        let amount = outer.constant(42.0);
        let _ = Contract::build(|inner| {
            let local = inner.constant(1.0);
            let sum = inner.add([amount, local]);
            Ok(inner.pay(Asset::Cash(Currency::JPY), sum, at))
        })?;
        Ok(outer.pay(Asset::Cash(Currency::JPY), amount, at))
    })?;
    Ok(())
}
`````

`````
crates/rhodium-products/tests/ui/cross_builder_expression_ok.rs
use rhodium_core::{
    Result,
    currency::{Asset, Currency},
};
use rhodium_products::ir::Contract;

fn main() -> Result<()> {
    let at = "2027-01-04T00:00:00Z".parse()?;
    let _ = Contract::build(|outer| {
        let amount = outer.constant(42.0);
        let _ = Contract::build(|inner| {
            let first = inner.constant(42.0);
            let second = inner.constant(1.0);
            let sum = inner.add([first, second]);
            Ok(inner.pay(Asset::Cash(Currency::JPY), sum, at))
        })?;
        Ok(outer.pay(Asset::Cash(Currency::JPY), amount, at))
    })?;
    Ok(())
}
`````

`````
crates/rhodium-risk/Cargo.toml
[package]
name = "rhodium-risk"
description = "価格の感応度、シナリオと規制計算"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
publish.workspace = true

[lints]
workspace = true

[dependencies]
rhodium-core.workspace = true
rhodium-market.workspace = true
rhodium-pricing.workspace = true
scioxide-ad.workspace = true
scioxide-core.workspace = true

[dev-dependencies]
proptest.workspace = true
rhodium-products.workspace = true
rhodium-testkit.workspace = true
serde.workspace = true
serde_json = "1"
toml = "1.2"
`````

`````
crates/rhodium-risk/src/adjoint.rs
//! 保存した校正の線形化から、自然な単位の quote 感応度を求める。

use std::{collections::BTreeMap, ops::Range};

use rhodium_core::{Error, Result, currency::Currency};
use rhodium_market::{
    quote::QuoteId,
    structure::{BlockId, ParamLayout},
};
use rhodium_pricing::calibration::{CalibratedMarket, GroupOutcome};
use scioxide_core::AdjointOperator;

/// 通貨ごとに集計した平坦な θ 勾配を、quote の自然な単位あたりの微分に戻す。
pub(crate) fn solve(
    cm: &CalibratedMarket,
    by_ccy: BTreeMap<Currency, Vec<f64>>,
) -> Result<Vec<(Currency, BTreeMap<QuoteId, f64>)>> {
    let layout = cm.market().structure().layout();
    let mut columns: Vec<_> = by_ccy
        .into_iter()
        .map(|(ccy, theta)| (ccy, theta, BTreeMap::new()))
        .collect();
    let mut rhs = Vec::new();
    let mut lambda = Vec::new();
    let mut contributions = Vec::new();
    let mut scratch = Vec::new();

    // 下流で生じた θ の寄与を、上流グループの右辺を解く前に足す。
    for group in cm.groups().iter().rev() {
        let own = ranges(layout, group.blocks())?;
        let upstream = ranges(layout, group.upstream())?;
        let n_upstream = upstream.iter().map(Range::len).sum::<usize>();
        for (ccy, theta, quotes) in &mut columns {
            rhs.clear();
            rhs.extend(own.iter().flat_map(|range| &theta[range.clone()]).copied());
            if rhs.iter().all(|&value| value == 0.0) {
                continue;
            }
            let GroupOutcome::Solved(solution) = group.outcome() else {
                return Err(Error::Internal(format!(
                    "nonzero {ccy} adjoint for unavailable calibration group {:?}",
                    group.blocks()
                )));
            };

            // λ = G_θ⁻ᵀ g。
            apply_adjoint(
                &solution.dg_dtheta().inverse_operator(),
                &rhs,
                &mut lambda,
                rhs.len(),
                &mut scratch,
                "calibration adjoint solve",
            )?;

            // ∂V/∂q = −G_qᵀλ。同じ quote を読む複数のグループの寄与を足す。
            apply_adjoint(
                solution.dg_dq(),
                &lambda,
                &mut contributions,
                group.quotes().len(),
                &mut scratch,
                "quote adjoint",
            )?;
            for (id, &value) in group.quotes().iter().zip(&contributions) {
                if value != 0.0 {
                    let total = quotes.entry(id.clone()).or_insert(0.0);
                    *total -= value;
                    if !total.is_finite() {
                        return Err(Error::Numerical(format!(
                            "non-finite adjoint for quote {id}"
                        )));
                    }
                }
            }

            if n_upstream != 0 {
                apply_adjoint(
                    solution.dg_dup(),
                    &lambda,
                    &mut contributions,
                    n_upstream,
                    &mut scratch,
                    "upstream adjoint",
                )?;
                let indices = upstream.iter().flat_map(Clone::clone);
                for (index, &value) in indices.zip(&contributions) {
                    theta[index] -= value;
                    if !theta[index].is_finite() {
                        return Err(Error::Numerical(format!(
                            "non-finite upstream adjoint at theta index {index}"
                        )));
                    }
                }
            }
        }
    }
    Ok(columns
        .into_iter()
        .map(|(ccy, _, quotes)| (ccy, quotes))
        .collect())
}

/// `op` の随伴 `opᵀ · input` を、長さ `len` にした `output` に書く。
///
/// `len` はグループの記録（ブロック、quote、上流の並び）から決めた出力の長さで、
/// 演算子の列数と食い違えば scioxide が次元の不一致として返す。
/// 作業領域の長さは分解や疎行列の実装ごとに違う（LU は n、LDL は 2n など）ため、
/// 推測せずに演算子に問い合わせる。
fn apply_adjoint<A: AdjointOperator<f64> + ?Sized>(
    op: &A,
    input: &[f64],
    output: &mut Vec<f64>,
    len: usize,
    scratch: &mut Vec<f64>,
    what: &str,
) -> Result<()> {
    let numerical = |error| Error::Numerical(format!("{what}: {error}"));
    output.resize(len, 0.0);
    scratch.resize(op.apply_adjoint_workspace_len().map_err(numerical)?, 0.0);
    op.apply_adjoint_into(input, output, scratch)
        .map_err(numerical)
}

fn ranges(layout: &ParamLayout, blocks: &[BlockId]) -> Result<Vec<Range<usize>>> {
    blocks
        .iter()
        .map(|block| {
            layout.range(block).ok_or_else(|| {
                Error::Internal(format!("missing theta range for curve {}", block.curve()))
            })
        })
        .collect()
}
`````

`````
crates/rhodium-risk/src/delta.rs
//! 受渡通貨別または報告通貨で集計した、校正入力 quote の感応度。

use std::collections::{BTreeMap, btree_map::Entry};

use rhodium_core::{
    Error, Result,
    currency::Currency,
    fx::CcyPair,
    money::MoneyBag,
    parallel::{self, ExecutionContext},
};
use rhodium_market::{
    quote::{QuoteId, QuoteKind},
    structure::ParamLayout,
};
use rhodium_pricing::{calibration::CalibratedMarket, prepare::Prepared};
