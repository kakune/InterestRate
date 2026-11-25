        &updated.calendars,
    )?;
    // BUS/252 ではどちらの区間も営業日が 1 日ずつ。
    assert_close(observations.fixings()[0].1, 1.0 / 252.0);
    assert_close(observations.fixings()[1].1, 1.0 / 252.0);
    assert_close(observations.annualization().value(), 2.0 / 252.0);
    assert_eq!(observations.annualization().basis(), &basis);
    Ok(())
}

#[test]
fn calendar_days_months_and_index_endpoints_follow_the_rule() -> Result<()> {
    let dates = Dates::new(&[])?;
    let schedule = dates.once("2027-01-08", "2027-01-11")?;
    let period = &schedule.periods()[0];
    let all_days = ObservationDays::AllDays {
        fixing_calendar: dates.spec(),
    };
    let rule = ObservationRule::new(
        all_days,
        ObservationMethod::Plain,
        WEIGHTED,
        ObservationRounding::default(),
    )?;
    // 金・土・日の行をそれぞれ残し、休日の行にも金曜の fixing を置く。
    check_daily(&observe(&dates, period, &rule)?, &[("2027-01-08", 1.0); 3])?;
    let mut month_rule = ObservationRule::new(
        ObservationDays::Months,
        ObservationMethod::Plain,
        WEIGHTED,
        ObservationRounding::default(),
    )?;
    let months = dates.once("2027-01-31", "2027-03-02")?;
    let monthly = observe(&dates, &months.periods()[0], &month_rule)?;
    // 1 月と 3 月は 1 日ずつ、2 月は 28 日を覆う。
    for ((point, weight), (month, days)) in
        monthly
            .fixings()
            .iter()
            .zip([(1, 1.0), (2, 28.0), (3, 1.0)])
    {
        assert_eq!(
            *point,
            ObservationPoint::Month(YearMonth::new(2027, month)?)
        );
        assert_close(*weight, days / 365.0);
    }
    assert_eq!(monthly.fixings().len(), 3);
    month_rule = ObservationRule::new(
        ObservationDays::Months,
        ObservationMethod::Plain,
        Aggregation::Average(AveragingWeights::Equal),
        ObservationRounding::default(),
    )?;
    let equal = observe(&dates, &months.periods()[0], &month_rule)?;
    for (_, weight) in equal.fixings() {
        assert_close(*weight, 1.0);
    }
    // 指数比はシフト後の木曜・金曜を使い、年率換算には 1 日分の係数を使う。
    let shifted = ObservationDays::Endpoints(EndpointDates::Shifted {
        calendar: dates.spec(),
        lag: NonZeroU8::MIN,
    });
    let index = ObservationRule::new(
        shifted,
        ObservationMethod::Plain,
        Aggregation::IndexRatio,
        ObservationRounding::default(),
    )?;
    let endpoints = observe(&dates, period, &index)?;
    assert_eq!(
        endpoints
            .fixings()
            .iter()
            .map(|&(point, _)| point)
            .collect::<Vec<_>>(),
        [
            ObservationPoint::Date("2027-01-07".parse()?),
            ObservationPoint::Date("2027-01-08".parse()?),
        ]
    );
    assert_close(endpoints.annualization().value(), 1.0 / 365.0);
    Ok(())
}

#[test]
fn incompatible_observation_rules_are_invalid_data() -> Result<()> {
    let dates = Dates::new(&[])?;
    let rounding = ObservationRounding {
        daily: Some(Rounding::new(1e-7, RoundingMode::HalfAwayFromZero)?),
        aggregated: None,
    };
    let plain = ObservationMethod::Plain;
    let none = ObservationRounding::default();
    // 日の選び方と方法、集約、丸めの許されない組をそれぞれ拒否する。
    for (days, method, aggregation, rounding) in [
        (
            ObservationDays::AllDays {
                fixing_calendar: dates.spec(),
            },
            plain,
            COMPOUND,
            none,
        ),
        (
            ObservationDays::Months,
            ObservationMethod::Lookback(NonZeroU8::MIN),
            WEIGHTED,
            none,
        ),
        (
            ObservationDays::Endpoints(EndpointDates::Accrual),
            plain,
            Aggregation::IndexRatio,
            rounding,
        ),
        (
            ObservationDays::BusinessDays(dates.spec()),
            plain,
            Aggregation::Compound {
                form: CompoundingForm::PowerFactor,
                spread: SpreadTreatment::Inclusive,
            },
            none,
        ),
    ] {
        assert!(matches!(
            ObservationRule::new(days, method, aggregation, rounding),
            Err(Error::InvalidData(_))
        ));
    }
    Ok(())
}

#[test]
fn all_day_lockout_and_accrual_endpoints_follow_the_rule() -> Result<()> {
    // 2027-01-05（火）を休日にし、月曜から金曜までの 5 暦日の期間を観測する。
    let dates = Dates::new(&["2027-01-05"])?;
    let schedule = dates.once("2027-01-04", "2027-01-09")?;
    let period = &schedule.periods()[0];
    let all_days = ObservationDays::AllDays {
        fixing_calendar: dates.spec(),
    };
    let lockout = |lag| {
        ObservationRule::new(
            all_days.clone(),
            ObservationMethod::Lockout(lag),
            Aggregation::Average(AveragingWeights::Equal),
            ObservationRounding::default(),
        )
    };
    let lag = NonZeroU8::new(2).ok_or_else(|| Error::Internal("zero lockout".into()))?;
    let observed = observe(&dates, period, &lockout(lag)?)?;
    // 暦日ごとの行は 月 04・火 05・水 06・木 07・金 08 で、休日の火曜は月曜の値を使う。
    // Lockout(2) は m = 5 の観測で締めの日 d_3（水 06）より後の 2 行を水曜の値にする。
    // 等重みの平均なので重みはすべて 1 で、各行は自分の暦日 1 日を覆う。
    let fixings: Vec<_> = observed.fixings().iter().map(|&(point, _)| point).collect();
    let expected = [
        "2027-01-04",
        "2027-01-04",
        "2027-01-06",
        "2027-01-06",
        "2027-01-06",
    ]
    .into_iter()
    .map(|date| Ok(ObservationPoint::Date(date.parse()?)))
    .collect::<Result<Vec<_>>>()?;
    assert_eq!(fixings, expected);
    for &(_, weight) in observed.fixings() {
        assert_close(weight, 1.0);
    }
    let covered = [
        "2027-01-04",
        "2027-01-05",
        "2027-01-06",
        "2027-01-07",
        "2027-01-08",
    ]
    .into_iter()
    .map(|date| {
        let date: rhodium_core::date::Date = date.parse()?;
        rhodium_core::date::DateInterval::new(date, date.checked_add_days(1)?)
    })
    .collect::<Result<Vec<_>>>()?;
    assert_eq!(observed.coverage(), covered.as_slice());
    assert_close(observed.annualization().value(), 5.0 / 365.0);

    // 締めの日が残らない短い期間（2 行で Lockout(2)）は適用できない。
    let short = dates.once("2027-01-07", "2027-01-09")?;
    assert!(matches!(
        observation_schedule(
            &short.periods()[0],
            &lockout(lag)?,
            &DayCount::Act365F,
            &AccrualContext::default(),
            &dates.calendars,
        ),
        Err(Error::InvalidData(_))
    ));

    // 複利指数の比は計算期間の両端の水準を読み、重みは 1、年率換算は期間全体の 3 日分。
    let schedule = dates.once("2027-01-08", "2027-01-11")?;
    let index = ObservationRule::new(
        ObservationDays::Endpoints(EndpointDates::Accrual),
        ObservationMethod::Plain,
        Aggregation::IndexRatio,
        ObservationRounding::default(),
    )?;
    let endpoints = observe(&dates, &schedule.periods()[0], &index)?;
    assert_eq!(
        endpoints
            .fixings()
            .iter()
            .map(|&(point, _)| point)
            .collect::<Vec<_>>(),
        [
            ObservationPoint::Date("2027-01-08".parse()?),
            ObservationPoint::Date("2027-01-11".parse()?),
        ]
    );
    for &(_, weight) in endpoints.fixings() {
        assert_close(weight, 1.0);
    }
    assert_close(endpoints.annualization().value(), 3.0 / 365.0);
    Ok(())
}

#[test]
fn factor_excess_keeps_the_digits_with_arithmetic_only() -> Result<()> {
    // 因子の超過 factor − 1 は、単利の因子では (r + spread)·w そのもので、1 を足して引く丸めも
    // 超越関数も使わない（確定した利息が環境の libm によらず同じ値になる）。冪の形は exp_m1(w·ln_1p(r))。
    let (rate, weight) = (0.005, 1.0 / 365.0);
    assert_eq!(
        CompoundingForm::SimpleFactor.factor_excess(rate, weight, 0.0)?,
        rate * weight
    );
    assert_eq!(
        CompoundingForm::SimpleFactor.factor_excess(0.01, 0.25, 0.002)?,
        (0.01 + 0.002) * 0.25
    );
    assert_close(
        CompoundingForm::PowerFactor.factor_excess(0.03, 0.25, 0.0)?,
        CompoundingForm::PowerFactor.factor(0.03, 0.25, 0.0)? - 1.0,
    );
    // 因子が正にならない値と非有限の値は、factor と同じく Numerical。
    for (form, rate, weight, spread) in [
        (CompoundingForm::SimpleFactor, -1.0, 1.0, 0.0),
        (CompoundingForm::SimpleFactor, 0.01, 1.0, -2.0),
        (CompoundingForm::SimpleFactor, f64::NAN, 1.0, 0.0),
        (CompoundingForm::PowerFactor, -1.5, 0.5, 0.0),
        (CompoundingForm::PowerFactor, f64::INFINITY, 0.5, 0.0),
    ] {
        assert!(matches!(
            form.factor_excess(rate, weight, spread),
            Err(Error::Numerical(_))
        ));
    }
    Ok(())
}

#[test]
fn log_factor_keeps_the_digits_that_one_plus_x_rounds_away() -> Result<()> {
    // 1 日の TONA の因子 1 + r·w（r = 0.5%、w = 1/365）は 1 との差が約 1.4e-5 で、
    // 因子を丸めると差の相対誤差が ulp(1)/1.4e-5 ≈ 1.6e-11 になる。
    // log_factor から exp_m1 で戻した超過は、r·w を倍精度の精度のまま再現する。
    let (rate, weight) = (0.005, 1.0 / 365.0);
    let increment = rate * weight;
    let log = CompoundingForm::SimpleFactor.log_factor(rate, weight, 0.0)?;
    let excess = log.exp_m1();
    assert!((excess - increment).abs() <= 2.0 * f64::EPSILON * increment);
    // 因子から 1 を引く書き方は、この例で相対 6.8e-12（ε の約 3 万倍）の誤差を持つ。
    let rounded = CompoundingForm::SimpleFactor.factor(rate, weight, 0.0)? - 1.0;
    assert!((rounded - increment).abs() > 1e4 * f64::EPSILON * increment);

    // スプレッドは SimpleFactor の (r + spread)·w に入り、PowerFactor は w·ln(1 + r) になる。
    let spread = CompoundingForm::SimpleFactor.log_factor(0.01, 0.5, 0.002)?;
    assert_close(spread, (0.012_f64 * 0.5).ln_1p());
    let power = CompoundingForm::PowerFactor.log_factor(0.03, 0.25, 0.0)?;
    assert_close(
        power,
        CompoundingForm::PowerFactor.factor(0.03, 0.25, 0.0)?.ln(),
    );

    // 因子が正にならない値と非有限の値は、factor と同じく Numerical。
    for (form, rate, weight, spread) in [
        (CompoundingForm::SimpleFactor, -1.0, 1.0, 0.0),
        (CompoundingForm::SimpleFactor, 0.01, 1.0, -2.0),
        (CompoundingForm::SimpleFactor, f64::NAN, 1.0, 0.0),
        (CompoundingForm::PowerFactor, -1.5, 0.5, 0.0),
        (CompoundingForm::PowerFactor, f64::INFINITY, 0.5, 0.0),
    ] {
        assert!(matches!(
            form.log_factor(rate, weight, spread),
            Err(Error::Numerical(_))
        ));
        assert!(matches!(
            form.factor(rate, weight, spread),
            Err(Error::Numerical(_))
        ));
    }
    Ok(())
}
`````

`````
crates/rhodium-core/tests/time_and_tenor.rs
//! 日付・時刻とテナーの値域、正準形、現地時刻の解決を確認する。

use chrono::{LocalResult, NaiveDateTime, Offset, TimeZone};
use rhodium_core::{
    Error, Result,
    date::{Date, DateInterval, Weekday},
    tenor::{DateOffset, Frequency, Tenor},
    time::{Instant, LocalDateTime, LocalEventTime, LocalTime, UtcOffset, ZoneRules, model_time},
};
use serde::{Deserialize, de::value::StrDeserializer};

#[test]
fn dates_and_times_check_precision_and_bounds() -> Result<()> {
    // Gregorian 暦では 1900 年は平年、2000 年は閏年。年と ISO 文字列の形も検証する。
    let leap = "2000-02-29".parse::<Date>()?;
    assert_eq!(leap.weekday(), Weekday::Tue);
    assert_eq!(leap.checked_add_days(1)?.ymd(), (2000, 3, 1));
    assert_eq!(leap.days_since("2000-03-01".parse()?), -1);
    assert_eq!(leap.to_string(), "2000-02-29");
    for invalid in [
        "1900-02-29",
        "2000-2-29",
        " 2000-02-29",
        "0000-01-01",
        "10000-01-01",
    ] {
        assert!(matches!(
            invalid.parse::<Date>(),
            Err(Error::InvalidInput(_))
        ));
        assert!(
            Date::deserialize(StrDeserializer::<serde::de::value::Error>::new(invalid)).is_err()
        );
    }
    let first: Date = "0001-01-01".parse()?;
    let last: Date = "9999-12-31".parse()?;
    assert!(first.checked_add_days(-1).is_err());
    assert!(last.checked_add_days(1).is_err());
    assert!(DateInterval::new(last, first).is_err());
    assert!(DateInterval::new(first, first)?.is_empty());

    // 明示オフセットを UTC に直して比較する。chrono が読み飛ばす桁でも 0 以外は拒否する。
    let at: Instant = "2026-11-18T17:00:00.123456000+09:00".parse()?;
    assert_eq!(at.to_string(), "2026-11-18T08:00:00.123456Z");
    assert_eq!(Instant::from_unix_micros(at.unix_micros())?, at);
    assert_eq!(
        Instant::from_unix_micros(-1)?.to_string(),
        "1969-12-31T23:59:59.999999Z"
    );
    for invalid in [
        "0001-01-01T00:00:00+00:01",
        "9999-12-31T23:59:59-00:01",
        "2027-01-01T00:00:60Z",
        "2027-01-01T00:00:00.1234560001Z",
    ] {
        assert!(matches!(
            invalid.parse::<Instant>(),
            Err(Error::InvalidInput(_))
        ));
        assert!(
            Instant::deserialize(StrDeserializer::<serde::de::value::Error>::new(invalid)).is_err()
        );
    }
    let local: LocalDateTime = "2026-11-18T17:00:00.123456".parse()?;
    assert_eq!(local.at_offset(UtcOffset::from_seconds(9 * 3_600)?)?, at);
    assert_eq!(local.time(), LocalTime::from_hms_micro(17, 0, 0, 123_456)?);
    assert!("00:00:60".parse::<LocalTime>().is_err());
    assert!("00:00:00.000000001".parse::<LocalTime>().is_err());
    assert!("0000-12-31T23:59:59".parse::<LocalDateTime>().is_err());
    assert!(LocalTime::from_hms_micro(23, 59, 59, 1_000_000).is_err());

    // 半日の実経過時間はモデル時間で 0.5 / 365 年。順序を逆にすると負になる。
    let noon: Instant = "2027-01-01T12:00:00Z".parse()?;
    let midnight: Instant = "2027-01-01T00:00:00Z".parse()?;
    for (actual, expected) in [
        (model_time(midnight, noon).years(), 0.5_f64 / 365.0),
        (model_time(noon, midnight).years(), -0.5 / 365.0),
    ] {
        let scale = actual.abs().max(expected.abs());
        assert!((actual - expected).abs() <= (1e-12 * scale).max(1e-14));
    }
    Ok(())
}

#[test]
fn local_time_text_requires_two_digit_fields() -> Result<()> {
    // chrono の書式解析は数値の前の空白と 1 桁の欄を許す。旧実装と同じく、現地時刻は
    // HH:MM:SS と任意の小数秒の形だけを受け付け、書き出しと往復できない綴りを拒否する。
    for invalid in [
        "1:02:03",
        "01:2:03",
        "01:02:3",
        " 01:02:03",
        "01: 02:03",
        "01:02",
        "01:02:03.",
        "01:02:03Z",
        "24:00:00",
    ] {
        assert!(
            matches!(invalid.parse::<LocalTime>(), Err(Error::InvalidInput(_))),
            "{invalid}"
        );
        assert!(
            LocalTime::deserialize(StrDeserializer::<serde::de::value::Error>::new(invalid))
                .is_err()
        );
    }
    // 日時の時刻部分にも同じ形を求める。
    assert!("2027-01-01T 1:02:03".parse::<LocalDateTime>().is_err());

    // 小数秒の末尾のゼロは 7 桁以上でも読め、書き出した綴りから同じ時刻に戻る。
    let time: LocalTime = "01:02:03.5000000".parse()?;
    assert_eq!(time, LocalTime::from_hms_micro(1, 2, 3, 500_000)?);
    assert_eq!(time.to_string().parse::<LocalTime>()?, time);
    Ok(())
}

#[test]
fn zone_rules_resolve_gaps_and_overlaps() -> Result<()> {
    let winter = UtcOffset::from_seconds(3_600)?;
    let summer = UtcOffset::from_seconds(7_200)?;
    let first: Instant = "2027-01-01T00:00:00Z".parse()?;
    let last: Instant = "2027-12-31T23:59:59Z".parse()?;
    let spring: Instant = "2027-03-28T01:00:00Z".parse()?;
    let autumn: Instant = "2027-10-31T01:00:00Z".parse()?;
    let berlin = ZoneRules::new(
        first..=last,
        winter,
        vec![(spring, summer), (autumn, winter)],
    )?;

    // 春の 02:30 は欠落。慣行で 03:30 に進めると、夏時間の +02:00 なので 01:30 UTC。
    let event = LocalEventTime {
        time: "02:30:00".parse()?,
        zone: "Europe/Berlin".parse()?,
    };
    let gap = event.on_date("2027-03-28".parse()?);
    assert!(matches!(
        berlin.resolve_input(gap, None),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        berlin.resolve_input(gap, Some(winter)),
        Err(Error::InvalidInput(_))
    ));
    assert_eq!(
        berlin.resolve_convention(gap)?,
        "2027-03-28T01:30:00Z".parse()?
    );
    assert_eq!(berlin.to_local(spring)?, "2027-03-28T03:00:00".parse()?);
    assert_eq!(
        berlin.offset_at(Instant::from_unix_micros(spring.unix_micros() - 1)?)?,
        winter
    );
    assert_eq!(berlin.offset_at(spring)?, summer);

    // 秋の 02:30 は +02:00 と +01:00 の二つ。慣行では早い 00:30 UTC を選ぶ。
    let overlap = event.on_date("2027-10-31".parse()?);
    assert!(matches!(
        berlin.resolve_input(overlap, None),
        Err(Error::InvalidInput(_))
    ));
    assert_eq!(
        berlin.resolve_input(overlap, Some(summer))?,
        "2027-10-31T00:30:00Z".parse()?
    );
    assert_eq!(
        berlin.resolve_input(overlap, Some(winter))?,
        "2027-10-31T01:30:00Z".parse()?
    );
    assert_eq!(
        berlin.resolve_convention(overlap)?,
        "2027-10-31T00:30:00Z".parse()?
    );
    let naive = "2027-10-31T02:30:00"
        .parse::<NaiveDateTime>()
        .map_err(|error| Error::Internal(error.to_string()))?;
    let LocalResult::Ambiguous(early, late) = berlin.from_local_datetime(&naive) else {
        return Err(Error::Internal("expected two chrono mappings".into()));
    };
    assert!(early < late);
    assert_eq!(early.offset().fix().local_minus_utc(), 7_200);
    assert_eq!(late.timezone().offset_at(last)?, winter);

    // 範囲の両端は含み、その直前・直後はギャップと異なる MissingData。
    assert_eq!(berlin.offset_at(first)?, winter);
    assert_eq!(berlin.offset_at(last)?, winter);
    let outside = Instant::from_unix_micros(last.unix_micros() + 1)?;
    assert!(matches!(
        berlin.to_local(outside),
        Err(Error::MissingData(_))
    ));
    assert!(matches!(
        berlin.resolve_input("2027-01-01T00:59:59".parse()?, None),
        Err(Error::MissingData(_))
    ));
    assert_eq!(
        berlin.resolve_input("2027-01-01T01:00:00".parse()?, None)?,
        first
    );

    // 遷移から離れた冬の正午に夏時間のオフセットを付けた入力は、規則と照合して拒否する。
    // 正しい冬時間のオフセットなら、その瞬間（11:00 UTC）に解決する。
    let noon = "2027-01-15T12:00:00".parse()?;
    assert!(matches!(
        berlin.resolve_input(noon, Some(summer)),
        Err(Error::InvalidInput(_))
    ));
    assert_eq!(
        berlin.resolve_input(noon, Some(winter))?,
        "2027-01-15T11:00:00Z".parse()?
    );
    assert_eq!(
        berlin.resolve_input(noon, None)?,
        "2027-01-15T11:00:00Z".parse()?
    );

    // 30 分だけ進む規則なら、同じ 02:10 の入力は 02:40 に進んで 01:10 UTC になる。
    let half_hour = ZoneRules::new(
        first..=last,
        winter,
        vec![(spring, UtcOffset::from_seconds(5_400)?)],
    )?;
    assert_eq!(
        half_hour.resolve_convention("2027-03-28T02:10:00".parse()?)?,
        "2027-03-28T01:10:00Z".parse()?
    );
    Ok(())
}

#[test]
fn zone_rules_reject_tables_that_cannot_be_mapped_by_chrono() -> Result<()> {
    let first: Instant = "2027-01-01T00:00:00Z".parse()?;
    let last: Instant = "2027-01-02T00:00:00Z".parse()?;
    let initial = UtcOffset::from_seconds(7_200)?;
    let middle = UtcOffset::from_seconds(3_600)?;
    let final_offset = UtcOffset::from_seconds(0)?;
    let early: Instant = "2027-01-01T01:00:00Z".parse()?;
    let late: Instant = "2027-01-01T01:30:00Z".parse()?;

    // 二度戻ると現地 02:15 に UTC の候補が三つでき、chrono の LocalResult では表せない。
    for transitions in [
        vec![(early, middle), (late, final_offset)],
        vec![(late, middle), (early, final_offset)],
        vec![(early, middle), (early, final_offset)],
        vec![(first, middle)],
    ] {
        assert!(matches!(
            ZoneRules::new(first..=last, initial, transitions),
            Err(Error::InvalidData(_))
        ));
    }
    Ok(())
}

#[test]
fn tenors_and_frequencies_preserve_contract_units() -> Result<()> {
    // 暦の長さは正準化するが、営業日と暦日は同じ長さとみなさない。
    for (text, canonical) in [
        ("14D", "2W"),
        ("12M", "1Y"),
        ("18M", "18M"),
        ("0Y", "0D"),
        ("0BD", "0BD"),
    ] {
        assert_eq!(text.parse::<Tenor>()?.to_string(), canonical);
    }
    assert_eq!(
        "1W".parse::<Tenor>()?.checked_add("1D".parse()?)?,
        "8D".parse()?
    );
    assert_eq!(
        "1Y".parse::<Tenor>()?.checked_add("6M".parse()?)?,
        "18M".parse()?
    );
    assert_eq!(
        "0D".parse::<Tenor>()?.checked_add("1M".parse()?)?,
        "1M".parse()?
    );
    assert_eq!("7D".parse::<Tenor>()?.times(2)?, "2W".parse()?);
    // 0 倍した年は暦のゼロ 0D になり、月・年の族の外の週とも足せる。
    assert_eq!(
        "1Y".parse::<Tenor>()?
            .times(0)?
            .checked_add("2W".parse()?)?,
        "14D".parse()?
    );
    // ON/TN/SN は長さではないので倍にできない。
    assert!(matches!(
        Tenor::Overnight.times(2),
        Err(Error::InvalidInput(_))
    ));
    assert!("65535Y".parse::<Tenor>()?.times(2).is_err());
    for (left, right) in [("30D", "1M"), ("0D", "2BD"), ("ON", "TN")] {
        assert!(matches!(
            left.parse::<Tenor>()?.checked_add(right.parse()?),
            Err(Error::InvalidInput(_))
        ));
    }
    assert!("1D".parse::<Tenor>()? < Tenor::Overnight);
    assert!("1W".parse::<Tenor>()? < "5BD".parse::<Tenor>()?);
    assert_eq!("-32768BD".parse::<DateOffset>()?.amount(), i16::MIN);
    assert!("-32769BD".parse::<DateOffset>().is_err());
    for invalid in ["1y", "-1D", "01D", "1.0M", "65536BD"] {
        assert!(matches!(
            invalid.parse::<Tenor>(),
            Err(Error::InvalidData(_))
        ));
    }
    for (text, periods) in [
        ("1M", 12),
        ("2M", 6),
        ("3M", 4),
        ("4M", 3),
        ("6M", 2),
        ("12M", 1),
    ] {
        assert_eq!(
            text.parse::<Frequency>()?
                .periods_per_year()
                .map(std::num::NonZeroU8::get),
            Some(periods)
        );
    }
    assert_eq!("4W".parse::<Frequency>()?.periods_per_year(), None);
    assert_eq!("ONCE".parse::<Frequency>()?, Frequency::ONCE);
    assert!("0D".parse::<Frequency>().is_err());
    assert!("ON".parse::<Frequency>().is_err());
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/clean_price_settlement.rs
use rhodium_core::{
    Result,
    currency::Currency,
    money::Money,
    quantity::{Clean, ParPrice},
};

fn main() -> Result<()> {
    // 経過利息を含める前の価格では決済額を求められない。
    let clean = ParPrice::<Clean>::new(0.98)?;
    let _ = clean.settlement(Money::new(Currency::JPY, 1000.0));
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/clean_price_settlement_ok.rs
use rhodium_core::{
    Result,
    currency::Currency,
    money::Money,
    quantity::{Clean, ParPrice},
};

fn main() -> Result<()> {
    // 経過利息を含む価格から決済額を求める。
    let clean = ParPrice::<Clean>::new(0.98)?;
    let _ = clean
        .with_accrued(0.03)?
        .settlement(Money::new(Currency::JPY, 1000.0))?;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/date_as_instant.rs
use rhodium_core::{
    date::Date,
    time::{Instant, model_time},
};

fn wrong(origin: Instant, at: Date) {
    let _ = model_time(origin, at);
}

fn main() {}
`````

`````
crates/rhodium-core/tests/ui/date_as_instant_ok.rs
use rhodium_core::{
    Result,
    time::{Instant, model_time},
};

fn main() -> Result<()> {
    let origin = "2027-01-01T00:00:00Z".parse::<Instant>()?;
    let at = "2027-01-02T00:00:00Z".parse::<Instant>()?;
    assert!(model_time(origin, at).years().is_sign_positive());
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/date_plus_tenor.rs
use rhodium_core::{date::Date, tenor::Tenor};

fn wrong(date: Date, tenor: Tenor) {
    let _ = date + tenor;
}

fn main() {}
`````

`````
crates/rhodium-core/tests/ui/date_plus_tenor_ok.rs
use rhodium_core::{Result, date::Date};

fn main() -> Result<()> {
    let date = Date::from_ymd(2027, 1, 1)?;
    assert_eq!(date.checked_add_days(1)?.ymd(), (2027, 1, 2));
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/futures_price_as_rate.rs
use rhodium_core::{
    Result,
    day_count::DayCount::Act360,
    quantity::{RateFuturesPrice, YearFraction, accrue},
};

fn main() -> Result<()> {
    // 先物価格はそのまま利息計算用の率として使えない。
    let price = RateFuturesPrice::new(0.96, Act360)?;
    let yf = YearFraction::new(0.5, Act360)?;
    let _ = accrue(&price, &yf);
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/futures_price_as_rate_ok.rs
use rhodium_core::{
    Result,
    day_count::DayCount::Act360,
    quantity::{RateFuturesPrice, YearFraction, accrue},
};

fn main() -> Result<()> {
    // 先物価格は implied_rate で率に変換してから使う。
    let price = RateFuturesPrice::new(0.96, Act360)?;
    let yf = YearFraction::new(0.5, Act360)?;
    let _ = accrue(&price.implied_rate()?, &yf)?;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/fx_spot_as_today.rs
use rhodium_core::{Result, currency::Currency, date::Date, fx::FxRate, money::Money};

fn main() -> Result<()> {
    // スポット受渡の FX レートでは現在価値を換算できない。
    let spot = FxRate::spot("USDJPY".parse()?, 150.0, Date::from_ymd(2026, 11, 24)?)?;
    let _ = spot.convert(Money::new(Currency::USD, 2.0));
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/fx_spot_as_today_ok.rs
use rhodium_core::{Result, currency::Currency, fx::FxRate, money::Money};

fn main() -> Result<()> {
    // 評価時点に受け渡す FX レートで現在価値を換算する。
    let today = FxRate::today("USDJPY".parse()?, 150.0)?;
    let _ = today.convert(Money::new(Currency::USD, 2.0))?;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/normal_plus_lognormal.rs
use rhodium_core::{
    Result,
    quantity::{LognormalVol, NormalVol},
};

fn main() -> Result<()> {
    // 原資産の単位を持つ正規ボラと対数正規ボラは比較できない。
    let normal = NormalVol::new(0.01)?;
    let lognormal = LognormalVol::new(0.20)?;
    let _ = normal < lognormal;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/normal_plus_lognormal_ok.rs
use rhodium_core::{Result, quantity::NormalVol};

fn main() -> Result<()> {
    // 正規ボラ同士なら比較できる。
    let _ = NormalVol::new(0.01)? < NormalVol::new(0.02)?;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/rate_kind_mismatch.rs
use rhodium_core::{
    Result,
    day_count::DayCount::Act360,
    quantity::{Continuous, Rate, Simple},
};

fn main() -> Result<()> {
    // 複利型が違う率は足せない。
    let simple = Rate::new(0.01, Act360, Simple)?;
    let continuous = Rate::new(0.02, Act360, Continuous)?;
    let _ = simple + continuous;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/rate_kind_mismatch_ok.rs
use rhodium_core::{
    Result,
    day_count::DayCount::Act360,
    quantity::{Rate, Simple},
};

fn main() -> Result<()> {
    // 同じ日数基準と複利型の率は足せる。
    let a = Rate::new(0.01, Act360, Simple)?;
    let b = Rate::new(0.02, Act360, Simple)?;
    let _ = (a + b)?;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/rate_times_year_fraction.rs
use rhodium_core::{
    Result,
    day_count::DayCount::Act360,
    quantity::{Rate, Simple, YearFraction},
};

fn main() -> Result<()> {
    // 率と日数係数は直接掛けられない。
    let rate = Rate::new(0.01, Act360, Simple)?;
    let yf = YearFraction::new(0.5, Act360)?;
    let _ = rate * yf;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/rate_times_year_fraction_ok.rs
use rhodium_core::{
    Result,
    day_count::DayCount::Act360,
    quantity::{Rate, Simple, YearFraction, accrue},
};

fn main() -> Result<()> {
    // accrue で日数基準を照合してから利息を求める。
    let rate = Rate::new(0.01, Act360, Simple)?;
    let yf = YearFraction::new(0.5, Act360)?;
    let _ = accrue(&rate, &yf)?;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/survival_as_discount.rs
use rhodium_core::{
    Result,
    quantity::{DiscountFactor, SurvivalProbability},
};

fn main() -> Result<()> {
    // 区間の成長因子 P(s)/P(e) は割引因子どうしからだけ作れ、生存確率は渡せない。
    let start = DiscountFactor::new(0.95)?;
    let survival = SurvivalProbability::new(0.9)?;
    let _ = start.growth_to(survival)?;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/survival_as_discount_ok.rs
use rhodium_core::{Result, quantity::DiscountFactor};

fn main() -> Result<()> {
    // 区間の成長因子 P(s)/P(e) は二つの割引因子から作る。
    let start = DiscountFactor::new(0.95)?;
    let end = DiscountFactor::new(0.9)?;
    let _ = start.growth_to(end)?;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/upfront_plus_running.rs
use rhodium_core::{
    Result,
    quantity::{RunningSpread, UpfrontFraction},
};

fn main() -> Result<()> {
    // アップフロントとランニング・スプレッドは足せない。
    let _ = UpfrontFraction::new(0.02)? + RunningSpread::new(0.01)?;
    Ok(())
}
`````

`````
crates/rhodium-core/tests/ui/upfront_plus_running_ok.rs
use rhodium_core::{
    Result,
    quantity::{RiskyAnnuity, RunningSpread},
};

fn main() -> Result<()> {
    // risky PV01 を使ってアップフロントに換算する。
    let _ = RunningSpread::new(0.01)?.to_upfront(RiskyAnnuity::new(2.0)?)?;
    Ok(())
}
`````

`````
crates/rhodium-market/Cargo.toml
[package]
name = "rhodium-market"
description = "市場の構造、quote、曲線とボラ面"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
publish.workspace = true

[lints]
workspace = true

[dependencies]
rhodium-core.workspace = true
scioxide-ad.workspace = true
serde.workspace = true

[dev-dependencies]
rhodium-testkit.workspace = true
scioxide-core.workspace = true
scioxide-interp.workspace = true
`````

`````
crates/rhodium-market/src/curve.rs
//! イールドカーブの識別子、役割、節点。

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
    str::FromStr,
    sync::Arc,
};

use rhodium_core::{
    Error, Result,
    currency::Currency,
    fx::CcyPair,
    time::{Instant, model_time},
};
use scioxide_ad::Real;
use serde::{Deserialize, Serialize};

use crate::{
    market::{Market, ThetaBlock},
    structure::{BlockGrid, BlockId, MarketStructure, ResolvedExpr},
};

/// 曲線の識別子。文字列から通貨や曲線の役割を読み取ることはない。
///
/// 1〜128 バイトで、英字、数字、`-`、`_`、`.`、`:`、`/`、`+` だけを使える
/// （例: `JPY-TONA`、`USD/csa:JPY-TONA`）。大小や区切りは正規化しない。
/// [`FromStr`] と serde は同じ文法で検査し、違反は `InvalidInput`。
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(try_from = "String", into = "String")]
pub struct CurveId(Arc<str>);

impl CurveId {
    /// ID の文字列。
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl FromStr for CurveId {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        if !(1..=128).contains(&text.len())
            || !text
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"-_.:/+".contains(&b))
        {
            return Err(Error::InvalidInput(format!("invalid curve ID: {text:?}")));
        }
        Ok(Self(text.into()))
    }
}

impl TryFrom<String> for CurveId {
    type Error = Error;

    fn try_from(value: String) -> Result<Self> {
        value.parse()
    }
}

impl From<CurveId> for String {
    fn from(id: CurveId) -> Self {
        id.as_str().to_owned()
    }
}

impl fmt::Display for CurveId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// イールドカーブが担う役割。一つの曲線に複数の役割を割り当てられる。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CurveRole {
    /// 支払いの現在価値を求める。
    Discount,
    /// 指標の将来の金利を予想する。
    Projection,
}

/// 最後のピラーより後の外挿の規則。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Extrapolation {
    /// 最後の区間の瞬間フォワードを保つ。
    #[default]
    FlatInstantaneousForward,
    /// 原点から最後のピラーまでの平均フォワードを保つ。
    FlatAverageForward,
    /// 最後のピラーより後の点を拒否する。
    Reject,
}

/// 割引因子または予想用の擬似割引因子の対数を作る部品。
///
/// ファイルの形式はこの型の serde の表現そのもので、瞬間は UTC の日時、`Flat` の率は
/// 接尾辞付きの文字列（`"25bp"`、`"0.25%"`）で書く。外貨担保曲線（`FxImplied`）は市場の定義の
/// FX から作るので、ファイルの曲線には書けない。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum CurveExpr {
    /// 原点を固定し、ln DF を線形補間する節点。
    Nodes {
        /// θ の部品名。省略すると `nodes`。
        #[serde(default = "nodes_part")]
        part: String,
        /// UTC のピラー。
        pillars: Vec<Instant>,
        /// 最後のピラー以後の規則。
        #[serde(default)]
        extrapolation: Extrapolation,
    },
    /// 同じ通貨の別の曲線の ln DF。
    Ref(CurveId),
    /// 指標のフォワードの積分に負号を付けた ln DF。付利の補正などでは別通貨の指標も参照できる。
    /// 参照先には [`CurveRole::Projection`] が必要。
    ForwardRef(CurveId),
    /// 各部品の ln DF の和。
    Sum(Vec<Self>),
    /// ln DF の定数倍。
    Scale {
        /// 有限な倍率。負の値とゼロも使える。
        factor: f64,
        /// 倍率を掛ける部品。
        body: Box<Self>,
    },
    /// 一定の瞬間フォワードで得る ln DF。
    Flat(#[serde(deserialize_with = "scaled_rate")] f64),
    /// 格子の区間ごとに候補の平均フォワードの最大値を取り、区間内ではその値を一定とする。
    /// 同点の候補には感応度を等分し、格子の最終時刻より先の点は拒否する。
    MaxForward {
        /// 一つ以上の候補。
        parts: Vec<Self>,
        /// 区間の右端となる UTC 時刻。評価時点は暗黙の原点で、ここには含めない。
        grid: Vec<Instant>,
    },
    /// 短期と長期を境界で連続につなぐ。
    Stitch {
        /// 接続点までの部品。
        short: Box<Self>,
        /// 接続点以後の傾きを与える部品。
        long: Box<Self>,
        /// 接続点。
        at: Instant,
    },
    /// 各会合で瞬間フォワードが切り替わる曲線。
    Steps {
        /// θ の部品名。要素は最初の会合前から始まる区間ごとの瞬間フォワード。
        part: String,
        /// 会合の発効時刻。
        meetings: Vec<Instant>,
    },
    /// 指定した瞬間から ln DF に負の跳びを加える部品。
    Jumps {
        /// θ の部品名。要素は各受渡期間のスプレッドと日数係数の積。
        part: String,
        /// 跳びが発効する瞬間。
        dates: Vec<Instant>,
    },
    /// 同じ担保の反対通貨の割引曲線と FX の対数フォワード比から作る曲線。
    #[serde(skip_deserializing)]
    FxImplied {
        /// θ の部品名。要素は 基準通貨と決済通貨の ln DF の差。
        part: String,
        /// 市場の標準の向き。
        pair: CcyPair,
        /// 反対通貨の割引曲線。
        base: CurveId,
        /// FX の対数フォワード比の UTC ピラー。
        pillars: Vec<Instant>,
        /// 最後のピラー以後の規則。
        extrapolation: Extrapolation,
    },
}

fn nodes_part() -> String {
    "nodes".into()
}

/// 率を接尾辞付きの文字列（`"25bp"`、`"0.25%"`）で読む。裸の数は単位の書き忘れとして拒否する。
fn scaled_rate<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> std::result::Result<f64, D::Error> {
    let text = <String as serde::Deserialize>::deserialize(deserializer)?;
    // 十進のまま尺度を変えてから一度だけ f64 に丸め、単位の違う等価な値を同じ値にする。
    let (number, scale) = text
        .strip_suffix("bp")
        .map(|number| (number, 4))
        .or_else(|| text.strip_suffix('%').map(|number| (number, 2)))
        .ok_or_else(|| serde::de::Error::custom(format!("rate {text:?} requires bp or %")))?;
    format!("{}e-{scale}", number.trim())
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .ok_or_else(|| serde::de::Error::custom(format!("invalid rate {text:?}")))
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum ParameterSpec<'a> {
    Nodes(&'a [Instant]),
    Steps(&'a [Instant]),
    Jumps(&'a [Instant]),
    FxImplied(&'a [Instant]),
}

impl ParameterSpec<'_> {
    pub(crate) fn len(self) -> usize {
        match self {
            Self::Nodes(dates) | Self::Jumps(dates) | Self::FxImplied(dates) => dates.len(),
            Self::Steps(dates) => dates.len() + 1,
        }
    }

    pub(crate) fn dates(&self) -> &[Instant] {
        match self {
            Self::Nodes(dates)
            | Self::Steps(dates)
            | Self::Jumps(dates)
            | Self::FxImplied(dates) => dates,
        }
    }
}

impl CurveExpr {
    fn visit_blocks<'a>(
        &'a self,
        visit: &mut impl FnMut(&'a str, ParameterSpec<'a>) -> Result<()>,
    ) -> Result<()> {
        match self {
            Self::Nodes { part, pillars, .. } => visit(part, ParameterSpec::Nodes(pillars)),
            Self::Steps { part, meetings } => visit(part, ParameterSpec::Steps(meetings)),
            Self::Jumps { part, dates } => visit(part, ParameterSpec::Jumps(dates)),
            Self::FxImplied { part, pillars, .. } => visit(part, ParameterSpec::FxImplied(pillars)),
            Self::Sum(parts) if !parts.is_empty() => {
                parts.iter().try_for_each(|part| part.visit_blocks(visit))
            }
            Self::Scale { factor, body } if factor.is_finite() => body.visit_blocks(visit),
            Self::MaxForward { parts, grid }
                if !parts.is_empty() && !grid.is_empty() && grid.is_sorted_by(|a, b| a < b) =>
            {
                parts.iter().try_for_each(|part| part.visit_blocks(visit))
            }
            Self::Stitch { short, long, .. } => {
                short.visit_blocks(visit)?;
                long.visit_blocks(visit)
            }
            Self::Flat(rate) if rate.is_finite() => Ok(()),
            Self::Ref(_) | Self::ForwardRef(_) => Ok(()),
            Self::Flat(_) | Self::Sum(_) => Err(Error::InvalidInput(
                "a curve sum must be nonempty and flat rates must be finite".into(),
            )),
            Self::Scale { .. } => Err(Error::InvalidInput(
                "a curve scale factor must be finite".into(),
            )),
            Self::MaxForward { .. } => Err(Error::InvalidInput(
                "max forward requires candidates and strictly increasing grid dates".into(),
            )),
        }
    }
}

/// 通貨と役割を持つイールドカーブの定義。
#[derive(Debug, Clone)]
pub struct CurveSpec {
    id: CurveId,
    currency: Currency,
    roles: BTreeSet<CurveRole>,
    pub(crate) expr: CurveExpr,
}

impl CurveSpec {
    /// 指定した通貨、役割、ピラー、外挿規則で曲線を定義する。
    ///
    /// # Errors
    ///
    /// 役割またはピラーが無い場合か、ピラーが厳密な昇順でなければ `InvalidInput`。
    pub fn new(
        id: CurveId,
        currency: Currency,
        roles: impl IntoIterator<Item = CurveRole>,
        pillars: impl IntoIterator<Item = Instant>,
        extrapolation: Extrapolation,
    ) -> Result<Self> {
        Self::from_expr(
            id,
            currency,
            roles,
            CurveExpr::Nodes {
                part: "nodes".into(),
                pillars: pillars.into_iter().collect(),
                extrapolation,
            },
        )
    }

    /// ln DF の式から曲線を定義する。
    ///
    /// # Errors
    ///
    /// 無効な部品、重複する部品名、昇順でない日付列、役割が無い場合は `InvalidInput`。
    pub fn from_expr(
        id: CurveId,
        currency: Currency,
        roles: impl IntoIterator<Item = CurveRole>,
        expr: CurveExpr,
    ) -> Result<Self> {
        let roles: BTreeSet<_> = roles.into_iter().collect();
        if roles.is_empty() {
            return Err(Error::InvalidInput(format!("curve {id} has no roles")));
        }
        let spec = Self {
            id,
            currency,
            roles,
            expr,
        };
        spec.collect_blocks()?;
        Ok(spec)
    }

    pub(crate) fn collect_blocks(&self) -> Result<BTreeMap<BlockId, ParameterSpec<'_>>> {
        let mut blocks = BTreeMap::new();
        self.expr.visit_blocks(&mut |part, spec| {
            let block = BlockId::new(self.id.clone(), part)?;
            if spec.len() == 0
                || !spec.dates().is_sorted_by(|a, b| a < b)
                || u32::try_from(spec.len() - 1).is_err()
            {
                return Err(Error::InvalidInput(format!(
                    "curve {} needs strictly increasing pillars or event dates for {part}",
                    self.id
                )));
            }
            if blocks.insert(block, spec).is_some() {
                return Err(Error::InvalidInput(format!(
                    "duplicate curve part {} {part}",
                    self.id
                )));
            }
            Ok(())
        })?;
        Ok(blocks)
    }

    /// この曲線が定義する会合日の階段（`Steps`）のブロック。参照先のブロックは含まない。
    pub fn step_blocks(&self) -> impl Iterator<Item = BlockId> {
        let mut blocks = Vec::new();
        // 構築時に部品名と日付を検査済み。
        let _: Result<()> = self.expr.visit_blocks(&mut |part, spec| {
            if matches!(spec, ParameterSpec::Steps(_)) {
                blocks.push(BlockId::for_part(self.id.clone(), part));
            }
            Ok(())
        });
        blocks.into_iter()
    }

    /// この曲線が定義する θ のブロックと要素数。参照先のブロックは含まない。
    pub fn parameter_blocks(&self) -> impl Iterator<Item = (BlockId, usize)> {
        let mut blocks = Vec::new();
        // 構築時に部品名と日付を検査済み。
        let _: Result<()> = self.expr.visit_blocks(&mut |part, spec| {
            blocks.push((BlockId::for_part(self.id.clone(), part), spec.len()));
            Ok(())
        });
        blocks.into_iter()
    }

    /// 曲線を作る部品の式。
    #[must_use]
    pub const fn expr(&self) -> &CurveExpr {
        &self.expr
    }

    /// 曲線の ID。
    #[must_use]
    pub const fn id(&self) -> &CurveId {
        &self.id
    }

    /// 曲線が表す通貨。
    #[must_use]
    pub const fn currency(&self) -> Currency {
        self.currency
    }

    /// 曲線が指定した役割を持つか。
    #[must_use]
    pub fn has_role(&self, role: CurveRole) -> bool {
        self.roles.contains(&role)
    }

    /// 単一の `Nodes` 曲線の UTC ピラー。合成曲線では空。
    #[must_use]
    pub fn pillars(&self) -> &[Instant] {
        match &self.expr {
            CurveExpr::Nodes { part, pillars, .. } if part == "nodes" => pillars,
            _ => &[],
        }
    }

    /// 単一の `Nodes` 曲線の外挿規則。合成曲線では部品ごとに指定する。
    #[must_use]
    pub const fn extrapolation(&self) -> Extrapolation {
        match &self.expr {
            CurveExpr::Nodes { extrapolation, .. } => *extrapolation,
            _ => Extrapolation::FlatInstantaneousForward,
        }
    }
}

/// 曲線と評価時点を指定して準備した点。
///
/// 区間と重みは構造だけで決まり、同じ構造の別の θ でも使える。
#[derive(Debug, Clone)]
pub struct CurvePoint {
    structure: Arc<MarketStructure>,
    curve: usize,
    at: Instant,
    weights: PointWeights,
}

#[derive(Debug, Clone)]
enum PointWeights {
    // 評価時点そのもの。L = 0 で θ を読まない。
    Origin,
    // 単一の Nodes 曲線。確保なしで一つの補間の重みを持つ。
    Nodes(WeightedNodes),
    // θ について線形な式。L = constant + Σ factor · (補間した節点)。
    Composite {
        constant: f64,
        terms: Box<[WeightedNodes]>,
    },
    // MaxForward を含み θ について非線形な式。どの候補が選ばれるかは評価時の θ で決まるので、
    // blocks には全候補が読むブロックを重複なしで持ち、利用可否の検査と依存の収集に使う。
    Nonlinear {
        expr: Box<PreparedExpr>,
        blocks: Box<[usize]>,
    },
}

// 区間 [start, end] の L の変化 ΔL = constant + Σ terms + Σ maxima を準備したもの。
// 評価時点から始まる区間なら、そのまま点の L になる。
#[derive(Debug, Clone, Default)]
struct PreparedExpr {
    constant: f64,
    terms: Vec<WeightedNodes>,
    maxima: Vec<MaxForwardInterval>,
}

// MaxForward の格子の一区間 [t_i, t_{i+1}] の寄与 factor · min_k ΔL_k。
// 各候補の ΔL_k は区間全体で準備し、読み出す点が区間の途中でも区間の平均フォワードを使う。
#[derive(Debug, Clone)]
struct MaxForwardInterval {
    // 外側の倍率に、準備する区間と格子の区間の重なりの割合 overlap / width を掛けた値。
    factor: f64,
    // 候補を first と rest に分け、候補が一つ以上あることを型で保つ。
    first: PreparedExpr,
    rest: Box<[PreparedExpr]>,
}

#[derive(Debug, Clone, Copy)]
struct WeightedNodes {
    block: usize,
    weights: NodeWeights,
    factor: f64,
}

// ln DF の L(t) を、曲線の節点 θ_i（ピラー t_i の L）の疎な一次結合として表す。
// 原点は L(0) = 0 の暗黙の定数なので、θ の添字を持たない。
#[derive(Debug, Clone, Copy)]
enum NodeWeights {
    // ピラー上: L = θ_node。
    Node(usize),
    // 原点を通る直線: L = weight · θ_node。
    // 最初のピラーより前（weight = t / t_0）と、原点からの平均フォワードを保つ外挿
    // （weight = t / t_n）の両方を表す。
    One {
        node: usize,
        weight: f64,
    },
    // ピラー間の線形補間（区間内のフォワードが一定）:
    // L = w_l · θ_left + w_r · θ_right、w_l = (t_r − t) / (t_r − t_l)、w_r = (t − t_l) / (t_r − t_l)。
    Two {
        left: usize,
        right: usize,
        weights: [f64; 2],
    },
    // 最終区間の瞬間フォワードを保つ外挿:
    // L = θ_right + (θ_right − θ_left) · extra、extra = (t − t_n) / (t_n − t_{n−1})。
    LastForward {
        left: usize,
        right: usize,
        extra: f64,
    },
}

impl CurvePoint {
    pub(crate) fn check(&self, structure: &Arc<MarketStructure>, curve: usize) -> Result<()> {
        if !Arc::ptr_eq(structure, &self.structure) || self.curve != curve {
            let expected = structure.curves[curve].spec.id();
            let actual = self.structure.curves[self.curve].spec.id();
            return Err(Error::InvalidInput(format!(
                "point on curve {actual} does not belong to curve {expected} in this structure"
            )));
        }
        Ok(())
    }

    pub(crate) const fn at(&self) -> Instant {
        self.at
    }

    /// この点の計算が読む θ ブロック。複数の部品で読むブロックは複数回現れうる。
    pub fn blocks(&self) -> impl Iterator<Item = &BlockId> {
        self.block_positions()
            .map(|block| &self.structure.layout().ids[block])
    }

    pub(crate) fn block_positions(&self) -> impl Iterator<Item = usize> + '_ {
        let (terms, blocks): (&[WeightedNodes], &[usize]) = match &self.weights {
            PointWeights::Origin => (&[], &[]),
            PointWeights::Nodes(term) => (std::slice::from_ref(term), &[]),
            PointWeights::Composite { terms, .. } => (terms, &[]),
            PointWeights::Nonlinear { blocks, .. } => (&[], blocks),
        };
        terms
            .iter()
            .map(|term| term.block)
            .chain(blocks.iter().copied())
    }

    pub(crate) fn log_df<R: Real>(&self, market: &Market<R>) -> R {
        match &self.weights {
            PointWeights::Origin => R::zero(),
            PointWeights::Nodes(node) => node.weights.eval(market.theta_block(node.block)),
            PointWeights::Composite { constant, terms } => eval_linear(*constant, terms, market),
            PointWeights::Nonlinear { expr, .. } => expr.eval(market),
        }
    }
}

fn eval_linear<R: Real>(constant: f64, terms: &[WeightedNodes], market: &Market<R>) -> R {
    terms.iter().fold(R::constant(constant), |value, term| {
        value + term.weights.eval(market.theta_block(term.block)) * term.factor
    })
}

impl PreparedExpr {
    fn collect_blocks(&self, blocks: &mut BTreeSet<usize>) {
        blocks.extend(self.terms.iter().map(|term| term.block));
        for interval in &self.maxima {
            interval.first.collect_blocks(blocks);
            for candidate in &interval.rest {
                candidate.collect_blocks(blocks);
            }
        }
    }

    fn eval<R: Real>(&self, market: &Market<R>) -> R {
        self.maxima.iter().fold(
            eval_linear(self.constant, &self.terms, market),
            |value, interval| value + interval.min_change(market) * interval.factor,
        )
    }
}

impl MaxForwardInterval {
    // 区間幅は候補で共通なので、平均フォワードの最大は ΔL の最小になる。
    // 同点の m 個には微分を 1/m ずつ配る。どれかの候補が有限でなければ NaN を返し、
    // 読み出しの側で Numerical にする。
    fn min_change<R: Real>(&self, market: &Market<R>) -> R {
        let first = self.first.eval(market);
        let mut valid = first.is_finite();
        let minimum = match self.rest.as_ref() {
            [] => first,
            // 2 個なら scioxide の min が同点に 1/2 ずつ配る（中点の規約）。
            [second] => {
                let second = second.eval(market);
                valid = valid & second.is_finite();
                first.min(second)
            }
            // 2 個ずつの min を重ねると同点の配分が 1/2、1/4、… に偏るので、
            // 同点の候補の微分の平均 mean と個数 count を走査しながら更新する。
            rest => {
                let mut best = first;
                let mut mean = first;
                let mut count = R::one();
                for part in rest {
                    let candidate = part.eval(market);
                    valid = valid & candidate.is_finite();
                    let lower = candidate.lt(best);
                    let tied = candidate.le(best) & candidate.ge(best);
                    let next_count = count + 1.0;
                    // 同点のときだけ新しい候補の微分を平均に足す。同点以外の値同士は引かない。
                    let tied_delta = R::select(tied, candidate, mean) - mean;
                    mean = R::select(lower, candidate, mean + tied_delta / next_count);
                    count = R::select(lower, R::one(), R::select(tied, next_count, count));
                    best = R::select(lower, candidate, best);
                }
                mean
            }
        };
        R::select(valid, minimum, R::constant(f64::NAN))
    }
}

impl NodeWeights {
    fn eval<R: Real>(self, nodes: ThetaBlock<'_, R>) -> R {
        match self {
            NodeWeights::Node(node) => nodes.at(node),
            NodeWeights::One { node, weight } => nodes.at(node) * weight,
            NodeWeights::Two {
                left,
                right,
                weights: [lw, rw],
            } => nodes.at(left) * lw + nodes.at(right) * rw,
            // 最後の傾きを先に求め、近いピラーを長く外挿するときの桁落ちを抑える。
            NodeWeights::LastForward { left, right, extra } => {
                let right = nodes.at(right);
                right + (right - nodes.at(left)) * extra
            }
        }
    }
}

impl MarketStructure {
    /// 曲線の節点を読み出すための点を、評価時点からのモデル時間で準備する。
    ///
    /// # Errors
    ///
    /// 曲線が無ければ `MissingData`。時刻が評価時点より前、または外挿が許されない
    /// 範囲（`MaxForward` の最終時刻より後など）なら `InvalidInput`。
    pub fn prepare_point(self: &Arc<Self>, id: &CurveId, at: Instant) -> Result<CurvePoint> {
        let curve = self
            .layout()
            .position(id)
            .ok_or_else(|| Error::MissingData(format!("unknown curve: {id}")))?;
        self.prepare_point_at(curve, at)
    }

    pub(crate) fn prepare_point_at(
        self: &Arc<Self>,
        curve: usize,
        at: Instant,
    ) -> Result<CurvePoint> {
        let definition = &self.curves[curve];
        let id = definition.spec.id();
        if at < self.as_of() {
            return Err(Error::InvalidInput(format!(
                "point before as-of on curve {id}"
            )));
        }

        let weights = if at == self.as_of() {
            PointWeights::Origin
        } else if let ResolvedExpr::Nodes {
            block,
            extrapolation,
        } = *definition.expr
        {
            // 単一の節点補間は確保なしで保持する。
            PointWeights::Nodes(self.prepare_nodes(block, extrapolation, at)?)
        } else {
            let mut expr = PreparedExpr::default();
            self.prepare_change(&definition.expr, self.as_of(), at, 1.0, &mut expr)?;
            if expr.maxima.is_empty() {
                PointWeights::Composite {
                    constant: expr.constant,
                    terms: expr.terms.into_boxed_slice(),
                }
            } else {
                let mut blocks = BTreeSet::new();
                expr.collect_blocks(&mut blocks);
                PointWeights::Nonlinear {
                    expr: Box::new(expr),
                    blocks: blocks.into_iter().collect(),
                }
            }
        };
        Ok(CurvePoint {
            structure: Arc::clone(self),
            curve,
            at,
            weights,
        })
    }

    fn prepare_nodes(
        &self,
        block: usize,
        extrapolation: Extrapolation,
        at: Instant,
    ) -> Result<WeightedNodes> {
        let BlockGrid { dates, times } = &self.grids[block];
        let last = times.len() - 1;
        let t = model_time(self.as_of(), at).years();
        // UTC で範囲を判定し、モデル時間の二分探索で補間の重みを決める。
        let weights = if at > dates[last] {
            let end = times[last].years();
            match (extrapolation, times[..last].last()) {
                (Extrapolation::Reject, _) => {
                    let id = &self.layout().ids[block];
                    return Err(Error::InvalidInput(format!(
                        "point beyond last pillar on curve {} part {}",
                        id.curve(),
                        id.part()
                    )));
                }
                (Extrapolation::FlatAverageForward, _)
                | (Extrapolation::FlatInstantaneousForward, None) => NodeWeights::One {
                    node: last,
                    weight: t / end,
                },
                (Extrapolation::FlatInstantaneousForward, Some(previous)) => {
                    let extra = (t - end) / (end - previous.years());
                    NodeWeights::LastForward {
                        left: last - 1,
                        right: last,
                        extra,
                    }
                }
            }
        } else {
            match times.binary_search_by(|node| node.years().total_cmp(&t)) {
                Ok(node) => NodeWeights::Node(node),
                Err(0) => NodeWeights::One {
                    node: 0,
                    weight: t / times[0].years(),
                },
                Err(right) => {
                    let left = right - 1;
                    let start = times[left].years();
                    let end = times[right].years();
                    NodeWeights::Two {
                        left,
                        right,
                        weights: [(end - t) / (end - start), (t - start) / (end - start)],
                    }
                }
            }
        };
        Ok(Self::weighted_nodes(block, weights, 1.0))
    }

    fn weighted_nodes(block: usize, weights: NodeWeights, factor: f64) -> WeightedNodes {
        WeightedNodes {
            block,
            weights,
            factor,
        }
    }

    // 区間内の L の変化を準備する。入れ子の MaxForward もこの区間だけを積分する。
    fn prepare_change(
        &self,
        expr: &ResolvedExpr,
        start: Instant,
        end: Instant,
        factor: f64,
        prepared: &mut PreparedExpr,
    ) -> Result<()> {
        if start == end {
            return Ok(());
        }
        match expr {
            ResolvedExpr::Nodes {
                block,
                extrapolation,
            } => {
                let mut last = self.prepare_nodes(*block, *extrapolation, end)?;
                last.factor = factor;
                prepared.terms.push(last);
                if start != self.as_of() {
                    let mut first = self.prepare_nodes(*block, *extrapolation, start)?;
                    first.factor = -factor;
                    prepared.terms.push(first);
                }
            }
            ResolvedExpr::Flat(rate) => {
                prepared.constant -= (rate * model_time(start, end).years()) * factor;
            }
            ResolvedExpr::Sum(parts) => {
                for part in parts {
                    self.prepare_change(part, start, end, factor, prepared)?;
                }
            }
            ResolvedExpr::Scale(scale, body) => {
                self.prepare_change(body, start, end, factor * scale, prepared)?;
            }
            ResolvedExpr::Stitch {
                short,
                long,
                at: boundary,
            } => {
                self.prepare_change(
                    short,
                    start.min(*boundary),
                    end.min(*boundary),
                    factor,
                    prepared,
                )?;
                self.prepare_change(
                    long,
                    start.max(*boundary),
                    end.max(*boundary),
                    factor,
                    prepared,
                )?;
            }
            ResolvedExpr::Steps(block) => {
                let dates = &self.grids[*block].dates;
                let first = dates.partition_point(|&date| date <= start);
                let mut left = start;
                // 最初に重なる区間から、終了時刻までのフォワードだけを積分する。
                for (node, right) in dates[first..]
                    .iter()
                    .copied()
                    .chain(std::iter::once(end))
                    .enumerate()
                {
                    let right = right.min(end);
                    prepared.terms.push(Self::weighted_nodes(
                        *block,
                        NodeWeights::One {
                            node: first + node,
                            weight: model_time(left, right).years(),
                        },
                        -factor,
                    ));
                    if right == end {
                        break;
                    }
                    left = right;
                }
            }
            ResolvedExpr::Jumps(block) => {
                let dates = &self.grids[*block].dates;
                let first = dates.partition_point(|&date| date <= start);
                let last = dates.partition_point(|&date| date <= end);
                prepared
                    .terms
                    .extend((first..last).map(|node| {
                        Self::weighted_nodes(*block, NodeWeights::Node(node), -factor)
                    }));
            }
            ResolvedExpr::MaxForward { parts, grid } => {
                self.prepare_max_forward(parts, grid, start, end, factor, &mut prepared.maxima)?;
            }
        }
        Ok(())
    }

    fn prepare_max_forward(
        &self,
        parts: &[Arc<ResolvedExpr>],
        grid: &[Instant],
        start: Instant,
        end: Instant,
        factor: f64,
        intervals: &mut Vec<MaxForwardInterval>,
    ) -> Result<()> {
        let last = grid[grid.len() - 1];
        if end > last {
            return Err(Error::InvalidInput(format!(
                "point {end} beyond last max forward grid date {last}"
            )));
        }
        let first = grid.partition_point(|&date| date <= start);
        let mut left = if first == 0 {
            self.as_of()
        } else {
            grid[first - 1]
        };
        for &right in &grid[first..] {
            let overlap_end = end.min(right);
            let width = model_time(left, right).years();
            let overlap = model_time(left.max(start), overlap_end).years();
            // 同じ区間幅で比べるので、割り算せず L の変化が最小の候補を選べる。
            // 区間の途中でも候補を比べるのは格子の両端。
            let prepare_candidate = |part: &Arc<ResolvedExpr>| {
                let mut candidate = PreparedExpr::default();
                self.prepare_change(part, left, right, 1.0, &mut candidate)?;
                Ok(candidate)
            };
            intervals.push(MaxForwardInterval {
                factor: factor * (overlap / width),
                first: prepare_candidate(&parts[0])?,
                rest: parts[1..]
                    .iter()
                    .map(prepare_candidate)
                    .collect::<Result<_>>()?,
            });
            if overlap_end == end {
                break;
            }
            left = right;
        }
        Ok(())
    }
}
`````

`````
crates/rhodium-market/src/fx.rs
//! 受渡時刻と担保を指定して読む FX フォワード。

use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use rhodium_core::{
    Error, Result, Warning,
    brand::Brand,
    currency::Currency,
    fx::{CcyPair, Forward, FxRate, Today},
    time::Instant,
};
use scioxide_ad::{Mask, Real};

use crate::{
    curve::{CurveId, CurvePoint, CurveRole},
    market::Market,
    structure::{BlockId, MarketStructure},
};

/// ある担保の下で、ペアの両通貨を割り引く曲線。
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FxCollateral {
    /// 担保区分の識別子。
    pub id: CurveId,
    /// ペアの base 通貨の割引曲線。
    pub base: CurveId,
    /// ペアの quote 通貨の割引曲線。
    pub quote: CurveId,
}

/// 評価時点の FX の入力。市場構造の構築時に検証する。
#[derive(Debug, Clone)]
pub enum FxAnchor {
    /// スポットと、スポットにつながる短期ポイントの θ を読む。
    Direct {
        /// スポット受渡の瞬間。
        spot: Instant,
        /// 短期区間の開始時刻。昇順で、各区間の終わりは次の開始時刻またはスポット。
        short_starts: Vec<Instant>,
        /// 短端を評価時点へ割り戻す担保区分。
        market_collateral: CurveId,
    },
    /// 宣言した中間通貨を経由する。
    Triangulated,
}

/// 評価時点の FX の起点と、担保ごとの割引曲線の対応。
#[derive(Debug, Clone)]
pub struct FxSpec {
    /// 市場の標準の向き。
    pub pair: CcyPair,
    /// 共通の X(0) の入力。
    pub anchor: FxAnchor,
    /// 担保別の割引曲線。
    pub collaterals: Vec<FxCollateral>,
    /// 三角経路。直接ペアにも、整合性の検査のために指定できる。
    pub via: Option<Currency>,
}

impl FxSpec {
    /// スポット、短期ポイント（レートの単位）の順に θ を持つブロック ID。
    ///
    /// # Errors
    ///
    /// ブロック ID に使えない通貨コードの場合は入力エラー。
    pub fn marks_block(pair: CcyPair) -> Result<BlockId> {
        BlockId::new(format!("FX/{pair}").parse()?, "marks")
    }

    pub(crate) fn mark_count(&self) -> Option<usize> {
        match &self.anchor {
            FxAnchor::Direct { short_starts, .. } => Some(1 + short_starts.len()),
            FxAnchor::Triangulated => None,
        }
    }
}

type FxCurves = [usize; 2];

#[derive(Debug)]
pub(crate) struct ResolvedFx {
    collaterals: BTreeMap<CurveId, FxCurves>,
    anchor: Option<ResolvedAnchor>,
    via: Option<Currency>,
}

#[derive(Debug)]
struct ResolvedAnchor {
    at: Instant,
    curves: FxCurves,
    marks: BlockId,
    collateral: CurveId,
}

pub(crate) fn resolve(
    structure: &MarketStructure,
    specs: Vec<FxSpec>,
) -> Result<BTreeMap<CcyPair, ResolvedFx>> {
    let mut result = BTreeMap::new();
    for spec in specs {
        if find(&result, spec.pair).is_some() {
            return Err(Error::InvalidData(format!(
                "duplicate FX pair {}",
                spec.pair
            )));
        }
        let position = |id: &CurveId, currency| {
            structure
                .layout()
                .position(id)
                .filter(|&position| {
                    let spec = &structure.curves[position].spec;
                    spec.currency() == currency && spec.has_role(CurveRole::Discount)
                })
                .ok_or_else(|| {
                    Error::InvalidData(format!(
                        "FX {} requires a {currency} discount curve for {id}",
                        spec.pair
                    ))
                })
        };
        let mut collaterals = BTreeMap::new();
        for collateral in spec.collaterals {
            let curves = [
                position(&collateral.base, spec.pair.base())?,
                position(&collateral.quote, spec.pair.quote())?,
            ];
            if collaterals.insert(collateral.id, curves).is_some() {
                return Err(Error::InvalidData(format!(
                    "duplicate FX collateral for {}",
                    spec.pair
                )));
            }
        }
        let anchor = match spec.anchor {
            FxAnchor::Direct {
                spot,
                short_starts,
                market_collateral,
            } => {
                let at = short_starts.first().copied().unwrap_or(spot);
                if at < structure.as_of()
                    || !short_starts.is_sorted_by(|a, b| a < b)
                    || short_starts.last().is_some_and(|last| *last >= spot)
                {
                    return Err(Error::InvalidData(format!(
                        "invalid FX short settlement chain for {}",
                        spec.pair
                    )));
                }
                let curves = *collaterals.get(&market_collateral).ok_or_else(|| {
                    Error::InvalidData(format!(
                        "missing market collateral {market_collateral} for {}",
                        spec.pair
                    ))
                })?;
                Some(ResolvedAnchor {
                    at,
                    curves,
                    marks: FxSpec::marks_block(spec.pair)?,
                    collateral: market_collateral,
                })
            }
            FxAnchor::Triangulated => None,
        };
        if spec
            .via
            .is_some_and(|via| via == spec.pair.base() || via == spec.pair.quote())
            || (anchor.is_none() && spec.via.is_none())
        {
            return Err(Error::InvalidData(format!(
                "invalid FX triangulation for {}",
                spec.pair
            )));
        }
        result.insert(
            spec.pair,
            ResolvedFx {
                collaterals,
                anchor,
                via: spec.via,
            },
        );
    }
    for (&pair, spec) in &result {
        if let Some(via) = spec.via {
            let mut visiting = BTreeSet::new();
            for leg in cross_legs(pair, via)? {
                validate_path(&result, leg, &mut visiting)?;
            }
        }
    }
    Ok(result)
}

fn find(pairs: &BTreeMap<CcyPair, ResolvedFx>, pair: CcyPair) -> Option<(CcyPair, &ResolvedFx)> {
    pairs
        .get_key_value(&pair)
        .map(|(&key, value)| (key, value))
        .or_else(|| {
            pairs
                .get_key_value(&pair.inverse())
                .map(|(&key, value)| (key, value))
        })
}

fn validate_path(
    pairs: &BTreeMap<CcyPair, ResolvedFx>,
    pair: CcyPair,
    visiting: &mut BTreeSet<CcyPair>,
) -> Result<()> {
    let (key, spec) = find(pairs, pair)
        .ok_or_else(|| Error::InvalidData(format!("undeclared FX triangulation leg {pair}")))?;
    if spec.anchor.is_none() {
        if !visiting.insert(key) {
            return Err(Error::InvalidData(format!(
                "cyclic FX triangulation: {key}"
            )));
        }
        let via = spec
            .via
            .ok_or_else(|| Error::Internal("missing FX triangulation".into()))?;
        for leg in cross_legs(pair, via)? {
            validate_path(pairs, leg, visiting)?;
        }
        visiting.remove(&key);
    }
    Ok(())
}

fn cross_legs(pair: CcyPair, via: Currency) -> Result<[CcyPair; 2]> {
    Ok([
        CcyPair::new(pair.base(), via)?,
        CcyPair::new(via, pair.quote())?,
    ])
}

/// 評価時点の FX を、同じ構造の市場から読む計画。
#[derive(Debug)]
pub struct FxTodayPlan {
    pair: CcyPair,
    source: FxTodaySource,
    structure: Arc<MarketStructure>,
}

#[derive(Debug)]
enum FxTodaySource {
    Direct {
        market_pair: CcyPair,
        anchor: Box<[CurvePoint; 2]>,
        marks: usize,
        block: BlockId,
    },
    Cross(Box<[FxTodayPlan; 2]>),
}

impl FxTodayPlan {
    /// 同じ構造の市場から FX を計算する。
    ///
    /// # Errors
    ///
    /// マークや曲線が利用できないか、FX の値が不正な場合。
    pub fn rate<R: Real>(&self, market: &Market<R>) -> Result<FxRate<Today, R>> {
        if !Arc::ptr_eq(&self.structure, market.structure()) {
            return Err(Error::InvalidInput(
                "FX plan belongs to another market structure".into(),
            ));
        }
        match &self.source {
            FxTodaySource::Direct {
                market_pair,
                anchor,
                marks,
                ..
            } => {
                if !market.available_at(*marks) {
                    return Err(Error::MissingData(format!(
                        "missing FX marks for {market_pair}"
                    )));
                }
                let marks = market.theta_block(*marks);
                let mut outright = marks.at(0);
                if !outright.gt(R::zero()).all() {
                    return Err(Error::InvalidData(format!(
                        "nonpositive FX spot for {market_pair}"
                    )));
                }
                for index in (1..marks.len()).rev() {
                    outright -= marks.at(index);
                    if !outright.is_finite().all() || !outright.gt(R::zero()).all() {
                        return Err(Error::InvalidData(format!(
                            "nonpositive or non-finite short FX outright for {market_pair}"
                        )));
                    }
                }
                let roll = -log_ratio(market, anchor)?;
                let rate = FxRate::today(*market_pair, outright * roll.exp()).map_err(|_| {
                    Error::Numerical(format!("invalid valuation FX rate for {market_pair}"))
                })?;
                if *market_pair == self.pair {
                    Ok(rate)
                } else {
                    rate.invert()
                }
            }
            FxTodaySource::Cross(legs) => legs[0].rate(market)?.cross(&legs[1].rate(market)?),
        }
    }

    fn collect_blocks(&self, blocks: &mut BTreeSet<BlockId>) {
        match &self.source {
            FxTodaySource::Direct { anchor, block, .. } => {
                blocks.insert(block.clone());
                blocks.extend(anchor.iter().flat_map(CurvePoint::blocks).cloned());
            }
            FxTodaySource::Cross(legs) => legs.iter().for_each(|leg| leg.collect_blocks(blocks)),
        }
    }

    /// 計算が読む θ ブロック。
    #[must_use]
    pub fn blocks(&self) -> BTreeSet<BlockId> {
        let mut blocks = BTreeSet::new();
        self.collect_blocks(&mut blocks);
        blocks
    }
}

impl MarketStructure {
    fn fx_spec(&self, pair: CcyPair) -> Result<(CcyPair, &ResolvedFx)> {
        find(&self.fx, pair).ok_or_else(|| Error::MissingData(format!("unknown FX pair {pair}")))
    }

    /// 評価時点の FX に必要なマークと曲線を解決する。
    ///
    /// # Errors
    ///
    /// ペアまたは三角経路が宣言されていない場合。
    pub fn fx_today_plan(self: &Arc<Self>, pair: CcyPair) -> Result<FxTodayPlan> {
        let (market_pair, spec) = self.fx_spec(pair)?;
        let source = if let Some(anchor) = &spec.anchor {
            FxTodaySource::Direct {
                market_pair,
                anchor: Box::new(self.fx_curve_points(anchor.curves, false, anchor.at)?),
                marks: self
                    .layout()
                    .block_position(&anchor.marks)
                    .ok_or_else(|| Error::Internal("no block for FX marks".into()))?,
                block: anchor.marks.clone(),
            }
        } else {
            let via = spec
                .via
                .ok_or_else(|| Error::Internal("missing FX triangulation".into()))?;
            let [first, second] = cross_legs(pair, via)?;
            FxTodaySource::Cross(Box::new([
                self.fx_today_plan(first)?,
                self.fx_today_plan(second)?,
            ]))
        };
        Ok(FxTodayPlan {
            pair,
            source,
            structure: Arc::clone(self),
        })
    }

    /// 支払曲線と同じ組を優先し、無ければ指定した担保またはペアの市場担保を選ぶ。
    /// 戻り値の真偽は支払曲線との一致を示す。
    ///
    /// # Errors
    ///
    /// 指定のペアや担保が無いか、三角ペアの既定の担保が一意でない場合。
    pub fn fx_projection_collateral(
        &self,
        pair: CcyPair,
        discounts: Option<[&CurveId; 2]>,
        preferred: Option<&CurveId>,
    ) -> Result<(&CurveId, bool)> {
        let (market_pair, spec) = self.fx_spec(pair)?;
        let matches = |[base, quote]: FxCurves| {
            let indices = if pair == market_pair {
                [base, quote]
            } else {
                [quote, base]
            };
            discounts.is_some_and(|ids| indices.map(|i| self.curves[i].spec.id()) == ids)
        };
        if let Some(preferred) = preferred {
            return spec
                .collaterals
                .get_key_value(preferred)
                .map(|(id, &curves)| (id, matches(curves)))
                .ok_or_else(|| {
                    Error::MissingData(format!("no {preferred} FX collateral for {pair}"))
                });
        }
        if let Some((id, _)) = spec
            .collaterals
            .iter()
            .find(|(_, curves)| matches(**curves))
        {
            return Ok((id, true));
        }
        if let Some(anchor) = &spec.anchor {
            Ok((&anchor.collateral, false))
        } else if spec.collaterals.len() == 1
            && let Some((id, curves)) = spec.collaterals.first_key_value()
        {
            Ok((id, matches(*curves)))
        } else {
            Err(Error::MissingData(format!(
                "no default FX collateral for {pair}"
            )))
        }
    }

    /// 通貨ペアに宣言された市場担保を取得する。
    ///
    /// # Errors
    ///
    /// ペアが未登録か、三角ペアの担保を一意に選べない場合。
    pub fn fx_market_collateral(&self, pair: CcyPair) -> Result<&CurveId> {
        self.fx_projection_collateral(pair, None, None)
            .map(|(id, _)| id)
    }

    /// 評価時点の FX と受渡時点の曲線を準備する。
    ///
    /// # Errors
    ///
    /// ペア、担保、曲線または受渡時刻を解決できない場合。
    pub fn fx_forward_plan(
        self: &Arc<Self>,
        pair: CcyPair,
        collateral: &CurveId,
        at: Instant,
    ) -> Result<FxForwardPlan> {
        let today = self.fx_today_plan(pair)?;
        let point = self.prepare_fx_point(pair, collateral, at)?;
        let mut blocks = today.blocks();
        blocks.extend(point.curves.iter().flat_map(CurvePoint::blocks).cloned());
        Ok(FxForwardPlan {
            today,
            point,
            blocks,
        })
    }

    fn fx_curve_points(
        self: &Arc<Self>,
        [base, quote]: FxCurves,
        inverse: bool,
        at: Instant,
    ) -> Result<[CurvePoint; 2]> {
        let [base, quote] = if inverse {
            [quote, base]
        } else {
            [base, quote]
        };
        Ok([
            self.prepare_point_at(base, at)?,
            self.prepare_point_at(quote, at)?,
        ])
    }

    /// FX フォワードを指定した UTC の受渡で読むための準備。
    ///
    /// # Errors
    ///
    /// ペア、担保または曲線が利用できないか、受渡の時刻が曲線の範囲外の場合。
    pub fn prepare_fx_point(
        self: &Arc<Self>,
        pair: CcyPair,
        collateral: &CurveId,
        at: Instant,
    ) -> Result<FxForwardPoint> {
        let (market_pair, spec) = self.fx_spec(pair)?;
        let curves = *spec.collaterals.get(collateral).ok_or_else(|| {
            Error::MissingData(format!("no {collateral} FX collateral for {pair}"))
        })?;
        Ok(FxForwardPoint {
            structure: Arc::clone(self),
            pair,
            collateral: collateral.clone(),
            at,
            curves: self.fx_curve_points(curves, market_pair != pair, at)?,
        })
    }

    /// 同じ担保の二つの受渡の FX ポイントを計算する準備。
    ///
    /// # Errors
    ///
    /// 受渡が昇順でないか、必要な曲線と時刻を解決できない場合。
    pub fn prepare_fx_interval(
        self: &Arc<Self>,
        pair: CcyPair,
        collateral: &CurveId,
        start: Instant,
        end: Instant,
    ) -> Result<FxForwardInterval> {
        if start >= end {
            return Err(Error::InvalidInput(
                "FX interval needs increasing settlements".into(),
            ));
        }
        let today = self.fx_today_plan(pair)?;
        let start = self.prepare_fx_point(pair, collateral, start)?;
        let end = self.prepare_fx_point(pair, collateral, end)?;
        let mut blocks = BTreeSet::new();
        today.collect_blocks(&mut blocks);
        blocks.extend(
            start
                .curves
                .iter()
                .chain(&end.curves)
                .flat_map(CurvePoint::blocks)
                .cloned(),
        );
        Ok(FxForwardInterval {
            today,
            start,
            end,
            blocks: blocks.into_iter().collect(),
        })
    }
}

/// θ から一つの FX フォワードを評価する準備済みの計画。
#[derive(Debug)]
pub struct FxForwardPlan {
    today: FxTodayPlan,
    point: FxForwardPoint,
    blocks: BTreeSet<BlockId>,
}

impl FxForwardPlan {
    /// 計算が読む θ ブロック。
    #[must_use]
    pub const fn blocks(&self) -> &BTreeSet<BlockId> {
        &self.blocks
    }

    pub(crate) fn rate<R: Real>(&self, market: &Market<R>) -> Result<R> {
        if !Arc::ptr_eq(&self.point.structure, market.structure()) {
            return Err(Error::InvalidInput(
                "FX plan belongs to another market structure".into(),
            ));
        }
        let value = self.today.rate(market)?.rate() * log_ratio(market, &self.point.curves)?.exp();
        Ok(FxRate::forward(self.point.pair, value, self.point.at)
            .map_err(|_| Error::Numerical(format!("invalid FX forward for {}", self.point.pair)))?
            .rate())
    }
}

fn log_ratio<R: Real>(market: &Market<R>, [base, quote]: &[CurvePoint; 2]) -> Result<R> {
    market.require_blocks(base.block_positions().chain(quote.block_positions()))?;
    Ok(base.log_df(market) - quote.log_df(market))
}

/// 通貨ペア、担保と受渡時刻に対応する準備済みの FX の点。
#[derive(Debug)]
pub struct FxForwardPoint {
    structure: Arc<MarketStructure>,
    pair: CcyPair,
    collateral: CurveId,
    at: Instant,
    curves: [CurvePoint; 2],
}

/// 受渡区間の FX ポイントの準備。同じ構造の別の θ に再利用できる。
#[derive(Debug)]
pub struct FxForwardInterval {
    today: FxTodayPlan,
    start: FxForwardPoint,
    end: FxForwardPoint,
    blocks: Vec<BlockId>,
}

impl FxForwardInterval {
    /// 計算が読む θ ブロック。
    #[must_use]
    pub fn blocks(&self) -> &[BlockId] {
        &self.blocks
    }

    /// 二つのフォワードの差をレートの単位で求める。
    ///
    /// # Errors
    ///
    /// 構造が異なるか、必要な値を利用できない場合。
    pub fn points<R: Real>(&self, market: &Market<R>) -> Result<R> {
        if !Arc::ptr_eq(&self.start.structure, market.structure()) {
            return Err(Error::InvalidInput(
                "FX interval belongs to another market structure".into(),
            ));
        }
        let start = log_ratio(market, &self.start.curves)?;
        let end = log_ratio(market, &self.end.curves)?;
        // 近接した受渡の大きなアウトライトを引き合わず、短期ポイントの桁を保つ。
        let points = self.today.rate(market)?.rate() * start.exp() * (end - start).exp_m1();
        if !points.is_finite().all() {
            return Err(Error::Numerical(format!(
                "non-finite FX points for {}",
                self.today.pair
            )));
        }
        Ok(points)
    }
}

/// 一つのビューで担保とペアを解決した FX ハンドル。
#[derive(Debug)]
pub struct FxRef<'v, R = f64> {
    market: Market<R>,
    collateral: CurveId,
    today: FxRate<Today, R>,
    _brand: Brand<'v>,
}

impl<'v, R: Real> FxRef<'v, R> {
    pub(crate) fn new(
        market: &Market<R>,
        pair: CcyPair,
        collateral: &CurveId,
        brand: Brand<'v>,
    ) -> Result<Self> {
        let (_, spec) = market.structure().fx_spec(pair)?;
        if !spec.collaterals.contains_key(collateral) {
            return Err(Error::MissingData(format!(
                "no {collateral} FX collateral for {pair}"
            )));
        }
        Ok(Self {
            market: market.clone(),
            collateral: collateral.clone(),
            today: market.fx_today(pair)?,
            _brand: brand,
        })
    }

    /// 担保共通の評価時点のレート。
    #[must_use]
    pub const fn today(&self) -> FxRate<Today, R> {
        self.today
    }

    /// 指定した時刻の FX フォワード。
    ///
    /// # Errors
    ///
    /// 受渡が曲線の範囲外、必要なブロックが欠損、または値を表現できない場合。
    pub fn forward(&self, at: Instant) -> Result<FxRate<Forward, R>> {
        let point =
            self.market
                .structure()
                .prepare_fx_point(self.today.pair(), &self.collateral, at)?;
        self.bind_forward(&point)
    }

    /// 準備済みの受渡点を市場に結合し、フォワードを計算する。
    ///
    /// # Errors
    ///
    /// ペア、担保、構造が点と異なるか、値を計算できない場合。
    pub fn bind_forward(&self, point: &FxForwardPoint) -> Result<FxRate<Forward, R>> {
        if !Arc::ptr_eq(self.market.structure(), &point.structure)
            || point.pair != self.today.pair()
            || point.collateral != self.collateral
        {
            return Err(Error::InvalidInput(
                "FX point belongs to another pair, collateral or market structure".into(),
            ));
        }
        let value = self.today.rate() * log_ratio(&self.market, &point.curves)?.exp();
        FxRate::forward(point.pair, value, point.at)
            .map_err(|_| Error::Numerical(format!("invalid FX forward for {}", point.pair)))
    }
}

impl<R: Real> Market<R> {
    pub(crate) fn fx_today(&self, pair: CcyPair) -> Result<FxRate<Today, R>> {
        self.structure().fx_today_plan(pair)?.rate(self)
    }
}

impl Market<f64> {
    /// 値が利用できる直接ペアと宣言された三角経路との不一致を報告する。
    #[must_use]
    pub fn fx_cross_warnings(&self) -> Vec<Warning> {
        self.structure()
            .fx
            .iter()
            .filter_map(|(&pair, spec)| {
                let via = spec.via?;
                spec.anchor.as_ref()?;
                let [first, second] = cross_legs(pair, via).ok()?;
                let direct = self.fx_today(pair).ok()?.rate();
                let cross = self
                    .fx_today(first)
                    .ok()?
                    .cross(&self.fx_today(second).ok()?)
                    .ok()?
                    .rate();
                ((direct - cross).abs() > (1e-12 * direct.abs().max(cross.abs())).max(1e-14))
                    .then_some(Warning::FxCrossInconsistent { pair })
            })
            .collect()
    }
}
`````

`````
crates/rhodium-market/src/lib.rs
//! 市場の構造、quote、曲線とボラ面。

#![forbid(unsafe_code)]

pub mod curve;
pub mod fx;
pub mod market;
pub mod quote;
pub mod structure;
pub mod view;
`````

`````
crates/rhodium-market/src/market.rs
//! 同じ構造で差し替えられる平坦な市場パラメータ θ。

use std::sync::{Arc, OnceLock};

use rhodium_core::{Error, Result};
use scioxide_ad::{Mask, Real, Var};

use crate::structure::{BlockId, MarketStructure};

/// 市場の構造と、明示的に設定した θ。
#[derive(Debug, Clone)]
pub struct Market<R = f64> {
    structure: Arc<MarketStructure>,
    state: Arc<MarketState<R>>,
}

// θ の持ち方。曲線や FX の読み出しは theta_block と available_at を通し、どちらの形かを意識しない。
#[derive(Debug)]
enum MarketState<R> {
    // θ の全体と、配置の位置ごとの利用可否を持つ。校正やバンプの市場はこの形。
    Dense {
        theta: Vec<R>,
        available: Vec<bool>,
    },
    // AAD の市場。校正済みの passive な市場を取引の間で共有し、取引が読むブロックだけを
    // tape の変数で持つ。ほかのブロックは読むたびに passive の値を定数にし、全 θ を写さない。
    Active {
        // 構造、値、利用可否の元になる校正済みの市場。
        passive: Market<f64>,
        // (配置の位置, 変数) を位置の昇順に並べ、theta_block が二分探索で引く。
        blocks: Vec<(usize, Box<[R]>)>,
        // 全要素のスライスを返す既存の API（theta）のために、呼ばれたときだけ作る写し。
        flat: OnceLock<Vec<R>>,
    },
}

// 一つの θ ブロックの読み出し。添字はブロックの先頭からの位置。
#[derive(Clone, Copy)]
pub(crate) enum ThetaBlock<'a, R> {
    // 市場の数値型の値。Dense の θ か、AAD の市場で tape の変数にしたブロック。
    Values(&'a [R]),
    // AAD の市場が共有する passive なブロックの値。読み出しのたびに R::constant にする。
    Constants(&'a [f64]),
}

impl<R: Real> ThetaBlock<'_, R> {
    pub(crate) fn at(&self, index: usize) -> R {
        match self {
            Self::Values(values) => values[index],
            Self::Constants(values) => R::constant(values[index]),
        }
    }

    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Values(values) => values.len(),
            Self::Constants(values) => values.len(),
        }
    }
}

impl Market<f64> {
    /// この市場を共有し、指定した θ ブロックを tape の変数に置き換える。
    ///
    /// ほかのブロックは読み出した値を定数にする。構造と利用可否は元の市場を引き継ぐ。
    ///
    /// # Errors
    ///
    /// ブロックが未知、利用不可、重複、または変数の長さや値が不正なら `InvalidInput`。
    #[expect(
        clippy::arc_with_non_send_sync,
        reason = "数値型を変えても Market と view で同じ状態を共有する。Var は tape の内側に閉じる"
    )]
    pub fn with_active_theta<'a, 't: 'a>(
        &self,
        active: impl IntoIterator<Item = (&'a BlockId, &'a [Var<'t, f64>])>,
    ) -> Result<Market<Var<'t, f64>>> {
        let layout = self.structure.layout();
        let mut blocks = active
            .into_iter()
            .map(|(block, values)| {
                let position = layout
                    .block_position(block)
                    .filter(|&position| self.available_at(position))
                    .ok_or_else(|| {
                        Error::InvalidInput(format!(
                            "active block is unavailable: {} {}",
                            block.curve(),
                            block.part()
                        ))
                    })?;
                if layout.ranges[position].len() != values.len()
                    || values.iter().any(|value| !value.is_finite())
                {
                    return Err(Error::InvalidInput(format!(
                        "invalid active theta in {} {}",
                        block.curve(),
                        block.part()
                    )));
                }
                Ok((position, Box::from(values)))
            })
            .collect::<Result<Vec<_>>>()?;
        blocks.sort_unstable_by_key(|(position, _)| *position);
        if blocks.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(Error::InvalidInput("duplicate active block".into()));
        }
        Ok(Market {
            structure: Arc::clone(&self.structure),
            state: Arc::new(MarketState::Active {
                passive: self.clone(),
                blocks,
                flat: OnceLock::new(),
            }),
        })
    }
}

impl<R: Real> Market<R> {
    /// θ の全ブロックが利用できる市場を作る。
    ///
    /// # Errors
    ///
    /// θ の長さが配置と一致しないか、値が有限でなければ `InvalidInput`。
    pub fn from_theta(structure: Arc<MarketStructure>, theta: Vec<R>) -> Result<Self> {
        let available = vec![true; structure.layout().ranges.len()];
        Self::build(structure, theta, available)
    }

    /// 指定したブロックだけ利用できる部分市場を作る。
    ///
    /// 利用できないブロックにも、θ の同じ範囲が割り当てられる。
    ///
    /// # Errors
    ///
    /// θ の長さが配置と一致しない、利用できるブロックが構造に無い、または
    /// 利用できる値が有限でなければ `InvalidInput`。
    pub fn from_partial_theta(
        structure: Arc<MarketStructure>,
        theta: Vec<R>,
        available: impl IntoIterator<Item = BlockId>,
    ) -> Result<Self> {
        let layout = structure.layout();
        let mut mask = vec![false; layout.ranges.len()];
        for block in available {
            let position = layout.block_position(&block).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "unknown available block {} {}",
                    block.curve(),
                    block.part()
                ))
            })?;
            mask[position] = true;
        }
        Self::build(structure, theta, mask)
    }

    fn build(structure: Arc<MarketStructure>, theta: Vec<R>, available: Vec<bool>) -> Result<Self> {
        let layout = structure.layout();
        if theta.len() != layout.len() {
            return Err(Error::InvalidInput(format!(
                "theta length {}, expected {}",
                theta.len(),
                layout.len()
            )));
        }
        // 利用できないブロックは読み出しの前に MissingData で止まるので、値を検査しない。
        for (position, range) in layout.ranges.iter().enumerate() {
            if available[position]
                && theta[range.clone()]
                    .iter()
                    .any(|value| !value.is_finite().all())
            {
                let block = &layout.ids[position];
                return Err(Error::InvalidInput(format!(
                    "non-finite theta in {} {}",
                    block.curve(),
                    block.part()
                )));
            }
        }
        Ok(Self {
            structure,
            state: Arc::new(MarketState::Dense { theta, available }),
        })
    }

    /// 準備済みの点やスロットと共有する市場の構造。
    #[must_use]
    pub const fn structure(&self) -> &Arc<MarketStructure> {
        &self.structure
    }

    /// θ の全要素。利用できないブロックの範囲も含む。
    ///
    /// AAD の市場はこのメソッドを呼んだときにだけ全要素を平坦な列として確保する。
    #[must_use]
    pub fn theta(&self) -> &[R] {
        match self.state.as_ref() {
            MarketState::Dense { theta, .. } => theta,
            MarketState::Active {
                passive,
                blocks,
                flat,
            } => flat.get_or_init(|| {
                let mut theta: Vec<_> = passive.theta().iter().copied().map(R::constant).collect();
                for (position, values) in blocks {
                    theta[self.structure.layout().ranges[*position].clone()]
                        .copy_from_slice(values);
                }
                theta
            }),
        }
    }

    // 配置の位置 position のブロックを、θ の全体を確保せずに読む。
    // 利用可否は確かめないので、呼び出し側が require_blocks などで先に確かめる。
    pub(crate) fn theta_block(&self, position: usize) -> ThetaBlock<'_, R> {
        let range = &self.structure.layout().ranges[position];
        match self.state.as_ref() {
            MarketState::Dense { theta, .. } => ThetaBlock::Values(&theta[range.clone()]),
            MarketState::Active {
                passive, blocks, ..
            } => match blocks.binary_search_by_key(&position, |(position, _)| *position) {
                Ok(index) => ThetaBlock::Values(&blocks[index].1),
                Err(_) => ThetaBlock::Constants(&passive.theta()[range.clone()]),
            },
        }
    }

    // 配置の位置 position のブロックを利用できるか。AAD の市場は元の市場の利用可否をそのまま使う
    // （with_active_theta は利用できないブロックを active にしない）。
    pub(crate) fn available_at(&self, position: usize) -> bool {
        match self.state.as_ref() {
            MarketState::Dense { available, .. } => available[position],
            MarketState::Active { passive, .. } => passive.available_at(position),
        }
    }

    /// 指定したブロックが存在し、利用できるか。
    #[must_use]
    pub fn block_is_available(&self, block: &BlockId) -> bool {
        self.structure
            .layout()
            .block_position(block)
            .is_some_and(|position| self.available_at(position))
    }

    pub(crate) fn require_curve(&self, position: usize) -> Result<()> {
        self.require_blocks(self.structure.curves[position].blocks.iter().copied())
    }

    pub(crate) fn require_blocks(&self, blocks: impl IntoIterator<Item = usize>) -> Result<()> {
        if let Some(position) = blocks
            .into_iter()
            .find(|&position| !self.available_at(position))
        {
            let block = &self.structure.layout().ids[position];
            return Err(Error::MissingData(format!(
                "market block is unavailable: {} {}",
                block.curve(),
                block.part()
            )));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use rhodium_core::currency::Currency;
    use rhodium_testkit::assert_close;

    use super::*;
    use crate::curve::{CurveId, CurveRole, CurveSpec, Extrapolation};

    #[test]
    fn active_market_reads_passive_blocks_without_materializing_theta() -> Result<()> {
        // 置き換えた値、元のまま読める値、校正に失敗したブロックを一つずつ作る。
        let [active, missing, passive] =
            ["ACTIVE", "MISSING", "PASSIVE"].map(str::parse::<CurveId>);
        let (active, missing, passive) = (active?, missing?, passive?);
        let pillar = "2027-01-04T00:00:00Z".parse()?;
        let curves = [&active, &missing, &passive]
            .into_iter()
            .map(|id| {
                CurveSpec::new(
                    id.clone(),
                    Currency::JPY,
                    [CurveRole::Discount],
                    [pillar],
                    Extrapolation::FlatAverageForward,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        let structure = Arc::new(MarketStructure::new(
            "2026-01-04T00:00:00Z".parse()?,
            curves,
        )?);
        let active = BlockId::nodes(active);
        let missing = BlockId::nodes(missing);
        let passive = BlockId::nodes(passive);
        let base = Market::from_partial_theta(
            structure,
            vec![-0.02, f64::NAN, -0.04],
            [active.clone(), passive.clone()],
        )?;
        let replacements = [Var::constant(-0.03)];
        let market = base.with_active_theta([(&active, replacements.as_slice())])?;
        assert!(market.block_is_available(&active));
        assert!(!market.block_is_available(&missing));
        assert!(market.block_is_available(&passive));
        let log_df = |block: &BlockId| {
            let point = market.structure().prepare_point(block.curve(), pillar)?;
            Ok::<_, Error>(point.log_df(&market).value())
        };
        assert_close!(log_df(&active)?, -0.03, 1e-12, 1e-14);
        assert_close!(log_df(&passive)?, -0.04, 1e-12, 1e-14);
        let MarketState::Active { flat, .. } = market.state.as_ref() else {
            return Err(Error::Internal("expected an active market".into()));
        };
        assert!(flat.get().is_none());
        assert_close!(market.theta()[2].value(), -0.04, 1e-12, 1e-14);
        assert!(market.theta()[1].value().is_nan());
        assert!(matches!(
            base.with_active_theta([(&missing, replacements.as_slice())]),
            Err(Error::InvalidInput(_))
        ));
        Ok(())
    }
}
`````

`````
crates/rhodium-market/src/quote.rs
//! 市場で観測した quote と、校正が参照する入力のマーク。

use std::{
    borrow::Cow,
    collections::{BTreeMap, BTreeSet},
    str::FromStr,
};

use rhodium_core::{Error, Result};
use serde::{Deserialize, Serialize};

pub use rhodium_core::id::QuoteId;

/// 評価時点の quote のスナップショット。金利は decimal（0.85% は 0.0085）、FX ポイントは pip 数。
#[derive(Debug, Clone, PartialEq, Default)]
pub struct QuoteSet(BTreeMap<QuoteId, f64>);

impl QuoteSet {
    /// 受け取った quote を ID 順に保持する。
    ///
    /// # Errors
    ///
    /// 値が有限でないか、同じ ID が複数回現れたときは `InvalidInput`。
    pub fn new(entries: impl IntoIterator<Item = (QuoteId, f64)>) -> Result<Self> {
        let mut quotes = BTreeMap::new();
        for (id, value) in entries {
            if !value.is_finite() {
                return Err(Error::InvalidInput(format!("non-finite quote: {id}")));
            }
            match quotes.entry(id) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(value);
                }
                std::collections::btree_map::Entry::Occupied(entry) => {
                    return Err(Error::InvalidInput(format!(
                        "duplicate quote: {}",
                        entry.key()
                    )));
                }
            }
        }
        Ok(Self(quotes))
    }

    /// ID に対応する値を返す。
    #[must_use]
    pub fn get(&self, id: &QuoteId) -> Option<f64> {
        self.0.get(id).copied()
    }

    /// スナップショットの値を ID の辞書順で返す。
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (&QuoteId, f64)> {
        self.0.iter().map(|(id, &value)| (id, value))
    }

    /// 追加・削除されたか、値の変わった quote の ID を返す。
    #[must_use]
    pub fn changed_ids(&self, other: &Self) -> BTreeSet<QuoteId> {
        self.0
            .iter()
            .filter(|(id, value)| other.0.get(*id) != Some(*value))
            .map(|(id, _)| id.clone())
            .chain(
                other
                    .0
                    .keys()
                    .filter(|id| !self.0.contains_key(*id))
                    .cloned(),
            )
            .collect()
    }

    /// 記述順に上書きを適用したスナップショットと、適用前後で値が変わった ID を返す。
    ///
    /// # Errors
    ///
    /// 現在のスナップショットにパターンが一致しない、登録されていない ID に一致する、
    /// 演算結果が非有限か、倍率を掛ける前後の値が正でない場合は `InvalidInput`。
    pub fn with_overrides(&self, overrides: &QuoteOverrides) -> Result<(Self, BTreeSet<QuoteId>)> {
        let mut quotes = self.clone();
        for operation in &overrides.operations {
            match operation {
                QuoteOverride::Set { id, value } => {
                    quotes.0.insert(id.clone(), *value);
                }
                QuoteOverride::Remove { id } => {
                    quotes.0.remove(id);
                }
                QuoteOverride::Add { pattern, shift } => {
                    quotes.apply_shift(pattern, &overrides.kinds, |value| Some(value + shift))?;
                }
                QuoteOverride::Scale { pattern, factor } => {
                    quotes.apply_shift(pattern, &overrides.kinds, |value| {
                        let next = value * factor;
                        (value > 0.0 && next > 0.0).then_some(next)
                    })?;
                }
            }
        }
        let changed = self.changed_ids(&quotes);
        Ok((quotes, changed))
    }

    fn apply_shift(
        &mut self,
        pattern: &QuotePattern,
        kinds: &BTreeMap<QuoteId, QuoteKind>,
        transform: impl Fn(f64) -> Option<f64>,
    ) -> Result<()> {
        let mut matching = self
            .0
            .iter_mut()
            .filter(|(id, _)| pattern.matches(id))
            .peekable();
        if matching.peek().is_none() {
            return Err(pattern.no_match());
        }
        for (id, value) in matching {
            require_quote(kinds, id)?;
            *value = transform(*value)
                .filter(|next| next.is_finite())
                .ok_or_else(|| Error::InvalidInput(format!("quote shift out of range: {id}")))?;
        }
        Ok(())
    }
}

/// quote ID を区間ごとに照合する検証済みのパターン。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuotePattern {
    text: String,
    segments: Vec<PatternSegment>,
    tail: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum PatternSegment {
    Any,
    OneOf(Vec<String>),
}

impl FromStr for QuotePattern {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        let invalid = || Error::InvalidInput(format!("invalid quote pattern: {text:?}"));
        let mut parts = text.split('.').peekable();
        let mut segments = Vec::new();
        let mut tail = false;
        while let Some(part) = parts.next() {
            match part {
                "*" => segments.push(PatternSegment::Any),
                "**" if parts.peek().is_none() => tail = true,
                _ => {
                    let alternatives = part
                        .strip_prefix('{')
                        .and_then(|part| part.strip_suffix('}'));
                    let names: Vec<_> =
                        alternatives.map_or_else(|| vec![part], |list| list.split(',').collect());
                    if !names.iter().all(|name| QuoteId::is_valid_segment(name)) {
                        return Err(invalid());
                    }
                    segments.push(PatternSegment::OneOf(
                        names.into_iter().map(str::to_owned).collect(),
                    ));
                }
            }
        }
        Ok(Self {
            text: text.into(),
            segments,
            tail,
        })
    }
}

impl QuotePattern {
    /// `*` は一つの区間、`{A,B}` は列挙した一つ、末尾の `**` は残りの 0 個以上に一致する。
    #[must_use]
    pub fn matches(&self, id: &QuoteId) -> bool {
        let mut parts = id.as_str().split('.');
        self.segments.iter().all(|segment| {
            parts.next().is_some_and(|part| match segment {
                PatternSegment::Any => true,
                PatternSegment::OneOf(names) => names.iter().any(|name| name == part),
            })
        }) && (self.tail || parts.next().is_none())
    }

    fn no_match(&self) -> Error {
        Error::InvalidInput(format!("quote pattern matches no quote: {}", self.text))
    }
}

/// 一つの quote 上書き。値と加算幅は入力 quote と同じ単位。
#[derive(Debug, Clone)]
pub enum QuoteOverride {
    /// 値を設定する。登録済みの ID はスナップショットに追加できる。
    Set {
        /// 対象の ID。
        id: QuoteId,
        /// 設定する値。
        value: f64,
    },
    /// パターンに一致する値に加算する。
    Add {
        /// 対象のパターン。
        pattern: QuotePattern,
        /// 加算幅。
        shift: f64,
    },
    /// パターンに一致する正の種類の値に、正の係数を掛ける。
    Scale {
        /// 対象のパターン。
        pattern: QuotePattern,
        /// 倍率。
        factor: f64,
    },
    /// スナップショットから値を削除する。
    Remove {
        /// 対象の ID。
        id: QuoteId,
    },
}

/// 市場の quote の種類に対して検証済みの、順序付きの上書き。
#[derive(Debug, Clone)]
pub struct QuoteOverrides {
    operations: Vec<QuoteOverride>,
    kinds: BTreeMap<QuoteId, QuoteKind>,
}

impl QuoteOverrides {
    /// 市場に登録された quote とその種類に照らして上書きを作る。
    ///
    /// # Errors
    ///
    /// 対象が未登録、パターンが登録済みの ID に一致しない、演算の引数が非有限か、
    /// 倍率が正でないか対象に倍率を掛けられない種類が含まれる場合は `InvalidInput`。
    pub fn new(
        operations: impl IntoIterator<Item = QuoteOverride>,
        kinds: &BTreeMap<QuoteId, QuoteKind>,
    ) -> Result<Self> {
        let operations: Vec<_> = operations.into_iter().collect();
        for operation in &operations {
            match operation {
                QuoteOverride::Set { id, value } => {
                    require_quote(kinds, id)?;
                    if !value.is_finite() {
                        return Err(Error::InvalidInput(format!("non-finite quote: {id}")));
                    }
                }
                QuoteOverride::Remove { id } => require_quote(kinds, id)?,
                QuoteOverride::Add { pattern, shift } => {
                    if !shift.is_finite() {
                        return Err(Error::InvalidInput(format!(
                            "non-finite quote shift: {}",
                            pattern.text
                        )));
                    }
                    check_pattern(kinds, pattern, false)?;
                }
                QuoteOverride::Scale { pattern, factor } => {
                    if !factor.is_finite() || *factor <= 0.0 {
                        return Err(Error::InvalidInput(format!(
                            "invalid quote scale factor: {factor}"
                        )));
                    }
                    check_pattern(kinds, pattern, true)?;
                }
            }
        }
        Ok(Self {
            operations,
            kinds: kinds.clone(),
        })
    }

    /// 別の市場の種類表に適用する場合は、その表に照らして上書きを検証する。
    ///
    /// # Errors
    ///
    /// [`Self::new`] と同じ。
    pub fn for_kinds(&self, kinds: &BTreeMap<QuoteId, QuoteKind>) -> Result<Cow<'_, Self>> {
        if &self.kinds == kinds {
            Ok(Cow::Borrowed(self))
        } else {
            Self::new(self.operations.iter().cloned(), kinds).map(Cow::Owned)
        }
    }
}

fn require_quote(kinds: &BTreeMap<QuoteId, QuoteKind>, id: &QuoteId) -> Result<()> {
    kinds
        .contains_key(id)
        .then_some(())
        .ok_or_else(|| Error::InvalidInput(format!("quote is not registered: {id}")))
}

fn check_pattern(
    kinds: &BTreeMap<QuoteId, QuoteKind>,
    pattern: &QuotePattern,
    scale: bool,
) -> Result<()> {
    let mut matching = kinds
        .iter()
        .filter(|(id, _)| pattern.matches(id))
        .peekable();
    if matching.peek().is_none() {
        return Err(pattern.no_match());
    }
    if scale && let Some((id, kind)) = matching.find(|&(_, kind)| !kind.can_scale()) {
        return Err(Error::InvalidInput(format!(
            "cannot scale quote {id} of kind {kind:?}"
        )));
    }
    Ok(())
}

/// 校正の binding が宣言する quote の種類。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum QuoteKind {
    /// パーレートなどの金利。
    Rate,
    /// 金利差などの符号付きスプレッド。
    Spread,
    /// 信用商品のランニングスプレッド。
    CreditSpread,
    /// 先物の価格。
    FuturesPrice,
    /// 価格。
    Price,
    /// 額面に対するアップフロントの比。
    Upfront,
    /// 正規ボラティリティ。
    NormalVol,
    /// 対数正規ボラティリティ。
    LognormalVol,
    /// シフト付き対数正規ボラティリティ。
    ShiftedLognormalVol,
    /// 正規ボラティリティの差。
    NormalVolDifference,
    /// 対数正規ボラティリティの差。
    LognormalVolDifference,
    /// FX スポット。
    FxSpot,
    /// FX フォワードポイント。
    FxForwardPoints,
    /// 資産のスポット価格。
    AssetSpot,
    /// 相関。
    Correlation,
    /// 回収率。
    Recovery,
    /// モデルなどの入力のマーク。
    Mark,
}

impl QuoteKind {
    const fn can_scale(self) -> bool {
        matches!(
            self,
            Self::NormalVol
                | Self::LognormalVol
                | Self::ShiftedLognormalVol
                | Self::CreditSpread
                | Self::FxSpot
                | Self::AssetSpot
                | Self::Price
        )
    }
}
`````

`````
crates/rhodium-market/src/structure.rs
//! Quote に依存しない市場の構造と、曲線パラメータ θ の配置。

use std::{
    borrow::Cow,
    collections::{BTreeMap, BTreeSet},
    ops::Range,
    sync::Arc,
};

use rhodium_core::{
    Error, Result,
    currency::Currency,
    fx::CcyPair,
    time::{Instant, ModelTime, model_time},
};

use crate::curve::{CurveExpr, CurveId, CurveRole, CurveSpec, Extrapolation};
use crate::fx::{self, FxSpec, ResolvedFx};

/// 曲線のパラメータブロックの識別子。
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BlockId(CurveId, Cow<'static, str>);

impl BlockId {
    /// 指定した曲線の節点ブロックを指す ID。
    #[must_use]
    pub const fn nodes(curve: CurveId) -> Self {
        Self(curve, Cow::Borrowed("nodes"))
    }

    /// 指定した曲線の名前付き部品を指す ID。
    ///
    /// # Errors
    ///
    /// 部品名が空か、英数字、`-`、`_`、`.` 以外を含む場合は `InvalidInput`。
    pub fn new(curve: CurveId, part: &str) -> Result<Self> {
        if part.is_empty()
            || !part
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"-_.".contains(&c))
        {
            return Err(Error::InvalidInput(format!(
                "invalid part name {part:?} on {curve}"
            )));
        }
        Ok(Self::for_part(curve, part))
    }

    pub(crate) fn for_part(curve: CurveId, part: &str) -> Self {
        if part == "nodes" {
            Self::nodes(curve)
        } else {
            Self(curve, Cow::Owned(part.into()))
        }
    }

    /// ブロックが属する曲線。
    #[must_use]
    pub const fn curve(&self) -> &CurveId {
        &self.0
    }

    /// 曲線の中でブロックを識別する部品名。
    #[must_use]
    pub fn part(&self) -> &str {
        &self.1
    }
}

/// θ の一つのパラメータの識別子。
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ParamId {
    block: BlockId,
    element: u32,
}

impl ParamId {
    /// ブロックとブロック内の 0 始まりの位置でパラメータを指定する。
    #[must_use]
    pub const fn new(block: BlockId, element: u32) -> Self {
        Self { block, element }
    }

    /// パラメータが属するブロック。
    #[must_use]
    pub const fn block(&self) -> &BlockId {
        &self.block
    }

    /// ブロック内の 0 始まりの位置。
    #[must_use]
    pub const fn element(&self) -> u32 {
        self.element
    }
}

/// 曲線 ID と部品名の辞書順に並べた θ の配置。原点は θ に含まない。
#[derive(Debug)]
pub struct ParamLayout {
    curve_positions: BTreeMap<CurveId, usize>,
    block_positions: BTreeMap<BlockId, usize>,
    pub(crate) ids: Vec<BlockId>,
    pub(crate) ranges: Vec<Range<usize>>,
    len: usize,
}

impl ParamLayout {
    /// 平坦な θ の要素数。
    #[must_use]
    pub const fn len(&self) -> usize {
        self.len
    }

    /// θ が空か。
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// 指定したブロックが占める θ の範囲。
    #[must_use]
    pub fn range(&self, block: &BlockId) -> Option<Range<usize>> {
        self.block_position(block)
            .map(|position| self.ranges[position].clone())
    }

    /// θ の配置順のブロック ID。
    #[must_use]
    pub fn blocks(&self) -> &[BlockId] {
        &self.ids
    }

    /// 指定したパラメータの θ の添字。対象のブロックか要素が無ければ `None`。
    #[must_use]
    pub fn index(&self, param: &ParamId) -> Option<usize> {
        let range = self.range(param.block())?;
        let element = usize::try_from(param.element()).ok()?;
        (element < range.len()).then(|| range.start + element)
    }

    pub(crate) fn position(&self, id: &CurveId) -> Option<usize> {
        self.curve_positions.get(id).copied()
    }

    pub(crate) fn block_position(&self, block: &BlockId) -> Option<usize> {
        self.block_positions.get(block).copied()
    }
}

#[derive(Debug)]
pub(crate) struct Curve {
    pub(crate) spec: CurveSpec,
    pub(crate) expr: Arc<ResolvedExpr>,
    pub(crate) blocks: Vec<usize>,
    /// 予想曲線の予想を定める担保の割引曲線（予想の測度）の位置。分からなければ `None`。
    pub(crate) forecast_measure: Option<usize>,
}

#[derive(Debug)]
pub(crate) enum ResolvedExpr {
    Nodes {
        block: usize,
        extrapolation: Extrapolation,
    },
    Sum(Vec<Arc<Self>>),
    Scale(f64, Arc<Self>),
    Flat(f64),
    MaxForward {
        parts: Box<[Arc<Self>]>,
        grid: Box<[Instant]>,
    },
    Stitch {
        short: Arc<Self>,
        long: Arc<Self>,
        at: Instant,
    },
    Steps(usize),
    Jumps(usize),
}

impl ResolvedExpr {
    fn collect_blocks(&self, blocks: &mut BTreeSet<usize>) {
        match self {
            Self::Nodes { block, .. } | Self::Steps(block) | Self::Jumps(block) => {
                blocks.insert(*block);
            }
            Self::Sum(parts) => parts.iter().for_each(|part| part.collect_blocks(blocks)),
            Self::MaxForward { parts, .. } => {
                parts.iter().for_each(|part| part.collect_blocks(blocks));
            }
            Self::Scale(_, body) => body.collect_blocks(blocks),
            Self::Stitch { short, long, .. } => {
                short.collect_blocks(blocks);
                long.collect_blocks(blocks);
            }
            Self::Flat(_) => {}
        }
    }
}

#[derive(Debug)]
pub(crate) struct BlockGrid {
    pub(crate) dates: Box<[Instant]>,
    pub(crate) times: Box<[ModelTime]>,
}

/// 評価の瞬間、曲線の定義、θ の配置を持つ検証済みの市場構造。
#[derive(Debug)]
pub struct MarketStructure {
    as_of: Instant,
    pub(crate) curves: Vec<Curve>,
    pub(crate) grids: Vec<BlockGrid>,
    pub(crate) fx: BTreeMap<CcyPair, ResolvedFx>,
    layout: ParamLayout,
}

impl MarketStructure {
    /// 曲線を指定した評価の瞬間に結び付け、節点ブロックを配置する。
    ///
    /// # Errors
    ///
    /// 曲線 ID の重複、ピラーや `MaxForward` の格子が評価時点以前にあるかモデル時間で
    /// 重複する場合、ブロックの要素が添字の範囲を超える場合は `InvalidInput`。
    pub fn new(as_of: Instant, curves: impl IntoIterator<Item = CurveSpec>) -> Result<Self> {
        Self::with_fx(as_of, curves, [])
    }

    /// 曲線、担保別の FX と θ の配置を評価の瞬間に結ぶ。
    ///
    /// # Errors
    ///
    /// 曲線の構造または FX の割引・受渡・三角経路の指定が矛盾する場合。
    pub fn with_fx(
        as_of: Instant,
        curves: impl IntoIterator<Item = CurveSpec>,
        fx: impl IntoIterator<Item = FxSpec>,
    ) -> Result<Self> {
        Self::with_fx_and_lgm(as_of, curves, fx, [])
    }

    /// FX と通貨別の LGM のパラメータを含む市場構造を作る。
    ///
    /// LGM は時刻の区切りの個数に 1 を足した区間ごとに、平均回帰とボラを持つ。
    /// θ では平均回帰を先に並べ、同じ区間順にボラを並べる。
    ///
    /// # Errors
    ///
    /// 市場構造の指定が不正か、LGM の区切りが評価時点以後で昇順でない場合。
    pub fn with_fx_and_lgm(
        as_of: Instant,
        curves: impl IntoIterator<Item = CurveSpec>,
        fx: impl IntoIterator<Item = FxSpec>,
        lgm: impl IntoIterator<Item = (Currency, Vec<Instant>)>,
    ) -> Result<Self> {
        let mut specs: Vec<_> = curves.into_iter().collect();
        // θ のブロック順を入力順から独立させるため、曲線 ID の辞書順に並べる。
        // 重複は並べた後の挿入で検出する。
        specs.sort_by(|a, b| a.id().cmp(b.id()));
        let fx: Vec<_> = fx.into_iter().collect();
        let (layout, grids) = build_layout(as_of, &specs, &fx, lgm)?;
        let mut resolver = CurveResolver {
            specs: &specs,
            layout: &layout,
            as_of,
            resolved: BTreeMap::new(),
            visiting: BTreeSet::new(),
        };
        let expressions = (0..specs.len())
            .map(|position| resolver.resolve(position))
            .collect::<Result<Vec<_>>>()?;
        let curves = specs
            .into_iter()
            .zip(expressions)
            .map(|(spec, expr)| {
                let mut blocks = BTreeSet::new();
                expr.collect_blocks(&mut blocks);
                Curve {
                    spec,
                    expr,
                    blocks: blocks.into_iter().collect(),
                    forecast_measure: None,
                }
            })
            .collect();
        let mut structure = Self {
            as_of,
            curves,
            grids,
            layout,
            fx: BTreeMap::new(),
        };
        structure.fx = fx::resolve(&structure, fx)?;
        Ok(structure)
    }

    /// 予想曲線ごとに、その予想を定める担保の割引曲線（予想の測度）を記録する。
    ///
    /// 予想曲線は指標ごとに一本で、その予想は校正商品の担保の下で定まる。
    /// 同じ担保で割り引く支払いでは予想と割引の担保が一致するので、曲線の ID が違っても
    /// `CrossCollateralProjection` の近似は要らない。記録の無い予想曲線は、
    /// 自身を割引にも使う場合だけを同じ担保とみなす。
    ///
    /// # Errors
    ///
    /// 曲線が無いか、予想曲線でないか、測度が同じ通貨の割引曲線でなければ `InvalidInput`。
    pub fn with_forecast_measures(
        mut self,
        measures: impl IntoIterator<Item = (CurveId, CurveId)>,
    ) -> Result<Self> {
        for (projection, discount) in measures {
            let position = self.role_position(&projection, CurveRole::Projection)?;
            let measure = self.role_position(&discount, CurveRole::Discount)?;
            let [projection_ccy, discount_ccy] =
                [position, measure].map(|index| self.curves[index].spec.currency());
            if projection_ccy != discount_ccy {
                return Err(Error::InvalidInput(format!(
                    "forecast measure {discount} of {projection} has currency {discount_ccy}, expected {projection_ccy}"
                )));
            }
            self.curves[position].forecast_measure = Some(measure);
        }
        Ok(self)
    }

    /// 予想曲線の予想を定める担保の割引曲線。記録が無ければ `None`。
    #[must_use]
    pub fn forecast_measure(&self, projection: &CurveId) -> Option<&CurveId> {
        let position = self.layout.position(projection)?;
        self.curves[position]
            .forecast_measure
            .map(|measure| self.curves[measure].spec.id())
    }

    fn role_position(&self, id: &CurveId, role: CurveRole) -> Result<usize> {
        let position = self
            .layout
            .position(id)
            .ok_or_else(|| Error::InvalidInput(format!("unknown curve {id}")))?;
        if self.curves[position].spec.has_role(role) {
            Ok(position)
        } else {
            Err(Error::InvalidInput(format!(
                "curve {id} has no {role:?} role"
            )))
        }
    }

    /// モデル時間の原点となる評価の瞬間。
    #[must_use]
    pub const fn as_of(&self) -> Instant {
        self.as_of
    }

    /// θ の配置。
    #[must_use]
    pub const fn layout(&self) -> &ParamLayout {
        &self.layout
    }

    /// 指定した曲線の定義。
    #[must_use]
    pub fn curve(&self, id: &CurveId) -> Option<&CurveSpec> {
        self.layout
            .position(id)
            .map(|position| &self.curves[position].spec)
    }

    /// 指定した曲線のピラーを評価時点からのモデル時間で返す。
    #[must_use]
    pub fn pillar_times(&self, id: &CurveId) -> Option<&[ModelTime]> {
        self.layout.position(id)?;
        Some(self.block_times(&BlockId::nodes(id.clone())).unwrap_or(&[]))
    }

    /// 指定したブロックのピラーまたはイベントのモデル時間。
    #[must_use]
    pub fn block_times(&self, block: &BlockId) -> Option<&[ModelTime]> {
        self.layout
            .block_position(block)
            .map(|position| self.grids[position].times.as_ref())
    }

    /// 通貨別 LGM の平均回帰とボラを持つ θ ブロックの ID。
    ///
    /// # Errors
    ///
    /// 通貨コードでブロック ID を構築できない場合。
    pub fn lgm_block(ccy: Currency) -> Result<BlockId> {
        BlockId::new(format!("LGM/{ccy}").parse()?, "marks")
    }
}

fn build_layout(
    as_of: Instant,
    specs: &[CurveSpec],
    fx: &[FxSpec],
    lgm: impl IntoIterator<Item = (Currency, Vec<Instant>)>,
) -> Result<(ParamLayout, Vec<BlockGrid>)> {
    let mut curve_positions = BTreeMap::new();
    let mut blocks = BTreeMap::new();
    for (position, spec) in specs.iter().enumerate() {
        if curve_positions
            .insert(spec.id().clone(), position)
            .is_some()
        {
            return Err(Error::InvalidInput(format!(
                "duplicate curve: {}",
                spec.id()
            )));
        }
        for (block, part) in spec.collect_blocks()? {
            let dates = Box::<[Instant]>::from(part.dates());
            let times: Box<_> = dates.iter().map(|&at| model_time(as_of, at)).collect();
            if times.first().is_some_and(|time| time.years() <= 0.0)
                || !times.is_sorted_by(|a, b| a < b)
            {
                return Err(Error::InvalidInput(format!(
                    "curve {} needs distinct pillar times after as-of for {}",
                    spec.id(),
                    block.part()
                )));
            }
            blocks.insert(block, (BlockGrid { dates, times }, part.len()));
        }
    }
    for spec in fx {
        if let Some(count) = spec.mark_count() {
            let block = FxSpec::marks_block(spec.pair)?;
            let grid = BlockGrid {
                dates: Box::new([]),
                times: Box::new([]),
            };
            if blocks.insert(block, (grid, count)).is_some() {
                return Err(Error::InvalidData(format!(
                    "duplicate FX marks for {}",
                    spec.pair
                )));
            }
        }
    }
    for (ccy, dates) in lgm {
        let count = (dates.len() + 1) * 2;
        let times: Box<_> = dates.iter().map(|&at| model_time(as_of, at)).collect();
        if times.first().is_some_and(|time| time.years() <= 0.0)
            || !times.is_sorted_by(|a, b| a < b)
        {
            return Err(Error::InvalidData(format!(
                "invalid LGM regime boundaries for {ccy}"
            )));
        }
        let grid = BlockGrid {
            dates: dates.into_boxed_slice(),
            times,
        };
        if blocks
            .insert(MarketStructure::lgm_block(ccy)?, (grid, count))
            .is_some()
        {
            return Err(Error::InvalidData(format!("duplicate LGM block for {ccy}")));
        }
    }
    let mut block_positions = BTreeMap::new();
    let mut grids = Vec::with_capacity(blocks.len());
    let mut ranges = Vec::with_capacity(blocks.len());
    let mut ids = Vec::with_capacity(blocks.len());
    let mut len = 0usize;
    for (block, (grid, count)) in blocks {
        let end = len
            .checked_add(count)
            .ok_or_else(|| Error::InvalidInput("theta layout is too large".into()))?;
        block_positions.insert(block.clone(), ids.len());
        ids.push(block);
        grids.push(grid);
        ranges.push(len..end);
        len = end;
    }
    Ok((
        ParamLayout {
            curve_positions,
            block_positions,
            ids,
            ranges,
            len,
        },
        grids,
    ))
}

struct CurveResolver<'a> {
    specs: &'a [CurveSpec],
    layout: &'a ParamLayout,
    as_of: Instant,
    resolved: BTreeMap<usize, Arc<ResolvedExpr>>,
    visiting: BTreeSet<usize>,
}

impl CurveResolver<'_> {
    fn resolve(&mut self, position: usize) -> Result<Arc<ResolvedExpr>> {
        if let Some(expr) = self.resolved.get(&position) {
            return Ok(Arc::clone(expr));
        }
        let spec = &self.specs[position];
        if !self.visiting.insert(position) {
            return Err(Error::InvalidData(format!(
                "cyclic curve reference: {}",
                spec.id()
            )));
        }
        let expr = self.resolve_part(spec, &spec.expr)?;
        self.visiting.remove(&position);
        self.resolved.insert(position, Arc::clone(&expr));
        Ok(expr)
    }

    fn position(&self, id: &CurveId) -> Result<usize> {
        self.layout
            .position(id)
            .ok_or_else(|| Error::InvalidData(format!("unknown curve reference: {id}")))
    }

    fn resolve_part(&mut self, spec: &CurveSpec, expr: &CurveExpr) -> Result<Arc<ResolvedExpr>> {
        let block = |part| self.layout.block_positions[&BlockId::for_part(spec.id().clone(), part)];
        let resolved = match expr {
            CurveExpr::Nodes {
                part,
                extrapolation,
                ..
            } => ResolvedExpr::Nodes {
                block: block(part),
                extrapolation: *extrapolation,
            },
            CurveExpr::Steps { part, .. } => ResolvedExpr::Steps(block(part)),
            CurveExpr::Jumps { part, .. } => ResolvedExpr::Jumps(block(part)),
            CurveExpr::Flat(rate) => ResolvedExpr::Flat(*rate),
            CurveExpr::Sum(parts) => ResolvedExpr::Sum(
                parts
                    .iter()
                    .map(|part| self.resolve_part(spec, part))
                    .collect::<Result<_>>()?,
            ),
            CurveExpr::Scale { factor, body } => {
                ResolvedExpr::Scale(*factor, self.resolve_part(spec, body)?)
            }
            CurveExpr::MaxForward { parts, grid } => self.resolve_max_forward(spec, parts, grid)?,
            CurveExpr::Stitch { short, long, at } => {
                if *at <= self.as_of {
                    return Err(Error::InvalidData(format!(
                        "stitch on {} must be after as-of",
                        spec.id()
                    )));
                }
                ResolvedExpr::Stitch {
                    short: self.resolve_part(spec, short)?,
                    long: self.resolve_part(spec, long)?,
                    at: *at,
                }
            }
            CurveExpr::Ref(id) => {
                let position = self.position(id)?;
                let source = &self.specs[position];
                if source.currency() != spec.currency()
                    || (spec.has_role(CurveRole::Discount) && !source.has_role(CurveRole::Discount))
                {
                    return Err(Error::InvalidData(format!(
                        "invalid curve reference from {} to {id}",
                        spec.id()
                    )));
                }
                return self.resolve(position);
            }
            CurveExpr::ForwardRef(id) => {
                let position = self.position(id)?;
                if !self.specs[position].has_role(CurveRole::Projection) {
                    return Err(Error::InvalidData(format!(
                        "forward reference from {} to {id} requires a projection curve",
                        spec.id()
                    )));
                }
                return self.resolve(position);
            }
            CurveExpr::FxImplied {
                part,
                pair,
                base,
                extrapolation,
                ..
            } => {
                let position = self.position(base)?;
                let source = &self.specs[position];
                let sign = match (spec.currency(), source.currency()) {
                    (ccy, other) if ccy == pair.base() && other == pair.quote() => 1.0,
                    (ccy, other) if ccy == pair.quote() && other == pair.base() => -1.0,
                    _ => {
                        return Err(Error::InvalidData(format!(
                            "currencies of {} and {base} do not match {pair}",
                            spec.id()
                        )));
                    }
                };
                if !spec.has_role(CurveRole::Discount) || !source.has_role(CurveRole::Discount) {
                    return Err(Error::InvalidData(format!(
                        "FX-implied curve {} requires discount curves",
                        spec.id()
                    )));
                }
                let basis = ResolvedExpr::Nodes {
                    block: block(part),
                    extrapolation: *extrapolation,
                };
                ResolvedExpr::Sum(vec![
                    self.resolve(position)?,
                    Arc::new(ResolvedExpr::Scale(sign, Arc::new(basis))),
                ])
            }
        };
        Ok(Arc::new(resolved))
    }

    fn resolve_max_forward(
        &mut self,
        spec: &CurveSpec,
        parts: &[CurveExpr],
        grid: &[Instant],
    ) -> Result<ResolvedExpr> {
        // 候補と格子の有無、UTC 時刻の順序は CurveSpec の構築時に検査済み。
        if grid[0] <= self.as_of
            || !grid
                .iter()
                .map(|&at| model_time(self.as_of, at))
                .is_sorted_by(|a, b| a < b)
        {
            return Err(Error::InvalidInput(format!(
                "max forward on {} needs distinct grid times after as-of",
                spec.id()
            )));
        }
        Ok(ResolvedExpr::MaxForward {
            parts: parts
                .iter()
                .map(|part| self.resolve_part(spec, part))
                .collect::<Result<_>>()?,
            grid: grid.into(),
        })
    }
}
`````

`````
crates/rhodium-market/src/view.rs
//! ブランド付きの市場読み出しと通貨スコープ。

use std::{
    marker::PhantomData,
    ops::{Add, Neg, Sub},
    sync::Arc,
};

use rhodium_core::{
    Error, Result,
    brand::{Brand, with_brand},
    currency::Currency,
    fx::{CcyPair, FxRate, Today},
    money::Money,
    quantity::{DiscountFactor, YearFraction},
    time::Instant,
};
use scioxide_ad::{Mask, Real};

use crate::{
    curve::{CurveId, CurvePoint, CurveRole, CurveSpec},
    fx::FxRef,
    market::{Market, ThetaBlock},
    structure::{BlockId, MarketStructure},
};

/// 割引用の曲線の役割。
#[derive(Debug, Clone, Copy)]
pub struct Discount;

/// 金利予想用の曲線の役割。
#[derive(Debug, Clone, Copy)]
pub struct Projection;

/// 市場構造に解決した、役割付きの曲線の位置。
#[derive(Debug, Clone)]
pub struct ResolvedCurve<K> {
    structure: Arc<MarketStructure>,
    position: usize,
    _role: PhantomData<K>,
}

impl<K> ResolvedCurve<K> {
    /// 解決した曲線の定義。
    #[must_use]
    pub fn spec(&self) -> &CurveSpec {
        &self.structure.curves[self.position].spec
    }

    /// 解決した曲線上の点を準備する。
    ///
    /// # Errors
    /// 評価時点より前、または曲線が外挿を許さない範囲なら `InvalidInput`。
    pub fn prepare_point(&self, at: Instant) -> Result<CurvePoint> {
        self.structure.prepare_point_at(self.position, at)
    }
}

impl ResolvedCurve<Projection> {
    /// この予想が `discount` と同じ担保の下で定まるか。
    ///
    /// 構造に予想の測度（[`MarketStructure::with_forecast_measures`]）があればそれと比べ、
    /// 無ければ同じ曲線を予想と割引に使う場合だけを同じ担保とみなす。
    /// 違う構造の曲線どうしは比べられないので、常に違う担保として扱う。
    #[must_use]
    pub fn is_forecast_under(&self, discount: &ResolvedCurve<Discount>) -> bool {
        Arc::ptr_eq(&self.structure, &discount.structure)
            && self.structure.curves[self.position]
                .forecast_measure
                .unwrap_or(self.position)
                == discount.position
    }
}

impl MarketStructure {
    /// 割引曲線の位置を解決する。
    ///
    /// # Errors
    /// 曲線が無ければ `MissingData`、役割が違えば `InvalidInput`。
    pub fn resolve_discount(self: &Arc<Self>, id: &CurveId) -> Result<ResolvedCurve<Discount>> {
        self.resolve_curve(id, CurveRole::Discount)
    }

    /// 予想曲線の位置を解決する。
    ///
    /// # Errors
    /// 曲線が無ければ `MissingData`、役割が違えば `InvalidInput`。
    pub fn resolve_projection(self: &Arc<Self>, id: &CurveId) -> Result<ResolvedCurve<Projection>> {
        self.resolve_curve(id, CurveRole::Projection)
    }

    fn resolve_curve<K>(
        self: &Arc<Self>,
        id: &CurveId,
        role: CurveRole,
    ) -> Result<ResolvedCurve<K>> {
        Ok(ResolvedCurve {
            structure: Arc::clone(self),
            position: curve_position(self, id, role)?,
            _role: PhantomData,
        })
    }
}

/// 一つの市場ビューの中で使える、役割付きの曲線ハンドル。
#[derive(Debug, Clone, Copy)]
pub struct CurveRef<'v, K> {
    position: usize,
    _brand: Brand<'v>,
    _role: PhantomData<K>,
}

/// 割引曲線と準備済みの点を結び付け、計算済みの割引因子を保持する。
#[derive(Debug, Clone, Copy)]
pub struct CurveAt<'v, K, R = f64> {
    value: DiscountFactor<R>,
    curve: CurveRef<'v, K>,
}

/// 予想曲線、期間の両端と日数係数を結び付け、単利の予想値を保持する。
#[derive(Debug, Clone, Copy)]
pub struct ForwardAt<'v, R = f64> {
    value: R,
    _curve: CurveRef<'v, Projection>,
}

/// 予想曲線と両端の点を結合した、複利の成長率。
#[derive(Debug, Clone, Copy)]
pub struct GrowthAt<'v, R = f64> {
    value: R,
    // 成長率の自然対数 ln P(開始) − ln P(終了)。超過成長 G − 1 を exp_m1 で求めるために持つ。
    log: R,
    _curve: CurveRef<'v, Projection>,
}

type ForwardPoints = (CurvePoint, CurvePoint, YearFraction);
type GrowthPoints = (CurvePoint, CurvePoint);

#[derive(Debug)]
struct CurveSlot<T> {
    position: usize,
    points: Vec<T>,
}

/// 準備済みの曲線スロットと点。同じ構造から作る市場の間で使い回せる。
#[derive(Debug)]
pub struct Slots {
    structure: Arc<MarketStructure>,
    discounts: Vec<CurveSlot<CurvePoint>>,
    projections: Vec<CurveSlot<ForwardPoints>>,
    growth: Vec<CurveSlot<GrowthPoints>>,
}

impl Slots {
    /// 割引に使う準備済みの曲線と支払時刻。
    #[must_use]
    pub fn discount_point(&self, slot: usize, point: usize) -> Option<(&CurveSpec, Instant)> {
        let slot = self.discounts.get(slot)?;
        Some((
            &self.structure.curves[slot.position].spec,
            slot.points.get(point)?.at(),
        ))
    }

    /// DF 比に使う準備済みの予想曲線と期間の始点・終点。
    #[must_use]
    pub fn growth_interval(
        &self,
        slot: usize,
        point: usize,
    ) -> Option<(&CurveSpec, Instant, Instant)> {
        let slot = self.growth.get(slot)?;
        let (start, end) = slot.points.get(point)?;
        Some((
            &self.structure.curves[slot.position].spec,
            start.at(),
            end.at(),
        ))
    }

    /// 役割別に曲線 ID を解決する。予想の点は開始、終了、日数係数の順。
    ///
    /// 入力順でスロットと各スロット内の点を保持する。曲線のハンドルだけが必要な
    /// スロットには空の点の列を渡せる。
    ///
    /// # Errors
    ///
    /// 曲線が無ければ `MissingData`、役割が違えば `InvalidInput`。
    pub fn new(
        structure: Arc<MarketStructure>,
        discounts: impl IntoIterator<Item = (CurveId, Vec<CurvePoint>)>,
        projections: impl IntoIterator<Item = (CurveId, Vec<ForwardPoints>)>,
    ) -> Result<Self> {
        let discounts = resolve_slots(&structure, CurveRole::Discount, discounts)?;
        let projections = resolve_slots(&structure, CurveRole::Projection, projections)?;
        Ok(Self {
            structure,
            discounts,
            projections,
            growth: Vec::new(),
        })
    }

    /// 解決した曲線と準備済みの点から、入力順にスロットを作る。
    ///
    /// # Errors
    /// 解決した曲線が別の構造に属する場合は `InvalidInput`。
    /// 点が対応する曲線に属するかは [`MarketView::bind_slots`] が検査する。
    pub fn from_resolved(
        structure: Arc<MarketStructure>,
        discounts: impl IntoIterator<Item = (ResolvedCurve<Discount>, Vec<CurvePoint>)>,
        projections: impl IntoIterator<Item = (ResolvedCurve<Projection>, Vec<ForwardPoints>)>,
    ) -> Result<Self> {
        let discounts = resolved_slots(&structure, discounts)?;
        let projections = resolved_slots(&structure, projections)?;
        Ok(Self {
            structure,
            discounts,
            projections,
            growth: Vec::new(),
        })
    }

    /// 予想曲線上の DF 比を読むスロットを、渡した順で追加する。
    ///
    /// # Errors
    /// 解決済みの曲線が別の構造に属するときは `InvalidInput`。
    pub fn with_projection_growth(
        mut self,
        growth: impl IntoIterator<Item = (ResolvedCurve<Projection>, Vec<GrowthPoints>)>,
    ) -> Result<Self> {
        self.growth.extend(resolved_slots(&self.structure, growth)?);
        Ok(self)
    }
}

fn resolved_slots<K, T>(
    structure: &Arc<MarketStructure>,
    specs: impl IntoIterator<Item = (ResolvedCurve<K>, Vec<T>)>,
) -> Result<Vec<CurveSlot<T>>> {
    specs
        .into_iter()
        .map(|(curve, points)| {
            if !Arc::ptr_eq(structure, &curve.structure) {
                return Err(Error::InvalidInput(
                    "resolved curve belongs to a different market structure".into(),
                ));
            }
            Ok(CurveSlot {
                position: curve.position,
                points,
            })
        })
        .collect()
}

fn resolve_slots<T>(
    structure: &MarketStructure,
    role: CurveRole,
    specs: impl IntoIterator<Item = (CurveId, Vec<T>)>,
) -> Result<Vec<CurveSlot<T>>> {
    specs
        .into_iter()
        .map(|(id, points)| {
            let position = curve_position(structure, &id, role)?;
            Ok(CurveSlot { position, points })
        })
        .collect()
}

fn curve_position(structure: &MarketStructure, id: &CurveId, role: CurveRole) -> Result<usize> {
    let position = structure
        .layout()
        .position(id)
        .ok_or_else(|| Error::MissingData(format!("unknown curve: {id}")))?;
    if !structure.curves[position].spec.has_role(role) {
        return Err(Error::InvalidInput(format!(
            "curve {id} does not have role {role:?}"
        )));
    }
    Ok(position)
}

/// 一つのビューに結合した曲線ハンドルと、その曲線の点。
#[derive(Debug)]
pub struct BoundCurve<'v, K, T> {
    curve: CurveRef<'v, K>,
    points: Vec<T>,
}

impl<'v, K: Copy, T> BoundCurve<'v, K, T> {
    /// 検証済みの曲線ハンドル。
    #[must_use]
    pub const fn curve(&self) -> CurveRef<'v, K> {
        self.curve
    }

    /// 結合した点を、準備時に指定した順で返す。
    #[must_use]
    pub fn points(&self) -> &[T] {
        &self.points
    }
}

/// 一つのビューに結合した役割別のスロット。
#[derive(Debug)]
pub struct BoundSlots<'v, R = f64> {
    discounts: Vec<BoundCurve<'v, Discount, CurveAt<'v, Discount, R>>>,
    projections: Vec<BoundCurve<'v, Projection, ForwardAt<'v, R>>>,
    growth: Vec<BoundCurve<'v, Projection, GrowthAt<'v, R>>>,
}

impl<'v, R> BoundSlots<'v, R> {
    /// 割引曲線のスロットを準備時に指定した順で返す。
    #[must_use]
    pub fn discounts(&self) -> &[BoundCurve<'v, Discount, CurveAt<'v, Discount, R>>] {
        &self.discounts
    }

    /// 予想曲線のスロットを準備時に指定した順で返す。
    #[must_use]
    pub fn projections(&self) -> &[BoundCurve<'v, Projection, ForwardAt<'v, R>>] {
        &self.projections
    }

    /// DF 比のスロットを、準備時に指定した順で返す。
    #[must_use]
    pub fn projection_growth(&self) -> &[BoundCurve<'v, Projection, GrowthAt<'v, R>>] {
        &self.growth
    }
}

/// 一つの市場の役割付きハンドルと点を作り、読み出すビュー。
#[derive(Debug)]
pub struct MarketView<'v, R = f64> {
    market: Market<R>,
    brand: Brand<'v>,
}

impl<R: Real> Market<R> {
    /// 新しいブランドでビューを開く。ハンドルや結合済みの点はビューの外へ出せない。
    pub fn with_view<T>(&self, f: impl for<'v> FnOnce(MarketView<'v, R>) -> T) -> T {
        with_brand(|brand| {
            f(MarketView {
                market: self.clone(),
                brand,
            })
        })
    }
}

impl<'v, R: Real> MarketView<'v, R> {
    /// 同じ構造で準備した LGM の θ ブロックを読む。
    ///
    /// # Errors
    ///
    /// 構造または LGM の指定が異なるか、ブロックを利用できない場合。
    pub fn lgm_parameters(
        &self,
        structure: &Arc<MarketStructure>,
        ccy: Currency,
        block: &BlockId,
    ) -> Result<&[R]> {
        if !Arc::ptr_eq(structure, self.market.structure())
            || MarketStructure::lgm_block(ccy)? != *block
        {
            return Err(Error::InvalidInput(
                "LGM parameters belong to another structure or currency".into(),
            ));
        }
        let position = structure
            .layout()
            .block_position(block)
            .filter(|&position| self.market.available_at(position))
            .ok_or_else(|| Error::MissingData(format!("missing LGM parameters for {ccy}")))?;
        Ok(match self.market.theta_block(position) {
            ThetaBlock::Values(values) => values,
            // 既存のスライス API を使い、passive なブロックを明示して読んだ場合に平坦化する。
            ThetaBlock::Constants(_) => {
                &self.market.theta()[structure.layout().ranges[position].clone()]
            }
        })
    }
    /// 準備済みの FX フォワードをこのビューの市場で評価する。
    ///
    /// # Errors
    ///
    /// 市場構造が異なるか、必要な値を表現できない場合。
    pub fn fx_forward_value(&self, plan: &crate::fx::FxForwardPlan) -> Result<R> {
        plan.rate(&self.market)
    }
    /// ペアと担保を解決し、共通の評価時点 FX を持つハンドルを作る。
    ///
    /// # Errors
    ///
    /// ペア、担保、起点の市場データが利用できないか、結果を表現できない場合。
    pub fn fx(&self, pair: CcyPair, collateral: &CurveId) -> Result<FxRef<'v, R>> {
        FxRef::new(&self.market, pair, collateral, self.brand)
    }

    /// 評価時点の通貨換算レート。ペアの向きに従って返す。
    ///
    /// # Errors
    ///
    /// ペアまたは起点の市場データが利用できないか、結果を表現できない場合。
    pub fn fx_today(&self, pair: CcyPair) -> Result<FxRate<Today, R>> {
        self.market.fx_today(pair)
    }

    /// ID から割引用の曲線ハンドルを得る。
    ///
    /// # Errors
    ///
    /// 曲線が無いか利用できなければ `MissingData`、役割が違えば `InvalidInput`。
    pub fn discount_curve(&self, id: &CurveId) -> Result<CurveRef<'v, Discount>> {
        self.curve_ref(id, CurveRole::Discount)
    }

    /// ID から予想用の曲線ハンドルを得る。
    ///
    /// # Errors
    ///
    /// 曲線が無いか利用できなければ `MissingData`、役割が違えば `InvalidInput`。
    pub fn projection_curve(&self, id: &CurveId) -> Result<CurveRef<'v, Projection>> {
        self.curve_ref(id, CurveRole::Projection)
    }

    fn curve_ref<K>(&self, id: &CurveId, role: CurveRole) -> Result<CurveRef<'v, K>> {
        let position = curve_position(self.market.structure(), id, role)?;
        self.available_ref(position)
    }

    fn available_ref<K>(&self, position: usize) -> Result<CurveRef<'v, K>> {
        self.market.require_curve(position)?;
        Ok(self.resolved_ref(position))
    }

    fn resolved_ref<K>(&self, position: usize) -> CurveRef<'v, K> {
        CurveRef {
            position,
            _brand: self.brand,
            _role: PhantomData,
        }
    }

    /// 準備済みのスロットを現在の市場に結び付け、読み出す値を計算する。
    ///
    /// # Errors
    ///
    /// 構造や点の曲線が違う、期間の順が逆か同時刻、日数係数が正でなければ
    /// `InvalidInput`。ブロックが利用できなければ `MissingData`。
    /// DF やフォワードを表現できなければ `Numerical`。
    pub fn bind_slots(&self, slots: &Slots) -> Result<BoundSlots<'v, R>> {
        if !Arc::ptr_eq(self.market.structure(), &slots.structure) {
            return Err(Error::InvalidInput(
                "slots belong to a different market structure".into(),
            ));
        }

        let discounts = slots
            .discounts
            .iter()
            .map(|slot| {
                if slot
                    .points
                    .iter()
                    .all(|point| point.at() == slots.structure.as_of())
                {
                    self.market.require_curve(slot.position)?;
                }
                let curve = self.resolved_ref(slot.position);
                let points = slot
                    .points
                    .iter()
                    .map(|point| {
                        let value = discount_at(&self.market, slot.position, point)?;
                        Ok(CurveAt { value, curve })
                    })
                    .collect::<Result<_>>()?;
                Ok(BoundCurve { curve, points })
            })
            .collect::<Result<_>>()?;
        let projections = slots
            .projections
            .iter()
            .map(|slot| {
                if slot.points.is_empty()
                    && !slots
                        .growth
                        .iter()
                        .any(|growth| growth.position == slot.position && !growth.points.is_empty())
                {
                    self.market.require_curve(slot.position)?;
                }
                let curve = self.resolved_ref(slot.position);
                let points = slot
                    .points
                    .iter()
                    .map(|(start, end, tau)| {
                        let value = forward_at(&self.market, slot.position, start, end, tau)?;
                        Ok(ForwardAt {
                            value,
                            _curve: curve,
                        })
                    })
                    .collect::<Result<_>>()?;
                Ok(BoundCurve { curve, points })
            })
            .collect::<Result<_>>()?;
        let growth = slots
            .growth
            .iter()
            .map(|slot| {
                if slot.points.is_empty()
                    && !slots.projections.iter().any(|projection| {
                        projection.position == slot.position && !projection.points.is_empty()
                    })
                {
                    self.market.require_curve(slot.position)?;
                }
                let curve = self.resolved_ref(slot.position);
                let points = slot
                    .points
                    .iter()
                    .map(|(start, end)| {
                        let delta = projection_log_ratio(&self.market, slot.position, start, end)?;
                        let value = delta.exp();
                        if !delta.is_finite().all()
                            || !value.is_finite().all()
                            || !value.gt(R::zero()).all()
                        {
                            let id = self.market.structure().curves[slot.position].spec.id();
                            return Err(Error::Numerical(format!(
                                "nonpositive or non-finite projection growth on curve {id}"
                            )));
                        }
                        Ok(GrowthAt {
                            value,
                            log: delta,
                            _curve: curve,
                        })
                    })
                    .collect::<Result<_>>()?;
                Ok(BoundCurve { curve, points })
            })
            .collect::<Result<_>>()?;
        Ok(BoundSlots {
            discounts,
            projections,
            growth,
        })
    }

    /// 結合済みの割引因子を読む。
    #[must_use]
    pub const fn df(&self, at: CurveAt<'v, Discount, R>) -> DiscountFactor<R> {
        at.value
    }

    /// 結合済みの単利フォワードを読む。
    #[must_use]
    pub const fn forward(&self, at: ForwardAt<'v, R>) -> R {
        at.value
    }

    /// 結合済みの予想曲線の DF 比 P(開始) / P(終了) を読む。
    #[must_use]
    pub const fn growth(&self, at: GrowthAt<'v, R>) -> R {
        at.value
    }

    /// 結合済みの予想曲線の DF 比の自然対数 ln P(開始) − ln P(終了) を読む。
    ///
    /// 短い期間の DF 比は 1 に近く、[`Self::growth`] から 1 を引くと ulp(1) の丸めで
    /// 超過の下位の桁を失う。超過成長は、この対数から `exp_m1` で求める。
    #[must_use]
    pub const fn log_growth(&self, at: GrowthAt<'v, R>) -> R {
        at.log
    }

    /// 指定した通貨の支払いを割り引くスコープを開く。
    ///
    /// # Errors
    ///
    /// 割引曲線の通貨がスコープの通貨と違えば `InvalidInput`。
    /// クロージャのエラーはそのまま返す。
    pub fn ccy_scope<T>(
        &self,
        ccy: Currency,
        discount: CurveRef<'v, Discount>,
        f: impl for<'c> FnOnce(CcyScope<'c, 'v, R>) -> Result<T>,
    ) -> Result<T> {
        let spec = &self.market.structure().curves[discount.position].spec;
        if spec.currency() != ccy {
            return Err(Error::InvalidInput(format!(
                "curve {} is in {}, not {ccy}",
                spec.id(),
                spec.currency()
            )));
        }
        with_brand(|brand| {
            f(CcyScope {
                ccy,
                market: self.market.clone(),
                discount,
                brand,
            })
        })
    }
}

fn discount_at<R: Real>(
    market: &Market<R>,
    curve: usize,
    point: &CurvePoint,
) -> Result<DiscountFactor<R>> {
    point.check(market.structure(), curve)?;
    // P = exp(L)。有限な θ でも exp は 0 や ∞ になり得るので、ここで一度だけ検査し、
    // 評価ループの読み出しから失敗を除く。
    market.require_blocks(point.block_positions())?;
    let df = point.log_df(market).exp();
    DiscountFactor::new(df).map_err(|_| {
        let id = market.structure().curves[curve].spec.id();
        Error::Numerical(format!(
            "discount factor on curve {id} is not finite and positive"
        ))
    })
}

fn forward_at<R: Real>(
    market: &Market<R>,
    curve: usize,
    start: &CurvePoint,
    end: &CurvePoint,
    tau: &YearFraction,
) -> Result<R> {
    if tau.value() <= 0.0 {
        return Err(Error::InvalidInput(
            "forward requires a positive year fraction".into(),
        ));
    }
    let delta = projection_log_ratio(market, curve, start, end)?;
    // exp_m1 は短い期間や小さい金利でも割引因子の比と 1 の差を保つ。
    let value = delta.exp_m1() / tau.value();
    if !delta.is_finite().all() || !value.is_finite().all() {
        let id = market.structure().curves[curve].spec.id();
        return Err(Error::Numerical(format!(
            "non-finite forward on curve {id}"
        )));
    }
    Ok(value)
}

fn projection_log_ratio<R: Real>(
    market: &Market<R>,
    curve: usize,
    start: &CurvePoint,
    end: &CurvePoint,
) -> Result<R> {
    let structure = market.structure();
    start.check(structure, curve)?;
    end.check(structure, curve)?;
    if start.at() >= end.at() {
        return Err(Error::InvalidInput(
            "projection requires increasing instants".into(),
        ));
    }
    market.require_blocks(start.block_positions().chain(end.block_positions()))?;
    Ok(start.log_df(market) - end.log_df(market))
}

/// 通貨スコープの割引曲線に結合した支払点。
#[derive(Debug, Clone, Copy)]
pub struct PayAt<'c, R = f64> {
    df: DiscountFactor<R>,
    _brand: Brand<'c>,
}

/// 指定した割引曲線による一つの通貨の現在価値を扱う。
#[derive(Debug)]
pub struct CcyScope<'c, 'v, R = f64> {
    ccy: Currency,
    market: Market<R>,
    discount: CurveRef<'v, Discount>,
    brand: Brand<'c>,
}

impl<'c, 'v, R: Real> CcyScope<'c, 'v, R> {
    /// 支払点をスコープの割引曲線に結び付ける。
    ///
    /// # Errors
    ///
    /// 点が別の曲線または別の構造のものなら `InvalidInput`。
    /// DF を正の有限値で表現できなければ `Numerical`。
    pub fn at(&self, pay: &CurvePoint) -> Result<PayAt<'c, R>> {
        let df = discount_at(&self.market, self.discount.position, pay)?;
        Ok(PayAt {
            df,
            _brand: self.brand,
        })
    }

    /// ビューで結合した割引点の値を、この通貨スコープに移す。
    ///
    /// # Errors
    /// スコープと別の割引曲線に結合した点なら `InvalidInput`。
    pub fn bind_pay(&self, pay: CurveAt<'v, Discount, R>) -> Result<PayAt<'c, R>> {
        if pay.curve.position != self.discount.position {
            return Err(Error::InvalidInput(
                "payment point belongs to a different discount curve".into(),
            ));
        }
        Ok(PayAt {
            df: pay.value,
            _brand: self.brand,
        })
    }

    /// 一つの支払いを結合済みの点で割り引く。
    #[must_use]
    pub fn pv(&self, amount: R, pay: PayAt<'c, R>) -> Amount<'c, R> {
        self.present(amount * pay.df.value())
    }

    /// 呼び出し側で計算済みの現在価値を、スコープに結び付ける。
    #[must_use]
    pub const fn present(&self, value: R) -> Amount<'c, R> {
        Amount {
            value,
            brand: self.brand,
        }
    }

    /// 別の通貨スコープの現在価値を評価時点の FX で換算する。
    ///
    /// # Errors
    ///
    /// ペアのレートが利用できないか、換算後の金額を表現できない場合。
    pub fn convert_from<'e>(
        &self,
        from: &CcyScope<'e, 'v, R>,
        amount: Amount<'e, R>,
    ) -> Result<Amount<'c, R>> {
        if from.ccy == self.ccy {
            return Ok(self.present(amount.value));
        }
        let pair = CcyPair::new(from.ccy, self.ccy)?;
        let converted = self
            .market
            .fx_today(pair)?
            .convert(from.into_money(amount))?;
        Ok(self.present(converted.amount()))
    }

    /// 現在価値を通貨タグ付きの金額として取り出す。
    #[must_use]
    pub const fn into_money(&self, amount: Amount<'c, R>) -> Money<R> {
        Money::new(self.ccy, amount.value)
    }
}

/// 同じ通貨スコープの内側で加減算できる現在価値。
#[derive(Debug, Clone, Copy)]
pub struct Amount<'c, R = f64> {
    value: R,
    brand: Brand<'c>,
}

impl<R: Real> Amount<'_, R> {
    /// 現在価値を倍率で拡大する。
    #[must_use]
    pub fn scale(self, factor: R) -> Self {
        Self {
            value: self.value * factor,
            brand: self.brand,
        }
    }
}

impl<R: Real> Add for Amount<'_, R> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            value: self.value + rhs.value,
            brand: self.brand,
        }
    }
}

impl<R: Real> Sub for Amount<'_, R> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            value: self.value - rhs.value,
            brand: self.brand,
        }
    }
}

impl<R: Real> Neg for Amount<'_, R> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            value: -self.value,
            brand: self.brand,
        }
    }
}
`````

`````
crates/rhodium-market/tests/compile_fail.rs
//! 曲線の役割、ビュー、通貨スコープの型の境界を検証する。

use rhodium_testkit::compile_fail_pair;

#[test]
fn curve_ref_cross_view() {
    // 同じ市場でも、一つのビューで解決した曲線は別のビューへ渡せない。
    compile_fail_pair!("curve_ref_cross_view");
}

#[test]
fn discount_point_cross_view() {
    // 割引点は、その市場ビューの値であることも型で保持する。
    compile_fail_pair!("discount_point_cross_view");
}

#[test]
fn forward_point_cross_view() {
    // 予想期間の結合済みの値も、結合したビューでだけ使える。
    compile_fail_pair!("forward_point_cross_view");
}

#[test]
fn growth_point_cross_view() {
    // DF 比も、予想曲線に結合した市場ビューの中でだけ読める。
    compile_fail_pair!("growth_point_cross_view");
}

#[test]
fn projection_curve_as_discount() {
    // 一つの曲線が両方の役割を持っても、割引のスコープには割引ハンドルが要る。
    compile_fail_pair!("projection_curve_as_discount");
}

#[test]
fn projection_position_as_discount_slot() {
    // 構造に解決した位置も、曲線の役割が異なるスロットには入れられない。
    compile_fail_pair!("projection_position_as_discount_slot");
}

#[test]
fn read_with_unbound_point() {
    // 準備しただけの点は、曲線と結び付けてから読む。
    compile_fail_pair!("read_with_unbound_point");
}

#[test]
fn pay_point_cross_scope() {
    // 同じ曲線の点でも、別の通貨スコープで結合した値では割り引けない。
    compile_fail_pair!("pay_point_cross_scope");
}

#[test]
fn amount_cross_scope() {
    // 異なるスコープの現在価値を、そのまま足すことはできない。
    compile_fail_pair!("amount_cross_scope");
}

#[test]
fn amount_plus_undiscounted() {
    // 割り引いていない将来の額は R のままなので、現在価値の Amount に直接は足せない。
    compile_fail_pair!("amount_plus_undiscounted");
}

#[test]
fn fx_ref_cannot_escape_view() {
    // FX ハンドルは結合先のビューの外へ持ち出せない。
    compile_fail_pair!("fx_ref_cannot_escape_view");
}
`````

`````
crates/rhodium-market/tests/curve_reads.rs
//! 曲線の数値と、市場ビューに結合したときの入力境界を検証する。

use std::{collections::BTreeSet, sync::Arc};

use rhodium_core::{
    Error, Result,
    currency::Currency,
    day_count::DayCount,
    quantity::YearFraction,
    time::{Instant, model_time},
};
use rhodium_market::{
    curve::{CurveExpr, CurveId, CurveRole, CurveSpec, Extrapolation},
    market::Market,
    structure::{BlockId, MarketStructure},
    view::Slots,
};
use rhodium_testkit::assert_close;
use scioxide_ad::{Dual, Real, Tape};
use scioxide_core::{ExecutionContext, InvalidInputReason, SciOxideError, argument};
use scioxide_interp::{BSpline, ExtrapolationPolicy, Interpolant};

const ORIGIN: &str = "2026-01-01T00:00:00Z";
const FIRST: &str = "2027-01-01T00:00:00Z";
const SECOND: &str = "2028-01-01T00:00:00Z";
const THIRD: &str = "2028-12-31T00:00:00Z";
const REL: f64 = 1e-12;
const ABS: f64 = 1e-14;

fn structure(extrapolation: Extrapolation) -> Result<(Arc<MarketStructure>, CurveId)> {
    let id: CurveId = "JPY-TONA".parse()?;
    let spec = CurveSpec::new(
        id.clone(),
        Currency::JPY,
        [CurveRole::Discount, CurveRole::Projection],
        [FIRST.parse()?, SECOND.parse()?],
        extrapolation,
    )?;
    Ok((Arc::new(MarketStructure::new(ORIGIN.parse()?, [spec])?), id))
}

fn two_curves() -> Result<(Arc<MarketStructure>, CurveId, CurveId)> {
    let discount: CurveId = "USD-SOFR-DISCOUNT".parse()?;
    let projection: CurveId = "USD-SOFR-PROJECTION".parse()?;
    let pillars = [FIRST.parse()?];
    let discount_spec = CurveSpec::new(
        discount.clone(),
        Currency::USD,
        [CurveRole::Discount],
        pillars,
        Extrapolation::default(),
    )?;
    let projection_spec = CurveSpec::new(
        projection.clone(),
        Currency::USD,
        [CurveRole::Projection],
        pillars,
        Extrapolation::default(),
    )?;
    let market = MarketStructure::new(ORIGIN.parse()?, [discount_spec, projection_spec])?;
    Ok((Arc::new(market), discount, projection))
}

fn discount_slots(structure: &Arc<MarketStructure>, id: &CurveId, at: &str) -> Result<Slots> {
    let point = structure.prepare_point(id, at.parse()?)?;
    Slots::new(Arc::clone(structure), [(id.clone(), vec![point])], [])
}

fn df<R: Real>(market: &Market<R>, slots: &Slots) -> Result<R> {
    market.with_view(|view| {
        let bound = view.bind_slots(slots)?;
        Ok(view.df(bound.discounts()[0].points()[0]).value())
    })
}

#[test]
fn max_forward_integrates_full_grid_intervals_and_rebinds() -> Result<()> {
    // 半年時点は B の方が高いが、最初の格子区間の平均は A の方が高い。
    // 後半は B を採り、2 年目までの合計はどちらの候補の曲線とも異なる。
    let a: CurveId = "A".parse()?;
    let b: CurveId = "B".parse()?;
    let maximum: CurveId = "MAX".parse()?;
    let half = "2026-07-02T12:00:00Z";
    let pillars = [half.parse()?, FIRST.parse()?, SECOND.parse()?];
    let candidate = |id| {
        CurveSpec::new(
            id,
            Currency::JPY,
            [CurveRole::Discount],
            pillars,
            Extrapolation::Reject,
        )
    };
    let composed = CurveSpec::from_expr(
        maximum.clone(),
        Currency::JPY,
        [CurveRole::Discount],
        CurveExpr::MaxForward {
            parts: vec![CurveExpr::Ref(a.clone()), CurveExpr::Ref(b.clone())],
            grid: pillars[1..].to_vec(),
        },
    )?;
    assert_eq!(composed.parameter_blocks().count(), 0);
    let structure = Arc::new(MarketStructure::new(
        ORIGIN.parse()?,
        [candidate(a.clone())?, candidate(b.clone())?, composed],
    )?);
    let market = Market::from_theta(
        Arc::clone(&structure),
        vec![-0.005, -0.06, -0.08, -0.02, -0.03, -0.10],
    )?;
    for (date, expected_log) in [
        (ORIGIN, 0.0_f64),
        (half, -0.03),
        (FIRST, -0.06),
        ("2027-07-02T12:00:00Z", -0.095),
        (SECOND, -0.13),
    ] {
        let slots = discount_slots(&structure, &maximum, date)?;
        assert_close!(df(&market, &slots)?, expected_log.exp(), REL, ABS);
    }
    assert!(matches!(
        structure.prepare_point(&maximum, "2028-01-01T00:00:00.000001Z".parse()?),
        Err(Error::InvalidInput(message)) if message.contains("max forward grid")
    ));

    // 候補のブロックはいずれも必要で、準備を再利用したまま θ で勝者が変わる。
    let point = structure.prepare_point(&maximum, SECOND.parse()?)?;
    let a = BlockId::nodes(a);
    let b = BlockId::nodes(b);
    assert_eq!(
        point.blocks().cloned().collect::<BTreeSet<_>>(),
        BTreeSet::from([a.clone(), b]),
    );
    let slots = Slots::new(Arc::clone(&structure), [(maximum, vec![point])], [])?;
    let rebased = Market::from_theta(
        Arc::clone(&structure),
        vec![-0.02, -0.02, -0.13, -0.025, -0.07, -0.08],
    )?;
    assert_close!(df(&rebased, &slots)?, (-0.18_f64).exp(), REL, ABS);
    let partial = Market::from_partial_theta(
        structure,
        vec![-0.02, -0.02, -0.13, f64::NAN, f64::NAN, f64::NAN],
        [a],
