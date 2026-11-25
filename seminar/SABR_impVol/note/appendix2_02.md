    method: ObservationMethod,
    aggregation: Aggregation,
    #[serde(default)]
    rounding: ObservationRounding,
}

impl TryFrom<ObservationRuleInput> for ObservationRule {
    type Error = Error;

    fn try_from(input: ObservationRuleInput) -> Result<Self> {
        Self::new(input.days, input.method, input.aggregation, input.rounding)
    }
}

impl ObservationRule {
    /// 観測の規則を組み立てる。
    ///
    /// # Errors
    ///
    /// 日の選び方、方法、集約または丸めの組が許されない場合は [`Error::InvalidData`]。
    pub fn new(
        days: ObservationDays,
        method: ObservationMethod,
        aggregation: Aggregation,
        rounding: ObservationRounding,
    ) -> Result<Self> {
        // 設計書の「組の正しさ」の表を写す。PowerFactor とスプレッドを複利の中に入れる組は参照の式を
        // 持たず、指数の比は公表値の比なので観測値ごとの丸めを許さない。
        let valid = match (&days, method, aggregation) {
            (ObservationDays::BusinessDays(_), _, Aggregation::Compound { form, spread }) => {
                form == CompoundingForm::SimpleFactor || spread == SpreadTreatment::Exclusive
            }
            (ObservationDays::BusinessDays(_), _, Aggregation::Average(_))
            | (
                ObservationDays::AllDays { .. },
                ObservationMethod::Plain | ObservationMethod::Lockout(_),
                Aggregation::Average(_),
            )
            | (ObservationDays::Months, ObservationMethod::Plain, Aggregation::Average(_)) => true,
            (ObservationDays::Endpoints(_), ObservationMethod::Plain, Aggregation::IndexRatio) => {
                rounding.daily.is_none()
            }
            _ => false,
        };
        if !valid {
            return Err(invalid(
                "incompatible observation days, method, aggregation or rounding",
            ));
        }
        Ok(Self {
            days,
            method,
            aggregation,
            rounding,
        })
    }

    /// 観測点を選ぶ規則。
    #[must_use]
    pub const fn days(&self) -> &ObservationDays {
        &self.days
    }

    /// 観測日を変換する方法。
    #[must_use]
    pub const fn method(&self) -> ObservationMethod {
        self.method
    }

    /// 観測値の集約方法。
    #[must_use]
    pub const fn aggregation(&self) -> Aggregation {
        self.aggregation
    }

    /// 観測値と年率の丸め。
    #[must_use]
    pub const fn rounding(&self) -> ObservationRounding {
        self.rounding
    }
}

/// 一つの計算期間で使う観測点・重みと年率換算の係数。
#[derive(Debug, Clone, PartialEq)]
pub struct ObservationSchedule {
    fixings: Vec<(ObservationPoint, f64)>,
    coverage: Vec<DateInterval>,
    annualization: YearFraction,
    rule: ObservationRule,
}

impl ObservationSchedule {
    /// この観測列を作った規則。
    #[must_use]
    pub const fn rule(&self) -> &ObservationRule {
        &self.rule
    }

    /// 集約順の観測点と重み。等重みと指数の比の重みは 1。
    #[must_use]
    pub fn fixings(&self) -> &[(ObservationPoint, f64)] {
        &self.fixings
    }

    /// 各観測について、重みが覆う日付の区間を集約順に返す。
    #[must_use]
    pub fn coverage(&self) -> &[DateInterval] {
        &self.coverage
    }

    /// 集約値を年率へ換算するときに使う日数係数。
    #[must_use]
    pub const fn annualization(&self) -> &YearFraction {
        &self.annualization
    }
}

struct Observation {
    point: ObservationPoint,
    cover: DateInterval,
}

impl Observation {
    fn daily(fixing: Date, start: Date, end: Date) -> Result<Self> {
        Ok(Self {
            point: ObservationPoint::Date(fixing),
            cover: DateInterval::new(start, end)?,
        })
    }
}

/// 観測点と重みを計算する。参照枠の指定がなければ計算期間の枠を使う。
/// シフト後の観測期間に別の参照枠が必要なときは、その枠を文脈に渡す。
/// BUS/252 のカレンダーも観測日と同じ参照データから解決する。
///
/// # Errors
///
/// 暦が解決できないか、日数基準の文脈が不足するか、観測規則を期間に適用できない場合。
pub fn observation_schedule(
    period: &Period,
    rule: &ObservationRule,
    day_count: &DayCount,
    context: &AccrualContext<'_>,
    calendars: &BTreeMap<CalendarId, BaseCalendar>,
) -> Result<(ObservationSchedule, Vec<Warning>)> {
    let context = AccrualContext {
        frame: context.frame.or_else(|| period.frame()),
        ..*context
    };
    observation_schedule_using(period.accrual(), rule, day_count, &context, |spec| {
        Calendar::compile(spec, calendars)
    })
}

/// コンパイル済みの暦を再利用しながら計算期間の観測列を作る。
///
/// # Errors
///
/// [`observation_schedule`] と同じ失敗。
pub fn observation_schedule_with_calendars(
    period: &Period,
    rule: &ObservationRule,
    day_count: &DayCount,
    context: &AccrualContext<'_>,
    calendars: &mut CalendarCache<'_>,
) -> Result<(ObservationSchedule, Vec<Warning>)> {
    let context = AccrualContext {
        frame: context.frame.or_else(|| period.frame()),
        ..*context
    };
    observation_schedule_using(period.accrual(), rule, day_count, &context, |spec| {
        calendars.resolve(spec)
    })
}

/// 預金期間に対するフォールバックの観測日と重みを求める。
///
/// # Errors
///
/// 暦、日数基準、または期間から観測列を作れない場合。
pub fn fallback_observation_schedule(
    period: DateInterval,
    replacement: &crate::refdata::FallbackReplacement,
    calendars: &mut CalendarCache<'_>,
) -> Result<(ObservationSchedule, Vec<Warning>)> {
    observation_schedule_using(
        period,
        &replacement.observation,
        &replacement.annualization,
        &AccrualContext {
            maturity: Some(period.end()),
            ..AccrualContext::default()
        },
        |spec| calendars.resolve(spec),
    )
}

fn observation_schedule_using(
    accrual: DateInterval,
    rule: &ObservationRule,
    day_count: &DayCount,
    context: &AccrualContext<'_>,
    mut resolve: impl FnMut(&CalendarSpec) -> Result<Calendar>,
) -> Result<(ObservationSchedule, Vec<Warning>)> {
    let mut warnings = Vec::new();
    let (observations, annualization, fixing_calendar) = match &rule.days {
        ObservationDays::BusinessDays(spec) => {
            let calendar = resolve(spec)?;
            let annualization = if let ObservationMethod::ObservationShift(lag) = rule.method {
                shifted_interval(accrual, lag, &calendar, &mut warnings)?
            } else {
                accrual
            };
            let days = business_observations(annualization, rule.method, &calendar, &mut warnings)?;
            (days, annualization, Some(calendar))
        }
        ObservationDays::AllDays { fixing_calendar } => {
            let calendar = resolve(fixing_calendar)?;
            let days = all_days_observations(accrual, rule.method, &calendar, &mut warnings)?;
            (days, accrual, Some(calendar))
        }
        ObservationDays::Months => (monthly_observations(accrual)?, accrual, None),
        ObservationDays::Endpoints(endpoints) => {
            let (annualization, calendar) = match endpoints {
                EndpointDates::Accrual => (accrual, None),
                EndpointDates::Shifted { calendar, lag } => {
                    let calendar = resolve(calendar)?;
                    let shifted = shifted_interval(accrual, *lag, &calendar, &mut warnings)?;
                    (shifted, Some(calendar))
                }
            };
            let points = [annualization.start(), annualization.end()]
                .into_iter()
                .map(|date| Observation::daily(date, date, date))
                .collect::<Result<_>>()?;
            (points, annualization, calendar)
        }
    };
    // BUS/252 の営業日数は観測日と同じ参照データから数える。観測カレンダーが同じ基本暦なら、
    // compile 済みの値を使い回す。
    let bus252 = match day_count {
        DayCount::Bus252 { calendar: id } => {
            let matching = fixing_calendar.filter(|cal| cal.base_id() == Some(id));
            Some(matching.map_or_else(|| resolve(&CalendarSpec::Base(id.clone())), Ok)?)
        }
        _ => None,
    };
    let context = AccrualContext {
        calendar: bus252.as_ref(),
        ..*context
    };
    let mut fraction = |interval: DateInterval| {
        let (value, notices) =
            year_fraction(day_count, interval.start(), interval.end(), &context)?;
        warnings.extend(notices);
        Ok(value)
    };
    let annualization = fraction(annualization)?;
    let (fixings, coverage) = observations
        .into_iter()
        .map(|obs| {
            let weight = match rule.aggregation {
                Aggregation::Average(AveragingWeights::Equal) | Aggregation::IndexRatio => 1.0,
                Aggregation::Average(AveragingWeights::CoverWeighted)
                | Aggregation::Compound { .. } => fraction(obs.cover)?.value(),
            };
            Ok(((obs.point, weight), obs.cover))
        })
        .collect::<Result<_>>()?;
    Ok((
        ObservationSchedule {
            fixings,
            coverage,
            annualization,
            rule: rule.clone(),
        },
        unique_warnings(warnings),
    ))
}

fn shifted_interval(
    accrual: DateInterval,
    lag: NonZeroU8,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<DateInterval> {
    let lag = -i32::from(lag.get());
    let start = shift_business(accrual.start(), lag, calendar, warnings)?;
    let end = shift_business(accrual.end(), lag, calendar, warnings)?;
    if start >= end {
        return Err(invalid("shifted observation period is empty"));
    }
    DateInterval::new(start, end)
}

fn business_observations(
    interval: DateInterval,
    method: ObservationMethod,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<Vec<Observation>> {
    let start = interval.start();
    let end = interval.end();
    let mut days = Vec::new();
    for date in dates(interval) {
        let date = date?;
        if calendar.lookup(date, warnings)? {
            days.push(date);
        }
    }
    // 期間の開始日 S が営業日でなければ、S より前で最後の営業日 d_0 を先頭に補い、[S, d_1) を覆わせる
    // （非営業日の率は直前の営業日の値）。観測日が一つも無い期間も d_0 だけで全体を覆う。
    let carry = if days.first() == Some(&start) {
        None
    } else {
        Some(shift_business(start, -1, calendar, warnings)?)
    };
    // 観測 i の重みの区間は [d_i, d_{i+1}) で、最後の観測は期間の終了日 E までを覆う。
    let ends = days
        .iter()
        .copied()
        .skip(usize::from(carry.is_none()))
        .chain([end]);
    let mut previous_fixing = None;
    let mut observations = carry
        .into_iter()
        .chain(days.iter().copied())
        .zip(ends)
        .map(|(day, next)| {
            let fixing = match method {
                ObservationMethod::Lookback(lag) => {
                    // 観測日は連続する営業日なので、ラグを一度だけ戻してから順に進める。
                    let fixing = if let Some(previous) = previous_fixing {
                        shift_business(previous, 1, calendar, warnings)?
                    } else {
                        shift_business(day, -i32::from(lag.get()), calendar, warnings)?
                    };
                    previous_fixing = Some(fixing);
                    fixing
                }
                _ => day,
            };
            Observation::daily(fixing, day.max(start), next)
        })
        .collect::<Result<Vec<_>>>()?;
    if let ObservationMethod::Lockout(lag) = method {
        // 計算期間より前にある補完用の観測日は締め日の候補に数えない。
        lockout(&mut observations[usize::from(carry.is_some())..], lag)?;
    }
    Ok(observations)
}

// 各暦日を 1 日分の重みで観測し、非営業日には直前の営業日の fixing を置く。
fn all_days_observations(
    interval: DateInterval,
    method: ObservationMethod,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<Vec<Observation>> {
    let mut last_business = None;
    let mut observations = dates(interval)
        .map(|date| {
            let date = date?;
            let fixing = match (calendar.lookup(date, warnings)?, last_business) {
                (true, _) => date,
                (false, Some(previous)) => previous,
                (false, None) => shift_business(date, -1, calendar, warnings)?,
            };
            last_business = Some(fixing);
            Observation::daily(fixing, date, date.checked_add_days(1)?)
        })
        .collect::<Result<Vec<_>>>()?;
    if let ObservationMethod::Lockout(lag) = method {
        lockout(&mut observations, lag)?;
    }
    Ok(observations)
}

// Lockout(n) は観測 d_1..d_m のうち i > m − n の fixing を締めの日 d_{m−n} の値にする。
// 0 始まりの添字では締めの日が m − n − 1 にあたり、m ≤ n なら締めの日が無い。
fn lockout(observations: &mut [Observation], lag: NonZeroU8) -> Result<()> {
    let cutoff = observations
        .len()
        .checked_sub(usize::from(lag.get()) + 1)
        .ok_or_else(|| invalid("lockout requires a cutoff observation"))?;
    let fixing = observations[cutoff].point;
    for observation in &mut observations[cutoff + 1..] {
        observation.point = fixing;
    }
    Ok(())
}

fn monthly_observations(interval: DateInterval) -> Result<Vec<Observation>> {
    let mut start = interval.start();
    let end = interval.end();
    let mut observations = Vec::new();
    while start < end {
        let (year, month, _) = start.ymd();
        let month_end = start.month_end()?;
        let next = if end <= month_end {
            end
        } else {
            month_end.checked_add_days(1)?
        };
        observations.push(Observation {
            point: ObservationPoint::Month(YearMonth::new(year, month)?),
            cover: DateInterval::new(start, next)?,
        });
        start = next;
    }
    Ok(observations)
}

fn dates(interval: DateInterval) -> impl Iterator<Item = Result<Date>> {
    let start = interval.start();
    (0..interval.end().days_since(start)).map(move |days| start.checked_add_days(days))
}

fn invalid(message: &str) -> Error {
    Error::InvalidData(message.into())
}
`````

`````
crates/rhodium-core/src/parallel.rs
//! 取引などの対象を一段だけ並列に分けて処理する。
//!
//! 並列の実行器は scioxide の [`ExecutionContext`] で、呼び出し側が引数で渡す。
//! [`ExecutionContext::sequential`] を渡せば逐次に処理する。
//! 対象は連続した区間に分け、区間ごとの結果を区間の順に返す。
//! 区間の中の処理は子の実行文脈を受け取り、入れ子の並列はしない。
//!
//! 実行器が取消と期限を確かめるのは区間の開始前だけなので、区間の中の処理は対象の境界で
//! [`check`] を呼ぶ。止まった実行は対象ごとの失敗にせず、呼び出し全体の `Execution` にする。

use std::{num::NonZeroUsize, ops::Range};

pub use scioxide_core::{ExecutionContext, ThreadExecutor};
use scioxide_core::{
    ExecutionFailureReason, ExecutionPolicy, PlanLimits, SciOxideError, SciOxideErrorView,
    backend::PartitionSlot,
};

use crate::{Error, Result};

/// 対象 `0..len` を連続した区間に分けて `run` を並列に呼び、区間の順に結果を返す。
///
/// 区間の数は実行器のワーカー数と対象の数の小さい方で、どの区間も 1 件以上の対象を持つ。
/// パスを持たない評価（リニア評価や準備）は常に対象で分ける。
/// `run` は区間ごとに一度だけ呼ばれ、区間の範囲と子の実行文脈を受け取る。
/// tape や作業領域は `run` の中で一つ作り、区間の対象に使い回す。
/// 区間ごとの結果を順に連結すれば対象の順になり、区間ごとの和を順に足せば
/// 実行器の完了の順によらない縮約になる（分け方が変われば和の順は変わる）。
///
/// 対象が 0 件でも、取消か期限で止まっていれば `Execution` を返す。
///
/// # Errors
///
/// 取り消しや期限で止まった（対象が 0 件の場合を含む）、ワーカーや領域を確保できなかった、
/// または実行器自身が失敗したら `Execution`。
/// `run` が panic したか、実行器が分割の約束を破ったら `Internal`。
pub fn map_partitions<T: Send>(
    execution: ExecutionContext<'_>,
    len: usize,
    run: impl for<'s> Fn(Range<usize>, ExecutionContext<'s>) -> T + Sync,
) -> Result<Vec<T>> {
    // 対象が無いと実行器は区間を作らず、区間の開始前の停止の確認も行わない。取り消した実行が
    // 空の結果で成功しないよう、ここで一度確かめる（対象があるときの確認の回数は変えない）。
    if len == 0 {
        check(execution)?;
    }
    // ワーカー数の上限は実行器の予算と対象の数が決める。区間の大きさは等しく分ける。
    let plan = execution.plan(
        ExecutionPolicy::Parallel(execution.worker_budget()),
        len,
        PlanLimits {
            auto_threshold: NonZeroUsize::MIN,
            work_per_item: NonZeroUsize::MIN,
            worker_limit: NonZeroUsize::MAX,
        },
    );
    let mut results: Vec<Option<T>> = std::iter::repeat_with(|| None)
        .take(plan.worker_count())
        .collect();
    {
        // 区間ごとに出力の場所を一つ渡し、区間の番号の順に結果を並べる。
        let slots: Vec<_> = results.iter_mut().map(PartitionSlot::new).collect();
        execution
            .execute_partitioned(plan, &slots, &|partition, slot, child| {
                *slot = Some(run(partition.range(), child));
                Ok(())
            })
            .map_err(execution_error)?;
    }
    results
        .into_iter()
        .map(|result| {
            result.ok_or_else(|| Error::Internal("executor left a partition without result".into()))
        })
        .collect()
}

/// 取消か期限で止まっていれば `Execution` を返す。
///
/// 区間の中の処理が対象（取引など）の境界で呼ぶ。実行器は区間の開始前にしか確かめないので、
/// 区間の開始後の取消や期限は、この確認か、内側の計算（AD の記録や逆伝播など）が見つける。
///
/// # Errors
///
/// 取消か期限で止まっているか、期限の時計が失敗したら `Execution`。
pub fn check(execution: ExecutionContext<'_>) -> Result<()> {
    execution.check().map_err(execution_error)
}

/// scioxide のエラーが実行の制御（取消、期限、時計、領域の確保、ワーカーの作成）によるものなら、
/// 呼び出し全体を止める `Execution` を返す。
///
/// 内側の計算が返した scioxide のエラーのうち、これに当たらないもの（数値や入力の失敗）は
/// `None` を返し、呼び出し側がその計算の文脈のエラーにする。
/// 実行の制御の失敗を対象ごとの失敗に変えると、止めた計算が一部の対象を除いた結果として
/// 成功してしまうため、この分類で分ける。
#[must_use]
pub fn control_failure(error: &SciOxideError) -> Option<Error> {
    match error.view() {
        // 呼び出し側が渡した取り消し・期限・時計と、実行環境の資源の不足。
        SciOxideErrorView::ExecutionStopped { .. }
        | SciOxideErrorView::CapabilityUnavailable { .. }
        | SciOxideErrorView::CapabilityFailure { .. }
        | SciOxideErrorView::AllocationFailure { .. }
        | SciOxideErrorView::ExecutionFailure {
            reason:
                ExecutionFailureReason::WorkerCreation { .. } | ExecutionFailureReason::ProviderFailure,
            ..
        } => Some(Error::Execution(format!("parallel execution: {error}"))),
        _ => None,
    }
}

/// scioxide の実行のエラーを、直す人で分けた rhodium のエラーにする。
#[expect(clippy::needless_pass_by_value, reason = "map_err の変換関数に使う")]
fn execution_error(error: SciOxideError) -> Error {
    // 実行の制御でなければ、コールバックの panic か分割の約束を破る実行器で、rhodium の誤り。
    control_failure(&error)
        .unwrap_or_else(|| Error::Internal(format!("parallel execution: {error}")))
}
`````

`````
crates/rhodium-core/src/quantity.rs
//! 金融計算で扱う、意味を持つ量。

mod prices;
mod rates;
mod volatility;

use scioxide_ad::{Mask, Real};

pub use prices::{
    Clean, Correlation, Dirty, HazardRate, ParPrice, Probability, RateFuturesPrice, RecoveryRate,
    RiskyAnnuity, RunningSpread, UpfrontFraction,
};
pub use rates::{Compounded, Continuous, Rate, RateKind, Simple, Spread, accrue};
pub use volatility::{CallPut, LognormalVol, NormalVol, ShiftedLognormalVol, VolKind};

use crate::{Error, Result, day_count::DayCount};

/// 日数基準のタグを持つ日数係数。
///
/// ファイルでは `{ value = 0.0164, basis = "ACT365F" }` と書き、読み込みで有限性を確かめる。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(try_from = "YearFractionInput")]
pub struct YearFraction {
    value: f64,
    basis: DayCount,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct YearFractionInput {
    value: f64,
    basis: DayCount,
}

impl TryFrom<YearFractionInput> for YearFraction {
    type Error = Error;

    fn try_from(input: YearFractionInput) -> Result<Self> {
        Self::new(input.value, input.basis)
    }
}

impl YearFraction {
    /// 有限な日数係数に日数基準を指定する。
    ///
    /// # Errors
    ///
    /// 係数が有限でなければ [`Error::InvalidInput`]。
    pub fn new(value: f64, basis: DayCount) -> Result<Self> {
        if !value.is_finite() {
            return Err(Error::InvalidInput("year fraction is not finite".into()));
        }
        Ok(Self { value, basis })
    }

    /// 係数の値。
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.value
    }

    /// 係数を求めるのに使った日数基準。
    #[must_use]
    pub const fn basis(&self) -> &DayCount {
        &self.basis
    }
}

/// 評価時点からの正の割引因子。生存確率と区別する。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiscountFactor<R = f64>(R);

impl<R: Real> DiscountFactor<R> {
    /// 正の有限な割引因子を作る。負金利のときは 1 を超えられる。
    ///
    /// # Errors
    ///
    /// 正の有限値でなければ `InvalidInput`。
    pub fn new(value: R) -> Result<Self> {
        if !value.is_finite().all() || !value.gt(R::zero()).all() {
            return Err(Error::InvalidInput(
                "discount factor must be finite and positive".into(),
            ));
        }
        Ok(Self(value))
    }

    /// 割引因子の値。
    #[must_use]
    pub const fn value(self) -> R {
        self.0
    }

    /// 区間の成長因子 P(s)/P(e) を求める。
    ///
    /// # Errors
    ///
    /// 比が正の有限値にならなければ `Numerical`。
    pub fn growth_to(self, end: Self) -> Result<R> {
        let factor = self.0 / end.0;
        if !factor.is_finite().all() || !factor.gt(R::zero()).all() {
            return Err(Error::Numerical(
                "interval discount-factor ratio is invalid".into(),
            ));
        }
        Ok(factor)
    }
}

/// 評価時点からの生存確率。割引因子と区別する。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurvivalProbability<R = f64>(R);

impl<R: Real> SurvivalProbability<R> {
    /// 0 以上 1 以下の有限な生存確率を作る。
    ///
    /// # Errors
    ///
    /// 範囲外か非有限の値なら `InvalidInput`。
    pub fn new(value: R) -> Result<Self> {
        if !value.is_finite().all() || !value.ge(R::zero()).all() || !R::one().ge(value).all() {
            return Err(Error::InvalidInput(
                "survival probability must be within [0, 1]".into(),
            ));
        }
        Ok(Self(value))
    }

    /// 生存確率の値。
    #[must_use]
    pub const fn value(self) -> R {
        self.0
    }

    /// s で生存していた条件下で、e まで生存する確率 Q(e)/Q(s)。
    ///
    /// # Errors
    ///
    /// 開始の確率がゼロか、終了の確率が開始より大きければ `InvalidInput`。
    /// 計算が非有限なら `Numerical`。
    pub fn survival_to(self, end: Self) -> Result<Self> {
        if !self.0.gt(R::zero()).all() || !self.0.ge(end.0).all() {
            return Err(Error::InvalidInput(
                "invalid interval survival probabilities".into(),
            ));
        }
        Self::new(end.0 / self.0)
            .map_err(|_| Error::Numerical("interval survival probability is invalid".into()))
    }
}

fn finite(value: f64, name: &str) -> Result<f64> {
    if !value.is_finite() {
        return Err(Error::InvalidInput(format!("{name} is not finite")));
    }
    Ok(value)
}

fn numerical(value: f64, name: &str) -> Result<f64> {
    if !value.is_finite() {
        return Err(Error::Numerical(format!("calculated {name} is not finite")));
    }
    Ok(value)
}
`````

`````
crates/rhodium-core/src/quantity/prices.rs
use std::marker::PhantomData;

use super::{Rate, Simple, finite, numerical};
use crate::{Error, Result, currency::Currency, day_count::DayCount, id::CreditKey, money::Money};

/// 経過利息を含まない価格の区分。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Clean;

/// 経過利息を含む価格の区分。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Dirty;

/// 額面 1 あたりの価格。経過利息の有無を型で区別する。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParPrice<P> {
    value: f64,
    basis: PhantomData<P>,
}

impl<P> ParPrice<P> {
    /// 有限な額面価格を作る。
    ///
    /// # Errors
    ///
    /// 非有限値は `InvalidInput`。
    pub fn new(value: f64) -> Result<Self> {
        Ok(Self {
            value: finite(value, "par price")?,
            basis: PhantomData,
        })
    }

    /// 額面 1 あたりの価格。
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.value
    }
}

impl ParPrice<Clean> {
    /// 額面 1 あたりの経過利息を加え、ダーティ価格を求める。
    ///
    /// # Errors
    ///
    /// 経過利息が非有限なら `InvalidInput`、ダーティ価格が非有限なら `Numerical`。
    pub fn with_accrued(self, accrued: f64) -> Result<ParPrice<Dirty>> {
        let accrued = finite(accrued, "accrued fraction")?;
        Ok(ParPrice {
            value: numerical(self.value + accrued, "dirty price")?,
            basis: PhantomData,
        })
    }
}

impl ParPrice<Dirty> {
    /// 額面金額に対する決済額。
    ///
    /// # Errors
    ///
    /// 結果が非有限なら `Numerical`。
    pub fn settlement(self, face: Money) -> Result<Money> {
        numerical(self.value * face.amount(), "settlement amount")
            .map(|amount| Money::new(face.ccy(), amount))
    }
}

/// 短期金利先物の価格。値は 1 − 率。
#[derive(Debug, Clone, PartialEq)]
pub struct RateFuturesPrice {
    value: f64,
    basis: DayCount,
}

impl RateFuturesPrice {
    /// 価格と原資産の率の日数基準を指定する。
    ///
    /// # Errors
    ///
    /// 非有限値は `InvalidInput`。
    pub fn new(value: f64, basis: DayCount) -> Result<Self> {
        Ok(Self {
            value: finite(value, "rate futures price")?,
            basis,
        })
    }

    /// 1 − 率の値。
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.value
    }

    /// 原資産の率の日数基準。
    #[must_use]
    pub const fn basis(&self) -> &DayCount {
        &self.basis
    }

    /// 先物価格から単利の率を求める。
    ///
    /// # Errors
    ///
    /// 換算結果が非有限なら `Numerical`。
    pub fn implied_rate(self) -> Result<Rate<Simple>> {
        Rate::new(1.0 - self.value, self.basis, Simple)
            .map_err(|_| Error::Numerical("calculated futures rate is non-finite".into()))
    }
}

macro_rules! scalar_quantity {
    ($(#[$doc:meta] $name:ident, $label:literal, $domain:expr;)+) => {
        $(
            #[$doc]
            #[derive(Debug, Clone, Copy, PartialEq)]
            pub struct $name(f64);

            impl $name {
                /// 定義域に収まる有限な値から作る。
                ///
                /// # Errors
                ///
                /// 非有限値または定義域外は `InvalidInput`。
                pub fn new(value: f64) -> Result<Self> {
                    let value = finite(value, $label)?;
                    if !($domain)(value) {
                        return Err(Error::InvalidInput(concat!($label, " is out of range").into()));
                    }
                    Ok(Self(value))
                }

                /// 小数表示の値。
                #[must_use]
                pub const fn value(self) -> f64 {
                    self.0
                }
            }
        )+
    };
}

scalar_quantity! {
    /// CDS のランニング・スプレッド。
    RunningSpread, "running spread", |_| true;
    /// CDS の初期に支払う、元本に対する比率。
    UpfrontFraction, "upfront fraction", |_| true;
    /// ISDA の標準モデルで計算した元本 1 あたりの risky PV01。
    RiskyAnnuity, "risky annuity", |value| value >= 0.0;
    /// 信用事由後の回収率。
    RecoveryRate, "recovery rate", |value| (0.0..=1.0).contains(&value);
    /// 事象が起きる確率。
    Probability, "probability", |value| (0.0..=1.0).contains(&value);
    /// 二つの量の相関係数。
    Correlation, "correlation", |value| (-1.0..=1.0).contains(&value);
}

impl RunningSpread {
    /// ISDA の標準モデルで評価済みの risky PV01 を掛け、元本あたりの額にする。
    ///
    /// # Errors
    ///
    /// 結果が非有限なら `Numerical`。
    pub fn to_upfront(self, risky_pv01: RiskyAnnuity) -> Result<UpfrontFraction> {
        numerical(self.0 * risky_pv01.value(), "CDS upfront fraction").map(UpfrontFraction)
    }
}

/// 信用主体と契約条件に対応する、非負のハザード率。
#[derive(Debug, Clone, PartialEq)]
pub struct HazardRate {
    value: f64,
    credit: CreditKey,
}

impl HazardRate {
    /// 信用参照とハザード率を結び付ける。
    ///
    /// # Errors
    ///
    /// 非有限または負の率は `InvalidInput`。
    pub fn new(value: f64, credit: CreditKey) -> Result<Self> {
        let value = finite(value, "hazard rate")?;
        if value < 0.0 {
            return Err(Error::InvalidInput("hazard rate is negative".into()));
        }
        Ok(Self { value, credit })
    }

    /// 小数表示の率。
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.value
    }

    /// 信用主体と契約条件。
    #[must_use]
    pub const fn credit(&self) -> &CreditKey {
        &self.credit
    }

    /// 契約の通貨。
    #[must_use]
    pub const fn ccy(&self) -> Currency {
        self.credit.ccy
    }
}
`````

`````
crates/rhodium-core/src/quantity/rates.rs
use std::{
    num::NonZeroU8,
    ops::{Add, Sub},
};

use super::{DiscountFactor, YearFraction, finite, numerical};
use crate::{Error, Result, day_count::DayCount, tenor::Frequency};

/// 単利。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Simple;

/// 期間内の複利。頻度は一年を整数個の期間に分割できる。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Compounded {
    frequency: Frequency,
    periods_per_year: NonZeroU8,
}

impl Compounded {
    /// 一年あたりの回数が決まる頻度から複利型を作る。
    ///
    /// # Errors
    ///
    /// 一年あたりの回数を決められなければ `InvalidInput`。
    pub fn new(frequency: Frequency) -> Result<Self> {
        let periods_per_year = frequency
            .periods_per_year()
            .ok_or_else(|| Error::InvalidInput("invalid compounding frequency".into()))?;
        Ok(Self {
            frequency,
            periods_per_year,
        })
    }

    fn log_growth(self, rate: f64, years: f64) -> f64 {
        let periods = f64::from(self.periods_per_year.get());
        let per_period = rate / periods;
        if per_period <= -1.0 {
            return f64::NAN;
        }
        periods * per_period.ln_1p() * years
    }

    /// 複利の頻度。
    #[must_use]
    pub const fn frequency(self) -> Frequency {
        self.frequency
    }
}

/// 連続複利。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Continuous;

mod sealed {
    pub trait Sealed {
        fn growth_excess(self, rate: f64, years: f64) -> f64;
        fn growth(self, rate: f64, years: f64) -> f64;
        fn discount(self, rate: f64, years: f64) -> f64;
    }
}

/// 率の複利型。単利、頻度付きの複利、連続複利の三つ。
pub trait RateKind: sealed::Sealed + Copy + Eq {}

impl RateKind for Simple {}
impl RateKind for Compounded {}
impl RateKind for Continuous {}

impl sealed::Sealed for Simple {
    fn growth_excess(self, rate: f64, years: f64) -> f64 {
        rate * years
    }

    fn growth(self, rate: f64, years: f64) -> f64 {
        rate.mul_add(years, 1.0)
    }

    fn discount(self, rate: f64, years: f64) -> f64 {
        self.growth(rate, years).recip()
    }
}

impl sealed::Sealed for Compounded {
    fn growth_excess(self, rate: f64, years: f64) -> f64 {
        self.log_growth(rate, years).exp_m1()
    }

    fn growth(self, rate: f64, years: f64) -> f64 {
        self.log_growth(rate, years).exp()
    }

    fn discount(self, rate: f64, years: f64) -> f64 {
        (-self.log_growth(rate, years)).exp()
    }
}

impl sealed::Sealed for Continuous {
    fn growth_excess(self, rate: f64, years: f64) -> f64 {
        (rate * years).exp_m1()
    }

    fn growth(self, rate: f64, years: f64) -> f64 {
        (rate * years).exp()
    }

    fn discount(self, rate: f64, years: f64) -> f64 {
        (-rate * years).exp()
    }
}

/// 複利の種類と日数基準を持つ年率。
#[derive(Debug, Clone, PartialEq)]
pub struct Rate<K> {
    value: f64,
    basis: DayCount,
    kind: K,
}

impl<K: RateKind> Rate<K> {
    /// 有限な率に日数基準と複利型を指定する。
    ///
    /// # Errors
    ///
    /// 非有限値は `InvalidInput`。
    pub fn new(value: f64, basis: DayCount, kind: K) -> Result<Self> {
        Ok(Self {
            value: finite(value, "rate")?,
            basis,
            kind,
        })
    }

    /// 小数表示の年率。
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.value
    }

    /// 率の日数基準。
    #[must_use]
    pub const fn basis(&self) -> &DayCount {
        &self.basis
    }

    /// 率の複利型。
    #[must_use]
    pub const fn kind(&self) -> K {
        self.kind
    }

    /// 同じ日数基準、複利型、頻度の率の和。
    ///
    /// # Errors
    ///
    /// 基準か頻度が異なれば `InvalidInput`、結果が非有限なら `Numerical`。
    pub fn checked_add(self, other: &Self) -> Result<Self> {
        let result = self.value + other.value;
        self.combine(other, result)
    }

    /// 同じ日数基準、複利型、頻度の率の差。
    ///
    /// # Errors
    ///
    /// 基準か頻度が異なれば `InvalidInput`、結果が非有限なら `Numerical`。
    pub fn checked_sub(self, other: &Self) -> Result<Self> {
        let result = self.value - other.value;
        self.combine(other, result)
    }

    fn combine(mut self, other: &Self, result: f64) -> Result<Self> {
        if self.basis != other.basis || self.kind != other.kind {
            return Err(Error::InvalidInput(
                "rate basis or compounding frequency differs".into(),
            ));
        }
        self.value = numerical(result, "rate")?;
        Ok(self)
    }

    /// 指定した日数係数に対する成長因子。
    ///
    /// # Errors
    ///
    /// 日数基準の不一致は `InvalidInput`、成長因子が正の有限値にならなければ `Numerical`。
    pub fn growth_factor(&self, fraction: &YearFraction) -> Result<f64> {
        self.check_basis(fraction)?;
        let growth = self.kind.growth(self.value, fraction.value());
        if !growth.is_finite() || growth <= 0.0 {
            return Err(Error::Numerical(
                "growth factor is not finite and positive".into(),
            ));
        }
        Ok(growth)
    }

    /// 指定した日数係数に対する割引因子。
    ///
    /// # Errors
    ///
    /// 日数基準の不一致は `InvalidInput`、割引因子が正の有限値にならなければ `Numerical`。
    pub fn discount_factor(&self, fraction: &YearFraction) -> Result<DiscountFactor> {
        self.check_basis(fraction)?;
        let discount = self.kind.discount(self.value, fraction.value());
        DiscountFactor::new(discount)
            .map_err(|_| Error::Numerical("calculated discount factor is invalid".into()))
    }

    fn check_basis(&self, fraction: &YearFraction) -> Result<()> {
        if self.basis() != fraction.basis() {
            return Err(Error::InvalidInput(format!(
                "rate basis {} differs from year fraction basis {}",
                self.basis(),
                fraction.basis()
            )));
        }
        Ok(())
    }
}

impl<K: RateKind> Add for Rate<K> {
    type Output = Result<Self>;

    fn add(self, rhs: Self) -> Self::Output {
        self.checked_add(&rhs)
    }
}

impl<K: RateKind> Sub for Rate<K> {
    type Output = Result<Self>;

    fn sub(self, rhs: Self) -> Self::Output {
        self.checked_sub(&rhs)
    }
}

/// 単利の率の差。値の単位は小数。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spread(f64);

impl std::str::FromStr for Spread {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        let (number, shift) = text
            .strip_suffix('%')
            .map(|number| (number, 2))
            .or_else(|| text.strip_suffix("bp").map(|number| (number, 4)))
            .ok_or_else(|| Error::InvalidData("rate must end with % or bp".into()))?;
        let unsigned = number.strip_prefix(['+', '-']).unwrap_or(number);
        let (whole, fraction) = unsigned.split_once('.').unwrap_or((unsigned, ""));
        if whole.is_empty() && fraction.is_empty()
            || !whole
                .bytes()
                .chain(fraction.bytes())
                .all(|b| b.is_ascii_digit())
        {
            return Err(Error::InvalidData("rate must be a scaled decimal".into()));
        }
        // 十進の指数を変えてから変換し、同じ値の bp と % を同じ f64 にする。
        let value = format!("{number}e-{shift}")
            .parse::<f64>()
            .map_err(|_| Error::InvalidData("invalid rate".into()))?;
        if !value.is_finite() {
            return Err(Error::InvalidData("rate must be finite".into()));
        }
        Ok(Self(value))
    }
}

impl Spread {
    /// 有限なスプレッドを作る。
    ///
    /// # Errors
    ///
    /// 非有限値は `InvalidInput`。
    pub fn new(value: f64) -> Result<Self> {
        finite(value, "spread").map(Self)
    }

    /// 小数表示のスプレッド。
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0
    }
}

impl Add for Spread {
    type Output = Result<Self>;

    fn add(self, rhs: Self) -> Self::Output {
        numerical(self.0 + rhs.0, "spread").map(Self)
    }
}

impl Add<Spread> for Rate<Simple> {
    type Output = Result<Self>;

    fn add(mut self, rhs: Spread) -> Self::Output {
        self.value = numerical(self.value + rhs.value(), "rate")?;
        Ok(self)
    }
}

/// 複利型に従い、元本 1 あたりの経過利息を求める。
///
/// # Errors
///
/// 率と係数の日数基準が異なれば `InvalidInput`、結果が非有限なら `Numerical`。
pub fn accrue<K: RateKind>(rate: &Rate<K>, fraction: &YearFraction) -> Result<f64> {
    rate.check_basis(fraction)?;
    numerical(
        rate.kind.growth_excess(rate.value, fraction.value()),
        "accrued interest",
    )
}
`````

`````
crates/rhodium-core/src/quantity/volatility.rs
use std::cmp::Ordering;

use super::{finite, numerical};
use crate::{Error, Result};

/// ボラの種類を選択するための分類。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VolKind {
    /// 正規モデルのボラ。
    Normal,
    /// 対数正規モデルのボラ。
    Lognormal,
    /// シフトを加えた対数正規モデルのボラ。
    ShiftedLognormal {
        /// 原資産に加えるシフト。
        shift: f64,
    },
}

/// オプションの権利の向き。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CallPut {
    /// 買う権利。
    Call,
    /// 売る権利。
    Put,
}

macro_rules! volatility {
    ($(#[$doc:meta] $name:ident, $label:literal;)+) => {
        $(
            #[$doc]
            #[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
            pub struct $name(f64);

            impl $name {
                /// 非負の有限なボラを作る。
                ///
                /// # Errors
                ///
                /// 負の値や非有限値なら `InvalidInput`。
                pub fn new(value: f64) -> Result<Self> {
                    nonnegative(value, $label).map(Self)
                }

                /// ボラの数値。
                #[must_use]
                pub const fn value(self) -> f64 {
                    self.0
                }

                /// 同じ種類のボラとの差。負の差も返す。
                #[must_use]
                pub fn difference(self, other: Self) -> f64 {
                    self.0 - other.0
                }
            }
        )+
    };
}

volatility! {
    /// 正規ボラ。原資産の単位で表す。
    NormalVol, "normal volatility";
    /// 対数正規ボラ。比率を小数で表す。
    LognormalVol, "lognormal volatility";
}

/// 原資産に加えるシフトを持つ対数正規ボラ。
#[derive(Debug, Clone, Copy)]
pub struct ShiftedLognormalVol {
    value: f64,
    shift: f64,
}

impl ShiftedLognormalVol {
    /// 非負の有限なボラと有限なシフトを結び付ける。
    ///
    /// # Errors
    ///
    /// ボラまたはシフトが定義域に収まらなければ `InvalidInput`。
    pub fn new(value: f64, shift: f64) -> Result<Self> {
        Ok(Self {
            value: nonnegative(value, "shifted lognormal volatility")?,
            shift: finite(shift, "volatility shift")?,
        })
    }

    /// ボラの数値。
    #[must_use]
    pub const fn value(self) -> f64 {
        self.value
    }

    /// 原資産に加えるシフト。
    #[must_use]
    pub const fn shift(self) -> f64 {
        self.shift
    }

    /// 同じシフトを持つボラを比べる。
    ///
    /// # Errors
    ///
    /// シフトが異なれば `InvalidInput`。
    pub fn compare(self, other: Self) -> Result<Ordering> {
        self.check_shift(other)?;
        // -0.0 と 0.0 は同じボラとして比較する。
        self.value
            .partial_cmp(&other.value)
            .ok_or_else(|| Error::Internal("validated volatility is not comparable".into()))
    }

    /// 同じシフトを持つボラとの差。負の差も返す。
    ///
    /// # Errors
    ///
    /// シフトが異なれば `InvalidInput`。
    pub fn difference(self, other: Self) -> Result<f64> {
        self.check_shift(other)?;
        numerical(self.value - other.value, "volatility difference")
    }

    fn check_shift(self, other: Self) -> Result<()> {
        if self.shift.partial_cmp(&other.shift) != Some(Ordering::Equal) {
            return Err(Error::InvalidInput("volatility shifts differ".into()));
        }
        Ok(())
    }
}

fn nonnegative(value: f64, name: &str) -> Result<f64> {
    let value = finite(value, name)?;
    if value < 0.0 {
        return Err(Error::InvalidInput(format!("{name} is negative")));
    }
    Ok(value)
}
`````

`````
crates/rhodium-core/src/refdata.rs
//! 通貨、地域、カレンダー、指標と慣行の参照データ。

mod fallback;
mod fx_convention;
mod index;
mod rate_future;

use std::{
    collections::{BTreeMap, btree_map::Entry},
    sync::Arc,
};

use crate::{
    Error, Result, Warning,
    calendar::{
        AdvanceRule, BaseCalendar, BusinessDayConvention, Calendar, CalendarSpec, EndOfMonthRule,
        money_market_period, shift_business, unique_warnings,
    },
    currency::{Asset, Currency},
    date::{Date, DateInterval},
    day_count::DayCount,
    fixing::{ObservationPeriod, ObservationPoint},
    fx::CcyPair,
    id::{CalendarId, CommodityId, ConventionId, EntityId, IndexId, SecurityId},
    observation::{EndpointDates, ObservationDays, ObservationRule},
    rounding::Rounding,
    schedule::PaymentRule,
    tenor::{DateOffset, Frequency},
    time::{Instant, LocalEventTime, TzId, ZoneRules},
};

pub use fallback::{
    FallbackApplicability, FallbackDefinition, FallbackReplacement, FallbackTrigger,
};
pub use fx_convention::{FxConvention, FxSpotLag};
pub use index::{
    DailyPublication, IndexDefinition, IndexKind, IndexStatus, InflationInterpolation, PriceKind,
    PublicationRule, Seasonality,
};
pub use rate_future::{
    RateFutureClearing, RateFutureContract, RateFutureConvention, RateFutureRate,
};

use index::{ResolvedPublication, credit_window_end, incompatible_point};

/// 現物決済の可否。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Deliverability {
    /// 現物を受け渡せる。
    Deliverable,
    /// FX fixing による差金で決済する。
    NonDeliverable,
}

/// 通貨の最小単位と支払の慣行。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CurrencyDefinition {
    /// この定義の通貨コード。
    pub code: Currency,
    /// 通貨の最小単位の十進桁数。
    pub minor_units: u8,
    /// 契約に従って支払額に適用する丸め。
    pub payment_rounding: Rounding,
    /// 現物決済の可否。
    pub deliverability: Deliverability,
    /// 決済の地域と現地時刻。
    pub settlement_time: LocalEventTime,
    /// 決済日の基本カレンダー。
    pub settlement_calendar: CalendarId,
}

/// 受け渡す証券の種類。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SecurityKind {
    /// 株式。
    Share,
    /// 債券。
    Bond,
    /// ファンドの持分。
    FundUnit,
    /// 株価指数の取引単位。
    IndexUnit,
    /// 上場先物。
    ListedFuture,
    /// 上場オプション。
    ListedOption,
}

/// 証券 1 単位の意味。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SecurityUnit {
    /// 1 株。
    Share,
    /// 額面の通貨単位 1。
    FaceValue,
    /// 1 持分。
    FundUnit,
    /// 先物やオプションの 1 枚。
    Contract,
    /// 指数の 1 ポイント。
    IndexPoint,
}

/// 保有資産の評価に用いる種類。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AssetValuation {
    /// 市場のフォワード。
    Forward,
    /// 証券自体の契約。
    Contract,
}

/// 証券の参照データ。
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SecurityDefinition {
    /// 証券の識別子。
    pub id: SecurityId,
    /// 証券の種類。
    pub kind: SecurityKind,
    /// 価格や額面の通貨。
    pub currency: Currency,
    /// 証券 1 単位の意味。
    pub unit: SecurityUnit,
    /// 発行体。債券では必須。
    pub issuer: Option<EntityId>,
    /// 証券の評価方法。
    pub valuation: AssetValuation,
}

/// 受渡場所が決まったコモディティの定義。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CommodityDefinition {
    /// コモディティの識別子。
    pub id: CommodityId,
    /// 受渡数量の物理単位。
    pub unit: String,
    /// 価格を表す通貨。
    pub currency: Currency,
}

/// 商品の条件に適用する型付きの慣行。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum Convention {
    /// 預金の慣行。
    Deposit(DepositConvention),
    /// OIS の慣行。
    Ois(Box<OisConvention>),
    /// 取引所上場の金利先物の慣行と限月。
    RateFuture(Box<RateFutureConvention>),
    /// 期首決済の FRA。
    Fra(FraConvention),
    /// 固定対ターム金利スワップ。
    Irs(Box<SwapConvention>),
    /// 同一通貨の変動対変動スワップ。
    BasisSwap(Box<SwapConvention>),
    /// クロスカレンシースワップ。
    XccySwap(Box<SwapConvention>),
    /// FX フォワードのスポット起点と受渡。
    FxForward(Box<FxForwardConvention>),
}

/// FRA の原指標とスポット起点。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FraConvention {
    /// ターム金利指標。
    pub index: IndexId,
    /// 約定日からスポット日までのずれ。
    pub spot_lag: DateOffset,
}

/// 確認書の元本・quote・満期を除く、標準スワップの条件。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SwapConvention {
    /// 約定日から開始日までのずれ。
    pub spot_lag: DateOffset,
    /// 払い脚と受け脚。ベーシスでは払い脚のスプレッドを quote する。
    pub legs: [SwapLegConvention; 2],
}

/// スワップの一脚の通貨、日程と観測の慣行。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SwapLegConvention {
    /// 支払通貨。
    pub currency: Currency,
    /// 頻度、日数基準および営業日規則。
    pub periods: FixedLegConvention,
    /// 支払い日の算出規則。
    pub payment: PaymentRule,
    /// クーポンの種類と指標。
    pub coupon: SwapCouponConvention,
    /// この脚の元本を反対脚の通貨から換算する FX fixing。
    pub fx_reset: Option<IndexId>,
}

/// 校正用の標準スワップで用いるクーポン。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum SwapCouponConvention {
    /// 単利の固定クーポン。
    Fixed,
    /// ターム金利。
    Term {
        /// ターム金利の指標。
        index: IndexId,
    },
    /// 複利する翌日物。
    Compounded(OvernightLegConvention),
}

/// FX フォワードのアウトライトの向きと満期の営業日規則。
///
/// スポットラグ、両通貨の休日の和の暦と受渡の時刻は、ペアの [`FxConvention`] から読む。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FxForwardConvention {
    /// レートの向き。
    pub pair: CcyPair,
    /// 満期に適用する営業日規則。
    pub bdc: BusinessDayConvention,
}

/// 預金の開始日・満期日と利息を決める慣行。
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DepositConvention {
    /// 約定日から開始日までのずれ。
    pub spot_lag: DateOffset,
    /// 利息の日数基準。
    pub day_count: DayCount,
    /// 開始と満期を調整する営業日規則。
    pub bdc: BusinessDayConvention,
    /// 満期日の月末規則。
    pub end_of_month: EndOfMonthRule,
    /// 開始と満期のカレンダー。
    pub calendar: CalendarSpec,
}

/// OIS の日程と二つの脚の慣行。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OisConvention {
    /// 約定日から開始日までのずれ。
    pub spot_lag: DateOffset,
    /// 固定脚の慣行。
    pub fixed: FixedLegConvention,
    /// 翌日物の脚の慣行。
    pub float: OvernightLegConvention,
    /// 支払い日の算出規則。
    pub payment: PaymentRule,
}

/// 固定脚の期間と利息の慣行。
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FixedLegConvention {
    /// 固定脚の頻度。
    pub frequency: Frequency,
    /// 利息の日数基準。
    pub day_count: DayCount,
    /// 期間の境界に適用する営業日調整。
    pub bdc: BusinessDayConvention,
    /// 固定脚のカレンダー。
    pub calendar: CalendarSpec,
}

/// 翌日物の脚で読む指標と観測の規則。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OvernightLegConvention {
    /// 翌日物金利の指標。
    pub index: IndexId,
    /// 期間中の観測と集約の規則。
    pub observation: ObservationRule,
}

/// 参照先を解決する前の参照データ。
#[derive(Debug, Clone, Default)]
pub struct RefDataParts {
    /// 通貨の定義。
    pub currencies: BTreeMap<Currency, CurrencyDefinition>,
    /// 基本カレンダーの定義。
    pub calendars: BTreeMap<CalendarId, BaseCalendar>,
    /// 地域別の UTC オフセットの遷移表。
    pub zones: BTreeMap<TzId, ZoneRules>,
    /// 指標の定義。
    pub indices: BTreeMap<IndexId, IndexDefinition>,
    /// 元の指標ごとのフォールバック定義。
    pub fallbacks: BTreeMap<IndexId, FallbackDefinition>,
    /// 商品の慣行。
    pub conventions: BTreeMap<ConventionId, Convention>,
    /// 通貨ペアの現物の受渡慣行。キーは慣行の `pair`（市場の標準の向き）。
    pub fx_conventions: BTreeMap<CcyPair, FxConvention>,
    /// 証券の定義。
    pub securities: BTreeMap<SecurityId, SecurityDefinition>,
    /// コモディティの定義。
    pub commodities: BTreeMap<CommodityId, CommodityDefinition>,
}

/// 参照を解決済みの参照データのスナップショット。
#[derive(Debug, Clone)]
pub struct RefData(Arc<Resolved>);

/// 単一の参照データに対して解決したカレンダーを、複数の日程や観測期間で共有する。
#[derive(Debug)]
pub struct CalendarCache<'a> {
    refdata: &'a RefData,
    calendars: BTreeMap<CalendarSpec, Calendar>,
}

impl CalendarCache<'_> {
    pub(crate) fn resolve(&mut self, spec: &CalendarSpec) -> Result<Calendar> {
        Ok(match self.calendars.entry(spec.clone()) {
            Entry::Occupied(entry) => entry.get().clone(),
            Entry::Vacant(entry) => entry.insert(self.refdata.calendar(spec)?).clone(),
        })
    }
}

#[derive(Debug)]
struct Resolved {
    parts: RefDataParts,
    calendars: BTreeMap<CalendarSpec, Calendar>,
    publications: BTreeMap<IndexId, ResolvedPublication>,
    rate_futures: BTreeMap<IndexId, (ConventionId, crate::fixing::YearMonth)>,
}

/// 参照データの照合でエラーが見つかった入力の項目。
#[derive(Debug)]
pub enum RefDataEntry {
    /// 通貨。
    Currency(Currency),
    /// 証券。
    Security(SecurityId),
    /// コモディティ。
    Commodity(CommodityId),
    /// 商品の慣行。
    Convention(ConventionId),
    /// 指標。
    Index(IndexId),
    /// フォールバックの元の指標。
    Fallback(IndexId),
    /// 通貨ペアの受渡慣行。
    FxConvention(CcyPair),
}

impl RefData {
    /// 各定義に必要な参照を解決する。同じカレンダー式の前計算は共有する。
    ///
    /// # Errors
    ///
    /// 参照先の欠損、種類の不一致、構成できないカレンダーは `InvalidData`。
    pub fn new(parts: RefDataParts) -> Result<Self> {
        Self::new_with_diagnostics(parts).map_err(|errors| {
            errors.into_iter().next().map_or_else(
                || {
                    Error::Internal(
                        "reference data construction failed without a diagnostic".into(),
                    )
                },
                |(_, error)| error,
            )
        })
    }

    /// 参照データを組み立て、照合に失敗した項目を入力の ID とともに返す。
    ///
    /// # Errors
    ///
    /// 参照先が無いか、定義同士が矛盾する場合。独立に照合できる各項目のエラーを返す。
    pub fn new_with_diagnostics(
        parts: RefDataParts,
    ) -> std::result::Result<Self, Vec<(RefDataEntry, Error)>> {
        let mut calendars = BTreeMap::new();
        let mut use_calendar = |spec: &CalendarSpec| -> Result<Calendar> {
            Ok(match calendars.entry(spec.clone()) {
                Entry::Vacant(entry) => entry
                    .insert(Calendar::compile(spec, &parts.calendars)?)
                    .clone(),
                Entry::Occupied(entry) => entry.get().clone(),
            })
        };
        let mut errors = Vec::new();
        for (&code, definition) in &parts.currencies {
            let result = (|| {
                check_key(&code, &definition.code, "currency")?;
                lookup(&parts.zones, &definition.settlement_time.zone, "zone")?;
                use_calendar(&CalendarSpec::Base(definition.settlement_calendar.clone()))?;
                Ok(())
            })();
            if let Err(error) = result {
                errors.push((RefDataEntry::Currency(code), error));
            }
        }
        validate_assets(&parts, &mut errors);
        for (&pair, convention) in &parts.fx_conventions {
            if let Err(error) = validate_fx_convention(pair, convention, &parts, &mut use_calendar)
            {
                errors.push((RefDataEntry::FxConvention(pair), error));
            }
        }
        let mut rate_futures = BTreeMap::new();
        for (id, convention) in &parts.conventions {
            if let Err(error) = validate_convention(convention, &parts, &mut use_calendar)
                .and_then(|()| rate_future::register(id, convention, &mut rate_futures))
            {
                errors.push((RefDataEntry::Convention(id.clone()), error));
            }
        }
        let mut publications = BTreeMap::new();
        for (id, definition) in &parts.indices {
            let result = (|| {
                check_key(id, &definition.id, "index")?;
                validate_index_kind(&definition.kind, &parts, &mut use_calendar)?;
                use_calendar(&definition.calendar)?;
                let publication = resolve_publication(definition, &parts.zones, &mut use_calendar)?;
                if let Some(time) = &definition.observation_time {
                    if definition.kind.is_credit()
                        || matches!(definition.publication, PublicationRule::Intraday { .. })
                    {
                        return Err(Error::InvalidData(format!(
                            "index {id} does not use a dated observation time"
                        )));
                    }
                    lookup(&parts.zones, &time.zone, "observation zone")?;
                }
                Ok(publication)
            })();
            match result {
                Ok(publication) => {
                    publications.insert(id.clone(), publication);
                }
                Err(error) => errors.push((RefDataEntry::Index(id.clone()), error)),
            }
        }
        for (id, definition) in &parts.fallbacks {
            let result = check_key(id, &definition.index, "fallback").and_then(|()| {
                fallback::validate(id, &definition.replacement, &parts, &mut use_calendar)
            });
            if let Err(error) = result {
                errors.push((RefDataEntry::Fallback(id.clone()), error));
            }
        }
        if !errors.is_empty() {
            return Err(errors);
        }
        Ok(Self(Arc::new(Resolved {
            parts,
            calendars,
            publications,
            rate_futures,
        })))
    }

    /// 通貨の定義を取得する。
    ///
    /// # Errors
    ///
    /// 通貨が未登録なら `InvalidData`。
    pub fn currency(&self, code: Currency) -> Result<&CurrencyDefinition> {
        lookup(&self.0.parts.currencies, &code, "currency")
    }

    /// 指標の定義を取得する。
    ///
    /// # Errors
    ///
    /// 指標が未登録なら `InvalidData`。
    pub fn index(&self, id: &IndexId) -> Result<&IndexDefinition> {
        lookup(&self.0.parts.indices, id, "index")
    }

    /// 指標に既定のフォールバックがあれば返す。
    #[must_use]
    pub fn fallback(&self, id: &IndexId) -> Option<&FallbackDefinition> {
        self.0.parts.fallbacks.get(id)
    }

    /// 取引の個別定義が参照する指標、観測規則、カレンダーを照合する。
    ///
    /// # Errors
    ///
    /// 指標や規則の組み合わせが不正か、参照先が無い場合は `InvalidData`。
    pub fn validate_fallback_replacement(
        &self,
        source: &IndexId,
        replacement: &FallbackReplacement,
    ) -> Result<()> {
        fallback::validate(source, replacement, &self.0.parts, &mut |spec| {
            self.calendar(spec)
        })
    }

    /// ターム金利の観測日が指す預金期間の現地日付。
    ///
    /// # Errors
    ///
    /// 指標の種類、カレンダー、テナーまたは日付が不正な場合。
    pub fn term_rate_dates(
        &self,
        id: &IndexId,
        date: Date,
    ) -> Result<(DateInterval, Vec<Warning>)> {
        let definition = self.index(id)?;
        let IndexKind::TermRate {
            tenor,
            value_offset,
            bdc,
            end_of_month,
            ..
        } = definition.kind
        else {
            return Err(Error::InvalidData(format!("{id} is not a term rate")));
        };
        let calendar = self.calendar(&definition.calendar)?;
        let rule = AdvanceRule {
            end_of_month,
            bdc,
            calendar: &calendar,
        };
        money_market_period(date, tenor, value_offset, &rule)
    }

    /// 証券の定義を取得する。
    ///
    /// # Errors
    ///
    /// 証券が未登録なら `InvalidData`。
    pub fn security(&self, id: &SecurityId) -> Result<&SecurityDefinition> {
        lookup(&self.0.parts.securities, id, "security")
    }

    /// コモディティの定義を取得する。
    ///
    /// # Errors
    ///
    /// コモディティが未登録なら `InvalidData`。
    pub fn commodity(&self, id: &CommodityId) -> Result<&CommodityDefinition> {
        lookup(&self.0.parts.commodities, id, "commodity")
    }

    /// 商品の慣行を取得する。
    ///
    /// # Errors
    ///
    /// 慣行が未登録なら `InvalidData`。
    pub fn convention(&self, id: &ConventionId) -> Result<&Convention> {
        lookup(&self.0.parts.conventions, id, "convention")
    }

    /// 定義済みの基本カレンダーを使って、式の休日集合を求める。
    ///
    /// # Errors
    ///
    /// 参照先の欠損または不正な式なら `InvalidData`。
    pub fn calendar(&self, spec: &CalendarSpec) -> Result<Calendar> {
        self.0.calendars.get(spec).map_or_else(
            || Calendar::compile(spec, &self.0.parts.calendars),
            |calendar| Ok(calendar.clone()),
        )
    }

    /// 商品ごとのカレンダー解決結果を再利用するためのキャッシュを開く。
    #[must_use]
    pub fn calendar_cache(&self) -> CalendarCache<'_> {
        CalendarCache {
            refdata: self,
            calendars: BTreeMap::new(),
        }
    }

    /// スケジュールや観測の計算に渡せる基本カレンダーの表。
    #[must_use]
    pub fn calendars(&self) -> &BTreeMap<CalendarId, BaseCalendar> {
        &self.0.parts.calendars
    }

    /// 現地時刻の解決に使う地域ごとの遷移表。
    #[must_use]
    pub fn zones(&self) -> &BTreeMap<TzId, ZoneRules> {
        &self.0.parts.zones
    }

    /// 地域の遷移表を返す。
    ///
    /// 地域を受け取る呼び出し（評価時点の解決、締め時刻、通貨の決済時刻、観測時刻）は、
    /// すべてここで引いて同じ種類のエラーにする。
    ///
    /// # Errors
    ///
    /// 登録されていない地域なら `InvalidData`。
    pub fn zone(&self, id: &TzId) -> Result<&ZoneRules> {
        lookup(&self.0.parts.zones, id, "zone")
    }

    /// 通貨ペアの受渡慣行を、向きによらずに返す。慣行の `pair` が市場の標準の向き。
    ///
    /// # Errors
    ///
    /// どちらの向きの慣行も無ければ `InvalidData`。
    pub fn fx_convention(&self, pair: CcyPair) -> Result<&FxConvention> {
        fx_convention(&self.0.parts, pair)
    }

    /// 通貨ペアの受渡日を数える、両通貨の決済カレンダーの和の式。向きによらず同じ式を返す。
    ///
    /// # Errors
    ///
    /// どちらかの通貨が未登録なら `InvalidData`。
    pub fn fx_calendar_spec(&self, pair: CcyPair) -> Result<CalendarSpec> {
        fx_calendar_spec(&self.0.parts, pair)
    }

    /// 指標の観測点が公表される瞬間と、カレンダー外挿の警告を返す。
    ///
    /// # Errors
    ///
    /// 指標や観測点の種類が不明か、観測点の公表時刻を決められない場合。
    pub fn publication(
        &self,
        id: &IndexId,
        point: ObservationPoint,
    ) -> Result<(Instant, Vec<Warning>)> {
        let definition = self.index(id)?;
        if definition.kind.is_credit() {
            return Err(incompatible_point(id));
        }
        lookup(&self.0.publications, id, "publication")?.at(definition, point)
    }

    /// 観測点が表す期間を UTC の瞬間に解決する。
    ///
    /// # Errors
    ///
    /// 未登録の指標、観測点の不一致、観測時刻の欠損、暦や地域の範囲外はエラー。
    pub fn observation_period(
        &self,
        id: &IndexId,
        point: ObservationPoint,
    ) -> Result<(ObservationPeriod, Vec<Warning>)> {
        self.validate_point(id, point)?;
        if let ObservationPoint::Instant(at) = point {
            return Ok((ObservationPeriod::at(at), Vec::new()));
        }
        let definition = self.index(id)?;
        let time = definition.observation_time.as_ref().ok_or_else(|| {
            Error::InvalidData(format!("missing observation time for index {id}"))
        })?;
        let zone = lookup(&self.0.parts.zones, &time.zone, "observation zone")?;
        let mut warnings = Vec::new();
        let (start, end) = match (point, &definition.kind) {
            (ObservationPoint::Date(date), IndexKind::OvernightRate { .. }) => {
                let calendar = self.calendar(&definition.calendar)?;
                (date, shift_business(date, 1, &calendar, &mut warnings)?)
            }
            (ObservationPoint::Date(date), IndexKind::TermRate { .. }) => {
                let (period, notices) = self.term_rate_dates(id, date)?;
                warnings.extend(notices);
                (period.start(), period.end())
            }
            (ObservationPoint::Date(date), _) => (date, date),
            (ObservationPoint::Month(month), _) => {
                let start = month.first_day();
                (start, start.month_end()?.checked_add_days(1)?)
            }
            (ObservationPoint::Instant(_), _) => return Err(incompatible_point(id)),
        };
        let start = zone.resolve_convention(time.on_date(start))?;
        let end = zone.resolve_convention(time.on_date(end))?;
        Ok((
            ObservationPeriod::new(start, end)?,
            unique_warnings(warnings),
        ))
    }

    /// 信用窓 (start, end] の値が判明する瞬間と、カレンダー外挿の警告を返す。
    ///
    /// # Errors
    ///
    /// 指標が信用事由・回収率の指標でないか、信用窓が空か、公表の瞬間を解決できない場合。
    pub fn credit_publication(
        &self,
        id: &IndexId,
        window: DateInterval,
    ) -> Result<(Instant, Vec<Warning>)> {
        if !self.index(id)?.kind.is_credit() {
            return Err(Error::InvalidData(format!(
                "index {id} does not observe a credit window"
            )));
        }
        lookup(&self.0.publications, id, "publication")?.on_date(credit_window_end(window)?)
    }

    pub(crate) fn validate_point(&self, id: &IndexId, point: ObservationPoint) -> Result<()> {
        if self.index(id)?.kind.is_credit()
            || !lookup(&self.0.publications, id, "publication")?.accepts(point)
        {
            return Err(incompatible_point(id));
        }
        Ok(())
    }
}

fn lookup<'a, K: Ord + std::fmt::Display, V>(
    map: &'a BTreeMap<K, V>,
    id: &K,
    kind: &str,
) -> Result<&'a V> {
    map.get(id)
        .ok_or_else(|| Error::InvalidData(format!("unknown {kind} {id}")))
}

/// 証券とコモディティの定義の参照を照合し、失敗した項目を `errors` に足す。
fn validate_assets(parts: &RefDataParts, errors: &mut Vec<(RefDataEntry, Error)>) {
    for (id, security) in &parts.securities {
        let result = (|| {
            check_key(id, &security.id, "security")?;
            lookup(&parts.currencies, &security.currency, "currency")?;
            if security.kind == SecurityKind::Bond && security.issuer.is_none() {
                return Err(Error::InvalidData(format!("bond {id} requires an issuer")));
            }
            Ok(())
        })();
        if let Err(error) = result {
            errors.push((RefDataEntry::Security(id.clone()), error));
        }
    }
    for (id, commodity) in &parts.commodities {
        let result = (|| {
            check_key(id, &commodity.id, "commodity")?;
            lookup(&parts.currencies, &commodity.currency, "currency")?;
            if commodity.unit.trim().is_empty() {
                return Err(Error::InvalidData(format!(
                    "commodity {id} requires a physical unit"
                )));
            }
            Ok(())
        })();
        if let Err(error) = result {
            errors.push((RefDataEntry::Commodity(id.clone()), error));
        }
    }
}

/// 表のキーと向き、逆向きの重複、受渡時刻の地域、両通貨の休日の和の暦を照合する。
fn validate_fx_convention(
    pair: CcyPair,
    convention: &FxConvention,
    parts: &RefDataParts,
    use_calendar: &mut impl FnMut(&CalendarSpec) -> Result<Calendar>,
) -> Result<()> {
    check_key(&pair, &convention.pair, "FX convention")?;
    // 同じペアの慣行を二つの向きで持つと、どちらを正本にするかが決まらない。
    if parts.fx_conventions.contains_key(&pair.inverse()) {
        return Err(Error::InvalidData(format!(
            "FX conventions for both {pair} and {}",
            pair.inverse()
        )));
    }
    lookup(
        &parts.zones,
        &convention.settlement.zone,
        "FX settlement zone",
    )?;
    use_calendar(&fx_calendar_spec(parts, pair)?)?;
    Ok(())
}

fn fx_convention(parts: &RefDataParts, pair: CcyPair) -> Result<&FxConvention> {
    parts
        .fx_conventions
        .get(&pair)
        .or_else(|| parts.fx_conventions.get(&pair.inverse()))
        .ok_or_else(|| Error::InvalidData(format!("unknown FX convention {pair}")))
}

fn fx_calendar_spec(parts: &RefDataParts, pair: CcyPair) -> Result<CalendarSpec> {
    let [base, quote] = [pair.base(), pair.quote()].map(|ccy| {
        lookup(&parts.currencies, &ccy, "currency")
            .map(|definition| &definition.settlement_calendar)
    });
    Ok(fx_convention::pair_calendar(base?, quote?))
}

fn check_key<T: Eq + std::fmt::Display>(key: &T, id: &T, kind: &str) -> Result<()> {
    if key != id {
        return Err(Error::InvalidData(format!(
            "{kind} key {key} does not match {id}"
        )));
    }
    Ok(())
}

fn validate_day_count(
    basis: &DayCount,
    use_calendar: &mut impl FnMut(&CalendarSpec) -> Result<Calendar>,
) -> Result<()> {
    if let DayCount::Bus252 { calendar } = basis {
        use_calendar(&CalendarSpec::Base(calendar.clone()))?;
    }
    Ok(())
}

fn validate_index_kind(
    kind: &IndexKind,
    parts: &RefDataParts,
    use_calendar: &mut impl FnMut(&CalendarSpec) -> Result<Calendar>,
) -> Result<()> {
    let currency = |code| lookup(&parts.currencies, &code, "currency").map(|_| ());
    let security = |id| lookup(&parts.securities, id, "security");
    match kind {
        IndexKind::TermRate { ccy, day_count, .. }
        | IndexKind::OvernightRate { ccy, day_count } => {
            currency(*ccy)?;
            validate_day_count(day_count, use_calendar)?;
        }
        IndexKind::CompoundedRateIndex { overnight } => validate_overnight(parts, overnight)?,
        IndexKind::SwapRate {
            ccy, convention, ..
        } => {
            currency(*ccy)?;
            let Convention::Ois(ois) = lookup(&parts.conventions, convention, "convention")? else {
                return Err(Error::InvalidData(format!(
                    "{convention} is not a swap convention"
                )));
            };
            let IndexKind::OvernightRate { ccy: float_ccy, .. } =
                &lookup(&parts.indices, &ois.float.index, "index")?.kind
            else {
                return Err(Error::InvalidData(format!(
                    "index {} is not an overnight rate",
                    ois.float.index
                )));
            };
            check_key(ccy, float_ccy, "swap convention currency")?;
        }
        IndexKind::FxFixing { pair } => {
            currency(pair.base())?;
            currency(pair.quote())?;
            // FX fixing の受渡日はペアのスポットラグで決まるので、慣行の無いペアは評価できない。
            fx_convention(parts, *pair)?;
        }
        IndexKind::AssetPrice { asset, ccy, .. } => {
            currency(*ccy)?;
            match asset {
                Asset::Cash(ccy) => currency(*ccy)?,
                Asset::Security(id) => {
                    security(id)?;
                }
                Asset::Commodity(id) => {
                    lookup(&parts.commodities, id, "commodity")?;
                }
            }
        }
        IndexKind::EquityIndexLevel { ccy, unit } => {
            currency(*ccy)?;
            validate_security_kind(security(unit)?, SecurityKind::IndexUnit)?;
        }
        IndexKind::FuturesSettlement { contract, ccy } => {
            currency(*ccy)?;
            validate_security_kind(security(contract)?, SecurityKind::ListedFuture)?;
        }
        IndexKind::RateFuturesSettlement {
            contract,
            underlying,
        } => {
            validate_security_kind(security(contract)?, SecurityKind::ListedFuture)?;
            if !matches!(
                lookup(&parts.indices, underlying, "index")?.kind,
                IndexKind::TermRate { .. } | IndexKind::OvernightRate { .. }
            ) {
                return Err(Error::InvalidData(format!(
                    "{underlying} is not an underlying rate index"
                )));
            }
        }
        IndexKind::InflationIndex { .. } => {}
        IndexKind::CreditEvent { credit } | IndexKind::Recovery { credit } => currency(credit.ccy)?,
        IndexKind::Dividend { asset, ccy } => {
            currency(*ccy)?;
            validate_security_kind(security(asset)?, SecurityKind::Share)?;
        }
    }
    Ok(())
}

fn validate_security_kind(definition: &SecurityDefinition, kind: SecurityKind) -> Result<()> {
    if definition.kind != kind {
        return Err(Error::InvalidData(format!(
            "unexpected security kind for {}",
            definition.id
        )));
    }
    Ok(())
}

fn validate_overnight(parts: &RefDataParts, id: &IndexId) -> Result<()> {
    if !matches!(
        lookup(&parts.indices, id, "index")?.kind,
        IndexKind::OvernightRate { .. }
    ) {
        return Err(Error::InvalidData(format!(
            "index {id} is not an overnight rate"
        )));
    }
    Ok(())
}

fn validate_convention(
    convention: &Convention,
    parts: &RefDataParts,
    use_calendar: &mut impl FnMut(&CalendarSpec) -> Result<Calendar>,
) -> Result<()> {
    match convention {
        Convention::Deposit(deposit) => {
            use_calendar(&deposit.calendar)?;
            validate_day_count(&deposit.day_count, use_calendar)?;
        }
        Convention::Ois(ois) => {
            use_calendar(&ois.fixed.calendar)?;
            use_calendar(&ois.payment.calendar)?;
            validate_day_count(&ois.fixed.day_count, use_calendar)?;
            if let Some(settlement) = &ois.payment.settlement {
                lookup(&parts.zones, &settlement.zone, "zone")?;
            }
            validate_overnight(parts, &ois.float.index)?;
            match ois.float.observation.days() {
                ObservationDays::BusinessDays(calendar)
                | ObservationDays::AllDays {
                    fixing_calendar: calendar,
                }
                | ObservationDays::Endpoints(EndpointDates::Shifted { calendar, .. }) => {
                    use_calendar(calendar)?;
                }
                ObservationDays::Months | ObservationDays::Endpoints(EndpointDates::Accrual) => {}
            }
        }
        Convention::RateFuture(future) => future.validate(parts, use_calendar)?,
        Convention::Fra(fra) => {
            if !matches!(
                lookup(&parts.indices, &fra.index, "FRA index")?.kind,
                IndexKind::TermRate { .. }
            ) {
                return Err(Error::InvalidData("FRA requires a term rate index".into()));
            }
        }
        Convention::Irs(swap) | Convention::BasisSwap(swap) | Convention::XccySwap(swap) => {
            validate_swap_convention(convention, swap, parts, use_calendar)?;
        }
        Convention::FxForward(fx) => {
            for ccy in [fx.pair.base(), fx.pair.quote()] {
                if lookup(&parts.currencies, &ccy, "FX currency")?.deliverability
                    != Deliverability::Deliverable
                {
                    return Err(Error::InvalidData(format!("{ccy} is not deliverable")));
                }
            }
            fx_convention(parts, fx.pair)?;
        }
    }
    Ok(())
}

fn validate_swap_convention(
    convention: &Convention,
    swap: &SwapConvention,
    parts: &RefDataParts,
    use_calendar: &mut impl FnMut(&CalendarSpec) -> Result<Calendar>,
) -> Result<()> {
    for leg in &swap.legs {
        lookup(&parts.currencies, &leg.currency, "swap leg currency")?;
        use_calendar(&leg.periods.calendar)?;
        use_calendar(&leg.payment.calendar)?;
        validate_day_count(&leg.periods.day_count, use_calendar)?;
        if let Some(settlement) = &leg.payment.settlement {
            lookup(&parts.zones, &settlement.zone, "swap payment zone")?;
        }
        match &leg.coupon {
            SwapCouponConvention::Fixed => {}
            SwapCouponConvention::Term { index } => {
                if !matches!(lookup(&parts.indices, index, "swap term index")?.kind, IndexKind::TermRate { ccy, .. } if ccy == leg.currency)
                {
                    return Err(Error::InvalidData(format!(
                        "term index {index} does not match swap leg"
                    )));
                }
            }
            SwapCouponConvention::Compounded(overnight) => {
                validate_overnight(parts, &overnight.index)?;
                if !matches!(lookup(&parts.indices, &overnight.index, "swap overnight index")?.kind, IndexKind::OvernightRate { ccy, .. } if ccy == leg.currency)
                {
                    return Err(Error::InvalidData(
                        "overnight index currency differs from swap leg currency".into(),
                    ));
                }
            }
        }
    }
    let [first, second] = &swap.legs;
    let same_currency = first.currency == second.currency;
    let first_fixed = matches!(first.coupon, SwapCouponConvention::Fixed);
    let second_fixed = matches!(second.coupon, SwapCouponConvention::Fixed);
    let resets = swap
        .legs
        .iter()
        .filter(|leg| leg.fx_reset.is_some())
        .count();
    let valid = match convention {
        Convention::Irs(_) => same_currency && first_fixed != second_fixed && resets == 0,
        Convention::BasisSwap(_) => same_currency && !first_fixed && !second_fixed && resets == 0,
        Convention::XccySwap(_) => !same_currency && !first_fixed && !second_fixed && resets == 1,
        _ => false,
    };
    if !valid {
        return Err(Error::InvalidData(
            "swap legs do not match the convention kind".into(),
        ));
    }
    for (leg, other) in [(first, second), (second, first)] {
        if let Some(index) = &leg.fx_reset {
            let pair = CcyPair::new(leg.currency, other.currency)?;
            if !matches!(lookup(&parts.indices, index, "FX reset index")?.kind, IndexKind::FxFixing { pair: fixing } if fixing == pair || fixing == pair.inverse())
            {
                return Err(Error::InvalidData(format!(
                    "FX reset {index} does not match the two swap legs"
                )));
            }
        }
    }
    Ok(())
}

fn resolve_publication(
    definition: &IndexDefinition,
    zones: &BTreeMap<TzId, ZoneRules>,
    use_calendar: &mut impl FnMut(&CalendarSpec) -> Result<Calendar>,
) -> Result<ResolvedPublication> {
    if !definition.kind.accepts_publication(&definition.publication) {
        return Err(Error::InvalidData(format!(
            "incompatible publication rule for {}",
            definition.id
        )));
    }
    Ok(match &definition.publication {
        PublicationRule::Daily(daily) => ResolvedPublication::Daily {
            business_days: daily.business_days,
            calendar: use_calendar(&daily.calendar)?,
            local_time: LocalEventTime {
                time: daily.time,
                zone: daily.zone.clone(),
            },
            zone: lookup(zones, &daily.zone, "zone")?.clone(),
        },
        PublicationRule::Intraday { delay_micros } => {
            ResolvedPublication::Intraday(i64::try_from(*delay_micros).map_err(|_| {
                Error::InvalidData(format!(
                    "publication delay for {} is out of range",
                    definition.id
                ))
            })?)
        }
        PublicationRule::Monthly { .. } => ResolvedPublication::Monthly,
    })
}
`````

`````
crates/rhodium-core/src/refdata/fallback.rs
use crate::{
    Error, Result,
    calendar::{Calendar, CalendarSpec},
    date::Date,
    day_count::DayCount,
    id::IndexId,
    observation::{Aggregation, ObservationDays, ObservationRule},
    quantity::Spread,
};

use super::{IndexKind, PublicationRule, RefDataParts, lookup, validate_day_count};

/// 指標に既定のフォールバックとその契約上の適用条件を与える。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FallbackDefinition {
    /// 置き換えるターム金利の指標。
    pub index: IndexId,
    /// フォールバックを開始する事由。
    pub trigger: FallbackTrigger,
    /// 代替の観測とスプレッド。
    pub replacement: FallbackReplacement,
    /// 取引が満たす条件の候補。
    pub applies: Vec<FallbackApplicability>,
}

/// 代替の翌日物金利からターム金利を求める条件。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FallbackReplacement {
    /// 複利にする翌日物金利の指標。
    pub index: IndexId,
    /// 観測日、複利と丸めの規則。
    pub observation: ObservationRule,
    /// 観測期間の年率化に使う日数基準。
    pub annualization: DayCount,
    /// 代替する率に適用するスプレッド。
    #[serde(deserialize_with = "scaled_spread")]
    pub spread: Spread,
}

/// 指標の状態に対する、契約上のフォールバック開始条件。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FallbackTrigger {
    /// 公表停止。
    Cessation,
    /// 公表停止または代表性の喪失。
    CessationOrNonRepresentative,
}

/// 取引に既定のフォールバックを適用する条件。
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum FallbackApplicability {
    /// 確認書が指定する定義書の版と適用開始日。
    Definitions {
        /// 定義書の版。
        version: String,
        /// この日以後に約定した取引に適用する。
        traded_on_or_after: Date,
    },
    /// 両当事者が批准したプロトコル。
    Protocol {
        /// プロトコルの識別子。
        protocol: String,
    },
    /// この指標についての個別合意。
    Bilateral,
}

fn scaled_spread<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> std::result::Result<Spread, D::Error> {
    use serde::Deserialize;
    String::deserialize(deserializer)?
        .parse()
        .map_err(serde::de::Error::custom)
}

pub(super) fn validate(
    source: &IndexId,
    replacement: &FallbackReplacement,
    parts: &RefDataParts,
    use_calendar: &mut impl FnMut(&CalendarSpec) -> Result<Calendar>,
) -> Result<()> {
    let original = lookup(&parts.indices, source, "index")?;
    let alternate = lookup(&parts.indices, &replacement.index, "index")?;
    match (&original.kind, &alternate.kind) {
        (IndexKind::TermRate { ccy: left, .. }, IndexKind::OvernightRate { ccy: right, .. })
            if left == right
                && matches!(original.publication, PublicationRule::Daily(_))
                && matches!(alternate.publication, PublicationRule::Daily(_))
                && alternate.observation_time.is_some() => {}
        _ => {
            return Err(Error::InvalidData(format!(
                "fallback for {source} requires daily term and overnight rates in the same currency"
            )));
        }
    }
    let (ObservationDays::BusinessDays(calendar), Aggregation::Compound { .. }) = (
        replacement.observation.days(),
        replacement.observation.aggregation(),
    ) else {
        return Err(Error::InvalidData(
            "fallback requires business-day compounding".into(),
        ));
    };
    use_calendar(calendar)?;
    validate_day_count(&replacement.annualization, use_calendar)
}
`````

`````
crates/rhodium-core/src/refdata/fx_convention.rs
//! 通貨ペアの現物の受渡慣行。

use crate::{
    calendar::CalendarSpec,
    fx::CcyPair,
    id::CalendarId,
    tenor::{DateOffset, TenorUnit},
    time::LocalEventTime,
};

/// 現物 FX のスポットラグ。両通貨の決済カレンダーの和の営業日で数える。
///
/// 閉じた enum にして、スポットまでの短期の quote の連鎖（T+1 は ON だけ、T+2 は TN と ON）を
/// 網羅して分岐し、市場に無い T+0 や T+5 の慣行を読み込みの境界で拒否する。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FxSpotLag {
    /// 翌営業日がスポット（USDCAD など）。
    OneBusinessDay,
    /// 二営業日後がスポット（USDJPY、EURUSD など）。
    TwoBusinessDays,
}

impl FxSpotLag {
    /// スポットまでの営業日数。
    #[must_use]
    pub const fn business_days(self) -> i16 {
        match self {
            Self::OneBusinessDay => 1,
            Self::TwoBusinessDays => 2,
        }
    }

    /// 約定日からスポット日へ進めるずれ。
    #[must_use]
    pub const fn offset(self) -> DateOffset {
        DateOffset::new(self.business_days(), TenorUnit::BusinessDay)
    }
}

/// 通貨ペアの現物の受渡慣行。
///
/// 評価時点の換算レート X(0) の短期の連鎖、FX fixing の受渡日、`MtM` XCCY のリセットの
/// 観測日、FX フォワードの校正商品の日程は、すべてこの慣行から決める（正本は一つ）。
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FxConvention {
    /// 市場の標準の向き（USDJPY、EURUSD）。pip とポイントの quote はこの向きで読む。
    pub pair: CcyPair,
    /// スポットラグ。
    pub spot_lag: FxSpotLag,
    /// 各受渡日の受渡の瞬間。この瞬間を過ぎた当日の受渡はできない（同日受渡の締め）。
    pub settlement: LocalEventTime,
    /// 1 pip がレートの 10 の何乗分の一か（USDJPY は 2、EURUSD は 4）。
    pub pip_exponent: u8,
}

impl FxConvention {
    /// 1 pip のレートの大きさ 10^(−`pip_exponent`)。
    #[must_use]
    pub fn pip(&self) -> f64 {
        10_f64.powi(-i32::from(self.pip_exponent))
    }
}

/// 両通貨の決済カレンダーの和。ID の順に並べ、向きの違うペアで同じ式（同じキャッシュ）にする。
pub(super) fn pair_calendar(base: &CalendarId, quote: &CalendarId) -> CalendarSpec {
    if base == quote {
        return CalendarSpec::Base(base.clone());
    }
    let (first, second) = if base < quote {
        (base, quote)
    } else {
        (quote, base)
    };
    CalendarSpec::Union(vec![
        CalendarSpec::Base(first.clone()),
        CalendarSpec::Base(second.clone()),
    ])
}
`````

`````
crates/rhodium-core/src/refdata/index.rs
use std::collections::BTreeMap;

use crate::{
    Error, Result, Warning,
    calendar::{BusinessDayConvention, EndOfMonthRule},
    calendar::{Calendar, CalendarSpec, shift_business, unique_warnings},
    currency::{Asset, Currency},
    date::{Date, DateInterval},
    day_count::DayCount,
    fixing::{ObservationPoint, YearMonth},
    fx::CcyPair,
    id::{ConventionId, CreditKey, IndexId, SecurityId},
    tenor::{DateOffset, Tenor},
    time::{Instant, LocalEventTime, LocalTime, TzId, ZoneRules},
};

/// 確認書が参照する指標の定義。参照先は `RefData` の構築時に解決する。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IndexDefinition {
    /// 指標の ID。
    pub id: IndexId,
    /// 値の意味。
    pub kind: IndexKind,
    /// 観測日のカレンダー。
    pub calendar: CalendarSpec,
    /// 指標の公表規則。
    pub publication: PublicationRule,
    /// 日次・月次の観測期間の境界を決める現地時刻。公表時刻とは独立に指定する。
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub observation_time: Option<LocalEventTime>,
    /// 指標の状態。
    #[serde(default)]
    pub status: IndexStatus,
}

/// 指標の代表性と公表の状態。
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum IndexStatus {
    /// 通常どおり公表される。
    #[default]
    Active,
    /// 指標としての代表性を失った。
    NonRepresentative {
        /// 代表性を失った日。
        since: Date,
    },
    /// 公表が終了した。
    Ceased {
        /// 公表を終了した日。
        since: Date,
    },
}

/// 指標が表す値の種類と、その意味を決める属性。
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum IndexKind {
    /// 預金期間の単利の率。
    TermRate {
        /// 預金の通貨。
        ccy: Currency,
        /// 率の日数基準。
        day_count: DayCount,
        /// 預金期間の長さ。
        tenor: Tenor,
        /// 観測日から預金開始日までのずれ。
        value_offset: DateOffset,
        /// 預金満期日の営業日調整。
        bdc: BusinessDayConvention,
        /// 預金満期日の月末規則。
        end_of_month: EndOfMonthRule,
    },
    /// 翌日物の率。
    OvernightRate {
        /// 預金の通貨。
        ccy: Currency,
        /// 率の日数基準。
        day_count: DayCount,
    },
    /// 翌日物を複利した指数の水準。
    CompoundedRateIndex {
        /// 元の翌日物指標。
        overnight: IndexId,
    },
    /// スワップのパーレート。
    SwapRate {
        /// スワップの通貨。
        ccy: Currency,
        /// スワップの条件を定める慣行。
        convention: ConventionId,
        /// スワップの期間。
        tenor: Tenor,
    },
    /// 指定した向きの FX のスポットレート。
    FxFixing {
        /// レートの通貨ペア。
        pair: CcyPair,
    },
    /// 資産 1 単位あたりの価格。
    AssetPrice {
        /// 価格の対象。
        asset: Asset,
        /// 価格を表す通貨。
        ccy: Currency,
        /// 価格の決め方。
        price: PriceKind,
    },
    /// 株価指数の水準。
    EquityIndexLevel {
        /// 水準を表す通貨。
        ccy: Currency,
        /// 指数の取引単位を表す証券。
        unit: SecurityId,
    },
    /// 先物の日々の清算価格。
    FuturesSettlement {
        /// 先物の識別子。
        contract: SecurityId,
        /// 清算の通貨。
        ccy: Currency,
    },
    /// 短期金利先物の日々の清算価格。
    RateFuturesSettlement {
        /// 先物の識別子。
        contract: SecurityId,
        /// 元の金利指標。
        underlying: IndexId,
    },
    /// 参照月の物価指数の水準。
    InflationIndex {
        /// 指標の対象地域。
        region: String,
        /// 市場の標準の参照月の遅れ。
        lag_months: u8,
        /// 取引に既定で適用する補間。
        interpolation: InflationInterpolation,
        /// 市場の季節性の扱い。
        seasonality: Seasonality,
    },
    /// 日付の窓における信用事由の有無。
    CreditEvent {
        /// 判定の対象となる信用参照。
        credit: CreditKey,
    },
    /// 日付の窓における最初の信用事由の回収率。
    Recovery {
        /// 判定の対象となる信用参照。
        credit: CreditKey,
    },
    /// 権利落ち日の 1 株あたりの配当。
    Dividend {
        /// 配当する証券。
        asset: SecurityId,
        /// 配当の通貨。
        ccy: Currency,
    },
}

/// 資産価格を決める方法。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PriceKind {
    /// 当日の終値。
    Close,
    /// 日次で公表される評価値。
    Assessment,
    /// 日中に観測される価格。
    Intraday,
}

/// 物価指数の参照月間の補間方法。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InflationInterpolation {
    /// 各参照月の値をそのまま使う。
    Flat,
    /// 参照月間を線形補間する。
    Linear,
}

/// 月別の季節性の扱い。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Seasonality {
    /// 季節性を調整しない。
    None,
    /// 月別の季節性を乗法で適用する。
    Multiplicative,
}

/// 観測点から一つの公表の瞬間を決める規則。
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum PublicationRule {
    /// 日中の観測から一定時間後に公表する。
    Intraday {
        /// 観測から公表までの非負のマイクロ秒数。
        delay_micros: u64,
    },
    /// 参照月別の公表日程を使う。
    Monthly {
        /// 参照月ごとの公表の瞬間。
        releases: BTreeMap<YearMonth, Instant>,
    },
    /// 日付を起点に営業日ラグを適用し、現地時刻を解決する。
    #[serde(untagged)]
    Daily(DailyPublication),
}

/// 日次の観測日または信用窓の終了日を起点にする公表規則。
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DailyPublication {
    /// 公表カレンダーで数える営業日数。ゼロは同日公表。
    pub business_days: u16,
    /// 観測日から公表日までの営業日を数えるカレンダー。
    pub calendar: CalendarSpec,
    /// 現地の公表時刻。
    #[serde(deserialize_with = "crate::time::deserialize_event_time")]
    pub time: LocalTime,
    /// 公表時刻を解決する地域。
    pub zone: TzId,
}

#[derive(Debug)]
pub(super) enum ResolvedPublication {
    Daily {
        business_days: u16,
        calendar: Calendar,
        local_time: LocalEventTime,
        zone: ZoneRules,
    },
    Intraday(i64),
    Monthly,
}

impl ResolvedPublication {
    pub(super) fn at(
        &self,
        definition: &IndexDefinition,
        point: ObservationPoint,
    ) -> Result<(Instant, Vec<Warning>)> {
        match (self, &definition.publication, point) {
            (Self::Daily { .. }, PublicationRule::Daily(_), ObservationPoint::Date(date)) => {
                self.on_date(date)
            }
            (
                Self::Intraday(delay),
                PublicationRule::Intraday { .. },
                ObservationPoint::Instant(at),
            ) => {
                let micros = at
                    .unix_micros()
                    .checked_add(*delay)
                    .ok_or_else(|| Error::InvalidInput("publication instant overflows".into()))?;
                Ok((Instant::from_unix_micros(micros)?, Vec::new()))
            }
            (
                Self::Monthly,
                PublicationRule::Monthly { releases },
                ObservationPoint::Month(month),
            ) => {
                let at = releases.get(&month).ok_or_else(|| {
                    Error::MissingData(format!(
                        "missing publication time for {} {month}",
                        definition.id
                    ))
                })?;
                Ok((*at, Vec::new()))
            }
            _ => Err(incompatible_point(&definition.id)),
        }
    }

    pub(super) fn on_date(&self, mut date: Date) -> Result<(Instant, Vec<Warning>)> {
        let Self::Daily {
            business_days,
            calendar,
            local_time,
            zone,
        } = self
        else {
            return Err(Error::Internal(
                "credit publication must be a daily rule".into(),
            ));
        };
        let mut warnings = Vec::new();
        // ラグがゼロなら、起点の日付をそのまま使う。
        if *business_days > 0 {
            date = shift_business(date, i32::from(*business_days), calendar, &mut warnings)?;
        }
        let instant = zone.resolve_convention(local_time.on_date(date))?;
        Ok((instant, unique_warnings(warnings)))
    }

    pub(super) fn accepts(&self, point: ObservationPoint) -> bool {
        matches!(
            (self, point),
            (Self::Daily { .. }, ObservationPoint::Date(_))
                | (Self::Intraday(_), ObservationPoint::Instant(_))
                | (Self::Monthly, ObservationPoint::Month(_))
        )
    }
}

impl IndexKind {
    pub(super) const fn accepts_publication(&self, rule: &PublicationRule) -> bool {
        match self {
            Self::InflationIndex { .. } => matches!(rule, PublicationRule::Monthly { .. }),
            Self::AssetPrice {
                price: PriceKind::Intraday,
                ..
            } => {
                matches!(rule, PublicationRule::Intraday { .. })
            }
            Self::FxFixing { .. } => !matches!(rule, PublicationRule::Monthly { .. }),
            _ => matches!(rule, PublicationRule::Daily(_)),
        }
    }

    pub(super) const fn is_credit(&self) -> bool {
        matches!(self, Self::CreditEvent { .. } | Self::Recovery { .. })
    }
}

pub(super) fn incompatible_point(id: &IndexId) -> Error {
    Error::InvalidData(format!("observation point does not match index {id}"))
}

pub(super) fn credit_window_end(window: DateInterval) -> Result<Date> {
    if window.is_empty() {
        return Err(Error::InvalidInput("credit event window is empty".into()));
    }
    Ok(window.end())
}
`````

`````
crates/rhodium-core/src/refdata/rate_future.rs
//! 限月別の金利先物の清算日程と最終決済の条件。

use std::collections::BTreeMap;

use crate::{
    Error, Result, Warning,
    calendar::{
        AdvanceRule, BusinessDayConvention, Calendar, CalendarSpec, EndOfMonthRule, advance,
    },
    currency::Currency,
    date::{Date, DateInterval},
    day_count::DayCount,
    fixing::{ObservationPoint, YearMonth},
    id::{ConventionId, IndexId},
    tenor::{DateOffset, TenorUnit},
    time::{Instant, LocalEventTime},
};

use super::{Convention, IndexKind, RefData, RefDataParts, lookup};

/// 最終清算値 1 − R の参照金利。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum RateFutureRate {
    /// 一回のターム金利。
    Term {
        /// 原指標の観測日。
        fixing: Date,
    },
    /// 参照期間の日次単利因子を複利する。
    Compounded(DateInterval),
    /// 参照期間の各翌日物を日数係数で加重平均する。
    Averaged(DateInterval),
}

/// 限月の清算指標と最終決済の条件。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RateFutureContract {
    /// 限月固有の日次清算指標。
    pub settlement_index: IndexId,
    /// 最後の清算日。
    pub final_clearing: Date,
    /// 最終清算値の参照金利。
    pub rate: RateFutureRate,
}

/// 取引所が定める日次値洗いの慣行。
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RateFutureConvention {
    /// 価格 1 単位あたりの一枚の受渡額。
    pub multiplier: f64,
    /// 清算日から値洗いの受渡日までのずれ。
    pub settlement_lag: DateOffset,
    /// 値洗いを受け渡す現地時刻。
    pub settlement_time: LocalEventTime,
    /// 取引所の値洗い口座の付利指標。省略時は支払の割引曲線で付利する。
    pub margin_index: Option<IndexId>,
    /// 月ごとの上場契約。
    pub contracts: BTreeMap<YearMonth, RateFutureContract>,
}

/// 日次の清算価格の観測と値洗いの受渡。
#[derive(Debug, Clone, Copy)]
pub struct RateFutureClearing {
    /// 清算指標の観測日。
    pub date: Date,
    /// 清算価格を決める瞬間。
    pub observed_at: Instant,
    /// 清算価格の公表の瞬間。
    pub published_at: Instant,
    /// 値洗いの受渡の瞬間。
    pub settlement: Instant,
}

impl RateFutureConvention {
    pub(super) fn validate(
        &self,
        parts: &RefDataParts,
        calendar: &mut impl FnMut(&CalendarSpec) -> Result<Calendar>,
    ) -> Result<()> {
        let invalid = |message: &str| Error::InvalidData(message.into());
        if !self.multiplier.is_finite() || self.multiplier <= 0.0 || self.contracts.is_empty() {
            return Err(invalid(
                "futures convention requires a positive multiplier and listed contracts",
            ));
        }
        if self.settlement_lag.amount() < 0
            || !matches!(
                self.settlement_lag.unit(),
                TenorUnit::Day | TenorUnit::BusinessDay
            )
        {
            return Err(invalid(
                "futures settlement lag must be a nonnegative day count",
            ));
        }
        lookup(
            &parts.zones,
            &self.settlement_time.zone,
            "futures settlement zone",
        )?;
        for contract in self.contracts.values() {
            let clearing = lookup(
                &parts.indices,
                &contract.settlement_index,
                "futures settlement index",
            )?;
            let IndexKind::RateFuturesSettlement {
                contract: security,
                underlying,
            } = &clearing.kind
            else {
                return Err(invalid("futures settlement index has the wrong kind"));
            };
            if clearing.observation_time.is_none() {
                return Err(invalid("futures settlement index requires a clearing time"));
            }
            let underlying = lookup(&parts.indices, underlying, "futures underlying")?;
            let ccy = match (&contract.rate, &underlying.kind) {
                (RateFutureRate::Term { .. }, IndexKind::TermRate { ccy, .. }) => *ccy,
                (
                    RateFutureRate::Compounded(period) | RateFutureRate::Averaged(period),
                    IndexKind::OvernightRate { ccy, .. },
                ) if !period.is_empty() => *ccy,
                _ => {
                    return Err(invalid(
                        "futures reference period does not match its underlying index",
                    ));
                }
            };
            if lookup(&parts.securities, security, "futures security")?.currency != ccy {
                return Err(invalid(
                    "futures underlying and settlement currencies differ",
                ));
            }
            if let Some(margin) = &self.margin_index
                && !matches!(lookup(&parts.indices, margin, "futures margin index")?.kind, IndexKind::OvernightRate { ccy: margin_ccy, .. } if margin_ccy == ccy)
            {
                return Err(invalid(
                    "futures margin index must be an overnight rate in the settlement currency",
                ));
            }
            if !calendar(&clearing.calendar)?
                .is_business_day(contract.final_clearing)?
                .0
            {
                return Err(invalid("final futures clearing must be on a business day"));
            }
        }
        Ok(())
    }
}

impl RefData {
    /// 清算指標から限月と取引所の慣行を取得する。
    ///
    /// # Errors
    ///
    /// 指標に対応する上場契約が無い場合。
    pub fn rate_future_contract(
        &self,
        index: &IndexId,
    ) -> Result<(&RateFutureConvention, &RateFutureContract)> {
        let (id, month) = lookup(&self.0.rate_futures, index, "futures settlement index")?;
        let Convention::RateFuture(convention) = self.convention(id)? else {
            return Err(Error::Internal(
                "futures index points to another kind of convention".into(),
            ));
        };
        Ok((convention, &convention.contracts[month]))
    }

    /// 指標から原資産、受渡通貨、率の日数基準を取得する。
    ///
    /// # Errors
    ///
    /// 清算指標が上場契約を参照しないか、原指標が金利でない場合。
    pub fn rate_future_underlying(
        &self,
        index: &IndexId,
    ) -> Result<(&IndexId, Currency, &DayCount)> {
        let IndexKind::RateFuturesSettlement { underlying, .. } = &self.index(index)?.kind else {
            return Err(Error::InvalidData(format!(
                "{index} is not a rate futures settlement index"
            )));
        };
        let (ccy, day_count) = match &self.index(underlying)?.kind {
            IndexKind::TermRate { ccy, day_count, .. }
            | IndexKind::OvernightRate { ccy, day_count } => (*ccy, day_count),
            _ => {
                return Err(Error::Internal(
                    "rate future has a non-rate underlying".into(),
                ));
            }
        };
        Ok((underlying, ccy, day_count))
    }

    /// 清算日から決定・公表・受渡の瞬間を求める。
    ///
    /// # Errors
    ///
    /// 暦や地域を解決できない、または時刻の順序が不正な場合。
    pub fn rate_future_clearing(
        &self,
        index: &IndexId,
        date: Date,
    ) -> Result<(RateFutureClearing, Vec<Warning>)> {
        let (convention, _) = self.rate_future_contract(index)?;
        let definition = self.index(index)?;
        let calendar = self.calendar(&definition.calendar)?;
        let point = ObservationPoint::Date(date);
        let (observed, mut warnings) = self.observation_period(index, point)?;
        let (published_at, notices) = self.publication(index, point)?;
        warnings.extend(notices);
        let (value_date, notices) = advance(
            date,
            convention.settlement_lag,
            &AdvanceRule {
                calendar: &calendar,
                bdc: BusinessDayConvention::Following,
                end_of_month: EndOfMonthRule::Ignore,
            },
        )?;
        warnings.extend(notices);
        let event = &convention.settlement_time;
        let settlement = self
            .zone(&event.zone)?
            .resolve_convention(event.on_date(value_date))?;
        if observed.start() > published_at || published_at > settlement {
            return Err(Error::InvalidData(format!(
                "clearing, publication and settlement are out of order for {index} on {date}"
            )));
        }
        Ok((
            RateFutureClearing {
                date,
                observed_at: observed.start(),
                published_at,
                settlement,
            },
            warnings,
        ))
    }

    /// 開始日から最終清算まで、取引所の営業日ごとの値洗いを列挙する。
    ///
    /// # Errors
    ///
    /// 限月・暦・時刻を解決できないか、清算日程が空または昇順でない場合。
    pub fn rate_future_schedule(
        &self,
        index: &IndexId,
        start: Date,
    ) -> Result<(Vec<RateFutureClearing>, Vec<Warning>)> {
        let (_, contract) = self.rate_future_contract(index)?;
        let calendar = self.calendar(&self.index(index)?.calendar)?;
        let mut date = start;
        let mut rows = Vec::new();
        let mut warnings = Vec::new();
        while date <= contract.final_clearing {
            let (open, notices) = calendar.is_business_day(date)?;
            warnings.extend(notices);
            if open {
                let (row, notices) = self.rate_future_clearing(index, date)?;
                warnings.extend(notices);
                rows.push(row);
            }
            if date == contract.final_clearing {
                break;
            }
            date = date.checked_add_days(1)?;
        }
        if rows.is_empty()
            || !rows
                .is_sorted_by(|a, b| a.observed_at < b.observed_at && a.settlement < b.settlement)
        {
            return Err(Error::InvalidData(format!(
                "invalid futures clearing schedule for {index} from {start}"
            )));
        }
        Ok((rows, crate::calendar::unique_warnings(warnings)))
    }
}

pub(super) fn register(
    id: &ConventionId,
    convention: &Convention,
    contracts: &mut BTreeMap<IndexId, (ConventionId, YearMonth)>,
) -> Result<()> {
    if let Convention::RateFuture(convention) = convention {
        for (&month, contract) in &convention.contracts {
            if contracts
                .insert(contract.settlement_index.clone(), (id.clone(), month))
                .is_some()
            {
                return Err(Error::InvalidData(format!(
                    "duplicate futures settlement index {}",
                    contract.settlement_index
                )));
            }
        }
    }
    Ok(())
}
`````

`````
crates/rhodium-core/src/rounding.rs
//! 確認書が定める丸めの条件。

use std::{
    cmp::Ordering,
    hash::{Hash, Hasher},
};

use scioxide_ad::{Mask, Real};
use scioxide_core::{NumericalIssue, SciOxideError};

use crate::{Error, Result};

/// 丸めの方向と、中間値の扱い。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RoundingMode {
    /// 中間値はゼロから遠い側。
    HalfAwayFromZero,
    /// 中間値は偶数倍の側。
    HalfEven,
    /// ゼロに近い側へ切り捨てる。
    TowardZero,
    /// ゼロから遠い側へ切り上げる。
    AwayFromZero,
    /// 負の無限大の方向。
    Floor,
    /// 正の無限大の方向。
    Ceiling,
}

/// 増分と丸め方の検証済みの組。
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(try_from = "RoundingInput")]
pub struct Rounding {
    increment: f64,
    mode: RoundingMode,
}

// 増分は有限かつ正なので、浮動小数点の等値とビット列の等値が一致する。
impl Eq for Rounding {}

impl Hash for Rounding {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.increment.to_bits().hash(state);
        self.mode.hash(state);
    }
}

impl Rounding {
    /// 丸めの条件を作る。
    ///
    /// # Errors
    ///
    /// 増分が有限かつ正の数でなければ [`Error::InvalidData`]。
    pub fn new(increment: f64, mode: RoundingMode) -> Result<Self> {
        if !increment.is_finite() || increment <= 0.0 {
            return Err(Error::InvalidData(
                "rounding increment must be finite and positive".into(),
            ));
        }
        Ok(Self { increment, mode })
    }

    /// 丸めの最小単位。
    #[must_use]
    pub const fn increment(self) -> f64 {
        self.increment
    }

    /// 丸め方。
    #[must_use]
    pub const fn mode(self) -> RoundingMode {
        self.mode
    }
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct RoundingInput {
    increment: f64,
    mode: RoundingMode,
}

impl TryFrom<RoundingInput> for Rounding {
    type Error = Error;

    fn try_from(input: RoundingInput) -> Result<Self> {
        Self::new(input.increment, input.mode)
    }
}

impl Rounding {
    /// 十進表記で増分の整数倍に丸め、AD では入力の微分を通す。
    ///
    /// # Errors
    ///
    /// 入出力が非有限か、整数倍の計算が 2^53 の範囲を超えれば `Numerical`。
    pub fn apply<R: Real>(&self, value: R) -> Result<R> {
        if !value.is_finite().all() {
            return Err(Error::Numerical("non-finite rounding input".into()));
        }
        let (coefficient, exponent) = decimal_parts(&format!("{:e}", self.increment))
            .ok_or_else(|| Error::Internal("invalid decimal rounding increment".into()))?;
        let rounded = value
            .try_map_scalar(|x| {
                let primal = round_primal(x, coefficient, exponent, self.mode)
                    .filter(|value| value.is_finite())
                    .ok_or_else(|| {
                        SciOxideError::numerical_instability(
                            NumericalIssue::Overflow,
                            "rounding result is not finite or the quotient exceeds 2^53",
                        )
                    })?;
                // 丸めるのは値だけ。一次微分は 1、丸めそのものの二次微分は 0 とする。
                Ok([primal, 1.0, 0.0])
            })
            .map_err(|error| Error::Numerical(error.to_string()))?;
        if !rounded.is_finite().all() {
            return Err(Error::Numerical("non-finite rounding result".into()));
        }
        Ok(rounded)
    }
}

fn decimal_parts(scientific: &str) -> Option<(u128, i32)> {
    let (digits, exponent) = scientific.split_once('e')?;
    let fractional_digits = digits.split_once('.').map_or(0, |(_, tail)| tail.len());
    let exponent = exponent
        .parse::<i32>()
        .ok()?
        .checked_sub(i32::try_from(fractional_digits).ok()?)?;
    let coefficient = digits.replace('.', "").parse().ok()?;
    Some((coefficient, exponent))
}

fn round_primal(x: f64, coefficient: u128, exponent: i32, mode: RoundingMode) -> Option<f64> {
    const MAX_TICKS: u128 = 1 << 53;
    if x == 0.0 {
        return Some(0.0);
    }

    let (digits, tail) = scaled_digits(x.abs(), exponent.checked_sub(1)?)?;
    // 1 桁細かく数えると、増分と半増分の境界がどちらも整数になる。
    let tick = coefficient.checked_mul(10)?;
    let whole = digits.checked_div(tick)?;
    let remainder = digits.checked_rem(tick)?;
    let on_tick = remainder == 0 && !tail;
    if whole > MAX_TICKS || (whole == MAX_TICKS && !on_tick) {
        return None;
    }

    let relation = match (remainder.cmp(&tick.checked_div(2)?), tail) {
        (Ordering::Equal, true) => Ordering::Greater,
        (relation, _) => relation,
    };
    let negative = x.is_sign_negative();
    let next = match mode {
        RoundingMode::TowardZero => false,
        RoundingMode::AwayFromZero => !on_tick,
        RoundingMode::Floor => negative && !on_tick,
        RoundingMode::Ceiling => !negative && !on_tick,
        RoundingMode::HalfAwayFromZero => relation != Ordering::Less,
        RoundingMode::HalfEven => {
            relation == Ordering::Greater
                || (relation == Ordering::Equal && !whole.is_multiple_of(2))
        }
    };
    let count = whole.checked_add(u128::from(next))?;
    let magnitude = count.checked_mul(coefficient)?;
    if magnitude == 0 {
        return Some(0.0);
    }
    let sign = if negative { "-" } else { "" };
    // 整数倍の十進値を作り、最後に一度だけ f64 に丸める。
    format!("{sign}{magnitude}e{exponent}").parse().ok()
}

fn scaled_digits(x: f64, place: i32) -> Option<(u128, bool)> {
    // 増分と同じく、元の f64 に戻る最短の十進表記を丸める。二進の厳密値を使うと、
    // 0.07 のように格子上にある値が格子の外に見え、Ceiling(0.07) が 0.08 になって冪等でなくなる。
    let shortest = format!("{x:e}");
    let (significand, exponent) = shortest.split_once('e')?;
    let count = exponent
        .parse::<i32>()
        .ok()?
        .checked_sub(place)?
        .checked_add(1)?;
    if count <= 0 {
        return Some((0, true));
    }
    // u128 に収まらない場合、商はすでに 2^53 を超えている。
    let count = usize::try_from(count).ok()?;
    let mut digits = significand.bytes().filter(u8::is_ascii_digit);
    let integer = (0..count).try_fold(0_u128, |value, _| {
        let digit = digits.next().unwrap_or(b'0').checked_sub(b'0')?;
        value.checked_mul(10)?.checked_add(u128::from(digit))
    })?;
    Some((integer, digits.any(|digit| digit != b'0')))
}
`````

`````
crates/rhodium-core/src/schedule.rs
//! 確認書の条件から計算期間、支払いとリセットの日程を導く。

mod grid;
mod wire;

use std::{collections::BTreeMap, num::NonZeroU8, ops::Range};

use crate::{
    Error, Result, Warning,
    calendar::{
        AdvanceRule, BaseCalendar, BusinessDayConvention, Calendar, CalendarSpec, EndOfMonthRule,
        adjust_date, advance_date, shift_calendar_date, unique_warnings,
    },
    date::{Date, DateInterval, Weekday},
    day_count::{AccrualContext, DayCount, MonthEndSchedule, PeriodFrame},
    id::CalendarId,
    refdata::CalendarCache,
    tenor::{DateOffset, Frequency, TenorLength, TenorUnit},
    time::{Instant, LocalEventTime, TzId, ZoneRules},
};

/// 一つの日程に含められる計算期間の上限。
pub const MAX_PERIODS: usize = 10_000;

/// 確認書が定める日程の条件。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ScheduleSpec {
    /// 日付の規則から導く。
    Generated {
        /// 計算期間を生成する条件。
        periods: GeneratedPeriods,
        /// 生成した期間を支払日に対応させる規則。
        payment: PaymentRule,
    },
    /// 確認書で計算期間と支払日を列挙する。
    Explicit(ExplicitPeriods),
}

/// 未調整の開始日を決める条件。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StartRule {
    /// 確認書の開始日。
    Date(Date),
    /// 約定日を計算期間のカレンダーでずらした日。
    FromTradeDate(DateOffset),
}

/// 未調整の終了日を決める条件。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndRule {
    /// 確認書の終了日。
    Date(Date),
    /// 未調整の開始日からのテナー。
    AfterStart(TenorLength),
}

/// スタブの位置と格子を作る向き。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum StubSpec {
    /// 開始日から進み、終了日と一致する必要がある。
    NoStub,
    /// 終了日から遡り、初回に余りを置く。
    Initial(StubLength),
    /// 開始日から進み、最後に余りを置く。
    Final(StubLength),
    /// 最初の通常期間の開始日を指定する。
    FirstRegular(Date),
    /// 最後の通常期間の終了日を指定する。
    LastRegular(Date),
    /// 両端の通常期間の境界を指定する。
    Both {
        /// 最初の通常期間の開始日。
        first_regular: Date,
        /// 最後の通常期間の終了日。
        last_regular: Date,
    },
}

/// 格子の余りの扱い。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum StubLength {
    /// 余りを独立した期間にする。
    Short,
    /// 余りを隣の通常期間と合わせる。
    Long,
}

/// 利息を数える区間の境界。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum AccrualDates {
    /// 営業日調整後の日付。
    Adjusted,
    /// 営業日調整前の日付。
    Unadjusted,
}

/// 1 日から 31 日までのロール日。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(try_from = "u8")]
pub struct DayOfMonth(u8);

impl DayOfMonth {
    /// 月内のロール日を作る。
    ///
    /// # Errors
    ///
    /// 1〜31 の範囲外なら [`Error::InvalidData`]。
    pub fn new(day: u8) -> Result<Self> {
        if (1..=31).contains(&day) {
            Ok(Self(day))
        } else {
            Err(invalid("roll day must be in 1..=31"))
        }
    }
}

impl TryFrom<u8> for DayOfMonth {
    type Error = Error;

    fn try_from(day: u8) -> Result<Self> {
        Self::new(day)
    }
}

/// 月内で曜日が現れる順番。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum WeekOfMonth {
    /// 第 1。
    First,
    /// 第 2。
    Second,
    /// 第 3。
    Third,
    /// 第 4。
    Fourth,
    /// 最後。
    Last,
}

/// 未調整の日付を月や週のどこに置くかを決める規則。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum RollConvention {
    /// 指定の日。短い月では月末。
    OnDay(DayOfMonth),
    /// 月の最終暦日。
    EndOfMonth,
    /// 月内で指定した順番の曜日。
    NthWeekday {
        /// 月内の順番。
        week: WeekOfMonth,
        /// 曜日。
        #[serde(with = "wire::weekday")]
        weekday: Weekday,
    },
    /// 週の指定の曜日。
    OnWeekday(#[serde(with = "wire::weekday")] Weekday),
    /// 直前の調整後の日付から進める FRN の慣行。
    Frn,
}

/// 検証済みの計算期間の生成条件。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneratedPeriods {
    start: StartRule,
    end: EndRule,
    step: Frequency,
    roll: Option<RollConvention>,
    stub: StubSpec,
    calendar: CalendarSpec,
    bdc: BusinessDayConvention,
    accrual: AccrualDates,
}

impl GeneratedPeriods {
    /// 期間の刻み、ロール、スタブを組み合わせる。ロールの省略時は格子の起点の日を使う。
    ///
    /// # Errors
    ///
    /// 規則の組が許されないか、明示された開始日が終了日より前でなければ
    /// [`Error::InvalidData`]。
    #[expect(
        clippy::too_many_arguments,
        reason = "確認書が独立に定める条件を一度に検証する"
    )]
    pub fn new(
        start: StartRule,
        end: EndRule,
        step: Frequency,
        roll: Option<RollConvention>,
        stub: StubSpec,
        calendar: CalendarSpec,
        bdc: BusinessDayConvention,
        accrual: AccrualDates,
    ) -> Result<Self> {
        let unit = step.step().map(TenorLength::unit);
        let monthly = matches!(unit, Some(TenorUnit::Month | TenorUnit::Year));
        // FRN は直前の調整後の日付から ModifiedFollowing で月を進める慣行なので、月か年の刻みで
        // 開始日から先へ生成するときだけ意味を持つ（終了日から遡る Initial と LastRegular は不可）。
        // 曜日のロールは週の刻み、月末と第 n 曜日のロールは月か年の刻み（一回限りならロールは使わない）に限る。
        let valid_roll = match roll {
            Some(RollConvention::Frn) => {
                monthly
                    && bdc == BusinessDayConvention::ModifiedFollowing
                    && !matches!(stub, StubSpec::Initial(_) | StubSpec::LastRegular(_))
            }
            Some(RollConvention::OnWeekday(_)) => unit == Some(TenorUnit::Week),
            Some(RollConvention::EndOfMonth | RollConvention::NthWeekday { .. }) => {
                unit.is_none() || monthly
            }
            Some(RollConvention::OnDay(_)) | None => true,
        };
        // 通常期間の境界を指定するスタブは、刻みの無い一回限りの日程では格子を作れない。
        if !valid_roll
            || (unit.is_none()
                && matches!(
                    stub,
                    StubSpec::FirstRegular(_) | StubSpec::LastRegular(_) | StubSpec::Both { .. }
                ))
        {
            return Err(invalid("incompatible schedule frequency, roll or stub"));
        }
        if let (StartRule::Date(start), EndRule::Date(end)) = (start, end)
            && start >= end
        {
            return Err(invalid("schedule start must precede its end"));
        }
        Ok(Self {
            start,
            end,
            step,
            roll,
            stub,
            calendar,
            bdc,
            accrual,
        })
    }
}

/// 確認書で指定した一つの期間。日付は調整後の日付で書く。
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExplicitPeriod {
    /// 利息を数える区間。
    pub accrual: DateInterval,
    /// 確認書の支払日。
    pub payment: Date,
    /// 確認書のリセット日。
    pub reset: Option<Date>,
    /// 決済の現地時刻と地域。
    pub settlement: Option<LocalEventTime>,
    /// 利息を数える区間を覆う、連続した通常の参照期間の列（ICMA の日数計算が使う）。
    ///
    /// 二つ以上の通常期間にまたがる長いスタブは、一つの期間に参照期間を並べて書く。
    /// 空なら、`accrual` が暦のちょうど一刻みのときだけ `accrual` 自身を参照期間にする。
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub references: Vec<DateInterval>,
}

/// 確認書の期間と支払日を検証済みの日程条件として持つ。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExplicitPeriods {
    periods: Vec<ExplicitPeriod>,
    on_holiday: ExplicitDatePolicy,
    combine: PeriodCombination,
    frequency: Option<Frequency>,
    /// 期間ごとに検証済みの参照枠。参照期間も頻度も無い期間は `None`。
    frames: Vec<Option<PeriodFrame>>,
}

impl ExplicitPeriods {
    /// 連続する明示期間と頻度を登録する。
    ///
    /// 参照期間は ICMA の日数計算が使う通常期間の連続した列で、利息を数える区間を覆い、
    /// それぞれが [`PeriodFrame::new`] と同じ規則で頻度のおよそ一刻み（ロールや営業日調整による
    /// ずれは許す）でなければならない。二つ以上の通常期間にまたがる長いスタブは、一つの期間に
    /// 参照期間を並べて書く。二つの期間に分けると、RFR の複利が分割点で切れる別のクーポンになる。
    ///
    /// # Errors
    ///
    /// 期間が空、重なるか連続しない、支払日が期首より前か減少する、参照期間が
    /// 連続しないか利息を数える区間を覆わない、頻度が無いのに参照期間がある、または
    /// 参照期間が頻度の一刻みでなければ [`Error::InvalidData`]。
    /// 期間数が上限を超える場合は [`Error::InvalidInput`]。
    pub fn new(
        periods: Vec<ExplicitPeriod>,
        on_holiday: ExplicitDatePolicy,
        combine: PeriodCombination,
        frequency: Option<Frequency>,
    ) -> Result<Self> {
        check_count(periods.len())?;
        if periods.is_empty()
            || periods
                .iter()
                .any(|p| p.accrual.is_empty() || p.payment < p.accrual.start())
            || periods
                .windows(2)
                .any(|p| p[0].accrual.end() != p[1].accrual.start() || p[0].payment > p[1].payment)
        {
            return Err(invalid("invalid explicit accrual or payment dates"));
        }
        let frames = periods
            .iter()
            .map(|p| explicit_frame(p, frequency))
            .collect::<Result<_>>()?;
        Ok(Self {
            periods,
            on_holiday,
            combine,
            frequency,
            frames,
        })
    }

    /// 確認書の順に並ぶ期間。
    #[must_use]
    pub fn periods(&self) -> &[ExplicitPeriod] {
        &self.periods
    }

    /// 通常の利払いの頻度。
    #[must_use]
    pub const fn frequency(&self) -> Option<Frequency> {
        self.frequency
    }
}

/// 確認書に列挙された支払日が休日の場合の扱い。
#[derive(Debug, Clone, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum ExplicitDatePolicy {
    /// 確認書の日付を使い、休日と判明した場合は警告する。
    #[default]
    KeepAsConfirmed,
    /// 確認書が定めた暦と規則で支払日を調整し直す。
    ReadjustWithRule {
        /// 営業日調整規則。
        bdc: BusinessDayConvention,
        /// 支払いのカレンダー。
        calendar: CalendarSpec,
    },
}

/// 支払いまたはリセットの起点。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum PeriodAnchor {
    /// 調整後の期首。
    PeriodStart,
    /// 調整後の期末。
    PeriodEnd,
}

impl PeriodAnchor {
    const fn date(self, interval: DateInterval) -> Date {
        match self {
            Self::PeriodStart => interval.start(),
            Self::PeriodEnd => interval.end(),
        }
    }
}

/// 計算期間の額を一つの支払いに合わせる規則。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum PeriodCombination {
    /// 期間の額の和。
    #[default]
    Sum,
    /// それまでの額をスプレッド込みの率で複利にする。
    StraightCompounding,
    /// それまでの額をスプレッドを除いた率で複利にする。
    FlatCompounding,
    /// 率の部分だけを複利、スプレッドを単利で数える。
    SpreadExclusiveCompounding,
}

/// 計算期間を支払いに割り当てる規則。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum PaymentGrouping {
    /// 計算期間ごとに支払う。
    #[default]
    EachPeriod,
    /// 所定の期間数ごとに支払う。最後には余りをまとめる。
    EveryN {
        /// 一回の支払いに含める期間数。
        periods: NonZeroU8,
        /// 期間の支払額の合わせ方。
        combine: PeriodCombination,
    },
}

/// 計算期間から支払日と決済時刻を決める規則。
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PaymentRule {
    /// 日付の起点。
    pub anchor: PeriodAnchor,
    /// 起点から支払いまでのオフセット。
    pub offset: DateOffset,
    /// 支払いのカレンダー。
    pub calendar: CalendarSpec,
    /// 支払いの営業日調整規則。
    pub bdc: BusinessDayConvention,
    /// 期間をまとめる規則。
    #[serde(default)]
    pub grouping: PaymentGrouping,
    /// 決済時刻の上書き。
    pub settlement: Option<LocalEventTime>,
}

/// リセット日を決める規則。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResetRule {
    /// 日付の起点。
    pub anchor: PeriodAnchor,
    /// 起点から fixing までのオフセット。
    pub offset: DateOffset,
    /// fixing のカレンダー。
    pub calendar: CalendarSpec,
    /// fixing の営業日調整規則。
    pub bdc: BusinessDayConvention,
}

/// 支払日を UTC の瞬間に解決するための参照データ。
#[derive(Debug, Clone, Copy)]
pub struct SettlementContext<'a> {
    /// 利用できる地域とそのゾーン規則。
    pub zones: &'a BTreeMap<TzId, ZoneRules>,
    /// 支払い通貨の決済時刻。明示された時刻が優先される。
    pub default_time: Option<&'a LocalEventTime>,
    /// 明示の支払日の休日判定に使う決済暦の式。確認書の日付を保持する場合に必要。
    pub calendar: Option<&'a CalendarSpec>,
}

impl SettlementContext<'_> {
    fn resolve(&self, date: Date, specified: Option<&LocalEventTime>) -> Result<Instant> {
        let event = specified
            .or(self.default_time)
            .ok_or_else(|| invalid("payment requires a settlement time"))?;
        let zone = self
            .zones
            .get(&event.zone)
            .ok_or_else(|| invalid(format!("unknown settlement zone {}", event.zone)))?;
        zone.resolve_convention(event.on_date(date))
    }
}

/// 計算期間の種類。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeriodKind {
    /// 通常の期間。
    Regular,
    /// 短いスタブ。
    ShortStub,
    /// 長いスタブ。
    LongStub,
    /// 確認書で明示された期間。
    Explicit,
}

/// 日程を構成する一つの計算期間。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Period {
    kind: PeriodKind,
    unadjusted: DateInterval,
    adjusted: DateInterval,
    accrual: DateInterval,
    frame: Option<PeriodFrame>,
    reset: Option<Date>,
}

impl Period {
    /// 期間の種類。
    #[must_use]
    pub const fn kind(&self) -> PeriodKind {
        self.kind
    }

    /// 営業日調整前の区間。
    #[must_use]
    pub const fn unadjusted(&self) -> DateInterval {
        self.unadjusted
    }

    /// 営業日調整後の区間。
    #[must_use]
    pub const fn adjusted(&self) -> DateInterval {
        self.adjusted
    }

    /// 利息を数える区間。
    #[must_use]
    pub const fn accrual(&self) -> DateInterval {
        self.accrual
    }

    /// 枠の要る日数計算のための参照期間と頻度。
    #[must_use]
    pub const fn frame(&self) -> Option<&PeriodFrame> {
        self.frame.as_ref()
    }

    /// 確認書が明示したリセット日。
    #[must_use]
    pub const fn reset(&self) -> Option<Date> {
        self.reset
    }
}

/// 連続する一つ以上の計算期間の支払い。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScheduledPayment {
    date: Date,
    at: Instant,
    periods: Range<usize>,
}

impl ScheduledPayment {
    /// 現地の支払日。
    #[must_use]
    pub const fn date(&self) -> Date {
        self.date
    }

    /// 決済の UTC 時刻。
    #[must_use]
    pub const fn at(&self) -> Instant {
        self.at
    }

    /// 対象となる計算期間の 0 始まりの半開範囲。
    #[must_use]
    pub fn periods(&self) -> Range<usize> {
        self.periods.clone()
    }
}

/// 検証済みの計算期間と支払日程。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Schedule {
    periods: Vec<Period>,
    payments: Vec<ScheduledPayment>,
    combine: PeriodCombination,
    month_end: MonthEndSchedule,
}

impl Schedule {
    /// 計算期間の順序付きの列。
    #[must_use]
    pub fn periods(&self) -> &[Period] {
        &self.periods
    }

    /// 計算期間の順に並ぶ支払い。
    #[must_use]
    pub fn payments(&self) -> &[ScheduledPayment] {
        &self.payments
    }

    /// 複数の計算期間の額をまとめる規則。
    #[must_use]
    pub const fn combine(&self) -> PeriodCombination {
        self.combine
    }

    /// 商品を構築するときに、日数基準を適用するための参照枠を確認する。
    ///
    /// # Errors
    ///
    /// ACT/ACT ICMA または ACT/365L に必要な参照枠か年あたりの頻度が
    /// 一つでも無ければ [`Error::InvalidData`]。
    pub fn validate_day_count(&self, basis: &DayCount) -> Result<()> {
        if matches!(basis, DayCount::ActActIcma | DayCount::Act365L)
            && self.periods.iter().any(|p| {
                p.frame()
                    .and_then(|f| f.frequency().periods_per_year())
                    .is_none()
            })
        {
            return Err(invalid(format!(
                "{basis} requires a reference period and frequency for every accrual"
            )));
        }
        Ok(())
    }

    /// 計算期間の参照枠、月末規則、満期日、および BUS/252 の暦を組み合わせる。
    /// 期間の添字が範囲外なら None。
    #[must_use]
    pub fn accrual_context<'a>(
        &'a self,
        index: usize,
        calendar: Option<&'a Calendar>,
    ) -> Option<AccrualContext<'a>> {
        Some(AccrualContext {
            frame: self.periods.get(index)?.frame(),
            month_end: self.month_end,
            maturity: self.periods.last().map(|p| p.accrual.end()),
            calendar,
        })
    }
}

/// 契約の計算期間と UTC の受渡時刻を導く。
///
/// # Errors
///
/// 期間・支払いの条件や参照先が無効、暦やゾーンの範囲が不足、または計算が値域を超えた場合。
pub fn generate(
    spec: &ScheduleSpec,
    trade_date: Date,
    calendars: &BTreeMap<CalendarId, BaseCalendar>,
    settlement: &SettlementContext<'_>,
) -> Result<(Schedule, Vec<Warning>)> {
    generate_using(spec, trade_date, settlement, |spec| {
        Calendar::compile(spec, calendars)
    })
}

/// コンパイル済みの暦を再利用しながら日程を生成する。
///
/// # Errors
///
/// [`generate`] と同じ失敗。
pub fn generate_with_calendars(
    spec: &ScheduleSpec,
    trade_date: Date,
    settlement: &SettlementContext<'_>,
    calendars: &mut CalendarCache<'_>,
) -> Result<(Schedule, Vec<Warning>)> {
    generate_using(spec, trade_date, settlement, |spec| calendars.resolve(spec))
}

fn generate_using(
    spec: &ScheduleSpec,
    trade_date: Date,
    settlement: &SettlementContext<'_>,
    mut resolve: impl FnMut(&CalendarSpec) -> Result<Calendar>,
) -> Result<(Schedule, Vec<Warning>)> {
    let mut warnings = Vec::new();
    let schedule = match spec {
        ScheduleSpec::Generated { periods, payment } => generate_regular(
            periods,
            payment,
            trade_date,
            &mut resolve,
            settlement,
            &mut warnings,
        )?,
        ScheduleSpec::Explicit(periods) => {
            generate_explicit(periods, &mut resolve, settlement, &mut warnings)?
        }
    };
    // 支払いのオフセットや明示の支払日の再調整の後でも、支払いは計算期間の順に並ぶことを保証する。
    if schedule.payments.windows(2).any(|p| p[0].date > p[1].date) {
        return Err(invalid("payment dates must be non-decreasing"));
    }
    Ok((schedule, unique_warnings(warnings)))
}

/// 各計算期間のリセット日を求める。確認書のリセット日があれば優先する。
///
/// # Errors
///
/// カレンダーを解決できないか、リセット日が参照データの範囲または日付の値域を超えた場合。
pub fn reset_dates(
    schedule: &Schedule,
    rule: &ResetRule,
    calendars: &BTreeMap<CalendarId, BaseCalendar>,
) -> Result<(Vec<Date>, Vec<Warning>)> {
    let calendar = schedule
        .periods
        .iter()
        .any(|p| p.reset.is_none())
        .then(|| Calendar::compile(&rule.calendar, calendars))
        .transpose()?;
    let mut warnings = Vec::new();
    let dates = schedule
        .periods
        .iter()
        .map(|p| match (p.reset, calendar.as_ref()) {
            (Some(date), _) => Ok(date),
            (None, Some(calendar)) => offset_date(
                rule.anchor.date(p.adjusted),
                rule.offset,
                rule.bdc,
                calendar,
                &mut warnings,
            ),
            (None, None) => Err(Error::Internal("reset calendar was not compiled".into())),
        })
        .collect::<Result<_>>()?;
    Ok((dates, unique_warnings(warnings)))
}

fn generate_regular(
    spec: &GeneratedPeriods,
    payment: &PaymentRule,
    trade_date: Date,
    resolve: &mut impl FnMut(&CalendarSpec) -> Result<Calendar>,
    settlement: &SettlementContext<'_>,
    warnings: &mut Vec<Warning>,
) -> Result<Schedule> {
    let calendar = resolve(&spec.calendar)?;
    let periods = grid::generate(spec, trade_date, &calendar, warnings)?;
    let payment_calendar = if spec.calendar == payment.calendar {
        calendar
    } else {
        resolve(&payment.calendar)?
    };
    // 期間を先頭から所定の数ずつ組にし、最後の組は余りでよい。一期間ずつなら合わせ方は和と同じ。
    let (size, combine) = match payment.grouping {
        PaymentGrouping::EachPeriod => (1, PeriodCombination::Sum),
        PaymentGrouping::EveryN { periods, combine } => (usize::from(periods.get()), combine),
    };
    let payments = periods
        .chunks(size)
        .enumerate()
        .map(|(group, periods)| {
            let anchor = match payment.anchor {
                PeriodAnchor::PeriodStart => periods[0].adjusted.start(),
                PeriodAnchor::PeriodEnd => periods[periods.len() - 1].adjusted.end(),
            };
            let date = offset_date(
                anchor,
                payment.offset,
                payment.bdc,
                &payment_calendar,
                warnings,
            )?;
            let at = settlement.resolve(date, payment.settlement.as_ref())?;
            let start = group * size;
            Ok(ScheduledPayment {
                date,
                at,
                periods: start..start + periods.len(),
            })
        })
        .collect::<Result<_>>()?;
    // 30/360 US の月末の慣行は、月末のロールを明示した日程だけに適用する（旧実装と同じ）。
    let month_end = if spec.roll == Some(RollConvention::EndOfMonth) {
        MonthEndSchedule::Aligned
    } else {
        MonthEndSchedule::NotAligned
    };
    Ok(Schedule {
        periods,
        payments,
        combine,
        month_end,
    })
}

fn generate_explicit(
    spec: &ExplicitPeriods,
    resolve: &mut impl FnMut(&CalendarSpec) -> Result<Calendar>,
    settlement: &SettlementContext<'_>,
    warnings: &mut Vec<Warning>,
) -> Result<Schedule> {
    let calendar_spec = match &spec.on_holiday {
        ExplicitDatePolicy::KeepAsConfirmed => settlement
            .calendar
            .ok_or_else(|| invalid("confirmed payment dates require a settlement calendar"))?,
        ExplicitDatePolicy::ReadjustWithRule { calendar, .. } => calendar,
    };
    let calendar = resolve(calendar_spec)?;
    let mut payments: Vec<ScheduledPayment> = Vec::new();
    let periods = spec
        .periods
        .iter()
        .enumerate()
        .map(|(index, p)| {
            let date = match spec.on_holiday {
                ExplicitDatePolicy::KeepAsConfirmed => {
                    match calendar.lookup(p.payment, warnings) {
                        Ok(false) => {
                            warnings.push(Warning::ExplicitDateIsHoliday { date: p.payment });
                        }
                        // 確認書の日付は暦の適用範囲外でも確定している。休日かどうかだけが不明。
                        Ok(true) | Err(Error::MissingData(_)) => {}
                        Err(error) => return Err(error),
                    }
                    p.payment
                }
                ExplicitDatePolicy::ReadjustWithRule { bdc, .. } => {
                    adjust_date(p.payment, bdc, &calendar, warnings)?
                }
            };
            if date < p.accrual.start() {
                return Err(invalid("adjusted payment precedes its accrual period"));
            }
            let at = settlement.resolve(date, p.settlement.as_ref())?;
            // 支払日と決済の瞬間がともに一致する連続した期間だけを一つの支払いにまとめる。
            // 同じ支払日でも確認書の時刻が違えば、別々の瞬間に払う別の支払いになる。
            if let Some(last) = payments.last_mut()
                && last.date == date
                && last.at == at
            {
                last.periods.end = index + 1;
            } else {
                payments.push(ScheduledPayment {
                    date,
                    at,
                    periods: index..index + 1,
                });
            }
            Ok(Period {
                kind: PeriodKind::Explicit,
                unadjusted: p.accrual,
                adjusted: p.accrual,
                accrual: p.accrual,
                frame: spec.frames[index].clone(),
                reset: p.reset,
            })
        })
        .collect::<Result<_>>()?;
    Ok(Schedule {
        periods,
        payments,
        combine: spec.combine,
        // 明示の期間はロールを持たないので、月末の慣行に従うかを日程から判定できない。
        month_end: MonthEndSchedule::NotAligned,
    })
}

/// 明示の期間の参照枠を作る。参照期間は、生成した日程と同じ規則（[`PeriodFrame::new`]）で検証する。
fn explicit_frame(
    period: &ExplicitPeriod,
    frequency: Option<Frequency>,
) -> Result<Option<PeriodFrame>> {
    let Some(frequency) = frequency else {
        return if period.references.is_empty() {
            Ok(None)
        } else {
            Err(invalid(format!(
                "explicit reference periods {:?} require a frequency",
                period.references
            )))
        };
    };
    let references = if period.references.is_empty() {
        // 暦上で正確に一刻みの区間だけ、自身を通常の参照期間として使える。
        // スタブや営業日調整でずれた期間は参照期間を一意に復元できないので、枠を作らない。
        let regular = frequency.step().is_some_and(|step| {
            shift_calendar_date(
                period.accrual.start(),
                i32::from(step.amount()),
                step.unit(),
                true,
            )
            .is_ok_and(|end| end == period.accrual.end())
        });
        if !regular {
            return Ok(None);
        }
        vec![period.accrual]
    } else {
        period.references.clone()
    };
    PeriodFrame::new(period.accrual, references, frequency)
        .map(Some)
        .map_err(|error| {
            invalid(format!(
                "explicit reference periods of {:?} are not regular periods: {error}",
                period.accrual
            ))
        })
}

fn offset_date(
    date: Date,
    offset: DateOffset,
    bdc: BusinessDayConvention,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<Date> {
    advance_date(
        date,
        i32::from(offset.amount()),
        offset.unit(),
        &AdvanceRule {
            end_of_month: EndOfMonthRule::Ignore,
            bdc,
            calendar,
        },
        warnings,
    )
}

fn check_count(count: usize) -> Result<()> {
    if count > MAX_PERIODS {
        Err(Error::InvalidInput(format!(
            "schedule exceeds {MAX_PERIODS} periods"
        )))
    } else {
        Ok(())
    }
}

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidData(message.into())
}
`````

`````
crates/rhodium-core/src/schedule/grid.rs
//! 取引の境界と通常の格子から、実際の期間とスタブの参照期間を作る。

use chrono::NaiveDate;

use crate::{
    Result, Warning,
    calendar::{
        AdvanceRule, BusinessDayConvention, Calendar, EndOfMonthRule, adjust_date, advance_date,
        is_last_business_day, shift_business, shift_calendar_date,
    },
    date::{Date, DateInterval, date_error},
    day_count::PeriodFrame,
    tenor::{TenorLength, TenorUnit},
};

use super::{
    AccrualDates, DayOfMonth, EndRule, GeneratedPeriods, Period, PeriodKind, RollConvention,
    StartRule, StubLength, StubSpec, WeekOfMonth, check_count, invalid, offset_date,
};

pub(super) fn generate(
    spec: &GeneratedPeriods,
    trade_date: Date,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<Vec<Period>> {
    let start = match spec.start {
        StartRule::Date(date) => date,
        StartRule::FromTradeDate(offset) => offset_date(
            trade_date,
            offset,
            BusinessDayConvention::Unadjusted,
            calendar,
            warnings,
        )?,
    };
    let end = match spec.end {
        EndRule::Date(date) => date,
        EndRule::AfterStart(length) => end_after_start(spec, start, length, calendar, warnings)?,
    };
    if start >= end {
        return Err(invalid("schedule start must precede its end"));
    }
    let periods = match spec.step.step() {
        Some(step) => build_grid(spec, start, end, step, calendar, warnings)?,
        None => vec![DraftPeriod::between(start.into(), end.into())],
    };
    periods
        .into_iter()
        .map(|p| p.resolve(spec, calendar, warnings))
        .collect()
}

fn end_after_start(
    spec: &GeneratedPeriods,
    start: Date,
    length: TenorLength,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<Date> {
    let rule = AdvanceRule {
        end_of_month: EndOfMonthRule::Ignore,
        bdc: BusinessDayConvention::Unadjusted,
        calendar,
    };
    let end = advance_date(
        start,
        i32::from(length.amount()),
        length.unit(),
        &rule,
        warnings,
    )?;
    let monthly = |unit| matches!(unit, TenorUnit::Month | TenorUnit::Year);
    // 起点が月の格子に乗っていれば、テナーで求めた終了日も同じ格子に置く。
    if monthly(length.unit())
        && spec.step.step().is_some_and(|step| monthly(step.unit()))
        && let Some(
            roll @ (RollConvention::OnDay(_)
            | RollConvention::EndOfMonth
            | RollConvention::NthWeekday { .. }),
        ) = spec.roll
        && apply_roll(start, roll, length.unit())? == start
    {
        return apply_roll(end, roll, length.unit());
    }
    Ok(end)
}

#[derive(Clone, Copy)]
struct GridDate {
    raw: Date,
    // FRN が次の刻みにも使う調整後の日付。
    adjusted: Option<Date>,
}

impl From<Date> for GridDate {
    fn from(raw: Date) -> Self {
        Self {
            raw,
            adjusted: None,
        }
    }
}

impl GridDate {
    fn adjust(
        self,
        bdc: BusinessDayConvention,
        calendar: &Calendar,
        warnings: &mut Vec<Warning>,
    ) -> Result<Date> {
        self.adjusted
            .map_or_else(|| adjust_date(self.raw, bdc, calendar, warnings), Ok)
    }
}

struct DraftPeriod {
    start: GridDate,
    end: GridDate,
    kind: PeriodKind,
    references: Vec<(GridDate, GridDate)>,
}

impl DraftPeriod {
    fn between(a: GridDate, b: GridDate) -> Self {
        let (start, end) = ordered(a, b);
        Self {
            start,
            end,
            kind: PeriodKind::Regular,
            references: Vec::new(),
        }
    }

    fn resolve(
        self,
        spec: &GeneratedPeriods,
        calendar: &Calendar,
        warnings: &mut Vec<Warning>,
    ) -> Result<Period> {
        let unadjusted = nonempty_interval(self.start.raw, self.end.raw)?;
        let adjusted = nonempty_interval(
            self.start.adjust(spec.bdc, calendar, warnings)?,
            self.end.adjust(spec.bdc, calendar, warnings)?,
        )?;
        let accrual = match spec.accrual {
            AccrualDates::Adjusted => adjusted,
            AccrualDates::Unadjusted => unadjusted,
        };
        // 設計の手順 6: 通常の期間はそれ自身が参照期間。スタブの仮想の通常期間は、
        // 利息を数える区間と同じ側（調整後か調整前）の日付で表し、枠が区間を覆うようにする。
        let frame = if self.references.is_empty() {
            // それ自身を参照期間にできるのは、通常期間の規則（頻度のおよそ一刻み）を満たすときだけ。
            // FRN の慣行で日付がずれ、確認書の満期で閉じた最後の期間は一刻みに満たないことがあり、
            // その期間は参照枠を持たない（枠の要る日数計算は構築のときに拒否し、他の日数基準は使える）。
            // 覆うことと連続することはそれ自身なら必ず満たすので、失敗するのはこの規則だけ。
            PeriodFrame::new(accrual, vec![accrual], spec.step).ok()
        } else {
            let regular = self
                .references
                .into_iter()
                .map(|(start, end)| {
                    let (start, end) = match spec.accrual {
                        AccrualDates::Adjusted => (
                            start.adjust(spec.bdc, calendar, warnings)?,
                            end.adjust(spec.bdc, calendar, warnings)?,
                        ),
                        AccrualDates::Unadjusted => (start.raw, end.raw),
                    };
                    nonempty_interval(start, end)
                })
                .collect::<Result<_>>()?;
            Some(
                PeriodFrame::new(accrual, regular, spec.step).map_err(|error| {
                    invalid(format!("invalid schedule reference periods: {error}"))
                })?,
            )
        };
        Ok(Period {
            kind: self.kind,
            unadjusted,
            adjusted,
            accrual,
            frame,
            reset: None,
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Direction {
    Forward,
    Backward,
}

impl Direction {
    const fn sign(self) -> i32 {
        match self {
            Self::Forward => 1,
            Self::Backward => -1,
        }
    }

    fn before(self, a: Date, b: Date) -> bool {
        match self {
            Self::Forward => a < b,
            Self::Backward => a > b,
        }
    }
}

struct Grid {
    anchor: Date,
    current: GridDate,
    step: TenorLength,
    roll: RollConvention,
    direction: Direction,
    index: i32,
}

impl Grid {
    fn new(
        anchor: GridDate,
        step: TenorLength,
        roll: Option<RollConvention>,
        direction: Direction,
    ) -> Self {
        let roll = roll.unwrap_or_else(|| RollConvention::OnDay(DayOfMonth(anchor.raw.ymd().2)));
        Self {
            anchor: anchor.raw,
            current: anchor,
            step,
            roll,
            direction,
            index: 0,
        }
    }

    fn next(&mut self, calendar: &Calendar, warnings: &mut Vec<Warning>) -> Result<GridDate> {
        self.index += self.direction.sign();
        let size = i32::from(self.step.amount());
        let unit = self.step.unit();
        let next = match (self.roll, unit) {
            (RollConvention::Frn, _) => frn_next(self.current, size, unit, calendar, warnings)?,
            // 同じ向きの営業日加算は直前の格子点から進め、済んだ暦を走査し直さない。
            (_, TenorUnit::BusinessDay) => shift_business(
                self.current.raw,
                size * self.direction.sign(),
                calendar,
                warnings,
            )?
            .into(),
            (roll, _) => apply_roll(self.unrolled()?, roll, unit)?.into(),
        };
        if !self.direction.before(self.current.raw, next.raw) {
            return Err(invalid("schedule grid did not advance"));
        }
        self.current = next;
        Ok(next)
    }

    fn unrolled(&self) -> Result<Date> {
        // index は期間数の上限、amount は u16 で制限される。
        shift_calendar_date(
            self.anchor,
            i32::from(self.step.amount()) * self.index,
            self.step.unit(),
            false,
        )
    }

    fn reaches_frn_end(&self, end: Date, next: GridDate) -> Result<bool> {
        if self.roll != RollConvention::Frn || self.direction != Direction::Forward {
            return Ok(false);
        }
        // 営業日調整で格子が後ろにずれても、ISDA の最終の計算期間は確認書の終了日で閉じる。
        Ok(next.adjusted == Some(end) || (end < next.raw && self.unrolled()? <= end))
    }
}

fn build_grid(
    spec: &GeneratedPeriods,
    start: Date,
    end: Date,
    step: TenorLength,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<Vec<DraftPeriod>> {
    // 設計の手順 2: スタブの指定から、通常の格子が覆う区間 [first, last)、生成の向き、
    // 向きの先に残る余りの扱いを決める。境界を指定するスタブの余りは、格子の外で後から作る。
    let (first, last, direction, residual) = match spec.stub {
        StubSpec::NoStub => (start, end, Direction::Forward, None),
        StubSpec::Initial(length) => (start, end, Direction::Backward, Some(length)),
        StubSpec::Final(length) => (start, end, Direction::Forward, Some(length)),
        StubSpec::FirstRegular(date) => (date, end, Direction::Forward, None),
        StubSpec::LastRegular(date) => (start, date, Direction::Backward, None),
        StubSpec::Both {
            first_regular,
            last_regular,
        } => (first_regular, last_regular, Direction::Forward, None),
    };
    if first < start || last > end || first >= last {
        return Err(invalid(
            "regular dates must form a nonempty interval inside the schedule",
        ));
    }
    let (anchor, limit) = match direction {
        Direction::Forward => (first, last),
        Direction::Backward => (last, first),
    };
    let mut grid = Grid::new(anchor.into(), step, spec.roll, direction);
    if apply_roll(anchor, grid.roll, step.unit())? != anchor {
        return Err(invalid(
            "regular grid anchor does not follow its roll convention",
        ));
    }
    // 設計の手順 3〜4: 格子点を一つずつ作り、向きの先の端に届いた時点で止める。
    // ちょうど端に着けば通常の期間、FRN は確認書の端で閉じ、越えたら余りをスタブにする。
    let mut periods = Vec::new();
    loop {
        let current = grid.current;
        let next = grid.next(calendar, warnings)?;
        if next.raw == limit {
            check_count(periods.len() + 1)?;
            periods.push(DraftPeriod::between(current, next));
            break;
        }
        if grid.reaches_frn_end(limit, next)? {
            check_count(periods.len() + 1)?;
            periods.push(DraftPeriod::between(current, limit.into()));
            break;
        }
        if direction.before(limit, next.raw) {
            let length =
                residual.ok_or_else(|| invalid("schedule does not end on its regular grid"))?;
            add_residual(&mut periods, current, next, limit.into(), length, direction)?;
            break;
        }
        check_count(periods.len() + 1)?;
        periods.push(DraftPeriod::between(current, next));
    }
    if direction == Direction::Backward {
        periods.reverse();
    }
    // FirstRegular、LastRegular、Both は、指定した通常期間の境界から外側へ格子を延ばしてスタブを作る。
    if start < first {
        check_count(periods.len() + 1)?;
        let stub = specified_stub(
            Some(grid.roll),
            step,
            periods[0].start,
            start,
            Direction::Backward,
            calendar,
            warnings,
        )?;
        periods.insert(0, stub);
    }
    if last < end {
        check_count(periods.len() + 1)?;
        let stub = specified_stub(
            Some(grid.roll),
            step,
            periods[periods.len() - 1].end,
            end,
            Direction::Forward,
            calendar,
            warnings,
        )?;
        periods.push(stub);
    }
    Ok(periods)
}

fn add_residual(
    periods: &mut Vec<DraftPeriod>,
    current: GridDate,
    next: GridDate,
    endpoint: GridDate,
    length: StubLength,
    direction: Direction,
) -> Result<()> {
    // 端を越えた格子点までの仮想の通常期間が、余りの参照期間になる。
    let reference = ordered(current, next);
    match length {
        StubLength::Short => {
            check_count(periods.len() + 1)?;
            let mut stub = DraftPeriod::between(current, endpoint);
            stub.kind = PeriodKind::ShortStub;
            stub.references.push(reference);
            periods.push(stub);
        }
        StubLength::Long => {
            let stub = periods
                .last_mut()
                .ok_or_else(|| invalid("a long stub requires an adjacent regular period"))?;
            // 隣の通常期間と仮想の通常期間を、日付の昇順に並べて参照枠にする。
            let regular = (stub.start, stub.end);
            stub.kind = PeriodKind::LongStub;
            match direction {
                Direction::Forward => {
                    stub.end = endpoint;
                    stub.references.extend([regular, reference]);
                }
                Direction::Backward => {
                    stub.start = endpoint;
                    stub.references.extend([reference, regular]);
                }
            }
        }
    }
    Ok(())
}

fn specified_stub(
    roll: Option<RollConvention>,
    step: TenorLength,
    anchor: GridDate,
    limit: Date,
    direction: Direction,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<DraftPeriod> {
    // FRN の初回スタブの参照枠は指定された通常開始日から通常のロールで遡る。
    let roll = roll.filter(|roll| *roll != RollConvention::Frn || direction == Direction::Forward);
    let mut grid = Grid::new(anchor, step, roll, direction);
    let mut stub = DraftPeriod::between(anchor, limit.into());
    while direction.before(grid.current.raw, limit) {
        check_count(stub.references.len() + 1)?;
        let current = grid.current;
        let next = grid.next(calendar, warnings)?;
        stub.references.push(ordered(current, next));
    }
    if direction == Direction::Backward {
        stub.references.reverse();
    }
    stub.kind = match stub.references.as_slice() {
        [(start, end)] if start.raw == stub.start.raw && end.raw == stub.end.raw => {
            PeriodKind::Regular
        }
        [_] => PeriodKind::ShortStub,
        _ => PeriodKind::LongStub,
    };
    Ok(stub)
}

fn apply_roll(date: Date, roll: RollConvention, unit: TenorUnit) -> Result<Date> {
    let (year, month, _) = date.ymd();
    match (roll, unit) {
        (RollConvention::OnDay(day), TenorUnit::Month | TenorUnit::Year) => {
            Date::from_ymd(year, month, day.0.min(date.month_end()?.ymd().2))
        }
        (RollConvention::EndOfMonth, TenorUnit::Month | TenorUnit::Year) => date.month_end(),
        (RollConvention::OnWeekday(weekday), TenorUnit::Week) => {
            date.checked_add_days(weekday_index(weekday) - weekday_index(date.weekday()))
        }
        (
            RollConvention::NthWeekday {
                week: WeekOfMonth::Last,
                weekday,
            },
            TenorUnit::Month | TenorUnit::Year,
        ) => {
            let last = date.month_end()?;
            last.checked_add_days(
                -(weekday_index(last.weekday()) - weekday_index(weekday)).rem_euclid(7),
            )
        }
        (RollConvention::NthWeekday { week, weekday }, TenorUnit::Month | TenorUnit::Year) => {
            let index = match week {
                WeekOfMonth::First | WeekOfMonth::Last => 1,
                WeekOfMonth::Second => 2,
                WeekOfMonth::Third => 3,
                WeekOfMonth::Fourth => 4,
            };
            Date::from_naive(
                NaiveDate::from_weekday_of_month_opt(year, u32::from(month), weekday, index)
                    .ok_or_else(date_error)?,
            )
        }
        _ => Ok(date),
    }
}

#[allow(clippy::cast_possible_wrap, reason = "曜日番号の範囲は 0〜6")]
fn weekday_index(weekday: chrono::Weekday) -> i32 {
    weekday.num_days_from_monday() as i32
}

fn frn_next(
    current: GridDate,
    size: i32,
    unit: TenorUnit,
    calendar: &Calendar,
    warnings: &mut Vec<Warning>,
) -> Result<GridDate> {
    let previous = current.adjusted.unwrap_or(current.raw);
    let next = shift_calendar_date(previous, size, unit, false)?;
    // 起点は未調整の effective date。以後は月の最終営業日に着いたらその慣行を引き継ぐ。
    if next.ymd().2 < previous.ymd().2
        || (current.adjusted.is_some() && is_last_business_day(previous, calendar, warnings)?)
    {
        let raw = next.month_end()?;
        let adjusted = adjust_date(raw, BusinessDayConvention::Preceding, calendar, warnings)?;
        if adjusted.ymd().1 != raw.ymd().1 {
            return Err(invalid("FRN target month has no business day"));
        }
        return Ok(GridDate {
            raw,
            adjusted: Some(adjusted),
        });
    }
    let adjusted = adjust_date(
        next,
        BusinessDayConvention::ModifiedFollowing,
        calendar,
        warnings,
    )?;
    Ok(GridDate {
        raw: next,
        adjusted: Some(adjusted),
    })
}

fn ordered(a: GridDate, b: GridDate) -> (GridDate, GridDate) {
    if a.raw <= b.raw { (a, b) } else { (b, a) }
}

fn nonempty_interval(start: Date, end: Date) -> Result<DateInterval> {
    if start >= end {
        return Err(invalid("period boundaries must strictly increase"));
    }
    DateInterval::new(start, end)
}
`````

`````
crates/rhodium-core/src/schedule/wire.rs
use serde::{Deserialize, Deserializer, Serialize, Serializer, de::Error as _};

use super::{
    AccrualDates, EndRule, ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, GeneratedPeriods,
    PaymentRule, PeriodCombination, RollConvention, ScheduleSpec, StartRule, StubSpec,
};
use crate::{
    Error, Result,
    calendar::CalendarSpec,
    date::Date,
    tenor::{DateOffset, Frequency, TenorLength},
};

#[derive(Serialize, Deserialize)]
#[serde(untagged)]
enum Start {
    Date(Date),
    Relative(RelativeStart),
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct RelativeStart {
    from_trade_date: DateOffset,
}

