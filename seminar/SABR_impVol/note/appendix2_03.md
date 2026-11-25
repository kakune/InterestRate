impl From<StartRule> for Start {
    fn from(rule: StartRule) -> Self {
        match rule {
            StartRule::Date(date) => Self::Date(date),
            StartRule::FromTradeDate(from_trade_date) => {
                Self::Relative(RelativeStart { from_trade_date })
            }
        }
    }
}

impl From<Start> for StartRule {
    fn from(start: Start) -> Self {
        match start {
            Start::Date(date) => Self::Date(date),
            Start::Relative(rule) => Self::FromTradeDate(rule.from_trade_date),
        }
    }
}

#[derive(Serialize, Deserialize)]
#[serde(untagged)]
enum End {
    Date(Date),
    Relative(RelativeEnd),
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct RelativeEnd {
    after_start: TenorLength,
}

impl From<EndRule> for End {
    fn from(rule: EndRule) -> Self {
        match rule {
            EndRule::Date(date) => Self::Date(date),
            EndRule::AfterStart(after_start) => Self::Relative(RelativeEnd { after_start }),
        }
    }
}

impl From<End> for EndRule {
    fn from(end: End) -> Self {
        match end {
            End::Date(date) => Self::Date(date),
            End::Relative(rule) => Self::AfterStart(rule.after_start),
        }
    }
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Generated {
    start: Start,
    end: End,
    step: Frequency,
    roll: Option<RollConvention>,
    stub: StubSpec,
    calendar: CalendarSpec,
    bdc: super::BusinessDayConvention,
    accrual: AccrualDates,
    payment: PaymentRule,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Explicit {
    periods: Vec<ExplicitPeriod>,
    #[serde(default)]
    on_holiday: ExplicitDatePolicy,
    #[serde(default)]
    combine: PeriodCombination,
    frequency: Option<Frequency>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ExplicitSpec {
    explicit: Explicit,
}

#[derive(Serialize, Deserialize)]
#[serde(untagged)]
enum Wire {
    Generated(Generated),
    Explicit(ExplicitSpec),
}

impl TryFrom<Wire> for ScheduleSpec {
    type Error = Error;

    fn try_from(wire: Wire) -> Result<Self> {
        match wire {
            Wire::Generated(value) => Ok(Self::Generated {
                periods: GeneratedPeriods::new(
                    value.start.into(),
                    value.end.into(),
                    value.step,
                    value.roll,
                    value.stub,
                    value.calendar,
                    value.bdc,
                    value.accrual,
                )?,
                payment: value.payment,
            }),
            Wire::Explicit(ExplicitSpec { explicit: value }) => ExplicitPeriods::new(
                value.periods,
                value.on_holiday,
                value.combine,
                value.frequency,
            )
            .map(Self::Explicit),
        }
    }
}

impl<'de> Deserialize<'de> for ScheduleSpec {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> std::result::Result<Self, D::Error> {
        Wire::deserialize(deserializer)?
            .try_into()
            .map_err(D::Error::custom)
    }
}

impl Serialize for ScheduleSpec {
    fn serialize<S: Serializer>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error> {
        let value = match self {
            Self::Generated { periods, payment } => Wire::Generated(Generated {
                start: periods.start.into(),
                end: periods.end.into(),
                step: periods.step,
                roll: periods.roll,
                stub: periods.stub,
                calendar: periods.calendar.clone(),
                bdc: periods.bdc,
                accrual: periods.accrual,
                payment: payment.clone(),
            }),
            Self::Explicit(periods) => Wire::Explicit(ExplicitSpec {
                explicit: Explicit {
                    periods: periods.periods.clone(),
                    on_holiday: periods.on_holiday.clone(),
                    combine: periods.combine,
                    frequency: periods.frequency,
                },
            }),
        };
        value.serialize(serializer)
    }
}

pub(super) mod weekday {
    use serde::{Deserialize, Deserializer, Serializer, de::Error as _};

    use crate::date::Weekday;

    #[expect(
        clippy::trivially_copy_pass_by_ref,
        reason = "serde の with が参照を渡す"
    )]
    pub(in crate::schedule) fn serialize<S: Serializer>(
        day: &Weekday,
        serializer: S,
    ) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&day.to_string().to_ascii_lowercase())
    }

    pub(in crate::schedule) fn deserialize<'de, D: Deserializer<'de>>(
        deserializer: D,
    ) -> Result<Weekday, D::Error> {
        String::deserialize(deserializer)?
            .parse()
            .map_err(D::Error::custom)
    }
}
`````

`````
crates/rhodium-core/src/tenor.rs
//! テナー、符号付きの日付オフセット、期の頻度。

use std::{cmp::Ordering, fmt, num::NonZeroU8, str::FromStr};

use crate::{Error, Result};

/// テナーと日付オフセットの単位。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TenorUnit {
    /// 暦日。
    Day,
    /// 週（7 暦日）。
    Week,
    /// 暦月。
    Month,
    /// 年（12 暦月）。
    Year,
    /// 営業日。
    BusinessDay,
}

impl TenorUnit {
    const fn suffix(self) -> &'static str {
        match self {
            Self::Day => "D",
            Self::Week => "W",
            Self::Month => "M",
            Self::Year => "Y",
            Self::BusinessDay => "BD",
        }
    }
}

/// 暦の長さか、短期の市場商品の期間。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Tenor {
    /// 取引日から翌営業日まで。
    Overnight,
    /// 取引日の翌営業日からその翌営業日まで。
    TomorrowNext,
    /// スポット日からその翌営業日まで。
    SpotNext,
    /// 正準化した量と単位。
    Length(TenorLength),
}

/// 数と単位の正準形。暦のゼロは `0D`、営業日のゼロは `0BD`。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TenorLength {
    amount: u16,
    unit: TenorUnit,
}

impl TenorLength {
    /// 量と単位を正準化する。
    #[must_use]
    pub const fn new(amount: u16, unit: TenorUnit) -> Self {
        let (amount, unit) = match unit {
            TenorUnit::BusinessDay => (amount, unit),
            _ if amount == 0 => (0, TenorUnit::Day),
            TenorUnit::Day if amount.is_multiple_of(7) => (amount / 7, TenorUnit::Week),
            TenorUnit::Month if amount.is_multiple_of(12) => (amount / 12, TenorUnit::Year),
            _ => (amount, unit),
        };
        Self { amount, unit }
    }

    /// 正準形での量。
    #[must_use]
    pub const fn amount(self) -> u16 {
        self.amount
    }

    /// 正準形での単位。
    #[must_use]
    pub const fn unit(self) -> TenorUnit {
        self.unit
    }

    fn expanded(self) -> (u32, TenorUnit) {
        let amount = u32::from(self.amount);
        match self.unit {
            TenorUnit::Week => (amount * 7, TenorUnit::Day),
            TenorUnit::Year => (amount * 12, TenorUnit::Month),
            unit => (amount, unit),
        }
    }

    fn from_expanded(amount: u32, unit: TenorUnit) -> Result<Self> {
        let (amount, unit) = match unit {
            TenorUnit::Day if amount.is_multiple_of(7) => (amount / 7, TenorUnit::Week),
            TenorUnit::Month if amount.is_multiple_of(12) => (amount / 12, TenorUnit::Year),
            _ => (amount, unit),
        };
        Ok(Self::new(
            u16::try_from(amount)
                .map_err(|_| Error::InvalidInput("tenor length overflow".into()))?,
            unit,
        ))
    }
}

impl Tenor {
    /// 同じ単位族のテナーを足し、結果を正準化する。`0D` は日・週と月・年で使える。
    ///
    /// # Errors
    ///
    /// 単位族が異なる、ON/TN/SN を含む、または結果が量の値域を超えると
    /// [`Error::InvalidInput`]。
    pub fn checked_add(self, other: Self) -> Result<Self> {
        let (Self::Length(left), Self::Length(right)) = (self, other) else {
            return Err(incompatible_tenors());
        };
        let (left, left_unit) = left.expanded();
        let (right, right_unit) = right.expanded();
        let unit = match (left, left_unit, right, right_unit) {
            (_, a, _, b) if a == b => a,
            (0, TenorUnit::Day, _, TenorUnit::Month) | (_, TenorUnit::Month, 0, TenorUnit::Day) => {
                TenorUnit::Month
            }
            _ => return Err(incompatible_tenors()),
        };
        // u16 の年を月に広げて足しても u32 に収まる。
        TenorLength::from_expanded(left + right, unit).map(Self::Length)
    }

    /// テナーの長さを整数倍して正準化する。
    ///
    /// # Errors
    ///
    /// ON/TN/SN を指定したか、結果が量の値域を超えると [`Error::InvalidInput`]。
    pub fn times(self, factor: u16) -> Result<Self> {
        let Self::Length(length) = self else {
            return Err(incompatible_tenors());
        };
        let (amount, unit) = length.expanded();
        let amount = amount
            .checked_mul(u32::from(factor))
            .ok_or_else(|| Error::InvalidInput("tenor length overflow".into()))?;
        TenorLength::from_expanded(amount, unit).map(Self::Length)
    }

    // 名目日数を 1/80 日単位の整数で比べる。
    fn nominal_eightieths(self) -> u32 {
        match self {
            Self::Overnight => 80,
            Self::TomorrowNext => 160,
            Self::SpotNext => 240,
            Self::Length(length) => {
                let scale = match length.unit {
                    TenorUnit::Day => 80,
                    TenorUnit::Week => 560,
                    TenorUnit::Month => 2_435,
                    TenorUnit::Year => 29_220,
                    TenorUnit::BusinessDay => 112,
                };
                u32::from(length.amount) * scale
            }
        }
    }
}

impl FromStr for Tenor {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        let tenor = match text {
            "ON" => Self::Overnight,
            "TN" => Self::TomorrowNext,
            "SN" => Self::SpotNext,
            _ => {
                let (digits, unit) = split_tenor(text).ok_or_else(invalid_tenor)?;
                if !unsigned_digits(digits) {
                    return Err(invalid_tenor());
                }
                let amount = digits.parse().map_err(|_| invalid_tenor())?;
                Self::Length(TenorLength::new(amount, unit))
            }
        };
        Ok(tenor)
    }
}

impl fmt::Display for Tenor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Overnight => f.write_str("ON"),
            Self::TomorrowNext => f.write_str("TN"),
            Self::SpotNext => f.write_str("SN"),
            Self::Length(length) => write!(f, "{}{}", length.amount, length.unit.suffix()),
        }
    }
}

impl Ord for Tenor {
    fn cmp(&self, other: &Self) -> Ordering {
        self.nominal_eightieths()
            .cmp(&other.nominal_eightieths())
            .then_with(|| self.to_string().cmp(&other.to_string()))
    }
}

impl PartialOrd for Tenor {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

serde_text!(Tenor);

impl FromStr for TenorLength {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        match text.parse()? {
            Tenor::Length(length) => Ok(length),
            _ => Err(invalid_tenor()),
        }
    }
}

impl fmt::Display for TenorLength {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        Tenor::Length(*self).fmt(f)
    }
}

serde_text!(TenorLength);

/// 日付を前後にずらす量と単位。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DateOffset {
    amount: i16,
    unit: TenorUnit,
}

impl DateOffset {
    /// 符号付きの量と単位を指定する。
    #[must_use]
    pub const fn new(amount: i16, unit: TenorUnit) -> Self {
        Self { amount, unit }
    }

    /// 符号付きの量。
    #[must_use]
    pub const fn amount(self) -> i16 {
        self.amount
    }

    /// 日付を進める単位。
    #[must_use]
    pub const fn unit(self) -> TenorUnit {
        self.unit
    }
}

impl FromStr for DateOffset {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        let invalid = || Error::InvalidData(format!("invalid date offset: {text}"));
        let (amount, unit) = split_tenor(text).ok_or_else(invalid)?;
        if !unsigned_digits(amount.strip_prefix(['+', '-']).unwrap_or(amount)) {
            return Err(invalid());
        }
        Ok(Self::new(amount.parse().map_err(|_| invalid())?, unit))
    }
}

impl fmt::Display for DateOffset {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.amount, self.unit.suffix())
    }
}

serde_text!(DateOffset);

/// 一回限りか、正のテナーを刻みとする頻度。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Frequency(Option<TenorLength>);

impl Frequency {
    /// 全期間で一回。
    pub const ONCE: Self = Self(None);
    /// 一年に一回。
    pub const ANNUAL: Self = Self(Some(TenorLength::new(1, TenorUnit::Year)));

    /// 正のテナーを刻みとする頻度を作る。
    ///
    /// # Errors
    ///
    /// ゼロまたは ON/TN/SN は [`Error::InvalidData`]。
    pub fn every(step: Tenor) -> Result<Self> {
        match step {
            Tenor::Length(length) if length.amount != 0 => Ok(Self(Some(length))),
            _ => Err(Error::InvalidData(
                "frequency requires a positive tenor length".into(),
            )),
        }
    }

    /// 繰り返すときの刻み。
    #[must_use]
    pub const fn step(self) -> Option<TenorLength> {
        self.0
    }

    /// 複利の頻度として使える場合の一年あたりの回数。
    #[must_use]
    pub fn periods_per_year(self) -> Option<NonZeroU8> {
        let step = self.step()?;
        let count = match (step.amount, step.unit) {
            (1, TenorUnit::Year) => 1,
            (amount, TenorUnit::Month) if 12_u16.is_multiple_of(amount) => 12 / amount,
            _ => return None,
        };
        NonZeroU8::new(u8::try_from(count).ok()?)
    }
}

impl FromStr for Frequency {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        if text == "ONCE" {
            Ok(Self::ONCE)
        } else {
            Self::every(text.parse()?)
        }
    }
}

impl fmt::Display for Frequency {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.step() {
            Some(step) => write!(f, "{}", Tenor::Length(step)),
            None => f.write_str("ONCE"),
        }
    }
}

serde_text!(Frequency);

fn unsigned_digits(text: &str) -> bool {
    !text.is_empty()
        && (text.len() == 1 || !text.starts_with('0'))
        && text.bytes().all(|byte| byte.is_ascii_digit())
}

fn split_tenor(text: &str) -> Option<(&str, TenorUnit)> {
    [
        TenorUnit::BusinessDay,
        TenorUnit::Day,
        TenorUnit::Week,
        TenorUnit::Month,
        TenorUnit::Year,
    ]
    .into_iter()
    .find_map(|unit| {
        text.strip_suffix(unit.suffix())
            .map(|amount| (amount, unit))
    })
}

fn invalid_tenor() -> Error {
    Error::InvalidData("invalid tenor".into())
}

fn incompatible_tenors() -> Error {
    Error::InvalidInput("incompatible tenor lengths".into())
}
`````

`````
crates/rhodium-core/src/time.rs
//! UTC の瞬間、現地日時、地域の規則、モデル時間。

mod zone;

use std::{fmt, str::FromStr, sync::Arc};

use chrono::{DateTime, FixedOffset, NaiveDateTime, NaiveTime, SecondsFormat, Timelike, Utc};

use crate::{Error, Result, date::Date};

pub use zone::{ZoneOffset, ZoneRules};

/// UTC の瞬間。UTC の年は 1〜9999、精度はマイクロ秒。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Instant(DateTime<Utc>);

impl Instant {
    /// Unix epoch からの符号付きマイクロ秒で瞬間を作る。
    ///
    /// # Errors
    ///
    /// UTC の年が値域を出る場合は [`Error::InvalidInput`]。
    pub fn from_unix_micros(micros: i64) -> Result<Self> {
        Self::from_utc(DateTime::from_timestamp_micros(micros).ok_or_else(time_error)?)
    }

    /// Unix epoch からの符号付きマイクロ秒。
    #[must_use]
    pub const fn unix_micros(self) -> i64 {
        self.0.timestamp_micros()
    }

    fn from_utc(datetime: DateTime<Utc>) -> Result<Self> {
        Date::from_naive(datetime.date_naive())?;
        check_precision(datetime.nanosecond())?;
        Ok(Self(datetime))
    }

    const fn naive_utc(self) -> NaiveDateTime {
        self.0.naive_utc()
    }
}

impl FromStr for Instant {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        check_fraction(text)?;
        let parsed = DateTime::parse_from_rfc3339(text).map_err(|_| time_error())?;
        UtcOffset::from_seconds(parsed.offset().local_minus_utc())?;
        Self::from_utc(parsed.to_utc())
    }
}

impl fmt::Display for Instant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0.to_rfc3339_opts(SecondsFormat::AutoSi, true))
    }
}

serde_text!(Instant);

/// 現地の一日の時刻。精度はマイクロ秒。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LocalTime(NaiveTime);

impl LocalTime {
    /// 時・分・秒・マイクロ秒から現地時刻を作る。
    ///
    /// # Errors
    ///
    /// 時刻が存在しないか、100 万マイクロ秒以上を指定すると [`Error::InvalidInput`]。
    pub fn from_hms_micro(hour: u32, minute: u32, second: u32, micros: u32) -> Result<Self> {
        Self::from_naive(
            NaiveTime::from_hms_micro_opt(hour, minute, second, micros).ok_or_else(time_error)?,
        )
    }

    fn from_naive(time: NaiveTime) -> Result<Self> {
        check_precision(time.nanosecond())?;
        Ok(Self(time))
    }
}

impl FromStr for LocalTime {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        check_clock_shape(text)?;
        check_fraction(text)?;
        Self::from_naive(NaiveTime::parse_from_str(text, "%H:%M:%S%.f").map_err(|_| time_error())?)
    }
}

impl fmt::Display for LocalTime {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

serde_text!(LocalTime);

/// 地域を解決する前の現地の日付と時刻。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LocalDateTime(NaiveDateTime);

impl LocalDateTime {
    /// 検証済みの日付と現地時刻を組み合わせる。
    #[must_use]
    pub const fn new(date: Date, time: LocalTime) -> Self {
        Self(date.naive().and_time(time.0))
    }

    /// 現地の日付。
    #[must_use]
    pub const fn date(self) -> Date {
        Date::from_valid_naive(self.0.date())
    }

    /// 現地の時刻。
    #[must_use]
    pub const fn time(self) -> LocalTime {
        LocalTime(self.0.time())
    }

    /// 明示された UTC オフセットから瞬間を求める。
    ///
    /// # Errors
    ///
    /// UTC の年が値域を出る場合は [`Error::InvalidInput`]。
    pub fn at_offset(self, offset: UtcOffset) -> Result<Instant> {
        let utc = self.0.checked_sub_offset(offset.0).ok_or_else(time_error)?;
        Instant::from_utc(utc.and_utc())
    }
}

impl FromStr for LocalDateTime {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        let (date, time) = text.split_once('T').ok_or_else(time_error)?;
        Ok(Self::new(date.parse()?, time.parse()?))
    }
}

impl fmt::Display for LocalDateTime {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}T{}", self.date(), self.time())
    }
}

serde_text!(LocalDateTime);

/// 地域の規則を参照する IANA タイムゾーン名。
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TzId(Arc<str>);

impl TzId {
    /// 地域名。
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl FromStr for TzId {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        if text.split('/').all(|part| {
            !matches!(part, "" | "." | "..")
                && part.bytes().all(|byte| {
                    byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'+' | b'.')
                })
        }) {
            Ok(Self(text.into()))
        } else {
            Err(Error::InvalidData("invalid time zone name".into()))
        }
    }
}

impl fmt::Display for TzId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

serde_text!(TzId);

/// 各営業日に繰り返す事象の時刻と地域。
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LocalEventTime {
    /// 現地時刻。
    #[serde(deserialize_with = "deserialize_event_time")]
    pub time: LocalTime,
    /// 解決に使う規則の地域。
    pub zone: TzId,
}

// 設定には HH:MM も書ける。LocalTime の文字列表現は秒までの形を保つ。
pub(crate) fn deserialize_event_time<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> std::result::Result<LocalTime, D::Error> {
    let text = <String as serde::Deserialize>::deserialize(deserializer)?;
    let parsed = if text.len() == 5 {
        format!("{text}:00").parse()
    } else {
        text.parse()
    };
    parsed.map_err(serde::de::Error::custom)
}

impl LocalEventTime {
    /// 営業日調整を済ませた日付と現地時刻を組み合わせる。
    #[must_use]
    pub const fn on_date(&self, date: Date) -> LocalDateTime {
        LocalDateTime::new(date, self.time)
    }
}

/// UTC から現地に進める秒数。範囲は ±18 時間。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct UtcOffset(FixedOffset);

impl UtcOffset {
    /// 符号付きの秒数でオフセットを作る。
    ///
    /// # Errors
    ///
    /// ±18 時間の範囲外は [`Error::InvalidInput`]。
    pub fn from_seconds(seconds: i32) -> Result<Self> {
        if !(-64_800..=64_800).contains(&seconds) {
            return Err(Error::InvalidInput("UTC offset exceeds ±18 hours".into()));
        }
        Ok(Self(FixedOffset::east_opt(seconds).ok_or_else(time_error)?))
    }

    /// UTC から現地に進める符号付き秒数。
    #[must_use]
    pub const fn seconds(self) -> i32 {
        self.0.local_minus_utc()
    }
}

/// 評価時点からの実経過時間を 365 日で割ったモデルの年数。
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct ModelTime(f64);

impl ModelTime {
    /// 符号付きのモデル年数。
    #[must_use]
    pub const fn years(self) -> f64 {
        self.0
    }
}

/// 評価時点からの実経過時間でモデルの時間軸を求める。
#[must_use]
pub fn model_time(origin: Instant, at: Instant) -> ModelTime {
    // UTC 年 1〜9999 の全期間でも差は i64 のマイクロ秒に収まる。
    ModelTime((at.unix_micros() - origin.unix_micros()) as f64 / 31_536_000_000_000.0)
}

// chrono の書式解析は数値の前の空白と 1 桁の欄を許すため、RFC 3339 の partial-time の形
// （2 桁ずつの時・分・秒と、任意の "." と 1 桁以上の小数秒）を先に確かめる。
fn check_clock_shape(text: &str) -> Result<()> {
    let (clock, fraction) = text.as_bytes().split_at_checked(8).ok_or_else(time_error)?;
    let clock_ok = clock.iter().enumerate().all(|(index, byte)| match index {
        2 | 5 => *byte == b':',
        _ => byte.is_ascii_digit(),
    });
    let fraction_ok = match fraction {
        [] => true,
        [b'.', digits @ ..] => !digits.is_empty() && digits.iter().all(u8::is_ascii_digit),
        _ => false,
    };
    if clock_ok && fraction_ok {
        Ok(())
    } else {
        Err(time_error())
    }
}

fn check_fraction(text: &str) -> Result<()> {
    // chrono が 10 桁目以降を切り捨てる前に、表せない値を検出する。
    if text.split_once('.').is_some_and(|(_, fraction)| {
        fraction
            .bytes()
            .take_while(u8::is_ascii_digit)
            .skip(6)
            .any(|byte| byte != b'0')
    }) {
        return Err(time_error());
    }
    Ok(())
}

fn check_precision(nanos: u32) -> Result<()> {
    if nanos >= 1_000_000_000 || !nanos.is_multiple_of(1_000) {
        return Err(time_error());
    }
    Ok(())
}

fn time_error() -> Error {
    Error::InvalidInput("invalid date/time or time not representable in microseconds".into())
}
`````

`````
crates/rhodium-core/src/time/zone.rs
use std::{cmp::Reverse, iter, ops::RangeInclusive, sync::Arc};

use chrono::{
    FixedOffset, LocalResult, NaiveDate, NaiveDateTime, NaiveTime, Offset, TimeDelta, TimeZone,
};

use super::{Instant, LocalDateTime, UtcOffset, time_error};
use crate::{Error, Result, date::Date};

/// 地域の UTC 適用範囲、初期オフセット、遷移表。
#[derive(Debug, Clone)]
pub struct ZoneRules(Arc<Rules>);

#[derive(Debug)]
struct Rules {
    coverage: RangeInclusive<NaiveDateTime>,
    initial: UtcOffset,
    transitions: Vec<(NaiveDateTime, UtcOffset)>,
    // 同じ現地日時に対して UTC が早い順に調べるため、オフセットの降順。
    offsets: Vec<UtcOffset>,
}

/// chrono の日時が保持するオフセットと、そのオフセットを導いた地域の規則。
#[derive(Debug, Clone)]
pub struct ZoneOffset {
    offset: UtcOffset,
    rules: ZoneRules,
}

impl Offset for ZoneOffset {
    fn fix(&self) -> FixedOffset {
        self.offset.0
    }
}

impl ZoneRules {
    /// UTC の両端を含む適用範囲と、UTC の昇順で並べた遷移から規則を作る。
    /// 初期オフセットは開始の瞬間から、各遷移のオフセットはその瞬間から有効。
    ///
    /// # Errors
    ///
    /// 適用範囲が空、遷移が範囲内の狭義昇順でない、開始の瞬間に重なる、
    /// または同じ現地日時に三つ以上の UTC 候補があると [`Error::InvalidData`]。
    pub fn new(
        coverage: RangeInclusive<Instant>,
        initial: UtcOffset,
        transitions: Vec<(Instant, UtcOffset)>,
    ) -> Result<Self> {
        if coverage.is_empty() {
            return Err(Error::InvalidData("zone coverage is empty".into()));
        }
        let start = coverage.start().naive_utc();
        let end = coverage.end().naive_utc();
        let mut previous = start;
        let transitions = transitions
            .into_iter()
            .map(|(at, offset)| {
                let at = at.naive_utc();
                if at <= previous || at > end {
                    return Err(Error::InvalidData("invalid zone transition time".into()));
                }
                previous = at;
                Ok((at, offset))
            })
            .collect::<Result<Vec<_>>>()?;
        check_local_overlap(start, end, initial, &transitions)?;
        let mut offsets: Vec<_> = iter::once(initial)
            .chain(transitions.iter().map(|&(_, offset)| offset))
            .collect();
        offsets.sort_unstable_by_key(|offset| Reverse(offset.seconds()));
        offsets.dedup();
        Ok(Self(Arc::new(Rules {
            coverage: start..=end,
            initial,
            transitions,
            offsets,
        })))
    }

    /// 適用範囲内の UTC の瞬間で有効なオフセット。
    ///
    /// # Errors
    ///
    /// 瞬間が適用範囲外なら [`Error::MissingData`]。
    pub fn offset_at(&self, at: Instant) -> Result<UtcOffset> {
        self.covered_offset(at.naive_utc()).ok_or_else(missing_zone)
    }

    /// 適用範囲内の瞬間を現地日時に変換する。
    ///
    /// # Errors
    ///
    /// 規則の適用範囲外なら [`Error::MissingData`]。
    /// 変換後の現地の年が値域外なら [`Error::InvalidInput`]。
    pub fn to_local(&self, at: Instant) -> Result<LocalDateTime> {
        let offset = self.offset_at(at)?;
        let local = at
            .naive_utc()
            .checked_add_offset(offset.0)
            .ok_or_else(time_error)?;
        Date::from_naive(local.date())?;
        Ok(LocalDateTime(local))
    }

    /// 利用者が指定した現地日時を地域の規則で解決する。
    ///
    /// # Errors
    ///
    /// 規則の適用範囲外なら [`Error::MissingData`]。欠落・重複する現地日時に
    /// オフセットを省いたか、明示したオフセットが規則と違うと [`Error::InvalidInput`]。
    pub fn resolve_input(
        &self,
        local: LocalDateTime,
        offset: Option<UtcOffset>,
    ) -> Result<Instant> {
        if let Some(offset) = offset {
            let at = local.at_offset(offset)?;
            return if self.offset_at(at)? == offset {
                Ok(at)
            } else {
                Err(Error::InvalidInput(
                    "UTC offset does not match zone rules".into(),
                ))
            };
        }
        match self.from_local_datetime(&local.0) {
            LocalResult::Single(at) => Instant::from_utc(at.to_utc()),
            LocalResult::Ambiguous(..) => Err(Error::InvalidInput(
                "ambiguous local time requires an explicit UTC offset".into(),
            )),
            LocalResult::None => Err(self.unmapped_local_error(local.0)),
        }
    }

    /// 慣行から得た現地日時を解決する。重複なら早い瞬間、欠落ならギャップ幅だけ進める。
    ///
    /// # Errors
    ///
    /// 解決した日時が適用範囲外なら [`Error::MissingData`]、日付の値域外や、
    /// ギャップを進んだ後も現地日時が存在しない場合は [`Error::InvalidInput`]。
    pub fn resolve_convention(&self, local: LocalDateTime) -> Result<Instant> {
        if let Some(at) = self.from_local_datetime(&local.0).earliest() {
            return Instant::from_utc(at.to_utc());
        }
        let width = self.gap_width(local.0).ok_or_else(missing_zone)?;
        let shifted = local.0.checked_add_signed(width).ok_or_else(time_error)?;
        Date::from_naive(shifted.date())?;
        let at = self
            .from_local_datetime(&shifted)
            .earliest()
            .ok_or_else(|| self.unmapped_local_error(shifted))?;
        Instant::from_utc(at.to_utc())
    }

    fn offset_for_utc(&self, utc: NaiveDateTime) -> UtcOffset {
        let index = self.0.transitions.partition_point(|&(at, _)| at <= utc);
        index
            .checked_sub(1)
            .and_then(|index| self.0.transitions.get(index))
            .map_or(self.0.initial, |&(_, offset)| offset)
    }

    fn covered_offset(&self, utc: NaiveDateTime) -> Option<UtcOffset> {
        self.0
            .coverage
            .contains(&utc)
            .then(|| self.offset_for_utc(utc))
    }

    // 現地日時を含む欠落の幅。オフセットが増える遷移では、遷移の瞬間を旧オフセットで表した
    // 現地日時から新オフセットで表した現地日時までが欠落する。減る遷移は start > end で一致しない。
    fn gap_width(&self, local: NaiveDateTime) -> Option<TimeDelta> {
        let old_offsets =
            iter::once(self.0.initial).chain(self.0.transitions.iter().map(|&(_, offset)| offset));
        self.0
            .transitions
            .iter()
            .zip(old_offsets)
            .find_map(|(&(at, new), old)| {
                let start = at.checked_add_offset(old.0)?;
                let end = at.checked_add_offset(new.0)?;
                (start <= local && local < end).then(|| end.signed_duration_since(start))
            })
    }

    fn unmapped_local_error(&self, local: NaiveDateTime) -> Error {
        if self.gap_width(local).is_some() {
            Error::InvalidInput("local time is in a zone transition gap".into())
        } else {
            missing_zone()
        }
    }

    fn with_rules(&self, offset: UtcOffset) -> ZoneOffset {
        ZoneOffset {
            offset,
            rules: self.clone(),
        }
    }
}

impl TimeZone for ZoneRules {
    type Offset = ZoneOffset;

    fn from_offset(offset: &Self::Offset) -> Self {
        offset.rules.clone()
    }

    fn offset_from_local_date(&self, local: &NaiveDate) -> LocalResult<Self::Offset> {
        self.offset_from_local_datetime(&local.and_time(NaiveTime::MIN))
    }

    fn offset_from_local_datetime(&self, local: &NaiveDateTime) -> LocalResult<Self::Offset> {
        // 現地日時からオフセットを引いた UTC で、そのオフセット自身が有効なものだけが候補になる。
        // 一つの UTC は一つの区間にしか属さないので、異なるオフセットごとに候補は高々一つで、
        // 構築時の検査により候補は二つ以下。オフセットの降順に調べるので、先に UTC が早い候補が出る。
        let mut matching = self.0.offsets.iter().copied().filter(|offset| {
            local
                .checked_sub_offset(offset.0)
                .and_then(|utc| self.covered_offset(utc))
                == Some(*offset)
        });
        match (matching.next(), matching.next()) {
            (Some(early), Some(late)) => {
                LocalResult::Ambiguous(self.with_rules(early), self.with_rules(late))
            }
            (Some(offset), None) => LocalResult::Single(self.with_rules(offset)),
            (None, _) => LocalResult::None,
        }
    }

    fn offset_from_utc_date(&self, utc: &NaiveDate) -> Self::Offset {
        self.offset_from_utc_datetime(&utc.and_time(NaiveTime::MIN))
    }

    fn offset_from_utc_datetime(&self, utc: &NaiveDateTime) -> Self::Offset {
        // chrono の UTC メソッドは失敗を返せない。範囲の検査が必要なら offset_at を使う。
        self.with_rules(self.offset_for_utc(*utc))
    }
}

fn check_local_overlap(
    start: NaiveDateTime,
    end: NaiveDateTime,
    initial: UtcOffset,
    transitions: &[(NaiveDateTime, UtcOffset)],
) -> Result<()> {
    let starts = iter::once((start, initial)).chain(transitions.iter().copied());
    // UTC の最後のナノ秒も含めた半開区間を走査する。同時刻の終端を始端より先に数える。
    let exclusive_end = end
        .checked_add_signed(TimeDelta::nanoseconds(1))
        .ok_or_else(time_error)?;
    let ends = transitions
        .iter()
        .map(|&(at, _)| at)
        .chain(iter::once(exclusive_end));
    let mut edges = Vec::with_capacity((transitions.len() + 1) * 2);
    for ((start, offset), end) in starts.zip(ends) {
        edges.push((
            start.checked_add_offset(offset.0).ok_or_else(time_error)?,
            true,
        ));
        edges.push((
            end.checked_add_offset(offset.0).ok_or_else(time_error)?,
            false,
        ));
    }
    edges.sort_unstable();
    let mut candidates = 0;
    for (_, begins) in edges {
        candidates += if begins { 1 } else { -1 };
        if candidates > 2 {
            return Err(Error::InvalidData(
                "local time has more than two UTC candidates".into(),
            ));
        }
    }
    Ok(())
}

fn missing_zone() -> Error {
    Error::MissingData("outside the coverage of the zone rules".into())
}
`````

`````
crates/rhodium-core/tests/calendar.rs
//! 休日の合成と、営業日を使った日付計算を確認する。

use std::collections::BTreeMap;

use rhodium_core::{
    Error, Result, Warning,
    calendar::{
        AdvanceRule, BaseCalendar, BusinessDayConvention as Bdc, Calendar, CalendarSpec,
        EndOfMonthRule, Extrapolation, WeekendRule, adjust, advance, money_market_period,
    },
    date::{Date, DateInterval, Weekday},
    id::CalendarId,
    tenor::Tenor,
};

fn base(
    known_start: &str,
    known_end: &str,
    holidays: &[&str],
    extrapolation: Extrapolation,
) -> Result<BaseCalendar> {
    BaseCalendar::new(
        known_start.parse()?..=known_end.parse()?,
        WeekendRule::new([Weekday::Sat, Weekday::Sun]),
        holidays
            .iter()
            .map(|date| date.parse())
            .collect::<Result<_>>()?,
        extrapolation,
    )
}

fn calendar(id: &str, base: BaseCalendar) -> Result<Calendar> {
    let id: CalendarId = id.parse()?;
    Calendar::compile(
        &CalendarSpec::Base(id.clone()),
        &BTreeMap::from([(id, base)]),
    )
}

#[test]
fn calendar_composition_and_extrapolation() -> Result<()> {
    let jpto: CalendarId = "JPTO".parse()?;
    let usny: CalendarId = "USNY".parse()?;
    let jp = CalendarSpec::Base(jpto.clone());
    let us = CalendarSpec::Base(usny.clone());
    let union = CalendarSpec::Union(vec![jp.clone(), us.clone()]);
    let intersection = CalendarSpec::Intersection(vec![jp.clone(), us.clone()]);
    let bases = BTreeMap::from([
        (
            jpto.clone(),
            base(
                "2027-01-01",
                "2027-01-22",
                &["2027-01-11"],
                Extrapolation::WeekendsOnly,
            )?,
        ),
        (
            usny,
            base(
                "2027-01-05",
                "2027-01-31",
                &["2027-01-18"],
                Extrapolation::Strict,
            )?,
        ),
    ]);
    let both_open = Calendar::compile(&union, &bases)?;
    let either_open = Calendar::compile(&intersection, &bases)?;

    // 片方だけの祝日は和なら休日、共通部分なら営業日。土曜はどちらも休日。
    for (date, business_in_union, business_in_intersection) in [
        ("2027-01-09", false, false),
        ("2027-01-11", false, true),
        ("2027-01-12", true, true),
        ("2027-01-18", false, true),
    ] {
        assert_eq!(
            both_open.is_business_day(date.parse()?)?,
            (business_in_union, vec![])
        );
        assert_eq!(
            either_open.is_business_day(date.parse()?)?,
            (business_in_intersection, vec![])
        );
    }

    // 祝日を一度営業日に変え、外側の式でまた休日にした。土曜の除去は週末より優先する。
    let monday = "2027-01-11".parse()?;
    let saturday = "2027-01-09".parse()?;
    let adjusted = CalendarSpec::Adjusted {
        inner: Box::new(CalendarSpec::Adjusted {
            inner: Box::new(union),
            add: vec![],
            remove: vec![monday, saturday],
        }),
        add: vec![monday],
        remove: vec![],
    };
    let adjusted = Calendar::compile(&adjusted, &bases)?;
    assert_eq!(adjusted.is_business_day(monday)?, (false, vec![]));
    assert_eq!(adjusted.is_business_day(saturday)?, (true, vec![]));

    // 既知範囲が異なる場合も各暦を照会する。二度参照した暦の外挿警告は一件になる。
    let repeated = Calendar::compile(
        &CalendarSpec::Union(vec![jp.clone(), us, jp.clone()]),
        &bases,
    )?;
    let outside = "2027-01-25".parse()?;
    assert_eq!(
        repeated.is_business_day(outside)?,
        (
            true,
            vec![Warning::CalendarExtrapolated {
                calendar: jpto,
                date: outside
            }]
        )
    );
    assert!(matches!(
        repeated.is_business_day("2027-02-01".parse()?),
        Err(Error::MissingData(_))
    ));

    // 不明な暦、要素が一つの和、休日の追加と除去の重複は式の解決時に拒否する。
    for invalid in [
        CalendarSpec::Base("NONE".parse()?),
        CalendarSpec::Union(vec![jp.clone()]),
        CalendarSpec::Adjusted {
            inner: Box::new(jp),
            add: vec![monday],
            remove: vec![monday],
        },
    ] {
        assert!(matches!(
            Calendar::compile(&invalid, &bases),
            Err(Error::InvalidData(_))
        ));
    }
    assert!(matches!(
        "日本".parse::<CalendarId>(),
        Err(Error::InvalidData(_))
    ));
    Ok(())
}

#[test]
fn calendar_algebra_laws_hold_day_by_day() -> Result<()> {
    // 月曜・火曜・水曜をそれぞれの祝日にした三つの暦で、休日集合の和と共通部分の法則が
    // 営業日の判定として成り立つことを、既知範囲の全日で確かめる。
    let known_start = "2026-12-01";
    let known_end = "2027-04-30";
    let [jpto, usny, gblo] = ["JPTO", "USNY", "GBLO"].map(str::parse::<CalendarId>);
    let (jpto, usny, gblo) = (jpto?, usny?, gblo?);
    let bases = BTreeMap::from([
        (
            jpto.clone(),
            base(
                known_start,
                known_end,
                &["2027-01-04"],
                Extrapolation::Strict,
            )?,
        ),
        (
            usny.clone(),
            base(
                known_start,
                known_end,
                &["2027-01-05"],
                Extrapolation::Strict,
            )?,
        ),
        (
            gblo.clone(),
            base(
                known_start,
                known_end,
                &["2027-01-06"],
                Extrapolation::Strict,
            )?,
        ),
    ]);
    let (a, b, c) = (
        CalendarSpec::Base(jpto),
        CalendarSpec::Base(usny),
        CalendarSpec::Base(gblo),
    );
    let start: Date = known_start.parse()?;
    let days = known_end.parse::<Date>()?.days_since(start);
    // 二つの式が既知範囲のすべての日で同じ営業日の判定を返すことを確かめる。
    let same_days = |left: &CalendarSpec, right: &CalendarSpec| -> Result<()> {
        let left = Calendar::compile(left, &bases)?;
        let right = Calendar::compile(right, &bases)?;
        for offset in 0..=days {
            let date = start.checked_add_days(offset)?;
            assert_eq!(
                left.is_business_day(date)?,
                right.is_business_day(date)?,
                "{date}"
            );
        }
        Ok(())
    };

    // 和と共通部分のそれぞれで、結合・交換をまとめた並べ替えと、冪等を確かめる。
    let groups: [fn(Vec<CalendarSpec>) -> CalendarSpec; 2] =
        [CalendarSpec::Union, CalendarSpec::Intersection];
    for group in groups {
        same_days(
            &group(vec![a.clone(), group(vec![b.clone(), c.clone()])]),
            &group(vec![group(vec![c.clone(), a.clone()]), b.clone()]),
        )?;
        same_days(&group(vec![a.clone(), a.clone()]), &a)?;
    }
    let union = CalendarSpec::Union;
    let intersection = CalendarSpec::Intersection;
    // 吸収: A ∪ (A ∩ B) = A、A ∩ (A ∪ B) = A。
    same_days(
        &union(vec![a.clone(), intersection(vec![a.clone(), b.clone()])]),
        &a,
    )?;
    same_days(
        &intersection(vec![a.clone(), union(vec![a.clone(), b.clone()])]),
        &a,
    )?;
    // 分配: A ∪ (B ∩ C) = (A ∪ B) ∩ (A ∪ C)。
    same_days(
        &union(vec![a.clone(), intersection(vec![b.clone(), c.clone()])]),
        &intersection(vec![
            union(vec![a.clone(), b.clone()]),
            union(vec![a.clone(), c]),
        ]),
    )?;

    // 営業日調整は和に分配しない。火曜から Preceding で、和の暦は月曜と火曜の祝日を
    // 飛ばして前週の金曜に戻るが、暦ごとに調整すると USNY の調整が JPTO の祝日の月曜に着く。
    let tuesday = "2027-01-05".parse()?;
    let joint = Calendar::compile(&union(vec![a.clone(), b.clone()]), &bases)?;
    assert_eq!(
        adjust(tuesday, Bdc::Preceding, &joint)?,
        ("2027-01-01".parse()?, vec![])
    );
    let jp = Calendar::compile(&a, &bases)?;
    let us = Calendar::compile(&b, &bases)?;
    let staggered = adjust(adjust(tuesday, Bdc::Preceding, &jp)?.0, Bdc::Preceding, &us)?.0;
    assert_eq!(staggered, "2027-01-04".parse()?);
    Ok(())
}

#[test]
fn business_day_count_respects_bit_boundaries_and_half_open_ranges() -> Result<()> {
    let cal = calendar(
        "JPTO",
        base(
            "2027-01-01",
            "2027-03-31",
            &["2027-01-11", "2027-02-11"],
            Extrapolation::Strict,
        )?,
    )?;
    // 70 日の区間は 10 週で平日 50 日。祝日を 2 日引くと 48 日。
    let start = "2027-01-03".parse()?;
    let end = "2027-03-14".parse()?;
    assert_eq!(cal.business_days_between(start, end)?, (48, vec![]));
    // 1 月 1 日から 64 日は 9 週と金曜 1 日。祝日を引くと 44 日。
    let first: Date = "2027-01-01".parse()?;
    assert_eq!(
        cal.business_days_between(first, first.checked_add_days(64)?)?,
        (44, vec![])
    );
    let last = "2027-03-31".parse()?;
    let after = "2027-04-01".parse()?;
    assert_eq!(cal.business_days_between(last, after)?, (1, vec![]));
    assert_eq!(cal.business_days_between(after, after)?, (0, vec![]));
    assert!(matches!(
        cal.business_days_between(first, "2027-04-02".parse()?),
        Err(Error::MissingData(_))
    ));
    assert!(matches!(
        cal.business_days_between(end, start),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn advance_respects_business_days_and_month_end() -> Result<()> {
    let cal = calendar(
        "JPTO",
        base(
            "2027-01-01",
            "2028-03-31",
            &["2027-07-14"],
            Extrapolation::Strict,
        )?,
    )?;
    // Nearest は曜日で向きを決める。水曜の祝日は翌日が営業日でも火曜に戻す。
    for (date, bdc, expected) in [
        ("2027-01-31", Bdc::Unadjusted, "2027-01-31"),
        ("2027-01-31", Bdc::Following, "2027-02-01"),
        ("2027-01-31", Bdc::ModifiedFollowing, "2027-01-29"),
        ("2027-05-01", Bdc::Preceding, "2027-04-30"),
        ("2027-05-01", Bdc::ModifiedPreceding, "2027-05-03"),
        ("2027-01-31", Bdc::Nearest, "2027-02-01"),
        ("2027-07-14", Bdc::Nearest, "2027-07-13"),
    ] {
        let (actual, warnings) = adjust(date.parse()?, bdc, &cal)?;
        assert_eq!(actual, expected.parse()?);
        assert_eq!(warnings, []);
        assert_eq!(adjust(actual, bdc, &cal)?.0, actual);
    }
    let rule = AdvanceRule {
        end_of_month: EndOfMonthRule::Ignore,
        bdc: Bdc::Unadjusted,
        calendar: &cal,
    };
    let plus = |date, offset: &str| advance(date, offset.parse()?, &rule).map(|(date, _)| date);
    let friday = "2027-01-29".parse()?;
    let month_end = "2027-01-31".parse()?;
    // 暦日と週は日数で足すので、分けて足しても一度に足しても同じ日になる。
    assert_eq!(plus(plus(friday, "1W")?, "2D")?, plus(friday, "9D")?);
    // 月の加算は結合的でなく、暦日と営業日を混ぜる場合は順序でも結果が違う。
    assert_eq!(plus(plus(month_end, "1M")?, "1M")?, "2027-03-28".parse()?);
    assert_eq!(plus(month_end, "2M")?, "2027-03-31".parse()?);
    assert_eq!(plus(plus(friday, "1D")?, "1BD")?, "2027-02-01".parse()?);
    assert_eq!(plus(plus(friday, "1BD")?, "1D")?, "2027-02-02".parse()?);
    assert_eq!(plus(plus(month_end, "1BD")?, "-1BD")?, friday);
    assert_eq!(plus(month_end, "0D")?, month_end);
    assert_eq!(plus(month_end, "0BD")?, "2027-02-01".parse()?);
    let modified = AdvanceRule {
        bdc: Bdc::Preceding,
        ..rule
    };
    assert_eq!(advance(month_end, "0D".parse()?, &modified)?.0, friday);
    assert_eq!(
        advance(month_end, "0BD".parse()?, &modified)?.0,
        "2027-02-01".parse()?
    );
    let last_biz = AdvanceRule {
        end_of_month: EndOfMonthRule::LastBusinessDay,
        ..rule
    };
    assert_eq!(
        advance(friday, "1M".parse()?, &last_biz)?.0,
        "2027-02-26".parse()?
    );
    let last_day = AdvanceRule {
        end_of_month: EndOfMonthRule::LastCalendarDay,
        ..rule
    };
    let february_end = "2027-02-28".parse()?;
    assert_eq!(
        advance(february_end, "1Y".parse()?, &last_day)?.0,
        "2028-02-29".parse()?
    );
    assert!(matches!(
        plus("9999-12-31".parse()?, "1M"),
        Err(Error::InvalidInput(_))
    ));

    // 金曜から 2BD のスポットは火曜。ON、TN、SN と長さの起点を区別する。
    for (tenor, start, end) in [
        (Tenor::Overnight, "2027-01-29", "2027-02-01"),
        (Tenor::TomorrowNext, "2027-02-01", "2027-02-02"),
        (Tenor::SpotNext, "2027-02-02", "2027-02-03"),
        ("1W".parse()?, "2027-02-02", "2027-02-09"),
    ] {
        assert_eq!(
            money_market_period(friday, tenor, "2BD".parse()?, &rule)?,
            (DateInterval::new(start.parse()?, end.parse()?)?, vec![])
        );
    }
    Ok(())
}

#[test]
fn adjustments_report_each_extrapolated_date_visited() -> Result<()> {
    let cal = calendar(
        "JPTO",
        base("2027-01-29", "2027-01-29", &[], Extrapolation::WeekendsOnly)?,
    )?;
    let saturday = "2027-01-30".parse()?;
    let sunday = "2027-01-31".parse()?;
    let monday = "2027-02-01".parse()?;
    let warn = |date| -> Result<_> {
        Ok(Warning::CalendarExtrapolated {
            calendar: "JPTO".parse()?,
            date,
        })
    };
    // ModifiedFollowing は月末の日曜から同じ月の金曜へ戻り、調べた土日を報告する。
    assert_eq!(
        adjust(sunday, Bdc::ModifiedFollowing, &cal)?,
        ("2027-01-29".parse()?, vec![warn(saturday)?, warn(sunday)?])
    );
    let rule = AdvanceRule {
        end_of_month: EndOfMonthRule::Ignore,
        bdc: Bdc::Unadjusted,
        calendar: &cal,
    };
    assert_eq!(
        advance(sunday, "0BD".parse()?, &rule)?,
        (monday, vec![warn(sunday)?, warn(monday)?])
    );
    Ok(())
}

#[test]
fn nearest_moves_sunday_and_monday_forward_and_other_days_back() -> Result<()> {
    // Nearest は距離でなく曜日で向きを決める（ISDA の Nearest）。
    // 2027-07-19 は月曜の祝日（海の日）とする。
    let cal = calendar(
        "JPTO",
        base(
            "2027-01-01",
            "2027-12-31",
            &["2027-07-19"],
            Extrapolation::Strict,
        )?,
    )?;
    for (date, expected, why) in [
        // 月曜の祝日は次の営業日（火曜）へ進む。
        ("2027-07-19", "2027-07-20", "Monday holiday moves forward"),
        // 日曜も次へ進み、翌日の月曜が祝日なら火曜まで進む。
        (
            "2027-07-18",
            "2027-07-20",
            "Sunday skips the Monday holiday",
        ),
        // 土曜は前の営業日（金曜）へ戻る。
        ("2027-01-30", "2027-01-29", "Saturday moves back"),
        (
            "2027-07-17",
            "2027-07-16",
            "Saturday before a holiday moves back",
        ),
    ] {
        let (actual, warnings) = adjust(date.parse()?, Bdc::Nearest, &cal)?;
        assert_eq!(actual, expected.parse::<Date>()?, "{why}");
        assert_eq!(warnings, []);
    }
    Ok(())
}
`````

`````
crates/rhodium-core/tests/compile_fail.rs
//! 日付と瞬間の違い、および日付の加算に規則が必要なことを確かめる。

use rhodium_testkit::compile_fail_pair;

#[test]
fn date_plus_tenor_requires_advance() {
    // 日付のずれは暦日数か、営業日規則を明示する advance で組み立てる。
    compile_fail_pair!("date_plus_tenor");
}

#[test]
fn model_time_requires_instants() {
    // モデル時間の計算には UTC の瞬間を渡し、現地の日付を取り違えない。
    compile_fail_pair!("date_as_instant");
}

#[test]
fn rate_kind_mismatch() {
    // 複利型が違う率は足せない。
    compile_fail_pair!("rate_kind_mismatch");
}

#[test]
fn rate_times_year_fraction() {
    // 率と日数係数は直接掛けられない。
    compile_fail_pair!("rate_times_year_fraction");
}

#[test]
fn survival_as_discount() {
    // 生存確率は割引因子として受け付けない。
    compile_fail_pair!("survival_as_discount");
}

#[test]
fn normal_plus_lognormal() {
    // 原資産の単位を持つ正規ボラと対数正規ボラは比較できない。
    compile_fail_pair!("normal_plus_lognormal");
}

#[test]
fn clean_price_settlement() {
    // 経過利息を含める前の価格では決済額を求められない。
    compile_fail_pair!("clean_price_settlement");
}

#[test]
fn fx_spot_as_today() {
    // スポット受渡の FX レートでは現在価値を換算できない。
    compile_fail_pair!("fx_spot_as_today");
}

#[test]
fn upfront_plus_running() {
    // アップフロントとランニング・スプレッドは足せない。
    compile_fail_pair!("upfront_plus_running");
}

#[test]
fn futures_price_as_rate() {
    // 先物価格はそのまま利息計算用の率として使えない。
    compile_fail_pair!("futures_price_as_rate");
}
`````

`````
crates/rhodium-core/tests/day_count.rs
//! 日数係数を独立に数えた分子・分母と照合する。

use std::collections::BTreeMap;

use rhodium_core::{
    Error, Result, Warning,
    calendar::{BaseCalendar, Calendar, CalendarSpec, Extrapolation, WeekendRule},
    date::{DateInterval, Weekday},
    day_count::{AccrualContext, DayCount, MonthEndSchedule, PeriodFrame, year_fraction},
    id::CalendarId,
    quantity::YearFraction,
    tenor::Frequency,
};

fn period(start: &str, end: &str) -> Result<DateInterval> {
    DateInterval::new(start.parse()?, end.parse()?)
}

fn calendar(id: &str, extrapolation: Extrapolation) -> Result<Calendar> {
    let id: CalendarId = id.parse()?;
    let base = BaseCalendar::new(
        "2028-02-01".parse()?..="2028-02-29".parse()?,
        WeekendRule::new([Weekday::Sat, Weekday::Sun]),
        vec!["2028-02-29".parse()?],
        extrapolation,
    )?;
    Calendar::compile(
        &CalendarSpec::Base(id.clone()),
        &BTreeMap::from([(id, base)]),
    )
}

fn assert_close(actual: f64, expected: f64, dc: &DayCount) {
    let scale = actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= (1e-12 * scale).max(1e-14),
        "{dc}: {actual} != {expected}"
    );
}

#[test]
#[expect(
    clippy::too_many_lines,
    reason = "15 種の字句と独立した手計算例を一つの表で照合する"
)]
fn all_day_counts_match_hand_calculations() -> Result<()> {
    let cal = calendar("JPTO", Extrapolation::Strict)?;
    let leap = period("2028-02-28", "2028-03-01")?;
    let february_to_august = period("2027-02-28", "2027-08-31")?;
    let st = period("2026-03-01", "2026-06-15")?;
    // 枠の通常期間は頻度のおよそ一刻みでなければならないので、閏日をまたぐ 2 日の区間は
    // それを含む年払いの通常期間で覆う。ACT/365L が読むのは利払いの区間と頻度だけ。
    let frames = [
        PeriodFrame::new(
            leap,
            vec![period("2027-03-01", "2028-03-01")?],
            Frequency::ANNUAL,
        )?,
        PeriodFrame::new(
            st,
            vec![period("2025-06-15", "2026-06-15")?],
            Frequency::ANNUAL,
        )?,
    ];
    // 閏日をまたぐ実日数は 2 日。NL/365 と、閏日を休日にした BUS/252 は 1 日を数える。
    // 年をまたぐ ISDA は 47 日と 135 日を各年の長さで割り、ICMA の短いスタブは 106 日を
    // 仮想の年払い期間 365 日で割る。30/360 の例は日数基準ごとの月末置き換えを反映する。
    let cases = [
        (DayCount::Act360, "ACT360", leap, 2.0 / 360.0, None),
        (DayCount::Act365F, "ACT365F", leap, 2.0 / 365.0, None),
        (
            DayCount::Act365L,
            "ACT365L",
            leap,
            2.0 / 366.0,
            Some(&frames[0]),
        ),
        (DayCount::Act364, "ACT364", leap, 2.0 / 364.0, None),
        (
            DayCount::ActActIsda,
            "ACTACT_ISDA",
            period("2027-11-15", "2028-05-15")?,
            47.0 / 365.0 + 135.0 / 366.0,
            None,
        ),
        (
            DayCount::ActActIcma,
            "ACTACT_ICMA",
            st,
            106.0 / 365.0,
            Some(&frames[1]),
        ),
        (
            DayCount::ActActAfb,
            "ACTACT_AFB",
            period("2027-02-27", "2028-02-29")?,
            1.0 + 1.0 / 365.0,
            None,
        ),
        (DayCount::Nl365, "NL365", leap, 1.0 / 365.0, None),
        (
            DayCount::Thirty360BondBasis,
            "THIRTY_360_BOND",
            february_to_august,
            183.0 / 360.0,
            None,
        ),
        (
            DayCount::Thirty360Us,
            "THIRTY_360_US",
            february_to_august,
            180.0 / 360.0,
            None,
        ),
        (
            DayCount::ThirtyE360,
            "THIRTY_E_360",
            february_to_august,
            182.0 / 360.0,
            None,
        ),
        (
            DayCount::ThirtyE360Isda,
            "THIRTY_E_360_ISDA",
            period("2027-08-31", "2028-02-29")?,
            179.0 / 360.0,
            None,
        ),
        (
            DayCount::ThirtyEPlus360,
            "THIRTY_E_PLUS_360",
            period("2027-01-31", "2027-03-31")?,
            61.0 / 360.0,
            None,
        ),
        (
            DayCount::Bus252 {
                calendar: "JPTO".parse()?,
            },
            "BUS252(JPTO)",
            leap,
            1.0 / 252.0,
            None,
        ),
        (DayCount::OneOne, "ONE_ONE", leap, 1.0, None),
    ];
    for (dc, text, accrual, expected, frame) in cases {
        assert_eq!(text.parse::<DayCount>()?, dc);
        assert_eq!(dc.to_string(), text);
        let ctx = AccrualContext {
            frame,
            month_end: MonthEndSchedule::Aligned,
            maturity: Some("2028-02-29".parse()?),
            calendar: Some(&cal),
        };
        let (fraction, warnings) = year_fraction(&dc, accrual.start(), accrual.end(), &ctx)?;
        assert_eq!(fraction.basis(), &dc);
        assert_close(fraction.value(), expected, &dc);
        assert_eq!(warnings, []);

        // 空区間では暦、満期、頻度が渡されていなくても 0 で、タグを保つ。
        let (empty, warnings) = year_fraction(
            &dc,
            accrual.start(),
            accrual.start(),
            &AccrualContext::default(),
        )?;
        assert_close(empty.value(), 0.0, &dc);
        assert_eq!(empty.basis(), &dc);
        assert_eq!(warnings, []);
    }
    Ok(())
}

#[test]
fn day_count_fractions_agree_with_published_oracles() -> Result<()> {
    // Strata の公開検算（固定コミット 910a144fc49532c1411ed95137556476d622c85c）の期間と値。
    // BUS/252 だけは QuantLib v1.40 の最初の例で、金曜から月曜までの区間で金曜だけを数える。
    let brsp: CalendarId = "BRSP".parse()?;
    let bus_calendar = Calendar::compile(
        &CalendarSpec::Base(brsp.clone()),
        &BTreeMap::from([(
            brsp.clone(),
            BaseCalendar::new(
                "2002-02-01".parse()?..="2002-02-03".parse()?,
                WeekendRule::new([Weekday::Sat, Weekday::Sun]),
                vec![],
                Extrapolation::Strict,
            )?,
        )]),
    )?;
    // 標準の区間は 62 日で、終了年の 2012 年は閏年。30 日系の区間は 1 年 3 か月で両端が 31 日。
    let standard = period("2011-12-28", "2012-02-28")?;
    let thirty_days = period("2012-05-31", "2013-08-31")?;
    // ACT/ACT ICMA の長いスタブは、91 日ずつの四半期の通常期間二つで覆う。
    let icma_frame = PeriodFrame::new(
        period("2011-10-01", "2012-02-29")?,
        vec![
            period("2011-08-31", "2011-11-30")?,
            period("2011-11-30", "2012-02-29")?,
        ],
        "3M".parse()?,
    )?;
    // 期待値の内訳:
    // - ISDA は 2011 年の 4 日と 2012 年の 60 日を各年の長さで割る。
    // - AFB は 1997-06-30 から 1994-06-30 まで 3 年さかのぼり、残り 140 日を 365 で割る。
    // - NL/365 は 63 日から終了日の 2 月 29 日を除く。
    // - 30/360 US は月末の日程で 2 月末の開始日を 30 日にする。
    // - 30E/360 ISDA は満期日の 2 月 29 日を置き換えない。
    // - 30E+/360 は 8 月 31 日を 9 月 1 日にする。
    let cases = [
        (DayCount::Act360, standard, 62.0 / 360.0),
        (DayCount::Act365F, standard, 62.0 / 365.0),
        (DayCount::Act365L, standard, 62.0 / 366.0),
        (
            DayCount::Act364,
            period("2012-02-29", "2012-03-28")?,
            28.0 / 364.0,
        ),
        (
            DayCount::ActActIsda,
            period("2011-12-28", "2012-03-01")?,
            4.0 / 365.0 + 60.0 / 366.0,
        ),
        (
            DayCount::ActActIcma,
            period("2011-10-01", "2012-01-12")?,
            60.0 / 364.0 + 43.0 / 364.0,
        ),
        (
            DayCount::ActActAfb,
            period("1994-02-10", "1997-06-30")?,
            3.0 + 140.0 / 365.0,
        ),
        (
            DayCount::Nl365,
            period("2011-12-28", "2012-02-29")?,
            62.0 / 365.0,
        ),
        (DayCount::Thirty360BondBasis, thirty_days, 450.0 / 360.0),
        (
            DayCount::Thirty360Us,
            period("2012-02-29", "2012-03-28")?,
            28.0 / 360.0,
        ),
        (DayCount::ThirtyE360, thirty_days, 450.0 / 360.0),
        (
            DayCount::ThirtyE360Isda,
            period("2011-02-28", "2012-02-29")?,
            359.0 / 360.0,
        ),
        (DayCount::ThirtyEPlus360, thirty_days, 451.0 / 360.0),
        (
            DayCount::Bus252 { calendar: brsp },
            period("2002-02-01", "2002-02-04")?,
            1.0 / 252.0,
        ),
        (DayCount::OneOne, standard, 1.0),
    ];
    // ACT/365L の 62 日の区間は半年払いの短いスタブとし、満期から 6 か月さかのぼった
    // 仮想の通常期間を参照枠にする（終了年の 2012 年が閏年なので分母は 366）。
    let short_stub = PeriodFrame::new(
        standard,
        vec![period("2011-08-28", "2012-02-28")?],
        "6M".parse()?,
    )?;
    for (dc, accrual, expected) in cases {
        // 枠を読むのは ACT/ACT ICMA と ACT/365L だけ。
        let frame = match dc {
            DayCount::ActActIcma => Some(&icma_frame),
            DayCount::Act365L => Some(&short_stub),
            _ => None,
        };
        let ctx = AccrualContext {
            frame,
            month_end: MonthEndSchedule::Aligned,
            maturity: Some("2012-02-29".parse()?),
            calendar: Some(&bus_calendar),
        };
        let (fraction, warnings) = year_fraction(&dc, accrual.start(), accrual.end(), &ctx)?;
        assert_close(fraction.value(), expected, &dc);
        assert_eq!(warnings, []);
    }
    Ok(())
}

#[test]
fn documented_day_count_examples_hold() -> Result<()> {
    // 設計書「日数計算」の計算の例の表を、期間・文脈・分数のまま単体テストにする。
    let default = AccrualContext::default();
    let icma_frame = PeriodFrame::new(
        period("2026-03-01", "2026-06-15")?,
        vec![period("2025-06-15", "2026-06-15")?],
        Frequency::ANNUAL,
    )?;
    let month_end = AccrualContext {
        month_end: MonthEndSchedule::Aligned,
        ..default
    };
    let february_maturity = AccrualContext {
        maturity: Some("2028-02-29".parse()?),
        ..default
    };
    let icma = AccrualContext {
        frame: Some(&icma_frame),
        ..default
    };
    for (dc, accrual, ctx, expected) in [
        (
            DayCount::Thirty360BondBasis,
            period("2027-01-31", "2027-03-31")?,
            default,
            60.0 / 360.0,
        ),
        (
            DayCount::ThirtyEPlus360,
            period("2027-01-31", "2027-03-31")?,
            default,
            61.0 / 360.0,
        ),
        (
            DayCount::Thirty360Us,
            period("2027-02-28", "2027-08-31")?,
            month_end,
            180.0 / 360.0,
        ),
        (
            DayCount::Thirty360BondBasis,
            period("2027-02-28", "2027-08-31")?,
            default,
            183.0 / 360.0,
        ),
        (
            DayCount::ThirtyE360Isda,
            period("2027-08-31", "2028-02-29")?,
            february_maturity,
            179.0 / 360.0,
        ),
        (
            DayCount::ActActIsda,
            period("2027-11-15", "2028-05-15")?,
            default,
            47.0 / 365.0 + 135.0 / 366.0,
        ),
        (
            DayCount::ActActIcma,
            period("2026-03-01", "2026-06-15")?,
            icma,
            106.0 / 365.0,
        ),
        (
            DayCount::Act365F,
            period("2027-04-30", "2027-05-07")?,
            default,
            7.0 / 365.0,
        ),
    ] {
        let (fraction, warnings) = year_fraction(&dc, accrual.start(), accrual.end(), &ctx)?;
        assert_close(fraction.value(), expected, &dc);
        assert_eq!(warnings, []);
    }
    Ok(())
}

#[test]
fn leap_day_boundaries_and_reference_periods_follow_each_basis() -> Result<()> {
    // NL/365 は終了日を含めて閏日を除く。AFB の余りは開始日を含めて閏日を探す。
    for (dc, start, end, expected) in [
        (DayCount::Nl365, "2028-02-28", "2028-02-29", 0.0),
        (DayCount::Nl365, "2028-02-29", "2028-03-01", 1.0 / 365.0),
        (
            DayCount::ActActAfb,
            "2028-02-01",
            "2028-02-29",
            28.0 / 365.0,
        ),
        (DayCount::ActActAfb, "2028-02-29", "2028-03-01", 1.0 / 366.0),
        (
            DayCount::ActActAfb,
            "2028-02-28",
            "2029-02-28",
            1.0 + 1.0 / 365.0,
        ),
        (DayCount::ActActIsda, "2027-01-01", "2030-01-01", 3.0),
    ] {
        let actual = year_fraction(
            &dc,
            start.parse()?,
            end.parse()?,
            &AccrualContext::default(),
        )?
        .0;
        assert_close(actual.value(), expected, &dc);
    }

    // 年払いは利払い区間が閏日を終端として含む場合に 366。年に複数回なら終了年が閏年かで決める。
    // 半年払いの 62 日の区間は短いスタブで、6 か月の仮想の通常期間を参照枠にする。
    for (start, end, frequency, reference, sampled, expected) in [
        (
            "2027-03-01",
            "2028-02-29",
            Frequency::ANNUAL,
            None,
            period("2028-02-28", "2028-02-29")?,
            1.0 / 366.0,
        ),
        (
            "2028-02-29",
            "2029-02-28",
            Frequency::ANNUAL,
            None,
            period("2028-02-29", "2028-03-01")?,
            1.0 / 365.0,
        ),
        (
            "2011-12-28",
            "2012-02-28",
            "6M".parse()?,
            Some(period("2011-08-28", "2012-02-28")?),
            period("2011-12-28", "2012-02-28")?,
            62.0 / 366.0,
        ),
    ] {
        let accrual = period(start, end)?;
        let frame = PeriodFrame::new(accrual, vec![reference.unwrap_or(accrual)], frequency)?;
        let ctx = AccrualContext {
            frame: Some(&frame),
            ..AccrualContext::default()
        };
        let dc = DayCount::Act365L;
        assert_close(
            year_fraction(&dc, sampled.start(), sampled.end(), &ctx)?
                .0
                .value(),
            expected,
            &dc,
        );
    }

    // 長いスタブは 91 日ある四半期枠に 60 日と 43 日だけ重なる。
    let accrued = period("2011-10-01", "2012-01-12")?;
    let frame = PeriodFrame::new(
        period("2011-10-01", "2012-02-29")?,
        vec![
            period("2011-08-31", "2011-11-30")?,
            period("2011-11-30", "2012-02-29")?,
        ],
        "3M".parse()?,
    )?;
    let dc = DayCount::ActActIcma;
    let ctx = AccrualContext {
        frame: Some(&frame),
        ..AccrualContext::default()
    };
    let actual = year_fraction(&dc, accrued.start(), accrued.end(), &ctx)?.0;
    assert_close(actual.value(), 60.0 / 364.0 + 43.0 / 364.0, &dc);
    Ok(())
}

#[test]
fn day_count_context_errors_and_calendar_warnings_propagate() -> Result<()> {
    let start = "2028-02-28".parse()?;
    let end = "2028-03-02".parse()?;
    let dc_bus = "BUS252(JPTO)".parse::<DayCount>()?;
    let default = AccrualContext::default();
    for dc in [
        &DayCount::Act365L,
        &DayCount::ActActIcma,
        &DayCount::ThirtyE360Isda,
        &dc_bus,
    ] {
        assert!(matches!(
            year_fraction(dc, start, end, &default),
            Err(Error::InvalidInput(_))
        ));
    }
    assert!(matches!(
        year_fraction(&DayCount::Act360, end, start, &default),
        Err(Error::InvalidInput(_))
    ));

    let first = period("2028-02-01", "2028-02-15")?;
    let second = period("2028-02-16", "2028-02-29")?;
    let regular = period("2028-02-01", "2028-02-29")?;
    assert!(PeriodFrame::new(regular, vec![first, second], Frequency::ANNUAL).is_err());
    assert!(PeriodFrame::new(regular, vec![first], Frequency::ANNUAL).is_err());
    let weekly = PeriodFrame::new(regular, vec![regular], "1W".parse()?)?;
    let yearly = PeriodFrame::new(
        regular,
        vec![period("2027-03-01", "2028-03-01")?],
        Frequency::ANNUAL,
    )?;
    let dc = DayCount::ActActIcma;
    let ctx = AccrualContext {
        frame: Some(&weekly),
        ..default
    };
    assert!(matches!(
        year_fraction(&dc, regular.start(), regular.end(), &ctx),
        Err(Error::InvalidInput(_))
    ));
    let ctx = AccrualContext {
        frame: Some(&yearly),
        ..default
    };
    assert!(matches!(
        year_fraction(&dc, start, end, &ctx),
        Err(Error::InvalidInput(_))
    ));

    let unrelated = calendar("USNY", Extrapolation::Strict)?;
    let strict = calendar("JPTO", Extrapolation::Strict)?;
    let extra = calendar("JPTO", Extrapolation::WeekendsOnly)?;
    let ctx = |calendar| AccrualContext {
        calendar: Some(calendar),
        ..default
    };
    assert!(matches!(
        year_fraction(&dc_bus, start, end, &ctx(&unrelated)),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        year_fraction(&dc_bus, start, end, &ctx(&strict)),
        Err(Error::MissingData(_))
    ));
    // 月曜と水曜の 2 日を数える。火曜は休日、水曜だけが外挿される。
    let (actual, warnings) = year_fraction(&dc_bus, start, end, &ctx(&extra))?;
    assert_close(actual.value(), 2.0 / 252.0, &dc_bus);
    assert_eq!(
        warnings,
        [Warning::CalendarExtrapolated {
            calendar: "JPTO".parse()?,
            date: "2028-03-01".parse()?
        }]
    );
    assert!("BUS252()".parse::<DayCount>().is_err());
    assert!(YearFraction::new(f64::NAN, DayCount::Act360).is_err());
    Ok(())
}

#[test]
fn thirty_day_counts_apply_february_maturity_and_december_rollover() -> Result<()> {
    let first = "2027-02-28".parse()?;
    let last = "2028-02-29".parse()?;
    // US 基準は EOM の両端を 30 日にする。EOM でなければ元の日付により +1 日。
    let aligned = AccrualContext {
        month_end: MonthEndSchedule::Aligned,
        ..AccrualContext::default()
    };
    for (ctx, expected) in [(aligned, 1.0), (AccrualContext::default(), 361.0 / 360.0)] {
        let dc = DayCount::Thirty360Us;
        assert_close(
            year_fraction(&dc, first, last, &ctx)?.0.value(),
            expected,
            &dc,
        );
    }
    // 2 月の満期だけ元の日を使い、それ以外の月末は 30 日に置き換える。
    let dc = DayCount::ThirtyE360Isda;
    for (maturity, expected) in [(last, 359.0 / 360.0), ("2028-08-31".parse()?, 1.0)] {
        let ctx = AccrualContext {
            maturity: Some(maturity),
            ..AccrualContext::default()
        };
        assert_close(
            year_fraction(&dc, first, last, &ctx)?.0.value(),
            expected,
            &dc,
        );
    }
    // 12 月 31 日は翌年 1 月 1 日に置き換えるので 30 + 1 日。
    let dc = DayCount::ThirtyEPlus360;
    let actual = year_fraction(
        &dc,
        "2027-11-30".parse()?,
        "2027-12-31".parse()?,
        &AccrualContext::default(),
    )?
    .0;
    assert_close(actual.value(), 31.0 / 360.0, &dc);
    Ok(())
}
`````

`````
crates/rhodium-core/tests/error.rs
//! 共有エラーの分類、表示、文脈の付与を確認する。

use std::{mem::discriminant, num::NonZeroU32};

use rhodium_core::{Error, error::Span};

#[test]
fn error_context_preserves_variant_and_prefixes_message() {
    // 呼び出し側がエラーの種類で分岐でき、表示から失敗した処理を特定できる。
    let cases = [
        (Error::InvalidInput("detail".into()), "invalid input"),
        (Error::InvalidData("detail".into()), "invalid data"),
        (Error::MissingData("detail".into()), "missing data"),
        (Error::Calibration("detail".into()), "calibration failed"),
        (Error::Numerical("detail".into()), "numerical failure"),
        (Error::Unsupported("detail".into()), "unsupported"),
        (Error::Execution("detail".into()), "execution failed"),
        (Error::Internal("detail".into()), "internal error"),
    ];

    for (error, label) in cases {
        let category = discriminant(&error);
        let error = error.context("pricing").context("trade T-001");
        assert_eq!(discriminant(&error), category);
        assert_eq!(
            error.to_string(),
            format!("{label}: trade T-001: pricing: detail")
        );
    }
}

#[test]
fn dsl_error_context_preserves_span_and_help() {
    // 取引の文脈を付けても、DSL を修正する人に位置と提案を渡せる。
    let span = Span {
        file: "swap.rh".into(),
        line: NonZeroU32::new(3).expect("positive line"),
        column: NonZeroU32::new(7).expect("positive column"),
    };
    let error = Error::Dsl {
        span: span.clone(),
        message: "unknown index".into(),
        help: Some("use TONA".into()),
    }
    .context("trade T-001");

    assert_eq!(error.to_string(), "trade T-001: unknown index");
    assert!(matches!(
        error,
        Error::Dsl { span: actual_span, help, .. }
            if actual_span == span && help.as_deref() == Some("use TONA")
    ));
}
`````

`````
crates/rhodium-core/tests/money_and_quantities.rs
//! 金額、量、FX、丸めの手計算と境界を確かめる。

use std::{cmp::Ordering, collections::BTreeSet};

use proptest::{
    prelude::*,
    test_runner::{Config, RngSeed},
};
use rhodium_core::{
    Error, Result,
    currency::Currency,
    date::Date,
    day_count::DayCount::{Act360, Act365F},
    fx::{CcyPair, FxForwardPoints, FxRate},
    id::{CalendarId, EntityId, IndexId},
    money::{Money, MoneyBag},
    quantity::{
        Clean, Compounded, Continuous, Correlation, DiscountFactor, ParPrice, Probability, Rate,
        RateFuturesPrice, RecoveryRate, RiskyAnnuity, RunningSpread, ShiftedLognormalVol, Simple,
        Spread, SurvivalProbability, YearFraction, accrue,
    },
    rounding::{Rounding, RoundingMode},
};
use rhodium_testkit::assert_close;
use scioxide_ad::{Dual, Lanes};

#[test]
fn money_currency_mismatch() {
    // 通貨が同じ額だけを加え、異なる通貨はエラーに両方のコードを示す。
    let yen = Money::new(Currency::JPY, 3.0);
    let dollar = Money::new(Currency::USD, 2.0);
    assert!(
        matches!(yen.checked_add(dollar), Err(Error::InvalidInput(message))
        if message.contains("JPY") && message.contains("USD"))
    );
}

#[test]
fn zero_money_bag_terms_keep_sensitivities() -> Result<()> {
    // 0 円に相殺されても、変数の微分 1 は残る。USD への換算後も 1/150。
    let mut bag = MoneyBag::new();
    bag.add(Money::new(Currency::JPY, Dual::variable(5.0)));
    bag.add(Money::new(Currency::JPY, Dual::new(-5.0, 0.0)));
    assert_eq!(bag.terms().count(), 1);
    let yen = bag.get(Currency::JPY);
    assert_close!(yen.value(), 0.0, 1e-12, 1e-6);
    assert_close!(yen.derivative(), 1.0, 1e-10, 1e-6);
    let rate = FxRate::today("USDJPY".parse()?, Dual::new(150.0, 0.0))?.invert()?;
    let converted = bag.convert(Currency::USD, &[rate])?.amount();
    assert_close!(converted.value(), 0.0, 1e-12, 1e-8);
    assert_close!(converted.derivative(), 1.0 / 150.0, 1e-10, 1e-8);

    // 値の等値では、ゼロ項の有無を区別しない。
    let mut zero = MoneyBag::new();
    zero.add(Money::new(Currency::JPY, 0.0));
    assert_eq!(zero, MoneyBag::new());
    Ok(())
}

#[test]
fn accrual_and_growth_follow_compounding_and_day_count() -> Result<()> {
    // 年率 4% を半年適用すると単利と半年複利は 2%、連続複利は exp(0.02) − 1。
    let half = YearFraction::new(0.5, Act360)?;
    let simple = Rate::new(0.04, Act360, Simple)?;
    let continuous = Rate::new(0.04, Act360, Continuous)?;
    let semiannual = Compounded::new("6M".parse()?)?;
    let compounded = Rate::new(0.04, Act360, semiannual)?;
    assert_close!(accrue(&simple, &half)?, 0.02, 1e-12, 1e-14);
    assert_close!(accrue(&compounded, &half)?, 0.02, 1e-12, 1e-14);
    assert_close!(
        accrue(&continuous, &half)?,
        0.020_201_340_026_755_81,
        1e-12,
        1e-14
    );
    assert_close!(
        simple.discount_factor(&half)?.value(),
        1.0 / 1.02,
        1e-12,
        1e-14
    );
    let with_spread = (simple + Spread::new(0.001)?)?;
    assert_close!(with_spread.value(), 0.041, 1e-12, 1e-14);
    let annual = Rate::new(0.04, Act360, Compounded::new("1Y".parse()?)?)?;
    assert!(matches!(compounded + annual, Err(Error::InvalidInput(_))));
    assert!(Compounded::new("5M".parse()?).is_err());
    Ok(())
}

#[test]
fn compounded_factors_keep_small_positive_results() -> Result<()> {
    // exp(-40) は有限で正。1 + expm1(-40) から求めると丸めでゼロになる。
    let year = YearFraction::new(1.0, Act360)?;
    let negative = Rate::new(-40.0, Act360, Continuous)?;
    let positive = Rate::new(40.0, Act360, Continuous)?;
    let expected = 4.248_354_255_291_589e-18;
    assert_close!(negative.growth_factor(&year)?, expected, 1e-12, 0.0);
    assert_close!(
        positive.discount_factor(&year)?.value(),
        expected,
        1e-12,
        0.0
    );
    // 年一回の複利で年率 -100% は対数の定義域を外れる。
    let invalid = Rate::new(-1.0, Act360, Compounded::new("1Y".parse()?)?)?;
    assert!(matches!(accrue(&invalid, &year), Err(Error::Numerical(_))));
    Ok(())
}

#[test]
fn accrue_day_count_mismatch() -> Result<()> {
    // ACT/360 の率と ACT/365F の日数係数を同じ式に渡すと基準の不一致を返す。
    let rate = Rate::new(0.03, Act360, Simple)?;
    let fraction = YearFraction::new(0.5, Act365F)?;
    assert!(
        matches!(accrue(&rate, &fraction), Err(Error::InvalidInput(message))
        if message.contains("ACT360") && message.contains("ACT365F"))
    );
    Ok(())
}

#[test]
fn rate_sum_rejects_a_different_day_count_basis() -> Result<()> {
    // 同じ単利でも ACT/360 と ACT/365F の率は足せない。日数基準は実行時のタグで守る。
    let act360 = Rate::new(0.03, Act360, Simple)?;
    let act365 = Rate::new(0.01, Act365F, Simple)?;
    assert!(matches!(
        act360.clone().checked_sub(&act365),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        act360.clone() + act365,
        Err(Error::InvalidInput(_))
    ));
    // 基準が同じなら和も差も同じ基準の率になる。
    let sum = (act360 + Rate::new(0.01, Act360, Simple)?)?;
    assert_close!(sum.value(), 0.04, 1e-12, 1e-14);
    assert_eq!(sum.basis(), &Act360);
    Ok(())
}

#[test]
fn prices_and_probabilities_use_their_own_operations() -> Result<()> {
    // 額面 1,000、クリーン価格 0.98、経過利息 0.03 の決済額は 1,010 円。
    let price = ParPrice::<Clean>::new(0.98)?.with_accrued(0.03)?;
    assert_close!(
        price
            .settlement(Money::new(Currency::JPY, 1000.0))?
            .amount(),
        1010.0,
        1e-12,
        1e-6
    );
    let future = RateFuturesPrice::new(0.96, Act360)?;
    assert_close!(future.implied_rate()?.value(), 0.04, 1e-12, 1e-14);
    let upfront = RunningSpread::new(0.01)?.to_upfront(RiskyAnnuity::new(3.0)?)?;
    assert_close!(upfront.value(), 0.03, 1e-12, 1e-14);

    // 区間の成長因子は P(s)/P(e)、区間の生存は Q(e)/Q(s)。
    let discount = DiscountFactor::new(0.96)?.growth_to(DiscountFactor::new(0.80)?)?;
    let survival = SurvivalProbability::new(0.80)?.survival_to(SurvivalProbability::new(0.60)?)?;
    assert_close!(discount, 1.2, 1e-12, 1e-14);
    assert_close!(survival.value(), 0.75, 1e-12, 1e-14);
    assert!(DiscountFactor::new(1.02).is_ok());
    assert!(DiscountFactor::new(0.0).is_err());
    assert!(SurvivalProbability::new(1.1).is_err());
    assert!(Probability::new(-0.1).is_err());
    assert!(RecoveryRate::new(f64::NAN).is_err());
    assert!(Correlation::new(-1.1).is_err());
    Ok(())
}

#[test]
fn shift_mismatch() -> Result<()> {
    // 値が同じでもシフトが違えば比較も差も失敗する。
    let a = ShiftedLognormalVol::new(0.20, 0.01)?;
    let b = ShiftedLognormalVol::new(0.20, 0.02)?;
    let c = ShiftedLognormalVol::new(0.21, 0.01)?;
    assert!(matches!(a.compare(b), Err(Error::InvalidInput(_))));
    assert!(matches!(a.difference(b), Err(Error::InvalidInput(_))));
    assert_eq!(a.compare(c)?, Ordering::Less);
    Ok(())
}

#[test]
fn fx_direction() -> Result<()> {
    // USDJPY は USD から JPY だけを換算する。逆向きの換算には反転したレートを渡す。
    let rate = FxRate::today("USDJPY".parse()?, 150.0)?;
    let yen = rate.convert(Money::new(Currency::USD, 2.0))?;
    assert_eq!(yen.ccy(), Currency::JPY);
    assert_close!(yen.amount(), 300.0, 1e-12, 1e-6);
    assert!(matches!(rate.convert(yen), Err(Error::InvalidInput(_))));
    assert_close!(rate.invert()?.convert(yen)?.amount(), 2.0, 1e-12, 1e-8);
    let mut bag = MoneyBag::new();
    bag.add(yen);
    assert!(matches!(
        bag.convert(Currency::USD, &[rate]),
        Err(Error::MissingData(_))
    ));
    let inverse = rate.invert()?;
    assert!(matches!(
        bag.convert(Currency::USD, &[inverse, inverse]),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn fx_cross_spot_dates() -> Result<()> {
    // 同じスポット日なら EURUSD 1.20 と USDJPY 150 の積は EURJPY 180。
    let date = Date::from_ymd(2026, 11, 24)?;
    let next = date.checked_add_days(1)?;
    let eur_usd = FxRate::spot("EURUSD".parse()?, 1.20, date)?;
    let usd_jpy = FxRate::spot("USDJPY".parse()?, 150.0, date)?;
    let crossed = eur_usd.cross(&usd_jpy)?;
    assert_eq!(crossed.pair(), "EURJPY".parse::<CcyPair>()?);
    assert_eq!(crossed.spot_date(), date);
    assert_close!(crossed.rate(), 180.0, 1e-12, 1e-14);
    assert!(matches!(
        eur_usd.cross(&FxRate::spot(usd_jpy.pair(), 150.0, next)?),
        Err(Error::InvalidInput(_))
    ));
    // USDJPY の 95.3 ポイントは 0.953。アウトライトは 150 − 0.953。
    let points = FxForwardPoints::new(usd_jpy.pair(), -95.3, 2, date, next)?;
    assert_close!(points.outright(&usd_jpy)?, 149.047, 1e-12, 1e-14);
    assert!(FxForwardPoints::new(usd_jpy.pair(), 1.0, 2, next, date).is_err());
    Ok(())
}

#[test]
fn fx_today_cross_and_inverse_keep_direction() -> Result<()> {
    // 評価時点の受渡どうしなら、EURUSD 2 と USDJPY 128 の積は EURJPY 256（2 の冪で割り切れる例）。
    let usd_jpy = FxRate::today("USDJPY".parse()?, 128.0)?;
    let eur_usd = FxRate::today("EURUSD".parse()?, 2.0)?;
    let crossed = eur_usd.cross(&usd_jpy)?;
    assert_eq!(crossed.pair(), "EURJPY".parse::<CcyPair>()?);
    assert_close!(crossed.rate(), 256.0, 1e-12, 1e-14);

    // 反転は向きだけを入れ替える。結果も FxRate<Today> なので、そのまま 256 JPY を 2 USD に換算できる。
    let inverse = usd_jpy.invert()?;
    assert_eq!(inverse.pair(), "JPYUSD".parse::<CcyPair>()?);
    assert_close!(inverse.rate(), 0.007_812_5, 1e-12, 1e-14);
    assert_close!(
        inverse.convert(Money::new(Currency::JPY, 256.0))?.amount(),
        2.0,
        1e-12,
        1e-8
    );
    // EURUSD の quote 通貨 USD と JPYUSD の base 通貨 JPY はつながらないので、積を作らない。
    assert!(matches!(
        eur_usd.cross(&inverse),
        Err(Error::InvalidInput(_))
    ));

    // フォワードポイントは、同じペアで自分のスポット日のスポットにだけ足せる。
    // USDJPY の 25 ポイントは pip 指数 2 で 0.25 なので、150 のスポットのアウトライトは 150.25。
    let spot_date = Date::from_ymd(2026, 11, 25)?;
    let points = FxForwardPoints::new(
        usd_jpy.pair(),
        25.0,
        2,
        spot_date,
        Date::from_ymd(2026, 12, 25)?,
    )?;
    let spot = FxRate::spot(usd_jpy.pair(), 150.0, spot_date)?;
    assert_close!(points.outright(&spot)?, 150.25, 1e-12, 1e-14);
    let earlier_spot = FxRate::spot(usd_jpy.pair(), 150.0, spot_date.checked_add_days(-1)?)?;
    let other_pair = FxRate::spot("EURJPY".parse()?, 150.0, spot_date)?;
    assert!(matches!(
        points.outright(&earlier_spot),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        points.outright(&other_pair),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn money_bag_conversion_uses_one_rate_table() -> Result<()> {
    // 4 USD と 2 JPY を USDJPY 128 で JPY に換算すると 4·128 + 2 = 514 JPY。
    let usd = Money::new(Currency::USD, 4.0);
    let yen = Money::new(Currency::JPY, 2.0);
    let rate = FxRate::today("USDJPY".parse()?, 128.0)?;
    let mut bag = MoneyBag::new();
    bag.add(usd);
    bag.add(yen);
    let reported = bag.convert(Currency::JPY, &[rate])?;
    assert_eq!(reported.ccy(), Currency::JPY);
    assert_close!(reported.amount(), 514.0, 1e-12, 1e-6);

    // 袋の換算は、同じ表で項ごとに換算してから足した値と一致する。
    let termwise = rate.convert(usd)?.checked_add(yen)?;
    assert_close!(reported.amount(), termwise.amount(), 1e-12, 1e-6);
    // 袋に無い通貨のレートは、表に含まれていても換算に使わず、結果を変えない。
    let unused = FxRate::today("GBPJPY".parse()?, 200.0)?;
    assert_close!(
        bag.convert(Currency::JPY, &[unused, rate])?.amount(),
        514.0,
        1e-12,
        1e-6
    );

    // 和が f64 の範囲を超える金額は、無限大を値として返さず Numerical にする。
    let largest = Money::new(Currency::USD, f64::MAX);
    assert!(matches!(
        largest.checked_add(largest),
        Err(Error::Numerical(_))
    ));
    Ok(())
}

#[test]
fn rounding_uses_decimal_ties_and_straight_through_lanes() -> Result<()> {
    // 0.07 は十進の格子上にあり、1.005 は 1.00 と 1.01 の中点。
    let ceiling = Rounding::new(0.01, RoundingMode::Ceiling)?;
    assert_close!(ceiling.apply(0.07)?, 0.07, 1e-12, 1e-8);
    let away = Rounding::new(0.01, RoundingMode::HalfAwayFromZero)?;
    let even = Rounding::new(0.01, RoundingMode::HalfEven)?;
    assert_close!(away.apply(1.005)?, 1.01, 1e-12, 1e-8);
    assert_close!(even.apply(1.005)?, 1.00, 1e-12, 1e-8);
    assert_close!(even.apply(-1.015)?, -1.02, 1e-12, 1e-8);
    let rounded = away.apply(Dual::new(
        Lanes::from_array([1.005, -1.005]),
        Lanes::from_array([1.0, 2.0]),
    ))?;
    let [positive, negative] = rounded.value().to_array();
    let [first, second] = rounded.derivative().to_array();
    assert_close!(positive, 1.01, 1e-12, 1e-8);
    assert_close!(negative, -1.01, 1e-12, 1e-8);
    assert_close!(first, 1.0, 1e-10, 1e-8);
    assert_close!(second, 2.0, 1e-10, 1e-8);
    assert!(matches!(
        away.apply(f64::INFINITY),
        Err(Error::Numerical(_))
    ));
    assert!(matches!(away.apply(1e20), Err(Error::Numerical(_))));
    Ok(())
}

// 丸めた値の比較では、許容誤差を増分（ここでは 1 セント）よりはるかに小さくし、
// 隣の格子点への取り違えを必ず失敗にする。

#[test]
fn rounding_modes_round_signed_midpoints() -> Result<()> {
    use RoundingMode::{AwayFromZero, Ceiling, Floor, HalfAwayFromZero, HalfEven, TowardZero};
    // 0.125 は 1 セントの格子点 0.12 と 0.13 のちょうど中点。最近接の二方式は中点の扱いだけが違い、
    // 方向付きの四方式は符号で行き先が分かれる（Floor は負側を、Ceiling は正側を 0 から遠ざける）。
    for (mode, positive, negative) in [
        (HalfAwayFromZero, 0.13, -0.13),
        (HalfEven, 0.12, -0.12),
        (TowardZero, 0.12, -0.12),
        (AwayFromZero, 0.13, -0.13),
        (Floor, 0.12, -0.13),
        (Ceiling, 0.13, -0.12),
    ] {
        let cent = Rounding::new(0.01, mode)?;
        assert_close!(cent.apply(0.125)?, positive, 1e-12, 1e-12);
        assert_close!(cent.apply(-0.125)?, negative, 1e-12, 1e-12);
    }

    // 中点のすぐ下の f64 は最短の十進表記も 0.12499… なので、中点として扱わず切り捨てる。
    let nearest = Rounding::new(0.01, HalfAwayFromZero)?;
    assert_close!(nearest.apply(0.125_f64.next_down())?, 0.12, 1e-12, 1e-12);
    // HalfEven の中点 0.375 は、37 と 38 のうち偶数の 38 セントへ丸める。
    let even = Rounding::new(0.01, HalfEven)?;
    assert_close!(even.apply(0.375)?, 0.38, 1e-12, 1e-12);
    // -0.001 は 0 に丸まり、支払額の表示で -0 にならないよう符号を持たないゼロを返す。
    let zero = nearest.apply(-0.001)?;
    assert_close!(zero, 0.0, 1e-12, 1e-12);
    assert!(zero.is_sign_positive());

    // 入れ子の Dual でも値だけを丸め、微分は丸めない値の微分を通す（straight-through）。
    // 丸めの二次微分は 0 とするので、外側の微分の二次成分は入力の 4 のまま残る。
    let nested = even.apply(Dual::new(Dual::new(0.125, 2.0), Dual::new(3.0, 4.0)))?;
    assert_close!(nested.value().value(), 0.12, 1e-12, 1e-12);
    assert_close!(nested.value().derivative(), 2.0, 1e-10, 1e-12);
    assert_close!(nested.derivative().value(), 3.0, 1e-10, 1e-12);
    assert_close!(nested.derivative().derivative(), 4.0, 1e-8, 1e-12);
    Ok(())
}

#[test]
fn rounding_keeps_grid_values_and_reports_range_limits() -> Result<()> {
    use RoundingMode::{AwayFromZero, Ceiling, Floor, HalfAwayFromZero, HalfEven, TowardZero};
    // 方向付きの方式は格子から外れた値を必ず動かす。f64 の 0.03 のように二進では格子から
    // わずかにずれる値も、最短の十進表記では格子上にあるので動かない。丸めた値の丸め直しも不変。
    for mode in [Floor, Ceiling, AwayFromZero, TowardZero] {
        let cent = Rounding::new(0.01, mode)?;
        for cents in -1000_i32..=1000 {
            let on_grid = f64::from(cents) / 100.0;
            assert_close!(cent.apply(on_grid)?, on_grid, 1e-12, 1e-12);
        }
        for x in [0.071, -0.071, 19.995, 1234.5649] {
            let once = cent.apply(x)?;
            assert_close!(cent.apply(once)?, once, 1e-12, 1e-12);
        }
    }

    // 二進では 2.67499… でも、十進の 2.675 は 2.67 と 2.68 の中点なので 0 から遠い 2.68 にする。
    let nearest = Rounding::new(0.01, HalfAwayFromZero)?;
    assert_close!(nearest.apply(2.675)?, 2.68, 1e-12, 1e-12);

    // 増分 1 で刻みの数がちょうど 2^53 の値は表せるので、同じ格子点（差が半刻み未満）に残る。
    // その次の f64 は刻みの数が 2^53 を超え、整数倍を正確に数えられないため Numerical。
    let unit = Rounding::new(1.0, HalfEven)?;
    let max_ticks = 9_007_199_254_740_992.0_f64;
    assert_close!(unit.apply(max_ticks)?, max_ticks, 0.0, 0.5);
    assert!(matches!(
        unit.apply(max_ticks.next_up()),
        Err(Error::Numerical(_))
    ));

    // 最小のサブノーマルを増分にしても、その値自身は 1 刻み目の格子点で動かない。
    let tiny = f64::from_bits(1);
    let subnormal = Rounding::new(tiny, HalfEven)?;
    assert_close!(subnormal.apply(tiny)? / tiny, 1.0, 1e-12, 0.0);

    // 1e308 刻みで 1.5e308 を切り上げると 2e308 で f64 の範囲を超え、NaN の入力も丸められない。
    let huge = Rounding::new(1.0e308, Ceiling)?;
    assert!(matches!(huge.apply(1.5e308), Err(Error::Numerical(_))));
    assert!(matches!(unit.apply(f64::NAN), Err(Error::Numerical(_))));
    Ok(())
}

#[test]
fn ids_preserve_opaque_labels_and_validate_lei() -> Result<()> {
    // CalendarId は + を含む 1 つの文字列。LEI: の場合は検査桁も照合する。
    let id: CalendarId = "JPTO+USNY".parse()?;
    assert_eq!(id.as_str(), "JPTO+USNY");
    assert!("JPY TONA".parse::<IndexId>().is_err());
    assert!("jpy".parse::<Currency>().is_err());
    assert!("LEI:5493001KJTIIGC8Y1R12".parse::<EntityId>().is_ok());
    assert!("LEI:5493001KJTIIGC8Y1R13".parse::<EntityId>().is_err());
    Ok(())
}

proptest! {
    #![proptest_config(Config { cases: 256, rng_seed: RngSeed::Fixed(158), ..Config::default() })]

    #[test]
    fn money_bag_sums_do_not_depend_on_order(
        rows in prop::collection::vec(
            (prop::sample::select(vec!["JPY", "USD", "EUR"]), -1.0e6..1.0e6_f64),
            0..24,
        ),
    ) {
        // 通貨ごとの和は加える順序に依らず（可換・結合）、袋は加えた通貨の項を一つも落とさない。
        // 浮動小数点の和は順序で丸めが変わるため、相対の尺度に項の絶対値の和 M を使い、
        // 品質の章の PV の許容誤差 max(1e-12·M, 1e-8) で比べる。
        let rows = rows
            .into_iter()
            .map(|(code, amount)| Ok((code.parse::<Currency>()?, amount)))
            .collect::<Result<Vec<_>>>()?;
        let mut forward = MoneyBag::new();
        let mut backward = MoneyBag::new();
        for &(ccy, amount) in &rows {
            forward.add(Money::new(ccy, amount));
        }
        for &(ccy, amount) in rows.iter().rev() {
            backward.add(Money::new(ccy, amount));
        }
        let currencies: BTreeSet<Currency> = rows.iter().map(|&(ccy, _)| ccy).collect();
        prop_assert!(forward.terms().map(Money::ccy).eq(currencies.iter().copied()));

        // 同じ額を符号を変えて足し戻すと値は 0 に戻るが、項は消えない（金額と量の D4）。
        let mut cancelled = forward.clone();
        for &(ccy, amount) in &rows {
            cancelled.add(Money::new(ccy, -amount));
        }
        prop_assert!(cancelled.terms().map(Money::ccy).eq(currencies.iter().copied()));

        for ccy in currencies {
            let (sum, scale) = rows
                .iter()
                .filter(|&&(code, _)| code == ccy)
                .fold((0.0, 0.0), |(sum, scale), &(_, amount)| (sum + amount, scale + amount.abs()));
            let tolerance = (1e-12 * scale).max(1e-8);
            prop_assert!((forward.get(ccy) - backward.get(ccy)).abs() <= tolerance);
            prop_assert!((forward.get(ccy) - sum).abs() <= tolerance);
            prop_assert!(cancelled.get(ccy).abs() <= 2.0 * tolerance);
        }
    }
}
`````

`````
crates/rhodium-core/tests/parallel.rs
//! 対象を区間に分ける並列の処理が、分け方によらず同じ対象を同じ順に覆うことを確かめる。

use std::{
    num::NonZeroUsize,
    sync::atomic::{AtomicUsize, Ordering},
};

use rhodium_core::{
    Error, Result,
    parallel::{ExecutionContext, ThreadExecutor, check, control_failure, map_partitions},
};
use rhodium_testkit::TickClock;
use scioxide_core::{
    CancellationToken, ClockBinding, Deadline, InvalidInputReason, SciOxideError, argument,
};

fn workers(count: usize) -> NonZeroUsize {
    NonZeroUsize::new(count).unwrap_or(NonZeroUsize::MIN)
}

#[test]
fn partitions_cover_every_item_once_in_order() -> Result<()> {
    // 対象の数とワーカー数の組ごとに、区間を順に連結すると 0..len になり、
    // 区間の数は min(ワーカー数, 対象の数)、どの区間も空でないことを確かめる。
    for len in [0, 1, 2, 5, 17, 100] {
        for count in [1, 2, 3, 8] {
            let executor = ThreadExecutor::new(workers(count));
            let ranges = map_partitions(ExecutionContext::new(&executor), len, |range, _| range)?;
            assert_eq!(ranges.len(), count.min(len), "len {len}, workers {count}");
            assert!(ranges.iter().all(|range| !range.is_empty()));
            let items: Vec<_> = ranges.into_iter().flatten().collect();
            assert_eq!(items, (0..len).collect::<Vec<_>>());
        }
    }
    // 逐次の実行文脈は全体を一つの区間で処理する。
    let ranges = map_partitions(ExecutionContext::sequential(), 7, |range, _| {
        (range.start, range.end)
    })?;
    assert_eq!(ranges, vec![(0, 7)]);
    Ok(())
}

#[test]
fn partition_work_runs_once_with_a_single_worker_budget() -> Result<()> {
    // 区間の中の処理は子の実行文脈を受け取り、その予算は 1 なので入れ子の並列にならない。
    let executor = ThreadExecutor::new(workers(4));
    let calls = AtomicUsize::new(0);
    let budgets = map_partitions(ExecutionContext::new(&executor), 10, |_, child| {
        calls.fetch_add(1, Ordering::Relaxed);
        child.worker_budget().get()
    })?;
    assert_eq!(budgets, [1, 1, 1, 1]);
    assert_eq!(calls.load(Ordering::Relaxed), 4);
    Ok(())
}

#[test]
fn stopped_or_panicking_execution_maps_to_rhodium_errors() {
    // 取り消された実行文脈は呼び出し側の設定による停止なので Execution。
    let executor = ThreadExecutor::new(workers(2));
    let token = CancellationToken::new();
    token.cancel();
    let cancelled = ExecutionContext::new(&executor).with_cancellation(&token);
    assert!(matches!(
        map_partitions(cancelled, 4, |range, _| range.len()),
        Err(Error::Execution(message)) if message.contains("cancel")
    ));
    // 区間の処理の panic はスレッドの実行器が捕まえ、rhodium の誤りとして Internal を返す。
    let result = map_partitions(ExecutionContext::new(&executor), 4, |range, _| {
        assert!(range.start == 0, "deliberate panic in the second partition");
        range.len()
    });
    assert!(matches!(result, Err(Error::Internal(message)) if message.contains("panicked")));
}

#[test]
fn empty_input_still_stops_when_cancelled() -> Result<()> {
    // 対象が 0 件だと実行器は区間を作らず、区間の開始前の停止の確認も行わない。
    // 取り消した実行が空の結果で成功すると、取引の無いポートフォリオで取消を見逃すので、
    // map_partitions の入口で確かめる。止まっていなければ空の結果を返す（対の正常例）。
    let executor = ThreadExecutor::new(workers(2));
    let token = CancellationToken::new();
    token.cancel();
    for execution in [
        ExecutionContext::sequential().with_cancellation(&token),
        ExecutionContext::new(&executor).with_cancellation(&token),
    ] {
        assert!(matches!(
            map_partitions(execution, 0, |range, _| range.len()),
            Err(Error::Execution(message)) if message.contains("cancel")
        ));
    }
    assert_eq!(
        map_partitions(ExecutionContext::sequential(), 0, |range, _| range.len())?,
        Vec::<usize>::new()
    );
    Ok(())
}

// 実行文脈の組み立ての失敗はテストの誤りなので、そのまま報告する。
#[expect(clippy::needless_pass_by_value, reason = "map_err の変換関数に使う")]
fn scioxide(error: SciOxideError) -> Error {
    Error::Internal(error.to_string())
}

#[test]
fn deadline_after_partitions_start_stops_at_the_next_item_boundary() -> Result<()> {
    // 実行器は区間の開始前にだけ期限を確かめる。区間の中で対象ごとに check を呼べば、
    // 区間の開始後に来た期限もその対象の境界で Execution として止まる。
    // 時計は読むたびに進むので、逐次なら止まる位置が決まる。開始前の確認が 0 回目、
    // 対象 i の前の check が i + 1 回目の読み出しで、期限を 3 回目に置くと対象 0 と 1 だけを処理する。
    let clock = TickClock::new();
    let binding = ClockBinding::new(&clock);
    let execution = ExecutionContext::sequential()
        .with_deadline(Deadline::at(&binding, TickClock::after(3)))
        .map_err(scioxide)?;
    let processed = AtomicUsize::new(0);
    let results = map_partitions(execution, 10, |range, child| {
        for _ in range {
            check(child)?;
            processed.fetch_add(1, Ordering::Relaxed);
        }
        Ok::<_, Error>(())
    })?;
    assert!(matches!(
        results.as_slice(),
        [Err(Error::Execution(message))] if message.contains("deadline")
    ));
    assert_eq!(processed.load(Ordering::Relaxed), 2);

    // 4 ワーカーでは読み出しの順がスレッドの進み方で変わる。遅れて始まった区間は開始前の確認で
    // 止まり（map_partitions 自体が Execution）、始まった区間は対象の境界で止まる。
    // どちらでも呼び出し全体は Execution になり、処理した対象は全体より少ない。
    let clock = TickClock::new();
    let binding = ClockBinding::new(&clock);
    let executor = ThreadExecutor::new(workers(4));
    let execution = ExecutionContext::new(&executor)
        .with_deadline(Deadline::at(&binding, TickClock::after(6)))
        .map_err(scioxide)?;
    let processed = AtomicUsize::new(0);
    let outcome = map_partitions(execution, 40, |range, child| {
        for _ in range {
            check(child)?;
            processed.fetch_add(1, Ordering::Relaxed);
        }
        Ok::<_, Error>(())
    })
    .and_then(|results| results.into_iter().collect::<Result<Vec<_>>>());
    assert!(matches!(outcome, Err(Error::Execution(message)) if message.contains("deadline")));
    assert!(processed.load(Ordering::Relaxed) < 40);
    Ok(())
}

#[test]
fn control_failure_separates_stops_from_numerical_errors() {
    // 取消で止まった実行のエラーは呼び出し全体の Execution、計算の入力や数値の失敗は None で、
    // 呼び出し側（AD の勾配など）が自分の文脈のエラーにする。
    let token = CancellationToken::new();
    token.cancel();
    let stopped = ExecutionContext::sequential()
        .with_cancellation(&token)
        .check()
        .err();
    assert!(matches!(
        stopped.as_ref().and_then(control_failure),
        Some(Error::Execution(message)) if message.contains("cancel")
    ));
    let invalid = SciOxideError::invalid_input(
        argument!("theta"),
        InvalidInputReason::Inconsistent,
        "valuation failed",
    );
    assert!(control_failure(&invalid).is_none());
}
`````

`````
crates/rhodium-core/tests/refdata_and_fixing.rs
//! 参照データと fixing の公表時刻・評価時点の境界を確かめる。

use std::collections::BTreeMap;

use rhodium_core::{
    Error, Result, Warning,
    calendar::{BaseCalendar, CalendarSpec, Extrapolation, WeekendRule},
    currency::{Asset, Currency},
    date::{Date, DateInterval, Weekday},
    day_count::DayCount,
    fixing::{
        CreditEventRecord, FixingResolution, FixingStore, FixingStoreParts, ObservationPeriod,
        ObservationPoint, YearMonth,
    },
    fx::CcyPair,
    id::{CalendarId, CommodityId, CreditKey, IndexId},
    quantity::RecoveryRate,
    refdata::{
        CommodityDefinition, CurrencyDefinition, DailyPublication, Deliverability, FxConvention,
        FxSpotLag, IndexDefinition, IndexKind, IndexStatus, InflationInterpolation, PriceKind,
        PublicationRule, RefData, RefDataEntry, RefDataParts, Seasonality,
    },
    rounding::{Rounding, RoundingMode},
    time::{Instant, LocalEventTime, UtcOffset, ZoneRules},
};
use rhodium_testkit::assert_close;

fn parts() -> Result<RefDataParts> {
    let calendar: CalendarId = "USNY".parse()?;
    let zone = "America/New_York".parse()?;
    let mut parts = RefDataParts::default();
    parts.calendars.insert(
        calendar.clone(),
        BaseCalendar::new(
            "2026-01-01".parse()?..="2026-12-31".parse()?,
            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
            vec!["2026-11-23".parse()?],
            Extrapolation::Strict,
        )?,
    );
    let est = UtcOffset::from_seconds(-5 * 3600)?;
    let edt = UtcOffset::from_seconds(-4 * 3600)?;
    parts.zones.insert(
        zone,
        ZoneRules::new(
            "2026-01-01T00:00:00Z".parse()?..="2026-12-31T23:59:59Z".parse()?,
            est,
            vec![
                ("2026-03-08T07:00:00Z".parse()?, edt),
                ("2026-11-01T06:00:00Z".parse()?, est),
            ],
        )?,
    );
    for code in [Currency::USD, Currency::JPY] {
        parts.currencies.insert(
            code,
            CurrencyDefinition {
                code,
                minor_units: 2,
                payment_rounding: Rounding::new(0.01, RoundingMode::HalfAwayFromZero)?,
                deliverability: Deliverability::Deliverable,
                settlement_time: LocalEventTime {
                    time: "17:00:00".parse()?,
                    zone: "America/New_York".parse()?,
                },
                settlement_calendar: calendar.clone(),
            },
        );
    }
    // USDJPY の FX fixing の受渡日は、ペアの慣行のスポットラグで決まる。
    add_usd_jpy(&mut parts)?;
    let calendar = CalendarSpec::Base(calendar);
    let definitions = [
        (
            "USD-SOFR",
            IndexKind::OvernightRate {
                ccy: Currency::USD,
                day_count: DayCount::Act360,
            },
            daily(calendar.clone(), 1, "08:00:00")?,
        ),
        (
            "USDJPY-WMR",
            IndexKind::FxFixing {
                pair: "USDJPY".parse()?,
            },
            PublicationRule::Intraday {
                delay_micros: 5_000_000,
            },
        ),
        (
            "US-CPI",
            IndexKind::InflationIndex {
                region: "US".into(),
                lag_months: 2,
                interpolation: InflationInterpolation::Flat,
                seasonality: Seasonality::None,
            },
            PublicationRule::Monthly {
                releases: BTreeMap::from([("2026-10".parse()?, "2026-11-12T13:30:00Z".parse()?)]),
            },
        ),
        (
            "USD-EVENT",
            IndexKind::CreditEvent {
                credit: credit_key()?,
            },
            daily(calendar.clone(), 1, "15:00:00")?,
        ),
    ];
    for (id, kind, publication) in definitions {
        let id: IndexId = id.parse()?;
        parts.indices.insert(
            id.clone(),
            IndexDefinition {
                id,
                kind,
                calendar: calendar.clone(),
                publication,
                observation_time: None,
                status: IndexStatus::Active,
            },
        );
    }
    Ok(parts)
}

/// USDJPY の受渡慣行（T+2、ニューヨーク 17 時、pip 0.01 円）を足す。
fn add_usd_jpy(parts: &mut RefDataParts) -> Result<()> {
    let usd_jpy: CcyPair = "USDJPY".parse()?;
    parts.fx_conventions.insert(
        usd_jpy,
        FxConvention {
            pair: usd_jpy,
            spot_lag: FxSpotLag::TwoBusinessDays,
            settlement: LocalEventTime {
                time: "17:00:00".parse()?,
                zone: "America/New_York".parse()?,
            },
            pip_exponent: 2,
        },
    );
    Ok(())
}

fn daily(calendar: CalendarSpec, business_days: u16, time: &str) -> Result<PublicationRule> {
    Ok(PublicationRule::Daily(DailyPublication {
        business_days,
        calendar,
        time: time.parse()?,
        zone: "America/New_York".parse()?,
    }))
}

fn credit_key() -> Result<CreditKey> {
    Ok(CreditKey {
        entity: "ACME".parse()?,
        seniority: "SENIOR".parse()?,
        doc_clause: "CR".parse()?,
        ccy: Currency::USD,
    })
}

#[test]
fn publication_uses_business_days_and_resolves_dst() -> Result<()> {
    // 架空の休業日が月曜にあるので、金曜分の翌営業日公表は火曜 08:00 EST。
    let id = "USD-SOFR".parse()?;
    let friday = ObservationPoint::Date("2026-11-20".parse()?);
    let refdata = RefData::new(parts()?)?;
    let (published_at, warnings) = refdata.publication(&id, friday)?;
    assert_eq!(published_at, "2026-11-24T13:00:00Z".parse()?);
    assert_eq!(warnings, []);

    // 春の 02:30 は欠落幅だけ進め、秋の 01:30 は先の UTC 候補を選ぶ。
    for (day, time, expected) in [
        ("2026-03-08", "02:30:00", "2026-03-08T07:30:00Z"),
        ("2026-11-01", "01:30:00", "2026-11-01T05:30:00Z"),
    ] {
        let mut parts = parts()?;
        let index = parts
            .indices
            .get_mut(&id)
            .ok_or_else(|| Error::Internal("missing test index".into()))?;
        index.publication = daily(index.calendar.clone(), 0, time)?;
        let (published_at, _) =
            RefData::new(parts)?.publication(&id, ObservationPoint::Date(day.parse()?))?;
        assert_eq!(published_at, expected.parse()?);
    }
    Ok(())
}

#[test]
fn intraday_monthly_and_credit_publications_have_distinct_anchors() -> Result<()> {
    let refdata = RefData::new(parts()?)?;
    // 5 秒の公表ラグは瞬間に足し、物価指数は参照月の表を引く。
    let intraday = ObservationPoint::Instant("2026-11-19T16:00:00Z".parse()?);
    let month = ObservationPoint::Month("2026-10".parse()?);
    assert_eq!(
        refdata.publication(&"USDJPY-WMR".parse()?, intraday)?.0,
        "2026-11-19T16:00:05Z".parse()?
    );
    assert_eq!(
        refdata.publication(&"US-CPI".parse()?, month)?.0,
        "2026-11-12T13:30:00Z".parse()?
    );
    assert!(matches!(
        refdata.publication(&"US-CPI".parse()?, intraday),
        Err(Error::InvalidData(_))
    ));
    let missing_month = ObservationPoint::Month(YearMonth::new(2026, 11)?);
    assert!(matches!(
        refdata.publication(&"US-CPI".parse()?, missing_month),
        Err(Error::MissingData(_))
    ));

    // 信用窓の終了日の翌営業日、ニューヨーク 15:00 に判定が公表される。
    let event_id = "USD-EVENT".parse()?;
    let friday = "2026-11-20".parse()?;
    let window = DateInterval::new("2026-11-19".parse()?, friday)?;
    assert_eq!(
        refdata.credit_publication(&event_id, window)?.0,
        "2026-11-24T20:00:00Z".parse()?
    );
    let key = credit_key()?;
    let event = CreditEventRecord {
        event_date: friday,
        kind: "failure_to_pay".into(),
    };
    let store = FixingStore::new(
        &refdata,
        FixingStoreParts {
            credit_events: vec![(key.clone(), event.clone())],
            recoveries: vec![(key.clone(), RecoveryRate::new(0.4)?)],
            ..FixingStoreParts::default()
        },
    )?;
    assert_eq!(store.credit_event(&key, window), Some(&event));
    assert!(
        store
            .credit_event(&key, DateInterval::new(friday, "2026-11-24".parse()?)?)
            .is_none()
    );
    let recovery = store
        .recovery(&key)
        .ok_or_else(|| Error::Internal("missing test recovery".into()))?;
    assert_close!(recovery.value(), 0.4, 1e-12, 1e-14);
    Ok(())
}

#[test]
fn fixing_switches_at_publication_and_preserves_snapshots() -> Result<()> {
    let refdata = RefData::new(parts()?)?;
    let id: IndexId = "USD-SOFR".parse()?;
    let point = ObservationPoint::Date("2026-11-20".parse()?);
    let (publication, _) = refdata.publication(&id, point)?;
    let just_before = Instant::from_unix_micros(publication.unix_micros() - 1)?;
    let just_after = Instant::from_unix_micros(publication.unix_micros() + 1)?;
    let period = ObservationPeriod::new(
        "2026-11-20T13:00:00Z".parse()?,
        "2026-11-24T12:00:00Z".parse()?,
    )?;
    let empty = FixingStore::new(&refdata, FixingStoreParts::default())?;
    let estimated = empty.with_estimates(&refdata, [(id.clone(), point, 0.038)])?;

    // 金曜分の期間は既に終了。公表の 1 μs 前は推定値だけを読み、警告が残る。
    assert!(matches!(
        empty.resolve(&id, point, just_before, publication, period),
        Err(Error::MissingData(_))
    ));
    let (value, warning) = estimated.resolve(&id, point, just_before, publication, period)?;
    assert!(
        matches!(warning, Some(Warning::FixingEstimateUsed { index, point: used }) if index == id && used == point)
    );
    assert!(matches!(value, FixingResolution::Estimated(v) if (v - 0.038).abs() < 1e-14));

    // 公表時刻の一致と通過後には推定値を使えず、公表値を入力したスナップショットだけ成功する。
    let updated = estimated.with_assumed_published(&refdata, [(id.clone(), point, 0.040)], [])?;
    for at in [publication, just_after] {
        assert!(matches!(
            estimated.resolve(&id, point, at, publication, period),
            Err(Error::MissingData(_))
        ));
        let (value, warning) = updated.resolve(&id, point, at, publication, period)?;
        assert!(matches!(value, FixingResolution::Published(v) if (v - 0.04).abs() < 1e-14));
        assert!(warning.is_none());
    }
    assert!(empty.estimate(&id, &point).is_none());
    assert!(estimated.published(&id, &point).is_none());
    Ok(())
}

#[test]
fn overnight_observation_period_uses_its_own_clock_and_business_calendar() -> Result<()> {
    let mut input = parts()?;
    let id: IndexId = "USD-SOFR".parse()?;
    let point = ObservationPoint::Date("2026-11-20".parse()?);
    assert!(matches!(
        RefData::new(input.clone())?.observation_period(&id, point),
        Err(Error::InvalidData(_))
    ));
    let definition = input
        .indices
        .get_mut(&id)
        .ok_or_else(|| Error::Internal("missing test index".into()))?;
    definition.observation_time = Some(LocalEventTime {
        time: "07:00:00".parse()?,
        zone: "America/New_York".parse()?,
    });
    let refdata = RefData::new(input)?;
    // 金曜の観測は休日の月曜を越えて火曜の現地 7 時まで。公表の時計は別に進む。
    let (period, warnings) = refdata.observation_period(&id, point)?;
    assert_eq!(period.start(), "2026-11-20T12:00:00Z".parse()?);
    assert_eq!(period.end(), "2026-11-24T12:00:00Z".parse()?);
    assert_ne!(refdata.publication(&id, point)?.0, period.end());
    assert_eq!(warnings, []);
    Ok(())
}

#[test]
fn unfinished_unpublished_period_uses_the_remaining_forward() -> Result<()> {
    // 評価時点までは推定値を使わず、期中に推定値が無ければ残存期間を予想する。
    let id: IndexId = "USD-SOFR".parse()?;
    let point = ObservationPoint::Date("2026-11-20".parse()?);
    let refdata = RefData::new(parts()?)?;
    let (publication, _) = refdata.publication(&id, point)?;
    let start = "2026-11-20T13:00:00Z".parse()?;
    let end = "2026-11-24T12:00:00Z".parse()?;
    let mid = "2026-11-21T00:00:00Z".parse()?;
    let period = ObservationPeriod::new(start, end)?;
    let empty = FixingStore::default();
    let estimated = empty.with_estimates(&refdata, [(id.clone(), point, 0.039)])?;
    assert_eq!(
        estimated.resolve(&id, point, start, publication, period)?,
        (FixingResolution::Forecast, None)
    );
    assert_eq!(
        empty.resolve(&id, point, mid, publication, period)?,
        (
            FixingResolution::ForecastFrom(mid),
            Some(Warning::UnpublishedFixingForecast {
                index: id.clone(),
                point
            }),
        )
    );
    assert!(
        matches!(estimated.resolve(&id, point, mid, publication, period)?.0,
        FixingResolution::Estimated(value) if (value - 0.039).abs() < 1e-14)
    );
    assert!(matches!(
        empty.resolve(&id, point, end, publication, period),
        Err(Error::MissingData(_))
    ));
    Ok(())
}

#[test]
fn refdata_and_fixing_input_reject_inconsistent_references() -> Result<()> {
    // 欠けた地域、指標と合わない公表規則、キーの種類の混同は読み込みで止める。
    let mut input = parts()?;
    input.zones.clear();
    assert!(
        matches!(RefData::new(input), Err(Error::InvalidData(message)) if message.contains("unknown zone"))
    );
    let mut input = parts()?;
    let id: IndexId = "USD-SOFR".parse()?;
    let index = input
        .indices
        .get_mut(&id)
        .ok_or_else(|| Error::Internal("missing test index".into()))?;
    index.publication = PublicationRule::Intraday { delay_micros: 1 };
    assert!(
        matches!(RefData::new(input), Err(Error::InvalidData(message)) if message.contains("incompatible publication"))
    );
    let refdata = RefData::new(parts()?)?;
    let point = ObservationPoint::Date("2026-11-20".parse::<Date>()?);
    for row in [
        ("USDJPY-WMR".parse()?, point, 150.0),
        ("US-CPI".parse()?, point, 300.0),
        ("UNKNOWN".parse()?, point, 0.01),
        (id.clone(), point, f64::NAN),
    ] {
        assert!(matches!(
            FixingStore::default().with_estimates(&refdata, [row]),
            Err(Error::InvalidData(_))
        ));
    }
    let row = (id, point, 0.03);
    assert!(
        matches!(FixingStore::default().with_estimates(&refdata, [row.clone(), row]),
        Err(Error::InvalidData(message)) if message.contains("duplicate fixing"))
    );
    Ok(())
}

// 既存の参照データに検査対象の指標を一つ足して構築する。日中の電力価格の指標に使う
// コモディティも登録し、種類と公表規則の組だけが構築の成否を分けるようにする。
fn refdata_with_index(kind: IndexKind, publication: PublicationRule) -> Result<RefData> {
    let mut parts = parts()?;
    let commodity: CommodityId = "PJM-WEST".parse()?;
    parts.commodities.insert(
        commodity.clone(),
        CommodityDefinition {
            id: commodity,
            unit: "MWh".into(),
            currency: Currency::USD,
        },
    );
    let id: IndexId = "TEST-INDEX".parse()?;
    parts.indices.insert(
        id.clone(),
        IndexDefinition {
            id,
            kind,
            calendar: CalendarSpec::Base("USNY".parse()?),
            publication,
            observation_time: None,
            status: IndexStatus::Active,
        },
    );
    RefData::new(parts)
}

fn power_price(price: PriceKind) -> Result<IndexKind> {
    Ok(IndexKind::AssetPrice {
        asset: Asset::Commodity("PJM-WEST".parse()?),
        ccy: Currency::USD,
        price,
    })
}

#[test]
fn index_kinds_accept_only_their_publication_rules() -> Result<()> {
    // 公表規則が観測点の種類（日次は観測日、日中は瞬間、月次は参照月）を決めるので、
    // 指標の種類ごとに許す規則を設計の表と照合する。物価指数は月次だけ、日中の価格は日中だけ、
    // FX の fixing は日次と日中、それ以外（終値、翌日物、信用事由など）は日次だけを許す。
    let rules = [
        daily(CalendarSpec::Base("USNY".parse()?), 1, "08:00:00")?,
        PublicationRule::Intraday {
            delay_micros: 5_000_000,
        },
        PublicationRule::Monthly {
            releases: BTreeMap::new(),
        },
    ];
    let inflation = IndexKind::InflationIndex {
        region: "US".into(),
        lag_months: 2,
        interpolation: InflationInterpolation::Flat,
        seasonality: Seasonality::None,
    };
    let fx = IndexKind::FxFixing {
        pair: "USDJPY".parse()?,
    };
    let overnight = IndexKind::OvernightRate {
        ccy: Currency::USD,
        day_count: DayCount::Act360,
    };
    let credit = IndexKind::CreditEvent {
        credit: credit_key()?,
    };
    // 期待値は [日次, 日中, 月次] の順に受理するか。
    for (kind, expected) in [
        (inflation, [false, false, true]),
        (fx, [true, true, false]),
        (power_price(PriceKind::Intraday)?, [false, true, false]),
        (power_price(PriceKind::Close)?, [true, false, false]),
        (overnight, [true, false, false]),
        (credit, [true, false, false]),
    ] {
        for (rule, accepted) in rules.iter().zip(expected) {
            match refdata_with_index(kind.clone(), rule.clone()) {
                Ok(_) => assert!(accepted, "{kind:?} accepted {rule:?}"),
                Err(Error::InvalidData(message))
                    if message.contains("incompatible publication") =>
                {
                    assert!(!accepted, "{kind:?} rejected {rule:?}");
                }
                Err(error) => return Err(error),
            }
        }
    }

    // 同じ FX の fixing でも、日次の規則なら観測日を、日中の規則なら瞬間を観測点として受け付ける。
    let id: IndexId = "TEST-INDEX".parse()?;
    let date = ObservationPoint::Date("2026-11-20".parse()?);
    let instant = ObservationPoint::Instant("2026-11-20T21:00:00Z".parse()?);
    let fx = IndexKind::FxFixing {
        pair: "USDJPY".parse()?,
    };
    let [daily_rule, intraday_rule, _] = rules;
    let daily_fx = refdata_with_index(fx.clone(), daily_rule)?;
    let intraday_fx = refdata_with_index(fx, intraday_rule)?;
    assert!(daily_fx.publication(&id, date).is_ok());
    assert!(matches!(
        daily_fx.publication(&id, instant),
        Err(Error::InvalidData(_))
    ));
    assert!(intraday_fx.publication(&id, instant).is_ok());
    assert!(matches!(
        intraday_fx.publication(&id, date),
        Err(Error::InvalidData(_))
    ));
    Ok(())
}

#[test]
fn compounded_rate_index_requires_an_overnight_rate() -> Result<()> {
    // 複利の指数の水準は元の翌日物の予想曲線から導くので、参照先は翌日物の指標に限る。
    let publication = daily(CalendarSpec::Base("USNY".parse()?), 1, "08:00:00")?;
    let sofr_index = IndexKind::CompoundedRateIndex {
        overnight: "USD-SOFR".parse()?,
    };
    assert!(refdata_with_index(sofr_index, publication.clone()).is_ok());
    let fx_index = IndexKind::CompoundedRateIndex {
        overnight: "USDJPY-WMR".parse()?,
    };
    assert!(matches!(
        refdata_with_index(fx_index, publication),
        Err(Error::InvalidData(message)) if message.contains("not an overnight rate")
    ));
    Ok(())
}

#[test]
fn zone_lookup_returns_registered_rules_and_rejects_unknown_zones() -> Result<()> {
    // 地域は RefData::zone で一度だけ引く。登録した遷移表は夏時間を含めて解決に使え、
    // 未登録の地域は参照の欠けとして InvalidData になる。
    let refdata = RefData::new(parts()?)?;
    let new_york = refdata.zone(&"America/New_York".parse()?)?;
    assert_eq!(
        new_york
            .to_local("2026-07-01T16:00:00Z".parse()?)?
            .to_string(),
        "2026-07-01T12:00:00"
    );
    assert!(matches!(
        refdata.zone(&"Asia/Tokyo".parse()?),
        Err(Error::InvalidData(message)) if message.contains("Asia/Tokyo")
    ));
    Ok(())
}

/// USDJPY の受渡慣行に、休日の違う決済カレンダーの USD と JPY を足した参照データ。
fn fx_parts() -> Result<RefDataParts> {
    let mut parts = parts()?;
    let tokyo: CalendarId = "JPTO".parse()?;
    parts.calendars.insert(
        tokyo.clone(),
        BaseCalendar::new(
            "2026-01-01".parse()?..="2026-12-31".parse()?,
            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
            vec!["2026-11-03".parse()?],
            Extrapolation::Strict,
        )?,
    );
    if let Some(jpy) = parts.currencies.get_mut(&Currency::JPY) {
        jpy.settlement_calendar = tokyo;
    }
    Ok(parts)
}

#[test]
fn fx_convention_lookup_accepts_either_direction() -> Result<()> {
    // 慣行は市場の標準の向き（USDJPY）で一つだけ持ち、JPYUSD で引いても同じ慣行を返す。
    // 受渡日を数える暦は両通貨の決済カレンダーの和で、向きによらず同じ式（同じキャッシュ）になる。
    let refdata = RefData::new(fx_parts()?)?;
    let usd_jpy: CcyPair = "USDJPY".parse()?;
    let convention = refdata.fx_convention(usd_jpy)?;
    assert!(ptr_eq(
        convention,
        refdata.fx_convention(usd_jpy.inverse())?
    ));
    assert_eq!(convention.pair, usd_jpy);
    assert_eq!(convention.spot_lag.business_days(), 2);
    assert_close!(convention.pip(), 0.01, 1e-15, 0.0);
    let spec = refdata.fx_calendar_spec(usd_jpy)?;
    assert_eq!(spec, refdata.fx_calendar_spec(usd_jpy.inverse())?);
    assert_eq!(
        spec,
        CalendarSpec::Union(vec![
            CalendarSpec::Base("JPTO".parse()?),
            CalendarSpec::Base("USNY".parse()?),
        ])
    );
    // 和の暦は、東京の祝日（11 月 3 日）とニューヨークの祝日（11 月 23 日）の両方で休む。
    let calendar = refdata.calendar(&spec)?;
    for holiday in ["2026-11-03", "2026-11-23"] {
        assert!(!calendar.is_business_day(holiday.parse()?)?.0);
    }
    assert!(matches!(
        refdata.fx_convention("EURUSD".parse()?),
        Err(Error::InvalidData(message)) if message.contains("unknown FX convention EURUSD")
    ));
    Ok(())
}

fn ptr_eq<T>(lhs: &T, rhs: &T) -> bool {
    std::ptr::eq(lhs, rhs)
}

#[test]
fn refdata_rejects_fx_conventions_in_both_directions() -> Result<()> {
    // 同じペアを両方の向きで持つと、pip とポイントをどちらの向きで読むかが決まらない。
    let mut parts = fx_parts()?;
    let usd_jpy: CcyPair = "USDJPY".parse()?;
    let mut inverse = parts.fx_conventions[&usd_jpy].clone();
    inverse.pair = usd_jpy.inverse();
    parts.fx_conventions.insert(inverse.pair, inverse);
    let Err(errors) = RefData::new_with_diagnostics(parts) else {
        return Err(Error::Internal("both directions were accepted".into()));
    };
    assert!(errors.iter().all(|(entry, error)| matches!(
        (entry, error),
        (RefDataEntry::FxConvention(_), Error::InvalidData(message))
            if message.contains("FX conventions for both")
    )));
    assert_eq!(errors.len(), 2);
    Ok(())
}

#[test]
fn refdata_requires_fx_convention_for_fx_fixing_index() -> Result<()> {
    // FX fixing の受渡日はペアのスポットラグで決まるので、慣行の無いペアの fixing 指標は
    // 参照データを作るときに拒否する。どちらの向きの慣行でもよい。
    let mut parts = fx_parts()?;
    parts.fx_conventions.clear();
    let Err(errors) = RefData::new_with_diagnostics(parts) else {
        return Err(Error::Internal(
            "FX fixing without a convention was accepted".into(),
        ));
    };
    assert!(errors.iter().any(|(entry, error)| matches!(
        (entry, error),
        (RefDataEntry::Index(id), Error::InvalidData(message))
            if id.as_str() == "USDJPY-WMR" && message.contains("unknown FX convention USDJPY")
    )));
    let mut parts = fx_parts()?;
    let usd_jpy: CcyPair = "USDJPY".parse()?;
    let mut convention = parts
        .fx_conventions
        .remove(&usd_jpy)
        .ok_or_else(|| Error::Internal("missing USDJPY convention".into()))?;
    convention.pair = usd_jpy.inverse();
    parts.fx_conventions.insert(convention.pair, convention);
    RefData::new(parts)?;
    Ok(())
}

#[test]
fn refdata_rejects_fx_convention_key_mismatch_and_unknown_zone() -> Result<()> {
    // 表のキーと慣行の向きが違う定義、受渡時刻の地域が未登録の定義を、慣行の項目として報告する。
    let usd_jpy: CcyPair = "USDJPY".parse()?;
    let mut parts = fx_parts()?;
    let convention = parts.fx_conventions[&usd_jpy].clone();
    parts.fx_conventions.insert("EURJPY".parse()?, convention);
    let Err(errors) = RefData::new_with_diagnostics(parts) else {
        return Err(Error::Internal("mismatched key was accepted".into()));
    };
    assert!(errors.iter().any(|(entry, _)| matches!(
        entry,
        RefDataEntry::FxConvention(pair) if pair.to_string() == "EURJPY"
    )));

    let mut parts = fx_parts()?;
    if let Some(convention) = parts.fx_conventions.get_mut(&usd_jpy) {
        convention.settlement.zone = "Asia/Tokyo".parse()?;
    }
    let Err(errors) = RefData::new_with_diagnostics(parts) else {
        return Err(Error::Internal("unknown zone was accepted".into()));
    };
    assert!(errors.iter().any(|(entry, error)| matches!(
        (entry, error),
        (RefDataEntry::FxConvention(pair), Error::InvalidData(message))
            if *pair == usd_jpy && message.contains("Asia/Tokyo")
    )));
    Ok(())
}
`````

`````
crates/rhodium-core/tests/schedules.rs
//! 計算期間と決済の瞬間を独立に求めた値と照合する。

#[path = "schedules/observation.rs"]
mod observation;

use std::collections::BTreeMap;

use proptest::{
    prelude::*,
    test_runner::{Config, RngSeed},
};
use rhodium_core::{
    Error, Result, Warning,
    calendar::{
        BaseCalendar, BusinessDayConvention as Bdc, CalendarSpec, Extrapolation, WeekendRule,
    },
    date::{Date, DateInterval, Weekday},
    day_count::{AccrualContext, DayCount, PeriodFrame, year_fraction},
    id::CalendarId,
    schedule::{
        AccrualDates, DayOfMonth, EndRule, ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods,
        GeneratedPeriods, MAX_PERIODS, PaymentGrouping, PaymentRule, Period, PeriodAnchor,
        PeriodCombination, PeriodKind, ResetRule, RollConvention, Schedule, ScheduleSpec,
        SettlementContext, StartRule, StubLength, StubSpec, WeekOfMonth, generate, reset_dates,
    },
    tenor::{Frequency, TenorLength, TenorUnit},
    time::{LocalEventTime, TzId, UtcOffset, ZoneRules},
};

struct Dates {
    id: CalendarId,
    calendars: BTreeMap<CalendarId, BaseCalendar>,
    calendar: CalendarSpec,
    settlement_time: LocalEventTime,
    zones: BTreeMap<TzId, ZoneRules>,
}

impl Dates {
    fn new(holidays: &[&str]) -> Result<Self> {
        let id: CalendarId = "JPTO".parse()?;
        let calendar = BaseCalendar::new(
            "2026-01-01".parse()?..="2040-12-31".parse()?,
            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
            holidays
                .iter()
                .map(|date| date.parse())
                .collect::<Result<_>>()?,
            Extrapolation::Strict,
        )?;
        let calendars = BTreeMap::from([(id.clone(), calendar)]);
        let calendar = CalendarSpec::Base(id.clone());
        let settlement_time = event("15:00:00", "Asia/Tokyo")?;
        let coverage = "2026-01-01T00:00:00Z".parse()?..="2100-12-31T23:59:59Z".parse()?;
        let zones = BTreeMap::from([
            (
                settlement_time.zone.clone(),
                ZoneRules::new(coverage.clone(), UtcOffset::from_seconds(9 * 3600)?, vec![])?,
            ),
            (
                "Europe/Berlin".parse()?,
                ZoneRules::new(
                    coverage,
                    UtcOffset::from_seconds(3600)?,
                    vec![
                        (
                            "2027-03-28T01:00:00Z".parse()?,
                            UtcOffset::from_seconds(7200)?,
                        ),
                        (
                            "2027-10-31T01:00:00Z".parse()?,
                            UtcOffset::from_seconds(3600)?,
                        ),
                    ],
                )?,
            ),
        ]);
        Ok(Self {
            id,
            calendars,
            calendar,
            settlement_time,
            zones,
        })
    }

    fn spec(&self) -> CalendarSpec {
        CalendarSpec::Base(self.id.clone())
    }

    fn settlement(&self) -> SettlementContext<'_> {
        SettlementContext {
            zones: &self.zones,
            default_time: Some(&self.settlement_time),
            calendar: Some(&self.calendar),
        }
    }

    fn generate(&self, spec: &ScheduleSpec) -> Result<(Schedule, Vec<Warning>)> {
        generate(
            spec,
            "2027-01-04".parse()?,
            &self.calendars,
            &self.settlement(),
        )
    }

    fn payment(&self) -> PaymentRule {
        PaymentRule {
            anchor: PeriodAnchor::PeriodEnd,
            offset: rhodium_core::tenor::DateOffset::new(0, TenorUnit::Day),
            calendar: self.spec(),
            bdc: Bdc::Unadjusted,
            grouping: PaymentGrouping::EachPeriod,
            settlement: None,
        }
    }
}

fn event(time: &str, zone: &str) -> Result<LocalEventTime> {
    Ok(LocalEventTime {
        time: time.parse()?,
        zone: zone.parse()?,
    })
}

fn interval(start: &str, end: &str) -> Result<DateInterval> {
    DateInterval::new(start.parse()?, end.parse()?)
}

fn generated(
    dates: &Dates,
    interval: DateInterval,
    frequency: Frequency,
    roll: Option<RollConvention>,
    stub: StubSpec,
    bdc: Bdc,
) -> Result<ScheduleSpec> {
    Ok(ScheduleSpec::Generated {
        periods: GeneratedPeriods::new(
            StartRule::Date(interval.start()),
            EndRule::Date(interval.end()),
            frequency,
            roll,
            stub,
            dates.spec(),
            bdc,
            AccrualDates::Adjusted,
        )?,
        payment: dates.payment(),
    })
}

fn assert_close(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= (1e-12 * scale).max(1e-14),
        "{actual} != {expected}"
    );
}

proptest! {
    #![proptest_config(Config { cases: 256, rng_seed: RngSeed::Fixed(157), ..Config::default() })]

    #[test]
    fn schedule_invariants(
        day in 1..=28u8,
        month in 1..=12u8,
        steps in 1..=12i32,
        months_per_step in prop::sample::select(vec![1, 3, 6, 12]),
        residual in prop_oneof![Just(0), 7..=20i32],
        stub_first in any::<bool>(),
        length in prop::sample::select(vec![StubLength::Short, StubLength::Long]),
        bdc in prop::sample::select(vec![
            Bdc::Unadjusted,
            Bdc::Following,
            Bdc::ModifiedFollowing,
            Bdc::Preceding,
        ]),
        accrual in prop::sample::select(vec![AccrualDates::Adjusted, AccrualDates::Unadjusted]),
    ) {
        // どのスタブ・刻み・調整規則でも、未調整の境界は格子と余りから手で数えた日付になり、
        // 未調整・調整後・利息計算の区間は連続して狭義に増え、支払いは各期間を一度ずつ覆う。
        // 境界は 7 日以上離れ、週末だけの暦の調整は 2 日以内なので、調整で境界が衝突することはない。
        let dates = Dates::new(&[])?;
        // 格子点 G_i = G_0 + i × 刻み（i = 0..=steps）。日は 28 日以下なので月末に詰めない。
        let grid = (0..=steps)
            .map(|i| months_after(month, day, i * months_per_step))
            .collect::<Result<Vec<_>>>()?;
        let (first, last) = (grid[0], grid[grid.len() - 1]);
        // Initial は終了日 G_k から遡り、Final は開始日 G_0 から進む。余りは生成の向きの先に置く。
        let (start, end, stub) = if stub_first {
            (first.checked_add_days(-residual)?, last, StubSpec::Initial(length))
        } else {
            (first, last.checked_add_days(residual)?, StubSpec::Final(length))
        };
        // 期待する境界: 余りが無ければ格子そのもの。Short は余りの側に端点を足し、Long は余りの側で
        // 隣の格子点を端点に置き換えて、最初か最後の通常期間に余りを合わせる。
        let mut boundaries = grid;
        let stub_kind = match (residual > 0, stub_first, length) {
            (false, _, _) => None,
            (true, true, StubLength::Short) => {
                boundaries.insert(0, start);
                Some(PeriodKind::ShortStub)
            }
            (true, true, StubLength::Long) => {
                boundaries[0] = start;
                Some(PeriodKind::LongStub)
            }
            (true, false, StubLength::Short) => {
                boundaries.push(end);
                Some(PeriodKind::ShortStub)
            }
            (true, false, StubLength::Long) => {
                let last = boundaries.len() - 1;
                boundaries[last] = end;
                Some(PeriodKind::LongStub)
            }
        };
        let spec = ScheduleSpec::Generated {
            periods: GeneratedPeriods::new(
                StartRule::Date(start),
                EndRule::Date(end),
                format!("{months_per_step}M").parse()?,
                None,
                stub,
                dates.spec(),
                bdc,
                accrual,
            )?,
            payment: dates.payment(),
        };
        let (schedule, _) = dates.generate(&spec)?;
        let periods = schedule.periods();
        prop_assert_eq!(
            periods.iter().map(|p| p.unadjusted().start()).chain([end]).collect::<Vec<_>>(),
            boundaries
        );
        for pair in periods.windows(2) {
            prop_assert_eq!(pair[0].unadjusted().end(), pair[1].unadjusted().start());
            prop_assert_eq!(pair[0].adjusted().end(), pair[1].adjusted().start());
            prop_assert_eq!(pair[0].accrual().end(), pair[1].accrual().start());
        }
        let stub_index = if stub_first { 0 } else { periods.len() - 1 };
        for (index, period) in periods.iter().enumerate() {
            prop_assert!(!period.unadjusted().is_empty() && !period.adjusted().is_empty());
            let expected_accrual = match accrual {
                AccrualDates::Adjusted => period.adjusted(),
                AccrualDates::Unadjusted => period.unadjusted(),
            };
            prop_assert_eq!(period.accrual(), expected_accrual);
            // 余りの側の期間だけがスタブで、その他は通常の期間。
            let kind = stub_kind.filter(|_| index == stub_index).unwrap_or(PeriodKind::Regular);
            prop_assert_eq!(period.kind(), kind);
            // 参照枠は利息計算の区間を覆い、長いスタブだけが仮想の通常期間を 2 つ持つ。
            let frame = period.frame().map_or_default(PeriodFrame::regular);
            prop_assert_eq!(frame.len(), if kind == PeriodKind::LongStub { 2 } else { 1 });
            prop_assert!(frame[0].start() <= period.accrual().start());
            prop_assert!(frame[frame.len() - 1].end() >= period.accrual().end());
            if kind == PeriodKind::Regular {
                prop_assert_eq!(frame, &[period.accrual()]);
            }
        }
        // 支払いは期間ごとに一つで、期末の 0 日後（調整なし）に払う。
        prop_assert_eq!(schedule.payments().len(), periods.len());
        for (index, (period, payment)) in periods.iter().zip(schedule.payments()).enumerate() {
            prop_assert_eq!(payment.periods(), index..index + 1);
            prop_assert_eq!(payment.date(), period.adjusted().end());
        }
    }
}

fn months_after(month: u8, day: u8, months: i32) -> Result<Date> {
    // 起点は 2027 年の month 月。年月は格子の実装に頼らず、月の通し番号から求める。
    let serial = i32::from(month) - 1 + months;
    let month = u8::try_from(serial.rem_euclid(12) + 1)
        .map_err(|error| Error::Internal(error.to_string()))?;
    Date::from_ymd(2027 + serial.div_euclid(12), month, day)
}

#[test]
fn stubs_preserve_the_regular_grid_and_group_payments() -> Result<()> {
    let dates = Dates::new(&[])?;
    let start = "2027-01-31".parse()?;
    let spec = ScheduleSpec::Generated {
        periods: GeneratedPeriods::new(
            StartRule::Date(start),
            EndRule::Date("2027-06-15".parse()?),
            "1M".parse()?,
            Some(RollConvention::EndOfMonth),
            StubSpec::Final(StubLength::Long),
            dates.spec(),
            Bdc::ModifiedFollowing,
            AccrualDates::Unadjusted,
        )?,
        payment: PaymentRule {
            offset: "2BD".parse()?,
            grouping: PaymentGrouping::EveryN {
                periods: 3
                    .try_into()
                    .map_err(|_| Error::Internal("invalid group size".into()))?,
                combine: PeriodCombination::StraightCompounding,
            },
            ..dates.payment()
        },
    };
    let (schedule, _) = dates.generate(&spec)?;
    // 5 月末で終わる通常期間に 6 月 15 日までの余りを足し、参照枠は 6 月末まで残す。
    let stub = &schedule.periods()[3];
    assert_eq!(stub.kind(), PeriodKind::LongStub);
    assert_eq!(stub.accrual(), interval("2027-04-30", "2027-06-15")?);
    assert_eq!(
        stub.frame().map(PeriodFrame::regular),
        Some(
            [
                interval("2027-04-30", "2027-05-31")?,
                interval("2027-05-31", "2027-06-30")?,
            ]
            .as_slice()
        )
    );
    assert_eq!(schedule.combine(), PeriodCombination::StraightCompounding);
    assert_eq!(schedule.payments()[0].periods(), 0..3);
    assert_eq!(schedule.payments()[0].date(), "2027-05-04".parse()?);
    assert_eq!(schedule.payments()[0].at(), "2027-05-04T06:00:00Z".parse()?);
    assert_eq!(schedule.payments()[1].periods(), 3..4);
    let reset = ResetRule {
        anchor: PeriodAnchor::PeriodStart,
        offset: "-2BD".parse()?,
        calendar: dates.spec(),
        bdc: Bdc::Preceding,
    };
    let (resets, _) = reset_dates(&schedule, &reset, &dates.calendars)?;
    // 1 月末の日曜は 29 日に前倒しされ、2 営業日前は 27 日。
    assert_eq!(resets[0], "2027-01-27".parse()?);
    Ok(())
}

#[test]
fn all_stub_positions_preserve_boundaries_and_references() -> Result<()> {
    use PeriodKind::{LongStub, Regular, ShortStub};
    let dates = Dates::new(&[])?;
    let jan = "2027-01-31".parse()?;
    let mar = "2027-03-31".parse()?;
    // 指定された通常期間を中心に、端の余りと仮想期間を手で数える。
    for (start, end, stub, kinds) in [
        (
            "2027-01-17",
            "2027-03-31",
            StubSpec::Initial(StubLength::Short),
            vec![ShortStub, Regular, Regular],
        ),
        (
            "2027-01-17",
            "2027-03-31",
            StubSpec::Initial(StubLength::Long),
            vec![LongStub, Regular],
        ),
        (
            "2027-01-17",
            "2027-03-31",
            StubSpec::FirstRegular(jan),
            vec![ShortStub, Regular, Regular],
        ),
        (
            "2027-01-31",
            "2027-04-17",
            StubSpec::LastRegular(mar),
            vec![Regular, Regular, ShortStub],
        ),
        (
            "2027-01-17",
            "2027-04-17",
            StubSpec::Both {
                first_regular: jan,
                last_regular: mar,
            },
            vec![ShortStub, Regular, Regular, ShortStub],
        ),
        (
            "2027-01-31",
            "2027-04-30",
            StubSpec::LastRegular(mar),
            vec![Regular, Regular, Regular],
        ),
        (
            "2027-01-31",
            "2027-04-17",
            StubSpec::Final(StubLength::Short),
            vec![Regular, Regular, ShortStub],
        ),
    ] {
        let (schedule, _) = dates.generate(&generated(
            &dates,
            interval(start, end)?,
            "1M".parse()?,
            Some(RollConvention::EndOfMonth),
            stub,
            Bdc::Unadjusted,
        )?)?;
        assert_eq!(
            schedule
                .periods()
                .iter()
                .map(Period::kind)
                .collect::<Vec<_>>(),
            kinds
        );
        let references: Vec<_> = schedule
            .periods()
            .iter()
            .map(|p| p.frame().map(|f| f.regular().len()))
            .collect();
        let expected: Vec<_> = kinds
            .iter()
            .map(|kind| Some(if *kind == LongStub { 2 } else { 1 }))
            .collect();
        assert_eq!(references, expected);
    }
    Ok(())
}

#[test]
fn invalid_grids_and_adjustment_collisions_are_rejected() -> Result<()> {
    let dates = Dates::new(&[])?;
    // 割り切れない通常格子、合算相手のない長いスタブ、週末の調整で同日にそろう境界を拒否する。
    for (start, end, step, stub, bdc) in [
        (
            "2027-01-29",
            "2027-02-17",
            "1M",
            StubSpec::NoStub,
            Bdc::Unadjusted,
        ),
        (
            "2027-01-29",
            "2027-02-17",
            "1M",
            StubSpec::Final(StubLength::Long),
            Bdc::Unadjusted,
        ),
        (
            "2027-01-29",
            "2027-02-01",
            "1D",
            StubSpec::NoStub,
            Bdc::Following,
        ),
    ] {
        let spec = generated(
            &dates,
            interval(start, end)?,
            step.parse()?,
            None,
            stub,
            bdc,
        )?;
        assert!(matches!(dates.generate(&spec), Err(Error::InvalidData(_))));
    }
    let start: Date = "2027-01-01".parse()?;
    let spec = generated(
        &dates,
        DateInterval::new(start, start.checked_add_days(10_001)?)?,
        "1D".parse()?,
        None,
        StubSpec::NoStub,
        Bdc::Unadjusted,
    )?;
    assert_eq!(MAX_PERIODS, 10_000);
    assert!(matches!(dates.generate(&spec), Err(Error::InvalidInput(_))));
    Ok(())
}

#[test]
fn tenor_end_on_the_roll_and_frn_drift_preserve_contract_dates() -> Result<()> {
    let dates = Dates::new(&[])?;
    // 月末の開始日から 1 年後は閏年の 2 月 29 日。
    let eom = ScheduleSpec::Generated {
        periods: GeneratedPeriods::new(
            StartRule::Date("2027-02-28".parse()?),
            EndRule::AfterStart(TenorLength::new(1, TenorUnit::Year)),
            Frequency::ANNUAL,
            Some(RollConvention::EndOfMonth),
            StubSpec::NoStub,
            dates.spec(),
            Bdc::Unadjusted,
            AccrualDates::Unadjusted,
        )?,
        payment: dates.payment(),
    };
    let (schedule, _) = dates.generate(&eom)?;
    assert_eq!(
        schedule.periods()[0].accrual(),
        interval("2027-02-28", "2028-02-29")?
    );
    // 5 月 15 日が土曜日なので 17 日に調整される。次の刻みは 6 月 17 日だが、確認書の満期 15 日で閉じる。
    let (frn, _) = dates.generate(&generated(
        &dates,
        interval("2027-01-15", "2027-06-15")?,
        "1M".parse()?,
        Some(RollConvention::Frn),
        StubSpec::NoStub,
        Bdc::ModifiedFollowing,
    )?)?;
    assert_eq!(
        frn.periods()[4].adjusted(),
        interval("2027-05-17", "2027-06-15")?
    );
    assert_eq!(frn.periods()[4].kind(), PeriodKind::Regular);
    // 2 月で月末に着くと 3 月も月末営業日を使う。
    let (sticky, _) = dates.generate(&generated(
        &dates,
        interval("2027-01-29", "2027-03-31")?,
        "1M".parse()?,
        Some(RollConvention::Frn),
        StubSpec::NoStub,
        Bdc::ModifiedFollowing,
    )?)?;
    assert_eq!(
        sticky.periods()[1].adjusted(),
        interval("2027-02-26", "2027-03-31")?
    );
    Ok(())
}

#[test]
fn frn_drift_closing_short_of_a_step_keeps_the_schedule_without_a_frame() -> Result<()> {
    // FRN の慣行では日付が調整後の日付から進むので、15 日始まりの月払いは週末のたびに後ろへずれ、
    // 月末に着くと月末営業日に張り付く。確認書の満期 15 日で閉じる最後の期間は、2029-12-31→2030-01-15 の
    // 15 日しかなく、それ自身を通常期間とみなせない（最も近い刻みの数が 0）。日程は作れ、
    // その期間だけ参照枠を持たないので、枠を読まない ACT/360 は使え、ICMA と ACT/365L は構築で拒否する
    // （以前は枠の検査で日程の生成ごと失敗した。さらに前は ICMA で黙って 1/12 になった）。
    let dates = Dates::new(&[])?;
    let (frn, _) = dates.generate(&generated(
        &dates,
        interval("2027-01-15", "2030-01-15")?,
        "1M".parse()?,
        Some(RollConvention::Frn),
        StubSpec::NoStub,
        Bdc::ModifiedFollowing,
    )?)?;
    let last = frn
        .periods()
        .last()
        .ok_or_else(|| Error::Internal("empty FRN schedule".into()))?;
    assert_eq!(last.adjusted(), interval("2029-12-31", "2030-01-15")?);
    assert_eq!(last.kind(), PeriodKind::Regular);
    assert!(last.frame().is_none());
    // ずれても一刻みに近い期間は、従来どおりそれ自身を参照枠に持つ。
    assert!(
        frn.periods()[..frn.periods().len() - 1]
            .iter()
            .all(|p| p.frame().map(PeriodFrame::regular) == Some(&[p.accrual()][..]))
    );
    frn.validate_day_count(&DayCount::Act360)?;
    for basis in [DayCount::ActActIcma, DayCount::Act365L] {
        assert!(matches!(
            frn.validate_day_count(&basis),
            Err(Error::InvalidData(_))
        ));
    }
    Ok(())
}

#[test]
fn trade_date_offset_and_imm_roll_produce_expected_dates() -> Result<()> {
    let dates = Dates::new(&["2027-01-05"])?;
    let spec = ScheduleSpec::Generated {
        periods: GeneratedPeriods::new(
            StartRule::FromTradeDate("2BD".parse()?),
            EndRule::AfterStart(TenorLength::new(1, TenorUnit::Week)),
            Frequency::ONCE,
            Some(RollConvention::OnDay(DayOfMonth::new(6)?)),
            StubSpec::NoStub,
            dates.spec(),
            Bdc::Unadjusted,
            AccrualDates::Adjusted,
        )?,
        payment: dates.payment(),
    };
    let (schedule, _) = dates.generate(&spec)?;
    // 月曜の約定から、火曜の休日を飛ばして木曜に開始する。
    assert_eq!(
        schedule.periods()[0].accrual(),
        interval("2027-01-07", "2027-01-14")?
    );
    let (imm, _) = dates.generate(&generated(
        &dates,
        interval("2027-03-17", "2027-09-15")?,
        "3M".parse()?,
        Some(RollConvention::NthWeekday {
            week: WeekOfMonth::Third,
            weekday: Weekday::Wed,
        }),
        StubSpec::NoStub,
        Bdc::Unadjusted,
    )?)?;
    assert_eq!(imm.periods()[0].unadjusted().end(), "2027-06-16".parse()?);
    Ok(())
}

fn explicit(start: &str, end: &str, payment: &str) -> Result<ExplicitPeriod> {
    Ok(ExplicitPeriod {
        accrual: interval(start, end)?,
        payment: payment.parse()?,
        settlement: None,
        references: Vec::new(),
        reset: None,
    })
}

#[test]
fn explicit_stub_reference_period() -> Result<()> {
    let dates = Dates::new(&[])?;
    // 半年払いの 1 か月のスタブを、頻度も参照期間も無い場合、頻度だけの場合、両方ある場合で比べる。
    // 参照期間を復元できない前の二つは枠の要る日数計算を InvalidData にし、ACT/365F は使えるままにする。
    // 両方あれば、accrual を通常期間とみなした 1/2 ではなく、参照期間で数えた ICMA の係数になる。
    let stub = explicit("2027-02-28", "2027-03-31", "2027-03-31")?;
    for (reference, frequency) in [
        (None, None),
        (None, Some("6M".parse()?)),
        (
            Some(interval("2027-02-28", "2027-08-31")?),
            Some("6M".parse()?),
        ),
    ] {
        let spec = ScheduleSpec::Explicit(ExplicitPeriods::new(
            vec![ExplicitPeriod {
                references: reference.into_iter().collect(),
                ..stub.clone()
            }],
            ExplicitDatePolicy::KeepAsConfirmed,
            PeriodCombination::Sum,
            frequency,
        )?);
        let (schedule, _) = dates.generate(&spec)?;
        let dc = DayCount::ActActIcma;
        if reference.is_none() {
            assert!(matches!(
                schedule.validate_day_count(&dc),
                Err(Error::InvalidData(_))
            ));
            assert!(matches!(
                schedule.validate_day_count(&DayCount::Act365L),
                Err(Error::InvalidData(_))
            ));
            schedule.validate_day_count(&DayCount::Act365F)?;
        } else {
            schedule.validate_day_count(&dc)?;
            let ctx = schedule
                .accrual_context(0, None)
                .ok_or_else(|| Error::Internal("missing accrual context".into()))?;
            let (fraction, _) = year_fraction(&dc, stub.accrual.start(), stub.accrual.end(), &ctx)?;
            // 2 月 28 日から 3 月 31 日まで 31 日。半年の参照期間は 184 日。
            assert_close(fraction.value(), 31.0 / (2.0 * 184.0));
        }
    }
    let regular = explicit("2028-02-29", "2028-08-31", "2028-08-31")?;
    let spec = ScheduleSpec::Explicit(ExplicitPeriods::new(
        vec![regular],
        ExplicitDatePolicy::KeepAsConfirmed,
        PeriodCombination::Sum,
        Some("6M".parse()?),
    )?);
    let (schedule, _) = dates.generate(&spec)?;
    // 開始日が月末なので、参照期間を省略しても 8 月 31 日が通常期間の終了日となる。
    schedule.validate_day_count(&DayCount::Act365L)?;
    let context = schedule
        .accrual_context(0, None)
        .ok_or_else(|| Error::Internal("missing accrual context".into()))?;
    let p = schedule.periods()[0].accrual();
    let (fraction, _) = year_fraction(&DayCount::Act365L, p.start(), p.end(), &context)?;
    assert_close(fraction.value(), 184.0 / 366.0);
    Ok(())
}

#[test]
fn explicit_long_stub_carries_several_reference_periods() -> Result<()> {
    // 半年払いの長い後ろのスタブ 2027-01-01→10-01 を、確認書の明示の期間で一つのクーポンとして書く。
    // 参照期間 01-01→07-01（181 日のうち 181 日、係数 1/2）と 07-01→2028-01-01（184 日のうち 92 日、
    // 係数 92/(2·184) = 1/4）を一つの期間に並べると 3/4 になり、規則から生成した長い後ろのスタブ
    // （同じ区間を一つの期間にし、二つの通常期間を枠に持つ）と期間の数も係数も一致する。
    // 二つの期間に分けて Sum で足しても固定の係数は 3/4 だが、RFR の複利は分割点で切れて
    // 別のクーポンになる（pricing の explicit_long_stub_matches_generated_rfr_coupon）。
    let dates = Dates::new(&[])?;
    let dc = DayCount::ActActIcma;
    let fractions = |schedule: &Schedule| -> Result<f64> {
        schedule
            .periods()
            .iter()
            .enumerate()
            .map(|(i, period)| {
                let context = schedule
                    .accrual_context(i, None)
                    .ok_or_else(|| Error::Internal("missing accrual context".into()))?;
                let accrual = period.accrual();
                Ok(
                    year_fraction(&dc, accrual.start(), accrual.end(), &context)?
                        .0
                        .value(),
                )
            })
            .sum()
    };
    let explicit_periods = |periods, frequency| {
        ExplicitPeriods::new(
            periods,
            ExplicitDatePolicy::KeepAsConfirmed,
            PeriodCombination::Sum,
            frequency,
        )
    };
    let semiannual = Some("6M".parse()?);
    let first = interval("2027-01-01", "2027-07-01")?;
    let second = interval("2027-07-01", "2028-01-01")?;
    let one_coupon = vec![ExplicitPeriod {
        references: vec![first, second],
        ..explicit("2027-01-01", "2027-10-01", "2027-10-01")?
    }];
    let (schedule, _) = dates.generate(&ScheduleSpec::Explicit(explicit_periods(
        one_coupon, semiannual,
    )?))?;
    let (rule, _) = dates.generate(&generated(
        &dates,
        interval("2027-01-01", "2027-10-01")?,
        "6M".parse()?,
        None,
        StubSpec::Final(StubLength::Long),
        Bdc::Unadjusted,
    )?)?;
    assert_eq!(schedule.periods().len(), 1);
    assert_eq!(rule.periods().len(), 1);
    assert_eq!(
        schedule.periods()[0].frame().map(PeriodFrame::regular),
        rule.periods()[0].frame().map(PeriodFrame::regular)
    );
    assert_close(fractions(&schedule)?, 0.75);
    assert_close(fractions(&rule)?, fractions(&schedule)?);

    // 二つの期間に分けた入力も書けるが、別の二つの計算期間になる（固定の係数の和は同じ 3/4）。
    let split = vec![
        explicit("2027-01-01", "2027-07-01", "2027-10-01")?,
        ExplicitPeriod {
            references: vec![second],
            ..explicit("2027-07-01", "2027-10-01", "2027-10-01")?
        },
    ];
    let (split, _) = dates.generate(&ScheduleSpec::Explicit(explicit_periods(
        split, semiannual,
    )?))?;
    assert_eq!(split.periods().len(), 2);
    assert_eq!(split.payments().len(), 1);
    assert_close(fractions(&split)?, 0.75);

    // 区間を覆わない、連続しない、頻度の無い参照期間の列は拒否する。
    let invalid = [
        (vec![first], semiannual),
        (
            vec![first, interval("2027-08-01", "2028-02-01")?],
            semiannual,
        ),
        (vec![first, second], None),
    ];
    for (references, frequency) in invalid {
        let period = ExplicitPeriod {
            references,
            ..explicit("2027-01-01", "2027-10-01", "2027-10-01")?
        };
        assert!(matches!(
            explicit_periods(vec![period], frequency),
            Err(Error::InvalidData(_))
        ));
    }
    // 月末から始まる参照期間は、月末に揃える慣行の有無のどちらの期末（08-28、08-31）も
    // 通常期間として受け付ける。どちらも利息を数える 02-28→03-31 を覆う。
    for end in ["2027-08-28", "2027-08-31"] {
        let stub = ExplicitPeriod {
            references: vec![interval("2027-02-28", end)?],
            ..explicit("2027-02-28", "2027-03-31", "2027-03-31")?
        };
        explicit_periods(vec![stub], semiannual)?;
    }
    Ok(())
}

#[test]
fn explicit_reference_accepts_rolled_and_adjusted_regular_periods() -> Result<()> {
    // 毎月第 3 水曜にロールする月払いの日程では、2026-11-18→12-16 の 28 日も通常期間になる。
    // 同じ区間を明示の参照期間に書いても受理され、生成した日程と同じ ICMA の係数 28/(12·28) = 1/12 になる。
    // 暦の上でちょうど一か月ではないので、以前の明示日程の判定（暦の一刻み）はこれを拒否していた。
    let dates = Dates::new(&[])?;
    let dc = DayCount::ActActIcma;
    let fraction = |schedule: &Schedule, index: usize| -> Result<f64> {
        let context = schedule
            .accrual_context(index, None)
            .ok_or_else(|| Error::Internal("missing accrual context".into()))?;
        let accrual = schedule.periods()[index].accrual();
        Ok(
            year_fraction(&dc, accrual.start(), accrual.end(), &context)?
                .0
                .value(),
        )
    };
    let (rule, _) = dates.generate(&generated(
        &dates,
        interval("2026-10-21", "2027-01-20")?,
        "1M".parse()?,
        Some(RollConvention::NthWeekday {
            week: WeekOfMonth::Third,
            weekday: Weekday::Wed,
        }),
        StubSpec::NoStub,
        Bdc::Unadjusted,
    )?)?;
    let rolled = interval("2026-11-18", "2026-12-16")?;
    assert_eq!(rule.periods()[1].accrual(), rolled);
    let explicit_schedule = |accrual: DateInterval, frequency: &str| -> Result<Schedule> {
        let period = ExplicitPeriod {
            references: vec![accrual],
            ..explicit(
                &accrual.start().to_string(),
                &accrual.end().to_string(),
                &accrual.end().to_string(),
            )?
        };
        let spec = ScheduleSpec::Explicit(ExplicitPeriods::new(
            vec![period],
            ExplicitDatePolicy::KeepAsConfirmed,
            PeriodCombination::Sum,
            Some(frequency.parse()?),
        )?);
        Ok(dates.generate(&spec)?.0)
    };
    let explicit_rolled = explicit_schedule(rolled, "1M")?;
    assert_close(fraction(&explicit_rolled, 0)?, fraction(&rule, 1)?);
    assert_close(fraction(&explicit_rolled, 0)?, 1.0 / 12.0);
    // 確認書が営業日調整後の日付で書いた半年の期間（01-04→07-01 の 178 日。暦の一刻みではない）も、
    // 通常期間として受理し、係数は 1/2 になる。
    let adjusted = explicit_schedule(interval("2027-01-04", "2027-07-01")?, "6M")?;
    assert_close(fraction(&adjusted, 0)?, 0.5);
    Ok(())
}

#[test]
fn period_frame_rejects_reference_periods_that_are_not_one_step() -> Result<()> {
    // 半年払いに一年の参照期間を渡すと、9 か月の期間の ICMA の係数が 3/4 ではなく
    // 273/(2·365) ≈ 0.374 になる。公開の PeriodFrame::new も明示日程も、頻度の一刻みに当たらない
    // （最も近い刻みの数が 1 でない）参照期間を拒否する。半刻みに満たない参照期間も同じく拒否する。
    let semiannual: Frequency = "6M".parse()?;
    let nine_months = interval("2027-01-01", "2027-10-01")?;
    let whole_year = interval("2027-01-01", "2028-01-01")?;
    assert!(matches!(
        PeriodFrame::new(nine_months, vec![whole_year], semiannual),
        Err(Error::InvalidInput(_))
    ));
    let one_month = interval("2027-01-01", "2027-02-01")?;
    assert!(matches!(
        PeriodFrame::new(one_month, vec![one_month], semiannual),
        Err(Error::InvalidInput(_))
    ));
    let explicit_year = ExplicitPeriod {
        references: vec![whole_year],
        ..explicit("2027-01-01", "2027-10-01", "2027-10-01")?
    };
    assert!(matches!(
        ExplicitPeriods::new(
            vec![explicit_year],
            ExplicitDatePolicy::KeepAsConfirmed,
            PeriodCombination::Sum,
            Some(semiannual),
        ),
        Err(Error::InvalidData(_))
    ));
    // 対になる正しい入力: 同じ 9 か月を二つの半年の通常期間で覆えば 1/2 + 92/(2·184) = 3/4。
    let frame = PeriodFrame::new(
        nine_months,
        vec![
            interval("2027-01-01", "2027-07-01")?,
            interval("2027-07-01", "2028-01-01")?,
        ],
        semiannual,
    )?;
    let context = AccrualContext {
        frame: Some(&frame),
        ..AccrualContext::default()
    };
    let (fraction, _) = year_fraction(
        &DayCount::ActActIcma,
        nine_months.start(),
        nine_months.end(),
        &context,
    )?;
    assert_close(fraction.value(), 0.75);
    Ok(())
}

#[test]
fn payment_uses_convention_time_through_dst_and_keeps_distinct_instants() -> Result<()> {
    let dates = Dates::new(&[])?;
    let berlin = event("02:30:00", "Europe/Berlin")?;
    let specs = ExplicitPeriods::new(
        vec![
            ExplicitPeriod {
                settlement: Some(berlin.clone()),
                reset: Some("2027-01-29".parse()?),
                ..explicit("2027-02-01", "2027-02-28", "2027-03-28")?
            },
            ExplicitPeriod {
                settlement: Some(berlin.clone()),
                ..explicit("2027-02-28", "2027-03-15", "2027-03-28")?
            },
            explicit("2027-03-15", "2027-03-28", "2027-03-28")?,
            ExplicitPeriod {
                settlement: Some(berlin),
                ..explicit("2027-03-28", "2027-10-31", "2027-10-31")?
            },
        ],
        ExplicitDatePolicy::KeepAsConfirmed,
        PeriodCombination::Sum,
        None,
    )?;
    let (schedule, warnings) = dates.generate(&ScheduleSpec::Explicit(specs))?;
    let payments = schedule.payments();
    assert_eq!(payments.len(), 3);
    assert_eq!(payments[0].periods(), 0..2);
    assert_eq!(payments[1].periods(), 2..3);
    assert_eq!(payments[0].date(), payments[1].date());
    // 3 月は 02:30 が無く、03:30 CEST に進む。10 月の 02:30 は早い方の CEST を使う。
    assert_eq!(payments[0].at(), "2027-03-28T01:30:00Z".parse()?);
    assert_eq!(payments[1].at(), "2027-03-28T06:00:00Z".parse()?);
    assert_eq!(payments[2].at(), "2027-10-31T00:30:00Z".parse()?);
    assert_eq!(
        warnings,
        [
            Warning::ExplicitDateIsHoliday {
                date: "2027-03-28".parse()?
            },
            Warning::ExplicitDateIsHoliday {
                date: "2027-10-31".parse()?
            },
        ]
    );
    let (resets, _) = reset_dates(
        &schedule,
        &ResetRule {
            anchor: PeriodAnchor::PeriodStart,
            offset: "0BD".parse()?,
            calendar: dates.spec(),
            bdc: Bdc::Following,
        },
        &dates.calendars,
    )?;
    assert_eq!(resets[0], "2027-01-29".parse()?);
    assert_eq!(resets[1], "2027-03-01".parse()?);
    Ok(())
}

#[test]
fn explicit_readjustment_applies_before_settlement_resolution() -> Result<()> {
    let original = Dates::new(&[])?;
    let updated = Dates::new(&["2027-05-06"])?;
    let spec = ScheduleSpec::Explicit(ExplicitPeriods::new(
        vec![explicit("2027-04-06", "2027-05-06", "2027-05-06")?],
        ExplicitDatePolicy::ReadjustWithRule {
            bdc: Bdc::ModifiedFollowing,
            calendar: original.spec(),
        },
        PeriodCombination::Sum,
        Some("1M".parse()?),
    )?);
    // 休日を追加したスナップショットでは支払日を金曜に動かしてから時刻を組み合わせる。
    assert_eq!(
        original.generate(&spec)?.0.payments()[0].at(),
        "2027-05-06T06:00:00Z".parse()?
    );
    let changed = updated.generate(&spec)?.0;
    assert_eq!(changed.payments()[0].at(), "2027-05-07T06:00:00Z".parse()?);
    assert_eq!(changed.periods()[0].accrual().end(), "2027-05-06".parse()?);
    Ok(())
}

#[test]
fn schedule_propagates_calendar_failures_and_deduplicates_warnings() -> Result<()> {
    let mut dates = Dates::new(&[])?;
    let spec = ScheduleSpec::Generated {
        periods: GeneratedPeriods::new(
            StartRule::Date("2027-01-02".parse()?),
            EndRule::Date("2027-01-04".parse()?),
            Frequency::ONCE,
            None,
            StubSpec::NoStub,
            dates.spec(),
            Bdc::Preceding,
            AccrualDates::Adjusted,
        )?,
        payment: PaymentRule {
            bdc: Bdc::Following,
            ..dates.payment()
        },
    };
    let id = dates.id.clone();
    let known = "2027-01-01".parse()?..="2027-01-03".parse()?;
    let weekend = WeekendRule::new([Weekday::Sat, Weekday::Sun]);
    // 計算期間と支払いが同じ範囲外の日を調べても、外挿の警告は 1 件。
    dates.calendars.insert(
        id.clone(),
        BaseCalendar::new(known.clone(), weekend, vec![], Extrapolation::WeekendsOnly)?,
    );
    let (schedule, warnings) = dates.generate(&spec)?;
    assert_eq!(
        schedule.periods()[0].adjusted(),
        interval("2027-01-01", "2027-01-04")?
    );
    assert_eq!(
        warnings,
        [Warning::CalendarExtrapolated {
            calendar: id.clone(),
            date: "2027-01-04".parse()?
        }]
    );
    // 範囲外の扱いを Strict に変えると MissingData、参照先そのものが無ければ InvalidData。
    dates.calendars.insert(
        id,
        BaseCalendar::new(known, weekend, vec![], Extrapolation::Strict)?,
    );
    assert!(matches!(dates.generate(&spec), Err(Error::MissingData(_))));
    dates.calendars.clear();
    assert!(matches!(dates.generate(&spec), Err(Error::InvalidData(_))));
    Ok(())
}

#[test]
fn weekly_schedules_roll_on_the_weekday_and_can_pay_in_advance() -> Result<()> {
    // 金曜に置く週次の格子。2027-01-15（金）を祝日にし、期間の境界は Following で 18 日（月）になる。
    let dates = Dates::new(&["2027-01-15"])?;
    let weekly = |start: &str, step: &str, end: &str, roll| -> Result<GeneratedPeriods> {
        GeneratedPeriods::new(
            StartRule::Date(start.parse()?),
            EndRule::Date(end.parse()?),
            step.parse()?,
            roll,
            StubSpec::NoStub,
            dates.spec(),
            Bdc::Following,
            AccrualDates::Adjusted,
        )
    };
    let friday = Some(RollConvention::OnWeekday(Weekday::Fri));
    let spec = ScheduleSpec::Generated {
        periods: weekly("2027-01-08", "1W", "2027-01-29", friday)?,
        // 期首払い: 調整後の期首の日に払う（前払いの利息や手数料の慣行）。
        payment: PaymentRule {
            anchor: PeriodAnchor::PeriodStart,
            bdc: Bdc::Following,
            ..dates.payment()
        },
    };
    let (schedule, _) = dates.generate(&spec)?;
    let unadjusted: Vec<_> = schedule.periods().iter().map(Period::unadjusted).collect();
    let accrual: Vec<_> = schedule.periods().iter().map(Period::accrual).collect();
    assert_eq!(
        unadjusted,
        [
            interval("2027-01-08", "2027-01-15")?,
            interval("2027-01-15", "2027-01-22")?,
            interval("2027-01-22", "2027-01-29")?,
        ]
    );
    assert_eq!(
        accrual,
        [
            interval("2027-01-08", "2027-01-18")?,
            interval("2027-01-18", "2027-01-22")?,
            interval("2027-01-22", "2027-01-29")?,
        ]
    );
    // 支払いは各期間の調整後の期首で、東京 15 時は UTC 06 時。
    let payments: Vec<_> = schedule
        .payments()
        .iter()
        .map(|payment| (payment.date(), payment.at(), payment.periods()))
        .collect();
    assert_eq!(
        payments,
        [
            ("2027-01-08".parse()?, "2027-01-08T06:00:00Z".parse()?, 0..1),
            ("2027-01-18".parse()?, "2027-01-18T06:00:00Z".parse()?, 1..2),
            ("2027-01-22".parse()?, "2027-01-22T06:00:00Z".parse()?, 2..3),
        ]
    );

    // 4 週の刻みは 28 日ごとに進む。ロールを省けば起点の曜日（月曜）のまま並ぶ。
    let spec = ScheduleSpec::Generated {
        periods: weekly("2027-01-04", "4W", "2027-03-29", None)?,
        payment: dates.payment(),
    };
    let (schedule, _) = dates.generate(&spec)?;
    let accrual: Vec<_> = schedule.periods().iter().map(Period::accrual).collect();
    assert_eq!(
        accrual,
        [
            interval("2027-01-04", "2027-02-01")?,
            interval("2027-02-01", "2027-03-01")?,
            interval("2027-03-01", "2027-03-29")?,
        ]
    );
    assert!(
        schedule
            .periods()
            .iter()
            .all(|period| period.kind() == PeriodKind::Regular)
    );

    // 曜日のロールは週の刻みだけに使え、格子の起点もその曜日でなければならない。
    assert!(matches!(
        weekly("2027-01-08", "1M", "2027-04-08", friday),
        Err(Error::InvalidData(_))
    ));
    let thursday = ScheduleSpec::Generated {
        periods: weekly("2027-01-07", "1W", "2027-01-28", friday)?,
        payment: dates.payment(),
    };
    assert!(matches!(
        dates.generate(&thursday),
        Err(Error::InvalidData(_))
    ));
    Ok(())
}
`````

`````
crates/rhodium-core/tests/schedules/observation.rs
//! RFR などの観測点・重み・年率換算の係数を、暦と日数基準から手で数えた値と照合する。

use std::num::NonZeroU8;

use rhodium_core::{
    Error, Result,
    calendar::BusinessDayConvention as Bdc,
    day_count::{AccrualContext, DayCount, year_fraction},
    fixing::{ObservationPoint, YearMonth},
    observation::{
        Aggregation, AveragingWeights, CompoundingForm, EndpointDates, ObservationDays,
        ObservationMethod, ObservationRounding, ObservationRule, ObservationSchedule,
        SpreadTreatment, observation_schedule,
    },
    rounding::{Rounding, RoundingMode},
    schedule::{Period, Schedule, StubSpec},
    tenor::Frequency,
};

use super::{Dates, assert_close, generated, interval};

impl Dates {
    fn once(&self, start: &str, end: &str) -> Result<Schedule> {
        self.generate(&generated(
            self,
            interval(start, end)?,
            Frequency::ONCE,
            None,
            StubSpec::NoStub,
            Bdc::Unadjusted,
        )?)
        .map(|(schedule, _)| schedule)
    }
}

const COMPOUND: Aggregation = Aggregation::Compound {
    form: CompoundingForm::SimpleFactor,
    spread: SpreadTreatment::Exclusive,
};
const WEIGHTED: Aggregation = Aggregation::Average(AveragingWeights::CoverWeighted);

fn compound(dates: &Dates, method: ObservationMethod) -> Result<ObservationRule> {
    ObservationRule::new(
        ObservationDays::BusinessDays(dates.spec()),
        method,
        COMPOUND,
        ObservationRounding::default(),
    )
}

fn observe(dates: &Dates, period: &Period, rule: &ObservationRule) -> Result<ObservationSchedule> {
    observation_schedule(
        period,
        rule,
        &DayCount::Act365F,
        &AccrualContext::default(),
        &dates.calendars,
    )
    .map(|(value, _)| value)
}

fn check_daily(schedule: &ObservationSchedule, days: &[(&str, f64)]) -> Result<()> {
    assert_eq!(schedule.fixings().len(), days.len());
    for (&(point, weight), &(date, days)) in schedule.fixings().iter().zip(days) {
        assert_eq!(point, ObservationPoint::Date(date.parse()?));
        assert_close(weight, days / 365.0);
    }
    Ok(())
}

#[test]
fn rfr_methods_move_fixing_dates_and_weights_independently() -> Result<()> {
    use ObservationMethod::{Lockout, Lookback, ObservationShift, Plain};
    let dates = Dates::new(&[])?;
    let schedule = dates.once("2027-01-08", "2027-01-13")?;
    let period = &schedule.periods()[0];
    let lag = NonZeroU8::MIN;
    // 金曜の率がそのままなら週末を含む 3 日を覆い、lookback では木曜の値だけが移る。
    for (method, weights) in [
        (
            Plain,
            [
                ("2027-01-08", 3.0),
                ("2027-01-11", 1.0),
                ("2027-01-12", 1.0),
            ],
        ),
        (
            Lookback(lag),
            [
                ("2027-01-07", 3.0),
                ("2027-01-08", 1.0),
                ("2027-01-11", 1.0),
            ],
        ),
        (
            ObservationShift(lag),
            [
                ("2027-01-07", 1.0),
                ("2027-01-08", 3.0),
                ("2027-01-11", 1.0),
            ],
        ),
        (
            Lockout(lag),
            [
                ("2027-01-08", 3.0),
                ("2027-01-11", 1.0),
                ("2027-01-11", 1.0),
            ],
        ),
    ] {
        let observations = observe(&dates, period, &compound(&dates, method)?)?;
        check_daily(&observations, &weights)?;
        assert_close(observations.annualization().value(), 5.0 / 365.0);
    }
    let weekend = dates.once("2027-01-08", "2027-01-11")?;
    let period = &weekend.periods()[0];
    let observations = observe(&dates, period, &compound(&dates, ObservationShift(lag))?)?;
    check_daily(&observations, &[("2027-01-07", 1.0)])?;
    let (coupon, _) = year_fraction(
        &DayCount::Act365F,
        period.accrual().start(),
        period.accrual().end(),
        &AccrualContext::default(),
    )?;
    // 観測は木曜から金曜までの 1 日、利息計算は金曜から月曜までの 3 日。
    assert_close(observations.annualization().value(), 1.0 / 365.0);
    assert_close(coupon.value(), 3.0 / 365.0);
    assert!(matches!(
        observe(&dates, period, &compound(&dates, Lockout(lag))?),
        Err(Error::InvalidData(_))
    ));
    Ok(())
}

#[test]
fn rfr_period_starting_on_holiday() -> Result<()> {
    let dates = Dates::new(&[])?;
    let schedule = dates.once("2027-01-09", "2027-01-12")?;
    let period = &schedule.periods()[0];
    // 土曜から月曜までを直前の金曜、月曜から火曜までを月曜の値で覆う。
    let observations = observe(&dates, period, &compound(&dates, ObservationMethod::Plain)?)?;
    check_daily(&observations, &[("2027-01-08", 2.0), ("2027-01-11", 1.0)])?;
    let coverage = [
        interval("2027-01-09", "2027-01-11")?,
        interval("2027-01-11", "2027-01-12")?,
    ];
    assert_eq!(observations.coverage(), &coverage);
    assert_close(observations.annualization().value(), 3.0 / 365.0);
    let lag = NonZeroU8::MIN;
    let lookback = observe(
        &dates,
        period,
        &compound(&dates, ObservationMethod::Lookback(lag))?,
    )?;
    check_daily(&lookback, &[("2027-01-07", 2.0), ("2027-01-08", 1.0)])?;
    assert_eq!(lookback.coverage(), &coverage);
    // 金曜は計算期間より前の補完なので、締めの日を決める観測数に含めない。
    assert!(matches!(
        observe(
            &dates,
            period,
            &compound(&dates, ObservationMethod::Lockout(lag))?
        ),
        Err(Error::InvalidData(_))
    ));
    let weekend = dates.once("2027-01-09", "2027-01-11")?;
    let observations = observe(
        &dates,
        &weekend.periods()[0],
        &compound(&dates, ObservationMethod::Plain)?,
    )?;
    check_daily(&observations, &[("2027-01-08", 2.0)])?;
    Ok(())
}

#[test]
fn observation_weights_use_the_passed_snapshot_and_day_count() -> Result<()> {
    let original = Dates::new(&["2027-05-03", "2027-05-04", "2027-05-05"])?;
    let updated = Dates::new(&["2027-05-03", "2027-05-04", "2027-05-05", "2027-05-06"])?;
    let schedule = original.once("2027-04-30", "2027-05-10")?;
    let p = &schedule.periods()[0];
    let rule = compound(&original, ObservationMethod::Plain)?;
    // 5 月 6 日を休日にすると 4 月 30 日の観測の重みは 6 日から 7 日になり、6 日の観測点が消える。
    check_daily(
        &observe(&original, p, &rule)?,
        &[
            ("2027-04-30", 6.0),
            ("2027-05-06", 1.0),
            ("2027-05-07", 3.0),
        ],
    )?;
    check_daily(
        &observe(&updated, p, &rule)?,
        &[("2027-04-30", 7.0), ("2027-05-07", 3.0)],
    )?;
    let basis = DayCount::Bus252 {
        calendar: updated.id.clone(),
    };
    let (observations, _) = observation_schedule(
        p,
        &rule,
        &basis,
        &AccrualContext::default(),
