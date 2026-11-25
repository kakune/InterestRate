pub(crate) const FINAL: &str = "2026-01-07";
pub(crate) const PILLAR: &str = "2026-05-08T16:00:00Z";
pub(crate) const CLEARING: &str = "USD-FUT-JAN";
pub(crate) const MULTIPLIER: f64 = 1_000.0;
pub(crate) const RATE: f64 = 0.04;

pub(crate) fn at(day: &str, time: &str) -> Result<Instant> {
    format!("{day}T{time}:00Z").parse()
}

pub(crate) fn refdata(delay: DateOffset, settlement_time: &str) -> Result<RefData> {
    RefData::new(reference_parts(delay, settlement_time)?)
}

#[expect(
    clippy::too_many_lines,
    reason = "暦、原指標、清算と限月のテスト用参照データをここにまとめる"
)]
pub(crate) fn reference_parts(delay: DateOffset, settlement_time: &str) -> Result<RefDataParts> {
    let calendar: CalendarId = "USD".parse()?;
    let zone: TzId = "Etc/UTC".parse()?;
    let time = "16:00:00".parse()?;
    let clock = LocalEventTime {
        time,
        zone: zone.clone(),
    };
    let mut parts = RefDataParts::default();
    parts.zones.insert(
        zone.clone(),
        ZoneRules::new(
            "2025-01-01T00:00:00Z".parse()?..="2027-01-01T00:00:00Z".parse()?,
            UtcOffset::from_seconds(0)?,
            vec![],
        )?,
    );
    parts.calendars.insert(
        calendar.clone(),
        BaseCalendar::new(
            "2025-01-01".parse()?..="2026-12-31".parse()?,
            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
            vec![],
            CalendarExtrapolation::Strict,
        )?,
    );
    parts.currencies.insert(
        Currency::USD,
        CurrencyDefinition {
            code: Currency::USD,
            minor_units: 2,
            payment_rounding: Rounding::new(0.01, RoundingMode::HalfEven)?,
            deliverability: Deliverability::Deliverable,
            settlement_time: clock.clone(),
            settlement_calendar: calendar.clone(),
        },
    );
    let security: SecurityId = CLEARING.parse()?;
    parts.securities.insert(
        security.clone(),
        SecurityDefinition {
            id: security.clone(),
            kind: SecurityKind::ListedFuture,
            currency: Currency::USD,
            unit: SecurityUnit::Contract,
            issuer: None,
            valuation: AssetValuation::Forward,
        },
    );
    let calendar = CalendarSpec::Base(calendar);
    for (name, kind) in [
        (
            "TERM",
            IndexKind::TermRate {
                ccy: Currency::USD,
                day_count: DayCount::Act365F,
                tenor: "1W".parse()?,
                value_offset: "0D".parse()?,
                bdc: BusinessDayConvention::Following,
                end_of_month: EndOfMonthRule::Ignore,
            },
        ),
        (
            CLEARING,
            IndexKind::RateFuturesSettlement {
                contract: security,
                underlying: "TERM".parse()?,
            },
        ),
    ] {
        let id: IndexId = name.parse()?;
        parts.indices.insert(
            id.clone(),
            IndexDefinition {
                id,
                kind,
                calendar: calendar.clone(),
                status: IndexStatus::Active,
                observation_time: Some(clock.clone()),
                publication: PublicationRule::Daily(DailyPublication {
                    business_days: 0,
                    calendar: calendar.clone(),
                    time,
                    zone: zone.clone(),
                }),
            },
        );
    }
    parts.conventions.insert(
        "USD-RATE-FUTURES".parse()?,
        Convention::RateFuture(Box::new(RateFutureConvention {
            multiplier: MULTIPLIER,
            settlement_lag: delay,
            settlement_time: LocalEventTime {
                time: settlement_time.parse()?,
                zone,
            },
            margin_index: None,
            contracts: BTreeMap::from([(
                "2026-01".parse()?,
                RateFutureContract {
                    settlement_index: CLEARING.parse()?,
                    final_clearing: FINAL.parse()?,
                    rate: RateFutureRate::Term {
                        fixing: FINAL.parse()?,
                    },
                },
            )]),
        })),
    );
    Ok(parts)
}

pub(crate) fn config(extra: &str) -> Result<PricingConfig> {
    toml::from_str(&format!(
        r#"
version = 1
reporting_currency = "USD"
[[discounting]]
currency = "USD"
collateral = "uncollateralized"
curve = "USD"
[forecasting]
TERM = "USD"
{extra}
"#
    ))
    .map_err(|error| Error::InvalidData(error.to_string()))
}

pub(crate) fn curve() -> Result<CurveSpec> {
    CurveSpec::new(
        "USD".parse()?,
        Currency::USD,
        [CurveRole::Discount, CurveRole::Projection],
        [PILLAR.parse()?],
        Extrapolation::Reject,
    )
}

pub(crate) fn market(as_of: Instant, volatility: f64) -> Result<Market<f64>> {
    let structure = Arc::new(MarketStructure::with_fx_and_lgm(
        as_of,
        [curve()?],
        [],
        [(Currency::USD, vec![])],
    )?);
    let curve = structure
        .layout()
        .range(&BlockId::nodes("USD".parse()?))
        .ok_or_else(|| Error::Internal("missing curve block".into()))?;
    let model = structure
        .layout()
        .range(&MarketStructure::lgm_block(Currency::USD)?)
        .ok_or_else(|| Error::Internal("missing model block".into()))?;
    let mut theta = vec![0.0; 3];
    theta[curve.start] = -RATE * model_time(as_of, PILLAR.parse()?).years();
    theta[model.start + 1] = volatility;
    Market::from_theta(structure, theta)
}

/// 2026 年 1 月限の金利先物の校正テンプレート。
pub(crate) fn futures_template() -> Result<InstrumentTemplate> {
    Ok(InstrumentTemplate::RateFuture {
        convention: "USD-RATE-FUTURES".parse()?,
        delivery_month: "2026-01".parse()?,
    })
}

/// USD 曲線の一つのピラーを先物の価格で校正し、凸性の補正に使う LGM の平均回帰とボラを
/// 直接マークにする市場の定義。
pub(crate) fn futures_definition() -> Result<MarketDefinition> {
    let binding = InstrumentBinding::new(
        "USD.FUT.2026-01".parse()?,
        QuoteKind::FuturesPrice,
        CalibrationTarget::Price(futures_template()?),
        QuoteTransform::Identity,
        CollateralKey::Uncollateralized,
        ElementSelector::Block(BlockId::nodes("USD".parse()?)),
    )?;
    MarketDefinition::new("USD".parse()?, "Etc/UTC".parse()?, [curve()?], [binding])?.with_lgm([(
        Currency::USD,
        LgmDefinition {
            breaks: vec![],
            mean_reversion: vec!["USD.LGM.MR".parse()?],
            volatility: vec!["USD.LGM.VOL".parse()?],
        },
    )])
}
`````

`````
crates/rhodium-pricing/tests/fixtures/multicurve.rs
//! OIS とターム指標のベーシスを校正する共通のテスト入力。

use std::fmt::Write as _;

use rhodium_core::{
    Error, Result,
    calendar::{BusinessDayConvention, EndOfMonthRule},
    currency::Currency,
    refdata::{
        Convention, DepositConvention, IndexKind, SwapConvention, SwapCouponConvention,
        SwapLegConvention,
    },
    tenor::{DateOffset, TenorUnit},
};
use rhodium_market::{
    curve::{CurveExpr, CurveRole, CurveSpec, Extrapolation},
    market::Market,
    quote::{QuoteId, QuoteKind, QuoteSet},
    structure::BlockId,
};
use rhodium_pricing::{
    context::Context,
    definition::{
        CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
    },
    prepare::{PrepareInput, Prepared, prepare},
    resolve::CollateralKey,
};
use rhodium_products::standard::InstrumentTemplate;
use serde::Deserialize;

use super::ois::{Case, Inputs, Node, Sample, cases, fixture_quotes, midnight};

const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;

#[derive(Deserialize)]
struct QuantLibBasis {
    overnight_nodes: Vec<Node>,
    overnight_samples: Vec<Sample>,
    term_index: String,
    term_tenor: String,
    nodes: Vec<Node>,
    samples: Vec<Sample>,
}

/// OIS とターム指標の校正商品、ピラーと比較点。
pub struct MulticurveCase {
    /// 割引と翌日物の予想を行う OIS の入力。
    pub overnight: Case,
    /// ターム指標とその予想曲線の ID。
    pub term_index: String,
    /// ターム指標のテナー（ベーシスの両脚の支払いの頻度）。
    pub term_tenor: String,
    /// ベーシス商品のテナーと spread、ターム予想曲線のピラー。
    pub nodes: Vec<Node>,
    /// ターム予想曲線のピラー以外の比較点。
    pub samples: Vec<Sample>,
}

/// JPY の TONA/TIBOR 6M と USD の SOFR/テスト用ターム指標の入力を読む。
pub(crate) fn multicurve_cases() -> Result<[MulticurveCase; 2]> {
    let [jpy, usd, _] = all_multicurve_cases()?;
    Ok([jpy, usd])
}

/// TIBOR 6M、USD のターム指標、TIBOR 3M の三つのベーシスの入力を読む。
pub(crate) fn all_multicurve_cases() -> Result<[MulticurveCase; 3]> {
    let [jpy, usd, jpy_3m]: [QuantLibBasis; 3] =
        serde_json::from_str(include_str!("quantlib/multicurve.json"))
            .map_err(|error| Error::InvalidData(error.to_string()))?;
    let [tona, sofr] = cases()?;
    Ok([
        jpy.combine(tona.clone()),
        usd.combine(sofr),
        jpy_3m.combine(tona),
    ])
}

impl QuantLibBasis {
    fn combine(self, mut overnight: Case) -> MulticurveCase {
        overnight.nodes = self.overnight_nodes;
        overnight.samples = self.overnight_samples;
        overnight.swaps.clear();
        MulticurveCase {
            overnight,
            term_index: self.term_index,
            term_tenor: self.term_tenor,
            nodes: self.nodes,
            samples: self.samples,
        }
    }
}

impl MulticurveCase {
    /// OIS にターム指標とベーシス商品の慣行を加える。
    pub(crate) fn inputs(&self) -> Result<Inputs> {
        self.inputs_with("")
    }

    /// [`Self::inputs`] の価格設定に、割引の行や近似の設定を加える。
    pub(crate) fn inputs_with(&self, extra_config: &str) -> Result<Inputs> {
        Inputs::with_forecasts(
            &[&self.overnight],
            self.overnight.as_of,
            extra_config,
            |parts, forecasting| {
                let overnight_id = self.overnight.index.parse()?;
                let term_id = self.term_index.parse()?;
                let Some(Convention::Ois(ois)) = parts.conventions.get(&overnight_id) else {
                    return Err(Error::Internal("missing OIS test convention".into()));
                };
                let mut term = parts.indices[&ois.float.index].clone();
                term.id = term_id;
                term.kind = IndexKind::TermRate {
                    ccy: self.overnight.currency,
                    day_count: self.overnight.basis(),
                    tenor: self.term_tenor.parse()?,
                    value_offset: DateOffset::new(2, TenorUnit::BusinessDay),
                    bdc: BusinessDayConvention::ModifiedFollowing,
                    end_of_month: EndOfMonthRule::Ignore,
                };
                let mut pay = SwapLegConvention {
                    currency: self.overnight.currency,
                    periods: ois.fixed.clone(),
                    payment: ois.payment.clone(),
                    coupon: SwapCouponConvention::Compounded(ois.float.clone()),
                    fx_reset: None,
                };
                // QuantLib のヘルパーと同じく、両脚をターム指標のテナーごとにラグなしで払う。
                pay.periods.frequency = self.term_tenor.parse()?;
                pay.payment.offset = DateOffset::new(0, TenorUnit::BusinessDay);
                let mut receive = pay.clone();
                receive.coupon = SwapCouponConvention::Term {
                    index: term.id.clone(),
                };
                // 固定対ターム金利の IRS（固定を払い、ターム金利を受ける。両脚はターム指標のテナーごと）。
                let mut fixed = pay.clone();
                fixed.coupon = SwapCouponConvention::Fixed;
                let spot_lag = ois.spot_lag;
                // スポット起点の預金（TONA の暦と日数基準）。
                let deposit = DepositConvention {
                    spot_lag,
                    day_count: self.overnight.basis(),
                    bdc: BusinessDayConvention::ModifiedFollowing,
                    end_of_month: EndOfMonthRule::Ignore,
                    calendar: ois.fixed.calendar.clone(),
                };
                parts.indices.insert(term.id.clone(), term);
                parts.conventions.insert(
                    self.term_index.parse()?,
                    Convention::BasisSwap(Box::new(SwapConvention {
                        spot_lag,
                        legs: [pay, receive.clone()],
                    })),
                );
                parts.conventions.insert(
                    self.irs_convention().parse()?,
                    Convention::Irs(Box::new(SwapConvention {
                        spot_lag,
                        legs: [fixed, receive],
                    })),
                );
                parts.conventions.insert(
                    self.deposit_convention().parse()?,
                    Convention::Deposit(deposit),
                );
                writeln!(
                    forecasting,
                    "\"{}\" = \"{}\"",
                    self.term_index, self.term_index
                )
                .map_err(|error| Error::Internal(error.to_string()))
            },
        )
    }

    /// 固定対ターム金利の IRS の慣行の ID。
    pub(crate) fn irs_convention(&self) -> String {
        format!("{}-IRS", self.term_index)
    }

    /// スポット起点の預金の慣行の ID。
    pub(crate) fn deposit_convention(&self) -> String {
        format!("{}-DEPO", self.overnight.index)
    }

    /// ターム予想曲線の校正対象。
    pub(crate) fn term_block(&self) -> Result<BlockId> {
        BlockId::new(self.term_index.parse()?, "basis")
    }

    /// ベーシス商品の quote ID。
    pub(crate) fn quote(&self, node: &Node) -> Result<QuoteId> {
        format!("{}.BASIS.{}", self.term_index, node.tenor).parse()
    }

    pub(crate) fn template(&self, node: &Node) -> Result<InstrumentTemplate> {
        Ok(InstrumentTemplate::BasisSwap {
            convention: self.term_index.parse()?,
            tenor: node.tenor,
        })
    }

    /// 共通の OIS 曲線と、それを上流に持つターム予想曲線を定義する。
    pub(crate) fn definition(&self) -> Result<MarketDefinition> {
        self.definition_with([])
    }

    /// [`Self::definition`] に、校正しない曲線（担保条件の違う割引曲線など）を加える。
    pub(crate) fn definition_with(
        &self,
        extra: impl IntoIterator<Item = CurveSpec>,
    ) -> Result<MarketDefinition> {
        self.definition_with_collateral(extra, |_| None)
    }

    /// 共通の OIS 曲線に、ターム予想曲線のブロックを当てる任意の binding を組み合わせる。
    pub(crate) fn definition_with_term_bindings(
        &self,
        term_bindings: impl IntoIterator<Item = InstrumentBinding>,
    ) -> Result<MarketDefinition> {
        let overnight_block = ElementSelector::Block(self.overnight.block()?);
        let bindings = term_bindings
            .into_iter()
            .map(Ok)
            .chain(self.overnight.nodes.iter().map(|node| {
                self.overnight
                    .binding(node, &COLLATERAL, overnight_block.clone())
            }))
            .collect::<Result<Vec<_>>>()?;
        MarketDefinition::new(
            self.overnight.currency.to_string().parse()?,
            "Etc/UTC".parse()?,
            [self.term_curve()?, self.overnight.curve()?],
            bindings,
        )
    }

    /// ターム予想曲線。JPY は OIS 曲線とベーシスの和、USD は独立した節点の曲線。
    pub(crate) fn term_curve(&self) -> Result<CurveSpec> {
        let term_nodes = CurveExpr::Nodes {
            part: "basis".into(),
            pillars: self
                .nodes
                .iter()
                .map(|node| midnight(node.date))
                .collect::<Result<_>>()?,
            extrapolation: Extrapolation::Reject,
        };
        let term_curve = match self.overnight.currency {
            Currency::JPY => CurveExpr::Sum(vec![
                CurveExpr::Ref(self.overnight.index.parse()?),
                term_nodes,
            ]),
            _ => term_nodes,
        };
        CurveSpec::from_expr(
            self.term_index.parse()?,
            self.overnight.currency,
            [CurveRole::Projection],
            term_curve,
        )
    }

    /// [`Self::definition_with`] で、`collateral` が返すベーシスの行（ノードの番号）だけ
    /// 校正商品の担保を変える。
    pub(crate) fn definition_with_collateral(
        &self,
        extra: impl IntoIterator<Item = CurveSpec>,
        collateral: impl Fn(usize) -> Option<CollateralKey>,
    ) -> Result<MarketDefinition> {
        let overnight = self.overnight.curve()?;
        let term_nodes = CurveExpr::Nodes {
            part: "basis".into(),
            pillars: self
                .nodes
                .iter()
                .map(|node| midnight(node.date))
                .collect::<Result<_>>()?,
            extrapolation: Extrapolation::Reject,
        };
        let term_curve = match self.overnight.currency {
            Currency::JPY => {
                CurveExpr::Sum(vec![CurveExpr::Ref(overnight.id().clone()), term_nodes])
            }
            _ => term_nodes,
        };
        let term = CurveSpec::from_expr(
            self.term_index.parse()?,
            self.overnight.currency,
            [CurveRole::Projection],
            term_curve,
        )?;
        let overnight_block = ElementSelector::Block(self.overnight.block()?);
        let term_block = ElementSelector::Block(self.term_block()?);
        // ベーシスの行を先に渡しても、足跡の依存関係で OIS が先に解かれる。
        let bindings = self
            .nodes
            .iter()
            .enumerate()
            .map(|(index, node)| {
                InstrumentBinding::new(
                    self.quote(node)?,
                    QuoteKind::Spread,
                    CalibrationTarget::ParRate(self.template(node)?),
                    QuoteTransform::Identity,
                    collateral(index).unwrap_or(COLLATERAL),
                    term_block.clone(),
                )
            })
            .chain(self.overnight.nodes.iter().map(|node| {
                self.overnight
                    .binding(node, &COLLATERAL, overnight_block.clone())
            }))
            .collect::<Result<Vec<_>>>()?;
        MarketDefinition::new(
            self.overnight.currency.to_string().parse()?,
            "Etc/UTC".parse()?,
            [term, overnight].into_iter().chain(extra),
            bindings,
        )
    }

    /// OIS とベーシスの入力 quote を作る。
    pub(crate) fn quotes(&self) -> Result<QuoteSet> {
        let term = self
            .nodes
            .iter()
            .map(|node| Ok((self.quote(node)?, node.rate)));
        QuoteSet::new(
            fixture_quotes(&self.overnight)?
                .into_iter()
                .chain(term.collect::<Result<Vec<_>>>()?),
        )
    }

    /// 指定のテナーと spread のベーシス取引を準備する。
    pub(crate) fn prepare_basis(
        &self,
        node: &Node,
        spread: f64,
        market: &Market<f64>,
        ctx: &Context<'_>,
    ) -> Result<Prepared> {
        let product = self.template(node)?.build(spread, ctx.refdata)?;
        prepare(
            PrepareInput::Standard {
                product: &product,
                trade_date: self.overnight.as_of,
            },
            &COLLATERAL,
            &[],
            market.structure(),
            ctx,
        )
    }
}
`````

`````
crates/rhodium-pricing/tests/fixtures/ois.rs
//! TONA と SOFR の校正テストで共用する市場入力。

use std::{collections::BTreeMap, fmt::Write as _};

use rhodium_core::{
    Error, Result,
    calendar::{
        BaseCalendar, BusinessDayConvention, CalendarSpec, Extrapolation as CalendarExtrapolation,
        WeekendRule,
    },
    currency::Currency,
    date::{Date, Weekday},
    day_count::DayCount,
    fixing::{FixingStore, FixingStoreParts, ObservationPoint},
    fx::CcyPair,
    id::{CalendarId, IndexId},
    observation::{
        Aggregation, CompoundingForm, ObservationDays, ObservationMethod, ObservationRounding,
        ObservationRule, SpreadTreatment,
    },
    refdata::{
        Convention, CurrencyDefinition, DailyPublication, Deliverability, FixedLegConvention,
        FxConvention, FxSpotLag, IndexDefinition, IndexKind, IndexStatus, OisConvention,
        OvernightLegConvention, PublicationRule, RefData, RefDataParts,
    },
    rounding::{Rounding, RoundingMode},
    schedule::{PaymentGrouping, PaymentRule, PeriodAnchor},
    tenor::{DateOffset, Frequency, TenorLength, TenorUnit},
    time::{Instant, LocalEventTime, UtcOffset, ZoneRules},
};
use rhodium_market::{
    curve::{CurveRole, CurveSpec, Extrapolation},
    market::Market,
    quote::{QuoteId, QuoteKind},
    structure::{BlockId, ParamId},
};
use rhodium_pricing::{
    config::PricingConfig,
    context::Context,
    definition::{
        CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
    },
    prepare::{PrepareInput, Prepared, prepare},
    resolve::CollateralKey,
};
use rhodium_products::standard::{InstrumentTemplate, PayReceive};
use serde::Deserialize;

#[derive(Clone, Deserialize)]
/// 一つの OIS 曲線と参照データのテスト入力。
pub struct Case {
    /// 曲線の通貨。
    pub currency: Currency,
    /// 指標、曲線、OIS 慣行に共通の ID。
    pub index: String,
    /// 評価日。
    pub as_of: Date,
    calendar_from: Date,
    calendar_through: Date,
    holidays: Vec<Date>,
    payment_lag: i16,
    /// テナー順のピラー。
    pub nodes: Vec<Node>,
    /// 評価日から最後のピラーまでの各区間に置いた、補間を照合する点。
    pub samples: Vec<Sample>,
    /// 評価日の前日までに公表された fixing。開始済みのスワップが読む。
    pub fixings: Vec<FixingRow>,
    /// 同じ曲線で `QuantLib` が評価した、慣行どおりの OIS。
    pub swaps: Vec<QuantLibOis>,
    /// OIS 慣行の観測の方法。JSON には無く、既定は `Plain`（`QuantLib` と同じ）。
    #[serde(skip)]
    pub observation_method: ObservationMethod,
    /// OIS 慣行の観測金利の丸め。JSON には無く、既定は丸めなし。
    #[serde(skip)]
    pub observation_rounding: ObservationRounding,
}

#[derive(Clone, Deserialize)]
/// ピラー以外の日付の `QuantLib` の割引係数。
pub struct Sample {
    /// 照合する日付。
    pub date: Date,
    /// その日付の割引係数。
    pub df: f64,
}

#[derive(Clone, Deserialize)]
/// 公表済みの日次の fixing。
pub struct FixingRow {
    /// 観測日。
    pub date: Date,
    /// 公表値。
    pub value: f64,
}

#[derive(Clone, Deserialize)]
/// `QuantLib` の `OvernightIndexedSwap` で評価したスワップと脚別の PV。
pub struct QuantLibOis {
    /// スポット起点（`spot`）か開始済み（`seasoned`）か。
    pub name: String,
    /// 計算期間の開始日。
    pub start: Date,
    /// 開始日から満期までの長さ。
    pub tenor: TenorLength,
    /// 固定金利。
    pub fixed_rate: f64,
    /// 複利の外で足す変動脚のスプレッド。
    pub spread: f64,
    /// 元本。
    pub notional: f64,
    /// 固定脚を払うか受けるか。
    pub fixed_side: PayReceive,
    /// 自社側の PV。
    pub pv: f64,
    /// 固定脚の PV。
    pub fixed_pv: f64,
    /// 変動脚の PV。
    pub floating_pv: f64,
}

#[derive(Clone, Deserialize)]
/// 一つの OIS 校正商品の入力とピラー。
pub struct Node {
    /// 校正商品のテナー。
    pub tenor: TenorLength,
    /// パーレート。
    pub rate: f64,
    /// 曲線のピラーの日付。
    pub date: Date,
    /// ピラーでの `QuantLib` の割引係数。
    pub df: f64,
}

/// TONA と SOFR のテスト入力を読み込む。
pub(crate) fn cases() -> Result<[Case; 2]> {
    serde_json::from_str(include_str!("quantlib/ois_curves.json"))
        .map_err(|error| Error::InvalidData(error.to_string()))
}

/// 日付を UTC の 00:00 に解決する。
pub(crate) fn midnight(date: Date) -> Result<Instant> {
    format!("{date}T00:00:00Z").parse()
}

impl Case {
    /// 固定脚と観測に使う日数基準。
    pub(crate) fn basis(&self) -> DayCount {
        if self.currency == Currency::JPY {
            DayCount::Act365F
        } else {
            DayCount::Act360
        }
    }

    /// 曲線の節点ブロック。
    pub(crate) fn block(&self) -> Result<BlockId> {
        Ok(BlockId::nodes(self.index.parse()?))
    }

    /// 割引と予想に使える曲線。
    pub(crate) fn curve(&self) -> Result<CurveSpec> {
        CurveSpec::new(
            self.index.parse()?,
            self.currency,
            [CurveRole::Discount, CurveRole::Projection],
            self.nodes
                .iter()
                .map(|node| midnight(node.date))
                .collect::<Result<Vec<_>>>()?,
            Extrapolation::Reject,
        )
    }

    /// ピラーに対応する慣行ベースの OIS。
    pub(crate) fn template(&self, node: &Node) -> Result<InstrumentTemplate> {
        Ok(InstrumentTemplate::Ois {
            convention: self.index.parse()?,
            tenor: node.tenor,
        })
    }

    /// ピラーの quote ID。
    pub(crate) fn quote(&self, node: &Node) -> Result<QuoteId> {
        format!("{}.OIS.{}", self.index, node.tenor).parse()
    }

    /// ピラーの校正条件。
    pub(crate) fn binding(
        &self,
        node: &Node,
        collateral: &CollateralKey,
        determines: ElementSelector,
    ) -> Result<InstrumentBinding> {
        InstrumentBinding::new(
            self.quote(node)?,
            QuoteKind::Rate,
            CalibrationTarget::ParRate(self.template(node)?),
            QuoteTransform::Identity,
            collateral.clone(),
            determines,
        )
    }

    fn add_refdata(&self, parts: &mut RefDataParts, clock: &LocalEventTime) -> Result<()> {
        let calendar_id: CalendarId = self.currency.to_string().parse()?;
        let calendar = CalendarSpec::Base(calendar_id.clone());
        parts.calendars.insert(
            calendar_id.clone(),
            BaseCalendar::new(
                self.calendar_from..=self.calendar_through,
                WeekendRule::new([Weekday::Sat, Weekday::Sun]),
                self.holidays.clone(),
                CalendarExtrapolation::Strict,
            )?,
        );
        let minor_units = if self.currency == Currency::JPY { 0 } else { 2 };
        parts.currencies.insert(
            self.currency,
            CurrencyDefinition {
                code: self.currency,
                minor_units,
                payment_rounding: Rounding::new(
                    10.0_f64.powi(-i32::from(minor_units)),
                    RoundingMode::HalfEven,
                )?,
                deliverability: Deliverability::Deliverable,
                settlement_time: clock.clone(),
                settlement_calendar: calendar_id.clone(),
            },
        );
        let index: IndexId = self.index.parse()?;
        parts.indices.insert(
            index.clone(),
            IndexDefinition {
                id: index.clone(),
                kind: IndexKind::OvernightRate {
                    ccy: self.currency,
                    day_count: self.basis(),
                },
                calendar: calendar.clone(),
                publication: PublicationRule::Daily(DailyPublication {
                    business_days: 1,
                    calendar: calendar.clone(),
                    time: clock.time,
                    zone: clock.zone.clone(),
                }),
                observation_time: Some(clock.clone()),
                status: IndexStatus::Active,
            },
        );
        let observation = ObservationRule::new(
            ObservationDays::BusinessDays(calendar.clone()),
            self.observation_method,
            Aggregation::Compound {
                form: CompoundingForm::SimpleFactor,
                spread: SpreadTreatment::Exclusive,
            },
            self.observation_rounding,
        )?;
        parts.conventions.insert(
            self.index.parse()?,
            Convention::Ois(Box::new(OisConvention {
                spot_lag: DateOffset::new(2, TenorUnit::BusinessDay),
                fixed: FixedLegConvention {
                    frequency: Frequency::ANNUAL,
                    day_count: self.basis(),
                    bdc: BusinessDayConvention::ModifiedFollowing,
                    calendar: calendar.clone(),
                },
                float: OvernightLegConvention { index, observation },
                payment: PaymentRule {
                    anchor: PeriodAnchor::PeriodEnd,
                    offset: DateOffset::new(self.payment_lag, TenorUnit::BusinessDay),
                    calendar,
                    bdc: BusinessDayConvention::Following,
                    grouping: PaymentGrouping::default(),
                    settlement: None,
                },
            })),
        );
        Ok(())
    }
}

/// テストの参照データ、fixing と価格設定。
pub struct Inputs {
    /// 評価の瞬間。
    pub as_of: Instant,
    refdata: RefData,
    fixings: FixingStore,
    config: PricingConfig,
}

impl Inputs {
    /// 同じ時刻に評価する曲線の入力をまとめる。
    pub(crate) fn new(cases: &[&Case], as_of: Date, extra_discounting: &str) -> Result<Self> {
        Self::with_forecasts(cases, as_of, extra_discounting, |_, _| Ok(()))
    }

    /// OIS の入力に、別のテストで使う参照データと予想曲線の結び付けを加える。
    pub(crate) fn with_forecasts(
        cases: &[&Case],
        as_of: Date,
        extra_discounting: &str,
        add_forecasts: impl FnOnce(&mut RefDataParts, &mut String) -> Result<()>,
    ) -> Result<Self> {
        let clock = LocalEventTime {
            time: "00:00:00".parse()?,
            zone: "Etc/UTC".parse()?,
        };
        let mut parts = RefDataParts::default();
        parts.zones.insert(
            clock.zone.clone(),
            ZoneRules::new(
                "2024-01-01T00:00:00Z".parse()?..="2027-02-15T00:00:00Z".parse()?,
                UtcOffset::from_seconds(0)?,
                Vec::new(),
            )?,
        );
        parts.zones.insert(
            "Asia/Tokyo".parse()?,
            ZoneRules::new(
                "2024-01-01T00:00:00Z".parse()?..="2027-02-15T00:00:00Z".parse()?,
                UtcOffset::from_seconds(9 * 3600)?,
                Vec::new(),
            )?,
        );
        let mut discounting = BTreeMap::new();
        let mut forecasting = String::from("[forecasting]\n");
        for case in cases {
            case.add_refdata(&mut parts, &clock)?;
            discounting.entry(case.currency).or_insert(&case.index);
            writeln!(forecasting, "\"{}\" = \"{}\"", case.index, case.index)
                .map_err(format_error)?;
        }
        // 両通貨の曲線があれば USDJPY の受渡慣行（T+2、16:00 UTC、pip は 0.01 円）を置き、
        // 評価時点の換算レートと FX のテストが参照データから読めるようにする。
        if [Currency::JPY, Currency::USD]
            .iter()
            .all(|ccy| parts.currencies.contains_key(ccy))
        {
            let usd_jpy: CcyPair = "USDJPY".parse()?;
            parts.fx_conventions.insert(
                usd_jpy,
                FxConvention {
                    pair: usd_jpy,
                    spot_lag: FxSpotLag::TwoBusinessDays,
                    settlement: LocalEventTime {
                        time: "16:00:00".parse()?,
                        zone: "Etc/UTC".parse()?,
                    },
                    pip_exponent: 2,
                },
            );
        }
        add_forecasts(&mut parts, &mut forecasting)?;
        let mut config = String::from("version = 1\nreporting_currency = \"JPY\"\n");
        for (ccy, curve) in discounting {
            writeln!(config, "[[discounting]]\ncurrency = \"{ccy}\"\ncollateral = \"uncollateralized\"\ncurve = \"{curve}\"")
                .map_err(format_error)?;
        }
        config.push_str(extra_discounting);
        config.push_str(&forecasting);
        let refdata = RefData::new(parts)?;
        // 各曲線の指標の公表済み fixing を入れる。スポット起点の商品は過去の観測を読まないので、
        // 開始済みのスワップを評価するときだけ値が使われる。
        let published = cases
            .iter()
            .flat_map(|case| {
                case.fixings.iter().map(|row| -> Result<_> {
                    Ok((
                        case.index.parse()?,
                        ObservationPoint::Date(row.date),
                        row.value,
                    ))
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            as_of: midnight(as_of)?,
            fixings: FixingStore::new(
                &refdata,
                FixingStoreParts {
                    published,
                    ..FixingStoreParts::default()
                },
            )?,
            refdata,
            config: toml::from_str(&config)
                .map_err(|error| Error::InvalidData(error.to_string()))?,
        })
    }

    /// 校正と準備に使うコンテキスト。
    pub(crate) const fn context(&self) -> Context<'_> {
        Context {
            as_of: self.as_of,
            refdata: &self.refdata,
            fixings: &self.fixings,
            config: &self.config,
        }
    }
}

fn format_error(error: std::fmt::Error) -> Error {
    Error::Internal(error.to_string())
}

/// 指定した OIS のピラーと校正商品から市場定義を作る。
pub(crate) fn definition(cases: &[(&Case, &CollateralKey)]) -> Result<MarketDefinition> {
    let mut curves = Vec::new();
    let mut bindings = Vec::new();
    for &(case, collateral) in cases {
        curves.push(case.curve()?);
        let block = case.block()?;
        for (index, node) in case.nodes.iter().enumerate() {
            let element =
                u32::try_from(index).map_err(|error| Error::InvalidInput(error.to_string()))?;
            bindings.push(case.binding(
                node,
                collateral,
                ElementSelector::Param(ParamId::new(block.clone(), element)),
            )?);
        }
    }
    MarketDefinition::new(
        cases[0].0.currency.to_string().parse()?,
        "Etc/UTC".parse()?,
        curves,
        bindings,
    )
}

/// ピラーに対応する入力 quote を返す。
pub(crate) fn fixture_quotes(case: &Case) -> Result<Vec<(QuoteId, f64)>> {
    case.nodes
        .iter()
        .map(|node| Ok((case.quote(node)?, node.rate)))
        .collect()
}

/// 指定のテナーと固定金利の OIS を準備する。
pub(crate) fn prepare_ois(
    case: &Case,
    node: &Node,
    rate: f64,
    collateral: &CollateralKey,
    market: &Market<f64>,
    ctx: &Context<'_>,
) -> Result<Prepared> {
    let product = case.template(node)?.build(rate, ctx.refdata)?;
    prepare(
        PrepareInput::Standard {
            product: &product,
            trade_date: case.as_of,
        },
        collateral,
        &[],
        market.structure(),
        ctx,
    )
}
`````

`````
crates/rhodium-pricing/tests/fixtures/quantlib/README.md
# QuantLib のスワップ PV

QuantLib 1.42 の Python バインディングで、固定払い・翌日物複利受けのスワップを評価した。
[QuantLib](https://github.com/lballabio/QuantLib/tree/v1.42) と [QuantLib-SWIG](https://github.com/lballabio/QuantLib-SWIG/tree/v1.42) は Modified BSD License で公開されている。
生成データとスクリプトはこのリポジトリと同じライセンスで配布する。

JPY では Japan カレンダー、ACT/365F の TONA を表す翌日物指標、同じ割引・予想曲線、支払ラグなしで評価する。
USD では QuantLib の SOFR 指標とその fixing カレンダー、ACT/360、独立した割引・予想曲線、2 営業日の支払ラグで評価する。
両ケースとも評価日は 2025 年 1 月 7 日で、fixing は全期間で将来にある。
曲線は評価日を DF 1 の起点にしたログ線形補間とし、入力の DF を直接与える。
日単位の QuantLib と同じ曲線時刻になるよう、rhodium では観測・決済・評価・ピラーを UTC の 00:00 にそろえる。
公表は翌営業日の UTC 00:00 とする。

再生成はリポジトリのルートから次を実行する。

```sh
uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-pricing/tests/fixtures/quantlib/generate.py > crates/rhodium-pricing/tests/fixtures/quantlib/swaps.json
```

Rust のテストは保存した PV と脚別 PV を照合する。
相対誤差は 1e-12、絶対誤差の下限は通貨の最小単位の 1e-6 倍とする。
ネット PV の比較では両脚の絶対値の和も相対誤差の尺度に含める。

## TONA と SOFR の OIS 曲線

`ois_curves.json` は、QuantLib 1.42 の `OISRateHelper` と
`PiecewiseLogLinearDiscount` で 1M・6M・1Y・18M・2Y の quote をブートストラップした結果である。
同じ曲線で割引と翌日物の予想を行う。
約定日は 2025 年 1 月 7 日、スポットラグは 2 営業日、固定脚は年払いとする。
両実装で満期から逆算して短い初回期間を作り、期間境界に Modified Following、支払いに Following を適用する。
JPY は Japan カレンダーと ACT/365F で支払ラグは 0、USD は SOFR の fixing カレンダーと ACT/360 で支払ラグは 2 営業日とする。
ピラーは QuantLib の `LastRelevantDate`（各商品の最終支払日）とする。
評価、観測と支払いを UTC 00:00 にそろえ、曲線の時間軸は両実装で ACT/365F とする。

```sh
uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py > crates/rhodium-pricing/tests/fixtures/quantlib/ois_curves.json
```

Rust では各ピラーの DF と、評価日から最後のピラーまでの 5 つの区間に一つずつ置いた照合点（2025-01-24、2025-04-15、2025-08-12、2026-04-15、2026-10-15）の DF を、相対 1e-12、絶対 1e-14 の下限で照合する。
両実装の停止判定と浮動小数点の演算順が異なるため、この許容誤差を使う。
quote を固定金利にした校正商品を取引と同じ準備で評価し、PV が年金現価 × 2 × 残差の許容誤差（1e-14）以内で 0 になることも確かめる。

同じ曲線で、QuantLib の `OvernightIndexedSwap`（日程は年払い、満期から逆算、期間境界に Modified Following、支払いに Following と各通貨の支払ラグ）を二つ評価し、PV と脚別の PV を `swaps` に置く。

| 名前 | 開始 | テナー | 固定金利 | スプレッド | 固定脚 | 元本 |
|---|---|---|---|---|---|---|
| `spot` | 約定日（評価日）の 2 営業日後 | 18M | 18M の quote + 25bp | 0 | 受け | JPY 1 億、USD 100 万 |
| `seasoned` | 2024-12-05 | 1Y | JPY 0.35%、USD 4.20% | 10bp | 払い | JPY 1 億、USD 100 万 |

開始済みのスワップは、2024-12-05 から評価日の前日までの各営業日の fixing（`fixings`、JPY は 0.23% から、USD は 4.55% から 1bp ずつ 3 日周期で変える値）を確定値として読み、評価日の観測は曲線から予想する。
rhodium の参照データでは fixing の公表を翌営業日の UTC 00:00 とするので、評価日の前日の値は評価の瞬間ちょうどに公表済みになり、評価日の観測期間は評価の瞬間に始まる。
QuantLib はスプレッドを日次の複利の外で足す（rhodium の `Exclusive`）。生成スクリプトは、スプレッドの有無の変動脚の差が Σ 元本 × スプレッド × 計算期間 × DF に一致することを確かめてから出力する。
Rust は QuantLib のノードの DF を θ に入れた市場で、同じ取引を `Swap::from_convention`（スポット起点は開始日を省いて慣行のスポットラグで求める）と同じ fixing で作り、PV、変動脚（固定金利 0 の PV）、固定脚（差）を照合する。
許容誤差は相対 1e-12、絶対の下限は通貨の最小単位の 1e-6 倍とし、相対の尺度には両脚の絶対値の和と、変動脚の二つの項（N·P(s)/P(e) と N）の大きさとして元本の 2 倍を含める。

## 会合日のステップと年末のターンを持つ TONA

`ois_turn.json` は同じ TONA の quote を、年末のターンを含む `PiecewiseLogLinearDiscount` に校正した結果である。
Japan カレンダーで期末の営業日から翌営業日までの日数係数を求め、50bp を掛けた値を ln DF の跳びにする。
QuantLib のジャンプは指定日より後に効くため、QuantLib には発効の前日を指定する。
rhodium は最初のピラーで会合のステップを切り替え、二つ目のピラーで長期の節点につなぐ。
両実装のピラーと中間日の DF を相対 1e-12、絶対 1e-14 で比べ、JPY と USD の OIS quote を再評価する。

```sh
uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py --turn > crates/rhodium-pricing/tests/fixtures/quantlib/ois_turn.json
```

## TONA と TIBOR 6M・3M、SOFR と USD ターム指標の複数曲線

`multicurve.json` は、QuantLib 1.42 の `OISRateHelper` で OIS 曲線を先に校正し、`OvernightIborBasisSwapRateHelper` で別のターム予想曲線を校正した結果である。
OIS とベーシスには 6M、1Y、18M、2Y を使う。
ベーシス商品は 2 営業日後に開始し、半年ごとにラグなしで支払う。
翌日物脚に複利の外で spread を足して払い、ターム脚を受ける。
ターム指標は JPY では QuantLib の TIBOR 6M と TIBOR 3M、USD では SOFR のカレンダーを使うテスト用の 6M 指標で、予想期間の端点は指標の fixing 規則で決める。
両脚の支払いの頻度はターム指標のテナー（TIBOR 3M は四半期ごと）で、出力の `term_tenor` に記録する。
配列の順は TIBOR 6M、USD のターム指標、TIBOR 3M（三つ目は M1 で足した）。
OIS 曲線を両脚の割引と翌日物の予想に使う。

JPY は割引とターム予想のピラーがそろっており、rhodium ではターム予想を TONA とベーシスの和で表す。
共通のピラーの間ではどちらの実装も対数 DF が線形になる。
USD の OIS は 2 営業日の支払ラグがあるため、ターム予想は独立した節点曲線で表す。
Rust のテストはピラーと補間点の割引係数および予想の擬似 DF を相対 1e-12、絶対 1e-14 で比較し、元本 1 の校正商品を quote どおりに評価した PV を残差の許容誤差から求めた幅で照合する。

```sh
uv run --no-project --python 3.13 --with QuantLib==1.42 python crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py --multicurve > crates/rhodium-pricing/tests/fixtures/quantlib/multicurve.json
```
`````

`````
crates/rhodium-pricing/tests/fixtures/quantlib/generate.py
"""QuantLib 1.42 で M0 のスワップ PV と同じ暦の日付を作る。"""

import json

import QuantLib as ql


def date(day, month):
    return ql.Date(day, month, 2025)


def case(currency, name, calendar, day_count, projection, discount, lag, rate):
    as_of = date(7, 1)
    schedule = [date(8, 1), date(10, 2), date(3, 3)]
    pillars = [date(22, 1), date(27, 2), date(14, 3)]

    def curve(discounts):
        return ql.YieldTermStructureHandle(
            ql.DiscountCurve([as_of, *pillars], [1.0, *discounts], ql.Actual365Fixed())
        )

    ql.Settings.instance().evaluationDate = as_of
    if currency == "JPY":
        index = ql.OvernightIndex("TONA", 0, ql.JPYCurrency(), calendar, day_count, curve(projection))
    else:
        index = ql.Sofr(curve(projection))

    notional = 1_000_000.0
    swap = ql.OvernightIndexedSwap(
        ql.OvernightIndexedSwap.Payer, notional, ql.Schedule(schedule, calendar),
        rate, day_count, index, 0.0, lag, ql.Following, calendar,
    )
    swap.setPricingEngine(ql.DiscountingSwapEngine(curve(discount)))
    return {
        "currency": currency,
        "index": name,
        "as_of": as_of.ISO(),
        "dates": [item.ISO() for item in schedule],
        "payments": [flow.date().ISO() for flow in swap.fixedLeg()],
        "pillars": [item.ISO() for item in pillars],
        "holidays": [item.ISO() for item in calendar.holidayList(date(1, 1), date(14, 3))],
        "projection_dfs": projection,
        "discount_dfs": discount,
        "single_curve": currency == "JPY",
        "payment_lag": lag,
        "fixed_rate": rate,
        "notional": notional,
        "pv": swap.NPV(),
        "fixed_pv": swap.fixedLegNPV(),
        "floating_pv": swap.overnightLegNPV(),
    }


assert ql.__version__ == "1.42"
values = [
    case("JPY", "JPY-TONA", ql.Japan(), ql.Actual365Fixed(),
         [0.99991, 0.99966, 0.99934], [0.99991, 0.99966, 0.99934], 0, 0.0125),
    case("USD", "USD-SOFR", ql.Sofr().fixingCalendar(), ql.Actual360(),
         [0.9983, 0.9935, 0.9911], [0.9985, 0.9941, 0.9917], 2, 0.036),
]
print(json.dumps(values, indent=2))
`````

`````
crates/rhodium-pricing/tests/fixtures/quantlib/generate_ois.py
"""QuantLib 1.42 で TONA と SOFR の単一 OIS 曲線をブートストラップし、同じ曲線でスワップを評価する。"""

import json
import math
import sys

import QuantLib as ql

AS_OF = ql.Date(7, 1, 2025)
TENORS = ["1M", "6M", "1Y", "18M", "2Y"]
COVERAGE = (ql.Date(1, 1, 2024), ql.Date(1, 2, 2027))
# 評価日から最後のピラーまでの各区間に一つずつ置く、ピラー以外の照合点。
SAMPLES = [
    ql.Date(24, 1, 2025),
    ql.Date(15, 4, 2025),
    ql.Date(12, 8, 2025),
    ql.Date(15, 4, 2026),
    ql.Date(15, 10, 2026),
]
# 開始済みのスワップの開始日。評価日の前日までの fixing を確定値として渡す。
SEASONED_START = ql.Date(5, 12, 2024)


def fixed_leg_day_count(currency):
    return ql.Actual365Fixed() if currency == "JPY" else ql.Actual360()


def fixings(index, base):
    """開始済みのスワップが読む、評価日の前日までの営業日の fixing を決まった規則で作る。"""
    calendar = index.fixingCalendar()
    values = []
    day = SEASONED_START
    step = 0
    while day < AS_OF:
        if calendar.isBusinessDay(day):
            values.append((day, base + 0.0001 * (step % 3)))
            step += 1
        day += 1
    return values


def swap_case(currency, index, curve, payment_lag, start, tenor, fixed_rate, spread, notional,
              payer):
    """固定払い（payer）か受けの OIS を、慣行どおりの日程で作って評価する。"""
    calendar = index.fixingCalendar()
    schedule = ql.Schedule(
        start, start + ql.Period(tenor), ql.Period(ql.Annual), calendar,
        ql.ModifiedFollowing, ql.ModifiedFollowing, ql.DateGeneration.Backward, False,
    )
    side = ql.Swap.Payer if payer else ql.Swap.Receiver
    # 同じ曲線で予想と割引を行う。fixing は指標の名前で共有されるので複製にも見える。
    handle = ql.YieldTermStructureHandle(curve)
    swap = ql.OvernightIndexedSwap(
        side, notional, schedule, fixed_rate, fixed_leg_day_count(currency),
        index.clone(handle), spread, payment_lag, ql.Following, calendar,
    )
    swap.setPricingEngine(ql.DiscountingSwapEngine(handle))
    return swap


def priced(swap):
    return {
        "pv": swap.NPV(),
        "fixed_pv": swap.fixedLegNPV(),
        "floating_pv": swap.overnightLegNPV(),
    }


def case(currency, name, index, payment_lag, rates, seasoned, turn=None, tenors=TENORS):
    calendar = index.fixingCalendar()
    helpers = [
        ql.OISRateHelper(
            2, ql.Period(tenor), rate, index,
            paymentLag=payment_lag, paymentCalendar=calendar, fixedCalendar=calendar,
            paymentFrequency=ql.Annual, paymentConvention=ql.Following,
            endOfMonth=False, convention=ql.ModifiedFollowing,
        )
        for tenor, rate in zip(tenors, rates, strict=True)
    ]
    bootstrap = ql.IterativeBootstrap(1.0e-14)
    if turn is None:
        curve = ql.PiecewiseLogLinearDiscount(AS_OF, helpers, ql.Actual365Fixed(), bootstrap)
    else:
        effective, log_jump = turn
        # QuantLib は jump 日の翌日から、rhodium は発効日の冒頭から跳びをかける。
        curve = ql.PiecewiseLogLinearDiscount(
            AS_OF, helpers, ql.Actual365Fixed(), [ql.QuoteHandle(ql.SimpleQuote(math.exp(-log_jump)))],
            [effective - 1], ql.LogLinear(), bootstrap,
        )
    nodes = curve.nodes()
    assert nodes[0] == (AS_OF, 1.0)
    assert [helper.pillarDate() for helper in helpers] == [day for day, _ in nodes[1:]]
    assert max(abs(helper.impliedQuote() - rate)
               for helper, rate in zip(helpers, rates, strict=True)) < 2.0e-14

    # 慣行どおりの OIS を二つ評価する。
    # spot: 約定日の 2 営業日後に始まる 18M の off-market の受け。
    # seasoned: 評価日より前に始まった 1Y の固定払いで、確定した fixing とスプレッドを持つ。
    notional = 100_000_000.0 if currency == "JPY" else 1_000_000.0
    spot = calendar.advance(AS_OF, 2, ql.Days)
    seasoned_fixed, seasoned_spread, seasoned_base = seasoned
    history = fixings(index, seasoned_base)
    index.clearFixings()
    for day, value in history:
        index.addFixing(day, value)
    swaps = []
    for label, start, tenor, fixed_rate, spread, payer in [
        ("spot", spot, "18M", rates[tenors.index("18M")] + 0.0025, 0.0, False),
        ("seasoned", SEASONED_START, "1Y", seasoned_fixed, seasoned_spread, True),
    ]:
        swap = swap_case(currency, index, curve, payment_lag, start, tenor, fixed_rate, spread,
                         notional, payer)
        # QuantLib の OIS は spread を複利の外で足す（rhodium の Exclusive）。差で確かめる。
        if spread != 0.0:
            plain = swap_case(currency, index, curve, payment_lag, start, tenor, fixed_rate, 0.0,
                              notional, payer)
            coupons = [ql.as_coupon(cf) for cf in swap.overnightLeg()]
            exclusive = sum(c.accrualPeriod() * spread * notional * curve.discount(c.date())
                            for c in coupons)
            assert abs((swap.overnightLegNPV() - plain.overnightLegNPV()) - exclusive) < 1e-6
        swaps.append({
            "name": label,
            "start": start.ISO(),
            "tenor": tenor,
            "fixed_rate": fixed_rate,
            "spread": spread,
            "notional": notional,
            "fixed_side": "pay" if payer else "receive",
            **priced(swap),
        })
    index.clearFixings()

    return {
        "currency": currency,
        "index": name,
        "as_of": AS_OF.ISO(),
        "calendar_from": COVERAGE[0].ISO(),
        "calendar_through": COVERAGE[1].ISO(),
        "holidays": [day.ISO() for day in calendar.holidayList(*COVERAGE)],
        "payment_lag": payment_lag,
        "nodes": [
            {"tenor": tenor, "rate": rate, "date": day.ISO(),
             "df": df if turn is None else curve.discount(day)}
            for tenor, rate, (day, df) in zip(tenors, rates, nodes[1:], strict=True)
        ],
        "samples": [{"date": day.ISO(), "df": curve.discount(day)} for day in SAMPLES],
        "fixings": [{"date": day.ISO(), "value": value} for day, value in history],
        "swaps": swaps,
    }, curve


def multicurve_case(currency, overnight_id, index, term_id, term, lag, rates, seasoned, spreads,
                    term_tenor="6M"):
    """OIS で割引と翌日物の予想を先に解き、ベーシスでターム指標の予想を解く。"""
    tenors = TENORS[1:]
    overnight, discount = case(currency, overnight_id, index, lag, rates, seasoned, tenors=tenors)
    calendar = index.fixingCalendar()
    handle = ql.YieldTermStructureHandle(discount)
    helpers = [
        ql.OvernightIborBasisSwapRateHelper(
            ql.QuoteHandle(ql.SimpleQuote(spread)), ql.Period(tenor), 2,
            calendar, ql.ModifiedFollowing, False, index.clone(handle), term, handle,
        )
        for tenor, spread in zip(tenors, spreads, strict=True)
    ]
    projection = ql.PiecewiseLogLinearDiscount(
        AS_OF, helpers, ql.Actual365Fixed(), ql.IterativeBootstrap(1.0e-14),
    )
    nodes = projection.nodes()[1:]
    assert max(abs(helper.impliedQuote() - spread)
               for helper, spread in zip(helpers, spreads, strict=True)) < 2.0e-14
    if currency == "JPY":
        # 共通の節点なら、対数 DF の差も同じ節点の間で線形になる。
        assert [day.ISO() for day, _ in nodes] == [row["date"] for row in overnight["nodes"]]
    return {
        "overnight_nodes": overnight["nodes"],
        "overnight_samples": overnight["samples"],
        "term_index": term_id,
        "term_tenor": term_tenor,
        "nodes": [
            {"tenor": tenor, "rate": spread, "date": day.ISO(), "df": df}
            for tenor, spread, (day, df) in zip(tenors, spreads, nodes, strict=True)
        ],
        "samples": [{"date": day.ISO(), "df": projection.discount(day)} for day in SAMPLES],
    }


assert ql.__version__ == "1.42"
ql.Settings.instance().evaluationDate = AS_OF
tona = ql.OvernightIndex("TONA", 0, ql.JPYCurrency(), ql.Japan(), ql.Actual365Fixed())
tona_rates = [-0.001, 0.002, 0.004, 0.005, 0.007]
tona_seasoned = (0.0035, 0.0010, 0.0023)
if "--turn" in sys.argv:
    calendar = tona.fixingCalendar()
    effective = calendar.adjust(ql.Date(1, 1, 2026), ql.Following)
    start = calendar.advance(effective, -1, ql.Days)
    fraction = ql.Actual365Fixed().yearFraction(start, effective)
    spread = 0.005
    data, _ = case("JPY", "JPY-TONA", tona, 0, tona_rates, tona_seasoned, (effective, spread * fraction))
    print(json.dumps({
        "effective": effective.ISO(),
        "spread": spread,
        "fraction": fraction,
        "nodes": data["nodes"],
        "samples": data["samples"],
    }, indent=2))
elif "--multicurve" in sys.argv:
    # fixing の開始と終了は指標のテナーで決める。両実装で支払期間と区別する。
    ql.IborCoupon.createIndexedCoupons()
    sofr = ql.Sofr()
    usd_term = ql.IborIndex(
        "TestUSDTerm", ql.Period("6M"), 2, ql.USDCurrency(), sofr.fixingCalendar(),
        ql.ModifiedFollowing, False, ql.Actual360(),
    )
    print(json.dumps([
        multicurve_case("JPY", "JPY-TONA", tona, "JPY-TIBOR-6M", ql.Tibor(ql.Period("6M")), 0,
                        tona_rates[1:], tona_seasoned, [0.0012, 0.0022, 0.0017, 0.0028]),
        multicurve_case("USD", "USD-SOFR", sofr, "USD-TERM-6M", usd_term, 2,
                        [0.042, 0.041, 0.040, 0.039], (0.0420, 0.0010, 0.0455),
                        [-0.0007, -0.0002, 0.0004, 0.0008]),
        # TIBOR 3M も TONA へのベーシスで、両脚を四半期ごとに払う。
        multicurve_case("JPY", "JPY-TONA", tona, "JPY-TIBOR-3M", ql.Tibor(ql.Period("3M")), 0,
                        tona_rates[1:], tona_seasoned, [0.0008, 0.0015, 0.0011, 0.0019],
                        term_tenor="3M"),
    ], indent=2))
else:
    print(json.dumps([
        case("JPY", "JPY-TONA", tona, 0, tona_rates, tona_seasoned)[0],
        case("USD", "USD-SOFR", ql.Sofr(), 2, [0.043, 0.042, 0.041, 0.040, 0.039],
             (0.0420, 0.0010, 0.0455))[0],
    ], indent=2))
`````

`````
crates/rhodium-pricing/tests/fixtures/quantlib/multicurve.json
[
  {
    "overnight_nodes": [
      {
        "tenor": "6M",
        "rate": 0.002,
        "date": "2025-07-09",
        "df": 0.9989982592720291
      },
      {
        "tenor": "1Y",
        "rate": 0.004,
        "date": "2026-01-09",
        "df": 0.9960050264807817
      },
      {
        "tenor": "18M",
        "rate": 0.005,
        "date": "2026-07-09",
        "df": 0.992549331641807
      },
      {
        "tenor": "2Y",
        "rate": 0.007,
        "date": "2027-01-12",
        "df": 0.9860578741851718
      }
    ],
    "overnight_samples": [
      {
        "date": "2025-01-24",
        "df": 0.9999068998111723
      },
      {
        "date": "2025-04-15",
        "df": 0.9994634238175699
      },
      {
        "date": "2025-08-12",
        "df": 0.9984444851883086
      },
      {
        "date": "2026-04-15",
        "df": 0.994170676035945
      },
      {
        "date": "2026-10-15",
        "df": 0.9891420796400195
      }
    ],
    "term_index": "JPY-TIBOR-6M",
    "term_tenor": "6M",
    "nodes": [
      {
        "tenor": "6M",
        "rate": 0.0012,
        "date": "2025-07-09",
        "df": 0.9983981723483519
      },
      {
        "tenor": "1Y",
        "rate": 0.0022,
        "date": "2026-01-09",
        "df": 0.9938150348912135
      },
      {
        "tenor": "18M",
        "rate": 0.0017,
        "date": "2026-07-09",
        "df": 0.9900309775934285
      },
      {
        "tenor": "2Y",
        "rate": 0.0028,
        "date": "2027-01-12",
        "df": 0.9805397593759969
      }
    ],
    "samples": [
      {
        "date": "2025-01-24",
        "df": 0.9998510881205138
      },
      {
        "date": "2025-04-15",
        "df": 0.9991418712937897
      },
      {
        "date": "2025-08-12",
        "df": 0.9975496992122486
      },
      {
        "date": "2026-04-15",
        "df": 0.9918062233474458
      },
      {
        "date": "2026-10-15",
        "df": 0.985045568902551
      }
    ]
  },
  {
    "overnight_nodes": [
      {
        "tenor": "6M",
        "rate": 0.042,
        "date": "2025-07-11",
        "df": 0.9788678723178676
      },
      {
        "tenor": "1Y",
        "rate": 0.041,
        "date": "2026-01-13",
        "df": 0.9594545329223071
      },
      {
        "tenor": "18M",
        "rate": 0.04,
        "date": "2026-07-13",
        "df": 0.9414867478891797
      },
      {
        "tenor": "2Y",
        "rate": 0.039,
        "date": "2027-01-13",
        "df": 0.9248568221589608
      }
    ],
    "overnight_samples": [
      {
        "date": "2025-01-24",
        "df": 0.9980392419565411
      },
      {
        "date": "2025-04-15",
        "df": 0.9887494765093302
      },
      {
        "date": "2025-08-12",
        "df": 0.9755001853706653
      },
      {
        "date": "2026-04-15",
        "df": 0.950279290958284
      },
      {
        "date": "2026-10-15",
        "df": 0.9329539991729342
      }
    ],
    "term_index": "USD-TERM-6M",
    "term_tenor": "6M",
    "nodes": [
      {
        "tenor": "6M",
        "rate": -0.0007,
        "date": "2025-07-09",
        "df": 0.9794352305461219
      },
      {
        "tenor": "1Y",
        "rate": -0.0002,
        "date": "2026-01-09",
        "df": 0.96005741554781
      },
      {
        "tenor": "18M",
        "rate": 0.0004,
        "date": "2026-07-09",
        "df": 0.9413020142845046
      },
      {
        "tenor": "2Y",
        "rate": 0.0008,
        "date": "2027-01-11",
        "df": 0.9235220706683834
      }
    ],
    "samples": [
      {
        "date": "2025-01-24",
        "df": 0.9980715565426879
      },
      {
        "date": "2025-04-15",
        "df": 0.9889340408949315
      },
      {
        "date": "2025-08-12",
        "df": 0.9758253206796759
      },
      {
        "date": "2026-04-15",
        "df": 0.9500637265357552
      },
      {
        "date": "2026-10-15",
        "df": 0.9318918355727929
      }
    ]
  },
  {
    "overnight_nodes": [
      {
        "tenor": "6M",
        "rate": 0.002,
        "date": "2025-07-09",
        "df": 0.9989982592720291
      },
      {
        "tenor": "1Y",
        "rate": 0.004,
        "date": "2026-01-09",
        "df": 0.9960050264807817
      },
      {
        "tenor": "18M",
        "rate": 0.005,
        "date": "2026-07-09",
        "df": 0.992549331641807
      },
      {
        "tenor": "2Y",
        "rate": 0.007,
        "date": "2027-01-12",
        "df": 0.9860578741851718
      }
    ],
    "overnight_samples": [
      {
        "date": "2025-01-24",
        "df": 0.9999068998111723
      },
      {
        "date": "2025-04-15",
        "df": 0.9994634238175699
      },
      {
        "date": "2025-08-12",
        "df": 0.9984444851883086
      },
      {
        "date": "2026-04-15",
        "df": 0.994170676035945
      },
      {
        "date": "2026-10-15",
        "df": 0.9891420796400195
      }
    ],
    "term_index": "JPY-TIBOR-3M",
    "term_tenor": "3M",
    "nodes": [
      {
        "tenor": "6M",
        "rate": 0.0008,
        "date": "2025-07-09",
        "df": 0.9985978836829034
      },
      {
        "tenor": "1Y",
        "rate": 0.0015,
        "date": "2026-01-09",
        "df": 0.9945091015080203
      },
      {
        "tenor": "18M",
        "rate": 0.0011,
        "date": "2026-07-09",
        "df": 0.9909165069568865
      },
      {
        "tenor": "2Y",
        "rate": 0.0019,
        "date": "2027-01-12",
        "df": 0.9822993956366051
      }
    ],
    "samples": [
      {
        "date": "2025-01-24",
        "df": 0.9998696658525201
      },
      {
        "date": "2025-04-15",
        "df": 0.9992488952271131
      },
      {
        "date": "2025-08-12",
        "df": 0.997841083751653
      },
      {
        "date": "2026-04-15",
        "df": 0.9926020178702298
      },
      {
        "date": "2026-10-15",
        "df": 0.9863912022439585
      }
    ]
  }
]
`````

`````
crates/rhodium-pricing/tests/fixtures/quantlib/ois_curves.json
[
  {
    "currency": "JPY",
    "index": "JPY-TONA",
    "as_of": "2025-01-07",
    "calendar_from": "2024-01-01",
    "calendar_through": "2027-02-01",
    "holidays": [
      "2024-01-01",
      "2024-01-02",
      "2024-01-03",
      "2024-01-08",
      "2024-02-12",
      "2024-02-23",
      "2024-03-20",
      "2024-04-29",
      "2024-05-03",
      "2024-05-06",
      "2024-07-15",
      "2024-08-12",
      "2024-09-16",
      "2024-09-23",
      "2024-10-14",
      "2024-11-04",
      "2024-12-31",
      "2025-01-01",
      "2025-01-02",
      "2025-01-03",
      "2025-01-13",
      "2025-02-11",
      "2025-02-24",
      "2025-03-20",
      "2025-04-29",
      "2025-05-05",
      "2025-05-06",
      "2025-07-21",
      "2025-08-11",
      "2025-09-15",
      "2025-09-23",
      "2025-10-13",
      "2025-11-03",
      "2025-11-24",
      "2025-12-31",
      "2026-01-01",
      "2026-01-02",
      "2026-01-12",
      "2026-02-11",
      "2026-02-23",
      "2026-03-20",
      "2026-04-29",
      "2026-05-04",
      "2026-05-05",
      "2026-05-06",
      "2026-07-20",
      "2026-08-11",
      "2026-09-21",
      "2026-09-22",
      "2026-09-23",
      "2026-10-12",
      "2026-11-03",
      "2026-11-23",
      "2026-12-31",
      "2027-01-01",
      "2027-01-11"
    ],
    "payment_lag": 0,
    "nodes": [
      {
        "tenor": "1M",
        "rate": -0.001,
        "date": "2025-02-10",
        "df": 1.0000931591075282
      },
      {
        "tenor": "6M",
        "rate": 0.002,
        "date": "2025-07-09",
        "df": 0.9990146761106965
      },
      {
        "tenor": "1Y",
        "rate": 0.004,
        "date": "2026-01-09",
        "df": 0.9960213941307547
      },
      {
        "tenor": "18M",
        "rate": 0.005,
        "date": "2026-07-09",
        "df": 0.9925656425033084
      },
      {
        "tenor": "2Y",
        "rate": 0.007,
        "date": "2027-01-12",
        "df": 0.9860740783706016
      }
    ],
    "samples": [
      {
        "date": "2025-01-24",
        "df": 1.0000465784689871
      },
      {
        "date": "2025-04-15",
        "df": 0.9996297755188744
      },
      {
        "date": "2025-08-12",
        "df": 0.9984608929266401
      },
      {
        "date": "2026-04-15",
        "df": 0.994187013541486
      },
      {
        "date": "2026-10-15",
        "df": 0.9891583345091246
      }
    ],
    "fixings": [
      {
        "date": "2024-12-05",
        "value": 0.0023
      },
      {
        "date": "2024-12-06",
        "value": 0.0024
      },
      {
        "date": "2024-12-09",
        "value": 0.0025
      },
      {
        "date": "2024-12-10",
        "value": 0.0023
      },
      {
        "date": "2024-12-11",
        "value": 0.0024
      },
      {
        "date": "2024-12-12",
        "value": 0.0025
      },
      {
        "date": "2024-12-13",
        "value": 0.0023
      },
      {
        "date": "2024-12-16",
        "value": 0.0024
      },
      {
        "date": "2024-12-17",
        "value": 0.0025
      },
      {
        "date": "2024-12-18",
        "value": 0.0023
      },
      {
        "date": "2024-12-19",
        "value": 0.0024
      },
      {
        "date": "2024-12-20",
        "value": 0.0025
      },
      {
        "date": "2024-12-23",
        "value": 0.0023
      },
      {
        "date": "2024-12-24",
        "value": 0.0024
      },
      {
        "date": "2024-12-25",
        "value": 0.0025
      },
      {
        "date": "2024-12-26",
        "value": 0.0023
      },
      {
        "date": "2024-12-27",
        "value": 0.0024
      },
      {
        "date": "2024-12-30",
        "value": 0.0025
      },
      {
        "date": "2025-01-06",
        "value": 0.0023
      }
    ],
    "swaps": [
      {
        "name": "spot",
        "start": "2025-01-09",
        "tenor": "18M",
        "fixed_rate": 0.0075,
        "spread": 0.0,
        "notional": 100000000.0,
        "fixed_side": "receive",
        "pv": 371991.86019845516,
        "fixed_pv": 1115975.5805953604,
        "floating_pv": -743983.7203969052
      },
      {
        "name": "seasoned",
        "start": "2024-12-05",
        "tenor": "1Y",
        "fixed_rate": 0.0035,
        "spread": 0.001,
        "notional": 100000000.0,
        "fixed_side": "pay",
        "pv": 113682.5723520185,
        "fixed_pv": -348806.5269017113,
        "floating_pv": 462489.0992537298
      }
    ]
  },
  {
    "currency": "USD",
    "index": "USD-SOFR",
    "as_of": "2025-01-07",
    "calendar_from": "2024-01-01",
    "calendar_through": "2027-02-01",
    "holidays": [
      "2024-01-01",
      "2024-01-15",
      "2024-02-19",
      "2024-03-29",
      "2024-05-27",
      "2024-06-19",
      "2024-07-04",
      "2024-09-02",
      "2024-10-14",
      "2024-11-11",
      "2024-11-28",
      "2024-12-25",
      "2025-01-01",
      "2025-01-20",
      "2025-02-17",
      "2025-04-18",
      "2025-05-26",
      "2025-06-19",
      "2025-07-04",
      "2025-09-01",
      "2025-10-13",
      "2025-11-11",
      "2025-11-27",
      "2025-12-25",
      "2026-01-01",
      "2026-01-19",
      "2026-02-16",
      "2026-04-03",
      "2026-05-25",
      "2026-06-19",
      "2026-07-03",
      "2026-09-07",
      "2026-10-12",
      "2026-11-11",
      "2026-11-26",
      "2026-12-25",
      "2027-01-01",
      "2027-01-18"
    ],
    "payment_lag": 2,
    "nodes": [
      {
        "tenor": "1M",
        "rate": 0.043,
        "date": "2025-02-12",
        "df": 0.9957173935230357
      },
      {
        "tenor": "6M",
        "rate": 0.042,
        "date": "2025-07-11",
        "df": 0.9788622065572212
      },
      {
        "tenor": "1Y",
        "rate": 0.041,
        "date": "2026-01-13",
        "df": 0.9594472718772433
      },
      {
        "tenor": "18M",
        "rate": 0.04,
        "date": "2026-07-13",
        "df": 0.941479661356336
      },
      {
        "tenor": "2Y",
        "rate": 0.039,
        "date": "2027-01-13",
        "df": 0.9248498582455866
      }
    ],
    "samples": [
      {
        "date": "2025-01-24",
        "df": 0.9979753675405209
      },
      {
        "date": "2025-04-15",
        "df": 0.9886688517214712
      },
      {
        "date": "2025-08-12",
        "df": 0.9754942403996766
      },
      {
        "date": "2026-04-15",
        "df": 0.950272119119701
      },
      {
        "date": "2026-10-15",
        "df": 0.932946975549917
      }
    ],
    "fixings": [
      {
        "date": "2024-12-05",
        "value": 0.0455
      },
      {
        "date": "2024-12-06",
        "value": 0.0456
      },
      {
        "date": "2024-12-09",
        "value": 0.0457
      },
      {
        "date": "2024-12-10",
        "value": 0.0455
      },
      {
        "date": "2024-12-11",
        "value": 0.0456
      },
      {
        "date": "2024-12-12",
        "value": 0.0457
      },
      {
        "date": "2024-12-13",
        "value": 0.0455
      },
      {
        "date": "2024-12-16",
        "value": 0.0456
      },
      {
        "date": "2024-12-17",
        "value": 0.0457
      },
      {
        "date": "2024-12-18",
        "value": 0.0455
      },
      {
        "date": "2024-12-19",
        "value": 0.0456
      },
      {
        "date": "2024-12-20",
        "value": 0.0457
      },
      {
        "date": "2024-12-23",
        "value": 0.0455
      },
      {
        "date": "2024-12-24",
        "value": 0.0456
      },
      {
        "date": "2024-12-26",
        "value": 0.0457
      },
      {
        "date": "2024-12-27",
        "value": 0.0455
      },
      {
        "date": "2024-12-30",
        "value": 0.0456
      },
      {
        "date": "2024-12-31",
        "value": 0.0457
      },
      {
        "date": "2025-01-02",
        "value": 0.0455
      },
      {
        "date": "2025-01-03",
        "value": 0.0456
      },
      {
        "date": "2025-01-06",
        "value": 0.0457
      }
    ],
    "swaps": [
      {
        "name": "spot",
        "start": "2025-01-09",
        "tenor": "18M",
        "fixed_rate": 0.0425,
        "spread": 0.0,
        "notional": 1000000.0,
        "fixed_side": "receive",
        "pv": 3616.7648318189604,
        "fixed_pv": 61485.00214092103,
        "floating_pv": -57868.23730910207
      },
      {
        "name": "seasoned",
        "start": "2024-12-05",
        "tenor": "1Y",
        "fixed_rate": 0.042,
        "spread": 0.001,
        "notional": 1000000.0,
        "fixed_side": "pay",
        "pv": 650.0919172484646,
        "fixed_pv": -41010.77227553883,
        "floating_pv": 41660.864192787296
      }
    ]
  }
]
`````

`````
crates/rhodium-pricing/tests/fixtures/quantlib/ois_turn.json
{
  "effective": "2026-01-05",
  "spread": 0.005,
  "fraction": 0.01643835616438356,
  "nodes": [
    {
      "tenor": "1M",
      "rate": -0.001,
      "date": "2025-02-10",
      "df": 1.0000931591075282
    },
    {
      "tenor": "6M",
      "rate": 0.002,
      "date": "2025-07-09",
      "df": 0.9990146761106965
    },
    {
      "tenor": "1Y",
      "rate": 0.004,
      "date": "2026-01-09",
      "df": 0.9960213941307549
    },
    {
      "tenor": "18M",
      "rate": 0.005,
      "date": "2026-07-09",
      "df": 0.9925656425033085
    },
    {
      "tenor": "2Y",
      "rate": 0.007,
      "date": "2027-01-12",
      "df": 0.9860740783706015
    }
  ],
  "samples": [
    {
      "date": "2025-01-24",
      "df": 1.0000465784689871
    },
    {
      "date": "2025-04-15",
      "df": 0.9996297755188744
    },
    {
      "date": "2025-08-12",
      "df": 0.998476057278108
    },
    {
      "date": "2026-04-15",
      "df": 0.9941870135414861
    },
    {
      "date": "2026-10-15",
      "df": 0.9891583345091246
    }
  ]
}
`````

`````
crates/rhodium-pricing/tests/fixtures/quantlib/swaps.json
[
  {
    "currency": "JPY",
    "index": "JPY-TONA",
    "as_of": "2025-01-07",
    "dates": [
      "2025-01-08",
      "2025-02-10",
      "2025-03-03"
    ],
    "payments": [
      "2025-02-10",
      "2025-03-03"
    ],
    "pillars": [
      "2025-01-22",
      "2025-02-27",
      "2025-03-14"
    ],
    "holidays": [
      "2025-01-01",
      "2025-01-02",
      "2025-01-03",
      "2025-01-13",
      "2025-02-11",
      "2025-02-24"
    ],
    "projection_dfs": [
      0.99991,
      0.99966,
      0.99934
    ],
    "discount_dfs": [
      0.99991,
      0.99966,
      0.99934
    ],
    "single_curve": true,
    "payment_lag": 0,
    "fixed_rate": 0.0125,
    "notional": 1000000.0,
    "pv": -1429.4152354267912,
    "fixed_pv": -1848.7583344482036,
    "floating_pv": 419.3430990214125
  },
  {
    "currency": "USD",
    "index": "USD-SOFR",
    "as_of": "2025-01-07",
    "dates": [
      "2025-01-08",
      "2025-02-10",
      "2025-03-03"
    ],
    "payments": [
      "2025-02-12",
      "2025-03-05"
    ],
    "pillars": [
      "2025-01-22",
      "2025-02-27",
      "2025-03-14"
    ],
    "holidays": [
      "2025-01-01",
      "2025-01-20",
      "2025-02-17"
    ],
    "projection_dfs": [
      0.9983,
      0.9935,
      0.9911
    ],
    "discount_dfs": [
      0.9985,
      0.9941,
      0.9917
    ],
    "single_curve": false,
    "payment_lag": 2,
    "fixed_rate": 0.036,
    "notional": 1000000.0,
    "pv": 1656.4904942248068,
    "fixed_pv": -5372.164745927916,
    "floating_pv": 7028.655240152722
  }
]
`````

`````
crates/rhodium-pricing/tests/fixtures/steps.rs
//! 定型テナーの OIS だけで会合日の階段を当てる、罰則付きの最小二乗の TONA 曲線（R4）。
//!
//! 曲線は Sum[Stitch{Steps(1 年までの 4 回の会合), Nodes(18M, 2Y), 1Y で接続}, Jumps(年末)]。
//! 1M、6M、1Y の OIS が 5 区間の階段を、18M と 2Y の OIS が長期の節点を当てる。
//! 階段は quote より未知数が多いので、隣り合う区間のフォワードの差の罰則で一意にする。
//! 罰則を強くすると OIS の quote を再現しきれず、残差の残る過剰決定の当てはめになる。

use rhodium_core::{Result, date::Date, quantity::YearFraction};
use rhodium_market::{
    curve::{CurveExpr, CurveId, CurveRole, CurveSpec, Extrapolation},
    quote::{QuoteKind, QuoteSet},
    structure::{BlockId, ParamId},
};
use rhodium_pricing::{
    definition::{
        CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
        StepPenalty,
    },
    resolve::CollateralKey,
};
use serde::Deserialize;

use super::ois::{Case, Node, Sample, cases, fixture_quotes, midnight};

/// `QuantLib` で年末のターンを含めて校正した TONA の quote と DF。
#[derive(Deserialize)]
pub struct QuantLibTurn {
    /// 年末の跳びが発効する日。
    pub effective: Date,
    /// 年末のスプレッド。
    pub spread: f64,
    /// 年末の受渡期間の日数係数。
    pub fraction: f64,
    /// OIS の quote とピラーの DF。
    pub nodes: Vec<Node>,
    /// ピラー以外の照合点の DF。
    pub samples: Vec<Sample>,
}

/// 年末のターンの quote の ID。
pub const TURN_QUOTE: &str = "JPY.TONA.TURN";
/// OIS の当てはめの許容幅 σ（0.5bp）。
pub const OIS_TOLERANCE: f64 = 0.5e-4;
/// 罰則の重み λ。
pub const PENALTY_WEIGHT: f64 = 1.0;
/// 罰則の参照の幅 s（1bp）。隣の区間との 1bp の差が、0.5bp の残差の 1/4 の重さになる。
pub const PENALTY_SCALE: f64 = 1e-4;
/// 1 年までの会合の発効日。
const MEETINGS: [&str; 4] = ["2025-03-19", "2025-06-17", "2025-09-19", "2025-12-19"];

/// 年末のターンの `QuantLib` の fixture を読む。
pub(crate) fn quantlib_turn() -> Result<QuantLibTurn> {
    serde_json::from_str(include_str!("quantlib/ois_turn.json"))
        .map_err(|error| rhodium_core::Error::InvalidData(error.to_string()))
}

/// 罰則付きの TONA の入力。`case.nodes` は `QuantLib` のターンの fixture の 5 テナー。
pub struct PenalizedTona {
    /// TONA の OIS の慣行と quote。
    pub case: Case,
    /// 年末のターン。
    pub turn: QuantLibTurn,
}

impl PenalizedTona {
    /// TONA の慣行と、ターンの fixture の quote を読む。
    pub(crate) fn new() -> Result<Self> {
        let [mut case, _] = cases()?;
        let turn = quantlib_turn()?;
        case.nodes.clone_from(&turn.nodes);
        Ok(Self { case, turn })
    }

    /// 曲線の ID。
    pub(crate) fn id(&self) -> Result<CurveId> {
        self.case.index.parse()
    }

    /// 会合日の階段のブロック。
    pub(crate) fn meetings(&self) -> Result<BlockId> {
        BlockId::new(self.id()?, "meetings")
    }

    /// 長期の節点のブロック。
    pub(crate) fn long(&self) -> Result<BlockId> {
        BlockId::new(self.id()?, "long")
    }

    /// 年末の跳びのブロック。
    pub(crate) fn turn_block(&self) -> Result<BlockId> {
        BlockId::new(self.id()?, "turn")
    }

    /// 階段、長期の節点、年末の跳びからなる TONA 曲線。
    pub(crate) fn curve(&self) -> Result<CurveSpec> {
        let nodes = &self.turn.nodes;
        CurveSpec::from_expr(
            self.id()?,
            self.case.currency,
            [CurveRole::Discount, CurveRole::Projection],
            CurveExpr::Sum(vec![
                CurveExpr::Stitch {
                    short: Box::new(CurveExpr::Steps {
                        part: "meetings".into(),
                        meetings: MEETINGS
                            .iter()
                            .map(|date| midnight(date.parse()?))
                            .collect::<Result<_>>()?,
                    }),
                    long: Box::new(CurveExpr::Nodes {
                        part: "long".into(),
                        pillars: nodes[3..]
                            .iter()
                            .map(|node| midnight(node.date))
                            .collect::<Result<_>>()?,
                        extrapolation: Extrapolation::Reject,
                    }),
                    at: midnight(nodes[2].date)?,
                },
                CurveExpr::Jumps {
                    part: "turn".into(),
                    dates: vec![midnight(self.turn.effective)?],
                },
            ]),
        )
    }

    /// 罰則付きの市場の定義。`gauss_newton` なら最小二乗の記録を Gauss–Newton で作る。
    pub(crate) fn definition(&self, gauss_newton: bool) -> Result<MarketDefinition> {
        let collateral = CollateralKey::Uncollateralized;
        let mut bindings = self
            .turn
            .nodes
            .iter()
            .enumerate()
            .map(|(position, node)| {
                let block = if position < 3 {
                    self.meetings()?
                } else {
                    self.long()?
                };
                self.case
                    .binding(node, &collateral, ElementSelector::Block(block))?
                    .with_tolerance(OIS_TOLERANCE)
            })
            .collect::<Result<Vec<_>>>()?;
        // 年末のスプレッドは直接マークで、変換 TurnSpreadToJump が s·τ_turn の跳びにする。
        bindings.push(InstrumentBinding::new(
            TURN_QUOTE.parse()?,
            QuoteKind::Spread,
            CalibrationTarget::DirectMark,
            QuoteTransform::TurnSpreadToJump(YearFraction::new(
                self.turn.fraction,
                self.case.basis(),
            )?),
            collateral,
            ElementSelector::Param(ParamId::new(self.turn_block()?, 0)),
        )?);
        let definition = MarketDefinition::new(
            self.case.currency.to_string().parse()?,
            "Etc/UTC".parse()?,
            [self.curve()?],
            bindings,
        )?
        .with_step_penalties([StepPenalty::new(
            self.meetings()?,
            PENALTY_WEIGHT,
            PENALTY_SCALE,
        )?])?;
        if gauss_newton {
            definition.with_gauss_newton([self.meetings()?, self.long()?])
        } else {
            Ok(definition)
        }
    }

    /// 会合 1 回（1M のピラー）で区切った 2 区間の階段を 1M と 6M の OIS が、1Y で接続した長期の
    /// 節点を 18M と 2Y の OIS が、年末の跳びを直接マークが決める正方系の市場の定義（R3）。
    pub(crate) fn square_definition(&self) -> Result<MarketDefinition> {
        let nodes = &self.turn.nodes;
        let curve = CurveSpec::from_expr(
            self.id()?,
            self.case.currency,
            [CurveRole::Discount, CurveRole::Projection],
            CurveExpr::Sum(vec![
                CurveExpr::Stitch {
                    short: Box::new(CurveExpr::Steps {
                        part: "meetings".into(),
                        meetings: vec![midnight(nodes[0].date)?],
                    }),
                    long: Box::new(CurveExpr::Nodes {
                        part: "long".into(),
                        pillars: nodes[2..]
                            .iter()
                            .map(|node| midnight(node.date))
                            .collect::<Result<_>>()?,
                        extrapolation: Extrapolation::Reject,
                    }),
                    at: midnight(nodes[1].date)?,
                },
                CurveExpr::Jumps {
                    part: "turn".into(),
                    dates: vec![midnight(self.turn.effective)?],
                },
            ]),
        )?;
        let collateral = CollateralKey::Uncollateralized;
        let mut bindings = nodes
            .iter()
            .enumerate()
            .map(|(position, node)| {
                let block = if position < 2 {
                    self.meetings()?
                } else {
                    self.long()?
                };
                self.case
                    .binding(node, &collateral, ElementSelector::Block(block))
            })
            .collect::<Result<Vec<_>>>()?;
        bindings.push(InstrumentBinding::new(
            TURN_QUOTE.parse()?,
            QuoteKind::Spread,
            CalibrationTarget::DirectMark,
            QuoteTransform::TurnSpreadToJump(YearFraction::new(
                self.turn.fraction,
                self.case.basis(),
            )?),
            collateral,
            ElementSelector::Param(ParamId::new(self.turn_block()?, 0)),
        )?);
        MarketDefinition::new(
            self.case.currency.to_string().parse()?,
            "Etc/UTC".parse()?,
            [curve],
            bindings,
        )
    }

    /// OIS とターンの quote。
    pub(crate) fn quotes(&self) -> Result<QuoteSet> {
        QuoteSet::new(
            fixture_quotes(&self.case)?
                .into_iter()
                .chain([(TURN_QUOTE.parse()?, self.turn.spread)]),
        )
    }
}
`````

`````
crates/rhodium-pricing/tests/fixtures/xccy.rs
//! TONA と SOFR の OIS、USDJPY の FX ポイントと `MtM` XCCY で校正する、JPY 現金担保の USD 曲線（R6）。
//!
//! 外貨担保曲線 `USD-under-JPY` は `FxImplied` で、短期は FX スワップのポイント（ON、TN と
//! 3M・6M の日付付きのポイント）、長期は TONA 対 SOFR の `MtM` XCCY のベーシス（1Y、18M）が決める。
//! XCCY は JPY の TONA 脚にベーシスを付けて払い、USD の SOFR 脚の元本を USDJPY の fixing で
//! 期首ごとにリセットする。JPY 現金担保の下で JPY は JPY-TONA、USD は外貨担保曲線で割り引く。

use std::fmt::Write as _;

use rhodium_core::{
    Error, Result,
    calendar::CalendarSpec,
    currency::Currency,
    fx::CcyPair,
    id::IndexId,
    quantity::Spread,
    refdata::{
        Convention, DailyPublication, IndexDefinition, IndexKind, IndexStatus, PublicationRule,
        SwapConvention, SwapCouponConvention, SwapLegConvention,
    },
    tenor::TenorLength,
    time::{Instant, LocalEventTime},
};
use rhodium_market::{
    curve::{CurveId, Extrapolation},
    fx::FxCollateral,
    quote::QuoteSet,
};
use rhodium_pricing::{
    definition::{
        DirectFx, FxDatedQuote, FxDefinition, FxImpliedCurve, FxImpliedLongEnd, FxQuotes,
        MarketDefinition,
    },
    resolve::CollateralKey,
};
use rhodium_products::{
    agreement::{CashCollateralTerms, Remuneration},
    standard::InstrumentTemplate,
};

use super::ois::{Case, Inputs, cases, definition, fixture_quotes, midnight};

/// JPY 現金担保の USD 割引曲線（外貨担保曲線）。
pub const IMPLIED: &str = "USD-under-JPY";
/// `MtM` XCCY の慣行。
pub const XCCY: &str = "USDJPY-XCCY";
/// 元本のリセットに使う USDJPY の fixing。
pub const FX_FIXING: &str = "USDJPY-WMR";
/// スポットのアウトライト。
pub const SPOT: f64 = 150.0;
/// 短期と日付付きの FX ポイント（pip、0.01 円）と、その受渡日。
pub const POINTS: [(&str, Option<&str>, f64); 4] = [
    ("USDJPY.ON", None, -2.0),
    ("USDJPY.TN", None, -2.5),
    ("USDJPY.3M", Some("2025-04-09"), -95.0),
    ("USDJPY.6M", Some("2025-07-09"), -190.0),
];
/// 長期の XCCY のテナー、quote（TONA 脚のベーシス）とピラー（最終の支払日）。
pub const XCCY_QUOTES: [(&str, &str, f64, &str); 2] = [
    ("USDJPY.XCCY.1Y", "1Y", -0.0035, "2026-01-13"),
    ("USDJPY.XCCY.18M", "18M", -0.0040, "2026-07-13"),
];

/// TONA と SOFR の慣行と quote。
pub struct XccyMarket {
    /// JPY の OIS。
    pub tona: Case,
    /// USD の OIS。
    pub sofr: Case,
}

/// JPY 現金担保（TONA 付利、スプレッド 0）。
pub(crate) fn jpy_cash() -> Result<CollateralKey> {
    Ok(CollateralKey::Cash(CashCollateralTerms {
        ccy: Currency::JPY,
        remuneration: Remuneration {
            index: "JPY-TONA".parse()?,
            spread: Spread::new(0.0)?,
            floor: None,
        },
    }))
}

fn write_error(error: std::fmt::Error) -> Error {
    Error::Internal(error.to_string())
}

impl XccyMarket {
    /// `QuantLib` の fixture の TONA と SOFR を読む。
    pub(crate) fn new() -> Result<Self> {
        let [tona, sofr] = cases()?;
        Ok(Self { tona, sofr })
    }

    /// 参照データ（FX fixing と XCCY の慣行を含む）と、JPY 現金担保の割引の行を持つ入力。
    pub(crate) fn inputs(&self, extra_config: &str) -> Result<Inputs> {
        let mut config = String::new();
        for (ccy, curve) in [("JPY", "JPY-TONA"), ("USD", IMPLIED)] {
            writeln!(
                config,
                "[[discounting]]\ncurrency = \"{ccy}\"\ncollateral = {{ cash = {{ ccy = \"JPY\", index = \"JPY-TONA\", spread = \"0bp\" }} }}\ncurve = \"{curve}\""
            )
            .map_err(write_error)?;
        }
        config.push_str(extra_config);
        Inputs::with_forecasts(
            &[&self.tona, &self.sofr],
            self.tona.as_of,
            &config,
            |parts, _| {
                let index: IndexId = FX_FIXING.parse()?;
                let calendar = CalendarSpec::Union(vec![
                    CalendarSpec::Base("JPY".parse()?),
                    CalendarSpec::Base("USD".parse()?),
                ]);
                let clock = LocalEventTime {
                    time: "16:00:00".parse()?,
                    zone: "Etc/UTC".parse()?,
                };
                parts.indices.insert(
                    index.clone(),
                    IndexDefinition {
                        id: index.clone(),
                        kind: IndexKind::FxFixing {
                            pair: "USDJPY".parse()?,
                        },
                        calendar: calendar.clone(),
                        publication: PublicationRule::Daily(DailyPublication {
                            business_days: 0,
                            calendar,
                            time: clock.time,
                            zone: clock.zone.clone(),
                        }),
                        observation_time: Some(clock),
                        status: IndexStatus::Active,
                    },
                );
                let leg = |case: &Case, fx_reset: Option<IndexId>| -> Result<SwapLegConvention> {
                    let Some(Convention::Ois(ois)) = parts.conventions.get(&case.index.parse()?)
                    else {
                        return Err(Error::Internal("missing OIS test convention".into()));
                    };
                    let mut periods = ois.fixed.clone();
                    periods.frequency = "3M".parse()?;
                    Ok(SwapLegConvention {
                        currency: case.currency,
                        periods,
                        payment: ois.payment.clone(),
                        coupon: SwapCouponConvention::Compounded(ois.float.clone()),
                        fx_reset,
                    })
                };
                // 払い脚は JPY の TONA（ベーシスを付ける）、受け脚は USD の SOFR（元本をリセットする）。
                let legs = [leg(&self.tona, None)?, leg(&self.sofr, Some(index))?];
                let Some(Convention::Ois(ois)) = parts.conventions.get(&self.tona.index.parse()?)
                else {
                    return Err(Error::Internal("missing OIS test convention".into()));
                };
                let spot_lag = ois.spot_lag;
                parts.conventions.insert(
                    XCCY.parse()?,
                    Convention::XccySwap(Box::new(SwapConvention { spot_lag, legs })),
                );
                Ok(())
            },
        )
    }

    /// 外貨担保曲線の ID。
    pub(crate) fn implied() -> Result<CurveId> {
        IMPLIED.parse()
    }

    /// USDJPY の FX ポイントと長期の XCCY。`long_end` を差し替えると長期の行を変えられる。
    pub(crate) fn fx(&self, long_end: Vec<FxImpliedLongEnd>) -> Result<FxDefinition> {
        let [on, tn, three, six] = POINTS;
        Ok(FxDefinition {
            pair: "USDJPY".parse()?,
            direct: Some(DirectFx {
                market_collateral: "market".parse()?,
                quotes: FxQuotes {
                    spot: "USDJPY.SPOT".parse()?,
                    on: Some(on.0.parse()?),
                    tn: Some(tn.0.parse()?),
                    sn: None,
                    forwards: [three, six]
                        .iter()
                        .map(|(quote, date, _)| {
                            Ok(FxDatedQuote {
                                settlement: date
                                    .ok_or_else(|| Error::Internal("dated FX point".into()))?
                                    .parse()?,
                                quote: quote.parse()?,
                            })
                        })
                        .collect::<Result<_>>()?,
                },
                implied: Some(FxImpliedCurve {
                    id: Self::implied()?,
                    extrapolation: Extrapolation::default(),
                    long_end,
                }),
            }),
            collaterals: vec![FxCollateral {
                id: "market".parse()?,
                base: IMPLIED.parse()?,
                quote: self.tona.index.parse()?,
            }],
            via: None,
        })
    }

    /// 1Y と 18M の XCCY の長期の行。
    pub(crate) fn long_end() -> Result<Vec<FxImpliedLongEnd>> {
        XCCY_QUOTES
            .iter()
            .map(|(quote, tenor, _, pillar)| {
                Ok(FxImpliedLongEnd {
                    pillar: midnight(pillar.parse()?)?,
                    quote: quote.parse()?,
                    convention: XCCY.parse()?,
                    tenor: tenor.parse()?,
                    collateral: jpy_cash()?,
                })
            })
            .collect()
    }

    /// TONA と SOFR の OIS 曲線と USDJPY の FX からなる市場の定義。
    pub(crate) fn definition(&self) -> Result<MarketDefinition> {
        let uncollateralized = CollateralKey::Uncollateralized;
        definition(&[
            (&self.tona, &uncollateralized),
            (&self.sofr, &uncollateralized),
        ])?
        .with_fx([self.fx(Self::long_end()?)?])
    }

    /// OIS、FX と XCCY の quote。
    pub(crate) fn quotes(&self) -> Result<QuoteSet> {
        QuoteSet::new(
            fixture_quotes(&self.tona)?
                .into_iter()
                .chain(fixture_quotes(&self.sofr)?)
                .chain([("USDJPY.SPOT".parse()?, SPOT)])
                .chain(
                    POINTS
                        .iter()
                        .map(|(quote, _, value)| Ok::<_, Error>((quote.parse()?, *value)))
                        .collect::<Result<Vec<_>>>()?,
                )
                .chain(
                    XCCY_QUOTES
                        .iter()
                        .map(|(quote, _, value, _)| Ok::<_, Error>((quote.parse()?, *value)))
                        .collect::<Result<Vec<_>>>()?,
                ),
        )
    }

    /// 長期の XCCY のテンプレート。
    pub(crate) fn template(tenor: &str) -> Result<InstrumentTemplate> {
        Ok(InstrumentTemplate::XccyBasis {
            convention: XCCY.parse()?,
            tenor: tenor.parse::<TenorLength>()?,
        })
    }

    /// 受渡日の 16:00 UTC（ペアの慣行の受渡時刻）。
    pub(crate) fn settlement(date: &str) -> Result<Instant> {
        format!("{date}T16:00:00Z").parse()
    }

    /// 通貨ペア。
    pub(crate) fn pair() -> Result<CcyPair> {
        "USDJPY".parse()
    }
}
`````

`````
crates/rhodium-pricing/tests/forecast_measure.rs
//! 予想と割引の担保の違い（`CrossCollateralProjection`）を、予想曲線の予想の測度で判定する。
//!
//! 予想曲線は指標ごとに一本で、その予想は曲線を当てる校正商品の担保の下で定まる
//! （設計書「担保による割引」）。同じ担保で割り引く取引は、予想曲線と割引曲線の ID が違っても
//! 予想の担保の調整を要しない。担保の違う取引だけが決定論的なベーシスの仮定に頼る。

use std::sync::Arc;

use rhodium_core::{Error, Result, Warning, currency::Currency, money::Money, quantity::Spread};
use rhodium_market::{
    curve::{CurveExpr, CurveId, CurveRole, CurveSpec},
    structure::MarketStructure,
};
use rhodium_pricing::{
    calibration::{CalibratedMarket, GroupOutcome, calibrate},
    context::Context,
    prepare::{PrepareInput, Prepared, prepare},
    resolve::CollateralKey,
};
use rhodium_products::{
    agreement::{CashCollateralTerms, Remuneration},
    standard::{Fra, Position, StandardProduct},
};

#[allow(
    dead_code,
    reason = "共有の fixture のうち、担保の判定に使わない関数がある"
)]
#[path = "fixtures/ois.rs"]
pub mod ois;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、担保の判定に使わない関数がある"
)]
#[path = "fixtures/multicurve.rs"]
pub mod multicurve;

use multicurve::{MulticurveCase, multicurve_cases};
use ois::{Inputs, prepare_ois};

/// 校正商品と同じ担保（無担保の行で TONA に割り引く）。
const CALIBRATION: CollateralKey = CollateralKey::Uncollateralized;

/// 25bp の付利スプレッドの JPY 現金担保で割り引く曲線（TONA + 25bp）。
const SPREAD_CURVE: &str = "JPY-TONA-25";

fn cash(ccy: Currency, index: &str, spread: f64) -> Result<CollateralKey> {
    Ok(CollateralKey::Cash(CashCollateralTerms {
        ccy,
        remuneration: Remuneration {
            index: index.parse()?,
            spread: Spread::new(spread)?,
            floor: None,
        },
    }))
}

/// 担保条件の違う割引曲線 `L(t) = L_base(t) − s·t`。校正しないので θ を持たない。
fn spread_curve(id: &str, ccy: Currency, base: &str, spread: f64) -> Result<CurveSpec> {
    CurveSpec::from_expr(
        id.parse()?,
        ccy,
        [CurveRole::Discount],
        CurveExpr::Sum(vec![CurveExpr::Ref(base.parse()?), CurveExpr::Flat(spread)]),
    )
}

fn calibrated(
    case: &MulticurveCase,
    extra: Vec<CurveSpec>,
    ctx: &Context<'_>,
) -> Result<CalibratedMarket> {
    let result = calibrate(&case.definition_with(extra)?, &case.quotes()?, ctx)?;
    for group in result.groups() {
        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
            return Err(Error::Calibration(format!("{group:?}")));
        }
    }
    Ok(result)
}

fn basis_trade(
    case: &MulticurveCase,
    collateral: &CollateralKey,
    market: &CalibratedMarket,
    ctx: &Context<'_>,
) -> Result<Prepared> {
    let node = &case.nodes[1];
    let product = case.template(node)?.build(node.rate + 0.001, ctx.refdata)?;
    prepare(
        PrepareInput::Standard {
            product: &product,
            trade_date: case.overnight.as_of,
        },
        collateral,
        &[],
        market.market().structure(),
        ctx,
    )
}

fn has(warnings: &[Warning], expected: &Warning) -> bool {
    warnings.contains(expected)
}

#[test]
fn tibor_swap_under_calibration_collateral_needs_no_cross_collateral_convention() -> Result<()> {
    // TIBOR 6M は TONA 担保のベーシスで校正するので、予想の測度は JPY-TONA。
    // 同じ担保で割り引く TIBOR の取引と校正商品は、予想曲線と割引曲線の ID が違っても
    // 担保の調整を要さず、警告を出さない。
    let [jpy, _] = multicurve_cases()?;
    let inputs = jpy.inputs()?;
    let ctx = inputs.context();
    let market = calibrated(&jpy, vec![], &ctx)?;
    let term: CurveId = jpy.term_index.parse()?;
    let overnight: CurveId = jpy.overnight.index.parse()?;
    assert_eq!(
        market.market().structure().forecast_measure(&term),
        Some(&overnight)
    );
    assert!(!has(market.warnings(), &Warning::CrossCollateralProjection));
    let trade = basis_trade(&jpy, &CALIBRATION, &market, &ctx)?;
    assert!(!has(trade.warnings(), &Warning::CrossCollateralProjection));

    // 近似を禁止しても、厳密に評価できる TIBOR のグループは解け、取引も準備できる。
    let strict = multicurve_inputs(
        &jpy,
        "[conventions]\ndisabled = [\"cross_collateral_projection\"]\n",
    )?;
    let strict_ctx = strict.context();
    let strict_market = calibrated(&jpy, vec![], &strict_ctx)?;
    basis_trade(&jpy, &CALIBRATION, &strict_market, &strict_ctx)?;
    Ok(())
}

#[test]
fn tibor_swap_under_a_spread_csa_records_cross_collateral_projection() -> Result<()> {
    // 25bp の付利の CSA は TONA + 25bp の曲線で割り引くので、TONA 担保で定まる TIBOR の予想と
    // 担保が違う。この取引だけが CrossCollateralProjection を記録し、禁止すると準備できない。
    let [jpy, _] = multicurve_cases()?;
    let rule = format!(
        "[[discounting]]\ncurrency = \"JPY\"\ncollateral = {{ cash = {{ ccy = \"JPY\", index = \"{}\", spread = \"25bp\" }} }}\ncurve = \"{SPREAD_CURVE}\"\n",
        jpy.overnight.index
    );
    let extra = || spread_curve(SPREAD_CURVE, Currency::JPY, &jpy.overnight.index, 0.0025);
    let inputs = multicurve_inputs(&jpy, &rule)?;
    let ctx = inputs.context();
    let market = calibrated(&jpy, vec![extra()?], &ctx)?;
    let collateral = cash(Currency::JPY, &jpy.overnight.index, 0.0025)?;
    let trade = basis_trade(&jpy, &collateral, &market, &ctx)?;
    assert!(has(trade.warnings(), &Warning::CrossCollateralProjection));
    // 校正商品の担保の取引は、同じ市場でも警告しない。
    let native = basis_trade(&jpy, &CALIBRATION, &market, &ctx)?;
    assert!(!has(native.warnings(), &Warning::CrossCollateralProjection));

    let strict = multicurve_inputs(
        &jpy,
        &format!("{rule}[conventions]\ndisabled = [\"cross_collateral_projection\"]\n"),
    )?;
    let strict_ctx = strict.context();
    let strict_market = calibrated(&jpy, vec![extra()?], &strict_ctx)?;
    assert!(matches!(
        basis_trade(&jpy, &collateral, &strict_market, &strict_ctx),
        Err(Error::Unsupported(_))
    ));
    Ok(())
}

#[test]
fn sofr_leg_under_jpy_collateral_records_cross_collateral_projection() -> Result<()> {
    // SOFR の予想は USD 担保の OIS で定まる。JPY 現金担保の USD の支払いは別の USD 曲線で
    // 割り引くので、SOFR の複利の予想と割引の担保が違う。
    let [_, usd] = multicurve_cases()?;
    let csa = "USD-JPYCSA";
    let rule = format!(
        "[[discounting]]\ncurrency = \"USD\"\ncollateral = {{ cash = {{ ccy = \"JPY\", index = \"JPY-TONA\", spread = \"0bp\" }} }}\ncurve = \"{csa}\"\n"
    );
    let inputs = multicurve_inputs(&usd, &rule)?;
    let ctx = inputs.context();
    let market = calibrated(
        &usd,
        vec![spread_curve(
            csa,
            Currency::USD,
            &usd.overnight.index,
            0.001,
        )?],
        &ctx,
    )?;
    let node = &usd.overnight.nodes[2];
    let foreign = cash(Currency::JPY, "JPY-TONA", 0.0)?;
    let trade = prepare_ois(
        &usd.overnight,
        node,
        node.rate,
        &foreign,
        market.market(),
        &ctx,
    )?;
    assert!(has(trade.warnings(), &Warning::CrossCollateralProjection));
    let native = prepare_ois(
        &usd.overnight,
        node,
        node.rate,
        &CALIBRATION,
        market.market(),
        &ctx,
    )?;
    assert!(!has(native.warnings(), &Warning::CrossCollateralProjection));
    Ok(())
}

#[test]
fn fra_under_calibration_collateral_records_only_fra_discounting() -> Result<()> {
    // FRA の決済は期首に P_d(s) で割り引くので、予想曲線と割引曲線が違えば同じ担保でも
    // FraDiscounting の凸性を無視する。担保は同じなので CrossCollateralProjection は記録しない。
    let [jpy, _] = multicurve_cases()?;
    let inputs = jpy.inputs()?;
    let ctx = inputs.context();
    let market = calibrated(&jpy, vec![], &ctx)?;
    let product = StandardProduct::Fra(Fra::new(
        Money::new(Currency::JPY, 100_000_000.0),
        jpy.term_index.parse()?,
        jpy.nodes[1].date,
        0.01,
        Position::Long,
    )?);
    let trade = prepare(
        PrepareInput::Standard {
            product: &product,
            trade_date: jpy.overnight.as_of,
        },
        &CALIBRATION,
        &[],
        market.market().structure(),
        &ctx,
    )?;
    assert!(has(trade.warnings(), &Warning::FraDiscounting));
    assert!(!has(trade.warnings(), &Warning::CrossCollateralProjection));
    Ok(())
}

#[test]
fn forecast_measure_is_unknown_when_bindings_disagree() -> Result<()> {
    // TIBOR のベーシスの一つを 25bp の CSA で校正すると、予想を定める担保が一つに決まらない。
    // 測度を記録せず、予想曲線と割引曲線が同じかどうかの比較に戻るので、TONA で割り引く
    // 取引も CrossCollateralProjection を記録する。
    let [jpy, _] = multicurve_cases()?;
    let rule = format!(
        "[[discounting]]\ncurrency = \"JPY\"\ncollateral = {{ cash = {{ ccy = \"JPY\", index = \"{}\", spread = \"25bp\" }} }}\ncurve = \"{SPREAD_CURVE}\"\n",
        jpy.overnight.index
    );
    let inputs = multicurve_inputs(&jpy, &rule)?;
    let ctx = inputs.context();
    let collateral = cash(Currency::JPY, &jpy.overnight.index, 0.0025)?;
    let definition = jpy.definition_with_collateral(
        [spread_curve(
            SPREAD_CURVE,
            Currency::JPY,
            &jpy.overnight.index,
            0.0025,
        )?],
        |index| (index == 0).then(|| collateral.clone()),
    )?;
    let market = calibrate(&definition, &jpy.quotes()?, &ctx)?;
    let term: CurveId = jpy.term_index.parse()?;
    assert_eq!(market.market().structure().forecast_measure(&term), None);
    let trade = basis_trade(&jpy, &CALIBRATION, &market, &ctx)?;
    assert!(has(trade.warnings(), &Warning::CrossCollateralProjection));
    Ok(())
}

#[test]
fn market_structure_rejects_inconsistent_forecast_measures() -> Result<()> {
    // 測度は同じ通貨の割引曲線でなければならず、予想の役割を持たない曲線には付けられない。
    let [jpy, usd] = multicurve_cases()?;
    let as_of = jpy.overnight.as_of;
    let structure = || {
        MarketStructure::new(
            ois::midnight(as_of)?,
            [
                jpy.overnight.curve()?,
                usd.overnight.curve()?,
                spread_curve(SPREAD_CURVE, Currency::JPY, &jpy.overnight.index, 0.0025)?,
            ],
        )
    };
    let tona: CurveId = jpy.overnight.index.parse()?;
    let sofr: CurveId = usd.overnight.index.parse()?;
    let spread: CurveId = SPREAD_CURVE.parse()?;
    for measures in [
        vec![(tona.clone(), sofr.clone())],
        vec![(spread.clone(), tona.clone())],
        vec![("UNKNOWN".parse()?, tona.clone())],
    ] {
        assert!(matches!(
            structure()?.with_forecast_measures(measures),
            Err(Error::InvalidInput(_))
        ));
    }
    let recorded = Arc::new(structure()?.with_forecast_measures([(tona.clone(), spread.clone())])?);
    assert_eq!(recorded.forecast_measure(&tona), Some(&spread));
    // 測度を付けた予想は測度の曲線と同じ担保で、予想曲線自身とは違う担保になる。
    let projection = recorded.resolve_projection(&tona)?;
    assert!(projection.is_forecast_under(&recorded.resolve_discount(&spread)?));
    assert!(!projection.is_forecast_under(&recorded.resolve_discount(&tona)?));
    // 記録の無い構造では、同じ曲線だけを同じ担保とみなす。
    let plain = Arc::new(structure()?);
    let projection = plain.resolve_projection(&tona)?;
    assert!(projection.is_forecast_under(&plain.resolve_discount(&tona)?));
    assert!(!projection.is_forecast_under(&plain.resolve_discount(&spread)?));
    Ok(())
}

/// multicurve の入力に、割引の行や近似の設定を加える。
fn multicurve_inputs(case: &MulticurveCase, extra: &str) -> Result<Inputs> {
    case.inputs_with(extra)
}
`````

`````
crates/rhodium-pricing/tests/fx_curves.rs
//! 担保別 FX フォワードを、実際の受渡と手計算のアウトライトで照合する。

use std::ptr;

use rhodium_core::{
    Error, Result, Warning,
    calendar::{BaseCalendar, Extrapolation as CalendarExtrapolation, WeekendRule},
    currency::Currency,
    date::Weekday,
    fixing::{FixingStore, FixingStoreParts},
    fx::CcyPair,
    refdata::{CurrencyDefinition, Deliverability, FxConvention, FxSpotLag, RefData, RefDataParts},
    rounding::{Rounding, RoundingMode},
    time::{Instant, LocalEventTime, UtcOffset, ZoneRules, model_time},
};
use rhodium_market::{
    curve::{CurveExpr, CurveRole, CurveSpec, Extrapolation},
    fx::FxCollateral,
    quote::{QuoteId, QuoteSet},
};
use rhodium_pricing::{
    calibration::{CalibratedMarket, GroupOutcome, calibrate},
    config::PricingConfig,
    context::Context,
    definition::{
        DirectFx, FxDatedQuote, FxDefinition, FxImpliedCurve, FxQuotes, MarketDefinition,
    },
};
use rhodium_testkit::assert_close;

const BEFORE: &str = "2026-01-02T09:00:00Z";
const JUST_BEFORE: &str = "2026-01-02T15:59:59Z";
const AT_CUTOFF: &str = "2026-01-02T16:00:00Z";
const AFTER: &str = "2026-01-02T16:00:01Z";
const USDJPY_TOMORROW: &str = "2026-01-07T16:00:00Z";
const USDJPY_SPOT: &str = "2026-01-08T16:00:00Z";
const USDCAD_SPOT: &str = "2026-01-06T16:00:00Z";
const USDCAD_SN: &str = "2026-01-07T16:00:00Z";
const LONG: &str = "2026-03-02T16:00:00Z";
const REL: f64 = 1e-12;
const ABS: f64 = 1e-14;

fn clock() -> Result<LocalEventTime> {
    Ok(LocalEventTime {
        time: "16:00:00".parse()?,
        zone: "Etc/UTC".parse()?,
    })
}

struct Inputs {
    refdata: RefData,
    fixings: FixingStore,
    config: PricingConfig,
}

/// 既定のペアの慣行（USDJPY は T+2、USDCAD と CADJPY は T+1）。
const CONVENTIONS: [(&str, FxSpotLag, u8); 3] = [
    ("USDJPY", FxSpotLag::TwoBusinessDays, 2),
    ("USDCAD", FxSpotLag::OneBusinessDay, 4),
    ("CADJPY", FxSpotLag::OneBusinessDay, 2),
];

impl Inputs {
    fn new() -> Result<Self> {
        Self::with_conventions(&CONVENTIONS)
    }

    fn with_conventions(conventions: &[(&str, FxSpotLag, u8)]) -> Result<Self> {
        let clock = clock()?;
        let mut parts = RefDataParts::default();
        parts.zones.insert(
            clock.zone.clone(),
            ZoneRules::new(
                "2025-12-01T00:00:00Z".parse()?..="2027-01-31T00:00:00Z".parse()?,
                UtcOffset::from_seconds(0)?,
                Vec::new(),
            )?,
        );
        // USD と JPY の異なる休日が USDJPY のスポットを木曜日にずらす。
        for (code, holidays) in [
            ("USD", vec!["2026-01-05".parse()?]),
            ("JPY", vec!["2026-01-06".parse()?]),
            ("CAD", vec![]),
        ] {
            let currency: Currency = code.parse()?;
            let calendar = code.parse()?;
            parts.calendars.insert(
                calendar,
                BaseCalendar::new(
                    "2026-01-01".parse()?..="2026-12-31".parse()?,
                    WeekendRule::new([Weekday::Sat, Weekday::Sun]),
                    holidays,
                    CalendarExtrapolation::Strict,
                )?,
            );
            parts.currencies.insert(
                currency,
                CurrencyDefinition {
                    code: currency,
                    minor_units: 2,
                    payment_rounding: Rounding::new(0.01, RoundingMode::HalfEven)?,
                    deliverability: Deliverability::Deliverable,
                    settlement_time: clock.clone(),
                    settlement_calendar: code.parse()?,
                },
            );
        }
        // スポットラグ、同日受渡の締めと pip の桁は参照データのペアの慣行が決める。
        for &(pair, spot_lag, pip_exponent) in conventions {
            let pair: CcyPair = pair.parse()?;
            parts.fx_conventions.insert(
                pair,
                FxConvention {
                    pair,
                    spot_lag,
                    settlement: clock.clone(),
                    pip_exponent,
                },
            );
        }
        let refdata = RefData::new(parts)?;
        let fixings = FixingStore::new(&refdata, FixingStoreParts::default())?;
        let config = toml::from_str("version = 1\nreporting_currency = 'USD'")
            .map_err(|error| Error::InvalidData(error.to_string()))?;
        Ok(Self {
            refdata,
            fixings,
            config,
        })
    }

    fn calibrate(
        &self,
        def: &MarketDefinition,
        quotes: &QuoteSet,
        at: &str,
    ) -> Result<CalibratedMarket> {
        calibrate(
            def,
            quotes,
            &Context {
                as_of: at.parse()?,
                refdata: &self.refdata,
                fixings: &self.fixings,
                config: &self.config,
            },
        )
    }
}

fn collateral(id: &str, base: &str, quote: &str) -> Result<FxCollateral> {
    Ok(FxCollateral {
        id: id.parse()?,
        base: base.parse()?,
        quote: quote.parse()?,
    })
}

fn pairs() -> Result<[FxDefinition; 2]> {
    let usd_jpy = FxDefinition {
        pair: "USDJPY".parse()?,
        direct: Some(DirectFx {
            market_collateral: "market".parse()?,
            quotes: FxQuotes {
                spot: "USDJPY.SPOT".parse()?,
                on: Some("USDJPY.ON".parse()?),
                tn: Some("USDJPY.TN".parse()?),
                sn: Some("USDJPY.SN".parse()?),
                forwards: vec![FxDatedQuote {
                    settlement: "2026-03-02".parse()?,
                    quote: "USDJPY.FWD".parse()?,
                }],
            },
            implied: Some(FxImpliedCurve {
                id: "USD-under-JPY".parse()?,
                extrapolation: Extrapolation::default(),
                long_end: vec![],
            }),
        }),
        collaterals: vec![
            collateral("market", "USD-under-JPY", "JPY")?,
            collateral("other", "USD", "JPY-under-USD")?,
        ],
        via: None,
    };
    let usd_cad = FxDefinition {
        pair: "USDCAD".parse()?,
        direct: Some(DirectFx {
            market_collateral: "market".parse()?,
            quotes: FxQuotes {
                spot: "USDCAD.SPOT".parse()?,
                on: Some("USDCAD.ON".parse()?),
                tn: None,
                sn: Some("USDCAD.SN".parse()?),
                forwards: vec![],
            },
            implied: Some(FxImpliedCurve {
                id: "CAD-under-USD".parse()?,
                extrapolation: Extrapolation::default(),
                long_end: vec![],
            }),
        }),
        collaterals: vec![collateral("market", "USD", "CAD-under-USD")?],
        via: None,
    };
    Ok([usd_jpy, usd_cad])
}

fn definition(pairs: impl IntoIterator<Item = FxDefinition>) -> Result<MarketDefinition> {
    let curves = [
        ("USD", "USD", 0.04),
        ("JPY", "JPY", 0.01),
        ("CAD", "CAD", 0.025),
        ("JPY-under-USD", "JPY", 0.008),
    ]
    .into_iter()
    .map(|(id, ccy, rate)| {
        CurveSpec::from_expr(
            id.parse()?,
            ccy.parse()?,
            [CurveRole::Discount],
            CurveExpr::Flat(rate),
        )
    })
    .collect::<Result<Vec<_>>>()?;
    MarketDefinition::new("USD".parse()?, "Etc/UTC".parse()?, curves, [])?.with_fx(pairs)
}

fn quotes(extra: &[(&str, f64)]) -> Result<QuoteSet> {
    QuoteSet::new(
        [
            ("USDJPY.SPOT", 150.0),
            ("USDJPY.ON", 2.0),
            ("USDJPY.TN", -1.0),
            ("USDJPY.SN", 3.0),
            ("USDJPY.FWD", 27.0),
            ("USDCAD.SPOT", 1.300),
            ("USDCAD.ON", 10.0),
            ("USDCAD.SN", 20.0),
        ]
        .into_iter()
        .chain(extra.iter().copied())
        .map(|(id, value)| Ok((id.parse()?, value)))
        .collect::<Result<Vec<_>>>()?,
    )
}

fn assert_solved(market: &CalibratedMarket) -> Result<()> {
    for group in market.groups() {
        let GroupOutcome::Solved(solution) = group.outcome() else {
            return Err(Error::Calibration(format!("{group:?}")));
        };
        assert!(solution.residual_norm() <= 1e-14);
    }
    Ok(())
}

fn hand_x0(as_of: Instant, start: &str, end: &str, first: f64, second: f64) -> Result<f64> {
    let start = model_time(as_of, start.parse()?).years();
    let end = model_time(as_of, end.parse()?).years();
    Ok(first * (-(second / first).ln() * start / (end - start)).exp())
}

#[test]
fn fx_x0_short_dated_chain() -> Result<()> {
    // T+2 は TN と ON、T+1 は ON。締め後は同日受渡の ON が不要になる。
    let inputs = Inputs::new()?;
    let def = definition(pairs()?)?;
    let quotes = quotes(&[])?;
    for as_of in [BEFORE, JUST_BEFORE, AT_CUTOFF, AFTER] {
        let result = inputs.calibrate(&def, &quotes, as_of)?;
        assert_solved(&result)?;
        let now = as_of.parse()?;
        result.market().with_view(|view| {
            let jpy: CcyPair = "USDJPY".parse()?;
            let cad: CcyPair = "USDCAD".parse()?;
            let market = "market".parse()?;
            let other = "other".parse()?;
            let jpy_market = view.fx(jpy, &market)?;
            let jpy_other = view.fx(jpy, &other)?;
            let cad_market = view.fx(cad, &market)?;
            let (jpy_x0, cad_x0) = if as_of == AFTER {
                (
                    hand_x0(now, USDJPY_TOMORROW, USDJPY_SPOT, 150.01, 150.0)?,
                    hand_x0(now, USDCAD_SPOT, USDCAD_SN, 1.300, 1.302)?,
                )
            } else {
                assert_close!(
                    jpy_market.forward(AT_CUTOFF.parse()?)?.rate(),
                    149.99,
                    REL,
                    ABS
                );
                assert_close!(
                    cad_market.forward(AT_CUTOFF.parse()?)?.rate(),
                    1.299,
                    REL,
                    ABS
                );
                (
                    hand_x0(now, AT_CUTOFF, USDJPY_TOMORROW, 149.99, 150.01)?,
                    hand_x0(now, AT_CUTOFF, USDCAD_SPOT, 1.299, 1.300)?,
                )
            };
            assert_close!(jpy_market.today().rate(), jpy_x0, REL, ABS);
            assert_close!(jpy_other.today().rate(), jpy_x0, REL, ABS);
            assert_close!(cad_market.today().rate(), cad_x0, REL, ABS);
            for (at, expected) in [
                (USDJPY_TOMORROW, 150.01),
                (USDJPY_SPOT, 150.0),
                ("2026-01-09T16:00:00Z", 150.03),
                (LONG, 150.27),
            ] {
                assert_close!(jpy_market.forward(at.parse()?)?.rate(), expected, REL, ABS);
            }
            assert_close!(
                cad_market.forward(USDCAD_SPOT.parse()?)?.rate(),
                1.300,
                REL,
                ABS
            );
            assert_close!(
                cad_market.forward(USDCAD_SN.parse()?)?.rate(),
                1.302,
                REL,
                ABS
            );
            let at = LONG.parse()?;
            let point = result
                .market()
                .structure()
                .prepare_fx_point(jpy, &other, at)?;
            assert_close!(
                jpy_other.bind_forward(&point)?.rate(),
                jpy_x0 * (-0.032 * model_time(now, at).years()).exp(),
                REL,
                ABS
            );
            assert!(matches!(
                jpy_market.bind_forward(&point),
                Err(Error::InvalidInput(_))
            ));
            Ok(())
        })?;
    }
    Ok(())
}

#[test]
fn fx_x0_on_tn_consistency() -> Result<()> {
    // 短期の quote が揃う定義には、スポットを含む各区間を再現する自由度が必要。
    // 締めの前後で成否が変わらないよう、評価時点を受け取る前の定義の作成で拒否する。
    let [mut jpy, mut cad] = pairs()?;
    let direct = jpy
        .direct
        .as_mut()
        .ok_or_else(|| Error::Internal("test needs direct FX".into()))?;
    direct.implied = None;
    jpy.collaterals[0].base = "USD".parse()?;
    assert!(
        matches!(definition([jpy]), Err(Error::InvalidData(message)) if message.contains("require an FX-implied curve"))
    );

    // T+1 のペアで外貨担保曲線を ON だけに校正すると、締め後に校正区間が残らない。
    // 午前は校正でき午後に校正全体が止まる定義を、作るときに拒否する。
    let direct = cad
        .direct
        .as_mut()
        .ok_or_else(|| Error::Internal("test needs direct FX".into()))?;
    direct.quotes.sn = None;
    assert!(
        matches!(definition([cad]), Err(Error::InvalidData(message)) if message.contains("same-day cutoff"))
    );
    Ok(())
}

#[test]
fn fx_cross_direct_vs_triangulated() -> Result<()> {
    // CADJPY は T+1 の逆向きと T+2 を経由する。直接のマークがあれば警告を出して採用する。
    let [jpy, cad] = pairs()?;
    let inputs = Inputs::new()?;
    let pair: CcyPair = "CADJPY".parse()?;
    let mut cross = FxDefinition {
        pair,
        direct: None,
        collaterals: vec![collateral("market", "CAD", "JPY")?],
        via: Some(Currency::USD),
    };
    let triangle = inputs.calibrate(
        &definition([jpy.clone(), cad.clone(), cross.clone()])?,
        &quotes(&[])?,
        BEFORE,
    )?;
    assert_solved(&triangle)?;
    let path = triangle.market().with_view(|view| {
        let usdjpy = view.fx_today("USDJPY".parse()?)?.rate();
        let usdcad = view.fx_today("USDCAD".parse()?)?.rate();
        let triangulated = view.fx_today(pair)?.rate();
        assert_close!(triangulated, usdjpy / usdcad, REL, ABS);
        view.ccy_scope(
            "CAD".parse()?,
            view.discount_curve(&"CAD".parse()?)?,
            |cad_scope| {
                view.ccy_scope(
                    Currency::JPY,
                    view.discount_curve(&"JPY".parse()?)?,
                    |jpy_scope| {
                        let value = jpy_scope.convert_from(&cad_scope, cad_scope.present(10.0))?;
                        assert_close!(
                            jpy_scope.into_money(value).amount(),
                            10.0 * triangulated,
                            REL,
                            ABS
                        );
                        Ok(())
                    },
                )
            },
        )?;
        Ok(triangulated)
    })?;
    cross.direct = Some(DirectFx {
        market_collateral: "market".parse()?,
        quotes: FxQuotes {
            spot: "CADJPY.SPOT".parse()?,
            on: None,
            tn: None,
            sn: None,
            forwards: vec![],
        },
        implied: None,
    });
    let direct = inputs.calibrate(
        &definition([jpy, cad, cross])?,
        &quotes(&[("CADJPY.SPOT", 111.0)])?,
        BEFORE,
    )?;
    assert!(
        direct
            .warnings()
            .contains(&Warning::FxCrossInconsistent { pair })
    );
    direct.market().with_view(|view| {
        let time = model_time(BEFORE.parse()?, "2026-01-05T16:00:00Z".parse()?).years();
        let actual = view.fx_today(pair)?.rate();
        assert_close!(actual, 111.0 * (0.015 * time).exp(), REL, ABS);
        assert!((actual - path).abs() > 1.0);
        Ok(())
    })
}

#[test]
fn fx_quote_update_matches_full_recalibration() -> Result<()> {
    // 円の TN の更新は円のマークと外貨担保曲線を更新する。CAD の解は共有する。
    let inputs = Inputs::new()?;
    let def = definition(pairs()?)?;
    let original = quotes(&[])?;
    let result = inputs.calibrate(&def, &original, BEFORE)?;
    let id: QuoteId = "USDJPY.TN".parse()?;
    let bumped = QuoteSet::new(original.iter().map(|(quote, value)| {
        (
            quote.clone(),
            if *quote == id { value + 0.5 } else { value },
        )
    }))?;
    let updated = result.update(&bumped)?;
    let fresh = inputs.calibrate(&def, &bumped, BEFORE)?;
    assert_solved(&updated)?;
    for (old, new) in result.groups().iter().zip(updated.groups()) {
        let cad = old
            .blocks()
            .iter()
            .any(|block| block.curve().as_str().contains("CAD"));
        assert_eq!(ptr::eq(old.outcome(), new.outcome()), cad);
    }
    for (&actual, &expected) in updated.market().theta().iter().zip(fresh.market().theta()) {
        assert_close!(actual, expected, REL, ABS);
    }
    Ok(())
}

#[test]
fn fx_missing_active_short_quote_only_blocks_its_pair() -> Result<()> {
    // 当日受渡を含む場合、ON の欠損は円の FX だけに伝わる。締め後は同じ入力で校正できる。
    let inputs = Inputs::new()?;
    let def = definition(pairs()?)?;
    let missing: QuoteId = "USDJPY.ON".parse()?;
    let quotes = QuoteSet::new(
        quotes(&[])?
            .iter()
            .filter(|&(id, _)| *id != missing)
            .map(|(id, value)| (id.clone(), value)),
    )?;
    let before = inputs.calibrate(&def, &quotes, BEFORE)?;
    assert!(
        before
            .groups()
            .iter()
            .any(|group| matches!(group.outcome(), GroupOutcome::Failed(Error::MissingData(_))))
    );
    assert!(
        before
            .groups()
            .iter()
            .any(|group| matches!(group.outcome(), GroupOutcome::Unavailable))
    );
    before.market().with_view(|view| {
        assert!(matches!(
            view.fx_today("USDJPY".parse()?),
            Err(Error::MissingData(_))
        ));
        assert!(view.fx_today("USDCAD".parse()?)?.rate() > 0.0);
        Ok(())
    })?;
    assert_solved(&inputs.calibrate(&def, &quotes, AFTER)?)
}

#[test]
fn fx_x0_chain_follows_refdata_spot_lag() -> Result<()> {
    // スポットまでの短期の連鎖は参照データのペアの慣行で決まる。USDCAD に TN を足すと、
    // T+1 の慣行ではスポットの後の区間になり連鎖にならないので、校正の入口で InvalidData になる。
    // 市場の定義を変えずに参照データの慣行だけを T+2 にすると、同じ定義で校正できる。
    let [_, mut cad] = pairs()?;
    let direct = cad
        .direct
        .as_mut()
        .ok_or_else(|| Error::Internal("test needs direct FX".into()))?;
    direct.quotes.tn = Some("USDCAD.TN".parse()?);
    let definition = definition([cad])?;
    let quotes = quotes(&[("USDCAD.TN", 5.0)])?;
    assert!(matches!(
        Inputs::new()?.calibrate(&definition, &quotes, BEFORE),
        Err(Error::InvalidData(message)) if message.contains("chain to spot")
    ));
    let t_plus_two = Inputs::with_conventions(&[("USDCAD", FxSpotLag::TwoBusinessDays, 4)])?;
    let market = t_plus_two.calibrate(&definition, &quotes, BEFORE)?;
    assert_solved(&market)?;
    Ok(())
}

#[test]
fn fx_definition_rejects_the_inverse_market_direction() -> Result<()> {
    // pip とポイントは慣行の標準の向き（USDJPY）で読むので、逆向きの JPYUSD で書いた市場の
    // quote は、ペアの慣行があっても校正の入口で拒否する。
    let [mut jpy, _] = pairs()?;
    jpy.pair = "JPYUSD".parse()?;
    let definition = definition([jpy])?;
    assert!(matches!(
        Inputs::new()?.calibrate(&definition, &quotes(&[])?, BEFORE),
        Err(Error::InvalidData(message)) if message.contains("convention direction USDJPY")
    ));
    Ok(())
}
`````

`````
crates/rhodium-pricing/tests/least_squares_calibration.rs
//! 最小二乗の校正グループ: 手計算の検算、当てはめの許容幅の警告、warm start、定義の検査、
//! 罰則付きの会合日の階段（R4）と Gauss–Newton の記録。

use std::sync::Arc;

use rhodium_core::{Error, Result, Warning, currency::Currency};
use rhodium_market::{
    curve::{CurveRole, CurveSpec, Extrapolation},
    market::Market,
    quote::{QuoteId, QuoteKind, QuoteSet},
    structure::{BlockId, ParamId},
};
use rhodium_pricing::{
    calibration::{CalibratedMarket, GroupOutcome, GroupSolution, HessianMode, calibrate},
    context::Context,
    definition::{
        CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
        StepPenalty,
    },
    prepare::evaluate,
    resolve::CollateralKey,
};
use rhodium_testkit::assert_close;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、最小二乗の検査に使わない関数がある"
)]
#[path = "fixtures/ois.rs"]
pub mod ois;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、最小二乗の検査に使わない関数がある"
)]
#[path = "fixtures/steps.rs"]
pub mod steps;

use ois::{Case, Inputs, cases, midnight, prepare_ois};
use steps::{OIS_TOLERANCE, PENALTY_SCALE, PENALTY_WEIGHT, PenalizedTona};

const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
const MARKS: [&str; 2] = ["MARK.1", "MARK.2"];

/// 一つのピラーの曲線の ln DF を、許容幅 σ の二つの直接マークで当てる市場の定義。
fn direct_marks(tona: &Case, sigmas: [f64; 2]) -> Result<MarketDefinition> {
    let id: rhodium_market::curve::CurveId = tona.index.parse()?;
    let curve = CurveSpec::new(
        tona.index.parse()?,
        Currency::JPY,
        [CurveRole::Discount, CurveRole::Projection],
        [midnight(tona.nodes[2].date)?],
        Extrapolation::Reject,
    )?;
    let bindings = MARKS
        .iter()
        .zip(sigmas)
        .map(|(quote, sigma)| {
            InstrumentBinding::new(
                quote.parse()?,
                QuoteKind::Mark,
                CalibrationTarget::DirectMark,
                QuoteTransform::Identity,
                COLLATERAL,
                ElementSelector::Param(ParamId::new(BlockId::nodes(id.clone()), 0)),
            )?
            .with_tolerance(sigma)
        })
        .collect::<Result<Vec<_>>>()?;
    MarketDefinition::new("JPY".parse()?, "Etc/UTC".parse()?, [curve], bindings)
}

fn marks(values: [f64; 2]) -> Result<QuoteSet> {
    QuoteSet::new(
        MARKS
            .iter()
            .zip(values)
            .map(|(quote, value)| Ok((quote.parse()?, value)))
            .collect::<Result<Vec<_>>>()?,
    )
}

fn solution(market: &CalibratedMarket, index: usize) -> Result<&GroupSolution> {
    match market.groups()[index].outcome() {
        GroupOutcome::Solved(solution) => Ok(solution),
        other => Err(Error::Calibration(format!("{other:?}"))),
    }
}

fn linear(error: impl std::fmt::Display) -> Error {
    Error::Numerical(error.to_string())
}

#[test]
fn weighted_direct_marks_solve_the_hand_least_squares() -> Result<()> {
    // 品質の手計算の例: min (θ − q₁)² + 4(θ − q₂)²（σ₁ = 1、σ₂ = 1/2）の解は θ = (q₁ + 4q₂)/5。
    // 停留条件 G = (θ − q₁) + 4(θ − q₂) の微分は ∂G/∂θ = 5、∂G/∂q = (−1, −4) なので、
    // quote への勾配は dθ/dq = −(∂G/∂θ)⁻¹·∂G/∂q = (1/5, 4/5)。
    let [tona, _] = cases()?;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let ctx = inputs.context();
    let [q1, q2] = [-0.01, -0.02];
    let result = calibrate(&direct_marks(&tona, [1.0, 0.5])?, &marks([q1, q2])?, &ctx)?;
    assert_eq!(result.groups().len(), 1);
    let solution = solution(&result, 0)?;
    assert_eq!(solution.hessian(), Some(HessianMode::Exact));
    let theta = result.market().theta()[0];
    assert_close!(theta, (q1 + 4.0 * q2) / 5.0, 1e-12, 1e-14);
    let inverse = solution
        .dg_dtheta()
        .solve(&[1.0][..])
        .map_err(linear)?
        .as_slice()[0];
    assert_close!(inverse, 0.2, 1e-10, 1e-12);
    // 列は quote の ID の昇順（MARK.1、MARK.2）。
    let dg_dq = [0, 1].map(|column| solution.dg_dq().get(0, column).map_err(linear));
    let [dg_dq1, dg_dq2] = [dg_dq[0].clone()?, dg_dq[1].clone()?];
    assert_close!(dg_dq1, -1.0, 1e-10, 1e-12);
    assert_close!(dg_dq2, -4.0, 1e-10, 1e-12);
    assert_close!(-inverse * dg_dq1, 0.2, 1e-10, 1e-12);
    assert_close!(-inverse * dg_dq2, 0.8, 1e-10, 1e-12);
    // 残差は重みを掛けない θ − q_i で、停止の誤差は停留点ではほぼ 0。
    assert_close!(solution.residuals()[0], theta - q1, 1e-12, 1e-14);
    assert_close!(solution.residuals()[1], theta - q2, 1e-12, 1e-14);
    assert!(solution.theta_error()[0].abs() <= 1e-15);
    Ok(())
}

#[test]
fn fit_outside_tolerance_warns_per_quote() -> Result<()> {
    // σ₁ = 1bp、σ₂ = 10bp で q₁ = −1%、q₂ = −2% を当てると、θ = (w₁²q₁ + w₂²q₂)/(w₁² + w₂²)
    // = −1.0099% で、残差は F₁ = −0.99bp（σ₁ の内）、F₂ = 99bp（σ₂ の外）。外れた quote だけを警告する。
    // 整合した quote に更新すると、再利用したグループの警告も解に合わせて消える。
    let [tona, _] = cases()?;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let ctx = inputs.context();
    let definition = direct_marks(&tona, [1e-4, 1e-3])?;
    let result = calibrate(&definition, &marks([-0.01, -0.02])?, &ctx)?;
    let outside: Vec<_> = result
        .warnings()
        .iter()
        .filter_map(|warning| match warning {
            Warning::CalibrationFitOutsideTolerance { quote } => Some(quote.as_str()),
            _ => None,
        })
        .collect();
    assert_eq!(outside, ["MARK.2"]);
    let consistent = result.update(&marks([-0.01, -0.01])?)?;
    assert!(
        !consistent
            .warnings()
            .iter()
            .any(|warning| matches!(warning, Warning::CalibrationFitOutsideTolerance { .. }))
    );
    Ok(())
}

#[test]
fn least_squares_update_warm_starts_a_strongly_convex_group() -> Result<()> {
    // 直接マークだけの最小二乗は θ の強凸な二次式で解が一意なので、前の θ から解き直す。
    // 更新の結果は、新しい quote から固定の初期値で解いた全体の校正と一致する。
    let [tona, _] = cases()?;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let ctx = inputs.context();
    let definition = direct_marks(&tona, [1.0, 0.5])?;
    let base = calibrate(&definition, &marks([-0.01, -0.02])?, &ctx)?;
    let quotes = marks([-0.01, -0.0195])?;
    let updated = base.update(&quotes)?;
    let full = calibrate(&definition, &quotes, &ctx)?;
    solution(&updated, 0)?;
    assert_close!(
        updated.market().theta()[0],
        full.market().theta()[0],
        1e-12,
        1e-14
    );
    assert_close!(
        updated.market().theta()[0],
        (-0.01 + 4.0 * -0.0195) / 5.0,
        1e-12,
        1e-14
    );
    Ok(())
}

#[test]
fn least_squares_definitions_reject_inconsistent_options() -> Result<()> {
    // 許容幅・罰則・Gauss–Newton の指定は、定義を作る境界と校正の入口で検査する。
    let [tona, _] = cases()?;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let ctx = inputs.context();
    let binding = || {
        tona.binding(
            &tona.nodes[0],
            &COLLATERAL,
            ElementSelector::Block(tona.block()?),
        )
    };
    for sigma in [0.0, -1e-4, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            binding()?.with_tolerance(sigma),
            Err(Error::InvalidInput(_))
        ));
    }
    for (weight, scale) in [(0.0, 1e-4), (1.0, f64::INFINITY), (-1.0, 1e-4)] {
        assert!(matches!(
            StepPenalty::new(tona.block()?, weight, scale),
            Err(Error::InvalidInput(_))
        ));
    }
    let square = ois::definition(&[(&tona, &COLLATERAL)])?;
    // 罰則は Steps の部品にだけ掛けられる。
    assert!(matches!(
        square
            .clone()
            .with_step_penalties([StepPenalty::new(tona.block()?, 1.0, 1e-4)?]),
        Err(Error::InvalidInput(message)) if message.contains("steps part")
    ));
    assert!(matches!(
        square
            .clone()
            .with_gauss_newton([BlockId::new(tona.index.parse()?, "missing")?]),
        Err(Error::InvalidInput(_))
    ));
    let quotes = QuoteSet::new(ois::fixture_quotes(&tona)?)?;
    // 正方系（許容幅の無い binding）に Gauss–Newton を指定するのは定義の誤り。
    let gauss_newton = square.clone().with_gauss_newton([tona.block()?])?;
    assert!(matches!(
        calibrate(&gauss_newton, &quotes, &ctx),
        Err(Error::InvalidInput(message)) if message.contains("no binding tolerances")
    ));
    // 許容幅を持つ binding と持たない binding が同じグループにあると、形が一つに決まらない。
    let block = tona.block()?;
    let mixed = tona
        .nodes
        .iter()
        .enumerate()
        .map(|(index, node)| {
            let element = u32::try_from(index).map_err(linear)?;
            let binding = tona.binding(
                node,
                &COLLATERAL,
                ElementSelector::Param(ParamId::new(block.clone(), element)),
            )?;
            if index == 0 {
                binding.with_tolerance(1e-4)
            } else {
                Ok(binding)
            }
        })
        .collect::<Result<Vec<_>>>()?;
    let mixed = MarketDefinition::new("JPY".parse()?, "Etc/UTC".parse()?, [tona.curve()?], mixed)?;
    assert!(matches!(
        calibrate(&mixed, &quotes, &ctx),
        Err(Error::InvalidInput(message)) if message.contains("mixes bindings")
    ));
    // 最小二乗でも残差と罰則の行が未知数より少なければ解が決まらない（5 要素を 4 本で当てる）。
    let short = tona.nodes[..4]
        .iter()
        .map(|node| {
            tona.binding(node, &COLLATERAL, ElementSelector::Block(block.clone()))?
                .with_tolerance(1e-4)
        })
        .collect::<Result<Vec<_>>>()?;
    let short = MarketDefinition::new("JPY".parse()?, "Etc/UTC".parse()?, [tona.curve()?], short)?;
    assert!(matches!(
        calibrate(&short, &quotes, &ctx),
        Err(Error::InvalidInput(message)) if message.contains("least-squares")
    ));
    Ok(())
}

fn penalized(gauss_newton: bool) -> Result<(PenalizedTona, Inputs, CalibratedMarket)> {
    let tona = PenalizedTona::new()?;
    let inputs = Inputs::new(&[&tona.case], tona.case.as_of, "")?;
    let result = calibrate(
        &tona.definition(gauss_newton)?,
        &tona.quotes()?,
        &inputs.context(),
    )?;
    Ok((tona, inputs, result))
}

/// 曲線の θ を置き換えた市場で、OIS のパーレートを求める。
fn par_rates(case: &Case, market: &Market<f64>, ctx: &Context<'_>) -> Result<Vec<f64>> {
    case.nodes
        .iter()
        .map(|node| {
            let value = |rate| {
                let plan = prepare_ois(case, node, rate, &COLLATERAL, market, ctx)?;
                Ok(evaluate(&plan, market)?.get(case.currency))
            };
            let zero = value(0.0)?;
            let one: f64 = value(1.0)?;
            Ok(-zero / (one - zero))
        })
        .collect()
}

#[test]
fn penalized_meeting_steps_fit_standard_tenor_ois() -> Result<()> {
    // R4: 1 年までの 4 回の会合（5 区間）を 1M、6M、1Y の OIS だけで当てる。未知数より quote が
    // 少ないので、隣り合う区間のフォワードの差の罰則 (λ/2)Σ((f_{i+1} − f_i)/s)² で一意にする。
    // 罰則を強くすると quote を再現しきれず、残差の残る当てはめになる。
    let (tona, inputs, result) = penalized(false)?;
    let ctx = inputs.context();
    let blocks: Vec<_> = result
        .groups()
        .iter()
        .map(|group| group.blocks().to_vec())
        .collect();
    // 年末の跳び（直接マーク）を上流に、階段、長期の節点の順に解く。
    assert_eq!(
        blocks,
        [
            vec![tona.turn_block()?],
            vec![tona.meetings()?],
            vec![tona.long()?]
        ]
    );
    let [turn, meetings, long] = [0, 1, 2].map(|index| solution(&result, index));
    let [turn, meetings, long] = [turn?, meetings?, long?];
    assert_eq!(turn.hessian(), None);
    assert_eq!(meetings.hessian(), Some(HessianMode::Exact));
    assert_eq!(long.hessian(), Some(HessianMode::Exact));
    assert!(
        meetings
            .residuals()
            .iter()
            .any(|residual| residual.abs() > 1e-6),
        "{:?}",
        meetings.residuals()
    );

    // 校正と取引の一致: OIS を取引として評価したパーは、quote に解の残差を足した値。
    let market = result.market();
    let pars = par_rates(&tona.case, market, &ctx)?;
    let residuals: Vec<f64> = meetings
        .residuals()
        .iter()
        .chain(long.residuals())
        .copied()
        .collect();
    for ((par, node), residual) in pars.iter().zip(&tona.case.nodes).zip(&residuals) {
        assert_close!(*par, node.rate + residual, 1e-10, 1e-12);
    }

    // 停留条件を独立に検算する。∂F_i/∂f_j を θ の中心差分で求め、
    // G_j = Σ_i (F_i/σ²)·∂F_i/∂f_j + (λ/s²)·((f_j − f_{j−1}) − (f_{j+1} − f_j)) が 0 になる。
    let range = market
        .structure()
        .layout()
        .range(&tona.meetings()?)
        .ok_or_else(|| Error::Internal("missing meetings block".into()))?;
    let forwards = &market.theta()[range.clone()];
    let h = 1e-7;
    let penalty = PENALTY_WEIGHT / (PENALTY_SCALE * PENALTY_SCALE);
    let last = forwards.len() - 1;
    for (j, index) in range.enumerate() {
        let bumped = |shift: f64| {
            let mut theta = market.theta().to_vec();
            theta[index] += shift;
            par_rates(
                &tona.case,
                &Market::from_theta(Arc::clone(market.structure()), theta)?,
                &ctx,
            )
        };
        let [up, down] = [bumped(h)?, bumped(-h)?];
        let data: Vec<f64> = (0..3)
            .map(|i| residuals[i] / (OIS_TOLERANCE * OIS_TOLERANCE) * (up[i] - down[i]) / (2.0 * h))
            .collect();
        let left = if j > 0 {
            forwards[j] - forwards[j - 1]
        } else {
            0.0
        };
        let right = if j < last {
            forwards[j + 1] - forwards[j]
        } else {
            0.0
        };
        let gradient = data.iter().sum::<f64>() + penalty * (left - right);
        let scale =
            data.iter().map(|term| term.abs()).sum::<f64>() + penalty * (left.abs() + right.abs());
        assert!(
            gradient.abs() <= 1e-6 * scale,
            "G[{j}] = {gradient:e} of {scale:e}"
        );
    }
    Ok(())
}

#[test]
fn gauss_newton_records_are_flagged_and_differ_from_exact() -> Result<()> {
    // Gauss–Newton の記録は残差の 2 階の項を省くので、残差の残る階段のグループでは ∂G/∂θ が
    // 厳密な記録と違う。解（θ）は同じで、近似の記録を解と警告の両方で報告する。
    let (_, _, exact) = penalized(false)?;
    let (tona, _, approximate) = penalized(true)?;
    for (lhs, rhs) in exact
        .market()
        .theta()
        .iter()
        .zip(approximate.market().theta())
    {
        assert_close!(*lhs, *rhs, 1e-10, 1e-12);
    }
    let gauss_newton = solution(&approximate, 1)?;
    assert_eq!(gauss_newton.hessian(), Some(HessianMode::GaussNewton));
    let flagged: Vec<QuoteId> = approximate
        .warnings()
        .iter()
        .filter_map(|warning| match warning {
            Warning::CalibrationGaussNewton { quote } => Some(quote.clone()),
            _ => None,
        })
        .collect();
    let expected: Vec<QuoteId> = tona
        .case
        .nodes
        .iter()
        .map(|node| tona.case.quote(node))
        .collect::<Result<_>>()?;
    assert_eq!(flagged.len(), expected.len());
    assert!(expected.iter().all(|quote| flagged.contains(quote)));
    // 同じ右辺の逆作用素の値が、2 階の項の分だけ違う。
    let rhs = vec![1.0; 5];
    let [exact_inverse, approximate_inverse] = [solution(&exact, 1)?, gauss_newton].map(|s| {
        s.dg_dtheta()
            .solve(rhs.as_slice())
            .map(|vector| vector.as_slice().to_vec())
            .map_err(linear)
    });
    let difference = exact_inverse?
        .iter()
        .zip(approximate_inverse?)
        .map(|(lhs, rhs)| (lhs - rhs).abs() / lhs.abs().max(1e-300))
        .fold(0.0_f64, f64::max);
    assert!(difference > 1e-6, "{difference:e}");
    Ok(())
}
`````

`````
crates/rhodium-pricing/tests/linear_pricing.rs
//! 標準商品と検証済み IR を、本番の準備・結合・評価の入口で照合する。

use std::{num::NonZeroU8, sync::Arc};

use rhodium_core::{
    Error, Result, Warning,
    calendar::{
        BaseCalendar, BusinessDayConvention, CalendarSpec, EndOfMonthRule,
        Extrapolation as CalendarExtrapolation, WeekendRule,
    },
    currency::{Asset, Currency},
    date::{Date, DateInterval, Weekday},
    day_count::{AccrualContext, DayCount},
    fixing::{FixingStore, FixingStoreParts, ObservationPoint},
    id::IndexId,
    money::Money,
    observation::{
        Aggregation, CompoundingForm, ObservationDays, ObservationMethod, ObservationRounding,
        ObservationRule, SpreadTreatment, observation_schedule,
    },
    quantity::{Rate, Simple, Spread},
    refdata::{
        CurrencyDefinition, DailyPublication, Deliverability, FallbackDefinition,
        FallbackReplacement, FallbackTrigger, IndexDefinition, IndexKind, IndexStatus,
        PublicationRule, RefData, RefDataParts,
    },
    rounding::{Rounding, RoundingMode},
    schedule::{
        AccrualDates, EndRule, ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods,
        GeneratedPeriods, PaymentGrouping, PaymentRule, PeriodAnchor, PeriodCombination,
        ScheduleSpec, SettlementContext, StartRule, StubLength, StubSpec, generate,
    },
    tenor::{DateOffset, Frequency, TenorUnit},
    time::{Instant, LocalEventTime, UtcOffset, ZoneRules},
};
use rhodium_market::{
    curve::{CurveRole, CurveSpec, Extrapolation},
    market::Market,
    structure::{BlockId, MarketStructure},
};
use rhodium_pricing::{
    config::{Method, PricingConfig},
    context::Context,
    prepare::{PrepareInput, Prepared, evaluate, prepare},
    resolve::CollateralKey,
};
use rhodium_products::{
    ir::{Contract, Product},
    standard::{Coupon, Leg, NotionalExchange, PayReceive, PaymentRounding, StandardProduct, Swap},
    trade::{FallbackSelection, FallbackTerms, ProductSpec, Trade},
};
use rhodium_testkit::assert_close;
use scioxide_ad::Tape;
use scioxide_core::{ExecutionContext, InvalidInputReason, SciOxideError, argument};

const INDEX: &str = "JPY-TONA";
const REL: f64 = 1e-12;
const JPY_ABS: f64 = 1e-6;
const NOTIONAL: f64 = 1_000_000.0;
const DISABLE_DAILY_FORWARD_PRODUCT: &str = r#"disabled = ["daily_forward_product"]"#;

fn at(date: Date) -> Result<Instant> {
    format!("{date}T00:00:00Z").parse()
}

fn calendar() -> Result<CalendarSpec> {
    Ok(CalendarSpec::Base("PRICING".parse()?))
}

fn refdata(ccy: Currency, index: &str, basis: DayCount, holidays: Vec<Date>) -> Result<RefData> {
    let utc = ZoneRules::new(
        "2024-10-01T00:00:00Z".parse()?..="2025-05-01T00:00:00Z".parse()?,
        UtcOffset::from_seconds(0)?,
        Vec::new(),
    )?;
    refdata_in_zone(ccy, index, basis, holidays, "Etc/UTC", utc)
}

/// 観測、公表、決済の時刻をすべて指定の zone の現地 00:00 にした参照データ。
fn refdata_in_zone(
    ccy: Currency,
    index: &str,
    basis: DayCount,
    holidays: Vec<Date>,
    zone: &str,
    rules: ZoneRules,
) -> Result<RefData> {
    RefData::new(refdata_parts(ccy, index, basis, holidays, zone, rules)?)
}

fn refdata_parts(
    ccy: Currency,
    index: &str,
    basis: DayCount,
    holidays: Vec<Date>,
    zone: &str,
    rules: ZoneRules,
) -> Result<RefDataParts> {
    let mut parts = RefDataParts::default();
    let calendar_id = "PRICING".parse()?;
    let clock = LocalEventTime {
        time: "00:00:00".parse()?,
        zone: zone.parse()?,
    };
    parts.calendars.insert(
        calendar_id,
        BaseCalendar::new(
            "2024-10-01".parse()?..="2025-04-30".parse()?,
            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
            holidays,
            CalendarExtrapolation::Strict,
        )?,
    );
    parts.zones.insert(clock.zone.clone(), rules);
    let minor_units = if ccy == Currency::JPY { 0 } else { 2 };
    parts.currencies.insert(
        ccy,
        CurrencyDefinition {
            code: ccy,
            minor_units,
            payment_rounding: Rounding::new(
                0.1_f64.powi(i32::from(minor_units)),
                RoundingMode::HalfEven,
            )?,
            deliverability: Deliverability::Deliverable,
            settlement_time: clock.clone(),
            settlement_calendar: "PRICING".parse()?,
        },
    );
    let index: IndexId = index.parse()?;
    parts.indices.insert(
        index.clone(),
        IndexDefinition {
            id: index,
            kind: IndexKind::OvernightRate {
                ccy,
                day_count: basis,
            },
            calendar: calendar()?,
            publication: PublicationRule::Daily(DailyPublication {
                business_days: 1,
                calendar: calendar()?,
                time: clock.time,
                zone: clock.zone.clone(),
            }),
            observation_time: Some(clock),
            status: IndexStatus::Active,
        },
    );
    Ok(parts)
}

fn config(
    ccy: Currency,
    index: &str,
    single_curve: bool,
    settlement: &str,
    conventions: &str,
) -> Result<PricingConfig> {
    let projection = if single_curve {
        "discount"
    } else {
        "projection"
    };
    toml::from_str(&format!(
        r#"
version = 1
reporting_currency = "{ccy}"
settlement = "{settlement}"

[conventions]
{conventions}

[[discounting]]
currency = "{ccy}"
collateral = "uncollateralized"
curve = "discount"

[forecasting]
"{index}" = "{projection}"
"#
    ))
    .map_err(|error| Error::InvalidData(error.to_string()))
}

fn structure(
    ccy: Currency,
    as_of: Instant,
    pillars: &[Instant],
    single_curve: bool,
) -> Result<Arc<MarketStructure>> {
    let discount_roles = if single_curve {
        vec![CurveRole::Discount, CurveRole::Projection]
    } else {
        vec![CurveRole::Discount]
    };
    let mut curves = vec![CurveSpec::new(
        "discount".parse()?,
        ccy,
        discount_roles,
        pillars.iter().copied(),
        Extrapolation::Reject,
    )?];
    if !single_curve {
        curves.push(CurveSpec::new(
            "projection".parse()?,
            ccy,
            [CurveRole::Projection],
            pillars.iter().copied(),
            Extrapolation::Reject,
        )?);
    }
    Ok(Arc::new(MarketStructure::new(as_of, curves)?))
}

fn schedule(dates: &[Date], pays: &[Date]) -> Result<ScheduleSpec> {
    let periods = dates
        .windows(2)
        .zip(pays)
        .map(|(dates, &payment)| {
            Ok(ExplicitPeriod {
                accrual: DateInterval::new(dates[0], dates[1])?,
                payment,
                reset: None,
                settlement: None,
                references: Vec::new(),
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

fn rule(
    method: ObservationMethod,
    form: CompoundingForm,
    spread: SpreadTreatment,
    rounding: ObservationRounding,
) -> Result<ObservationRule> {
    ObservationRule::new(
        ObservationDays::BusinessDays(calendar()?),
        method,
        Aggregation::Compound { form, spread },
        rounding,
    )
}

fn plain() -> Result<ObservationRule> {
    rule(
        ObservationMethod::Plain,
        CompoundingForm::SimpleFactor,
        SpreadTreatment::Exclusive,
        ObservationRounding::default(),
    )
}

fn swap(
    money: Money,
    schedule: &ScheduleSpec,
    basis: &DayCount,
    index: &str,
    fixed_rate: f64,
    float_spread: f64,
    rule: ObservationRule,
) -> Result<StandardProduct> {
    let leg = |side, coupon| {
        Leg::new(
            side,
            money,
            schedule.clone(),
            basis.clone(),
            coupon,
            NotionalExchange::None,
            PaymentRounding::None,
        )
    };
    Ok(StandardProduct::Swap(Swap::from_legs([
        leg(
            PayReceive::Pay,
            Coupon::Fixed {
                rate: Rate::new(fixed_rate, basis.clone(), Simple)?,
            },
        )?,
        leg(
            PayReceive::Receive,
            Coupon::Compounded {
                index: index.parse()?,
                rule,
                spread: Spread::new(float_spread)?,
            },
        )?,
    ])?))
}

fn prepare_standard(
    product: &StandardProduct,
    date: Date,
    structure: &Arc<MarketStructure>,
    context: &Context<'_>,
) -> Result<Prepared> {
    prepare(
        PrepareInput::Standard {
            product,
            trade_date: date,
        },
        &CollateralKey::Uncollateralized,
        &[],
        structure,
        context,
    )
}

#[derive(serde::Deserialize)]
struct QuantLibSwap {
    currency: Currency,
    index: String,
    as_of: Date,
    dates: Vec<Date>,
    payments: Vec<Date>,
    pillars: Vec<Date>,
    holidays: Vec<Date>,
    projection_dfs: Vec<f64>,
    discount_dfs: Vec<f64>,
    single_curve: bool,
    payment_lag: u32,
    fixed_rate: f64,
    notional: f64,
    pv: f64,
    fixed_pv: f64,
    floating_pv: f64,
}

fn compare_quantlib(currency: Currency) -> Result<()> {
    let fixtures: Vec<QuantLibSwap> =
        serde_json::from_str(include_str!("fixtures/quantlib/swaps.json"))
            .map_err(|error| Error::InvalidData(error.to_string()))?;
    let fixture = fixtures
        .into_iter()
        .find(|case| case.currency == currency)
        .ok_or_else(|| Error::InvalidData(format!("missing fixture for {currency}")))?;
    let basis = if currency == Currency::JPY {
        DayCount::Act365F
    } else {
        DayCount::Act360
    };
    let as_of = at(fixture.as_of)?;
    let pillars = fixture
        .pillars
        .into_iter()
        .map(at)
        .collect::<Result<Vec<_>>>()?;
    let structure = structure(currency, as_of, &pillars, fixture.single_curve)?;
    let refdata = refdata(currency, &fixture.index, basis.clone(), fixture.holidays)?;
    let config = config(
        currency,
        &fixture.index,
        fixture.single_curve,
        "exclude",
        "",
    )?;
    let context = Context {
        as_of,
        refdata: &refdata,
        config: &config,
        fixings: &FixingStore::default(),
    };
    let schedule = schedule(&fixture.dates, &fixture.payments)?;
    let notional = Money::new(currency, fixture.notional);
    let product = swap(
        notional,
        &schedule,
        &basis,
        &fixture.index,
        fixture.fixed_rate,
        0.0,
        plain()?,
    )?;
    let prepared = prepare_standard(&product, fixture.as_of, &structure, &context)?;
    let mut theta: Vec<_> = fixture.discount_dfs.into_iter().map(f64::ln).collect();
    if !fixture.single_curve {
        theta.extend(fixture.projection_dfs.into_iter().map(f64::ln));
    }

    // 両脚の絶対値の和を、差し引いた PV の相対誤差の尺度に使う。
    let minor_units = refdata.currency(currency)?.minor_units;
    let abs = 0.1_f64.powi(i32::from(minor_units)) * 1e-6;
    let scale = fixture.fixed_pv.abs() + fixture.floating_pv.abs();
    let market = Market::from_theta(structure.clone(), theta)?;
    let pv = evaluate(&prepared, &market)?.get(currency);
    assert_close!(pv, fixture.pv, REL, abs.max(REL * scale));

    // 固定金利が 0 の同じ契約で変動脚の PV を照合し、差から固定脚も確かめる。
    let float_only = swap(
        notional,
        &schedule,
        &basis,
        &fixture.index,
        0.0,
        0.0,
        plain()?,
    )?;
    let float_only = prepare_standard(&float_only, fixture.as_of, &structure, &context)?;
    let floating = evaluate(&float_only, &market)?.get(currency);
    assert_close!(floating, fixture.floating_pv, REL, abs);
    assert_close!(pv - floating, fixture.fixed_pv, REL, abs);
    assert_eq!(prepared.method(), Method::Linear);
    assert_eq!(
        prepared.blocks().len(),
        if fixture.single_curve { 1 } else { 2 }
    );
    assert_eq!(
        prepared.warnings().contains(&Warning::PaymentDelayIgnored),
        fixture.payment_lag > 0
    );
    assert_eq!(
        prepared
            .warnings()
            .contains(&Warning::CrossCollateralProjection),
        !fixture.single_curve
    );
    assert!(!prepared.warnings().contains(&Warning::DailyForwardProduct));
    Ok(())
}

#[test]
fn tona_ois_pv_matches_quantlib() -> Result<()> {
    // 同じ TONA 曲線で割引と予想を行い、支払ラグなしの OIS を照合する。
    compare_quantlib(Currency::JPY)
}

#[test]
fn sofr_fixed_rfr_swap_pv_matches_quantlib() -> Result<()> {
    // 予想曲線と割引曲線を分け、支払日は観測期間の末尾から 2 営業日後に置く。
    compare_quantlib(Currency::USD)
}

fn jpy() -> Result<RefData> {
    refdata(Currency::JPY, INDEX, DayCount::Act365F, Vec::new())
}

fn fixings(
    data: &RefData,
    published: &[(&str, f64)],
    estimates: &[(&str, f64)],
) -> Result<FixingStore> {
    let index: IndexId = INDEX.parse()?;
    let rows = |items: &[(&str, f64)]| {
        items
            .iter()
            .map(|&(day, rate)| Ok((index.clone(), ObservationPoint::Date(day.parse()?), rate)))
            .collect::<Result<Vec<_>>>()
    };
    FixingStore::new(
        data,
        FixingStoreParts {
            published: rows(published)?,
            estimates: rows(estimates)?,
            ..FixingStoreParts::default()
        },
    )
}

fn jpy_swap(
    dates: &[&str],
    pays: &[&str],
    fixed: f64,
    spread: f64,
    rule: ObservationRule,
) -> Result<StandardProduct> {
    let dates = dates
        .iter()
        .map(|date| date.parse())
        .collect::<Result<Vec<_>>>()?;
    let pays = pays
        .iter()
        .map(|date| date.parse())
        .collect::<Result<Vec<_>>>()?;
    swap(
        Money::new(Currency::JPY, NOTIONAL),
        &schedule(&dates, &pays)?,
        &DayCount::Act365F,
        INDEX,
        fixed,
        spread,
        rule,
    )
}

fn linear_df(as_of: Instant, pillar: Instant, theta: f64, at: Instant) -> f64 {
    let elapsed = at.unix_micros() - as_of.unix_micros();
    let total = pillar.unix_micros() - as_of.unix_micros();
    (theta * elapsed as f64 / total as f64).exp()
}

fn jpy_pv_and_gradient(
    structure: &Arc<MarketStructure>,
    prepared: &Prepared,
    theta: f64,
) -> std::result::Result<(f64, f64), SciOxideError> {
    let mut gradient = [0.0];
    let pv = Tape::default().gradient_into(
        |values| {
            Market::from_theta(structure.clone(), values.to_vec())
                .and_then(|market| evaluate(prepared, &market))
                .map(|pv| pv.get(Currency::JPY))
                .map_err(|error| {
                    SciOxideError::invalid_input(
                        argument!("market"),
                        InvalidInputReason::Inconsistent,
                        error.to_string(),
                    )
                })
        },
        &[theta],
        &mut gradient,
        ExecutionContext::sequential(),
    )?;
    Ok((pv, gradient[0]))
}

#[test]
fn fixed_prefix_and_intraday_forward_use_original_coupon_weights()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // 木曜分は確定済み、金曜の観測の途中で評価する。月曜分以降だけ端点の比に畳む。
    let as_of = "2025-01-10T12:00:00Z".parse()?;
    let pillar = "2025-02-10T00:00:00Z".parse()?;
    let data = jpy()?;
    let fixings = fixings(&data, &[("2025-01-09", 0.04)], &[])?;
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
    let structure = structure(Currency::JPY, as_of, &[pillar], true)?;
    let product = jpy_swap(
        &["2025-01-09", "2025-01-14"],
        &["2025-01-14"],
        0.0,
        0.015,
        plain()?,
    )?;
    let prepared = prepare_standard(&product, "2025-01-08".parse()?, &structure, &context)?;

    // 週末をまたぐ金曜の率は残り 2.5/365 を分母とし、元の 3/365 の重みに掛ける。
    // 月曜・火曜の DF は、評価時点から補間のノードまでの 30.5 日に対する経過日数で求める。
    let monday_time = 2.5 / 30.5;
    let tuesday_time = 3.5 / 30.5;
    let full_to_remaining = 3.0 / 2.5;
    let fixed_factor = 1.0 + 0.04 / 365.0;
    let spread_accrual = 0.015 * 5.0 / 365.0;
    let expected = |theta: f64| {
        let monday = (theta * monday_time).exp();
        let tuesday = (theta * tuesday_time).exp();
        let friday_factor = 1.0 + full_to_remaining * (monday.recip() - 1.0);
        let growth = fixed_factor * friday_factor * monday / tuesday;
        NOTIONAL * (growth - 1.0 + spread_accrual) * tuesday
    };
    let theta = -0.0018;
    let market = Market::from_theta(structure.clone(), vec![theta])?;
    let pv = evaluate(&prepared, &market)?.get(Currency::JPY);
    assert_close!(pv, expected(theta), REL, JPY_ABS);

    // 同じ準備で θ を変えても、金曜の因子と月曜以降の比の両方が価格と感応度に反映される。
    let theta = -0.0042;
    let (pv, gradient) = jpy_pv_and_gradient(&structure, &prepared, theta)?;
    assert_close!(pv, expected(theta), REL, JPY_ABS);
    // 金曜の因子と月曜の DF の積は q + (1 − q)·DF(月曜)。q は元の重みと残存期間の比。
    let derivative = NOTIONAL
        * (fixed_factor * (1.0 - full_to_remaining) * monday_time * (theta * monday_time).exp()
            + (spread_accrual - 1.0) * tuesday_time * (theta * tuesday_time).exp());
    assert_close!(gradient, derivative, 1e-10, JPY_ABS);
    assert!(
        prepared
            .warnings()
            .iter()
            .any(|warning| matches!(warning, Warning::UnpublishedFixingForecast { .. }))
    );
    Ok(())
}

#[test]
fn partially_fixed_excess_keeps_the_published_digits() -> Result<()> {
    // 公表済みの 1 日分（0.5%）と予想の残り（θ = 0 なので予想の率は 0）からなる 4 日の期間で、
    // 元本 73,109,500 円の利息はちょうど 1,001.5 円、円単位の HalfEven では 1,002 円になる。
    // 簡約が確定因子を積 1 + 0.005/365 に丸めると、予想の側で exp_m1 を使っても超過の桁が戻らず
    // 1,001.4999999932 円（1,001 円）になる。端点の比に畳む経路（Plain）と日次積の経路（lookback）の
    // どちらでも、確定因子の超過 e が評価まで届き、予想の超過 f と e + f + e·f で合わさることを確かめる
    // （f = exp_m1(0) = 0 なので、利息は e = r·w から求めた値そのものになる）。
    let data = jpy()?;
    let pillar = "2025-02-10T00:00:00Z".parse()?;
    let notional = Money::new(Currency::JPY, 73_109_500.0);
    let lookback = rule(
        ObservationMethod::Lookback(NonZeroU8::MIN),
        CompoundingForm::SimpleFactor,
        SpreadTreatment::Exclusive,
        ObservationRounding::default(),
    )?;
    // Plain は木曜の公表後（金曜 0 時）に、lookback は水曜の公表後で木曜の公表前（木曜の昼）に評価する。
    // どちらも確定した最初の要素の率が 0.5%、残りの要素は予想になる。
    for (rule, as_of, published) in [
        (plain()?, "2025-01-10T00:00:00Z", "2025-01-09"),
        (lookback, "2025-01-09T12:00:00Z", "2025-01-08"),
    ] {
        let as_of: Instant = as_of.parse()?;
        let fixings = fixings(&data, &[(published, 0.005)], &[])?;
        let config = config(Currency::JPY, INDEX, true, "exclude", "")?;
        let context = Context {
            as_of,
            refdata: &data,
            fixings: &fixings,
            config: &config,
        };
        let schedule = schedule(
            &["2025-01-09".parse()?, "2025-01-13".parse()?],
            &["2025-01-13".parse()?],
        )?;
        let leg = |side, coupon| {
            Leg::new(
                side,
                notional,
                schedule.clone(),
                DayCount::Act365F,
                coupon,
                NotionalExchange::None,
                PaymentRounding::CurrencyConvention,
            )
        };
        let product = StandardProduct::Swap(Swap::from_legs([
            leg(
                PayReceive::Pay,
                Coupon::Fixed {
                    rate: Rate::new(0.0, DayCount::Act365F, Simple)?,
                },
            )?,
            leg(
                PayReceive::Receive,
                Coupon::Compounded {
                    index: INDEX.parse()?,
                    rule,
                    spread: Spread::new(0.0)?,
                },
            )?,
        ])?);
        let structure = structure(Currency::JPY, as_of, &[pillar], true)?;
        let prepared = prepare_standard(&product, "2025-01-07".parse()?, &structure, &context)?;
        // θ = 0 では割引因子が 1 なので、PV は丸めた受渡額そのもの。
        let pv =
            evaluate(&prepared, &Market::from_theta(structure, vec![0.0])?)?.get(Currency::JPY);
        assert_close!(pv, 1002.0, REL, JPY_ABS);
    }
    Ok(())
}

#[test]
fn growth_and_excess_of_one_compounding_share_its_plan() -> Result<()> {
    // 同じ Compounded C を `C − 1` と別の親（M·C）の両方が読むと、超過の命令と成長率の命令は
    // 一つの計画を共有し、どちらも確定部分の超過 e を使う（超過だけが読む場合は成長率の命令を作らない）。
    // 1 日分が公表済み（0.5%）、残りが予想（θ = 0 なので率 0）の 4 日の期間で、
    // N·(C − 1) + M·C の PV は N·e + M·(1 + e)（e = 0.005/365、割引因子は 1）になる。
    let data = jpy()?;
    let as_of: Instant = "2025-01-10T00:00:00Z".parse()?;
    let fixings = fixings(&data, &[("2025-01-09", 0.005)], &[])?;
    let config = config(Currency::JPY, INDEX, true, "exclude", "")?;
    let context = Context {
        as_of,
        refdata: &data,
        fixings: &fixings,
        config: &config,
    };
    let currency = data.currency(Currency::JPY)?;
    let settlement_calendar = CalendarSpec::Base(currency.settlement_calendar.clone());
    let spec = schedule(
        &["2025-01-09".parse()?, "2025-01-13".parse()?],
        &["2025-01-13".parse()?],
    )?;
    let (generated, _) = generate(
        &spec,
        "2025-01-07".parse()?,
        data.calendars(),
        &SettlementContext {
            zones: data.zones(),
            default_time: Some(&currency.settlement_time),
            calendar: Some(&settlement_calendar),
        },
    )?;
    let accrual = generated
        .accrual_context(0, None)
        .ok_or_else(|| Error::Internal("missing accrual context".into()))?;
    let rule = plain()?;
    let (observations, _) = observation_schedule(
        &generated.periods()[0],
        &rule,
        &DayCount::Act365F,
        &accrual,
        data.calendars(),
    )?;
    let index: IndexId = INDEX.parse()?;
    let (n, m) = (100_000_000.0, 7.0);
    let product = Product::new(
        Contract::build(|b| {
            let growth = b.compounded(&index, &observations, &rule, Spread::new(0.0)?)?;
            let one = b.constant(1.0);
            let excess = b.sub(growth, one);
            let notional = b.constant(n);
            let scale = b.constant(m);
            let interest = b.mul([notional, excess]);
            let level = b.mul([scale, growth]);
            let amount = b.add([interest, level]);
            Ok(b.pay(
                Asset::Cash(Currency::JPY),
                amount,
                generated.payments()[0].at(),
            ))
        })?,
        &data,
    )?;
    let market_structure = structure(
        Currency::JPY,
        as_of,
        &["2025-02-10T00:00:00Z".parse()?],
        true,
    )?;
    let prepared = prepare(
        PrepareInput::Ir(&product),
        &CollateralKey::Uncollateralized,
        &[],
        &market_structure,
        &context,
    )?;
    let pv =
        evaluate(&prepared, &Market::from_theta(market_structure, vec![0.0])?)?.get(Currency::JPY);
    let excess = 0.005 / 365.0;
    assert_close!(pv, n * excess + m * (1.0 + excess), REL, JPY_ABS);
    Ok(())
}

#[test]
fn explicit_long_stub_matches_generated_rfr_coupon() -> Result<()> {
    // 月払いの 2024-12-02→2025-01-15 は、規則から生成すると長い後ろのスタブ一つ（通常期間
    // 12-02→01-02 と仮想の 01-02→02-02 を参照枠に持つ）で、TONA を 45 日通して複利にする。
    // 確認書の明示の期間でも、一つの期間に二つの参照期間を並べれば同じクーポンで、PV が一致する。
    // 二つの期間に分けて同じ支払日に Sum で足すと、複利が 01-02 で切れて
    // (Π₁ − 1) + (Π₂ − 1) になり、交差項 (Π₁ − 1)(Π₂ − 1) の分だけ別の額になる（元本 1 億円で数百円）。
    let as_of: Instant = "2024-11-29T00:00:00Z".parse()?;
    let pillar = "2025-02-10T00:00:00Z".parse()?;
    let data = jpy()?;
    let fixings = fixings(&data, &[], &[])?;
    let config = config(Currency::JPY, INDEX, true, "exclude", "")?;
    let context = Context {
        as_of,
        refdata: &data,
        fixings: &fixings,
        config: &config,
    };
    let structure = structure(Currency::JPY, as_of, &[pillar], true)?;
    // 約 5% の翌日物金利になる傾き。
    let market = Market::from_theta(structure.clone(), vec![-0.01])?;
    let monthly: Frequency = "1M".parse()?;
    let (start, end) = ("2024-12-02".parse()?, "2025-01-15".parse()?);
    let generated = ScheduleSpec::Generated {
        periods: GeneratedPeriods::new(
            StartRule::Date(start),
            EndRule::Date(end),
            monthly,
            None,
            StubSpec::Final(StubLength::Long),
            calendar()?,
            BusinessDayConvention::ModifiedFollowing,
            AccrualDates::Adjusted,
        )?,
        payment: PaymentRule {
            anchor: PeriodAnchor::PeriodEnd,
            offset: DateOffset::new(0, TenorUnit::Day),
            calendar: calendar()?,
            bdc: BusinessDayConvention::Following,
            grouping: PaymentGrouping::EachPeriod,
            settlement: None,
        },
    };
    // どの行も満期の 2025-01-15 に払う。
    let end_payment = end;
    let period = |start: &str, end: &str, references: Vec<DateInterval>| -> Result<_> {
        Ok(ExplicitPeriod {
            accrual: DateInterval::new(start.parse()?, end.parse()?)?,
            payment: end_payment,
            reset: None,
            settlement: None,
            references,
        })
    };
    let regular = DateInterval::new(start, "2025-01-02".parse()?)?;
    let virtual_regular = DateInterval::new("2025-01-02".parse()?, "2025-02-02".parse()?)?;
    let explicit = |periods| -> Result<ScheduleSpec> {
        Ok(ScheduleSpec::Explicit(ExplicitPeriods::new(
            periods,
            ExplicitDatePolicy::KeepAsConfirmed,
            PeriodCombination::Sum,
