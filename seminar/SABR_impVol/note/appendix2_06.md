    /// 同じ quote を読む binding の種類が異なる場合。
    pub fn new(
        valuation_calendar: CalendarId,
        valuation_zone: TzId,
        curves: impl IntoIterator<Item = CurveSpec>,
        bindings: impl IntoIterator<Item = InstrumentBinding>,
    ) -> Result<Self> {
        let mut assignment = BTreeMap::new();
        let mut curves: Vec<_> = curves.into_iter().collect();
        curves.sort_by(|lhs, rhs| lhs.id().cmp(rhs.id()));
        if let Some(pair) = curves.windows(2).find(|pair| pair[0].id() == pair[1].id()) {
            return Err(Error::InvalidInput(format!(
                "duplicate curve {}",
                pair[0].id()
            )));
        }
        assignment.extend(
            curves
                .iter()
                .flat_map(CurveSpec::parameter_blocks)
                .map(|(block, count)| (block, vec![None; count])),
        );
        let bindings: Vec<_> = bindings.into_iter().collect();
        let mut quote_kinds = BTreeMap::new();
        for binding in &bindings {
            let kind = binding.kind();
            if quote_kinds
                .insert(binding.quote.clone(), kind)
                .is_some_and(|previous| previous != kind)
            {
                return Err(Error::InvalidInput(format!(
                    "inconsistent kinds for quote {}",
                    binding.quote
                )));
            }
            let block = binding.determines.block();
            let assigned = assignment.get_mut(block).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "quote {} targets unknown curve {}",
                    binding.quote,
                    block.curve()
                ))
            })?;
            // 直接マークは要素の値を決める。許容幅を持つ直接マークどうしだけは同じ要素を重ねて当て、
            // 最小二乗で一つの値に決める（手計算の検算 min (θ − q₁)² + 4(θ − q₂)²）。
            let current = if matches!(binding.target, CalibrationTarget::DirectMark) {
                Assigned::Direct {
                    fitted: binding.tolerance.is_some(),
                }
            } else {
                Assigned::Fit
            };
            let selected = match &binding.determines {
                ElementSelector::Block(_) => assigned.as_mut_slice(),
                ElementSelector::Param(param) => {
                    let index = usize::try_from(param.element())
                        .ok()
                        .and_then(|element| assigned.get_mut(element))
                        .ok_or_else(|| {
                            Error::InvalidInput(format!(
                                "quote {} targets unknown element {} of curve {}",
                                binding.quote,
                                param.element(),
                                block.curve()
                            ))
                        })?;
                    std::slice::from_mut(index)
                }
            };
            for status in selected {
                match (*status, current) {
                    (None, _) => *status = Some(current),
                    (Some(Assigned::Fit), Assigned::Fit)
                    | (
                        Some(Assigned::Direct { fitted: true }),
                        Assigned::Direct { fitted: true },
                    ) => {}
                    (Some(Assigned::Fit | Assigned::Direct { .. }), _) => {
                        return Err(Error::InvalidInput(format!(
                            "conflicting direct mark for {} {}",
                            block.curve(),
                            block.part()
                        )));
                    }
                }
            }
        }
        check_assigned(&assignment)?;
        Ok(Self {
            valuation_calendar,
            valuation_zone,
            curves,
            bindings,
            fx: Vec::new(),
            lgm: BTreeMap::new(),
            penalties: BTreeMap::new(),
            gauss_newton: BTreeSet::new(),
            quote_kinds,
        })
    }

    /// 会合日の階段のブロックに罰則を掛け、そのブロックを最小二乗で当てる。
    ///
    /// # Errors
    ///
    /// ブロックが市場の定義の曲線の `Steps` の部品でないか、区間が二つ未満で罰則の行が無いか、
    /// 同じブロックに罰則を二度掛けた場合は `InvalidInput`。
    pub fn with_step_penalties(
        mut self,
        penalties: impl IntoIterator<Item = StepPenalty>,
    ) -> Result<Self> {
        for penalty in penalties {
            let steps = self
                .curves
                .iter()
                .find(|curve| curve.id() == penalty.block.curve())
                .and_then(|curve| {
                    curve
                        .parameter_blocks()
                        .find(|(block, _)| *block == penalty.block)
                        .filter(|_| curve.step_blocks().any(|block| block == penalty.block))
                });
            if steps.is_none_or(|(_, count)| count < 2) {
                return Err(Error::InvalidInput(format!(
                    "step penalty needs a steps part with two or more intervals: {} {}",
                    penalty.block.curve(),
                    penalty.block.part()
                )));
            }
            if self
                .penalties
                .insert(penalty.block.clone(), penalty)
                .is_some()
            {
                return Err(Error::InvalidInput("duplicate step penalty".into()));
            }
        }
        Ok(self)
    }

    /// 最小二乗のグループの解の記録を Gauss–Newton（残差の 2 階の項を省く）で作るブロック。
    ///
    /// 既定は 2 階の項を含む厳密な記録で、このときだけ 1 次の quote 感応度が厳密になる。
    /// 要素の多い当てはめで記録の費用を下げるために選べ、選んだことを解の記録と警告に残す。
    ///
    /// # Errors
    ///
    /// ブロックが市場の定義の曲線に無ければ `InvalidInput`。
    pub fn with_gauss_newton(mut self, blocks: impl IntoIterator<Item = BlockId>) -> Result<Self> {
        for block in blocks {
            if !self
                .curves
                .iter()
                .flat_map(CurveSpec::parameter_blocks)
                .any(|(known, _)| known == block)
            {
                return Err(Error::InvalidInput(format!(
                    "Gauss-Newton block {} {} is not in the market definition",
                    block.curve(),
                    block.part()
                )));
            }
            self.gauss_newton.insert(block);
        }
        Ok(self)
    }

    /// LGM の区間と quote ID を市場定義に追加する。
    ///
    /// # Errors
    ///
    /// 通貨の重複、区間とマーク数の不一致、または quote の種類の衝突がある場合。
    pub fn with_lgm(
        mut self,
        models: impl IntoIterator<Item = (Currency, LgmDefinition)>,
    ) -> Result<Self> {
        for (ccy, model) in models {
            let count = model.breaks.len() + 1;
            if model.mean_reversion.len() != count
                || model.volatility.len() != count
                || !model.breaks.is_sorted_by(|a, b| a < b)
            {
                return Err(Error::InvalidData(format!("invalid LGM marks for {ccy}")));
            }
            for quote in model.mean_reversion.iter().chain(&model.volatility) {
                if self
                    .quote_kinds
                    .insert(quote.clone(), QuoteKind::Mark)
                    .is_some_and(|kind| kind != QuoteKind::Mark)
                {
                    return Err(Error::InvalidInput(format!(
                        "inconsistent kinds for quote {quote}"
                    )));
                }
            }
            if self.lgm.insert(ccy, model).is_some() {
                return Err(Error::InvalidData(format!("duplicate LGM model for {ccy}")));
            }
        }
        Ok(self)
    }

    /// 通貨ペアの受渡慣行、quotes、担保別曲線を市場定義に追加する。
    ///
    /// # Errors
    ///
    /// 短期の quote がスポットへの連鎖にならない場合、FX ポイントがあるのに外貨担保曲線
    /// （`FxImplied`）が無い場合、外貨担保曲線を当日受渡の締め後も校正できる quote が無い場合は
    /// `InvalidData`。同じ quote ID を異なる種類で使う場合は `InvalidInput`。
    /// 受渡日と曲線の照合は、評価時点と参照データを受け取る校正時に行う。
    pub fn with_fx(mut self, definitions: impl IntoIterator<Item = FxDefinition>) -> Result<Self> {
        for definition in definitions {
            definition.validate()?;
            // 外貨担保曲線の長期の XCCY は、その曲線の basis ブロックを当てる binding になる。
            if let Some(implied) = definition
                .direct
                .as_ref()
                .and_then(|direct| direct.implied.as_ref())
            {
                self.bindings.extend(implied.long_end_bindings()?);
            }
            for (quote, kind) in definition.quote_kinds() {
                if self
                    .quote_kinds
                    .insert(quote.clone(), kind)
                    .is_some_and(|previous| previous != kind)
                {
                    return Err(Error::InvalidInput(format!(
                        "inconsistent kinds for quote {quote}"
                    )));
                }
            }
            self.fx.push(definition);
        }
        Ok(self)
    }

    pub(crate) fn contains_curve(&self, id: &CurveId) -> bool {
        self.curves
            .binary_search_by(|curve| curve.id().cmp(id))
            .is_ok()
            || self
                .fx
                .iter()
                .filter_map(|pair| pair.direct.as_ref()?.implied.as_ref())
                .any(|implied| implied.id == *id)
    }

    /// この市場で binding が読む quote とその種類。
    #[must_use]
    pub const fn quote_kinds(&self) -> &BTreeMap<QuoteId, QuoteKind> {
        &self.quote_kinds
    }
}
`````

`````
crates/rhodium-pricing/src/definition/fx.rs
//! 市場の FX quote と受渡の定義。

use rhodium_core::{
    Error, Result, Warning,
    calendar::{AdvanceRule, BusinessDayConvention, Calendar, EndOfMonthRule, adjust, advance},
    currency::Currency,
    date::Date,
    fx::CcyPair,
    id::ConventionId,
    refdata::{FxSpotLag, RefData},
    tenor::{DateOffset, TenorLength, TenorUnit},
    time::{Instant, LocalEventTime, ZoneRules},
};
use rhodium_market::{
    curve::{CurveExpr, CurveId, CurveRole, CurveSpec, Extrapolation},
    fx::{FxAnchor, FxCollateral, FxSpec},
    quote::{QuoteId, QuoteKind},
    structure::BlockId,
};
use rhodium_products::standard::InstrumentTemplate;

use super::{CalibrationTarget, ElementSelector, InstrumentBinding, QuoteTransform};
use crate::resolve::CollateralKey;

/// スポットから受渡日までの FX ポイント。
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FxDatedQuote {
    /// 両通貨の営業日で指定する受渡日。
    pub settlement: Date,
    /// pip 単位の quote。
    pub quote: QuoteId,
}

/// 一つの直接ペアの市場入力。
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FxQuotes {
    /// スポットのアウトライト。
    pub spot: QuoteId,
    /// 約定日から次の営業日までの pip。
    pub on: Option<QuoteId>,
    /// 翌営業日から次の営業日までの pip。
    pub tn: Option<QuoteId>,
    /// スポットから次の営業日までの pip。
    pub sn: Option<QuoteId>,
    /// スポットより後の受渡に対する pip。
    #[serde(default)]
    pub forwards: Vec<FxDatedQuote>,
}

/// FX ポイントと `MtM` XCCY から作る外貨担保割引曲線。
///
/// 短期は FX スワップのポイント（スポットと短期の quote、日付付きのポイント）、長期は `MtM` XCCY の
/// ベーシスに校正する。どちらの行も同じ `basis` ブロックの節点を決め、一つのグループで解く。
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FxImpliedCurve {
    /// 市場の担保の割引の組にある、今回作る側の曲線 ID。
    pub id: CurveId,
    /// 最後のピラー以後の補間規則。
    #[serde(default)]
    pub extrapolation: Extrapolation,
    /// FX ポイントより先の、`MtM` XCCY のベーシスで決めるピラー。ピラーの昇順。
    #[serde(default)]
    pub long_end: Vec<FxImpliedLongEnd>,
}

/// 外貨担保曲線の長期の一つのピラーと、それを決める `MtM` XCCY のベーシスの quote。
///
/// XCCY の quote は外貨担保曲線の `basis` ブロックを当てる。ブロックの要素の番号は当日受渡の締めの
/// 前後で（ON の区間が抜けて）ずれるので、一般の binding から要素を指さず、ここで宣言する。
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FxImpliedLongEnd {
    /// ピラーの UTC の瞬間。最後の FX ポイントの受渡より後。
    pub pillar: Instant,
    /// ベーシスの quote（スプレッド）。
    pub quote: QuoteId,
    /// `XccySwap` の慣行。払い脚にベーシスを付け、もう一方の脚の元本を FX でリセットする。
    pub convention: ConventionId,
    /// スポットからの満期。
    pub tenor: TenorLength,
    /// 校正商品の担保。両通貨の割引は価格設定の割引の表で解決する。
    pub collateral: CollateralKey,
}

impl FxImpliedCurve {
    /// FX ポイントと XCCY の行が当てる θ ブロック。
    ///
    /// # Errors
    ///
    /// 曲線 ID からブロック ID を作れない場合。
    pub fn block(&self) -> Result<BlockId> {
        BlockId::new(self.id.clone(), "basis")
    }

    /// 長期の XCCY の行を、外貨担保曲線のブロックを当てる校正の binding にする。
    pub(crate) fn long_end_bindings(&self) -> Result<Vec<InstrumentBinding>> {
        let block = self.block()?;
        self.long_end
            .iter()
            .map(|row| {
                InstrumentBinding::new(
                    row.quote.clone(),
                    QuoteKind::Spread,
                    CalibrationTarget::ParRate(InstrumentTemplate::XccyBasis {
                        convention: row.convention.clone(),
                        tenor: row.tenor,
                    }),
                    QuoteTransform::Identity,
                    row.collateral.clone(),
                    ElementSelector::Block(block.clone()),
                )
            })
            .collect()
    }
}

/// 一つの直接ペアの起点と、その市場担保の校正。
///
/// スポットラグ、同日受渡の締めと pip の桁は、参照データのペアの
/// [`FxConvention`](rhodium_core::refdata::FxConvention) から読む。
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DirectFx {
    /// ペアの市場の担保区分。
    pub market_collateral: CurveId,
    /// 市場の FX quote。
    pub quotes: FxQuotes,
    /// FX ポイントを再現する外貨担保割引曲線。スポットのみなら省略できる。
    pub implied: Option<FxImpliedCurve>,
}

/// 評価時点の FX を定める一つの通貨ペア。
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FxDefinition {
    /// 市場の標準の向き。
    pub pair: CcyPair,
    /// 直接ペアの市場入力。省略する場合は三角経路が必要。
    pub direct: Option<DirectFx>,
    /// 当該ペアを読む担保別の割引曲線。
    pub collaterals: Vec<FxCollateral>,
    /// 宣言した三角経路。直接ペアでは整合性の検査に使う。
    pub via: Option<Currency>,
}

pub(crate) struct ResolvedFxDefinition {
    pub(crate) spec: FxSpec,
    pub(crate) marks: Vec<FxMark>,
    pub(crate) fit: Option<FxFit>,
    pub(crate) warnings: Vec<Warning>,
}

pub(crate) struct FxMark {
    pub(crate) quote: QuoteId,
    pub(crate) scale: f64,
    pub(crate) positive: bool,
}

pub(crate) struct FxFit {
    pub(crate) curve: CurveSpec,
    pub(crate) target: BlockId,
    pub(crate) collateral: CurveId,
    pub(crate) intervals: Vec<FxQuoteInterval>,
    pub(crate) pip: f64,
}

pub(crate) struct FxQuoteInterval {
    pub(crate) start: Instant,
    pub(crate) end: Instant,
    pub(crate) quote: QuoteId,
}

impl FxDefinition {
    /// 参照データのペアの慣行と、宣言された短期の quote が合うかを確かめる。
    ///
    /// quote は慣行の標準の向きで読むので、逆向きのペアの定義は拒否する。
    /// T+1 のペアは ON だけ、T+2 のペアは TN（と ON）でスポットへ連鎖する。
    ///
    /// # Errors
    ///
    /// ペアの慣行が無いか、向きが慣行と違うか、短期の quote がスポットへ連鎖しなければ
    /// `InvalidData`。
    pub fn validate_references(&self, refdata: &RefData) -> Result<()> {
        let Some(direct) = &self.direct else {
            return Ok(());
        };
        let convention = refdata.fx_convention(self.pair)?;
        if convention.pair != self.pair {
            return Err(Error::InvalidData(format!(
                "market FX quotes for {} must use the convention direction {}",
                self.pair, convention.pair
            )));
        }
        let quotes = &direct.quotes;
        let chain = match convention.spot_lag {
            FxSpotLag::OneBusinessDay => quotes.tn.is_none(),
            FxSpotLag::TwoBusinessDays => quotes.on.is_none() || quotes.tn.is_some(),
        };
        if chain {
            Ok(())
        } else {
            Err(Error::InvalidData(format!(
                "short FX quotes do not form a chain to spot for {}",
                self.pair
            )))
        }
    }

    /// 宣言された quote だけで決まる定義の整合性を、評価時点と参照データによらずに確かめる。
    pub(crate) fn validate(&self) -> Result<()> {
        let Some(direct) = &self.direct else {
            return Ok(());
        };
        if let Some(implied) = &direct.implied
            && !implied
                .long_end
                .windows(2)
                .all(|pair| pair[0].pillar < pair[1].pillar)
        {
            return Err(Error::InvalidData(format!(
                "long-end pillars of {} must be strictly increasing",
                implied.id
            )));
        }
        let quotes = &direct.quotes;
        // ON は当日受渡の締めを過ぎると連鎖から外れる。締めの前後で定義の成否が変わらないよう、
        // 外貨担保曲線の要否は使う区間ではなく宣言された quote で決める。
        // TN の区間は翌営業日に始まり、SN と日付付きの区間はスポットに始まるので、締め後も残る。
        let after_cutoff =
            quotes.tn.is_some() || quotes.sn.is_some() || !quotes.forwards.is_empty();
        let points = quotes.on.is_some() || after_cutoff;
        match (&direct.implied, points, after_cutoff) {
            // ポイントを市場の担保のフォワードで再現する自由度が無い（設計の fx_x0_on_tn_consistency）。
            (None, true, _) => Err(Error::InvalidData(format!(
                "market FX points for {} require an FX-implied curve",
                self.pair
            ))),
            // T+1 のペアの ON だけでは、締め後に外貨担保曲線の校正区間が一つも残らない。
            (Some(_), _, false) => Err(Error::InvalidData(format!(
                "FX-implied curve for {} needs a TN, SN or dated quote that remains after the same-day cutoff",
                self.pair
            ))),
            _ => Ok(()),
        }
    }

    pub(crate) fn quote_kinds(&self) -> impl Iterator<Item = (&QuoteId, QuoteKind)> {
        self.direct.iter().flat_map(|direct| {
            let quotes = &direct.quotes;
            let long_end = direct
                .implied
                .iter()
                .flat_map(|implied| &implied.long_end)
                .map(|row| (&row.quote, QuoteKind::Spread));
            std::iter::once((&quotes.spot, QuoteKind::FxSpot))
                .chain(
                    [&quotes.on, &quotes.tn, &quotes.sn]
                        .into_iter()
                        .filter_map(Option::as_ref)
                        .chain(quotes.forwards.iter().map(|row| &row.quote))
                        .map(|quote| (quote, QuoteKind::FxForwardPoints)),
                )
                .chain(long_end)
        })
    }

    pub(crate) fn resolve(
        &self,
        as_of: Instant,
        refdata: &RefData,
    ) -> Result<ResolvedFxDefinition> {
        let mut spec = FxSpec {
            pair: self.pair,
            anchor: FxAnchor::Triangulated,
            collaterals: self.collaterals.clone(),
            via: self.via,
        };
        let Some(direct) = &self.direct else {
            return Ok(ResolvedFxDefinition {
                spec,
                marks: Vec::new(),
                fit: None,
                warnings: Vec::new(),
            });
        };
        self.validate_references(refdata)?;
        let mut warnings = Vec::new();
        let pair = self.pair;
        let convention = refdata.fx_convention(pair)?;
        let pip = convention.pip();
        let calendar = FxCalendar::new(pair, &convention.settlement, refdata)?;
        let local_date = calendar.zone.to_local(as_of)?.date();
        let (trade, adjust_warnings) = adjust(
            local_date,
            BusinessDayConvention::Following,
            &calendar.joint,
        )?;
        warnings.extend(adjust_warnings);
        let trade_at = calendar.instant(trade)?;
        let tomorrow = calendar.next(trade, &mut warnings)?;
        let tomorrow_at = calendar.instant(tomorrow)?;
        let spot_date = match convention.spot_lag {
            FxSpotLag::OneBusinessDay => tomorrow,
            FxSpotLag::TwoBusinessDays => calendar.next(tomorrow, &mut warnings)?,
        };
        let spot = calendar.instant(spot_date)?;
        let mut intervals = Vec::new();
        let mut marks = vec![FxMark {
            quote: direct.quotes.spot.clone(),
            scale: 1.0,
            positive: true,
        }];
        let short = [
            (trade_at, tomorrow_at, direct.quotes.on.as_ref()),
            (tomorrow_at, spot, direct.quotes.tn.as_ref()),
        ];
        for (start, end, quote) in short {
            if let Some(quote) = quote.filter(|_| start >= as_of) {
                marks.push(FxMark {
                    quote: quote.clone(),
                    scale: pip,
                    positive: false,
                });
                intervals.push(FxQuoteInterval {
                    start,
                    end,
                    quote: quote.clone(),
                });
            }
        }
        spec.anchor = FxAnchor::Direct {
            spot,
            short_starts: intervals.iter().map(|interval| interval.start).collect(),
            market_collateral: direct.market_collateral.clone(),
        };
        if let Some(sn) = &direct.quotes.sn {
            intervals.push(FxQuoteInterval {
                start: spot,
                end: calendar.instant(calendar.next(spot_date, &mut warnings)?)?,
                quote: sn.clone(),
            });
        }
        for row in &direct.quotes.forwards {
            let (valid, notices) = calendar.joint.is_business_day(row.settlement)?;
            warnings.extend(notices);
            if !valid || row.settlement <= spot_date {
                return Err(Error::InvalidData(format!(
                    "FX quote {} needs a joint business day after spot",
                    row.quote
                )));
            }
            intervals.push(FxQuoteInterval {
                start: spot,
                end: calendar.instant(row.settlement)?,
                quote: row.quote.clone(),
            });
        }
        intervals.sort_by_key(|interval| interval.end);
        if intervals.windows(2).any(|pair| pair[0].end == pair[1].end) {
            return Err(Error::InvalidData(format!(
                "multiple FX targets at the same settlement for {pair}"
            )));
        }
        let fit = self.resolve_fit(direct, intervals, pip)?;
        Ok(ResolvedFxDefinition {
            spec,
            marks,
            fit,
            warnings,
        })
    }

    fn resolve_fit(
        &self,
        direct: &DirectFx,
        intervals: Vec<FxQuoteInterval>,
        pip: f64,
    ) -> Result<Option<FxFit>> {
        let pair = self.pair;
        match (&direct.implied, intervals.is_empty()) {
            (None, true) => Ok(None),
            (Some(implied), false) => {
                let collateral = self
                    .collaterals
                    .iter()
                    .find(|collateral| collateral.id == direct.market_collateral)
                    .ok_or_else(|| {
                        Error::InvalidData(format!("missing FX market collateral for {pair}"))
                    })?;
                let (ccy, opposite) = match (&collateral.base, &collateral.quote) {
                    (base, quote) if *base == implied.id => (pair.base(), quote),
                    (base, quote) if *quote == implied.id => (pair.quote(), base),
                    _ => {
                        return Err(Error::InvalidData(format!(
                            "FX-implied curve {} is outside the market collateral of {pair}",
                            implied.id
                        )));
                    }
                };
                let target = implied.block()?;
                // 短期は FX ポイントの受渡、長期は XCCY のピラー。ピラーの順に一つずつの行が決める。
                let last = intervals.iter().map(|interval| interval.end).max();
                if let (Some(last), Some(first)) = (last, implied.long_end.first())
                    && first.pillar <= last
                {
                    return Err(Error::InvalidData(format!(
                        "long-end pillar of {} must follow the last FX point settlement",
                        implied.id
                    )));
                }
                let curve = CurveSpec::from_expr(
                    implied.id.clone(),
                    ccy,
                    [CurveRole::Discount],
                    CurveExpr::FxImplied {
                        part: target.part().into(),
                        pair,
                        base: opposite.clone(),
                        pillars: intervals
                            .iter()
                            .map(|interval| interval.end)
                            .chain(implied.long_end.iter().map(|row| row.pillar))
                            .collect(),
                        extrapolation: implied.extrapolation,
                    },
                )?;
                Ok(Some(FxFit {
                    curve,
                    target,
                    collateral: direct.market_collateral.clone(),
                    intervals,
                    pip,
                }))
            }
            // `validate` が宣言された quote で拒否するので、次の二つは定義の検査を経ない場合の防御。
            (None, false) => Err(Error::InvalidData(format!(
                "market FX points for {pair} require an FX-implied curve"
            ))),
            (Some(_), true) => Err(Error::InvalidData(format!(
                "FX-implied curve for {pair} has no eligible calibration interval"
            ))),
        }
    }
}

struct FxCalendar<'a> {
    joint: Calendar,
    settlement: &'a LocalEventTime,
    zone: &'a ZoneRules,
}

impl<'a> FxCalendar<'a> {
    fn new(pair: CcyPair, settlement: &'a LocalEventTime, refdata: &'a RefData) -> Result<Self> {
        Ok(Self {
            joint: refdata.calendar(&refdata.fx_calendar_spec(pair)?)?,
            settlement,
            zone: refdata.zone(&settlement.zone)?,
        })
    }

    fn instant(&self, date: Date) -> Result<Instant> {
        self.zone.resolve_convention(self.settlement.on_date(date))
    }

    fn next(&self, date: Date, warnings: &mut Vec<Warning>) -> Result<Date> {
        let rule = AdvanceRule {
            calendar: &self.joint,
            bdc: BusinessDayConvention::Following,
            end_of_month: EndOfMonthRule::Ignore,
        };
        let (date, notices) = advance(date, DateOffset::new(1, TenorUnit::BusinessDay), &rule)?;
        warnings.extend(notices);
        Ok(date)
    }
}
`````

`````
crates/rhodium-pricing/src/lib.rs
//! 市場の校正、評価の準備と価格計算。

#![forbid(unsafe_code)]

// 校正の単体テストでも統合テストと同じ市場 fixture を使う。
#[cfg(test)]
extern crate self as rhodium_pricing;

pub mod calibration;
pub mod config;
pub mod context;
pub mod definition;
mod linear;
pub mod par;
pub mod prepare;
pub mod resolve;
pub mod valuation;
`````

`````
crates/rhodium-pricing/src/linear.rs
//! リニア方式の cashflow の計画と評価。

mod future;
mod fx;
mod rfr;

use std::{
    collections::{BTreeMap, BTreeSet, HashMap, btree_map::Entry, hash_map},
    sync::Arc,
};

use rhodium_core::{
    Error, Result, Warning,
    currency::{Asset, Currency},
    fixing::ObservationPoint,
    fx::CcyPair,
    id::IndexId,
    money::{Money, MoneyBag},
    rounding::Rounding,
    time::Instant,
};
use rhodium_market::{
    curve::{CurvePoint, CurveSpec},
    fx::FxForwardPlan,
    structure::{BlockId, MarketStructure},
    view::{BoundSlots, Discount, MarketView, ResolvedCurve, Slots},
};
use rhodium_products::ir::{
    Product,
    analysis::{Cashflow, ExprRef, Expression, FxFixing, cashflows},
};
use rhodium_products::{
    standard::{PayReceive, RateFutureMargin},
    trade::Fee,
};
use scioxide_ad::{Mask, Real};

use crate::{
    config::LinearConvention,
    context::Context,
    resolve::{CollateralKey, Requirement, ResolvedRequirement, resolve},
};
use future::FuturePrice;
use rfr::{CompoundPlan, Forecast};

#[derive(Debug)]
pub(super) struct CashflowPlan {
    groups: Vec<PaymentGroup>,
    fx: Vec<FxForwardPlan>,
    futures: Vec<FuturePrice>,
}

/// OIS の方程式が、このピラーで一意な根を持ちうることを示す。
#[derive(Debug, Clone, Copy)]
pub(super) struct OisWarmStart {
    pub(super) pillar: usize,
    pub(super) requires_nonnegative_quote: bool,
}

struct OisFloatingPayment<'a> {
    curve: &'a CurveSpec,
    discount: &'a CurveSpec,
    start: Instant,
    end: Instant,
    pay: Instant,
    coefficient: f64,
    fixed_excess: f64,
}

#[derive(Debug)]
struct PaymentGroup {
    currency: Currency,
    payments: Vec<Payment>,
}

#[derive(Debug)]
struct Payment {
    point: usize,
    amount: Vec<Instruction>,
    /// 金額の式が読む複利の計画。`Compound` と `CompoundExcess` が添字で共有する。
    compounds: Vec<CompoundPlan>,
}

#[derive(Debug)]
enum Instruction {
    Constant(f64),
    Negate(usize),
    Round(Rounding, usize),
    Add(Vec<usize>),
    Subtract(usize, usize),
    Multiply(Vec<usize>),
    Divide(usize, usize),
    FraDivide(usize, usize),
    FxForward(usize),
    RateFuture(usize),
    /// ターム金利の予想値を読むスロットと点。
    Forward {
        slot: usize,
        point: usize,
    },
    /// 複利の成長率 G。値は `Payment::compounds` の計画で求める。
    Compound(usize),
    /// 超過成長 G − 1。値は G と同じ計画（観測の行ごとの解決）から、予想する部分は対数の和と
    /// `exp_m1`、確定部分は簡約が残した超過で求める。
    CompoundExcess(usize),
}

impl CashflowPlan {
    /// 単位固定金利で準備したスポット起点の標準 OIS の、単調な方程式を見分ける。
    pub(super) fn ois_warm_start(&self, slots: &Slots, target: &BlockId) -> Option<OisWarmStart> {
        let [group] = self.groups.as_slice() else {
            return None;
        };
        let mut fixed = Vec::new();
        let mut floating = Vec::new();
        let mut values = Vec::new();
        for payment in &group.payments {
            let (discount, at) = slots.discount_point(0, payment.point)?;
            // 固定脚は負の定数、変動脚は正の係数を掛けた DF 比の超過になる。
            let amount = payment
                .evaluate_amount(
                    &mut values,
                    |_, excess| {
                        if excess {
                            Ok(1.0_f64)
                        } else {
                            Err(Error::Internal("expected excess growth".into()))
                        }
                    },
                    |_, _| Err(Error::Internal("unexpected term rate in OIS".into())),
                    |_| Err(Error::Internal("unexpected FX rate in OIS".into())),
                    |_| Err(Error::Internal("unexpected futures price in OIS".into())),
                )
                .ok()?;
            if !amount.is_finite() {
                return None;
            }
            match payment.compounds.as_slice() {
                [] if amount < 0.0 => fixed.push(at),
                [compound] if amount > 0.0 => {
                    let (slot, point, fixed_excess) = compound.pure_growth()?;
                    let (curve, start, end) = slots.growth_interval(slot, point)?;
                    if target.part() != "nodes" || curve.id() != target.curve() || start >= end {
                        return None;
                    }
                    floating.push(OisFloatingPayment {
                        curve,
                        discount,
                        start,
                        end,
                        pay: at,
                        coefficient: amount,
                        fixed_excess,
                    });
                }
                _ => return None,
            }
        }
        floating.sort_by_key(|period| period.start);
        let first = floating.first()?;
        let last = floating.last()?;
        let pillars = last.curve.pillars();
        let pillar = pillars.partition_point(|&time| time < last.end);
        if pillar == pillars.len() {
            return None;
        }

        // 同時刻の単一支払なら DF は分子・分母で消える。ln DF の線形補間では、
        // 最新ピラーの係数は期間の終点のほうが始点より大きいため、パーは厳密に単調。
        if floating.len() == 1 && fixed.as_slice() == [last.pay] {
            return Some(OisWarmStart {
                pillar,
                requires_nonnegative_quote: false,
            });
        }
        // DF 比が連続し同じ曲線で各終了時刻に割り引けば、変動脚は c(P(start) − P(end))。
        // P(start) が既知で q >= 0 なら、固定脚を引いた残差は最新ピラーに対して厳密に減少する。
        let previous_pillar = pillars.get(pillar.checked_sub(1)?)?;
        (first.start <= *previous_pillar
            && !fixed.is_empty()
            && fixed.iter().all(|&pay| pay <= last.end)
            && floating.windows(2).all(|pair| pair[0].end == pair[1].start)
            && floating.iter().all(|period| {
                period.end == period.pay
                    && period.discount.id() == target.curve()
                    && period.coefficient.total_cmp(&first.coefficient).is_eq()
                    && period.fixed_excess == 0.0
            }))
        .then_some(OisWarmStart {
            pillar,
            requires_nonnegative_quote: true,
        })
    }

    pub(super) fn evaluate<'v, R: Real>(
        &self,
        view: &MarketView<'v, R>,
        slots: &BoundSlots<'v, R>,
    ) -> Result<MoneyBag<R>> {
        let mut result = MoneyBag::new();
        let mut values = Vec::new();
        let fx = self
            .fx
            .iter()
            .map(|plan| view.fx_forward_value(plan))
            .collect::<Result<Vec<_>>>()?;
        let futures = self
            .futures
            .iter()
            .map(|plan| plan.evaluate(view))
            .collect::<Result<Vec<_>>>()?;
        for (group, discount) in self.groups.iter().zip(slots.discounts()) {
            let money = view.ccy_scope(group.currency, discount.curve(), |scope| {
                let pay_points = discount
                    .points()
                    .iter()
                    .map(|&point| scope.bind_pay(point))
                    .collect::<Result<Vec<_>>>()?;
                let mut total = scope.present(R::zero());
                for payment in &group.payments {
                    let amount = payment.evaluate(view, slots, &mut values, &fx, &futures)?;
                    total = total + scope.pv(amount, pay_points[payment.point]);
                }
                Ok(scope.into_money(total))
            })?;
            result.add(checked_money(money)?);
        }
        Ok(result)
    }
}

/// 複利から 1 を引く `Compounded − 1` なら、その複利の式を返す。
///
/// 簡約は確定因子を `Compounded` の超過 `fixed_excess` に畳むので、一部が確定した期間もこの形のままになる。
fn excess_growth(expr: ExprRef<'_>) -> Option<ExprRef<'_>> {
    if !matches!(expr.expression(), Expression::Subtract) {
        return None;
    }
    let mut children = expr.children();
    let (growth, one) = (children.next()?, children.next()?);
    (matches!(growth.expression(), Expression::Compounded(_))
        && matches!(one.expression(), Expression::Constant(value) if value.to_bits() == 1.0_f64.to_bits()))
    .then_some(growth)
}

fn inverse_fx(expr: ExprRef<'_>) -> Option<ExprRef<'_>> {
    if !matches!(expr.expression(), Expression::Divide) {
        return None;
    }
    let denominator = expr.children().next_back()?;
    matches!(denominator.expression(), Expression::FxFixing(_)).then_some(denominator)
}

fn checked_money<R: Real>(money: Money<R>) -> Result<Money<R>> {
    if money.amount().is_finite().all() {
        Ok(money)
    } else {
        Err(Error::Numerical(format!(
            "non-finite PV in {}",
            money.ccy()
        )))
    }
}

impl Payment {
    fn evaluate<'v, R: Real>(
        &self,
        view: &MarketView<'v, R>,
        slots: &BoundSlots<'v, R>,
        values: &mut Vec<R>,
        fx: &[R],
        futures: &[R],
    ) -> Result<R> {
        self.evaluate_amount(
            values,
            |compound, excess| {
                if excess {
                    compound.evaluate_excess(view, slots)
                } else {
                    compound.evaluate(view, slots)
                }
            },
            |slot, point| Ok(view.forward(slots.projections()[slot].points()[point])),
            |point| Ok(fx[point]),
            |point| Ok(futures[point]),
        )
    }

    fn evaluate_amount<R: Real>(
        &self,
        values: &mut Vec<R>,
        compound_value: impl Fn(&CompoundPlan, bool) -> Result<R>,
        forward_value: impl Fn(usize, usize) -> Result<R>,
        fx_value: impl Fn(usize) -> Result<R>,
        futures_value: impl Fn(usize) -> Result<R>,
    ) -> Result<R> {
        values.clear();
        for instruction in &self.amount {
            let value = match instruction {
                Instruction::Constant(value) => R::constant(*value),
                Instruction::Negate(x) => -values[*x],
                Instruction::Round(rounding, x) => rounding.apply(values[*x])?,
                Instruction::Add(xs) => xs.iter().fold(R::zero(), |sum, &x| sum + values[x]),
                Instruction::Subtract(x, y) => values[*x] - values[*y],
                Instruction::Multiply(xs) => {
                    xs.iter().fold(R::one(), |product, &x| product * values[x])
                }
                Instruction::Divide(x, y) => values[*x] / values[*y],
                Instruction::FraDivide(x, y) => {
                    if !values[*y].gt(R::zero()).all() {
                        return Err(Error::Numerical(
                            "nonpositive projected FRA settlement factor".into(),
                        ));
                    }
                    values[*x] / values[*y]
                }
                Instruction::FxForward(point) => fx_value(*point)?,
                Instruction::RateFuture(point) => futures_value(*point)?,
                Instruction::Forward { slot, point } => forward_value(*slot, *point)?,
                Instruction::Compound(plan) | Instruction::CompoundExcess(plan) => compound_value(
                    &self.compounds[*plan],
                    matches!(instruction, Instruction::CompoundExcess(_)),
                )?,
            };
            values.push(value);
        }
        values
            .last()
            .copied()
            .ok_or_else(|| Error::Internal("empty cashflow expression".into()))
    }
}

/// 一つの支払いの金額の式が読む複利の計画と、複利の式の番号からその添字への表。
#[derive(Default)]
struct CompoundPlans {
    plans: Vec<CompoundPlan>,
    by_expr: HashMap<usize, usize>,
}

struct DiscountSlot {
    curve: ResolvedCurve<Discount>,
    points: Vec<CurvePoint>,
    times: HashMap<Instant, usize>,
    group: PaymentGroup,
}

pub(super) struct LinearPreparation {
    pub(super) plan: CashflowPlan,
    pub(super) slots: Slots,
    pub(super) blocks: Vec<BlockId>,
    pub(super) warnings: Vec<Warning>,
}

pub(super) fn prepare(
    product: &Product,
    collateral: &CollateralKey,
    structure: &Arc<MarketStructure>,
    ctx: &Context<'_>,
) -> Result<LinearPreparation> {
    let mut builder = PlanBuilder::new(
        collateral,
        structure,
        ctx,
        product.warnings().iter().cloned(),
    );
    for payment in cashflows(product) {
        builder.payment(product, &payment)?;
    }
    builder.finish()
}

pub(super) fn prepare_rate_future(
    margin: &RateFutureMargin<'_>,
    fees: &[Fee],
    warnings: Vec<Warning>,
    collateral: &CollateralKey,
    structure: &Arc<MarketStructure>,
    ctx: &Context<'_>,
) -> Result<LinearPreparation> {
    use rhodium_core::fixing::{FixingResolution, ObservationPeriod};

    let mut builder = PlanBuilder::new(collateral, structure, ctx, warnings);
    // 取引の手数料は値洗いと独立の受渡で、IR の評価と同じく決済の境界より後のものだけを数える。
    for fee in fees
        .iter()
        .filter(|fee| ctx.config.settlement().includes(ctx.as_of, fee.at))
    {
        let signed = match fee.side {
            PayReceive::Pay => -fee.amount.amount(),
            PayReceive::Receive => fee.amount.amount(),
        };
        let discount = builder.discount(fee.amount.ccy())?;
        let point = builder.payment_point(discount, fee.at)?;
        builder.discounts[discount].group.payments.push(Payment {
            point,
            amount: vec![Instruction::Constant(signed)],
            compounds: Vec::new(),
        });
    }
    let first = margin
        .clearings
        .partition_point(|row| !ctx.config.settlement().includes(ctx.as_of, row.settlement));
    if first == margin.clearings.len() {
        return builder.finish();
    }
    let resolve = |row: &rhodium_core::refdata::RateFutureClearing,
                   warnings: &mut BTreeSet<Warning>|
     -> Result<Option<f64>> {
        let (fixing, warning) = ctx.fixings.resolve(
            margin.index,
            ObservationPoint::Date(row.date),
            ctx.as_of,
            row.published_at,
            ObservationPeriod::at(row.observed_at),
        )?;
        warnings.extend(warning);
        Ok(match fixing {
            FixingResolution::Published(value) | FixingResolution::Estimated(value) => Some(value),
            FixingResolution::Forecast | FixingResolution::ForecastFrom(_) => None,
        })
    };
    let mut previous = if first == 0 {
        Some(margin.trade_price)
    } else {
        resolve(&margin.clearings[first - 1], &mut builder.warnings)?
    };
    for row in &margin.clearings[first..] {
        let price = resolve(row, &mut builder.warnings)?;
        if price.is_none() {
            future::clearing_delay(*row, margin.currency, ctx, &mut builder.warnings)?;
        }
        if price.is_some() || previous.is_some() {
            let discount = builder.discount(margin.currency)?;
            let point = builder.payment_point(discount, row.settlement)?;
            let amount = match (price, previous) {
                (Some(price), Some(previous)) => vec![Instruction::Constant(
                    margin.signed_multiplier * (price - previous),
                )],
                (price, previous) => {
                    let forecast = builder.future_price(margin.index, discount)?;
                    let (known, sign) = match (price, previous) {
                        (None, Some(value)) => (value, 1.0),
                        (Some(value), None) => (value, -1.0),
                        _ => return Err(Error::Internal("expected one known margin price".into())),
                    };
                    vec![
                        Instruction::RateFuture(forecast),
                        Instruction::Constant(known),
                        Instruction::Subtract(0, 1),
                        Instruction::Constant(sign * margin.signed_multiplier),
                        Instruction::Multiply(vec![2, 3]),
                    ]
                }
            };
            builder.discounts[discount].group.payments.push(Payment {
                point,
                amount,
                compounds: Vec::new(),
            });
        }
        previous = price;
    }
    builder.finish()
}

struct PlanBuilder<'a> {
    collateral: &'a CollateralKey,
    structure: &'a Arc<MarketStructure>,
    ctx: &'a Context<'a>,
    discount_indices: HashMap<Currency, usize>,
    discounts: Vec<DiscountSlot>,
    forecast_indices: HashMap<IndexId, usize>,
    forecasts: Vec<Forecast<'a>>,
    fx_indices: HashMap<(IndexId, Instant, CcyPair), usize>,
    fx: Vec<FxForwardPlan>,
    futures_indices: HashMap<IndexId, usize>,
    futures: Vec<FuturePrice>,
    warnings: BTreeSet<Warning>,
}

impl<'a> PlanBuilder<'a> {
    fn new(
        collateral: &'a CollateralKey,
        structure: &'a Arc<MarketStructure>,
        ctx: &'a Context<'a>,
        warnings: impl IntoIterator<Item = Warning>,
    ) -> Self {
        Self {
            collateral,
            structure,
            ctx,
            discount_indices: HashMap::new(),
            discounts: Vec::new(),
            forecast_indices: HashMap::new(),
            forecasts: Vec::new(),
            fx_indices: HashMap::new(),
            fx: Vec::new(),
            futures_indices: HashMap::new(),
            futures: Vec::new(),
            warnings: warnings.into_iter().collect(),
        }
    }

    fn payment(&mut self, product: &Product, payment: &Cashflow<'_>) -> Result<()> {
        let Asset::Cash(currency) = *payment.asset else {
            return Err(Error::Internal(
                "linear plan requires a currency payment".into(),
            ));
        };
        let position = self.discount(currency)?;
        let point = self.payment_point(position, payment.at)?;
        let (amount, compounds) = self.amount(product, payment.amount, position, payment.at)?;
        self.discounts[position].group.payments.push(Payment {
            point,
            amount,
            compounds,
        });
        Ok(())
    }

    fn payment_point(&mut self, slot: usize, at: Instant) -> Result<usize> {
        let discount = &mut self.discounts[slot];
        if let Some(&point) = discount.times.get(&at) {
            return Ok(point);
        }
        let point = discount.points.len();
        discount.points.push(discount.curve.prepare_point(at)?);
        discount.times.insert(at, point);
        Ok(point)
    }

    fn future_price(&mut self, index: &IndexId, discount: usize) -> Result<usize> {
        if let Some(&slot) = self.futures_indices.get(index) {
            return Ok(slot);
        }
        let slot = self.futures.len();
        let discount = self.discounts[discount].curve.spec().id();
        self.futures.push(FuturePrice::prepare(
            index,
            discount,
            self.structure,
            self.ctx,
            &mut self.warnings,
        )?);
        self.futures_indices.insert(index.clone(), slot);
        Ok(slot)
    }

    fn fx_point(
        &mut self,
        fixing: FxFixing<'_>,
        currency: Currency,
        pay: Instant,
    ) -> Result<usize> {
        let value = fx::value_time(fixing, self.ctx, &mut self.warnings)?;
        // 予想はスポットの受渡で読むので、その受渡と実際に受け渡す瞬間の遅れを照合する。
        // 通常の fixing は支払いの瞬間、元本リセットは新しい元本を交換する瞬間と比べる。リセットの
        // 元本に掛かる期末のクーポンの遅れは、金利の予想との積の規約 MtmFxForwardProduct が扱う。
        let exchanged = fixing.value_at.unwrap_or(pay);
        if value != exchanged {
            let terms = self.ctx.refdata.currency(currency)?;
            let date = self
                .ctx
                .refdata
                .zone(&terms.settlement_time.zone)?
                .to_local(value)?
                .date();
            let warning =
                rfr::payment_delay(self.ctx, currency, date, exchanged, &mut self.warnings)?;
            self.warnings.extend(warning);
        }
        match self
            .fx_indices
            .entry((fixing.index.clone(), value, fixing.pair))
        {
            hash_map::Entry::Occupied(entry) => Ok(*entry.get()),
            hash_map::Entry::Vacant(entry) => {
                let point = self.fx.len();
                self.fx.push(fx::prepare(
                    fixing,
                    value,
                    self.collateral,
                    self.structure,
                    self.ctx,
                    &mut self.warnings,
                )?);
                entry.insert(point);
                Ok(point)
            }
        }
    }

    fn discount(&mut self, currency: Currency) -> Result<usize> {
        if let Some(&position) = self.discount_indices.get(&currency) {
            return Ok(position);
        }
        let request = Requirement::Discounting {
            ccy: currency,
            collateral: self.collateral.clone(),
        };
        let ResolvedRequirement::Discounting(curve) =
            resolve(&request, self.ctx.config, self.structure, self.ctx.refdata)?
        else {
            return Err(Error::Internal("expected a discount curve".into()));
        };
        let position = self.discounts.len();
        self.discounts.push(DiscountSlot {
            curve,
            points: Vec::new(),
            times: HashMap::new(),
            group: PaymentGroup {
                currency,
                payments: Vec::new(),
            },
        });
        self.discount_indices.insert(currency, position);
        Ok(position)
    }

    fn forecast(&mut self, index: &IndexId) -> Result<usize> {
        if let Some(&position) = self.forecast_indices.get(index) {
            return Ok(position);
        }
        let request = Requirement::Forecast {
            index: index.clone(),
        };
        let ResolvedRequirement::Forecast(curve) =
            resolve(&request, self.ctx.config, self.structure, self.ctx.refdata)?
        else {
            return Err(Error::Internal("expected a projection curve".into()));
        };
        let position = self.forecasts.len();
        self.forecasts
            .push(Forecast::new(index, curve, self.ctx.refdata, position)?);
        self.forecast_indices.insert(index.clone(), position);
        Ok(position)
    }

    #[expect(
        clippy::too_many_lines,
        reason = "式の種類ごとの命令を一箇所で組み立てる"
    )]
    fn amount(
        &mut self,
        product: &Product,
        root: ExprRef<'_>,
        discount: usize,
        at: Instant,
    ) -> Result<(Vec<Instruction>, Vec<CompoundPlan>)> {
        const RATE: u8 = 1;
        const FX: u8 = 2;
        const CLEARING: u8 = 4;

        // 子の番号は親より小さい。共有した式は一度だけ計画へ写す。
        // 超過成長 `Compounded − 1` の命令は子の複利を自分で評価するので、子をたどらない。
        // 成長率の命令も作ると、日次積の計画では同じ観測の列を二度評価してしまう。
        // 同じ複利を他の親も読むときだけ、その親からたどって成長率の命令を作り、計画を共有する。
        let mut reached = BTreeMap::new();
        let currency = self.discounts[discount].curve.spec().currency();
        let mut work = vec![(root, false)];
        while let Some((expr, inverse)) = work.pop() {
            if let Expression::FxFixing(fixing) = expr.expression() {
                let expected = if inverse {
                    fixing.pair.base()
                } else {
                    fixing.pair.quote()
                };
                if expected != currency {
                    return Err(Error::Unsupported(
                        "FX observation has the wrong orientation for the payment currency".into(),
                    ));
                }
            }
            if let Entry::Vacant(entry) = reached.entry(expr.index()) {
                entry.insert(expr);
                if excess_growth(expr).is_none() {
                    if let Some(inverse) = inverse_fx(expr) {
                        work.extend(expr.children().take(1).map(|child| (child, false)));
                        work.push((inverse, true));
                    } else {
                        work.extend(expr.children().map(|child| (child, false)));
                    }
                }
            }
        }
        let positions: HashMap<_, _> = reached
            .keys()
            .enumerate()
            .map(|(pos, &ix)| (ix, pos))
            .collect();
        let mut instructions = Vec::with_capacity(reached.len());
        let mut dynamic = Vec::with_capacity(reached.len());
        let mut compounds = CompoundPlans::default();
        let is_fra = reached.values().any(|expr| expr.fra_settlement().is_some());
        for expr in reached.into_values() {
            if let Some(compound) = excess_growth(expr) {
                // 複利から 1 を引く式は、成長率を丸めてから引くと短い期間の超過の桁を落とすので、
                // 超過のまま求める命令にする（CompoundPlan::evaluate_excess）。値は θ に依存する。
                let plan = self.compound_plan(product, compound, &mut compounds, discount, at)?;
                instructions.push(Instruction::CompoundExcess(plan));
                dynamic.push(RATE);
                continue;
            }
            let children: Vec<_> = expr
                .children()
                .map(|child| positions[&child.index()])
                .collect();
            let expression = expr.expression();
            let from_children = children
                .iter()
                .fold(0, |mask, &child| mask | dynamic[child]);
            let projected = match expression {
                Expression::Compounded(_) | Expression::TermFixing(_) => RATE,
                Expression::FxFixing(_) => FX,
                Expression::RateFutureFixing(_) => CLEARING,
                _ => from_children,
            };
            let instruction = match (expression, children.as_slice()) {
                (Expression::Constant(value), []) => Instruction::Constant(value),
                (Expression::Negate, [x]) => Instruction::Negate(*x),
                (Expression::Round(rounding), [x]) => {
                    if projected != 0 {
                        self.warnings.insert(
                            self.ctx
                                .config
                                .conventions()
                                .require(LinearConvention::RoundingAtForward)?,
                        );
                    }
                    Instruction::Round(rounding, *x)
                }
                (Expression::Add, _) => Instruction::Add(children),
                (Expression::Subtract, [x, y]) => Instruction::Subtract(*x, *y),
                (Expression::Multiply, _) => {
                    if from_children == RATE | FX
                        && children.iter().any(|&child| dynamic[child] == RATE)
                        && children.iter().any(|&child| dynamic[child] == FX)
                    {
                        self.warnings.insert(
                            self.ctx
                                .config
                                .conventions()
                                .require(LinearConvention::MtmFxForwardProduct)?,
                        );
                    }
                    Instruction::Multiply(children)
                }
                (Expression::Divide, [x, y]) if expr.fra_settlement().is_some() => {
                    Instruction::FraDivide(*x, *y)
                }
                (Expression::Divide, [x, y]) => Instruction::Divide(*x, *y),
                (Expression::TermFixing(fixing), []) => {
                    let slot = self.forecast(fixing.index)?;
                    let forecast = &mut self.forecasts[slot];
                    let discount = &self.discounts[discount].curve;
                    let point = if is_fra {
                        forecast.fra_point(fixing, discount, at, self.ctx, &mut self.warnings)?
                    } else {
                        forecast.term_point(fixing, discount, at, self.ctx, &mut self.warnings)?
                    };
                    Instruction::Forward { slot, point }
                }
                (Expression::FxFixing(fixing), []) => {
                    Instruction::FxForward(self.fx_point(fixing, currency, at)?)
                }
                (Expression::RateFutureFixing(fixing), []) => {
                    let (_, settlement_currency, _) =
                        self.ctx.refdata.rate_future_underlying(fixing.index)?;
                    if settlement_currency != currency {
                        return Err(Error::Unsupported(
                            "rate futures fixing in another payment currency".into(),
                        ));
                    }
                    if let ObservationPoint::Date(date) = fixing.point {
                        let (row, notices) =
                            self.ctx.refdata.rate_future_clearing(fixing.index, date)?;
                        self.warnings.extend(notices);
                        if row.settlement == at {
                            future::clearing_delay(row, currency, self.ctx, &mut self.warnings)?;
                        }
                    }
                    Instruction::RateFuture(self.future_price(fixing.index, discount)?)
                }
                (Expression::Compounded(_), []) => Instruction::Compound(self.compound_plan(
                    product,
                    expr,
                    &mut compounds,
                    discount,
                    at,
                )?),
                _ => {
                    return Err(Error::Internal(
                        "expression does not match the linear plan".into(),
                    ));
                }
            };
            instructions.push(instruction);
            dynamic.push(projected);
        }
        Ok((instructions, compounds.plans))
    }

    /// 複利の式の計画の添字を返す。同じ式の成長率と超過成長の命令は、一つの計画を共有する。
    fn compound_plan(
        &mut self,
        product: &Product,
        expr: ExprRef<'_>,
        compounds: &mut CompoundPlans,
        discount: usize,
        at: Instant,
    ) -> Result<usize> {
        if let Some(&plan) = compounds.by_expr.get(&expr.index()) {
            return Ok(plan);
        }
        let Expression::Compounded(compound) = expr.expression() else {
            return Err(Error::Internal("expected a compounded expression".into()));
        };
        let forecast = self.forecast(compound.index)?;
        let plan = self.forecasts[forecast].plan(
            product,
            compound,
            &self.discounts[discount].curve,
            at,
            self.ctx,
            &mut self.warnings,
        )?;
        let position = compounds.plans.len();
        compounds.plans.push(plan);
        compounds.by_expr.insert(expr.index(), position);
        Ok(position)
    }

    fn finish(self) -> Result<LinearPreparation> {
        let mut blocks = BTreeSet::new();
        blocks.extend(self.fx.iter().flat_map(FxForwardPlan::blocks).cloned());
        blocks.extend(self.futures.iter().flat_map(FuturePrice::blocks).cloned());
        let mut groups = Vec::with_capacity(self.discounts.len());
        let discounts = self
            .discounts
            .into_iter()
            .map(|slot| {
                blocks.extend(slot.points.iter().flat_map(CurvePoint::blocks).cloned());
                groups.push(slot.group);
                (slot.curve, slot.points)
            })
            .collect::<Vec<_>>();
        let mut projections = Vec::with_capacity(self.forecasts.len());
        let mut growth = Vec::with_capacity(self.forecasts.len());
        // 予想の slot は CompoundPlan が持つ番号の順に並べる。点の無い slot も番号をずらさないよう残す。
        for forecast in self.forecasts {
            blocks.extend(
                forecast
                    .rates
                    .iter()
                    .flat_map(|(start, end, _)| start.blocks().chain(end.blocks()))
                    .cloned(),
            );
            blocks.extend(
                forecast
                    .growth
                    .iter()
                    .flat_map(|(start, end)| start.blocks().chain(end.blocks()))
                    .cloned(),
            );
            projections.push((forecast.curve.clone(), forecast.rates));
            growth.push((forecast.curve, forecast.growth));
        }
        let slots = Slots::from_resolved(self.structure.clone(), discounts, projections)?
            .with_projection_growth(growth)?;
        Ok(LinearPreparation {
            plan: CashflowPlan {
                groups,
                fx: self.fx,
                futures: self.futures,
            },
            slots,
            blocks: blocks.into_iter().collect(),
            warnings: self.warnings.into_iter().collect(),
        })
    }
}
`````

`````
crates/rhodium-pricing/src/linear/future.rs
//! 日次の値洗い口座を数値基準とする金利先物価格。

use std::{collections::BTreeSet, sync::Arc};

use rhodium_core::{
    Error, Result, Warning,
    calendar::{AdvanceRule, BusinessDayConvention, EndOfMonthRule, advance},
    currency::Currency,
    date::{Date, DateInterval},
    day_count::{AccrualContext, DayCount, year_fraction},
    fixing::{FixingResolution, ObservationPoint},
    id::IndexId,
    refdata::{RateFutureClearing, RateFutureRate},
    tenor::{DateOffset, TenorUnit},
    time::{Instant, model_time},
};
use rhodium_market::{
    curve::CurveId,
    structure::{BlockId, MarketStructure},
    view::{MarketView, Slots},
};
use scioxide_ad::{Mask, Real};

use crate::{
    config::LinearConvention,
    context::Context,
    resolve::{Requirement, ResolvedRequirement, resolve},
};

#[derive(Debug)]
pub(super) enum FuturePrice {
    Fixed(f64),
    Projected(Box<ProjectedPrice>),
}

#[derive(Debug)]
pub(super) struct ProjectedPrice {
    structure: Arc<MarketStructure>,
    block: BlockId,
    currency: Currency,
    slots: Slots,
    blocks: BTreeSet<BlockId>,
    steps: Vec<Step>,
    margin: Vec<usize>,
    legs: Vec<RateLeg>,
    rate: RatePlan,
}

#[derive(Clone, Copy, Debug)]
struct Step {
    dt: f64,
    regime: usize,
}

#[derive(Clone, Copy, Debug)]
struct RateLeg {
    fixing: usize,
    start: usize,
    end: usize,
    weight: f64,
    span: Span,
}

/// 予想する 1 本の率が覆う観測の期間。
#[derive(Clone, Copy, Debug)]
enum Span {
    /// 観測の期間の全体の成長 G(開始, 終了) から率を作る。
    Whole,
    /// 評価時点 A が期間をまたぐ。[A, 終了) の成長を残りの日数係数 `τ_rem` で年率化し、
    /// その率を元の期間の日数係数 τ に当てる（参照データの章の fixing の規則）。
    /// `scale` は `τ/τ_rem` で、元の期間の増分は `τ·(G − 1)/τ_rem = scale·(G − 1)`。
    /// 成長は評価時点で決まるので、fixing の瞬間を評価時点に置く。
    FromAsOf { scale: f64 },
}

#[derive(Debug, Clone, Copy)]
enum RatePlan {
    Compound { log_fixed: f64, tau: f64 },
    Sum { fixed: f64 },
}

#[derive(Debug)]
struct ForecastLeg {
    fixing: Instant,
    start: Instant,
    end: Instant,
    weight: f64,
    span: Span,
}

impl FuturePrice {
    pub(super) fn prepare(
        index: &IndexId,
        discount: &CurveId,
        structure: &Arc<MarketStructure>,
        ctx: &Context<'_>,
        warnings: &mut BTreeSet<Warning>,
    ) -> Result<Self> {
        let (convention, contract) = ctx.refdata.rate_future_contract(index)?;
        let (underlying, currency, basis) = ctx.refdata.rate_future_underlying(index)?;
        if let Some(margin_index) = &convention.margin_index {
            let margin_curve = ctx.config.forecast_curve(margin_index).ok_or_else(|| {
                Error::MissingData(format!("no curve for margin index {margin_index}"))
            })?;
            if margin_curve != discount {
                warnings.insert(
                    ctx.config
                        .conventions()
                        .require(LinearConvention::CrossCollateralProjection)?,
                );
            }
        }
        let mut legs = Vec::new();
        let rate = match contract.rate {
            RateFutureRate::Term { fixing } => {
                let (period, notices) = ctx.refdata.term_rate_dates(underlying, fixing)?;
                warnings.extend(notices);
                let tau = fraction(basis, period, ctx, warnings)?;
                let read = RateRead {
                    weight: tau.recip(),
                    accrual: tau,
                    basis,
                };
                let fixed =
                    read_rate(underlying, fixing, read, &mut legs, ctx, warnings)?.unwrap_or(0.0);
                RatePlan::Sum { fixed }
            }
            RateFutureRate::Compounded(period) | RateFutureRate::Averaged(period) => {
                let compound = matches!(contract.rate, RateFutureRate::Compounded(_));
                overnight_plan(
                    underlying, basis, period, compound, &mut legs, ctx, warnings,
                )?
            }
        };
        if legs.is_empty() {
            return Ok(Self::Fixed(
                1.0 - match rate {
                    RatePlan::Compound { log_fixed, tau } => log_fixed.exp_m1() / tau,
                    RatePlan::Sum { fixed } => fixed,
                },
            ));
        }

        let zone = ctx.refdata.zone(&convention.settlement_time.zone)?;
        let local_as_of = zone.to_local(ctx.as_of)?.date();
        let calendar = ctx.refdata.calendar(&ctx.refdata.index(index)?.calendar)?;
        let lag = convention.settlement_lag;
        let (first, notices) = advance(
            local_as_of,
            DateOffset::new(-lag.amount(), lag.unit()),
            &AdvanceRule {
                calendar: &calendar,
                bdc: BusinessDayConvention::Preceding,
                end_of_month: EndOfMonthRule::Ignore,
            },
        )?;
        warnings.extend(notices);
        let (schedule, notices) = ctx
            .refdata
            .rate_future_schedule(index, first.min(contract.final_clearing))?;
        warnings.extend(notices);
        let margin: Vec<_> = std::iter::once(ctx.as_of)
            .chain(
                schedule
                    .into_iter()
                    .map(|row| row.settlement)
                    .filter(|&at| at > ctx.as_of),
            )
            .collect();
        if margin.len() < 2 {
            return Err(Error::InvalidData(
                "forecast futures price requires a future margin settlement".into(),
            ));
        }
        Ok(Self::Projected(Box::new(ProjectedPrice::new(
            underlying, currency, structure, ctx, legs, margin, rate,
        )?)))
    }

    pub(super) fn blocks(&self) -> impl Iterator<Item = &BlockId> {
        let blocks = match self {
            Self::Fixed(_) => None,
            Self::Projected(plan) => Some(&plan.blocks),
        };
        blocks.into_iter().flatten()
    }

    pub(super) fn evaluate<R: Real>(&self, view: &MarketView<'_, R>) -> Result<R> {
        let price = match self {
            Self::Fixed(price) => R::constant(*price),
            Self::Projected(plan) => plan.evaluate(view)?,
        };
        if !price.is_finite().all() {
            return Err(Error::Numerical("non-finite futures price".into()));
        }
        Ok(price)
    }
}

impl ProjectedPrice {
    fn new(
        underlying: &IndexId,
        currency: Currency,
        structure: &Arc<MarketStructure>,
        ctx: &Context<'_>,
        legs: Vec<ForecastLeg>,
        margin_dates: Vec<Instant>,
        rate: RatePlan,
    ) -> Result<Self> {
        let block = MarketStructure::lgm_block(currency)?;
        let breaks = structure.block_times(&block).ok_or_else(|| {
            Error::MissingData(format!("missing LGM specification for {currency}"))
        })?;
        let ResolvedRequirement::Forecast(projection) = resolve(
            &Requirement::Forecast {
                index: underlying.clone(),
            },
            ctx.config,
            structure,
            ctx.refdata,
        )?
        else {
            return Err(Error::Internal(
                "expected futures underlying projection".into(),
            ));
        };
        let growth = legs
            .iter()
            .map(|leg| {
                Ok((
                    projection.prepare_point(leg.start)?,
                    projection.prepare_point(leg.end)?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut blocks = BTreeSet::from([block.clone()]);
        blocks.extend(
            growth
                .iter()
                .flat_map(|(start, end)| start.blocks().chain(end.blocks()))
                .cloned(),
        );
        let slots = Slots::from_resolved(structure.clone(), [], [])?
            .with_projection_growth([(projection, growth)])?;
        let time = |at| model_time(ctx.as_of, at).years();
        let mut times: Vec<_> = margin_dates
            .iter()
            .copied()
            .chain(legs.iter().flat_map(|leg| [leg.fixing, leg.start, leg.end]))
            .map(time)
            .collect();
        let horizon = times.iter().copied().fold(0.0, f64::max);
        times.extend(
            breaks
                .iter()
                .map(|t| t.years())
                .take_while(|&t| t < horizon),
        );
        times.sort_by(f64::total_cmp);
        times.dedup_by(|a, b| a.total_cmp(b).is_eq());
        let position = |at| {
            times
                .binary_search_by(|t| t.total_cmp(&time(at)))
                .map_err(|_| Error::Internal("missing futures model time".into()))
        };
        let margin = margin_dates
            .into_iter()
            .map(position)
            .collect::<Result<_>>()?;
        let legs = legs
            .into_iter()
            .map(|leg| {
                Ok(RateLeg {
                    fixing: position(leg.fixing)?,
                    start: position(leg.start)?,
                    end: position(leg.end)?,
                    weight: leg.weight,
                    span: leg.span,
                })
            })
            .collect::<Result<_>>()?;
        let steps = times
            .windows(2)
            .map(|pair| Step {
                dt: pair[1] - pair[0],
                regime: breaks.partition_point(|t| t.years() <= pair[0]),
            })
            .collect();
        Ok(Self {
            structure: structure.clone(),
            block,
            currency,
            slots,
            blocks,
            steps,
            margin,
            legs,
            rate,
        })
    }

    fn evaluate<R: Real>(&self, view: &MarketView<'_, R>) -> Result<R> {
        let model = view.lgm_parameters(&self.structure, self.currency, &self.block)?;
        let (reversion, vol) = model.split_at(model.len() / 2);
        if vol.iter().any(|&value| !value.ge(R::zero()).all()) {
            return Err(Error::InvalidData(format!(
                "negative LGM volatility for {}",
                self.currency
            )));
        }
        let n = self.steps.len() + 1;
        let mut h = Vec::with_capacity(n);
        let mut zeta = Vec::with_capacity(n);
        h.push(R::zero());
        zeta.push(R::zero());
        let mut derivative = R::one();
        for (i, step) in self.steps.iter().enumerate() {
            let a_dt = reversion[step.regime] * step.dt;
            h.push(h[i] + derivative * integration_factor(a_dt) * step.dt);
            derivative *= (-a_dt).exp();
            let volatility = vol[step.regime];
            zeta.push(zeta[i] + volatility * volatility * step.dt);
        }

        // log(B_M / N_LGM) のガウス部分。各区間の付利は期首に決まり、
        // 終端で LGM の数値基準を割る。指数傾斜による x の平均の移動を一回の走査で求める。
        let mut density = vec![R::zero(); n];
        for pair in self.margin.windows(2) {
            density[pair[0]] += h[pair[1]] - h[pair[0]];
        }
        let terminal = self.margin[self.margin.len() - 1];
        density[terminal] -= h[terminal];
        let mut suffix = R::zero();
        for weight in density.iter_mut().rev() {
            suffix += *weight;
            *weight = suffix;
        }
        let mut tilt = vec![R::zero(); n];
        for i in 1..n {
            tilt[i] = tilt[i - 1] + density[i] * (zeta[i] - zeta[i - 1]);
        }

        let slots = view.bind_slots(&self.slots)?;
        let growth = slots.projection_growth()[0].points();
        let log_expected = |leg: RateLeg, position: usize| {
            let loading = h[leg.end] - h[leg.start];
            let drift = loading * (h[leg.end] + h[leg.start]) * zeta[leg.fixing] * 0.5;
            (
                view.log_growth(growth[position]) + drift + loading * tilt[leg.fixing],
                loading,
            )
        };
        let rate =
            match self.rate {
                // 期間をまたぐ脚は fixing が評価時点（ζ = 0、傾き 0）なので成長は決定論的で、
                // 期待値は G そのもの。増分 G − 1 を scale 倍して元の期間の増分にする。
                RatePlan::Sum { fixed } => self.legs.iter().copied().enumerate().fold(
                    R::constant(fixed),
                    |rate, (i, leg)| {
                        let (log, loading) = log_expected(leg, i);
                        let scale = match leg.span {
                            Span::Whole => 1.0,
                            Span::FromAsOf { scale } => scale,
                        };
                        rate + (log + loading * loading * zeta[leg.fixing] * 0.5).exp_m1()
                            * (leg.weight * scale)
                    },
                ),
                RatePlan::Compound { log_fixed, tau } => {
                    let mut weights = vec![R::zero(); n];
                    let mut log = R::constant(log_fixed);
                    for (i, leg) in self.legs.iter().copied().enumerate() {
                        let (component, loading) = log_expected(leg, i);
                        match leg.span {
                            Span::Whole => {
                                log += component;
                                weights[leg.fixing] += loading;
                            }
                            // 元の期間の複利因子 1 + τ·r = 1 + scale·(G − 1)。対数正規ではないが、
                            // 決定論的なので共分散の項に加えない。
                            Span::FromAsOf { scale } => {
                                log += (component.exp_m1() * scale).ln_1p();
                            }
                        }
                    }
                    let mut suffix = R::zero();
                    for i in (1..n).rev() {
                        suffix += weights[i];
                        log += suffix * suffix * (zeta[i] - zeta[i - 1]) * 0.5;
                    }
                    log.exp_m1() / tau
                }
            };
        Ok(R::one() - rate)
    }
}

fn overnight_plan(
    underlying: &IndexId,
    basis: &DayCount,
    period: DateInterval,
    compound: bool,
    legs: &mut Vec<ForecastLeg>,
    ctx: &Context<'_>,
    warnings: &mut BTreeSet<Warning>,
) -> Result<RatePlan> {
    let tau = fraction(basis, period, ctx, warnings)?;
    let cal = ctx
        .refdata
        .calendar(&ctx.refdata.index(underlying)?.calendar)?;
    let preceding = AdvanceRule {
        end_of_month: EndOfMonthRule::Ignore,
        bdc: BusinessDayConvention::Preceding,
        calendar: &cal,
    };
    let (mut start, notices) = advance(
        period.start(),
        DateOffset::new(0, TenorUnit::Day),
        &preceding,
    )?;
    warnings.extend(notices);
    let mut fixed = 0.0;
    while start < period.end() {
        let (end, notices) = advance(
            start,
            DateOffset::new(1, TenorUnit::BusinessDay),
            &preceding,
        )?;
        warnings.extend(notices);
        let coverage = DateInterval::new(start.max(period.start()), end.min(period.end()))?;
        let covered = fraction(basis, coverage, ctx, warnings)?;
        let actual = fraction(basis, DateInterval::new(start, end)?, ctx, warnings)?;
        if compound && coverage != DateInterval::new(start, end)? {
            return Err(Error::InvalidData(
                "compounded futures reference dates must align with overnight fixing periods"
                    .into(),
            ));
        }
        let read = RateRead {
            weight: if compound {
                1.0
            } else {
                covered / (actual * tau)
            },
            accrual: actual,
            basis,
        };
        if let Some(value) = read_rate(underlying, start, read, legs, ctx, warnings)? {
            if compound {
                let simple = actual * value;
                if simple <= -1.0 {
                    return Err(Error::Numerical(
                        "nonpositive fixed futures accrual factor".into(),
                    ));
                }
                fixed += simple.ln_1p();
            } else {
                fixed += covered / tau * value;
            }
        }
        start = end;
    }
    Ok(if compound {
        RatePlan::Compound {
            log_fixed: fixed,
            tau,
        }
    } else {
        RatePlan::Sum { fixed }
    })
}

fn integration_factor<R: Real>(value: R) -> R {
    let small = value.abs().le(R::constant(1e-5));
    let safe = R::select(small, R::one(), value);
    let series = R::one()
        - value * (R::constant(0.5) - value * (R::constant(1.0 / 6.0) - value * (1.0 / 24.0)));
    R::select(small, series, -(-value).exp_m1() / safe)
}

fn fraction(
    basis: &DayCount,
    interval: DateInterval,
    ctx: &Context<'_>,
    warnings: &mut BTreeSet<Warning>,
) -> Result<f64> {
    let calendar = match basis {
        DayCount::Bus252 { calendar } => Some(ctx.refdata.calendar(
            &rhodium_core::calendar::CalendarSpec::Base(calendar.clone()),
        )?),
        _ => None,
    };
    let accrual = AccrualContext {
        calendar: calendar.as_ref(),
        maturity: Some(interval.end()),
        ..AccrualContext::default()
    };
    let (fraction, notices) = year_fraction(basis, interval.start(), interval.end(), &accrual)?;
    warnings.extend(notices);
    if fraction.value() <= 0.0 {
        return Err(Error::InvalidData(
            "futures requires a positive accrual fraction".into(),
        ));
    }
    Ok(fraction.value())
}

/// 1 本の率の読み方。
#[derive(Clone, Copy)]
struct RateRead<'a> {
    /// 率の増分に掛ける重み。
    weight: f64,
    /// 率を当てる観測の期間の全体の日数係数 τ。
    accrual: f64,
    /// 日数係数の基準。
    basis: &'a DayCount,
}

fn read_rate(
    index: &IndexId,
    date: Date,
    read: RateRead<'_>,
    legs: &mut Vec<ForecastLeg>,
    ctx: &Context<'_>,
    warnings: &mut BTreeSet<Warning>,
) -> Result<Option<f64>> {
    let point = ObservationPoint::Date(date);
    let definition = ctx.refdata.index(index)?;
    let event = definition
        .observation_time
        .as_ref()
        .ok_or_else(|| Error::InvalidData(format!("missing observation time for {index}")))?;
    let fixing = ctx
        .refdata
        .zone(&event.zone)?
        .resolve_convention(event.on_date(date))?;
    let (period, notices) = ctx.refdata.observation_period(index, point)?;
    warnings.extend(notices);
    let (publication, notices) = ctx.refdata.publication(index, point)?;
    warnings.extend(notices);
    let (resolution, warning) =
        ctx.fixings
            .resolve(index, point, ctx.as_of, publication, period)?;
    warnings.extend(warning);
    match resolution {
        FixingResolution::Published(rate) | FixingResolution::Estimated(rate) => Ok(Some(rate)),
        FixingResolution::Forecast => {
            legs.push(ForecastLeg {
                fixing: fixing.max(ctx.as_of),
                start: period.start(),
                end: period.end(),
                weight: read.weight,
                span: Span::Whole,
            });
            Ok(None)
        }
        // 評価時点が期間をまたぐ。[A, 終了) の forward の率を期間全体の率とみなす（RFR の
        // `insert_forward` と同じ規則）。残りの日数係数は観測時刻の間を時間で按分して求める。
        FixingResolution::ForecastFrom(at) => {
            let zone = ctx.refdata.zone(&event.zone)?;
            let end = zone.to_local(period.end())?.date();
            let remaining = super::rfr::remaining_fraction(at, end, zone, event, |from, to| {
                if from >= to {
                    Ok(0.0)
                } else {
                    fraction(read.basis, DateInterval::new(from, to)?, ctx, warnings)
                }
            })?;
            if remaining <= 0.0 {
                return Err(Error::Numerical(
                    "futures forecast requires a positive remaining fraction".into(),
                ));
            }
            legs.push(ForecastLeg {
                fixing: at,
                start: at,
                end: period.end(),
                weight: read.weight,
                span: Span::FromAsOf {
                    scale: read.accrual / remaining,
                },
            });
            Ok(None)
        }
    }
}

pub(super) fn clearing_delay(
    row: RateFutureClearing,
    currency: Currency,
    ctx: &Context<'_>,
    warnings: &mut BTreeSet<Warning>,
) -> Result<()> {
    if row.observed_at == row.settlement {
        return Ok(());
    }
    let terms = ctx.refdata.currency(currency)?;
    let zone = ctx.refdata.zone(&terms.settlement_time.zone)?;
    let start = zone.to_local(row.observed_at)?.date();
    let warning = super::rfr::payment_delay(ctx, currency, start, row.settlement, warnings)?;
    warnings.extend(warning);
    warnings.insert(
        ctx.config
            .conventions()
            .require(LinearConvention::PaymentDelayIgnored)?,
    );
    Ok(())
}
`````

`````
crates/rhodium-pricing/src/linear/fx.rs
//! FX リセットに使う担保別フォワードを準備する。

use std::{collections::BTreeSet, sync::Arc};

use rhodium_core::{
    Error, Result, Warning,
    calendar::{AdvanceRule, BusinessDayConvention, EndOfMonthRule, advance},
    fixing::ObservationPoint,
    time::Instant,
};
use rhodium_market::{fx::FxForwardPlan, structure::MarketStructure};
use rhodium_products::ir::analysis::FxFixing;

use crate::{config::LinearConvention, context::Context, resolve::CollateralKey};

pub(super) fn prepare(
    fixing: FxFixing<'_>,
    value: Instant,
    collateral: &CollateralKey,
    structure: &Arc<MarketStructure>,
    ctx: &Context<'_>,
    warnings: &mut BTreeSet<Warning>,
) -> Result<FxForwardPlan> {
    let pair = fixing.pair;
    let discounts = ctx
        .config
        .discount_curve(pair.base(), collateral)
        .zip(ctx.config.discount_curve(pair.quote(), collateral))
        .map(|(base, quote)| [base, quote]);
    let uncollateralized = matches!(collateral, CollateralKey::Uncollateralized);
    let preferred = if uncollateralized {
        Some(match ctx.config.fx_collateral(pair) {
            Some(id) => id,
            None => structure.fx_market_collateral(pair)?,
        })
    } else {
        None
    };
    let (id, same_collateral) = structure.fx_projection_collateral(pair, discounts, preferred)?;
    if !same_collateral {
        warnings.insert(
            ctx.config
                .conventions()
                .require(LinearConvention::CrossCollateralProjection)?,
        );
    }
    structure.fx_forward_plan(pair, id, value)
}

/// FX fixing の予想を読む受渡の瞬間。fixing 日からペアのスポットラグを両通貨の休日の和で
/// 進めた日の、ペアの慣行の受渡時刻とする。
///
/// 同じ指標・観測日の公表値は一つなので、元本リセットの fixing もこの瞬間で読む。新しい元本を
/// 受け渡す瞬間（[`FxFixing::value_at`]）は脚の通貨の決済時刻で、予想の受渡とは別に扱う。
pub(super) fn value_time(
    fixing: FxFixing<'_>,
    ctx: &Context<'_>,
    warnings: &mut BTreeSet<Warning>,
) -> Result<Instant> {
    let pair = fixing.pair;
    let convention = ctx.refdata.fx_convention(pair)?;
    let calendar = ctx.refdata.calendar(&ctx.refdata.fx_calendar_spec(pair)?)?;
    let event = &convention.settlement;
    let zone = ctx.refdata.zone(&event.zone)?;
    let date = match fixing.point {
        ObservationPoint::Date(date) => date,
        ObservationPoint::Instant(at) => zone.to_local(at)?.date(),
        ObservationPoint::Month(_) => return Err(Error::Internal("monthly FX fixing".into())),
    };
    let (value_date, notices) = advance(
        date,
        convention.spot_lag.offset(),
        &AdvanceRule {
            calendar: &calendar,
            bdc: BusinessDayConvention::Following,
            end_of_month: EndOfMonthRule::Ignore,
        },
    )?;
    warnings.extend(notices);
    zone.resolve_convention(event.on_date(value_date))
}
`````

`````
crates/rhodium-pricing/src/linear/rfr.rs
//! 翌日物予想と、端点の比に畳む条件。

use std::collections::{BTreeSet, HashMap};

use rhodium_core::{
    Error, Result, Warning,
    calendar::{Calendar, CalendarSpec},
    date::{Date, DateInterval},
    day_count::{AccrualContext, DayCount, year_fraction},
    fixing::{ObservationPeriod, ObservationPoint},
    id::IndexId,
    observation::{Aggregation, CompoundingForm, ObservationMethod, SpreadTreatment},
    quantity::YearFraction,
    refdata::{IndexKind, RefData},
    rounding::Rounding,
    time::{Instant, LocalEventTime, ZoneRules},
};
use rhodium_market::{
    curve::CurvePoint,
    view::{BoundSlots, Discount, MarketView, Projection, ResolvedCurve},
};
use rhodium_products::ir::{
    Product,
    analysis::{Compounded, RfrObservation, TermFixing},
};
use scioxide_ad::{Mask, Real};
use scioxide_core::Tolerance;

use crate::{config::LinearConvention, context::Context};

#[derive(Debug)]
pub(super) struct CompoundPlan {
    slot: usize,
    growth: Option<usize>,
    daily: Vec<DailyFactor>,
    form: CompoundingForm,
    spread: f64,
    rounding: Option<Rounding>,
    /// 簡約が畳んだ確定因子の積から 1 を引いた値（IR の `Compounded::fixed_excess`）。
    fixed_excess: f64,
}

#[derive(Debug)]
struct DailyFactor {
    point: usize,
    weight: f64,
}

impl CompoundPlan {
    /// 日次の因子を持たない成長のスロット、点、既知の成長の超過。
    pub(super) fn pure_growth(&self) -> Option<(usize, usize, f64)> {
        self.growth
            .filter(|_| self.daily.is_empty() && self.fixed_excess > -1.0)
            .map(|point| (self.slot, point, self.fixed_excess))
    }

    pub(super) fn evaluate<'v, R: Real>(
        &self,
        view: &MarketView<'v, R>,
        slots: &BoundSlots<'v, R>,
    ) -> Result<R> {
        let suffix = self.growth.map_or_else(R::one, |point| {
            view.growth(slots.projection_growth()[self.slot].points()[point])
        });
        // 確定因子の積 1 + e は θ によらない定数。
        let suffix = suffix * (1.0 + self.fixed_excess);
        if self.daily.is_empty() {
            return Ok(suffix);
        }

        let rates = slots.projections()[self.slot].points();
        let growth = self.daily.iter().try_fold(suffix, |growth, item| {
            let rate = view.forward(rates[item.point]);
            let rate = self
                .rounding
                .map_or(Ok(rate), |rounding| rounding.apply(rate))?;
            let factor = self.form.factor(rate, item.weight, self.spread)?;
            Ok(growth * factor)
        })?;
        if growth.is_finite().all() && growth.gt(R::zero()).all() {
            Ok(growth)
        } else {
            Err(Error::Numerical(
                "compounded growth is nonpositive or non-finite".into(),
            ))
        }
    }

    /// 超過成長 G − 1 を、予想する部分の対数の和と確定部分の超過から桁を落とさずに求める。
    ///
    /// G を丸めてから 1 を引くと、短い期間（1 日の期間で G − 1 ≈ 1.4e-5）の超過は
    /// ulp(1) ≈ 2.2e-16 の刻みになり、年率化したパーが 1e-14 の校正の許容誤差に届かない。
    /// 予想する部分は、DF 比に畳んだ部分の ln P(開始) − ln P(終了) と日次の因子の `log_factor` を足し、
    /// `exp_m1` で超過 f にする。確定部分の超過 e（簡約が残した `fixed_excess`）と合わせた
    /// (1 + e)(1 + f) − 1 = e + f + e·f の相対誤差は、倍精度の ε 程度に保たれる。
    pub(super) fn evaluate_excess<'v, R: Real>(
        &self,
        view: &MarketView<'v, R>,
        slots: &BoundSlots<'v, R>,
    ) -> Result<R> {
        let suffix = self.growth.map_or_else(R::zero, |point| {
            view.log_growth(slots.projection_growth()[self.slot].points()[point])
        });
        let rates = slots.projections()[self.slot].points();
        let log = self.daily.iter().try_fold(suffix, |log, item| {
            let rate = view.forward(rates[item.point]);
            let rate = self
                .rounding
                .map_or(Ok(rate), |rounding| rounding.apply(rate))?;
            Ok::<_, Error>(log + self.form.log_factor(rate, item.weight, self.spread)?)
        })?;
        let forecast = log.exp_m1();
        let excess = forecast + self.fixed_excess + forecast * self.fixed_excess;
        if excess.is_finite().all() {
            Ok(excess)
        } else {
            Err(Error::Numerical(
                "compounded growth is nonpositive or non-finite".into(),
            ))
        }
    }
}

/// ターム金利の観測を受け渡す時点。
#[derive(Debug, Clone, Copy)]
enum TermPayment {
    /// 観測期間の終わりに払うクーポン（`TermForward`）。
    PeriodEnd,
    /// 期首に割り引いて払う FRA の決済（`FraSettlement`）。
    FraStart,
}

pub(super) struct Forecast<'a> {
    pub(super) curve: ResolvedCurve<Projection>,
    pub(super) rates: Vec<(CurvePoint, CurvePoint, YearFraction)>,
    pub(super) growth: Vec<(CurvePoint, CurvePoint)>,
    slot: usize,
    rate_points: HashMap<ObservationPoint, usize>,
    growth_points: HashMap<(Instant, Instant), usize>,
    points: HashMap<Instant, CurvePoint>,
    day_count: &'a DayCount,
    observation_time: &'a LocalEventTime,
    zone: &'a ZoneRules,
    calendar: Option<Calendar>,
}

struct ProjectionRow {
    original: RfrObservation,
    dates: DateInterval,
    partial_from: Option<Instant>,
    tau: YearFraction,
}

impl<'a> Forecast<'a> {
    pub(super) fn new(
        index: &IndexId,
        curve: ResolvedCurve<Projection>,
        refdata: &'a RefData,
        slot: usize,
    ) -> Result<Self> {
        let definition = refdata.index(index)?;
        let (IndexKind::OvernightRate { day_count, .. } | IndexKind::TermRate { day_count, .. }) =
            &definition.kind
        else {
            return Err(Error::InvalidData(format!(
                "forward requires a dated rate: {index}"
            )));
        };
        let observation_time = definition
            .observation_time
            .as_ref()
            .ok_or_else(|| Error::InvalidData(format!("missing observation time for {index}")))?;
        let zone = refdata.zone(&observation_time.zone)?;
        let calendar = match day_count {
            DayCount::Bus252 { calendar } => {
                Some(refdata.calendar(&CalendarSpec::Base(calendar.clone()))?)
            }
            _ => None,
        };
        Ok(Self {
            curve,
            rates: Vec::new(),
            growth: Vec::new(),
            slot,
            rate_points: HashMap::new(),
            growth_points: HashMap::new(),
            points: HashMap::new(),
            day_count,
            observation_time,
            zone,
            calendar,
        })
    }

    pub(super) fn term_point(
        &mut self,
        fixing: TermFixing<'_>,
        discount: &ResolvedCurve<Discount>,
        pay: Instant,
        ctx: &Context<'_>,
        warnings: &mut BTreeSet<Warning>,
    ) -> Result<usize> {
        self.term_point_paid(fixing, discount, pay, ctx, warnings, TermPayment::PeriodEnd)
    }

    pub(super) fn fra_point(
        &mut self,
        fixing: TermFixing<'_>,
        discount: &ResolvedCurve<Discount>,
        pay: Instant,
        ctx: &Context<'_>,
        warnings: &mut BTreeSet<Warning>,
    ) -> Result<usize> {
        self.term_point_paid(fixing, discount, pay, ctx, warnings, TermPayment::FraStart)
    }

    fn term_point_paid(
        &mut self,
        fixing: TermFixing<'_>,
        discount: &ResolvedCurve<Discount>,
        pay: Instant,
        ctx: &Context<'_>,
        warnings: &mut BTreeSet<Warning>,
        payment: TermPayment,
    ) -> Result<usize> {
        if fixing.publication <= ctx.as_of || fixing.period.end() <= ctx.as_of {
            return Err(Error::MissingData(format!(
                "missing fixing {} at {:?}",
                fixing.index, fixing.point
            )));
        }
        if self.curve.spec().currency() != discount.spec().currency() {
            return Err(Error::Unsupported(
                "term projection in another payment currency".into(),
            ));
        }
        let conventions = ctx.config.conventions();
        // 予想の担保（予想の測度）と支払いの担保が違うときだけ、担保の違いによる予想の調整を無視する。
        if !self.curve.is_forecast_under(discount) {
            warnings.insert(conventions.require(LinearConvention::CrossCollateralProjection)?);
        }
        let natural = match payment {
            TermPayment::PeriodEnd => fixing.period.end(),
            // FRA の決済 N·(1 − (1 + τK)·P_p(e)/P_p(s))·P_d(s) は P_p = P_d のときだけ厳密で、
            // 同じ担保でも予想曲線と割引曲線が違えば凸性を無視する。
            TermPayment::FraStart => {
                if self.curve.spec().id() != discount.spec().id() {
                    warnings.insert(conventions.require(LinearConvention::FraDiscounting)?);
                }
                fixing.period.start()
            }
        };
        let date = self.zone.to_local(natural)?.date();
        let delay = payment_delay(ctx, discount.spec().currency(), date, pay, warnings)?;
        warnings.extend(delay);
        self.forward_point(fixing.point, fixing.period, fixing.forecast_from, warnings)
    }

    pub(super) fn plan(
        &mut self,
        product: &Product,
        compound: Compounded<'_>,
        discount: &ResolvedCurve<Discount>,
        pay: Instant,
        ctx: &Context<'_>,
        warnings: &mut BTreeSet<Warning>,
    ) -> Result<CompoundPlan> {
        let Aggregation::Compound { form, spread } = compound.rule.aggregation() else {
            return Err(Error::Internal(
                "expected a compound observation rule".into(),
            ));
        };
        if self.curve.spec().currency() != discount.spec().currency() {
            return Err(Error::Unsupported(
                "RFR projection in another payment currency".into(),
            ));
        }

        let rows = compound
            .rows()
            .map(|row| self.row(product, compound.index, row, ctx.as_of, warnings))
            .collect::<Result<Vec<_>>>()?;
        let rounding = compound.rule.rounding().daily;
        let mut plan = CompoundPlan {
            slot: self.slot,
            growth: None,
            daily: Vec::new(),
            form,
            spread: compound.daily_spread,
            rounding,
            fixed_excess: compound.fixed_excess,
        };
        let Some(first) = rows.first() else {
            return Ok(plan);
        };
        let conventions = ctx.config.conventions();
        // 予想の担保（予想の測度）と支払いの担保が違うときだけ、担保の違いによる予想の調整を無視する。
        if !self.curve.is_forecast_under(discount) {
            warnings.insert(conventions.require(LinearConvention::CrossCollateralProjection)?);
        }
        if rounding.is_some() {
            warnings.insert(conventions.require(LinearConvention::RoundingAtForward)?);
        }
        if let Some(last) = compound.last_period() {
            let natural = self.zone.to_local(last.end())?.date();
            if let Some(warning) =
                payment_delay(ctx, discount.spec().currency(), natural, pay, warnings)?
            {
                warnings.insert(warning);
            }
        }

        // 設計書「RFR の telescoping」の 5 条件を確かめる。
        // 条件 4（予想の率を予想曲線の DF 比で定義する）は Forecast の forward の定義で常に満たす。
        let can_telescope =
            // 条件 1: 日次丸めのない SimpleFactor で、スプレッドは Exclusive。各因子が 1 + w·r になる。
            form == CompoundingForm::SimpleFactor
            && spread == SpreadTreatment::Exclusive
            && rounding.is_none()
            // 条件 2 の前提: lookback と lockout は金利の期間と重みの期間がずれる。
            && matches!(
                compound.rule.method(),
                ObservationMethod::Plain | ObservationMethod::ObservationShift(_)
            )
            // 条件 3: 確定した因子は元の列の先頭に連続していた。
            && compound.has_fixed_prefix()
            // 条件 5: 重みの日数基準が予想の日数基準と一致する。
            && compound.basis() == self.day_count
            // 条件 2 と 3: 期間が隙間なく続き、重みが金利の期間の日数係数に等しく、
            // 評価時点をまたぐ要素は残りの列の先頭にだけある。
            && is_contiguous(&rows)?;
        let daily = if can_telescope {
            // またぐ要素は厳密な因子 Π_A = 1 + w_A·r_A として日次で残し、その次の要素から端点の比に畳む。
            let n = usize::from(first.partial_from.is_some());
            let suffix = &rows[n..];
            if let Some((start, end)) = suffix.first().zip(suffix.last()) {
                plan.growth = Some(
                    self.growth_point(start.original.period.start(), end.original.period.end())?,
                );
            }
            &rows[..n]
        } else {
            if matches!(compound.rule.method(), ObservationMethod::Lookback(lag) | ObservationMethod::Lockout(lag)
                if lag.get() > conventions.max_lookback_lockout_business_days())
            {
                return Err(Error::Unsupported(
                    "lookback or lockout exceeds the daily forward product limit".into(),
                ));
            }
            warnings.insert(conventions.require(LinearConvention::DailyForwardProduct)?);
            &rows[..]
        };
        plan.daily = daily
            .iter()
            .map(|row| {
                Ok(DailyFactor {
                    point: self.rate_point(row, warnings)?,
                    weight: row.original.weight,
                })
            })
            .collect::<Result<_>>()?;
        Ok(plan)
    }

    fn row(
        &self,
        product: &Product,
        index: &IndexId,
        original: RfrObservation,
        as_of: Instant,
        warnings: &mut BTreeSet<Warning>,
    ) -> Result<ProjectionRow> {
        if original.publication <= as_of || original.period.end() <= as_of {
            return Err(Error::MissingData(format!(
                "missing fixing {index} at {:?}",
                original.point
            )));
        }
        let dates = DateInterval::new(
            self.zone.to_local(original.period.start())?.date(),
            self.zone.to_local(original.period.end())?.date(),
        )?;
        let tau = self.fraction(dates.start(), dates.end(), warnings)?;
        Ok(ProjectionRow {
            original,
            dates,
            partial_from: product.forecast_from(index, original.point),
            tau,
        })
    }

    fn fraction(
        &self,
        start: Date,
        end: Date,
        warnings: &mut BTreeSet<Warning>,
    ) -> Result<YearFraction> {
        let context = AccrualContext {
            calendar: self.calendar.as_ref(),
            ..AccrualContext::default()
        };
        let (fraction, notices) = year_fraction(self.day_count, start, end, &context)?;
        warnings.extend(notices);
        Ok(fraction)
    }

    fn remaining_fraction(
        &self,
        at: Instant,
        end: Date,
        warnings: &mut BTreeSet<Warning>,
    ) -> Result<YearFraction> {
        let value = remaining_fraction(at, end, self.zone, self.observation_time, |from, to| {
            Ok(self.fraction(from, to, warnings)?.value())
        })?;
        YearFraction::new(value, self.day_count.clone())
    }

    fn rate_point(
        &mut self,
        row: &ProjectionRow,
        warnings: &mut BTreeSet<Warning>,
    ) -> Result<usize> {
        self.insert_forward(
            row.original.point,
            row.original.period,
            row.partial_from,
            row.tau.clone(),
            warnings,
        )
    }

    fn forward_point(
        &mut self,
        point: ObservationPoint,
        period: ObservationPeriod,
        partial: Option<Instant>,
        warnings: &mut BTreeSet<Warning>,
    ) -> Result<usize> {
        if let Some(&pos) = self.rate_points.get(&point) {
            return Ok(pos);
        }
        let start = self.zone.to_local(period.start())?.date();
        let end = self.zone.to_local(period.end())?.date();
        let tau = self.fraction(start, end, warnings)?;
        self.insert_forward(point, period, partial, tau, warnings)
    }

    fn insert_forward(
        &mut self,
        point: ObservationPoint,
        period: ObservationPeriod,
        partial: Option<Instant>,
        full: YearFraction,
        warnings: &mut BTreeSet<Warning>,
    ) -> Result<usize> {
        if let Some(&pos) = self.rate_points.get(&point) {
            return Ok(pos);
        }
        let (start, tau) = match partial {
            Some(at) => (
                at,
                self.remaining_fraction(at, self.zone.to_local(period.end())?.date(), warnings)?,
            ),
            None => (period.start(), full),
        };
        if tau.value() <= 0.0 {
            return Err(Error::Numerical(
                "forward requires a positive year fraction".into(),
            ));
        }
        let start = self.curve_point(start)?;
        let end = self.curve_point(period.end())?;
        let pos = self.rates.len();
        self.rates.push((start, end, tau));
        self.rate_points.insert(point, pos);
        Ok(pos)
    }

    fn growth_point(&mut self, start: Instant, end: Instant) -> Result<usize> {
        if let Some(&pos) = self.growth_points.get(&(start, end)) {
            return Ok(pos);
        }
        let begin = self.curve_point(start)?;
        let finish = self.curve_point(end)?;
        let pos = self.growth.len();
        self.growth.push((begin, finish));
        self.growth_points.insert((start, end), pos);
        Ok(pos)
    }

    fn curve_point(&mut self, at: Instant) -> Result<CurvePoint> {
        use std::collections::hash_map::Entry;
        match self.points.entry(at) {
            Entry::Occupied(entry) => Ok(entry.get().clone()),
            Entry::Vacant(entry) => Ok(entry.insert(self.curve.prepare_point(at)?).clone()),
        }
    }
}

fn is_contiguous(rows: &[ProjectionRow]) -> Result<bool> {
    // 期間と日数基準は構造として照合済みなので、重みと τ の差は演算の順序の誤差だけになる。
    // P8 に従いビット単位の一致は求めず、品質の表の market の値の許容誤差（ε = 1e-12、a = 1e-14）で比べる。
    let tolerance = Tolerance::new(1e-14, 1e-12).map_err(numerical)?;
    let mut previous_end = None;
    for row in rows {
        // 同じ重みが偶然得られても、金利の期間と重みの期間が違えば畳まない。
        let follows = previous_end
            .is_none_or(|end| end == row.original.period.start() && row.partial_from.is_none());
        previous_end = Some(row.original.period.end());
        let same_weight = row.tau.value() > 0.0
            && tolerance
                .is_converged(row.original.weight, row.tau.value())
                .map_err(numerical)?;
        if !follows || row.original.coverage != row.dates || !same_weight {
            return Ok(false);
        }
    }
    Ok(true)
}

fn numerical(error: impl std::fmt::Display) -> Error {
    Error::Numerical(error.to_string())
}

/// 評価時点 `at` から観測の期間の終了日 `end` までの残りの日数係数。
///
/// `at` を含む暦日の観測時刻から翌日の観測時刻までを時間で按分し、その前後の暦日から `end` までの
/// 日数係数を補間する。観測時刻から翌日の観測時刻までの長さを使うので、夏時間の開始・終了日
/// （24 時間と限らない日）でも按分の比は 0 から 1 に収まる。`fraction(start, end)` は日付の間の
/// 日数係数で、`start == end` のときは 0 を返すこと。
pub(super) fn remaining_fraction(
    at: Instant,
    end: Date,
    zone: &ZoneRules,
    observation_time: &LocalEventTime,
    mut fraction: impl FnMut(Date, Date) -> Result<f64>,
) -> Result<f64> {
    let day = zone.to_local(at)?.date();
    let cutoff = |date| zone.resolve_convention(observation_time.on_date(date));
    let start = if at < cutoff(day)? {
        day.checked_add_days(-1)?
    } else {
        day
    };
    let next = start.checked_add_days(1)?;
    let from = cutoff(start)?.unix_micros();
    let to = cutoff(next)?.unix_micros();
    let remaining = (to - at.unix_micros()) as f64 / (to - from) as f64;
    let first = fraction(start, end)?;
    let second = fraction(next, end)?;
    Ok(second + remaining * (first - second))
}

pub(super) fn payment_delay(
    ctx: &Context<'_>,
    currency: rhodium_core::currency::Currency,
    natural: Date,
    pay: Instant,
    warnings: &mut BTreeSet<Warning>,
) -> Result<Option<Warning>> {
    let definition = ctx.refdata.currency(currency)?;
    let zone = ctx.refdata.zone(&definition.settlement_time.zone)?;
    let pay = zone.to_local(pay)?.date();
    if natural == pay {
        return Ok(None);
    }
    let calendar = ctx
        .refdata
        .calendar(&CalendarSpec::Base(definition.settlement_calendar.clone()))?;
    let (delay, notices) = calendar.business_days_between(natural.min(pay), natural.max(pay))?;
    warnings.extend(notices);
    let conventions = ctx.config.conventions();
    if delay > conventions.max_payment_delay_business_days() {
        return Err(Error::Unsupported(format!(
            "payment delay of {delay} business days exceeds the linear limit"
        )));
    }
    conventions
        .require(LinearConvention::PaymentDelayIgnored)
        .map(Some)
}
`````

`````
crates/rhodium-pricing/src/par.rs
//! 一つの商品の PV をゼロにする値の計算。

use rhodium_core::{Error, Result};
use scioxide_ad::{BaseReal, Dual, Mask, Real, Residual, implicit_solution_into};
use scioxide_core::{ResidualTolerance, SciOxideError};
use scioxide_rootfind::{
    SafeguardedInitialGuess, SafeguardedNewton, ScalarFn, ScalarRootOptions, ScalarRootSolver,
};

pub(crate) fn affine<R: Real>(pv_zero: R, pv_one: R) -> Result<R> {
    let slope = pv_one - pv_zero;
    if !slope.is_finite().all() || !slope.abs().gt(R::zero()).all() {
        return Err(Error::Calibration(
            "par value requires a finite nonzero PV slope".into(),
        ));
    }
    let par = -pv_zero / slope;
    if !par.is_finite().all() {
        return Err(Error::Calibration("non-finite par value".into()));
    }
    Ok(par)
}

/// 一つの未知数を持つ非線形の PV 残差を解き、パーを求める。
///
/// 残差は第一引数の長さ 1 の配列からパーの候補を読み、
/// 第二引数の市場パラメータに対する PV を出力の一要素に書く。
/// 両端が根を挟む区間とその内側の初期値を渡す。
/// PV の絶対許容誤差は 1e-14。入れ子の AD 型では収束点の陰関数の微分を返す。
///
/// # Errors
///
/// 残差評価、求根、陰関数の復元が失敗するか、根での PV の傾きがゼロの場合は
/// [`Error::Calibration`]。
pub fn solve_nonlinear<R: Real<f64, Mask = bool>>(
    pv: &impl Residual,
    parameters: &[R],
    initial: SafeguardedInitialGuess<f64>,
) -> Result<R> {
    let tolerance = ResidualTolerance::new(1e-14, 0.0).map_err(par_error)?;
    let base: Vec<_> = parameters
        .iter()
        .map(|value| value.base().lane(0).unwrap_or(f64::NAN))
        .collect();
    let dual: Vec<_> = base.iter().copied().map(Dual::lift).collect();
    let derivative = |rate| {
        let mut out = [Dual::lift(0.0)];
        pv.residual(&[Dual::new(rate, 1.0)], &dual, &mut out)?;
        Ok(out[0].derivative())
    };
    let function = ScalarFn::new(|rate| {
        let mut out = [0.0];
        pv.residual(&[rate], &base, &mut out)?;
        Ok(out[0])
    })
    .with_derivative(derivative);
    let solver = SafeguardedNewton::default();
    let options = ScalarRootOptions::default().with_residual_tolerance(tolerance);
    let root = solver
        .solve(&function, initial, &options)
        .map_err(par_error)?
        .point();
    let slope = derivative(root).map_err(par_error)?;
    if !slope.is_finite() || slope == 0.0 {
        return Err(Error::Calibration(
            "nonlinear par value has a singular PV slope".into(),
        ));
    }
    // 収束点を指定し、入れ子の AD の次数まで陰関数の微分を復元する。
    let mut out = [R::zero()];
    implicit_solution_into(pv, &[root], parameters, tolerance, &mut out).map_err(par_error)?;
    Ok(out[0])
}

#[expect(clippy::needless_pass_by_value, reason = "map_err の変換関数に使う")]
fn par_error(error: SciOxideError) -> Error {
    Error::Calibration(format!("nonlinear par: {error}"))
}
`````

`````
crates/rhodium-pricing/src/prepare.rs
//! 市場構造に対する準備と、特定の市場での評価。

use std::sync::Arc;

use rhodium_core::{Error, Result, Warning, date::Date, money::MoneyBag, time::Instant};
use rhodium_market::{
    market::Market,
    structure::{BlockId, MarketStructure},
    view::{BoundSlots, MarketView, Slots},
};
use rhodium_products::{
    agreement::MarginSet,
    ir::{
        Product,
        analysis::{Linearity, linearity},
        reduce,
    },
    standard::{ProductKind, StandardProduct},
    trade::{ExerciseRecord, ProductSpec, Trade},
};
use scioxide_ad::Real;

use crate::{
    config::{AnalyticPricer, Method},
    context::Context,
    linear::{self, CashflowPlan},
    resolve::CollateralKey,
};

/// 準備する商品。方式の表は商品の種類（標準商品の `ProductKind`）で引く。
#[derive(Debug, Clone, Copy)]
pub enum PrepareInput<'a> {
    /// 標準商品と、日程の解決に使う約定日。
    Standard {
        /// 商品の条件。
        product: &'a StandardProduct,
        /// 日程を生成する起点の約定日。
        trade_date: Date,
    },
    /// 簡約前の検証済み IR。構築時と同じ参照データを評価にも使う。種類を持たない。
    Ir(&'a Product),
    /// 取引と、同じ参照データでその取引を下ろした検証済み IR（`Trade::to_product` か
    /// `TradeProductResolver::resolve` の結果）。
    ///
    /// 下ろした IR を作り直さずに使い、方式の表は取引の商品の種類で引く。解析プライサーは
    /// 取引の標準商品、約定日と未決済の手数料を読む。
    Trade {
        /// 取引の条件。
        trade: &'a Trade,
        /// `trade` を下ろした IR。
        product: &'a Product,
    },
}

impl PrepareInput<'_> {
    /// 方式の表が照合する商品の種類。IR だけを渡した場合と、変更した取引は持たない。
    fn kind(&self) -> Option<ProductKind> {
        match self {
            Self::Standard { product, .. } => Some(product.kind()),
            Self::Ir(_) => None,
            Self::Trade { trade, .. } => trade.product_kind(),
        }
    }
}

/// 市場のパラメータ θ に依存しない、評価に必要な計画。
#[derive(Debug)]
pub struct Prepared {
    as_of: Instant,
    method: Method,
    slots: Slots,
    plan: CashflowPlan,
    blocks: Vec<BlockId>,
    warnings: Vec<Warning>,
}

impl Prepared {
    pub(crate) fn ois_warm_start(&self, target: &BlockId) -> Option<linear::OisWarmStart> {
        self.plan.ois_warm_start(&self.slots, target)
    }

    /// fixing と決済を簡約した評価の瞬間。
    #[must_use]
    pub const fn as_of(&self) -> Instant {
        self.as_of
    }

    /// 使用する評価方式。
    #[must_use]
    pub const fn method(&self) -> Method {
        self.method
    }

    /// 評価が読む市場パラメータのブロック。
    #[must_use]
    pub fn blocks(&self) -> &[BlockId] {
        &self.blocks
    }

    /// 簡約と計画作成で生じた警告。
    #[must_use]
    pub fn warnings(&self) -> &[Warning] {
        &self.warnings
    }

    /// 計画の曲線と点を、評価する市場ビューに結び付ける。
    ///
    /// # Errors
    /// 構造が異なれば `InvalidInput`、必要なブロックが使えなければ
    /// `MissingData`、曲線の値を表現できなければ `Numerical`。
    pub fn bind<'p, 'v, R: Real>(
        &'p self,
        view: &'p MarketView<'v, R>,
    ) -> Result<BoundPrepared<'p, 'v, R>> {
        Ok(BoundPrepared {
            plan: &self.plan,
            slots: view.bind_slots(&self.slots)?,
            view,
        })
    }
}

/// 一つの市場ビューに結び付けた評価計画。
#[derive(Debug)]
pub struct BoundPrepared<'p, 'v, R = f64> {
    plan: &'p CashflowPlan,
    slots: BoundSlots<'v, R>,
    view: &'p MarketView<'v, R>,
}

impl<R: Real> BoundPrepared<'_, '_, R> {
    /// 結び付けた値から通貨ごとの PV を求める。
    ///
    /// # Errors
    /// 複利因子、丸めまたは PV を表現できなければ `Numerical`。
    pub fn evaluate(&self) -> Result<MoneyBag<R>> {
        self.plan.evaluate(self.view, &self.slots)
    }
}

/// fixing・行使・決済を簡約し、市場の位置と補間点を準備する。
///
/// # Errors
/// 構造の評価時点、入力の参照データ、日程、行使記録または設定が不正か、必要な
/// fixing や曲線が無い場合。リニア方式で扱えない商品は `Unsupported`。
pub fn prepare(
    input: PrepareInput<'_>,
    collateral: &CollateralKey,
    exercises: &[ExerciseRecord],
    structure: &Arc<MarketStructure>,
    ctx: &Context<'_>,
) -> Result<Prepared> {
    if structure.as_of() != ctx.as_of {
        return Err(Error::InvalidInput(
            "market structure and context have different valuation instants".into(),
        ));
    }
    let kind = input.kind();
    // 方式の表で解析プライサーを選んだ標準商品は、IR を簡約せずに型付きの商品から準備する。
    // 表の照合はリニア評価の可否によらないので、ここでは可否を問わずに引く。
    if kind.is_some()
        && let Method::Analytic(pricer) = ctx.config.select_method(kind, Linearity::Linear)?
    {
        return prepare_analytic(pricer, input, collateral, structure, ctx);
    }
    let lowered;
    let product = match input {
        PrepareInput::Standard {
            product,
            trade_date,
        } => {
            lowered = product.to_product(trade_date, ctx.refdata)?;
            &lowered
        }
        PrepareInput::Ir(product) | PrepareInput::Trade { product, .. } => product,
    };
    let simplified = reduce(
        product,
        ctx.as_of,
        ctx.fixings,
        exercises,
        ctx.config.settlement(),
    )?;
    let applicable = linearity(&simplified);
    let method = ctx
        .config
        .select_method(kind, applicable)
        .map_err(|error| match applicable {
            Linearity::NonLinear(reason) => error.context(reason),
            Linearity::Linear => error,
        })?;
    // 方式を足したときに、その計画を作らないままリニアで評価しないよう網羅して分岐する。
    let linear = match method {
        Method::Linear => linear::prepare(&simplified, collateral, structure, ctx)?,
        Method::Analytic(_) => {
            return Err(Error::Unsupported(
                "analytic evaluation requires a supported standard product".into(),
            ));
        }
    };
    Ok(Prepared {
        as_of: ctx.as_of,
        method,
        slots: linear.slots,
        plan: linear.plan,
        blocks: linear.blocks,
        warnings: linear.warnings,
    })
}

/// 解析プライサーの入力を型付きの標準商品から作る。プライサーを網羅して分岐し、
/// 足したプライサーの入力を決めないままリニアの計画に落とさない。
///
/// 取引の入力では、約定日と未決済の手数料も読む（IR の評価と同じ受渡を数える）。
fn prepare_analytic(
    pricer: AnalyticPricer,
    input: PrepareInput<'_>,
    collateral: &CollateralKey,
    structure: &Arc<MarketStructure>,
    ctx: &Context<'_>,
) -> Result<Prepared> {
    match pricer {
        AnalyticPricer::RateFuture => {
            let (future, trade_date, fees) = match input {
                PrepareInput::Standard {
                    product: StandardProduct::RateFuture(future),
                    trade_date,
                } => (future, trade_date, &[][..]),
                PrepareInput::Trade { trade, .. } => match &trade.product {
                    ProductSpec::Standard(StandardProduct::RateFuture(future)) => {
                        (future, trade.trade_date, trade.fees.as_slice())
                    }
                    _ => return Err(rate_future_input()),
                },
                PrepareInput::Standard { .. } | PrepareInput::Ir(_) => {
                    return Err(rate_future_input());
                }
            };
            let (margin, warnings) = future.margin(trade_date, ctx.refdata)?;
            let plan =
                linear::prepare_rate_future(&margin, fees, warnings, collateral, structure, ctx)?;
            Ok(Prepared {
                as_of: ctx.as_of,
                method: Method::Analytic(pricer),
                slots: plan.slots,
                plan: plan.plan,
                blocks: plan.blocks,
                warnings: plan.warnings,
            })
        }
    }
}

fn rate_future_input() -> Error {
    Error::Unsupported("the rate-future analytic pricer requires a standard rate future".into())
}

/// マージンセットの VM 条件から担保の鍵を導いて、価格評価を準備する。
///
/// 支払いごとの割引規則は、支払い通貨と導出した担保の鍵で照合する。
///
/// # Errors
///
/// 担保の鍵を導けない場合は `Unsupported`。準備に必要なデータなどのエラーは
/// [`prepare`] と同じ。
pub fn prepare_with_margin_set(
    input: PrepareInput<'_>,
    margin_set: &MarginSet,
    exercises: &[ExerciseRecord],
    structure: &Arc<MarketStructure>,
    ctx: &Context<'_>,
) -> Result<Prepared> {
    let collateral = CollateralKey::from_margin_set(margin_set)?;
    prepare(input, &collateral, exercises, structure, ctx)
}

/// 準備を現在の市場に結び付けて、通貨ごとの PV を評価する。
///
/// # Errors
/// 構造が異なる、必要なブロックが無い、または数値計算が失敗した場合。
pub fn evaluate<R: Real>(prepared: &Prepared, market: &Market<R>) -> Result<MoneyBag<R>> {
    market.with_view(|view| prepared.bind(&view)?.evaluate())
}
`````

`````
crates/rhodium-pricing/src/resolve.rs
//! 価格設定で市場の要求を構造の位置に解決する。

use std::{cmp::Ordering, sync::Arc};

use rhodium_core::{
    Error, Result,
    currency::Currency,
    day_count::DayCount,
    id::IndexId,
    quantity::{Rate, Simple, Spread},
    refdata::{IndexKind, RefData},
};
use rhodium_market::{
    curve::CurveSpec,
    structure::MarketStructure,
    view::{Discount, Projection, ResolvedCurve},
};
use rhodium_products::agreement::{
    CashCollateralTerms, CollateralAsset, Eligibility, MarginSet, PerParty, Remuneration,
};
use serde::{Deserialize, Deserializer, de};

use crate::config::{FundingTreatment, PricingConfig};

/// 支払いの割引規則を照合するための担保条件。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum CollateralKey {
    /// 現金担保の付利条件。
    Cash(#[serde(deserialize_with = "deserialize_cash")] CashCollateralTerms),
    /// 最割安担保（CTD）の候補。
    CheapestToDeliver(CheapestToDeliverTerms),
    /// 担保なし。
    Uncollateralized,
}

impl CollateralKey {
    /// マージンセットの VM 条件から、割引に使う担保の鍵を導く。
    ///
    /// # Errors
    ///
    /// VM の CSA があり、両当事者が差し入れられる現金担保が無ければ `Unsupported`。
    pub fn from_margin_set(margin_set: &MarginSet) -> Result<Self> {
        let Some(vm) = &margin_set.vm else {
            return Ok(Self::Uncollateralized);
        };
        // 資産の種類を網羅して分岐し、担保資産を足したときに割引での扱いを決めずに
        // 黙って候補から外さないようにする。
        let candidates = vm
            .eligible()
            .iter()
            .filter_map(|entry| match (entry.asset(), entry.posters()) {
                (
                    CollateralAsset::Cash(terms),
                    PerParty {
                        own: Eligibility::Eligible,
                        cpty: Eligibility::Eligible,
                    },
                ) => Some(terms.clone()),
                // 片側だけが差し入れられる現金は割引の候補にしない（担保の選択の価値は
                // M12 で扱う）。証券担保の割引は M11 で CollateralKey::Securities として扱う。
                (CollateralAsset::Cash(_) | CollateralAsset::Security(_), _) => None,
            })
            .collect();
        let candidates = normalize_candidates(candidates);
        match candidates.as_slice() {
            [] => Err(Error::Unsupported(format!(
                "no cash collateral eligible for both parties in margin set {}",
                margin_set.id
            ))),
            [terms] => Ok(Self::Cash(terms.clone())),
            _ => Ok(Self::CheapestToDeliver(CheapestToDeliverTerms(
                candidates.into(),
            ))),
        }
    }

    /// 現金担保の全候補について、付利に必要な参照先を照合する。
    ///
    /// # Errors
    ///
    /// 通貨、指標またはフロアに必要なカレンダーが無ければ `InvalidData`。
    pub fn validate_references(&self, refdata: &RefData) -> Result<()> {
        match self {
            Self::Cash(terms) => terms.validate_references(refdata),
            Self::CheapestToDeliver(terms) => terms
                .candidates()
                .iter()
                .try_for_each(|terms| terms.validate_references(refdata)),
            Self::Uncollateralized => Ok(()),
        }
    }
}

/// 付利条件まで含めて順序をそろえた、2 種類以上の現金担保の集合。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(try_from = "Vec<CashInput>")]
pub struct CheapestToDeliverTerms(Arc<[CashCollateralTerms]>);

impl CheapestToDeliverTerms {
    /// 候補の重複を除き、CTD の条件を作る。
    ///
    /// # Errors
    ///
    /// 異なる候補が 2 種類未満なら `InvalidData`。
    pub fn new(candidates: Vec<CashCollateralTerms>) -> Result<Self> {
        let candidates = normalize_candidates(candidates);
        if candidates.len() < 2 {
            return Err(Error::InvalidData(
                "cheapest_to_deliver requires at least two distinct cash collateral terms".into(),
            ));
        }
        Ok(Self(candidates.into()))
    }

    /// 照合に使う順序に並べた候補。
    #[must_use]
    pub fn candidates(&self) -> &[CashCollateralTerms] {
        &self.0
    }
}

// CTD の候補を付利条件の全項目で並べ、重複を除く。設定から読んだ鍵とマージンセットから
// 導いた鍵が、候補の順序や同じ条件の繰り返し（ヘアカット違いなど）によらず `==` で
// 一致するようにする。並べ替えのキーは `PartialEq` と同じ項目を使い、等しい候補を
// 隣り合わせて `dedup` で除けるようにする。
fn normalize_candidates(mut candidates: Vec<CashCollateralTerms>) -> Vec<CashCollateralTerms> {
    candidates.sort_by(|a, b| {
        let (left, right) = (&a.remuneration, &b.remuneration);
        a.ccy
            .cmp(&b.ccy)
            .then_with(|| left.index.cmp(&right.index))
            .then_with(|| compare_rates(left.spread.value(), right.spread.value()))
            .then_with(|| match (&left.floor, &right.floor) {
                (Some(left), Some(right)) => compare_rates(left.value(), right.value())
                    .then_with(|| left.basis().to_string().cmp(&right.basis().to_string())),
                (Some(_), None) => Ordering::Greater,
                (None, Some(_)) => Ordering::Less,
                (None, None) => Ordering::Equal,
            })
    });
    candidates.dedup();
    candidates
}

fn compare_rates(left: f64, right: f64) -> Ordering {
    // 有限な率だけを比べる。±0 の並び順を PartialEq とそろえる。
    (left + 0.0).total_cmp(&(right + 0.0))
}

// 価格設定の現金担保の書式。`Remuneration` の入れ子にせず、付利の項目を通貨と同じ
// 階層に平らに書く。率は bp か % の文字列で受け取り、`into_terms` で数値にする。
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CashInput {
    ccy: Currency,
    index: IndexId,
    spread: String,
    floor: Option<FloorInput>,
}

// 付利のフロア。単利の率なので、率と日数基準を組で書かせる。
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FloorInput {
    rate: String,
    day_count: DayCount,
}

impl CashInput {
    fn into_terms(self) -> Result<CashCollateralTerms> {
        let spread = parse_collateral_rate(&self.spread)
            .and_then(Spread::new)
            .map_err(|error| error.context("cash collateral spread"))?;
        let floor = self
            .floor
            .map(|floor| Rate::new(parse_collateral_rate(&floor.rate)?, floor.day_count, Simple))
            .transpose()
            .map_err(|error| error.context("cash collateral floor"))?;
        Ok(CashCollateralTerms {
            ccy: self.ccy,
            remuneration: Remuneration {
                index: self.index,
                spread,
                floor,
            },
        })
    }
}

impl TryFrom<Vec<CashInput>> for CheapestToDeliverTerms {
    type Error = Error;

    fn try_from(input: Vec<CashInput>) -> Result<Self> {
        let candidates = input
            .into_iter()
            .map(CashInput::into_terms)
            .collect::<Result<_>>()?;
        Self::new(candidates)
    }
}

fn deserialize_cash<'de, D: Deserializer<'de>>(
    deserializer: D,
) -> std::result::Result<CashCollateralTerms, D::Error> {
    CashInput::deserialize(deserializer)?
        .into_terms()
        .map_err(de::Error::custom)
}

fn parse_collateral_rate(text: &str) -> Result<f64> {
    let (number, scale) = text
        .strip_suffix("bp")
        .map(|number| (number, 4))
        .or_else(|| text.strip_suffix('%').map(|number| (number, 2)))
        .ok_or_else(|| Error::InvalidData("cash collateral rate requires bp or %".into()))?;
    // 十進のまま尺度を変え、単位の違う等価な値を同じ f64 に丸める。
    format!("{number}e-{scale}")
        .parse()
        .map_err(|_| Error::InvalidData(format!("invalid cash collateral rate {text}")))
}

/// 商品の評価に必要な市場の要求。
#[derive(Debug, Clone, PartialEq)]
pub enum Requirement {
    /// 支払いの割引。
    Discounting {
        /// 支払いの通貨。
        ccy: Currency,
        /// 支払いに適用する担保条件。
        collateral: CollateralKey,
    },
    /// 指標の将来の fixing の予想。
    Forecast {
        /// 予想する指標。
        index: IndexId,
    },
}

/// 市場構造の位置に解決した要求。曲線の点の準備とスロットの作成に使う。
#[derive(Debug)]
pub enum ResolvedRequirement {
    /// 割引曲線。
    Discounting(ResolvedCurve<Discount>),
    /// 予想曲線。
    Forecast(ResolvedCurve<Projection>),
}

/// 市場の要求に対応する曲線を選び、構造の位置まで解決する。
///
/// # Errors
/// 割引・予想の規則、必要な資金調達曲線の指定、または指定先の曲線が無ければ
/// `MissingData`。曲線の役割や通貨、無担保の資金調達の指定が矛盾すれば
/// `InvalidInput`。指標定義が無ければ `InvalidData`、予想が未対応の指標なら
/// `Unsupported`。
pub fn resolve(
    requirement: &Requirement,
    config: &PricingConfig,
    structure: &Arc<MarketStructure>,
    refdata: &RefData,
) -> Result<ResolvedRequirement> {
    match requirement {
        Requirement::Discounting { ccy, collateral } => {
            let id = config.discount_curve(*ccy, collateral).ok_or_else(|| {
                Error::MissingData(format!(
                    "no discounting rule for {ccy} under {collateral:?}"
                ))
            })?;
            // InDiscounting では無担保の支払いを自社の資金調達曲線で割り引く。XVA の
            // FCA/FBA と同じ曲線を読むよう、割引行が借入・運用と同じ ID を指すことを求める。
            // 借入と運用が違う曲線では、一つの割引曲線で資金調達を表せないので拒否する。
            if config.funding() == FundingTreatment::InDiscounting
                && matches!(collateral, CollateralKey::Uncollateralized)
            {
                let funding = config
                    .funding_curves()
                    .get(ccy)
                    .ok_or_else(|| Error::MissingData(format!("no funding curves for {ccy}")))?;
                if funding.borrow != funding.lend || funding.borrow != *id {
                    return Err(Error::InvalidInput(format!(
                        "{ccy} uncollateralized discounting requires the same curve for discounting ({id}), borrowing ({}) and lending ({})",
                        funding.borrow, funding.lend
                    )));
                }
            }
            let curve = structure.resolve_discount(id)?;
            check_currency(curve.spec(), *ccy)?;
            Ok(ResolvedRequirement::Discounting(curve))
        }
        Requirement::Forecast { index } => {
            let id = config
                .forecast_curve(index)
                .ok_or_else(|| Error::MissingData(format!("no forecasting rule for {index}")))?;
            let curve = structure.resolve_projection(id)?;
            // 予想の通貨は指標 ID の文字列から推測せず、参照データの定義から取る。
            // 金利曲線から予想できるのは金利の指標だけで、FX や資産の指標は扱わない。
            let ccy = match &refdata.index(index)?.kind {
                IndexKind::OvernightRate { ccy, .. } | IndexKind::TermRate { ccy, .. } => *ccy,
                _ => {
                    return Err(Error::Unsupported(format!(
                        "rate curve forecasting for {index}"
                    )));
                }
            };
            check_currency(curve.spec(), ccy)?;
            Ok(ResolvedRequirement::Forecast(curve))
        }
    }
}

fn check_currency(curve: &CurveSpec, expected: Currency) -> Result<()> {
    if curve.currency() != expected {
        return Err(Error::InvalidInput(format!(
            "curve {} has currency {}, expected {expected}",
            curve.id(),
            curve.currency()
        )));
    }
    Ok(())
}
`````

`````
crates/rhodium-pricing/src/valuation.rs
//! 取引の評価結果と、報告通貨でのポートフォリオ集計。

use std::{borrow::Borrow, collections::BTreeMap};

use rhodium_core::{
    Error, Result, Warning,
    currency::Currency,
    fx::{CcyPair, FxRate, Today},
    id::TradeId,
    money::{Money, MoneyBag},
    parallel::{self, ExecutionContext, map_partitions},
};
use rhodium_market::view::MarketView;
use scioxide_core::backend::RealValueSum;

use crate::{calibration::CalibratedMarket, config::Method, context::Context, prepare::Prepared};

/// 一件の取引の評価結果。
#[derive(Debug)]
pub struct Valuation {
    /// 取引の ID。
    pub trade: TradeId,
    /// 受渡通貨ごとの PV。ゼロの項も残す。
    pub pv: MoneyBag,
    /// 評価時点の FX で報告通貨に換算した PV。
    pub reporting: Money,
    /// 使用した評価方式。
    pub method: Method,
    /// 受渡通貨ごとの標準誤差。現行の評価方式では `None`。
    pub std_error: Option<MoneyBag>,
    /// 取引の準備で発生した警告。
    pub warnings: Vec<Warning>,
}

/// 取引ごとの結果と、成功した取引の報告通貨での合計。
#[derive(Debug)]
pub struct PortfolioValuation {
    /// 入力順の結果。取引の準備、評価、換算の失敗もその取引の行に残す。
    pub results: Vec<(TradeId, Result<Valuation>)>,
    /// 成功した取引の PV の合計。合計を有限値で表せないときは `Numerical`。
    pub total: Result<Money>,
    /// 合計から除いた取引の ID。除外理由は `results` に含む。
    pub excluded: Vec<TradeId>,
}

/// 準備済みの取引を評価し、`ctx` の価格設定が指定する報告通貨で集計する。
///
/// 各取引は `market` と同じ市場構造に対して準備する。評価で必要な校正グループや
/// 換算レートを利用できない場合は、該当する取引を除外する。取引は実行器のワーカーに
/// 連続した区間で分け、区間内では同じ評価時点の FX レートを使い回す。
/// 成功した取引が 0 件なら、合計は報告通貨のゼロとなる。
///
/// # Errors
///
/// 取消、期限、並列実行器の失敗は呼び出し全体の `Execution` または `Internal`。
/// 取引の失敗と合計の数値計算の失敗は [`PortfolioValuation`] の各結果に残す。
pub fn value_portfolio<P: Borrow<Prepared> + Sync>(
    prepared: &[(TradeId, Result<P>)],
    market: &CalibratedMarket,
    ctx: &Context<'_>,
    execution: ExecutionContext<'_>,
) -> Result<PortfolioValuation> {
    let reporting = ctx.config.reporting_currency();
    let parts = map_partitions(execution, prepared.len(), |range, child| {
        market.market().with_view(|view| {
            let mut rates = BTreeMap::new();
            prepared[range]
                .iter()
                .map(|(trade, result)| {
                    parallel::check(child)?;
                    let value = result.as_ref().map_err(Clone::clone).and_then(|plan| {
                        value_trade(trade, plan.borrow(), &view, reporting, &mut rates)
                    });
                    match value {
                        Err(error @ Error::Execution(_)) => Err(error),
                        result => Ok((trade.clone(), result)),
                    }
                })
                .collect::<Result<Vec<_>>>()
        })
    })?;
    // 取引が 0 件のときの停止は map_partitions が、それ以外は取引の境界の check が確かめる。
    // 区間は連続しているので、区間の順に連結すれば入力順になる。
    let results: Vec<_> = parts
        .into_iter()
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect();
    let excluded = results
        .iter()
        .filter(|(_, result)| result.is_err())
        .map(|(trade, _)| trade.clone())
        .collect();
    // 途中で f64 の範囲を超える組み合わせでも、相殺後の合計を求める。
    let mut sum = RealValueSum::new();
    let total = results
        .iter()
        .filter_map(|(_, result)| result.as_ref().ok())
        .try_for_each(|value| sum.add_value(value.reporting.amount()))
        .and_then(|()| sum.finish())
        .map(|amount| Money::new(reporting, amount))
        .map_err(|error| Error::Numerical(format!("portfolio PV in {reporting}: {error}")));
    Ok(PortfolioValuation {
        results,
        total,
        excluded,
    })
}

/// 一件の取引を評価し、受渡通貨ごとの PV を評価時点の FX で報告通貨に換算する。
///
/// `rates` は区間内で共有する換算レートの控え。取得に失敗したペアもエラーごと控え、
/// 同じ通貨を持つ後続の取引で取得をやり直さない。
fn value_trade(
    trade: &TradeId,
    plan: &Prepared,
    view: &MarketView<'_>,
    reporting: Currency,
    rates: &mut BTreeMap<Currency, Result<FxRate<Today>>>,
) -> Result<Valuation> {
    let pv = plan.bind(view)?.evaluate()?;
    // この取引が持つ報告通貨以外の通貨だけについてレートを集め、換算は MoneyBag::convert に任せる。
    let direct = pv
        .terms()
        .map(Money::ccy)
        .filter(|&ccy| ccy != reporting)
        .map(|ccy| {
            rates
                .entry(ccy)
                .or_insert_with(|| view.fx_today(CcyPair::new(ccy, reporting)?))
                .clone()
        })
        .collect::<Result<Vec<_>>>()?;
    let reporting = pv.convert(reporting, &direct)?;
    Ok(Valuation {
        trade: trade.clone(),
        pv,
        reporting,
        method: plan.method(),
        std_error: None,
        warnings: plan.warnings().to_vec(),
    })
}
`````

`````
crates/rhodium-pricing/tests/calibration_instruments.rs
//! IRS と預金のテンプレートを校正商品に使い、校正と取引の評価が一致することを確かめる。
//!
//! ベーシスで校正した市場で IRS（預金）のパーを求め、それを quote にして曲線を校正し直すと、
//! 同じ θ に戻る（解が一意なので、元の曲線が新しい方程式も満たす）。校正した IRS と預金を
//! quote で約定した取引として評価すると、PV は校正の残差の許容誤差の内でゼロになる。

use rhodium_core::{Error, Result, currency::Currency};
use rhodium_market::{
    market::Market,
    quote::{QuoteId, QuoteKind, QuoteSet},
};
use rhodium_pricing::{
    calibration::{CalibratedMarket, GroupOutcome, calibrate},
    context::Context,
    definition::{CalibrationTarget, ElementSelector, InstrumentBinding, QuoteTransform},
    prepare::{PrepareInput, evaluate, prepare},
    resolve::CollateralKey,
};
use rhodium_products::standard::InstrumentTemplate;
use rhodium_testkit::assert_close;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、校正商品の検査に使わない関数がある"
)]
#[path = "fixtures/ois.rs"]
pub mod ois;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、校正商品の検査に使わない関数がある"
)]
#[path = "fixtures/multicurve.rs"]
pub mod multicurve;

use multicurve::{MulticurveCase, all_multicurve_cases};
use ois::{Node, fixture_quotes};

const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
const TOLERANCE: f64 = 1e-14;

fn solved(market: &CalibratedMarket) -> Result<()> {
    for group in market.groups() {
        let GroupOutcome::Solved(solution) = group.outcome() else {
            return Err(Error::Calibration(format!("{group:?}")));
        };
        assert!(solution.residual_norm() <= TOLERANCE);
    }
    Ok(())
}

/// テンプレートの PV が quote K について affine なので、K = 0 と 1 の PV からパーと傾きを求める。
fn par_and_slope(
    template: &InstrumentTemplate,
    case: &MulticurveCase,
    market: &Market<f64>,
    ctx: &Context<'_>,
) -> Result<(f64, f64)> {
    let pv = |rate| {
        let product = template.build(rate, ctx.refdata)?;
        let plan = prepare(
            PrepareInput::Standard {
                product: &product,
                trade_date: case.overnight.as_of,
            },
            &COLLATERAL,
            &[],
            market.structure(),
            ctx,
        )?;
        Ok::<_, Error>(evaluate(&plan, market)?.get(case.overnight.currency))
    };
    let zero = pv(0.0)?;
    let slope = pv(1.0)? - zero;
    Ok((-zero / slope, slope))
}

fn irs(case: &MulticurveCase, node: &Node) -> Result<InstrumentTemplate> {
    Ok(InstrumentTemplate::Swap {
        convention: case.irs_convention().parse()?,
        tenor: node.tenor,
    })
}

fn assert_same_theta(lhs: &CalibratedMarket, rhs: &CalibratedMarket) {
    // 二つの校正はそれぞれ残差 1e-14 で止まるので、θ はその誤差を伝えた幅で一致する。
    for (lhs, rhs) in lhs.market().theta().iter().zip(rhs.market().theta()) {
        assert_close!(*lhs, *rhs, 1e-10, 1e-12);
    }
}

#[test]
fn irs_quotes_calibrate_the_tibor_curve_and_reprice_at_par() -> Result<()> {
    // JPY の TIBOR 6M と 3M（TONA へのベーシス）と USD のターム指標（独立の節点）で、ベーシスの代わりに
    // 固定対ターム金利の IRS のパーレートでターム予想曲線を当てる。
    for case in all_multicurve_cases()? {
        let inputs = case.inputs()?;
        let ctx = inputs.context();
        let basis = calibrate(&case.definition()?, &case.quotes()?, &ctx)?;
        solved(&basis)?;
        let term_block = ElementSelector::Block(case.term_block()?);
        let mut bindings = Vec::new();
        let mut quotes = fixture_quotes(&case.overnight)?;
        for node in &case.nodes {
            let template = irs(&case, node)?;
            let (par, _) = par_and_slope(&template, &case, basis.market(), &ctx)?;
            let quote: QuoteId = format!("{}.IRS.{}", case.term_index, node.tenor).parse()?;
            bindings.push(InstrumentBinding::new(
                quote.clone(),
                QuoteKind::Rate,
                CalibrationTarget::ParRate(template),
                QuoteTransform::Identity,
                COLLATERAL,
                term_block.clone(),
            )?);
            quotes.push((quote, par));
        }
        let quotes = QuoteSet::new(quotes)?;
        let swaps = calibrate(
            &case.definition_with_term_bindings(bindings)?,
            &quotes,
            &ctx,
        )?;
        solved(&swaps)?;
        assert_same_theta(&basis, &swaps);
        // 校正した IRS を quote の固定金利で約定した取引の PV は、残差の許容誤差 × 傾きの幅でゼロ。
        for node in &case.nodes {
            let template = irs(&case, node)?;
            let id = format!("{}.IRS.{}", case.term_index, node.tenor).parse()?;
            let quote = quotes
                .get(&id)
                .ok_or_else(|| Error::Internal("missing IRS quote".into()))?;
            let (par, slope) = par_and_slope(&template, &case, swaps.market(), &ctx)?;
            assert_close!(par, quote, 0.0, TOLERANCE);
            assert_close!(
                slope * (quote - par),
                0.0,
                0.0,
                2.0 * slope.abs() * TOLERANCE
            );
        }
    }
    Ok(())
}

#[test]
fn deposit_template_calibrates_and_reprices_at_par() -> Result<()> {
    // TONA の最初のピラー（6M）を OIS の代わりにスポット起点の 6M の預金で当てる。
    // 預金の満期は OIS の最終の支払日と同じなので、預金のパーを quote にすると同じ曲線に戻る。
    let [case, _, _] = all_multicurve_cases()?;
    let inputs = case.inputs()?;
    let ctx = inputs.context();
    let original = calibrate(&case.definition()?, &case.quotes()?, &ctx)?;
    let first = &case.overnight.nodes[0];
    let template = InstrumentTemplate::Deposit {
        convention: case.deposit_convention().parse()?,
        tenor: first.tenor.to_string().parse()?,
        currency: Currency::JPY,
    };
    let (par, _) = par_and_slope(&template, &case, original.market(), &ctx)?;
    let deposit: QuoteId = "JPY.DEPO.6M".parse()?;
    let overnight_block = ElementSelector::Block(case.overnight.block()?);
    let term_block = ElementSelector::Block(case.term_block()?);
    let bindings = std::iter::once(InstrumentBinding::new(
        deposit.clone(),
        QuoteKind::Rate,
        CalibrationTarget::ParRate(template.clone()),
        QuoteTransform::Identity,
        COLLATERAL,
        overnight_block.clone(),
    ))
    .chain(case.overnight.nodes[1..].iter().map(|node| {
        case.overnight
            .binding(node, &COLLATERAL, overnight_block.clone())
    }))
    .chain(case.nodes.iter().map(|node| {
        InstrumentBinding::new(
            case.quote(node)?,
            QuoteKind::Spread,
            CalibrationTarget::ParRate(case.template(node)?),
            QuoteTransform::Identity,
            COLLATERAL,
            term_block.clone(),
        )
    }))
    .collect::<Result<Vec<_>>>()?;
    let definition = rhodium_pricing::definition::MarketDefinition::new(
        "JPY".parse()?,
        "Etc/UTC".parse()?,
        [case.term_curve()?, case.overnight.curve()?],
        bindings,
    )?;
    let first_quote = case.overnight.quote(first)?;
    let quotes = QuoteSet::new(
        case.quotes()?
            .iter()
            .filter(|(id, _)| **id != first_quote)
            .map(|(id, value)| (id.clone(), value))
            .chain([(deposit, par)]),
    )?;
    let deposits = calibrate(&definition, &quotes, &ctx)?;
    solved(&deposits)?;
    assert_same_theta(&original, &deposits);
    let (repriced, _) = par_and_slope(&template, &case, deposits.market(), &ctx)?;
    assert_close!(repriced, par, 0.0, TOLERANCE);
    Ok(())
}
`````

`````
crates/rhodium-pricing/tests/config_and_resolve.rs
//! 価格設定から市場ビューへの曲線の解決を確かめる。

use std::sync::Arc;

use rhodium_core::{
    Error, Result,
    calendar::{BaseCalendar, CalendarSpec, Extrapolation as CalendarExtrapolation, WeekendRule},
    currency::{Asset, Currency},
    date::Weekday,
    day_count::DayCount,
    fixing::FixingStore,
    fx::CcyPair,
    id::{CalendarId, IndexId, SecurityId},
    quantity::{Rate, Simple, Spread, YearFraction},
    refdata::{
        AssetValuation, CurrencyDefinition, DailyPublication, Deliverability, FxConvention,
        FxSpotLag, IndexDefinition, IndexKind, IndexStatus, PublicationRule, RefData, RefDataParts,
        SecurityDefinition, SecurityKind, SecurityUnit,
    },
    rounding::{Rounding, RoundingMode},
    time::{Instant, LocalEventTime, UtcOffset, ZoneRules},
};
use rhodium_market::{
    curve::{CurveExpr, CurveRole, CurveSpec, Extrapolation},
    market::Market,
    structure::MarketStructure,
    view::{Discount, Projection, ResolvedCurve, Slots},
};
use rhodium_pricing::{
    config::{FundingTreatment, Method, PricingConfig},
    context::Context,
    prepare::{PrepareInput, evaluate, prepare_with_margin_set},
    resolve::{CollateralKey, Requirement, ResolvedRequirement, resolve},
};
use rhodium_products::{
    agreement::{
        CashCollateralTerms, CollateralAsset, CsaTerms, Eligibility, EligibleCollateral, MarginSet,
        PerParty, Remuneration,
    },
    ir::{Contract, Product, SettlementBoundary, analysis::Linearity},
    standard::ProductKind,
};
use rhodium_testkit::assert_close;

const START: &str = "2026-01-01T00:00:00Z";
const PILLAR: &str = "2027-01-01T00:00:00Z";
const SECURITY: &str = "TEST-SHARE";
const CONFIG: &str = r#"
version = 1
reporting_currency = "JPY"
valuation_close = { time = "17:00", zone = "Asia/Tokyo" }
settlement = "include"

[[discounting]]
currency = "JPY"
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
curve = "JPY-DISCOUNT"

[[discounting]]
currency = "JPY"
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "25bp" } }
curve = "JPY-CSA"

[[discounting]]
currency = "USD"
collateral = "uncollateralized"
curve = "USD-SOFR"

[[methods]]
product = "swap"
method = "linear"

[forecasting]
"JPY-TONA" = "JPY-PROJECTION"
"USD-SOFR" = "USD-SOFR"
"#;

fn parse(text: &str) -> Result<PricingConfig> {
    toml::from_str(text).map_err(|error| Error::InvalidData(error.to_string()))
}

fn fixture_refdata() -> Result<RefData> {
    let mut parts = RefDataParts::default();
    let calendar_id: CalendarId = "WEEKDAY".parse()?;
    let calendar = CalendarSpec::Base(calendar_id.clone());
    parts.calendars.insert(
        calendar_id.clone(),
        BaseCalendar::new(
            "2026-01-01".parse()?..="2027-01-01".parse()?,
            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
            vec![],
            CalendarExtrapolation::Strict,
        )?,
    );
    let close = LocalEventTime {
        time: "17:00:00".parse()?,
        zone: "Asia/Tokyo".parse()?,
    };
    parts.zones.insert(
        close.zone.clone(),
        ZoneRules::new(
            START.parse()?..=PILLAR.parse()?,
            UtcOffset::from_seconds(9 * 3600)?,
            vec![],
        )?,
    );
    for (ccy, minor_units, unit) in [(Currency::JPY, 0, 1.0), (Currency::USD, 2, 0.01)] {
        parts.currencies.insert(
            ccy,
            CurrencyDefinition {
                code: ccy,
                minor_units,
                payment_rounding: Rounding::new(unit, RoundingMode::HalfAwayFromZero)?,
                deliverability: Deliverability::Deliverable,
                settlement_time: close.clone(),
                settlement_calendar: calendar_id.clone(),
            },
        );
    }
    let daily = PublicationRule::Daily(DailyPublication {
        business_days: 1,
        calendar: calendar.clone(),
        time: close.time,
        zone: close.zone,
    });
    for (name, ccy, day_count) in [
        ("JPY-TONA", Currency::JPY, DayCount::Act365F),
        ("USD-SOFR", Currency::USD, DayCount::Act360),
    ] {
        let id: IndexId = name.parse()?;
        parts.indices.insert(
            id.clone(),
            IndexDefinition {
                id,
                kind: IndexKind::OvernightRate { ccy, day_count },
                calendar: calendar.clone(),
                publication: daily.clone(),
                observation_time: None,
                status: IndexStatus::Active,
            },
        );
    }
    // USDJPY の FX fixing の受渡日は、ペアの慣行のスポットラグで決まる。
    let usd_jpy: CcyPair = "USDJPY".parse()?;
    parts.fx_conventions.insert(
        usd_jpy,
        FxConvention {
            pair: usd_jpy,
            spot_lag: FxSpotLag::TwoBusinessDays,
            settlement: parts.currencies[&Currency::JPY].settlement_time.clone(),
            pip_exponent: 2,
        },
    );
    let fx: IndexId = "USDJPY-WMR".parse()?;
    parts.indices.insert(
        fx.clone(),
        IndexDefinition {
            id: fx,
            kind: IndexKind::FxFixing {
                pair: "USDJPY".parse()?,
            },
            calendar,
            publication: PublicationRule::Intraday { delay_micros: 0 },
            observation_time: None,
            status: IndexStatus::Active,
        },
    );
    let security: SecurityId = SECURITY.parse()?;
    parts.securities.insert(
        security.clone(),
        SecurityDefinition {
            id: security,
            kind: SecurityKind::Share,
            currency: Currency::JPY,
            unit: SecurityUnit::Share,
            issuer: None,
            valuation: AssetValuation::Forward,
        },
    );
    RefData::new(parts)
}

fn fixture() -> Result<(Arc<MarketStructure>, RefData)> {
    let [discount, projection] = [CurveRole::Discount, CurveRole::Projection];
    let specs = [
        ("USD-SOFR", Currency::USD, vec![discount, projection]),
        ("JPY-DISCOUNT", Currency::JPY, vec![discount]),
        ("JPY-CSA", Currency::JPY, vec![discount]),
        ("JPY-PROJECTION", Currency::JPY, vec![projection]),
    ]
    .into_iter()
    .map(|(id, ccy, roles)| {
        CurveSpec::new(
            id.parse()?,
            ccy,
            roles,
            [PILLAR.parse()?],
            Extrapolation::default(),
        )
    })
    .collect::<Result<Vec<_>>>()?;
    Ok((
        Arc::new(MarketStructure::new(START.parse()?, specs)?),
        fixture_refdata()?,
    ))
}

fn cash(spread: &str) -> Result<CollateralKey> {
    serde_json::from_value(serde_json::json!({
        "cash": { "ccy": "JPY", "index": "JPY-TONA", "spread": spread }
    }))
    .map_err(|error| Error::InvalidData(error.to_string()))
}

const BOTH: PerParty<Eligibility> = PerParty {
    own: Eligibility::Eligible,
    cpty: Eligibility::Eligible,
};

fn cash_terms(ccy: Currency, index: &str, spread: f64) -> Result<CashCollateralTerms> {
    Ok(CashCollateralTerms {
        ccy,
        remuneration: Remuneration {
            index: index.parse()?,
            spread: Spread::new(spread)?,
            floor: None,
        },
    })
}

fn cash_candidate(terms: CashCollateralTerms) -> Result<EligibleCollateral> {
    EligibleCollateral::new(CollateralAsset::Cash(terms), 0.0, BOTH)
}

fn margin_set(eligible: Vec<EligibleCollateral>, refdata: &RefData) -> Result<MarginSet> {
    Ok(MarginSet {
        id: "MS-TEST".parse()?,
        vm: Some(CsaTerms::new(eligible, refdata)?),
    })
}

fn discounting(ccy: Currency, collateral: CollateralKey) -> Requirement {
    Requirement::Discounting { ccy, collateral }
}

fn forecasting(name: &str) -> Result<Requirement> {
    Ok(Requirement::Forecast {
        index: name.parse()?,
    })
}

fn discount(
    requirement: &Requirement,
    config: &PricingConfig,
    structure: &Arc<MarketStructure>,
    refdata: &RefData,
) -> Result<ResolvedCurve<Discount>> {
    match resolve(requirement, config, structure, refdata)? {
        ResolvedRequirement::Discounting(curve) => Ok(curve),
        ResolvedRequirement::Forecast(_) => Err(Error::Internal("expected discount curve".into())),
    }
}

fn projection(
    index: &str,
    config: &PricingConfig,
    structure: &Arc<MarketStructure>,
    refdata: &RefData,
) -> Result<ResolvedCurve<Projection>> {
    match resolve(&forecasting(index)?, config, structure, refdata)? {
        ResolvedRequirement::Forecast(curve) => Ok(curve),
        ResolvedRequirement::Discounting(_) => {
            Err(Error::Internal("expected projection curve".into()))
        }
    }
}

#[test]
fn pricing_config_reads_defaults_and_explicit_settings() -> Result<()> {
    // 設定の版と報告通貨だけを指定したとき、他の項目には設計上の既定値が入る。
    let json = r#"{"version":1,"reporting_currency":"USD"}"#;
    let defaults: PricingConfig =
        serde_json::from_str(json).map_err(|error| Error::InvalidData(error.to_string()))?;
    assert_eq!(defaults.reporting_currency(), Currency::USD);
    assert_eq!(defaults.settlement(), SettlementBoundary::Exclude);
    assert_eq!(defaults.funding(), FundingTreatment::InXva);
    assert!(defaults.funding_curves().is_empty());
    assert!(defaults.valuation_close().is_none());

    // 明示された現地の締め時刻と受渡境界を同じ型から取得できる。
    let config = parse(CONFIG)?;
    let close = config
        .valuation_close()
        .ok_or_else(|| Error::Internal("missing close".into()))?;
    assert_eq!(close.time, "17:00:00".parse()?);
    assert_eq!(close.zone, "Asia/Tokyo".parse()?);
    let at = START.parse()?;
    assert!(config.settlement().includes(at, at));
    assert!(!defaults.settlement().includes(at, at));
    Ok(())
}

#[test]
fn pricing_config_rejects_invalid_documents() {
    // 誤記や等価な担保条件の重複を読み込み時に報告する。
    for (source, message) in [
        ("version = 1", "reporting_currency"),
        ("reporting_currency = \"JPY\"", "version"),
        ("version = 2\nreporting_currency = \"JPY\"", "version 2"),
        (
            "version = 1\nreporting_currency = \"JPY\"\nsettlemant = \"include\"",
            "unknown field",
        ),
    ] {
        assert!(matches!(parse(source), Err(Error::InvalidData(error)) if error.contains(message)));
    }
    let duplicate = CONFIG.replace("25bp", "0%");
    assert!(
        matches!(parse(&duplicate), Err(Error::InvalidData(error)) if error.contains("duplicate discounting rule"))
    );
    let duplicate = format!("{CONFIG}\n[[methods]]\nproduct = \"swap\"\nmethod = \"linear\"");
    assert!(
        matches!(parse(&duplicate), Err(Error::InvalidData(error)) if error.contains("duplicate method rule"))
    );
    assert!(parse(&CONFIG.replace("\"25bp\"", "25")).is_err());
    assert!(cash("NaNbp").is_err());
    assert!(parse(&CONFIG.replace("method = \"linear\"", "method = \"analytic\"")).is_err());

    // JSON の表でも、同じキーを複数回書いた設定は拒否する。
    let duplicate = r#"{"version":1,"reporting_currency":"JPY","forecasting":{"JPY-TONA":"FIRST","JPY-TONA":"SECOND"}}"#;
    assert!(
        matches!(serde_json::from_str::<PricingConfig>(duplicate), Err(error) if error.to_string().contains("duplicate rule for JPY-TONA"))
    );
}

#[test]
fn method_selection_checks_explicit_and_default_linear_applicability() -> Result<()> {
    // スワップには明示行があり、預金や標準商品に分類されない IR では既定を使う。
    let config = parse(CONFIG)?;
    assert_eq!(
        config.select_method(Some(ProductKind::Swap), Linearity::Linear)?,
        Method::Linear
    );
    assert_eq!(
        config.select_method(Some(ProductKind::Deposit), Linearity::Linear)?,
        Method::Linear
    );
    assert_eq!(
        config.select_method(None, Linearity::Linear)?,
        Method::Linear
    );
    assert!(
        matches!(config.select_method(Some(ProductKind::Swap), Linearity::NonLinear("choice")), Err(Error::Unsupported(error)) if error.contains("configured"))
    );
    assert!(
        matches!(config.select_method(Some(ProductKind::Deposit), Linearity::NonLinear("choice")), Err(Error::Unsupported(error)) if error.contains("default"))
    );
    Ok(())
}

#[test]
fn resolved_jpy_and_usd_curves_bind_with_their_own_points() -> Result<()> {
    // 決められた曲線から割引因子と単利の予想を取り出し、節点の DF から手計算する。
    let config = parse(CONFIG)?;
    let (structure, refdata) = fixture()?;
    let start = START.parse()?;
    let end = PILLAR.parse()?;
    let jpy = discount(
        &discounting(Currency::JPY, cash("0.25%")?),
        &config,
        &structure,
        &refdata,
    )?;
    let usd = discount(
        &discounting(Currency::USD, CollateralKey::Uncollateralized),
        &config,
        &structure,
        &refdata,
    )?;
    let tona = projection("JPY-TONA", &config, &structure, &refdata)?;
    let sofr = projection("USD-SOFR", &config, &structure, &refdata)?;
    assert_eq!(jpy.spec().id().as_str(), "JPY-CSA");
    assert_eq!(usd.spec().id(), sofr.spec().id());
    assert_eq!(tona.spec().id().as_str(), "JPY-PROJECTION");
    let jpy_pay = jpy.prepare_point(end)?;
    let usd_pay = usd.prepare_point(end)?;
    let tona_points = (
        tona.prepare_point(start)?,
        tona.prepare_point(end)?,
        YearFraction::new(1.0, DayCount::Act365F)?,
    );
    let sofr_tau = 365.0 / 360.0;
    let sofr_points = (
        sofr.prepare_point(start)?,
        sofr.prepare_point(end)?,
        YearFraction::new(sofr_tau, DayCount::Act360)?,
    );
    let slots = Slots::from_resolved(
        Arc::clone(&structure),
        [(jpy, vec![jpy_pay]), (usd, vec![usd_pay])],
        [(tona, vec![tona_points]), (sofr, vec![sofr_points])],
    )?;
    // 構造内の順序は JPY-CSA、JPY-DISCOUNT、JPY-PROJECTION、USD-SOFR。
    let market = Market::from_theta(
        structure,
        vec![0.95_f64.ln(), 0.96_f64.ln(), 0.98_f64.ln(), 0.94_f64.ln()],
    )?;
    market.with_view(|view| {
        let bound = view.bind_slots(&slots)?;
        assert_close!(
            view.df(bound.discounts()[0].points()[0]).value(),
            0.95,
            1e-12,
            1e-14
        );
        assert_close!(
            view.df(bound.discounts()[1].points()[0]).value(),
            0.94,
            1e-12,
            1e-14
        );
        assert_close!(
            view.forward(bound.projections()[0].points()[0]),
            1.0 / 0.98 - 1.0,
            1e-12,
            1e-14
        );
        assert_close!(
            view.forward(bound.projections()[1].points()[0]),
            (1.0 / 0.94 - 1.0) / sofr_tau,
            1e-12,
            1e-14
        );
        Ok(())
    })
}

#[test]
fn missing_rules_do_not_fall_back_to_other_curves() -> Result<()> {
    // 同じ通貨の他の担保行や、同名の市場曲線が存在しても、必要な表の行を要求する。
    let (structure, refdata) = fixture()?;
    let config = parse(CONFIG)?;
    for requirement in [
        discounting(Currency::JPY, cash("5bp")?),
        discounting(Currency::JPY, CollateralKey::Uncollateralized),
    ] {
        assert!(matches!(
            resolve(&requirement, &config, &structure, &refdata),
            Err(Error::MissingData(_))
        ));
    }
    let missing_curve =
        parse(&CONFIG.replace("\"USD-SOFR\" = \"USD-SOFR\"", "\"USD-SOFR\" = \"NO-CURVE\""))?;
    let sofr = forecasting("USD-SOFR")?;
    assert!(matches!(
        resolve(&sofr, &missing_curve, &structure, &refdata),
        Err(Error::MissingData(_))
    ));
    let empty = parse("version = 1\nreporting_currency = \"JPY\"")?;
    assert!(matches!(
        resolve(&sofr, &empty, &structure, &refdata),
        Err(Error::MissingData(_))
    ));
    Ok(())
}

#[test]
fn resolve_rejects_wrong_curve_roles_and_currencies() -> Result<()> {
    // 予想と割引それぞれで、市場の役割と通貨が要求に合うことを検証する。
    let (structure, refdata) = fixture()?;
    let sofr = forecasting("USD-SOFR")?;
    for curve in ["JPY-DISCOUNT", "JPY-PROJECTION"] {
        let wrong = parse(&CONFIG.replace(
            "\"USD-SOFR\" = \"USD-SOFR\"",
            &format!("\"USD-SOFR\" = \"{curve}\""),
        ))?;
        assert!(matches!(
            resolve(&sofr, &wrong, &structure, &refdata),
            Err(Error::InvalidInput(_))
        ));
    }
    let config = parse(&format!(
        "{CONFIG}\n\"USDJPY-WMR\" = \"USD-SOFR\"\n\"UNKNOWN-INDEX\" = \"USD-SOFR\""
    ))?;
    let cash_discount = discounting(Currency::JPY, cash("0bp")?);
    for curve in ["JPY-PROJECTION", "USD-SOFR"] {
        let wrong =
            parse(&CONFIG.replace("curve = \"JPY-DISCOUNT\"", &format!("curve = \"{curve}\"")))?;
        assert!(matches!(
            resolve(&cash_discount, &wrong, &structure, &refdata),
            Err(Error::InvalidInput(_))
        ));
    }
    assert!(matches!(
        resolve(&forecasting("USDJPY-WMR")?, &config, &structure, &refdata),
        Err(Error::Unsupported(_))
    ));
    assert!(matches!(
        resolve(
            &forecasting("UNKNOWN-INDEX")?,
            &config,
            &structure,
            &refdata
        ),
        Err(Error::InvalidData(_))
    ));
    Ok(())
}

#[test]
fn funding_curves_single_source() -> Result<()> {
    // 無担保の割引と資金調達が同じ曲線を読むことを、評価の準備時に照合する。
    let (structure, refdata) = fixture()?;
    let source = format!("funding = \"in_discounting\"\n{CONFIG}");
    let requirement = discounting(Currency::USD, CollateralKey::Uncollateralized);
    let with_funding = |borrow, lend| {
        parse(&format!(
            "{source}\n[funding_curves.USD]\nborrow = \"{borrow}\"\nlend = \"{lend}\""
        ))
    };
    let config = with_funding("USD-SOFR", "USD-SOFR")?;
    let funding_curve = &config.funding_curves()[&Currency::USD].borrow;
    assert_eq!(
        discount(&requirement, &config, &structure, &refdata)?
            .spec()
            .id(),
        funding_curve
    );
    for (borrow, lend) in [("USD-FUNDING", "USD-FUNDING"), ("USD-SOFR", "USD-LEND")] {
        let wrong = with_funding(borrow, lend)?;
        assert!(matches!(
            resolve(&requirement, &wrong, &structure, &refdata),
            Err(Error::InvalidInput(_))
        ));
    }

    // 無担保の行が無い場合と、参照する通貨に資金調達の表が無い場合を区別する。
    let missing_rule = parse(
        "version = 1\nreporting_currency = \"JPY\"\nfunding = \"in_discounting\"\n[funding_curves.USD]\nborrow = \"USD-SOFR\"\nlend = \"USD-SOFR\"",
    )?;
    assert!(
        matches!(resolve(&requirement, &missing_rule, &structure, &refdata), Err(Error::MissingData(message)) if message.contains("discounting"))
    );
    let missing_funding = parse(&source)?;
    assert!(
        matches!(resolve(&requirement, &missing_funding, &structure, &refdata), Err(Error::MissingData(message)) if message.contains("funding"))
    );
    let secured = discounting(Currency::JPY, cash("0bp")?);
    assert!(resolve(&secured, &missing_funding, &structure, &refdata).is_ok());
    Ok(())
}

fn csa_market(as_of: Instant, first_year: Instant, second_year: Instant) -> Result<Market> {
    let curves = [
        ("JPY-TONA", Currency::JPY, CurveExpr::Flat(0.01)),
        (
            "USD-TONA",
            Currency::USD,
            CurveExpr::Stitch {
                short: Box::new(CurveExpr::Flat(0.05)),
                long: Box::new(CurveExpr::Flat(0.01)),
                at: first_year,
            },
        ),
        (
            "USD-SOFR",
            Currency::USD,
            CurveExpr::Stitch {
                short: Box::new(CurveExpr::Flat(0.02)),
                long: Box::new(CurveExpr::Flat(0.04)),
                at: first_year,
            },
        ),
        ("USD-UNSECURED", Currency::USD, CurveExpr::Flat(0.035)),
        (
            "USD-CTD",
            Currency::USD,
            CurveExpr::MaxForward {
                parts: vec![
                    CurveExpr::Ref("USD-TONA".parse()?),
                    CurveExpr::Ref("USD-SOFR".parse()?),
                ],
                grid: vec![first_year, second_year],
            },
        ),
    ]
    .into_iter()
    .map(|(id, ccy, expr)| CurveSpec::from_expr(id.parse()?, ccy, [CurveRole::Discount], expr))
    .collect::<Result<Vec<_>>>()?;
    let structure = Arc::new(MarketStructure::new(as_of, curves)?);
    Market::from_theta(structure, vec![])
}

const CSA_CONFIG: &str = r#"
version = 1
reporting_currency = "JPY"

[[discounting]]
currency = "JPY"
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
curve = "JPY-TONA"

[[discounting]]
currency = "USD"
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
curve = "USD-TONA"

[[discounting]]
currency = "USD"
collateral = { cheapest_to_deliver = [{ ccy = "JPY", index = "JPY-TONA", spread = "0bp" }, { ccy = "USD", index = "USD-SOFR", spread = "0bp" }] }
curve = "USD-CTD"

[[discounting]]
currency = "USD"
collateral = "uncollateralized"
curve = "USD-UNSECURED"
"#;

fn cash_payment(currencies: &[Currency], at: Instant, refdata: &RefData) -> Result<Product> {
    Product::new(
        Contract::build(|builder| {
            let amount = builder.constant(100.0);
            let legs = currencies
                .iter()
                .map(|&ccy| builder.pay(Asset::Cash(ccy), amount, at))
                .collect::<Vec<_>>();
            Ok(builder.all(legs))
        })?,
        refdata,
    )
}

#[test]
fn csa_discounting_matches_manual_pv_for_four_cases() -> Result<()> {
    // CSA の 4 つの場合（自国通貨、外貨、CTD、無担保）で、マージンセットから導いた鍵が
    // それぞれの割引行を引き、PV が曲線から手計算した値と一致することを確かめる。
    // まず同じ TONA 担保で JPY と USD の支払いを評価し、支払い通貨で割引行が変わることを見る。
    let refdata = fixture_refdata()?;
    let as_of = START.parse()?;
    let second_year = "2028-01-01T00:00:00Z".parse()?;
    let market = csa_market(as_of, PILLAR.parse()?, second_year)?;
    let config = parse(CSA_CONFIG)?;
    let ctx = Context {
        as_of,
        refdata: &refdata,
        fixings: &FixingStore::default(),
        config: &config,
    };
    let terms = cash_terms(Currency::JPY, "JPY-TONA", 0.0)?;
    let tona = cash_candidate(terms)?;
    let sofr = cash_candidate(cash_terms(Currency::USD, "USD-SOFR", 0.0)?)?;
    let single = margin_set(vec![tona.clone()], &refdata)?;
    let ctd = margin_set(vec![sofr, tona.clone(), tona], &refdata)?;
    let unsecured = MarginSet {
        id: "MS-UNSECURED".parse()?,
        vm: None,
    };
    let pv = |product, margin_set| {
        let prepared = prepare_with_margin_set(
            PrepareInput::Ir(product),
            margin_set,
            &[],
            market.structure(),
            &ctx,
        )?;
        evaluate(&prepared, &market)
    };
    let payment_both = cash_payment(&[Currency::JPY, Currency::USD], second_year, &refdata)?;
    let both_pv = pv(&payment_both, &single)?;
    assert_close!(
        both_pv.get(Currency::JPY),
        100.0 * (-0.02_f64).exp(),
        1e-12,
        1e-6
    );
    let tona_forward = [0.05_f64, 0.01];
    let sofr_forward = [0.02_f64, 0.04];
    assert_close!(
        both_pv.get(Currency::USD),
        100.0 * (-tona_forward.iter().sum::<f64>()).exp(),
        1e-12,
        1e-8
    );

    // CSA の候補に重複があっても USD-CTD の行を引く。最初の年は TONA 担保、
    // 次の年は SOFR 担保のフォワードで積算する。無担保は別の曲線で割り引く。
    let ctd_forward = tona_forward[0].max(sofr_forward[0]) + tona_forward[1].max(sofr_forward[1]);
    let payment_usd = cash_payment(&[Currency::USD], second_year, &refdata)?;
    assert_close!(
        pv(&payment_usd, &ctd)?.get(Currency::USD),
        100.0 * (-ctd_forward).exp(),
        1e-12,
        1e-8
    );
    assert_close!(
        pv(&payment_usd, &unsecured)?.get(Currency::USD),
        100.0 * (-0.07_f64).exp(),
        1e-12,
        1e-8
    );

    // USD の無担保行と CTD 行があっても、単独の SOFR 担保には対応する割引行が要る。
    let sofr_only = margin_set(
        vec![cash_candidate(cash_terms(Currency::USD, "USD-SOFR", 0.0)?)?],
        &refdata,
    )?;
    assert!(matches!(
        pv(&payment_usd, &sofr_only),
        Err(Error::MissingData(message)) if message.contains("discounting rule for USD")
    ));
    Ok(())
}

#[test]
fn csa_remuneration_floor_matches_manual_domestic_and_foreign_pv() -> Result<()> {
    // TONA は −2%、3%。50 bp の付利スプレッドと 0% のフロアでは、初年は
    // 2%、翌年は 0.5% を基準曲線に足す。外貨の支払いにも同じ補正を使う。
    let as_of = START.parse()?;
    let first = PILLAR.parse()?;
    let second = "2028-01-01T00:00:00Z".parse()?;
    let steps = |short, long| CurveExpr::Stitch {
        short: Box::new(CurveExpr::Flat(short)),
        long: Box::new(CurveExpr::Flat(long)),
        at: first,
    };
    let floor = |base: &str| -> Result<_> {
        Ok(CurveExpr::Sum(vec![
            CurveExpr::Ref(base.parse()?),
            CurveExpr::MaxForward {
                parts: vec![
                    CurveExpr::Flat(0.005),
                    CurveExpr::Sum(vec![
                        CurveExpr::Flat(0.0),
                        CurveExpr::Scale {
                            factor: -1.0,
                            body: Box::new(CurveExpr::ForwardRef("JPY-TONA".parse()?)),
                        },
                    ]),
                ],
                grid: vec![first, second],
            },
        ]))
    };
    let curves = [
        CurveSpec::from_expr(
            "JPY-TONA".parse()?,
            Currency::JPY,
            [CurveRole::Discount, CurveRole::Projection],
            steps(-0.02, 0.03),
        )?,
        CurveSpec::from_expr(
            "USD-TONA".parse()?,
            Currency::USD,
            [CurveRole::Discount],
            steps(0.05, 0.01),
        )?,
        CurveSpec::from_expr(
            "JPY-FLOOR".parse()?,
            Currency::JPY,
            [CurveRole::Discount],
            floor("JPY-TONA")?,
        )?,
        CurveSpec::from_expr(
            "USD-FLOOR".parse()?,
            Currency::USD,
            [CurveRole::Discount],
            floor("USD-TONA")?,
        )?,
    ];
    let structure = Arc::new(MarketStructure::new(as_of, curves)?);
    let market = Market::from_theta(structure, vec![])?;
    let config = parse(
        r#"
version = 1
reporting_currency = "JPY"

[[discounting]]
currency = "JPY"
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "50bp", floor = { rate = "0bp", day_count = "ACT365F" } } }
curve = "JPY-FLOOR"

[[discounting]]
currency = "USD"
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "50bp", floor = { rate = "0bp", day_count = "ACT365F" } } }
curve = "USD-FLOOR"
"#,
    )?;
    let refdata = fixture_refdata()?;
    let ctx = Context {
        as_of,
        refdata: &refdata,
        fixings: &FixingStore::default(),
        config: &config,
    };
    let mut terms = cash_terms(Currency::JPY, "JPY-TONA", 0.005)?;
    terms.remuneration.floor = Some(Rate::new(0.0, DayCount::Act365F, Simple)?);
    let margin = margin_set(vec![cash_candidate(terms)?], &refdata)?;
    for (at, jpy_log, usd_log) in [(first, 0.0_f64, -0.07_f64), (second, -0.035, -0.085)] {
        let product = cash_payment(&[Currency::JPY, Currency::USD], at, &refdata)?;
        let prepared = prepare_with_margin_set(
            PrepareInput::Ir(&product),
            &margin,
            &[],
            market.structure(),
            &ctx,
        )?;
        let pv = evaluate(&prepared, &market)?;
        assert_close!(pv.get(Currency::JPY), 100.0 * jpy_log.exp(), 1e-12, 1e-6);
        assert_close!(pv.get(Currency::USD), 100.0 * usd_log.exp(), 1e-12, 1e-8);
    }
    Ok(())
}

#[test]
fn csa_cash_key_matches_complete_remuneration() -> Result<()> {
    // 同じ通貨・指標でもスプレッドとフロアが違う条件は別の割引行を使う。
    let (structure, refdata) = fixture()?;
    let config = parse(&CONFIG.replace(
        "[[methods]]",
        r#"[[discounting]]
currency = "JPY"
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0%", floor = { rate = "0bp", day_count = "ACT365F" } } }
curve = "JPY-CSA"

[[methods]]"#,
    ))?;
    let plain = cash_terms(Currency::JPY, "JPY-TONA", 0.0)?;
    let with_spread = cash_terms(Currency::JPY, "JPY-TONA", 0.0025)?;
    let mut with_floor = plain.clone();
    with_floor.remuneration.floor = Some(Rate::new(0.0, DayCount::Act365F, Simple)?);
    let request = |terms| {
        let margin = margin_set(vec![cash_candidate(terms)?], &refdata)?;
        Ok(discounting(
            Currency::JPY,
            CollateralKey::from_margin_set(&margin)?,
        ))
    };
    for (terms, curve) in [
        (plain, "JPY-DISCOUNT"),
        (with_spread, "JPY-CSA"),
        (with_floor.clone(), "JPY-CSA"),
    ] {
        assert_eq!(
            discount(&request(terms)?, &config, &structure, &refdata)?
                .spec()
                .id()
                .as_str(),
            curve
        );
    }
    with_floor.remuneration.floor = Some(Rate::new(0.0, DayCount::Act360, Simple)?);
    assert!(matches!(
        resolve(&request(with_floor)?, &config, &structure, &refdata),
        Err(Error::MissingData(_))
    ));
    Ok(())
}

#[test]
fn csa_excludes_asymmetric_and_security_collateral() -> Result<()> {
    // 証券と片側だけが差し入れられる現金を除外し、ヘアカットの違いで同じ現金条件を二重に数えない。
    let refdata = fixture_refdata()?;
    let bilateral = cash_terms(Currency::JPY, "JPY-TONA", 0.0)?;
    let unilateral = cash_terms(Currency::USD, "USD-SOFR", 0.0)?;
    let mut ignored = vec![EligibleCollateral::new(
        CollateralAsset::Security(SECURITY.parse()?),
        0.1,
        BOTH,
    )?];
    for posters in [
        PerParty {
            cpty: Eligibility::Ineligible,
            ..BOTH
        },
        PerParty {
            own: Eligibility::Ineligible,
            ..BOTH
        },
    ] {
        ignored.push(EligibleCollateral::new(
            CollateralAsset::Cash(unilateral.clone()),
            0.0,
            posters,
        )?);
    }
    let no_cash = margin_set(ignored.clone(), &refdata)?;
    assert!(matches!(
        CollateralKey::from_margin_set(&no_cash),
        Err(Error::Unsupported(message)) if message.contains("MS-TEST")
    ));

    ignored.push(cash_candidate(bilateral.clone())?);
    ignored.push(EligibleCollateral::new(
        CollateralAsset::Cash(bilateral.clone()),
        0.1,
        BOTH,
    )?);
    let margin = margin_set(ignored, &refdata)?;
    assert_eq!(
        CollateralKey::from_margin_set(&margin)?,
        CollateralKey::Cash(bilateral)
    );
    Ok(())
}

#[test]
fn pricing_config_normalizes_ctd_candidates() {
    // 単位・候補順・同じ条件の繰り返しは CTD の集合を変えず、同じ割引行は拒否する。
    let jpy = r#"{ ccy = "JPY", index = "JPY-TONA", spread = "0bp" }"#;
    let jpy_zero = r#"{ ccy = "JPY", index = "JPY-TONA", spread = "-0%" }"#;
    let usd = r#"{ ccy = "USD", index = "USD-SOFR", spread = "25bp", floor = { rate = "10bp", day_count = "ACT360" } }"#;
    let usd_percent = r#"{ ccy = "USD", index = "USD-SOFR", spread = "0.25%", floor = { rate = "0.10%", day_count = "ACT360" } }"#;
    let rule = |candidates: &str| {
        format!(
            r#"
[[discounting]]
currency = "USD"
collateral = {{ cheapest_to_deliver = [{candidates}] }}
curve = "USD-SOFR"
"#
        )
    };
    let base = "version = 1\nreporting_currency = \"JPY\"\n";
    let first = rule(&format!("{jpy}, {usd}"));
    let other = rule(&format!("{usd_percent}, {jpy_zero}, {usd}"));
    assert!(parse(&format!("{base}{first}")).is_ok());
    assert!(parse(&format!("{base}{other}")).is_ok());
    assert!(matches!(
        parse(&format!("{base}{first}{other}")),
        Err(Error::InvalidData(message)) if message.contains("duplicate discounting rule")
    ));
    for candidates in [String::new(), jpy.to_owned(), format!("{jpy}, {jpy_zero}")] {
        assert!(matches!(
            parse(&format!("{base}{}", rule(&candidates))),
            Err(Error::InvalidData(message)) if message.contains("at least two distinct")
        ));
    }
}

#[test]
fn csa_construction_and_ctd_check_references() -> Result<()> {
    // 契約の入力は構築時に止める。CTD の参照検証では並びの後ろにある候補も照合する。
    let refdata = fixture_refdata()?;
    assert!(matches!(
        CsaTerms::new(vec![], &refdata),
        Err(Error::InvalidData(_))
    ));
    let asset = CollateralAsset::Cash(cash_terms(Currency::JPY, "JPY-TONA", 0.0)?);
    for haircut in [-0.01, 1.0, f64::NAN] {
        assert!(matches!(
            EligibleCollateral::new(asset.clone(), haircut, BOTH),
            Err(Error::InvalidData(_))
        ));
    }
    let unknown = cash_terms(Currency::USD, "USD-UNKNOWN", 0.0)?;
    assert!(matches!(
        margin_set(vec![cash_candidate(unknown)?], &refdata),
        Err(Error::InvalidData(message)) if message.contains("USD-UNKNOWN")
    ));
    let ctd: CollateralKey = serde_json::from_str(
        r#"{"cheapest_to_deliver":[{"ccy":"JPY","index":"JPY-TONA","spread":"0bp"},{"ccy":"USD","index":"USD-UNKNOWN","spread":"0bp"}]}"#,
    )
    .map_err(|error| Error::InvalidData(error.to_string()))?;
    assert!(matches!(
        ctd.validate_references(&refdata),
        Err(Error::InvalidData(message)) if message.contains("USD-UNKNOWN")
    ));
    Ok(())
}
`````

`````
crates/rhodium-pricing/tests/fixtures/futures.rs
//! 一つの USD 曲線、金利先物の清算指標と限月、LGM のマークからなる M1 の先物と FX の試験用の市場。
//!
//! 週末だけの暦（2025〜2026 年）、UTC の 16 時の公表と受渡、1W の USD ターム指標と、その指標を
//! 原資産とする 2026 年 1 月限の先物（乗数 1,000、最終清算 2026-01-07）を持つ。

use std::{collections::BTreeMap, sync::Arc};

use rhodium_core::{
    Error, Result,
    calendar::{
        BaseCalendar, BusinessDayConvention, CalendarSpec, EndOfMonthRule,
        Extrapolation as CalendarExtrapolation, WeekendRule,
    },
    currency::Currency,
    date::Weekday,
    day_count::DayCount,
    id::{CalendarId, IndexId, SecurityId},
    refdata::{
        AssetValuation, Convention, CurrencyDefinition, DailyPublication, Deliverability,
        IndexDefinition, IndexKind, IndexStatus, PublicationRule, RateFutureContract,
        RateFutureConvention, RateFutureRate, RefData, RefDataParts, SecurityDefinition,
        SecurityKind, SecurityUnit,
    },
    rounding::{Rounding, RoundingMode},
    tenor::DateOffset,
    time::{Instant, LocalEventTime, TzId, UtcOffset, ZoneRules, model_time},
};
use rhodium_market::{
    curve::{CurveRole, CurveSpec, Extrapolation},
    market::Market,
    quote::QuoteKind,
    structure::{BlockId, MarketStructure},
};
use rhodium_pricing::{
    config::PricingConfig,
    definition::{
        CalibrationTarget, ElementSelector, InstrumentBinding, LgmDefinition, MarketDefinition,
        QuoteTransform,
    },
    resolve::CollateralKey,
};
use rhodium_products::standard::InstrumentTemplate;

pub(crate) const START: &str = "2026-01-05";
pub(crate) const AS_OF: &str = "2026-01-05T00:00:00Z";
