
use crate::{
    adjoint,
    gradient::{self, ParamGradient},
    vector::RiskVector,
};

/// quote デルタと、集計から除いた取引の評価エラー。
#[derive(Debug)]
pub struct DeltaReport {
    by_ccy: Vec<RiskVector<QuoteId>>,
    failures: Vec<(usize, Error)>,
}

impl DeltaReport {
    /// 通貨順に並んだ quote デルタ。報告通貨を指定した集計では、その通貨のベクトルを返す。
    ///
    /// 値は quote の種類が定める変化幅あたりの金額。金利とベーシスは絶対 1 bp、
    /// FX スポットは相対 1%、FX フォワードポイントは 1 pip。
    #[must_use]
    pub fn by_ccy(&self) -> &[RiskVector<QuoteId>] {
        &self.by_ccy
    }

    /// 評価に失敗した取引の、入力スライスでの 0 始まりの位置とエラー。
    #[must_use]
    pub fn failures(&self) -> &[(usize, Error)] {
        &self.failures
    }
}

/// 成功した取引の θ 勾配を通貨別に合計し、quote デルタを求める。
///
/// 各取引は `cm` と同じ市場構造に対して準備する。取引ごとの評価エラーは
/// [`DeltaReport::failures`] に残す。取引が全件決済済みなら結果は空になる。
///
/// 取引は `execution` のワーカーに連続した区間で分ける。各ワーカーは tape を一つ作って
/// 区間の取引に使い回し、通貨別の θ 勾配を足す。区間ごとの和を区間の順に足してから、
/// 校正の随伴を通貨ごとに一度だけ解く（右辺の数は取引の数によらない）。
/// 和の順はワーカーの数で変わるため、結果はスレッド数によらず許容誤差の範囲で一致する。
///
/// 取消と期限は取引の境界と AD の記録・逆伝播の中で確かめる。止まった実行は取引の失敗に
/// 数えず、呼び出し全体の `Execution` にする（一部の取引を除いた結果を成功として返さない）。
///
/// # Errors
///
/// θ の集計や校正の随伴に失敗した場合は呼び出し全体のエラーを返す。
/// 並列の実行が取消か期限で止まったか、実行器が失敗した場合は `Execution` か `Internal`。
pub fn quote_delta(
    prepared: &[Prepared],
    cm: &CalibratedMarket,
    execution: ExecutionContext<'_>,
) -> Result<DeltaReport> {
    calculate(prepared, cm, execution, None)
}

/// 成功した取引の quote デルタを、指定した報告通貨で合計する。
///
/// 評価時点の FX で換算した現地通貨の勾配と、評価時点の FX 自身の勾配を合わせてから
/// 校正の随伴を一度解く。全件が決済済みか評価に失敗した場合、ベクトルの列は空になる。
///
/// # Errors
///
/// [`quote_delta`] のエラーに加え、換算に必要な FX 経路や市場データの欠損は
/// `MissingData`、換算やその勾配が非有限なら `Numerical`。
pub fn quote_delta_in(
    prepared: &[Prepared],
    cm: &CalibratedMarket,
    reporting: Currency,
    execution: ExecutionContext<'_>,
) -> Result<DeltaReport> {
    calculate(prepared, cm, execution, Some(reporting))
}

fn calculate(
    prepared: &[Prepared],
    cm: &CalibratedMarket,
    execution: ExecutionContext<'_>,
    reporting: Option<Currency>,
) -> Result<DeltaReport> {
    let layout = cm.market().structure().layout();
    let partials = parallel::map_partitions(execution, prepared.len(), |range, child| {
        let mut tape = gradient::new_tape();
        let mut sums = BTreeMap::new();
        let mut pv = MoneyBag::new();
        let mut failures = Vec::new();
        for position in range {
            parallel::check(child)?;
            match gradient::with_tape(&prepared[position], cm, &mut tape, child) {
                Ok(gradient) => {
                    accumulate(&mut sums, &gradient, layout)?;
                    if let Some(reporting) = reporting {
                        for term in gradient.pv().terms().filter(|term| term.ccy() != reporting) {
                            pv.add(term);
                        }
                    }
                }
                // 実行の制御の失敗は取引の性質ではないので、区間ごと止める。
                Err(error @ Error::Execution(_)) => return Err(error),
                Err(error) => failures.push((position, error)),
            }
        }
        Ok::<_, Error>((sums, pv, failures))
    })?;

    // 区間の順に足す。失敗の位置は区間が連続しているので、連結しても昇順のまま。
    let mut by_ccy: BTreeMap<Currency, Vec<f64>> = BTreeMap::new();
    let mut pv = MoneyBag::new();
    let mut failures = Vec::new();
    for partial in partials {
        let (sums, partial_pv, partial_failures) = partial?;
        for term in partial_pv.terms() {
            pv.add(term);
        }
        for (ccy, theta) in sums {
            match by_ccy.entry(ccy) {
                Entry::Vacant(entry) => {
                    entry.insert(theta);
                }
                Entry::Occupied(mut entry) => {
                    for (total, value) in entry.get_mut().iter_mut().zip(theta) {
                        *total += value;
                    }
                }
            }
        }
        failures.extend(partial_failures);
    }
    if let Some(reporting) = reporting {
        parallel::check(execution)?;
        by_ccy = translate(by_ccy, &pv, cm, reporting, execution)?;
    }
    for (ccy, theta) in &by_ccy {
        if let Some(index) = theta.iter().position(|value| !value.is_finite()) {
            return Err(Error::Numerical(format!(
                "non-finite {ccy} gradient sum at theta index {index}"
            )));
        }
    }

    let by_ccy = adjoint::solve(cm, by_ccy)?
        .into_iter()
        .map(|(ccy, quotes)| {
            let entries = quotes
                .into_iter()
                .map(|(id, value)| {
                    let scaled = quote_unit(cm, &id)? * value;
                    if !scaled.is_finite() {
                        return Err(Error::Numerical(format!("non-finite {ccy} delta for {id}")));
                    }
                    Ok((id, scaled))
                })
                .collect::<Result<_>>()?;
            Ok(RiskVector::from_sorted(ccy, entries))
        })
        .collect::<Result<_>>()?;
    Ok(DeltaReport { by_ccy, failures })
}

/// 受渡通貨別の θ 勾配を、評価時点の FX で報告通貨 r の一本の勾配にまとめる。
///
/// 報告通貨の PV と、その θ 勾配は次のとおり。
///
/// ```text
/// V_r = Σ_c X_{c→r}(0)·V_c
/// g_r = Σ_c X_{c→r}(0)·g_c + Σ_c V_c·∇_θ X_{c→r}(0)
/// ```
///
/// 第 1 項は現地の勾配の換算（r 自身は X = 1）、第 2 項は集計した PV を種にした
/// FX 自身の勾配で、FX のマークと X(0) を決める曲線の θ に入る。どちらも θ について線形なので、
/// 足してから校正の随伴を一度だけ解けば、FX と上流の quote にも戻る。
/// 全通貨が r なら換算せずに返す。
fn translate(
    mut by_ccy: BTreeMap<Currency, Vec<f64>>,
    pv: &MoneyBag,
    cm: &CalibratedMarket,
    reporting: Currency,
    execution: ExecutionContext<'_>,
) -> Result<BTreeMap<Currency, Vec<f64>>> {
    if by_ccy.keys().all(|&ccy| ccy == reporting) {
        return Ok(by_ccy);
    }
    let reporting_gradient = by_ccy.remove(&reporting);
    let structure = cm.market().structure();
    let (plans, pvs): (Vec<_>, Vec<_>) = by_ccy
        .keys()
        .map(|&ccy| {
            let amount = pv.get(ccy);
            if !amount.is_finite() {
                return Err(Error::Numerical(format!("non-finite {ccy} PV sum")));
            }
            Ok((
                structure.fx_today_plan(CcyPair::new(ccy, reporting)?)?,
                amount,
            ))
        })
        .collect::<Result<_>>()?;
    let (rates, mut theta) = gradient::fx_translation(cm, &plans, &pvs, execution)?;
    // 換算レート自身の勾配に、換算先と各受渡通貨の現地勾配を足す。
    let gradients = rates.into_iter().zip(by_ccy.into_values());
    for (rate, gradient) in reporting_gradient
        .into_iter()
        .map(|g| (1.0, g))
        .chain(gradients)
    {
        for (total, value) in theta.iter_mut().zip(gradient) {
            *total += rate * value;
        }
    }
    Ok(BTreeMap::from([(reporting, theta)]))
}

/// quote の種類が定める報告の単位 u。相対の単位（FX スポット、株価）では u·q を返す。
///
/// 随伴の結果は quote の値 1 あたりなので、最後に一度だけこの値を掛ける。
/// 表は設計書「quote」の感応度の単位（`docs/design/06-market/01-quotes.md`）に従う。
fn quote_unit(cm: &CalibratedMarket, id: &QuoteId) -> Result<f64> {
    let kind = cm
        .quote_kind(id)
        .ok_or_else(|| Error::Internal(format!("unknown quote kind for {id}")))?;
    Ok(match kind {
        QuoteKind::Rate
        | QuoteKind::Spread
        | QuoteKind::CreditSpread
        | QuoteKind::FuturesPrice
        | QuoteKind::Upfront
        | QuoteKind::NormalVol
        | QuoteKind::NormalVolDifference => 1e-4,
        QuoteKind::Price
        | QuoteKind::LognormalVol
        | QuoteKind::ShiftedLognormalVol
        | QuoteKind::LognormalVolDifference
        | QuoteKind::Correlation
        | QuoteKind::Recovery => 0.01,
        QuoteKind::FxSpot | QuoteKind::AssetSpot => {
            0.01 * cm
                .quotes()
                .get(id)
                .ok_or_else(|| Error::Internal(format!("missing calibrated quote {id}")))?
        }
        // FX ポイントの入力は pip 数。ペアの pip 幅は校正残差の quote 微分に入っている。
        QuoteKind::FxForwardPoints | QuoteKind::Mark => 1.0,
    })
}

// 一つの取引の通貨別の疎な θ 勾配を、区間の通貨別の平坦な和に足す。
fn accumulate(
    sums: &mut BTreeMap<Currency, Vec<f64>>,
    gradient: &ParamGradient,
    layout: &ParamLayout,
) -> Result<()> {
    for vector in gradient.by_ccy() {
        let total = sums
            .entry(vector.ccy())
            .or_insert_with(|| vec![0.0; layout.len()]);
        for (id, value) in vector.entries() {
            let index = layout.index(id).ok_or_else(|| {
                Error::Internal(format!("unknown parameter in trade gradient: {id:?}"))
            })?;
            total[index] += value;
            if !total[index].is_finite() {
                return Err(Error::Numerical(format!(
                    "non-finite {} gradient sum for {id:?}",
                    vector.ccy()
                )));
            }
        }
    }
    Ok(())
}
`````

`````
crates/rhodium-risk/src/gradient.rs
//! 準備済み商品の PV と、市場パラメータに対する AAD 勾配。

use std::{
    cell::{Cell, RefCell},
    collections::BTreeSet,
    ops::Range,
};

use rhodium_core::{
    Error, Result,
    money::{Money, MoneyBag},
    parallel,
};
use rhodium_market::{
    fx::FxTodayPlan,
    market::Market,
    structure::{BlockId, ParamId},
};
use rhodium_pricing::{
    calibration::CalibratedMarket,
    prepare::{Prepared, evaluate},
};
use scioxide_ad::{Linearization, NonFinitePolicy, Tape, TapeOptions, Var};
use scioxide_core::{ExecutionContext, InvalidInputReason, SciOxideError, argument};

use crate::vector::RiskVector;

/// 通貨ごとの PV と疎な θ 勾配。θ の自然な単位あたりの金額で表す。
#[derive(Debug, Clone)]
pub struct ParamGradient {
    pv: MoneyBag<f64>,
    by_ccy: Vec<RiskVector<ParamId>>,
}

impl ParamGradient {
    /// 評価で得た通貨ごとの PV。
    #[must_use]
    pub const fn pv(&self) -> &MoneyBag<f64> {
        &self.pv
    }

    /// 通貨順に並んだ θ 勾配。
    #[must_use]
    pub fn by_ccy(&self) -> &[RiskVector<ParamId>] {
        &self.by_ccy
    }
}

/// 準備済み商品を評価し、通貨ごとの PV と θ 勾配を求める。
///
/// # Errors
///
/// 準備時と異なる市場構造、利用できない曲線、評価または AAD の失敗を返す。
pub fn param_gradient(prepared: &Prepared, cm: &CalibratedMarket) -> Result<ParamGradient> {
    with_tape(
        prepared,
        cm,
        &mut new_tape(),
        ExecutionContext::sequential(),
    )
}

/// 取引の評価に使う tape を作る。
///
/// 評価の途中で NaN や無限大が生じたら、勾配に黙って伝えずにその取引のエラーにする。
pub(crate) fn new_tape() -> Tape<f64> {
    Tape::new(TapeOptions::new().with_non_finite_policy(NonFinitePolicy::Error))
}

/// 呼び出し側の tape を使って [`param_gradient`] を計算する。
///
/// [`crate::delta::quote_delta`] はワーカーごとに一つの tape を取引間で使い回し、
/// 取引ごとに tape の領域を確保し直す費用を避ける。`execution` はワーカーの子の
/// 実行文脈で、記録と逆伝播の中では並列にしない。
///
/// 記録と逆伝播は `execution` の取消と期限を確かめる。止まった場合は取引の失敗ではなく
/// `Execution` を返し、呼び出し側が集計全体を止められるようにする。
pub(crate) fn with_tape(
    prepared: &Prepared,
    cm: &CalibratedMarket,
    tape: &mut Tape<f64>,
    execution: ExecutionContext<'_>,
) -> Result<ParamGradient> {
    let active = active_blocks(cm.market(), prepared.blocks())?;
    let currencies = RefCell::new(Vec::new());
    let mut evaluation = linearize_market(cm, tape, &active, execution, |market, outputs| {
        let pv = evaluate(prepared, market)?;
        // tape の容量拡張で再記録されても、通貨と出力の対応を一度だけ記録する。
        let mut currencies = currencies.borrow_mut();
        currencies.clear();
        for term in pv.terms() {
            currencies.push(term.ccy());
            outputs.push(term.amount());
        }
        Ok(())
    })?;

    let currencies = currencies.into_inner();
    let mut pv = MoneyBag::new();
    let mut by_ccy = Vec::with_capacity(currencies.len());
    let mut seed = vec![0.0; currencies.len()];
    let mut gradient = vec![0.0; evaluation.inputs()];
    for (row, ccy) in currencies.into_iter().enumerate() {
        pv.add(Money::new(ccy, evaluation.values()[row]));
        seed[row] = 1.0;
        evaluation
            .vjp_into(&seed, &mut gradient, execution)
            .map_err(|error| {
                parallel::control_failure(&error)
                    .unwrap_or_else(|| Error::Numerical(format!("AAD gradient for {ccy}: {error}")))
            })?;
        seed[row] = 0.0;

        let positions = active.iter().flat_map(|&(block, ref range)| {
            // 市場構造はブロック内の添字が u32 に収まることを保証する。
            (0..=u32::MAX).take(range.len()).map(move |i| (block, i))
        });
        let entries = positions
            .zip(gradient.iter().copied())
            .filter(|&(_, value)| value != 0.0)
            .map(|((block, i), value)| (ParamId::new(block.clone(), i), value))
            .collect();
        by_ccy.push(RiskVector::from_sorted(ccy, entries));
    }
    Ok(ParamGradient { pv, by_ccy })
}

/// 評価時点の FX と、各通貨の PV を固定した換算部分の θ 勾配を返す。
pub(crate) fn fx_translation(
    cm: &CalibratedMarket,
    plans: &[FxTodayPlan],
    pvs: &[f64],
    execution: ExecutionContext<'_>,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let blocks: BTreeSet<_> = plans.iter().flat_map(FxTodayPlan::blocks).collect();
    let active = active_blocks(cm.market(), &blocks)?;
    let mut tape = new_tape();
    let mut rates = linearize_market(cm, &mut tape, &active, execution, |market, outputs| {
        for plan in plans {
            outputs.push(plan.rate(market)?.rate());
        }
        Ok(())
    })?;
    // 通貨ごとの集計 PV を種にして、一度の逆伝播で Σ V_c ∂X_c/∂θ を求める。
    let mut gradient = vec![0.0; rates.inputs()];
    rates
        .vjp_into(pvs, &mut gradient, execution)
        .map_err(|error| {
            parallel::control_failure(&error)
                .unwrap_or_else(|| Error::Numerical(format!("AAD FX translation: {error}")))
        })?;
    let mut theta = vec![0.0; cm.market().theta().len()];
    for (index, value) in active
        .iter()
        .flat_map(|(_, range)| range.clone())
        .zip(gradient)
    {
        theta[index] = value;
    }
    Ok((rates.values().to_vec(), theta))
}

fn active_blocks<'a>(
    market: &Market<f64>,
    blocks: impl IntoIterator<Item = &'a BlockId>,
) -> Result<Vec<(&'a BlockId, Range<usize>)>> {
    let layout = market.structure().layout();
    blocks
        .into_iter()
        .filter(|block| market.block_is_available(block))
        .map(|block| {
            let range = layout.range(block).ok_or_else(|| {
                Error::Internal(format!("missing theta range for curve {}", block.curve()))
            })?;
            Ok((block, range))
        })
        .collect()
}

fn linearize_market<'a>(
    cm: &CalibratedMarket,
    tape: &'a mut Tape<f64>,
    active: &[(&BlockId, Range<usize>)],
    execution: ExecutionContext<'_>,
    evaluate: impl for<'t> Fn(&Market<Var<'t, f64>>, &mut Vec<Var<'t, f64>>) -> Result<()>,
) -> Result<Linearization<'a, f64>> {
    let market = cm.market();
    let point: Vec<_> = active
        .iter()
        .flat_map(|(_, range)| market.theta()[range.clone()].iter().copied())
        .collect();
    let failure = Cell::new(None);
    tape.linearize(
        |inputs, outputs| {
            // 元の市場の利用可否を保ったまま、読み出す θ だけを tape の変数にする。
            let mut input = inputs;
            let blocks = active.iter().map(|&(block, ref range)| {
                let (values, rest) = input.split_at(range.len());
                input = rest;
                (block, values)
            });
            market
                .with_active_theta(blocks)
                .and_then(|market| evaluate(&market, outputs))
                .map_err(|error| {
                    failure.set(Some(error));
                    SciOxideError::invalid_input(
                        argument!("theta"),
                        InvalidInputReason::Inconsistent,
                        "AAD market evaluation failed",
                    )
                })
        },
        &point,
        execution,
    )
    .map_err(|error| {
        failure.take().unwrap_or_else(|| {
            parallel::control_failure(&error)
                .unwrap_or_else(|| Error::Numerical(format!("AAD market evaluation: {error}")))
        })
    })
}
`````

`````
crates/rhodium-risk/src/lib.rs
//! 価格の感応度、シナリオと規制計算。

#![forbid(unsafe_code)]

mod adjoint;
pub mod delta;
pub mod gradient;
pub mod vector;
`````

`````
crates/rhodium-risk/src/vector.rs
//! 通貨と感応度のキー空間を区別する疎ベクトル。

use std::cmp::Ordering;

use rhodium_core::{Error, Result, currency::Currency};

/// キーごとの感応度。値は指定通貨の金額で、各キーの単位あたり。
#[derive(Debug, Clone)]
pub struct RiskVector<K> {
    ccy: Currency,
    entries: Vec<(K, f64)>,
}

impl<K> RiskVector<K> {
    /// 有限な値をキーの昇順で一つずつ受け取る。ゼロの要素は省く。
    pub(crate) fn from_sorted(ccy: Currency, mut entries: Vec<(K, f64)>) -> Self {
        entries.retain(|&(_, value)| value != 0.0);
        Self { ccy, entries }
    }

    /// 感応度の金額の通貨。
    #[must_use]
    pub const fn ccy(&self) -> Currency {
        self.ccy
    }

    /// 非ゼロの感応度をキーの昇順で返す。
    #[must_use]
    pub fn entries(&self) -> &[(K, f64)] {
        &self.entries
    }
}

impl<K: Ord + Clone> RiskVector<K> {
    /// 同じ通貨の感応度をキーごとに足す。
    ///
    /// # Errors
    ///
    /// 通貨が異なれば `InvalidInput`、合計が非有限なら `Numerical`。
    pub fn checked_add(&self, other: &Self) -> Result<Self> {
        if self.ccy != other.ccy {
            return Err(Error::InvalidInput(format!(
                "cannot add {} and {} sensitivities",
                self.ccy, other.ccy
            )));
        }

        let mut left = self.entries.as_slice();
        let mut right = other.entries.as_slice();
        let mut entries = Vec::with_capacity(left.len().max(right.len()));
        while let (Some((a, left_tail)), Some((b, right_tail))) =
            (left.split_first(), right.split_first())
        {
            match a.0.cmp(&b.0) {
                Ordering::Less => {
                    entries.push(a.clone());
                    left = left_tail;
                }
                Ordering::Greater => {
                    entries.push(b.clone());
                    right = right_tail;
                }
                Ordering::Equal => {
                    let sum = a.1 + b.1;
                    if !sum.is_finite() {
                        return Err(Error::Numerical(format!(
                            "non-finite {} sensitivity sum",
                            self.ccy
                        )));
                    }
                    entries.push((a.0.clone(), sum));
                    left = left_tail;
                    right = right_tail;
                }
            }
        }
        entries.extend_from_slice(left);
        entries.extend_from_slice(right);
        Ok(Self::from_sorted(self.ccy, entries))
    }
}
`````

`````
crates/rhodium-risk/tests/compile_fail.rs
//! 感応度のキー型が異なるベクトルを足せないことを確認する。

use rhodium_testkit::compile_fail_pair;

#[test]
fn parameter_gradient_cannot_be_added_to_quote_delta() {
    // θ と quote をそれぞれ正しい空間のベクトルとして扱う。
    compile_fail_pair!("quote_delta_with_param_gradient");
}

#[test]
fn tape_escape() {
    // AD の変数はブランド 't を持ち、記録の関数の外へ持ち出せない。
    compile_fail_pair!("tape_escape");
}
`````

`````
crates/rhodium-risk/tests/least_squares.rs
//! 最小二乗の校正グループを通る quote デルタを、手計算と再校正の中心差分で照合する。

use rhodium_core::{
    Error, Result,
    currency::{Asset, Currency},
    date::Date,
    money::Money,
};
use rhodium_market::{
    curve::{CurveId, CurveRole, CurveSpec, Extrapolation},
    quote::{QuoteId, QuoteKind, QuoteSet},
    structure::{BlockId, ParamId},
};
use rhodium_pricing::{
    calibration::{CalibratedMarket, GroupOutcome, calibrate, calibrate_with_residual_tolerance},
    context::Context,
    definition::{
        CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
    },
    prepare::{PrepareInput, Prepared, evaluate, prepare},
    resolve::CollateralKey,
};
use rhodium_products::{
    ir::{Contract, Product},
    standard::{Deposit, DepositDates, PaymentRounding, Position, StandardProduct},
};
use rhodium_risk::{delta::quote_delta, gradient::param_gradient, vector::RiskVector};
use rhodium_testkit::{assert_close, check_aad_vs_bump};
use scioxide_core::ExecutionContext;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、最小二乗の照合に使わない関数がある"
)]
#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
pub mod ois;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、最小二乗の照合に使わない関数がある"
)]
#[path = "../../rhodium-pricing/tests/fixtures/steps.rs"]
pub mod steps;

use ois::{Inputs, Node, cases, midnight, prepare_ois};
use steps::{PenalizedTona, TURN_QUOTE};

const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
const NOTIONAL: f64 = 1_000_000.0;
const BP: f64 = 1e-4;

fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
    vector
        .entries()
        .binary_search_by(|(key, _)| key.cmp(id))
        .ok()
        .map_or(0.0, |index| vector.entries()[index].1)
}

fn jpy(report: &rhodium_risk::delta::DeltaReport) -> Result<&RiskVector<QuoteId>> {
    report
        .by_ccy()
        .iter()
        .find(|vector| vector.ccy() == Currency::JPY)
        .ok_or_else(|| Error::Internal("no JPY quote deltas".into()))
}

#[test]
fn lsq_adjoint_oracle() -> Result<()> {
    // 感応度の検証表の例: min (θ − q₁)² + 4(θ − q₂)² の解 θ = (q₁ + 4q₂)/5 で、
    // quote への勾配 (1/5, 4/5) を校正の随伴が再現する。θ はピラーの ln DF なので、
    // ピラーで 100 万円を受け取る取引の V = N·exp(θ) の θ 勾配 V を、二つのマークに 1:4 で配る。
    let [tona, _] = cases()?;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let ctx = inputs.context();
    let id: CurveId = tona.index.parse()?;
    let pillar = midnight(tona.nodes[2].date)?;
    let curve = CurveSpec::new(
        id.clone(),
        Currency::JPY,
        [CurveRole::Discount, CurveRole::Projection],
        [pillar],
        Extrapolation::Reject,
    )?;
    let marks: [QuoteId; 2] = ["MARK.1".parse()?, "MARK.2".parse()?];
    let bindings = marks
        .iter()
        .zip([1.0, 0.5])
        .map(|(quote, sigma)| {
            InstrumentBinding::new(
                quote.clone(),
                QuoteKind::Mark,
                CalibrationTarget::DirectMark,
                QuoteTransform::Identity,
                COLLATERAL,
                ElementSelector::Param(ParamId::new(BlockId::nodes(id.clone()), 0)),
            )?
            .with_tolerance(sigma)
        })
        .collect::<Result<Vec<_>>>()?;
    let definition = MarketDefinition::new("JPY".parse()?, "Etc/UTC".parse()?, [curve], bindings)?;
    let quotes = QuoteSet::new([(marks[0].clone(), -0.01), (marks[1].clone(), -0.02)])?;
    let market = calibrate(&definition, &quotes, &ctx)?;
    let product = Product::new(
        Contract::build(|builder| {
            let amount = builder.constant(NOTIONAL);
            Ok(builder.pay(Asset::Cash(Currency::JPY), amount, pillar))
        })?,
        ctx.refdata,
    )?;
    let trade = prepare(
        PrepareInput::Ir(&product),
        &COLLATERAL,
        &[],
        market.market().structure(),
        &ctx,
    )?;
    let gradient = param_gradient(&trade, &market)?;
    let [theta_gradient] = gradient.by_ccy()[0].entries() else {
        return Err(Error::Internal("expected one theta gradient".into()));
    };
    // V = N·exp(θ) なので ∂V/∂θ = V。
    let pv = evaluate(&trade, market.market())?.get(Currency::JPY);
    assert_close!(theta_gradient.1, pv, 1e-12, 1e-6);
    let report = quote_delta(std::slice::from_ref(&trade), &market, SEQUENTIAL)?;
    let vector = jpy(&report)?;
    // マークの感応度の単位は自然な単位 1。
    assert_close!(sensitivity(vector, &marks[0]), pv / 5.0, 1e-10, 1e-6);
    assert_close!(sensitivity(vector, &marks[1]), 4.0 * pv / 5.0, 1e-10, 1e-6);
    Ok(())
}

/// バンプの再校正の許容誤差。最小二乗のグループは停留条件の残差を θ に戻した各行の動きで判定する。
const BUMP_TOLERANCE: f64 = 1e-15;

fn calibrated(
    definition: &MarketDefinition,
    quotes: &QuoteSet,
    ctx: &Context<'_>,
    tolerance: f64,
) -> Result<CalibratedMarket> {
    let market = calibrate_with_residual_tolerance(definition, quotes, ctx, tolerance)?;
    for group in market.groups() {
        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
            return Err(Error::Calibration(format!("{group:?}")));
        }
    }
    Ok(market)
}

fn bumped(quotes: &QuoteSet, id: &QuoteId, bump: f64) -> Result<QuoteSet> {
    QuoteSet::new(
        quotes
            .iter()
            .map(|(key, value)| (key.clone(), if key == id { value + bump } else { value })),
    )
}

/// 階段の区間をまたぐ 9M と、接続点をまたぐ 18M の OIS（quote から外した固定金利）と、
/// 年末をまたぐ預金。
fn trades(
    tona: &PenalizedTona,
    market: &CalibratedMarket,
    ctx: &Context<'_>,
) -> Result<Vec<Prepared>> {
    let node = |tenor: &str| -> Result<Node> {
        Ok(Node {
            tenor: tenor.parse()?,
            rate: 0.0,
            date: tona.case.as_of,
            df: 1.0,
        })
    };
    let mut prepared = [("9M", 0.003), ("18M", 0.0065)]
        .into_iter()
        .map(|(tenor, rate)| {
            prepare_ois(
                &tona.case,
                &node(tenor)?,
                rate,
                &COLLATERAL,
                market.market(),
                ctx,
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let start: Date = "2025-07-09".parse()?;
    let deposit = StandardProduct::Deposit(Deposit::new(
        Money::new(Currency::JPY, NOTIONAL),
        0.004,
        DepositDates::Confirmed {
            start,
            maturity: "2026-04-09".parse()?,
        },
        None,
        Some(tona.case.basis()),
        Position::Long,
        PaymentRounding::None,
    )?);
    prepared.push(prepare(
        PrepareInput::Standard {
            product: &deposit,
            trade_date: tona.case.as_of,
        },
        &COLLATERAL,
        &[],
        market.market().structure(),
        ctx,
    )?);
    Ok(prepared)
}

#[test]
fn penalized_meeting_steps_quote_delta_matches_bump() -> Result<()> {
    // R4 の網羅表の行: 罰則付きの最小二乗で当てた会合日の階段を通る quote デルタ（AAD と
    // 2 階の項を含む校正の随伴）が、quote ごとに全体を再校正した中心差分と一致する。
    // 残差が残る当てはめなので、随伴が残差の 2 階の項を落とすと一致しない（Gauss–Newton の記録）。
    let tona = PenalizedTona::new()?;
    let inputs = Inputs::new(&[&tona.case], tona.case.as_of, "")?;
    let ctx = inputs.context();
    let definition = tona.definition(false)?;
    let quotes = tona.quotes()?;
    let base = calibrated(&definition, &quotes, &ctx, BUMP_TOLERANCE)?;
    let ids: Vec<QuoteId> = tona
        .case
        .nodes
        .iter()
        .map(|node| tona.case.quote(node))
        .chain([TURN_QUOTE.parse()])
        .collect::<Result<_>>()?;
    let prepared = trades(&tona, &base, &ctx)?;
    for (index, trade) in prepared.iter().enumerate() {
        let report = quote_delta(std::slice::from_ref(trade), &base, SEQUENTIAL)?;
        let vector = jpy(&report)?;
        // 停留条件の誤差は各行を自然な単位（quote）で高々 BUMP_TOLERANCE だけ動かすので、
        // quote の感応度の和に掛けて再校正の PV の誤差にする。
        let scale = vector
            .entries()
            .iter()
            .map(|(_, delta)| delta.abs() / BP)
            .sum::<f64>();
        for id in &ids {
            check_aad_vs_bump(
                |h| {
                    let market =
                        calibrated(&definition, &bumped(&quotes, id, h)?, &ctx, BUMP_TOLERANCE)?;
                    let trade = &trades(&tona, &market, &ctx)?[index];
                    Ok::<_, Error>(evaluate(trade, market.market())?.get(Currency::JPY))
                },
                sensitivity(vector, id),
                BP,
                BUMP_TOLERANCE * scale,
                1e-6,
                &format!("trade {index} {id}"),
            )?;
        }
    }
    // 対照: Gauss–Newton の記録の随伴は残差の 2 階の項を落とすので、同じ取引のデルタが厳密な値から
    // 比べる許容誤差（相対 1e-10）よりはるかに大きくずれる。上の照合がこの差を検出できることを示す。
    let approximate = calibrated(&tona.definition(true)?, &quotes, &ctx, BUMP_TOLERANCE)?;
    let trade = &trades(&tona, &approximate, &ctx)?[0];
    let exact = quote_delta(std::slice::from_ref(&prepared[0]), &base, SEQUENTIAL)?;
    let gauss_newton = quote_delta(std::slice::from_ref(trade), &approximate, SEQUENTIAL)?;
    let deviation = ids
        .iter()
        .map(|id| {
            let lhs = sensitivity(jpy(&exact)?, id);
            let rhs = sensitivity(jpy(&gauss_newton)?, id);
            Ok::<_, Error>((lhs - rhs).abs() / lhs.abs().max(1e-6))
        })
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .fold(0.0_f64, f64::max);
    assert!(deviation > 1e-6, "{deviation:e}");
    Ok(())
}
`````

`````
crates/rhodium-risk/tests/m1_quote_delta.rs
//! AAD とバンプの網羅表の M1 の行: 商品の族と曲線の組み方ごとに、quote デルタ（AAD と校正の随伴）を
//! quote ごとに全体を再校正した中心差分と照合する。
//!
//! - R2・R5: OIS の割引とターム指標の予想（TONA へのベーシスの TIBOR 6M・3M、USD の独立したターム曲線）
//!   の上の IRS、FRA、ベーシススワップ。
//! - R1 と LGM のマーク: 金利先物（解析プライサー、凸性の補正）の先物価格と LGM のマークへの感応度。
//! - R3: 会合日の階段と年末の跳びを持つ TONA の上の OIS（lookback、observation shift、lockout、
//!   支払いの遅れ）と、年末をまたぐ預金。
//!
//! 再校正した市場は quote とバンプの組ごとに一度だけ作り、同じ市場で全取引を評価する。

use std::{
    cell::RefCell,
    collections::HashMap,
    num::{NonZeroU8, NonZeroU32},
};

use rhodium_core::{
    Error, Result,
    calendar::CalendarSpec,
    currency::Currency,
    date::{Date, DateInterval},
    day_count::DayCount,
    fixing::FixingStore,
    money::Money,
    observation::{
        Aggregation, CompoundingForm, ObservationDays, ObservationMethod, ObservationRounding,
        ObservationRule, SpreadTreatment,
    },
    quantity::{Rate, RateFuturesPrice, Simple, Spread},
    schedule::{
        ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination, ScheduleSpec,
    },
};
use rhodium_market::quote::{QuoteId, QuoteKind, QuoteSet};
use rhodium_pricing::{
    calibration::{CalibratedMarket, GroupOutcome, calibrate_with_residual_tolerance},
    context::Context,
    definition::MarketDefinition,
    prepare::{PrepareInput, Prepared, evaluate, prepare},
    resolve::CollateralKey,
};
use rhodium_products::standard::{
    Coupon, Deposit, DepositDates, Fra, Leg, NotionalExchange, PayReceive, PaymentRounding,
    Position, RateFuture, StandardProduct, Swap,
};
use rhodium_risk::{delta::quote_delta, vector::RiskVector};
use rhodium_testkit::check_aad_vs_bump;
use scioxide_core::ExecutionContext;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
)]
#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
pub mod ois;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
)]
#[path = "../../rhodium-pricing/tests/fixtures/multicurve.rs"]
pub mod multicurve;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
)]
#[path = "../../rhodium-pricing/tests/fixtures/steps.rs"]
pub mod steps;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、網羅表の照合に使わない関数がある"
)]
#[path = "../../rhodium-pricing/tests/fixtures/futures.rs"]
pub mod futures;

use multicurve::{MulticurveCase, all_multicurve_cases};
use ois::{Inputs, Node, prepare_ois};
use steps::PenalizedTona;

const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
const BP: f64 = 1e-4;
const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
/// バンプの再校正の残差の許容誤差（通常の 1e-14 の 1/100）。
const BUMP_TOLERANCE: f64 = 1e-16;

/// 市場から照合する取引を準備する関数。取引の名前と準備の組を返す。
type Trades<'a> = dyn Fn(&CalibratedMarket) -> Result<Vec<(String, Prepared)>> + 'a;

fn calibrated(
    definition: &MarketDefinition,
    quotes: &QuoteSet,
    ctx: &Context<'_>,
) -> Result<CalibratedMarket> {
    let market = calibrate_with_residual_tolerance(definition, quotes, ctx, BUMP_TOLERANCE)?;
    for group in market.groups() {
        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
            return Err(Error::Calibration(format!("{group:?}")));
        }
    }
    Ok(market)
}

fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
    vector
        .entries()
        .binary_search_by(|(key, _)| key.cmp(id))
        .ok()
        .map_or(0.0, |index| vector.entries()[index].1)
}

/// quote の報告単位（デルタの値の単位）と、バンプの幅。マークは報告単位 1 だが、1 だけ動かすと
/// 非線形の項が大きいので、1bp の幅で照合し、AAD を幅あたりに換算する。
fn units(market: &CalibratedMarket, id: &QuoteId) -> Result<(f64, f64)> {
    Ok(
        match market
            .quote_kind(id)
            .ok_or_else(|| Error::Internal(format!("unknown kind of {id}")))?
        {
            QuoteKind::Mark => (1.0, BP),
            _ => (BP, BP),
        },
    )
}

/// 全取引と全 quote で、AAD の quote デルタを再校正の中心差分と照合する。
fn check_all(
    definition: &MarketDefinition,
    quotes: &QuoteSet,
    ctx: &Context<'_>,
    ccy: Currency,
    floor: f64,
    trades: &Trades<'_>,
) -> Result<()> {
    let base = calibrated(definition, quotes, ctx)?;
    // 同じ (quote, 幅) の再校正は取引の間で共有する。
    let cache: RefCell<HashMap<(QuoteId, u64), CalibratedMarket>> = RefCell::new(HashMap::new());
    let reprice = |id: &QuoteId, bump: f64, index: usize| -> Result<f64> {
        let key = (id.clone(), bump.to_bits());
        if !cache.borrow().contains_key(&key) {
            let bumped =
                QuoteSet::new(quotes.iter().map(|(key, value)| {
                    (key.clone(), if key == id { value + bump } else { value })
                }))?;
            let market = calibrated(definition, &bumped, ctx)?;
            cache.borrow_mut().insert(key.clone(), market);
        }
        let cache = cache.borrow();
        let market = &cache[&key];
        let (_, prepared) = &trades(market)?[index];
        Ok(evaluate(prepared, market.market())?.get(ccy))
    };
    for (index, (name, prepared)) in trades(&base)?.iter().enumerate() {
        let report = quote_delta(std::slice::from_ref(prepared), &base, SEQUENTIAL)?;
        assert!(report.failures().is_empty(), "{name}");
        let vector = report
            .by_ccy()
            .iter()
            .find(|vector| vector.ccy() == ccy)
            .ok_or_else(|| Error::Internal(format!("{name}: no {ccy} deltas")))?;
        // 再校正の残差の許容誤差を、quote ごとの感応度（自然な単位）の和で PV の誤差に伝える。
        let scale = vector
            .entries()
            .iter()
            .map(|(id, delta)| Ok(delta.abs() / units(&base, id)?.0))
            .sum::<Result<f64>>()?;
        assert!(scale > 0.0, "{name} has no quote deltas");
        for (id, _) in quotes.iter() {
            let (report_unit, width) = units(&base, id)?;
            check_aad_vs_bump(
                |h| reprice(id, h, index),
                sensitivity(vector, id) * width / report_unit,
                width,
                BUMP_TOLERANCE * scale,
                floor,
                &format!("{name} {id}"),
            )?;
        }
    }
    Ok(())
}

fn standard(
    product: &StandardProduct,
    trade_date: Date,
    market: &CalibratedMarket,
    ctx: &Context<'_>,
) -> Result<Prepared> {
    prepare(
        PrepareInput::Standard {
            product,
            trade_date,
        },
        &COLLATERAL,
        &[],
        market.market().structure(),
        ctx,
    )
}

/// 固定 1% の IRS、ターム指標の FRA、quote から 10bp 外したベーシススワップ。元本は 1。
fn term_products(
    case: &MulticurveCase,
    ctx: &Context<'_>,
) -> Result<Vec<(String, StandardProduct)>> {
    let node = &case.nodes[2];
    let irs = rhodium_products::standard::InstrumentTemplate::Swap {
        convention: case.irs_convention().parse()?,
        tenor: node.tenor,
    }
    .build(0.01, ctx.refdata)?;
    let fra = StandardProduct::Fra(Fra::new(
        Money::new(case.overnight.currency, 1.0),
        case.term_index.parse()?,
        "2025-07-07".parse()?,
        0.01,
        Position::Long,
    )?);
    let basis = case.template(node)?.build(node.rate + 0.001, ctx.refdata)?;
    Ok(vec![
        (format!("{} IRS", case.term_index), irs),
        (format!("{} FRA", case.term_index), fra),
        (format!("{} basis", case.term_index), basis),
    ])
}

#[test]
fn term_products_quote_deltas_match_bumps_on_term_curves() -> Result<()> {
    // R2・R5: TIBOR 6M と 3M は TONA の曲線とベーシスの和、USD のターム曲線は独立の節点。
    // どれも OIS の曲線で割り引き、ターム曲線のリスクは随伴でベーシスと OIS の quote に戻る。
    for case in all_multicurve_cases()? {
        let inputs = case.inputs()?;
        let ctx = inputs.context();
        let products = term_products(&case, &ctx)?;
        let trade_date = case.overnight.as_of;
        let trades = |market: &CalibratedMarket| {
            products
                .iter()
                .map(|(name, product)| {
                    Ok((name.clone(), standard(product, trade_date, market, &ctx)?))
                })
                .collect()
        };
        check_all(
            &case.definition()?,
            &case.quotes()?,
            &ctx,
            case.overnight.currency,
            1e-14,
            &trades,
        )?;
    }
    Ok(())
}

#[test]
fn rate_future_quote_deltas_match_bumps() -> Result<()> {
    // 金利先物は解析プライサー（LGM の凸性の補正）で評価する。曲線を先物の価格で校正し、
    // 凸性に効く LGM の平均回帰とボラは直接マーク。先物価格への感応度は曲線の随伴で、
    // マークへの感応度は凸性の補正で、取引から quote へ戻る。
    let as_of = futures::AS_OF.parse()?;
    let reference = futures::refdata("0D".parse()?, "16:00:00")?;
    let config = futures::config("")?;
    let fixings = FixingStore::default();
    let ctx = Context {
        as_of,
        refdata: &reference,
        config: &config,
        fixings: &fixings,
    };
    let quotes = QuoteSet::new([
        ("USD.FUT.2026-01".parse()?, 0.96),
        ("USD.LGM.MR".parse()?, 0.03),
        ("USD.LGM.VOL".parse()?, 0.01),
    ])?;
    let future = StandardProduct::RateFuture(RateFuture::new(
        futures::CLEARING.parse()?,
        NonZeroU32::MIN,
        RateFuturesPrice::new(0.955, DayCount::Act365F)?,
        Position::Long,
        None,
    ));
    let trade_date = futures::START.parse()?;
    let trades = |market: &CalibratedMarket| {
        let prepared = standard(&future, trade_date, market, &ctx)?;
        Ok(vec![("rate future".to_owned(), prepared)])
    };
    check_all(
        &futures::futures_definition()?,
        &quotes,
        &ctx,
        Currency::USD,
        1e-8,
        &trades,
    )
}

/// 期間 2025-06-10〜2026-03-10（年末の跳びをまたぐ）の固定 1% 払い・TONA 複利受けの OIS で、
/// 観測の規約 `method` と、期間の終わりから 2 営業日後の支払い（支払いの遅れ）を持つ。
fn rfr_swap(case: &ois::Case, method: ObservationMethod) -> Result<StandardProduct> {
    let calendar = CalendarSpec::Base(case.currency.to_string().parse()?);
    let rule = ObservationRule::new(
        ObservationDays::BusinessDays(calendar),
        method,
        Aggregation::Compound {
            form: CompoundingForm::SimpleFactor,
            spread: SpreadTreatment::Exclusive,
        },
        ObservationRounding::default(),
    )?;
    let schedule = ScheduleSpec::Explicit(ExplicitPeriods::new(
        vec![ExplicitPeriod {
            accrual: DateInterval::new("2025-06-10".parse()?, "2026-03-10".parse()?)?,
            payment: "2026-03-12".parse()?,
            reset: None,
            settlement: None,
            references: Vec::new(),
        }],
        ExplicitDatePolicy::KeepAsConfirmed,
        PeriodCombination::Sum,
        None,
    )?);
    let notional = Money::new(case.currency, 1_000_000.0);
    let leg = |side, coupon| {
        Leg::new(
            side,
            notional,
            schedule.clone(),
            case.basis(),
            coupon,
            NotionalExchange::None,
            PaymentRounding::None,
        )
    };
    Ok(StandardProduct::Swap(Swap::from_legs([
        leg(
            PayReceive::Pay,
            Coupon::Fixed {
                rate: Rate::new(0.01, case.basis(), Simple)?,
            },
        )?,
        leg(
            PayReceive::Receive,
            Coupon::Compounded {
                index: case.index.parse()?,
                rule,
                spread: Spread::new(0.0)?,
            },
        )?,
    ])?))
}

#[test]
fn meeting_steps_and_turn_quote_deltas_match_bumps() -> Result<()> {
    // R3: 会合日の階段と年末の跳びの TONA。OIS の観測の規約 4 種（Plain は慣行の OIS）と支払いの遅れ、
    // 年末をまたぐ預金のデルタが、年末のスプレッドを含む全 quote の再校正と一致する。
    // lookback と lockout は日次積（DailyForwardProduct）で、observation shift は端点の比で評価する。
    let tona = PenalizedTona::new()?;
    let inputs = Inputs::new(&[&tona.case], tona.case.as_of, "")?;
    let ctx = inputs.context();
    let node = |tenor: &str| -> Result<Node> {
        Ok(Node {
            tenor: tenor.parse()?,
            rate: 0.0,
            date: tona.case.as_of,
            df: 1.0,
        })
    };
    let lag = |days: u8| NonZeroU8::new(days).ok_or_else(|| Error::Internal("zero lag".into()));
    let swaps = [
        ("lookback", ObservationMethod::Lookback(lag(2)?)),
        (
            "observation shift",
            ObservationMethod::ObservationShift(lag(2)?),
        ),
        ("lockout", ObservationMethod::Lockout(lag(1)?)),
    ]
    .into_iter()
    .map(|(name, method)| Ok((name.to_owned(), rfr_swap(&tona.case, method)?)))
    .collect::<Result<Vec<_>>>()?;
    let deposit = StandardProduct::Deposit(Deposit::new(
        Money::new(Currency::JPY, 1_000_000.0),
        0.004,
        DepositDates::Confirmed {
            start: "2025-07-09".parse()?,
            maturity: "2026-04-09".parse()?,
        },
        None,
        Some(tona.case.basis()),
        Position::Long,
        PaymentRounding::None,
    )?);
    let trade_date = tona.case.as_of;
    let trades = |market: &CalibratedMarket| {
        let mut prepared = [("9M", 0.003), ("18M", 0.0065)]
            .into_iter()
            .map(|(tenor, rate)| {
                Ok((
                    format!("OIS {tenor}"),
                    prepare_ois(
                        &tona.case,
                        &node(tenor)?,
                        rate,
                        &COLLATERAL,
                        market.market(),
                        &ctx,
                    )?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        for (name, product) in &swaps {
            prepared.push((name.clone(), standard(product, trade_date, market, &ctx)?));
        }
        prepared.push((
            "deposit".to_owned(),
            standard(&deposit, trade_date, market, &ctx)?,
        ));
        Ok(prepared)
    };
    check_all(
        &tona.square_definition()?,
        &tona.quotes()?,
        &ctx,
        Currency::JPY,
        1e-6,
        &trades,
    )
}
`````

`````
crates/rhodium-risk/tests/quote_delta.rs
//! M0 の市場と商品で、quote デルタを全体再校正の中心差分と照合する。

use std::{collections::BTreeSet, num::NonZeroUsize, sync::Arc};

use proptest::{
    prelude::*,
    test_runner::{Config, RngSeed, TestCaseError},
};
use rhodium_core::{
    Error, Result,
    currency::{Asset, Currency},
    date::{Date, DateInterval},
    fx::CcyPair,
    money::{Money, MoneyBag},
    parallel::{ExecutionContext, ThreadExecutor},
    quantity::{Rate, Simple, Spread},
    refdata::Convention,
    schedule::{
        EndRule, ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination,
        ScheduleSpec,
    },
};
use rhodium_market::{
    curve::{CurveId, Extrapolation},
    fx::FxCollateral,
    market::Market,
    quote::{QuoteId, QuoteSet},
};
use rhodium_pricing::{
    calibration::{CalibratedMarket, GroupOutcome, calibrate, calibrate_with_residual_tolerance},
    context::Context,
    definition::{DirectFx, FxDefinition, FxImpliedCurve, FxQuotes, MarketDefinition},
    prepare::{PrepareInput, Prepared, evaluate, prepare},
    resolve::CollateralKey,
};
use rhodium_products::{
    agreement::{CashCollateralTerms, Remuneration},
    ir::{Contract, Product, analysis::cashflows},
    standard::{
        Coupon, Deposit, DepositDates, Leg, NotionalExchange, PayReceive, PaymentRounding,
        Position, StandardProduct, Swap,
    },
};
use rhodium_risk::{
    delta::{DeltaReport, quote_delta, quote_delta_in},
    gradient::param_gradient,
    vector::RiskVector,
};
use rhodium_testkit::{TickClock, assert_close, check_aad_vs_bump};
use scioxide_core::{ClockBinding, Deadline};

#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
pub mod ois;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、このテストで使わない関数がある"
)]
#[path = "../../rhodium-pricing/tests/fixtures/multicurve.rs"]
pub mod multicurve;

use ois::{Case, Inputs, cases, definition, fixture_quotes, prepare_ois};

const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
/// 一つの取引のデルタは逐次に求める。
const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
const BP: f64 = 1e-4;
const NOTIONAL: f64 = 1_000_000.0;
// バンプの再校正の残差の許容誤差。品質の照合手順どおり、通常の校正の許容誤差（1e-14）の 1/100。
// 超過成長を対数から exp_m1 で求めるので、1M の OIS でも残差は約 1e-17 まで下がり、この値に届く。
// この残差を δ_V に伝える。δ_V / h₂ は 1 bp あたりの感応度の約 2e-12 倍で、比べる相対誤差 1e-10 より小さい。
const BUMP_TOLERANCE: f64 = 1e-16;

/// 両通貨の現金担保の割引を、それぞれの OIS 曲線に結ぶ価格設定。
const CASH_DISCOUNTING: &str = r#"
[[discounting]]
currency = "JPY"
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
curve = "JPY-TONA"

[[discounting]]
currency = "USD"
collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
curve = "USD-SOFR"
"#;

fn prepare_trade(
    product: &StandardProduct,
    date: Date,
    collateral: &CollateralKey,
    cm: &CalibratedMarket,
    ctx: &Context<'_>,
) -> Result<Prepared> {
    prepare(
        PrepareInput::Standard {
            product,
            trade_date: date,
        },
        collateral,
        &[],
        cm.market().structure(),
        ctx,
    )
}

fn cash(case: &Case) -> Result<CollateralKey> {
    Ok(CollateralKey::Cash(CashCollateralTerms {
        ccy: case.currency,
        remuneration: Remuneration {
            index: case.index.parse()?,
            spread: Spread::new(0.0)?,
            floor: None,
        },
    }))
}

/// 慣行（ID は指標と同じ）どおりの OIS。開始日を省けば約定日からスポットラグで決める。
fn convention_swap(
    case: &Case,
    start: Option<Date>,
    tenor: &str,
    fixed_rate: f64,
    side: PayReceive,
    spread: f64,
) -> Result<StandardProduct> {
    Ok(StandardProduct::Swap(Swap::from_convention(
        case.index.parse()?,
        start,
        EndRule::AfterStart(tenor.parse()?),
        Money::new(case.currency, NOTIONAL),
        fixed_rate,
        side,
        Spread::new(spread)?,
        PaymentRounding::None,
    )?))
}

fn deposit(
    case: &Case,
    start: Date,
    maturity: Date,
    position: Position,
) -> Result<StandardProduct> {
    Ok(StandardProduct::Deposit(Deposit::new(
        Money::new(case.currency, NOTIONAL),
        case.nodes[0].rate + 0.003,
        DepositDates::Confirmed { start, maturity },
        None,
        Some(case.basis()),
        position,
        PaymentRounding::None,
    )?))
}

fn fixed_rfr_swap(case: &Case, ctx: &Context<'_>) -> Result<StandardProduct> {
    let Convention::Ois(convention) = ctx.refdata.convention(&case.index.parse()?)? else {
        return Err(Error::Internal("expected OIS test convention".into()));
    };
    let date = "2026-07-09".parse()?;
    let schedule = ScheduleSpec::Explicit(ExplicitPeriods::new(
        vec![ExplicitPeriod {
            accrual: DateInterval::new("2025-07-09".parse()?, date)?,
            payment: date,
            reset: None,
            settlement: None,
            references: Vec::new(),
        }],
        ExplicitDatePolicy::KeepAsConfirmed,
        PeriodCombination::Sum,
        None,
    )?);
    let leg = |side, coupon| {
        Leg::new(
            side,
            Money::new(case.currency, NOTIONAL),
            schedule.clone(),
            case.basis(),
            coupon,
            NotionalExchange::None,
            PaymentRounding::None,
        )
    };
    Ok(StandardProduct::Swap(Swap::from_legs([
        leg(
            PayReceive::Receive,
            Coupon::Fixed {
                rate: Rate::new(case.nodes[1].rate + 0.004, case.basis(), Simple)?,
            },
        )?,
        leg(
            PayReceive::Pay,
            Coupon::Compounded {
                index: convention.float.index.clone(),
                rule: convention.float.observation.clone(),
                spread: Spread::new(0.001)?,
            },
        )?,
    ])?))
}

fn calibrated(
    def: &MarketDefinition,
    quotes: &QuoteSet,
    ctx: &Context<'_>,
    tolerance: f64,
) -> Result<CalibratedMarket> {
    let cm = calibrate_with_residual_tolerance(def, quotes, ctx, tolerance)?;
    for group in cm.groups() {
        let GroupOutcome::Solved(solution) = group.outcome() else {
            return Err(Error::Calibration(format!("risk test fixture: {group:?}")));
        };
        assert!(solution.residual_norm() <= tolerance);
    }
    Ok(cm)
}

fn bumped(quotes: &QuoteSet, id: &QuoteId, bump: f64) -> Result<QuoteSet> {
    QuoteSet::new(
        quotes
            .iter()
            .map(|(key, value)| (key.clone(), if key == id { value + bump } else { value })),
    )
}

fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
    vector
        .entries()
        .binary_search_by(|(key, _)| key.cmp(id))
        .ok()
        .map_or(0.0, |index| vector.entries()[index].1)
}

fn vector_for(vectors: &[RiskVector<QuoteId>], ccy: Currency) -> Result<&RiskVector<QuoteId>> {
    vectors
        .iter()
        .find(|vector| vector.ccy() == ccy)
        .ok_or_else(|| Error::Internal(format!("no quote deltas for {ccy}")))
}

fn absolute_tolerance(ctx: &Context<'_>, ccy: Currency) -> Result<f64> {
    let digits = ctx.refdata.currency(ccy)?.minor_units;
    Ok(1e-6 * 10.0_f64.powi(-i32::from(digits)))
}

/// バンプで照合する取引の名前、約定日、担保、商品。
type BumpTrade = (&'static str, Date, CollateralKey, StandardProduct);

// 一つの通貨で照合する 5 種類の取引を作る。
fn bump_trades(case: &Case, ctx: &Context<'_>) -> Result<Vec<BumpTrade>> {
    let seasoned = case
        .swaps
        .iter()
        .find(|swap| swap.name == "seasoned")
        .ok_or_else(|| Error::Internal("missing seasoned swap fixture".into()))?;
    let deposit = deposit(
        case,
        "2025-07-09".parse()?,
        "2026-07-09".parse()?,
        Position::Long,
    )?;
    let six_months = case.nodes[1].rate + 0.003;
    let spot = convention_swap(case, None, "6M", six_months, PayReceive::Receive, 0.0)?;
    let started = convention_swap(
        case,
        Some(seasoned.start),
        "1Y",
        seasoned.fixed_rate,
        PayReceive::Pay,
        seasoned.spread,
    )?;
    let cash_trade = convention_swap(case, None, "18M", case.nodes[3].rate, PayReceive::Pay, 0.0)?;
    Ok(vec![
        ("deposit", case.as_of, COLLATERAL, deposit),
        ("OIS", case.as_of, COLLATERAL, spot),
        (
            "fixed/RFR swap",
            case.as_of,
            COLLATERAL,
            fixed_rfr_swap(case, ctx)?,
        ),
        ("seasoned OIS", seasoned.start, COLLATERAL, started),
        (
            "cash-collateralized OIS",
            case.as_of,
            cash(case)?,
            cash_trade,
        ),
    ])
}

#[test]
fn quote_delta_vs_bump() -> Result<()> {
    // R1 の両通貨（QuantLib と同じ 1M〜2Y の 5 ピラー）で、預金、慣行どおりの OIS、明示日程の
    // 固定対 RFR スワップ、fixing とスプレッドを持つ開始済みの OIS、現金担保で割り引く OIS を照合する。
    // 独立した通貨の quote も動かし、他通貨の PV に感応度が出ないことを確認する。
    let [tona, sofr] = cases()?;
    let inputs = Inputs::new(&[&tona, &sofr], tona.as_of, CASH_DISCOUNTING)?;
    let ctx = inputs.context();
    let def = definition(&[(&tona, &COLLATERAL), (&sofr, &COLLATERAL)])?;
    let quotes = QuoteSet::new(
        fixture_quotes(&tona)?
            .into_iter()
            .chain(fixture_quotes(&sofr)?),
    )?;
    let base = calibrate(&def, &quotes, &ctx)?;

    for case in [&tona, &sofr] {
        let absolute = absolute_tolerance(&ctx, case.currency)?;
        let block = case.block()?;
        for (name, trade_date, collateral, product) in bump_trades(case, &ctx)? {
            let prepare_one =
                |cm: &CalibratedMarket| prepare_trade(&product, trade_date, &collateral, cm, &ctx);
            let plan = prepare_one(&base)?;
            let gradient = param_gradient(&plan, &base)?;
            assert_close!(
                gradient.pv().get(case.currency),
                evaluate(&plan, base.market())?.get(case.currency),
                1e-12,
                absolute
            );
            let [active] = gradient.by_ccy() else {
                return Err(Error::Internal("expected one currency in gradient".into()));
            };
            assert_ne!(active.entries(), []);
            assert!(active.entries().iter().all(|(id, _)| id.block() == &block));

            let delta = quote_delta(&[plan], &base, SEQUENTIAL)?;
            assert!(delta.failures().is_empty());
            let vector = vector_for(delta.by_ccy(), case.currency)?;
            assert_ne!(vector.entries(), []);
            let scale: f64 = vector
                .entries()
                .iter()
                .map(|(_, value)| value.abs() / BP)
                .sum();
            for (id, _) in quotes.iter() {
                check_aad_vs_bump(
                    |bump| {
                        let cm =
                            calibrated(&def, &bumped(&quotes, id, bump)?, &ctx, BUMP_TOLERANCE)?;
                        Ok::<_, Error>(
                            evaluate(&prepare_one(&cm)?, cm.market())?.get(case.currency),
                        )
                    },
                    sensitivity(vector, id),
                    BP,
                    BUMP_TOLERANCE * scale,
                    absolute,
                    &format!("{name} {} {id}", case.currency),
                )?;
            }
        }
    }
    Ok(())
}

#[test]
fn rounded_payments_report_the_delta_of_the_unrounded_amount() -> Result<()> {
    // 受渡額を円単位で丸める元本 100 円の 6M OIS は、quote を ±1bp 動かしても丸めた PV が変わらない
    // （中心差分は 0）。丸めの微分を 1 とする設計どおり、AAD は丸めた額 round(x) の微分を dx とし、
    // ∂PV = DF·∂x + round(x₀)·∂DF を返す。報告する感応度の意味に合わせ、丸め付きの取引は、
    // 基準点の丸め差 c = round(x₀) − x₀ を固定した再評価 PV_c(θ) = Σ (x(θ) + c)·DF(θ) と照合する
    // （品質の「AAD とバンプの照合」）。両脚の受渡は同じ支払日 t なので、
    // PV_c(θ) = PV_丸めなし(θ) + c·DF_t(θ) で、DF_t は t に 1 円を払う商品の PV から求める。
    let [tona, _] = cases()?;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let ctx = inputs.context();
    let def = definition(&[(&tona, &COLLATERAL)])?;
    let quotes = QuoteSet::new(fixture_quotes(&tona)?)?;
    let base = calibrate(&def, &quotes, &ctx)?;
    let swap = |rounding| -> Result<StandardProduct> {
        Ok(StandardProduct::Swap(Swap::from_convention(
            tona.index.parse()?,
            None,
            EndRule::AfterStart("6M".parse()?),
            Money::new(tona.currency, 100.0),
            tona.nodes[1].rate + 0.003,
            PayReceive::Receive,
            Spread::new(0.0)?,
            rounding,
        )?))
    };
    let rounded = swap(PaymentRounding::CurrencyConvention)?;
    let unrounded = swap(PaymentRounding::None)?;
    let ir = unrounded.to_product(tona.as_of, ctx.refdata)?;
    let pays: BTreeSet<_> = cashflows(&ir).iter().map(|cashflow| cashflow.at).collect();
    let [pay] = pays.into_iter().collect::<Vec<_>>()[..] else {
        return Err(Error::Internal("expected one payment instant".into()));
    };
    let unit = Product::new(
        Contract::build(|builder| {
            let amount = builder.constant(1.0);
            Ok(builder.pay(Asset::Cash(tona.currency), amount, pay))
        })?,
        ctx.refdata,
    )?;
    let pv = |product: &StandardProduct, cm: &CalibratedMarket| -> Result<f64> {
        let plan = prepare_trade(product, tona.as_of, &COLLATERAL, cm, &ctx)?;
        Ok(evaluate(&plan, cm.market())?.get(tona.currency))
    };
    let discount = |cm: &CalibratedMarket| -> Result<f64> {
        let plan = prepare(
            PrepareInput::Ir(&unit),
            &COLLATERAL,
            &[],
            cm.market().structure(),
            &ctx,
        )?;
        Ok(evaluate(&plan, cm.market())?.get(tona.currency))
    };
    let difference = (pv(&rounded, &base)? - pv(&unrounded, &base)?) / discount(&base)?;
    // 丸めは実際に働いている（基準点で 1 円未満の差がある）。
    assert!(difference != 0.0 && difference.abs() < 1.0, "{difference}");

    let plan = prepare_trade(&rounded, tona.as_of, &COLLATERAL, &base, &ctx)?;
    let report = quote_delta(&[plan], &base, SEQUENTIAL)?;
    assert!(report.failures().is_empty());
    let vector = vector_for(report.by_ccy(), tona.currency)?;
    let absolute = absolute_tolerance(&ctx, tona.currency)?;
    let scale: f64 = vector
        .entries()
        .iter()
        .map(|(_, value)| value.abs() / BP)
        .sum();
    // 丸めた PV は 6M の quote を動かしても変わらないが、AAD は 0 でない。
    let six_months = tona.quote(&tona.nodes[1])?;
    let rounded_pv = |bump: f64| -> Result<f64> {
        let cm = calibrate(&def, &bumped(&quotes, &six_months, bump)?, &ctx)?;
        pv(&rounded, &cm)
    };
    assert_eq!(rounded_pv(BP)?, rounded_pv(-BP)?);
    assert!(sensitivity(vector, &six_months).abs() > absolute);
    for (id, _) in quotes.iter() {
        check_aad_vs_bump(
            |bump| {
                let cm = calibrated(&def, &bumped(&quotes, id, bump)?, &ctx, BUMP_TOLERANCE)?;
                Ok::<_, Error>(pv(&unrounded, &cm)? + difference * discount(&cm)?)
            },
            sensitivity(vector, id),
            BP,
            BUMP_TOLERANCE * scale,
            absolute,
            &format!("rounded OIS {id}"),
        )?;
    }
    Ok(())
}

#[test]
fn short_tenor_quote_deltas_match_bumps() -> Result<()> {
    // 1D・1W・1M の TONA の OIS だけの曲線で、quote が 0.5% 付近と 0 付近の場合に、
    // 短い OIS の quote デルタ（AAD と校正の随伴）がバンプと再校正の値と一致することを確かめる。
    // 超過成長を丸めた DF 比から 1 を引いて求めると、1D の校正が 1e-14 に届かず、
    // 1/100 の再校正の許容誤差でも解けない（短期の桁落ちの回帰テスト）。
    let [tona, _] = cases()?;
    let node = |tenor: &str, rate: f64, date: &str| -> Result<ois::Node> {
        Ok(ois::Node {
            tenor: tenor.parse()?,
            rate,
            date: date.parse()?,
            df: f64::NAN,
        })
    };
    for rates in [[0.005, 0.0051, 0.0052], [1e-7, 0.0, -1e-7]] {
        // 評価日 2025-01-07（火）のスポットは 01-09 で、ピラーは各 OIS の最終支払日。
        let mut case = tona.clone();
        case.nodes = vec![
            node("1D", rates[0], "2025-01-10")?,
            node("1W", rates[1], "2025-01-16")?,
            node("1M", rates[2], "2025-02-10")?,
        ];
        let inputs = Inputs::new(&[&case], case.as_of, "")?;
        let ctx = inputs.context();
        let def = definition(&[(&case, &COLLATERAL)])?;
        let quotes = QuoteSet::new(fixture_quotes(&case)?)?;
        let base = calibrated(&def, &quotes, &ctx, 1e-14)?;
        let absolute = absolute_tolerance(&ctx, case.currency)?;
        // ピラーちょうどの 1W と、ピラーの間で終わる 2W（スプレッド付き）の OIS。
        for (name, product) in [
            (
                "1W OIS",
                convention_swap(
                    &case,
                    None,
                    "1W",
                    rates[1] + 0.001,
                    PayReceive::Receive,
                    0.0,
                )?,
            ),
            (
                "2W OIS",
                convention_swap(&case, None, "2W", rates[2], PayReceive::Pay, 0.0005)?,
            ),
        ] {
            let prepare_one =
                |cm: &CalibratedMarket| prepare_trade(&product, case.as_of, &COLLATERAL, cm, &ctx);
            let delta = quote_delta(&[prepare_one(&base)?], &base, SEQUENTIAL)?;
            assert!(delta.failures().is_empty());
            let vector = vector_for(delta.by_ccy(), case.currency)?;
            assert_ne!(vector.entries(), []);
            let scale: f64 = vector
                .entries()
                .iter()
                .map(|(_, value)| value.abs() / BP)
                .sum();
            for (id, _) in quotes.iter() {
                check_aad_vs_bump(
                    |bump| {
                        let cm =
                            calibrated(&def, &bumped(&quotes, id, bump)?, &ctx, BUMP_TOLERANCE)?;
                        Ok::<_, Error>(
                            evaluate(&prepare_one(&cm)?, cm.market())?.get(case.currency),
                        )
                    },
                    sensitivity(vector, id),
                    BP,
                    BUMP_TOLERANCE * scale,
                    absolute,
                    &format!("{name} {rates:?} {id}"),
                )?;
            }
        }
    }
    Ok(())
}

// 読むたびに進む時計の `reads` 回の読み出しの後に期限を置いて、quote デルタを求める。
fn delta_with_deadline(
    trades: &[Prepared],
    cm: &CalibratedMarket,
    executor: Option<&ThreadExecutor>,
    reads: u64,
) -> Result<DeltaReport> {
    let clock = TickClock::new();
    let binding = ClockBinding::new(&clock);
    let execution = match executor {
        Some(executor) => ExecutionContext::new(executor),
        None => ExecutionContext::sequential(),
    }
    .with_deadline(Deadline::at(&binding, TickClock::after(reads)))
    .map_err(|error| Error::Internal(error.to_string()))?;
    quote_delta(trades, cm, execution)
}

#[test]
fn quote_delta_stops_instead_of_dropping_trades_when_the_deadline_passes() -> Result<()> {
    // 実行器は区間の開始前にだけ期限を確かめ、区間の中では取引の境界と AD の記録・逆伝播が確かめる。
    // 区間の開始後に期限が来ても、止まった取引を「失敗した取引」として除いた結果を成功として
    // 返してはいけない。期限を 1 回目から全体の読み出しの回数まで 1 回ずつずらし、
    // 結果が「全件を集計した成功」か Execution のどちらかだけであることを確かめる。
    let [tona, _] = cases()?;
    let inputs = Inputs::new(&[&tona], tona.as_of, CASH_DISCOUNTING)?;
    let ctx = inputs.context();
    let cm = calibrate(
        &definition(&[(&tona, &COLLATERAL)])?,
        &QuoteSet::new(fixture_quotes(&tona)?)?,
        &ctx,
    )?;
    let trades = bump_trades(&tona, &ctx)?
        .into_iter()
        .map(|(_, date, collateral, product)| prepare_trade(&product, date, &collateral, &cm, &ctx))
        .collect::<Result<Vec<_>>>()?;
    let full = quote_delta(&trades, &cm, SEQUENTIAL)?;
    assert!(full.failures().is_empty());
    // 期限の無い時計で、逐次の全体が何回時計を読むかを数える。
    let clock = TickClock::new();
    let binding = ClockBinding::new(&clock);
    quote_delta(
        &trades,
        &cm,
        SEQUENTIAL
            .with_deadline(Deadline::at(&binding, TickClock::after(u64::MAX)))
            .map_err(|error| Error::Internal(error.to_string()))?,
    )?;
    let total = clock.reads();
    assert!(total > trades.len() as u64);

    let executor = ThreadExecutor::new(NonZeroUsize::MIN.saturating_add(1));
    let mut stopped = 0;
    for reads in 1..=total {
        for executor in [None, Some(&executor)] {
            match delta_with_deadline(&trades, &cm, executor, reads) {
                Err(Error::Execution(message)) => {
                    assert!(message.contains("deadline"), "{message}");
                    stopped += 1;
                }
                Ok(report) => {
                    // 期限より前に終われば、全件を集計した結果と同じになる。
                    assert!(
                        report.failures().is_empty(),
                        "{reads}: {:?}",
                        report.failures()
                    );
                    assert_eq!(report.by_ccy().len(), full.by_ccy().len());
                }
                Err(error) => return Err(error),
            }
        }
    }
    // 逐次では全体の読み出しの回数より前の期限はすべて止まる。
    assert!(stopped >= total);
    Ok(())
}

/// 線形性の性質テストで作る取引。通貨、種類、向き、固定金利のずらし幅を持つ。
#[derive(Debug, Clone, Copy)]
struct TradeSpec {
    usd: bool,
    kind: u8,
    long: bool,
    offset_bp: i32,
}

fn trade_spec() -> impl Strategy<Value = TradeSpec> {
    (any::<bool>(), 0_u8..3, any::<bool>(), -50_i32..=50).prop_map(
        |(usd, kind, long, offset_bp)| TradeSpec {
            usd,
            kind,
            long,
            offset_bp,
        },
    )
}

// 取引ごとに quote デルタを求めて足した値と、通貨別に θ 勾配を足してから一度だけ随伴を解いた値を比べる。
// 随伴は θ 勾配について線形なので、集計の順によらず許容誤差で一致する。
fn check_linearity(
    specs: &[TradeSpec],
    cases: &[Case; 2],
    cm: &CalibratedMarket,
    ctx: &Context<'_>,
) -> Result<()> {
    let trades = specs
        .iter()
        .map(|spec| {
            let case = &cases[usize::from(spec.usd)];
            let rate = case.nodes[1].rate + f64::from(spec.offset_bp) * BP;
            let side = if spec.long {
                PayReceive::Receive
            } else {
                PayReceive::Pay
            };
            let product = match spec.kind {
                0 => deposit(
                    case,
                    "2025-07-09".parse()?,
                    "2026-07-09".parse()?,
                    if spec.long {
                        Position::Long
                    } else {
                        Position::Short
                    },
                )?,
                1 => convention_swap(case, None, "1Y", rate, side, 0.0)?,
                _ => convention_swap(case, None, "2Y", rate, side, 0.001)?,
            };
            prepare_trade(&product, case.as_of, &COLLATERAL, cm, ctx)
        })
        .collect::<Result<Vec<_>>>()?;
    // 合計は 3 ワーカーで区間ごとの和を足し、単独の値は逐次で求める。
    let executor = ThreadExecutor::new(NonZeroUsize::MIN.saturating_add(2));
    let combined = quote_delta(&trades, cm, ExecutionContext::new(&executor))?;
    assert!(combined.failures().is_empty());
    let singles = trades
        .iter()
        .map(|plan| quote_delta(std::slice::from_ref(plan), cm, SEQUENTIAL))
        .collect::<Result<Vec<_>>>()?;
    // 単独の結果も除外ゼロでなければ、合計との比較が失敗の取引を見逃す。
    assert!(singles.iter().all(|single| single.failures().is_empty()));
    for ccy in [Currency::JPY, Currency::USD] {
        let parts: Vec<_> = singles
            .iter()
            .filter_map(|part| part.by_ccy().iter().find(|vector| vector.ccy() == ccy))
            .collect();
        let total = combined.by_ccy().iter().find(|vector| vector.ccy() == ccy);
        assert_eq!(total.is_some(), !parts.is_empty());
        let Some(total) = total else { continue };
        let floor = absolute_tolerance(ctx, ccy)?;
        // 合計と単独の結果のどちらかに現れる quote をすべて比べる。片方にしか無い quote は
        // もう片方を 0 とみなすので、合計から非ゼロの quote が抜けても見逃さない。
        let ids: BTreeSet<&QuoteId> = total
            .entries()
            .iter()
            .chain(parts.iter().flat_map(|part| part.entries()))
            .map(|(id, _)| id)
            .collect();
        for id in ids {
            let sum: f64 = parts.iter().map(|part| sensitivity(part, id)).sum();
            let gross: f64 = parts.iter().map(|part| sensitivity(part, id).abs()).sum();
            assert_close!(sensitivity(total, id), sum, 1e-10, floor, gross);
        }
    }
    // 通貨の違う感応度は足せない。
    if let [Some(jpy), Some(usd)] = [Currency::JPY, Currency::USD]
        .map(|ccy| combined.by_ccy().iter().find(|vector| vector.ccy() == ccy))
    {
        assert!(matches!(jpy.checked_add(usd), Err(Error::InvalidInput(_))));
    }
    Ok(())
}

proptest! {
    #![proptest_config(Config { cases: 24, rng_seed: RngSeed::Fixed(28), ..Config::default() })]

    #[test]
    fn adjoint_linearity(specs in prop::collection::vec(trade_spec(), 1..6)) {
        // 基準の市場は一度だけ校正し、乱数で選んだ取引の組で集計の線形性を確かめる。
        let fixture = || -> Result<()> {
            let cases = cases()?;
            let inputs = Inputs::new(&[&cases[0], &cases[1]], cases[0].as_of, "")?;
            let ctx = inputs.context();
            let def = definition(&[(&cases[0], &COLLATERAL), (&cases[1], &COLLATERAL)])?;
            let quotes = QuoteSet::new(
                fixture_quotes(&cases[0])?
                    .into_iter()
                    .chain(fixture_quotes(&cases[1])?),
            )?;
            let cm = calibrate(&def, &quotes, &ctx)?;
            check_linearity(&specs, &cases, &cm, &ctx)
        };
        fixture().map_err(|error| TestCaseError::fail(error.to_string()))?;
    }
}

#[test]
fn quote_delta_propagates_to_upstream_quotes() -> Result<()> {
    // 2 本の OIS 曲線で、予想曲線が割引曲線の変化を受ける場合も再校正と一致する。
    let [tona, _] = cases()?;
    let mut forecast = tona.clone();
    forecast.index = "JPY-FORECAST".into();
    for node in &mut forecast.nodes {
        node.rate += 0.002;
    }
    let inputs = Inputs::new(&[&tona, &forecast], tona.as_of, "")?;
    let ctx = inputs.context();
    let def = definition(&[(&forecast, &COLLATERAL), (&tona, &COLLATERAL)])?;
    let quotes = QuoteSet::new(
        fixture_quotes(&tona)?
            .into_iter()
            .chain(fixture_quotes(&forecast)?),
    )?;
    // 基準は通常の許容誤差で解き、上流と下流の両グループが解けたことを確かめる。
    // バンプの再校正は他のテストと同じく、校正の許容誤差の 1/100 の BUMP_TOLERANCE で解く。
    let cm = calibrated(&def, &quotes, &ctx, 1e-14)?;
    assert_eq!(cm.groups()[1].upstream(), &[tona.block()?]);
    let node = &forecast.nodes[1];
    let prepare_one = |cm: &CalibratedMarket| {
        prepare_ois(
            &forecast,
            node,
            node.rate + 0.003,
            &COLLATERAL,
            cm.market(),
            &ctx,
        )
    };
    let report = quote_delta(&[prepare_one(&cm)?], &cm, SEQUENTIAL)?;
    assert!(report.failures().is_empty());
    let vector = vector_for(report.by_ccy(), Currency::JPY)?;
    let scale: f64 = vector
        .entries()
        .iter()
        .map(|(_, value)| value.abs() / BP)
        .sum();
    for (id, _) in fixture_quotes(&tona)? {
        check_aad_vs_bump(
            |bump| {
                let cm = calibrated(&def, &bumped(&quotes, &id, bump)?, &ctx, BUMP_TOLERANCE)?;
                Ok::<_, Error>(evaluate(&prepare_one(&cm)?, cm.market())?.get(Currency::JPY))
            },
            sensitivity(vector, &id),
            BP,
            BUMP_TOLERANCE * scale,
            absolute_tolerance(&ctx, Currency::JPY)?,
            id.as_str(),
        )?;
    }
    Ok(())
}

#[test]
fn tibor_basis_quote_delta_matches_full_recalibration() -> Result<()> {
    // 市場 quote と異なる spread の取引で、ベーシスと TONA の全 quote を再校正で照合する。
    let [jpy, _] = multicurve::multicurve_cases()?;
    let inputs = jpy.inputs()?;
    let ctx = inputs.context();
    let definition = jpy.definition()?;
    let quotes = jpy.quotes()?;
    let base = calibrated(&definition, &quotes, &ctx, 1e-14)?;
    let node = &jpy.nodes[2];
    let trade = |market: &Market<f64>| jpy.prepare_basis(node, node.rate + 0.0025, market, &ctx);
    let report = quote_delta(&[trade(base.market())?], &base, SEQUENTIAL)?;
    assert!(report.failures().is_empty());
    let vector = vector_for(report.by_ccy(), Currency::JPY)?;
    let scale = vector
        .entries()
        .iter()
        .map(|(_, delta)| delta.abs() / BP)
        .sum::<f64>();
    // テンプレートの元本は 1。元本 100 万の円取引の許容誤差を元本 1 に換算する。
    let floor = absolute_tolerance(&ctx, Currency::JPY)? / NOTIONAL;
    let reprice = |market: &Market<f64>| {
        Ok::<_, Error>(evaluate(&trade(market)?, market)?.get(Currency::JPY))
    };
    for (id, _) in quotes.iter() {
        check_aad_vs_bump(
            |h| {
                let market =
                    calibrated(&definition, &bumped(&quotes, id, h)?, &ctx, BUMP_TOLERANCE)?;
                reprice(market.market())
            },
            sensitivity(vector, id),
            BP,
            BUMP_TOLERANCE * scale,
            floor,
            id.as_str(),
        )?;
    }
    assert_ne!(sensitivity(vector, &jpy.quote(node)?), 0.0);

    // TONA を動かしたとき、基準のベーシス θ を固定すると下流の再校正を経由する寄与が消える。
    let basis = base
        .market()
        .structure()
        .layout()
        .range(&jpy.term_block()?)
        .ok_or_else(|| Error::Internal("missing TIBOR basis block".into()))?;
    let mut max_indirect: f64 = 0.0;
    for (id, _) in fixture_quotes(&jpy.overnight)? {
        let direct_pv = |h| {
            let bumped = calibrated(&definition, &bumped(&quotes, &id, h)?, &ctx, BUMP_TOLERANCE)?;
            let mut theta = bumped.market().theta().to_vec();
            theta[basis.clone()].copy_from_slice(&base.market().theta()[basis.clone()]);
            let market = Market::from_theta(Arc::clone(bumped.market().structure()), theta)?;
            reprice(&market)
        };
        let direct_delta = direct_pv(BP / 2.0)? - direct_pv(-BP / 2.0)?;
        max_indirect = max_indirect.max((sensitivity(vector, &id) - direct_delta).abs());
    }
    // 比較の単位は元本 1。丸めや中心差分の打ち切りを超える間接寄与がある。
    assert!(
        max_indirect > 100.0 * floor,
        "indirect TONA delta: {max_indirect:e}"
    );
    Ok(())
}

fn usd_jpy(tona: &Case, sofr: &Case, implied: Option<FxImpliedCurve>) -> Result<FxDefinition> {
    let collateral: CurveId = "usd-jpy".parse()?;
    let base = implied
        .as_ref()
        .map_or_else(|| sofr.index.parse(), |curve| Ok(curve.id.clone()))?;
    let short = |id: &str| implied.as_ref().map(|_| id.parse()).transpose();
    Ok(FxDefinition {
        pair: "USDJPY".parse()?,
        direct: Some(DirectFx {
            market_collateral: collateral.clone(),
            quotes: FxQuotes {
                spot: "USDJPY.SPOT".parse()?,
                on: short("USDJPY.ON")?,
                tn: short("USDJPY.TN")?,
                sn: None,
                forwards: vec![],
            },
            implied,
        }),
        collaterals: vec![FxCollateral {
            id: collateral,
            base,
            quote: tona.index.parse()?,
        }],
        via: None,
    })
}

#[test]
fn translation_fx() -> Result<()> {
    // スポットだけなら割引曲線で、短期ポイントがあれば連鎖した最初の受渡のレートから X(0) を求める。
    let [tona, sofr] = cases()?;
    let inputs = Inputs::new(&[&tona, &sofr], tona.as_of, "")?;
    let ctx = inputs.context();
    let definition = definition(&[(&tona, &COLLATERAL), (&sofr, &COLLATERAL)])?;
    let ois_quotes: Vec<_> = fixture_quotes(&tona)?
        .into_iter()
        .chain(fixture_quotes(&sofr)?)
        .collect();
    let spot = "USDJPY.SPOT".parse::<QuoteId>()?;
    let on = "USDJPY.ON".parse::<QuoteId>()?;
    let tn = "USDJPY.TN".parse::<QuoteId>()?;
    let [jpy, usd] = [Currency::JPY, Currency::USD];
    let jpy_swap = convention_swap(&tona, None, "1Y", 0.005, PayReceive::Receive, 0.0)?;
    let usd_swap = convention_swap(&sofr, None, "1Y", 0.045, PayReceive::Receive, 0.0)?;
    let prepare = |market: &CalibratedMarket| {
        [&jpy_swap, &usd_swap]
            .map(|swap| prepare_trade(swap, tona.as_of, &COLLATERAL, market, &ctx))
    };
    let front_jpy = tona.quote(&tona.nodes[0])?;
    let usd_rate = sofr.quote(&sofr.nodes[2])?;

    for implied in [
        None,
        Some(FxImpliedCurve {
            id: "USD-under-JPY".parse()?,
            extrapolation: Extrapolation::default(),
            long_end: vec![],
        }),
    ] {
        let short = implied.is_some();
        let definition = definition
            .clone()
            .with_fx([usd_jpy(&tona, &sofr, implied)?])?;
        let quotes = QuoteSet::new(
            ois_quotes
                .iter()
                .cloned()
                .chain([(spot.clone(), 150.0)])
                .chain(
                    short
                        .then(|| [(on.clone(), -2.0), (tn.clone(), 3.0)])
                        .into_iter()
                        .flatten(),
                ),
        )?;
        let base = calibrated(&definition, &quotes, &ctx, 1e-14)?;
        let [first, second] = prepare(&base);
        let prepared = [first?, second?];
        let mut targets = vec![(&spot, 1.5), (&front_jpy, BP), (&usd_rate, BP)];
        if short {
            targets.extend([(&on, 1.0), (&tn, 1.0)]);
        }
        let local = quote_delta(&prepared, &base, SEQUENTIAL)?;
        assert!(local.failures().is_empty());
        for (reporting, foreign) in [(jpy, usd), (usd, jpy)] {
            let pair = CcyPair::new(foreign, reporting)?;
            let floor = absolute_tolerance(&ctx, reporting)?;
            // 取引は市場ごとに同じ条件で準備し、PV の換算は型付きの FX レートで別途行う。
            let reprice = |id: &QuoteId, bump| {
                let bumped = calibrated(
                    &definition,
                    &bumped(&quotes, id, bump)?,
                    &ctx,
                    BUMP_TOLERANCE,
                )?;
                converted_pv(prepare(&bumped), bumped.market(), pair)
            };
            let report =
                check_translation_bumps(&prepared, &base, reporting, floor, &targets, reprice)?;
            let vector = vector_for(report.by_ccy(), reporting)?;
            if !short {
                let x0 = base
                    .market()
                    .structure()
                    .fx_today_plan(pair)?
                    .rate(base.market())?
                    .rate();
                let foreign_trade = if foreign == jpy {
                    &prepared[0]
                } else {
                    &prepared[1]
                };
                let foreign_pv = evaluate(foreign_trade, base.market())?.get(foreign);
                let sign = if foreign == usd { 1.0 } else { -1.0 };
                assert_close!(
                    sensitivity(vector, &spot),
                    sign * 0.01 * x0 * foreign_pv,
                    1e-10,
                    floor
                );
                let frozen_fx = sensitivity(vector_for(local.by_ccy(), reporting)?, &front_jpy)
                    + x0 * sensitivity(vector_for(local.by_ccy(), foreign)?, &front_jpy);
                // スポット受渡日までの JPY の割引を経由する項が、現地デルタの換算に加わる。
                assert!((sensitivity(vector, &front_jpy) - frozen_fx).abs() > floor);
            }
        }
    }
    Ok(())
}

fn converted_pv(trades: [Result<Prepared>; 2], market: &Market<f64>, pair: CcyPair) -> Result<f64> {
    let mut pv = MoneyBag::new();
    for trade in trades {
        for term in evaluate(&trade?, market)?.terms() {
            pv.add(term);
        }
    }
    let rate = market.structure().fx_today_plan(pair)?.rate(market)?;
    Ok(pv.convert(pair.quote(), &[rate])?.amount())
}

fn check_translation_bumps(
    prepared: &[Prepared],
    base: &CalibratedMarket,
    reporting: Currency,
    floor: f64,
    targets: &[(&QuoteId, f64)],
    mut reprice: impl FnMut(&QuoteId, f64) -> Result<f64>,
) -> Result<DeltaReport> {
    let report = quote_delta_in(prepared, base, reporting, SEQUENTIAL)?;
    let executor = ThreadExecutor::new(NonZeroUsize::MIN.saturating_add(1));
    let concurrent = quote_delta_in(prepared, base, reporting, ExecutionContext::new(&executor))?;
    assert!(report.failures().is_empty());
    assert!(concurrent.failures().is_empty());
    assert_eq!(report.by_ccy().len(), 1);
    assert_eq!(concurrent.by_ccy().len(), 1);
    let vector = vector_for(report.by_ccy(), reporting)?;
    let parallel = vector_for(concurrent.by_ccy(), reporting)?;
    for (id, _) in base.quotes().iter() {
        assert_close!(
            sensitivity(vector, id),
            sensitivity(parallel, id),
            1e-10,
            floor
        );
    }
    // 入力 quote の自然な単位に戻して、ソルバーの残差が PV に伝わる幅を求める。
    let scale = vector
        .entries()
        .iter()
        .map(|(id, value)| {
            let unit = targets
                .iter()
                .find_map(|&(target, unit)| (id == target).then_some(unit))
                .unwrap_or(BP);
            value.abs() / unit
        })
        .sum::<f64>();
    for &(id, unit) in targets {
        let aad = sensitivity(vector, id);
        assert!(aad.abs() > floor, "{reporting} {id}: {aad:e}");
        check_aad_vs_bump(
            |h| reprice(id, h),
            aad,
            unit,
            BUMP_TOLERANCE * scale,
            floor,
            id.as_str(),
        )?;
    }
    Ok(report)
}

#[test]
fn quote_delta_reporting_requires_fx_for_foreign_pv() -> Result<()> {
    // FX の宣言が無い場合と、宣言があってスポットが欠けた場合を照合する。
    let [tona, sofr] = cases()?;
    let inputs = Inputs::new(&[&tona, &sofr], tona.as_of, "")?;
    let ctx = inputs.context();
    let definition = definition(&[(&tona, &COLLATERAL), (&sofr, &COLLATERAL)])?;
    let quotes = QuoteSet::new(
        fixture_quotes(&tona)?
            .into_iter()
            .chain(fixture_quotes(&sofr)?),
    )?;
    let with_fx = definition.clone().with_fx([usd_jpy(&tona, &sofr, None)?])?;
    for def in [&definition, &with_fx] {
        let market = calibrate(def, &quotes, &ctx)?;
        let trade = prepare_ois(
            &sofr,
            &sofr.nodes[2],
            0.045,
            &COLLATERAL,
            market.market(),
            &ctx,
        )?;
        let trades = [trade];
        assert!(matches!(
            quote_delta_in(&trades, &market, Currency::JPY, SEQUENTIAL),
            Err(Error::MissingData(_))
        ));
        let native = quote_delta(&trades, &market, SEQUENTIAL)?;
        let same = quote_delta_in(&trades, &market, Currency::USD, SEQUENTIAL)?;
        assert!(same.failures().is_empty());
        let lhs = vector_for(native.by_ccy(), Currency::USD)?;
        let rhs = vector_for(same.by_ccy(), Currency::USD)?;
        for (id, _) in quotes.iter() {
            assert_close!(
                sensitivity(lhs, id),
                sensitivity(rhs, id),
                1e-10,
                absolute_tolerance(&ctx, Currency::USD)?
            );
        }
        let empty = quote_delta_in(&[], &market, Currency::JPY, SEQUENTIAL)?;
        assert!(empty.by_ccy().is_empty() && empty.failures().is_empty());
    }
    Ok(())
}

#[test]
fn quote_delta_keeps_successful_trades_when_calibration_fails() -> Result<()> {
    // TONA の欠損は下流の予想曲線にも伝わるが、独立した SOFR の取引は集計される。
    let [tona, sofr] = cases()?;
    let mut forecast = tona.clone();
    forecast.index = "JPY-FORECAST".into();
    let inputs = Inputs::new(&[&tona, &forecast, &sofr], tona.as_of, "")?;
    let ctx = inputs.context();
    let def = definition(&[
        (&tona, &COLLATERAL),
        (&forecast, &COLLATERAL),
        (&sofr, &COLLATERAL),
    ])?;
    let quotes = QuoteSet::new(
        fixture_quotes(&forecast)?
            .into_iter()
            .chain(fixture_quotes(&sofr)?),
    )?;
    let cm = calibrate(&def, &quotes, &ctx)?;
    assert!(matches!(
        cm.groups()[1].outcome(),
        GroupOutcome::Unavailable
    ));
    let trades = [&tona, &sofr, &forecast]
        .into_iter()
        .map(|case| prepare_ois(case, &case.nodes[1], 0.01, &COLLATERAL, cm.market(), &ctx))
        .collect::<Result<Vec<_>>>()?;
    let report = quote_delta(&trades, &cm, SEQUENTIAL)?;
    assert!(matches!(
        report.failures(),
        [(0, Error::MissingData(_)), (2, Error::MissingData(_))]
    ));
    assert!(matches!(report.by_ccy(), [vector] if vector.ccy() == Currency::USD));
    assert_ne!(report.by_ccy()[0].entries(), []);
    Ok(())
}

#[test]
fn settlement_at_valuation_instant_has_pv_without_delta() -> Result<()> {
    // 評価時点に限って残る支払いは、結合に曲線が必要でも θ を読み出さない。
    let [tona, _] = cases()?;
    let inputs = Inputs::new(&[&tona], tona.as_of, "")?;
    let ctx = inputs.context();
    let def = definition(&[(&tona, &COLLATERAL)])?;
    let cm = calibrate(&def, &QuoteSet::new(fixture_quotes(&tona)?)?, &ctx)?;
    let product = deposit(&tona, "2024-12-16".parse()?, tona.as_of, Position::Long)?;
    let settled = prepare_trade(&product, tona.as_of, &COLLATERAL, &cm, &ctx)?;
    let report = quote_delta(&[settled], &cm, SEQUENTIAL)?;
    assert!(report.by_ccy().is_empty() && report.failures().is_empty());

    let config = toml::from_str(
        r#"
version = 1
reporting_currency = "JPY"
settlement = "include"
[[discounting]]
currency = "JPY"
collateral = "uncollateralized"
curve = "JPY-TONA"
"#,
    )
    .map_err(|error| Error::InvalidData(error.to_string()))?;
    let ctx = Context {
        config: &config,
        ..ctx
    };
    let due_today = prepare_trade(&product, tona.as_of, &COLLATERAL, &cm, &ctx)?;
    assert_eq!(due_today.blocks(), []);
    let gradient = param_gradient(&due_today, &cm)?;
    assert!(gradient.pv().get(Currency::JPY) > 0.0);
    assert!(matches!(gradient.by_ccy(), [vector] if vector.entries().is_empty()));
    let report = quote_delta(&[due_today], &cm, SEQUENTIAL)?;
    assert!(report.failures().is_empty());
    assert!(matches!(report.by_ccy(), [vector] if vector.entries().is_empty()));
    Ok(())
}
`````

`````
crates/rhodium-risk/tests/ui/quote_delta_with_param_gradient.rs
use rhodium_market::{quote::QuoteId, structure::ParamId};
use rhodium_risk::vector::RiskVector;

fn add(delta: &RiskVector<QuoteId>, gradient: &RiskVector<ParamId>) {
    // 同じ通貨であっても、quote と θ の感応度には異なるキーを使う。
    let _ = delta.checked_add(gradient);
}

fn main() {}
`````

`````
crates/rhodium-risk/tests/ui/quote_delta_with_param_gradient_ok.rs
use rhodium_core::Result;
use rhodium_market::quote::QuoteId;
use rhodium_risk::vector::RiskVector;

// quote 空間の感応度どうしは足せる。
fn add(left: &RiskVector<QuoteId>, right: &RiskVector<QuoteId>) -> Result<RiskVector<QuoteId>> {
    left.checked_add(right)
}

fn main() {
    let _ = add;
}
`````

`````
crates/rhodium-risk/tests/ui/tape_escape.rs
use std::cell::Cell;

use scioxide_ad::{Tape, TapeOptions, Var};
use scioxide_core::{ExecutionContext, SciOxideError};

fn main() -> Result<(), SciOxideError> {
    let mut tape = Tape::<f64>::new(TapeOptions::new());
    let escaped = Cell::new(None::<Var<'_, f64>>);
    let _ = tape.linearize(
        |inputs, outputs| {
            // 記録中の変数を外へ保存すると、tape を使い回した後にも読めてしまう。
            escaped.set(inputs.first().copied());
            outputs.extend_from_slice(inputs);
            Ok(())
        },
        &[1.0],
        ExecutionContext::sequential(),
    )?;
    Ok(())
}
`````

`````
crates/rhodium-risk/tests/ui/tape_escape_ok.rs
use scioxide_ad::{Tape, TapeOptions};
use scioxide_core::{ExecutionContext, SciOxideError};

fn main() -> Result<(), SciOxideError> {
    let mut tape = Tape::<f64>::new(TapeOptions::new());
    let linearization = tape.linearize(
        |inputs, outputs| {
            // 変数は記録の中だけで使い、値と微分は記録の結果から f64 で受け取る。
            outputs.extend(inputs.iter().map(|&x| x * x));
            Ok(())
        },
        &[3.0],
        ExecutionContext::sequential(),
    )?;
    let _value: f64 = linearization.values()[0];
    Ok(())
}
`````

`````
crates/rhodium-risk/tests/xccy_quote_delta.rs
//! R6 の網羅表の行: 外貨担保曲線（FX ポイントと `MtM` XCCY）を通る `MtM` XCCY、FX フォワード、
//! FX スワップの quote デルタ（Δq と Δx）を、quote ごとに全体を再校正した中心差分と照合する。

use rhodium_core::{Error, Result, currency::Currency, date::Date, fx::FxRate, money::Money};
use rhodium_market::quote::{QuoteId, QuoteKind, QuoteSet};
use rhodium_pricing::{
    calibration::{CalibratedMarket, GroupOutcome, calibrate_with_residual_tolerance},
    context::Context,
    definition::MarketDefinition,
    prepare::{PrepareInput, Prepared, evaluate, prepare},
};
use rhodium_products::standard::{FxForward, FxSwap, Position, StandardProduct};
use rhodium_risk::{delta::quote_delta_in, vector::RiskVector};
use rhodium_testkit::check_aad_vs_bump;
use scioxide_core::ExecutionContext;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、XCCY のデルタに使わない関数がある"
)]
#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
pub mod ois;

#[allow(
    dead_code,
    reason = "共有の fixture のうち、XCCY のデルタに使わない関数がある"
)]
#[path = "../../rhodium-pricing/tests/fixtures/xccy.rs"]
pub mod xccy;

use xccy::{SPOT, XccyMarket, jpy_cash};

const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();
const BP: f64 = 1e-4;
/// バンプの再校正の残差の許容誤差。品質の照合手順は通常の 1e-14 の 1/100 を求めるが、XCCY のパーの
/// 残差は倍精度の分解能で約 1.5e-16 までしか下がらない。届く最小の値で解き、`δ_V` に伝える。
const BUMP_TOLERANCE: f64 = 1e-15;

fn calibrated(
    definition: &MarketDefinition,
    quotes: &QuoteSet,
    ctx: &Context<'_>,
) -> Result<CalibratedMarket> {
    let market = calibrate_with_residual_tolerance(definition, quotes, ctx, BUMP_TOLERANCE)?;
    for group in market.groups() {
        if !matches!(group.outcome(), GroupOutcome::Solved(_)) {
            return Err(Error::Calibration(format!("{group:?}")));
        }
    }
    Ok(market)
}

fn bumped(quotes: &QuoteSet, id: &QuoteId, bump: f64) -> Result<QuoteSet> {
    QuoteSet::new(
        quotes
            .iter()
            .map(|(key, value)| (key.clone(), if key == id { value + bump } else { value })),
    )
}

fn sensitivity(vector: &RiskVector<QuoteId>, id: &QuoteId) -> f64 {
    vector
        .entries()
        .binary_search_by(|(key, _)| key.cmp(id))
        .ok()
        .map_or(0.0, |index| vector.entries()[index].1)
}

/// 照合する取引と、その PV の絶対の下限（JPY）。
fn products() -> Result<Vec<(&'static str, StandardProduct, f64)>> {
    let pair = XccyMarket::pair()?;
    let notional = Money::new(Currency::USD, 1_000_000.0);
    // quote から 10bp 外した 18M の `MtM` XCCY。テンプレートの元本は JPY 1 円。
    let xccy = XccyMarket::template("18M")?.build_for_trade(
        -0.0030,
        trade_date()?,
        &XccyMarket::new()?.inputs("")?.context().refdata.clone(),
    )?;
    let forward = StandardProduct::FxForward(FxForward::new(
        notional,
        FxRate::forward(pair, 147.0, XccyMarket::settlement("2025-10-09")?)?,
        Position::Long,
    )?);
    let swap = StandardProduct::FxSwap(FxSwap::new(
        notional,
        FxRate::forward(pair, SPOT, XccyMarket::settlement("2025-01-09")?)?,
        FxRate::forward(pair, 145.0, XccyMarket::settlement("2026-01-09")?)?,
        Position::Long,
    )?);
    Ok(vec![
        ("MtM XCCY", xccy, 1e-14),
        ("FX forward", forward, 1e-6),
        ("FX swap", swap, 1e-6),
    ])
}

fn trade_date() -> Result<Date> {
    "2025-01-07".parse()
}

fn prepare_trade(
    product: &StandardProduct,
    market: &CalibratedMarket,
    ctx: &Context<'_>,
) -> Result<Prepared> {
    prepare(
        PrepareInput::Standard {
            product,
            trade_date: trade_date()?,
        },
        &jpy_cash()?,
        &[],
        market.market().structure(),
        ctx,
    )
}

/// 評価時点の X(0) で USD の受渡を JPY に換算した PV。
fn jpy_pv(prepared: &Prepared, market: &CalibratedMarket) -> Result<f64> {
    let pv = evaluate(prepared, market.market())?;
    let x0 = market
        .market()
        .structure()
        .fx_today_plan(XccyMarket::pair()?)?
        .rate(market.market())?
        .rate();
    Ok(pv.get(Currency::JPY) + pv.get(Currency::USD) * x0)
}

/// quote の種類の感応度の単位（金利とスプレッドは 1bp、FX スポットは相対 1%、FX ポイントは 1 pip）。
fn unit(market: &CalibratedMarket, id: &QuoteId) -> Result<f64> {
    Ok(
        match market
            .quote_kind(id)
            .ok_or_else(|| Error::Internal(format!("unknown kind of {id}")))?
        {
            QuoteKind::FxSpot => 0.01 * SPOT,
            QuoteKind::FxForwardPoints => 1.0,
            _ => BP,
        },
    )
}

#[test]
fn xccy_fx_forward_and_fx_swap_quote_deltas_match_bumps() -> Result<()> {
    // JPY 現金担保の下で、USD の受渡は外貨担保曲線、JPY の受渡は TONA で割り引き、報告通貨 JPY に
    // 評価時点の X(0) で換算する。Δq は TONA、SOFR の OIS と XCCY のベーシスへ、Δx は USDJPY の
    // スポット（相対 1%）と ON、TN、日付付きのポイント（1 pip）へ、校正の随伴で戻る。
    let fixture = XccyMarket::new()?;
    let inputs = fixture.inputs("")?;
    let ctx = inputs.context();
    let definition = fixture.definition()?;
    let quotes = fixture.quotes()?;
    let base = calibrated(&definition, &quotes, &ctx)?;
    for (name, product, floor) in products()? {
        let prepared = prepare_trade(&product, &base, &ctx)?;
        let report = quote_delta_in(
            std::slice::from_ref(&prepared),
            &base,
            Currency::JPY,
            SEQUENTIAL,
        )?;
        assert!(report.failures().is_empty());
        let [vector] = report.by_ccy() else {
            return Err(Error::Internal("expected one reporting currency".into()));
        };
        // 再校正の残差の許容誤差を、quote ごとの感応度（自然な単位）の和で PV の誤差に伝える。
        let scale = vector
            .entries()
            .iter()
            .map(|(id, delta)| Ok(delta.abs() / unit(&base, id)?))
            .sum::<Result<f64>>()?;
        for (id, _) in quotes.iter() {
            let unit = unit(&base, id)?;
            check_aad_vs_bump(
                |h| {
                    let market = calibrated(&definition, &bumped(&quotes, id, h)?, &ctx)?;
                    jpy_pv(&prepare_trade(&product, &market, &ctx)?, &market)
                },
                sensitivity(vector, id),
                unit,
                BUMP_TOLERANCE * scale,
                floor,
                &format!("{name} {id}"),
            )?;
        }
        // 照合が空振りでないことを確かめる。`MtM` XCCY はベーシスに、FX の取引はスポットに感応する
        // （`MtM` XCCY は元本のリセットで FX の露出がほぼ打ち消されるので、スポットは確かめない）。
        let expected = if name == "MtM XCCY" {
            "USDJPY.XCCY.18M"
        } else {
            "USDJPY.SPOT"
        };
        assert!(
            sensitivity(vector, &expected.parse()?).abs() > floor,
            "{name} {expected}"
        );
    }
    Ok(())
}
`````

`````
crates/rhodium-testkit/Cargo.toml
[package]
name = "rhodium-testkit"
description = "テストの補助と独立なオラクルを持つ開発専用の crate"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
publish.workspace = true

[lints]
workspace = true

[dependencies]
scioxide-core.workspace = true
trybuild.workspace = true
`````

`````
crates/rhodium-testkit/src/lib.rs
//! テストの補助と独立なオラクルを持つ開発専用の crate。

/// 呼び出し元の `tests/ui/` にある `<case>_ok.rs` を実行し、`<case>.rs` の診断を照合する。
///
/// ケース名は ASCII 小文字・数字・アンダースコアからなる空でない文字列リテラルで指定する。
#[macro_export]
macro_rules! compile_fail_pair {
    ($case:literal) => {{
        const CASE: &str = {
            // 登録検査と同じく、通常の引用符で書いたケース名を受け付ける。
            let literal = stringify!($case).as_bytes();
            assert!(
                literal.len() > 2 && literal[0] == b'"' && literal[literal.len() - 1] == b'"',
                "expected a nonempty quoted UI case name",
            );
            let mut index = 1;
            while index < literal.len() - 1 {
                let byte = literal[index];
                assert!(
                    byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'_',
                    "expected a snake_case UI case name",
                );
                index += 1;
            }
            $case
        };
        $crate::run_compile_fail_pair(CASE);
    }};
}

#[doc(hidden)]
pub fn run_compile_fail_pair(case: &str) {
    // trybuild はテストを実行する crate のマニフェストを起点にパスを解決する。
    let cases = trybuild::TestCases::new();
    cases.pass(format!("tests/ui/{case}_ok.rs"));
    cases.compile_fail(format!("tests/ui/{case}.rs"));
}

/// 数値を相対誤差と絶対誤差の大きい方で照合する。
///
/// `|x − y| ≤ max(rel · M, abs)` を確かめる。M は `max(|x|, |y|)` で、5 番目の引数を渡すと
/// それも M に含める。打ち消しの大きい和（ネットした PV やデルタ）は、和の項の絶対値の和を
/// 5 番目に渡す（品質の許容誤差の M = max(|x|, |y|, Σ|項|)）。
#[macro_export]
macro_rules! assert_close {
    ($actual:expr, $expected:expr, $relative:expr, $absolute:expr $(,)?) => {
        $crate::assert_close!($actual, $expected, $relative, $absolute, 0.0)
    };
    ($actual:expr, $expected:expr, $relative:expr, $absolute:expr, $scale:expr $(,)?) => {{
        let actual: f64 = $actual;
        let expected: f64 = $expected;
        let relative: f64 = $relative;
        let absolute: f64 = $absolute;
        let scale: f64 = $scale;
        let magnitude = actual.abs().max(expected.abs()).max(scale);
        let tolerance = (relative * magnitude).max(absolute);
        assert!(
            actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= tolerance,
            "actual {actual}, expected {expected}, tolerance {tolerance}"
        );
    }};
}

/// 全体を再校正する中心差分と、報告単位あたりの一次感応度を照合する。
///
/// `reprice` は自然な単位で指定した quote の増分に対する PV を返す。
/// `unit` は報告単位の幅、`pv_error` は一つの PV への校正誤差の伝播幅。
/// `absolute` は報告単位あたりの感応度に対する絶対誤差の下限。
///
/// # Errors
///
/// 校正または再評価で `reprice` が返したエラーを伝える。
///
/// # Panics
///
/// Richardson の打ち切り誤差と数値誤差を加えた許容幅で一致しない場合。
pub fn check_aad_vs_bump<E>(
    mut reprice: impl FnMut(f64) -> Result<f64, E>,
    aad: f64,
    unit: f64,
    pv_error: f64,
    absolute: f64,
    label: &str,
) -> Result<(), E> {
    let half = unit / 2.0;
    let mut central = |h| Ok((reprice(h)? - reprice(-h)?) * (unit / (2.0 * h)));
    let coarse = central(unit)?;
    let fine = central(half)?;
    let truncation = (coarse - fine).abs() / 3.0;
    let solver = pv_error * (unit / half);
    let rounding = (1e-10 * aad.abs().max(fine.abs())).max(absolute);
    let tolerance = 2.0 * truncation + solver + rounding;
    assert!(
        aad.is_finite() && fine.is_finite() && (aad - fine).abs() <= tolerance,
        "{label}: AAD {aad:e}, central difference {fine:e}, tolerance {tolerance:e}"
    );
    Ok(())
}

/// 読むたびに 1 ns 進む決定的な時計。
///
/// 期限を何回目の読み出しに置くかで、実行が止まる位置をスレッドの速さや実際の時刻によらず
/// 決められる。区間の開始後に期限が来る場合の取消の伝わり方を確かめるテストに使う。
/// `ClockBinding::new(&clock)` で実行文脈に結び、[`TickClock::after`] の瞬間を期限にする。
#[derive(Debug, Default)]
pub struct TickClock(std::sync::atomic::AtomicU64);

impl TickClock {
    /// まだ一度も読まれていない時計。
    #[must_use]
    pub const fn new() -> Self {
        Self(std::sync::atomic::AtomicU64::new(0))
    }

    /// これまでの読み出しの回数。
    #[must_use]
    pub fn reads(&self) -> u64 {
        self.0.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// `reads` 回の読み出しの後、次の読み出しで達する瞬間。
    ///
    /// 期限 `at` は読み出しの値が `at` 以上で達するので、0 回目から数えて `reads` 回目の
    /// 読み出し（値 `reads`）で初めて止まる。
    pub const fn after(reads: u64) -> scioxide_core::MonotonicInstant {
        scioxide_core::MonotonicInstant::from_duration_since_epoch(std::time::Duration::from_nanos(
            reads,
        ))
    }
}

impl scioxide_core::Clock for TickClock {
    fn now(&self) -> Result<scioxide_core::MonotonicInstant, scioxide_core::CapabilityError> {
        let tick = self.0.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        Ok(Self::after(tick))
    }
}
`````

`````
crates/rhodium-testkit/tests/assert_close.rs
//! 数値の照合が、品質の許容誤差の式 |x − y| ≤ max(ε·M, a) に従うことを確かめる。

use rhodium_testkit::assert_close;

#[test]
fn scale_widens_the_relative_tolerance_for_cancelling_sums() {
    // 1 + 1e-11 と 1 の差 1e-11 は、M = max(|x|, |y|) = 1 なら相対 1e-12 を超える。
    // 打ち消しの大きい和の項の絶対値の和 100 を尺度に入れれば、許容誤差は 1e-10 になる。
    assert_close!(1.0 + 1e-11, 1.0, 1e-12, 0.0, 100.0);
    // 尺度を渡さない形は M = max(|x|, |y|) を使う。
    assert_close!(1.0 + 1e-13, 1.0, 1e-12, 0.0);
    // ゼロに近い値は絶対の下限で比べる。
    assert_close!(1e-16, 0.0, 1e-12, 1e-14);
}

#[test]
#[should_panic(expected = "tolerance")]
fn difference_beyond_the_tolerance_fails() {
    // 尺度を渡さなければ、同じ差 1e-11 は相対 1e-12 の許容誤差を超える。
    assert_close!(1.0 + 1e-11, 1.0, 1e-12, 0.0);
}

#[test]
#[should_panic(expected = "tolerance")]
fn non_finite_values_never_match() {
    // NaN は許容誤差によらず一致しない。
    assert_close!(f64::NAN, f64::NAN, 1.0, 1.0);
}
`````

`````
crates/rhodium-testkit/tests/compile_fail_pairs.rs
//! ワークスペースの型テストが対になっていることを確認する。

use std::{
    collections::BTreeSet,
    fs, io,
    path::{Path, PathBuf},
};

// 各ケースが成功例と stderr を登録していること、登録されていないケースがないことを確認する。
// また、compile-fail の例を対の登録の外に置く書き方（doctest の compile_fail と trybuild の直接呼び出し）を拒否する。
// 各 crate のケースの実行と診断の照合は compile_fail_pair! に任せる。
#[test]
#[expect(
    clippy::disallowed_macros,
    reason = "開発用の統合テストで、Cargo が指定するソースツリーを検査する"
)]
fn each_compile_failure_has_an_executed_positive_case() -> io::Result<()> {
    let crates = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    for directory in children(&crates)? {
        if !directory.join("Cargo.toml").is_file() {
            continue;
        }
        check_ui(&directory)?;
        for root in ["src", "tests", "examples", "benches", "build.rs"] {
            check_sources(&directory.join(root))?;
        }
    }
    Ok(())
}

fn check_ui(directory: &Path) -> io::Result<()> {
    let tests = directory.join("tests");
    let ui = tests.join("ui");
    let harness = tests.join("compile_fail.rs");
    let mut expected = BTreeSet::new();
    if harness.is_file() {
        let source = fs::read_to_string(&harness)?;
        let mut cases = BTreeSet::new();
        for name in source.lines().filter_map(|line| {
            line.trim()
                .strip_prefix("compile_fail_pair!(\"")
                .and_then(|tail| tail.strip_suffix("\");"))
        }) {
            assert!(
                cases.insert(name),
                "{}: duplicate case {name}",
                harness.display()
            );
            expected.extend([
                ui.join(format!("{name}.rs")),
                ui.join(format!("{name}_ok.rs")),
                ui.join(format!("{name}.stderr")),
            ]);
        }
        assert!(
            !cases.is_empty(),
            "{}: missing registrations",
            harness.display()
        );
    }

    let actual = if ui.is_dir() {
        children(&ui)?
            .into_iter()
            .filter(|path| {
                path.is_file()
                    && matches!(
                        path.extension().and_then(|ext| ext.to_str()),
                        Some("rs" | "stderr")
                    )
            })
            .collect()
    } else {
        BTreeSet::new()
    };
    let missing = expected.difference(&actual).collect::<Vec<_>>();
    let unregistered = actual.difference(&expected).collect::<Vec<_>>();
    assert!(
        missing.is_empty() && unregistered.is_empty(),
        "{}: missing {missing:?}, unregistered {unregistered:?}",
        directory.display(),
    );
    Ok(())
}

// 対の登録を通らない compile-fail の例を拒否する。doctest の compile_fail は成功例と診断を
// 対にできず、trybuild の直接呼び出しは tests/ui の外のファイルを検査から外せるため。
fn check_sources(path: &Path) -> io::Result<()> {
    if path.is_dir() {
        for child in children(path)? {
            check_sources(&child)?;
        }
    } else if path.is_file()
        && path.extension().is_some_and(|ext| ext == "rs")
        // この検査自身は探す文字列をソースに含むので対象から外す。
        && !path.ends_with("rhodium-testkit/tests/compile_fail_pairs.rs")
    {
        // trybuild の直接呼び出しはマクロの実装ファイルだけに許可する。
        let allows_trybuild = path.ends_with("rhodium-testkit/src/lib.rs");
        let source = fs::read_to_string(path)?;
        for (number, line) in (1_u64..).zip(source.lines()) {
            let line = line.trim();
            assert!(
                allows_trybuild
                    || !(line.contains("trybuild::")
                        || line.contains("TestCases")
                        || line.contains("run_compile_fail_pair")),
                "{}:{number}: register the case with compile_fail_pair! instead of calling trybuild",
                path.display(),
            );
            let body = line
                .strip_prefix("///")
                .or_else(|| line.strip_prefix("//!"))
                .or_else(|| line.strip_prefix('*'))
                .map(str::trim_start);
            let fence = body.and_then(|body| {
                body.strip_prefix("```")
                    .or_else(|| body.strip_prefix("~~~"))
            });
            assert!(
                !fence.is_some_and(|header| header
                    .split([',', ' ', '\t'])
                    .any(|flag| flag == "compile_fail")),
                "{}:{number}: move the compile-fail example to a paired trybuild case",
                path.display(),
            );
        }
    }
    Ok(())
}

fn children(directory: &Path) -> io::Result<BTreeSet<PathBuf>> {
    fs::read_dir(directory)?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect()
}
`````

`````
crates/rhodium-xva/Cargo.toml
[package]
name = "rhodium-xva"
description = "エクスポージャー、担保、XVA の指標と感応度"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
publish.workspace = true

[lints]
workspace = true

[dependencies]
rhodium-risk.workspace = true
`````

`````
crates/rhodium-xva/src/lib.rs
//! エクスポージャー、担保、XVA の指標と感応度。

#![forbid(unsafe_code)]

// XVA を実装するまで、依存先を明示的に参照する。
use rhodium_risk as _;
`````

`````
crates/rhodium/Cargo.toml
[package]
name = "rhodium"
description = "rhodium の公開ライブラリ API"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true
publish.workspace = true

[lints]
workspace = true

[dependencies]
csv.workspace = true
rhodium-core.workspace = true
rhodium-products.workspace = true
rhodium-market.workspace = true
rhodium-pricing.workspace = true
rhodium-risk.workspace = true
rhodium-xva.workspace = true
serde.workspace = true
serde_json = "1"
toml = "1.2"

[dev-dependencies]
criterion.workspace = true
rhodium-testkit.workspace = true
scioxide-core.workspace = true
tempfile = "3.27"

[[bench]]
name = "m0"
harness = false

[[bench]]
name = "active_market"
harness = false

[[bench]]
name = "m1"
harness = false
`````

`````
crates/rhodium/benches/active_market.rs
//! 固定した 1,000 件の OIS の quote デルタを、取引が読まない曲線の本数だけ変えて測る。

#[allow(dead_code, reason = "校正用 fixture の定義関数は未使用曲線を足せない")]
#[path = "../../rhodium-pricing/tests/fixtures/ois.rs"]
pub mod ois;

use std::{hint::black_box, time::Duration};

use criterion::{BenchmarkId, Criterion};
use rhodium_core::{Result, parallel::ExecutionContext, time::Instant};
use rhodium_market::{
    curve::{CurveRole, CurveSpec, Extrapolation},
    quote::{QuoteId, QuoteKind, QuoteSet},
    structure::{BlockId, ParamId},
};
use rhodium_pricing::{
    calibration::{CalibratedMarket, GroupOutcome, calibrate},
    definition::{
        CalibrationTarget, ElementSelector, InstrumentBinding, MarketDefinition, QuoteTransform,
    },
    prepare::Prepared,
    resolve::CollateralKey,
};
use rhodium_risk::delta::quote_delta;

const TRADES: usize = 1_000;
const COLLATERAL: CollateralKey = CollateralKey::Uncollateralized;
const SEQUENTIAL: ExecutionContext<'static> = ExecutionContext::sequential();

fn portfolio(extra_curves: usize) -> Result<(CalibratedMarket, Vec<Prepared>)> {
    let [case, _] = ois::cases()?;
    let inputs = ois::Inputs::new(&[&case], case.as_of, "")?;
    let ctx = inputs.context();
    let mut curves = vec![case.curve()?];
    let mut bindings = Vec::new();
    let mut quotes = ois::fixture_quotes(&case)?;
    let active = case.block()?;
    for (element, node) in (0..).zip(&case.nodes) {
        bindings.push(case.binding(
            node,
            &COLLATERAL,
            ElementSelector::Param(ParamId::new(active.clone(), element)),
        )?);
    }
    let pillars = (2027..2042)
        .map(|year| format!("{year}-12-31T00:00:00Z").parse::<Instant>())
        .collect::<Result<Vec<_>>>()?;
    for curve in 0..extra_curves {
        let name = format!("PASSIVE-{curve:04}");
        let id = name.parse()?;
        curves.push(CurveSpec::new(
            id,
            case.currency,
            [CurveRole::Discount],
            pillars.clone(),
            Extrapolation::FlatAverageForward,
        )?);
        let block = BlockId::nodes(name.parse()?);
        // 直接マークで独立に校正し、解いた曲線や依存グループを取引側に加えない。
        for element in 0..15 {
            let quote: QuoteId = format!("{name}.N{element}").parse()?;
            bindings.push(InstrumentBinding::new(
                quote.clone(),
                QuoteKind::Rate,
                CalibrationTarget::DirectMark,
                QuoteTransform::Identity,
                COLLATERAL.clone(),
                ElementSelector::Param(ParamId::new(block.clone(), element)),
            )?);
            quotes.push((quote, -0.01 * f64::from(element + 1)));
        }
    }
    let definition = MarketDefinition::new(
        case.currency.to_string().parse()?,
        "Etc/UTC".parse()?,
        curves,
        bindings,
    )?;
    let market = calibrate(&definition, &QuoteSet::new(quotes)?, &ctx)?;
    assert!(
        market
            .groups()
            .iter()
            .all(|group| matches!(group.outcome(), GroupOutcome::Solved(_)))
    );
    let trades = case
        .nodes
        .iter()
        .cycle()
        .take(TRADES)
        .map(|node| {
            ois::prepare_ois(
                &case,
                node,
                node.rate + 0.001,
                &COLLATERAL,
                market.market(),
                &ctx,
            )
        })
        .collect::<Result<Vec<_>>>()?;
    assert!(
        trades
            .iter()
            .all(|trade| trade.blocks() == std::slice::from_ref(&active))
    );
    Ok((market, trades))
}

#[expect(clippy::expect_used, reason = "計測が失敗した場合はそこで停止する")]
fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    let mut group = criterion.benchmark_group("active_market/quote_delta_1000");
    group
        .sample_size(10)
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3));
    for extra_curves in [0, 32, 128, 512] {
        let (market, trades) = portfolio(extra_curves).expect("synthetic portfolio");
        let report = quote_delta(&trades, &market, SEQUENTIAL).expect("reference delta");
        assert!(report.failures().is_empty(), "{:?}", report.failures());
        assert!(
            report
                .by_ccy()
                .iter()
                .any(|column| !column.entries().is_empty()),
            "expected nonzero quote delta"
        );
        group.bench_function(BenchmarkId::from_parameter(extra_curves), |b| {
            b.iter(|| black_box(quote_delta(&trades, &market, SEQUENTIAL).expect("quote delta")));
        });
    }
    group.finish();
    criterion.final_summary();
}
`````

`````
crates/rhodium/benches/m0.rs
//! M0-9 の所要時間: 校正と、1 万件のスワップのデスクの読み込み、準備、PV、quote デルタ。
//!
//! `cargo bench -p rhodium --bench m0` で実行する。入力は `tests/support/desk.rs` が一時
//! ディレクトリに書く合成のデスク（JPY と USD の 1M〜30Y の OIS 曲線と慣行どおりの OIS）。
//! 結果は日付、コミット、マシンとともに `docs/benchmarks/` に残す。

#[allow(
    dead_code,
    reason = "テストと共有する生成器のうち、計測で使わない定数がある"
)]
#[path = "../tests/support/desk.rs"]
mod desk;

use std::{hint::black_box, num::NonZeroUsize, path::Path};

use criterion::{BatchSize, Criterion};
use rhodium::{
    prelude::{ExecutionContext, ThreadExecutor},
    session::{GroupOutcome, Session},
    workspace::Workspace,
};
use rhodium_core::time::Instant;
use rhodium_market::quote::QuoteOverride;
use tempfile::TempDir;

/// 計測するデスクの取引の数。
const TRADES: usize = 10_000;

// 逐次（1）、一部のコア（4）、使えるすべてのスレッドの三つで測る。
fn thread_counts() -> Vec<NonZeroUsize> {
    let available = std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN);
    let mut counts = vec![
        NonZeroUsize::MIN,
        NonZeroUsize::MIN.saturating_add(3),
        available,
    ];
    counts.sort_unstable();
    counts.dedup();
    counts
}

#[expect(
    clippy::expect_used,
    reason = "計測の入力を作れなければ、計測を続けられない"
)]
fn desk(trades: usize) -> TempDir {
    let dir = tempfile::tempdir().expect("temporary directory for the desk");
    desk::write(dir.path(), trades).expect("desk workspace");
    dir
}

#[expect(
    clippy::expect_used,
    reason = "計測の入力を読めなければ、計測を続けられない"
)]
fn load(dir: &Path) -> (Workspace, Instant) {
    let workspace = Workspace::load(dir).expect("desk workspace loads");
    let as_of = workspace
        .resolve_date(desk::VALUATION_DATE.parse().expect("valuation date"))
        .expect("valuation close");
    (workspace, as_of)
}

/// 計測の前に、校正の全グループが解け、全取引の準備と PV が成功し、quote デルタが取引を
/// 除いていないことを確かめる。
///
/// 失敗した計算は速く終わるので、確かめずに測ると失敗を高速化として記録してしまう。
#[expect(
    clippy::expect_used,
    reason = "計測の前提が崩れていれば、計測を続けられない"
)]
fn validate(session: &Session, execution: ExecutionContext<'_>, trades: usize) {
    let unsolved: Vec<_> = session
        .calibrated_market()
        .groups()
        .iter()
        .filter(|group| !matches!(group.outcome(), GroupOutcome::Solved(_)))
        .collect();
    assert!(
        unsolved.is_empty(),
        "unsolved calibration groups: {unsolved:?}"
    );
    let pv = session.pv(execution).expect("pv");
    let failed: Vec<_> = pv.iter().filter(|(_, result)| result.is_err()).collect();
    assert!(failed.is_empty(), "failed trades: {failed:?}");
    assert_eq!(pv.len(), trades);
    if trades > 0 {
        let delta = session.quote_delta(execution).expect("quote delta");
        assert!(
            delta.failures().is_empty(),
            "excluded trades: {:?}",
            delta.failures()
        );
    }
}

#[expect(
    clippy::expect_used,
    reason = "計測の入力から市場を作れなければ、計測を続けられない"
)]
fn calibration(c: &mut Criterion) {
    // 取引の無いデスクで Session::new を測ると、校正（JPY と USD の 15 ピラーの曲線）だけになる。
    let dir = desk(0);
    let (workspace, as_of) = load(dir.path());
    let session =
        Session::new(workspace, as_of, ExecutionContext::sequential()).expect("calibration");
    validate(&session, ExecutionContext::sequential(), 0);
    drop(session);
    c.bench_function("m0/calibration", |b| {
        b.iter_batched(
            || load(dir.path()),
            |(workspace, as_of)| {
                Session::new(workspace, as_of, ExecutionContext::sequential()).expect("calibration")
            },
            BatchSize::PerIteration,
        );
    });
}

#[expect(
    clippy::expect_used,
    reason = "計測の入力を読めなければ、計測を続けられない"
)]
fn workspace_load(c: &mut Criterion) {
    // ファイルの解析と、1 万件の取引の IR への変換（観測の列と公表時刻の解決を含む）。
    // 読み込みは誤りがあれば全体が失敗する。全件を読めたことは、準備と PV の検査が確かめる。
    let dir = desk(TRADES);
    c.bench_function("m0/workspace_load_10k", |b| {
        b.iter(|| Workspace::load(dir.path()).expect("desk workspace loads"));
    });
}

#[expect(
    clippy::expect_used,
    reason = "計測の入力から市場を作れなければ、計測を続けられない"
)]
fn session_new(c: &mut Criterion) {
    // 校正と、1 万件の取引の準備（fixing と決済の簡約、cashflow の計画）。読み込みは測らない。
    let dir = desk(TRADES);
    let mut group = c.benchmark_group("m0/session_new_10k");
    for threads in thread_counts() {
        let executor = ThreadExecutor::new(threads);
        // 同じスレッド数で作ったセッションが全件を準備できることを、測る前に一度確かめる。
        let (workspace, as_of) = load(dir.path());
        let session =
            Session::new(workspace, as_of, ExecutionContext::new(&executor)).expect("session");
        validate(&session, ExecutionContext::new(&executor), TRADES);
        drop(session);
        group.bench_function(format!("threads={threads}"), |b| {
            b.iter_batched(
                || load(dir.path()),
                |(workspace, as_of)| {
                    Session::new(workspace, as_of, ExecutionContext::new(&executor))
                        .expect("session")
                },
                BatchSize::PerIteration,
            );
        });
    }
    group.finish();
}

#[expect(
    clippy::expect_used,
    reason = "計測の入力から市場を作れなければ、計測を続けられない"
)]
fn pv_and_delta(c: &mut Criterion) {
    // 準備済みのセッションを一度だけ作り、PV と quote デルタの再計算を測る。
    let dir = desk(TRADES);
    let (workspace, as_of) = load(dir.path());
    let session =
        Session::new(workspace, as_of, ExecutionContext::sequential()).expect("desk session");
    // 各スレッド数で、PV が全件成功し quote デルタが取引を除かないことを、測る前に確かめる。
    for threads in thread_counts() {
        let executor = ThreadExecutor::new(threads);
        validate(&session, ExecutionContext::new(&executor), TRADES);
    }
    let mut group = c.benchmark_group("m0/pv_10k");
    for threads in thread_counts() {
        let executor = ThreadExecutor::new(threads);
        group.bench_function(format!("threads={threads}"), |b| {
            b.iter(|| black_box(session.pv(ExecutionContext::new(&executor)).expect("pv")));
        });
    }
    group.finish();
    let mut group = c.benchmark_group("m0/quote_delta_10k");
    for threads in thread_counts() {
        let executor = ThreadExecutor::new(threads);
        group.bench_function(format!("threads={threads}"), |b| {
            b.iter(|| {
                black_box(
                    session
                        .quote_delta(ExecutionContext::new(&executor))
                        .expect("quote delta"),
                )
            });
        });
    }
    group.finish();
}

#[expect(
    clippy::expect_used,
    reason = "計測の前提が崩れていれば、計測を続けられない"
)]
fn quote_update(c: &mut Criterion) {
    // 既存の 1 万件デスクで、変更 quote 1 本の再校正から PV とデルタまで通して測る。
    let dir = desk(TRADES);
    let (workspace, as_of) = load(dir.path());
    let threads = std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN);
    let executor = ThreadExecutor::new(threads);
    let execution = ExecutionContext::new(&executor);
    let session = Session::new(workspace, as_of, execution).expect("desk session");
    let overrides = session
        .quote_overrides([QuoteOverride::Set {
            id: "JPY.TONA.OIS.10Y".parse().expect("quote ID"),
            value: 0.014,
        }])
        .expect("quote override");
    validate(
        &session
            .with_quotes(&overrides, execution)
            .expect("quote update"),
        execution,
        TRADES,
    );
    c.bench_function("m1/quote_update_10k", |b| {
        b.iter(|| {
            let bumped = session
                .with_quotes(black_box(&overrides), execution)
                .expect("quote update");
            black_box(bumped.pv(execution).expect("pv"));
            black_box(bumped.quote_delta(execution).expect("quote delta"));
        });
    });
}

// criterion の生成する関数を公開しないよう、非公開のモジュールに置く。
mod harness {
    use criterion::Criterion;

    criterion::criterion_group! {
        name = benches;
        config = Criterion::default().sample_size(10);
        targets = super::calibration, super::workspace_load, super::session_new, super::pv_and_delta, super::quote_update
    }
}

criterion::criterion_main!(harness::benches);
`````

`````
crates/rhodium/benches/m1.rs
//! M1 の所要時間: 最初の検証データセットの校正と、M1 の商品を混ぜたデスクの準備、PV、quote デルタ、
//! quote の更新。
//!
//! `cargo bench -p rhodium --bench m1` で実行する。入力は `examples/m1-rates` で、デスクはその
//! 18 件の取引を ID を変えて繰り返したもの（同じ商品の繰り返しなので、観測の列の共有が最大に効く）。
//! 結果は日付、コミット、マシンとともに `docs/benchmarks/m1.md` に残す。
//! 計測の名前は `m1_rates/` で始め、M0 のデスクで M1-7 を測る `m0` の `m1/quote_update_10k` と区別する。

use std::{
    fs,
    hint::black_box,
    num::NonZeroUsize,
    path::{Path, PathBuf},
};

use criterion::{BatchSize, Criterion};
use rhodium::{
    prelude::{ExecutionContext, ThreadExecutor},
    session::{GroupOutcome, Session},
    workspace::Workspace,
};
use rhodium_core::time::Instant;
use rhodium_market::quote::QuoteOverride;
use tempfile::TempDir;

/// 1 万件のデスク。
const DESK: usize = 10_000;
/// 10 万件のデスク（性能の目安の「線形な取引 10 万件の PV と quote デルタ」）。
const LARGE_DESK: usize = 100_000;

#[expect(
    clippy::disallowed_macros,
    reason = "リポジトリのサンプルの場所をベンチのソースの位置から決める"
)]
fn sample() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/m1-rates")
}

#[expect(
    clippy::expect_used,
    reason = "計測の入力を作れなければ、計測を続けられない"
)]
fn copy_dir(from: &Path, to: &Path) {
    for entry in fs::read_dir(from).expect("sample directory") {
        let entry = entry.expect("sample entry");
        let target = to.join(entry.file_name());
        if entry.file_type().expect("entry type").is_dir() {
            fs::create_dir(&target).expect("directory");
            copy_dir(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), target).expect("file");
        }
    }
}

/// サンプルの取引ファイルの `[[trade]]` ごとの塊（文書の `version` の行を除く）。
#[expect(
    clippy::expect_used,
    reason = "計測の入力を作れなければ、計測を続けられない"
)]
fn trade_blocks() -> Vec<String> {
    let mut paths: Vec<_> = fs::read_dir(sample().join("trades"))
        .expect("trades directory")
        .map(|entry| entry.expect("trade file").path())
        .collect();
    paths.sort();
    paths
        .iter()
        .flat_map(|path| {
            let text = fs::read_to_string(path).expect("trade file");
            text.split("\n[[trade]]\n")
                .skip(1)
                .map(|block| format!("\n[[trade]]\n{block}"))
                .collect::<Vec<_>>()
        })
        .collect()
}

/// サンプルの市場と参照データに、取引を `trades` 件まで繰り返したデスクを置く。
#[expect(
    clippy::expect_used,
    reason = "計測の入力を作れなければ、計測を続けられない"
)]
fn desk(trades: usize) -> TempDir {
    let dir = tempfile::tempdir().expect("temporary directory");
    copy_dir(&sample(), dir.path());
    let trades_dir = dir.path().join("trades");
    fs::remove_dir_all(&trades_dir).expect("remove sample trades");
    fs::create_dir(&trades_dir).expect("trades directory");
    let blocks = trade_blocks();
    let mut text = String::from("version = 1\n");
    for (index, block) in blocks.iter().cycle().take(trades).enumerate() {
        // ID の行だけを変え、条件は同じ取引を繰り返す。
        let renamed = block.replacen("\nid = \"", &format!("\nid = \"D{index:06}-"), 1);
        text.push_str(&renamed);
    }
    fs::write(trades_dir.join("desk.toml"), text).expect("desk trades");
    dir
}

#[expect(
    clippy::expect_used,
    reason = "計測の入力を読めなければ、計測を続けられない"
)]
fn load(dir: &Path) -> (Workspace, Instant) {
    let workspace = Workspace::load(dir).expect("m1 workspace loads");
    let as_of = workspace
        .resolve_date("2026-11-24".parse().expect("valuation date"))
        .expect("valuation close");
    (workspace, as_of)
}

fn threads() -> NonZeroUsize {
    std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN)
}

/// 計測の前に、校正の全グループが解け、全取引の PV が成功し、デルタが取引を除かないことを確かめる。
#[expect(
    clippy::expect_used,
    reason = "計測の前提が崩れていれば、計測を続けられない"
)]
fn validate(session: &Session, execution: ExecutionContext<'_>, trades: usize) {
    assert!(
        session
            .calibrated_market()
            .groups()
            .iter()
            .all(|group| matches!(group.outcome(), GroupOutcome::Solved(_)))
    );
    let pv = session.pv(execution).expect("pv");
    assert_eq!(pv.len(), trades);
    assert!(pv.values().all(Result::is_ok));
    if trades > 0 {
        let delta = session.quote_delta(execution).expect("quote delta");
        assert!(delta.failures().is_empty(), "{:?}", delta.failures());
    }
}

#[expect(
    clippy::expect_used,
    reason = "計測の入力から市場を作れなければ、計測を続けられない"
)]
fn calibration(c: &mut Criterion) {
    // 取引の無いサンプルで Session::new を 1 スレッドで測ると、最初の検証データセットの全曲線の
    // 校正（9 グループ、TONA の階段と跳び、TIBOR 6M・3M、SOFR、FX ポイントと XCCY）だけになる。
    let dir = desk(0);
    let (workspace, as_of) = load(dir.path());
    let session = Session::new(workspace, as_of, ExecutionContext::sequential()).expect("session");
    validate(&session, ExecutionContext::sequential(), 0);
    c.bench_function("m1_rates/calibration", |b| {
        b.iter_batched(
            || load(dir.path()),
            |(workspace, as_of)| {
                Session::new(workspace, as_of, ExecutionContext::sequential()).expect("session")
            },
            BatchSize::PerIteration,
        );
    });
}

#[expect(
    clippy::expect_used,
    reason = "計測の入力から市場を作れなければ、計測を続けられない"
)]
fn desk_benchmarks(c: &mut Criterion, trades: usize, name: &str) {
    // 読み込み、準備（校正込み）、PV、quote デルタ、TONA の 10Y の変更からの再計算を、使える全スレッドで測る。
    let dir = desk(trades);
    let executor = ThreadExecutor::new(threads());
    let execution = ExecutionContext::new(&executor);
    c.bench_function(&format!("m1_rates/workspace_load_{name}"), |b| {
        b.iter(|| black_box(load(dir.path())));
    });
    c.bench_function(&format!("m1_rates/session_new_{name}"), |b| {
        b.iter_batched(
            || load(dir.path()),
            |(workspace, as_of)| Session::new(workspace, as_of, execution).expect("session"),
            BatchSize::PerIteration,
        );
    });
    let (workspace, as_of) = load(dir.path());
    let session = Session::new(workspace, as_of, execution).expect("session");
    validate(&session, execution, trades);
    c.bench_function(&format!("m1_rates/pv_{name}"), |b| {
        b.iter(|| black_box(session.value_portfolio(execution).expect("portfolio")));
    });
    c.bench_function(&format!("m1_rates/quote_delta_{name}"), |b| {
        b.iter(|| {
            black_box(
                session
                    .quote_delta_reporting(execution)
                    .expect("quote delta"),
            )
        });
    });
    let overrides = session
        .quote_overrides([QuoteOverride::Set {
            id: "JPY.TONA.OIS.10Y".parse().expect("quote ID"),
            value: 0.0133,
        }])
        .expect("quote override");
    validate(
        &session.with_quotes(&overrides, execution).expect("update"),
        execution,
        trades,
    );
    c.bench_function(&format!("m1_rates/quote_update_{name}"), |b| {
        b.iter(|| {
            let bumped = session
                .with_quotes(black_box(&overrides), execution)
                .expect("quote update");
            black_box(bumped.value_portfolio(execution).expect("portfolio"));
            black_box(
                bumped
                    .quote_delta_reporting(execution)
                    .expect("quote delta"),
            );
        });
    });
}

fn desk_10k(c: &mut Criterion) {
    desk_benchmarks(c, DESK, "10k");
}

fn desk_100k(c: &mut Criterion) {
    desk_benchmarks(c, LARGE_DESK, "100k");
}

// criterion の生成する関数を公開しないよう、非公開のモジュールに置く。
mod harness {
    use criterion::Criterion;

    criterion::criterion_group! {
        name = benches;
        config = Criterion::default().sample_size(10);
        targets = super::calibration, super::desk_10k, super::desk_100k
    }
}

criterion::criterion_main!(harness::benches);
`````

`````
crates/rhodium/examples/issue_180_workspace.rs
//! Issue #180 のデスク生成、読み込み・準備・最大 RSS の計測、前後比較のための数値出力。

#[allow(dead_code, reason = "計測では desk の締め日の定数を使わない")]
#[path = "../tests/support/desk.rs"]
mod desk;

use std::{hint::black_box, num::NonZeroUsize, path::Path, time::Instant};

use rhodium::{
    prelude::{Error, ExecutionContext, Result, ThreadExecutor},
    session::{GroupOutcome, Session},
    workspace::Workspace,
};

#[expect(
    clippy::disallowed_methods,
    clippy::print_stderr,
    reason = "計測バイナリの入口で引数を読み、工程ごとの実時間を測って出力する"
)]
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    let [mode, dir, count, threads] = args.as_slice() else {
        return Err(Error::InvalidInput(
            "usage: issue_180_workspace <generate|load|prepare|values> <directory> <trades> <threads>".into(),
        ));
    };
    let count: usize = count
        .to_str()
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| Error::InvalidInput("invalid trade count".into()))?;
    let threads: NonZeroUsize = threads
        .to_str()
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| Error::InvalidInput("invalid thread count".into()))?;
    let dir = Path::new(dir);
    let mode = mode.to_str().unwrap_or_default();
    if mode == "generate" {
        return desk::write(dir, count);
    }
    let executor = ThreadExecutor::new(threads);
    let execution = ExecutionContext::new(&executor);
    let start = Instant::now();
    let workspace = Workspace::load_with_execution(dir, execution)?;
    eprintln!("load_s\t{:.6}", start.elapsed().as_secs_f64());
    match mode {
        "load" => {
            black_box(workspace);
        }
        "prepare" | "values" => {
            let as_of = workspace.resolve_as_of(desk::AS_OF)?;
            let start = Instant::now();
            let session = Session::new(workspace, as_of, execution)?;
            eprintln!("prepare_s\t{:.6}", start.elapsed().as_secs_f64());
            assert_eq!(session.preparation_warnings().count(), count);
            assert!(
                session
                    .calibrated_market()
                    .groups()
                    .iter()
                    .all(|group| { matches!(group.outcome(), GroupOutcome::Solved(_)) })
            );
            if mode == "values" {
                print_values(&session, execution)?;
            }
            black_box(session);
        }
        _ => return Err(Error::InvalidInput("unknown benchmark mode".into())),
    }
    Ok(())
}

#[expect(
    clippy::print_stdout,
    reason = "変更前後で照合する数値を TSV で出力する"
)]
fn print_values(session: &Session, execution: ExecutionContext<'_>) -> Result<()> {
    for (id, pv) in session.pv(execution)? {
        for term in pv?.terms() {
            println!("pv\t{id}\t{}\t{:.17e}", term.ccy(), term.amount());
        }
    }
    let report = session.quote_delta(execution)?;
    assert!(report.failures().is_empty(), "{:?}", report.failures());
    for vector in report.by_ccy() {
        for (id, value) in vector.entries() {
            println!("delta\t{id}\t{}\t{value:.17e}", vector.ccy());
        }
    }
    Ok(())
}
`````

`````
crates/rhodium/src/io.rs
//! ファイル形式から、検証済みの値に変換する。

use rhodium_core::{
    Error, Result,
    calendar::{Calendar, CalendarSpec},
    refdata::RefData,
};
use rhodium_products::standard::StandardProduct;

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct CalendarDocument {
    version: u32,
    calendar: CalendarSpec,
}

impl CalendarDocument {
    fn resolve(self, refdata: &RefData) -> Result<Calendar> {
        if self.version != 1 {
            return Err(Error::InvalidData(format!(
                "unsupported calendar document version {}",
                self.version
            )));
        }
        refdata.calendar(&self.calendar)
    }
}

/// TOML のカレンダー指定を読み、参照データに照らして解決する。
///
/// # Errors
///
/// 形式、版、カレンダー式またはカレンダーへの参照が不正なら `InvalidData`。
pub fn read_calendar_toml(text: &str, refdata: &RefData) -> Result<Calendar> {
    toml::from_str::<CalendarDocument>(text)
        .map_err(|error| Error::InvalidData(error.to_string()))?
        .resolve(refdata)
}

/// JSON のカレンダー指定を読み、参照データに照らして解決する。
///
/// # Errors
///
/// 形式、版、カレンダー式またはカレンダーへの参照が不正なら `InvalidData`。
pub fn read_calendar_json(text: &str, refdata: &RefData) -> Result<Calendar> {
    // serde の構造体には配列でも値を渡せるため、ここでは文書の先頭を確認する。
    if !text.trim_start().starts_with('{') {
        return Err(Error::InvalidData(
            "calendar JSON document must be an object".into(),
        ));
    }
    serde_json::from_str::<CalendarDocument>(text)
        .map_err(|error| Error::InvalidData(error.to_string()))?
        .resolve(refdata)
}

/// TOML の商品テーブルを検証済みの標準商品に変換する。
///
/// # Errors
///
/// TOML の形式または商品の条件が不正なら `InvalidData`。
pub fn read_standard_product_toml(text: &str) -> Result<StandardProduct> {
    let mut value: toml::Value =
        toml::from_str(text).map_err(|error| Error::InvalidData(error.to_string()))?;
    normalize_datetimes(&mut value, "", false)?;
    value
        .try_into()
        .map_err(|error: toml::de::Error| Error::InvalidData(error.to_string()))
}

pub(crate) fn normalize_workspace_datetimes(value: &mut toml::Value) -> Result<()> {
    normalize_datetimes(value, "", true)
}

// ドメイン型の共通表現は文字列。日付と瞬間のフィールドだけでリテラルを変換する。
fn normalize_datetimes(value: &mut toml::Value, field: &str, workspace: bool) -> Result<()> {
    match value {
        toml::Value::Datetime(datetime) => {
            let date_field = matches!(
                field,
                "start"
                    | "end"
                    | "maturity"
                    | "payment"
                    | "reset"
                    | "fixing"
                    | "first_clearing"
                    | "final_clearing"
                    | "first_regular"
                    | "last_regular"
                    | "add"
                    | "remove"
            ) || (workspace
                && matches!(
                    field,
                    "trade_date"
                        | "from"
                        | "through"
                        | "since"
                        | "traded_on_or_after"
                        | "settlement"
                ));
            let instant_field = field == "settlement"
                || workspace
                    && matches!(
                        field,
                        "from"
                            | "through"
                            | "at"
                            | "effective_at"
                            | "pillars"
                            | "pillar"
                            | "breaks"
                            | "grid"
                            | "meetings"
                            | "dates"
                    );
            let text = match (datetime.date, datetime.time, datetime.offset) {
                (Some(date), None, None) if date_field => date.to_string(),
                (Some(_), Some(_), Some(_)) if instant_field => datetime.to_string(),
                _ => {
                    return Err(Error::InvalidData(format!(
                        "invalid TOML date or instant literal in {field}"
                    )));
                }
            };
            *value = toml::Value::String(text);
            Ok(())
        }
        toml::Value::Array(items) => items
            .iter_mut()
            .try_for_each(|item| normalize_datetimes(item, field, workspace)),
        toml::Value::Table(table) => table
            .iter_mut()
            .try_for_each(|(key, value)| normalize_datetimes(value, key, workspace)),
        _ => Ok(()),
    }
}
`````

`````
crates/rhodium/src/lib.rs
//! rhodium の公開ライブラリ API。

#![forbid(unsafe_code)]

// 公開 API を追加するまで、依存先を明示的に参照する。
use rhodium_xva as _;

pub mod io;
pub mod session;
pub mod workspace;

/// ライブラリ全体で共通に使う型。
///
/// 並列の実行器は呼び出し側が引数で渡す。逐次なら
/// [`ExecutionContext::sequential`](rhodium_core::parallel::ExecutionContext::sequential)、
/// 並列なら `ExecutionContext::new(&ThreadExecutor::new(workers))` を渡す。
pub mod prelude {
    pub use rhodium_core::{
        Error, Result,
        parallel::{ExecutionContext, ThreadExecutor},
    };
}
`````

`````
crates/rhodium/src/session.rs
//! 同じ評価時点で市場の校正と取引の準備を使い回す。

use std::{collections::BTreeMap, sync::Arc};

use rhodium_core::{
    Error, Result, Warning,
    id::TradeId,
    money::MoneyBag,
    parallel::{self, ExecutionContext, map_partitions},
    time::Instant,
};
use rhodium_market::quote::{QuoteId, QuoteOverride, QuoteOverrides};
use rhodium_pricing::{
    calibration::{CalibratedMarket, calibrate},
    context::Context,
    prepare::{PrepareInput, Prepared, evaluate, prepare},
    valuation,
};
use rhodium_risk::{
    delta::{self, DeltaReport},
    vector::RiskVector,
};

use crate::workspace::Workspace;

pub use rhodium_pricing::{
    calibration::GroupOutcome,
    valuation::{PortfolioValuation, Valuation},
};

/// ワークスペースの取引を一つの評価時点と校正済みの市場に結び付けた不変の値。
#[derive(Debug)]
pub struct Session {
    workspace: Workspace,
    calibrated: CalibratedMarket,
    trades: Arc<PreparedTrades>,
}

#[derive(Debug)]
struct PreparedTrades {
    prepared: Vec<Prepared>,
    // デルタが入力スライス内の位置で失敗を返すため、成功した取引の ID を同じ順で持つ。
    ids: Vec<TradeId>,
    failures: BTreeMap<TradeId, Error>,
}

/// quote デルタと、集計から除いた取引。
#[derive(Debug)]
pub struct QuoteDeltaReport {
    delta: DeltaReport,
    failures: BTreeMap<TradeId, Error>,
}

impl QuoteDeltaReport {
    /// 受渡通貨別、または報告通貨の quote デルタ。値は quote ごとの報告単位あたりの変化額。
    #[must_use]
    pub fn by_ccy(&self) -> &[RiskVector<QuoteId>] {
        self.delta.by_ccy()
    }

    /// 準備またはデルタの評価に失敗し、集計から除いた取引。
    #[must_use]
    pub const fn failures(&self) -> &BTreeMap<TradeId, Error> {
        &self.failures
    }
}

impl Session {
    /// 指定の瞬間で校正と準備を行う。取引ごとの準備エラーは結果に残す。
    ///
    /// 評価時点、fixing、担保条件を変更するときは、変更後の入力を渡してセッションを作る。
    ///
    /// 校正の後、取引の準備（fixing と決済の簡約、cashflow の計画）を `execution` の
    /// ワーカーに取引の区間で分けて行う。準備は取引ごとに独立なので、結果はスレッド数によらない。
    /// 取消と期限は校正の前と取引の境界で確かめ、止まれば取引の失敗にせず全体を `Execution` にする。
    /// 取引が 0 件でも、止まった実行は成功にしない。
    ///
    /// # Errors
    ///
    /// 市場構造や校正条件を組み立てられない場合に失敗する。
    /// 校正グループ単位の失敗は [`Self::calibrated_market`] で確認できる。
    /// 並列の実行が取消か期限で止まったか、実行器が失敗した場合は `Execution` か `Internal`。
    pub fn new(
        workspace: Workspace,
        as_of: Instant,
        execution: ExecutionContext<'_>,
    ) -> Result<Self> {
        let context = context(&workspace, as_of);
        // 校正は実行器を使わないので、取り消した実行で校正を始めないよう先に確かめる。
        parallel::check(execution)?;
        let calibrated = calibrate(&workspace.inputs.market, workspace.quotes(), &context)?;
        let structure = calibrated.market().structure();
        let trades: Vec<_> = workspace.trades().collect();
        let plans = map_partitions(execution, trades.len(), |range, child| {
            trades[range]
                .iter()
                .map(|&(_, trade, product, collateral)| {
                    parallel::check(child)?;
                    // 下ろした IR を使い回し、方式の表は取引の商品の種類で引く。
                    Ok(prepare(
                        PrepareInput::Trade { trade, product },
                        collateral,
                        &[],
                        structure,
                        &context,
                    ))
                })
                .collect::<Result<Vec<_>>>()
        })?
        .into_iter()
        .collect::<Result<Vec<_>>>()?;
        let mut prepared = Vec::new();
        let mut prepared_ids = Vec::new();
        let mut preparation_failures = BTreeMap::new();
        for (&(id, _, _, _), plan) in trades.iter().zip(plans.into_iter().flatten()) {
            match plan {
                Ok(plan) => {
                    prepared.push(plan);
                    prepared_ids.push(id.clone());
                }
                Err(error) => {
                    preparation_failures.insert(id.clone(), error);
                }
            }
        }
        Ok(Self {
            workspace,
            calibrated,
            trades: Arc::new(PreparedTrades {
                prepared,
                ids: prepared_ids,
                failures: preparation_failures,
            }),
        })
    }

    /// このセッションの市場の binding に照らして quote の上書きを作る。
    ///
    /// # Errors
    ///
    /// [`QuoteOverrides::new`] と同じ。
    pub fn quote_overrides(
        &self,
        operations: impl IntoIterator<Item = QuoteOverride>,
    ) -> Result<QuoteOverrides> {
        QuoteOverrides::new(operations, self.workspace.inputs.market.quote_kinds())
    }

    /// quote を差し替えて影響を受けるグループを再校正し、取引の準備を共有したセッションを返す。
    ///
    /// # Errors
    ///
    /// 上書きがこの市場やスナップショットに適用できない場合、市場を組み立てられない場合、
    /// 取消または期限で停止した場合。グループの失敗は [`Self::calibrated_market`] で確認できる。
    pub fn with_quotes(
        &self,
        overrides: &QuoteOverrides,
        execution: ExecutionContext<'_>,
    ) -> Result<Self> {
        parallel::check(execution)?;
        let overrides = overrides.for_kinds(self.workspace.inputs.market.quote_kinds())?;
        let (quotes, _) = self.workspace.quotes().with_overrides(&overrides)?;
        let calibrated = self.calibrated.update(&quotes)?;
        parallel::check(execution)?;
        Ok(Self {
            workspace: self.workspace.with_quote_set(quotes),
            calibrated,
            trades: Arc::clone(&self.trades),
        })
    }

    /// 市場と取引を準備した評価の瞬間。
    #[must_use]
    pub fn as_of(&self) -> Instant {
        self.calibrated.market().structure().as_of()
    }

    /// このセッションの入力スナップショット。
    #[must_use]
    pub const fn workspace(&self) -> &Workspace {
        &self.workspace
    }

    /// 利用できる曲線とグループごとの校正結果。
    #[must_use]
    pub const fn calibrated_market(&self) -> &CalibratedMarket {
        &self.calibrated
    }

    /// 取引の準備で発生した警告を取引 ID の順に返す。
    pub fn preparation_warnings(&self) -> impl Iterator<Item = (&TradeId, &[Warning])> {
        self.trades
            .ids
            .iter()
            .zip(self.trades.prepared.iter().map(Prepared::warnings))
    }

    /// 取引ごとの PV を受渡通貨別に返す。準備と評価の失敗も取引 ID で報告する。
    ///
    /// 取引を `execution` のワーカーに区間で分けて評価する。取引ごとの PV は分け方によらない。
    /// 取消と期限は取引の境界で確かめ、止まれば取引の失敗にせず全体を `Execution` にする。
    ///
    /// # Errors
    ///
    /// 並列の実行が取消か期限で止まったか、実行器が失敗した場合は `Execution` か `Internal`。
    /// 取引ごとの評価の失敗は呼び出し全体のエラーにせず、その取引の結果に残す。
    pub fn pv(
        &self,
        execution: ExecutionContext<'_>,
    ) -> Result<BTreeMap<TradeId, Result<MoneyBag>>> {
        let market = self.calibrated.market();
        let values = map_partitions(execution, self.trades.prepared.len(), |range, child| {
            self.trades.prepared[range]
                .iter()
                .map(|plan| {
                    parallel::check(child)?;
                    Ok(evaluate(plan, market))
                })
                .collect::<Result<Vec<_>>>()
        })?
        .into_iter()
        .collect::<Result<Vec<_>>>()?;
        Ok(self
            .trades
            .ids
            .iter()
            .cloned()
            .zip(values.into_iter().flatten())
            .chain(
                self.trades
                    .failures
                    .iter()
                    .map(|(id, error)| (id.clone(), Err(error.clone()))),
            )
            .collect())
    }

    /// 取引 ID 順の評価結果、報告通貨での合計、除外した取引を返す。
    ///
    /// 取引の準備と受渡通貨別の評価は [`Self::pv`] と同じ計画を使う。報告通貨への
    /// 換算に失敗した取引も、除外一覧と取引ごとの結果に含める。
    ///
    /// # Errors
    ///
    /// [`valuation::value_portfolio`] と同じ。
    pub fn value_portfolio(&self, execution: ExecutionContext<'_>) -> Result<PortfolioValuation> {
        let mut prepared: Vec<_> = self
            .trades
            .ids
            .iter()
            .cloned()
            .zip(self.trades.prepared.iter().map(Ok))
            .chain(
                self.trades
                    .failures
                    .iter()
                    .map(|(id, error)| (id.clone(), Err(error.clone()))),
            )
            .collect();
        if !self.trades.failures.is_empty() {
            prepared.sort_by(|(left, _), (right, _)| left.cmp(right));
        }
        valuation::value_portfolio(
            &prepared,
            &self.calibrated,
            &context(&self.workspace, self.as_of()),
            execution,
        )
    }

    /// 全取引の quote デルタを受渡通貨別に集計する。失敗した取引は集計から除く。
    ///
    /// 取引の θ 勾配は `execution` のワーカーごとに一つの tape で求め、通貨別に足してから
    /// 校正の随伴を一度だけ解く（[`delta::quote_delta`]）。
    ///
    /// # Errors
    ///
    /// 勾配の集計や校正の随伴に失敗した場合。
    /// 並列の実行が取消か期限で止まったか、実行器が失敗した場合は `Execution` か `Internal`
    /// （止まった取引を集計から除いた結果は返さない）。
    pub fn quote_delta(&self, execution: ExecutionContext<'_>) -> Result<QuoteDeltaReport> {
        let delta = delta::quote_delta(&self.trades.prepared, &self.calibrated, execution)?;
        self.report_quote_delta(delta)
    }

    /// 価格設定の報告通貨で全取引の quote デルタを集計する。
    ///
    /// 評価時点の FX で換算し、その FX 自身の quote 感応度も含める。
    /// 準備または評価に失敗した取引は取引 ID 付きで返す。
    ///
    /// # Errors
    ///
    /// [`Self::quote_delta`] のエラーに加え、換算に必要な FX 経路かデータの欠損は
    /// `MissingData`、換算の数値計算の失敗は `Numerical`。
    pub fn quote_delta_reporting(
        &self,
        execution: ExecutionContext<'_>,
    ) -> Result<QuoteDeltaReport> {
        let delta = delta::quote_delta_in(
            &self.trades.prepared,
            &self.calibrated,
            self.workspace.inputs.config.reporting_currency(),
            execution,
        )?;
        self.report_quote_delta(delta)
    }

    fn report_quote_delta(&self, delta: DeltaReport) -> Result<QuoteDeltaReport> {
        let mut failures = self.trades.failures.clone();
        for (index, error) in delta.failures() {
            let id = self.trades.ids.get(*index).ok_or_else(|| {
                Error::Internal(format!("unknown prepared trade position {index}"))
            })?;
            failures.insert(id.clone(), error.clone());
        }
        Ok(QuoteDeltaReport { delta, failures })
    }
}

fn context(workspace: &Workspace, as_of: Instant) -> Context<'_> {
    Context {
        as_of,
        refdata: &workspace.inputs.refdata,
        fixings: &workspace.inputs.fixings,
        config: &workspace.inputs.config,
    }
}
`````

`````
crates/rhodium/src/workspace.rs
//! JPY・USD 金利評価に使うワークスペースの入力。

mod csv_files;
mod input;
mod market_file;
mod refdata_files;
mod trade_files;

use std::{collections::BTreeMap, path::Path, sync::Arc};

use rhodium_core::{
    Error, Result,
    date::Date,
    fixing::FixingStore,
    id::{BookId, TradeId},
    parallel::{self, ExecutionContext},
    refdata::RefData,
    time::{Instant, LocalDateTime},
};
use rhodium_market::quote::QuoteSet;
use rhodium_pricing::{
    config::PricingConfig, definition::MarketDefinition, resolve::CollateralKey,
};
use rhodium_products::{ir::Product, trade::Trade};

use self::{input::Source, trade_files::LoadedTrade};

/// 一つの評価に使う、参照先を解決済みの入力スナップショット。
#[derive(Debug)]
pub struct Workspace {
    pub(crate) inputs: Arc<WorkspaceInputs>,
    quotes: QuoteSet,
}

/// quote を差し替えても共有できる、解決済みのワークスペース入力。
#[derive(Debug)]
pub(crate) struct WorkspaceInputs {
    pub(crate) refdata: RefData,
    pub(crate) fixings: FixingStore,
    pub(crate) market: MarketDefinition,
    pub(crate) config: PricingConfig,
    trades: BTreeMap<TradeId, LoadedTrade>,
}

impl Workspace {
    /// ディレクトリにある TOML と CSV を読み、参照データと取引の参照先を解決する。
    ///
    /// # Errors
    ///
    /// 必須ファイルが読めないか、形式、内容、データ間の参照が不正なら `InvalidData`。
    pub fn load(dir: impl AsRef<Path>) -> Result<Self> {
        Self::load_with_execution(dir, ExecutionContext::sequential())
    }

    /// ワークスペースを読み、指定した実行器で取引を検証済み IR に変換する。
    ///
    /// # Errors
    ///
    /// [`Self::load`] の入力エラーに加え、取消や期限切れなどの実行エラーを返す。
    pub fn load_with_execution(
        dir: impl AsRef<Path>,
        execution: ExecutionContext<'_>,
    ) -> Result<Self> {
        parallel::check(execution)?;
        let dir = dir.as_ref();
        let refdata = refdata_files::load(&dir.join("refdata"))?;
        let market = market_file::load(dir.join("market/definition.toml"), &refdata)?;
        let source = Source::read(dir.join("config/pricing.toml"))?;
        let config: PricingConfig = source.decode()?;
        config
            .validate_references(&refdata, &market)
            .map_err(|error| source.invalid(error))?;
        let quotes = csv_files::quotes(&dir.join("quotes/snapshot.csv"))?;
        let fixings = csv_files::fixings(&dir.join("fixings"), &refdata)?;
        let trades = trade_files::load(dir, &refdata, execution)?;
        Ok(Self {
            inputs: Arc::new(WorkspaceInputs {
                refdata,
                fixings,
                market,
                config,
                trades,
            }),
            quotes,
        })
    }

    /// このワークスペースの quote スナップショット。
    #[must_use]
    pub const fn quotes(&self) -> &QuoteSet {
        &self.quotes
    }

    pub(crate) fn with_quote_set(&self, quotes: QuoteSet) -> Self {
        Self {
            inputs: Arc::clone(&self.inputs),
            quotes,
        }
    }

    /// オフセット付き日時、または末尾に `[地域名]` を指定した日時を UTC に解決する。
    ///
    /// 地域名とオフセットを併記した場合は、指定のオフセットと地域の規則を照合する。
    ///
    /// # Errors
    ///
    /// 日時や地域が不正、オフセットが規則と異なる、または瞬間を一意に決められない場合。
    pub fn resolve_as_of(&self, text: &str) -> Result<Instant> {
        resolve_as_of(text, &self.inputs.refdata)
    }

    /// 日付を価格設定の締め時刻と地域で UTC に解決する。
    ///
    /// # Errors
    ///
    /// 締め時刻が未設定か、その日の現地時刻を一意に解決できない場合。
    pub fn resolve_date(&self, date: Date) -> Result<Instant> {
        let close = self.inputs.config.valuation_close().ok_or_else(|| {
            Error::InvalidInput("--date requires valuation_close in config/pricing.toml".into())
        })?;
        self.inputs
            .refdata
            .zone(&close.zone)?
            .resolve_input(close.on_date(date), None)
    }

    /// 取引のブックを返す。指定の取引が無ければ `None`。
    #[must_use]
    pub fn book(&self, trade: &TradeId) -> Option<&BookId> {
        self.inputs
            .trades
            .get(trade)
            .map(|loaded| &loaded.trade.book)
    }

    pub(crate) fn trades(
        &self,
    ) -> impl Iterator<Item = (&TradeId, &Trade, &Product, &CollateralKey)> {
        self.inputs
            .trades
            .iter()
            .map(|(id, loaded)| (id, &loaded.trade, &loaded.product, &loaded.collateral))
    }
}

fn resolve_as_of(text: &str, refdata: &RefData) -> Result<Instant> {
    let Some((datetime, id)) = text
        .strip_suffix(']')
        .and_then(|text| text.rsplit_once('['))
    else {
        return text.parse();
    };
    let zone = refdata.zone(&id.parse()?)?;
    if let Ok(at) = datetime.parse::<Instant>() {
        // オフセット付きの文字列を解決した瞬間を現地日時に戻し、元の表示と照合する。
        let (date, time) = datetime
            .split_once(['T', 't', ' '])
            .ok_or_else(|| Error::InvalidInput("expected a date and time".into()))?;
        let time = time.split(['+', '-', 'Z', 'z']).next().unwrap_or_default();
        let local = LocalDateTime::new(date.parse()?, time.parse()?);
        if zone.to_local(at)? != local {
            return Err(Error::InvalidInput(
                "UTC offset does not match zone rules".into(),
            ));
        }
        return Ok(at);
    }
    zone.resolve_input(datetime.parse()?, None)
}
`````

`````
crates/rhodium/src/workspace/csv_files.rs
use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};

use csv::StringRecord;
use rhodium_core::{
    Error, Result,
    date::Date,
    fixing::{FixingStore, FixingStoreParts, ObservationPoint},
    id::IndexId,
    refdata::{IndexKind, PublicationRule, RefData},
};
use rhodium_market::quote::{QuoteId, QuoteSet};

use super::{
    input::{add_unique, file_id, files, finish, invalid},
    resolve_as_of,
};

pub(super) fn holidays(path: &Path) -> Result<Vec<Date>> {
    let mut holidays = BTreeSet::new();
    rows(path, &["date"], |record| {
        let date = record[0]
            .parse()
            .map_err(|error: Error| error.context("date"))?;
        if !holidays.insert(date) {
            return Err(Error::InvalidData(format!(
                "date: duplicate holiday {date}"
            )));
        }
        Ok(())
    })?;
    Ok(holidays.into_iter().collect())
}

pub(super) fn quotes(path: &Path) -> Result<QuoteSet> {
    let mut quotes = BTreeMap::new();
    rows(path, &["quote_id", "value"], |record| {
        let id: QuoteId = record[0]
            .parse()
            .map_err(|error: Error| error.context("quote_id"))?;
        let value = parse_value(&record[1])?;
        add_unique(&mut quotes, id, value).map_err(|error| error.context("quote_id"))
    })?;
    QuoteSet::new(quotes).map_err(|error| invalid(path, error))
}

// 指標ごとのファイルを読む。列名と観測点の型（日付、瞬間、月）は指標の公表規則で決まる。
pub(super) fn fixings(dir: &Path, refdata: &RefData) -> Result<FixingStore> {
    let mut parts = FixingStoreParts::default();
    for path in files(dir, "csv")? {
        let index: IndexId = file_id(&path)?;
        let definition = refdata
            .index(&index)
            .map_err(|error| invalid(&path, error))?;
        if matches!(
            definition.kind,
            IndexKind::CreditEvent { .. } | IndexKind::Recovery { .. }
        ) {
            return Err(invalid(
                &path,
                format!("{index} does not have numeric fixings"),
            ));
        }
        let field = match definition.publication {
            PublicationRule::Daily(_) => "date",
            PublicationRule::Intraday { .. } => "at",
            PublicationRule::Monthly { .. } => "month",
        };
        let mut points = BTreeSet::new();
        rows(&path, &[field, "value"], |record| {
            let point = match definition.publication {
                PublicationRule::Daily(_) => record[0].parse().map(ObservationPoint::Date),
                PublicationRule::Intraday { .. } => {
                    resolve_as_of(&record[0], refdata).map(ObservationPoint::Instant)
                }
                PublicationRule::Monthly { .. } => record[0].parse().map(ObservationPoint::Month),
            }
            .map_err(|error: Error| error.context(field))?;
            let value = parse_value(&record[1])?;
            if !points.insert(point) {
                return Err(Error::InvalidData(format!(
                    "{field}: duplicate fixing for {index} on {}",
                    &record[0]
                )));
            }
            parts.published.push((index.clone(), point, value));
            Ok(())
        })?;
    }
    FixingStore::new(refdata, parts).map_err(|error| invalid(dir, error))
}

// ヘッダーの完全一致を求めてから各レコードを解析する。列数の違いと値の誤りは
// 後続の行も調べ、行番号付きでまとめて返す。
fn rows(
    path: &Path,
    headers: &[&str],
    mut parse: impl FnMut(&StringRecord) -> Result<()>,
) -> Result<()> {
    let mut reader = csv::Reader::from_path(path).map_err(|error| invalid(path, error))?;
    let header = reader.headers().map_err(|error| invalid(path, error))?;
    if !header.iter().eq(headers.iter().copied()) {
        return Err(invalid(
            path,
            format!("line 1: expected CSV header {}", headers.join(",")),
        ));
    }
    let mut errors = Vec::new();
    for row in reader.records() {
        let record = match row {
            Ok(record) => record,
            Err(error) => {
                errors.push(invalid(path, error));
                continue;
            }
        };
        if let Err(error) = parse(&record) {
            let line = record.position().map_or(0, csv::Position::line);
            errors.push(invalid(path, format!("line {line}: {error}")));
        }
    }
    finish(errors)
}

pub(super) fn parse_value(text: &str) -> Result<f64> {
    // 接尾辞付きの値も十進で尺度を変えてから一度だけ f64 に丸める。
    let scaled;
    let decimal = if let Some(number) = text.strip_suffix("bp") {
        scaled = format!("{number}e-4");
        &scaled
    } else if let Some(number) = text.strip_suffix('%') {
        scaled = format!("{number}e-2");
        &scaled
    } else {
        text
    };
    decimal
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .ok_or_else(|| {
            Error::InvalidData(format!("value: expected a finite decimal, got {text:?}"))
        })
}
`````

`````
crates/rhodium/src/workspace/input.rs
use std::{
    collections::{BTreeMap, btree_map::Entry},
    fmt::Display,
    fs,
    path::{Path, PathBuf},
};

use rhodium_core::{Error, Result};
use serde::de::DeserializeOwned;
use toml::{Spanned, Value};

use crate::io::normalize_workspace_datetimes;

// 一つの TOML 文書の本文とパス。行の位置付きの診断を作るため、復号の後も本文を保持する。
pub(super) struct Source {
    path: PathBuf,
    text: String,
}

impl Source {
    pub(super) fn read(path: PathBuf) -> Result<Self> {
        let text = fs::read_to_string(&path).map_err(|error| invalid(&path, error))?;
        Ok(Self { path, text })
    }

    pub(super) fn decode<T: DeserializeOwned>(&self) -> Result<T> {
        toml::from_str(&self.text).map_err(|error| self.invalid(error))
    }

    // 文書の版を照合する。1.0 までの文書の版は 1 だけ（設計書「ファイルの形式」）。
    pub(super) fn version(&self, version: Spanned<u32>) -> Result<()> {
        let start = version.span().start;
        let number = version.into_inner();
        if number != 1 {
            return Err(self.at(start, "version", format!("unsupported version {number}")));
        }
        Ok(())
    }

    pub(super) fn rows<T: DeserializeOwned>(
        &self,
        field: &str,
        rows: Vec<Spanned<Value>>,
        mut insert: impl FnMut(T) -> Result<()>,
    ) -> Result<()> {
        self.entries(
            field,
            rows.into_iter()
                .enumerate()
                .map(|(i, row)| (i.to_string(), row)),
            |_, row| insert(row),
        )
    }

    // 各行を独立に復号して登録する。失敗した行があっても残りの行を調べ、
    // 同じ文書のエラーを行の位置と項目名付きでまとめて返す。
    pub(super) fn entries<T: DeserializeOwned>(
        &self,
        field: &str,
        rows: impl IntoIterator<Item = (String, Spanned<Value>)>,
        mut insert: impl FnMut(&str, T) -> Result<()>,
    ) -> Result<()> {
        let mut errors = Vec::new();
        for (key, row) in rows {
            let start = row.span().start;
            let decoded = Self::row(row.into_inner()).and_then(|row| insert(&key, row));
            if let Err(error) = decoded {
                errors.push(self.at(start, &format!("{field}[{key}]"), error));
            }
        }
        finish(errors)
    }

    // ドメイン型は日付と瞬間を文字列で受け取るため、TOML の日時リテラルを先に文字列へ直す。
    pub(super) fn row<T: DeserializeOwned>(mut value: Value) -> Result<T> {
        normalize_workspace_datetimes(&mut value)?;
        value
            .try_into()
            .map_err(|error| Error::InvalidData(error.to_string()))
    }

    pub(super) fn invalid(&self, error: impl Display) -> Error {
        invalid(&self.path, error)
    }

    // toml の span はバイト位置なので、1 始まりの行と文字単位の列に直してエディタで辿れる形にする。
    pub(super) fn at(&self, offset: usize, field: &str, error: impl Display) -> Error {
        let before = &self.text[..offset];
        let line = before.bytes().filter(|&byte| byte == b'\n').count() + 1;
        let column = before
            .rsplit('\n')
            .next()
            .map_or(1, |line| line.chars().count() + 1);
        Error::InvalidData(format!(
            "{}:{line}:{column}: {field}: {error}",
            self.path.display()
        ))
    }
}

pub(super) fn invalid(path: &Path, error: impl Display) -> Error {
    Error::InvalidData(format!("{}: {error}", path.display()))
}

// 集めたエラーを改行区切りの一つの InvalidData にする。InvalidData の接頭辞は外側で一度だけ付ける。
pub(super) fn finish(errors: Vec<Error>) -> Result<()> {
    if errors.is_empty() {
        return Ok(());
    }
    Err(Error::InvalidData(
        errors
            .into_iter()
            .map(|error| match error {
                Error::InvalidData(message) => message,
                error => error.to_string(),
            })
            .collect::<Vec<_>>()
            .join("\n"),
    ))
}

// 同じ種類のデータがファイルをまたいで同じ ID を持つことを拒否する。
pub(super) fn add_unique<K: Ord + Display, V>(
    map: &mut BTreeMap<K, V>,
    id: K,
    value: V,
) -> Result<()> {
    match map.entry(id) {
        Entry::Vacant(entry) => {
            entry.insert(value);
            Ok(())
        }
        Entry::Occupied(entry) => Err(Error::InvalidData(format!("duplicate ID {}", entry.key()))),
    }
}

// 任意のディレクトリ（trades/、fixings/ など）は無ければ空として扱う。
// 読み込み順とエラーの順を環境に依らず決めるため、パスの順に並べる。
pub(super) fn files(dir: &Path, extension: &str) -> Result<Vec<PathBuf>> {
    let entries = match fs::read_dir(dir) {
        Ok(entries) => entries,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(Vec::new()),
        Err(error) => return Err(invalid(dir, error)),
    };
    let mut paths = entries
        .map(|entry| {
            entry
                .map(|entry| entry.path())
                .map_err(|error| invalid(dir, error))
        })
        .collect::<Result<Vec<_>>>()?;
    paths.retain(|path| path.extension().is_some_and(|ext| ext == extension));
    paths.sort_unstable();
    Ok(paths)
}

// 休日や fixing のファイルは、拡張子を除いたファイル名を対象の ID として読む。
pub(super) fn file_id<T: std::str::FromStr>(path: &Path) -> Result<T>
where
    T::Err: Display,
{
    path.file_stem()
        .and_then(|name| name.to_str())
        .ok_or_else(|| invalid(path, "file name must contain a UTF-8 ID"))?
        .parse()
        .map_err(|error| invalid(path, error))
}
`````

`````
crates/rhodium/src/workspace/market_file.rs
use std::{collections::BTreeMap, path::PathBuf};

use rhodium_core::{
    Error, Result,
    calendar::CalendarSpec,
    currency::Currency,
    id::CalendarId,
    refdata::RefData,
    time::{Instant, TzId},
};
use rhodium_market::{
    curve::{CurveExpr, CurveId, CurveRole, CurveSpec, Extrapolation},
    quote::{QuoteId, QuoteKind},
    structure::{BlockId, ParamId},
};
use rhodium_pricing::{
    definition::{
        CalibrationTarget, ElementSelector, FxDefinition, InstrumentBinding, LgmDefinition,
        MarketDefinition, QuoteTransform, StepPenalty,
    },
    resolve::CollateralKey,
};
use rhodium_products::standard::InstrumentTemplate;
use serde::Deserialize;
use toml::{Spanned, Value};

use super::{
    csv_files,
    input::{Source, finish},
};

// 曲線と binding は行ごとの診断を返せるよう、位置付きの値として受けてから個別に復号する。
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MarketDocument {
    version: Spanned<u32>,
    valuation_calendar: CalendarId,
    valuation_zone: TzId,
    curves: BTreeMap<String, Spanned<Value>>,
    instruments: Vec<Spanned<Value>>,
    #[serde(default)]
    fx: Vec<Spanned<Value>>,
    #[serde(default)]
    lgm: BTreeMap<String, Spanned<Value>>,
    #[serde(default)]
    step_penalties: Vec<Spanned<Value>>,
    #[serde(default)]
    gauss_newton: Vec<Spanned<Value>>,
}

// 曲線は明示したピラーの節点（`pillars`）か、部品の式（`expr`、`CurveExpr` の serde の表現）で書く。
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CurveInput {
    currency: Currency,
    roles: Vec<CurveRole>,
    #[serde(default)]
    pillars: Option<Vec<Instant>>,
    #[serde(default)]
    extrapolation: Option<Extrapolation>,
    #[serde(default)]
    expr: Option<CurveExpr>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InstrumentInput {
    quote: QuoteId,
    kind: QuoteKind,
    target: Target,
    #[serde(default)]
    template: Option<InstrumentTemplate>,
    #[serde(default)]
    transform: QuoteTransform,
    #[serde(default)]
    collateral: Option<CollateralKey>,
    determines: Determines,
    /// 最小二乗の当てはめの許容幅（quote の単位。`"0.5bp"` などの接尾辞付きか十進）。
    #[serde(default)]
    tolerance: Option<Value>,
}

#[derive(Deserialize)]
#[serde(rename_all = "snake_case")]
enum Target {
    ParRate,
    Price,
    DirectMark,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Determines {
