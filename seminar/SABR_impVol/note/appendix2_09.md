
    /// この商品を解決または簡約するときに生じた警告。
    #[must_use]
    pub fn warnings(&self) -> &[Warning] {
        &self.warnings
    }

    /// 簡約で観測期間の一部を市場から予想する場合に、その開始の瞬間を返す。
    #[must_use]
    pub fn forecast_from(&self, index: &IndexId, point: ObservationPoint) -> Option<Instant> {
        self.forecasts.get(index)?.get(&point).copied()
    }
}

fn reachable_contracts(nodes: &[ContractNode], root: Ix) -> Vec<bool> {
    let mut reached = vec![false; root as usize + 1];
    reached[root as usize] = true;
    for pos in (0..reached.len()).rev() {
        if reached[pos] {
            nodes[pos].visit_children(|ix| reached[ix as usize] = true);
        }
    }
    reached
}

/// 受渡の瞬間が評価時点と等しいときに、受渡を現在価値に含めるか。
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SettlementBoundary {
    /// 決済済みとして、評価時点ちょうどの受渡を除く。
    #[default]
    Exclude,
    /// 評価時点ちょうどの受渡を含める。
    Include,
}

impl SettlementBoundary {
    /// この規則と受渡の瞬間から、評価に含めるかを決める。
    #[must_use]
    pub fn includes(self, as_of: Instant, settlement: Instant) -> bool {
        settlement > as_of || (settlement == as_of && self == Self::Include)
    }
}
`````

`````
crates/rhodium-products/src/ir/amendment.rs
use std::collections::{HashMap, hash_map::Entry};

use rhodium_core::{Error, Result, time::Instant};

use super::{
    Contract,
    nodes::{ContractNode, ExprNode, Ix, MAX_NODES, Nary, Number, node_limit},
    reachable_contracts,
};
use crate::trade::ExerciseRecord;

/// 取引の行使記録を `(label, at)` で引く表。
///
/// 値は選んだ選択肢の番号と、到達可能な `Choose` に一致したかの印。`ProductSpec` を IR にするとき、
/// 変更前の商品の判断を確定させるのに使い、最後にどの判断にも一致しなかった記録を `InvalidData` にする。
pub(crate) struct Exercises<'a> {
    records: HashMap<&'a str, HashMap<Instant, (usize, bool)>>,
}

impl<'a> Exercises<'a> {
    pub(crate) fn new(records: &'a [ExerciseRecord]) -> Result<Self> {
        let mut result = Self {
            records: HashMap::new(),
        };
        for record in records {
            if record.effective_at > record.at
                || result
                    .records
                    .entry(&record.label)
                    .or_default()
                    .insert(record.at, (record.option, false))
                    .is_some()
            {
                return Err(invalid_record(&record.label, record.at));
            }
        }
        Ok(result)
    }

    /// `(label, at)` の判断で記録が選んだ選択肢の番号。
    fn option(&self, label: &str, at: Instant) -> Option<usize> {
        self.records.get(label)?.get(&at).map(|&(option, _)| option)
    }

    /// 契約の到達可能な `Choose` に一致した記録へ印を付け、選択肢の番号が範囲内かを確かめる。
    pub(crate) fn observe(&mut self, contract: &Contract) -> Result<()> {
        if self.records.is_empty() {
            return Ok(());
        }
        for (node, needed) in contract
            .contracts
            .iter()
            .zip(reachable_contracts(&contract.contracts, contract.root))
        {
            if let (
                ContractNode::Choose {
                    label, at, options, ..
                },
                true,
            ) = (node, needed)
                && let Some((option, matched)) = self
                    .records
                    .get_mut(label.as_str())
                    .and_then(|times| times.get_mut(at))
            {
                if *option >= options.len() {
                    return Err(invalid_record(label, *at));
                }
                *matched = true;
            }
        }
        Ok(())
    }

    /// どの商品の判断にも一致しなかった記録があれば、取引のデータの誤りとして返す。
    pub(crate) fn finish(&self) -> Result<()> {
        // Map の順序によらず、同じ不正な記録を報告する。
        if let Some((label, at)) = self
            .records
            .iter()
            .flat_map(|(&label, times)| {
                times
                    .iter()
                    .filter_map(move |(&at, &(_, matched))| (!matched).then_some((label, at)))
            })
            .min()
        {
            return Err(invalid_record(label, at));
        }
        Ok(())
    }
}

fn invalid_record(label: &str, at: Instant) -> Error {
    Error::InvalidData(format!("invalid exercise record {label:?} at {at}"))
}

fn missing_record(from: Instant, label: &str, at: Instant) -> Error {
    Error::Unsupported(format!(
        "replacement at {from} requires an exercise record for {label:?} at {at}"
    ))
}

impl Contract {
    /// 瞬間 `from` 以後に到達しうるすべての `Pay` の最終額に `factor` を掛ける（一部解約）。
    ///
    /// 丸めや上限を適用した後の額を縮めるので、元本を先に縮めた場合とは結果が異なりうる。
    /// `factor` が 0 の全解約では、`from` 以後の受渡と判断を契約から取り除く。
    pub(crate) fn scale_from(mut self, from: Instant, factor: f64) -> Result<Self> {
        if !factor.is_finite() || !(0.0..1.0).contains(&factor) {
            return Err(Error::InvalidData(
                "amendment scale must be in [0, 1)".into(),
            ));
        }
        if factor == 0.0 {
            return self.retain(from, Side::Before, None);
        }
        // 係数の定数を一つ足し、同じ額の式を共有する受渡には同じ積の式を共有させる。
        // 式の arena は子が親より前に並ぶので、末尾に足した積はその順序を保つ。
        let constant = Ix::try_from(self.exprs.len()).map_err(|_| node_limit())?;
        self.exprs.push(ExprNode::Const(Number::new(factor)));
        let mut amounts = HashMap::new();
        let reached = reachable_contracts(&self.contracts, self.root);
        for (node, needed) in self.contracts.iter_mut().zip(reached) {
            if let (ContractNode::Pay { amount, at, .. }, true) = (node, needed)
                && *at >= from
            {
                *amount = match amounts.entry(*amount) {
                    Entry::Occupied(entry) => *entry.get(),
                    Entry::Vacant(entry) => {
                        let ix = Ix::try_from(self.exprs.len()).map_err(|_| node_limit())?;
                        self.exprs
                            .push(ExprNode::Nary(Nary::Mul, vec![*amount, constant]));
                        *entry.insert(ix)
                    }
                };
            }
        }
        self.compact()
    }

    /// `from` より前は `self`、`from` 以後は `new` に従う契約を作る（条件変更）。
    ///
    /// 境界をまたぐ旧契約の判断は、分ける前に行使記録で選択肢を確定させる。
    /// 記録が無ければ、どちらの受渡を残すか決まらないので `Unsupported` を返す。
    pub(crate) fn replace_from(
        self,
        new: Self,
        from: Instant,
        exercises: &Exercises<'_>,
    ) -> Result<Self> {
        let before = self.retain(from, Side::Before, Some(exercises))?;
        let after = new.retain(from, Side::After, Some(exercises))?;
        before.combine(after)
    }

    /// `side` の側の受渡だけを残す。`Before` は `from` より前、`After` は `from` 以後の受渡を残す。
    ///
    /// `exercises` を渡すと、`from` より前の判断で `from` 以後の受渡に関わるものを記録の選択肢に畳む。
    fn retain(self, from: Instant, side: Side, exercises: Option<&Exercises<'_>>) -> Result<Self> {
        // 各ノードの下に `from` 以後の受渡があるか。arena は子が親より前に並ぶので、前から一度で求まる。
        let reachable = reachable_contracts(&self.contracts, self.root);
        let mut has_late_pay = vec![false; reachable.len()];
        for (ix, (node, needed)) in self.contracts.iter().zip(reachable).enumerate() {
            if needed {
                let mut late = matches!(node, ContractNode::Pay { at, .. } if *at >= from);
                node.visit_children(|child| late |= has_late_pay[child as usize]);
                has_late_pay[ix] = late;
            }
        }

        // 先に根から選ばれた部分をたどる。記録で退けた選択肢では内側の判断は不要。
        // `Before` では `from` 以後の判断の下をたどらない（因果性により、その下の受渡は `from` 以後）。
        // `After` では `from` 以後の受渡を持たない部分の下をたどらない。
        let mut retained = vec![false; has_late_pay.len()];
        let mut choices = HashMap::new();
        let mut pending = vec![self.root];
        while let Some(ix) = pending.pop() {
            if std::mem::replace(&mut retained[ix as usize], true) {
                continue;
            }
            let node = &self.contracts[ix as usize];
            match node {
                ContractNode::Cond { at, .. } | ContractNode::Choose { at, .. }
                    if side == Side::Before && *at >= from => {}
                ContractNode::Choose {
                    label, at, options, ..
                } if *at < from && has_late_pay[ix as usize] && exercises.is_some() => {
                    // 境界をまたぐ判断は記録の選択肢だけを残す。番号の範囲は observe で確かめてある。
                    let option = exercises
                        .and_then(|records| records.option(label, *at))
                        .ok_or_else(|| missing_record(from, label, *at))?;
                    let child = options[option];
                    choices.insert(ix, child);
                    pending.push(child);
                }
                _ if side == Side::After && !has_late_pay[ix as usize] => {}
                _ => node.visit_children(|child| pending.push(child)),
            }
        }

        // 添字 0 は空の `All`（何も受け渡さない契約）にする。取り除いたノードは 0 に写し、
        // `Choose` では 0 の選択肢として残して、行使記録が指す選択肢の番号を変えない。
        let mut contracts = vec![ContractNode::All(Vec::new())];
        let mut mapped = vec![0; retained.len()];
        for (ix, (node, needed)) in (0..).zip(self.contracts.into_iter().zip(retained)) {
            if !needed {
                continue;
            }
            let child = |old: Ix| mapped[old as usize];
            let node = match node {
                ContractNode::Pay { at, .. } if (at < from) == (side == Side::Before) => node,
                ContractNode::Pay { .. } => continue,
                ContractNode::Cond { at, .. } | ContractNode::Choose { at, .. }
                    if side == Side::Before && at >= from =>
                {
                    continue;
                }
                ContractNode::All(children) => {
                    let children: Vec<_> = children
                        .into_iter()
                        .map(child)
                        .filter(|&ix| ix != 0)
                        .collect();
                    if children.is_empty() {
                        continue;
                    }
                    ContractNode::All(children)
                }
                ContractNode::Cond { cond, at, yes, no } => {
                    let (yes, no) = (child(yes), child(no));
                    if yes == 0 && no == 0 {
                        continue;
                    }
                    ContractNode::Cond { cond, at, yes, no }
                }
                ContractNode::Choose {
                    by,
                    label,
                    at,
                    options,
                } => {
                    // 記録で確定した判断は、選んだ選択肢の写し先をそのまま使う。
                    if let Some(&selected) = choices.get(&ix) {
                        mapped[ix as usize] = child(selected);
                        continue;
                    }
                    let options: Vec<_> = options.into_iter().map(child).collect();
                    if options.iter().all(|&ix| ix == 0) {
                        continue;
                    }
                    ContractNode::Choose {
                        by,
                        label,
                        at,
                        options,
                    }
                }
            };
            mapped[ix as usize] = Ix::try_from(contracts.len()).map_err(|_| node_limit())?;
            contracts.push(node);
        }
        Self {
            contracts,
            exprs: self.exprs,
            root: mapped[self.root as usize],
        }
        .compact()
    }

    /// 二つの契約を `All` で並べる。`other` の式と契約の添字は `self` の後ろへずらす。
    pub(crate) fn combine(mut self, other: Self) -> Result<Self> {
        let other = other.compact()?;
        self = self.compact()?;
        if self.contracts.len() + other.contracts.len() + 1 + self.exprs.len() + other.exprs.len()
            > MAX_NODES
        {
            return Err(node_limit());
        }
        let expr_offset = Ix::try_from(self.exprs.len()).map_err(|_| node_limit())?;
        let offset = Ix::try_from(self.contracts.len()).map_err(|_| node_limit())?;
        for mut node in other.exprs {
            node.try_map_children(|child| Ok(child + expr_offset))?;
            self.exprs.push(node);
        }
        self.contracts
            .extend(other.contracts.into_iter().map(|mut node| {
                match &mut node {
                    ContractNode::Pay { amount, .. } => *amount += expr_offset,
                    ContractNode::Cond { cond, yes, no, .. } => {
                        *cond += expr_offset;
                        *yes += offset;
                        *no += offset;
                    }
                    ContractNode::Choose { options, .. } | ContractNode::All(options) => {
                        for ix in options {
                            *ix += offset;
                        }
                    }
                }
                node
            }));
        let root = Ix::try_from(self.contracts.len()).map_err(|_| node_limit())?;
        self.contracts
            .push(ContractNode::All(vec![self.root, other.root + offset]));
        self.root = root;
        Ok(self)
    }
}

/// 変更の境界で残す側。
#[derive(Clone, Copy, PartialEq, Eq)]
enum Side {
    Before,
    After,
}
`````

`````
crates/rhodium-products/src/ir/analysis.rs
//! 検証済み IR の受渡、観測と、リニア方式への適用可否。

use std::collections::HashSet;

use rhodium_core::{
    currency::Asset,
    date::DateInterval,
    day_count::DayCount,
    fixing::{ObservationPeriod, ObservationPoint},
    fx::CcyPair,
    id::IndexId,
    observation::ObservationRule,
    rounding::Rounding,
    time::Instant,
};

use super::{
    FixingKind, Ix, Party, Product, ResolvedCompounding,
    nodes::{Binary, CompoundingDates, CompoundingRows, ContractNode, ExprNode, Nary, Unary},
    reachable_contracts,
};

/// 根から到達する一つの受渡。
#[derive(Debug)]
pub struct Cashflow<'a> {
    /// 受け渡す資産。
    pub asset: &'a Asset,
    /// 受渡の瞬間。
    pub at: Instant,
    /// 符号を含む受渡量の式。
    pub amount: ExprRef<'a>,
    /// 外側から内側へ並べた、この受渡を囲む分岐。
    pub branches: Vec<Branch<'a>>,
}

/// 受渡が属する条件の枝または行使の選択肢。
#[derive(Debug, Clone, Copy)]
pub enum Branch<'a> {
    /// 条件による分岐。
    Condition {
        /// 判定の瞬間。
        at: Instant,
        /// 判断する条件式。
        condition: ExprRef<'a>,
        /// 選ぶ条件の値。
        when: bool,
    },
    /// 行使する当事者による選択。
    Choice {
        /// 判断する当事者。
        by: Party,
        /// 選択の識別子。
        label: &'a str,
        /// 判断の締切。
        at: Instant,
        /// 選択肢のゼロ始まりの位置。
        option: usize,
    },
}

/// 受渡を契約の順で列挙する。共有した支払ノードも出現ごとに数える。
#[must_use]
pub fn cashflows(product: &Product) -> Vec<Cashflow<'_>> {
    enum Work<'a> {
        Visit(Ix),
        Enter(Branch<'a>, Ix),
        Leave,
    }

    let contract = &product.contract;
    let mut work = vec![Work::Visit(contract.root)];
    let mut branches = Vec::new();
    let mut result = Vec::new();
    while let Some(step) = work.pop() {
        let pos = match step {
            Work::Visit(pos) => pos,
            Work::Enter(branch, pos) => {
                branches.push(branch);
                work.push(Work::Leave);
                pos
            }
            Work::Leave => {
                branches.pop();
                continue;
            }
        };
        match &contract.contracts[pos as usize] {
            ContractNode::Pay { asset, amount, at } => result.push(Cashflow {
                asset,
                at: *at,
                amount: ExprRef {
                    product,
                    ix: *amount,
                },
                branches: branches.clone(),
            }),
            ContractNode::All(children) => {
                work.extend(children.iter().rev().copied().map(Work::Visit));
            }
            ContractNode::Cond { cond, at, yes, no } => {
                let condition = ExprRef { product, ix: *cond };
                for (when, child) in [(false, *no), (true, *yes)] {
                    work.push(Work::Enter(
                        Branch::Condition {
                            at: *at,
                            condition,
                            when,
                        },
                        child,
                    ));
                }
            }
            ContractNode::Choose {
                by,
                label,
                at,
                options,
            } => {
                work.extend(options.iter().enumerate().rev().map(|(option, &child)| {
                    Work::Enter(
                        Branch::Choice {
                            by: *by,
                            label,
                            at: *at,
                            option,
                        },
                        child,
                    )
                }));
            }
        }
    }
    result
}

/// 一つの必要な観測。信用事由と回収率は別の観測として扱う。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Observation<'a> {
    /// 指標の fixing。
    Fixing {
        /// 指標。
        index: &'a IndexId,
        /// 値を参照する観測点。
        point: ObservationPoint,
    },
    /// 信用事由の有無。
    Defaulted {
        /// 信用事由の指標。
        index: &'a IndexId,
        /// 信用事由を判定する窓。
        window: DateInterval,
    },
    /// 信用事由の回収率。
    Recovery {
        /// 回収率の指標。
        index: &'a IndexId,
        /// 信用事由を判定する窓。
        window: DateInterval,
    },
}

/// 到達しうる式が読む観測を、同じ観測を一度ずつ含む一覧にする。
#[must_use]
pub fn observations(product: &Product) -> Vec<Observation<'_>> {
    let contract = &product.contract;
    let reached = reachable_contracts(&contract.contracts, contract.root);
    let mut seen = HashSet::new();
    let mut result = Vec::new();
    let mut add = |observation| {
        if seen.insert(observation) {
            result.push(observation);
        }
    };
    for (node, reached) in contract
        .exprs
        .iter()
        .zip(contract.reachable_exprs(&reached))
    {
        if !reached {
            continue;
        }
        match node {
            ExprNode::Fixing { index, point, .. } => add(Observation::Fixing {
                index,
                point: *point,
            }),
            ExprNode::Compounded { index, rows, .. } => {
                for position in rows.positions() {
                    let ((point, _), _) = rows.source.entries[position];
                    add(Observation::Fixing { index, point });
                }
            }
            ExprNode::Average { index, rows, .. } => {
                for &(point, _) in rows {
                    add(Observation::Fixing { index, point });
                }
            }
            ExprNode::Defaulted { index, window } => add(Observation::Defaulted {
                index,
                window: *window,
            }),
            ExprNode::Recovery { index, window } => add(Observation::Recovery {
                index,
                window: *window,
            }),
            _ => {}
        }
    }
    result
}

/// 検証済みの商品に属する式への読み取り専用の参照。
#[derive(Debug, Clone, Copy)]
pub struct ExprRef<'a> {
    product: &'a Product,
    ix: Ix,
}

impl<'a> ExprRef<'a> {
    /// 同じ商品内の式を識別する番号。子の番号は親より小さい。
    #[must_use]
    pub const fn index(self) -> usize {
        self.ix as usize
    }

    /// 式の演算と、商品内の子への参照を返す。
    #[must_use]
    pub fn expression(self) -> Expression<'a> {
        match self.node() {
            ExprNode::Const(value) => Expression::Constant(value.value()),
            ExprNode::Unary(Unary::Neg, _) => Expression::Negate,
            ExprNode::Unary(Unary::Round(rounding), _) => Expression::Round(*rounding),
            ExprNode::Binary(Binary::Sub, _) => Expression::Subtract,
            ExprNode::Binary(Binary::Div, _) => Expression::Divide,
            ExprNode::Nary(Nary::Add, _) => Expression::Add,
            ExprNode::Nary(Nary::Mul, _) => Expression::Multiply,
            ExprNode::Fixing {
                index,
                point,
                fx_reset_at,
            } => {
                let reference = &self.product.references.fixings[index][point];
                match reference.kind {
                    FixingKind::TermRate => Expression::TermFixing(TermFixing {
                        index,
                        point: *point,
                        period: reference.period,
                        publication: reference.publication,
                        forecast_from: self.product.forecast_from(index, *point),
                    }),
                    FixingKind::Fx(pair) => Expression::FxFixing(FxFixing {
                        index,
                        point: *point,
                        pair,
                        observed_at: reference.period.start(),
                        value_at: *fx_reset_at,
                    }),
                    FixingKind::RateFuture => Expression::RateFutureFixing(RateFutureFixing {
                        index,
                        point: *point,
                        observed_at: reference.period.start(),
                    }),
                    FixingKind::Other => Expression::Other,
                }
            }
            ExprNode::Compounded {
                index,
                rows,
                dates,
                rule,
                daily_spread,
                fixed_excess,
            } => Expression::Compounded(Compounded {
                resolved: &self.product.references.compounded[index][&rows.source],
                index,
                rows,
                dates,
                rule,
                daily_spread: daily_spread.value(),
                fixed_excess: fixed_excess.value(),
            }),
            _ => Expression::Other,
        }
    }

    /// 子を演算の引数の順で返す。
    pub fn children(self) -> impl DoubleEndedIterator<Item = Self> + 'a {
        self.node().children().iter().map(move |&ix| Self {
            product: self.product,
            ix,
        })
    }

    /// 分母が 1 + τL で、分子が同じ L と固定率の差の定数倍なら、L と τ を返す。
    #[must_use]
    pub fn fra_settlement(self) -> Option<(TermFixing<'a>, f64)> {
        let [numerator, denominator] = self.binary_children(Expression::Divide)?;
        let factors = |expr: Self| {
            let [left, right] = expr.binary_children(Expression::Multiply)?;
            match (left.expression(), right.expression()) {
                (Expression::Constant(coefficient), _) => Some((coefficient, right)),
                (_, Expression::Constant(coefficient)) => Some((coefficient, left)),
                _ => None,
            }
        };
        let (_, difference) = factors(numerator)?;
        let [observed, strike] = difference.binary_children(Expression::Subtract)?;
        if !matches!(strike.expression(), Expression::Constant(_)) {
            return None;
        }
        let Expression::TermFixing(fixing) = observed.expression() else {
            return None;
        };
        let [left, right] = denominator.binary_children(Expression::Add)?;
        let interest = match (left.expression(), right.expression()) {
            (Expression::Constant(one), _) if one.to_bits() == 1.0_f64.to_bits() => right,
            (_, Expression::Constant(one)) if one.to_bits() == 1.0_f64.to_bits() => left,
            _ => return None,
        };
        let (tau, rate) = factors(interest)?;
        (tau > 0.0 && rate.index() == observed.index()).then_some((fixing, tau))
    }

    fn binary_children(self, operator: Expression<'_>) -> Option<[Self; 2]> {
        if std::mem::discriminant(&self.expression()) != std::mem::discriminant(&operator) {
            return None;
        }
        let mut children = self.children();
        let left = children.next()?;
        let right = children.next()?;
        children.next().is_none().then_some([left, right])
    }

    fn node(self) -> &'a ExprNode {
        &self.product.contract.exprs[self.index()]
    }
}

/// リニア評価で照合する式の種類。
#[derive(Debug, Clone, Copy)]
pub enum Expression<'a> {
    /// 定数。
    Constant(f64),
    /// 符号の反転。
    Negate,
    /// 契約の丸め。
    Round(Rounding),
    /// 和。
    Add,
    /// 左から右を引く。
    Subtract,
    /// 積。
    Multiply,
    /// 左を右で割る。
    Divide,
    /// 予想するターム金利。
    TermFixing(TermFixing<'a>),
    /// 将来の FX fixing。
    FxFixing(FxFixing<'a>),
    /// 将来の金利先物の清算価格。
    RateFutureFixing(RateFutureFixing<'a>),
    /// 翌日物の複利。
    Compounded(Compounded<'a>),
    /// リニア方式が扱う額の演算以外。
    Other,
}

/// 簡約後も予想が必要な FX の観測。
#[derive(Debug, Clone, Copy)]
pub struct FxFixing<'a> {
    /// 観測する指標。
    pub index: &'a IndexId,
    /// 指標の観測点。
    pub point: ObservationPoint,
    /// スポットを観測する瞬間。
    pub observed_at: Instant,
    /// 元本リセットで、この fixing で換算した新しい元本を交換する瞬間（脚の通貨の決済時刻）。
    ///
    /// fixing の予想を読む受渡ではない。同じ指標・観測日の公表値は一つなので、予想は通常の
    /// fixing と同じくペアのスポットの受渡で読み、この瞬間とのずれは支払いの遅れとして扱う。
    pub value_at: Option<Instant>,
    /// 観測のペアの向き。
    pub pair: CcyPair,
}

/// 簡約後も予想が必要な日次の清算価格。
#[derive(Debug, Clone, Copy)]
pub struct RateFutureFixing<'a> {
    /// 限月固有の指標。
    pub index: &'a IndexId,
    /// 指標の観測点。
    pub point: ObservationPoint,
    /// 清算価格を決める瞬間。
    pub observed_at: Instant,
}

/// 簡約後も予想が必要なターム金利の観測。
#[derive(Debug, Clone, Copy)]
pub struct TermFixing<'a> {
    /// 元の指標。
    pub index: &'a IndexId,
    /// 指標の観測日。
    pub point: ObservationPoint,
    /// 指標の預金期間。
    pub period: ObservationPeriod,
    /// 公表の瞬間。
    pub publication: Instant,
    /// 期間の途中から予想する場合の開始の瞬間。
    pub forecast_from: Option<Instant>,
}

/// 簡約済み IR に残る翌日物の複利。
#[derive(Debug, Clone, Copy)]
pub struct Compounded<'a> {
    resolved: &'a ResolvedCompounding,
    rows: &'a CompoundingRows,
    dates: &'a CompoundingDates,
    /// 翌日物の指標。
    pub index: &'a IndexId,
    /// 複利と観測の規則。
    pub rule: &'a ObservationRule,
    /// 日次因子に含めるスプレッド。
    pub daily_spread: f64,
    /// 簡約が畳んだ確定因子の積から 1 を引いた値 e。成長率は (1 + e)·Π（残りの因子）。
    ///
    /// 超過成長 Π − 1 は、残りの因子の超過 f と合わせて e + f + e·f で求めると桁を落とさない。
    pub fixed_excess: f64,
}

impl<'a> Compounded<'a> {
    /// 契約の観測の重みを求めた日数基準。
    #[must_use]
    pub const fn basis(self) -> &'a DayCount {
        &self.dates.basis
    }

    /// 確定済みの因子が元の列の先頭に連続していたか。
    #[must_use]
    pub const fn has_fixed_prefix(self) -> bool {
        self.dates.fixed_prefix
    }

    /// 簡約前の最後の観測の期間。
    #[must_use]
    pub fn last_period(self) -> Option<ObservationPeriod> {
        self.resolved
            .observations
            .last()
            .map(|fixing| fixing.period)
    }

    /// まだ値の必要な観測を、元の重みとともに返す。
    pub fn rows(self) -> impl ExactSizeIterator<Item = RfrObservation> + 'a {
        self.rows.positions().map(move |position| {
            let ((point, weight), coverage) = self.rows.source.entries[position];
            let fixing = &self.resolved.observations[position];
            RfrObservation {
                point,
                weight: weight.value(),
                coverage,
                period: fixing.period,
                publication: fixing.publication,
            }
        })
    }
}

/// RFR 複利に残る一つの観測。
#[derive(Debug, Clone, Copy)]
pub struct RfrObservation {
    /// 指標の観測点。
    pub point: ObservationPoint,
    /// 契約で用いる重み。
    pub weight: f64,
    /// 契約の重みが覆う現地日の区間。
    pub coverage: DateInterval,
    /// 指標が表す金利の期間。
    pub period: ObservationPeriod,
    /// 指標が確定する瞬間。
    pub publication: Instant,
}

/// リニア評価への適用可否。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Linearity {
    /// リニア評価で扱える式。
    Linear,
    /// リニア評価を適用できない理由。
    NonLinear(&'static str),
}

/// 簡約後の契約にリニア方式を適用できるかを判定する。
///
/// 受渡額が定数、ターム金利と翌日物の複利、FX の fixing（その逆数を含む）、金利先物の清算価格の
/// 線形結合で、予想の積は FX と金利の一つずつ（MtM XCCY の元本リセット）、予想での割り算は FX の
/// 逆数と FRA の決済だけのときにリニアとする。
#[must_use]
pub fn linearity(product: &Product) -> Linearity {
    // 受渡額が含む予想の種類を追い、FX と金利の一つずつの積だけを近似の対象にする。
    const RATE: u8 = 1;
    const FX: u8 = 2;
    const CLEARING: u8 = 4;
    let contract = &product.contract;
    let reached = reachable_contracts(&contract.contracts, contract.root);
    for (node, &reached) in contract.contracts.iter().zip(&reached) {
        if !reached {
            continue;
        }
        match node {
            ContractNode::Choose { .. } => {
                return Linearity::NonLinear("unresolved exercise choice");
            }
            ContractNode::Cond { .. } => {
                return Linearity::NonLinear("unresolved cashflow condition");
            }
            ContractNode::Pay { asset, .. } if !matches!(asset, Asset::Cash(_)) => {
                return Linearity::NonLinear("non-currency delivery");
            }
            _ => {}
        }
    }

    let reached = contract.reachable_exprs(&reached);
    let mut projected = vec![0; reached.len()];
    for (ix, reached) in (0..).zip(reached) {
        if !reached {
            continue;
        }
        let pos = ix as usize;
        let reference = ExprRef { product, ix };
        let mut dynamics = reference.children().map(|child| projected[child.index()]);
        projected[pos] = match reference.expression() {
            Expression::Constant(_) => 0,
            Expression::Compounded(_) | Expression::TermFixing(_) => RATE,
            Expression::FxFixing(_) => FX,
            Expression::RateFutureFixing(_) => CLEARING,
            Expression::Negate | Expression::Round(_) | Expression::Add | Expression::Subtract => {
                dynamics.fold(0, |kinds, kind| kinds | kind)
            }
            Expression::Multiply => {
                let mut dynamic = dynamics.filter(|&item| item != 0);
                match (dynamic.next(), dynamic.next(), dynamic.next()) {
                    (None, None, None) => 0,
                    (Some(kind), None, None) => kind,
                    (Some(FX), Some(RATE), None) | (Some(RATE), Some(FX), None) => FX | RATE,
                    _ => return Linearity::NonLinear("multiplication of two projected quantities"),
                }
            }
            Expression::Divide => {
                let denominator = dynamics.next_back().unwrap_or(0);
                let numerator = dynamics.next().unwrap_or(0);
                match denominator {
                    0 => numerator,
                    FX if numerator == 0
                        && reference.children().next_back().is_some_and(|child| {
                            matches!(child.expression(), Expression::FxFixing(_))
                        }) =>
                    {
                        FX
                    }
                    RATE if reference.fra_settlement().is_some() => RATE,
                    _ => return Linearity::NonLinear("division by a projected quantity"),
                }
            }
            Expression::Other => return Linearity::NonLinear("unsupported cashflow expression"),
        };
    }
    Linearity::Linear
}
`````

`````
crates/rhodium-products/src/ir/builder.rs
use std::marker::PhantomData;

use rhodium_core::{
    Error, Result,
    brand::Brand,
    currency::Asset,
    date::DateInterval,
    fixing::ObservationPoint,
    id::IndexId,
    observation::{Aggregation, ObservationRule, ObservationSchedule},
    quantity::Spread,
    rounding::Rounding,
    time::Instant,
};
use scioxide_core::checked_sum;

use super::{
    Bool, Contract, ContractIx, ExprIx, Party,
    nodes::{
        Arena, Binary, ContractNode, ExprNode, MAX_NODES, Nary, Number, Row, Unary, node_limit,
    },
};

/// 契約と式を子から親の順に組み立て、同じノードを共有する。
#[derive(Debug)]
pub struct Builder<'b> {
    brand: Brand<'b>,
    contracts: Arena<ContractNode>,
    exprs: Arena<ExprNode>,
    error: Option<Error>,
}

macro_rules! unary {
    ($name:ident, $op:ident, $doc:literal) => {
        #[doc = $doc]
        pub fn $name(&mut self, x: ExprIx<'b>) -> ExprIx<'b> {
            self.expr(ExprNode::Unary(Unary::$op, [x.ix]))
        }
    };
}

macro_rules! binary {
    ($name:ident, $op:ident, $out:ty, $doc:literal) => {
        #[doc = $doc]
        pub fn $name(&mut self, left: ExprIx<'b>, right: ExprIx<'b>) -> ExprIx<'b, $out> {
            self.expr(ExprNode::Binary(Binary::$op, [left.ix, right.ix]))
        }
    };
}

macro_rules! nary {
    ($name:ident, $op:ident, $sort:ty, $doc:literal) => {
        #[doc = $doc]
        pub fn $name(
            &mut self,
            terms: impl IntoIterator<Item = ExprIx<'b, $sort>>,
        ) -> ExprIx<'b, $sort> {
            self.expr(ExprNode::Nary(
                Nary::$op,
                terms.into_iter().map(|term| term.ix).collect(),
            ))
        }
    };
}

impl<'b> Builder<'b> {
    pub(super) fn new(brand: Brand<'b>) -> Self {
        Self {
            brand,
            contracts: Arena::default(),
            exprs: Arena::default(),
            error: None,
        }
    }

    pub(super) fn finish(self, root: ContractIx<'b>) -> Result<Contract> {
        if let Some(error) = self.error {
            return Err(error);
        }
        Ok(Contract {
            contracts: self.contracts.nodes,
            exprs: self.exprs.nodes,
            root: root.ix,
        })
    }

    fn check_budget(&mut self) {
        if self.contracts.nodes.len() + self.exprs.nodes.len() > MAX_NODES {
            self.error = Some(node_limit());
        }
    }

    fn expr<T>(&mut self, node: ExprNode) -> ExprIx<'b, T> {
        let ix = if self.error.is_some() {
            0
        } else {
            match self.exprs.insert(node) {
                Ok(ix) => ix,
                Err(error) => {
                    self.error = Some(error);
                    0
                }
            }
        };
        self.check_budget();
        ExprIx {
            ix,
            _brand: self.brand,
            sort: PhantomData,
        }
    }

    fn contract(&mut self, node: ContractNode) -> ContractIx<'b> {
        let ix = if self.error.is_some() {
            0
        } else {
            match self.contracts.insert(node) {
                Ok(ix) => ix,
                Err(error) => {
                    self.error = Some(error);
                    0
                }
            }
        };
        self.check_budget();
        ContractIx {
            ix,
            _brand: self.brand,
        }
    }

    /// 数の定数。非有限値を含む商品は検査時に拒否する。
    pub fn constant(&mut self, value: f64) -> ExprIx<'b> {
        self.expr(ExprNode::Const(Number::new(value)))
    }

    /// 真偽値の定数。
    pub fn boolean(&mut self, value: bool) -> ExprIx<'b, Bool> {
        self.expr(ExprNode::Bool(value))
    }

    /// 一つの観測点の fixing。
    pub fn fixing(&mut self, index: &IndexId, point: ObservationPoint) -> ExprIx<'b> {
        self.expr(ExprNode::Fixing {
            index: index.clone(),
            point,
            fx_reset_at: None,
        })
    }

    /// 元本換算に使う FX fixing と、新しい元本の受渡の瞬間を指定する。
    ///
    /// `value_at` は元本を交換する瞬間で、fixing の予想はペアのスポットの受渡で読む。
    pub fn fx_reset(
        &mut self,
        index: &IndexId,
        point: ObservationPoint,
        value_at: Instant,
    ) -> ExprIx<'b> {
        self.expr(ExprNode::Fixing {
            index: index.clone(),
            point,
            fx_reset_at: Some(value_at),
        })
    }

    /// 複利の観測列から、年率化する前の成長率を作る。
    ///
    /// # Errors
    ///
    /// 規則が複利用でないか、観測列を生成した規則と異なるときは `InvalidInput`。
    pub fn compounded(
        &mut self,
        index: &IndexId,
        observations: &ObservationSchedule,
        rule: &ObservationRule,
        spread: Spread,
    ) -> Result<ExprIx<'b>> {
        let Aggregation::Compound { .. } = rule.aggregation() else {
            return Err(Error::InvalidInput(
                "compounded requires a compounding rule".into(),
            ));
        };
        if observations.rule() != rule {
            return Err(Error::InvalidInput(
                "compounding rule does not match observation schedule".into(),
            ));
        }
        Ok(self.expr(ExprNode::compounded(index, observations, spread)))
    }

    /// 平均用の観測列の重みを一度だけ正規化し、加重平均を作る。
    ///
    /// # Errors
    ///
    /// 平均の規則が無いか、重みの和が正の有限値にならないときは `InvalidInput`。
    pub fn average(
        &mut self,
        index: &IndexId,
        observations: &ObservationSchedule,
    ) -> Result<ExprIx<'b>> {
        let rule = observations.rule();
        if !matches!(rule.aggregation(), Aggregation::Average(_)) {
            return Err(Error::InvalidInput(
                "average requires an averaging rule".into(),
            ));
        }
        let weights: Vec<_> = observations.fixings().iter().map(|&(_, w)| w).collect();
        let total = checked_sum(&weights)
            .ok()
            .filter(|total| *total > 0.0)
            .ok_or_else(|| Error::InvalidInput("average requires positive total weight".into()))?;
        let rows: Vec<Row> = observations
            .fixings()
            .iter()
            .map(|&(point, weight)| (point, Number::new(weight / total)))
            .collect();
        Ok(self.expr(ExprNode::Average {
            index: index.clone(),
            rows,
            daily_rounding: rule.rounding().daily,
        }))
    }

    /// 窓 (開始日, 終了日] に信用事由が起きたか。
    pub fn defaulted(&mut self, index: &IndexId, window: DateInterval) -> ExprIx<'b, Bool> {
        self.expr(ExprNode::Defaulted {
            index: index.clone(),
            window,
        })
    }

    /// 窓 (開始日, 終了日] に起きた信用事由の回収率。
    pub fn recovery(&mut self, index: &IndexId, window: DateInterval) -> ExprIx<'b> {
        self.expr(ExprNode::Recovery {
            index: index.clone(),
            window,
        })
    }

    unary!(neg, Neg, "符号の反転。");
    unary!(abs, Abs, "絶対値。");
    unary!(exp, Exp, "指数関数。");
    unary!(ln, Ln, "自然対数。");
    unary!(sqrt, Sqrt, "平方根。");

    /// 契約が定める丸め。
    pub fn round(&mut self, x: ExprIx<'b>, rounding: Rounding) -> ExprIx<'b> {
        self.expr(ExprNode::Unary(Unary::Round(rounding), [x.ix]))
    }

    /// 論理否定。
    pub fn not(&mut self, x: ExprIx<'b, Bool>) -> ExprIx<'b, Bool> {
        self.expr(ExprNode::Unary(Unary::Not, [x.ix]))
    }

    binary!(sub, Sub, super::Num, "二項の差。");
    binary!(div, Div, super::Num, "二項の商。");
    binary!(pow, Pow, super::Num, "底を指数で累乗する。");
    binary!(lt, Lt, Bool, "左辺が右辺より小さいか。");
    binary!(le, Le, Bool, "左辺が右辺以下か。");
    binary!(gt, Gt, Bool, "左辺が右辺より大きいか。");
    binary!(ge, Ge, Bool, "左辺が右辺以上か。");

    nary!(add, Add, super::Num, "項の和。空の和は 0。");
    nary!(mul, Mul, super::Num, "項の積。空の積は 1。");
    nary!(min, Min, super::Num, "項の最小値。少なくとも一項を要する。");
    nary!(max, Max, super::Num, "項の最大値。少なくとも一項を要する。");
    nary!(and, And, Bool, "順に短絡する論理積。空なら真。");
    nary!(or, Or, Bool, "順に短絡する論理和。空なら偽。");

    /// 真偽値に従い、選んだ数の式だけを評価する。
    pub fn select(
        &mut self,
        cond: ExprIx<'b, Bool>,
        yes: ExprIx<'b>,
        no: ExprIx<'b>,
    ) -> ExprIx<'b> {
        self.expr(ExprNode::Select([cond.ix, yes.ix, no.ix]))
    }

    /// 受渡の瞬間に指定した量の資産を受け取る。負の量なら支払う。
    pub fn pay(&mut self, asset: Asset, amount: ExprIx<'b>, at: Instant) -> ContractIx<'b> {
        self.contract(ContractNode::Pay {
            asset,
            amount: amount.ix,
            at,
        })
    }

    /// すべての子契約を受け取る。
    pub fn all(&mut self, children: impl IntoIterator<Item = ContractIx<'b>>) -> ContractIx<'b> {
        self.contract(ContractNode::All(
            children.into_iter().map(|child| child.ix).collect(),
        ))
    }

    /// 受渡のない契約。
    pub fn zero(&mut self) -> ContractIx<'b> {
        self.all([])
    }

    /// 判定の瞬間の条件に従い、一つの枝を受け取る。
    pub fn cond(
        &mut self,
        cond: ExprIx<'b, Bool>,
        at: Instant,
        yes: ContractIx<'b>,
        no: ContractIx<'b>,
    ) -> ContractIx<'b> {
        self.contract(ContractNode::Cond {
            cond: cond.ix,
            at,
            yes: yes.ix,
            no: no.ix,
        })
    }

    /// 判断の瞬間に、指定した当事者が選択肢の一つを選ぶ。
    pub fn choose(
        &mut self,
        by: Party,
        label: &str,
        at: Instant,
        options: impl IntoIterator<Item = ContractIx<'b>>,
    ) -> ContractIx<'b> {
        self.contract(ContractNode::Choose {
            by,
            label: label.to_owned(),
            at,
            options: options.into_iter().map(|option| option.ix).collect(),
        })
    }

    /// 支払額の符号と選択の当事者を入れ替える。
    pub fn give(&mut self, contract: ContractIx<'b>) -> ContractIx<'b> {
        if self.error.is_some() {
            return contract;
        }
        let needed = super::reachable_contracts(&self.contracts.nodes, contract.ix);
        let mut reversed = vec![0; needed.len()];
        for (pos, needed) in needed.into_iter().enumerate() {
            if !needed {
                continue;
            }
            let node = match self.contracts.nodes[pos].clone() {
                ContractNode::Pay { asset, amount, at } => {
                    let amount: ExprIx<'b> = self.expr(ExprNode::Unary(Unary::Neg, [amount]));
                    ContractNode::Pay {
                        asset,
                        amount: amount.ix,
                        at,
                    }
                }
                ContractNode::All(children) => ContractNode::All(
                    children
                        .into_iter()
                        .map(|ix| reversed[ix as usize])
                        .collect(),
                ),
                ContractNode::Cond { cond, at, yes, no } => ContractNode::Cond {
                    cond,
                    at,
                    yes: reversed[yes as usize],
                    no: reversed[no as usize],
                },
                ContractNode::Choose {
                    by,
                    label,
                    at,
                    options,
                } => ContractNode::Choose {
                    by: by.other(),
                    label,
                    at,
                    options: options
                        .into_iter()
                        .map(|ix| reversed[ix as usize])
                        .collect(),
                },
            };
            reversed[pos] = self.contract(node).ix;
        }
        ContractIx {
            ix: reversed[contract.ix as usize],
            _brand: self.brand,
        }
    }
}
`````

`````
crates/rhodium-products/src/ir/fallback.rs
use rhodium_core::{
    Error, Result, Warning,
    date::Date,
    fixing::ObservationPoint,
    id::IndexId,
    observation::{Aggregation, SpreadTreatment, fallback_observation_schedule},
    refdata::{CalendarCache, FallbackReplacement, FallbackTrigger, IndexStatus, RefData},
};

use super::{
    Contract,
    nodes::{Arena, Binary, ContractNode, ExprNode, Ix, Nary, Number, Unary},
    reachable_contracts,
};
use crate::trade::FallbackTerms;

pub(super) fn lower(
    mut contract: Contract,
    refdata: &RefData,
    terms: &FallbackTerms,
    trade_date: Date,
) -> Result<Contract> {
    if terms.by_index.is_empty()
        && terms.definitions.is_none()
        && terms.protocols.is_empty()
        && terms.bilateral.is_empty()
    {
        return Ok(contract);
    }
    let select = |node: &ExprNode| {
        let ExprNode::Fixing {
            index,
            point: ObservationPoint::Date(date),
            fx_reset_at: None,
        } = node
        else {
            return Ok(None);
        };
        let Some((trigger, replacement)) = terms.select(index, trade_date, refdata) else {
            return Ok(None);
        };
        let applies = match refdata.index(index)?.status {
            IndexStatus::Active => false,
            IndexStatus::Ceased { since } => *date >= since,
            IndexStatus::NonRepresentative { since } => {
                *date >= since && trigger == FallbackTrigger::CessationOrNonRepresentative
            }
        };
        Ok(applies.then(|| (index.clone(), *date, replacement)))
    };

    // 通常の商品の式をコピーせず、必要な代替がある場合だけ同じ DAG を書き直す。
    if !contract.exprs.iter().any(|node| {
        matches!(node, ExprNode::Fixing { index, point: ObservationPoint::Date(_), fx_reset_at: None }
            if terms.select(index, trade_date, refdata).is_some())
    }) {
        return Ok(contract);
    }
    contract = contract.compact()?;
    let mut replacements = contract
        .exprs
        .iter()
        .map(select)
        .collect::<Result<Vec<_>>>()?;
    if replacements.iter().all(Option::is_none) {
        return Ok(contract);
    }
    let mut calendars = refdata.calendar_cache();
    let reached = reachable_contracts(&contract.contracts, contract.root);
    let mut exprs = Arena::default();
    let mut ids = Vec::with_capacity(contract.exprs.len());
    for (mut node, replacement) in contract.exprs.into_iter().zip(replacements.drain(..)) {
        let ix = if let Some((index, date, replacement)) = replacement {
            build_rate(
                &mut exprs,
                &index,
                date,
                replacement,
                refdata,
                &mut calendars,
            )?
        } else {
            node.try_map_children(|child| Ok(ids[child as usize]))?;
            exprs.insert(node)?
        };
        ids.push(ix);
    }
    for (node, needed) in contract.contracts.iter_mut().zip(reached) {
        if needed {
            match node {
                ContractNode::Pay { amount, .. } => *amount = ids[*amount as usize],
                ContractNode::Cond { cond, .. } => *cond = ids[*cond as usize],
                ContractNode::All(_) | ContractNode::Choose { .. } => {}
            }
        }
    }
    contract.exprs = exprs.nodes;
    Ok(contract)
}

fn build_rate(
    exprs: &mut Arena<ExprNode>,
    index: &IndexId,
    date: Date,
    replacement: &FallbackReplacement,
    refdata: &RefData,
    calendars: &mut CalendarCache<'_>,
) -> Result<Ix> {
    let (period, mut warnings) = refdata.term_rate_dates(index, date)?;
    let (observations, notices) = fallback_observation_schedule(period, replacement, calendars)?;
    warnings.extend(notices);
    let tau = observations.annualization().value();
    if tau <= 0.0 {
        return Err(Error::InvalidData(format!(
            "fallback for {index} requires positive annualization"
        )));
    }
    let growth = exprs.insert(ExprNode::compounded(
        &replacement.index,
        &observations,
        replacement.spread,
    ))?;
    let one = exprs.insert(ExprNode::Const(Number::new(1.0)))?;
    let excess = exprs.insert(ExprNode::Binary(Binary::Sub, [growth, one]))?;
    let tau = exprs.insert(ExprNode::Const(Number::new(tau)))?;
    let mut rate = exprs.insert(ExprNode::Binary(Binary::Div, [excess, tau]))?;
    if matches!(
        replacement.observation.aggregation(),
        Aggregation::Compound {
            spread: SpreadTreatment::Exclusive,
            ..
        }
    ) && replacement.spread.value() != 0.0
    {
        let spread = exprs.insert(ExprNode::Const(Number::new(replacement.spread.value())))?;
        rate = exprs.insert(ExprNode::Nary(Nary::Add, vec![rate, spread]))?;
    }
    if let Some(rounding) = replacement.observation.rounding().aggregated {
        rate = exprs.insert(ExprNode::Unary(Unary::Round(rounding), [rate]))?;
    }
    warnings.push(Warning::FallbackApplied {
        index: index.clone(),
        replacement: replacement.index.clone(),
    });
    warnings.sort_unstable();
    warnings.dedup();
    exprs.insert(ExprNode::WithWarnings(warnings.into_boxed_slice(), [rate]))
}
`````

`````
crates/rhodium-products/src/ir/nodes.rs
use std::{
    collections::HashMap,
    hash::{BuildHasher, DefaultHasher, Hash, Hasher},
    sync::Arc,
};

use rhodium_core::{
    Error, Result, Warning,
    currency::Asset,
    date::DateInterval,
    day_count::DayCount,
    fixing::ObservationPoint,
    id::IndexId,
    observation::{Aggregation, ObservationRule, ObservationSchedule, SpreadTreatment},
    quantity::Spread,
    rounding::Rounding,
    time::Instant,
};

use super::Party;

pub(super) const MAX_NODES: usize = 1_000_000;
pub(super) type Ix = u32;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(super) enum ContractNode {
    Pay {
        asset: Asset,
        amount: Ix,
        at: Instant,
    },
    All(Vec<Ix>),
    Cond {
        cond: Ix,
        at: Instant,
        yes: Ix,
        no: Ix,
    },
    Choose {
        by: Party,
        label: String,
        at: Instant,
        options: Vec<Ix>,
    },
}

impl ContractNode {
    pub(super) fn children(&self) -> &[Ix] {
        match self {
            Self::Pay { .. } | Self::Cond { .. } => &[],
            Self::All(children)
            | Self::Choose {
                options: children, ..
            } => children,
        }
    }

    pub(super) fn visit_children(&self, mut visit: impl FnMut(Ix)) {
        match self {
            Self::Cond { yes, no, .. } => {
                visit(*yes);
                visit(*no);
            }
            _ => self.children().iter().copied().for_each(visit),
        }
    }
}

/// ビット列で共有を判定し、負のゼロも元の値のまま保持する。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) struct Number(u64);

impl Number {
    pub(super) const fn new(value: f64) -> Self {
        Self(value.to_bits())
    }

    pub(super) const fn value(self) -> f64 {
        f64::from_bits(self.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) enum Unary {
    Neg,
    Abs,
    Exp,
    Ln,
    Sqrt,
    Round(Rounding),
    Not,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) enum Binary {
    Sub,
    Div,
    Pow,
    Lt,
    Le,
    Gt,
    Ge,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) enum Nary {
    Add,
    Mul,
    Min,
    Max,
    And,
    Or,
}

pub(super) type Row = (ObservationPoint, Number);

/// 重みと対象期間を合わせて共有する。arena と参照表は計算済みのハッシュを使う。
#[derive(Debug)]
pub(super) struct RfrRows {
    pub(super) entries: Box<[(Row, DateInterval)]>,
    hash: u64,
}

impl PartialEq for RfrRows {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self, other) || (self.hash == other.hash && self.entries == other.entries)
    }
}

impl Eq for RfrRows {}

impl Hash for RfrRows {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.hash);
    }
}

/// `Compounded` が残す観測の列。行の値は写さず、共有した元の列を位置で参照する。
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(super) struct CompoundingRows {
    /// 簡約前の全行。検証で同じ内容の列を一つの `Arc` にまとめ、`TradeProductResolver` を使えば取引の間でも共有する。
    /// `References::compounded` の解決もこの `Arc` を鍵にして引く。
    pub(super) source: Arc<RfrRows>,
    /// まだ値の必要な行の、`source` での位置（昇順）。`None` は全行が残っていることを表す。
    /// 簡約し直したときも現在の列の添字ではなく元の位置を持つので、解決と対象期間を同じ位置で引ける。
    pending: Option<Arc<[usize]>>,
}

impl CompoundingRows {
    fn new(observations: &ObservationSchedule) -> Self {
        let entries: Box<[_]> = observations
            .fixings()
            .iter()
            .zip(observations.coverage())
            .map(|(&(point, weight), &coverage)| ((point, Number::new(weight)), coverage))
            .collect();
        let mut hasher = DefaultHasher::new();
        entries.hash(&mut hasher);
        Self {
            source: Arc::new(RfrRows {
                entries,
                hash: hasher.finish(),
            }),
            pending: None,
        }
    }

    pub(super) fn len(&self) -> usize {
        self.pending
            .as_ref()
            .map_or(self.source.entries.len(), |pending| pending.len())
    }

    pub(super) fn positions(&self) -> impl ExactSizeIterator<Item = usize> + '_ {
        (0..self.len()).map(|pos| self.pending.as_ref().map_or(pos, |pending| pending[pos]))
    }

    pub(super) fn with_pending(&self, pending: Vec<usize>) -> Self {
        if pending.len() == self.len() {
            return self.clone();
        }
        Self {
            source: Arc::clone(&self.source),
            pending: Some(pending.into()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(super) struct CompoundingDates {
    pub(super) basis: DayCount,
    pub(super) fixed_prefix: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(super) enum ExprNode {
    Const(Number),
    Bool(bool),
    Fixing {
        index: IndexId,
        point: ObservationPoint,
        fx_reset_at: Option<Instant>,
    },
    Compounded {
        index: IndexId,
        rows: CompoundingRows,
        dates: CompoundingDates,
        rule: ObservationRule,
        daily_spread: Number,
        /// 簡約が畳んだ確定因子の積から 1 を引いた値（積の超過）。builder は 0（因子を畳んでいない）を入れる。
        ///
        /// 確定因子の積を丸めて持つと、短い期間の超過成長 Π − 1 の下位の桁を ulp(1) で失うので、
        /// 超過のまま持つ。因子の超過 x を e ← e + x + e·x で畳むので、四則演算だけで環境によらず同じ値になる。
        fixed_excess: Number,
    },
    Average {
        index: IndexId,
        rows: Vec<Row>,
        daily_rounding: Option<Rounding>,
    },
    Defaulted {
        index: IndexId,
        window: DateInterval,
    },
    Recovery {
        index: IndexId,
        window: DateInterval,
    },
    Unary(Unary, [Ix; 1]),
    Binary(Binary, [Ix; 2]),
    Nary(Nary, Vec<Ix>),
    Select([Ix; 3]),
    /// 代替の率を観測したときに加える診断。簡約で中の式だけを残す。
    WithWarnings(Box<[Warning]>, [Ix; 1]),
    /// 簡約が未確定の枝に残した、確定した観測から値を作れなかった集約（理由を持つ）。
    ///
    /// `Compounded` と `Average` の確定部分（丸め、因子、超過や項の和）の失敗は残りの観測の予想によらず
    /// 決まるので、公表値を残さず失敗の理由だけを持ち、値が必要になった評価で報告させる。
    /// 簡約だけが作り、builder と `Product::new` の入力には現れない。
    Undefined(Box<str>),
}

impl ExprNode {
    pub(super) fn compounded(
        index: &IndexId,
        observations: &ObservationSchedule,
        spread: Spread,
    ) -> Self {
        let daily_spread = if matches!(
            observations.rule().aggregation(),
            Aggregation::Compound {
                spread: SpreadTreatment::Inclusive,
                ..
            }
        ) {
            spread.value()
        } else {
            0.0
        };
        Self::Compounded {
            index: index.clone(),
            rows: CompoundingRows::new(observations),
            dates: CompoundingDates {
                basis: observations.annualization().basis().clone(),
                fixed_prefix: true,
            },
            rule: observations.rule().clone(),
            daily_spread: Number::new(daily_spread),
            fixed_excess: Number::new(0.0),
        }
    }

    pub(super) fn children(&self) -> &[Ix] {
        match self {
            Self::Unary(_, xs) | Self::WithWarnings(_, xs) => xs,
            Self::Binary(_, xs) => xs,
            Self::Nary(_, xs) => xs,
            Self::Select(xs) => xs,
            _ => &[],
        }
    }

    pub(super) fn try_map_children(&mut self, mut map: impl FnMut(Ix) -> Result<Ix>) -> Result<()> {
        let children: &mut [Ix] = match self {
            Self::Unary(_, xs) | Self::WithWarnings(_, xs) => xs,
            Self::Binary(_, xs) => xs,
            Self::Nary(_, xs) => xs,
            Self::Select(xs) => xs,
            _ => return Ok(()),
        };
        for child in children {
            *child = map(*child)?;
        }
        Ok(())
    }
}

#[derive(Debug)]
pub(super) struct Arena<N> {
    pub(super) nodes: Vec<N>,
    heads: HashMap<u64, Ix>,
    collisions: Vec<Option<Ix>>,
}

impl<N> Default for Arena<N> {
    fn default() -> Self {
        Self {
            nodes: Vec::new(),
            heads: HashMap::new(),
            collisions: Vec::new(),
        }
    }
}

impl<N: Eq + Hash> Arena<N> {
    pub(super) fn insert(&mut self, node: N) -> Result<Ix> {
        let hash = self.heads.hasher().hash_one(&node);
        let head = self.heads.get(&hash).copied();
        let mut next = head;
        while let Some(ix) = next {
            if self.nodes[ix as usize] == node {
                return Ok(ix);
            }
            next = self.collisions[ix as usize];
        }
        let ix = Ix::try_from(self.nodes.len()).map_err(|_| node_limit())?;
        self.nodes.push(node);
        self.collisions.push(head);
        self.heads.insert(hash, ix);
        Ok(ix)
    }
}

pub(super) fn node_limit() -> Error {
    Error::InvalidInput(format!("contract exceeds the {MAX_NODES} node limit"))
}
`````

`````
crates/rhodium-products/src/ir/reduce.rs
use std::collections::{BTreeSet, HashMap, HashSet};

use rhodium_core::{
    Error, Result, Warning,
    date::DateInterval,
    fixing::{FixingResolution, FixingStore, ObservationPoint},
    id::IndexId,
    observation::Aggregation,
    rounding::Rounding,
    time::Instant,
};
use scioxide_core::{
    checked_float_product, checked_float_product_factors, checked_float_quotient, checked_sum,
};

use super::{
    ByChoice, ByIndex, Contract, Product, ResolvedFixing, SettlementBoundary,
    nodes::{
        Arena, Binary, CompoundingDates, ContractNode, ExprNode, Ix, Nary, Number, Row, Unary,
    },
};
use crate::trade::ExerciseRecord;

mod fold;

/// 評価時点で残る受渡と判断を求め、到達した式を確定値で簡約する。
///
/// 評価時点やスナップショットを変える場合は、簡約前の商品を使う。
///
/// # Errors
///
/// 不正な行使記録は `InvalidData`。必要な fixing または行使記録の欠損は
/// `MissingData`。到達した定数演算の失敗は `Numerical`。
pub fn reduce(
    product: &Product,
    as_of: Instant,
    fixings: &FixingStore,
    exercises: &[ExerciseRecord],
    settlement: SettlementBoundary,
) -> Result<Product> {
    let exercises = effective_exercises(product, as_of, exercises)?;
    let mut reducer = Reducer {
        product,
        as_of,
        fixings,
        exercises,
        settlement,
        contracts: Arena::default(),
        exprs: Arena::default(),
        contract_ids: Copies::new(product.contract.contracts.len()),
        expr_ids: Copies::new(product.contract.exprs.len()),
        observations: HashMap::new(),
        forecasts: HashMap::new(),
        warnings: product.construction_warnings.iter().cloned().collect(),
    };
    let root = reducer.contract()?;
    let contract = Contract {
        contracts: reducer.contracts.nodes,
        exprs: reducer.exprs.nodes,
        root,
    }
    .compact()?;
    Ok(Product {
        contract,
        references: product.references.clone(),
        forecasts: reducer.forecasts,
        construction_warnings: product.construction_warnings.clone(),
        warnings: reducer.warnings.into_iter().collect(),
        trade_exercises: product.trade_exercises.clone(),
    })
}

fn effective_exercises(
    product: &Product,
    as_of: Instant,
    records: &[ExerciseRecord],
) -> Result<ByChoice<usize>> {
    let mut seen = HashSet::new();
    let mut effective: ByChoice<usize> = HashMap::new();
    for record in product.trade_exercises.iter() {
        seen.insert((&record.label, &record.at));
        if record.effective_at <= as_of {
            effective
                .entry(record.label.clone())
                .or_default()
                .insert(record.at, record.option);
        }
    }
    for record in records {
        let ExerciseRecord {
            label,
            at,
            effective_at,
            option,
        } = record;
        let count = product
            .references
            .choices
            .get(label)
            .and_then(|times| times.get(at));
        if count.is_none_or(|count| option >= count)
            || effective_at > at
            || !seen.insert((label, at))
        {
            return Err(Error::InvalidData(format!(
                "invalid exercise record {label:?} at {at}"
            )));
        }
        if *effective_at <= as_of {
            effective
                .entry(label.clone())
                .or_default()
                .insert(*at, *option);
        }
    }
    Ok(effective)
}

/// 簡約した値を使う文脈。
///
/// どちらの文脈でも、公表済みの fixing、推定値、判明した信用事由、効力のある行使記録を同じ規則で埋め込む。
/// 違いは失敗の扱いだけである。到達した計算（`Required`）は、欠損と定数演算の失敗をそのまま返す。
/// まだ選ばれるか分からない枝（`Deferred`）は、失敗したノードを残し、値が必要になった評価で報告させる。
#[derive(Clone, Copy)]
enum Reach {
    Required,
    Deferred,
}

/// 未確定の枝の欠損だけをノードのまま残す合図に変える。
///
/// 到達した計算の欠損と、欠損以外の失敗は呼び出し側に返す。
fn deferrable<T>(reach: Reach, result: Result<T>) -> Result<Option<T>> {
    match (reach, result) {
        (_, Ok(value)) => Ok(Some(value)),
        (Reach::Deferred, Err(Error::MissingData(_))) => Ok(None),
        (_, Err(error)) => Err(error),
    }
}

/// 未確定の枝の数値の失敗だけを、値が必要になった評価で報告する理由に変える。
///
/// 到達した計算の失敗と、数値以外の失敗（データの誤りなど）は呼び出し側に返す。
fn deferrable_numerical<T>(
    reach: Reach,
    result: Result<T>,
) -> Result<std::result::Result<T, Box<str>>> {
    match (reach, result) {
        (_, Ok(value)) => Ok(Ok(value)),
        (Reach::Deferred, Err(Error::Numerical(reason))) => Ok(Err(reason.into())),
        (_, Err(error)) => Err(error),
    }
}

/// 同じノードでも、到達した場合と将来の枝に残す場合で簡約結果が異なる。
struct Copies {
    required: Vec<Option<Ix>>,
    deferred: Vec<Option<Ix>>,
}

impl Copies {
    fn new(len: usize) -> Self {
        Self {
            required: vec![None; len],
            deferred: vec![None; len],
        }
    }

    fn get(&self, reach: Reach, ix: Ix) -> Option<Ix> {
        match reach {
            Reach::Required => self.required[ix as usize],
            Reach::Deferred => self.deferred[ix as usize],
        }
    }

    fn set(&mut self, reach: Reach, ix: Ix, result: Ix) {
        match reach {
            Reach::Required => self.required[ix as usize] = Some(result),
            Reach::Deferred => self.deferred[ix as usize] = Some(result),
        }
    }

    fn mapped(&self, reach: Reach, ix: Ix) -> Result<Ix> {
        self.get(reach, ix)
            .ok_or_else(|| Error::Internal(format!("unvisited IR node {ix}")))
    }
}

struct Reducer<'a> {
    product: &'a Product,
    as_of: Instant,
    fixings: &'a FixingStore,
    exercises: ByChoice<usize>,
    settlement: SettlementBoundary,
    contracts: Arena<ContractNode>,
    exprs: Arena<ExprNode>,
    contract_ids: Copies,
    expr_ids: Copies,
    // 鍵は元の商品の中の指標 ID への参照。ID は全取引が共有する Arc<str> なので、日次の行ごとに
    // 複製すると並列の準備でワーカーどうしが同じ参照数を更新し合う。参照なら複製しない。
    observations: HashMap<&'a IndexId, HashMap<ObservationPoint, FixingResolution>>,
    forecasts: ByIndex<ObservationPoint, Instant>,
    warnings: BTreeSet<Warning>,
}

enum ContractWork {
    Visit(Ix, Reach),
    Write(Ix, Reach),
    /// 実現した条件または効力のある行使記録が選んだ枝を、元のノードの簡約結果にする。
    Alias {
        ix: Ix,
        selected: Ix,
        reach: Reach,
    },
}

// 深い DAG でもスタックの深さに依存しないよう、式の簡約は一つの作業スタックで進める。
// 未確定の腕や項も同じスタックに `Deferred` として積み、入れ子の簡約を再帰で呼ばない。
enum ExprWork {
    Visit(Ix, Reach),
    Fold(Ix, Reach),
    /// 条件を簡約した後に、選ばれた腕に進むか、両腕を残すかを決める。
    Select(Ix, Reach),
    /// 決まらない `Select` の両腕を `Deferred` で簡約した後に組み直す。
    Branches(Ix, Reach),
    Logic {
        ix: Ix,
        pos: usize,
        reach: Reach,
    },
    /// 結果の決まらない項より後ろを `Deferred` で簡約した後に、論理演算を組み直す。
    Rest {
        ix: Ix,
        pos: usize,
        reach: Reach,
    },
    Alias {
        ix: Ix,
        selected: Ix,
        reach: Reach,
    },
}

impl<'a> Reducer<'a> {
    fn contract(&mut self) -> Result<Ix> {
        let root = self.product.contract.root;
        let nodes = &self.product.contract.contracts;
        let mut work = vec![ContractWork::Visit(root, Reach::Required)];
        while let Some(task) = work.pop() {
            match task {
                ContractWork::Visit(ix, reach) => {
                    if self.contract_ids.get(reach, ix).is_none() {
                        self.visit_contract(&nodes[ix as usize], ix, reach, &mut work)?;
                    }
                }
                ContractWork::Write(ix, reach) => {
                    let new = self.write_contract(&nodes[ix as usize], reach)?;
                    self.contract_ids.set(reach, ix, new);
                }
                ContractWork::Alias {
                    ix,
                    selected,
                    reach,
                } => {
                    let branch = self.contract_ids.mapped(reach, selected)?;
                    self.contract_ids.set(reach, ix, branch);
                }
            }
        }
        self.contract_ids.mapped(Reach::Required, root)
    }

    fn visit_contract(
        &mut self,
        node: &ContractNode,
        ix: Ix,
        reach: Reach,
        work: &mut Vec<ContractWork>,
    ) -> Result<()> {
        let select = |selected, work: &mut Vec<_>| {
            work.push(ContractWork::Alias {
                ix,
                selected,
                reach,
            });
            work.push(ContractWork::Visit(selected, reach));
        };
        // 未確定の枝の受渡と判断は、因果性の検査により外側の判断より後、つまり評価時点より後にある。
        // そのため決済済みの除外、判定済みの Cond、締切を過ぎた Choose は到達した経路でだけ起きる。
        let child_reach = match node {
            ContractNode::Pay { at, amount, .. } => {
                if !self.settlement.includes(self.as_of, *at) {
                    // 決済済みの受渡は額を読まずに除く。
                    let zero = self.contracts.insert(ContractNode::All(Vec::new()))?;
                    self.contract_ids.set(reach, ix, zero);
                    return Ok(());
                }
                self.expr(*amount, reach)?;
                reach
            }
            ContractNode::Cond { cond, at, yes, no } => {
                // 条件は Cond と同じ文脈で簡約する。到達した Cond の条件は判定時刻に必ず評価するので、
                // 判定が評価時点より後でも Required のまま、欠損をエラーにする。
                let cond = self.expr(*cond, reach)?;
                if *at <= self.as_of {
                    let ExprNode::Bool(value) = self.exprs.nodes[cond as usize] else {
                        return Err(Error::Internal(format!(
                            "condition at {at} remains unresolved"
                        )));
                    };
                    select(if value { *yes } else { *no }, work);
                    return Ok(());
                }
                Reach::Deferred
            }
            ContractNode::Choose {
                label, at, options, ..
            } => {
                // 効力のある記録は、外側の判断が未確定でもこの判断の結果を決める。
                if let Some(option) = self.exercises.get(label).and_then(|times| times.get(at)) {
                    select(options[*option], work);
                    return Ok(());
                }
                if *at <= self.as_of {
                    return Err(Error::MissingData(format!(
                        "missing exercise record {label:?} at {at}"
                    )));
                }
                Reach::Deferred
            }
            ContractNode::All(_) => reach,
        };
        work.push(ContractWork::Write(ix, reach));
        // 兄弟は元の順でたどる。選択肢の位置は行使記録と対応する。
        match node {
            ContractNode::Cond { yes, no, .. } => {
                work.push(ContractWork::Visit(*no, child_reach));
                work.push(ContractWork::Visit(*yes, child_reach));
            }
            _ => work.extend(
                node.children()
                    .iter()
                    .rev()
                    .map(|&child| ContractWork::Visit(child, child_reach)),
            ),
        }
        Ok(())
    }

    fn write_contract(&mut self, node: &ContractNode, reach: Reach) -> Result<Ix> {
        // 未確定の判断の枝と選択肢は、どれを受け取るか分からないので Deferred の簡約結果を使う。
        let branch = |ix| self.contract_ids.mapped(Reach::Deferred, ix);
        let node = match node {
            ContractNode::Pay { asset, amount, at } => ContractNode::Pay {
                asset: asset.clone(),
                amount: self.expr_ids.mapped(reach, *amount)?,
                at: *at,
            },
            ContractNode::All(children) => {
                let mut children = children
                    .iter()
                    .map(|&ix| self.contract_ids.mapped(reach, ix))
                    .collect::<Result<Vec<_>>>()?;
                if matches!(reach, Reach::Required) {
                    children.retain(|&child| !matches!(&self.contracts.nodes[child as usize], ContractNode::All(xs) if xs.is_empty()));
                    if let [only] = children.as_slice() {
                        return Ok(*only);
                    }
                }
                ContractNode::All(children)
            }
            ContractNode::Cond { cond, at, yes, no } => ContractNode::Cond {
                cond: self.expr_ids.mapped(reach, *cond)?,
                at: *at,
                yes: branch(*yes)?,
                no: branch(*no)?,
            },
            ContractNode::Choose {
                by,
                label,
                at,
                options,
            } => ContractNode::Choose {
                by: *by,
                label: label.clone(),
                at: *at,
                options: options
                    .iter()
                    .map(|&ix| branch(ix))
                    .collect::<Result<_>>()?,
            },
        };
        self.contracts.insert(node)
    }

    fn expr(&mut self, root: Ix, root_reach: Reach) -> Result<Ix> {
        // 元の商品のノードを 'a のまま借り、観測の表の鍵に ID の参照を使えるようにする。
        let source: &'a [ExprNode] = &self.product.contract.exprs;
        let mut work = vec![ExprWork::Visit(root, root_reach)];
        while let Some(task) = work.pop() {
            let (ix, reach, result) = match task {
                ExprWork::Visit(ix, reach) => {
                    if self.expr_ids.get(reach, ix).is_none() {
                        self.visit_expr(&source[ix as usize], ix, reach, &mut work)?;
                    }
                    continue;
                }
                ExprWork::Fold(ix, reach) => {
                    let mut node = source[ix as usize].clone();
                    node.try_map_children(|child| self.expr_ids.mapped(reach, child))?;
                    (ix, reach, self.fold_expr(node, reach)?)
                }
                ExprWork::Select(ix, reach) => {
                    let [cond, yes, no] = select_children(&source[ix as usize])?;
                    let cond = self.expr_ids.mapped(reach, cond)?;
                    if let ExprNode::Bool(value) = self.exprs.nodes[cond as usize] {
                        let selected = if value { yes } else { no };
                        work.push(ExprWork::Alias {
                            ix,
                            selected,
                            reach,
                        });
                        work.push(ExprWork::Visit(selected, reach));
                    } else {
                        // どちらの腕も選ばれうるので、欠損と定義域外の演算を残せる Deferred で簡約する。
                        work.push(ExprWork::Branches(ix, reach));
                        work.push(ExprWork::Visit(no, Reach::Deferred));
                        work.push(ExprWork::Visit(yes, Reach::Deferred));
                    }
                    continue;
                }
                ExprWork::Branches(ix, reach) => {
                    let [cond, yes, no] = select_children(&source[ix as usize])?;
                    let node = ExprNode::Select([
                        self.expr_ids.mapped(reach, cond)?,
                        self.expr_ids.mapped(Reach::Deferred, yes)?,
                        self.expr_ids.mapped(Reach::Deferred, no)?,
                    ]);
                    (ix, reach, self.exprs.insert(node)?)
                }
                ExprWork::Logic { ix, pos, reach } => {
                    let Some(result) =
                        self.logic(&source[ix as usize], ix, pos, reach, &mut work)?
                    else {
                        continue;
                    };
                    (ix, reach, result)
                }
                ExprWork::Rest { ix, pos, reach } => {
                    let ExprNode::Nary(op, xs) = &source[ix as usize] else {
                        return Err(Error::Internal(
                            "logic work requires a logical expression".into(),
                        ));
                    };
                    let mut children = Vec::with_capacity(xs.len() - pos);
                    children.push(self.expr_ids.mapped(reach, xs[pos])?);
                    for &child in &xs[pos + 1..] {
                        children.push(self.expr_ids.mapped(Reach::Deferred, child)?);
                    }
                    (ix, reach, self.exprs.insert(ExprNode::Nary(*op, children))?)
                }
                ExprWork::Alias {
                    ix,
                    selected,
                    reach,
                } => (ix, reach, self.expr_ids.mapped(reach, selected)?),
            };
            self.expr_ids.set(reach, ix, result);
        }
        self.expr_ids.mapped(root_reach, root)
    }

    fn visit_expr(
        &mut self,
        node: &'a ExprNode,
        ix: Ix,
        reach: Reach,
        work: &mut Vec<ExprWork>,
    ) -> Result<()> {
        let source: &'a [ExprNode] = &self.product.contract.exprs;
        match node {
            ExprNode::WithWarnings(warnings, [child]) => {
                self.warnings.extend(warnings.iter().cloned());
                work.push(ExprWork::Alias {
                    ix,
                    selected: *child,
                    reach,
                });
                work.push(ExprWork::Visit(*child, reach));
            }
            ExprNode::Select([cond, ..]) => {
                work.push(ExprWork::Select(ix, reach));
                work.push(ExprWork::Visit(*cond, reach));
            }
            ExprNode::Nary(Nary::And | Nary::Or, _) => {
                work.push(ExprWork::Logic { ix, pos: 0, reach });
            }
            ExprNode::Fixing { index, point, .. } => {
                let resolved = deferrable(reach, self.fixing(index, *point, None))?.flatten();
                let node = resolved
                    .map_or_else(|| node.clone(), |value| ExprNode::Const(Number::new(value)));
                self.store_expr(ix, reach, node)?;
            }
            ExprNode::Compounded { .. } => {
                let result = match self.fold_compounded(node, reach)? {
                    FoldedCompound::Known(excess) => known_constant(reach, 1.0 + excess)?,
                    FoldedCompound::Residual(residual) => *residual,
                    FoldedCompound::Undefined(reason) => ExprNode::Undefined(reason),
                };
                self.store_expr(ix, reach, result)?;
            }
            // 複利の超過成長 Compounded − 1 は、成長率を丸めてから 1 を引くと短い期間の下位の桁を
            // ulp(1) で失う（1 日の 0.5% で相対 7e-12 ほど。円の丸めの境界で 1 円ずれる）。
            // 確定因子は超過のまま畳んで定数にし、残りがあれば確定部分の超過を Compounded に持たせて評価に渡す。
            ExprNode::Binary(Binary::Sub, [growth, one])
                if matches!(source[*growth as usize], ExprNode::Compounded { .. })
                    && source[*one as usize] == ExprNode::Const(Number::new(1.0)) =>
            {
                let result = match self.fold_compounded(&source[*growth as usize], reach)? {
                    FoldedCompound::Known(excess) => known_constant(reach, excess)?,
                    FoldedCompound::Residual(residual) => {
                        let growth = self.exprs.insert(*residual)?;
                        let one = self.exprs.insert(ExprNode::Const(Number::new(1.0)))?;
                        ExprNode::Binary(Binary::Sub, [growth, one])
                    }
                    FoldedCompound::Undefined(reason) => ExprNode::Undefined(reason),
                };
                self.store_expr(ix, reach, result)?;
            }
            ExprNode::Average {
                index,
                rows,
                daily_rounding,
            } => {
                let folded = self
                    .partition_rows(
                        index,
                        rows.iter()
                            .copied()
                            .enumerate()
                            .map(|(pos, row)| (pos, row, None)),
                        *daily_rounding,
                        reach,
                        |value, weight| {
                            checked_float_product(value, weight).map_err(fold::numerical)
                        },
                    )
                    .and_then(|(terms, pending)| {
                        Ok((checked_sum(&terms).map_err(fold::numerical)?, pending))
                    });
                let (fixed, pending) = match deferrable_numerical(reach, folded)? {
                    Ok(folded) => folded,
                    Err(reason) => return self.store_expr(ix, reach, ExprNode::Undefined(reason)),
                };
                let residual = ExprNode::Average {
                    index: index.clone(),
                    rows: pending.iter().map(|&pos| rows[pos]).collect(),
                    daily_rounding: *daily_rounding,
                };
                self.combine_average(ix, reach, fixed, residual)?;
            }
            ExprNode::Defaulted { index, window } | ExprNode::Recovery { index, window } => {
                let result = deferrable(reach, self.credit(node, index, *window))?
                    .unwrap_or_else(|| node.clone());
                self.store_expr(ix, reach, result)?;
            }
            _ => {
                work.push(ExprWork::Fold(ix, reach));
                work.extend(
                    node.children()
                        .iter()
                        .rev()
                        .map(|&child| ExprWork::Visit(child, reach)),
                );
            }
        }
        Ok(())
    }

    fn logic(
        &mut self,
        node: &ExprNode,
        ix: Ix,
        pos: usize,
        reach: Reach,
        work: &mut Vec<ExprWork>,
    ) -> Result<Option<Ix>> {
        let ExprNode::Nary(op, xs) = node else {
            return Err(Error::Internal(
                "logic work requires a logical expression".into(),
            ));
        };
        let identity = *op == Nary::And;
        let Some(&child) = xs.get(pos) else {
            return Ok(Some(self.exprs.insert(ExprNode::Bool(identity))?));
        };
        let Some(result) = self.expr_ids.get(reach, child) else {
            work.push(ExprWork::Logic { ix, pos, reach });
            work.push(ExprWork::Visit(child, reach));
            return Ok(None);
        };
        match self.exprs.nodes[result as usize] {
            ExprNode::Bool(value) if value == identity => {
                work.push(ExprWork::Logic {
                    ix,
                    pos: pos + 1,
                    reach,
                });
                Ok(None)
            }
            ExprNode::Bool(_) => Ok(Some(result)),
            _ if pos + 1 == xs.len() => Ok(Some(result)),
            _ => {
                // ここで止まる可能性があるので、後ろの項は欠損を残せる Deferred で簡約する。
                work.push(ExprWork::Rest { ix, pos, reach });
                work.extend(
                    xs[pos + 1..]
                        .iter()
                        .rev()
                        .map(|&child| ExprWork::Visit(child, Reach::Deferred)),
                );
                Ok(None)
            }
        }
    }

    fn store_expr(&mut self, ix: Ix, reach: Reach, node: ExprNode) -> Result<()> {
        let copy = self.exprs.insert(node)?;
        self.expr_ids.set(reach, ix, copy);
        Ok(())
    }

    fn fold_expr(&mut self, node: ExprNode, reach: Reach) -> Result<Ix> {
        let constant = |ix| match self.exprs.nodes[ix as usize] {
            ExprNode::Const(value) => Some(value.value()),
            _ => None,
        };
        let folded = match &node {
            ExprNode::Unary(Unary::Not, [x]) => Ok(match self.exprs.nodes[*x as usize] {
                ExprNode::Bool(value) => Some(ExprNode::Bool(!value)),
                _ => None,
            }),
            ExprNode::Unary(op, [x]) => constant(*x).map(|x| fold::unary(*op, x)).transpose(),
            ExprNode::Binary(op, [left, right]) => constant(*left)
                .zip(constant(*right))
                .map(|(left, right)| fold::binary(*op, left, right))
                .transpose(),
            ExprNode::Nary(op, xs) => xs
                .iter()
                .map(|&ix| constant(ix))
                .collect::<Option<Vec<_>>>()
                .map(|xs| fold::nary(*op, &xs))
                .transpose(),
            _ => Ok(None),
        };
        let folded = match (reach, folded) {
            (_, Ok(folded)) => folded,
            // 選ばれるか分からない枝の部分関数（負の数の平方根など）は、その枝を選んだ評価で失敗させる。
            // 子は定数化済みなので、畳まずに残しても観測の情報は失わない。
            (Reach::Deferred, Err(Error::Numerical(_))) => None,
            (_, Err(error)) => return Err(error),
        };
        self.exprs.insert(folded.unwrap_or(node))
    }

    fn fixing(
        &mut self,
        index: &'a IndexId,
        point: ObservationPoint,
        reference: Option<&ResolvedFixing>,
    ) -> Result<Option<f64>> {
        let mut resolve = |reference: &ResolvedFixing| {
            if reference.ceased {
                return Err(Error::MissingData(format!(
                    "no fallback for ceased index {index} {point:?}"
                )));
            }
            if reference.non_representative {
                self.warnings.insert(Warning::NonRepresentativeIndex {
                    index: index.clone(),
                });
            }
            self.warnings.extend(reference.warnings.iter().cloned());
            let (resolution, warning) = self.fixings.resolve(
                index,
                point,
                self.as_of,
                reference.publication,
                reference.period,
            )?;
            self.warnings.extend(warning);
            Ok(resolution)
        };
        let resolution = if let Some(reference) = reference {
            resolve(reference)?
        } else {
            let lookup = self.observations.entry(index).or_default();
            match lookup.entry(point) {
                std::collections::hash_map::Entry::Occupied(entry) => *entry.get(),
                std::collections::hash_map::Entry::Vacant(entry) => {
                    let reference = &self.product.references.fixings[index][&point];
                    *entry.insert(resolve(reference)?)
                }
            }
        };
        Ok(match resolution {
            FixingResolution::Published(value) | FixingResolution::Estimated(value) => Some(value),
            FixingResolution::Forecast => None,
            FixingResolution::ForecastFrom(from) => {
                self.forecasts
                    .entry(index.clone())
                    .or_default()
                    .insert(point, from);
                None
            }
        })
    }

    fn partition_rows(
        &mut self,
        index: &'a IndexId,
        rows: impl IntoIterator<Item = (usize, Row, Option<&'a ResolvedFixing>)>,
        rounding: Option<Rounding>,
        reach: Reach,
        term: impl Fn(f64, f64) -> Result<f64>,
    ) -> Result<(Vec<f64>, Vec<usize>)> {
        let mut fixed = Vec::new();
        let mut pending = Vec::new();
        for (pos, (point, weight), reference) in rows {
            // 未確定の枝で欠けた行は、予想する行と同じく元の重みのまま残りの観測にする。
            // 確定値から項を作れない数値の失敗（丸めのあふれ、契約の Inclusive のスプレッドで
            // 非正になる因子など）はここで返し、未確定の枝なら呼び出し側が集約を Undefined にする。
            match deferrable(reach, self.fixing(index, point, reference))?.flatten() {
                Some(value) => {
                    let value = rounding.map_or(Ok(value), |rounding| rounding.apply(value))?;
                    fixed.push(term(value, weight.value())?);
                }
                None => pending.push(pos),
            }
        }
        Ok((fixed, pending))
    }

    /// 確定項の和と残りの観測の平均を組み合わせる。残りが無ければ定数、和が 0 なら残りだけにする。
    fn combine_average(
        &mut self,
        ix: Ix,
        reach: Reach,
        fixed: f64,
        residual: ExprNode,
    ) -> Result<()> {
        let ExprNode::Average { rows, .. } = &residual else {
            return Err(Error::Internal("expected an observation average".into()));
        };
        let value = Number::new(fixed);
        if rows.is_empty() {
            self.store_expr(ix, reach, ExprNode::Const(value))?;
        } else if value == Number::new(0.0) {
            self.store_expr(ix, reach, residual)?;
        } else {
            let known = self.exprs.insert(ExprNode::Const(value))?;
            let pending = self.exprs.insert(residual)?;
            self.store_expr(ix, reach, ExprNode::Nary(Nary::Add, vec![known, pending]))?;
        }
        Ok(())
    }

    /// `Compounded` の確定した観測を、積の超過（確定因子の積から 1 を引いた値）として畳む。
    ///
    /// 因子 1 + x を丸めてから掛けると x の下位の桁を失うので、因子の超過 x（`factor_excess`）を
    /// (1 + e)(1 + x) − 1 = e + x + e·x で順に畳む。単利の因子では四則演算だけなので、確定した受渡額は
    /// 環境（libm）によらず同じ値になり、丸めの境界の判定も環境で変わらない。
    /// 確定部分の失敗（丸めのあふれ、非正の因子、積のあふれ）は、到達した計算なら返し、
    /// 未確定の枝なら理由に変える。
    fn fold_compounded(&mut self, node: &'a ExprNode, reach: Reach) -> Result<FoldedCompound> {
        let ExprNode::Compounded {
            index,
            rows,
            dates,
            rule,
            daily_spread,
            fixed_excess,
        } = node
        else {
            return Err(Error::Internal("expected a compounded observation".into()));
        };
        let Aggregation::Compound { form, .. } = rule.aggregation() else {
            return Err(Error::Internal("invalid compounding rule".into()));
        };
        let resolved = &self.product.references.compounded[index][&rows.source];
        let observations = rows.positions().map(|pos| {
            (
                pos,
                rows.source.entries[pos].0,
                Some(&resolved.observations[pos]),
            )
        });
        let folded = self
            .partition_rows(
                index,
                observations,
                rule.rounding().daily,
                reach,
                |rate, weight| form.factor_excess(rate, weight, daily_spread.value()),
            )
            .and_then(|(increments, pending)| {
                // 先に畳んだ確定因子（同じ商品を簡約し直した場合）の超過から始める。
                let excess =
                    increments
                        .iter()
                        .try_fold(fixed_excess.value(), |excess, &increment| {
                            let next = excess + increment + excess * increment;
                            if next.is_finite() && next > -1.0 {
                                Ok(next)
                            } else {
                                Err(fold::numerical(
                                    "published compounding factors overflow or are not positive",
                                ))
                            }
                        })?;
                Ok((excess, increments.len(), pending))
            });
        let (excess, known, pending) = match deferrable_numerical(reach, folded)? {
            Ok(folded) => folded,
            Err(reason) => return Ok(FoldedCompound::Undefined(reason)),
        };
        if pending.is_empty() {
            return Ok(FoldedCompound::Known(excess));
        }
        // 残った最初の位置が、今回確定した先頭の行の直後かを確かめる。
        let fixed_prefix =
            dates.fixed_prefix && rows.positions().nth(known) == pending.first().copied();
        Ok(FoldedCompound::Residual(Box::new(ExprNode::Compounded {
            index: index.clone(),
            rows: rows.with_pending(pending),
            dates: CompoundingDates {
                basis: dates.basis.clone(),
                fixed_prefix,
            },
            rule: rule.clone(),
            daily_spread: *daily_spread,
            fixed_excess: Number::new(excess),
        })))
    }

    fn credit(
        &mut self,
        node: &ExprNode,
        index: &IndexId,
        window: DateInterval,
    ) -> Result<ExprNode> {
        let reference = &self.product.references.credit[index][&window];
        self.warnings.extend(reference.warnings.iter().cloned());
        let confirmed = self.as_of >= reference.publication;
        let event = self.fixings.credit_event(&reference.key, window);
        match node {
            ExprNode::Defaulted { .. } if event.is_some() || confirmed => {
                Ok(ExprNode::Bool(event.is_some()))
            }
            ExprNode::Recovery { .. } if confirmed => event
                .and_then(|_| self.fixings.recovery(&reference.key))
                .map(|recovery| ExprNode::Const(Number::new(recovery.value())))
                .ok_or_else(|| {
                    Error::MissingData(format!("missing recovery {index} for {window:?}"))
                }),
            _ => Ok(node.clone()),
        }
    }
}

/// 確定した観測を畳んだ `Compounded` の結果。
enum FoldedCompound {
    /// すべての観測が確定した。値は確定因子の積から 1 を引いた値。
    Known(f64),
    /// 未確定の観測が残る。確定因子の超過は `fixed_excess` に入れてある。
    /// 他の場合より大きいので、箱に入れて返り値を小さく保つ。
    Residual(Box<ExprNode>),
    /// 未確定の枝で、確定部分から値を作れなかった（理由を持つ）。
    Undefined(Box<str>),
}

/// 確定した値の定数ノードを作る。値域の外は、未確定の枝なら理由付きの `Undefined` にする。
fn known_constant(reach: Reach, value: f64) -> Result<ExprNode> {
    Ok(match deferrable_numerical(reach, fold::constant(value))? {
        Ok(node) => node,
        Err(reason) => ExprNode::Undefined(reason),
    })
}

fn select_children(node: &ExprNode) -> Result<[Ix; 3]> {
    match node {
        ExprNode::Select(children) => Ok(*children),
        _ => Err(Error::Internal("select work requires a select node".into())),
    }
}
`````

`````
crates/rhodium-products/src/ir/reduce/fold.rs
use rhodium_core::{Error, Result};

use super::{
    Binary, ExprNode, Nary, Number, Unary, checked_float_product_factors, checked_float_quotient,
    checked_sum,
};

pub(super) fn unary(op: Unary, x: f64) -> Result<ExprNode> {
    let result = match op {
        Unary::Neg => -x,
        Unary::Abs => x.abs(),
        Unary::Exp => x.exp(),
        Unary::Ln => x.ln(),
        Unary::Sqrt => x.sqrt(),
        Unary::Round(rounding) => rounding.apply(x)?,
        Unary::Not => return Err(Error::Internal("expected numeric unary operator".into())),
    };
    constant(result)
}

pub(super) fn binary(op: Binary, left: f64, right: f64) -> Result<ExprNode> {
    let result = match op {
        Binary::Sub => left - right,
        Binary::Div => checked_float_quotient(left, right).map_err(numerical)?,
        Binary::Pow => left.powf(right),
        Binary::Lt => return Ok(ExprNode::Bool(left < right)),
        Binary::Le => return Ok(ExprNode::Bool(left <= right)),
        Binary::Gt => return Ok(ExprNode::Bool(left > right)),
        Binary::Ge => return Ok(ExprNode::Bool(left >= right)),
    };
    constant(result)
}

pub(super) fn nary(op: Nary, values: &[f64]) -> Result<ExprNode> {
    let result = match op {
        Nary::Add => checked_sum(values).map_err(numerical)?,
        Nary::Mul => checked_float_product_factors(values).map_err(numerical)?,
        Nary::Min => values.iter().copied().fold(f64::INFINITY, f64::min),
        Nary::Max => values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        Nary::And | Nary::Or => {
            return Err(Error::Internal("expected numeric n-ary operator".into()));
        }
    };
    constant(result)
}

pub(super) fn constant(value: f64) -> Result<ExprNode> {
    if value.is_finite() {
        Ok(ExprNode::Const(Number::new(value)))
    } else {
        Err(numerical("constant expression has a nonfinite result"))
    }
}

pub(super) fn numerical(error: impl std::fmt::Display) -> Error {
    Error::Numerical(error.to_string())
}
`````

`````
crates/rhodium-products/src/ir/tests.rs
use std::num::NonZeroU8;

use rhodium_core::{
    Error, Result, Warning,
    calendar::{BaseCalendar, CalendarSpec, Extrapolation, WeekendRule},
    currency::{Asset, Currency},
    date::{DateInterval, Weekday},
    day_count::{AccrualContext, DayCount},
    fixing::{CreditEventRecord, FixingStore, FixingStoreParts, ObservationPoint},
    id::{CommodityId, CreditKey, IndexId},
    observation::{
        Aggregation, AveragingWeights, CompoundingForm, ObservationDays, ObservationMethod,
        ObservationRounding, ObservationRule, ObservationSchedule, SpreadTreatment,
        observation_schedule,
    },
    quantity::{RecoveryRate, Spread},
    refdata::{
        CommodityDefinition, CurrencyDefinition, DailyPublication, Deliverability, IndexDefinition,
        IndexKind, IndexStatus, PublicationRule, RefData, RefDataParts,
    },
    rounding::{Rounding, RoundingMode},
    schedule::{
        ExplicitDatePolicy, ExplicitPeriod, ExplicitPeriods, PeriodCombination, ScheduleSpec,
        SettlementContext, generate,
    },
    time::{Instant, LocalEventTime, UtcOffset, ZoneRules},
};
use rhodium_testkit::assert_close;

use super::*;
use crate::trade::ExerciseRecord;

const CASH: Asset = Asset::Cash(Currency::JPY);

fn credit_key() -> Result<CreditKey> {
    Ok(CreditKey {
        entity: "EXAMPLE".parse()?,
        seniority: "SENIOR".parse()?,
        doc_clause: "CR".parse()?,
        ccy: Currency::JPY,
    })
}

fn refdata() -> Result<RefData> {
    RefData::new(refdata_parts()?)
}

/// USDJPY の受渡慣行（T+2、pip 0.01 円）を、通貨と同じ決済時刻で足す。
fn add_usd_jpy(parts: &mut RefDataParts, clock: &LocalEventTime) -> Result<()> {
    let usd_jpy: rhodium_core::fx::CcyPair = "USDJPY".parse()?;
    parts.fx_conventions.insert(
        usd_jpy,
        rhodium_core::refdata::FxConvention {
            pair: usd_jpy,
            spot_lag: rhodium_core::refdata::FxSpotLag::TwoBusinessDays,
            settlement: clock.clone(),
            pip_exponent: 2,
        },
    );
    Ok(())
}

fn refdata_parts() -> Result<RefDataParts> {
    let mut parts = RefDataParts::default();
    let calendar = "JPTO".parse()?;
    let clock = LocalEventTime {
        time: "09:00:00".parse()?,
        zone: "Asia/Tokyo".parse()?,
    };
    parts.calendars.insert(
        calendar,
        BaseCalendar::new(
            "2027-01-01".parse()?..="2027-12-31".parse()?,
            WeekendRule::new([Weekday::Sat, Weekday::Sun]),
            Vec::new(),
            Extrapolation::Strict,
        )?,
    );
    parts.zones.insert(
        clock.zone.clone(),
        ZoneRules::new(
            "2027-01-01T00:00:00Z".parse()?..="2027-12-31T23:59:59Z".parse()?,
            UtcOffset::from_seconds(9 * 3600)?,
            Vec::new(),
        )?,
    );
    for code in [Currency::JPY, Currency::USD] {
        parts.currencies.insert(
            code,
            CurrencyDefinition {
                code,
                minor_units: u8::from(code == Currency::USD) * 2,
                payment_rounding: Rounding::new(1.0, RoundingMode::HalfEven)?,
                deliverability: Deliverability::Deliverable,
                settlement_time: clock.clone(),
                settlement_calendar: "JPTO".parse()?,
            },
        );
    }
    let calendar = CalendarSpec::Base("JPTO".parse()?);
    let daily = |business_days| -> Result<_> {
        Ok(PublicationRule::Daily(DailyPublication {
            business_days,
            calendar: calendar.clone(),
            time: "10:00:00".parse()?,
            zone: clock.zone.clone(),
        }))
    };
    // USDJPY の FX fixing の受渡日は、ペアの慣行のスポットラグで決まる。
    add_usd_jpy(&mut parts, &clock)?;
    for (name, kind, publication, observation_time) in [
        (
            "JPY-TONA",
            IndexKind::OvernightRate {
                ccy: Currency::JPY,
                day_count: DayCount::Act365F,
            },
            daily(1)?,
            Some(clock.clone()),
        ),
        (
            "USDJPY-SPOT",
            IndexKind::FxFixing {
                pair: "USDJPY".parse()?,
            },
            PublicationRule::Intraday { delay_micros: 0 },
            None,
        ),
        (
            "JP-EVENT",
            IndexKind::CreditEvent {
                credit: credit_key()?,
            },
            daily(0)?,
            None,
        ),
        (
            "JP-RECOVERY",
            IndexKind::Recovery {
                credit: credit_key()?,
            },
            daily(1)?,
            None,
        ),
    ] {
        let id: IndexId = name.parse()?;
        parts.indices.insert(
            id.clone(),
            IndexDefinition {
                id,
                kind,
                calendar: calendar.clone(),
                publication,
                observation_time,
                status: IndexStatus::Active,
            },
        );
    }
    insert_commodity(&mut parts)?;
    Ok(parts)
}

// 現物受渡の例が使うコモディティ。
fn insert_commodity(parts: &mut RefDataParts) -> Result<()> {
    let oil: CommodityId = "WTI-CUSHING".parse()?;
    parts.commodities.insert(
        oil.clone(),
        CommodityDefinition {
            id: oil,
            unit: "bbl".into(),
            currency: Currency::USD,
        },
    );
    Ok(())
}

fn observation_rows(
    refdata: &RefData,
    method: ObservationMethod,
    aggregation: Aggregation,
    rounding: ObservationRounding,
) -> Result<ObservationSchedule> {
    let accrual = DateInterval::new("2027-01-08".parse()?, "2027-01-12".parse()?)?;
    let spec = ScheduleSpec::Explicit(ExplicitPeriods::new(
        vec![ExplicitPeriod {
            accrual,
            payment: accrual.end(),
            reset: None,
            settlement: None,
            references: Vec::new(),
        }],
        ExplicitDatePolicy::KeepAsConfirmed,
        PeriodCombination::Sum,
        None,
    )?);
    let calendar = CalendarSpec::Base("JPTO".parse()?);
    let (schedule, _) = generate(
        &spec,
        accrual.start(),
        refdata.calendars(),
        &SettlementContext {
            zones: refdata.zones(),
            default_time: Some(&refdata.currency(Currency::JPY)?.settlement_time),
            calendar: Some(&calendar),
        },
    )?;
    let rule = ObservationRule::new(
        ObservationDays::BusinessDays(CalendarSpec::Base("JPTO".parse()?)),
        method,
        aggregation,
        rounding,
    )?;
    observation_schedule(
        &schedule.periods()[0],
        &rule,
        &DayCount::Act365F,
        &AccrualContext::default(),
        refdata.calendars(),
    )
    .map(|result| result.0)
}

fn root(product: &Product) -> &ContractNode {
    &product.contract.contracts[product.contract.root as usize]
}

fn amount_expr(product: &Product) -> Result<&ExprNode> {
    let ContractNode::Pay { amount, .. } = root(product) else {
        return Err(Error::Internal("expected one payment".into()));
    };
    Ok(&product.contract.exprs[*amount as usize])
}

fn const_amount(product: &Product) -> Result<f64> {
    let ExprNode::Const(value) = amount_expr(product)? else {
        return Err(Error::Internal("expected a fixed payment".into()));
    };
    Ok(value.value())
}

fn cond_branches(product: &Product) -> Result<(Ix, Ix)> {
    let ContractNode::Cond { yes, no, .. } = root(product) else {
        return Err(Error::Internal("expected a pending condition".into()));
    };
    Ok((*yes, *no))
}

fn pay_amount(product: &Product, contract: Ix) -> Result<&ExprNode> {
    let ContractNode::Pay { amount, .. } = product.contract.contracts[contract as usize] else {
        return Err(Error::Internal("expected a payment".into()));
    };
    Ok(&product.contract.exprs[amount as usize])
}

fn points(rows: &nodes::CompoundingRows) -> Vec<ObservationPoint> {
    rows.positions()
        .map(|position| {
            let ((point, _), _) = rows.source.entries[position];
            point
        })
        .collect()
}

#[test]
fn matching_rfr_schedules_share_rows_and_resolved_fixings() -> Result<()> {
    // 同じ規約と期間の列を独立に作り、支払額を変えても列と公表時刻を共有する。
    // 観測規約を変えた契約では元の観測列を参照しない。
    let refdata = refdata()?;
    let cache = ObservationCache::default();
    let index: IndexId = "JPY-TONA".parse()?;
    let at = "2027-01-15T03:00:00Z".parse()?;
    let product = |method, notional| {
        let observations = observation_rows(
            &refdata,
            method,
            Aggregation::Compound {
                form: CompoundingForm::SimpleFactor,
                spread: SpreadTreatment::Exclusive,
            },
            ObservationRounding::default(),
        )?;
        let contract = Contract::build(|b| {
            let growth = b.compounded(
                &index,
                &observations,
                observations.rule(),
                Spread::new(0.0)?,
            )?;
            let notional = b.constant(notional);
            let amount = b.mul([growth, notional]);
            Ok(b.pay(CASH, amount, at))
        })?;
        Product::validated(contract, &refdata, Vec::new(), Some(&cache))
    };
    let left = product(ObservationMethod::Plain, 1_000_000.0)?;
    let right = product(ObservationMethod::Plain, 2_000_000.0)?;
    let lookback = product(ObservationMethod::Lookback(NonZeroU8::MIN), 1_000_000.0)?;
    let entries = |product: &Product| {
        let (source, resolved) = product.references.compounded[&index]
            .iter()
            .next()
            .ok_or_else(|| Error::Internal("missing RFR observations".into()))?;
        Ok::<_, Error>((Arc::clone(source), Arc::clone(resolved)))
    };
    let (rows, resolved) = entries(&left)?;
    let (same_rows, same_resolved) = entries(&right)?;
    let (different_rows, different_resolved) = entries(&lookback)?;
    assert!(Arc::ptr_eq(&rows, &same_rows));
    assert!(Arc::ptr_eq(&resolved, &same_resolved));
    assert!(!Arc::ptr_eq(&rows, &different_rows));
    assert!(!Arc::ptr_eq(&resolved, &different_resolved));
    Ok(())
}

#[test]
fn fixed_leg_keeps_shared_payments_and_respects_settlement_boundary() -> Result<()> {
    let at = "2027-01-08T03:00:00Z".parse::<Instant>()?;
    let later = "2027-07-08T03:00:00Z".parse()?;
    let refdata = refdata()?;
    let contract = Contract::build(|b| {
        let factors = [100_000.0, 0.03, 0.5].map(|value| b.constant(value));
        let amount = b.mul(factors);
        let pay = b.pay(CASH, amount, at);
        let same = b.pay(CASH, amount, at);
        assert_eq!(pay.ix, same.ix);
        let received = b.pay(CASH, amount, later);
        let paid = b.give(received);
        Ok(b.all([pay, same, paid]))
    })?;
    let product = Product::new(contract, &refdata)?;
    let fixings = FixingStore::default();
    // 受取は 100,000 円 × 3% × 1/2 を 2 回、後日の支払は同額を 1 回。
    let before = Instant::from_unix_micros(at.unix_micros() - 1)?;
    for (as_of, boundary) in [
        (before, SettlementBoundary::Exclude),
        (at, SettlementBoundary::Include),
    ] {
        let included = reduce(&product, as_of, &fixings, &[], boundary)?;
        let ContractNode::All(xs) = root(&included) else {
            return Err(Error::Internal("expected fixed leg".into()));
        };
        assert_eq!(xs.len(), 3);
        assert_eq!(xs[0], xs[1]);
        let total = xs
            .iter()
            .map(|&ix| {
                let ContractNode::Pay { amount, .. } = included.contract.contracts[ix as usize]
                else {
                    return Err(Error::Internal("expected payment".into()));
                };
                let ExprNode::Const(amount) = included.contract.exprs[amount as usize] else {
                    return Err(Error::Internal("expected each coupon to be folded".into()));
                };
                Ok(amount.value())
            })
            .sum::<Result<f64>>()?;
        assert_close!(total, 1_500.0, 1e-12, 1e-6);
    }
    let after = Instant::from_unix_micros(at.unix_micros() + 1)?;
    for (as_of, boundary) in [
        (at, SettlementBoundary::Exclude),
        (after, SettlementBoundary::Include),
    ] {
        let excluded = reduce(&product, as_of, &fixings, &[], boundary)?;
        assert_close!(const_amount(&excluded)?, -1_500.0, 1e-12, 1e-6);
    }
    Ok(())
}

#[test]
fn tona_compounding_keeps_weekend_weights_and_folds_the_published_factors() -> Result<()> {
    let refdata = refdata()?;
    let id: IndexId = "JPY-TONA".parse()?;
    let at = "2027-01-12T03:00:00Z".parse()?;
    let observations = observation_rows(
        &refdata,
        ObservationMethod::Lookback(NonZeroU8::MIN),
        Aggregation::Compound {
            form: CompoundingForm::SimpleFactor,
            spread: SpreadTreatment::Exclusive,
        },
        ObservationRounding::default(),
    )?;
    let thursday = ObservationPoint::Date("2027-01-07".parse()?);
    let friday = ObservationPoint::Date("2027-01-08".parse()?);
    let product = Product::new(
        Contract::build(|b| {
            let amount =
                b.compounded(&id, &observations, observations.rule(), Spread::new(0.005)?)?;
            Ok(b.pay(CASH, amount, at))
        })?,
        &refdata,
    )?;
    let fixings = FixingStore::new(
        &refdata,
        FixingStoreParts {
            published: vec![(id.clone(), thursday, 0.1), (id.clone(), friday, 0.2)],
            ..FixingStoreParts::default()
        },
    )?;
    // 金曜公表済みの木曜分は週末の 3 日。金曜の観測は期間中で 1 日分の重みを保つ。
    let as_of = "2027-01-08T01:00:00Z".parse()?;
    let partial = reduce(&product, as_of, &fixings, &[], SettlementBoundary::Exclude)?;
    // 確定した木曜分は積に丸めず、因子の超過 0.1·3/365 として残りの複利に持たせる。
    let ExprNode::Compounded {
        rows, fixed_excess, ..
    } = amount_expr(&partial)?
    else {
        return Err(Error::Internal(
            "expected residual compounding with the folded excess".into(),
        ));
    };
    assert_close!(fixed_excess.value(), 0.1 * 3.0 / 365.0, 1e-15, 0.0);
    assert_eq!(rows.len(), 1);
    let position = rows
        .positions()
        .next()
        .ok_or_else(|| Error::Internal("missing remaining observation".into()))?;
    let ((point, weight), _) = rows.source.entries[position];
    assert_eq!(point, friday);
    assert_close!(weight.value(), 1.0 / 365.0, 1e-12, 1e-14);
    assert_eq!(partial.forecast_from(&id, friday), Some(as_of));
    assert_eq!(
        partial.warnings(),
        &[Warning::UnpublishedFixingForecast {
            index: id,
            point: friday
        }]
    );

    let full = reduce(
        &product,
        "2027-01-11T01:00:00Z".parse()?,
        &fixings,
        &[],
        SettlementBoundary::Exclude,
    )?;
    let expected = (1.0 + 0.1 * 3.0 / 365.0) * (1.0 + 0.2 / 365.0);
    // 元本 3,650 円の利息としても比較する。Exclusive のスプレッドは日次因子に入らない。
    assert_close!(
        3650.0 * (const_amount(&full)? - 1.0),
        3650.0 * (expected - 1.0),
        1e-12,
        1e-6
    );
    assert_eq!(full.warnings(), []);
    Ok(())
}

#[test]
fn average_preserves_the_original_weights_after_partial_reduction() -> Result<()> {
    let refdata = refdata()?;
    let id: IndexId = "JPY-TONA".parse()?;
    let observations = observation_rows(
        &refdata,
        ObservationMethod::Lookback(NonZeroU8::MIN),
        Aggregation::Average(AveragingWeights::CoverWeighted),
        ObservationRounding::default(),
    )?;
    let payment = "2027-01-12T03:00:00Z".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let avg = b.average(&id, &observations)?;
            Ok(b.pay(CASH, avg, payment))
        })?,
        &refdata,
    )?;
    let thursday = ObservationPoint::Date("2027-01-07".parse()?);
    let friday = ObservationPoint::Date("2027-01-08".parse()?);
    let fixings = FixingStore::new(
        &refdata,
        FixingStoreParts {
            published: vec![(id.clone(), thursday, 0.1), (id, friday, 0.3)],
            ..FixingStoreParts::default()
        },
    )?;
    let partial = reduce(
        &product,
        "2027-01-08T01:00:00Z".parse()?,
        &fixings,
        &[],
        SettlementBoundary::Exclude,
    )?;
    let ExprNode::Nary(nodes::Nary::Add, xs) = amount_expr(&partial)? else {
        return Err(Error::Internal("expected partial average".into()));
    };
    let ExprNode::Average { rows, .. } = &partial.contract.exprs[xs[1] as usize] else {
        return Err(Error::Internal("expected remaining term".into()));
    };
    // 4 日中 1 日分のまま残り、残りの 1 観測だけで重みを振り直さない。
    assert_eq!(rows.len(), 1);
    assert_close!(rows[0].1.value(), 0.25, 1e-12, 1e-14);
    let full = reduce(
        &product,
        "2027-01-11T01:00:00Z".parse()?,
        &fixings,
        &[],
        SettlementBoundary::Exclude,
    )?;
    assert_close!(const_amount(&full)?, 0.1 * 0.75 + 0.3 * 0.25, 1e-12, 1e-14);
    Ok(())
}

#[test]
fn fixing_switches_at_publication_and_past_unpublished_requires_an_estimate() -> Result<()> {
    let refdata = refdata()?;
    let id: IndexId = "USDJPY-SPOT".parse()?;
    let pub_at = "2027-01-08T01:00:00Z".parse::<Instant>()?;
    let before = Instant::from_unix_micros(pub_at.unix_micros() - 1)?;
    let after = Instant::from_unix_micros(pub_at.unix_micros() + 1)?;
    let at = "2027-01-12T03:00:00Z".parse()?;
    let point = ObservationPoint::Instant(pub_at);
    let product = Product::new(
        Contract::build(|b| {
            let obs = b.fixing(&id, point);
            Ok(b.pay(CASH, obs, at))
        })?,
        &refdata,
    )?;
    let actual = FixingStore::new(
        &refdata,
        FixingStoreParts {
            published: vec![(id, point, 150.0)],
            ..FixingStoreParts::default()
        },
    )?;
    let empty = FixingStore::default();
    // 公表値を先に渡しても観測前は予想を残す。公表の瞬間から実績が必要。
    let forecast = reduce(&product, before, &actual, &[], SettlementBoundary::Exclude)?;
    assert!(matches!(amount_expr(&forecast)?, ExprNode::Fixing { .. }));
    for as_of in [pub_at, after] {
        assert!(matches!(
            reduce(&product, as_of, &empty, &[], SettlementBoundary::Exclude),
            Err(Error::MissingData(_))
        ));
        let fixed = reduce(&product, as_of, &actual, &[], SettlementBoundary::Exclude)?;
        assert_close!(const_amount(&fixed)?, 150.0, 1e-12, 1e-6);
    }

    let id: IndexId = "JPY-TONA".parse()?;
    let friday = ObservationPoint::Date("2027-01-08".parse()?);
    let before = "2027-01-11T00:59:59.999999Z".parse()?;
    let publication = "2027-01-11T01:00:00Z".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let obs = b.fixing(&id, friday);
            Ok(b.pay(CASH, obs, at))
        })?,
        &refdata,
    )?;
    assert!(matches!(
        reduce(&product, before, &empty, &[], SettlementBoundary::Exclude),
        Err(Error::MissingData(_))
    ));
    let estimates = empty.with_estimates(&refdata, [(id.clone(), friday, 0.03)])?;
    let estimated = reduce(
        &product,
        before,
        &estimates,
        &[],
        SettlementBoundary::Exclude,
    )?;
    assert_close!(const_amount(&estimated)?, 0.03, 1e-12, 1e-14);
    assert_eq!(
        estimated.warnings(),
        &[Warning::FixingEstimateUsed {
            index: id,
            point: friday
        }]
    );
    assert!(matches!(
        reduce(
            &product,
            publication,
            &estimates,
            &[],
            SettlementBoundary::Exclude
        ),
        Err(Error::MissingData(_))
    ));
    Ok(())
}

#[test]
fn causality_uses_publication_instants_and_checks_each_reachable_branch() -> Result<()> {
    let refdata = refdata()?;
    let id: IndexId = "USDJPY-SPOT".parse()?;
    let publication = "2027-01-08T01:00:00Z".parse::<Instant>()?;
    let before = Instant::from_unix_micros(publication.unix_micros() - 1)?;
    let at = "2027-01-12T03:00:00Z".parse()?;
    for decision in [before, publication] {
        let contract = Contract::build(|b| {
            let spot = b.fixing(&id, ObservationPoint::Instant(publication));
            let strike = b.constant(100.0);
            let condition = b.gt(spot, strike);
            let pay = b.pay(CASH, strike, at);
            let zero = b.zero();
            Ok(b.cond(condition, decision, pay, zero))
        })?;
        assert_eq!(
            Product::new(contract, &refdata).is_ok(),
            decision == publication
        );
    }
    let invalid_amount = Contract::build(|b| {
        let spot = b.fixing(&id, ObservationPoint::Instant(publication));
        Ok(b.pay(CASH, spot, before))
    })?;
    assert!(matches!(
        Product::new(invalid_amount, &refdata),
        Err(Error::InvalidInput(_))
    ));
    for child_is_early in [true, false] {
        let contract = Contract::build(|b| {
            let spot = b.fixing(&id, ObservationPoint::Instant(publication));
            let pay = b.pay(CASH, spot, if child_is_early { publication } else { at });
            let zero = b.zero();
            // 支払額の公表を判断より後にしても、支払時刻までに公表されればよい。
            let decision = if child_is_early { at } else { before };
            Ok(b.choose(Party::Us, "exercise", decision, [pay, zero]))
        })?;
        assert_eq!(Product::new(contract, &refdata).is_ok(), !child_is_early);
    }
    Ok(())
}

#[test]
fn exercise_effectiveness_and_settled_payments_control_which_fixings_are_read() -> Result<()> {
    let refdata = refdata()?;
    let fixing = "2027-01-07T01:00:00Z".parse()?;
    let cutoff = "2027-01-08T01:00:00Z".parse::<Instant>()?;
    let future = "2027-01-12T03:00:00Z".parse()?;
    let id = "USDJPY-SPOT".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let missing = b.fixing(&id, ObservationPoint::Instant(fixing));
            let expired = b.pay(CASH, missing, cutoff);
            let bad = b.pay(CASH, missing, future);
            let amount = b.constant(110.0);
            let good = b.pay(CASH, amount, future);
            let choice = b.choose(Party::Us, "exercise", cutoff, [bad, good]);
            let sold = b.give(choice);
            Ok(b.all([expired, sold]))
        })?,
        &refdata,
    )?;
    let empty = FixingStore::default();
    let record = ExerciseRecord {
        label: "exercise".into(),
        at: cutoff,
        effective_at: cutoff,
        option: 1,
    };
    let before = Instant::from_unix_micros(cutoff.unix_micros() - 1)?;
    // 締切前には同時点の受渡が生きているため、まずその fixing が必要。
    let fixings = empty.with_assumed_published(
        &refdata,
        [(id, ObservationPoint::Instant(fixing), 100.0)],
        [],
    )?;
    let pending = reduce(
        &product,
        before,
        &fixings,
        std::slice::from_ref(&record),
        SettlementBoundary::Exclude,
    )?;
    let ContractNode::All(xs) = root(&pending) else {
        return Err(Error::Internal(
            "expected scheduled payment and pending choice".into(),
        ));
    };
    assert!(matches!(
        pending.contract.contracts[xs[1] as usize],
        ContractNode::Choose {
            by: Party::Counterparty,
            ..
        }
    ));
    assert!(matches!(
        reduce(&product, cutoff, &empty, &[], SettlementBoundary::Exclude),
        Err(Error::MissingData(_))
    ));
    let chosen = reduce(
        &product,
        cutoff,
        &empty,
        std::slice::from_ref(&record),
        SettlementBoundary::Exclude,
    )?;
    assert_close!(const_amount(&chosen)?, -110.0, 1e-12, 1e-6);
    assert!(matches!(
        reduce(
            &product,
            cutoff,
            &empty,
            std::slice::from_ref(&record),
            SettlementBoundary::Include
        ),
        Err(Error::MissingData(_))
    ));
    let other = ExerciseRecord {
        option: 0,
        ..record.clone()
    };
    assert!(matches!(
        reduce(
            &product,
            cutoff,
            &empty,
            &[other],
            SettlementBoundary::Exclude
        ),
        Err(Error::MissingData(_))
    ));
    Ok(())
}

#[test]
fn exercise_records_must_match_the_product_and_deadline() -> Result<()> {
    let refdata = refdata()?;
    let cutoff = "2027-01-08T01:00:00Z".parse()?;
    let future = "2027-01-12T03:00:00Z".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let value = b.constant(110.0);
            let pay = b.pay(CASH, value, future);
            let zero = b.zero();
            Ok(b.choose(Party::Us, "exercise", cutoff, [pay, zero]))
        })?,
        &refdata,
    )?;
    let empty = FixingStore::default();
    // キー、効力時刻、選択肢を商品に照合する。重複は選択肢が同じでも拒否する。
    let record = ExerciseRecord {
        label: "exercise".into(),
        at: cutoff,
        effective_at: cutoff,
        option: 1,
    };
    let invalid = ExerciseRecord {
        option: 2,
        ..record.clone()
    };
    let late = ExerciseRecord {
        effective_at: future,
        ..record.clone()
    };
    let unknown = ExerciseRecord {
        label: "other".into(),
        ..record.clone()
    };
    for records in [
        vec![invalid],
        vec![late],
        vec![unknown],
        vec![record.clone(), record],
    ] {
        assert!(matches!(
            reduce(
                &product,
                cutoff,
                &empty,
                &records,
                SettlementBoundary::Exclude
            ),
            Err(Error::InvalidData(_))
        ));
    }
    Ok(())
}

#[test]
fn cond_select_and_logic_skip_unreached_fixings_and_arithmetic() -> Result<()> {
    let refdata = refdata()?;
    let id: IndexId = "USDJPY-SPOT".parse()?;
    let today = "2027-01-08T01:00:00Z".parse()?;
    let future = "2027-01-11T01:00:00Z".parse()?;
    let payment = "2027-01-12T03:00:00Z".parse()?;
    let empty = FixingStore::default();
    for gate in [true, false] {
        let contract = Contract::build(|b| {
            let zero = b.constant(0.0);
            let missing = b.fixing(&id, ObservationPoint::Instant(today));
            let required = b.gt(missing, zero);
            let yes = b.boolean(true);
            let no = b.boolean(false);
            let short_and = b.and([no, required]);
            let short_or = b.or([yes, required]);
            let no_pay = b.pay(CASH, missing, payment);
            let value = b.constant(42.0);
            let negative = b.constant(-1.0);
            let undefined = b.sqrt(negative);
            let value = b.select(short_or, value, undefined);
            let conditional = b.select(short_and, missing, value);
            let pay = b.pay(CASH, conditional, payment);
            let gate = b.boolean(gate);
            Ok(b.cond(gate, today, pay, no_pay))
        })?;
        let product = Product::new(contract, &refdata)?;
        let reduced = reduce(&product, today, &empty, &[], SettlementBoundary::Exclude);
        if gate {
            assert_close!(const_amount(&reduced?)?, 42.0, 1e-12, 1e-6);
        } else {
            assert!(matches!(reduced, Err(Error::MissingData(_))));
        }
    }
    let product = Product::new(
        Contract::build(|b| {
            let late = b.fixing(&id, ObservationPoint::Instant(future));
            let missing = b.fixing(&id, ObservationPoint::Instant(today));
            let zero = b.constant(0.0);
            let condition = b.gt(late, zero);
            let amount = b.select(condition, missing, zero);
            Ok(b.pay(CASH, amount, payment))
        })?,
        &refdata,
    )?;
    let reduced = reduce(&product, today, &empty, &[], SettlementBoundary::Exclude)?;
    assert!(matches!(amount_expr(&reduced)?, ExprNode::Select(_)));
    for (spot, should_succeed) in [(-1.0, true), (1.0, false)] {
        let store = empty.with_assumed_published(
            &refdata,
            [(id.clone(), ObservationPoint::Instant(future), spot)],
            [],
        )?;
        assert_eq!(
            reduce(&product, future, &store, &[], SettlementBoundary::Exclude).is_ok(),
            should_succeed
        );
    }
    Ok(())
}

#[test]
fn undecided_branches_fold_published_fixings_and_keep_missing_ones() -> Result<()> {
    // 判定が評価時点より後の Cond の枝でも、公表済みの fixing と部分確定した複利を定数にする。
    // 価格計算は残った観測を予想するので、ここで定数にしないと過去の値を予想してしまう。
    // 欠損だけは枝に残し、判定でその枝を受け取ると決まった時点で報告する。
    let refdata = refdata()?;
    let spot: IndexId = "USDJPY-SPOT".parse()?;
    let tona: IndexId = "JPY-TONA".parse()?;
    let today = "2027-01-08T01:00:00Z".parse::<Instant>()?;
    let gate_at = "2027-01-11T01:00:00Z".parse::<Instant>()?;
    let payment = "2027-01-12T03:00:00Z".parse()?;
    let observations = observation_rows(
        &refdata,
        ObservationMethod::Lookback(NonZeroU8::MIN),
        Aggregation::Compound {
            form: CompoundingForm::SimpleFactor,
            spread: SpreadTreatment::Exclusive,
        },
        ObservationRounding::default(),
    )?;
    let thursday = ObservationPoint::Date("2027-01-07".parse()?);
    let friday = ObservationPoint::Date("2027-01-08".parse()?);
    let product = Product::new(
        Contract::build(|b| {
            let later = b.fixing(&spot, ObservationPoint::Instant(gate_at));
            let zero = b.constant(0.0);
            let gate = b.gt(later, zero);
            let now = b.fixing(&spot, ObservationPoint::Instant(today));
            let growth =
                b.compounded(&tona, &observations, observations.rule(), Spread::new(0.0)?)?;
            let yes = b.pay(CASH, now, payment);
            let no = b.pay(CASH, growth, payment);
            Ok(b.cond(gate, gate_at, yes, no))
        })?,
        &refdata,
    )?;
    let published = FixingStore::new(
        &refdata,
        FixingStoreParts {
            published: vec![
                (spot.clone(), ObservationPoint::Instant(today), 150.0),
                (tona.clone(), thursday, 0.1),
            ],
            ..FixingStoreParts::default()
        },
    )?;
    let reduced = reduce(
        &product,
        today,
        &published,
        &[],
        SettlementBoundary::Exclude,
    )?;
    let (yes, no) = cond_branches(&reduced)?;
    // 受け取るか未定の枝でも、公表済みの 150 円は定数になる。
    let ExprNode::Const(value) = pay_amount(&reduced, yes)? else {
        return Err(Error::Internal("expected the published spot".into()));
    };
    assert_close!(value.value(), 150.0, 1e-12, 1e-6);
    // 木曜分の因子だけを超過に畳み、期間中の金曜分は残りの複利として予想の開始時刻を記録する。
    let ExprNode::Compounded {
        rows, fixed_excess, ..
    } = pay_amount(&reduced, no)?
    else {
        return Err(Error::Internal(
            "expected residual compounding with the folded excess".into(),
        ));
    };
    assert_close!(fixed_excess.value(), 0.1 * 3.0 / 365.0, 1e-15, 0.0);
    assert_eq!(points(rows), [friday]);
    assert_eq!(reduced.forecast_from(&tona, friday), Some(today));

    // 値の無い過去の観測は、未定の枝の中ではエラーにせず観測のまま残す。
    let empty = FixingStore::default();
    let deferred = reduce(&product, today, &empty, &[], SettlementBoundary::Exclude)?;
    let (yes, no) = cond_branches(&deferred)?;
    assert!(matches!(
        pay_amount(&deferred, yes)?,
        ExprNode::Fixing { .. }
    ));
    let ExprNode::Compounded { rows, .. } = pay_amount(&deferred, no)? else {
        return Err(Error::Internal("expected the whole compounding".into()));
    };
    assert_eq!(points(rows), [thursday, friday]);

    // 判定で yes の枝を受け取ると決まると、その枝の欠損を報告する。
    let decided = empty.with_assumed_published(
        &refdata,
        [(spot, ObservationPoint::Instant(gate_at), 1.0)],
        [],
    )?;
    assert!(matches!(
        reduce(
            &product,
            gate_at,
            &decided,
            &[],
            SettlementBoundary::Exclude
        ),
        Err(Error::MissingData(_))
    ));
    Ok(())
}

#[test]
fn undecided_branches_defer_numerical_failures_of_published_aggregates() -> Result<()> {
    // 判定前の Cond の枝にある集約で、公表済みの値から作る確定部分が定義域の外になる場合、
    // 簡約は失敗せず、その集約を理由付きの Undefined にして値が必要になるまで残す。
    // 確定部分の失敗は残りの予想によらず決まるので、公表値を枝に残す必要はない。
    // 判定で枝を受け取らないと決まれば簡約は成功し、受け取ると決まれば Numerical を返す。
    let refdata = refdata()?;
    let spot: IndexId = "USDJPY-SPOT".parse()?;
    let tona: IndexId = "JPY-TONA".parse()?;
    let today = "2027-01-08T01:00:00Z".parse::<Instant>()?;
    let gate_at = "2027-01-11T01:00:00Z".parse::<Instant>()?;
    let payment = "2027-01-12T03:00:00Z".parse()?;
    let thursday = ObservationPoint::Date("2027-01-07".parse()?);
    // 複利: 契約の Inclusive のスプレッド −200（−20,000%）で、木曜の公表値 0.1% の因子
    // 1 + (0.001 − 200)·3/365 が負になる。
    let compounding = observation_rows(
        &refdata,
        ObservationMethod::Lookback(NonZeroU8::MIN),
        Aggregation::Compound {
            form: CompoundingForm::SimpleFactor,
            spread: SpreadTreatment::Inclusive,
        },
        ObservationRounding::default(),
    )?;
    // 平均: 日次の丸め 1e-6 で、木曜の公表値 1e12 は商が 2^53 を超えて丸められない。
    let averaging = observation_rows(
        &refdata,
        ObservationMethod::Lookback(NonZeroU8::MIN),
        Aggregation::Average(AveragingWeights::CoverWeighted),
        ObservationRounding {
            daily: Some(Rounding::new(1e-6, RoundingMode::HalfEven)?),
            aggregated: None,
        },
    )?;
    for (name, value) in [("compounded", 0.001), ("average", 1e12)] {
        let product = Product::new(
            Contract::build(|b| {
                let later = b.fixing(&spot, ObservationPoint::Instant(gate_at));
                let zero = b.constant(0.0);
                let gate = b.gt(later, zero);
                let aggregate = if name == "compounded" {
                    b.compounded(
                        &tona,
                        &compounding,
                        compounding.rule(),
                        Spread::new(-200.0)?,
                    )?
                } else {
                    b.average(&tona, &averaging)?
                };
                let one = b.constant(1.0);
                let yes = b.pay(CASH, one, payment);
                let no = b.pay(CASH, aggregate, payment);
                Ok(b.cond(gate, gate_at, yes, no))
            })?,
            &refdata,
        )?;
        let published = FixingStore::new(
            &refdata,
            FixingStoreParts {
                published: vec![(tona.clone(), thursday, value)],
                ..FixingStoreParts::default()
            },
        )?;
        let reduced = reduce(
            &product,
            today,
            &published,
            &[],
            SettlementBoundary::Exclude,
        )?;
        let (_, no) = cond_branches(&reduced)?;
        let ExprNode::Undefined(reason) = pay_amount(&reduced, no)? else {
            return Err(Error::Internal(format!(
                "{name}: expected a deferred failure"
            )));
        };
        assert!(!reason.is_empty());

        // 判定の値で受け取る枝を決め、判定時刻の簡約で結果を確かめる。
        for (gate, should_fail) in [(1.0, false), (-1.0, true)] {
            let decided = published.with_assumed_published(
                &refdata,
                [(spot.clone(), ObservationPoint::Instant(gate_at), gate)],
                [],
            )?;
            let result = reduce(
                &product,
                gate_at,
                &decided,
                &[],
                SettlementBoundary::Exclude,
            );
            if should_fail {
                assert!(
                    matches!(result, Err(Error::Numerical(_))),
                    "{name}: {result:?}"
                );
            } else {
                assert_close!(const_amount(&result?)?, 1.0, 1e-12, 1e-6);
            }
        }
    }
    Ok(())
}

#[test]
fn pending_condition_requires_its_published_fixings() -> Result<()> {
    // 到達した Cond の条件は判定時刻に必ず評価するので、判定が評価時点より後でも、
    // 公表済みの値が欠ければエラーにする。分かる項は畳み、まだ公表されない項だけを残す。
    let refdata = refdata()?;
    let spot: IndexId = "USDJPY-SPOT".parse()?;
    let today = "2027-01-08T01:00:00Z".parse::<Instant>()?;
    let gate_at = "2027-01-11T01:00:00Z".parse::<Instant>()?;
    let payment = "2027-01-12T03:00:00Z".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let now = b.fixing(&spot, ObservationPoint::Instant(today));
            let later = b.fixing(&spot, ObservationPoint::Instant(gate_at));
            let strike = b.constant(100.0);
            let zero = b.constant(0.0);
            let above = b.gt(now, strike);
            let positive = b.gt(later, zero);
            let gate = b.and([above, positive]);
            let one = b.constant(1.0);
            let pay = b.pay(CASH, one, payment);
            let nothing = b.zero();
            Ok(b.cond(gate, gate_at, pay, nothing))
        })?,
        &refdata,
    )?;
    let empty = FixingStore::default();
    assert!(matches!(
        reduce(&product, today, &empty, &[], SettlementBoundary::Exclude),
        Err(Error::MissingData(_))
    ));
    for (now, pending) in [(150.0, true), (50.0, false)] {
        let store = empty.with_assumed_published(
            &refdata,
            [(spot.clone(), ObservationPoint::Instant(today), now)],
            [],
        )?;
        let reduced = reduce(&product, today, &store, &[], SettlementBoundary::Exclude)?;
        let ContractNode::Cond { cond, .. } = root(&reduced) else {
            return Err(Error::Internal("expected a pending condition".into()));
        };
        match &reduced.contract.exprs[*cond as usize] {
            // 150 > 100 は真なので、論理積はまだ公表されない比較だけになる。
            ExprNode::Binary(nodes::Binary::Gt, [later, _]) => {
                assert!(pending);
                assert!(matches!(
                    reduced.contract.exprs[*later as usize],
                    ExprNode::Fixing {
                        point: ObservationPoint::Instant(at),
                        ..
                    } if at == gate_at
                ));
            }
            // 50 > 100 は偽なので論理積は偽に決まる。判定時刻は先なので Cond は残す。
            ExprNode::Bool(false) => assert!(!pending),
            other => {
                return Err(Error::Internal(format!("unexpected condition {other:?}")));
            }
        }
    }
    Ok(())
}

#[test]
fn effective_exercise_inside_an_undecided_branch_selects_its_option() -> Result<()> {
    // 外側の Cond が未定でも、効力のある行使記録は内側の Choose の結果を決める。
    // 記録を無視すると、価格計算は既に決まった枝ではなく有利な枝の価値を返してしまう。
    let refdata = refdata()?;
    let spot: IndexId = "USDJPY-SPOT".parse()?;
    let today = "2027-01-08T01:00:00Z".parse::<Instant>()?;
    let gate_at = "2027-01-08T05:00:00Z".parse::<Instant>()?;
    let deadline = "2027-01-11T01:00:00Z".parse::<Instant>()?;
    let payment = "2027-01-12T03:00:00Z".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let observed = b.fixing(&spot, ObservationPoint::Instant(gate_at));
            let zero = b.constant(0.0);
            let gate = b.gt(observed, zero);
            let low = b.constant(100.0);
            let high = b.constant(110.0);
            let first = b.pay(CASH, low, payment);
            let second = b.pay(CASH, high, payment);
            let call = b.choose(Party::Us, "call", deadline, [first, second]);
            let nothing = b.zero();
            Ok(b.cond(gate, gate_at, call, nothing))
        })?,
        &refdata,
    )?;
    let empty = FixingStore::default();
    let record = ExerciseRecord {
        label: "call".into(),
        at: deadline,
        effective_at: today,
        option: 1,
    };
    let exercised = reduce(
        &product,
        today,
        &empty,
        std::slice::from_ref(&record),
        SettlementBoundary::Exclude,
    )?;
    let (yes, _) = cond_branches(&exercised)?;
    let ExprNode::Const(value) = pay_amount(&exercised, yes)? else {
        return Err(Error::Internal("expected the exercised payment".into()));
    };
    assert_close!(value.value(), 110.0, 1e-12, 1e-6);
    // 記録が無ければ、判断は締切まで選択肢を残す。
    let open = reduce(&product, today, &empty, &[], SettlementBoundary::Exclude)?;
    let (yes, _) = cond_branches(&open)?;
    assert!(matches!(
        open.contract.contracts[yes as usize],
        ContractNode::Choose { .. }
    ));
    Ok(())
}

#[test]
fn deep_undecided_selects_are_reduced_without_recursion() -> Result<()> {
    // 条件の決まらない Select の腕は Deferred として同じ作業スタックで簡約する。
    // 1 万 2 千段の入れ子でもスタックの深さに依存せず、共有した部分式を木に展開しない。
    const DEPTH: usize = 12_000;
    let refdata = refdata()?;
    let spot: IndexId = "USDJPY-SPOT".parse()?;
    let today = "2027-01-08T01:00:00Z".parse()?;
    let gate_at = "2027-01-11T01:00:00Z".parse()?;
    let payment = "2027-01-12T03:00:00Z".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let later = b.fixing(&spot, ObservationPoint::Instant(gate_at));
            let zero = b.constant(0.0);
            let undecided = b.gt(later, zero);
            let step = b.constant(0.25);
            let mut amount = step;
            for _ in 0..DEPTH {
                let next = b.add([step, amount]);
                amount = b.select(undecided, next, amount);
            }
            Ok(b.pay(CASH, amount, payment))
        })?,
        &refdata,
    )?;
    let reduced = reduce(
        &product,
        today,
        &FixingStore::default(),
        &[],
        SettlementBoundary::Exclude,
    )?;
    assert!(matches!(amount_expr(&reduced)?, ExprNode::Select(_)));
    assert!(reduced.contract.exprs.len() < 3 * DEPTH);
    Ok(())
}

#[test]
fn recovery_is_required_only_on_the_default_branch() -> Result<()> {
    let refdata = refdata()?;
    let event = "JP-EVENT".parse()?;
    let recovery = "JP-RECOVERY".parse()?;
    let window = DateInterval::new("2027-01-07".parse()?, "2027-01-08".parse()?)?;
    let payment = "2027-01-12T03:00:00Z".parse()?;
    let as_of = "2027-01-11T05:00:00Z".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let event = b.defaulted(&event, window);
            let recovery = b.recovery(&recovery, window);
            let one = b.constant(1.0);
            let zero = b.constant(0.0);
            let loss = b.sub(one, recovery);
            let amount = b.select(event, loss, zero);
            Ok(b.pay(CASH, amount, payment))
        })?,
        &refdata,
    )?;
    let empty = FixingStore::default();
    assert_close!(
        const_amount(&reduce(
            &product,
            as_of,
            &empty,
            &[],
            SettlementBoundary::Exclude
        )?)?,
        0.0,
        1e-12,
        1e-14
    );
    let events = [(
        credit_key()?,
        CreditEventRecord {
            event_date: window.end(),
            kind: "failure_to_pay".into(),
        },
    )];
    let incomplete = empty.with_assumed_published(&refdata, [], events.clone())?;
    assert!(matches!(
        reduce(
            &product,
            as_of,
            &incomplete,
            &[],
            SettlementBoundary::Exclude
        ),
        Err(Error::MissingData(_))
    ));
    let resolved = FixingStore::new(
        &refdata,
        FixingStoreParts {
            credit_events: events.into(),
            recoveries: vec![(credit_key()?, RecoveryRate::new(0.4)?)],
            ..FixingStoreParts::default()
        },
    )?;
    let reduced = reduce(&product, as_of, &resolved, &[], SettlementBoundary::Exclude)?;
    assert_close!(const_amount(&reduced)?, 0.6, 1e-12, 1e-14);
    Ok(())
}

#[test]
fn deep_shared_dag_is_validated_and_reduced_without_recursion() -> Result<()> {
    let refdata = refdata()?;
    let at = "2027-01-12T03:00:00Z".parse()?;
    let as_of = "2027-01-11T05:00:00Z".parse()?;
    let contract = Contract::build(|b| {
        let step = b.constant(0.25);
        let mut amount = step;
        for _ in 0..12_000 {
            amount = b.add([step, amount]);
        }
        let mut root = b.pay(CASH, amount, at);
        for _ in 0..12_000 {
            root = b.all([root]);
        }
        Ok(root)
    })?;
    let product = Product::new(contract, &refdata)?;
    let reduced = reduce(
        &product,
        as_of,
        &FixingStore::default(),
        &[],
        SettlementBoundary::Exclude,
    )?;
    assert_close!(const_amount(&reduced)?, 12_001.0 * 0.25, 1e-12, 1e-6);
    assert_eq!(reduced.contract.contracts.len(), 1);
    assert_eq!(reduced.contract.exprs.len(), 1);
    Ok(())
}

#[test]
fn compound_factors_apply_daily_rounding_and_the_selected_spread_rule() -> Result<()> {
    let refdata = refdata()?;
    let id: IndexId = "JPY-TONA".parse()?;
    let payment = "2027-01-12T03:00:00Z".parse()?;
    let as_of = "2027-01-11T01:00:00Z".parse()?;
    let fixings = FixingStore::new(
        &refdata,
        FixingStoreParts {
            published: vec![
                (
                    id.clone(),
                    ObservationPoint::Date("2027-01-07".parse()?),
                    0.104,
                ),
                (
                    id.clone(),
                    ObservationPoint::Date("2027-01-08".parse()?),
                    0.206,
                ),
            ],
            ..FixingStoreParts::default()
        },
    )?;
    // 丸めた 0.10 と 0.21 を、単利では日次スプレッド込み、冪ではスプレッドなしで使う。
    for (form, spread, expected) in [
        (
            CompoundingForm::SimpleFactor,
            SpreadTreatment::Inclusive,
            (1.0 + 0.12 * 3.0 / 365.0) * (1.0 + 0.23 / 365.0),
        ),
        (
            CompoundingForm::PowerFactor,
            SpreadTreatment::Exclusive,
            1.1_f64.powf(3.0 / 365.0) * 1.21_f64.powf(1.0 / 365.0),
        ),
    ] {
        let observations = observation_rows(
            &refdata,
            ObservationMethod::Lookback(NonZeroU8::MIN),
            Aggregation::Compound { form, spread },
            ObservationRounding {
                daily: Some(Rounding::new(0.01, RoundingMode::HalfEven)?),
                aggregated: None,
            },
        )?;
        let product = Product::new(
            Contract::build(|b| {
                let gain =
                    b.compounded(&id, &observations, observations.rule(), Spread::new(0.02)?)?;
                Ok(b.pay(CASH, gain, payment))
            })?,
            &refdata,
        )?;
        let reduced = reduce(&product, as_of, &fixings, &[], SettlementBoundary::Exclude)?;
        assert_close!(const_amount(&reduced)?, expected, 1e-12, 1e-14);
    }
    Ok(())
}

#[test]
fn product_rejects_distinct_choices_with_the_same_key_and_empty_choices() -> Result<()> {
    let refdata = refdata()?;
    let gate = "2027-01-08T01:00:00Z".parse()?;
    let at = "2027-01-11T01:00:00Z".parse()?;
    let duplicate = Contract::build(|b| {
        let x = b.constant(1.0);
        let pay = b.pay(CASH, x, at);
        let zero = b.zero();
        let first = b.choose(Party::Us, "exercise", gate, [pay, zero]);
        let second = b.choose(Party::Us, "exercise", gate, [zero, pay]);
        Ok(b.all([first, second]))
    })?;
    // 同じ行使記録が異なる選択に一致する構成と、選べる枝のない構成を拒否する。
    assert!(matches!(
        Product::new(duplicate, &refdata),
        Err(Error::InvalidInput(_))
    ));
    let empty = Contract::build(|b| Ok(b.choose(Party::Us, "exercise", gate, [])))?;
    assert!(matches!(
        Product::new(empty, &refdata),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn give_flips_signs_and_parties_by_parity() -> Result<()> {
    // 旧 IR の give の例を移す。give の偶奇は額の符号と選択の当事者だけに作用し、
    // 符号付きゼロも符号を反転した値として保つ（0 は −0、−0 は +0 になる）。
    let refdata = refdata()?;
    let at = "2027-01-12T03:00:00Z".parse()?;
    let deadline = "2027-01-11T06:00:00Z".parse()?;
    let as_of = "2027-01-11T05:00:00Z".parse()?;
    let empty = FixingStore::default();
    for (value, gifts, expected) in [
        (100.0, 1, -100.0),
        (100.0, 2, 100.0),
        (0.0, 1, -0.0),
        (-0.0, 1, 0.0),
    ] {
        let product = Product::new(
            Contract::build(|b| {
                let amount = b.constant(value);
                let pay = b.pay(CASH, amount, at);
                Ok((0..gifts).fold(pay, |contract, _| b.give(contract)))
            })?,
            &refdata,
        )?;
        let amount = const_amount(&reduce(
            &product,
            as_of,
            &empty,
            &[],
            SettlementBoundary::Exclude,
        )?)?;
        assert_close!(amount, expected, 1e-12, 1e-6);
        assert_eq!(amount.is_sign_negative(), expected.is_sign_negative());
    }
    // 受渡のない契約は give しても受渡がない。
    let nothing = Product::new(
        Contract::build(|b| {
            let zero = b.zero();
            Ok(b.give(zero))
        })?,
        &refdata,
    )?;
    let reduced = reduce(&nothing, as_of, &empty, &[], SettlementBoundary::Exclude)?;
    assert!(matches!(root(&reduced), ContractNode::All(xs) if xs.is_empty()));
    // 選択の当事者も同じ偶奇で入れ替わる。
    for (gifts, party) in [(0, Party::Us), (1, Party::Counterparty), (2, Party::Us)] {
        let product = Product::new(
            Contract::build(|b| {
                let amount = b.constant(1.0);
                let pay = b.pay(CASH, amount, at);
                let zero = b.zero();
                let choice = b.choose(Party::Us, "exercise", deadline, [pay, zero]);
                Ok((0..gifts).fold(choice, |contract, _| b.give(contract)))
            })?,
            &refdata,
        )?;
        let reduced = reduce(&product, as_of, &empty, &[], SettlementBoundary::Exclude)?;
        assert!(matches!(root(&reduced), ContractNode::Choose { by, .. } if *by == party));
    }
    Ok(())
}

#[test]
fn physical_delivery_keeps_the_asset_and_each_shared_quantity() -> Result<()> {
    // 旧 IR の現物受渡の例を移す。同じ資産・時刻の受渡を二度持てば二度受け取り、
    // 係数倍は額に掛かり、give した 0 単位は −0 として残る。資産と受渡の瞬間は変わらない。
    let refdata = refdata()?;
    let oil = Asset::Commodity("WTI-CUSHING".parse()?);
    let first_at = "2027-01-12T03:00:00Z".parse::<Instant>()?;
    let later_at = "2027-01-13T03:00:00Z".parse::<Instant>()?;
    let as_of = "2027-01-11T05:00:00Z".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let factor = b.constant(3.0);
            let lot = b.constant(2.0);
            let none = b.constant(0.0);
            let quantity = b.mul([lot, factor]);
            let empty = b.mul([none, factor]);
            let first = b.pay(oil.clone(), quantity, first_at);
            let second = b.pay(oil.clone(), quantity, first_at);
            let later = b.pay(oil.clone(), empty, later_at);
            let returned = b.give(later);
            Ok(b.all([first, second, returned]))
        })?,
        &refdata,
    )?;
    let reduced = reduce(
        &product,
        as_of,
        &FixingStore::default(),
        &[],
        SettlementBoundary::Exclude,
    )?;
    let ContractNode::All(xs) = root(&reduced) else {
        return Err(Error::Internal("expected three deliveries".into()));
    };
    // 同じ受渡は一つのノードを共有したまま二度数える。
    assert_eq!(xs.len(), 3);
    assert_eq!(xs[0], xs[1]);
    for (&ix, (expected_at, expected)) in
        xs.iter()
            .zip([(first_at, 6.0), (first_at, 6.0), (later_at, -0.0)])
    {
        let ContractNode::Pay { asset, at, .. } = &reduced.contract.contracts[ix as usize] else {
            return Err(Error::Internal("expected a delivery".into()));
        };
        let ExprNode::Const(quantity) = pay_amount(&reduced, ix)? else {
            return Err(Error::Internal("expected a fixed quantity".into()));
        };
        assert_eq!(asset, &oil);
        assert_eq!(*at, expected_at);
        assert_close!(quantity.value(), expected, 1e-12, 1e-12);
        assert_eq!(
            quantity.value().is_sign_negative(),
            expected.is_sign_negative()
        );
    }
    // 参照データに無い資産の受渡は検査で拒否する。
    let unknown = Contract::build(|b| {
        let one = b.constant(1.0);
        Ok(b.pay(Asset::Commodity("BRENT".parse()?), one, first_at))
    })?;
    assert!(matches!(
        Product::new(unknown, &refdata),
        Err(Error::InvalidData(_))
    ));
    Ok(())
}

#[test]
fn contract_rounding_folds_half_even() -> Result<()> {
    // 旧 IR の丸めの例を移す。契約が定める丸めは定数の式として畳み、
    // 0.01 単位の HalfEven は中点を偶数側に寄せる（0.125 → 0.12、0.135 → 0.14）。
    let refdata = refdata()?;
    let at = "2027-01-12T03:00:00Z".parse()?;
    let as_of = "2027-01-11T05:00:00Z".parse()?;
    let rounding = Rounding::new(0.01, RoundingMode::HalfEven)?;
    for (value, expected) in [(0.125, 0.12), (-0.125, -0.12), (0.135, 0.14)] {
        let product = Product::new(
            Contract::build(|b| {
                let raw = b.constant(value);
                let rounded = b.round(raw, rounding);
                Ok(b.pay(CASH, rounded, at))
            })?,
            &refdata,
        )?;
        let reduced = reduce(
            &product,
            as_of,
            &FixingStore::default(),
            &[],
            SettlementBoundary::Exclude,
        )?;
        assert_close!(const_amount(&reduced)?, expected, 1e-12, 1e-14);
    }
    Ok(())
}

#[test]
fn node_budget_rejects_contracts_over_the_limit() -> Result<()> {
    // 旧 IR の予算の例を移す。商品の誤りによるノード数の爆発を、上限ちょうどまでは許し、
    // 一つでも超えたら builder が InvalidInput で止める。
    let at = "2027-01-12T03:00:00Z".parse()?;
    let build = |nodes: usize| {
        Contract::build(|b| {
            // 値の異なる定数は共有されないので、定数 nodes − 1 個と支払い 1 個でちょうど nodes 個になる。
            let amount = (1..nodes - 1).fold(b.constant(0.0), |_, value| b.constant(value as f64));
            Ok(b.pay(CASH, amount, at))
        })
    };
    assert!(build(MAX_NODES).is_ok());
    assert!(matches!(build(MAX_NODES + 1), Err(Error::InvalidInput(_))));
    Ok(())
}

#[test]
fn observation_builders_reject_rules_that_did_not_make_the_schedule() -> Result<()> {
    // 旧 IR の観測列の検査の例のうち、新しい builder に当たるものを移す。
    // 複利は観測列を作った複利の規則とだけ、平均は平均の規則で作った観測列とだけ組める。
    let refdata = refdata()?;
    let id: IndexId = "JPY-TONA".parse()?;
    let at = "2027-01-12T03:00:00Z".parse()?;
    let compound = Aggregation::Compound {
        form: CompoundingForm::SimpleFactor,
        spread: SpreadTreatment::Exclusive,
    };
    let lookback = observation_rows(
        &refdata,
        ObservationMethod::Lookback(NonZeroU8::MIN),
        compound,
        ObservationRounding::default(),
    )?;
    let plain = observation_rows(
        &refdata,
        ObservationMethod::Plain,
        compound,
        ObservationRounding::default(),
    )?;
    let average = observation_rows(
        &refdata,
        ObservationMethod::Lookback(NonZeroU8::MIN),
        Aggregation::Average(AveragingWeights::Equal),
        ObservationRounding::default(),
    )?;
    for (observations, rule) in [(&lookback, plain.rule()), (&average, average.rule())] {
        let result = Contract::build(|b| {
            let gain = b.compounded(&id, observations, rule, Spread::new(0.0)?)?;
            Ok(b.pay(CASH, gain, at))
        });
        assert!(matches!(result, Err(Error::InvalidInput(_))));
    }
    let result = Contract::build(|b| {
        let mean = b.average(&id, &lookback)?;
        Ok(b.pay(CASH, mean, at))
    });
    assert!(matches!(result, Err(Error::InvalidInput(_))));
    Ok(())
}

mod lifecycle;
mod standard;

#[test]
fn analysis_preserves_payment_occurrences_and_branch_paths() -> Result<()> {
    // 同じ支払ノードが繰り返し現れても支払回数を保つ。式の観測だけを重複排除する。
    use analysis::{Branch, Linearity, Observation, cashflows, linearity, observations};

    let data = refdata()?;
    let spot: IndexId = "USDJPY-SPOT".parse()?;
    let observed = "2027-02-01T00:00:00Z".parse()?;
    let decided = "2027-02-02T00:00:00Z".parse()?;
    let paid = "2027-02-03T00:00:00Z".parse()?;
    let point = ObservationPoint::Instant(observed);
    let contract = Contract::build(|builder| {
        let fixing = builder.fixing(&spot, point);
        let one = builder.constant(1.0);
        let pay = builder.pay(CASH, one, paid);
        let conditional_amount = builder.pay(CASH, fixing, paid);
        let choice = builder.choose(Party::Us, "exercise", decided, [conditional_amount, pay]);
        let condition = builder.gt(fixing, one);
        let contingent = builder.cond(condition, decided, choice, pay);
        Ok(builder.all([pay, pay, contingent]))
    })?;
    let product = Product::new(contract, &data)?;
    let flows = cashflows(&product);
    assert_eq!(flows.len(), 5);
    assert_eq!(flows[0].branches.len(), 0);
    assert_eq!(flows[1].amount.index(), flows[0].amount.index());
    assert!(matches!(
        flows[2].branches.as_slice(),
        [
            Branch::Condition { when: true, .. },
            Branch::Choice {
                label: "exercise",
                option: 0,
                ..
            }
        ]
    ));
    assert!(matches!(
        flows[4].branches.as_slice(),
        [Branch::Condition { when: false, .. }]
    ));
    assert_eq!(
        observations(&product),
        [Observation::Fixing {
            index: &spot,
            point
        }]
    );
    assert!(matches!(linearity(&product), Linearity::NonLinear(_)));
    Ok(())
}

/// 予想する複利 g と定数 2 から支払額を組み、リニア方式への適用可否を判定する。
fn linearity_of(
    refdata: &RefData,
    observations: &ObservationSchedule,
    amount: impl for<'b> FnOnce(&mut Builder<'b>, ExprIx<'b>, ExprIx<'b>) -> ExprIx<'b>,
) -> Result<analysis::Linearity> {
    let id: IndexId = "JPY-TONA".parse()?;
    let at = "2027-01-12T03:00:00Z".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let growth = b.compounded(&id, observations, observations.rule(), Spread::new(0.0)?)?;
            let constant = b.constant(2.0);
            let amount = amount(b, growth, constant);
            Ok(b.pay(CASH, amount, at))
        })?,
        refdata,
    )?;
    Ok(analysis::linearity(&product))
}

#[test]
fn linearity_rejects_products_and_quotients_of_projected_growth() -> Result<()> {
    // 予想する複利 g の定数倍、定数での割り算、丸め、符号の反転と差は、g について線形なので
    // リニア評価で扱える。g 同士の積、g による割り算、g を含む最大値は線形でないため、
    // 方式の選択がその理由とともに Unsupported を返せるよう NonLinear にする。
    use analysis::Linearity;

    let refdata = refdata()?;
    let observations = observation_rows(
        &refdata,
        ObservationMethod::Lookback(NonZeroU8::MIN),
        Aggregation::Compound {
            form: CompoundingForm::SimpleFactor,
            spread: SpreadTreatment::Exclusive,
        },
        ObservationRounding::default(),
    )?;
    let rounding = Rounding::new(0.0001, RoundingMode::HalfEven)?;

    let scaled = linearity_of(&refdata, &observations, |b, g, c| {
        let notional = b.constant(1_000.0);
        let term = b.mul([g, notional]);
        b.add([term, c])
    })?;
    assert_eq!(scaled, Linearity::Linear);
    let divided = linearity_of(&refdata, &observations, |b, g, _| {
        let denominator = b.constant(4.0);
        b.div(g, denominator)
    })?;
    assert_eq!(divided, Linearity::Linear);
    let rounded = linearity_of(&refdata, &observations, |b, g, _| b.round(g, rounding))?;
    assert_eq!(rounded, Linearity::Linear);
    let negated = linearity_of(&refdata, &observations, |b, g, c| {
        let negative = b.neg(g);
        b.sub(negative, c)
    })?;
    assert_eq!(negated, Linearity::Linear);

    let squared = linearity_of(&refdata, &observations, |b, g, _| b.mul([g, g]))?;
    assert_eq!(
        squared,
        Linearity::NonLinear("multiplication of two projected quantities")
    );
    let reciprocal = linearity_of(&refdata, &observations, |b, g, c| b.div(c, g))?;
    assert_eq!(
        reciprocal,
        Linearity::NonLinear("division by a projected quantity")
    );
    let floored = linearity_of(&refdata, &observations, |b, g, c| b.max([g, c]))?;
    assert_eq!(
        floored,
        Linearity::NonLinear("unsupported cashflow expression")
    );
    Ok(())
}

#[test]
fn cond_rejects_a_branch_paying_before_its_decision() -> Result<()> {
    // Cond は判断の瞬間より前に受け渡す枝を持てない（Choose と同じ因果性の規則）。
    let refdata = refdata()?;
    let id: IndexId = "USDJPY-SPOT".parse()?;
    let publication = "2027-01-08T01:00:00Z".parse::<Instant>()?;
    let before = Instant::from_unix_micros(publication.unix_micros() - 1)?;
    for (paid, accepted) in [(before, false), (publication, true)] {
        let contract = Contract::build(|b| {
            let spot = b.fixing(&id, ObservationPoint::Instant(publication));
            let strike = b.constant(100.0);
            let condition = b.gt(spot, strike);
            let pay = b.pay(CASH, strike, paid);
            let zero = b.zero();
            Ok(b.cond(condition, publication, pay, zero))
        })?;
        let product = Product::new(contract, &refdata);
        if accepted {
            // 判断と同じ瞬間の受渡は許す。
            product?;
        } else {
            assert!(
                matches!(&product, Err(Error::InvalidInput(message)) if message.contains("precedes decision")),
                "{product:?}"
            );
        }
    }
    Ok(())
}

// Product::new が契約を拒否し、エラーの種類が期待どおりであることを確かめる。
fn assert_rejected(
    refdata: &RefData,
    what: &str,
    expected: fn(&Error) -> bool,
    contract: Contract,
) {
    let result = Product::new(contract, refdata);
    assert!(
        matches!(&result, Err(error) if expected(error)),
        "{what}: {result:?}"
    );
}

fn is_input(error: &Error) -> bool {
    matches!(error, Error::InvalidInput(_))
}

#[test]
fn product_new_rejects_malformed_expressions() -> Result<()> {
    // 非有限の定数と、項の無い min・max は builder で組めても検査で拒否する。
    let refdata = refdata()?;
    let at = "2027-01-12T03:00:00Z".parse()?;
    assert_rejected(
        &refdata,
        "non-finite constant",
        is_input,
        Contract::build(|b| {
            let nan = b.constant(f64::NAN);
            Ok(b.pay(CASH, nan, at))
        })?,
    );
    assert_rejected(
        &refdata,
        "empty min",
        is_input,
        Contract::build(|b| {
            let min = b.min([]);
            Ok(b.pay(CASH, min, at))
        })?,
    );
    assert_rejected(
        &refdata,
        "empty max",
        is_input,
        Contract::build(|b| {
            let max = b.max([]);
            Ok(b.pay(CASH, max, at))
        })?,
    );
    Ok(())
}

#[test]
fn product_new_rejects_index_misuse() -> Result<()> {
    // 指標の種類に合わない読み方を、参照データと照合して一度だけ拒否する。
    let refdata = refdata()?;
    let at = "2027-01-12T03:00:00Z".parse()?;
    let event: IndexId = "JP-EVENT".parse()?;
    let recovery: IndexId = "JP-RECOVERY".parse()?;
    let day = ObservationPoint::Date("2027-01-08".parse()?);
    let window = DateInterval::new("2027-01-08".parse()?, "2027-01-11".parse()?)?;
    let empty = DateInterval::new("2027-01-08".parse()?, "2027-01-08".parse()?)?;
    let rows = observation_rows(
        &refdata,
        ObservationMethod::Plain,
        Aggregation::Compound {
            form: CompoundingForm::SimpleFactor,
            spread: SpreadTreatment::Exclusive,
        },
        ObservationRounding::default(),
    )?;
    // 信用事由の指標は数の fixing として読めず、信用窓が要る。
    assert_rejected(
        &refdata,
        "credit index as a numeric fixing",
        is_input,
        Contract::build(|b| {
            let value = b.fixing(&event, day);
            Ok(b.pay(CASH, value, at))
        })?,
    );
    // 複利は翌日物の指標だけに使える。FX の fixing は数でも複利にできない。
    let compounded = |index: &IndexId| {
        Contract::build(|b| {
            let growth = b.compounded(index, &rows, rows.rule(), Spread::new(0.0)?)?;
            Ok(b.pay(CASH, growth, at))
        })
    };
    assert_rejected(
        &refdata,
        "compounding a non-overnight index",
        is_input,
        compounded(&"USDJPY-SPOT".parse()?)?,
    );
    Product::new(compounded(&"JPY-TONA".parse()?)?, &refdata)?;
    // 事由の判定と回収率は、それぞれ対応する種類の指標と空でない窓で読む。
    let defaulted = |index: &IndexId, window| {
        Contract::build(|b| {
            let defaulted = b.defaulted(index, window);
            let one = b.constant(1.0);
            let zero = b.constant(0.0);
            let amount = b.select(defaulted, one, zero);
            Ok(b.pay(CASH, amount, at))
        })
    };
    assert_rejected(
        &refdata,
        "defaulted on a recovery index",
        is_input,
        defaulted(&recovery, window)?,
    );
    assert_rejected(
        &refdata,
        "empty credit window",
        is_input,
        defaulted(&event, empty)?,
    );
    assert_rejected(
        &refdata,
        "recovery on an event index",
        is_input,
        Contract::build(|b| {
            let rate = b.recovery(&event, window);
            Ok(b.pay(CASH, rate, at))
        })?,
    );
    // 参照データに無い指標は、データの作り手が直す参照の欠け。
    assert_rejected(
        &refdata,
        "unknown index",
        |error| matches!(error, Error::InvalidData(_)),
        Contract::build(|b| {
            let value = b.fixing(&"EUR-ESTR".parse()?, day);
            Ok(b.pay(CASH, value, at))
        })?,
    );
    Ok(())
}

#[test]
fn tona_fixing_switches_at_publication_and_forecasts_a_spanning_period() -> Result<()> {
    // 金曜 2027-01-08 の TONA。観測期間は金曜 09:00 から月曜 09:00（JST）まで、
    // 公表は 1 営業日後の月曜 10:00 JST（01:00 UTC）。
    let refdata = refdata()?;
    let id: IndexId = "JPY-TONA".parse()?;
    let friday = ObservationPoint::Date("2027-01-08".parse()?);
    let paid = "2027-01-12T03:00:00Z".parse()?;
    let product = Product::new(
        Contract::build(|b| {
            let obs = b.fixing(&id, friday);
            Ok(b.pay(CASH, obs, paid))
        })?,
        &refdata,
    )?;
    let published = FixingStore::new(
        &refdata,
        FixingStoreParts {
            published: vec![(id.clone(), friday, 0.0046)],
            ..FixingStoreParts::default()
        },
    )?;
    let publication = "2027-01-11T01:00:00Z".parse::<Instant>()?;
    let before = Instant::from_unix_micros(publication.unix_micros() - 1)?;
    let after = Instant::from_unix_micros(publication.unix_micros() + 1)?;

    // 公表の 1µs 前は、ストアに値があってもまだ使えない。観測期間は終わっているので推定値が要る。
    assert!(matches!(
        reduce(
            &product,
            before,
            &published,
            &[],
            SettlementBoundary::Exclude
        ),
        Err(Error::MissingData(_))
    ));
    // 公表ちょうどと 1µs 後は公表値の定数になり、予想の開始時刻も警告も残らない。
    for as_of in [publication, after] {
        let fixed = reduce(
            &product,
            as_of,
            &published,
            &[],
            SettlementBoundary::Exclude,
        )?;
        assert_close!(const_amount(&fixed)?, 0.0046, 1e-12, 1e-14);
        assert_eq!(fixed.forecast_from(&id, friday), None);
        assert_eq!(fixed.warnings(), &[]);
    }

    // 観測期間の途中（金曜 14:00 JST）で推定値が無ければ、評価時点から期間の終わりまでの
    // 市場予想を期間全体の率とみなし、予想の開始時刻と警告を残す。値は評価のときに読む。
    let spanning = "2027-01-08T05:00:00Z".parse()?;
    let forecast = reduce(
        &product,
        spanning,
        &FixingStore::default(),
        &[],
        SettlementBoundary::Exclude,
    )?;
    assert!(matches!(amount_expr(&forecast)?, ExprNode::Fixing { .. }));
    assert_eq!(forecast.forecast_from(&id, friday), Some(spanning));
    assert_eq!(
        forecast.warnings(),
        &[Warning::UnpublishedFixingForecast {
            index: id.clone(),
            point: friday
        }]
    );
    // 期間が始まる前（金曜 08:59 JST）は、期間全体を市場から予想し、警告は出ない。
    let ahead = "2027-01-07T23:59:00Z".parse()?;
    let ahead = reduce(
        &product,
        ahead,
        &FixingStore::default(),
        &[],
        SettlementBoundary::Exclude,
    )?;
    assert!(matches!(amount_expr(&ahead)?, ExprNode::Fixing { .. }));
    assert_eq!(ahead.forecast_from(&id, friday), None);
    assert_eq!(ahead.warnings(), &[]);
    Ok(())
}
`````

`````
crates/rhodium-products/src/ir/tests/lifecycle.rs
use std::collections::BTreeMap;

use rhodium_core::{
    calendar::{BusinessDayConvention, EndOfMonthRule},
    date::Date,
    money::Money,
    refdata::{FallbackApplicability, FallbackDefinition, FallbackReplacement, FallbackTrigger},
};

use super::*;
use crate::{
    ir::analysis::{Observation, cashflows, observations},
    standard::PayReceive,
    trade::{FallbackSelection, FallbackTerms, Fee, ProductSpec, Trade},
};

fn trade(contract: Contract) -> Result<Trade> {
    Ok(Trade {
        id: "T-001".parse()?,
        margin_set: "JPY-CASH".parse()?,
        book: "BOOK".parse()?,
        trade_date: "2027-01-04".parse()?,
        product: ProductSpec::Contract(contract),
        attributes: BTreeMap::new(),
        exercises: Vec::new(),
        fees: Vec::new(),
        fallback: FallbackTerms::default(),
    })
}

fn amounts(product: &Product) -> Result<Vec<f64>> {
    cashflows(product)
        .iter()
        .map(|payment| match payment.amount.expression() {
            analysis::Expression::Constant(value) => Ok(value),
            _ => Err(Error::Internal("expected a resolved amount".into())),
        })
        .collect()
}

#[test]
fn partial_terminations_scale_final_amounts_at_the_cutoff_and_keep_fees() -> Result<()> {
    // 変動クーポンに上限をかけた後で解約する。元本を縮めてから上限をかけると 50 になる。
    let from: Instant = "2027-01-08T08:00:00Z".parse()?;
    let before = Instant::from_unix_micros(from.unix_micros() - 1)?;
    let after = Instant::from_unix_micros(from.unix_micros() + 1)?;
    let data = refdata()?;
    let original = trade(Contract::build(|b| {
        let principal = b.constant(100.0);
        let limit = b.constant(60.0);
        let amount = b.min([principal, limit]);
        let payments = [before, from, after].map(|at| b.pay(CASH, amount, at));
        Ok(b.all(payments))
    })?)?;
    let terminated = original
        .partially_terminated(
            from,
            0.5,
            Some(Fee {
                amount: Money::new(Currency::JPY, 7.0),
                side: PayReceive::Pay,
                at: from,
            }),
        )?
        .partially_terminated(after, 0.5, None)?;
    let fixings = FixingStore::default();
    let simplified = reduce(
        &terminated.to_product(&data)?,
        before,
        &fixings,
        &[],
        SettlementBoundary::Include,
    )?;
    let actual = amounts(&simplified)?;
    assert_eq!(actual.len(), 4);
    for (a, e) in actual.into_iter().zip([60.0, 30.0, 15.0, -7.0]) {
        assert_close!(a, e, 1e-12, 1e-6);
    }
    let original = reduce(
        &original.to_product(&data)?,
        from,
        &fixings,
        &[],
        SettlementBoundary::Exclude,
    )?;
    assert_close!(const_amount(&original)?, 60.0, 1e-12, 1e-6);
    Ok(())
}

#[test]
fn full_termination_drops_fixings_needed_only_by_cancelled_payments() -> Result<()> {
    // 解約後の分岐では、公表済みの fixing も締切済みの行使も読まない。
    let data = refdata()?;
    let index = "USDJPY-SPOT".parse()?;
    let from = "2027-01-08T08:00:00Z".parse()?;
    let observed = "2027-01-07T00:00:00Z".parse()?;
    let decided = "2027-01-07T08:00:00Z".parse()?;
    let original = trade(Contract::build(|b| {
        let fixing = b.fixing(&index, ObservationPoint::Instant(observed));
        let payment = b.pay(CASH, fixing, from);
        Ok(b.choose(Party::Us, "choice", decided, [payment]))
    })?)?;
    let product = original.terminated(from, None).to_product(&data)?;
    let simplified = reduce(
        &product,
        from,
        &FixingStore::default(),
        &[],
        SettlementBoundary::Include,
    )?;
    assert!(cashflows(&simplified).is_empty());
    Ok(())
}

#[test]
fn replacement_applies_pre_cutoff_exercise_before_splitting() -> Result<()> {
    // 旧契約の 100 円か境界ちょうどの 110 円かを選んでから切る。新契約の境界前の 999 円は除く。
    let data = refdata()?;
    let at = "2027-01-07T08:00:00Z".parse()?;
    let early = "2027-01-08T07:00:00Z".parse()?;
    let from = "2027-01-08T08:00:00Z".parse()?;
    let old = trade(Contract::build(|b| {
        let amount = b.constant(100.0);
        let first = b.pay(CASH, amount, early);
        let amount = b.constant(110.0);
        let second = b.pay(CASH, amount, from);
        Ok(b.choose(Party::Us, "old", at, [first, second]))
    })?)?;
    let new = ProductSpec::Contract(Contract::build(|b| {
        let wrong = b.constant(999.0);
        let dropped = b.pay(CASH, wrong, early);
        let amount = b.constant(30.0);
        let received = b.pay(CASH, amount, from);
        Ok(b.all([dropped, received]))
    })?);
    let amended = old.amended(from, new);
    assert!(matches!(
        amended.to_product(&data),
        Err(Error::Unsupported(_))
    ));
    for (option, expected) in [(0, vec![100.0, 30.0]), (1, vec![30.0])] {
        let record = ExerciseRecord {
            label: "old".into(),
            at,
            effective_at: at,
            option,
        };
        let selected = amended.exercised(record.clone());
        let simplified = reduce(
            &selected.to_product(&data)?,
            at,
            &FixingStore::default(),
            &[],
            SettlementBoundary::Exclude,
        )?;
        let actual = amounts(&simplified)?;
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.into_iter().zip(expected) {
            assert_close!(a, e, 1e-12, 1e-6);
        }
        assert!(matches!(
            selected.exercised(record).to_product(&data),
            Err(Error::InvalidData(_))
        ));
    }
    for (label, effective_at, option) in [("unknown", at, 0), ("old", early, 0), ("old", at, 2)] {
        assert!(matches!(
            old.exercised(ExerciseRecord {
                label: label.into(),
                at,
                effective_at,
                option
            })
            .to_product(&data),
            Err(Error::InvalidData(_))
        ));
    }
    Ok(())
}

fn fallback_data(status: IndexStatus) -> Result<(RefData, IndexId, IndexId)> {
    let mut parts = refdata_parts()?;
    let term: IndexId = "JPY-TERM".parse()?;
    let overnight = "JPY-TONA".parse()?;
    let mut definition = parts.indices[&overnight].clone();
    definition.id = term.clone();
    definition.kind = IndexKind::TermRate {
        ccy: Currency::JPY,
        day_count: DayCount::Act365F,
        tenor: "2BD".parse()?,
        value_offset: "0BD".parse()?,
        bdc: BusinessDayConvention::Following,
        end_of_month: EndOfMonthRule::Ignore,
    };
    definition.status = status;
    let PublicationRule::Daily(publication) = &mut definition.publication else {
        return Err(Error::Internal("expected a daily test index".into()));
    };
    publication.business_days = 0;
    parts.indices.insert(term.clone(), definition);
    parts.fallbacks.insert(
        term.clone(),
        FallbackDefinition {
            index: term.clone(),
            trigger: FallbackTrigger::CessationOrNonRepresentative,
            replacement: FallbackReplacement {
                index: overnight.clone(),
                observation: ObservationRule::new(
                    ObservationDays::BusinessDays(CalendarSpec::Base("JPTO".parse()?)),
                    ObservationMethod::Plain,
                    Aggregation::Compound {
                        form: CompoundingForm::SimpleFactor,
                        spread: SpreadTreatment::Exclusive,
                    },
                    ObservationRounding::default(),
                )?,
                annualization: DayCount::Act365F,
                spread: "10bp".parse()?,
            },
            applies: vec![
                FallbackApplicability::Definitions {
                    version: "2021".into(),
                    traded_on_or_after: "2027-01-05".parse()?,
                },
                FallbackApplicability::Protocol {
                    protocol: "p".into(),
                },
                FallbackApplicability::Bilateral,
            ],
        },
