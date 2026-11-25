2026-11-06,3.75%
2026-11-09,3.75%
2026-11-10,3.75%
2026-11-12,3.75%
2026-11-13,3.75%
2026-11-16,3.75%
2026-11-17,3.75%
2026-11-18,3.75%
2026-11-19,3.75%
2026-11-20,3.75%
`````

`````
examples/m1-rates/fixings/USDJPY-TKY-1500.csv
date,value
2026-11-24,152.35
`````

`````
examples/m1-rates/market/definition.toml
version = 1
valuation_calendar = "JPTO"
valuation_zone = "Asia/Tokyo"

# TONA: 1 年までの 3 回の会合で区切った 4 区間の階段、1 年で接続する長期の節点、2026 年末の跳び。
[curves.JPY-TONA]
currency = "JPY"
roles = ["discount", "projection"]

[[curves.JPY-TONA.expr.sum]]
stitch.at = 2027-12-03T00:00:00Z
stitch.short.steps = { part = "meetings", meetings = [2026-12-18T00:00:00Z, 2027-01-22T00:00:00Z, 2027-03-18T00:00:00Z] }
stitch.long.nodes = { part = "long", pillars = [2028-12-04T00:00:00Z, 2029-12-03T00:00:00Z, 2031-12-03T00:00:00Z, 2033-12-05T00:00:00Z, 2036-12-03T00:00:00Z] }

[[curves.JPY-TONA.expr.sum]]
jumps = { part = "turns", dates = [2027-01-03T15:00:00Z] }

# TIBOR 6M: TONA の曲線に、TONA 担保のベーシスで決める節点を足す。
[curves.JPY-TIBOR-6M]
currency = "JPY"
roles = ["projection"]
expr.sum = [{ ref = "JPY-TONA" }, { nodes = { part = "basis", pillars = [2027-12-03T00:00:00Z, 2028-12-04T00:00:00Z, 2031-12-03T00:00:00Z, 2036-12-03T00:00:00Z] } }]

# TIBOR 3M: TONA の曲線に、TONA 担保のベーシスで決める節点を足す。
[curves.JPY-TIBOR-3M]
currency = "JPY"
roles = ["projection"]
expr.sum = [{ ref = "JPY-TONA" }, { nodes = { part = "basis", pillars = [2027-12-03T00:00:00Z, 2028-12-04T00:00:00Z, 2031-12-03T00:00:00Z, 2036-12-03T00:00:00Z] } }]

[curves.USD-SOFR]
currency = "USD"
roles = ["discount", "projection"]
pillars = [2027-01-10T00:00:00Z, 2027-03-13T00:00:00Z, 2027-06-10T00:00:00Z, 2027-12-11T00:00:00Z, 2028-12-11T00:00:00Z, 2031-12-11T00:00:00Z, 2036-12-11T00:00:00Z]

# 25bp の付利の JPY 現金担保で割り引く曲線。
[curves.JPY-TONA-25]
currency = "JPY"
roles = ["discount"]
expr.sum = [{ ref = "JPY-TONA" }, { flat = "25bp" }]

# USD の支払いの、JPY と USD の現金担保の CTD（各年の平均フォワードの大きい方）。
[curves.USD-CTD]
currency = "USD"
roles = ["discount"]
expr.max_forward = { parts = [{ ref = "USD-SOFR" }, { ref = "USD-under-JPY" }], grid = [2027-11-27T00:00:00Z, 2028-11-27T00:00:00Z, 2029-11-27T00:00:00Z, 2030-11-27T00:00:00Z, 2031-11-27T00:00:00Z, 2032-11-27T00:00:00Z, 2033-11-27T00:00:00Z, 2034-11-27T00:00:00Z, 2035-11-27T00:00:00Z, 2036-11-27T00:00:00Z] }

[[instruments]]
quote = "JPY.TONA.OIS.1M"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "JPY-TONA-OIS", tenor = "1M" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TONA", part = "meetings" }

[[instruments]]
quote = "JPY.TONA.OIS.3M"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "JPY-TONA-OIS", tenor = "3M" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TONA", part = "meetings" }

[[instruments]]
quote = "JPY.TONA.OIS.6M"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "JPY-TONA-OIS", tenor = "6M" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TONA", part = "meetings" }

[[instruments]]
quote = "JPY.TONA.OIS.1Y"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "JPY-TONA-OIS", tenor = "1Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TONA", part = "meetings" }

[[instruments]]
quote = "JPY.TONA.OIS.2Y"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "JPY-TONA-OIS", tenor = "2Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TONA", part = "long" }

[[instruments]]
quote = "JPY.TONA.OIS.3Y"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "JPY-TONA-OIS", tenor = "3Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TONA", part = "long" }

[[instruments]]
quote = "JPY.TONA.OIS.5Y"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "JPY-TONA-OIS", tenor = "5Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TONA", part = "long" }

[[instruments]]
quote = "JPY.TONA.OIS.7Y"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "JPY-TONA-OIS", tenor = "7Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TONA", part = "long" }

[[instruments]]
quote = "JPY.TONA.OIS.10Y"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "JPY-TONA-OIS", tenor = "10Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TONA", part = "long" }

[[instruments]]
quote = "JPY.TONA.TURN.2026"
kind = "spread"
target = "direct_mark"
determines = { curve = "JPY-TONA", part = "turns", element = 0 }
transform = { turn_spread_to_jump = { value = 0.0136986301369863, basis = "ACT365F" } }

[[instruments]]
quote = "JPY.TIBOR6M-TONA.BASIS.1Y"
kind = "spread"
target = "par_rate"
template = { basis_swap = { convention = "JPY-TIBOR6M-TONA-BASIS", tenor = "1Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TIBOR-6M", part = "basis" }

[[instruments]]
quote = "JPY.TIBOR6M-TONA.BASIS.2Y"
kind = "spread"
target = "par_rate"
template = { basis_swap = { convention = "JPY-TIBOR6M-TONA-BASIS", tenor = "2Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TIBOR-6M", part = "basis" }

[[instruments]]
quote = "JPY.TIBOR6M-TONA.BASIS.5Y"
kind = "spread"
target = "par_rate"
template = { basis_swap = { convention = "JPY-TIBOR6M-TONA-BASIS", tenor = "5Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TIBOR-6M", part = "basis" }

[[instruments]]
quote = "JPY.TIBOR6M-TONA.BASIS.10Y"
kind = "spread"
target = "par_rate"
template = { basis_swap = { convention = "JPY-TIBOR6M-TONA-BASIS", tenor = "10Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TIBOR-6M", part = "basis" }

[[instruments]]
quote = "JPY.TIBOR3M-TONA.BASIS.1Y"
kind = "spread"
target = "par_rate"
template = { basis_swap = { convention = "JPY-TIBOR3M-TONA-BASIS", tenor = "1Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TIBOR-3M", part = "basis" }

[[instruments]]
quote = "JPY.TIBOR3M-TONA.BASIS.2Y"
kind = "spread"
target = "par_rate"
template = { basis_swap = { convention = "JPY-TIBOR3M-TONA-BASIS", tenor = "2Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TIBOR-3M", part = "basis" }

[[instruments]]
quote = "JPY.TIBOR3M-TONA.BASIS.5Y"
kind = "spread"
target = "par_rate"
template = { basis_swap = { convention = "JPY-TIBOR3M-TONA-BASIS", tenor = "5Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TIBOR-3M", part = "basis" }

[[instruments]]
quote = "JPY.TIBOR3M-TONA.BASIS.10Y"
kind = "spread"
target = "par_rate"
template = { basis_swap = { convention = "JPY-TIBOR3M-TONA-BASIS", tenor = "10Y" } }
collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } }
determines = { curve = "JPY-TIBOR-3M", part = "basis" }

[[instruments]]
quote = "USD.SOFR.OIS.1M"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "USD-SOFR-OIS", tenor = "1M" } }
collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
determines = { curve = "USD-SOFR", element = 0 }

[[instruments]]
quote = "USD.SOFR.OIS.3M"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "USD-SOFR-OIS", tenor = "3M" } }
collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
determines = { curve = "USD-SOFR", element = 1 }

[[instruments]]
quote = "USD.SOFR.OIS.6M"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "USD-SOFR-OIS", tenor = "6M" } }
collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
determines = { curve = "USD-SOFR", element = 2 }

[[instruments]]
quote = "USD.SOFR.OIS.1Y"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "USD-SOFR-OIS", tenor = "1Y" } }
collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
determines = { curve = "USD-SOFR", element = 3 }

[[instruments]]
quote = "USD.SOFR.OIS.2Y"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "USD-SOFR-OIS", tenor = "2Y" } }
collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
determines = { curve = "USD-SOFR", element = 4 }

[[instruments]]
quote = "USD.SOFR.OIS.5Y"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "USD-SOFR-OIS", tenor = "5Y" } }
collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
determines = { curve = "USD-SOFR", element = 5 }

[[instruments]]
quote = "USD.SOFR.OIS.10Y"
kind = "rate"
target = "par_rate"
template = { ois = { convention = "USD-SOFR-OIS", tenor = "10Y" } }
collateral = { cash = { ccy = "USD", index = "USD-SOFR", spread = "0bp" } }
determines = { curve = "USD-SOFR", element = 6 }

# USDJPY: スポット、TN、3M と 6M のポイントと、外貨担保曲線の長期の MtM XCCY。
# 評価の 17 時は当日の受渡の締め（15 時）を過ぎているので、ON は連鎖に使わない。
[[fx]]
pair = "USDJPY"
collaterals = [{ id = "market", base = "USD-under-JPY", quote = "JPY-TONA" }]

[fx.direct]
market_collateral = "market"
quotes = { spot = "USDJPY.SPOT", on = "USDJPY.ON", tn = "USDJPY.TN", forwards = [{ settlement = 2027-02-26, quote = "USDJPY.3M" }, { settlement = 2027-05-27, quote = "USDJPY.6M" }] }

[fx.direct.implied]
id = "USD-under-JPY"
long_end = [
    { pillar = 2027-12-11T00:00:00Z, quote = "USDJPY.XCCY.1Y", convention = "USDJPY-MTM-XCCY", tenor = "1Y", collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } } },
    { pillar = 2028-12-11T00:00:00Z, quote = "USDJPY.XCCY.2Y", convention = "USDJPY-MTM-XCCY", tenor = "2Y", collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } } },
    { pillar = 2031-12-11T00:00:00Z, quote = "USDJPY.XCCY.5Y", convention = "USDJPY-MTM-XCCY", tenor = "5Y", collateral = { cash = { ccy = "JPY", index = "JPY-TONA", spread = "0bp" } } },
]

# 金利先物の凸性の補正に使う JPY の LGM のマーク。
[lgm.JPY]
mean_reversion = ["JPY.LGM.MR"]
volatility = ["JPY.LGM.VOL"]
`````

`````
examples/m1-rates/quotes/snapshot.csv
quote_id,value
JPY.TONA.OIS.1M,0.48%
JPY.TONA.OIS.3M,0.52%
JPY.TONA.OIS.6M,0.58%
JPY.TONA.OIS.1Y,0.65%
JPY.TONA.OIS.2Y,0.78%
JPY.TONA.OIS.3Y,0.88%
JPY.TONA.OIS.5Y,1.02%
JPY.TONA.OIS.7Y,1.15%
JPY.TONA.OIS.10Y,1.32%
JPY.TONA.TURN.2026,50bp
JPY.TIBOR6M-TONA.BASIS.1Y,18bp
JPY.TIBOR6M-TONA.BASIS.2Y,20bp
JPY.TIBOR6M-TONA.BASIS.5Y,24bp
JPY.TIBOR6M-TONA.BASIS.10Y,27bp
JPY.TIBOR3M-TONA.BASIS.1Y,12bp
JPY.TIBOR3M-TONA.BASIS.2Y,14bp
JPY.TIBOR3M-TONA.BASIS.5Y,17bp
JPY.TIBOR3M-TONA.BASIS.10Y,19bp
USD.SOFR.OIS.1M,3.70%
USD.SOFR.OIS.3M,3.62%
USD.SOFR.OIS.6M,3.50%
USD.SOFR.OIS.1Y,3.35%
USD.SOFR.OIS.2Y,3.25%
USD.SOFR.OIS.5Y,3.35%
USD.SOFR.OIS.10Y,3.55%
USDJPY.SPOT,152.40
USDJPY.ON,-2.5
USDJPY.TN,-2.4
USDJPY.3M,-112
USDJPY.6M,-218
USDJPY.XCCY.1Y,-35bp
USDJPY.XCCY.2Y,-40bp
USDJPY.XCCY.5Y,-45bp
JPY.LGM.MR,0.03
JPY.LGM.VOL,0.006
`````

`````
examples/m1-rates/refdata/calendars.toml
version = 1

# 2026 年 11 月から 2027 年 12 月までは主な祝日を、それ以後は週末だけを休日にするサンプルの暦。
[[calendar]]
id = "JPTO"
from = 2025-01-01
through = 2037-12-31
weekend = ["Sat", "Sun"]

[[calendar]]
id = "USNY"
from = 2025-01-01
through = 2037-12-31
weekend = ["Sat", "Sun"]
`````

`````
examples/m1-rates/refdata/conventions.toml
version = 1

[[convention]]
id = "JPY-TONA-OIS"

[convention.kind.ois]
spot_lag = "2BD"
fixed = { frequency = "1Y", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }

[convention.kind.ois.float]
index = "JPY-TONA"
observation = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } }

[[convention]]
id = "USD-SOFR-OIS"

[convention.kind.ois]
spot_lag = "2BD"
fixed = { frequency = "1Y", day_count = "ACT360", bdc = "modified_following", calendar = "USNY" }
payment = { anchor = "period_end", offset = "2BD", calendar = "USNY", bdc = "following" }

[convention.kind.ois.float]
index = "USD-SOFR"
observation = { days = { business_days = "USNY" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } }

[[convention]]
id = "JPY-TIBOR6M-TONA-BASIS"

[convention.kind.basis_swap]
spot_lag = "2BD"

[[convention.kind.basis_swap.legs]]
currency = "JPY"
periods = { frequency = "6M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
coupon = { compounded = { index = "JPY-TONA", observation = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }

[[convention.kind.basis_swap.legs]]
currency = "JPY"
periods = { frequency = "6M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
coupon = { term = { index = "JPY-TIBOR-6M" } }

[[convention]]
id = "JPY-TIBOR3M-TONA-BASIS"

[convention.kind.basis_swap]
spot_lag = "2BD"

[[convention.kind.basis_swap.legs]]
currency = "JPY"
periods = { frequency = "3M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
coupon = { compounded = { index = "JPY-TONA", observation = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }

[[convention.kind.basis_swap.legs]]
currency = "JPY"
periods = { frequency = "3M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
coupon = { term = { index = "JPY-TIBOR-3M" } }

[[convention]]
id = "JPY-TIBOR6M-IRS"

[convention.kind.irs]
spot_lag = "2BD"

[[convention.kind.irs.legs]]
currency = "JPY"
periods = { frequency = "6M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
coupon = "fixed"

[[convention.kind.irs.legs]]
currency = "JPY"
periods = { frequency = "6M", day_count = "ACT365F", bdc = "modified_following", calendar = "JPTO" }
payment = { anchor = "period_end", offset = "0BD", calendar = "JPTO", bdc = "following" }
coupon = { term = { index = "JPY-TIBOR-6M" } }

[[convention]]
id = "USDJPY-MTM-XCCY"

[convention.kind.xccy_swap]
spot_lag = "2BD"

[[convention.kind.xccy_swap.legs]]
currency = "JPY"
periods = { frequency = "3M", day_count = "ACT365F", bdc = "modified_following", calendar = { union = ["JPTO", "USNY"] } }
payment = { anchor = "period_end", offset = "0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" }
coupon = { compounded = { index = "JPY-TONA", observation = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }

[[convention.kind.xccy_swap.legs]]
currency = "USD"
periods = { frequency = "3M", day_count = "ACT360", bdc = "modified_following", calendar = { union = ["JPTO", "USNY"] } }
payment = { anchor = "period_end", offset = "0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" }
coupon = { compounded = { index = "USD-SOFR", observation = { days = { business_days = "USNY" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }
fx_reset = "USDJPY-TKY-1500"

[[convention]]
id = "JPY-TIBOR6M-FRA"
kind.fra = { index = "JPY-TIBOR-6M", spot_lag = "2BD" }

[[convention]]
id = "USDJPY-FWD"
kind.fx_forward = { pair = "USDJPY", bdc = "modified_following" }

[[convention]]
id = "JPY-TONA3M-FUT"

[convention.kind.rate_future]
multiplier = 25_000_000.0
settlement_lag = "1BD"
settlement_time = { time = "10:00", zone = "Asia/Tokyo" }
contracts = { "2027-03" = { settlement_index = "JPY-TONA3M-FUT-2027-03", final_clearing = 2027-06-15, rate = { compounded = { start = 2027-03-17, end = 2027-06-16 } } } }

# USDJPY は T+2、東京 15 時の受渡、pip は 0.01 円。評価の 17 時には当日の受渡の締めを過ぎている。
[[fx]]
pair = "USDJPY"
spot_lag = "two_business_days"
settlement = { time = "15:00", zone = "Asia/Tokyo" }
pip_exponent = 2
`````

`````
examples/m1-rates/refdata/currencies.toml
version = 1

[[currency]]
code = "JPY"
minor_units = 0
payment_rounding = { increment = 1.0, mode = "half_even" }
deliverability = "deliverable"
settlement_time = { time = "17:00", zone = "Asia/Tokyo" }
settlement_calendar = "JPTO"

[[currency]]
code = "USD"
minor_units = 2
payment_rounding = { increment = 0.01, mode = "half_even" }
deliverability = "deliverable"
settlement_time = { time = "03:00", zone = "America/New_York" }
settlement_calendar = "USNY"
`````

`````
examples/m1-rates/refdata/holidays/JPTO.csv
date
2026-11-03
2026-11-23
2026-12-31
2027-01-01
2027-01-11
2027-02-11
2027-02-23
2027-03-22
2027-04-29
2027-05-03
2027-05-04
2027-05-05
2027-07-19
2027-08-11
2027-09-20
2027-09-23
2027-10-11
2027-11-03
2027-11-23
2027-12-31
`````

`````
examples/m1-rates/refdata/holidays/USNY.csv
date
2026-11-11
2026-11-26
2026-12-25
2027-01-01
2027-01-18
2027-02-15
2027-05-31
2027-06-18
2027-07-05
2027-09-06
2027-10-11
2027-11-11
2027-11-25
2027-12-24
`````

`````
examples/m1-rates/refdata/indices.toml
version = 1

[[index]]
id = "JPY-TONA"
calendar = "JPTO"
publication = { business_days = 1, calendar = "JPTO", time = "17:00", zone = "Asia/Tokyo" }
observation_time = { time = "17:00", zone = "Asia/Tokyo" }
kind = { overnight_rate = { ccy = "JPY", day_count = "ACT365F" } }

[[index]]
id = "USD-SOFR"
calendar = "USNY"
publication = { business_days = 1, calendar = "USNY", time = "03:00", zone = "America/New_York" }
observation_time = { time = "03:00", zone = "America/New_York" }
kind = { overnight_rate = { ccy = "USD", day_count = "ACT360" } }

[[index]]
id = "JPY-TIBOR-6M"
calendar = "JPTO"
publication = { business_days = 0, calendar = "JPTO", time = "11:00", zone = "Asia/Tokyo" }
observation_time = { time = "11:00", zone = "Asia/Tokyo" }
kind = { term_rate = { ccy = "JPY", day_count = "ACT365F", tenor = "6M", value_offset = "2BD", bdc = "modified_following", end_of_month = "ignore" } }

[[index]]
id = "JPY-TIBOR-3M"
calendar = "JPTO"
publication = { business_days = 0, calendar = "JPTO", time = "11:00", zone = "Asia/Tokyo" }
observation_time = { time = "11:00", zone = "Asia/Tokyo" }
kind = { term_rate = { ccy = "JPY", day_count = "ACT365F", tenor = "3M", value_offset = "2BD", bdc = "modified_following", end_of_month = "ignore" } }

# 東京 15 時の USDJPY の公表値。MtM XCCY の元本のリセットに使う。
[[index]]
id = "USDJPY-TKY-1500"
calendar = { union = ["JPTO", "USNY"] }
publication = { business_days = 0, calendar = { union = ["JPTO", "USNY"] }, time = "15:00", zone = "Asia/Tokyo" }
observation_time = { time = "15:00", zone = "Asia/Tokyo" }
kind = { fx_fixing = { pair = "USDJPY" } }

# 2027 年 3 月限の 3 か月 TONA 先物の日次清算価格。
[[index]]
id = "JPY-TONA3M-FUT-2027-03"
calendar = "JPTO"
publication = { business_days = 0, calendar = "JPTO", time = "15:30", zone = "Asia/Tokyo" }
observation_time = { time = "15:30", zone = "Asia/Tokyo" }
kind = { rate_futures_settlement = { contract = "JPY-TONA3M-2027-03", underlying = "JPY-TONA" } }
`````

`````
examples/m1-rates/refdata/securities.toml
version = 1

[[security]]
id = "JPY-TONA3M-2027-03"
kind = "listed_future"
currency = "JPY"
unit = "contract"
valuation = "contract"
`````

`````
examples/m1-rates/refdata/zones.toml
version = 1

[[zone]]
id = "Asia/Tokyo"
from = 2025-01-01T00:00:00Z
through = 2038-01-01T00:00:00Z
initial_offset_seconds = 32400

[[zone]]
id = "America/New_York"
from = 2025-01-01T00:00:00Z
through = 2038-01-01T00:00:00Z
initial_offset_seconds = -18000
transitions = [
    { at = 2025-03-09T07:00:00Z, offset_seconds = -14400 },
    { at = 2025-11-02T06:00:00Z, offset_seconds = -18000 },
    { at = 2026-03-08T07:00:00Z, offset_seconds = -14400 },
    { at = 2026-11-01T06:00:00Z, offset_seconds = -18000 },
    { at = 2027-03-14T07:00:00Z, offset_seconds = -14400 },
    { at = 2027-11-07T06:00:00Z, offset_seconds = -18000 },
    { at = 2028-03-12T07:00:00Z, offset_seconds = -14400 },
    { at = 2028-11-05T06:00:00Z, offset_seconds = -18000 },
    { at = 2029-03-11T07:00:00Z, offset_seconds = -14400 },
    { at = 2029-11-04T06:00:00Z, offset_seconds = -18000 },
    { at = 2030-03-10T07:00:00Z, offset_seconds = -14400 },
    { at = 2030-11-03T06:00:00Z, offset_seconds = -18000 },
    { at = 2031-03-09T07:00:00Z, offset_seconds = -14400 },
    { at = 2031-11-02T06:00:00Z, offset_seconds = -18000 },
    { at = 2032-03-14T07:00:00Z, offset_seconds = -14400 },
    { at = 2032-11-07T06:00:00Z, offset_seconds = -18000 },
    { at = 2033-03-13T07:00:00Z, offset_seconds = -14400 },
    { at = 2033-11-06T06:00:00Z, offset_seconds = -18000 },
    { at = 2034-03-12T07:00:00Z, offset_seconds = -14400 },
    { at = 2034-11-05T06:00:00Z, offset_seconds = -18000 },
    { at = 2035-03-11T07:00:00Z, offset_seconds = -14400 },
    { at = 2035-11-04T06:00:00Z, offset_seconds = -18000 },
    { at = 2036-03-09T07:00:00Z, offset_seconds = -14400 },
    { at = 2036-11-02T06:00:00Z, offset_seconds = -18000 },
    { at = 2037-03-08T07:00:00Z, offset_seconds = -14400 },
    { at = 2037-11-01T06:00:00Z, offset_seconds = -18000 },
]
`````

`````
examples/m1-rates/trades/csa.toml
version = 1

[[trade]]
id = "USD-OIS-JPY-CASH"
margin_set = "JPY-CASH"
book = "CSA"
trade_date = 2026-11-24

[trade.product]
type = "swap"
convention = "USD-SOFR-OIS"
tenor = "2Y"
notional = "10_000_000 USD"
fixed_rate = "3.40%"
fixed_side = "receive"

[[trade]]
id = "USD-OIS-USD-CASH"
margin_set = "USD-CASH"
book = "CSA"
trade_date = 2026-11-24

[trade.product]
type = "swap"
convention = "USD-SOFR-OIS"
tenor = "2Y"
notional = "10_000_000 USD"
fixed_rate = "3.40%"
fixed_side = "receive"

[[trade]]
id = "USD-OIS-CTD"
margin_set = "CTD-JPY-USD"
book = "CSA"
trade_date = 2026-11-24

[trade.product]
type = "swap"
convention = "USD-SOFR-OIS"
tenor = "2Y"
notional = "10_000_000 USD"
fixed_rate = "3.40%"
fixed_side = "receive"

[[trade]]
id = "USD-OIS-UNCOLLATERALIZED"
margin_set = "UNCOLLATERALIZED"
book = "CSA"
trade_date = 2026-11-24

[trade.product]
type = "swap"
convention = "USD-SOFR-OIS"
tenor = "2Y"
notional = "10_000_000 USD"
fixed_rate = "3.40%"
fixed_side = "receive"

[[trade]]
id = "JPY-IRS-5Y-CASH-25"
margin_set = "JPY-CASH-25"
book = "CSA"
trade_date = 2026-11-24

[trade.product]
type = "swap"
legs = [
  { side = "pay", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "1.10%" } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
  { side = "receive", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { term = { index = "JPY-TIBOR-6M" } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
]
`````

`````
examples/m1-rates/trades/fx.toml
version = 1

[[trade]]
id = "USDJPY-FWD-1Y"
margin_set = "JPY-CASH"
book = "FX"
trade_date = 2026-11-24

[trade.product]
type = "fx_forward"
pair = "USDJPY"
notional = "10_000_000 USD"
rate = 148.0
settlement = 2027-11-29T06:00:00Z

[[trade]]
id = "USDJPY-FXSWAP-6M"
margin_set = "JPY-CASH"
book = "FX"
trade_date = 2026-11-24

[trade.product]
type = "fx_swap"
pair = "USDJPY"
notional = "10_000_000 USD"
near = { rate = 152.40, settlement = 2026-11-27T06:00:00Z }
far = { rate = 150.20, settlement = 2027-05-27T06:00:00Z }
`````

`````
examples/m1-rates/trades/rates.toml
version = 1

[[trade]]
id = "JPY-DEPO-3M"
margin_set = "JPY-CASH"
book = "JPY-MM"
trade_date = 2026-11-24

[trade.product]
type = "deposit"
start = 2026-11-26
maturity = 2027-02-26
notional = "1_000_000_000 JPY"
rate = "0.50%"
day_count = "ACT365F"

[[trade]]
id = "JPY-FRA-6X12"
margin_set = "JPY-CASH"
book = "JPY-MM"
trade_date = 2026-11-24

[trade.product]
type = "fra"
notional = "1_000_000_000 JPY"
index = "JPY-TIBOR-6M"
fixing = 2027-05-24
fixed_rate = "0.75%"

[[trade]]
id = "JPY-TONA-FUT"
margin_set = "JPY-CASH"
book = "JPY-FUTURES"
trade_date = 2026-11-24

[trade.product]
type = "rate_future"
settlement_index = "JPY-TONA3M-FUT-2027-03"
contracts = 10
trade_price = 0.9935
day_count = "ACT365F"

[[trade]]
id = "JPY-IRS-5Y"
margin_set = "JPY-CASH"
book = "JPY-SWAPS"
trade_date = 2026-11-24

[trade.product]
type = "swap"
legs = [
  { side = "pay", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "1.10%" } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
  { side = "receive", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { term = { index = "JPY-TIBOR-6M" } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
]

[[trade]]
id = "JPY-BASIS-3M-5Y"
margin_set = "JPY-CASH"
book = "JPY-SWAPS"
trade_date = 2026-11-24

[trade.product]
type = "swap"
legs = [
  { side = "pay", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", spread = "20bp", rule = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "3M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
  { side = "receive", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { term = { index = "JPY-TIBOR-3M" } }, schedule = { start = 2026-11-26, end = { after_start = "5Y" }, step = "3M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
]

[[trade]]
id = "JPY-OIS-LOOKBACK"
margin_set = "JPY-CASH"
book = "JPY-SWAPS"
trade_date = 2026-11-24

[trade.product]
type = "swap"
legs = [
  { side = "receive", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "0.70%" } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
  { side = "pay", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", rule = { days = { business_days = "JPTO" }, method = { lookback = 2 }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
]

[[trade]]
id = "JPY-OIS-SHIFT"
margin_set = "JPY-CASH"
book = "JPY-SWAPS"
trade_date = 2026-11-24

[trade.product]
type = "swap"
legs = [
  { side = "receive", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "0.70%" } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
  { side = "pay", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", rule = { days = { business_days = "JPTO" }, method = { observation_shift = 2 }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
]

[[trade]]
id = "JPY-OIS-LOCKOUT"
margin_set = "JPY-CASH"
book = "JPY-SWAPS"
trade_date = 2026-11-24

[trade.product]
type = "swap"
legs = [
  { side = "receive", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "0.70%" } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
  { side = "pay", notional = { fixed = "5_000_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", rule = { days = { business_days = "JPTO" }, method = { lockout = 1 }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-26, end = { after_start = "2Y" }, step = "1Y", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+2BD", calendar = "JPTO", bdc = "following" } } },
]

[[trade]]
id = "JPY-IRS-AMENDED"
margin_set = "JPY-CASH"
book = "JPY-SWAPS"
trade_date = 2026-05-28
fees = [{ amount = "2_000_000 JPY", side = "receive", at = 2026-11-26T08:00:00Z }]

[trade.product]
type = "amended"
from = 2026-11-24T08:00:00Z
after = { scale = 0.4 }

[trade.product.before]
type = "swap"
legs = [
  { side = "pay", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { fixed = { rate = "0.90%" } }, schedule = { start = 2026-06-01, end = { after_start = "3Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
  { side = "receive", notional = { fixed = "10_000_000_000 JPY" }, day_count = "ACT365F", coupon = { term = { index = "JPY-TIBOR-6M" } }, schedule = { start = 2026-06-01, end = { after_start = "3Y" }, step = "6M", stub = "no_stub", calendar = "JPTO", bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = "JPTO", bdc = "following" } } },
]
`````

`````
examples/m1-rates/trades/xccy.toml
version = 1

[[trade]]
id = "USDJPY-MTM-XCCY-2Y"
margin_set = "JPY-CASH"
book = "XCCY"
trade_date = 2026-11-24

[trade.product]
type = "swap"
legs = [
  { side = "pay", notional = { fixed = "1_524_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", spread = "-38bp", rule = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-27, end = { after_start = "2Y" }, step = "3M", stub = "no_stub", calendar = { union = ["JPTO", "USNY"] }, bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" } }, exchange = "initial_and_final" },
  { side = "receive", notional = { fx_reset = { base = "1_524_000_000 JPY", fixing = "USDJPY-TKY-1500" } }, day_count = "ACT360", coupon = { compounded = { index = "USD-SOFR", rule = { days = { business_days = "USNY" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-27, end = { after_start = "2Y" }, step = "3M", stub = "no_stub", calendar = { union = ["JPTO", "USNY"] }, bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" } }, exchange = "initial_and_final" },
]

[[trade]]
id = "USDJPY-XCCY-2Y"
margin_set = "JPY-CASH"
book = "XCCY"
trade_date = 2026-11-24

[trade.product]
type = "swap"
legs = [
  { side = "pay", notional = { fixed = "1_524_000_000 JPY" }, day_count = "ACT365F", coupon = { compounded = { index = "JPY-TONA", spread = "-40bp", rule = { days = { business_days = "JPTO" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-27, end = { after_start = "2Y" }, step = "3M", stub = "no_stub", calendar = { union = ["JPTO", "USNY"] }, bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" } }, exchange = "initial_and_final" },
  { side = "receive", notional = { fixed = "10_000_000 USD" }, day_count = "ACT360", coupon = { compounded = { index = "USD-SOFR", rule = { days = { business_days = "USNY" }, aggregation = { compound = { form = "simple_factor", spread = "exclusive" } } } } }, schedule = { start = 2026-11-27, end = { after_start = "2Y" }, step = "3M", stub = "no_stub", calendar = { union = ["JPTO", "USNY"] }, bdc = "modified_following", accrual = "adjusted", payment = { anchor = "period_end", offset = "+0BD", calendar = { union = ["JPTO", "USNY"] }, bdc = "following" } }, exchange = "initial_and_final" },
]
`````

`````
rust-toolchain.toml
[toolchain]
channel = "1.100.0"
profile = "minimal"
components = ["clippy", "rust-analyzer", "rust-src", "rustfmt"]
`````
