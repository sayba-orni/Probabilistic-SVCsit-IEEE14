# ProbabilisticSVCsitIEEE14
We study where to place a single Static VAR Compensator (SVC) on the IEEE-14 bus system
when solar and wind injections, loads, and N-1 line outages are uncertain. We generate hundreds
of scenarios via Monte Carlo sampling of weather, load noise, and random outages; solve a
robust Newton–Raphson (NR) power flow with an SVC model (either at a bus or mid-line); and
rank candidates using a risk-aware objective: active-power losses plus a weighted voltage-
violation penalty. We summarize risk with CVaR0.90 (tail average of the worst 10% outcomes).

Electricity systems don’t behave the same every day:
• Solar power depends on sunshine.
• Wind power depends on wind speed.
• Loads (demand) vary randomly.
• Transmission lines can fail (N-1 outage = one line out at a time).
All this uncertainty means voltages may go too high/low and power losses may increase.
An SVC can help by injecting or absorbing reactive power to stabilize voltages. But its impact
depends a lot on where you place it in the network.
Instead of assuming one “average” day, we create hundreds of possible futures (scenarios)
using Monte Carlo sampling:
1. Pick random weather (irradiance for solar, wind speed for wind farm).
2. Add random noise to loads (like ±5%).
3. Randomly remove one line with small probability (outage).
Each scenario is one possible “state of the world.”
