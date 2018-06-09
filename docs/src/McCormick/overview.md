## **McCormick Operator Capabilities**

EAGO provides a library of McCormick relaxations in native Julia code. It supports
relaxing functions using both **nonsmooth McCormick relaxations** (Mitsos2009), **smooth McCormick relaxations** (Khan2017), **multi-variant McCormick relaxations** (Tsoukalas2014), as well
as **subgradient-based interval refinement** (Najman2017). For functions with
arbitrarily differentiable relaxations, the differentiable can be modified by adjusting a constant value. Additionally, validated interval bounds are supported via **ValidatedNumerics.jl** and nonvalidated interval operators are available through use of the **MCInterval{T}** object included in EAGO (which is essentially just a copy of the **IntervalArithmetic.jl** and **IntervalContractor.jl** library with corrected rounding features removed).
