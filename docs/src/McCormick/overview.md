# Overview

EAGO provides a library of McCormick relaxations in native Julia code. It supports
relaxing functions using both **nonsmooth McCormick relaxations** (Mitsos2009), **smooth McCormick relaxations** (Khan2017), **multi-variant McCormick relaxations** (Tsoukalas2014), as well
as **subgradient-based interval refinement** (Najman2017). For functions with
arbitrarily differentiable relaxations, the differentiable can be modified by adjusting a constant value. Additionally, and nonvalidated validated interval bounds are supported via **ValidatedNumerics.jl**.

## Guard Operators
By default, the McCormick operators used by EAGO perform are **unguarded**. When
the interval domain is extensive enough to cause a domain violation an error is
thrown. If models are known to be well-posed overloading results in a domain violation due to expansiveness of the interval bounds e.g. x/y on Y = [-2, 2], a small region is cut out around the point of domain violation. The relaxation of x/y is then computed as the union of the
relaxation of the X/Y1 and X/Y2. 
