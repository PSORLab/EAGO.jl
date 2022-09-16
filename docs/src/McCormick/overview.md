# Overview

EAGO provides a library of McCormick relaxations in native Julia code. The EAGO optimizer supports
relaxing functions using **nonsmooth McCormick relaxations** ([Mitsos2009](https://epubs.siam.org/doi/abs/10.1137/080717341), [Scott2011](https://link.springer.com/article/10.1007/s10898-011-9664-7)), **smooth McCormick relaxations** ([Khan2016](https://link.springer.com/article/10.1007/s10898-016-0440-6), [Khan2018](https://link.springer.com/article/10.1007/s10898-017-0601-2), [Khan2019](https://www.tandfonline.com/doi/abs/10.1080/02331934.2018.1534108)), and **multi-variant McCormick relaxations** ([Tsoukalas2014](https://link.springer.com/article/10.1007/s10898-014-0176-0); a variant of **subgradient-based interval refinement** ([Najman2017](https://link.springer.com/article/10.1007/s10898-016-0470-0))). For functions with arbitrarily differentiable relaxations, the differentiable constant μ can be modified by adjusting a constant value in the package. Additionally, validated and nonvalidated interval bounds are supported via [**IntervalArithmetic.jl**](https://github.com/JuliaIntervals/IntervalArithmetic.jl) which is reexported. The basic McCormick operator and reverse McCormick operator ([Wechsung2015](https://link.springer.com/article/10.1007/s10898-015-0303-6)) libraries are included in two dependent subpackages which can loaded and used independently:
- **[McCormick.jl](https://github.com/PSORLab/McCormick.jl)**: A library of forward-mode and implciit McCormick operators.
- **[ReverseMcCormick.jl](https://github.com/PSORLab/ReverseMcCormick.jl)**: A reverse-mode McCormick operator library.

## NaN Numerics
When a relaxation is computed at an undefined point or over an unbounded domain, the resulting relaxation is defined as "not a number" (`NaN`) rather than throwing an error. This allows algorithms to check for these cases without resorting to `try-catch` statements. Moreover, when the interval domain is extensive enough to cause a domain violation, an `x::MC` structure is returned that satisfies `isnan(x) === true`.

- **Khan KA, Watson HAJ, Barton PI (2017).** Differentiable McCormick relaxations. *Journal of Global Optimization*, 67(4): 687-729.
- **Khan KA, Wilhelm ME, Stuber MD, Cao H, Watson HAJ, Barton PI (2018).** Corrections to: Differentiable McCormick relaxations. *Journal of Global Optimization*, 70(3): 705-706.
- **Khan KA (2019).** Whitney differentiability of optimal-value functions for bound-constrained convex programming problems. *Optimization*, 68(2-3): 691-711
- **Mitsos A, Chachuat B, and Barton PI. (2009).** McCormick-based relaxations of algorithms. *SIAM Journal on Optimization*, 20(2): 573–601.
- **Najman J, Bongratz D, Tsoukalas A, and Mitsos A (2017).** Erratum to: Multivariate McCormick relaxations. *Journal of Global Optimization*, 68: 219-225.
- **Scott JK,  Stuber MD, and Barton PI. (2011).** Generalized McCormick relaxations. *Journal of Global Optimization*, 51(4): 569–606.
- **Stuber MD, Scott JK, Barton PI (2015).** Convex and concave relaxations of implicit functions. *Optim. Methods Softw.*, 30(3): 424–460
- **Tsoukalas A and Mitsos A (2014).** Multivariate McCormick Relaxations. *Journal of Global Optimization*, 59: 633–662.
- **Wechsung A, Scott JK, Watson HAJ, and Barton PI. (2015).** Reverse propagation of McCormick relaxations. *Journal of Global Optimization*, 63(1): 1-36.
