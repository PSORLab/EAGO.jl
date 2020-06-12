# **Currently supported operators**

The operators currently supported are listed below. The operators with a check box
have been subject to a large degree of scrutiny and have been implemented for
both forward and reverse McCormick relaxations ([Wechsung2015](https://link.springer.com/article/10.1007/s10898-015-0303-6)). Each McCormick object is associated with a
parameter `T <: RelaxTag` which is either `NS` for nonsmooth relaxations ([Mitsos2009](https://epubs.siam.org/doi/abs/10.1137/080717341), [Scott2011](https://link.springer.com/article/10.1007/s10898-011-9664-7)), `MV` for multivariate relaxations ([Tsoukalas2014](https://link.springer.com/article/10.1007/s10898-014-0176-0), [Najman2017](https://link.springer.com/article/10.1007/s10898-016-0470-0)),
and `Diff` for differentiable relaxations ([Khan2016](https://link.springer.com/article/10.1007/s10898-016-0440-6), [Khan2018](https://link.springer.com/article/10.1007/s10898-017-0601-2), [Khan2019](https://www.tandfonline.com/doi/abs/10.1080/02331934.2018.1534108)). Conversion between `MV`, `NS`, and `Diff` relax tags are not currently supported. Convex and concave envelopes are used to compute relaxations of univariate functions.

### **Univariate McCormick Operators**

Arbitrarily differentiable relaxations can be constructed for the following operators:

- [x] **Inverse** (inv)
- [x] **Logarithms** (log, log2, log10)
- [x] **Exponential Functions** (exp, exp2, exp10)
- [x] **Square Root** (sqrt)
- [x] **Absolute Value** (abs)

Both nonsmooth and Whitney-1 (once differentiable) relaxations are supported:

- [x] **Step Functions** (step, sign)
- [x] **Trignometric Functions** (sin, cos, tan)
- [x] **Inverse Trignometric Functions** (asin, acos, atan)
- [x] **Hyperbolic Functions** (sinh, cosh, tanh)
- [x] **Inverse Hyperbolic Functions** (asinh, acosh, atanh)

### **Bivariate Operators: McCormick & McCormick**

The following bivariant operators are supported for two **MC** objects. Both nonsmooth and Whitney-1 (once differentiable) relaxations are supported.

- [x] **multiplication** (\*)
- [x] **division** (/)

Arbitrarily differentiable relaxations can be constructed for the following operators:

- [x] **addition** (+)
- [x] **subtraction** (-)
- [x] **minimization** (min)
- [x] **maximization** (max)

### References
- **Khan KA, Watson HAJ, Barton PI (2017).** Differentiable McCormick relaxations. *Journal of Global Optimization*, 67(4):687-729.
- **Khan KA, Wilhelm ME, Stuber MD, Cao H, Watson HAJ, Barton PI (2018).** Corrections to: Differentiable McCormick relaxations. *Journal of Global Optimization*, 70(3):705-706.
- **Khan KA (2019).** Whitney differentiability of optimal-value functions for bound-constrained convex programming problems. *Optimization* 68(2-3): 691-711
- **Mitsos A, Chachuat B, and Barton PI. (2009).** McCormick-based relaxations of algorithms. *SIAM Journal on Optimization*, 20(2):573–601.
- **Najman J, Bongratz D, Tsoukalas A, and Mitsos A (2017).** Erratum to: Multivariate McCormick relaxations. *Journal of Global Optimization*, 68:219-225.
- **Scott JK,  Stuber MD, and Barton PI. (2011).** Generalized McCormick relaxations. *Journal of Global Optimization*, 51(4):569–606.
- **Stuber MD, Scott JK, Barton PI (2015).** Convex and concave relaxations of implicit functions. *Optim. Methods Softw.* 30(3), 424–460
- **Tsoukalas A and Mitsos A (2014).** Multivariate McCormick Relaxations. *Journal of Global Optimization*, 59:633–662.
- **Wechsung A, Scott JK, Watson HAJ, and Barton PI. (2015).** Reverse propagation of McCormick relaxations. *Journal of Global Optimization* 63(1):1-36.
