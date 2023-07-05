# Currently Supported Operators

The operators currently supported are listed below. The operators with a check box have been subject to a large degree of scrutiny and have been implemented for both forward and reverse McCormick relaxations ([Wechsung2015](https://link.springer.com/article/10.1007/s10898-015-0303-6)). Each McCormick object is associated with a parameter `T <: RelaxTag` which is either `NS` for nonsmooth relaxations ([Mitsos2009](https://epubs.siam.org/doi/abs/10.1137/080717341), [Scott2011](https://link.springer.com/article/10.1007/s10898-011-9664-7)), `MV` for multivariate relaxations ([Tsoukalas2014](https://link.springer.com/article/10.1007/s10898-014-0176-0), [Najman2017](https://link.springer.com/article/10.1007/s10898-016-0470-0)), or `Diff` for differentiable relaxations ([Khan2016](https://link.springer.com/article/10.1007/s10898-016-0440-6), [Khan2018](https://link.springer.com/article/10.1007/s10898-017-0601-2), [Khan2019](https://www.tandfonline.com/doi/abs/10.1080/02331934.2018.1534108)). Conversion between `NS`, `MV`, and `Diff` relax tags is not currently supported. Convex and concave envelopes are used to compute relaxations of univariate functions.

## Univariate McCormick Operators

Arbitrarily differentiable relaxations can be constructed for the following operators:

- **Inverse** (`inv`)
- **Logarithms** (`log`, `log2`, `log10`)
- **Exponential Functions** (`exp`, `exp2`, `exp10`)
- **Square Root** (`sqrt`)
- **Absolute Value** (`abs`)

Both nonsmooth and Whitney-1 (once differentiable) relaxations are supported for the following operators:

- **Step Functions** (`step`, `sign`)
- **Trigonometric Functions** (`sin`, `cos`, `tan`)
- **Inverse Trigonometric Functions** (`asin`, `acos`, `atan`)
- **Hyperbolic Functions** (`sinh`, `cosh`, `tanh`)
- **Inverse Hyperbolic Functions** (`asinh`, `acosh`, `atanh`)
- **Common Activation Functions** (`relu`, `leaky_relu`, `param_relu`, `sigmoid`, `bisigmoid`,
                                       `softsign`, `softplu`s, `maxtanh`, `pentanh`,
                                       `gelu`, `elu`, `selu`, `swish`)
- **Special Functions** (`erf`)

## Bivariate McCormick Operators

The following bivariate operators are supported for two [`MC`](@ref McCormick.MC) objects. Both nonsmooth and Whitney-1 (once differentiable) relaxations are supported.

- **Multiplication** (`*`)
- **Division** (`/`)

Arbitrarily differentiable relaxations can be constructed for the following operators:

- **Addition** (`+`)
- **Subtraction** (`-`)
- **Minimization** (`min`)
- **Maximization** (`max`)

## Common Subexpressions

The following functions can be used in place of common subexpressions encountered in optimization and will result in improved performance (in each case, the standard McCormick composition rules are often more expansive).

```@docs
xexpax
arh
xlogx
mm
```

## Bound Setting Functions

The following functions are used to specify that known bounds on a subexpression exist and that the relaxation/interval bounds propagated should make use of this information. The utility functions can be helpful in avoiding domain violations that arise due to the overly expansive nature of composite relaxations. Improper use of these functions may lead to cases in which the resulting relaxations are empty, so the user is encouraged to use discretion.

```@docs
positive
negative
lower_bnd
upper_bnd
bnd
```

## Specialized Activation Functions

```@docs
    pentanh
    leaky_relu
    param_relu
    maxtanh
```

## References

- **Khan KA, Watson HAJ, Barton PI (2017).** Differentiable McCormick relaxations. *Journal of Global Optimization*, 67(4): 687-729.
- **Khan KA, Wilhelm ME, Stuber MD, Cao H, Watson HAJ, Barton PI (2018).** Corrections to: Differentiable McCormick relaxations. *Journal of Global Optimization*, 70(3): 705-706.
- **Khan KA (2019).** Whitney differentiability of optimal-value functions for bound-constrained convex programming problems. *Optimization*, 68(2-3): 691-711
- **Mitsos A, Chachuat B, and Barton PI. (2009).** McCormick-based relaxations of algorithms. *SIAM Journal on Optimization*, 20(2): 573–601.
- **Najman J, Bongratz D, Tsoukalas A, and Mitsos A (2017).** Erratum to: Multivariate McCormick relaxations. *Journal of Global Optimization*, 68: 219-225.
- **Scott JK,  Stuber MD, and Barton PI. (2011).** Generalized McCormick relaxations. *Journal of Global Optimization*, 51(4): 569–606.
- **Stuber MD, Scott JK, Barton PI (2015).** Convex and concave relaxations of implicit functions. *Optim. Methods Softw.*, 30(3): 424–460
- **Tsoukalas A and Mitsos A (2014).** Multivariate McCormick Relaxations. *Journal of Global Optimization*, 59:633–662.
- **Wechsung A, Scott JK, Watson HAJ, and Barton PI. (2015).** Reverse propagation of McCormick relaxations. *Journal of Global Optimization*, 63(1): 1-36.
