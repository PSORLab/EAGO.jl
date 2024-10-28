<img src="https://github.com/PSORLab/EAGO.jl/blob/master/docs/src/assets/logo.png" width="75%" height="75%">

# EAGO - Easy Advanced Global Optimization

EAGO is an open-source development environment for **robust and global optimization** in Julia.

| **PSOR Lab** | **Build Status**                                                                                |
|:------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/Developed_by-PSOR_Lab-342674)](https://psor.uconn.edu/) | [![Build Status](https://github.com/PSORLab/EAGO.jl/workflows/CI/badge.svg?branch=master)](https://github.com/PSORLab/EAGO.jl/actions?query=workflow%3ACI) [![codecov](https://codecov.io/gh/PSORLab/EAGO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PSORLab/EAGO.jl)|

| **Documentation**                                                 | **Persistent DOI**                                                                              |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://PSORLab.github.io/EAGO.jl/stable) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://PSORLab.github.io/EAGO.jl/dev) | [![DOI](https://zenodo.org/badge/108954118.svg)](https://zenodo.org/badge/latestdoi/108954118) |

| **Chat**     |
|:------------:|
| [![Join the chat at https://gitter.im/EAGODevelopment](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/EAGODevelopment/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link)

## EAGO's Optimizer Capabilities

EAGO is a deterministic global optimizer designed to address a wide variety of optimization problems, emphasizing nonlinear programs (NLPs), by propagating McCormick relaxations along the factorable structure of each expression in the NLP. Most operators supported by modern automatic differentiation (AD) packages (e.g., `+`, `sin`, `cosh`) are supported by EAGO and a number utilities for sanitizing native Julia code and generating relaxations on a wide variety of user-defined functions have been included. Currently, EAGO supports problems that have a priori variable bounds defined and have differentiable constraints. That is, problems should be specified in the generic form below:

$$
\begin{align*}
f^{\*} = & \min_{\mathbf y \in Y \subset \mathbb R^{n_{y}}} f(\mathbf y) \\
{\rm s.t.} \\;\\; & \mathbf h(\mathbf y) = \mathbf 0 \\
& \mathbf g(\mathbf y) \leq \mathbf 0 \\
& Y = [\mathbf y^{\mathbf L}, \mathbf y^{\mathbf U}] \in \mathbb{IR}^{n} \\
& \qquad \mathbf y^{\mathbf L}, \mathbf y^{\mathbf U} \in \mathbb R^{n}
\end{align*}
$$

## EAGO's Relaxations

For each nonlinear term, EAGO makes use of factorable representations to construct bounds and relaxations. In the case of $f(x) = x (x - 5) \sin(x)$, a list is generated and rules for constructing McCormick relaxations are used to formulate relaxations in the original decision space, $X$ [[1](#references)]:

- $v_{1} = x$
- $v_{2} = v_{1} - 5$
- $v_{3} = \sin(v_{1})$
- $v_{4} = v_{1} v_{2}$
- $v_{5} = v_{4} v_{3}$
- $f(x) = v_{5}$

<p align="center">
<img src="https://github.com/PSORLab/EAGO.jl/blob/master/docs/src/mccormick/Figure_1.png" width="60%" height="60%">

Either these original relaxations, differentiable McCormick relaxations [[2](#references)], or affine relaxations thereof can be used to construct relaxations of optimization problems useful in branch and bound routines for global optimization. Utilities are included to combine these with algorithms for relaxing implicit functions [[3](#references)] and forward-reverse propagation of McCormick arithmetic [[4](#references)].

## Sample Usage

EAGO makes use of the JuMP algebraic modeling language to improve the user's experience in setting up optimization models. Consider the familiar "process" problem instance [[5](#references)]:

$$
\begin{align*}
& \max_{\mathbf x \in X} 0.063 x_{4} x_{7} - 5.04 x_{1} - 0.035 x_{2} - 10 x_{3} - 3.36 x_{2} \\
{\rm s.t.} \\;\\; & x_{1} (1.12 + 0.13167 x_{8} - 0.00667 x_{8}^{2}) + x_{4} = 0 \\
& -0.001 x_{4} x_{9} x_{6} / (98 - x_{6}) + x_{3} = 0 \\
& -(1.098 x_{8} - 0.038 x_{8}^{2}) - 0.325 x_{6} + x_{7} = 0 \\
& -(x_{2} + x_{5}) / x_{1} + x_{8} = 0 \\
& -x_{1} + 1.22 x_{4} - x_{5} = 0 \\
& x_{9} + 0.222 x_{10} - 35.82 = 0 \\
& -3.0 x_{7} + x_{10} + 133.0 = 0 \\
& X = [10, 2000] \times [0, 16000] \times [0, 120] \times [0, 5000] \\
& \qquad \times [0, 2000] \times [85, 93] \times [90,9 5] \times [3, 12] \times [1.2, 4] \times [145, 162]
\end{align*}
$$


This model can be formulated using JuMP code as:

```julia
using JuMP, EAGO

model = Model(EAGO.Optimizer)

# Define bounded variables
xL = [10.0; 0.0; 0.0; 0.0; 0.0; 85.0; 90.0; 3.0; 1.2; 145.0]
xU = [2000.0; 16000.0; 120.0; 5000.0; 2000.0; 93.0; 95.0; 12.0; 4.0; 162.0]
@variable(model, xL[i] <= x[i=1:10] <= xU[i])

# Define constraints
@constraint(model, e1, -x[1]*(1.12 + 0.13167*x[8] - 0.00667*(x[8])^2) + x[4] == 0.0)
@constraint(model, e2, -x[1] + 1.22*x[4] - x[5] == 0.0)
@constraint(model, e3, -0.001*x[4]*x[9]*x[6]/(98.0 - x[6]) + x[3] == 0.0)
@constraint(model, e4, -(1.098*x[8] - 0.038*(x[8])^2) - 0.325*x[6] + x[7] == 57.425)
@constraint(model, e5, -(x[2] + x[5])/x[1] + x[8] == 0.0)
@constraint(model, e6, x[9] + 0.222*x[10] == 35.82)
@constraint(model, e7, -3.0*x[7] + x[10] == -133.0)

# Define objective
@objective(model, Max, 0.063*x[4]*x[7] - 5.04*x[1] - 0.035*x[2] - 10*x[3] - 3.36*x[5])

# Solve the optimization problem
JuMP.optimize!(model)
```

## A Cautionary Note on Global Optimization

As a global optimization platform, EAGO's solvers can be used to find solutions of general nonconvex problems with a guaranteed certificate of optimality. However, global solvers suffer from the curse of dimensionality and therefore their performance is outstripped by convex/local solvers. For users interested in large-scale applications, be warned that problems generally larger than a few variables may prove challenging for certain types of global optimization problems.

## Package Capabilities

The EAGO package has numerous features: a solver accessible from JuMP/MathOptInterface (MOI), domain reduction routines, McCormick relaxations, and specialized nonconvex semi-infinite program solvers. A full description of all EAGO features is available on the [documentation website](https://psorlab.github.io/EAGO.jl/dev/). A series of example have been provided in the form of Jupyter Notebooks in the separate [EAGO-notebooks](https://github.com/PSORLab/EAGO-notebooks) repository.

## Recent News

### [v0.8.2](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.8.2) (October 27, 2024)
- Added support for `MOI.ScalarNonlinearFunction`.
  - Users can now define all constraints using `@constraint` instead of needing to use `@NLconstraint`. This applies to `@objective` as well.
- Added support for variable names.
- Updated display.

### [v0.8.1](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.8.1) (June 15, 2023)

- Resolved an issue where integer and binary variables would sometimes throw a `MathOptInterface.UpperBoundAlreadySet` error.
- Added the function `unbounded_check!` which warns users if they are missing variable bounds and sets them to +/- 1E10 by default.
  - Added an EAGO parameter `unbounded_check` which defaults to `true` and enables `unbounded_check!`.
- Bumped requirement for PrettyTables.jl to v2+ to accommodate the latest version of DataFrames.jl.

### [v0.8.0](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.8.0) (June 12, 2023)

- Updated EAGO for compatibility with the nonlinear expression API changes introduced in JuMP v1.2: https://discourse.julialang.org/t/ann-upcoming-refactoring-of-jumps-nonlinear-api/83052.
  - EAGO now uses the `MOI.Nonlinear` submodule instead of `JuMP._Derivatives`.
  - Models, nodes, expressions, constraints, and operators are now compatible with MOI.
- Added logic and comparison operators to `EAGO.OperatorRegistry`.

For a full list of EAGO release news, click [here](https://github.com/PSORLab/EAGO.jl/blob/master/NEWS.md).

## Installing EAGO

EAGO is a registered Julia package and it can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the Package manager (Pkg) mode and run the following command:

```jldoctest
pkg> add EAGO
```

Currently, EAGO is compatible with version 1.12 of JuMP. This allows a replication of some of the internal features shared by EAGO and JuMP's AD scheme, e.g., generation of Wengert Tapes, passing evaluators between JuMP and EAGO, etc.

```jldoctest
pkg> add JuMP
```

EAGO v0.8.1 is the current tagged version and requires Julia 1.6+ for full functionality (however Julia 1.0+ versions support partial functionality). Use with version 1.8 is recommended as the majority of in-house testing has occurred using this version of Julia. The user is directed to the [High-Performance Configuration](https://psorlab.github.io/EAGO.jl/optimizer/high_performance/) for instructions on how to install a high performance version of EAGO (rather than the basic entirely open-source version). If any issues are encountered when loading EAGO (or when using it), please submit an issue using the GitHub [issue tracker](https://github.com/PSORLab/EAGO.jl/issues).

## Bug Reporting, Support, and Feature Requests

Please report bugs or feature requests by opening an issue using the GitHub [issue tracker](https://github.com/PSORLab/EAGO.jl/issues). All manners of feedback are encouraged.

## Current Limitations

- Nonlinear handling assumes that box-constraints of nonlinear terms are available or can be inferred from bounds-tightening.
- Only currently supports continuous functions. Support for mixed-integer problems is forthcoming.

## Work in Progress

- Extensions for nonconvex dynamic global & robust optimization.
- Provide support for mixed-integer problems.
- Update EAGO to support nonsmooth problems (requires: a nonsmooth local nlp optimizer or lexicographic AD, support for relaxations is already included).
- Performance assessment of nonlinear (differentiable) relaxations and incorporation into main EAGO routine.
- Evaluation and incorporation of implicit relaxation routines in basic solver.

## Citing EAGO

Please cite the following paper when using EAGO. In plain text form this is:

```
Wilhelm, M.E. and Stuber, M.D. EAGO.jl: easy advanced global optimization in Julia.
Optimization Methods and Software. 37(2): 425-450 (2022). DOI: 10.1080/10556788.2020.1786566
```

A BibTeX entry is given below and a corresponding .bib file is given in citation.bib.

```bibtex
@article{doi:10.1080/10556788.2020.1786566,
author = {Wilhelm, M.E. and Stuber, M.D.},
title = {EAGO.jl: easy advanced global optimization in Julia},
journal = {Optimization Methods and Software},
volume = {37},
number = {2},
pages = {425-450},
year  = {2022},
publisher = {Taylor & Francis},
doi = {10.1080/10556788.2020.1786566},
URL = {https://doi.org/10.1080/10556788.2020.1786566},
eprint = {https://doi.org/10.1080/10556788.2020.1786566}
}
```

## Related Packages

- [ValidatedNumerics.jl](https://github.com/JuliaIntervals/ValidatedNumerics.jl): A Julia library for validated interval calculations, including basic interval extensions, constraint programming, and interval contractors
- [MAiNGO](https://avt-svt.pages.rwth-aachen.de/public/maingo/): An open-source mixed-integer nonlinear programming package in C++ that utilizes MC++ for relaxations
- [MC++](https://github.com/coin-or/MCpp): A mature McCormick relaxation package in C++ that also includes McCormick-Taylor, Chebyshev Polyhedral, and Ellipsoidal arithmetics

## References

1. Mitsos, A., Chachuat, B., and Barton, P.I. **McCormick-based relaxations of algorithms.** *SIAM Journal on Optimization*. 20(2): 573–601 (2009).
2. Khan, K.A., Watson, H.A.J., and Barton, P.I. **Differentiable McCormick relaxations.** *Journal of Global Optimization*. 67(4): 687-729 (2017).
3. Stuber, M.D., Scott, J.K., and Barton, P.I.: **Convex and concave relaxations of implicit functions.** *Optimization Methods and Software* 30(3): 424–460 (2015).
4. Wechsung, A., Scott, J.K., Watson, H.A.J., and Barton, P.I. **Reverse propagation of McCormick relaxations.** *Journal of Global Optimization* 63(1): 1-36 (2015).
5. Bracken, J., and McCormick, G.P. *Selected Applications of Nonlinear Programming.* John Wiley and Sons, New York (1968).