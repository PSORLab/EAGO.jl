<img src="https://github.com/PSORLab/EAGO.jl/blob/master/docs/src/full_Logo1.png" width="75%" height="75%">

# EAGO: Easy-Advanced Global Optimization
EAGO is an open-source development environment for **robust and global optimization** in Julia.

| **Documentation**                                                | **Linux/OS/Windows**                                                                     | **Persistent DOI**                                                                     |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://PSORLab.github.io/EAGO.jl/dev) | [![Build Status](https://github.com/PSORLab/EAGO.jl/workflows/CI/badge.svg?branch=master)](https://github.com/PSORLab/EAGO.jl/actions?query=workflow%3ACI) | [![DOI](https://zenodo.org/badge/108954118.svg)](https://zenodo.org/badge/latestdoi/108954118) |

| **Coverage** | **Chat** |
|:------------:|:------------:|
| [![codecov](https://codecov.io/gh/PSORLab/EAGO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PSORLab/EAGO.jl) | [![Join the chat at https://gitter.im/EAGODevelopment](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/EAGODevelopment/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link)

## EAGO's Optimizer Capabilities

EAGO is a deterministic global optimizer designed to address a wide variety of optimization problems by propagating McCormick relaxations along the factorable structure of each expression in the NLP. Most operators supported by modern AD packages (e.g. **+**, **sin**, **cosh**) are supported by EAGO and a number utilities for sanitizing native Julia code and generating relaxations on a wide variety of user-defined functions have been included. Currently, EAGO supports problems that have aprior variable bounds defined and have differentiable constraints. That is problems should be specified in the generic form below:

<p align="center">
<img src="https://github.com/PSORLab/EAGO.jl/blob/master/docs/readme/OptForm.svg" width="30%" height="30%">

## EAGO's Relaxations

For each nonlinear term EAGO makes use of factorable representation to construct bounds and relaxations. In the case of F = y(y-1)sin(y), a list is generated and rules for constructing McCormick relaxations are used to formulate relaxations in the original Y decision space<sup>1</sup>:

- *v*<sub>1</sub> = y
- *v*<sub>2</sub> = *v*<sub>1</sub> - 5
- *v*<sub>3</sub> = sin(*v*<sub>1</sub>)
- *v*<sub>4</sub> = *v*<sub>1</sub>*v*<sub>2</sub>
- *v*<sub>5</sub> = *v*<sub>4</sub>*v*<sub>3</sub>
- F = *v*<sub>5</sub>

<p align="center">
<img src="https://github.com/PSORLab/EAGO.jl/blob/master/docs/src/McCormick/Figure_1.png" width="60%" height="60%">

Either these original relaxations, differentiable McCormick relaxations<sup>2</sup>, or affine relaxations thereof can be used to construct relaxations of optimization problems useful in branch and bound routines for global optimization. Utilities are included to combine these with algorithms for relaxing implicit functions<sup>3</sup> and forward-reverse propagation of McCormick arithmetic<sup>4</sup>.

## Sample Usage

EAGO makes use of the JuMP modeling language to. Consider the familiar "process" problem instance<sup>5</sup>:

<p align="center">
<img src="https://github.com/PSORLab/EAGO.jl/blob/master/docs/readme/ProcessFormulation.svg" width="60%" height="60%">


This model can be formulated using JuMP code as

```julia
using JuMP, EAGO

m = Model(EAGO.Optimizer)

# Define bounded variables
xL = [10.0; 0.0; 0.0; 0.0; 0.0; 85.0; 90.0; 3.0; 1.2; 145.0]
xU = [2000.0; 16000.0; 120.0; 5000.0; 2000.0; 93.0; 95.0; 12.0; 4.0; 162.0]
@variable(m, xL[i] <= x[i=1:10] <= xU[i])

# Define nonlinear constraints
@NLconstraint(m, e1, -x[1]*(1.12+0.13167*x[8]-0.00667* (x[8])^2)+x[4] == 0.0)
@NLconstraint(m, e3, -0.001*x[4]*x[9]*x[6]/(98-x[6])+x[3] == 0.0)
@NLconstraint(m, e4, -(1.098*x[8]-0.038* (x[8])^2)-0.325*x[6]+x[7] == 57.425)
@NLconstraint(m, e5, -(x[2]+x[5])/x[1]+x[8] == 0.0)

# Define linear constraints
@constraint(m, e2, -x[1]+1.22*x[4]-x[5] == 0.0)
@constraint(m, e6, x[9]+0.222*x[10] == 35.82)
@constraint(m, e7, -3*x[7]+x[10] == -133.0)

# Define nonlinear objective
@NLobjective(m, Max, 0.063*x[4]*x[7] - 5.04*x[1] - 0.035*x[2] - 10*x[3] - 3.36*x[5])

# Solve the optimization problem
JuMP.optimize!(m)
```
<!--
$
\begin{align}
\max_{x \in X}{0.063x_4x_7 - 5.04x_1 - 0.035x_2 - 10x_3 - 3.36x_2} & \\
s.t. \; x_1(1.12 + 0.13167x_8-0.00667x_8^2) + x_4 &= 0 \\
-0.001x_4x_9x_6/(98-x_6)+x_3 &= 0 \\
-(1.098x_8 - 0.038x_8^2) - 0.325x_6 + x_7 &= 0 \\
-(x_2 + x_5)/x_1 + x_8 &= 0 \\
-x_1 + 1.22x_4 - x_5 &= 0 \\
x_9 + 0.222x_{10} - 35.82 &= 0 \\
-3.0x_7 + x_{10} + 133.0 &= 0
\end{align}
$
where $X= [10, 2000] \times[0, 16000] \times[0, 120] \times[0, 5000]$
                 $\qquad \times[0, 2000] \times[85, 93] \times [90,95] \times [3,12] \times [1.2, 4] \times [145, 162]$
-->

Special handling has been included for linear/quadratic functions defined using the `@constraint` macro in JuMP and these can generally be expected to perform better than specifying quadratic or linear terms with the `@NLconstraint` macro.

## A Cautionary Note on Global Optimization

As a global optimization platform, EAGO's solvers can be used to find solutions of general nonconvex problems with a guaranteed certificate of optimality. However, global solvers suffer from the curse of dimensionality and therefore their performance is  outstripped by convex solvers. For users interested in large-scale applications, be warned that problems generally larger than a few variables may prove challenging for certain types of global optimization problems.

## Package Capabilities

The EAGO package has numerous features: a solver accessible from JuMP/MathOptInterface, domain reduction routines, McCormick relaxations, and specialized non-convex semi-infinite program solvers. A full description of all EAGO features is available in the [**documentation website**](https://psorlab.github.io/EAGO.jl/dev/). A series of example have been provided in the form of Jupyter notebooks in the separate [**EAGO-notebooks**](https://github.com/PSORLab/EAGO-notebooks) repository

## Recent News

- 6/7/2020: [**EAGO v0.4.0 has been tagged**](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.4.0).
  - Support for new MOI/JuMP `RawParameter` input and a number of new attributes.
  - Separates McCormick and ReverseMcCormick libraries (now [McCormick.jl](https://github.com/PSORLab/McCormick.jl) and [ReverseMcCormick.jl](https://github.com/PSORLab/ReverseMcCormick.jl))
    from main package.  McCormick.jl is reexported.
  - Relaxation calculations now return NaN values on a domain violation.
  - Tolerance based validation of cuts has been added to generate numerically safe cuts.
  - Significantly simplify internal codebase for `EAGO.Optimizer` (no changes to API): fully decouples input problem specifications from the formulation used internally, stack only stores variables that are branched on, and a number of internal rearrangements to clearly delineate different routines.
  - Add problem classification preprocessing that throws to simpler routines if LP problem types are detected (enables future support for SOCP, MILP, MISOCP, and Convex forms).
  - Fix multiple bugs and add more transparent error codes.

- 06/17/2020: [**EAGO v0.4.1 has been tagged**](https://github.com/PSORLab/EAGO.jl/commit/9c1bcf024a19840a0ac49c8c6da13619a5f3845f#comments) Contains minor bug releases.
- 08/29/2020: [**EAGO v0.4.2 has been tagged**](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.4.2) Support for Julia v1.5.
- 11/18/2020: [**EAGO v0.5 has been tagged**](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.5.0)
  - Introduces the `register_eago_operators!(m::JuMP.Model)` which can be used
    to register all nonstandard nonlinear terms used in EAGO in any JuMP model.
  - Introduces `positive`, `negative`, `lower_bnd`, `upper_bnd`, and `bnd`
    functions which can be used to enforce bounds on intermediate terms in
    nonlinear expressions (EAGO.Optimizer only).
  - Adds envelopes: `abs2`, `sinpi`, `cospi`, `fma`, `cbrt`.
  - Adds envelopes and functions: `xlogx`
  - Adds envelopes of special functions: `erf`, `erfc`, `erfinv`, `erfcinv`
  - Adds envelopes of activation functions: `relu`, `gelu`, `elu`, `selu`, `swish1`,
                                             `sigmoid`, `softsign`, `softplus`,
                                             `bisigmoid`, `pentanh`, `leaky_relu`,
                                             `param_relu`.
  - Error messages in `sip_explicit` have been made more transparent.
  - Fixes some issues with documentation image rendering and links.
  - Drops appveyor CI and Travis CI in favor of Github Actions.


For a full list of EAGO release news, see click [**here**](https://github.com/PSORLab/EAGO.jl/releases)

## Installing EAGO

EAGO is registered Julia package. It can be installed using the Julia package manager.
From the Julia REPL, type ] to enter the Pkg REPL mode and run the following command

```julia
pkg> add EAGO
```

Currently, EAGO is tied to a 0.19+ or greater version of JuMP. This allows a replication
of some of the internal features shared by EAGO and JuMP's AD scheme aka
generation of Wergert Tapes pass evaluators between JuMP and EAGO etc.

```julia
pkg> add JuMP
```

EAGO v0.5.1 is the current tagged version and requires Julia 1.2+ for full functionality (however Julia 1.0+ versions support partial functionality). Use with version 1.5 is recommended as the majority of in-house testing has occurred using this version of Julia. The user is directed to the [**High-Performance Configuration**](https://psorlab.github.io/EAGO.jl/Optimizer/high_performance/)
for instructions on how to install a high performance version of EAGO (rather than the basic entirely open-source version).
If any issues are encountered when loading EAGO (or when using it), please submit an issue using the Github [**issue tracker**](https://github.com/PSORLab/EAGO.jl/issues).

## Bug reporting, support and feature requests

Please report bugs or feature requests by opening an issue using the Github [**issue tracker**](https://github.com/PSORLab/EAGO.jl/issues). All manners of feedback are encouraged.

## Current limitations
- Nonlinear handling assumes that box-constraints of nonlinear terms are available or can be inferred from bounds-tightening.
- Only currently supports continuous functions. Support for mixed-integer problems is forthcoming.

## Work In Progress
- Extensions for nonconvex dynamic global & robust optimization.
- Provide support for mixed-integer problems.
- Update EAGO to support nonsmooth problems (requires: a nonsmooth local nlp optimizer or lexiographic AD, support for relaxations is already included).
- Performance assessment of nonlinear (differentiable) relaxations and incorporation into main EAGO routine.
- Evaluation and incorporation of implicit relaxation routines in basic solver.

## Citing EAGO

Please cite the following paper when using EAGO. In plain tex form this is:
```
 M. E. Wilhelm & M. D. Stuber (2020) EAGO.jl: easy advanced global optimization in Julia,
 Optimization Methods and Software, DOI: 10.1080/10556788.2020.1786566
```

A corresponding bibtex entry text is given below and a corresponding .bib file is given in citation.bib.
```
@article{doi:10.1080/10556788.2020.1786566,
author = { M. E.   Wilhelm  and  M. D.   Stuber },
title = {EAGO.jl: easy advanced global optimization in Julia},
journal = {Optimization Methods and Software},
volume = {0},
number = {0},
pages = {1-26},
year  = {2020},
publisher = {Taylor & Francis},
doi = {10.1080/10556788.2020.1786566},
URL = {https://doi.org/10.1080/10556788.2020.1786566},
eprint = {https://doi.org/10.1080/10556788.2020.1786566}
}
```

## Related Packages

- [**ValidatedNumerics.jl**](https://github.com/JuliaIntervals/ValidatedNumerics.jl), a Julia library for validated interval calculations, including basic interval extensions, constraint programming, and interval contactors   
- [**MC++**](https://omega-icl.github.io/mcpp/): A mature McCormick relaxation package in C++ that also includes McCormick-Taylor, Chebyshev
Polyhedral and Ellipsoidal arithmetics.

## References
1. A. Mitsos, B. Chachuat, and P. I. Barton. **McCormick-based relaxations of algorithms.** *SIAM Journal on Optimization*, 20(2):573–601, 2009.
2. K.A. Khan, HAJ Watson, P.I. Barton. **Differentiable McCormick relaxations.** *Journal of Global Optimization*, 67(4):687-729 (2017).
3. Stuber, M.D., Scott, J.K., Barton, P.I.: **Convex and concave relaxations of implicit functions.** *Optim. Methods Softw.* 30(3), 424–460 (2015)
4. A., Wechsung JK Scott, HAJ Watson, and PI Barton. **Reverse propagation of McCormick relaxations.** *Journal of Global Optimization* 63(1):1-36 (2015).
5. Bracken, Jerome and McCormick, Garth P. **Selected Applications of Nonlinear Programming**, John Wiley and Sons, New York, 1968.
