<img src="https://github.com/PSORLab/EAGO.jl/blob/master/docs/src/full_Logo1.png" width="75%" height="75%">

# EAGO: Easy-Advanced Global Optimization
EAGO is an open-source development environment for **robust and global optimization** in Julia.

| **Documentation**                                                               | **Linux/OS**                                                                     | **Windows**                                                                     |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://PSORLab.github.io/EAGO.jl/stable) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://PSORLab.github.io/EAGO.jl/latest) | [![Build Status](https://travis-ci.org/PSORLab/EAGO.jl.svg?branch=master)](https://travis-ci.org/PSORLab/EAGO.jl) | [![Build status](https://ci.appveyor.com/api/projects/status/oc8ql7kdugc19skf?svg=true)](https://ci.appveyor.com/project/mewilhel/eago-jl) |

| **Coverage** | **Chat** |
|:------------:|:------------:|
|[![Coverage Status](https://coveralls.io/repos/github/PSORLab/EAGO.jl/badge.svg?branch=master)](https://coveralls.io/github/PSORLab/EAGO.jl?branch=master) [![codecov](https://codecov.io/gh/PSORLab/EAGO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PSORLab/EAGO.jl) | [![Join the chat at https://gitter.im/EAGODevelopment](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/EAGODevelopment/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link)

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

m = Model(with_optimizer(EAGO.Optimizer))

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

Special handling has been included for linear/quadratic functions defined using the `@constraint` macro in JuMP and these can generally be expected to perform better than specify quadratic or linear terms with the `@NLconstraint` macro.

## A Cautionary Note on Global Optimization

As a global optimization platform, EAGO's solvers can be used to find solutions of general nonconvex problems with a guaranteed certificate of optimality. However, global solvers suffer from the curse of dimensionality and therefore their performance is  outstripped by convex solvers. For users interested in large-scale applications, be warned that problems generally larger than a few variables may prove challenging for certain types of global optimization problems.

## Package Capabilities

The EAGO package has numerous features: a solver accessible from JuMP/MathOptInterface, domain reduction routines, McCormick relaxations, and specialized non-convex semi-infinite program solvers. A full description of all EAGO features is available in the [**documentation website**](https://psorlab.github.io/EAGO.jl/stable/). A series of example have been provided in the form of Jupyter notebooks in the separate [**EAGO-notebooks**](https://github.com/PSORLab/EAGO-notebooks) repository

## Recent News

- 6/14/2019: [**EAGO v0.2.0 has been tagged**](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.2.0). This update creates a number of breaking changes to the EAGO API. Please review the use cases provided in the documentation to update examples.
  - Updated to support Julia 1.0+, MathOptInterface (MOI), and MOI construction of subproblems.
  - Additional domain reduction routines available.
  - Support for specialized handling of linear and quadratic terms.
  - Significant performance improvements due to pre-allocation of Wengert tapes and MOI support.
  - A more intuitive API for McCormick relaxation construction.
- 7/5/2019: [**EAGO v0.2.1 has been tagged**](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.2.1). This contains fixes for a few minor issues.
  - Bug fix for explicit SIP solving routine that occurred for uncertainty sets of dimension greater than 1.
  - Bug fix for Max objective sense.
<!--
- 11/1/2019: [**EAGO v0.3.0 has been tagged**](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.3.0): This update is intended to be the last to create a large number of breaking changes to the EAGO API. Please review the use cases provided in the documentation to update examples.
  - A number of performance improvements have been made to the underlying McCormick relaxation library.
  - The optimizer used to construct relaxations is now modified in place.
  - All subproblem storage has been moved to the Optimizer object and storage types (e.g. LowerInfo) have been removed.
  - A MinMax heap structure is now used to store nodes.
  - Speed and aesthetics for logging and printing utilities have been updated.
  - Subroutines are now customized by creating a subtype of 'ExtensionType' and defining subroutines which dispatch on this new structure.
  - Parametric interval methods and the Implicit optimizer have been move to a separate package (to be tagged shortly.)
  - JIT compilation time has been reduced substantially.
  - Support for silent tag and time limits.
-->

For a full list of EAGO release news, see click [**here**](https://github.com/PSORLab/EAGO.jl/releases)

## Bug reporting and support

Please report bugs or feature requests by opening an issue using the Github [**issue tracker**](https://github.com/PSORLab/EAGO.jl/issues). All manners of feedback are encouraged.

## Work In Progress
- Extensions for nonconvex dynamic global & robust optimization.
- Provide support for mixed-integer problems.
- Update EAGO to support nonsmooth problems (requires: a nonsmooth local nlp optimizer or lexiographic AD, support for relaxations is already included).
- Performance assessment of nonlinear (differentiable) relaxations and incorporation into main EAGO routine.
- Evaluation and incorporation of implicit relaxation routines in basic solver.

## Citing EAGO

```
Wilhelm, Matthew; Stuber, Matthew (October 2017) Easy Advanced Global
Optimization (EAGO): An Open-Source Platform for Robust and Global Optimization
in Julia. Presented at the AIChE Annual Meeting in Minneapolis, MN.
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
