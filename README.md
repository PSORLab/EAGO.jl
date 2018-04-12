# EAGO
EAGO is an open-source development environment for **robust and global optimization** in Julia. 

[![Build Status](https://travis-ci.org/MatthewStuber/EAGO.jl.svg?branch=master)](https://travis-ci.org/MatthewStuber/EAGO.jl)
[![Coverage Status](https://coveralls.io/repos/github/MatthewStuber/EAGO.jl/badge.svg?branch=master)](https://coveralls.io/github/MatthewStuber/EAGO.jl?branch=master)
[![codecov.io](http://codecov.io/github/MatthewStuber/EAGO.jl/coverage.svg?branch=master)](http://codecov.io/github/MatthewStuber/EAGO.jl?branch=master)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://MatthewStuber.github.io/EAGO.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://MatthewStuber.github.io/EAGO.jl/latest)

## Global Optimization

As a global optimization platform, EAGO's solvers can be used to find solutions of general nonconvex problems with a guaranteed certificate of optimality. However, global solvers suffer from the curse of dimensionality and therefore their performance is  outstripped by convex solvers. For users interested in large-scale applications, be warned that problems generally larger than a few variables may prove challenging for certain types of global optimization problems. 

## So what differentiates the EAGO's global solver from currently available offerings (BARON, Antigone, etc.)?

- **McCormick Relaxations**:
First, the global **EAGONLP_Solver** supports the use of McCormick relaxations to construct the convex lower-bounding problem. These relaxations have some appealing properties (2nd order convergence, low computational cost) and can be evaluated directly via operator overloading schemes. As such, the solver supports standard arithmetic operators, trignometric operators (**sin**, **cos**, **tan**, **asin**, **acos**, **atan**), hyperbolic operators (**sinh**, **cosh**, ...), and a set of nonsmooth operators (**abs**, **min**, **max**, **step**, **sign**). It can easily be extended to additional univariate functions provided that explicit forms of the convex/concave relaxations and interval bounds are known.

- **Relaxation of Implicit Functions** As the evaluation of relaxations proceeds in an operator overloading fashion, functions implicitly defined via equality constraints h(x,p) = 0 can be relaxed in the p dimension only. This can substantially reduce the problem size and provide a significant performance benefit relative to existing solvers for problems with large embedded implicit functions.

- **Expression Handling**:
Second, the majority of the code is written directly in Julia allow for easy manipulation of expressions using an AST. This allows for facile identification of implicit functions and

## What's unique about the robust optimization capabilities?

- **Nonconvex Semi-infinite Programming**: An implementation of this is available in GAMS using BARON. However, this is the first nonconvex robust optimizer available in Julia. Second, the majority of the code is written directly in Julia allowing for easy manipulation of expressions using an abstract syntax tree (AST). This allows for facile identification of implicit functions.

## What's unique about the robust optimization capabilities?

- **Nonconvex Semi-infinite Programming**: An implementation of this in GAMS is available using the BARON solver. However, this is the first nonconvex robust optimizer available in Julia. 
- **Semi-Infinite Equality Constraints**: Nonconvex semi-infinite programming with equality constraints represents an extremely challenging problem. In principle, the problem can be decomposed into a series of semi-infinite constraints and handled with a restriction of the right-hand side technique. However, the nonexistence of Slater points in the explicit formulation is likely and  thus the above algorithm will lack a guarantee of convergence. By solving the problem using the reduced-space implicit approach, Slater points can be shown to trivially exist and we recover the convergence guarantee.

## Capabilities as a development platform

- **McCormick Relaxation Library**
  * *Standard McCormick Relaxations*: 
  * *Multivariate McCormick Relaxations*: 
  * *Differentiable McCormick Relaxations*: 
  * *Subgradient Interval Bound Tightening*:
  * *Implicit Fixed-Point Relaxations*:
- **Domain Reduction Techniques**:
  * *Forward-Backward Interval Contractor*: Runs forward and reverse interval contractor propagation.
  * *Standard Range Reduction*: Range Reduction using McCormick relaxations to compute linear relaxtions of function and bounds.
  * *Implicit Range Reduction*: Range Reduction on decision space only using McCormick relaxations to compute linear relaxations of implicit bounds.
  * *Probing on Variable Bounds:* Uses McCormick Relaxation to compute linear relaxations and subproblems.
- **Native Lower Bounding Problem**:
  * *Interval Arithmetic*: Calculated using IntervalArithmetic.jl.
  * *LP of McCormick Relaxation*: Supports use of any LP solver using the MathProgBases interface.
  * *Convex NLP of Differentiable McCormick Relaxation*: Currently, supports using either SNOPT or Ipopt.
  * *Convex NLP AlphaBB*: Currently, solved using Ipopt.
  * *LP of McCormick Relaxation (with Implicit Bounding Routine)*: Supports use of any LP solver using the MathProgBases interface.
  * *Convex NLP of Differentiable McCormick Relaxation (with Implicit Bounding Routine)*: Currently, supports using either SNOPT or Ipopt.
- **Native Upper Bounding Problem**:
  * *Interval Arithmetic*: Calculated using IntervalArithmetic.jl.
  * *LP of McCormick Relaxation*: Supports use of any LP solver using the MathProgBases interface.
  * *Local NLP Problem*: Currently, supports using either SNOPT or Ipopt.
- **Optimization Interface**
  * *MathProgBase support*
  * *JuMP interface support*
  * *Benchmarking through JuMP interface*

## Description
This package contains a global solver accessible from JuMP and subroutines for solving semi-infinite programs. The EAGO solver is usable in the JuMP enviroment.

## Subpackages
While the EAGO package is independent of the below subpackages. Their functionality can be accessed independently using the below packages:
- [EAGOBranchBound](https://github.com/MatthewStuber/EAGOSmoothMcCormick): A fully customizable branch and bound library. Supports standard bisection, node selection, visualization, and timing.
- [EAGOParametricInterval](https://github.com/MatthewStuber/EAGOParametricInterval): Contains a series of parametric interval contractors which provide interval bounds on implicit functions. Also contains various tests for the existence and uniqueness of implicit functions inside an interval bound and a generalized bisection routine for partitioning an interval box into sub-boxes that contain a unique implicit function.
- [EAGOSmoothMcCormick](https://github.com/MatthewStuber/EAGOSmoothMcCormick): A library that supports not only smooth McCormick relaxations but allows for calculation of gradients using ForwardDiff.jl/ReverseDiff.jl.
- [EAGOSmoothMcCormickGrad](https://github.com/MatthewStuber/EAGOSmoothMcCormickGrad): A library of smooth McCormick relaxations with an embedded vector storage used for subgradient and gradient calculations. Supports standard, nonsmooth multivariate, and smooth McCormick relaxations. Supports interval bound tightening via using subgradients. Also, contains subroutines for generating relaxations of implicit functions using McCormick relaxations.
- [EAGODomainReduction](https://github.com/MatthewStuber/EAGODomainReduction): A library for domain reduction techniques used in global optimization. Currently includes: forward-reverse interval constraint propagation, standard range reduction (using McCormick relaxations), duality-based bound tightening, and probing techniques. Further updates planned.

## News

Trying to switch to the latest version of ForwardDiff? See our [upgrade guide](http://www.juliadiff.org/ForwardDiff.jl/stable/user/upgrade.html) for details regarding user-facing changes between releases.

- 4/2018: [EAGO v0.1.0 has been tagged](https://github.com/JuliaLang/METADATA.jl/pull/14218). We're currently moving all the EAGO subpackages `EAGODomainReduction` into the main `EAGO.jl` package. Once this is finalized and the documentations is complete we'll take EAGO v0.2.0. 

## Further Work
We're currently improving the EAGO platform in multiple ways:
- Adding support for MathOptInterface.
- Incorporate subroutines to automatically detect implicit functions and reformulate explicitly defined problems.
- Support for MINLP problem types as well as mixed-integer nonconvex SIPs.
- A GUI interface to handle flowsheeting type problems and data-fitting problems.
- Optimizing default parameter settings.

## How to Contribute

- Feedback on applications is always appreciated.
- Add new relaxation

## Related Packages

- **ValidatedNumerics.jl**: Provides the interval arithmetic backbone for the SmoothMcCormickGrad package. It allow contains a series of subroutine that allow for locating roots and minimal of unconstrained functions via interval arithmetic. EAGO allows provides an implementation for interval optimization that can be called through JuMP and MathProgBase. In general, McCormick relaxations exhibit superior convergence properties which may help mitigate clustering but will not guarantee correctly rounded bounds.    
- JuMP.jl

## References
See subpackage pages.
