# EAGO
A is an open-source development environment for robust and global optimization in Julia. Included in EAGO 

## Global Optimization

As a global optimization platform, EAGO's solvers can be used to find the solution to general nonconvex problems with a guaranteed certificate of optimality. However, the performance of global solvers is generally outstripped by convex solvers so for users interested in large applications be warned that problems generally larger a few variables may prove challenging. For certain types of global optimization problems 

## So what differentiates the EAGO's global solver from currently available offerings (BARON, Antigone, etc.)?

- **McCormick Relaxations**:
First, the global **EAGONLP_Solver** supports the use of McCormick relaxations to construct the convex lower bounding problem. These relaxations have some appealing properties (2nd order convergence, low computational cost) and can be evaluated directly via overloading schemes. As such, the solver supports standard arithmetic operators, trignometric operators (sin,cos,tan,asin,acos,atan), hyperbolic operators (sinh,cosh,...), and a set of nonsmooth operators (abs,min,max,step,sign). It can easily be extended to additional univariant functions provided that explicit forms of the convex, concave relaxations and interval bounds are known.

- **Relaxation of Implicit Functions** As the evaluation of relaxations proceeds in an operator overloading fashion. Functions implicitly defined by equality constraints h(x,p) = 0 can be relaxed in the p dimension only. This can substantially reduces the problem size and provide a significant perform benefit relative to existing solvers for with large embedded implicit functions.

- **Expression Handling**:
Second, the majority of the code is written directly in Julia allow for easy manipulation of expressions using an AST. This allows for facile identification of implicit functions and 

## What's unique about the robust optimization capabilities?

- **Nonconvex Semi-infinite Programming**: An implementation of this by is available BARON. However, this is the first nonconvex robust optimizer available in Julia. 

- **Semi-Infinite Equality Constraints**: Nonconvex semi-infinite programming with equality constraints represents an extremely challenging problem. In principle, the problem can be decomposed into a series of semi-infinite constraints and handled with a restriction of the rhs technique. However, the nonexistence of Slater points in the explicit formulation is likely and the thus the above algorithm will lack a guarantee of convergence. By solve, the problem in the implicit space, p, Slater points can be shown to trivially exist and we recover the guarantee of convergence.

## Capabilities as a development platform

- **Native Lower and Upper Bounding Problems**
- **Domain Reduction Techniques**
- **McCormick Relaxation Library**
- **Optimization Interface**
- **Benchmarking**

## Description
This package contains a global solver accessible from JuMP and subroutines for solving semi-infinite programs. The EAGO solver is usable in the JuMP enviroment.

## Subpackages
The currently available subpackages are:
- [EAGOBranchBound](https://github.com/MatthewStuber/EAGOSmoothMcCormick): A fully customizable branch and bound library. Supports standard bisection, node selection, visualization, and timing.
- [EAGOParametricInterval](https://github.com/MatthewStuber/EAGOParametricInterval): Contains a series of parametric interval contractors which provide interval bounds on implicit functions. Also contains various tests for the existence and uniqueness of implicit functions inside an interval bound and a generalized bisection routine for partitioning an interval box into sub-boxes that contain a unique implicit function. 
- [EAGOSmoothMcCormick](https://github.com/MatthewStuber/EAGOSmoothMcCormick): A library that supports only of smooth McCormick relaxations but allows for calculation of gradients using ForwardDiff.jl/ReverseDiff.jl.
- [EAGOSmoothMcCormickGrad](https://github.com/MatthewStuber/EAGOSmoothMcCormickGrad): A library of smooth McCormick relaxations with an imbedded vector storage used for subgradient and gradient calculations. Supports standard, nonsmooth multivariant, and smooth McCormick relaxations. Supports interval bound tightening via using subgradients. Also, contains subroutines for generating relaxations of implicit functions using McCormick relaxations.
- [EAGODomainReduction](https://github.com/MatthewStuber/EAGODomainReduction): A library for domain reduction techniques used in global optimization. Currently includes: forward-reverse interval constraint propagation, standard range reduction (using McCormick relaxations), duality-based bound tightening, and probing techniques. Further updates planned.

## Further Work
We're currently improving the EAGO platform in multiple ways:
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

