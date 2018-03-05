# EAGO
A development environment for robust and global optimization.

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
