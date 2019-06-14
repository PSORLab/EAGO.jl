<img src="https://github.com/PSORLab/EAGO.jl/blob/master/docs/src/full_Logo1.png" width="75%" height="75%">

# EAGO: Easy-Advanced Global Optimization
EAGO is an open-source development environment for **robust and global optimization** in Julia.


| **Documentation**                                                               | **PackageEvaluator**                                            | **Linux/OS**                                                                     | **Windows 32/64** |
|:-------------------------------------------------------------------------------:|:---------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://PSORLab.github.io/EAGO.jl/stable) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://PSORLab.github.io/EAGO.jl/latest) | Coming Soon | [![Build Status](https://travis-ci.org/PSORLab/EAGO.jl.svg?branch=master)](https://travis-ci.org/PSORLab/EAGO.jl) | [![Build status](https://ci.appveyor.com/api/projects/status/edwwjgvbkdsqcr1t?svg=true)](https://ci.appveyor.com/project/MatthewStuber/eago-jl)

| **Coverage** | **Chat** |
|:------------:|:------------:|
|[![Coverage Status](https://coveralls.io/repos/github/PSORLab/EAGO.jl/badge.svg?branch=master)](https://coveralls.io/github/PSORLab/EAGO.jl?branch=master) [![codecov](https://codecov.io/gh/PSORLab/EAGO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PSORLab/EAGO.jl) | [![Join the chat at https://gitter.im/EAGODevelopment](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/EAGODevelopment/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link)

## Global Optimization

As a global optimization platform, EAGO's solvers can be used to find solutions of general nonconvex problems with a guaranteed certificate of optimality. However, global solvers suffer from the curse of dimensionality and therefore their performance is  outstripped by convex solvers. For users interested in large-scale applications, be warned that problems generally larger than a few variables may prove challenging for certain types of global optimization problems.

## Package Capabilities

The EAGO package has numerous features: a solver accessible from JuMP/MathOptInterface, domain reduction routines, McCormick relaxations, and specialized non-convex semi-infinite program solvers. A full description of all EAGO features is available in the documentation website.

## News

- 4/12/2018: Initial release of combined EAGO packages.
- 6/20/2018: [EAGO v0.1.2 has been tagged](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.1.2). Significant speed and functionality updates.
- 6/14/2018: [EAGO v0.2.0 has been tagged](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.2.0). This update creates a number of breaking changes to the EAGO API. Please review the use cases provided in the documentation to update examples.
  - Updated to support Julia 1.0, MathOptInterface (MOI), and MOI construction of subproblems.
  - Additional domain reduction routines available.
  - Support for specialized handling of linear and quadratic terms.
  - Significant performance improvements due to pre-allocation of Wengert tapes and MOI support.
  - A more intuitive API for McCormick relaxation construction.

## Related Packages

- [**ValidatedNumerics.jl**](https://github.com/JuliaIntervals/ValidatedNumerics.jl), a Julia library for validated interval calculations, including basic interval extensions, constraint programming, and interval contactors   
- [**MC++**](https://omega-icl.github.io/mcpp/): A mature McCormick relaxation package in C++ that also includes McCormick-Taylor, Chebyshev
Polyhedral and Ellipsoidal arithmetics.
