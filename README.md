# EAGO: Easy-Advanced Global Optimization
EAGO is an open-source development environment for **robust and global optimization** in Julia.


| **Documentation**                                                               | **PackageEvaluator**                                            | **Linux/OS**                                                                     | **Windows 32/64** |
|:-------------------------------------------------------------------------------:|:---------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://psorlab.github.io/EAGO.jl/) | Coming Soon | [![Build Status](https://travis-ci.org/PSORLab/EAGO.jl.svg?branch=master)](https://travis-ci.org/PSORLab/EAGO.jl) | [![Build status](https://ci.appveyor.com/api/projects/status/edwwjgvbkdsqcr1t?svg=true)](https://ci.appveyor.com/project/MatthewStuber/eago-jl)

| **Coverage** | **Chat** |
|:------------:|:------------:|
|[![Coverage Status](https://coveralls.io/repos/github/PSORLab/EAGO.jl/badge.svg?branch=master)](https://coveralls.io/github/PSORLab/EAGO.jl?branch=master) [![codecov](https://codecov.io/gh/PSORLab/EAGO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PSORLab/EAGO.jl) | [![Join the chat at https://gitter.im/EAGODevelopment](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/EAGODevelopment/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link)

## Global Optimization

As a global optimization platform, EAGO's solvers can be used to find solutions of general nonconvex problems with a guaranteed certificate of optimality. However, global solvers suffer from the curse of dimensionality and therefore their performance is  outstripped by convex solvers. For users interested in large-scale applications, be warned that problems generally larger than a few variables may prove challenging for certain types of global optimization problems.

## Package Capabilities

The EAGO package has numerous features: a solver accessible from JuMP/MathProgBase, domain reduction routines, McCormick relaxations, and specialized non-convex semi-infinite program solvers. A full description of all EAGO features is available in the documentation website: https://psorlab.github.io/EAGO.jl/.

## News

- 4/12/2018: [EAGO v0.1.0 has been tagged](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.1.0). Initial release of combined EAGO packages.
- 6/08/2018: [EAGO v0.1.1 has been tagged](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.1.1). Significant speed and functionality updates.

## Related Packages

- [**ValidatedNumerics.jl**](https://github.com/JuliaIntervals/ValidatedNumerics.jl), a Julia library for validated interval calculations, including basic interval extensions, constraint programming, and interval contactors   
- [**MC++**](https://omega-icl.github.io/mcpp/): A mature McCormick relaxation package in C++ that also includes McCormick-Taylor, Chebyshev
Polyhedral and Ellipsoidal arithmetics.
