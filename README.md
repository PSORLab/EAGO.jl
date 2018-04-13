# EAGO
EAGO is an open-source development environment for **robust and global optimization** in Julia. 

[![Build Status](https://travis-ci.org/MatthewStuber/EAGO.jl.svg?branch=master)](https://travis-ci.org/MatthewStuber/EAGO.jl)
[![Coverage Status](https://coveralls.io/repos/github/MatthewStuber/EAGO.jl/badge.svg?branch=master)](https://coveralls.io/github/MatthewStuber/EAGO.jl?branch=master)
[![codecov.io](http://codecov.io/github/MatthewStuber/EAGO.jl/coverage.svg?branch=master)](http://codecov.io/github/MatthewStuber/EAGO.jl?branch=master)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://MatthewStuber.github.io/EAGO.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://MatthewStuber.github.io/EAGO.jl/latest)

## Global Optimization

As a global optimization platform, EAGO's solvers can be used to find solutions of general nonconvex problems with a guaranteed certificate of optimality. However, global solvers suffer from the curse of dimensionality and therefore their performance is  outstripped by convex solvers. For users interested in large-scale applications, be warned that problems generally larger than a few variables may prove challenging for certain types of global optimization problems. 

## Package Capabilities

The EAGO package has numerous features: a solver accessible from JuMP/MathProgBase, domain reduction routines, McCormick relaxations, and specialized non-convex semi-infinite program solvers. A full description of all EAGO features is available in the documentation website: https://matthewstuber.github.io/EAGO.jl/.

## News

- 4/12/2018: [EAGO v0.1.0 has been tagged](https://github.com/JuliaLang/METADATA.jl/pull/14218). We're currently moving all the EAGO subpackages (e.g. `EAGODomainReduction.jl`) into the main `EAGO.jl` package. Once this is finalized and the documentations is complete we'll take EAGO v0.2.0. 


## Related Packages

- **ValidatedNumerics.jl**: Provides the interval arithmetic backbone for the SmoothMcCormickGrad package. It allow contains a series of subroutine that allow for locating roots and minimal of unconstrained functions via interval arithmetic. EAGO allows provides an implementation for interval optimization that can be called through JuMP and MathProgBase. In general, McCormick relaxations exhibit superior convergence properties which may help mitigate clustering but will not guarantee correctly rounded bounds.    
- JuMP.jl
