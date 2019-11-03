# News for EAGO releases

## v0.1.1
- 4/12/2018: Initial release of combined EAGO packages v0.1.1.

## v0.1.2
- 6/20/2018: [EAGO v0.1.2 has been tagged](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.1.2). Significant speed and functionality updates.

## v0.2.0
- 6/14/2019: [EAGO v0.2.0 has been tagged](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.2.0). This update creates a number of breaking changes to the EAGO API. Please review the use cases provided in the documentation to update examples.
  - Updated to support Julia 1.0+, MathOptInterface (MOI), and MOI construction of subproblems.
  - Additional domain reduction routines available.
  - Support for specialized handling of linear and quadratic terms.
  - Significant performance improvements due to pre-allocation of Wengert tapes and MOI support.
  - A more intuitive API for McCormick relaxation construction.

## v0.2.1
- 7/5/2019: []**EAGO v0.2.1 has been tagged**](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.2.1). This contains fixes for a few minor issues.
  - Bug fix for explicit SIP solving routine that occurred for uncertainty sets of dimension greater than 1.
  - Bug fix for Max objective sense.

## v0.3.0
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
