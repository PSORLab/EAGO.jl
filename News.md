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
- 7/5/2019: **EAGO v0.2.1 has been tagged**. This contains fixes for a few minor issues.
  - Bug fix for explicit SIP solving routine that occurred for uncertainty sets of dimension greater than 1.
  - Bug fix for Max objective sense.
