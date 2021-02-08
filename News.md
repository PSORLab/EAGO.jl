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
- 7/5/2019: [**EAGO v0.2.1 has been tagged**](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.2.1). This contains fixes for a few minor issues.
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

## v0.4.0
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

## v0.5.0
- 11/18/2020: [**EAGO v0.5.0 has been tagged**](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.5.0)
    - Introduces the `register_eago_operators!(m::JuMP.Model)` which can be used
        to register all nonstandard nonlinear terms used in EAGO in any JuMP model.
      - Introduces `positive`, `negative`, `lower_bnd`, `upper_bnd`, and `bnd`
        functions which can be used to enforce bounds on intermediate terms in
        nonlinear expressions (EAGO.Optimizer only).
      - Adds envelopes: `abs2`, `sinpi`, `cospi`, `fma`, `cbrt`.
      - Adds envelopes and functions: `xlogx`
      - Adds envelopes of special functions: `erf`, `erfc`, `erfinv`, `erfcinv`
      - Adds envelopes of activation functions: `relu`, `gelu`, `elu`, `selu`, `swish1`,`sigmoid`, `softsign`, `softplus`,`bisigmoid`, `pentanh`, `leaky_relu`, `param_relu`.
      - Error messages in `sip_explicit` have been made more transparent.
      - Fixes some issues with documentation image rendering and links.
      - Drops appveyor CI and Travis CI in favor of Github Actions.
- 11/18/2020 [**EAGO v0.5.1 has been tagged**](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.5.1)
     - Support for Julia ~1 (with limited functionality for Julia 1.0, 1.1).
- 11/18/2020 **EAGO v0.5.2 has been tagged**
     - Fix user specified branching variables.
