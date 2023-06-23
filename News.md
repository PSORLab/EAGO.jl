# News for EAGO Releases

## [v0.8.1](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.8.1) (June 15, 2023)

- Resolved an issue where integer and binary variables would sometimes throw a `MathOptInterface.UpperBoundAlreadySet` error.
- Added the function `unbounded_check!` which warns users if they are missing variable bounds and sets them to +/- 1E10 by default.
  - Added an EAGO parameter `unbounded_check` which defaults to `true` and enables `unbounded_check!`.
- Bumped requirement for PrettyTables.jl to v2+ to accommodate the latest version of DataFrames.jl.

## [v0.8.0](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.8.0) (June 12, 2023)

- Updated EAGO for compatibility with the nonlinear expression API changes introduced in JuMP v1.2: https://discourse.julialang.org/t/ann-upcoming-refactoring-of-jumps-nonlinear-api/83052.
  - EAGO now uses the `MOI.Nonlinear` submodule instead of `JuMP._Derivatives`.
  - Models, nodes, expressions, constraints, and operators are now compatible with MOI.
- Added logic and comparison operators to `EAGO.OperatorRegistry`.

## [v0.7.3](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.7.3) (April 11, 2023)

- Bumped DocStringExtensions.jl compatibility.

## [v0.7.2](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.7.2) (November 22, 2022)

- Added support for Julia 1.7.
- Bumped NaNMath.jl compatibility.
- Added `help?` information for various functions and structures.
- Updated documentation and some formatting.

## [v0.7.1](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.7.1) (June 26, 2022)

- Added the function `print_problem_summary`, an internal script used to display all constraints, objectives in a linear program which is added to functions for debug purposes while writing code.
- Adjusted default `EAGOParameters`.
  - `branch_cvx_factor`: 0.5 => 0.25
  - `branch_offset`: 0.2 => 0.15
  - `time_limit` and `_time_left`: 1000.0 => 3600.0
  - `obbt_depth`: 0 => 6
  - `obbt_repetitions`: 1 => 3
  - `cut_tolerance_rel`: 1E-2 => 1E-3
- Adjusted `Ipopt.Optimizer` attributes.
  - `max_iter`: 20000 => 10000
  - `acceptable_iter`: 10000 => 1000
- Excluded `test_quadratic_nonconvex_constraint_basic` from MOI tests.
- Restricted JuMP compatibility to 1.0.0 - 1.1.1.

## [v0.7.0](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.7.0) (March 28, 2022)

- Added envelopes of activation functions: `xabsx`, `logcosh`
- Added `estimator_extrema`, `estimator_under`, and `estimator_over` functions for McCormick relaxations.
- Moved various functions and related structures to new files.
- Added `RelaxCache` structure to hold relaxed problem information.
- Updated forward and reverse propagation.
- Added PrettyTables.jl.
- Added test examples.
- Added a memory allocation analysis.
- Updated documentation.

## [v0.6.1](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.6.1) (March 4, 2021)

- Minor update to tests.

## [v0.6.0](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.6.0) (February 19, 2021)

- License changed from CC BY-NC-SA 4.0 to MIT.
- Fix deprecated Ipopt constructor.
- Fix discrepancy between the returned objective value and the objective evaluated at the solution.
- Dramatically decrease allocates and first-run performance of SIP routines.
- Add two algorithms which modify `SIPRes` detailed in Djelassi, H. and Mitsos A. 2017.
- Fix objective interval fallback function.
- New SIP interface with extendable subroutines.
- Fix x^y relaxation bug.
- Add issues template.
- Add SIP subroutine documentation.

## [v0.5.2](https://github.com/PSORLab/EAGO.jl/commit/bc59c5a8a5e26960c159e06e7b26e2e5c2472956) (November 18, 2020)

- Fix user specified branching variables.

## [v0.5.1](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.5.1) (November 18, 2020)

- Support for Julia ~1 (with limited functionality for Julia 1.0 and 1.1).

## [v0.5.0](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.5.0) (November 18, 2020)

- Introduces the `register_eago_operators!(m::JuMP.Model)` which can be used to register all nonstandard nonlinear terms used in EAGO in any JuMP model.
- Introduces `positive`, `negative`, `lower_bnd`, `upper_bnd`, and `bnd` functions which can be used to enforce bounds on intermediate terms in nonlinear expressions (`EAGO.Optimizer` only).
- Adds envelopes: `abs2`, `sinpi`, `cospi`, `fma`, `cbrt`
- Adds envelopes and functions: `xlogx`
- Adds envelopes of special functions: `erf`, `erfc`, `erfinv`, `erfcinv`
- Adds envelopes of activation functions: `relu`, `gelu`, `elu`, `selu`, `swish`, `sigmoid`, `softsign`, `softplus`, `bisigmoid`, `pentanh`, `leaky_relu`, `param_relu`
- Error messages in `sip_explicit` have been made more transparent.
- Fixes some issues with documentation image rendering and links.
- Drops appveyor CI and Travis CI in favor of GitHub Actions.

## [v0.4.2](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.4.2) (August 28, 2020)

- Support for Julia 1.5.
  
## [v0.4.1](https://github.com/PSORLab/EAGO.jl/commit/9c1bcf024a19840a0ac49c8c6da13619a5f3845f#comments) (June 17, 2020)

- Minor bug fixes.

## [v0.4.0](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.4.0) (June 12, 2020)

- Support for new MOI/JuMP `RawParameter` input and a number of new attributes.
- Separates McCormick and ReverseMcCormick libraries (now [McCormick.jl](https://github.com/PSORLab/McCormick.jl) and [ReverseMcCormick.jl](https://github.com/PSORLab/ReverseMcCormick.jl)) from main package. McCormick.jl is reexported.
- Relaxation calculations now return NaN values on a domain violation.
- Tolerance based validation of cuts has been added to generate numerically safe cuts.
- Significantly simplify internal codebase for `EAGO.Optimizer` (no changes to API): fully decouples input problem specifications from the formulation used internally, stack only stores variables that are branched on, and a number of internal rearrangements to clearly delineate different routines.
- Add problem classification preprocessing that throws to simpler routines if LP problem types are detected (enables future support for SOCP, MILP, MISOCP, and Convex forms).
- Fix multiple bugs and add more transparent error codes.

## [v0.3.1](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.3.1) (January 29, 2020)

- Add unit tests.
- Support for Julia 1.3.
- Fix IntervalContractors.jl dependency issue.

## [v0.3.0](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.3.0) (November 5, 2019)

This update is intended to be the last to create a large number of breaking changes to the EAGO API. Please review the use cases provided in the documentation to update examples.
- A number of performance improvements have been made to the underlying McCormick relaxation library.
- The optimizer used to construct relaxations is now modified in place.
- All subproblem storage has been moved to the `Optimizer` object and storage types (e.g. `LowerInfo`) have been removed.
- A `BinaryMinMaxHeap` structure is now used to store nodes.
- Speed and aesthetics for logging and printing utilities have been updated.
- Subroutines are now customized by creating a subtype of `ExtensionType` and defining subroutines which dispatch on this new structure.
- Parametric interval methods and the Implicit optimizer have been move to a separate package (to be tagged shortly).
- JIT compilation time has been reduced substantially.
- Support for silent tag and time limits.

## [v0.2.1](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.2.1) (July 7, 2019)

- Bug fix for explicit SIP solving routine that occurred for uncertainty sets of dimension greater than 1.
- Bug fix for `MOI.MAX_SENSE` (max objective sense).

## [v0.2.0](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.2.0) (June 14, 2019)

This update creates a number of breaking changes to the EAGO API. Please review the use cases provided in the documentation to update examples.
- Updated to support Julia 1.0+, MathOptInterface (MOI), and MOI construction of subproblems.
- Additional domain reduction routines available.
- Support for specialized handling of linear and quadratic terms.
- Significant performance improvements due to pre-allocation of Wengert tapes and MOI support.
- A more intuitive API for McCormick relaxation construction.

## [v0.1.2](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.1.2) (June 20, 2018)

- Significant speed and functionality updates.

## [v0.1.1](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.1.1) (June 7, 2018)

- Initial release of combined EAGO packages.

## [v0.1.0](https://github.com/PSORLab/EAGO.jl/releases/tag/v0.1.0) (April 10, 2018)

- Main global solver release.
