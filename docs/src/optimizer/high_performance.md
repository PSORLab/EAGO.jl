# High-Performance Configuration

## Linear Programming Solver Selection

By default, EAGO uses GLPK for solving linear subproblems introduced. Using a commercial linear programming (LP) solver is highly recommended, such as Gurobi, CPLEX, or XPRESS. Both Gurobi and CPLEX are free for academics and installation information can be found on the [Gurobi website](http://www.gurobi.com/academia/academia-center) and the [IBM website](https://www.ibm.com/developerworks/community/blogs/jfp/entry/CPLEX_Is_Free_For_Students?lang=en), respectively.

A non-default LP solver can then be selected by the user via a series of keyword argument inputs as illustrated in the code snippet below. The `relaxed_optimizer` contains an instance optimizer with valid relaxations that are made at the root node and is updated with affine relaxations in place.

```julia
# Create opt EAGO Optimizer with Gurobi as a lower subsolver
eago_factory = () -> EAGO.Optimizer(SubSolvers(; r = Gurobi.Optimizer()))
m = Model(eago_factory)
```

## Rounding Mode

The [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) package supports a number of different directed rounding modes. The default directed rounding mode is `:tight`. It is recommended that the user specify that `:accurate` directed rounding mode be used as it may results in a significant performance improvement. Setting a rounding mode can requires the redefinition of a number of functions. As a result, this should only be done at the top-level by the user (rather than by using keyword arguments). To set the rounding mode to `:accurate` using the following syntax when loading the EAGO package initially:

```julia
using IntervalArithmetic; setrounding(Interval, :accurate)
using EAGO
# REST OF CODE
```

## Ipopt Build

Ipopt is the recommended solver for upper-bounding problems. Ipopt's performance is highly dependent on the linear algebra package used (up to 30x). By default MUMPS is used. It is recommended that you either compile Ipopt with HSL MA57 or the Pardiso linear algebra packages with a machine specific Blas library (for Intel users the JuliaPro MKL version is recommended). For information on this, see the links below:

- [Compiling Ipopt](https://www.coin-or.org/Ipopt/documentation/node13.html)
- Julia specifics:
   - Pointing Ipopt to a compiled version:
      - [Ipopt Package Info](https://github.com/JuliaOpt/Ipopt.jl)
      - [Discourse discussion](https://discourse.julialang.org/t/use-ipopt-with-custom-version/9176)
   - Issues using Pardiso:
      - [Ubuntu](https://github.com/JuliaOpt/Ipopt.jl/issues/106)
      - [Windows](https://github.com/JuliaOpt/Ipopt.jl/issues/83)
- [HSL website](http://www.hsl.rl.ac.uk/ipopt/)
- [Pardiso website](https://pardiso-project.org/)
