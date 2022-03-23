# High-Performance Configuration

## LP Solver Selection

By default, EAGO uses GLPK for solving linear subproblems introduced. Using a
commercial linear solver is highly recommended such as Gurobi, CPLEX, or XPRESS
is highly recommended. Both Gurobi and CPLEX are free for academics and
installation information can be found through http://www.gurobi.com/academia/academia-center and
https://www.ibm.com/developerworks/community/blogs/jfp/entry/CPLEX_Is_Free_For_Students?lang=en, respectively.

A non-default LP solver can then be selected by the user via a series of keyword argument inputs as illustrated in the code snippet below. The `relaxed_optimizer` contains an instance optimizer with valid relaxations that are made at the root node and is updated with affine relaxations in place.

```julia

# Create opt EAGO Optimizer with Gurobi as a lower subsolver
subsolver_config = SubSolvers(relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0))
eago_factory = () -> EAGO.Optimizer(subsolvers = subsolver_config)
m = Model(eago_factory)
```

## Rounding Mode

The `IntervalArithmetic.jl` package supports a number of different directed rounding
modes. The default directed rounding mode is `:tight`. It is recommended that the
user specify that `:accurate` directed rounding mode be used as it may results
in a significant performance  improvement. Setting a rounding mode can requires
the redefinition of a number of functions. As a result, this should only be done
at the top-level by the user (rather than by using keyword arguments). To set the
rounding mode to `:accurate` using the following syntax when loading the EAGO package
initially:

```julia
using IntervalArithmetic; setrounding(Interval, :accurate)
using EAGO
# REST OF CODE
```


## Ipopt Build

Ipopt is the recommended solver for upper bounding problems. Ipopt's performance is highly
dependent on the linear algebra package used (up to 30x). By default MUMPS is used.
It's recommended that you either compile Ipopt with HSL MA57 or the Pardiso linear
algebra packages with a machine specific Blas library (for Intel users the JuliaPro
MKL version is recommended). For information on this, see the below links:

- Compiling Ipopt: https://www.coin-or.org/Ipopt/documentation/node13.html
- Julia Specifics:
   - Pointing Ipopt to a compiled version:
      - Ipopt Package Info: https://github.com/JuliaOpt/Ipopt.jl
      - Discourse discussion: https://discourse.julialang.org/t/use-ipopt-with-custom-version/9176
   - Issues using Pardiso:
      - Ubuntu: https://github.com/JuliaOpt/Ipopt.jl/issues/106
      - Windows: https://github.com/JuliaOpt/Ipopt.jl/issues/83
- HSL Website: http://www.hsl.rl.ac.uk/ipopt/
- Pardiso Website: https://pardiso-project.org/
