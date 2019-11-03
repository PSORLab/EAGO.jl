# High-Performance Configuration

## LP Solver Selection

By default, EAGO uses GLPK for solving linear subproblems introduced. Using a
commercial linear solver is highly recommended such as Gurobi, CPLEX, or XPRESS
is highly recommended. Both Gurobi and CPLEX are free for academics and
installation information can be found through http://www.gurobi.com/academia/academia-center and
https://www.ibm.com/developerworks/community/blogs/jfp/entry/CPLEX_Is_Free_For_Students?lang=en, respectively.  

A non-default LP solver can then be selected by the user via a series of keyword argument inputs as illustrated in the code snippet below. The `relaxed_optimizer` contains an instance optimizer with valid relaxations that are made at the root node and is updated with affine relaxations in place. Options can be passed to this optimizer using keyword arguments when initializing EAGO using the with_optimizer syntax in JuMP by
defining an `Iterators.Pairs` structure assigning it to the `relaxed_optimizer_kwargs` keyword argument.
MOI.

```julia

# Create opt EAGO Optimizer with CPLEX for use with MOI routines
opt = EAGO.Optimizer(relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0))

# Create the same model m using an options dictionary in JuMP
relaxed_optimizer_kwargs = Dict{Symbol, Any}()
opt_dict[:relaxed_optimizer] = Gurobi.Optimizer()
opt_dict[:relaxed_optimizer_kwargs] = Iterators.Pairs([:OutputFlag], [0])

m = JuMP.Model(with_optimizer(EAGO.Optimizer; opt_dict...))

# Create the same model m is keyword arguments in JuMP
m = JuMP.Model(with_optimizer(EAGO.Optimizer; relaxed_optimizer = Gurobi.Optimizer(),
                                              relaxed_optimizer_kwargs = Iterators.Pairs([:OutputFlag], [0])))
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
