# High-Performance Configuration

### Solver Parameters

EAGO defaults to using reasonable estimates of performant solver parameters. Parameter considerations for explicit optimization:
- **Validation**: The validated interval arithmetic option comes with a
                  significant performance decrease but can be useful for some
                  problems.
- **OBBT**: Selecting an arbitrary high depth for range reduction for constrained problems may significantly improve performance on some problems and dramatically reduce solution times of others.
- **DBBT**: Recommend selecting an arbitrary high depth for duality-based tightening.
- **Constraint Propagation**: Recommended using for problems with highly nonlinear and complex constraints.

Parameter considerations specifically for implicit optimization:
- **Interval Contractor**: Run roughly 5-10 interval iterations per McCormick contractor iteration. Recommend starting with 10.
- **McCormick Contractor**: Limit iterations to three or fewer.

### Lower Bounding Problem

Currently, two (McCormick-based) lower bounding problems are available: one generates
an affine relaxation of the McCormick relaxation while the other provides the relaxation
directly to the solver for use. In general, the affine relaxation can be solved using
an LP and may solve quickly as a result. If the McCormick relaxation is used, then
either a nonsmooth local nlp solver should be used or local

### LP Solver Selection

By default, EAGO uses GLPK for solving linear subproblems introduced. Using a
commercial linear solver is highly recommended such as Gurobi, CPLEX, or XPRESS
is highly recommended. Both Gurobi and CPLEX are free for academics and
installation information can be found through http://www.gurobi.com/academia/academia-center and
https://www.ibm.com/developerworks/community/blogs/jfp/entry/CPLEX_Is_Free_For_Students?lang=en, respectively.  

A non-default LP solver can then be selected by the user via a series of keyword argument inputs as illustrated in the code snippet below. The `initial_relaxed_optimizer` contains an optimizer with valid relaxations that are made at the root node. This object can be copied to the `working_relaxed_optimizer` which contains all valid relaxations of the current node and will solve the lower subproblem. If the Optimizer does not support resolving then a new structure is built using the `lower_factory` each time. The `lower_optimizer_options` contain keyword arguments for options that are passed to the linear optimizer when constructed using the factory option.

```julia

# Create opt EAGO Optimizer with CPLEX
opt = EAGO.Optimizer(initial_relaxed_optimizer = CPLEX.Optimizer(), working_relaxed_optimizer = CPLEX.Optimizer(), lower_factory = JuMP.with_optimizer(CPLEX.Optimizer), linear_optimizer = CPLEX.Optimizer(), lower_optimizer_options = Dict{Symbol,Any}())

# Create the same model m using an options dictionary in JuMP
opt_dict = Dict{Symbol, Any}()
opt_dict[:initial_relaxed_optimizer] = CPLEX.Optimizer()
opt_dict[:working_relaxed_optimizer] = CPLEX.Optimizer()
opt_dict[:lower_factory] = JuMP.with_optimizer(CPLEX.Optimizer)
opt_dict[:linear_optimizer] = CPLEX.Optimizer()
opt_dict[:lower_optimizer_options] = Dict{Symbol,Any}()

m = JuMP.Model(with_optimizer(EAGO.Optimizer; opt_dict...))
```

### Ipopt Build

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

### Implicit Relaxation Notes
When the implicit routines are used, it's recommended that the user make use of the
provided midpoint Newton Gauss-Siedel bound provided. Selectively branching only on
the $p$ variables but solving the problem in the full $y$ space can prevent the
upper bound from ever converging to the solution for some nlp solvers.
