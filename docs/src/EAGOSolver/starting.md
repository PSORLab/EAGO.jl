## Solvers with a guarantee of global optimality
For unconstrained problems, the following upper bounding problem modes will
provide a solution that is globally optimal to within an epsilon tolerance:
- Interval
- LP (Explicit)

The following upper bounding problems solvers can also provide a valid solution to constrained problems:
- SNOPT
- Ipopt
- MPBNonlinear

The following lower bounding solvers are currently under-construction and will likely furnish incorrect answers/errors.
- AlphaBB
- Quadratic
- Ipopt

Lower bounding problem options are:
- Interval
- SNOPT
- LP

The solver contains a hook into JuMP that can be used to solve simple explicit problems
as shown below: 

## Solving a basic problem problem
```julia
jumpmodel4 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBDsolvertype = "LP",
                                         probe_depth = -1,
                                         variable_depth = 1000,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel4, -200 <= x <= -100)
@variable(jumpmodel4, 200 <= y <= 400)
@constraint(jumpmodel4, -500 <= x+2y <= 400)
@NLobjective(jumpmodel4, Min, x*y)
status2 = solve(jumpmodel4)
```
