## Solvers with a guarantee of global optimality
For unconstrained problems, the following lower/upper bounding problem modes will
provide a solution that is globally optimal to within an epsilon tolerance:
- Explicit LP
-
-
-

The following solvers can provide a valid solution to constrained problems:
-
-
-
-

The following solvers are currently under-construction and will likely furnish incorrect answers/errors.
- AlphaBB:
- QCQP McCormick:


## Solving a basic problem problem
```julia
println("Test Problem 5 (Mult):")
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
