  workspace()
  using EAGO
  using JuMP
  using MathProgBase

  println("test jumpmodel 9")
  jumpmodel9 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Diff1-MV",
                                         LBDsolvertype = "Ipopt",
                                         probe_depth = -1,
                                         variable_depth = 1000,
                                         DAG_depth = 10,
                                         STD_RR_depth = -1,
                                         validated = true))
  @variable(jumpmodel9, -5 <= x9 <= 5)
  @variable(jumpmodel9, -5 <= y9 <= 5)
  @NLobjective(jumpmodel9, Min, 2*x9^2-1.05*x9^4+(x9^6)/6+x9*y9+y9^2)
  status6 = solve(jumpmodel9)
