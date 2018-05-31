workspace()
using EAGO
using JuMP
using MathProgBase

jumpmodel4 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                         LBDsolvertype = "LP",
                                         probe_depth = -1,
                                         variable_depth = 1000,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel4, -200 <= x <= -100)
@variable(jumpmodel4, 200 <= y <= 400)
@constraint(jumpmodel4, -500 <= x+2y <= 400)
@NLobjective(jumpmodel4, Min, x*y)
status4 = solve(jumpmodel4)
objval4 = getobjectivevalue(jumpmodel4)
Xval4 = getvalue(x)
Yval4 = getvalue(y)

jumpmodel6 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                         LBDsolvertype = "LP",
                                         probe_depth = -1,
                                         variable_depth = 1000,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel6, -5 <= x1 <= 5)
@variable(jumpmodel6, -5 <= y1 <= 5)
@NLobjective(jumpmodel6, Min, 2*x1^2-1.05*x1^4+(x1^6)/6+x1*y1+y1^2)
status6 = solve(jumpmodel6)
objval6 = getobjectivevalue(jumpmodel6)
Xval6 = getvalue(x1)
Yval6 = getvalue(y1)
#=
println("Test Problem 6 (Matyas):")
jumpmodel5 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBD_problem_relax = "LP",
                                         LBD_problem_solver = "Clp",
                                         UBD_func_relax = "NS-STD-OFF",
                                         UBD_problem_relax = "LP",
                                         UBD_problem_solver = "Clp",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel5, -10 <= x <= 10)
@variable(jumpmodel5, -10 <= y <= 10)
@NLobjective(jumpmodel5, Min, 0.26*(x^2+y^2)-0.48*x*y)
status3 = solve(jumpmodel5)



println("Test Problem 7 (Three-Hump Camel):")
jumpmodel6 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBD_problem_relax = "LP",
                                         LBD_problem_solver = "Clp",
                                         UBD_func_relax = "NS-STD-OFF",
                                         UBD_problem_relax = "LP",
                                         UBD_problem_solver = "Clp",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel6, -5 <= x <= 5)
@variable(jumpmodel6, -5 <= y <= 5)
@NLobjective(jumpmodel6, Min, 2*x^2-1.05*x^4+(x^6)/6+x*y+y^2)
status4 = solve(jumpmodel6)



println("Test Problem 8 (SQRT):")
jumpmodel7 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBD_problem_relax = "LP",
                                         LBD_problem_solver = "Clp",
                                         UBD_func_relax = "NS-STD-OFF",
                                         UBD_problem_relax = "LP",
                                         UBD_problem_solver = "Clp",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel7, 0.1 <= x <= 2)
@variable(jumpmodel7, 0.1 <= y <= 2)
@NLconstraint(jumpmodel7, (2x)^3-y <= 0)
@NLconstraint(jumpmodel7, (-x+1)^3-y <= 0)
@NLobjective(jumpmodel7, Min, sqrt(y))
status5 = solve(jumpmodel7)



println("Test Problem 8 (McCormick Function): ")

jumpmodel8 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBD_problem_relax = "LP",
                                         LBD_problem_solver = "Clp",
                                         UBD_func_relax = "NS-STD-OFF",
                                         UBD_problem_relax = "LP",
                                         UBD_problem_solver = "Clp",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel8, -1.5 <= x <= 5)
@variable(jumpmodel8, -3 <= y <= 4)
@NLobjective(jumpmodel8, Min, sin(x+y)+(x-y)^2-1.5x+2.5y+1)
status6 = solve(jumpmodel8)



println("Test Problem 9 (Easom Function): ")

jumpmodel9 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBD_problem_relax = "LP",
                                         LBD_problem_solver = "Clp",
                                         UBD_func_relax = "NS-STD-OFF",
                                         UBD_problem_relax = "LP",
                                         UBD_problem_solver = "Clp",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel9, -10 <= x <= 10)
@variable(jumpmodel9, -10 <= y <= 10)
@NLobjective(jumpmodel9, Min, -cos(x)*cos(y)*exp(-(x-3.14)^2-(y-3.14)^2))
status7 = solve(jumpmodel9)



println("Test Problem 10 (N2 Schaffer Function): ")
jumpmodel10 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBD_problem_relax = "LP",
                                         LBD_problem_solver = "Clp",
                                         UBD_func_relax = "NS-STD-OFF",
                                         UBD_problem_relax = "LP",
                                         UBD_problem_solver = "Clp",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel10, -15 <= x <= 10)
@variable(jumpmodel10, -15 <= y <= 10)
@NLobjective(jumpmodel10, Min, 0.5 + (sin(x^2-y^2)^2-0.5)/(1+0.001*(x^2+y^2))^2)
status8 = solve(jumpmodel10)



println("Test Problem 11 (GR & Lee): ")

jumpmodel11 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBD_problem_relax = "LP",
                                         LBD_problem_solver = "Clp",
                                         UBD_func_relax = "NS-STD-OFF",
                                         UBD_problem_relax = "LP",
                                         UBD_problem_solver = "Clp",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel11, 0.5 <= x <= 2.5)
@NLobjective(jumpmodel11, Min, sin(10*3.14*x)/(2x)+(x-1)^4)
status9 = solve(jumpmodel11)



println("Test Problem 12 (Mishra's Bird): ")

jumpmodel12 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBD_problem_relax = "LP",
                                         LBD_problem_solver = "Clp",
                                         UBD_func_relax = "NS-STD-OFF",
                                         UBD_problem_relax = "LP",
                                         UBD_problem_solver = "Clp",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel12, -10 <= x <= 0)
@variable(jumpmodel12, -6.5 <= y <= 0)
@NLconstraint(jumpmodel12, (x+5)^2+(y+5)^2 <= 25)
@NLobjective(jumpmodel12, Min, sin(y)*exp((1-cos(x))^2)+cos(x)*exp((1-sin(y))^2)+(x-y)^2)
status10 = solve(jumpmodel12)
=#

jumpmodel7 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                         LBDsolvertype = "LP",
                                         probe_depth = -1,
                                         variable_depth = 1000,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel7, -5 <= x7 <= 5)
@variable(jumpmodel7, -5 <= y7 <= 5)
@NLobjective(jumpmodel7, Min, 2*x7^2-1.05*x7^4+(x7^6)/6+x7*y7+y7^2)
