workspace()
using EAGO
using JuMP
using MathProgBase

#=
Solves the optimization problem f(x,y) = (x-5)^2 + (y-3)^2 on the domain
X = [0,10], Y = [0,10] s.t x-y<1, 0<x-y using the MathProgBase interface and
the interval solver.
=#

#=
println("Test Problem 1")
f(x) = (x[1]-5)^2 + (x[2]-3)^2
g(x) = [x[1] - x[2]]
s1 = EAGO_NLPSolver(probe_depth = -1,
                    variable_depth = -1,
                    DAG_depth = -1,
                    STD_RR_depth = -1)
tl = []
m1 = MathProgBase.NonlinearModel(s1)
MathProgBase.loadproblem!(m1, 2, 1, [0.0, 0.0], [10.0, 10.0],
            [0.0], [1.0], :Min, f, g)
MathProgBase.optimize!(m1)
=#


#=
Solves the optimization problem f(x,y) = (x-5)^2 + (y-3)^2 on the domain
X = [0,10], Y = [0,10] s.t x-y<1, 0<x-y using the JuMP and the interval solver.
=#
#=
println("Test Problem 2")
jumpmodel = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Interval",
                                         LBD_problem_relax = "Interval",
                                         LBD_problem_solver = "Interval",
                                         UBD_func_relax = "Interval",
                                         UBD_problem_relax = "Interval",
                                         UBD_problem_solver = "Interval",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel, 0.0 <= x <= 10.0)
@variable(jumpmodel, 0.0 <= y <= 10.0)
@constraint(jumpmodel, 0.0 <= x - y <= 1.0 )
@NLobjective(jumpmodel, Min, (x-5)^2 + (y-3)^2)
status = solve(jumpmodel)
=#

#=
Solves the optimization problem f(x,y) = (x-5)^2 + (y-3)^2 on the domain
X = [0,10], Y = [0,10] s.t x-y<1, 0<x-y using the JuMP and the LP relaxation solver
and Ipopt upper bounding problem with dual-based bound tightening and constraint
propagation enabled.
=#

#=
println("Test Problem 3")
jumpmodel1 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBD_problem_relax = "LP",
                                         LBD_problem_solver = "Clp",
                                         UBD_func_relax = "Original",
                                         UBD_problem_relax = "NLP2",
                                         UBD_problem_solver = "Ipopt",
                                         UBD_full_depth = 2,
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel1, 0.0 <= x <= 10.0)
@variable(jumpmodel1, 0.0 <= y <= 10.0)
@constraint(jumpmodel1, 0.0 <= x - y <= 1.0 )
@NLobjective(jumpmodel1, Min, (x-5)^2 + (y-3)^2)
status1 = solve(jumpmodel1)

println("Test Problem 3a")
jumpmodel1a = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBD_problem_relax = "LP",
                                         LBD_problem_solver = "Clp",
                                         UBD_func_relax = "Original",
                                         UBD_problem_relax = "NLP2",
                                         UBD_problem_solver = "Ipopt",
                                         UBD_full_depth = 2,
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         STD_RR_depth = -1))
@variable(jumpmodel1a, 0.0 <= x <= 10.0)
@variable(jumpmodel1a, 0.0 <= y <= 10.0)
@constraint(jumpmodel1a, 0.0 <= x - y <= 1.0 )
@NLobjective(jumpmodel1a, Min, (x-5)^2 + (y-3)^2)
status1a = solve(jumpmodel1a)
=#

#=
Solves the optimization problem f(x,y) = (x-5)^2 + (y-3)^2 on the domain
X = [0,10], Y = [0,10] s.t x-y<1, 0<x-y using the JuMP and the NLP relaxation solver.
=#

#=
println("Test Problem 4")
jumpmodel2 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBD_problem_relax = "LP",
                                         LBD_problem_solver = "Clp",
                                         UBD_func_relax = "NS-STD-OFF",
                                         UBD_problem_relax = "LP",
                                         UBD_problem_solver = "Clp",
                                         probe_depth = -1,
                                         variable_depth = -1
                                         ))
@variable(jumpmodel2, 0.0 <= x <= 10.0)
@variable(jumpmodel2, 0.0 <= y <= 10.0)
@constraint(jumpmodel2, 0.0 <= x - y <= 1.0 )
@NLobjective(jumpmodel2, Min, (x-5)^2 + (y-3)^2)
status1 = solve(jumpmodel2)
=#
#=
println("Test Problem 5a (Mult):")
jumpmodel4 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Diff2-MV-OFF",
                                         LBD_problem_relax = "NLP2",
                                         LBD_problem_solver = "Ipopt",
                                         UBD_func_relax = "Original",
                                         UBD_problem_relax = "NLP2",
                                         UBD_problem_solver = "Ipopt",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1,
                                         verbosity = "Full"))
@variable(jumpmodel4, 0 <= x <= 400)
@variable(jumpmodel4, 0 <= y <= 200)
@constraint(jumpmodel4, x+2y == 500)
@NLobjective(jumpmodel4, Min, x*y)
status2 = solve(jumpmodel4)
=#
#=
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
=#

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
status4 = solve(jumpmodel4)
objval4 = getobjectivevalue(jumpmodel4)
Xval4 = getvalue(x)
Yval4 = getvalue(y)

jumpmodel6 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
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
