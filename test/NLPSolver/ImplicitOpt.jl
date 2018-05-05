module Test_Explicit_NLP

using Compat
using Compat.Test
using EAGO
using JuMP
using MathProgBase

@testset "JuMP Implicit (LP)" begin

# Solves the quadratically constrained problem (Example 5.1, Stuber2015)
h1(x,p) = [x[1]^2 + p[1]*x[1] + 4.0]
hj1(x,p) = [2.0*x[1] + p[1]]
f1(x,p) = x[1]
LBD1a_func(i) = (i==1) ? (-0.78) : (6.0)
UBD1a_func(i) = (i==1) ? (-0.4) : (9.0)
LBD1b_func(i) = (i==1) ? (-10.0) : (6.0)
UBD1b_func(i) = (i==1) ? (-5.0) : (9.0)

jm1a = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                   LBDsolvertype = "LP",
                                   probe_depth = -1,
                                   variable_depth = 1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Normal",
                                   validated = false))
xa = @variable(jm1a, [i=1:2], lowerbound=LBD1a_func(i), upperbound=UBD1a_func(i))
@NLconstraint(jm1a, xa[1]^2 + xa[2]*xa[1] + 4 == 0.0 )
@NLobjective(jm1a, Min, xa[1])
status1a = Solve_Implicit(jm1a,f1,h1,hj1,x->[],1)


jm1b = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                   LBDsolvertype = "LP",
                                   probe_depth = -1,
                                   variable_depth = 1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Normal",
                                   validated = false))
xb = @variable(jm1b, [i=1:2], lowerbound=LBD1b_func(i), upperbound=UBD1b_func(i))
@NLconstraint(jm1b, xb[1]^2 + xb[2]*xb[1] + 4 == 0.0 )
@NLobjective(jm1b, Min, xb[1])
status1b = Solve_Implicit(jm1b,f1,h1,hj1,x->[],1)


sol1a = getsolution(internalmodel(jm1a))
@test status1a == :Optimal
@test isapprox(sol1a[1],-0.7639320302435194,atol=1E-3)
@test isapprox(sol1a[2],6.00,atol=1E-3)
@test isapprox(getobjval(internalmodel(jm1a)),-0.7639320302435194,atol=1E-3)

sol1b = getsolution(internalmodel(jm1b))
@test status1b == :Optimal
@test isapprox(sol1b[1],-8.531128966881054,atol=1E-3)
@test isapprox(sol1b[2],9.00,atol=1E-3)
@test isapprox(getobjval(internalmodel(jm1b)),-8.531128966881054,atol=1E-3)

jm1c = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                   LBDsolvertype = "LP",
                                   probe_depth = -1,
                                   variable_depth = 1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Normal",
                                   validated = true))
xc = @variable(jm1c, [i=1:2], lowerbound=LBD1a_func(i), upperbound=UBD1a_func(i))
@NLconstraint(jm1c, xc[1]^2 + xc[2]*xc[1] + 4 == 0.0 )
@NLobjective(jm1c, Min, xc[1])
status1c = Solve_Implicit(jm1c,f1,h1,hj1,x->[],1)


jm1d = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                   LBDsolvertype = "LP",
                                   probe_depth = -1,
                                   variable_depth = 1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Normal",
                                   validated = true))
xd = @variable(jm1d, [i=1:2], lowerbound=LBD1b_func(i), upperbound=UBD1b_func(i))
@NLconstraint(jm1d, xd[1]^2 + xd[2]*xd[1] + 4 == 0.0 )
@NLobjective(jm1d, Min, xd[1])
status1d = Solve_Implicit(jm1d,f1,h1,hj1,x->[],1)


sol1c = getsolution(internalmodel(jm1a))
@test status1c == :Optimal
@test isapprox(sol1c[1],-0.7639320302435194,atol=1E-3)
@test isapprox(sol1c[2],6.00,atol=1E-3)
@test isapprox(getobjval(internalmodel(jm1c)),-0.7639320302435194,atol=1E-3)

sol1d = getsolution(internalmodel(jm1d))
@test status1d == :Optimal
@test isapprox(sol1d[1],-8.531128966881054,atol=1E-3)
@test isapprox(sol1d[2],9.00,atol=1E-3)
@test isapprox(getobjval(internalmodel(jm1d)),-8.531128966881054,atol=1E-3)

end
end
