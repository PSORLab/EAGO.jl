module Test_Explicit_NLP

using Compat
using Compat.Test
using EAGO
using JuMP
using MathProgBase

@testset "JuMP Implicit (LP)" begin
###############################################################################
# TEST NONVALIDATED CALCULATIONS WITHOUT INEQUALITY CONSTRAINTS
###############################################################################
# Solves the quadratically constrained problem (Example 5.1, Stuber2015)
h1(x,p) = [x[1]^2 + p[1]*x[1] + 4.0]
hj1(x,p) = [2.0*x[1] + p[1]]
f1(x,p) = x[1]
LBD1a_func(i) = (i==1) ? (-0.78) : (6.0)
UBD1a_func(i) = (i==1) ? (-0.4) : (9.0)
LBD1b_func(i) = (i==1) ? (-10.0) : (6.0)
UBD1b_func(i) = (i==1) ? (-5.0) : (9.0)

jm1a = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
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


jm1b = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
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


###############################################################################
# TEST VALIDATED CALCULATIONS WITHOUT INEQUALITY CONSTRAINTS
###############################################################################
jm1c = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
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


jm1d = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
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

###############################################################################
# TEST VALIDATED CALCULATIONS WITH INEQUALITY CONSTRAINTS
###############################################################################
jm1e = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                   LBDsolvertype = "LP",
                                   probe_depth = -1,
                                   variable_depth = 1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Normal",
                                   validated = true))
xe = @variable(jm1e, [i=1:2], lowerbound=LBD1b_func(i), upperbound=UBD1b_func(i))
@NLconstraint(jm1e, xe[1]^2 + xe[2]*xe[1] + 4 == 0.0 )
@NLconstraint(jm1e, -100 <= -xe[1] <= 8.0 )
@NLobjective(jm1e, Min, xe[1])
status1e = Solve_Implicit(jm1e,f1,h1,hj1,(x,p)->[-x[1]],1,Imp_gL_Loc = [Int64(1)],
                             Imp_gU_Loc = [Int64(1)],
                             Imp_gL = [Float64(-100.0)],
                             Imp_gU = [Float64(8.0)],
                             Imp_nCons = Int64(1))
sol1e = getsolution(internalmodel(jm1e))
@test status1e == :Optimal
@test isapprox(sol1e[1],-8.0,atol=1E-3)
@test isapprox(sol1e[2],8.50,atol=1E-3)
@test isapprox(getobjval(internalmodel(jm1e)),-8.0,atol=1E-3) # getobjval returns void???

###############################################################################
# TEST NONVALIDATED CALCULATIONS WITH INEQUALITY CONSTRAINTS
###############################################################################

jm1f = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                   LBDsolvertype = "LP",
                                   probe_depth = -1,
                                   variable_depth = 1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Normal",
                                   validated = false))
xf = @variable(jm1f, [i=1:2], lowerbound=LBD1b_func(i), upperbound=UBD1b_func(i))
@NLconstraint(jm1f, xf[1]^2 + xf[2]*xf[1] + 4 == 0.0 )
@NLconstraint(jm1f, -100 <= -xf[1] <= 8.0 )
@NLobjective(jm1f, Min, xf[1])
status1f = Solve_Implicit(jm1f,f1,h1,hj1,(x,p)->[-x[1]],1,Imp_gL_Loc = [Int64(1)],
                             Imp_gU_Loc = [Int64(1)],
                             Imp_gL = [Float64(-100.0)],
                             Imp_gU = [Float64(8.0)],
                             Imp_nCons = Int64(1))
sol1f = getsolution(internalmodel(jm1f))
@test status1f == :Optimal
@test isapprox(sol1f[1],-8.0,atol=1E-3)
@test isapprox(sol1f[2],8.50,atol=1E-3)
@test isapprox(getobjval(internalmodel(jm1f)),-8.0,atol=1E-3) # getobjval returns void???

###############################################################################
# TEST NEWTON METHOD (1D) WITH INEQUALITY CONSTRANTS
###############################################################################
MCOPT = mc_opts()
set_default!(MCOPT)
MCOPT.CTyp = :Newton

jm1g = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                   LBDsolvertype = "LP",
                                   probe_depth = -1,
                                   variable_depth = 1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Full",
                                   validated = true,
                                   PSmcOpt = MCOPT))
xg = @variable(jm1g, [i=1:2], lowerbound=LBD1b_func(i), upperbound=UBD1b_func(i))
@NLconstraint(jm1g, xg[1]^2 + xg[2]*xg[1] + 4 == 0.0 )
@NLconstraint(jm1g, -100 <= -xg[1] <= 8.0 )
@NLobjective(jm1g, Min, xg[1])
status1g = Solve_Implicit(jm1g,f1,h1,hj1,(x,p)->[-x[1]],1,Imp_gL_Loc = [Int64(1)],
                             Imp_gU_Loc = [Int64(1)],
                             Imp_gL = [Float64(-100.0)],
                             Imp_gU = [Float64(8.0)],
                             Imp_nCons = Int64(1))
sol1g = getsolution(internalmodel(jm1g))
@test status1g == :Optimal
@test isapprox(sol1g[1],-8.0,atol=1E-3)
@test isapprox(sol1g[2],8.50,atol=1E-3)
@test isapprox(getobjval(internalmodel(jm1g)),-8.0,atol=1E-3) # getobjval returns void???

###############################################################################
# TEST NEWTON METHOD (1D) WITHOUT CONSTRANTS (Using Midpoint Upperbound)
###############################################################################
jm1z = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                   LBDsolvertype = "LP",
                                   UBDsolvertype = "Interval",
                                   probe_depth = -1,
                                   variable_depth = 1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Normal",
                                   validated = true))
xz = @variable(jm1z, [i=1:2], lowerbound=LBD1b_func(i), upperbound=UBD1b_func(i))
@NLconstraint(jm1z, xz[1]^2 + xz[2]*xz[1] + 4 == 0.0 )
@NLobjective(jm1z, Min, xz[1])
status1dn = Solve_Implicit(jm1z,f1,h1,hj1,x->[],1)

sol1dn = getsolution(internalmodel(jm1z))
@test status1dn == :Optimal
@test isapprox(sol1dn[1],-8.531128966881054,atol=1E-3)
@test isapprox(sol1dn[2],9.00,atol=1E-3)
@test isapprox(getobjval(internalmodel(jm1d)),-8.531128966881054,atol=1E-3)

end
end
