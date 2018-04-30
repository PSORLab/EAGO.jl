#workspace()

using EAGO
using Ipopt
using JuMP

# Haverly Pooling Problem Implicit Model
f(x,p) = 6.0*p[1] + 16.0*p[2] + 10.0*p[3] - 9.0*x[2] + 10.0*p[4] - 15.0*x[4]
h(x,p) = [p[1]+p[2]-9.0*x[1]+x[3];
          0.03*p[1]+0.01*p[2]-(x[1]+x[3])*p[5];
          x[1]+p[3]-x[2];
          x[3]+p[4]-x[4]]
hj(x,p) = [-9.0  0.0  1.0  0.0;
           -p[5] 0.0 -p[5] 0.0;
            1.0 -1.0  0.0  0.0;
            0.0  0.0  1.0 -1.0]
g(x,p) = [x[1]*p[5]+0.02*p[3]-0.025*p[5];
          x[3]*p[5]+0.02*p[4]-0.015*x[4]]

# JuMP Model Setup

jm1 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBDsolvertype = "LP",
                                         probe_depth = -1,
                                         variable_depth = 1000,
                                         DAG_depth = -1,
                                         STD_RR_depth = 1000,
                                         ImplicitFlag = true,
                                         verbosity = "Full",
                                         validated = true))

#jm1 = Model(solver=IpoptSolver())

@variable(jm1, 0.001 <= x1 <= 300)
@variable(jm1, 0.001 <= x2 <= 300)
@variable(jm1, 0.001 <= x3 <= 100)
@variable(jm1, 0.001 <= x4 <= 100)
@variable(jm1, 0.001 <= x5 <= 100)
@variable(jm1, 0.001 <= x6 <= 200)
@variable(jm1, 0.001 <= x7 <= 200)
@variable(jm1, 0.001 <= x8 <= 200)
@variable(jm1, 0.01 <= x9 <= 0.03)

@NLconstraint(jm1, x1 + x2 - 9.0*x3 + x6 == 0.0)
@NLconstraint(jm1, 0.03*x1 + 0.01*x2 - x3*x9 - x6*x9 == 0.0)
@NLconstraint(jm1, x3 + x4 - x5 == 0.0)
@NLconstraint(jm1, x6 + x7 - x8 == 0.0)
@NLconstraint(jm1, x3*x8 + 0.02*x4 - 0.025*x5 <= 0)
@NLconstraint(jm1, x6*x8 + 0.02*x7 - 0.015*x8 <= 0)

@NLobjective(jm1, Min, 6.0*x1 + 16.0*x2 + 10.0*x4 - 9.0*x5 + 10.0*x7 - 15.0*x8)
#solve(jm1)

status = Solve_Implicit(jm1,f,h,hj,g,4)
#status4 = solve(jm1)
#objval4 = getobjectivevalue(jm1)
#Xval4 = getvalue(x)
#Yval4 = getvalue(y)
