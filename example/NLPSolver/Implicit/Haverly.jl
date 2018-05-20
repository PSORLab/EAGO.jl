#workspace()

using EAGO
using Ipopt
using JuMP

# Haverly Pooling Problem Implicit Model
function f(x,p)
          return 6.0*x[1] + 16.0*x[2] + 10.0*p[2] - 9.0*p[1] + 10.0*p[3] - 15.0*p[4]
end
function h(x,p)
          [x[1] + x[2] - x[3] - x[4];
          0.03*x[1] + 0.01*x[2] - (x[3] + x[4])*p[5];
          x[3] + p[2] - p[1];
          x[4] + p[3] - p[4]]
end
hj(x,p) = [1.0  1.0  -1.0  -1.0;
           0.03 0.01 -p[5] -p[5];
            0.0 0.0  1.0  0.0;
            0.0  0.0  0.0 1.0]
g(x,p) = [x[3]*p[5] + 0.02*p[2] - 0.025*p[1];
          x[4]*p[5] + 0.02*p[3] - 0.015*p[4]]

# JuMP Model Setup
implicit = false
if (~implicit)
          println("Explicit Model")
          jm1 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                  LBDsolvertype = "LP",
                                  probe_depth = -1,
                                  variable_depth = 1000,
                                  DAG_depth = -1,
                                  STD_RR_depth = -1,
                                  ImplicitFlag = false,
                                  verbosity = "Full",
                                  validated = true))
else
          println("Implicit Model")
          jm1 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                  LBDsolvertype = "LP",
                                  probe_depth = -1,
                                  variable_depth = -1,
                                  DAG_depth = -1,
                                  STD_RR_depth = -1,
                                  ImplicitFlag = true,
                                  verbosity = "Full",
                                  validated = true))
end

#jm1 = Model(solver=IpoptSolver())

@variable(jm1, 0.001 <= x1 <= 300)
@variable(jm1, 0.001 <= x2 <= 300)
@variable(jm1, 0.001 <= x3 <= 100)
@variable(jm1, 0.001 <= x4 <= 200)
@variable(jm1, 0.001 <= x5 <= 100)
@variable(jm1, 0.001 <= x6 <= 100)
@variable(jm1, 0.001 <= x7 <= 200)
@variable(jm1, 0.001 <= x8 <= 200)
@variable(jm1, 0.01 <= x9 <= 0.03)

@NLconstraint(jm1, x1 + x2 - x3 - x4 == 0.0)
@NLconstraint(jm1, 0.03*x1 + 0.01*x2 - x3*x9 - x4*x9 == 0.0)
@NLconstraint(jm1, x3 + x6 - x5 == 0.0)
@NLconstraint(jm1, x4 + x7 - x8 == 0.0)
@NLconstraint(jm1, x3*x9 + 0.02*x6 - 0.025*x5 <= 0)
@NLconstraint(jm1, x4*x9 + 0.02*x7 - 0.015*x8 <= 0)

@NLobjective(jm1, Min, 6.0*x1 + 16.0*x2 + 10.0*x6 - 9.0*x5 + 10.0*x7 - 15.0*x8)

if (~implicit)
          println("Explicit Solve")
          status4 = solve(jm1)
else
          println("Implicit Solve")
          status = Solve_Implicit(jm1,f,h,hj,g, 4, Imp_gL_Loc = [],
                                         Imp_gU_Loc = [Int64(1); Int64(2)],
                                         Imp_gL = [-Inf; -Inf],
                                         Imp_gU = [Float64(0.0); Float64(0.0)],
                                         Imp_nCons = Int64(2))
end
#objval4 = getobjectivevalue(jm1)
#Xval4 = getvalue(x)
#Yval4 = getvalue(y)
