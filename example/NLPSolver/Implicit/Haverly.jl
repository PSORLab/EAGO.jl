using EAGO

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
                                         STD_RR_depth = -1,
                                         ImplicitFlag = true))

set_Implicit_Model!(jm1,f,h,hj,g)

@variable(jm1, -200 <= x1 <= -100)
@variable(jm1, 200 <= x2 <= 400)
@variable(jm1, 200 <= x3 <= 400)
@variable(jm1, 200 <= x4 <= 400)
@variable(jm1, 200 <= x5 <= 400)
@variable(jm1, 200 <= x6 <= 400)
@variable(jm1, 200 <= x7 <= 400)
@variable(jm1, 200 <= x8 <= 400)
@variable(jm1, 200 <= x9 <= 400)

@constraint(jm1, x+2y == 0.0)
@constraint(jm1, x+2y == 0.0)
@constraint(jm1, -x+2y == 0.0)
@constraint(jm1, x+2y == 0.0)
@constraint(jm1, x+2y <= 0)
@constraint(jm1, x+2y <= 0)

@NLobjective(jm1, Min, x*y)

status4 = solve(jm1)
#objval4 = getobjectivevalue(jm1)
#Xval4 = getvalue(x)
#Yval4 = getvalue(y)
