using JuMP

m = Model()

# ----- Variables ----- #
x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
@variable(m, x[x_Idx])
setlowerbound(x[11], 0.0)
setupperbound(x[11], 10.0)


# ----- Constraints ----- #
@NLconstraint(m, e2, x[10]* (x[11])^4-x[8]* (x[11])^2+x[6] == 0.0)
@NLconstraint(m, e3, x[9]* (x[11])^2-x[7] == 0.0)
@constraint(m, e4, -x[1]-x[12] <= -10.0)
@constraint(m, e5, x[1]-x[12] <= 10.0)
@constraint(m, e6, x[2]-0.1*x[12] <= 1.0)
@constraint(m, e7, -x[2]-0.1*x[12] <= -1.0)
@constraint(m, e8, -x[3]-0.1*x[12] <= -1.0)
@constraint(m, e9, x[3]-0.1*x[12] <= 1.0)
@constraint(m, e10, -x[4]-0.01*x[12] <= -0.2)
@constraint(m, e11, x[4]-0.01*x[12] <= 0.2)
@constraint(m, e12, -x[5]-0.005*x[12] <= -0.05)
@constraint(m, e13, x[5]-0.005*x[12] <= 0.05)
@NLconstraint(m, e14, -54.387*x[3]*x[2]+x[6] == 0.0)
@NLconstraint(m, e15, -0.2*(1364.67*x[3]*x[2]-147.15*x[4]*x[3]*x[2])+5.544*x[5]+x[7] == 0.0)
@NLconstraint(m, e16, -3*(-9.81*x[3]* (x[2])^2-9.81*x[3]*x[1]*x[2]-4.312* (x[3])^2*x[2]+264.896*x[3]*x[2]+x[4]*x[5]-9.274*x[5])+x[8] == 0.0)
@NLconstraint(m, e17, -(7*x[4]* (x[3])^2*x[2]-64.918* (x[3])^2*x[2]+380.067*x[3]*x[2]+3*x[5]*x[2]+3*x[5]*x[1])+x[9] == 0.0)
@NLconstraint(m, e18, - (x[3])^2*x[2]*(7*x[1]+4*x[2])+x[10] == 0.0)


# ----- Objective ----- #
@objective(m, Min, x[12])

m = m 		 # model get returned when including this script.
