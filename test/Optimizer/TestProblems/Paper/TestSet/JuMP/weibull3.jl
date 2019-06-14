using JuMP

m = Model()

# ----- Variables ----- #
x_Idx = Any[1, 2, 3, 4, 5]
@variable(m, x[x_Idx])
setlowerbound(x[1], 0.0)
setlowerbound(x[4], 0.0)
setlowerbound(x[2], 0.0)
setlowerbound(x[3], 0.0)
setupperbound(x[1], 100.0)
setupperbound(x[2], 100.0)
setupperbound(x[3], 1.0)
setupperbound(x[4], 10.0)


# ----- Constraints ----- #
@NLconstraint(m, e1, -( (1.23+exp(-0^x[4]*x[3])*x[2]-x[1])^2+ (1.52+exp(-1^x[4]*x[3])*x[2]-x[1])^2+ (2.95+exp(-2^x[4]*x[3])*x[2]-x[1])^2+ (4.34+exp(-3^x[4]*x[3])*x[2]-x[1])^2+ (5.26+exp(-4^x[4]*x[3])*x[2]-x[1])^2+ (5.84+exp(-5^x[4]*x[3])*x[2]-x[1])^2+ (6.21+exp(-6^x[4]*x[3])*x[2]-x[1])^2+ (6.5+exp(-8^x[4]*x[3])*x[2]-x[1])^2+ (6.83+exp(-10^x[4]*x[3])*x[2]-x[1])^2)+x[5] == 0.0)


# ----- Objective ----- #
@objective(m, Min, x[5])

m = m 		 # model get returned when including this script. 
