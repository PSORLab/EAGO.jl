using JuMP

m = Model()

# ----- Variables ----- #
@variable(m, objvar)
x_Idx = Any[1, 2]
@variable(m, x[x_Idx])
setlowerbound(x[1], -10.0)
setupperbound(x[1], 10.0)
setlowerbound(x[2], -10.0)
setupperbound(x[2], 10.0)

# ----- Objective ----- #
@NLobjective(m, Min, ((-5+sqrt( (x[1])^2))*(-5+sqrt( (x[1])^2))+(-5+sqrt( (x[2])^2))*(-5+sqrt( (x[2])^2))))

m = m 		 # model get returned when including this script.
