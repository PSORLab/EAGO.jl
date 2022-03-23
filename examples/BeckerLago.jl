using JuMP

m = Model()

# ----- Variables ----- #
@variable(m, objvar)
x_Idx = Any[1, 2]
@variable(m, x[x_Idx])
set_lower_bound(x[1], -10.0)
set_upper_bound(x[1], 10.0)
set_lower_bound(x[2], -10.0)
set_upper_bound(x[2], 10.0)


# ----- Constraints ----- #
@NLconstraint(m, e1, -((-5+sqrt( (x[1])^2))*(-5+sqrt( (x[1])^2))+(-5+sqrt( (x[2])^2))*(-5+sqrt( (x[2])^2)))+objvar == 0.0)


# ----- Objective ----- #
@objective(m, Min, objvar)

m = m 		 # model get returned when including this script. 
