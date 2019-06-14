using JuMP

m = Model()

# ----- Variables ----- #
x_Idx = Any[1, 2, 3, 4, 5, 6]
@variable(m, x[x_Idx])
setlowerbound(x[5], 0.0)
setlowerbound(x[1], 0.0)
setlowerbound(x[4], 0.0)
setlowerbound(x[2], 0.0)
setlowerbound(x[6], 0.0)
setlowerbound(x[3], 0.0)
setupperbound(x[1], 1.0)
setupperbound(x[2], 1.0)
setupperbound(x[3], 1.0)
setupperbound(x[4], 1.0)
setupperbound(x[5], 1.0)
setupperbound(x[6], 1.0)


# ----- Objective ----- #
@NLobjective(m, Min, -exp(-(10* ((-0.1312)+x[1])^2+0.05* ((-0.1696)+x[2])^2+17* ((-0.5569)+x[3])^2+3.5* ((-0.0124)+x[4])^2+1.7* ((-0.8283)+x[5])^2+8* ((-0.5886)+x[6])^2))+1.2*exp(-(0.05* ((-0.2329)+x[1])^2+10* ((-0.4135)+x[2])^2+17* ((-0.8307)+x[3])^2+0.1* ((-0.3736)+x[4])^2+8* ((-0.1004)+x[5])^2+14* ((-0.9991)+x[6])^2))+3*exp(-(3* ((-0.2348)+x[1])^2+3.5* ((-0.1451)+x[2])^2+1.7* ((-0.3522)+x[3])^2+10* ((-0.2883)+x[4])^2+17* ((-0.3047)+x[5])^2+8* ((-0.665)+x[6])^2))+3.2*exp(-(17* ((-0.4047)+x[1])^2+8* ((-0.8828)+x[2])^2+0.05* ((-0.8732)+x[3])^2+10* ((-0.5743)+x[4])^2+0.1* ((-0.1091)+x[5])^2+14* ((-0.0381)+x[6])^2)))

m = m 		 # model get returned when including this script.
