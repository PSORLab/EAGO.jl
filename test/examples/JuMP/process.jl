using JuMP, EAGO, Gurobi

m = Model(with_optimizer(EAGO.Optimizer, relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))


# ----- Variables ----- #
x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[5], 0.0)
JuMP.set_lower_bound(x[4], 0.0)
JuMP.set_lower_bound(x[2], 0.0)
JuMP.set_lower_bound(x[3], 0.0)
JuMP.set_lower_bound(x[1], 10.0)
JuMP.set_upper_bound(x[1], 2000.0)
JuMP.set_upper_bound(x[2], 16000.0)
JuMP.set_upper_bound(x[3], 120.0)
JuMP.set_upper_bound(x[4], 5000.0)
JuMP.set_upper_bound(x[5], 2000.0)
JuMP.set_lower_bound(x[6], 85.0)
JuMP.set_upper_bound(x[6], 93.0)
JuMP.set_lower_bound(x[7], 90.0)
JuMP.set_upper_bound(x[7], 95.0)
JuMP.set_lower_bound(x[8], 3.0)
JuMP.set_upper_bound(x[8], 12.0)
JuMP.set_lower_bound(x[9], 1.2)
JuMP.set_upper_bound(x[9], 4.0)
JuMP.set_lower_bound(x[10], 145.0)
JuMP.set_upper_bound(x[10], 162.0)

# ----- Constraints ----- #
@NLconstraint(m, e1, -x[1]*(1.12+0.13167*x[8]-0.00667* (x[8])^2)+x[4] == 0.0)
@constraint(m, e2, -x[1]+1.22*x[4]-x[5] == 0.0)
@NLconstraint(m, e3, -0.001*x[4]*x[9]*x[6]/(98-x[6])+x[3] == 0.0)
@NLconstraint(m, e4, -(1.098*x[8]-0.038* (x[8])^2)-0.325*x[6]+x[7] == 57.425)
@NLconstraint(m, e5, -(x[2]+x[5])/x[1]+x[8] == 0.0)
@constraint(m, e6, x[9]+0.222*x[10] == 35.82)
@constraint(m, e7, -3*x[7]+x[10] == -133.0)

# ----- Objective ----- #
@NLobjective(m, Max, 0.063*x[4]*x[7] - 5.04*x[1] - 0.035*x[2] - 10*x[3] - 3.36*x[5]) # CAN BE QUADRATIC

JuMP.optimize!(m)

run_time = backend(m).optimizer.model.optimizer._run_time
#objective_value =
#objective_bound =
#println("objective between: ")
println("run time: $run_time")
