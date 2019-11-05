using JuMP, EAGO, Gurobi, CSV, DataFrames, Ipopt

m = Model(with_optimizer(EAGO.Optimizer, relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))

# ----- Variables ----- #
x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[1], 100.0)
JuMP.set_upper_bound(x[1], 10000.0)
JuMP.set_lower_bound(x[2], 1000.0)
JuMP.set_upper_bound(x[2], 10000.0)
JuMP.set_lower_bound(x[3], 1000.0)
JuMP.set_upper_bound(x[3], 10000.0)
JuMP.set_lower_bound(x[4], 10.0)
JuMP.set_upper_bound(x[4], 1000.0)
JuMP.set_lower_bound(x[5], 10.0)
JuMP.set_upper_bound(x[5], 1000.0)
JuMP.set_lower_bound(x[6], 10.0)
JuMP.set_upper_bound(x[6], 1000.0)
JuMP.set_lower_bound(x[7], 10.0)
JuMP.set_upper_bound(x[7], 1000.0)
JuMP.set_lower_bound(x[8], 10.0)
JuMP.set_upper_bound(x[8], 1000.0)


# ----- Constraints ----- #
@NLconstraint(m, e2, 833.33252*x[4]/x[1]/x[6]+100/x[6]-83333.333/(x[1]*x[6]) <= 1.0)
@NLconstraint(m, e3, 1250*x[5]/x[2]/x[7]+x[4]/x[7]-1250*x[4]/x[2]/x[7] <= 1.0)
@NLconstraint(m, e4, 1250000/(x[3]*x[8])+x[5]/x[8]-2500*x[5]/x[3]/x[8] <= 1.0)
@constraint(m, e5, 0.0025*x[4]+0.0025*x[6] <= 1.0)
@constraint(m, e6, -0.0025*x[4]+0.0025*x[5]+0.0025*x[7] <= 1.0)
@constraint(m, e7, -0.01*x[5]+0.01*x[8] <= 1.0)


# ----- Objective ----- #
@objective(m, Min, x[1] + x[2] + x[3])

JuMP.optimize!(m)

backend_opt = backend(m).optimizer.model.optimizer
last_relaxed_opt = backend_opt.relaxed_optimizer
run_time = backend_opt._run_time
println("run time: $run_time")
