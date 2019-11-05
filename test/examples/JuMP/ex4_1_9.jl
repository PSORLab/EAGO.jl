using JuMP, EAGO, Gurobi

m = Model(with_optimizer(EAGO.Optimizer, relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))

# ----- Variables ----- #
x_Idx = Any[1, 2]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[1], 0.0)
JuMP.set_lower_bound(x[2], 0.0)
JuMP.set_upper_bound(x[1], 3.0)
JuMP.set_upper_bound(x[2], 4.0)

# ----- Constraints ----- #
@NLconstraint(m, e2, 8* (x[1])^3-2* (x[1])^4-8* (x[1])^2+x[2] <= 2.0)
@NLconstraint(m, e3, 32* (x[1])^3-4* (x[1])^4-88* (x[1])^2+96*x[1]+x[2] <= 36.0)

# ----- Objective ----- #
@objective(m, Min, -x[1] -x[2])

JuMP.optimize!(m)
run_time = backend(m).optimizer.model.optimizer._run_time
objval = objective_value(m)
println("run time: $run_time")
