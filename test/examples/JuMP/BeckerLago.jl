using JuMP, EAGO, Gurobi, CSV, DataFrames, Ipopt

m = Model(with_optimizer(
    EAGO.Optimizer,
    relaxed_optimizer = Gurobi.Optimizer(OutputFlag = 0),
))


# ----- Variables ----- #
x_Idx = Any[1, 2]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[1], -10.0)
JuMP.set_upper_bound(x[1], 10.0)
JuMP.set_lower_bound(x[2], -10.0)
JuMP.set_upper_bound(x[2], 10.0)

# ----- Objective ----- #
@NLobjective(
    m,
    Min,
    ((-5 + sqrt((x[1])^2)) * (-5 + sqrt((x[1])^2)) +
     (-5 + sqrt((x[2])^2)) * (-5 + sqrt((x[2])^2)))
)

JuMP.optimize!(m)
val = objective_value(m)
bnd = objective_bound(m)
gap = bnd - val
backend_opt = backend(m).optimizer.model.optimizer
last_relaxed_opt = backend_opt.relaxed_optimizer
run_time = backend_opt._run_time
println("run time: $run_time")
