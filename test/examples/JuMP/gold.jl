using JuMP, EAGO, Gurobi, CSV, DataFrames, Ipopt
#m = Model(with_optimizer(Ipopt.Optimizer))

m = Model(with_optimizer(EAGO.Optimizer, relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))
# ----- Variables ----- #
x_Idx = Any[1, 2]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[1], -2.0)
JuMP.set_upper_bound(x[1], 2.0)
JuMP.set_lower_bound(x[2], -2.0)
JuMP.set_upper_bound(x[2], 2.0)


# ----- Objective ----- #
@NLobjective(m, Min, (1+ (1+x[1]+x[2])^2*(19+3* (x[1])^2-14*x[1]+6*x[1]*x[2]-14*x[2]+3* (x[2])^2))*(30+ (2*x[1]-3*x[2])^2*(18+12* (x[1])^2-32*x[1]-36*x[1]*x[2]+48*x[2]+27* (x[2])^2)))

JuMP.optimize!(m)
run_time = backend(m).optimizer.model.optimizer._run_time
println("run time: $run_time")
