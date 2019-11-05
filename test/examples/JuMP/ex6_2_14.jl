using JuMP, EAGO, Gurobi, CSV, DataFrames, Ipopt
#m = Model(with_optimizer(Ipopt.Optimizer))

m = Model(with_optimizer(EAGO.Optimizer, relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))

# ----- Variables ----- #
x_Idx = Any[2, 3, 4, 5]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[2], 1.0e-7)
JuMP.set_upper_bound(x[2], 0.5)
JuMP.set_lower_bound(x[3], 1.0e-7)
JuMP.set_upper_bound(x[3], 0.5)
JuMP.set_lower_bound(x[4], 1.0e-7)
JuMP.set_upper_bound(x[4], 0.5)
JuMP.set_lower_bound(x[5], 1.0e-7)
JuMP.set_upper_bound(x[5], 0.5)

# ----- Constraints ----- #
@constraint(m, e2, x[2]+x[3] == 0.5)
@constraint(m, e3, x[4]+x[5] == 0.5)

# ----- Objective ----- #
@NLobjective(m, Min, (x[2]*(log(x[2]/(x[2]+x[4]))+log(x[2]/(x[2]+0.095173*x[4])))+x[4]*(log(x[4]/(x[2]+x[4]))+log(x[4]/(0.30384*x[2]+x[4])))+(x[2]+2.6738*x[4])*log(x[2]+2.6738*x[4])+(0.374*x[2]+x[4])*log(0.374*x[2]+x[4])+2.6738*x[4]*log(x[4]/(x[2]+2.6738*x[4]))+0.374*x[2]*log(x[2]/(0.374*x[2]+x[4]))+x[3]*(log(x[3]/(x[3]+x[5]))+log(x[3]/(x[3]+0.095173*x[5])))+x[5]*(log(x[5]/(x[3]+x[5]))+log(x[5]/(0.30384*x[3]+x[5])))+(x[3]+2.6738*x[5])*log(x[3]+2.6738*x[5])+(0.374*x[3]+x[5])*log(0.374*x[3]+x[5])+2.6738*x[5]*log(x[5]/(x[3]+2.6738*x[5]))+0.374*x[3]*log(x[3]/(0.374*x[3]+x[5]))-3.6838*x[2]*log(x[2])-1.59549*x[4]*log(x[4])-3.6838*x[3]*log(x[3])-1.59549*x[5]*log(x[5])))

JuMP.optimize!(m)

backend_opt = backend(m).optimizer.model.optimizer
last_relaxed_opt = backend_opt.relaxed_optimizer
run_time = backend_opt._run_time
println("run time: $run_time")
