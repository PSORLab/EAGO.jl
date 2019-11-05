using JuMP, EAGO, Gurobi

m = Model(with_optimizer(EAGO.Optimizer, relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))

# ----- Variables ----- #
x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[21], 0.0)
JuMP.set_lower_bound(x[1], -0.5)
JuMP.set_upper_bound(x[1], 0.5)
JuMP.set_lower_bound(x[2], 5.4)
JuMP.set_upper_bound(x[2], 6.4)
JuMP.set_lower_bound(x[3], 0.4)
JuMP.set_upper_bound(x[3], 1.4)
JuMP.set_lower_bound(x[4], 4.9)
JuMP.set_upper_bound(x[4], 5.9)
JuMP.set_lower_bound(x[5], 1.3)
JuMP.set_upper_bound(x[5], 2.3)
JuMP.set_lower_bound(x[6], 3.9)
JuMP.set_upper_bound(x[6], 4.9)
JuMP.set_lower_bound(x[7], 2.1)
JuMP.set_upper_bound(x[7], 3.1)
JuMP.set_lower_bound(x[8], 4.1)
JuMP.set_upper_bound(x[8], 5.1)
JuMP.set_lower_bound(x[9], 2.8)
JuMP.set_upper_bound(x[9], 3.8)
JuMP.set_lower_bound(x[10], 3.0)
JuMP.set_upper_bound(x[10], 4.0)
JuMP.set_lower_bound(x[11], 3.9)
JuMP.set_upper_bound(x[11], 4.9)
JuMP.set_lower_bound(x[12], 3.2)
JuMP.set_upper_bound(x[12], 4.2)
JuMP.set_lower_bound(x[13], 4.7)
JuMP.set_upper_bound(x[13], 5.7)
JuMP.set_lower_bound(x[14], 2.3)
JuMP.set_upper_bound(x[14], 3.3)
JuMP.set_lower_bound(x[15], 5.6)
JuMP.set_upper_bound(x[15], 6.6)
JuMP.set_lower_bound(x[16], 2.3)
JuMP.set_upper_bound(x[16], 3.3)
JuMP.set_lower_bound(x[17], 6.0)
JuMP.set_upper_bound(x[17], 7.0)
JuMP.set_lower_bound(x[18], 1.9)
JuMP.set_upper_bound(x[18], 2.9)
JuMP.set_lower_bound(x[19], 6.9)
JuMP.set_upper_bound(x[19], 7.9)
JuMP.set_lower_bound(x[20], 1.0)
JuMP.set_upper_bound(x[20], 2.0)
JuMP.set_upper_bound(x[21], 10.0)
JuMP.set_lower_bound(x[22], -2.0)
JuMP.set_upper_bound(x[22], 2.0)


# ----- Constraints ----- #
@NLconstraint(m, e2, x[22]*x[1]-x[2]+x[21] == 0.0)
@NLconstraint(m, e3, x[22]*x[3]-x[4]+x[21] == 0.0)
@NLconstraint(m, e4, x[22]*x[5]-x[6]+x[21] == 0.0)
@NLconstraint(m, e5, x[22]*x[7]-x[8]+x[21] == 0.0)
@NLconstraint(m, e6, x[22]*x[9]-x[10]+x[21] == 0.0)
@NLconstraint(m, e7, x[22]*x[11]-x[12]+x[21] == 0.0)
@NLconstraint(m, e8, x[22]*x[13]-x[14]+x[21] == 0.0)
@NLconstraint(m, e9, x[22]*x[15]-x[16]+x[21] == 0.0)
@NLconstraint(m, e10, x[22]*x[17]-x[18]+x[21] == 0.0)
@NLconstraint(m, e11, x[22]*x[19]-x[20]+x[21] == 0.0)


# ----- Objective ----- #
@NLobjective(m, Min, ( (x[1])^2+ (x[2]-5.9)^2+ (x[3]-0.9)^2+ (x[4]-5.4)^2+ (x[5]-1.8)^2+ (x[6]-4.4)^2+ (x[7]-2.6)^2+ (x[8]-4.6)^2+ (x[9]-3.3)^2+ (x[10]-3.5)^2+ (x[11]-4.4)^2+ (x[12]-3.7)^2+ (x[13]-5.2)^2+ (x[14]-2.8)^2+ (x[15]-6.1)^2+ (x[16]-2.8)^2+ (x[17]-6.5)^2+ (x[18]-2.4)^2+ (x[19]-7.4)^2+ (x[20]-1.5)^2))

JuMP.optimize!(m)

backend_opt = backend(m).optimizer.model.optimizer
last_relaxed_opt = backend_opt.relaxed_optimizer
run_time = backend_opt._run_time
println("run time: $run_time")
