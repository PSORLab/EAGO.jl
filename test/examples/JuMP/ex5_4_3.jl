using JuMP, EAGO, Gurobi

m = Model(with_optimizer(EAGO.Optimizer, relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))

# ----- Variables ----- #
x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[5], 0.0)
JuMP.set_lower_bound(x[9], 0.0)
JuMP.set_lower_bound(x[8], 0.0)
JuMP.set_lower_bound(x[7], 0.0)
JuMP.set_lower_bound(x[6], 0.0)
JuMP.set_lower_bound(x[11], 0.0)
JuMP.set_lower_bound(x[10], 0.0)
JuMP.set_lower_bound(x[12], 0.0)
JuMP.set_lower_bound(x[1], 10.0)
JuMP.set_upper_bound(x[1], 350.0)
JuMP.set_lower_bound(x[2], 10.0)
JuMP.set_upper_bound(x[2], 350.0)
JuMP.set_lower_bound(x[3], 10.0)
JuMP.set_upper_bound(x[3], 200.0)
JuMP.set_lower_bound(x[4], 10.0)
JuMP.set_upper_bound(x[4], 200.0)
JuMP.set_upper_bound(x[5], 10.0)
JuMP.set_upper_bound(x[6], 10.0)
JuMP.set_upper_bound(x[7], 10.0)
JuMP.set_upper_bound(x[8], 10.0)
JuMP.set_upper_bound(x[9], 10.0)
JuMP.set_upper_bound(x[10], 10.0)
JuMP.set_upper_bound(x[11], 10.0)
JuMP.set_upper_bound(x[12], 10.0)
JuMP.set_lower_bound(x[13], 150.0)
JuMP.set_upper_bound(x[13], 310.0)
JuMP.set_lower_bound(x[14], 150.0)
JuMP.set_upper_bound(x[14], 310.0)
JuMP.set_lower_bound(x[15], 150.0)
JuMP.set_upper_bound(x[15], 310.0)
JuMP.set_lower_bound(x[16], 150.0)
JuMP.set_upper_bound(x[16], 310.0)


# ----- Constraints ----- #
@constraint(m, e1, x[5]+x[9] == 10.0)
@constraint(m, e2, x[5]-x[6]+x[11] == 0.0)
@constraint(m, e3, x[7]+x[9]-x[10] == 0.0)
@constraint(m, e4, -x[6]+x[7]+x[8] == 0.0)
@constraint(m, e5, -x[10]+x[11]+x[12] == 0.0)
@NLconstraint(m, e6, x[16]*x[11]-x[13]*x[6]+150*x[5] == 0.0)
@NLconstraint(m, e7, x[15]*x[7]-x[14]*x[10]+150*x[9] == 0.0)
@NLconstraint(m, e8, x[6]*x[15]-x[6]*x[13] == 1000.0)
@NLconstraint(m, e9, x[10]*x[16]-x[10]*x[14] == 600.0)
@constraint(m, e10, x[1]+x[15] == 500.0)
@constraint(m, e11, x[2]+x[13] == 250.0)
@constraint(m, e12, x[3]+x[16] == 350.0)
@constraint(m, e13, x[4]+x[14] == 200.0)


# ----- Objective ----- #
@NLobjective(m, Min, (1300*(1000/(0.0333333333333333*x[1]*x[2]+0.166666666666667*x[1]+0.166666666666667*x[2]))^0.6+1300*(600/(0.0333333333333333*x[3]*x[4]+0.166666666666667*x[3]+0.166666666666667*x[4]))^0.6))


JuMP.optimize!(m)
run_time = backend(m).optimizer.model.optimizer._run_time
objval = objective_value(m)
println("run time: $run_time")
