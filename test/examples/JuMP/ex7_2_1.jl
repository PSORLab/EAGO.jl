using JuMP, EAGO, Gurobi, CSV, DataFrames, Ipopt

m = Model(with_optimizer(EAGO.Optimizer, relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))

# ----- Variables ----- #
x_Idx = Any[1, 2, 3, 4, 5, 6, 7]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[1], 1500.0)
JuMP.set_upper_bound(x[1], 2000.0)
JuMP.set_lower_bound(x[2], 1.0)
JuMP.set_upper_bound(x[2], 120.0)
JuMP.set_lower_bound(x[3], 3000.0)
JuMP.set_upper_bound(x[3], 3500.0)
JuMP.set_lower_bound(x[4], 85.0)
JuMP.set_upper_bound(x[4], 93.0)
JuMP.set_lower_bound(x[5], 90.0)
JuMP.set_upper_bound(x[5], 95.0)
JuMP.set_lower_bound(x[6], 3.0)
JuMP.set_upper_bound(x[6], 12.0)
JuMP.set_lower_bound(x[7], 145.0)
JuMP.set_upper_bound(x[7], 162.0)


# ----- Constraints ----- #
@NLconstraint(m, e2, 0.0059553571* (x[6])^2+0.88392857*x[3]/x[1]-0.1175625*x[6] <= 1.0)
@NLconstraint(m, e3, 1.1088*x[1]/x[3]+0.1303533*x[1]/x[3]*x[6]-0.0066033*x[1]/x[3]* (x[6])^2 <= 1.0)
@NLconstraint(m, e4, 0.00066173269* (x[6])^2-0.019120592*x[6]-0.0056595559*x[4]+0.017239878*x[5] <= 1.0)
@NLconstraint(m, e5, 56.85075/x[5]+1.08702*x[6]/x[5]+0.32175*x[4]/x[5]-0.03762* (x[6])^2/x[5] <= 1.0)
@NLconstraint(m, e6, 2462.3121*x[2]/x[3]/x[4]-25.125634*x[2]/x[3]+0.006198*x[7] <= 1.0)
@NLconstraint(m, e7, 161.18996/x[7]+5000*x[2]/x[3]/x[7]-489510*x[2]/x[3]/x[4]/x[7] <= 1.0)
@NLconstraint(m, e8, 44.333333/x[5]+0.33*x[7]/x[5] <= 1.0)
@constraint(m, e9, 0.022556*x[5]-0.007595*x[7] <= 1.0)
@constraint(m, e10, -0.0005*x[1]+0.00061*x[3] <= 1.0)
@NLconstraint(m, e11, 0.819672*x[1]/x[3]+0.819672/x[3] <= 1.0)
@NLconstraint(m, e12, 24500*x[2]/x[3]/x[4]-250*x[2]/x[3] <= 1.0)
@NLconstraint(m, e13, 1.2244898e-5*x[3]/x[2]*x[4]+0.010204082*x[4] <= 1.0)
@NLconstraint(m, e14, 6.25e-5*x[1]*x[6]+6.25e-5*x[1]-7.625E-5*x[3] <= 1.0)
@NLconstraint(m, e15, 1.22*x[3]/x[1]+1/x[1]-x[6] <= 1.0)

# ----- Objective ----- #
@NLobjective(m, Min, 3000.0 + 0.035*x[1]*x[6]-0.063*x[3]*x[5]+1.715*x[1]+4.0565*x[3] + 10*x[2])

JuMP.optimize!(m)

backend_opt = backend(m).optimizer.model.optimizer
last_relaxed_opt = backend_opt.relaxed_optimizer
run_time = backend_opt._run_time
println("run time: $run_time")
