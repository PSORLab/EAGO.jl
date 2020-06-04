#using JuMP, EAGO, Gurobi, CSV, DataFrames, Ipopt
#m = Model(with_optimizer(Ipopt.Optimizer))

using Revise
using JuMP, EAGO, Ipopt, SCIP
#m = Model(SCIP.Optimizer)
#m = Model(Ipopt.Optimizer)

m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 4,
                                                    "output_iterations" => 1,
                                                    "cp_repetitions" => -1,
                                                    "iteration_limit" => 500,
                                                    "cut_min_iterations" => 0,
                                                    "cut_max_iterations" => 0,
                                                    "objective_cut_on" => false,
                                                    "subgrad_tighten" => false,
                                                    "obbt_depth" => -1))

# cp_repetitions
# obbt_repetitions

# ----- Variables ----- #
x_Idx = Any[1, 2]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[1], -2.0)
JuMP.set_upper_bound(x[1], 2.0)
JuMP.set_lower_bound(x[2], -2.0)
JuMP.set_upper_bound(x[2], 2.0)

# FALSE FATHOM.... sqr, mult, +, -

# ----- Objective ----- #
@NLobjective(m, Min, exp((x[1] + 1.0)^3))
#@NLobjective(m, Min, (1+ (1+x[1]+x[2])^2*(19+3* (x[1])^2-14*x[1]+6*x[1]*x[2]-14*x[2]+3* (x[2])^2))*
#(30+ (2*x[1]-3*x[2])^2*(18+12* (x[1])^2-32*x[1]-36*x[1]*x[2]+48*x[2]+27* (x[2])^2)))




JuMP.optimize!(m)
run_time = backend(m).optimizer.model.optimizer._run_time
println("run time: $run_time")
