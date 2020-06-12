using Revise
using JuMP, EAGO

m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 1,
                                                    "output_iterations" => 1000,
                                                    "iteration_limit" => 100000,

                                                    "cp_depth" => -1,
                                                    "cp_repetitions" => -1,
                                                    "cp_forward_reverse_limit" => 2,

                                                    "cut_min_iterations" => 3,
                                                    "cut_max_iterations" => 3,
                                                    "objective_cut_on" => true,
                                                    "subgrad_tighten" => true,

                                                    "obbt_depth" => 4,
                                                    "obbt_repetitions" => 4,
                                                    "obbt_aggressive_on" => true,

                                                    "upper_bounding_depth" => 4,

                                                    "fbbt_lp_depth" => 100000,
                                                    "fbbt_lp_repetitions" => 3))

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
#@NLobjective(m, Min, 3*x[1]^2 + 6*x[1]*x[2] + 3*x[2]^2 - 14*x[1] - 14*x[2])
#@NLobjective(m, Min, 19+3*(x[1])^2-14*x[1]+6*x[1]*x[2]-14*x[2]+3* (x[2])^2)
@NLobjective(m, Min, (1+ (1+x[1]+x[2])^2*(19+3* (x[1])^2-14*x[1]+6*x[1]*x[2]-14*x[2]+3* (x[2])^2))*
(30+ (2*x[1]-3*x[2])^2*(18+12* (x[1])^2-32*x[1]-36*x[1]*x[2]+48*x[2]+27* (x[2])^2)))




JuMP.optimize!(m)
run_time = backend(m).optimizer.model.optimizer._run_time
println("run time: $run_time")
