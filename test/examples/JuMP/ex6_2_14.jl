using Revise
using JuMP, EAGO

m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 1,
                                                    "output_iterations" => 1000,
                                                    "iteration_limit" => 100000,

                                                    "cp_depth" => -1,
                                                    "cp_repetitions" => -1,
                                                    "cp_forward_reverse_limit" => 2,

                                                    "cut_min_iterations" => 2,
                                                    "cut_max_iterations" => 2,
                                                    "objective_cut_on" => true,
                                                    "subgrad_tighten" => true,

                                                    "obbt_depth" => 4,
                                                    "obbt_repetitions" => 4,
                                                    "obbt_aggressive_on" => true,

                                                    "upper_bounding_depth" => 4,

                                                    "fbbt_lp_depth" => 100000,
                                                    "fbbt_lp_repetitions" => 3))

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
