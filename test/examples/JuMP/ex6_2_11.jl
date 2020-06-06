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
# ----- Variables ----- #
x_Idx = Any[2, 3, 4]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[2], 1.0e-6)
JuMP.set_upper_bound(x[2], 1.0)
JuMP.set_lower_bound(x[3], 1.0e-6)
JuMP.set_upper_bound(x[3], 1.0)
JuMP.set_lower_bound(x[4], 1.0e-6)
JuMP.set_upper_bound(x[4], 1.0)

# ----- Constraints ----- #
@constraint(m, e2, x[2]+x[3]+x[4] == 1.0)

# ----- Objective ----- #
@NLobjective(m, Min, ((15.3261663216011*x[2]+23.2043471859416*x[3]+6.69678129464404*x[4])*log(2.1055*x[2]+3.1878*x[3]+0.92*x[4])+1.04055250396734*x[2]-2.24199441248417*x[3]+3.1618173099828*x[4]+6.4661663216011*x[2]*log(x[2]/(2.1055*x[2]+3.1878*x[3]+0.92*x[4]))+12.2043471859416*x[3]*log(x[3]/(2.1055*x[2]+3.1878*x[3]+0.92*x[4]))+0.696781294644034*x[4]*log(x[4]/(2.1055*x[2]+3.1878*x[3]+0.92*x[4]))+9.86*x[2]*log(x[2]/(1.972*x[2]+2.4*x[3]+1.4*x[4]))+12*x[3]*log(x[3]/(1.972*x[2]+2.4*x[3]+1.4*x[4]))+7*x[4]*log(x[4]/(1.972*x[2]+2.4*x[3]+1.4*x[4]))+(1.972*x[2]+2.4*x[3]+1.4*x[4])*log(1.972*x[2]+2.4*x[3]+1.4*x[4])+1.972*x[2]*log(x[2]/(1.972*x[2]+0.283910843616504*x[3]+3.02002220174195*x[4]))+2.4*x[3]*log(x[3]/(1.45991339466884*x[2]+2.4*x[3]+0.415073537580851*x[4]))+1.4*x[4]*log(x[4]/(0.602183324335333*x[2]+0.115623371371275*x[3]+1.4*x[4]))-17.2981663216011*x[2]*log(x[2])-25.6043471859416*x[3]*log(x[3])-8.09678129464404*x[4]*log(x[4])))

optimize!(m)
