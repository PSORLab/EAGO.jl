using JuMP, EAGO

opt = with_optimizer(EAGO.Optimizer, cp_depth = 0, cp_reptitions = 0,
                                     obbt_depth = 10,
                                     absolute_tolerance = 1E-3,
                                     relative_tolerance = 1E-3,
                                     subgrad_tighten = false,
                                     obbt_aggressive_on = false,
                                     dbbt_depth = 1000,
                                     reform_epigraph_flag = false,
                                     reform_cse_flag = false,
                                     reform_flatten_flag = false,
                                     poor_man_lp_depth = 0,
                                     poor_man_lp_reptitions = 10,
                                     verbosity = 1,
                                     header_iterations = 20,
                                     output_iterations = 1,
                                     cut_max_iterations = 4,
                                     upper_bounding_interval = 2)

m = Model(opt)

# ----- Variables ----- #
x_Idx = Any[2, 3, 4, 5, 6, 7]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[2], 1.0e-7)
JuMP.set_upper_bound(x[2], 0.2)
JuMP.set_lower_bound(x[3], 1.0e-7)
JuMP.set_upper_bound(x[3], 0.2)
JuMP.set_lower_bound(x[4], 1.0e-7)
JuMP.set_upper_bound(x[4], 0.4)
JuMP.set_lower_bound(x[5], 1.0e-7)
JuMP.set_upper_bound(x[5], 0.4)
JuMP.set_lower_bound(x[6], 1.0e-7)
JuMP.set_upper_bound(x[6], 0.4)
JuMP.set_lower_bound(x[7], 1.0e-7)
JuMP.set_upper_bound(x[7], 0.4)


# ----- Constraints ----- #
@constraint(m, e2, x[2]+x[3] == 0.2)
@constraint(m, e3, x[4]+x[5] == 0.4)
@constraint(m, e4, x[6]+x[7] == 0.4)


# ----- Objective ----- #
@NLobjective(m, Min, ((15.3261663216011*x[2]+23.2043471859416*x[4]+6.69678129464404*x[6])*log(2.1055*x[2]+3.1878*x[4]+0.92*x[6])-2.46348749603266*x[2]-4.33155441248417*x[4]-0.626542690017204*x[6]+6.4661663216011*x[2]*log(x[2]/(2.1055*x[2]+3.1878*x[4]+0.92*x[6]))+12.2043471859416*x[4]*log(x[4]/(2.1055*x[2]+3.1878*x[4]+0.92*x[6]))+0.696781294644034*x[6]*log(x[6]/(2.1055*x[2]+3.1878*x[4]+0.92*x[6]))+9.86*x[2]*log(x[2]/(1.972*x[2]+2.4*x[4]+1.4*x[6]))+12*x[4]*log(x[4]/(1.972*x[2]+2.4*x[4]+1.4*x[6]))+7*x[6]*log(x[6]/(1.972*x[2]+2.4*x[4]+1.4*x[6]))+(1.972*x[2]+2.4*x[4]+1.4*x[6])*log(1.972*x[2]+2.4*x[4]+1.4*x[6])+1.972*x[2]*log(x[2]/(1.972*x[2]+0.283910843616504*x[4]+3.02002220174195*x[6]))+2.4*x[4]*log(x[4]/(1.45991339466884*x[2]+2.4*x[4]+0.415073537580851*x[6]))+1.4*x[6]*log(x[6]/(0.602183324335333*x[2]+0.115623371371275*x[4]+1.4*x[6]))+(15.3261663216011*x[3]+23.2043471859416*x[5]+6.69678129464404*x[7])*log(2.1055*x[3]+3.1878*x[5]+0.92*x[7])-2.46348749603266*x[3]-4.33155441248417*x[5]-0.626542690017204*x[7]+6.4661663216011*x[3]*log(x[3]/(2.1055*x[3]+3.1878*x[5]+0.92*x[7]))+12.2043471859416*x[5]*log(x[5]/(2.1055*x[3]+3.1878*x[5]+0.92*x[7]))+0.696781294644034*x[7]*log(x[7]/(2.1055*x[3]+3.1878*x[5]+0.92*x[7]))+9.86*x[3]*log(x[3]/(1.972*x[3]+2.4*x[5]+1.4*x[7]))+12*x[5]*log(x[5]/(1.972*x[3]+2.4*x[5]+1.4*x[7]))+7*x[7]*log(x[7]/(1.972*x[3]+2.4*x[5]+1.4*x[7]))+(1.972*x[3]+2.4*x[5]+1.4*x[7])*log(1.972*x[3]+2.4*x[5]+1.4*x[7])+1.972*x[3]*log(x[3]/(1.972*x[3]+0.283910843616504*x[5]+3.02002220174195*x[7]))+2.4*x[5]*log(x[5]/(1.45991339466884*x[3]+2.4*x[5]+0.415073537580851*x[7]))+1.4*x[7]*log(x[7]/(0.602183324335333*x[3]+0.115623371371275*x[5]+1.4*x[7]))-17.2981663216011*x[2]*log(x[2])-25.6043471859416*x[4]*log(x[4])-8.09678129464404*x[6]*log(x[6])-17.2981663216011*x[3]*log(x[3])-25.6043471859416*x[5]*log(x[5])-8.09678129464404*x[7]*log(x[7])))

JuMP.optimize!(m)

fval = JuMP.objective_value(m)
TermStatus = JuMP.termination_status(m)
PrimStatus = JuMP.primal_status(m)
psol = JuMP.value.(x)
println("Alg. terminated with a status of $TermStatus and a result code of $PrimStatus")
println("The optimal value is f*=$fval, the solution found is p*=$psol.")

history = backend(m).optimizer.model.optimizer.history
optimizer = backend(m).optimizer.model.optimizer
evaluator_obj = optimizer.working_evaluator_block.evaluator.objective
evaluator_constraints = optimizer.working_evaluator_block.evaluator.constraints
evaluator_constraints_local = optimizer.nlp_data.evaluator.constraints
