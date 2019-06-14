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
x_Idx = Any[3, 5, 6]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[3], 1.0e-7)
JuMP.set_upper_bound(x[3], 0.2)
JuMP.set_lower_bound(x[5], 1.0e-7)
JuMP.set_upper_bound(x[5], 0.4)
JuMP.set_lower_bound(x[6], 1.0e-7)
JuMP.set_upper_bound(x[6], 0.4)


# ----- Objective ----- #
@NLobjective(m, Min, ((15.3261663216011*(0.2 - x[3])+23.2043471859416*(0.4 - x[5])+6.69678129464404*x[6])*log(2.1055*(0.2 - x[3])+3.1878*(0.4 - x[5])+0.92*x[6])-2.46348749603266*(0.2 - x[3])-4.33155441248417*(0.4 - x[5])-0.626542690017204*x[6]+6.4661663216011*(0.2 - x[3])*log((0.2 - x[3])/(2.1055*(0.2 - x[3])+3.1878*(0.4 - x[5])+0.92*x[6]))+12.2043471859416*(0.4 - x[5])*log((0.4 - x[5])/(2.1055*(0.2 - x[3])+3.1878*(0.4 - x[5])+0.92*x[6]))+0.696781294644034*x[6]*log(x[6]/(2.1055*(0.2 - x[3])+3.1878*(0.4 - x[5])+0.92*x[6]))+9.86*(0.2 - x[3])*log((0.2 - x[3])/(1.972*(0.2 - x[3])+2.4*(0.4 - x[5])+1.4*x[6]))+12*(0.4 - x[5])*log((0.4 - x[5])/(1.972*(0.2 - x[3])+2.4*(0.4 - x[5])+1.4*x[6]))+7*x[6]*log(x[6]/(1.972*(0.2 - x[3])+2.4*(0.4 - x[5])+1.4*x[6]))+(1.972*(0.2 - x[3])+2.4*(0.4 - x[5])+1.4*x[6])*log(1.972*(0.2 - x[3])+2.4*(0.4 - x[5])+1.4*x[6])+1.972*(0.2 - x[3])*log((0.2 - x[3])/(1.972*(0.2 - x[3])+0.283910843616504*(0.4 - x[5])+3.02002220174195*x[6]))+2.4*(0.4 - x[5])*log((0.4 - x[5])/(1.45991339466884*(0.2 - x[3])+2.4*(0.4 - x[5])+0.415073537580851*x[6]))+1.4*x[6]*log(x[6]/(0.602183324335333*(0.2 - x[3])+0.115623371371275*(0.4 - x[5])+1.4*x[6]))+(15.3261663216011*x[3]+23.2043471859416*x[5]+6.69678129464404*(0.4 - x[6]))*log(2.1055*x[3]+3.1878*x[5]+0.92*(0.4 - x[6]))-2.46348749603266*x[3]-4.33155441248417*x[5]-0.626542690017204*(0.4 - x[6])+6.4661663216011*x[3]*log(x[3]/(2.1055*x[3]+3.1878*x[5]+0.92*(0.4 - x[6])))+12.2043471859416*x[5]*log(x[5]/(2.1055*x[3]+3.1878*x[5]+0.92*(0.4 - x[6])))+0.696781294644034*(0.4 - x[6])*log((0.4 - x[6])/(2.1055*x[3]+3.1878*x[5]+0.92*(0.4 - x[6])))+9.86*x[3]*log(x[3]/(1.972*x[3]+2.4*x[5]+1.4*(0.4 - x[6])))+12*x[5]*log(x[5]/(1.972*x[3]+2.4*x[5]+1.4*(0.4 - x[6])))+7*(0.4 - x[6])*log((0.4 - x[6])/(1.972*x[3]+2.4*x[5]+1.4*(0.4 - x[6])))+(1.972*x[3]+2.4*x[5]+1.4*(0.4 - x[6]))*log(1.972*x[3]+2.4*x[5]+1.4*(0.4 - x[6]))+1.972*x[3]*log(x[3]/(1.972*x[3]+0.283910843616504*x[5]+3.02002220174195*(0.4 - x[6])))+2.4*x[5]*log(x[5]/(1.45991339466884*x[3]+2.4*x[5]+0.415073537580851*(0.4 - x[6])))+1.4*(0.4 - x[6])*log((0.4 - x[6])/(0.602183324335333*x[3]+0.115623371371275*x[5]+1.4*(0.4 - x[6])))-17.2981663216011*(0.2 - x[3])*log((0.2 - x[3]))-25.6043471859416*(0.4 - x[5])*log((0.4 - x[5]))-8.09678129464404*x[6]*log(x[6])-17.2981663216011*x[3]*log(x[3])-25.6043471859416*x[5]*log(x[5])-8.09678129464404*(0.4 - x[6])*log((0.4 - x[6]))))

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
