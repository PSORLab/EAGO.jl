using JuMP, EAGO

opt = with_optimizer(EAGO.Optimizer, cp_depth = 0, cp_reptitions = 0,
                                     obbt_depth = 50,
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
JuMP.set_upper_bound(x[2], 0.08)
JuMP.set_lower_bound(x[3], 1.0e-7)
JuMP.set_upper_bound(x[3], 0.08)
JuMP.set_lower_bound(x[4], 1.0e-7)
JuMP.set_upper_bound(x[4], 0.3)
JuMP.set_lower_bound(x[5], 1.0e-7)
JuMP.set_upper_bound(x[5], 0.3)
JuMP.set_lower_bound(x[6], 1.0e-7)
JuMP.set_upper_bound(x[6], 0.62)
JuMP.set_lower_bound(x[7], 1.0e-7)
JuMP.set_upper_bound(x[7], 0.62)


# ----- Constraints ----- #
@constraint(m, e2, x[2]+x[3] == 0.08)
@constraint(m, e3, x[4]+x[5] == 0.3)
@constraint(m, e4, x[6]+x[7] == 0.62)


# ----- Objective ----- #
@NLobjective(m, Min, (x[2]*log(x[2]/(3*x[2]+6*x[4]+x[6]))+x[4]*log(x[4]/(3*x[2]+6*x[4]+x[6]))+x[6]*log(x[6]/(3*x[2]+6*x[4]+x[6]))-0.80323071133189*x[2]+1.79175946922805*x[4]+0.752006*x[6]+(3*x[2]+6*x[4]+1.6*x[6])*log(3*x[2]+6*x[4]+1.6*x[6])+2*x[2]*log(x[2]/(2.00000019368913*x[2]+4.64593*x[4]+0.480353*x[6]))+x[2]*log(x[2]/(1.00772874182154*x[2]+0.724703350369523*x[4]+0.947722362492017*x[6]))+6*x[4]*log(x[4]/(3.36359157977228*x[2]+6*x[4]+1.13841069150863*x[6]))+1.6*x[6]*log(x[6]/(1.6359356134845*x[2]+3.39220996773471*x[4]+1.6*x[6]))+x[3]*log(x[3]/(3*x[3]+6*x[5]+x[7]))+x[5]*log(x[5]/(3*x[3]+6*x[5]+x[7]))+x[7]*log(x[7]/(3*x[3]+6*x[5]+x[7]))-0.80323071133189*x[3]+1.79175946922805*x[5]+0.752006*x[7]+(3*x[3]+6*x[5]+1.6*x[7])*log(3*x[3]+6*x[5]+1.6*x[7])+2*x[3]*log(x[3]/(2.00000019368913*x[3]+4.64593*x[5]+0.480353*x[7]))+x[3]*log(x[3]/(1.00772874182154*x[3]+0.724703350369523*x[5]+0.947722362492017*x[7]))+6*x[5]*log(x[5]/(3.36359157977228*x[3]+6*x[5]+1.13841069150863*x[7]))+1.6*x[7]*log(x[7]/(1.6359356134845*x[3]+3.39220996773471*x[5]+1.6*x[7]))-3*x[2]*log(x[2])-6*x[4]*log(x[4])-1.6*x[6]*log(x[6])-3*x[3]*log(x[3])-6*x[5]*log(x[5])-1.6*x[7]*log(x[7])))

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
