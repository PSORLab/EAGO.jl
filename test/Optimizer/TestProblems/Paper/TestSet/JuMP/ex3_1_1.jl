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
                                     cut_max_iterations = 10,
                                     upper_bounding_interval = 2)

m = Model(opt)

# ----- Variables ----- #
x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[1], 100.0)
JuMP.set_upper_bound(x[1], 10000.0)
JuMP.set_lower_bound(x[2], 1000.0)
JuMP.set_upper_bound(x[2], 10000.0)
JuMP.set_lower_bound(x[3], 1000.0)
JuMP.set_upper_bound(x[3], 10000.0)
JuMP.set_lower_bound(x[4], 10.0)
JuMP.set_upper_bound(x[4], 1000.0)
JuMP.set_lower_bound(x[5], 10.0)
JuMP.set_upper_bound(x[5], 1000.0)
JuMP.set_lower_bound(x[6], 10.0)
JuMP.set_upper_bound(x[6], 1000.0)
JuMP.set_lower_bound(x[7], 10.0)
JuMP.set_upper_bound(x[7], 1000.0)
JuMP.set_lower_bound(x[8], 10.0)
JuMP.set_upper_bound(x[8], 1000.0)


# ----- Constraints ----- #
@constraint(m, e2, 0.0025*x[4]+0.0025*x[6] <= 1.0)
@constraint(m, e3, -0.0025*x[4]+0.0025*x[5]+0.0025*x[7] <= 1.0)
@constraint(m, e4, -0.01*x[5]+0.01*x[8] <= 1.0)
@NLconstraint(m, e5, 100*x[1]-x[1]*x[6]+833.33252*x[4] <= 83333.333)
@NLconstraint(m, e6, x[2]*x[4]-x[2]*x[7]-1250*x[4]+1250*x[5] <= 0.0)
@NLconstraint(m, e7, x[3]*x[5]-x[3]*x[8]-2500*x[5] <= -1.25e6)


# ----- Objective ----- #
@objective(m, Min, x[1]+x[2]+x[3])

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
