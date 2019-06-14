using JuMP, EAGO

opt = with_optimizer(EAGO.Optimizer, cp_depth = 0, cp_reptitions = 0,
                                     obbt_depth = 0,
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
                                     cut_max_iterations = 1,
                                     upper_bounding_interval = 2)

m = Model(opt)

# ----- Variables ----- #
x_Idx = Any[1, 2, 3]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[1], -2.0)
JuMP.set_upper_bound(x[1], 2.0)
JuMP.set_lower_bound(x[2], -2.0)
JuMP.set_upper_bound(x[2], 2.0)


# ----- Objective ----- #
@NLobjective(m, Min, (1+ (1+x[1]+x[2])^2*(19+3* (x[1])^2-14*x[1]+6*x[1]*x[2]-14*x[2]+3* (x[2])^2))*(30+ (2*x[1]-3*x[2])^2*(18+12* (x[1])^2-32*x[1]-36*x[1]*x[2]+48*x[2]+27* (x[2])^2)))

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
