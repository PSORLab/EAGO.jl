using JuMP, EAGO

opt = with_optimizer(EAGO.Optimizer, cp_depth = 0, cp_reptitions = 0,
                                     obbt_depth = 2,
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
x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[5], 0.0)
JuMP.set_lower_bound(x[4], 0.0)
JuMP.set_lower_bound(x[2], 0.0)
JuMP.set_lower_bound(x[3], 0.0)
JuMP.set_lower_bound(x[1], 10.0)
JuMP.set_upper_bound(x[1], 2000.0)
JuMP.set_upper_bound(x[2], 16000.0)
JuMP.set_upper_bound(x[3], 120.0)
JuMP.set_upper_bound(x[4], 5000.0)
JuMP.set_upper_bound(x[5], 2000.0)
JuMP.set_lower_bound(x[6], 85.0)
JuMP.set_upper_bound(x[6], 93.0)
JuMP.set_lower_bound(x[7], 90.0)
JuMP.set_upper_bound(x[7], 95.0)
JuMP.set_lower_bound(x[8], 3.0)
JuMP.set_upper_bound(x[8], 12.0)
JuMP.set_lower_bound(x[9], 1.2)
JuMP.set_upper_bound(x[9], 4.0)
JuMP.set_lower_bound(x[10], 145.0)
JuMP.set_upper_bound(x[10], 162.0)

# ----- Constraints ----- #
@NLconstraint(m, e1, -x[1]*(1.12+0.13167*x[8]-0.00667* (x[8])^2)+x[4] == 0.0)
@constraint(m, e2, -x[1]+1.22*x[4]-x[5] == 0.0)
@NLconstraint(m, e3, -0.001*x[4]*x[9]*x[6]/(98-x[6])+x[3] == 0.0)
@NLconstraint(m, e4, -(1.098*x[8]-0.038* (x[8])^2)-0.325*x[6]+x[7] == 57.425)
@NLconstraint(m, e5, -(x[2]+x[5])/x[1]+x[8] == 0.0)
@constraint(m, e6, x[9]+0.222*x[10] == 35.82)
@constraint(m, e7, -3*x[7]+x[10] == -133.0)

# ----- Objective ----- #
@NLobjective(m, Max, 0.063*x[4]*x[7] - 5.04*x[1] - 0.035*x[2] - 10*x[3] - 3.36*x[5]) # CAN BE QUADRATIC

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

#lower_time =
#upper_time
