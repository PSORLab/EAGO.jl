using JuMP, EAGO

opt = with_optimizer(EAGO.Optimizer,
                     cp_depth = 5,
                     cp_interval_reptitions = 2,
                     obbt_depth = 0,
                     absolute_tolerance = 1E-3,
                     relative_tolerance = 1E-3,
                     subgrad_tighten = true,
                     obbt_aggressive_on = false,
                     dbbt_depth = 1000,
                     reform_epigraph_flag = false,
                     reform_cse_flag = false,
                     reform_flatten_flag = false,
                     poor_man_lp_depth = 0,
                     poor_man_lp_reptitions = 10,
                     verbosity = 4,
                     header_iterations = 10,
                     output_iterations = 1,
                     cut_max_iterations = 3,
                     upper_bounding_interval = 2,
                     iteration_limit = 2)

m = Model(opt)

# ----- Variables ----- #
x_Idx = Any[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[5], 0.0)
JuMP.set_lower_bound(x[4], 0.0)
JuMP.set_lower_bound(x[2], 0.0)
JuMP.set_lower_bound(x[6], 0.0)
JuMP.set_lower_bound(x[3], 0.0)
JuMP.set_upper_bound(x[2], 2.0)
JuMP.set_upper_bound(x[3], 1.6)
JuMP.set_upper_bound(x[4], 1.2)
JuMP.set_upper_bound(x[5], 5.0)
JuMP.set_upper_bound(x[6], 2.0)
JuMP.set_lower_bound(x[7], 0.85)
JuMP.set_upper_bound(x[7], 0.93)
JuMP.set_lower_bound(x[8], 0.9)
JuMP.set_upper_bound(x[8], 0.95)
JuMP.set_lower_bound(x[9], 3.0)
JuMP.set_upper_bound(x[9], 12.0)
JuMP.set_lower_bound(x[10], 1.2)
JuMP.set_upper_bound(x[10], 4.0)
JuMP.set_lower_bound(x[11], 1.45)
JuMP.set_upper_bound(x[11], 1.62)
JuMP.set_lower_bound(x[12], 0.99)
JuMP.set_upper_bound(x[12], 1.01010101010101)
JuMP.set_lower_bound(x[13], 0.99)
JuMP.set_upper_bound(x[13], 1.01010101010101)
JuMP.set_lower_bound(x[14], 0.9)
JuMP.set_upper_bound(x[14], 1.11111111111111)
JuMP.set_lower_bound(x[15], 0.99)
JuMP.set_upper_bound(x[15], 1.01010101010101)


# ----- Constraints ----- #
@constraint(m, e2, -0.819672131147541*x[2]+x[5]-0.819672131147541*x[6] == 0.0)
#@NLconstraint(m, e3, 0.98*x[4]-x[7]*(0.01*x[5]*x[10]+x[4]) == 0.0)
@NLconstraint(m, e4, -x[2]*x[9]+10*x[3]+x[6] == 1.0)
#@NLconstraint(m, e5, x[5]*x[12]-x[2]*(1.12+0.13167*x[9]-0.0067*x[9]*x[9]) == 0.0)
#@NLconstraint(m, e6, x[8]*x[13]-0.01*(1.098*x[9]-0.038*x[9]*x[9])-0.325*x[7] == 0.57425)
#@NLconstraint(m, e7, x[10]*x[14]+22.2*x[11] == 35.82)
#@NLconstraint(m, e8, x[11]*x[15]-3*x[8] == -1.33)


# ----- Objective ----- #
#@NLobjective(m, Min, 5.04*x[2] + 0.35*x[3] + x[4] + 3.36*x[6] - 6.3*x[5]*x[8])
@objective(m, Min, x[2])
println("start opt")
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

final_time = get_solution_time(history)[end]
println("final time: $final_time")
