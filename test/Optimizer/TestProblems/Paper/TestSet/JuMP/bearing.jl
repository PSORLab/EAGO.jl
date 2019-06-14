using JuMP, EAGO

opt = with_optimizer(EAGO.Optimizer, cp_depth = 0, cp_reptitions = 0,
                                     obbt_depth = 4,
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
                                     verbosity = 0,
                                     header_iterations = 1000,
                                     output_iterations = 200,
                                     cut_max_iterations = 3,
                                     upper_bounding_interval = 2)

m = Model(opt)

# ----- Variables ----- #
x_Idx = Any[1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14]
#x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]

@variable(m, x[x_Idx])

JuMP.set_lower_bound(x[1], 1.0)
JuMP.set_upper_bound(x[1], 16.0)

JuMP.set_lower_bound(x[2], 1.0)
JuMP.set_upper_bound(x[2], 16.0)

JuMP.set_lower_bound(x[3], 1.0)
JuMP.set_upper_bound(x[3], 16.0)

JuMP.set_lower_bound(x[4], 1.0)
JuMP.set_upper_bound(x[4], 16.0)

JuMP.set_lower_bound(x[6], 1.0)
JuMP.set_upper_bound(x[6], 1000.0)            # stored in 5th position

JuMP.set_lower_bound(x[7], 0.0001)
JuMP.set_upper_bound(x[7], 2.0000)            # stored in 6th position

JuMP.set_lower_bound(x[8], 0.0001)
JuMP.set_upper_bound(x[8], 2.0000)            # stored in 7th position

JuMP.set_lower_bound(x[9], 1.0)
JuMP.set_upper_bound(x[9], 2.0)               # stored in 8th position

JuMP.set_lower_bound(x[10], 0.0)
JuMP.set_upper_bound(x[10], 50.0)

JuMP.set_lower_bound(x[11], 550.0)
JuMP.set_upper_bound(x[11], 600.0)

JuMP.set_lower_bound(x[12], 1.0)
JuMP.set_upper_bound(x[12], 3.0)

JuMP.set_lower_bound(x[13], 0.0001)
JuMP.set_upper_bound(x[13], 3.0)

JuMP.set_lower_bound(x[14], 0.01)
JuMP.set_upper_bound(x[14], 10.0)


# ----- Constraints ----- #

@NLconstraint(m, e2, -1.42857142857143*x[4]*x[6] + 10000*x[8] == 0.0)                      # multiplication
#@NLconstraint(m, e3, 10*x[7]*x[9] - 0.00968946189201592*(x[1]^4 - x[2]^4)*x[3] == 0.0)     # multiplication, ^4 (Ill conditioned)
#@NLconstraint(m, e4, 143.3076*x[10]*x[4] - 10000*x[7] == 0.0)                              # multiplication
#@NLconstraint(m, e5, 3.1415927*(0.001*x[9])^3*x[6] - 6e-6*x[3]*x[4]*x[13] == 0.0)          # multiplication, ^3
#@NLconstraint(m, e6, 101000*x[12]*x[13] - 1.57079635*x[6]*x[14] == 0.0)                    # multiplication
#@NLconstraint(m, e7, log10(0.8 + 8.112*x[3]) - 10964781961.4318*x[11]^(-3.55) == 0.0)      # multiplication, log10, ^(float)
@constraint(m, e8, -0.5*x[10] + x[11] == 560.0)                                            # linear
@constraint(m, e9, x[1] - x[2] >= 0.0)                                                     # linear
@NLconstraint(m, e10, 0.0307*(x[4])^2 - 0.3864*(0.0062831854*x[1]*x[9])^2*x[6] <= 0.0)     # multiplication, ^2
@constraint(m, e11, 101000*x[12] - 15707.9635*x[14] <= 0.0)                                # linear
@NLconstraint(m, e12, -(log(x[1]) - log(x[2])) + x[13] == 0.0)                             # log
@NLconstraint(m, e13, -( (x[1])^2 - (x[2])^2) + x[14] == 0.0)                              # ^2


# ----- Objective ----- #
@objective(m, Min, x[7] + x[8])

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


#runtime = history.lower_time[11] + history.upper_time[11] + history.preprocess_time[11] + history.postprocess_time[11]
#println("Run time is $runtime seconds.")
