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
x_Idx = Any[2, 3, 4, 5, 6, 7, 8, 9]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[5], 0.0)
JuMP.set_lower_bound(x[9], 0.0)
JuMP.set_lower_bound(x[8], 0.0)
JuMP.set_lower_bound(x[7], 0.0)
JuMP.set_lower_bound(x[4], 0.0)
JuMP.set_lower_bound(x[6], 0.0)
JuMP.set_lower_bound(x[3], 0.0)
JuMP.set_lower_bound(x[2], 0.115)
JuMP.set_upper_bound(x[2], 0.115)
JuMP.set_upper_bound(x[3], 1.0)
JuMP.set_upper_bound(x[4], 1.0)
JuMP.set_upper_bound(x[5], 1.0)
JuMP.set_upper_bound(x[6], 1.0)
JuMP.set_upper_bound(x[7], 1.0)
JuMP.set_upper_bound(x[8], 1.0)
JuMP.set_upper_bound(x[9], 1.0)


# ----- Constraints ----- #
@constraint(m, e2, x[2]-0.1287*x[3]-0.1096*x[4]-0.0501*x[5]-0.1524*x[6]-0.0763*x[7]-0.1854*x[8]-0.062*x[9] == 0.0)
@constraint(m, e3, x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9] == 1.0)


# ----- Objective ----- #
@NLobjective(m, Min, 0.5*(42.18*x[3]*x[3]+20.18*x[3]*x[4]+10.88*x[3]*x[5]+5.3*x[3]*x[6]+12.32*x[3]*x[7]+23.84*x[3]*x[8]+17.41*x[3]*x[9]+20.18*x[4]*x[3]+70.89*x[4]*x[4]+21.58*x[4]*x[5]+15.41*x[4]*x[6]+23.24*x[4]*x[7]+23.8*x[4]*x[8]+12.62*x[4]*x[9]+10.88*x[5]*x[3]+21.58*x[5]*x[4]+25.51*x[5]*x[5]+9.6*x[5]*x[6]+22.63*x[5]*x[7]+13.22*x[5]*x[8]+4.7*x[5]*x[9]+5.3*x[6]*x[3]+15.41*x[6]*x[4]+9.6*x[6]*x[5]+22.33*x[6]*x[6]+10.32*x[6]*x[7]+10.46*x[6]*x[8]+x[6]*x[9]+12.32*x[7]*x[3]+23.24*x[7]*x[4]+22.63*x[7]*x[5]+10.32*x[7]*x[6]+30.01*x[7]*x[7]+16.36*x[7]*x[8]+7.2*x[7]*x[9]+23.84*x[8]*x[3]+23.8*x[8]*x[4]+13.22*x[8]*x[5]+10.46*x[8]*x[6]+16.36*x[8]*x[7]+42.23*x[8]*x[8]+9.9*x[8]*x[9]+17.41*x[9]*x[3]+12.62*x[9]*x[4]+4.7*x[9]*x[5]+x[9]*x[6]+7.2*x[9]*x[7]+9.9*x[9]*x[8]+16.42*x[9]*x[9]))


println("start optimize")
JuMP.optimize!(m)
println("finish optimize")

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
evaluator_constraints_local = optimizer.nlp_data.evaluator.constraints		 # model get returned when including this script.
