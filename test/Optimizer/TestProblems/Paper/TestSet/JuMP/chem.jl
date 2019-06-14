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
                                     cut_max_iterations = 0,
                                     upper_bounding_interval = 2)

m = Model(opt)

# ----- Variables ----- #
x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
@variable(m, x[x_Idx])
setlowerbound(x[1], 0.001)
setlowerbound(x[2], 0.001)
setlowerbound(x[3], 0.001)
setlowerbound(x[4], 0.001)
setlowerbound(x[5], 0.001)
setlowerbound(x[6], 0.001)
setlowerbound(x[7], 0.001)
setlowerbound(x[8], 0.001)
setlowerbound(x[9], 0.001)
setlowerbound(x[10], 0.001)
setlowerbound(x[11], 0.01)


# ----- Constraints ----- #
@constraint(m, e1, x[1]+2*x[2]+2*x[3]+x[6]+x[10] == 2.0)
@constraint(m, e2, x[4]+2*x[5]+x[6]+x[7] == 1.0)
@constraint(m, e3, x[3]+x[7]+x[8]+2*x[9]+x[10] == 1.0)
@constraint(m, e5, -x[1]-x[2]-x[3]-x[4]-x[5]-x[6]-x[7]-x[8]-x[9]-x[10]+x[11] == 0.0)


# ----- Objective ----- #
@NLobjective(m, Min, (x[1]*(log(x[1]/x[11])-6.05576803624071)+x[2]*(log(x[2]/x[11])-17.1307680362407)+x[3]*(log(x[3]/x[11])-34.0207680362407)+x[4]*(log(x[4]/x[11])-5.88076803624071)+x[5]*(log(x[5]/x[11])-24.6877680362407)+x[6]*(log(x[6]/x[11])-14.9527680362407)+x[7]*(log(x[7]/x[11])-24.0667680362407)+x[8]*(log(x[8]/x[11])-10.6747680362407)+x[9]*(log(x[9]/x[11])-26.6287680362407)+x[10]*(log(x[10]/x[11])-22.1447680362407)))

m = m 		 # model get returned when including this script.
