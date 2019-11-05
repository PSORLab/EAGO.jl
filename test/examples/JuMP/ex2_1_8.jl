
m = Model(with_optimizer(EAGO.Optimizer,
                             lp_depth = 100000,
                             lp_reptitions = 3,
                             quad_uni_depth = -1,
                             obbt_depth = 6,
                             cp_reptitions = 3,
                             cp_depth = 1000,
                             iteration_limit = 10000000,
                             verbosity = 0,
                             output_iterations = 1,
                             header_iterations = 400000,
                             relative_tolerance = 1E-3,
                             absolute_tolerance = 1E-3,
                             dbbt_depth = 100000000,
                             subgrad_tighten = true, #true,
                             objective_cut_on = true,
                             cut_max_iterations = 3,
                             upper_bounding_depth = 4,
                             time_limit = 1000.0,
                             relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))

                             x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
                             @variable(m, 0 <= x[x_Idx] <= 100)

                             # ----- Constraints ----- #
                             @constraint(m, e2, x[1]+x[2]+x[3]+x[4] == 8.0)
                             @constraint(m, e3, x[5]+x[6]+x[7]+x[8] == 24.0)
                             @constraint(m, e4, x[9]+x[10]+x[11]+x[12] == 20.0)
                             @constraint(m, e5, x[13]+x[14]+x[15]+x[16] == 24.0)
                             @constraint(m, e6, x[17]+x[18]+x[19]+x[20] == 16.0)
                             @constraint(m, e7, x[21]+x[22]+x[23]+x[24] == 12.0)
                             @constraint(m, e8, x[1]+x[5]+x[9]+x[13]+x[17]+x[21] == 29.0)
                             @constraint(m, e9, x[2]+x[6]+x[10]+x[14]+x[18]+x[22] == 41.0)
                             @constraint(m, e10, x[3]+x[7]+x[11]+x[15]+x[19]+x[23] == 13.0)
                             @constraint(m, e11, x[4]+x[8]+x[12]+x[16]+x[20]+x[24] == 21.0)

@NLobjective(m, Min, (300*x[1]-7*x[1]*x[1]-4*x[2]*x[2]+270*x[2]-6*x[3]*x[3]+460*x[3]-8*x[4]*x[4]+800*x[4]-12*x[5]*x[5]+740*x[5]-9*x[6]*x[6]+600*x[6]-14*x[7]*x[7]+540*x[7]-7*x[8]*x[8]+380*x[8]-13*x[9]*x[9]+300*x[9]-12*x[10]*x[10]+490*x[10]-8*x[11]*x[11]+380*x[11]-4*x[12]*x[12]+760*x[12]-7*x[13]*x[13]+430*x[13]-9*x[14]*x[14]+250*x[14]-16*x[15]*x[15]+390*x[15]-8*x[16]*x[16]+600*x[16]-4*x[17]*x[17]+210*x[17]-10*x[18]*x[18]+830*x[18]-21*x[19]*x[19]+470*x[19]-13*x[20]*x[20]+680*x[20]-17*x[21]*x[21]+360*x[21]-9*x[22]*x[22]+290*x[22]-8*x[23]*x[23]+400*x[23]-4*x[24]*x[24]+310*x[24]))

JuMP.optimize!(m)
run_time = backend(m).optimizer.model.optimizer._run_time
#objective_value =
#objective_bound =
#println("objective between: ")
println("run time: $run_time")
