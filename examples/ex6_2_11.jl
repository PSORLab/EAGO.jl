
println("start run")

using JuMP, EAGO

m = Model(EAGO.Optimizer)
set_optimizer_attribute(m, "mul_relax_style", 1)
set_optimizer_attribute(m, "verbosity", 1)
set_optimizer_attribute(m, "output_iterations", 1000)
set_optimizer_attribute(m, "iteration_limit", 10000000)
set_optimizer_attribute(m, "cut_max_iterations", 2)

# OBBT depth 0 -> 20... increases number of iterations...
set_optimizer_attribute(m, "obbt_depth", 8)
set_optimizer_attribute(m, "obbt_repetitions", 2)

# ----- Variables ----- #
x_Idx = Any[2, 3, 4]
@variable(m, x[x_Idx])
set_lower_bound(x[2], 1.0e-6)
set_upper_bound(x[2], 1.0)
set_lower_bound(x[3], 1.0e-6)
set_upper_bound(x[3], 1.0)
set_lower_bound(x[4], 1.0e-6)
set_upper_bound(x[4], 1.0)


# ----- Constraints ----- #
s = time()
#=
@NLobjective(m, Min, (15.3261663216011*x[2]+23.2043471859416*x[3]+6.69678129464404*x[4])*log(2.1055*x[2]+3.1878*x[3] + 0.92*x[4]) + 
                      1.04055250396734*x[2]-2.24199441248417*x[3]+3.1618173099828*x[4]+6.4661663216011*x[2]*log(x[2]/(2.1055*x[2]+3.1878*x[3]+0.92*x[4])) + 
                      12.2043471859416*x[3]*log(x[3]/(2.1055*x[2]+3.1878*x[3]+0.92*x[4]))+0.696781294644034*x[4]*log(x[4]/(2.1055*x[2]+3.1878*x[3]+0.92*x[4])) + 
                      9.86*x[2]*log(x[2]/(1.972*x[2]+2.4*x[3]+1.4*x[4]))+12*x[3]*log(x[3]/(1.972*x[2]+2.4*x[3]+1.4*x[4]))+7*x[4]*log(x[4]/(1.972*x[2]+2.4*x[3] + 
                      1.4*x[4]))+(1.972*x[2]+2.4*x[3]+1.4*x[4])*log(1.972*x[2]+2.4*x[3]+1.4*x[4])+1.972*x[2]*log(x[2]/(1.972*x[2]+0.283910843616504*x[3] + 
                      3.02002220174195*x[4]))+2.4*x[3]*log(x[3]/(1.45991339466884*x[2]+2.4*x[3]+0.415073537580851*x[4]))+1.4*x[4]*log(x[4]/(0.602183324335333*x[2] + 
                      0.115623371371275*x[3]+1.4*x[4]))-17.2981663216011*x[2]*log(x[2])-25.6043471859416*x[3]*log(x[3])-8.09678129464404*x[4]*log(x[4]))
=#
#@NLobjective(m, Min, (x[3]^3 - x[3]^3)*(x[3]^3 - x[3]^3))
@NLobjective(m, Min, x[2]*x[3]*x[4] + x[2]*x[4])
@constraint(m, e2, x[2]+x[3]+x[4] == 1.0)

m = m 		 # model get returned when including this script. 

optimize!(m)
@show termination_status(m)
sout = s - time()
@show solve_time(m)
@show sout