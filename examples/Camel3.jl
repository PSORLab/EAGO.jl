using JuMP, EAGO

m = Model(EAGO.Optimizer)
set_optimizer_attribute(m, "mul_relax_style", 1)
set_optimizer_attribute(m, "verbosity", 4)
set_optimizer_attribute(m, "output_iterations", 1)
set_optimizer_attribute(m, "iteration_limit", 1000)
set_optimizer_attribute(m, "cut_max_iterations", 2)
set_optimizer_attribute(m, "subgrad_tighten", false)

# OBBT depth 0 -> 20... increases number of iterations...
set_optimizer_attribute(m, "obbt_depth", 8)
set_optimizer_attribute(m, "obbt_repetitions", 2)

# ----- Variables ----- #
x_Idx = Any[1,2]
@variable(m, x[x_Idx])
set_lower_bound(x[1], 0.1)
set_upper_bound(x[1], 0.9)

set_lower_bound(x[2], 0.1)
set_upper_bound(x[2], 0.9)

# ----- Constraints ----- #

@NLobjective(m, Min, (x[1]^2 - x[1])*(x[2]^2 - x[2]))
s = time()
optimize!(m)
@show termination_status(m)
sout = s - time()
@show solve_time(m)
@show sout
