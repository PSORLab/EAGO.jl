
using JuMP, EAGO, Cbc

function make_model()
    sb = SubSolvers(;relaxed_optimizer= Cbc.Optimizer())
    m = Model(() -> EAGO.Optimizer(subsolver_block = sb))
    set_optimizer_attribute(m, "output_iterations", 1)
    set_optimizer_attribute(m, "iteration_limit", 5)
    set_optimizer_attribute(m, "verbosity", 0)

    # Define bounded variables
    xL = [10.0; 0.0; 0.0; 0.0; 0.0; 85.0; 90.0; 3.0; 1.2; 145.0]
    xU = [2000.0; 16000.0; 120.0; 5000.0; 2000.0; 93.0; 95.0; 12.0; 4.0; 162.0]
    @variable(m, xL[1] <= x[1] <= xU[1])
    @NLobjective(m, Max, -2.0 * x[1])
    return m
end

function f(m)
    JuMP.optimize!(m)
end

m = make_model()
m1 = make_model()

println("EAGO first call ")
@time f(m)
@time f(m1)

#=
m1 = Model(EAGO.Optimizer)
set_optimizer_attribute(m1, "output_iterations", 1)
set_optimizer_attribute(m1, "iteration_limit", 5)
set_optimizer_attribute(m1, "verbosity", 0)

# Define bounded variables
xL = [10.0; 0.0; 0.0; 0.0; 0.0; 85.0; 90.0; 3.0; 1.2; 145.0]
xU = [2000.0; 16000.0; 120.0; 5000.0; 2000.0; 93.0; 95.0; 12.0; 4.0; 162.0]
@variable(m1, xL[1] <= x[1] <= xU[1])

@constraint(m, e2, -x[1]+1.22*x[4]-x[5] == 0.0)
@constraint(m, e6, x[9]+0.222*x[10] == 35.82)
@constraint(m, e7, -3*x[7]+x[10] == -133.0)

# Define nonlinear constraints
@NLconstraint(m, e1, -x[1]*(1.12+0.13167*x[8]-0.00667* (x[8])^2)+x[4] == 0.0)
@NLconstraint(m, e3, -0.001*x[4]*x[9]*x[6]/(98-x[6])+x[3] == 0.0)

###Original constraints
@NLconstraint(m, e4, -(1.098*x[8]-0.038* (x[8])^2)-0.325*x[6]+x[7] == 57.425)
@NLconstraint(m, e5, -(x[2]+x[5])/x[1]+x[8] == 0.0)

@NLobjective(m1, Max, -2.0 * x[1])
JuMP.optimize!(m1)
=#
#=
=#
# Define nonlinear constraints
#@NLconstraint(m, e1, -x[1]*(1.12+0.13167*x[8]-0.00667* (x[8])^2)+x[4] == 0.0)
#@NLconstraint(m, e3, -0.001*x[4]*x[9]*x[6]/(98-x[6])+x[3] == 0.0)

###Original constraints
#@NLconstraint(m, e4, -(1.098*x[8]-0.038* (x[8])^2)-0.325*x[6]+x[7] == 57.425)
#@NLconstraint(m, e5, -(x[2]+x[5])/x[1]+x[8] == 0.0)

###Rewritten constraints
#=
@NLexpressions(m, begin
 ex1, -(x[2]+x[5])/x[1]+x[8]
 ex2, -(1.098*x[8]-0.038* (x[8])^2)-0.325*x[6]+x[7]
end)
@NLconstraint(m, e4, ex2 == 57.425)
@NLconstraint(m, e5, ex1 == 0.0)
=#

# Define linear constraints
#@constraint(m, e2, -x[1]+1.22*x[4]-x[5] == 0.0)
#@constraint(m, e6, x[9]+0.222*x[10] == 35.82)
#@constraint(m, e7, -3*x[7]+x[10] == -133.0)

# Define nonlinear objective
#@NLobjective(m, Max, 0.063*x[4]*x[7] - 5.04*x[1] - 0.035*x[2] - 10*x[3] - 3.36*x[5])
#@NLobjective(m, Min, -(x[4]*x[7]))
#@NLobjective(m, Min, (0.063*x[4]*x[7]^2 - 5.04*x[1]^2 - 0.035*x[2] - 10*x[3] - 3.36*x[5]))

#=
@NLobjective(m, Max, -2.0*x[1])


JuMP.optimize!(m)
# Solve the optimization problem
#using ProfileView
#@profview f(m)
#@profview f(m)

@show objective_value(m)
@show objective_bound(m)
@show value(x[1])
=#
