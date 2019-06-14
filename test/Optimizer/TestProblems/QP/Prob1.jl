#m = Model(with_optimizer(Ipopt.Optimizer))
m = Model(with_optimizer(EAGO.Optimizer))

# Need nonnegativity for (rotated) second-order cone
@variable(m, 0 <= x <= 1)
@variable(m, 0 <= y <= 1)
@variable(m, 0 <= z <= 1)
#@variable(m, 0.496094 <= x <= 1)
#@variable(m, 0 <= y <= 1)
#@variable(m, 0 <= z <= 1)



# Maximize x
@objective(m, Max, x)
#@objective(m, Min, -x)


# Subject to 1 linear and 2 nonlinear constraints
#@constraint(m, x + y + z == 1)
#@constraint(m, x*x + y*y - z*z <= 0)
#@constraint(m, x*x - y*z <= 0)
@constraint(m, x + y + z == 1)
@constraint(m, x*x + y*y - z*z <= 0)
@constraint(m, x*x - y*z <= 0)

# Solve with Gurobi
JuMP.optimize!(m)

# Solution
println("Objective Value: ", JuMP.objective_value(m))
println("Termination Status: ", JuMP.termination_status(m))
println("Primal Status: ", JuMP.primal_status(m))
println("x = ", JuMP.value(x))
println("y = ", JuMP.value(y))
println("y = ", JuMP.value(z))
#=
@testset "QP Problem #1" begin

    m = Model(with_optimizer(EAGO.Optimizer))

    # Need nonnegativity for (rotated) second-order cone
    @variable(m, x)
    @variable(m, y >= 0)
    @variable(m, z >= 0)

    # Maximize x
    @objective(m, Max, x)

    # Subject to 1 linear and 2 nonlinear constraints
    @constraint(m, x + y + z == 1)
    @constraint(m, x*x + y*y - z*z <= 0)
    @constraint(m, x*x - y*z <= 0)

    # Print the model to check correctness
    print(m)

    # Solve with Gurobi
    JuMP.optimize!(m)

    # Solution
    println("Objective value: ", JuMP.objective_value(m))
    println("x = ", JuMP.value(x))
    println("y = ", JuMP.value(y))

end
=#
