@testset "Set Objective" begin
    model = EAGO.Optimizer()

    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test model.optimization_sense == MOI.MIN_SENSE

    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    @test model.optimization_sense == MOI.MAX_SENSE

    MOI.set(model, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    @test model.optimization_sense == MOI.FEASIBILITY_SENSE
end

@testset "Get Termination Code " begin

    model = EAGO.Optimizer()

    # Termination Status Code Checks
    model.termination_status_code = MOI.OPTIMAL
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == MOI.OPTIMAL

    model.termination_status_code = MOI.ITERATION_LIMIT
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == MOI.ITERATION_LIMIT
end

@testset "Add Variable, Get Number/Index" begin

    model = EAGO.Optimizer()

    # Add variable
    MOI.add_variable(model)
    nVar = MOI.get(model, MOI.NumberOfVariables())
    @test nVar == 1

    # Add second variable
    MOI.add_variables(model,3)
    nVar = MOI.get(model, MOI.NumberOfVariables())
    @test nVar == 4

    # Get variable indices
    indx = MOI.get(model, MOI.ListOfVariableIndices())
    @test indx ==  MOI.VariableIndex[MOI.VariableIndex(1), MOI.VariableIndex(2),
                                    MOI.VariableIndex(3), MOI.VariableIndex(4)]

    @test_nowarn EAGO.check_inbounds(model, MOI.VariableIndex(1))
    @test_throws ErrorException EAGO.check_inbounds(model,MOI.VariableIndex(6))

end

@testset "Add Variable Bounds" begin
    model = EAGO.Optimizer()

    x = MOI.add_variables(model,3)
    z = MOI.add_variable(model)

    MOI.add_constraint(model, MOI.SingleVariable(x[1]), MOI.GreaterThan(-1.0))
    MOI.add_constraint(model, MOI.SingleVariable(x[2]), MOI.LessThan(-1.0))
    MOI.add_constraint(model, MOI.SingleVariable(x[3]), MOI.EqualTo(2.0))
    MOI.add_constraint(model, MOI.SingleVariable(z), MOI.ZeroOne())

    @test model.variable_info[1].is_integer == false
    @test model.variable_info[1].lower_bound == -1.0
    @test model.variable_info[1].has_lower_bound == true
    @test model.variable_info[1].upper_bound == Inf
    @test model.variable_info[1].has_upper_bound == false
    @test model.variable_info[1].is_fixed == false

    @test model.variable_info[2].is_integer == false
    @test model.variable_info[2].lower_bound == -Inf
    @test model.variable_info[2].has_lower_bound == false
    @test model.variable_info[2].upper_bound == -1.0
    @test model.variable_info[2].has_upper_bound == true
    @test model.variable_info[2].is_fixed == false

    @test model.variable_info[3].is_integer == false
    @test model.variable_info[3].lower_bound == 2.0
    @test model.variable_info[3].has_lower_bound == true
    @test model.variable_info[3].upper_bound == 2.0
    @test model.variable_info[3].has_upper_bound == true
    @test model.variable_info[3].is_fixed == true

    @test model.variable_info[4].is_integer == true
    @test model.variable_info[4].lower_bound == 0.0
    @test model.variable_info[4].has_lower_bound == true
    @test model.variable_info[4].upper_bound == 1.0
    @test model.variable_info[4].has_upper_bound == true
    @test model.variable_info[4].is_fixed == false
end

@testset "Add Linear Constraint " begin

    model = EAGO.Optimizer()

    x = MOI.add_variables(model,3)

    func1 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([5.0,-2.3],[x[1],x[2]]),2.0)
    func2 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([4.0,-2.2],[x[2],x[3]]),2.1)
    func3 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([3.0,-3.3],[x[1],x[3]]),2.2)

    set1 = MOI.LessThan{Float64}(1.0)
    set2 = MOI.GreaterThan{Float64}(2.0)
    set3 = MOI.EqualTo{Float64}(3.0)

    MOI.add_constraint(model, func1, set1)
    MOI.add_constraint(model, func2, set2)
    MOI.add_constraint(model, func3, set3)

    @test model.linear_leq_constraints[1][1].constant == 2.0
    @test model.linear_geq_constraints[1][1].constant == 2.1
    @test model.linear_eq_constraints[1][1].constant == 2.2
    @test model.linear_leq_constraints[1][1].terms[1].coefficient == 5.0
    @test model.linear_geq_constraints[1][1].terms[1].coefficient == 4.0
    @test model.linear_eq_constraints[1][1].terms[1].coefficient == 3.0
    @test model.linear_leq_constraints[1][1].terms[2].coefficient == -2.3
    @test model.linear_geq_constraints[1][1].terms[2].coefficient == -2.2
    @test model.linear_eq_constraints[1][1].terms[2].coefficient == -3.3
    @test model.linear_leq_constraints[1][1].terms[1].variable_index.value == 1
    @test model.linear_geq_constraints[1][1].terms[1].variable_index.value == 2
    @test model.linear_eq_constraints[1][1].terms[1].variable_index.value == 1
    @test model.linear_leq_constraints[1][1].terms[2].variable_index.value == 2
    @test model.linear_geq_constraints[1][1].terms[2].variable_index.value == 3
    @test model.linear_eq_constraints[1][1].terms[2].variable_index.value == 3
    @test MOI.LessThan{Float64}(1.0) == model.linear_leq_constraints[1][2]
    @test MOI.GreaterThan{Float64}(2.0) == model.linear_geq_constraints[1][2]
    @test MOI.EqualTo{Float64}(3.0) == model.linear_eq_constraints[1][2]
    @test model.linear_leq_constraints[1][3] == 2
    @test model.linear_geq_constraints[1][3] == 2
    @test model.linear_eq_constraints[1][3] == 2
end

@testset "Add Quadratic Constraint " begin

    model = EAGO.Optimizer()

    x = MOI.add_variables(model,3)

    func1 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(5.0,x[1])],
                                                 [MOI.ScalarQuadraticTerm{Float64}(2.5,x[2],x[2])],2.0)
    func2 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(4.0,x[2])],
                                                 [MOI.ScalarQuadraticTerm{Float64}(2.2,x[1],x[2])],2.1)
    func3 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(3.0,x[3])],
                                                 [MOI.ScalarQuadraticTerm{Float64}(2.1,x[1],x[1])],2.2)

    set1 = MOI.LessThan{Float64}(1.0)
    set2 = MOI.GreaterThan{Float64}(2.0)
    set3 = MOI.EqualTo{Float64}(3.0)

    MOI.add_constraint(model, func1, set1)
    MOI.add_constraint(model, func2, set2)
    MOI.add_constraint(model, func3, set3)

    @test model.quadratic_leq_constraints[1][1].constant == 2.0
    @test model.quadratic_geq_constraints[1][1].constant == 2.1
    @test model.quadratic_eq_constraints[1][1].constant == 2.2
    @test model.quadratic_leq_constraints[1][1].quadratic_terms[1].coefficient == 2.5
    @test model.quadratic_geq_constraints[1][1].quadratic_terms[1].coefficient == 2.2
    @test model.quadratic_eq_constraints[1][1].quadratic_terms[1].coefficient == 2.1
    @test model.quadratic_leq_constraints[1][1].affine_terms[1].coefficient == 5.0
    @test model.quadratic_geq_constraints[1][1].affine_terms[1].coefficient == 4.0
    @test model.quadratic_eq_constraints[1][1].affine_terms[1].coefficient == 3.0
    @test model.quadratic_leq_constraints[1][1].quadratic_terms[1].variable_index_1.value == 2
    @test model.quadratic_geq_constraints[1][1].quadratic_terms[1].variable_index_1.value == 1
    @test model.quadratic_eq_constraints[1][1].quadratic_terms[1].variable_index_1.value == 1
    @test model.quadratic_leq_constraints[1][1].quadratic_terms[1].variable_index_2.value == 2
    @test model.quadratic_geq_constraints[1][1].quadratic_terms[1].variable_index_2.value == 2
    @test model.quadratic_eq_constraints[1][1].quadratic_terms[1].variable_index_2.value == 1
    @test model.quadratic_leq_constraints[1][1].affine_terms[1].variable_index.value == 1
    @test model.quadratic_geq_constraints[1][1].affine_terms[1].variable_index.value == 2
    @test model.quadratic_eq_constraints[1][1].affine_terms[1].variable_index.value == 3
    @test MOI.LessThan{Float64}(1.0) == model.quadratic_leq_constraints[1][2]
    @test MOI.GreaterThan{Float64}(2.0) == model.quadratic_geq_constraints[1][2]
    @test MOI.EqualTo{Float64}(3.0) == model.quadratic_eq_constraints[1][2]
    @test model.quadratic_leq_constraints[1][3] == 1
    @test model.quadratic_geq_constraints[1][3] == 1
    @test model.quadratic_eq_constraints[1][3] == 1
end

@testset "Empty/Isempty, EaGO Model " begin
    model = EAGO.Optimizer()
    @test MOI.is_empty(model)
end

println("BEGIN TESTING QUADRATIC RELAXATIONS...")
#include("quadratic_relaxation.jl")
println("TESTING QUADRATIC RELAXATIONS COMPLETE.")

println("BEGIN TESTING STANDARD EVALUATOR...")
include("standard_evaluator.jl")
println("TESTING STANDARD EVALUATOR COMPLETE.")

println("BEGIN TESTING IMPLICIT EVALUATOR...")
#include("implicit_optimizer.jl")
println("TESTING IMPLICIT EVALUATOR COMPLETE.")

@testset "LP Problems" begin
    m = Model(with_optimizer(EAGO.Optimizer))
    #m = Model(with_optimizer(Clp.Optimizer))

    @variable(m, 1 <= x <= 3)
    @variable(m, 1 <= y <= 3)

    @NLobjective(m, Min, x + y)

    @NLconstraint(m, x + y <= 10)
    @NLconstraint(m, x - y <= 10)
    @NLconstraint(m, y >= 0)

    JuMP.optimize!(m)

    xval = JuMP.value(x)
    yval = JuMP.value(y)

    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(xval,1.0,atol=1E-4)
    @test isapprox(yval,1.0,atol=1E-4)
    @test isapprox(fval,2.0,atol=1E-4)
    @test status_term == MOI.OPTIMAL
    @test status_prim == MOI.FEASIBLE_POINT

    m = Model(with_optimizer(EAGO.Optimizer))

    @variable(m, -3 <= x <= -1)
    @variable(m, -2 <= y <= 2)
    @variable(m, 1 <= z <= 3)

    @NLobjective(m, Min, x - y + 2z)

    @NLconstraint(m, x + 2y >= -10)
    @NLconstraint(m, z - 2y <= 2)
    @NLconstraint(m, y >= 0)

    JuMP.optimize!(m)

    xval = JuMP.value(x)
    yval = JuMP.value(y)
    zval = JuMP.value(z)
    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(xval,-3.0,atol=1E-4)
    @test isapprox(yval,2.0,atol=1E-4)
    @test isapprox(zval,1.0,atol=1E-4)
    @test isapprox(fval,-3.0,atol=1E-4)
    @test status_term == MOI.OPTIMAL
    @test status_prim == MOI.FEASIBLE_POINT

    m = Model(with_optimizer(EAGO.Optimizer))
    #m = Model(with_optimizer(Clp.Optimizer))

    @variable(m, -3 <= x <= -1)
    @variable(m, -2 <= y <= 2)
    @variable(m, 1 <= z <= 3)
    @variable(m, -10 <= q <= 9)

    @NLobjective(m, Min, 2x - 3y + 2z)

    @NLconstraint(m, x + 2y >= -10)
    @NLconstraint(m, z - 2y <= 2)
    @NLconstraint(m, y >= 0)
    @NLconstraint(m, q-3*z-y >= 0)

    JuMP.optimize!(m)

    xval = JuMP.value(x)
    yval = JuMP.value(y)
    zval = JuMP.value(z)
    qval = JuMP.value(q)

    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(xval,-3.0,atol=1E-4)
    @test isapprox(yval,2.0,atol=1E-4)
    @test isapprox(zval,1.0,atol=1E-4)
    @test isapprox(fval,-10.0,atol=1E-4)
    @test status_term == MOI.OPTIMAL
    @test status_prim == MOI.FEASIBLE_POINT

    m = Model(with_optimizer(EAGO.Optimizer))
    #m = Model(with_optimizer(Clp.Optimizer))

    @variable(m, -3 <= x <= -1)
    @variable(m, -2 <= y <= 2)
    @variable(m, 1 <= z <= 3)
    @variable(m, -10 <= q <= 9)

    @NLobjective(m, Min, 2x - 3y + 2z)

    @NLconstraint(m, x + 2y >= -10)
    @NLconstraint(m, z - 2y <= 2)
    @NLconstraint(m, y >= 4)
    @NLconstraint(m, q-3*z-y >= 0)

    JuMP.optimize!(m)

    xval = JuMP.value(x)
    yval = JuMP.value(y)
    zval = JuMP.value(z)
    qval = JuMP.value(q)

    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(xval,0.0,atol=1E-4)
    @test isapprox(yval,0.0,atol=1E-4)
    @test isapprox(zval,0.0,atol=1E-4)
    @test isapprox(qval,0.0,atol=1E-4)
    @test isapprox(fval,Inf,atol=1E-4)
    @test status_term == MOI.INFEASIBLE
    @test status_prim == MOI.INFEASIBILITY_CERTIFICATE
end


println("BEGIN NLP TEST PROBLEMS...")
for i in 1:5
   include("TestProblems/NLP/Prob$i.jl")
end
println("NLP TEST COMPLETE.")
