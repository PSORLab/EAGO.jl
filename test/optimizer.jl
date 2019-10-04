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

@testset "Empty/Isempty, EAGO Model " begin
    model = EAGO.Optimizer()
    @test MOI.is_empty(model)
end

#=
function test_forward_evaluator(m::JuMP.Model, expr::Expr, xpnt::Vector{Float64},
                                xL::Vector{Float64}, xU::Vector{Float64},
                                tfunc::Function; atol = 0.0)
    n = length(xpnt)
    JuMP.set_NL_objective(m, MOI.MIN_SENSE, expr)
    JuMP.add_NL_constraint(m, :($expr <= 0.0))
    source_evaluator = JuMP.NLPEvaluator(m)
    MOI.initialize(source_evaluator , Symbol[:Grad, :Jac])

    opt = m.moi_backend.optimizer.model.optimizer
    built_evaluator = EAGO.build_nlp_evaluator(2, source_evaluator, opt, true)

    # Add current node and define point
    built_evaluator.current_node = EAGO.NodeBB(xL, xU, -Inf, Inf, 2, 1, true)
    user_operators = built_evaluator.m.nlp_data.user_operators
    nlp_data = built_evaluator.m.nlp_data
    user_input_buffer = built_evaluator.jac_storage

    EAGO.forward_eval_all(built_evaluator, xpnt)

    cv_val_eval = MOI.eval_objective(built_evaluator, xpnt)
    cv_grad_eval = zeros(Float64, (n,))
    MOI.eval_objective_gradient(built_evaluator, cv_grad_eval, xpnt)
    Jstor = zeros(Float64, (1,n))
    MOI.eval_constraint_jacobian(built_evaluator, Jstor, xpnt)

    xMC = MC{n}.(xpnt, IntervalType.(xL,xU), [i for i in 1:n])
    zMC = tfunc(xMC)


    pass_flag::Bool = true
    descr = "Failing Components: ("
    ~isapprox(zMC.cv, cv_val_eval, atol = atol) &&  (descr = descr*" (obj CV, $(zMC.cv) !== $(cv_val_eval))"; pass_flag = false)
    for i in 1:n
        ~isapprox(zMC.cv_grad[i], cv_grad_eval[i], atol = atol) &&  (descr = descr*" (obj CVgrad[$i], $(zMC.cv_grad[i]) !== $(cv_grad_eval[i]))"; pass_flag = false)
    end
    for i in 1:n
        ~isapprox(zMC.cv_grad[i], Jstor[1,i], atol = atol) &&  (descr = descr*" (jac CVgrad[$i], $(zMC.cv_grad[i]) !== $(Jstor[1,i]))"; pass_flag = false)
    end
    (descr !== "Failing Components: (") && println(descr*")")
    pass_flag
end


@testset "NLP Evaluator" begin

        m = Model(with_optimizer(EAGO.Optimizer))
        @variable(m, x[1:2])
        expr = :($(x[1]) + $(x[2])^2)
        test_func = x -> (x[1]) + (x[2])^2
        xL = [1.0; 2.0]
        xU = [5.0; 6.0]
        midx = (xL + xU)/2.0
        @test test_forward_evaluator(m, expr, midx, xL, xU, test_func)

        m = Model(with_optimizer(EAGO.Optimizer))
        @variable(m, x[1:2])
        expr = :(exp(1.0 - $(x[1]))^2 + 100.0*($(x[2]) - $(x[1])^2)^2)
        test_func = x -> exp(1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2
        xL = [1.0; 5.0]
        xU = [2.0; 6.0]
        midx = (xL + xU)/2.0
        @test test_forward_evaluator(m, expr, midx, xL, xU, test_func)

        m = Model(with_optimizer(EAGO.Optimizer))
        @variable(m, x[1:2])
        expr = :($(x[2])/1.0)
        test_func = x -> x[2]/1.0
        xL = [1.0; 5.0]
        xU = [2.0; 6.0]
        midx = (xL + xU)/2.0
        @test test_forward_evaluator(m, expr, midx, xL, xU, test_func)

        m = Model(with_optimizer(EAGO.Optimizer))
        @variable(m, x[1:2])
        expr = :($(x[2])/$(x[1]))
        test_func = x -> x[2]/x[1]
        xL = [1.0; 5.0]
        xU = [2.0; 6.0]
        midx = (xL + xU)/2.0
        @test test_forward_evaluator(m, expr, midx, xL, xU, test_func)

        m = Model(with_optimizer(EAGO.Optimizer))
        @variable(m, x[1:2])
        expr = :($(x[2]) - $(x[1]))
        test_func = x -> x[2] - x[1]
        xL = [1.0; 5.0]
        xU = [2.0; 6.0]
        midx = (xL + xU)/2.0
        @test test_forward_evaluator(m, expr, midx, xL, xU, test_func)

        m = Model(with_optimizer(EAGO.Optimizer))
        @variable(m, x[1:2])
        expr = :(sin($(x[2])) - abs($(x[1])))
        test_func = x -> sin(x[2]) - abs(x[1])
        xL = [1.0; 5.0]
        xU = [2.0; 6.0]
        midx = (xL + xU)/2.0
        @test test_forward_evaluator(m, expr, midx, xL, xU, test_func)

        m = Model(with_optimizer(EAGO.Optimizer))
        @variable(m, x[1:2])
        expr = :(($(x[2]))^3 - ($(x[1]))^6)
        test_func = x -> x[2]^3 - x[1]^6
        xL = [1.0; 5.0]
        xU = [2.0; 6.0]
        midx = (xL + xU)/2.0
        @test test_forward_evaluator(m, expr, midx, xL, xU, test_func)

        m = Model(with_optimizer(EAGO.Optimizer))
        @variable(m, x[1:2])
        expr = :($(x[2])*$(x[1]))
        test_func = x -> x[2]*x[1]
        xL = [1.0; 5.0]
        xU = [2.0; 6.0]
        midx = (xL + xU)/2.0
        @test test_forward_evaluator(m, expr, midx, xL, xU, test_func, atol = 1E-8)

        m = Model(with_optimizer(EAGO.Optimizer))
        @variable(m, x[1:2])
        expr = :($(x[2])*$(x[1])*1.1*$(x[1]))
        test_func = x -> x[2]*x[1]*1.1*x[1]
        xL = [1.0; 5.0]
        xU = [2.0; 6.0]
        midx = (xL + xU)/2.0
        @test test_forward_evaluator(m, expr, midx, xL, xU, test_func, atol = 1E-8)

        m = Model(with_optimizer(EAGO.Optimizer))
        @variable(m, x[1:2])
        expr = :($(x[2])+$(x[1])+$(x[2])+1.2)
        test_func = x -> x[2]+x[1]+x[2]+1.2
        xL = [1.0; 5.0]
        xU = [2.0; 6.0]
        midx = (xL + xU)/2.0
        @test test_forward_evaluator(m, expr, midx, xL, xU, test_func, atol = 1E-8)
end
=#

@testset "LP Problems" begin
    m = Model(with_optimizer(EAGO.Optimizer,
                             udf_scrubber_flag = false,
                             udf_to_JuMP_flag = false,
                             verbosity = 0))

    @variable(m, 1 <= x <= 3)
    @variable(m, 1 <= y <= 3)

    @NLobjective(m, Min, x + y)

    @NLconstraint(m, x + y <= 10)
    @NLconstraint(m, x - y <= 10)
    @NLconstraint(m, y >= 0)

    JuMP.optimize!(m)

    @test isapprox(JuMP.value(x), 1.0, atol=1E-4)
    @test isapprox(JuMP.value(y), 1.0, atol=1E-4)
    @test isapprox(JuMP.objective_value(m), 2.0, atol=1E-4)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

    m = Model(with_optimizer(EAGO.Optimizer,
                             udf_scrubber_flag = false,
                             udf_to_JuMP_flag = false,
                             verbosity = 0))

    @variable(m, -3 <= x <= -1)
    @variable(m, -2 <= y <= 2)
    @variable(m, 1 <= z <= 3)

    @NLobjective(m, Min, x - y + 2z)

    @NLconstraint(m, x + 2y >= -10)
    @NLconstraint(m, z - 2y <= 2)
    @NLconstraint(m, y >= 0)

    JuMP.optimize!(m)

    @test isapprox(JuMP.value(x), -3.0, atol=1E-4)
    @test isapprox(JuMP.value(y), 2.0, atol=1E-4)
    @test isapprox(JuMP.value(z), 1.0, atol=1E-4)
    @test isapprox(JuMP.objective_value(m), -3.0, atol=1E-4)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

    m = Model(with_optimizer(EAGO.Optimizer,
                             udf_scrubber_flag = false,
                             udf_to_JuMP_flag = false,
                             verbosity = 0))

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

    @test isapprox(JuMP.value(x), -3.0, atol=1E-4)
    @test isapprox(JuMP.value(y), 2.0, atol=1E-4)
    @test isapprox(JuMP.value(z), 1.0, atol=1E-4)
    @test isapprox(JuMP.objective_value(m), -10.0, atol=1E-4)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

    m = Model(with_optimizer(EAGO.Optimizer,
                             udf_scrubber_flag = false,
                             udf_to_JuMP_flag = false,
                             verbosity = 0))

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

    @test isapprox(JuMP.value(x), 0.0, atol=1E-4)
    @test isapprox(JuMP.value(y), 0.0, atol=1E-4)
    @test isapprox(JuMP.value(z), 0.0, atol=1E-4)
    @test isapprox(JuMP.value(q), 0.0, atol=1E-4)
    @test isapprox(JuMP.objective_value(m), Inf, atol=1E-4)
    @test JuMP.termination_status(m) == MOI.INFEASIBLE
    @test JuMP.primal_status(m) == MOI.INFEASIBILITY_CERTIFICATE
end

#=
@testset "NLP Problems" begin
    m = Model(with_optimizer(EAGO.Optimizer,
                             udf_scrubber_flag = false,
                             udf_to_JuMP_flag = false,
                             verbosity = 0))
    @variable(m, -200 <= x <= -100)
    @variable(m, 200 <= y <= 400)
    @constraint(m, -500 <= x+2y <= 400)
    @NLobjective(m, Min, x*y)
    JuMP.optimize!(m)

    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.value(x), -200.0, atol=1E-5)
    @test isapprox(JuMP.value(y), 300.0, atol=1E-5)
    @test isapprox(JuMP.objective_value(m), -60000.00119999499, atol=2.0)

    m = Model(with_optimizer(EAGO.Optimizer,
                             udf_scrubber_flag = false,
                             udf_to_JuMP_flag = false,
                             verbosity = 0))
    x_Idx = Any[1, 2, 3]
    @variable(m, x[x_Idx])
    JuMP.set_lower_bound(x[1], -5.0)
    JuMP.set_upper_bound(x[1], 5.0)
    JuMP.set_lower_bound(x[2], -5.0)
    JuMP.set_upper_bound(x[2], 5.0)
    JuMP.set_lower_bound(x[3], -100.0)
    JuMP.set_upper_bound(x[3], 100.0)

    @NLconstraint(m, e2, 2* (x[2])^2+4*x[1]*x[2]-42*x[1]+4* (x[1])^3-x[3] <= 14.0)
    @NLconstraint(m, e3, (-2* (x[2])^2)-4*x[1]*x[2]+42*x[1]-4* (x[1])^3-x[3] <= -14.0)
    @NLconstraint(m, e4, 2* (x[1])^2+4*x[1]*x[2]-26*x[2]+4* (x[2])^3-x[3] <= 22.0)
    @NLconstraint(m, e5, (-2* (x[1])^2)-4*x[1]*x[2]+26*x[2]-4* (x[2])^3-x[3] <= -22.0)
    @NLobjective(m, Min, x[3])
    JuMP.optimize!(m)

    @test isapprox(JuMP.objective_value(m), 0.000, atol=1E-3)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT


    m = Model(with_optimizer(EAGO.Optimizer,
                             udf_scrubber_flag = false,
                             udf_to_JuMP_flag = false,
                             verbosity = 0))

    @variable(m, y)
    x_Idx = Any[1, 2]
    @variable(m, x[x_Idx])
    JuMP.set_lower_bound(x[1], 0.0)
    JuMP.set_lower_bound(x[2], 0.0)
    JuMP.set_lower_bound(y, -20.0)
    JuMP.set_upper_bound(x[1], 3.0)
    JuMP.set_upper_bound(x[2], 4.0)
    JuMP.set_upper_bound(y, 20.0)

    @constraint(m, e1, x[1] + x[2] + y == 0.0)
    @NLconstraint(m, e2, x[2] - 8*(x[1])^2 + 8*(x[1])^3-2*(x[1])^4 <= 2.0)
    @NLconstraint(m, e3, 32*(x[1])^3-4*(x[1])^4-88*(x[1])^2+96*x[1]+x[2] <= 36.0)

    @objective(m, Min, y)
    JuMP.optimize!(m)

    @test isapprox(JuMP.objective_value(m), -5.5080, atol=1E-3)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
end

#include("Optimizer/quadratic_relaxation.jl")
#include("implicit_optimizer.jl")
=#
