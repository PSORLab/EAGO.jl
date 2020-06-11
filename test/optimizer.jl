@testset "Set/Get Attributes" begin

    m = EAGO.Optimizer()
    @test MOI.get(m, MOI.SolverName()) === "EAGO: Easy Advanced Global Optimization"

    m._maximum_node_id = 55
    @test MOI.get(m, MOI.NodeCount()) === 55

    m._result_status_code = MOI.FEASIBLE_POINT
    @test MOI.get(m, MOI.ResultCount()) === 1

    m._result_status_code = MOI.OTHER_RESULT_STATUS
    @test MOI.get(m, MOI.ResultCount()) === 0

    m._global_lower_bound = 4.0
    m._global_upper_bound = 6.0
    m._input_problem._optimization_sense = MOI.MIN_SENSE
    @test isapprox(MOI.get(m, MOI.RelativeGap()), 0.33333333, atol=1E-5)

    m._input_problem._optimization_sense = MOI.MAX_SENSE
    @test MOI.get(m, MOI.RelativeGap()) === 0.5
    @test MOI.get(m, MOI.ObjectiveBound()) === -4.0

    m._parameters.verbosity = 2
    m._parameters.log_on = true
    MOI.set(m, MOI.Silent(), 1)
    @test m._parameters.verbosity === 0
    @test m._parameters.log_on === false

    @test MOI.supports(m, MOI.ObjectiveSense())
end

@testset "Get Termination Code " begin

    model = EAGO.Optimizer()

    # Termination Status Code Checks
    model._termination_status_code = MOI.OPTIMAL
    status = @inferred MOI.get(model, MOI.TerminationStatus())
    @test status == MOI.OPTIMAL

    model._termination_status_code = MOI.ITERATION_LIMIT
    status = @inferred MOI.get(model, MOI.TerminationStatus())
    @test status == MOI.ITERATION_LIMIT
end

@testset "Add Variable, Get Number/Index" begin

    model = EAGO.Optimizer()

    # Add variable
    @inferred MOI.add_variable(model)
    nVar = @inferred MOI.get(model, MOI.NumberOfVariables())
    @test nVar == 1

    # Add second variable
    MOI.add_variables(model,3)
    nVar = @inferred MOI.get(model, MOI.NumberOfVariables())
    @test nVar == 4

    # Get variable indices
    indx = @inferred MOI.get(model, MOI.ListOfVariableIndices())
    @test indx == @inferred MOI.VariableIndex[MOI.VariableIndex(1), MOI.VariableIndex(2),
                                    MOI.VariableIndex(3), MOI.VariableIndex(4)]

    @test_nowarn @inferred EAGO.check_inbounds!(model, MOI.VariableIndex(1))
    @test_throws ErrorException @inferred EAGO.check_inbounds!(model,MOI.VariableIndex(6))
end

@testset "Add Variable Bounds" begin
    model = EAGO.Optimizer()

    @test MOI.supports_constraint(model, MOI.SingleVariable, MOI.LessThan{Float64})
    @test MOI.supports_constraint(model, MOI.SingleVariable, MOI.GreaterThan{Float64})
    @test MOI.supports_constraint(model, MOI.SingleVariable, MOI.EqualTo{Float64})

    x = MOI.add_variables(model,3)
    z = MOI.add_variable(model)

    @inferred MOI.add_constraint(model, MOI.SingleVariable(x[1]), MOI.GreaterThan(-1.0))
    @inferred MOI.add_constraint(model, MOI.SingleVariable(x[2]), MOI.LessThan(-1.0))
    @inferred MOI.add_constraint(model, MOI.SingleVariable(x[3]), MOI.EqualTo(2.0))

    @test model._input_problem._variable_info[1].is_integer == false
    @test model._input_problem._variable_info[1].lower_bound == -1.0
    @test model._input_problem._variable_info[1].has_lower_bound == true
    @test model._input_problem._variable_info[1].upper_bound == Inf
    @test model._input_problem._variable_info[1].has_upper_bound == false
    @test model._input_problem._variable_info[1].is_fixed == false

    @test model._input_problem._variable_info[2].is_integer == false
    @test model._input_problem._variable_info[2].lower_bound == -Inf
    @test model._input_problem._variable_info[2].has_lower_bound == false
    @test model._input_problem._variable_info[2].upper_bound == -1.0
    @test model._input_problem._variable_info[2].has_upper_bound == true
    @test model._input_problem._variable_info[2].is_fixed == false

    @test model._input_problem._variable_info[3].is_integer == false
    @test model._input_problem._variable_info[3].lower_bound == 2.0
    @test model._input_problem._variable_info[3].has_lower_bound == true
    @test model._input_problem._variable_info[3].upper_bound == 2.0
    @test model._input_problem._variable_info[3].has_upper_bound == true
    @test model._input_problem._variable_info[3].is_fixed == true
end

@testset "Add Linear Constraint " begin

    model = EAGO.Optimizer()

    @test MOI.supports_constraint(model, MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64})
    @test MOI.supports_constraint(model, MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64})
    @test MOI.supports_constraint(model, MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64})

    x = MOI.add_variables(model,3)

    func1 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.(Float64[5.0,-2.3],[x[1],x[2]]),2.0)
    func2 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.(Float64[4.0,-2.2],[x[2],x[3]]),2.1)
    func3 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.(Float64[3.0,-3.3],[x[1],x[3]]),2.2)

    set1 = MOI.LessThan{Float64}(1.0)
    set2 = MOI.GreaterThan{Float64}(2.0)
    set3 = MOI.EqualTo{Float64}(3.0)

    @inferred MOI.add_constraint(model, func1, set1)
    @inferred MOI.add_constraint(model, func2, set2)
    @inferred MOI.add_constraint(model, func3, set3)

    @test model._input_problem._linear_leq_constraints[1][1].constant == 2.0
    @test model._input_problem._linear_geq_constraints[1][1].constant == 2.1
    @test model._input_problem._linear_eq_constraints[1][1].constant == 2.2
    @test model._input_problem._linear_leq_constraints[1][1].terms[1].coefficient == 5.0
    @test model._input_problem._linear_geq_constraints[1][1].terms[1].coefficient == 4.0
    @test model._input_problem._linear_eq_constraints[1][1].terms[1].coefficient == 3.0
    @test model._input_problem._linear_leq_constraints[1][1].terms[2].coefficient == -2.3
    @test model._input_problem._linear_geq_constraints[1][1].terms[2].coefficient == -2.2
    @test model._input_problem._linear_eq_constraints[1][1].terms[2].coefficient == -3.3
    @test model._input_problem._linear_leq_constraints[1][1].terms[1].variable_index.value == 1
    @test model._input_problem._linear_geq_constraints[1][1].terms[1].variable_index.value == 2
    @test model._input_problem._linear_eq_constraints[1][1].terms[1].variable_index.value == 1
    @test model._input_problem._linear_leq_constraints[1][1].terms[2].variable_index.value == 2
    @test model._input_problem._linear_geq_constraints[1][1].terms[2].variable_index.value == 3
    @test model._input_problem._linear_eq_constraints[1][1].terms[2].variable_index.value == 3
    @test MOI.LessThan{Float64}(1.0) == model._input_problem._linear_leq_constraints[1][2]
    @test MOI.GreaterThan{Float64}(2.0) == model._input_problem._linear_geq_constraints[1][2]
    @test MOI.EqualTo{Float64}(3.0) == model._input_problem._linear_eq_constraints[1][2]
end

@testset "Add Quadratic Constraint " begin

    model = EAGO.Optimizer()

    @test MOI.supports_constraint(model, MOI.ScalarQuadraticFunction{Float64}, MOI.LessThan{Float64})
    @test MOI.supports_constraint(model, MOI.ScalarQuadraticFunction{Float64}, MOI.GreaterThan{Float64})
    @test MOI.supports_constraint(model, MOI.ScalarQuadraticFunction{Float64}, MOI.EqualTo{Float64})

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

    @inferred MOI.add_constraint(model, func1, set1)
    @inferred MOI.add_constraint(model, func2, set2)
    @inferred MOI.add_constraint(model, func3, set3)

    @test model._input_problem._quadratic_leq_constraints[1][1].constant == 2.0
    @test model._input_problem._quadratic_geq_constraints[1][1].constant == 2.1
    @test model._input_problem._quadratic_eq_constraints[1][1].constant == 2.2
    @test model._input_problem._quadratic_leq_constraints[1][1].quadratic_terms[1].coefficient == 2.5
    @test model._input_problem._quadratic_geq_constraints[1][1].quadratic_terms[1].coefficient == 2.2
    @test model._input_problem._quadratic_eq_constraints[1][1].quadratic_terms[1].coefficient == 2.1
    @test model._input_problem._quadratic_leq_constraints[1][1].affine_terms[1].coefficient == 5.0
    @test model._input_problem._quadratic_geq_constraints[1][1].affine_terms[1].coefficient == 4.0
    @test model._input_problem._quadratic_eq_constraints[1][1].affine_terms[1].coefficient == 3.0
    @test model._input_problem._quadratic_leq_constraints[1][1].quadratic_terms[1].variable_index_1.value == 2
    @test model._input_problem._quadratic_geq_constraints[1][1].quadratic_terms[1].variable_index_1.value == 1
    @test model._input_problem._quadratic_eq_constraints[1][1].quadratic_terms[1].variable_index_1.value == 1
    @test model._input_problem._quadratic_leq_constraints[1][1].quadratic_terms[1].variable_index_2.value == 2
    @test model._input_problem._quadratic_geq_constraints[1][1].quadratic_terms[1].variable_index_2.value == 2
    @test model._input_problem._quadratic_eq_constraints[1][1].quadratic_terms[1].variable_index_2.value == 1
    @test model._input_problem._quadratic_leq_constraints[1][1].affine_terms[1].variable_index.value == 1
    @test model._input_problem._quadratic_geq_constraints[1][1].affine_terms[1].variable_index.value == 2
    @test model._input_problem._quadratic_eq_constraints[1][1].affine_terms[1].variable_index.value == 3
    @test MOI.LessThan{Float64}(1.0) == model._input_problem._quadratic_leq_constraints[1][2]
    @test MOI.GreaterThan{Float64}(2.0) == model._input_problem._quadratic_geq_constraints[1][2]
    @test MOI.EqualTo{Float64}(3.0) == model._input_problem._quadratic_eq_constraints[1][2]
end

@testset "Set Objective" begin
    model = @inferred EAGO.Optimizer()

    x = MOI.add_variables(model,3)

    @inferred MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test model._input_problem._optimization_sense == MOI.MIN_SENSE

    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    @test model._input_problem._optimization_sense == MOI.MAX_SENSE

    MOI.set(model, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    @test model._input_problem._optimization_sense == MOI.FEASIBILITY_SENSE

    @test MOI.supports(model, MOI.ObjectiveFunction{MOI.SingleVariable}())
    @test MOI.supports(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    @test MOI.supports(model, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}())

    x = MOI.add_variables(model,3)

    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(MOI.VariableIndex(2)))
    @test model._input_problem._objective_type == EAGO.SINGLE_VARIABLE
    @test model._input_problem._objective_sv == MOI.SingleVariable(MOI.VariableIndex(2))

    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.(Float64[5.0,-2.3],[x[1],x[2]]),2.0))
    @test model._input_problem._objective_type == EAGO.SCALAR_AFFINE
    @test model._input_problem._objective_saf.constant == 2.0

    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(),
                                         MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(5.0,x[1])],
                                        [MOI.ScalarQuadraticTerm{Float64}(2.5,x[2],x[2])],3.0))
    @test model._input_problem._objective_type == EAGO.SCALAR_QUADRATIC
    @test model._input_problem._objective_sqf.constant == 3.0
end

@testset "Empty/Isempty, EAGO Model, Single Storage, Optimize Hook " begin
    model = EAGO.Optimizer()
    @test @inferred MOI.is_empty(model)

    t = EAGO.DefaultExt()
    model._current_node = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1)
    model._lower_objective_value = -3.0
    model._upper_objective_value = 0.0
    @test_nowarn EAGO.single_storage!(t, model)
    new_node = pop!(model._stack)
    @test new_node.lower_bound == -3.0
    @test new_node.upper_bound == 0.0

    @test_nowarn EAGO.single_storage!(model)

    @test_nowarn EAGO.optimize_hook!(EAGO.DefaultExt(), model)
    @test_nowarn EAGO.throw_optimize_hook!(model)
end
#=
@testset "Fallback Interval Bounds" begin

    # test linear expression interval fallback
    X = Interval(1.0, 2.0)
    Y = Interval(-4.0, 6.0)
    Z = Interval(-6.0,-5.0)

    ylower = Float64[X.lo; Y.lo; Z.lo]
    yupper = Float64[X.hi; Y.hi; Z.hi]

    n = EAGO.NodeBB(ylower, yupper, -Inf, Inf, 3, 2)

    coeff = Float64[3.0; 4.0; -7.0]
    cons = 3.2
    fintv_saf = coeff[1]*X + coeff[2]*Y + coeff[3]*Z + cons
    x = MOI.VariableIndex.(Int64[1;2;3])
    saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(coeff, x), cons)
    saf_lo = @inferred EAGO.interval_bound(saf, n, true)
    saf_hi = @inferred EAGO.interval_bound(saf, n, false)
    @test isapprox(saf_lo, fintv_saf.lo, atol = 1E-7)
    @test isapprox(saf_hi, fintv_saf.hi, atol = 1E-7)

    # test quadratic expression interval fallback
    vi1 = [MOI.VariableIndex(1); MOI.VariableIndex(2); MOI.VariableIndex(3);
           MOI.VariableIndex(1); MOI.VariableIndex(2)];
    vi2 = [MOI.VariableIndex(1); MOI.VariableIndex(2); MOI.VariableIndex(3);
           MOI.VariableIndex(2); MOI.VariableIndex(3)];
    coef1 = Float64[1.0; 4.0; 1.0; 3.0; -5.5]
    sat_arr = MOI.ScalarAffineTerm{Float64}[MOI.ScalarAffineTerm{Float64}(4.0,MOI.VariableIndex(2))]
    sqf =  MOI.ScalarQuadraticFunction(sat_arr,
                                       MOI.ScalarQuadraticTerm.(coef1, vi1, vi2),
                                       cons)
    fintv_sqf = X^2 + 4.0*Y^2 + Z^2 + 3.0*X*Y - 5.5*Y*Z + 4.0*Y + 3.2
    sqf_lo = @inferred EAGO.interval_bound(sqf, n, true)
    sqf_hi = @inferred EAGO.interval_bound(sqf, n, false)
    @test isapprox(sqf_lo, fintv_sqf.lo, atol = 1E-7)
    @test isapprox(sqf_hi, fintv_sqf.hi, atol = 1E-7)
end

@testset "Fallback Interval Bound 2" begin
    # test linear expression interval fallback
    n = EAGO.NodeBB([-1.0; -1.0], [2.0; 2.0], -Inf, Inf, 3, 2)

    m1 = Model(with_optimizer(EAGO.Optimizer, verbosity = 4, cp_depth = -1))
    @variable(m1, -1.0 <= x <= 2.0)
    @variable(m1, -1.0 <= y <= 2.0)
    @constraint(m1, x^2 + y + x*y <= 50.0)
    @constraint(m1, x^2 + y >= -50.0)
    @constraint(m1, 2x + y <= 50.0)
    @constraint(m1, 2x - y >= -50.0)
    @NLconstraint(m1, x^2 + sin(x) <= 50.0)
    @NLconstraint(m1, x^2 + sin(x) >= -50.0)
    @NLobjective(m1, Min, cos(x)*x)
    optimize!(m1)
    obj_value = objective_value(m1)
    @test isapprox(obj_value, -0.8322104859004729, atol=1E-3)

    b = backend(m1).optimizer.model.optimizer
    EAGO.interval_lower_bound!(b, n)
    #@test isapprox(b._lower_objective_value, -0.1513, atol=1E-3)
end
=#

@testset "LP Problems" begin
    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0,
                                        "presolve_scrubber_flag" => false,
                                        "presolve_to_JuMP_flag" => false))

    @variable(m, 1 <= x <= 3)
    @variable(m, 1 <= y <= 3)

    @objective(m, Min, x + y)

    @constraint(m, x + y <= 10)
    @constraint(m, x - y <= 10)
    @constraint(m, y >= 0)

    JuMP.optimize!(m)

    @test isapprox(JuMP.value(x), 1.0, atol=1E-4)
    @test isapprox(JuMP.value(y), 1.0, atol=1E-4)
    @test isapprox(JuMP.objective_value(m), 2.0, atol=1E-4)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0,
                                        "presolve_scrubber_flag" => false,
                                        "presolve_to_JuMP_flag" => false))

    @variable(m, -3 <= x <= -1)
    @variable(m, -2 <= y <= 2)
    @variable(m, 1 <= z <= 3)

    @objective(m, Min, x - y + 2z)

    @constraint(m, x + 2y >= -10)
    @constraint(m, z - 2y <= 2)
    @constraint(m, y >= 0)

    JuMP.optimize!(m)

    @test isapprox(JuMP.value(x), -3.0, atol=1E-4)
    @test isapprox(JuMP.value(y), 2.0, atol=1E-4)
    @test isapprox(JuMP.value(z), 1.0, atol=1E-4)
    @test isapprox(JuMP.objective_value(m), -3.0, atol=1E-4)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0,
                                        "presolve_scrubber_flag" => false,
                                        "presolve_to_JuMP_flag" => false))

    @variable(m, -3 <= x <= -1)
    @variable(m, -2 <= y <= 2)
    @variable(m, 1 <= z <= 3)
    @variable(m, -10 <= q <= 9)

    @objective(m, Min, 2x - 3y + 2z)

    @constraint(m, x + 2y >= -10)
    @constraint(m, z - 2y <= 2)
    @constraint(m, y >= 0)
    @constraint(m, q-3*z-y >= 0)

    JuMP.optimize!(m)

    @test isapprox(JuMP.value(x), -3.0, atol=1E-4)
    @test isapprox(JuMP.value(y), 2.0, atol=1E-4)
    @test isapprox(JuMP.value(z), 1.0, atol=1E-4)
    @test isapprox(JuMP.objective_value(m), -10.0, atol=1E-4)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0,
                                        "presolve_scrubber_flag" => false,
                                        "presolve_to_JuMP_flag" => false))

    @variable(m, -3 <= x <= -1)
    @variable(m, -2 <= y <= 2)
    @variable(m, 1 <= z <= 3)
    @variable(m, -10 <= q <= 9)

    @NLobjective(m, Min, 2x - 3y + 2z)

    @constraint(m, x + 2y >= -10)
    @constraint(m, z - 2y <= 2)
    @constraint(m, y >= 4)
    @constraint(m, q-3*z-y >= 0)

    JuMP.optimize!(m)

    @test isapprox(JuMP.value(x), 0.0, atol=1E-4)
    @test isapprox(JuMP.value(y), 0.0, atol=1E-4)
    @test isapprox(JuMP.value(z), 0.0, atol=1E-4)
    @test isapprox(JuMP.value(q), 0.0, atol=1E-4)
    @test isapprox(JuMP.objective_value(m), Inf, atol=1E-4)
    @test JuMP.termination_status(m) == MOI.INFEASIBLE
    @test JuMP.primal_status(m) == MOI.INFEASIBILITY_CERTIFICATE

    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0,
                                        "presolve_scrubber_flag" => false,
                                        "presolve_to_JuMP_flag" => false))

    @variable(m, -3 <= x <= -1)
    @variable(m, -2 <= y <= 2)
    @variable(m, 1 <= z <= 3)
    @variable(m, q == 1)

    @objective(m, Max, x)

    @constraint(m, x + 2y >= -10)
    @constraint(m, z - 2y <= 2)
    @constraint(m, y >= 0)

    JuMP.optimize!(m)

    @test isapprox(JuMP.value(x), -1.0, atol=1E-4)
    @test isapprox(JuMP.value(q), 1.0, atol=1E-4)
    @test isapprox(JuMP.objective_value(m), -1.0, atol=1E-4)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0,
                                        "presolve_scrubber_flag" => false,
                                        "presolve_to_JuMP_flag" => false))

    @variable(m, -3 <= x <= -1)
    @variable(m, -2 <= y <= 2)
    @variable(m, 1 <= z <= 3)
    @variable(m, -10 <= q <= 9)

    @objective(m, Max, 2x - 3y + 2z)

    @constraint(m, x + 2y >= -10)
    @constraint(m, z - 2y <= 2)
    @constraint(m, y >= 0)
    @constraint(m, q-3*z-y >= 0)

    JuMP.optimize!(m)

    @test isapprox(JuMP.value(x), -1.0, atol=1E-4)
    @test isapprox(JuMP.value(y), 0.428571, atol=1E-4)
    @test isapprox(JuMP.value(z), 2.857142, atol=1E-4)
    @test isapprox(JuMP.objective_value(m), 2.428571, atol=1E-4)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
end

@testset "NLP Problems" begin
    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0))

    @variable(m, -200 <= x <= -100)
    @variable(m, 200 <= y <= 400)
    @constraint(m, -500 <= x + 2y <= 400)
    @NLobjective(m, Min, x*y)
    JuMP.optimize!(m)

    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.value(x), -200.0, atol=1E-5)
    @test isapprox(JuMP.value(y), 300.0, atol=1E-5)
    @test isapprox(JuMP.objective_value(m), -59999.4011899692, atol=2.0)

    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0))
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


    #=
    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0))
    # ----- Variables ----- #
    x_Idx = Any[2, 3, 4]
    @variable(m, x[x_Idx])
    JuMP.set_lower_bound(x[2], 1.0e-6)
    JuMP.set_upper_bound(x[2], 1.0)
    JuMP.set_lower_bound(x[3], 1.0e-6)
    JuMP.set_upper_bound(x[3], 1.0)
    JuMP.set_lower_bound(x[4], 1.0e-6)
    JuMP.set_upper_bound(x[4], 1.0)

    # ----- Constraints ----- #
    @constraint(m, e2, x[2]+x[3]+x[4] == 1.0)

    # ----- Objective ----- #
    @NLobjective(m, Min, ((15.3261663216011*x[2]+23.2043471859416*x[3]+6.69678129464404*x[4])*log(2.1055*x[2]+3.1878*x[3]+0.92*x[4])+1.04055250396734*x[2]-2.24199441248417*x[3]+3.1618173099828*x[4]+6.4661663216011*x[2]*log(x[2]/(2.1055*x[2]+3.1878*x[3]+0.92*x[4]))+12.2043471859416*x[3]*log(x[3]/(2.1055*x[2]+3.1878*x[3]+0.92*x[4]))+0.696781294644034*x[4]*log(x[4]/(2.1055*x[2]+3.1878*x[3]+0.92*x[4]))+9.86*x[2]*log(x[2]/(1.972*x[2]+2.4*x[3]+1.4*x[4]))+12*x[3]*log(x[3]/(1.972*x[2]+2.4*x[3]+1.4*x[4]))+7*x[4]*log(x[4]/(1.972*x[2]+2.4*x[3]+1.4*x[4]))+(1.972*x[2]+2.4*x[3]+1.4*x[4])*log(1.972*x[2]+2.4*x[3]+1.4*x[4])+1.972*x[2]*log(x[2]/(1.972*x[2]+0.283910843616504*x[3]+3.02002220174195*x[4]))+2.4*x[3]*log(x[3]/(1.45991339466884*x[2]+2.4*x[3]+0.415073537580851*x[4]))+1.4*x[4]*log(x[4]/(0.602183324335333*x[2]+0.115623371371275*x[3]+1.4*x[4]))-17.2981663216011*x[2]*log(x[2])-25.6043471859416*x[3]*log(x[3])-8.09678129464404*x[4]*log(x[4])))

    JuMP.optimize!(m)
    @test isapprox(JuMP.objective_value(m), 0.000, atol=1E-3)
    =#

    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0, "absolute_tolerance" => 1.0E-2))
    # ----- Variables ----- #
    x_Idx = Any[2, 3, 4]
    @variable(m, x[x_Idx])
    JuMP.set_lower_bound(x[2], 1.0e-6)
    JuMP.set_upper_bound(x[2], 1.0)
    JuMP.set_lower_bound(x[3], 1.0e-6)
    JuMP.set_upper_bound(x[3], 1.0)
    JuMP.set_lower_bound(x[4], 1.0e-6)
    JuMP.set_upper_bound(x[4], 1.0)

    # ----- Constraints ----- #
    @constraint(m, e2, x[2]+x[3]+x[4] == 1.0)

    # ----- Objective ----- #
    @NLobjective(m, Max, (15.3261663216011*x[2]+23.2043471859416*x[3]+6.69678129464404*x[4])*log(2.1055*x[2]+3.1878*x[3]+0.92*x[4]) - (1.0-2.0) + (2.0*3.0) + (4.1+6.2) + (4.1/6.2) + (3.1^2))

    JuMP.optimize!(m)
    @test isapprox(objective_value(m), 54.0, atol=1E-0)

    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0, "output_iterations" => 0, "absolute_tolerance" => 1.0E-2))
    # ----- Variables ----- #
    x_Idx = Any[2, 3, 4]
    @variable(m, x[x_Idx])
    JuMP.set_lower_bound(x[2], 1.0e-6)
    JuMP.set_upper_bound(x[2], 1.0)
    JuMP.set_lower_bound(x[3], 1.0e-6)
    JuMP.set_upper_bound(x[3], 1.0)
    JuMP.set_lower_bound(x[4], 1.0e-6)
    JuMP.set_upper_bound(x[4], 1.0)

    @NLexpression(m, nl_expr, x[2]*x[3] + cos(x[4]) + tanh(x[2]))
    @NLconstraint(m, nl_expr >= -3.0)
    @NLobjective(m, Min, nl_expr)
    JuMP.optimize!(m)

    @test isapprox(JuMP.objective_value(m), 0.5404086991071391, atol=1E-6)
    @test JuMP.termination_status(m) === MOI.OPTIMAL
    @test JuMP.primal_status(m) === MOI.FEASIBLE_POINT

    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0, "output_iterations" => 0, "absolute_tolerance" => 1.0E-2))

    # ----- Variables ----- #
    xL = [500.0 1300.0 5000.0 100.0 200.0 200.0 200.0 300.0 6900.0]
    xU = [600.0 1500.0 6000.0 200.0 300.0 300.0 400.0 500.0 7100.0]
    @variable(m, xL[i] <= x[i=1:9] <= xU[i])

    # ----- Constraints ----- #
    @constraint(m, e1, x[1] + x[2]+ x[3] - x[9] == 0.0)
    @constraint(m, e2, 0.0025*x[4]+0.0025*x[6] <= 1.0)
    @constraint(m, e3, -0.0025*x[4]+0.0025*x[5]+0.0025*x[7] <= 1.0)
    @constraint(m, e4, -0.01*x[5]+0.01*x[8] <= 1.0)
    @constraint(m, e5, 100*x[1]-x[1]*x[6]+833.33252*x[4] <= 83333.333)
    @constraint(m, e6, -x[2]*x[4]+x[2]*x[7]+1250*x[4]-1250*x[5] >= 0.0)
    @constraint(m, e7, x[3]*x[5]-x[3]*x[8]-2500*x[5] <= -1.25e6)

    # ----- Objective ----- #
    @objective(m, Min, x[9])
    optimize!(m)
    @test isapprox(objective_value(m), 7049.31835811113, atol=1E-2)

    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0, "output_iterations" => 0, "absolute_tolerance" => 1.0E-2))

    xL = [-2.0 0.0]; xU = [2.0 4.0]
    @variable(m, xL[i] <= x[i=1:2] <= xU[i])
    @constraint(m, -x[1]^2 + x[2]^2 - 1.0 == 0.0)
    @constraint(m, 0.4*(x[1]-3.0)^2 + 0.2*x[2]^2 <= 3.0)
    @objective(m, Min, x[2])
    optimize!(m)
    @test isapprox(JuMP.objective_value(m), 1.0652212400578724, atol=1E-3)

    m = Model(with_optimizer(EAGO.Optimizer))
    xL = [-2.0 0.0]; xU = [2.0 4.0]
    @variable(m, xL[i] <= x[i=1:2] <= xU[i])
    @NLobjective(m, Max, x[2]^2 + x[1]^2 + x[1]*x[2])
    optimize!(m)
    @test isapprox(JuMP.objective_value(m), 27.99971055498271, atol=1E-3)
end

@testset "Empty Evaluator" begin
    x = EAGO.EmptyNLPEvaluator()
    n = EAGO.NodeBB()
    @test_nowarn EAGO.set_current_node!(x,n)

    fa = MOI.features_available(x)
    @test fa[1] === :Grad
    @test fa[2] === :Jac
    @test fa[3] === :Hess

    MOI.initialize(x, [:Grad, :Jac, :Hess]) === nothing
    @test MOI.eval_objective(x, 0.0) === NaN

    @test_throws AssertionError MOI.eval_constraint(x, [0.0], 0.0)
    @test_throws AssertionError MOI.eval_constraint_jacobian(x, [0.0], 0.0)
    MOI.eval_objective_gradient(x, [0.0], 0.0) === nothing
    MOI.jacobian_structure(x) === nothing
    MOI.hessian_lagrangian_structure(x) === nothing
    @test_throws AssertionError MOI.eval_hessian_lagrangian(x, [0.0], 0.0, 0.0, 0.0)
    MOI.eval_hessian_lagrangian(x, [], 0.0, 0.0, 0.0) === nothing
end
#=
@testset "User Defined Function Scrubber" begin
    gamma1_x1(z) = z[1]*(1253/z[3])/(1 + 2.62*(z[1]/z[2]))^2
    gamma2_x2(z) = z[2]*(479/z[3])/(1 + 0.382*(z[2]/z[1]))^2

    function cons_1ex(z...)
        temp = ForwardDiff.gradient(z -> log(gamma1_x1(z)), collect(z))
        temp[1]
    end

    unity(x) = x::Float64
    function cons_2ex(z...)
        temp = ForwardDiff.gradient(z -> log(gamma2_x2(z)), collect(z))
        unity(temp[2])
    end

    # Define the JuMP model and solve
    m = Model(with_optimizer(EAGO.Optimizer, presolve_scrubber_flag = true))
    register(m, :cons_1ex, 3, cons_1ex, autodiff = true)
    register(m, :cons_2ex, 3, cons_2ex, autodiff = true)
    @variable(m, 0.01 <= x1 <= 0.99)
    @variable(m, 0.01 <= x2 <= 0.99)
    @variable(m, 363.15 <= T <= 398.15)
    @constraint(m, x1 + x2 == 1.0)
    @NLexpression(m, P1, 1.33*exp(11.83572 - 4169.84/(T - 17.665)))
    @NLexpression(m, P2, 1.33*exp(11.33986 - 3724.523/(T - 69.854)))

    @NLconstraint(m, cons_1ex(x1,x2,T)*P1 + cons_2ex(x1,x2,T)*P2 == 1.02)
    @NLconstraint(m, cons_1ex(x1,x2,T)*P1/1.02 >= 0.95)
    @NLconstraint(m, cons1, cons_1ex(x1, x2, T) >= 0.001)
    @NLconstraint(m, cons2, cons_2ex(x1, x2, T) >= 0.001)
    @NLobjective(m, Min, T)
    optimize!(m)
    @test MOI.INFEASIBILITY_CERTIFICATE === primal_status(m)
end
=#
