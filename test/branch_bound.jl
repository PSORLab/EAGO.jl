@testset "Test Continuous Branch Rules" begin
    B = EAGO.Optimizer()
    B.variable_number = 2
    B.fixed_variable[1]  = false
    B.fixed_variable[2]  = false
    B.bisection_variable[1] = true
    B.bisection_variable[2] = true
    B.variable_info = [EAGO.VariableInfo(false,1.0,false,2.0,false,false),
                       EAGO.VariableInfo(false,2.0,false,6.0,false,false)]
    B.current_lower_info.solution = [1.4, 5.3]
    S = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -Inf, Inf, 2, 1, true)
    X1,X2 = EAGO.continuous_relative_bisect(B,S)
    @test isapprox(X1.lower_variable_bounds[1], 1.0; atol = 1E-4)
    @test isapprox(X1.upper_variable_bounds[1], 1.475; atol = 1E-2)
    @test isapprox(X2.lower_variable_bounds[1], 1.475; atol = 1E-2)
    @test isapprox(X2.upper_variable_bounds[1], 2.0; atol = 1E-4)
end

@testset "Test Implicit Branch Rules" begin
    B = EAGO.Optimizer()
    B.variable_number = 2
    B.fixed_variable[1]  = false
    B.fixed_variable[2]  = false
    B.fixed_variable[3]  = false
    B.bisection_variable[1] = true
    B.bisection_variable[2] = true
    B.bisection_variable[3] = true
    B.variable_info = [EAGO.VariableInfo(false,1.0,false,2.0,false,false),
                       EAGO.VariableInfo(false,2.0,false,6.0,false,false),
                       EAGO.VariableInfo(false,2.0,false,6.0,false,false)]
    S = EAGO.NodeBB(Float64[1.0,2.0,2.0], Float64[1.5,5.0,5.5], -Inf, Inf, 2, 1, true)
    B.current_lower_info.solution = [1.25, 3.5, 4.0]

    B.working_evaluator_block = JuMP._NLPData()
    B.working_evaluator_block.evaluator = EAGO.ImplicitLowerEvaluator{2}()
    X1,X2 = EAGO.implicit_bisection(B,S)

    @test isapprox(X1.lower_variable_bounds[2], 2.0; atol = 1E-4)
    @test isapprox(X1.upper_variable_bounds[2], 3.48829; atol = 1E-2)
    @test isapprox(X2.lower_variable_bounds[2], 3.48829; atol = 1E-2)
    @test isapprox(X2.upper_variable_bounds[2], 5.0; atol = 1E-4)
end

@testset "Test B&B Checks" begin
    B = EAGO.Optimizer()
    B.variable_number = 2
    B.variable_info = [EAGO.VariableInfo(false,1.0,false,2.0,false,false),
                      EAGO.VariableInfo(false,2.0,false,6.0,false,false)]
    S = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -Inf, Inf, 2, 1, true)

    @test EAGO.default_repeat_check(B,S) == false
    @test EAGO.default_termination_check(B) == false

    B.stack[1] = S; B.iteration_limit = -1; B.current_iteration_count = 2
    @test EAGO.default_termination_check(B) == false

    B.iteration_limit = 1E8; B.node_limit = -1;
    @test EAGO.default_termination_check(B) == false

    B.node_limit = 1E8; B.current_lower_info.value = 1.1;
    B.current_upper_info.value = 1.1 + 1.0E-6; B.absolute_tolerance = 1.0E-4
    @test EAGO.default_termination_check(B) == true

    B.absolute_tolerance = 1.0E-1; B.relative_tolerance = 1.0E-12
    @test EAGO.default_termination_check(B) == true

    B.current_lower_info.value = -Inf;  B.current_upper_info.value = Inf
    B.relative_tolerance = 1.0E10;
    @test EAGO.default_termination_check(B) == true

    B.current_lower_info.value = 2.1;  B.current_upper_info.value = 2.1+1E-9
    B.relative_tolerance = 1.0E10; B.absolute_tolerance = 1.0E-6
    @test EAGO.default_termination_check(B) == true

    B.relative_tolerance = 1.0E-6; B.absolute_tolerance = 1.0E10
    @test EAGO.default_termination_check(B) == true

    B.relative_tolerance = 1.0E-20; B.absolute_tolerance = 1.0E-20
    @test EAGO.default_termination_check(B) == true
end

@testset "Find Lower Bound" begin
    B = EAGO.Optimizer()
    B.global_upper_bound = -4.5
    B.stack[1] = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1, true)
    B.stack[2] = EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1, true)
    B.stack[3] = EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1, true)
    Lower = EAGO.find_lower_bound(B)

    @test Lower == -5.0
end

@testset "Test Fathom!" begin
    B = EAGO.Optimizer()
    B.global_upper_bound = -4.5
    B.stack[1] = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1, true)
    B.stack[2] = EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1, true)
    B.stack[3] = EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1, true)
    EAGO.fathom!(B)

    @test length(B.stack) == 1
    @test B.stack[2].lower_bound == -5.0
end

@testset "Node Selection" begin
    B = EAGO.Optimizer()
    B.global_upper_bound = -4.5
    B.stack[1] = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1, true)
    B.stack[2] = EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1, true)
    B.stack[3] = EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1, true)
    key,node = EAGO.node_select_best!(B)

    @test key == 2
    @test node.lower_bound == -5.0
end

@testset "Node Storage" begin
    B = EAGO.Optimizer()
    y = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1, true)
    y1 = EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1, true)
    y2 = EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1, true)
    EAGO.single_storage!(B,y)
    @test B.maximum_node_id == 1

    EAGO.default_storage!(B,y1,y2)
    @test B.stack[1].lower_bound == -4.0
    @test B.stack[2].lower_bound == -5.0
    @test B.stack[3].lower_bound == -2.0
    @test B.maximum_node_id == 3
end

@testset "Node Access Functions" begin
    x = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -3.4, 2.1, 2, 1, true)

    @test EAGO.lower_variable_bounds(x) == Float64[1.0,5.0]
    @test EAGO.upper_variable_bounds(x) == Float64[2.0,6.0]
    @test EAGO.lower_bound(x) == -3.4
    @test EAGO.upper_bound(x) == 2.1
    @test EAGO.depth(x) == 2
    @test EAGO.last_branch(x) == 1
    @test EAGO.branch_direction(x) == true
end
