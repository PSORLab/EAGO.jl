@testset "Test Continuous Branch Rules" begin
    B = EAGO.Optimizer(verbosity = 0)
    B._variable_number = 2
    B._fixed_variable = [false, false]
    B.branch_variable = [true, true]
    B._variable_info = [EAGO.VariableInfo(false,1.0,false,2.0,false,false),
                       EAGO.VariableInfo(false,2.0,false,6.0,false,false)]
    B._lower_solution = [1.4, 5.3]
    S = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -Inf, Inf, 2, 1)
    EAGO.branch_node!(B)
    #@test isapprox(X1.lower_variable_bounds[1], 1.0; atol = 1E-4)
    #@test isapprox(X1.upper_variable_bounds[1], 1.475; atol = 1E-2)
    #@test isapprox(X2.lower_variable_bounds[1], 1.475; atol = 1E-2)
    #@test isapprox(X2.upper_variable_bounds[1], 2.0; atol = 1E-4)
end

@testset "Test Implicit Branch Rules" begin
    B = EAGO.Optimizer(verbosity = 0)
    B._variable_number = 2
    B._fixed_variable  = [false, false, false]
    B.branch_variable = [true, true, true]
    B._variable_info = [EAGO.VariableInfo(false,1.0,false,2.0,false,false),
                       EAGO.VariableInfo(false,2.0,false,6.0,false,false),
                       EAGO.VariableInfo(false,2.0,false,6.0,false,false)]
    S = EAGO.NodeBB(Float64[1.0,2.0,2.0], Float64[1.5,5.0,5.5], -Inf, Inf, 2, 1)
    B._lower_solution = [1.25, 3.5, 4.0]
end

@testset "Test B&B Checks" begin
    B = EAGO.Optimizer(verbosity = 0)
    B._variable_number = 2
    B._variable_info = [EAGO.VariableInfo(false,1.0,false,2.0,false,false),
                      EAGO.VariableInfo(false,2.0,false,6.0,false,false)]
    S = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -Inf, Inf, 2, 1)

    @test EAGO.repeat_check(B) == false

    #=
    @test EAGO.termination_check(B) == false

    push!(B._stack, S); B.iteration_limit = -1; B._iteration_count = 2
    @test EAGO.termination_check(B) == false

    B.iteration_limit = 1E8; B.node_limit = -1;
    @test EAGO.termination_check(B) == false

    B.node_limit = 1E8; B._lower_objective_value = 1.1;
    B._upper_objective_value = 1.1 + 1.0E-6; B.absolute_tolerance = 1.0E-4
    @test EAGO.termination_check(B) == true

    B.absolute_tolerance = 1.0E-1; B.relative_tolerance = 1.0E-12
    @test EAGO.termination_check(B) == true

    B._lower_objective_value = -Inf;  B._upper_objective_value = Inf
    B.relative_tolerance = 1.0E10;
    @test EAGO.termination_check(B) == true

    B._lower_objective_value = 2.1;  B._upper_objective_valuee = 2.1+1E-9
    B.relative_tolerance = 1.0E10; B.absolute_tolerance = 1.0E-6
    @test EAGO.termination_check(B) == true

    B.relative_tolerance = 1.0E-6; B.absolute_tolerance = 1.0E10
    @test EAGO.termination_check(B) == true

    B.relative_tolerance = 1.0E-20; B.absolute_tolerance = 1.0E-20
    @test EAGO.termination_check(B) == true
    =#
end

@testset "Find Lower Bound" begin
    B = EAGO.Optimizer(verbosity = 0)
    B._global_upper_bound = -4.5
    push!(B._stack, EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1))
    push!(B._stack, EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1))
    push!(B._stack, EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1))
    EAGO.set_global_lower_bound!(B)

    @test B._global_lower_bound == -5.0
end

@testset "Test Fathom!" begin
    B = EAGO.Optimizer(verbosity = 0)
    B._global_upper_bound = -4.5
    push!(B._stack, EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1))
    push!(B._stack, EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1))
    push!(B._stack, EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1))
    EAGO.fathom!(B)

    @test length(B._stack) == 1
end

@testset "Node Selection" begin
    B = EAGO.Optimizer(verbosity = 0)
    B._global_upper_bound = -4.5
    push!(B._stack, EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1))
    push!(B._stack, EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1))
    push!(B._stack, EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1))
    EAGO.node_selection!(B)

    @test B._current_node.lower_bound == -5.0
end

@testset "Node Storage" begin
    B = EAGO.Optimizer(verbosity = 0)
    y = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1)
    y1 = EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1)
    y2 = EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1)
    EAGO.single_storage!(B)
    @test B.maximum_node_id == 1
end

@testset "Node Access Functions" begin
    x = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -3.4, 2.1, 2, 1)

    @test EAGO.lower_variable_bounds(x) == Float64[1.0,5.0]
    @test EAGO.upper_variable_bounds(x) == Float64[2.0,6.0]
    @test EAGO.lower_bound(x) == -3.4
    @test EAGO.upper_bound(x) == 2.1
    @test EAGO.depth(x) == 2
    @test EAGO.last_branch(x) == 1
    @test EAGO.branch_direction(x) == true
end
