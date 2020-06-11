@testset "Test Stack Management Functions" begin
    x = EAGO.Optimizer()
    x._parameters.verbosity = 0
    x._branch_variable_count = 2
    x._branch_to_sol_map = [1; 2]
    x._fixed_variable = [false, false]
    x._parameters.branch_variable = [true, true]
    x._working_problem._variable_info = [EAGO.VariableInfo(false,1.0,false,2.0,false,false,EAGO.BRANCH),
                       EAGO.VariableInfo(false,2.0,false,6.0,false,false,EAGO.BRANCH)]
    x._lower_solution = [1.4, 5.3]
    y = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -Inf, Inf, 2, 1)
    x._current_node = y
    EAGO.branch_node!(x)

    EAGO.node_selection!(x.ext_type, x)
    @test isapprox(x._current_node.lower_variable_bounds[1], 1.0; atol = 1E-4)
    @test isapprox(x._current_node.upper_variable_bounds[1], 1.475; atol = 1E-2)

    EAGO.node_selection!(x.ext_type, x)
    @test isapprox(x._current_node.lower_variable_bounds[1], 1.475; atol = 1E-2)
    @test isapprox(x._current_node.upper_variable_bounds[1], 2.0; atol = 1E-4)

    x._global_upper_bound = -4.5
    empty!(x._stack)
    push!(x._stack, EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1))
    push!(x._stack, EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1))
    push!(x._stack, EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1))
    EAGO.node_selection!(x)
    @test x._current_node.lower_bound == -5.0

    x._global_upper_bound = -4.5
    push!(x._stack, EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1))
    push!(x._stack, EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1))
    push!(x._stack, EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1))
    EAGO.set_global_lower_bound!(x)
    @test x._global_lower_bound == -5.0

    empty!(x._stack)
    x._global_upper_bound = -4.5
    push!(x._stack, EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1))
    push!(x._stack, EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1))
    push!(x._stack, EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1))
    EAGO.fathom!(x)

    @test length(x._stack) == 1

    @inferred EAGO.branch_node!(x)
    @inferred EAGO.node_selection!(x.ext_type, x)
    @inferred EAGO.set_global_lower_bound!(x)
    @inferred EAGO.fathom!(x)

    @test ~isless(EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1), EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1))
    @test length(EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1)) == 2
    diamx = EAGO.diam(EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,7.0], -4.0, 1.0, 2, 1))
    @test diamx[1] == 1.0
    @test diamx[2] == 2.0
    #create_initial_node!(x)
end

@testset "Test B&B Checks" begin
    x = EAGO.Optimizer()
    x._parameters.verbosity = 0
    x._branch_variable_count = 2
    x._working_problem._variable_info = [EAGO.VariableInfo(false,1.0,false,2.0,false,false,EAGO.BRANCH),
                      EAGO.VariableInfo(false,2.0,false,6.0,false,false,EAGO.BRANCH)]
    y = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -Inf, Inf, 2, 1)

    @test EAGO.repeat_check(EAGO.DefaultExt(), x) == false

    @test EAGO.termination_check(EAGO.DefaultExt(), x) == true

    push!(x._stack, y); x._parameters.iteration_limit = -1; x._iteration_count = 2
    @test EAGO.termination_check(EAGO.DefaultExt(), x) == true

    x._parameters.iteration_limit = 1E8; x._parameters.node_limit = -1;
    @test EAGO.termination_check(EAGO.DefaultExt(), x) == true

    x._parameters.node_limit = 1E8; x._lower_objective_value = 1.1;
    x._upper_objective_value = 1.1 + 1.0E-6; x._parameters.absolute_tolerance = 1.0E-4
    @test EAGO.termination_check(EAGO.DefaultExt(), x) == false

    x._parameters.absolute_tolerance = 1.0E-1; x._parameters.relative_tolerance = 1.0E-12
    @test EAGO.termination_check(EAGO.DefaultExt(), x) == false

    x._lower_objective_value = -Inf;  x._upper_objective_value = Inf
    x._parameters.relative_tolerance = 1.0E10;
    @test EAGO.termination_check(EAGO.DefaultExt(), x) == false

    x._lower_objective_value = 2.1;  x._upper_objective_value = 2.1+1E-9
    x._parameters.relative_tolerance = 1.0E10; x._parameters.absolute_tolerance = 1.0E-6
    @test EAGO.termination_check(EAGO.DefaultExt(), x) == false

    x._parameters.relative_tolerance = 1.0E-30; x._parameters.absolute_tolerance = 1.0E10
    @test EAGO.termination_check(EAGO.DefaultExt(), x) == false

    x._parameters.relative_tolerance = 1.0E-20; x._parameters.absolute_tolerance = 1.0E-20
    @test EAGO.termination_check(EAGO.DefaultExt(),x) == false

    @inferred EAGO.repeat_check(x)
    @inferred EAGO.termination_check(x)

    x._run_time = 100.0
    x._parameters.time_limit = 10.0
    @test EAGO.termination_check(EAGO.DefaultExt(),x)
    @test x._termination_status_code === MOI.TIME_LIMIT

    @test @inferred EAGO.is_feasible_solution(MOI.OPTIMAL, MOI.FEASIBLE_POINT)
    @test EAGO.is_feasible_solution(MOI.LOCALLY_SOLVED, MOI.FEASIBLE_POINT)
    @test EAGO.is_feasible_solution(MOI.ALMOST_LOCALLY_SOLVED, MOI.NEARLY_FEASIBLE_POINT)
    @test ~EAGO.is_feasible_solution(MOI.INFEASIBLE, MOI.FEASIBLE_POINT)

    valid, feas = @inferred EAGO.is_globally_optimal(MOI.INFEASIBLE, MOI.INFEASIBILITY_CERTIFICATE)
    @test (valid & ~feas)
    valid, feas = EAGO.is_globally_optimal(MOI.INFEASIBLE, MOI.INFEASIBILITY_CERTIFICATE)
    @test (valid & ~feas)
    valid, feas = EAGO.is_globally_optimal(MOI.INFEASIBLE, MOI.NO_SOLUTION)
    @test (valid & ~feas)
    valid, feas = EAGO.is_globally_optimal(MOI.INFEASIBLE, MOI.UNKNOWN_RESULT_STATUS)
    @test (valid & ~feas)
    valid, feas = EAGO.is_globally_optimal(MOI.OPTIMAL, MOI.FEASIBLE_POINT)
    @test (valid & feas)
    valid, feas = EAGO.is_globally_optimal(MOI.INFEASIBLE_OR_UNBOUNDED, MOI.NO_SOLUTION)
    @test (valid & ~feas)
    valid, feas = EAGO.is_globally_optimal(MOI.SLOW_PROGRESS, MOI.REDUCTION_CERTIFICATE)
    @test ~valid
end

@testset "Node Access Functions" begin
    x = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -3.4, 2.1, 2, 1)
    emptyx = EAGO.NodeBB(Float64[7.0,5.0], Float64[1.0,6.0], -3.4, 2.1, 2, 1)

    @test isempty(emptyx)

    @test EAGO.lower_variable_bounds(x) == Float64[1.0,5.0]
    @test EAGO.upper_variable_bounds(x) == Float64[2.0,6.0]
    @test EAGO.lower_bound(x) == -3.4
    @test EAGO.upper_bound(x) == 2.1
    @test EAGO.depth(x) == 2

    @test_nowarn @inferred EAGO.lower_variable_bounds(x)
    @test_nowarn @inferred EAGO.upper_variable_bounds(x)
    @test_nowarn @inferred EAGO.lower_bound(x)
    @test_nowarn @inferred EAGO.upper_bound(x)
    @test_nowarn @inferred EAGO.depth(x)

    x1 = EAGO.NodeBB(Float64[1.0,5.0], Float64[3.0,6.0], -3.4, 2.1, 2, 1)
    x2 = EAGO.NodeBB(Float64[0.9,5.0], Float64[2.0,6.0], -3.4, 2.1, 2, 1)
    x3 = copy(x)
    x4 = EAGO.NodeBB(x)
    @test ~EAGO.same_box(x, x1, 0.0)
    @test ~EAGO.same_box(x, x2, 0.0)
    @test EAGO.same_box(x, x, 0.0)
    @test EAGO.same_box(x, x3, 0.0)
    @test EAGO.same_box(x, x4, 0.0)

    lower = EAGO.lower_variable_bounds(x, 1, 2)
    upper = EAGO.upper_variable_bounds(x, 1, 2)
    @test lower == Float64[1.0,5.0]
    @test upper == Float64[2.0,6.0]
end
