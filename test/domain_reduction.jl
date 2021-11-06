@testset "Duality-Based Bound Tightening" begin

    ylower = Float64[1.0, 1.0, 1.0, 1.0]
    yupper = Float64[4.0, 4.0, 4.0, 4.0]
    ymult_lo = Float64[50.0, 0.0, 1.0, 0.0]
    ymult_hi = Float64[0.0, 0.0, 0.8, 3.0]
    isint = Bool[false, false]
    n = EAGO.NodeBB(ylower, yupper, isint, true, -Inf, Inf, 2, 1, 1, EAGO.BD_NONE, 1, 0.1)
    @inferred EAGO.variable_dbbt!(n, ymult_lo, ymult_hi, 1.0, 3.0, 4)
    lvb = n.lower_variable_bounds
    uvb = n.upper_variable_bounds

    @test lvb[1] == 1.0; @test uvb[1] == 1.04
    @test lvb[2] == 1.0; @test uvb[2] == 4.0
    @test lvb[3] == 1.0; @test uvb[3] == 3.0
    # lvb[3] doesn't tighten since dbbt assumes variables aren't fixed...
    @test isapprox(lvb[4], 3.33333333, atol= 1E-5); @test uvb[4] == 4.0
end

#=
@testset "Classify Quadratic Types" begin

    @test EAGO.get_value(MOI.LessThan{Float64}(1.1)) == 1.1
    @test EAGO.get_value(MOI.GreaterThan{Float64}(2.1)) == 2.1
    @test EAGO.get_value(MOI.EqualTo{Float64}(1.3)) == 1.3

    opt1 = EAGO.Optimizer()
    opt1._parameters.verbosity = 0
    x = MOI.add_variables(opt1, 3)

    MOI.add_constraint(opt1,MOI.SingleVariable(x[1]), MOI.LessThan(4.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[1]), MOI.GreaterThan(-2.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[2]), MOI.LessThan(4.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[2]), MOI.GreaterThan(0.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[3]), MOI.LessThan(1.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[3]), MOI.GreaterThan(-1.0))

    # adds four univariate quadratic constraints
    qterm1a = [MOI.ScalarQuadraticTerm{Float64}(1.5, MOI.VariableIndex(1), MOI.VariableIndex(1))]
    qterm2a = [MOI.ScalarQuadraticTerm{Float64}(1.5, MOI.VariableIndex(1), MOI.VariableIndex(1))]
    qterm3a = [MOI.ScalarQuadraticTerm{Float64}(1.5, MOI.VariableIndex(2), MOI.VariableIndex(2))]
    qterm4a = [MOI.ScalarQuadraticTerm{Float64}(1.5, MOI.VariableIndex(2), MOI.VariableIndex(2))]
    MOI.add_constraint(opt1,MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm(1.1, x[1])], qterm1a, 2.0), MOI.GreaterThan(2.0))
    MOI.add_constraint(opt1,MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm(1.2, x[1])], qterm2a, 1.0), MOI.LessThan(0.0))
    MOI.add_constraint(opt1,MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm(1.3, x[2])], qterm3a, -2.0), MOI.GreaterThan(2.0))
    MOI.add_constraint(opt1,MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm(1.4, x[2])], qterm4a, -1.0), MOI.EqualTo(0.0))

    # adds four bivariate quadratic constraints
    qterm1b = [MOI.ScalarQuadraticTerm{Float64}(1.5, MOI.VariableIndex(1), MOI.VariableIndex(2))]
    qterm2b = [MOI.ScalarQuadraticTerm{Float64}(1.5, MOI.VariableIndex(2), MOI.VariableIndex(2))]
    qterm3b = [MOI.ScalarQuadraticTerm{Float64}(1.5, MOI.VariableIndex(1), MOI.VariableIndex(1))]
    qterm4b = [MOI.ScalarQuadraticTerm{Float64}(1.5, MOI.VariableIndex(1), MOI.VariableIndex(2))]
    MOI.add_constraint(opt1,MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm(2.1, x[1])], qterm1b, -2.0), MOI.GreaterThan(2.0))
    MOI.add_constraint(opt1,MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm(2.2, x[1])], qterm2b, -1.0), MOI.LessThan(0.0))
    MOI.add_constraint(opt1,MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm(2.3, x[2])], qterm3b, -2.0), MOI.GreaterThan(2.0))
    MOI.add_constraint(opt1,MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm(2.4, x[2])], qterm4b, -1.0), MOI.EqualTo(0.0))
end
=#

#=
@testset "Quadratic Domain Reduction (Univariate)" begin
    m = EAGO.Optimizer(verbosity = 0)

    m._current_node = NodeBB([-10.0, -10.0, -10.0], [10.0, 10.0, 10.0], -Inf, Inf, 1, 1)
    a1, b1, c1 = 2.0, -4.0, 3.0
    a2, b2, c2 = 2.0, -4.0, 3.0
    a3, b3, c3 = 2.0, -4.0, 3.0
    push!(m._univariate_quadratic_geq_constraints,(a1,b1,c1,1))
    push!(m._univariate_quadratic_leq_constraints,(a2,b2,c2,2))
    push!(m._univariate_quadratic_eq_constraints,(a3,b3,c3,3))
    feas = EAGO.univariate_quadratic(m)

    @test feas == true
    @test isapprox(m._current_node.lower_variable_bounds[1],-10.0,atol=1E-3)
    @test isapprox(m._current_node.lower_variable_bounds[2],-0.58113,atol=1E-3)
    @test isapprox(m._current_node.lower_variable_bounds[3],-0.58113,atol=1E-3)
    @test isapprox(m._current_node.upper_variable_bounds[1],10.0,atol=1E-3)
    @test isapprox(m._current_node.upper_variable_bounds[2],2.5811,atol=1E-3)
    @test isapprox(m._current_node.upper_variable_bounds[3],2.5811,atol=1E-3)
end
=#

#=
@testset "Quadratic Domain Reduction (Bivariate)" begin
end

@testset "Interval CSP" begin
    # Vigerske 2017 example
    opt1 = EAGO.Optimizer()
    x = MOI.add_variable(opt1)
    y = MOI.add_variable(opt1)
    z = MOI.add_variable(opt1)

    MOI.add_constraint(opt1,MOI.SingleVariable(x), MOI.LessThan(Inf))
    MOI.add_constraint(opt1,MOI.SingleVariable(x), MOI.GreaterThan(-Inf))
    MOI.add_constraint(opt1,MOI.SingleVariable(y), MOI.LessThan(Inf))
    MOI.add_constraint(opt1,MOI.SingleVariable(y), MOI.GreaterThan(-Inf))
    MOI.add_constraint(opt1,MOI.SingleVariable(z), MOI.LessThan(Inf))
    MOI.add_constraint(opt1,MOI.SingleVariable(z), MOI.GreaterThan(-Inf))

    MOI.add_constraint(opt1,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, x[1]), MOI.ScalarAffineTerm(1.0, x[2])], 0.0), MOI.GreaterThan(4.0))
    MOI.add_constraint(opt1,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, x[1]), MOI.ScalarAffineTerm(1.0, x[2])], 0.0), MOI.LessThan(1.0))

    n1 = NodeData()
    n1.lower_variable_bounds = [-Inf, -Inf, 0.0]
    n1.upper_variable_bounds = [Inf, Inf, 1.0]
    n1.lower_bound = -Inf
    n1.upper_bound = Inf
    n1.depth = 3

    # C1 = exp(x^2 + y^2) - z <= 0
    # C2 = z*sqrt(x^2+y^2) <= 1

    feas1 = EAGO.IntervalCSP(opt1,n1)
    @test (n1.lower_variable_bounds[1] == 0.0)
    @test (n1.lower_variable_bounds[2] == 0.0)
    @test (n1.lower_variable_bounds[3] == 1.0)
    @test (n1.upper_variable_bounds[1] == 0.0)
    @test (n1.upper_variable_bounds[2] == 0.0)
    @test (n1.upper_variable_bounds[3] == 1.0)
    @test feas1 == true
end
=#
#=
@testset "Optimization-Based Bound Tightening (Linear)" begin
    m = Model(optimizer_with_attributes(EAGO.Optimizer, "verbosity" => 0))
    xL = [-4.0; -2.0]
    xU = [4.0; 2.0]

    @variable(m, xL[i] <= x[i=1:2] <= xU[i])

    @objective(m, Min, x[2] - x[1])
    @constraint(m, o1, x[1] >= -1.0)
    @constraint(m, o2, x[2] <= 1.0)
    @constraint(m, e1, -0.1*x[1] + x[2] >= -1.6)
    @constraint(m, e2, -3.7*x[1] + x[2] >= -12.8)
    @constraint(m, e3, 0.1*x[1] - x[2] >= -1.6)
    @constraint(m, e4, 3.7*x[1] - x[2] >= -12.8)

    m.moi_backend.optimizer.model.optimizer._global_lower_bound = Inf
    m.moi_backend.optimizer.model.optimizer._global_upper_bound = -Inf

    optimize!(m)
    opt = m.moi_backend.optimizer.model.optimizer
    m.moi_backend.optimizer.model.optimizer._branch_variables[1] = true
    m.moi_backend.optimizer.model.optimizer._branch_variables[2] = true
    m.moi_backend.optimizer.model.optimizer._global_lower_bound = -12.0
    m.moi_backend.optimizer.model.optimizer._global_upper_bound = 0.0
    opt._current_node = NodeBB(xL, xU,-12.0,0.0,1,1)
    opt.obbt_variable_values = [true; true]
    feas = EAGO.obbt!(opt)

    @test feas
    @test isapprox(opt._current_node.lower_variable_bounds[1], -1.0, atol = 1E-6)
    @test isapprox(opt._current_node.lower_variable_bounds[2], -1.7, atol = 1E-6)
    @test isapprox(opt._current_node.upper_variable_bounds[1], 3.729729, atol = 1E-6)
    @test isapprox(opt._current_node.upper_variable_bounds[2], 1.0, atol = 1E-6)
end

@testset "Optimization-Based Bound Tightening (Nonlinear)" begin
    m = EAGO.Optimizer()
    m._parameters.verbosity = 0

    xL = [0.0; -2.0]
    xU = [4.0; 4.0]

    @variable(m, xL[i] <= x[i=1:2] <= xU[i])
    @objective(m, Min, x[2] - x[1])
    @NLconstraint(m, e1, exp(x[1]/3) - x[2] == 0)
    @NLconstraint(m, e3, x[2] <= 3*log(x[1]+4))

    m.moi_backend.optimizer.model.optimizer._global_lower_bound = Inf
    m.moi_backend.optimizer.model.optimizer._global_upper_bound = -Inf

    optimize!(m)
    opt = m.moi_backend.optimizer.model.optimizer
    m.moi_backend.optimizer.model.optimizer.branch_variable[1] = true
    m.moi_backend.optimizer.model.optimizer.branch_variable[2] = true
    m.moi_backend.optimizer.model.optimizer._global_lower_bound = -10.0
    m.moi_backend.optimizer.model.optimizer._global_upper_bound = 10.0
    oldn = pop!(opt._stack)
    opt._current_node = NodeBB(oldn.lower_variable_bounds, oldn.upper_variable_bounds,-10.0,10.0,1,1)
    opt.obbt_variable_values[1] = true
    opt.obbt_variable_values[2] = true
    feas = EAGO.obbt(opt)

    @test feas
    @test opt._current_node.lower_variable_bounds[1] == 0.0
    @test isapprox(opt._current_node.upper_variable_bounds[1], 4.0, atol = 1E-6)
    @test isapprox(opt._current_node.lower_variable_bounds[2], 0.6492446803515586, atol = 1E-6)
    @test isapprox(opt._current_node.upper_variable_bounds[2], 3.7936678946831783, atol = 1E-6)
end
=#
