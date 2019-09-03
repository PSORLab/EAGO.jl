@testset "Duality-Based Bound Tightening" begin

    ylower = [1.0, 1.0, 1.0, 1.0]
    yupper = [4.0, 4.0, 4.0, 4.0]
    ymult_lo = [50, 1.0, 2.0, 3.0]
    ymult_hi = [0, 1.0, 2.0, 3.0]
    yLBD = 1.0
    yUBD = 3.0
    n = EAGO.NodeBB(ylower,yupper,yLBD,yUBD,3,2,false)

    EAGO.variable_duality_based_tightening!(n,ymult_lo,ymult_hi,yLBD,yUBD)
    @test 3.95999-1E-4 <= n.lower_variable_bounds[1] <= 3.95999+1E-4
    @test 4.0-1E-4 <= n.upper_variable_bounds[1] <= 4.0+1E-4
    @test 2.0-1E-4 <= n.lower_variable_bounds[2] <= 2.0+1E-4
    @test 4.0-1E-4 <= n.upper_variable_bounds[2] <= 4.0+1E-4
    @test 3.0-1E-4 <= n.lower_variable_bounds[3] <= 3.0+1E-4
    @test 4.0-1E-4 <= n.upper_variable_bounds[3] <= 4.0+1E-4
    @test 3.33333-1E-4 <= n.lower_variable_bounds[4] <= 3.33333+1E-4
    @test 4.0-1E-4 <= n.upper_variable_bounds[4] <= 4.0+1E-4
end

@testset "Poor Man's LP" begin
    # Puranik 2017 example
    opt1 = EAGO.Optimizer()
    x = MOI.add_variables(opt1,3)

    MOI.add_constraint(opt1,MOI.SingleVariable(x[1]), MOI.LessThan(4.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[1]), MOI.GreaterThan(-2.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[2]), MOI.LessThan(4.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[2]), MOI.GreaterThan(0.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[3]), MOI.LessThan(1.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[3]), MOI.GreaterThan(-1.0))

    MOI.add_constraint(opt1,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, x[1]), MOI.ScalarAffineTerm(1.0, x[2])], -2.0), MOI.GreaterThan(2.0))
    MOI.add_constraint(opt1,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, x[2]), MOI.ScalarAffineTerm(1.0, x[3])], -1.0), MOI.LessThan(0.0))

    n1 = EAGO.NodeBB()
    n1.lower_variable_bounds = [-2.0, 0.0, -1.0]
    n1.upper_variable_bounds = [4.0, 4.0, 1.0]
    n1.lower_bound = -Inf
    n1.upper_bound = Inf
    n1.depth = 3

    opt1.variable_index_to_storage = Dict{Int,Int}(1 => 1, 2 => 2, 3 => 3)
    opt1.bisection_variable = Dict{Int,Int}(1 => true, 2 => true, 3 => true)

    feas1 = EAGO.poor_man_lp(opt1,n1)
    feas1 = EAGO.poor_man_lp(opt1,n1)
    @test (n1.lower_variable_bounds == [2.0, 0.0, -1.0])
    @test (n1.upper_variable_bounds == [4.0, 2.0, 1.0])
    @test (feas1 == true)
end

@testset "Classify Quadratic Types" begin
    opt1 = EAGO.Optimizer()
    x = MOI.add_variables(opt1,3)

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

    # allocates special storage for these forms
    EAGO.classify_quadratics!(opt1)

    # checks to see if classification was correct
    @test opt1.univariate_quadratic_eq_constraints[1][1] == 1.5
    @test opt1.univariate_quadratic_eq_constraints[1][2] == 1.4
    @test opt1.univariate_quadratic_eq_constraints[1][3] == 1.0
    @test opt1.univariate_quadratic_eq_constraints[1][4] == 2

    @test opt1.univariate_quadratic_leq_constraints[1][1] == -1.5
    @test opt1.univariate_quadratic_leq_constraints[1][2] == -1.2
    @test opt1.univariate_quadratic_leq_constraints[1][3] == 1.0
    @test opt1.univariate_quadratic_leq_constraints[1][4] == 1

    @test opt1.univariate_quadratic_geq_constraints[1][1] == 1.5
    @test opt1.univariate_quadratic_geq_constraints[1][2] == 1.1
    @test opt1.univariate_quadratic_geq_constraints[1][3] == 0.0
    @test opt1.univariate_quadratic_geq_constraints[1][4] == 1

    @test opt1.univariate_quadratic_geq_constraints[2][1] == 1.5
    @test opt1.univariate_quadratic_geq_constraints[2][2] == 1.3
    @test opt1.univariate_quadratic_geq_constraints[2][3] == 4.0
    @test opt1.univariate_quadratic_geq_constraints[2][4] == 2
end

@testset "Quadratic Domain Reduction (Univariate)" begin
    m = EAGO.Optimizer()

    n2 = NodeBB()
    n2.lower_variable_bounds = [-10.0, -10.0, -10.0]
    n2.upper_variable_bounds = [10.0, 10.0, 10.0]
    n2.lower_bound = -Inf
    n2.upper_bound = Inf
    n2.depth = 1
    a1, b1, c1 = 2.0, -4.0, 3.0
    a2, b2, c2 = 2.0, -4.0, 3.0
    a3, b3, c3 = 2.0, -4.0, 3.0
    push!(m.univariate_quadratic_geq_constraints,(a1,b1,c1,1))
    push!(m.univariate_quadratic_leq_constraints,(a2,b2,c2,2))
    push!(m.univariate_quadratic_eq_constraints,(a3,b3,c3,3))
    feas = EAGO.univariate_quadratic(m,n2)

    @test feas == true
    @test isapprox(n2.lower_variable_bounds[1],-10.0,atol=1E-3)
    @test isapprox(n2.lower_variable_bounds[2],-0.58113,atol=1E-3)
    @test isapprox(n2.lower_variable_bounds[3],-0.58113,atol=1E-3)
    @test isapprox(n2.upper_variable_bounds[1],10.0,atol=1E-3)
    @test isapprox(n2.upper_variable_bounds[2],2.5811,atol=1E-3)
    @test isapprox(n2.upper_variable_bounds[3],2.5811,atol=1E-3)
end

@testset "Quadratic Domain Reduction (Bivariate)" begin
end

#=
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

@testset "Optimization-Based Bound Tightening (Linear)" begin
    m = Model(with_optimizer(EAGO.Optimizer))
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

    m.moi_backend.optimizer.model.optimizer.global_lower_bound = Inf
    m.moi_backend.optimizer.model.optimizer.global_upper_bound = -Inf

    optimize!(m)
    opt = m.moi_backend.optimizer.model.optimizer
    m.moi_backend.optimizer.model.optimizer.bisection_variable[1] = true
    m.moi_backend.optimizer.model.optimizer.bisection_variable[2] = true
    m.moi_backend.optimizer.model.optimizer.obbt_variables = [MOI.VariableIndex(1), MOI.VariableIndex(2)]
    m.moi_backend.optimizer.model.optimizer.global_lower_bound = -12.0
    m.moi_backend.optimizer.model.optimizer.global_upper_bound = 0.0
    n = opt.stack[1]
    n.lower_bound = -12.0
    n.upper_bound = 0.0
    feas = EAGO.obbt(opt, n)

    @test feas
    @test n.lower_variable_bounds[1] == -1.0
    @test isapprox(n.lower_variable_bounds[2], -1.7, atol = 1E-6)
    @test isapprox(n.upper_variable_bounds[1], 3.729729, atol = 1E-6)
    @test n.upper_variable_bounds[2] == 1.0
end

@testset "Optimization-Based Bound Tightening (Nonlinear)" begin
    m = Model(with_optimizer(EAGO.Optimizer,
                             udf_scrubber_flag = false,
                             udf_to_JuMP_flag = false))
    xL = [0.0; -2.0]
    xU = [4.0; 4.0]

    @variable(m, xL[i] <= x[i=1:2] <= xU[i])
    @objective(m, Min, x[2] - x[1])
    @NLconstraint(m, e1, exp(x[1]/3) - x[2] == 0)
    @NLconstraint(m, e3, x[2] <= 3*log(x[1]+4))

    m.moi_backend.optimizer.model.optimizer.global_lower_bound = Inf
    m.moi_backend.optimizer.model.optimizer.global_upper_bound = -Inf

    optimize!(m)
    opt = m.moi_backend.optimizer.model.optimizer
    m.moi_backend.optimizer.model.optimizer.bisection_variable[1] = true
    m.moi_backend.optimizer.model.optimizer.bisection_variable[2] = true
    m.moi_backend.optimizer.model.optimizer.obbt_variables = [MOI.VariableIndex(1), MOI.VariableIndex(2)]
    m.moi_backend.optimizer.model.optimizer.global_lower_bound = -10.0
    m.moi_backend.optimizer.model.optimizer.global_upper_bound = 10.0
    n = opt.stack[1]
    n.lower_bound = -10.0
    n.upper_bound = 10.0
    feas = EAGO.obbt(opt, n)

    @test feas

    @test n.lower_variable_bounds[1] == 0.0
    @test isapprox(n.upper_variable_bounds[1], 4.0, atol = 1E-6)

    @test isapprox(n.lower_variable_bounds[2], 0.6492446803515586, atol = 1E-6)
    @test isapprox(n.upper_variable_bounds[2], 3.7936678946831783, atol = 1E-6)
end
