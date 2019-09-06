@testset "Linear Relaxations" begin
    src = EAGO.Optimizer()
    model = EAGO.Optimizer()

    x = MOI.add_variables(model,3)
    x = MOI.add_variables(src,3)

    func1 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([5.0,-2.3],[x[1],x[2]]),2.0)
    func2 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([4.0,-2.2],[x[2],x[3]]),2.1)
    func3 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([3.0,-3.3],[x[1],x[3]]),2.2)

    set1 = MOI.LessThan{Float64}(1.0)
    set2 = MOI.GreaterThan{Float64}(2.0)
    set3 = MOI.EqualTo{Float64}(3.0)

    MOI.add_constraint(src, func1, set1)
    MOI.add_constraint(src, func2, set2)
    MOI.add_constraint(src, func3, set3)

    EAGO.relax_linear!(src,model)

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
    @test model.linear_eq_constraints[1][3] == 2
    @test model.linear_eq_constraints[1][3] == 2
    @test model.linear_eq_constraints[1][3] == 2
end

@testset "Quadratic Classification" begin
    vi1 = MOI.VariableIndex(1)
    vi2 = MOI.VariableIndex(2)
    vi3 = MOI.VariableIndex(3)

    f1a_terms = [MOI.ScalarAffineTerm{Float64}(1.0,vi1)]
    f2a_terms = [MOI.ScalarAffineTerm{Float64}(1.0,vi2)]
    f3a_terms = MOI.ScalarAffineTerm{Float64}[]
    f4a_terms = MOI.ScalarAffineTerm{Float64}[]
    f5a_terms = MOI.ScalarAffineTerm{Float64}[]

    f1_quad = [MOI.ScalarQuadraticTerm{Float64}(1.0,vi1,vi1)]
    f2_quad = [MOI.ScalarQuadraticTerm{Float64}(-2.0,vi2,vi2)]
    f3_quad = [MOI.ScalarQuadraticTerm{Float64}(2.0,vi1,vi1), MOI.ScalarQuadraticTerm{Float64}(1.0,vi2,vi2)]
    f4_quad = [MOI.ScalarQuadraticTerm{Float64}(-2.0,vi1,vi1), MOI.ScalarQuadraticTerm{Float64}(-1.0,vi2,vi2)]
    f5_quad = [MOI.ScalarQuadraticTerm{Float64}(-2.0,vi1,vi1), MOI.ScalarQuadraticTerm{Float64}(1.0,vi2,vi2)]

    func1 = MOI.ScalarQuadraticFunction{Float64}(f1a_terms,f1_quad,2.0)     # convex (univariate)
    func2 = MOI.ScalarQuadraticFunction{Float64}(f2a_terms,f2_quad,-2.0)    # concave (univariate)
    func3 = MOI.ScalarQuadraticFunction{Float64}(f3a_terms,f3_quad,-1.0)    # convex (multivariate)
    func4 = MOI.ScalarQuadraticFunction{Float64}(f4a_terms,f4_quad,3.0)     # concave (multivariate)
    func5 = MOI.ScalarQuadraticFunction{Float64}(f5a_terms,f5_quad,-4.0)    # saddel point (multivariate)

    func1_pos = EAGO.is_convex_quadratic(func1,1,1.0)
    func1_neg = EAGO.is_convex_quadratic(func1,1,-1.0)
    func2_pos = EAGO.is_convex_quadratic(func2,2,1.0)
    func2_neg = EAGO.is_convex_quadratic(func2,2,-1.0)
    func3_pos = EAGO.is_convex_quadratic(func3,2,1.0)
    func3_neg = EAGO.is_convex_quadratic(func3,2,-1.0)
    func4_pos = EAGO.is_convex_quadratic(func4,2,1.0)
    func4_neg = EAGO.is_convex_quadratic(func4,2,-1.0)
    func5_pos = EAGO.is_convex_quadratic(func5,2,1.0)
    func5_neg = EAGO.is_convex_quadratic(func5,2,-1.0)

    @test func1_pos == true
    @test func1_neg == false
    @test func2_pos == false
    @test func2_neg == true
    @test func3_pos == true
    @test func3_neg == false
    @test func4_pos == false
    @test func4_neg == true
    @test func5_pos == false
    @test func5_neg == false

    # Classify Quadratics in model model
    m1 = EAGO.Optimizer()
    MOI.add_variables(m1,2)
    c1i1 = MOI.add_constraint(m1, func1, MOI.LessThan{Float64}(1.0))
    c1i2 = MOI.add_constraint(m1, func2, MOI.GreaterThan{Float64}(2.0))
    c1i3 = MOI.add_constraint(m1, func3, MOI.LessThan{Float64}(3.0))
    c1i4 = MOI.add_constraint(m1, func4, MOI.GreaterThan{Float64}(-4.0))
    c1i5 = MOI.add_constraint(m1, func5, MOI.EqualTo{Float64}(-5.0))
    EAGO.quadratic_convexity!(m1)
    const11_convex = m1.constraint_convexity[c1i1]
    const12_convex = m1.constraint_convexity[c1i2]
    const13_convex = m1.constraint_convexity[c1i3]
    const14_convex = m1.constraint_convexity[c1i4]
    const15_convex = m1.constraint_convexity[c1i5]

    # Classify Quadratics on re-ordered model
    m2 = EAGO.Optimizer()
    MOI.add_variables(m2,2)
    c2i1 = MOI.add_constraint(m2, func5, MOI.LessThan{Float64}(1.0))
    c2i2 = MOI.add_constraint(m2, func4, MOI.GreaterThan{Float64}(2.0))
    c2i3 = MOI.add_constraint(m2, func3, MOI.LessThan{Float64}(3.0))
    c2i4 = MOI.add_constraint(m2, func1, MOI.GreaterThan{Float64}(-4.0))
    c2i5 = MOI.add_constraint(m2, func2, MOI.EqualTo{Float64}(-5.0))
    EAGO.quadratic_convexity!(m2)
    const21_convex = m2.constraint_convexity[c2i1]
    const22_convex = m2.constraint_convexity[c2i2]
    const23_convex = m2.constraint_convexity[c2i3]
    const24_convex = m2.constraint_convexity[c2i4]
    const25_convex = m2.constraint_convexity[c2i5]

    @test const11_convex
    @test const12_convex
    @test const13_convex
    @test const14_convex
    @test ~const15_convex
    @test ~const21_convex
    @test const22_convex
    @test const23_convex
    @test ~const24_convex
    @test ~const25_convex
end

#=
@testset "Quadratic Relaxations" begin
    model = EAGO.Optimizer()
    x = MOI.add_variables(model,3)

        func1 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(5.0,x[1])],
                                                 [MOI.ScalarQuadraticTerm{Float64}(2.5,x[2],x[2])],2.0)
        func2 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(4.0,x[2])],
                                                 [MOI.ScalarQuadraticTerm{Float64}(-2.2,x[1],x[2])],2.1)
        func3 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(3.0,x[3])],
                                                 [MOI.ScalarQuadraticTerm{Float64}(2.1,x[1],x[1])],2.2)
        func4 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(3.0,x[3])],
                                                 [MOI.ScalarQuadraticTerm{Float64}(-2.1,x[1],x[1])],2.2)
        func5 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(3.0,x[3])],
                                                 [MOI.ScalarQuadraticTerm{Float64}(2.1,x[1],x[1])],2.2)

        set1 = MOI.LessThan{Float64}(1.0)
        set2 = MOI.GreaterThan{Float64}(2.0)
        set3 = MOI.EqualTo{Float64}(3.0)

        MOI.add_constraint(model, func1, set1)
        MOI.add_constraint(model, func2, set1)
        MOI.add_constraint(model, func3, set2)
        MOI.add_constraint(model, func4, set2)
        MOI.add_constraint(model, func5, set3)

        EAGO.Quadratic_Convexity!(model)

        @test model.constraint_convexity[MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.GreaterThan{Float64}}(1)] == false
        @test model.constraint_convexity[MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.GreaterThan{Float64}}(2)]
        @test model.constraint_convexity[MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.LessThan{Float64}}(1)]
        @test model.constraint_convexity[MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.LessThan{Float64}}(2)] == false
        @test model.constraint_convexity[MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.EqualTo{Float64}}(1)] == false

        target = EAGO.Optimizer()
        x = MOI.add_variables(target,3)

        model.VItoSto[1] = 1; model.VItoSto[2] = 2; model.VItoSto[3] = 3;
        model.StoToVI = EAGO.ReverseDict(model.VItoSto)

        n = EAGO.NodeBB(Float64[1.0,5.0,8.0], Float64[2.0,6.0,9.0], -Inf, Inf, 2, 1, true)

        r = EAGO.DefaultRelaxationScheme()

        EAGO.RelaxQuadratic!(target, model, n, r)

        @test target.linear_leq_constraints[1][1].terms[1].coefficient == 5.0
        @test target.linear_leq_constraints[1][1].terms[1].variable_index.value == 1
        @test target.linear_leq_constraints[1][1].terms[2].coefficient == 2.5
        @test target.linear_leq_constraints[1][1].terms[2].variable_index.value == 2
        @test target.linear_leq_constraints[1][1].constant == 1.375
        @test target.linear_leq_constraints[1][2].upper == 1.0
        @test target.linear_leq_constraints[1][3] == 2

        @test target.linear_leq_constraints[2][1].terms[1].coefficient == -11.0
        @test target.linear_leq_constraints[2][1].terms[1].variable_index.value == 1
        @test isapprox(target.linear_leq_constraints[2][1].terms[2].coefficient, -0.4, atol=1E7)
        @test target.linear_leq_constraints[2][1].terms[2].variable_index.value == 2
        @test target.linear_leq_constraints[2][1].constant == 24.1
        @test target.linear_leq_constraints[2][2].upper == 1.0
        @test target.linear_leq_constraints[2][3] == 2

        @test isapprox(target.linear_leq_constraints[3][1].terms[1].coefficient, -0.63, atol=1E7)
        @test target.linear_leq_constraints[3][1].terms[1].variable_index.value == 1
        @test target.linear_leq_constraints[3][1].terms[2].coefficient == -3.0
        @test target.linear_leq_constraints[3][1].terms[2].variable_index.value == 3
        @test target.linear_leq_constraints[3][1].constant == 2.0
        @test target.linear_leq_constraints[3][2].upper == -2.0
        @test target.linear_leq_constraints[3][3] == 2

        @test target.linear_leq_constraints[4][1].terms[1].coefficient == 2.1
        @test target.linear_leq_constraints[4][1].terms[1].variable_index.value == 1
        @test target.linear_leq_constraints[4][1].terms[2].coefficient == -3.0
        @test target.linear_leq_constraints[4][1].terms[2].variable_index.value == 3
        @test target.linear_leq_constraints[4][1].constant == -2.725
        @test target.linear_leq_constraints[4][2].upper == -2.0
        @test target.linear_leq_constraints[4][3] == 2

        @test isapprox(target.linear_leq_constraints[5][1].terms[1].coefficient, 6.3, atol=1E7)
        @test target.linear_leq_constraints[5][1].terms[1].variable_index.value == 1
        @test target.linear_leq_constraints[5][1].terms[2].coefficient == 3.0
        @test target.linear_leq_constraints[5][1].terms[2].variable_index.value == 3
        @test target.linear_leq_constraints[5][1].constant == -2.0
        @test target.linear_leq_constraints[5][2].upper == 3.0
        @test target.linear_leq_constraints[5][3] == 2

        @test isapprox(target.linear_leq_constraints[6][1].terms[1].coefficient, -6.3, atol=1E7)
        @test target.linear_leq_constraints[6][1].terms[1].variable_index.value == 1
        @test target.linear_leq_constraints[6][1].terms[2].coefficient == -3.0
        @test target.linear_leq_constraints[6][1].terms[2].variable_index.value == 3
        @test target.linear_leq_constraints[6][1].constant == 2.0
        @test target.linear_leq_constraints[6][2].upper == -3.0
        @test target.linear_leq_constraints[6][3] == 2
    end
    =#
