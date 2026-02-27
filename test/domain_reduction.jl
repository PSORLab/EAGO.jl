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
    # lvb[3] doesn't tighten since dbbt assumes variables aren't fixed
    @test isapprox(lvb[4], 3.33333333, atol= 1E-5); @test uvb[4] == 4.0
end

@testset "Constraint Propagation" begin

    # SQRT
    model = JuMP.Model(EAGO.Optimizer)

    @variable(model, 0 <= x )
    @variable(model, 0 <= y )
    @constraint(model, sqrt(x) + y <= 4)

    @objective(model, Min, x^2 + y^2)
    JuMP.optimize!(model)

    @test isapprox(JuMP.value(x), 0.000297, atol=1E-6)
    @test isapprox(JuMP.value(y), 0.000297, atol=1E-6)
    
    # POW
    model = JuMP.Model(EAGO.Optimizer)
    @variable(model, x)
    @constraint(model, x^1.852 <= 1)
    JuMP.optimize!(model)

    @test 0.0 <= JuMP.value(x) <= 1.0

    # LOG
    model = JuMP.Model(EAGO.Optimizer)
    @variable(model, x <= 2)
    @variable(model, t)
    @constraint(model, log(x) >= t)
    @objective(model, Max, t)
    JuMP.optimize!(model)

    @test isapprox(JuMP.value(x), 2, atol=1E-3)
    @test isapprox(JuMP.value(t), log(2), atol=1E-3)

    # SIN 
    model = JuMP.Model(EAGO.Optimizer)
    @variable(model, 0 <= x <= 2pi)
    @variable(model, y <= 3)
    @constraint(model, sin(x) + y >= 1.5)
    @objective(model, Min, x + y)
    JuMP.optimize!(model)

    @test isapprox(JuMP.value(x), 0.00567, atol=1E-3)
    @test isapprox(JuMP.value(y), 1.4943, atol=1E-3)

    # ACOS
    model = JuMP.Model(EAGO.Optimizer)
    @variable(model, -2 <= x <= 2)
    @constraint(model, acos(x) <= 1)
    JuMP.optimize!(model)

    @test 0.540302 <= JuMP.value(x) <= 1.0
end