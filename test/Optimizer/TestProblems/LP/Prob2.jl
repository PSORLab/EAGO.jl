# Test LP #2
@testset "LP Problem #2" begin
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
end
