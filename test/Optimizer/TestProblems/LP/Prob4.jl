@testset "LP Problem #4" begin
    # Test LP #4
    println("----- Test Example 4 -----")
    m = Model(with_optimizer(EAGO.Optimizer))
    #m = Model(with_optimizer(Clp.Optimizer))

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

    xval = JuMP.value(x)
    yval = JuMP.value(y)
    zval = JuMP.value(z)
    qval = JuMP.value(q)

    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(xval,0.0,atol=1E-4)
    @test isapprox(yval,0.0,atol=1E-4)
    @test isapprox(zval,0.0,atol=1E-4)
    @test isapprox(qval,0.0,atol=1E-4)
    @test isapprox(fval,Inf,atol=1E-4)
    @test status_term == MOI.INFEASIBLE
    @test status_prim == MOI.INFEASIBILITY_CERTIFICATE
end
