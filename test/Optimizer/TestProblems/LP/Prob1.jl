@testset "LP Problem #1" begin
    m = Model(with_optimizer(EAGO.Optimizer))
    #m = Model(with_optimizer(Clp.Optimizer))

    @variable(m, 1 <= x <= 3)
    @variable(m, 1 <= y <= 3)

    @NLobjective(m, Min, x + y)

    @NLconstraint(m, x + y <= 10)
    @NLconstraint(m, x - y <= 10)
    @NLconstraint(m, y >= 0)

    JuMP.optimize!(m)

    xval = JuMP.value(x)
    yval = JuMP.value(y)

    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(xval,1.0,atol=1E-4)
    @test isapprox(yval,1.0,atol=1E-4)
    @test isapprox(fval,2.0,atol=1E-4)
    @test status_term == MOI.OPTIMAL
    @test status_prim == MOI.FEASIBLE_POINT
end
println("end problem 1")
