println("----- Test Example 3 -----")
m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))

@variable(m, -3 <= x <= -1)
@variable(m, -2 <= y <= 2)
@variable(m, 1 <= z <= 3)
@variable(m, -10 <= q <= 9)

@NLobjective(m, Min, 2x - 3y + 2z)

@NLconstraint(m, x + 2y >= -10)
@NLconstraint(m, z - 2y <= 2)
@NLconstraint(m, y >= 0)
@NLconstraint(m, q-3*z-y >= 0)

JuMP.optimize!(m)

xval = JuMP.value(x)
yval = JuMP.value(y)
zval = JuMP.value(z)
qval = JuMP.value(q)

fval = JuMP.objective_value(m)
status_term = JuMP.termination_status(m)
status_prim = JuMP.primal_status(m)

println("xval: $xval")
println("yval: $yval")
println("zval: $zval")
println("qval: $qval")
println("fval: $fval")
println("status_term: $status_term")
println("status_prim: $status_prim")

@testset "LP Problem #3" begin
    println("----- Test Example 3 -----")
    m = Model(with_optimizer(EAGO.Optimizer))
    #m = Model(with_optimizer(Clp.Optimizer))

    @variable(m, -3 <= x <= -1)
    @variable(m, -2 <= y <= 2)
    @variable(m, 1 <= z <= 3)
    @variable(m, -10 <= q <= 9)

    @NLobjective(m, Min, 2x - 3y + 2z)

    @NLconstraint(m, x + 2y >= -10)
    @NLconstraint(m, z - 2y <= 2)
    @NLconstraint(m, y >= 0)
    @NLconstraint(m, q-3*z-y >= 0)

    JuMP.optimize!(m)

    xval = JuMP.value(x)
    yval = JuMP.value(y)
    zval = JuMP.value(z)
    qval = JuMP.value(q)

    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(xval,-3.0,atol=1E-4)
    @test isapprox(yval,2.0,atol=1E-4)
    @test isapprox(zval,1.0,atol=1E-4)
    @test isapprox(fval,-10.0,atol=1E-4)
    @test status_term == MOI.OPTIMAL
    @test status_prim == MOI.FEASIBLE_POINT
end
