#jumpmodel6 = Model(with_optimizer(EAGO.Optimizer))
#@variable(jumpmodel6, -5 <= x1 <= 5)
#@variable(jumpmodel6, -5 <= y1 <= 5)
#@NLobjective(jumpmodel6, Min, 2*x1^2-1.05*x1^4+(x1^6)/6+x1*y1+y1^2)
#status6 = JuMP.optimize!(jumpmodel6)

#=
@testset "NLP Problem #2" begin
    jumpmodel6 = Model(with_optimizer(EAGO.Optimizer))
    @variable(jumpmodel6, -5 <= x1 <= 5)
    @variable(jumpmodel6, -5 <= y1 <= 5)
    @NLobjective(jumpmodel6, Min, 2*x1^2-1.05*x1^4+(x1^6)/6+x1*y1+y1^2)
    status6 = JuMP.optimize!(jumpmodel6)

    @test isapprox(getvalue(x1),0.0,atol=1E-3)
    @test isapprox(getvalue(y1),0.0,atol=1E-3)
    @test isapprox(getobjectivevalue(jumpmodel6),0.0,atol=1E-6)
    @test status6 == :Optimal
end
=#
