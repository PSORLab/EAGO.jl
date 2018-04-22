module Test_Explicit_NLP

using Compat
using Compat.Test
using EAGO
using JuMP
using MathProgBase


@testset "JuMP Interface Explicit (LP)" begin

  jumpmodel4 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                           LBDsolvertype = "LP",
                                           probe_depth = -1,
                                           variable_depth = 1000,
                                           DAG_depth = -1,
                                           STD_RR_depth = -1))
  @variable(jumpmodel4, -200 <= x <= -100)
  @variable(jumpmodel4, 200 <= y <= 400)
  @constraint(jumpmodel4, -500 <= x+2y <= 400)
  @NLobjective(jumpmodel4, Min, x*y)
  status4 = solve(jumpmodel4)

  @test status4 == :Optimal
  @test isapprox(getvalue(x),-200.0,atol=1E-6)
  @test isapprox(getvalue(y),300.0,atol=1E-6)
  @test isapprox(getobjectivevalue(jumpmodel4),-60000.00119999499,atol=2.0)

  jumpmodel5 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                           LBDsolvertype = "LP",
                                           probe_depth = -1,
                                           variable_depth = 1000,
                                           DAG_depth = -1,
                                           STD_RR_depth = -1,
                                           validated = false))
  @variable(jumpmodel5, -200 <= x2 <= -100)
  @variable(jumpmodel5, 200 <= y2 <= 400)
  @constraint(jumpmodel5, -500 <= x2+2y2 <= 400)
  @NLobjective(jumpmodel5, Min, x2*y2)
  status5 = solve(jumpmodel5)

  @test status5 == :Optimal
  @test isapprox(getvalue(x2),-200.0,atol=1E-6)
  @test isapprox(getvalue(y2),300.0,atol=1E-6)
  @test isapprox(getobjectivevalue(jumpmodel5),-60000.00119999499,atol=2.0)

  jumpmodel6 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                         LBDsolvertype = "LP",
                                         probe_depth = -1,
                                         variable_depth = 1000,
                                         DAG_depth = -1,
                                         STD_RR_depth = -1))
  @variable(jumpmodel6, -5 <= x1 <= 5)
  @variable(jumpmodel6, -5 <= y1 <= 5)
  @NLobjective(jumpmodel6, Min, 2*x1^2-1.05*x1^4+(x1^6)/6+x1*y1+y1^2)
  status6 = solve(jumpmodel6)

  @test isapprox(getvalue(x1),0.0,atol=1E-3)
  @test isapprox(getvalue(y1),0.0,atol=1E-3)
  @test isapprox(getobjectivevalue(jumpmodel6),0.0,atol=1E-6)
  @test status6 == :Optimal

  end

#=
 @testset "JuMP Interface Explicit (SNOPT)" begin

    jumpmodel4 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Diff1-MV-OFF",
                                             LBDsolvertype = "SNOPT",
                                             probe_depth = -1,
                                             variable_depth = -1,
                                             DAG_depth = -1,
                                             STD_RR_depth = 10))
    @variable(jumpmodel4, -200 <= x <= -100)
    @variable(jumpmodel4, 200 <= y <= 400)
    @constraint(jumpmodel4, -500 <= x+2y <= 400)
    @NLobjective(jumpmodel4, Min, x*y)
    status4 = solve(jumpmodel4)

    @test status4 == :Optimal
    @test isapprox(getvalue(x),-200.0,atol=1E-6)
    @test isapprox(getvalue(y),300.0,atol=1E-6)
    @test isapprox(getobjectivevalue(jumpmodel4),-60000.00119999499,atol=2.0)

    end
=#

    #=
@testset "JuMP Interface Explicit (SNOPT)" begin
jumpmodel5 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Diff1-MV-OFF",
                                         LBDsolvertype = "SNOPT",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = 10,
                                         STD_RR_depth = 10))
@variable(jumpmodel5, 0 <= x <= 400)
@variable(jumpmodel5, 0 <= y <= 200)
@constraint(jumpmodel5, x+2y == 500)
@NLobjective(jumpmodel5, Min, x*y)
status5 = solve(jumpmodel5)

@test status5 == :Optimal

jumpmodel6 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Diff1-MV-OFF",
                                         LBDsolvertype = "SNOPT",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = 10,
                                         STD_RR_depth = 10))
@variable(jumpmodel6, -5 <= x <= 5)
@variable(jumpmodel6, -5 <= y <= 5)
@NLobjective(jumpmodel6, Min, 2*x^2-1.05*x^4+(x^6)/6+x*y+y^2)
status4 = solve(jumpmodel6)

@test status6 == :Optimal
end

@testset "JuMP Interface Explicit (Interval)" begin
jumpmodel5 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Interval",
                                         LBDsolvertype = "Interva",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = 10,
                                         STD_RR_depth = 10))
@variable(jumpmodel5, 0 <= x <= 400)
@variable(jumpmodel5, 0 <= y <= 200)
@constraint(jumpmodel5, x+2y == 500)
@NLobjective(jumpmodel5, Min, x*y)
status5 = solve(jumpmodel5)

jumpmodel6 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Interval",
                                         LBDsolvertype = "Interva",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = 10,
                                         STD_RR_depth = 10))
@variable(jumpmodel6, -5 <= x <= 5)
@variable(jumpmodel6, -5 <= y <= 5)
@NLobjective(jumpmodel6, Min, 2*x^2-1.05*x^4+(x^6)/6+x*y+y^2)
status4 = solve(jumpmodel6)
end

@testset "MathProgBase Interface Explicit" begin
end
=#
end
