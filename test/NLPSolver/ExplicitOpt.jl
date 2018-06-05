module Test_Explicit_NLP

using Compat
using Compat.Test
using EAGO
using JuMP
using MathProgBase


@testset "JuMP Interface Explicit (LP)" begin

  jumpmodel4 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                           LBDsolvertype = "LP",
                                           probe_depth = -1,
                                           variable_depth = 1000,
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

  jumpmodel5 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
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

  jumpmodel5a = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                           LBDsolvertype = "LP",
                                           probe_depth = -1,
                                           variable_depth = 1000,
                                           DAG_depth = 10,
                                           STD_RR_depth = 1,
                                           validated = true))
  @variable(jumpmodel5a, -200 <= x2a <= -100)
  @variable(jumpmodel5a, 200 <= y2a <= 400)
  @constraint(jumpmodel5a, -500 <= x2a+2y2a <= 400)
  @NLobjective(jumpmodel5a, Min, x2a*y2a)
  status5 = solve(jumpmodel5a)


  @test status5 == :Optimal
  @test isapprox(getvalue(x2),-200.0,atol=1E-6)
  @test isapprox(getvalue(y2),300.0,atol=1E-6)
  @test isapprox(getobjectivevalue(jumpmodel5),-60000.00119999499,atol=2.0)

  println("test jumpmodel 6")
  jumpmodel6 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                         LBDsolvertype = "LP",
                                         probe_depth = -1,
                                         variable_depth = 1000,
                                         DAG_depth = 10,
                                         STD_RR_depth = -1,
                                         validated = true))
  @variable(jumpmodel6, -5 <= x1 <= 5)
  @variable(jumpmodel6, -5 <= y1 <= 5)
  @NLobjective(jumpmodel6, Min, 2*x1^2-1.05*x1^4+(x1^6)/6+x1*y1+y1^2)
  status6 = solve(jumpmodel6)

  @test isapprox(getvalue(x1),0.0,atol=1E-3)
  @test isapprox(getvalue(y1),0.0,atol=1E-3)
  @test isapprox(getobjectivevalue(jumpmodel6),0.0,atol=1E-6)
  @test status6 == :Optimal

  println("test jumpmodel 6b")
  jumpmodel6b = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Interval",
                                         LBDsolvertype = "Interval",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = 10,
                                         STD_RR_depth = -1,
                                         validated = true,
                                         atol=1E-1))
  @variable(jumpmodel6b, -5 <= x1b <= 5)
  @variable(jumpmodel6b, -5 <= y1b <= 5)
  @NLobjective(jumpmodel6b, Min, 2*x1b^2-1.05*x1b^4+(x1b^6)/6+x1b*y1b+y1b^2)
  status6b = solve(jumpmodel6b)

  @test isapprox(getvalue(x1b),0.0,atol=1E0)
  @test isapprox(getvalue(y1b),0.0,atol=1E0)
  @test isapprox(getobjectivevalue(jumpmodel6b),0.0,atol=1E-1)
  @test status6b == :Optimal

  println("test jumpmodel 8")
  jumpmodel8 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Interval",
                                         LBDsolvertype = "Interval",
                                           probe_depth = -1,
                                           variable_depth = -1,
                                           DAG_depth = -1,
                                           STD_RR_depth = -1,
                                           atol=1E-1,
                                           validated = true))
  @variable(jumpmodel8, -200 <= x2b <= -100)
  @variable(jumpmodel8, 200 <= y2b <= 400)
  @constraint(jumpmodel8, -500 <= x2b+2y2b <= 400)
  @NLobjective(jumpmodel8, Min, x2b*y2b)
  status4 = solve(jumpmodel8)

  @test status4 == :Optimal
  @test isapprox(getvalue(x),-200.0,atol=1E-1)
  @test isapprox(getvalue(y),300.0,atol=1E-1)
  @test isapprox(getobjectivevalue(jumpmodel8),-60000.00119999499,atol=2.0)

  println("test jumpmodel 8a")
  jumpmodel8a = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                         LBDsolvertype = "LP",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = 10,
                                         STD_RR_depth = -1,
                                         validated = true,
                                         atol=1E-1))
  @variable(jumpmodel8a, -5 <= x3b <= 5)
  @variable(jumpmodel8a, -5 <= y3b <= 5)
  @NLobjective(jumpmodel8a, Min, 2*x3b^2-1.05*x3b^4+(x3b^6)/6+x3b*y3b+y3b^2)
  status6b = solve(jumpmodel8a)
  @test isapprox(getvalue(x3b),0.0,atol=1E0)
  @test isapprox(getvalue(y3b),0.0,atol=1E0)
  @test isapprox(getobjectivevalue(jumpmodel8a),0.0,atol=1E-1)
  @test status6b == :Optimal

  println("test jumpmodel 8b")
  jumpmodel8b = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                         LBDsolvertype = "LP",
                                         UBDsolvertype = "Interval",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         DAG_depth = 10,
                                         STD_RR_depth = -1,
                                         validated = true,
                                         atol=1E-1))
  @variable(jumpmodel8b, -5 <= x3c <= 5)
  @variable(jumpmodel8b, -5 <= y3c <= 5)
  @NLobjective(jumpmodel8b, Min, 2*x3c^2-1.05*x3c^4+(x3c^6)/6+x3c*y3c+y3c^2)
  status6b = solve(jumpmodel8b)
  @test isapprox(getvalue(x3c),0.0,atol=1E0)
  @test isapprox(getvalue(y3c),0.0,atol=1E0)
  @test isapprox(getobjectivevalue(jumpmodel8b),0.0,atol=1E-1)
  @test status6b == :Optimal

  println("test jumpmodel 8c")
  jumpmodel8c = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                           LBDsolvertype = "LP",
                                           probe_depth = -1,
                                           variable_depth = 1000,
                                           DAG_depth = 10,
                                           STD_RR_depth = 1,
                                           validated = true))
  @variable(jumpmodel8c, -200 <= x2b <= -100)
  @variable(jumpmodel8c, 200 <= y2b <= 400)
  @constraint(jumpmodel8c, -500 <= x2b+2y2b <= 400)
  @NLobjective(jumpmodel8c, Min, x2b*y2b)
  status5 = solve(jumpmodel8c)


  println("test jumpmodel 6")
  jumpmodel9 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Diff1-MV",
                                         LBDsolvertype = "Ipopt",
                                         probe_depth = -1,
                                         variable_depth = 1000,
                                         DAG_depth = 10,
                                         STD_RR_depth = -1,
                                         validated = true))
  @variable(jumpmodel9, -5 <= x9 <= 5)
  @variable(jumpmodel9, -5 <= y9 <= 5)
  @NLobjective(jumpmodel9, Min, 2*x9^2-1.05*x9^4+(x9^6)/6+x9*y9+y9^2)
  status9 = solve(jumpmodel9)
  @test isapprox(getvalue(x9),0.0,atol=1E0)
  @test isapprox(getvalue(y9),0.0,atol=1E0)
  @test isapprox(getobjectivevalue(jumpmodel9),0.0,atol=1E-1)
  @test status9 == :Optimal

#=
  println("test jumpmodel 10")
  jumpmodel10 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "Diff1-MV",
                                           LBDsolvertype = "Ipopt",
                                           probe_depth = -1,
                                           variable_depth = 1000,
                                           DAG_depth = -1,
                                           STD_RR_depth = -1,
                                           atol=1E-1,
                                           validated = true))
  @variable(jumpmodel10, -200 <= x10 <= -100)
  @variable(jumpmodel10, 200 <= y10 <= 400)
  @constraint(jumpmodel10, -500 <= x10+2y10 <= 400)
  @NLobjective(jumpmodel10, Min, x10*y10)
  status4 = solve(jumpmodel10)

  @test status10 == :Optimal
  @test isapprox(getvalue(x10),-200.0,atol=1E-1)
  @test isapprox(getvalue(y10),300.0,atol=1E-1)
  @test isapprox(getobjectivevalue(jumpmodel10),-60000.00119999499,atol=2.0)
=#
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
