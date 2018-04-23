module IntervalArithmetic_Check

using Compat
using Compat.Test
using EAGO

a = MCInterval(2.0,3.0)
b = MCInterval(4.0,6.0)

@test a != b
@test a <= b
@test a < b
@test ~(a >= b)
@test ~(a > b)

@testset "infty" begin
    @test EAGO.infty(Float64) == Inf64
    @test EAGO.ninfty(Float64) == -Inf64
    @test EAGO.infty(Float32) == Inf32
    @test EAGO.ninfty(Float32) == -Inf32
    @test EAGO.infty(Float16) == Inf16
    @test EAGO.ninfty(Float16) == -Inf16
end

@testset "add/sub/div/mult" begin
    a = MCInterval(0.1, 1.1)
    b = MCInterval(0.9, 2.0)
    c = MCInterval(0.25, 4.0)


    ## Basic arithmetic
    @test a == MCInterval(0.1, 1.1)
    @test +a == a
    @test a+b == MCInterval(1.0, 3.1)
    @test -a == MCInterval(-1.1, -0.1)
    @test a-b == MCInterval(-1.9, 0.20000000000000007)

    @test MCInterval(-30.0,-15.0) / MCInterval(-5.0,-3.0) == MCInterval(3.0, 10.0)
    @test b/a == MCInterval(0.8181818181818181, 20.0)
    @test a/c == MCInterval(0.025, 4.4)
    @test c/4.0 == MCInterval(6.25e-02, 1e+00)
    @test c/zero(c) == emptyMCinterval(c)

    @test MCInterval(0.0, 1.0)/MCInterval(0.0,1.0) == MCInterval(0.0, Inf)
    @test MCInterval(-1.0, 1.0)/MCInterval(0.0,1.0) == entireMCinterval(c)
    @test MCInterval(-1.0, 1.0)/MCInterval(-1.0,1.0) == entireMCinterval(c)

    x = MCInterval(1.0,2.0)

    @test 0.1 + x == MCInterval(1.1, 2.1)
    @test 3.0 - x == MCInterval(1.0, 2.0)
    @test 3.1 - x == MCInterval(1.1, 2.1)
    @test 0.1 * MCInterval(1.0, 1.0) == MCInterval(0.1, 0.1)
    @test (-0.1) * MCInterval(1.0, 1.0) == MCInterval(-0.1, -0.1)
    @test MCInterval(1.0, 1.0) * 0.1 == MCInterval(0.1, 0.1)
    @test MCInterval(1.0, 1.0) * 1 == MCInterval(1.0, 1.0)
    @test MCInterval(1.0, 1.0) / 10.0 == MCInterval(0.1, 0.1)
end


@testset "Power tests" begin
    @test MCInterval(0,3) ^ -2 == MCInterval(1/9, Inf)
    @test MCInterval(-3,0) ^ -2 == MCInterval(1/9, Inf)
    @test MCInterval(-3,2) ^ -2 == MCInterval(1/9, Inf)
    @test MCInterval(2,3) ^ -2 == MCInterval(1/9, 1/4)
    @test MCInterval(1,2) ^ -3 == MCInterval(1/8, 1.0)
    @test MCInterval(0,3) ^ -3 == MCInterval(1/27, Inf)
    @test MCInterval(-1,2) ^ -3 == entireMCinterval(Float64)
    @test MCInterval(2.5,2.5)^3 == MCInterval(15.625, 15.625)
end

@testset "Exp and log tests" begin
    @test exp(MCInterval(1.1051709180756477, 1.1051709180756477)) == MCInterval(3.019740552945523, 3.019740552945523)
    @test diam(exp(MCInterval(0.1,0.1))) == 0.0
    @test log(MCInterval(0.1,0.1)) == MCInterval(-2.3025850929940455, -2.3025850929940455)
    @test diam(log(MCInterval(0.1,0.1))) == 0.0
    @test exp2(MCInterval(1024.0,1024.0)) == MCInterval(Inf, Inf)
    @test exp10(MCInterval(308.5,308.5)) == MCInterval(Inf, Inf)
    @test log2(MCInterval(0.25, 0.5)) == MCInterval(-2.0, -1.0)
    @test log10(MCInterval(0.01, 0.1)) == MCInterval(log10(0.01), log10(0.1))
end

@testset "Min/Max" begin
xf = MCInterval(1.0,3.0)
xe = emptyMCinterval(Float64)
minf = min(xf,2.0)
mine = min(xe,2.0)
maxf = max(xf,2.0)
maxe = max(xe,2.0)

@test mine == xe
@test maxe == xe
@test minf == MCInterval(1.0,2.0)
@test maxf == MCInterval(2.0,3.0)

xf = MCInterval(1.0,3.0)
xe = emptyMCinterval(Float64)

minf = min(2.0,xf)
mine = min(2.0,xe)
maxf = max(2.0,xf)
maxe = max(2.0,xe)

@test mine == xe
@test maxe == xe
@test minf == MCInterval(1.0,2.0)
@test maxf == MCInterval(2.0,3.0)

@test inf(xf) == 1.0
@test sup(xf) == 3.0
@test real(xf) == xf
@test abs(xf) == xf
@test abs(-xf) == xf
end

#=
t1 = zero(a)
t2 = one(a)
t3 = zero(MCInterval{Float64})
t4 = one(MCInterval{Float64})
t5 = +a
t6 = -a
t7 = b-a
t8 = a+b
t9 = a*b
t10 = inv(a)
t11 = a/b
t12 = fma(a,a,b)
t13 = EAGOIntervalArithmetic.mag(a)
t14 = EAGOIntervalArithmetic.mig(a)
t15 = EAGOIntervalArithmetic.inf(a)
t16 = EAGOIntervalArithmetic.sup(a)
t17 = real(a)
t18 = abs(a)
t20 = min(a,b)
t21 = max(a,b)
t22 = EAGOIntervalArithmetic.dist(a,b)
t23 = eps(a)
t24 = floor(a)
t25 = ceil(a)
t26 = trunc(a)
t27 = sign(a)
t28 = EAGOIntervalArithmetic.mid(a,0.4)
t29 = EAGOIntervalArithmetic.mid(a)
t30 = EAGOIntervalArithmetic.diam(a)
t31 = EAGOIntervalArithmetic.radius(a)
=#

end
