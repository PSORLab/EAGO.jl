module IntervalArithmetic_Check

using Compat
using Compat.Test
using EAGO

@testset "boolean check" begin
a = MCInterval(2.0,3.0)
b = MCInterval(4.0,6.0)

@test a != b
@test a <= b
@test a < b
@test ~(a >= b)
@test ~(a > b)
end

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

    @test (-MCInterval(-30.0,-15.0)) / MCInterval(-5.0,-3.0) == MCInterval(-10.0, -3.0)
    @test (-b)/a == MCInterval(-20.0, -0.8181818181818181)
    @test (-a)/c == MCInterval(-4.4, -0.025)
    @test (-c)/4.0 == MCInterval(-1.0, -0.0625)
    @test (-c)/zero(c) == emptyMCinterval(c)

    @test MCInterval(0.0, 1.0)/MCInterval(0.0,1.0) == MCInterval(0.0, Inf)
    @test MCInterval(-1.0, 1.0)/MCInterval(0.0,1.0) == entireMCinterval(c)
    @test MCInterval(-1.0, 1.0)/MCInterval(-1.0,1.0) == entireMCinterval(c)

    xz1 = MCInterval(-0.9, 2.0)
    xz2 = MCInterval(0.0, 2.0)
    xz3 = MCInterval(-0.9, 0.0)

    @test 4.0/xz1 == MCInterval(-Inf, Inf)
    @test 4.0/xz2 == MCInterval(2.0, Inf)
    @test 4.0/xz2 == MCInterval(2.0, Inf)
    @test -4.0/xz1 == MCInterval(-Inf, Inf)
    @test -4.0/xz2 == MCInterval(-Inf, -2.0)
    @test -4.0/xz2 == MCInterval(-Inf, -2.0)
    @test 0.0/xz1 == MCInterval(0.0, 0.0)
    @test -0.0/xz2 == MCInterval(0.0, 0.0)
    @test 0.0/xz2 == MCInterval(0.0, 0.0)

    x = MCInterval(1.0, 2.0)

    @test 0.1 + x == MCInterval(1.1, 2.1)
    @test 3.0 - x == MCInterval(1.0, 2.0)
    @test 3.1 - x == MCInterval(1.1, 2.1)
    @test 0.1 * MCInterval(1.0, 1.0) == MCInterval(0.1, 0.1)
    @test (-0.1) * MCInterval(1.0, 1.0) == MCInterval(-0.1, -0.1)
    @test (-0.1) * MCInterval(-2.0, 1.0) == MCInterval(-0.1, 0.2)
    @test MCInterval(1.0, 1.0) * 0.1 == MCInterval(0.1, 0.1)
    @test MCInterval(1.0, 1.0) * 1 == MCInterval(1.0, 1.0)
    @test 1 * MCInterval(1.0, 1.0) == MCInterval(1.0, 1.0)
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
    @test MCInterval{Float64}(2.5,2.5)^3 == MCInterval{Float64}(15.625, 15.625)
    @test MCInterval{Float64}(-3.5,-2.5)^-3 == MCInterval{Float64}(-0.064, -0.023323615160349854)

    @test pow(emptyMCinterval(Float64),2) == emptyMCinterval(Float64)
    @test pow(MCInterval{Float64}(-3,2),2) == EAGO.MCInterval{Float64}(0.0, 9.0)
    @test pow(MCInterval{Float64}(-3,2),3) == EAGO.MCInterval{Float64}(-27.0, 8.0)
    @test pow(emptyMCinterval(Float64),2.1) == emptyMCinterval(Float64)
    @test pow(MCInterval{Float64}(1,10),2.1) == EAGO.MCInterval{Float64}(1.0, 125.89254117941677)
    #=
    ADD CONVERSION METHOD TO COVER NEGATIVE POWERS
    @test MCInterval(0,3) ^ -2 == MCInterval(1/9, Inf)
    @test MCInterval(-3,0) ^ -2 == MCInterval(1/9, Inf)
    @test MCInterval(-3,2) ^ -2 == MCInterval(1/9, Inf)
    @test MCInterval(2,3) ^ -2 == MCInterval(1/9, 1/4)
    @test MCInterval(1,2) ^ -3 == MCInterval(1/8, 1.0)
    @test MCInterval(0,3) ^ -3 == MCInterval(1/27, Inf)
    @test MCInterval(-1,2) ^ -3 == entireMCinterval(Float64)
    =#

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

@test min(MCInterval(1.0,3.0),MCInterval(2.0,4.0)) == MCInterval(1.0, 3.0)
@test max(MCInterval(1.0,3.0),MCInterval(2.0,4.0)) == MCInterval(2.0,4.0)

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

@testset "Sqrt" begin
xf = MCInterval(1.0,3.0)
xe = emptyMCinterval(Float64)
sf = sqrt(xf)
se = sqrt(xe)
@test se == xe
@test sf == EAGO.MCInterval(1.0, 1.7320508075688772)
end


@testset "Utilities" begin

@test zero(MCInterval{Float64}) == MCInterval(zero(Float64))
@test one(zero(MCInterval{Float64})) == MCInterval(one(Float64))

@test floor(emptyMCinterval(Float64)) == emptyMCinterval(Float64)
@test floor(MCInterval(1.1,2.3)) == MCInterval(1.0, 2.0)

@test ceil(emptyMCinterval(Float64)) == emptyMCinterval(Float64)
@test ceil(MCInterval(1.1,2.3)) == MCInterval(2.0, 3.0)

@test trunc(emptyMCinterval(Float64)) == emptyMCinterval(Float64)
@test trunc(MCInterval(1.1,2.3)) == MCInterval(1.0, 2.0)

@test step(emptyMCinterval(Float64)) == emptyMCinterval(Float64)
@test step(MCInterval(1.1,2.3)) == MCInterval(1.0, 1.0)

@test sign(emptyMCinterval(Float64)) == emptyMCinterval(Float64)
@test sign(MCInterval(-1.1,2.3)) == MCInterval(-1.0, 1.0)

@test dist(MCInterval(-1.1,2.3),MCInterval(-1.1,2.3)) == 0.0
@test eps(MCInterval(-1.1,2.3)) == 4.440892098500626e-16
@test radius(MCInterval(-1.1,2.3)) == 1.7
@test mid(MCInterval(-1.1,2.3),0.6) == 0.94
end

end
