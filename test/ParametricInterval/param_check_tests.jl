module ParamChk_Tests

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

# Tests strict inclusion procedure for interval vectors
@testset "Strictly In" begin
    Y = [IntervalType(0,5),IntervalType(0,5),IntervalType(0,5)]
    X1 = [IntervalType(1,2),IntervalType(1,2),IntervalType(1,2)]
    X2 = [IntervalType(1,2),IntervalType(-10,10),IntervalType(1,2)]
    X3 = [IntervalType(1,2),IntervalType(1,2),IntervalType(0,5)]
    flag1 = EAGO.strict_x_in_y(X1,Y)
    flag2 = EAGO.strict_x_in_y(X2,Y)
    flag3 = EAGO.strict_x_in_y(X3,Y)
    @test flag1 == true
    @test flag2 == false
    @test flag3 == false
end

# Checks extended division routine
@testset "Extended Division" begin
    A1 = IntervalType(0)
    A2 = IntervalType(0,3)
    A3 = IntervalType(-2,0)
    A4 = IntervalType(-3,2)
    ind1,B1,C1 = EAGO.extended_divide(A1)
    ind2,B2,C2 = EAGO.extended_divide(A2)
    ind3,B3,C3 = EAGO.extended_divide(A3)
    ind4,B4,C4 = EAGO.extended_divide(A4)
    @test ind1 == 0
    @test ind2 == 1
    @test ind3 == 2
    @test ind4 == 3
    @test B1 == IntervalType(-Inf,Inf)
    @test 0.33333 - 1E-4 <= B2.lo <= 0.33333 + 1E-4
    @test B2.hi == Inf
    @test B3 == IntervalType(-Inf,-0.5)
    @test B4.lo == -Inf
    @test -0.33333 - 1E-4 <= B4.hi <= -0.33333 + 1E-4
    @test C1 == IntervalType(-Inf,Inf)
    @test C2 == IntervalType(Inf,Inf)
    @test C3 == IntervalType(-Inf,-Inf)
    @test C4 == IntervalType(0.5,Inf)
end

@testset "Extended Process" begin
    N =  IntervalType(-5,5)
    X = IntervalType(-5,5)
    Mii = IntervalType(-5,5)
    S1 = IntervalType(-5,5)
    S2 = IntervalType(-5,5)
    B = IntervalType(-5,5)
    rtol = 1E-4
    indx1,box11,box12 = EAGO.extended_process(N,X,Mii,S1,S2,B,rtol)
    Miib = IntervalType(0,5)
    S1b = IntervalType(1,5)
    S2b = IntervalType(1,5)
    Bb = IntervalType(1,5)
    indx2,box21,box22 = EAGO.extended_process(N,X,Miib,S1b,S2b,Bb,rtol)
    Miic = IntervalType(-5,0)
    S1c = IntervalType(1,5)
    S2c = IntervalType(1,5)
    Bc = IntervalType(1,5)
    indx3,box31,box32 = EAGO.extended_process(N,X,Miic,S1c,S2c,Bc,rtol)
    Miia = IntervalType(1,5)
    S1a = IntervalType(1,5)
    S2a = IntervalType(1,5)
    Ba = IntervalType(1,5)
    indx6,box61,box62 = EAGO.extended_process(N,X,Miia,S1a,S2a,Ba,rtol)
    Miid = IntervalType(0,0)
    S1d = IntervalType(1,5)
    S2d = IntervalType(1,5)
    Bd = IntervalType(1,5)
    indx8,box81,box82 = EAGO.extended_process(N,X,Miid,S1d,S2d,Bd,rtol)

    @test indx1 == 0
    @test box11 == IntervalType(-Inf,Inf)
    @test box12 == IntervalType(-5,5)

    @test indx2 == 0
    @test box21.hi > -Inf
    @test box22 == IntervalType(-5,5)

    @test indx3 == 0
    @test box31.lo < Inf
    @test box32 == IntervalType(-5,5)

    @test indx6 == 1
    @test -15.0002 - 1E-4 <= box61.lo <= -15.0002 + 1E-4
    @test box62.lo == -Inf
    @test box61.hi == Inf
    @test -0.599979 - 1E-4 <= box62.hi <= -0.599979 + 1E-4

    @test indx8 == 0
    @test box81 == IntervalType(-5,5)
    @test box82 == IntervalType(-5,5)
end

end
