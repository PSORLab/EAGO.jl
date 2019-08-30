#=
function MC_1_is_equal(y, x, tol)
    bool1 = isapprox(y.cc,x.cc,atol=tol)
    bool2 = isapprox(y.cv,x.cv,atol=tol)
    bool3 = isapprox(y.cv_grad[1], x.cv_grad[1], atol=tol)
    bool4 = isapprox(y.cc_grad[1], x.cc_grad[1], atol=tol)
    bool5 = isapprox(y.Intv.lo, x.Intv.lo, atol=tol)
    bool6 = isapprox(y.Intv.hi, x.Intv.hi, atol=tol)
    return (bool1 && bool2 && bool3 && bool4 && bool5 && bool6)
end
=#
#=
a = MC{1}(1.0,EAGO.IntervalType(0.4,3.0),1)
a1 = MC{1}(-7.0,EAGO.IntervalType(-12.0,-4.0),1)
b = MC{1}(EAGO.IntervalType(-10.0,-1.0))
c = MC{1}(2.0,EAGO.IntervalType(1.1,4.5),1)
aout1, bout1, cout1 = mul_rev(a,b,c)
aout2, bout2, cout2 = div_rev(a,b,c)
=#

#=
@testset "Reverse Multiplication" begin

    # THE BINARY OPERATOR
    a = MC{1}(1.0, EAGO.IntervalType(0.4,3.0), 1)
    b = MC{1}(EAGO.IntervalType(-10.0,-1.0))
    c = MC{1}(2.0, EAGO.IntervalType(1.1,4.5), 1)

    aout1, bout1, cout1 = mul_rev(a,b,c)

    @test bout1.Intv.lo == Inf
    @test bout1.Intv.hi == -Inf
    @test cout1.Intv.lo == Inf
    @test cout1.Intv.hi == -Inf

    bout1 = MC{1}(1.0, EAGO.IntervalType(0.4,3.0), 1)
    cout1 = MC{1}(EAGO.IntervalType(-10.0,-1.0))
    aout1 = bout1*cout1

    aout1_a, bout1_a, cout1_a = mul_rev(aout1, bout1, cout1)

    MC_1_is_equal(aout1_a, aout1, 0.00001)
    MC_1_is_equal(bout1_a, bout1, 0.00001)
    MC_1_is_equal(cout1_a, cout1, 0.00001)

    bout2 = MC{1}(1.0, EAGO.IntervalType(0.4,3.0), 1)
    cout2 = MC{1}(EAGO.IntervalType(-10.0,-1.0))
    aout2 = 0.3*bout1*cout1+1.0

    aout2_a, bout2_a, cout2_a = mul_rev(aout2, bout2, cout2)

    MC_1_is_equal(aout2_a, aout2, 0.00001)
    MC_1_is_equal(bout2_a, bout2, 0.00001)
    @test cout2_a.Intv.lo == -10.0
    @test cout2_a.Intv.hi == -1.0
    @test cout2_a.cv == -6.32
    @test cout2_a.cc == -1.0
    @test cout2_a.cv_grad[1] == -8.38
    @test cout2_a.cc_grad[1] == 0.0

    # WITH FLOAT

end
=#

#=
@testset "Reverse Addition" begin
end
=#

#=
@testset "Reverse Division" begin
end
=#

#=
@testset "Reverse Subtraction" begin
end
=#

#=
@testset "Reverse Exponential" begin
    a = MC{1}(1.0,EAGO.IntervalType(0.4,3.0),1)
    expa = exp(a)*1.1
    y,x = exp_rev(expa,a)

    @test MC_1_is_equal(expa, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test x.Intv.lo == 0.49531
    @test x.Intv.hi == 3.0
    @test x.cc_grad == 1.0
    @test x.cv_grad == 1.0

    exp2a = exp2(a)*1.1
    y,x = exp2_rev(exp2a,a)

    @test MC_1_is_equal(exp2a, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test x.Intv.lo == 0.53753
    @test x.Intv.hi == 3.0
    @test x.cc_grad == 1.0
    @test x.cv_grad == 1.0

    exp10a = exp10(a)*1.1
    y,x = exp10_rev(exp10a,a)

    @test MC_1_is_equal(exp10a, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test x.Intv.lo == 0.441392
    @test x.Intv.hi == 3.0
    @test x.cc_grad == 1.0
    @test x.cv_grad == 1.0

    expm1a = expm1(a)*1.1
    y,x = expm1_rev(expm1a,a)

    @test MC_1_is_equal(expm1a, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test x.Intv.lo == 0.432436
    @test x.Intv.hi == 3.0
    @test x.cc_grad == 1.0
    @test x.cv_grad == 1.0
end
=#

#=
@testset "Reverse Logarithm" begin

    a = MC{1}(9.6,EAGO.IntervalType(9.4,10.0),1)
    a1 = MC{1}(1.0,EAGO.IntervalType(0.2,5.0),1)
    loga = log(a)*5.1
    y,x = log_rev(loga,a)

    @test x.Intv.lo == Inf
    @test x.Intv.hi == -Inf

end
=#

#=
# BROKEN FLOAT REVERSE
bout1 = MC{1}(1.0,EAGO.IntervalType(0.4,3.0),1)
cout1 = -3.0
aout1 = bout1*cout1
aout1_a, bout1_a, cout1_a = mul_rev(aout1,bout1,cout1)
=#

a = MC{1}(9.6,EAGO.IntervalType(9.4,10.0),1)
a1 = MC{1}(1.0,EAGO.IntervalType(0.2,5.0),1)
loga = log(a)*5.1
y,x = log_rev(loga,a)

#a0 = MC{1}(7.0,EAGO.IntervalType(4.5,12.0),1)
#b0 = MC{1}(EAGO.IntervalType(6.0,9.0))
#c0 = MC{1}(5.0,EAGO.IntervalType(4.1,9.5),1)
#aout3, bout3, cout3 = pow_rev(a0,b0,c0)
