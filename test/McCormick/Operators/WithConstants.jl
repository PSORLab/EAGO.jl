module Test_WithConstant

using Compat
using Compat.Test
using EAGO
using IntervalArithmetic
using StaticArrays

EAGO.set_diff_relax(0)

@testset "Addition/Subtraction Constant" begin
    EAGO.set_diff_relax(0)
    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0);Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)

    X1 = X + 2.1
    X2 = 2.3 + X
    X3 = X + 2
    X4 = 3 + X
    @test +X == X
    @test X1.cv == 6.6
    @test X1.cc == 6.6
    @test X2.cv == 6.8
    @test X2.cc == 6.8
    @test X3.cv == 6.5
    @test X3.cc == 6.5
    @test X4.cv == 7.5
    @test X4.cc == 7.5

    X1 = X + Float16(2.1)
    X2 = Float16(2.3) + X
    X3 = X + Int16(2)
    X4 =  Int16(3) + X
    @test +X == X
    @test X1.cv == 6.599609375
    @test X1.cc == 6.599609375
    @test X2.cv == 6.80078125
    @test X2.cc == 6.80078125
    @test X3.cv == 6.5
    @test X3.cc == 6.5
    @test X4.cv == 7.5
    @test X4.cc == 7.5

    X1n = X - 2.1
    X2n = 2.3 - X
    X3n = X - 2
    X4n = 3 - X
    @test X1n.cv == 2.4
    @test X1n.cc == 2.4
    @test X2n.cv == -2.2
    @test X2n.cc == -2.2
    @test X3n.cv == 2.5
    @test X3n.cc == 2.5
    @test X4n.cv == -1.5
    @test X4n.cc == -1.5

    X1n = X - Float16(2.1)
    X2n = Float16(2.3) - X
    X3n = X - Int16(2)
    X4n = Int16(3) - X
    @test X1n.cv == 2.400390625
    @test X1n.cc == 2.400390625
    @test X2n.cv == -2.19921875
    @test X2n.cc == -2.19921875
    @test X3n.cv == 2.5
    @test X3n.cc == 2.5
    @test X4n.cv == -1.5
    @test X4n.cc == -1.5
end

@testset "Multiplication/Division Constant" begin

    EAGO.set_diff_relax(0)
    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0);Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)

    X1 = X * 2.1
    X2 = 2.3 * X
    X3 = X * 2
    X4 = 3 * X
    @test X1.cv == 9.450000000000001
    @test X1.cc == 9.450000000000001
    @test X2.cv == 10.35
    @test X2.cc == 10.35
    @test X3.cv == 9.0
    @test X3.cc == 9.0
    @test X4.cv == 13.5
    @test X4.cc == 13.5

    X1 = X * Float16(2.1)
    X2 = Float16(2.3) * X
    X3 = X * Int16(2)
    X4 =  Int16(3) * X
    @test X1.cv == 9.4482421875
    @test X1.cc == 9.4482421875
    @test X2.cv == 10.353515625
    @test X2.cc == 10.353515625
    @test X3.cv == 9.0
    @test X3.cc == 9.0
    @test X4.cv == 13.5
    @test X4.cc == 13.5

    X1 = X * (-2.1)
    X2 = (-2.3) * X
    X3 = X * (-2)
    X4 = (-3) * X
    @test X1.cv == -9.450000000000001
    @test X1.cc == -9.450000000000001
    @test X2.cv == -10.35
    @test X2.cc == -10.35
    @test X3.cv == -9.0
    @test X3.cc == -9.0
    @test X4.cv == -13.5
    @test X4.cc == -13.5

    X1 = X * Float16(-2.1)
    X2 = Float16(-2.3) * X
    X3 = X * Int16(-2)
    X4 =  Int16(-3) * X
    @test X1.cv == -9.4482421875
    @test X1.cc == -9.4482421875
    @test X2.cv == -10.353515625
    @test X2.cc == -10.353515625
    @test X3.cv == -9.0
    @test X3.cc == -9.0
    @test X4.cv == -13.5
    @test X4.cc == -13.5

    a = seed_g(Float64,1,2)
    b = seed_g(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,4.0);Interval(-5.0,-3.0)])
    mBox = mid.(xIBox)
    Y = SMCg{2,Interval{Float64},Float64}(-4.0,-4.0,b,b,xIBox[2],false)
    out = 1.0/Y
    out1 = 1/Y
    @test isapprox(out.cc,-0.25,atol=1E-6)
    @test isapprox(out.cv,-0.266666666,atol=1E-6)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-6)
    @test isapprox(out.cc_grad[2],-0.0625,atol=1E-6)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-6)
    @test isapprox(out.cv_grad[2],-0.0666667,atol=1E-6)
    @test isapprox(out.Intv.lo,-0.33333333,atol=1E-5)
    @test isapprox(out.Intv.hi,-0.199999,atol=1E-5)

    @test isapprox(out1.cc,-0.25,atol=1E-6)
    @test isapprox(out1.cv,-0.266666666,atol=1E-6)
    @test isapprox(out1.cc_grad[1],0.0,atol=1E-6)
    @test isapprox(out1.cc_grad[2],-0.0625,atol=1E-6)
    @test isapprox(out1.cv_grad[1],0.0,atol=1E-6)
    @test isapprox(out1.cv_grad[2],-0.0666667,atol=1E-6)
    @test isapprox(out1.Intv.lo,-0.33333333,atol=1E-5)
    @test isapprox(out1.Intv.hi,-0.199999,atol=1E-5)
end

@testset "Minimum/Maximum Constant" begin
    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)
    out = min(3,X)
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    out = min(X,3)
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)
    out = max(5,X)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)
    out = max(X,5)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)
    out = min(3.0,X)
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    out = min(X,3.0)
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)
    out = max(5.0,X)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)
    out = max(X,5.0)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)
    out = min(Float16(3.0),X)
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    out = min(X,Float16(3.0))
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)
    out = max(Float16(5.0),X)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)
    out = max(X,Float16(5.0))
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)
end

@testset "Conversion" begin
    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false)
    X1 = convert(SMCg{2,Interval{Float64},Float64},1)
    X2 = convert(SMCg{2,Interval{Float64},Float64},1.1)
    X3 = convert(SMCg{2,Interval{Float64},Float64},Interval(2.1,4.3))
    @test X1.cc == 1.0
    @test X1.cv == 1.0
    @test X2.cc == 1.1
    @test X2.cv == 1.1
    @test X3.cc == 4.3
    @test X3.cv == 2.1
end

end
