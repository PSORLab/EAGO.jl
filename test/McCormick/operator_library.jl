set_mc_differentiability!# DONE
@testset "Test Exponentials" begin

   mctol = 1E-4
   m = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0),
             seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   y = McCormick.exp(m)
   yref = MC{2}(15.154262241479262, 37.30486063158251, EAGO.IntervalType(2.71828, 54.5982), SVector{2,Float64}([0.0, 0.0]), SVector{2,Float64}([17.2933, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = McCormick.exp2(m)
   yref = MC{2}(4.0, 11.33333333333333, EAGO.IntervalType(1.99999, 16.0001), SVector{2,Float64}([2.77259, 0.0]), SVector{2,Float64}([4.66667, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = McCormick.exp10(m)
   yref = MC{2}(100.0, 6670.0000000000055, EAGO.IntervalType(9.999999999999999999, 10000.00000000001), SVector{2,Float64}([2302.5850929940457, 0.0]), SVector{2,Float64}([3330.0, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = McCormick.expm1(m)
   yref = MC{2}(6.38905609893065, 36.304860631582514, EAGO.IntervalType(1.71828, 53.5982), SVector{2,Float64}([7.38906, 0.0]), SVector{2,Float64}([17.2933, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)
end

# DONE
@testset "Test Logarithms" begin

   mctol = 1E-4
   m = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0),
             seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   y = log(m)
   yref = MC{2}(0.46209812037329695, 0.6931471805599453, EAGO.IntervalType(0, 1.3863), SVector{2,Float64}([0.462098, 0.0]), SVector{2,Float64}([0.5, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = log2(m)
   yref = MC{2}(0.6666666666666666, 1.0, EAGO.IntervalType(0, 2), SVector{2,Float64}([0.666667, 0.0]), SVector{2,Float64}([0.721348, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = log10(m)
   yref = MC{2}(0.20068666377598746, 0.3010299956639812, EAGO.IntervalType(0, 0.60206), SVector{2,Float64}([0.200687, 0.0]), SVector{2,Float64}([0.217147, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = log1p(m)
   yref = MC{2}(0.998577424517997, 1.0986122886681098, EAGO.IntervalType(0.693147, 1.60944), SVector{2,Float64}([0.30543, 0.0]), SVector{2,Float64}([0.333333, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)
end

# DONE
@testset "Addition/Subtraction Constant" begin
    EAGO.set_mc_differentiability!(0)
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,EAGO.IntervalType}([EAGO.IntervalType(-3.0,8.0);EAGO.IntervalType(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)

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

# DOUBLE CHECK DIVISION
@testset "Multiplication/Division Constant" begin

    EAGO.set_mc_differentiability!(0)
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0);Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)

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

    #=
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,4.0);Interval(-5.0,-3.0)])
    Y = MC{2}(-4.0,-4.0,xIBox[2],b,b,false)
    out = 1.0/Y
    out1 = 1/Y
    @test isapprox(out.cc,-0.25,atol=1E-6)                    # FAIL
    @test isapprox(out.cv,-0.266666666,atol=1E-6)             # FAIL
    @test isapprox(out.cc_grad[1],0.0,atol=1E-6)
    @test isapprox(out.cc_grad[2],-0.0625,atol=1E-6)          # FAIL
    @test isapprox(out.cv_grad[1],0.0,atol=1E-6)
    @test isapprox(out.cv_grad[2],-0.0666667,atol=1E-6)       # FAIL
    @test isapprox(out.Intv.lo,-0.33333333,atol=1E-5)         # FAIL
    @test isapprox(out.Intv.hi,-0.199999,atol=1E-5)           # FAIL

    @test isapprox(out1.cc,-0.25,atol=1E-6)                   # FAIL
    @test isapprox(out1.cv,-0.266666666,atol=1E-6)            # FAIL
    @test isapprox(out1.cc_grad[1],0.0,atol=1E-6)
    @test isapprox(out1.cc_grad[2],-0.0625,atol=1E-6)         # FAIL
    @test isapprox(out1.cv_grad[1],0.0,atol=1E-6)
    @test isapprox(out1.cv_grad[2],-0.0666667,atol=1E-6)      # FAIL
    @test isapprox(out1.Intv.lo,-0.33333333,atol=1E-5)        # FAIL
    @test isapprox(out1.Intv.hi,-0.199999,atol=1E-5)          # FAIL
    =#
end

# DONE
@testset "Minimum/Maximum Constant" begin
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
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

    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = max(5,X)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = max(X,5)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
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

    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = max(5.0,X)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = max(X,5.0)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
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

    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = max(Float16(5.0),X)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
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

# DONE
@testset "Conversion" begin
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    X1 = convert(MC{2},1)
    X2 = convert(MC{2},1.1)
    X3 = convert(MC{2},Interval(2.1,4.3))
    @test X1.cc == 1.0
    @test X1.cv == 1.0
    @test X2.cc == 1.1
    @test X2.cv == 1.1
    @test X3.cc == 4.3
    @test X3.cv == 2.1
end

# DOUBLE CHECK A CASE
@testset "Multiplication Operator" begin


    ##############      Testing for Nonsmooth Standard Mult           ##############
    EAGO.set_mc_differentiability!(0)

    ################### Test Nonsmooth Zero in Both Case (Failing)   ###############
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-2.0,1.0);Interval(-1.0,2.0)])
    X = MC{2}(0.0,0.0,xIBox[1],a,a,false)
    Y = MC{2}(1.0,1.0,xIBox[2],b,b,false)
    out = X*Y

    @test out.cc == 2.0
    @test out.cv == -1.0
    @test out.cc_grad[1] == 2.0
    @test out.cc_grad[2] == -2.0
    @test out.cv_grad[1] == 2.0
    @test out.cv_grad[2] == 1.0
    @test out.Intv.lo == -4.0
    @test out.Intv.hi == 2.0

    ###################### Test Nonsmooth X1.l>0   (Passing)  ######################
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(1.0,5.0);Interval(-1.0,2.0)])
    X = MC{2}(3.0,3.0,xIBox[1],a,a,false)
    Y = MC{2}(1.0,1.0,xIBox[2],b,b,false)
    out = X*Y

    @test out.cc == 5.0
    @test out.cv == 1.0
    @test out.cc_grad[1] == 2.0
    @test out.cc_grad[2] == 1.0
    @test out.cv_grad[1] == 2.0
    @test out.cv_grad[2] == 5.0
    @test out.Intv.lo == -5.0
    @test out.Intv.hi == 10.0

    ############## Test Nonsmooth X1.h<0  &&  X2.l>0 (Passing)  ######################
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-6.0,-2.0);Interval(1.0,3.0)])
    X = MC{2}(-4.0,-4.0,xIBox[1],a,a,false)
    Y = MC{2}(2.0,2.0,xIBox[2],b,b,false)
    out = X*Y

    @test out.cc == -6.0
    @test out.cv == -10.0
    @test out.cc_grad[1] == 1.0
    @test out.cc_grad[2] == -2.0
    @test out.cv_grad[1] == 3.0
    @test out.cv_grad[2] == -2.0
    @test out.Intv.lo == -18.0
    @test out.Intv.hi == -2.0

    ############## Test Nonsmooth X1.h<0  &&  X2.h<0 (Passing)  ######################
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-6.0,-2.0);Interval(-7.0,-3.0)])
    X = MC{2}(-4.0,-4.0,xIBox[1],a,a,false)
    Y = MC{2}(-5.0,-5.0,xIBox[2],b,b,false)
    out = X*Y

    @test out.cc == 24.0
    @test out.cv == 16.0
    @test out.cc_grad[1] == -7.0
    @test out.cc_grad[2] == -2.0
    @test out.cv_grad[1] == -3.0
    @test out.cv_grad[2] == -2.0
    @test out.Intv.lo == 6.0
    @test out.Intv.hi == 42.0

    ############## Test Nonsmooth X1.h<0  &&  0 in X2 (Passing)  ###################
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-6.0,-2.0);Interval(-7.0,4.0)])
    X = MC{2}(-4.0,-4.0,xIBox[1],a,a,false)
    Y = MC{2}(-5.0,-5.0,xIBox[2],b,b,false)
    out = X*Y

    @test out.cc == 24.0
    @test out.cv == 16.0
    @test out.cc_grad[1] == -7.0
    @test out.cc_grad[2] == -2.0
    @test out.cv_grad[1] == -7.0
    @test out.cv_grad[2] == -6.0
    @test out.Intv.lo == -24.0
    @test out.Intv.hi == 42.0

    ############## Test Nonsmooth 0 in X1  &&  X2.l > 0 ()  #################
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,4.0);Interval(1.0,4.0)])
    X = MC{2}(-2.0,-2.0,xIBox[1],a,a,false)
    Y = MC{2}(3.0,3.0,xIBox[2],b,b,false)
    out = X*Y

    @test out.cc == -5.0
    @test out.cv == -8.0
    @test out.cc_grad[1] == 4.0
    @test out.cc_grad[2] == -3.0
    @test out.cv_grad[1] == 1.0
    @test out.cv_grad[2] == -3.0
    @test out.Intv.lo == -12.0
    @test out.Intv.hi == 16.0

    ############## Test Nonsmooth 0 in X1  &&  X2.h < 0 ()         #################
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,4.0);Interval(-5.0,-3.0)])
    X = MC{2}(-2.0,-2.0,xIBox[1],a,a,false)
    Y = MC{2}(-4.0,-4.0,xIBox[2],b,b,false)
    out = X*Y

    @test out.cc == 9.0
    @test out.cv == 7.0
    @test out.cc_grad[1] == -3.0
    @test out.cc_grad[2] == -3.0
    @test out.cv_grad[1] == -5.0
    @test out.cv_grad[2] == -3.0
    @test out.Intv.lo == -20.0
    @test out.Intv.hi == 15.0

    ##############      Testing for Smooth Standard Mult           ##############
    EAGO.set_mc_differentiability!(1)

    seed1 = seed_gradient(Float64,1,2)
    seed2 = seed_gradient(Float64,2,2)
    x1 = MC{2}(0.0,0.0,Interval(-200.0,200.0),seed1,seed1,false)
    y1 = MC{2}(200.0,200.0,Interval(0.0,400.0),seed2,seed2,false)
    z1 = x1*y1
    @test isapprox(z1.cc,40000,atol=1E-4)
    @test isapprox(z1.cv,-40000,atol=1E-4)

    x2 = MC{2}(170.0,170.0,Interval(100.0,240.0),seed1,seed1,false)
    y2 = MC{2}(250.0,250.0,Interval(100.0,400.0),seed2,seed2,false)
    z2 = x2*y2
    @test isapprox(z2.cc,53000,atol=1E-4)
    @test isapprox(z2.cv,32000,atol=1E-4)

    x3 = MC{2}(-200.0,-200.0,Interval(-300.0,-100.0),seed1,seed1,false)
    y3 = MC{2}(-300.0,-300.0,Interval(-400.0,-200.0),seed2,seed2,false)
    z3 = x3*y3
    @test isapprox(z3.cc,70000,atol=1E-4)
    @test isapprox(z3.cv,50000,atol=1E-4)

    # CHECK ME AGAIN???? -47187.5 new, -47460.9375 old
    x4 = MC{2}(150.0,150.0,Interval(100.0,200.0),seed1,seed1,false)
    y4 = MC{2}(-250.0,-250.0,Interval(-500.0,-100.0),seed2,seed2,false)
    z4 = x4*y4
    @test isapprox(z4.cc,-30000,atol=1E-3)
    @test isapprox(z4.cv,-47187.5,atol=1E-3)

    x5 = MC{2}(-150.0,-150.0,Interval(-200.0,-100.0),seed1,seed1,false)
    y5 = MC{2}(300.0,300.0,Interval(200.0,400.0),seed2,seed2,false)
    z5 = x5*y5
    @test isapprox(z5.cv,-50000,atol=1E-4)
    @test isapprox(z5.cc,-40000,atol=1E-4)
end

# MOSTLY DONE!
@testset "Power" begin

    # tests powers (square)
    println("test 1")
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0); Interval(3.0,7.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    out = X^2
    @test isapprox(out.cc,19,atol=1E-8)
    @test isapprox(out.cv,16,atol=1E-8)
    @test isapprox(out.cc_grad[1],10.0,atol=1E-1)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],8.0,atol=1E-5)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out.Intv.lo,9,atol=1E-4)
    @test isapprox(out.Intv.hi,49,atol=1E-4)

    # tests powers (^-2 on positive domain)
    println("test 2")
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,7.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    out = X^(-2)
    @test isapprox(out.cc,0.08843537414965986,atol=1E-8)
    @test isapprox(out.cv,0.0625,atol=1E-8)
    @test isapprox(out.cc_grad[1],-0.0226757,atol=1E-1)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[1],-0.03125,atol=1E-5)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-4)
    @test isapprox(out.Intv.lo,0.0204081,atol=1E-4)
    @test isapprox(out.Intv.hi,0.111112,atol=1E-4)

    # tests powers (^-2 on negative domain)
    println("test 3")
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-8.0,-3.0);Interval(-8.0,-3.0)])
    X = MC{2}(-4.5,-4.5,xIBox[1],a,a,false)
    out = X^(-2)
    @test isapprox(out.cc,0.08246527777777776,atol=1E-8)
    @test isapprox(out.cv,0.04938271604938271,atol=1E-8)
    @test isapprox(out.cc_grad[1],0.0190972,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0219479,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,0.015625,atol=1E-4)
    @test isapprox(out.Intv.hi,0.111112,atol=1E-4)

    # tests powers (^1)
    println("test 4")
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-8.0,-3.0),Interval(-8.0,-3.0)])
    X = MC{2}(-4.5,-4.5,xIBox[1],a,a,false)
    out = X^(1)
    @test isapprox(out.cc,-4.5,atol=1E-8)
    @test isapprox(out.cv,-4.5,atol=1E-8)
    @test isapprox(out.cc_grad[1],1.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],1.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-8.0,atol=1E-4)
    @test isapprox(out.Intv.hi,-3.0,atol=1E-4)

    # tests powers (^2)
    println("test 5")
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-8.0,-3.0),Interval(-8.0,-3.0)])
    X = MC{2}(-4.5,-4.5,xIBox[1],a,a,false)
    out = X^(2)
    @test isapprox(out.cc,25.5,atol=1E-8)
    @test isapprox(out.cv,20.25,atol=1E-8)
    @test isapprox(out.cc_grad[1],-11.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],-9.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,9.0,atol=1E-4)
    @test isapprox(out.Intv.hi,64.0,atol=1E-4)

    # tests powers (^3)
    println("test 6")
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-8.0,-3.0);Interval(-8.0,-3.0)])
    X = MC{2}(-4.5,-4.5,xIBox[1],a,a,false)
    out = X^(3)
    @test isapprox(out.cc,-91.125,atol=1E-8)
    @test isapprox(out.cv,-172.5,atol=1E-8)
    @test isapprox(out.cc_grad[1],60.75,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],97.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-512,atol=1E-4)
    @test isapprox(out.Intv.hi,-27,atol=1E-4)

    # tests powers (^4)
    println("test 7")
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-8.0,-3.0),Interval(-8.0,-3.0)])
    X = MC{2}(-4.5,-4.5,xIBox[1],a,a,false)
    out = X^(4)
    @test isapprox(out.cc,1285.5,atol=1E-8)
    @test isapprox(out.cv,410.0625,atol=1E-8)
    @test isapprox(out.cc_grad[1],-803.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],-364.5,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,81,atol=1E-4)
    @test isapprox(out.Intv.hi,4096,atol=1E-4)

    # tests powers (^3 greater than zero ISSUE WITH CC)
    println("test 8")
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,8.0),Interval(3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = X^(3)
    @test isapprox(out.cc,172.5,atol=1E-8)
    @test isapprox(out.cv,91.125,atol=1E-8)
    @test isapprox(out.cc_grad[1],97.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],60.75,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,27,atol=1E-4)
    @test isapprox(out.Intv.hi,512,atol=1E-4)

    # tests powers (^4 greater than zero)
    println("test 9")
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,8.0),Interval(3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = X^(4)
    @test isapprox(out.cc,1285.5,atol=1E-1)
    @test isapprox(out.cv,410.0625,atol=1E-1)
    @test isapprox(out.cc_grad[1],803.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],364.5,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,81,atol=1E-4)
    @test isapprox(out.Intv.hi,4096,atol=1E-4)

    # tests powers (^4 zero in range)
    println("test 10")
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-8.0,3.0),Interval(-8.0,3.0)])
    X = MC{2}(-4.5,-4.5,xIBox[1],a,a,false)
    out = X^(4)
    @test isapprox(out.cc,2818.5,atol=1)
    @test isapprox(out.cv,410.0625,atol=1)
    @test isapprox(out.cc_grad[1],-365.0,atol=1)
    @test isapprox(out.cc_grad[2],0.0,atol=1)
    @test isapprox(out.cv_grad[1],-364.5,atol=1)
    @test isapprox(out.cv_grad[2],0.0,atol=1)
    @test isapprox(out.Intv.lo,0,atol=1)
    @test isapprox(out.Intv.hi,4096,atol=1)

    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    Y = MC{2}(7.0,7.0,xIBox[2],b,b,false)
    Xn = MC{2}(-4.0,-4.0,-xIBox[1],a,a,false)
    Xz = MC{2}(-2.0,-2.0,Interval(-3.0,1.0),a,a,false)

    println("test 11")
    out14 = inv(X)
    @test isapprox(out14.cc,0.2857142857142857,atol=1E-5)
    @test isapprox(out14.cv,0.25,atol=1E-5)
    @test isapprox(out14.cc_grad[1],-0.047619,atol=1E-5)
    @test isapprox(out14.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out14.cv_grad[1],-0.0625,atol=1E-5)
    @test isapprox(out14.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out14.Intv.lo,0.142857,atol=1E-5)
    @test isapprox(out14.Intv.hi,0.333334,atol=1E-5)
    out14a = inv(Xn)

    println("test 12")
    out23 = pow(X,2)
    @test isapprox(out23.cc,19.0,atol=1E-5)
    @test isapprox(out23.cv,16.0,atol=1E-5)
    @test isapprox(out23.cc_grad[1],10.0,atol=1E-5)
    @test isapprox(out23.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out23.cv_grad[1],8.0,atol=1E-5)
    @test isapprox(out23.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out23.Intv.lo,9.0,atol=1E-5)
    @test isapprox(out23.Intv.hi,49.0,atol=1E-5)

    println("test 13")
    out23a = pow(Xn,2)
    @test isapprox(out23a.cc,19.0,atol=1E-5)
    @test isapprox(out23a.cv,16.0,atol=1E-5)
    @test isapprox(out23a.cc_grad[1],-10.0,atol=1E-5)
    @test isapprox(out23a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out23a.cv_grad[1],-8.0,atol=1E-5)
    @test isapprox(out23a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out23a.Intv.lo,9.0,atol=1E-5)
    @test isapprox(out23a.Intv.hi,49.0,atol=1E-5)

    println("test 14")
    out23b = pow(Xz,2)
    @test isapprox(out23b.cc,7.0,atol=1E-5)
    @test isapprox(out23b.cv,2.66666666666666,atol=1E-5) #  double check me 4.0
    @test isapprox(out23b.cc_grad[1],-2.0,atol=1E-5)
    @test isapprox(out23b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out23b.cv_grad[1],-4.0,atol=1E-5) # double check me -4.0
    @test isapprox(out23b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out23b.Intv.lo,0.0,atol=1E-5)
    @test isapprox(out23b.Intv.hi,9.0,atol=1E-5)

    println("test 15")
    out1a = pow(X,3)
    @test isapprox(out1a.cc,106.0,atol=1E-5)
    @test isapprox(out1a.cv,64.0,atol=1E-5)
    @test isapprox(out1a.cc_grad[1],79.0,atol=1E-5)
    @test isapprox(out1a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out1a.cv_grad[1],48.0,atol=1E-5)
    @test isapprox(out1a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out1a.Intv.lo,27.0,atol=1E-5)
    @test isapprox(out1a.Intv.hi,343.0,atol=1E-5)

    println("test 16")
    out1b = pow(Xn,3)
    @test isapprox(out1b.cc,-64.0,atol=1E-5)
    @test isapprox(out1b.cv,-106.0,atol=1E-5)
    @test isapprox(out1b.cc_grad[1],48.0,atol=1E-5)
    @test isapprox(out1b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out1b.cv_grad[1],79,atol=1E-5)
    @test isapprox(out1b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out1b.Intv.lo,-343.0,atol=1E-5)
    @test isapprox(out1b.Intv.hi,-27.0,atol=1E-5)

    println("test 17")
    out1c = pow(Xz,3)
    @test isapprox(out1c.cc,-7.75,atol=1E-5)
    @test isapprox(out1c.cv,-20.25,atol=1E-5)
    @test isapprox(out1c.cc_grad[1],12.25,atol=1E-5)
    @test isapprox(out1c.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out1c.cv_grad[1],6.75,atol=1E-5)
    @test isapprox(out1c.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out1c.Intv.lo,-27.0,atol=1E-5)
    @test isapprox(out1c.Intv.hi,1.0,atol=1E-5)

    println("test 18")
    out2a = pow(X,-3)
    @test isapprox(out2a.cc,0.02850664075153871,atol=1E-5)
    @test isapprox(out2a.cv,0.015625,atol=1E-5)
    @test isapprox(out2a.cc_grad[1],-0.0085304,atol=1E-5)
    @test isapprox(out2a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out2a.cv_grad[1],-0.0117188,atol=1E-5)
    @test isapprox(out2a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out2a.Intv.lo,0.00291545,atol=1E-5)
    @test isapprox(out2a.Intv.hi,0.0370371,atol=1E-5)

    println("test 19")
    out2b = pow(Xn,-3)
    @test isapprox(out2b.cc,-0.015625,atol=1E-5)
    @test isapprox(out2b.cv,-0.02850664075153871,atol=1E-5)
    @test isapprox(out2b.cc_grad[1],-0.0117188,atol=1E-5)
    @test isapprox(out2b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out2b.cv_grad[1],-0.0085304,atol=1E-5)
    @test isapprox(out2b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out2b.Intv.lo,-0.0370371,atol=1E-5)
    @test isapprox(out2b.Intv.hi,-0.00291545,atol=1E-5)

    println("test 20")
    out3a = pow(X,-4)
    @test isapprox(out3a.cc,0.009363382541225106,atol=1E-5)
    @test isapprox(out3a.cv,0.00390625,atol=1E-5)
    @test isapprox(out3a.cc_grad[1],-0.0029823,atol=1E-5)
    @test isapprox(out3a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out3a.cv_grad[1],-0.00390625,atol=1E-5)
    @test isapprox(out3a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out3a.Intv.lo,0.000416493,atol=1E-5)
    @test isapprox(out3a.Intv.hi,0.0123457,atol=1E-5)

    println("test 21")
    out3b = pow(Xn,-4)
    @test isapprox(out3b.cc,0.009363382541225106,atol=1E-5)
    @test isapprox(out3b.cv,0.00390625,atol=1E-5)
    @test isapprox(out3b.cc_grad[1],0.0029823,atol=1E-5)
    @test isapprox(out3b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out3b.cv_grad[1],0.00390625,atol=1E-5)
    @test isapprox(out3b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out3b.Intv.lo,0.000416493,atol=1E-5)
    @test isapprox(out3b.Intv.hi,0.0123457,atol=1E-5)

    println("test 22")
    out4 = pow(X,4)
    @test isapprox(out4.cc,661.0,atol=1E-5)
    @test isapprox(out4.cv,256.0,atol=1E-5)
    @test isapprox(out4.cc_grad[1],580.0,atol=1E-5)
    @test isapprox(out4.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out4.cv_grad[1],256,atol=1E-5)
    @test isapprox(out4.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out4.Intv.lo,81.0,atol=1E-5)
    @test isapprox(out4.Intv.hi,2401.0,atol=1E-5)

    println("test 23")
    out4a = pow(Xn,4)
    @test isapprox(out4a.cc,661.0,atol=1E-5)
    @test isapprox(out4a.cv,256.0,atol=1E-5)
    @test isapprox(out4a.cc_grad[1],-580.0,atol=1E-5)
    @test isapprox(out4a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out4a.cv_grad[1],-256,atol=1E-5)
    @test isapprox(out4a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out4a.Intv.lo,81.0,atol=1E-5)
    @test isapprox(out4a.Intv.hi,2401.0,atol=1E-5)

    println("test 24")
    out4b = pow(Xz,4)
    @test isapprox(out4b.cc,61.0,atol=1E-5)
    @test isapprox(out4b.cv,16.0,atol=1E-5)
    @test isapprox(out4b.cc_grad[1],-20.0,atol=1E-5)
    @test isapprox(out4b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out4b.cv_grad[1],-32.0,atol=1E-5)
    @test isapprox(out4b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out4b.Intv.lo,0.0,atol=1E-5)
    @test isapprox(out4b.Intv.hi,81.0,atol=1E-5)
end

#=
# MOSTLY DONE! (NEED TO FIX THIS)
@testset "Square Root" begin
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,9.0),Interval(3.0,9.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)

    out = sqrt(X)
    @test isapprox(out.cc,2.1213203435596424,atol=1E-8)
    @test isapprox(out.cv,2.049038105676658,atol=1E-8)
    @test isapprox(out.cc_grad[1],0.235702,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.211325,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,1.73205,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)

    out1 = sqrt(X)
    @test isapprox(out1.cc,2.0,atol=1E-5)
    @test isapprox(out1.cv,1.9604759334428057,atol=1E-5)
    @test isapprox(out1.cc_grad[1],0.25,atol=1E-5)
    @test isapprox(out1.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out1.cv_grad[1],0.228425,atol=1E-5)
    @test isapprox(out1.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out1.Intv.lo,1.73205,atol=1E-5)
    @test isapprox(out1.Intv.hi,2.64576,atol=1E-5)
end
=#

# MOSTLY DONE!
@testset "Division" begin
    EAGO.set_mc_differentiability!(0)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,4.0);Interval(-5.0,-3.0)])
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    X = MC{2}(-2.0,-2.0,xIBox[1],a,a,false)
    Y = MC{2}(-4.0,-4.0,xIBox[2],b,b,false)
    out = X/Y
    @test isapprox(out.cc,0.6,atol=1E-6)
    @test isapprox(out.cv,0.41666666,atol=1E-6)
    @test isapprox(out.cc_grad[1],-0.2,atol=1E-6)
    @test isapprox(out.cc_grad[2],0.2,atol=1E-6)
    @test isapprox(out.cv_grad[1],-0.333333,atol=1E-6)
    @test isapprox(out.cv_grad[2],0.1875,atol=1E-6)
    @test isapprox(out.Intv.lo,-1.33333333,atol=1E-6)
    @test isapprox(out.Intv.hi,1.0,atol=1E-6)
end

# MOSTLY DONE!
@testset "Step Function" begin
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    Xn = MC{2}(-4.0,-4.0,-xIBox[1],a,a,false)
    Xz = MC{2}(-2.0,-2.0,Interval(-3.0,1.0),a,a,false)
    Xzp = MC{2}(0.5,0.5,Interval(-3.0,1.0),a,a,false)

    EAGO.set_mc_differentiability!(0)

    out21 = step(X)
    @test isapprox(out21.cc,1.0,atol=1E-5)
    @test isapprox(out21.cv,1.0,atol=1E-5)
    @test isapprox(out21.cc_grad[1],0.0,atol=1E-5)
    @test isapprox(out21.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out21.cv_grad[1],0.0,atol=1E-5)
    @test isapprox(out21.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out21.Intv.lo,1.0,atol=1E-5)
    @test isapprox(out21.Intv.hi,1.0,atol=1E-5)

    out21a = step(Xn)
    @test isapprox(out21a.cc,0.0,atol=1E-5)
    @test isapprox(out21a.cv,0.0,atol=1E-5)
    @test isapprox(out21a.cc_grad[1],0.0,atol=1E-5)
    @test isapprox(out21a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out21a.cv_grad[1],0.0,atol=1E-5)
    @test isapprox(out21a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out21a.Intv.lo,0.0,atol=1E-5)
    @test isapprox(out21a.Intv.hi,0.0,atol=1E-5)

    out21b = step(Xz)
    @test isapprox(out21b.cc,0.3333333333333333,atol=1E-5)
    @test isapprox(out21b.cv,0.0,atol=1E-5)
    @test isapprox(out21b.cc_grad[1],-0.6666666666666666,atol=1E-5) #0.3333333333333333
    @test isapprox(out21b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out21b.cv_grad[1],0.0,atol=1E-5)
    @test isapprox(out21b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out21b.Intv.lo,0.0,atol=1E-5)
    @test isapprox(out21b.Intv.hi,1.0,atol=1E-5)

    out21b = step(Xzp)
    @test isapprox(out21b.cc,1.0,atol=1E-5)
    @test isapprox(out21b.cv,0.5,atol=1E-5)
    @test isapprox(out21b.cc_grad[1],0.0,atol=1E-5)
    @test isapprox(out21b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out21b.cv_grad[1],1.0,atol=1E-5)
    @test isapprox(out21b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out21b.Intv.lo,0.0,atol=1E-5)
    @test isapprox(out21b.Intv.hi,1.0,atol=1E-5)

    EAGO.set_mc_differentiability!(1)

    out21 = step(X)
    @test isapprox(out21.cc,1.0,atol=1E-5)
    @test isapprox(out21.cv,1.0,atol=1E-5)
    @test isapprox(out21.cc_grad[1],0.0,atol=1E-5)
    @test isapprox(out21.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out21.cv_grad[1],0.0,atol=1E-5)
    @test isapprox(out21.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out21.Intv.lo,1.0,atol=1E-5)
    @test isapprox(out21.Intv.hi,1.0,atol=1E-5)

    out21a = step(Xn)
    @test isapprox(out21a.cc,0.0,atol=1E-5)
    @test isapprox(out21a.cv,0.0,atol=1E-5)
    @test isapprox(out21a.cc_grad[1],0.0,atol=1E-5)
    @test isapprox(out21a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out21a.cv_grad[1],0.0,atol=1E-5)
    @test isapprox(out21a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out21a.Intv.lo,0.0,atol=1E-5)
    @test isapprox(out21a.Intv.hi,0.0,atol=1E-5)

    out21b = step(Xz)
    @test isapprox(out21b.cc,0.5555555555555556,atol=1E-5)
    @test isapprox(out21b.cv,0.0,atol=1E-5)
    @test isapprox(out21b.cc_grad[1],0.444444,atol=1E-5)
    @test isapprox(out21b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out21b.cv_grad[1],0.0,atol=1E-5)
    @test isapprox(out21b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out21b.Intv.lo,0.0,atol=1E-5)
    @test isapprox(out21b.Intv.hi,1.0,atol=1E-5)
end

# MOSTLY DONE!
@testset "Sign Function" begin

    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    Xn = MC{2}(-4.0,-4.0,-xIBox[1],a,a,false)
    Xz = MC{2}(-2.0,-2.0,Interval(-3.0,1.0),a,a,false)

    out22 = sign(X)
    @test isapprox(out22.cc,1.0,atol=1E-5)
    @test isapprox(out22.cv,1.0,atol=1E-5)
    @test isapprox(out22.cc_grad[1],0.0,atol=1E-5)
    @test isapprox(out22.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out22.cv_grad[1],0.0,atol=1E-5)
    @test isapprox(out22.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out22.Intv.lo,1.0,atol=1E-5)
    @test isapprox(out22.Intv.hi,1.0,atol=1E-5)

    out22a = sign(Xn)
    @test isapprox(out22a.cc,-1.0,atol=1E-5)
    @test isapprox(out22a.cv,-1.0,atol=1E-5)
    @test isapprox(out22a.cc_grad[1],0.0,atol=1E-5)
    @test isapprox(out22a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out22a.cv_grad[1],0.0,atol=1E-5)
    @test isapprox(out22a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out22a.Intv.lo,-1.0,atol=1E-5)
    @test isapprox(out22a.Intv.hi,-1.0,atol=1E-5)

    out22b = sign(Xz)
    @test isapprox(out22b.cc,0.11111111111111116,atol=1E-5)
    @test isapprox(out22b.cv,-1.0,atol=1E-5)
    @test isapprox(out22b.cc_grad[1],0.888889,atol=1E-5)
    @test isapprox(out22b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out22b.cv_grad[1],0.0,atol=1E-5)
    @test isapprox(out22b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out22b.Intv.lo,-1.0,atol=1E-5)
    @test isapprox(out22b.Intv.hi,1.0,atol=1E-5)
end

# MOSTLY DONE!
#=
@testset "Absolute Value" begin
    EAGO.set_mc_differentiability!(0)
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0);Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)

    out = abs(X)
    @test isapprox(out.cc,6.409090909090908,atol=1E-1)
    @test isapprox(out.cv,4.5,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.454545,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],1.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,0,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    Y = MC{2}(7.0,7.0,xIBox[2],b,b,false)
    Xn = MC{2}(-4.0,-4.0,-xIBox[1],a,a,false)
    Xz = MC{2}(-2.0,-2.0,Interval(-3.0,1.0),a,a,false)

    out7 = abs(X)
    @test isapprox(out7.cc,atol=4.0,1E-5)
    @test isapprox(out7.cv,1.3061224489795915,atol=1E-5)
    @test isapprox(out7.cc_grad[1],1.0,atol=1E-5)
    @test isapprox(out7.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out7.cv_grad[1],0.979592,atol=1E-5)
    @test isapprox(out7.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out7.Intv.lo,3.0,atol=1E-5)
    @test isapprox(out7.Intv.hi,7.0,atol=1E-5)
end
=#
# MOSTLY DONE!
@testset "Sine" begin

    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    Y = MC{2}(7.0,7.0,xIBox[2],b,b,false)
    Xn = MC{2}(-4.0,-4.0,-xIBox[1],a,a,false)
    Xz = MC{2}(-2.0,-2.0,Interval(-3.0,1.0),a,a,false)

    out17 = sin(X)
    @test isapprox(out17.cc,0.2700866557245978,atol=1E-5)
    @test isapprox(out17.cv,-0.7568024953079283,atol=1E-5)
    @test isapprox(out17.cc_grad[1],0.128967,atol=1E-5)
    @test isapprox(out17.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out17.cv_grad[1],-0.653644,atol=1E-5)
    @test isapprox(out17.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out17.Intv.lo,-1,atol=1E-2)
    @test isapprox(out17.Intv.hi,0.656987,atol=1E-5)

    out17a = sin(Xn)
    @test isapprox(out17a.cc,0.7568024953079283,atol=1E-5)
    @test isapprox(out17a.cv,-0.2700866557245979,atol=1E-5)
    @test isapprox(out17a.cc_grad[1],-0.653644,atol=1E-5)
    @test isapprox(out17a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out17a.cv_grad[1],0.128967,atol=1E-5)
    @test isapprox(out17a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out17a.Intv.lo,-0.656987,atol=1E-2)
    @test isapprox(out17a.Intv.hi,1.0,atol=1E-5)

    out17b = sin(Xz)
    @test isapprox(out17b.cc,0.10452774015707458,atol=1E-5)
    @test isapprox(out17b.cv,-0.9092974268256817,atol=1E-5)
    @test isapprox(out17b.cc_grad[1],0.245648,atol=1E-5)
    @test isapprox(out17b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out17b.cv_grad[1],-0.416147,atol=1E-5)
    @test isapprox(out17b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out17b.Intv.lo,-1,atol=1E-2)
    @test isapprox(out17b.Intv.hi,0.841471,atol=1E-5)

    EAGO.set_mc_differentiability!(0)

    out17b = sin(Xz)
    @test isapprox(out17b.cc,0.10452774015707458,atol=1E-5)
    @test isapprox(out17b.cv,-0.9092974268256817,atol=1E-5)
    @test isapprox(out17b.cc_grad[1],0.245648,atol=1E-5)
    @test isapprox(out17b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out17b.cv_grad[1],-0.416147,atol=1E-5)
    @test isapprox(out17b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out17b.Intv.lo,-1,atol=1E-2)
    @test isapprox(out17b.Intv.hi,0.841471,atol=1E-5)
end

# MOSTLY DONE!
@testset "Cosine" begin

    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    Y = MC{2}(7.0,7.0,xIBox[2],b,b,false)
    Xn = MC{2}(-4.0,-4.0,-xIBox[1],a,a,false)
    Xz = MC{2}(-2.0,-2.0,Interval(-3.0,1.0),a,a,false)

    out19 = cos(X)
    @test isapprox(out19.cc,-0.31034065427934965,atol=1E-5)
    @test isapprox(out19.cv,-0.703492113936536,atol=1E-5)
    @test isapprox(out19.cc_grad[1],0.679652,atol=1E-5)
    @test isapprox(out19.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out19.cv_grad[1],0.485798,atol=1E-5)
    @test isapprox(out19.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out19.Intv.lo,-1.0,atol=1E-5)
    @test isapprox(out19.Intv.hi,1.0,atol=1E-5)

    out19a = cos(Xn)
    @test isapprox(out19a.cc,-0.31034065427934965,atol=1E-5)
    @test isapprox(out19a.cv,-0.703492113936536,atol=1E-5)
    @test isapprox(out19a.cc_grad[1],-0.679652,atol=1E-5)
    @test isapprox(out19a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out19a.cv_grad[1],-0.485798,atol=1E-5)
    @test isapprox(out19a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out19a.Intv.lo,-1.0,atol=1E-5)
    @test isapprox(out19a.Intv.hi,1.0,atol=1E-5)

    out19b = cos(Xz)
    @test isapprox(out19b.cc,-0.222468094224762,atol=1E-5)
    @test isapprox(out19b.cv,-0.6314158569813042,atol=1E-5)
    @test isapprox(out19b.cc_grad[1],0.76752,atol=1E-5)
    @test isapprox(out19b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out19b.cv_grad[1],0.390573,atol=1E-5)
    @test isapprox(out19b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out19b.Intv.lo,-0.989993,atol=1E-5)
    @test isapprox(out19b.Intv.hi,1.0,atol=1E-5)

    EAGO.set_mc_differentiability!(0)

    out19b = cos(Xz)
    @test isapprox(out19b.cc,-0.222468094224762,atol=1E-5)
    @test isapprox(out19b.cv,-0.6314158569813042,atol=1E-5)
    @test isapprox(out19b.cc_grad[1],0.76752,atol=1E-5)
    @test isapprox(out19b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out19b.cv_grad[1],0.390573,atol=1E-5)
    @test isapprox(out19b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out19b.Intv.lo,-0.989993,atol=1E-5)
    @test isapprox(out19b.Intv.hi,1.0,atol=1E-5)
end

#=
# MOSTLY DONE!
@testset "Tangent" begin

    EAGO.set_mc_differentiability!(0)
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(0.5,1.0);Interval(-0.5,0.5)])
    X = MC{2}(0.6,0.6,xIBox[1],a,a,false)
    Xn = MC{2}(-0.8,-0.8,-xIBox[1],a,a,false)
    Xz = MC{2}(-0.3,-0.3,xIBox[2],a,a,false)
    Xerr = MC{2}(0.6,0.6,Interval(-4.5,5.0),a,a,false)

    out19 = tan(X)
    @test isapprox(out19.cc,0.7485235368060128,atol=1E-5)
    @test isapprox(out19.cv,0.6841368083416923,atol=1E-5)
    @test isapprox(out19.cc_grad[1],2.02221,atol=1E-5)
    @test isapprox(out19.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out19.cv_grad[1],1.46804,atol=1E-5)
    @test isapprox(out19.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out19.Intv.lo,0.546302,atol=1E-5)
    @test isapprox(out19.Intv.hi,1.55741,atol=1E-5)

    out19a = tan(Xn)
    @test isapprox(out19a.cc,-1.0296385570503641,atol=1E-5)
    @test isapprox(out19a.cv,-1.1529656307304577,atol=1E-5)
    @test isapprox(out19a.cc_grad[1],2.06016,atol=1E-5)
    @test isapprox(out19a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out19a.cv_grad[1],2.02221,atol=1E-5)
    @test isapprox(out19a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out19a.Intv.lo,-1.55741,atol=1E-5)
    @test isapprox(out19a.Intv.hi,-0.546302,atol=1E-5)

    out19b = tan(Xz)
    @test isapprox(out19b.cc,-0.30933624960962325,atol=1E-5)
    @test isapprox(out19b.cv,-0.332534,atol=1E-5)
    @test isapprox(out19b.cc_grad[1],1.09569,atol=1E-5)
    @test isapprox(out19b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out19b.cv_grad[1],1.06884,atol=1E-5)
    @test isapprox(out19b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out19b.Intv.lo,-0.546303,atol=1E-5)
    @test isapprox(out19b.Intv.hi,0.546303,atol=1E-5)

    EAGO.set_mc_differentiability!(1)

    out19c = tan(Xz)
    @test isapprox(out19c.cc,-0.309336,atol=1E-5)
    @test isapprox(out19c.cv,-0.332534,atol=1E-5)
    @test isapprox(out19c.cc_grad[1],1.09569,atol=1E-5)
    @test isapprox(out19c.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out19c.cv_grad[1],1.06884,atol=1E-5)
    @test isapprox(out19c.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out19c.Intv.lo,-0.546303,atol=1E-5)
    @test isapprox(out19c.Intv.hi,0.546303,atol=1E-5)

    @test_throws ErrorException tan(Xerr)
end

# MOSTLY DONE!
@testset "Inverse Sine" begin

    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-0.9,-0.5);Interval(-0.5,0.5)])
    X = MC{2}(-0.7,-0.7,xIBox[1],a,a,false)
    Xn = MC{2}(0.7,0.7,-xIBox[1],a,a,false)
    Xz = MC{2}(-0.1,-0.1,xIBox[2],a,a,false)

    out17 = asin(X)
    @test isapprox(out17.cc,-0.775397496610753,atol=1E-5)
    @test isapprox(out17.cv,-0.8216841452984665,atol=1E-5)
    @test isapprox(out17.cc_grad[1],1.40028,atol=1E-5)
    @test isapprox(out17.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out17.cv_grad[1],1.49043,atol=1E-5)
    @test isapprox(out17.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out17.Intv.lo,-1.11977,atol=1E-2)
    @test isapprox(out17.Intv.hi,-0.523598,atol=1E-5)

    out17a = asin(Xn)
    @test isapprox(out17a.cc,0.8216841452984665,atol=1E-5)
    @test isapprox(out17a.cv,0.775397496610753,atol=1E-5)
    @test isapprox(out17a.cc_grad[1],1.49043,atol=1E-5)
    @test isapprox(out17a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out17a.cv_grad[1],1.40028,atol=1E-5)
    @test isapprox(out17a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out17a.Intv.lo,0.523598,atol=1E-2)
    @test isapprox(out17a.Intv.hi,1.11977,atol=1E-5)

    out17b = asin(Xz)
    @test isapprox(out17b.cc,-0.0974173098978382,atol=1E-5)
    @test isapprox(out17b.cv,-0.10958805193420748,atol=1E-5)
    @test isapprox(out17b.cc_grad[1],1.03503,atol=1E-5)
    @test isapprox(out17b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out17b.cv_grad[1],1.03503,atol=1E-5)
    @test isapprox(out17b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out17b.Intv.lo,-0.523599,atol=1E-2)
    @test isapprox(out17b.Intv.hi,0.523599,atol=1E-5)

    EAGO.set_mc_differentiability!(0)
    out17b = asin(Xz)
    @test isapprox(out17b.cc,-0.0974173098978382,atol=1E-5)
    @test isapprox(out17b.cv,-0.10958805193420748,atol=1E-5)
    @test isapprox(out17b.cc_grad[1],1.03503,atol=1E-5)
    @test isapprox(out17b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out17b.cv_grad[1],1.03503,atol=1E-5)
    @test isapprox(out17b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out17b.Intv.lo,-0.523599,atol=1E-2)
    @test isapprox(out17b.Intv.hi,0.523599,atol=1E-5)
end

#=
@testset "Inverse Cosine" begin
end
=#

# MOSTLY DONE!
@testset "Inverse Tangent" begin
    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)

    out16 = atan(X)
    @test isapprox(out16.cc,1.3258176636680326,atol=1E-5)
    @test isapprox(out16.cv,1.294009147346374,atol=1E-5)
    @test isapprox(out16.cc_grad[1],0.0588235,atol=1E-5)
    @test isapprox(out16.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out16.cv_grad[1],0.0449634,atol=1E-5)
    @test isapprox(out16.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out16.Intv.lo,1.24904,atol=1E-3)
    @test isapprox(out16.Intv.hi,1.4289,atol=1E-3)

    out16a = atan(Xn)
    @test isapprox(out16a.cc,-1.294009147346374,atol=1E-5)
    @test isapprox(out16a.cv,-1.3258176636680326,atol=1E-5)
    @test isapprox(out16a.cc_grad[1],0.0449634,atol=1E-5)
    @test isapprox(out16a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out16a.cv_grad[1],.0588235,atol=1E-5)
    @test isapprox(out16a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out16a.Intv.lo,-1.4289,atol=1E-3)
    @test isapprox(out16a.Intv.hi,-1.24904,atol=1E-3)

    out16b = atan(Xz)
    @test isapprox(out16b.cc,-0.7404162771337869,atol=1E-5)
    @test isapprox(out16b.cv,-1.1071487177940904,atol=1E-5)
    @test isapprox(out16b.cc_grad[1],0.508629,atol=1E-5)
    @test isapprox(out16b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out16b.cv_grad[1],0.2,atol=1E-5)
    @test isapprox(out16b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out16b.Intv.lo,-1.24905,atol=1E-3)
    @test isapprox(out16b.Intv.hi,0.785399,atol=1E-3)
end

# MOSTLY DONE!
@testset "Hyperbolic Sine" begin
    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    Y = MC{2}(7.0,7.0,xIBox[2],b,b,false)
    Xn = MC{2}(-4.0,-4.0,-xIBox[1],a,a,false)
    Xz = MC{2}(-2.0,-2.0,Interval(-3.0,1.0),a,a,false)

    out10 = sinh(X)
    @test isapprox(out10.cc,144.59243701386904,atol=1E-5)
    @test isapprox(out10.cv,27.28991719712775,atol=1E-5)
    @test isapprox(out10.cc_grad[1],134.575,atol=1E-2)
    @test isapprox(out10.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10.cv_grad[1],27.3082,atol=1E-2)
    @test isapprox(out10.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10.Intv.lo,10.0178,atol=1E-2)
    @test isapprox(out10.Intv.hi,548.317,atol=1E-2)

    out10a = sinh(Xn)
    @test isapprox(out10a.cc,-27.28991719712775,atol=1E-5)
    @test isapprox(out10a.cv,-144.59243701386904,atol=1E-5)
    @test isapprox(out10a.cc_grad[1],27.3082,atol=1E-2)
    @test isapprox(out10a.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10a.cv_grad[1],134.575,atol=1E-2)
    @test isapprox(out10a.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10a.Intv.lo,-548.317,atol=1E-2)
    @test isapprox(out10a.Intv.hi,-10.0178,atol=1E-2)

    out10b = sinh(Xz)
    @test isapprox(out10b.cc,-3.626860407847019,atol=1E-5)
    @test isapprox(out10b.cv,-7.219605897146477,atol=1E-5)
    @test isapprox(out10b.cc_grad[1],3.7622,atol=1E-2)
    @test isapprox(out10b.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10b.cv_grad[1],2.79827,atol=1E-2)
    @test isapprox(out10b.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10b.Intv.lo,-10.0179,atol=1E-2)
    @test isapprox(out10b.Intv.hi,1.17521,atol=1E-2)

    EAGO.set_mc_differentiability!(0)
    out10 = sinh(X)
    @test isapprox(out10.cc,144.59243701386904,atol=1E-5)
    @test isapprox(out10.cv,27.28991719712775,atol=1E-5)
    @test isapprox(out10.cc_grad[1],134.575,atol=1E-2)
    @test isapprox(out10.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10.cv_grad[1],27.3082,atol=1E-2)
    @test isapprox(out10.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10.Intv.lo,10.0178,atol=1E-2)
    @test isapprox(out10.Intv.hi,548.317,atol=1E-2)
end

# MOSTLY DONE!
@testset "Hyperbolic Cosine" begin
    EAGO.set_mc_differentiability!(0)
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    Y = MC{2}(7.0,7.0,xIBox[2],b,b,false)

    out8 = cosh(X)
    @test isapprox(out8.cc,144.63000528563632,atol=1E-5)
    @test isapprox(out8.cv,27.308232836016487,atol=1E-5)
    @test isapprox(out8.cc_grad[1],134.562,atol=1E-2)
    @test isapprox(out8.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out8.cv_grad[1],-27.2899,atol=1E-3)
    @test isapprox(out8.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out8.Intv.lo,10.0676,atol=1E-3)
    @test isapprox(out8.Intv.hi,548.318,atol=1E-3)

    EAGO.set_mc_differentiability!(1)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    mBox = mid.(xIBox)
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    Y = MC{2}(7.0,7.0,xIBox[2],b,b,false)
    Xn = MC{2}(-4.0,-4.0,-xIBox[1],a,a,false)
    Xz = MC{2}(-2.0,-2.0,Interval(-3.0,1.0),a,a,false)

    out8 = cosh(X)
    @test isapprox(out8.cc,144.63000528563632,atol=1E-5)
    @test isapprox(out8.cv,27.308232836016487,atol=1E-5)
    @test isapprox(out8.cc_grad[1],134.562,atol=1E-2)
    @test isapprox(out8.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out8.cv_grad[1],-27.2899,atol=1E-3)
    @test isapprox(out8.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out8.Intv.lo,10.0676,atol=1E-3)
    @test isapprox(out8.Intv.hi,548.318,atol=1E-3)
end

# MOSTLY DONE!
@testset "Hyperbolic Tangent" begin

    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    Y = MC{2}(7.0,7.0,xIBox[2],b,b,false)
    Xn = MC{2}(-4.0,-4.0,-xIBox[1],a,a,false)
    Xz = MC{2}(-2.0,-2.0,Interval(-3.0,1.0),a,a,false)
    Xz1 = MC{2}(2.0,2.0,Interval(-1.0,3.0),a,a,false)

    out12 = tanh(X)
    @test isapprox(out12.cc,0.999329299739067,atol=1E-5)
    @test isapprox(out12.cv,0.996290649501034,atol=1E-5)
    @test isapprox(out12.cc_grad[1],0.00134095,atol=1E-5)
    @test isapprox(out12.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out12.cv_grad[1],0.0012359,atol=1E-5)
    @test isapprox(out12.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out12.Intv.lo,0.995054,atol=1E-5)
    @test isapprox(out12.Intv.hi,0.999999,atol=1E-5)

    out12a = tanh(Xn)
    @test isapprox(out12a.cc,-0.996290649501034,atol=1E-5)
    @test isapprox(out12a.cv,-0.999329299739067,atol=1E-5)
    @test isapprox(out12a.cc_grad[1],0.0012359,atol=1E-5)
    @test isapprox(out12a.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out12a.cv_grad[1],0.00134095,atol=1E-5)
    @test isapprox(out12a.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out12a.Intv.lo,-0.999999,atol=1E-5)
    @test isapprox(out12a.Intv.hi,-0.995054,atol=1E-5)

    out12b = tanh(Xz)
    @test isapprox(out12b.cc,-0.5558207301372651,atol=1E-5)
    @test isapprox(out12b.cv,-0.9640275800758169,atol=1E-5)
    @test isapprox(out12b.cc_grad[1],0.439234,atol=1E-5)
    @test isapprox(out12b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out12b.cv_grad[1],0.0706508,atol=1E-5)
    @test isapprox(out12b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out12b.Intv.lo,-0.995055,atol=1E-5)
    @test isapprox(out12b.Intv.hi,0.761595,atol=1E-5)

    out12b = tanh(Xz1)
    @test isapprox(out12b.cc,0.9640275800758169,atol=1E-5)
    @test isapprox(out12b.cv,0.5558207301372651,atol=1E-5)
    @test isapprox(out12b.cc_grad[1],0.0706508,atol=1E-5)
    @test isapprox(out12b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out12b.cv_grad[1],0.439234,atol=1E-5)
    @test isapprox(out12b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out12b.Intv.lo,-0.761595,atol=1E-5)
    @test isapprox(out12b.Intv.hi,0.995055,atol=1E-5)

    EAGO.set_mc_differentiability!(0)
    out12b = tanh(Xz)
    @test isapprox(out12b.cc,-0.5558207301372651,atol=1E-5)
    @test isapprox(out12b.cv,-0.9640275800758169,atol=1E-5)
    @test isapprox(out12b.cc_grad[1],0.439234,atol=1E-5)
    @test isapprox(out12b.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out12b.cv_grad[1],0.0706508,atol=1E-5)
    @test isapprox(out12b.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out12b.Intv.lo,-0.995055,atol=1E-5)
    @test isapprox(out12b.Intv.hi,0.761595,atol=1E-5)
end

# MOSTLY DONE!
@testset "Inverse Hyperbolic Sine" begin

    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(0.1,0.7);Interval(-3,7)])
    X = MC{2}(0.3,0.3,xIBox[1],a,a,false)
    Xn = MC{2}(-0.3,-0.3,-xIBox[1],a,a,false)
    Xz = MC{2}(2.0,2.0,xIBox[2],a,a,false)
    Xz1 = MC{2}(-2.0,-2.0,-xIBox[2],a,a,false)

    out10 = asinh(X)
    @test isapprox(out10.cc,0.29567304756342244,atol=1E-5)
    @test isapprox(out10.cv,0.2841115746269236,atol=1E-5)
    @test isapprox(out10.cc_grad[1],0.957826,atol=1E-2)
    @test isapprox(out10.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10.cv_grad[1],0.921387,atol=1E-2)
    @test isapprox(out10.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10.Intv.lo,0.099834,atol=1E-2)
    @test isapprox(out10.Intv.hi,0.652667,atol=1E-2)

    out10a = asinh(Xn)
    @test isapprox(out10a.cc,-0.2841115746269236,atol=1E-5)
    @test isapprox(out10a.cv,-0.29567304756342244,atol=1E-5)
    @test isapprox(out10a.cc_grad[1],0.921387,atol=1E-2)
    @test isapprox(out10a.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10a.cv_grad[1],0.957826,atol=1E-2)
    @test isapprox(out10a.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10a.Intv.lo,-0.652667,atol=1E-2)
    @test isapprox(out10a.Intv.hi,-0.099834,atol=1E-2)

    out10b = asinh(Xz)
    @test isapprox(out10b.cc,1.4436354751788103,atol=1E-5)
    @test isapprox(out10b.cv,0.3730697449603356,atol=1E-5)
    @test isapprox(out10b.cc_grad[1],0.447214,atol=1E-2)
    @test isapprox(out10b.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10b.cv_grad[1],0.45421,atol=1E-2)
    @test isapprox(out10b.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10b.Intv.lo,-1.81845,atol=1E-2)
    @test isapprox(out10b.Intv.hi,2.64413,atol=1E-2)

    out10c = asinh(Xz1)
    @test isapprox(out10c.cc,-0.3730697449603356,atol=1E-5)
    @test isapprox(out10c.cv,-1.4436354751788103,atol=1E-5)
    @test isapprox(out10c.cc_grad[1],0.45421,atol=1E-2)
    @test isapprox(out10c.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10c.cv_grad[1],0.447214,atol=1E-2)
    @test isapprox(out10c.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10c.Intv.lo,-2.64413,atol=1E-2)
    @test isapprox(out10c.Intv.hi,1.81845,atol=1E-2)

    EAGO.set_mc_differentiability!(0)

    out10d = asinh(X)
    @test isapprox(out10d.cc,0.29567304756342244,atol=1E-5)
    @test isapprox(out10d.cv,0.2841115746269236,atol=1E-5)
    @test isapprox(out10d.cc_grad[1],0.957826,atol=1E-2)
    @test isapprox(out10d.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10d.cv_grad[1],0.921387,atol=1E-2)
    @test isapprox(out10d.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10d.Intv.lo,0.099834,atol=1E-2)
    @test isapprox(out10d.Intv.hi,0.652667,atol=1E-2)
end

# MOSTLY DONE!
@testset "Inverse Hyperbolic Cosine" begin
    a = seed_gradient(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    X = MC{2}(4.0,4.0,xIBox[1],a,a,false)
    out9 = acosh(X)

    @test isapprox(out9.cc,2.0634370688955608,atol=1E-5)
    @test isapprox(out9.cv,1.9805393289917226,atol=1E-5)
    @test isapprox(out9.cc_grad[1],0.258199,atol=1E-5)
    @test isapprox(out9.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out9.cv_grad[1],0.217792,atol=1E-5)
    @test isapprox(out9.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out9.Intv.lo,1.76274,atol=1E-5)
    @test isapprox(out9.Intv.hi,2.63392,atol=1E-5)
end

# MOSTLY DONE!
@testset "Inverse Hyperbolic Tangent" begin

    EAGO.set_mc_differentiability!(1)
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(0.1,0.7);Interval(-0.3,0.7)])
    X = MC{2}(0.3,0.3,xIBox[1],a,a,false)
    Xn = MC{2}(-0.3,-0.3,-xIBox[1],a,a,false)
    Xz = MC{2}(0.2,0.2,xIBox[2],a,a,false)
    Xz1 = MC{2}(-0.2,-0.2,-xIBox[2],a,a,false)

    out10 = atanh(X)
    @test isapprox(out10.cc,0.3559904077187347,atol=1E-5)
    @test isapprox(out10.cv,0.30951960420311175,atol=1E-5)
    @test isapprox(out10.cc_grad[1],1.27828,atol=1E-2)
    @test isapprox(out10.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10.cv_grad[1],1.0989,atol=1E-2)
    @test isapprox(out10.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10.Intv.lo,0.100335,atol=1E-2)
    @test isapprox(out10.Intv.hi,0.867301,atol=1E-2)

    out10a = atanh(Xn)
    @test isapprox(out10a.cc,-0.30951960420311175,atol=1E-5)
    @test isapprox(out10a.cv,-0.3559904077187347,atol=1E-5)
    @test isapprox(out10a.cc_grad[1],1.0989,atol=1E-2)
    @test isapprox(out10a.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10a.cv_grad[1],1.27828,atol=1E-2)
    @test isapprox(out10a.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10a.Intv.lo,-0.867301,atol=1E-2)
    @test isapprox(out10a.Intv.hi,-0.100335,atol=1E-2)

    out10b = atanh(Xz)
    @test isapprox(out10b.cc,0.2788904617454707,atol=1E-5)
    @test isapprox(out10b.cv,0.2027325540540822,atol=1E-5)
    @test isapprox(out10b.cc_grad[1],1.17682,atol=1E-2)
    @test isapprox(out10b.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10b.cv_grad[1],1.04167,atol=1E-2)
    @test isapprox(out10b.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10b.Intv.lo,-0.30952,atol=1E-2)
    @test isapprox(out10b.Intv.hi,0.867301,atol=1E-2)

    out10c = atanh(Xz1)
    @test isapprox(out10c.cc,-0.2027325540540822,atol=1E-5)
    @test isapprox(out10c.cv,-0.2788904617454707,atol=1E-5)
    @test isapprox(out10c.cc_grad[1],1.04167,atol=1E-2)
    @test isapprox(out10c.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10c.cv_grad[1],1.17682,atol=1E-2)
    @test isapprox(out10c.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10c.Intv.lo,-0.867301,atol=1E-2)
    @test isapprox(out10c.Intv.hi,0.30952,atol=1E-2)

    EAGO.set_mc_differentiability!(0)

    out10d = atanh(X)
    @test isapprox(out10d.cc,0.3559904077187347,atol=1E-5)
    @test isapprox(out10d.cv,0.30951960420311175,atol=1E-5)
    @test isapprox(out10d.cc_grad[1],1.27828,atol=1E-2)
    @test isapprox(out10d.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out10d.cv_grad[1],1.0989,atol=1E-2)
    @test isapprox(out10d.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out10d.Intv.lo,0.100335,atol=1E-2)
    @test isapprox(out10d.Intv.hi,0.867301,atol=1E-2)
end
=#
