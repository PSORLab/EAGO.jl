module Op_InvHyperbolic

using Compat
using Compat.Test
using EAGO
using IntervalArithmetic
using StaticArrays

@testset "Test Asinh" begin

    EAGO.set_diff_relax(1)
    a = seed_g(Float64,1,2)
    b = seed_g(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(0.1,0.7);Interval(-3,7)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(0.3,0.3,a,a,xIBox[1],false,xIBox,mBox)
    Xn = SMCg{2,Interval{Float64},Float64}(-0.3,-0.3,a,a,-xIBox[1],false,xIBox,mBox)
    Xz = SMCg{2,Interval{Float64},Float64}(2.0,2.0,a,a,xIBox[2],false,xIBox,mBox)
    Xz1 = SMCg{2,Interval{Float64},Float64}(-2.0,-2.0,a,a,-xIBox[2],false,xIBox,mBox)

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

    EAGO.set_diff_relax(0)

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

@testset "Test Acosh" begin
end

@testset "Test Atanh" begin
end

end
