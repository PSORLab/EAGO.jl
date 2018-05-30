module Op_Hyperbolic

using Compat
using Compat.Test
using EAGO
using IntervalArithmetic
using StaticArrays

@testset "Test Sinh" begin
    # ADD nonsmooth test
    EAGO.set_diff_relax(1)
    a = seed_g(Float64,1,2)
    b = seed_g(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.0,4.0,a,a,xIBox[1],false)
    Y = SMCg{2,Interval{Float64},Float64}(7.0,7.0,b,b,xIBox[2],false)
    Xn = SMCg{2,Interval{Float64},Float64}(-4.0,-4.0,a,a,-xIBox[1],false)
    Xz = SMCg{2,Interval{Float64},Float64}(-2.0,-2.0,a,a,Interval(-3.0,1.0),false)

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

    EAGO.set_diff_relax(0)
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

@testset "Test Cosh" begin

    EAGO.set_diff_relax(0)
    a = seed_g(Float64,1,2)
    b = seed_g(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.0,4.0,a,a,xIBox[1],false)
    Y = SMCg{2,Interval{Float64},Float64}(7.0,7.0,b,b,xIBox[2],false)

    out8 = cosh(X)
    @test isapprox(out8.cc,144.63000528563632,atol=1E-5)
    @test isapprox(out8.cv,27.308232836016487,atol=1E-5)
    @test isapprox(out8.cc_grad[1],134.562,atol=1E-2)
    @test isapprox(out8.cc_grad[2],0.0,atol=1E-5)
    @test isapprox(out8.cv_grad[1],-27.2899,atol=1E-3)
    @test isapprox(out8.cv_grad[2],0.0,atol=1E-5)
    @test isapprox(out8.Intv.lo,10.0676,atol=1E-3)
    @test isapprox(out8.Intv.hi,548.318,atol=1E-3)

    EAGO.set_diff_relax(1)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.0,4.0,a,a,xIBox[1],false)
    Y = SMCg{2,Interval{Float64},Float64}(7.0,7.0,b,b,xIBox[2],false)
    Xn = SMCg{2,Interval{Float64},Float64}(-4.0,-4.0,a,a,-xIBox[1],false)
    Xz = SMCg{2,Interval{Float64},Float64}(-2.0,-2.0,a,a,Interval(-3.0,1.0),false)

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

@testset "Test Tanh" begin

    EAGO.set_diff_relax(1)
    a = seed_g(Float64,1,2)
    b = seed_g(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.0,4.0,a,a,xIBox[1],false)
    Y = SMCg{2,Interval{Float64},Float64}(7.0,7.0,b,b,xIBox[2],false)
    Xn = SMCg{2,Interval{Float64},Float64}(-4.0,-4.0,a,a,-xIBox[1],false)
    Xz = SMCg{2,Interval{Float64},Float64}(-2.0,-2.0,a,a,Interval(-3.0,1.0),false)
    Xz1 = SMCg{2,Interval{Float64},Float64}(2.0,2.0,a,a,Interval(-1.0,3.0),false)

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

    EAGO.set_diff_relax(0)
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

end
