module Op_Trignometric

using Compat
using Compat.Test
using EAGO
using IntervalArithmetic
using StaticArrays

@testset "Test Sin" begin

    EAGO.set_diff_relax(1)
    a = seed_g(Float64,1,2)
    b = seed_g(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.0,4.0,a,a,xIBox[1],false,xIBox,mBox)
    Y = SMCg{2,Interval{Float64},Float64}(7.0,7.0,b,b,xIBox[2],false,xIBox,mBox)
    Xn = SMCg{2,Interval{Float64},Float64}(-4.0,-4.0,a,a,-xIBox[1],false,xIBox,mBox)
    Xz = SMCg{2,Interval{Float64},Float64}(-2.0,-2.0,a,a,Interval(-3.0,1.0),false,xIBox,mBox)

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

    EAGO.set_diff_relax(0)

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

@testset "Test Cos" begin

    EAGO.set_diff_relax(1)
    a = seed_g(Float64,1,2)
    b = seed_g(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0);Interval(3.0,9.0)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(4.0,4.0,a,a,xIBox[1],false,xIBox,mBox)
    Y = SMCg{2,Interval{Float64},Float64}(7.0,7.0,b,b,xIBox[2],false,xIBox,mBox)
    Xn = SMCg{2,Interval{Float64},Float64}(-4.0,-4.0,a,a,-xIBox[1],false,xIBox,mBox)
    Xz = SMCg{2,Interval{Float64},Float64}(-2.0,-2.0,a,a,Interval(-3.0,1.0),false,xIBox,mBox)

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

    EAGO.set_diff_relax(0)

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

@testset "Test Tan" begin

    EAGO.set_diff_relax(0)
    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(0.5,1.0);Interval(-0.5,0.5)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(0.6,0.6,a,a,xIBox[1],false,xIBox,mBox)
    Xn = SMCg{2,Interval{Float64},Float64}(-0.8,-0.8,a,a,-xIBox[1],false,xIBox,mBox)
    Xz = SMCg{2,Interval{Float64},Float64}(-0.3,-0.3,a,a,xIBox[2],false,xIBox,mBox)
    Xerr = SMCg{2,Interval{Float64},Float64}(0.6,0.6,a,a,Interval(-4.5,5.0),false,xIBox,mBox)

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

    EAGO.set_diff_relax(1)

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

end
