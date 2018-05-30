module Op_InvTrignometric

using Compat
using Compat.Test
using EAGO
using IntervalArithmetic
using StaticArrays


# FIX CROSS OVER
@testset "Test Asin" begin

    EAGO.set_diff_relax(1)
    a = seed_g(Float64,1,2)
    b = seed_g(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-0.9,-0.5);Interval(-0.5,0.5)])
    mBox = mid.(xIBox)
    X = SMCg{2,Interval{Float64},Float64}(-0.7,-0.7,a,a,xIBox[1],false)
    Xn = SMCg{2,Interval{Float64},Float64}(0.7,0.7,a,a,-xIBox[1],false)
    Xz = SMCg{2,Interval{Float64},Float64}(-0.1,-0.1,a,a,xIBox[2],false)

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

    EAGO.set_diff_relax(0)
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

@testset "Test Acos" begin
end

@testset "Test Atan" begin
end

end
