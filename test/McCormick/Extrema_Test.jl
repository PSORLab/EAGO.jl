module Extrema_Test

using Compat
using Compat.Test
using EAGO
using IntervalArithmetic
using StaticArrays


@testset "Test Multivariate Max/Min" begin

    EAGO.set_diff_relax(0)

    seed1 = seed_g(Float64,1,2)
    seed2 = seed_g(Float64,2,2)


    X = SMCg{2,Interval{Float64},Float64}(129.625,129.625,seed1,seed1,Interval(109.349, 149.901),false)
    Y = SMCg{2,Interval{Float64},Float64}(124.25,124.25,seed2,seed2,Interval(120.5,139.0),false)
    out = max(X,Y)

    @test isapprox(out.cc,141.689,atol=1E-1)
    @test isapprox(out.cv,129.625,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.497883,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.502117,atol=1E-4)
    @test isapprox(out.cv_grad[1],1.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,120.5,atol=1E-2)
    @test isapprox(out.Intv.hi,149.901,atol=1E-2)

    EAGO.set_diff_relax(1)

    X = SMCg{2,Interval{Float64},Float64}(129.625,129.625,seed1,seed1,Interval(109.349, 149.901),false)
    Y = SMCg{2,Interval{Float64},Float64}(124.25,124.25,seed2,seed2,Interval(120.5, 139.0),false)
    out = max(X,Y)

    @test isapprox(out.cc,144.4505,atol=1E-3)
    @test isapprox(out.cv,124.42964337332246,atol=1E-3)
    @test isapprox(out.cc_grad[1],0.3897811551934656,atol=1E-3)
    @test isapprox(out.cv_grad[1],0.10026606883114281,atol=1E-3)
    @test isapprox(out.cc_grad[2],0.2651570489408331,atol=1E-3)
    @test isapprox(out.cv_grad[2],0.8997339311688572,atol=1E-3)
    @test isapprox(out.Intv.lo,120.5,atol=1E-3)
    @test isapprox(out.Intv.hi,149.902,atol=1E-3)

    EAGO.set_diff_relax(0)

end
#=
a = seed_g(Float64,Int64(1),Int64(2))
xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
mBox = mid.(xIBox)
X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false,xIBox,mBox)
out = min(3,X)
println("out 3: $out")
@test about(out.cc,3.0,1E-1)
@test about(out.cv,1.0909090909090908,1E-1)
@test about(out.cc_grad[1],0.0,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],0.545455,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,-3,1E-4)
@test about(out.Intv.hi,3,1E-4)

a = seed_g(Float64,Int64(1),Int64(2))
xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
mBox = mid.(xIBox)
X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false,xIBox,mBox)
out = max(5,X)
println("out 4: $out")
@test about(out.cc,7.045454545454545,1E-1)
@test about(out.cv,5.0,1E-1)
@test about(out.cc_grad[1],0.272727,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],0.0,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,5,1E-4)
@test about(out.Intv.hi,8,1E-4)
=#
end
