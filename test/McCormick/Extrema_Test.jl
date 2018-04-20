module Extrema_Test

using Compat
using Compat.Test
using EAGO
using IntervalArithmetic
using StaticArrays

function about(calc,val,tol)
    return (val - tol <= calc <= val + tol)
end

EAGO.set_diff_relax(0)

a = seed_g(Float64,Int64(1),Int64(2))
xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
mBox = mid.(xIBox)
X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,a,a,xIBox[1],false,xIBox,mBox)
out = min(3,X)
println("out 1: $out")
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
println("out 2: $out")
@test about(out.cc,7.045454545454545,1E-1)
@test about(out.cv,5.0,1E-1)
@test about(out.cc_grad[1],0.272727,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],0.0,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,5,1E-4)
@test about(out.Intv.hi,8,1E-4)


EAGO.set_diff_relax(1)
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
