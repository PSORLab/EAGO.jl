module Utilities

using Compat
using Compat.Test
using EAGO
using IntervalArithmetic
using StaticArrays

# create seed gradient
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)

# create SmoothMcCormick seed object for x1 = 2.0 on [1.0,3.0] for relaxing
# a function f(x1,x2) on the interval box xIbox using mBox as a reference point
x = 2.0
xIntv1 = Interval(1.0,3.0)
xIBox =  SVector{2,Interval{Float64}}([xIntv1,xIntv1])
mBox = mid.(xIBox)
SMCg1 = SMCg{2,Interval{Float64},Float64}(x,x,a,a,xIntv1,false,xIBox,mBox)

# resets gradient to seed gradient with seed at j=2
grad(SMCg1,2)
@test SMCg1.cv_grad == SVector{2,Float64}([1,0])
@test SMCg1.cc_grad == SVector{2,Float64}([1,0])

# sets gradients to zero
SMCg1 = zgrad(SMCg1)
@test SMCg1.cv_grad == SVector{2,Float64}([0,0])
@test SMCg1.cc_grad == SVector{2,Float64}([0,0])

# checks that mid3 returns midpoint
mid_return = mid3(9.0,-3.0,5.0)
@test mid_return[1] == 5.0
@test mid_return[2] == 3

# checks that mid_grad is working
mg_return1 = mid_grad(a, b, 1)
mg_return2 = mid_grad(a, b, 2)
mg_return3 = mid_grad(a, b, 3)
@test mg_return1 == SVector{2,Float64}([1,0])
@test mg_return2 == SVector{2,Float64}([0,1])
@test mg_return3 == SVector{2,Float64}([0,0])

# generates line segment value and derivative
y0 = line_seg(2.0,1.0,0.0,6.0,10.0)
dy0 = dline_seg(2.0,1.0,0.0,6.0,10.0,1.0)
@test y0 == 2.0
@test dy0 == 2.0

# converts float/integer to SMCg
#=
promoted1 = promote_rule(EAGOSmoothMcCormickGrad.SMCg{3,Float64}, Int64)
promoted2 = promote_rule(EAGOSmoothMcCormickGrad.SMCg{3,Float32}, Float64)
promoted3 = promote_rule(EAGOSmoothMcCormickGrad.SMCg{3,Float64}, Interval{Float64})
promoted4 = promote_rule(EAGOSmoothMcCormickGrad.SMCg{3,Float64}, Real)
@test promoted1 == EAGOSmoothMcCormickGrad.SMCg{3,Float64}
@test promoted2 == EAGOSmoothMcCormickGrad.SMCg{3,Float32}
@test promoted3 == EAGOSmoothMcCormickGrad.SMCg{3,Float64}
@test promoted4 == EAGOSmoothMcCormickGrad.SMCg{3,Float64}
=#
#=
conv_float = convert(EAGOSmoothMcCormickGrad.SMCg{3,Float64},1.0)
conv_int = convert(EAGOSmoothMcCormickGrad.SMCg{3,Float64},3)

@test conv_float == EAGOSmoothMcCormickGrad.SMCg{3,Float64}(1.0, 1.0, SVector{3,Float64}([0.0, 0.0, 0.0]), SVector{3,Float64}([0.0, 0.0, 0.0]), Interval(1, 1), false, IntervalArithmetic.Interval{Float64}[∅], [0.0])
@test conv_int == EAGOSmoothMcCormickGrad.SMCg{3,Float64}(3.0, 3.0, SVector{3,Float64}([0.0, 0.0, 0.0]), SVector{3,Float64}([0.0, 0.0, 0.0]), Interval(3, 3), false, IntervalArithmetic.Interval{Float64}[∅], [0.0])
=#

#=
# test subgradient refinement
Intv1 = tighten_subgrad(2.0,1.0,b,b,Interval(-100.0,100.0),[Interval(-1.0,1.0)],[0.0])
Intv2 = tighten_subgrad(2.0,1.0,b,b,Interval(-100.0,100.0),[∅],[0.0])
@test Intv1 == Interval(1.0, 2.0)
@test Intv2 == Interval(-100.0, 100.0)

# tests outer rounding
rnd_chk = Interval(-1.0,1.0)
rnd_chk_out = outer_rnd(rnd_chk)
@test rnd_chk == rnd_chk_out
=#

end
