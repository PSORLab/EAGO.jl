workspace()

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
xIBox = [xIntv1;xIntv1]
mBox = mid.(xIBox)
SMCg = SMCg{2,Float64}(x,x,a,a,xIntv1,false)

# resets gradient to seed gradient with seed at j=2
grad(SMCg,2)
println("convex gradient: ", SMCg.cv_grad)
println("concave gradient: ", SMCg.cc_grad)

# sets gradients to zero
zgrad(SMCg)
println("convex gradient: ", SMCg.cv_grad)
println("concave gradient: ", SMCg.cc_grad)

# checks that mid3 returns midpoint
mid_return = mid3(9.0,-3.0,5.0)

# checks that mid_grad is working
mg_return1 = mid_grad(a, b, 1)
mg_return2 = mid_grad(a, b, 2)
mg_return3 = mid_grad(a, b, 3)

# generates line segment value and derivative
y0 = line_seg(2,1,0,6,10)
dy0 = dline_seg(2,1,0,6,10)

# converts float/integer to SMCg
promoted1 = promote_rule(EAGO.SMCg{3,Float64}, Int64)
promoted2 = promote_rule(EAGO.SMCg{3,Int64}, Float64)
promoted3 = promote_rule(EAGO.SMCg{3,Float64}, Interval{Float64})
promoted4 = promote_rule(EAGO.SMCg{3,Float64}, Real)

conv_float = convert(EAGO.SMCg{3,Float64},1.0)
conv_int = convert(EAGO.SMCg{3,Float64},3)

# test subgradient refinement
Intv1 = tighten_subgrad(2.0,1.0,b,b,Interval(-100.0,100.0),[Interval(-1.0,1.0)],[0.0])
Intv2 = tighten_subgrad(2.0,1.0,b,b,Interval(-100.0,100.0),[∅],[0.0])

# tests outer rounding
rnd_chk = Interval(-1.0,1.0)
rnd_chk_out = outer_rnd(rnd_chk)
