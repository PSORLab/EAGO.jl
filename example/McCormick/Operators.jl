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
xIBox = SVector{2,Interval{Float64}}([xIntv1,xIntv1])
mBox = mid.(xIBox)
SMCg = SMCg{2,Interval{Float64},Float64}(x,x,a,a,xIntv1,false,xIBox,mBox)

# calculates relaxations & gradients of (exp(x1)*x1-x1^3)*sin(x1) for x1 = 2.0 on [1.0,3.0]
relaxed_f = (exp(SMCg)*SMCg-SMCg^3)*sin(SMCg)
println("concave relaxaxation value: ", relaxed_f.cc)
println("convex relaxaxation value: ", relaxed_f.cv)
println("concave (sub)gradient: ", relaxed_f.cc_grad)
println("convex (sub)gradient: ", relaxed_f.cv_grad)
println("Interval Bounds: ", relaxed_f.Intv)

SMCglist = []
