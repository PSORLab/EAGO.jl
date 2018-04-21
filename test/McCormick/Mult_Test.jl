module Mult_Test

using Compat
using Compat.Test
using IntervalArithmetic
using StaticArrays
using EAGO

function XaboutY(x,y,tol)
    return abs(x-y) <= tol
end

################################################################################
################################################################################
##############      Testing for Nonsmooth Standard Mult           ##############
################################################################################
################################################################################
EAGO.set_diff_relax(0)

################################################################################
################### Test Nonsmooth Zero in Both Case (Failing)   ###############
################################################################################
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = SVector{2,Interval{Float64}}([Interval(-2.0,1.0);Interval(-1.0,2.0)])
mBox = mid.(xIBox)
X = SMCg{2,Interval{Float64},Float64}(0.0,0.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Interval{Float64},Float64}(1.0,1.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == 2.0
@test out.cv == -1.0
@test out.cc_grad[1] == 2.0
@test out.cc_grad[2] == -2.0
@test out.cv_grad[1] == 2.0
@test out.cv_grad[2] == 1.0
@test out.Intv.lo == -4.0
@test out.Intv.hi == 2.0


################################################################################
###################### Test Nonsmooth X1.l>0   (Passing)  ######################
################################################################################
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = SVector{2,Interval{Float64}}([Interval(1.0,5.0);Interval(-1.0,2.0)])
mBox = mid.(xIBox)
X = SMCg{2,Interval{Float64},Float64}(3.0,3.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Interval{Float64},Float64}(1.0,1.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == 5.0
@test out.cv == 1.0
@test out.cc_grad[1] == 2.0
@test out.cc_grad[2] == 1.0
@test out.cv_grad[1] == 2.0
@test out.cv_grad[2] == 5.0
@test out.Intv.lo == -5.0
@test out.Intv.hi == 10.0

################################################################################
############## Test Nonsmooth X1.h<0  &&  X2.l>0 (Passing)  ######################
################################################################################

a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = SVector{2,Interval{Float64}}([Interval(-6.0,-2.0);Interval(1.0,3.0)])
mBox = mid.(xIBox)
X = SMCg{2,Interval{Float64},Float64}(-4.0,-4.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Interval{Float64},Float64}(2.0,2.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == -6.0
@test out.cv == -10.0
@test out.cc_grad[1] == 1.0
@test out.cc_grad[2] == -2.0
@test out.cv_grad[1] == 3.0
@test out.cv_grad[2] == -2.0
@test out.Intv.lo == -18.0
@test out.Intv.hi == -2.0

################################################################################
############## Test Nonsmooth X1.h<0  &&  X2.h<0 (Passing)  ######################
################################################################################

a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = SVector{2,Interval{Float64}}([Interval(-6.0,-2.0);Interval(-7.0,-3.0)])
mBox = mid.(xIBox)
X = SMCg{2,Interval{Float64},Float64}(-4.0,-4.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Interval{Float64},Float64}(-5.0,-5.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == 24.0
@test out.cv == 16.0
@test out.cc_grad[1] == -7.0
@test out.cc_grad[2] == -2.0
@test out.cv_grad[1] == -3.0
@test out.cv_grad[2] == -2.0
@test out.Intv.lo == 6.0
@test out.Intv.hi == 42.0

################################################################################
############## Test Nonsmooth X1.h<0  &&  0 in X2 (Passing)  ###################
################################################################################

a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = SVector{2,Interval{Float64}}([Interval(-6.0,-2.0);Interval(-7.0,4.0)])
mBox = mid.(xIBox)
X = SMCg{2,Interval{Float64},Float64}(-4.0,-4.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Interval{Float64},Float64}(-5.0,-5.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == 24.0
@test out.cv == 16.0
@test out.cc_grad[1] == -7.0
@test out.cc_grad[2] == -2.0
@test out.cv_grad[1] == -7.0
@test out.cv_grad[2] == -6.0
@test out.Intv.lo == -24.0
@test out.Intv.hi == 42.0

################################################################################
############## Test Nonsmooth 0 in X1  &&  X2.l > 0 ()  #################
################################################################################
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,4.0);Interval(1.0,4.0)])
mBox = mid.(xIBox)
X = SMCg{2,Interval{Float64},Float64}(-2.0,-2.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Interval{Float64},Float64}(3.0,3.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == -5.0
@test out.cv == -8.0
@test out.cc_grad[1] == 4.0
@test out.cc_grad[2] == -3.0
@test out.cv_grad[1] == 1.0
@test out.cv_grad[2] == -3.0
@test out.Intv.lo == -12.0
@test out.Intv.hi == 16.0

################################################################################
############## Test Nonsmooth 0 in X1  &&  X2.h < 0 ()         #################
################################################################################
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,4.0);Interval(-5.0,-3.0)])
mBox = mid.(xIBox)
X = SMCg{2,Interval{Float64},Float64}(-2.0,-2.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Interval{Float64},Float64}(-4.0,-4.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == 9.0
@test out.cv == 7.0
@test out.cc_grad[1] == -3.0
@test out.cc_grad[2] == -3.0
@test out.cv_grad[1] == -5.0
@test out.cv_grad[2] == -3.0
@test out.Intv.lo == -20.0
@test out.Intv.hi == 15.0


################################################################################
################################################################################
##############      Testing for Smooth Standard Mult           ##############
################################################################################
################################################################################
EAGO.set_diff_relax(1)

seed1 = seed_g(Float64,1,2)
seed2 = seed_g(Float64,2,2)
x1 = SMCg{2,Interval{Float64},Float64}(0.0,0.0,seed1,seed1,Interval(-200.0,200.0),false, SVector{2,Interval{Float64}}([Interval(-200.0,200.0),Interval(0.0,400.0)]),SVector{2,Float64}([0.0,200.0]))
y1 = SMCg{2,Interval{Float64},Float64}(200.0,200.0,seed2,seed2,Interval(0.0,400.0),false, SVector{2,Interval{Float64}}([Interval(-200.0,200.0),Interval(0.0,400.0)]),SVector{2,Float64}([0.0,200.0]))
z1 = x1*y1
@test XaboutY(z1.cc,40000,1E-4)
@test XaboutY(z1.cv,-40000,1E-4)

x2 = SMCg{2,Interval{Float64},Float64}(170.0,170.0,seed1,seed1,Interval(100.0,240.0),false, SVector{2,Interval{Float64}}([Interval(100.0,240.0),Interval(100.0,400.0)]),SVector{2,Float64}([170.0,250.0]))
y2 = SMCg{2,Interval{Float64},Float64}(250.0,250.0,seed2,seed2,Interval(100.0,400.0),false, SVector{2,Interval{Float64}}([Interval(100.0,240.0),Interval(100.0,400.0)]),SVector{2,Float64}([170.0,250.0]))
z2 = x2*y2
@test XaboutY(z2.cc,53000,1E-4)
@test XaboutY(z2.cv,32000,1E-4)

x3 = SMCg{2,Interval{Float64},Float64}(-200.0,-200.0,seed1,seed1,Interval(-300.0,-100.0),false, SVector{2,Interval{Float64}}([Interval(-300.0,-100.0),Interval(-400.0,-200.0)]),SVector{2,Float64}([-200.0,-300.0]))
y3 = SMCg{2,Interval{Float64},Float64}(-300.0,-300.0,seed2,seed2,Interval(-400.0,-200.0),false, SVector{2,Interval{Float64}}([Interval(-300.0,-100.0),Interval(-400.0,-200.0)]),SVector{2,Float64}([-200.0,-300.0]))
z3 = x3*y3
@test XaboutY(z3.cc,70000,1E-4)
@test XaboutY(z3.cv,50000,1E-4)

x4 = SMCg{2,Interval{Float64},Float64}(150.0,150.0,seed1,seed1,Interval(100.0,200.0),false, SVector{2,Interval{Float64}}([Interval(100.0,200.0),Interval(-500.0,-100.0)]),SVector{2,Float64}([150.0,-300.0]))
y4 = SMCg{2,Interval{Float64},Float64}(-250.0,-250.0,seed2,seed2,Interval(-500.0,-100.0),false, SVector{2,Interval{Float64}}([Interval(100.0,200.0),Interval(-500.0,-100.0)]),SVector{2,Float64}([150.0,-300.0]))
z4 = x4*y4
@test XaboutY(z4.cc,-30000,1E-3)
@test XaboutY(z4.cv,-47460.9375,1E-3)

x5 = SMCg{2,Interval{Float64},Float64}(-150.0,-150.0,seed1,seed1,Interval(-200.0,-100.0),false, SVector{2,Interval{Float64}}([Interval(-200.0,-100.0),Interval(200.0,400.0)]),SVector{2,Float64}([-150.0,300.0]))
y5 = SMCg{2,Interval{Float64},Float64}(300.0,300.0,seed2,seed2,Interval(200.0,400.0),false, SVector{2,Interval{Float64}}([Interval(-200.0,-100.0),Interval(200.0,400.0)]),SVector{2,Float64}([-150.0,300.0]))
z5 = x5*y5
@test XaboutY(z5.cv,-50000,1E-4)
@test XaboutY(z5.cc,-40000,1E-4)

end
