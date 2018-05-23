#workspace()

# Just a collection of tests that had been used while debugging
using EAGO
using IntervalArithmetic
using StaticArrays

F(X,Y) = X-(Y-pow(Y,3)/6+pow(Y,5)/120)/sqrt(X)-80; #pow((Y-3.5),3) #pow((Y-3.5),4) - 5.0*pow((Y-3.5),3) - 2.0*pow((Y-3.5),2) + 15.0*(Y-3.5)
EAGO.set_diff_relax(0)
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = SVector{2,Interval{Float64}}([Interval(3.0,7.0),Interval(3.0,7.0)])
mBox = mid.(xIBox)
xI1 = Interval(68.7999,109.35)
xI2 = Interval(109.349,149.901)
xI3 = Interval(78.9375,84.0063)
xI4 = Interval(73.8687,76.4032)

pI1 = Interval(0.5,8.0)
pI2 = Interval(0.5,8.0)
pI3 = Interval(0.5,4.25)
pI4 = Interval(6.125,8.0)

X1 = SMCg{2,Interval{Float64},Float64}(mid(xI1),mid(xI1),a,a,xI1,false,xIBox,mBox)
X2 = SMCg{2,Interval{Float64},Float64}(mid(xI2),mid(xI2),a,a,xI2,false,xIBox,mBox)
X3 = SMCg{2,Interval{Float64},Float64}(mid(xI3),mid(xI3),a,a,xI3,false,xIBox,mBox)
X4 = SMCg{2,Interval{Float64},Float64}(mid(xI4),mid(xI4),a,a,xI4,false,xIBox,mBox)

P1 = SMCg{2,Interval{Float64},Float64}(mid(pI1),mid(pI1),b,b,pI1,false,xIBox,mBox)
P2 = SMCg{2,Interval{Float64},Float64}(mid(pI2),mid(pI2),b,b,pI2,false,xIBox,mBox)
P3 = SMCg{2,Interval{Float64},Float64}(mid(pI3),mid(pI3),b,b,pI3,false,xIBox,mBox)
P4 = SMCg{2,Interval{Float64},Float64}(mid(pI4),mid(pI4),b,b,pI4,false,xIBox,mBox)

F1 = F(X1,P1)
F2 = F(X2,P2)
F3 = F(X3,P3)
F4 = F(X4,P4)

F1F(X,Y) = (Y-pow(Y,3)/6+pow(Y,5)/120)
F11 = F1F(X1,P1)
F12 = F1F(X2,P2)
#F5 = F(X5)
