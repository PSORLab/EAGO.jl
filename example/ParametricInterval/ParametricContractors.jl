#workspace()
using BenchmarkTools
using IntervalArithmetic
using ForwardDiff
using EAGO


# Test simple contractor cases from Stuber2011
opt1 = PIntvParams(:Dense,:Krawczyk,1E-30,1E-6,2,0,100)
opt2 = PIntvParams(:Dense,:Newton,1E-30,1E-6,2,0,100)

P1 = [MCInterval(5.0,7.0),MCInterval(5.0,7.0)]
Z1 = [MCInterval(-1.5, 0.0),MCInterval(0.0, 0.5)]
Z1c = copy(Z1)
P1a= [Interval(5.0,7.0),Interval(5.0,7.0)]
Z1a = [Interval(-1.5, 0.0),Interval(0.0, 0.5)]
Z1ca = copy(Z1a)

h1(z,p) = [z[1]^2+z[2]^2+p[1]*z[1]+4.0;
           z[1]+p[2]*z[2]]
hj1(z,p) = [(2.0*z[1]+p[1]) (2.0*z[2]);
              one(p[1])              p[2]]
Eflag = false
Iflag = false
eDflag = false
H = h1(mid.(Z1),P1)
J = hj1(Z1,P1)
Y = mid.(J)
F = lufact(Y)
YH = F\H
YJ = F\J

output1 = Param_Intv_Contractor(h1,hj1,Z1,P1,Eflag,Iflag,eDflag,opt1)
output2 = Param_Intv_Contractor(h1,hj1,Z1c,P1,Eflag,Iflag,eDflag,opt2) # mutates Z1 in place
output1a = Param_Intv_Contractor(h1,hj1,Z1a,P1a,Eflag,Iflag,eDflag,opt1)
output2a = Param_Intv_Contractor(h1,hj1,Z1ca,P1a,Eflag,Iflag,eDflag,opt2)

# Test simple contractor cases from Stuber2015 (1-D example)
function hcos(z,y)
    [z[1]-(y[1]-(y[1]^3)/6+(y[1]^5)/120)/sqrt(z[1])-y[2]]
end
function hjcos(z,y)
    [(y[1]^5-20*y[1]^3+120*y[1])/(240*sqrt(z[1]^3))+1]
end
yBnds = [Interval(0.5,8.0),Interval(80,120)]
xBnds = [Interval(50.8,200.9)]
xBndsc = [Interval(50.8,200.9)]
opt1a = PIntvParams(:Dense,:Krawczyk,1E-30,1E-6,1,0,100)
opt2a = PIntvParams(:Dense,:Newton,1E-30,1E-6,1,0,100)
Eflag = false
Iflag = false
eDflag = false
output3a = Param_Intv_Contractor(hcos,hjcos,xBnds,yBnds,Eflag,Iflag,eDflag,opt1a)
Eflag = false
Iflag = false
eDflag = false
output3b = Param_Intv_Contractor(hcos,hjcos,xBndsc,yBnds,Eflag,Iflag,eDflag,opt2a)

# Stuber SIP contractor example 2
g2(x,p) = x[1] + cos(p[1]-80/90) - 80
function h2(y,x)
    [y[1]-(x[1]-(x[1]^3)/6+(x[1]^5)/120)/sqrt(y[1])-80]
end
function hj2(y,x)
    [1.0+(x[1]-(x[1]^3)/6+(x[1]^5)/120)/(2.0*sqrt(y[1]^3))]
end
Xlast = [Interval(80.2046, 80.5215)]
Plast = [Interval(2.60937, 2.84375)]
h2val = h2(Xlast,Plast)

Eflag = false
Iflag = false
eDflag = false
outputlasta = Param_Intv_Contractor(h2,hj2,Xlast,Plast,Eflag,Iflag,eDflag,opt1a)
Eflag = false
Iflag = false
eDflag = false
outputlastb = Param_Intv_Contractor(h2,hj2,Xlast,Plast,Eflag,Iflag,eDflag,opt2a)
