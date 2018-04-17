#workspace()
using BenchmarkTools
using IntervalArithmetic
using EAGO

opt1 = PIntvParams(:Dense,:Krawczyk,1E-30,1E-6,2,0,100)
opt2 = PIntvParams(:Dense,:Newton,1E-30,1E-6,2,0,100)


P1 = [MCInterval(5.0,7.0),MCInterval(5.0,7.0)]
Z1 = [MCInterval(-1.5, 0.0),MCInterval(0.0, 0.5)]
Z1c = copy(Z1)
h1(z,p) = [z[1]^2+z[2]^2+p[1]*z[1]+4.0;
           z[1]+p[2]*z[2]]
hj1(z,p) = [(z[1]+p[1]) (z[2]);
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

#@benchmark Param_Intv_Contractor($h1,$hj1,$Z1,$P1,$Eflag,$Iflag,$eDflag,$opt1)
# @code_warntype Param_Intv_Contractor(h1,hj1,Z1,P1,Eflag,Iflag,eDflag,opt1)

#=
# sets options for contractor
opt1 = Any[50    #Int64: Number of iterations
       1.0E-30 #Float64: Tolerance for equality of
       1.0E-10 #Float64: Add Interval(1E-8,1E8) to add to M[i,i] when
              #         processing extended interval division.
      ]
################################################################################
#########################  Test Problem #1 #####################################
################################################################################

P1 = [Interval(5.0,7.0),Interval(5.0,7.0)]
Z1 = [Interval(-1.5, 0.0),Interval(0.0, 0.5)]
h1(z,p) = [z[1]^2+z[2]^2+p[1]*z[1]+4;
           z[1]+p[2]*z[2]]
hj1(z,p) = [(2*z[1]+p[1]) (2*z[2]);
              1              p[2]]

Eflag = false
Iflag = false
eDflag = false
newtonGS1 = PI_NewtonGS(Z1,P1,hj1,h1,opt1,Eflag,Iflag,eDflag)
Eflag = false
Iflag = false
krawczykCW1 = PI_KrawczykCW(Z1,P1,hj1,h1,opt1,Eflag,Iflag)

function h1!(hh,z,p)
       hh[1] = z[1]^2+z[2]^2+p[1]*z[1]+4
       hh[2] = z[1]+p[2]*z[2]
end
function hj1!(hh,z,p)
       hh[1,1] = 2*z[1]+p[1]
       hh[1,2] = 2*z[2]
       hh[2,1] = 1
       hh[2,2] = p[2]
end

Eflag = false
Iflag = false
eDflag = false
InnewtonGS1 = PIn_NewtonGS(Z1,P1,hj1!,h1!,opt1,Eflag,Iflag,eDflag)
Eflag = false
Iflag = false
InkrawczykCW1 = PIn_KrawczykCW(Z1,P1,hj1!,h1!,opt1,Eflag,Iflag)

################################################################################
#########################  Test Problem #2 #####################################
################################################################################

function LBD_func(i)
          YL = [-30, -30, -30, 1800, 900]
          return YL[i]
end
function UBD_func(i)
          YU = [30, 30, 30, 2200, 1100]
          return YU[i]
end

P2 = [Interval(1800,2200),Interval(900,1100)]
Z2 = [Interval(-30, 30),Interval(-30, 30),Interval(-30, 30)]
function h2(x,p)
       return [(3.25-x[1])/p[1]-x[3]; x[1]/p[2]-x[3]; x[2]-(x[1]^2)/(1+x[1]^2)]
end
function hj2(x,p)
       return [-1/p[1] 0 -1;
       1/p[2] 0 -1;
       -2*x[1]/(1+x[1]^2)^2 1 0]
end

Eflag = false
Iflag = false
eDflag = false
newtonGS2 = PI_NewtonGS(Z2,P2,hj2,h2,opt1,Eflag,Iflag,eDflag)
P2 = [Interval(1800,2200),Interval(900,1100)]
Z2 = [Interval(-30, 30),Interval(-30, 30),Interval(-30, 30)]
Eflag = false
Iflag = false
krawczykCW2 = PI_KrawczykCW(Z2,P2,hj2,h2,opt1,Eflag,Iflag)

function h2!(hout,x,p)
       hout[1] = (3.25-x[1])/p[1]-x[3]
       hout[2] = x[1]/p[2]-x[3]
       hout[3] = x[2]-(x[1]^2)/(1+x[1]^2)
end
function hj2!(hout,x,p)
       hout[1,1] = -1/p[1]
       hout[1,2] = 0
       hout[1,3] = -1
       hout[2,1] = 1/p[2]
       hout[2,2] = 0
       hout[2,3] = -1
       hout[3,1] = -2*x[1]/(1+x[1]^2)^2
       hout[3,2] = 1
       hout[3,3] = 0
end
P2 = [Interval(1800,2200),Interval(900,1100)]
Z2 = [Interval(-30, 30),Interval(-30, 30),Interval(-30, 30)]
Eflag = false
Iflag = false
eDflag = false
InnewtonGS2 = PIn_NewtonGS(Z2,P2,hj2!,h2!,opt1,Eflag,Iflag,eDflag)
P2 = [Interval(1800,2200),Interval(900,1100)]
Z2 = [Interval(-30, 30),Interval(-30, 30),Interval(-30, 30)]
Eflag = false
Iflag = false
InkrawczykCW2 = PIn_KrawczykCW(Z2,P2,hj2!,h2!,opt1,Eflag,Iflag)





# sets options for contractor
opt1 = Any[50    #Int64: Number of iterations
       1.0E-30 #Float64: Tolerance for equality of
       1.0E-10 #Float64: Add Interval(1E-8,1E8) to add to M[i,i] when
              #         processing extended interval division.
      ]
################################################################################
#########################  Test Problem #1 #####################################
################################################################################

P1 = [MCInterval(5.0,7.0),MCInterval(5.0,7.0)]
Z1 = [MCInterval(-1.5, 0.0),MCInterval(0.0, 0.5)]
h1(z,p) = [z[1]^2+z[2]^2+p[1]*z[1]+4;
           z[1]+p[2]*z[2]]
hj1(z,p) = [(2*z[1]+p[1]) (2*z[2]);
              1              p[2]]

Eflag = false
Iflag = false
eDflag = false
newtonGS1 = PI_NewtonGS(Z1,P1,hj1,h1,opt1,Eflag,Iflag,eDflag)
Eflag = false
Iflag = false
krawczykCW1 = PI_KrawczykCW(Z1,P1,hj1,h1,opt1,Eflag,Iflag)

function h1!(hh,z,p)
       hh[1] = z[1]^2+z[2]^2+p[1]*z[1]+4
       hh[2] = z[1]+p[2]*z[2]
end
function hj1!(hh,z,p)
       hh[1,1] = 2*z[1]+p[1]
       hh[1,2] = 2*z[2]
       hh[2,1] = 1
       hh[2,2] = p[2]
end

Eflag = false
Iflag = false
eDflag = false
InnewtonGS1 = PIn_NewtonGS(Z1,P1,hj1!,h1!,opt1,Eflag,Iflag,eDflag)
Eflag = false
Iflag = false
InkrawczykCW1 = PIn_KrawczykCW(Z1,P1,hj1!,h1!,opt1,Eflag,Iflag)

################################################################################
#########################  Test Problem #2 #####################################
################################################################################

function LBD_func(i)
          YL = [-30, -30, -30, 1800, 900]
          return YL[i]
end
function UBD_func(i)
          YU = [30, 30, 30, 2200, 1100]
          return YU[i]
end

P2 = [MCInterval(1800,2200),MCInterval(900,1100)]
Z2 = [MCInterval(-30, 30),MCInterval(-30, 30),MCInterval(-30, 30)]
function h2(x,p)
       return [(3.25-x[1])/p[1]-x[3]; x[1]/p[2]-x[3]; x[2]-(x[1]^2)/(1+x[1]^2)]
end
function hj2(x,p)
       return [-1/p[1] 0 -1;
       1/p[2] 0 -1;
       -2*x[1]/(1+x[1]^2)^2 1 0]
end

Eflag = false
Iflag = false
eDflag = false
newtonGS2 = PI_NewtonGS(Z2,P2,hj2,h2,opt1,Eflag,Iflag,eDflag)
P2 = [MCInterval(1800,2200),MCInterval(900,1100)]
Z2 = [MCInterval(-30, 30),MCInterval(-30, 30),MCInterval(-30, 30)]
Eflag = false
Iflag = false
krawczykCW2 = PI_KrawczykCW(Z2,P2,hj2,h2,opt1,Eflag,Iflag)

function h2!(hout,x,p)
       hout[1] = (3.25-x[1])/p[1]-x[3]
       hout[2] = x[1]/p[2]-x[3]
       hout[3] = x[2]-(x[1]^2)/(1+x[1]^2)
end
function hj2!(hout,x,p)
       hout[1,1] = -1/p[1]
       hout[1,2] = 0
       hout[1,3] = -1
       hout[2,1] = 1/p[2]
       hout[2,2] = 0
       hout[2,3] = -1
       hout[3,1] = -2*x[1]/(1+x[1]^2)^2
       hout[3,2] = 1
       hout[3,3] = 0
end
P2 = [MCInterval(1800,2200),MCInterval(900,1100)]
Z2 = [MCInterval(-30, 30),MCInterval(-30, 30),MCInterval(-30, 30)]
Eflag = false
Iflag = false
eDflag = false
InnewtonGS2 = PIn_NewtonGS(Z2,P2,hj2!,h2!,opt1,Eflag,Iflag,eDflag)
P2 = [MCInterval(1800,2200),MCInterval(900,1100)]
Z2 = [MCInterval(-30, 30),MCInterval(-30, 30),MCInterval(-30, 30)]
Eflag = false
Iflag = false
InkrawczykCW2 = PIn_KrawczykCW(Z2,P2,hj2!,h2!,opt1,Eflag,Iflag)
=#
