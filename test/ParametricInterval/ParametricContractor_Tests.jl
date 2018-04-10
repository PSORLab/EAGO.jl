module ParametricContractor

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

# sets options for contractor
opt1 = Any[100    #Int64: Number of iterations
       1.0E-30 #Float64: Tolerance for equality of
       1.0E-6 #Float64: Add Interval(1E-8,1E8) to add to M[i,i] when
              #         processing extended interval division.
      ]

# Test Problem #1 Out of place
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
krawczykCW1 = PI_KrawczykCW(Z1,P1,hj1,h1,opt1,Eflag,Iflag)

@test -1E-4 <= newtonGS1[1][1].lo + 1.04243 <= 1E-4
@test -1E-4 <= newtonGS1[1][1].hi + 0.492759 <= 1E-4
@test -1E-4 <= newtonGS1[1][2].lo - 0.0473789 <= 1E-4
@test -1E-4 <= newtonGS1[1][2].hi - 0.208486 <= 1E-4
@test -1E-4 <= krawczykCW1[1][1].lo + 1.04243 <= 1E-4
@test -1E-4 <= krawczykCW1[1][1].hi + 0.492759 <= 1E-4
@test -1E-4 <= krawczykCW1[1][2].lo - 0.0473789 <= 1E-4
@test -1E-4 <= krawczykCW1[1][2].hi - 0.208486 <= 1E-4


# Test Problem #1 In place
function h1!(h,z,p)
   h[1] = z[1]^2+z[2]^2+p[1]*z[1]+4
   h[2] = z[1]+p[2]*z[2]
end

function hj1!(h,z,p)
   h[1,1] = (2*z[1]+p[1])
   h[1,2] = 2*z[2]
   h[2,1] = 1
   h[2,2] = p[2]
end

Eflag = false
Iflag = false
eDflag = false
newtonGS2 = PIn_NewtonGS(Z1,P1,hj1!,h1!,opt1,Eflag,Iflag,eDflag)
krawczykCW2 = PIn_KrawczykCW(Z1,P1,hj1!,h1!,opt1,Eflag,Iflag)

@test -1E-4 <= newtonGS2[1][1].lo + 1.04243 <= 1E-4
@test -1E-4 <= newtonGS2[1][1].hi + 0.492759 <= 1E-4
@test -1E-4 <= newtonGS2[1][2].lo - 0.0473789 <= 1E-4
@test -1E-4 <= newtonGS2[1][2].hi - 0.208486 <= 1E-4
@test -1E-4 <= krawczykCW2[1][1].lo + 1.04243 <= 1E-4
@test -1E-4 <= krawczykCW2[1][1].hi + 0.492759 <= 1E-4
@test -1E-4 <= krawczykCW2[1][2].lo - 0.0473789 <= 1E-4
@test -1E-4 <= krawczykCW2[1][2].hi - 0.208486 <= 1E-4

# Starts testing problem # 2
P2 = [Interval(6.0,9.0)]
Z2 = [Interval(-1.78, -0.1)]

# Test Problem #2 Out of Place
h2(x,p) = [x[1]^2 + p[1]*x[1] + 4.0]
hj2(x,p) = [2*x[1]+p[1]]
Eflag = false
Iflag = false
eDflag = false
newtonGS3 = PI_NewtonGS(Z2,P2,hj2,h2,opt1,Eflag,Iflag,eDflag)
krawczykCW3 = PI_KrawczykCW(Z2,P2,hj2,h2,opt1,Eflag,Iflag)
@test -0.772413 - 1E-4 <= krawczykCW3[1][1].lo <= -0.772413 + 1E-4
@test -0.383299 - 1E-4 <= krawczykCW3[1][1].hi <= -0.383299 + 1E-4
@test krawczykCW3[2] == false
@test krawczykCW3[3] == true
@test -0.769225 - 1E-4 <= newtonGS3[1][1].lo <= -0.769225 + 1E-4
@test -0.461725 - 1E-4 <= newtonGS3[1][1].hi <= -0.461725 + 1E-4
@test newtonGS3[3] == false
@test newtonGS3[4] == true

# Test Problem #2 Inplace
function h2!(hout,x,p)
   hout[:] = [x[1]^2 + p[1]*x[1] + 4.0]
end
function hj2!(hout,x,p)
   hout[:] = [2*x[1]+p[1]]
end
Eflag = false
Iflag = false
eDflag = false
newtonGS4 = PIn_NewtonGS(Z2,P2,hj2!,h2!,opt1,Eflag,Iflag,eDflag)
krawczykCW4 = PIn_KrawczykCW(Z2,P2,hj2!,h2!,opt1,Eflag,Iflag)
@test -0.772413 - 1E-4 <= krawczykCW4[1][1].lo <= -0.772413 + 1E-4
@test -0.383299 - 1E-4 <= krawczykCW4[1][1].hi <= -0.383299 + 1E-4
@test krawczykCW4[2] == false
@test krawczykCW4[3] == true
@test -0.769225 - 1E-4 <= newtonGS4[1][1].lo <= -0.769225 + 1E-4
@test -0.461725 - 1E-4 <= newtonGS4[1][1].hi <= -0.461725 + 1E-4
@test newtonGS4[3] == false
@test newtonGS4[4] == true

# Test Problem #3 (Exclusion)
Z2ex = [Interval(-30, -25)]
Eflag = false
Iflag = false
eDflag = false
newtonGS5 = PI_NewtonGS(Z2ex,P2,hj2,h2,opt1,Eflag,Iflag,eDflag)
krawczykCW5 = PI_KrawczykCW(Z2ex,P2,hj2,h2,opt1,Eflag,Iflag)
@test newtonGS5[3] == true
@test krawczykCW5[2] == true

Eflag = false
Iflag = false
eDflag = false
newtonGS6 = PIn_NewtonGS(Z2ex,P2,hj2!,h2!,opt1,Eflag,Iflag,eDflag)
krawczykCW6 = PIn_KrawczykCW(Z2ex,P2,hj2!,h2!,opt1,Eflag,Iflag)
@test newtonGS6[3] == true
@test krawczykCW6[2] == true

# Test Problem #4 (Extended Division)
Z3 = [Interval(-10, 10)]
newtonGS7 = PI_NewtonGS(Z3,P2,hj2,h2,opt1,Eflag,Iflag,eDflag)
newtonGS8 = PIn_NewtonGS(Z3,P2,hj2!,h2!,opt1,Eflag,Iflag,eDflag)
@test 0.285713 - 1E-4 <= newtonGS7[1][1].lo <= 0.285713 + 1E-4
@test 10 - 1E-4 <= newtonGS7[1][1].hi <= 10 + 1E-4
@test -10 - 1E-4 <= newtonGS7[2][1].lo <= -10 + 1E-4
@test -0.13793 - 1E-4 <= newtonGS7[2][1].hi <= -0.13793 + 1E-4
@test newtonGS7[3] == false
@test newtonGS7[4] == false
@test newtonGS7[5] == true
@test 0.285713 - 1E-4 <= newtonGS8[1][1].lo <= 0.285713 + 1E-4
@test 10 - 1E-4 <= newtonGS8[1][1].hi <= 10 + 1E-4
@test -10 - 1E-4 <= newtonGS8[2][1].lo <= -10 + 1E-4
@test -0.13793 - 1E-4 <= newtonGS8[2][1].hi <= -0.13793 + 1E-4
@test newtonGS8[3] == false
@test newtonGS8[4] == false
@test newtonGS8[5] == true

#=
Eflag = false
Iflag = false
eDflag = false
newtonGS9 = PId_NewtonGS(Z1,P1,hj1!,h1!,opt1,Eflag,Iflag,eDflag)
Eflag = false
Iflag = false
krawczykCW9 = PId_KrawczykCW(Z1,P1,hj1!,h1!,opt1,Eflag,Iflag)

@test -1E-4 <= newtonGS9[1][1].lo + 1.04243 <= 1E-4
@test -1E-4 <= newtonGS9[1][1].hi + 0.492759 <= 1E-4
@test -1E-4 <= newtonGS9[1][2].lo - 0.0473789 <= 1E-4
@test -1E-4 <= newtonGS9[1][2].hi - 0.208486 <= 1E-4
@test -1E-4 <= krawczykCW9[1][1].lo + 1.04243 <= 1E-4
@test -1E-4 <= krawczykCW9[1][1].hi + 0.492759 <= 1E-4
@test -1E-4 <= krawczykCW9[1][2].lo - 0.0473789 <= 1E-4
@test -1E-4 <= krawczykCW9[1][2].hi - 0.208486 <= 1E-4
=#

end
