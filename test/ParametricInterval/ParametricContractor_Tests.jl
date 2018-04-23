module ParametricContractor

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

opt1 = PIntvParams(:Dense,:Newton,1E-30,1E-6,Int64(2),Int64(2),100)
opt2 = PIntvParams(:Dense,:Krawczyk,1E-30,1E-6,Int64(2),Int64(2),100)

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
newtonGS1 = Param_Intv_Contractor(h1,hj1,Z1,P1,Eflag,Iflag,eDflag,opt1)
krawczykCW1 = Param_Intv_Contractor(h1,hj1,Z1,P1,Eflag,Iflag,eDflag,opt2)

@test -1E-4 <= newtonGS1[1][1].lo + 1.04243 <= 1E-4
@test -1E-4 <= newtonGS1[1][1].hi + 0.492759 <= 1E-4
@test -1E-4 <= newtonGS1[1][2].lo - 0.0473789 <= 1E-4
@test -1E-4 <= newtonGS1[1][2].hi - 0.208486 <= 1E-4
@test -1E-4 <= krawczykCW1[1][1].lo + 1.04243 <= 1E-4
@test -1E-4 <= krawczykCW1[1][1].hi + 0.492759 <= 1E-4
@test -1E-4 <= krawczykCW1[1][2].lo - 0.0473789 <= 1E-4
@test -1E-4 <= krawczykCW1[1][2].hi - 0.208486 <= 1E-4

# Starts testing problem # 2
opt1a = PIntvParams(:Dense,:Newton,1E-30,1E-6,1,1,100)
opt2a = PIntvParams(:Dense,:Krawczyk,1E-30,1E-6,1,1,100)
P2 = [Interval(6.0,9.0)]
Z2 = [Interval(-1.78, -0.1)]

# Test Problem #2 Out of Place
h2(x,p) = [x[1]^2 + p[1]*x[1] + 4.0]
hj2(x,p) = [2*x[1]+p[1]]
Eflag = false
Iflag = false
eDflag = false
Z2c = copy(Z2)
newtonGS3 = Param_Intv_Contractor(h2,hj2,Z2,P2,Eflag,Iflag,eDflag,opt1a)
krawczykCW3 = Param_Intv_Contractor(h2,hj2,Z2c,P2,Eflag,Iflag,eDflag,opt2a)
@test -0.772413 - 1E-4 <= krawczykCW3[1][1].lo <= -0.772413 + 1E-4
@test -0.383299 - 1E-4 <= krawczykCW3[1][1].hi <= -0.383299 + 1E-4
@test krawczykCW3[3] == false
@test krawczykCW3[4] == true
@test -0.769225 - 1E-4 <= newtonGS3[1][1].lo <= -0.769225 + 1E-4
@test -0.461725 - 1E-4 <= newtonGS3[1][1].hi <= -0.461725 + 1E-4
@test newtonGS3[3] == false
@test newtonGS3[4] == true

# Test Problem #3 (Exclusion)
Z2ex = [Interval(-30, -25)]
Eflag = false
Iflag = false
eDflag = false
Z2exc = copy(Z2ex)
Z2exc1 = copy(Z2ex)
Z2exc2 = copy(Z2ex)
newtonGS5 = Param_Intv_Contractor(h2,hj2,Z2ex,P2,Eflag,Iflag,eDflag,opt1a)
krawczykCW5 = Param_Intv_Contractor(h2,hj2,Z2exc,P2,Eflag,Iflag,eDflag,opt2a)
@test newtonGS5[3] == true
@test krawczykCW5[3] == true

Eflag = false
Iflag = false
eDflag = false
newtonGS6 = Param_Intv_Contractor(h2,hj2,Z2exc1,P2,Eflag,Iflag,eDflag,opt1a)
krawczykCW6 = Param_Intv_Contractor(h2,hj2,Z2exc2,P2,Eflag,Iflag,eDflag,opt2a)
@test newtonGS6[3] == true
@test krawczykCW6[3] == true

# Test Problem #4 (Extended Division)
Z3 = [Interval(-10, 10)]
Z3c = copy(Z3)
newtonGS7 = Param_Intv_Contractor(h2,hj2,Z3,P2,Eflag,Iflag,eDflag,opt1a)
newtonGS8 = Param_Intv_Contractor(h2,hj2,Z3c,P2,Eflag,Iflag,eDflag,opt1a)
@test 0.285713 - 1E-4 <= newtonGS7[1][1].lo <= 0.285713 + 1E-4
@test 0.548861585037655 - 1E-4 <= newtonGS7[1][1].hi <= 0.548861585037655 + 1E-4
@test 0.28571374999999993 - 1E-4 <= newtonGS7[2][1].lo <= 0.28571374999999993 + 1E-4
@test 0.548861585037655 - 1E-4 <= newtonGS7[2][1].hi <= 0.548861585037655 + 1E-4
@test newtonGS7[3] == true
@test newtonGS7[4] == false
@test newtonGS7[5] == false
@test 0.285713 - 1E-4 <= newtonGS8[1][1].lo <= 0.285713 + 1E-4
@test 0.548861585037655 - 1E-4 <= newtonGS8[1][1].hi <= 0.548861585037655 + 1E-4
@test 0.28571374999999993 - 1E-4 <= newtonGS8[2][1].lo <= 0.28571374999999993 + 1E-4
@test 0.548861585037655 - 1E-4 <= newtonGS8[2][1].hi <= 0.548861585037655 + 1E-4
@test newtonGS8[3] == true
@test newtonGS8[4] == false
@test newtonGS8[5] == false

@testset "Test Contractor Warnings" begin
opt3a = PIntvParams(:asdada,:Krawczyk,1E-30,1E-6,1,1,100)
opt3b = PIntvParams(:asdada,:Newton,1E-30,1E-6,1,1,100)
opt3c = PIntvParams(:asdada,:asdasdasd,1E-30,1E-6,1,1,100)

P1 = [Interval(5.0,7.0),Interval(5.0,7.0)]
Z1 = [Interval(-1.5, 0.0),Interval(0.0, 0.5)]
h1(z,p) = [z[1]^2+z[2]^2+p[1]*z[1]+4;
           z[1]+p[2]*z[2]]
hj1(z,p) = [(2*z[1]+p[1]) (2*z[2]);
              1              p[2]]
Eflag = false
Iflag = false
eDflag = false
@test_throws ErrorException Param_Intv_Contractor(h1,hj1,Z1,P1,Eflag,Iflag,eDflag,opt3a)
@test_throws ErrorException Param_Intv_Contractor(h1,hj1,Z1,P1,Eflag,Iflag,eDflag,opt3b)
@test_throws ErrorException Param_Intv_Contractor(h1,hj1,Z1,P1,Eflag,Iflag,eDflag,opt3c)
end
end
