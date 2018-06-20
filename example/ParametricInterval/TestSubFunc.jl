
# Will be moved to test functions

#=
workspace()

using IntervalArithmetic
using EAGO

# Miranda Passed
X = [Interval(1,3),Interval(1,3)]
P = [Interval(1,3),Interval(1,3)]
h(z,p) = [p[1]-z[1]; p[1]+z[1]]
h1(z,p) = [p[1]-z[1]; p[1]-z[1]]
test_mir1 = Miranda(h,X,P)
test_mir2 = Miranda(h1,X,P)
out_mir = h(X,P)

# Miranda Exclusion Test Passed
Eflag_in = false
ex_mir1 = MirandaExc(h,X,P,Eflag_in)
ex_mir2 = MirandaExc(h1,X,P,Eflag_in)
ex_mir_flag_out = Eflag_in

#
Y1 = [Interval(1,3),Interval(1,3)]
Y2 = [Interval(1.5,2.5),Interval(2.5,2.5)]
Y3 = [Interval(-1,2),Interval(-2,2)]
Y4 = [Interval(2,4),Interval(2,4)]
Y5 = [Interval(1.5,2.5),Interval(2.3,2.5)]
Y6 = [Interval(-1,4),Interval(-1,4)]
Y7 = [Interval(1.5,2.5),emptyinterval()]
Y8 = [emptyinterval(),emptyinterval()]
Y9 = [Interval(4,5),Interval(4,5)]
Y10 = [Interval(-2,-1),Interval(-2,-1)]

sXY_flag1 = Strict_XinY(Y1,X)
sXY_flag2 = Strict_XinY(Y2,X)
sXY_flag3 = Strict_XinY(Y3,X)
sXY_flag4 = Strict_XinY(Y4,X)
sXY_flag5 = Strict_XinY(Y5,X)
sXY_flag6 = Strict_XinY(Y6,X)
sXY_flag7 = Strict_XinY(Y7,X)
sXY_flag8 = Strict_XinY(Y8,X)
sXY_flag9 = Strict_XinY(Y9,X)
sXY_flag10 = Strict_XinY(Y10,X)

#=
partialIncTop(h,X,P,PIflag,incHigh)

partialIncBot(h,X,P,PIflag,incLow)

Strict_XinY(X,Y)

isEqual(X1,X2,atol)

extDivide(A,B,C)

extProcess(N,X,Mii,S1,S2,B,rtol)
=#
