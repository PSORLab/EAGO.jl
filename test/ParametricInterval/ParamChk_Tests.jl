module ParamChk_Tests

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

# Tests strict inclusion procedure for interval vectors
Y = [Interval(0,5),Interval(0,5),Interval(0,5)]
X1 = [Interval(1,2),Interval(1,2),Interval(1,2)]
X2 = [Interval(1,2),Interval(-10,10),Interval(1,2)]
X3 = [Interval(1,2),Interval(1,2),Interval(0,5)]
flag1 = EAGO.Strict_XinY(X1,Y)
flag2 = EAGO.Strict_XinY(X2,Y)
flag3 = EAGO.Strict_XinY(X3,Y)
@test flag1 == true
@test flag2 == false
@test flag3 == false

# Checks extended division routine
A1 = Interval(0)
A2 = Interval(0,3)
A3 = Interval(-2,0)
A4 = Interval(-3,2)
ind1,B1,C1 = EAGO.extDivide(A1)
ind2,B2,C2 = EAGO.extDivide(A2)
ind3,B3,C3 = EAGO.extDivide(A3)
ind4,B4,C4 = EAGO.extDivide(A4)
@test ind1 == 0
@test ind2 == 1
@test ind3 == 2
@test ind4 == 3
@test B1 == Interval(-Inf,Inf)
@test 0.33333 - 1E-4 <= B2.lo <= 0.33333 + 1E-4
@test B2.hi == Inf
@test B3 == Interval(-Inf,-0.5)
@test B4.lo == -Inf
@test -0.33333 - 1E-4 <= B4.hi <= -0.33333 + 1E-4
@test C1 == Interval(-Inf,Inf)
@test C2 == Interval(Inf,Inf)
@test C3 == Interval(-Inf,-Inf)
@test C4 == Interval(0.5,Inf)

#param = EAGO.Param_Bisect_Opts()
#@test param.DAGflag == false

N =  Interval(-5,5)
X = Interval(-5,5)
Mii = Interval(-5,5)
S1 = Interval(-5,5)
S2 = Interval(-5,5)
B = Interval(-5,5)
rtol = 1E-4
indx1,box11,box12 = EAGO.extProcess(N,X,Mii,S1,S2,B,rtol)
Miib = Interval(0,5)
S1b = Interval(1,5)
S2b = Interval(1,5)
Bb = Interval(1,5)
indx2,box21,box22 = EAGO.extProcess(N,X,Miib,S1b,S2b,Bb,rtol)
Miic = Interval(-5,0)
S1c = Interval(1,5)
S2c = Interval(1,5)
Bc = Interval(1,5)
indx3,box31,box32 = EAGO.extProcess(N,X,Miic,S1c,S2c,Bc,rtol)
Miia = Interval(1,5)
S1a = Interval(1,5)
S2a = Interval(1,5)
Ba = Interval(1,5)
indx6,box61,box62 = EAGO.extProcess(N,X,Miia,S1a,S2a,Ba,rtol)
Miid = Interval(0,0)
S1d = Interval(1,5)
S2d = Interval(1,5)
Bd = Interval(1,5)
indx8,box81,box82 = EAGO.extProcess(N,X,Miid,S1d,S2d,Bd,rtol)

@test indx1 == 0
@test box11 == Interval(-Inf,Inf)
@test box12 == Interval(-5,5)

@test indx2 == 0
@test box21.hi > -Inf
@test box22 == Interval(-5,5)

@test indx3 == 0
@test box31.lo < Inf
@test box32 == Interval(-5,5)

@test indx6 == 1
@test -15.0002 - 1E-4 <= box61.lo <= -15.0002 + 1E-4
@test box62.lo == -Inf
@test box61.hi == Inf
@test -0.599979 - 1E-4 <= box62.hi <= -0.599979 + 1E-4

@test indx8 == 0
@test box81 == Interval(-5,5)
@test box82 == Interval(-5,5)


# Tests strict inclusion procedure for MCMCInterval vectors
Y = [MCInterval(0,5),MCInterval(0,5),MCInterval(0,5)]
X1 = [MCInterval(1,2),MCInterval(1,2),MCInterval(1,2)]
X2 = [MCInterval(1,2),MCInterval(-10,10),MCInterval(1,2)]
X3 = [MCInterval(1,2),MCInterval(1,2),MCInterval(0,5)]
flag1 = EAGO.Strict_XinY(X1,Y)
flag2 = EAGO.Strict_XinY(X2,Y)
flag3 = EAGO.Strict_XinY(X3,Y)
@test flag1 == true
@test flag2 == false
@test flag3 == false

# Checks extended division routine
A1 = MCInterval(0,0)
A2 = MCInterval(0,3)
A3 = MCInterval(-2,0)
A4 = MCInterval(-3,2)
ind1,B1,C1 = EAGO.extDivide(A1)
ind2,B2,C2 = EAGO.extDivide(A2)
ind3,B3,C3 = EAGO.extDivide(A3)
ind4,B4,C4 = EAGO.extDivide(A4)
@test ind1 == 0
@test ind2 == 1
@test ind3 == 2
@test ind4 == 3
@test B1 == MCInterval(-Inf,Inf)
@test 0.33333 - 1E-4 <= B2.lo <= 0.33333 + 1E-4
@test B2.hi == Inf
@test B3 == MCInterval(-Inf,-0.5)
@test B4.lo == -Inf
@test -0.33333 - 1E-4 <= B4.hi <= -0.33333 + 1E-4
@test C1 == MCInterval(-Inf,Inf)
@test C2 == MCInterval(Inf,Inf)
@test C3 == MCInterval(-Inf,-Inf)
@test C4 == MCInterval(0.5,Inf)

#param = EAGO.Param_Bisect_Opts()
#@test param.DAGflag == false

N =  MCInterval(-5,5)
X = MCInterval(-5,5)
Mii = MCInterval(-5,5)
S1 = MCInterval(-5,5)
S2 = MCInterval(-5,5)
B = MCInterval(-5,5)
rtol = 1E-4
indx1,box11,box12 = EAGO.extProcess(N,X,Mii,S1,S2,B,rtol)
Miib = MCInterval(0,5)
S1b = MCInterval(1,5)
S2b = MCInterval(1,5)
Bb = MCInterval(1,5)
indx2,box21,box22 = EAGO.extProcess(N,X,Miib,S1b,S2b,Bb,rtol)
Miic = MCInterval(-5,0)
S1c = MCInterval(1,5)
S2c = MCInterval(1,5)
Bc = MCInterval(1,5)
indx3,box31,box32 = EAGO.extProcess(N,X,Miic,S1c,S2c,Bc,rtol)
Miia = MCInterval(1,5)
S1a = MCInterval(1,5)
S2a = MCInterval(1,5)
Ba = MCInterval(1,5)
indx6,box61,box62 = EAGO.extProcess(N,X,Miia,S1a,S2a,Ba,rtol)
Miid = MCInterval(0,0)
S1d = MCInterval(1,5)
S2d = MCInterval(1,5)
Bd = MCInterval(1,5)
indx8,box81,box82 = EAGO.extProcess(N,X,Miid,S1d,S2d,Bd,rtol)

@test indx1 == 0
@test box11 == MCInterval(-Inf,Inf)
@test box12 == MCInterval(-5,5)

@test indx2 == 0
@test box21.hi > -Inf
@test box22 == MCInterval(-5,5)

#@test indx3 == 0
#@test box31.lo < Inf
#@test box32 == MCInterval(-5,5)

@test indx6 == 1
@test -15.0002 - 1E-4 <= box61.lo <= -15.0002 + 1E-4
@test box62.lo == -Inf
@test box61.hi == Inf
@test -0.599979 - 1E-4 <= box62.hi <= -0.599979 + 1E-4

@test indx8 == 0
@test box81 == MCInterval(-5,5)
@test box82 == MCInterval(-5,5)
end
