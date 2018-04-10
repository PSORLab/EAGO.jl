module Access_Test

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

S = BnBSolver()
B = BnBModel([Interval(1.0,2.0)])
#B1 = BnBModel()
tL = 3.0
tU = 5.0
X1 = [Interval(1.0,2.0)]
X2 = [Interval(3.0,5.0)]
X = [Interval(3.0,5.0)]
pos = 4
EAGO.BM_breadth!(S,B,tL,tU,X1,X2,pos)
EAGO.BM_Single!(S,B,tL,tU,X,pos)
B.soln = [7.0]
B.UBDg = 14.7

#@test B1.Init_Box == [Interval(0,1)]
@test getsolution(B) == [7.0]
@test getobjval(B) == 14.7
@test getobjbound(B) == 14.7
@test getfeasibility(B) == false
@test LBDtime(B) == 0.0
@test UBDtime(B) == 0.0
end
