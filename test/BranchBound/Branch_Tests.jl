module Branch_Test

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

S = BnBSolver()
B = BnBModel([Interval(1.0,2.0)])
tL = 3.0
tU = 5.0
X1 = [Interval(1.0,2.0)]
X2 = [Interval(3.0,5.0)]
X = [Interval(3.0,5.0)]
pos = Int64(4)
EAGO.BM_breadth!(S,B,tL,tU,X1,X2,pos)
@test (B.pos[1:3] == [5; 5; 1])
@test (B.LBD[1:3] == [3.0; 3.0; -Inf])
@test (B.UBD[1:3] == [5.0; 5.0; Inf])
@test (B.id[1:3] == [2; 3; 1])

EAGO.BM_Single!(S,B,tL,tU,X,pos)
@test (B.pos[4] == 4)
@test (B.LBD[4] == 3.0)
@test (B.UBD[4] == 5.0)
@test (B.id[4] == 3)

end
