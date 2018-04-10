module Bisect_Tests

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

S = BnBSolver()
B = BnBModel([Interval(1.0,2.0),Interval(3.0,4.0),Interval(3.0,5.0),Interval(3.0,6.0)])
N = [Interval(1.0,2.0),Interval(3.0,4.0),Interval(4.0,5.0),Interval(4.0,6.0)]
nx = 2

X1,X2 = EAGO.Bisect_Rel(S,B,N)
Y1,Y2 = EAGO.Bisect_Abs_Imp(S,B,N,nx)
Z1,Z2 = EAGO.Bisect_Rel_Imp(S,B,N,nx)
@test X1[1] == Interval(1.0,1.5)
@test X2[1] == Interval(1.5,2.0)
@test X2[2:4] == [Interval(3.0,4.0),Interval(4.0,5.0),Interval(4.0,6.0)]
@test Y1[1:3] == [Interval(1.0,2.0),Interval(3.0,4.0),Interval(4.0,5.0)]
@test Y1[4] == Interval(4.0,5.0)
@test Y2[4] == Interval(5.0,6.0)
@test Z1[1:3] == [Interval(1.0,2.0),Interval(3.0,4.0),Interval(4.0,5.0)]
@test Z1[4] == Interval(4.0,5.0)
@test Z2[4] == Interval(5.0,6.0)

end
