module Node_Test

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

B = BnBModel([Interval(1.0,2.0),Interval(3.0,4.0)])
A1,A2,A3,A4,A5 = EAGO.NS_depth_breadth(B)
@test A1 == [Interval(1.0,2.0),Interval(3.0,4.0)]
@test A2 == -Inf
@test A3 == Inf
@test A4 == 1
@test A5 == 1

end
