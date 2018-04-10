module Checks_Test

using Compat
using Compat.Test
using IntervalArithmetic
using EAGOBranchBound

S = BnBSolver()
B = BnBModel([Interval(1.0,2.0)])

B.LBD = Float64[]
B.first_num = 100
@test EAGO.Term_Check(S,B,10) == false

B.first_num = 0
@test EAGO.Term_Check(S,B,10) == false

S.max_nodes = -1
B.LBD = Float64[1.0]
@test EAGO.Term_Check(S,B,10) == false

S.max_nodes = 100
S.iter_lim = true
S.max_iter = 5
@test EAGO.Term_Check(S,B,10) == false
end
