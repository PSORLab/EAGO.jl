module Checks_Test

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

@testset "Test Termination Checks" begin
S = BnBSolver()
B = BnBModel([Interval(1.0,2.0)])

B.LBD = Float64[]
B.first_num = 100
@test EAGO.Term_Check(S,B,Int64(10)) == false

B.first_num = 0
@test EAGO.Term_Check(S,B,Int64(10)) == false

S.max_nodes = -1
B.LBD = Float64[1.0]
@test EAGO.Term_Check(S,B,Int64(10)) == false

S.max_nodes = 100
S.iter_lim = true
S.max_iter = 5
@test EAGO.Term_Check(S,B,Int64(10)) == false
end

@testset "Test Default Pre/Post" begin
XM = [MCInterval(1,5)]
S = BnBSolver()
model = BnBModel(XM)
set_to_default!(S)
fo1,Xo1 = default_pre(true,XM,Float64(1.0),Int64(2),Int64(3),[],Float64(1.0),Float64(1.0),S,model)
fo2,Xo2 = default_post(true,XM,Int64(2),Int64(3),[],1.0,2.0,2.0,3.0)
@test Xo1[1] == XM[1]
@test Xo2[1] == XM[1]
@test fo1
@test fo2
end
end
