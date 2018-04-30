module TestDAGcntr

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

@testset "Test Tape Generation" begin
single_expr = :(x[1]+x[2]^2-exp(x[1]))
tape_out = Generate_Tape(single_expr,2,-2,4,Interval{Float64})
dag_out = EAGO.getDAG(single_expr,[1,2],Interval{Float64})
@test dag_out[1] == Int64[1; 2; 3; 4; 5; 6; 7]
@test dag_out[3][4] == Int64[4; 5]

mult_expr = [:(x[1]+x[2]^2)
             :(x[1]+x[5])
            ]
X = [Interval(-10.0,20.0) for i=1:5]
ftapelist_out = Generate_Fixed_TapeList(mult_expr,4,[-8.0,0.5],[2.0,4.5],[[5.0]],Interval{Float64})
DAGContractor!(X,ftapelist_out,6)
@test -10-1E-4 <= X[1].lo <= -10+1E-4
@test 0-1E-4 <= X[2].lo <= 0+1E-4
@test -10-1E-4 <= X[3].lo <= -10+1E-4
@test -10-1E-4 <= X[4].lo <= -10+1E-4
@test 5-1E-4 <= X[5].lo <= 5+1E-4
@test 20-1E-4 <= X[1].hi <= 20+1E-4
@test 20-1E-4 <= X[2].hi <= 20+1E-4
@test 20-1E-4 <= X[3].hi <= 20+1E-4
@test 20-1E-4 <= X[4].hi <= 20+1E-4
@test 5-1E-4 <= X[5].hi <= 5+1E-4

@test Tape() == Tape(MCInterval{Float64})
@test TapeList().sto == TapeList([]).sto
end

@testset "Test Node Finder Construction" begin
a = NodeFinder(1)
@test exp2(a) == EAGO.NodeFinder(7)
@test exp10(a) == EAGO.NodeFinder(8)
@test log(a) == EAGO.NodeFinder(9)
@test log2(a) == EAGO.NodeFinder(10)
@test log10(a) == EAGO.NodeFinder(11)
@test acosh(a) == EAGO.NodeFinder(12)
@test cosh(a) == EAGO.NodeFinder(13)
@test sqrt(a) == EAGO.NodeFinder(14)
@test asin(a) == EAGO.NodeFinder(15)
@test sinh(a) == EAGO.NodeFinder(16)
@test atanh(a) == EAGO.NodeFinder(17)
@test tan(a) == EAGO.NodeFinder(18)
@test atan(a) == EAGO.NodeFinder(19)
@test acos(a) == EAGO.NodeFinder(20)
@test tanh(a) == EAGO.NodeFinder(21)
@test asinh(a) == EAGO.NodeFinder(22)
@test sin(a) == EAGO.NodeFinder(23)
@test cos(a) == EAGO.NodeFinder(24)
@test abs(a) == EAGO.NodeFinder(25)
@test step(a) == EAGO.NodeFinder(26)
@test sign(a) == EAGO.NodeFinder(27)
@test inv(a) == EAGO.NodeFinder(28)
@test -(a) == EAGO.NodeFinder(29)
@test one(a) == EAGO.NodeFinder(30)
@test zero(a) == EAGO.NodeFinder(31)
@test /(a,a) == EAGO.NodeFinder(32)
@test -(a,a) == EAGO.NodeFinder(33)
@test +(a,a) == EAGO.NodeFinder(34)
@test *(a,a) == EAGO.NodeFinder(35)
@test min(a,a) == EAGO.NodeFinder(36)
@test max(a,a) == EAGO.NodeFinder(37)
@test ^(a,a) == EAGO.NodeFinder(38)
end

@testset "Test Reverse MCInterval/Interval" begin
a = Interval(1.0,2.0)
b = Interval(1.0,20.0)
b1 = Interval(-1.0,3.0)
c = Interval(1.0,5.0)
c1 = Interval(-1.0,5.0)
e = emptyinterval()

ma = MCInterval(1.0,2.0)
mb = MCInterval(1.0,20.0)
mb1 = MCInterval(-1.0,3.0)
mc = MCInterval(1.0,5.0)
mc1 = MCInterval(-1.0,5.0)
me = EAGO.emptyMCinterval(Float64)

o1,o2,o3 = EAGO.mul_revDR(a,b1,c)
@test o1 == Interval(1.0,2.0)
@test isapprox(o2.lo,0.199999,atol=1E-4)
@test isapprox(o2.hi,2.0,atol=1E-4)
@test o3 == Interval(1.0,5.0)

o1,o2,o3 = EAGO.mul_revDR(a,b,c1)
@test o1 == Interval(1.0,2.0)
@test o2 == Interval(1.0,20.0)
@test isapprox(o3.lo,0.0499999,atol=1E-4)
@test isapprox(o3.hi,2.0,atol=1E-4)

o1,o2,o3 = EAGO.mul_revDR(ma,mb1,mc)
@test o1 == MCInterval(1.0,2.0)
@test isapprox(o2.lo,0.199999,atol=1E-4)
@test isapprox(o2.hi,2.0,atol=1E-4)
@test o3 == MCInterval(1.0,5.0)

o1,o2,o3 = EAGO.mul_revDR(ma,mb,mc1)
@test o1 == MCInterval(1.0,2.0)
@test o2 == MCInterval(1.0,20.0)
@test isapprox(o3.lo,0.0499999,atol=1E-4)
@test isapprox(o3.hi,2.0,atol=1E-4)


o1,o2 = EAGO.exp_revDR(a,b)
@test o1 == Interval(1.0,2.0)
@test o2 == ∅

o1,o2 = EAGO.exp_revDR(mb,ma)
@test o1 == MCInterval(1.0,20.0)
@test o2 == MCInterval(1.0,2.0)

o1,o2 = EAGO.acos_rev(a,b)
@test o1 == Interval(1.0,2.0)
@test o2 == ∅

o1,o2 = EAGO.atan_rev(a,b)
@test isapprox(o1.hi,1.5708,atol=1E-4)
@test isapprox(o2.lo,1.5574,atol=1E-4)

o1,o2 = EAGO.atan_rev(ma,mb)
@test isapprox(o1.hi,1.5708,atol=1E-4)
@test isapprox(o2.lo,1.5574,atol=1E-4)

o1,o2 = EAGO.sinh_rev(a,b)
@test isapprox(o1.hi,2,atol=1E-4)
@test isapprox(o2.hi,1.44364,atol=1E-4)

o1,o2 = EAGO.sinh_rev(ma,mb)
@test isapprox(o1.hi,2,atol=1E-4)
@test isapprox(o2.hi,1.44364,atol=1E-4)

o1,o2 = EAGO.cosh_rev(a,b)
@test isapprox(o1.hi,2,atol=1E-4)
@test isapprox(o2.hi,1.3169578969248166,atol=1E-4)

o1,o2 = EAGO.cosh_rev(ma,mb)
@test isapprox(o1.hi,2,atol=1E-4)
@test isapprox(o2.hi,1.3169578969248166,atol=1E-4)

o1,o2 = EAGO.tanh_rev(a,b)
@test isapprox(o1.hi,1,atol=1E-4)
@test isapprox(o2.hi,-Inf,atol=1E-4)

o1,o2 = EAGO.tanh_rev(ma,mb)
@test isapprox(o1.hi,1,atol=1E-4)
@test isapprox(o2.hi,-Inf,atol=1E-4)

o1,o2 = EAGO.asinh_rev(a,b)
@test isapprox(o2.lo,1.1752,atol=1E-4)
@test isapprox(o2.hi,3.62687,atol=1E-4)

o1,o2 = EAGO.asinh_rev(ma,mb)
@test isapprox(o2.lo,1.175201193643801,atol=1E-4)
@test isapprox(o2.hi,3.626860407847019,atol=1E-4)

o1,o2 = EAGO.acosh_rev(a,b)
@test isapprox(o2.lo,1.54308,atol=1E-4)
@test isapprox(o2.hi,3.7622,atol=1E-4)

o1,o2 = EAGO.acosh_rev(ma,mb)
@test isapprox(o2.lo,1.54308,atol=1E-4)
@test isapprox(o2.hi,3.7622,atol=1E-4)

o1,o2 = EAGO.atanh_rev(a,b)
@test isapprox(o1.hi,2,atol=1E-4)
@test isapprox(o2.hi,-Inf,atol=1E-4)

o1,o2 = EAGO.atanh_rev(ma,mb)
@test isapprox(o1.hi,2,atol=1E-4)
@test isapprox(o2.hi,-Inf,atol=1E-4)

o1,o2 = EAGO.step_rev(a,b)
@test isapprox(o1.hi,1,atol=1E-4)
@test isapprox(o2.hi,20,atol=1E-4)

o1,o2 = EAGO.step_rev(ma,mb)
@test isapprox(o1.hi,1,atol=1E-4)
@test isapprox(o2.hi,20,atol=1E-4)

o1,o2 = EAGO.sign_rev(a,b)
@test isapprox(o1.hi,1,atol=1E-4)
@test isapprox(o2.hi,20,atol=1E-4)

o1,o2 = EAGO.sign_rev(ma,mb)
@test isapprox(o1.hi,1,atol=1E-4)
@test isapprox(o2.hi,20,atol=1E-4)

o1,o2 = EAGO.exp2_rev(a,b)
@test isapprox(o1.hi,2,atol=1E-4)
@test isapprox(o2.hi,1,atol=1E-4)

o1,o2 = EAGO.exp2_rev(ma,mb)
@test isapprox(o1.hi,2,atol=1E-4)
@test isapprox(o2.hi,1,atol=1E-4)

o1,o2 = EAGO.exp10_rev(a,b)
@test isapprox(o1.hi,2,atol=1E-4)
@test isapprox(o2.hi,-Inf,atol=1E-4)

o1,o2 = EAGO.exp10_rev(ma,mb)
@test isapprox(o1.hi,2,atol=1E-4)
@test isapprox(o2.hi,-Inf,atol=1E-4)

o1,o2 = EAGO.log2_rev(a,b)
@test isapprox(o2.lo,2,atol=1E-4)
@test isapprox(o2.hi,4,atol=1E-4)

o1,o2 = EAGO.log2_rev(ma,mb)
@test isapprox(o2.lo,2,atol=1E-4)
@test isapprox(o2.hi,4,atol=1E-4)

o1,o2 = EAGO.log10_rev(a,b)
@test isapprox(o2.lo,10,atol=1E-4)
@test isapprox(o2.hi,20,atol=1E-4)

o1,o2 = EAGO.log10_rev(ma,mb)
@test isapprox(o2.lo,10,atol=1E-4)
@test isapprox(o2.hi,20,atol=1E-4)

o1,o2 = EAGO.one_rev(a,b)
@test isapprox(o2.lo,1,atol=1E-4)
@test isapprox(o2.hi,20,atol=1E-4)

o1,o2 = EAGO.one_rev(ma,mb)
@test isapprox(o1.hi,1,atol=1E-4)
@test isapprox(o2.hi,20,atol=1E-4)

o1,o2 = EAGO.zero_rev(a,b)
@test isapprox(o1.hi,-Inf,atol=1E-4)
@test isapprox(o2.hi,20,atol=1E-4)

o1,o2 = EAGO.zero_rev(ma,mb)
@test isapprox(o1.hi,-Inf,atol=1E-4)
@test isapprox(o2.hi,20,atol=1E-4)

@test EAGO.mul_revDR(a,b,c) == EAGO.mul_revDR(promote(a,b,c)...)
@test EAGO.exp_revDR(a,b) == EAGO.exp_revDR(promote(a,b)...)
@test EAGO.sinh_rev(a,b) == EAGO.sinh_rev(promote(a,b)...)
@test EAGO.cosh_rev(a,b) == EAGO.cosh_rev(promote(a,b)...)
@test EAGO.tanh_rev(a,b) == EAGO.tanh_rev(promote(a,b)...)
@test EAGO.asinh_rev(a,b) == EAGO.asinh_rev(promote(a,b)...)
@test EAGO.acosh_rev(a,b) == EAGO.acosh_rev(promote(a,b)...)
@test EAGO.atanh_rev(a,b) == EAGO.atanh_rev(promote(a,b)...)
@test EAGO.step_rev(a,b) == EAGO.step_rev(promote(a,b)...)
@test EAGO.sign_rev(a,b) == EAGO.sign_rev(promote(a,b)...)
@test EAGO.exp2_rev(a,b) == EAGO.exp2_rev(promote(a,b)...)
@test EAGO.exp10_rev(a,b) == EAGO.exp10_rev(promote(a,b)...)
@test EAGO.log2_rev(a,b) == EAGO.log2_rev(promote(a,b)...)
@test EAGO.log10_rev(a,b) == EAGO.log10_rev(promote(a,b)...)
@test EAGO.one_rev(a,b) == EAGO.one_rev(promote(a,b)...)
@test EAGO.zero_rev(a,b) == EAGO.zero_rev(promote(a,b)...)


#=
o1,o2,o3 = EAGO.min_rev(a,b,c)
println("mult ")
println("o1: $(o1)")
println("o2: $(o2)")
println("o3: $(o3)")

o1,o2,o3 = EAGO.max_rev(a,b,c)
println("mult ")
println("o1: $(o1)")
println("o2: $(o2)")
println("o3: $(o3)")
=#

end

end
