module ParamIntvPrecond_Tests

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

@testset "Preconditioners" begin
opt1 = PIntvParams(:None,:Newton,1E-30,1E-6,Int64(2),Int64(2),100)
h(x,p) = [p[1]]
hj(x,p) = [p[1]]
X = [MCInterval(1.0,2.0),MCInterval(1.0,2.0)]
P = [MCInterval(1.0,2.0),MCInterval(1.0,3.0)]

X1,P2 = EAGO.Precondition(h,hj,X,P,opt1)
@test X == X
@test P == P

opt2 = PIntvParams(:XYZYZYYA,:Newton,1E-30,1E-6,Int64(2),Int64(2),100)
@test_throws ErrorException EAGO.Precondition(h,hj,X,P,opt2)
end

@testset "Preconditioner Utilities" begin
A = speye(spzeros(Float64,9,9))
A[2,3] = 2.1
A[8,7] = 3.1
kl,ku = EAGO.SparseBandwidth(A)
@test kl == 1
@test ku == 1
end

end
