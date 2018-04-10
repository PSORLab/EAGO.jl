module SparseCntr

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

nx = 3
A = spzeros(Float64,nx,nx)
c = [3.0, 0.5, 1.0]
P = spzeros(Float64,nx,nx)
A[1,1] = P[1,1] = 0.5
A[1,2] = P[1,2] = 0.2
A[1,3] = P[1,3] = 0.1
A[2,1] = P[2,1] = 0.3
A[2,2] = P[2,2] = 0.5
A[2,3] = P[2,3] = 0.5
A[3,1] = P[3,1] = 0.4
A[3,2] = P[3,2] = 0.5
A[3,3] = P[3,3] = 0.3

lufct = lufact(P)
At = copy(A)
ct = copy(c)
Pt = copy(P)
preA = inv(lufct)*A
prec = inv(lufct)*c
SSto = EAGO.SparseInSto()
SSto.Xh = zeros(Float64,nx)
SSto.Yh = zeros(Float64,nx)
SSto.Zh = zeros(Float64,nx)
SSto.Xj = spzeros(Float64,nx,nx)
SSto.Yj = spzeros(Float64,nx,nx)
SSto.Zj = spzeros(Float64,nx,nx)
SSto.nx = nx
EAGO.Sparse_Precondition!(c,A,P,SSto)

@test (7.7272-10E-4 <= c[1] <= 7.7272+10E-4)
@test (-5.0-10E-4 <= c[2] <= -5.0+10E-4)
@test (1.3636-10E-4 <= c[3] <= 1.3636+10E-4)

end
