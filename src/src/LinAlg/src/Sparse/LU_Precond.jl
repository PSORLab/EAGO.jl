"""
     SparseInSto

Storage type for in-place sparse calculations
"""
type SparseInSto
     Xh
     Xj
     Yh
     Yj
     Zh
     Zj
     nx
end
SparseInSto() = SparseInSto([],[],[],[],[],[],[])

"""
     Sparse_Precondition!(H::Vector{Ta},J::SparseMatrixCSC{Ta,Ti},
                          P::SparseMatrixCSC{Tp,Ti},st::SparseInSto)

Preconditions the H & J to inv(P)H and inv(P)J using a sparse LU factorization
method with full pivoting. J and P must be of size nx-by-nx and H must be of
size nx. st is the inplace storage type.
"""
function Sparse_Precondition!(H::Vector{Ta},J::SparseMatrixCSC{Ta,Ti},
                              P::SparseMatrixCSC{Tp,Ti},st::SparseInSto) where {Ta,Tp,Ti}

     # generate LU-PDQ factorization
     lu = lufact(P)

     # solves Lz = PDH for z
     Sparse_Forward_Elimination!(st.Zh,lu[:L],(lu[:Rs].*H)[lu[:p]],st.nx)
     Sparse_Forward_Elimination!(st.Zj,lu[:L],(lu[:Rs].*J)[lu[:p],:],st.nx)

     # solves Uy = z for y
     Sparse_Back_Substitution!(st.Yh,lu[:U],st.Zh,st.nx)
     Sparse_Back_Substitution!(st.Yj,lu[:U],st.Zj,st.nx)

     # solves x = Qy
     st.Xh = st.Yh[lu[:q]]
     st.Xj = st.Yj[lu[:q],:]

     # stores the preconditioned matrices back in place
     H[:] = st.Xh[:,1]
     J[:] = st.Xj[:,1:(st.nx)]
end
