"""
     Sparse_Forward_Elimination!(x::Vector{Tq},L::SparseMatrixCSC{Tv,Ti},
                                 b::Vector{Tq},nx::Int64) where {Tv,Tq,Ti}

Solves for `Lx=b` via forward elimination. A must be a lower triangular sparse
matrix of CSC format.
"""
function Sparse_Forward_Elimination!(x::Vector{Tq},L::SparseMatrixCSC{Tv,Ti},
                                     b::Vector{Tq},nx::Int64) where {Tv,Tq,Ti}
     # converts to CSR (expensive if dense ... much cheaper calc if sparse)
     A::SparseMatrixCSC{Tv,Ti} = L'
     # standard row-oriented forward elimination
     x[1,:] = b[1,:]/A.nzval[1]
     for i=2:(nx)
          for k=(A.colptr[i]):(A.colptr[i+1]-2)
               b[i,:] = b[i,:] - A.nzval[k]*x[A.rowval[k],:]
          end
          x[i,:] = b[i,:]/A.nzval[A.colptr[i+1]-1]
     end
end

"""
     Sparse_Forward_Elimination!(x::SparseMatrixCSC{Tq,Ti},U::SparseMatrixCSC{Tv,Ti},
                                  b::SparseMatrixCSC{Tq,Ti},nx::Int64)

Solves for `Lx=b` via forward elimination. A must be a lower triangular sparse
matrix of CSC format.
"""
# TO DO: Improve row access for b
function Sparse_Forward_Elimination!(x::SparseMatrixCSC{Tq,Ti},L::SparseMatrixCSC{Tv,Ti},
                                     b::SparseMatrixCSC{Tq,Ti},nx::Int64) where {Tv,Tq,Ti}
     # converts to CSR (expensive if dense ... much cheaper calc if sparse)
     A::SparseMatrixCSC{Tv,Ti} = L'
     # standard row-oriented forward elimination
     x[1,:] = b[1,:]/A.nzval[1]
     for i=2:(nx)
          for k=(A.colptr[i]):(A.colptr[i+1]-2)
               b[i,:] = b[i,:] - A.nzval[k]*x[A.rowval[k],:]
          end
          x[i,:] = b[i,:]/A.nzval[A.colptr[i+1]-1]
     end
end

"""
     Sparse_Back_Substitution!(x::Vector{Tq},U::SparseMatrixCSC{Tv,Ti},
                              b::Vector{Tq},nx::Int64)

Solves for `Ux=b` via backsubstitution. A must be a upper triangular sparse
matrix of CSC format.
"""
function Sparse_Back_Substitution!(x::Vector{Tq},U::SparseMatrixCSC{Tv,Ti},
                                   b::Vector{Tq},nx::Int64) where {Tv,Tq,Ti}
     # converts to CSR (expensive if dense ... much cheaper calc if sparse)
     A::SparseMatrixCSC{Tv,Ti} = U'
     # standard row-oriented back substituion
     x[end,:] = b[end,:]/A.nzval[end]
     for i=(nx-1):-1:1
          for k=(A.colptr[i+1]-1):-1:(A.colptr[i]+1)
               b[i,:] = b[i,:] - A.nzval[k]*x[A.rowval[k],:]
          end
          x[i,:] = b[i,:]/A.nzval[A.colptr[i]]
     end
end

"""
     Sparse_Back_Substitution!(x::SparseMatrixCSC{Tq,Ti},U::SparseMatrixCSC{Tv,Ti},
                               b::SparseMatrixCSC{Tq,Ti},nx::Int64)

Solves for `Ux=b` via backsubstitution. A must be a upper triangular sparse
matrix of CSC format.
"""
# TO DO: Improve row access for b
function Sparse_Back_Substitution!(x::SparseMatrixCSC{Tq,Ti},U::SparseMatrixCSC{Tv,Ti},
                                   b::SparseMatrixCSC{Tq,Ti},nx::Int64) where {Tv,Tq,Ti}
     # converts to CSR (expensive if dense ... much cheaper calc if sparse)
     A::SparseMatrixCSC{Tv,Ti} = U'
     # standard row-oriented back substituion
     x[end,:] = b[end,:]/A.nzval[end]
     for i=(nx-1):-1:1
          for k=(A.colptr[i+1]-1):-1:(A.colptr[i]+1)
               b[i,:] = b[i,:] - A.nzval[k]*x[A.rowval[k],:]
          end
          x[i,:] = b[i,:]/A.nzval[A.colptr[i]]
     end
end
