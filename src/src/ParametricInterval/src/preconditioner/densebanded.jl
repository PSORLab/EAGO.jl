#=
function VecSwap!(n::Int64,incx::Int64,incy::Int64,dx::Array,dy::Array) where {T<:AbstractFloat}
     temp = copy(dx)
     dx = dy
     dy = temp
end

function plus_xyT!(m,n,X,Y,incx,incy,A)
     if (n == 0) or (m == 0)
          error("array length wrong....")
     end
     jy = ((incy > 0) ? 1 : (n-1)*incy)
     if (incx == 1)
          for i=1:n
               if (Y[jy] != 0)
                    temp = -Y[jy]
                    for i=1:m
                         A[i,j] = A[i,j] + X[i]*temp
                    end
               end
          end
     else
          kx = (incx>0) ?  1 : (1-(m-1)*incx)
          for j=1:n
               if (Y[jy] != 0)
               temp = -Y[jy]
               ix = kx
               for i=1:m
                    A[i,j] = A[i,j] + X[ix]*temp
                    ix = ix + incx
               end
          end
          jy = jy +  incy
     end
end
=#
"""
N - matrix order
K - number of superdiagonals
A - leading (k+1)-by-n part is upper triangular band, supplied column by column
    the following program converts matrix from regular to band storage
        for j = 1:N
             M = K + 1 - j
             for i = MAX( 1, j - K ):j
                  A[M+i,j] = matrix[i,j]
             end
        end
lda - leading dimension of A
X -  Dimension at least ( 1 + ( n - 1 )*abs( INCX ) ). Before entry, the incremented
     array X must contain the n element right-hand side vector b.
     On exit, X is overwritten with the solution vector x.
IncX -  increment for the elements of x.
"""
function Usolver_Band!(n,k,A,lda,X,INCX)

     if (incx <= 0)
          kx = 1 - (n-1)*incx
     elseif  (incx != 1)
          kx = 1
     end

     kplus1 = k + 1
     if (incx == 1)
          for j = n:-1:1
               if (X[j] != 0)
                    l = kplus1 - j
                    X[j] = X[j]/A[kplus1,j]
                    temp = x[j]
                    for i = (j-1):-1:max(1,j-k)
                         x[i] = x[i] - temp*A[l+i,j]
                    end
               end
          end
     end
     kx = kx + (n-1)*incx
     jx = kx
     for j = n:-1:1
         kx = kx - incx
         if (X[jx] != 0)
             ix = kx
             l = kplus1 - j
             X[jx] = X[jx]/A[kplus1,j]
             temp = x[jx]
             for i = (j-1):-1:max(1,j-k)
                 X[ix] = X[ix] - temp*A[l+i,j]
                 ix = ix - incx
            end
          end
         jx = jx - incx
    end
end

function LUsolveB_Band(LU::BandedMatrices.BandedLU{T},B::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti<:Integer,T<:AbstractFloat}
     N = B.n
     nrhs = B.m
     kl = LU.l
     ku = LU.u
     data = LU.data
     ipiv = LU.ipiv
     #=
     if (kl>0) # matrix not upper triangular (condition via lower triangular)
          for j=1:(N-1)                                                            # for every row
               lm = min(kl,N-j)                                                    #
               l = ipiv[j]                                                         # get pivot index for row j
               if (l!=j)                                                           #
                    #VecSwap!(nrhs, N, N, b( l, 1 ), b( j, 1 ))
                    #plus_xyT!(lm,nrhs,ab(kd+1,j),b(j,1),1,N,b(j+1,1))
                    VecSwap!(nrhs, N, N, view(B,l,:), view(B,j,:))
                    plus_xyT!(lm,nrhs,data(kd+1,j),view(B,j,:),1,N,view(B,j+1,:))
               end
          end
     end
     =#
     for i=1:nrhs
          #Usolver_Band!(N, kl+ku, ab, ldab, b( 1, i ), 1 )
          Usolver_Band!(N, kl+ku, ab, ldab, view(B,:,i), 1 ) # dimensions check out
     end
end

"""
    DenseBand_Precondition!(H,J,Y,opt)
"""
function DenseBand_Precondition!(H::Union{Array{Th,2},Vector{Th}},
                                 J::BandedMatrix{Tj},
                                 Y::BandedMatrix{Ty},
                                 opt::PIntvParams) where {Th,Tj,Ty}
    if (opt.nx == 1)
        H[:] = H/Y[i,i]
        J[:] = J/Y[i,i]
    else
        F = lufact(Y)
        A_ldivB!(H,F,H)
        A_ldivB!(J,F,J)
    end
end
