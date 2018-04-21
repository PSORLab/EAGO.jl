#=

function LowMatMult!(a::Array{T,2},b::Array{T,2},m::Int64,n::Int64) where {T<:AbstractFloat}
    zeroT::T = zero(T)
    @inbounds  for i=1:m
        @inbounds for j = 1:n
            if (b[j,i] != zeroT)
                temp::T = b[j,i]
                @inbounds for k=(j+1):(n)
                    b[k,i] = b[k,i]-a[k,j]*temp
                end
            end
        end
    end
end

function UppMatMult!(a::Array{T,2},b::Array{T,2},m::Int64,n::Int64) where {T<:AbstractFloat}
    zeroT::T = zero(T)
    @inbounds for i = 1:m
        @inbounds for j = n:-1:1
            if (b[j,i] != zeroT)
                temp::T = b[j,i]
                b[j,i] = temp/a[j,j]
                @inbounds for k=1:(j-1)
                    b[k,i] = b[k,i]-a[k,j]*temp
                end
            end
        end
    end
end

function Dense_STD_Precondition!(H::Vector{T},J::Array{T,2},Y::Array{T,2},nx::Int64) where {T<:AbstractFloat}
    lu = lufact(Y)
    HJ::Array{T,2} = [H J][lu[:p],:]
    LowMatMult!(lu[:L],HJ,nx+1,nx)
    UppMatMult!(lu[:U],HJ,nx+1,nx)
    H[:] = HJ[:,1]
    J[:] = HJ[:,2:(nx+1)]
end


function LowBandMatMult!(a::Array{T,2},x::Array{T,1},n::Int64,k::Int64) where {T<:AbstractFloat}

    if (incx <= 0)
        kx = 1 - (n-1)*incx
    else (incx != 1)
        kx = 1
    end

    if incx == 1
        for j = 1:n
            if (x[j] != zero(T))
                l = 1 - j
                x[j] = x[j]/a[1,j]
                temp = x[j]
                for i = (j + 1):(min(n,j+k))
                    x[i] = x[i] - temp*a[l+i,j]
                end
            end
        end
    else
        jx = kx
        for j = 1:n
            kx = kx + incx
            if (x[jx] != 0)
                ix = kx
                l = 1 - j
                x[jx] = x[jx]/a[1,j]
                temp = x[jx]
                for i = (j + 1):(min(n,j+k))
                    x[ix] = x[ix] - temp*a[l+i,j]
                    ix = ix + incx
                end
            end
            jx = jx + incx
        end
    end
end

function UppBandMatMult!(a::Array{T,2},x::Array{T,1},n::Int64,k::Int64,incx::Int64) where {T<:AbstractFloat}

    if (incx <= 0)
        kx = 1 - (n-1)*incx
    else (incx != 1)
        kx = 1
    end

    kplus1 = k + 1
    if incx == 1
        for j = n:-1:1
             if (x[j] != zero(T))
                 l = kplus1 - j
                 x[j] = x[j]/a[kplus1,j]
                 temp = x[j]
                 for i = (j-1):-1:max(1,j-k)
                     x[i] = x[i] - temp*a[l+i,j]
                 end
             end
        end
    else
        kx = kx + (n-1)*incx
        jx = kx
        for j = n:-1:1
            kx = kx - incx
            if (x[jx] != zero(T))
                ix = kx
                l = kplus1 - j
                x[jx] = x[jx]/a[kplus1,j]
                temp = x[jx]
                for i = (j-1):-1:(max(1,j-k))
                    x[ix] = x[ix] - temp*a[l+i,j]
                    ix = ix - incx
                end
            jx = jx - incx
            end
        end
    end
end

function ConvLowBand(L::Array{T,2},nx::Int64,K::Int64) where {T<:AbstractFloat}
    A::Array{T,2} = zeros(L)
    oneT::Int64 = 1
    @inbounds for j = 1:nx
        m::Int64 = oneT - j
        @inbounds for i = j:min(nx,j+K)
            A[m+i,j] = L[i,j]
        end
    end
    return A
end

function ConvUppBand(U::Array{T,2},nx::Int64,K::Int64) where {T<:AbstractFloat}
    A::Array{T,2} = zeros(U)
    oneT::Int64 = 1
    @inbounds for j = 1:nx
        m::Int64 = K + oneT - j
        @inbounds for i = max(oneT,j-K):(j)
            A[m+i,j] = U[I,J]
        end
    end
    return A
end

function VecInter!(x::Vector{T},y::Vector{T},incx::Int64,incy::Int64,n::Int64) where {T<:AbstractFloat}
    if (n<0)
        return
    else (incx == 1) && (incy == 1)
          m::Int64 = mod(n,3)
          if (m != 0)
              for i = 1:m
                  dtemp = dx[i]
                  dx[i] = dy[i]
                  dy[i] = dtemp
              end
             if (n < 3)
                 return
             end
         end
         mp1::Int64 = m + 1
         for i = mp1:3:n
             dtemp = dx[i]
             dx[i] = dy[i]
             dy[i] = dtemp
             dtemp = dx[i+1]
             dx[i+1] = dy[i+1]
             dy[i+1] = dtemp
             dtemp = dx[i+2]
             dx[i+2] = dy[i+2]
             dy[i+2] = dtemp
         end
     else
         ix::Int64 = 1
         iy::Int64 = 1
         (incx < 0) && ix = (-n+1)*incx + 1
         (incy < 0) && iy = (-n+1)*incy + 1
         if i = 1:n
             dtemp = dx[ix]
             dx[ix] = dy[iy]
             dy[iy] = dtemp
             ix = ix + incx
             iy = iy + incy
         end
     end
end

function xyT_pA!(m,n,x,incx,y,incy,A,LDA)
    if (incy > zero(T))
     jy = 1
    else
     jy = 1 - (n-1)*incy
    end
    if (incx == 1)
        for j = 1:n
            if (y[jy] != zero(T))
                temp = y[jy]
            end
            for i = 1:m
                a[i,j] = a[i,j] + x[i]*temp
            end
         jy = jy + incy
     else
     if (incx > 0)
         kx = 1
     else
         kx = 1 - (m-1)*incx
     end
     for j = 1:n
         if (y[jy] != zero(T))
             temp = y[jy]
             ix = kx
             for i = 1:m
                 a[i,j] = a[i,j] + x[ix]*temp
                 ix = ix + incx
             end
         end
         jy = jy + incy
     end
 end
end

function Dense_Band_Precondition!(H::Vector{T},J::Array{T,2},Y::Array{T,2},nx::Int64,kL::Int64,kU::Int64) where {T<:AbstractFloat}
    # factor and convert to band ordering
    lu = lufact(Y)
    L = ConvLowBand(lu[:L],nx,kL)
    U = ConvUppBand(lu[:U],nx,kU)
    HJ::Array{T,2} = [H J][lu[:p],:]

    # banded matrix preconditioning by column
    nrhs::Int64 = nx + 1
    kd::Int64 = ku + kl + 1
    lnoti::Bool = (kL>0)
    ldb =

    if (lnoti)
        for j = 1:(nx-1)
            lm = min(kL, nx-j)
            l = lu[:p][j]
            if (l != j)
                VecInter!(b[l,1],b[j,1],nx,nx,nx+1)
                xyT_pA!(lm,nx+1,-one(T),incx,y,incy,A,LDA)
                #CALL dger(lm, nx+1, -one, ab( kd+1, j ), 1, b( j, 1 ),nx, b( j+1, 1 ), nx )
            end
        end
    end

    for i = 1:(nx+1)
        UppBandMatMult!(U,HJ[],nx,kL+kU,1)
        #CALL dtbsv(nx, kL+kL,ab, ldab, b( 1, i ), 1 )
    end

    LowBandMatMult!(L,HJ,nx+1,nx,kL)
    UppBandMatMult!(U,HJ,nx+1,nx,kU)
    H[:] = HJ[:,1]
    J[:] = HJ[:,2:(nx+1)]
end

function Dense_Precondition!(H::Vector{T},J::Array{T,2},Y::Array{T,2},nx::Int64,kL::Int64,kU::Int64)
    if kL<0
        Dense_STD_Precondition!(H,J,Y,nx)
    else
        Dense_Band_Precondition!(H,J,Y,nx,kL,kU)
    end
end

nx = 1000
B = rand(Float64,nx,nx)
C = rand(Float64,nx)
A = zeros(Float64,nx,nx)
for i = 1:nx
    for j = 1:nx
        if abs(i-j) <= 1
            A[i,j] = rand(Float64)
        end
    end
end
nx1 = copy(nx)
B1 = copy(B)
B2 = copy(B)
B3 = copy(B)
C1 = copy(C)
C2 = copy(C)
A1 = copy(A)
A2 = copy(A)
A3 = copy(A)

#@code_warntype Dense_Precondition!(C1,B1,A1,nx1)
@time Dense_Precondition!(C1,A1,B1,nx1)
@time Dense_Precondition!(C2,A2,B2,nx1)
lu = @time lufact(A3)
@time lu[:U]\(lu[:L]\A3[lu[:p],:])
@time ConvLowBand(lu[:L],nx,5)
@time ConvLowBand(lu[:L],nx,5)
for i=1:2
    println("typeof(i): $(typeof(i))")
end

#@code_warntype ConvLowBand(lu[:L],nx,5)

=#

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
