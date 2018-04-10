"""
     Dense_Precondition!(H::Vector{Ta},J::Array{Ta,2},P::Array{Tp,2})

Preconditions the H & J to inv(P)H and inv(P)J using a sparse LU factorization
method with full pivoting. J and P must be of size nx-by-nx and H must be of
size nx. st is the inplace storage type.
"""

function LowSolveMult!(a::Array{U,2},b::Array{T,2},m::Int64,n::Int64) where {T,U<:AbstractFloat}
    # Add bounds check!
    @inbounds for i=1:m
        @inbounds for j = 1:n
            if (b[j,i] != zero(T))
                temp::T = b[j,i]
                @inbounds for k=(j+1):(n)
                    b[k,i] = b[k,i]-A[k,j]*temp
                end
            end
        end
    end
end

function UppSolveMult!(a::Array{U,2},b::Array{T,2},m::Int64,n::Int64) where {T,U<:AbstractFloat}
     # Add bounds check!
     @inbounds for i = 1:m
          @inbounds for j = n:-1:1
               if (b[j,i] != zero(T))
                    temp::T = b[j,i]
                    b[j,i] = temp/A[j,j]
                    @inbounds for k=1:(j-1)
                    b[k,i] = b[k,i]-A[k,j]*temp
                    end
               end
          end
     end
end
