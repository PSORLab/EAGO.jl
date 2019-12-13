"""
$(FUNCTIONNAME)

Dense LU preconditioner for implicit McCormick relaxation.
"""
mutable struct DenseMidInv{T<:RelaxTag} <: AbstractPrecondionerMC
    Y::Array{Float64,2}
    F::LU{Float64,Array{Float64,2}}
    nx::Int
    np::Int
    tag::T
end
function precondition!(d::DenseMidInv{T}, H::Vector{MC{N,T}}, J::Array{MC{N,T},2}) where {N, T <: RelaxTag}
    @. d.Y = 0.5*(lo(J) + hi(J))
    d.F = lu!(d.Y)
    H[:] = d.F\H
    J[:,:] = d.F\J
    return
end
function (d::DenseMidInv{T})(h!::FH, hj!::FJ, nx::Int, np::Int, tag::T) where {FH <: Function, FJ <: Function, T<:RelaxTag}
    return DenseMidInv{T}(zeros(nx,nx), lu!(diagm(fill(1.0,(nx,)))), nx, np, tag)
end
function preconditioner_storage(d::DenseMidInv{T}) where T <: RelaxTag
    zeros(MC{d.np, T}, (d.nx, d.nx))
end
