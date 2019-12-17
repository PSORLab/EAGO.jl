"""
$(FUNCTIONNAME)

Dense LU preconditioner for implicit McCormick relaxation.
"""
struct DenseMidInv <: AbstractPrecondionerMC
    Y::Array{Float64,2}
    YInterval::Vector{Interval{Float64}}
    nx::Int
    np::Int
end
DenseMidInv(nx::Int, np::Int) = DenseMidInv(zeros(Float64,nx,nx), zeros(Interval{Float64},1), nx, np)
function precondition!(d::DenseMidInv, H::Vector{MC{N,T}}, J::Array{MC{N,T},2}) where {N, T <: RelaxTag}
    for i in eachindex(J)
        @inbounds d.YInterval[1] = J[i].Intv
        @inbounds d.Y[i] = 0.5*(d.YInterval[1].lo + d.YInterval[1].hi)
    end
    F = lu!(d.Y)
    H .= F\H
    J .= F\J
    return
end
function (d::DenseMidInv)(h!::FH, hj!::FJ, nx::Int, np::Int) where {FH <: Function, FJ <: Function}
    return DenseMidInv(zeros(nx,nx), zeros(Interval{Float64},1), nx, np)
end
function preconditioner_storage(d::DenseMidInv, tag::T) where T <: RelaxTag
    zeros(MC{d.np, T}, (d.nx, d.nx))
end
