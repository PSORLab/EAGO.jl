#=
# TODO: Infer dimensionality of user function from context type
struct RelaxMeta{N,T} 
    p::SVector{N,T}
    pL::SVector{N,T}
    pU::SVector{N,T}
    v::VariableValues{Float64}
    s::Vector{Int}
    ϵ::Float64
    subgradient_refinement::Bool
end

@context RelaxCtx

function ctx_subgradient_refine(c::RelaxCtx{RelaxMeta{N,T},Nothing}, z::MC{N,T}) where {N,T}
    m = ctx.metadata
    v = ctx.metadata.v
    s = ctx.metadata.s
    ϵ = ctx.metadata.ϵ
    return m.subgradient_refinement ? set_value_post(z, m.v, m.s, m.ϵ), z
end

for k in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[k]
    @eval function Cassette.overdub(ctx::RelaxCtx{RelaxMeta{N,T},Nothing}, ::typeof($f), x::MC{N,T}) where {N,T}
        return ctx_subgradient_refine(ctx, ($f)(x))
    end
end

for k in BIVARIATE_ATOM_TYPES
    f = BIVARIATE_ATOM_DICT[k]
    @eval function Cassette.overdub(ctx::RelaxCtx{RelaxMeta{N,T},Nothing}, ::typeof($f), x::MC{N,T}, y::NumberNotRelax) where {N,T}
        return ctx_subgradient_refine(ctx, ($f)(x, y))
    end
    @eval function Cassette.overdub(ctx::RelaxCtx{RelaxMeta{N,T},Nothing}, ::typeof($f), y::NumberNotRelax, x::MC{N,T}) where {N,T}
        return ctx_subgradient_refine(ctx, ($f)(y, x))
    end
end

#=
for k in NARITY_ATOM_TYPES
    f = NARITY_ATOM_DICT[k]
    @eval function Cassette.overdub(ctx::RelaxCtx{RelaxMeta{N,T},Nothing}, ::typeof($f), ...) where {N,T}
    end
end
=#

# commonly used storage objects automatically get promoted to use the set-valued storage
Cassette.overdub(::F_RELAX_CTX{T}, ::typeof(zeros), n::Int) where {T<:Number} = zeros(T, n)
Cassette.overdub(::F_RELAX_CTX{T}, ::typeof(zeros), dims...) where {T<:Number} = zeros(T, dims...)
Cassette.overdub(::F_RELAX_CTX{T}, ::typeof(fill), v::Any, dims::Int...) where {T<:Number} = fill(convert(T,v), dims)
function Cassette.overdub(::F_RELAX_CTX{T}, ::typeof(hcat), A::AbstractArray, Bs::AbstractArray...) where {T<:Number}
    vA = hcat(A,Bs...)
    vR = zeros(T, size(vA)...)
    return copyto!(vR,vA)
end
function Cassette.overdub(::F_RELAX_CTX{T}, ::typeof(vcat), A::AbstractArray, Bs::AbstractArray...) where {T<:Number}
    vA = vcat(A,Bs...)
    vR = zeros(T, size(vA)...)
    return copyto!(vR,vA)
end
=#
#=
# solution object
function Cassette.overdub(::F_RELAX_CTX{T}, ::typeof(ODEProblem), f, u0, tspan) where {T<:Number}
end
function Cassette.overdub(::F_RELAX_CTX{T}, ::typeof(solve), prob, alg, kwags...) where {T<:Number}
end
=#
