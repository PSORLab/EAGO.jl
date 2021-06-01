
# TODO: Infer dimensionality of user function from context type
struct RelaxMeta{T} end
const F_RELAX_CTX{T} = RelaxCtx{RelaxMeta{T}}

for k in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[k]
    eval(quote Cassette.overdub(F_RELAX_CTX{T}, ::typeof($f), x::T) where {T<:Number} = $f(x) end)
end
Cassette.overdub(c::F_RELAX_CTX{T}, ::typeof($f), x::T) where {T<:Number} = $f(x)

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

#=
# solution object
function Cassette.overdub(::F_RELAX_CTX{T}, ::typeof(ODEProblem), f, u0, tspan) where {T<:Number}
end
function Cassette.overdub(::F_RELAX_CTX{T}, ::typeof(solve), prob, alg, kwags...) where {T<:Number}
end
=#
