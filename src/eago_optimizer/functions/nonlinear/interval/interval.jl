struct RelaxInterval <: AbstractCacheAttribute end

Base.@kwdef mutable struct IntervalCache{T<:Real} <: AbstractCache
    ic::IntervalCache                          = IntervalCache{T}()
    v::VariableValues{Float64}                 = VariableValues{T}()
    _set::Vector{Interval{T}}                  = Interval{T}[]
    _subexpression_set::Dict{Int,Interval{T}}  = Dict{Int,Interval{T}}()
    _set_mv_buffer::Vector{Interval{T}}        = Interval{T}[]
    first_eval::Bool                           = true
end
function IntervalCache{T<:Real}(::RelaxInterval, n::Int, m::Int, p::Int) where T<:Real
    IntervalCache{T}(_set = zeros(Interval{T}, n), _set_mv_buffer = zeros(Interval{T}, p))
end
function initialize!(c::IntervalCache{T}, g::DirectedTree) where T<:Real
    c._set                 = zeros(Interval{T}, node_count(g))
    c._set_mv_buffer       = zeros(Interval{T}, length(sparsity(g, 1)))
    return
end

set(::RelaxInterval, b::IntervalCache{T}, i) where T<:Real = b._set[i]
num(::RelaxInterval, b::IntervalCache{T}, i) where T<:Real = b._set[i].lo
is_num(::RelaxInterval, b::IntervalCache{T}, i) where T<:Real = b._set[i]
subexpression_set(::RelaxInterval, b::IntervalCache{T}, i) where T<:Real = b._subexpression_set[i]
set_mv_buffer(::RelaxInterval, b::IntervalCache{T}, i) where T<:Real = b._set_mv_buffer[i]

function Base.setindex!(b::IntervalCache{T}, v::Interval{T}, i::Int) where T<:Real
    if first_eval
        b._set[i] = v
    else
        b._set[i] = b._set[i] ∩ v
    end
    nothing
end

include(joinpath(@__DIR__, "forward_interval.jl"))