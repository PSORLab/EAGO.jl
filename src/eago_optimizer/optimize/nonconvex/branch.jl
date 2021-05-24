abstract type AbstractBranchCost end
struct BranchCostInfeasible <: AbstractBranchCost end
struct BranchCostInterval <: AbstractBranchCost end
struct BranchCostIntervalRev <: AbstractBranchCost end
struct BranchCostIntervalLP <: AbstractBranchCost end
struct BranchCostIntervalLPRev <: AbstractBranchCost end

abstract type AbstractBranchAlgorithm end
struct PseudoCostBraching <: AbstractBranchAlgorithm end

Base.@kwdef mutable struct BranchOracle{T<:Real}
    strategy::BRANCH_WEIGHT_STRATEGY = BW_INTERVAL_LP
    𝛹n::Vector{T}            = T[]
    𝛹p::Vector{T}            = T[]
    δn::Vector{T}            = T[]
    δp::Vector{T}            = T[]
    ηn::Vector{T}            = T[]
    ηp::Vector{T}            = T[]
    μ1::T                    = 0.1
    μ2::T                    = 1.3
    μ3::T                    = 0.8
    β::T                     = 0.05
    μ_score::T               = 0.15
end
function BranchOracle(T::Type{<:AbstractFloat}, n)
    BranchOracle{T}(𝛹n = ones(T,n),
                    𝛹p = ones(T,n),
                    δn = zeros(T,n),
                    δp = zeros(T,n),
                    ηn = zeros(T,n),
                    ηp = zeros(T,n)
                    )
end

function _variable_infeasibility(m::Optimizer, d::BranchOracle{T}, i) where T<:Real
    tsum = zero(T); tmin = typemax(T); tmax = typemin(T)
    for j in _sparsity(m, i)
        v = m._constraint_infeasiblity[j]
        tsum += v
        (v > tmax) && (tmax = v)
        (v < tmin) && (tmin = v)
    end
    return d.μ1*tsum + d.μ2*tmin + d.μ3*tmax
end

function _set_δ!(::Val{BW_INFEASIBLE}, m::Optimizer, d::BranchOracle, i) where T<:Real
    v = _variable_infeasibility(m, d, i)
    m._variable_infeasibility[i] = v
    d.δn[i] = v
    d.δp[i] = v
    return
end

function _set_δ!(::Val{BW_INTERVAL_BRANCH}, m::GlobalOptimizer, d::BranchOracle, i) where {T<:Real}
    l = _lower_bound(BranchVar, m, i)
    u = _upper_bound(BranchVar, m, i)
    d.δn[i] = isfinite(l) ? (xb - l) : _variable_infeasibility(m, d, i)
    d.δp[i] = isfinite(u) ? (u - xb) : _variable_infeasibility(m, d, i)
    return
end

function _set_δ!(::Val{BW_INTERVAL_BRANCH_REV}, m::GlobalOptimizer, d::BranchOracle, i) where {T<:Real}
    l = _lower_bound(BranchVar, m, i)
    u = _upper_bound(BranchVar, m, i)
    d.δn[i] = isfinite(l) ? (u - xb) : _variable_infeasibility(m, d, i)
    d.δp[i] = isfinite(u) ? (xb - l) : _variable_infeasibility(m, d, i)
    return
end

function _set_δ!(::Val{BW_INTERVAL_LP}, m::GlobalOptimizer, d::BranchOracle, i) where {T<:Real}
    l = _lower_bound(BranchVar, m, i)
    u = _upper_bound(BranchVar, m, i)
    ρ = d.β*(u - l)
    xlp_adj = max(min(xlp, u - ρ), l + ρ)
    d.δn[i] = isfinite(l) ? (xlp_adj - l) : _variable_infeasibility(m, d, i)
    d.δp[i] = isfinite(u) ? (u - xlp_adj) : _variable_infeasibility(m, d, i)
    return
end

function _set_δ!(::Val{BW_INTERVAL_LP_REV}, m::GlobalOptimizer, d::BranchOracle, i) where {T<:Real}
    l = _lower_bound(BranchVar, m, i)
    u = _upper_bound(BranchVar, m, i)
    ρ = d.β*(u - l)
    xlp_adj = max(min(xlp, u - ρ), l + ρ)
    d.δn[i] = isfinite(l) ? (u - xlp_adj) : _variable_infeasibility(m, d, i)
    d.δp[i] = isfinite(u) ? (xlp_adj - l) : _variable_infeasibility(m, d, i)
    return
end

function _store_pseudocosts!(m::Optimizer, b::BranchOracle{T}, n::NodeBB)
    k = n.last_branch
    Δunit = (m._lower_objective_value - n.lower_bound)/n.branch_extent
    if n.branch_direction == BD_POS
        d.ηp[k] += 1
        d.δp[k] = Δunit
        d.𝛹p[k] = d.δp[k]/d.ηp[k]
    elseif n.branch_direction == BD_NEG
        d.ηn[k] += 1
        d.δn[k] = Δunit
        d.𝛹n[k] =  d.δn[k]/d.ηn[k]
    end
    return
end

@inline _score(x::T, y::T, μ::T) where T<:Real = (one(T) - μ)*min(x, y) + max(x, y)
@inline _score(d::BranchOracle{T}, i) where T<:Real
    _score(d.𝛹n[i]*d.δn[i], d.𝛹p[i]*d.δp[i], d.μ_score)
end
function _select_branch_variable!(m, b::BranchOracle{T}) where T <: Real
    j = 1
    s = typemin(T)
    for i = 1: _branch_variable_num(m)
        v = _score(m._branch_oracle, i)
        if v > s
            v = s
            j = i
        end
    end
    return j
end
function _select_branch_variable!(t::ExtensionType, m::Optimizer) where T<:Real
    _select_branch_variable!(m, m.branch_oracle)
end


function _select_branch_point(m::Optimizer, i)
    l = _lower_bound(BranchVar, m, i)
    u = _upper_bound(BranchVar, m, i)
    s = _lower_solution(BranchVar, m, i)
    α = _branch_cvx_α(m)
    b = _branch_offset_β(m)*(u - l)
    return max(l + b, min(u - b, α*s + (one(T) - α)*_mid(BranchVar, m, i)))
end

"""
$(SIGNATURES)

Creates two nodes from `current_node` using information available the `x`
and stores them to the stack. By default, relative width bisection is perfomed
at a point `branch_pnt` which is a convex combination
(parameter: `branch_cvx_factor`) of the solution to the relaxation and
the midpoint of the node. If this solution lies within `branch_offset/width` of
a bound then the branch point is moved to a distance of `branch_offset/width`
from the bound.
"""
function branch_node!(t::ExtensionType, m::Optimizer)

    k = _select_branch_variable(t, m)
    x = _select_branch_point!(t, m, k)
    n = _current_node(m)

    isfinite(n.last_branch) && _store_pseudocosts!(m, n)

    n.lower_bound = max(n.lower_bound, m._lower_objective_value)
    n.upper_bound = min(n.upper_bound, m._upper_objective_value)
    n.last_branch = k
    n.depth += 1

    l = NodeBB(n);  l.id += 1;  l.branch_direction = BD_NEG
    u = NodeBB(n);  u.id += 2;  u.branch_direction = BD_POS
    l.branch_extent = x - _lower_bound(BranchVar, m, k)
    u.branch_extent = _upper_bound(BranchVar, m, k) - x

    is_integer_flag = _is_integer(BranchVar, m, k)
    if is_integer_flag
        l.is_integer[k] = floor(x) != l.lower_variable_bound[k]
        l.continuous = !any(l.is_integer)
        u.is_integer[k] = ceil(x) != u.upper_variable_bound[k]
        u.continuous = !any(u.is_integer)
    end
    u.lower_variable_bound[k] = is_integer_flag ? ceil(x)  : x
    l.upper_variable_bound[k] = is_integer_flag ? floor(x) : x

    push!(m._stack, l, u)
    m._node_repetitions = 1
    m._maximum_node_id += 2
    m._node_count += 2
    return
end
