@enum(BranchCost, BC_INFEASIBLE, BC_INTERVAL, BC_INTERVAL_REV, BC_INTERVAL_LP, BC_INTERVAL_LP_REV)

Base.@kwdef mutable struct BranchOracle{T<:Real}
    strategy::BranchCost = BW_INTERVAL_LP
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
function BranchOracle{T}(n::Int) where T <:AbstractFloat
    BranchOracle{T}(𝛹n = ones(T,n),  𝛹p = ones(T,n),
                    δn = zeros(T,n),  δp = zeros(T,n),
                    ηn = zeros(T,n),  ηp = zeros(T,n))
end

function _variable_infeasibility(m::Optimizer, i) where T<:Real
    tsum = zero(T); tmin = typemax(T); tmax = typemin(T)
    d = m._branch_oracle
    for j in _sparsity(m, i)
        v = m._constraint_infeasiblity[j]
        tsum += v
        (v > tmax) && (tmax = v)
        (v < tmin) && (tmin = v)
    end
    return d.μ1*tsum + d.μ2*tmin + d.μ3*tmax
end

function _store_pseudocosts!(m::Optimizer, n::NodeBB)
    k = n.last_branch
    d = m._branch_oracle
    Δunit = (n.lower_bound - m._lower_objective_value)
    if n.branch_direction == BD_POS
        d.ηp[k] += 1
        d.𝛹p[k] = Δunit
        d.δp[k] = n.branch_extent
    elseif n.branch_direction == BD_NEG
        d.ηn[k] += 1
        d.𝛹n[k] = Δunit
        d.δn[k] = n.branch_extent
    end
    return
end

function _lo_extent(m, xb, k)
    !isfinite(l) && return _variable_infeasibility(m, i)

    c _branch_cost(m)
    (c == BC_INFEASIBLE)   && return _variable_infeasibility(m, i)
    l = _lower_bound(BranchVar, m, k)
    u = _upper_bound(BranchVar, m, k)
    (c == BC_INTERVAL)     && return xb - l
    (c == BC_INTERVAL_REV) && return u - xb

    xlp = _lower_solution(BranchVar, m, k)
    ρ = _cost_offset_β(m)*(u - l)
    y = max(min(xlp, u - ρ), l + ρ)
    return (c == BC_INTERVAL_LP) ? (y - l) : (u - y)
end
function _hi_extent(m, xb, k)
    !isfinite(u) && return _variable_infeasibility(m, i)

    c _branch_cost(m)
    (c == BC_INFEASIBLE)   && return _variable_infeasibility(m, i)
    l = _lower_bound(BranchVar, m, k)
    u = _upper_bound(BranchVar, m, k)
    (c == BC_INTERVAL)     && return u - xb
    (c == BC_INTERVAL_REV) && return xb - l

    xlp = _lower_solution(BranchVar, m, k)
    ρ = _cost_offset_β(m)*(u - l)
    y = max(min(xlp, u - ρ), l + ρ)
    return (c == BC_INTERVAL_LP) ? (u - y) : (y - l)
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

    is_integer_flag = _is_integer(BranchVar, m, k)
    if is_integer_flag
        l.is_integer[k] = floor(x) != l.lower_variable_bound[k]
        l.continuous = !any(l.is_integer)
        u.is_integer[k] = ceil(x) != u.upper_variable_bound[k]
        u.continuous = !any(u.is_integer)
    end
    lx = is_integer_flag ? floor(x) : x
    ux = is_integer_flag ? ceil(x)  : x
    l.upper_variable_bound[k] = lx
    u.lower_variable_bound[k] = ux
    l.branch_extent = _lo_extent(m, lx, k)
    u.branch_extent = _hi_extent(m, ux, k)

    push!(m._stack, l, u)
    m._node_repetitions = 1
    m._maximum_node_id += 2
    m._node_count += 2
    return
end
