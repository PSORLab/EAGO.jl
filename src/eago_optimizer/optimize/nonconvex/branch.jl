function _variable_infeasibility(m::Optimizer, i::Int)
    tsum = zero(Float64); tmin = typemax(Float64); tmax = typemin(Float64)
    d = m._branch_cost
    for j in _sparsity(BranchVar(), m, i)
        v = m._constraint_infeasiblity[j]
        tsum += v
        (v > tmax) && (tmax = v)
        (v < tmin) && (tmin = v)
    end
    return d.μ1*tsum + d.μ2*tmin + d.μ3*tmax
end

function _store_pseudocosts!(m::Optimizer, n::NodeBB)
    k = n.last_branch
    d = m._branch_cost
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

function _lo_extent(m::Optimizer, xb::Float64, k::Int)
    !isfinite(l) && return _variable_infeasibility(m, k)

    c = _branch_cost(m)
    (c == BC_INFEASIBLE)   && return _variable_infeasibility(m, k)

    l = _lower_bound(BranchVar(), m, k)
    u = _upper_bound(BranchVar(), m, k)
    (c == BC_INTERVAL)     && return xb - l
    (c == BC_INTERVAL_REV) && return u - xb

    xlp = _lower_solution(BranchVar(), m, k)
    ρ = _cost_offset_β(m)*(u - l)
    y = max(min(xlp, u - ρ), l + ρ)
    return (c == BC_INTERVAL_LP) ? (y - l) : (u - y)
end

function _hi_extent(m::Optimizer, xb::Float64, k::Int)
    !isfinite(u) && return _variable_infeasibility(m, k)

    c = _branch_cost(m)
    (c == BC_INFEASIBLE)   && return _variable_infeasibility(m, k)

    l = _lower_bound(BranchVar(), m, k)
    u = _upper_bound(BranchVar(), m, k)
    (c == BC_INTERVAL)     && return u - xb
    (c == BC_INTERVAL_REV) && return xb - l

    xlp = _lower_solution(BranchVar, m, k)
    ρ = _cost_offset_β(m)*(u - l)
    y = max(min(xlp, u - ρ), l + ρ)
    return (c == BC_INTERVAL_LP) ? (u - y) : (y - l)
end

@inline _score(x::T, y::T, μ::T) where T<:Real = (one(T) - μ)*min(x, y) + max(x, y)
@inline function _score(d::BranchCostStorage{T}, i::Int) where T<:Real
    _score(d.𝛹n[i]*d.δn[i], d.𝛹p[i]*d.δp[i], d.μ_score)
end

function _select_branch_variable_cost(m::Optimizer)
    return map_argmax(i -> score(m.branch_cost, i), 1:_branch_variable_num(m))
end

function _select_branch_variable_width(m::Optimizer)
    function branch_goal(x)
        v = m._working_problem._variable_info[_bv(m,x)]
        return (v.branch_on === BRANCH) ? _diam(BranchVar(),m,x)/diam(v) : -Inf
    end
    return map_argmax(branch_goal, 1:_branch_variable_num(m))
end

"""
$(SIGNATURES)

Selects the variable to branch on psuedocost branching is used if
(parameter: `branch_pseudocost_on` = true).
"""
function select_branch_variable(t::ExtensionType, m::Optimizer)
    _branch_pseudocost_on(m) && return _select_branch_variable_cost(m)
    return _select_branch_variable_width(m)
end

"""
$(SIGNATURES)

Selects a point `xb` which is a convex combination (parameter:
`branch_cvx_factor`) of the solution to the relaxation and the midpoint of the
node. If this solution lies within (parameter: `branch_offset`) of a bound then
the branch point is moved to a distance of `branch_offset` from the bound.
"""
function select_branch_point(t::ExtensionType, m::Optimizer, i)
    l = _lower_bound(BranchVar(), m, i)
    u = _upper_bound(BranchVar(), m, i)
    s = _lower_solution(BranchVar(), m, i)
    α = _branch_cvx_α(m)
    b = _branch_offset_β(m)*(u - l)
    return max(l + b, min(u - b, α*s + (one(T) - α)*_mid(BranchVar(), m, i)))
end

"""
$(SIGNATURES)

Creates two nodes from `current_node` and stores them to the stack. Calls
`select_branch_variable(t, m)` and `select_branch_point(t, m, k)`.
"""
function branch_node!(t::ExtensionType, m::Optimizer)

    k = select_branch_variable(t, m)
    x = select_branch_point(t, m, k)
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
