# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimize/nonconvex_branch/stack_management.jl
# Contains the subroutines used for stack management. Namely, initialize_stack!
# select_branch_variable!, select_branch_point!, branch_node!, and fathom!.
#############################################################################
#=
function _variable_infeasibility(m::GlobalOptimizer, i::Int)
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

function _store_pseudocosts!(m::GlobalOptimizer, n::NodeBB)
    k = n.last_branch
    d = m._branch_cost
    if n.branch_direction == BD_POS
        d.ηp[k] += 1
        d.𝛹p[k] = n.lower_bound - m._lower_objective_value
        d.δp[k] = n.branch_extent
    elseif n.branch_direction == BD_NEG
        d.ηn[k] += 1
        d.𝛹n[k] = n.lower_bound - m._lower_objective_value
        d.δn[k] = n.branch_extent
    end
    return
end

function _lo_extent(m::GlobalOptimizer, xb::Float64, k::Int)
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

function _hi_extent(m::GlobalOptimizer, xb::Float64, k::Int)
    c = _branch_cost(m)
    (c == BC_INFEASIBLE)   && return _variable_infeasibility(m, k)

    l = _lower_bound(BranchVar(), m, k)
    u = _upper_bound(BranchVar(), m, k)
    (c == BC_INTERVAL)     && return u - xb
    (c == BC_INTERVAL_REV) && return xb - l

    xlp = _lower_solution(BranchVar(), m, k)
    ρ = _cost_offset_β(m)*(u - l)
    y = max(min(xlp, u - ρ), l + ρ)
    return (c == BC_INTERVAL_LP) ? (u - y) : (y - l)
end

@inline _score(x::T, y::T, μ::T) where T<:Real = (one(T) - μ)*min(x, y) + max(x, y)
@inline function _score(d::BranchCostStorage{T}, i::Int) where T<:Real
    _score(d.𝛹n[i]*d.δn[i], d.𝛹p[i]*d.δp[i], d.μ_score)
end

function _select_branch_variable_cost(m::GlobalOptimizer)
    return map_argmax(i -> score(m.branch_cost, i), 1:_variable_num(BranchVar(),m))
end
=#

# In the case of diam(X) = Inf variables rel_diam may equal NaN, this should be set to
# a zero relative diameter to prevent "branching" on these variables for a functionally infinite
# amount of time
"""
$(TYPEDSIGNATURES)

Return the relative diameter of a variable. In the case of diam(X)=Inf, `rel_diam`
returns 0.0 to prevent "branching" on this variable endlessly.
"""
function rel_diam(m::GlobalOptimizer, i::Int)
    rd = _diam(BranchVar(), m, i)/diam(_working_variable_info(m, _bvi(m, i)))
    return isnan(rd) ? 0.0 : rd
end
function _select_branch_variable_width(m::GlobalOptimizer)
    map_argmax(i -> rel_diam(m,i), 1:_variable_num(BranchVar(), m))
end

"""
$(TYPEDSIGNATURES)

Choose a variable to branch on. A maximum relative width branching rule is used
by default.
"""
function select_branch_variable(t::ExtensionType, m::GlobalOptimizer)
    #_branch_pseudocost_on(m) && return _select_branch_variable_cost(m)
    return _select_branch_variable_width(m)
end

"""
$(TYPEDSIGNATURES)

Select a point `xb` within the domain of the `i`th branching variable. By default,
this point is a convex combination of the solution to the relaxation and the midpoint
of the node (`branch_cvx_factor*xmid + (1-branch_cvx_factor)*xsol`). If the solution
lies within `branch_offset` of a bound, then the branch point is moved to a distance
of `branch_offset` from that bound.
"""
function select_branch_point(t::ExtensionType, m::GlobalOptimizer, i)
    l = _lower_bound(BranchVar(), m, i)
    u = _upper_bound(BranchVar(), m, i)
    s = _lower_solution(BranchVar(), m, i)
    α = _branch_cvx_α(m)
    b = _branch_offset_β(m)*(u - l)
    return max(l + b, min(u - b, α*_mid(BranchVar(), m, i) + (one(Float64) - α)*s))
end

"""
$(TYPEDSIGNATURES)

Create two nodes from `current_node` and store them on the stack. Call
`select_branch_variable(t, m)` and `select_branch_point(t, m, k)` to determine
the variable that should be branched on and the point at which branching should
occur, respectively.
"""
function branch_node!(t::ExtensionType, m::GlobalOptimizer)

    k = select_branch_variable(t, m)
    x = select_branch_point(t, m, k)
    n = m._current_node

    #isfinite(n.last_branch) && _store_pseudocosts!(m, n)

    l_bound = max(n.lower_bound, m._lower_objective_value)
    u_bound = min(n.upper_bound, m._upper_objective_value)

    l_lbd = copy(n.lower_variable_bounds);  u_lbd = copy(n.lower_variable_bounds)
    l_ubd = copy(n.upper_variable_bounds);  u_ubd = copy(n.upper_variable_bounds)
    l_int = copy(n.is_integer);             u_int = copy(n.is_integer)

    flag = is_integer(BranchVar(), m, k)
    if flag
        l_int[k] = floor(x) != n.lower_variable_bounds[k]
        u_int[k] = ceil(x)  != n.upper_variable_bounds[k]
    end
    l_cont = flag ? !any(l_int) : true
    u_cont = flag ? !any(u_int) : true
    lx = flag ? floor(x) : x
    ux = flag ? ceil(x)  : x
    l_ubd[k] = lx
    u_lbd[k] = ux

    #psuedo_cost_flag = _branch_pseudocost_on(m)
    l_ext = zero(Float64) #psuedo_cost_flag ? _lo_extent(m, lx, k) : zero(Float64)
    u_ext = zero(Float64) #psuedo_cost_flag ? _hi_extent(m, ux, k) : zero(Float64)

    l_idepth = (flag && l_cont) ? n.depth + 1 : n.cont_depth
    u_idepth = (flag && u_cont) ? n.depth + 1 : n.cont_depth

    n1 = NodeBB(l_lbd, l_ubd, l_int, l_cont, l_bound, u_bound, n.depth + 1, l_idepth, n.id + 1, BD_NEG, k, l_ext)
    n2 = NodeBB(u_lbd, u_ubd, u_int, u_cont, l_bound, u_bound, n.depth + 1, u_idepth, n.id + 2, BD_POS, k, u_ext)
    push!(m._stack, n1)
    push!(m._stack, n2)

    m._node_repetitions = 1
    m._maximum_node_id += 2
    m._node_count += 2
    return
end
branch_node!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = branch_node!(_ext(m), m)

"""
$(TYPEDSIGNATURES)

Store the current node to the stack, without branching, after updating lower/upper bounds.
"""
function single_storage!(t::ExtensionType, m::GlobalOptimizer)
    y = m._current_node
    m._node_repetitions += 1
    m._node_count += 1
    lower_bound = max(y.lower_bound, m._lower_objective_value)
    upper_bound = min(y.upper_bound, m._upper_objective_value)
    push!(m._stack, NodeBB(y.lower_variable_bounds, y.upper_variable_bounds, y.is_integer, y.continuous,
                           lower_bound, upper_bound, y.depth, y.cont_depth, y.id, y.branch_direction, y.last_branch, y.branch_extent))
    return
end
single_storage!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = single_storage!(_ext(m), m)

"""
$(TYPEDSIGNATURES)

Select the next node in the stack to evaluate. By default, perform best-first 
node selection (select the node with the lowest lower bound in the stack).
"""
function node_selection!(t::ExtensionType, m::GlobalOptimizer)
    m._node_count -= 1
    m._current_node = popmin!(m._stack)
    return
end
node_selection!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = node_selection!(_ext(m), m)

"""
$(TYPEDSIGNATURES)

Remove nodes from the stack. By default, delete nodes from the stack if their
lower bounds are greater than the current global upper bound.
"""
function fathom!(t::ExtensionType, m::GlobalOptimizer)
    u = m._global_upper_bound
    continue_flag = !isempty(m._stack)
    while continue_flag
        n = maximum(m._stack)
        max_check = n.lower_bound > u
        if max_check
            popmax!(m._stack)
            m._node_count -= 1
        end
        continue_flag = !isempty(m._stack) && max_check
    end
    return
end
fathom!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = fathom!(_ext(m), m)

"""
$(TYPEDSIGNATURES)

Check the optimization problem for unbounded branching variables, which would interfere
with EAGO's branch-and-bound routine since there are no well-defined branching rules
for cases where the interval bounds contain `-Inf` or `Inf`. If any branching variables
are missing bounds, add the missing bound at +/- 1E10 and warn the user.
"""
function unbounded_check!(m::GlobalOptimizer)
    if m._parameters.unbounded_check
        unbounded_flag = false
        wp = m._working_problem
        epigraph_flag = _variable_num(FullVar(), m) != m._input_problem._variable_count
        for i = 1:_variable_num(BranchVar(), m) - epigraph_flag #Not including epigraph reformulation variable
            if !wp._variable_info[i].has_lower_bound
                unbounded_flag = true
                wp._variable_info[i] = VariableInfo(wp._variable_info[i], GT(-1E10))
            end
            if !wp._variable_info[i].has_upper_bound
                unbounded_flag = true
                wp._variable_info[i] = VariableInfo(wp._variable_info[i], LT(1E10))
            end
        end
        unbounded_flag && @warn("""
        At least one branching variable is unbounded. This will interfere with EAGO's global
        optimization routine and may cause unexpected results. Bounds have been automatically
        generated at +/- 1E10 for all unbounded variables, but tighter user-defined bounds are
        highly recommended. To disable this warning and the automatic generation of bounds, use
        the option `unbounded_check = false`.""")
    end
end

"""
$(TYPEDSIGNATURES)

Prepare the stack for the branch-and-bound routine. By default, create an
initial node with the variable bounds as box constraints and add it to the stack.
"""
function initialize_stack!(t::ExtensionType, m::GlobalOptimizer)
    unbounded_check!(m)
    d = _working_variable_info.(m, m._branch_to_sol_map)
    push!(m._stack, NodeBB(lower_bound.(d), upper_bound.(d), is_integer.(d)))
    m._node_count = 1
    m._maximum_node_id += 1
    return
end
initialize_stack!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = initialize_stack!(_ext(m), m)
