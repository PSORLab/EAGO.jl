function is_integer_feasible(m::GlobalOptimizer)
    bool = true
    atol = _integer_abs_tol(m)
    rtol = _integer_rel_tol(m)
    for i = 1:_variable_num(BranchVar(), m)
        if is_integer(BranchVar(), m, i)
            xsol = _lower_solution(BranchVar(), m, i)
            if isapprox(floor(xsol), xsol; atol = atol, rtol = rtol)
                continue
            elseif isapprox(ceil(xsol), xsol; atol = atol, rtol = rtol)
                continue
            else
                bool &= false
                break
            end
        end
    end
    return bool
end

"""
$(SIGNATURES)

Default check to see if the upper bounding problem should be run. By default,
The upper bounding problem is run on every node up to depth `upper_bounding_depth`
and is triggered with a probability of `0.5^(depth - upper_bounding_depth)`
afterwards.
"""
function default_nlp_heurestic(m::GlobalOptimizer)
    bool = false
    ubd_limit = m._parameters.upper_bounding_depth
    n = _current_node(m)
    if is_integer_feasible(m)
        Δdepth = n.depth - n.cont_depth
        bool |= (Δdepth <= ubd_limit)
        bool |= (rand() < 0.5^(Δdepth - ubd_limit))
    end
    return bool
end

"""
$(SIGNATURES)

Default upper bounding problem which simply calls `solve_local_nlp!` to solve
the nlp locally.
"""
function upper_problem!(t::ExtensionType, m::GlobalOptimizer)
    if !default_nlp_heurestic(m)
        m._upper_feasibility = false
        m._upper_objective_value = Inf
    else
        solve_local_nlp!(m)
    end
    return
end
upper_problem!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = upper_problem!(_ext(m), m)
