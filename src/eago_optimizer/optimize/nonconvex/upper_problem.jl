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
    depth = m._current_node.depth
    bool |= (depth <= ubd_limit)
    bool |= (rand() < 0.5^(depth - ubd_limit))
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
upper_problem!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = upper_problem!(_ext_typ(m), m)
