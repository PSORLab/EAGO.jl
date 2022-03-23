
"""
$(SIGNATURES)

Default postprocess perfoms duality-based bound tightening on the `y`.
"""
function postprocess!(t::ExtensionType, m::GlobalOptimizer)
    if m._parameters.dbbt_depth > m._iteration_count
        variable_dbbt!(m._current_node, m._lower_lvd, m._lower_uvd,
                       m._lower_objective_value, m._global_upper_bound,
                       m._branch_variable_count)
    end
    return
end
postprocess!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = postprocess!(_ext(m), m)

"""
$(SIGNATURES)

Checks to see if current node should be reprocessed.
"""
repeat_check(t::ExtensionType, m::GlobalOptimizer) = false
repeat_check(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = repeat_check(_ext(m), m)
