# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimize/nonconvex/postprocess.jl
# Routine which calls duality-based bound tightening after solving the upper
# bounding problem.
#############################################################################

"""
$(TYPEDSIGNATURES)

Default postprocess perfoms duality-based bound tightening (Tawarmalani, M., 
Sahinidis, N.V.: Global optimization of mixed-integer nonlinear programs: 
a theoretical and computational study. Math. Progr. 99, 563–591 (2004).) up to
an iteration limit set by `dbbt_depth`.
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
$(TYPEDSIGNATURES)

Check to see if current node should be reprocessed. Without any custom extension,
return `false` by default.
"""
repeat_check(t::ExtensionType, m::GlobalOptimizer) = false
repeat_check(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = repeat_check(_ext(m), m)
