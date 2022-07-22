# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimize/nonconvex.jl
# Functions which determine when the upper bounding (local NLP) problem should
# be solved as well as routines used to call the upper bounding problem.
#############################################################################

"""
$(TYPEDSIGNATURES)

Check that the solution of the lower (relaxed problem) is integer feasible to
within tolerances specified by the parameters: `integer_abs_tol` (absolute tolerance)
and `integer_rel_tol` (relative tolerance).
"""
function is_integer_feasible_relaxed(m::GlobalOptimizer)
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
$(TYPEDSIGNATURES)

Default check to see if the upper bounding problem should be run. By default,
The upper bounding problem is run on every node up to depth `upper_bounding_depth`
and is triggered with a probability of `0.5^(depth - upper_bounding_depth)`
afterwards for continuous problems. For integral problems, the `upper_bounding_depth`
approach is used as well as running on every node up to depth 
`upper_bounding_depth + cont_depth` with another trigger of probability
`0.5^(depth - upper_bounding_depth - cont_depth)`.
"""
function default_upper_heurestic(m::GlobalOptimizer)
    bool = false
    ubd_limit = m._parameters.upper_bounding_depth
    n = _current_node(m)
    if is_integer_feasible_relaxed(m)
        Δdepth = n.depth - n.cont_depth
        bool |= (Δdepth <= ubd_limit)
        bool |= (rand() < 0.5^(Δdepth - ubd_limit))
    end
    bool |= (n.depth <= ubd_limit)
    bool |= (rand() < 0.5^(n.depth - ubd_limit))
    return bool
end

"""
$(TYPEDSIGNATURES)

Default upper bounding problem which simply calls `solve_local_nlp!` to solve
the NLP locally.
"""
function upper_problem!(t::ExtensionType, m::GlobalOptimizer)
    if !default_upper_heurestic(m)
        m._upper_feasibility = false # Ensures that global upper bound not updated
        m._upper_objective_value = Inf
    else
        solve_local_nlp!(m)
    end
    return
end
upper_problem!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = upper_problem!(_ext(m), m)
