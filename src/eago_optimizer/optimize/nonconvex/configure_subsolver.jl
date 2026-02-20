# Copyright (c) 2018 Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Alireza Miraliakbar, Matthew Stuber, and the University of Connecticut (UConn)
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/optimize/nonconvex/configure_subsolver.jl
# Contains subroutines used to set default configuration for select supported
# solvers along with routines needed to adjust tolerances to mirror tolerance
# adjustments in the global solve.
################################################################################

function set_default_config_udf!(s, m::MOI.AbstractOptimizer, verbosity::Int)
    if verbosity > 0
        @warn("""
        EAGO lacks a specialized configuration routine for the subsolver ($(MOI.get(m, MOI.SolverName())))
        you selected. As a result, EAGO cannot set the subsolver tolerances based on the
        absolute_tolerance, relative_tolerance, and absolute_constraint_feas_tolerance
        parameters passed to the EAGO optimizer. Consequently, you need to ensure that the tolerances
        set in the provided subsolver are appropriate (for instance if the absolute_tolerance = 1E-3
        then the absolute tolerance for a subsolver should be < 1E-4 and any feasibility tolerances
        should be as conservative as the absolute_constraint_feas_tolerance). If you see this message
        please submit an issue at https://github.com/PSORLab/EAGO.jl/issues/new/choose requesting
        that a configuration routine be added for this subsolver.""")
    end
    return
end

function set_default_config!(ext, d::GlobalOptimizer, m::MOI.AbstractOptimizer, local_solver)
    set_default_config_udf!(MOI.get(m, MOI.SolverName()), m, _verbosity(d))
end

function set_default_subsolver_config!(ext::DefaultExt, d::GlobalOptimizer,  m::T, local_solver::Bool) where T
    !_user_solver_config(d) && set_default_config!(ext, d, m, local_solver)
    MOI.set(m, MOI.Silent(), true)
    nothing
end

function set_default_config!(ext::DefaultExt, m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} 
    set_default_subsolver_config!(DefaultExt(), m, m._subsolvers.relaxed_optimizer, false)
    set_default_subsolver_config!(DefaultExt(), m, m._subsolvers.upper_optimizer, true)
end

"""
    set_default_config!

Configures subsolver tolerances based on tolerance parameters provided to
EAGO (provided that a specialized subsolver configuration routine has been
provided and `m.user_solver_config = false`).
"""
set_default_config!(ext::ExtensionType, m::GlobalOptimizer) = set_default_config!(DefaultExt(), m)
set_default_config!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = set_default_config!(_ext(m), m)