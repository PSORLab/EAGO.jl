# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/subsolver_config/config/ipopt.jl
# Configuration adjustment subroutines for Ipopt.
#############################################################################

function set_default_config!(::Val{:ipopt}, ext::ExtensionType, m::T) where {T <: MOI.AbstractOptimizer}

    MOI.set(m, MOI.RawParameter("max_iter"),3000)
    MOI.set(m, MOI.RawParameter("acceptable_tol"), 1E30)
    MOI.set(m, MOI.RawParameter("acceptable_iter"), 300)
    MOI.set(m, MOI.RawParameter("constr_viol_tol"), 0.000001)
    MOI.set(m, MOI.RawParameter("acceptable_compl_inf_tol"), 0.000001)
    MOI.set(m, MOI.RawParameter("acceptable_dual_inf_tol"), 1.0)
    MOI.set(m, MOI.RawParameter("acceptable_constr_viol_tol"), 0.000001)
    MOI.set(m, MOI.RawParameter("print_level"), 0)

    return nothing
end
