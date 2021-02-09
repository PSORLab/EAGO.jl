# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/subsolver_config/config/clp.jl
# Configuration adjustment subroutines for Clp.
#############################################################################

function set_default_config!(::Val{:clp}, ext::ExtensionType, m::T) where T <: MOI.AbstractOptimizer

    MOI.set(m, MOI.RawParameter("PrimalTolerance"),       1E-7)
    MOI.set(m, MOI.RawParameter("DualTolerance"),         1E-7)
    MOI.set(m, MOI.RawParameter("DualObjectiveLimit"),    1e308)
    MOI.set(m, MOI.RawParameter("MaximumIterations"),     2147483647)
    MOI.set(m, MOI.RawParameter("PresolveType"),          0)
    MOI.set(m, MOI.RawParameter("SolveType"),             5)
    MOI.set(m, MOI.RawParameter("InfeasibleReturn"),      1)
    MOI.set(m, MOI.RawParameter("Scaling"),               3)
    MOI.set(m, MOI.RawParameter("Perturbation"),          100)

    return nothing
end
