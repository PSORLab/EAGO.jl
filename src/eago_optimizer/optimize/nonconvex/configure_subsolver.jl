# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/subsolver_config/config.jl
# Contains subroutines used to set default configuration for select supported
# solvers along with routines needed to adjust tolerances to mirror tolerance
# adjustments in the global solve.
#############################################################################

const SUPPORTED_SUBSOLVER_NAMES = Dict{String,Symbol}("GLPK" => :glpk,
                                                      "Clp"  => :clp,
                                                      "COIN Branch-and-Cut (Cbc)" => :cbc,
                                                      "Ipopt" => :ipopt)

function set_default_config_udf!(s::String, m::Optimizer)
    if _verbosity(m) > 0
        println("EAGO lacks a specialized configuration routine for the subsolver ($s)")
        println("you selected. As a result, EAGO cannot set the subsolver tolerances based on the")
        println("absolute_tolerance, relative tolerance, and absolute_constraint_feas_tolerance")
        println("parameters passed to the EAGO optimizer. Consequently, need to ensure that the tolerances")
        println("set in the provided subsolver are appropriate (for instance if the absolute_tolerance = 1E-3")
        println("then the absolute tolerance for a subsolver should be < 1E-4 and any feasibility tolerances")
        println(" should be as conservative as the absolute_constraint_feas_tolerance). If you see this message")
        println("please submit an issue at https://github.com/PSORLab/EAGO.jl/issues/new/choose requesting")
        println("that a configuration routine be added for this subsolver.")
    end
    return
end

function set_default_config!(::Val{:clp}, ext::ExtensionType, d::Optimizer, m::T, local_solver::Bool) where T
    MOI.set(m, MOI.RawParameter("PrimalTolerance"),       _absolute_tol(d)*1E-2)
    MOI.set(m, MOI.RawParameter("DualTolerance"),         _absolute_tol(d)*1E-2)
    MOI.set(m, MOI.RawParameter("DualObjectiveLimit"),    1e308)
    MOI.set(m, MOI.RawParameter("MaximumIterations"),     2147483647)
    MOI.set(m, MOI.RawParameter("PresolveType"),          0)
    MOI.set(m, MOI.RawParameter("SolveType"),             5)
    MOI.set(m, MOI.RawParameter("InfeasibleReturn"),      1)
    MOI.set(m, MOI.RawParameter("Scaling"),               3)
    MOI.set(m, MOI.RawParameter("Perturbation"),          100)
    return
end

function set_default_config!(::Val{:cbc}, ext::ExtensionType, d::Optimizer, m::T, local_solver::Bool) where T
    MOI.set(m, MOI.RawParameter("allowableGap"), _absolute_tol(d)*1E-2)
    MOI.set(m, MOI.RawParameter("ratioGap"),     _absolute_tol(d)*1E-2)
    #MOI.set(m, MOI.RawParameter("threads"), Threads.nthreads())
    return
end

function set_default_config!(::Val{:glpk}, ext::ExtensionType, d::Optimizer, m::T, local_solver::Bool) where T
    return
end

function set_default_config!(::Val{:ipopt}, ext::ExtensionType, d::Optimizer, m::T, local_solver::Bool) where T
    c_tol = _constraint_tol(d)*1E-3
    MOI.set(m, MOI.RawParameter("tol"),_absolute_tol(d)*1E-3)
    MOI.set(m, MOI.RawParameter("print_level"), 0)
    MOI.set(m, MOI.RawParameter("constr_viol_tol"), c_tol)
    if local_solver
        MOI.set(m, MOI.RawParameter("max_iter"),3000)
        MOI.set(m, MOI.RawParameter("acceptable_tol"), 1E30)
        MOI.set(m, MOI.RawParameter("acceptable_iter"), 300)
        MOI.set(m, MOI.RawParameter("acceptable_compl_inf_tol"), c_tol)
        MOI.set(m, MOI.RawParameter("acceptable_dual_inf_tol"), 1.0)
        MOI.set(m, MOI.RawParameter("acceptable_constr_viol_tol"), c_tol)
    else
        MOI.set(m, MOI.RawParameter("max_iter"), 1E5)
        MOI.set(m, MOI.RawParameter("acceptable_iter"), 1E5+1)
    end
    return
end

function set_default_subsolver_config!(ext::DefaultExt, d::Optimizer,  m::T, local_solver::Bool) where T
    sname = MOI.get(m, MOI.SolverName())
    if (sname in keys(SUPPORTED_SUBSOLVER_NAMES)) && !_user_solver_config(d)
        sym = SUPPORTED_SUBSOLVER_NAMES[sname]
        set_default_config!(Val(sym), ext, d, m, local_solver)
    elseif !_user_solver_config(d)
        set_default_config_udf!(sname, m)
    end
    MOI.set(d, MOI.Silent())
    return
end
function set_default_config!(ext::DefaultExt, m::T) where T
    set_default_subsolver_config!(DefaultExt(), m, m.relaxed_optimizer, false)
    set_default_subsolver_config!(DefaultExt(), m, m.upper_optimizer, true)
end

"""
    set_default_config!

Configures subsolver tolerances based on tolerance parameters provided to
EAGO (provided that a specialized subsolver configuration routine has been
provided and `m.user_solver_config = false`).
"""
function set_default_config!(ext::ExtensionType, m::T) where T
    set_default_config!(DefaultExt(), m)
end
function set_default_config!(m::T) where T
    set_default_config!(m.ext_type, m)
end
