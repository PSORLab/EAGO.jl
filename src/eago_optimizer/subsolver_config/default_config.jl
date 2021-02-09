# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/subsolver_config/default_config.jl
# Contains subroutines used to set default configuration for select supported solvers.
#############################################################################

function set_default_config_udf!(ext::ExtensionType, m::MOI.AbstractOptimizer) end

include("default/cbc.jl")
include("default/clp.jl")
include("default/cosmo.jl")
include("default/cplex.jl")
include("default/ecos.jl")
include("default/glpk.jl")
include("default/gurobi.jl")
include("default/ipopt.jl")
include("default/knitro.jl")
include("default/mosek.jl")
include("default/scs.jl")
include("default/tulip.jl")

function set_default_config!(ext::DefaultExt, m::T) where {T <: MOI.AbstractOptimizer}
    solver_name = MOI.get(m, MOI.SolverName())
    if solver_name == "Ipopt"
        set_default_config!(Val{:ipopt}(), ext, m)
    elseif solver_name == "GLPK"
        set_default_config!(Val{:glpk}(), ext, m)
    elseif solver_name == "CPLEX"
        set_default_config!(Val{:cplex}(), ext, m)
    elseif solver_name == "Clp"
        set_default_config!(Val{:clp}(), ext, m)
    elseif solver_name == "COIN Branch-and-Cut (Cbc)"
        set_default_config!(Val{:cbc}(), ext, m)
    elseif solver_name == "SCS"
        set_default_config!(Val{:scs}(), ext, m)
    elseif solver_name == "Mosek"
        set_default_config!(Val{:mosek}(), ext, m)
    elseif solver_name == "ECOS"
        set_default_config!(Val{:ecos}(), ext, m)
    elseif solver_name == "COSMO"
        set_default_config!(Val{:cosmo}(), ext, m)
    elseif solver_name == "Tulip"
        set_default_config!(Val{:tulip}(), ext, m)
    elseif solver_name == "Gurobi"
        set_default_config!(Val{:gurobi}(), ext, m)
    elseif solver_name == "Knitro"
        set_default_config!(Val{:knitro}(), ext, m)
    else
        set_default_config_udf!(ext, m)
    end
    return nothing
end

function set_default_config!(ext::ExtensionType, m::T) where {T <: MOI.AbstractOptimizer}
    set_default_config!(DefaultExt(), m::T)
end
