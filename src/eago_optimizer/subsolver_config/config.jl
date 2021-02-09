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

function set_default_config_udf!(ext::ExtensionType, m::MOI.AbstractOptimizer) end

include("config/cbc.jl")
include("config/clp.jl")
include("config/cosmo.jl")
include("config/cplex.jl")
include("config/ecos.jl")
include("config/glpk.jl")
include("config/gurobi.jl")
include("config/ipopt.jl")
include("config/knitro.jl")
include("config/mosek.jl")
include("config/scs.jl")
include("config/tulip.jl")

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
