# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/subsolver_config/config/cplex.jl
# Configuration adjustment subroutines for CPLEX.
#############################################################################

function set_default_config!(::Val{:cplex}, ext::ExtensionType, m::T) where {T <: MOI.AbstractOptimizer}

    return nothing
end
