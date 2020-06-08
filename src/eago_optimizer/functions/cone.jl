# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/functions/nonlinear/cone.jl
# Defines buffered structures to store quadratic functions:
# BufferedSOC as well as the  lower_interval_bound, interval_bound, and
# eliminate_fixed_variables! functions associated with each structure.
#############################################################################

###
### Structure definitions
###
"""
$(TYPEDEF)

Stores a second-order cone with a buffer.
"""
mutable struct BufferedSOC <: AbstractEAGOConstraint
    variables::VECOFVAR
    buffer::Dict{Int, Float64}
    saf::SAF
    len::Int
end

###
### Constructor definitions
###
function BufferedSOC(func::VECOFVAR, set::SOC)
    len = length(func.variables)
    buffer = Dict{Int, Float64}([(variable.value, 0.0) for variable in func.variables])
    saf = SAF(fill(SAT(0.0, VI(1)), len), 0.0)
    return BufferedSOC(copy(func), buffer, saf, len)
end

###
### Interval bounding definitions
###

###
### Parsing definitions
###
