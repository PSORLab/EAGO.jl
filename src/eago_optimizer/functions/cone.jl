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
# TODO
#############################################################################

###
### Structure definitions
###
"""
$(TYPEDEF)

Stores a second-order conen with a buffer.
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
    len = length(func)
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
