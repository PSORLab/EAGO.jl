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
# src/eago_optimizer/functions/nonlinear/quadratic.jl
# Defines buffered structures to store quadratic functions:
# BufferedQuadraticIneq, BufferedQuadraticEq, as well as the
# lower_interval_bound, interval_bound, and eliminate_fixed_variables!
# functions associated with each structure.
#############################################################################

###
### Structure definitions
###

"""
$(TYPEDEF)

Stores a general quadratic inequality constraint with a buffer.
"""
mutable struct BufferedQuadraticIneq <: AbstractEAGOConstraint
    func::SQF
    buffer::Dict{Int, Float64}
    saf::SAF
    len::Int
end

"""
$(TYPEDEF)

Stores a general quadratic equality constraint with a buffer.
"""
mutable struct BufferedQuadraticEq <: AbstractEAGOConstraint
    func::SQF
    minus_func::SQF
    buffer::Dict{Int, Float64}
    saf::SAF
    len::Int
end

#=
mutable struct BufferedConvexQuadratic <: AbstractEAGOConstraint
    func::SQF
    buffer::Dict{Int, Float64}
    saf::SAF
    len::Int
end
=#

###
### Constructor definitions
###

function create_buffer_dict(func::SQF)

    buffer = Dict{Int, Float64}()

    for term in func.quadratic_terms
        buffer[term.variable_index_1.value] = 0.0
        buffer[term.variable_index_2.value] = 0.0
    end

    for term in func.affine_terms
        buffer[term.variable_index.value] = 0.0
    end

    return buffer
end

BufferedQuadraticIneq() = BufferedQuadraticIneq(SQF(SQT[], SAT[], 0.0), Dict{Int, Float64}(), SAF(SAT[], 0.0), 0)

function BufferedQuadraticIneq(func::SQF, set::LT)

    buffer = create_buffer_dict(func)
    saf = SAF([SAT(0.0, VI(k)) for k in keys(buffer)], 0.0)
    len = length(buffer)
    cfunc = copy(func)
    cfunc.constant -= set.upper

    return BufferedQuadraticIneq(cfunc, buffer, saf, len)
end

function BufferedQuadraticIneq(func::SQF, set::GT)

    buffer = create_buffer_dict(func)
    saf = SAF([SAT(0.0, VI(k)) for k in keys(buffer)], 0.0)
    len = length(buffer)
    cfunc = MOIU.operate(-, Float64, func)
    cfunc.constant += set.lower

    return BufferedQuadraticIneq(cfunc, buffer, saf, len)
end

BufferedQuadraticEq() = BufferedQuadraticEq(SQF(SQT[], SAT[], 0.0), SQF(SQT[], SAT[], 0.0), Dict{Int, Float64}(), SAF(SAT[], 0.0), 0)

function BufferedQuadraticEq(func::SQF, set::ET)

    buffer = create_buffer_dict(func)
    saf = SAF([SAT(0.0, VI(k)) for k in keys(buffer)], 0.0)
    len = length(buffer)
    cfunc1 = copy(func)
    cfunc1.constant -= set.value
    cfunc2 = MOIU.operate(-, Float64, func)
    cfunc2.constant += set.value

    return BufferedQuadraticEq(cfunc1, cfunc2, buffer, saf, len)
end

#=
function BufferedConvexQuadratic(f::BufferedQuadraticIneq)
    BufferedConvexQuadratic(copy(f.func), copy(f.buffer), copy(f.saf), f.len)
end
=#

###
### Parsing definitions
###

function eliminate_fixed_variables!(f::T, v::Vector{VariableInfo}) where T <: Union{BufferedQuadraticIneq,
                                                                                    BufferedQuadraticIneq}
    deleted_count = 0
    index = 1
    while i + deleted_count <= f.len
        term = @inbounds f.sqf.terms[i]
        variable_info_1 = @inbounds v[term.variable_index_1.value]
        variable_info_2 = @inbounds v[term.variable_index_2.value]
        if variable_info_1.is_fixed && variable_index_2.is_fixed
            f.sqf.constant += coeff*variable_info_1.lower_bound*variable_index_2.lower_bound
            deleteat!(f.sqf.terms, i)
            deleted_count += 1
        else
            i += 1
        end
    end
    f.len -= deleted_count

    return nothing
end
