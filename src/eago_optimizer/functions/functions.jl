# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/functions/nonlinear/functions.jl
# Defines variable info and function types.
#############################################################################

#include(joinpath(@__DIR__, "nonlinear\\auxiliary_variables.jl"))

"""
$(TYPEDEF)

An abstract super-type used for representing constraints built by EAGO's backend.
"""
abstract type AbstractEAGOConstraint end

"""

Computes the lower interval bound for `AbstractEAGOConstraint` representing an
inequality constraint.
"""
function lower_interval_bound end

"""

Computes a tuple representing the lower and upper interval bounds for a
`AbstractEAGOConstraint` representing an equality constraint.
"""
function interval_bound end

"""

Eliminate fixed variables by rearrangment or restructuring of `AbstractEAGOConstraint`.
"""
function eliminate_fixed_variables! end


###
###
### Affine function
###
###

"""
$(TYPEDEF)

Current only used for bound tightening. Stores a representation
of an affine inequality.
"""
mutable struct AffineFunctionIneq <: AbstractEAGOConstraint
    terms::Vector{Tuple{Float64,Int}}
    constant::Float64
    len::Int
end
const AFI = AffineFunctionIneq

AffineFunctionIneq() = AffineFunctionIneq(Tuple{Float64,Int}[], 0.0, 0)
function AffineFunctionIneq(f::SAF, s::LT)
    terms = map(x -> (x.coefficient, x.variable.value), f.terms)
    AffineFunctionIneq(terms, f.constant - s.upper, length(f.terms))
end
function AffineFunctionIneq(f::SAF, s::GT)
    terms = map(x -> (-x.coefficient, x.variable.value), f.terms)
    AffineFunctionIneq(terms, s.lower - f.constant, length(f.terms))
end
function AffineFunctionIneq(f::VI; is_max = false)
    if is_max
        return AffineFunctionIneq(Tuple{Float64,Int}[(-1.0, f.value)], 0.0, 1)
    end
    AffineFunctionIneq(Tuple{Float64,Int}[(1.0, f.value)], 0.0, 1)
end


"""
$(TYPEDEF)

Current only used for bound tightening. Stores a representation
of an affine equality.
"""
mutable struct AffineFunctionEq <: AbstractEAGOConstraint
    terms::Vector{Tuple{Float64,Int}}
    constant::Float64
    len::Int
end
const AFE = AffineFunctionEq

AffineFunctionEq() = AffineFunctionEq(Tuple{Float64,Int}[], 0.0, 0)
function AffineFunctionEq(func::SAF, set::ET)
    terms = map(x -> (x.coefficient, x.variable_index.value), func.terms)
    return AffineFunctionEq(terms, func.constant - set.value, length(func.terms))
end

###
### Parsing definitions
###

function eliminate_fixed_variables!(f::T, v::Vector{VariableInfo}) where T <: Union{AffineFunctionIneq,
                                                                                    AffineFunctionEq}
    deleted_count = 0
    i = 1
    while i + deleted_count <= f.len
        coeff, indx = @inbounds f.terms[i]
        variable_info = @inbounds v[indx]
        if variable_info.is_fixed
            f.constant += coeff*variable_info.lower_bound
            deleteat!(f.terms, i)
            deleted_count += 1
        else
            i += 1
        end
    end
    f.len -= deleted_count
    return nothing
end


###
###
### Quadratic function
###
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
const BQI = BufferedQuadraticIneq

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
const BQE = BufferedQuadraticEq

#=
mutable struct BufferedConvexQuadratic <: AbstractEAGOConstraint
    func::SQF
    buffer::Dict{Int, Float64}
    saf::SAF
    len::Int
end
=#

function create_buffer_dict(func::SQF)

    buffer = Dict{Int, Float64}()

    for term in func.quadratic_terms
        buffer[term.variable_1.value] = 0.0
        buffer[term.variable_2.value] = 0.0
    end

    for term in func.affine_terms
        buffer[term.variable.value] = 0.0
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

function eliminate_fixed_variables!(f::T, v::Vector{VariableInfo}) where T <: Union{BufferedQuadraticIneq,
                                                                                    BufferedQuadraticIneq}
    deleted_count = 0
    i = 1
    while i + deleted_count <= f.len
        term = f.sqf.terms[i]
        variable_info_1 = v[term.variable_1.value]
        variable_info_2 = v[term.variable_2.value]
        if variable_info_1.is_fixed && variable_info_2.is_fixed
            f.sqf.constant += coeff*variable_info_1.lower_bound*variable_info_2.lower_bound
            deleteat!(f.sqf.terms, i)
            deleted_count += 1
        else
            i += 1
        end
    end
    f.len -= deleted_count

    return nothing
end


###
###
### Cone 
###
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

function BufferedSOC(func::VECOFVAR, set::SOC)
    len = length(func.variables)
    buffer = Dict{Int, Float64}([(variable.value, 0.0) for variable in func.variables])
    saf = SAF(fill(SAT(0.0, VI(1)), len), 0.0)
    return BufferedSOC(copy(func), buffer, saf, len)
end

include("nonlinear/nonlinear.jl")
