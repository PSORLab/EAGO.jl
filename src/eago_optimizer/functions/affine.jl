# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/functions/nonlinear/affine.jl
# Defines buffered structures to store quadratic functions:
# AffineFunctionIneq, AffineFunctionEq, as well as the
# lower_interval_bound, interval_bound, and eliminate_fixed_variables!
# functions associated with each structure.
#############################################################################

###
### Structure definitions
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

AffineFunctionIneq() = AffineFunctionIneq(Tuple{Float64,Int}[], 0.0, 0)
function AffineFunctionIneq(f::SAF, s::LT)
    terms = map(x -> (x.coefficient, x.variable_index.value), f.terms)
    AffineFunctionIneq(terms, f.constant - s.upper, length(f.terms))
end
function AffineFunctionIneq(f::SAF, s::GT)
    terms = map(x -> (-x.coefficient, x.variable_index.value), f.terms)
    AffineFunctionIneq(terms, s.lower - f.constant, length(f.terms))
end
function AffineFunctionIneq(f::SV; is_max = false)
    if is_max
        return AffineFunctionIneq(Tuple{Float64,Int}[(-1.0, f.variable.value)], 0.0, 1)
    end
    AffineFunctionIneq(Tuple{Float64,Int}[(1.0, f.variable.value)], 0.0, 1)
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
