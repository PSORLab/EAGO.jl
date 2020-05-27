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

Stores a general quadratic function with a buffer.
"""
mutable struct BufferedQuadraticIneq <: AbstractEAGOConstraint
    func::SQF
    buffer::Dict{Int, Float64}
    saf::SAF
    len::Int
end

"""
$(TYPEDEF)

Stores a general quadratic function with a buffer.
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
    return BufferQuadraticIneq(cfunc, buffer, saf, len)
end

function BufferedQuadraticIneq(func::SQF, set::GT)
    buffer = create_buffer_dict(func)
    saf = SAF([SAT(0.0, VI(k)) for k in keys(buffer)], 0.0)
    len = length(buffer)
    cfunc = MOIU.operate(-, Float64, func)
    cfunc.constant += set.lower
    BufferQuadraticIneq(cfunc, buffer, saf, len)
end

BufferedQuadraticEq() = BufferedQuadraticEq(SQF(SQT[], SAT[], 0.0), Dict{Int, Float64}(), SAF(SAT[], 0.0), 0)

function BufferedQuadraticEq(func::SQF, set::ET)
    buffer = create_buffer_dict(func)
    saf = SAF([SAT(0.0, VI(k)) for k in keys(buffer)], 0.0)
    len = length(buffer)
    cfunc1 = copy(func)
    cfunc1.constant -= set.value
    cfunc2 = MOIU.operate(-, Float64, func)
    cfunc2.constant += set.value
    BufferQuadratic(cfunc1, cfunc2, buffer, saf, len)
end

#=
function BufferedConvexQuadratic(f::BufferedQuadraticIneq)
    BufferedConvexQuadratic(copy(f.func), copy(f.buffer), copy(f.saf), f.len)
end
=#

###
### Interval bounding definitions
###

function lower_interval_bound(f::BufferedQuadraticIneq, n::NodeBB)

    lo_bnds = n.lower_variable_bounds
    up_bnds = n.upper_variable_bounds
    lower_interval_bound = Interval{Float64}(f.func.constant)

    for aff_term in f.func.affine_terms
        coeff = aff_term.coefficient
        vi = aff_term.variable_index.value
        @inbounds xL = lo_bnds[vi]
        @inbounds xU = up_bnds[vi]
        lower_interval_bound += coeff > 0.0 ? coeff*xL : coeff*xU
    end

    for quad_term in f.func.quadratic_terms
        coeff = quad_term.coefficient
        vi1 = quad_term.variable_index_1.value
        vi2 = quad_term.variable_index_2.value
        @inbounds xL = lo_bnds[vi1]
        @inbounds xU = up_bnds[vi1]
        if vi1 === vi2
            if coeff > 0.0
                lower_interval_bound += (0.0 < xL) ? coeff*xL*xL : ((xU <= 0.0) ? coeff*xU*xU : 0.0)
            else
                lower_interval_bound += (xL < xU) ? coeff*xU*xU : coeff*xL*xL
            end
        else
            @inbounds il2b = lo_bnds[vi2]
            @inbounds iu2b = up_bnds[vi2]
            lower_interval_bound += coeff*Interval{Float64}(xL, xU)*Interval{Float64}(il2b, iu2b)
        end
    end

    return val_intv.lo, val_intv.hi
end

function interval_bound(f::BufferedQuadraticEq, n::NodeBB)

    lo_bnds = n.lower_variable_bounds
    up_bnds = n.upper_variable_bounds
    val_intv = Interval(f.func.constant)

    for aff_term in f.func.affine_terms
        coeff = aff_term.coefficient
        vi = aff_term.variable_index.value
        @inbounds il1b = lo_bnds[vi]
        @inbounds iu1b = up_bnds[vi]
        val_intv += coeff*Interval(il1b, iu1b)
    end

    for quad_term in f.func.quadratic_terms
        coeff = quad_term.coefficient
        vi1 = quad_term.variable_index_1.value
        vi2 = quad_term.variable_index_2.value
        @inbounds il1b = lo_bnds[vi1]
        @inbounds iu1b = up_bnds[vi1]
        if vi1 === vi2
            val_intv += coeff*pow(Interval(il1b, iu1b), 2)
        else
            @inbounds il2b = lo_bnds[vi2]
            @inbounds iu2b = up_bnds[vi2]
            val_intv += coeff*Interval(il1b, iu1b)*Interval(il2b, iu2b)
        end
    end

    return val_intv.lo, val_intv.hi
end

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
