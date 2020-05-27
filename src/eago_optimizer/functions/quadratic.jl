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

"""
$(TYPEDEF)

Stores a general quadratic function with a buffer.
"""
mutable struct BufferedQuadraticIneq <: AbstractEAGOConstraint
    func::SQF
    buffer::OrderedDict{Int, Float64}
    saf::SAF
    len::Int
end

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

"""
$(TYPEDEF)

Stores a general quadratic function with a buffer.
"""
mutable struct BufferedQuadraticEq <: AbstractEAGOConstraint
    func::SQF
    buffer::OrderedDict{Int, Float64}
    saf::SAF
    len::Int
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
