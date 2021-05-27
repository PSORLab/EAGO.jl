# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/bound.jl
# Computes interval bounds of various functions.
#############################################################################

###
### AFFINE FUNCTIONS
###
function lower_interval_bound(m::Optimizer, f::AffineFunctionIneq)
    fL = f.constant
    for (c, j) in f.terms
        xL = _lower_bound(FullVar(), m, j)
        xU = _upper_bound(FullVar(), m, j)
        fL += (c > 0.0) ? c*xL : c*xU
    end
    return fL
end

function interval_bound(m::Optimizer, f::AffineFunctionEq, y::NodeBB)
    fL = fU = f.constant
    for (c, j) in f.terms
        xL = _lower_bound(FullVar(), m, j)
        xU = _upper_bound(FullVar(), m, j)
        if c > 0.0
            fL += c*xL
            fU += c*xU
        else
            fL += c*xU
            fU += c*xL
        end
    end
    return fL, fU
end

###
### QUADRATIC FUNCTIONS
###

function lower_interval_bound(m::Optimizer, f::BufferedQuadraticIneq, n::NodeBB)

    sol_branch_map = m._sol_to_branch_map
    lo_bnds = n.lower_variable_bounds
    up_bnds = n.upper_variable_bounds
    lower_interval_bound = Interval{Float64}(f.func.constant)

    for aff_term in f.func.affine_terms

        coeff = aff_term.coefficient

        vi = aff_term.variable_index.value
        if m._branch_variables[vi]
            mapped_vi = @inbounds sol_branch_map[vi]
            xL = @inbounds lo_bnds[mapped_vi]
            xU = @inbounds up_bnds[mapped_vi]
        else
            xL = @inbounds m._working_problem._variable_info[vi].lower_bound
            xU = @inbounds m._working_problem._variable_info[vi].upper_bound
        end

        lower_interval_bound += coeff > 0.0 ? coeff*xL : coeff*xU
    end

    # assumes branching on all quadratic terms, otherwise we'd need to distinguish
    # between look ups to the node bounds and lookups to the variable info in the
    # working problem.
    for quad_term in f.func.quadratic_terms

        coeff = quad_term.coefficient

        vi1 = quad_term.variable_index_1.value
        vi2 = quad_term.variable_index_2.value
        mapped_vi1 = sol_branch_map[vi1]

        xL = @inbounds lo_bnds[mapped_vi1]
        xU = @inbounds up_bnds[mapped_vi1]

        if vi1 === vi2
            if coeff > 0.0
                lower_interval_bound += (0.0 < xL) ? 0.5*coeff*xL*xL : ((xU <= 0.0) ? 0.5*coeff*xU*xU : 0.0)
            else
                lower_interval_bound += (xL < xU) ? 0.5*coeff*xU*xU : 0.5*coeff*xL*xL
            end
        else
            mapped_vi2 = sol_branch_map[vi2]

            il2b = @inbounds lo_bnds[mapped_vi2]
            iu2b = @inbounds up_bnds[mapped_vi2]
            lower_interval_bound += coeff*Interval{Float64}(xL, xU)*Interval{Float64}(il2b, iu2b)

        end
    end

    return lower_interval_bound.lo
end

function interval_bound(m::Optimizer, f::BufferedQuadraticEq, n::NodeBB)

    sol_branch_map = m._sol_to_branch_map
    lo_bnds = n.lower_variable_bounds
    up_bnds = n.upper_variable_bounds
    val_intv = Interval(f.func.constant)

    for aff_term in f.func.affine_terms

        coeff = aff_term.coefficient

        vi = aff_term.variable_index.value
        if m._branch_variables[vi]
            mapped_vi = @inbounds sol_branch_map[vi]
            xL = @inbounds lo_bnds[mapped_vi]
            xU = @inbounds up_bnds[mapped_vi]
        else
            xL = @inbounds m._working_problem._variable_info[vi].lower_bound
            xU = @inbounds m._working_problem._variable_info[vi].upper_bound
        end

        val_intv += coeff*Interval(xL, xU)
    end

    # assumes branching on all quadratic terms, otherwise we'd need to distinguish
    # between look ups to the node bounds and lookups to the variable info in the
    # working problem.
    for quad_term in f.func.quadratic_terms

        coeff = quad_term.coefficient

        vi1 = quad_term.variable_index_1.value
        vi2 = quad_term.variable_index_2.value

        mapped_vi1 = @inbounds sol_branch_map[vi1]

        xL = @inbounds lo_bnds[mapped_vi1]
        xU = @inbounds up_bnds[mapped_vi1]

        if vi1 === vi2
            val_intv += 0.5*coeff*pow(Interval(xL, xU), 2)
        else
            mapped_vi2 = @inbounds sol_branch_map[vi2]
            @inbounds il2b = lo_bnds[mapped_vi2]
            @inbounds iu2b = up_bnds[mapped_vi2]
            val_intv += coeff*Interval(xL, xU)*Interval(il2b, iu2b)
        end
    end

    return val_intv.lo, val_intv.hi
end

###
### SECOND-ORDER CONE
###
function lower_interval_bound(m::Optimizer, d::BufferedSOC, n::NodeBB)

    sol_branch_map = m._sol_to_branch_map
    lo_bnds = n.lower_variable_bounds
    up_bnds = n.upper_variable_bounds
    vec_of_vi = d.variables.variables

    norm_bound = Interval(0.0)
    for i = 2:length(vec_of_vi)
        mapped_vi = @inbounds sol_branch_map[vec_of_vi[i].value]
        x = Interval{Float64}(lo_bnds[mapped_vi], up_bnds[mapped_vi])
        norm_bound += pow(x, 2)
    end
    norm_bound = sqrt(norm_bound)

    mapped_vi = @inbounds sol_branch_map[vec_of_vi[1].value]
    lower_bound = norm_bound.lo -(@inbounds up_bnds[mapped_vi])

    return lower_bound
end


###
### NONLINEAR FUNCTIONS
###
function lower_interval_bound(m::Optimizer, d::BufferedNonlinearFunction{V,S}, n::NodeBB) where {V,S<:Real}
    !_has_value(d) && forward_pass!(m._working_problem._relaxed_evaluator, d)
    _is_num(d) ? _num(d) : _interval(d).lo
end
function interval_bound(m::Optimizer, d::BufferedNonlinearFunction{V,S}, n::NodeBB) where {V,S<:Real}
    !_has_value(d) && forward_pass!(d.evaluator, d)
    v = _is_num(d) ? Interval{S}(_num(d)) : _interval(d)
    return v.lo, v.hi
end

function is_feasible(m::Optimizer, f::Union{AffineFunctionIneq,BufferedQuadraticIneq}, y::NodeBB)
    lower_interval_bound(m, f, y) <= 0.0
end
function is_feasible(m::Optimizer, f::Union{AffineFunctionEq,BufferedQuadraticEq}, y::NodeBB)
    l, u = lower_interval_bound(m, f, y)
    l <= 0.0 <= u
end
function is_feasible(m::Optimizer, f::BufferedNonlinearFunction{V,S}, y::NodeBB) where {V,S<:Real}
    l, u = lower_interval_bound(m, f, y)
    feasible_flag = (upper_value < _lower_bound(f))
    feasible_flag && (lower_value > _upper_bound(f))
end

bound_objective(m::Optimizer, f::BufferedNonlinearFunction, n::NodeBB) = lower_interval_bound(m, f, n)
bound_objective(m::Optimizer, f::AffineFunctionIneq, n::NodeBB) = lower_interval_bound(m, f, n)
bound_objective(m::Optimizer, f::BufferedQuadraticIneq, n::NodeBB) = lower_interval_bound(m, f, n)
bound_objective(m::Optimizer, f::SV, n::NodeBB) = _lower_bound(FullVar(), m, f.variable.value)
function bound_objective(t::ExtensionType, m::Optimizer)
    bound_objective(m, m._working_problem._objective, m._current_node)
end
bound_objective(m::Optimizer) = bound_objective(m.ext_type, m)
