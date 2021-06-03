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
function lower_interval_bound(m::GlobalOptimizer, f::AffineFunctionIneq)
    fL = f.constant
    for (c, j) in f.terms
        xL = _lower_bound(FullVar(), m, j)
        xU = _upper_bound(FullVar(), m, j)
        fL += (c > 0.0) ? c*xL : c*xU
    end
    return fL
end

function interval_bound(m::GlobalOptimizer, f::AffineFunctionEq)
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

function lower_interval_bound(m::GlobalOptimizer, f::BufferedQuadraticIneq, n::NodeBB)
    fval = Interval{Float64}(f.func.constant)
    for t in f.func.affine_terms
        c = t.coefficient
        j = aff_term.variable_index.value
        xL = _lower_bound(FullVar(), m, j)
        xU = _upper_bound(FullVar(), m, j)
        fval += c > 0.0 ? c*xL : c*xU
    end
    for t in f.func.quadratic_terms
        c = t.coefficient
        i = t.variable_index_1.value
        j = t.variable_index_2.value
        xL = _lower_bound(FullVar(), m, i)
        xU = _upper_bound(FullVar(), m, i)
        if i == j
            if c > 0.0
                fval += (0.0 < xL) ? 0.5*c*xL*xL : ((xU <= 0.0) ? 0.5*c*xU*xU : 0.0)
            else
                fval += (xL < xU) ? 0.5*c*xU*xU : 0.5*c*xL*xL
            end
        else
            yL = _lower_bound(FullVar(), m, j)
            yU = _upper_bound(FullVar(), m, j)
            fval += c*Interval{Float64}(xL, xU)*Interval{Float64}(yL, yU)
        end
    end
    return fval.lo
end
function interval_bound(m::GlobalOptimizer, f::BufferedQuadraticIneq, n::NodeBB)

    vi = m._input_problem._variable_info
    fval = Interval{Float64}(f.func.constant)
    for t in f.func.affine_terms
        c = t.coefficient
        j = t.variable_index.value
        xL = lower_bound(vi[j])
        xU = upper_bound(vi[j])
        fval += c*Interval(xL, xU)
    end
    for t in f.func.quadratic_terms
        c = t.coefficient
        i = t.variable_index_1.value
        j = t.variable_index_2.value
        xL = lower_bound(vi[i])
        xU = upper_bound(vi[i])
        if i == j
            fval += 0.5*c*pow(Interval(xL, xU), 2)
        else
            yL = lower_bound(vi[j])
            yU = upper_bound(vi[j])
            fval += c*Interval(xL, xU)*Interval(yL, yU)
        end
    end
    return fval.lo, fval.hi
end

function interval_bound(m::GlobalOptimizer, f::BufferedQuadraticEq, n::NodeBB)

    fval = Interval{Float64}(f.func.constant)
    for t in f.func.affine_terms
        c = t.coefficient
        j = t.variable_index.value
        xL = _lower_bound(FullVar(), m, j)
        xU = _upper_bound(FullVar(), m, j)
        fval += c*Interval(xL, xU)
    end
    for t in f.func.quadratic_terms
        c = t.coefficient
        i = t.variable_index_1.value
        j = t.variable_index_2.value
        xL = _lower_bound(FullVar(), m, i)
        xU = _upper_bound(FullVar(), m, i)
        if i == j
            fval += 0.5*c*pow(Interval(xL, xU), 2)
        else
            yL = _lower_bound(FullVar(), m, j)
            yU = _upper_bound(FullVar(), m, j)
            fval += c*Interval(xL, xU)*Interval(yL, yU)
        end
    end
    return fval.lo, fval.hi
end

###
### SECOND-ORDER CONE
###
function lower_interval_bound(m::GlobalOptimizer, d::BufferedSOC, n::NodeBB)

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
function lower_interval_bound(m::GlobalOptimizer, d::BufferedNonlinearFunction{V}, n::NodeBB) where V
    !_has_value(d) && forward_pass!(m._working_problem._relaxed_evaluator, d)
    _is_num(d) ? _num(d) : _interval(d).lo
end
function interval_bound(m::GlobalOptimizer, d::BufferedNonlinearFunction{V}, n::NodeBB) where V
    !_has_value(d) && forward_pass!(m._working_problem._relaxed_evaluator, d)
    v = _is_num(d) ? Interval{Float64}(_num(d)) : _interval(d)
    return v.lo, v.hi
end

function is_feasible(m::GlobalOptimizer, f::Union{AffineFunctionIneq,BufferedQuadraticIneq}, y::NodeBB)
    lower_interval_bound(m, f, y) <= 0.0
end
function is_feasible(m::GlobalOptimizer, f::Union{AffineFunctionEq,BufferedQuadraticEq}, y::NodeBB)
    l, u = lower_interval_bound(m, f, y)
    l <= 0.0 <= u
end
function is_feasible(m::GlobalOptimizer, f::BufferedNonlinearFunction{V}, y::NodeBB) where V
    l, u = lower_interval_bound(m, f, y)
    feasible_flag = (u < _lower_bound(f))
    feasible_flag && (l > _upper_bound(f))
end

bound_objective(m::GlobalOptimizer, f::BufferedNonlinearFunction, n::NodeBB) = interval_bound(m, f, n)
bound_objective(m::GlobalOptimizer, f::AffineFunctionIneq, n::NodeBB) = interval_bound(m, f)
bound_objective(m::GlobalOptimizer, f::BufferedQuadraticIneq, n::NodeBB) = interval_bound(m, f)
function bound_objective(m::GlobalOptimizer, f::SV, n::NodeBB)
    l = _lower_bound(FullVar(), m, f.variable.value)
    u = _lower_bound(FullVar(), m, f.variable.value)
    return l, u
end
function bound_objective(t::ExtensionType, m::GlobalOptimizer)
    bound_objective(m, m._working_problem._objective, m._current_node)
end
bound_objective(m::GlobalOptimizer{R,Q,S}) where {R,Q,S<:ExtensionType} = bound_objective(_ext_typ(m), m)
