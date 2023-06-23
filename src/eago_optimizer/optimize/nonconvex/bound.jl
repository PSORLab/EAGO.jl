# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/optimize/nonconvex/bound.jl
# Computes interval bounds of various functions.
################################################################################

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

function interval_bound(m::GlobalOptimizer, f::Union{AffineFunctionEq,AffineFunctionIneq})
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

function lower_interval_bound(m::GlobalOptimizer, f::BufferedQuadraticIneq)
    fval = Interval{Float64}(f.func.constant)
    for t in f.func.affine_terms
        c = t.coefficient
        j = t.variable.value
        xL = _lower_bound(FullVar(), m, j)
        xU = _upper_bound(FullVar(), m, j)
        fval += c > 0.0 ? c*xL : c*xU
    end
    for t in f.func.quadratic_terms
        c = t.coefficient
        i = t.variable_1.value
        j = t.variable_2.value
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
function interval_bound(m::GlobalOptimizer, f::BufferedQuadraticIneq)

    fval = Interval{Float64}(f.func.constant)
    for t in f.func.affine_terms
        c = t.coefficient
        j = t.variable.value
        xL = _lower_bound(FullVar(), m, j)
        xU = _upper_bound(FullVar(), m, j)
        fval += c*Interval(xL, xU)
    end
    for t in f.func.quadratic_terms
        c = t.coefficient
        i = t.variable_1.value
        j = t.variable_2.value
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

function interval_bound(m::GlobalOptimizer, f::BufferedQuadraticEq)

    fval = Interval{Float64}(f.func.constant)
    for t in f.func.affine_terms
        c = t.coefficient
        j = t.variable.value
        xL = _lower_bound(FullVar(), m, j)
        xU = _upper_bound(FullVar(), m, j)
        fval += c*Interval(xL, xU)
    end
    for t in f.func.quadratic_terms
        c = t.coefficient
        i = t.variable_1.value
        j = t.variable_2.value
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
function lower_interval_bound(m::GlobalOptimizer, d::BufferedSOC)

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
function lower_interval_bound(m::GlobalOptimizer, d::BufferedNonlinearFunction{V,N,T}) where {V,N,T}
    !has_value(d) && forward_pass!(m._working_problem._relaxed_evaluator, d)
    is_num(d) ? num(d) : interval(d).lo
end
function interval_bound(m::GlobalOptimizer, d::BufferedNonlinearFunction{V,N,T}) where {V,N,T}
    !has_value(d) && forward_pass!(m._working_problem._relaxed_evaluator, d)
    v = is_num(d) ? Interval{Float64}(num(d)) : interval(d)
    return v.lo, v.hi
end

"""
    is_feasible(::GlobalOptimizer, ::T)

Check if a given bound is feasible.

# Options for T (all are subtypes of AbstractEAGOConstraint):
- AffineFunctionIneq
- AffineFunctionEq
- BufferedQuadraticIneq
- BufferedQuadraticEq
- BufferedNonlinearFunction{V,N,T}
"""
is_feasible(m::GlobalOptimizer, f::Union{AFI,BQI}) = lower_interval_bound(m, f) <= 0.0
function is_feasible(m::GlobalOptimizer, f::Union{AFE,BQE})
    l, u = interval_bound(m, f)
    l <= 0.0 <= u
end
function is_feasible(m::GlobalOptimizer, f::BufferedNonlinearFunction{V,N,T}) where {V,N,T}
    l, u = interval_bound(m, f)
    feasible_flag = (u >= lower_bound(f))
    feasible_flag && (l <= upper_bound(f))
end

"""
    bound_objective(::GlobalOptimizer, ::T)
    bound_objective(::ExtensionType, ::GlobalOptimizer)

Compute a tuple representing the lower and upper bounds for an objective function.
Note: `bound_objective(::GlobalOptimizer)` dispatches to 
`bound_objective(::ExtensionType, ::GlobalOptimizer)`.

# Options for T:
- AffineFunctionIneq
- BufferedNonlinearFunction
- BufferedQuadraticIneq
- MOI.VariableIndex
"""
bound_objective(m::GlobalOptimizer, f::BufferedNonlinearFunction) = interval_bound(m, f)
bound_objective(m::GlobalOptimizer, f::AffineFunctionIneq) = interval_bound(m, f)
bound_objective(m::GlobalOptimizer, f::BufferedQuadraticIneq) = interval_bound(m, f)
function bound_objective(m::GlobalOptimizer, f::VI)
    vval = f.value
    l = lower_bound(FullVar(), m, vval)
    u = upper_bound(FullVar(), m, vval)
    return l, u
end

bound_objective(t::ExtensionType, m::GlobalOptimizer) = bound_objective(m, m._working_problem._objective)
bound_objective(m::GlobalOptimizer{R,Q,S}) where {R,Q,S<:ExtensionType} = bound_objective(_ext(m), m)
