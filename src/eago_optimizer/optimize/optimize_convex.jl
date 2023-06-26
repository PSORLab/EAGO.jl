# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/optimize/optimize_convex.jl
# Contains the solve_local_nlp! routine which computes the optimal value
# of a convex function. This is used to compute the upper bound in the
# branch and bound routine. A number of utility functions required for
# solve_local_nlp! are also included.
################################################################################

"""
$(TYPEDSIGNATURES)

Checks that the solution of a local solve is integer feasible to within
the tolerances specified by `integer_abs_tol` and `integer_rel_tol`.
"""
function is_integer_feasible_local(m::GlobalOptimizer, d)
    bool = true
    atol = _integer_abs_tol(m)
    rtol = _integer_rel_tol(m)
    for i = 1:_variable_num(BranchVar(), m)
        if is_integer(BranchVar(), m, i)
            xsol = MOI.get(d, MOI.VariablePrimal(), m._upper_variables[i])
            if isapprox(floor(xsol), xsol; atol = atol, rtol = rtol)
                continue
            elseif isapprox(ceil(xsol), xsol; atol = atol, rtol = rtol)
                continue
            else
                bool &= false
                break
            end
        end
    end
    return bool
end

"""
$(SIGNATURES)

Shifts the resulting local nlp objective value `f*` by `(1.0 + relative_tolerance/100.0)*f* + absolute_tolerance/100.0`.
This assumes that the local solvers relative tolerance and absolute tolerance is significantly lower than the global
tolerance (local problem is minimum).
"""
function stored_adjusted_upper_bound!(d::GlobalOptimizer, v::Float64)
    adj_atol = d._parameters.absolute_tolerance/100.0
    adj_rtol = d._parameters.relative_tolerance/100.0
    if v > 0.0
        d._upper_objective_value = v*(1.0 + adj_rtol) + adj_atol
    else
        d._upper_objective_value = v*(1.0 - adj_rtol) + adj_atol
    end
    return nothing
end

function _update_upper_variables!(d, m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    for i = 1:_variable_num(FullVar(), m)
        v = m._upper_variables[i]
        l  = _lower_bound(FullVar(), m, i)
        u  = _upper_bound(FullVar(), m, i)
        if is_integer(FullVar(), m, i)
            l = ceil(l)
            u = floor(u)
        end
        is_fixed_int = l == u
        vi = _working_variable_info(m,i)
        if is_fixed(vi) || is_fixed_int
            MOI.add_constraint(d, v, ET(l))
        elseif is_less_than(vi)
            MOI.add_constraint(d, v, LT(u))
        elseif is_greater_than(vi)
            MOI.add_constraint(d, v, GT(l))
        elseif is_real_interval(vi)
            MOI.add_constraint(d, v, LT(u))
            MOI.add_constraint(d, v, GT(l))
        end
    end
    return
end

function _finite_mid(l::T, u::T) where T
    (isfinite(l) && isfinite(u)) && return 0.5*(l + u)
    isfinite(l) ? l : (isfinite(u) ? u : zero(T))
end
function _set_starting_point!(d, m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    for i = 1:_variable_num(FullVar(), m)
        l  = _lower_bound(FullVar(), m, i)
        u  = _upper_bound(FullVar(), m, i)
        v = m._upper_variables[i]
        MOI.set(d, MOI.VariablePrimalStart(), v, _finite_mid(l, u))
    end
    return
end

"""
    LocalResultStatus

Status code used internally to determine how to interpret the results from the
solution of a local problem solve.
"""
@enum(LocalResultStatus, LRS_FEASIBLE, LRS_OTHER)

"""
$(SIGNATURES)

Takes an `MOI.TerminationStatusCode` and a `MOI.ResultStatusCode` and returns `true`
if this corresponds to a solution that is proven to be feasible.
Returns `false` otherwise.
"""
function local_problem_status(t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)

    if (t == MOI.OPTIMAL) && (r == MOI.FEASIBLE_POINT)
        return LRS_FEASIBLE
    elseif (t == MOI.LOCALLY_SOLVED) && (r == MOI.FEASIBLE_POINT)
        return LRS_FEASIBLE
    # This is default solver specific. The acceptable constraint tolerances
    # are set to the same values as the basic tolerance. As a result, an
    # acceptably solved solution is feasible but non necessarily optimal
    # so it should be treated as a feasible point.
    elseif (t == MOI.ALMOST_LOCALLY_SOLVED) && (r == MOI.NEARLY_FEASIBLE_POINT)
        return LRS_FEASIBLE
    end
    return LRS_OTHER
end

function _unpack_local_nlp_solve!(m::GlobalOptimizer, d::T) where T

    tstatus = MOI.get(d, MOI.TerminationStatus())
    pstatus = MOI.get(d, MOI.PrimalStatus())
    m._upper_termination_status = tstatus
    m._upper_result_status = pstatus

    if local_problem_status(tstatus, pstatus) == LRS_FEASIBLE
        if is_integer_feasible_local(m, d)

            m._upper_feasibility = true
            obj_val = MOI.get(d, MOI.ObjectiveValue())
            stored_adjusted_upper_bound!(m, obj_val)
            m._best_upper_value = min(obj_val, m._best_upper_value)
            m._upper_solution .= MOI.get(d, MOI.VariablePrimal(), m._upper_variables)
            
            ip = m._input_problem
            _extract_primal_linear!(d, ip)
            _extract_primal_quadratic!(d, ip)
        end
    else
        m._upper_feasibility = false
        m._upper_objective_value = Inf
    end
    return
end

"""

Constructs and solves the problem locally on node `y` updated the upper
solution informaton in the optimizer.
"""
function solve_local_nlp!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}

    upper_optimizer = _upper_optimizer(m)
    MOI.empty!(upper_optimizer)

    for i = 1:m._working_problem._variable_count
        m._upper_variables[i] = MOI.add_variable(upper_optimizer)
    end
    _update_upper_variables!(upper_optimizer, m)
    _set_starting_point!(upper_optimizer, m)

    # Add constraints
    ip = m._input_problem
    _add_constraint_store_ci_linear!(upper_optimizer, ip)
    _add_constraint_store_ci_quadratic!(upper_optimizer, ip)
     #add_soc_constraints!(m, upper_optimizer)
  
    # Add nonlinear evaluation block
    MOI.set(upper_optimizer, MOI.NLPBlock(), m._working_problem._nlp_data)
    MOI.set(upper_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(upper_optimizer, MOI.ObjectiveFunction{SAF}(), m._working_problem._objective_saf)

    # Optimizes the object
    MOI.optimize!(upper_optimizer)
    _unpack_local_nlp_solve!(m, upper_optimizer)
end

optimize!(::DIFF_CVX, m::Optimizer) = solve_local_nlp!(m._global_optimizer)
optimize!(::DIFF_CVX, m::GlobalOptimizer) = solve_local_nlp!(m)