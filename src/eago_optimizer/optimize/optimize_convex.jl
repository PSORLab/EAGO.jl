# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimize/optimize_convex.jl
# Contains the solve_local_nlp! routine which computes the optimal value
# of a convex function. This is used to compute the upper bound in the
# branch and bound routine. A number of utility function required for
# solve_local_nlp! are also included.
#############################################################################

"""

Shifts the resulting local nlp objective value `f*` by `(1.0 + relative_tolerance/100.0)*f* + absolute_tolerance/100.0`.
This assumes that the local solvers relative tolerance and absolute tolerance is significantly lower than the global
tolerance (local problem is minimum).
"""
function stored_adjusted_upper_bound!(d::Optimizer, v::Float64)
    adj_atol = d._parameters.absolute_tolerance/100.0
    adj_rtol = d._parameters.relative_tolerance/100.0
    if v > 0.0
        d._upper_objective_value = v*(1.0 + adj_rtol) + adj_atol
    else
        d._upper_objective_value = v*(1.0 - adj_rtol) + adj_atol
    end

    return nothing
end

revert_adjusted_upper_bound!(t::ExtensionType, d::Optimizer) = nothing

function revert_adjusted_upper_bound!(t::DefaultExt, d::Optimizer)

    adj_atol = d._parameters.absolute_tolerance/100.0
    adj_rtol = d._parameters.relative_tolerance/100.0

    adj_objective_value = d._global_upper_bound
    adj_objective_value -= adj_atol
    if adj_objective_value > 0.0
        adj_objective_value /= (1.0 + adj_rtol)
    else
        adj_objective_value /= (1.0 - adj_rtol)
    end
    d._global_upper_bound = adj_objective_value

    return nothing
end

# translates quadratic cone
function add_soc_constraints_as_quad!(m::Optimizer, opt::T) where T

    for (func, set) in m._input_problem._conic_second_order
        # quadratic cone implies variable[1] >= 0.0, bounds contracted accordingly in initial_parse!
        quad_terms = SQT[SQT((), func.variables[i], func.variables[i]) for i = 1:length(func.variables)]
        sqf = SQF(SQT[], SAF[], 0.0)
        MOI.add_constraint(opt, sqf, LT_ZERO)
    end

    return nothing
end

function _update_branch_variables!(d, m::Optimizer)
    for i = 1:_variable_num(FullVar(), m)
        upper_variable_index = m._upper_variables[i]
        single_variable = MOI.SingleVariable(upper_variable_index)
        if !_is_integer(FullVar(), m, i)
            vinfo = _variable_info(m,i)
            lvb  = _lower_bound(FullVar(), m, i)
            uvb  = _upper_bound(FullVar(), m, i)
            is_fixed(vinfo)        && MOI.add_constraint(d, single_variable, ET(lvb))
            is_less_than(vinfo)    && MOI.add_constraint(d, single_variable, LT(uvb))
            is_greater_than(vinfo) && MOI.add_constraint(d, single_variable, GT(lvb))
            if is_real_interval(vinfo)
                MOI.add_constraint(d, single_variable, LT(uvb))
                MOI.add_constraint(d, single_variable, GT(lvb))
            end
        end
    end
    return
end

function _finite_mid(l::T, u::T) where T
    (isfinite(l) && isfinite(u)) && return 0.5*(l + u)
    isfinite(l) ? l : (isfinite(u) ? u : zero(T))
end
function _set_starting_point!(d, m::Optimizer)
    for i = 1:_variable_num(FullVar(), m)
        l  = _lower_bound(FullVar(), m, i)
        u  = _upper_bound(FullVar(), m, i)
        v = m._upper_variables[i]
        MOI.set(d, MOI.VariablePrimalStart(), v, _finite_mid(l, u))
    end
    return
end

"""
$(SIGNATURES)

Takes an `MOI.TerminationStatusCode` and a `MOI.ResultStatusCode` and returns `true`
if this corresponds to a solution that is proven to be feasible.
Returns `false` otherwise.
"""
function is_feasible_solution(t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)

    termination_flag = false
    result_flag = false

    (t == MOI.OPTIMAL) && (termination_flag = true)
    (t == MOI.LOCALLY_SOLVED) && (termination_flag = true)

    # This is default solver specific... the acceptable constraint tolerances
    # are set to the same values as the basic tolerance. As a result, an
    # acceptably solved solution is feasible but non necessarily optimal
    # so it should be treated as a feasible point
    if (t == MOI.ALMOST_LOCALLY_SOLVED) && (r == MOI.NEARLY_FEASIBLE_POINT)
        termination_flag = true
        result_flag = true
    end
    (r == MOI.FEASIBLE_POINT) && (result_flag = true)

    return (termination_flag && result_flag)
end

function _unpack_local_nlp_solve!(m::Optimizer, opt::T) where T
    m._upper_termination_status = MOI.get(opt, MOI.TerminationStatus())
    m._upper_result_status = MOI.get(opt, MOI.PrimalStatus())
    if is_feasible_solution(m._upper_termination_status, m._upper_result_status)
        m._upper_feasibility = true
        obj_val = MOI.get(opt, MOI.ObjectiveValue())
        stored_adjusted_upper_bound!(m, obj_val)
        m._best_upper_value = min(obj_val, m._best_upper_value)
        m._upper_solution .= MOI.get(opt, MOI.VariablePrimal(), m._upper_variables)
    else
        m._upper_feasibility = false
        m._upper_objective_value = Inf
    end
    return
end

"""

Constructs and solves the problem locally on on node `y` updated the upper
solution informaton in the optimizer.
"""
function solve_local_nlp!(m::Optimizer)

    upper_optimizer = m.upper_optimizer
    MOI.empty!(upper_optimizer)

    for i = 1:m._working_problem._variable_count
        m._upper_variables[i] = MOI.add_variable(upper_optimizer)
    end
    _update_branch_variables!(upper_optimizer, m)
    _set_starting_point!(upper_optimizer, m)

    # Add linear and quadratic constraints to model
    add_linear_constraints!(m, upper_optimizer)

    for (func, set) in m._input_problem._quadratic_leq_constraints
        MOI.add_constraint(upper_optimizer, func, set)
    end
    for (func, set) in m._input_problem._quadratic_geq_constraints
        MOI.add_constraint(upper_optimizer, func, set)
    end
    for (func, set) in m._input_problem._quadratic_eq_constraints
        MOI.add_constraint(upper_optimizer, func, set)
    end

    if MOI.supports_constraint(upper_optimizer, VECOFVAR, SOC)
        add_soc_constraints!(m, upper_optimizer)
    else
        add_soc_constraints_as_quad!(m, upper_optimizer)
    end

    # Add nonlinear evaluation block
    MOI.set(upper_optimizer, MOI.NLPBlock(), m._working_problem._nlp_data)
    MOI.set(upper_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    # set objective as NECESSARY
    add_sv_or_aff_obj!(m, upper_optimizer)
    if m._input_problem._objective_type === SCALAR_QUADRATIC
        MOI.set(upper_optimizer, MOI.ObjectiveFunction{SQF}(), m._input_problem._objective_sqf)
    end

    # Optimizes the object
    MOI.optimize!(upper_optimizer)
    _unpack_local_nlp_solve!(m, upper_optimizer)
end

optimize!(::Val{DIFF_CVX}, m::Optimizer) = solve_local_nlp!(m)
