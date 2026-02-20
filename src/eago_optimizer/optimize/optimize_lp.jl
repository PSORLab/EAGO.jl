# Copyright (c) 2018 Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Alireza Miraliakbar, Matthew Stuber, and the University of Connecticut (UConn)
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/optimize/optimize_lp.jl
# Contains the optimize! routines for LP and MILP type problems. This also
# includes functions to add variables, linear constraints, soc constraints, 
# and unpack solutions.
################################################################################

function add_variables(m::GlobalOptimizer, d)
    n = m._input_problem._variable_count
    z = fill(VI(1), n)
    for i = 1:n
        z[i] = MOI.add_variable(d)
        vi = m._working_problem._variable_info[i]
        if is_fixed(vi)
            MOI.add_constraint(d, z[i], ET(vi))
        elseif is_interval(vi)
            MOI.add_constraint(d, z[i], IT(vi))
        elseif is_greater_than(vi)
            MOI.add_constraint(d, z[i], GT(vi))
        elseif is_less_than(vi)
            MOI.add_constraint(d, z[i], LT(vi))
        end
        if is_integer(vi)
            MOI.add_constraint(d, z[i], INT())
        end
    end
    return z
end

lp_obj!(m::GlobalOptimizer, d, f::Nothing) = false
function lp_obj!(m::GlobalOptimizer, d, f::VI)
    MOI.set(d, MOI.ObjectiveFunction{VI}(), f)
    MOI.set(d, MOI.ObjectiveSense(), m._input_problem._optimization_sense)
    return false
end
function lp_obj!(m::GlobalOptimizer, d, f::SAF)
    MOI.set(d, MOI.ObjectiveFunction{SAF}(), f)
    MOI.set(d, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    return m._input_problem._optimization_sense == MOI.MAX_SENSE
end

function optimize!(::LP, m::Optimizer{Q,S,T}) where {Q,S,T}

    d = m._global_optimizer
    ip = d._input_problem
    r = _relaxed_optimizer(d)
    MOI.empty!(r)

    d._relaxed_variable_index = add_variables(d, r)
    _add_constraint_store_ci_linear!(r, ip)
    
    min_to_max = lp_obj!(d, r, ip._objective)
    if ip._optimization_sense == MOI.FEASIBILITY_SENSE
        MOI.set(r, MOI.ObjectiveSense(), ip._optimization_sense)
    end

    (d._parameters.verbosity < 5) && MOI.set(r, MOI.Silent(), true)
    d._parse_time = time() - d._start_time

    MOI.optimize!(r)

    m._termination_status_code = MOI.get(r, MOI.TerminationStatus())
    m._result_status_code = MOI.get(r, MOI.PrimalStatus())

    if MOI.get(r, MOI.ResultCount()) > 0

        obj_val = MOI.get(r, MOI.ObjectiveValue())
        if min_to_max
            obj_val *= -1.0
        end
        d._global_lower_bound = obj_val
        d._global_upper_bound = obj_val
        d._best_upper_value = obj_val
        d._solution_value = obj_val
        m._objective_value = obj_val
        m._objective_bound = obj_val

        d._continuous_solution = zeros(d._input_problem._variable_count)
        for i = 1:d._input_problem._variable_count
            d._continuous_solution[i] = MOI.get(r, MOI.VariablePrimal(),  d._relaxed_variable_index[i])
        end

        _extract_primal_linear!(r, ip)
    end
    d._run_time = time() - d._start_time
    return
end

optimize!(::MILP, m::Optimizer{Q,S,T}) where {Q,S,T} = optimize!(LP(), m)