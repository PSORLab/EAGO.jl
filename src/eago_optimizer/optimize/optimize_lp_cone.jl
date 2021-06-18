# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eeago_optimizer/optimize/optimize_lp_cone.jl
# Contains the optimize! routines for LP, SOCP, (and in the future MILP and
# MISOCP) type problems. This also includes functions to add variables,
# linear constraints, soc constraints, and unpack solutions.
#############################################################################

function add_variables(m::GlobalOptimizer, optimizer, variable_number::Int)

    variable_index = fill(VI(1), variable_number)
    for i = 1:variable_number
        @inbounds variable_index[i] = MOI.add_variable(optimizer)
        relaxed_variable = SV(@inbounds variable_index[i])
        v_info = @inbounds m._working_problem._variable_info[i]
        if v_info.is_integer && v_info.is_fixed
            MOI.add_constraint(optimizer, relaxed_variable, ET(v_info.lower_bound))
        elseif v_info.is_integer
            MOI.add_constraint(optimizer, relaxed_variable, ZO())
        elseif v_info.is_fixed
            MOI.add_constraint(optimizer, relaxed_variable, ET(v_info.lower_bound))
        elseif v_info.has_lower_bound && v_info.has_upper_bound
            MOI.add_constraint(optimizer, relaxed_variable, IT(v_info.lower_bound, v_info.upper_bound))
        elseif v_info.has_lower_bound
            MOI.add_constraint(optimizer, relaxed_variable, GT(v_info.lower_bound))
        elseif v_info.has_upper_bound
            MOI.add_constraint(optimizer, relaxed_variable, LT(v_info.upper_bound))
        end
    end

    return variable_index
end

### LP and MILP routines
function add_linear_constraints!(m::GlobalOptimizer, d::T) where T

    # add linear constraints
    for (f, leq) in m._input_problem._linear_leq_constraints
        MOI.add_constraint(d, f, leq)
    end
    for (f, geq) in m._input_problem._linear_geq_constraints
        MOI.add_constraint(d, f, geq)
    end
    for (f, eq) in m._input_problem._linear_eq_constraints
        MOI.add_constraint(d, f, eq)
    end

    return nothing
end

### LP and MILP routines
function add_soc_constraints!(m::GlobalOptimizer, opt::T) where T

    for (func, set) in m._input_problem._conic_second_order
         MOI.add_constraint(opt, func, set)
    end

    return nothing
end

add_sv_or_aff_obj!(m::GlobalOptimizer, d::T, f::Nothing) where T = nothing
function add_sv_or_aff_obj!(m::GlobalOptimizer, d::T, f::F) where {T,F<:Union{SV,SAF}}
    MOI.set(d, MOI.ObjectiveFunction{F}(), f)
    return
end

function unpack_local_solve!(m::GlobalOptimizer, opt::T) where T

    m._maximum_node_id = 0

    m._termination_status_code = MOI.get(opt, MOI.TerminationStatus())
    m._result_status_code = MOI.get(opt, MOI.PrimalStatus())
    m._feasible_solution_found = m._result_status_code === MOI.FEASIBLE_POINT

    if MOI.get(opt, MOI.ResultCount()) > 0

        objective_value = MOI.get(opt, MOI.ObjectiveValue())

        # corrects for standard printing multiplier
        if !_is_input_min(m)
            objective_value *= -1.0
        end

        m._global_lower_bound = objective_value
        m._global_upper_bound = objective_value
        m._best_upper_value = objective_value
        m._solution_value = objective_value
    end

    m._continuous_solution = zeros(m._input_problem._variable_count)
    for i = 1:m._input_problem._variable_count
         m._continuous_solution[i] = MOI.get(opt, MOI.VariablePrimal(),  m._relaxed_variable_index[i])
    end

    m._run_time = time() - m._start_time

    return nothing
end

function optimize!(::LP, m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}

    relaxed_optimizer = m.relaxed_optimizer
    MOI.empty!(relaxed_optimizer)

    m._relaxed_variable_index = add_variables(m, relaxed_optimizer, m._input_problem._variable_count)
    add_linear_constraints!(m, relaxed_optimizer)
    add_sv_or_aff_obj!(m, relaxed_optimizer,  m._input_problem._objective)
    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), m._input_problem._optimization_sense)

    if m._parameters.verbosity < 5
        MOI.set(relaxed_optimizer, MOI.Silent(), true)
    end
    m._parse_time = time() - m._start_time

    MOI.optimize!(relaxed_optimizer)

    unpack_local_solve!(m, relaxed_optimizer)

    return nothing
end

optimize!(::MILP, m::GlobalOptimizer) = optimize!(LP(), m)

function optimize!(::SOCP, m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}

    relaxed_optimizer = m.relaxed_optimizer
    MOI.empty!(relaxed_optimizer)

    m._relaxed_variable_index = add_variables(m, relaxed_optimizer, m._input_problem._variable_count)
    add_linear_constraints!(m, relaxed_optimizer)
    add_soc_constraints!(m, relaxed_optimizer)
    add_sv_or_aff_obj!(m, relaxed_optimizer,  m._input_problem._objective)
    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), m._input_problem._optimization_sense)

    if m._parameters.verbosity < 5
        MOI.set(relaxed_optimizer, MOI.Silent(), true)
    end
    m._parse_time = time() - m._start_time

    MOI.optimize!(relaxed_optimizer)

    unpack_local_solve!(m, relaxed_optimizer)

    return
end

optimize!(::MISOCP, m::GlobalOptimizer) = optimize!(SOCP(), m)
