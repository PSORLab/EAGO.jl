function interval_lower_bound!(x::Optimizer, y::NodeBB, flag::Bool)

    feas = true

    d = x.working_evaluator_block.evaluator
    if x.objective === nothing
        objective_lo = get_node_lower(d.objective, 1)
    #elseif isa(x.objective)
    #    objective_lo
    #elseif isa(x.objective)
#        objective_lo
#    elseif isa(x.objective)
#        objective_lo
    end
    constraints_intv_lo = get_node_lower.(d.constraints, 1)
    constraints_intv_hi = get_node_upper.(d.constraints, 1)
    constraints_bnd_lo = d.constraints_lbd
    constrains_bnd_hi = d.constraints_ubd

    for i in 1:d.constraint_number
        if (constraints_bnd_lo[i] > constraints_intv_hi[i]) ||
           (constrains_bnd_hi[i] < constraints_intv_lo[i])
            feas = false
            break
        end
    end

    x.current_lower_info.feasibility = feas
    flag && (objective_lo = max(x.current_lower_info.value, objective_lo))
    x.current_lower_info.value = feas ? objective_lo : -Inf

end

"""
    default_lower_bounding!

Constructs and solves the relaxation using the default EAGO relaxation scheme
and optimizer on node `y`.
"""
function lower_problem!(x::Optimizer, y::NodeBB)

    # Copies initial model into working model (if initial model isn't dummy)
    # A dummy model is left iff all terms are relaxed
    if x.use_lower_factory
        factory = x.lower_factory(; x.lower_optimizer_options...)     # Should accept keyword arguments
        x.working_relaxed_optimizer = factory
        MOI.add_variables(x.working_relaxed_optimizer, x.variable_number)
    else
        if (x.initial_relaxed_optimizer != DummyOptimizer()) && (x.obbt_performed_flag)
            x.working_relaxed_optimizer = deepcopy(x.initial_relaxed_optimizer)           # TODO: Nix deepcopy here.
        end
    end

    xmid = 0.5*(lower_variable_bounds(y) + upper_variable_bounds(y))
    update_lower_variable_bounds1!(x, y, x.working_relaxed_optimizer)
    relax_function!(x, x.working_relaxed_optimizer, y, x.relaxation, xmid, load = true)
    relax_function!(x, x.working_relaxed_optimizer, y, x.relaxation, xmid, load = false)

    # Optimizes the object
    MOI.optimize!(x.working_relaxed_optimizer)

    # Process output info and save to CurrentUpperInfo object
    termination_status = MOI.get(x.working_relaxed_optimizer, MOI.TerminationStatus())
    result_status_code = MOI.get(x.working_relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(termination_status, result_status_code)
    solution = MOI.get(x.working_relaxed_optimizer, MOI.VariablePrimal(), x.lower_variables)

    if valid_flag
        if feasible_flag
            x.current_lower_info.feasibility = true
            x.current_lower_info.value = MOI.get(x.working_relaxed_optimizer, MOI.ObjectiveValue())
            vprimal_solution = MOI.get(x.working_relaxed_optimizer, MOI.VariablePrimal(), x.lower_variables)
            x.current_lower_info.solution[1:end] = vprimal_solution
            x._cut_add_flag = x.current_lower_info.feasibility
            set_dual!(x)
        else
            x._cut_add_flag = false
            x.current_lower_info.feasibility = false
            x.current_lower_info.value = -Inf
        end
    else
        interval_lower_bound!(x, y, true)
        x._cut_add_flag = false
    end
end
