# is root node? is at iteration number? did last bound improve? did last bound land in current domain? can omit is_integer_feasible(x) since still an NLP solver
function default_nlp_heurestic(x::Optimizer, y::NodeBB)
    bool = false
    bool |= (y.depth <= x.upper_bounding_depth)
    bool |= (~y.branch_direction & (rand() < 0.5^(y.depth-x.upper_bounding_depth)))   # last attempt improved and no solution has been found in y
    return bool
end

"""
    default_upper_bounding!

Constructs and solves the problem locally on on node `y` and saves upper
bounding info to `x.current_upper_info`.
"""
function default_upper_bounding!(x::Optimizer,y::NodeBB)
    if default_nlp_heurestic(x,y)
        if x.use_upper_factory
            factory = x.upper_factory()
            x.initial_upper_optimizer = factory
            x.upper_variables = MOI.add_variables(x.initial_upper_optimizer, x.variable_number)
            set_local_nlp!(x)
            x.working_upper_optimizer = x.initial_upper_optimizer
        else
            if x.initial_upper_optimizer != DummyOptimizer()
                x.working_upper_optimizer = deepcopy(x.initial_upper_optimizer)
            end
        end
        update_upper_variable_bounds!(x,y,x.working_upper_optimizer)

        if x.upper_has_node
            set_current_node!(x.working_upper_optimizer.nlp_data.evaluator,y)
        end

        # Optimizes the object
        MOI.optimize!(x.working_upper_optimizer)

        # Process output info and save to CurrentUpperInfo object
        termination_status = MOI.get(x.working_upper_optimizer, MOI.TerminationStatus())
        result_status = MOI.get(x.working_upper_optimizer, MOI.PrimalStatus())
        solution = MOI.get(x.working_upper_optimizer, MOI.VariablePrimal(), x.upper_variables)

        if is_feasible_solution(termination_status, result_status)
            x.current_upper_info.feasibility = true
            mult = (x.optimization_sense == MOI.MAX_SENSE && x.objective == nothing) ? -1.0 : 1.0
            x.current_upper_info.value = mult*MOI.get(x.working_upper_optimizer, MOI.ObjectiveValue())
            x.current_upper_info.solution[1:end] = MOI.get(x.working_upper_optimizer, MOI.VariablePrimal(), x.upper_variables)
        else
            x.current_upper_info.feasibility = false
            x.current_upper_info.value = Inf
        end
    else
        x.current_upper_info.feasibility = false
        x.current_upper_info.value = Inf
    end
end
