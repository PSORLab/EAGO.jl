function triv_function(x) end

linear_solve!(m::Optimizer) = (println("Linear solve to be implemented. Recommend using linear solver such as GLPK directly. EAGO will continue..."))

function initialize_evaluators!(m::Optimizer, flag::Bool)

    num_nlp_constraints = length(m.nlp_data.constraint_bounds)

    # Build the JuMP NLP evaluator
    evaluator = m.nlp_data.evaluator
    features = MOI.features_available(evaluator)
    has_hessian = (:Hess in features)
    init_feat = [:Grad, :ExprGraph]
    num_nlp_constraints > 0 && push!(init_feat, :Jac)
    MOI.initialize(evaluator, init_feat)

    # Scrub user-defined functions
    if ~isa(evaluator, EAGO.EmptyNLPEvaluator)
        m.udf_scrubber_flag && Script.scrub!(evaluator.m.nlp_data)
        if m.udf_to_JuMP_flag
            Script.udf_loader!(m)
        end
    end
    #m.nlp_data.evaluator = evaluator #TODO: Rebuilt entire nlp_block...

    # Transform UDFs to JuMP ASTs
    #m.udf_to_JuMP_flag && Script.udf_loader!(m)
    # Rebuild the nlp-evaluator with udfs -> JuMP expressions

    #####
    #Script to unpack UDF from here
    #unpacked_evaluator = script_to_dag()
    #m.nlp.evaluator = unpacked_evaluator
    #####

    # Creates initial EAGO nlp evaluator for relaxations
    m.working_evaluator_block = m.nlp_data
    if ~isa(m.nlp_data.evaluator, EAGO.EmptyNLPEvaluator) || false #flag
        built_evaluator = build_nlp_evaluator(m.variable_number, m.nlp_data.evaluator, m, false)
        (m.optimization_sense == MOI.MAX_SENSE) && neg_objective!(built_evaluator)
        m.working_evaluator_block = MOI.NLPBlockData(m.nlp_data.constraint_bounds, built_evaluator, m.nlp_data.has_objective)
    end
    evaluator
end

function label_fixed_variables!(m::Optimizer)
    lbd = -Inf
    ubd = Inf
    for i in 1:m.variable_number
        lbd = m.variable_info[i].lower_bound
        ubd = m.variable_info[i].upper_bound
        if (lbd == ubd)
            m.variable_info[i].is_fixed = true
            m.fixed_variable[i] = true
        end
    end
end

function label_obbt_variables!(m::Optimizer)

    if ~in(true, values(m.bisection_variable))
        linear_solve!(m)
    end

    for i=1:length(m.variable_info)
        if m.bisection_variable[i]
            push!(m.obbt_variables, MOI.VariableIndex(i))
        end
    end

    no_constraints = true
    if length(m.nlp_data.constraint_bounds) > 0
        no_constraints &= false
    end
    if length(m.quadratic_leq_constraints) > 0
        no_constraints &= false
    end
    if length(m.quadratic_geq_constraints) > 0
        no_constraints &= false
    end
    if length(m.quadratic_eq_constraints) > 0
        no_constraints &= false
    end
    if length(m.quadratic_interval_constraints) > 0
        no_constraints &= false
    end
    if no_constraints
        empty!(m.obbt_variables)
        m.obbt_depth = 0
        m.cp_depth = 0
    end

end

function bld_user_upper_fact!(m::Optimizer)
    MOI.add_variables(m.initial_upper_optimizer, m.variable_number)
    m.upper_variables = MOI.add_variables(m.working_upper_optimizer, m.variable_number)
    set_local_nlp!(m)
end

function user_reformed_optimizer(m::Optimizer)
    deepcopy(backend(m.nlp_data.evaluator.m).optimizer.model.optimizer)
end

function local_solve!(x::Optimizer)

  key_pair::Pair{Int,NodeBB} = node_selection(x)
  current_node::NodeBB = kn_pair[2]
  upper_time = @elapsed solve_local_nlp!(x, current_node)

  if x.current_upper_info.feasibility
    x.feasible_solution_found = true
    x.first_solution_node = x.maximum_node_id
    x.solution_value = x.current_upper_info.value
    x.continuous_solution[:] = x.current_upper_info.solution
    x.termination_status_code = MOI.LOCALLY_SOLVED
    x.result_status_code = MOI.FEASIBLE_POINT
  else
    x.feasible_solution_found = false
    x.first_solution_node = x.maximum_node_id
    x.termination_status_code = MOI.LOCALLY_INFEASIBLE
    x.result_status_code = MOI.INFEASIBLE_POINT
  end
  x.history.upper_time[1] = upper_time

end

function MOI.optimize!(m::Optimizer)
#function MOI.optimize!(m::Optimizer; ignore_optimize_hook = (model.optimize_hook === nothing))


    setrounding(Interval, m.rounding_mode)

    ########### Reformulate DAG using auxilliary variables ###########
    _variable_len = length(m.variable_info)
    m.continuous_variable_number = _variable_len
    m.variable_number = _variable_len

    ########### Set Correct Size for Problem Storage #########
    m.current_lower_info.solution = fill(0.0, _variable_len)
    m.current_lower_info.lower_variable_dual = fill(0.0, _variable_len)
    m.current_lower_info.upper_variable_dual = fill(0.0, _variable_len)
    m.current_upper_info.solution = fill(0.0, _variable_len)

    # loads variables into the working model
    m.variable_index_to_storage = Dict{Int,Int}()
    for i=1:_variable_len
        m.variable_index_to_storage[i] = i
    end
    m.storage_index_to_variable = ReverseDict(m.variable_index_to_storage)

    # Get various other sizes
    m.continuous_solution = zeros(Float64, _variable_len)

    set_to_default!(m)                             # Sets any unset functions to default values
    initialize_evaluators!(m, false)               # initializes the EAGO and JuMP NLP evaluators


    m.reform_epigraph_flag && reform_epigraph!(m)  # perform epigraph rearrangement
    #m.reform_cse_flag && dag_cse_simplify!(m)      #
    #m.reform_flatten_flag && dag_flattening!(m)

    #m = user_reformed_optimizer(m)
    #m.debug1 = initialize_evaluators!(m, true)                      # re-initializes evaluators after reformulations are performed

    create_initial_node!(m)                        # Create initial node and add it to the stack

    m.lower_variables = MOI.VariableIndex[MOI.VariableIndex(i) for i in 1:_variable_len]
    m.upper_variables = MOI.add_variables(m.initial_relaxed_optimizer, _variable_len)

    ###### OBBT Setup #####
    label_fixed_variables!(m)                                   # label variable fixed to a value
    _nlpdata = m.nlp_data
    _evaluator = _nlpdata.evaluator::MOI.AbstractNLPEvaluator
    label_nonlinear_variables!(m, _evaluator)
    label_obbt_variables!(m)

    # Relax initial model terms
    xmid = 0.5*(lower_variable_bounds(m.stack[1]) + upper_variable_bounds(m.stack[1]))
    relax_function!(m, m.initial_relaxed_optimizer, m.stack[1], m.relaxation, xmid, load = true)

    # Allow for a hook to modify branch and bound routine
    ignore_optimize_hook = (m.optimize_hook === nothing)
    if !ignore_optimize_hook
        return m.optimize_hook(m; kwargs...)
    end

    (~m.use_upper_factory) && bld_user_upper_fact!(m) # if optimizer type is supplied for upper, build factory

    # Runs the branch and bound routine
    if m.local_solve_only
        local_solve!(m)
    else
        global_solve!(m)
    end
end
