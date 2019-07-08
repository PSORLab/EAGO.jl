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
    MOI.initialize(evaluator,init_feat)

    # Creates initial EAGO nlp evaluator for relaxations
    m.working_evaluator_block = m.nlp_data
    if ~isa(m.nlp_data.evaluator, EAGO.EmptyNLPEvaluator) || false #flag
        built_evaluator = build_nlp_evaluator(MC{m.variable_number}, m.nlp_data.evaluator, m, false)
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

    if ~in(true, values(m.nonlinear_variable))
        linear_solve!(m)
    end

    for i=1:length(m.variable_info)
        if m.nonlinear_variable[i]
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

function MOI.optimize!(m::Optimizer; custom_mod! = triv_function, custom_mod_args = (1,))


    setrounding(Interval, m.rounding_mode)

    ########### Reformulate DAG using auxilliary variables ###########
    NewVariableSize = length(m.variable_info)
    m.continuous_variable_number = NewVariableSize
    m.variable_number = NewVariableSize

    ########### Set Correct Size for Problem Storage #########
    m.current_lower_info.solution = Float64[0.0 for i=1:NewVariableSize]
    m.current_lower_info.lower_variable_dual = Float64[0.0 for i=1:NewVariableSize]
    m.current_lower_info.upper_variable_dual = Float64[0.0 for i=1:NewVariableSize]
    m.current_upper_info.solution = Float64[0.0 for i=1:NewVariableSize]

    # loads variables into the working model
    m.variable_index_to_storage = Dict{Int,Int}()
    for i=1:NewVariableSize
        m.variable_index_to_storage[i] = i
    end
    m.storage_index_to_variable = ReverseDict(m.variable_index_to_storage)

    # Get various other sizes
    m.continuous_solution = zeros(Float64,NewVariableSize)

    set_to_default!(m)                             # Sets any unset functions to default values
    initialize_evaluators!(m, false)               # initializes the EAGO and JuMP NLP evaluators

    m.reform_epigraph_flag && reform_epigraph!(m)  # perform epigraph rearrangement
    #m.reform_cse_flag && dag_cse_simplify!(m)      #
    #m.reform_flatten_flag && dag_flattening!(m)

    #m = user_reformed_optimizer(m)
    #m.debug1 = initialize_evaluators!(m, true)                      # re-initializes evaluators after reformulations are performed

    create_initial_node!(m)                        # Create initial node and add it to the stack

    m.upper_variables = MOI.add_variables(m.initial_relaxed_optimizer, m.variable_number)
    m.lower_variables = MOI.VariableIndex.(1:m.variable_number)

    ###### OBBT Setup #####
    label_fixed_variables!(m)                                   # label variable fixed to a value
    label_nonlinear_variables!(m, m.nlp_data.evaluator)
    label_obbt_variables!(m)

    # Relax initial model terms
    xmid = 0.5*(lower_variable_bounds(m.stack[1]) + upper_variable_bounds(m.stack[1]))
    relax_model!(m, m.initial_relaxed_optimizer, m.stack[1], m.relaxation, xmid, load = true)

    # Runs a customized function if one is provided
    m.custom_mod_flag = (custom_mod! != triv_function)
    if m.custom_mod_flag
        custom_mod!(m, custom_mod_args)
    end

    (~m.use_upper_factory) && bld_user_upper_fact!(m) # if optimizer type is supplied for upper, build factory

    # Runs the branch and bound routine
    solve_nlp!(m)
end
