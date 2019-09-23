"""
    gen_quad_sparsity

Generate the list of variables (1,...,n) that are included in the func.
"""
function gen_quad_sparsity(func::SQF)
    list = Int64[]
    for term in func.affine_terms
        val = term.variable_index.value
        push!(list, val)
    end
    for term in func.quadratic_terms
        val = term.variable_index_1.value
        valn = term.variable_index_1.value
        push!(list, val, valn)
    end
    sort!(list); unique!(list)
    return list
end

"""
    gen_quadratic_storage!

Generate the storage for the affine relaxations of quadratic functions.
"""
function gen_quadratic_storage!(x::Optimizer)

    opt = x.relaxed_optimizer
    variable_index = x._lower_variable_index

    for (func, set, ind) in x._quadratic_leq_constraints

        v = gen_quad_sparsity(func)
        nv = length(v)
        push!(x._quadratic_leq_sparsity, v)
        push!(x._quadratic_leq_gradnz, nv)

        d = ImmutableDict{Int64,Int64}()
        for (i, val) in enumerate(v)
            ImmutableDict(d, val => i)
        end
        push!(x._quadratic_leq_dict, d)

        @inbounds nzvar = variable_index[v]
        func = SAF(SAT.(zeros(nv), nzvar), 0.0)
        ci = MOI.add_constraint(opt, func, LT(0.0))
        push!(x._quadratic_ci_leq, ci)

    end

    for (func, set, ind) in x._quadratic_geq_constraints

        v = gen_quad_sparsity(func)
        nv = length(v)
        push!(x._quadratic_geq_sparsity, v)
        push!(x._quadratic_geq_gradnz, nv)

        d = ImmutableDict{Int64,Int64}()
        for (i, val) in enumerate(v)
            ImmutableDict(d, val => i)
        end
        push!(x._quadratic_geq_dict, d)

        @inbounds nzvar = variable_index[v]
        func = SAF(SAT.(zeros(nv), nzvar), 0.0)
        ci = MOI.add_constraint(opt, func, LT(0.0))
        push!(x._quadratic_ci_geq, ci)

    end

    for (func, set, ind) in x._quadratic_eq_constraints

        v = gen_quad_sparsity(func)
        nv = length(v)
        push!(x._quadratic_eq_sparsity, v)
        push!(x._quadratic_eq_gradnz, nv)

        d = ImmutableDict{Int64,Int64}()
        for (i, val) in enumerate(v)
            ImmutableDict(d, val => i)
        end
        push!(x._quadratic_eq_dict, d)

        @inbounds nzvar = variable_index[v]
        func = SAF(SAT.(zeros(nv), nzvar), 0.0)
        c1 = MOI.add_constraint(opt, func, LT(0.0))
        c2 = MOI.add_constraint(opt, func, LT(0.0))
        push!(x._quadratic_ci_eq, c1, c2)

    end
    return
end

"""
    load_relaxed_problem!

Loads variables, linear constraints, and empty storage for first nlp and
quadratic cut.
"""
function load_relaxed_problem!(x::Optimizer)
    opt = x.relaxed_optimizer

    # add variables and indices
    variable_index = MOI.add_variables(opt, x._variable_number)
    variable = SV.(variable_index)
    append!(x._lower_variable_index, variable_index)
    append!(x._lower_variable, variable)

    for i in 1:x._variable_number

        @inbounds var = x._variable_info[i]
        @inbounds variable_i = variable[i]

        if var.is_integer
        else
            if var.is_fixed
                @inbounds bnd = var.lower_bound
                ci1 = MOI.add_constraint(opt, variable_i, ET(bnd))
                push!(x._lower_variable_et, ci1)
                push!(x._lower_variable_style, 1)

            elseif var.has_lower_bound && var.has_upper_bound
                @inbounds bnd = var.upper_bound
                ci2 = MOI.add_constraint(opt, variable_i, LT(bnd))
                push!(x._lower_variable_lt, ci2)
                @inbounds bnd = var.lower_bound
                ci3 = MOI.add_constraint(opt, variable_i, GT(bnd))
                push!(x._lower_variable_gt, ci3)
                push!(x._lower_variable_style, 2)

            elseif var.has_lower_bound
                @inbounds bnd = var.lower_bound
                ci4 = MOI.add_constraint(opt, variable_i, GT(bnd))
                push!(x._lower_variable_gt, ci4)
                push!(x._lower_variable_style, 3)

            elseif var.has_upper_bound
                @inbounds bnd = var.upper_bound
                ci5 = MOI.add_constraint(opt, var_xi, LT(bnd))
                push!(x._lower_variable_lt, ci5)
                push!(x._lower_variable_style, 4)
            end
        end
    end

    # add linear constraints
    for (func, set) in x._linear_leq_constraints
         MOI.add_constraint(opt, func, set)
    end
    for (func, set) in x._linear_geq_constraints
        MOI.add_constraint(opt, func, set)
    end
    for (func, set) in x._linear_eq_constraints
        MOI.add_constraint(opt, func, set)
    end

    gen_quadratic_storage!(x)

    # add a linear constraint for each nonlinear constraint
    eval_block = x._working_evaluator_block
    evaluator = eval_block.evaluator
    count_lower = 0
    count_upper = 0
    for (j, bns) in enumerate(eval_block.constraint_bounds)
        @inbounds nzidx = evaluator.constraints[j].grad_sparsity
        @inbounds nzvar = variable_index[nzidx]
        func = SAF(SAT.(zeros(length(nzidx)), nzvar), 0.0)
        if !(bns.upper == Inf)
            ci = MOI.add_constraint(opt, func, LT(0.0))
            push!(x._upper_nlp_affine, ci)
            count_upper += 1
        elseif !(bns.lower == Inf)
            ci = MOI.add_constraint(opt, func, LT(0.0))
            push!(x._lower_nlp_affine, ci)
            count_lower += 1
        end
    end

    # add an empty objective
    MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    return
end

is_lp(m::Optimizer) = ~in(true, m.branch_variable)
function linear_solve!(m::Optimizer)
    opt = m.relaxed_optimizer
    # TODO: Nonlinear terms which are actually linear
    if isa(m._objective, SV) || isa(m._objective, SAF)
        MOI.set(opt, MOI.ObjectiveFunction(), m._objective)
        MOI.optimize!(opt)
        m._solution_value = MOI.get(opt, MOI.ObjectiveValue())
        m._termination_status_code = MOI.get(opt, MOI.TerminationStatusCode())
        m._result_status_code = MOI.get(opt, MOI.ResultStatusCode())
        m._continuous_solution = MOI.get.(opt, MOI.VariablePrimal(), m._lower_variable_index)
        #m._run_time = MOI.get(opt, MOI.SolveTime())
    end
    return
end

"""
    convert_to_min!

Converts MOI.MAX_SENSE objective to equivalent MOI.MIN_SENSE objective
max(f) = - min(-f).
"""
function convert_to_min!(x::Optimizer)
    if x._optimization_sense === MOI.MAX_SENSE
        if isa(x._objective, SV)
            x._objective = SAF(SAT[SAT(-1.0, x._objective.variable)], 0.0)
        elseif isa(x._objective, SAF)
            @inbounds x._objective.terms[:] = SAT.(-getfield.(x._objective.terms, :coefficient),
                                                        getfield.(x._objective.terms, :variable_index))
            x._objective.constant *= -1.0
        elseif isa(x._objective, SQF)
            @inbounds x._objective.affine_terms[:] = SAT.(-getfield.(x._objective.affine_terms, :coefficient),
                                                               getfield.(x._objective.affine_terms, :variable_index))
            @inbounds x._objective.quadratic_terms[:] = SQT.(-getfield.(x._objective.quadratic_terms, :coefficient),
                                                                  getfield.(x._objective.quadratic_terms, :variable_index_1),
                                                                  getfield.(x._objective.quadratic_terms, :variable_index_2))
            x._objective.constant *= -1.0
        else
            nd = x._nlp_data.nlobj.nd
            pushfirst!(nd, NodeData(JuMP._Derivatives.CALLUNIVAR, 2, -1))
            for i in 2:length(nd)
                @inbounds nd[i] = NodeData(nd[i].nodetype, nd[i].index, nd[i].parent + 1)
            end
        end
    end
    return
end

# TODO
function initialize_evaluators!(m::Optimizer, flag::Bool)

    num_nlp_constraints = length(m._nlp_data.constraint_bounds)

    # Build the JuMP NLP evaluator
    evaluator = m._nlp_data.evaluator
    features = MOI.features_available(evaluator)
    has_hessian = (:Hess in features)
    init_feat = [:Grad, :ExprGraph]
    num_nlp_constraints > 0 && push!(init_feat, :Jac)
    MOI.initialize(evaluator, init_feat)

    # Scrub user-defined functions
    if ~isa(evaluator, EAGO.EmptyNLPEvaluator)
        m.presolve_scrubber_flag && Script.scrub!(evaluator.m.nlp_data)
        if m.presolve_to_JuMP_flag
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
    m._working_evaluator_block = m._nlp_data
    if ~isa(m._nlp_data.evaluator, EAGO.EmptyNLPEvaluator) || false #flag
        built_evaluator = build_nlp_evaluator(m._variable_number, m._nlp_data.evaluator, m, false)
        (m._optimization_sense == MOI.MAX_SENSE) && neg_objective!(built_evaluator)
        m._working_evaluator_block = MOI.NLPBlockData(m._nlp_data.constraint_bounds,
                                                      built_evaluator,
                                                      m._nlp_data.has_objective)
    end
    return evaluator
end

# DONE
"""
    label_fixed_variables!

Detects any variables set to a fixed value by equality or inequality constraints
and populates the _fixed_variable storage array.
"""
function label_fixed_variables!(m::Optimizer)
    lbd = 0.0
    ubd = 0.0
    @simd for i in 1:m._variable_number
        @inbounds lbd = m._variable_info[i].lower_bound
        @inbounds ubd = m._variable_info[i].upper_bound
        if lbd == ubd
            @inbounds m._variable_info[i].is_fixed = true
            @inbounds m._fixed_variable[i] = true
        end
    end
    return
end

"""
    is_convex_quadratic

Returns true if `func` < 0  based on eigenvalue tests, false otherwise.
"""
function is_convex_quadratic(func::SQF, NumVar::Int, mult::Float64)
    # Riguous Convexity Test
    Q = spzeros(NumVar, NumVar)
    for term in func.quadratic_terms
        if term.coefficient != 0.0
            Q[term.variable_index_1.value, term.variable_index_2.value] = mult*term.coefficient
        end
    end
    if length(Q.nzval) > 1
        eigval = eigmin(Array(Q))
        if (eigval) > 0.0
            return true
        end
    else
        if Q.nzval[1] > 0.0
            return true
        else
            return false
        end
    end
    return false
end


"""
    quadratic_convexity!

Assigns boolean value to constraint_convexity dictionary entry corresponding to
constraint index that is true if constraint is shown to be convex and false
otherwise.
"""
function quadratic_convexity!(x::Optimizer)
    for (func, set, ind) in x._quadratic_leq_constraints
        MOIindx = CI{SQF,LT}(ind)
        x._constraint_convexity[MOIindx] = is_convex_quadratic(func, x._variable_number, 1.0)
    end
    for (func, set, ind) in x._quadratic_geq_constraints
        MOIindx = CI{SQF, GT}(ind)
        x._constraint_convexity[MOIindx] = is_convex_quadratic(func, x.variable_number, -1.0)
    end
    for (func, set, ind) in x._quadratic_eq_constraints
        MOIindx = CI{SQF, ET}(ind)
        x.constraint_convexity[MOIindx] = false
    end
    return
end

"""
     local_solve!

Performs a single local solve of problem.
"""
function local_solve!(x::Optimizer)

    node_selection!(x)
    solve_local_nlp!(x)

    if x._upper_feasibility
        x._feasible_solution_found = true
        x._first_solution_node = x._maximum_node_id
        x._solution_value = x._upper_objective_value
        @inbounds x._continuous_solution[:] = x._upper_solution
        x._termination_status_code = MOI.LOCALLY_SOLVED
        x._result_status_code = MOI.FEASIBLE_POINT
    else
        x._feasible_solution_found = false
        x._first_solution_node = x._maximum_node_id
        x._termination_status_code = MOI.LOCALLY_INFEASIBLE
        x._result_status_code = MOI.INFEASIBLE_POINT
    end
    if m.log_on
        m._log[:total_time] = time() - m._start_time
    end
    return
end

"""
    build_nlp_evaluator

Builds the evaluator used to generate relaxations of the nonlinear equations
and constraints from a source model.
"""
function build_nlp_evaluator(N::Int, src::T, x::Optimizer, bool_flag::Bool) where {T<:MOI.AbstractNLPEvaluator}

    # Checks to see if nlp data block evaluator is present
    if ~isa(src, EAGO.EmptyNLPEvaluator)

        # Creates the empty evaluator
        d = Evaluator{N}(src.m)

        num_variables_ = JuMP.num_variables(d.m)
        d.variable_number = num_variables_
        nldata::JuMP._NLPData = deepcopy(d.m.nlp_data)

        # Copy state of user-defined multivariable functions
        d.has_user_mv_operator = src.disable_2ndorder
        d.parameter_values = nldata.nlparamvalues
        d.last_x = fill(NaN, d.variable_number)
        d.last_node = NodeBB()

        # Set valued operators must contain a (sub)gradient and/or (sub)hessian field to support 1st/2nd order eval
        d.disable_1storder = false
        d.disable_2ndorder = true

        # Add objective functions, constraints, subexpressions
        d.has_nlobj = isa(nldata.nlobj, JuMP._NonlinearExprData)
        if (src.has_nlobj)
            d.objective = copy_to_function(N, src.objective)
        end

        for i in 1:length(src.constraints)
            push!(d.constraints, copy_to_function(N, src.constraints[i]))
        end

        d.subexpression_order = src.subexpression_order
        d.subexpression_linearity = src.subexpression_linearity
        d.subexpressions_as_julia_expressions = Any[]
        if isdefined(src, :subexpressions_as_julia_expressions)
            d.subexpressions_as_julia_expressions = src.subexpressions_as_julia_expressions
        end

        d.subexpression_values_set = MC{N}[]
        d.subexpression_values_flt = Float64[]
        d.subexpressions = SubexpressionSetStorage{N}[]
        for i in 1:length(src.subexpressions)
            temp = copy_to_subexpr(N, src.subexpressions[i])
            push!(d.subexpressions,temp)
        end
        d.subexpression_values_set = fill(NaN, length(d.subexpressions))
        d.subexpression_values_flt = fill(NaN, length(d.subexpressions))

        # Add bounds to evaluator
        for bnds in x._nlp_data.constraint_bounds
            push!(d.constraints_lbd, bnds.lower)
            push!(d.constraints_ubd, bnds.upper)
        end

        # USER OUTPUT BUFFERS??????
        d.cp_tolerance = x.cp_interval_tolerance
        d.cp_reptitions = x.cp_interval_reptitions
        d.has_reverse = x._cp_evaluation_reverse
        d.subgrad_tighten = x.subgrad_tighten
        d.subgrad_tighten_reverse = x.subgrad_tighten_reverse
        d.jac_storage = Array{Float64}(undef, max(num_variables_, d.m.nlp_data.largest_user_input_dimension))

        d.constraint_number = length(d.constraints)
        d.subexpression_number = length(d.subexpressions)

        # calculate an index for each variable via search on
        d.index_to_variable = fill((-1,-1,-1), (d.variable_number,))
        for (oindx,node) in enumerate(d.objective.nd)
            if (node.nodetype == JuMP._Derivatives.VARIABLE)
                current_value = d.index_to_variable[node.index]
                if (current_value[1] == current_value[2] == current_value[3] == -1)
                    d.index_to_variable[node.index] = (oindx, 1, 1)
                end
            end
        end
        for (cindx,constraint) in enumerate(d.constraints)
            for (indx,node) in enumerate(constraint.nd)
                if (node.nodetype == JuMP._Derivatives.VARIABLE)
                    current_value = d.index_to_variable[node.index]
                    if (current_value[1] == current_value[2] == current_value[3] == -1)
                        d.index_to_variable[node.index] = (indx, cindx, 2)
                    end
                end
            end
        end
        for (cindx,subexpress) in enumerate(d.subexpressions)
            for (indx,node) in enumerate(subexpress.nd)
                if (node.nodetype == JuMP._Derivatives.VARIABLE)
                    current_value = d.index_to_variable[node.index]
                    if (current_value[1] == current_value[2] == current_value[3] == -1)
                        d.index_to_variable[node.index] = (indx, cindx, 3)
                    end
                end
            end
        end

        d.subexpression_isnum = fill(true, (d.subexpression_number,))

        return d #deepcopy(d)
    end
end

function parse_problem!(m::Optimizer)
    setrounding(Interval, m.rounding_mode)

    ########### Reformulate DAG using auxilliary variables ###########
    _variable_len = length(m._variable_info)
    m._continuous_variable_number = _variable_len
    m._variable_number = _variable_len

    ########### Set Correct Size for Problem Storage #########
    m._lower_solution = fill(0.0, _variable_len)
    m._upper_solution = fill(0.0, _variable_len)
    m._lower_lvd = fill(0.0, _variable_len)
    m._lower_uvd = fill(0.0, _variable_len)

    # Get various other sizes
    m._continuous_solution = zeros(Float64, _variable_len)

    convert_to_min!(m)
    initialize_evaluators!(m, false)               # initializes the EAGO and JuMP NLP evaluators
    return
end

"""
    create_initial_node!

Creates an initial node with initial box constraints and adds it to the stack.
"""
function create_initial_node!(m::Optimizer)
    @inbounds lower_variable_bounds = lower_bound.(m._variable_info)
    @inbounds upper_variable_bounds = upper_bound.(m._variable_info)
    n = NodeBB(lower_variable_bounds, upper_variable_bounds, -Inf, Inf, 1, 1)
    push!(m._stack, n)
    m._node_count = 1
    m._maximum_node_id += 1
    return
end

"""
    check_disable_fbbt!

Disables bound tightening for problems lacking the appropriate constraints.
"""
function check_disable_fbbt!(m::Optimizer)

    no_constraints = true
    if length(m._nlp_data.constraint_bounds) > 0
        no_constraints &= false
    end

    no_quad_constraints = true
    if length(m._quadratic_leq_constraints) > 0
        no_constraints &= false
        no_quad_constraints &= false
    end
    if length(m._quadratic_geq_constraints) > 0
        no_constraints &= false
        no_quad_constraints &= false
    end
    if length(m._quadratic_eq_constraints) > 0
        no_constraints &= false
        no_quad_constraints &= false
    end
    if no_quad_constraints
        m.quad_uni_depth = 0
        m.quad_bi_depth = 0
    end

    no_lin_constraints = true
    if length(m._linear_leq_constraints) > 0
        no_constraints &= false
        no_lin_constraints &= false
    end
    if length(m._linear_geq_constraints) > 0
        no_constraints &= false
        no_lin_constraints &= false
    end
    if length(m._linear_eq_constraints) > 0
        no_constraints &= false
        no_lin_constraints &= false
    end
    if no_lin_constraints
        m.lp_depth = 0
    end

    if no_constraints
        m.obbt_depth = 0
        m.cp_depth = 0
    end

    return
end

function presolve_problem!(m::Optimizer)

    m.presolve_epigraph_flag && reform_epigraph!(m)  # perform epigraph rearrangement
    #m.presolve_cse_flag && dag_cse_simplify!(m)      #
    #m.presolve_flatten_flag && dag_flattening!(m)

    #m = user_reformed_optimizer(m)
    #m.debug1 = initialize_evaluators!(m, true)                      # re-initializes evaluators after reformulations are performed

    create_initial_node!(m)                        # Create initial node and add it to the stack

    #m.lower_variables = MOI.VariableIndex[MOI.VariableIndex(i) for i in 1:_variable_len]
    #m.upper_variables = MOI.add_variables(m.initial_relaxed_optimizer, _variable_len)

    ###### OBBT Setup #####
    label_fixed_variables!(m)                                   # label variable fixed to a value
    _nlpdata = m._nlp_data
    _evaluator = _nlpdata.evaluator::MOI.AbstractNLPEvaluator
    label_nonlinear_variables!(m, _evaluator)

    check_disable_fbbt!(m)

    load_relaxed_problem!(m)

    return
end

function store_candidate_solution!(x::Optimizer)
    flag = false
    if x._upper_feasibility
        if x._upper_objective_value < x._global_upper_bound
            x._feasible_solution_found = true
            x._first_solution_node = x._maximum_node_id
            x._solution_value = x._upper_objective_value
            @inbounds x._continuous_solution[:] = x._upper_solution
            if x._optimization_sense == MOI.FEASIBILITY_SENSE
                if ~x.feasible_local_continue || x.local_solve_only
                    flag = true
                end
            end
        end
    end
    return flag
end

"""
    global_solve!

Solves the branch and bound problem with the input EAGO optimizer object.
"""
function global_solve!(x::Optimizer)

    ext_type = x.ext_type

    x._iteration_count = 0
    x._node_count = 1

    # terminates when max nodes or iteration is reach, or when node stack is empty
    while ~termination_check(ext_type, x)

        # Selects node, deletes it from stack, prints based on verbosity
        node_selection!(ext_type, x)
        (x.verbosity >= 3) && print_node!(x)

        # Performs prepocessing and times
        x.log_on && (start_time = time())
        preprocess!(ext_type, x)
        if x.log_on
            x._last_preprocess_time = time() - start_time
        end

        if x._preprocess_feasibility

            # solves & times lower bounding problem
            x.log_on && (start_time = time())
            x._cut_iterations = 0
            lower_problem!(ext_type, x)
            while cut_condition(ext_type, x)
                add_cut!(ext_type, x)
            end
            if x.log_on
                x._last_lower_problem_time = time() - start_time
            end
            print_results!(x, true)
            print_results_post_cut!(x)

            # checks for infeasibility stores solution
            if x._lower_feasibility

                if ~convergence_check(ext_type, x)

                    x.log_on && (start_time = time())
                    upper_problem!(ext_type, x)
                    if x.log_on
                        x._last_upper_problem_time = time() - start_time
                    end
                    print_results!(x, false)
                    store_candidate_solution!(x) && break

                    # Performs and times post processing
                    x.log_on && (start_time = time())
                    postprocess!(ext_type, x)
                    if x.log_on
                        x._last_postprocessing_time = time() - start_time
                    end

                    # Checks to see if the node
                    if (x._postprocess_feasibility)
                        if repeat_check(ext_type, x)
                            single_storage!(ext_type, x)
                        else
                            branch_node!(ext_type, x)
                        end
                    end
                end
            end
            fathom!(ext_type, x)
        else
            x._lower_objective_value = -Inf
            x._lower_feasibility = false
            x._upper_feasibility = false
        end
        x._run_time = time() - x._start_time
        x._iteration_count += 1
        log_iteration!(x)
        print_iteration!(x)
    end

    x._objective_value = x._global_upper_bound

    # Prints the solution
    print_solution!(x)
end

function MOI.optimize!(m::Optimizer)

    m._start_time = time()
    parse_problem!(m)
    if m.log_on
        m._log[:parse_time] = time() - start_time
        m._log[:total_time] = m._log[:parse_time]
    end

    presolve_problem!(m)
    if m.log_on
        m._log[:presolve_time] = time() - m._log[:parse_time]
    end

    # Allow for a hook to modify branch and bound routine
    if m.enable_optimize_hook
        return optimize_hook(m.ext_type, m)
    end

    # Runs the branch and bound routine
    if is_lp(m)
        linear_solve!(m)
    elseif m.local_solve_only
        local_solve!(m)
    else
        global_solve!(m)
    end
    return
end
