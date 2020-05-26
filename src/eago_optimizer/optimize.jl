function create_inner_functions!(m::Optimizer)

    # add buffered quadratic constraints
    for i = 1:m._input_problem._quadratic_leq_count
        sqf, lt = m._input_problem._quadratic_leq_constraints
        add_buffered_quadratic!(m, sqf, lt)
    end

    # add buffered quadratic constraints
    for i = 1:m._input_problem._quadratic_geq_count
        sqf, gt = m._input_problem._quadratic_geq_constraints
        add_buffered_quadratic!(m, sqf, gt)
    end

    # add buffered quadratic constraints
    for i = 1:m._input_problem._quadratic_eq_count
        sqf, et = m._input_problem._quadratic_eq_constraints
        add_buffered_quadratic!(m, sqf, et)
    end

    nothing
end

"""
$(TYPEDSIGNATURES)

Loads variables, linear constraints, and empty storage for first nlp and
quadratic cut.
"""
function load_relaxed_problem!(m::Optimizer)
    relaxed_optimizer = m.relaxed_optimizer

    # add variables and indices
    variable_count = m._working_problem._variable_count
    for i = 1:variable_count
        push!(x._working_problem._relaxed_variable_index, MOI.add_variable(opt))
    end

    # add variables
    for i = 1:variable_count
        @inbounds vinfo = x._variable_info[i]
        single_variable = SV(x._working_problem._relaxed_variable_index[i])
        if vinfo.is_integer
        elseif vinfo.is_fixed
            ci_sv_et = MOI.add_constraint(relaxed_optimizer, single_variable, ET(vinfo.lower_bound))
            m._relaxed_variable_node_map[ci_sv_et] = i
        else
            if vinfo.has_lower_bound
                ci_sv_gt = MOI.add_constraint(relaxed_optimizer, single_variable, GT(vinfo.lower_bound))
                m._relaxed_variable_node_map[ci_sv_gt] = i
            end
            if vinfo.has_upper_bound
                ci_sv_lt = MOI.add_constraint(relaxed_optimizer, single_variable, LT(vinfo.upper_bound))
                m._relaxed_variable_node_map[ci_sv_lt] = i
            end
        end
    end

    # add linear constraints
    for (func, set) in m._input_problem._linear_leq_constraints
         MOI.add_constraint(opt, func, set)
    end
    for (func, set) in m._input_problem._linear_geq_constraints
        MOI.add_constraint(opt, func, set)
    end
    for (func, set) in m._input_problem._linear_eq_constraints
        MOI.add_constraint(opt, func, set)
    end

    MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    return
end

"""
$(TYPEDSIGNATURES)

Converts MOI.MAX_SENSE objective to equivalent MOI.MIN_SENSE objective
max(f) = - min(-f).
"""
function convert_to_min!(x::Optimizer)

    if x._input_problem._optimization_sense === MOI.MAX_SENSE

        if x._objective_type === SINGLE_VARIABLE
            x._objective_type = SCALAR_AFFINE
            x._objective_saf = SAF(SAT[SAT(-1.0, x._objective_sv.variable)], 0.0)

        elseif x._objective_type === SCALAR_AFFINE
            @__dot__ x._objective_saf.terms = SAT(-getfield(x._objective_saf.terms, :coefficient),
                                                   getfield(x._objective_saf.terms, :variable_index))
            x._objective_saf.constant *= -1.0

        elseif x._objective_type === SCALAR_QUADRATIC
            @__dot__ x._objective_sqf.affine_terms = SAT(-getfield(x._objective_sqf.affine_terms, :coefficient),
                                                          getfield(x._objective_sqf.affine_terms, :variable_index))
            @__dot__ x._objective_sqf.quadratic_terms = SQT(-getfield(x._objective_sqf.quadratic_terms, :coefficient),
                                                             getfield(x._objective_sqf.quadratic_terms, :variable_index_1),
                                                             getfield(x._objective_sqf.quadratic_terms, :variable_index_2))
            x._objective_sqf.constant *= -1.0

        else
            #=
            nd = x._nlp_data.evaluator.m.nlp_data.nlobj.nd
            pushfirst!(nd, NodeData(JuMP._Derivatives.CALLUNIVAR, 2, -1))
            nd[2] = NodeData(nd[2].nodetype, nd[2].index, 1)
            for i = 3:length(nd)
                @inbounds nd[i] = NodeData(nd[i].nodetype, nd[i].index, nd[i].parent + 1)
            end
            =#
        end
    end
    return
end

"""
$(TYPEDSIGNATURES)

Detects any variables set to a fixed value by equality or inequality constraints
and populates the _fixed_variable storage array.
"""
function label_fixed_variables!(m::Optimizer)
    map!(m -> (m.is_fixed |= (m.lower_bound == m.upper_bound)), m._fixed_variable, m._variable_info)
end

"""
$(TYPEDSIGNATURES)

Performs a single local solve of problem.
"""
function local_solve!(m::Optimizer)

    node_selection!(m)
    solve_local_nlp!(m)

    if m._upper_feasibility
        m._feasible_solution_found = true
        m._first_solution_node = m._maximum_node_id
        m._solution_value = m._upper_objective_value
        m._continuous_solution .= m._upper_solution
        m._objective_value = m._upper_objective_value
        m._termination_status_code = MOI.LOCALLY_SOLVED
        m._result_status_code = MOI.FEASIBLE_POINT

    else
        m._feasible_solution_found = false
        m._first_solution_node = m._maximum_node_id
        m._termination_status_code = MOI.LOCALLY_INFEASIBLE
        m._result_status_code = MOI.INFEASIBLE_POINT

    end

    return nothing
end

"""
Translates input problem to working problem. Routines and checks and optional manipulation is left to the presolve stage.
"""
function parse_problem!(m::Optimizer)

    m._user_branch_variables = ~isempty(m.branch_variable)
    m._time_left = m.time_limit

    ip = m._input_problem
    for i = 1:ip._linear_leq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._linear_leq_constraints[i])
    end
    for i = 1:ip._linear_geq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._linear_geq_constraints[i])
    end
    for i = 1:ip._linear_eq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._linear_eq_constraints[i])
    end

    for i = 1:ip._quadratic_leq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._quadratic_leq_constraints[i])
    end
    for i = 1:ip._quadratic_geq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._quadratic_geq_constraints[i])
    end
    for i = 1:ip._quadratic_eq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._quadratic_eq_constraints[i])
    end

    convert_to_min!(m)
    label_fixed_variables!(m)
    label_nonlinear_variables!(m)

    new_time = time() - m._start_time
    m._parse_time = new_time
    m._run_time = new_time

    return nothing
end

"""
$(TYPEDSIGNATURES)

Creates an initial node with initial box constraints and adds it to the stack.
"""
function create_initial_node!(m::Optimizer)
    n = NodeBB(lower_bound.(m._variable_info), upper_bound.(m._variable_info), -Inf, Inf, 1, 1)
    push!(m._stack, n)
    m._node_count = 1
    m._maximum_node_id += 1
    return nothing
end

function presolve_problem!(m::Optimizer)

    ########### Reformulate DAG using auxilliary variables ###########
    _variable_len = length(m._variable_info)
    m._variable_number = _variable_len

    ########### Set Correct Size for Problem Storage #########
    m._current_xref = fill(0.0, _variable_len)
    m._cut_solution = fill(0.0, _variable_len)
    m._lower_solution = fill(0.0, _variable_len)
    m._upper_solution = fill(0.0, _variable_len)
    m._lower_lvd = fill(0.0, _variable_len)
    m._lower_uvd = fill(0.0, _variable_len)
    m._continuous_solution = zeros(Float64, _variable_len)

    m.presolve_epigraph_flag && reform_epigraph!(m)  # perform epigraph rearrangement
    m.presolve_cse_flag && dag_cse_simplify!(m)      #

    create_initial_node!(m)                          # Create initial node and add it to the stack
    load_relaxed_problem!(m)

    m._presolve_time = time() - m._parse_time

    return
end

function store_candidate_solution!(m::Optimizer)
    if m._upper_feasibility && (m._upper_objective_value < m._global_upper_bound)
        m._feasible_solution_found = true
        m._first_solution_node = m._maximum_node_id
        m._solution_value = m._upper_objective_value
        m._global_upper_bound = m._upper_objective_value
        @__dot__ m._continuous_solution = m._upper_solution
    end
    return
end

function set_global_lower_bound!(m::Optimizer)
    if !isempty(m._stack)
        min_node = minimum(m._stack)
        lower_bound = min_node.lower_bound
        if m._global_lower_bound < lower_bound
            m._global_lower_bound = lower_bound
        end
    end
    return
end

# wraps subroutine call to isolate ExtensionType
termination_check(m::Optimizer) = termination_check(m.ext_type, m)
cut_condition(m::Optimizer) = cut_condition(m.ext_type, m)
convergence_check(m::Optimizer) = convergence_check(m.ext_type, m)
repeat_check(m::Optimizer) = repeat_check(m.ext_type, m)
node_selection!(m::Optimizer) = node_selection!(m.ext_type, m)
preprocess!(m::Optimizer) = preprocess!(m.ext_type, m)
lower_problem!(m::Optimizer) = lower_problem!(m.ext_type, m)
add_cut!(m::Optimizer) = add_cut!(m.ext_type, m)
upper_problem!(m::Optimizer) = upper_problem!(m.ext_type, m)
postprocess!(m::Optimizer) = postprocess!(m.ext_type, m)
single_storage!(m::Optimizer) = single_storage!(m.ext_type, m)
branch_node!(m::Optimizer) = branch_node!(m.ext_type, m)
fathom!(m::Optimizer) = fathom!(m.ext_type, m)

"""
$(TYPEDSIGNATURES)

Solves the branch and bound problem with the input EAGO optimizer object.
"""
function global_solve!(x::Optimizer)

    x._iteration_count = 1
    x._node_count = 1

    # terminates when max nodes or iteration is reach, or when node stack is empty
    while !termination_check(x)

        # Selects node, deletes it from stack, prints based on verbosity
        node_selection!(x)
        (x.verbosity >= 3) && print_node!(x)

        # Performs prepocessing and times
        x.log_on && (start_time = time())
        preprocess!(x)
        if x.log_on
            x._last_preprocess_time = time() - start_time
        end

        if x._preprocess_feasibility

            # solves & times lower bounding problem
            x.log_on && (start_time = time())
            x._cut_iterations = 1
            lower_problem!(x)
            while cut_condition(x)
                add_cut!(x)
            end
            if x.log_on
                x._last_lower_problem_time = time() - start_time
            end
            print_results!(x, true)
            print_results_post_cut!(x)

            # checks for infeasibility stores solution
            if x._lower_feasibility
                if ~convergence_check(x)

                    x.log_on && (start_time = time())
                    upper_problem!(x)
                    if x.log_on
                        x._last_upper_problem_time = time() - start_time
                    end
                    print_results!(x, false)
                    store_candidate_solution!(x)
                    if x._input_problem._optimization_sense === MOI.FEASIBILITY_SENSE
                        if ~x.feasible_local_continue || x.local_solve_only
                            break
                        end
                    end

                    # Performs and times post processing
                    x.log_on && (start_time = time())
                    postprocess!(x)
                    if x.log_on
                        x._last_postprocessing_time = time() - start_time
                    end

                    # Checks to see if the node
                    if (x._postprocess_feasibility)
                        if repeat_check(x)
                            single_storage!(x)
                        else
                            branch_node!(x)
                        end
                    end
                else
                    #x._global_lower_bound = x._lower_objective_value
                end
            end
            fathom!(x)
        else
            x._lower_objective_value = -Inf
            x._lower_feasibility = false
            x._upper_feasibility = false
        end
        set_global_lower_bound!(x)
        x._run_time = time() - x._start_time
        x._time_left = x.time_limit - x._run_time
        log_iteration!(x)
        print_iteration!(x)
        x._iteration_count += 1
    end

    x._objective_value = x._global_upper_bound

    # Prints the solution
    print_solution!(x)
    return
end

throw_optimize_hook!(m::Optimizer) = optimize_hook!(m.ext_type, m)
function MOI.optimize!(m::Optimizer)

    m._start_time = time()
    presolve_problem!(m)
    parse_problem!(m)

    # Runs the branch and bound routine
    if !m.enable_optimize_hook
        if is_lp(m)
            linear_solve!(m)
        elseif m.local_solve_only
            local_solve!(m)
        else
            global_solve!(m)
        end
    else
        throw_optimize_hook!(m)
    end

    return nothing
end

#=
function label_nonlinear_variables!(m::Optimizer)
    _nlpdata = m._nlp_data
    x = _nlpdata.evaluator

    # scans subexpressions, objective, and constraints for nonlinear terms
    if ~isa(x, EmptyNLPEvaluator)
        if x.has_nlobj
            if (x.objective.linearity != JuMP._Derivatives.LINEAR) &&
               (x.objective.linearity != JuMP._Derivatives.CONSTANT)
                for i in 1:length(x.objective.nd)
                    nd = x.objective.nd[i]
                    if (nd.nodetype == JuMP._Derivatives.VARIABLE)
                        m.branch_variable[nd.index] = true
                        m.obbt_variable_values[nd.index] = true
                    end
                end
            end
        end
        for i in 1:length(x.constraints)
            if (x.constraints[i].linearity != JuMP._Derivatives.LINEAR) &&
               (x.constraints[i].linearity != JuMP._Derivatives.CONSTANT)
                for j in 1:length(x.constraints[i].nd)
                    nd = x.constraints[i].nd[j]
                    bool1 = (nd.nodetype == JuMP._Derivatives.VARIABLE)
                    if (nd.nodetype == JuMP._Derivatives.VARIABLE)
                        m.branch_variable[nd.index] = true
                        m.obbt_variable_values[nd.index] = true
                    end
                end
            end
        end
        for i in 1:length(x.subexpressions)
            if (x.subexpressions[i].linearity != JuMP._Derivatives.LINEAR) &&
               (x.subexpressions[i].linearity != JuMP._Derivatives.CONSTANT)
                for j in 1:length(x.subexpressions[i].nd)
                    nd = x.subexpressions[i].nd[j]
                    bool1 = (nd.nodetype == JuMP._Derivatives.VARIABLE)
                    if (nd.nodetype == JuMP._Derivatives.VARIABLE)
                        m.branch_variable[nd.index] = true
                        m.obbt_variable_values[nd.index] = true
                    end
                end
            end
        end
    end
    return
end
=#

#=
function build_nlp_kernel!(d::Evaluator{N,T}, src::JuMP.NLPEvaluator, x::Optimizer, bool_flag::Bool) where {N,T<:RelaxTag}

    m = src.m::Model
    num_variables_ = JuMP.num_variables(m)
    d.variable_number = num_variables_
    #nldata = m.nlp_data::JuMP._NLPData
    nldata = deepcopy(m.nlp_data)
    parameter_values = nldata.nlparamvalues

    # Copy state of user-defined multivariable functions
    d.has_user_mv_operator = src.disable_2ndorder
    d.last_x = fill(NaN, d.variable_number)
    d.last_node = NodeBB()

    # Set valued operators must contain a (sub)gradient and/or (sub)hessian field to support 1st/2nd order eval
    d.disable_1storder = false

    # Add objective functions, constraints, subexpressions
    d.has_nlobj = src.has_nlobj
    if src.has_nlobj
        copy_to_function!(d, 1, src.objective)
    end

    for i in 1:length(src.constraints)
        copy_to_function!(d, i + 1, src.constraints[i])
    end

    d.subexpression_order = src.subexpression_order
    d.subexpression_linearity = src.subexpression_linearity

    d.subexpression_values_set = MC{N,T}[]
    d.subexpression_values_flt = Float64[]
    d.subexpressions = SubexpressionSetStorage{N}[]
    for i in 1:length(src.subexpressions)
        copy_to_subexpr!(d, src.subexpressions[i])
    end
    d.subexpression_values_set = fill(zero(MC{N,T}), length(d.subexpressions))
    d.subexpression_values_flt = fill(NaN, length(d.subexpressions))

    # Add bounds to evaluator
    for bnds in x._nlp_data.constraint_bounds
        push!(d.constraints_lbd, bnds.lower)
        push!(d.constraints_ubd, bnds.upper)
    end

    d.cp_tolerance = x.cp_tolerance
    d.cp_repetitions = x.cp_repetitions
    d.has_reverse = x._cp_evaluation_reverse
    d.subgrad_tighten = x.subgrad_tighten
    d.subgrad_tighten_reverse = x.subgrad_tighten_reverse
    d.jac_storage = fill(zero(MC{N,T}), max(num_variables_, nldata.largest_user_input_dimension))
    d.flt_jac_storage = fill(0.0, max(num_variables_, nldata.largest_user_input_dimension))

    d.constraint_number = length(d.constraints)
    d.subexpression_number = length(d.subexpressions)

    # calculate an index for each variable via search on
    unvisited_tuple = (-1,-1,-1)
    d.index_to_variable = fill(unvisited_tuple, (d.variable_number,))
    for (indx, node) in enumerate(d.objective.nd)
        op = node.index
        ntype = node.nodetype
        if (ntype == JuMP._Derivatives.VARIABLE)
            current_value = d.index_to_variable[op]
            if (current_value == unvisited_tuple)
                d.index_to_variable[op] = (indx, 1, 1)
            end
            @inbounds d.objective.numvalued[indx] = false
        elseif ntype == JuMP._Derivatives.VALUE
            @inbounds d.objective.numberstorage[indx] = d.objective.const_values[op]
            @inbounds d.objective.numvalued[indx] = true
        elseif ntype == JuMP._Derivatives.PARAMETER
            @inbounds d.objective.numberstorage[indx] = parameter_values[op]
            @inbounds d.objective.numvalued[indx] = true
        end
    end
    for (cindx,constraint) in enumerate(d.constraints)
        for (indx,node) in enumerate(constraint.nd)
            op = node.index
            if node.nodetype == JuMP._Derivatives.VARIABLE
                current_value = d.index_to_variable[op]
                if (current_value == unvisited_tuple)
                    d.index_to_variable[op] = (indx, cindx, 2)
                end
                @inbounds constraint.numvalued[indx] = false
            elseif node.nodetype == JuMP._Derivatives.VALUE
                @inbounds constraint.numberstorage[indx] = constraint.const_values[op]
                @inbounds constraint.numvalued[indx] = true
            elseif node.nodetype == JuMP._Derivatives.PARAMETER
                @inbounds constraint.numberstorage[indx] = parameter_values[op]
                @inbounds constraint.numvalued[indx] = true
            end
        end
    end
    for (cindx,subexpress) in enumerate(d.subexpressions)
        for (indx,node) in enumerate(subexpress.nd)
            op = node.index
            if node.nodetype == JuMP._Derivatives.VARIABLE
                current_value = d.index_to_variable[op]
                if (current_value == unvisited_tuple)
                    d.index_to_variable[op] = (indx, cindx, 3)
                end
                @inbounds subexpress.numvalued[indx] = false
            elseif node.nodetype == JuMP._Derivatives.VALUE
                @inbounds subexpress.numberstorage[indx] = subexpress.const_values[op]
                @inbounds subexpress.numvalued[indx] = true
            elseif node.nodetype == JuMP._Derivatives.PARAMETER
                @inbounds subexpress.numberstorage[indx] = parameter_values[op]
                @inbounds subexpress.numvalued[indx] = true
            end
        end
    end

    d.subexpression_isnum = fill(true, (d.subexpression_number,))

    d.user_operators = nldata.user_operators
    return
end

"""
$(TYPEDSIGNATURES)

Builds the evaluator used to generate relaxations of the nonlinear equations
and constraints from a source model.
"""
function build_nlp_evaluator(N::Int64, s::T, src::JuMP.NLPEvaluator, x::Optimizer, bool_flag::Bool) where {T<:RelaxTag}

    # Creates the empty evaluator
    d::Evaluator{N,T} = Evaluator{N,T}()
    build_nlp_kernel!(d, src, x, bool_flag)

    return d
end
=#

#=
"""
$(TYPEDSIGNATURES)

Disables bound tightening for problems lacking the appropriate constraints.
"""
function check_disable_fbbt!(m::Optimizer)

    no_constraints = true
    if length(m._nlp_data.constraint_bounds) > 0
        no_constraints &= false
    end
    if no_constraints
        m.cp_depth = -1
    end

    no_quad_constraints = true
    if m._input_problem._quadratic_leq_count !== 0
        no_constraints &= false
        no_quad_constraints &= false
    end
    if m._input_problem._quadratic_geq_count !== 0
        no_constraints &= false
        no_quad_constraints &= false
    end
    if m._input_problem._quadratic_eq_count !== 0
        no_constraints &= false
        no_quad_constraints &= false
    end
    if no_quad_constraints
        m.quad_uni_depth = -1
        m.quad_bi_depth = -1
    end

    only_lin_constraints = no_constraints
    no_lin_constraints = true
    if m._input_problem._linear_leq_count !== 0
        no_lin_constraints &= false
    end
    if m._input_problem._linear_geq_count !== 0
        no_lin_constraints &= false
    end
    if m._input_problem._linear_eq_count !== 0
        no_lin_constraints &= false
    end
    if no_lin_constraints
        m.lp_depth = -1
    end

    if only_lin_constraints
        m.obbt_depth = -1
    end

    return
end
=#

#=
function has_evaluator(x::MOI.NLPBlockData)
    flag = x.evaluator !== nothing
    flag &= ~isa(x.evaluator, EmptyNLPEvaluator)
    return flag
end

function initialize_scrub!(m::Optimizer, y::JuMP.NLPEvaluator)
    m.presolve_scrubber_flag && Script.scrub!(y.m.nlp_data)
    if m.presolve_to_JuMP_flag
        Script.udf_loader!(m)
    end
    return
end

function initialize_evaluators!(m::Optimizer, flag::Bool)

    nlp_data = deepcopy(m._nlp_data)

    has_eval = has_evaluator(nlp_data)
    if has_evaluator(nlp_data)

        # Build the JuMP NLP evaluator
        evaluator = nlp_data.evaluator::JuMP.NLPEvaluator
        num_nlp_constraints = length(nlp_data.constraint_bounds)
        features = MOI.features_available(evaluator)
        has_hessian = (:Hess in features)
        init_feat = [:Grad, :ExprGraph]
        num_nlp_constraints > 0 && push!(init_feat, :Jac)
        MOI.initialize(evaluator, init_feat)
        m._nlp_data = nlp_data

        # Scrub user-defined functions
        initialize_scrub!(m, evaluator)

        built_evaluator = build_nlp_evaluator(m._variable_number, NS(), deepcopy(evaluator), m, false)
        m._relaxed_evaluator = built_evaluator
        m._relaxed_eval_has_objective = m._nlp_data.has_objective
        append!(m._relaxed_constraint_bounds, m._nlp_data.constraint_bounds)

        # add info to guard context
        m._relaxed_evaluator.ctx = GuardCtx(metadata = GuardTracker(m.domain_violation_ϵ))
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

    return
end
=#

is_lp(m::Optimizer) = ~in(true, m.branch_variable)

function linear_solve!(m::Optimizer)

    opt = m.relaxed_optimizer
    if m._objective_type === SINGLE_VARIABLE
        MOI.set(opt, MOI.ObjectiveFunction{SV}(), m._objective_sv)
    elseif  m._objective_type === SCALAR_AFFINE
        MOI.set(opt, MOI.ObjectiveFunction{SAF}(), m._objective_saf)
    end

    MOI.optimize!(opt)
    m._objective_value = MOI.get(opt, MOI.ObjectiveValue())
    m._solution_value = MOI.get(opt, MOI.ObjectiveValue())
    m._global_lower_bound = MOI.get(opt, MOI.ObjectiveValue())
    m._global_upper_bound = MOI.get(opt, MOI.ObjectiveValue())
    m._termination_status_code = MOI.get(opt, MOI.TerminationStatus())
    m._result_status_code = MOI.get(opt, MOI.PrimalStatus())
    m._continuous_solution = MOI.get.(opt, MOI.VariablePrimal(), m._lower_variable_index)
    #m._run_time = MOI.get(opt, MOI.SolveTime())

    return
end
