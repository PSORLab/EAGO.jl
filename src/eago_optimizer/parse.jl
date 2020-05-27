function create_working_functions!(m::Optimizer)

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

Converts MOI.MAX_SENSE objective to equivalent MOI.MIN_SENSE objective
max(f) = - min(-f).
"""
function convert_to_min!(m::Optimizer)

    if m._input_problem._optimization_sense === MOI.MAX_SENSE

        if m._objective_type === SINGLE_VARIABLE
            m._objective_type = SCALAR_AFFINE
            m._objective_saf = SAF(SAT[SAT(-1.0, m._objective_sv.variable)], 0.0)

        elseif m._objective_type === SCALAR_AFFINE
            @__dot__ m._objective_saf.terms = SAT(-getfield(m._objective_saf.terms, :coefficient),
                                                   getfield(m._objective_saf.terms, :variable_index))
            m._objective_saf.constant *= -1.0

        elseif m._objective_type === SCALAR_QUADRATIC
            @__dot__ m._objective_sqf.affine_terms = SAT(-getfield(m._objective_sqf.affine_terms, :coefficient),
                                                          getfield(m._objective_sqf.affine_terms, :variable_index))
            @__dot__ m._objective_sqf.quadratic_terms = SQT(-getfield(m._objective_sqf.quadratic_terms, :coefficient),
                                                             getfield(m._objective_sqf.quadratic_terms, :variable_index_1),
                                                             getfield(m._objective_sqf.quadratic_terms, :variable_index_2))
            m._objective_sqf.constant *= -1.0

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

function check_set_is_fixed(v::VariableInfo)
    v.is_fixed && return true
    v.is_fixed = x.lower_bound === x.upper_bound
    return v.is_fixed
end

"""
$(TYPEDSIGNATURES)

Detects any variables set to a fixed value by equality or inequality constraints
and populates the _fixed_variable storage array.
"""
function label_fixed_variables!(m::Optimizer)
    map!(x -> check_set_is_fixed(x), m._fixed_variable, m._variable_info)
end

function label_branch_variables!(m::Optimizer)

    m._user_branch_variables = !isempty(m.parameters.branch_variable)
    if m._user_branch_variables
        copyto!(m._branch_variables, m.parameters.branch_variable)
    end

    # adds nonlinear terms in quadratic constraints
    sqf_leq = m._working_problem._sqf_leq
    for i = 1:m._working_problem._sqf_leq_count
        quad_ineq = @inbounds sqf_leq[i]
        for term in quad_ineq.sqf
            variable_index_1 = term.variable_index_1.value
            variable_index_2 = term.variable_index_2.value
            @inbounds m._branch_variables[variable_index_1] = true
            @inbounds m._branch_variables[variable_index_2] = true
        end
    end

    sqf_eq = m._working_problem._sqf_eq
    for i = 1:m._working_problem._sqf_eq_count
        quad_eq = @inbounds sqf_eq[i]
        for term in quad_eq.sqf
            variable_index_1 = term.variable_index_1.value
            variable_index_2 = term.variable_index_2.value
            @inbounds m._branch_variables[variable_index_1] = true
            @inbounds m._branch_variables[variable_index_2] = true
        end
    end

    # adds nonlinear terms in objectives if
    obj_type = m._working_problem._objective_type
    if obj_type === SCALAR_QUADRATIC
        for term in m._working_problem._objective_sqf
            variable_index_1 = term.variable_index_1.value
            variable_index_2 = term.variable_index_2.value
            @inbounds m._branch_variables[variable_index_1] = true
            @inbounds m._branch_variables[variable_index_2] = true
        end
    end

    return
end

"""
Translates input problem to working problem. Routines and checks and optional manipulation is left to the presolve stage.
"""
function parse_problem!(m::Optimizer)

    m._time_left = m.time_limit
    ip = m._input_problem

    # add linear constraints to the working problem
    for i = 1:ip._linear_leq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._linear_leq_constraints[i])
    end
    for i = 1:ip._linear_geq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._linear_geq_constraints[i])
    end
    for i = 1:ip._linear_eq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._linear_eq_constraints[i])
    end

    # add quadratic constraints to the working problem
    for i = 1:ip._quadratic_leq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._quadratic_leq_constraints[i])
    end
    for i = 1:ip._quadratic_geq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._quadratic_geq_constraints[i])
    end
    for i = 1:ip._quadratic_eq_count
        add_to_working_problem!(m._working_problem, @inbounds ip._quadratic_eq_constraints[i])
    end

    # converts a maximum problem to a minimum problem (internally) if necessary
    convert_to_min!(m)

    # labels the variable info and the _fixed_variable vector for each fixed variable
    label_fixed_variables!(m)

    # labels variables to branch on
    label_branch_variables!(m)

    new_time = time() - m._start_time
    m._parse_time = new_time
    m._run_time = new_time

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

### Routines for parsing the full nonconvex problem
"""
Reformulates quadratic terms in SOC constraints if possible. For <= or >=,
the quadratic term is deleted if an SOCP is detected. For ==, the SOC check
is done for each >= and <=, the convex constraint is reformulated to a SOC,
the concave constraint is keep as a quadratic.
"""
function parse_classify_quadratic!(m::Optimizer)
    #=
    for (id, cinfo) in m._quadratic_constraint
        is_soc, add_concave, cfunc, cset, qfunc, qset = check_convexity(cinfo.func, cinfo.set)
        if is_soc
            MOI.add_constraint(m, cfunc, cset)
            deleteat!(m._quadratic_constraint, id)
            if add_concave
                MOI.add_constraint(m, qfunc, qset)
            end
        end
    end
    =#
    nothing
end

"""
"""
function parse_classify_nlp(m)
    nothing
end

"""
Classifies the problem type
"""
function parse_classify_problem!(m::Optimizer)

    ip = m._input_problem
    integer_variable_number = count(is_integer.(ip._variable_info))
    cone_constraint_number = ip._conic_second_order_count
    quad_constraint_number = ip._quadratic_leq_count + ip._quadratic_geq_count + ip._quadratic_eq_count

    if integer_variable_number === 0
        if cone_constraint_number === 0 && quad_constraint_number === 0
            # && iszero(m._input_nonlinear_constraint_number)
            m._working_problem._problem_type = LP
        elseif quad_constraint_number === 0
            # && iszero(m._input_nonlinear_constraint_number)
            m._working_problem = SOCP
        else
            #parse_classify_quadratic!(m)
            #if iszero(m._input_nonlinear_constraint_number)
            #    if isempty(m._quadratic_constraint)
            #        m._problem_type = SOCP
            #    end
            #else
            #    # Check if DIFF_CVX, NS_CVX, DIFF_NCVX, OR NS_NCVX
            #    m._problem_type = parse_classify_nlp(m)
            #end
            m._problem_type = NS_NCVX
        end
    else
        if cone_constraint_number === 0 && quad_constraint_number === 0
            m._problem_type = MILP
        elseif quad_constraint_number === 0
            m._problem_type = MISOCP
        else
            parse_classify_quadratic!(m)
            #=
            if iszero(m._nonlinear_constraint_number)
                if iszero(m._quadratic_constraint_number)
                    m._problem_type = MISOCP
                end
            else
                # Performs parsing
                _ = parse_classify_nlp(m)
            end
            =#
            m._problem_type = MINCVX
        end
    end

    return nothing
end

#=
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
