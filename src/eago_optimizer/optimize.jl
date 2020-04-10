"""
    gen_quad_sparsity

Generate the list of variables (1,...,n) that are included in the func.
"""
function gen_quad_vals(func::SQF)
    list = Int64[]
    for term in func.affine_terms
        val = term.variable_index.value
        push!(list, val)
    end
    for term in func.quadratic_terms
        val1 = term.variable_index_1.value
        val2 = term.variable_index_2.value
        push!(list, val1, val2)
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

    temp_leq_ci = CI{SAF,LT}[]
    for (func, set, ind) in x._quadratic_leq_constraints

        v = gen_quad_vals(func)
        nv = length(v)
        @inbounds nzvar = variable_index[v]
        push!(x._quadratic_leq_sparsity, nzvar)
        push!(x._quadratic_leq_gradnz, nv)

        d = ImmutableDict{Int64,Int64}()
        for (i, val) in enumerate(v)
            d = ImmutableDict(d, val => i)
        end
        push!(x._quadratic_leq_dict, d)

        func = SAF(SAT.(zeros(nv), nzvar), 0.0)
        ci = MOI.add_constraint(opt, func, LT(0.0))
        push!(temp_leq_ci, ci)
    end

    len_leq_ci = length(temp_leq_ci)
    push!(x._quadratic_ci_leq, temp_leq_ci)
    for i in 2:x.cut_max_iterations
        push!(x._quadratic_ci_leq, fill(CI{SAF,LT}(-1), (len_leq_ci,)))
    end

    temp_geq_ci = CI{SAF,LT}[]
    for (func, set, ind) in x._quadratic_geq_constraints

        v = gen_quad_vals(func)
        nv = length(v)
        @inbounds nzvar = variable_index[v]
        push!(x._quadratic_geq_sparsity, nzvar)
        push!(x._quadratic_geq_gradnz, nv)

        d = ImmutableDict{Int64,Int64}()
        for (i, val) in enumerate(v)
            d = ImmutableDict(d, val => i)
        end
        push!(x._quadratic_geq_dict, d)

        func = SAF(SAT.(zeros(nv), nzvar), 0.0)
        ci = MOI.add_constraint(opt, func, LT(0.0))
        push!(temp_geq_ci, ci)
    end

    len_geq_ci = length(temp_geq_ci)
    push!(x._quadratic_ci_geq, temp_geq_ci)
    for i in 2:x.cut_max_iterations
        push!(x._quadratic_ci_geq, fill(CI{SAF,LT}(-1), (len_geq_ci,)))
    end

    temp_eq_ci = Tuple{CI{SAF,LT},CI{SAF,LT}}[]
    for (func, set, ind) in x._quadratic_eq_constraints

        v = gen_quad_vals(func)
        nv = length(v)
        @inbounds nzvar = variable_index[v]
        push!(x._quadratic_eq_sparsity, nzvar)
        push!(x._quadratic_eq_gradnz, nv)

        d = ImmutableDict{Int64,Int64}()
        for (i, val) in enumerate(v)
            d = ImmutableDict(d, val => i)
        end
        push!(x._quadratic_eq_dict, d)

        func = SAF(SAT.(zeros(nv), nzvar), 0.0)
        c1 = MOI.add_constraint(opt, func, LT(0.0))
        c2 = MOI.add_constraint(opt, func, LT(0.0))
        push!(temp_eq_ci, (c1, c2))
    end

    len_eq_ci = length(temp_eq_ci)
    push!(x._quadratic_ci_eq, temp_eq_ci)
    for i in 2:x.cut_max_iterations
        push!(x._quadratic_ci_eq, fill((CI{SAF,LT}(-1), CI{SAF,LT}(-1)), (len_eq_ci,)))
    end

    if x._objective_is_sqf
        v = gen_quad_vals(x._objective_sqf)
        d = ImmutableDict{Int64,Int64}()
        for (i, val) in enumerate(v)
            d = ImmutableDict(d, val => i)
        end
        x._quadratic_obj_dict = d
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
    variable_number = x._variable_number
    for i in 1:variable_number
        variable_index = MOI.add_variable(opt)
        single_variable = SV(variable_index)
        push!(x._lower_variable_index, variable_index)
        push!(x._lower_variable, single_variable) # FIX ME WHEN
    end

    # add variables
    variable = x._lower_variable
    for i in 1:variable_number

        @inbounds var = x._variable_info[i]
        @inbounds variable_i = variable[i]

        if var.is_integer
        else
            if var.is_fixed
                @inbounds bnd = var.lower_bound
                ci1 = MOI.add_constraint(opt, variable_i, ET(bnd))
                push!(x._lower_variable_et, ci1)
                push!(x._lower_variable_et_indx, i)
            else
                if var.has_lower_bound
                    @inbounds bnd = var.lower_bound
                    ci4 = MOI.add_constraint(opt, variable_i, GT(bnd))
                    push!(x._lower_variable_gt, ci4)
                    push!(x._lower_variable_gt_indx, i)
                end
                if var.has_upper_bound
                    @inbounds bnd = var.upper_bound
                    ci5 = MOI.add_constraint(opt, variable_i, LT(bnd))
                    push!(x._lower_variable_lt, ci5)
                    push!(x._lower_variable_lt_indx, i)
                end
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
    evaluator = x._relaxed_evaluator
    constraint_bounds = x._relaxed_constraint_bounds
    temp_la = CI{SAF,LT}[]
    temp_ua = CI{SAF,LT}[]
    for (j, bns) in enumerate(constraint_bounds)
        sparsity = grad_sparsity(evaluator, j+1)
        @inbounds nzvar = x._lower_variable_index[sparsity]
        func = SAF(SAT.(zeros(length(sparsity)), nzvar), 0.0)
        if !(bns.upper == Inf)
            ci = MOI.add_constraint(opt, func, LT(0.0))
            push!(temp_la, ci)
            push!(x._lower_nlp_affine_indx, j)
            push!(x._lower_nlp_sparsity, sparsity)
        end
        if !(bns.lower == -Inf)
            ci = MOI.add_constraint(opt, func, LT(0.0))
            push!(temp_ua, ci)
            push!(x._upper_nlp_affine_indx, j)
            push!(x._upper_nlp_sparsity, sparsity)
        end
    end

    len_la = length(temp_la)
    len_ua = length(temp_ua)
    push!(x._lower_nlp_affine, temp_la)
    push!(x._upper_nlp_affine, temp_ua)
    for i in 2:x.cut_max_iterations
        push!(x._lower_nlp_affine, fill(CI{SAF,LT}(-1), (len_la,)))
        push!(x._upper_nlp_affine, fill(CI{SAF,LT}(-1), (len_ua,)))
    end

    # only solves box constrained problems (+ other constraints) so single
    # variable objective implies that a bound has already been set and it
    # is not fixed
    if x._objective_is_sv
        for (i, z) in enumerate(x._lower_variable_lt)
            if x._objective_sv.variable.value == x._lower_variable_lt_indx[i]
                x._objective_cut_ci_sv = z
                break
            end
        end
    end

    for i in 1:x.cut_max_iterations
        push!(x._objective_cut_ci_saf, CI{SAF,LT}(-1))
    end
    MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    return
end

is_lp(m::Optimizer) = ~in(true, m.branch_variable)

function linear_solve!(m::Optimizer)

    opt = m.relaxed_optimizer

    # TODO: Add check for nonlinear terms which are actually linear
    if m._objective_is_sv
        MOI.set(opt, MOI.ObjectiveFunction{SV}(), m._objective_sv)
    elseif m._objective_is_saf
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

"""
    convert_to_min!

Converts MOI.MAX_SENSE objective to equivalent MOI.MIN_SENSE objective
max(f) = - min(-f).
"""
function convert_to_min!(x::Optimizer)
    if x._optimization_sense === MOI.MAX_SENSE
        if x._objective_is_sv
            x._objective_is_saf = true
            x._objective_is_sv = false
            x._objective_saf = SAF(SAT[SAT(-1.0, x._objective_sv.variable)], 0.0)
        elseif x._objective_is_saf
            @inbounds x._objective_saf.terms[:] = SAT.(-getfield.(x._objective_saf.terms, :coefficient),
                                                        getfield.(x._objective_saf.terms, :variable_index))
            x._objective_saf.constant *= -1.0
        elseif x._objective_is_sqf
            @inbounds x._objective_sqf.affine_terms[:] = SAT.(-getfield.(x._objective_sqf.affine_terms, :coefficient),
                                                               getfield.(x._objective_sqf.affine_terms, :variable_index))
            @inbounds x._objective_sqf.quadratic_terms[:] = SQT.(-getfield.(x._objective_sqf.quadratic_terms, :coefficient),
                                                                  getfield.(x._objective_sqf.quadratic_terms, :variable_index_1),
                                                                  getfield.(x._objective_sqf.quadratic_terms, :variable_index_2))
            x._objective_sqf.constant *= -1.0
        else
            nd = x._nlp_data.evaluator.m.nlp_data.nlobj.nd
            pushfirst!(nd, NodeData(JuMP._Derivatives.CALLUNIVAR, 2, -1))
            nd[2] = NodeData(nd[2].nodetype, nd[2].index, 1)
            for i in 3:length(nd)
                @inbounds nd[i] = NodeData(nd[i].nodetype, nd[i].index, nd[i].parent + 1)
            end
        end
    end
    return
end

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

# DONE
"""
    label_fixed_variables!

Detects any variables set to a fixed value by equality or inequality constraints
and populates the _fixed_variable storage array.
"""
function label_fixed_variables!(m::Optimizer)
    lbd = 0.0
    ubd = 0.0
    for i in 1:m._variable_number
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
function is_convex_quadratic(func::SQF, mult::Float64, cvx_dict::ImmutableDict{Int64,Int64})
    # Riguous Convexity Test
    flag = false
    row = Int64[]
    column = Int64[]
    value = Float64[]
    for term in func.quadratic_terms
        coeff = term.coefficient
        if coeff != 0.0
            value1 = cvx_dict[term.variable_index_1.value]
            value2 = cvx_dict[term.variable_index_2.value]
            mcoeff = mult*coeff
            push!(row, value1)
            push!(column, value2)
            push!(value, mcoeff)
            push!(row, value2)
            push!(column, value1)
            push!(value, mcoeff)
        end
    end
    Q = sparse(row, column, value)
    s1, s2 = size(Q)
    if length(Q.nzval) > 1
        eigval = eigmin(Array(Q))
        if (eigval) > 0.0
            flag = true
        end
    elseif (length(Q.nzval) == 1) && (Q.nzval[1] > 0.0)
        flag =  true
    else
        flag =  true
    end
    return flag
end


"""
    label_quadratic_convexity!

Assigns boolean value to constraint_convexity dictionary entry corresponding to
constraint index that is true if constraint is shown to be convex and false
otherwise.
"""
function label_quadratic_convexity!(x::Optimizer)

    for i in 1:length(x._quadratic_leq_constraints)
        @inbounds func, set, ind = x._quadratic_leq_constraints[i]
        @inbounds cvx_dict = x._quadratic_leq_dict[i]
        push!(x._quadratic_leq_convexity, is_convex_quadratic(func, 1.0, cvx_dict))
    end

    for i in 1:length(x._quadratic_geq_constraints)
        @inbounds func, set, ind = x._quadratic_geq_constraints[i]
        @inbounds cvx_dict = x._quadratic_geq_dict[i]
        push!(x._quadratic_geq_convexity, is_convex_quadratic(func, -1.0, cvx_dict))
    end

    for i in 1:length(x._quadratic_eq_constraints)
        @inbounds func, set, ind = x._quadratic_eq_constraints[i]
        @inbounds cvx_dict = x._quadratic_eq_dict[i]
        push!(x._quadratic_eq_convexity_1, is_convex_quadratic(func, 1.0, cvx_dict))
        push!(x._quadratic_eq_convexity_2, is_convex_quadratic(func, -1.0, cvx_dict))
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
        x._continuous_solution .= x._upper_solution
        x._objective_value = x._upper_objective_value
        x._termination_status_code = MOI.LOCALLY_SOLVED
        x._result_status_code = MOI.FEASIBLE_POINT
    else
        x._feasible_solution_found = false
        x._first_solution_node = x._maximum_node_id
        x._termination_status_code = MOI.LOCALLY_INFEASIBLE
        x._result_status_code = MOI.INFEASIBLE_POINT
    end
    return
end

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
    build_nlp_evaluator

Builds the evaluator used to generate relaxations of the nonlinear equations
and constraints from a source model.
"""
function build_nlp_evaluator(N::Int64, s::T, src::JuMP.NLPEvaluator, x::Optimizer, bool_flag::Bool) where {T<:RelaxTag}

    # Creates the empty evaluator
    d::Evaluator{N,T} = Evaluator{N,T}()
    build_nlp_kernel!(d, src, x, bool_flag)

    return d
end

function parse_problem!(m::Optimizer)

    ########### Reformulate DAG using auxilliary variables ###########
    _variable_len = length(m._variable_info)
    m._continuous_variable_number = _variable_len
    m._variable_number = _variable_len

    ########### Set Correct Size for Problem Storage #########
    m._current_xref = fill(0.0, _variable_len)
    m._cut_solution = fill(0.0, _variable_len)
    m._lower_solution = fill(0.0, _variable_len)
    m._upper_solution = fill(0.0, _variable_len)
    m._lower_lvd = fill(0.0, _variable_len)
    m._lower_uvd = fill(0.0, _variable_len)

    # Get various other sizes
    m._continuous_solution = zeros(Float64, _variable_len)

    m.presolve_flatten_flag && Script.dag_flattening!(m)
    convert_to_min!(m)
    initialize_evaluators!(m, false)               # initializes the EAGO and JuMP NLP evaluators

    new_time = time() - m._start_time
    m._parse_time = new_time
    m._run_time = new_time

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
    if no_constraints
        m.cp_depth = -1
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
        m.quad_uni_depth = -1
        m.quad_bi_depth = -1
    end

    only_lin_constraints = no_constraints
    no_lin_constraints = true
    if length(m._linear_leq_constraints) > 0
        no_lin_constraints &= false
    end
    if length(m._linear_geq_constraints) > 0
        no_lin_constraints &= false
    end
    if length(m._linear_eq_constraints) > 0
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

function presolve_problem!(m::Optimizer)

    m.presolve_epigraph_flag && reform_epigraph!(m)  # perform epigraph rearrangement
    m.presolve_cse_flag && dag_cse_simplify!(m)      #

    #m = user_reformed_optimizer(m)
    create_initial_node!(m)                        # Create initial node and add it to the stack
    label_fixed_variables!(m)
    label_nonlinear_variables!(m)
    check_disable_fbbt!(m)
    load_relaxed_problem!(m)
    label_quadratic_convexity!(m)

    m._presolve_time = time() - m._parse_time

    return
end

function store_candidate_solution!(x::Optimizer)
    if x._upper_feasibility
        if x._upper_objective_value < x._global_upper_bound
            x._feasible_solution_found = true
            x._first_solution_node = x._maximum_node_id
            x._solution_value = x._upper_objective_value
            x._global_upper_bound = x._upper_objective_value
            @inbounds x._continuous_solution[:] = x._upper_solution
        end
    end
    return
end

function set_global_lower_bound!(x::Optimizer)
    if ~isempty(x._stack)
        min_node = minimum(x._stack)
        lower_bound = min_node.lower_bound
        if x._global_lower_bound < lower_bound
            x._global_lower_bound = lower_bound
        end
    end
    return
end

# wraps subroutine call to isolate ExtensionType
termination_check(x::Optimizer) = termination_check(x.ext_type, x)
cut_condition(x::Optimizer) = cut_condition(x.ext_type, x)
convergence_check(x::Optimizer) = convergence_check(x.ext_type, x)
repeat_check(x::Optimizer) = repeat_check(x.ext_type, x)
node_selection!(x::Optimizer) = node_selection!(x.ext_type, x)
preprocess!(x::Optimizer) = preprocess!(x.ext_type, x)
lower_problem!(x::Optimizer) = lower_problem!(x.ext_type, x)
add_cut!(x::Optimizer) = add_cut!(x.ext_type, x)
upper_problem!(x::Optimizer) = upper_problem!(x.ext_type, x)
postprocess!(x::Optimizer) = postprocess!(x.ext_type, x)
single_storage!(x::Optimizer) = single_storage!(x.ext_type, x)
branch_node!(x::Optimizer) = branch_node!(x.ext_type, x)
fathom!(x::Optimizer) = fathom!(x.ext_type, x)
"""
    global_solve!

Solves the branch and bound problem with the input EAGO optimizer object.
"""
function global_solve!(x::Optimizer)

    x._iteration_count = 1
    x._node_count = 1

    # terminates when max nodes or iteration is reach, or when node stack is empty
    while ~termination_check(x)

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
                    if x._optimization_sense == MOI.FEASIBILITY_SENSE
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

function throw_optimize_hook!(m::Optimizer)
    optimize_hook!(m.ext_type, m)
end

function MOI.optimize!(m::Optimizer)

    m._start_time = time()
    parse_problem!(m)
    presolve_problem!(m)

    # Runs the branch and bound routine
    if ~m.enable_optimize_hook
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

    return
end
