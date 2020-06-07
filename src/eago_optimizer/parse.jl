# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# Defines functions used to parse the input optimization problem into
# a solvable form including routines used to classify input problems as a
# LP, SOCP, MILP, MISOCP, and convex problem types.
#############################################################################

"""
$(TYPEDSIGNATURES)

Converts MOI.MAX_SENSE objective to equivalent MOI.MIN_SENSE objective
max(f) = - min(-f).
"""
function convert_to_min!(m::Optimizer)

    m._working_problem._optimization_sense = MOI.MIN_SENSE
    if m._input_problem._optimization_sense === MOI.MAX_SENSE

        obj_type = m._input_problem._objective_type
        if obj_type === SINGLE_VARIABLE
            m._working_problem._objective_type = SCALAR_AFFINE
            m._working_problem._objective_saf = MOIU.operate(-, Float64, m._working_problem._objective_sv)
            m._working_problem._objective_saf_parsed = AffineFunctionIneq(m._working_problem._objective_saf, LT_ZERO)

        elseif obj_type === SCALAR_AFFINE
            m._working_problem._objective_saf = MOIU.operate(-, Float64, m._working_problem._objective_saf)
            m._working_problem._objective_saf_parsed = AffineFunctionIneq(m._working_problem._objective_saf, LT_ZERO)

        elseif obj_type === SCALAR_QUADRATIC
            sqf = m._working_problem._objective_sqf.sqf
            m._working_problem._objective_sqf.sqf = MOIU.operate(-, Float64, sqf)

        elseif obj_type === NONLINEAR

            # updates tape for nlp_data block (used by local optimizer)
            nd = m._working_problem._nlp_data.evaluator.m.nlp_data.nlobj.nd
            pushfirst!(nd, NodeData(JuMP._Derivatives.CALLUNIVAR, 2, -1))
            nd[2] = NodeData(nd[2].nodetype, nd[2].index, 1)
            for i = 3:length(nd)
                @inbounds nd[i] = NodeData(nd[i].nodetype, nd[i].index, nd[i].parent + 1)
            end

            # updates tape used by evaluator for the nonlinear objective (used by the relaxed optimizer)
            nd = m._working_problem._objective_nl.expr.nd
            pushfirst!(nd, NodeData(JuMP._Derivatives.CALLUNIVAR, 2, -1))
            nd[2] = NodeData(nd[2].nodetype, nd[2].index, 1)
            for i = 3:length(nd)
                @inbounds nd[i] = NodeData(nd[i].nodetype, nd[i].index, nd[i].parent + 1)
            end
        end
    end

    return nothing
end

"""

Performs an epigraph reformulation assuming the working_problem is a minimization problem.
"""
function reform_epigraph!(m::Optimizer)

    if m._parameters.presolve_epigraph_flag
        #=
        # add epigraph variable
        obj_variable_index = MOI.add_variable(m)

        # converts ax + b objective to ax - y <= -b constraint with y objective
        obj_type = m._working_problem._objective_type
        if obj_type === SCALAR_AFFINE

            # update unparsed expression
            objective_saf = m._working_problem._objective_saf
            push!(objective_saf, SAT(-1.0, obj_variable_index))
            obj_ci = MOI.add_constraint(m, saf, LT(-objective_saf.constant))

            # update parsed expression (needed for interval bounds)

        # converts ax + b objective to ax - y <= -b constraint with y objective
        elseif obj_type === SCALAR_QUADRATIC

            # update parsed expression
            objective_sqf = m._working_problem._objective_sqf
            obj_ci = MOI.add_constraint(m, saf, LT())

        elseif obj_type === NONLINEAR

            # updated parsed expressions
            objective_nl = m._working_problem._objective_nl

        end

        MOI.set(m, MOI.ObjectiveFunction{SV}(), SV(obj_variable_index))
        =#
    end

    return nothing
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
    map!(x -> check_set_is_fixed(x), m._fixed_variable, m._working_problem._variable_info)
end

"""
$(TYPEDSIGNATURES)

Detects any variables participating in nonconvex terms and populates the
_branch_variables storage array.
"""
function label_branch_variables!(m::Optimizer)

    m._user_branch_variables = !isempty(m._branch_variables)
    m._user_branch_variables && (return nothing)

    append!(m._branch_variables, fill(false, m._working_problem._variable_count))

    # adds nonlinear terms in quadratic constraints
    sqf_leq = m._working_problem._sqf_leq
    for i = 1:m._working_problem._sqf_leq_count
        quad_ineq = @inbounds sqf_leq[i]
        for term in quad_ineq.func.quadratic_terms
            variable_index_1 = term.variable_index_1.value
            variable_index_2 = term.variable_index_2.value
            @inbounds m._branch_variables[variable_index_1] = true
            @inbounds m._branch_variables[variable_index_2] = true
        end
    end

    sqf_eq = m._working_problem._sqf_eq
    for i = 1:m._working_problem._sqf_eq_count
        quad_eq = @inbounds sqf_eq[i]
        for term in quad_eq.func.quadratic_terms
            variable_index_1 = term.variable_index_1.value
            variable_index_2 = term.variable_index_2.value
            @inbounds m._branch_variables[variable_index_1] = true
            @inbounds m._branch_variables[variable_index_2] = true
        end
    end

    obj_type = m._working_problem._objective_type
    if obj_type === SCALAR_QUADRATIC
        for term in m._working_problem._objective_sqf.func.quadratic_terms
            variable_index_1 = term.variable_index_1.value
            variable_index_2 = term.variable_index_2.value
            @inbounds m._branch_variables[variable_index_1] = true
            @inbounds m._branch_variables[variable_index_2] = true
        end
    end

    # label nonlinear branch variables (assumes affine terms have been extracted)
    nl_constr = m._working_problem._nonlinear_constr
    for i = 1:m._working_problem._nonlinear_count
        nl_constr_eq = @inbounds nl_constr[i]
        node_list = nl_constr_eq.expr.nd
        for node in node_list
            if node.nodetype === JuMP._Derivatives.VARIABLE
                @inbounds m._branch_variables[node.index] = true
            end
        end
    end

    if obj_type === NONLINEAR
        node_list = m._working_problem._objective_nl.expr.nd
        for node in node_list
            if node.nodetype === JuMP._Derivatives.VARIABLE
                @inbounds m._branch_variables[node.index] = true
            end
        end
    end

    # add a map of branch/node index to variables in the continuous solution
    for i = 1:m._working_problem._variable_count
        if m._working_problem._variable_info[i].is_fixed
            m._branch_variables[i] = false
            continue
        end
        if m._branch_variables[i]
            push!(m._branch_to_sol_map, i)
        end
    end

    # creates reverse map
    m._sol_to_branch_map = zeros(m._working_problem._variable_count)
    for i = 1:length(m._branch_to_sol_map)
        j = m._branch_to_sol_map[i]
        m._sol_to_branch_map[j] = i
    end

    # adds branch solution to branch map to evaluator
    m._working_problem._relaxed_evaluator.node_to_variable_map = m._branch_to_sol_map
    m._working_problem._relaxed_evaluator.variable_to_node_map = m._sol_to_branch_map
    m._working_problem._relaxed_evaluator.node_count = length(m._branch_to_sol_map)

    return nothing
end


add_nonlinear_functions!(m::Optimizer) = add_nonlinear_functions!(m, m._input_problem._nlp_data.evaluator)

add_nonlinear_functions!(m::Optimizer, evaluator::Nothing) = nothing
add_nonlinear_functions!(m::Optimizer, evaluator::EmptyNLPEvaluator) = nothing
function add_nonlinear_functions!(m::Optimizer, evaluator::JuMP.NLPEvaluator)

    nlp_data = m._input_problem._nlp_data
    MOI.initialize(evaluator, Symbol[:Grad, :ExprGraph])

    # set nlp data structure
    m._working_problem._nlp_data = nlp_data


    #println("typeof(nlp_data) = $(nlp_data)")

    # scrubs udf functions using Cassette to remove odd data structures...
    # alternatively convert udfs to JuMP scripts...
    m._parameters.presolve_scrubber_flag && Script.scrub!(m._working_problem._nlp_data)
    if m._parameters.presolve_to_JuMP_flag
        Script.udf_loader!(m)
    end

    parameter_values = copy(evaluator.parameter_values)

    # add nonlinear objective
    if evaluator.has_nlobj
        m._working_problem._objective_nl = BufferedNonlinearFunction(evaluator.objective, MOI.NLPBoundsPair(-Inf, Inf),
                                                                     evaluator.subexpression_linearity,
                                                                     m._parameters.relax_tag)
    end

    # add nonlinear constraints
    constraint_bounds = m._working_problem._nlp_data.constraint_bounds
    for i = 1:length(evaluator.constraints)
        constraint = evaluator.constraints[i]
        bnds = constraint_bounds[i]
        push!(m._working_problem._nonlinear_constr, BufferedNonlinearFunction(constraint, bnds,
                                                                              evaluator.subexpression_linearity,
                                                                              m._parameters.relax_tag))
    end
     m._working_problem._nonlinear_count =  length(m._working_problem._nonlinear_constr)
    return nothing
end

function add_nonlinear_evaluator!(m::Optimizer)
    evaluator = m._input_problem._nlp_data.evaluator
    add_nonlinear_evaluator!(m, evaluator)
    return nothing
end

add_nonlinear_evaluator!(m::Optimizer, evaluator::Nothing) = nothing
add_nonlinear_evaluator!(m::Optimizer, evaluator::EmptyNLPEvaluator) = nothing
function add_nonlinear_evaluator!(m::Optimizer, evaluator::JuMP.NLPEvaluator)
    m._working_problem._relaxed_evaluator = Evaluator()

    relax_evaluator = m._working_problem._relaxed_evaluator
    has_subexpressions = length(evaluator.m.nlp_data.nlexpr) > 0
    if has_subexpressions
        for i = 1:length(evaluator.subexpressions)
            subexpr = evaluator.subexpressions[i]
            push!(relax_evaluator._nonlinear_subexpr, BufferedNonlinearSubexpr(subexpr,
                                                                               evaluator.subexpression_linearity,
                                                                               m._parameters.relax_tag))
        end
    end

    println("length(m._working_problem._variable_info) = $(length(m._working_problem._variable_info))")
    relax_evaluator.variable_count = length(m._working_problem._variable_info)
    relax_evaluator.user_operators = evaluator.m.nlp_data.user_operators

    relax_evaluator.lower_variable_bounds = zeros(relax_evaluator.variable_count)
    relax_evaluator.upper_variable_bounds = zeros(relax_evaluator.variable_count)
    relax_evaluator.x                     = zeros(relax_evaluator.variable_count)
    relax_evaluator.num_mv_buffer         = zeros(relax_evaluator.variable_count)
    relax_evaluator.cv_grad_buffer        = zeros(relax_evaluator.variable_count)
    relax_evaluator.cc_grad_buffer        = zeros(relax_evaluator.variable_count)
    relax_evaluator.treat_x_as_number     = fill(false, relax_evaluator.variable_count)
    relax_evaluator.ctx       = GuardCtx(metadata = GuardTracker(m._parameters.domain_violation_ϵ))
    relax_evaluator.subgrad_tol       = m._parameters.subgrad_tol

    m._nonlinear_evaluator_created = true

    return nothing
end

"""
Translates input problem to working problem. Routines and checks and optional manipulation is left to the presolve stage.
"""
function initial_parse!(m::Optimizer)

    # reset initial time and solution statistics
    m._time_left = m._parameters.time_limit

    # add variables to working model
    ip = m._input_problem
    append!(m._working_problem._variable_info, ip._variable_info)
    m._working_problem._variable_count = ip._variable_count

    # add linear constraints to the working problem
    linear_leq = ip._linear_leq_constraints
    for i = 1:ip._linear_leq_count
        linear_func, leq_set = @inbounds linear_leq[i]
        push!(m._working_problem._saf_leq, AffineFunctionIneq(linear_func, leq_set))
        m._working_problem._saf_leq_count += 1
    end

    linear_geq = ip._linear_geq_constraints
    for i = 1:ip._linear_geq_count
        linear_func, geq_set = @inbounds linear_geq[i]
        push!(m._working_problem._saf_leq, AffineFunctionEq(linear_func, geq_set))
        m._working_problem._saf_leq_count += 1
    end

    linear_eq = ip._linear_eq_constraints
    for i = 1:ip._linear_eq_count
        linear_func, eq_set = @inbounds linear_eq[i]
        push!(m._working_problem._saf_eq, AffineFunctionEq(linear_func, eq_set))
        m._working_problem._saf_eq_count += 1
    end

    # add quadratic constraints to the working problem
    quad_leq = ip._quadratic_leq_constraints
    for i = 1:ip._quadratic_leq_count
        quad_func, leq_set = @inbounds quad_leq[i]
        push!(m._working_problem._sqf_leq, BufferedQuadraticIneq(quad_func, leq_set))
        m._working_problem._sqf_leq_count += 1
    end

    quad_geq = ip._quadratic_geq_constraints
    for i = 1:ip._quadratic_geq_count
        quad_func, geq_set = @inbounds quad_geq[i]
        push!(m._working_problem._sqf_leq, BufferedQuadraticIneq(quad_func, geq_set))
        m._working_problem._sqf_leq_count += 1
    end

    quad_eq = ip._quadratic_eq_constraints
    for i = 1:ip._quadratic_eq_count
        quad_func, eq_set = @inbounds quad_eq[i]
        push!(m._working_problem._sqf_eq, BufferedQuadraticEq(quad_func, eq_set))
        m._working_problem._sqf_eq_count += 1
    end

    # add conic constraints to the working problem
    soc_vec = m._input_problem._conic_second_order
    for i = 1:ip._conic_second_order_count
        soc_func, soc_set = @inbounds soc_vec[i]
        first_variable_loc = soc_func.variables[1].value
        prior_lbnd = m._working_problem._variable_info[first_variable_loc].lower_bound
        m._working_problem._variable_info[first_variable_loc].lower_bound = max(prior_lbnd, 0.0)
        push!(m._working_problem._conic_second_order, BufferedSOC(soc_func, soc_set))
        m._working_problem._conic_second_order_count += 1
    end

    # set objective function
    m._working_problem._objective_type = ip._objective_type
    m._working_problem._objective_sv = ip._objective_sv
    m._working_problem._objective_saf = ip._objective_saf
    m._working_problem._objective_saf_parsed = AffineFunctionIneq(ip._objective_saf, LT_ZERO)
    m._working_problem._objective_sqf = BufferedQuadraticIneq(ip._objective_sqf, LT_ZERO)

    # add nonlinear constraints
    # the nonlinear evaluator loads with populated subexpressions which are then used
    # to asssess the linearity of subexpressions
    add_nonlinear_evaluator!(m)
    add_nonlinear_functions!(m)

    # converts a maximum problem to a minimum problem (internally) if necessary
    # this is placed after adding nonlinear functions as this prior routine
    # copies the nlp_block from the input_problem to the working problem
    convert_to_min!(m)
    reform_epigraph!(m)

    # labels the variable info and the _fixed_variable vector for each fixed variable
    label_fixed_variables!(m)

    # labels variables to branch on
    label_branch_variables!(m)

    # updates run and parse times
    new_time = time() - m._start_time
    m._parse_time = new_time
    m._run_time = new_time

    return nothing
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

    nl_constraint_number = ip._nonlinear_leq_count + ip._nonlinear_eq_count
    cone_constraint_number = ip._conic_second_order_count
    quad_constraint_number = ip._quadratic_leq_count + ip._quadratic_geq_count + ip._quadratic_eq_count

    linear_or_sv_objective = (ip._objective_type === SINGLE_VARIABLE || ip._objective_type === SCALAR_AFFINE)
    relaxed_supports_soc = MOI.supports_constraint(m.relaxed_optimizer, VECOFVAR, SOC)

    if integer_variable_number === 0

        if cone_constraint_number === 0 && quad_constraint_number === 0 &&
            nl_constraint_number === 0 && linear_or_sv_objective
            m._working_problem._problem_type = LP

            println(" was lp")
        elseif quad_constraint_number === 0 && relaxed_supports_soc &&
               nl_constraint_number === 0 && linear_or_sv_objective
            m._working_problem._problem_type = SOCP
            println(" was socp")
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
            m._working_problem._problem_type = MINCVX
            println(" was minlp")

        end
    else
        #=
        if cone_constraint_number === 0 && quad_constraint_number === 0 && linear_or_sv_objective
        elseif quad_constraint_number === 0 && relaxed_supports_soc && linear_or_sv_objective
            m._working_problem._problem_type = MISOCP
        else
            #parse_classify_quadratic!(m)
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
        =#
    end

    return nothing
end

"""

Basic parsing for global solutions (no extensive manipulation)
"""
function parse_global!(t::ExtensionType, m::Optimizer)
end
