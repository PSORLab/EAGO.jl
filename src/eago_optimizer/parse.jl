# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# Defines functions used to parse the input optimization problem into
# a solvable form including routines used to classify input problems as a
# LP, SOCP, MILP, MISOCP, and convex problem types.
#############################################################################

function reform_epigraph_min!(m::ParsedProblem, f::SV)
    ip = m._input_problem
    wp = m._working_problem
    flag = ip._optimization_sense == MOI.MAX_SENSE
    wp._objective = AffineFunctionIneq(f, is_max = flag)
    return
end
function reform_epigraph_min!(m::ParsedProblem, f::AffineFunctionIneq)
    ip = m._input_problem
    wp = m._working_problem
    if ip._optimization_sense == MOI.MAX_SENSE
        wp._objective_saf = MOIU.operate(-, Float64, wp._objective_saf)
        wp._objective = AffineFunctionIneq(f, GT_ZERO)
    end
    return
end
function convert_to_min!(m::ParsedProblem, f::BufferedQuadraticIneq)
    wp = m._working_problem
    ηi = m._variable_count + 1


    f.buffer[ηi] = 0.0                            # update buffered function
    f.len += 1
    MOIU.operate!(+, Float64, f.sqf, SV(VI(ηi)))
    if ip._optimization_sense == MOI.MAX_SENSE
        MOIU.operate!(-, Float64, f.sqf)
    end

    return
end
function convert_to_min!(m::ParsedProblem, f::BufferedNonlinearFunction)
    # updates tape for nlp_data block (used by local optimizer)
    wp = m._working_problem
    nd = wp._nlp_data.evaluator.m.nlp_data.nlobj.nd
    pushfirst!(nd, NodeData(JuMP._Derivatives.CALLUNIVAR, 2, -1))
    nd[2] = NodeData(nd[2].nodetype, nd[2].index, 1)
    for i = 3:length(nd)
        @inbounds nd[i] = NodeData(nd[i].nodetype, nd[i].index, nd[i].parent + 1)
    end

    # updates tape used by evaluator for the nonlinear objective (used by the relaxed optimizer)
    _negate!(wp._objective_nl.ex.g)
    pushfirst!(_set(wp._objective_nl.ex.relax_cache), _set(wp._objective_nl))
    pushfirst!(_num(wp._objective_nl.ex.relax_cache), 0.0)
    pushfirst!(_is_num(wp._objective_nl.ex.relax_cache), false)
    return
end

"""
$(TYPEDSIGNATURES)

Converts `MOI.MAX_SENSE` objective to equivalent `MOI.MIN_SENSE` objective
`max(f) = -min(-f)`.
"""
function convert_to_min!(m::Optimizer)
    wp = m._working_problem
    if m._input_problem._optimization_sense == MOI.MAX_SENSE
        convert_to_min!(wp, wp._objective)
    end
    return
end

"""

Performs an epigraph reformulation assuming the working_problem is a minimization problem.
"""
function reform_epigraph_min!(m::Optimizer)

    # add epigraph variable
    obj_variable_index = MOI.add_variable(m)

    # converts ax + b objective to ax - y <= -b constraint with y objective
    obj_type = m._working_problem._objective_type
    if obj_type == SCALAR_AFFINE

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
and populates the `_fixed_variable` storage array.
"""
function label_fixed_variables!(m::Optimizer)
    map!(x -> check_set_is_fixed(x), m._fixed_variable, m._working_problem._variable_info)
end

"""
$(TYPEDSIGNATURES)

Detects any variables participating in nonconvex terms and populates the
`_branch_variables` storage array.
"""
function label_branch_variables!(m::Optimizer)

    wp = m._working_problem
    m._user_branch_variables = !isempty(m._parameters.branch_variable)
    if m._user_branch_variables
        m._branch_variables = m._parameters.branch_variable
    else
        m._branch_variables = fill(false, m._working_problem._variable_count)
        for f in wp._sqf_leq, t in f.func.quadratic_terms
            m._branch_variables[t.variable_index_1.value] = true
            m._branch_variables[t.variable_index_2.value] = true
        end
        for f in wp._sqf_eq, t in f.func.quadratic_terms
            m._branch_variables[t.variable_index_1.value] = true
            m._branch_variables[t.variable_index_2.value] = true
        end
        for f in wp._nonlinear_constr, i in _sparsity(f)
            m._branch_variables[i] = true
        end
        if wp._objective isa BufferedQuadraticIneq
            for t in wp._objective.func.quadratic_terms
                m._branch_variables[t.variable_index_1.value] = true
                m._branch_variables[t.variable_index_2.value] = true
            end
        elseif wp._objective isa BufferedNonlinearFunction
            for i in _sparsity(wp._objective)
                m._branch_variables[i] = true
            end
        end
    end

    # add a map of branch/node index to variables in the continuous solution
    for i = 1:wp._variable_count
        if is_fixed(wp._variable_info[i])
            m._branch_variables[i] = false
        elseif m._branch_variables[i]
            push!(m._branch_to_sol_map, i)
        end
    end

    # creates reverse map
    m._sol_to_branch_map = zeros(wp._variable_count)
    for i = 1:length(m._branch_to_sol_map)
        j = m._branch_to_sol_map[i]
        m._sol_to_branch_map[j] = i
    end

    # adds branch solution to branch map to evaluator
    vnum = wp._variable_count
    initialize!(m._branch_cost, length(m._branch_to_sol_map))
    v = VariableValues{Float64}(x = zeros(vnum),
                                lower_variable_bounds = zeros(vnum),
                                upper_variable_bounds = zeros(vnum),
                                node_to_variable_map = m._branch_to_sol_map,
                                variable_to_node_map = m._sol_to_branch_map)
    wp._relaxed_evaluator.variable_values = v
    # sets variable_values in nonlinear expressions/functions to v
    # add nonlinear objective
    (wp._objective isa BufferedNonlinearFunction) && _set_variable_storage!(wp._objective, v)
    foreach(i -> _set_variable_storage!(i, v), wp._nonlinear_constr)
    foreach(i -> _set_variable_storage!(i, v), wp._relaxed_evaluator.subexpressions)
    return
end

add_nonlinear_evaluator!(m::Optimizer, evaluator::Nothing) = nothing
add_nonlinear_evaluator!(m::Optimizer, evaluator::EmptyNLPEvaluator) = nothing
function add_nonlinear_evaluator!(m::Optimizer, d::JuMP.NLPEvaluator)
    m._working_problem._relaxed_evaluator = Evaluator()

    variable_count = length(m._working_problem._variable_info)
    relax_eval = m._working_problem._relaxed_evaluator
    relax_eval.user_operators = OperatorRegistry(d.m.nlp_data.user_operators)
    relax_eval.ctx            = GuardCtx(metadata = GuardTracker(m._parameters.domain_violation_ϵ,
                                                                 m._parameters.domain_violation_guard_on))
    relax_eval.subgrad_tol           = m._parameters.subgrad_tol
    m._nonlinear_evaluator_created = true
    return
end
function add_nonlinear_evaluator!(m::Optimizer)
    add_nonlinear_evaluator!(m, m._input_problem._nlp_data.evaluator)
end

add_nonlinear!(m::Optimizer, evaluator::Nothing) = nothing
add_nonlinear!(m::Optimizer, evaluator::EmptyNLPEvaluator) = nothing
function add_nonlinear!(m::Optimizer, evaluator::JuMP.NLPEvaluator)

    add_nonlinear_evaluator!(m, m._input_problem._nlp_data.evaluator)

    nlp_data = m._input_problem._nlp_data
    MOI.initialize(evaluator, Symbol[:Grad, :ExprGraph])
    user_operator_registry = OperatorRegistry(evaluator.m.nlp_data.user_operators)

    # set nlp data structure
    m._working_problem._nlp_data = nlp_data

    # add subexpressions (assumes they are already ordered by JuMP)
    # creates a dictionary that lists the subexpression sparsity
    # by search each node for variables dict[2] = [2,3] indicates
    # that subexpression 2 depends on variables 2 and 3
    # this is referenced when subexpressions are called by other
    # subexpressions or functions to determine overall sparsity
    # the sparsity of a function is the collection of indices
    # in all participating subexpressions and the function itself
    # it is necessary to define this as such to enable reverse
    # McCormick constraint propagation
    relax_evaluator = m._working_problem._relaxed_evaluator
    dict_sparsity = Dict{Int,Vector{Int}}()
    if length(evaluator.m.nlp_data.nlexpr) > 0      # should check for nonlinear objective, constraint
        for i = 1:length(evaluator.subexpressions)
            subexpr = evaluator.subexpressions[i]
            push!(relax_evaluator.subexpressions, NonlinearExpression!(subexpr, dict_sparsity, i,
                                                                      evaluator.subexpression_linearity,
                                                                      user_operator_registry,
                                                                      m._parameters.relax_tag,
                                                                      evaluator.m.nlp_data.nlparamvalues))
        end
    end

    # scrubs udf functions using Cassette to remove odd data structures...
    # alternatively convert udfs to JuMP scripts...
    m._parameters.presolve_scrubber_flag && Script.scrub!(m._working_problem._nlp_data)
    if m._parameters.presolve_to_JuMP_flag
        Script.udf_loader!(m)
    end

    parameter_values = copy(evaluator.parameter_values)

    # add nonlinear objective
    if evaluator.has_nlobj
        m._working_problem._objective = BufferedNonlinearFunction(evaluator.objective, MOI.NLPBoundsPair(-Inf, Inf),
                                                                 dict_sparsity, evaluator.subexpression_linearity,
                                                                 user_operator_registry,
                                                                 evaluator.m.nlp_data.nlparamvalues,
                                                                 m._parameters.relax_tag)
    end

    # add nonlinear constraints
    constraint_bounds = m._working_problem._nlp_data.constraint_bounds
    for i = 1:length(evaluator.constraints)
        constraint = evaluator.constraints[i]
        bnds = constraint_bounds[i]
        push!(m._working_problem._nonlinear_constr, BufferedNonlinearFunction(constraint, bnds, dict_sparsity,
                                                                              evaluator.subexpression_linearity,
                                                                              user_operator_registry,
                                                                              evaluator.m.nlp_data.nlparamvalues,
                                                                              m._parameters.relax_tag))
    end

    m._input_problem._nonlinear_count = length(m._working_problem._nonlinear_constr)
    m._working_problem._nonlinear_count = length(m._working_problem._nonlinear_constr)
    relax_evaluator.subexpressions_eval = fill(false, length(relax_evaluator.subexpressions))
    return
end
add_nonlinear!(m::Optimizer) = add_nonlinear_functions!(m, m._input_problem._nlp_data.evaluator)

add_objective!(m::ParsedProblem, f) = nothing
add_objective!(m::ParsedProblem, f::SV) = (m._objective = f; return )
function add_objective!(m::ParsedProblem, f::SAF)
    m._objective = AffineFunctionIneq(f, LT_ZERO)
    return
end
function add_objective!(m::ParsedProblem, f::SQF)
    m._objective = BufferedQuadraticIneq(f, LT_ZERO)
end

"""
Translates input problem to working problem. Routines and checks and optional manipulation is left to the presolve stage.
"""
function initial_parse!(m::Optimizer)

    # reset initial time and solution statistics
    m._time_left = m._parameters.time_limit

    # add variables to working model
    ip = m._input_problem
    wp = m._working_problem
    append!(wp._variable_info, ip._variable_info)
    wp._variable_count = ip._variable_count

    # add linear constraints to the working problem
    append!(wp._saf_leq, [AffineFunctionIneq(c[1], c[2]) for c in ip._linear_leq_constraints])
    append!(wp._sqf_leq, [AffineFunctionIneq(c[1], c[2]) for c in ip._linear_geq_constraints])
    wp._saf_eq  = [AffineFunctionEq(c[1], c[2]) for c in ip._linear_eq_constraints]

    append!(wp._sqf_leq, [BufferedQuadraticIneq(c[1], c[2]) for c in ip._quadratic_leq_constraints])
    append!(wp._sqf_leq, [BufferedQuadraticIneq(c[1], c[2]) for c in ip._quadratic_geq_constraints])
    wp._sqf_eq  = [BufferedQuadraticEq(c[1], c[2]) for c in ip._quadratic_eq_constraints]

    # set objective function
    add_objective!(m._working_problem, ip._objective)

    # add nonlinear constraints, evaluator, subexpressions
    add_nonlinear!(m)

    # converts a maximum problem to a minimum problem (internally) if necessary
    # this is placed after adding nonlinear functions as this prior routine
    # copies the nlp_block from the input_problem to the working problem
    reform_epigraph_min!(m)

    label_fixed_variables!(m)
    label_branch_variables!(m)

    # updates run and parse times
    new_time = time() - m._start_time
    m._parse_time = new_time
    m._run_time = new_time

    return
end

"""
Classifies the problem type
"""
function parse_classify_problem!(m::Optimizer)

    ip = m._input_problem
    integer_variable_number = count(is_integer.(ip._variable_info))

    nl_expr_number = ip._objective == nothing ? 1 : 0
    nl_expr_number += length(m._working_problem._nonlinear_constr)
    cone_constraint_number = length(ip._conic_second_order)
    quad_constraint_number = length(ip._quadratic_leq_constraints) +
                             length(ip._quadratic_geq_constraints) +
                             length(ip._quadratic_eq_constraints)

    linear_or_sv_objective = (ip._objective isa SV || ip._objective isa SAF)
    relaxed_supports_soc = false

    if cone_constraint_number === 0 && quad_constraint_number === 0 &&
        nl_expr_number === 0 && linear_or_sv_objective
        m._working_problem._problem_type = LP
    elseif quad_constraint_number === 0 && relaxed_supports_soc &&
        nl_expr_number === 0 && linear_or_sv_objective
        m._working_problem._problem_type = SOCP
    else
        m._working_problem._problem_type = MINCVX
    end
    return
end

"""

Basic parsing for global solutions (no extensive manipulation)
"""
parse_global!(t::ExtensionType, m::Optimizer) = nothing
parse_global!(m::Optimizer) = parse_global!(m.ext_type, m)
