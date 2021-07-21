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

"""

Adds objective function (if any) to the parsed problem.
"""
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

Adds a Evaluator structure if nonlinear terms are attached.
"""
add_nonlinear_evaluator!(m::GlobalOptimizer, evaluator::Nothing) = nothing
function add_nonlinear_evaluator!(m::GlobalOptimizer, d::JuMP.NLPEvaluator)
    m._working_problem._relaxed_evaluator = Evaluator()
    relax_eval = m._working_problem._relaxed_evaluator
    relax_eval.user_operators = OperatorRegistry(d.m.nlp_data.user_operators)
    relax_eval.subgrad_tol    = m._parameters.subgrad_tol
    m._nonlinear_evaluator_created = true
    return
end
add_nonlinear_evaluator!(m::GlobalOptimizer, nldata) = add_nonlinear_evaluator!(m, nldata.evaluator)
add_nonlinear_evaluator!(m::GlobalOptimizer) = add_nonlinear_evaluator!(m, m._input_problem._nlp_data)

"""

Adds a Evaluator, nonlinear functions, and populates each appropriately.
"""
add_nonlinear!(m::GlobalOptimizer, evaluator::Nothing) = nothing
function add_nonlinear!(m::GlobalOptimizer, evaluator::JuMP.NLPEvaluator)

    add_nonlinear_evaluator!(m, m._input_problem._nlp_data.evaluator)

    nlp_data = m._input_problem._nlp_data
    MOI.initialize(evaluator, Symbol[:Grad, :Jac, :ExprGraph])
    user_operator_registry = OperatorRegistry(evaluator.m.nlp_data.user_operators)

    # set nlp data structure
    m._working_problem._nlp_data = nlp_data
    mul_relax = m._parameters.mul_relax_style
    rtype = (mul_relax == 1) ? ((mul_relax_style == 2) ? RelaxAA() : RelaxMulEnum()) : Relax()
    renum = (mul_relax == 1) ? ((mul_relax_style == 2) ? MC_AFF_RELAX : MC_ENUM_RELAX) : STD_RELAX
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
    relax_evaluator.relax_type = renum
    dict_sparsity = Dict{Int,Vector{Int}}()
    if length(evaluator.m.nlp_data.nlexpr) > 0      # should check for nonlinear objective, constraint
        for i = 1:length(evaluator.subexpressions)
            subexpr = evaluator.subexpressions[i]
            nlexpr = NonlinearExpression!(rtype, subexpr, dict_sparsity, i, evaluator.subexpression_linearity, 
                                          user_operator_registry, m._parameters.relax_tag, evaluator.m.nlp_data.nlparamvalues)
            push!(relax_evaluator.subexpressions, nlexpr)
        end
    end

    # scrubs udf functions using Cassette to remove odd data structures...
    # alternatively convert udfs to JuMP scripts...
    m._parameters.presolve_scrubber_flag && Script.scrub!(m._working_problem._nlp_data)
    if m._parameters.presolve_to_JuMP_flag
        Script.udf_loader!(m)
    end

    # add nonlinear objective
    if evaluator.has_nlobj
        m._working_problem._objective = BufferedNonlinearFunction(rtype, evaluator.objective, MOI.NLPBoundsPair(-Inf, Inf),
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
        push!(m._working_problem._nonlinear_constr, BufferedNonlinearFunction(rtype, constraint, bnds, dict_sparsity,
                                                                              evaluator.subexpression_linearity,
                                                                              user_operator_registry,
                                                                              evaluator.m.nlp_data.nlparamvalues,
                                                                              m._parameters.relax_tag))
    end

    relax_evaluator.subexpressions_eval = fill(false, length(relax_evaluator.subexpressions))
    return
end
add_nonlinear!(m::GlobalOptimizer, nldata) = add_nonlinear!(m, nldata.evaluator)
add_nonlinear!(m::GlobalOptimizer) = add_nonlinear!(m, m._input_problem._nlp_data)

function reform_epigraph_min!(m::GlobalOptimizer, d::ParsedProblem, f::SV)
    flag = m._input_problem._optimization_sense == MOI.MAX_SENSE
    d._objective = AffineFunctionIneq(f, is_max = flag)
    d._objective_saf = SAF([SAT(flag ? -1.0 : 1.0, VI(f.variable.value))], 0.0)
    return
end
function reform_epigraph_min!(m::GlobalOptimizer, d::ParsedProblem, f::AffineFunctionIneq)
    d._objective_saf = m._input_problem._objective
    if m._input_problem._optimization_sense == MOI.MAX_SENSE
        MOIU.operate!(-, Float64, d._objective_saf)
        d._objective = AffineFunctionIneq([(-i,j) for (i,j) in f.terms], -f.constant, f.len)
    else
        d._objective = f
    end
    return
end

function add_η!(m::ParsedProblem, l::Float64, u::Float64)
    m._variable_count += 1
    push!(m._variable_info, VariableInfo(MOI.Interval{Float64}(l, u)))
    return m._variable_count
end

reform_epigraph_min!(m::GlobalOptimizer, d::ParsedProblem, f::Nothing) = nothing
function reform_epigraph_min!(m::GlobalOptimizer, d::ParsedProblem, f::BufferedQuadraticIneq)
    ip = m._input_problem

    vi = d._variable_info
    n = NodeBB(lower_bound.(vi), upper_bound.(vi), is_integer.(vi))
    m._current_node = n

    l, u = interval_bound(m, f)
    ηi = add_η!(d, l, u)
    m._global_lower_bound = l
    m._global_upper_bound = u
    d._objective_saf = SAF([SAT(1.0, VI(ηi))], 0.0)

    f.buffer[ηi] = 0.0
    f.len += 1
    if ip._optimization_sense == MOI.MAX_SENSE
        MOIU.operate!(-, Float64, f.func)
    end
    MOIU.operate!(-, Float64, f.func, SV(VI(ηi)))
    push!(f.saf.terms, SAT(0.0, VI(ηi)))
    m._obj_var_slack_added = true

    return
end
function reform_epigraph_min!(m::GlobalOptimizer, d::ParsedProblem, f::BufferedNonlinearFunction)
    ip = m._input_problem
    wp = m._working_problem

    vi = ip._variable_info
    vi_mid = mid.(vi)
    vi_lo = lower_bound.(vi)
    vi_hi = upper_bound.(vi)

    q = _variable_num(FullVar(), m)
    v = VariableValues{Float64}(x = vi_mid,
                                lower_variable_bounds = vi_lo,
                                upper_variable_bounds = vi_hi,
                                node_to_variable_map = [i for i in 1:q],
                                variable_to_node_map = [i for i in 1:q])
    wp._relaxed_evaluator.variable_values = v
    _set_variable_storage!(f,v)

    n = NodeBB(vi_lo, vi_hi, is_integer.(vi))
    m._current_node = n
    set_node!(wp._relaxed_evaluator, n)

    l, u = interval_bound(m, f)
    if !_is_input_min(m)
        l, u = -u, -l
    end
    ηi = add_η!(d, l, u)
    @variable(ip._nlp_data.evaluator.m, l <= η <= u)
    m._global_lower_bound = l
    m._global_upper_bound = u
    wp._objective_saf = SAF([SAT(1.0, VI(ηi))], 0.0)

    nd = ip._nlp_data.evaluator.m.nlp_data.nlobj.nd
    if !_is_input_min(m)
        pushfirst!(nd, NodeData(JuMP._Derivatives.CALLUNIVAR, 2, 1))
        pushfirst!(nd, NodeData(JuMP._Derivatives.CALL, 2, -1))
        nd[3] = NodeData(nd[3].nodetype, nd[3].index, 2)
        for i = 4:length(nd)
            @inbounds nd[i] = NodeData(nd[i].nodetype, nd[i].index, nd[i].parent + 2)
        end
    else
        pushfirst!(nd, NodeData(JuMP._Derivatives.CALL, 2, -1))
        nd[2] = NodeData(nd[2].nodetype, nd[2].index, 1)
        for i = 3:length(nd)
            @inbounds nd[i] = NodeData(nd[i].nodetype, nd[i].index, nd[i].parent + 1)
        end
    end
    push!(nd, NodeData(JuMP._Derivatives.VARIABLE, ηi, 1))
    nlobj = ip._nlp_data.evaluator.m.nlp_data.nlobj
    nlexpr = JuMP._NonlinearExprData(copy(nlobj.nd), copy(nlobj.const_values))
    nlcons = JuMP._NonlinearConstraint(nlexpr, -Inf, 0.0)

    ip._nlp_data.evaluator.m.nlp_data.nlobj = nothing
    ip._nlp_data.evaluator.has_nlobj = false
    push!(ip._nlp_data.evaluator.m.nlp_data.nlconstr, nlcons)
    constraint_bounds = ip._nlp_data.constraint_bounds
    push!(constraint_bounds, MOI.NLPBoundsPair(-Inf, 0.0))
    ip._nlp_data = MOI.NLPBlockData(constraint_bounds, ip._nlp_data.evaluator, false)

    empty!(d._relaxed_evaluator.subexpressions)
    empty!(d._nonlinear_constr)
    add_nonlinear!(m)
    m._obj_var_slack_added = true
    return
end

"""

Performs an epigraph reformulation assuming the working_problem is a minimization problem.
"""
function reform_epigraph_min!(m::GlobalOptimizer)
    ip = m._input_problem
    m._branch_variables = fill(false, m._working_problem._variable_count)
    m._obj_mult = (ip._optimization_sense == MOI.MAX_SENSE) ? -1.0 : 1.0
    reform_epigraph_min!(m, m._working_problem, m._working_problem._objective)
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
function label_fixed_variables!(m::GlobalOptimizer)
    map!(x -> check_set_is_fixed(x), m._fixed_variable, m._working_problem._variable_info)
end

"""
$(TYPEDSIGNATURES)

Detects any variables participating in nonconvex terms and populates the
`_branch_variables` storage array.
"""
function label_branch_variables!(m::GlobalOptimizer)

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
        elseif m._branch_variables[i] || wp._variable_info[i].is_integer 
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
    (wp._objective isa BufferedNonlinearFunction) && _set_variable_storage!(wp._objective, v)
    foreach(i -> _set_variable_storage!(i, v), wp._nonlinear_constr)
    foreach(i -> _set_variable_storage!(i, v), wp._relaxed_evaluator.subexpressions)
    return
end

"""
Translates input problem to working problem. Routines and checks and optional manipulation is left to the presolve stage.
"""
function initial_parse!(m::Optimizer{R,S,T}) where {R,S,T}

    # reset initial time and solution statistics
    m._global_optimizer._time_left = m._parameters.time_limit

    # add variables to working model
    ip = m._global_optimizer._input_problem = m._input_problem
    wp = m._global_optimizer._working_problem = m._working_problem
    m._global_optimizer._parameters = m._parameters
    append!(wp._variable_info, ip._variable_info)
    wp._variable_count = ip._variable_count

    # add linear constraints to the working problem
    append!(wp._saf_leq, [AffineFunctionIneq(c[1], c[2]) for c in ip._linear_leq_constraints])
    append!(wp._saf_leq, [AffineFunctionIneq(c[1], c[2]) for c in ip._linear_geq_constraints])
    wp._saf_eq  = [AffineFunctionEq(c[1], c[2]) for c in ip._linear_eq_constraints]

    append!(wp._sqf_leq, [BufferedQuadraticIneq(c[1], c[2]) for c in ip._quadratic_leq_constraints])
    append!(wp._sqf_leq, [BufferedQuadraticIneq(c[1], c[2]) for c in ip._quadratic_geq_constraints])
    wp._sqf_eq  = [BufferedQuadraticEq(c[1], c[2]) for c in ip._quadratic_eq_constraints]

    add_objective!(wp, ip._objective)    # set objective function
    add_nonlinear!(m._global_optimizer)         # add nonlinear constraints, evaluator, subexpressions

    # converts a maximum problem to a minimum problem (internally) if necessary
    # this is placed after adding nonlinear functions as this prior routine
    # copies the nlp_block from the input_problem to the working problem
    reform_epigraph_min!(m._global_optimizer)

    label_fixed_variables!(m._global_optimizer)
    label_branch_variables!(m._global_optimizer)

    # updates run and parse times
    new_time = time() - m._global_optimizer._start_time
    m._global_optimizer._parse_time = new_time
    m._global_optimizer._run_time = new_time
    return
end

"""
Classifies the problem type
"""
function parse_classify_problem!(m::GlobalOptimizer)

    ip = m._input_problem

    nl_expr_num = 0
    if (ip._objective === nothing) && (ip._nlp_data !== nothing)
        if ip._nlp_data.has_objective
            nl_expr_num += 1
            has_objective = true
        else
            has_objective = false
        end
    elseif (ip._objective === nothing) && (ip._nlp_data === nothing)
        has_objective = false
    end
    nl_expr_num += length(m._working_problem._nonlinear_constr)
    cone_constr_num = length(ip._conic_second_order)
    quad_constr_num = length(ip._quadratic_leq_constraints) +
                      length(ip._quadratic_geq_constraints) +
                      length(ip._quadratic_eq_constraints)

    if (ip._objective !== nothing)
        has_objective = true
    end
    lin_or_sv_obj = (ip._objective isa SV || ip._objective isa SAF || !has_objective)
    relaxed_supports_soc = false

    if (cone_constr_num == 0) && (quad_constr_num == 0) && (nl_expr_num == 0) && lin_or_sv_obj
        m._working_problem._problem_type = LP()
    elseif (quad_constr_num == 0) && relaxed_supports_soc && (nl_expr_num == 0) && lin_or_sv_obj
        m._working_problem._problem_type = SOCP()
    else
        m._working_problem._problem_type = MINCVX()
    end
    return
end