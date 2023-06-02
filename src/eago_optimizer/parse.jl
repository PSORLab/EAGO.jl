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
    add_objective!

Add objective function (if any) to the parsed problem.
"""
add_objective!(m::ParsedProblem, f) = nothing
add_objective!(m::ParsedProblem, f::VI) = (m._objective = f; nothing)
add_objective!(m::ParsedProblem, f::SAF) = (m._objective = AffineFunctionIneq(f, LT_ZERO); nothing)
add_objective!(m::ParsedProblem, f::SQF) = (m._objective = BufferedQuadraticIneq(f, LT_ZERO); nothing)

"""
    add_nonlinear_evaluator!

Add an Evaluator structure if nonlinear terms are attached.
"""
add_nonlinear_evaluator!(m::GlobalOptimizer, evaluator::Nothing) = nothing
function add_nonlinear_evaluator!(m::GlobalOptimizer, d::MOI.AbstractNLPEvaluator)
    m._working_problem._relaxed_evaluator = Evaluator()
    relax_eval = m._working_problem._relaxed_evaluator
    relax_eval.user_operators = OperatorRegistry(d.model.operators)
    relax_eval.subgrad_tol    = m._parameters.subgrad_tol
    m._nonlinear_evaluator_created = true
    return
end
add_nonlinear_evaluator!(m::GlobalOptimizer, nldata) = add_nonlinear_evaluator!(m, nldata.evaluator)
add_nonlinear_evaluator!(m::GlobalOptimizer) = add_nonlinear_evaluator!(m, m._input_problem._nlp_data)


function link_subexpression_dicts!(m::GlobalOptimizer)
    evaluator = m._working_problem._relaxed_evaluator
    n_subexpr = length(evaluator.subexpressions)

    dn = Dict{Int,Float64}()
    din = Dict{Int,Bool}()
    for i = 1:n_subexpr
        dn[i] = 0.0
        din[i] = false
    end
    for ex in evaluator.subexpressions
        mctyp = mc_type(ex)
        ds = Dict{Int,mctyp}()
        di = Dict{Int,mctyp}()
        for i = 1:n_subexpr
            ds[i] = zero(mctyp)
            di[i] = zero(mctyp)
        end
        copy_subexpr!(ex.relax_cache, ds, dn, din, di)
    end
    if m._working_problem._objective isa BufferedNonlinearFunction
        ex = m._working_problem._objective.ex
        mctyp = mc_type(ex)
        ds = Dict{Int,mctyp}()
        di = Dict{Int,mctyp}()
        for i = 1:n_subexpr
            ds[i] = zero(mctyp)
            di[i] = zero(mctyp)
        end
        copy_subexpr!(ex.relax_cache, ds, dn, din, di)
    end
    for f in m._working_problem._nonlinear_constr
        mctyp = mc_type(f.ex)
        ds = Dict{Int,mctyp}()
        di = Dict{Int,mctyp}()
        for i = 1:n_subexpr
            ds[i] = zero(mctyp)
            di[i] = zero(mctyp)
        end
        copy_subexpr!(f.ex.relax_cache, ds, dn, din, di)
    end
    return
end

add_nonlinear!(m::GlobalOptimizer, evaluator::Nothing) = nothing
function add_nonlinear!(m::GlobalOptimizer, evaluator::MOI.AbstractNLPEvaluator)

    add_nonlinear_evaluator!(m, m._input_problem._nlp_data.evaluator)

    nlp_data = m._input_problem._nlp_data

    MOI.initialize(evaluator, Symbol[:Grad, :Jac, :ExprGraph])

    user_operator_registry = OperatorRegistry(evaluator.model.operators)

    # set nlp data structure
    m._working_problem._nlp_data = nlp_data
    mul_relax = m._parameters.mul_relax_style
    if mul_relax == 1
        rtype = Relax()
        renum = STD_RELAX
        ruse_apriori = true
    elseif mul_relax == 2
        rtype = RelaxAA()
        renum = MC_AFF_RELAX
        ruse_apriori = true
    elseif mul_relax == 3
        rtype = RelaxMulEnum()
        renum = MC_ENUM_RELAX
        ruse_apriori = true
    else
        rtype = Relax()
        renum = STD_RELAX
        ruse_apriori = false
    end

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
    if length(evaluator.model.expressions) > 0      # should check for nonlinear objective, constraint
        for i = 1:length(evaluator.backend.subexpressions)
            subexpr = evaluator.backend.subexpressions[i]
            nlexpr = NonlinearExpression!(m._auxiliary_variable_info, rtype, subexpr, MOI.NLPBoundsPair(-Inf, Inf),
                                          dict_sparsity, i, evaluator.backend.subexpression_linearity, 
                                          user_operator_registry, evaluator.model.parameters,
                                          m._parameters.relax_tag, ruse_apriori; is_sub = true)
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
    if evaluator.model.objective !== nothing
        m._working_problem._objective = BufferedNonlinearFunction(m._auxiliary_variable_info, rtype, evaluator.backend.objective, MOI.NLPBoundsPair(-Inf, Inf),
                                                                 dict_sparsity, evaluator.backend.subexpression_linearity,
                                                                 user_operator_registry,
                                                                 evaluator.model.parameters,
                                                                 m._parameters.relax_tag, ruse_apriori)
    end

    # add nonlinear constraints
    constraint_bounds = m._working_problem._nlp_data.constraint_bounds
    for i = 1:length(evaluator.model.constraints)
        constraint = evaluator.backend.constraints[i]
        bnds = constraint_bounds[i]
        push!(m._working_problem._nonlinear_constr, BufferedNonlinearFunction(m._auxiliary_variable_info, rtype, constraint, bnds, dict_sparsity,
                                                                              evaluator.backend.subexpression_linearity,
                                                                              user_operator_registry,
                                                                              evaluator.model.parameters,
                                                                              m._parameters.relax_tag, ruse_apriori))
    end

    relax_evaluator.subexpressions_eval = fill(false, length(relax_evaluator.subexpressions))
    return
end
add_nonlinear!(m::GlobalOptimizer, nldata) = add_nonlinear!(m, nldata.evaluator)

"""
$(TYPEDSIGNATURES)

Add an Evaluator and nonlinear functions and populate each appropriately.
"""
add_nonlinear!(m::GlobalOptimizer) = add_nonlinear!(m, m._input_problem._nlp_data)

function reform_epigraph_min!(m::GlobalOptimizer, d::ParsedProblem, f::VI)
    flag = m._input_problem._optimization_sense == MOI.MAX_SENSE
    d._objective = AffineFunctionIneq(f, is_max = flag)
    d._objective_saf = SAF([SAT(flag ? -1.0 : 1.0, VI(f.value))], 0.0)
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

function add_η!(m::ParsedProblem)
    m._variable_count += 1
    push!(m._variable_info, VariableInfo{Float64}())
    return m._variable_count
end

set_variable_values!(wp, v) = set_variable_values!(wp._relaxed_evaluator, v)
reform_epigraph_min!(m::GlobalOptimizer, d::ParsedProblem, f::Nothing) = nothing
function reform_epigraph_min!(m::GlobalOptimizer, d::ParsedProblem, f::BufferedQuadraticIneq)
    ip = m._input_problem
    wp = m._working_problem

    vi = wp._variable_info
    q = _variable_num(FullVar(), m)
    v = VariableValues{Float64}(x = mid.(vi), x0 = mid.(vi),
                                lower_variable_bounds = lower_bound.(vi),
                                upper_variable_bounds = upper_bound.(vi),
                                node_to_variable_map = [i for i in 1:q],
                                variable_to_node_map = [i for i in 1:q])
    set_variable_values!(wp, v)
    
    n = NodeBB(lower_bound.(vi), upper_bound.(vi), is_integer.(vi))
    m._current_node = n
    set_node!(wp._relaxed_evaluator, n)

    ηi = add_η!(d)
    if !isnothing(ip._nlp_data)
        # add variable
        η = MOI.VariableIndex(ip._variable_count + 1)
        push!(ip._nlp_data.evaluator.backend.ordered_variables, η)
        # add objective
        MOINL.set_objective(ip._nlp_data.evaluator.model, :($η))
    end

    sqf_obj = copy(m._input_problem._objective)
    if !_is_input_min(m)
        MOIU.operate!(-, Float64, sqf_obj)
    end
    d._objective_saf = SAF([SAT(1.0, VI(ηi))], 0.0)    
    push!(sqf_obj.affine_terms, SAT(-1.0, VI(ηi)))
    _constraints(ip, SQF, LT)[CI{SQF,LT}(ip._constraint_count += 1)] = (sqf_obj, LT(0.0))
    push!(wp._sqf_leq, BufferedQuadraticIneq(sqf_obj, LT(0.0)))

    f.buffer[ηi] = 0.0
    f.len += 1
    if !_is_input_min(m)
        MOIU.operate!(-, Float64, f.func)
    end
    MOIU.operate!(-, Float64, f.func, VI(ηi))
    push!(f.saf.terms, SAT(-1.0, VI(ηi)))
    m._obj_var_slack_added = true
    return
end
function reform_epigraph_min!(m::GlobalOptimizer, d::ParsedProblem, f::BufferedNonlinearFunction)
    ip = m._input_problem
    wp = m._working_problem

    vi = wp._variable_info
    q = _variable_num(FullVar(), m)
    v = VariableValues{Float64}(x = mid.(vi), x0 = mid.(vi),
                                lower_variable_bounds = lower_bound.(vi),
                                upper_variable_bounds = upper_bound.(vi),
                                node_to_variable_map = [i for i in 1:q],
                                variable_to_node_map = [i for i in 1:q]) 
   
    ηi = add_η!(d)
    η = MOI.VariableIndex(ip._variable_count + 1)
    push!(ip._nlp_data.evaluator.backend.ordered_variables, η)
    expr = MOINL.convert_to_expr(ip._nlp_data.evaluator.model,ip._nlp_data.evaluator.model.objective)
    MOINL.set_objective(ip._nlp_data.evaluator.model, :($η))
    wp._objective_saf = SAF([SAT(1.0, VI(ηi))], 0.0)

    # check if input problem is min or max
    if !_is_input_min(m)
        cons = MOINL.add_expression(ip._nlp_data.evaluator.model, :(-$expr - $η))
        MOINL.add_constraint(ip._nlp_data.evaluator.model, :($cons), MOI.LessThan(0.0))
    else
        cons = MOINL.add_expression(ip._nlp_data.evaluator.model, :($expr - $η))
        MOINL.add_constraint(ip._nlp_data.evaluator.model, :($cons), MOI.LessThan(0.0))
    end

    ip._nlp_data.evaluator.model.objective = nothing
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
$(TYPEDSIGNATURES)

Perform an epigraph reformulation assuming the working_problem is a minimization problem.
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

Detect any variables set to a fixed value by equality or inequality constraints
and populate the `_fixed_variable` storage array.
"""
label_fixed_variables!(m::GlobalOptimizer) = map!(x -> check_set_is_fixed(x), m._fixed_variable, m._working_problem._variable_info)

"""
$(TYPEDSIGNATURES)

Detect any variables participating in nonconvex terms and populate the
`_branch_variables` storage array.
"""
function label_branch_variables!(m::GlobalOptimizer)

    wp = m._working_problem
    m._user_branch_variables = !isempty(m._parameters.branch_variable)
    if m._user_branch_variables
        m._branch_variables = m._parameters.branch_variable
        if length(wp._variable_info) > length(m._branch_variables) #Should only need 1
            push!(m._parameters.branch_variable, true)
            push!(m._branch_variables, true)
        end
    else
        m._branch_variables = fill(false, m._working_problem._variable_count)
        for f in wp._sqf_leq, t in f.func.quadratic_terms
            m._branch_variables[t.variable_1.value] = true
            m._branch_variables[t.variable_2.value] = true
        end
        for f in wp._sqf_eq, t in f.func.quadratic_terms
            m._branch_variables[t.variable_1.value] = true
            m._branch_variables[t.variable_2.value] = true
        end
        for f in wp._nonlinear_constr, i in sparsity(f)
             m._branch_variables[i] = true
        end
        if wp._objective isa BufferedQuadraticIneq
            for t in wp._objective.func.quadratic_terms
                m._branch_variables[t.variable_1.value] = true
                m._branch_variables[t.variable_2.value] = true
            end
        elseif wp._objective isa BufferedNonlinearFunction
            for i in sparsity(wp._objective)
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
        elseif i == wp._variable_count
            m._branch_variables[i] = false
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
    l = lower_bound.(m._working_problem._variable_info)
    u = upper_bound.(m._working_problem._variable_info)
    v = VariableValues{Float64}(x = zeros(vnum),
                                x0 = zeros(vnum),
                                lower_variable_bounds = l,
                                upper_variable_bounds = u,
                                node_to_variable_map = m._branch_to_sol_map,
                                variable_to_node_map = m._sol_to_branch_map)

    wp._relaxed_evaluator.variable_values = v
    (wp._objective isa BufferedNonlinearFunction) && set_variable_storage!(wp._objective, v)
    foreach(i -> set_variable_storage!(i, v), wp._nonlinear_constr)
    foreach(i -> set_variable_storage!(i, v), wp._relaxed_evaluator.subexpressions)
    return
end

"""
    variable_load_parse!

Parse constraint information of a given constraint type to fill in related variable information.
"""
function variable_load_parse!(m::Optimizer, ::Type{VI}, ::Type{T}) where T
    wp = m._global_optimizer._working_problem = m._working_problem
    for (i, v) in enumerate(values(_constraints(m, VI, T)))
        wp._variable_info[v[1].value] = VariableInfo(wp._variable_info[v[1].value], v[2])
    end
    return
end

"""
$(TYPEDSIGNATURES)

Translate the input problem to the working problem. Any checks or optional
manipulation are left to the presolve stage.
"""
function initial_parse!(m::Optimizer{R,S,T}) where {R,S,T}

    # Copy time information to the global optimizer
    m._global_optimizer._time_left = m._parameters.time_limit

    # Copy input and working problem info to the global optimizer
    ip = m._global_optimizer._input_problem = m._input_problem
    wp = m._global_optimizer._working_problem = m._working_problem
    m._global_optimizer._parameters = m._parameters

    # Parse different constraint types to add information to variables in
    # the working problem
    wp._variable_info = VariableInfo{Float64}[VariableInfo{Float64}() for i=1:ip._variable_count]
    variable_load_parse!(m, VI, LT)
    variable_load_parse!(m, VI, GT)
    variable_load_parse!(m, VI, ET)
    variable_load_parse!(m, VI, ZO)
    variable_load_parse!(m, VI, MOI.Integer)
    wp._variable_count = ip._variable_count

    # Copy all constraint dictionaries from the input problem to the working problem
    for (f, s) in values(ip._linear_leq_constraints)
        push!(wp._saf_leq, AffineFunctionIneq(f, s))
    end
    for (f, s) in values(ip._linear_geq_constraints)
        push!(wp._saf_leq, AffineFunctionIneq(f, s))
    end
    for (f, s) in values(ip._linear_eq_constraints)
        push!(wp._saf_eq, AffineFunctionEq(f, s))
    end

    for (f, s) in values(ip._quadratic_leq_constraints)
        push!(wp._sqf_leq, BufferedQuadraticIneq(f, s))
    end
    for (f, s) in values(ip._quadratic_geq_constraints)
        push!(wp._sqf_leq, BufferedQuadraticIneq(f, s))
    end
    for (f, s) in values(ip._quadratic_eq_constraints)
        push!(wp._sqf_eq, BufferedQuadraticEq(f, s))
    end

    # Copy the input problem objective to the working problem
    add_objective!(wp, ip._objective)

    # Add an evaluator, nonlinear constraints, and subexpressions
    add_nonlinear!(m._global_optimizer)

    # Convert maximum problem to a minimum problem (internally) if necessary.
    # This is done after adding nonlinear functions as `add_nonlinear!` copies
    # the nlp_block from the input problem to the working problem
    reform_epigraph_min!(m._global_optimizer)

    # Identify any fixed variables and variables to branch on
    label_fixed_variables!(m._global_optimizer)
    label_branch_variables!(m._global_optimizer)

    link_subexpression_dicts!(m._global_optimizer)

    # Update run and parse times
    new_time = time() - m._global_optimizer._start_time
    m._global_optimizer._run_time = new_time
    m._global_optimizer._parse_time = new_time
    return
end

"""
    parse_classify_problem!(m::GlobalOptimizer)

Interprets the type/number of constraints and the type of objective function
to infer a problem type. Current possible types include:

- `LP`: Linear program; sent to `optimize_lp.jl`
- 'MILP`: Mixed integer linear program; sent to `optimize_lp.jl`
- `SOCP`: Second-order cone program; sent to `optimize_conic.jl`
- `MINCVX`: Mixed-integer nonconvex; sent to `optimize_nonconvex.jl`

If the `force_global_solve` parameter is set to `true`, `parse_classify_problem!`
will set the problem type to `MINCVX` to pass the problem to `optimize_nonconvex.jl`.
"""
function parse_classify_problem!(m::GlobalOptimizer)

    if !m._parameters.force_global_solve
        ip = m._input_problem

        nl_expr_num = 0
        if isnothing(ip._objective) && !isnothing(ip._nlp_data)
            if ip._nlp_data.has_objective
                nl_expr_num += 1
                has_objective = true
            else
                has_objective = false
            end
        elseif isnothing(ip._objective) && isnothing(ip._nlp_data)
            has_objective = false
        end
        nl_expr_num += length(m._working_problem._nonlinear_constr)
        cone_constr_num = length(ip._conic_second_order)
        quad_constr_num = length(ip._quadratic_leq_constraints) +
                        length(ip._quadratic_geq_constraints) +
                        length(ip._quadratic_eq_constraints)

        if !isnothing(ip._objective)
            has_objective = true
        end
        has_int_var = !iszero(length(ip._vi_zo_constraints) + length(ip._vi_int_constraints))

        lin_or_sv_obj = (ip._objective isa VI || ip._objective isa SAF || !has_objective)
        relaxed_supports_soc = false

        if (cone_constr_num == 0) && (quad_constr_num == 0) && (nl_expr_num == 0) && lin_or_sv_obj && !has_int_var
            m._working_problem._problem_type = LP()
        elseif (cone_constr_num == 0) && (quad_constr_num == 0) && (nl_expr_num == 0) && lin_or_sv_obj
            m._working_problem._problem_type = MILP()
        elseif (quad_constr_num == 0) && relaxed_supports_soc && (nl_expr_num == 0) && lin_or_sv_obj && !has_int_var
            m._working_problem._problem_type = SOCP()
        else
            m._working_problem._problem_type = MINCVX()
        end
    else
        m._working_problem._problem_type = MINCVX()
    end

    return
end