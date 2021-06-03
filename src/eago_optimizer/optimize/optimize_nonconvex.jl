# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimize/optimize_nonconvex.jl
# Contains the optimize! routine and subroutines needed in the branch and
# bound routine called by EAGO.
#############################################################################

include(joinpath(@__DIR__,"nonconvex","stack_management.jl"))
include(joinpath(@__DIR__,"nonconvex","lower_problem.jl"))
include(joinpath(@__DIR__,"nonconvex","upper_problem.jl"))
include(joinpath(@__DIR__,"nonconvex","postprocess.jl"))
include(joinpath(@__DIR__,"nonconvex","log_iteration.jl"))
include(joinpath(@__DIR__,"nonconvex","display.jl"))
include(joinpath(@__DIR__,"nonconvex","configure_subsolver.jl"))

"""

Basic parsing for global solutions (no extensive manipulation)
"""
parse_global!(t::ExtensionType, m::GlobalOptimizer) = nothing
parse_global!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}  = parse_global!(_ext_typ(m), m)

"""
$(TYPEDSIGNATURES)

Loads variables, linear constraints, and empty storage for first nlp and
quadratic cut.
"""
function load_relaxed_problem!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    relaxed_optimizer = _relaxed_optimizer(m)

    # add variables and indices and constraints
    wp = m._working_problem
    branch_variable_count = 0

    for i = 1:_variable_num(FullVar(), m)

        relaxed_variable_indx = MOI.add_variable(relaxed_optimizer)
        relaxed_variable = SV(relaxed_variable_indx)
        push!(m._relaxed_variable_index, relaxed_variable_indx)

        vinfo =  wp._variable_info[i]

        is_branch_variable =  m._branch_variables[i]
        is_branch_variable && (branch_variable_count += 1)

        if vinfo.is_integer
        elseif vinfo.is_fixed
            ci_sv_et = MOI.add_constraint(relaxed_optimizer, relaxed_variable, ET(vinfo.lower_bound))
            if is_branch_variable
                push!(m._relaxed_variable_eq, (ci_sv_et, branch_variable_count))
            end
        else
            if vinfo.has_lower_bound
                ci_sv_gt = MOI.add_constraint(relaxed_optimizer, relaxed_variable, GT(vinfo.lower_bound))
                if is_branch_variable
                    push!(m._relaxed_variable_gt, (ci_sv_gt, branch_variable_count))
                end
            end
            if vinfo.has_upper_bound
                ci_sv_lt = MOI.add_constraint(relaxed_optimizer, relaxed_variable, LT(vinfo.upper_bound))
                if is_branch_variable
                    push!(m._relaxed_variable_lt, (ci_sv_lt, branch_variable_count))
                end
            end
        end
    end

    # set node index to single variable constraint index maps
    m._node_to_sv_leq_ci = fill(CI{SV,LT}(-1), branch_variable_count)
    m._node_to_sv_geq_ci = fill(CI{SV,GT}(-1), branch_variable_count)
    for v in m._relaxed_variable_lt
        m._node_to_sv_leq_ci[v[2]] = v[1]
    end
    for v in m._relaxed_variable_gt
        m._node_to_sv_geq_ci[v[2]] = v[1]
    end

    # set number of variables to branch on
    m._branch_variable_count = branch_variable_count

    # add linear constraints
    add_linear_constraints!(m, relaxed_optimizer)

    # sets relaxed problem objective sense to Min as all problems
    # are internally converted in Min problems in EAGO
    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), wp._objective_saf)

    return
end

function presolve_global!(t::ExtensionType, m::GlobalOptimizer)

    set_default_config!(m)
    load_relaxed_problem!(m)
    initialize_stack!(m)

    wp = m._working_problem
    branch_variable_count = m._branch_variable_count

    m._current_xref             = fill(0.0, branch_variable_count)
    m._candidate_xref           = fill(0.0, branch_variable_count)
    m._current_objective_xref   = fill(0.0, branch_variable_count)
    m._prior_objective_xref     = fill(0.0, branch_variable_count)
    m._lower_lvd                = fill(0.0, branch_variable_count)
    m._lower_uvd                = fill(0.0, branch_variable_count)

    # populate in full space until local MOI nlp solves support constraint deletion
    # uses input model for local nlp solves... may adjust this if a convincing reason
    # to use a reformulated upper problem presents itself
    m._lower_solution      = zeros(Float64, wp._variable_count)
    m._continuous_solution = zeros(Float64, wp._variable_count)
    m._upper_solution      = zeros(Float64, wp._variable_count)
    m._upper_variables     = fill(VI(-1), wp._variable_count)

    # add storage for fbbt
    m._lower_fbbt_buffer   = zeros(Float64, wp._variable_count)
    m._upper_fbbt_buffer   = zeros(Float64, wp._variable_count)

    # add storage for obbt ( perform obbt on all relaxed variables, potentially)
    m._obbt_working_lower_index = fill(false, branch_variable_count)
    m._obbt_working_upper_index = fill(false, branch_variable_count)
    m._old_low_index            = fill(false, branch_variable_count)
    m._old_upp_index            = fill(false, branch_variable_count)
    m._new_low_index            = fill(false, branch_variable_count)
    m._new_upp_index            = fill(false, branch_variable_count)
    m._lower_indx_diff          = fill(false, branch_variable_count)
    m._upper_indx_diff          = fill(false, branch_variable_count)
    m._obbt_variable_count      = branch_variable_count

    # set subgradient refinement flag
    wp._relaxed_evaluator.is_post = m._parameters.subgrad_tighten
    wp._relaxed_evaluator.subgrad_tighten = m._parameters.subgrad_tighten
    wp._relaxed_evaluator.reverse_subgrad_tighten =  m._parameters.reverse_subgrad_tighten

    m._presolve_time = time() - m._parse_time
    return
end
presolve_global!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = presolve_global!(_ext_typ(m), m)

"""
$(SIGNATURES)

Checks for termination of algorithm due to satisfying absolute or relative
tolerance, infeasibility, or a specified limit, returns a boolean valued true
if algorithm should continue.
"""
function termination_check(t::ExtensionType, m::GlobalOptimizer)
    nlen = length(m._stack)
    L = m._global_lower_bound
    U = m._global_upper_bound
    if nlen == 0 && m._first_solution_node > 0
        m._end_state = GS_OPTIMAL
    elseif nlen == 0 && !(m._first_solution_node > 0)
        m._end_state = GS_INFEASIBLE
    elseif nlen >= m._parameters.node_limit
        m._end_state = GS_NODE_LIMIT
    elseif m._iteration_count >= m._parameters.iteration_limit
        m._end_state = GS_ITERATION_LIMIT
    elseif !relative_tolerance(L, U, m._parameters.relative_tolerance)
        m._end_state = GS_RELATIVE_TOL
    elseif (U - L) < m._parameters.absolute_tolerance
        m._end_state = GS_ABSOLUTE_TOL
    elseif m._run_time > m._parameters.time_limit
        m._end_state = GS_TIME_LIMIT
    else
        return false
    end
    return true
end
termination_check(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = termination_check(_ext_typ(m), m)

const GLOBALEND_TSTATUS = Dict{GlobalEndState, MOI.TerminationStatusCode}(
        GS_OPTIMAL => MOI.OPTIMAL,
        GS_INFEASIBLE => MOI.INFEASIBLE,
        GS_NODE_LIMIT => MOI.NODE_LIMIT,
        GS_ITERATION_LIMIT => MOI.ITERATION_LIMIT,
        GS_RELATIVE_TOL => MOI.OPTIMAL,
        GS_ABSOLUTE_TOL => MOI.OPTIMAL,
        GS_TIME_LIMIT => MOI.TIME_LIMIT
        )

function set_termination_status!(m::GlobalOptimizer)
    m._termination_status_code = GLOBALEND_TSTATUS[m._end_state]
    return
end

const GLOBALEND_PSTATUS = Dict{GlobalEndState, MOI.ResultStatusCode}(
        GS_OPTIMAL => MOI.FEASIBLE_POINT,
        GS_INFEASIBLE => MOI.INFEASIBILITY_CERTIFICATE,
        GS_NODE_LIMIT => MOI.UNKNOWN_RESULT_STATUS,
        GS_ITERATION_LIMIT => MOI.UNKNOWN_RESULT_STATUS,
        GS_RELATIVE_TOL => MOI.FEASIBLE_POINT,
        GS_ABSOLUTE_TOL => MOI.FEASIBLE_POINT,
        GS_TIME_LIMIT => MOI.UNKNOWN_RESULT_STATUS
        )

function set_result_status!(m::GlobalOptimizer)
    m._result_status_code = GLOBALEND_PSTATUS[m._end_state]
    return
end

"""
$(SIGNATURES)

Checks for convergence of algorithm with respect to absolute and/or relative
tolerances.
"""
function convergence_check(t::ExtensionType, m::GlobalOptimizer)

  L = m._lower_objective_value
  U = m._global_upper_bound
  t = (U - L) <= m._parameters.absolute_tolerance
  if (U < Inf) && (L > Inf)
      t |= (abs(U - L)/(max(abs(L), abs(U))) <= m._parameters.relative_tolerance)
  end

  if t && m._min_converged_value < Inf
      m._min_converged_value = min(m._min_converged_value, L)
  else
      m._min_converged_value = L
  end

  return t
end
convergence_check(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = convergence_check(_ext_typ(m), m)

"""
$(SIGNATURES)

Provides a hook for extensions to EAGO as opposed to standard global, local,
or linear solvers.
"""
optimize_hook!(t::ExtensionType, m::GlobalOptimizer) = nothing

function store_candidate_solution!(m::GlobalOptimizer)
    if m._upper_feasibility && (m._upper_objective_value < m._global_upper_bound)
        m._feasible_solution_found = true
        m._first_solution_node = m._maximum_node_id
        m._global_upper_bound = m._upper_objective_value
        @__dot__ m._continuous_solution = m._upper_solution
    end
    return
end

function set_global_lower_bound!(m::GlobalOptimizer)
    if !isempty(m._stack)
        n = minimum(m._stack)
        lower_bound = n.lower_bound
        if m._global_lower_bound < lower_bound
            m._global_lower_bound = lower_bound
        end
    end
    return
end

"""
$(TYPEDSIGNATURES)

Solves the branch and bound problem with the input EAGO optimizer object.
"""
function global_solve!(m::GlobalOptimizer)

    m._iteration_count = 1
    m._node_count = 1

    parse_global!(m)
    presolve_global!(m)

    print_preamble!(m)

    # terminates when max nodes or iteration is reach, or when node stack is empty
    while !termination_check(m)

        # Selects node, deletes it from stack, prints based on verbosity
        node_selection!(m)
        print_node!(m)

        # Performs prepocessing and times
        m._last_preprocess_time += @elapsed preprocess!(m)

        if m._preprocess_feasibility

            # solves & times lower bounding problem
            m._last_lower_problem_time += @elapsed lower_problem!(m)
            print_results!(m, true)

            # checks for infeasibility stores solution
            if m._lower_feasibility && !convergence_check(m)

                # Solves upper problem
                m._last_upper_problem_time += @elapsed upper_problem!(m)
                print_results!(m, false)
                store_candidate_solution!(m)

                # Performs post processing
                m._last_postprocessing_time += @elapsed postprocess!(m)

                # Checks to see if the node
                if m._postprocess_feasibility
                    repeat_check(m) ? single_storage!(m) : branch_node!(m)
                end
            end
            fathom!(m)
        else
            m._lower_objective_value = -Inf
            m._lower_feasibility = false
            m._upper_feasibility = false
        end
        set_global_lower_bound!(m)
        m._run_time = time() - m._start_time
        m._time_left = m._parameters.time_limit - m._run_time
        log_iteration!(m)
        print_iteration!(m)
        m._iteration_count += 1
    end

    set_termination_status!(m)
    set_result_status!(m)

    print_solution!(m)
end

function unpack_global_solution!(m::Optimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    g = m._global_optimizer
    
    m._termination_status_code = g._termination_status_code
    m._result_status_code      = g._result_status_code

    m._run_time = g._run_time
    m._node_count = g._maximum_node_id


    if g._input_problem._optimization_sense == MOI.MIN_SENSE
        m._objective_bound = g._global_lower_bound
        m._objective_value = g._global_upper_bound
    end
    m._objective_bound = -g._global_upper_bound
    m._objective_value = -g._global_lower_bound
    return
end

function optimize!(::MINCVX, m::Optimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    global_solve!(m._global_optimizer)
    unpack_global_solution!(m)
    return
end
