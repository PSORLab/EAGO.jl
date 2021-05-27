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

"""
$(TYPEDSIGNATURES)

Loads variables, linear constraints, and empty storage for first nlp and
quadratic cut.
"""
function load_relaxed_problem!(m::Optimizer)

    relaxed_optimizer = m.relaxed_optimizer

    # add variables and indices and constraints
    wp = m._working_problem
    branch_variable_count = 0

    variable_count = wp._variable_count
    for i = 1:variable_count

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
                wp._var_eq_count += 1
            end

        else
            if vinfo.has_lower_bound
                ci_sv_gt = MOI.add_constraint(relaxed_optimizer, relaxed_variable, GT(vinfo.lower_bound))
                if is_branch_variable
                    push!(m._relaxed_variable_gt, (ci_sv_gt, branch_variable_count))
                    wp._var_geq_count += 1
                end
            end

            if vinfo.has_upper_bound
                ci_sv_lt = MOI.add_constraint(relaxed_optimizer, relaxed_variable, LT(vinfo.upper_bound))
                if is_branch_variable
                    push!(m._relaxed_variable_lt, (ci_sv_lt, branch_variable_count))
                    wp._var_leq_count += 1
                end
            end
        end
    end

    # set node index to single variable constraint index maps
    m._node_to_sv_leq_ci = fill(CI{SV,LT}(-1), branch_variable_count)
    m._node_to_sv_geq_ci = fill(CI{SV,GT}(-1), branch_variable_count)
    for i = 1:wp._var_leq_count
        ci_sv_lt, branch_index = m._relaxed_variable_lt[i]
        m._node_to_sv_leq_ci[branch_index] = ci_sv_lt
    end
    for i = 1:wp._var_geq_count
        ci_sv_gt, branch_index = m._relaxed_variable_gt[i]
        m._node_to_sv_geq_ci[branch_index] = ci_sv_gt
    end

    # set number of variables to branch on
    m._branch_variable_count = branch_variable_count

    # add linear constraints
    add_linear_constraints!(m, relaxed_optimizer)

    # sets relaxed problem objective sense to Min as all problems
    # are internally converted in Min problems in EAGO
    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    return
end

function presolve_global!(t::ExtensionType, m::Optimizer)

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

    # add storage for objective cut if quadratic or nonlinear
    wp._objective_saf.terms = copy(wp._objective.saf.terms)

    # set subgradient refinement flag
    wp._relaxed_evaluator.is_post = m._parameters.subgrad_tighten
    wp._relaxed_evaluator.subgrad_tighten = m._parameters.subgrad_tighten
    wp._relaxed_evaluator.reverse_subgrad_tighten =  m._parameters.reverse_subgrad_tighten

    m._presolve_time = time() - m._parse_time
    return
end
presolve_global!(m::Optimizer) = presolve_global!(m.ext_type, m)

"""
$(SIGNATURES)

Checks to see if current node should be reprocessed.
"""
repeat_check(t::ExtensionType, m::Optimizer) = false
repeat_check(m::Optimizer) = repeat_check(m.ext_type, m)

relative_gap(L::Float64, U::Float64) = ((L > -Inf) && (U < Inf)) ?  abs(U - L)/(max(abs(L), abs(U))) : Inf
relative_tolerance(L::Float64, U::Float64, tol::Float64) = relative_gap(L, U)  > tol || ~(L > -Inf)

"""
$(SIGNATURES)

Checks for termination of algorithm due to satisfying absolute or relative
tolerance, infeasibility, or a specified limit, returns a boolean valued true
if algorithm should continue.
"""
function termination_check(t::ExtensionType, m::Optimizer)
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
termination_check(m::Optimizer) = termination_check(m.ext_type, m)

const GLOBALEND_TSTATUS = Dict{GlobalEndState, MOI.TerminationStatusCode}(
        GS_OPTIMAL => MOI.OPTIMAL,
        GS_INFEASIBLE => MOI.INFEASIBLE,
        GS_NODE_LIMIT => MOI.NODE_LIMIT,
        GS_ITERATION_LIMIT => MOI.ITERATION_LIMIT,
        GS_RELATIVE_TOL => MOI.OPTIMAL,
        GS_ABSOLUTE_TOL => MOI.OPTIMAL,
        GS_TIME_LIMIT => MOI.TIME_LIMIT
        )

function set_termination_status!(m::Optimizer)
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

function set_result_status!(m::Optimizer)
    m._result_status_code = GLOBALEND_PSTATUS[m._end_state]
    return
end

"""
$(SIGNATURES)

Checks for convergence of algorithm with respect to absolute and/or relative
tolerances.
"""
function convergence_check(t::ExtensionType, m::Optimizer)

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
convergence_check(m::Optimizer) = convergence_check(m.ext_type, m)

"""
$(SIGNATURES)

Default postprocess perfoms duality-based bound tightening on the `y`.
"""
function postprocess!(t::ExtensionType, m::Optimizer)
    if m._parameters.dbbt_depth > m._iteration_count
        variable_dbbt!(m._current_node, m._lower_lvd, m._lower_uvd,
                       m._lower_objective_value, m._global_upper_bound,
                       m._branch_variable_count)
    end
    return
end
postprocess!(m::Optimizer) = postprocess!(m.ext_type, m)

"""
$(SIGNATURES)

Provides a hook for extensions to EAGO as opposed to standard global, local,
or linear solvers.
"""
optimize_hook!(t::ExtensionType, m::Optimizer) = nothing

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
function global_solve!(m::Optimizer)

    m._iteration_count = 1
    m._node_count = 1

    parse_global!(m)
    presolve_global!(m)

    logging_on = m._parameters.log_on
    verbosity = m._parameters.verbosity

    # terminates when max nodes or iteration is reach, or when node stack is empty
    while !termination_check(m)

        # Selects node, deletes it from stack, prints based on verbosity
        node_selection!(m)
        (verbosity >= 3) && print_node!(m)

        # Performs prepocessing and times
        logging_on && (start_time = time())
        preprocess!(m)
        if logging_on
            m._last_preprocess_time = time() - start_time
        end

        if m._preprocess_feasibility

            # solves & times lower bounding problem
            logging_on && (start_time = time())
            lower_problem!(m)
            if logging_on
                m._last_lower_problem_time = time() - start_time
            end
            print_results!(m, true)

            # checks for infeasibility stores solution
            if m._lower_feasibility && !convergence_check(m)

                logging_on && (start_time = time())
                upper_problem!(m)
                if logging_on
                    m._last_upper_problem_time = time() - start_time
                end
                print_results!(m, false)
                store_candidate_solution!(m)
                if m._input_problem._optimization_sense === MOI.FEASIBILITY_SENSE
                    if !m.feasible_local_continue || m.local_solve_only
                        break
                    end
                end

                # Performs and times post processing
                logging_on && (start_time = time())
                postprocess!(m)
                if logging_on
                    m._last_postprocessing_time = time() - start_time
                end

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

    revert_adjusted_upper_bound!(m)
    m._objective_value = m._global_upper_bound
    set_termination_status!(m)
    set_primal_status!(m)

    print_solution!(m)
end

optimize!(::Val{MINCVX}, m::Optimizer) = global_solve!(m)
