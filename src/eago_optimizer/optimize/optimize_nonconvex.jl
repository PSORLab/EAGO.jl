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
Basic parsing for global solutions (no extensive manipulation). By default,
does nothing.
"""
parse_global!(t::ExtensionType, m::GlobalOptimizer) = nothing
parse_global!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}  = parse_global!(_ext(m), m)

"""
$(TYPEDSIGNATURES)

Load variables, linear constraints, and empty storage space for the first NLP
and quadratic cut into the relaxed optimizer.
"""
function load_relaxed_problem!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    d = _relaxed_optimizer(m)

    # add variables and indices and constraints
    wp = m._working_problem
    branch_variable_count = 0

    full_var_num = _variable_num(FullVar(), m)
    relaxed_index_new = length(m._relaxed_variable_index) != full_var_num
    for i = 1:full_var_num
        v = MOI.add_variable(d)
        if relaxed_index_new
            push!(m._relaxed_variable_index, v)
        else
            m._relaxed_variable_index[i] = v
        end

        is_branch_variable =  m._branch_variables[i]
        is_branch_variable && (branch_variable_count += 1)

        vi = wp._variable_info[i]
        if !is_branch_variable
            if is_fixed(vi)
                MOI.add_constraint(d, v, ET(vi))
            elseif is_interval(vi)
                MOI.add_constraint(d, v, IT(vi))
            elseif is_greater_than(vi)
                MOI.add_constraint(d, v, GT(vi))
            elseif is_less_than(vi)
                MOI.add_constraint(d, v, LT(vi))
            end
        end
    end

    # TODO: Remove when upstream Cbc issue https://github.com/jump-dev/Cbc.jl/issues/168 is fixed
    # Add extra binary variable `issue_var` fixed to zero to prevent Cbc from displaying even though 
    # silent is set to off. Sets `issue_var` to zero. 
    issue_var = MOI.add_variable(d)
    MOI.add_constraint(d, issue_var, ZO())
    MOI.add_constraint(d, issue_var, ET(0.0))


    # set number of variables to branch on
    m._branch_variable_count = branch_variable_count

    # add linear constraints
    for (f, s) in collect(values(m._input_problem._linear_leq_constraints))
        MOI.add_constraint(d, f, s)
    end
    for (f, s) in collect(values(m._input_problem._linear_geq_constraints))
        MOI.add_constraint(d, f, s)
    end
    for (f, s) in collect(values(m._input_problem._linear_eq_constraints))
        MOI.add_constraint(d, f, s)
    end

    # sets relaxed problem objective sense to Min as all problems
    # are internally converted in Min problems in EAGO
    MOI.set(d, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(d, MOI.ObjectiveFunction{SAF}(), wp._objective_saf)

    return
end

"""
$(TYPEDSIGNATURES)

Perform any necessary work prior to running branch-and-bound. 

- Set subsolver configs using values in `EAGOParameters`;
- Load variables, linear constraints, and empty storage space into the relaxed optimizer;
- Prepare the stack for the start of branch-and-bound;
- Fill fields of the `GlobalOptimizer` with zeros of the proper dimensions;
- Pass necessary flags from the `GlobalOptimizer` to the working problem.
"""
function presolve_global!(t::ExtensionType, m::GlobalOptimizer)

    set_default_config!(m)
    load_relaxed_problem!(m)
    initialize_stack!(m)

    wp = m._working_problem
    branch_variable_count = m._branch_variable_count

    m._current_xref             = fill(0.0, branch_variable_count) #Note: Unused?
    m._candidate_xref           = fill(0.0, branch_variable_count) #Note: Unused?
    m._current_objective_xref   = fill(0.0, branch_variable_count) #Note: Unused?
    m._prior_objective_xref     = fill(0.0, branch_variable_count) #Note: Unused?
    m._lower_lvd                = fill(0.0, branch_variable_count)
    m._lower_uvd                = fill(0.0, branch_variable_count)

    # Populate in full space until local MOI NLP solves support constraint deletion.
    # Uses input model for local NLP solves... may adjust this if there's ever a 
    # convincing reason to use a reformulated upper problem
    m._lower_solution      = zeros(Float64, wp._variable_count)
    m._continuous_solution = zeros(Float64, wp._variable_count)
    m._upper_solution      = zeros(Float64, wp._variable_count)
    m._upper_variables     = fill(VI(-1), wp._variable_count)

    # Add storage for fbbt
    m._lower_fbbt_buffer   = zeros(Float64, wp._variable_count)
    m._upper_fbbt_buffer   = zeros(Float64, wp._variable_count)

    # Add storage for obbt (perform obbt on all relaxed variables, potentially)
    m._obbt_working_lower_index = fill(false, branch_variable_count)
    m._obbt_working_upper_index = fill(false, branch_variable_count)
    m._old_low_index            = fill(false, branch_variable_count)
    m._old_upp_index            = fill(false, branch_variable_count)
    m._new_low_index            = fill(false, branch_variable_count)
    m._new_upp_index            = fill(false, branch_variable_count)
    m._lower_indx_diff          = fill(false, branch_variable_count)
    m._upper_indx_diff          = fill(false, branch_variable_count)
    m._obbt_variable_count      = branch_variable_count

    # Set subgradient refinement flag
    wp._relaxed_evaluator.is_post = m._parameters.subgrad_tighten
    wp._relaxed_evaluator.subgrad_tighten = m._parameters.subgrad_tighten
    wp._relaxed_evaluator.reverse_subgrad_tighten =  m._parameters.reverse_subgrad_tighten

    m._presolve_time = time() - m._parse_time #TODO check that this works as expected. Isn't time() >>> _parse_time?
    return
end
presolve_global!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = presolve_global!(_ext(m), m)

"""
    termination_check(m::GlobalOptimizer)
    termination_check(t::ExtensionType, m::GlobalOptimizer) -> Bool

Check for termination of the branch-and-bound algorithm.

If only the `GlobalOptimizer` is given as an argument, `termination_check` dispatches
to the other form using the `ExtensionType` given in the `SubSolvers`. If there is no
user-defined extension, then by default, this will check for satisfaction of absolute
or relative tolerances, solution infeasibility, and other specified limits. Returns 
`true` if any conditions are met and branch-and-bound should end, and `false` otherwise.
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
    elseif abs(U - L) < m._parameters.absolute_tolerance
        m._end_state = GS_ABSOLUTE_TOL
    elseif m._run_time > m._parameters.time_limit
        m._end_state = GS_TIME_LIMIT
    else
        return false
    end
    return true
end
termination_check(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = termination_check(_ext(m), m)

const GLOBALEND_TSTATUS = Dict{GlobalEndState, MOI.TerminationStatusCode}(
        GS_OPTIMAL => MOI.OPTIMAL,
        GS_INFEASIBLE => MOI.INFEASIBLE,
        GS_NODE_LIMIT => MOI.NODE_LIMIT,
        GS_ITERATION_LIMIT => MOI.ITERATION_LIMIT,
        GS_RELATIVE_TOL => MOI.OPTIMAL,
        GS_ABSOLUTE_TOL => MOI.OPTIMAL,
        GS_TIME_LIMIT => MOI.TIME_LIMIT
        )

"""
$(TYPEDSIGNATURES)

Convert EAGO's ending status code into an `MOI.TerminationStatusCode`.
"""
function set_termination_status!(m::GlobalOptimizer)
    m._termination_status_code = GLOBALEND_TSTATUS[m._end_state]
    return
end

const GLOBALEND_PSTATUS = Dict{GlobalEndState, MOI.ResultStatusCode}(
        GS_OPTIMAL => MOI.FEASIBLE_POINT,
        GS_INFEASIBLE => MOI.NO_SOLUTION,                 # Proof of infeasibility implies not solution found
        GS_NODE_LIMIT => MOI.UNKNOWN_RESULT_STATUS,
        GS_ITERATION_LIMIT => MOI.UNKNOWN_RESULT_STATUS,
        GS_RELATIVE_TOL => MOI.FEASIBLE_POINT,
        GS_ABSOLUTE_TOL => MOI.FEASIBLE_POINT,
        GS_TIME_LIMIT => MOI.UNKNOWN_RESULT_STATUS
        )

"""
$(TYPEDSIGNATURES)

Convert EAGO's ending status code into an `MOI.ResultStatusCode`.
"""
function set_result_status!(m::GlobalOptimizer)
    m._result_status_code = GLOBALEND_PSTATUS[m._end_state]
    return
end

"""
$(TYPEDSIGNATURES)

Check for problem convergence.
    
By default, check if the lower and upper bounds have converged to within absolute
and/or relative tolerances.
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
convergence_check(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = convergence_check(_ext(m), m)

"""
    optimize_hook!(t::ExtensionType, m::Optimizer)

Provide a hook for extensions to EAGO. 

The user-defined extension of `optimize_hook!` is used in EAGO's overloading of
`MOI.optimize!` (see `EAGO.jl/src/eago_optimizer/optimize/optimize.jl`). Without
the `optimize_hook!` specified, EAGO will run `initial_parse!`, `parse_classify_problem!`,
and then `optimize!` using the parsed problem type. The user-specified `optimize_hook!` 
should thus take the new extension and `Optimizer` as inputs and will execute when the 
user writes `optimize!(model)`.

# Example
Here, `optimize_hook!` is used to bypass EAGO's problem parsing and
treat every problem using its branch-and-bound routine. This is done
in this example by telling EAGO to treat the problem as a mixed integer
non-convex problem, which normally dispatches to branch-and-bound.
```julia-repl
struct MyNewExtension <: EAGO.ExtensionType end
import EAGO: optimize_hook!
function EAGO.optimize_hook!(t::MyNewExtension, m::Optimizer)
    initial_parse!(m)
    optimize!(EAGO.MINCVX(), m)
end
```

The same functionality could be accomplished by setting the `EAGOParameter` field
`force_global_solve` to be true.
"""
optimize_hook!(t::ExtensionType, m::Optimizer) = nothing

"""
$(TYPEDSIGNATURES)

If the most recent upper problem returned a feasible result, and the upper
objective value is less than the previous best-known global upper bound,
set the most recent upper problem result to be the new global upper bound.
Update the `_feasible_solution_found`, `_first_solution_node`, 
`_global_upper_bound`, and `_continuous_solution` fields of the `GlobalOptimizer`
accordingly.
"""
function store_candidate_solution!(m::GlobalOptimizer)
    if m._upper_feasibility && (m._upper_objective_value < m._global_upper_bound)
        m._feasible_solution_found = true
        m._first_solution_node = m._maximum_node_id
        m._global_upper_bound = m._upper_objective_value
        @__dot__ m._continuous_solution = m._upper_solution
    end
    return
end

"""
$(TYPEDSIGNATURES)

If the previous best-known global lower bound is lower than the lowest
lower bound in the stack, set the global lower bound equal to the lowest
lower bound in the stack.
"""
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

Solves the branch-and-bound problem with the input `EAGO.GlobalOptimizer` object.

Pseudocode description of the algorithm, as implemented here:

-I) Prepare optimizers and stack for branch-and-bound

-II) While no reason to terminate the algorithm has occurred:

---II.A) Fathom nodes from the stack

---II.B) Select the new "current node" from the stack

---II.C) Perform preprocessing on current node

---II.D) If preprocessing result is feasible:

-----II.D.1) Solve lower problem for current node

-----II.D.2) If lower problem result is feasible and lower/upper bounds have not converged:

-------II.D.2.a) Solve upper problem for current node

-------II.D.2.b) Update the global upper bound if necessary

-------II.D.2.c) Perform postprocessing

-------II.D.2.d) If postprocessing result is feasible:

---------II.D.2.d.α) Branch and add the nodes back to the stack

---II.E) Update the global lower bound if necessary

---II.F) Update log information

-III) Set termination and result statuses

-IV) Print solution

"""
function global_solve!(m::GlobalOptimizer)

    # Set counts to 1
    m._iteration_count = 1
    m._node_count = 1

    # Prepare to run branch-and-bound
    parse_global!(m)
    presolve_global!(m)
    print_preamble!(m)

    # Run branch and bound; terminate when the stack is empty or when some
    # tolerance or limit is hit
    while !termination_check(m)

        # Fathom nodes from the stack, then pick a node and temporarily remove
        # it from the stack
        fathom!(m)
        node_selection!(m)
        print_node!(m)

        # Perform prepocessing and log the time
        m._last_preprocess_time += @elapsed preprocess!(m)

        # Continue if the node has not been proven infeasible
        if m._preprocess_feasibility

            # Solve the lower bounding problem and log the time
            m._last_lower_problem_time += @elapsed lower_problem!(m)
            print_results!(m, true)

            # Continue if lower problem is not infeasible and problem
            # problem has not yet converged
            if m._lower_feasibility && !convergence_check(m)

                # Solve the upper bounding problem and log the time
                m._last_upper_problem_time += @elapsed upper_problem!(m)
                print_results!(m, false)

                # Update the global upper bound if necessary
                store_candidate_solution!(m)

                # Perform post processing and log the time
                m._last_postprocessing_time += @elapsed postprocess!(m)

                # Continue if the node is not infeasible after postprocessing
                if m._postprocess_feasibility

                    # If the node is to be repeatedly evaluated, add it back
                    # onto the stack. If not, branch on the node and add the
                    # two new nodes to the stack
                    repeat_check(m) ? single_storage!(m) : branch_node!(m)
                end
            end
        else
            # "Disqualify" the node (TODO: Not strictly necessary, since it's
            # not added back to the stack? This node will simply be overwritten
            # in `node_selection!`)
            m._lower_objective_value = -Inf
            m._lower_feasibility = false
            m._upper_feasibility = false
        end

        # Update the global lower bound if necessary
        set_global_lower_bound!(m)

        # Adjust algorithm run times
        m._run_time = time() - m._start_time
        m._time_left = m._parameters.time_limit - m._run_time

        # Log and print information as needed and update the iteration counter
        log_iteration!(m)
        print_iteration!(m)
        m._iteration_count += 1
    end

    # Since the algorithm has terminated, convert EAGO's end status into
    # MOI.TerminationStatusCode and MOI.ResultStatusCode
    set_termination_status!(m)
    set_result_status!(m)

    # Print final information about the solution
    print_solution!(m)
end

"""
$(TYPEDSIGNATURES)

If global optimization was performed, much of the work happened within the
`_global_optimizer::GlobalOptimizer`. The `unpack_global_solution!` function
extracts results from the `GlobalOptimizer` and puts them in the correct fields
of the `Optimizer`.
"""
function unpack_global_solution!(m::Optimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    g = m._global_optimizer
    
    m._termination_status_code = g._termination_status_code
    m._result_status_code      = g._result_status_code

    m._run_time = g._run_time
    m._node_count = g._maximum_node_id

    # evaluate objective (so there isn't a small difference in f(x) and objective_value)
    # local solvers that solve to feasibility may result in a slightly lower than true solve...
    # TODO
    
    # Store objective value and objective bound 
    if _is_input_min(g)
        m._objective_bound = g._global_lower_bound
        m._objective_value = g._global_upper_bound
    else
        m._objective_bound = -g._global_lower_bound
        m._objective_value = -g._global_upper_bound
    end

    return
end

function optimize!(::MINCVX, m::Optimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    global_solve!(m._global_optimizer)
    unpack_global_solution!(m)
    return
end
