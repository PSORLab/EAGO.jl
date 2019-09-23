"""
    initialize_log!

Creates the required storage variables for logging data at various iterations.
"""
function initialize_log!(m::Optimizer)

    log = Dict{Symbol,Any}()

    if x.log_on

        if x.log_subproblem_info
            log[:current_lower_bound] = Float64[]
            log[:current_upper_bound] = Float64[]

            log[:preprocessing_time] = Float64[]
            log[:lower_problem_time] = Float64[]
            log[:upper_problem_time] = Float64[]
            log[:postprocessing_time] = Float64[]

            log[:preprocessing_feas] = Bool[]
            log[:lower_problem_feas] = Bool[]
            log[:upper_problem_feas] = Bool[]
            log[:postprocessing_feas] = Bool[]
        end

        log[:global_lower_bound] = Float64[]
        log[:global_upper_bound] = Float64[]
        log[:node_count] = Float64[]
        log[:run_time] = Float64[]

        log[:parse_time] = 0.0
        log[:presolve_time] = 0.0

    end

    m.log = log
    return
end

"""
    log_iteration!

If 'logging_on' is true, the 'global_lower_bound', 'global_upper_bound',
'run_time', and 'node_count' are stored every 'log_interval'. If
'log_subproblem_info' then the lower bound, feasibility and run times of the
subproblems are logged every 'log_interval'.
"""
function log_iteration!(x::Optimizer)

    if x.log_on

        if (mod(x._iteration_count, x.log_interval) == 0 || x._iteration_count == 1)

            if x.log_subproblem_info
                push!(log[:current_lower_bound], x._lower_objective_value)
                push!(log[:current_upper_bound], x._upper_objective_value)

                push!(log[:preprocessing_time], x._last_preprocess_time)
                push!(log[:lower_problem_time], x._last_lower_problem_time)
                push!(log[:upper_problem_time], x._last_upper_problem_time)
                push!(log[:postprocessing_time], x._last_postprocessing_time)

                push!(log[:preprocess_feasibility], x._preprocess_feasibility)
                push!(log[:lower_problem_feasibility], x._lower_feasibility)
                push!(log[:upper_problem_feasibility], x._upper_feasibility)
                push!(log[:postprocess_feasibility], x._postprocess_feasibility)
            end

            push!(log[:global_lower_bound], x._global_lower_bound)
            push!(log[:global_upper_bound], x._global_upper_bound)
            push!(log[:run_time], x._run_time)
            push!(log[:node_count], x._node_count)

        end
    end
    return
end
