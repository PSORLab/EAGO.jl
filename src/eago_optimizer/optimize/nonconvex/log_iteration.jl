# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/logging/log_iteration.jl
# Defines all routine to store information at each iteration.
#############################################################################

"""
$(TYPEDSIGNATURES)

If 'logging_on' is true, the 'global_lower_bound', 'global_upper_bound',
'run_time', and 'node_count' are stored every 'log_interval'. If
'log_subproblem_info' then the lower bound, feasibility and run times of the
subproblems are logged every 'log_interval'.
"""
function log_iteration!(x::GlobalOptimizer)

    if x._parameters.log_on
        log = x._log
        if (mod(x._iteration_count, x._parameters.log_interval) == 0 || x._iteration_count == 1)
            if x._parameters.log_subproblem_info
                if _is_input_min(m)
                    push!(log.current_lower_bound, x._lower_objective_value)
                    push!(log.current_upper_bound, x._upper_objective_value)
                else
                    push!(log.current_lower_bound, -x._upper_objective_value)
                    push!(log.current_upper_bound, -x._lower_objective_value)
                end

                push!(log.preprocessing_time, x._last_preprocess_time)
                push!(log.lower_problem_time, x._last_lower_problem_time)
                push!(log.upper_problem_time, x._last_upper_problem_time)
                push!(log.postprocessing_time, x._last_postprocessing_time)

                push!(log.preprocess_feasibility, x._preprocess_feasibility)
                push!(log.lower_problem_feasibility, x._lower_feasibility)
                push!(log.upper_problem_feasibility, x._upper_feasibility)
                push!(log.postprocess_feasibility, x._postprocess_feasibility)
            end

            if _is_input_min(x)
                push!(log.global_lower_bound, x._global_lower_bound)
                push!(log.global_upper_bound, x._global_upper_bound)
            else
                push!(log.global_lower_bound, -x._global_upper_bound)
                push!(log.global_upper_bound, -x._global_lower_bound)
            end
            push!(log.run_time, x._run_time)
            push!(log.node_count, x._node_count)
        end
    end
    return
end
