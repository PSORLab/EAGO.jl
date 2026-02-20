# Copyright (c) 2018 Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Alireza Miraliakbar, Matthew Stuber, and the University of Connecticut (UConn)
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/optimize/nonconvex/log_iteration.jl
# Defines all routine to store information at each iteration.
################################################################################

"""
$(TYPEDSIGNATURES)

If `log_on` is true, the `global_lower_bound`, `global_upper_bound`,
`run_time`, and `node_count` are stored every `log_interval`. If
`log_subproblem_info` then the lower bound, feasibility and run times of the
subproblems are logged every `log_interval`.
"""
function log_iteration!(m::GlobalOptimizer)
    if m._parameters.log_on
        log = m._log
        if (mod(m._iteration_count, m._parameters.log_interval) == 0 || m._iteration_count == 1)
            if m._parameters.log_subproblem_info
                if _is_input_min(m)
                    push!(log.current_lower_bound, m._lower_objective_value)
                    push!(log.current_upper_bound, m._upper_objective_value)
                else
                    push!(log.current_lower_bound, -m._upper_objective_value)
                    push!(log.current_upper_bound, -m._lower_objective_value)
                end

                push!(log.preprocessing_time, m._last_preprocess_time)
                push!(log.lower_problem_time, m._last_lower_problem_time)
                push!(log.upper_problem_time, m._last_upper_problem_time)
                push!(log.postprocessing_time, m._last_postprocessing_time)

                push!(log.preprocess_feasibility, m._preprocess_feasibility)
                push!(log.lower_problem_feasibility, m._lower_feasibility)
                push!(log.upper_problem_feasibility, m._upper_feasibility)
                push!(log.postprocess_feasibility, m._postprocess_feasibility)
            end

            if _is_input_min(m)
                push!(log.global_lower_bound, m._global_lower_bound)
                push!(log.global_upper_bound, m._global_upper_bound)
            else
                push!(log.global_lower_bound, -m._global_upper_bound)
                push!(log.global_upper_bound, -m._global_lower_bound)
            end
            push!(log.run_time, m._run_time)
            push!(log.node_count, m._node_count)
        end
    end
    return
end
