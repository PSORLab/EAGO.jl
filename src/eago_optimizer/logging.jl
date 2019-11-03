"""
    Log

A structure used to store information on the history of the solution procedure
for generating convergence plots and other analysis.

- `current_lower_bound::Vector{Float64}`: Storage for lower bound calculated for current node.
- `current_upper_bound::Vector{Float64}`: Storage for upper bound calculated for current node.
- `preprocessing_time::Vector{Float64}`: Storage for preprocessing time of each iteration.
- `lower_problem_time::Vector{Float64}`: Storage for lower bounding time of each iteration.
- `upper_problem_time::Vector{Float64}`: Storage for upper bounding time of each iteration.
- `postprocessing_time::Vector{Float64}`: Storage for postprocessing time of each iteration.
- `preprocessing_feas::Vector{Bool}`: Storage for preprocessing feasibility of each iteration.
- `lower_problem_feas::Vector{Bool}`: Storage for lower bounding feasibility of each iteration.
- `upper_problem_feas::Vector{Bool}`: Storage for upper bounding feasibility of each iteration.
- `postprocessing_feas::Vector{Bool}`: Storage for postprocessing feasibility of each iteration.
- `global_lower_bound::Vector{Float64}`: Storage for best (global) lower bound at each iteration.
- `global_upper_bound::Vector{Float64}`: Storage for best (global) upper bound at each iteration.
"""
struct Log
    current_lower_bound::Vector{Float64}
    current_upper_bound::Vector{Float64}
    preprocessing_time::Vector{Float64}
    lower_problem_time::Vector{Float64}
    upper_problem_time::Vector{Float64}
    postprocessing_time::Vector{Float64}
    preprocessing_feas::Vector{Bool}
    lower_problem_feas::Vector{Bool}
    upper_problem_feas::Vector{Bool}
    postprocessing_feas::Vector{Bool}
    global_lower_bound::Vector{Float64}
    global_upper_bound::Vector{Float64}
    node_count::Vector{Int64}
    run_time::Vector{Float64}
end

Log() = Log(Float64[], Float64[],
            Float64[], Float64[],
            Float64[], Float64[],
            Bool[],    Bool[],
            Bool[],    Bool[],
            Float64[], Float64[],
            Int64[], Float64[])
