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
