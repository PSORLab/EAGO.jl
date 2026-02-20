# Copyright (c) 2018 Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Alireza Miraliakbar, Matthew Stuber, and the University of Connecticut (UConn)
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/types/log.jl
# Defines all type used to store solution information at given iterations.
################################################################################

"""
    $(TYPEDEF)

A structure used to store information on the history of the solution procedure
for generating convergence plots and other analyses.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Log
    "Storage for lower bound calculated for current node."
    current_lower_bound::Vector{Float64} = Float64[]
    "Storage for upper bound calculated for current node."
    current_upper_bound::Vector{Float64} = Float64[]
    "Storage for preprocessing time of each iteration."
    preprocessing_time::Vector{Float64} = Float64[]
    "Storage for lower bounding time of each iteration."
    lower_problem_time::Vector{Float64} = Float64[]
    "Storage for upper bounding time of each iteration."
    upper_problem_time::Vector{Float64} = Float64[]
    "Storage for postprocessing time of each iteration."
    postprocessing_time::Vector{Float64} = Float64[]
    "Storage for preprocessing feasibility of each iteration."
    preprocessing_feas::Vector{Bool} = Bool[]
    "Storage for lower bounding feasibility of each iteration."
    lower_problem_feas::Vector{Bool} = Bool[]
    "Storage for upper bounding feasibility of each iteration."
    upper_problem_feas::Vector{Bool} = Bool[]
    "Storage for postprocessing feasibility of each iteration."
    postprocessing_feas::Vector{Bool} = Bool[]
    "Storage for best (global) lower bound at each iteration."
    global_lower_bound::Vector{Float64} = Float64[]
    "Storage for best (global) upper bound at each iteration."
    global_upper_bound::Vector{Float64} = Float64[]
    "Number of nodes at each iteration."
    node_count::Vector{Int64} = Int64[]
    "Run time at each iteration."
    run_time::Vector{Float64} = Float64[]
end

function Base.isempty(x::Log)
    for field in fieldnames(Log)
        if !isempty(getfield(x, field))
            return false
        end
    end
    return true
end
