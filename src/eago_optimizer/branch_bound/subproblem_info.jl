abstract type SubProblemInfo <: Any end

"""
    LowerInfo

Storage type for lower problem results.
"""
mutable struct LowerInfo <: SubProblemInfo
    feasibility::Bool
    value::Float64
    solution::Vector{Float64}
    lower_variable_dual::Vector{Float64}
    upper_variable_dual::Vector{Float64}
end
LowerInfo() = LowerInfo(true,-Inf,Float64[],Float64[],Float64[])

"""
    UpperInfo

Storage type for upper problem results.
"""
mutable struct UpperInfo <: SubProblemInfo
    feasibility::Bool
    value::Float64
    solution::Vector{Float64}
end
UpperInfo() = UpperInfo(true,Inf,Float64[])

"""
    PreprocessInfo

Storage type for preprocessing results.
"""
mutable struct PreprocessInfo <: SubProblemInfo
    feasibility::Bool
end
PreprocessInfo() = PreprocessInfo(true)

"""
    PreprocessInfo

Storage type for postprocessing results.
"""
mutable struct PostprocessInfo <: SubProblemInfo
    feasibility::Bool
end
PostprocessInfo() = PostprocessInfo(true)
