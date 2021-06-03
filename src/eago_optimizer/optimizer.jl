# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimizer.jl
# Defines optimizer structure used by EAGO. Namely, ObjectiveType, ProblemType
# EAGOParameters, InputProblem, ParsedProblem, and Optimizer.
#############################################################################

export Optimizer
"""
$(TYPEDEF)

The main optimizer object used by EAGO to solve problems during the optimization
routine. The following commonly used options are described below and can be set
via keyword arguments in the JuMP/MOI model. The raw parameter interface however
is likely preferable. The Optimizer is organized in the following manner. Parameters
which are expected to be constant over the entire solve are stored in
`_parameters::EAGOParameters` field. User-facing keywords not in EAGOParameters field:
- `relaxed_optimizer::MOI.AbstractOptimizer`: An instance of the optimizer used to solve the relaxed subproblems (default = GLPK.Optimizer())
- `obbt_variable_values::Vector{Bool}`: Variables to perform OBBT on (default: all variables in nonlinear expressions).
- `upper_optimizer::MOI.AbstractOptimizer`: Optimizer used to solve upper bounding problems. (default = Ipopt.Optimizer)
- `enable_optimize_hook::Bool`: Specifies that the optimize_hook! function should be called rather than throw the problem to the standard B&B routine (default = false).
- `ext::Dict{Symbol, Any}`: Holds additional storage needed for constructing extensions to EAGO (default = Dict{Symbol,Any}).
- `ext_type::ExtensionType`: Holds an instance of a subtype of `EAGO.ExtensionType` used to define new custom subroutines (default = DefaultExt()).
"""
mutable struct Optimizer{Q,S,T} <: MOI.AbstractOptimizer

    subsolver_block::SubSolvers{Q,S,T}
    enable_optimize_hook::Bool
    ext::Dict{Symbol, Any}
  
    _global_optimizer::GlobalOptimizer{Q,S,T}
    _input_problem::InputProblem
    _working_problem::ParsedProblem

       # set as user-specified option
    _parameters::EAGOParameters

    _termination_status_code::MOI.TerminationStatusCode
    _result_status_code::MOI.ResultStatusCode

    # set constructor reset on empty! and  to zero in initial parse! in parse.jl
    _run_time::Float64

    _objective_value::Float64
    _objective_bound::Float64
    _relative_gap::Float64
    _iteration_count::Int
    _node_count::Int
end
function Optimizer{Q,S,T}(; subsolver_block::SubSolvers{Q,S,T} = SubSolvers{Q,S,T}(),
                            enable_optimize_hook::Bool = false,
                            ext::Dict{Symbol, Any} = Dict{Symbol,Any}(),
                            _global_optimizer::GlobalOptimizer{Q,S,T} = GlobalOptimizer{Q,S,T}(),
                            _input_problem::InputProblem = InputProblem(),
                            _working_problem::ParsedProblem = ParsedProblem(),
                            _parameters::EAGOParameters = EAGOParameters(),
                            _termination_status_code::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED,
                            _result_status_code::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS,
                            _run_time::Float64   = 0.0,
                            _objective_value::Float64 = -Inf,
                            _objective_bound::Float64 =  Inf,
                            _relative_gap::Float64    = Inf,
                            _iteration_count::Int     = 0,
                            _node_count::Int          = 0) where {Q,S,T}

    Optimizer{Q,S,T}(subsolver_block,
                     enable_optimize_hook,
                     ext, 
                     _global_optimizer,
                     _input_problem,
                     _working_problem,
                     _parameters,
                     _termination_status_code,
                     _result_status_code,
                     _run_time,
                     _objective_value,
                     _objective_bound,
                     _relative_gap ,
                     _iteration_count,
                     _node_count)
end

function Optimizer(; subsolver_block::SubSolvers{Q,S,T} = SubSolvers(), kwargs...) where {Q,S,T}
    sb = SubSolvers{Incremental{Q},Incremental{S},T}(Incremental(subsolver_block.relaxed_optimizer), 
                                                     Incremental(subsolver_block.upper_optimizer), 
                                                     subsolver_block.ext_typ )
    gopt = GlobalOptimizer{Incremental{Q},Incremental{S},T}(; _subsolvers = sb)
    return Optimizer{Incremental{Q},Incremental{S},T}(; subsolver_block = sb, _global_optimizer = gopt)
end