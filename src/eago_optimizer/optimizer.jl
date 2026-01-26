## Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/optimizer.jl
# Defines optimizer structure used by EAGO. Namely, ObjectiveType, ProblemType
# EAGOParameters, InputProblem, ParsedProblem, and Optimizer.
################################################################################

export Optimizer
"""
$(TYPEDEF)

The highest level optimizer object used by EAGO to solve problems during the optimization
routine. Additional options and temporary storage are located in the
`_global_optimizer::GlobalOptimizer{Q,S,T}` field. Parameters which are expected to
be constant over the entire solve are stored in the `_parameters::EAGOParameters` field. 
Some user-facing keywords not in the `EAGOParameters` field include:
- `relaxed_optimizer::MOI.AbstractOptimizer`: An instance of the optimizer used to solve 
    the relaxed subproblems (default = `GLPK.Optimizer()`). Located in `subsolver_block::SubSolvers{Q,S,T}`.
- `upper_optimizer::MOI.AbstractOptimizer`: Optimizer used to solve upper bounding problems 
    (default = `Ipopt.Optimizer()`). Located in `subsolver_block::SubSolvers{Q,S,T}`.
- `ext::ExtensionType`: Holds an instance of a subtype of `EAGO.ExtensionType`, used to define
    new custom subroutines (default = `DefaultExt()`). Located in `subsolver_block::SubSolvers{Q,S,T}`.
- `enable_optimize_hook::Bool`: Specifies that the user-defined `optimize_hook!` function should
    be called rather than use the standard EAGO optimization routines. Located in `Optimizer`
    and `_global_optimizer::GlobalOptimizer{Q,S,T}`.
- `obbt_variable_values::Vector{Bool}`: Variables to perform OBBT on (default: all variables in nonlinear
    expressions). Located in `_global_optimizer::GlobalOptimizer{Q,S,T}`.

Descriptions of all `Optimizer` fields available in extended help.

# Extended Help
$(TYPEDFIELDS)
"""
mutable struct Optimizer{Q,S,T} <: MOI.AbstractOptimizer
    "Holds definitions of the relaxed and upper optimizers, as well as any user-defined extension types"
    subsolver_block::SubSolvers{Q,S,T}
    "Specifies that the optimize_hook! function should be called rather than throw the
    problem to the standard routine"
    enable_optimize_hook::Bool
    "(Deprecated, use `subsolver_block` instead) Storage for custom extension types" 
    ext::Union{Nothing,T}
  
    "Information on any auxiliary variables"
    _auxiliary_variable_info::Union{Nothing,_AuxVarData}
    "Additional options and temporary storage for solving optimization problems"
    _global_optimizer::GlobalOptimizer{Q,S,T}
    "Expressions and constraints added to the EAGO model (not directly
    used for relaxations)"
    _input_problem::InputProblem
    "Expressions and problem descriptions that EAGO uses to formulate
    relaxed problems"
    _working_problem::ParsedProblem

    # Set as user-specified option
    "Parameters that do not change during a global solve"
    _parameters::EAGOParameters
    "Set of optimizer attributes"
    _optimizer_attributes_set::Vector{MOI.AbstractOptimizerAttribute}

    "The MathOptInterface-compliant completion status code"
    _termination_status_code::MOI.TerminationStatusCode
    "Value indicating the feasibility status of the result"
    _result_status_code::MOI.ResultStatusCode

    # Set constructor reset on empty! and to zero in initial parse! in parse.jl
    "Optimization run time"
    _run_time::Float64

    "The objective value of the primal solution"
    _objective_value::Float64
    "The best-known bound on the optimal objective value"
    _objective_bound::Float64
    "The gap between the upper and lower bound, relative to the bound with the larger magnitude"
    _relative_gap::Float64
    "The number of iterations the branch-and-bound algorithm has completed"
    _iteration_count::Int
    "The number of nodes in the stack"
    _node_count::Int
end
function Optimizer{Q,S,T}(sb::SubSolvers{Q,S,T}) where {Q,S,T}
    return Optimizer{Q,S,T}(sb, false, nothing, nothing, GlobalOptimizer{Q,S,T}(_subsolvers = sb, ext = _ext(sb)),
                     InputProblem(), ParsedProblem(), EAGOParameters(), MOI.AbstractOptimizerAttribute[],
                     MOI.OPTIMIZE_NOT_CALLED, MOI.OTHER_RESULT_STATUS,
                     0.0, -Inf, Inf, Inf, 0, 0)
end
function Optimizer(subsolver_block::SubSolvers{Q,S,T} = SubSolvers(); kwargs...) where {Q,S,T}
    if length(kwargs) > 0
        error("""Passing optimizer attributes as keyword arguments to `EAGO.Optimizer` is deprecated. 
                 Use MOI.set(model, MOI.RawParameter("key"), value) or 
                 JuMP.set_optimizer_attribute(model, "key", value) instead.""")
    end
    sb = SubSolvers{Incremental{Q},Incremental{S},T}(Incremental(subsolver_block.relaxed_optimizer), 
                                                     Incremental(subsolver_block.upper_optimizer),  
                                                     subsolver_block.ext)
    m = Optimizer{Incremental{Q},Incremental{S},T}(sb)
    m._global_optimizer = GlobalOptimizer{Incremental{Q},Incremental{S},T}(; _subsolvers = sb, ext = _ext(sb))
    return m
end

"""
    _constraints

Helper function which simplifies finding constraints of different types. See also: [`_constraint_primal`](@ref)
"""
_constraints(m::Optimizer, ::Type{VI}, ::Type{LT}) = m._input_problem._vi_leq_constraints
_constraints(m::Optimizer, ::Type{VI}, ::Type{GT}) = m._input_problem._vi_geq_constraints
_constraints(m::Optimizer, ::Type{VI}, ::Type{ET}) = m._input_problem._vi_eq_constraints
_constraints(m::Optimizer, ::Type{VI}, ::Type{IT}) = m._input_problem._vi_it_constraints
_constraints(m::Optimizer, ::Type{VI}, ::Type{ZO}) = m._input_problem._vi_zo_constraints
_constraints(m::Optimizer, ::Type{VI}, ::Type{INT}) = m._input_problem._vi_int_constraints

_constraints(m::Optimizer, ::Type{VI}, ::Type{MOI.Parameter{Float64}}) = m._input_problem._parameter_constraints

_constraints(m::Optimizer, ::Type{SAF}, ::Type{LT}) = _constraints(m._input_problem, SAF, LT)
_constraints(m::Optimizer, ::Type{SAF}, ::Type{GT}) = _constraints(m._input_problem, SAF, GT)
_constraints(m::Optimizer, ::Type{SAF}, ::Type{ET}) = _constraints(m._input_problem, SAF, ET)

_constraints(m::Optimizer, ::Type{SQF}, ::Type{LT}) = _constraints(m._input_problem, SQF, LT)
_constraints(m::Optimizer, ::Type{SQF}, ::Type{GT}) = _constraints(m._input_problem, SQF, GT)
_constraints(m::Optimizer, ::Type{SQF}, ::Type{ET}) = _constraints(m._input_problem, SQF, ET)

_ext(m::Optimizer) = _ext(m._global_optimizer)