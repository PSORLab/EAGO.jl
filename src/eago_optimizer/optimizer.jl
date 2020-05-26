# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimizer.jl
# Defines optimizer structure used by EAGO, storage functions, and access
# functions.
#############################################################################

"""
$(TYPEDEF)

A structure used to store information related to the bounds assigned to each
variable.

$(TYPEDFIELDS)
"""
mutable struct VariableInfo
    "Is the variable integer valued?"
    is_integer::Bool
    "Lower bounds. May be -Inf."
    lower_bound::Float64
    "Boolean indicating whether finite lower bound exists."
    has_lower_bound::Bool
    "Upper bounds. May be Inf."
    upper_bound::Float64
    "Boolean indicating whether finite upper bound exists."
    has_upper_bound::Bool
    "Boolean indicating variable is fixed to a finite value."
    is_fixed::Bool
end
VariableInfo() = VariableInfo(false,-Inf, false, Inf, false, false)
lower_bound(x::VariableInfo) = x.lower_bound
upper_bound(x::VariableInfo) = x.upper_bound

"""
$(TYPEDEF)

An abstract type the subtypes of which are associated with functions method
overloaded for for new extensions. An instance of the `DefaultExt <:ExtensionType`
structure to the `Optimizer` in the `ext_type` field.
"""
abstract type ExtensionType end
struct DefaultExt <: ExtensionType end

@enum(ObjectiveType, UNSET, SINGLE_VARIABLE, SCALAR_AFFINE, SCALAR_QUADRATIC, NONLINEAR)

"""
$(TYPEDEF)

Storage for parameters that do not change during a global solve.

$(TYPEDFIELDS)
"""
Base.@kwdef struct EAGOParameters


    # Presolving options
    presolve_scrubber_flag::Bool = false
    "Create and use DAG representations of user-defined function (default = false)."
    presolve_to_JuMP_flag::Bool = false
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Apply the epigraph reformulation
    to the problem (default = false)."
    presolve_epigraph_flag::Bool = false
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Enable common subexpression
    elimination for DAG (default = false)."
    presolve_cse_flag::Bool = false
    "Rerranges the DAG using registered transformations (default = false)"
    presolve_flatten_flag::Bool = false

    # Conic reformulations
    "Attempt to bridge convex constraint to second order cone"
    conic_convert_quadratic::Bool = false

    # Iteration logging options
    "Turns logging on records global bounds, node count and run time. Additional
     options are available for recording information specific to subproblems (default = false)."
    log_on::Bool = false
    "Turns on logging of times and feasibility of subproblems (default = false)"
    log_subproblem_info::Bool = false
    "Log data every `log_interval` iterations (default = 1)."
    log_interval::Int64 = 1

    # Optimizer display options
    "The amount of information that should be printed to console while solving
    values range from 0 - 4: 0 is silent, 1 shows iteration summary statistics
    only, 2-4 show varying degrees of details about calculations within each
    iteration (default = 1)."
    verbosity::Int64 = 1
    "Display summary of iteration to console every `output_iterations` (default = 10)"
    output_iterations::Int64 = 1000
    "Display header for summary to console every `output_iterations` (default = 100)"
    header_iterations::Int64 = 10000

    # Node branching options
    "Convex coefficient used to select branch point. Branch point is given by
    `branch_cvx_factor*xmid + (1-branch_cvx_factor)*xsol` (default = 0.25)"
    branch_cvx_factor::Float64 = 0.25
    "Minimum distance from bound to have branch point normalized by width of
    dimension to branch on (default = 0.15)"
    branch_offset::Float64 = 0.15
    "Variables to branch on (default is all nonlinear)."
    branch_variable::Vector{Bool} = Bool[]
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Number of times repeat node
    processing priorto branching (default = 4)."
    branch_max_repetitions::Int64 = 4
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Volume ratio tolerance required
    to repeat processing the current node (default = 0.9)"
    branch_repetition_tol::Float64 = 0.9

    # Termination limits
    "Maximum number of nodes (default = 1E-7)"
    node_limit::Int64 = 1*10^7
    "Maximum CPU time in seconds (default = 1000)"
    time_limit::Float64 = 1000.0
    "Maximum number of iterations (default 3E6)"
    iteration_limit::Int64 = 3*10^6
    "Absolute tolerance for termination (default = 1E-3)"
    absolute_tolerance::Float64 = 1E-3
    "Relative tolerance for termination (default = 1E-3)"
    relative_tolerance::Float64 = 1E-3
    "Perform only a local solve of the problem (default = false)."
    local_solve_only::Bool = false
    feasible_local_continue::Bool = false

    # Options for constraint propagation
    "Depth in B&B tree above which constraint propagation should be disabled (default = 1000)"
    cp_depth::Int64 = 1000
    "Number of repetitions of forward-reverse passes to perform in constraint propagation (default = 3)"
    cp_repetitions::Int64 = 4
    "Disable constraint propagation if the ratio of new node volume to beginning node volume exceeds
    this number (default = 0.99)"
    cp_tolerance::Float64 = 0.99
    "Use only valid interval bounds during constraint propagation (default = false)"
    cp_interval_only::Bool = false

    # obbt options
    "Depth in B&B tree above which OBBT should be disabled (default = 6)"
    obbt_depth::Int64 = 6
    "Number of repetitions of OBBT to perform in preprocessing (default = 3)"
    obbt_repetitions::Int64 = 4
    "Turn aggresive OBBT on (default = false)"
    obbt_aggressive_on::Bool = true
    "Maximum iteration to perform aggresive OBBT (default = 2)"
    obbt_aggressive_max_iteration::Int64 = 2
    "Minimum dimension to perform aggresive OBBT (default = 2)"
    obbt_aggressive_min_dimension::Int64 = 2
    "Tolerance to consider bounds equal (default = 1E-9)"
    obbt_tolerance::Float64 = 1E-9

    # Options for linear bound tightening
    "Depth in B&B tree above which linear FBBT should be disabled (default = 1000)"
    lp_depth::Int64  = 100000
    "Number of repetitions of linear FBBT to perform in preprocessing (default = 3)"
    lp_repetitions::Int64  = 3

    # Options for quadratic bound tightening
    "Depth in B&B tree above which univariate quadratic FBBT should be disabled (default = -1)"
    quad_uni_depth::Int64 = -1
    "Number of repetitions of univariate quadratic FBBT to perform in preprocessing (default = 2)"
    quad_uni_repetitions::Int64 = 2
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Depth in B&B tree above which bivariate
    quadratic FBBT should be disabled (default = -1)"
    quad_bi_depth::Int64 = -1
    "Number of repetitions of bivariate quadratic FBBT to perform in preprocessing (default = 2)."
    quad_bi_repetitions::Int64 = 2

    # Duality-based bound tightening (DBBT) options
    "Depth in B&B tree above which duality-based bound tightening should be disabled (default = 1E10)"
    dbbt_depth::Int64 = 10^10
    "New bound is considered equal to the prior bound if within dbbt_tolerance (default = 1E-9)."
    dbbt_tolerance::Float64 = 1E-8

    # Subgradient tightening flag
    "Perform tightening of interval bounds using subgradients at each factor in
    each nonlinear tape during a forward-reverse pass (default = true)."
    subgrad_tighten::Bool = true
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Used to enable/disable subgradient
    tightening of interval bounds on the reverse pass (default = true)"
    subgrad_tighten_reverse::Bool = false

    # Tolerance to add cuts and max number of cuts
    cut_max_iterations::Int64 = 1
    "Convex coefficient used to select point for new added cuts. Branch point is
    given by `(1-cut_cvx)*xmid + cut_cvx*xsol` (default = 0.9)."
    cut_cvx::Float64 = 0.9
    "Add cut if the L1 distance from the prior cutting point to the new cutting
    point normalized by the box volume is greater than the tolerance (default = 0.05)."
    cut_tolerance::Float64 = 0.05
    "Adds an objective cut to the relaxed problem (default = true)."
    objective_cut_on::Bool = true
    "Lower tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_l::Float64 = 1E-8
    "Upper tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_u::Float64 = 1E8
    "Constant tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_b::Float64 = 1E9

    "Solve upper problem for every node with depth less than `upper_bounding_depth`
    and with a probabilityof (1/2)^(depth-upper_bounding_depth) otherwise (default = 6)"
    upper_bounding_depth::Int64 = 6

    # handling for domain violations
    domain_violation_ϵ::Float64 = 1E-9
end

export Optimizer
"""
$(TYPEDEF)

The main optimizer object used by EAGO to solve problems during the optimization
routine. The following commonly used options are described below and can be set
via keyword arguments in the JuMP/MOI model. The raw parameter interface however
is likely preferable. The Optimizer is organized in the following manner. Parameters
which are expected to be constant over the entire solve are stored in
`_parameters::EAGOParameters` field,

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Optimizer <: MOI.AbstractOptimizer

    # Options for optimality-based bound tightening
    "An instance of the optimizer used to solve the relaxed subproblems (default = GLPK.Optimizer())"
    relaxed_optimizer::MOI.AbstractOptimizer = GLPK.Optimizer()
    "Keyword arguments for the relaxed optimizer."
    relaxed_optimizer_kwargs::Base.Iterators.Pairs = Base.Iterators.Pairs(NamedTuple(),())

    "Variables to perform OBBT on (default: all variables in nonlinear expressions)."
    obbt_variable_values::Vector{Bool} = Bool[]

    # Upper bounding options
    upper_optimizer::MOI.AbstractOptimizer = Ipopt.Optimizer(max_iter = 3000, acceptable_tol = 1E30,
                                                             acceptable_iter = 300, constr_viol_tol = 1E-8,
                                                             acceptable_constr_viol_tol = 1E-8, print_level = 0)
    upper_factory::JuMP.OptimizerFactory = with_optimizer(Ipopt.Optimizer, max_iter = 3000, acceptable_tol = 1E30,
                                                          acceptable_iter = 300, constr_viol_tol = 1E-8,
                                                          acceptable_constr_viol_tol = 1E-8, print_level = 0)

    # Extensions
    "Specifies that the optimize_hook! function should be called rather than
    throw the problem to the standard B&B routine (default = false)."
    enable_optimize_hook::Bool = false
    "Holds additional storage needed for constructing extensions to EAGO
    (default = Dict{Symbol,Any})."
    ext::Dict{Symbol, Any} = Dict{Symbol,Any}()
    "Holds an instance of a subtype of `EAGO.ExtensionType` used to define
    new custom subroutines (default = DefaultExt())."
    ext_type::ExtensionType = DefaultExt()

    _parameters::EAGOParameters = EAGOParameters()
    _input_problem::InputProblem = InputProblem()
    _working_problem::ParsedProblem = ParsedProblem()

    _termination_status_code::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _result_status_code::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS

    _current_node::NodeBB = NodeBB()
    _current_xref::Vector{Float64} = Float64[]
    _stack::BinaryMinMaxHeap{NodeBB} = BinaryMinMaxHeap{NodeBB}()

    _user_branch_variables::Bool = false
    _fixed_variable::Vector{Bool} = Bool[]
    _continuous_solution::Vector{Float64} = Float64[]
    _upper_variables::Vector{VI} =  VI[]

    _relaxed_variable::Vector{SV} = SV[]
    _relaxed_variable_index::Vector{VI} = VI[]
    _relaxed_variable_et::Vector{CI{SV, ET}} = CI{SV, ET}[]
    _relaxed_variable_lt::Vector{CI{SV, LT}} = CI{SV, LT}[]
    _relaxed_variable_gt::Vector{CI{SV, GT}} = CI{SV, GT}[]
    _relaxed_variable_et_indx::Vector{Int64} = Int64[]
    _relaxed_variable_lt_indx::Vector{Int64} = Int64[]
    _relaxed_variable_gt_indx::Vector{Int64} = Int64[]

    _preprocess_feasibility::Bool = true
    _preprocess_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _preprocess_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED

    _lower_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _lower_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _lower_feasibility::Bool = false
    _lower_objective_value::Float64 = -Inf
    _lower_solution::Vector{Float64} = Float64[]
    _lower_lvd::Vector{Float64} = Float64[]
    _lower_uvd::Vector{Float64} = Float64[]

    _cut_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _cut_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _cut_solution::Vector{Float64} = Float64[]
    _cut_objective_value::Float64 = -Inf
    _cut_feasibility::Bool = false

    _upper_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _upper_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _upper_feasibility::Bool = false
    _upper_objective_value::Float64 = Inf
    _upper_solution::Vector{Float64} = Float64[]

    _postprocess_feasibility::Bool = false

    _start_time::Float64 = 0.0
    _run_time::Float64 = 0.0
    _time_left::Float64 = 1000.0
    _parse_time::Float64 = 0.0
    _presolve_time::Float64 = 0.0
    _last_preprocess_time::Float64 = 0.0
    _last_lower_problem_time::Float64 = 0.0
    _last_upper_problem_time::Float64 = 0.0
    _last_postprocessing_time::Float64 = 0.0

    _objective_cut_ci_sv::CI{SV,LT} = CI{SV,LT}(-1.0)
    _objective_cut_ci_saf::Vector{CI{SAF,LT}} = CI{SAF,LT}[]

    _global_lower_bound::Float64 = -Inf
    _global_upper_bound::Float64 = Inf
    _maximum_node_id::Int64 = 0
    _iteration_count::Int64 = 0
    _node_count::Int64 = 0

    # Storage for output
    _solution_value::Float64 = 0.0
    _feasible_solution_found::Bool = false
    _first_solution_node::Int64 = -1
    _objective_value::Float64 = -Inf
    _best_upper_value::Float64 = Inf

    # Optimality-Based Bound Tightening (OBBT) Options
    _obbt_working_lower_index::Vector{Bool} = Bool[]
    _obbt_working_upper_index::Vector{Bool} = Bool[]
    _lower_indx_diff::Vector{Bool} = Bool[]
    _upper_indx_diff::Vector{Bool} = Bool[]
    _old_low_index::Vector{Bool} = Bool[]
    _old_upp_index::Vector{Bool} = Bool[]
    _new_low_index::Vector{Bool} = Bool[]
    _new_upp_index::Vector{Bool} = Bool[]
    _obbt_variables::Vector{VI} = VI[]
    _obbt_performed_flag::Bool = false

    # Feasibility-Based Bound Tightening Options
    _cp_improvement::Float64 = 0.0
    _cp_evaluation_reverse::Bool = false

    _cut_iterations::Int64 = 0
    _cut_add_flag::Bool = false

    # Options for Repetition (If DBBT Performed Well)
    _node_repetitions::Int64 = 0
    _initial_volume::Float64 = 0.0
    _final_volume::Float64 = 0.0

    # Log
    _log::Log = Log()

    _relaxed_evaluator::Evaluator = Evaluator{1,NS}()
    _relaxed_constraint_bounds::Vector{MOI.NLPBoundsPair} = Vector{MOI.NLPBoundsPair}[]
    _relaxed_eval_has_objective::Bool = false
end

function MOI.empty!(m::Optimizer)
    m = Optimizer()
    nothing
end

function MOI.is_empty(m::Optimizer)

    flag = true
    flag &= isempty(m._variable_info)
    flag &= m._input_problem._optimization_sense === MOI.MIN_SENSE
    flag &= m._termination_status_code === MOI.OPTIMIZE_NOT_CALLED

    return flag
end

##### Utilities for checking that JuMP model contains variables used in expression
function check_inbounds!(m::Optimizer, vi::VI)
    if !(1 <= vi.value <= m._variable_number)
        error("Invalid variable index $vi. ($(m._variable_numbe) variables in the model.)")
    end
    return
end
check_inbounds!(m::Optimizer, var::SV) = check_inbounds!(m, var.variable)

function check_inbounds!(m::Optimizer, aff::SAF)
    for term in aff.terms
        check_inbounds!(m, term.variable_index)
    end
    return
end
function check_inbounds!(m::Optimizer, quad::SQF)
    for term in quad.affine_terms
        check_inbounds!(m, term.variable_index)
    end
    for term in quad.quadratic_terms
        check_inbounds!(m, term.variable_index_1)
        check_inbounds!(m, term.variable_index_2)
    end
    return
end
#=
function check_inbounds!(m::Optimizer, vov::VECOFVAR)
    for vi in vov.variables
        check_inbounds!(m, vi)
    end
    return
end
=#

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; copy_names = false)
    return MOI.Utilities.default_copy_to(model, src, copy_names)
end

function MOI.set(m::Optimizer, ::MOI.Silent, value)
     m.verbosity = 0
     m.log_on = false
     return
end

#=
function MOI.set(m::Optimizer, ::MOI.TimeLimitSec, value::Nothing)
    m.time_limit = Inf
    return
end

function MOI.set(m::Optimizer, ::MOI.TimeLimitSec, value::Float64)
    m.time_limit = value
    return
end
=#

MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices) = [MOI.VariableIndex(i) for i = 1:length(m._variable_info)]
function MOI.get(m::Optimizer, ::MOI.ObjectiveValue)
    mult = 1.0
    if m._input_problem._optimization_sense === MOI.MAX_SENSE
        mult *= -1.0
    end
    return mult*m._objective_value
end
MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = m._variable_number
function MOI.get(m::Optimizer, ::MOI.ObjectiveBound)
    if m._input_problem._optimization_sense === MOI.MAX_SENSE
        bound = -m._global_lower_bound
    else
        bound = m._global_upper_bound
    end
    return bound
end
function MOI.get(m::Optimizer, ::MOI.RelativeGap)
    LBD = m._global_lower_bound
    UBD = m._global_upper_bound
    if m._input_problem._optimization_sense === MOI.MAX_SENSE
        gap = abs(UBD - LBD)/abs(LBD)
    else
        gap = abs(UBD - LBD)/abs(UBD)
    end
    return gap
end
MOI.get(m::Optimizer, ::MOI.SolverName) = "EAGO: Easy Advanced Global Optimization"
MOI.get(m::Optimizer, ::MOI.TerminationStatus) = m._termination_status_code
MOI.get(m::Optimizer, ::MOI.PrimalStatus) = m._result_status_code
MOI.get(m::Optimizer, ::MOI.SolveTime) = m._run_time
MOI.get(m::Optimizer, ::MOI.NodeCount) = m._maximum_node_id
MOI.get(m::Optimizer, ::MOI.ResultCount) = (m._result_status_code === MOI.FEASIBLE_POINT) ? 1 : 0

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    check_inbounds!(model, vi)
    return model._continuous_solution[vi.value]
end

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.set(m::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    if nlp_data.has_objective
        m._input_problem._objective_type = NONLINEAR
    end
    m._input_problem._nlp_data = nlp_data
    return
end

##### Support, set, and evaluate objective functions
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{F}) where {F <: Union{SV,SAF,SQF}} = true
function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction{SV}, func::SV)
    check_inbounds!(m, func)
    m._input_problem._objective_sv = func
    m._input_problem._objective_type = SINGLE_VARIABLE
    return
end
function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction{SAF}, func::SAF)
    check_inbounds!(m, func)
    m._input_problem._objective_saf = func
    m._input_problem._objective_type = SCALAR_AFFINE
    return
end
function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction{SQF}, func::SQF)
    check_inbounds!(m, func)
    m._input_problem._objective_sqf = func
    m._input_problem._objective_type = SCALAR_QUADRATIC
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    m._input_problem._optimization_sense = sense
    return
end
