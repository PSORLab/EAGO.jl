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

@enum(BranchCost, BC_INFEASIBLE, BC_INTERVAL, BC_INTERVAL_REV, BC_INTERVAL_LP, BC_INTERVAL_LP_REV)

Base.@kwdef mutable struct BranchCostStorage{T<:Real}
    cost::BranchCost         = BC_INTERVAL
    𝛹n::Vector{T}           = T[]
    𝛹p::Vector{T}           = T[]
    δn::Vector{T}            = T[]
    δp::Vector{T}            = T[]
    ηn::Vector{T}            = T[]
    ηp::Vector{T}            = T[]
    μ1::T                    = 0.1
    μ2::T                    = 1.3
    μ3::T                    = 0.8
    β::T                     = 0.05
    μ_score::T               = 0.15
end
function initialize!(d::BranchCostStorage{T}, n::Int) where T <:AbstractFloat
    append!(d.𝛹n, ones(T,n));  append!(d.𝛹p, ones(T,n))
    append!(d.δn, zeros(T,n));  append!(d.δp, zeros(T,n))
    append!(d.ηn, zeros(T,n));  append!(d.ηp, zeros(T,n))
    return
end

@enum(ProblemType, UNCLASSIFIED, LP, MILP, SOCP, MISOCP, DIFF_CVX, MINCVX)
@enum(GlobalEndState, GS_OPTIMAL, GS_INFEASIBLE, GS_NODE_LIMIT,
                      GS_ITERATION_LIMIT, GS_RELATIVE_TOL,
                      GS_ABSOLUTE_TOL, GS_TIME_LIMIT, GS_UNSET)

export EAGOParameters
"""
$(TYPEDEF)

Storage for parameters that do not change during a global solve.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct EAGOParameters

    # Presolving options
    "Should EAGO attempt to remove type assert issues for user-defined functions (default = false)"
    presolve_scrubber_flag::Bool = false
    "Create and use DAG representations of user-defined function (default = false)."
    presolve_to_JuMP_flag::Bool = false
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
    log_interval::Int = 1

    # Optimizer display options
    "The amount of information that should be printed to console while solving
    values range from 0 - 4: 0 is silent, 1 shows iteration summary statistics
    only, 2-4 show varying degrees of details about calculations within each
    iteration (default = 1)."
    verbosity::Int = 1
    "Display summary of iteration to console every `output_iterations` (default = 10)"
    output_iterations::Int = 1000
    "Display header for summary to console every `output_iterations` (default = 100)"
    header_iterations::Int = 10000

    # Node branching options
    "Convex coefficient used to select branch point. Branch point is given by
    `branch_cvx_factor*xmid + (1-branch_cvx_factor)*xsol` (default = 0.25)"
    branch_cvx_factor::Float64 = 0.25
    "Minimum distance from bound to have branch point normalized by width of
    dimension to branch on (default = 0.15)"
    branch_offset::Float64 = 0.15
    "Indicates that pseudocost branching should be used"
    branch_pseudocost_on::Bool = false
    "Variables to branch on (default is all nonlinear)."
    branch_variable::Vector{Bool} = Bool[]
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Number of times repeat node
    processing priorto branching (default = 4)."
    branch_max_repetitions::Int = 4
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Volume ratio tolerance required
    to repeat processing the current node (default = 0.9)"
    branch_repetition_tol::Float64 = 0.9

    # Termination limits
    "Maximum number of nodes (default = 1E-7)"
    node_limit::Int = 1*10^7
    "Maximum CPU time in seconds (default = 1000)"
    time_limit::Float64 = 1000.0
    "Maximum number of iterations (default 3E6)"
    iteration_limit::Int = 3*10^6
    "Absolute tolerance for termination (default = 1E-3)"
    absolute_tolerance::Float64 = 1E-3
    "Relative tolerance for termination (default = 1E-3)"
    relative_tolerance::Float64 = 1E-3
    "Absolute constraint feasibility tolerance"
    absolute_constraint_feas_tolerance::Float64 = 1E-6
    "Perform only a local solve of the problem (default = false)."
    local_solve_only::Bool = false
    "[TO BE REMOVED] Flag stops B&B loop if feasible point found."
    feasible_local_continue::Bool = false

    # Options for constraint propagation
    "Depth in B&B tree above which constraint propagation should be disabled (default = 1000)"
    cp_depth::Int = 20
    "Number of times to repeat forward-reverse pass routine (default = 3)"
    cp_repetitions::Int = 3
    "Disable constraint propagation if the ratio of new node volume to beginning node volume exceeds
    this number (default = 0.99)"
    cp_tolerance::Float64 = 0.99
    "Use only valid interval bounds during constraint propagation (default = false)"
    cp_interval_only::Bool = false

    # obbt options
    "Depth in B&B tree above which OBBT should be disabled (default = 6)"
    obbt_depth::Int = 4
    "Number of repetitions of OBBT to perform in preprocessing (default = 3)"
    obbt_repetitions::Int = 20
    "Turn aggresive OBBT on (default = false)"
    obbt_aggressive_on::Bool = true
    "Maximum iteration to perform aggresive OBBT (default = 2)"
    obbt_aggressive_max_iteration::Int = 2
    "Minimum dimension to perform aggresive OBBT (default = 2)"
    obbt_aggressive_min_dimension::Int = 2
    "Tolerance to consider bounds equal (default = 1E-9)"
    obbt_tolerance::Float64 = 1E-9

    # Options for linear bound tightening
    "Depth in B&B tree above which linear FBBT should be disabled (default = 1000)"
    fbbt_lp_depth::Int  = 1000
    "Number of repetitions of linear FBBT to perform in preprocessing (default = 3)"
    fbbt_lp_repetitions::Int  = 3

    # Options for quadratic bound tightening
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Depth in B&B tree above which univariate quadratic FBBT should be disabled (default = -1)"
    quad_uni_depth::Int = -1
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Number of repetitions of univariate quadratic FBBT to perform in preprocessing (default = 2)"
    quad_uni_repetitions::Int = 2
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Depth in B&B tree above which bivariate
    quadratic FBBT should be disabled (default = -1)"
    quad_bi_depth::Int = -1
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Number of repetitions of bivariate quadratic FBBT to perform in preprocessing (default = 2)."
    quad_bi_repetitions::Int = 2

    # Duality-based bound tightening (DBBT) options
    "Depth in B&B tree above which duality-based bound tightening should be disabled (default = 1E10)"
    dbbt_depth::Int = 10^10
    "New bound is considered equal to the prior bound if within dbbt_tolerance (default = 1E-9)."
    dbbt_tolerance::Float64 = 1E-8

    # Subgradient tightening flag
    "Relax Tag used to specify type of McCormick operator"
    relax_tag::RelaxTag = NS()
    "Perform tightening of interval bounds using subgradients at each factor in
    each nonlinear tape during a forward pass (default = true)."
    subgrad_tighten::Bool = true
    "Perform tightening of interval bounds using subgradients at each factor in
    each nonlinear tape during a reverse pass (default = false)."
    reverse_subgrad_tighten::Bool = false
    "Outer round computed subgradient bounds by this amount"
    subgrad_tol::Float64 = 1E-10

    # Tolerance to add cuts and max number of cuts
    "Minimum number of cuts at each node to attempt (unsafe cuts not necessarily added)"
    cut_min_iterations::Int = 1
    "Maximum number of cuts at each node to attempt"
    cut_max_iterations::Int = 5
    "Absolute tolerance checked for continuing cut"
    cut_tolerance_abs::Float64 = 1E-3
    "Relative tolerance checked for continuing cut"
    cut_tolerance_rel::Float64 = 1E-3

    "Use tolerances to determine safe cuts in a Khajavirad 2018 manner"
    cut_safe_on::Bool = true
    "Lower tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_l::Float64 = 1E-8
    "Upper tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_u::Float64 = 1E8
    "Constant tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_b::Float64 = 1E9

    "Solve upper problem for every node with depth less than `upper_bounding_depth`
    and with a probabilityof (1/2)^(depth-upper_bounding_depth) otherwise (default = 6)"
    upper_bounding_depth::Int = 8

    # handling for domain violations
    "Amount about a domain violation to ignore when propagating bounds."
    domain_violation_guard_on::Bool = false
    "Amount about a domain violation to ignore when propagating bounds."
    domain_violation_ϵ::Float64 = 1E-9
end

"""
$(TYPEDEF)

A structure used to hold objectives and constraints added to EAGO model.
The constraints generally aren't used for relaxations.
"""
Base.@kwdef mutable struct InputProblem

    # variables (set by MOI.add_variable in variables.jl)
    _variable_info::Vector{VariableInfo{Float64}} = VariableInfo{Float64}[]
    _variable_count::Int = 0

    # last constraint index added
    _last_constraint_index::Int = 0

    # linear constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _linear_leq_constraints::Vector{Tuple{SAF, LT}} = Tuple{SAF, LT}[]
    _linear_geq_constraints::Vector{Tuple{SAF, GT}} = Tuple{SAF, GT}[]
    _linear_eq_constraints::Vector{Tuple{SAF, ET}} = Tuple{SAF, ET}[]

    # quadratic constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _quadratic_leq_constraints::Vector{Tuple{SQF, LT}} = Tuple{SQF, LT}[]
    _quadratic_geq_constraints::Vector{Tuple{SQF, GT}} = Tuple{SQF, GT}[]
    _quadratic_eq_constraints::Vector{Tuple{SQF, ET}} = Tuple{SQF, ET}[]

    # conic constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _conic_second_order::Vector{Tuple{VECOFVAR, MOI.SecondOrderCone}} = Tuple{VECOFVAR, MOI.SecondOrderCone}[]

    # nonlinear constraint storage
    _objective::Union{SV,SAF,SQF,Nothing} = nothing

    # nlp constraints (set by MOI.set(m, ::NLPBlockData...) in optimizer.jl)
    _nlp_data::MOI.NLPBlockData = empty_nlp_data()

    # objective sense information (set by MOI.set(m, ::ObjectiveSense...) in optimizer.jl)
    _optimization_sense::MOI.OptimizationSense = MOI.MIN_SENSE
end

function Base.isempty(x::InputProblem)

    is_empty_flag = true
    new_input_problem = InputProblem()

    for field in fieldnames(InputProblem)
        field_value = getfield(x, field)
        if field_value isa Array
            if !isempty(field_value)
                is_empty_flag = false
                break
            end
        elseif field_value isa Number
            if getfield(new_input_problem, field) !== field_value
                is_empty_flag = false
                break
            end
        end
    end

    is_empty_flag &= x._nlp_data.evaluator isa EmptyNLPEvaluator
    is_empty_flag &= !x._nlp_data.has_objective
    is_empty_flag &= isempty(x._nlp_data.constraint_bounds)
    is_empty_flag &= x._objective == nothing

    return is_empty_flag
end

"""
$(TYPEDEF)

A structure used to expressions and problem descriptions EAGO uses to formulate
relaxed problems.
"""
Base.@kwdef mutable struct ParsedProblem

    # Problem classification (set in parse_classify_problem!)
    _problem_type::ProblemType = UNCLASSIFIED

    "_objective_saf stores the objective and is used for constructing linear affine cuts"
    _objective_saf::SAF = SAF(SAT[], 0.0)
    _objective::Union{SV,AffineFunctionIneq,BufferedQuadraticIneq,BufferedNonlinearFunction,Nothing} = nothing

    # objective sense information (set by convert_to_min in parse.jl)
    _optimization_sense::MOI.OptimizationSense = MOI.MIN_SENSE

    # non-single variable constraints (set in initial_parse)
    _saf_leq::Vector{AffineFunctionIneq} = AffineFunctionIneq[]
    _saf_eq::Vector{AffineFunctionEq} = AffineFunctionEq[]
    _sqf_leq::Vector{BufferedQuadraticIneq} = BufferedQuadraticIneq[]
    _sqf_eq::Vector{BufferedQuadraticEq} = BufferedQuadraticEq[]
    _conic_second_order::Vector{BufferedSOC} = BufferedSOC[]

    # nlp constraints (set in initial_parse)
    _nlp_data::MOI.NLPBlockData = empty_nlp_data()

    # storage for nonlinear functions
    _nonlinear_constr::Vector{BufferedNonlinearFunction} = BufferedNonlinearFunction[]

    # nonlinear evaluator
    _relaxed_evaluator = Evaluator()

    # variables (set in initial_parse)
    _variable_info::Vector{VariableInfo{Float64}} = VariableInfo{Float64}[]
    _variable_count::Int = 0
end

function Base.isempty(x::ParsedProblem)

    is_empty_flag = true

    new_input_problem = ParsedProblem()
    for field in fieldnames(ParsedProblem)

        field_value = getfield(x, field)
        if field_value isa Array
            if !isempty(field_value)
                is_empty_flag = false
                break
            end

        elseif field_value isa Number
            if getfield(new_input_problem, field) !== field_value
                is_empty_flag = false
                break
            end
        end
    end

    is_empty_flag &= x._nlp_data.evaluator isa EmptyNLPEvaluator
    is_empty_flag &= !x._nlp_data.has_objective
    is_empty_flag &= isempty(x._nlp_data.constraint_bounds)

    is_empty_flag &= isempty(x._objective_saf.terms)
    is_empty_flag &= x._objective_saf.constant === 0.0
    is_empty_flag &= x._objective === nothing

    return is_empty_flag
end


function default_nlp_solver()

    upper_optimizer = Ipopt.Optimizer()

    MOI.set(upper_optimizer, MOI.RawParameter("max_iter"),5000)
    MOI.set(upper_optimizer, MOI.RawParameter("acceptable_tol"), 1E30)
    MOI.set(upper_optimizer, MOI.RawParameter("acceptable_iter"), 500)
    MOI.set(upper_optimizer, MOI.RawParameter("constr_viol_tol"), 0.000001)
    MOI.set(upper_optimizer, MOI.RawParameter("acceptable_compl_inf_tol"), 0.000001)
    MOI.set(upper_optimizer, MOI.RawParameter("acceptable_dual_inf_tol"), 1.0)
    MOI.set(upper_optimizer, MOI.RawParameter("acceptable_constr_viol_tol"), 0.000001)
    MOI.set(upper_optimizer, MOI.RawParameter("print_level"), 0)

    return upper_optimizer
end

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
Base.@kwdef mutable struct Optimizer <: MOI.AbstractOptimizer

    # Options for optimality-based bound tightening
    # set as a user-specified option
    relaxed_optimizer::MOI.AbstractOptimizer = GLPK.Optimizer()

    # set as a user-specified option (if empty set to all nonlinear by TODO in TODO)
    obbt_variable_values::Vector{Bool} = Bool[]

    # Upper bounding options (set as a user-specified option)
    upper_optimizer::MOI.AbstractOptimizer = default_nlp_solver()

    # Extensions (set as user-specified option)
    enable_optimize_hook::Bool = false
    ext::Dict{Symbol, Any} = Dict{Symbol,Any}()
    ext_type::ExtensionType = DefaultExt()

    # set as user-specified option
    _parameters::EAGOParameters = EAGOParameters()

    # set by MOI manipulations (see Input problem structure)
    _input_problem::InputProblem = InputProblem()

    # loaded from _input_problem by TODO
    _working_problem::ParsedProblem = ParsedProblem()

    _end_state::GlobalEndState = GS_UNSET
    _termination_status_code::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _result_status_code::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _obj_mult::Float64 = 1.0

    _stack::BinaryMinMaxHeap{NodeBB} = BinaryMinMaxHeap{NodeBB}()

    # set in node_selection!
    _current_node::NodeBB = NodeBB()

    _first_relax_point_set::Bool = false
    _current_xref::Vector{Float64} = Float64[]
    _candidate_xref::Vector{Float64} = Float64[]

    _use_prior_objective_xref::Bool = false
    _current_objective_xref::Vector{Float64} = Float64[]
    _prior_objective_xref::Vector{Float64} = Float64[]

    # set in label_branch_variables! and label_fixed_variables! respectively in parse.jl
    _user_branch_variables::Bool = false
    _fixed_variable::Vector{Bool} = Bool[]
    _branch_variable_count::Int = 0
    _branch_to_sol_map::Vector{Int} = Int[]
    _sol_to_branch_map::Vector{Int} = Int[]

    _continuous_solution::Vector{Float64} = Float64[]

    # all subproblem immutable subproblem status are set in global_solve in corresponding routines
    # in optimize_nonconvex.jl
    _preprocess_feasibility::Bool = true
    _preprocess_primal_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _preprocess_dual_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _preprocess_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED

    _lower_primal_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _lower_dual_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _lower_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _lower_feasibility::Bool = true
    _lower_objective_value::Float64 = -Inf

    # set in TODO
    _lower_solution::Vector{Float64} = Float64[]
    _lower_lvd::Vector{Float64} = Float64[]
    _lower_uvd::Vector{Float64} = Float64[]

    _last_cut_objective::Float64 = -Inf

    _upper_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _upper_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _upper_feasibility::Bool = true
    _upper_objective_value::Float64 = Inf

    # array is initialized to correct size in TODO, reset in single_nlp_solve! in optimize_convex.jl
    _upper_variables::Vector{VI} =  VI[]

    # set in TODO
    _upper_solution::Vector{Float64} = Float64[]

    _postprocess_feasibility::Bool = true

    # set to time limit in initial_parse! in parse.jl, decremented throughout global_solve in optimize_nonconvex.jl
    _time_left::Float64 = 1000.0

    # set constructor reset on empty! and  to zero in initial parse! in parse.jl
    _start_time::Float64 = 0.0
    _run_time::Float64 = 0.0
    _parse_time::Float64 = 0.0
    _presolve_time::Float64 = 0.0
    _last_preprocess_time::Float64 = 0.0
    _last_lower_problem_time::Float64 = 0.0
    _last_upper_problem_time::Float64 = 0.0
    _last_postprocessing_time::Float64 = 0.0

    # reset in initial_parse! in parse.jl
    _min_converged_value::Float64 = Inf
    _global_lower_bound::Float64 = -Inf
    _global_upper_bound::Float64 = Inf
    _maximum_node_id::Int = 0
    _iteration_count::Int = 0
    _node_count::Int = 0

    # Storage for output, reset in initial_parse! in parse.jl
    _solution_value::Float64 = 0.0
    _feasible_solution_found::Bool = false
    _first_solution_node::Int = -1
    _objective_value::Float64 = -Inf
    _best_upper_value::Float64 = Inf

    # Optimality-Based Bound Tightening (OBBT) Options
    # set in TODO
    _obbt_working_lower_index::Vector{Bool} = Bool[]
    _obbt_working_upper_index::Vector{Bool} = Bool[]
    _lower_indx_diff::Vector{Bool} = Bool[]
    _upper_indx_diff::Vector{Bool} = Bool[]
    _old_low_index::Vector{Bool} = Bool[]
    _old_upp_index::Vector{Bool} = Bool[]
    _new_low_index::Vector{Bool} = Bool[]
    _new_upp_index::Vector{Bool} = Bool[]
    _obbt_variables::Vector{VI} = VI[]
    _obbt_variable_count::Int = 0
    _obbt_performed_flag::Bool = false

    # Buffers for fbbt, set in presolve, used in preprocess
    _lower_fbbt_buffer::Vector{Float64} = Float64[]
    _upper_fbbt_buffer::Vector{Float64} = Float64[]

    # Feasibility-Based Bound Tightening Options
    # set in set_constraint_propagation_fbbt in domain_reduction.jl
    _cp_improvement::Float64 = 0.0
    _cp_evaluation_reverse::Bool = false

    _cut_iterations::Int = 0
    _cut_add_flag::Bool = false

    # Options for Repetition (If DBBT Performed Well)
    # set in within preprocess in optimize_nonconvex.jl
    _node_repetitions::Int = 0
    _initial_volume::Float64 = 0.0
    _final_volume::Float64 = 0.0

    # Log
    _log::Log = Log()

    _affine_relax_ci = CI{SAF,LT}[]
    _affine_objective_cut_ci::Union{CI{SV,LT},CI{SAF,LT},Nothing} = nothing

    # need to retreive primal _relaxed_variable_index
    # set in TODO
    #"Number of variables actively branched on in B&B routine (excludes linear and fixed)"
    _relaxed_variable_number::Int = 0
    _relaxed_variable_index::Vector{VI} = VI[]
    _relaxed_variable_eq::Vector{Tuple{CI{SV, ET}, Int}} = Tuple{CI{SV, ET}, Int}[]
    _relaxed_variable_lt::Vector{Tuple{CI{SV, LT}, Int}} = Tuple{CI{SV, LT}, Int}[]
    _relaxed_variable_gt::Vector{Tuple{CI{SV, GT}, Int}} = Tuple{CI{SV, GT}, Int}[]

    # set as user-input
    _branch_variables::Vector{Bool} = Bool[]

    _new_eval_constraint::Bool = false
    _new_eval_objective::Bool = false

    _node_to_sv_leq_ci::Vector{CI{SV,LT}} = CI{SV,LT}[]
    _node_to_sv_geq_ci::Vector{CI{SV,GT}} = CI{SV,GT}[]

    #"Set to true if a nonlinear evaluator was created (NLconstraint or NLobjective specified)"
    _nonlinear_evaluator_created::Bool = false

    _branch_cost::BranchCostStorage{Float64} = BranchCostStorage{Float64}()
    _branch_variable_sparsity::SparseMatrixCSC{Bool,Int} = spzeros(Bool,1,1)
    _constraint_infeasiblity::Vector{Float64} = Float64[]

    #_relaxed_evaluator::Evaluator = Evaluator{1,NS}()
    #_relaxed_constraint_bounds::Vector{MOI.NLPBoundsPair} = Vector{MOI.NLPBoundsPair}[]
end

#####
#####
##### General MOI utilities required
#####
#####

const EAGO_OPTIMIZER_ATTRIBUTES = Symbol[:relaxed_optimizer, :relaxed_optimizer_kwargs, :upper_optimizer,
                                         :enable_optimize_hook, :ext, :ext_type, :_parameters]
const EAGO_MODEL_STRUCT_ATTRIBUTES = Symbol[:_stack, :_log, :_current_node, :_working_problem, :_input_problem, :_branch_cost]
const EAGO_MODEL_NOT_STRUCT_ATTRIBUTES = setdiff(fieldnames(Optimizer), union(EAGO_OPTIMIZER_ATTRIBUTES,
                                                                              EAGO_MODEL_STRUCT_ATTRIBUTES))

function MOI.empty!(m::Optimizer)

    # create a new empty optimizer and copy fields to m
    new_optimizer = Optimizer()
    for field in union(EAGO_MODEL_STRUCT_ATTRIBUTES, EAGO_MODEL_NOT_STRUCT_ATTRIBUTES)
        setfield!(m, field, getfield(new_optimizer, field))
    end

    return nothing
end

function MOI.is_empty(m::Optimizer)

    is_empty_flag = uninitialized(_current_node(m))
    is_empty_flag &= isempty(m._stack)
    is_empty_flag &= isempty(m._log)
    is_empty_flag &= isempty(m._input_problem)
    is_empty_flag &= isempty(m._working_problem)

    new_optimizer = Optimizer()
    for field in EAGO_MODEL_NOT_STRUCT_ATTRIBUTES
        if getfield(m, field) != getfield(new_optimizer, field)
            is_empty_flag = false
            break
        end
    end

    return is_empty_flag
end

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; copy_names = false)
    return MOI.Utilities.default_copy_to(model, src, copy_names)
end

#####
#####
##### Set & get attributes of model
#####
#####

function MOI.set(m::Optimizer, ::MOI.Silent, value)

     m._parameters.verbosity = 0
     m._parameters.log_on = false
     return nothing

end

function MOI.set(m::Optimizer, ::MOI.TimeLimitSec, value::Nothing)
    m._parameters.time_limit = Inf
    return
end

function MOI.set(m::Optimizer, ::MOI.TimeLimitSec, value::Float64)
    m._parameters.time_limit = value
    return
end

function MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices)
    return [MOI.VariableIndex(i) for i = 1:length(m._input_problem._variable_info)]
end

function MOI.get(m::Optimizer, ::MOI.ObjectiveValue)
    mult = 1.0
    if m._input_problem._optimization_sense === MOI.MAX_SENSE
        mult *= -1.0
    end
    return mult*m._objective_value
end

MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = m._input_problem._variable_count

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
MOI.get(m::Optimizer, ::MOI.TimeLimitSec) = m.time_limit

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    check_inbounds!(model, vi)
    return model._continuous_solution[vi.value]
end

const EAGO_PARAMETERS = fieldnames(EAGOParameters)

_to_sym(d) = error("EAGO only supports raw parameters with Symbol or String names.")
_to_sym(d::String) = Symbol(d)
_to_sym(d::Symbol) = d

function MOI.get(m::Optimizer, p::MOI.RawParameter)
    s = _to_sym(p.name)
    s in EAGO_PARAMETERS ? getfield(m._parameters, s) : getfield(m, s)
end
function MOI.set(m::Optimizer, p::MOI.RawParameter, x)
    s = _to_sym(p.name)
    s in EAGO_PARAMETERS ? setfield!(m._parameters, s, x) : setfield!(m, s, x)
    return
end

#####
#####
##### Support, set, and evaluate objective functions
#####
#####
MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{F}) where {F <: Union{SV, SAF, SQF}} = true

function MOI.set(m::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    if nlp_data.has_objective
        m._input_problem._objective = nothing
    end
    m._input_problem._nlp_data = nlp_data
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction{T}, f::T) where T <: Union{SV,SAF,SQF}
    check_inbounds!(m, f)
    m._input_problem._objective = f
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveSense, s::MOI.OptimizationSense)
    m._input_problem._optimization_sense = s
    return
end

@inline _branch_cost(m::Optimizer) = m._branch_cost.cost
@inline _cost_offset_β(m::Optimizer) = m._branch_cost.β
@inline _branch_cvx_α(m::Optimizer) =  m._parameters.branch_cvx_factor
@inline _branch_offset_β(m::Optimizer) = m._parameters.branch_offset
@inline _branch_pseudocost_on(m::Optimizer) = m._parameters.branch_pseudocost_on
@inline _cut_ϵ_abs(m::Optimizer) = m._parameters.cut_tolerance_abs
@inline _cut_ϵ_rel(m::Optimizer) = m._parameters.cut_tolerance_rel
@inline _cut_max_iterations(m::Optimizer) = m._parameters.cut_max_iterations
@inline _constraint_tol(m::Optimizer) = m._parameters.absolute_constraint_feas_tolerance

@inline _bvi(m::Optimizer, i::Int) = m._branch_to_sol_map[i]
@inline _svi(m::Optimizer, i::Int) = m._sol_to_branch_map[i]

@inline _is_branch_var(m::Optimizer, i) = m._branch_variables[i]

@inline _current_node(m::Optimizer) = m._current_node
@inline _variable_info(m::Optimizer, i::Int) = m._input_problem._variable_info[i]
@inline _working_variable_info(m::Optimizer, i::Int) = m._working_problem._variable_info[i]

@inline _sparsity(::BranchVar, m::Optimizer, i::Int) = view(m._branch_variable_sparsity,1,:)

@inline _variable_num(::BranchVar, m::Optimizer) = m._branch_variable_count
@inline _variable_num(::FullVar, m::Optimizer) = m._working_problem._variable_count

@inline _is_integer(::BranchVar, m::Optimizer, i::Int) = is_integer(_current_node(m), i)
@inline function _is_integer(::FullVar, m::Optimizer, i::Int)
     if _is_branch_var(m,i)
        return is_integer(_current_node(m), _svi(m, i))
     end
    is_integer(_working_variable_info(m,i))
end

@inline _lower_bound(::BranchVar, m::Optimizer, i::Int) = lower_variable_bounds(_current_node(m), i)
@inline function _lower_bound(::FullVar, m::Optimizer, i::Int)
     if _is_branch_var(m,i)
         return lower_variable_bounds(_current_node(m), _svi(m, i))
     end
    lower_bound(_working_variable_info(m,i))
end

@inline _upper_bound(::BranchVar, m::Optimizer, i::Int) = upper_variable_bounds(_current_node(m), i)
@inline function _upper_bound(::FullVar, m::Optimizer, i::Int)
    if _is_branch_var(m,i)
        return upper_variable_bounds(_current_node(m), _svi(m, i))
    end
    upper_bound(_working_variable_info(m,i))
end

@inline _mid(::BranchVar, m::Optimizer, i::Int) = mid(_current_node(m), i)
@inline function _mid(::FullVar, m::Optimizer, i::Int)
    if _is_branch_var(m,i)
        return mid(_current_node(m), _svi(m, i))
    end
    mid(_working_variable_info(m,i))
end

@inline _diam(::BranchVar, m::Optimizer, i::Int) = diam(_current_node(m), i)
@inline function _diam(::FullVar, m::Optimizer, i::Int)
    if _is_branch_var(m,i)
        return diam(_current_node(m), _svi(m, i))
    end
    diam(_working_variable_info(m,i))
end

@inline _lower_solution(::BranchVar, m::Optimizer, i::Int) = m._lower_solution[_bvi(m, i)]
@inline _lower_solution(::FullVar, m::Optimizer, i::Int) = m._lower_solution[i]

@inline function _set_lower_solution!(::BranchVar, m::Optimizer, v::Float64, i::Int)
    m._lower_solution[_bvi(m, i)] = v
    return
end
