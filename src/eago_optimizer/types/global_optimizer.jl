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

#=
LP          -> COPY TO RELAXED SOLVER AND SOLVE
MILP        -> COPY TO RELAXED SOLVER AND SOLVE
SOCP        -> COPY TO RELAXED SOLVER AND SOLVE
MISOCP      -> COPY TO RELAXED SOLVER AND SOLVE
DIFF_CVX    -> COPY TO NLP SOLVER AND SOLVE (POTENTIAL MULTISTART)
NS_CVX      -> COPY TO NLP SOLVER AND SOLVE (POTENTIAL MULTISTART)
DIFF_NCVX   -> APPLY GLOBAL SOLVER (UNLESS USER REQUEST LOCAL SOLVE THEN NLP)
NS_NCVX     -> APPLY GLOBAL SOLVER (UNLESS USER REQUEST LOCAL SOLVE THEN NLP)
MINCVX      -> APPLY GLOBAL SOLVER (LOCAL SOLVE OPTION FUTURE FEATURE)
=#

abstract type AbstractProblemType end
struct LP <: AbstractProblemType end
struct MILP <: AbstractProblemType end
struct SOCP <: AbstractProblemType end
struct MISOCP <: AbstractProblemType end
struct DIFF_CVX <: AbstractProblemType end
struct MINCVX <: AbstractProblemType end

const ANY_PROBLEM_TYPE = Union{Nothing, LP, MILP, SOCP, MISOCP, DIFF_CVX, MINCVX}

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
    "Display summary of iteration to console every `output_iterations` (default = 1000)"
    output_iterations::Int = 1
    "Display header for summary to console every `output_iterations` (default = 10000)"
    header_iterations::Int = 100000

    # Node branching options
    "Convex coefficient used to select branch point. Branch point is given by
    `branch_cvx_factor*xmid + (1-branch_cvx_factor)*xsol` (default = 0.5)"
    branch_cvx_factor::Float64 = 0.5
    "Minimum distance from bound to have branch point normalized by width of
    dimension to branch on (default = 0.2)"
    branch_offset::Float64 = 0.2
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
    time_limit::Float64 = 100.0
    "Maximum number of iterations (default 3E6)"
    iteration_limit::Int = 1E9 #2*10^5
    "Absolute tolerance for termination (default = 1E-3)"
    absolute_tolerance::Float64 = 1E-4
    "Relative tolerance for termination (default = 1E-3)"
    relative_tolerance::Float64 = 1E-4
    "Absolute constraint feasibility tolerance"
    absolute_constraint_feas_tolerance::Float64 = 1E-7

    # Options for constraint propagation
    "Depth in B&B tree above which constraint propagation should be disabled (default = 1000)"
    cp_depth::Int = 0
    "Number of times to repeat forward-reverse pass routine (default = 3)"
    cp_repetitions::Int = 0
    "Disable constraint propagation if the ratio of new node volume to beginning node volume exceeds
    this number (default = 0.99)"
    cp_tolerance::Float64 = 0.99
    "Use only valid interval bounds during constraint propagation (default = false)"
    cp_interval_only::Bool = false

    # obbt options
    "Depth in B&B tree above which OBBT should be disabled (default = 6)"
    obbt_depth::Int = 4
    "Number of repetitions of OBBT to perform in preprocessing (default = 3)"
    obbt_repetitions::Int = 10
    "Turn aggresive OBBT on (default = false)"
    obbt_aggressive_on::Bool = true
    "Maximum iteration to perform aggresive OBBT (default = 2)"
    obbt_aggressive_max_iteration::Int = 2
    "Minimum dimension to perform aggresive OBBT (default = 2)"
    obbt_aggressive_min_dimension::Int = 2
    "Tolerance to consider bounds equal (default = 1E-10)"
    obbt_tolerance::Float64 = 1E-10

    # Options for linear bound tightening
    "Depth in B&B tree above which linear FBBT should be disabled (default = 1000)"
    fbbt_lp_depth::Int  = 1000
    "Number of repetitions of linear FBBT to perform in preprocessing (default = 3)"
    fbbt_lp_repetitions::Int  = 3

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
    mul_relax_style::Int = 1

    # Tolerance to add cuts and max number of cuts
    "Minimum number of cuts at each node to attempt (unsafe cuts not necessarily added)"
    cut_min_iterations::Int = 2
    "Maximum number of cuts at each node to attempt"
    cut_max_iterations::Int = 8
    "Absolute tolerance checked for continuing cut"
    cut_tolerance_abs::Float64 = 1E-6
    "Relative tolerance checked for continuing cut"
    cut_tolerance_rel::Float64 = 1E-2

    "Use tolerances to determine safe cuts in a Khajavirad 2018 manner"
    cut_safe_on::Bool = true
    "Lower tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_l::Float64 = 1E-8
    "Upper tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_u::Float64 = 1E8
    "Constant tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_b::Float64 = 1E9

    "Solve upper problem for every node with depth less than `upper_bounding_depth`
    and with a probabilityof (1/2)^(depth-upper_bounding_depth) otherwise (default = 8)"
    upper_bounding_depth::Int = 8

    # handling for domain violations
    "Amount about a domain violation to ignore when propagating bounds."
    domain_violation_guard_on::Bool = false
    "Amount about a domain violation to ignore when propagating bounds."
    domain_violation_ϵ::Float64 = 1E-9

    "If true, then EAGO forgos its default configuration process for subsolvers"
    user_solver_config::Bool = false

    integer_abs_tol::Float64 = 1E-9
    integer_rel_tol::Float64 = 1E-9
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
    _linear_leq_constraints::Vector{Tuple{SAF, LT, Int}} = Tuple{SAF, LT, Int}[]
    _linear_geq_constraints::Vector{Tuple{SAF, GT, Int}} = Tuple{SAF, GT, Int}[]
    _linear_eq_constraints::Vector{Tuple{SAF, ET, Int}} = Tuple{SAF, ET, Int}[]

    # quadratic constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _quadratic_leq_constraints::Vector{Tuple{SQF, LT, Int}} = Tuple{SQF, LT, Int}[]
    _quadratic_geq_constraints::Vector{Tuple{SQF, GT, Int}} = Tuple{SQF, GT, Int}[]
    _quadratic_eq_constraints::Vector{Tuple{SQF, ET, Int}} = Tuple{SQF, ET, Int}[]

    # dictionary mapping ci for input problem to local nlp problem
    _linear_leq_ci_dict::Dict{Int,CI{SAF,LT}} = Dict{Int,CI{SAF,LT}}()
    _linear_geq_ci_dict::Dict{Int,CI{SAF,GT}} = Dict{Int,CI{SAF,GT}}()
    _linear_eq_ci_dict::Dict{Int,CI{SAF,ET}} = Dict{Int,CI{SAF,ET}}()
    _quadratic_leq_ci_dict::Dict{Int,CI{SQF,LT}} = Dict{Int,CI{SQF,LT}}()
    _quadratic_geq_ci_dict::Dict{Int,CI{SQF,GT}} = Dict{Int,CI{SQF,GT}}()
    _quadratic_eq_ci_dict::Dict{Int,CI{SQF,ET}} = Dict{Int,CI{SQF,ET}}()

    _constraint_primal::Dict{Int,Float64} = Dict{Int,Float64}()

    # conic constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _conic_second_order::Vector{Tuple{VECOFVAR, MOI.SecondOrderCone}} = Tuple{VECOFVAR, MOI.SecondOrderCone}[]

    # nonlinear constraint storage
    _objective::Union{SV,SAF,SQF,Nothing} = nothing

    # nlp constraints (set by MOI.set(m, ::NLPBlockData...) in optimizer.jl)
    _nlp_data::Union{MOI.NLPBlockData,Nothing} = nothing

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

    is_empty_flag &= x._nlp_data === nothing
    is_empty_flag &= x._objective === nothing

    return is_empty_flag
end

"""
$(TYPEDEF)

A structure used to expressions and problem descriptions EAGO uses to formulate
relaxed problems.
"""
Base.@kwdef mutable struct ParsedProblem

    # Problem classification (set in parse_classify_problem!)
    _problem_type::ANY_PROBLEM_TYPE = nothing

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

    # nlp constraints
    _nlp_data::Union{MOI.NLPBlockData,Nothing} = nothing
    _nonlinear_constr::Vector{BufferedNonlinearFunction} = BufferedNonlinearFunction[]
    _relaxed_evaluator::Evaluator = Evaluator()

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

    is_empty_flag &= x._nlp_data === nothing
    is_empty_flag &= isempty(x._objective_saf.terms)
    is_empty_flag &= x._objective_saf.constant === 0.0
    is_empty_flag &= x._objective === nothing

    return is_empty_flag
end

Base.@kwdef mutable struct GlobalOptimizer{R,Q,S<:ExtensionType} <: MOI.AbstractOptimizer
    
    _subsolvers::SubSolvers{R,Q,S}
    _parameters::EAGOParameters = EAGOParameters()
    _input_problem::InputProblem = InputProblem()
    _working_problem::ParsedProblem = ParsedProblem()

    obbt_variable_values::Vector{Bool} = Bool[]

    enable_optimize_hook::Bool = false
    ext::Dict{Symbol, Any} = Dict{Symbol,Any}()

    # set as user-specified option
    _end_state::GlobalEndState = GS_UNSET
    _termination_status_code::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _result_status_code::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _obj_mult::Float64 = 1.0
    _obj_var_slack_added::Bool = false

    _stack::BinaryMinMaxHeap{NodeBB} = BinaryMinMaxHeap{NodeBB}()
    _current_node::NodeBB = NodeBB()

    _first_relax_point_set::Bool = false
    _current_xref::Vector{Float64} = Float64[]
    _candidate_xref::Vector{Float64} = Float64[]

    _use_prior_objective_xref::Bool = false
    _current_objective_xref::Vector{Float64} = Float64[]
    _prior_objective_xref::Vector{Float64} = Float64[]

    _user_branch_variables::Bool = false
    _fixed_variable::Vector{Bool} = Bool[]
    _branch_variable_count::Int = 0
    _branch_to_sol_map::Vector{Int} = Int[]
    _sol_to_branch_map::Vector{Int} = Int[]

    _continuous_solution::Vector{Float64} = Float64[]
    _constraint_primal::Dict{Int,Float64} = Dict{Int,Float64}()

    _preprocess_feasibility::Bool = true
    _preprocess_primal_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _preprocess_dual_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _preprocess_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED

    _lower_primal_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _lower_dual_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _lower_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _lower_feasibility::Bool = true
    _lower_objective_value::Float64 = -Inf
    _lower_solution::Vector{Float64} = Float64[]
    _lower_lvd::Vector{Float64} = Float64[]
    _lower_uvd::Vector{Float64} = Float64[]

    _last_cut_objective::Float64 = -Inf

    _upper_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    _upper_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _upper_feasibility::Bool = true
    _upper_objective_value::Float64 = Inf
    _upper_variables::Vector{VI} =  VI[]
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
    _obbt_variable_count::Int = 0
    _obbt_performed_flag::Bool = false

    # Buffers for fbbt, set in presolve, used in preprocess
    _lower_fbbt_buffer::Vector{Float64} = Float64[]
    _upper_fbbt_buffer::Vector{Float64} = Float64[]
    _cp_improvement::Float64 = 0.0
    _cp_evaluation_reverse::Bool = false

    _cut_iterations::Int = 0
    _cut_add_flag::Bool = false
    _node_repetitions::Int = 0

    _log::Log = Log()

    _affine_relax_ci::Vector{CI{SAF,LT}} = CI{SAF,LT}[]
    _affine_objective_cut_ci::Union{CI{SV,LT},CI{SAF,LT},Nothing} = nothing

    _relaxed_variable_number::Int = 0
    _relaxed_variable_index::Vector{VI} = VI[]
    _relaxed_variable_et::Vector{CI{SV, ET}} = CI{SV, ET}[]
    _relaxed_variable_lt::Vector{Tuple{CI{SV, LT}, Int}} = Tuple{CI{SV, LT}, Int}[]
    _relaxed_variable_gt::Vector{Tuple{CI{SV, GT}, Int}} = Tuple{CI{SV, GT}, Int}[]
    _relaxed_variable_integer::Vector{CI{SV, MOI.Integer}} = CI{SV, MOI.Integer}[]

    _branch_variables::Vector{Bool} = Bool[]
    _nonbranching_int::Bool = false

    _new_eval_constraint::Bool = false
    _new_eval_objective::Bool = false

    _node_to_sv_leq_ci::Dict{Int,CI{SV,LT}} = Dict{Int,CI{SV,LT}}()
    _node_to_sv_geq_ci::Dict{Int,CI{SV,GT}} = Dict{Int,CI{SV,GT}}()
    _nonlinear_evaluator_created::Bool = false

    _branch_cost::BranchCostStorage{Float64} = BranchCostStorage{Float64}()
    _branch_variable_sparsity::SparseMatrixCSC{Bool,Int} = spzeros(Bool,1,1)
    _constraint_infeasiblity::Vector{Float64} = Float64[]
end

const EAGO_OPTIMIZER_ATTRIBUTES = Symbol[:relaxed_optimizer, :upper_optimizer,
                                         :enable_optimize_hook, :ext, :ext_type, :_parameters]
const EAGO_MODEL_STRUCT_ATTRIBUTES = Symbol[:_stack, :_log, :_current_node, :_working_problem, :_input_problem, :_branch_cost]
const EAGO_MODEL_NOT_STRUCT_ATTRIBUTES = setdiff(fieldnames(GlobalOptimizer), union(EAGO_OPTIMIZER_ATTRIBUTES, EAGO_MODEL_STRUCT_ATTRIBUTES))
                           
function MOI.empty!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q}
                           
    # create a new empty optimizer and copy fields to m
    new_optimizer = GlobalOptimizer{R,S,Q}(_subsolvers = m._subsolvers,
                                           _parameters = m._parameters, 
                                           _input_problem = m._input_problem,
                                           _working_problem = m._working_problem)

    MOI.empty!(new_optimizer._subsolvers)
   # MOI.empty!(new_optimizer._input_problem)
   # MOI.empty!(new_optimizer._working_problem)
    for field in union(EAGO_MODEL_STRUCT_ATTRIBUTES, EAGO_MODEL_NOT_STRUCT_ATTRIBUTES)
        setfield!(m, field, getfield(new_optimizer, field))
    end
                           
    return nothing
end
                           
function MOI.is_empty(m::GlobalOptimizer{R,S,Q}) where {R,S,Q}
                           
    is_empty_flag = uninitialized(_current_node(m))
    is_empty_flag &= isempty(m._stack)
    is_empty_flag &= isempty(m._log)
    is_empty_flag &= isempty(m._input_problem)
    is_empty_flag &= isempty(m._working_problem)
                           
    new_optimizer = GlobalOptimizer{R,S,Q}(_subsolvers = m._subsolvers,
                                           _parameters = m._parameters, 
                                           _input_problem = m._input_problem,
                                           _working_problem = m._working_problem)
    for field in EAGO_MODEL_NOT_STRUCT_ATTRIBUTES
        if getfield(m, field) != getfield(new_optimizer, field)
            is_empty_flag = false
            break
        end
    end
                           
    return is_empty_flag
end

MOI.get(m::GlobalOptimizer, ::MOI.ObjectiveBound) = _is_input_min(m) ? m._global_lower_bound : -m._global_upper_bound
MOI.get(m::GlobalOptimizer, ::MOI.ObjectiveValue) = _is_input_min(m) ? m._global_upper_bound : -m._global_lower_bound

_relaxed_optimizer(m::GlobalOptimizer{R,S,Q}) where {R,S,Q} = m._subsolvers.relaxed_optimizer
_upper_optimizer(m::GlobalOptimizer{R,S,Q})   where {R,S,Q} = m._subsolvers.upper_optimizer
_ext_typ(m::GlobalOptimizer{R,S,Q})           where {R,S,Q} = m._subsolvers.ext_typ

@inline _is_input_min(m::GlobalOptimizer) = m._input_problem._optimization_sense == MOI.MIN_SENSE
@inline _branch_cost(m::GlobalOptimizer) = m._branch_cost.cost
@inline _cost_offset_β(m::GlobalOptimizer) = m._branch_cost.β
@inline _branch_cvx_α(m::GlobalOptimizer) =  m._parameters.branch_cvx_factor
@inline _branch_offset_β(m::GlobalOptimizer) = m._parameters.branch_offset
@inline _branch_pseudocost_on(m::GlobalOptimizer) = m._parameters.branch_pseudocost_on
@inline _cut_ϵ_abs(m::GlobalOptimizer) = m._parameters.cut_tolerance_abs
@inline _cut_ϵ_rel(m::GlobalOptimizer) = m._parameters.cut_tolerance_rel
@inline _cut_max_iterations(m::GlobalOptimizer) = m._parameters.cut_max_iterations

@inline _integer_abs_tol(m::GlobalOptimizer) = m._parameters.integer_abs_tol
@inline _integer_rel_tol(m::GlobalOptimizer) = m._parameters.integer_rel_tol

@inline _absolute_tol(m::GlobalOptimizer) = m._parameters.absolute_tolerance
@inline _relative_tol(m::GlobalOptimizer) = m._parameters.relative_tolerance
@inline _constraint_tol(m::GlobalOptimizer) = m._parameters.absolute_constraint_feas_tolerance

@inline _fbbt_lp_depth(m::GlobalOptimizer) = m._parameters.fbbt_lp_depth
@inline _fbbt_lp_repetitions(m::GlobalOptimizer) = m._parameters.fbbt_lp_repetitions

@inline _obbt_depth(m::GlobalOptimizer) = m._parameters.obbt_depth
@inline _obbt_repetitions(m::GlobalOptimizer) = m._parameters.obbt_repetitions
@inline _obbt_tolerance(m::GlobalOptimizer) = m._parameters.obbt_repetitions
@inline _obbt_aggressive_on(m::GlobalOptimizer) = m._parameters.obbt_aggressive_on
@inline _obbt_aggressive_max_iteration(m::GlobalOptimizer) = m._parameters.obbt_aggressive_max_iteration

@inline _user_solver_config(m::GlobalOptimizer) = m._parameters.user_solver_config
@inline _verbosity(m::GlobalOptimizer) = m._parameters.verbosity

@inline _cp_depth(m::GlobalOptimizer) = m._parameters.cp_depth
@inline _cp_repetitions(m::GlobalOptimizer) = m._parameters.cp_repetitions

@inline _iteration_count(m::GlobalOptimizer) = m._iteration_count
@inline _obj_var_slack_added(m::GlobalOptimizer) = m._obj_var_slack_added

@inline _bvi(m::GlobalOptimizer, i::Int) = m._branch_to_sol_map[i]
@inline _svi(m::GlobalOptimizer, i::Int) = m._sol_to_branch_map[i]

@inline _is_branch_var(m::GlobalOptimizer, i) = m._branch_variables[i]

@inline _current_node(m::GlobalOptimizer) = m._current_node
@inline _variable_info(m::GlobalOptimizer, i::Int) = m._input_problem._variable_info[i]
@inline _working_variable_info(m::GlobalOptimizer, i::Int) = m._working_problem._variable_info[i]

@inline _sparsity(::BranchVar, m::GlobalOptimizer, i::Int) = view(m._branch_variable_sparsity,1,:)

@inline _variable_num(::BranchVar, m::GlobalOptimizer) = m._branch_variable_count
@inline _variable_num(::FullVar, m::GlobalOptimizer) = m._working_problem._variable_count

@inline is_integer(::BranchVar, m::GlobalOptimizer, i::Int) = is_integer(_current_node(m), i)
@inline function is_integer(::FullVar, m::GlobalOptimizer, i::Int)
    if _is_branch_var(m,i)
        return is_integer(_current_node(m), _svi(m, i))
    end
    is_integer(_working_variable_info(m,i))
end

@inline _lower_bound(::BranchVar, m::GlobalOptimizer, i::Int) = lower_variable_bounds(_current_node(m), i)
@inline function _lower_bound(::FullVar, m::GlobalOptimizer, i::Int)
     if _is_branch_var(m,i)
         return lower_variable_bounds(_current_node(m), _svi(m, i))
     end
    lower_bound(_working_variable_info(m,i))
end

@inline _upper_bound(::BranchVar, m::GlobalOptimizer, i::Int) = upper_variable_bounds(_current_node(m), i)
@inline function _upper_bound(::FullVar, m::GlobalOptimizer, i::Int)
    if _is_branch_var(m,i)
        return upper_variable_bounds(_current_node(m), _svi(m, i))
    end
    upper_bound(_working_variable_info(m,i))
end

@inline _mid(::BranchVar, m::GlobalOptimizer, i::Int) = mid(_current_node(m), i)
@inline function _mid(::FullVar, m::GlobalOptimizer, i::Int)
    if _is_branch_var(m,i)
        return mid(_current_node(m), _svi(m, i))
    end
    mid(_working_variable_info(m,i))
end

@inline _diam(::BranchVar, m::GlobalOptimizer, i::Int) = diam(_current_node(m), i)
@inline function _diam(::FullVar, m::GlobalOptimizer, i::Int)
    if _is_branch_var(m,i)
        return diam(_current_node(m), _svi(m, i))
    end
    diam(_working_variable_info(m,i))
end

@inline _lower_solution(::BranchVar, m::GlobalOptimizer, i::Int) = m._lower_solution[_bvi(m, i)]
@inline _lower_solution(::FullVar, m::GlobalOptimizer, i::Int) = m._lower_solution[i]

@inline function _set_lower_solution!(::BranchVar, m::GlobalOptimizer, v::Float64, i::Int)
    m._lower_solution[_bvi(m, i)] = v
    return
end
@inline function _set_lower_solution!(::FullVar, m::GlobalOptimizer, v::Float64, i::Int)
    m._lower_solution[i] = v
    return
end