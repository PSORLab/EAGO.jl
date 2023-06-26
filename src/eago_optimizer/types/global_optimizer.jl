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

"""
$(TYPEDEF)

An Enum of possible values for EAGO's termination status. This attribute is used
by EAGO to explain why the optimizer stopped executing in the most recent call
to `optimize!`. See also [`MathOptInterface.TerminationStatusCode`](@ref).

If no call has been made to `optimize!`, the `GlobalEndState` value is:
- `GS_UNSET`: The optimization algorithm has not stated.

# OK
- `GS_OPTIMAL`: A globally optimal solution was found.
- `GS_INFEASIBLE`: The algorithm concluded that no feasible solution exists.

# Limits reached
- `GS_NODE_LIMIT`: The branch-and-bound algorithm stopped because it reached the
        user-set maximum number of nodes in the branch-and-bound tree.
- `GS_ITERATION_LIMIT`: The maximum number of iterations was reached.
- `GS_RELATIVE_TOL`: The gap between the lower and upper bounds, relative 
        to the bound with the larger magnitude, is within the user-set relative
        tolerance.
- `GS_ABSOLUTE_TOL`: The gap between the lower and upper bounds is within the user-set 
        absolute tolerance.
- `GS_TIME_LIMIT`: The algorithm stopped after the user-specified time limit was reached.

"""
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
    "Should EAGO attempt to remove type-assert issues for user-defined functions (default = false)"
    presolve_scrubber_flag::Bool = false
    "Create and use DAG representations of user-defined functions (default = false)"
    presolve_to_JuMP_flag::Bool = false
    "Rerrange the DAG using registered transformations (default = false)"
    presolve_flatten_flag::Bool = false

    # Conic reformulations
    "Attempt to bridge convex constraint to second-order cone (default = false)"
    conic_convert_quadratic::Bool = false

    # Iteration logging options
    "Turn logging on; record global bounds, node count, and run time. Additional
     options are available for recording information specific to subproblems (default = false)"
    log_on::Bool = false
    "Turn on logging of times and feasibility of subproblems (default = false)"
    log_subproblem_info::Bool = false
    "Log data every `log_interval` iterations (default = 1)"
    log_interval::Int = 1

    # Optimizer display options
    "The amount of information that should be printed to console while solving.
    Values range from 0 - 4: 0 is silent, 1 shows iteration summary statistics
    only, 2-4 show varying degrees of detail about calculations within each
    iteration (default = 1)"
    verbosity::Int = 1
    "Display summary of iteration to console every `output_iterations` (default = 1000)"
    output_iterations::Int = 1000
    "Display header for summary to console every `output_iterations` (default = 100000)"
    header_iterations::Int = 100000

    # Node branching options
    "Convex coefficient used to select branch point. Branch point is given by
    `branch_cvx_factor*xmid + (1-branch_cvx_factor)*xsol` (default = 0.25)"
    branch_cvx_factor::Float64 = 0.25
    "Minimum distance from bound to have branch point, normalized by width of
    dimension to branch on (default = 0.15)"
    branch_offset::Float64 = 0.15
    "Indicate that pseudocost branching should be used (default = false)"
    branch_pseudocost_on::Bool = false
    "Variables to branch on (default is all nonlinear)"
    branch_variable::Vector{Bool} = Bool[]
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Number of times to repeat node
    processing prior to branching (default = 4)"
    branch_max_repetitions::Int = 4
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Volume ratio tolerance required
    to repeat processing the current node (default = 0.9)"
    branch_repetition_tol::Float64 = 0.9

    # Termination limits
    "Maximum number of nodes (default = 1E7)"
    node_limit::Int = 1E7
    "Maximum CPU time in seconds (default = 3600)"
    time_limit::Float64 = 3600.0
    "Maximum number of iterations (default 1E9)"
    iteration_limit::Int = 1E9
    "Absolute tolerance for termination (default = 1E-3)"
    absolute_tolerance::Float64 = 1E-3
    "Relative tolerance for termination (default = 1E-3)"
    relative_tolerance::Float64 = 1E-3
    "Absolute constraint feasibility tolerance (default = 1E-8)"
    absolute_constraint_feas_tolerance::Float64 = 1E-8

    # Options for constraint propagation
    "Depth in B&B tree above which constraint propagation should be disabled (default = 0)"
    cp_depth::Int = 0
    "Number of times to repeat forward-reverse pass routine (default = 0)"
    cp_repetitions::Int = 0
    "Disable constraint propagation if the ratio of new node volume to beginning node volume exceeds
    this number (default = 0.99)"
    cp_tolerance::Float64 = 0.99
    "Use only valid interval bounds during constraint propagation (default = false)"
    cp_interval_only::Bool = false

    # Optimization-Based Bound Tightening (OBBT) options
    "Depth in B&B tree above which OBBT should be disabled (default = 6)"
    obbt_depth::Int = 6
    "Number of repetitions of OBBT to perform in preprocessing (default = 3)"
    obbt_repetitions::Int = 3
    "Turn on aggresive OBBT (default = true)"
    obbt_aggressive_on::Bool = true
    "Maximum iteration to perform aggresive OBBT (default = 2)"
    obbt_aggressive_max_iteration::Int = 2
    "Minimum dimension to perform aggresive OBBT (default = 2)"
    obbt_aggressive_min_dimension::Int = 2
    "Tolerance to consider bounds equal (default = 1E-10)"
    obbt_tolerance::Float64 = 1E-10

    # Feasibility-Based Bound Tightening (FBBT) [linear bound tightening] options
    "Depth in B&B tree above which linear FBBT should be disabled (default = 1000)"
    fbbt_lp_depth::Int = 1000
    "Number of repetitions of linear FBBT to perform in preprocessing (default = 3)"
    fbbt_lp_repetitions::Int = 3

    # Duality-Based Bound Tightening (DBBT) options
    "Depth in B&B tree above which duality-based bound tightening should be disabled (default = 1E10)"
    dbbt_depth::Int = 1E10
    "New bound is considered equal to the prior bound if within dbbt_tolerance (default = 1E-8)"
    dbbt_tolerance::Float64 = 1E-8

    # Subgradient tightening options
    "RelaxTag used to specify type of McCormick operator (default = NS())"
    relax_tag::RelaxTag = NS()
    "Perform tightening of interval bounds using subgradients at each factor in
    each nonlinear tape during a forward pass (default = true)"
    subgrad_tighten::Bool = true
    "Perform tightening of interval bounds using subgradients at each factor in
    each nonlinear tape during a reverse pass (default = false)"
    reverse_subgrad_tighten::Bool = false
    "Outer-round computed subgradient bounds by this amount (default = 1E-10)"
    subgrad_tol::Float64 = 1E-10
    "Select the type of relaxation to use for the bilinear term (multiplication): 0 corresponds to
     a standard McCormick arithmetic approach. Settings 1-3 augment the standard McCormick relaxation
     with implied apriori relaxations: (1) corresponds to a subgradient-based apriori relaxation approach; (2) 
     corresponds to an affine arithmetic-based apriori approach; and (3) corresponds to a enumerative apriori
     relaxation-based approach (default = 0)"
    mul_relax_style::Int = 0

    # Tolerance to add cuts and max number of cuts
    "Minimum number of cuts at each node to attempt (unsafe cuts not necessarily added) (default = 2)"
    cut_min_iterations::Int = 2
    "Maximum number of cuts at each node to attempt (default = 8)"
    cut_max_iterations::Int = 8
    "Absolute tolerance checked for continuing cut (default = 1E-6)"
    cut_tolerance_abs::Float64 = 1E-6
    "Relative tolerance checked for continuing cut (default = 1E-3)"
    cut_tolerance_rel::Float64 = 1E-3

    "Use tolerances to determine safe cuts in a Khajavirad 2018 manner (default = true)"
    cut_safe_on::Bool = true
    "Lower tolerance for safe-lp cut, Khajavirad 2018 (default = 1E-7)"
    cut_safe_l::Float64 = 1E-7
    "Upper tolerance for safe-lp cut, Khajavirad 2018 (default = 1E7)"
    cut_safe_u::Float64 = 1E7
    "Constant tolerance for safe-lp cut, Khajavirad 2018 (default = 1E9)"
    cut_safe_b::Float64 = 1E9

    "Solve upper problem for every node with depth less than `upper_bounding_depth`,
    and otherwise solve upper problems with a probability of `(1/2)^(depth-upper_bounding_depth)`
    (default = 8)"
    upper_bounding_depth::Int = 8

    # Handling for domain violations
    "(Unused) Protect against domain violation (default = false)"
    domain_violation_guard_on::Bool = false
    "(Unused) Amount about a domain violation to ignore when propagating bounds (default = 1E-9)"
    domain_violation_ϵ::Float64 = 1E-9

    "If true, EAGO forgoes its default configuration process for subsolvers (default = false)"
    user_solver_config::Bool = false

    "Absolute tolerance used to check for integrality of decision variables (default = 1E-9)"
    integer_abs_tol::Float64 = 1E-9
    "Relative tolerance used to check for integrality of decision variables (default = 1E-9)"
    integer_rel_tol::Float64 = 1E-9

    # Other forcing options
    "Ignore EAGO's ability to parse problem types and force it to run global optimization (default = false)"
    force_global_solve::Bool = false
    "Check that all branching variables have finite bounds and set them to +/- 1E10 if not (default = true)"
    unbounded_check::Bool = true
end
const EAGO_PARAMETERS = fieldnames(EAGOParameters)

"""
    $(TYPEDEF)

A structure used to hold objectives and constraints added to the EAGO model.
The constraints generally aren't used for relaxations.

All field information available in extended help.

# Extended Help
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct InputProblem

    "Count for the number of variables"
    _variable_count::Int = 0
    "Count for the number of constraints"
    _constraint_count::Int = 0

    # Dictionaries of [ConstraintIndex{F,S} => Tuple{Function, Set}]
    # The function `F` is the type of function in the constraint;
    # The set `S` is the type of set in the constraint
    _vi_leq_constraints::Dict{CI{VI,LT}, Tuple{VI,LT}} = Dict{CI{VI,LT}, Tuple{VI,LT}}()
    _vi_geq_constraints::Dict{CI{VI,GT}, Tuple{VI,GT}} = Dict{CI{VI,GT}, Tuple{VI,GT}}()
    _vi_eq_constraints::Dict{CI{VI,ET}, Tuple{VI,ET}} = Dict{CI{VI,ET}, Tuple{VI,ET}}()
    _vi_it_constraints::Dict{CI{VI,IT}, Tuple{VI,IT}} = Dict{CI{VI,IT}, Tuple{VI,IT}}()
    _vi_zo_constraints::Dict{CI{VI,ZO}, Tuple{VI,ZO}} = Dict{CI{VI,ZO}, Tuple{VI,ZO}}()
    _vi_int_constraints::Dict{CI{VI,MOI.Integer}, Tuple{VI,MOI.Integer}} = Dict{CI{VI,MOI.Integer}, Tuple{VI,MOI.Integer}}()

    _linear_leq_constraints::Dict{CI{SAF,LT}, Tuple{SAF,LT}} = Dict{CI{SAF,LT}, Tuple{SAF,LT}}()
    _linear_geq_constraints::Dict{CI{SAF,GT}, Tuple{SAF,GT}} = Dict{CI{SAF,GT}, Tuple{SAF,GT}}()
    _linear_eq_constraints::Dict{CI{SAF,ET}, Tuple{SAF,ET}} = Dict{CI{SAF,ET}, Tuple{SAF,ET}}()

    _quadratic_leq_constraints::Dict{CI{SQF,LT}, Tuple{SQF,LT}} = Dict{CI{SQF,LT}, Tuple{SQF,LT}}()
    _quadratic_geq_constraints::Dict{CI{SQF,GT}, Tuple{SQF,GT}} = Dict{CI{SQF,GT}, Tuple{SQF,GT}}()
    _quadratic_eq_constraints::Dict{CI{SQF,ET}, Tuple{SQF,ET}} = Dict{CI{SQF,ET}, Tuple{SQF,ET}}()

    _conic_second_order::Dict{CI{VECOFVAR,SOC}, Tuple{VECOFVAR,SOC}} = Dict{CI{VECOFVAR,SOC}, Tuple{VECOFVAR,SOC}}()

    # Dictionaries of [ConstraintIndex{F,S} => Primal Storage]
    # The function `F` is the type of function in the constraint;
    # The set `S` is the type of set in the constraint
    _linear_leq_primal::Dict{CI{SAF,LT},Float64} = Dict{CI{SAF,LT},Float64}()
    _linear_geq_primal::Dict{CI{SAF,GT},Float64} = Dict{CI{SAF,GT},Float64}()
    _linear_eq_primal::Dict{CI{SAF,ET},Float64} = Dict{CI{SAF,ET},Float64}()

    _quadratic_leq_primal::Dict{CI{SQF,LT},Float64} = Dict{CI{SQF,LT},Float64}()
    _quadratic_geq_primal::Dict{CI{SQF,GT},Float64} = Dict{CI{SQF,GT},Float64}()
    _quadratic_eq_primal::Dict{CI{SQF,ET},Float64} = Dict{CI{SQF,ET},Float64}()

    _linear_leq_prob_to_ip::Dict{CI{SAF,LT},CI{SAF,LT}} = Dict{CI{SAF,LT},CI{SAF,LT}}()
    _linear_geq_prob_to_ip::Dict{CI{SAF,GT},CI{SAF,GT}} = Dict{CI{SAF,GT},CI{SAF,GT}}()
    _linear_eq_prob_to_ip::Dict{CI{SAF,ET},CI{SAF,ET}} = Dict{CI{SAF,ET},CI{SAF,ET}}()

    _quadratic_leq_prob_to_ip::Dict{CI{SQF,LT},CI{SQF,LT}} = Dict{CI{SQF,LT},CI{SQF,LT}}()
    _quadratic_geq_prob_to_ip::Dict{CI{SQF,GT},CI{SQF,GT}} = Dict{CI{SQF,GT},CI{SQF,GT}}()
    _quadratic_eq_prob_to_ip::Dict{CI{SQF,ET},CI{SQF,ET}} = Dict{CI{SQF,ET},CI{SQF,ET}}()

    "Storage for the objective function"
    _objective::Union{VI,SAF,SQF,Nothing} = nothing

    "Storage for NLP constraints (set by `MOI.set(m, ::NLPBlockData...)` in `moi_wrapper.jl`)"
    _nlp_data::Union{MOI.NLPBlockData,Nothing} = nothing

    "Objective sense information (set by `MOI.set(m, ::ObjectiveSense...)`)"
    _optimization_sense::MOI.OptimizationSense = MOI.MIN_SENSE
end

function MOI.empty!(ip::InputProblem)

    # Empty all arrays and dictionaries in input problem
    for field in fieldnames(InputProblem)
        field_value = getfield(ip, field)
        if (field_value isa Array) || (field_value isa Dict)
            empty!(field_value)
        end
    end

    # Reset all non-array fields and dictionaries in input problem to their original values
    ip._variable_count = 0
    ip._constraint_count = 0
    ip._objective = nothing
    ip._nlp_data = nothing
    ip._optimization_sense = MOI.MIN_SENSE
    return
end

function Base.isempty(x::InputProblem)

    is_empty_flag = true
    new_input_problem = InputProblem()

    for field in fieldnames(InputProblem)
        field_value = getfield(x, field)
        if (field_value isa Array) || (field_value isa Dict)
            if !isempty(field_value)
                is_empty_flag = false
                break
            end
        elseif field_value isa Number
            if getfield(new_input_problem, field) != field_value
                is_empty_flag = false
                break
            end
        end
    end
    is_empty_flag &= isnothing(x._nlp_data)
    is_empty_flag &= isnothing(x._objective)
    return is_empty_flag
end

# Helper functions which simplify finding constraints of different types
_constraints(m::InputProblem, ::Type{VI}, ::Type{LT}) = m._vi_leq_constraints
_constraints(m::InputProblem, ::Type{VI}, ::Type{GT}) = m._vi_geq_constraints
_constraints(m::InputProblem, ::Type{VI}, ::Type{ET}) = m._vi_eq_constraints
_constraints(m::InputProblem, ::Type{VI}, ::Type{IT}) = m._vi_it_constraints
_constraints(m::InputProblem, ::Type{VI}, ::Type{ZO}) = m._vi_zo_constraints
_constraints(m::InputProblem, ::Type{VI}, ::Type{MOI.Integer}) = m._vi_int_constraints

_constraints(m::InputProblem, ::Type{SAF}, ::Type{LT}) = m._linear_leq_constraints
_constraints(m::InputProblem, ::Type{SAF}, ::Type{GT}) = m._linear_geq_constraints
_constraints(m::InputProblem, ::Type{SAF}, ::Type{ET}) = m._linear_eq_constraints

_constraints(m::InputProblem, ::Type{SQF}, ::Type{LT}) = m._quadratic_leq_constraints
_constraints(m::InputProblem, ::Type{SQF}, ::Type{GT}) = m._quadratic_geq_constraints
_constraints(m::InputProblem, ::Type{SQF}, ::Type{ET}) = m._quadratic_eq_constraints

_constraint_primal(m::InputProblem, ::Type{SAF}, ::Type{LT}) = m._linear_leq_primal
_constraint_primal(m::InputProblem, ::Type{SAF}, ::Type{GT}) = m._linear_geq_primal
_constraint_primal(m::InputProblem, ::Type{SAF}, ::Type{ET}) = m._linear_eq_primal

_constraint_primal(m::InputProblem, ::Type{SQF}, ::Type{LT}) = m._quadratic_leq_primal
_constraint_primal(m::InputProblem, ::Type{SQF}, ::Type{GT}) = m._quadratic_geq_primal
_constraint_primal(m::InputProblem, ::Type{SQF}, ::Type{ET}) = m._quadratic_eq_primal

"""
    $(TYPEDSIGNATURES)

Extract primal constraint value from the local problem and save the result to the appropriate
field of the optimizer.
"""
function _extract_primal!(d, m::InputProblem, ::Type{F}, ::Type{S}) where {F,S}
    for (k, v) in _constraint_index_to_ip(m, F, S)
        _constraint_primal(m, F, S)[v] = MOI.get(d, MOI.ConstraintPrimal(), k)
    end
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Extract linear primal constraint values from the local problem and save the result to the
appropriate field of the optimizer.
"""
function _extract_primal_linear!(d, ip::InputProblem)
    _extract_primal!(d, ip, SAF, LT)
    _extract_primal!(d, ip, SAF, GT)
    _extract_primal!(d, ip, SAF, ET)
end
"""
    $(TYPEDSIGNATURES)

Extract quadratic primal constraint values from the local problem and save the result to the
appropriate field of the optimizer.
"""
function _extract_primal_quadratic!(d, ip::InputProblem)
    _extract_primal!(d, ip, SQF, LT)
    _extract_primal!(d, ip, SQF, GT)
    _extract_primal!(d, ip, SQF, ET)
end

_constraint_index_to_ip(m::InputProblem, ::Type{SAF}, ::Type{LT}) = m._linear_leq_prob_to_ip
_constraint_index_to_ip(m::InputProblem, ::Type{SAF}, ::Type{GT}) = m._linear_geq_prob_to_ip
_constraint_index_to_ip(m::InputProblem, ::Type{SAF}, ::Type{ET}) = m._linear_eq_prob_to_ip

_constraint_index_to_ip(m::InputProblem, ::Type{SQF}, ::Type{LT}) = m._quadratic_leq_prob_to_ip
_constraint_index_to_ip(m::InputProblem, ::Type{SQF}, ::Type{GT}) = m._quadratic_geq_prob_to_ip
_constraint_index_to_ip(m::InputProblem, ::Type{SQF}, ::Type{ET}) = m._quadratic_eq_prob_to_ip

"""
    $(TYPEDSIGNATURES)

Add a constraint to the local problem, storing the new constraint index and the associated
index in the input problem.
"""
function _add_constraint_store_ci!(d, m::InputProblem, ::Type{F}, ::Type{S}) where {F,S}
    for (ci_ip, fs) in _constraints(m, F, S)
        ci_wp = MOI.add_constraint(d, fs[1], fs[2])
        _constraint_index_to_ip(m, F, S)[ci_wp] = ci_ip
    end
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Add linear constraints to the local problem, storing the new constraint indices and the associated
indices in the input problem.
"""
function _add_constraint_store_ci_linear!(d, ip::InputProblem)
    _add_constraint_store_ci!(d, ip, SAF, LT)
    _add_constraint_store_ci!(d, ip, SAF, GT)
    _add_constraint_store_ci!(d, ip, SAF, ET)
end

"""
    $(TYPEDSIGNATURES)

Add quadratic constraints to the local problem, storing the new constraint indices and the associated
indices in the input problem.
"""
function _add_constraint_store_ci_quadratic!(d, ip::InputProblem)
    _add_constraint_store_ci!(d, ip, SQF, LT)
    _add_constraint_store_ci!(d, ip, SQF, GT)
    _add_constraint_store_ci!(d, ip, SQF, ET)
end

"""
    $(TYPEDEF)

A structure used to store expressions and problem descriptions EAGO uses to formulate
relaxed problems.

All field information available in extended help.

# Extended Help
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct ParsedProblem

    "Problem classification (set in `parse_classify_problem!`)"
    _problem_type::ANY_PROBLEM_TYPE = nothing

    "Stores the objective and is used for constructing linear affine cuts"
    _objective_saf::SAF = SAF(SAT[], 0.0)
    "Storage for the objective function"
    _objective::Union{VI,AffineFunctionIneq,BufferedQuadraticIneq,BufferedNonlinearFunction,Nothing} = nothing

    "Objective sense information (set by `MOI.set(m, ::ObjectiveSense...)`)"
    _optimization_sense::MOI.OptimizationSense = MOI.MIN_SENSE

    # Non-single variable constraints (set in `initial_parse!`)
    _saf_leq::Vector{AffineFunctionIneq} = AffineFunctionIneq[]
    _saf_eq::Vector{AffineFunctionEq} = AffineFunctionEq[]
    _sqf_leq::Vector{BufferedQuadraticIneq} = BufferedQuadraticIneq[]
    _sqf_eq::Vector{BufferedQuadraticEq} = BufferedQuadraticEq[]
    _conic_second_order::Vector{BufferedSOC} = BufferedSOC[]

    # NLP constraints
    _nlp_data::Union{MOI.NLPBlockData,Nothing} = nothing
    _nonlinear_constr::Vector{BufferedNonlinearFunction} = BufferedNonlinearFunction[]
    _relaxed_evaluator::Evaluator = Evaluator()

    "Variable information (set in `initial_parse!`)"
    _variable_info::Vector{VariableInfo{Float64}} = VariableInfo{Float64}[]
    "Count for the number of variables"
    _variable_count::Int = 0
end

function MOI.empty!(x::ParsedProblem)

    for field in fieldnames(ParsedProblem)
        field_value = getfield(x, field)
        if (field_value isa Array) || (field_value isa Dict)
            empty!(field_value)
        end
    end

    x._objective    = nothing
    x._problem_type = nothing
    x._nlp_data     = nothing

    x._optimization_sense = MOI.MIN_SENSE
    x._relaxed_evaluator = Evaluator()
    x._objective_saf = SAF(SAT[], 0.0)
    x._variable_count = 0
    return
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
            if getfield(new_input_problem, field) != field_value
                is_empty_flag = false
                break
            end
        end
    end

    is_empty_flag &= isempty(x._objective_saf.terms)
    is_empty_flag &= iszero(x._objective_saf.constant)
    is_empty_flag &= isnothing(x._objective)
    is_empty_flag &= isnothing(x._nlp_data)

    return is_empty_flag
end

"""
    $(TYPEDEF)

Optimizer internal to EAGO which holds information used to perform
branch-and-bound in order to solve nonconvex MINLPs.

Descriptions of all fields available in extended help.

# Extended Help
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct GlobalOptimizer{Q,S,T<:ExtensionType} <: MOI.AbstractOptimizer
    
    "Storage for relaxed and upper optimizers to use, and any custom extensions"
    _subsolvers::SubSolvers{Q,S,T} = SubSolvers{Q,S,T}()
    "Parameters that do not change during a global solve"
    _parameters::EAGOParameters = EAGOParameters()
    "Expressions and constraints added to the EAGO model (not directly
    used for relaxations)"
    _input_problem::InputProblem = InputProblem()
    "Expressions and problem descriptions that EAGO uses to formulate
    relaxed problems"
    _working_problem::ParsedProblem = ParsedProblem()
    "Information on any auxiliary variables"
    _auxiliary_variable_info::Union{Nothing,_AuxVarData} = nothing
    "Variables to perform OBBT on (default: all variables in nonlinear expressions)"
    obbt_variable_values::Vector{Bool} = Bool[]
    "Specifies that the optimize_hook! function should be called rather than throw the
    problem to the standard routine"
    enable_optimize_hook::Bool = false
    "(Deprecated, use _subsolvers instead) Storage for custom extension types" 
    ext

    # Set as user-specified option
    "The completion status code for the branch-and-bound algorithm"
    _end_state::GlobalEndState = GS_UNSET
    "The MathOptInterface-compliant completion status code"
    _termination_status_code::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    "Value indicating the feasibility status of the result"
    _result_status_code::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    "Multiplier used internally to convert objective sense from Max to Min. Only
    takes on values of {-1.0, 1.0}"
    _obj_mult::Float64 = 1.0
    "Flag to indicate if a slack variable was added for the objective function. This
    is done in some epigraph reformulations (see [`reform_epigraph_min!`](@ref))"
    _obj_var_slack_added::Bool = false

    "A heap of all nodes in the branch-and-bound tree"
    _stack::BinaryMinMaxHeap{NodeBB} = BinaryMinMaxHeap{NodeBB}()
    "The individual node being examined at any particular time. Nodes are removed from
    the stack and placed here, evaluated, and then sent back to the stack"
    _current_node::NodeBB = NodeBB()

    "(Unused) Flag for relaxation points"
    _first_relax_point_set::Bool = false
    "(Unused) Variable values of a particular point"
    _current_xref::Vector{Float64} = Float64[]
    "(Unused) Variable values of a candidate point"
    _candidate_xref::Vector{Float64} = Float64[]

    "(Unused) Flag to use variable values from previous evaluation on the current step"
    _use_prior_objective_xref::Bool = false
    "(Unused) Variable values for objective evaluation"
    _current_objective_xref::Vector{Float64} = Float64[]
    "(Unused) Variable values for previous objective evaluation"
    _prior_objective_xref::Vector{Float64} = Float64[]

    "Flag for if the user has specified branch variables (see [`label_branch_variables!`](@ref))"
    _user_branch_variables::Bool = false
    "Variables that are fixed in place"
    _fixed_variable::Vector{Bool} = Bool[]
    "Number of variables that can be branched on"
    _branch_variable_count::Int = 0
    "Mapping from the branch variables to the full set of variables in the problem"
    _branch_to_sol_map::Vector{Int} = Int[]
    "Mapping from the full set of variables in the problem to the branch variables"
    _sol_to_branch_map::Vector{Int} = Int[]

    "The final (or intermediate) variable values of the solution"
    _continuous_solution::Vector{Float64} = Float64[]

    "Flag to ensure preprocessing result is feasible"
    _preprocess_feasibility::Bool = true
    "Status codes for use in bounds tightening"
    _preprocess_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    "Status codes for use in bounds tightening"
    _preprocess_primal_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    "Status codes for use in bounds tightening"
    _preprocess_dual_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS

    "Primal status of the lower problem"
    _lower_primal_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    "Dual status of the lower problem"
    _lower_dual_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    "Termination status of the lower problem"
    _lower_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    "Flag for lower problem feasibility"
    _lower_feasibility::Bool = true
    "Objective value result from the lower problem"
    _lower_objective_value::Float64 = -Inf
    "Variable values of the lower problem solution"
    _lower_solution::Vector{Float64} = Float64[]
    "Lower variable duals for use in duality-based bound tightening"
    _lower_lvd::Vector{Float64} = Float64[]
    "Upper variable duals for use in duality-based bound tightening"
    _lower_uvd::Vector{Float64} = Float64[]

    "Objective value associated with the previous cut in the cutting planes algorithm"
    _last_cut_objective::Float64 = -Inf

    "Primal status of the upper problem"
    _upper_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
    "Termination status of the upper problem"
    _upper_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    "Flag for upper problem feasibility"
    _upper_feasibility::Bool = true
    "Objective value result from the upper problem"
    _upper_objective_value::Float64 = Inf
    ""
    _upper_variables::Vector{VI} =  VI[]
    ""
    _upper_solution::Vector{Float64} = Float64[]

    "(Unused) Flag to ensure postprocessing result is feasible"
    _postprocess_feasibility::Bool = true

    # Set to time limit in initial_parse! in parse.jl, decremented throughout global_solve in optimize_nonconvex.jl
    "Time remaining for the optimization algorithm. This is set in `initial_parse!` to the user-defined
    time limit and is decremented throughout `global_solve!`"
    _time_left::Float64 = 3600.0

    # Set constructor reset on empty! and to zero in initial parse! in parse.jl
    "Storage for the `time()` when optimization began"
    _start_time::Float64 = 0.0
    "Current run time, incremented using `time()-_start_time`"
    _run_time::Float64 = 0.0
    "A field to keep track of time spent on initial problem parsing"
    _parse_time::Float64 = 0.0
    "Used in `optimize_nonconvex.jl` to track how long the presolve step takes"
    _presolve_time::Float64 = 0.0
    "Updated each iteration to track the time of the preprocess step"
    _last_preprocess_time::Float64 = 0.0
    "Updated each iteration to track the time of the lower problem step"
    _last_lower_problem_time::Float64 = 0.0
    "Updated each iteration to track the time of the upper problem step"
    _last_upper_problem_time::Float64 = 0.0
    "Updated each iteration to track the time of the postprocess step"
    _last_postprocessing_time::Float64 = 0.0

    # Reset in initial_parse! in parse.jl
    "A field to track convergence progress across iterations"
    _min_converged_value::Float64 = Inf
    "The best-known lower bound"
    _global_lower_bound::Float64 = -Inf
    "The best-known upper bound"
    _global_upper_bound::Float64 = Inf
    "The total number of nodes that have been created"
    _maximum_node_id::Int = 0
    "The number of iterations the branch-and-bound algorithm has completed"
    _iteration_count::Int = 0
    "The number of nodes in the stack"
    _node_count::Int = 0

    # Storage for output, reset in initial_parse! in parse.jl
    "(Unused) The best-known solution value"
    _solution_value::Float64 = 0.0
    "A flag for if a feasible solution was identified. Updated if preprocessing,
    lower problem, and upper problem all return feasible values"
    _feasible_solution_found::Bool = false
    "The node ID of the best-known feasible upper problem solution (default = -1, if no feasible solution is found)"
    _solution_node::Int = -1
    "The best-known upper bound"
    _best_upper_value::Float64 = Inf #TODO: Duplicate of _global_upper_bound? Why is this ever used instead?

    # Optimality-Based Bound Tightening (OBBT) options
    "Indices of variables to perform OBBT on"
    _obbt_working_lower_index::Vector{Bool} = Bool[]
    "Indices of variables to perform OBBT on"
    _obbt_working_upper_index::Vector{Bool} = Bool[]
    "Tracker for changes in _obbt_working_lower_index across iterations"
    _lower_indx_diff::Vector{Bool} = Bool[]
    "Tracker for changes in _obbt_working_upper_index across iterations"
    _upper_indx_diff::Vector{Bool} = Bool[]
    "Storage for indices prior to OBBT step"
    _old_low_index::Vector{Bool} = Bool[]
    "Storage for indices prior to OBBT step"
    _old_upp_index::Vector{Bool} = Bool[]
    "New indices following OBBT step; compared with `_old_low_index`"
    _new_low_index::Vector{Bool} = Bool[]
    "New indices following OBBT step; compared with `_old_upp_index`"
    _new_upp_index::Vector{Bool} = Bool[]
    "(Deprecated) Variables to perform OBBT on. Replaced by `_obbt_working_lower_index` and `_obbt_working_upper_index`"
    _obbt_variables::Vector{VI} = VI[]
    "The number of variables to perform OBBT on"
    _obbt_variable_count::Int = 0
    "(Unused) Flag to indicate whether OBBT has been performed"
    _obbt_performed_flag::Bool = false

    # Buffers for FBBT, set in presolve, used in preprocess
    "Buffer for FBBT lower bounds. Set in presolve, used in preprocess"
    _lower_fbbt_buffer::Vector{Float64} = Float64[]
    "Buffer for FBBT upper bounds. Set in presolve, used in preprocess"
    _upper_fbbt_buffer::Vector{Float64} = Float64[]
    "(Unused) Improvement in constraint propagation"
    _cp_improvement::Float64 = 0.0
    "(Unused) Flag for if constraint propagation results need to be reversed"
    _cp_evaluation_reverse::Bool = false

    "Iterations of the cutting planes algorithm completed"
    _cut_iterations::Int = 0
    "(Unused) Flag to check if cuts should be added"
    _cut_add_flag::Bool = false
    "Counter for number of times a node is evaluated. If the `repeat_check` function
    is overloaded to return `true`, a node will not be branched on, but will instead
    be added back into the stack using `single_storage!`. In this case, `_node_repetitions`
    is incremented"
    _node_repetitions::Int = 0

    "Storage for logging information during a branch-and-bound run"
    _log::Log = Log()

    "Storage for affine constraints"
    _affine_relax_ci::Vector{CI{SAF,LT}} = CI{SAF,LT}[]
    "Storage for a linear objective cut constraint"
    _affine_objective_cut_ci::Union{CI{VI,LT},CI{SAF,LT},Nothing} = nothing

    "(Unused) Number of relaxed variables"
    _relaxed_variable_number::Int = 0
    "Indices of relaxed variables"
    _relaxed_variable_index::Vector{VI} = VI[]
    "Stored EqualTo constraints"
    _relaxed_variable_et::Vector{Tuple{CI{VI, ET}, Int}} = Tuple{CI{VI, ET}, Int}[]
    "Stored LessThan constraints"
    _relaxed_variable_lt::Vector{Tuple{CI{VI, LT}, Int}} = Tuple{CI{VI, LT}, Int}[]
    "Stored GreaterThan constraints"
    _relaxed_variable_gt::Vector{Tuple{CI{VI, GT}, Int}} = Tuple{CI{VI, GT}, Int}[]
    "Stored integer constraints"
    _relaxed_variable_integer::Vector{CI{VI, MOI.Integer}} = CI{VI, MOI.Integer}[]

    "List of variables that can be branched on. If not user-specified, branch variables
    are identified in `label_branch_variables!`"
    _branch_variables::Vector{Bool} = Bool[]
    "(Unused) Flag for non-branching integers"
    _nonbranching_int::Bool = false

    "Flag indicating if an initial evaluation of the constraints has occurred"
    _new_eval_constraint::Bool = false
    "Flag indicating if the objective expression was evaluated"
    _new_eval_objective::Bool = false

    "Storage for carrying LessThan constraint information. Used in `obbt!` and
    `update_relaxed_problem_box!`"
    _node_to_sv_leq_ci::Dict{Int,CI{VI,LT}} = Dict{Int,CI{VI,LT}}()
    "Storage for carrying GreaterThan constraint information. Used in `obbt!` and
    `update_relaxed_problem_box!`"
    _node_to_sv_geq_ci::Dict{Int,CI{VI,GT}} = Dict{Int,CI{VI,GT}}()
    "Flag to check for nonlinear evaluators. Set to `true` in `add_nonlinear_evaluator!`"
    _nonlinear_evaluator_created::Bool = false

    "(FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED) Storage for pseudocost branching"
    _branch_cost::BranchCostStorage{Float64} = BranchCostStorage{Float64}()
    "(FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED) Sparsity information of the
    branch variables"
    _branch_variable_sparsity::SparseMatrixCSC{Bool,Int} = spzeros(Bool,1,1)
    "(FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED) Information on the infeasibility
    of each constraint"
    _constraint_infeasiblity::Vector{Float64} = Float64[]
end

const EAGO_OPTIMIZER_ATTRIBUTES = Symbol[:relaxed_optimizer, :upper_optimizer,
                                         :enable_optimize_hook, :ext, :_parameters]
const EAGO_MODEL_STRUCT_ATTRIBUTES = Symbol[:_stack, :_log, :_current_node, :_working_problem, :_input_problem, :_branch_cost]
const EAGO_MODEL_NOT_STRUCT_ATTRIBUTES = setdiff(fieldnames(GlobalOptimizer), union(EAGO_OPTIMIZER_ATTRIBUTES, EAGO_MODEL_STRUCT_ATTRIBUTES))
                           
function MOI.empty!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q}
                           
    # Create a new empty optimizer and copy fields to m
    new_optimizer = GlobalOptimizer{R,S,Q}(_subsolvers = m._subsolvers,
                                           _parameters = m._parameters, 
                                           _input_problem = m._input_problem,
                                           _working_problem = m._working_problem,
                                           ext = nothing)

    MOI.empty!(new_optimizer._subsolvers)
    MOI.empty!(new_optimizer._input_problem)
    MOI.empty!(new_optimizer._working_problem)
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
                                           _working_problem = m._working_problem,
                                           ext = nothing)
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
_ext(m::GlobalOptimizer{R,S,Q})               where {R,S,Q} = m._subsolvers.ext

"""
$(TYPEDSIGNATURES)

Check to see if the sense of the optimization problem is Min (`true`)
or Max (`false`).
"""
@inline function _is_input_min(m::GlobalOptimizer)
    return m._obj_mult == 1.0
end
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

"""
$(TYPEDSIGNATURES)

Return the branch-to-sol mapping for the i'th variable.
"""
@inline _bvi(m::GlobalOptimizer, i::Int) = m._branch_to_sol_map[i]
"""
$(TYPEDSIGNATURES)

Return the sol-to-branch mapping for the i'th variable.
"""
@inline _svi(m::GlobalOptimizer, i::Int) = m._sol_to_branch_map[i]

"""
$(TYPEDSIGNATURES)

Check if `m._branch_variables[i]` is true or false.
"""
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

"""
    _diam(::BranchVar, ::GlobalOptimizer, ::Int)
    _diam(::FullVar, ::GlobalOptimizer, ::Int)

Return the diameter of a variable (upper bound - lower bound).
"""
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

"""
    _constraint_primal

Helper function which simplifies finding constraints of different types. See also: [`_constraints`](@ref)
"""
_constraint_primal(m::GlobalOptimizer, ::Type{F}, ::Type{S}) where {F,S} = _constraint_primal(m._input_problem, F, S)