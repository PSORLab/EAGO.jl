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

struct EmptyNLPEvaluator <: MOI.AbstractNLPEvaluator
    _current_node::NodeBB
    has_nlobj::Bool
end
EmptyNLPEvaluator() = EmptyNLPEvaluator(NodeBB(),false)
set_current_node!(x::EmptyNLPEvaluator, n::NodeBB) = ()

MOI.features_available(::EmptyNLPEvaluator) = [:Grad, :Jac, :Hess]
MOI.initialize(::EmptyNLPEvaluator, features) = nothing
MOI.eval_objective(::EmptyNLPEvaluator, x) = NaN
function MOI.eval_constraint(::EmptyNLPEvaluator, g, x)
    @assert length(g) == 0
    return
end
MOI.eval_objective_gradient(::EmptyNLPEvaluator, g, x) = nothing
MOI.jacobian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
function MOI.eval_constraint_jacobian(::EmptyNLPEvaluator, J, x)
    @assert length(J) == 0
    return
end
function MOI.eval_hessian_lagrangian(::EmptyNLPEvaluator, H, x, σ, μ)
    @assert length(H) == 0
    return
end

empty_nlp_data() = MOI.NLPBlockData([], EmptyNLPEvaluator(), false)

#=
struct BufferedSQF{S}
    sqf::SQF
    set::S
    buffer::SAF
    quad1_to_indx::Vector{Int}
    quad2_to_indx::Vector{Int}
    affine_to_indx::Vector{Int}
end
=#


export Optimizer
"""
$(TYPEDEF)

The main optimizer object used by EAGO to solve problems during the optimization
routine. The following commonly used options are described below and can be set
via keyword arguments in the JuMP/MOI model:

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Optimizer{S<:MOI.AbstractOptimizer, T<:MOI.AbstractOptimizer} <: MOI.AbstractOptimizer

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
    # "Attempt to bridge convex constraint to second order cone"
    # conic_convert_quadratic::Bool

    # Options for constraint propagation
    "Depth in B&B tree above which constraint propagation should be disabled (default = 1000)"
    cp_depth::Int64 = 1000
    "Number of repetitions of forward-reverse passes to perform in constraint propagation (default = 3)"
    cp_repetitions::Int64 = 3
    "Disable constraint propagation if the ratio of new node volume to beginning node volume exceeds
    this number (default = 0.99)"
    cp_tolerance::Float64 = 0.99
    "Use only valid interval bounds during constraint propagation (default = false)"
    cp_interval_only::Bool = false

    # Options for optimality-based bound tightening
    "An instance of the optimizer used to solve the relaxed subproblems (default = GLPK.Optimizer())"
    relaxed_optimizer::S = GLPK.Optimizer()
    "Keyword arguments for the relaxed optimizer."
    relaxed_optimizer_kwargs::Base.Iterators.Pairs = Base.Iterators.Pairs(NamedTuple(),())
    relaxed_inplace_mod::Bool = true

    "Depth in B&B tree above which OBBT should be disabled (default = 6)"
    obbt_depth::Int64 = 6
    "Number of repetitions of OBBT to perform in preprocessing (default = 3)"
    obbt_repetitions::Int64 = 4
    "Turn aggresive OBBT on (default = false)"
    obbt_aggressive_on::Bool = false
    "Maximum iteration to perform aggresive OBBT (default = 2)"
    obbt_aggressive_max_iteration::Int64 = 2
    "Minimum dimension to perform aggresive OBBT (default = 2)"
    obbt_aggressive_min_dimension::Int64 = 2
    "Tolerance to consider bounds equal (default = 1E-9)"
    obbt_tolerance::Float64 = 1E-9
    "Variables to perform OBBT on (default: all variables in nonlinear expressions)."
    obbt_variable_values::Vector{Bool} = Bool[]

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

    # Subgradient tightening flag
    "Perform tightening of interval bounds using subgradients at each factor in
    each nonlinear tape during a forward-reverse pass (default = true)."
    subgrad_tighten::Bool = true
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Used to enable/disable subgradient
    tightening of interval bounds on the reverse pass (default = true)"
    subgrad_tighten_reverse::Bool = false

    # Tolerance to add cuts and max number of cuts
    cut_max_iterations::Int64 = 3
    "Convex coefficient used to select point for new added cuts. Branch point is
    given by `(1-cut_cvx)*xmid + cut_cvx*xsol` (default = 0.9)."
    cut_cvx::Float64 = 0.9
    "Add cut if the L1 distance from the prior cutting point to the new cutting
    point normalized by the box volume is greater than the tolerance (default = 0.05)."
    cut_tolerance::Float64 = 0.05
    "Adds an objective cut to the relaxed problem (default = true)."
    objective_cut_on::Bool = true
    cut_safe_l::Float64 = 0.00001
    cut_safe_u::Float64 = 10000.0
    cut_safe_b::Float64 = 100000.0

    # Upper bounding options
    upper_optimizer::T = Ipopt.Optimizer()
    upper_factory::JuMP.OptimizerFactory = with_optimizer(Ipopt.Optimizer)
    "Solve upper problem for every node with depth less than `upper_bounding_depth`
    and with a probabilityof (1/2)^(depth-upper_bounding_depth) otherwise (default = 4)"
    upper_bounding_depth::Int64 = 4

    # Duality-based bound tightening (DBBT) options
    "Depth in B&B tree above which duality-based bound tightening should be disabled (default = 1E10)"
    dbbt_depth::Int64 = 10^10
    "New bound is considered equal to the prior bound if within dbbt_tolerance (default = 1E-9)."
    dbbt_tolerance::Float64 = 1E-8

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
    node_limit::Int64 = 1E-7
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

    # Debug
    "Specifies that the optimize_hook! function should be called rather than
    throw the problem to the standard B&B routine (default = false)."
    enable_optimize_hook::Bool = false
    "Holds additional storage needed for constructing extensions to EAGO
    (default = Dict{Symbol,Any})."
    ext::Dict{Symbol, Any} = Dict{Symbol,Any}()
    "Holds an instance of a subtype of `EAGO.ExtensionType` used to define
    new custom subroutines (default = DefaultExt())."
    ext_type::ExtensionType = DefaultExt()

    # handling for domain violations
    domain_violation_ϵ::Float64 = 1E-9

    _current_node::NodeBB = NodeBB()
    _current_xref::Vector{Float64} = Float64[]
    _sense_multiplier::Float64

    _variable_number::Int64 = 0
    _state_variables::Int64
    _continuous_variable_number::Int64 = 0
    _integer_variable_number::Int64

    _user_branch_variables::Bool
    _fixed_variable::Vector{Bool} = Bool[]

    _continuous_solution::Vector{Float64} = Float64[]

    _integer_variables::Vector{Int64} = Int64[]
    _variable_info::Vector{VariableInfo} = VariableInfo[]
    _upper_variables::Vector{VI} =  VI[]

    _stack::BinaryMinMaxHeap{NodeBB} = BinaryMinMaxHeap{NodeBB}()

    _lower_variable::Vector{SV} = SV[]
    _lower_variable_index::Vector{VI} = VI[]
    _lower_variable_et::Vector{CI{SV, ET}} = CI{SV, ET}[]
    _lower_variable_lt::Vector{CI{SV, LT}} = CI{SV, LT}[]
    _lower_variable_gt::Vector{CI{SV, GT}} = CI{SV, GT}[]
    _lower_variable_et_indx::Vector{Int64} = Int64[]
    _lower_variable_lt_indx::Vector{Int64} = Int64[]
    _lower_variable_gt_indx::Vector{Int64} = Int64[]

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

    _best_upper_value::Float64 = Inf

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

    _objective::Union{Nothing, SV, SAF, SQF} = nothing

    _objective_convexity::Bool = false

    _objective_cut_set::Int64 = -1
    _objective_cut_ci_sv::CI{SV,LT} = CI{SV,LT}(-1.0)
    _objective_cut_ci_saf::Vector{CI{SAF,LT}} = CI{SAF,LT}[]

    #_conic_norm_infinity::Vector{Tuple{VecOfVar, NormInfinityCone}} = Tuple{VecOfVar, NormInfinityCone}[]
    #_conic_norm_one::Vector{Tuple{VecOfVar, NormOneCone}} = Tuple{VecOfVar, NormOneCone}[]
    #_conic_second_order::Vector{Tuple{VecOfVar, SecondOrderCone}} = Tuple{VecOfVar, SecondOrderCone}[]
    #_conic_rotated_second_order::Vector{Tuple{VecOfVar, RotatedSecondOrderCone}} = Tuple{VecOfVar, RotatedSecondOrderCone}[]
    #_conic_geometric_mean::Vector{Tuple{VecOfVar, GeometricMeanCone}} = Tuple{VecOfVar, GeometricMeanCone}[]
    #_conic_exponential::Vector{Tuple{VecOfVar, ExponentialCone}} = Tuple{VecOfVar, ExponentialCone}[]
    #_conic_dual_exponential::Vector{Tuple{VecOfVar, DualExponentialCone}} = Tuple{VecOfVar, DualExponentialCone}[]
    #_conic_power_cone::Vector{Tuple{VecOfVar, PowerCone}} = Tuple{VecOfVar, PowerCone}[]
    #_conic_dual_power::Vector{Tuple{VecOfVar, DualPowerCone}} = Tuple{VecOfVar, DualPowerCone}[]
    #_conic_relative_entropy::Vector{Tuple{VecOfVar, RelativeEntropyCone}} = Tuple{VecOfVar, RelativeEntropyCone}[]
    #_conic_norm_spectral::Vector{Tuple{VecOfVar, NormSpectralCone}} = Tuple{VecOfVar, NormSpectralCone}[]
    #_conic_norm_nuclear::Vector{Tuple{VecOfVar, NormNuclearCone}} = Tuple{VecOfVar, NormNuclearCone}[]

    _linear_leq_constraints::Vector{Tuple{SAF, LT, Int64}} = Tuple{SAF, LT, Int64}[]
    _linear_geq_constraints::Vector{Tuple{SAF, GT, Int64}} = Tuple{SAF, GT, Int64}[]
    _linear_eq_constraints::Vector{Tuple{SAF, ET, Int64}} = Tuple{SAF, ET, Int64}[]

    _quadratic_leq_constraints::Vector{Tuple{SQF, LT, Int64}} = Tuple{SQF, LT, Int64}[]
    _quadratic_geq_constraints::Vector{Tuple{SQF, GT, Int64}} = Tuple{SQF, GT, Int64}[]
    _quadratic_eq_constraints::Vector{Tuple{SQF, ET, Int64}} = Tuple{SQF, ET, Int64}[]

    _quadratic_leq_dict::Vector{ImmutableDict{Int64,Int64}} = ImmutableDict{Int64,Int64}[]
    _quadratic_geq_dict::Vector{ImmutableDict{Int64,Int64}} = ImmutableDict{Int64,Int64}[]
    _quadratic_eq_dict::Vector{ImmutableDict{Int64,Int64}} = ImmutableDict{Int64,Int64}[]
    _quadratic_obj_dict::ImmutableDict{Int64,Int64} = ImmutableDict{Int64,Int64}()

    _quadratic_ci_leq::Vector{Vector{CI{SAF,LT}}} = CI{SAF,LT}[]
    _quadratic_ci_geq::Vector{Vector{CI{SAF,LT}}} = CI{SAF,LT}[]
    _quadratic_ci_eq::Vector{Vector{Tuple{CI{SAF,LT},CI{SAF,LT}}}} = Tuple{CI{SAF,LT},CI{SAF,LT}}[]

    _quadratic_leq_sparsity::Vector{Vector{VI}} = Vector{VI}[]
    _quadratic_geq_sparsity::Vector{Vector{VI}} = Vector{VI}[]
    _quadratic_eq_sparsity::Vector{Vector{VI}} = Vector{VI}[]

    _quadratic_leq_gradnz::Vector{Int64} = Int64[]
    _quadratic_geq_gradnz::Vector{Int64} = Int64[]
    _quadratic_eq_gradnz::Vector{Int64} = Int64[]

    _quadratic_leq_convexity::Vector{Bool} = Bool[]
    _quadratic_geq_convexity::Vector{Bool} = Bool[]
    _quadratic_eq_convexity_1::Vector{Bool} = Bool[]
    _quadratic_eq_convexity_2::Vector{Bool} = Bool[]

    _lower_nlp_affine::Vector{Vector{CI{SAF,LT}}} = Vector{CI{SAF,LT}}[]
    _upper_nlp_affine::Vector{Vector{CI{SAF,LT}}} = Vector{CI{SAF,LT}}[]

    _lower_nlp_affine_indx::Vector{Int64} = Int64[]
    _upper_nlp_affine_indx::Vector{Int64} = Int64[]

    _lower_nlp_sparsity::Vector{Vector{Int64}} = Vector{Int64}[]
    _upper_nlp_sparsity::Vector{Vector{Int64}} = Vector{Int64}[]

    _univariate_quadratic_leq_constraints::Vector{Tuple{Float64,Float64,Float64,Int64}} = Tuple{Float64,Float64,Float64,Int64}[]
    _univariate_quadratic_geq_constraints::Vector{Tuple{Float64,Float64,Float64,Int64}} = Tuple{Float64,Float64,Float64,Int64}[]
    _univariate_quadratic_eq_constraints::Vector{Tuple{Float64,Float64,Float64,Int64}} = Tuple{Float64,Float64,Float64,Int64}[]
    _bivariate_quadratic_leq_constraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}} = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}[]
    _bivariate_quadratic_geq_constraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}} = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}[]
    _bivariate_quadratic_eq_constraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}} = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}[]

    _global_lower_bound::Float64 = -Inf
    _global_upper_bound::Float64 = Inf
    _maximum_node_id::Int64 = 0
    _iteration_count::Int64 = 0
    _node_count::Int64 = 0

    # Storage for output
    _solution_value::Float64 = 0.0
    _feasible_solution_found::Bool = false
    _first_solution_node::Int64 = -1
    _optimization_sense::MOI.OptimizationSense = MOI.MIN_SENSE
    _objective_value::Float64 = -Inf
    _termination_status_code::MOI.TerminationStatusCode = MOI.OTHER_RESULT_STATUS
    _result_status_code::MOI.ResultStatusCode = MOI.OPTIMIZE_NOT_CALLED

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

    _nlp_data::MOI.NLPBlockData = empty_nlp_data()

    _relaxed_evaluator::Evaluator = Evaluator{1,NS}()
    _relaxed_constraint_bounds::Vector{MOI.NLPBoundsPair} = Vector{MOI.NLPBoundsPair}[]
    _relaxed_eval_has_objective::Bool = false

    function Optimizer{S, T}(;options...) where {S <: MOI.AbstractOptimizer,
                                                 T <: MOI.AbstractOptimizer}

        m = new()

        # checks keywords supplied are field names and throws error otherwise
        allowed_kwargs = fieldnames(Optimizer{S,T})
        disallowed_kwargs = setdiff(collect(keys(options)), allowed_kwargs)
        if ~isempty(disallowed_kwargs)
            error("The following keyword arguments are not recognized by the
                   EAGO optimizer: $(disallowed_kwargs). Please consult
                   the documentation for allowed arguments.")
        end

        for i in keys(options)
            haskey(options,i) && setfield!(m, i, options[i])
        end

        m._time_left = m.time_limit
        m._user_branch_variables = ~isempty(m.branch_variable)

        return m
    end
end
function Optimizer(;options...)
    rtype = haskey(options, :relaxed_optimizer) ? typeof(options[:relaxed_optimizer]) : GLPK.Optimizer
    ropts = haskey(options, :relaxed_optimizer_kwargs) ? haskey(options, :relaxed_optimizer_kwargs) : Base.Iterators.Pairs(NamedTuple(),())
    utype = haskey(options, :upper_optimizer) ? typeof(options[:upper_optimizer]) : Ipopt.Optimizer

    opt = Optimizer{rtype, utype}(;options...)
    relax_fact = with_optimizer(rtype; ropts...)
    opt.relaxed_optimizer = relax_fact()
    if MOI.supports(opt.relaxed_optimizer, MOI.Silent())
        MOI.set(opt.relaxed_optimizer, MOI.Silent(), true)
    end
    return opt
end

function MOI.empty!(m::Optimizer)
    m = Optimizer()
end

function MOI.is_empty(m::Optimizer)

    flag = true
    flag &= isempty(m._variable_info)
    flag &= m._optimization_sense == MOI.MIN_SENSE
    flag &= m._termination_status_code == MOI.OPTIMIZE_NOT_CALLED

    return flag
end

function check_inbounds!(m::Optimizer, vi::MOI.VariableIndex)
    num_variables = length(m._variable_info)
    if !(1 <= vi.value <= num_variables)
        error("Invalid variable index $vi. ($num_variables variables in the model.)")
    end
    return
end

check_inbounds!(m::Optimizer, var::MOI.SingleVariable) = check_inbounds!(m, var.variable)

function check_inbounds!(m::Optimizer, aff::MOI.ScalarAffineFunction)
    for term in aff.terms
        check_inbounds!(m, term.variable_index)
    end
    return
end

function check_inbounds!(m::Optimizer, quad::MOI.ScalarQuadraticFunction)
    for term in quad.affine_terms
        check_inbounds!(m, term.variable_index)
    end
    for term in quad.quadratic_terms
        check_inbounds!(m, term.variable_index_1)
        check_inbounds!(m, term.variable_index_2)
    end
    return
end

function has_upper_bound(m::Optimizer, vi::MOI.VariableIndex)
    @inbounds val = m._variable_info[vi.value]
    return val.has_upper_bound
end

function has_lower_bound(m::Optimizer, vi::MOI.VariableIndex)
    @inbounds val = m._variable_info[vi.value]
    return val.has_lower_bound
end

function is_fixed(m::Optimizer, vi::MOI.VariableIndex)
    @inbounds val = m._variable_info[vi.value]
    return val.is_fixed
end

function is_integer_feasible(m::Optimizer)
    flag = true
    for var in m._integer_variables
        @inbounds val = m._lower_solution[var]
        if (0.0 < val < 1.0)
            flag = false
            break
        end
    end
    return flag
end

is_integer_variable(m::Optimizer, i::Int64) = m._variable_info[i].is_integer

function ReverseDict(dict)
    Dict(value => key for (key, value) in dict)
end

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; copy_names = false)
    return MOI.Utilities.default_copy_to(model, src, copy_names)
end

function label_nonlinear_variables!(m::Optimizer)
    _nlpdata = m._nlp_data
    x = _nlpdata.evaluator

    # scans subexpressions, objective, and constraints for nonlinear terms
    if ~isa(x, EmptyNLPEvaluator)
        if x.has_nlobj
            if (x.objective.linearity != JuMP._Derivatives.LINEAR) &&
               (x.objective.linearity != JuMP._Derivatives.CONSTANT)
                for i in 1:length(x.objective.nd)
                    nd = x.objective.nd[i]
                    if (nd.nodetype == JuMP._Derivatives.VARIABLE)
                        m.branch_variable[nd.index] = true
                        m.obbt_variable_values[nd.index] = true
                    end
                end
            end
        end
        for i in 1:length(x.constraints)
            if (x.constraints[i].linearity != JuMP._Derivatives.LINEAR) &&
               (x.constraints[i].linearity != JuMP._Derivatives.CONSTANT)
                for j in 1:length(x.constraints[i].nd)
                    nd = x.constraints[i].nd[j]
                    bool1 = (nd.nodetype == JuMP._Derivatives.VARIABLE)
                    if (nd.nodetype == JuMP._Derivatives.VARIABLE)
                        m.branch_variable[nd.index] = true
                        m.obbt_variable_values[nd.index] = true
                    end
                end
            end
        end
        for i in 1:length(x.subexpressions)
            if (x.subexpressions[i].linearity != JuMP._Derivatives.LINEAR) &&
               (x.subexpressions[i].linearity != JuMP._Derivatives.CONSTANT)
                for j in 1:length(x.subexpressions[i].nd)
                    nd = x.subexpressions[i].nd[j]
                    bool1 = (nd.nodetype == JuMP._Derivatives.VARIABLE)
                    if (nd.nodetype == JuMP._Derivatives.VARIABLE)
                        m.branch_variable[nd.index] = true
                        m.obbt_variable_values[nd.index] = true
                    end
                end
            end
        end
    end
    return
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

MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices) = [MOI.VariableIndex(i) for i in 1:length(m._variable_info)]
function MOI.get(m::Optimizer, ::MOI.ObjectiveValue)
    mult = 1.0
    if m._optimization_sense === MOI.MAX_SENSE
        mult *= -1.0
    end
    return mult*m._objective_value
end
MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = m._variable_number
function MOI.get(m::Optimizer, ::MOI.ObjectiveBound)
    if m._optimization_sense === MOI.MAX_SENSE
        bound = -m._global_lower_bound
    else
        bound = m._global_upper_bound
    end
    return bound
end
function MOI.get(m::Optimizer, ::MOI.RelativeGap)
    LBD = m._global_lower_bound
    UBD = m._global_upper_bound
    if m._optimization_sense === MOI.MAX_SENSE
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
function MOI.get(m::Optimizer, ::MOI.ResultCount)
    (m._result_status_code === MOI.FEASIBLE_POINT) ? 1 : 0
end

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    check_inbounds!(model, vi)
    return model._continuous_solution[vi.value]
end

const SCAL_FUNCS = Union{SV, SAF, SQF}
const INEQ_SETS = Union{LT, GT, ET}
MOI.supports_constraint(::Optimizer, ::Type{F}, ::Type{S}) where {F <: SCAL_FUNCS, S <: INEQ_SETS} = true

#=
const CONE_SETS = Union{NormInfinityCone, NormOneCone, SecondOrderCone, RotatedSecondOrderCone,
                        GeometricMeanCone, ExponentialCone, DualExponentialCone, PowerCone,
                        DualPowerCone, RelativeEntropyCone, NormSpectralCone, NormNuclearCone}
MOI.supports_constraint(::Optimizer, ::Type{VecOfVar}, ::Type{S}) where {S <: CONE_SETS} = true
=#

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.add_variable(m::Optimizer)
    m._variable_number += 1
    if ~m._user_branch_variables
        push!(m.branch_variable, false)
    end
    push!(m.obbt_variable_values, false)
    push!(m._obbt_working_lower_index, false)
    push!(m._obbt_working_upper_index, false)
    push!(m._lower_indx_diff, false)
    push!(m._upper_indx_diff, false)
    push!(m._old_low_index, false)
    push!(m._old_upp_index, false)
    push!(m._new_low_index, false)
    push!(m._new_upp_index, false)
    push!(m._fixed_variable, false)
    push!(m._variable_info, VariableInfo())
    return VI(m._variable_number)
end
MOI.add_variables(m::Optimizer, n::Int) = [MOI.add_variable(m) for i in 1:n]

#=
function MOI.add_constraint(m::Optimizer, v::SV, zo::MOI.ZO)
    vi = v.variable
    check_inbounds!(m, vi)
    has_upper_bound(m, vi) && error("Upper bound on variable $vi already exists.")
    has_lower_bound(m, vi) && error("Lower bound on variable $vi already exists.")
    is_fixed(m, vi) && error("Variable $vi is fixed. Cannot also set upper bound.")
    m._variable_info[vi.value].lower_bound = 0.0
    m._variable_info[vi.value].upper_bound = 1.0
    m._variable_info[vi.value].has_lower_bound = true
    m._variable_info[vi.value].has_upper_bound = true
    m._variable_info[vi.value].is_integer = true
    return CI{SV, MOI.ZO}(vi.value)
end
=#

function MOI.add_constraint(m::Optimizer, v::SV, lt::LT)
    vi = v.variable
    check_inbounds!(m, vi)
    if isnan(lt.upper)
        error("Invalid upper bound value $(lt.upper).")
    end
    if has_upper_bound(m, vi)
        error("Upper bound on variable $vi already exists.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is fixed. Cannot also set upper bound.")
    end
    m._variable_info[vi.value].upper_bound = lt.upper
    m._variable_info[vi.value].has_upper_bound = true
    return CI{SV, LT}(vi.value)
end

function MOI.add_constraint(m::Optimizer, v::SV, gt::GT)
    vi = v.variable
    check_inbounds!(m, vi)
    if isnan(gt.lower)
        error("Invalid lower bound value $(gt.lower).")
    end
    if has_lower_bound(m, vi)
        error("Lower bound on variable $vi already exists.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is fixed. Cannot also set lower bound.")
    end
    m._variable_info[vi.value].lower_bound = gt.lower
    m._variable_info[vi.value].has_lower_bound = true
    return CI{SV, GT}(vi.value)
end

function MOI.add_constraint(m::Optimizer, v::SV, eq::ET)
    vi = v.variable
    check_inbounds!(m, vi)
    if isnan(eq.value)
        error("Invalid fixed value $(gt.lower).")
    end
    if has_lower_bound(m, vi)
        error("Variable $vi has a lower bound. Cannot be fixed.")
    end
    if has_upper_bound(m, vi)
        error("Variable $vi has an upper bound. Cannot be fixed.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is already fixed.")
    end
    m._variable_info[vi.value].lower_bound = eq.value
    m._variable_info[vi.value].upper_bound = eq.value
    m._variable_info[vi.value].has_lower_bound = true
    m._variable_info[vi.value].has_upper_bound = true
    m._variable_info[vi.value].is_fixed = true
    return CI{SV, ET}(vi.value)
end

#=
function MOI.add_constraint(m::Optimizer, v::SV, eq::MOI.IT)
    vi = v.variable
    check_inbounds!(m, vi)
    if isnan(eq.lower)
        error("Invalid fixed value $(gt.lower).")
    end
    if isnan(eq.upper)
        error("Invalid fixed value $(gt.upper).")
    end
    if has_lower_bound(m, vi)
        error("Lower bound on variable $vi already exists. Cannot also set interval bounds.")
    end
    if has_upper_bound(m, vi)
        error("Upper bound on variable $vi already exists. Cannot also set interval bounds.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is fixed. Cannot also set interval bounds.")
    end
    m._variable_info[vi.value].lower_bound = eq.lower
    m._variable_info[vi.value].upper_bound = eq.upper
    m._variable_info[vi.value].has_lower_bound = true
    m._variable_info[vi.value].has_upper_bound = true
    return CI{SV, MOI.IT}(vi.value)
end
=#

macro define_addconstraint_linear(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            push!(m.$(array_name), (func, set, length(func.terms)))
            indx = CI{$function_type, $set_type}(length(m.$(array_name)))
            return indx
        end
    end
end

@define_addconstraint_linear SAF LT _linear_leq_constraints
@define_addconstraint_linear SAF GT _linear_geq_constraints
@define_addconstraint_linear SAF ET _linear_eq_constraints

macro define_addconstraint_quadratic(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            for i in func.affine_terms m.branch_variable[i.variable_index.value] = true end
            for i in func.quadratic_terms
                m.branch_variable[i.variable_index_1.value] = true
                m.branch_variable[i.variable_index_2.value] = true
            end
            push!(m.$(array_name), (func, set, length(m.$(array_name))+1))
            indx = CI{$function_type, $set_type}(length(m.$(array_name)))
            return indx
        end
    end
end

@define_addconstraint_quadratic SQF LT _quadratic_leq_constraints
@define_addconstraint_quadratic SQF GT _quadratic_geq_constraints
@define_addconstraint_quadratic SQF ET _quadratic_eq_constraints

#=
macro define_addconstraint_cone(set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::VecOfVar, set::$set_type)
            if length(func.variables) != dimension(s)
                error("Dimension of $(s) does not match number of terms in $(f)")
            end
            check_inbounds!(m, func)
            push!(m.$(array_name), (func, set)))
            return CI{$function_type, $set_type}(Number)
        end
    end
end

@define_addconstraint_cone NormInfinityCone _conic_norm_infinity
@define_addconstraint_cone NormOneCone _conic_norm_one
@define_addconstraint_cone SecondOrderCone _conic_second_order
@define_addconstraint_cone RotatedSecondOrderCone _conic_rotated_second_order
@define_addconstraint_cone GeometricMeanCone _conic_geometric_mean
@define_addconstraint_cone ExponentialCone _conic_exponential
@define_addconstraint_cone DualExponentialCone _conic_dual_exponential
@define_addconstraint_cone PowerCone _conic_power_cone
@define_addconstraint_cone DualPowerCone _conic_dual_power
@define_addconstraint_cone RelativeEntropyCone _conic_relative_entropy
@define_addconstraint_cone NormSpectralCone _conic_norm_spectral
@define_addconstraint_cone NormNuclearCone _conic_norm_nuclear
=#

function MOI.set(m::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    m._nlp_data = nlp_data
    return
end

MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{SV}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{SAF}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{SQF}) = true

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction, func::SV)
    check_inbounds!(m, func)
    m._objective_sv = func
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction, func::SAF)
    check_inbounds!(m, func)
    m._objective_saf = func
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction, func::SQF)
    check_inbounds!(m, func)
    m._objective_sqf = func
    for term in func.quadratic_terms
        @inbounds m.branch_variable[term.variable_index_1.value] = true
        @inbounds m.branch_variable[term.variable_index_1.value] = true
    end
    #for term in func.affine_terms
    #    m.branch_variable[term.variable_index] = true
    #end
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    m._optimization_sense = sense
    return
end

# Defines single variable objective function
eval_function(var::SV, x) = x[var.variable.value]

# Defines affine objective function
function eval_function(aff::SAF, x)
    function_value = aff.constant
    for term in aff.terms
        @inbounds function_value += term.coefficient*x[term.variable_index.value]
    end
    return function_value
end

# Defines quadratic objective function
function eval_function(quad::SQF, x)
    function_value = quad.constant
    for term in quad.affine_terms
        @inbounds function_value += term.coefficient*x[term.variable_index.value]
    end
    for term in quad.quadratic_terms
        row_idx = term.variable_index_1
        col_idx = term.variable_index_2
        coefficient = term.coefficient
        function_value += 0.5*coefficient*x[row_idx.value]*x[col_idx.value]
    end
    return function_value
end

# Defines evaluation function for objective
function eval_objective(m::Optimizer, x)
    @assert !(m._nlp_data.has_objective && m._objective !== nothing)
    if m._nlp_data.has_objective
        return MOI.eval_objective(m._nlp_data.evaluator, x)
    elseif m._objective !== nothing
        return eval_function(m._objective, x)
    else
        return 0.0
    end
end


"""
    log_iteration!(x::Optimizer)

If 'logging_on' is true, the 'global_lower_bound', 'global_upper_bound',
'run_time', and 'node_count' are stored every 'log_interval'. If
'log_subproblem_info' then the lower bound, feasibility and run times of the
subproblems are logged every 'log_interval'.
"""
function log_iteration!(x::Optimizer)

    if x.log_on

        log = x._log
        if (mod(x._iteration_count, x.log_interval) == 0 || x._iteration_count == 1)

            if x.log_subproblem_info
                if x._optimization_sense === MOI.MIN_SENSE
                    push!(log.current_lower_bound, x._lower_objective_value)
                    push!(log.current_upper_bound, x._upper_objective_value)
                else
                    push!(log.current_lower_bound, -x._upper_objective_value)
                    push!(log.current_upper_bound, -x._lower_objective_value)
                end

                push!(log.preprocessing_time, x._last_preprocess_time)
                push!(log.lower_problem_time, x._last_lower_problem_time)
                push!(log.upper_problem_time, x._last_upper_problem_time)
                push!(log.postprocessing_time, x._last_postprocessing_time)

                push!(log.preprocess_feasibility, x._preprocess_feasibility)
                push!(log.lower_problem_feasibility, x._lower_feasibility)
                push!(log.upper_problem_feasibility, x._upper_feasibility)
                push!(log.postprocess_feasibility, x._postprocess_feasibility)
            end

            if x._optimization_sense === MOI.MIN_SENSE
                push!(log.global_lower_bound, x._global_lower_bound)
                push!(log.global_upper_bound, x._global_upper_bound)
            else
                push!(log.global_lower_bound, -x._global_upper_bound)
                push!(log.global_upper_bound, -x._global_lower_bound)
            end
            push!(log.run_time, x._run_time)
            push!(log.node_count, x._node_count)

        end
    end
    return
end
