"""
    VariableInfo

A structure used to store information related to the bounds assigned to each
variable.
- `is_integer::Bool`:      Is the variable integer valued?
- `lower_bound::Float64`:  May be -Inf even if has_lower_bound == true
- `has_lower_bound::Bool`: Implies lower_bound == Inf
- `upper_bound::Float64`:  May be Inf even if has_upper_bound == true
- `has_upper_bound::Bool`: Implies upper_bound == Inf
- `is_fixed::Bool`:        Implies lower_bound == upper_bound and
                           !has_lower_bound and !has_upper_bound.
"""
mutable struct VariableInfo
    is_integer::Bool
    lower_bound::Float64           # May be -Inf even if has_lower_bound == true
    has_lower_bound::Bool          # Implies lower_bound == Inf
    upper_bound::Float64           # May be Inf even if has_upper_bound == true
    has_upper_bound::Bool          # Implies upper_bound == Inf
    is_fixed::Bool                 # Implies lower_bound == upper_bound and
                                   # !has_lower_bound and !has_upper_bound.
end
VariableInfo() = VariableInfo(false,-Inf, false, Inf, false, false)
lower_bound(x::VariableInfo) = x.lower_bound
upper_bound(x::VariableInfo) = x.upper_bound

"""
    ExtensionType

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

export Optimizer
"""
    Optimizer

The main optimizer object used by EAGO to solve problems during the optimization
routine. The following commonly used options are described below and can be set
via keyword arguments in the JuMP/MOI model:
- `presolve_scrubber_flag::Bool`: Replace code in user-defined functions which
                                  may prevent method overloading on Real subtypes
                                  (default = false).
- `presolve_to_JuMP_flag::Bool`: Create and use DAG representations of user-defined
                                 function (default = false).
- `presolve_epigraph_flag::Bool`: [FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED]
                                  Apply the epigraph reformulation to the problem
                                  (default = false).
- `presolve_cse_flag::Bool`: [FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Enable
                             common subexpression elimination for DAG (default = false).
- `presolve_flatten_flag::Bool`: Rerranges the DAG using registered transformations
                                 (default = false)
- `cp_depth::Int64`: Depth in B&B tree above which constraint propagation should
                     be disabled (default = 1000).
- `cp_repetitions::Int64`: Number of repetitions of forward-reverse passes to perform in
                          constraint propagation (default = 3).
- `cp_tolerance::Float64`: Disable constraint propagation if the ratio of new node
                           volume to beginning node volume exceeds this number
                           (default = 0.99).
- `cp_interval_only::Bool`: Use only valid interval bounds during constraint
                            propagation (default = false).
- `relaxed_optimizer::S`: An instance of the optimizer used to solve the relaxed
                          subproblems (default = GLPK.Optimizer()).
- `relaxed_optimizer_kwargs::Base.Iterators.Pairs`: Keyword arguments for the
                                                    relaxed optimizer.
- `obbt_depth::Int64`: Depth in B&B tree above which OBBT should
                     be disabled (default = 1000).
- `obbt_repetitions::Int64`: Number of repetitions of OBBT to perform in
                            preprocessing (default = 3).
- `obbt_aggressive_on::Bool`: Turn aggresive OBBT on (default = false).
- `obbt_aggressive_max_iteration::Int64`: Maximum iteration to perform aggresive
                                          OBBT (default = 2)
- `obbt_aggressive_min_dimension::Int64`: Minimum dimension to perform aggresive
                                          OBBT (default = 2)
- `obbt_tolerance::Float64`: Tolerance to consider bounds equal (default = 1E-9).
- `obbt_variable_values::Vector{Bool}`: Variables to perform OBBT on
                                        (default: all variables in nonlinear expressions).
- `lp_depth::Int64`: Depth in B&B tree above which linear FBBT should
                     be disabled (default = 1000).
- `lp_repetitions::Int64`: Number of repetitions of linear FBBT to perform in
                            preprocessing (default = 3).
- `quad_uni_depth::Int64`: Depth in B&B tree above which univariate quadratic
                           FBBT should be disabled (default = -1).
- `quad_uni_repetitions::Int64`: Number of repetitions of univariate quadratic FBBT
                                 to perform in preprocessing (default = 2).
- `quad_bi_depth::Int64`: [FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Depth in
                          B&B tree above which bivariate quadratic FBBT should
                          be disabled (default = -1).
- `quad_bi_repetitions::Int64`: Number of repetitions of bivariate quadratic FBBT
                                 to perform in preprocessing (default = 2).
- `subgrad_tighten::Bool`: Perform tightening of interval bounds using subgradients
                           at each factor in each nonlinear tape during a forward-reverse
                           pass (default = true).
- `subgrad_tighten_reverse::Bool`: [FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Used to
                                   enable/disable subgradient tightening of interval
                                   bounds on the reverse pass (default = true).
- `cut_max_iterations::Int64`
- `cut_cvx::Float64`: Convex coefficient used to select point for new added cuts.
                      Branch point is given by `(1-cut_cvx)*xmid + cut_cvx*xsol`
                      (default = 0.9).
- `cut_tolerance::Float64`: Add cut if the L1 distance from the prior cutting point
                            to the new cutting point normalized by the box volume
                            is greater than the tolerance (default = 0.05).
- `objective_cut_on::Bool`: Adds an objective cut to the relaxed problem (default = true).
- `upper_optimizer::T`: Optimizer used to solve upper bounding problem (default = Ipopt.Optimizer())
- `upper_factory::JuMP.OptimizerFactory`: OptimizerFactory used to build optimizer that
                                          solves the upper bounding problem
                                          (default = with_optimizer(Ipopt.Optimizer, kwargs),
                                          check Optimizer constructor for kwargs used).
- `upper_bounding_depth::Int64`: Solve upper problem for every node with depth
                                 less than `upper_bounding_depth` and with a probability
                                 of (1/2)^(depth-upper_bounding_depth) otherwise
                                 (default = 4).
- `dbbt_depth::Int64`: Depth in B&B tree above which duality-based bound tightening
                       should be disabled (default = 1E10).
- `dbbt_tolerance::Float64`: New bound is considered equal to the prior bound
                             if within dbbt_tolerance (default = 1E-9).
- `branch_cvx_factor::Float64`: Convex coefficient used to select branch point.
                                Branch point is given by
                                `branch_cvx_factor*xmid + (1-branch_cvx_factor)*xsol`
                                (default = 0.25)
- `branch_offset::Float64`: Minimum distance from bound to have branch point
                            normalized by width of dimension to branch on
                            (default = 0.15)
- `branch_variable::Vector{Bool}`: Variables to branch on (default is all nonlinear).
- `branch_max_repetitions::Int64`: [FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED]
                                   Number of times repeat node processing prior
                                   to branching (default = 4).
- `branch_repetition_tol::Float64`: [FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED]
                                    Volume ratio tolerance required to repeat
                                    processing the current node (default = 0.9).
- `rounding_mode::Symbol`: Interval rounding mode to use (default = :accurate)
- `node_limit::Int64`: Node limit  (default = 10^7).
- `time_limit::Float64`: Time limit in seconds (default = 3600).
- `iteration_limit::Int64`: Iteration limit (default = 3000000).
- `absolute_tolerance::Float64`: Absolute tolerance for terminatin (default = 1E-3).
- `relative_tolerance::Float64`: Relative tolerance for terminatin (default = 1E-3).
- `local_solve_only::Bool`: Perform only a local solve of the problem (default = false).
- `log_on::Bool`: Turns logging on records global bounds, node count and run time.
                  Additional options are available for recording information specific
                  to subproblems (default = false).
- `log_subproblem_info::Bool`: Turns on logging of times and feasibility of
                               subproblems (default = false).
- `log_interval::Bool`: Log data every `log_interval` iterations (default = 1).
- `verbosity::Int64`: The amount of information that should be printed to console
                      while solving values range from 0 - 4: 0 is silent, 1 shows
                      iteration summary statistics only, 2-4 show varying degrees
                      of details about calculations within each iteration
                      (default = 1).
- `output_iterations::Int64`: Display summary of iteration to console every
                             `output_iterations` (default = 10).
- `header_iterations::Int64`: Display header for summary to console every
                             `output_iterations` (default = 100).
- `enable_optimize_hook::Bool`: Specifies that the optimize_hook! function should
                                be called rather than throw the problem to the
                                standard B&B routine (default = false).
- `ext::Dict{Symbol, Any}`: Holds additional storage needed for constructing
                            extensions to EAGO (default = Dict{Symbol,Any}).
- `ext_type::ExtensionType`: Holds an instance of a subtype of `EAGO.ExtensionType`
                             used to define new custom subroutines
                             (default = DefaultExt()).
"""
mutable struct Optimizer{S<:MOI.AbstractOptimizer, T<:MOI.AbstractOptimizer} <: MOI.AbstractOptimizer

    # Presolving options
    presolve_scrubber_flag::Bool
    presolve_to_JuMP_flag::Bool
    presolve_epigraph_flag::Bool
    presolve_cse_flag::Bool
    presolve_flatten_flag::Bool

    # Options for constraint propagation
    cp_depth::Int64
    cp_improvement::Float64
    cp_repetitions::Int64
    cp_tolerance::Float64
    cp_interval_only::Bool

    # Options for optimality-based bound tightening
    relaxed_optimizer::S
    relaxed_optimizer_kwargs::Base.Iterators.Pairs
    relaxed_inplace_mod::Bool
    obbt_depth::Int64
    obbt_repetitions::Int64
    obbt_aggressive_on::Bool
    obbt_aggressive_max_iteration::Int64
    obbt_aggressive_min_dimension::Int64
    obbt_tolerance::Float64
    obbt_variable_values::Vector{Bool}

    # Options for linear bound tightening
    lp_depth::Int64
    lp_repetitions::Int64

    # Options for quadratic bound tightening
    quad_uni_depth::Int64
    quad_uni_repetitions::Int64
    quad_bi_depth::Int64
    quad_bi_repetitions::Int64

    # Subgradient tightening flag
    subgrad_tighten::Bool
    subgrad_tighten_reverse::Bool

    # Tolerance to add cuts and max number of cuts
    cut_max_iterations::Int64
    cut_cvx::Float64
    cut_tolerance::Float64
    objective_cut_on::Bool

    # Upper bounding options
    upper_optimizer::T
    upper_factory::JuMP.OptimizerFactory
    upper_bounding_depth::Int64

    # Duality-based bound tightening (DBBT) options
    dbbt_depth::Int64
    dbbt_tolerance::Float64

    # Node branching options
    branch_cvx_factor::Float64
    branch_offset::Float64
    branch_variable::Vector{Bool}
    branch_max_repetitions::Int64
    branch_repetition_tol::Float64

    # Rounding mode used with interval arithmetic
    rounding_mode::Symbol

    # Termination limits
    node_limit::Int64
    time_limit::Float64
    iteration_limit::Int64
    absolute_tolerance::Float64
    relative_tolerance::Float64
    local_solve_only::Bool
    feasible_local_continue::Bool

    # Iteration logging options
    log_on::Bool
    log_subproblem_info::Bool
    log_interval::Bool

    # Optimizer display options
    verbosity::Int64
    output_iterations::Int64
    header_iterations::Int64

    # Debug
    enable_optimize_hook::Bool
    ext::Dict{Symbol, Any}
    ext_type::ExtensionType

    _current_node::NodeBB
    _current_xref::Vector{Float64}
    _sense_multiplier::Float64

    _variable_number::Int64
    _state_variables::Int64
    _continuous_variable_number::Int64
    _integer_variable_number::Int64

    _user_branch_variables::Bool
    _fixed_variable::Vector{Bool}
    _constraint_convexity::Dict{CI, Bool}

    _continuous_solution::Vector{Float64}

    _integer_variables::Vector{Int64}
    _variable_info::Vector{VariableInfo}
    _upper_variables::Vector{VI}

    _stack::BinaryMinMaxHeap{NodeBB}

    _lower_variable::Vector{SV}
    _lower_variable_index::Vector{VI}
    _lower_variable_values::Vector{Int64}
    _lower_variable_et::Vector{CI{SV, ET}}
    _lower_variable_lt::Vector{CI{SV, LT}}
    _lower_variable_gt::Vector{CI{SV, GT}}
    _lower_variable_et_indx::Vector{Int64}
    _lower_variable_lt_indx::Vector{Int64}
    _lower_variable_gt_indx::Vector{Int64}

    _preprocess_feasibility::Bool
    _preprocess_result_status::MOI.ResultStatusCode
    _preprocess_termination_status::MOI.TerminationStatusCode

    _lower_result_status::MOI.ResultStatusCode
    _lower_termination_status::MOI.TerminationStatusCode
    _lower_feasibility::Bool
    _lower_objective_value::Float64
    _lower_solution::Vector{Float64}
    _lower_lvd::Vector{Float64}
    _lower_uvd::Vector{Float64}

    _cut_result_status::MOI.ResultStatusCode
    _cut_termination_status::MOI.TerminationStatusCode
    _cut_solution::Vector{Float64}
    _cut_objective_value::Float64
    _cut_feasibility::Bool

    _upper_result_status::MOI.ResultStatusCode
    _upper_termination_status::MOI.TerminationStatusCode
    _upper_feasibility::Bool
    _upper_objective_value::Float64
    _upper_solution::Vector{Float64}

    _best_upper_value::Float64

    _postprocess_feasibility::Bool

    _start_time::Float64
    _run_time::Float64
    _time_left::Float64
    _parse_time::Float64
    _presolve_time::Float64
    _last_preprocess_time::Float64
    _last_lower_problem_time::Float64
    _last_upper_problem_time::Float64
    _last_postprocessing_time::Float64

    _objective_sv::SV
    _objective_saf::SAF
    _objective_sqf::SQF
    _objective_is_sv::Bool
    _objective_is_saf::Bool
    _objective_is_sqf::Bool
    _objective_is_nlp::Bool

    _objective_convexity::Bool

    _objective_cut_set::Int64
    _objective_cut_ci_sv::CI{SV,LT}
    _objective_cut_ci_saf::Vector{CI{SAF,LT}}

    _linear_leq_constraints::Vector{Tuple{SAF, LT, Int64}}
    _linear_geq_constraints::Vector{Tuple{SAF, GT, Int64}}
    _linear_eq_constraints::Vector{Tuple{SAF, ET, Int64}}

    _quadratic_leq_constraints::Vector{Tuple{SQF, LT, Int64}}
    _quadratic_geq_constraints::Vector{Tuple{SQF, GT, Int64}}
    _quadratic_eq_constraints::Vector{Tuple{SQF, ET, Int64}}

    _quadratic_leq_dict::Vector{ImmutableDict{Int64,Int64}}
    _quadratic_geq_dict::Vector{ImmutableDict{Int64,Int64}}
    _quadratic_eq_dict::Vector{ImmutableDict{Int64,Int64}}

    _quadratic_ci_leq::Vector{Vector{CI{SAF,LT}}}
    _quadratic_ci_geq::Vector{Vector{CI{SAF,LT}}}
    _quadratic_ci_eq::Vector{Vector{Tuple{CI{SAF,LT},CI{SAF,LT}}}}

    _quadratic_leq_sparsity::Vector{Vector{VI}}
    _quadratic_geq_sparsity::Vector{Vector{VI}}
    _quadratic_eq_sparsity::Vector{Vector{VI}}

    _quadratic_leq_gradnz::Vector{Int64}
    _quadratic_geq_gradnz::Vector{Int64}
    _quadratic_eq_gradnz::Vector{Int64}

    _quadratic_leq_convexity::Vector{Bool}
    _quadratic_geq_convexity::Vector{Bool}
    _quadratic_eq_convexity_1::Vector{Bool}
    _quadratic_eq_convexity_2::Vector{Bool}

    _lower_nlp_affine::Vector{Vector{CI{SAF,LT}}}
    _upper_nlp_affine::Vector{Vector{CI{SAF,LT}}}

    _lower_nlp_affine_indx::Vector{Int64}
    _upper_nlp_affine_indx::Vector{Int64}

    _lower_nlp_sparsity::Vector{Vector{Int64}}
    _upper_nlp_sparsity::Vector{Vector{Int64}}

    _univariate_quadratic_leq_constraints::Vector{Tuple{Float64,Float64,Float64,Int64}}
    _univariate_quadratic_geq_constraints::Vector{Tuple{Float64,Float64,Float64,Int64}}
    _univariate_quadratic_eq_constraints::Vector{Tuple{Float64,Float64,Float64,Int64}}
    _bivariate_quadratic_leq_constraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}}
    _bivariate_quadratic_geq_constraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}}
    _bivariate_quadratic_eq_constraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}}

    _global_lower_bound::Float64
    _global_upper_bound::Float64
    _maximum_node_id::Int64
    _iteration_count::Int64
    _node_count::Int64

    # Storage for output
    _solution_value::Float64
    _feasible_solution_found::Bool
    _first_solution_node::Int64
    _optimization_sense::MOI.OptimizationSense
    _objective_value::Float64
    _termination_status_code::MOI.TerminationStatusCode
    _result_status_code::MOI.ResultStatusCode

    # Optimality-Based Bound Tightening (OBBT) Options
    _obbt_working_lower_index::Vector{Bool}
    _obbt_working_upper_index::Vector{Bool}
    _lower_indx_diff::Vector{Bool}
    _upper_indx_diff::Vector{Bool}
    _old_low_index::Vector{Bool}
    _old_upp_index::Vector{Bool}
    _new_low_index::Vector{Bool}
    _new_upp_index::Vector{Bool}
    _obbt_variables::Vector{VI}
    _obbt_performed_flag::Bool

    # Feasibility-Based Bound Tightening Options
    _cp_improvement::Float64
    _cp_evaluation_reverse::Bool

    _cut_iterations::Int64
    _cut_add_flag::Bool

    # Options for Repetition (If DBBT Performed Well)
    _node_repetitions::Int64
    _initial_volume::Float64
    _final_volume::Float64

    # Log
    _log::Log

    _nlp_data::MOI.NLPBlockData

    _relaxed_evaluator::Evaluator
    _relaxed_constraint_bounds::Vector{MOI.NLPBoundsPair}
    _relaxed_eval_has_objective::Bool

    function Optimizer{S, T}(;options...) where {S <: MOI.AbstractOptimizer,
                                                 T <: MOI.AbstractOptimizer}

        m = new()

        # checks that all keywords supplied to the optimizers are field names
        # throws error otherwise
        allowed_kwargs = fieldnames(Optimizer{S,T})
        disallowed_kwargs = setdiff(collect(keys(options)), allowed_kwargs)
        if ~isempty(disallowed_kwargs)
            error("The following keyword arguments are not recognized by the
                   EAGO optimizer: $(disallowed_kwargs). Please consult
                   the documentation for allowed arguments.")
        end

        default_opt_dict = Dict{Symbol,Any}()

        # Presolving options
        default_opt_dict[:presolve_scrubber_flag] = false
        default_opt_dict[:presolve_to_JuMP_flag] = false
        default_opt_dict[:presolve_epigraph_flag] = false
        default_opt_dict[:presolve_cse_flag] = false
        default_opt_dict[:presolve_flatten_flag] = false

        # Options for constraint propagation
        default_opt_dict[:cp_depth] = 1000
        default_opt_dict[:cp_repetitions] = 3
        default_opt_dict[:cp_tolerance] = 0.99
        default_opt_dict[:cp_interval_only] = false

        # Options for optimality-based bound tightening
        default_opt_dict[:obbt_depth] = 6
        default_opt_dict[:obbt_repetitions] = 20
        default_opt_dict[:obbt_aggressive_on] = false
        default_opt_dict[:obbt_aggressive_max_iteration] = 2
        default_opt_dict[:obbt_aggressive_min_dimension] = 2
        default_opt_dict[:obbt_tolerance] = 1E-9
        default_opt_dict[:obbt_variable_values] = Bool[]

        # Options for linear bound tightening
        default_opt_dict[:lp_depth] = 100000
        default_opt_dict[:lp_repetitions] = 3

        # Options for quadratic bound tightening
        default_opt_dict[:quad_uni_depth] = -1
        default_opt_dict[:quad_uni_repetitions] = 2
        default_opt_dict[:quad_bi_depth] = -1
        default_opt_dict[:quad_bi_repetitions] = 2

        # Subgradient tightening flags for evaluation
        default_opt_dict[:subgrad_tighten] = true
        default_opt_dict[:subgrad_tighten_reverse] = false

        # Tolerance to add cuts and max number of cuts
        default_opt_dict[:objective_cut_on] = true
        default_opt_dict[:cut_max_iterations] = 3
        default_opt_dict[:cut_cvx] = 0.9
        default_opt_dict[:cut_tolerance] = 0.05

        # Upper bounding options
        default_opt_dict[:upper_bounding_depth] = 4

        # Duality-based bound tightening (DBBT) options
        default_opt_dict[:dbbt_depth] = 10^10
        default_opt_dict[:dbbt_tolerance] = 1E-8

        # Node branching options
        default_opt_dict[:branch_cvx_factor] = 0.25
        default_opt_dict[:branch_offset] = 0.15
        default_opt_dict[:branch_variable] = Bool[]
        default_opt_dict[:branch_max_repetitions] = 4
        default_opt_dict[:branch_repetition_tol] = 0.9

        # Rounding mode used with interval arithmetic
        default_opt_dict[:rounding_mode] = :accurate

        # Termination limits
        default_opt_dict[:node_limit] = 10^7
        default_opt_dict[:time_limit] = 1000.0
        default_opt_dict[:iteration_limit] = 3000000
        default_opt_dict[:absolute_tolerance] = 1E-3
        default_opt_dict[:relative_tolerance] = 1E-3
        default_opt_dict[:local_solve_only] = false
        default_opt_dict[:feasible_local_continue] = false

        # Iteration logging options
        default_opt_dict[:log_on] = false
        default_opt_dict[:log_subproblem_info] = false
        default_opt_dict[:log_interval] = 1

        # Optimizer display options
        default_opt_dict[:verbosity] = 1
        default_opt_dict[:output_iterations] = 10
        default_opt_dict[:header_iterations] = 100

        # Extension options
        default_opt_dict[:enable_optimize_hook] = false
        default_opt_dict[:ext] = Dict{Symbol, Any}()
        default_opt_dict[:ext_type] = DefaultExt()

        default_opt_dict[:relaxed_optimizer] = GLPK.Optimizer()
        default_opt_dict[:relaxed_optimizer_kwargs] = Base.Iterators.Pairs(NamedTuple(),())
        default_opt_dict[:relaxed_inplace_mod] = true
        default_opt_dict[:upper_optimizer] = Ipopt.Optimizer()
        #=
        fac = with_optimizer(Ipopt.Optimizer, max_iter = 1000, acceptable_tol = 1E30,
                             acceptable_iter = 100, constr_viol_tol = 0.00001,
                             acceptable_constr_viol_tol = 1E-6, print_level = 0)
        =#
        fac = with_optimizer(Ipopt.Optimizer, print_level = 0,
                             acceptable_tol = 1E30,
                             max_iter = 1000000,
                             acceptable_iter = 50000,
                             constr_viol_tol = 0.0000001,
                             acceptable_constr_viol_tol = 0.0000001,
                             acceptable_dual_inf_tol = 1.0,
                             acceptable_compl_inf_tol = 0.0000001)
        default_opt_dict[:upper_factory] = fac

        for i in keys(default_opt_dict)
            if (haskey(options,i))
                setfield!(m, i, options[i])
            else
                setfield!(m, i, default_opt_dict[i])
            end
        end

        m._global_lower_bound = -Inf
        m._global_upper_bound = Inf
        m._maximum_node_id = 0
        m._iteration_count = 0
        m._node_count = 0

        m._stack = BinaryMinMaxHeap{NodeBB}()
        m._current_node = NodeBB()
        m._current_xref = Float64[]

        m._variable_info = VariableInfo[]
        m._variable_number = 0
        m._continuous_variable_number = 0

        m._user_branch_variables = ~isempty(m.branch_variable)
        m._fixed_variable = Bool[]
        m._constraint_convexity = Dict{CI,Bool}()

        m._continuous_solution = Float64[]

        m._integer_variables = Int64[]
        m._upper_variables = VI[]

        m._lower_variable = SV[]
        m._lower_variable_index = VI[]
        m._lower_variable_et = CI{SV, ET}[]
        m._lower_variable_lt = CI{SV, LT}[]
        m._lower_variable_gt = CI{SV, GT}[]
        m._lower_variable_et_indx = Int64[]
        m._lower_variable_lt_indx = Int64[]
        m._lower_variable_gt_indx = Int64[]

        m._preprocess_feasibility = true
        m._preprocess_result_status = MOI.OTHER_RESULT_STATUS
        m._preprocess_termination_status = MOI.OPTIMIZE_NOT_CALLED

        m._lower_result_status = MOI.OTHER_RESULT_STATUS
        m._lower_termination_status = MOI.OPTIMIZE_NOT_CALLED
        m._lower_feasibility = false
        m._lower_objective_value = -Inf
        m._lower_solution = Float64[]
        m._lower_lvd = Float64[]
        m._lower_uvd = Float64[]

        m._cut_result_status = MOI.OTHER_RESULT_STATUS
        m._cut_termination_status = MOI.OPTIMIZE_NOT_CALLED
        m._cut_solution = Float64[]
        m._cut_objective_value = -Inf
        m._cut_feasibility = false

        m._upper_result_status = MOI.OTHER_RESULT_STATUS
        m._upper_termination_status = MOI.OPTIMIZE_NOT_CALLED
        m._upper_feasibility = false
        m._upper_objective_value = Inf
        m._upper_solution = Float64[]

        m._best_upper_value = Inf

        m._postprocess_feasibility = false

        m._start_time = 0.0
        m._run_time = 0.0
        m._time_left = m.time_limit
        m._parse_time = 0.0
        m._presolve_time = 0.0
        m._last_preprocess_time = 0.0
        m._last_lower_problem_time = 0.0
        m._last_upper_problem_time = 0.0
        m._last_postprocessing_time = 0.0

        m._objective_sv = SV(VI(1))
        m._objective_saf = SAF(SAT.([0.0],[VI(1)]), 0.0)
        m._objective_sqf = SQF(SAT.([0.0],[VI(1)]), SQT.([0.0],[VI(1)],[VI(1)]), 0.0)
        m._objective_is_sv = false
        m._objective_is_saf = false
        m._objective_is_sqf = false
        m._objective_is_nlp = false
        m._objective_cut_set = -1
        m._objective_cut_ci_sv = CI{SV,LT}(-1.0)
        m._objective_cut_ci_saf = CI{SAF,LT}[]

        m._objective_convexity = false

        m._linear_leq_constraints = Tuple{SAF, LT, Int64}[]
        m._linear_geq_constraints = Tuple{SAF, GT, Int64}[]
        m._linear_eq_constraints = Tuple{SAF, ET, Int64}[]

        m._quadratic_leq_constraints = Tuple{SQF, LT, Int64}[]
        m._quadratic_geq_constraints = Tuple{SQF, GT, Int64}[]
        m._quadratic_eq_constraints = Tuple{SQF, ET, Int64}[]

        m._quadratic_leq_dict = ImmutableDict{Int64,Int64}[]
        m._quadratic_geq_dict = ImmutableDict{Int64,Int64}[]
        m._quadratic_eq_dict = ImmutableDict{Int64,Int64}[]

        m._quadratic_leq_sparsity = Vector{VI}[]
        m._quadratic_geq_sparsity = Vector{VI}[]
        m._quadratic_eq_sparsity = Vector{VI}[]

        m._quadratic_ci_leq = CI{SAF,LT}[]
        m._quadratic_ci_geq = CI{SAF,LT}[]
        m._quadratic_ci_eq = Tuple{CI{SAF,LT},CI{SAF,LT}}[]

        m._quadratic_leq_gradnz = Int64[]
        m._quadratic_geq_gradnz = Int64[]
        m._quadratic_eq_gradnz = Int64[]

        m._quadratic_leq_convexity = Bool[]
        m._quadratic_geq_convexity = Bool[]
        m._quadratic_eq_convexity_1 = Bool[]
        m._quadratic_eq_convexity_2 = Bool[]


        m._lower_nlp_affine = Vector{CI{SAF,LT}}[]
        m._upper_nlp_affine = Vector{CI{SAF,LT}}[]

        m._lower_nlp_affine_indx = Int64[]
        m._upper_nlp_affine_indx = Int64[]

        m._lower_nlp_sparsity = Vector{Int64}[]
        m._upper_nlp_sparsity = Vector{Int64}[]

        m._univariate_quadratic_leq_constraints = Tuple{Float64,Float64,Float64,Int64}[]
        m._univariate_quadratic_geq_constraints = Tuple{Float64,Float64,Float64,Int64}[]
        m._univariate_quadratic_eq_constraints = Tuple{Float64,Float64,Float64,Int64}[]
        m._bivariate_quadratic_leq_constraints = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}[]
        m._bivariate_quadratic_geq_constraints = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}[]
        m._bivariate_quadratic_eq_constraints = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Int64}[]

        # Storage for output
        m._solution_value = 0.0
        m._feasible_solution_found = false
        m._first_solution_node = -1
        m._objective_value = -Inf
        m._optimization_sense = MOI.MIN_SENSE
        m._result_status_code = MOI.OTHER_RESULT_STATUS
        m._termination_status_code = MOI.OPTIMIZE_NOT_CALLED

        # Optimality-based bound tightening (OBBT) options
        m._obbt_working_lower_index = Bool[]
        m._obbt_working_upper_index = Bool[]
        m._lower_indx_diff = Bool[]
        m._upper_indx_diff = Bool[]
        m._old_low_index = Bool[]
        m._old_upp_index = Bool[]
        m._new_low_index = Bool[]
        m._new_upp_index = Bool[]
        m._obbt_variables = VI[]
        m._obbt_performed_flag = false

        # Feasibility-based bound tightening options
        m._cp_improvement = 0.0
        m._cp_evaluation_reverse = false

        # Options for adding additional cuts
        m._cut_iterations = 0
        m._cut_add_flag = false

        # Options for repetition
        m._node_repetitions = 0
        m._initial_volume = 0.0
        m._final_volume = 0.0

        # Log
        m._log = Log()

        m._nlp_data = empty_nlp_data()

        m._relaxed_evaluator = Evaluator{1,NS}()
        m._relaxed_constraint_bounds = Vector{MOI.NLPBoundsPair}[]
        m._relaxed_eval_has_objective = false

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

function check_inbounds(m::Optimizer, vi::MOI.VariableIndex)
    num_variables = length(m._variable_info)
    if !(1 <= vi.value <= num_variables)
        error("Invalid variable index $vi. ($num_variables variables in the model.)")
    end
end

check_inbounds(m::Optimizer, var::MOI.SingleVariable) = check_inbounds(m, var.variable)

function check_inbounds(m::Optimizer, aff::MOI.ScalarAffineFunction)
    for term in aff.terms
        check_inbounds(m, term.variable_index)
    end
end

function check_inbounds(m::Optimizer, quad::MOI.ScalarQuadraticFunction)
    for term in quad.affine_terms
        check_inbounds(m, term.variable_index)
    end
    for term in quad.quadratic_terms
        check_inbounds(m, term.variable_index_1)
        check_inbounds(m, term.variable_index_2)
    end
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
    for var in m.integer_variables
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
     m.logging_on = false
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
    check_inbounds(model, vi)
    return model._continuous_solution[vi.value]
end

MOI.supports_constraint(::Optimizer, ::Type{SV}, ::Type{LT}) = true
MOI.supports_constraint(::Optimizer, ::Type{SV}, ::Type{GT}) = true
MOI.supports_constraint(::Optimizer, ::Type{SV}, ::Type{ET}) = true
#MOI.supports_constraint(::Optimizer, ::Type{SV}, ::Type{IT}) = true
#MOI.supports_constraint(::Optimizer, ::Type{SV}, ::Type{MOI.ZO}) = true

MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{LT}) = true
MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{GT}) = true
MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{ET}) = true
#MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{IT}) = true

MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{LT}) = true
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{GT}) = true
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{ET}) = true
#MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{IT}) = true

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
    check_inbounds(m, vi)
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
    check_inbounds(m, vi)
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
    check_inbounds(m, vi)
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
    check_inbounds(m, vi)
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
    check_inbounds(m, vi)
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
            check_inbounds(m, func)
            push!(m.$(array_name), (func, set, length(func.terms)))
            indx = CI{$function_type, $set_type}(length(m.$(array_name)))
            m._constraint_convexity[indx] = true
            return indx
        end
    end
end

macro define_addconstraint_quadratic(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds(m, func)
            for i in func.affine_terms m.branch_variable[i.variable_index.value] = true end
            for i in func.quadratic_terms
                m.branch_variable[i.variable_index_1.value] = true
                m.branch_variable[i.variable_index_2.value] = true
            end
            push!(m.$(array_name), (func, set, length(m.$(array_name))+1))
            indx = CI{$function_type, $set_type}(length(m.$(array_name)))
            m._constraint_convexity[indx] = false
            return indx
        end
    end
end

@define_addconstraint_linear SAF LT _linear_leq_constraints
@define_addconstraint_linear SAF GT _linear_geq_constraints
@define_addconstraint_linear SAF ET _linear_eq_constraints

@define_addconstraint_quadratic SQF LT _quadratic_leq_constraints
@define_addconstraint_quadratic SQF GT _quadratic_geq_constraints
@define_addconstraint_quadratic SQF ET _quadratic_eq_constraints

function MOI.set(m::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    m._nlp_data = nlp_data
    if nlp_data.has_objective
        m._objective_is_sv = false
        m._objective_is_saf = false
        m._objective_is_sqf = false
        m._objective_is_nlp = true
    end
    return
end

MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{SV}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{SAF}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{SQF}) = true

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction, func::SV)
    check_inbounds(m, func)
    m._objective_is_sv = true
    m._objective_is_saf = false
    m._objective_is_sqf = false
    m._objective_is_nlp = false
    m._objective_sv = func
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction, func::SAF)
    check_inbounds(m, func)
    m._objective_is_sv = false
    m._objective_is_saf = true
    m._objective_is_sqf = false
    m._objective_is_nlp = false
    m._objective_saf = func
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction, func::SQF)
    check_inbounds(m, func)
    m._objective_is_sv = false
    m._objective_is_saf = false
    m._objective_is_sqf = true
    m._objective_is_nlp = false
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
function eval_function(var::SV, x)
    return x[var.variable.value]
end

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
        if row_idx == col_idx
            function_value += 0.5*coefficient*x[row_idx.value]*x[col_idx.value]
        else
            function_value += coefficient*x[row_idx.value]*x[col_idx.value]
        end
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
