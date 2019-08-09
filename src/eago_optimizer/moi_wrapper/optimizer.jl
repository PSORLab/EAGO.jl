mutable struct VariableInfo
    is_integer::Bool
    lower_bound::Float64  # May be -Inf even if has_lower_bound == true
    has_lower_bound::Bool # Implies lower_bound == Inf
    upper_bound::Float64  # May be Inf even if has_upper_bound == true
    has_upper_bound::Bool # Implies upper_bound == Inf
    is_fixed::Bool        # Implies lower_bound == upper_bound and !has_lower_bound and !has_upper_bound.
end
VariableInfo() = VariableInfo(false,-Inf, false, Inf, false, false)

lower_bound(x::VariableInfo) = x.lower_bound
upper_bound(x::VariableInfo) = x.upper_bound

@enum FailureLocation LOWER_SOLVER_FAILED UPPER_SOLVER_FAILED PREPROCESSING_FAILED POSTPROCESSING_FAILED NO_FAILURE
@enum OptimizerType LP MILP NLP MINLP

dummy_function() = nothing

export Optimizer
"""
    Optimizer

The main optimizer object used by EAGO to solve problems during the optimization routine. The following
commonly used options are described below and can be set via keyword arguments in the JuMP/MOI model:

* `lower_problem!::Function` - Lower problem function
* `upper_problem!::Function` - Upper problem function
* `preprocess!::Function` - Preprocessing function
* `postprocess!::Function` - Postprocessing function
* `single_check::Function` - Check if single node should be stored
* `convergence_check::Function` - Convergence criterion function
* `termination_check::Function` - Termination check function
* `node_storage!::Function` - Function defining manner in which node is stored to stack
* `node_selection::Function` - Function which selects node from stack
* `bisection_function::Function` - Bisection function
* `cut_condition::Function` - Condition for adding cutting plane
* `add_cut!::Function` - Function for adding additional cutting plane
* `relax_function!::Function`- Function used to relax the model

# Output specification fields
* `verbosity::Int`- Specifies detail of output to console 0 (none) - 4 (node level information)
* `output_iterations::Int`- Number of iterations to skip between printing iteration summary
* `header_iterations::Int`- Number of iterations to skip between printing heade
* `digits_displayed::Int`- Digits displayed before decimal

# Termination Limits
* `iteration_limit::Int`- Maximum number of iterations
* `node_limit::Int`- Maximum number of nodes to store in stack
* `absolute_tolerance::Float64`- Absolute tolerance
* `relative_tolerance::Float64`- Relative tolerance (UBD-LBD)/MIN(|UBD|,|LBD|)

# Optimality-Based Bound Tightening (OBBT) Options
* `obbt_depth::Int`- Depth in B&B tree to perform OBBT on every node to. Probability 2^{depth - obbt_depth} otherwise.
* `obbt_reptitions::Int`- Maximum number of OBBT repetitions to perform.
* `obbt_aggressive_on::Bool`- Use aggressive bound tightening.
* `obbt_aggressive_max_iteration::Int`- Maximum number of iterations for aggressive OBBT.
* `obbt_aggressive_min_dimension::Int`- Minimum dimensions for aggressive OBBT.
* `obbt_tolerance::Float64`- Continue obbt if (Final Volume)/(Initial Volume) < obbt_tolerance

# Duality-Based Bound Tightening (DBBT) Options
* `dbbt_depth::Int`- Depth in B&B tree to perform DBBT on every node to.
* `dbbt_tolerance::Float64`- Tolerance to determine if solution point lies on variable bound.

# Feasibility-Based Bound Tightening Options
* `cp_depth::Int`- Depth in B&B tree to perform constraint propagation on every node to.
* `cp_interval_reptitions::Int`- Maximum number of repetitions for interval constraint propagation.
* `cp_interval_tolerance::Float64`- Continue constraint propagation if (Final Volume)/(Initial Volume) < cp_interval_tolerance.
* `cp_mccormick_reptitions::Int`- Maximum number of repetitions for McCormick constraint propagation.
* `cp_mccormick_tolerance::Float64`- Continue constraint propagation if (Final Volume)/(Initial Volume) < cp_mccormick_tolerance.
* `evaluation_reverse::Bool`- Perform a reverse evaluation pass when evaluating the relaxation.

# Rounding mode (interval arithmetic options)
* `rounding_mode::Symbol`- Rounding mode for interval arithmetic.
* `treat_as_nonlinear::Vector{Int}`- Specifies variable index to treat as nonlinear.
* `cut_max_iterations::Int`-
* `cut_tolerance::Float64`-

# Subgradient Tightening Flag
* `subgrad_tighten::Bool`- Perform an affine contraction on forward McCormick calculations.
* `subgrad_tighten_reverse::Bool`- Perform an affine contraction on reverse McCormick calculations.

# Branching Options
* `mid_cvx_factor::Float64`- Convex factor to sum solution point and midpoint of largest relative dimension.

# Options for Poor Man's LP
* `poor_man_lp_depth::Int`- Depth in B&B tree to perform constraint propagation on every node to
* `poor_man_lp_reptitions::Int`- Maximum number of repetitions for interval constraint propagation

# Options for Quadratic Bounding Tightening
* `univariate_quadratic_depth::Int`- Depth in B&B tree to perform univariate quadratic contraction on every node to
* `univariate_quadratic_reptitions::Int`- Maximum number of repetitions for univariate quadratic contraction
* `bivariate_quadratic_depth::Int`- Depth in B&B tree to perform bivariate quadratic contraction on every node to
* `bivariate_quadratic_reptitions::Int`- Maximum number of repetitions for bivariate quadratic contraction

# Options for Repetition (If DBBT Performed Well)
* `maximum_repetitions::Int`- Maximum number of repetitions to perform all subproblems on
* `repetition_volume_tolerance::Float64`- Tolerance below which to repeat repetitions

# Upper bounding options
* `upper_bounding_depth::Int`- Depth in B&B tree to perform local nlp solve on every node to. Probability 2^{depth - upper_bounding_depth} otherwise.

# Options for specifying Optimizers
* `initial_relaxed_optimizer::MOI.AbstractOptimizer`- Optimizer to use when solving relaxations, initial
* `working_relaxed_optimizer::MOI.AbstractOptimizer`- Optimizer to use when solving relaxations, working
* `initial_upper_optimizer::MOI.AbstractOptimizer`- Optimizer to use when local nlp, initial
* `working_upper_optimizer::MOI.AbstractOptimizer`- Optimizer to use when local nlp, working
* `lower_factory::JuMP.OptimizerFactory`- Optimizer factory to use for creating the optimizer, relaxations
* `upper_factory::JuMP.OptimizerFactory`- Optimizer factory to use for creating the optimizer, local nlp
* `lower_optimizer_options::Dict{Symbol,Any}`- Options to use when constructing problem by factory, relaxations
* `upper_optimizer_options::Dict{Symbol,Any}`- Options to use when constructing problem by factory, local nlp
* `use_lower_factory::Bool`- Use an optimizer factory versus modifying a single optimizer instance, relaxations
* `use_upper_factory::Bool`- Use an optimizer factory versus modifying a single optimizer instance, local nlp

"""
mutable struct Optimizer <: MOI.AbstractOptimizer

    #input_model::Any
    integer_variables::Vector{Int}
    variable_info::Vector{VariableInfo}
    lower_variables::Vector{MOI.VariableIndex}
    upper_variables::Vector{MOI.VariableIndex}
    lower_variable_index::Vector{Tuple{MOI.ConstraintIndex,MOI.ConstraintIndex,Int}}
    upper_variable_index::Vector{Tuple{MOI.ConstraintIndex,MOI.ConstraintIndex,Int}}
    objective_constraint_index::Vector{MOI.ConstraintIndex}
    initial_continuous_values::IntervalBox          # Interval box constraints
    initial_integer_values::Vector{Int}               # Potential Integer Values

    bisection_variable::Dict{Int,Bool}
    fixed_variable::Dict{Int,Bool}
    #=
    pseudo_cost_lower::Vector{Float64}
    pseudo_cost_upper::Vector{Float64}
    prob_count_lower::Vector{Float64}
    prob_count_upper::Vector{Float64}
    =#
    variable_index_to_storage::Dict{Int,Int}
    storage_index_to_variable::Dict{Int,Int}
    constraint_convexity::Dict{MOI.ConstraintIndex,Bool}
    constraint_label::Dict{Int,Symbol}

    continuous_solution::Vector{Float64}             # Stores a point in the IntervalSolutionBox
    integer_solution::Vector{Bool}                   # Stores the integer solution point
    stack::Dict{Int,NodeBB}                       # Map of Node ID to NodeData

    nlp_data
    working_evaluator_block
    variable_number::Int
    state_variables::Int
    continuous_variable_number::Int
    integer_variable_number::Int
    constraint_number::Int
    linear_number::Int
    quadratic_number::Int

    current_lower_info::LowerInfo                       # Problem solution info for lower bounding program
    current_upper_info::UpperInfo                       # Problem solution info for upper bounding program
    current_preprocess_info::PreprocessInfo             # Problem solution info for preprocessing step
    current_postprocess_info::PostprocessInfo           # Problem solution info for postprocessing step

    objective::Union{MOI.SingleVariable,MOI.ScalarAffineFunction{Float64},MOI.ScalarQuadraticFunction{Float64},Nothing}
    objective_convexity::Bool
    custom_mod_flag::Bool

    linear_leq_constraints::Vector{Tuple{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64},Int}}
    linear_geq_constraints::Vector{Tuple{MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64},Int}}
    linear_eq_constraints::Vector{Tuple{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64},Int}}
    linear_interval_constraints::Vector{Tuple{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64},Int}}
    quadratic_leq_constraints::Vector{Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.LessThan{Float64},Int}}
    quadratic_geq_constraints::Vector{Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.GreaterThan{Float64},Int}}
    quadratic_eq_constraints::Vector{Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.EqualTo{Float64},Int}}
    quadratic_interval_constraints::Vector{Tuple{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64},Int}}
    quadratic_convexity::Vector{Bool}

    univariate_quadratic_leq_constraints::Vector{Tuple{Float64,Float64,Float64,Int}}
    univariate_quadratic_geq_constraints::Vector{Tuple{Float64,Float64,Float64,Int}}
    univariate_quadratic_eq_constraints::Vector{Tuple{Float64,Float64,Float64,Int}}
    bivariate_quadratic_leq_constraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}}
    bivariate_quadratic_geq_constraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}}
    bivariate_quadratic_eq_constraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}}

    linear_optimizer::MOI.AbstractOptimizer
    nlp_optimizer::MOI.AbstractOptimizer

    initial_relaxed_optimizer::MOI.AbstractOptimizer
    working_relaxed_optimizer::MOI.AbstractOptimizer
    initial_upper_optimizer::MOI.AbstractOptimizer
    working_upper_optimizer::MOI.AbstractOptimizer
    lower_factory::JuMP.OptimizerFactory
    upper_factory::JuMP.OptimizerFactory
    lower_optimizer_options::Dict{Symbol,Any}
    upper_optimizer_options::Dict{Symbol,Any}
    use_lower_factory::Bool
    use_upper_factory::Bool

    relaxation::RelaxationScheme

    # Stores functions for implementing Branch and Bound
    lower_problem!::Function                                                 # Stores lower problem function
    upper_problem!::Function                                                 # Stores upper problem function
    preprocess!::Function                                                   # Preprocessing function
    postprocess!::Function                                                  # Post processing function
    single_check::Function                                                  # Repeation check
    convergence_check::Function                                             # convergence criterion
    termination_check::Function                                             # Stores termination check function
    node_storage!::Function                                                  # Stores branching function
    node_selection::Function                                                # Stores node selection function
    bisection_function::Function                                            #
    cut_condition::Function                                                 #
    add_cut!::Function                                                       #
    relax_function!::Function                                                # Stores code used to relax the model

    global_lower_bound::Float64                                              # Global Lower Bound
    global_upper_bound::Float64                                              # Global Upper Bound
    maximum_node_id::Int
    history::NodeHistory
    current_iteration_count::Int
    current_node_count::Int

    # Storage for output
    solution_value::Float64                                                 # Value of the solution
    first_found::Bool                                                       #
    feasible_solution_found::Bool                                           #
    first_solution_node::Int                                                #
    last_gap::Float64                                                       #
    optimization_sense::MOI.OptimizationSense                               #
    objective_value::Float64
    termination_status_code::MOI.TerminationStatusCode                      #
    result_status_code::MOI.ResultStatusCode
    started_solve::Bool                                                     #
    failed_solver::FailureLocation                                          # Stores branching function

    # Output specification fields
    verbosity::Int                                                         # Stores output selection
    warm_start::Bool                                                       # Boolean
    output_iterations::Int                                                 # Number of iterations to skip between printing iteration summary
    header_iterations::Int                                                 # Number of iterations to skip between printing heade
    digits_displayed::Int                                                  # digits displayed before decimal
    return_history::Bool                                                   # returns LBD, UBD array and time vector
    flag_subsolver_errors::Bool                                            # If a subsolver has a problem termination code then stop the algorithm
                                                                           # and record it
    # Termination Limits
    iteration_limit::Int                                                   # Maximum number of iterations
    node_limit::Int                                                        # Maximum number of nodes to store in memory
    absolute_tolerance::Float64                                            # Absolute tolerance for BnB
    relative_tolerance::Float64                                             # Relative tolerance for BnB
    exhaustive_search::Bool                                                 # Exhaustive search: find all solns or find first
    local_solve_only::Bool
    feasible_local_continue::Bool

    # Optimality-Based Bound Tightening (OBBT) Options
    obbt_variables::Vector{MOI.VariableIndex}
    obbt_depth::Int
    obbt_reptitions::Int
    obbt_aggressive_on::Bool
    obbt_aggressive_max_iteration::Int
    obbt_aggressive_min_dimension::Int
    obbt_tolerance::Float64
    obbt_working_lower_index::Vector{MOI.VariableIndex}
    obbt_working_upper_index::Vector{MOI.VariableIndex}
    obbt_active_current::Bool
    obbt_performed_flag::Bool

    # Duality-Based Bound Tightening (DBBT) Options
    dbbt_depth::Int
    dbbt_tolerance::Float64

    # Feasibility-Based Bound Tightening Options
    cp_depth::Int
    cp_improvement::Float64
    cp_interval_reptitions::Int
    cp_interval_tolerance::Float64
    cp_mccormick_reptitions::Int
    cp_mccormick_tolerance::Float64
    evaluation_reverse::Bool

    # Rounding mode (interval arithmetic options)
    rounding_mode::Symbol
    treat_as_nonlinear::Vector{Int}
    cut_max_iterations::Int
    cut_tolerance::Float64

    # Subgradient Tightening Flag
    subgrad_tighten::Bool
    subgrad_tighten_reverse::Bool

    # Options for Poor Man's LP
    poor_man_lp_depth::Int
    poor_man_lp_reptitions::Int

    # Options for Quadratic Bounding Tightening
    univariate_quadratic_depth::Int
    univariate_quadratic_reptitions::Int
    bivariate_quadratic_depth::Int
    bivariate_quadratic_reptitions::Int

    # Options for Repetition (If DBBT Performed Well)
    node_repetitions::Int
    maximum_repetitions::Int
    initial_volume::Float64
    final_volume::Float64
    repetition_volume_tolerance::Float64

    # Upper bounding nodes skipped
    upper_bounding_interval::Int #Not used
    upper_bounding_depth::Int
    upper_bnd_ni_cnt::Int  #Not used
    upper_bnd_ni_tol::Int  #Not used
    upper_bnd_this_int::Bool  #Not used

    # Cutting Plane Options
    cut_iterations::Int
    cut_add_flag::Bool
    mid_cvx_factor::Float64

    # Status flags
    first_relaxed::Bool
    upper_has_node::Bool

    # UDF options
    udf_scrubber_flag::Bool
    udf_to_JuMP_flag::Bool

    # Problem reformulation options
    reform_epigraph_flag::Bool
    reform_cse_flag::Bool
    reform_flatten_flag::Bool

    # Debug
    ext

    function Optimizer(;options...)

        m = new()

        default_opt_dict = Dict{Symbol,Any}()

        # set fallback for potentially user defined functions
        for i in (:lower_problem!, :upper_problem!, :preprocess!, :postprocess!, :single_check,
                  :convergence_check, :termination_check, :node_storage!, :node_selection,
                  :bisection_function, :cut_condition, :add_cut!, :relax_function!)
                  default_opt_dict[i] = dummy_function
        end

        # set fallback for optimizers
        for i in (:linear_optimizer, :nlp_optimizer,
                  :initial_relaxed_optimizer, :working_relaxed_optimizer,
                  :initial_upper_optimizer, :working_upper_optimizer)
            default_opt_dict[i] = DummyOptimizer()
        end

        default_opt_dict[:use_upper_factory] = true
        default_opt_dict[:use_lower_factory] = true

        default_opt_dict[:relaxation] = default_relaxation_scheme()
        default_opt_dict[:global_lower_bound] = -Inf
        default_opt_dict[:global_upper_bound] = Inf

        # Output specification fields
        default_opt_dict[:verbosity] = 4
        default_opt_dict[:warm_start] = false
        default_opt_dict[:output_iterations] = 1
        default_opt_dict[:header_iterations] = 10
        default_opt_dict[:digits_displayed] = 3
        default_opt_dict[:return_history] = false
        default_opt_dict[:flag_subsolver_errors] = true

        # Duality-based bound tightening parameters
        default_opt_dict[:dbbt_depth] = Int(1E6)
        default_opt_dict[:dbbt_tolerance] = 1E-8

        # Optimality-based bound tightening parameters
        default_opt_dict[:obbt_depth] = 0
        default_opt_dict[:obbt_aggressive_on] = false
        default_opt_dict[:obbt_aggressive_max_iteration] = 2
        default_opt_dict[:obbt_aggressive_min_dimension] = 2
        default_opt_dict[:obbt_tolerance] = 1E-9
        default_opt_dict[:obbt_reptitions] = 20

        # Feasibility-Based Bound Tightening Options
        default_opt_dict[:cp_depth] = 0
        default_opt_dict[:cp_interval_reptitions] = 0
        default_opt_dict[:cp_interval_tolerance] = 0.99
        default_opt_dict[:cp_mccormick_reptitions] = 0
        default_opt_dict[:cp_mccormick_tolerance] = 0.99
        default_opt_dict[:evaluation_reverse] = false

        # Cutting Options
        default_opt_dict[:cut_max_iterations] = 0

        # Interval options
        default_opt_dict[:rounding_mode] = :accurate

        # Subgradient tightening options
        default_opt_dict[:subgrad_tighten] = true
        default_opt_dict[:subgrad_tighten_reverse] = false

        # Feasibility-Based Bound Tightening for Quadratics
        default_opt_dict[:univariate_quadratic_depth] = 0
        default_opt_dict[:univariate_quadratic_reptitions] = 2
        default_opt_dict[:bivariate_quadratic_depth] = 1000
        default_opt_dict[:bivariate_quadratic_reptitions] = 2

        # Options for Repetition (If DBBT Performed Well)
        default_opt_dict[:maximum_repetitions] = 4
        default_opt_dict[:repetition_volume_tolerance] = 0.9

        # Poor Man's Options
        default_opt_dict[:poor_man_lp_depth] = 0
        default_opt_dict[:poor_man_lp_reptitions] = 1

        # Upper Bounding Interval
        default_opt_dict[:upper_bounding_interval] = 10
        default_opt_dict[:upper_bounding_depth] = 10
        default_opt_dict[:upper_bnd_ni_tol] = 3

        # Termination Limits
        default_opt_dict[:iteration_limit] = 100000 #Int(1E6)
        default_opt_dict[:node_limit] = Int(1E6)
        default_opt_dict[:absolute_tolerance] = 1E-3
        default_opt_dict[:relative_tolerance] = 1E-3
        default_opt_dict[:exhaustive_search] = false
        default_opt_dict[:first_relaxed] = false
        default_opt_dict[:upper_has_node] = false
        default_opt_dict[:local_solve_only] = false
        default_opt_dict[:feasible_local_continue] = false

        default_opt_dict[:mid_cvx_factor] = 0.25

        # UDF options
        default_opt_dict[:udf_scrubber_flag] = true
        default_opt_dict[:udf_to_JuMP_flag] = true

        # Problem Reformulation Options
        default_opt_dict[:reform_epigraph_flag] = false
        default_opt_dict[:reform_cse_flag] = false
        default_opt_dict[:reform_flatten_flag] = false

        default_opt_dict[:treat_as_nonlinear] = Int[]

        for i in keys(default_opt_dict)
            if (haskey(options,i))
                setfield!(m, i, options[i])
            else
                setfield!(m, i, default_opt_dict[i])
            end
        end

        if haskey(options, :lower_optimizer_options)
            lower_options = options[:lower_optimizer_options]
        else
            lower_options = Dict{Symbol,Any}()
        end
        m.lower_optimizer_options = lower_options
        m.lower_factory = JuMP.OptimizerFactory(DummyOptimizer(),(),lower_options)

        if haskey(options, :upper_optimizer_options)
            upper_options = options[:upper_optimizer_options]
        else
            upper_options = Dict{Symbol,Any}()
            # specifies tolerance to force ipopt to terminate after
            # a small number of feasible iterations, or fail after a few
            # iterations. Scales ipopt tolerance to overall tolerance.
            if haskey(options, :absolute_tolerance)
                upper_options[:tol] = options[:absolute_tolerance]/100.0
            else
                upper_options[:tol] = default_opt_dict[:absolute_tolerance]/100.0
            end
            upper_options[:constr_viol_tol] = 0.00001
            upper_options[:max_iter] = 1000
            upper_options[:acceptable_tol] = 1E30
            upper_options[:acceptable_iter] = 100
            upper_options[:acceptable_constr_viol_tol] = 0.0000001
            upper_options[:print_level] = 0
        end
        m.upper_optimizer_options = upper_options
        m.upper_factory = JuMP.OptimizerFactory(DummyOptimizer(),(),upper_options)

        # Set storage used for interior calculations

        m.objective_value = -Inf

        m.ext = []
        #m.input_model = 0
        m.integer_variables = Int[]
        m.variable_info = VariableInfo[]
        m.lower_variables = MOI.VariableIndex[]
        m.upper_variables = MOI.VariableIndex[]
        m.lower_variable_index = Tuple{MOI.ConstraintIndex,MOI.ConstraintIndex,Int}[]
        m.upper_variable_index = Tuple{MOI.ConstraintIndex,MOI.ConstraintIndex,Int}[]
        m.objective_constraint_index = MOI.ConstraintIndex[]
        m.initial_continuous_values = IntervalBox(Interval(0.0))      # Interval box constraints
        m.initial_integer_values = Vector{Int}[]                      # Potential Integer Values
        m.nlp_data = empty_nlp_data()
        m.working_evaluator_block = empty_nlp_data()

        m.linear_number = 0
        m.quadratic_number = 0
        m.cut_iterations = 0
        m.cut_add_flag = true

        m.constraint_convexity = Dict{MOI.ConstraintIndex,Bool}()
        m.variable_index_to_storage = Dict{Int,Int}()
        m.storage_index_to_variable = Dict{Int,Int}()
        m.bisection_variable = Dict{Int,Bool}()
        m.fixed_variable =  Dict{Int,Bool}()
        m.quadratic_convexity = Bool[]
        m.constraint_label = Dict{Int,Symbol}()
        #=
        m.PseudoCostLower = Float64[]
        m.PseudoCostUpper = Float64[]
        m.ProbCountLower = Float64[]
        m.ProbCountUpper = Float64[]
        =#
        m.continuous_solution = Float64[]                            # Stores a point in the IntervalSolutionBox
        m.integer_solution = Bool[]                                  # Stores the integer solution point
        m.stack = Dict{Int,NodeBB}()                              # Map of Node ID to NodeData
        m.variable_number = 0
        m.state_variables = 0
        m.continuous_variable_number = 0
        m.integer_variable_number = 0
        m.constraint_number = 0

        m.linear_leq_constraints = Tuple{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}, Int}[]
        m.linear_geq_constraints = Tuple{MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}, Int}[]
        m.linear_eq_constraints = Tuple{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}, Int}[]
        m.linear_interval_constraints = Tuple{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}, Int}[]

        m.quadratic_leq_constraints = Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.LessThan{Float64}, Int}[]
        m.quadratic_geq_constraints = Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.GreaterThan{Float64}, Int}[]
        m.quadratic_eq_constraints = Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.EqualTo{Float64}, Int}[]
        m.quadratic_interval_constraints = Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.Interval{Float64}, Int}[]

        m.univariate_quadratic_leq_constraints = Tuple{Float64,Float64,Float64,Int}[]
        m.univariate_quadratic_geq_constraints = Tuple{Float64,Float64,Float64,Int}[]
        m.univariate_quadratic_eq_constraints = Tuple{Float64,Float64,Float64,Int}[]
        m.bivariate_quadratic_leq_constraints = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}[]
        m.bivariate_quadratic_geq_constraints = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}[]
        m.bivariate_quadratic_eq_constraints = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}[]

        m.current_lower_info = LowerInfo()
        m.current_upper_info = UpperInfo()
        m.current_preprocess_info = PreprocessInfo()
        m.current_postprocess_info = PostprocessInfo()

        m.objective = nothing
        m.objective_convexity = false
        m.custom_mod_flag = false

        # Historical Information
        m.history = NodeHistory()
        m.current_iteration_count = 0
        m.current_node_count = 0
        m.maximum_node_id = 0

        m.node_repetitions = 0
        m.initial_volume = 0.0
        m.final_volume = 0.0

        # Output for Solution Storage
        m.solution_value = -Inf
        m.first_found = false
        m.feasible_solution_found = false
        m.first_solution_node = -1
        m.last_gap = -Inf
        m.optimization_sense = MOI.FEASIBILITY_SENSE
        m.termination_status_code = MOI.OPTIMIZE_NOT_CALLED #MOI.OptimizeNotCalled
        m.result_status_code = MOI.OTHER_RESULT_STATUS
        m.started_solve = false
        m.failed_solver = NO_FAILURE

        # Optimality-Based Bound Tightening (OBBT) Storage
        m.obbt_performed_flag = false
        m.obbt_variables = MOI.VariableIndex[]
        m.obbt_working_lower_index = MOI.VariableIndex[]
        m.obbt_working_upper_index = MOI.VariableIndex[]

        # NLP Upper Bound Heurestics
        m.upper_bnd_ni_cnt = 0
        m.upper_bnd_this_int = false

        return m
    end
end
