Base.@kwdef mutable struct GlobalOptimizer{R,Q,S<:ExtensionType} <: MOI.AbstractOptimizer
    relaxed_optimizer::R = Incremental(GLPK.Optimizer())
    upper_optimizer::Q = Incremental(Ipopt.Optimizer())

    _termination_status_code::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _result_status_code::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS

    _working_problem::ParsedProblem = ParsedProblem()

    _end_state::GlobalEndState = GS_UNSET
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

    # Storage for output, reset in initial_parse! in parse.jl
    _solution_value::Float64 = 0.0
    _feasible_solution_found::Bool = false
    _first_solution_node::Int = -1
    _best_upper_value::Float64 = Inf
   
    _cut_iterations::Int = 0
    _cut_add_flag::Bool = false
    _node_repetitions::Int = 0
   
    _log::Log = Log()
   
    _affine_relax_ci::Vector{CI{SAF,LT}} = CI{SAF,LT}[]
    _affine_objective_cut_ci::Union{CI{SV,LT},CI{SAF,LT},Nothing} = nothing
   
    _relaxed_variable_number::Int = 0
    _relaxed_variable_index::Vector{VI} = VI[]
    _relaxed_variable_eq::Vector{Tuple{CI{SV, ET}, Int}} = Tuple{CI{SV, ET}, Int}[]
    _relaxed_variable_lt::Vector{Tuple{CI{SV, LT}, Int}} = Tuple{CI{SV, LT}, Int}[]
    _relaxed_variable_gt::Vector{Tuple{CI{SV, GT}, Int}} = Tuple{CI{SV, GT}, Int}[]
   
    _branch_variables::Vector{Bool} = Bool[]
   
    _new_eval_constraint::Bool = false
    _new_eval_objective::Bool = false
   
    _node_to_sv_leq_ci::Vector{CI{SV,LT}} = CI{SV,LT}[]
    _node_to_sv_geq_ci::Vector{CI{SV,GT}} = CI{SV,GT}[]
    _nonlinear_evaluator_created::Bool = false
   
    _branch_cost::BranchCostStorage{Float64} = BranchCostStorage{Float64}()
    _branch_variable_sparsity::SparseMatrixCSC{Bool,Int} = spzeros(Bool,1,1)
    _constraint_infeasiblity::Vector{Float64} = Float64[]

    # reset in initial_parse! in parse.jl
    _min_converged_value::Float64 = Inf
    _global_lower_bound::Float64 = -Inf
    _global_upper_bound::Float64 = Inf
    _maximum_node_id::Int = 0
    _iteration_count::Int = 0
    _node_count::Int = 0
end