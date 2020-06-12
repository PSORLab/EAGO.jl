function build_nlp_kernel!(d::Evaluator{N,T}, src::JuMP.NLPEvaluator, x::Optimizer, bool_flag::Bool) where {N,T<:RelaxTag}

    m = src.m::Model
    num_variables_ = JuMP.num_variables(m)
    d.variable_number = num_variables_

    nldata = deepcopy(m.nlp_data)
    parameter_values = nldata.nlparamvalues

    # Copy state of user-defined multivariable functions
    d.has_user_mv_operator = src.disable_2ndorder
    d.last_x = fill(NaN, d.variable_number)
    d.last_node = NodeBB()

    # Set valued operators must contain a (sub)gradient and/or (sub)hessian field to support 1st/2nd order eval
    d.disable_1storder = false

    # Add objective functions, constraints, subexpressions
    d.has_nlobj = src.has_nlobj
    if src.has_nlobj
        copy_to_function!(d, 1, src.objective)
    end

    for i in 1:length(src.constraints)
        copy_to_function!(d, i + 1, src.constraints[i])
    end

    d.subexpression_order = src.subexpression_order
    d.subexpression_linearity = src.subexpression_linearity

    d.subexpression_values_set = MC{N,T}[]
    d.subexpression_values_flt = Float64[]
    d.subexpressions = SubexpressionSetStorage{N}[]
    for i in 1:length(src.subexpressions)
        copy_to_subexpr!(d, src.subexpressions[i])
    end
    d.subexpression_values_set = fill(zero(MC{N,T}), length(d.subexpressions))
    d.subexpression_values_flt = fill(NaN, length(d.subexpressions))

    # Add bounds to evaluator
    for bnds in x._nlp_data.constraint_bounds
        push!(d.constraints_lbd, bnds.lower)
        push!(d.constraints_ubd, bnds.upper)
    end

    d.cp_tolerance = x.cp_tolerance
    d.cp_repetitions = x.cp_repetitions
    d.has_reverse = x._cp_evaluation_reverse
    d.subgrad_tighten = x.subgrad_tighten
    d.subgrad_tighten_reverse = x.subgrad_tighten_reverse
    d.jac_storage = fill(zero(MC{N,T}), max(num_variables_, nldata.largest_user_input_dimension))
    d.flt_jac_storage = fill(0.0, max(num_variables_, nldata.largest_user_input_dimension))

    d.constraint_number = length(d.constraints)
    d.subexpression_number = length(d.subexpressions)

    # calculate an index for each variable via search on
    unvisited_tuple = (-1,-1,-1)
    d.index_to_variable = fill(unvisited_tuple, (d.variable_number,))
    for (indx, node) in enumerate(d.objective.nd)
        op = node.index
        ntype = node.nodetype
        if (ntype == JuMP._Derivatives.VARIABLE)
            current_value = d.index_to_variable[op]
            if (current_value == unvisited_tuple)
                d.index_to_variable[op] = (indx, 1, 1)
            end
            @inbounds d.objective.numvalued[indx] = false
        elseif ntype == JuMP._Derivatives.VALUE
            @inbounds d.objective.numberstorage[indx] = d.objective.const_values[op]
            @inbounds d.objective.numvalued[indx] = true
        elseif ntype == JuMP._Derivatives.PARAMETER
            @inbounds d.objective.numberstorage[indx] = parameter_values[op]
            @inbounds d.objective.numvalued[indx] = true
        end
    end
    for (cindx,constraint) in enumerate(d.constraints)
        for (indx,node) in enumerate(constraint.nd)
            op = node.index
            if node.nodetype == JuMP._Derivatives.VARIABLE
                current_value = d.index_to_variable[op]
                if (current_value == unvisited_tuple)
                    d.index_to_variable[op] = (indx, cindx, 2)
                end
                @inbounds constraint.numvalued[indx] = false
            elseif node.nodetype == JuMP._Derivatives.VALUE
                @inbounds constraint.numberstorage[indx] = constraint.const_values[op]
                @inbounds constraint.numvalued[indx] = true
            elseif node.nodetype == JuMP._Derivatives.PARAMETER
                @inbounds constraint.numberstorage[indx] = parameter_values[op]
                @inbounds constraint.numvalued[indx] = true
            end
        end
    end
    for (cindx,subexpress) in enumerate(d.subexpressions)
        for (indx,node) in enumerate(subexpress.nd)
            op = node.index
            if node.nodetype == JuMP._Derivatives.VARIABLE
                current_value = d.index_to_variable[op]
                if (current_value == unvisited_tuple)
                    d.index_to_variable[op] = (indx, cindx, 3)
                end
                @inbounds subexpress.numvalued[indx] = false
            elseif node.nodetype == JuMP._Derivatives.VALUE
                @inbounds subexpress.numberstorage[indx] = subexpress.const_values[op]
                @inbounds subexpress.numvalued[indx] = true
            elseif node.nodetype == JuMP._Derivatives.PARAMETER
                @inbounds subexpress.numberstorage[indx] = parameter_values[op]
                @inbounds subexpress.numvalued[indx] = true
            end
        end
    end

    d.subexpression_isnum = fill(true, (d.subexpression_number,))

    d.user_operators = nldata.user_operators
    return
end

"""
$(TYPEDSIGNATURES)

Builds the evaluator used to generate relaxations of the nonlinear equations
and constraints from a source model.
"""
function build_nlp_evaluator(N::Int64, s::T, src::JuMP.NLPEvaluator, x::Optimizer, bool_flag::Bool) where {T<:RelaxTag}

    # Creates the empty evaluator
    d::Evaluator{N,T} = Evaluator{N,T}()
    build_nlp_kernel!(d, src, x, bool_flag)

    return d
end
=#

#=
"""
$(TYPEDSIGNATURES)

Disables bound tightening for problems lacking the appropriate constraints.
"""
function check_disable_fbbt!(m::Optimizer)

    no_constraints = true
    if length(m._nlp_data.constraint_bounds) > 0
        no_constraints &= false
    end
    if no_constraints
        m.cp_depth = -1
    end

    no_quad_constraints = true
    if m._input_problem._quadratic_leq_count !== 0
        no_constraints &= false
        no_quad_constraints &= false
    end
    if m._input_problem._quadratic_geq_count !== 0
        no_constraints &= false
        no_quad_constraints &= false
    end
    if m._input_problem._quadratic_eq_count !== 0
        no_constraints &= false
        no_quad_constraints &= false
    end
    if no_quad_constraints
        m.quad_uni_depth = -1
        m.quad_bi_depth = -1
    end

    only_lin_constraints = no_constraints
    no_lin_constraints = true
    if m._input_problem._linear_leq_count !== 0
        no_lin_constraints &= false
    end
    if m._input_problem._linear_geq_count !== 0
        no_lin_constraints &= false
    end
    if m._input_problem._linear_eq_count !== 0
        no_lin_constraints &= false
    end
    if no_lin_constraints
        m.lp_depth = -1
    end

    if only_lin_constraints
        m.obbt_depth = -1
    end

    return
end
=#

#=
function has_evaluator(x::MOI.NLPBlockData)
    flag = x.evaluator !== nothing
    flag &= ~isa(x.evaluator, EmptyNLPEvaluator)
    return flag
end

function initialize_scrub!(m::Optimizer, y::JuMP.NLPEvaluator)
    m.presolve_scrubber_flag && Script.scrub!(y.m.nlp_data)
    if m.presolve_to_JuMP_flag
        Script.udf_loader!(m)
    end
    return
end

function initialize_evaluators!(m::Optimizer, flag::Bool)

    nlp_data = deepcopy(m._nlp_data)

    has_eval = has_evaluator(nlp_data)
    if has_evaluator(nlp_data)

        # Build the JuMP NLP evaluator
        evaluator = nlp_data.evaluator::JuMP.NLPEvaluator
        num_nlp_constraints = length(nlp_data.constraint_bounds)
        features = MOI.features_available(evaluator)
        has_hessian = (:Hess in features)
        init_feat = [:Grad, :ExprGraph]
        num_nlp_constraints > 0 && push!(init_feat, :Jac)
        MOI.initialize(evaluator, init_feat)
        m._nlp_data = nlp_data

        # Scrub user-defined functions
        initialize_scrub!(m, evaluator)

        built_evaluator = build_nlp_evaluator(m._variable_number, NS(), deepcopy(evaluator), m, false)
        m._relaxed_evaluator = built_evaluator
        m._relaxed_eval_has_objective = m._nlp_data.has_objective
        append!(m._relaxed_constraint_bounds, m._nlp_data.constraint_bounds)

        # add info to guard context
        m._relaxed_evaluator.ctx = GuardCtx(metadata = GuardTracker(m.domain_violation_ϵ))
    end



    #m.nlp_data.evaluator = evaluator #TODO: Rebuilt entire nlp_block...

    # Transform UDFs to JuMP ASTs
    #m.udf_to_JuMP_flag && Script.udf_loader!(m)
    # Rebuild the nlp-evaluator with udfs -> JuMP expressions

    #####
    #Script to unpack UDF from here
    #unpacked_evaluator = script_to_dag()
    #m.nlp.evaluator = unpacked_evaluator
    #####

    return
end

#=
is_lp(m::Optimizer) = ~in(true, m.branch_variable)

function linear_solve!(m::Optimizer)

    opt = m.relaxed_optimizer
    if m._objective_type === SINGLE_VARIABLE
        MOI.set(opt, MOI.ObjectiveFunction{SV}(), m._objective_sv)
    elseif  m._objective_type === SCALAR_AFFINE
        MOI.set(opt, MOI.ObjectiveFunction{SAF}(), m._objective_saf)
    end

    MOI.optimize!(opt)
    m._objective_value = MOI.get(opt, MOI.ObjectiveValue())
    m._solution_value = MOI.get(opt, MOI.ObjectiveValue())
    m._global_lower_bound = MOI.get(opt, MOI.ObjectiveValue())
    m._global_upper_bound = MOI.get(opt, MOI.ObjectiveValue())
    m._termination_status_code = MOI.get(opt, MOI.TerminationStatus())
    m._result_status_code = MOI.get(opt, MOI.PrimalStatus())
    m._continuous_solution = MOI.get.(opt, MOI.VariablePrimal(), m._lower_variable_index)
    #m._run_time = MOI.get(opt, MOI.SolveTime())

    return
end

=#

#=
"""
$(FUNCTIONNAME)

A rountine that relaxes all nonlinear constraints excluding
constraints specified as quadratic.
"""
function relax_nlp!(x::Optimizer, v::Vector{Float64}, q::Int64)

    evaluator = x._relaxed_evaluator

    if ~isempty(x.branch_variable)

        if MOI.supports(x.relaxed_optimizer, MOI.NLPBlock())

            _nlp_data = MOI.NLPBlockData(x._nlp_data.constraint_bounds,
                                         evaluator,
                                         x._nlp_data.has_objective)
            MOI.set(x._relaxed_optimizer, MOI.NLPBlock(), _nlp_data)

        else
            # Add other affine constraints
            constraint_bounds = x._relaxed_constraint_bounds
            leng = length(constraint_bounds)
            if leng > 0
                nx = x._variable_number
                vi = x._lower_variable_index

                g = zeros(leng)
                dg = zeros(leng, nx)

                g_cc = zeros(leng)
                dg_cc = zeros(leng, nx)

                MOI.eval_constraint(evaluator, g, v)
                MOI.eval_constraint_jacobian(evaluator, dg, v)

                eval_constraint_cc(evaluator, g_cc, v)
                eval_constraint_cc_grad(evaluator, dg_cc, v)

                lower_nlp_affine = x._lower_nlp_affine[q]
                upper_nlp_affine = x._upper_nlp_affine[q]
                lower_nlp_sparsity = x._lower_nlp_sparsity
                upper_nlp_sparsity = x._upper_nlp_sparsity
                lower_nlp_affine_indx = x._lower_nlp_affine_indx
                upper_nlp_affine_indx = x._upper_nlp_affine_indx
                if (q == 1) & x.relaxed_inplace_mod
                    for i = 1:length(lower_nlp_affine_indx)
                        @inbounds g_indx = lower_nlp_affine_indx[i]
                        @inbounds aff_ci = lower_nlp_affine[i]
                        @inbounds nzidx = lower_nlp_sparsity[i]
                        @inbounds nzvar = vi[nzidx]
                        @inbounds constant = g[g_indx]
                        dg_cv_val = 0.0
                        for j in nzidx
                            @inbounds dg_cv_val = dg[i,j]
                            @inbounds vindx = vi[j]
                            @inbounds constant -= v[j]*dg_cv_val
                            MOI.modify(x.relaxed_optimizer, aff_ci, SCoefC(vindx, dg_cv_val))
                        end
                        set = LT(-constant)
                        MOI.set(x.relaxed_optimizer, MOI.ConstraintSet(), aff_ci, set)
                    end
                    for i = 1:length(upper_nlp_affine_indx)
                        @inbounds g_indx = upper_nlp_affine_indx[i]
                        @inbounds aff_ci = upper_nlp_affine[i]
                        @inbounds nzidx = upper_nlp_sparsity[i]
                        @inbounds nzvar = vi[nzidx]
                        @inbounds constant = g_cc[g_indx]
                        dg_cc_val = 0.0
                        for j in nzidx
                            @inbounds dg_cc_val = -dg_cc[i,j]
                            @inbounds vindx = vi[j]
                            @inbounds constant += v[j]*dg_cc_val
                            MOI.modify(x.relaxed_optimizer, aff_ci, SCoefC(vindx, dg_cc_val))
                        end
                        set = LT(constant)
                        MOI.set(x.relaxed_optimizer, MOI.ConstraintSet(), aff_ci, set)
                    end
                else
                    for i = 1:length(lower_nlp_affine_indx)
                        @inbounds g_indx = lower_nlp_affine_indx[i]
                        @inbounds aff_ci = lower_nlp_affine[i]
                        @inbounds nzidx = lower_nlp_sparsity[i]
                        @inbounds nzvar = vi[nzidx]
                        @inbounds constant = g[g_indx]
                        dg_cv_val = 0.0
                        coeff = zeros(Float64,length(nzidx))
                        vindices = vi[nzidx]
                        for j in 1:length(nzidx)
                            @inbounds indx = nzidx[j]
                            @inbounds coeff[j] = dg[i,indx]
                            @inbounds constant -= v[indx]*coeff[j]
                        end
                        set = LT(-constant)
                        saf = SAF(SAT.(coeff,vindices), 0.0)
                        x._lower_nlp_affine[q][i] = MOI.add_constraint(x.relaxed_optimizer,
                                                                   saf, set)
                    end
                    for i = 1:length(upper_nlp_affine_indx)
                        @inbounds g_indx = upper_nlp_affine_indx[i]
                        @inbounds aff_ci = upper_nlp_affine[i]
                        @inbounds nzidx = upper_nlp_sparsity[i]
                        @inbounds nzvar = vi[nzidx]
                        @inbounds constant = g_cc[g_indx]
                        dg_cc_val = 0.0
                        coeff = zeros(Float64,length(nzidx))
                        @inbounds vindices = vi[nzidx]
                        for j in 1:length(nzidx)
                            @inbounds indx = nzidx[j]
                            @inbounds dg_cc_val = -dg_cc[i,indx]
                            @inbounds coeff[j] = dg_cc_val
                            @inbounds constant += v[indx]*dg_cc_val
                        end
                        set = LT(constant)
                        saf = SAF(SAT.(coeff,vindices), 0.0)
                        x._upper_nlp_affine[q][i] = MOI.add_constraint(x.relaxed_optimizer,
                                                                   saf, set)
                    end
                end
            end
        end
    end
    return nothing
end
=#


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
# src/eago_optimizer/evaluator/load.jl
# Create FunctionSetStorage and SubexpressionSetStorage from JuMP nonlinear
# backend. We'll evenetually want to build this from an expression drag instead.
#############################################################################


function copy_to_function!(d::Evaluator{N,T}, i::Int64, y::JuMP._FunctionStorage) where {N,T<:RelaxTag}
    x = deepcopy(y)
    lenx = length(x.nd)
    temp_set = fill(MC{N,T}(Interval(-Inf, Inf)), (lenx,))
    temp_flt = Array{Float64}(undef, lenx)
    temp_bool = Array{Bool}(undef, lenx)

    tpdict = Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}()
    tp1_count = 0
    tp2_count = 0
    tp3_count = 0
    tp4_count = 0
    for i = 1:lenx
        op = x.nd[i].index
        if double_tp(op)
            tp1_count += 1
            tp2_count += 1
            tp3_count += 1
            tp4_count += 1
            tpdict[i] = (tp1_count, tp2_count, tp3_count, tp4_count)
        elseif single_tp(op)
            tp1_count += 1
            tp2_count += 1
            tpdict[i] = (tp1_count, tp2_count, -1, -1)
        end
    end
    tp1storage = zeros(tp1_count)
    tp2storage = zeros(tp2_count)
    tp3storage = zeros(tp3_count)
    tp4storage = zeros(tp4_count)

    sto = FunctionSetStorage{N,T}(x.nd, x.adj, x.const_values, temp_set, temp_flt,
                                temp_bool, tp1storage, tp2storage, tp3storage, tp4storage, tpdict,
                                x.grad_sparsity, x.dependent_subexpressions)
    if i == 1
        d.objective = sto
    else
        push!(d.constraints, sto)
    end
    return
end
function copy_to_subexpr!(d::Evaluator{N,T}, y::JuMP._SubexpressionStorage) where {N,T<:RelaxTag}
    x = deepcopy(y)
    lenx = length(x.nd)
    temp_set = fill(MC{N,T}(Interval(-Inf, Inf)), (lenx,))
    temp_flt = Array{Float64}(undef, lenx)
    temp_bool = Array{Bool}(undef, lenx)

    tpdict = Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}()
    tp1_count = 0
    tp2_count = 0
    tp3_count = 0
    tp4_count = 0
    for i = 1:lenx
        op = x.nd[i].index
        if double_tp(op)
            tp1_count += 1
            tp2_count += 1
            tp3_count += 1
            tp4_count += 1
            tpdict[i] = (tp1_count, tp2_count, tp3_count, tp4_count)
        elseif single_tp(op)
            tp1_count += 1
            tp2_count += 1
            tpdict[i] = (tp1_count, tp2_count, -1, -1)
        end
    end
    tp1storage = zeros(tp1_count)
    tp2storage = zeros(tp2_count)
    tp3storage = zeros(tp3_count)
    tp4storage = zeros(tp4_count)

    sto = SubexpressionSetStorage{N,T}(x.nd, x.adj, x.const_values, temp_set, temp_flt,
                                     temp_bool, tp1storage, tp2storage, tp3storage, tp4storage, tpdict,
                                     x.linearity)
    push!(d.subexpressions, sto)
end


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
# src/eago_optimizer/evaluator/evaluator.jl
# Structures used store nonlinear functions used in computing relaxations.
#############################################################################

"""
$(TYPEDEF)

A storage object for both set and number valued data required to
compute relaxations which contains the tape used to compute a nonlinear function.
The object is parameterized by a `{N,T<:RelaxTag}` where N corresponds the
subgradient size used in the MC object.

$(TYPEDFIELDS)
"""
mutable struct FunctionSetStorage{N, T<:RelaxTag}
    "List of nodes in nonlinear expression"
    #nd::Vector{JuMP.NodeData}
    "Adjacency Matrix for the expression"
    adj::SparseMatrixCSC{Bool,Int64}
    #const_values::Vector{Float64}
    #setstorage::Vector{MC{N,T}}
    numberstorage::Vector{Float64}
    numvalued::Vector{Bool}
    tp1storage::Vector{Float64}
    tp2storage::Vector{Float64}
    tp3storage::Vector{Float64}
    tp4storage::Vector{Float64}
    tpdict::Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}
    grad_sparsity::Vector{Int64}
    #dependent_subexpressions::Vector{Int64}
end

FunctionSetStorage(N,T) = FunctionSetStorage{N,T}(JuMP.NodeData[],spzeros(Bool,1),
                                           Float64[],MC{N,T}[],Float64[], Bool[],
                                           Float64[], Float64[], Float64[], Float64[],
                                           Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}(), Int[],Int[])

"""
$(TYPEDEF)

A storage object for both set and number valued data required to
compute relaxations  which contains the tape used to compute a nonlinear
subexpression. The object is parameterized by a `{N,T<:RelaxTag}` where
N corresponds the the subgradient size used in the MC object.

$(TYPEDFIELDS)
"""
mutable struct SubexpressionSetStorage{N,T<:RelaxTag}
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int64}
    const_values::Vector{Float64}
    setstorage::Vector{MC{N,T}}
    numberstorage::Vector{Float64}
    numvalued::Vector{Bool}
    tp1storage::Vector{Float64}
    tp2storage::Vector{Float64}
    tp3storage::Vector{Float64}
    tp4storage::Vector{Float64}
    tpdict::Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}
    linearity::JuMP._Derivatives.Linearity
end

"""
$(TYPEDEF)

MOI.AbstractNLPEvaluator for calculating relaxations of nonlinear terms.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Evaluator <: MOI.AbstractNLPEvaluator

    user_operators::JuMP._Derivatives.UserOperatorRegistry = JuMP._Derivatives.UserOperatorRegistry()
    has_user_mv_operator::Bool = false
    parameter_values::Vector{Float64} = Float64[]

    current_node::NodeBB = NodeBB()
    lower_variable_bounds::Vector{Float64} = Float64[]
    upper_variable_bounds::Vector{Float64} = Float64[]
    x_value::Vector{Float64} = Float64[]

    "Context used to guard against domain violations & branch on these violations if necessary"
    subgrad_tighten::Bool = false
    subgrad_tighten_reverse::Bool = false
    ctx::GuardCtx = GuardCtx()

    subexpressions::Vector{SubexpressionSetStorage{N,T}}
end

"""
$(FUNCTIONNAME)

Sets the current node in the Evaluator structure.
"""
function set_node!(x::Evaluator, n::NodeBB)
    x.current_node = NodeBB(n)
    return nothing
end

# TODO: Unpacks variable bounds....
get_node(d::Evaluator) = d.current_node

include("univariate.jl")
include("passes.jl")
include("get_info.jl")
include("load.jl")

#=
# WORK ON NEW EVALUATOR
struct NonlinearFunction
    d::Evaluator
    indx::
end

struct EvaluatorParams
    parameter_values::Vector{Float64}
    has_user_mv_operator::Bool
    has_nlobj::Bool
    has_reverse::Bool
    subgrad_tighten::Bool
    subgrad_tighten_reverse::Bool
    cp_repetitions::Int64
    cp_tolerance::Float64
    "Context used to guard against domain violations & branch on these violations if necessary"
    ctx::GuardCtx
end
=#

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
# src/eago_optimizer/evaluator/get_info.jl
# Access functions for information from evaluator.
#############################################################################

function MOI.eval_objective(d::Evaluator, x::Vector{Float64})
    forward_reverse_pass(d, x)
    val = 0.0
    if d.has_nlobj
        if d.objective.numvalued[1]
            val = d.objective.numberstorage[1]
        else
            val = d.objective.setstorage[1].cv
        end
    else
        error("No nonlinear objective.")
    end
    return val
end

"""
$(FUNCTIONNAME)

Retrieves the lower bound of the objective.
"""
function eval_objective_lo(d::Evaluator)
    val = 0.0
    if d.has_nlobj
        if d.objective.numvalued[1]
            val = d.objective.numberstorage[1]
        else
            val = d.objective.setstorage[1].Intv.lo
        end
    else
        error("No nonlinear objective.")
    end
    return val
end

function MOI.eval_constraint(d::Evaluator, g::Vector{Float64}, x::Vector{Float64})
    forward_reverse_pass(d, x)

    n = d.current_node

    print_info = false
    print_info && println("n = $n")
    #objvar = MC{11,NS}(x[1], Interval{Float64}(n.lower_variable_bounds[1], n.upper_variable_bounds[1]), 1)
    xMC = [MC{11,NS}(x[i], Interval{Float64}(n.lower_variable_bounds[i], n.upper_variable_bounds[i]), i) for i in 1:10]

    print_info && println("d.constraints[1].setstorage: $(length(d.constraints[1].nd))")
    for i = 1:length(d.constraints)
        if i === 1
            print_info && println(" ")
            print_info && println("i = $i, setstorage = $(d.constraints[i].setstorage[1])")
            #val = -xMC[1]*(1.12 + 0.13167*xMC[8] - 0.00667*(xMC[8])^2)+xMC[4] # interval and cc are different
            val9 = 0.0
            print_info && println("val9 = $(val9)")
            val8 = xMC[1]
            print_info && println("val8 = $(val8), type = $(typeof(val8))")
            val7 = xMC[8]
            print_info && println("val7 = $(val7)")
            val6 = 0.13167
            print_info && println("val6 = $(val6)")
            val5 = val6*val7
            print_info && println("val5 = $(val5)")
            val4 = 1.12
            print_info && println("val4 = $(val4)")
            val3 = val4 - val5
            print_info && println("val3 = $(val3), type = $(typeof(val3))")
            val2 = val8*val3
            print_info && println("val2 = $(val2), type = $(typeof(val2))")
            val1 = xMC[1]*(1.12 - 0.13167*xMC[8])
            print_info && println("val1 = $(val1)")

            #=
        elseif i === 2
            val = -0.001*xMC[4]*xMC[9]*xMC[6]/(98 - xMC[6]) + xMC[3]
        elseif i === 3
            #val = -(1.098*xMC[8] - 0.038*(xMC[8])^2) - 0.325*xMC[6] + xMC[7] - 57.425 # INTERVAL BOUNDS ARE DIFFERENT
            val = -1.098*xMC[8] + xMC[8]*xMC[8] - 0.325*xMC[6] + xMC[7] - 57.425
        elseif i === 4
            val = -(xMC[2] + xMC[5])/xMC[1] + xMC[8]
        elseif i === 5
            val = -0.063*xMC[4]*xMC[7] + 5.04*xMC[1] + 0.035*xMC[2] + 10*xMC[3] + 3.36*xMC[5] - objvar
            =#
        end
        #println("i = $i, val = $(val), val0 = $(val0)")
        print_info && println(" ")
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].cv
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Populates `g` with the concave relaxations of the constraints of `d` evaluated
at `y`.
"""
function eval_constraint_cc(d::Evaluator, g::Vector{Float64}, y::Vector{Float64})
    forward_reverse_pass(d,y)
    for i = 1:length(d.constraints)
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].cc
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Populates `g` with the lower bounds of the constraints of `d`.
"""
function eval_constraint_lo!(d::Evaluator, g::Vector{Float64})
    for i = 1:length(d.constraints)
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].Intv.lo
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Populates `g` with the upper bounds of the constraints of `d`.
"""
function eval_constraint_hi!(d::Evaluator, g::Vector{Float64})
    for i = 1:length(d.constraints)
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].Intv.hi
        end
    end
    return
end


function MOI.eval_objective_gradient(d::Evaluator, df::Vector{Float64}, x::Vector{Float64})
    forward_reverse_pass(d,x)
    if d.has_nlobj
        if ~d.objective.numvalued[1]
            for j = 1:length(d.objective.setstorage[1].cv_grad)
                df[j] = d.objective.setstorage[1].cv_grad[j]
            end
        end
    else
        error("No nonlinear objective.")
    end
    return
end

function MOI.eval_constraint_jacobian(d::Evaluator,g,x)
    forward_reverse_pass(d,x)
    for i = 1:length(d.constraints)
        if ~d.constraints[i].numvalued[1]
            for j = 1:d.variable_number
                g[i,j] = d.constraints[i].setstorage[1].cv_grad[j]
            end
        else
            for j = 1:length(d.constraints[i].setstorage[1].cv_grad)
                g[i,j] = 0.0
            end
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Populates `g` with the subgradients of the constraints of `d` evaluated at `y`.
"""
function eval_constraint_cc_grad(d::Evaluator, g, y)
        forward_reverse_pass(d,y)
        for i = 1:length(d.constraints)
            if ~d.constraints[i].numvalued[1]
                for j = 1:d.variable_number
                    g[i,j] = d.constraints[i].setstorage[1].cc_grad[j]
                end
            else
                for j = 1:d.variable_number
                    g[i,j] = 0.0
                end
            end
        end
    return
end

# looks good
function MOI.features_available(d::Evaluator)
    features = Symbol[]
    if !d.disable_1storder
        push!(features,:Grad)
        push!(features,:Jac)
    end
    return features
end

# looks good, doesn't do anything, EAGO builds the evaluator and attaches it to lower problems
function MOI.initialize(d::Evaluator, requested_features::Vector{Symbol})
end

# TO DO (CHECK GRADIENT DIMS)
function MOI.eval_constraint_jacobian_product(d::Evaluator, y, x, w)
    if !d.disable_1storder
        forward_reverse_pass(d,x)
        t = typeof(d.constraints[1].setstorage[1])
        y = zeros(t,length(d.constraints[1].setstorage[1].cv_grad),length(d.constraints))
        for i = 1:length(d.constraints)
            if ~d.constraints[i].numvalued[1]
                for j = 1:d.variable_number
                    y[i] += d.constraints[i].setstorage[1].cv_grad[j]*w[j]
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# TO DO
function MOI.eval_constraint_jacobian_transpose_product(d::Evaluator, y, x, w)
    if !d.disable_1storder
        forward_reverse_pass(d,x)
        y = zeros(Float64,length(d.constraints[1].setstorage[1].cv_grad),length(d.constraints))
        for i in 1:length(d.constraints)
            if ~d.constraints[i].numvalued[1]
                for j = 1:d.variable_number
                    y[i] += d.constraints[i].setstorage[1].cv_grad[j]*w[j]
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

function MOI.jacobian_structure(d::Evaluator)
    jacobian_sparsity = Tuple{Int64,Int64}[]
    for row = 1:length(d.constraints)
        row_sparsity = d.constraints[row].grad_sparsity
        for idx in row_sparsity
            push!(jacobian_sparsity, (row, idx))
        end
    end
    return jacobian_sparsity
end

function grad_sparsity(d::Evaluator, j::Int64)
    sparsity = Int64[]
    if j == 1
        sparsity = d.objective.grad_sparsity
    else
        sparsity = d.constraints[j-1].grad_sparsity
    end
    return sparsity
end

#=
const SCALAR_ACCESS_FUNCTIONS = Union{cv, cc, lo, hi}
function eval_objective(::typeof{f}, d::NonlinearFunction{N}, x::Vector{Float64}) where {f <: SCALAR_ACCESS_FUNCTIONS, N}
    @assert d.has_nlobj "No nonlinear objective."
    forward_reverse_pass(d, x, 0)
    d.numvalued[1] ? d.numberstorage[1] : f(d.setstorage[1])::Float64
end
function eval_constraint(::typeof{f}, d::NonlinearFunction{N}, x::Vector{Float64}, i::Int) where {f <: SCALAR_ACCESS_FUNCTIONS, N}
    forward_reverse_pass(d, x, i)
    d.numvalued[1] ? d.numberstorage[1] : f(d.setstorage[1])
end

const VECTOR_ACCESS_FUNCTIONS = Union{cv_grad, cc_grad}
function eval_objective(::typeof{f}, d::NonlinearFunction{N}, df::Vector, x::Vector{Float64}) where {f <: VECTOR_ACCESS_FUNCTIONS, N}
    @assert d.has_nlobj "No nonlinear objective."
    forward_reverse_pass(d, x, 0)
    if d.numvalued[1]
        for i = 1:d.variable_number
            df[i] = f(d.setstorage[1])[i]
        end
    else
        fill!(df, 0.0)
    end
    nothing
end
function eval_objective(::typeof{f}, d::NonlinearFunction{N}, df::Vector, x::Vector{Float64}) where {f <: VECTOR_ACCESS_FUNCTIONS, N}
    forward_reverse_pass(d, x, 0)
    if d.numvalued[1]
        for i = 1:d.variable_number
            df[i] = f(d.setstorage[1])[i]
        end
    else
        fill!(df, 0.0)
    end
    nothing
end
grad_sparsity(d::NonlinearFunction) = d.grad_sparsity
=#


# looks good
function reverse_eval_all(d::Evaluator, x::Vector{Float64})
    #println("ran reverse evals")
    subexpr_values_set = d.subexpression_values_set
    subexpr_isnum = d.subexpression_isnum
    feas = true
    subgrad_tighten = d.subgrad_tighten

    if d.has_nlobj
        # Cut Objective at upper bound
        ex = d.objective
        reverse_updated_mc = ex.setstorage[1] ∩ Interval{Float64}(-Inf, d.objective_ubd)
        if isnan(reverse_updated_mc)
            reverse_updated_mc = interval_MC(reverse_updated_mc)
        end
        ex.setstorage[1] = reverse_updated_mc
        feas &= reverse_eval(ex.setstorage, ex.numberstorage, ex.numvalued, subexpr_isnum,
                             subexpr_values_set, ex.nd, ex.adj, x, d.current_node, subgrad_tighten)
    end
    for i = 1:length(d.constraints)
        # Cut constraints on constraint bounds & reverse
        if feas
            ex = d.constraints[i]
            reverse_updated_mc = ex.setstorage[1] ∩ Interval{Float64}(d.constraints_lbd[i], d.constraints_ubd[i])
            if isnan(reverse_updated_mc)
                reverse_updated_mc = interval_MC(reverse_updated_mc)
            end
            ex.setstorage[1] = reverse_updated_mc
            feas &= reverse_eval(ex.setstorage, ex.numberstorage, ex.numvalued, subexpr_isnum,
                                subexpr_values_set, ex.nd, ex.adj, x, d.current_node, subgrad_tighten)
        else
            break
        end
    end
    for k = 1:length(d.subexpression_order)
        if feas
            ex = d.subexpressions[d.subexpression_order[k]]
            ex.setstorage[1] = subexpr_values_set[d.subexpression_order[k]]
            feas &= reverse_eval(ex.setstorage, ex.numberstorage, ex.numvalued, subexpr_isnum,
                                  subexpr_values_set, ex.nd, ex.adj, x, d.current_node, subgrad_tighten)
        else
            break
        end
    end
    copyto!(d.last_x, x)

    return feas
end

"""
$(FUNCTIONNAME)

Performs a `d.cp_repetitions` number of forward passes of the set-value evaluator each followed
by a reverse pass if `d.has_reverse` as long as the node between passes differs
by more that `d.fw_atol` at each iteration.
"""
function forward_reverse_pass(d::Evaluator, x::Vector{Float64})
    flag = true
    converged_flag = false
    if d.last_x != x
        if d.has_reverse
            for i = 1:d.cp_repetitions
                d.first_eval_flag = (i == 1)
                if flag
                    forward_eval_all(d, x)
                    flag = reverse_eval_all(d, x)
                    ~flag && break
                    converged_flag = same_box(d.current_node, get_node(d), d.cp_tolerance)
                    converged_flag && break
                end
            end
            flag && forward_eval_all(d, x)
        else
            d.first_eval_flag = true
            forward_eval_all(d, x)
        end
    end

    d.last_x .= x

     return flag
end


function forward_eval_obj(d::Evaluator, x::Vector{Float64})
    subexpr_values_flt = d.subexpression_values_flt
    subexpr_values_set = d.subexpression_values_set
    user_operators = d.user_operators
    user_input_buffer = d.jac_storage
    flt_user_input_buffer = d.flt_jac_storage
    subgrad_tighten = d.subgrad_tighten
    first_eval_flag = d.first_eval_flag
    seeds = d.seeds

    seeds = d.seeds
    for (ind, k) in enumerate(reverse(d.subexpression_order))
        subex = d.subexpressions[k]
        temp = forward_eval(subex.setstorage, subex.numberstorage, subex.numvalued,
                    subex.nd, subex.adj, d.current_node,
                    x, subexpr_values_flt, subexpr_values_set, d.subexpression_isnum,
                    user_input_buffer, flt_user_input_buffer, subgrad_tighten, subex.tpdict,
                    subex.tp1storage, subex.tp2storage, subex.tp3storage,
                    subex.tp4storage, first_eval_flag, user_operators, seeds, d.ctx)
        d.subexpression_isnum[ind] = subex.numvalued[1]
        if d.subexpression_isnum[ind]
            d.subexpression_values_flt[k] = temp
        else
            d.subexpression_values_set[k] = temp
        end
    end

    if d.has_nlobj
        ex = d.objective
        forward_eval(ex.setstorage, ex.numberstorage, ex.numvalued,
                     ex.nd, ex.adj, d.current_node,
                     x, subexpr_values_flt, subexpr_values_set,
                     d.subexpression_isnum, user_input_buffer, flt_user_input_buffer, subgrad_tighten, ex.tpdict,
                     ex.tp1storage, ex.tp2storage, ex.tp3storage, ex.tp4storage,
                     first_eval_flag, user_operators, seeds, d.ctx)
    end
    return
end

"""
$(TYPEDSIGNATURES)
"""
function forward_eval_all(d::Evaluator, x::Vector{Float64})

    subexpr_values_flt = d.subexpression_values_flt
    subexpr_values_set = d.subexpression_values_set
    user_operators = d.user_operators
    user_input_buffer = d.jac_storage
    flt_user_input_buffer = d.flt_jac_storage
    subgrad_tighten = d.subgrad_tighten
    first_eval_flag = d.first_eval_flag
    seeds = d.seeds

    for (ind, k) in enumerate(reverse(d.subexpression_order))
        subex = d.subexpressions[k]
        forward_eval(subex.setstorage, subex.numberstorage, subex.numvalued,
                    subex.nd, subex.adj, d.current_node,
                    x, subexpr_values_flt, subexpr_values_set, d.subexpression_isnum,
                    user_input_buffer, flt_user_input_buffer, subgrad_tighten, subex.tpdict,
                    subex.tp1storage, subex.tp2storage, subex.tp3storage,
                    subex.tp4storage, first_eval_flag, user_operators, seeds, d.ctx)

        d.subexpression_isnum[ind] = subex.numvalued[1]
        if d.subexpression_isnum[ind]
            d.subexpression_values_flt[k] = subex.numberstorage[1]
        else
            d.subexpression_values_set[k] = subex.setstorage[1]
        end
    end

    if d.has_nlobj
        ex = d.objective
        forward_eval(ex.setstorage, ex.numberstorage, ex.numvalued,
                     ex.nd, ex.adj, d.current_node,
                     x, subexpr_values_flt, subexpr_values_set,
                     d.subexpression_isnum, user_input_buffer, flt_user_input_buffer, subgrad_tighten, ex.tpdict,
                     ex.tp1storage, ex.tp2storage, ex.tp3storage, ex.tp4storage,
                     first_eval_flag, user_operators, seeds, d.ctx)
    end

    for (ind,ex) in enumerate(d.constraints)
        #println(" ")
        #println(" ")
        #println("CONSTRAINT # = $ind")
        #println(" ")
        #println(" ")
        forward_eval(ex.setstorage, ex.numberstorage, ex.numvalued,
                     ex.nd, ex.adj, d.current_node,
                     x, subexpr_values_flt, subexpr_values_set,
                     d.subexpression_isnum, user_input_buffer, flt_user_input_buffer, subgrad_tighten, ex.tpdict,
                     ex.tp1storage, ex.tp2storage, ex.tp3storage, ex.tp4storage,
                     first_eval_flag, user_operators, seeds, d.ctx)
    end

    return
end
