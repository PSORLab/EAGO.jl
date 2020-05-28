function build_nlp_kernel!(d::Evaluator{N,T}, src::JuMP.NLPEvaluator, x::Optimizer, bool_flag::Bool) where {N,T<:RelaxTag}

    m = src.m::Model
    num_variables_ = JuMP.num_variables(m)
    d.variable_number = num_variables_
    #nldata = m.nlp_data::JuMP._NLPData
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
