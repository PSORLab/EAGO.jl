"""
    relax_convex_kernel

Stores the kernel of the calculation required to relax convex quadratic constraints.
"""
function relax_convex_kernel(func::SQF, vi::Vector{VI}, nx::Int64, x0::Vector{Float64})

    _variable_number = x._variable_numer
    terms_coeff = zeros(_variable_number)

    term_constant = func.constant
    for term in func.quadratic_terms
        @inbounds temp_constant = x0[term.variable_index_1.value]
        @inbounds temp_constant *= x0[term.variable_index_2.value]
        term_constant -= term.coefficient*temp_constant
    end

    # Adds affine terms
    for term in func.affine_terms
        @inbounds terms_coeff[term.variable_index.value] += term.coefficient
    end

    # Adds quadratic terms
    for term in func.quadratic_terms
        @inbounds terms_coeff[term.variable_index_1.value] += 2.0*term.coefficient*x0[term.variable_index_1.value]
    end

    saf = SAF(SAT.(terms_coeff, vi), term_constant)

    return saf
end

function relax_nonconvex_kernel(func::SQF, vi::Vector{VI},
                                x::Optimizer, nx::Int64, x0::Vector{Float64})
    n = x._current_node
    quadratic_coefficient = zeros(nx)

    quadratic_constant = func.constant
    for term in func.quadratic_terms
        a = term.coefficient
        NodeIndx1 = term.variable_index_1.value
        NodeIndx2 = term.variable_index_2.value
        @inbounds xL1 = n.lower_variable_bounds[NodeIndx1]
        @inbounds xU1 = n.upper_variable_bounds[NodeIndx1]
        @inbounds x01 = x0[NodeIndx1]
        if NodeIndx1 == NodeIndx2               # quadratic terms
            if (a > 0.0)
                @inbounds quadratic_coefficient[NodeIndx1] = 2.0*a*x01
                quadratic_constant -= 2.0*x01*x01
            else
                @inbounds quadratic_coefficient[NodeIndx1] = a*(xL1 + xU1)
                quadratic_constant -= a*xL1*xU1
            end
        else
            @inbounds xL2 = n.lower_variable_bounds[NodeIndx2]
            @inbounds xU2 = n.upper_variable_bounds[NodeIndx2]
            if (a > 0.0)
                @inbounds check_ref = (xU1 - xL1)*x0[NodeIndx2] + (xU2 - xL2)*x0[NodeIndx1]
                if (check_ref - xU1*xU2 + xL1*xL2) <= 0.0
                    @inbounds quadratic_coefficient[NodeIndx1] += a*xL2
                    @inbounds quadratic_coefficient[NodeIndx2] += a*xL1
                    quadratic_constant -= a*xL1*xL2
                else
                    @inbounds quadratic_coefficient[NodeIndx1] += a*xU2
                    @inbounds quadratic_coefficient[NodeIndx2] += a*xU1
                    quadratic_constant -= a*xU1*xU2
                end
            else
                @inbounds check_ref = (xU1 - xL1)*x0[NodeIndx2] + (xU2 - xL2)*x0[NodeIndx1]
                if check_ref <= (xU1*xL2 - xL1*xU2)
                    @inbounds quadratic_coefficient[NodeIndx1] += a*xL2
                    @inbounds quadratic_coefficient[NodeIndx2] += a*xU1
                    quadratic_constant -= a*xU1*xL2
                else
                    @inbounds quadratic_coefficient[NodeIndx1] += a*xU2
                    @inbounds quadratic_coefficient[NodeIndx2] += a*xL1
                    quadratic_constant -= a*xL1*xU2
                end
            end
        end
    end

    for term in func.affine_terms
        NodeIndx3 = term.variable_index.value
        @inbounds quadratic_coefficient[NodeIndx3] += term.coefficient
    end

    saf = SAF(SAT.(quadratic_coefficient, vi), quadratic_constant)

    return saf
end

"""
    relax_objective!

A rountine that only relaxes the objective.
"""
function relax_objective!(t::ExtensionType, x::Optimizer, x0::Vector{Float64})

    nx = x._variable_number
    @inbounds vi = x._lower_variable_index[1:nx]

    # Add objective
    if isa(x._objective, SV)

        MOI.set(x._relaxed_optimizer, MOI.ObjectiveFunction{SV}(), x._objective)
        MOI.set(x._relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    elseif isa(x._objective, SAF)

        MOI.set(x._relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), x._objective)
        MOI.set(x._relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    elseif isa(x._objective, SQF)

        if (x._objective_convexity)
            saf = relax_convex_kernel(x._objective, vi, nx, x0)
        else
            saf = relax_nonconvex_kernel(x._objective, vi, x,  nx, x0)
        end

        MOI.set(x._relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), saf)
        MOI.set(x._relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    elseif eval_block.has_objective

        eval_block = x._working_evaluator_block

        # Calculates convex relaxation
        f = MOI.eval_objective(eval_block.evaluator, x0)

        # calculates the convex relaxation subgradient
        df = zeros(nx)
        MOI.eval_objective_gradient(eval_block.evaluator, df, x0)

        # Add objective relaxation to model
        saf_const = f
        grad_c = 0.0
        x0_c = 0.0
        @simd for i in 1:nx
            @inbounds grad_c = df[i]
            @inbounds x0_c = x0[i]
            @inbounds vindx = vi[i]
            saf_const -= xpnt_c*grad_c
            MOI.modify(x._relaxed_optimizer,  MOI.ObjectiveFunction{SAF}(), MOI.SCoef(vindx, grad_c))
        end
        MOI.modify(x._relaxed_optimizer,  MOI.ObjectiveFunction{SAF}(), MOI.SConsC(saf_const))
        MOI.set(x._relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    end
    return
end

"""
    relax_quadratic_inner!

Default routine for relaxing nonconvex quadratic constraint `lower` < `func` < `upper`
on node `n`. Takes affine bounds of convex part at point `x0` and secant line
bounds on concave parts.
"""
function relax_quadratic_inner!(flag::Bool, func::MOI.ScalarQuadraticFunction{Float64},
                                x::Optimizer, lower::Float64, upper::Float64,
                                vi::Vector{VI}, x0::Vector{Float64}, nx::Int,
                                c1::CI{SAF, LT},
                                c2::CI{SAF, LT})


    if flag
        saf = relax_convex_kernel(func, vi, nx, x0)
    else
        saf = relax_nonconvex_kernel(func, vi, x, nx, x0)
    end

    if (lower == upper)

        for (i, term) in enumerate(saf.terms)
            MOI.modify(x._relaxed_optimizer, c1, MOI.SCC(term.variable_index, term.coefficent))
        end
        MOI.set(x._relaxed_optimizer, MOI.ConstraintSet(), c1, LT(upper - saf.constant))

        saf1 = copy(saf)
        for (i, term) in enumerate(saf1.terms)
            MOI.modify(x._relaxed_optimizer, c2, MOI.SCC(term.variable_index, -1.0*term.coefficent))
        end
        MOI.set(x._relaxed_optimizer, MOI.ConstraintSet(), c2, LT(-lower + saf1.constant))

    elseif (lower != -Inf)

        for (i, term) in enumerate(saf.terms)
            MOI.modify(x._relaxed_optimizer, c1, MOI.SCC(term.variable_index, -1.0*term.coefficent))
        end
        MOI.set(x._relaxed_optimizer, MOI.ConstraintSet(), c1, LT(-lower + saf.constant))

    else

        for (i, term) in enumerate(saf.terms)
            MOI.modify(x._relaxed_optimizer, c1, MOI.SCC(term.variable_index, term.coefficent))
        end
        MOI.set(x._relaxed_optimizer, MOI.ConstraintSet(), c1, LT(upper - saf.constant))

    end

    return
end

"""
    relax_quadratic!

Relaxes all quadratic constraints in `x` optimizer.
"""
function relax_quadratic!(x::Optimizer, x0::Vector{Float64})

    nx = x._variable_number
    vi = x._lower_variable_index

    # Relax Convex Constraint Terms
    for (func, set, ci1, ci2, i) in x._quadratic_leq_constraints
        flag = x._constraint_convexity[i]
        relax_quadratic_inner!(flag, func, x, -Inf, set.upper, n, x0, ci1, ci2)
    end
    for (func, set, ci1, ci2, i) in x._quadratic_geq_constraints
        flag = x._constraint_convexity[i]
        relax_quadratic_inner!(flag, func, x, set.lower, Inf, vi, x0, nx, ci1, ci2)
    end
    for (func, set, ci1, ci2, i) in x._quadratic_eq_constraints
        flag =  x._constraint_convexity[i]
        relax_quadratic_inner!(flag, func, x, set.value, set.valu, vi, x0, nx, ci1, ci2)
    end
    return
end

function relax_nlp!(x::Optimizer, v::Vector{Float64})

    eval_block = x._working_evaluator_block

    if ~isempty(x._branch_variable)

        if MOI.supports(x._relaxed_optimizer, MOI.NLPBlock())

            _nlp_data = MOI.NLPBlockData(x._nlp_data.constraint_bounds,
                                         eval_block.evaluator,
                                         x._nlp_data.has_objective)
            MOI.set(x._relaxed_optimizer, MOI.NLPBlock(), _nlp_data)

        else

            # Add other affine constraints
            if length(eval_block.constraint_bounds) > 0

                nx = x._variable_number
                vi = x._lower_variable_index

                leng = length(eval_block.constraint_bounds)
                g = zeros(leng)
                dg = zeros(leng, nx)

                g_cc = zeros(leng)
                dg_cc = zeros(leng, nx)

                MOI.eval_constraint(eval_block.evaluator, g, v)
                MOI.eval_constraint_jacobian(eval_block.evaluator, dg, v)

                eval_constraint_cc(eval_block.evaluator, g_cc, v)
                eval_constraint_cc_grad(eval_block.evaluator, dg_cc, v)

                count_lower = 0
                count_upper = 0
                for (j, bns) in enumerate(eval_block.constraint_bounds)
                    @inbounds nzidx = eval.constraints[j].grad_sparsity
                    @inbounds nzvar = vi[nzidx]
                    if bns.upper != Inf
                        @inbounds cu = x._upper_nlp_affine[count_upper]
                        @inbounds constant = g[j]
                        dg_cv_val = 0.0
                        for i in nzidx
                            @inbounds dg_cv_val = dg[j,i]
                            @inbounds vindx = vi[i]
                            constant -= vindx*dg_cv_val
                            MOI.modify(x._relaxed_optimizer, cu, MOI.SCC(vindx, dg_cv_val))
                        end
                        set = LT(bns.upper-constant)
                        MOI.set(x._relaxed_optimizer, MOI.ConstraintSet(), cu, set)
                        count_upper += 1
                    end
                    if bns.lower > -Inf
                        @inbounds cl = x._upper_nlp_affine[count_lower]
                        @inbounds constant = g_cc[j]
                        dg_cc_val = 0.0
                        for i in nzidx
                            @inbounds dg_cc_val = dg_cc[j,i]
                            @inbounds vindx = vi[i]
                            constant -= vindx*dg_cc_val
                            MOI.modify(x._relaxed_optimizer, cl, MOI.SCC(vindx, -dg_cc_val))
                        end
                        set = LT(constant)
                        MOI.set(x._relaxed_optimizer, MOI.ConstraintSet(), cl, set)
                        count_lower += 1
                    end
                end
            end
        end
    end
    return
end

function relax_problem!(t::ExtensionType, x::Optimizer, v::Vector{Float64})

    eval_block = x._working_evaluator_block
    set_current_node!(eval_block.evaluator, x._current_node)
    relax_quadratic!(x, v)
    relax_nlp!(x,v)

    return
end

"""
    objective_cut_linear!

Adds linear objective cut constraint to `trg` optimizer.
"""
function objective_cut_linear!(x::Optimizer)
    if x._global_upper_bound < Inf
        if x._objective_cut_ci === nothing
            set = LT(x._global_upper_bound)
            if isa(x._objective, SV)
                ci = MOI.add_constraint(x._relaxed_optimizer, x._objective, set)
            elseif isa(x._objective, SAF)
                ci = MOI.add_constraint(x._relaxed_optimizer, x._objective, set)
            elseif isa(x._objective, SQF) || x._working_evaluator_block.evaluator.has_nlobj
                saf = MOI.get(x._relaxed_optimizer, MOI.ObjectiveFunction{SAF}())
                ci = MOI.add_constraint(x._relaxed_optimizer, saf, set)
            end
            x._objective_cut_ci = ci
        else
            ci = x._objective_cut_ci
            if ~isa(x._objective, SV)
                saf = MOI.get(x._relaxed_optimizer, MOI.ObjectiveFunction{SAF}())
                for term in saf.terms
                    MOI.modify(x._relaxed_optimizer, ci, MOI.SCC(term.variable_index, term.coefficient))
                end
            end
            MOI.set(x._relaxed_optimizer, MOI.ConstraintSet(), ci, set)
        end
    end
    return
end
