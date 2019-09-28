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
    opt = x.relaxed_optimizer
    @inbounds vi = x._lower_variable_index[1:nx]

    # Add objective
    if x._objective_is_sv

        MOI.set(opt, MOI.ObjectiveFunction{SV}(), x._objective_sv)
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    elseif x._objective_is_saf

        MOI.set(opt, MOI.ObjectiveFunction{SAF}(), x._objective_saf)
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    elseif x._objective_is_sqf

        if (x._objective_convexity)
            saf = relax_convex_kernel(x._objective, vi, nx, x0)
        else
            saf = relax_nonconvex_kernel(x._objective, vi, x,  nx, x0)
        end

        MOI.set(opt, MOI.ObjectiveFunction{SAF}(), saf)
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    else
        if x._objective_is_nlp

            evaluator = x._relaxed_evaluator

            # Calculates convex relaxation
            f = MOI.eval_objective(evaluator, x0)

            # calculates the convex relaxation subgradient
            df = zeros(nx)
            MOI.eval_objective_gradient(evaluator, df, x0)

            # Add objective relaxation to model
            saf_const = f
            grad_c = 0.0
            x0_c = 0.0
            for i in 1:nx
                @inbounds grad_c = df[i]
                @inbounds x0_c = x0[i]
                @inbounds vindx = vi[i]
                saf_const -= x0_c*grad_c
                MOI.modify(opt,  MOI.ObjectiveFunction{SAF}(), SCoefC(vindx, grad_c))
            end
            MOI.modify(opt,  MOI.ObjectiveFunction{SAF}(), SConsC(saf_const))
            MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        end
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


    opt = x.relaxed_optimizer
    if flag
        saf = relax_convex_kernel(func, vi, nx, x0)
    else
        saf = relax_nonconvex_kernel(func, vi, x, nx, x0)
    end

    if (lower == upper)

        for (i, term) in enumerate(saf.terms)
            MOI.modify(opt, c1, MOI.SCC(term.variable_index, term.coefficent))
        end
        MOI.set(opt, MOI.ConstraintSet(), c1, LT(upper - saf.constant))

        saf1 = copy(saf)
        for (i, term) in enumerate(saf1.terms)
            MOI.modify(opt, c2, MOI.SCC(term.variable_index, -1.0*term.coefficent))
        end
        MOI.set(opt, MOI.ConstraintSet(), c2, LT(-lower + saf1.constant))

    elseif (lower != -Inf)

        for (i, term) in enumerate(saf.terms)
            MOI.modify(opt, c1, MOI.SCC(term.variable_index, -1.0*term.coefficent))
        end
        MOI.set(opt, MOI.ConstraintSet(), c1, LT(-lower + saf.constant))

    else

        for (i, term) in enumerate(saf.terms)
            MOI.modify(opt, c1, MOI.SCC(term.variable_index, term.coefficent))
        end
        MOI.set(opt, MOI.ConstraintSet(), c1, LT(upper - saf.constant))

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

                lower_nlp_sparsity = x._lower_nlp_sparsity
                lower_nlp_affine = x._lower_nlp_affine
                lower_nlp_affine_indx = x._lower_nlp_affine_indx
                for i in 1:length(lower_nlp_affine_indx)
                    @inbounds g_indx = lower_nlp_affine_indx[i]
                    @inbounds aff_ci = lower_nlp_affine[i]
                    @inbounds nzidx = lower_nlp_sparsity[i]
                    @inbounds nzvar = vi[nzidx]
                    @inbounds constant = g[g_indx]
                    dg_cv_val = 0.0
                    for j in nzidx
                        @inbounds dg_cv_val = dg[i,j]
                        @inbounds vindx = vi[i]
                        @inbounds constant -= v[i]*dg_cv_val
                        MOI.modify(x.relaxed_optimizer, aff_ci, SCoefC(vindx, dg_cv_val))
                    end
                    @inbounds bns = constraint_bounds[i]
                    set = LT(bns.upper-constant)
                    MOI.set(x.relaxed_optimizer, MOI.ConstraintSet(), aff_ci, set)
                end

                upper_nlp_sparsity = x._upper_nlp_sparsity
                upper_nlp_affine = x._upper_nlp_affine
                upper_nlp_affine_indx = x._upper_nlp_affine_indx
                for i in 1:length(upper_nlp_affine_indx)
                    @inbounds g_indx = upper_nlp_affine_indx[i]
                    @inbounds aff_ci = upper_nlp_affine[i]
                    @inbounds nzidx = upper_nlp_sparsity[i]
                    @inbounds nzvar = vi[nzidx]
                    @inbounds constant = g_cc[g_indx]
                    dg_cc_val = 0.0
                    for j in nzidx
                        @inbounds dg_cc_val = dg[i,j]
                        @inbounds vindx = vi[i]
                        @inbounds constant -= v[i]*dg_cc_val
                        MOI.modify(x.relaxed_optimizer, aff_ci, SCoefC(vindx, dg_cc_val))
                    end
                    @inbounds bns = constraint_bounds[i]
                    set = LT(constant - bns.lower)
                    MOI.set(x.relaxed_optimizer, MOI.ConstraintSet(), aff_ci, set)
                end
            end
        end
    end
    return
end

function relax_problem!(t::ExtensionType, x::Optimizer, v::Vector{Float64})

    evaluator = x._relaxed_evaluator
    set_current_node!(evaluator, x._current_node)
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
        if ~x._objective_cut_set
            set = LT(x._global_upper_bound)
            if x._objective_is_sv
                ci_sv = MOI.add_constraint(x_relaxed_optimizer, x._objective_sv, set)
                x._objective_cut_ci_sv = ci
            elseif x._objective_is_saf
                ci_saf = MOI.add_constraint(x.relaxed_optimizer, x._objective_saf, set)
                x._objective_cut_ci_saf = ci_saf
            elseif x._objective_is_sqf || x._objective_is_nlp
                saf = MOI.get(x.relaxed_optimizer, MOI.ObjectiveFunction{SAF}())
                set = LT(x._global_upper_bound - saf.constant)
                saf.constant = 0.0
                ci_saf = MOI.add_constraint(x.relaxed_optimizer, saf, set)
                x._objective_cut_ci_saf = ci_saf
            end
        else
            if ~x._objective_is_sv
                ci_saf = x._objective_cut_ci_saf
                saf = MOI.get(x.relaxed_optimizer, MOI.ObjectiveFunction{SAF}())
                for term in saf.terms
                    MOI.modify(x.relaxed_optimizer, ci_saf, MOI.SCC(term.variable_index, term.coefficient))
                end
                MOI.set(x.relaxed_optimizer, MOI.ConstraintSet(), ci_saf, set)
            else
                ci_sv = x._objective_cut_ci_sv
                MOI.set(x.relaxed_optimizer, MOI.ConstraintSet(), ci_sv, set)
            end
        end
    end
    return
end
