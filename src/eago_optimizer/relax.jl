"""
    relax_convex_kernel

Stores the kernel of the calculation required to relax convex quadratic
constraints using the immutable dictionary to label terms.
"""
function relax_convex_kernel(func::SQF, vi::Vector{VI}, cvx_dict::ImmutableDict{Int64,Int64},
                             nx::Int64, x0::Vector{Float64})

    # initialize storage terms
    terms_coeff = zeros(nx)
    term_constant = func.constant

    # iterate over quadratic terms to build convex quadratic cut
    for term in func.quadratic_terms

        coeff = term.coefficient
        indx1 = cvx_dict[term.variable_index_1.value]
        indx2 = cvx_dict[term.variable_index_2.value]
        @inbounds temp_constant = x0[indx1]
        @inbounds temp_constant *= x0[indx2]
        term_constant -= coeff*temp_constant

        @inbounds terms_coeff[indx1] += 2.0*coeff*x0[indx1]
    end

    # Adds affine terms
    for term in func.affine_terms
        coeff = term.coefficient
        indx1 = cvx_dict[term.variable_index.value]
        @inbounds terms_coeff[indx1] += coeff
    end

    saf = SAF(SAT.(terms_coeff, vi), term_constant)

    return saf
end

function relax_nonconvex_kernel(func::SQF, vi::Vector{VI}, n::NodeBB,
                                cvx_dict::ImmutableDict{Int64,Int64}, nx::Int64, x0::Vector{Float64})

    terms_coeff = zeros(nx)

    lvb = n.lower_variable_bounds
    uvb = n.upper_variable_bounds
    quadratic_constant = func.constant
    for term in func.quadratic_terms
        a = term.coefficient
        nidx1 = term.variable_index_1.value
        nidx2 = term.variable_index_2.value
        vidx1 = cvx_dict[nidx1]
        vidx2 = cvx_dict[nidx2]
        @inbounds xL1 = lvb[nidx1]
        @inbounds xU1 = uvb[nidx2]
        @inbounds x01 = x0[nidx1]
        if nidx1 == nidx2               # quadratic terms
            if (a > 0.0)
                @inbounds terms_coeff[vidx1] = 2.0*a*x01
                quadratic_constant -= 2.0*x01*x01
            else
                @inbounds terms_coeff[vidx1] = a*(xL1 + xU1)
                quadratic_constant -= a*xL1*xU1
            end
        else
            @inbounds xL2 = lvb[nidx2]
            @inbounds xU2 = uvb[nidx2]
            @inbounds x02 = x0[nidx2]
            if (a > 0.0)
                check_ref = (xU1 - xL1)*x02 + (xU2 - xL2)*x01
                if (check_ref - xU1*xU2 + xL1*xL2) <= 0.0
                    @inbounds terms_coeff[vidx1] += a*xL2
                    @inbounds terms_coeff[vidx2] += a*xL1
                    quadratic_constant -= a*xL1*xL2
                else
                    @inbounds terms_coeff[vidx1] += a*xU2
                    @inbounds terms_coeff[vidx2] += a*xU1
                    quadratic_constant -= a*xU1*xU2
                end
            else
                check_ref = (xU1 - xL1)*x02 + (xU2 - xL2)*x01
                if check_ref <= (xU1*xL2 - xL1*xU2)
                    @inbounds terms_coeff[vidx1] += a*xL2
                    @inbounds terms_coeff[vidx2] += a*xU1
                    quadratic_constant -= a*xU1*xL2
                else
                    @inbounds terms_coeff[vidx1] += a*xU2
                    @inbounds terms_coeff[vidx2] += a*xL1
                    quadratic_constant -= a*xL1*xU2
                end
            end
        end
    end

    for term in func.affine_terms
        nidx1 = term.variable_index.value
        vidx1 = cvx_dict[nidx1]
        @inbounds terms_coeff[vidx1] += term.coefficient
    end

    saf = SAF(SAT.(terms_coeff, vi), quadratic_constant)

    return saf
end

"""
    relax_quadratic_inner!

Default routine for relaxing nonconvex quadratic constraint `lower` < `func` < `upper`
on node `n`. Takes affine bounds of convex part at point `x0` and secant line
bounds on concave parts.
"""
function relax_quadratic_inner!(eqflag::Bool, flag1::Bool, flag2::Bool,
                                func::MOI.ScalarQuadraticFunction{Float64},
                                x::Optimizer, lower::Float64, upper::Float64,
                                vi::Vector{VI}, x0::Vector{Float64}, nx::Int,
                                c1::CI{SAF, LT}, c2::CI{SAF, LT},
                                n::NodeBB, cvx_dict::ImmutableDict{Int64,Int64})


    opt = x.relaxed_optimizer
    if flag1
        saf1 = relax_convex_kernel(func, vi, cvx_dict, nx, x0)
    else
        saf1 = relax_nonconvex_kernel(func, vi, n, cvx_dict, nx, x0)
    end

    if eqflag

        for (i, term) in enumerate(saf1.terms)
            MOI.modify(opt, c1, SCoefC(term.variable_index, term.coefficient))
        end
        MOI.set(opt, MOI.ConstraintSet(), c1, LT(upper - saf1.constant))


        if flag2
            saf2 = relax_convex_kernel(func, vi, cvx_dict, nx, x0)
        else
            saf2 = relax_nonconvex_kernel(func, vi, n, cvx_dict, nx, x0)
        end

        for (i, term) in enumerate(saf2.terms)
            MOI.modify(opt, c2, SCoefC(term.variable_index, -1.0*term.coefficient))
        end
        MOI.set(opt, MOI.ConstraintSet(), c2, LT(-lower + saf2.constant))

    elseif (lower != -Inf)

        for (i, term) in enumerate(saf1.terms)
            MOI.modify(opt, c1, SCoefC(term.variable_index, -1.0*term.coefficient))
        end
        MOI.set(opt, MOI.ConstraintSet(), c1, LT(-lower + saf1.constant))

    else

        for (i, term) in enumerate(saf1.terms)
            MOI.modify(opt, c1, SCoefC(term.variable_index, term.coefficient))
        end
        MOI.set(opt, MOI.ConstraintSet(), c1, LT(upper - saf1.constant))

    end

    return
end

"""
    relax_quadratic!

Relaxes all quadratic constraints in `x` optimizer.
"""
function relax_quadratic!(x::Optimizer, x0::Vector{Float64})

    n = x._current_node

    # Relax Convex Constraint Terms TODO: place all quadratic info into one vector of tuples?
    for i in 1:length(x._quadratic_leq_constraints)
        @inbounds func, set, j = x._quadratic_leq_constraints[i]
        @inbounds cvx_dict = x._quadratic_leq_dict[i]
        @inbounds vi = x._quadratic_leq_sparsity[i]
        @inbounds nz = x._quadratic_leq_gradnz[i]
        @inbounds ci1 = x._quadratic_ci_leq[i]
        @inbounds flag1 = x._quadratic_leq_convexity[i]
        relax_quadratic_inner!(false, flag1, false, func, x, -Inf, set.upper, vi, x0,
                               nz, ci1, ci1,  n, cvx_dict)
    end

    for i in 1:length(x._quadratic_geq_constraints)
        @inbounds func, set, j = x._quadratic_geq_constraints[i]
        @inbounds cvx_dict = x._quadratic_geq_dict[i]
        @inbounds vi = x._quadratic_geq_sparsity[i]
        @inbounds nz = x._quadratic_geq_gradnz[i]
        @inbounds ci1 = x._quadratic_ci_geq[i]
        @inbounds flag1 = x._quadratic_geq_convexity[i]
        relax_quadratic_inner!(false, flag1, false, func, x, set.lower, Inf, vi, x0,
                               nz, ci1, ci1,  n, cvx_dict)
    end

    for i in 1:length(x._quadratic_eq_constraints)
        @inbounds func, set, j = x._quadratic_eq_constraints[i]
        @inbounds x._quadratic_eq_convexity_1[i] = is_convex_quadratic(func, 1.0)
        @inbounds x._quadratic_eq_convexity_2[i] = is_convex_quadratic(func, -1.0)
        @inbounds cvx_dict = x._quadratic_eq_dict[i]
        @inbounds vi = x._quadratic_eq_sparsity[i]
        @inbounds nz = x._quadratic_eq_gradnz[i]
        @inbounds ci1, ci2 = x._quadratic_ci_eq[i]
        @inbounds flag1 = x._quadratic_eq_convexity_1[i]
        @inbounds flag2 = x._quadratic_eq_convexity_2[i]
        relax_quadratic_inner!(true, flag1, flag2, func, x, set.value, set.value, vi, x0,
                               nz, ci1, ci2, n, cvx_dict)
    end

    return
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
                x_val = x._current_xref

                f_grad = zeros(x._variable_number)
                MOI.eval_objective_gradient(x._relaxed_evaluator, f_grad, x_val)

                c_term = 0.0
                for term in saf.terms
                    term_coeff = term.coefficient
                    term_indx = term.variable_index
                    MOI.modify(x.relaxed_optimizer, ci_saf, MOI.SCC(term_indx, term_coeff))
                    @inbounds c_term += x_val[term_indx]*f_grad[term_indx]
                end
                fstar = MOI.eval_objective(x._relaxed_evaluator, x_val)
                set = LT(x._global_upper_bound - fstar + cterm)
                MOI.set(x.relaxed_optimizer, MOI.ConstraintSet(), ci_saf, set)
            else
                ci_sv = x._objective_cut_ci_sv
                set = LT(x._global_upper_bound)
                MOI.set(x.relaxed_optimizer, MOI.ConstraintSet(), ci_sv, set)
            end
        end
    end
    return
end

relax_problem!(x::Optimizer, v::Vector{Float64}) = relax_problem!(x.ext_type, x, v)
relax_objective!(x::Optimizer, v::Vector{Float64}) = relax_objective!(x.ext_type, x, v)
