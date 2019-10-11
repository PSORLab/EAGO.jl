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
        @inbounds xU1 = uvb[nidx1]
        @inbounds x01 = x0[nidx1]
        if nidx1 == nidx2               # quadratic terms
            if (a > 0.0)
                @inbounds terms_coeff[vidx1] = 2.0*a*x01
                quadratic_constant -= x01*x01
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
                check_ref = (xU1 - xL1)*x02 - (xU2 - xL2)*x01
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
function relax_quadratic_gen_saf(func::SQF, vi::Vector{VI}, n::NodeBB,
                                 nx::Int64, x0::Vector{Float64},
                                 cvx_dict::ImmutableDict{Int64,Int64}, flag::Bool)
    if flag
        saf = relax_convex_kernel(func, vi, cvx_dict, nx, x0)
    else
        saf = relax_nonconvex_kernel(func, vi, n, cvx_dict, nx, x0)
    end
    return saf
end
function store_ge_quadratic!(x::Optimizer, ci::CI{SAF,LT}, saf::SAF,
                             lower::Float64, i::Int64, q::Int64)
    opt = x.relaxed_optimizer
    if (q == 1) & x.relaxed_inplace_mod
        for (i, term) in enumerate(saf.terms)
            MOI.modify(opt, ci, SCoefC(term.variable_index, -1.0*term.coefficient))
        end
        MOI.set(opt, MOI.ConstraintSet(), ci, LT(-lower + saf.constant))
    else
        c = saf.constant
        saf.constant = 0.0
        x._quadratic_ci_geq[q][i] = MOI.add_constraint(opt, saf, LT(-lower + c))
    end
    return
end
function store_le_quadratic!(x::Optimizer, ci::CI{SAF,LT}, saf::SAF,
                            upper::Float64, i::Int64, q::Int64)
    opt = x.relaxed_optimizer
    if (q == 1) & x.relaxed_inplace_mod
        opt = x.relaxed_optimizer
        for (i, term) in enumerate(saf.terms)
            MOI.modify(opt, ci, SCoefC(term.variable_index, term.coefficient))
        end
        MOI.set(opt, MOI.ConstraintSet(), ci, LT(upper - saf.constant))

    else
        c = saf.constant
        saf.constant = 0.0
        #println("c: $c")
        new_set = LT(upper - c)
        #println("c: $c")
        cindxo = MOI.add_constraint(opt, saf, new_set)
        x._quadratic_ci_leq[q][i] = cindxo
    end
    return
end
function store_eq_quadratic!(x::Optimizer, ci1::CI{SAF,LT}, ci2::CI{SAF,LT},
                            saf1::SAF, saf2::SAF, value::Float64, i::Int64,
                            q::Int64)
    opt = x.relaxed_optimizer
    if (q == 1) & x.relaxed_inplace_mod
        store_ge_quadratic!(x, ci1, saf1, value, q)
        store_le_quadratic!(x, ci2, saf2, value, q)
    else
        c1 = saf1.constant
        c2 = saf2.constant
        saf1.constant = 0.0
        saf2.constant = 0.0
        c1 = MOI.add_constraint(opt, saf1, LT(value - c1))
        c2 = MOI.add_constraint(opt, saf2, LT(-value + c2))
        x._quadratic_ci_eq[q][i] = (c1,c2)
    end
    return
end

"""
    relax_quadratic!

Relaxes all quadratic constraints in `x` optimizer.
"""
function relax_quadratic!(x::Optimizer, x0::Vector{Float64}, q::Int64)

    n = x._current_node

    # Relax Convex Constraint Terms TODO: place all quadratic info into one vector of tuples?
    for i in 1:length(x._quadratic_leq_constraints)
        func, set, j = x._quadratic_leq_constraints[i]
        cvx_dict = x._quadratic_leq_dict[i]
        vi = x._quadratic_leq_sparsity[i]
        nz = x._quadratic_leq_gradnz[i]
        ci1 = x._quadratic_ci_leq[q][i]
        flag1 = x._quadratic_leq_convexity[i]

        saf = relax_quadratic_gen_saf(func, vi, n, nz, x0, cvx_dict, flag1)
        store_le_quadratic!(x, ci1, saf, set.upper, i, q)
    end

    for i in 1:length(x._quadratic_geq_constraints)
        func, set, j = x._quadratic_geq_constraints[i]
        cvx_dict = x._quadratic_geq_dict[i]
        vi = x._quadratic_geq_sparsity[i]
        nz = x._quadratic_geq_gradnz[i]
        ci1 = x._quadratic_ci_geq[q][i]
        flag1 = x._quadratic_geq_convexity[i]
        vec_sat = SAT[SAT(-t.coefficient, t.variable_index) for t in func.affine_terms]
        vec_sqt = SQT[SQT(-t.coefficient, t.variable_index_1, t.variable_index_2) for t in func.quadratic_terms]
        func_minus = SQF(vec_sat, vec_sqt, -func.constant)
        saf = relax_quadratic_gen_saf(func_minus, vi, n, nz, x0, cvx_dict, flag1)
        store_ge_quadratic!(x, ci1, saf, set.lower, i, q)
    end

    for i in 1:length(x._quadratic_eq_constraints)
        func, set, j = x._quadratic_eq_constraints[i]
        cvx_dict = x._quadratic_eq_dict[i]
        vi = x._quadratic_eq_sparsity[i]
        nz = x._quadratic_eq_gradnz[i]
        ci1, ci2 = x._quadratic_ci_eq[q][i]
        flag1 = x._quadratic_eq_convexity_1[i]
        flag2 = x._quadratic_eq_convexity_2[i]
        saf1 = relax_quadratic_gen_saf(func, vi, n, nz, x0, cvx_dict, flag1)
        vec_sat = SAT[SAT(-t.coefficient, t.variable_index) for t in func.affine_terms]
        vec_sqt = SQT[SQT(-t.coefficient, t.variable_index_1, t.variable_index_2) for t in func.quadratic_terms]
        func_minus = SQF(vec_sat, vec_sqt, -func.constant)
        saf2 = relax_quadratic_gen_saf(func_minus, vi, n, nz, x0, cvx_dict, flag2)
        store_eq_quadratic!(x, ci1, ci2, saf1, saf2, set.value, i, q)
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
                    for i in 1:length(lower_nlp_affine_indx)
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
                    for i in 1:length(upper_nlp_affine_indx)
                        @inbounds g_indx = upper_nlp_affine_indx[i]
                        @inbounds aff_ci = upper_nlp_affine[i]
                        @inbounds nzidx = upper_nlp_sparsity[i]
                        @inbounds nzvar = vi[nzidx]
                        @inbounds constant = g_cc[g_indx]
                        dg_cc_val = 0.0
                        for j in nzidx
                            @inbounds dg_cc_val = dg[i,j]
                            @inbounds vindx = vi[j]
                            @inbounds constant -= v[j]*dg_cc_val
                            MOI.modify(x.relaxed_optimizer, aff_ci, SCoefC(vindx, dg_cc_val))
                        end
                        @inbounds bns = constraint_bounds[i]
                        set = LT(constant - bns.lower)
                        MOI.set(x.relaxed_optimizer, MOI.ConstraintSet(), aff_ci, set)
                    end
                else
                    for i in 1:length(lower_nlp_affine_indx)
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
                    for i in 1:length(upper_nlp_affine_indx)
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
                            @inbounds dg_cc_val = dg[i,indx]
                            @inbounds coeff[j] = dg_cc_val
                            @inbounds constant -= v[indx]*dg_cc_val
                        end
                        @inbounds bns = constraint_bounds[i]
                        set = LT(constant - bns.lower)
                        saf = SAF(SAT.(coeff,vindices), 0.0)
                        x._upper_nlp_affine[q][i] = MOI.add_constraint(x.relaxed_optimizer,
                                                                   saf, set)
                    end
                end
            end
        end
    end
    return
end

function relax_problem!(t::ExtensionType, x::Optimizer, v::Vector{Float64}, q::Int64)

    evaluator = x._relaxed_evaluator
    if q == 1
        set_current_node!(evaluator, x._current_node)
    end
    relax_quadratic!(x, v, q)
    relax_nlp!(x, v, q)

    return
end

"""
    objective_cut_linear!

Adds linear objective cut constraint to `trg` optimizer.
"""
function objective_cut_linear!(x::Optimizer, q::Int64)
    if x._global_upper_bound < Inf
        if (x._objective_cut_set == -1) || (q > 1) ||  ~x.relaxed_inplace_mod
            z = (x._objective_cut_set == -1) ? 1 : q
            set = LT(x._global_upper_bound)
            if x._objective_is_sv
                ci_sv = x._objective_cut_ci_sv
                MOI.set(x.relaxed_optimizer, MOI.ConstraintSet(), ci_sv, set)
            elseif x._objective_is_saf
                ci_saf = MOI.add_constraint(x.relaxed_optimizer, x._objective_saf, set)
                x._objective_cut_ci_saf[z] = ci_saf
            elseif x._objective_is_sqf || x._objective_is_nlp
                saf = MOI.get(x.relaxed_optimizer, MOI.ObjectiveFunction{SAF}())
                set = LT(x._global_upper_bound - saf.constant)
                saf.constant = 0.0
                ci_saf = MOI.add_constraint(x.relaxed_optimizer, saf, set)
                x._objective_cut_ci_saf[z] = ci_saf
            end
            x._objective_cut_set = x._iteration_count
        elseif q == 1
            if ~x._objective_is_sv
                ci_saf = x._objective_cut_ci_saf[1]
                saf = MOI.get(x.relaxed_optimizer, MOI.ObjectiveFunction{SAF}())
                for term in saf.terms
                    term_coeff = term.coefficient
                    term_indx = term.variable_index
                    MOI.modify(x.relaxed_optimizer, ci_saf, SCoefC(term_indx, term_coeff))
                end
                #MOI.modify(x.relaxed_optimizer, ci_saf, SConsC(0.0))
                set = LT(x._global_upper_bound - saf.constant)
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

relax_problem!(x::Optimizer, v::Vector{Float64}, q::Int64) = relax_problem!(x.ext_type, x, v, q)
relax_objective!(x::Optimizer, v::Vector{Float64}) = relax_objective!(x.ext_type, x, v)
