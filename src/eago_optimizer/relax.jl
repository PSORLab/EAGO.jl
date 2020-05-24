"""
$(FUNCTIONNAME)

Default routine for relaxing quadratic constraint `lower` < `func` < `upper`
on node `n`. Takes affine bounds of convex part at point `x0` and secant line
bounds on concave parts.
"""
function affine_relax_quadratic!(func::SQF, buffer::OrderedDict{Int,Float64}, saf::SAF, n::NodeBB, x::Vector{Float64}, p1::Bool)

    lower_bounds = n.lower_variable_bounds
    upper_bounds = n.upper_variable_bounds

    quadratic_constant = func.constant

    for term in func.quadratic_terms

        a = p1 ? term.coefficient : -term.coefficient
        idx1 = term.variable_index_1.value
        idx2 = term.variable_index_2.value

        x0_1 = @inbounds x0[idx1]
        xL_1 = @inbounds lower_bounds[idx1]
        xU_1 = @inbounds upper_bounds[idx1]

        if idx1 === idx2

            @inbounds buffer[idx1] = (a > 0.0) ? 2.0*a*x0_1 : a*(xL_1 + xU_1)
            quadratic_constant -= (a > 0.0) ? x0_1*x0_1 : a*xL_1*xU_1

        else
            x0_2 = @inbounds x0[idx2]
            xL_2 = @inbounds lower_bounds[idx2]
            xU_2 = @inbounds upper_bounds[idx2]

            if a > 0.0
                check_ref = (xU_1 - xL_1)*x0_2 + (xU_2 - xL_2)*x0_1
                if check_ref <= xU_1*xU_2 - xL_1*xL_2
                    @inbounds buffer[idx1] += a*xL_2
                    @inbounds buffer[idx2] += a*xL_1
                    quadratic_constant -= a*xL_1*xL_2
                else
                    @inbounds buffer[idx1] += a*xU_2
                    @inbounds buffer[idx2] += a*xU_1
                    quadratic_constant -= a*xU_1*xU_2
                end
            else
                check_ref = (xU_1 - xL_1)*x0_2 - (xU_2 - xL_2)*x0_1
                if check_ref <= xU_1*xL_2 - xL_1*xU_2
                    @inbounds buffer[idx1] += a*xL_2
                    @inbounds buffer[idx2] += a*xU_1
                    quadratic_constant -= a*xU_1*xL_2
                else
                    @inbounds buffer[idx1] += a*xU_2
                    @inbounds buffer[idx2] += a*xL_1
                    quadratic_constant -= a*xL_1*xU_2
                end
            end
        end
    end

    for term in func.affine_terms
        a = p1 ? term.coefficient : -term.coefficient
        idx = term.variable_index.value
        @inbounds buffer[idx] += a
    end

    count = 1
    for (key, value) in buffer
        @inbounds saf.terms[count] = SAT(value, VI(key))
        count += 1
    end
    saf.constant = quadratic_constant

    return
end

affine_relax!(f::BufferedQuadratic{LT}, n::NodeBB, x::Vector{Float64}) = affine_relax_quadratic!(f.func, f.buffer, f.saf1, n, x, true)
affine_relax!(f::BufferedQuadratic{GT}, n::NodeBB, x::Vector{Float64}) = affine_relax_quadratic!(f.func, f.buffer, f.saf1, n, x, false)
function affine_relax!(f::BufferedQuadratic{ET}, n::NodeBB, x::Vector{Float64})
    affine_relax_quadratic!(f.func, f.buffer, f.saf1, n, x, true)
    affine_relax_quadratic!(f.func, f.buffer, f.saf2, n, x, false)
    nothing
end

function relax!(m::Optimizer, f::BufferedQuadratic{LT}, is_first::Bool)
    affine_relax!(f, m._current_node, x)
    if is_first && x.relaxed_inplace_mod
        for term in saf.terms
            MOI.modify(opt, ci, SCoefC(term.variable_index, term.coefficient))
        end
        MOI.set(opt, MOI.ConstraintSet(), ci, LT(set.upper - f.saf1.constant))
    else
        x._quadratic_ci_leq[q][i] = MOI.add_constraint(opt, f.saf1, LT(set.upper  - f.saf1.constant))
    end
    nothing
end

function relax!(m::Optimizer, f::BufferedQuadratic{ET}, k::Int)
    opt = x.relaxed_optimizer
    if is_first && x.relaxed_inplace_mod
        store_ge_quadratic!(x, ci1, saf1, value, i, q)
        store_le_quadratic!(x, ci2, saf2, value, i, q)
    else
        saf1.constant = 0.0
        saf2.constant = 0.0
        cindx1 = MOI.add_constraint(opt, saf1, LT(value))
        cindx2 = MOI.add_constraint(opt, saf2, LT(-value))
    end
    nothing
end
#=
function relax_quadratic_gen_saf(full_buffer::Vector{Float64}, full_active::Vector{Bool}, func::SQF, n::NodeBB, x0::Vector{Float64}, pos::Bool, nx::Int)

    fill!(full_buffer, 0.0)
    fill!(full_active, false)
    multiplier = pos ? 1.0 : - 1.0

    lower_bounds = n.lower_variable_bounds
    upper_bounds = n.upper_variable_bounds

    quadratic_constant = func.constant

    for term in func.quadratic_terms

        a = multiplier*term.coefficient
        idx1 = term.variable_index_1.value
        idx2 = term.variable_index_2.value

        x0_1 = @inbounds x0[idx1]
        xL_1 = @inbounds lower_bounds[idx1]
        xU_1 = @inbounds upper_bounds[idx1]

        if idx1 === idx2

            @inbounds full_active[idx1] = true
            @inbounds full_buffer[idx1] = (a > 0.0) ? 2.0*a*x0_1 : a*(xL_1 + xU_1)
            quadratic_constant -= (a > 0.0) ? x0_1*x0_1 : a*xL_1*xU_1

        else
            x0_2 = @inbounds x0[idx2]
            xL_2 = @inbounds lower_bounds[idx2]
            xU_2 = @inbounds upper_bounds[idx2]

            @inbounds full_active[idx1] = true
            @inbounds full_active[idx2] = true

            if a > 0.0
                check_ref = (xU_1 - xL_1)*x0_2 + (xU_2 - xL_2)*x0_1
                if check_ref <= xU_1*xU_2 - xL_1*xL_2
                    @inbounds full_buffer[idx1] += a*xL_2
                    @inbounds full_buffer[idx2] += a*xL_1
                    quadratic_constant -= a*xL_1*xL_2
                else
                    @inbounds full_buffer[idx1] += a*xU_2
                    @inbounds full_buffer[idx2] += a*xU_1
                    quadratic_constant -= a*xU_1*xU_2
                end
            else
                check_ref = (xU_1 - xL_1)*x0_2 - (xU_2 - xL_2)*x0_1
                if check_ref <= xU_1*xL_2 - xL_1*xU_2
                    @inbounds full_buffer[idx1] += a*xL_2
                    @inbounds full_buffer[idx2] += a*xU_1
                    quadratic_constant -= a*xU_1*xL_2
                else
                    @inbounds full_buffer[idx1] += a*xU_2
                    @inbounds full_buffer[idx2] += a*xL_1
                    quadratic_constant -= a*xL_1*xU_2
                end
            end
        end
    end

    for term in func.affine_terms
        a = multiplier*term.coefficient
        idx = term.variable_index.value
        @inbounds full_buffer[idx] += a
        @inbounds full_active[idx] = true
    end

    #saf = SAF(SAT.(terms_coeff, vi), quadratic_constant)

    return saf
end

function store_ge_quadratic!(x::Optimizer, ci::CI{SAF,LT}, saf::SAF, lower::Float64, i::Int64, q::Int64)
    opt = x.relaxed_optimizer
    if (q == 1) & x.relaxed_inplace_mod
        for (i, term) in enumerate(saf.terms)
            MOI.modify(opt, ci, SCoefC(term.variable_index, term.coefficient))
        end
        MOI.set(opt, MOI.ConstraintSet(), ci, LT(-lower - saf.constant))
    else
        saf.constant = 0.0
        x._quadratic_ci_geq[q][i] = MOI.add_constraint(opt, saf, LT(-lower - saf.constant))
    end
    return
end

function store_le_quadratic!(x::Optimizer, ci::CI{SAF,LT}, saf::SAF, upper::Float64, i::Int64, q::Int64)
    opt = x.relaxed_optimizer
    if (q == 1) & x.relaxed_inplace_mod
        opt = x.relaxed_optimizer
        for (i, term) in enumerate(saf.terms)
            MOI.modify(opt, ci, SCoefC(term.variable_index, term.coefficient))
        end
        MOI.set(opt, MOI.ConstraintSet(), ci, LT(upper - saf.constant))

    else
        saf.constant = 0.0
        cindxo = MOI.add_constraint(opt, saf, LT(upper - saf.constant))
        x._quadratic_ci_leq[q][i] = cindxo
    end
    return
end

function store_eq_quadratic!(x::Optimizer, ci1::CI{SAF,LT}, ci2::CI{SAF,LT}, saf1::SAF, saf2::SAF, value::Float64, i::Int64, q::Int64)
    opt = x.relaxed_optimizer
    if (q == 1) & x.relaxed_inplace_mod
        store_ge_quadratic!(x, ci1, saf1, value, i, q)
        store_le_quadratic!(x, ci2, saf2, value, i, q)
    else
        saf1.constant = 0.0
        saf2.constant = 0.0
        cindx1 = MOI.add_constraint(opt, saf1, LT(value - saf1.constant))
        cindx2 = MOI.add_constraint(opt, saf2, LT(-value + saf2.constant))
        x._quadratic_ci_eq[q][i] = (cindx1, cindx2)
    end
    return
end

"""
$(FUNCTIONNAME)

Relaxes all quadratic constraints in `x` optimizer.
"""
function relax_quadratic!(x::Optimizer, x0::Vector{Float64}, q::Int64)

    n = x._current_node

    # Relax Convex Constraint Terms TODO: place all quadratic info into one vector of tuples?
    for i = 1:length(x._quadratic_leq_constraints)
        func, set = x._quadratic_leq_constraints[i]
        vi = x._quadratic_leq_sparsity[i]
        nz = x._quadratic_leq_gradnz[i]
        ci1 = x._quadratic_ci_leq[q][i]

        saf = relax_quadratic_gen_saf(func, vi, n, nz, x0)
        store_le_quadratic!(x, ci1, saf, set.upper, i, q)
    end

    for i = 1:length(x._quadratic_geq_constraints)
        func, set = x._quadratic_geq_constraints[i]
        vi = x._quadratic_geq_sparsity[i]
        nz = x._quadratic_geq_gradnz[i]
        ci1 = x._quadratic_ci_geq[q][i]

        vec_sat = SAT[SAT(-t.coefficient, t.variable_index) for t in func.affine_terms]
        vec_sqt = SQT[SQT(-t.coefficient, t.variable_index_1, t.variable_index_2) for t in func.quadratic_terms]
        func_minus = SQF(vec_sat, vec_sqt, -func.constant)
        saf = relax_quadratic_gen_saf(func_minus, vi, n, nz, x0)
        store_ge_quadratic!(x, ci1, saf, set.lower, i, q)
    end

    for i = 1:length(x._quadratic_eq_constraints)
        func, set = x._quadratic_eq_constraints[i]
        vi = x._quadratic_eq_sparsity[i]
        nz = x._quadratic_eq_gradnz[i]
        ci1, ci2 = x._quadratic_ci_eq[q][i]

        saf1 = relax_quadratic_gen_saf(func, vi, n, nz, x0)

        vec_sat = SAT[SAT(-t.coefficient, t.variable_index) for t in func.affine_terms]
        vec_sqt = SQT[SQT(-t.coefficient, t.variable_index_1, t.variable_index_2) for t in func.quadratic_terms]
        func_minus = SQF(vec_sat, vec_sqt, -func.constant)
        saf2 = relax_quadratic_gen_saf(func_minus, vi, n, nz, x0)
        store_eq_quadratic!(x, ci1, ci2, saf1, saf2, set.value, i, q)
    end

    return
end
=#

"""
$(TYPEDSIGNATURES)

A rountine that only relaxes the objective.
"""
function relax_objective!(t::ExtensionType, x::Optimizer, x0::Vector{Float64})

    nx = x._variable_number
    opt = x.relaxed_optimizer
    @inbounds vi = x._lower_variable_index[1:nx]

    # Add objective
    if x._input_problem._objective_type === SINGLE_VARIABLE
        MOI.set(opt, MOI.ObjectiveFunction{SV}(), x._input_problem._objective_sv)
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    elseif x._input_problem._objective_type === SCALAR_AFFINE
        MOI.set(opt, MOI.ObjectiveFunction{SAF}(), x._input_problem._objective_saf)
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    elseif x._input_problem._objective_type === SCALAR_QUADRATIC

        # TODO: ADD QUADRATIC RELAXATION HERE
        MOI.set(opt, MOI.ObjectiveFunction{SAF}(), saf)
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    elseif x._objective_type === NONLINEAR

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
    return
end

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
    return
end

"""
$(TYPEDSIGNATURES)

A rountine that updates the current node for the `Evaluator` and relaxes all
nonlinear constraints and quadratic constraints.
"""
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
$(FUNCTIONNAME)

Adds linear objective cut constraint to the `x.relaxed_optimizer`.
"""
function objective_cut_linear!(m::Optimizer, q::Int64)
    if m._global_upper_bound < Inf
        if (m._objective_cut_set == -1) || (q > 1) ||  ~m.relaxed_inplace_mod
            z = (m._objective_cut_set == -1) ? 1 : q
            set = LT(m._global_upper_bound)
            if m._objective_type === SINGLE_VARIABLE
                ci_sv = m._objective_cut_ci_sv
                MOI.set(m.relaxed_optimizer, MOI.ConstraintSet(), ci_sv, set)
            elseif m._objective_type === SCALAR_AFFINE
                ci_saf = MOI.add_constraint(m.relaxed_optimizer, m._objective_saf, set)
                m._objective_cut_ci_saf[z] = ci_saf
            elseif m._objective_type === SCALAR_QUADRATIC || m._objective_type === NONLINEAR
                saf = MOI.get(m.relaxed_optimizer, MOI.ObjectiveFunction{SAF}())
                set = LT(m._global_upper_bound - saf.constant)
                saf.constant = 0.0
                ci_saf = MOI.add_constraint(m.relaxed_optimizer, saf, set)
                m._objective_cut_ci_saf[z] = ci_saf
            end
            m._objective_cut_set = m._iteration_count
        elseif q == 1
            if m._objective_type !== SINGLE_VARIABLE
                ci_saf = m._objective_cut_ci_saf[1]
                saf = MOI.get(m.relaxed_optimizer, MOI.ObjectiveFunction{SAF}())
                for term in saf.terms
                    term_coeff = term.coefficient
                    term_indx = term.variable_index
                    MOI.modify(m.relaxed_optimizer, ci_saf, SCoefC(term_indx, term_coeff))
                end
                #MOI.modify(x.relaxed_optimizer, ci_saf, SConsC(0.0))
                set = LT(m._global_upper_bound - saf.constant)
                MOI.set(m.relaxed_optimizer, MOI.ConstraintSet(), ci_saf, set)
            else
                ci_sv = m._objective_cut_ci_sv
                set = LT(m._global_upper_bound)
                MOI.set(m.relaxed_optimizer, MOI.ConstraintSet(), ci_sv, set)
            end
        end
    end
    return
end

relax_problem!(m::Optimizer, x::Vector{Float64}, q::Int64) = relax_problem!(m.ext_type, m, x, q)
relax_objective!(m::Optimizer, x::Vector{Float64}) = relax_objective!(m.ext_type, m, x)

#=
work on new constraint relaxation
"""
Takes the optimizer and constraint index and computes a relaxation inplace if
possible (no prior relaxation) and no set.
"""

function copy_add_from_buffer!(x::Optimizer, c::CI{})
end
function relax_to_buffer!(x, c)
end
function is_safe_relax(x, c)
    buffer = x.buffer[x.buffer_indx[c]]
    coeffs = buffer.coeffs
    flag = abs(buffer.b) < x.cut_safe_b
    ~flag && (return flag)
    for i=1:length(coeffs)
        ai = coeffs[i]
        if ~iszero(ai)
            (x.cut_safe_l > abs(ai)) && (flag = false; break)
            (x.cut_safe_u < abs(ai)) && (flag = false; break)
            for j=1:length(coeffs)
                aj = coeffs[j]
                if ~iszero(coeffs[j])
                    d = abs(ai/aj)
                    (x.cut_safe_l > d) && (flag = false; break)
                    (x.cut_safe_u < d) && (flag = false; break)
                end
            end
        end
    end
    return flag
end
function relax_expr!(x::Optimizer, c::CI)
    relax_to_buffer!(x, c)
    is_safe_relax(x, c) && copy_add_from_buffer!(x, c)
    x._cut_number[x] += 1
    nothing
end
=#
