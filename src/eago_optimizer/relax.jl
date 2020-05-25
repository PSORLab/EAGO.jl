"""
$(FUNCTIONNAME)

Applies the safe cut checks detailed in Khajavirad, 2018 [Khajavirad, Aida, and Nikolaos V. Sahinidis. "A hybrid LP/NLP
paradigm for global optimization relaxations." Mathematical Programming Computation 10.3 (2018): 383-421] to ensure
that only numerically safe affine relaxations are added. Checks that i) |b| <= safe_b, ii) safe_l <= abs(ai) <= safe_u,
and iii) violates safe_l <= abs(ai/aj) <= safe_u.
"""
function is_safe_cut!(m::Optimizer, f::SAF)

    safe_l = m._parameters.cut_safe_l
    safe_u = m._parameters.cut_safe_u
    safe_b = m._parameters.cut_safe_b

    (abs(f.constant) > safe_b) && return false                          # violates |b| <= safe_b

    term_count = length(f.terms)
    for i = 1:term_count

        ai = @inbounds f.terms[i]
        if ai !== 0.0

            ai_abs = abs(ai)
            !(safe_l <= abs(ai) <= safe_u) && return false              # violates safe_l <= abs(ai) <= safe_u

            for j = i:term_count
                aj = @inbounds f.terms[j]
                if aj !== 0.0
                    !(safe_l <= abs(ai/aj) <= safe_u) && return false   # violates safe_l <= abs(ai/aj) <= safe_u
                end
            end
        end
    end

    return true
end

"""
Relaxs the constraint adding an affine function.
"""
function relax! end

"""
$(FUNCTIONNAME)

Default routine for relaxing quadratic constraint `func` < `0.0` on node `n`. Takes affine bounds of convex part at
point `x0` and secant line bounds on concave parts.
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

const LT_ZERO = LT(0.0)
macro define_relax_ineq(function_type, ci_array_name)
    quote
        function relax!(m::Optimizer, f::$function_type, indx::Int, check_safe::Bool)

            affine_relax_quadratic!(f.func, f.buffer, f.saf, m._current_node, m._current_xref, true)
            if check_safe && is_safe_cut!(m, f.saf)
                ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, LT_ZERO)
                push!(m.$(ci_array_name), ci)
            end
            m.relaxed_to_problem_map[ci] = indx

            nothing
        end
    end
end

@define_relax_ineq BufferedQuadraticIneq _buffered_quadratic_ineq_ci

macro define_relax_eq(function_type, ci_array_name)
    quote
        function relax!(m::Optimizer, f::$function_type, indx::Int, check_safe::Bool)

            affine_relax_quadratic!(f.func, f.buffer, f.saf1, m._current_node, m._current_xref, true)
            if check_safe && is_safe_cut!(m, f.saf1)
                ci = MOI.add_constraint(m.relaxed_optimizer, f.saf1, LT_ZERO)
                push!(m.$(ci_array_name), ci)
            end
            m.relaxed_to_problem_map[ci] = indx

            affine_relax_quadratic!(f.function, f.buffer, f.saf2, m._current_node, m._current_xref, false)
            if check_safe && is_safe_cut!(m, f.saf2)
                ci = MOI.add_constraint(m.relaxed_optimizer, f.saf2, LT_ZERO)
                push!(m.$(ci_array_name), ci)
            end
            m.relaxed_to_problem_map[ci] = indx

            nothing
        end
    end
end

@define_relax_ineq BufferedQuadraticEq _buffered_quadratic_eq_ci

function relax_midpoint!(m::Optimizer)

    # set midpoint
    set_reference_midpoint!(m)

    # adds a cut for each scalar quadratic function at the midpoint
    for i = 1:m._sqf_leq_count
        relax_midpoint!(m, m._sqf_leq, i)
    end

        relax_midpoint!.(m, m_sqf_eq)

    # adds a cut for each nonlinear function at the midpoint
    relax!.(m._nl_leq)
    relax!.(m._nl_geq)
end

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

function delete_nl_relaxations!(m::Optimizer)

    # delete affine relaxations added from quadatic inequality
    for ci in m._buffered_quadratic_ineq_ci
        MOI.delete(m, ci)
    end

    # delete affine relaxations added from quadratic equality
    for ci in m._buffered_quadratic_eq_ci
        MOI.delete(m, ci)
    end
end
