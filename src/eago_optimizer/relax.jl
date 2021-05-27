# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/relax.jl
# Defines routines used construct the relaxed subproblem.
#############################################################################

"""
$(FUNCTIONNAME)

Applies the safe cut checks detailed in Khajavirad, 2018 [Khajavirad, Aida,
and Nikolaos V. Sahinidis. "A hybrid LP/NLP paradigm for global optimization
relaxations." Mathematical Programming Computation 10.3 (2018): 383-421] to
ensure that only numerically safe affine relaxations are added. Checks that
i) ``|b| <= safe b`, ii) `safe_l <= abs(ai) <= safe u`, and iii) violates
`safe_l <= abs(ai/aj) <= safe_u`.
"""
function is_safe_cut!(m::Optimizer, f::SAF)

    safe_l = m._parameters.cut_safe_l
    safe_u = m._parameters.cut_safe_u
    safe_b = m._parameters.cut_safe_b

    # violates |b| <= safe_b
    (abs(f.constant) > safe_b) && return false

    term_count = length(f.terms)
    for i = 1:term_count

        ai = (@inbounds f.terms[i]).coefficient
        if ai !== 0.0

            # violates safe_l <= abs(ai) <= safe_u
            ai_abs = abs(ai)
            !(safe_l <= abs(ai) <= safe_u) && return false

            # violates safe_l <= abs(ai/aj) <= safe_u
            for j = i:term_count
                aj = (@inbounds f.terms[j]).coefficient
                if aj !== 0.0
                    !(safe_l <= abs(ai/aj) <= safe_u) && return false
                end
            end
        end
    end

    return true
end

"""
$(FUNCTIONNAME)

Relaxs the constraint by adding an affine constraint to the model.
"""
function relax! end

"""
$(FUNCTIONNAME)

Default routine for relaxing quadratic constraint `func < 0.0` on node `n`.
Takes affine bounds of convex part at point `x0` and secant line bounds on
concave parts.
"""
function affine_relax_quadratic!(func::SQF, buffer::Dict{Int,Float64}, saf::SAF,
                                 n::NodeBB, sol_to_branch_map::Vector{Int},
                                 x::Vector{Float64})

    lower_bounds = n.lower_variable_bounds
    upper_bounds = n.upper_variable_bounds
    quadratic_constant = func.constant

    # Affine terms only contribute coefficients, so the respective
    # values do not contribute to the cut. Since all quadratic terms
    # are considered to be branch variables we exclude any potential
    # need to retrieve variable bounds from locations other than
    # the node.
    for term in func.quadratic_terms

        a = term.coefficient
        idx1 = term.variable_index_1.value
        idx2 = term.variable_index_2.value
        sol_idx1 = sol_to_branch_map[idx1]
        sol_idx2 = sol_to_branch_map[idx2]
        x0_1 = x[sol_idx1]
        xL_1 = lower_bounds[sol_idx1]
        xU_1 = upper_bounds[sol_idx1]

        if idx1 === idx2

            if a > 0.0
                buffer[idx1] += a*x0_1
                quadratic_constant -= 0.5*a*x0_1*x0_1

            else
                if !isinf(xL_1) && !isinf(xU_1)
                    buffer[idx1] += 0.5*a*(xL_1 + xU_1)
                    quadratic_constant -= 0.5*a*xL_1*xU_1
                else
                    return false
                end
            end

        else
            x0_2 = x[sol_idx2]
            xL_2 = lower_bounds[sol_idx2]
            xU_2 = upper_bounds[sol_idx2]

            if a > 0.0
                if (!isinf(xL_1) && !isinf(xL_2)) &&
                   ((xU_1 - xL_1)*x0_2 + (xU_2 - xL_2)*x0_1 <= xU_1*xU_2 - xL_1*xL_2)
                    buffer[idx1] += a*xL_2
                    buffer[idx2] += a*xL_1
                    quadratic_constant -= a*xL_1*xL_2

                elseif !isinf(xU_1) && !isinf(xU_2)
                    buffer[idx1] += a*xU_2
                    buffer[idx2] += a*xU_1
                    quadratic_constant -= a*xU_1*xU_2

                else
                    return false

                end
            else
                if (!isinf(xU_1) && !isinf(xL_2)) &&
                   ((xU_1 - xL_1)*x0_2 - (xU_2 - xL_2)*x0_1 <= xU_1*xL_2 - xL_1*xU_2)

                    buffer[idx1] += a*xL_2
                    buffer[idx2] += a*xU_1
                    quadratic_constant -= a*xU_1*xL_2

                elseif !isinf(xL_1) && !isinf(xU_2)
                    buffer[idx1] += a*xU_2
                    buffer[idx2] += a*xL_1
                    quadratic_constant -= a*xL_1*xU_2

                else
                    return false
                end
            end
        end
    end

    for term in func.affine_terms
        a0 = term.coefficient
        idx = term.variable_index.value
        buffer[idx] += a0
    end

    count = 1
    for (key, value) in buffer
        saf.terms[count] = SAT(value, VI(key))
        buffer[key] = 0.0
        count += 1
    end
    saf.constant = quadratic_constant

    return true
end

"""
$(TYPEDSIGNATURES)
"""
function relax!(m::Optimizer, f::BufferedQuadraticIneq, indx::Int, check_safe::Bool)

    constraint_tol = m._parameters.absolute_constraint_feas_tolerance
    finite_cut_generated = affine_relax_quadratic!(f.func, f.buffer, f.saf, m._current_node, m._sol_to_branch_map, m._current_xref)
    if finite_cut_generated
        if !check_safe || is_safe_cut!(m, f.saf)
            lt = LT(-f.saf.constant + constraint_tol)
            f.saf.constant = 0.0
            ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, lt)
            push!(m._buffered_quadratic_ineq_ci, ci)
        end
    end
    #m.relaxed_to_problem_map[ci] = indx

    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function relax!(m::Optimizer, f::BufferedQuadraticEq, indx::Int, check_safe::Bool)

    constraint_tol = m._parameters.absolute_constraint_feas_tolerance
    finite_cut_generated = affine_relax_quadratic!(f.func, f.buffer, f.saf, m._current_node, m._sol_to_branch_map, m._current_xref)
    if finite_cut_generated
        if !check_safe || is_safe_cut!(m, f.saf)
            lt = LT(-f.saf.constant + constraint_tol)
            f.saf.constant = 0.0
            ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, lt)
            push!(m._buffered_quadratic_eq_ci, ci)
        end
    end
    #m.relaxed_to_problem_map[ci] = indx

    finite_cut_generated = affine_relax_quadratic!(f.minus_func, f.buffer, f.saf, m._current_node, m._sol_to_branch_map, m._current_xref)
    if finite_cut_generated
        if !check_safe || is_safe_cut!(m, f.saf)
            lt = LT(-f.saf.constant + constraint_tol)
            f.saf.constant = 0.0
            ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, lt)
            push!(m._buffered_quadratic_eq_ci, ci)
        end
    end
    #m.relaxed_to_problem_map[ci] = indx

    return nothing
end

"""
$(FUNCTIONNAME)
"""
function affine_relax_nonlinear!(f::BufferedNonlinearFunction{MC{N,T}}, evaluator::Evaluator,
                                 use_cvx::Bool, new_pass::Bool, is_constraint::Bool) where {N,T<:RelaxTag}

    new_pass && forward_pass!(evaluator, f)
    x = evaluator.variable_values.x
    finite_cut = true

    grad_sparsity = _sparsity(f)
    if _is_num(f)
        f.saf.constant = _num(f)
        for i = 1:N
            vval = @inbounds grad_sparsity[i]
            f.saf.terms[i] = SAT(0.0, VI(vval))
        end

    else
        setvalue = _set(f)
        finite_cut &= !(isempty(setvalue) || isnan(setvalue))

        if finite_cut
            value = _set(f)
            f.saf.constant = use_cvx ? value.cv : -value.cc
            for i = 1:N
                vval = @inbounds grad_sparsity[i]
                if use_cvx
                    coef = @inbounds value.cv_grad[i]
                else
                    coef = @inbounds -value.cc_grad[i]
                end
                f.saf.terms[i] = SAT(coef, VI(vval))
                xv = @inbounds x[vval]
                f.saf.constant = sub_round(f.saf.constant , mul_round(coef, xv, RoundUp), RoundDown)
            end
            if is_constraint
                bnd_used =  use_cvx ? -_upper_bound(f) : _lower_bound(f)
                f.saf.constant = add_round(f.saf.constant, bnd_used, RoundDown)
            end
        end
    end

    return finite_cut
end

"""
$(TYPEDSIGNATURES)
"""
function check_set_affine_nl!(m::Optimizer, f::BufferedNonlinearFunction{MC{N,T}}, finite_cut_generated::Bool, check_safe::Bool) where {N,T<:RelaxTag}

    constraint_tol = m._parameters.absolute_constraint_feas_tolerance
    if finite_cut_generated
        if !check_safe || is_safe_cut!(m, f.saf)
            lt = LT(-f.saf.constant + constraint_tol)
            f.saf.constant = 0.0
            ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, lt)
            push!(m._buffered_nonlinear_ci, ci)
        end
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function relax!(m::Optimizer, f::BufferedNonlinearFunction{MC{N,T}}, indx::Int, check_safe::Bool) where {N,T<:RelaxTag}
    evaluator = m._working_problem._relaxed_evaluator

    finite_cut_generated = affine_relax_nonlinear!(f, evaluator, true, true, true)
    check_set_affine_nl!(m, f, finite_cut_generated, check_safe)

    finite_cut_generated = affine_relax_nonlinear!(f, evaluator, false, false, true)
    check_set_affine_nl!(m, f, finite_cut_generated, check_safe)

    return nothing
end

bound_objective(m::Optimizer, f::BufferedNonlinearFunction, n::NodeBB) = lower_interval_bound(m, f, n)
bound_objective(m::Optimizer, f::AffineFunctionIneq, n::NodeBB) = lower_interval_bound(m, f, n)
bound_objective(m::Optimizer, f::BufferedQuadraticIneq, n::NodeBB) = lower_interval_bound(m, f, n)
bound_objective(m::Optimizer, f::SV, n::NodeBB) = _lower_bound(FullVar(), m, f.variable.value)
function bound_objective(t::ExtensionType, m::Optimizer)
    bound_objective(m, m._working_problem._objective, m._current_node)
end
bound_objective(m::Optimizer) = bound_objective(m.ext_type, m)

function relax_objective!(m, f::T, check_safe) where T<:Union{SV,SAF}
    MOI.set(m.relaxed_optimizer, MOI.ObjectiveFunction{T}(), f)
    return
end
function relax_objective!(m, f::BufferedQuadraticIneq, check_safe)
    wp = m._working_problem
    finite_cut_generated = affine_relax_quadratic!(f.func, f.buffer, f.saf,
                                                   m._current_node, m._sol_to_branch_map,
                                                   m._current_xref)
    if finite_cut_generated && (!check_safe || is_safe_cut!(m, f.saf))
        copyto!(wp._objective_saf.terms, f.saf.terms)
        wp._objective_saf.constant = f.saf.constant
        MOI.set(m.relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), wp._objective_saf)
    end
    return
end
function relax_objective!(m, f::BufferedNonlinearFunction, check_safe)
    wp = m._working_problem
    new_flag = m._new_eval_objective
    wp._relaxed_evaluator.is_first_eval = new_flag
    finite_cut_generated = affine_relax_nonlinear!(f, wp._relaxed_evaluator, true, new_flag, false)
    wp._relaxed_evaluator.is_first_eval = false

    if finite_cut_generated && (!check_safe || is_safe_cut!(m, f.saf))
        copyto!(wp._objective_saf.terms, f.saf.terms)
        wp._objective_saf.constant = f.saf.constant
        MOI.set(m.relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), wp._objective_saf)
    end
    return
end

"""
$(TYPEDSIGNATURES)

Triggers an evaluation of the objective function and then updates
the affine relaxation of the objective function.
"""
function relax_objective!(t::ExtensionType, m::Optimizer, q::Int64)
    check_safe = (q === 1) ? false : m._parameters.cut_safe_on
    relax_objective!(m, m._working_problem._objective, check_safe)
    m._new_eval_objective = false
    return
end
relax_objective!(m::Optimizer, q::Int64) = relax_objective!(m.ext_type, m, q)

function objective_cut!(m::Optimizer, f::SV, check_safe::Bool, ϵ)
    if !isinf(m._global_upper_bound) && (m._objective_cut_ci_sv.value === -1)
        m._objective_cut_ci_sv = CI{SV,LT}(f.variable.value)
    end
    MOI.set(m.relaxed_optimizer, MOI.ConstraintSet(), m._objective_cut_ci_sv, LT(UBD))
    return
end
function objective_cut!(m::Optimizer, f::AffineFunctionIneq, check_safe::Bool, ϵ)
    formulated_constant = wp._objective_saf.constant
    wp._objective_saf.constant = 0.0
    if check_safe && is_safe_cut!(m, wp._objective_saf)
        ci_saf = MOI.add_constraint(m.relaxed_optimizer, wp._objective_saf, LT(UBD - wp._objective_saf.constant + ϵ))
        push!(m._objective_cut_ci_saf, ci_saf)
    end
    wp._objective_saf.constant = formulated_constant
    return
end
function objective_cut!(m::Optimizer, f::BufferedQuadraticIneq, check_safe::Bool)
    finite_cut_generated = affine_relax_quadratic!(f.func, f.buffer,
                                                   f.saf, m._current_node, m._sol_to_branch_map,
                                                   m._current_xref)
    if finite_cut_generated && (!check_safe || is_safe_cut!(m, f.saf))
        copyto!(wp._objective_saf.terms, f.saf.terms)
        wp._objective_saf.constant = 0.0
        ci_saf = MOI.add_constraint(m.relaxed_optimizer, wp._objective_saf, LT(UBD - buffered_sqf.saf.constant + ϵ))
        push!(m._objective_cut_ci_saf, ci_saf)
    end
    return
end
function objective_cut!(m::Optimizer, f::BufferedNonlinearFunction, check_safe::Bool, ϵ)

    UBD = m._global_upper_bound
    relaxed_optimizer = m.relaxed_optimizer
    wp = m._working_problem
    relaxed_evaluator = wp._relaxed_evaluator

    # if the objective cut is the first evaluation of the objective expression
    # then perform a a forward pass
    new_flag = m._new_eval_objective
    relaxed_evaluator.is_first_eval = new_flag
    finite_cut_generated = affine_relax_nonlinear!(f, relaxed_evaluator, true, new_flag, false)

    constraint_tol = m._parameters.absolute_constraint_feas_tolerance
    if finite_cut_generated
        copyto!(wp._objective_saf.terms, f.saf.terms)
        wp._objective_saf.constant = 0.0
        if !check_safe || is_safe_cut!(m,  f.saf)
            # TODO: When we introduce numerically safe McCormick operators we'll need to replace
            # the UBD - buffered_nl.saf.constant with a correctly rounded version. For now,
            # a small factor is added to the UBD calculation initially which should be sufficient.
            ci_saf = MOI.add_constraint(m.relaxed_optimizer, wp._objective_saf, LT(UBD - f.saf.constant + ϵ))
            push!(m._objective_cut_ci_saf, ci_saf)
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Adds linear objective cut constraint to the `x.relaxed_optimizer`.
"""
function objective_cut!(m::Optimizer, check_safe::Bool)
    ϵ = m._parameters.absolute_constraint_feas_tolerance
    if m._parameters.objective_cut_on && m._global_upper_bound < Inf
        objective_cut!(m, m._working_problem._objective, check_safe, ϵ)
        m._new_eval_objective = false
    end
    return
end

"""
$(TYPEDSIGNATURES)

A routine that adds relaxations for all nonlinear constraints and quadratic constraints
corresponding to the current node to the relaxed problem. This adds an objective cut
(if specified by `objective_cut_on`) and then sets the `_new_eval_constraint` flag
to false indicating that an initial evaluation of the constraints has occurred. If
the `objective_cut_on` flag is `true` then the `_new_eval_objective` flag is also
set to `false` indicating that the objective expression was evaluated.
"""
function relax_all_constraints!(t::ExtensionType, m::Optimizer, q::Int64)

    check_safe = (q === 1) ? false : m._parameters.cut_safe_on
    m._working_problem._relaxed_evaluator.is_first_eval = m._new_eval_constraint

    sqf_leq_list = m._working_problem._sqf_leq
    for i = 1:m._working_problem._sqf_leq_count
        sqf_leq = @inbounds sqf_leq_list[i]
        relax!(m, sqf_leq, i, check_safe)
    end

    sqf_eq_list = m._working_problem._sqf_eq
    for i = 1:m._working_problem._sqf_eq_count
        sqf_eq = @inbounds sqf_eq_list[i]
        relax!(m, sqf_eq, i, check_safe)
    end

    nl_list = m._working_problem._nonlinear_constr
    for i = 1:m._working_problem._nonlinear_count
        nl = @inbounds nl_list[i]
        relax!(m, nl, i, check_safe)
    end

    m._new_eval_constraint = false

    objective_cut!(m, check_safe)

    return
end
relax_constraints!(t::ExtensionType, m::Optimizer, q::Int64) = relax_all_constraints!(t, m, q)
relax_constraints!(m::Optimizer, q::Int64) = relax_constraints!(m.ext_type, m, q)

"""
$(FUNCTIONNAME)

Deletes all nonlinear constraints added to the relaxed optimizer.
"""
function delete_nl_constraints!(m::Optimizer)

    # delete affine relaxations added from quadratic inequality
    for ci in m._buffered_quadratic_ineq_ci
        MOI.delete(m.relaxed_optimizer, ci)
    end
    empty!(m._buffered_quadratic_ineq_ci)

    # delete affine relaxations added from quadratic equality
    for ci in m._buffered_quadratic_eq_ci
        MOI.delete(m.relaxed_optimizer, ci)
    end
    empty!(m._buffered_quadratic_eq_ci)

    # delete affine relaxations added from nonlinear inequality
    for ci in m._buffered_nonlinear_ci
        MOI.delete(m.relaxed_optimizer, ci)
    end
    empty!(m._buffered_nonlinear_ci)

    return
end

"""
$(FUNCTIONNAME)

Deletes all scalar-affine objective cuts added to the relaxed optimizer.
"""
function delete_objective_cuts!(m::Optimizer)
    for ci in m._objective_cut_ci_saf
        MOI.delete(m.relaxed_optimizer, ci)
    end
    empty!(m._objective_cut_ci_saf)
    return
end

"""
$(FUNCTIONNAME)

"""
function set_first_relax_point!(m::Optimizer)
    m._working_problem._relaxed_evaluator.is_first_eval = true
    m._new_eval_constraint = true
    m._new_eval_objective = true
    n = m._current_node

    m._current_xref .= mid(n)
    unsafe_check_fill!(isnan, m._current_xref, 0.0, length(m._current_xref))
    return
end
