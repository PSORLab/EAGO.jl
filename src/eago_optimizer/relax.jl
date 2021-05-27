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
    return
end

"""
$(TYPEDSIGNATURES)
"""
function relax!(m::Optimizer, f::BufferedNonlinearFunction{MC{N,T}}, k::Int, check_safe::Bool) where {N,T<:RelaxTag}
    d = m._working_problem._relaxed_evaluator
    check_set_affine_nl!(m, f, affine_relax_nonlinear!(f, d, true, true, true), check_safe)
    check_set_affine_nl!(m, f, affine_relax_nonlinear!(f, d, false, false, true), check_safe)
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

Deletes all constraints corresponding to relaxations of nonlinear terms
added to the relaxed optimizer and clear buffers of constraint indicies.
"""
function delete_nl_constraints!(m::Optimizer)
    foreach(c -> MOI.delete(m.relaxed_optimizer, c), m._buffered_quadratic_ineq_ci)
    foreach(c -> MOI.delete(m.relaxed_optimizer, c), m._buffered_quadratic_eq_ci)
    foreach(c -> MOI.delete(m.relaxed_optimizer, c), m._buffered_nonlinear_ci)
    empty!(m._buffered_quadratic_ineq_ci)
    empty!(m._buffered_quadratic_eq_ci)
    empty!(m._buffered_nonlinear_ci)
    return
end

"""
$(FUNCTIONNAME)

Deletes all scalar-affine objective cuts added to the relaxed optimizer.
"""
function delete_objective_cuts!(m::Optimizer)
    foreach(c -> MOI.delete(m.relaxed_optimizer, c), m._objective_cut_ci_saf)
    empty!(m._objective_cut_ci_saf)
    return
end
