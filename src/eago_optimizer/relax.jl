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
# TODO
#############################################################################

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

        ai = (@inbounds f.terms[i]).coefficient
        if ai !== 0.0

            ai_abs = abs(ai)
            !(safe_l <= abs(ai) <= safe_u) && return false              # violates safe_l <= abs(ai) <= safe_u

            for j = i:term_count
                aj = (@inbounds f.terms[j]).coefficient
                if aj !== 0.0
                    !(safe_l <= abs(ai/aj) <= safe_u) && return false   # violates safe_l <= abs(ai/aj) <= safe_u
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

Default routine for relaxing quadratic constraint `func` < `0.0` on node `n`. Takes affine bounds of convex part at
point `x0` and secant line bounds on concave parts.
"""
function affine_relax_quadratic!(func::SQF, buffer::Dict{Int,Float64}, saf::SAF,
                                 n::NodeBB, sol_to_branch_map::Vector{Int},
                                 x::Vector{Float64})

    lower_bounds = n.lower_variable_bounds
    upper_bounds = n.upper_variable_bounds
    quadratic_constant = func.constant

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
                buffer[idx1] += 2.0*a*x0_1
                quadratic_constant -= a*x0_1*x0_1

            else
                if !isinf(xL_1) && !isinf(xU_1)
                    buffer[idx1] += a*(xL_1 + xU_1)
                    quadratic_constant -= a*xL_1*xU_1
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

function relax!(m::Optimizer, f::BufferedQuadraticIneq, indx::Int, check_safe::Bool)

    finite_cut_generated = affine_relax_quadratic!(f.func, f.buffer, f.saf, m._current_node, m._sol_to_branch_map, m._current_xref)
    if finite_cut_generated
        if !check_safe || is_safe_cut!(m, f.saf)
            lt = LT(-f.saf.constant)
            f.saf.constant = 0.0
            ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, lt)
            push!(m._buffered_quadratic_ineq_ci, ci)
        end
    end
    #m.relaxed_to_problem_map[ci] = indx

    return nothing
end

function relax!(m::Optimizer, f::BufferedQuadraticEq, indx::Int, check_safe::Bool)

    finite_cut_generated = affine_relax_quadratic!(f.func, f.buffer, f.saf, m._current_node, m._sol_to_branch_map, m._current_xref)
    if finite_cut_generated
        if !check_safe || is_safe_cut!(m, f.saf)
            lt = LT(-f.saf.constant)
            f.saf.constant = 0.0
            ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, lt)
            push!(m._buffered_quadratic_eq_ci, ci)
        end
    end
    #m.relaxed_to_problem_map[ci] = indx

    finite_cut_generated = affine_relax_quadratic!(f.minus_func, f.buffer, f.saf, m._current_node, m._sol_to_branch_map, m._current_xref)
    if finite_cut_generated
        if !check_safe || is_safe_cut!(m, f.saf)
            lt = LT(-f.saf.constant)
            f.saf.constant = 0.0
            ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, lt)
            push!(m._buffered_quadratic_eq_ci, ci)
        end
    end
    #m.relaxed_to_problem_map[ci] = indx

    return nothing
end


function bound_objective(t::ExtensionType, m::Optimizer)

    n = m._current_node
    sb_map = m._sol_to_branch_map
    wp = m._working_problem
    obj_type = wp._objective_type

    if obj_type === NONLINEAR
        #=
        objective_lo = eval_objective_lo(d)
        constraints = d.constraints
        constr_num = d.constraint_number
        constraints_intv_lo = zeros(Float64, constr_num)
        constraints_intv_hi = zeros(Float64, constr_num)
        eval_constraint_lo!(d, constraints_intv_lo)
        eval_constraint_hi!(d, constraints_intv_hi)
        constraints_bnd_lo = d.constraints_lbd
        constraints_bnd_hi = d.constraints_ubd

        for i = 1:d.constraint_number
            @inbounds constraints_intv_lo = constraints_bnd_lo[i]
            @inbounds constraints_intv_hi = constraints_bnd_hi[i]
            if (constraints_intv_lo > constraints_intv_hi) || (constraints_intv_hi < constraints_intv_lo)
                feas = false
                break
            end
        end
        =#

    elseif obj_type === SINGLE_VARIABLE
        obj_indx = @inbounds sb_map[wp._objective_sv.variable.value]
        objective_lo = @inbounds n.lower_variable_bounds[obj_indx]

    elseif obj_type === SCALAR_AFFINE
        objective_lo = lower_interval_bound(wp._objective_saf_parsed, n)

    elseif obj_type === SCALAR_QUADRATIC
        objective_lo = lower_interval_bound(wp._objective_sqf, n)

    end

    return objective_lo
end
bound_objective(m::Optimizer) = bound_objective(m.ext_type, m)

"""
$(TYPEDSIGNATURES)

A rountine that only relaxes the objective.
"""
function relax_objective!(t::ExtensionType, m::Optimizer, q::Int64)

    relaxed_optimizer = m.relaxed_optimizer

    # Add objective
    wp = m._working_problem
    obj_type = wp._objective_type
    check_safe = (q === 1) ? false : m._parameters.cut_safe_on

    if obj_type === SINGLE_VARIABLE
        MOI.set(relaxed_optimizer, MOI.ObjectiveFunction{SV}(), wp._objective_sv)

    elseif obj_type === SCALAR_AFFINE
        MOI.set(relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), wp._objective_saf)

    elseif obj_type === SCALAR_QUADRATIC
        buffered_sqf = wp._objective_sqf
        finite_cut_generated = affine_relax_quadratic!(buffered_sqf.func, buffered_sqf.buffer, buffered_sqf.saf,
                                m._current_node, m._sol_to_branch_map, m._current_xref)
        if finite_cut_generated
            if !check_safe || is_safe_cut!(m, buffered_sqf.saf)
                copyto!(wp._objective_saf.terms, buffered_sqf.saf.terms)
                wp._objective_saf.constant = buffered_sqf.saf.constant
                MOI.set(relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), wp._objective_saf)
            end
        end

    #=
    elseif obj_type === NONLINEAR


    =#
    end
    return nothing
end
relax_objective!(m::Optimizer, q::Int64) = relax_objective!(m.ext_type, m, q)

"""
$(FUNCTIONNAME)

Adds linear objective cut constraint to the `x.relaxed_optimizer`.
"""
function objective_cut!(m::Optimizer, check_safe::Bool)

    UBD = m._global_upper_bound
    if m._parameters.objective_cut_on && m._global_upper_bound < Inf

        wp = m._working_problem
        obj_type = wp._objective_type

        if obj_type === SINGLE_VARIABLE
            if !isinf(UBD) && (m._objective_cut_ci_sv.value === -1)
                m._objective_cut_ci_sv = MOI.add_constraint(m.relaxed_optimizer, wp._objective_sv, LT_ZERO)
            else
                MOI.set(m.relaxed_optimizer, MOI.ConstraintSet(), m._objective_cut_ci_sv, LT(UBD))
            end
        elseif obj_type === SCALAR_AFFINE
            wp._objective_saf.constant -= UBD
            relax!(wp._objective_saf)
            if check_safe && is_safe_cut!(m, wp._objective_saf)
                ci_saf = MOI.add_constraint(m.relaxed_optimizer, wp._objective_saf, LT_ZERO)
                push!(m._objective_cut_ci_saf, ci_saf)
            end
            wp._objective_saf.constant += UBD

        elseif obj_type === SCALAR_QUADRATIC
            buffered_sqf = wp._objective_sqf
            finite_cut_generated =  affine_relax_quadratic!(buffered_sqf.func, buffered_sqf.buffer,
                                                            buffered_sqf.saf, m._current_node, m._sol_to_branch_map,
                                                            m._current_xref)

            if finite_cut_generated
                if !check_safe || is_safe_cut!(m, buffered_sqf.saf)
                    copyto!(wp._objective_saf.terms, buffered_sqf.saf.terms)
                    wp._objective_saf.constant = 0.0
                    ci_saf = MOI.add_constraint(m.relaxed_optimizer, wp._objective_saf, LT(UBD - buffered_sqf.saf.constant))
                    push!(m._objective_cut_ci_saf, ci_saf)
                end
            end
        #=
        elseif obj_type === NONLINEAR
            buffered_nl = wp._objective_nl
            finite_cut_generated =  affine_relax_nl!(TODO)
            if finite_cut_generated
                copyto!(wp._objective_saf.terms, buffered_nl.saf.terms)
                wp._objective_saf.constant = 0.0
                if !check_safe || is_safe_cut!(m,  buffered_nl.saf)
                    ci_saf = MOI.add_constraint(m.relaxed_optimizer, wp._objective_saf, LT(UBD - buffered_nl.saf.constant))
                    push!(m._objective_cut_ci_saf, ci_saf)
                end
            end
        =#
        end
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)

A rountine that updates the current node for the `Evaluator` and relaxes all
nonlinear constraints and quadratic constraints.
"""
function relax_all_constraints!(t::ExtensionType, m::Optimizer, q::Int64)

    check_safe = (q === 1) ? false : m._parameters.cut_safe_on

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

    #=
    nl_leq_list = m._working_problem._nl_leq
    for i = 1:m._working_problem._nl_leq_count
        nl_leq = @inbounds nl_leq_list[i]
        relax!(m, nl_leq, i, check_safe)
    end

    nl_eq_list = m._working_problem._nl_eq
    for i = 1:m._working_problem._nl_eq_count
        nl_eq = @inbounds nl_eq_list[i]
        relax!(m, nl_eq, i, check_safe)
    end
    =#

    objective_cut!(m, check_safe)

    return nothing
end
relax_constraints!(t::ExtensionType, m::Optimizer, q::Int64) = relax_all_constraints!(t, m, q)
relax_constraints!(m::Optimizer, q::Int64) = relax_constraints!(m.ext_type, m, q)

"""

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

    #=
    # delete affine relaxations added from nonlinear inequality
    for ci in m._buffered_nonlinear_ineq_ci
        MOI.delete(m.relaxed_optimizer, ci)
    end
    empty!(m._buffered_nonlinear_ineq_ci)

    # delete affine relaxations added from nonlinear equality
    for ci in m._buffered_nonlinear_eq_ci
        MOI.delete(m.relaxed_optimizer, ci)
    end
    empty!(m._buffered_nonlinear_eq_ci)
    =#

    return nothing
end

"""
Deletes all scalar-affine objective cuts added to the relaxed optimizer.
"""
function delete_objective_cuts!(m::Optimizer)

    for ci in m._objective_cut_ci_saf
        MOI.delete(m.relaxed_optimizer, ci)
    end
    empty!(m._objective_cut_ci_saf)

    return nothing
end


function set_first_relax_point!(m::Optimizer)

    if !m._first_relax_point_set
        m._first_relax_point_set = true

        n = m._current_node
        @__dot__ m._current_xref = 0.5*(n.upper_variable_bounds + n.lower_variable_bounds)
        unsafe_check_fill!(isnan, m._current_xref, 0.0, length(m._current_xref))

    end

    return nothing
end
