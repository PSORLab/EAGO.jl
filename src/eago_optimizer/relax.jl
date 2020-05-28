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
$(FUNCTIONNAME)

Relaxs the constraint by adding an affine constraint to the model.
"""
function relax! end

"""
$(FUNCTIONNAME)

Default routine for relaxing quadratic constraint `func` < `0.0` on node `n`. Takes affine bounds of convex part at
point `x0` and secant line bounds on concave parts.
"""
function affine_relax_quadratic!(func::SQF, buffer::Dict{Int,Float64}, saf::SAF, n::NodeBB, x::Vector{Float64}, p1::Bool)

    lower_bounds = n.lower_variable_bounds
    upper_bounds = n.upper_variable_bounds

    quadratic_constant = func.constant

    for term in func.quadratic_terms

        a = p1 ? term.coefficient : -term.coefficient
        idx1 = term.variable_index_1.value
        idx2 = term.variable_index_2.value

        x0_1 = @inbounds x[idx1]
        xL_1 = @inbounds lower_bounds[idx1]
        xU_1 = @inbounds upper_bounds[idx1]

        if idx1 === idx2

            @inbounds buffer[idx1] = (a > 0.0) ? 2.0*a*x0_1 : a*(xL_1 + xU_1)
            quadratic_constant -= (a > 0.0) ? x0_1*x0_1 : a*xL_1*xU_1

        else
            x0_2 = @inbounds x[idx2]
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

function relax!(m::Optimizer, f::BufferedQuadraticIneq, indx::Int, check_safe::Bool)

    affine_relax_quadratic!(f.func, f.buffer, f.saf, m._current_node, m._current_xref, true)
    if check_safe && is_safe_cut!(m, f.saf)
        ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, LT_ZERO)
        push!(m._buffered_quadratic_ineq_ci, ci)
    end
    m.relaxed_to_problem_map[ci] = indx

    return nothing
end

function relax!(m::Optimizer, f::BufferedQuadraticEq, indx::Int, check_safe::Bool)

    affine_relax_quadratic!(f.func, f.buffer, f.saf, m._current_node, m._current_xref, true)
    if check_safe && is_safe_cut!(m, f.saf)
        ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, LT_ZERO)
        push!(m._buffered_quadratic_eq_ci, ci)
    end
    m.relaxed_to_problem_map[ci] = indx

    affine_relax_quadratic!(f.minus_func, f.buffer, f.saf, m._current_node, m._current_xref, false)
    if check_safe && is_safe_cut!(m, f.saf)
        ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, LT_ZERO)
        push!(m._buffered_quadratic_eq_ci, ci)
    end
    m.relaxed_to_problem_map[ci] = indx

    return nothing
end


function bound_objective!(t::ExtensionType, m::Optimizer, x0::Vector{Float64})

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
        obj_indx = wp._objective_sv.variable.value
        objective_lo = @inbounds n.lower_variable_bounds[obj_indx]

    elseif obj_type === SCALAR_AFFINE
        objective_lo = lower_interval_bound(wp._objective_saf_parsed, n)

    elseif obj_type === SCALAR_QUADRATIC
        objective_lo = lower_interval_bound(wp._objective_sqf, n)

    end

    if objective_lo > x._lower_objective_value
        m._lower_objective_value = objective_lo
    end

    return nothing
end
bound_objective!(m::Optimizer, x::Vector{Float64}) = bound_objective!(m.ext_type, m, x)

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
        affine_relax_quadratic!(buffered_sqf.func, buffered_sqf.buffer, buffered_sqf.saf,
                                m._current_node, m._current_xref, true)

        println("wp._objective_saf.terms: $(wp._objective_saf.terms)")
        println("buffered_sqf.saf.terms: $(buffered_sqf.saf.terms)")

        copyto!(wp._objective_saf.terms, buffered_sqf.saf.terms)
        if check_safe && is_safe_cut!(m, wp._objective_saf)
            MOI.set(relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), wp._objective_saf)
        end

    #=
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
            MOI.set(m.relaxed_optimizer, MOI.ConstraintSet(), wp._objective_cut_ci_sv, LT(UBD))

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
            affine_relax_quadratic!(buffered_sqf.func, buffered_sqf.buffer, buffered_sqf.saf,
                                    m._current_node, m._current_xref, true)
            copyto!(wp._objective_saf.terms, buffered_sqf.saf.terms)
            m._objective_saf.constant = buffered_sqf.saf.constant - UBD
            if check_safe && is_safe_cut!(m, wp._objective_saf)
                ci_saf = MOI.add_constraint(m.relaxed_optimizer, wp._objective_saf, LT_ZERO)
                push!(m._objective_cut_ci_saf, ci_saf)
            end
        #=
        elseif obj_type === NONLINEAR
            relax(XXX)
            copyto!(m._objective_saf.terms, XXX)
            m._objective_saf.constant = XXX - UBD
            if check_safe && is_safe_cut!(m, wp._objective_saf)
                ci_saf = MOI.add_constraint(m.relaxed_optimizer, wp._objective_saf, LT_ZERO)
                push!(m._objective_cut_ci_saf, ci_saf)
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
    for i = 1:m._working_problem._saf_leq_count
        sqf_leq = @inbounds sqf_leq_list[i]
        relax!(m, sqf_leq, i, check_safe)
    end

    sqf_eq_list = m._working_problem._sqf_eq
    for i = 1:m._working_problem._saf_eq_count
        sqf_eq = @inbounds sqf_eq_list[i]
        relax!(m, sqf_eq, i, check_safe)
    end

    objective_cut!(m, check_safe)

    return nothing
end
relax_constraints!(t::ExtensionType, m::Optimizer, q::Int64) = relax_all_constraints!(t, m, q)
relax_constraints!(m::Optimizer, q::Int64) = relax_constraints!(m.ext_type, m, q)

"""

Deletes all nonlinear constraints added to the relaxed optimizer.
"""
function delete_nl_constraints!(m::Optimizer)

    # delete affine relaxations added from quadatic inequality
    for ci in m._buffered_quadratic_ineq_ci
        MOI.delete(m.relaxed_optimizer, ci)
    end

    # delete affine relaxations added from quadratic equality
    for ci in m._buffered_quadratic_eq_ci
        MOI.delete(m.relaxed_optimizer, ci)
    end

    return nothing
end

"""
Deletes all scalar-affine objective cuts added to the relaxed optimizer.
"""
function delete_objective_cuts!(m::Optimizer)
    for ci in m._objective_cut_ci_saf
        MOI.delete(m.relaxed_optimizer, ci)
    end
    return nothing
end
