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
function is_safe_cut!(m::GlobalOptimizer, f::SAF)

    safe_l = m._parameters.cut_safe_l
    safe_u = m._parameters.cut_safe_u
    safe_b = m._parameters.cut_safe_b
    (abs(f.constant) > safe_b) && return false # violates |b| <= safe_b

    term_count = length(f.terms)
    @inbounds for i = 1:term_count
        ai = f.terms[i].coefficient
        if ai !== 0.0
            if !(safe_l <= abs(ai) <= safe_u)  # violates safe_l <= abs(ai) <= safe_u
                return false
            end
            @inbounds for j = i:term_count     # violates safe_l <= abs(ai/aj) <= safe_u
                aj = f.terms[j].coefficient
                if aj !== 0.0
                    if !(safe_l <= abs(ai/aj) <= safe_u)
                        return false
                    end
                end
            end
        end
    end
    return true
end

function add_affine_relaxation!(m::GlobalOptimizer{R,S,Q}, f::SAF, check_safe::Bool) where {R,S,Q<:ExtensionType}
    valid_cut_flag = !check_safe || is_safe_cut!(m, f)
    if valid_cut_flag
        s = LT(-f.constant + _constraint_tol(m))
        f.constant = 0.0
        ci = MOI.add_constraint(_relaxed_optimizer(m), f, s)::CI{SAF,LT}
        push!(m._affine_relax_ci, ci)
    end
    return valid_cut_flag
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
function affine_relax_quadratic!(m::GlobalOptimizer, func::SQF, buffer::Dict{Int,Float64}, saf::SAF) where {R,S,Q<:ExtensionType}

    quadratic_constant = func.constant

    # Affine terms only contribute coefficients, so the respective
    # values do not contribute to the cut. Since all quadratic terms
    # are considered to be branch variables we exclude any potential
    # need to retrieve variable bounds from locations other than
    # the node.
    for term in func.quadratic_terms
        a = term.coefficient
        i = term.variable_index_1.value
        j = term.variable_index_2.value
        x0 = _lower_solution(FullVar(), m, i)
        xL = _lower_bound(FullVar(), m, i)
        xU = _upper_bound(FullVar(), m, i)
        if i == j
            if a > 0.0
                buffer[i] += a*x0
                quadratic_constant -= 0.5*a*x0*x0
            else
                if !isinf(xL) && !isinf(xU)
                    buffer[i] += 0.5*a*(xL + xU)
                    quadratic_constant -= 0.5*a*xL*xU
                else
                    return false
                end
            end
        else
            y0 = _lower_solution(FullVar(), m, j)
            yL = _lower_bound(FullVar(), m, j)
            yU = _upper_bound(FullVar(), m, j)
            if a > 0.0
                if (!isinf(xL) && !isinf(yL)) && ((xU - xL)*y0 + (yU - yL)*x0 <= xU*yU - xL*yL)
                    buffer[i] += a*yL
                    buffer[j] += a*xL
                    quadratic_constant -= a*xL*yL

                elseif !isinf(xU) && !isinf(yU)
                    buffer[i] += a*yU
                    buffer[j] += a*xU
                    quadratic_constant -= a*xU*yU
                else
                    return false
                end
            else
                if (!isinf(xU) && !isinf(yL)) && ((xU - xL)*y0 - (yU - yL)*x0 <= xU*yL - xL*yU)
                    buffer[i] += a*yL
                    buffer[j] += a*xU
                    quadratic_constant -= a*xU*yL
                elseif !isinf(xL) && !isinf(yU)
                    buffer[i] += a*yU
                    buffer[j] += a*xL
                    quadratic_constant -= a*xL*yU
                else
                    return false
                end
            end
        end
    end

    for t in func.affine_terms
        buffer[t.variable_index.value] += t.coefficient
    end
    count = 1
    for (key, value) in buffer
        saf.terms[count] = SAT(value, VI(key))
        buffer[key] = 0.0
        count += 1
    end
    saf.constant = quadratic_constant
    return
end

"""
$(TYPEDSIGNATURES)
"""
function relax!(m::GlobalOptimizer, f::BufferedQuadraticIneq, k::Int, check_safe::Bool)
    affine_relax_quadratic!(m, f.func, f.buffer, f.saf)
    valid_cut_flag = add_affine_relaxation!(m, f.saf, check_safe)
    return valid_cut_flag
end

"""
$(TYPEDSIGNATURES)
"""
function relax!(m::GlobalOptimizer, f::BufferedQuadraticEq, indx::Int, check_safe::Bool)

    affine_relax_quadratic!(m, f.func, f.buffer, f.saf)
    valid_cut_flag = add_affine_relaxation!(m, f.saf, check_safe)

    affine_relax_quadratic!(m, f.minus_func, f.buffer, f.saf)
    valid_cut_flag &= add_affine_relaxation!(m, f.saf, check_safe)
    return valid_cut_flag
end

"""
$(TYPEDSIGNATURES)
"""
function check_set_affine_nl!(m::GlobalOptimizer{R,S,Q}, f::BufferedNonlinearFunction{V,N,T}, finite_cut::Bool, check_safe::Bool) where {V,R,S,N,T<:RelaxTag,Q<:ExtensionType}
    valid_cut_flag = finite_cut && (!check_safe || is_safe_cut!(m, f.saf))
    if valid_cut_flag
        lt = LT(_constraint_tol(m) - f.saf.constant)
        f.saf.constant = 0.0
        push!(m._affine_relax_ci, MOI.add_constraint(_relaxed_optimizer(m), f.saf, lt))
    end
    return valid_cut_flag
end

"""
$(TYPEDSIGNATURES)
"""
function relax!(m::GlobalOptimizer{R,S,Q}, f::BufferedNonlinearFunction{V,N,T}, k::Int, check_safe::Bool) where {V,R,S,N,T<:RelaxTag,Q<:ExtensionType}
    d = m._working_problem._relaxed_evaluator
    x = d.variable_values.x
    forward_pass!(d, f)
    valid_cut_flag = true

    grad_sparsity = _sparsity(f)
    if _is_num(f)
        f.saf.constant = _num(f)
        for i = 1:length(grad_sparsity)
            f.saf.terms[i] = SAT(0.0, VI(grad_sparsity[i]))
        end
    else
        v = _set(f)
        if !isempty(v)
            # if has less than or equal to bound (<=)
            if isfinite(_upper_bound(f))
                lower_cut_valid = !isnan(v.cv) && isfinite(v.cv)
                if lower_cut_valid
                    f.saf.constant = v.cv - _upper_bound(f)
                    for (i, k) in enumerate(grad_sparsity)
                        c = v.cv_grad[i]
                        f.saf.terms[i] = SAT(c, VI(k))
                        f.saf.constant -= c*x[k]
                    end
                    valid_cut_flag = check_set_affine_nl!(m, f, lower_cut_valid, check_safe)
                end
            end
            # if has greater than or equal to bound (>=)
            if isfinite(_lower_bound(f))
                upper_cut_valid = !isnan(v.cc) && isfinite(v.cc)
                if upper_cut_valid
                    f.saf.constant = -v.cc +  _lower_bound(f) 
                    for (i, k) in enumerate(grad_sparsity)
                        c = -v.cc_grad[i]
                        f.saf.terms[i] = SAT(c, VI(k))
                        f.saf.constant -= c*x[k]
                    end
                    valid_cut_flag &= check_set_affine_nl!(m, f, upper_cut_valid, check_safe)
                end
            end
        end
    end
    return valid_cut_flag
end

relax!(m::GlobalOptimizer, f::Union{Nothing, SV, AffineFunctionIneq}, k::Int, b::Bool) = true
