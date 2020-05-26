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
# src/eago_optimizer/problem_types.jl
# Defines constraints supported Optimizer and how to store them.
#############################################################################

#=
LP          -> COPY TO RELAXED SOLVER AND SOLVE
MILP        -> COPY TO RELAXED SOLVER AND SOLVE
SOCP        -> COPY TO RELAXED SOLVER AND SOLVE
MISOCP      -> COPY TO RELAXED SOLVER AND SOLVE
DIFF_CVX    -> COPY TO NLP SOLVER AND SOLVE (POTENTIAL MULTISTART)
NS_CVX      -> COPY TO NLP SOLVER AND SOLVE (POTENTIAL MULTISTART)
DIFF_NCVX   -> APPLY GLOBAL SOLVER (UNLESS USER REQUEST LOCAL SOLVE THEN NLP)
NS_NCVX     -> APPLY GLOBAL SOLVER (UNLESS USER REQUEST LOCAL SOLVE THEN NLP)
MINCVX      -> APPLY GLOBAL SOLVER (LOCAL SOLVE OPTION FUTURE FEATURE)
=#

function add_variables(m::Optimizer, opt::T, nvar::Int) where T

    variable_index = fill(VI(1), nvar)
    for i = 1:nvar
        @inbounds variable_index[i] = MOI.add_variable(opt)
        relaxed_variable = SV(@inbounds variable_index[i])
        v_info = @inbounds m._variable_info[i]
        if v_info.is_integer && v_info.is_fixed
            MOI.add_constraint(opt, relaxed_variable, ET(v_info.lower_bound))
        elseif v_info.is_integer
            MOI.add_constraint(opt, relaxed_variable, ZO())
        elseif v_info.is_fixed
            MOI.add_constraint(opt, relaxed_variable, ET(v_info.lower_bound))
        elseif v_info.has_lower_bound && v_info.has_upper_bound
            MOI.add_constraint(opt, relaxed_variable, IT(v_info.lower_bound, v_info.upper_bound))
        elseif v_info.has_lower_bound
            MOI.add_constraint(opt, relaxed_variable, GT(v_info.lower_bound))
        elseif v_info.has_upper_bound
            MOI.add_constraint(opt, relaxed_variable, LT(v_info.upper_bound))
        end
    end
    variable_index
end

### LP and MILP routines
function add_linear_constraints!(m::Optimizer, opt::T) where T
    for (id, cinfo) in m._linear_constraint
        MOI.add_constraint(opt, cinfo.func, cinfo.set)
    end
    nothing
end

function optimize!(::Val{LP}, m::Optimizer)

    relaxed_optimizer = m.relaxed_optimizer()
    m._relaxed_variable_index = add_variables(m, relaxed_optimizer, m._input_variable_number)
    add_linear_constraints!(m, relaxed_optimizer)
    add_sv_obj_relax!(m, relaxed_optimizer)

    (m.verbosity < 5) && MOI.set(relaxed_optimizer, MOI.Silent(), true)
    MOI.optimize!(relaxed_optimizer)

    unpack_local_solve!(m, relaxed_optimizer)
end
#optimize!(::Val{MILP}, m::Optimizer) = optimize!(Val{LP}(), m::Optimizer)

function optimize!(::Val{SOCP}, m::Optimizer)

    relaxed_optimizer = m.relaxed_optimizer()
    m._relaxed_variable_index = add_variables(m, relaxed_optimizer, m._input_variable_number)
    add_linear_constraints!(m, relaxed_optimizer)
    add_soc_constraints!(m, relaxed_optimizer)
    add_sv_obj_relax!(m, relaxed_optimizer)

    (m.verbosity < 5) && MOI.set(relaxed_optimizer, MOI.Silent(), true)
    MOI.optimize!(relaxed_optimizer)

    unpack_local_solve!(m, relaxed_optimizer)
end
#optimize!(::Val{MISOCP}, m::Optimizer) = optimize!(Val{SOCP}(), m)

function optimize!(::Val{DIFF_CVX}, m::Optimizer)
    single_nlp_solve!(m)
    nothing
end
optimize!(::Val{MINCVX}, m::Optimizer) = global_solve!(m)

### Routines for parsing the full nonconvex problem
"""
Reformulates quadratic terms in SOC constraints if possible. For <= or >=,
the quadratic term is deleted if an SOCP is detected. For ==, the SOC check
is done for each >= and <=, the convex constraint is reformulated to a SOC,
the concave constraint is keep as a quadratic.
"""
function parse_classify_quadratic!(m::Optimizer)
    #=
    for (id, cinfo) in m._quadratic_constraint
        is_soc, add_concave, cfunc, cset, qfunc, qset = check_convexity(cinfo.func, cinfo.set)
        if is_soc
            MOI.add_constraint(m, cfunc, cset)
            deleteat!(m._quadratic_constraint, id)
            if add_concave
                MOI.add_constraint(m, qfunc, qset)
            end
        end
    end
    =#
    nothing
end

"""
"""
function parse_classify_nlp(m)
    nothing
end

"""
Classifies the problem type
"""
function parse_classify_problem!(m::Optimizer)

    ip = m._input_problem
    integer_variable_number = count(is_integer.(ip._variable_info))
    cone_constraint_number = ip._conic_second_order_count
    quad_constraint_number = ip._quadratic_leq_count + ip._quadratic_geq_count + ip._quadratic_eq_count

    if integer_variable_number === 0
        if cone_constraint_number === 0 && quad_constraint_number === 0
            # && iszero(m._input_nonlinear_constraint_number)
            m._problem_type = LP
        elseif quad_constraint_number === 0
            # && iszero(m._input_nonlinear_constraint_number)
            m._problem_type = SOCP
        else
            #parse_classify_quadratic!(m)
            #if iszero(m._input_nonlinear_constraint_number)
            #    if isempty(m._quadratic_constraint)
            #        m._problem_type = SOCP
            #    end
            #else
            #    # Check if DIFF_CVX, NS_CVX, DIFF_NCVX, OR NS_NCVX
            #    m._problem_type = parse_classify_nlp(m)
            #end
            m._problem_type = NS_NCVX
        end
    else
        if cone_constraint_number === 0 && quad_constraint_number === 0
            m._problem_type = MILP
        elseif quad_constraint_number === 0
            m._problem_type = MISOCP
        else
            parse_classify_quadratic!(m)
            #=
            if iszero(m._nonlinear_constraint_number)
                if iszero(m._quadratic_constraint_number)
                    m._problem_type = MISOCP
                end
            else
                # Performs parsing
                _ = parse_classify_nlp(m)
            end
            =#
            m._problem_type = MINCVX
        end
    end

    return nothing
end
