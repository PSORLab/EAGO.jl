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

function add_variables(m::Optimizer, optimizer::T, variable_number::Int) where T

    variable_index = fill(VI(1), variable_number)
    for i = 1:variable_number
        @inbounds variable_index[i] = MOI.add_variable(optimizer)
        relaxed_variable = SV(@inbounds variable_index[i])
        v_info = @inbounds m._variable_info[i]
        if v_info.is_integer && v_info.is_fixed
            MOI.add_constraint(optimizer, relaxed_variable, ET(v_info.lower_bound))
        elseif v_info.is_integer
            MOI.add_constraint(optimizer, relaxed_variable, ZO())
        elseif v_info.is_fixed
            MOI.add_constraint(optimizer, relaxed_variable, ET(v_info.lower_bound))
        elseif v_info.has_lower_bound && v_info.has_upper_bound
            MOI.add_constraint(optimizer, relaxed_variable, IT(v_info.lower_bound, v_info.upper_bound))
        elseif v_info.has_lower_bound
            MOI.add_constraint(optimizer, relaxed_variable, GT(v_info.lower_bound))
        elseif v_info.has_upper_bound
            MOI.add_constraint(optimizer, relaxed_variable, LT(v_info.upper_bound))
        end
    end

    return variable_index
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
optimize!(::Val{MILP}, m::Optimizer) = optimize!(Val{LP}(), m::Optimizer)

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
optimize!(::Val{MISOCP}, m::Optimizer) = optimize!(Val{SOCP}(), m)

function optimize!(::Val{DIFF_CVX}, m::Optimizer)
    single_nlp_solve!(m)
    return nothing
end
