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
        v_info = @inbounds m._working_problem._variable_info[i]
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

    # add linear constraints
    for (func, set) in m._input_problem._linear_leq_constraints
         MOI.add_constraint(opt, func, set)
    end
    for (func, set) in m._input_problem._linear_geq_constraints
        MOI.add_constraint(opt, func, set)
    end
    for (func, set) in m._input_problem._linear_eq_constraints
        MOI.add_constraint(opt, func, set)
    end
    nothing
end

### LP and MILP routines
function add_soc_constraints!(m::Optimizer, opt::T) where T

    for (func, set) in m._input_problem._conic_second_order
         MOI.add_constraint(opt, func, set)
    end

    nothing
end

function add_sv_or_aff_obj!(m::Optimizer, opt::T) where T

    if m._input_problem._objective_type === SINGLE_VARIABLE
        MOI.set(opt, MOI.ObjectiveFunction{SV}(), m._input_problem._objective_sv)
    elseif m._input_problem._objective_type === SCALAR_AFFINE
        MOI.set(opt, MOI.ObjectiveFunction{SAF}(), m._input_problem._objective_saf)
    end

    return nothing
end

function unpack_local_solve!(m::Optimizer, opt::T) where T

    m._maximum_node_id = 0

    m._termination_status_code = MOI.get(opt, MOI.TerminationStatus())
    m._result_status_code = MOI.get(opt, MOI.PrimalStatus())
    m._feasible_solution_found = m._result_status_code === MOI.FEASIBLE_POINT

    if MOI.get(opt, MOI.ResultCount()) > 0

        objective_value = MOI.get(opt, MOI.ObjectiveValue())

        # corrects for standard printing multiplier
        if m._input_problem._optimization_sense === MOI.MAX_SENSE
            objective_value *= -1.0
        end

        m._global_lower_bound = objective_value
        m._global_upper_bound = objective_value
        m._objective_value = objective_value
        m._best_upper_value = objective_value
        m._solution_value = objective_value
    end

    m._continuous_solution = zeros(m._input_problem._variable_count)
    for i = 1:m._input_problem._variable_count
         m._continuous_solution[i] = MOI.get(opt, MOI.VariablePrimal(),  m._relaxed_variable_index[i])
    end

    m._run_time = time() - m._start_time

    return nothing
end

function optimize!(::Val{LP}, m::Optimizer)

    relaxed_optimizer = m.relaxed_optimizer
    MOI.empty!(relaxed_optimizer)

    m._relaxed_variable_index = add_variables(m, relaxed_optimizer, m._input_problem._variable_count)
    add_linear_constraints!(m, relaxed_optimizer)
    add_sv_or_aff_obj!(m, relaxed_optimizer)
    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), m._input_problem._optimization_sense)

    if m._parameters.verbosity < 5
        MOI.set(relaxed_optimizer, MOI.Silent(), true)
    end
    m._parse_time = time() - m._start_time

    MOI.optimize!(relaxed_optimizer)

    unpack_local_solve!(m, relaxed_optimizer)

    return nothing
end

optimize!(::Val{MILP}, m::Optimizer) = optimize!(Val{LP}(), m)

function optimize!(::Val{SOCP}, m::Optimizer)

    relaxed_optimizer = m.relaxed_optimizer
    MOI.empty!(relaxed_optimizer)

    m._relaxed_variable_index = add_variables(m, relaxed_optimizer, m._input_problem._variable_count)
    add_linear_constraints!(m, relaxed_optimizer)
    add_soc_constraints!(m, relaxed_optimizer)
    add_sv_or_aff_obj!(m, relaxed_optimizer)
    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), m._input_problem._optimization_sense)

    if m._parameters.verbosity < 5
        MOI.set(relaxed_optimizer, MOI.Silent(), true)
    end
    m._parse_time = time() - m._start_time

    MOI.optimize!(relaxed_optimizer)

    unpack_local_solve!(m, relaxed_optimizer)

    return nothing
end

optimize!(::Val{MISOCP}, m::Optimizer) = optimize!(Val{SOCP}(), m)

function optimize!(::Val{DIFF_CVX}, m::Optimizer)

    single_nlp_solve!(m)

    return nothing
end
