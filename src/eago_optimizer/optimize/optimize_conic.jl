### LP and MILP routines
function add_soc_constraints!(m::GlobalOptimizer, opt::T) where T

    for (func, set) in m._input_problem._conic_second_order
         MOI.add_constraint(opt, func, set)
    end

    return nothing
end

function optimize!(::SOCP, m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}

    relaxed_optimizer = m.relaxed_optimizer
    MOI.empty!(relaxed_optimizer)

    m._relaxed_variable_index = add_variables(m, relaxed_optimizer, m._input_problem._variable_count)
    add_linear_constraints!(m, relaxed_optimizer)
    add_soc_constraints!(m, relaxed_optimizer)
    add_sv_or_aff_obj!(m, relaxed_optimizer,  m._input_problem._objective)
    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), m._input_problem._optimization_sense)

    if m._parameters.verbosity < 5
        MOI.set(relaxed_optimizer, MOI.Silent(), true)
    end
    m._parse_time = time() - m._start_time

    MOI.optimize!(relaxed_optimizer)

    #unpack_local_solve!(m, relaxed_optimizer)

    return
end

optimize!(::MISOCP, m::GlobalOptimizer) = optimize!(SOCP(), m)