_is_incremental(::Type{IpoptMathOptInterfaceExt.Optimizer}) = true 

function set_default_config!(ext::ExtensionType, d::GlobalOptimizer, m::Incremental{IpoptMathOptInterfaceExt.Optimizer}, local_solver::Bool)
    c_tol = _constraint_tol(d)*1E-3
    MOI.set(m, MOI.RawOptimizerAttribute("tol"), _absolute_tol(d)*1E-3)
    MOI.set(m, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.set(m, MOI.RawOptimizerAttribute("constr_viol_tol"), c_tol)
    if local_solver
        MOI.set(m, MOI.RawOptimizerAttribute("max_iter"), 10000)
        MOI.set(m, MOI.RawOptimizerAttribute("acceptable_tol"), 1E30)
        MOI.set(m, MOI.RawOptimizerAttribute("acceptable_iter"), 1000)
        MOI.set(m, MOI.RawOptimizerAttribute("acceptable_compl_inf_tol"), c_tol)
        MOI.set(m, MOI.RawOptimizerAttribute("acceptable_dual_inf_tol"), 1.0)
        MOI.set(m, MOI.RawOptimizerAttribute("acceptable_constr_viol_tol"), c_tol)
    else
        MOI.set(m, MOI.RawOptimizerAttribute("max_iter"), 1E5)
        MOI.set(m, MOI.RawOptimizerAttribute("acceptable_iter"), 1E5+1)
    end
    return
end