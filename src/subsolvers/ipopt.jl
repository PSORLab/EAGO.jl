_is_incremental(::Type{Ipopt.Optimizer}) = true 

function set_default_config!(ext::ExtensionType, d::GlobalOptimizer, m::Ipopt.Optimizer, local_solver::Bool)
    c_tol = _constraint_tol(d)*1E-3
    MOI.set(m, MOI.RawParameter("tol"),_absolute_tol(d)*1E-3)
    MOI.set(m, MOI.RawParameter("print_level"), 0)
    MOI.set(m, MOI.RawParameter("constr_viol_tol"), c_tol)
    if local_solver
        MOI.set(m, MOI.RawParameter("max_iter"),3000)
        MOI.set(m, MOI.RawParameter("acceptable_tol"), 1E30)
        MOI.set(m, MOI.RawParameter("acceptable_iter"), 300)
        MOI.set(m, MOI.RawParameter("acceptable_compl_inf_tol"), c_tol)
        MOI.set(m, MOI.RawParameter("acceptable_dual_inf_tol"), 1.0)
        MOI.set(m, MOI.RawParameter("acceptable_constr_viol_tol"), c_tol)
    else
        MOI.set(m, MOI.RawParameter("max_iter"), 1E5)
        MOI.set(m, MOI.RawParameter("acceptable_iter"), 1E5+1)
    end
    return
end