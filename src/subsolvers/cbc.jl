_is_incremental(::Type{Cbc.Optimizer}) = false

function set_default_config!(ext::ExtensionType, d::GlobalOptimizer, m::Cbc.Optimizer, local_solver::Bool)
    MOI.set(m, MOI.RawParameter("allowableGap"), _absolute_tol(d)*1E-2)
    MOI.set(m, MOI.RawParameter("ratioGap"),     _absolute_tol(d)*1E-2)
    MOI.set(m, MOI.RawParameter("logLevel"),     0)
    #MOI.set(m, MOI.RawParameter("threads"), Threads.nthreads())
    return
end