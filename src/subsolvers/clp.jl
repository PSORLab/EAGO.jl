_is_incremental(::Type{Clp.Optimizer}) = false

function set_default_config!(ext::ExtensionType, d::GlobalOptimizer, m::Clp.Optimizer, local_solver::Bool)
    MOI.set(m, MOI.RawParameter("PrimalTolerance"),       _absolute_tol(d)*1E-2)
    MOI.set(m, MOI.RawParameter("DualTolerance"),         _absolute_tol(d)*1E-2)
    MOI.set(m, MOI.RawParameter("DualObjectiveLimit"),    1e308)
    MOI.set(m, MOI.RawParameter("MaximumIterations"),     2147483647)
    MOI.set(m, MOI.RawParameter("PresolveType"),          0)
    MOI.set(m, MOI.RawParameter("SolveType"),             5)
    MOI.set(m, MOI.RawParameter("InfeasibleReturn"),      1)
    MOI.set(m, MOI.RawParameter("Scaling"),               3)
    MOI.set(m, MOI.RawParameter("Perturbation"),          100)
    return
end