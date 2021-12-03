_is_incremental(::Type{Clp.Optimizer}) = false

function set_default_config!(ext::ExtensionType, d::GlobalOptimizer, m::Clp.Optimizer, local_solver::Bool)
    MOI.set(m, MOI.RawOptimizerAttribute("PrimalTolerance"),       _absolute_tol(d)*1E-2)
    MOI.set(m, MOI.RawOptimizerAttribute("DualTolerance"),         _absolute_tol(d)*1E-2)
    MOI.set(m, MOI.RawOptimizerAttribute("DualObjectiveLimit"),    1e308)
    MOI.set(m, MOI.RawOptimizerAttribute("MaximumIterations"),     2147483647)
    MOI.set(m, MOI.RawOptimizerAttribute("PresolveType"),          0)
    MOI.set(m, MOI.RawOptimizerAttribute("SolveType"),             5)
    MOI.set(m, MOI.RawOptimizerAttribute("InfeasibleReturn"),      1)
    MOI.set(m, MOI.RawOptimizerAttribute("Scaling"),               3)
    MOI.set(m, MOI.RawOptimizerAttribute("Perturbation"),          100)
    return
end