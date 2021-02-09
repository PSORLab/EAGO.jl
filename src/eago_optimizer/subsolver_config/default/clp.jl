
function set_default_config!(::Val{:clp}, ext::ExtensionType, m::T) where T <: MOI.AbstractOptimizer

    MOI.set(m, MOI.RawParameter("PrimalTolerance"),       1E-7)
    MOI.set(m, MOI.RawParameter("DualTolerance"),         1E-7)
    MOI.set(m, MOI.RawParameter("DualObjectiveLimit"),    1e308)
    MOI.set(m, MOI.RawParameter("MaximumIterations"),     2147483647)
    MOI.set(m, MOI.RawParameter("PresolveType"),          0)
    MOI.set(m, MOI.RawParameter("SolveType"),             5)
    MOI.set(m, MOI.RawParameter("InfeasibleReturn"),      1)
    MOI.set(m, MOI.RawParameter("Scaling"),               3)
    MOI.set(m, MOI.RawParameter("Perturbation"),          100)

    return nothing
end
