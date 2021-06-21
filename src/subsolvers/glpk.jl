_is_incremental(::Type{GLPK.Optimizer}) = true

function set_default_config!(ext::ExtensionType, d::GlobalOptimizer, m::Incremental{GLPK.Optimizer}, local_solver::Bool)
    return
end