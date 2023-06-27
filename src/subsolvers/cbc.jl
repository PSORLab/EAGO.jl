_is_incremental(::Type{Cbc.Optimizer}) = false

function set_default_config!(ext::ExtensionType, d::GlobalOptimizer, m::Incremental{Cbc.Optimizer}, local_solver::Bool)
    #MOI.set(m, MOI.RawOptimizerAttribute("allowableGap"), _absolute_tol(d)*1E-2)
    #MOI.set(m, MOI.RawOptimizerAttribute("ratioGap"),     _absolute_tol(d)*1E-2)
    #MOI.set(m, MOI.RawOptimizerAttribute("logLevel"),     0)
    #MOI.set(m, MOI.RawOptimizerAttribute("threads"), Threads.nthreads())
    return
end

function capture_stdout(f)
    open("nul", "w") do out
        redirect_stdout(out) do
            f()        
        end
    end
end

# Used to prevent Cbc from printing to console even when set to silent on continuous problems. 
function MOI.optimize!(::Type{Cbc.Optimizer}, d::Incremental{Cbc.Optimizer})
    continuous_problem = true
    for (F,S) in MOI.get(d, MOI.ListOfConstraintTypesPresent())
        if S isa Union{INT, MOI.SOS1, MOI.SOS2, ZO}
            continuous_problem = false
            break
        end 
    end
    if continuous_problem
        capture_stdout() do
            MOI.optimize!(_get_storage(d))
        end
    else
        MOI.optimize!(_get_storage(d))
    end
end