using Distributed, DataFrames, CSV
using MINLPLib, JuMP, SCIP

new_lib = "ImprovedCompMult"

# Build library
build_lib = false
instance_names = ["bearing", "ex6_2_10", "ex6_2_10", "ex6_2_11", "ex6_2_12", "ex6_2_13", "ex7_2_4", "ex7_3_1", "ex7_3_2", "ex14_1_8", "ex14_1_9"]
if build_lib
    source_lib = "global"
    for n in instance_names
        MINLPLib.add_to_lib(new_lib, source_lib, n)
    end
end

# Solution handler
struct TrialSummary
    was_error::Bool
    solver_name::String
    instance_name::String
    term_status::MOI.TerminationStatusCode
    primal_status::MOI.ResultStatusCode
    dual_status::MOI.ResultStatusCode
    obj_bound::Float64
    obj_value::Float64
    solve_time::Float64
end
TrialSummary() = TrialSummary(true, "", "", MOI.OTHER_ERROR, MOI.OTHER_RESULT_STATUS, MOI.OTHER_RESULT_STATUS, -Inf, -Inf, 0.0)

function benchmark_problem(d, lib, sname::String)
    m = fetch_model(lib, sname)
    set_optimizer(m, d)

    set_time_limit_sec(m, 1000.0)
    set_silent(m)
    
    JuMP.optimize!(m)
    p = primal_status(m)
    ds = dual_status(m)
    t = termination_status(m)
    v = objective_value(m)
    b = objective_bound(m)
    n = solve_name(m)
    s = solve_time(m)
    return TrialSummary(false,n,sname,t,p,ds,v,b,s)
end


trial_summaries = TrialSummary[]

# Need to assign tolerance here (since tolerances aren't standardized among solvers)
function build_scip()
    m = SCIP.Optimizer(limits_gap=1E-3,    # Absolute tolerance
                       limits_absgap=1E-3  # Relative tolerance
                       )
    return m
end
append!(trial_summaries, pmap(x->benchmark_problem(build_scip, new_lib, x), instance_names; on_error=ex->TrialSummary()))