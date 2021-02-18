
#=
Base.@kwdef mutable struct SIPOptimizer <: MOI.AbstractOptimizer
    objective_sense::MOI.OptimizationSense = MOI.FEASIBLE_SENSE
    time_limit::Float64 = 1000.0
    is_silent::Bool = false
    x::Vector{Float64} = Float64[]
    obj_bound::Float64 = 0.0
    obj_value::Float64 = 0.0
    termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    result_status_code::MOI.ResultStatusCode = MOI.UNKNOWN_RESULT_STATUS
    run_time::Float64 = 0.0
end

function MOI.empty!(m::SIPOptimizer)
end
function MOI.is_empty(m::Optimizer)
end

function MOI.set(m::SIPOptimizer, ::MOI.Silent, value)
    m.is_silent = value
    return nothing
end

function MOI.set(m::SIPOptimizer, ::MOI.TimeLimitSec, value::Nothing)
    m.time_limit = Inf
    return nothing
end
function MOI.set(m::SIPOptimizer, ::MOI.TimeLimitSec, value::Float64)
    m.time_limit = value
    return nothing
end

MOI.set(m::SIPOptimizer, p::MOI.RawParameter, value) = error("The SIPOptimizer does not have any raw parameters.")

function MOI.get(m::SIPOptimizer, ::MOI.ListOfVariableIndices)
end
function MOI.get(m::SIPOptimizer, ::MOI.NumberOfVariables)
end

function MOI.get(m::SIPOptimizer, ::MOI.ObjectiveValue)

end
function MOI.get(m::SIPOptimizer, ::MOI.ObjectiveBound)

end

MOI.get(m::SIPOptimizer, p::MOI.RawParameter) = error("The SIPOptimizer does not have any raw parameters.")


function MOI.get(m::SIPOptimizer, ::MOI.RelativeGap)
    abs(objective_bound() - objective_value())/objective_value()
end
MOI.get(m::SIPOptimizer, ::MOI.SolverName) = "EAGO: SIPOptimizer"
MOI.get(m::SIPOptimizer, ::MOI.TerminationStatus) = m._termination_status_code
MOI.get(m::SIPOptimizer, ::MOI.PrimalStatus) = m.result_status_code
MOI.get(m::SIPOptimizer, ::MOI.SolveTime) = m.run_time
MOI.get(m::SIPOptimizer, ::MOI.ResultCount) = (m.result_status_code === MOI.FEASIBLE_POINT) ? 1 : 0
MOI.get(m::SIPOptimizer, ::MOI.TimeLimitSec) = m.time_limit

MOI.supports(::SIPOptimizer, ::MOI.TimeLimitSec) = true
MOI.supports(::SIPOptimizer, ::MOI.ObjectiveSense) = true
function MOI.set(m::SIPOptimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    m.objective_sense = sense
    return nothing
end


function MOI.optimize!(m::SIPOptimizer)
end
=#
