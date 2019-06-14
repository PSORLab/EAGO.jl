#=
MOI.canget(::EAGOOptimizer, ::MOI.NumberOfVariables) = true
MOI.canget(::EAGOOptimizer, ::MOI.ListOfVariableIndices) = true
MOI.canget(::EAGOOptimizer, ::MOI.SolverName) = true
MOI.canget(m::EAGOOptimizer, ::MOI.TerminationStatus) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.ObjectiveValue) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.ObjectiveBound) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.RelativeGap) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.SolveTime) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.NodeCount) = m.started_solve
=#

MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices) = [MOI.VariableIndex(i) for i in 1:length(m.variable_info)]
MOI.get(m::Optimizer, ::MOI.ObjectiveValue) = m.objective_value
MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = length(m.variable_info)
MOI.get(m::Optimizer, ::MOI.ObjectiveBound) = m.global_upper_bound
MOI.get(m::Optimizer, ::MOI.RelativeGap) = m.global_upper_bound
MOI.get(m::Optimizer, ::MOI.SolverName) = "EAGO: Easy Advanced Global Optimization"
MOI.get(m::Optimizer, ::MOI.TerminationStatus) = m.termination_status_code
MOI.get(m::Optimizer, ::MOI.PrimalStatus) = m.result_status_code

function MOI.get(m::Optimizer, ::MOI.SolveTime)
    if m.current_iteration_count > 0
        return m.history.preprocess_time[m.current_iteration_count-1] +
               m.history.postprocess_time[m.current_iteration_count-1] +
               m.history.lower_time[m.current_iteration_count-1] +
               m.history.upper_time[m.current_iteration_count-1]
     else
         return 0.0
     end
 end
#MOI.get(m::Optimizer, ::MOI.NodeCount) = length(m.MaximumNodeID)

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    check_inbounds(model, vi)
    return model.continuous_solution[vi.value]
end

get_history(x::Optimizer) = x.history
function get_history(x::JuMP.Model)
    inner_optimizer = x.moi_backend.optimizer.model.optimizer
    if isa(inner_optimizer,EAGO.Optimizer)
        return get_history(inner_optimizer)
    else
        error("get_history(x::JuMP.Model) requires that the JuMP model use a EAGO.Optimizer.")
    end
end
