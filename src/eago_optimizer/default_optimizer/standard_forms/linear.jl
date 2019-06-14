"""
    relax_linear!

Adds linear constraints to `trg` optimizer.
"""
function relax_linear!(src::Optimizer,trg::T) where {T<:MOI.AbstractOptimizer}
    for (func,set,ind) in src.linear_leq_constraints
        MOI.add_constraint(trg, func, set)
    end
    for (func,set,ind) in src.linear_geq_constraints
        MOI.add_constraint(trg, func, set)
    end
    for (func,set,ind) in src.linear_eq_constraints
        MOI.add_constraint(trg, func, set)
    end
    for (func,set,ind) in src.linear_interval_constraints
        MOI.add_constraint(trg, func, MOI.GreaterThan{Float64}(set.lower))
        MOI.add_constraint(trg, func, MOI.LessThan{Float64}(set.upper))
    end

    if isa(src.objective, MOI.SingleVariable)
        if (src.optimization_sense == MOI.MIN_SENSE)
            MOI.set(trg, MOI.ObjectiveFunction{MOI.SingleVariable}(), src.objective)
        elseif (src.optimization_sense == MOI.MAX_SENSE)
            neg_obj_var = MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm{Float64}(-1.0, src.objective.variable)], 0.0)
            MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), neg_obj_var)
        else
            error("Objective sense must be MOI.MinSense or MOI.MaxSense")
        end
    elseif isa(src.objective, MOI.ScalarAffineFunction{Float64})
        if (src.optimization_sense == MOI.MIN_SENSE)
            MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), src.objective)
        elseif (src.optimization_sense == MOI.MAX_SENSE)
            neg_obj_aff_terms = []
            for term in src.objective.terms
                push!(neg_obj_aff_terms,MOI.ScalarAffineTerm{Float64}(-term.coefficient,term.variable_index))
            end
            neg_obj_aff = MOI.ScalarAffineFunction{Float64}(neg_obj_aff_terms, -src.objective.constant)
            MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), neg_obj_aff)
        else
            error("Objective sense must be MOI.MIN_SENSE or MOI.MAX_SENSE")
        end
    end
end

"""
    objective_cut_linear!

Adds linear objective cut constraint to `trg` optimizer.
"""
function objective_cut_linear!(src::Optimizer,trg) where {T<:MOI.AbstractOptimizer}
    if isa(src.objective, MOI.SingleVariable)
        obj_set = MOI.LessThan(src.global_upper_bound)
        if (src.optimization_sense == MOI.MIN_SENSE)
            MOI.add_constraint(trg, src.objective, obj_set)
        elseif (src.optimization_sense == MOI.MAX_SENSE)
            neg_obj_var = MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm{Float64}(-1.0, src.objective.variable)], 0.0)
            MOI.add_constraint(trg, neg_obj_var, obj_set)
        else
            error("Objective sense must be MOI.MinSense or MOI.MaxSense")
        end
    elseif isa(src.objective, MOI.ScalarAffineFunction{Float64})
        obj_set = MOI.LessThan(src.global_upper_bound)
        if (src.optimization_sense == MOI.MIN_SENSE)
            MOI.add_constraint(trg, src.objective, obj_set)
        elseif (src.optimization_sense == MOI.MAX_SENSE)
            neg_obj_aff_terms = []
            for term in src.objective.terms
                push!(neg_obj_aff_terms,MOI.ScalarAffineTerm{Float64}(-term.coefficient,term.variable_index))
            end
            neg_obj_aff = MOI.ScalarAffineFunction{Float64}(neg_obj_aff_terms, -src.objective.constant)
            MOI.add_constraint(trg, neg_obj_aff, obj_set)
        else
            error("Objective sense must be MOI.MinSense or MOI.MaxSense")
        end
    end
end
