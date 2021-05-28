
"""

A type-stable cache used to wrapper for an optimizer that enables incremental
modification of solvers that don't inherently suppport this. Explicitly checks
support of MOI functionality used in EAGO.
"""
mutable struct Incremental{Q,S} <: MOI.AbstractOptimizer
    optimizer::S
    cache::MOIB.LazyBridgeOptimizer{MOIU.GenericModel{Float64,MOIU.ModelFunctionConstraints{Float64}}}
end
function Incremental(m::S) where S
    is_incremental = MOI.supports_incremental_interface(m, false)
    cache = MOIB.full_bridge_optimizer(MOIU.Model{Float64}(), Float64)
    b = MOIB.full_bridge_optimizer(m, Float64)
    return Incremental{Val(is_incremental),typeof(b)}(b, cache)
end

_get_storage(d::Incremental{Val{:true},S}) where S <: MOI.AbstractOptimizer = d.optimizer
_get_storage(d::Incremental{Val{:false},S}) where S <: MOI.AbstractOptimizer = d.cache

# Set attributes
for F in (SAF, SV)
    @eval function MOI.set(d::Incremental{Q,S}, ::MOI.ObjectiveFunction{F}, f::F) where {Q, S <: MOI.AbstractOptimizer}
        MOI.set(_get_storage(d), MOI.ObjectiveFunction{F}(), f)
    end
end
function MOI.set(d::Incremental{Q,S}, ::MOI.ObjectiveSense, s) where {Q, S <: MOI.AbstractOptimizer}
    MOI.set(_get_storage(d), MOI.ObjectiveSense(), s)
end

# Add variable/constraint
function MOI.add_variable(d::Incremental{Q,S}) where {Q, S <: MOI.AbstractOptimizer}
    MOI.add_variable(_get_storage(d))
end
function MOI.add_constraint(d::Incremental{Q,S}, f::SV, s::Union{LT,GT,ET}) where {Q, S <: MOI.AbstractOptimizer}
    MOI.add_constraint(_get_storage(d), f, s)
end

function MOI.add_constraint(d::Incremental{Q,S}, f::SAF, s::LT) where {Q, S <: MOI.AbstractOptimizer}
    MOI.add_constraint(_get_storage(d), f, s)
end

# Delete
function MOI.delete(d::Incremental{Q,S}, ci::CI{SAF,LT}) where {Q, S <: MOI.AbstractOptimizer}
     MOI.delete(_get_storage(d), ci)
end

# Set modifications
function MOI.set(d::Incremental{Q,S}, ::MOI.ConstraintSet, ci::CI{SV,T}, s::T) where  {T <: Union{LT,GT,ET}, Q, S <: MOI.AbstractOptimizer}
     MOI.set(_get_storage(d), MOI.ConstraintSet(), ci, s)
end

# Get attributes
for attr in (MOI.TerminationStatus, MOI.PrimalStatus, MOI.DualStatus, MOI.ObjectiveValue,
             MOI.DualObjectiveValue, MOI.ListOfConstraintIndices{SAF,LT})
    @eval function MOI.get(d::Incremental{Q,T}, ::$attr) where {Q,T}
        MOI.get(d.optimizer, $attr())
    end
end
function MOI.get(d::Incremental{Q,T}, ::MOI.VariablePrimal, vi::VI) where {Q,T}
    MOI.get(d.optimizer, MOI.VariablePrimal(), vi)
end
function MOI.get(d::Incremental{Q,T}, ::MOI.ConstraintDual, ci::Union{CI{SV,LT},CI{SV,GT}}) where {Q,T}
    MOI.get(d.optimizer, MOI.ConstraintDual(), ci)
end

# define optimize!
function MOI.optimize!(d::Incremental{Val{:true},S}) where S <: MOI.AbstractOptimizer
    MOI.optimize!(d.optimizer)
    return
end
function MOI.optimize!(d::Incremental{Val{:false},S}) where S <: MOI.AbstractOptimizer
    MOIU.copy_to(d.optimizer, d.cache)
    MOI.optimize!(d.optimizer)
    return
end
