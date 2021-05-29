# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/types/incremental.jl
# A type-stable wrapper for optimizers used by EAGO to enable bridging and
# incremental loading.
#############################################################################

"""
    Incremental{Q,S}

A type-stable cache used to wrapper for an optimizer that enables incremental
modification of solvers that don't inherently suppport this. Explicitly checks
support of MOI functionality used in EAGO.
"""
mutable struct Incremental{Q,S} <: MOI.AbstractOptimizer
    optimizer::S
    cache::MOIB.LazyBridgeOptimizer{MOIU.GenericModel{Float64,MOIU.ModelFunctionConstraints{Float64}}}
end
function Incremental(m::S) where S
    is_incremental = MOIU.supports_default_copy_to(m, false)
    cache = MOIB.full_bridge_optimizer(MOIU.Model{Float64}(), Float64)
    b = MOIB.full_bridge_optimizer(m, Float64)
    return Incremental{Val{is_incremental},typeof(b)}(b, cache)
end

_get_storage(d::Incremental{Val{true},S}) where S = d.optimizer
_get_storage(d::Incremental{Val{false},S}) where S = d.cache

# Set attributes
for F in (SV, SAF, SQF)
    @eval function MOI.set(d::Incremental{Q,S}, ::MOI.ObjectiveFunction{$F}, f::$F) where {Q, S}
        MOI.set(_get_storage(d), MOI.ObjectiveFunction{$F}(), f)
    end
end
function MOI.set(d::Incremental{Q,S}, ::MOI.ObjectiveSense, s) where {Q, S}
    MOI.set(_get_storage(d), MOI.ObjectiveSense(), s)
end
function MOI.set(d::Incremental{Q,S}, ::MOI.VariablePrimalStart, v::VI, x) where {Q, S}
    MOI.set(_get_storage(d), MOI.VariablePrimalStart(), v, x)
end

function MOI.set(d::Incremental{Q,S}, ::MOI.NLPBlock, x) where {Q, S}
    MOI.set(_get_storage(d), MOI.NLPBlock(), x)
end

# Add variable/constraint
function MOI.add_variable(d::Incremental{Q,S}) where {Q, S}
    MOI.add_variable(_get_storage(d))
end
function MOI.add_constraint(d::Incremental{Q,S}, f::SV, s::Union{LT,GT,ET}) where {Q, S}
    MOI.add_constraint(_get_storage(d), f, s)
end

function MOI.add_constraint(d::Incremental{Q,S}, f::Union{SAF,SQF}, s::Union{LT,GT,ET}) where {Q, S}
    MOI.add_constraint(_get_storage(d), f, s)
end

# Delete
function MOI.delete(d::Incremental{Q,S}, ci::CI{SAF,LT}) where {Q, S}
     MOI.delete(_get_storage(d), ci)
end

# Set modifications
function MOI.set(d::Incremental{Q,S}, ::MOI.ConstraintSet, ci::CI{SV,T}, s::T) where  {T <: Union{LT,GT,ET}, Q, S}
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
function MOI.optimize!(d::Incremental{Val{:true},S}) where S
    MOI.optimize!(d.optimizer)
    return
end
function MOI.optimize!(d::Incremental{Val{:false},S}) where S
    MOIU.copy_to(d.optimizer, d.cache)
    MOI.optimize!(d.optimizer)
    return
end

MOI.empty!(d::Incremental{Q,S}) where {Q, S} = MOI.empty!(_get_storage(d))
