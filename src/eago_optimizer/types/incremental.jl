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
support of MOI functionality used in EAGO. For `Q = Val{true}`, the subsolver
supports incremental loading. For `Q = Val{false}`, the subsolver does not.
"""
mutable struct Incremental{Q,S <: MOI.AbstractOptimizer} <: MOI.AbstractOptimizer
    optimizer::MOIB.LazyBridgeOptimizer{S}
    cache::MOIB.LazyBridgeOptimizer{MOIU.CachingOptimizer{S,MOIU.GenericModel{Float64,MOIU.ModelFunctionConstraints{Float64}}}}
end
function Incremental(m::S) where S <: MOI.AbstractOptimizer
    is_incremental = MOIU.supports_default_copy_to(m, false)
    b = MOIB.full_bridge_optimizer(m, Float64)
    cache = MOIB.full_bridge_optimizer(MOIU.CachingOptimizer(MOIU.Model{Float64}(), m), Float64)
    return Incremental{Val{is_incremental},typeof(m)}(b, cache)
end

_get_storage(d::Incremental{Val{true},S}) where S = d.optimizer
_get_storage(d::Incremental{Val{false},S}) where S = d.cache

# Set attributes
for F in (SV, SAF, SQF)
    @eval function MOI.set(d::Incremental{Q,S}, ::MOI.ObjectiveFunction{$F}, f::$F) where {Q, S}
        MOI.set(_get_storage(d), MOI.ObjectiveFunction{$F}(), f)
        return
    end
end
function MOI.set(d::Incremental{Q,S}, ::MOI.ObjectiveSense, s) where {Q, S}
    MOI.set(_get_storage(d), MOI.ObjectiveSense(), s)
    return
end
function MOI.set(d::Incremental{Q,S}, ::MOI.Silent, s) where {Q, S}
    MOI.set(_get_storage(d), MOI.Silent(), s)
    return
end

function MOI.set(d::Incremental{Q,S}, ::MOI.VariablePrimalStart, v::VI, x) where {Q, S}
    MOI.set(_get_storage(d), MOI.VariablePrimalStart(), v, x)
    return
end

function MOI.set(d::Incremental{Q,S}, ::MOI.NLPBlock, s) where {Q, S}
    MOI.set(_get_storage(d), MOI.NLPBlock(), s)
    return
end

function MOI.set(d::Incremental{Q,S}, p::MOI.RawParameter, s) where {Q, S}
    MOI.set(_get_storage(d), p, s)
    return
end

# Add variable/constraint
function MOI.add_variable(d::Incremental{Q,S}) where {Q, S}
    MOI.add_variable(_get_storage(d))::VI
end
function MOI.add_constraint(d::Incremental{Q,S}, f::SV, s::T) where {Q, S, T<:Union{LT,GT,ET,IT}}
    MOI.add_constraint(_get_storage(d), f, s)::CI{SV,T}
end

function MOI.add_constraint(d::Incremental{Q,S}, f::R, s::T) where {Q, S, R<:Union{SAF,SQF}, T<:Union{LT,GT,ET,IT}}
    MOI.add_constraint(_get_storage(d), f, s)::CI{R,T}
end

# Delete
function MOI.delete(d::Incremental{Q,S}, ci::CI{SAF,LT}) where {Q, S}
     MOI.delete(_get_storage(d), ci)
     return
end

# Set modifications
function MOI.set(d::Incremental{Q,S}, ::MOI.ConstraintSet, ci::CI{SV,T}, s::T) where  {T <: Union{LT,GT,ET,IT}, Q, S}
     MOI.set(_get_storage(d), MOI.ConstraintSet(), ci, s)
     return
end

# Get attributes
function MOI.get(d::Incremental{Q,T}, ::MOI.TerminationStatus) where {Q,T}
    MOI.get(d.optimizer, MOI.TerminationStatus())::MOI.TerminationStatusCode
end
function MOI.get(d::Incremental{Q,T}, ::MOI.PrimalStatus) where {Q,T}
    MOI.get(d.optimizer, MOI.PrimalStatus())::MOI.ResultStatusCode
end
function MOI.get(d::Incremental{Q,T}, ::MOI.DualStatus) where {Q,T}
    MOI.get(d.optimizer, MOI.DualStatus())::MOI.ResultStatusCode
end
function MOI.get(d::Incremental{Q,T}, ::MOI.ObjectiveValue) where {Q,T}
    MOI.get(d.optimizer, MOI.ObjectiveValue())::Float64
end
function MOI.get(d::Incremental{Q,T}, ::MOI.DualObjectiveValue) where {Q,T}
    MOI.get(d.optimizer, MOI.DualObjectiveValue())::Float64
end

function MOI.get(d::Incremental{Q,T}, ::MOI.VariablePrimal, vi::VI) where {Q,T}
    MOI.get(d.optimizer, MOI.VariablePrimal(), vi)::Float64
end
function MOI.get(d::Incremental{Q,T}, ::MOI.ConstraintDual, ci::Union{CI{SV,LT},CI{SV,GT}}) where {Q,T}
    MOI.get(d.optimizer, MOI.ConstraintDual(), ci)::Float64
end
function MOI.get(d::Incremental{Q,T}, ::MOI.ResultCount) where {Q,T}
    MOI.get(d.optimizer, MOI.ResultCount())::Int
end
function MOI.get(d::Incremental{Q,S}, n::MOI.SolverName) where {Q, S}
    MOI.get(d.optimizer, n)::String
end

# define optimize!
function MOI.optimize!(d::Incremental{Val{:true},S}) where S
    MOI.optimize!(d.optimizer)
    return
end
function MOI.optimize!(d::Incremental{Val{:false},S}) where S
    MOI.optimize!(d.cache)
    return
end

function MOI.empty!(d::Incremental{Q,S}) where {Q, S}
    MOI.empty!(_get_storage(d))
    return
end
