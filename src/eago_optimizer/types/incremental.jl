# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/types/incremental.jl
# A type-stable wrapper with for optimizers used by EAGO to enable bridging and
# incremental loading. This is taylored to the internal routines used by EAGO.jl
# so methods may be specialized by optimizer types and error checking is often
# avoided.
#############################################################################

#=
mutable struct IncrementalCache{S <: MOI.AbstractOptimizer} <: MOI.AbstractOptimizer}
end
=#

"""
    Incremental{Q,S}

A type-stable cache used to wrapper for an optimizer that enables incremental
modification of solvers that don't inherently suppport this. Explicitly checks
support of MOI functionality used in EAGO. For `Q = Val{true}`, the subsolver
supports incremental loading. For `Q = Val{false}`, the subsolver does not.
"""
mutable struct Incremental{S <: MOI.AbstractOptimizer} <: MOI.AbstractOptimizer
    optimizer::MOIB.LazyBridgeOptimizer{S}
    cache::MOIB.LazyBridgeOptimizer{MOIU.CachingOptimizer{S,MOIU.Model{Float64}}}
end
function Incremental(m::S) where {S <: MOI.AbstractOptimizer}
    b = MOIB.full_bridge_optimizer(m, Float64)
    cache = MOIB.full_bridge_optimizer(MOIU.CachingOptimizer(MOIU.Model{Float64}(), m), Float64)
    return Incremental{S}(b, cache)
end

function MOI.copy_to(model::Incremental{S}, src::MOI.ModelLike) where S <: MOI.AbstractOptimizer
    return MOI.Utilities.default_copy_to(model, src, copy_names)
end

_is_incremental(x) = false
_get_storage(d::Incremental{S}) where S = _is_incremental(S) ? d.optimizer : d.cache

# Set attributes
MOI.set(d::Incremental, ::MOI.ObjectiveFunction{VI}, f::VI)   = (MOI.set(_get_storage(d), MOI.ObjectiveFunction{VI}(), f); nothing)
MOI.set(d::Incremental, ::MOI.ObjectiveFunction{SAF}, f::SAF) = (MOI.set(_get_storage(d), MOI.ObjectiveFunction{SAF}(), f); nothing)
MOI.set(d::Incremental, ::MOI.ObjectiveFunction{SQF}, f::SQF) = (MOI.set(_get_storage(d), MOI.ObjectiveFunction{SQF}(), f); nothing)

function MOI.set(d::Incremental, ::MOI.ObjectiveSense, s)
    MOI.set(_get_storage(d), MOI.ObjectiveSense(), s)
    return
end
function MOI.set(d::Incremental{S}, ::MOI.Silent, s) where S <: MOI.AbstractOptimizer
    if _is_incremental(S)
        MOI.set(d.optimizer, MOI.Silent(), s)
    else
        MOI.set(d.cache.model.optimizer, MOI.Silent(), s)
    end
    return
end

function MOI.set(d::Incremental, ::MOI.VariablePrimalStart, v::VI, x)
    MOI.set(_get_storage(d), MOI.VariablePrimalStart(), v, x)
    return
end

MOI.set(d::Incremental, ::MOI.NLPBlock, ::Nothing) = nothing
function MOI.set(d::Incremental, ::MOI.NLPBlock, s)
    MOI.set(_get_storage(d), MOI.NLPBlock(), s)
    return
end

function MOI.set(d::Incremental, p::MOI.RawOptimizerAttribute, s)
    MOI.set(_get_storage(d), p, s)
    return
end

# Add variable/constraint
function MOI.add_variable(d::Incremental)
    MOI.add_variable(_get_storage(d))::VI
end

MOI.add_constraint(d::Incremental, f::VI, s::LT) = MOI.add_constraint(_get_storage(d), f, s)::CI{VI,LT}
MOI.add_constraint(d::Incremental, f::VI, s::GT) = MOI.add_constraint(_get_storage(d), f, s)::CI{VI,GT}
MOI.add_constraint(d::Incremental, f::VI, s::ET) = MOI.add_constraint(_get_storage(d), f, s)::CI{VI,ET}
MOI.add_constraint(d::Incremental, f::VI, s::IT) = MOI.add_constraint(_get_storage(d), f, s)::CI{VI,IT}
MOI.add_constraint(d::Incremental, f::VI, s::ZO) = MOI.add_constraint(_get_storage(d), f, s)::CI{VI,ZO}
MOI.add_constraint(d::Incremental, f::VI, s::MOI.Integer) = MOI.add_constraint(_get_storage(d), f, s)::CI{VI,MOI.Integer}

MOI.add_constraint(d::Incremental, f::SAF, s::LT) = MOI.add_constraint(_get_storage(d), f, s)::CI{SAF,LT}
MOI.add_constraint(d::Incremental, f::SAF, s::GT) = MOI.add_constraint(_get_storage(d), f, s)::CI{SAF,GT}
MOI.add_constraint(d::Incremental, f::SAF, s::ET) = MOI.add_constraint(_get_storage(d), f, s)::CI{SAF,ET}
MOI.add_constraint(d::Incremental, f::SAF, s::IT) = MOI.add_constraint(_get_storage(d), f, s)::CI{SAF,IT}

MOI.add_constraint(d::Incremental, f::SQF, s::LT) = MOI.add_constraint(_get_storage(d), f, s)::CI{SQF,LT}
MOI.add_constraint(d::Incremental, f::SQF, s::GT) = MOI.add_constraint(_get_storage(d), f, s)::CI{SQF,GT}
MOI.add_constraint(d::Incremental, f::SQF, s::ET) = MOI.add_constraint(_get_storage(d), f, s)::CI{SQF,ET}
MOI.add_constraint(d::Incremental, f::SQF, s::IT) = MOI.add_constraint(_get_storage(d), f, s)::CI{SQF,IT}

# Delete
function MOI.delete(d::Incremental, ci::CI{VI,T}) where T <: Union{LT,GT,ET,IT,MOI.Integer}
    MOI.delete(_get_storage(d), ci)
    return
end
function MOI.delete(d::Incremental, ci::CI{SAF,LT})
     MOI.delete(_get_storage(d), ci)
     return
end

# Set modifications
function MOI.set(d::Incremental, ::MOI.ConstraintSet, ci::CI{VI,T}, s::T) where T <: Union{LT,GT,ET,IT}
     MOI.set(_get_storage(d), MOI.ConstraintSet(), ci, s)
     return
end

# Get attributes
function MOI.get(d::Incremental{S}, ::MOI.TerminationStatus) where S
    MOI.get(d.optimizer, MOI.TerminationStatus())::MOI.TerminationStatusCode
end
function MOI.get(d::Incremental{S}, ::MOI.PrimalStatus) where S
    MOI.get(d.optimizer, MOI.PrimalStatus())::MOI.ResultStatusCode
end
function MOI.get(d::Incremental{S}, ::MOI.DualStatus) where S
    MOI.get(d.optimizer, MOI.DualStatus())::MOI.ResultStatusCode
end
function MOI.get(d::Incremental{S}, ::MOI.RawStatusString) where S
    MOI.get(d.optimizer, MOI.RawStatusString())::MOI.String
end

function MOI.get(d::Incremental{S}, ::MOI.ObjectiveBound) where S
    MOI.get(d.optimizer, MOI.ObjectiveBound())::Float64
end
function MOI.get(d::Incremental{S}, ::MOI.ObjectiveValue) where S
    MOI.get(d.optimizer, MOI.ObjectiveValue())::Float64
end
function MOI.get(d::Incremental{S}, ::MOI.DualObjectiveValue) where S
    MOI.get(d.optimizer, MOI.DualObjectiveValue())::Float64
end

function MOI.get(d::Incremental{S}, ::MOI.VariablePrimal, vi::VI) where S
    MOI.get(d.optimizer, MOI.VariablePrimal(), vi)::Float64
end

const SAF_CI_TYPES = Union{CI{SAF,LT},CI{SAF,GT},CI{SAF,ET}}
function MOI.get(d::Incremental{S}, ::MOI.ConstraintPrimal, ci::SAF_CI_TYPES) where S
    MOI.get(d.optimizer, MOI.ConstraintPrimal(), ci)::Float64
end

const SQF_CI_TYPES = Union{CI{SQF,LT},CI{SQF,GT},CI{SQF,ET}}
function MOI.get(d::Incremental{S}, ::MOI.ConstraintPrimal, ci::SQF_CI_TYPES) where S
    MOI.get(d.optimizer, MOI.ConstraintPrimal(), ci)::Float64
end

function MOI.get(d::Incremental{S}, ::MOI.ConstraintDual, ci::Union{CI{VI,LT},CI{VI,GT}}) where S
    MOI.get(d.optimizer, MOI.ConstraintDual(), ci)::Float64
end
function MOI.get(d::Incremental{S}, ::MOI.ResultCount) where S
    MOI.get(d.optimizer, MOI.ResultCount())::Int
end
function MOI.get(d::Incremental{S}, n::MOI.SolverName) where S
    MOI.get(d.optimizer, n)::String
end

# define optimize!
function MOI.optimize!(d::Incremental{S}) where S
    MOI.optimize!(_get_storage(d))
    return
end

function MOI.empty!(d::Incremental{S}) where S
    MOI.empty!(_get_storage(d))
    return
end

MOI.is_empty(d::Incremental{S}) where S = MOI.is_empty(_get_storage(d))