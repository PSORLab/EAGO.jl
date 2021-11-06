"""

A structure containing the 
"""
mutable struct SubSolvers{Q<:MOI.AbstractOptimizer, S<:MOI.AbstractOptimizer, T<:ExtensionType}
    relaxed_optimizer::Q
    upper_optimizer::S
    ext_typ::T
end
SubSolvers{Q,S,T}(; r::Q = GLPK.Optimizer(), u::S = Ipopt.Optimizer(), t::T = DefaultExt()) where {Q,S,T} = SubSolvers{Q,S,T}(r, u, t)
SubSolvers(; r::Q = GLPK.Optimizer(), u::S = Ipopt.Optimizer(), t::T = DefaultExt()) where {Q,S,T} = SubSolvers{Q,S,T}(r,u,t)

function _relaxed_optimizer(d::SubSolvers{Q,S,T}) where {Q <: MOI.AbstractOptimizer, 
                                                         S <: MOI.AbstractOptimizer, 
                                                         T <: ExtensionType}
    return d.relaxed_optimizer
end

function _upper_optimizer(d::SubSolvers{Q,S,T}) where {Q <: MOI.AbstractOptimizer, 
                                                       S <: MOI.AbstractOptimizer, 
                                                       T <: ExtensionType}
    return d.upper_optimizer
end

function _ext_type(d::SubSolvers{Q,S,T}) where {Q <: MOI.AbstractOptimizer, 
                                                S <: MOI.AbstractOptimizer, 
                                                T <: ExtensionType}
    return d.ext_typ
end

function isempty(d::SubSolvers{Q,S,T}) where {Q,S,T}
    MOI.is_empty(d.relaxed_optimizer) && MOI.is_empty(d.upper_optimizer)
end
function MOI.empty!(d::SubSolvers{Q,S,T}) where {Q,S,T}
    MOI.empty!(d.relaxed_optimizer)
    MOI.empty!(d.upper_optimizer)
    return 
end