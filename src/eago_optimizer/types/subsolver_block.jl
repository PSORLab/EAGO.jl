"""

A structure containing the 
"""
mutable struct SubSolvers{Q<:MOI.AbstractOptimizer, S<:MOI.AbstractOptimizer, T<:ExtensionType}
    relaxed_optimizer::Q
    upper_optimizer::S
    ext_typ::T
end
function SubSolvers{Q,S,T}(; relaxed_optimizer::Q = GLPK.Optimizer(), upper_optimizer::S = Ipopt.Optimizer(), 
                             ext_typ::T = DefaultExt()) where {Q,S,T}
    SubSolvers{Q,S,T}(relaxed_optimizer, upper_optimizer, ext_typ)
end
function SubSolvers(; relaxed_optimizer::Q = GLPK.Optimizer(), upper_optimizer::S = Ipopt.Optimizer(), 
                      ext_typ::T = DefaultExt()) where {Q,S,T}
    return SubSolvers{Q,S,T}(relaxed_optimizer = relaxed_optimizer, upper_optimizer = upper_optimizer, ext_typ = ext_typ)
end

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