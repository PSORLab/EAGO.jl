"""
    $(TYPEDEF)

A structure containing the relaxed and upper optimizers to be used, as well as any user-defined
extension.

$(TYPEDFIELDS)
"""
mutable struct SubSolvers{Q<:MOI.AbstractOptimizer, S<:MOI.AbstractOptimizer, T<:ExtensionType}
    "Optimizer used to solve relaxed subproblems. Set using `r = [...]` (<: `MOI.AbstractOptimizer`) 
        (default = `Cbc.Optimizer()`)"
    relaxed_optimizer::Q
    "Optimizer used to solve upper bounding problems. Set using `u = [...]` (<: `MOI.AbstractOptimizer`) 
        (default = `Ipopt.Optimizer()`)"
    upper_optimizer::S
    "User-defined extension to use. Set using `t = [...]`(<: `EAGO.ExtensionType`)"
    ext::T
end
SubSolvers{Q,S,T}(; r::Q = Cbc.Optimizer(), u::S = IpoptMathOptInterfaceExt.Optimizer(), t::T = DefaultExt()) where {Q,S,T} = SubSolvers{Q,S,T}(r, u, t)
SubSolvers(; r::Q = Cbc.Optimizer(), u::S = IpoptMathOptInterfaceExt.Optimizer(), t::T = DefaultExt()) where {Q,S,T} = SubSolvers{Q,S,T}(r,u,t)

"""
    _relaxed_optimizer(::GlobalOptimizer)
    _relaxed_optimizer(::SubSolvers)

Return the relaxed optimizer (<: `MOI.AbstractOptimizer`) from the `SubSolvers` field or object.
"""
function _relaxed_optimizer(d::SubSolvers{Q,S,T}) where {Q <: MOI.AbstractOptimizer, 
                                                         S <: MOI.AbstractOptimizer, 
                                                         T <: ExtensionType}
    return d.relaxed_optimizer
end

"""
    _upper_optimizer(::GlobalOptimizer)
    _upper_optimizer(::SubSolvers)

Return the upper optimizer (<: `MOI.AbstractOptimizer`) from the `SubSolvers` field or object.
"""
function _upper_optimizer(d::SubSolvers{Q,S,T}) where {Q <: MOI.AbstractOptimizer, 
                                                       S <: MOI.AbstractOptimizer, 
                                                       T <: ExtensionType}
    return d.upper_optimizer
end

"""
    _ext(::Optimizer)
    _ext(::GlobalOptimizer)
    _ext(::SubSolvers)

Return the extension (<: `ExtensionType`) from the `SubSolvers` (sub)field or object.
"""
function _ext(d::SubSolvers{Q,S,T}) where {Q <: MOI.AbstractOptimizer, 
                                                S <: MOI.AbstractOptimizer, 
                                                T <: ExtensionType}
    return d.ext
end

function isempty(d::SubSolvers{Q,S,T}) where {Q,S,T}
    MOI.is_empty(d.relaxed_optimizer) && MOI.is_empty(d.upper_optimizer)
end
function MOI.empty!(d::SubSolvers{Q,S,T}) where {Q,S,T}
    MOI.empty!(d.relaxed_optimizer)
    MOI.empty!(d.upper_optimizer)
    return 
end