# Defines type
"""
    HybridMC

Defines the hybridMC type used for constructing nonstandard McCormick relaxations
holds the SMC type.
"""
struct HybridMC{N,V,T} <: Real
    SMC::SMCg{N,V,T}
end

############################## constructors ####################################

#################### conversion/promotion functions ############################

function convert(::Type{HybridMC{N,V,T}},x::S) where {S<:Integer,N,V<:AbstractInterval,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          HybridMC{N,V,T}(SMCg{N,V,T}(convert(T,x),convert(T,x),seed,seed,V(convert(V,x)),false))
end
function convert(::Type{HybridMC{N,V,T}},x::S) where {S<:AbstractFloat,N,V<:AbstractInterval,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          HybridMC{N,V,T}(SMCg{N,V,T}(convert(T,x),convert(T,x),seed,seed,V(convert(V,x)),false))
end
function convert(::Type{HybridMC{N,V,T}},x::S) where {S<:Interval,N,V<:AbstractInterval,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          HybridMC{N,V,T}(SMCg{N,V,T}(convert(T,x.hi),convert(T,x.lo),seed,seed,convert(V,x),false))
end

promote_rule(::Type{HybridMC{N,V,T}}, ::Type{S}) where {S<:Integer,N,V<:AbstractInterval,T<:AbstractFloat} = HybridMC{N,V,T}
promote_rule(::Type{HybridMC{N,V,T}}, ::Type{S}) where {S<:AbstractFloat,N,V<:AbstractInterval,T<:AbstractFloat} = HybridMC{N,V,T}
promote_rule(::Type{HybridMC{N,V,T}}, ::Type{S}) where {S<:Interval,N,V<:AbstractInterval,T<:AbstractFloat} = HybridMC{N,V,T}
promote_rule(::Type{HybridMC{N,V,T}}, ::Type{S}) where {S<:Real,N,V<:AbstractInterval,T<:AbstractFloat} = HybridMC{N,V,T}

###########################  Access functions  #################################
"""
    cc
Access concave relaxation of hybrid McCormick object.
"""
cc(x::HybridMC) = x.SMC.cc

"""
    cv

Access convex relaxation of hybrid McCormick object.
"""
cv(x::HybridMC) = x.SMC.cv

"""
    cc_grad
Access concave (sub)gradient of hybrid McCormick object.
"""
cc_grad(x::HybridMC) = x.SMC.cc_grad

"""
    cv_grad
Access convex (sub)gradient of hybrid McCormick object.
"""
cv_grad(x::HybridMC) = x.SMC.cv_grad

"""
    Intv
Access the interval bounds of the hybrid McCormick object.
"""
Intv(x::HybridMC) = x.SMC.Intv
