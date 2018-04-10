function default_pre(feas::Bool,X::Vector{Interval{T}},UBD::Float64,
                     k::Int64,pos::Int64,opt,LBDn::Float64,UBDn::Float64,x::BnBSolver,
                     y::BnBModel{Interval{T}}) where {T<:AbstractFloat}
   return feas,X
end

function default_pre(feas::Bool,X::Vector{MCInterval{T}},UBD::Float64,
                     k::Int64,pos::Int64,opt,LBDn::Float64,UBDn::Float64,x::BnBSolver,
                     y::BnBModel{MCInterval{T}}) where {T<:AbstractFloat}
   return feas,X
end

function default_post(feas::Bool,X::Vector{Interval{T}},k::Int64,
                     pos::Int64,opt,tempL,tempU,LBD,UBD) where {T<:AbstractFloat}
    return feas,X
end

function default_post(feas::Bool,X::Vector{MCInterval{T}},k::Int64,
                     pos::Int64,opt,tempL,tempU,LBD,UBD) where {T<:AbstractFloat}
    return feas,X
end
