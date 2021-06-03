"""
$(TYPEDEF)

A structure used to store information related to the bounds assigned to each
variable.

$(TYPEDFIELDS)
"""
Base.@kwdef struct VariableInfo{T<:AbstractFloat}
    "Is the variable integer valued?"
    is_integer::Bool                                            = false
    "Boolean indicating whether finite lower bound exists."
    has_lower_bound::Bool                                       = false
    "Boolean indicating whether finite upper bound exists."
    has_upper_bound::Bool                                       = false
    "Boolean indicating variable is fixed to a finite value."
    is_fixed::Bool                                              = false
    "Indicates that constraints have been set"
    has_constraints::Bool                                       = false
    "Lower bounds. May be -Inf."
    lower_bound::T                                              = typemin(T)
    "Upper bounds. May be Inf."
    upper_bound::T                                              = typemax(T)
end
is_integer(x::VariableInfo) = x.is_integer
has_lower_bound(x::VariableInfo) = x.has_lower_bound
has_upper_bound(x::VariableInfo) = x.has_upper_bound
lower_bound(x::VariableInfo{T}) where {T <: AbstractFloat} = x.lower_bound
upper_bound(x::VariableInfo{T}) where {T <: AbstractFloat} = x.upper_bound

is_fixed(x::VariableInfo) = x.is_fixed
function is_less_than(x::VariableInfo)
    flag = x.has_upper_bound
    flag &= !x.has_lower_bound
    flag &= !x.is_integer
    return flag
end
function is_greater_than(x::VariableInfo)
    flag = x.has_lower_bound
    flag &= !x.has_upper_bound
    flag &= !x.is_integer
    return flag
end
function is_int_interval(x::VariableInfo)
    flag = x.has_lower_bound
    flag &= x.has_upper_bound
    flag &= x.is_integer
    return flag
end
function is_real_interval(x::VariableInfo)
    flag = x.has_lower_bound
    flag &= x.has_upper_bound
    flag &= !x.is_integer
    return flag
end
function is_zero_one(x::VariableInfo{T}) where {T <: AbstractFloat}
    flag = iszero(x.lower_bound)
    flag &= isone(x.upper_bound)
    flag &= x.is_integer
    return flag
end

mid(x::VariableInfo{T}) where {T <: AbstractFloat} = 0.5*(upper_bound(x) - lower_bound(x))
diam(x::VariableInfo{T}) where {T <: AbstractFloat} = 0.5*(upper_bound(x) - lower_bound(x))
empty_variable_info(::Type{T}) where T = VariableInfo{T}(lower_bound = Inf,
                                                         upper_bound = -Inf)

Base.isempty(v::VariableInfo{T}) where {T <: AbstractFloat} = lower_bound(v) > upper_bound(v)
function check_isempty(l, u, b)
    flag = l < u
    if b
        flag &= (l <= 0.0) | (u >= 1.0)
    end
    return !flag
end

function VariableInfo(::Type{T}, ::ZO) where {T <: AbstractFloat}
    return VariableInfo{T}(is_integer = true,
                           has_lower_bound = true,
                           has_upper_bound = true,
                           has_constraints = true,
                           lower_bound = zero(T),
                           upper_bound = one(T))
end

function VariableInfo(it::MOI.Interval{T}) where {T <: AbstractFloat}
    VariableInfo{T}(has_lower_bound = !isinf(it.lower),
                    has_upper_bound = !isinf(it.upper),
                    has_constraints = !isinf(it.lower) | !isinf(it.upper),
                    is_fixed = it.lower == it.upper,
                    lower_bound = it.lower,
                    upper_bound = it.upper)
end

function VariableInfo(v::VariableInfo{T}, ::ZO) where {T <: AbstractFloat}
    isempty(v) && (return v)
    l = max(zero(T), lower_bound(v))
    u = min(one(T), upper_bound(v))
    check_isempty(l, u, is_integer(v)) && return empty_variable_info(T)
    return VariableInfo{T}(is_integer = true,
                        has_lower_bound = true,
                        has_upper_bound = true,
                        has_constraints = true,
                        is_fixed = l == u,
                        lower_bound = l,
                        upper_bound = u)
end

function VariableInfo(v::VariableInfo{T}, it::MOI.Interval{T}) where {T <: AbstractFloat}
    isempty(v) && return v
    l = max(it.lower, lower_bound(v))
    u = min(it.upper, upper_bound(v))
    check_isempty(l, u, is_integer(v)) && return empty_variable_info(T)
    return VariableInfo(is_integer = is_integer(v),
                        has_lower_bound = !isinf(l),
                        has_upper_bound = !isinf(u),
                        has_constraints = !isinf(l) | !isinf(u),
                        is_fixed = l == u,
                        lower_bound = l,
                        upper_bound = u)
end

function VariableInfo(v::VariableInfo{T}, gt::MOI.GreaterThan{T}) where {T <: AbstractFloat}
    isempty(v) && return v
    l = max(gt.lower, lower_bound(v))
    u = upper_bound(v)
    check_isempty(l, u, is_integer(v)) && return empty_variable_info(T)
    return VariableInfo(is_integer = is_integer(v),
                        has_lower_bound = !isinf(l),
                        has_upper_bound = !isinf(u),
                        has_constraints = !isinf(l) | !isinf(u),
                        is_fixed = l == u,
                        lower_bound = l,
                        upper_bound = u)
end
function VariableInfo(v::VariableInfo{T}, lt::MOI.LessThan{T}) where {T <: AbstractFloat}
    isempty(v) && return v
    l = lower_bound(v)
    u = min(lt.upper, upper_bound(v))
    check_isempty(l, u, is_integer(v)) && return empty_variable_info(T)
    return VariableInfo(is_integer = is_integer(v),
                        has_lower_bound = !isinf(l),
                        has_upper_bound = !isinf(u),
                        has_constraints = !isinf(l) | !isinf(u),
                        is_fixed = l == u,
                        lower_bound = l,
                        upper_bound = u)
end

ZO(v::VariableInfo)  = MOI.ZeroOne()
ET(v::VariableInfo{T}) where {T <: AbstractFloat}  = MOI.EqualTo{T}(v.lower_bound)
IT(v::VariableInfo{T}) where {T <: AbstractFloat}  = MOI.Interval{T}(v.lower_bound)
GT(v::VariableInfo{T}) where {T <: AbstractFloat}  = MOI.GreaterThan{T}(v.lower_bound)
LT(v::VariableInfo{T}) where {T <: AbstractFloat}  = MOI.LessThan{T}(v.upper_bound)
INT(v::VariableInfo{T}) where {T <: AbstractFloat} = MOI.Semiinteger{T}(v.lower_bound, v.upper_bound)
