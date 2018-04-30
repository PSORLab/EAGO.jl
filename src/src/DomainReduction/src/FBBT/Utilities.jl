function mul_revDR(a::T,b::T,c::T) where {T<:AbstractInterval}
    if (in(zero(T),c) && ~in(zero(T),b))
        c = c ∩ (a / b)
    elseif (~in(zero(T),c) && in(zero(T),b))
        b = b ∩ (a / c)
    elseif (~in(zero(T),c) && ~in(zero(T),b))
        b = b ∩ (a / c)
        c = c ∩ (a / b)
    end
    return a, b, c
end

function exp_revDR(y::Interval{T}, x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(zero(T),Inf)
    x_new = x ∩ log(y_new)
    return y_new, x_new
end
function exp_revDR(y::MCInterval{T}, x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(zero(T),Inf)
    x_new = x ∩ log(y_new)
    return y_new, x_new
end

function acos_rev(y::Interval{T}, x::Interval{T}) where {T<:AbstractFloat}
        y_new = y ∩ Interval{T}(zero(T),2.0*pi)
        x_new = x ∩ cos(y_new)

        return y_new, x_new
end
function acos_rev(y::MCInterval{T}, x::MCInterval{T}) where {T<:AbstractFloat}
        y_new = y ∩ MCInterval{T}(zero(T),2.0*pi)
        x_new = x ∩ cos(y_new)

        return y_new, x_new
end

function atan_rev(y::Interval{T}, x::Interval{T}) where {T<:AbstractFloat}
        y_new = y ∩ Interval{T}(-pi/2.0,pi/2.0)
        x_new = x ∩ tan(y_new)

        return y_new, x_new
end
function atan_rev(y::MCInterval{T}, x::MCInterval{T}) where {T<:AbstractFloat}
        y_new = y ∩ MCInterval{T}(-pi/2.0,pi/2.0)
        x_new = x ∩ tan(y_new)

        return y_new, x_new
end

function sinh_rev(y::T,x::T) where {T<:AbstractInterval}
    x_new = x ∩ asinh(y)

    return y, x_new
end

function cosh_rev(y::Interval{T},x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(one(T),∞)
    x_new = x ∩ acosh(y_new)

    return y_new, x_new
end

function cosh_rev(y::MCInterval{T},x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(one(T),∞)
    x_new = x ∩ acosh(y_new)

    return y_new, x_new
end

function tanh_rev(y::Interval{T},x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(-one(T),one(T))
    x_new = x ∩ atanh(y_new)

    return y_new, x_new
end

function tanh_rev(y::MCInterval{T},x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(-one(T),one(T))
    x_new = x ∩ atanh(y_new)

    return y_new, x_new
end

function asinh_rev(y::T,x::T) where {T<:AbstractInterval}
    x_new = x ∩ sinh(y)

    return y, x_new
end

function acosh_rev(y::Interval{T},x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(zero(T),∞)
    x_new = x ∩ cosh(y_new)

    return y_new, x_new
end

function acosh_rev(y::MCInterval{T},x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(zero(T),∞)
    x_new = x ∩ cosh(y_new)

    return y_new, x_new
end

function atanh_rev(y::T,x::T) where {T<:AbstractInterval}
    x_new = x ∩ tanh(y)

    return y, x_new
end

function step_rev(y::Interval{T},x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(zero(T),one(T))

    return y_new, x
end

function step_rev(y::MCInterval{T},x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(zero(T),one(T))

    return y_new, x
end

function sign_rev(y::Interval{T},x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(-one(T),one(T))

    return y_new, x
end

function sign_rev(y::MCInterval{T},x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(-one(T),one(T))

    return y_new, x
end

function exp2_rev(y::Interval{T},x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(zero(T),∞)
    x_new = x ∩ log2(y_new)

    return y_new, x_new
end

function exp2_rev(y::MCInterval{T},x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(zero(T),∞)
    x_new = x ∩ log2(y_new)

    return y_new, x_new
end

function exp10_rev(y::Interval{T},x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(zero(T),∞)
    x_new = x ∩ log10(y_new)

    return y_new, x_new
end

function exp10_rev(y::MCInterval{T},x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(zero(T),∞)
    x_new = x ∩ log10(y_new)

    return y_new, x_new
end

function log2_rev(y::Interval{T},x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(zero(T),∞)
    x_new = x ∩ exp2(y_new)

    return y_new, x_new
end

function log2_rev(y::MCInterval{T},x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(zero(T),∞)
    x_new = x ∩ exp2(y_new)

    return y_new, x_new
end

function log10_rev(y::Interval{T},x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(zero(T),∞)
    x_new = x ∩ exp10(y_new)

    return y_new, x_new
end

function log10_rev(y::MCInterval{T},x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(zero(T),∞)
    x_new = x ∩ exp10(y_new)

    return y_new, x_new
end

function one_rev(y::Interval{T},x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(one(T),one(T))

    return y_new, x
end

function one_rev(y::MCInterval{T},x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(one(T),one(T))

    return y_new, x
end

function zero_rev(y::Interval{T},x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(zero(T),zero(T))

    return y_new, x
end

function zero_rev(y::MCInterval{T},x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(zero(T),zero(T))

    return y_new, x
end
#=
function min_rev(z::Interval{T},x::Interval{T},y::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(zero(T),∞)
    x_new = x ∩ log(y_new)

    return y_new, x_new
end

function min_rev(z::MCInterval{T},x::MCInterval{T},y::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(zero(T),∞)
    x_new = x ∩ log(y_new)

    return y_new, x_new
end

function max_rev(z::Interval{T},x::Interval{T},y::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(zero(T),∞)
    x_new = x ∩ log(y_new)

    return y_new, x_new
end

function max_rev(z::MCInterval{T},x::MCInterval{T},y::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(zero(T),∞)
    x_new = x ∩ log(y_new)

    return y_new, x_new
end
=#
mul_revDR(a,b,c) = mul_revDR(promote(a,b,c)...)
exp_revDR(a,b) = exp_revDR(promote(a,b)...)
sinh_rev(a,b) = sinh_rev(promote(a,b)...)
cosh_rev(a,b) = cosh_rev(promote(a,b)...)
tanh_rev(a,b) = tanh_rev(promote(a,b)...)
asinh_rev(a,b) = asinh_rev(promote(a,b)...)
acosh_rev(a,b) = acosh_rev(promote(a,b)...)
atanh_rev(a,b) = atanh_rev(promote(a,b)...)
step_rev(a,b) = step_rev(promote(a,b)...)
sign_rev(a,b) = sign_rev(promote(a,b)...)
exp2_rev(a,b) = exp2_rev(promote(a,b)...)
exp10_rev(a,b) = exp10_rev(promote(a,b)...)
log2_rev(a,b) = log2_rev(promote(a,b)...)
log10_rev(a,b) = log10_rev(promote(a,b)...)
one_rev(a,b) = one_rev(promote(a,b)...)
zero_rev(a,b) = zero_rev(promote(a,b)...)
