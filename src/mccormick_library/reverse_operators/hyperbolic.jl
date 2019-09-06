"""
    sinh_rev!

Reverse McCormick operator for `sinh`.
"""
function sinh_rev(y::MC, x::MC)
    x = x ∩ asinh(y)
    y,x
end

"""
    cosh_rev!

Reverse McCormick operator for `cosh`.
"""
function cosh_rev(y::MC, x::MC)
    y = y ∩ IntervalType(1.0, Inf)
    if ~isempty(y)
        x = x ∩ acosh(y)
    end
    y,x
end

"""
    tanh_rev!

Reverse McCormick operator for `tanh`.
"""
function tanh_rev(y::MC, x::MC)
    y = y ∩ IntervalType(-1.0, 1.0)
    if ~isempty(y)
        x = x ∩ atanh(y)
    end
    y,x
end
