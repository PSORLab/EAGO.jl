"""
$(FUNCTIONNAME)

Reverse McCormick operator for `sinh`.
"""
function sinh_rev(y::MC, x::MC)
    x = x ∩ asinh(y)
    y,x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `cosh`.
"""
function cosh_rev(y::MC, x::MC)
    y = y ∩ Interval{Float64}(1.0, Inf)
    if ~isempty(y)
        x = x ∩ acosh(y)
    end
    y,x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `tanh`.
"""
function tanh_rev(y::MC, x::MC)
    y = y ∩ Interval{Float64}(-1.0, 1.0)
    if ~isempty(y)
        x = x ∩ atanh(y)
    end
    y,x
end
