"""
$(FUNCTIONNAME)

Reverse McCormick operator for `asin`.
"""
function asin_rev(y::MC, x::MC)  # y = asin(x)
    h = lo(half_pi(Float64))
    y = y ∩ Interval{Float64}(-h, h)
    if ~isempty(y)
        x = sin(y)
    end
    y,x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `acos`.
"""
function acos_rev(y::MC, x::MC)
    y = y ∩ Interval{Float64}(0.0, hi(half_pi(Float64)))
    if ~isempty(y)
        x = x ∩ cos(y)
    end
    y,x
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `atan`.
"""
function atan_rev(y::MC, x::MC)
    y = y ∩ Interval{Float64}(-lo(half_pi(Float64)), hi(half_pi(Float64)))
    if ~isempty(y)
        x = x ∩ tan(y)
    end
    y,x
end
