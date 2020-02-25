function IntervalContractors.cos_main(X::IntervalBox)

    x, y = X

    x_range = Interval(0, Interval{Float64}(π).lo)
    y_range = -1..1

    x = x ∩ x_range
    y = y ∩ y_range

    isempty(IntervalBox(x, y)) && return IntervalBox(x, y)

    y = y ∩ cos(x)
    x = x ∩ acos(y)

    return IntervalBox(x, y)

end

function IntervalContractors.tan_main(X::IntervalBox)

    x, y = X

    x_range = Interval(-(Interval{Float64}(π)/2).hi, (Interval{Float64}(π)/2).hi)

    x = x ∩ x_range

    isempty(x) && return IntervalBox(x, y)

    y = y ∩ tan(x)
    x = x ∩ atan(y)

    return IntervalBox(x, y)

end

function IntervalContractors.sin_main(X::IntervalBox)

    x, y = X

    x_range = Interval(-(Interval{Float64}(π)/2).hi, (Interval{Float64}(π)/2).hi)
    y_range = -1..1

    x = x ∩ x_range
    y = y ∩ y_range

    isempty(IntervalBox(x, y)) && return IntervalBox(x, y)

    y = y ∩ sin(x)
    x = x ∩ asin(y)

    return IntervalBox(x, y)

end
