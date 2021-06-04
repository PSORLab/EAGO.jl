
# TODO: NEED VARIABLE SUBTYPES...

is_convex(::typeof(VARIABLE), x::Interval{T}) where T = true
is_concave(::typeof(VARIABLE), x::Interval{T}) where T = true

is_loglogconvex(::typeof(VARIABLE), x::Interval{T}) where T = true
is_loglogconcave(::typeof(VARIABLE), x::Interval{T}) where T = true

is_increasing(::typeof(VARIABLE), x::Interval{T}) where T = true
is_decreasing(::typeof(VARIABLE), x::Interval{T}) where T = false

is_positive(::typeof(VARIABLE), x::Interval{T}) where T = x >= zero(T)
is_negative(::typeof(VARIABLE), x::Interval{T}) where T = x <= zero(T)
is_locked(::typeof(VARIABLE), x::Interval{T}) where T = one(T) ∉ x
