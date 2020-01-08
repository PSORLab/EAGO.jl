# Intersect Interval then promote and intersect (can do full McCormick later)

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `sin`.
"""
function sin_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    aintv, bintv = sin_rev(a.Intv, b.Intv)
    a = a ∩ MC{N,T}(aintv)
    b = b ∩ MC{N,T}(bintv)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `cos`.
"""
function cos_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    aintv, bintv = cos_rev(a.Intv, b.Intv)
    a = a ∩ MC{N,T}(aintv)
    b = b ∩ MC{N,T}(bintv)
    a, b
end

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `tan`.
"""
function tan_rev(a::MC{N,T}, b::MC{N,T}) where {N, T<:RelaxTag}
    aintv, bintv = tan_rev(a.Intv, b.Intv)
    a = a ∩ MC{N,T}(aintv)
    b = b ∩ MC{N,T}(bintv)
    a, b
end
