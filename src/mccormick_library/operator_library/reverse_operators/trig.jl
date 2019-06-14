# Intersect Interval then promote and intersect (can do full McCormick later)

function sin_rev(a::MC{N}, b::MC{N}) where N
    aintv, bintv = sin_rev(a.Intv, b.Intv)
    a = a ∩ MC{N}(aintv)
    b = b ∩ MC{N}(bintv)
    a, b
end
function cos_rev(a::MC{N}, b::MC{N}) where N
    aintv, bintv = cos_rev(a.Intv, b.Intv)
    a = a ∩ MC{N}(aintv)
    b = b ∩ MC{N}(bintv)
    a, b
end
function tan_rev(a::MC{N}, b::MC{N}) where N
    aintv, bintv = tan_rev(a.Intv, b.Intv)
    a = a ∩ MC{N}(aintv)
    b = b ∩ MC{N}(bintv)
    a, b
end
