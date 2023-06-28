
function estimator_extrema(x::MC{N,T}, y::MC{N,T}, s, dP) where {N,T}
    xcv = x.cv;   xcvg = x.cv_grad
    ycv = y.cv;   ycvg = y.cv_grad
    xccn = -x.cc; xccgn = -x.cc_grad
    yccn = -y.cc; yccgn = -y.cc_grad
    t3 = affine_expand_del(dP, xcv,  xcvg,  s)
    t4 = affine_expand_del(dP, ycv,  ycvg,  s)
    s3 = affine_expand_del(dP, xccn, xccgn, s)
    s4 = affine_expand_del(dP, yccn, yccgn, s)
    return t3, t4, s3, s4
end

function estimator_under(x::MC{N,T}, y::MC{N,T}, s, dp) where {N,T}
    xcv = x.cv;   xcvg = x.cv_grad
    ycv = y.cv;   ycvg = y.cv_grad
    u1cv = affine_expand_del(dp, xcv, xcvg, s)
    u2cv = affine_expand_del(dp, ycv, ycvg, s)
    return u1cv, u2cv, xcvg, ycvg
end

function estimator_over(x::MC{N,T}, y::MC{N,T}, s, dp) where {N,T}
    xccn = -x.cc; xccgn = -x.cc_grad
    yccn = -y.cc; yccgn = -y.cc_grad
    v1ccn = affine_expand_del(dp, xccn, xccgn, s)
    v2ccn = affine_expand_del(dp, yccn, yccgn, s)
    return v1ccn, v2ccn, xccgn, yccgn
end


include(joinpath(@__DIR__, "enumeration.jl"))
include(joinpath(@__DIR__, "affine_arithmetic.jl"))