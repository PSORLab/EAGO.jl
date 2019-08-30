@inline div_alphaxy(es::Float64, nu::Float64, x::Interval{Float64}, y::Interval{Float64}) = (es/y.hi) + (x.lo/(y.lo*y.hi))*(y.hi-nu)
@inline div_gammay(omega::Float64, y::Interval{Float64}) = (y.lo*(max(0.0, omega))^2)/(y.hi - omega*(y.hi-y.lo))
@inline div_deltaxy(omega::Float64, x::Interval{Float64}, y::Interval{Float64}) = (1.0/(y.hi*y.lo))*(x.hi - x.lo)*(y.hi - y.lo)*div_gammay(omega, y)
@inline div_psixy(es::Float64, nu::Float64, x::Interval{Float64}, y::Interval{Float64}) = div_alphaxy(es, nu, x, y) + div_deltaxy(((es - x.lo)/(x.hi - x.lo))-((nu - y.lo)/(y.hi - y.lo)), x, y)
@inline div_omegaxy(x::Interval{Float64}, y::Interval{Float64}) = (y.hi/(y.hi-y.lo))*(1.0 - sqrt((y.lo*(x.hi-x.lo))/((-x.lo)*(y.hi-y.lo)+(y.lo)*(x.hi-x.lo))))
@inline div_lambdaxy(es::Float64, nu::Float64, x::Interval{Float64}) = (((es + sqrt(x.lo*x.hi))/(sqrt(x.lo) + sqrt(x.hi)))^2)/nu
@inline div_nuline(x::Interval{Float64}, y::Interval{Float64}, z::Float64) = y.lo + (y.hi - y.lo)*(z - x.lo)/(x.hi - x.lo)
@inline function div_diffcv(x::MC, y::MC)
    nu_bar = div_nuline(x.Intv, y.Intv, x.cv)
    (0.0 <= x.Intv.lo) && (return div_lambdaxy(x.cv, y.cc, x.Intv))
    ((x.Intv.lo < 0.0) && (nu_bar <= y.cv)) && (return div_alphaxy(x.cv, y.cv, x.Intv, y.Intv))
    if (x.Intv.lo < 0.0) && (nu_bar > y.cv)
        return div_psixy(x.cv, mid(y.cv, y.cc, nu_bar - (y.Intv.hi - y.Intv.lo)*div_omegaxy(x.Intv, y.Intv)), x.Intv, y.Intv)
    end
end

@inline function div_MV(x::MC, y::MC)
    if (0.0 < y.Intv.lo)
        cv = div_diffcv(x, y)
        cc = -div_diffcv(-x, y)
        #cv_grad =
        #cc_grad =
    elseif (y.Intv.hi < 0.0)
        cv = div_diffcv(-x, -y)
        cc = -div_diffcv(-x, y)
        #cv_grad =
        #cc_grad =
    else
        error("Division (x/y) is unbounded on intervals y containing 0.")
    end
    return cv, cc, cv_grad, cc_grad
end

@inline function div_kernel(x::MC{N}, y::MC{N}, z::Interval{Float64}) where N
    pos_orth::Bool = (x.Intv.lo >= 0) && (y.Intv.lo >= 0)
    degen1 = ((x.Intv.hi - x.Intv.lo) == 0.0)
    degen2 = ((y.Intv.hi - y.Intv.lo) == 0.0)
    if (x === y)
        cv, cc, cv_grad, cc_grad = one_kernel(x)
    elseif  (MC_param.mu >= 1 && ~(degen1||degen2))
        cv, cc, cv_grad, cc_grad = div_MV(x, y, z)
    elseif (MC_param.multivar_refine) && (MC_param.mu < 1) && pos_orth
        # TODO ADD Multivariant Nonsmooth Division
    else
        q = inv(x)
        cv, cc, cv_grad, cc_grad = mult_kernel(x, q, z)
        cv, cc, cv_grad, cc_grad = cut(z.lo, z.hi, cv, cc, cv_grad, cc_grad)
    end
    return MC{N}(cv, cc, z, cv_grad, cc_grad, (x.cnst && y.cnst))
end

@inline function /(x::MC, y::MC)
    @assert ~(y.lo <= 0.0 <= y.hi)
    return div_kernel(x, y, x.Intv/y.Intv)
end
