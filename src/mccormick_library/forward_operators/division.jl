function div_alphaxy(es::T, nu::T, x::Interval{Float64}, y::Interval{Float64}) where T <: Real
    return (es/y.hi) + (x.lo/(y.lo*y.hi))*(y.hi-nu)
end
function div_gammay(omega::T, y::Interval{Float64}) where T <: Real
    return (y.lo*(max(0.0, omega))^2)/(y.hi - omega*(y.hi-y.lo))
end
function div_deltaxy(omega::T, x::Interval{Float64}, y::Interval{Float64}) where T <: Real
    return (1.0/(y.hi*y.lo))*(x.hi - x.lo)*(y.hi - y.lo)*div_gammay(omega, y)
end
function div_psixy(es::T, nu::T, x::Interval{Float64}, y::Interval{Float64}) where T <: Real
    return div_alphaxy(es, nu, x, y) + div_deltaxy(((es - x.lo)/(x.hi - x.lo))-((nu - y.lo)/(y.hi - y.lo)), x, y)
end
function div_omegaxy(x::Interval{Float64}, y::Interval{Float64}) where T <: Real
    return (y.hi/(y.hi-y.lo))*(1.0 - sqrt((y.lo*(x.hi-x.lo))/((-x.lo)*(y.hi-y.lo)+(y.lo)*(x.hi-x.lo))))
end
function div_lambdaxy(es::T, nu::T, x::Interval{Float64}) where T <: Real
    return (((es + sqrt(x.lo*x.hi))/(sqrt(x.lo) + sqrt(x.hi)))^2)/nu
end
function div_nuline(x::Interval{Float64}, y::Interval{Float64}, z::T) where T <: Real
    return y.lo + (y.hi - y.lo)*(z - x.lo)/(x.hi - x.lo)
end

function mid3(x::T, y::T, z::T) where T<:Real
  (((x>=y)&&(y>=z))||((z>=y)&&(y>=x))) && (return y,2)
  (((y>=x)&&(x>=z))||((z>=x)&&(x>=y))) && (return x,1)
  return z,3
end

struct DiffMCTag end

function div_diffcv(x::MC{N,T}, y::MC{N,T}) where {N, T<:RelaxTag}
    dualx_cv = Dual{DiffMCTag}(x.cv, Partials{N,Float64}(NTuple{N,Float64}(x.cv_grad)))
    dualy_cv = Dual{DiffMCTag}(y.cv, Partials{N,Float64}(NTuple{N,Float64}(x.cv_grad)))
    dualy_cc = Dual{DiffMCTag}(y.cc, Partials{N,Float64}(NTuple{N,Float64}(y.cc_grad)))
    nu_bar = div_nuline(x.Intv, y.Intv, dualx_cv)
    if (0.0 <= x.Intv.lo)
        dual_div = div_lambdaxy(dualx_cv, dualy_cc, x.Intv)
    elseif (x.Intv.lo < 0.0) && (nu_bar <= dualy_cv)
        dual_div = div_alphaxy(dualx_cv, dualy_cv, x.Intv, y.Intv)
    elseif (x.Intv.lo < 0.0) && (nu_bar > dualy_cv)
        dual_div = div_psixy(dualx_cv, mid(dualy_cv, dualy_cc, nu_bar - (y.Intv.hi - y.Intv.lo)*div_omegaxy(x.Intv, y.Intv)),
                         x.Intv, y.Intv)
    end
    val, grad = dual_div.value, SVector{N,Float64}(dual_div.partials)
    return val, grad
end

function div_MV(x::MC{N,Diff}, y::MC{N,Diff}, z::Interval{Float64}) where N
    if (0.0 < y.Intv.lo)
        cv, cv_grad = div_diffcv(x, y)
        cc, cc_grad = div_diffcv(-x, y)
        cc = -cc
        cc_grad = -cc_grad
    elseif (y.Intv.hi < 0.0)
        cv, cv_grad = div_diffcv(-x, -y)
        cc, cc_grad = div_diffcv(-x, y)
        cc = -cc
        cc_grad = -cc_grad
    else
        error("Division (x/y) is unbounded on intervals y containing 0.")
    end
    return MC{N,Diff}(cv, cc, z, cv_grad, cc_grad, x.cnst && y.cnst)
end

function div_kernel(x::MC{N,NS}, y::MC{N,NS}, z::Interval{Float64}) where N
    if (x === y)
        zMC = one(x)
    else
        q = inv(y)
        zMC = mult_kernel(x, q, z)
    end
    return zMC
end
#=
function div_kernel(x::MC{N,MV}, y::MC{N,MV}, z::Interval{Float64}) where N
    pos_orth::Bool = (x.Intv.lo >= 0) && (y.Intv.lo >= 0)
    degen1 = ((x.Intv.hi - x.Intv.lo) == 0.0)
    degen2 = ((y.Intv.hi - y.Intv.lo) == 0.0)
    error("To be implemented")
end
=#
function div_kernel(x::MC{N,Diff}, y::MC{N,Diff}, z::Interval{Float64}) where {N}
    degen1 = ((x.Intv.hi - x.Intv.lo) == 0.0)
    degen2 = ((y.Intv.hi - y.Intv.lo) == 0.0)
    if (x === y)
        zMC = one(x)
    elseif  ~(degen1||degen2)
        zMC = div_MV(x, y, z)
    else
        q = inv(y)
        zMC = mult_kernel(x, q, z)
    end
    return zMC
end

function /(x::MC{N,T}, y::MC{N,T}) where {N,T<:RelaxTag}
    @assert ~(y.Intv.lo <= 0.0 <= y.Intv.hi)
    return div_kernel(x, y, x.Intv/y.Intv)
end
