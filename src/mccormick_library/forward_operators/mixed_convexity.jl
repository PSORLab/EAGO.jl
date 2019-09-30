# convex relaxation (envelope) of cos function
@inline function cv_cos(x::Float64, xL::Float64, xU::Float64, tp1::Float64, tp2::Float64)
  r = 0.0
  kL = Base.ceil((-0.5-xL/(2.0*pi)))
  if (x<=(pi-2.0*pi*kL))
    xL1 = xL+2.0*pi*kL
    if (xL1 >= pi/2.0)
      return cos(x), -sin(x), tp1, tp2
    end
    xU1 = min(xU+2.0*pi*kL,pi)
    if ((xL1 >= (-pi/2)) && (xU1 <= (pi/2)))
      if (abs(xL-xU) < MC_ENV_TOL)
        r = 0.0
      else
        r = (cos(xU) - cos(xL))/(xU - xL)
      end
      return cos(xL)+r*(x - xL), r, tp1, tp2
    end
    val, dval, tp1 = cv_cosin(x+(2.0*pi)*kL, xL1, xU1, tp1)
    return val, dval, tp1, tp2
  end
  kU = Base.floor((0.5-xU/(2.0*pi)))
  if (x >= (-pi-2.0*pi*kU))
    xU2 = xU+2.0*pi*kU
    if (xU2 <=-pi/2.0)
      return cos(x), -sin(x), tp1, tp2
    end
    val, dval, tp2 = cv_cosin(x+2.0*pi*kU, max(xL+2.0*pi*kU,-pi), xU2, tp2)
    return val, dval, tp1, tp2
  end
  return -1.0, 0.0, tp1, tp2
end
# function for computing convex relaxation over nonconvex and nonconcave regions
@inline function cv_cosin(x::Float64, xL::Float64, xU::Float64, xj::Float64)
  if (abs(xL) <= abs(xU))
    left = false
    x0 = xU
    xm = xL
  else
    left = true
    x0 = xL
    xm = xU
  end
  if xj === Inf
    xj, flag = newton(x0, xL, xU, cv_cosenv, dcv_cosenv, xm, 0.0)
    if flag
        xj = golden_section(xL, xU, cv_cosenv, xm, 0.0)
      end
  end
  if ((left && x <= xj)||((~left) && x >= xj))
    return cos(x), -sin(x), xj
  else
    if abs(xm - xj) < MC_ENV_TOL
      r = 0.0
    else
      r = (cos(xm) - cos(xj))/(xm - xj)
    end
    return cos(xm) + r*(x - xm), r, xj
  end
end
@inline cv_cosenv(x::Float64, y::Float64, z::Float64) = (x - y)*sin(x) + cos(x) - cos(y)
@inline dcv_cosenv(x::Float64, y::Float64, z::Float64) = (x - y)*cos(x)
# concave relaxation (envelope) of cos function
@inline function cc_cos(x::Float64, xL::Float64, xU::Float64, tp1::Float64, tp2::Float64)
  temp = cv_cos(x - pi, xL - pi, xU - pi, tp1, tp2)
  return -temp[1],-temp[2], temp[3], temp[4]
end
@inline function cos_arg(xL::Float64, xU::Float64)
  kL = Base.ceil(-0.5-xL/(2.0*pi))
  xL1 = xL + 2.0*pi*kL
  xU1 = xU + 2.0*pi*kL
  if ~( (xL1 >= -pi) && (xL1 <= pi) )
    error("Cosine Argument Calculation: xL out of bounds.")
  end
  if (xL1 <= 0.0)
    if (xU1 <= 0.0)
      arg1 = xL
      arg2 = xU
    elseif (xU1 >= pi)
      arg1 = pi-2.0*pi*kL
      arg2 = -2.0*pi*kL
    else
      arg1 = (cos(xL1) <= cos(xU1)) ? xL : xU
      arg2 = -2.0*pi*kL
    end
  end
  if (xU1 <= pi)
    arg1 = xU
    arg2 = xL
  elseif (xU1 >= (2.0*pi))
    arg1 = pi-2.0*pi*kL
    arg2 = 2.0*pi-2.0*pi*kL
  else
    arg1 = pi-2.0*pi*kL
    arg2 = (cos(xL1) >= cos(xU1)) ? xL : xU
  end
  return arg1,arg2
end
@inline function cos_kernel(x::MC{N, Diff}, y::Interval{Float64}, cv_tp1::Float64,
                              cv_tp2::Float64, cc_tp1::Float64, cc_tp2::Float64) where N
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = y.lo
  xUc = y.hi
  eps_max, eps_min = cos_arg(x.Intv.lo, x.Intv.hi)
  midcc,cc_id = mid3(x.cc, x.cv, eps_max)
  midcv,cv_id = mid3(x.cc, x.cv, eps_min)
  cc, dcc, cc_tp1, cc_tp2 = cc_cos(midcc, x.Intv.lo, x.Intv.hi, cc_tp1, cc_tp2)
  cv, dcv, cv_tp1, cv_tp2 = cv_cos(midcv, x.Intv.lo, x.Intv.hi, cv_tp1, cv_tp2)
  gcc1, gdcc1, cc_tp1, cc_tp2 = cc_cos(x.cv, x.Intv.lo, x.Intv.hi, cc_tp1, cc_tp2)
  gcv1, gdcv1, cv_tp1, cv_tp2 = cv_cos(x.cv, x.Intv.lo, x.Intv.hi, cv_tp1, cv_tp2)
  gcc2, gdcc2, cc_tp1, cc_tp2 = cc_cos(x.cc, x.Intv.lo, x.Intv.hi, cc_tp1, cc_tp2)
  gcv2, gdcv2, cv_tp1, cv_tp2 = cv_cos(x.cc, x.Intv.lo, x.Intv.hi, cv_tp1, cv_tp2)
  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  return MC{N, Diff}(cv, cc, y, cv_grad, cc_grad, x.cnst), cv_tp1, cv_tp2, cc_tp1, cc_tp2
end
@inline function cos_kernel(x::MC{N, NS}, y::Interval{Float64}, cv_tp1::Float64,
                              cv_tp2::Float64, cc_tp1::Float64, cc_tp2::Float64) where N
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = y.lo
  xUc = y.hi
  eps_max, eps_min = cos_arg(x.Intv.lo, x.Intv.hi)
  midcc,cc_id = mid3(x.cc, x.cv, eps_max)
  midcv,cv_id = mid3(x.cc, x.cv, eps_min)
  cc, dcc, cc_tp1, cc_tp2 = cc_cos(midcc, x.Intv.lo, x.Intv.hi, cc_tp1, cc_tp2)
  cv, dcv, cv_tp1, cv_tp2 = cv_cos(midcv, x.Intv.lo, x.Intv.hi, cv_tp1, cv_tp2)
  cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
  cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  return MC{N,NS}(cv, cc, y, cv_grad, cc_grad, x.cnst), cv_tp1, cv_tp2, cc_tp1, cc_tp2
end
@inline function cos(x::MC)
  y, tp1, tp2, tp3, tp4 = cos_kernel(x, cos(x.Intv), Inf, Inf, Inf, Inf)
  return y
end
@inline function sin_kernel(x::MC, y::Interval{Float64}, cv_tp1::Float64,
                            cv_tp2::Float64, cc_tp1::Float64, cc_tp2::Float64)
    cos_kernel(x-pi/2.0, y, cv_tp1, cv_tp2, cc_tp1, cc_tp2)
end
@inline function sin(x::MC)
  y, tp1, tp2, tp3, tp4 = sin_kernel(x, sin(x.Intv), Inf, Inf, Inf, Inf)
  return y
end

@inline sinh_deriv(x::Float64) = cosh(x)
@inline sinh_env(x::Float64, y::Float64, z::Float64) = (x-y)*cosh(x)-(sinh(x)-sinh(y))
@inline function cv_sinh(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return sinh(x), sinh_deriv(x), p)
  (xU <= 0.0) && (return dline_seg(sinh, sinh_deriv, x, xL, xU)..., p)
  if p === Inf
      p, flag = secant(xL, 0.0, xL, 0.0, sinh_env, xU, 0.0)
      if flag
          p = golden_section(xL, 0.0, sinh_env, xU, 0.0)
      end
  end
  (x <= p) && (return dline_seg(sinh, sinh_deriv, x, xL, p)..., p)
  return sinh(x), sinh_deriv(x), p
end
@inline function cc_sinh(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return dline_seg(sinh, sinh_deriv, x, xL, xU)..., p)
  (xU <= 0.0) && (return sinh(x), sinh_deriv(x), p)
  if p === Inf
      p, flag = secant(0.0, xU, 0.0, xU, sinh_env, xL, 0.0)
      if flag
          p = golden_section(0.0, xU, sinh_env, xL, 0.0)
      end
  end
  (x <= p) && (return sinh(x), sinh_deriv(x), p)
  return dline_seg(sinh, sinh_deriv, x, p, xU)..., p
end

@inline atanh_deriv(x::Float64) = 1.0/(1.0 - x^2)
@inline atanh_env(x::Float64, y::Float64, z::Float64) = (atanh(x) - atanh(y))*(1.0 - x^2) - x + y
@inline function cv_atanh(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return atanh(x), atanh_deriv(x), p)
  (xU <= 0.0) && (return dline_seg(atanh, atanh_deriv, x, xL, xU)..., p)
  if p === Inf
      p, flag = secant(xL, 0.0, xL, 0.0, atanh_env, xU, 0.0)
      if flag
          p = golden_section(xL, 0.0, atanh_env, xU, 0.0)
      end
  end
  (x <= p) && (return dline_seg(atanh, atanh_deriv, x, xL, p)..., p)
  return atanh(x), atanh_deriv(x), p
end
@inline function cc_atanh(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return dline_seg(atanh, atanh_deriv, x, xL, xU)..., p)
  (xU <= 0.0) && (return atanh(x), atanh_deriv(x), p)
  if p === Inf
      p, flag = secant(0.0, xU, 0.0, xU, atanh_env, xL, 0.0)
      if flag
          p = golden_section(0.0, xU, atanh_env, xL, 0.0)
      end
  end
  (x <= p) && (return atanh(x), atanh_deriv(x), p)
  return dline_seg(atanh, atanh_deriv, x, p, xU)..., p
end

@inline tanh_deriv(x::Float64, y::Float64, z::Float64) = sech(x)^2
@inline tanh_env(x::Float64, y::Float64, z::Float64) = (x-y)-(tanh(x)-tanh(y))/(1.0-tanh(x)^2)
@inline tanh_envd(x::Float64, y::Float64, z::Float64)= 2*tanh(x)/(1.0-tanh(x)^2)*(tanh(x)-tanh(y))
@inline function cv_tanh(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return dline_seg(tanh, tanh_deriv, x, xL, xU)..., p)
  (xU <= 0.0) && (return tanh(x), sech(x)^2, p)
  if p === Inf
    p, flag = secant(xL, 0.0, xL, 0.0, tanh_env, xU, 0.0)
    if flag
        p = golden_section(xL, 0.0, tanh_env, xU, 0.0)
    end
   end
   (x <= p) && (return tanh(x), sech(x)^2, p)
   return dline_seg(tanh, tanh_deriv, x, p, xU)..., p
end
@inline function cc_tanh(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return tanh(x), sech(x)^2, p)
  (xU <= 0.0) && (return dline_seg(tanh, tanh_deriv, x, xL, xU)..., p)
  if p === Inf
    p, flag = secant(0.0, xU, 0.0, xU, tanh_env, xL, 0.0)
    if flag
      p = golden_section(0.0, xU, tanh_env, xL, 0.0)
    end
   end
   (x <= p) && (return dline_seg(tanh, tanh_deriv, x, xL, p)..., p)
   return tanh(x), sech(x)^2, p
end

@inline atan_deriv(x::Float64, y::Float64, z::Float64) = 1.0/(1.0+x^2)
@inline atan_env(x::Float64, y::Float64, z::Float64) = (x-y)-(1.0+sqr(x))*(atan(x)-atan(y))
@inline function cv_atan(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return dline_seg(atan, atan_deriv, x, xL, xU)..., p)
  (xU <= 0.0) && (return atan(x), 1.0/(1.0+x^2), p)
  if p === Inf
    p,flag = secant(xL, 0.0, xL, 0.0, atan_env, xU, 0.0)
    if flag
        p = golden_section(xL, 0.0, atan_env, xU, 0.0)
    end
  end
  (x <= p) && (return atan(x), 1.0/(1.0+x^2), p)
  return dline_seg(atan, atan_deriv, x, p, xU)..., p
end
@inline function cc_atan(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return atan(x), 1.0/(1.0+x^2), p)
  (xU <= 0.0) && (return dline_seg(atan, atan_deriv, x, xL, xU)..., p)
  if p === Inf
    p,flag = secant(0.0,xU,0.0,xU,atan_env,xL,0.0)
    if flag
      p = golden_section(0.0,xU,atan_env,xL,0.0)
    end
  end
  (x <= p) && (return dline_seg(atan, atan_deriv, x, xL, p)..., p)
  return atan(x), 1.0/(1.0+x^2), p
end

@inline asin_deriv(x::Float64, y::Float64, z::Float64) = 1.0/sqrt(1.0-x^2)
@inline asin_env(x::Float64, y::Float64, z::Float64) = (asin(x) - asin(y))*sqrt(1.0-x^2) - x + y
@inline function cv_asin(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return asin(x), 1.0/sqrt(1.0-x^2), p)
  (xU <= 0.0) && (return dline_seg(asin, asin_deriv, x, xL, xU)..., p)
  if p === Inf
    p,flag = secant(0.0, xU, 0.0, xU, asin_env, xL, 0.0)
    if flag
      p = golden_section(0.0, xU, asin_env, xL, 0.0)
    end
  end
  (x <= p) && (return dline_seg(asin, asin_deriv, x, xL, p)..., p)
  return asin(x), 1.0/sqrt(1.0-x^2, p)
end
@inline function cc_asin(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return dline_seg(asin, asin_deriv, x, xL, xU)..., p)
  (xU <= 0.0) && (return asin(x), 1.0/sqrt(1.0-x^2), p)
  if p === Inf
    p, flag = secant(xL, 0.0, xL, 0.0, asin_env, xU, 0.0)
    if flag
        p = golden_section(xL, 0.0, asin_env, xU, 0.0)
    end
  end
  (x <= p) && (return asin(x), 1.0/sqrt(1.0-x^2), p)
  return dline_seg(asin, asin_deriv, x, p, xU)..., p
end

@inline tan_deriv(x::Float64, y::Float64, z::Float64) = sec(x)^2
@inline tan_env(x::Float64, y::Float64, z::Float64) = (x-y)-(tan(x)-tan(y))/(1.0+tan(x)^2)
@inline tan_envd(x::Float64, y::Float64, z::Float64)= 2*tan(x)/(1.0+tan(x)^2)*(tan(x)-tan(y))
@inline function cv_tan(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return tan(x), sec(x)^2, p)
  (xU <= 0.0) && (return dline_seg(tan, tan_deriv, x, xL, xU)..., p)
  if p === Inf
    p, flag = secant(0.0, xU, 0.0, xU, tan_env, xL, 0.0)
    if flag
        p = golden_section(0.0, xU, tan_env, xL, 0.0)
    end
  end
  (x <= p) && (return dline_seg(tan, tan_deriv, x, xL, p)..., p)
  return tan(x), sec(x)^2, p
end
@inline function cc_tan(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return dline_seg(tan, tan_deriv, x, xL, xU)..., p)
  (xU <= 0.0) && (return tan(x), sec(x)^2, p)
  if p === Inf
    p, flag = secant(0.0, xL, xL, 0.0, tan_env, xU, 0.0)
    if flag
      p = golden_section(xL, 0.0, tan_env, xU, 0.0)
    end
  end
  (x <= p) && (return tan(x), sec(x)^2, p)
  return dline_seg(tan, tan_deriv, x, p, xU)..., p
end

@inline acos_deriv(x::Float64, y::Float64, z::Float64) = -1.0/sqrt(1.0-x^2)
@inline acos_env(x::Float64, y::Float64, z::Float64) = -(acos(x) - acos(y))*sqrt(1-x^2) - x + y
@inline function cv_acos(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return acos(x), -1.0/sqrt(1.0-x^2), p)
  (xU <= 0.0) && (return dline_seg(acos, acos_deriv, x, xL, xU)..., p)
  if p === Inf
    p, flag = secant(0.0, xU, 0.0, xU, acos_env, xL, 0.0)
    if flag
      p = golden_section(0.0, xU, acos_env, xL, 0.0)
    end
  end
  (x <= p) && (return dline_seg(acos, acos_deriv, x, xL, p)..., p)
  return acos(x), -1.0/sqrt(1.0-x^2), p
end
@inline function cc_acos(x::Float64, xL::Float64, xU::Float64, p::Float64)
    (xL >= 0.0) && (return dline_seg(acos, acos_deriv, x, xL, xU)..., p)
    (xU <= 0.0) && (return acos(x), -1.0/sqrt(1.0-x^2), p)
    if p === Inf
        p, flag = secant(xL, 0.0, xL, 0.0, acos_env, xU, 0.0)
        if flag
            p = golden_section(xL, 0.0, acos_env, xU, 0.0)
        end
      end
    (x <= p) && (return acos(x), -1.0/sqrt(1.0-x^2), p)
    return dline_seg(acos, acos_deriv, x, p, xU)..., p
end

@inline asinh_deriv(x::Float64, y::Float64, z::Float64) = 1.0/sqrt(1.0+x^2)
@inline asinh_env(x::Float64, y::Float64, z::Float64) = (x-y)-sqrt(1.0+sqr(x))*(asinh(x)-asinh(y))
@inline function cv_asinh(x::Float64, xL::Float64, xU::Float64, p::Float64)
    (xL >= 0.0) && (return dline_seg(asinh, asinh_deriv, x, xL, xU)..., p)
    (xU <= 0.0) && (return asinh(x), 1.0/sqrt(1.0 + x^2), p)
    if p === Inf
        p, flag = secant(xL, 0.0, xL, 0.0, asinh_env, xU, 0.0)
        if flag
            p = golden_section(xL, 0.0, asinh_env, xU, 0.0)
        end
    end
    (x <= p) && (return asinh(x), 1.0/sqrt(1.0+x^2), p)
    return dline_seg(asinh, asinh_deriv, x, p, xU)..., p
end
@inline function cc_asinh(x::Float64, xL::Float64, xU::Float64, p::Float64)
  (xL >= 0.0) && (return asinh(x), 1.0/sqrt(1.0 + x^2), p)
  (xU <= 0.0) && (return dline_seg(asinh, asinh_deriv, x, xL, xU)..., p)
  if p === Inf
    p, flag = secant(0.0, xU, 0.0, xU, asinh_env, xL, 0.0)
    if flag
      p = golden_section(0.0, xU, asinh_env, xL, 0.0)
    end
  end
  (x <= p) && (return dline_seg(asinh, asinh_deriv, x, xL, p)..., p)
  return asinh(x), 1.0/sqrt(1.0+x^2), p
end

# basic method overloading operator (sinh, tanh, atanh, asinh), convexoconcave or concavoconvex
eps_min_dict = Dict{Symbol,Symbol}(:sinh => :xL, :tanh => :xL, :asinh => :xL,
                                 :atanh => :xL, :tan => :xL, :acos => :xU,
                                 :asin => :xL, :atan => :xL)
eps_max_dict = Dict{Symbol,Symbol}(:sinh => :xU, :tanh => :xU, :asinh => :xU,
                                 :atanh => :xU, :tan => :xU, :acos => :xL,
                                 :asin => :xU, :atan => :xU)
for expri in (:sinh, :tanh, :asinh, :atanh, :tan, :acos, :asin, :atan)
    expri_cv = Symbol("cv_"*String(expri))
    expri_cc = Symbol("cc_"*String(expri))
    expri_kernel = Symbol(String(expri)*"_kernel")
    eps_min = eps_min_dict[expri]
    eps_max = eps_max_dict[expri]
    @eval @inline function ($expri_kernel)(x::MC{N, NS}, y::Interval{Float64},
                            cv_p::Float64, cc_p::Float64) where N
        if (y.lo == -Inf) || (y.hi == Inf)
            error("Function unbounded on this domain")
        end
        xL = x.Intv.lo
        xU = x.Intv.hi
        if (x.cc >= x.cv >= $eps_min) || ($eps_min >= x.cv >= x.cc)
            midcv::Float64 = x.cv
            cv_id::Int = 1
        elseif (x.cv >= x.cc >= $eps_min) || ($eps_min >= x.cc >= x.cv)
            midcv = x.cc
            cv_id = 2
        else
            midcv = $eps_min
            cv_id = 3
        end
        if (x.cc >= x.cv >= $eps_max) || ($eps_max >= x.cv >= x.cc)
            midcc::Float64 = x.cv
            cc_id::Int = 1
        elseif (x.cv >= x.cc >= $eps_max) || ($eps_max >= x.cc >= x.cv)
            midcc = x.cc
            cc_id = 2
        else
            midcc = $eps_max
            cc_id = 3
        end
        cv, dcv, cv_p = $(expri_cv)(midcv, xL, xU, cv_p)
        cc, dcc, cc_p = $(expri_cc)(midcc, xL, xU, cc_p)
        (cv_id == 1) && (cv_grad = x.cc_grad*dcv)
        (cv_id == 2) && (cv_grad = x.cv_grad*dcv)
        (cv_id == 3) && (cv_grad = zeros(SVector{N,Float64}))
        (cc_id == 1) && (cc_grad = x.cc_grad*dcc)
        (cc_id == 2) && (cc_grad = x.cv_grad*dcc)
        (cc_id == 3) && (cc_grad = zeros(SVector{N,Float64}))
        cv, cc, cv_grad, cc_grad = cut(y.lo, y.hi, cv, cc, cv_grad, cc_grad)
        return MC{N, NS}(cv, cc, y, cv_grad, cc_grad, x.cnst), cv_p, cc_p
    end
    @eval @inline function ($expri_kernel)(x::MC{N, Diff}, y::Interval{Float64},
                            cv_p::Float64, cc_p::Float64) where N
        if (y.lo == -Inf) || (y.hi == Inf)
            error("Function unbounded on this domain")
        end
        xL = x.Intv.lo
        xU = x.Intv.hi
        if (x.cc >= x.cv >= $eps_min) || ($eps_min >= x.cv >= x.cc)
            midcv::Float64 = x.cv
            cv_id::Int64 = 1
        elseif (x.cv >= x.cc >= $eps_min) || ($eps_min >= x.cc >= x.cv)
            midcv = x.cc
            cv_id = 2
        else
            midcv = $eps_min
            cv_id = 3
        end
        if (x.cc >= x.cv >= $eps_max) || ($eps_max >= x.cv >= x.cc)
            midcc::Float64 = x.cv
            cc_id::Int64 = 1
        elseif (x.cv >= x.cc >= $eps_max) || ($eps_max >= x.cc >= x.cv)
            midcc = x.cc
            cc_id = 2
        else
            midcc = $eps_max
            cc_id = 3
        end
        cv, dcv, cv_p = $(expri_cv)(midcv, xL, xU, cv_p)
        cc, dcc, cc_p = $(expri_cc)(midcc, xL, xU, cc_p)
        gcv1, gdcv1, cv_p = $(expri_cv)(x.cv, xL, xU, cv_p)
        gcc1, gdcc1, cc_p = $(expri_cc)(x.cv, xL, xU, cc_p)
        gcv2, gdcv2, cv_p = $(expri_cv)(x.cc, xL, xU, cv_p)
        gcc2, gdcc2, cc_p = $(expri_cc)(x.cc, xL, xU, cc_p)
        cv_grad = max(0.0, gdcv1)*x.cv_grad + min(0.0, gdcv2)*x.cc_grad
        cc_grad = min(0.0, gdcc1)*x.cv_grad + max(0.0, gdcc2)*x.cc_grad
        return MC{N,Diff}(cv, cc, y, cv_grad, cc_grad, x.cnst), cv_p, cc_p
    end
    @eval @inline function ($expri)(x::MC)
        z, tp1, tp2 = ($expri_kernel)(x, ($expri)(x.Intv), Inf, Inf)
        return z
    end
end


# cosh convex
@inline cv_cosh(x::Float64, xL::Float64, xU::Float64) = cosh(x), sinh(x)
@inline cc_cosh(x::Float64, xL::Float64, xU::Float64) = dline_seg(cosh, sinh, x, xL, xU)
@inline function cosh_kernel(x::MC{N, NS}, yintv::Interval{Float64}) where N
  if (yintv.lo == -Inf) || (yintv.hi == Inf)
    error("Function unbounded on this domain")
  end
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = yintv.lo
  xUc = yintv.hi
  eps_max = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.hi : x.Intv.lo
  if (x.Intv.lo < 0.0 < x.Intv.hi)
    eps_min = 0.0
  else
    eps_min = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.lo : x.Intv.hi
  end
  midcc,cc_id = mid3(x.cc, x.cv, eps_max)
  midcv,cv_id = mid3(x.cc, x.cv, eps_min)
  cc,dcc = cc_cosh(midcc, x.Intv.lo, x.Intv.hi)
  cv,dcv = cv_cosh(midcv, x.Intv.lo, x.Intv.hi)
  cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
  cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  cv, cc, cv_grad, cc_grad = cut(xLc, xUc, cv, cc, cv_grad, cc_grad)
  return MC{N}(cv, cc, yintv, cv_grad, cc_grad, x.cnst)
end
@inline function cosh_kernel(x::MC{N,Diff}, yintv::Interval{Float64}) where N
  if (yintv.lo == -Inf) || (yintv.hi == Inf)
    error("Function unbounded on this domain")
  end
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = yintv.lo
  xUc = yintv.hi
  eps_max = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.hi : x.Intv.lo
  if (x.Intv.lo < 0.0 < x.Intv.hi)
    eps_min = 0.0
  else
    eps_min = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.lo : x.Intv.hi
  end
  midcc,cc_id = mid3(x.cc, x.cv, eps_max)
  midcv,cv_id = mid3(x.cc, x.cv, eps_min)
  cc,dcc = cc_cosh(midcc, x.Intv.lo, x.Intv.hi)
  cv,dcv = cv_cosh(midcv, x.Intv.lo, x.Intv.hi)
  gcc1,gdcc1 = cc_cosh(x.cv, x.Intv.lo, x.Intv.hi)
  gcv1,gdcv1 = cv_cosh(x.cv, x.Intv.lo, x.Intv.hi)
  gcc2,gdcc2 = cc_cosh(x.cc, x.Intv.lo, x.Intv.hi)
  gcv2,gdcv2 = cv_cosh(x.cc, x.Intv.lo, x.Intv.hi)
  cv_grad = max(0.0, gdcv1)*x.cv_grad + min(0.0, gdcv2)*x.cc_grad
  cc_grad = min(0.0, gdcc1)*x.cv_grad + max(0.0, gdcc2)*x.cc_grad
  return MC{N}(cv, cc, yintv, cv_grad, cc_grad, x.cnst)
end
@inline cosh(x::MC) = cosh_kernel(x, cosh(x.Intv))

@inline deg2rad(x::MC, y::Interval{Float64}) = mult_kernel(x, pi/180.0, y)
@inline rad2deg(x::MC, y::Interval{Float64}) = mult_kernel(x, 180.0/pi, y)
@inline deg2rad(x::MC) = deg2rad(x, (pi/180.0)*x.Intv)
@inline rad2deg(x::MC) = rad2deg(x, (180.0/pi)*x.Intv)

# TODO: ADD efficient kernels for below (if applicable)
@inline sec(x::MC)= inv(cos(x))
@inline csc(x::MC)= inv(sin(x))
@inline cot(x::MC)= inv(tan(x))

@inline asec(x::MC) = acos(inv(x))
@inline acsc(x::MC) = asin(inv(x))
@inline acot(x::MC) = atan(inv(x))

@inline sech(x::MC) = inv(cosh(x))
@inline csch(x::MC) = inv(sinh(x))
@inline coth(x::MC) = inv(tanh(x))

@inline acsch(x::MC) = log(sqrt(1.0+inv(sqr(x)))+inv(x))
@inline acoth(x::MC) = 0.5*(log(1.0+inv(x))-log(1.0-inv(x)))

@inline sind(x::MC) = sin(deg2rad(x))
@inline cosd(x::MC) = cos(deg2rad(x))
@inline tand(x::MC) = tan(deg2rad(x))
@inline secd(x::MC) = inv(cosd(x))
@inline cscd(x::MC) = inv(sind(x))
@inline cotd(x::MC) = inv(tand(x))

@inline asind(x::MC) = rad2deg(asin(x))
@inline acosd(x::MC) = rad2deg(acos(x))
@inline atand(x::MC) = rad2deg(atan(x))
@inline asecd(x::MC) = rad2deg(asec(x))
@inline acscd(x::MC) = rad2deg(acsc(x))
@inline acotd(x::MC) = rad2deg(acot(x))
