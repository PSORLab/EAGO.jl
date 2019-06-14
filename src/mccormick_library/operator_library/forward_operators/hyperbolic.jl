# sinh concave... convex
sinh_env(x::Float64, y::Float64, z::Float64) = (x-y)*cosh(x)-(sinh(x)-sinh(y))
function cv_sinh(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return line_seg(x, xL, sinh(xL), xU, sinh(xU)), dline_seg(x, xL, sinh(xL), xU, sinh(xU), cosh(x))
  elseif (xU <= 0.0)
    return sinh(x), cosh(x)
  else
    p,flag = secant(xL, 0.0, xL, 0.0, sinh_env, xU, 0.0)
    if flag
        p = golden_section(xL, 0.0, sinh_env, xU, 0.0)
    end
    if (x <= p)
      return sinh(x), cosh(x)
    else
      return line_seg(x, p, sinh(p), xU, sinh(xU)), dline_seg(x, p, sinh(p), xU, sinh(xU), cosh(x))
    end
  end
end
function cc_sinh(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return sinh(x), cosh(x)
  elseif (xU <= 0.0)
    return line_seg(x, xL, sinh(xL), xU, sinh(xU)), dline_seg(x, xL, sinh(xL), xU, sinh(xU), cosh(x))
  else
    p,flag = secant(0.0,xU,0.0,xU,sinh_env,xL,0.0)
    if flag
      p = golden_section(0.0,xU,sinh_env,xL,0.0)
    end
    if (x <= p)
       return line_seg(x, xL, sinh(xL), p, sinh(p)), dline_seg(x, xL, sinh(xL), p, sinh(p), cosh(x))
    else
       return sinh(x), cosh(x)
     end
  end
end
function sinh(x::MC{N}) where N
    Intv = sinh(x.Intv)
    xL = x.Intv.lo
    xU = x.Intv.hi
    xLc = Intv.lo
    xUc = Intv.hi
    midcc,cc_id = mid3(x.cc,x.cv,xU)
    midcv,cv_id = mid3(x.cc,x.cv,xL)
    cc,dcc = cc_sinh(midcc,x.Intv.lo,x.Intv.hi)
    cv,dcv = cv_sinh(midcv,x.Intv.lo,x.Intv.hi)
    if ((Intv.lo==-Inf)||(Intv.hi==Inf))
      error("Function unbounded on domain")
    end
    if (MC_param.mu >= 1)
      gcc1,gdcc1 = cc_sinh(x.cv,x.Intv.lo,x.Intv.hi)
      gcv1,gdcv1 = cv_sinh(x.cv,x.Intv.lo,x.Intv.hi)
      gcc2,gdcc2 = cc_sinh(x.cc,x.Intv.lo,x.Intv.hi)
      gcv2,gdcv2 = cv_sinh(x.cc,x.Intv.lo,x.Intv.hi)
      cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
      cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
    else
      cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
      cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
      cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
    end
    return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

# cosh convex
cv_cosh(x::Float64, xL::Float64, xU::Float64) = cosh(x), sinh(x)
function cc_cosh(x::Float64, xL::Float64, xU::Float64)
  fxL = cosh(xL)
  fxU = cosh(xU)
  line_seg(x, xL, fxL, xU, fxU), dline_seg(x, xL, fxL, xU, fxU, sinh(x))
end
function cosh(x::MC{N}) where N
  Intv = cosh(x.Intv)
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = Intv.lo
  xUc = Intv.hi
  eps_max = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.hi : x.Intv.lo
  if (x.Intv.lo < 0.0 < x.Intv.hi)
    eps_min = 0.0
  else
    eps_min = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.lo : x.Intv.hi
  end
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_cosh(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_cosh(midcv,x.Intv.lo,x.Intv.hi)
  if ((Intv.lo==-Inf)||(Intv.hi==Inf))
    error("Function unbounded on domain")
  end
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_cosh(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1,gdcv1 = cv_cosh(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2,gdcc2 = cc_cosh(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2,gdcv2 = cv_cosh(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

#####
# tanh convex... concave
tanh_env(x::Float64, y::Float64, z::Float64) = (x-y)-(tanh(x)-tanh(y))/(1.0-tanh(x)^2)
tanh_envd(x::Float64, y::Float64, z::Float64)= 2*tanh(x)/(1.0-tanh(x)^2)*(tanh(x)-tanh(y))
function cv_tanh(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL>=0.0)
    return line_seg(x,xL,tanh(xL),xU,tanh(xU)),dline_seg(x,xL,tanh(xL),xU,tanh(xU),sech(x)^2)
  elseif (xU<=0.0)
    return tanh(x),sech(x)^2
  else
    p,flag = secant(xL,0.0,xL,0.0,tanh_env,xU,0.0)
    if flag
        p = golden_section(xL,0.0,tanh_env,xU,0.0)
    end
    if (x<=p)
      return tanh(x),sech(x)^2
    else
      return line_seg(x,p,tanh(p),xU,tanh(xU)),dline_seg(x,p,tanh(p),xU,tanh(xU),sech(x)^2)
    end
  end
end
function cc_tanh(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return tanh(x),sech(x)^2
  elseif (xU <= 0.0)
    return line_seg(x,xL,tanh(xL),xU,tanh(xU)),dline_seg(x,xL,tanh(xL),xU,tanh(xU),sech(x)^2)
  else
    p,flag = secant(0.0,xU,0.0,xU,tanh_env,xL,0.0)
    if flag
      p = golden_section(0.0,xU,tanh_env,xL,0.0)
    end
    if (x <= p)
       return line_seg(x,xL,tanh(xL),p,tanh(p)),dline_seg(x,xL,tanh(xL),p,tanh(p),sech(x)^2)
    else
       return tanh(x),sech(x)^2
     end
  end
end
function tanh(x::MC{N}) where N
    Intv = tanh(x.Intv)
    xL = x.Intv.lo
    xU = x.Intv.hi
    xLc = Intv.lo
    xUc = Intv.hi
    midcc,cc_id = mid3(x.cc,x.cv,xU)
    midcv,cv_id = mid3(x.cc,x.cv,xL)
    cc,dcc = cc_tanh(midcc,x.Intv.lo,x.Intv.hi)
    cv,dcv = cv_tanh(midcv,x.Intv.lo,x.Intv.hi)
    if ((Intv.lo==-Inf)||(Intv.hi==Inf))
      error("Function unbounded on domain")
    end
    if (MC_param.mu >= 1)
      gcc1,gdcc1 = cc_tanh(x.cv,x.Intv.lo,x.Intv.hi)
      gcv1,gdcv1 = cv_tanh(x.cv,x.Intv.lo,x.Intv.hi)
      gcc2,gdcc2 = cc_tanh(x.cc,x.Intv.lo,x.Intv.hi)
      gcv2,gdcv2 = cv_tanh(x.cc,x.Intv.lo,x.Intv.hi)
      cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
      cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
    else
      cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
      cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
      cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
    end
    return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

# asinh concave... convex
asinh_env(x::Float64, y::Float64, z::Float64) = (x-y)-sqrt(1.0+sqr(x))*(asinh(x)-asinh(y))
function cv_asinh(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return line_seg(x, xL, asinh(xL), xU, asinh(xU)), dline_seg(x, xL, asinh(xL), xU, asinh(xU), 1.0/sqrt(1.0+x^2))
  elseif (xU <= 0.0)
    return asinh(x), 1.0/sqrt(1.0+x^2)
  else
    p,flag = secant(xL, 0.0, xL, 0.0, asinh_env, xU, 0.0)
    if flag
        p = golden_section(xL, 0.0, asinh_env, xU, 0.0)
    end
    if (x <= p)
      return asinh(x), 1.0/sqrt(1.0+x^2)
    else
      return line_seg(x, p, asinh(p), xU, asinh(xU)), dline_seg(x, p, asinh(p), xU, asinh(xU), 1.0/sqrt(1.0+x^2))
    end
  end
end
function cc_asinh(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return asinh(x), 1.0/sqrt(1.0+x^2)
  elseif (xU <= 0.0)
    return line_seg(x, xL, asinh(xL), xU, asinh(xU)), dline_seg(x, xL, asinh(xL), xU, asinh(xU), 1.0/sqrt(1.0+x^2))
  else
    p,flag = secant(0.0,xU,0.0,xU,asinh_env,xL,0.0)
    if flag
      p = golden_section(0.0,xU,asinh_env,xL,0.0)
    end
    if (x <= p)
       return line_seg(x, xL, asinh(xL), p, asinh(p)), dline_seg(x, xL, asinh(xL), p, asinh(p), 1.0/sqrt(1.0+x^2))
    else
       return asinh(x), 1.0/sqrt(1.0+x^2)
     end
  end
end
function asinh(x::MC{N}) where N
    Intv = asinh(x.Intv)
    xL = x.Intv.lo
    xU = x.Intv.hi
    xLc = Intv.lo
    xUc = Intv.hi
    midcc,cc_id = mid3(x.cc,x.cv,xU)
    midcv,cv_id = mid3(x.cc,x.cv,xL)
    cc,dcc = cc_asinh(midcc,x.Intv.lo,x.Intv.hi)
    cv,dcv = cv_asinh(midcv,x.Intv.lo,x.Intv.hi)
    if ((Intv.lo==-Inf)||(Intv.hi==Inf))
      error("Function unbounded on domain")
    end
    if (MC_param.mu >= 1)
      gcc1,gdcc1 = cc_asinh(x.cv,x.Intv.lo,x.Intv.hi)
      gcv1,gdcv1 = cv_asinh(x.cv,x.Intv.lo,x.Intv.hi)
      gcc2,gdcc2 = cc_asinh(x.cc,x.Intv.lo,x.Intv.hi)
      gcv2,gdcv2 = cv_asinh(x.cc,x.Intv.lo,x.Intv.hi)
      cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
      cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
    else
      cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
      cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
      cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
    end
    return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

# atanh convex... concave
atanh_env(x::Float64, y::Float64, z::Float64) = (atanh(x) - atanh(y))*(1-x^2) - x + y
function cv_atanh(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return atanh(x), 1.0/(1.0-x^2)
  elseif (xU <= 0.0)
    return line_seg(x, xL, atanh(xL), xU, atanh(xU)), dline_seg(x, xL, atanh(xL), xU, atanh(xU), 1.0/(1.0-x^2))
  else
    p,flag = secant(0.0, xU, 0.0, xU, atanh_env, xL, 0.0)
    if flag
      p = golden_section(0.0, xU, atanh_env, xL, 0.0)
    end
    if (x <= p)
       return line_seg(x, xL, atanh(xL), p, atanh(p)), dline_seg(x, xL, atanh(xL), p, atanh(p), 1.0/(1.0-x^2))
    else
       return atanh(x), 1.0/(1.0-x^2)
    end
  end
end
function cc_atanh(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return line_seg(x, xL, atanh(xL), xU, atanh(xU)), dline_seg(x, xL, atanh(xL), xU, atanh(xU), 1.0/(1.0-x^2))
  elseif (xU <= 0.0)
    return atanh(x), 1.0/(1.0-x^2)
  else
    p,flag = secant(xL, 0.0, xL, 0.0, atanh_env, xU, 0.0)
    if flag
        p = golden_section(xL, 0.0, atanh_env, xU, 0.0)
    end
    if (x <= p)
      return atanh(x), 1.0/(1.0-x^2)
    else
      return line_seg(x, p, atanh(p), xU, atanh(xU)), dline_seg(x, p, atanh(p), xU, atanh(xU), 1.0/(1.0-x^2))
    end
  end
end
function atanh(x::MC{N}) where N
  Intv = atanh(x.Intv)
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = Intv.lo
  xUc = Intv.hi
  midcc,cc_id = mid3(x.cc,x.cv,xU)
  midcv,cv_id = mid3(x.cc,x.cv,xL)
  cc,dcc = cc_atanh(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_atanh(midcv,x.Intv.lo,x.Intv.hi)
  if ((Intv.lo==-Inf)||(Intv.hi==Inf))
    error("Function unbounded on domain")
  end
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_atanh(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1,gdcv1 = cv_atanh(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2,gdcc2 = cc_atanh(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2,gdcv2 = cv_atanh(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end
