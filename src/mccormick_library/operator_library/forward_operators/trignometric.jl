# convex relaxation (envelope) of cos function
function cv_cos(x,xL,xU)
  r = 0.0
  kL = Base.ceil((-0.5-xL/(2.0*pi)))
  if (x<=(pi-2.0*pi*kL))
    xL1 = xL+2.0*pi*kL
    if (xL1 >= pi/2.0)
      return cos(x),-sin(x)
    end
    xU1 = min(xU+2.0*pi*kL,pi)
    if ((xL1>=(-pi/2))&&(xU1<=(pi/2)))
      if (abs(xL-xU)<MC_param.env_tol)
        r = 0.0
      else
        r = (cos(xU)-cos(xL))/(xU-xL)
      end
      return cos(xL)+r*(x-xL),r
    end
    return cv_cosin(x+(2.0*pi)*kL,xL1,xU1)
  end
  kU = Base.floor((0.5-xU/(2.0*pi)))
  if (x>=(-pi-2.0*pi*kU))
    xU2 = xU+2.0*pi*kU
    if (xU2<=-pi/2.0)
      return cos(x),-sin(x)
    end
    return cv_cosin(x+2.0*pi*kU,max(xL+2.0*pi*kU,-pi),xU2)
  end
  return -1.0,0.0
end
# function for computing convex relaxation over nonconvex and nonconcave regions
function cv_cosin(x,xL,xU)
  xj = -Inf
  if (abs(xL)<=abs(xU))
    left = false
    x0 = xU
    xm = xL
  else
    left = true
    x0 = xL
    xm = xU
  end
  xj,flag = newton(x0,xL,xU,cv_cosenv,dcv_cosenv,xm,0.0)
  if flag
      xj = golden_section(xL,xU,cv_cosenv,xm,0.0)
  end
  if ((left && x<=xj)||((~left) && x>=xj))
    return cos(x),-sin(x)
  else
    if abs(xm-xj)<MC_param.env_tol
      r = 0.0
    else
      r = (cos(xm)-cos(xj))/(xm-xj)
    end
    return cos(xm)+r*(x-xm),r
  end
end

# pivot point calculation function for convex relaxation of cosine
cv_cosenv(x,y,z) = (x-y)*sin(x)+cos(x)-cos(y)
dcv_cosenv(x,y,z) = (x-y)*cos(x)

# concave relaxation (envelope) of cos function
function cc_cos(x,xL,xU)
  temp = cv_cos(x-pi,xL-pi,xU-pi)
  return -temp[1],-temp[2]
end
function cos_arg(xL,xU)
  kL = Base.ceil(-0.5-xL/(2.0*pi))
  xL1 = xL+2.0*pi*kL
  xU1 = xU+2.0*pi*kL
  if ~((xL1>=-pi)&&(xL1<=pi))
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

function cos(x::MC{N}) where N
  Intv = cos(x.Intv)
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = Intv.lo
  xUc = Intv.hi
  eps_max,eps_min = cos_arg(x.Intv.lo,x.Intv.hi)
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_cos(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_cos(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu>=1)
    gcc1,gdcc1 = cc_cos(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1,gdcv1 = cv_cos(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2,gdcc2 = cc_cos(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2,gdcv2 = cv_cos(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

sin(x::MC) = cos(x-pi/2.0)

# pivot point calculation function for convex relaxation of complimentary error function
tan_env(x,y,z) = (x-y)-(tan(x)-tan(y))/(1.0+tan(x)^2)
# derivative of pivot point calculation function for convex relaxation of complimentary error function
tan_envd(x,y,z)= 2*tan(x)/(1.0+tan(x)^2)*(tan(x)-tan(y))

# convex relaxation (envelope) of tangent function
function cv_tan(x,xL,xU)
  p = 0.0
  if (xL>=0.0)
    return tan(x),sec(x)^2
  elseif (xU<=0.0)
    return line_seg(x,xL,tan(xL),xU,tan(xU)),dline_seg(x,xL,tan(xL),xU,tan(xU),sec(x)^2)
  else
    p,flag = secant(0.0,xU,0.0,xU,tan_env,xL,0.0)
    if flag
        p = golden_section(0.0,xU,tan_env,xL,0.0)
    end
    if (x<=p)
      return line_seg(x,xL,tan(xL),p,tan(p)),dline_seg(x,xL,tan(xL),p,tan(p),sec(x)^2)
    else
      return tan(x),sec(x)^2
    end
  end
end
# concave relaxation (envelope) of tangent function
function cc_tan(x,xL,xU)
  p = 0.0
  if (xL>=0.0)
    return line_seg(x,xL,tan(xL),xU,tan(xU)),dline_seg(x,xL,tan(xL),xU,tan(xU),sec(x)^2)
  elseif (xU<=0.0)
    return tan(x),sec(x)^2
  else
    p,flag = secant(0.0,xL,xL,0.0,tan_env,xU,0.0)
    if flag
      p = golden_section(xL,0.0,tan_env,xU,0.0)
    end
    if (x<=p)
       return tan(x),sec(x)^2
    else
       return line_seg(x,p,tan(p),xU,tan(xU)),dline_seg(x,p,tan(p),xU,tan(xU),sec(x)^2)
     end
  end
end
function tan(x::MC{N}) where N
  Intv = tan(x.Intv)
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = Intv.lo
  xUc = Intv.hi
  midcc,cc_id = mid3(x.cc,x.cv,xU)
  midcv,cv_id = mid3(x.cc,x.cv,xL)
  cc,dcc = cc_tan(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_tan(midcv,x.Intv.lo,x.Intv.hi)
  if ((Intv.lo==-Inf)||(Intv.hi==Inf))
    error("Function unbounded on domain")
  end
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_tan(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1,gdcv1 = cv_tan(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2,gdcc2 = cc_tan(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2,gdcv2 = cv_tan(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

######

acos_env(x::Float64, y::Float64, z::Float64) = -(acos(x) - acos(y))*sqrt(1-x^2) - x + y
function cv_acos(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return acos(x), -1.0/sqrt(1.0-x^2)
  elseif (xU <= 0.0)
    return line_seg(x, xL, acos(xL), xU, acos(xU)), dline_seg(x, xL, acos(xL), xU, acos(xU), -1.0/sqrt(1.0-x^2))
  else
    p,flag = secant(0.0, xU, 0.0, xU, acos_env, xL, 0.0)
    if flag
      p = golden_section(0.0, xU, acos_env, xL, 0.0)
    end
    if (x <= p)
       return line_seg(x, xL, acos(xL), p, acos(p)), dline_seg(x, xL, acos(xL), p, acos(p), -1.0/sqrt(1.0-x^2))
    else
       return acos(x), -1.0/sqrt(1.0-x^2)
    end
  end
end
function cc_acos(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return line_seg(x, xL, acos(xL), xU, acos(xU)), dline_seg(x, xL, acos(xL), xU, acos(xU), -1.0/sqrt(1.0-x^2))
  elseif (xU <= 0.0)
    return acos(x), -1.0/sqrt(1.0-x^2)
  else
    p,flag = secant(xL, 0.0, xL, 0.0, acos_env, xU, 0.0)
    if flag
        p = golden_section(xL, 0.0, acos_env, xU, 0.0)
    end
    if (x <= p)
      return acos(x), -1.0/sqrt(1.0-x^2)
    else
      return line_seg(x, p, acos(p), xU, acos(xU)), dline_seg(x, p, acos(p), xU, acos(xU), -1.0/sqrt(1.0-x^2))
    end
  end
end
function acos(x::MC{N}) where N
  Intv = acos(x.Intv)
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = Intv.lo
  xUc = Intv.hi
  midcc,cc_id = mid3(x.cc,x.cv,xU)
  midcv,cv_id = mid3(x.cc,x.cv,xL)
  cc,dcc = cc_acos(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_acos(midcv,x.Intv.lo,x.Intv.hi)
  if ((Intv.lo==-Inf)||(Intv.hi==Inf))
    error("Function unbounded on domain")
  end
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_acos(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1,gdcv1 = cv_acos(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2,gdcc2 = cc_acos(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2,gdcv2 = cv_acos(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

######

asin_env(x::Float64, y::Float64, z::Float64) = (asin(x) - asin(y))*sqrt(1-x^2) - x + y
function cv_asin(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return asin(x), 1.0/sqrt(1.0-x^2)
  elseif (xU <= 0.0)
    return line_seg(x, xL, asin(xL), xU, asin(xU)), dline_seg(x, xL, asin(xL), xU, asin(xU), 1.0/sqrt(1.0-x^2))
  else
    p,flag = secant(0.0, xU, 0.0, xU, asin_env, xL, 0.0)
    if flag
      p = golden_section(0.0, xU, asin_env, xL, 0.0)
    end
    if (x <= p)
       return line_seg(x, xL, asin(xL), p, asin(p)), dline_seg(x, xL, asin(xL), p, asin(p), 1.0/sqrt(1.0-x^2))
    else
       return asin(x), 1.0/sqrt(1.0-x^2)
    end
  end
end
function cc_asin(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return line_seg(x, xL, asin(xL), xU, asin(xU)), dline_seg(x, xL, asin(xL), xU, asin(xU), 1.0/sqrt(1.0-x^2))
  elseif (xU <= 0.0)
    return asin(x), 1.0/sqrt(1.0-x^2)
  else
    p,flag = secant(xL, 0.0, xL, 0.0, asin_env, xU, 0.0)
    if flag
        p = golden_section(xL, 0.0, asin_env, xU, 0.0)
    end
    if (x <= p)
      return asin(x), 1.0/sqrt(1.0-x^2)
    else
      return line_seg(x, p, asin(p), xU, asin(xU)), dline_seg(x, p, asin(p), xU, asin(xU), 1.0/sqrt(1.0-x^2))
    end
  end
end
function asin(x::MC{N}) where N
  Intv = asin(x.Intv)
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = Intv.lo
  xUc = Intv.hi
  midcc,cc_id = mid3(x.cc,x.cv,xU)
  midcv,cv_id = mid3(x.cc,x.cv,xL)
  cc,dcc = cc_asin(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_asin(midcv,x.Intv.lo,x.Intv.hi)
  if ((Intv.lo==-Inf)||(Intv.hi==Inf))
    error("Function unbounded on domain")
  end
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_asin(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1,gdcv1 = cv_asin(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2,gdcc2 = cc_asin(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2,gdcv2 = cv_asin(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

####

atan_env(x::Float64, y::Float64, z::Float64) = (x-y)-(1.0+sqr(x))*(atan(x)-atan(y))
function cv_atan(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return line_seg(x, xL, atan(xL), xU, atan(xU)), dline_seg(x, xL, atan(xL), xU, atan(xU), 1.0/(1.0+x^2))
  elseif (xU <= 0.0)
    return atan(x), 1.0/(1.0+x^2)
  else
    p,flag = secant(xL, 0.0, xL, 0.0, atan_env, xU, 0.0)
    if flag
        p = golden_section(xL, 0.0, atan_env, xU, 0.0)
    end
    if (x <= p)
      return atan(x), 1.0/(1.0+x^2)
    else
      return line_seg(x, p, atan(p), xU, atan(xU)), dline_seg(x, p, atan(p), xU, atan(xU), 1.0/(1.0+x^2))
    end
  end
end
function cc_atan(x::Float64, xL::Float64, xU::Float64)
  p = 0.0
  if (xL >= 0.0)
    return atan(x), 1.0/(1.0+x^2)
  elseif (xU <= 0.0)
    return line_seg(x, xL, atan(xL), xU, atan(xU)), dline_seg(x, xL, atan(xL), xU, atan(xU), 1.0/(1.0+x^2))
  else
    p,flag = secant(0.0,xU,0.0,xU,atan_env,xL,0.0)
    if flag
      p = golden_section(0.0,xU,atan_env,xL,0.0)
    end
    if (x <= p)
       return line_seg(x, xL, atan(xL), p, atan(p)), dline_seg(x, xL, atan(xL), p, atan(p), 1.0/(1.0+x^2))
    else
       return atan(x), 1.0/(1.0+x^2)
     end
  end
end
function atan(x::MC{N}) where N
    Intv = atan(x.Intv)
    xL = x.Intv.lo
    xU = x.Intv.hi
    xLc = Intv.lo
    xUc = Intv.hi
    midcc,cc_id = mid3(x.cc,x.cv,xU)
    midcv,cv_id = mid3(x.cc,x.cv,xL)
    cc,dcc = cc_atan(midcc,x.Intv.lo,x.Intv.hi)
    cv,dcv = cv_atan(midcv,x.Intv.lo,x.Intv.hi)
    if ((Intv.lo==-Inf)||(Intv.hi==Inf))
      error("Function unbounded on domain")
    end
    if (MC_param.mu >= 1)
      gcc1,gdcc1 = cc_atan(x.cv,x.Intv.lo,x.Intv.hi)
      gcv1,gdcv1 = cv_atan(x.cv,x.Intv.lo,x.Intv.hi)
      gcc2,gdcc2 = cc_atan(x.cc,x.Intv.lo,x.Intv.hi)
      gcv2,gdcv2 = cv_atan(x.cc,x.Intv.lo,x.Intv.hi)
      cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
      cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
    else
      cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
      cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
      cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
    end
    return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end


deg2rad(x::MC) = (pi/180.0)*x
rad2deg(x::MC) = (180.0/pi)*x

sec(x::MC)= inv(cos(x))
csc(x::MC)= inv(sin(x))
cot(x::MC)= inv(tan(x))

asec(x::MC) = acos(inv(x))
acsc(x::MC) = asin(inv(x))
acot(x::MC) = atan(inv(x))

sech(x::MC) = inv(cosh(x))
csch(x::MC) = inv(sinh(x))
coth(x::MC) = inv(tanh(x))

acsch(x::MC) = log(sqrt(1.0+inv(sqr(x)))+inv(x))
acoth(x::MC) = 0.5*(log(1.0+inv(x))-log(1.0-inv(x)))

sind(x::MC) = sin(deg2rad(x))
cosd(x::MC) = cos(deg2rad(x))
tand(x::MC) = tan(deg2rad(x))
secd(x::MC) = inv(cosd(x))
cscd(x::MC) = inv(sind(x))
cotd(x::MC) = inv(tand(x))

asind(x::MC) = rad2deg(asin(x))
acosd(x::MC) = rad2deg(acos(x))
atand(x::MC) = rad2deg(atan(x))
asecd(x::MC) = rad2deg(asec(x))
acscd(x::MC) = rad2deg(acsc(x))
acotd(x::MC) = rad2deg(acot(x))
