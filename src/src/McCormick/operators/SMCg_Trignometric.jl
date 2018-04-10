# convex relaxation (envelope) of cos function
function cv_cos(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  r::T = zero(T)
  kL::Int64 = Base.ceil((-0.5-xL/(2.0*pi)))
  if (x<=(pi-2.0*pi*kL))
    xL1::T = convert(T,xL+2.0*pi*kL)
    if (xL1 >= pi/2.0)
      return cos(x),-sin(x)
    end
    xU1::T = convert(T,min(xU+2.0*pi*kL,pi))
    if ((xL1>=(-pi/2))&&(xU1<=(pi/2)))
      if (abs(xL-xU)<MC_param.env_tol)
        r = zero(T)
      else
        r = (cos(xU)-cos(xL))/(xU-xL)
      end
      return cos(xL)+r*(x-xL),r
    end
    return cv_cosin(convert(T,x+(2.0*pi)*kL),xL1,xU1)
  end
  kU::Int64 = Base.floor((0.5-xU/(2.0*pi)))
  if (x>=(-pi-2.0*pi*kU))
    xU2::T = convert(T,xU+2.0*pi*kU)
    if (xU2<=-pi/2.0)
      return cos(x),-sin(x)
    end
    return cv_cosin(convert(T,x+2.0*pi*kU),convert(T,max(xL+2.0*pi*kU,-pi)),xU2)
  end
  return -one(T),zero(T)
end
# function for computing convex relaxation over nonconvex and nonconcave regions
function cv_cosin(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  xj::T = -Inf
  if (abs(xL)<=abs(xU))
    left::Bool = false
    x0::T = xU
    xm::T = xL
  else
    left = true
    x0 = xL
    xm = xU
  end
  try
    xj = newton(x0,xL,xU,cv_cosenv,dcv_cosenv,xm,zero(x0))
  catch e
    if isa(e, ErrorException)
      xj = golden_section(xL,xU,cv_cosenv,xm,zero(x0))
    end
  end
  if ((left && x<=xj)||((~left) && x>=xj))
    return cos(x),-sin(x)
  else
    if abs(xm-xj)<MC_param.env_tol
      r::T = zero(T)
    else
      r = (cos(xm)-cos(xj))/(xm-xj)
    end
    return cos(xm)+r*(x-xm),r
  end
end
# pivot point calculation function for convex relaxation of cosine
function cv_cosenv(x::T,y::T,z::T) where {T<:AbstractFloat}
  return (x-y)*sin(x)+cos(x)-cos(y)
end
function dcv_cosenv(x::T,y::T,z::T) where {T<:AbstractFloat}
  return (x-y)*cos(x)
end
# concave relaxation (envelope) of cos function
function cc_cos(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  temp = cv_cos(convert(T,x-pi),convert(T,xL-pi),convert(T,xU-pi))
  return -temp[1],-temp[2]
end
function cos_arg(xL::T,xU::T) where {T<:AbstractFloat}
  kL::Int64 = Base.ceil(-0.5-xL/(2.0*pi))
  xL1::T = convert(T,xL+2.0*pi*kL)
  xU1::T = convert(T,xU+2.0*pi*kL)
  if ~((xL1>=-pi)&&(xL1<=pi))
    error("Cosine Argument Calculation: xL out of bounds.")
  end
  if (xL1<=zero(T))
    if (xU1<=zero(T))
      arg1::T = xL
      arg2::T = xU
    elseif (xU1>=pi)
      arg1 = convert(T,pi-2.0*pi*kL)
      arg2 = convert(T,-2.0*pi*kL)
    else
      arg1 = (cos(xL1)<=cos(xU1)) ? xL : xU
      arg2 = convert(T,-2.0*pi*kL)
    end
  end
  if (xU1<=pi)
    arg1 = xU
    arg2 = xL
  elseif (xU1>=(2.0*pi))
    arg1 = convert(T,pi-2.0*pi*kL)
    arg2 = convert(T,2.0*pi-2.0*pi*kL)
  else
    arg1 = convert(T,pi-2.0*pi*kL)
    arg2 = (cos(xL1)>=cos(xU1)) ? xL : xU
  end
  return [arg1,arg2]
end

function cos(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = cos(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T,eps_min::T = cos_arg(x.Intv.lo,x.Intv.hi)
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = cc_cos(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = cv_cos(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu>=1)
    gcc1::T,gdcc1::T = cc_cos(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1::T,gdcv1::T = cv_cos(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2::T,gdcc2::T = cc_cos(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2::T,gdcv2::T = cv_cos(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst, x.IntvBox, x.xref)
end

function sin(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  cos(x-pi/2.0)
end

# pivot point calculation function for convex relaxation of complimentary error function
@inline function tan_env(x::T,y::T,z::T) where {T<:AbstractFloat}
  return (x-y)-(tan(x)-tan(y))/(one(T)+tan(x)^2)
end
# derivative of pivot point calculation function for convex relaxation of complimentary error function
@inline function tan_envd(x::T,y::T,z::T) where {T<:AbstractFloat}
  return 2*tan(x)/(one(T)+tan(x)^2)*(tan(x)-tan(y))
end
# convex relaxation (envelope) of tangent function
@inline function cv_tan(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  p::T = zero(T)
  if (xL>=zero(T))
    return tan(x),sec(x)^2
  elseif (xU<=zero(T))
    return line_seg(x,xL,tan(xL),xU,tan(xU)),dline_seg(x,xL,tan(xL),xU,tan(xU),sec(x)^2)
  else
    try
      p = secant(zero(T),xU,zero(T),xU,tan_env,xL,zero(T))
    catch e
      if isa(e, ErrorException)
        p = golden_section(zero(T),xU,tan_env,xL,zero(T))
      end
    end
    if (x<=p)
      return line_seg(x,xL,tan(xL),p,tan(p)),dline_seg(x,xL,tan(xL),p,tan(p),sec(x)^2)
    else
      return tan(x),sec(x)^2
    end
  end
end
# concave relaxation (envelope) of tangent function
@inline function cc_tan(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  p::T = zero(T)
  if (xL>=zero(T))
    return line_seg(x,xL,tan(xL),xU,tan(xU)),dline_seg(x,xL,tan(xL),xU,tan(xU),sec(x)^2)
  elseif (xU<=zero(T))
    return tan(x),sec(x)^2
  else
    try
      p = secant(zero(T),xL,xL,zero(T),tan_env,xU,zero(T))
    catch e
      if isa(e, ErrorException)
        p = golden_section(xL,zero(T),tan_env,xU,zero(T))
      end
    end
    if (x<=p)
       return tan(x),sec(x)^2
    else
       return line_seg(x,p,tan(p),xU,tan(xU)),dline_seg(x,p,tan(p),xU,tan(xU),sec(x)^2)
     end
  end
end
@inline function tan(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = tan(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T = x.Intv.hi
  eps_min::T = x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = cc_tan(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = cv_tan(midcv,x.Intv.lo,x.Intv.hi)
  if ((x.Intv.lo==-Inf)||(x.Intv.hi==Inf))
    error("Function unbounded on domain")
  end
  if (MC_param.mu >= 1)
    gcc1::T,gdcc1::T = cc_tan(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1::T,gdcv1::T = cv_tan(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2::T,gdcc2::T = cc_tan(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2::T,gdcv2::T = cv_tan(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst, x.IntvBox, x.xref)
end

@inline function acos(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  return asin(-x)+pi/2.0
end

# pivot point calculation function for convex relaxation of arcsine
@inline function asin_env(x::T,y::T,z::T) where {T<:AbstractFloat}
  return (asin(x)-asin(y))/(x-y)-one(T)/sqrt(one(T)-x^2)
end
# derivative of pivot point calculation function for convex relaxation of arcsine
@inline function asin_envd(x::T,y::T,z::T) where {T<:AbstractFloat}
  return one(T)/((x-y)*sqrt(one(T)-x^2))-x/(one(T)-x^2)^(3/2)-(asin(x)-asin(y))/(x-y)^2
end
# convex relaxation (envelope) of arcsine function
@inline function cv_asin(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  p::T = zero(T)
  if (xL>=zero(T))
    return asin(x),one(T)/sqrt(one(T)-x^2)
  elseif (xU<=zero(T))
    return line_seg(x,xL,asin(xL),xU,asin(xU)),dline_seg(x,xL,asin(xL),xU,asin(xU),one(T)/sqrt(one(T)-x^2))
  else
    try
      p = newton(xU,xL,xU,asin_env,asin_envd,xL,xU)
    catch e
      if isa(e, ErrorException)
        p = golden_section(zero(x),xU,asin_env,xL,zero(T))
      end
    end
    if (x<=p)
      return line_seg(x,xL,asin(xL),p,asin(p)),dline_seg(x,xL,asin(xL),p,asin(p),one(T)/sqrt(one(T)-x^2))
    else
      return asin(x),one(T)/sqrt(one(T)-x^2)
    end
  end
end
# concave relaxation (envelope) of arcsine function
@inline function cc_asin(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  p::T = zero(T)
  if (xL>=zero(T))
    return line_seg(x,xL,asin(xL),xU,asin(xU)),dline_seg(x,xL,asin(xL),xU,asin(xU),one(T)/sqrt(one(T)-x^2))
  elseif (xU<=zero(T))
    return asin(x),one(T)/sqrt(one(T)-x^2)
  else
    try
      p = secant(zero(T),xL,xL,zero(T),asin_env,xU,zero(T))
    catch e
      if isa(e, ErrorException)
        p = golden_section(xL,zero(T),asin_env,xU,zero(T))
      end
    end
    if (x<=p)
      return asin(x),one(T)/sqrt(one(T)-x^2)
    else
      return line_seg(x,p,asin(p),xU,asin(xU)),dline_seg(x,p,asin(p),xU,asin(xU),one(T)/sqrt(one(T)-x^2))
    end
  end
end
@inline function asin(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = asin(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T = x.Intv.hi
  eps_min::T = x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = cc_asin(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = cv_asin(midcv,x.Intv.lo,x.Intv.hi)
  if ((x.Intv.lo==-Inf)||(x.Intv.hi==Inf))
    error("Function unbounded on domain")
  end
  if (MC_param.mu >= 1)
    gcc1::T,gdcc1::T = cc_asin(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1::T,gdcv1::T = cv_asin(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2::T,gdcc2::T = cc_asin(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2::T,gdcv2::T = cv_asin(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst, x.IntvBox, x.xref)
end

# pivot point calculation function for convex relaxation of arctangent
@inline function atan_env(x::T,y::T,z::T) where {T<:AbstractFloat}
  return (x-y)-(one(T)+x^2)*(atan(x)-atan(y))
end
# derivative of pivot point calculation function for convex relaxation of arctangent
@inline function atan_envd(x::T,y::T,z::T) where {T<:AbstractFloat}
  return -2*x*(atan(x)-atan(y))
end
# convex relaxation (envelope) of arctan function
@inline function cv_atan(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  p::T = zero(T)
  if (xL>=zero(T))
    return line_seg(x,xL,atan(xL),xU,atan(xU)),dline_seg(x,xL,atan(xL),xU,atan(xU),one(T)/(x^2+one(T)))
  elseif (xU<=zero(T))
    return atan(x),one(T)/(x^2+one(T))
  else
    try
      p = newton(xL,xL,zero(T),atan_env,atan_envd,xU,zero(T))
    catch e
      if isa(e, ErrorException)
        p = golden_section(xL,zero(T),atan_env,xU,zero(T))
      end
    end
    if (x<=p)
      return atan(x),one(T)/(x^2+one(T))
    else
      return line_seg(x,p,atan(p),xU,atan(xU)),dline_seg(x,p,atan(p),xU,atan(xU),one(T)/(x^2+one(T)))
    end
  end
end
# concave relaxation (envelope) of arctan function
@inline function cc_atan(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  p::T = zero(T)
  if (xL>=zero(T))
    return atan(x),one(T)/(x^2+one(T))
  elseif (xU<=zero(T))
    return line_seg(x,xL,atan(xL),xU,atan(xU)),dline_seg(x,xL,atan(xL),xU,atan(xU),one(T)/(x^2+one(T)))
  else
    try
      p = newton(xU,zero(T),xU,atan_env,atan_envd,xL,zero(T))
    catch e
      if isa(e, ErrorException)
        p = golden_section(zero(T),xU,atan_env,xL,zero(T))
      end
    end
    if (x<=p)
      return line_seg(x,xL,atan(xL),p,atan(p)),dline_seg(x,xL,atan(xL),p,atan(p),one(T)/(x^2+one(T)))
    else
      return atan(x),one(T)/(x^2+one(T))
    end
  end
end
@inline function atan(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = atan(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T = x.Intv.hi
  eps_min::T = x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = cc_atan(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = cv_atan(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1::T,gdcc1::T = cc_atan(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1::T,gdcv1::T = cv_atan(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2::T,gdcc2::T = cc_atan(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2::T,gdcv2::T = cv_atan(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst, x.IntvBox, x.xref)
end
