########### Defines differentiable step relaxations
@inline function cv_step(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  if (xU<=zero(T))
    return zero(T),zero(T)
  elseif (xL>=zero(xL))
    return one(T),zero(T)
  elseif (x>=zero(x))
    return (x/xU)^2,2*x/xU^2
  else
    return zero(T),zero(T)
  end
end
@inline function cc_step(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  if (xU<=zero(T))
     return zero(T),zero(T)
  elseif (xL>=zero(T))
     return one(T),zero(T)
  elseif (x>=zero(T))
     return one(T),zero(T)
  else
    return one(T)-(x/xL)^2,-2*x/xL^2
  end
end
@inline function cv_step_NS(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  if (xU<=zero(T))
    return zero(T),zero(T)
  elseif (xL>=zero(T))
    return one(T),zero(T)
  elseif (x>=zero(T))
    return line_seg(x,zero(T),zero(T),xU,one(T)),dline_seg(x,zero(T),zero(T),xU,one(T),1/x)
  else
    return zero(T),zero(T)
  end
end
@inline function cc_step_NS(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  if (xU<=zero(T))
     return zero(T),zero(T)
  elseif (xL>=zero(T))
     return one(T),zero(T)
  elseif (x>=zero(x))
     return one(T),zero(T)
  else
    return line_seg(x,xL,zero(T),zero(T),one(T)),dline_seg(x,xL,zero(T),zero(T),one(T),1/x)
  end
end
@inline function step(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = step(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T = x.Intv.hi
  eps_min::T = x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  if (MC_param.mu >= 1)
    cc::T,dcc::T = cc_step(midcc,x.Intv.lo,x.Intv.hi)
    cv::T,dcv::T = cv_step(midcv,x.Intv.lo,x.Intv.hi)
    gcc1::T,gdcc1::T = cc_step(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1::T,gdcv1::T = cv_step(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2::T,gdcc2::T = cc_step(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2::T,gdcv2::T = cv_step(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
	  cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc,dcc = cc_step_NS(midcc,x.Intv.lo,x.Intv.hi)
    cv,dcv = cv_step_NS(midcv,x.Intv.lo,x.Intv.hi)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst)
end

function step(x::IntervalArithmetic.Interval{T}) where {T}
      isempty(x) && return emptyinterval(x)
      xmin::T = ((x.lo)<zero(T)) ? zero(T) : one(T)
      xmax::T = ((x.hi)>=zero(T)) ? one(T) : zero(T)
      return Interval{T}(xmin,xmax)
end

########### Defines sign
@inline sign(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = -step(-x)+step(x)
