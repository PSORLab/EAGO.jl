function exp(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}

  # initialize reused terms
  Intv::V = exp(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  expxL::T = Intv.lo
  expxU::T = Intv.hi

  if (MC_param.mu >= 1)
    if (x.cc <= xL)
      cv::T = exp(x.cc)
      cv_grad::SVector{N,T} = exp(x.cv)*x.cv_grad
    elseif (xL >= x.cv)
      cv = expxL
      cv_grad = exp(x.cv)*x.cv_grad
    else
      cv = exp(x.cv)
      cv_grad = cv*x.cv_grad
    end
    # calc concave term
    if (neq(xU,xL))
      dcc::T = (expxU-expxL)/(xU-xL)
      cc_grad::SVector{N,T} = dcc*x.cc_grad
      if (x.cc <= xU)
        cc::T = expxL + dcc*(x.cc - xL)
      elseif (xU >= x.cv)
        cc = expxL + dcc*(xU - xL)
      else
        cc = expxL + dcc*(x.cv - xL)
      end
    else
      cc = expxU
      cc_grad = zeros(SVector{N,T})
    end
  else

    if (x.cc <= xL)
      cv = exp(x.cc)
      cv_grad = cv*x.cc_grad
    elseif (xL >= x.cv)
      cv = expxL
      cv_grad = zeros(SVector{N,T})
    else
      cv = exp(x.cv)
      cv_grad = cv*x.cv_grad
    end

    # calc concave term
    if (neq(xU,xL))
      dcc = (expxU-expxL)/(xU-xL)
      if (x.cc <= xU)
        cc = expxL + dcc*(x.cc - xL)
        cc_grad = dcc*x.cc_grad
      elseif (xU >= x.cv)
        cc = expxL + dcc*(xU - xL)
        cc_grad = zeros(SVector{N,T})
      else
        cc = expxL + dcc*(x.cv - xL)
        cc_grad = dcc*x.cv_grad
      end
    else
      cc = expxU
      cc_grad = zeros(SVector{N,T})
    end

    # applies cut operator
    cv,cc,cv_grad,cc_grad = cut(expxL,expxU,cv,cc,cv_grad,cc_grad)
  end

  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, exp(x.Intv),x.cnst, x.IntvBox, x.xref)
end

function exp2_cc(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return line_seg(x,xL,exp2(xL),xU,exp2(xU)),dline_seg(x,xL,exp2(xL),xU,exp2(xU),exp2(x)*log(2))
end
function exp2_cv(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return exp2(x),exp2(x)*log(2)
end
function exp2(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = exp2(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T = x.Intv.hi
  eps_min::T = x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = exp2_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = exp2_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1::T,gdcc1::T = exp2_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1::T,gdcv1::T = exp2_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2::T,gdcc2::T = exp2_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2::T,gdcv2::T = exp2_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
	  cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst,x.IntvBox, x.xref)
end

function exp10_cc(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return line_seg(x,xL,exp10(xL),xU,exp10(xU)),dline_seg(x,xL,exp10(xL),xU,exp10(xU),exp10(x)*log(10))
end
function exp10_cv(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return exp10(x),exp10(x)*log(10)
end
function exp10(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = exp10(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T = x.Intv.hi
  eps_min::T = x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = exp10_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = exp10_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1::T,gdcc1::T = exp10_cc(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1::T,gdcv1::T = exp10_cv(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2::T,gdcc2::T = exp10_cc(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2::T,gdcv2::T = exp10_cv(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst,x.IntvBox, x.xref)
end

function sqrt_cc(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return sqrt(x),one(T)/(2*sqrt(x))
end
function sqrt_cv(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return line_seg(x,xL,sqrt(xL),xU,sqrt(xU)),dline_seg(x,xL,sqrt(xL),xU,sqrt(xU),one(T)/(2*sqrt(x)))
end
function sqrt(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = sqrt(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T = x.Intv.hi
  eps_min::T = x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = sqrt_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = sqrt_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1::T,gdcc1::T = sqrt_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1::T,gdcv1::T = sqrt_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2::T,gdcc2::T = sqrt_cc(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2::T,gdcv2::T = sqrt_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
	  cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst,x.IntvBox, x.xref)
end

function log_cc(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return log(x),one(T)/x
end
function log_cv(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return line_seg(x,xL,log(xL),xU,log(xU)),dline_seg(x,xL,log(xL),xU,log(xU),one(T)/x)
end
function log(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = log(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T = x.Intv.hi
  eps_min::T = x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = log_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = log_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = log_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1,gdcv1 = log_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2,gdcc2 = log_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2,gdcv2 = log_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst,x.IntvBox, x.xref)
end

function log2_cc(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return log2(x),one(T)/(log(2)*x)
end
function log2_cv(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return line_seg(x,xL,log2(xL),xU,log2(xU)),dline_seg(x,xL,log2(xL),xU,log2(xU),one(T)/(log(2)*x))
end
function log2(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = log2(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T = x.Intv.hi
  eps_min::T = x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = log2_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = log2_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1::T,gdcc1::T = log2_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1::T,gdcv1::T = log2_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2::T,gdcc2::T = log2_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2::T,gdcv2::T = log2_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
	  cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst,x.IntvBox, x.xref)
end

function log10_cc(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return log10(x),one(T)/(log(10)*x)
end
function log10_cv(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return line_seg(x,xL,log10(xL),xU,log10(xU)),dline_seg(x,xL,log10(xL),xU,log10(xU),one(T)/(log(10)*x))
end
function log10(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = log10(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T = x.Intv.hi
  eps_min::T = x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = log10_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = log10_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1::T,gdcc1::T = log10_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1::T,gdcv1::T = log10_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2::T,gdcc2::T = log10_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2::T,gdcv2::T = log10_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
	  cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst,x.IntvBox, x.xref)
end

function acosh_cc(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return acosh(x),one(T)/sqrt(x^2 - one(T))
end
function acosh_cv(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return line_seg(x,xL,acosh(xL),xU,acosh(xU)),dline_seg(x,xL,acosh(xL),xU,acosh(xU),one(T)/sqrt(x^2 - one(T)))
end
function acosh(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = acosh(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_max::T = x.Intv.hi
  eps_min::T = x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = acosh_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = acosh_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1::T,gdcc1::T = acosh_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1::T,gdcv1::T = acosh_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2::T,gdcc2::T = acosh_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2::T,gdcv2::T = acosh_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
	  cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst,x.IntvBox, x.xref)
end

function abs_cv(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  if (x>=zero(x))
    return xU*(x/xU)^(MC_param.mu+1),(MC_param.mu+1)*(x/xU)^MC_param.mu
  else
    return -xL*(x/xL)^(MC_param.mu+1), -(MC_param.mu+1)*(x/xL)^MC_param.mu
  end
end
function abs_cc(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  if (x>=zero(x))
    d = (MC_param.mu+1)*(x/xU)^MC_param.mu
  else
    d = -(MC_param.mu+1)*(x/xL)^MC_param.mu
  end
  return line_seg(x,xL,abs(xL),xU,abs(xU)),dline_seg(x,xL,abs(xL),xU,abs(xU),d)
end
function abs_cc_NS(x::T,lo::T,hi::T) where {T<:AbstractFloat}
  return line_seg(x,lo,abs(lo),hi,abs(hi)),dline_seg(x,lo,abs(lo),hi,abs(hi),sign(x))
end
function abs_cv_NS(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  return abs(x),sign(x)
end
function abs(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = abs(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_min::T,blank = mid3(x.Intv.lo,x.Intv.hi,zero(x.Intv.lo))
  eps_max::T = (abs(x.Intv.hi)>=abs(x.Intv.lo)) ? x.Intv.hi : x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  if (MC_param.mu >= 1)
    cc::T,dcc::T = abs_cc(midcc,x.Intv.lo,x.Intv.hi)
    cv::T,dcv::T = abs_cv(midcv,x.Intv.lo,x.Intv.hi)
    gcc1::T,gdcc1::T = abs_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1::T,gdcv1::T = abs_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2::T,gdcc2::T = abs_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2::T,gdcv2::T = abs_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
	  cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    #println("abs x in: ", x)
    cc,dcc = abs_cc_NS(midcc,x.Intv.lo,x.Intv.hi)
    cv,dcv = abs_cv_NS(midcv,x.Intv.lo,x.Intv.hi)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst ,x.IntvBox, x.xref)
end

function cosh_cv(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  cosh(x),-sinh(x)
end
function cosh_cc(x::T,xL::T,xU::T) where {T<:AbstractFloat}
  line_seg(x,xL,cosh(xL),xU,cosh(xU)),dline_seg(x,xL,cosh(xL),xU,cosh(xU),-sinh(x))
end
function cosh(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  Intv::V = cosh(x.Intv)
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  eps_min::T,blank = mid3(x.Intv.lo,x.Intv.hi,zero(x.Intv.lo))
  eps_max::T = (abs(x.Intv.hi)>=abs(x.Intv.lo)) ? x.Intv.hi : x.Intv.lo
  midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
  midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  cc::T,dcc::T = cosh_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv::T,dcv::T = cosh_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1::T,gdcc1::T = cosh_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1::T,gdcv1::T = cosh_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2::T,gdcc2::T = cosh_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2::T,gdcv2::T = cosh_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
	  cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv, x.cnst ,x.IntvBox, x.xref)
end
