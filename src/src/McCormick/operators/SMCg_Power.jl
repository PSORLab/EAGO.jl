# convex/concave relaxation (Khan 3.1-3.2) of integer powers of 1/x for positive reals
function cv_negpowpos(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  return x^n,n*x^(n-1)
end
function cc_negpowpos(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
end

# convex/concave relaxation of integer powers of 1/x for negative reals
@inline function cv_negpowneg(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  if (isodd(n))
    return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
  else
    return x^n,n*x^(n-1)
  end
end
@inline function cc_negpowneg(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  if (isodd(n))
    return x^n,n*x^(n-1)
  else
    return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
  end
end
# convex/concave relaxation of even powers greater than or equal to 4
@inline function cv_pow4(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  return x^n,n*x^(n-1)
end
@inline function cc_pow4(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
end
# convex/concave relaxation of odd powers
@inline function cv_powodd(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
          if (xU<=zero(T))
             return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
          elseif (zero(T)<=xL)
            return x^n,n*x^(n-1)
          else
            val::T = (xL^n)*(xU-x)/(xU-xL)+(max(zero(T),x))^n
            dval::T = -(xL^n)/(xU-xL)+n*(max(zero(T),x))^(n-1)
            return val,dval
          end
end
@inline function cc_powodd(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  if (xU<=zero(T))
    return x^n,n*x^(n-1)
  elseif (zero(T)<=xL)
    return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
  else
    val::T = (xU^n)*(x-xL)/(xU-xL)+(min(zero(T),x))^n
    dval::T = (xU^n)/(xU-xL)+n*(min(zero(T),x))^(n-1)
    return val,dval
  end
end
function pos_odd(x::SMCg{N,V,T},c::Integer) where {N,V,T<:AbstractFloat}
  intv::V = pow(x.Intv,c)
  xLc::T = intv.lo
  xUc::T = intv.hi
  if (MC_param.mu >= 1)
    eps_max::T = x.Intv.hi
    eps_min::T = x.Intv.lo
    midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
    midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
    cc,dcc::T = cc_powodd(midcc,x.Intv.lo,x.Intv.hi,c)
    cv,dcv::T = cv_powodd(midcv,x.Intv.lo,x.Intv.hi,c)
    gcc1::T,gdcc1::T = cc_powodd(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1::T,gdcv1::T = cv_powodd(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2::T,gdcc2::T = cc_powodd(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2::T,gdcv2::T = cv_powodd(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    eps_max = x.Intv.hi
    eps_min = x.Intv.lo
    midcc,cc_id = mid3(x.cc,x.cv,eps_max)
    midcv,cv_id = mid3(x.cc,x.cv,eps_min)
    cc,dcc = cc_powodd(midcc,x.Intv.lo,x.Intv.hi,c)
    cv,dcv = cv_powodd(midcv,x.Intv.lo,x.Intv.hi,c)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, intv,x.cnst,x.IntvBox,x.xref)
end

function pos_even(x::SMCg{N,V,T},c::Integer) where {N,V,T<:AbstractFloat}
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  intv::V = x.Intv^c
  xLc::T = intv.lo
  xUc::T = intv.hi
  if (x.Intv.hi<zero(x.cc))
    eps_min = x.Intv.hi
    eps_max = x.Intv.lo
  elseif (x.Intv.lo>zero(x.cc))
    eps_min = x.Intv.lo
    eps_max = x.Intv.hi
  else
    eps_min = zero(x.cc)
    eps_max = (abs(x.Intv.lo)>=abs(x.Intv.hi))? x.Intv.lo : x.Intv.hi
  end
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_pow4(midcc,x.Intv.lo,x.Intv.hi,c)
  cv,dcv = cv_pow4(midcv,x.Intv.lo,x.Intv.hi,c)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_pow4(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_pow4(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_pow4(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_pow4(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    cc_grad = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, intv,x.cnst,x.IntvBox,x.xref)
end

# TO DO: CUT CORRECTED
function neg_powneg_odd(x::SMCg{N,V,T},c::Integer) where {N,V,T<:AbstractFloat}
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  INTV::V = x.Intv^c
  xLc::T = INTV.lo
  xUc::T = INTV.hi
  eps_max = INTV.hi
  eps_min = INTV.lo
  if (MC_param.mu >= 1)
    midcc,cc_id = mid3(x.cc,x.cv,eps_max)
    midcv,cv_id = mid3(x.cc,x.cv,eps_min)
    cc,dcc = cc_negpowneg(midcc,x.Intv.lo,x.Intv.hi,c)
    cv,dcv = cv_negpowneg(midcv,x.Intv.lo,x.Intv.hi,c)
    gcc1,gdcc1 = cc_negpowneg(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_negpowneg(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_negpowneg(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_negpowneg(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    cc_grad = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
    return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, ((V<:AbstractMCInterval) ? V(xUc,xLc) : x.Intv^c), x.cnst,x.IntvBox,x.xref)
  else
    # calc cc
    if (xL < x.cv)
      cc = x.cv^c
      cc_grad = (c*x.cv^(c-1))*x.cv_grad
    elseif (xL > x.cc)
      cc = x.cc^c
      cc_grad = (c*x.cc^(c-1))*x.cc_grad
    else
      cc = xL^c
      cc_grad = zero(SVector{N,T})
    end

    if (xU == xL)
      cv = xLc
      cv_grad = zero(SVector{N,T})
      cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
      return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, ((V<:AbstractMCInterval) ? V(xUc,xLc) : x.Intv^c),x.cnst,x.IntvBox,x.xref)
    else
      dcv = (xU^c-xL^c)/(xU-xL) # function decreasing
      if (xU < x.cv)
        cv = xUc + dcv*(x.cv - xL)
        cv_grad = dcv*x.cv_grad
      elseif (xU > x.cc)
        cv = xUc + dcv*(x.cc - xL)
        cv_grad = dcv*x.cc_grad
      else
        cv = xUc + dcv*(xU - xL)
        cv_grad = zero(SVector{N,T})
      end
      cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
      return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, ((V<:AbstractMCInterval) ? V(xUc,xLc) : x.Intv^c),x.cnst,x.IntvBox,x.xref)
    end
  end
end

# TO DO: lATER, CUT CORRECTED
function neg_powneg_even(x::SMCg{N,V,T},c::Integer) where {N,V,T<:AbstractFloat}
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  intv = x.Intv^c
  xLc::T = intv.lo
  xUc::T = intv.hi
  eps_min = x.Intv.lo
  eps_max = x.Intv.hi
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_negpowneg(midcc,x.Intv.lo,x.Intv.hi,c)
  cv,dcv = cv_negpowneg(midcv,x.Intv.lo,x.Intv.hi,c)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_negpowneg(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_negpowneg(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_negpowneg(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_negpowneg(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    cc_grad = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, intv,x.cnst,x.IntvBox,x.xref)
end

# TO DO: CHECK CUT CORRECTED
function neg_powpos(x::SMCg{N,V,T},c::Int64) where {N,V,T<:AbstractFloat}
  xL::T = x.Intv.lo
  xU::T = x.Intv.hi
  xLc::T = xL^c
  xUc::T = xU^c
  eps_min::T = x.Intv.hi
  eps_max::T = x.Intv.lo
  if (MC_param.mu >= 1)
    # midcv,cv_id = mid3(x.cc,x.cv,eps_min)
    if (eps_min < x.cv)
      cv::T = x.cv^c
    elseif (eps_min>x.cc)
      cv = x.cc^c
    else
      cv = eps_min^c
    end

    # midcc,cc_id = mid3(x.cc,x.cv,eps_max)
    m::T = (xUc-xLc)/(xU-xL)
    if (eps_max < x.cv)
      cc::T = xLc + m*(x.cv-xL)
    elseif (eps_max > x.cc)
      cc = xLc + m*(x.cc-xL)
    else
      cc = xLc + m*(eps_max-xL)
    end
    cv_grad::SVector{N,T} =(min(zero(T),c*x.cc^(c-1)))*x.cc_grad
    cc_grad::SVector{N,T} = m*x.cv_grad
    return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, ((V<:AbstractMCInterval) ? V(xUc,xLc) : x.Intv^c), x.cnst,x.IntvBox,x.xref)
  else
    xcvc::T = x.cv^c
    xccc::T = x.cc^c
    # cv calculation
    if (xU < x.cv)
      cv = xcvc
      cv_grad = (c*x.cv^(c-1))*x.cv_grad
    elseif (xU > x.cc)
      cv = xccc
      cv_grad = (c*x.cc^(c-1))*x.cc_grad
    else
      cv = xU^c
      cv_grad = zero(SVector{N,T})
    end
    # cc calculation
    if (xU == xL)
      cc = xU
      cc_grad = zero(SVector{N,T})
      cv,cc,cv_grad,cc_grad = cut(xUc,xLc,cv,cc,cv_grad,cc_grad)
      return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, ((V<:AbstractMCInterval) ? V(xUc,xLc) : x.Intv^c), x.cnst,x.IntvBox,x.xref)
    else
      dcc = (xUc-xLc)/(xU-xL)
      if (xL < x.cv)
        cc = xLc + dcc*(x.cv - xL)
        cc_grad = dcc*x.cv_grad
      elseif (xL > x.cc)
        cc = xLc + dcc*(x.cc - xL)
        cc_grad = dcc*x.cc_grad
      else
        cc = xLc
        cc_grad = zero(SVector{N,T})
      end
      cv,cc,cv_grad,cc_grad = cut(xUc,xLc,cv,cc,cv_grad,cc_grad)
      return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, ((V<:AbstractMCInterval) ? V(xUc,xLc) : x.Intv^c), x.cnst, x.IntvBox, x.xref)
    end
  end
end

# DONE EXCEPT FOR SUBFUNCTIONS (neg_powneg_even,pos_odd,pos_even,sqr)
function pow(x::SMCg{N,V,T},c::Q) where {N,V,Q<:Integer,T<:AbstractFloat}
  if (c==0)
    return one(x) # DONE
  elseif (c>0)
    if (c==2)
      return sqr(x)
    elseif (c==1)
      return x
    elseif (isodd(c))
      return pos_odd(x,c)
    else
      return pos_even(x,c)
    end
  else
    if (x.Intv.hi<zero(T))
      if (isodd(c))
        return neg_powneg_odd(x,c)
      else
        return neg_powneg_even(x,c) # PRIORITY 2
      end
    elseif (x.Intv.lo>zero(T))
      return neg_powpos(x,c) # MOSTLY DONE
    else
      error("Function unbounded on domain")
    end
  end
end
########### power of a generalized McCormick object raised to c
function  (^)(x::SMCg{N,V,T},c::Q) where {N,V,Q<:Integer,T<:AbstractFloat}
  pow(x,c)
end
function  (^)(x::SMCg{N,V,T},c::Float64) where {N,V,T<:AbstractFloat}
  if (isinteger(c))
    return pow(x,Int64(c))
  else
    error("Relaxation of non-integer powers are unavailable.")
  end
end
########### Defines inverse
function inv(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  if (x.Intv.hi<zero(T))
    return neg_powneg_odd(x,-1)
  elseif (x.Intv.lo>zero(T))
    return neg_powpos(x,-1)
  else
    error("Function unbounded on domain")
  end
end
