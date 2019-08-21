#=

# defines square operator
sqr(x::T) where T = x*x
cv_sqr_NS(x,xL,xU) = x^2
dcv_sqr_NS(x,xL,xU) = 2.0*x

cc_sqr(x,xL,xU) = (xU>xL) ? xL^2 + (xL+xU)*(x-xL) : xU^2
dcc_sqr(x,xL,xU) = (xU>xL) ? (xL+xU) : 0.0
function cv_sqr(x,xL,xU)
  (0.0 <= xL || xU <= 0.0) && return x^2
	((xL < 0.0) && (0.0 <= xU) && (0.0 <= x)) && return (x^3)/xU
	return (x^3)/xL
end
function dcv_sqr(x,xL,xU)
  (0.0 <= xL || xU <= 0.0) && return 2.0*x
	((xL < 0.0) && (0.0 <= xU) && (0.0 <= x)) && (3.0*x^2)/xU
	return (3.0*x^2)/xL
end

cv_negpowpos(x,xL,xU,n) = x^n
dcv_negpowpos(x,xL,xU,n) = n*x^(n-1)
cc_negpowpos(x,xL,xU,n) = (xU == xL) ? xU : ((xU-x)/(xU-xL))*xL^n +((x-xL)/(xU-xL))*xU^n
dcc_negpowpos(x,xL,xU,n) = (xU == xL) ? 0.0 : (xU^n-xL^n)/(xU-xL)

cv_pow4(x,xL,xU,n) = x^n
dcv_pow4(x,xL,xU,n) = n*x^(n-1)
cc_pow4(x,xL,xU,n) = (xU == xL) ? xU : ((xU-x)/(xU-xL))*xL^n +((x-xL)/(xU-xL))*xU^n
dcc_pow4(x,xL,xU,n) = (xU == xL) ? 0.0 : (xU^n-xL^n)/(xU-xL)

# convex/concave relaxation of integer powers of 1/x for negative reals
cv_negpowneg(x,xL,xU,n::Int) = isodd(n) ? cc_negpowpos(x,xL,xU,n) : x^n
dcv_negpowneg(x,xL,xU,n::Int) = isodd(n) ? dcc_negpowpos(x,xL,xU,n) : n*x^(n-1)
cc_negpowneg(x,xL,xU,n::Int) = isodd(n) ? x^n : (xU == xL) ? xU : ((xU-x)/(xU-xL))*xL^n +((x-xL)/(xU-xL))*xU^n
dcc_negpowneg(x,xL,xU,n::Int) = isodd(n) ? n*x^(n-1) : (xU == xL) ? 0.0 : (xU^n-xL^n)/(xU-xL)

cv_negpowneg(x,xL,xU,n::Float64) = isodd(Int(n)) ? cc_negpowpos(x,xL,xU,n) : x^n
dcv_negpowneg(x,xL,xU,n::Float64) = isodd(Int(n)) ? dcc_negpowpos(x,xL,xU,n) : n*x^(n-1)
cc_negpowneg(x,xL,xU,n::Float64) = isodd(Int(n)) ? x^n : (xU == xL) ? xU : ((xU-x)/(xU-xL))*xL^n +((x-xL)/(xU-xL))*xU^n
dcc_negpowneg(x,xL,xU,n::Float64) = isodd(Int(n)) ? n*x^(n-1) : (xU == xL) ? 0.0 : (xU^n-xL^n)/(xU-xL)

function cv_powodd(x,xL,xU,n)
  (xU <= 0.0) && return cc_pow4(x,xL,xU,n)
  (0.0 <= xL) && return x^n
  return (xL^n)*(xU-x)/(xU-xL)+(max(0.0,x))^n
end
function dcv_powodd(x,xL,xU,n)
  (xU <= 0.0) && return dcc_pow4(x,xL,xU,n)
  (0.0 <= xL) && return n*x^(n-1)
  return -(xL^n)/(xU-xL)+n*(max(0.0,x))^(n-1)
end
function cc_powodd(x,xL,xU,n)
  println("ran 3")
  println("x: $x, xL: $xL, xU: $xU, c: $n")
  println("x^n: $(x^n)")
  (xU <= 0.0) && return x^n
  println("ran 1")
  (0.0 <= xL) && return cc_pow4(x,xL,xU,n)
  println("ran 2")
  return (xU^n)*(x-xL)/(xU-xL)+(min(0.0,x))^n
end
function dcc_powodd(x,xL,xU,n)
  (xU <= 0.0) && return (n-1)*x^n
  (0.0 <= xL) && return dcc_pow4(x,xL,xU,n)
  return (xU^n)/(xU-xL)+n*(min(0.0,x))^(n-1)
end

function cv_pow(x,xL,xU,c)
  if isinteger(c)
    if (c > 0)
     if (isodd(Int(c)))
        return cv_powodd(x,xL,xU,c)
      else
        return cv_pow4(x,xL,xU,c)
      end
    else
      if (xU < 0.0)
        if (isodd(Int(c)))
          return cv_powodd(x,xL,xU,c)
        else
          return cv_negpowneg(x,xL,xU,c)
        end
      elseif (xL > 0.0)
        return cv_negpowpos(x,xL,xU,c)
      else
        error("Function unbounded on domain")
      end
    end
  else
    if (xL >= 0.0)
      (c > 1.0) && return x^c
      (0.0 < c < 1.0) && return (xL^c)*(xU-x)/(xU-xL)+(max(0.0,x))^c # Concave
      (c < 0.0) && return x^c # Convex
    end
  end
end

function dcv_pow(x,xL,xU,c)
  if isinteger(c)
    if (c > 0)
     if (isodd(Int(c)))
        return dcv_powodd(x,xL,xU,c)
      else
        return dcv_pow4(x,xL,xU,c)
      end
    else
      if (xU < 0.0)
        if (isodd(Int(c)))
          return dcv_powodd(x,xL,xU,c)
        else
          return dcv_negpowneg(x,xL,xU,c)
        end
      elseif (xL > 0.0)
        return dcv_negpowpos(x,xL,xU,c)
      else
        error("Function unbounded on domain")
      end
    end
  else
    if (xL >= 0.0)
      (c > 1.0) && return c*x^(c-1.0)
      (0.0 < c < 1.0) && return -(xL^c)/(xU-xL)+c*(max(0.0,x))^(c-1.0) # Concave
      (c < 0.0) && return c*x^(c-1.0) # Convex
    end
  end
end

function cc_pow(x,xL,xU,c)
  if isinteger(c)
    if (c > 0)
     if (isodd(Int(c)))
		println("cc pow-pos odd")
        return cc_powodd(x,xL,xU,c)
      else
		 println("cc pow-pos even")
        return cc_pow4(x,xL,xU,c)
      end
    else
      if (xU < 0.0)
        if (isodd(Int(c)))
		  println("cc pow-neg dom-neg odd")
		  println("x: $x, xL: $xL, xU: $xU, c: $c")
          return cc_powodd(x,xL,xU,c)
        else
		  println("cc pow-neg dom-neg even")
          return cc_negpowneg(x,xL,xU,c)
        end
      elseif (xL > 0.0)
		 println("cc pow-neg dom-pos")
        return cc_negpowpos(x,xL,xU,c)
      else
        error("Function unbounded on domain")
      end
    end
  else
    if (xL >= 0.0)
      (c > 1.0) && return (xL^c)*(xU-x)/(xU-xL)+(max(0.0,x))^c
      (0.0 < c < 1.0) && return x^c # Concave
      (c < 0.0) && return (xL^c)*(xU-x)/(xU-xL)+(max(0.0,x))^c # Convex
    end
  end
end

function dcc_pow(x,xL,xU,c)
  if isinteger(c)
    if (c > 0)
     if (isodd(Int(c)))
        return dcc_powodd(x,xL,xU,c)
      else
        return dcc_pow4(x,xL,xU,c)
      end
    else
      if (xU < 0.0)
        if (isodd(Int(c)))
          return dcc_powodd(x,xL,xU,c)
        else
          return dcc_negpowneg(x,xL,xU,c)
        end
      elseif (xL > 0.0)
        return dcc_negpowpos(x,xL,xU,c)
      else
        error("Function unbounded on domain")
      end
    end
  else
    if (xL >= 0.0)
      (c > 1.0) && return -(xL^c)/(xU-xL)+c*(max(0.0,x))^(c-1.0)
      (0.0 < c < 1.0) && return c*x^(c-1.0) # Concave
      (c < 0.0) && return -(xL^c)/(xU-xL)+c*(max(0.0,x))^(c-1.0) # Convex
    end
  end
end

function sqr(x::MC{N}) where N
	Intv::IntervalType = x.Intv^2
  xL::Float64 = x.Intv.lo
  xU::Float64 = x.Intv.hi
  xLc::Float64 = Intv.lo
  xUc::Float64 = Intv.hi
  eps_max::Float64 = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.hi : x.Intv.lo
	if (x.Intv.lo < zero(Float64) < x.Intv.hi)
		eps_min::Float64 = zero(Float64)
	else
		eps_min = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.lo : x.Intv.hi
	end
  midcc::Float64,cc_id::Int = mid3(x.cc,x.cv,eps_max)
  midcv::Float64,cv_id::Int = mid3(x.cc,x.cv,eps_min)
  cc::Float64 = cc_sqr(midcc,x.Intv.lo,x.Intv.hi)
  dcc::Float64 = dcc_sqr(midcc,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
	   cv::Float64 = cv_sqr(midcv,x.Intv.lo,x.Intv.hi)
     dcv::Float64 = dcv_sqr(midcv,x.Intv.lo,x.Intv.hi)
     gdcc1::Float64 = dcc_sqr(x.cv,x.Intv.lo,x.Intv.hi)
     gdcv1::Float64 = dcv_sqr(x.cv,x.Intv.lo,x.Intv.hi)
     gdcc2::Float64 = dcc_sqr(x.cc,x.Intv.lo,x.Intv.hi)
     gdcv2::Float64 = dcv_sqr(x.cc,x.Intv.lo,x.Intv.hi)
     cv_grad::SVector{N,Float64} = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
     cc_grad::SVector{N,Float64} = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
	   cv = cv_sqr_NS(midcv,x.Intv.lo,x.Intv.hi)
     dcv = dcv_sqr_NS(midcv,x.Intv.lo,x.Intv.hi)
     cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
     cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
     cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

function pow_max(Intv::IntervalType,c::Float64)
  if isinteger(c)
    ci  = Int(c)
    if isodd(ci)
      return hi(Intv)
    elseif iseven(ci)
      if (abs(lo(Intv)) < abs(hi(Intv)))
        (lo(Intv) < 0.0) && (return 0.0)
        return (hi(Intv))
      else (abs(lo(Intv)) > abs(hi(Intv)))
        (hi(Intv) > 0.0) && (return 0.0)
        return (lo(Intv))
      end
    end
  end
  return hi(Intv)
end

function pow_min(Intv::IntervalType,c::Float64)
  if isinteger(c)
    ci  = Int(c)
    if isodd(ci)
      return lo(Intv)
    elseif iseven(ci)
      if (abs(lo(Intv)) < abs(hi(Intv)))
        (lo(Intv) < 0.0) && (return 0.0)
        return (lo(Intv))
      else (abs(lo(Intv)) > abs(hi(Intv)))
        (hi(Intv) > 0.0) && (return 0.0)
        return (hi(Intv))
      end
    end
  end
  return lo(Intv)
end

function pow(x::MC{N},c::Float64) where N
  println("x: $x")
  println("c: $c")
  if (c == 0)
	  println("ran zero arc")
    return one(x)
  elseif (c == 1)
	  println("ran one arc")
    return x
  elseif (c == 2)
	  println("ran two arc")
    return sqr(x)
  elseif (!isinteger(c) && lo(x.Intv) <= 0.0)
	  println("ran float arc")
    return exp(c*log(x))
  else
	println("ran other integer arc")
    Intv = x.Intv^c
    xL = x.Intv.lo
    xU = x.Intv.hi
    xLc = Intv.lo
    xUc = Intv.hi
    eps_max = pow_max(x.Intv,c)
    eps_min = pow_min(x.Intv,c)
    midcc,cc_id = mid3(x.cc,x.cv,eps_max)
    midcv,cv_id = mid3(x.cc,x.cv,eps_min)
    cc = cc_pow(midcc,x.Intv.lo,x.Intv.hi,c)
	println("cc asigned: $cc")
    cv = cv_pow(midcv,x.Intv.lo,x.Intv.hi,c)
	println("cv asigned: $cv")
    if (MC_param.mu >= 1)
       gdcc1 = dcc_pow(x.cv,x.Intv.lo,x.Intv.hi,c)
       gdcv1 = dcv_pow(x.cv,x.Intv.lo,x.Intv.hi,c)
       gdcc2 = dcc_pow(x.cc,x.Intv.lo,x.Intv.hi,c)
       gdcv2 = dcv_pow(x.cc,x.Intv.lo,x.Intv.hi,c)
       cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
       cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
    else
       dcv = dcv_pow(midcv,x.Intv.lo,x.Intv.hi,c)
       dcc = dcc_pow(midcc,x.Intv.lo,x.Intv.hi,c)
       cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
       cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
       cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
    end
	temp = MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
	println("temp: $temp")
    return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
  end
end


(^)(x::MC,c::Int) = pow(x,Float64(c))
pow(x::MC,c::Int) = pow(x,Float64(c))
inv(x::MC) = pow(x,-1)

=#

# defines square operator
sqr(x::T) where T = x*x
cv_sqr_NS(x,xL,xU) = x^2
dcv_sqr_NS(x,xL,xU) = 2.0*x

cc_sqr(x,xL,xU) = (xU>xL) ? xL^2 + (xL+xU)*(x-xL) : xU^2
dcc_sqr(x,xL,xU) = (xU>xL) ? (xL+xU) : 0.0
function cv_sqr(x,xL,xU)
  (0.0 <= xL || xU <= 0.0) && return x^2
	((xL < 0.0) && (0.0 <= xU) && (0.0 <= x)) && return (x^3)/xU
	return (x^3)/xL
end
function dcv_sqr(x,xL,xU)
  (0.0 <= xL || xU <= 0.0) && return 2.0*x
	((xL < 0.0) && (0.0 <= xU) && (0.0 <= x)) && (3.0*x^2)/xU
	return (3.0*x^2)/xL
end

function sqr(x::MC{N}) where N
	Intv = x.Intv^2
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
  cc = cc_sqr(midcc,x.Intv.lo,x.Intv.hi)
  dcc = dcc_sqr(midcc,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
	   cv = cv_sqr(midcv,x.Intv.lo,x.Intv.hi)
     dcv = dcv_sqr(midcv,x.Intv.lo,x.Intv.hi)
     gdcc1 = dcc_sqr(x.cv,x.Intv.lo,x.Intv.hi)
     gdcv1 = dcv_sqr(x.cv,x.Intv.lo,x.Intv.hi)
     gdcc2 = dcc_sqr(x.cc,x.Intv.lo,x.Intv.hi)
     gdcv2 = dcv_sqr(x.cc,x.Intv.lo,x.Intv.hi)
     cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
     cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
	   cv = cv_sqr_NS(midcv,x.Intv.lo,x.Intv.hi)
     dcv = dcv_sqr_NS(midcv,x.Intv.lo,x.Intv.hi)
     cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
     cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
     cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

# convex/concave relaxation (Khan 3.1-3.2) of integer powers of 1/x for positive reals
function cv_negpowpos(x,xL,xU,n::Integer)
  return x^n,n*x^(n-1)
end
function cc_negpowpos(x,xL,xU,n::Integer)
  return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
end

# convex/concave relaxation of integer powers of 1/x for negative reals
@inline function cv_negpowneg(x,xL,xU,n::Integer)
  if (isodd(n))
    return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
  else
    return x^n,n*x^(n-1)
  end
end
@inline function cc_negpowneg(x,xL,xU,n::Integer)
  if (isodd(n))
    return x^n,n*x^(n-1)
  else
    return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
  end
end
# convex/concave relaxation of even powers greater than or equal to 4
@inline function cv_pow4(x,xL,xU,n::Integer)
  return x^n,n*x^(n-1)
end
@inline function cc_pow4(x,xL,xU,n::Integer)
  return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
end
# convex/concave relaxation of odd powers
@inline function cv_powodd(x,xL,xU,n::Integer)
          if (xU <= 0.0)
             return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
          elseif (0.0 <= xL)
            return x^n,n*x^(n-1)
          else
            val = (xL^n)*(xU-x)/(xU-xL)+(max(0.0,x))^n
            dval = -(xL^n)/(xU-xL)+n*(max(0.0,x))^(n-1)
            return val,dval
          end
end
@inline function cc_powodd(x,xL,xU,n::Integer)
  if (xU <= 0.0)
    return x^n,n*x^(n-1)
  elseif (0.0 <= xL)
    return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
  else
    val = (xU^n)*(x-xL)/(xU-xL)+(min(0.0,x))^n
    dval = (xU^n)/(xU-xL)+n*(min(0.0,x))^(n-1)
    return val,dval
  end
end
function pos_odd(x::MC{N},c::Integer) where {N}
  intv = pow(x.Intv,c)
  xLc = intv.lo
  xUc = intv.hi
  if (MC_param.mu >= 1)
    eps_max = x.Intv.hi
    eps_min = x.Intv.lo
    midcc,cc_id = mid3(x.cc,x.cv,eps_max)
    midcv,cv_id = mid3(x.cc,x.cv,eps_min)
    cc,dcc = cc_powodd(midcc,x.Intv.lo,x.Intv.hi,c)
    cv,dcv = cv_powodd(midcv,x.Intv.lo,x.Intv.hi,c)
    gcc1,gdcc1 = cc_powodd(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_powodd(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_powodd(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_powodd(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
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
  return MC{N}(cv, cc,  intv, cv_grad, cc_grad, x.cnst)
end

function pos_even(x::MC{N},c::Integer) where {N}
  xL = x.Intv.lo
  xU = x.Intv.hi
  intv = x.Intv^c
  xLc = intv.lo
  xUc = intv.hi
  if (x.Intv.hi<0.0)
    eps_min = x.Intv.hi
    eps_max = x.Intv.lo
  elseif (x.Intv.lo>0.0)
    eps_min = x.Intv.lo
    eps_max = x.Intv.hi
  else
    eps_min = 0.0
    eps_max = (abs(x.Intv.lo)>=abs(x.Intv.hi)) ? x.Intv.lo : x.Intv.hi
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
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, intv, cv_grad, cc_grad, x.cnst)
end

# TO DO: CUT CORRECTED
function neg_powneg_odd(x::MC{N},c::Integer) where {N}
  xL = x.Intv.lo
  xU = x.Intv.hi
  INTV = x.Intv^c
  xLc = INTV.lo
  xUc = INTV.hi
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
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
		xt = MC{N}(cv, cc, x.Intv^c, cv_grad, cc_grad, x.cnst)
		return affine_intv_contract(xt)
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
      cc_grad = zero(x.cc_grad)
    end

    if (xU == xL)
      cv = xLc
      cv_grad = zero(x.cv_grad)
      cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
			xt = MC{N}(cv, cc, x.Intv^c, cv_grad, cc_grad, x.cnst)
			return affine_intv_contract(xt)
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
        cv_grad = zero(x.cv_grad)
      end
      cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
			xt = MC{N}(cv, cc, x.Intv^c, cv_grad, cc_grad, x.cnst)
			return affine_intv_contract(xt)
    end
  end
end

# TO DO: lATER, CUT CORRECTED
function neg_powneg_even(x::MC{N},c::Integer) where {N}
  xL = x.Intv.lo
  xU = x.Intv.hi
  intv = x.Intv^c
  xLc = intv.lo
  xUc = intv.hi
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
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, intv, cv_grad, cc_grad, x.cnst)
end

# TO DO: CHECK CUT CORRECTED
function neg_powpos(x::MC{N},c::Integer) where {N}
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = xL^c
  xUc = xU^c
  eps_min = x.Intv.hi
  eps_max = x.Intv.lo
  if (MC_param.mu >= 1)
    # midcv,cv_id = mid3(x.cc,x.cv,eps_min)
    if (eps_min < x.cv)
      cv = x.cv^c
    elseif (eps_min>x.cc)
      cv = x.cc^c
    else
      cv = eps_min^c
    end

    # midcc,cc_id = mid3(x.cc,x.cv,eps_max)
    m = (xUc-xLc)/(xU-xL)
    if (eps_max < x.cv)
      cc = xLc + m*(x.cv-xL)
    elseif (eps_max > x.cc)
      cc = xLc + m*(x.cc-xL)
    else
      cc = xLc + m*(eps_max-xL)
    end
    cv_grad =(min(0.0,c*x.cc^(c-1)))*x.cc_grad
    cc_grad = m*x.cv_grad
    xt = MC{N}(cv, cc, x.Intv^c, cv_grad, cc_grad, x.cnst)
    return affine_intv_contract(xt)
  else
    xcvc = x.cv^c
    xccc = x.cc^c
    # cv calculation
    if (xU < x.cv)
      cv = xcvc
      cv_grad = (c*x.cv^(c-1))*x.cv_grad
    elseif (xU > x.cc)
      cv = xccc
      cv_grad = (c*x.cc^(c-1))*x.cc_grad
    else
      cv = xU^c
      cv_grad = zero(x.cv_grad)
    end
    # cc calculation
    if (xU == xL)
      cc = xU
      cc_grad = zero(x.cc_grad)
      cv,cc,cv_grad,cc_grad = cut(xUc,xLc,cv,cc,cv_grad,cc_grad)
			xt = MC{N}(cv, cc, x.Intv^c, cv_grad, cc_grad, x.cnst)
      return affine_intv_contract(xt)
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
        cc_grad = zero(x.cc_grad)
      end
      cv,cc,cv_grad,cc_grad = cut(xUc,xLc,cv,cc,cv_grad,cc_grad)
			xt = MC{N}(cv, cc, x.Intv^c, cv_grad, cc_grad, x.cnst)
      return affine_intv_contract(xt)
    end
  end
end

# DONE EXCEPT FOR SUBFUNCTIONS (neg_powneg_even,pos_odd,pos_even,sqr)
function pow(x::MC{N},c::Q) where {N,Q<:Integer}
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
    if (x.Intv.hi < 0.0)
      if (isodd(c))
        return neg_powneg_odd(x,c)
      else
        return neg_powneg_even(x,c) # PRIORITY 2
      end
    elseif (x.Intv.lo > 0.0)
      return neg_powpos(x,c) # MOSTLY DONE
    else
      error("Function unbounded on domain")
    end
  end
end
########### power of a generalized McCormick object raised to c
function pow(x::MC{N},c::Float64) where N
	return x^c
end


# convex/concave relaxation of integer powers of 1/x for negative reals
@inline function cv_flt_pow_1(x::T,xL::T,xU::T,n::Float64) where {T<:AbstractFloat}
    return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
end
@inline function cc_flt_pow_1(x::T,xL::T,xU::T,n::Float64) where {T<:AbstractFloat}
    return x^n,n*x^(n-1)
end
function flt_pow_1(x::MC{N},c::Float64) where N
  intv = pow(x.Intv,c)
  xLc = intv.lo
  xUc = intv.hi
  if (MC_param.mu >= 1)
    eps_max = x.Intv.hi
    eps_min = x.Intv.lo
    midcc,cc_id = mid3(x.cc,x.cv,eps_max)
    midcv,cv_id = mid3(x.cc,x.cv,eps_min)
    cc,dcc = cc_flt_pow_1(midcc,x.Intv.lo,x.Intv.hi,c)
    cv,dcv = cv_flt_pow_1(midcv,x.Intv.lo,x.Intv.hi,c)
    gcc1,gdcc1 = cc_flt_pow_1(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_flt_pow_1(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_flt_pow_1(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_flt_pow_1(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    eps_max = x.Intv.hi
    eps_min = x.Intv.lo
    midcc,cc_id = mid3(x.cc,x.cv,eps_max)
    midcv,cv_id = mid3(x.cc,x.cv,eps_min)
    cc,dcc = cc_flt_pow_1(midcc,x.Intv.lo,x.Intv.hi,c)
    cv,dcv = cv_flt_pow_1(midcv,x.Intv.lo,x.Intv.hi,c)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, intv, cv_grad, cc_grad, x.cnst)
end
function  (^)(x::MC{N},c::Q) where {N,Q<:Integer}
  pow(x,c)
end
function  (^)(x::MC{N},c::Float64) where N
  if (isinteger(c))
    return pow(x,Int(c))
  elseif ((x.Intv.lo >= 0) && (0.0 < c < 1.0))
    return flt_pow_1(x,c)
  else
    return exp(c*log(x))
  end
end

(^)(x::MC{N}, c::MC{N}) where N = exp(c*log(x))
(^)(x::MC{N}, c::Float32) where N = x^Float64(c)
(^)(x::MC{N}, c::Float16) where N = x^Float64(c)

# Define powers to MC of floating point number
function pow(b::Float64, x::MC{N}) where N
	(b <= 0.0) && error("Relaxations of a^x where a<=0 not currently defined in library.
		                   Functions of this type may prevent convergences in global
			                 optimization algorithm as they may be discontinuous.")
	exp(x*log(b))
end
^(b::Float64, x::MC{N}) where N = pow(b::Float64, x::MC{N})

########### Defines inverse
function inv(x::MC{N}) where N
  if (x.Intv.hi < 0.0)
    return neg_powneg_odd(x,-1)
  elseif (x.Intv.lo > 0.0)
    return neg_powpos(x,-1)
  else
    error("Function unbounded on domain: $(x.Intv)")
  end
end
