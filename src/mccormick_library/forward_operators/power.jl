# defines square operator
sqr(x::Float64) = x*x
function cv_sqr_NS(x::Float64, xL::Float64, xU::Float64)
	x^2
end
dcv_sqr_NS(x::Float64, xL::Float64, xU::Float64) = 2.0*x
function cc_sqr(x::Float64, xL::Float64, xU::Float64)
	if (xU > xL)
		cc = xL^2 + (xL + xU)*(x - xL)
	else
		cc = xU^2
	end
	return cc
end
dcc_sqr(x::Float64, xL::Float64, xU::Float64) = (xU > xL) ? (xL + xU) : 0.0
function cv_sqr(x::Float64, xL::Float64, xU::Float64)
    (0.0 <= xL || xU <= 0.0) && return x^2
	((xL < 0.0) && (0.0 <= xU) && (0.0 <= x)) && return (x^3)/xU
	return (x^3)/xL
end
function dcv_sqr(x::Float64, xL::Float64, xU::Float64)
    (0.0 <= xL || xU <= 0.0) && return 2.0*x
	((xL < 0.0) && (0.0 <= xU) && (0.0 <= x)) && (3.0*x^2)/xU
	return (3.0*x^2)/xL
end
function sqr_kernel(x::MC{N,NS}, y::Interval{Float64}) where N
    eps_min = y.lo
    eps_max = y.hi
	midcc, cc_id = mid3(x.cc, x.cv, eps_max)
	midcv, cv_id = mid3(x.cc, x.cv, eps_min)
	cc = cc_sqr(midcc, x.Intv.lo, x.Intv.hi)
	dcc = dcc_sqr(midcc, x.Intv.lo, x.Intv.hi)
	cv = cv_sqr_NS(midcv, x.Intv.lo, x.Intv.hi)
	dcv = dcv_sqr_NS(midcv, x.Intv.lo, x.Intv.hi)
	cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
	cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
	#cv, cc, cv_grad, cc_grad = cut(y.lo, y.hi, cv, cc, cv_grad, cc_grad)
	return MC{N,NS}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end
function sqr_kernel(x::MC{N,Diff}, y::Interval{Float64}) where N
	if (x.Intv.hi < 0.0)
      eps_min = x.Intv.hi
      eps_max = x.Intv.lo
    elseif (x.Intv.lo > 0.0)
      eps_min = x.Intv.lo
      eps_max = x.Intv.hi
    else
      eps_min = 0.0
      eps_max = (abs(x.Intv.lo) >= abs(x.Intv.hi)) ? x.Intv.lo : x.Intv.hi
    end
	midcc, cc_id = mid3(x.cc, x.cv, eps_max)
	midcv, cv_id = mid3(x.cc, x.cv, eps_min)
	cc = cc_sqr(midcc, x.Intv.lo, x.Intv.hi)
	dcc = dcc_sqr(midcc, x.Intv.lo, x.Intv.hi)
    cv = cv_sqr(midcv, x.Intv.lo, x.Intv.hi)
    dcv = dcv_sqr(midcv, x.Intv.lo, x.Intv.hi)
    gdcc1 = dcc_sqr(x.cv, x.Intv.lo, x.Intv.hi)
    gdcv1 = dcv_sqr(x.cv, x.Intv.lo, x.Intv.hi)
    gdcc2 = dcc_sqr(x.cc, x.Intv.lo, x.Intv.hi)
    gdcv2 = dcv_sqr(x.cc, x.Intv.lo, x.Intv.hi)
    cv_grad = max(0.0, gdcv1)*x.cv_grad + min(0.0, gdcv2)*x.cc_grad
    cc_grad = min(0.0, gdcc1)*x.cv_grad + max(0.0, gdcc2)*x.cc_grad
    return MC{N,Diff}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end
sqr(x::MC) = sqr_kernel(x, (x.Intv)^2)

# convex/concave relaxation (Khan 3.1-3.2) of integer powers of 1/x for positive reals
pow_deriv(x::Float64, n::Z) where {Z <: Integer} = n*x^(n-1)
function cv_npp_or_pow4(x::Float64, xL::Float64, xU::Float64, n::Z) where {Z <: Integer}
	x^n, n*x^(n-1)
end
function cc_npp_or_pow4(x::Float64, xL::Float64, xU::Float64, n::Z) where {Z <: Integer}
	dline_seg(^, pow_deriv, x, xL, xU, n)
end
# convex/concave relaxation of integer powers of 1/x for negative reals
function cv_negpowneg(x::Float64, xL::Float64, xU::Float64, n::Z) where {Z <: Integer}
  isodd(n) && (return dline_seg(^, pow_deriv, x, xL, xU, n))
  return x^n, n*x^(n-1)
end
function cc_negpowneg(x::Float64, xL::Float64, xU::Float64, n::Z) where {Z <: Integer}
  isodd(n) && (return x^n,n*x^(n-1))
  return dline_seg(^, pow_deriv, x, xL, xU, n)
end
# convex/concave relaxation of odd powers
function cv_powodd(x::Float64, xL::Float64, xU::Float64, n::Z) where {Z <: Integer}
    (xU <= 0.0) && (return dline_seg(^, pow_deriv, x, xL, xU, n))
    (0.0 <= xL) && (return x^n, n*x^(n - 1))
    val = (xL^n)*(xU - x)/(xU - xL) + (max(0.0, x))^n
    dval = -(xL^n)/(xU - xL) + n*(max(0.0, x))^(n-1)
    return val, dval
end
function cc_powodd(x::Float64, xL::Float64, xU::Float64, n::Z) where {Z <: Integer}
    (xU <= 0.0) && (return x^n, n*x^(n - 1))
    (0.0 <= xL) && (return dline_seg(^, pow_deriv, x, xL, xU, n))
    val = (xU^n)*(x - xL)/(xU - xL) + (min(0.0, x))^n
    dval = (xU^n)/(xU - xL) + n*(min(0.0, x))^(n-1)
    return val, dval
end

function npp_or_pow4(x::MC{N,NS}, c::Z, y::Interval{Float64}) where {N, Z<:Integer}
  if (x.Intv.hi < 0.0)
    eps_min = x.Intv.hi
    eps_max = x.Intv.lo
  elseif (x.Intv.lo > 0.0)
    eps_min = x.Intv.lo
    eps_max = x.Intv.hi
  else
    eps_min = 0.0
    eps_max = (abs(x.Intv.lo) >= abs(x.Intv.hi)) ? x.Intv.lo : x.Intv.hi
  end
  midcc, cc_id = mid3(x.cc, x.cv, eps_max)
  midcv, cv_id = mid3(x.cc, x.cv, eps_min)
  cc, dcc = cc_npp_or_pow4(midcc, x.Intv.lo, x.Intv.hi, c)
  cv, dcv = cv_npp_or_pow4(midcv, x.Intv.lo, x.Intv.hi, c)
  cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
  cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
	cv, cc, cv_grad, cc_grad = cut(y.lo, y.hi, cv, cc, cv_grad, cc_grad)
  return MC{N,NS}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end
function npp_or_pow4(x::MC{N,Diff}, c::Z, y::Interval{Float64}) where {N, Z<:Integer}
  if (x.Intv.hi < 0.0)
    eps_min = x.Intv.hi
    eps_max = x.Intv.lo
  elseif (x.Intv.lo > 0.0)
    eps_min = x.Intv.lo
    eps_max = x.Intv.hi
  else
    eps_min = 0.0
    eps_max = (abs(x.Intv.lo) >= abs(x.Intv.hi)) ? x.Intv.lo : x.Intv.hi
  end
  midcc, cc_id = mid3(x.cc, x.cv, eps_max)
  midcv, cv_id = mid3(x.cc, x.cv, eps_min)
  cc, dcc = cc_npp_or_pow4(midcc, x.Intv.lo, x.Intv.hi, c)
  cv, dcv = cv_npp_or_pow4(midcv, x.Intv.lo, x.Intv.hi, c)
  gcc1, gdcc1 = cc_npp_or_pow4(x.cv, x.Intv.lo, x.Intv.hi, c)
  gcv1, gdcv1 = cv_npp_or_pow4(x.cv, x.Intv.lo, x.Intv.hi, c)
  gcc2, gdcc2 = cc_npp_or_pow4(x.cc, x.Intv.lo, x.Intv.hi, c)
  gcv2, gdcv2 = cv_npp_or_pow4(x.cc, x.Intv.lo, x.Intv.hi, c)
  cv_grad = max(0.0, gdcv1)*x.cv_grad + min(0.0, gdcv2)*x.cc_grad
  cc_grad = min(0.0, gdcc1)*x.cv_grad + max(0.0, gdcc2)*x.cc_grad
  return MC{N,Diff}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end

function pos_odd(x::MC{N,NS}, c::Z, y::Interval{Float64}) where {N, Z<:Integer}
    midcc, cc_id = mid3(x.cc, x.cv, x.Intv.hi)
    midcv, cv_id = mid3(x.cc, x.cv, x.Intv.lo)
    cc, dcc = cc_powodd(midcc, x.Intv.lo, x.Intv.hi, c)
    cv, dcv = cv_powodd(midcv, x.Intv.lo, x.Intv.hi, c)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
		cv, cc, cv_grad, cc_grad = cut(y.lo, y.hi, cv, cc, cv_grad, cc_grad)
    return MC{N,NS}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end
function pos_odd(x::MC{N,Diff}, c::Z, y::Interval{Float64}) where {N, Z<:Integer}
    midcc, cc_id = mid3(x.cc, x.cv, x.Intv.hi)
    midcv, cv_id = mid3(x.cc, x.cv, x.Intv.lo)
    cc, dcc = cc_powodd(midcc, x.Intv.lo, x.Intv.hi, c)
    cv, dcv = cv_powodd(midcv, x.Intv.lo, x.Intv.hi, c)
    gcc1, gdcc1 = cc_powodd(x.cv, x.Intv.lo, x.Intv.hi, c)
    gcv1, gdcv1 = cv_powodd(x.cv, x.Intv.lo, x.Intv.hi, c)
    gcc2, gdcc2 = cc_powodd(x.cc, x.Intv.lo, x.Intv.hi, c)
    gcv2, gdcv2 = cv_powodd(x.cc, x.Intv.lo, x.Intv.hi, c)
    cv_grad = max(0.0, gdcv1)*x.cv_grad + min(0.0, gdcv2)*x.cc_grad
    cc_grad = min(0.0, gdcc1)*x.cv_grad + max(0.0, gdcc2)*x.cc_grad
    return MC{N,Diff}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end

function neg_powneg_odd(x::MC{N,NS}, c::Z, y::Interval{Float64}) where {N, Z<:Integer}
  xL = x.Intv.lo
  xU = x.Intv.hi
  xUc = y.hi
  eps_max = xU
  eps_min = xL
  if (xL < x.cv)
    cc = x.cv^c
    cc_grad = (c*x.cv^(c-1))*x.cv_grad
  elseif (xL > x.cc)
    cc = x.cc^c
    cc_grad = (c*x.cc^(c-1))*x.cc_grad
  else
    cc = xL^c
    cc_grad = zeros(SVector{N,Float64})
  end
  if (xU == xL)
    cv = xLc
    cv_grad = zeros(SVector{N,Float64})
  else
    dcv = (xU^c - xL^c)/(xU - xL) # function decreasing
    if (xU < x.cv)
      cv = xUc + dcv*(x.cv - xL)
      cv_grad = dcv*x.cv_grad
    elseif (xU > x.cc)
      cv = xUc + dcv*(x.cc - xL)
      cv_grad = dcv*x.cc_grad
    else
      cv = xUc + dcv*(xU - xL)
      cv_grad = zeros(SVector{N,Float64})
    end
  end
	cv, cc, cv_grad, cc_grad = cut(y.lo,xUc, cv, cc, cv_grad, cc_grad)
	return MC{N,NS}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end
function neg_powneg_odd(x::MC{N,Diff}, c::Z, y::Interval{Float64}) where {N, Z<:Integer}
  xL = x.Intv.lo
  xU = x.Intv.hi
  xUc = y.hi
  eps_max = xU
  eps_min = xL
  midcc, cc_id = mid3(x.cc, x.cv, eps_max)
  midcv, cv_id = mid3(x.cc, x.cv, eps_min)
  cc, dcc = cc_negpowneg(midcc, xL, xU, c)
  cv, dcv = cv_negpowneg(midcv, xL, xU, c)
  gcc1, gdcc1 = cc_negpowneg(x.cv, xL, xU, c)
  gcv1, gdcv1 = cv_negpowneg(x.cv, xL, xU, c)
  gcc2, gdcc2 = cc_negpowneg(x.cc, xL, xU, c)
  gcv2, gdcv2 = cv_negpowneg(x.cc, xL, xU, c)
  cv_grad = max(0.0, gdcv1)*x.cv_grad + min(0.0, gdcv2)*x.cc_grad
  cc_grad = min(0.0, gdcc1)*x.cv_grad + max(0.0, gdcc2)*x.cc_grad
	return MC{N,Diff}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end

function neg_powneg_even(x::MC{N,NS}, c::Z, y::Interval{Float64}) where {N, Z<:Integer}
  midcc, cc_id = mid3(x.cc, x.cv, x.Intv.hi)
  midcv, cv_id = mid3(x.cc, x.cv, x.Intv.lo)
  cc, dcc = cc_negpowneg(midcc, x.Intv.lo, x.Intv.hi, c)
  cv, dcv = cv_negpowneg(midcv, x.Intv.lo, x.Intv.hi, c)
  cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
  cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
	cv, cc, cv_grad, cc_grad = cut(y.lo, y.hi, cv, cc, cv_grad, cc_grad)
  return MC{N,NS}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end
function neg_powneg_even(x::MC{N,Diff}, c::Z, y::Interval{Float64}) where {N, Z<:Integer}
  midcc, cc_id = mid3(x.cc, x.cv, x.Intv.hi)
  midcv, cv_id = mid3(x.cc, x.cv, x.Intv.lo)
  cc, dcc = cc_negpowneg(midcc, x.Intv.lo, x.Intv.hi, c)
  cv, dcv = cv_negpowneg(midcv, x.Intv.lo, x.Intv.hi, c)
  gcc1, gdcc1 = cc_negpowneg(x.cv, x.Intv.lo, x.Intv.hi, c)
  gcv1, gdcv1 = cv_negpowneg(x.cv, x.Intv.lo, x.Intv.hi, c)
  gcc2, gdcc2 = cc_negpowneg(x.cc, x.Intv.lo, x.Intv.hi, c)
  gcv2, gdcv2 = cv_negpowneg(x.cc, x.Intv.lo, x.Intv.hi, c)
  cv_grad = max(0.0, gdcv1)*x.cv_grad + min(0.0, gdcv2)*x.cc_grad
  cc_grad = min(0.0, gdcc1)*x.cv_grad + max(0.0, gdcc2)*x.cc_grad
  return MC{N,Diff}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end

function pow_kernel(x::MC, c::Z, y::Interval{Float64}) where {Z<:Integer}
    if (c == 0)
		    z = one(x)
	  elseif (c == 1)
		    z = x
	  elseif (c > 0)
        if (c == 2)
			#println("x = $x")
			#println("y = $y")
			      z = sqr_kernel(x, y)
				elseif isodd(c)
						z = pos_odd(x, c, y)
				else
						z =	npp_or_pow4(x, c, y)
				end
    else
        if (x.Intv.hi < 0.0)
        		if isodd(c)
						  	z = neg_powneg_odd(x, c, y)
						else
        	  		z = neg_powneg_even(x, c, y)
						end
        elseif (x.Intv.lo > 0.0)
					  z = npp_or_pow4(x, c, y)
				else
					error("Envelope not defined.")
				end
    end
	  return z
end
function pow(x::MC, c::Z) where {Z<:Integer}
	if (x.Intv.lo <= 0.0 <= x.Intv.hi) && (c < 0)
		error("Function unbounded on this domain")
	end
	return pow_kernel(x, c, x.Intv^c)
end
(^)(x::MC, c::Z) where {Z <: Integer} = pow(x,c)

# Power of MC to float
cv_flt_pow_1(x::Float64, xL::Float64, xU::Float64, n::Float64) = dline_seg(^, pow_deriv, x, xL, xU, n)
cc_flt_pow_1(x::Float64, xL::Float64, xU::Float64, n::Float64) = x^n, n*x^(n-1)
function flt_pow_1(x::MC{N,NS}, c::Float64, y::Interval{Float64}) where N
	midcc, cc_id = mid3(x.cc, x.cv, x.Intv.hi)
	midcv, cv_id = mid3(x.cc, x.cv, x.Intv.lo)
	cc, dcc = cc_flt_pow_1(midcc, x.Intv.lo, x.Intv.hi, c)
	cv, dcv = cv_flt_pow_1(midcv, x.Intv.lo, x.Intv.hi, c)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv, cc, cv_grad, cc_grad = cut(y.lo, y.hi, cv, cc, cv_grad, cc_grad)
    return MC{N,NS}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end
function flt_pow_1(x::MC{N,Diff}, c::Float64, y::Interval{Float64}) where N
	midcc, cc_id = mid3(x.cc, x.cv, x.Intv.hi)
	midcv, cv_id = mid3(x.cc, x.cv, x.Intv.lo)
	cc, dcc = cc_flt_pow_1(midcc, x.Intv.lo, x.Intv.hi, c)
	cv, dcv = cv_flt_pow_1(midcv, x.Intv.lo, x.Intv.hi, c)
    gcc1, gdcc1 = cc_flt_pow_1(x.cv ,x.Intv.lo, x.Intv.hi, c)
    gcv1, gdcv1 = cv_flt_pow_1(x.cv, x.Intv.lo, x.Intv.hi, c)
    gcc2, gdcc2 = cc_flt_pow_1(x.cc, x.Intv.lo, x.Intv.hi, c)
    gcv2, gdcv2 = cv_flt_pow_1(x.cc, x.Intv.lo, x.Intv.hi, c)
    cv_grad = max(0.0, gdcv1)*x.cv_grad + min(0.0, gdcv2)*x.cc_grad
    cc_grad = min(0.0, gdcc1)*x.cv_grad + max(0.0, gdcc2)*x.cc_grad
    return MC{N,Diff}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end

function (^)(x::MC{N,NS}, c::Float64, y::Interval{Float64}) where N
    isinteger(c) && (return pow_kernel(x, Int(c), y))
    ((x.Intv.lo >= 0) && (0.0 < c < 1.0)) && (return flt_pow_1(x, c, y))
	z = exp(c*log(x))
    return MC{N,NS}(z.cv, z.cc, y, z.cv_grad, z.cc_grad, x.cnst)
end
function (^)(x::MC{N,Diff}, c::Float64, y::Interval{Float64}) where N
    isinteger(c) && (return pow_kernel(x, Int(c), y))
    ((x.Intv.lo >= 0) && (0.0 < c < 1.0)) && (return flt_pow_1(x, c, y))
	z = exp(c*log(x))
    return MC{N,Diff}(z.cv, z.cc, y, z.cv_grad, z.cc_grad, x.cnst)
end

(^)(x::MC, c::Float32, y::Interval{Float64}) = (^)(x, Float32(c), y)
(^)(x::MC, c::Float16, y::Interval{Float64}) = (^)(x, Float16(c), y)
(^)(x::MC, c::Float64) = (^)(x, c, x.Intv^c)
(^)(x::MC, c::Float32) = x^Float64(c) # DONE
(^)(x::MC, c::Float16) = x^Float64(c) # DONE
(^)(x::MC, c::MC) = exp(c*log(x)) # DONE (no kernel)
pow(x::MC, c::F) where {F <: AbstractFloat} = x^c

# Define powers to MC of floating point number
function pow(b::Float64, x::MC) # DONE (no kernel)
	(b <= 0.0) && error("Relaxations of a^x where a<=0 not currently defined in library.
		                   Functions of this type may prevent convergences in global
			                 optimization algorithm as they may be discontinuous.")
	exp(x*log(b))
end
^(b::Float64, x::MC) = pow(b, x) # DONE (no kernel)

########### Defines inverse
function cc_inv1(x::Float64, xL::Float64, xU::Float64)
	t = (xL*xU)
	cc = (xU + xL - x)/t
	dcc = -1.0/t
	return cc, dcc
end
function cv_inv1(x::Float64, xL::Float64, xU::Float64)
	cv = 1.0/x
	dcv = -1.0/(x*x)
	return cv, dcv
end
function inv1(x::MC{N,NS}, y::Interval{Float64}) where N
  midcc, cc_id = mid3(x.cc, x.cv, x.Intv.lo)
  midcv, cv_id = mid3(x.cc, x.cv, x.Intv.hi)
  cc, dcc = cc_inv1(midcc, x.Intv.lo, x.Intv.hi)
  cv, dcv = cv_inv1(midcv, x.Intv.lo, x.Intv.hi)
  cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
  cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  cv, cc, cv_grad, cc_grad = cut(y.lo, y.hi, cv, cc, cv_grad, cc_grad)
  return MC{N,NS}(cv, cc, y, cv_grad, cc_grad, x.cnst)
end
function inv_kernel(x::MC{N,T}, y::Interval{Float64}) where {N,T<:RelaxTag}
	if (x.Intv.lo <= 0.0 <= x.Intv.hi)
		error("Function unbounded on domain: $(x.Intv)")
	end
	if (x.Intv.hi < 0.0)
		x = neg_powneg_odd(x, -1, y)
  	elseif (x.Intv.lo > 0.0)
		x = inv1(x, y)
	end
	return x
end
inv(x::MC{N,T}) where {N,T<:RelaxTag} = inv_kernel(x, (x.Intv)^(-1))
