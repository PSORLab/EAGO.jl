
function interval_MC(x::MC{S,T}) where {S,T<:RelaxTag}
	MC{S,T}(x.Intv)
end

@inline isnan(x::MC) = isnan(x.cc) || isnan(x.cv)

########### Defines differentiable step relaxations
@inline function cv_step(x::Float64, xL::Float64, xU::Float64)
	  (xU <= 0.0) && (return 0.0, 0.0)
	  (xL >= 0.0) && (return 1.0, 0.0)
	  (x >= 0.0) ? ((x/xU)^2, 2.0*x/xU^2) : (0.0, 0.0)
end
@inline function cc_step(x::Float64, xL::Float64, xU::Float64)
	  (xU <= 0.0) && (return 0.0, 0.0)
	  (xL >= 0.0) && (return 1.0, 0.0)
	  (x >= 0.0) ? (1.0, 0.0) : (1.0-(x/xL)^2, -2.0*x/xL^2)
end
@inline function cv_step_NS(x::Float64, xL::Float64, xU::Float64)
	  (xU <= 0.0) && (return 0.0, 0.0)
	  (xL >= 0.0) && (return 1.0, 0.0)
	  (x > 0.0) ? (x/xU, 1.0/xU) : (0.0, 0.0)
end
@inline function cc_step_NS(x::Float64, xL::Float64, xU::Float64)
	  (xU <= 0.0) && (return 0.0, 0.0)
	  (xL >= 0.0) && (return 1.0, 0.0)
	  (x >= 0.0) ? (1.0, 0.0) : ((1.0 - (x/xL)), (-x/xL))
end
@inline function step_kernel(x::MC{N, NS}, z::Interval{Float64}) where N
	xL = x.Intv.lo
	xU = x.Intv.hi
	xLc = z.lo
	xUc = z.hi
	midcc, cc_id = mid3(x.cc, x.cv, xU)
	midcv, cv_id = mid3(x.cc, x.cv, xL)
	cc, dcc = cc_step_NS(midcc,x.Intv.lo,x.Intv.hi)
	cv, dcv = cv_step_NS(midcv,x.Intv.lo,x.Intv.hi)
	cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
	cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
	cv, cc, cv_grad, cc_grad = cut(xLc, xUc, cv, cc, cv_grad, cc_grad)
	return MC{N,NS}(cv, cc, z, cv_grad, cc_grad, x.cnst)
end
@inline function step_kernel(x::MC{N, Diff}, z::Interval{Float64}) where N
	 xL = x.Intv.lo
	 xU = x.Intv.hi
	 xLc = z.lo
	 xUc = z.hi
	 midcc, cc_id = mid3(x.cc, x.cv, xU)
	 midcv, cv_id = mid3(x.cc, x.cv, xL)
	 cc, dcc = cc_step(midcc, x.Intv.lo, x.Intv.hi)
	 cv, dcv = cv_step(midcv, x.Intv.lo, x.Intv.hi)
	 cc1, gdcc1 = cc_step(x.cv, x.Intv.lo, x.Intv.hi)
	 cv1, gdcv1 = cv_step(x.cv, x.Intv.lo, x.Intv.hi)
	 cc2, gdcc2 = cc_step(x.cc, x.Intv.lo, x.Intv.hi)
	 cv2, gdcv2 = cv_step(x.cc, x.Intv.lo, x.Intv.hi)
	 cv_grad = max(0.0, gdcv1)*x.cv_grad + min(0.0, gdcv2)*x.cc_grad
	 cc_grad = min(0.0, gdcc1)*x.cv_grad + max(0.0, gdcc2)*x.cc_grad
	 return MC{N,Diff}(cv, cc, z, cv_grad, cc_grad, x.cnst)
end
@inline step(x::MC) = step_kernel(x, step(x.Intv))

########### Defines sign
@inline function sign_kernel(x::MC{N, T}, z::Interval{Float64}) where {N, T<:RelaxTag}
	zMC = -step(-x) + step(x)
	return MC{N,T}(zMC.cv, zMC.cc, z, zMC.cv_grad, zMC.cc_grad, zMC.cnst)
end
@inline sign(x::MC) = sign_kernel(x, sign(x.Intv))

@inline function cv_abs(x::Float64, xL::Float64, xU::Float64)
	 (xL >= 0.0) && (return x, 1.0)
	 (xU <= 0.0) && (return -x, -1.0)
	 if (x >= 0.0)
	 	return xU*(x/xU)^(MC_DIFF_MU+1), (MC_DIFF_MU+1)*(x/xU)^MC_DIFF_MU
	 else
		return -xL*(x/xL)^(MC_DIFF_MU+1), -(MC_DIFF_MU+1)*(x/xL)^MC_DIFF_MU
	 end
 end
@inline cc_abs(x::Float64,xL::Float64,xU::Float64) = dline_seg(abs, sign, x, xL, xU)
@inline cv_abs_NS(x::Float64,xL::Float64,xU::Float64) = abs(x), sign(x)

@inline function abs_kernel(x::MC{N, NS}, z::Interval{Float64}) where N
    xL = x.Intv.lo
	xU = x.Intv.hi
	xLc = z.lo
	xUc = z.hi
	eps_min, blank = mid3(x.Intv.lo,x.Intv.hi,0.0)
	eps_max = (abs(x.Intv.hi) >= abs(x.Intv.lo)) ? x.Intv.hi : x.Intv.lo
	midcc,cc_id = mid3(x.cc, x.cv, eps_max)
	midcv,cv_id = mid3(x.cc, x.cv, eps_min)
	cc, dcc = cc_abs(midcc, x.Intv.lo, x.Intv.hi)
	cv, dcv = cv_abs_NS(midcv, x.Intv.lo, x.Intv.hi)
	cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
	cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
	cv, cc, cv_grad, cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
	return MC{N,NS}(cv, cc, z, cv_grad, cc_grad, x.cnst)
end
@inline function abs_kernel(x::MC{N, Diff}, z::Interval{Float64}) where N
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = z.lo
  xUc = z.hi
  eps_min, blank = mid3(x.Intv.lo,x.Intv.hi,0.0)
  eps_max = (abs(x.Intv.hi) >= abs(x.Intv.lo)) ? x.Intv.hi : x.Intv.lo
  midcc,cc_id = mid3(x.cc, x.cv, eps_max)
  midcv,cv_id = mid3(x.cc, x.cv, eps_min)
  cc, dcc = cc_abs(midcc, x.Intv.lo, x.Intv.hi)
  cv, dcv = cv_abs(midcv, x.Intv.lo, x.Intv.hi)
  cc1, gdcc1 = cc_abs(x.cv,x.Intv.lo,x.Intv.hi)
  cv1, gdcv1 = cv_abs(x.cv,x.Intv.lo,x.Intv.hi)
  cc2, gdcc2 = cc_abs(x.cc,x.Intv.lo,x.Intv.hi)
  cv2, gdcv2 = cv_abs(x.cc,x.Intv.lo,x.Intv.hi)
  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  return MC{N, Diff}(cv, cc, z, cv_grad, cc_grad, x.cnst)
end
@inline abs(x::MC) = abs_kernel(x, abs(x.Intv))

@inline function intersect(x::MC{N, NS}, y::MC{N, NS}) where N
    if (x.cv >= y.cv)
      cv = x.cv
      cv_grad = x.cv_grad
    else
      cv = y.cv
      cv_grad = y.cv_grad
    end
    if (x.cc <= y.cc)
      cc = x.cc
      cc_grad = x.cc_grad
    else
      cc = y.cc
      cc_grad = y.cc_grad
    end
	# ADD NAN CHECK HERE
    return MC{N, NS}(cv, cc, intersect(x.Intv,y.Intv), cv_grad, cc_grad, (x.cnst && y.cnst))
end
@inline function intersect(x::MC{N, Diff}, y::MC{N, Diff}) where N
    max_MC = x - max(x - y, 0.0)   # used for convex
    min_MC = y - max(y - x, 0.0)   # used for concave
    return MC{N,Diff}(max_MC.cv, min_MC.cc, intersect(x.Intv,y.Intv), max_MC.cv_grad, min_MC.cc_grad, (x.cnst && y.cnst))
end

@inline function intersect(x::MC{N, NS}, y::Interval{Float64}) where N

	if (x.cv >= y.lo)
  		cv = x.cv
  		cv_grad = x.cv_grad
	else
  		cv = y.lo
  		cv_grad = zero(SVector{N,Float64})
	end
	if (x.cc <= y.hi)
  		cc = x.cc
  		cc_grad = x.cc_grad
	else
  		cc = y.hi
  		cc_grad = zero(SVector{N,Float64})
	end
	if (cc < y.lo)
		cv = NaN
		cc = NaN
	end
	if (y.hi < cv)
		cv = NaN
		cc = NaN
	end
    return MC{N, NS}(cv, cc, intersect(x.Intv,y), cv_grad, cc_grad, (x.cnst && y.cnst))
end
@inline function intersect(x::MC{N, Diff}, y::Interval{Float64}) where N
     max_MC = x - max(x - y, 0.0)   # used for convex
     min_MC = y - max(y - x, 0.0)   # used for concave
     return MC{N, Diff}(max_MC.cv, min_MC.cc, intersect(x.Intv,y),
	                    max_MC.cv_grad, min_MC.cc_grad, (x.cnst && y.cnst))
end

@inline in(x::MC) = in(x.Intv)
@inline isempty(x::MC) = isempty(x.Intv)
