########### Defines differentiable step relaxations
function cv_step(x::Float64,xL::Float64,xU::Float64)
  (xU <= 0.0) && (return 0.0)
  (xL >= 0.0) && (return 1.0)
  (x >= 0.0) ? (x/xU)^2 : 0.0
end

function dcv_step(x::Float64,xL::Float64,xU::Float64)
  (xU <= 0.0) && (return 0.0)
  (xL >= 0.0) && (return 0.0)
  (x >= 0.0) ? 2.0*x/xU^2 : 0.0
end

function cc_step(x::Float64,xL::Float64,xU::Float64)
  (xU <= 0.0) && (return 0.0)
  (xL >= 0.0) && (return 1.0)
  (x >= 0.0) ? 1.0 : 1.0-(x/xL)^2
end

function dcc_step(x::Float64,xL::Float64,xU::Float64)
  (xU <= 0.0) && (return 0.0)
  (xL >= 0.0) && (return 0.0)
  (x >= 0.0) ? 0.0 : -2.0*x/xL^2
end

function cv_step_NS(x::Float64,xL::Float64,xU::Float64)
  (xU <= 0.0) && (return 0.0)
  (xL >= 0.0) && (return 1.0)
  (x > 0.0) ? (x/xU) : 0.0
end

function dcv_step_NS(x::Float64,xL::Float64,xU::Float64)
  (xU <= 0.0) && (return 0.0)
  (xL >= 0.0) && (return 0.0)
  (x >= 0.0) ? 1.0/xU : 0.0
end

function cc_step_NS(x::Float64,xL::Float64,xU::Float64)
  (xU <= 0.0) && (return 0.0)
  (xL >= 0.0) && (return 1.0)
  (x >= 0.0) ? 1.0 : (1.0 - (x/xL))
end

function dcc_step_NS(x::Float64,xL::Float64,xU::Float64)
  (xU <= 0.0) && (return 0.0)
  (xL >= 0.0) && (return 0.0)
  (x >= 0.0) ? 0.0 : (-x/xL)
end

function step(x::MC{N}) where N
  Intv = step(x.Intv)
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = Intv.lo
  xUc = Intv.hi
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  if (MC_param.mu >= 1)
    cc = cc_step(midcc,x.Intv.lo,x.Intv.hi)
    dcc = dcc_step(midcc,x.Intv.lo,x.Intv.hi)
    cv = cv_step(midcv,x.Intv.lo,x.Intv.hi)
    dcv = dcv_step(midcv,x.Intv.lo,x.Intv.hi)
    gdcc1 = dcc_step(x.cv,x.Intv.lo,x.Intv.hi)
	  gdcv1 = dcv_step(x.cv,x.Intv.lo,x.Intv.hi)
	  gdcc2 = dcc_step(x.cc,x.Intv.lo,x.Intv.hi)
	  gdcv2 = dcv_step(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
	  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc = cc_step_NS(midcc,x.Intv.lo,x.Intv.hi)
    dcc = dcc_step_NS(midcc,x.Intv.lo,x.Intv.hi)
    cv = cv_step_NS(midcv,x.Intv.lo,x.Intv.hi)
    dcv = dcv_step_NS(midcv,x.Intv.lo,x.Intv.hi)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

########### Defines sign
sign(x::MC) = -step(-x) + step(x)

function cv_abs(x,xL,xU)
	 (x >= 0.0) ? xU*(x/xU)^(MC_param.mu+1) : -xL*(x/xL)^(MC_param.mu+1)
 end
function dcv_abs(x,xL,xU)
	 (x >= 0.0) ? (MC_param.mu+1)*(x/xU)^MC_param.mu : -(MC_param.mu+1)*(x/xL)^MC_param.mu
 end
function cc_abs(x::Float64,xL::Float64,xU::Float64)
	(x >= 0.0) ? (x/xU)^(MC_param.mu+1) : -(x/xL)^(MC_param.mu+1)
end
function dcc_abs(x::Float64,xL::Float64,xU::Float64)
	(MC_param.mu+1)*(x/xU)^MC_param.mu : -(MC_param.mu+1)*(x/xL)^MC_param.mu
end
cc_abs_NS(x::Float64,xL::Float64,xU::Float64) = line_seg(x,xL,abs(xL),xU,abs(xU))[1]
dcc_abs_NS(x::Float64,xL::Float64,xU::Float64) = dline_seg(x,xL,abs(xL),xU,abs(xU),sign(x))[1]
cv_abs_NS(x::Float64,xL::Float64,xU::Float64) = abs(x)
dcv_abs_NS(x::Float64,xL::Float64,xU::Float64) = sign(x)

function abs(x::MC{N}) where N
  Intv = abs(x.Intv)
  xL = x.Intv.lo
  xU = x.Intv.hi
  xLc = Intv.lo
  xUc = Intv.hi
  eps_min,blank = mid3(x.Intv.lo,x.Intv.hi,0.0)
  eps_max = (abs(x.Intv.hi)>=abs(x.Intv.lo)) ? x.Intv.hi : x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  if (MC_param.mu >= 1)
    cc = cc_abs(midcc,x.Intv.lo,x.Intv.hi)
    dcc = dcc_abs(midcc,x.Intv.lo,x.Intv.hi)
    cv = cv_abs(midcv,x.Intv.lo,x.Intv.hi)
    dcv = dcv_abs(midcv,x.Intv.lo,x.Intv.hi)
    gdcc1 = cc_abs(x.cv,x.Intv.lo,x.Intv.hi)
	gdcv1 = cv_abs(x.cv,x.Intv.lo,x.Intv.hi)
	gdcc2 = cc_abs(x.cc,x.Intv.lo,x.Intv.hi)
	gdcv2 = cv_abs(x.cc,x.Intv.lo,x.Intv.hi)
	cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
	cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc = cc_abs_NS(midcc, x.Intv.lo, x.Intv.hi)
    dcc = dcc_abs_NS(midcc, x.Intv.lo, x.Intv.hi)
    cv = cv_abs_NS(midcv, x.Intv.lo, x.Intv.hi)
    dcv = dcv_abs_NS(midcv, x.Intv.lo, x.Intv.hi)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

function intersect(x::MC{N}, y::MC{N}) where N
  if (MC_param.mu >= 1)
    max_MC = x - max(x - y, 0.0)   # used for convex
    min_MC = y - max(y - x, 0.0)   # used for concave
    return MC{N}(max_MC.cv, min_MC.cc, intersect(x.Intv,y.Intv), max_MC.cv_grad, min_MC.cc_grad, (x.cnst && y.cnst))
  else
    # Convex terms
    if (x.cv >= y.cv)
      cv = x.cv
      cv_grad = x.cv_grad
    else
      cv = y.cv
      cv_grad = y.cv_grad
			if (x.cc < y.cv)
				cc = y.cv
				cc_grad = zero(x.cc_grad)
			end
    end
    # Concave terms
    if (x.cc <= y.cc)
      cc = x.cc
      cc_grad = x.cc_grad
    else
      cc = y.cc
      cc_grad = y.cc_grad
			if (x.cv > y.cc)
				cv = y.cc
				cv_grad = zero(x.cc_grad)
			end
    end
    return MC{N}(cv, cc, intersect(x.Intv,y.Intv), cv_grad, cc_grad, (x.cnst && y.cnst))
  end
end

function intersect(x::MC{N}, y::IntervalType) where N
    if (MC_param.mu >= 1)
        max_MC = x - max(x - y, 0.0)   # used for convex
    		min_MC = y - max(y - x, 0.0)   # used for concave
    		return MC{N}(max_MC.cv, min_MC.cc, intersect(Intv(x),y), max_MC.cv_grad, min_MC.cc_grad, (x.cnst && y.cnst))
  	else
    	#println("ran intersect 1")
		#	println("x: $x")
		#	println("y: $y")
    	# Convex terms
    	if (x.cv >= lo(y))
      	cv = x.cv
      	cv_grad = x.cv_grad
    	else
      	cv = lo(y)
	  		cc = x.cc
      	cv_grad = zero(x.cv_grad)
	  		if (cc < lo(y))
		  		cc = lo(y)
		  		cc_grad = zero(x.cc_grad)
	  		end
    	end
    	# Concave terms
    	if (x.cc <= hi(y))
      	cc = x.cc
      	cc_grad = x.cc_grad
    	else
	  		cv = x.cv
      	cc = hi(y)
      	cc_grad = zero(x.cc_grad)
	  		if (hi(y) < cv)
		  		cv = hi(y)
		  		cv_grad = zero(x.cc_grad)
	  		end
    	end
		#println("ran intersect 1: $(MC{N}(cv, cc, intersect(Intv(x),y), cv_grad, cc_grad, (x.cnst && y.cnst)))")
    return MC{N}(cv, cc, intersect(Intv(x),y), cv_grad, cc_grad, (x.cnst && y.cnst))
  end
end
#=
function hull(x::MC{N}, y::MC{N}) where N
  #return MC{N}(cv, cc, intersect(Intv(a),Intv(b)), cv_grad, cc_grad, (x.cnst && y.cnst))
end
=#
in(x::MC) = in(Intv(x))
isempty(x::MC) = isempty(Intv(x))
