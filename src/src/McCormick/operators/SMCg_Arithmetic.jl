function +(x::SMCg{N,V,T},y::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	return SMCg{N,V,T}(x.cc+y.cc, x.cv+y.cv, x.cc_grad+y.cc_grad, x.cv_grad+y.cv_grad,
						 (x.Intv+y.Intv),(x.cnst && y.cnst),x.IntvBox,x.xref)
end
function -(x::SMCg{N,V,T},c::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	return x + (-c)
end
function -(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	return SMCg{N,V,T}(-x.cv, -x.cc, -x.cv_grad, -x.cc_grad, -x.Intv,
							 x.cnst, x.IntvBox, x.xref)
end
function /(x::SMCg{N,V,T},y::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	pos_orth::Bool = (x.Intv.lo >= 0) && (y.Intv.lo >= 0)
	if (x==y)
		#println("ran me to 1")
		return one(T)
	#=
	elseif (MC_param.multivar_refine) && (~ (MC_param.mu >= 1)) && (pos_orth)
		#println("ran me to 2")
		Intv::V = x.Intv/y.Intv
		xLc::T = Intv.lo
		xUc::T = Intv.hi
		cv1::T,pos1::Int64 = mid3(x.cv,x.cc,x.Intv.lo)
		cv1t::T = (cv1+sqrt(x.intv.lo*x.intv.hi))
				  /(sqrt(x.intv.lo)+sqrt(x.intv.hi))
		cv2t::T,pos2::Int64 = mid3(y.cv,y.cc,y.Intv.hi)
		cv::T = sqr(cv1t)/cv2t
		cv_grad::SVector{N,T} = 2.0*(cv1t/cv2t)/(sqrt(x.Intv.lo)+sqrt(x.Intv.hi))*
								  mid_grad(x.cc_grad, x.cv_grad, pos1) - ((cv1t/cv2t)^2)*
									mid_grad(y.cc_grad, y.cv_grad, pos2)
		cc1::T,pos1 = mid3(x.cv,x.cc,x.Intv.hi)
		cc2::T,pos2 = mid3(y.cv,y.cc,y.Intv.lo)
		gcc1::T = y.Intv.hi*cc1 - x.Intv.lo*cc2 + x.Intv.lo*y.Intv.lo
		gcc2::T = y.Intv.lo*cc1 - x.Intv.hi*cc2 + x.Intv.hi*y.Intv.hi
		if gcc1 <= gcc2
			cc::T = gcc1/(y.Intv.lo*y.Intv.hi)
			cc_grad::SVector{N,T} = one(T)/y.Intv.lo*mid_grad(x.cc_grad, x.cv_grad, pos1)
							- x.Intv.lo/(y.Intv.lo*y.Intv.hi)*mid_grad(y.cc_grad, y.cv_grad, pos2)
		else
			cc = gcc2/(y.Intv.lo*y.Intv.hi)
			cc_grad = one(T)/y.Intv.hi*mid_grad(x.cc_grad, x.cv_grad, pos1)
							- x.Intv.hi/(y.Intv.lo*y.Intv.hi)*mid_grad(y.cc_grad, y.cv_grad, pos2)
		end
		cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
		cnst::Bool = y.cnst ? x.cnst : (x.cnst ? y.cnst : x.cnst || y.cnst)
		return SMCg{T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x.IntvBox,x.xref)
	=#
	else
		#println("ran me to 3")
		return x*inv(y)
	end
end


function sqr_cv_NS(x::T,xL::T,xU::T) where {T<:AbstractFloat}
	return x^2,two(T)*x
end
function sqr_cc(x::T,xL::T,xU::T) where {T<:AbstractFloat}
	if (xU>xL)
		return xL^2 + (xL+xU)*(x-xL),(xL+xU)
	else
		return xU^2,zero(T)
	end
end
function sqr_cv(x::T,xL::T,xU::T) where {T<:AbstractFloat}
	if (zero(x) <= xL || xU <= zero(x))
		return x^2,two(T)*x
	elseif ((xL<zero(x)) && (zero(x) <= xU))
		return (x^3)/xU,(three(T)*x^2)/xU
	else
		return (x^3)/xL,(three(T)*x^2)/xL
	end
end
function sqr(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	Intv::V = x.Intv^2
    xL::T = x.Intv.lo
    xU::T = x.Intv.hi
    xLc::T = Intv.lo
    xUc::T = Intv.hi
    eps_max::T = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.hi : x.Intv.lo
    eps_min::T = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.lo : x.Intv.hi
    midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
    midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
    cc::T,dcc::T = sqr_cc(midcc,x.Intv.lo,x.Intv.hi)
    cv::T,dcv::T = sqr_cv_NS(midcv,x.Intv.lo,x.Intv.hi)
    if (MC_param.mu >= 1)
      gcc1::T,gdcc1::T = sqr_cc(x.cv,x.Intv.lo,x.Intv.hi)
  	  gcv1::T,gdcv1::T = sqr_cv(x.cv,x.Intv.lo,x.Intv.hi)
  	  gcc2::T,gdcc2::T = sqr_cc(x.cc,x.Intv.lo,x.Intv.hi)
  	  gcv2::T,gdcv2::T = sqr_cv(x.cc,x.Intv.lo,x.Intv.hi)
  	  cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
  	  cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
    else
      cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
      cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
      cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
    end
    return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv,x.cnst,x.IntvBox, x.xref)
end

########### Defines functions required for linear algebra packages

one(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = SMCg{N,V,T}(one(T),one(T),zeros(SVector{N,T}),zeros(SVector{N,T}),V(one(T)),x.cnst,x.IntvBox,x.xref)

zero(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = SMCg{N,V,T}(zero(x.cc),zero(x.cv),zeros(SVector{N,T}),zeros(SVector{N,T}),V(zero(T)),x.cnst,x.IntvBox,x.xref)

real(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = x

dist(x1::SMCg{N,V,T}, x2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = max(abs(x1.cc-x2.cc), abs(x1.cv-x2.cv))

eps(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = max(eps(x.cc), eps(x.cv))

mid(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = mid(Intv(x.Intv))


#=
###### Defines boolean operators

"""==(x::SMC,y::SMC) defines == for SMC type
"""
function ==(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv == y.Intv && x.cv == y.cv && x.cc == y.cc
end
"""!=(x::SMC,y::SMC) defines != for SMC type
"""
function !=(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv != y.Intv || x.cv != y.cv || x.cc != y.cc
end
"""<=(x::SMC,y::SMC) defines <= for SMC type
"""
function <=(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv <= y.Intv && x.cv <= y.cv && x.cc <= y.cc
end
"""=>(x::SMC,y::SMC) defines => for SMC type
"""
function >=(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv >= y.Intv && x.cv >= y.cv && x.cc >= y.cc
end
""">(x::SMC,y::SMC) defines > for SMC type
"""
function >(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv > y.Intv && x.cv > y.cv && x.cc > y.cc
end
"""<(x::SMC,y::SMC) defines < for SMC type
"""
function <(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv < y.Intv && x.cv < y.cv && x.cc < y.cc
end
=#
