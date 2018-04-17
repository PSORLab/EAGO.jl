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

function sqr(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	Intv::V = x.Intv^2
	xL::T = x.Intv.lo
	xU::T = x.Intv.hi
	xLc::T = xL*xL
	xUc::T = xU*xU
	zT::T = zero(T)

	#calculates: eps_min = mid3(xLc,xUc,zT)
	if (xL>=xU)
		if (xU<=zT)
			if (xL>=zT)
				eps_min::T = zT
			else
				eps_min = xL
			end
		else
			eps_min = xU
		end
	elseif (xU>=zT)
		if (xL<=zT)
			eps_min = zT
		else
			eps_min = xL
		end
	else
		eps_min = xU
	end
	comp_LU::Bool = (abs(xLc)>=abs(xUc))
  	eps_max::T = ifelse(comp_LU,xL,xU)

	if (MC_param.mu >= 1)

		# calculate mid3(x.cc,x.cv,eps_max)
		if (eps_max<x.cv)
			cc::T = xLc+(xL+xU)*(x.cc-xL)
		elseif (eps_max>x.cc)
			cc = xLc+(xL+xU)*(x.cv-xL)
		else
			cc = xLc+(xL+xU)*(eps_max-xL)
		end
		if comp_LU
			cc_grad::SVector{N,T} = (xL+xU)*x.cc_grad
		else
			cc_grad = (xL+xU)*x.cv_grad
		end

		# calculate mid3(x.cc,x.cv,min)
		if (eps_min<x.cv)
			midcv::T = x.cv
		elseif (eps_min>x.cc)
			midcv = x.cc
		else
			midcv = eps_min
		end
		if (zT<=xL) 		# increasing
			cv::T = midcv^2
			cv_grad::SVector{N,T} = (two(T)*midcv)*x.cv_grad
		elseif (xU<=zT) 	# decreasing
			cv = midcv^2
			cv_grad = (two(T)*midcv)*x.cc_grad
		elseif (zT<=midcv)  # increasing
			cv = (midcv^3)/xU
			cv_grad = ((three(T)/xU)*midcv^2)*x.cv_grad
		elseif (zT>midcv)   # decreasing
			cv = (midcv^3)/xL
			cv_grad = ((three(T)/xL)*midcv^2)*x.cc_grad
		end
	else
		# calculate mid3(x.cc,x.cv,min)
		if (eps_min<x.cv)
			cv = x.cv^2
			cv_grad = (two(T)*x.cv)*x.cv_grad
		elseif (eps_min>x.cc)
			cv = x.cc^2
			cv_grad = (two(T)*x.cc)*x.cc_grad
		else
			cv = eps_min^2
			cv_grad = zeros(SVector{N,T})
		end
		# calculate mid3(x.cc,x.cv,eps_max)
		m::T = (xUc-xLc)/(xU-xL)
		if (eps_max<x.cv)
			cc = xLc + m*(x.cc-xL)
			cc_grad = m*x.cc_grad
		elseif (eps_max>x.cc)
			cc = xLc + m*(x.cv-xL)
			cc_grad = m*x.cv_grad
		else
			cc = xLc + m*(eps_max-xL)
			cc_grad = zeros(SVector{N,T})
		end
		cv,cc,cv_grad,cc_grad = cut(Intv.lo,Intv.hi,cv,cc,cv_grad,cc_grad)
	end
  return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, x.Intv^2, x.cnst, x.IntvBox,x.xref)
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
