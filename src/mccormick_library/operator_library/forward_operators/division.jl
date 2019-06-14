div_alphaxy(es, nu, x, y) = (es/y.hi) + (x.lo/(y.lo*y.hi))*(y.hi-nu)
div_gammay(omega, y) = (y.lo*(max(0.0, omega))^2)/(y.hi - omega*(y.hi-y.lo))
div_deltaxy(omega, x, y) = (1.0/(y.hi*y.lo))*(x.hi - x.lo)*(y.hi - y.lo)*div_gammay(omega, y)
div_psixy(es, nu, x, y) = div_alphaxy(es, nu, x, y) + div_deltaxy(((es - x.lo)/(x.hi - x.lo))-((nu - y.lo)/(y.hi - y.lo)), x, y)
div_lambdaxy(es, nu, x) = (((eps + sqrt(x.lo*x.hi))/(sqrt(x.lo) + sqrt(x.hi)))^2)/nu
div_omegaxy(x, y) = (y.hi/(y.hi-y.lo))*(1.0 - sqrt((y.lo*(x.hi-x.lo))/((-x.lo)*(y.hi-y.lo)+(y.lo)*(x.hi-x.lo))))
div_nuline(x, y, z) = y.lo + (y.hi - y.lo)*(z - x.lo)/(x.hi - x.lo)

function div_diffcv(x::MC, y::MC)
	nu_bar = div_nuline(x.Intv, y.Intv, x.cv)
	if (0.0 <= x.Intv.lo)
		return div_lambdaxy(x.cv, y.cc, x)
	elseif (x.Intv.lo < 0.0) && (nu_bar <= y.cv)
		return div_alphaxy(x.cv, y.cv, x.Intv, y.Intv)
	elseif (x.Intv.lo < 0.0) && (nu_bar > y.cv)
		return div_psixy(x.cv, mid(y.cv, y.cc, nu_bar-(y.Intv.hi - y.Intv.lo)*div_omegaxy(x.Intv, y.Intv)), x.Intv, y.Intv)
	end
end

function div_MV(x::MC{N},y::MC{N}) where N
	if (0.0 < y.Intv.lo)
		cv = div_diffcv(x, y)
		cc = -div_diffcv(-x, y)
	elseif (y.Intv.hi < 0.0)
		cv = div_diffcv(-x, -y)
		cc = -div_diffcv(-x, y)
	else
		DomainError()
	end
	MC{N}(cv, cc, x.Intv/y.Intv, cv_grad, cc_grad, (x.cnst && y.cnst))
end

function /(x::MC{N},y::MC{N}) where N
	pos_orth::Bool = (x.Intv.lo >= 0) && (y.Intv.lo >= 0)
	degen1 = ((x.Intv.hi - x.Intv.lo) == 0.0)
	degen2 = ((y.Intv.hi - y.Intv.lo) == 0.0)
	if (x==y)
		return one(Float64)
	#elseif  (MC_param.mu >= 1 && ~(degen1||degen2))
	#	println("MC_param.mu: $(MC_param.mu)")
	#	return div_MV(x,y)
	#=
	elseif (MC_param.multivar_refine) && (~ (MC_param.mu >= 1)) && (pos_orth)
		#println("ran me to 2")
		Intv::V = x.Intv/y.Intv
		xLc::T = .lo
		xUc::T = .hi
		cv1::T,pos1::Int64 = mid3(x.cv,x.cc,x.lo)
		cv1t::T = (cv1+sqrt(x.lo*x.hi))
				  /(sqrt(x.lo)+sqrt(x.hi))
		cv2t::T,pos2::Int64 = mid3(y.cv,y.cc,y.hi)
		cv::T = sqr(cv1t)/cv2t
		cv_grad::SVector{N,T} = 2.0*(cv1t/cv2t)/(sqrt(x.lo)+sqrt(x.hi))*
								  mid_grad(x.cc_grad, x.cv_grad, pos1) - ((cv1t/cv2t)^2)*
									mid_grad(y.cc_grad, y.cv_grad, pos2)
		cc1::T,pos1 = mid3(x.cv,x.cc,x.hi)
		cc2::T,pos2 = mid3(y.cv,y.cc,y.lo)
		gcc1::T = y.hi*cc1 - x.lo*cc2 + x.lo*y.lo
		gcc2::T = y.lo*cc1 - x.hi*cc2 + x.hi*y.hi
		if gcc1 <= gcc2
			cc::T = gcc1/(y.lo*y.hi)
			cc_grad::SVector{N,T} = one(T)/y.lo*mid_grad(x.cc_grad, x.cv_grad, pos1)
							- x.lo/(y.lo*y.hi)*mid_grad(y.cc_grad, y.cv_grad, pos2)
		else
			cc = gcc2/(y.lo*y.hi)
			cc_grad = one(T)/y.hi*mid_grad(x.cc_grad, x.cv_grad, pos1)
							- x.hi/(y.lo*y.hi)*mid_grad(y.cc_grad, y.cv_grad, pos2)
		end
		cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
		cnst::Bool = y.cnst ? x.cnst : (x.cnst ? y.cnst : x.cnst || y.cnst)
		return SMCg{T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x.IntvBox,x.xref)
	=#
	else
		return x*inv(y)
	end
end
