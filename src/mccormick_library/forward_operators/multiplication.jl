########### Defines multiplicaiton of SMC and SMC
#=
function gCxAcv(alpha, beta ,lambda ,nu ,x1::MC{N},x2::MC{N}) where N

		# gCxA pre-terms
		alplo = alpha.lo
		alphi = alpha.hi
		betlo = beta.lo
		bethi = beta.hi
		LmdDel = lambda.hi-lambda.lo
		NuDel = nu.hi-nu.lo
		LmdSum = lambda.lo+lambda.hi
		NuSum = nu.lo+nu.hi
		muT =  MC_param.mu_flt
		mu1 = MC_param.mu + 1
		mu1T = MC_param.mu_flt + 1.0
		mu1n = MC_param.mu - 1

		xslo = lambda.lo+LmdDel*(((nu.hi-betlo)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT))
		xshi = lambda.lo+LmdDel*(((nu.hi-bethi)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT))
		yslo = nu.lo+NuDel*(((lambda.hi-alplo)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))
		yshi = nu.lo+NuDel*(((lambda.hi-alphi)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))

		# calculates term 1
		if (xslo <= alplo)
			term1 = alplo
		elseif (alphi <= xslo)
			term1 = alphi
		else
			term1 = xslo
		end

		# calculates term 2
		if (xshi <= alplo)
			term2 = alplo
		elseif (alphi <= xshi)
			term2 = alphi
		else
			term2 = xshi
		end

		# calculates term 3
		if (yslo <= betlo)
			term3 = betlo
		elseif (bethi <= yslo)
			term3 = bethi
		else
			term3 = yslo
		end

		# calculates term 4
		if (yshi <= betlo)
			term4 = betlo
		elseif (bethi <= yshi)
			term4 = bethi
		else
			term4 = yshi
		end

		# calculates the convex relaxation
		tempGxA1a = LmdDel*NuDel*abs((betlo-nu.lo)/NuDel-(lambda.hi-term1)/LmdDel)^mu1
		tempGxA1 = 0.5*(term1*NuSum+betlo*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA1a)
		tempGxA2a = LmdDel*NuDel*abs((bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel)^mu1
		tempGxA2 = 0.5*(term2*NuSum+bethi*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA2a)
		tempGxA3a = LmdDel*NuDel*abs((term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel)^mu1
		tempGxA3 = 0.5*(alplo*NuSum+term3*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA3a)
		tempGxA4a = LmdDel*NuDel*abs((term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel)^mu1
		tempGxA4 = 0.5*(alphi*NuSum+term4*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA4a)

		# gets minima which is cv/cc/Intv depending on arguments
		a = min(tempGxA1,tempGxA2,tempGxA3,tempGxA4)

		if (a == tempGxA1)
				psi_mul = (betlo-nu.lo)/NuDel - (lambda.hi-term1)/LmdDel
				psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1 = 0.5*(LmdSum+mu1T*(LmdDel)*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA2)
				psi_mul = (bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel
				psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA3)
				psi_mul = (term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel
				psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		else
				psi_mul = (term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel
				psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		end

    	grad = max(0.0,psi_mult_x1)*x1.cv_grad+
			   min(0.0,psi_mult_x1)*x1.cc_grad+
			   max(0.0,psi_mult_y1)*x2.cv_grad+
			   min(0.0,psi_mult_y1)*x2.cc_grad

    return a,grad
end

function gCxAcc(alpha, beta, lambda, nu, x1::MC{N}, x2::MC{N}) where N

		# gCxA pre-terms
		alplo = alpha.lo
		alphi = alpha.hi
		betlo = beta.lo
		bethi = beta.hi
		LmdDel = lambda.hi-lambda.lo
		NuDel = nu.hi-nu.lo
		LmdSum = lambda.lo+lambda.hi
		NuSum = nu.lo+nu.hi
		NuDotLmd = lambda.lo*nu.lo+lambda.hi*nu.hi
		muT = MC_param.mu_flt
		mu1 = MC_param.mu + 1
		mu1n = MC_param.mu - 1
		mu1T = MC_param.mu + 1.0

		xslo = lambda.lo+LmdDel*(((nu.hi-betlo)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT))
		xshi = lambda.lo+LmdDel*(((nu.hi-bethi)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT))
		yslo = nu.lo+NuDel*(((lambda.hi-alplo)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))
		yshi = nu.lo+NuDel*(((lambda.hi-alphi)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))

		# calculates term 1
		if (xslo <= alplo)
			term1 = alplo
		elseif (alphi <= xslo)
			term1 = alphi
		else
			term1 = xslo
		end

		# calculates term 2
		if (xshi <= alplo)
			term2 = alplo
		elseif (alphi <= xshi)
			term2 = alphi
		else
			term2 = xshi
		end

		# calculates term 3
		if (yslo <= betlo)
			term3 = betlo
		elseif (bethi <= yslo)
			term3 = bethi
		else
			term3 = yslo
		end

		# calculates term 4
		if (yshi <= betlo)
			term4 = betlo
		elseif (bethi <= yshi)
			term4 = bethi
		else
			term4 = yshi
		end

		# calculates the convex relaxation
		tempGxA1a = LmdDel*NuDel*abs((betlo-nu.lo)/NuDel-(lambda.hi-term1)/LmdDel)^mu1
		tempGxA1 = 0.5*(term1*NuSum+betlo*LmdSum-NuDotLmd+tempGxA1a)
		tempGxA2a = LmdDel*NuDel*abs((bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel)^mu1
		tempGxA2 = 0.5*(term2*NuSum+bethi*LmdSum-NuDotLmd+tempGxA2a)
		tempGxA3a = LmdDel*NuDel*abs((term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel)^mu1
		tempGxA3 = 0.5*(alplo*NuSum+term3*LmdSum-NuDotLmd+tempGxA3a)
		tempGxA4a = LmdDel*NuDel*abs((term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel)^mu1
		tempGxA4 = 0.5*(alphi*NuSum+term4*LmdSum-NuDotLmd+tempGxA4a)

		# gets minima which is cv/cc/Intv depending on arguments
		a = min(tempGxA1,tempGxA2,tempGxA3,tempGxA4)

	   if (a == tempGxA1)
		   psi_mul = (betlo-nu.lo)/NuDel - (lambda.hi-term1)/LmdDel
		   psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
		   psi_mult_y1 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA2)
			psi_mul = (bethi-nu.lo)/NuDel - (lambda.hi-term2)/LmdDel
			psi_mult_x1 = 0.5*(NuSum+mu1*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
			psi_mult_y1 = 0.5*(LmdSum+mu1*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA3)
			psi_mul = (term3-nu.lo)/NuDel - (lambda.hi-alplo)/LmdDel
			psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
			psi_mult_y1 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		else
			psi_mul= (term4-nu.lo)/NuDel - (lambda.hi-alphi)/LmdDel
			psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
			psi_mult_y1 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		end
		grad2::SVector{N,Float64} = max(0.0,psi_mult_x1)*x1.cc_grad+
							        min(0.0,psi_mult_x1)*x1.cv_grad-
							        max(0.0,psi_mult_y1)*x2.cc_grad-
							        min(0.0,psi_mult_y1)*x2.cv_grad
    return a,grad2
end

function gCxAIntv(alpha,beta,lambda,nu,x1::MC{N},x2::MC{N}) where N

		# gCxA pre-terms
		alplo = alpha.lo
		alphi = alpha.hi
		betlo = beta.lo
		bethi = beta.hi
		LmdDel = lambda.hi-lambda.lo
		NuDel = nu.hi-nu.lo
		LmdSum = lambda.lo+lambda.hi
		NuSum = nu.lo+nu.hi
		NuLmd = lambda.lo*nu.lo+lambda.hi*nu.hi
		mu1 = MC_param.mu + 1
		mu1n = MC_param.mu - 1
		mu1T = MC_param.mu_flt + 1.0

		xslo = lambda.lo+LmdDel*(((nu.hi-betlo)/NuDel)+sigu(-(NuSum)/(mu1*(NuDel)),mu1T))
		xshi = lambda.lo+LmdDel*(((nu.hi-bethi)/NuDel)+sigu(-(NuSum)/(mu1*(NuDel)),mu1T))
		yslo = nu.lo+NuDel*(((lambda.hi-alplo)/LmdDel)-sigu((LmdSum)/(mu1*LmdDel),mu1T))
		yshi = nu.lo+NuDel*(((lambda.hi-alphi)/LmdDel)-sigu((LmdSum)/(mu1*LmdDel),mu1T))

		# calculates term 1
		if (xslo <= alplo)
			term1 = alplo
		elseif (alphi <= xslo)
			term1 = alphi
		else
			term1 = xslo
		end

		# calculates term 2
		if (xshi <= alplo)
			term2 = alplo
		elseif (alphi <= xshi)
			term2 = alphi
		else
			term2 = xshi
		end

		# calculates term 3
		if (yslo <= betlo)
			term3 = betlo
		elseif (bethi <= yslo)
			term3 = bethi
		else
			term3 = yslo
		end

		# calculates term 4
		if (yshi <= betlo)
			term4 = betlo
		elseif (bethi <= yshi)
			term4 = bethi
		else
			term4 = yshi
		end

		# calculates the convex relaxation
		tempGxA1a = LmdDel*NuDel*abs((betlo-nu.lo)/NuDel-(lambda.hi-term1)/LmdDel)^mu1
		tempGxA1 =  0.5*(term1*NuSum+betlo*LmdSum-NuLmd+tempGxA1a)
		tempGxA2a = LmdDel*NuDel*abs((bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel)^mu1
		tempGxA2 =  0.5*(term2*NuSum+bethi*LmdSum-NuLmd+tempGxA2a)
		tempGxA3a = LmdDel*NuDel*abs((term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel)^mu1
		tempGxA3 =  0.5*(alplo*NuSum+term3*LmdSum-NuLmd+tempGxA3a)
		tempGxA4a = LmdDel*NuDel*abs((term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel)^mu1
		tempGxA4 = 0.5*(alphi*NuSum+term4*LmdSum-NuLmd+tempGxA4a)

		# gets minima which is cv/cc/Intv depending on arguments
		return min(tempGxA1,tempGxA2,tempGxA3,tempGxA4)
end

function multiply_MV(x1::MC{N},x2::MC{N}) where N

	x1cv = x1.cv
	x2cv = x2.cv
	x1cc = x1.cc
	x2cc = x2.cc
	x1lo = x1.Intv.lo
	x2lo = x2.Intv.lo
	x1hi = x1.Intv.hi
	x2hi = x2.Intv.hi
	x1h_l = x1hi-x1lo
	x1h_v = x1hi-x1cv
	x2h_l = x2hi-x2lo
	x2v_l = x2cv-x2lo

	cv = 0.0
	cc = 0.0
	mu1 = MC_param.mu + 1
	muf1 = MC_param.mu_flt + 1.0
	alpha0 = IntervalType(x1cv,x1cc)
	beta0 = IntervalType(x2cv,x2cc)
	if (0.0 <= x1lo) && (0.0 <= x2lo)
		b1_temp1 = max(0.0,((x2v_l/x2h_l)-(x1h_v/x1h_l)))
		b1_term2 = x2lo+muf1*x2h_l*b1_temp1^MC_param.mu
		b1_term3 = x1lo+muf1*x1h_l*b1_temp1^MC_param.mu
		cv = x1cv*x2lo + x1lo*x2v_l + x1h_l*x2h_l*b1_temp1^mu1
		cv_grad::SVector{N,Float64} = b1_term2*x1.cv_grad + b1_term3*x2.cv_grad
	elseif ((x1hi <= 0.0)) && (x2hi <= 0.0)
		b1_temp = max(0.0,(x2hi - x2cc)/x2h_l - (x1cc-x1lo)/(x1h_l))
		cv = x1cc*x2hi + x1hi*(x2cc - x2hi) + (x1h_l)*(x2h_l)*b1_temp^mu1
		b1_term2 = -x2hi+muf1*(x2h_l)*max(0.0, b1_temp)^MC_param.mu
		b1_term3 = -x1hi+muf1*(x1h_l)*max(0.0, b1_temp)^MC_param.mu
		cv_grad = (-b1_term2)*x1.cc_grad - b1_term3*x2.cc_grad
	else
		cv,cv_grad = gCxAcv(alpha0, beta0, x1.Intv, x2.Intv, x1, x2)
	end
	if ((x1hi <= 0.0) && (0.0 <= x2lo))
		btemp1 = max(0.0,x2v_l/x2h_l - (x1cc-x1lo)/(x1h_l))
		btemp2 = x2lo+muf1*x2h_l*btemp1^MC_param.mu
		btemp3 = -x1hi+muf1*x1h_l*btemp1^MC_param.mu
		cc = x1cc*x2lo-x1hi*(x2lo-x2cv)-x1h_l*x2h_l*btemp1^mu1
		cc_grad::SVector{N,Float64} = btemp2*x1.cc_grad - btemp3*x2.cv_grad
	elseif ((0.0<=x1.Intv.lo) && (x2.Intv.hi<=0.0))
		btemp1 = max(0.0,(x2hi-x2cc)/x2h_l - x1h_v/x1h_l)
		btemp2 = -x2hi+muf1*(x2h_l)*btemp1^MC_param.mu
		btemp3 = x1lo+muf1*x1h_l*btemp1^MC_param.mu
		cc = x1cv*x2hi+x1lo*(x2cc-x2hi)-x1h_l*x2h_l*btemp1^mu1
		cc_grad = -btemp2*x1.cv_grad + btemp3*x2.cc_grad
	else
		cct,cc_gradt::SVector{N,Float64} = gCxAcc(-alpha0,beta0,-x1.Intv,x2.Intv,x1,x2)
		cc = -cct
		cc_grad = cc_gradt
	end
	if (min(x1lo,x2lo) < 0.0 < max(x1hi,x2hi))
		lo_Intv_calc = gCxAIntv(x1.Intv,x2.Intv,x1.Intv,x2.Intv,x1,x2)
		hi_Intv_calc = -gCxAIntv(-x1.Intv,x2.Intv,-x1.Intv,x2.Intv,x1,x2)
		Intv_calc = IntervalType(lo_Intv_calc,hi_Intv_calc)
	else
		Intv_calc = x1.Intv*x2.Intv
	end
	cnst = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	return MC{N}(cv, cc, Intv_calc, cv_grad, cc_grad, cnst)
end
=#

function sigu(x::Float64, mu1T::Float64)
	 (0.0 <= x) ? x^(1.0/mu1T) : -abs(x)^(1.0/mu1T)
 end

function gCxAcv(alpha::IntervalType,beta::IntervalType,lambda::IntervalType,nu::IntervalType,x1::MC{N},x2::MC{N}) where N
		# gCxA pre-terms
		alplo::Float64 = alpha.lo
		alphi::Float64 = alpha.hi
		betlo::Float64 = beta.lo
		bethi::Float64 = beta.hi
		LmdDel::Float64 = lambda.hi-lambda.lo
		NuDel::Float64 = nu.hi-nu.lo
		LmdSum::Float64 = lambda.lo+lambda.hi
		NuSum::Float64 = nu.lo+nu.hi
		muT::Float64 = convert(Float64,MC_param.mu)
		mu1::Int = MC_param.mu+1
		mu1T::Float64 = convert(Float64,MC_param.mu+1)
		mu1n::Int = MC_param.mu-1

		xslo::Float64 = lambda.lo+LmdDel*(((nu.hi-betlo)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT)) # correct to minus?, correct internal minus? exponent to mu
		xshi::Float64 = lambda.lo+LmdDel*(((nu.hi-bethi)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT)) # correct to minus?
		yslo::Float64 = nu.lo+NuDel*(((lambda.hi-alplo)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))
		yshi::Float64 = nu.lo+NuDel*(((lambda.hi-alphi)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))

		# calculates term 1
		if (xslo <= alplo)
			term1::Float64 = alplo
		elseif (alphi <= xslo)
			term1 = alphi
		else
			term1 = xslo
		end

		# calculates term 2
		if (xshi <= alplo)
			term2::Float64 = alplo
		elseif (alphi <= xshi)
			term2 = alphi
		else
			term2 = xshi
		end

		# calculates term 3
		if (yslo <= betlo)
			term3::Float64 = betlo
		elseif (bethi <= yslo)
			term3 = bethi
		else
			term3 = yslo
		end

		# calculates term 4
		if (yshi <= betlo)
			term4::Float64 = betlo
		elseif (bethi <= yshi)
			term4 = bethi
		else
			term4 = yshi
		end

		# calculates the convex relaxation
		tempGxA1a::Float64 = LmdDel*NuDel*abs((betlo-nu.lo)/NuDel-(lambda.hi-term1)/LmdDel)^mu1
		tempGxA1::Float64 = 0.5*(term1*NuSum+betlo*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA1a)
		tempGxA2a::Float64 = LmdDel*NuDel*abs((bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel)^mu1
		tempGxA2::Float64 = 0.5*(term2*NuSum+bethi*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA2a)
		tempGxA3a::Float64 = LmdDel*NuDel*abs((term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel)^mu1
		tempGxA3::Float64 = 0.5*(alplo*NuSum+term3*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA3a)
		tempGxA4a::Float64 = LmdDel*NuDel*abs((term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel)^mu1
		tempGxA4::Float64 = 0.5*(alphi*NuSum+term4*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA4a)

		# gets minima which is cv/cc/Intv depending on arguments
		a = min(tempGxA1,tempGxA2,tempGxA3,tempGxA4)

		if (a == tempGxA1)
				psi_mul::Float64 = (betlo-nu.lo)/NuDel - (lambda.hi-term1)/LmdDel
				psi_mult_x1::Float64 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1::Float64 = 0.5*(LmdSum+mu1T*(LmdDel)*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA2)
				psi_mul = (bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel
				psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA3)
				psi_mul = (term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel
				psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		else
				psi_mul = (term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel
				psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		end

    	grad = max(0.0, psi_mult_x1)*x1.cv_grad+
			   min(0.0, psi_mult_x1)*x1.cc_grad+
			   max(0.0, psi_mult_y1)*x2.cv_grad+
			   min(0.0, psi_mult_y1)*x2.cc_grad

    return a,grad
end

function gCxAcc(alpha::IntervalType,beta::IntervalType,lambda::IntervalType,nu::IntervalType,x1::MC{N},x2::MC{N}) where N
		# gCxA pre-terms
		alplo::Float64 = alpha.lo
		alphi::Float64 = alpha.hi
		betlo::Float64 = beta.lo
		bethi::Float64 = beta.hi
		LmdDel::Float64 = lambda.hi-lambda.lo
		NuDel::Float64 = nu.hi-nu.lo
		LmdSum::Float64 = lambda.lo+lambda.hi
		NuSum::Float64 = nu.lo+nu.hi
		NuDotLmd::Float64 = lambda.lo*nu.lo+lambda.hi*nu.hi
		muT::Float64 = convert(Float64,MC_param.mu)
		mu1::Int = MC_param.mu+1
		mu1n::Int = MC_param.mu-1
		mu1T::Float64 = convert(Float64,MC_param.mu+1)

		xslo::Float64 = lambda.lo+LmdDel*(((nu.hi-betlo)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT))
		xshi::Float64 = lambda.lo+LmdDel*(((nu.hi-bethi)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT))
		yslo::Float64 = nu.lo+NuDel*(((lambda.hi-alplo)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))
		yshi::Float64 = nu.lo+NuDel*(((lambda.hi-alphi)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))

		# calculates term 1
		if (xslo <= alplo)
			term1::Float64 = alplo
		elseif (alphi <= xslo)
			term1 = alphi
		else
			term1 = xslo
		end

		# calculates term 2
		if (xshi <= alplo)
			term2::Float64 = alplo
		elseif (alphi <= xshi)
			term2 = alphi
		else
			term2 = xshi
		end

		# calculates term 3
		if (yslo <= betlo)
			term3::Float64 = betlo
		elseif (bethi <= yslo)
			term3 = bethi
		else
			term3 = yslo
		end

		# calculates term 4
		if (yshi <= betlo)
			term4::Float64 = betlo
		elseif (bethi <= yshi)
			term4 = bethi
		else
			term4 = yshi
		end

		# calculates the convex relaxation
		tempGxA1a::Float64 = LmdDel*NuDel*abs((betlo-nu.lo)/NuDel-(lambda.hi-term1)/LmdDel)^mu1
		tempGxA1::Float64 = 0.5*(term1*NuSum+betlo*LmdSum-NuDotLmd+tempGxA1a)
		tempGxA2a::Float64 = LmdDel*NuDel*abs((bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel)^mu1
		tempGxA2::Float64 = 0.5*(term2*NuSum+bethi*LmdSum-NuDotLmd+tempGxA2a)
		tempGxA3a::Float64 = LmdDel*NuDel*abs((term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel)^mu1
		tempGxA3::Float64 = 0.5*(alplo*NuSum+term3*LmdSum-NuDotLmd+tempGxA3a)
		tempGxA4a::Float64 = LmdDel*NuDel*abs((term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel)^mu1
		tempGxA4::Float64 = 0.5*(alphi*NuSum+term4*LmdSum-NuDotLmd+tempGxA4a)

		# gets minima which is cv/cc/Intv depending on arguments
		a::Float64 = min(tempGxA1,tempGxA2,tempGxA3,tempGxA4)

	   if (a == tempGxA1)
		   psi_mul::Float64 = (betlo-nu.lo)/NuDel - (lambda.hi-term1)/LmdDel
		   psi_mult_x1::Float64 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
		   psi_mult_y1::Float64 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA2)
			psi_mul = (bethi-nu.lo)/NuDel - (lambda.hi-term2)/LmdDel
			psi_mult_x1 = 0.5*(NuSum+mu1*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
			psi_mult_y1 = 0.5*(LmdSum+mu1*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA3)
			psi_mul = (term3-nu.lo)/NuDel - (lambda.hi-alplo)/LmdDel
			psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
			psi_mult_y1 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		else
			psi_mul= (term4-nu.lo)/NuDel - (lambda.hi-alphi)/LmdDel
			psi_mult_x1 = 0.5*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
			psi_mult_y1 = 0.5*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		end
		grad2 = max(0.0,psi_mult_x1)*x1.cc_grad+
				min(0.0,psi_mult_x1)*x1.cv_grad-
				max(0.0,psi_mult_y1)*x2.cc_grad-
				min(0.0,psi_mult_y1)*x2.cv_grad
    return a,grad2
end

function gCxAIntv(alpha::IntervalType,beta::IntervalType,lambda::IntervalType,nu::IntervalType,x1::MC{N},x2::MC{N}) where N
		# gCxA pre-terms
		alplo::Float64 = alpha.lo
		alphi::Float64 = alpha.hi
		betlo::Float64 = beta.lo
		bethi::Float64 = beta.hi
		LmdDel::Float64 = lambda.hi-lambda.lo
		NuDel::Float64 = nu.hi-nu.lo
		LmdSum::Float64 = lambda.lo+lambda.hi
		NuSum::Float64 = nu.lo+nu.hi
		NuLmd::Float64 = lambda.lo*nu.lo+lambda.hi*nu.hi
		mu1::Int = MC_param.mu+1
		mu1n::Int = MC_param.mu-1
		mu1T::Float64 = convert(Float64,MC_param.mu+1)

		xslo::Float64 = lambda.lo+LmdDel*(((nu.hi-betlo)/NuDel)+sigu(-(NuSum)/(mu1*(NuDel)),mu1T))
		xshi::Float64 = lambda.lo+LmdDel*(((nu.hi-bethi)/NuDel)+sigu(-(NuSum)/(mu1*(NuDel)),mu1T))
		yslo::Float64 = nu.lo+NuDel*(((lambda.hi-alplo)/LmdDel)-sigu((LmdSum)/(mu1*LmdDel),mu1T))
		yshi::Float64 = nu.lo+NuDel*(((lambda.hi-alphi)/LmdDel)-sigu((LmdSum)/(mu1*LmdDel),mu1T))

		# calculates term 1
		if (xslo <= alplo)
			term1::Float64 = alplo
		elseif (alphi <= xslo)
			term1 = alphi
		else
			term1 = xslo
		end

		# calculates term 2
		if (xshi <= alplo)
			term2::Float64 = alplo
		elseif (alphi <= xshi)
			term2 = alphi
		else
			term2 = xshi
		end

		# calculates term 3
		if (yslo <= betlo)
			term3::Float64 = betlo
		elseif (bethi <= yslo)
			term3 = bethi
		else
			term3 = yslo
		end

		# calculates term 4
		if (yshi <= betlo)
			term4::Float64 = betlo
		elseif (bethi <= yshi)
			term4 = bethi
		else
			term4 = yshi
		end

		# calculates the convex relaxation
		tempGxA1a::Float64 = LmdDel*NuDel*abs((betlo-nu.lo)/NuDel-(lambda.hi-term1)/LmdDel)^mu1
		tempGxA1::Float64 =  0.5*(term1*NuSum+betlo*LmdSum-NuLmd+tempGxA1a)
		tempGxA2a::Float64 = LmdDel*NuDel*abs((bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel)^mu1
		tempGxA2::Float64 =  0.5*(term2*NuSum+bethi*LmdSum-NuLmd+tempGxA2a)
		tempGxA3a::Float64 = LmdDel*NuDel*abs((term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel)^mu1
		tempGxA3::Float64 =  0.5*(alplo*NuSum+term3*LmdSum-NuLmd+tempGxA3a)
		tempGxA4a::Float64 = LmdDel*NuDel*abs((term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel)^mu1
		tempGxA4::Float64 = 0.5*(alphi*NuSum+term4*LmdSum-NuLmd+tempGxA4a)

		# gets minima which is cv/cc/Intv depending on arguments
		return min(tempGxA1,tempGxA2,tempGxA3,tempGxA4)
end

function multiply_MV(x1::MC{N}, x2::MC{N}, z::Interval{Float64}) where {N}
	x1cv::Float64 = x1.cv
	x2cv::Float64 = x2.cv
	x1cc::Float64 = x1.cc
	x2cc::Float64 = x2.cc
	x1lo::Float64 = x1.Intv.lo
	x2lo::Float64 = x2.Intv.lo
	x1hi::Float64 = x1.Intv.hi
	x2hi::Float64 = x2.Intv.hi
	x1h_l::Float64 = x1hi-x1lo
	x1h_v::Float64 = x1hi-x1cv
	x2h_l::Float64 = x2hi-x2lo
	x2v_l::Float64 = x2cv-x2lo
	cv::Float64 = 0.0
	cc::Float64 = 0.0
	mu1::Int = MC_param.mu+1
	muf1::Float64 = convert(Float64,MC_param.mu+1)
	alpha0::IntervalType = IntervalType(x1cv,x1cc)
	beta0::IntervalType = IntervalType(x2cv,x2cc)
	if (0.0 <= x1lo) && (0.0 <= x2lo)
		b1_temp1::Float64 = max(0.0,((x2v_l/x2h_l)-(x1h_v/x1h_l)))
		b1_term2::Float64 = x2lo+muf1*x2h_l*b1_temp1^MC_param.mu
		b1_term3::Float64 = x1lo+muf1*x1h_l*b1_temp1^MC_param.mu
		cv = x1cv*x2lo + x1lo*x2v_l + x1h_l*x2h_l*b1_temp1^mu1
		cv_grad = b1_term2*x1.cv_grad + b1_term3*x2.cv_grad
	elseif ((x1hi <= 0.0)) && (x2hi <= 0.0)
		b1_temp = max(0.0,(x2hi - x2cc)/x2h_l - (x1cc-x1lo)/(x1h_l))
		cv = x1cc*x2hi + x1hi*(x2cc - x2hi) + (x1h_l)*(x2h_l)*b1_temp^mu1
		b1_term2 = -x2hi+muf1*(x2h_l)*max(0.0,b1_temp)^MC_param.mu
		b1_term3 = -x1hi+muf1*(x1h_l)*max(0.0,b1_temp)^MC_param.mu
		cv_grad = (-b1_term2)*x1.cc_grad - b1_term3*x2.cc_grad
	else
		cv,cv_grad = gCxAcv(alpha0,beta0,x1.Intv,x2.Intv,x1,x2)
	end
	if ((x1hi <= 0.0)) && ((0.0) <= x2lo)
		btemp1 = max(0.0,x2v_l/x2h_l - (x1cc-x1lo)/(x1h_l))
		btemp2 = x2lo+muf1*x2h_l*btemp1^MC_param.mu
		btemp3 = -x1hi+muf1*x1h_l*btemp1^MC_param.mu
		cc = x1cc*x2lo-x1hi*(x2lo-x2cv)-x1h_l*x2h_l*btemp1^mu1
		cc_grad = btemp2*x1.cc_grad - btemp3*x2.cv_grad
	elseif ((0.0)<=x1.Intv.lo) && (x2.Intv.hi<=0.0)
		btemp1 = max(0.0,(x2hi-x2cc)/x2h_l - x1h_v/x1h_l)
		btemp2 = -x2hi+muf1*(x2h_l)*btemp1^MC_param.mu
		btemp3 = x1lo+muf1*x1h_l*btemp1^MC_param.mu
		cc = x1cv*x2hi+x1lo*(x2cc-x2hi)-x1h_l*x2h_l*btemp1^mu1
		cc_grad = -btemp2*x1.cv_grad + btemp3*x2.cc_grad
	else
		cct::Float64,cc_gradt = gCxAcc(-alpha0,beta0,-x1.Intv,x2.Intv,x1,x2)
		cc = -cct
		cc_grad = cc_gradt
	end
	cnst = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	return MC{N}(cv, cc, z, cv_grad, cc_grad, cnst)
end

function mul1_u1pos_u2pos(x1::MC{N}, x2::MC{N}, z::Interval{Float64}, cnst::Bool) where N
  xLc = z.lo
  xUc = z.hi
  cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2 = x2.Intv.lo*x1.cv + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
  if (cv1 > cv2)
    cv = cv1
    cv_grad = x2.Intv.hi*x1.cv_grad
  else
    cv = cv2
    cv_grad = x2.Intv.lo*x1.cv_grad
  end

  cc1 = x2.Intv.lo*x1.cc + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
  if (cc1 < cc2)
    cc = cc1
    cc_grad = x2.Intv.lo*x1.cc_grad
  else
    cc = cc2
    cc_grad = x2.Intv.hi*x1.cc_grad
  end
  cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  return MC{N}(cv,cc,z,cv_grad,cc_grad,cnst)
end

function mul1_u1pos_u2mix(x1::MC{N}, x2::MC{N}, z::Interval{Float64}, cnst::Bool) where N
  xLc = z.lo
  xUc = z.hi
  cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2 = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
  if (cv1 > cv2)
    cv = cv1
    cv_grad = x2.Intv.hi*x1.cv_grad
  else
    cv = cv2
    cv_grad = x2.Intv.lo*x1.cc_grad
  end

  cc1 = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
  if (cc1 < cc2)
    cc = cc1
    cc_grad = x2.Intv.lo*x1.cv_grad
  else
    cc = cc2
    cc_grad = x2.Intv.hi*x1.cc_grad
  end
  cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  return MC{N}(cv, cc, z, cv_grad, cc_grad, cnst)
end
function mul1_u1mix_u2mix(x1::MC{N}, x2::MC{N}, z::Interval{Float64}, cnst::Bool) where {N}
  xLc = z.lo
  xUc = z.hi
  cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2 = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.lo
  cv = cv1
  if (cv1 > cv2)
    cv = cv1
    cv_grad = x2.Intv.hi*x1.cv_grad
  else
    cv = cv2
    cv_grad = x2.Intv.lo*x1.cc_grad
  end
  cc1 = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.hi
  cc = cc1
  if (cc1 < cc2)
    cc = cc1
    cc_grad = x2.Intv.lo*x1.cv_grad
  else
    cc = cc2
    cc_grad = x2.Intv.hi*x1.cc_grad
  end
  cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  return MC{N}(cv, cc, z, cv_grad, cc_grad, cnst)
end
function mul2_u1pos_u2pos(x1::MC{N}, x2::MC{N}, z::Interval{Float64}) where N
	xLc = z.lo
	xUc = z.hi
	cnst = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
	cv2 = x2.Intv.lo*x1.cv + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
	if (cv1 > cv2)
		cv = cv1
		cv_grad = x2.Intv.hi*x1.cv_grad + x1.Intv.hi*x2.cv_grad
	else
		cv = cv2
		cv_grad = x2.Intv.lo*x1.cv_grad + x1.Intv.lo*x2.cv_grad
	end
	cc1 = x2.Intv.lo*x1.cc + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
	cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
	if (cc1 < cc2)
		cc = cc1
		cc_grad = x2.Intv.lo*x1.cc_grad + x1.Intv.hi*x2.cc_grad
	else
		cc = cc2
		cc_grad = x2.Intv.hi*x1.cc_grad + x1.Intv.lo*x2.cc_grad
	end
	cv, cc, cv_grad, cc_grad = cut(xLc, xUc, cv, cc, cv_grad, cc_grad)
	return MC{N}(cv, cc, z, cv_grad, cc_grad, cnst)
end
function mul2_u1pos_u2mix(x1::MC{N}, x2::MC{N}, z::Interval{Float64}, cnst::Bool) where N
  xLc = z.lo
  xUc = z.hi
  cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2 = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
  if (cv1 > cv2)
    cv = cv1
    cv_grad = x1.Intv.hi*x2.cv_grad
  else
    cv = cv2
    cv_grad = x1.Intv.lo*x2.cc_grad
  end
  cc1 = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
  if (cc1 < cc2)
    cc = cc1
    cc_grad = x1.Intv.hi*x2.cc_grad
  else
    cc = cc2
    cc_grad = x1.Intv.lo*x2.cc_grad
  end
  cv, cc, cv_grad, cc_grad = cut(xLc, xUc, cv, cc, cv_grad, cc_grad)
  return MC{N}(cv, cc, z, cv_grad, cc_grad, cnst)
end
function mul2_u1mix_u2mix(x1::MC{N}, x2::MC{N}, z::Interval{Float64}) where N
	xLc = z.lo
	xUc = z.hi
  	cnst = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
	cv2 = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.lo
	if (cv1 > cv2)
		cv = cv1
		cv_grad = x2.Intv.hi*x1.cv_grad + x1.Intv.hi*x2.cv_grad
	else
		cv = cv2
		cv_grad = x2.Intv.lo*x1.cc_grad + x1.Intv.lo*x2.cc_grad
	end
	cc1 = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
	cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.hi
	if (cc1 < cc2)
		cc = cc1
		cc_grad = x2.Intv.lo*x1.cv_grad + x1.Intv.hi*x2.cc_grad
	else
		cc = cc2
		cc_grad = x2.Intv.hi*x1.cc_grad + x1.Intv.lo*x2.cv_grad
	end
	cv, cc, cv_grad, cc_grad = cut(xLc, xUc, cv, cc, cv_grad, cc_grad)
	return MC{N}(cv, cc, z, cv_grad, cc_grad, cnst)
end
function mul3_u1pos_u2mix(x1::MC{N}, x2::MC{N}, z::Interval{Float64}) where N
	xLc = z.lo
	xUc = z.hi
    cnst = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
	cv2 = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
	if (cv1 > cv2)
		cv = cv1
		cv_grad = x2.Intv.hi*x1.cv_grad + x1.Intv.hi*x2.cv_grad
	else
		cv = cv2
		cv_grad = x2.Intv.lo*x1.cc_grad + x1.Intv.lo*x2.cv_grad
	end
	cc1 = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
	cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
	if (cc1 < cc2)
		cc = cc1
		cc_grad = x2.Intv.lo*x1.cv_grad + x1.Intv.hi*x2.cc_grad
	else
		cc = cc2
		cc_grad = x2.Intv.hi*x1.cc_grad + x1.Intv.lo*x2.cc_grad
	end
	cv, cc, cv_grad, cc_grad = cut(xLc, xUc, cv, cc, cv_grad, cc_grad)
	return MC{N}(cv, cc, z, cv_grad, cc_grad, cnst)
end

mul_MV_ns1cv(x1::Float64,x2::Float64,MC1::MC,MC2::MC) = MC2.Intv.hi*x1+MC1.Intv.hi*x2-MC2.Intv.hi*MC1.Intv.hi
mul_MV_ns2cv(x1::Float64,x2::Float64,MC1::MC,MC2::MC) = MC2.Intv.lo*x1+MC1.Intv.lo*x2-MC2.Intv.lo*MC1.Intv.lo
mul_MV_ns3cv(x1::Float64,x2::Float64,MC1::MC,MC2::MC) = max(mul_MV_ns1cv(x1,x2,MC1,MC2),mul_MV_ns2cv(x1,x2,MC1,MC2))
mul_MV_ns1cc(x1::Float64,x2::Float64,MC1::MC,MC2::MC) = MC2.Intv.lo*x1+MC1.Intv.hi*x2-MC2.Intv.lo*MC1.Intv.hi
mul_MV_ns2cc(x1::Float64,x2::Float64,MC1::MC,MC2::MC) = MC2.Intv.hi*x1+MC1.Intv.lo*x2-MC2.Intv.hi*MC1.Intv.lo
mul_MV_ns3cc(x1::Float64,x2::Float64,MC1::MC,MC2::MC) = min(mul_MV_ns1cc(x1,x2,MC1,MC2),mul_MV_ns2cc(x1,x2,MC1,MC2))

function multiply_MV_NS(x1::MC{N},x2::MC{N},ngrad::Int,cnst::Bool) where N

 k = diam(x2.Intv)/diam(x1.Intv)
 z = (x1.Intv.hi*x2.Intv.hi - x1.Intv.lo*x2.Intv.lo)/diam(x1.Intv)
 x1vta,blank = mid3(x1.cv,x1.cc,(x2.cv-z)/k)
 x1vtb,blank = mid3(x1.cv,x1.cc,(x2.cc-z)/k)
 x2vta,blank = mid3(x2.cv,x2.cc, k*x1.cv+z)
 x2vtb,blank = mid3(x2.cv,x2.cc, k*x1.cc+z)
 x1vt = [x1.cv, x1.cc, x1vta, x1vtb, x1.cv, x1.cc]
 x2vt = [x2vta, x2vtb, x2.cv, x2.cc, x2.cv, x2.cc]
 vt  = [mul_MV_ns3cv(x1vt[1],x2vt[1],x1,x2), mul_MV_ns3cv(x1vt[2],x2vt[2],x1,x2),
        mul_MV_ns3cv(x1vt[3],x2vt[3],x1,x2), mul_MV_ns3cv(x1vt[4],x2vt[4],x1,x2),
				mul_MV_ns3cv(x1vt[5],x2vt[5],x1,x2), mul_MV_ns3cv(x1vt[6],x2vt[6],x1,x2)]
 cv,cvind = findmax(vt)

 if (ngrad>0)
	if isequal(mul_MV_ns1cv(x1vt[cvind],x2vt[cvind],x1,x2),
						 mul_MV_ns2cv(x1vt[cvind],x2vt[cvind],x1,x2),
						 MC_param.mv_tol,MC_param.mv_tol)
 		alph = [0.0,1.0]

 		MC1thin::Bool = isequal(x1.cv,x1.cc,MC_param.mv_tol,MC_param.mv_tol)
 		if ((~MC1thin) && (x1vt[cvind] > x1.cv))
	 		if (~isequal(x1vt[cvind],x1.cv,MC_param.mv_tol,MC_param.mv_tol))
		 		alph[2] = min(alph[2],-x2.Intv.lo/diam(x2.Intv))
	 		end
 		end
 		if ((~MC1thin) && (x1vt[cvind] < x1.cc))
			if (~isequal(x1vt[cvind],x1.cc,MC_param.mv_tol,MC_param.mv_tol))
				alph[1] = max(alph[1],-x2.Intv.lo/diam(x2.Intv))
			end
 		end

 		MC2thin::Bool = isequal(x2.cv,x2.cc,MC_param.mv_tol,MC_param.mv_tol)
 		if ((~MC2thin) && (x2vt[cvind] > x2.cv))
			if (~isequal(x2vt[cvind],x2.cv,MC_param.mv_tol,MC_param.mv_tol))
				alph[2] = min(alph[2],-x1.Intv.lo/diam(x1.Intv))
			end
 		end
 		if ((~MC2thin) && (x2vt[cvind] < x2.cc))
 			if (~isequal(x2vt[cvind],x2.cc,MC_param.mv_tol,MC_param.mv_tol))
	 			alph[1] = max(alph[1],-x1.Intv.lo/diam(x1.Intv))
 			end
 		end

 		alphthin = isequal(alph[1],alph[2],MC_param.mv_tol,MC_param.mv_tol)
 		if (~alphthin && (alph[1]>alph[2]))
	 		error("Multivariant mult error alphaL = alphaU")
 		end
 		myalph = (alph[1]+alph[2])/2.0
 	elseif (mul_MV_ns1cv(x1vt[cvind],x2vt[cvind],x1,x2) >
			    mul_MV_ns2cv(x1vt[cvind],x2vt[cvind],x1,x2))
		myalph = 1.0
	else
		myalph = 0.0
	end
	sigma_cv1 = x2.Intv.lo + myalph*diam(x2.Intv)
	sigma_cv2 = x1.Intv.lo + myalph*diam(x1.Intv)
	if (x1.cnst)
		term1 = @SVector zeros(Float64,N)
	elseif (sigma_cv1>=0.0)
		term1 = x1.cv_grad
	else
		term1 = x1.cc_grad
	end
	if (x2.cnst)
		term2 = @SVector zeros(Float64,N)
	elseif (sigma_cv1 >= 0.0)
		term2 = x2.cv_grad
	else
		term2 = x2.cc_grad
	end
	cv_grad = term1*sigma_cv1 + term2*sigma_cv2
 end

 z = (x1.Intv.hi*x2.Intv.lo - x1.Intv.lo*x2.Intv.hi)/diam(x1.Intv)
 x1vta,blank = mid3(x1.cv,x1.cc,(x2.cv-z)/k)
 x1vtb,blank = mid3(x1.cv,x1.cc,(x2.cc-z)/k)
 x2vta,blank = mid3(x2.cv,x2.cc, k*x1.cv+z)
 x2vtb,blank = mid3(x2.cv,x2.cc, k*x1.cc+z)
 x1vt = [x1.cv, x1.cc, x1vta, x1vtb, x1.cv, x1.cc]
 x2vt = [x2vta, x2vtb, x2.cv, x2.cc, x2.cc, x2.cv]
 vt  = [mul_MV_ns3cc(x1vt[1],x2vt[1],x1,x2), mul_MV_ns3cc(x1vt[2],x2vt[2],x1,x2),
				mul_MV_ns3cc(x1vt[3],x2vt[3],x1,x2), mul_MV_ns3cc(x1vt[4],x2vt[4],x1,x2),
			  mul_MV_ns3cc(x1vt[5],x2vt[5],x1,x2), mul_MV_ns3cc(x1vt[6],x2vt[6],x1,x2)]
 cc,ccind = findmax(vt)

 if (ngrad>0)
 	if isequal(mul_MV_ns1cc(x1vt[cvind],x2vt[cvind],x1,x2),
						 mul_MV_ns2cc(x1vt[cvind],x2vt[cvind],x1,x2),
						 MC_param.mv_tol,MC_param.mv_tol)
		 alph = [0.0,1.0]

		 MC1thin = isequal(x1.cv,x1.cc,MC_param.mv_tol,MC_param.mv_tol)
		 if ((~MC1thin) && (x1vt[cvind] > x1.cv))
		 	if (~isequal(x1vt[cvind],x1.cv,MC_param.mv_tol,MC_param.mv_tol))
			 	alph[1] = max(alph[1],-x2.Intv.lo/diam(x2.Intv))
		 	end
		 end
		 if ((~MC1thin) && (x1vt[cvind] < x1.cc))
		 	if (~isequal(x1vt[cvind],x1.cc,MC_param.mv_tol,MC_param.mv_tol))
			 	alph[2] = min(alph[2],-x2.Intv.lo/diam(x2.Intv))
		 	end
		 end

		 MC2thin = isequal(x2.cv,x2.cc,MC_param.mv_tol,MC_param.mv_tol)
		 if ((~MC2thin) && (x2vt[cvind] > x2.cv))
		 	if (~isequal(x2vt[cvind],x2.cv,MC_param.mv_tol,MC_param.mv_tol))
			 	alph[2] = min(alph[2],x1.Intv.hi/diam(x1.Intv))
		 	end
		 end
		 if ((~MC2thin) && (x2vt[cvind] < x2.cc))
			if (~isequal(x2vt[cvind],x2.cc,MC_param.mv_tol,MC_param.mv_tol))
			 	alph[1] = max(alph[1],x1.Intv.hi/diam(x1.Intv))
			end
		 end

		 alphthin = isequal(alph[1],alph[2],MC_param.mv_tol,MC_param.mv_tol)
		 if (~alphthin && (alph[1]>alph[2]))
		 	error("Multivariant mult error alphaL = alphaU")
		 end
		 myalph = (alph[1]+alph[2])/2
	elseif (mul_MV_ns1cv(x1vt[cvind],x2vt[cvind],x1,x2) >
				  mul_MV_ns2cv(x1vt[cvind],x2vt[cvind],x1,x2))
	   myalph = 1.0
  	else
		myalph = 0.0
  	end
	sigma_cc1 = x2.Intv.lo + myalph*diam(x2.Intv)
	sigma_cc2 = x1.Intv.hi - myalph*diam(x1.Intv)
	if (x1.cnst)
		term1 = @SVector zeros(Float64,N)
	elseif (sigma_cc1>=0.0)
		term1 = x1.cc_grad
	else
		term1 = x1.cv_grad
	end
	if (x2.cnst)
		term2 =  @SVector zeros(Float64,N)
	elseif (sigma_cc1>= 0.0)
		term2 = x2.cc_grad
	else
		term2 = x2.cv_grad
	end
	cc_grad = term1*sigma_cc1 + term2*sigma_cc2
 end
 return MC{N}(cv, cc, x1.Intv*x2.Intv, cv_grad, cc_grad, cnst)
end

function multiply_STD_NS(x1::MC{N}, x2::MC{N}, y::Interval{Float64}) where N
	if (x2.Intv.lo >= 0.0)
    	(x2.cnst) && (return mul1_u1pos_u2pos(x1, x2, y, x1.cnst))
    	(x1.cnst) && (return mul1_u1pos_u2pos(x2, x1, y, x2.cnst))
  		return mul2_u1pos_u2pos(x1, x2, y)
	elseif (x2.Intv.hi <= 0.0)
		return -mult_kernel(x1, -x2, -y)
	else
    	(x2.cnst) && (return mul1_u1pos_u2mix(x1, x2, y, x1.cnst))
    	(x1.cnst) && (return mul2_u1pos_u2mix(x1, x2, y, x2.cnst))
  		return mul3_u1pos_u2mix(x1, x2, y)
	end
end

function STD_NS_ALT_kernel(x::MC{N}, y::MC{N}, z::Interval{Float64}) where N
	alpha1 = min( y.Intv.lo*x.cv,  y.Intv.lo*x.cc )
	alpha2 = min( x.Intv.lo*y.cv,  x.Intv.lo*y.cc )
	beta1  = min( y.Intv.hi*x.cv,  y.Intv.hi*x.cc )
	beta2  = min( x.Intv.hi*y.cv,  x.Intv.hi*y.cc )
	gamma1 = max( y.Intv.lo*x.cv,  y.Intv.lo*x.cc )
	gamma2 = max( x.Intv.hi*y.cv,  x.Intv.hi*y.cc )
	delta1 = max( y.Intv.hi*x.cv,  y.Intv.hi*x.cc )
	delta2 = max( x.Intv.lo*y.cv,  x.Intv.lo*y.cc )

	cv1 = alpha1 + alpha2 - x.Intv.lo*y.Intv.lo
	cv2 = beta1  + beta2  - x.Intv.hi*y.Intv.hi
	cc1 = gamma1 + gamma2 - x.Intv.hi*y.Intv.lo
	cc2 = delta1 + delta2 - x.Intv.lo*y.Intv.hi

	s_alpha1 = (y.Intv.lo >= 0.0) ? y.Intv.lo*x.cv_grad : y.Intv.lo*x.cc_grad
	s_alpha2 = (x.Intv.lo >= 0.0) ? x.Intv.lo*x.cv_grad : x.Intv.lo*x.cc_grad
	s_beta1  = (y.Intv.hi >= 0.0) ? y.Intv.hi*x.cv_grad : y.Intv.hi*x.cc_grad
	s_beta2  = (x.Intv.hi >= 0.0) ? x.Intv.hi*x.cv_grad : x.Intv.hi*x.cc_grad
	s_gamma1 = (y.Intv.lo >= 0.0) ? y.Intv.lo*x.cc_grad : y.Intv.lo*x.cv_grad
	s_gamma2 = (x.Intv.hi >= 0.0) ? x.Intv.hi*x.cc_grad : x.Intv.hi*x.cv_grad
	s_delta1 = (y.Intv.hi >= 0.0) ? y.Intv.hi*x.cc_grad : y.Intv.hi*x.cv_grad
	s_delta2 = (x.Intv.lo >= 0.0) ? x.Intv.lo*x.cc_grad : x.Intv.lo*x.cv_grad

	if (cv1 >= cv2)
	    cv = cv1
	    cv_grad = s_alpha1 + s_alpha2
	else
	    cv = cv2
	    cv_grad = s_beta1 + s_beta2
	end

	if (cc1 <= cc2)
	    cc = cc1
	    cc_grad = s_gamma1 + s_gamma2
	else
	    cc = cc2
	    cc_grad = s_delta1 + s_delta2
	end
	return MC{N}(cv, cc, z, cv_grad, cc_grad, (x.cnst && y.cnst))
end

function mult_kernel(x1::MC, x2::MC, y::Interval{Float64})

	degen1 = ((x1.Intv.hi - x1.Intv.lo) == 0.0)
	degen2 = ((x2.Intv.hi - x2.Intv.lo) == 0.0)
	if (MC_param.mu >= 1 && ~(degen1||degen2))
		return multiply_MV(x1, x2, y)
		#=
	elseif (MC_param.multivar_refine && ~(degen1||degen2))
		#println("NS MV mult trace 1")
		if (x2.cnst)
			cnst = x1.cnst
		elseif (x1.cnst)
			cnst = x2.cnst
		elseif (length(x1.cc_grad) != length(x2.cc_grad))
			error("Unequal subgradients")
		else
			cnst = (x1.cnst||x2.cnst)
		end
		return multiply_MV_NS(x1,x2,N,cnst) # DONE (minus gradients & case handling?)
		=#
	elseif (x1.Intv.lo >= 0.0)
		return multiply_STD_NS(x1, x2, y)
		#return STD_NS_ALT(x1,x2)
	elseif (x1.Intv.hi <= 0.0)
		(x2.Intv.lo >= 0.0) && (return -mult_kernel(-x1, x2, -y))
	    (x2.Intv.hi <= 0.0) && (return mult_kernel(-x1, -x2, y))
		return -mult_kernel(x2, -x1, -y)
	elseif (x2.Intv.lo >= 0.0)
		return mult_kernel(x2, x1, y)
	elseif (x2.Intv.hi <= 0.0)
		return -mult_kernel(-x2, x1, -y)
	else
    	(x2.cnst) && (return STD_NS_ALT_kernel(x1, x2, y)) #return mul1_u1mix_u2mix(x1,x2,x1.cnst)
    	(x1.cnst) && (return STD_NS_ALT_kernel(x1, x2, y)) #return mul1_u1mix_u2mix(x2,x1,x2.cnst)
  		return mul2_u1mix_u2mix(x1, x2, y) #return STD_NS_ALT(x1,x2)
	end
end
function *(x1::MC, x2::MC)
	degen1 = ((x1.Intv.hi - x1.Intv.lo) == 0.0)
	degen2 = ((x2.Intv.hi - x2.Intv.lo) == 0.0)
	if (MC_param.mu >= 1 && ~(degen1||degen2))
		if (min(x1.Intv.lo, x2.Intv.lo) < 0.0 < max(x1.Intv.hi, x2.Intv.hi))
			lo_Intv_calc::Float64 = gCxAIntv(x1.Intv, x2.Intv, x1.Intv, x2.Intv, x1, x2)
			hi_Intv_calc::Float64 = -gCxAIntv(-x1.Intv, x2.Intv, -x1.Intv, x2.Intv, x1, x2)
			z = Interval{Float64}(lo_Intv_calc, hi_Intv_calc)
		else
			z = x1.Intv*x2.Intv
		end
	else
		z = x1.Intv*x2.Intv
	end
	return mult_kernel(x1, x2, z)
end
