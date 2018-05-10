########### Defines multiplicaiton of SMC and SMC
function sigu(x::T,mu1T::T) where {T<:AbstractFloat}
	 (zero(T)<=x) ? x^(one(T)/mu1T) : -abs(x)^(one(T)/mu1T)
 end

function gCxAcv(alpha::V,beta::V,lambda::V,nu::V,x1::SMCg{N,V,T},x2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
		# gCxA pre-terms
		alplo::T = alpha.lo
		alphi::T = alpha.hi
		betlo::T = beta.lo
		bethi::T = beta.hi
		LmdDel::T = lambda.hi-lambda.lo
		NuDel::T = nu.hi-nu.lo
		LmdSum::T = lambda.lo+lambda.hi
		NuSum::T = nu.lo+nu.hi
		muT::T = convert(T,MC_param.mu)
		mu1::Int64 = MC_param.mu+1
		mu1T::T = convert(T,MC_param.mu+1)
		mu1n::Int64 = MC_param.mu-1

		xslo::T = lambda.lo+LmdDel*(((nu.hi-betlo)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT)) # correct to minus?, correct internal minus? exponent to mu
		xshi::T = lambda.lo+LmdDel*(((nu.hi-bethi)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT)) # correct to minus?
		yslo::T = nu.lo+NuDel*(((lambda.hi-alplo)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))
		yshi::T = nu.lo+NuDel*(((lambda.hi-alphi)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))

		# calculates term 1
		if (xslo <= alplo)
			term1::T = alplo
		elseif (alphi <= xslo)
			term1 = alphi
		else
			term1 = xslo
		end

		# calculates term 2
		if (xshi <= alplo)
			term2::T = alplo
		elseif (alphi <= xshi)
			term2 = alphi
		else
			term2 = xshi
		end

		# calculates term 3
		if (yslo <= betlo)
			term3::T = betlo
		elseif (bethi <= yslo)
			term3 = bethi
		else
			term3 = yslo
		end

		# calculates term 4
		if (yshi <= betlo)
			term4::T = betlo
		elseif (bethi <= yshi)
			term4 = bethi
		else
			term4 = yshi
		end

		# calculates the convex relaxation
		tempGxA1a::T = LmdDel*NuDel*abs((betlo-nu.lo)/NuDel-(lambda.hi-term1)/LmdDel)^mu1
		tempGxA1::T = half(T)*(term1*NuSum+betlo*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA1a)
		tempGxA2a::T = LmdDel*NuDel*abs((bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel)^mu1
		tempGxA2::T = half(T)*(term2*NuSum+bethi*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA2a)
		tempGxA3a::T = LmdDel*NuDel*abs((term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel)^mu1
		tempGxA3::T = half(T)*(alplo*NuSum+term3*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA3a)
		tempGxA4a::T = LmdDel*NuDel*abs((term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel)^mu1
		tempGxA4::T = half(T)*(alphi*NuSum+term4*LmdSum-(lambda.lo*nu.lo+lambda.hi*nu.hi)+tempGxA4a)

		# gets minima which is cv/cc/Intv depending on arguments
		a = min(tempGxA1,tempGxA2,tempGxA3,tempGxA4)

		if (a == tempGxA1)
				psi_mul::T = (betlo-nu.lo)/NuDel - (lambda.hi-term1)/LmdDel
				psi_mult_x1::T = half(T)*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1::T = half(T)*(LmdSum+mu1T*(LmdDel)*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA2)
				psi_mul = (bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel
				psi_mult_x1 = half(T)*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1 = half(T)*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA3)
				psi_mul = (term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel
				psi_mult_x1 = half(T)*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1 = half(T)*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		else
				psi_mul = (term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel
				psi_mult_x1 = half(T)*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
				psi_mult_y1 = half(T)*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		end

    	grad::SVector{N,T} = max(zero(T),psi_mult_x1)*x1.cv_grad+
							 min(zero(T),psi_mult_x1)*x1.cc_grad+
							 max(zero(T),psi_mult_y1)*x2.cv_grad+
							 min(zero(T),psi_mult_y1)*x2.cc_grad

    return a,grad
end

function gCxAcc(alpha::V,beta::V,lambda::V,nu::V,x1::SMCg{N,V,T},x2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
		# gCxA pre-terms
		alplo::T = alpha.lo
		alphi::T = alpha.hi
		betlo::T = beta.lo
		bethi::T = beta.hi
		LmdDel::T = lambda.hi-lambda.lo
		NuDel::T = nu.hi-nu.lo
		LmdSum::T = lambda.lo+lambda.hi
		NuSum::T = nu.lo+nu.hi
		NuDotLmd::T = lambda.lo*nu.lo+lambda.hi*nu.hi
		muT::T = convert(T,MC_param.mu)
		mu1::Int64 = MC_param.mu+1
		mu1n::Int64 = MC_param.mu-1
		mu1T::T = convert(T,MC_param.mu+1)

		xslo::T = lambda.lo+LmdDel*(((nu.hi-betlo)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT))
		xshi::T = lambda.lo+LmdDel*(((nu.hi-bethi)/NuDel)-sigu((NuSum)/(mu1T*(NuDel)),muT))
		yslo::T = nu.lo+NuDel*(((lambda.hi-alplo)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))
		yshi::T = nu.lo+NuDel*(((lambda.hi-alphi)/LmdDel)-sigu((LmdSum)/(mu1T*LmdDel),muT))

		# calculates term 1
		if (xslo <= alplo)
			term1::T = alplo
		elseif (alphi <= xslo)
			term1 = alphi
		else
			term1 = xslo
		end

		# calculates term 2
		if (xshi <= alplo)
			term2::T = alplo
		elseif (alphi <= xshi)
			term2 = alphi
		else
			term2 = xshi
		end

		# calculates term 3
		if (yslo <= betlo)
			term3::T = betlo
		elseif (bethi <= yslo)
			term3 = bethi
		else
			term3 = yslo
		end

		# calculates term 4
		if (yshi <= betlo)
			term4::T = betlo
		elseif (bethi <= yshi)
			term4 = bethi
		else
			term4 = yshi
		end

		# calculates the convex relaxation
		tempGxA1a::T = LmdDel*NuDel*abs((betlo-nu.lo)/NuDel-(lambda.hi-term1)/LmdDel)^mu1
		tempGxA1::T = half(T)*(term1*NuSum+betlo*LmdSum-NuDotLmd+tempGxA1a)
		tempGxA2a::T = LmdDel*NuDel*abs((bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel)^mu1
		tempGxA2::T = half(T)*(term2*NuSum+bethi*LmdSum-NuDotLmd+tempGxA2a)
		tempGxA3a::T = LmdDel*NuDel*abs((term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel)^mu1
		tempGxA3::T = half(T)*(alplo*NuSum+term3*LmdSum-NuDotLmd+tempGxA3a)
		tempGxA4a::T = LmdDel*NuDel*abs((term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel)^mu1
		tempGxA4::T = half(T)*(alphi*NuSum+term4*LmdSum-NuDotLmd+tempGxA4a)

		# gets minima which is cv/cc/Intv depending on arguments
		a::T = min(tempGxA1,tempGxA2,tempGxA3,tempGxA4)

	   if (a == tempGxA1)
		   psi_mul::T = (betlo-nu.lo)/NuDel - (lambda.hi-term1)/LmdDel
		   psi_mult_x1::T = half(T)*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
		   psi_mult_y1::T = half(T)*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA2)
			psi_mul = (bethi-nu.lo)/NuDel - (lambda.hi-term2)/LmdDel
			psi_mult_x1 = half(T)*(NuSum+mu1*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
			psi_mult_y1 = half(T)*(LmdSum+mu1*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		elseif (a == tempGxA3)
			psi_mul = (term3-nu.lo)/NuDel - (lambda.hi-alplo)/LmdDel
			psi_mult_x1 = half(T)*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
			psi_mult_y1 = half(T)*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		else
			psi_mul= (term4-nu.lo)/NuDel - (lambda.hi-alphi)/LmdDel
			psi_mult_x1 = half(T)*(NuSum+mu1T*NuDel*(psi_mul)*abs(psi_mul)^mu1n)
			psi_mult_y1 = half(T)*(LmdSum+mu1T*LmdDel*(psi_mul)*abs(psi_mul)^mu1n)
		end
		grad2::SVector{N,T} = max(zero(T),psi_mult_x1)*x1.cc_grad+
							  min(zero(T),psi_mult_x1)*x1.cv_grad-
							  max(zero(T),psi_mult_y1)*x2.cc_grad-
							  min(zero(T),psi_mult_y1)*x2.cv_grad
    return a,grad2
end

function gCxAIntv(alpha::V,beta::V,lambda::V,nu::V,x1::SMCg{N,V,T},x2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
		# gCxA pre-terms
		alplo::T = alpha.lo
		alphi::T = alpha.hi
		betlo::T = beta.lo
		bethi::T = beta.hi
		LmdDel::T = lambda.hi-lambda.lo
		NuDel::T = nu.hi-nu.lo
		LmdSum::T = lambda.lo+lambda.hi
		NuSum::T = nu.lo+nu.hi
		NuLmd::T = lambda.lo*nu.lo+lambda.hi*nu.hi
		mu1::Int64 = MC_param.mu+1
		mu1n::Int64 = MC_param.mu-1
		mu1T::T = convert(T,MC_param.mu+1)

		xslo::T = lambda.lo+LmdDel*(((nu.hi-betlo)/NuDel)+sigu(-(NuSum)/(mu1*(NuDel)),mu1T))
		xshi::T = lambda.lo+LmdDel*(((nu.hi-bethi)/NuDel)+sigu(-(NuSum)/(mu1*(NuDel)),mu1T))
		yslo::T = nu.lo+NuDel*(((lambda.hi-alplo)/LmdDel)-sigu((LmdSum)/(mu1*LmdDel),mu1T))
		yshi::T = nu.lo+NuDel*(((lambda.hi-alphi)/LmdDel)-sigu((LmdSum)/(mu1*LmdDel),mu1T))

		# calculates term 1
		if (xslo <= alplo)
			term1::T = alplo
		elseif (alphi <= xslo)
			term1 = alphi
		else
			term1 = xslo
		end

		# calculates term 2
		if (xshi <= alplo)
			term2::T = alplo
		elseif (alphi <= xshi)
			term2 = alphi
		else
			term2 = xshi
		end

		# calculates term 3
		if (yslo <= betlo)
			term3::T = betlo
		elseif (bethi <= yslo)
			term3 = bethi
		else
			term3 = yslo
		end

		# calculates term 4
		if (yshi <= betlo)
			term4::T = betlo
		elseif (bethi <= yshi)
			term4 = bethi
		else
			term4 = yshi
		end

		# calculates the convex relaxation
		tempGxA1a::T = LmdDel*NuDel*abs((betlo-nu.lo)/NuDel-(lambda.hi-term1)/LmdDel)^mu1
		tempGxA1::T =  half(T)*(term1*NuSum+betlo*LmdSum-NuLmd+tempGxA1a)
		tempGxA2a::T = LmdDel*NuDel*abs((bethi-nu.lo)/NuDel-(lambda.hi-term2)/LmdDel)^mu1
		tempGxA2::T =  half(T)*(term2*NuSum+bethi*LmdSum-NuLmd+tempGxA2a)
		tempGxA3a::T = LmdDel*NuDel*abs((term3-nu.lo)/NuDel-(lambda.hi-alplo)/LmdDel)^mu1
		tempGxA3::T =  half(T)*(alplo*NuSum+term3*LmdSum-NuLmd+tempGxA3a)
		tempGxA4a::T = LmdDel*NuDel*abs((term4-nu.lo)/NuDel-(lambda.hi-alphi)/LmdDel)^mu1
		tempGxA4::T = half(T)*(alphi*NuSum+term4*LmdSum-NuLmd+tempGxA4a)

		# gets minima which is cv/cc/Intv depending on arguments
		return min(tempGxA1,tempGxA2,tempGxA3,tempGxA4)
end

function multiply_MV(x1::SMCg{N,V,T},x2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	x1cv::T = x1.cv
	x2cv::T = x2.cv
	x1cc::T = x1.cc
	x2cc::T = x2.cc
	x1lo::T = x1.Intv.lo
	x2lo::T = x2.Intv.lo
	x1hi::T = x1.Intv.hi
	x2hi::T = x2.Intv.hi
	x1h_l::T = x1hi-x1lo
	x1h_v::T = x1hi-x1cv
	x2h_l::T = x2hi-x2lo
	x2v_l::T = x2cv-x2lo
	cv::T = zero(T)
	cc::T = zero(T)
	mu1::Int64 = MC_param.mu+1
	muf1::T = convert(T,MC_param.mu+1)
	alpha0::V = V(x1cv,x1cc)
	beta0::V = V(x2cv,x2cc)
	if (zero(T)<=x1lo) && (zero(T)<=x2lo)
		b1_temp1::T = max(zero(T),((x2v_l/x2h_l)-(x1h_v/x1h_l)))
		b1_term2::T = x2lo+muf1*x2h_l*b1_temp1^MC_param.mu
		b1_term3::T = x1lo+muf1*x1h_l*b1_temp1^MC_param.mu
		cv = x1cv*x2lo + x1lo*x2v_l + x1h_l*x2h_l*b1_temp1^mu1
		cv_grad::SVector{N,T} = b1_term2*x1.cv_grad + b1_term3*x2.cv_grad
	elseif ((x1hi<=zero(T))) && (x2hi<=zero(T))
		b1_temp = max(zero(T),(x2hi - x2cc)/x2h_l - (x1cc-x1lo)/(x1h_l))
		cv = x1cc*x2hi + x1hi*(x2cc - x2hi) + (x1h_l)*(x2h_l)*b1_temp^mu1
		b1_term2 = -x2hi+muf1*(x2h_l)*max(zero(T),b1_temp)^MC_param.mu
		b1_term3 = -x1hi+muf1*(x1h_l)*max(zero(T),b1_temp)^MC_param.mu
		cv_grad = (-b1_term2)*x1.cc_grad - b1_term3*x2.cc_grad
	else
		cv,cv_grad = gCxAcv(alpha0,beta0,x1.Intv,x2.Intv,x1,x2)
	end
	if ((x1hi<=zero(T))) && ((zero(T))<=x2lo)
		btemp1 = max(zero(T),x2v_l/x2h_l - (x1cc-x1lo)/(x1h_l))
		btemp2 = x2lo+muf1*x2h_l*btemp1^MC_param.mu
		btemp3 = -x1hi+muf1*x1h_l*btemp1^MC_param.mu
		cc = x1cc*x2lo-x1hi*(x2lo-x2cv)-x1h_l*x2h_l*btemp1^mu1
		cc_grad::SVector{N,T} = btemp2*x1.cc_grad - btemp3*x2.cv_grad
	elseif ((zero(T))<=x1.Intv.lo) && (x2.Intv.hi<=zero(T))
		btemp1 = max(zero(T),(x2hi-x2cc)/x2h_l - x1h_v/x1h_l)
		btemp2 = -x2hi+muf1*(x2h_l)*btemp1^MC_param.mu
		btemp3 = x1lo+muf1*x1h_l*btemp1^MC_param.mu
		cc = x1cv*x2hi+x1lo*(x2cc-x2hi)-x1h_l*x2h_l*btemp1^mu1
		cc_grad = -btemp2*x1.cv_grad + btemp3*x2.cc_grad
	else
		cct::T,cc_gradt::SVector{N,T} = gCxAcc(-alpha0,beta0,-x1.Intv,x2.Intv,x1,x2)
		cc = -cct
		cc_grad = cc_gradt
	end
	if (min(x1lo,x2lo)<zero(T)<max(x1hi,x2hi))
		lo_Intv_calc::T = gCxAIntv(x1.Intv,x2.Intv,x1.Intv,x2.Intv,x1,x2)
		hi_Intv_calc::T = -gCxAIntv(-x1.Intv,x2.Intv,-x1.Intv,x2.Intv,x1,x2)
		Intv_calc::V = V(lo_Intv_calc,hi_Intv_calc)
	else
		Intv_calc = x1.Intv*x2.Intv
	end
	cnst = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, Intv_calc, cnst, x1.IntvBox,x1.xref)
end

function mul1_u1pos_u2pos(x1::SMCg{N,V,T},x2::SMCg{N,V,T},cnst::Bool) where {N,V,T<:AbstractFloat}
  Intv::V = x1.Intv*x2.Intv
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  cv1::T = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2::T = x2.Intv.lo*x1.cv + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
  if (cv1 > cv2)
    cv::T = cv1
    cv_grad::SVector{N,T} = x2.Intv.hi*x1.cv_grad
  else
    cv = cv2
    cv_grad = x2.Intv.lo*x1.cv_grad
  end

  cc1::T = x2.Intv.lo*x1.cc + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2::T = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
  if (cc1 < cc2)
    cc::T = cc1
    cc_grad::SVector{N,T} = x2.Intv.lo*x1.cc_grad
  else
    cc = cc2
    cc_grad = x2.Intv.hi*x1.cc_grad
  end
  cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  return SMCg{N,V,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end
function mul1_u1pos_u2mix(x1::SMCg{N,V,T},x2::SMCg{N,V,T},cnst::Bool) where {N,V,T<:AbstractFloat}
  Intv::V = x1.Intv*x2.Intv
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  cv1::T = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2::T = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
  if (cv1 > cv2)
    cv::T = cv1
    cv_grad::SVector{N,T} = x2.Intv.hi*x1.cv_grad
  else
    cv = cv2
    cv_grad = x2.Intv.lo*x1.cc_grad
  end

  cc1::T = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2::T = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
  if (cc1 < cc2)
    cc::T = cc1
    cc_grad::SVector{N,T} = x2.Intv.lo*x1.cv_grad
  else
    cc = cc2
    cc_grad = x2.Intv.hi*x1.cc_grad
  end
  cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  return SMCg{N,V,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end
function mul1_u1mix_u2mix(x1::SMCg{N,V,T},x2::SMCg{N,V,T},cnst::Bool) where {N,V,T<:AbstractFloat}
  Intv::V = x1.Intv*x2.Intv
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  cv1::T = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2::T = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.lo
  cv::T = cv1
  if (cv1 > cv2)
    cv = cv1
    cv_grad::SVector{N,T} = x2.Intv.hi*x1.cv_grad
  else
    cv = cv2
    cv_grad = x2.Intv.lo*x1.cc_grad
  end
  cc1::T = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2::T = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.hi
  cc::T = cc1
  if (cc1 < cc2)
    cc = cc1
    cc_grad::SVector{N,T} = x2.Intv.lo*x1.cv_grad
  else
    cc = cc2
    cc_grad = x2.Intv.hi*x1.cc_grad
  end
  cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  return SMCg{N,V,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end
function mul2_u1pos_u2pos(x1::SMCg{N,V,T},x2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	Intv::V = x1.Intv*x2.Intv
	xLc::T = Intv.lo
	xUc::T = Intv.hi
	cnst::Bool = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	cv1::T = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
	cv2::T = x2.Intv.lo*x1.cv + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
	if (cv1 > cv2)
		cv::T = cv1
		cv_grad::SVector{N,T} = x2.Intv.hi*x1.cv_grad + x1.Intv.hi*x2.cv_grad
	else
		cv = cv2
		cv_grad = x2.Intv.lo*x1.cv_grad + x1.Intv.lo*x2.cv_grad
	end

	cc1::T = x2.Intv.lo*x1.cc + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
	cc2::T = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
	if (cc1 < cc2)
		cc::T = cc1
		cc_grad::SVector{N,T} = x2.Intv.lo*x1.cc_grad + x1.Intv.hi*x2.cc_grad
	else
		cc = cc2
		cc_grad = x2.Intv.hi*x1.cc_grad + x1.Intv.lo*x2.cc_grad
	end
	cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
	return SMCg{N,V,T}(cc,cv,cc_grad,cv_grad,Intv, cnst,x1.IntvBox,x1.xref)
end
function mul2_u1pos_u2mix(x1::SMCg{N,V,T},x2::SMCg{N,V,T},cnst::Bool) where {N,V,T<:AbstractFloat}
  Intv::V = x1.Intv*x2.Intv
  xLc::T = Intv.lo
  xUc::T = Intv.hi
  cv1::T = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2::T = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
  if (cv1 > cv2)
    cv::T = cv1
    cv_grad::SVector{N,T} = x1.Intv.hi*x2.cv_grad
  else
    cv = cv2
    cv_grad = x1.Intv.lo*x2.cc_grad
  end

  cc1::T = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2::T = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
  if (cc1 < cc2)
    cc::T = cc1
    cc_grad::SVector{N,T} = x1.Intv.hi*x2.cc_grad
  else
    cc = cc2
    cc_grad = x1.Intv.lo*x2.cc_grad
  end
  cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  return SMCg{N,V,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end
function mul2_u1mix_u2mix(x1::SMCg{N,V,T},x2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	Intv::V = x1.Intv*x2.Intv
	xLc::T = Intv.lo
	xUc::T = Intv.hi
  	cnst::Bool = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	cv1::T = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
	cv2::T = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.lo
	if (cv1 > cv2)
		cv::T = cv1
		cv_grad::SVector{N,T} = x2.Intv.hi*x1.cv_grad + x1.Intv.hi*x2.cv_grad
	else
		cv = cv2
		cv_grad = x2.Intv.lo*x1.cc_grad + x1.Intv.lo*x2.cc_grad
	end

	cc1::T = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
	cc2::T = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.hi
	if (cc1 < cc2)
		cc::T = cc1
		cc_grad::SVector{N,T} = x2.Intv.lo*x1.cv_grad + x1.Intv.hi*x2.cc_grad
	else
		cc = cc2
		cc_grad = x2.Intv.hi*x1.cc_grad + x1.Intv.lo*x2.cv_grad
	end
	cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
	return SMCg{N,V,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end
function mul3_u1pos_u2mix(x1::SMCg{N,V,T},x2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	Intv::V = x1.Intv*x2.Intv
	xLc::T = Intv.lo
	xUc::T = Intv.hi
    cnst::Bool = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	cv1::T = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
	cv2::T = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
	if (cv1 > cv2)
		cv::T = cv1
		cv_grad::SVector{N,T} = x2.Intv.hi*x1.cv_grad + x1.Intv.hi*x2.cv_grad
	else
		cv = cv2
		cv_grad = x2.Intv.lo*x1.cc_grad + x1.Intv.lo*x2.cv_grad
	end

	cc1::T = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
	cc2::T = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
	if (cc1 < cc2)
		cc::T = cc1
		cc_grad::SVector{N,T} = x2.Intv.lo*x1.cv_grad + x1.Intv.hi*x2.cc_grad
	else
		cc = cc2
		cc_grad = x2.Intv.hi*x1.cc_grad + x1.Intv.lo*x2.cc_grad
	end
	cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
	return SMCg{N,V,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end

function mul_MV_ns1cv(x1::T,x2::T,MC1::SMCg{N,V,T},MC2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	return MC2.Intv.hi*x1+MC1.Intv.hi*x2-MC2.Intv.hi*MC1.Intv.hi
end
function mul_MV_ns2cv(x1::T,x2::T,MC1::SMCg{N,V,T},MC2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	return MC2.Intv.lo*x1+MC1.Intv.lo*x2-MC2.Intv.lo*MC1.Intv.lo
end
function mul_MV_ns3cv(x1::T,x2::T,MC1::SMCg{N,V,T},MC2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	return max(mul_MV_ns1cv(x1,x2,MC1,MC2),mul_MV_ns2cv(x1,x2,MC1,MC2))
end
function mul_MV_ns1cc(x1::T,x2::T,MC1::SMCg{N,V,T},MC2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	return MC2.Intv.lo*x1+MC1.Intv.hi*x2-MC2.Intv.lo*MC1.Intv.hi
end
function mul_MV_ns2cc(x1::T,x2::T,MC1::SMCg{N,V,T},MC2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	return MC2.Intv.hi*x1+MC1.Intv.lo*x2-MC2.Intv.hi*MC1.Intv.lo
end
function mul_MV_ns3cc(x1::T,x2::T,MC1::SMCg{N,V,T},MC2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	return min(mul_MV_ns1cc(x1,x2,MC1,MC2),mul_MV_ns2cc(x1,x2,MC1,MC2))
end
function multiply_MV_NS(x1::SMCg{N,V,T},x2::SMCg{N,V,T},ngrad::Int64,cnst::Bool) where {N,V,T<:AbstractFloat}

 k::T = diam(x2.Intv)/diam(x1.Intv)
 z::T = (x1.Intv.hi*x2.Intv.hi - x1.Intv.lo*x2.Intv.lo)/diam(x1.Intv)
 x1vta::T,blank::Int64 = mid3(x1.cv,x1.cc,(x2.cv-z)/k)
 x1vtb::T,blank = mid3(x1.cv,x1.cc,(x2.cc-z)/k)
 x2vta::T,blank = mid3(x2.cv,x2.cc, k*x1.cv+z)
 x2vtb::T,blank = mid3(x2.cv,x2.cc, k*x1.cc+z)
 x1vt = [x1.cv, x1.cc, x1vta, x1vtb, x1.cv, x1.cc]
 x2vt = [x2vta, x2vtb, x2.cv, x2.cc, x2.cv, x2.cc]
 vt  = [mul_MV_ns3cv(x1vt[1],x2vt[1],x1,x2), mul_MV_ns3cv(x1vt[2],x2vt[2],x1,x2),
        mul_MV_ns3cv(x1vt[3],x2vt[3],x1,x2), mul_MV_ns3cv(x1vt[4],x2vt[4],x1,x2),
				mul_MV_ns3cv(x1vt[5],x2vt[5],x1,x2), mul_MV_ns3cv(x1vt[6],x2vt[6],x1,x2)]
 cv::T,cvind::Int64 = findmax(vt)

 if (ngrad>0)
	if isequal(mul_MV_ns1cv(x1vt[cvind],x2vt[cvind],x1,x2),
						 mul_MV_ns2cv(x1vt[cvind],x2vt[cvind],x1,x2),
						 MC_param.mv_tol,MC_param.mv_tol)
 		alph = [zero(T),one(T)]

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

 		alphthin::Bool = isequal(alph[1],alph[2],MC_param.mv_tol,MC_param.mv_tol)
 		if (~alphthin && (alph[1]>alph[2]))
	 		error("Multivariant mult error alphaL = alphaU")
 		end
 		myalph::T = (alph[1]+alph[2])/2
 	elseif (mul_MV_ns1cv(x1vt[cvind],x2vt[cvind],x1,x2) >
			    mul_MV_ns2cv(x1vt[cvind],x2vt[cvind],x1,x2))
		myalph = one(T)
	else
		myalph = zero(T)
	end
	sigma_cv1::T = x2.Intv.lo + myalph*diam(x2.Intv)
	sigma_cv2::T = x1.Intv.lo + myalph*diam(x1.Intv)
	if (x1.cnst)
		term1::SVector{N,T} = @SVector zeros(T,N)
	elseif (sigma_cv1>=zero(T))
		term1 = x1.cv_grad
	else
		term1 = x1.cc_grad
	end
	if (x2.cnst)
		term2::SVector{N,T} = @SVector zeros(T,N)
	elseif (sigma_cv1>=zero(T))
		term2 = x2.cv_grad
	else
		term2 = x2.cc_grad
	end
	cv_grad::SVector{N,T} = term1*sigma_cv1 + term2*sigma_cv2
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
 cc::T,ccind::Int64 = findmax(vt)

 if (ngrad>0)
 	if isequal(mul_MV_ns1cc(x1vt[cvind],x2vt[cvind],x1,x2),
						 mul_MV_ns2cc(x1vt[cvind],x2vt[cvind],x1,x2),
						 MC_param.mv_tol,MC_param.mv_tol)
		 alph = [zero(T),one(T)]

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
	   myalph = one(T)
  	else
		myalph = zero(T)
  	end
	sigma_cc1::T = x2.Intv.lo + myalph*diam(x2.Intv)
	sigma_cc2::T = x1.Intv.hi - myalph*diam(x1.Intv)
	if (x1.cnst)
		term1 = @SVector zeros(T,N)
	elseif (sigma_cc1>=zero(T))
		term1 = x1.cc_grad
	else
		term1 = x1.cv_grad
	end
	if (x2.cnst)
		term2 =  @SVector zeros(T,N)
	elseif (sigma_cc1>=zero(T))
		term2 = x2.cc_grad
	else
		term2 = x2.cv_grad
	end
	cc_grad = term1*sigma_cc1 + term2*sigma_cc2
 end
 return SMCg{N,V,T}(cc,cv,cc_grad,cv_grad,x1.Intv*x2.Intv,cnst,x1.IntvBox,x1.xref)
end

function multiply_STD_NS(x1::SMCg{N,V,T},x2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	if (x2.Intv.lo >= zero(T))
    	if (x2.cnst)
      		return mul1_u1pos_u2pos(x1,x2,x1.cnst) # cnst to do
    	elseif (x1.cnst)
      		return mul1_u1pos_u2pos(x2,x1,x2.cnst) # cnst to do
    	else
      		return mul2_u1pos_u2pos(x1,x2) # DONE
    	end
	elseif (x2.Intv.hi <= zero(T))
		return -(x1*(-x2))
	else
    	if (x2.cnst)
      		return mul1_u1pos_u2mix(x1,x2,x1.cnst) # cnst to do
    	elseif (x1.cnst)
      		return mul2_u1pos_u2mix(x1,x2,x2.cnst) # cnst to do
    	else
	  		return mul3_u1pos_u2mix(x1,x2) # cnst to do
    	end
	end
end

function STD_NS_ALT(x::SMCg{N,V,T},y::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	alpha1::T = min( y.Intv.lo*x.cv,  y.Intv.lo*x.cc )
	alpha2::T = min( x.Intv.lo*y.cv,  x.Intv.lo*y.cc )
	beta1::T  = min( y.Intv.hi*x.cv,  y.Intv.hi*x.cc )
	beta2::T  = min( x.Intv.hi*y.cv,  x.Intv.hi*y.cc )
	gamma1::T = max( y.Intv.lo*x.cv,  y.Intv.lo*x.cc )
	gamma2::T = max( x.Intv.hi*y.cv,  x.Intv.hi*y.cc )
	delta1::T = max( y.Intv.hi*x.cv,  y.Intv.hi*x.cc )
	delta2::T = max( x.Intv.lo*y.cv,  x.Intv.lo*y.cc )

	cv1::T = alpha1 + alpha2 - x.Intv.lo*y.Intv.lo
	cv2::T = beta1  + beta2  - x.Intv.hi*y.Intv.hi
	cc1::T = gamma1 + gamma2 - x.Intv.hi*y.Intv.lo
	cc2::T = delta1 + delta2 - x.Intv.lo*y.Intv.hi

	s_alpha1::SVector{N,T} = (y.Intv.lo >= zero(T)) ? y.Intv.lo*x.cv_grad : y.Intv.lo*x.cc_grad
	s_alpha2::SVector{N,T} = (x.Intv.lo >= zero(T)) ? x.Intv.lo*x.cv_grad : x.Intv.lo*x.cc_grad
	s_beta1::SVector{N,T}  = (y.Intv.hi >= zero(T)) ? y.Intv.hi*x.cv_grad : y.Intv.hi*x.cc_grad
	s_beta2::SVector{N,T}  = (x.Intv.hi >= zero(T)) ? x.Intv.hi*x.cv_grad : x.Intv.hi*x.cc_grad
	s_gamma1::SVector{N,T} = (y.Intv.lo >= zero(T)) ? y.Intv.lo*x.cc_grad : y.Intv.lo*x.cv_grad
	s_gamma2::SVector{N,T} = (x.Intv.hi >= zero(T)) ? x.Intv.hi*x.cc_grad : x.Intv.hi*x.cv_grad
	s_delta1::SVector{N,T} = (y.Intv.hi >= zero(T)) ? y.Intv.hi*x.cc_grad : y.Intv.hi*x.cv_grad
	s_delta2::SVector{N,T} = (x.Intv.lo >= zero(T)) ? x.Intv.lo*x.cc_grad : x.Intv.lo*x.cv_grad

	if (cv1 >= cv2)
	    cv::T = cv1
	    cv_grad::SVector{N,T} = s_alpha1 + s_alpha2
	else
	    cv = cv2
	    cv_grad = s_beta1 + s_beta2
	end

	if (cc1 <= cc2)
	    cc::T = cc1
	    cc_grad::SVector{N,T} = s_gamma1 + s_gamma2
	else
	    cc = cc2
	    cc_grad = s_delta1 + s_delta2
	end
	Intv::V = x.Intv*y.Intv
	return SMCg{N,V,T}(cc,cv,cc_grad,cv_grad,Intv,(x.cnst && y.cnst),x.IntvBox,x.xref)
end

function *(x1::SMCg{N,V,T},x2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	if x1 == x2
		#println("sqr trace")
		return sqr(x1)
	end

	degen1::Bool = ((x1.Intv.hi - x1.Intv.lo) == zero(T))
	degen2::Bool = ((x2.Intv.hi - x2.Intv.lo) == zero(T))

	if (MC_param.mu >= 1 && ~(degen1||degen2))
		#println("traces multi")
		return multiply_MV(x1,x2)
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
	elseif (x1.Intv.lo >= zero(T))
		return multiply_STD_NS(x1,x2)
		#return STD_NS_ALT(x1,x2)
	elseif (x1.Intv.hi <= zero(T))
		if (x2.Intv.lo >= zero(T))
			return -((-x1)*x2)
		elseif (x2.Intv.hi <= zero(T))
			return (-x1)*(-x2)
		else
			return -(x2*(-x1))
		end
	elseif (x2.Intv.lo >= zero(T))
		return x2*x1
	elseif (x2.Intv.hi <= zero(T))
		return -((-x2)*x1)
	else
    	if (x2.cnst)
			#println("NS mult trace 5")
	  		return STD_NS_ALT(x1,x2)
      		#return mul1_u1mix_u2mix(x1,x2,x1.cnst)
    	elseif (x1.cnst)
	  		#println("NS mult trace 6")
	  		return STD_NS_ALT(x1,x2)
      		#return mul1_u1mix_u2mix(x2,x1,x2.cnst)
    	else
	  		#return STD_NS_ALT(x1,x2)
	  		return mul2_u1mix_u2mix(x1,x2)
    	end
	end
end
