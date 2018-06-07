#=
function cv_max(x::T,xL::T,xU::T,ca::T) where {T<:AbstractFloat}
        if (xU<=ca)
          return ca, zero(T)
        elseif (ca<=xL)
          return x, one(T)
        else
          term::T = max(zero(T),(x-ca)/(xU-ca))
          val::T = ca + (xU-ca)*(term)^(MC_param.mu+1)
          dval::T = (MC_param.mu+1)*(term)^(MC_param.mu)
          return val, dval
        end
end
=#
function max(x::SMCg{N,V,T},c::T) where {N,V,T<:AbstractFloat}
    xL::T = x.Intv.lo
    xU::T = x.Intv.hi
    maxxL::T = max(xL,c)
    maxxU::T = max(xU,c)
    if (MC_param.mu >= 1)
      if (x.cc <= xL)
        cv::T,dcv::T = cv_max(x.cc,xL,xU,c)
        cv_grad::SVector{N,T} = dcv*x.cv_grad
      elseif (x.Intv.lo >= x.cv)
        cv,dcv = cv_max(xL,xL,xU,c)
        cv_grad = dcv*x.cv_grad
      else
        cv,dcv = cv_max(x.cv,xL,xU,c)
        cv_grad = cv*x.cv_grad
      end

      # calc concave term
      if (neq(xU,xL))
        dcc::T = (maxxU-maxxL)/(xU-xL)
        cc_grad::SVector{N,T} = dcc*x.cc_grad
        if (x.cc <= xU)
          cc::T = maxxL + dcc*(x.cc - xL)
        elseif (xU >= x.cv)
          cc = maxxL + dcc*(xU - xL)
        else
          cc = maxxL + dcc*(x.cv - xL)
        end
      else
        cc = maxxU
        cc_grad = zeros(SVector{N,T})
      end
    else
      # calc convex term
      if (xL < x.cv)
        cv = max(x.cv,c)
        cv_grad = (x.cv > c) ? x.cv_grad : zeros(SVector{N,T})
      elseif (xL > x.cc)
        cv = max(x.cc,c)
        cv_grad = (x.cc > c) ? x.cc_grad : zeros(SVector{N,T})
      else
        cv = max(xL,c)
        cv_grad = zeros(SVector{N,T})
      end

      # calc concave term
      if (neq(xU,xL))
          dcc = (maxxU-maxxL)/(xU-xL)
          if (xU < x.cv)
            cc = maxxL + dcc*(x.cv - xL)
            cc_grad = dcc*x.cv_grad
          elseif (xU > x.cc)
            cc = maxxL + dcc*(x.cc - xL)
            cc_grad = dcc*x.cc_grad
          else
            cc = maxxL + dcc*(xU - xL)
            cc_grad = zeros(SVector{N,T})
          end
      else
          cc = xU
          cc_grad = zeros(SVector{N,T})
      end
      # applies cut operator
      cv,cc,cv_grad,cc_grad = cut(maxxL,maxxU,cv,cc,cv_grad,cc_grad)
    end
    # ((V<:AbstractMCInterval) ? V(maxxL,maxxU) : max(x.Intv,c))
    return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, max(x.Intv,c), x.cnst)
end

# defines functions on which bivariant maximum mapping from Khan 2016
@inline function psil_max(x::T,y::T,lambda::V,nu::V,
                          f1::SMCg{N,V,T},f2::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
   if (nu.hi<=lambda.lo)
     val::T = x
   elseif (lambda.hi<=nu.lo)
     val = y
   elseif ((nu.lo<=lambda.lo)&&(lambda.lo<nu.hi))
     val = x+(nu.hi-lambda.lo)*max(zero(T),((y-x)/(nu.hi-lambda.lo)))^(MC_param.mu+1)
   else
     val =  y + (lambda.hi-nu.lo)*max(zero(T),(x-y)/(lambda.hi-nu.lo))^(MC_param.mu+1)
   end
   if (nu.hi <= lambda.lo)
     grad_val::SVector{N,T} = f1.cv_grad
   elseif (lambda.hi <= nu.lo)
     grad_val = f1.cc_grad
   else
     grad_val = max(zero(T),psil_max_dx(x,y,lambda,nu))*f1.cv_grad +
                min(zero(T),psil_max_dx(x,y,lambda,nu))*f1.cc_grad +
                max(zero(T),psil_max_dy(x,y,lambda,nu))*f2.cv_grad +
                min(zero(T),psil_max_dy(x,y,lambda,nu))*f2.cc_grad
   end
   return val,grad_val
end
@inline function thetar(x::T,y::T,lambda::V,nu::V) where {V,T<:AbstractFloat}
    return (max(lambda.lo,nu.lo) + max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi) +
    max(lambda.hi,nu.lo))*max(zero(T),((lambda.hi-x)/(lambda.hi-lambda.lo)-(y-nu.lo)/(nu.hi-nu.lo)))^(MC_param.mu+1)
end
function psil_max_dx(x::T,y::T,lambda::V,nu::V) where {V,T<:AbstractFloat}
  if (nu.lo <= lambda.lo < nu.hi)
    return one(T)-(MC_param.mu+1)*max(zero(T),(y-x)/(nu.hi-lambda.lo))^MC_param.mu
  else
    return (MC_param.mu+1)*max(zero(T),(x-y)/(lambda.hi-nu.lo))^MC_param.mu
  end
end
function psil_max_dy(x::T,y::T,lambda::V,nu::V) where {V,T<:AbstractFloat}
  if (nu.lo <= lambda.lo < nu.hi)
    return (MC_param.mu+1)*max(zero(T),(y-x)/(nu.hi-lambda.lo))^MC_param.mu
  else
    return one(T)-(MC_param.mu+1)*max(zero(T),(x-y)/(lambda.hi-nu.lo))^MC_param.mu
  end
end

@inline function psir_max(x::T,y::T,xgrad::SVector{N,T},ygrad::SVector{N,T},
                          lambda::V,nu::V) where {N,V,T<:AbstractFloat}
    if (nu.hi<=lambda.lo)
      return x,xgrad
    elseif (lambda.hi<=nu.lo)
      return y,ygrad
    else
      val::T = max(lambda.hi,nu.hi)-(max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi))*
          ((lambda.hi-x)/(lambda.hi-lambda.lo))-
          (max(lambda.hi,nu.hi)-max(lambda.hi,nu.lo))*((nu.hi-y)/(nu.hi-nu.lo))
          +thetar(x,y,nu,lambda)
      coeff = [(max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi))/(lambda.hi-lambda.lo)
               (max(lambda.hi,nu.hi)-max(nu.lo,lambda.hi))/(nu.hi-nu.lo)
               (MC_param.mu+1)*(max(lambda.hi,nu.lo)+max(lambda.lo,nu.hi)-max(lambda.lo,nu.lo)-max(lambda.hi,nu.hi))
               max(0,((lambda.hi-x)/(lambda.hi-lambda.lo))-((y-nu.lo)/(nu.hi-nu.lo)))^MC_param.mu
               1/(lambda.hi-lambda.lo)
               1/(nu.hi-nu.lo)
               ]
      grad_val::SVector{N,T} = coeff[1]*xgrad + coeff[2]*ygrad +
                 coeff[3]*coeff[4]*(coeff[5]*xgrad+coeff[6]*ygrad)
      return val,grad_val
    end
end

function max(x::SMCg{N,V,T},y::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
    maxIntv = max(x.Intv,y.Intv)
    maxxL = maxIntv.lo
    maxxU = maxIntv.hi
    if (MC_param.mu >= 1)
      cc::T = zero(T)
      cv::T = zero(T)
      temp_mid::T = zero(T)
      if ((y.Intv.hi<=x.Intv.lo)||(x.Intv.hi<=y.Intv.lo))
        cv,cv_grad::SVector{N,T} = psil_max(x.cv,y.cv,x.Intv,y.Intv,x,y)
      elseif ((y.Intv.lo<=x.Intv.lo) & (x.Intv.lo<y.Intv.hi))
        temp_mid,blank = mid3(x.cv,x.cc,y.cv-(y.Intv.hi-x.Intv.lo)*(MC_param.mu+1)^(-1/MC_param.mu))
        cv,cv_grad = psil_max(temp_mid,y.cv,x.Intv,y.Intv,x,y)
      elseif ((x.Intv.lo<y.Intv.lo) & (y.Intv.lo<x.Intv.hi))
        temp_mid,blank = mid3(y.cv,y.cc,x.cv-(x.Intv.hi-y.Intv.lo)*(MC_param.mu+1)^(-1/MC_param.mu))
        cv,cv_grad = psil_max(x.cv,temp_mid,x.Intv,y.Intv,x,y)
      end
      cc,cc_grad::SVector{N,T} = psir_max(x.cc,y.cc,x.cc_grad,y.cv_grad,x.Intv,y.Intv)
      return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, max(x.Intv,y.Intv),(x.cnst && y.cnst))
    elseif (x.Intv.hi <= y.Intv.lo)
      cc = y.cc
      cc_grad = y.cnst ? zeros(y.cc_grad) : y.cc_grad
    elseif (x.Intv.lo >= y.Intv.hi)
      cc = x.cc
      cc_grad = x.cnst ? zeros(x.cc_grad) : x.cc_grad
      #=
    elseif (MC_param.multivar_refine)
      maxLL::T = max(x.Intv.lo,y.Intv.lo)
      maxLU::T = max(x.Intv.lo,y.Intv.hi)
      maxUL::T = max(x.Intv.hi,y.Intv.lo)
      maxUU::T = max(x.Intv.hi,y.Intv.hi)
      thin1::Bool = (diam(x.Intv) == zero(T))
      thin2::Bool = (diam(y.Intv) == zero(T))
      r11::T = thin1 ? zero(T) : (maxUL-maxLL)/diam(x.Intv)
      r21::T = thin1 ? zero(T) : (maxLU-maxUU)/diam(x.Intv)
      r12::T = thin2 ? zero(T) : (maxLU-maxLL)/diam(y.Intv)
      r22::T = thin2 ? zero(T) : (maxUL-maxUU)/diam(y.Intv)
      cc1::T = maxLL + r11*(x.cc-x.Intv.lo) + r12*(y.cc-y.Intv.lo)
      cc2::T = maxUU - r21*(x.cc-x.Intv.hi) - r22*(y.cc-y.Intv.hi)
      if (cc1 <= cc2)
        cc = cc1
        cc_grad =  (x.cnst ? zeros(y.cc_grad) : r11*x.cc_grad) + (y.cnst ? zeros(x.cc_grad) : r12*y.cc_grad)
      else
        cc = cc2
        cc_grad = -(x.cnst ? zeros(y.cc_grad) : r21*x.cc_grad) - (y.cnst ? zeros(x.cc_grad) : r22*y.cc_grad)
      end
      cv,cc,cv_grad,cc_grad = cut(maxxL,maxxU,cv,cc,cv_grad,cc_grad)
      =#
    else
      ccMC::SMCg{N,V,T} = (x+y+abs(x-y))/2
      cc = ccMC.cc
      cc_grad = ccMC.cc_grad
      cv = max(x.cv,y.cv)
      cv_grad = (x.cv > y.cv) ? (x.cnst ? zeros(x.cv_grad): x.cv_grad) :
                                (y.cnst ? zeros(y.cv_grad): y.cv_grad)
      cv,cc,cv_grad,cc_grad = cut(maxxL,maxxU,cv,cc,cv_grad,cc_grad)
      cnst = y.cnst ? x.cnst : (x.cnst ? y.cnst : (x.cnst || y.cnst) )

      return SMCg{N,V,T}(cc, cv, cc_grad, cv_grad, max(x.Intv,y.Intv),cnst)
    end
end

@inline function maxcv(x::SMCg{N,V,T},y::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
        cv::T = zero(T)
        temp_mid::T = zero(T)
        if ((y.Intv.hi<=x.Intv.lo)||(x.Intv.hi<=y.Intv.lo))
            cv = psil_max(x.cv,y.cv,x.Intv,y.Intv)
        elseif ((y.Intv.lo<=x.Intv.lo) & (x.Intv.lo<y.Intv.hi))
          temp_mid = mid3(x.cv,x.cc,y.cv-(y.Intv.hi-x.Intv.lo)*(MC_param.mu+1)^(-1.0/MC_param.mu))
          cv = psil_max(temp_mid,y.cv,x.Intv,y.Intv)
        elseif ((x.Intv.lo<y.Intv.lo) & (y.Intv.lo<x.Intv.hi))
          temp_mid = mid3(y.cv,y.cc,x.cv-(x.Intv.hi-y.Intv.lo)*(MC_param.mu+1)^(-1.0/MC_param.mu))
          cv = psil_max(x.cv,temp_mid,x.Intv,y.Intv)
        end
        return cv
end
@inline function maxcc(x,y::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
        cc::T = psir_max(x.cc,y.cc,x.Intv,y.Intv)
end
@inline mincv(x,y::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = - maxcc(-x,-y)
@inline min(x::SMCg{N,V,T},y::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = -max(-x,-y)
#=
@inline max(x::SMCg{N,V,T},y::V) where {N,V,T<:AbstractFloat} = max(x,SMCg{N,V,T}(y))
@inline max(y::V,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = max(x,SMCg{N,V,T}(y))
@inline min(x::SMCg{N,V,T},y::V) where {N,V,T<:AbstractFloat} = min(x,SMCg{N,V,T}(y))
@inline min(y::V,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = min(x,SMCg{N,V,T}(y))
=#
