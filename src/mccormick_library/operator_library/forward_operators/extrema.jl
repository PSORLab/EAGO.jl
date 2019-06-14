function cv_max(x::Float64, xL::Float64, xU::Float64, ca::Float64)
    if (xU<=ca)
        return ca, zero(Float64)
    elseif (ca<=xL)
        return x, one(Float64)
    else
        term::Float64 = max(zero(Float64),(x-ca)/(xU-ca))
        val::Float64 = ca + (xU-ca)*(term)^(MC_param.mu+1)
        dval::Float64 = (MC_param.mu+1)*(term)^(MC_param.mu)
        return val, dval
    end
end

function max(x::MC{N}, c::Float64) where N
    xL = x.Intv.lo
    xU = x.Intv.hi
    maxxL = max(xL,c)
    maxxU = max(xU,c)
    if (MC_param.mu >= 1)
      if (x.cc <= xL)
        cv, dcv = cv_max(x.cc, xL, xU, c)
        cv_grad = dcv*x.cv_grad
      elseif (x.Intv.lo >= x.cv)
        cv, dcv = cv_max(xL, xL, xU, c)
        cv_grad = dcv*x.cv_grad
      else
        cv, dcv = cv_max(x.cv, xL, xU, c)
        cv_grad = cv*x.cv_grad
      end

      # calc concave term
      if (neq(xU,xL))
        dcc = (maxxU-maxxL)/(xU-xL)
        cc_grad = dcc*x.cc_grad
        if (x.cc <= xU)
          cc = maxxL + dcc*(x.cc - xL)
        elseif (xU >= x.cv)
          cc = maxxL + dcc*(xU - xL)
        else
          cc = maxxL + dcc*(x.cv - xL)
        end
      else
        cc = maxxU
        cc_grad = zeros(SVector{N,Float64})
      end
    else
      # calc convex term
      if (xL < x.cv)
        cv = max(x.cv,c)
        cv_grad = (x.cv > c) ? x.cv_grad : zeros(SVector{N,Float64})
      elseif (xL > x.cc)
        cv = max(x.cc,c)
        cv_grad = (x.cc > c) ? x.cc_grad : zeros(SVector{N,Float64})
      else
        cv = max(xL,c)
        cv_grad = zeros(SVector{N,Float64})
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
            cc_grad = zeros(SVector{N,Float64})
          end
      else
          cc = xU
          cc_grad = zeros(SVector{N,Float64})
      end
      # applies cut operator
      cv,cc,cv_grad,cc_grad = cut(maxxL,maxxU,cv,cc,cv_grad,cc_grad)
    end
    return MC{N}(cv, cc, max(x.Intv,c), cv_grad, cc_grad, x.cnst)
end

# defines functions on which bivariant maximum mapping from Khan 2016
function psil_max(x,y,lambda,nu,f1::MC{N},f2::MC{N}) where N
   if (nu.hi<=lambda.lo)
     val = x
   elseif (lambda.hi<=nu.lo)
     val = y
   elseif ((nu.lo<=lambda.lo)&&(lambda.lo<nu.hi))
     val = x+(nu.hi-lambda.lo)*max(0.0,((y-x)/(nu.hi-lambda.lo)))^(MC_param.mu+1)
   else
     val =  y + (lambda.hi-nu.lo)*max(0.0,(x-y)/(lambda.hi-nu.lo))^(MC_param.mu+1)
   end
   if (nu.hi <= lambda.lo)
     grad_val = f1.cv_grad
   elseif (lambda.hi <= nu.lo)
     grad_val = f1.cc_grad
   else
     grad_val = max(0.0,psil_max_dx(x,y,lambda,nu))*f1.cv_grad +
                min(0.0,psil_max_dx(x,y,lambda,nu))*f1.cc_grad +
                max(0.0,psil_max_dy(x,y,lambda,nu))*f2.cv_grad +
                min(0.0,psil_max_dy(x,y,lambda,nu))*f2.cc_grad
   end
   return val,grad_val
end
function thetar(x,y,lambda,nu)
    return (max(lambda.lo,nu.lo) + max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi) +
    max(lambda.hi,nu.lo))*max(0.0,((lambda.hi-x)/(lambda.hi-lambda.lo)-(y-nu.lo)/(nu.hi-nu.lo)))^(MC_param.mu+1)
end
function psil_max_dx(x,y,lambda,nu)
  if (nu.lo <= lambda.lo < nu.hi)
    return 1.0-(MC_param.mu+1)*max(0.0,(y-x)/(nu.hi-lambda.lo))^MC_param.mu
  else
    return (MC_param.mu+1)*max(0.0,(x-y)/(lambda.hi-nu.lo))^MC_param.mu
  end
end
function psil_max_dy(x,y,lambda,nu)
  if (nu.lo <= lambda.lo < nu.hi)
    return (MC_param.mu+1)*max(0.0,(y-x)/(nu.hi-lambda.lo))^MC_param.mu
  else
    return 1.0-(MC_param.mu+1)*max(0.0,(x-y)/(lambda.hi-nu.lo))^MC_param.mu
  end
end

function psir_max(x,y,xgrad,ygrad,lambda,nu)
    if (nu.hi<=lambda.lo)
      return x,xgrad
    elseif (lambda.hi<=nu.lo)
      return y,ygrad
    else
      val = max(lambda.hi,nu.hi)-(max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi))*
          ((lambda.hi-x)/(lambda.hi-lambda.lo))-
          (max(lambda.hi,nu.hi)-max(lambda.hi,nu.lo))*((nu.hi-y)/(nu.hi-nu.lo))+thetar(x,y,nu,lambda)
      grad_val = (max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi))/(lambda.hi-lambda.lo)*xgrad +
                 (max(lambda.hi,nu.hi)-max(nu.lo,lambda.hi))/(nu.hi-nu.lo)*ygrad +
                 (MC_param.mu+1)*(max(lambda.hi,nu.lo)+max(lambda.lo,nu.hi)-max(lambda.lo,nu.lo)-max(lambda.hi,nu.hi))*
                 max(0.0,((lambda.hi-x)/(lambda.hi-lambda.lo))-((y-nu.lo)/(nu.hi-nu.lo)))^MC_param.mu*
                 ((1.0/(lambda.hi-lambda.lo))*xgrad+(1.0/(nu.hi-nu.lo))*ygrad)
      return val,grad_val
    end
end

function max(x::MC{N},y::MC{N}) where N
    maxIntv = max(x.Intv,y.Intv)
    maxxL = maxIntv.lo
    maxxU = maxIntv.hi
    if (MC_param.mu >= 1)
      cc = 0.0
      cv = 0.0
      temp_mid = 0.0
      if ((y.Intv.hi<=x.Intv.lo)||(x.Intv.hi<=y.Intv.lo))
        cv,cv_grad = psil_max(x.cv,y.cv,x.Intv,y.Intv,x,y)
      elseif ((y.Intv.lo<=x.Intv.lo) & (x.Intv.lo<y.Intv.hi))
        temp_mid,blank = mid3(x.cv,x.cc,y.cv-(y.Intv.hi-x.Intv.lo)*(MC_param.mu+1)^(-1/MC_param.mu))
        cv,cv_grad = psil_max(temp_mid,y.cv,x.Intv,y.Intv,x,y)
      elseif ((x.Intv.lo<y.Intv.lo) & (y.Intv.lo<x.Intv.hi))
        temp_mid,blank = mid3(y.cv,y.cc,x.cv-(x.Intv.hi-y.Intv.lo)*(MC_param.mu+1)^(-1/MC_param.mu))
        cv,cv_grad = psil_max(x.cv,temp_mid,x.Intv,y.Intv,x,y)
      end
      cc,cc_grad = psir_max(x.cc,y.cc,x.cc_grad,y.cv_grad,x.Intv,y.Intv)
      return MC{N}(cv, cc, maxIntv, cv_grad, cc_grad, (x.cnst && y.cnst))
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
      ccMC = (x+y+abs(x-y))/2.0
      cc = ccMC.cc
      cc_grad = ccMC.cc_grad
    end
      cv = max(x.cv,y.cv)
      cv_grad = (x.cv > y.cv) ? (x.cnst ? zeros(x.cv_grad) : x.cv_grad) :
                                (y.cnst ? zeros(y.cv_grad) : y.cv_grad)
      cv,cc,cv_grad,cc_grad = cut(maxxL,maxxU,cv,cc,cv_grad,cc_grad)
      cnst = y.cnst ? x.cnst : (x.cnst ? y.cnst : (x.cnst || y.cnst) )
      temp = MC{N}(cv, cc, maxIntv, cv_grad, cc_grad, cnst)
      return temp
end

function maxcv(x::MC,y::MC)
        cv = 0.0
        temp_mid = 0.0
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
function maxcc(x::MC,y::MC)
        psir_max(x.cc,y.cc,x.Intv,y.Intv)
end
mincv(x::MC,y::MC) = - maxcc(-x,-y)
min(x::MC,y::MC) = -max(-x,-y)
#=
@inline max(x::SMCg{N,V,T},y::V) where {N,V,T<:AbstractFloat} = max(x,SMCg{N,V,T}(y))
@inline max(y::V,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = max(x,SMCg{N,V,T}(y))
@inline min(x::SMCg{N,V,T},y::V) where {N,V,T<:AbstractFloat} = min(x,SMCg{N,V,T}(y))
@inline min(y::V,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = min(x,SMCg{N,V,T}(y))
=#
