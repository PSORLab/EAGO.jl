@inline deriv_max(x::Float64, c::Float64) = (x < c) ? 0.0 : 1.0
@inline function cv_max(x::Float64, xL::Float64, xU::Float64, a::Float64)
    (xU <= a) && (return a, 0.0)
    (a <= xL) && (return x, 1.0)
    term::Float64 = max(0.0, (x - a)/(xU - a))
    val::Float64 = a + (xU - a)*(term)^(MC_DIFF_MU + 1)
    dval::Float64 = (MC_DIFF_MU + 1)*(term)^(MC_DIFF_MU)
    return val, dval
end
@inline cc_max(x::Float64, xL::Float64, xU::Float64, a::Float64) = dline_seg(max, deriv_max, x, xL, xU, a)
@inline function max_kernel(x::MC{N, Diff}, c::Float64, z::Interval{Float64}) where N
    midcv, cv_id = mid3(x.cc, x.cv, x.Intv.lo)
    midcc, cc_id = mid3(x.cc, x.cv, x.Intv.hi)
    cv, dcv = cv_max(midcv, x.Intv.lo, x.Intv.hi, c)
    cc, dcc = cc_max(midcc, x.Intv.lo, x.Intv.hi, c)
    gcc1,gdcc1 = cc_max(x.cv, x.Intv.lo, x.Intv.hi, c)
    gcv1,gdcv1 = cv_max(x.cv, x.Intv.lo, x.Intv.hi, c)
    gcc2,gdcc2 = cc_max(x.cc, x.Intv.lo, x.Intv.hi, c)
    gcv2,gdcv2 = cv_max(x.cc, x.Intv.lo, x.Intv.hi, c)
    cv_grad = max(0.0, gdcv1)*x.cv_grad + min(0.0, gdcv2)*x.cc_grad
    cc_grad = min(0.0, gdcc1)*x.cv_grad + max(0.0, gdcc2)*x.cc_grad
    return MC{N,Diff}(cv, cc, z, cv_grad, cc_grad, x.cnst)
end
@inline function max_kernel(x::MC{N, NS}, c::Float64, z::Interval{Float64}) where N
    midcv, cv_id = mid3(x.cc, x.cv, x.Intv.lo)
    midcc, cc_id = mid3(x.cc, x.cv, x.Intv.hi)
    cv, dcv = cv_max(midcv, x.Intv.lo, x.Intv.hi, c)
    cc, dcc = cc_max(midcc, x.Intv.lo, x.Intv.hi, c)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    cv, cc, cv_grad, cc_grad = cut(z.lo, z.hi, cv, cc, cv_grad, cc_grad)
    return MC{N,NS}(cv, cc, z, cv_grad, cc_grad, x.cnst)
end
@inline max(x::MC, c::Float64) = max_kernel(x, c, max(x.Intv, c))
@inline min_kernel(x::MC, c::Float64, z::Interval{Float64}) = -max(-x,-c)
@inline min(x::MC, c::Float64) = -max(-x,-c)

# defines functions on which bivariant maximum mapping from Khan 2016
@inline function psil_max(x::Float64, y::Float64, lambda::Interval{Float64}, nu::Interval{Float64}, f1::MC{N, Diff}, f2::MC{N, Diff}) where N
   if (nu.hi <= lambda.lo)
     val = x
   elseif (lambda.hi <= nu.lo)
     val = y
   elseif ((nu.lo<=lambda.lo)&&(lambda.lo<nu.hi))
     val = x + (nu.hi-lambda.lo)*max(0.0,((y-x)/(nu.hi-lambda.lo)))^(MC_DIFF_MU+1)
   else
     val =  y + (lambda.hi-nu.lo)*max(0.0,(x-y)/(lambda.hi-nu.lo))^(MC_DIFF_MU+1)
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
@inline function thetar(x::Float64, y::Float64, lambda::Interval{Float64}, nu::Interval{Float64})
    return (max(lambda.lo,nu.lo) + max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi) +
    max(lambda.hi,nu.lo))*max(0.0,((lambda.hi-x)/(lambda.hi-lambda.lo)-(y-nu.lo)/(nu.hi-nu.lo)))^(MC_DIFF_MU+1)
end
@inline function psil_max_dx(x::Float64, y::Float64, lambda::Interval{Float64}, nu::Interval{Float64})
  (nu.lo <= lambda.lo < nu.hi) && (1.0-(MC_DIFF_MU+1)*max(0.0,(y-x)/(nu.hi-lambda.lo))^MC_DIFF_MU)
  return (MC_DIFF_MU+1)*max(0.0,(x-y)/(lambda.hi-nu.lo))^MC_DIFF_MU
end
@inline function psil_max_dy(x::Float64, y::Float64, lambda::Interval{Float64}, nu::Interval{Float64})
  (nu.lo <= lambda.lo < nu.hi) && (return (MC_DIFF_MU+1)*max(0.0,(y-x)/(nu.hi-lambda.lo))^MC_DIFF_MU)
  return 1.0-(MC_DIFF_MU+1)*max(0.0,(x-y)/(lambda.hi-nu.lo))^MC_DIFF_MU
end
@inline function psir_max(x::Float64, y::Float64, xgrad::SVector{N,Float64}, ygrad::SVector{N,Float64},
                  lambda::Interval{Float64}, nu::Interval{Float64}) where N
    if (nu.hi <= lambda.lo)
      return x,xgrad
    elseif (lambda.hi <= nu.lo)
      return y,ygrad
    else
      val = max(lambda.hi,nu.hi)-(max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi))*
          ((lambda.hi-x)/(lambda.hi-lambda.lo))-
          (max(lambda.hi,nu.hi)-max(lambda.hi,nu.lo))*((nu.hi-y)/(nu.hi-nu.lo))+thetar(x,y,nu,lambda)
      grad_val = (max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi))/(lambda.hi-lambda.lo)*xgrad +
                 (max(lambda.hi,nu.hi)-max(nu.lo,lambda.hi))/(nu.hi-nu.lo)*ygrad +
                 (MC_DIFF_MU+1)*(max(lambda.hi,nu.lo)+max(lambda.lo,nu.hi)-max(lambda.lo,nu.lo)-max(lambda.hi,nu.hi))*
                 max(0.0,((lambda.hi-x)/(lambda.hi-lambda.lo))-((y-nu.lo)/(nu.hi-nu.lo)))^MC_DIFF_MU*
                 ((1.0/(lambda.hi-lambda.lo))*xgrad+(1.0/(nu.hi-nu.lo))*ygrad)
      return val, grad_val
    end
end

@inline function max_kernel(x::MC{N, NS}, y::MC{N, NS}, z::Interval{Float64}) where N
    maxxL = z.lo
    maxxU = z.hi
    if (x.Intv.hi <= y.Intv.lo)
        cc = y.cc
        cc_grad = y.cnst ? zero(SVector{N,Float64}) : y.cc_grad
    elseif (x.Intv.lo >= y.Intv.hi)
        cc = x.cc
        cc_grad = x.cnst ? zero(SVector{N,Float64}) : x.cc_grad
    else
        ccMC = (x + y + abs(x-y))/2.0
        cc = ccMC.cc
        cc_grad = ccMC.cc_grad
    end
    cv = max(x.cv,y.cv)
    cv_grad = (x.cv > y.cv) ? (x.cnst ? zero(SVector{N,Float64}) : x.cv_grad) :
                              (y.cnst ? zero(SVector{N,Float64}) : y.cv_grad)
    cv, cc, cv_grad, cc_grad = cut(maxxL, maxxU, cv,cc,cv_grad,cc_grad)
    return MC{N, NS}(cv, cc, z, cv_grad, cc_grad, y.cnst ? x.cnst : (x.cnst ? y.cnst : (x.cnst || y.cnst) ))
end
@inline function max_kernel(x::MC{N, Diff}, y::MC{N, Diff}, z::Interval{Float64}) where N
    maxxL = z.lo
    maxxU = z.hi
    cc = 0.0
    cv = 0.0
    temp_mid = 0.0
    if (y.Intv.hi <= x.Intv.lo) || (x.Intv.hi <= y.Intv.lo)
        cv,cv_grad = psil_max(x.cv, y.cv, x.Intv, y.Intv, x, y)
    elseif (y.Intv.lo <= x.Intv.lo) & (x.Intv.lo < y.Intv.hi)
        temp_mid, blank = mid3(x.cv, x.cc, y.cv - (y.Intv.hi - x.Intv.lo)*(MC_DIFF_MU+1)^(-1/MC_DIFF_MU))
        cv,cv_grad = psil_max(temp_mid, y.cv, x.Intv, y.Intv, x, y)
    elseif (x.Intv.lo < y.Intv.lo) & (y.Intv.lo < x.Intv.hi)
        temp_mid,blank = mid3(y.cv, y.cc, x.cv - (x.Intv.hi - y.Intv.lo)*(MC_DIFF_MU+1)^(-1/MC_DIFF_MU))
        cv,cv_grad = psil_max(x.cv,temp_mid,x.Intv,y.Intv,x,y)
    end
    cc,cc_grad = psir_max(x.cc, y.cc, x.cc_grad, y.cv_grad, x.Intv, y.Intv)
    return MC{N, Diff}(cv, cc, max(x.Intv, y.Intv), cv_grad, cc_grad, (x.cnst && y.cnst))
end
@inline max(x::MC, y::MC) = max_kernel(x, y, max(x.Intv, y.Intv))

@inline min_kernel(x::MC, y::MC, z::Interval{Float64}) = -max(-x,-y)
@inline min(x::MC,y::MC) = -max(-x,-y)
