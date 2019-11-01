"""
    affine_exp!(x,p::Vector{MC{N}},p_ref::Vector{MC{N}},
               xa::Vector{MC{N}},xA::Vector{MC{N}},z::Vector{MC{N}},
               opt::Array{Any})

Computates the affine relaxations of the state variable. Inputs are:
* `x::Vector{MC{N,T}}`: State variable relaxation
* `p::Vector{MC{N,T}}`: Decision variable relaxation
* `p_ref::Vector{MC{N,T}}`: Reference variable relaxation
* `xa::Vector{MC{N,T}}`: Lower affine relaxation of the state variable
* `xA::Vector{MC{N,T}}`: Upper affine relaxation of the state variable
* `z::Vector{MC{N,T}}`: Affine function in `X`
* `opt::Array{Any,1}`: `[np,nx,lambda]` values for relaxation
Returns the tuple `(xa,xA,z)`:
* `xa::Vector{MC{N,T}}`: Lower affine relaxation of the state variable
* `xA::Vector{MC{N,T}}`: Upper affine relaxation of the state variable
* `z::Vector{MC{N,T}}`: Affine function in X
--------------------------------------------------------------------------------
"""
function affine_exp!(x, p::Vector{MC{N,T}},
                     p_ref::Vector{MC{N,T}}, xa::Vector{MC{N,T}},
                     xA::Vector{MC{N,T}}, z::Vector{MC{N,T}},
                     nx::Int64, lambda::Float64) where {N, T<:RelaxTag}

  S1 = zero(MC{N,T})
  S2 = zero(MC{N,T})
  S3 = zero(MC{N,T})
  for i = 1:nx
   S1 = zero(MC{N,T})
   S2 = zero(MC{N,T})
   S3 = zero(MC{N,T})
   for j = 1:N
      @inbounds S1 = S1 + (p[j]-p_ref[j])*x[i].cv_grad[j]
      @inbounds S2 = S2 + (p[j]-p_ref[j])*x[i].cc_grad[j]
      @inbounds S3 = S3 + (lambda*x[i].cv_grad[j]+(1.0-lambda)*x[i].cc_grad[j])*(p[j]-p_ref[j])
   end
   @inbounds temp1 = x[i].cv + S1
   @inbounds temp2 = x[i].cc + S2
   @inbounds temp3 = lambda*x[i].cv+(1.0-lambda)*x[i].cc+S3
   @inbounds xa[i] = MC{N,T}(temp1.cv,temp1.cc,temp1.Intv,x[i].cv_grad,x[i].cv_grad,S1.cnst)
   @inbounds xA[i] = MC{N,T}(temp2.cv,temp2.cc,temp2.Intv,x[i].cc_grad,x[i].cc_grad,S2.cnst)
   @inbounds z[i] = MC{N,T}(temp3.cv,temp3.cc,temp3.Intv,
                            lambda*x[i].cv_grad+(1.0-lambda)*x[i].cc_grad,
                            lambda*x[i].cv_grad+(1.0-lambda)*x[i].cc_grad,
                            S3.cnst)
  end
  return
end

"""
    correct_exp!(z_mc::Vector{MC{N}},x_mc::Vector{MC{N}},
                 X::Vector{IntervalType},nx::Int,np::Int,
                 epsv::Float64) where N

Corrects the relaxation of the state variable `x_mc` if the affine relaxation,
'z_mc', exceeds the interval bounds `xL` or `xU`.
* `z_mc::Vector{MC{N}}`: Affine relaxation
* `x_mc::Vector{MC{N}}`: Relaxation of state variable
* `X::Vector{IntervalType}`: Lower bound on state vector
* `nx::Int64`: Size of the state vector
* `np::Int64`: Size of the decision vector
* `epsv::Float64`: Tolerance for checking that subgradient exceeds bound
"""
function correct_exp!(xa_mc::Vector{MC{N,T}}, xA_mc::Vector{MC{N,T}},
                      z_mc::Vector{MC{N,T}}, x_mc, X::Vector{Interval{Float64}},
                      nx::Int64, epsv::Float64) where {N, T<:RelaxTag}
  zero_grad = zeros(SVector{N,Float64})
  for i = 1:nx
    @inbounds z_mc_term = z_mc[i]
    z_mc_lo = z_mc_term.Intv.lo
    z_mc_hi = z_mc_term.Intv.hi
    @inbounds X_term = X[i]
    X_lo = X_term.lo
    X_hi = X_term.hi
    if (z_mc_lo - epsv < X_lo) && (z_mc_hi + epsv > X_hi)
      @inbounds x_mc[i] = MC{N,T}(X_lo, X_hi, X_term, zero_grad, zero_grad,true)
    end
    if (z_mc_lo - epsv < X_lo)
      @inbounds x_mc[i] = MC{N,T}(X_lo, x_mc[i].cc, Interval(X_lo, x_mc[i].Intv.hi),
                                  zero_grad, x_mc[i].cc_grad, x_mc[i].cnst)
    end
    if (z_mc_hi + epsv > X_hi)
      @inbounds x_mc[i] = MC{N,T}(x_mc[i].cv, X_hi,Interval(x_mc[i].Intv.lo, X_hi),
                                  x_mc[i].cv_grad, zero_grad, x_mc[i].cnst)
    end
  end
  return
end
