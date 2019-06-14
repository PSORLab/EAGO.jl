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
function affine_exp!(x, p::Vector{MC{N}},
                     p_ref::Vector{MC{N}}, xa::Vector{MC{N}},
                     xA::Vector{MC{N}}, z::Vector{MC{N}},
                     nx::Int, lambda::Float64) where N

  S1::MC{N},S2::MC{N},S3::MC{N} = zero(x[1]),zero(x[1]),zero(x[1])
  for i = 1:nx
   S1 = zero(x[1])
   S2 = zero(x[1])
   S3 = zero(x[1])
   for j = 1:N
      S1 = S1 + (p[j]-p_ref[j])*x[i].cv_grad[j]
      S2 = S2 + (p[j]-p_ref[j])*x[i].cc_grad[j]
      S3 = S3 + (lambda*x[i].cv_grad[j]+(1.0-lambda)*x[i].cc_grad[j])*(p[j]-p_ref[j])
   end
   temp1 = x[i].cv + S1
   temp2 = x[i].cc + S2
   temp3 = lambda*x[i].cv+(1.0-lambda)*x[i].cc+S3
   xa[i] = MC{N}(temp1.cv,temp1.cc,temp1.Intv,x[i].cv_grad,x[i].cv_grad,S1.cnst)
   xA[i] = MC{N}(temp2.cv,temp2.cc,temp2.Intv,x[i].cc_grad,x[i].cc_grad,S2.cnst)
   z[i] = MC{N}(temp3.cv,temp3.cc,temp3.Intv,
                lambda*x[i].cv_grad+(1.0-lambda)*x[i].cc_grad,
                lambda*x[i].cv_grad+(1.0-lambda)*x[i].cc_grad,
                S3.cnst)
  end
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
function correct_exp!(xa_mc::Vector{MC{N}},xA_mc::Vector{MC{N}},z_mc::Vector{MC{N}},x_mc,
                     X::Vector{IntervalType},nx::Int,
                     epsv::Float64) where N
  zero_grad::SVector{N,Float64} = @SVector zeros(Float64,N)
  for i = 1:nx
    if (z_mc[i].Intv.lo-epsv < X[i].lo) && (z_mc[i].Intv.hi+epsv > X[i].hi)
      x_mc[i] = MC{N}(X[i].lo,X[i].hi,X[i],zero_grad,zero_grad,true)
      #xa_mc[i] = MC{N}(IntervalType(X[i].lo))
      #xA_mc[i] = MC{N}(IntervalType(X[i].hi))
    end
    if (z_mc[i].Intv.lo-epsv < X[i].lo)
      x_mc[i] = MC{N}(X[i].lo,x_mc[i].cc,IntervalType(X[i].lo,x_mc[i].Intv.hi),zero_grad,x_mc[i].cc_grad,x_mc[i].cnst)
      #xa_mc[i] = MC{N}(IntervalType(X[i].lo))
    end
    if (z_mc[i].Intv.hi+epsv > X[i].hi)
      x_mc[i] = MC{N}(x_mc[i].cv,X[i].hi,IntervalType(x_mc[i].Intv.lo,X[i].hi),x_mc[i].cv_grad,zero_grad,x_mc[i].cnst)
      #xA_mc[i] = MC{N}(IntervalType(X[i].hi))
    end
  end
end
