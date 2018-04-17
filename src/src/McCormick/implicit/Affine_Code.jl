"""
    Affine_Exp!(x::Vector{SMCg{N,T}},p::Vector{SMCg{N,T}},p_ref::Vector{SMCg{N,T}},
               xa::Vector{SMCg{N,T}},xA::Vector{SMCg{N,T}},z::Vector{SMCg{N,T}},
               opt::Array{Any})

Computates the affine relaxations of the state variable. Inputs are:
* `x::Vector{SMCg{N,T}}`: State variable relaxation
* `p::Vector{SMCg{N,T}}`: Decision variable relaxation
* `p_ref::Vector{SMCg{N,T}}`: Reference variable relaxation
* `xa::Vector{SMCg{N,T}}`: Lower affine relaxation of the state variable
* `xA::Vector{SMCg{N,T}}`: Upper affine relaxation of the state variable
* `z::Vector{SMCg{N,T}}`: Affine function in `X`
* `opt::Array{Any,1}`: `[np,nx,lambda]` values for relaxation
Returns the tuple `(xa,xA,z)`:
* `xa::Vector{SMCg{N,T}}`: Lower affine relaxation of the state variable
* `xA::Vector{SMCg{N,T}}`: Upper affine relaxation of the state variable
* `z::Vector{SMCg{N,T}}`: Affine function in X
--------------------------------------------------------------------------------
"""
function Affine_Exp!(x::Vector{SMCg{N,V,T}}, p::Vector{SMCg{N,V,T}}, p_ref,
                     xa::Vector{SMCg{N,V,T}}, xA::Vector{SMCg{N,V,T}}, z::Vector{SMCg{N,V,T}},
                     mc_opts::mc_opts{T}) where {N,V,T<:AbstractFloat}

  nx::Int64 = mc_opts.nx
  lambda::T = mc_opts.lambda
  S1::SMCg{N,V,T},S2::SMCg{N,V,T},S3::SMCg{N,V,T} = zero(x[1]),zero(x[1]),zero(x[1])
  for i = 1:nx
   S1 = zero(x[1])
   S2 = zero(x[1])
   S3 = zero(x[1])
   for j = 1:N
      S1 = S1 + (p[j]-p_ref[j])*x[i].cv_grad[j]
      S2 = S2 + (p[j]-p_ref[j])*x[i].cc_grad[j]
      S3 = S3 + (lambda*x[i].cv_grad[j]+(one(T)-lambda)*x[i].cc_grad[j])*(p[j]-p_ref[j])
   end
   temp1::SMCg{N,V,T} = x[i].cv + S1
   temp2::SMCg{N,V,T} = x[i].cc + S2
   temp3::SMCg{N,V,T} = lambda*x[i].cv+(one(T)-lambda)*x[i].cc+S3
   xa[i] = SMCg{N,V,T}(temp1.cc,temp1.cv,x[i].cv_grad,x[i].cv_grad,temp1.Intv,S1.cnst,x[i].IntvBox,x[i].xref)
   xA[i] = SMCg{N,V,T}(temp2.cc,temp2.cv,x[i].cc_grad,x[i].cc_grad,temp2.Intv,S2.cnst,x[i].IntvBox,x[i].xref)
   z[i] = SMCg{N,V,T}(temp3.cc,temp3.cv,
                    lambda*x[i].cv_grad+(one(T)-lambda)*x[i].cc_grad,
                    lambda*x[i].cv_grad+(one(T)-lambda)*x[i].cc_grad,
                    temp3.Intv,S3.cnst,x[i].IntvBox,x[i].xref)
  end
end

"""
    Correct_Exp!(z_mc::Vector{SMCg{N,T}},x_mc::Vector{SMCg{N,T}},xL::Vector{Float64},
                 xU::Vector{Float64},nx::Int64,np::Int64,epsv::Float64)

Corrects the relaxation of the state variable `x_mc` if the affine relaxation,
'z_mc', exceeds the interval bounds `xL` or `xU`.
* `z_mc::Vector{SMCg{N,T}}`: Affine relaxation
* `x_mc::Vector{SMCg{N,T}}`: Relaxation of state variable
* `xL::Vector{Float64}`: Lower bound on state vector
* `xU::Vector{Float64}`: Upper bound on state vector
* `nx::Int64`: Size of the state vector
* `np::Int64`: Size of the decision vector
* `epsv::Float64`: Tolerance for checking that subgradient exceeds bound
"""
function Correct_Exp!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                     X::Vector{V},nx::Int64,np::Int64,
                     epsv::Float64) where {N,V,T<:AbstractFloat}
  zero_grad::SVector{N,T} = @SVector zeros(T,N)
  for i = 1:nx
    if (z_mc[i].Intv.lo-epsv < X[i].lo)
      x_mc[i] = SMCg{N,V,T}(x_mc[i].cc,X[i].lo,x_mc[i].cc_grad,zero_grad,x_mc[i].Intv,true,x_mc[i].IntvBox,x_mc[i].xref)
    end
    if (z_mc[i].Intv.hi+epsv > X[i].hi)
      x_mc[i] = SMCg{N,V,T}(X[i].hi,x_mc[i].cv,zero_grad,x_mc[i].cv_grad,x_mc[i].Intv,true,x_mc[i].IntvBox,x_mc[i].xref)
    end
  end
end
