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
function Affine_Exp!(x::Vector{SMCg{N,V,T}}, p::Vector{SMCg{N,V,T}},
                     p_ref::Vector{SMCg{N,V,T}}, xa::Vector{SMCg{N,V,T}},
                     xA::Vector{SMCg{N,V,T}}, z::Vector{SMCg{N,V,T}},
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
   xa[i] = SMCg{N,V,T}(temp1.cc,temp1.cv,x[i].cv_grad,x[i].cv_grad,temp1.Intv,S1.cnst)
   xA[i] = SMCg{N,V,T}(temp2.cc,temp2.cv,x[i].cc_grad,x[i].cc_grad,temp2.Intv,S2.cnst)
   z[i] = SMCg{N,V,T}(temp3.cc,temp3.cv,
                    lambda*x[i].cv_grad+(one(T)-lambda)*x[i].cc_grad,
                    lambda*x[i].cv_grad+(one(T)-lambda)*x[i].cc_grad,
                    temp3.Intv,S3.cnst)
  end
end

#=
"""
    Affine_Exp!(x::Vector{SMCg{N,T}},p::Vector{SMCg{N,T}},p_ref::Vector{SMCg{N,T}},
               xa::Vector{SMCg{N,T}},xA::Vector{SMCg{N,T}},z::Vector{SMCg{N,T}},
               opt::Array{Any})

Computates the affine relaxations of the state variable. Inputs are:
* `x::Vector{HybridMC{N,V,T}}`: State variable relaxation
* `p::Vector{HybridMC{N,V,T}}`: Decision variable relaxation
* `p_ref::Vector{HybridMC{N,V,T}}`: Reference variable relaxation
* `xa::Vector{HybridMC{N,V,T}}`: Lower affine relaxation of the state variable
* `xA::Vector{HybridMC{N,V,T}}`: Upper affine relaxation of the state variable
* `z::Vector{HybridMC{N,V,T}}`: Affine function in `X`
* `opt::Array{Any,1}`: `[np,nx,lambda]` values for relaxation
--------------------------------------------------------------------------------
"""
function Affine_Exp!(x::Vector{HybridMC{N,V,T}}, p::Vector{HybridMC{N,V,T}},
                     p_ref::Vector{HybridMC{N,V,T}}, xa::Vector{HybridMC{N,V,T}},
                     xA::Vector{HybridMC{N,V,T}}, z::Vector{HybridMC{N,V,T}},
                     mc_opts::mc_opts{T}) where {N,V,T<:AbstractFloat}

  nx::Int64 = mc_opts.nx
  lambda::T = mc_opts.lambda
  S1::HybridMC{N,V,T},S2::HybridMC{N,V,T},S3::HybridMC{N,V,T} = zero(x[1]),zero(x[1]),zero(x[1])
  for i = 1:nx
   S1 = zero(x[1])
   S2 = zero(x[1])
   S3 = zero(x[1])
   for j = 1:N
      S1 = S1 + (p[j]-p_ref[j])*x[i].cv_grad[j]
      S2 = S2 + (p[j]-p_ref[j])*x[i].cc_grad[j]
      S3 = S3 + (lambda*x[i].cv_grad[j]+(one(T)-lambda)*x[i].cc_grad[j])*(p[j]-p_ref[j])
   end
   temp1::HybridMC{N,V,T} = x[i].cv + S1
   temp2::HybridMC{N,V,T} = x[i].cc + S2
   temp3::HybridMC{N,V,T} = lambda*x[i].cv+(one(T)-lambda)*x[i].cc+S3
   xa[i] = HybridMC{N,V,T}(SMCg{N,V,T}(cc(temp1),cv(temp1),cv_grad(x[i]),cv_grad(x[i]),Intv(temp1),S1.SMC.cnst))
   xA[i] = SMCg{N,V,T}(cc(temp2),cv(temp2),cc_grad(x[i]),cc_grad(x[i]),Intv(temp2),S2.SMC.cnst)
   z[i] = SMCg{N,V,T}(cc(temp3),cv(temp3),
                    lambda*cv_grad(x[i])+(one(T)-lambda)*cc_grad(x[i]),
                    lambda*cv_grad(x[i])+(one(T)-lambda)*cc_grad(x[i]),
                    Intv(temp3),S3.SMC.cnst)
  end
end
=#

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
    if (z_mc[i].Intv.lo-epsv < X[i].lo) && (z_mc[i].Intv.hi+epsv > X[i].hi)
      x_mc[i] = SMCg{N,V,T}(X[i].hi,X[i].lo,zero_grad,zero_grad,X[i],true)
    end
    if (z_mc[i].Intv.lo-epsv < X[i].lo)
      x_mc[i] = SMCg{N,V,T}(x_mc[i].cc,X[i].lo,x_mc[i].cc_grad,zero_grad,V(X[i].lo,x_mc[i].Intv.hi),true)
    end
    if (z_mc[i].Intv.hi+epsv > X[i].hi)
      x_mc[i] = SMCg{N,V,T}(X[i].hi,x_mc[i].cv,zero_grad,x_mc[i].cv_grad,V(x_mc[i].Intv.lo,X[i].hi),true)
    end
  end
end

#=
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
function Correct_Exp!(z_mc::Vector{HybridMC{N,V,T}},x_mc::Vector{HybridMC{N,V,T}},
                     X::Vector{V},nx::Int64,np::Int64,
                     epsv::Float64) where {N,V,T<:AbstractFloat}
  zero_grad::SVector{N,T} = @SVector zeros(T,N)
  for i = 1:nx
    if (Intv(z_mc[i]).lo-epsv < X[i].lo) && (Intv(z_mc[i]).hi+epsv > X[i].hi)
      x_mc[i] = HybridMC{N,V,T}(SMCg{N,V,T}(X[i].hi,X[i].lo,zero_grad,zero_grad,X[i],true))
    end
    if (Intv(z_mc[i]).lo-epsv < X[i].lo)
      x_mc[i] = HybridMC{N,V,T}(SMCg{N,V,T}(cc(x_mc[i]),X[i].lo,cc_grad(x_mc[i]),zero_grad,V(X[i].lo,Intv(x_mc[i]).hi),true))
    end
    if (Intv(z_mc[i]).hi+epsv > X[i].hi)
      x_mc[i] = HybridMC{N,V,T}(SMCg{N,V,T}(X[i].hi,cv(x_mc[i]),zero_grad,cv_grad(x_mc[i]),V(Intv(x_mc[i]).lo,X[i].hi),true))
    end
  end
end
=#
