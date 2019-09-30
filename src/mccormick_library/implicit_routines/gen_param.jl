function pmc_kernel_std!(h::Function,hj::Function,z_mc::Vector{MC{N,T}},
                       aff_mc::Vector{MC{N,T}},p_mc::Vector{MC{N,T}},
                       x_mc::Vector{MC{N,T}},opt::mc_opts) where {N, T<:RelaxTag}
  H,J = precondition_mc(h,hj,z_mc,aff_mc,p_mc,opt)
  if (opt.contractor_type == :Newton)
      mc_dense_newton_gs!(z_mc,x_mc,J,H,opt.nx)
  elseif (opt.contractor_type == :Krawczyk)
      mc_dense_krawczyk_cw!(z_mc,x_mc,J,H,opt.nx)
  else
      error("The contractor type $(opt.contractor_type) is not currently supported. The
             contractors :Newton and :Krawczyk are currently supported.")
  end
  return
end

"""
    pmc_kernel!

Peforms the following steps in sequence:
- Evaluates the function h!(H, x, xp, p, t) in place with x = z_mc, p = p_mc,
xp = xp_mc, t = flt_param and preconditions H using an interval midpoint
preconditioner if precond = true.
- Evaluates the function hj!(J, x, xp, p, t) in place with x = aff_mc, p = p_mc,
xp = xp_mc, t = flt_param and preconditions J using an interval midpoint
preconditioner if precond = true.
- Lastly, applies a Newton-type contractor method. The parametric GS Newton
contractor if cntr = :Newton and the componentwise Krawczyk contractor otherwise.
"""
function pmc_kernel!(h!::Function, hj!::Function, H, J,
                     Y::Array{Float64,2}, z_mc::Vector{MC{N,T}}, aff_mc::Vector{MC{N,T}},
                     p_mc::Vector{MC{N,T}}, x_mc, xa_mc, xA_mc, cntr::Symbol,
                     nx::Int, xp_mc::Vector{MC{N,T}},
                     flt_param::Vector{Float64}, precond::Bool) where {N, T<:RelaxTag}

  for i in 1:nx
    aff_mc[i] = MC{N,T}(xa_mc[i].cv, xA_mc[i].cc)
  end

  if nx == 1
    mc_precondition_1!(h!, hj!, z_mc, aff_mc, p_mc, J, H, Y, xp_mc, flt_param, precond)
  else
    mc_denseband_precondition!(h!, hj!, z_mc, aff_mc, p_mc, J, H, Y, xp_mc, flt_param, precond, nx)
  end

  if (cntr == :Newton)
      mc_dense_newton_gs!(z_mc, x_mc, J, H, nx)
  elseif (cntr == :Krawczyk)
      mc_dense_krawczyk_cw!(z_mc, x_mc, J, H, nx)
  else
      error("The contractor type $(cntr) is not currently supported. The
             contractors :Newton and :Krawczyk are currently supported.")
  end
  return
end

function gen_expansion_params(h::Function, hj::Function,
                      X::Vector{Interval{Float64}},
                      P::Vector{Interval{Float64}},
                      pmid::Vector{Float64},
                      mc_opts::mc_opts,
                      T::RelaxTag)

  nx::Int = length(X)
  np::Int = length(P)
  lambda::Float64 = mc_opts.lambda
  szero::SVector{np,Float64} = zeros(SVector{np,Float64})
  sone::SVector{np,Float64} = ones(SVector{np,Float64})

  x_mc::Vector{MC{np,T}} = MC{np}[MC{np,T}(X[i].lo,X[i].hi,Interval{Float64}(X[i].lo,X[i].hi),szero,szero,false) for i=1:nx]
  xa_mc::Vector{MC{np,T}} = MC{np}[MC{np,T}(X[i].lo,X[i].lo,Interval{Float64}(X[i].lo,X[i].lo),szero,szero,false) for i=1:nx]
  xA_mc::Vector{MC{np,T}} = MC{np}[MC{np,T}(X[i].hi,X[i].hi,Interval{Float64}(X[i].hi,X[i].hi),szero,szero,false) for i=1:nx]
  z_mct::Vector{MC{np,T}} = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc
  z_mc = rnd_out_z_all(z_mct,mc_opts.aff_correct_eps)

  p_mc::Vector{MC{np,T}} = MC{np}[MC{np,T}(pmid[i],pmid[i],Interval{Float64}(P[i].lo,P[i].hi),sone,sone,false) for i=1:np]
  pref_mc::Vector{MC{np,T}} = copy(p_mc)
  aff_mc::Vector{MC{np,T}} = MC{np}[MC{np,T}(xa_mc[i].cv,xA_mc[i].cc,Interval{Float64}(xa_mc[i].Intv.lo,xA_mc[i].Intv.hi),szero,szero,false) for i=1:nx]
  sto_out::Vector{Vector{MC{np,T}}} = Vector{MC{np}}[x_mc for j=1:(mc_opts.kmax+1)]
  sto_out[1] = copy(x_mc)
  optc = Any[szero,sone]
  for k=1:mc_opts.kmax
    pmc_kernel_std!(h,hj,z_mc,aff_mc,p_mc,x_mc,mc_opts)
    affine_exp!(x_mc,p_mc,p_mc,xa_mc,xA_mc,z_mc,nx,lambda)
    z_mc = rnd_out_z_all(z_mc,mc_opts.aff_correct_eps)
    correct_exp!(xa_mc,xA_mc,z_mc,x_mc,X,nx,mc_opts.aff_correct_eps)
    aff_mc = MC{np}[MC{np}(xa_mc[i].cv,xA_mc[i].cc,Interval{Float64}(xA_mc[i].Intv.lo,
                           xA_mc[i].Intv.hi),xa_mc[i].cv_grad,xA_mc[i].cc_grad,false) for i=1:nx]
    # store relaxation
    sto_out[k+1] = copy(x_mc)
    end
  return sto_out
end

function gen_expansion_params!(h!::Function, hj!::Function, pref_mc::Vector{MC{N,T}}, xp_mc::Vector{MC{N,T}},
                               x_mc, xa_mc::Vector{MC{N,T}}, xA_mc::Vector{MC{N,T}}, z_mc::Vector{MC{N,T}},
                               aff_mc::Vector{MC{N,T}}, X::Vector{Interval{Float64}}, P::Vector{Interval{Float64}},
                               opts::mc_opts, sto_out, H::Vector{MC{N,T}},
                               J::Array{MC{N,T},2}, Y::Array{Float64,2}, interval_bnds::Bool,
                               flt_param::Vector{Float64}, precond::Bool; subgrad_cntr::Bool = false) where {N, T<:RelaxTag}

  nx::Int = length(X)
  kmax::Int = opts.kmax
  lambda::Float64 = opts.lambda
  aff_eps::Float64 = opts.aff_correct_eps
  cntr::Symbol = opts.contractor_type

  for i in 1:nx
    if interval_bnds
      x_mc[i] = MC{N,T}(X[i].lo, X[i].hi)
      xa_mc[i] = MC{N,T}(X[i].lo, X[i].lo)
      xA_mc[i] = MC{N,T}(X[i].hi, X[i].hi)
    end
  end
  z_mc[:] = lambda*xa_mc[:] + (1.0 - lambda)*xA_mc[:]
  sto_out[:,1:1] .= x_mc

  subgrad_cntr && set_reference!(cv.(pref_mc),P,subgrad_cntr)
  for k=1:(kmax-1)
    pmc_kernel!(h!, hj!, H, J, Y, z_mc, aff_mc, pref_mc, x_mc, xa_mc, xA_mc,
                cntr, nx, xp_mc, flt_param, precond)
    affine_exp!(x_mc, pref_mc, pref_mc, xa_mc, xA_mc, z_mc, nx, lambda)
    correct_exp!(xa_mc, xA_mc, z_mc, x_mc, X, nx, aff_eps)
    sto_out[:,(k+1):(k+1)] .= x_mc
  end
  subgrad_cntr && set_reference!(cv.(pref_mc), P, false)
  return
end
