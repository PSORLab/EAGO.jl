"""
    implicit_relax_h

Relaxes the implicit function determined by `h(x,p)` with `x` in `X` and `p` in
`P`. The reference point for the affine relaxations is `pmid`. The parameters
generated from the relaxation at `pmid` are `param` and the basic parameters of the
fixed point method are `mc_opt`.
"""
function implicit_relax_h(h::Function, hj::Function, p::Vector{MC{N,T}}, pmid::Vector{Float64},
                     X::Vector{Interval{Float64}}, P::Vector{Interval{Float64}},
                     mc_opts::mc_opts, param::Vector{Vector{MC{N,T}}}) where {N, T<:RelaxTag}

    nx::Int = mc_opts.nx
    np::Int = mc_opts.np
    lambda::Float64 = mc_opts.lambda
    szero::SVector{np,Float64} = zeros(SVector{np,Float64})
    sone::SVector{np,Float64} = ones(SVector{np,Float64})

    x_mc::Vector{MC{np,T}} = [MC{np,T}(X[i].lo,X[i].hi,Interval{Float64}(X[i].lo,X[i].hi),szero,szero,true) for i=1:nx]
    xa_mc::Vector{MC{np,T}} = [MC{np,T}(X[i].lo,X[i].lo,Interval{Float64}(X[i].lo,X[i].lo),szero,szero,true) for i=1:nx]
    xA_mc::Vector{MC{np,T}} = [MC{np,T}(X[i].hi,X[i].hi,Interval{Float64}(X[i].hi,X[i].hi),szero,szero,true) for i=1:nx]
    z_mc::Vector{MC{np,T}} = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc

    p_mc::Vector{MC{np,T}} = copy(p)
    pref_mc::Vector{MC{np,T}} = [MC{np,T}(pmid[i],pmid[i],Interval{Float64}(P[i].lo,P[i].hi),sone,sone,false) for i=1:np]
    aff_mc::Vector{MC{np,T}} = [MC{np,T}(xa_mc[i].cv,xA_mc[i].cc,Interval{Float64}(xa_mc[i].cv,xA_mc[i].cc),szero,szero,false) for i=1:nx]

    # Begins loop to generate parameters
    for k=1:mc_opts.kmax
      affine_exp!(param[k], p_mc, pref_mc, xa_mc, xA_mc, z_mc, nx, lambda)
      aff_mc = MC{np}[MC{np,T}(x_mc[i].cv,x_mc[i].cc,Interval{Float64}(x_mc[i].cv,x_mc[i].cc),
                             szero,szero,false) for i=1:nx]
      pmc_kernel_std!(h,hj,z_mc,aff_mc,p_mc,x_mc,mc_opts)
    end
    return x_mc
end

function implicit_relax_h!(h::Function, hj::Function, p_mc::Vector{MC{N,T}}, pref_mc::Vector{MC{N,T}},
                     xp_mc::Vector{MC{N,T}}, x_mc, xa_mc::Vector{MC{N,T}}, xA_mc::Vector{MC{N,T}}, z_mc::Vector{MC{N,T}}, aff_mc::Vector{MC{N,T}},
                     X::Vector{Interval{Float64}}, P::Vector{Interval{Float64}}, mc_opts::mc_opts, param,
                     H, J, Y, interval_bnds::Bool, flt_param::Vector{Float64}, precond::Bool; subgrad_cntr::Bool = false) where {N, T<:RelaxTag}

    nx::Int = mc_opts.nx
    kmax::Int = mc_opts.kmax
    lambda::Float64 = mc_opts.lambda
    cntr::Symbol = mc_opts.contractor_type

    for i in 1:nx
      if interval_bnds
        x_mc[i] = MC{N,T}(X[i].lo, X[i].hi)
        xa_mc[i] = MC{N,T}(X[i].lo, X[i].lo)
        xA_mc[i] = MC{N,T}(X[i].hi, X[i].hi)
      end
    end
    z_mc[:] = lambda*xa_mc[:] + (1.0 - lambda)*xA_mc[:]

    # Begins loop to generate parameters
    #subgrad_cntr && set_reference!(cv.(pref_mc),P,subgrad_cntr)
    for k=1:(kmax)
      affine_exp!(param[:,k], p_mc, pref_mc, xa_mc, xA_mc, z_mc, nx, lambda)
      pmc_kernel!(h, hj, H, J, Y, z_mc, aff_mc, p_mc, x_mc, xa_mc, xA_mc,
                 cntr, nx, xp_mc, flt_param, precond)
    end
    #subgrad_cntr && set_reference!(cv.(pref_mc), P, false)
end

"""
    implicit_relax_f

Relaxes the function `f(x,p)` by relaxation the state variable `x` using the implicit
function determined by `h(x,p)` with `x` in `X` and `p` in `P`. The reference
point for the affine relaxations is `pmid`. The parameters generated from the
relaxation at `pmid` are `param` and the basic parameters of the fixed point
method are `mc_opt`.
"""
function implicit_relax_f(f::Function,h::Function,hj::Function,X::Vector{Interval{Float64}},
                    P::Vector{Interval{Float64}},p::Vector{Float64},pmid::Vector{Float64},
                    mc_opt::mc_opts,param::Vector{Vector{MC{N,T}}}) where {N, T<:RelaxTag}
  np::Int = length(P)
  sone::SVector{np,Float64} = ones(SVector{np,Float64})
  p_mc::Vector{MC{np,T}} = [MC{np,T}(p[i],p[i],Interval{Float64}(P[i].lo,P[i].hi),sone,sone,false) for i=1:np]
  xpMC::Vector{MC{np,T}} = implicit_relax_h(h,hj,p_mc,pmid,X,P,mc_opt,param)
  return f(xpMC,p_mc)
end

"""
    implicit_relax_fg

Relaxes the functions `f(x,p)` and `g(x,p)` by relaxation the state variable `x`
using the implicit function determined by `h(x,p)` with `x` in `X` and `p` in `P`.
The reference point for the affine relaxations is `pmid`. The parameters generated
from the relaxation at `pmid` are `param` and the basic parameters of the fixed point
method are `mc_opt`.
"""
function implicit_relax_fg(f::Function,g::Function,h::Function,hj::Function,
                     X::Vector{Interval{Float64}},P::Vector{Interval{Float64}},
                     p::Vector{Float64},pmid::Vector{Float64},
                     mc_opt::mc_opts,param::Vector{Vector{MC{N,T}}}) where {N, T<:RelaxTag}
  np::Int = length(P)
  sone::SVector{np,Float64} = ones(SVector{np,Float64})
  p_mc::Vector{MC{np,T}} = [MC{np,T}(p[i],p[i],Interval{Float64}(P[i].lo,P[i].hi),sone,sone,false) for i=1:np]
  xpMC::Vector{MC{np,T}} = implicit_relax_h(h,hj,p_mc,pmid,X,P,mc_opt,param)
  return f(xpMC,p_mc),g(xpMC,p_mc)
end
