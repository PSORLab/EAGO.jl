"""
    impRelax_f(f::Function,h::Function,hj::Function,X::Vector{Interval{T}},
               P::Vector{Interval{T}},p::Vector{T},pmid::Vector{T},
               mc_opt::mc_opts{T},param::Vector{Vector{SMCg{N,T}}})

Relaxes the function `f(x,p)` by relaxation the state variable `x` using the implicit
function determined by `h(x,p)` with `x` in `X` and `p` in `P`. The reference
point for the affine relaxations is `pmid`. The parameters generated from the
relaxation at `pmid` are `param` and the basic parameters of the fixed point
method are `mc_opt`.
"""
function impRelax_f(f::Function,h::Function,hj::Function,X::Vector{V},
                    P::Vector{V},p::Vector{T},pmid::Vector{T},
                    mc_opt::mc_opts{T},param::Vector{Vector{SMCg{N,V,T}}}) where {N,V,T<:AbstractFloat}
  np::Int64 = length(P)
  sone::SVector{np,T} = ones(SVector{np,T})
  p_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(p[i],p[i],sone,sone,V(P[i].lo,P[i].hi),false) for i=1:np]
  xpMC::Vector{SMCg{np,V,T}} = MC_impRelax(h,hj,p_mc,pmid,X,P,mc_opt,param)
  return f(xpMC,p_mc)
end

"""
    impRelax_fg(f::Function,g::Function,h::Function,hj::Function,
                X::Vector{Interval{T}},P::Vector{Interval{T}},
                p::Vector{T},pmid::Vector{T},
                mc_opt::mc_opts,param::Vector{Vector{SMCg{N,T}}})

Relaxes the functions `f(x,p)` and `g(x,p)` by relaxation the state variable `x`
using the implicit function determined by `h(x,p)` with `x` in `X` and `p` in `P`.
The reference point for the affine relaxations is `pmid`. The parameters generated
from the relaxation at `pmid` are `param` and the basic parameters of the fixed point
method are `mc_opt`.
"""
function impRelax_fg(f::Function,g::Function,h::Function,hj::Function,
                     X::Vector{V},P::Vector{V},
                     p::Vector{T},pmid::Vector{T},
                     mc_opt::mc_opts{T},param::Vector{Vector{SMCg{N,V,T}}}) where {N,V,T<:AbstractFloat}
  np::Int64 = length(P)
  sone::SVector{np,T} = ones(SVector{np,T})
  p_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(p[i],p[i],sone,sone,V(P[i].lo,P[i].hi),false) for i=1:np]
  xpMC::Vector{SMCg{np,V,T}} = MC_impRelax(h,hj,p_mc,pmid,X,P,mc_opt,param)
  return f(xpMC,p_mc),g(xpMC,p_mc)
end
