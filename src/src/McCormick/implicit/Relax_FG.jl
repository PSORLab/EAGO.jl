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
  sone::SVector{np,T} = @SVector ones(np)
  p_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(p[i],p[i],sone,sone,V(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
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
  sone::SVector{np,T} = @SVector ones(np)
  p_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(p[i],p[i],sone,sone,V(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  xpMC::Vector{SMCg{np,V,T}} = MC_impRelax(h,hj,p_mc,pmid,X,P,mc_opt,param)
  return f(xpMC,p_mc),g(xpMC,p_mc)
end

"""
    NimpRelax_f(f::Function,h::Function,hj::Function,X::Vector{Interval{T}},
               P::Vector{Interval{T}},p::Vector{T},pmid::Vector{T},
               mc_opt::mc_opts{T},param::Vector{Vector{SMCg{N,T}}})

Relaxes the function `f(x,p)` by relaxation the state variable `x` using the implicit
function determined by `h(x,p)` with `x` in `X` and `p` in `P`. The reference
point for the affine relaxations is `pmid`. The parameters generated from the
relaxation at `pmid` are `param` and the basic parameters of the fixed point
method are `mc_opt`. Preconditioning is using a sparse LDU with full pivoting approach,
the jacobian with respect to X and the preconditioner are stored in a sparse format,
sparse contractors are applied, and h & hj are evaluated in place.
"""
function NimpRelax_f(f::Function,h!::Function,hj!::Function,X::Vector{Interval{T}},
                    P::Vector{Interval{T}},p::Vector{T},pmid::Vector{T},
                    mc_opt::mc_opts{T},param::Vector{Vector{SMCg{N,V,T}}}) where {N,V,T<:AbstractFloat}
  np::Int64 = length(P)
  sone::SVector{np,V,T} = @SVector ones(np)
  p_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(p[i],p[i],sone,sone,V(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  xpMC::Vector{SMCg{np,V,T}} = MC_NimpRelax(h!,hj!,p_mc,pmid,X,P,mc_opt,param)
  return f(xpMC,p_mc)
end

"""
    NimpRelax_fg(f::Function,g::Function,h::Function,hj::Function,
                X::Vector{Interval{T}},P::Vector{Interval{T}},
                p::Vector{T},pmid::Vector{T},
                mc_opt::mc_opts,param::Vector{Vector{SMCg{N,T}}})

Relaxes the functions `f(x,p)` and `g(x,p)` by relaxation the state variable `x`
using the implicit function determined by `h(x,p)` with `x` in `X` and `p` in `P`.
The reference point for the affine relaxations is `pmid`. The parameters generated
from the relaxation at `pmid` are `param` and the basic parameters of the fixed point
method are `mc_opt`. Preconditioning is using a sparse LDU with full pivoting approach,
the jacobian with respect to X and the preconditioner are stored in a sparse format,
sparse contractors are applied, and h & hj are evaluated in place.
"""
function NimpRelax_fg(f::Function,g::Function,h!::Function,hj!::Function,
                     X::Vector{V},P::Vector{V},
                     p::Vector{T},pmid::Vector{T},
                     mc_opt::mc_opts{T},param::Vector{Vector{SMCg{N,V,T}}}) where {N,V,T<:AbstractFloat}
  np::Int64 = length(P)
  sone::SVector{np,T} = @SVector ones(np)
  p_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(p[i],p[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  xpMC::Vector{SMCg{np,V,T}} = MC_NimpRelax(h!,hj!,p_mc,pmid,X,P,mc_opt,param)
  return f(xpMC,p_mc),g(xpMC,p_mc)
end

"""
    NdimpRelax_f(f::Function,h::Function,hj::Function,X::Vector{Interval{T}},
               P::Vector{Interval{T}},p::Vector{T},pmid::Vector{T},
               mc_opt::mc_opts{T},param::Vector{Vector{SMCg{N,T}}})

Relaxes the function `f(x,p)` by relaxation the state variable `x` using the implicit
function determined by `h(x,p)` with `x` in `X` and `p` in `P`. The reference
point for the affine relaxations is `pmid`. The parameters generated from the
relaxation at `pmid` are `param` and the basic parameters of the fixed point
method are `mc_opt`. Preconditioning is using a LU with full pivoting approach,
the jacobian with respect to X and the preconditioner are stored, and h & hj are evaluated in place.
"""
function NdimpRelax_f(f::Function,h!::Function,hj!::Function,X::Vector{V},
                    P::Vector{V},p::Vector{T},pmid::Vector{T},
                    mc_opt::mc_opts{T},param::Vector{Vector{SMCg{N,V,T}}}) where {N,V,T<:AbstractFloat}
  np::Int64 = length(P)
  sone::SVector{np,T} = @SVector ones(np)
  p_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(p[i],p[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  xpMC::Vector{SMCg{np,V,T}} = MC_NdimpRelax(h!,hj!,p_mc,pmid,X,P,mc_opt,param)
  return f(xpMC,p_mc)
end

"""
    NdimpRelax_fg(f::Function,g::Function,h::Function,hj::Function,
                X::Vector{Interval{T}},P::Vector{Interval{T}},
                p::Vector{T},pmid::Vector{T},
                mc_opt::mc_opts,param::Vector{Vector{SMCg{N,T}}})

Relaxes the functions `f(x,p)` and `g(x,p)` by relaxation the state variable `x`
using the implicit function determined by `h(x,p)` with `x` in `X` and `p` in `P`.
The reference point for the affine relaxations is `pmid`. The parameters generated from the
relaxation at `pmid` are `param` and the basic parameters of the fixed point
method are `mc_opt`. Preconditioning is using a LU with full pivoting approach,
the jacobian with respect to X and the preconditioner are stored, and h & hj are evaluated in place.
"""
function NdimpRelax_fg(f::Function,g::Function,h!::Function,hj!::Function,
                     X::Vector{V},P::Vector{V},
                     p::Vector{T},pmid::Vector{T},
                     mc_opt::mc_opts{T},param::Vector{Vector{SMCg{N,V,T}}}) where {N,V,T<:AbstractFloat}
  np::Int64 = length(P)
  sone::SVector{np,T} = @SVector ones(np)
  p_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(p[i],p[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  xpMC::Vector{SMCg{np,V,T}} = MC_NdimpRelax(h!,hj!,p_mc,pmid,X,P,mc_opt,param)
  return f(xpMC,p_mc),g(xpMC,p_mc)
end
