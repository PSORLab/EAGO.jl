
"""
    MC_impRelax(h::Function, hj::Function, p::Vector{SMCg{N,T}}, pmid::Vector{T},
                X::Vector{Interval{T}}, P::Vector{Interval{T}},mc_opts::mc_opts{T},
                param::Vector{Vector{SMCg{N,T}}})

Relaxes the implicit function determined by `h(x,p)` with `x` in `X` and `p` in
`P`. The reference point for the affine relaxations is `pmid`. The parameters
generated from the relaxation at `pmid` are `param` and the basic parameters of the
fixed point method are `mc_opt`.
"""
function MC_impRelax(h::Function, hj::Function, p::Vector{SMCg{N,V,T}}, pmid::Vector{T},
                     X::Vector{V}, P::Vector{V},
                     mc_opts::mc_opts{T},param::Vector{Vector{SMCg{N,V,T}}}) where {N,V,T<:AbstractFloat}

    nx::Int64 = mc_opts.nx
    np::Int64 = mc_opts.np
    szero::SVector{np,T} = zeros(SVector{np,T})
    sone::SVector{np,T} = ones(SVector{np,T})

    x_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),true) for i=1:nx]
    xa_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),true) for i=1:nx]
    xA_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),true) for i=1:nx]
    z_mc::Vector{SMCg{np,V,T}} = mc_opts.lambda*xa_mc+(one(T)-mc_opts.lambda)*xA_mc

    p_mc::Vector{SMCg{np,V,T}} = copy(p)
    pref_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(pmid[i],pmid[i],sone,sone,V(P[i].lo,P[i].hi),false) for i=1:np]
    aff_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,Float64}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,V(xa_mc[i].cv,xA_mc[i].cc),false) for i=1:nx]

    # Begins loop to generate parameters
    for k=1:mc_opts.kmax
      Affine_Exp!(param[k],p_mc,pmid,xa_mc,xA_mc,z_mc,mc_opts)
      aff_mc = SMCg{np,V,Float64}[SMCg{np,V,Float64}(x_mc[i].cc,x_mc[i].cv,szero,szero,
                     V(x_mc[i].cv,x_mc[i].cc),false) for i=1:nx]
      PSMCg_Kernel!(h,hj,z_mc,aff_mc,p_mc,x_mc,mc_opts)

    end
    return x_mc
end
