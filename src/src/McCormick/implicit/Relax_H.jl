
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

    SP = SVector{np,V}(P)
    SP0 = SVector{np,T}(pmid)

    x_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),true,SP,SP0) for i=1:nx]
    xa_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),true,SP,SP0) for i=1:nx]
    xA_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),true,SP,SP0) for i=1:nx]
    z_mc::Vector{SMCg{np,V,T}} = mc_opts.lambda*xa_mc+(one(T)-mc_opts.lambda)*xA_mc

    p_mc::Vector{SMCg{np,V,T}} = copy(p)
    pref_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(pmid[i],pmid[i],sone,sone,V(P[i].lo,P[i].hi),false,SP,SP0) for i=1:np]
    aff_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,Float64}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,V(xa_mc[i].cv,xA_mc[i].cc),false,SP,SP0) for i=1:nx]

    # Begins loop to generate parameters
    for k=1:mc_opts.kmax
      Affine_Exp!(param[k],p_mc,pmid,xa_mc,xA_mc,z_mc,mc_opts)
      aff_mc = SMCg{np,V,Float64}[SMCg{np,V,Float64}(x_mc[i].cc,x_mc[i].cv,szero,szero,
                     V(x_mc[i].cv,x_mc[i].cc),false,SP,SP0) for i=1:nx]
      PSMCg_Kernel!(h,hj,z_mc,aff_mc,p_mc,x_mc,mc_opts)

    end
    return x_mc
end


#=
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
    exp_opt::Vector{Any} = Any[nx,np,mc_opts.lambda]

    x_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),true,[∅],[1.0]) for i=1:nx]
    xa_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),true,[∅],[1.0]) for i=1:nx]
    xA_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),true,[∅],[1.0]) for i=1:nx]
    z_mc::Vector{SMCg{np,V,T}} = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc

    if mc_opts.z_rnd_all == true
      z_mc::Vector{SMCg{np,T}} = Rnd_Out_Z_All(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_all_eps)
    elseif mc_opts.z_rnd == true
      z_mc = Rnd_Out_Z_Intv(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_eps)
    else
      z_mc = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc
    end

    p_mc::Vector{SMCg{np,V,T}} = copy(p)
    pref_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(pmid[i],pmid[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
    aff_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,Float64}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,@interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]

    H_mc::Vector{SMCg{np,V,T}} = h(z_mc,p_mc)
    dH_mc::Union{Vector{SMCg{np,V,T}},Array{SMCg{np,V,T},2}} = hj(aff_mc,p_mc)
    Y::Union{Vector{T},Array{T,2}} = (nx == 1) ? [one(T)/mid(dH_mc[1].Intv)] : mid.(Intv.(dH_mc))

    # stores some things for repeated us
    optc = Any[szero,sone]

    # Begins loop to generate parameters
    for k=1:mc_opts.kmax
      Affine_Exp!(param[k],p_mc,pref_mc,xa_mc,xA_mc,z_mc,exp_opt)

      # sets up affine relaxations and preconditioning
      aff_mc = [SMCg(xA_mc[i].cc,xa_mc[i].cv,xA_mc[i].cc_grad,xa_mc[i].cv_grad,
                     @interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
      H_mc,dH_mc = h(z_mc,p_mc),hj(aff_mc,p_mc)
      Y = (nx == 1) ? [one(T)/mid(dH_mc[1].Intv)] : mid.(Intv.(dH_mc))
      # applies preconditioner
      Precondition!(H_mc,dH_mc,Y,nx)

    end
    return x_mc
end
=#
