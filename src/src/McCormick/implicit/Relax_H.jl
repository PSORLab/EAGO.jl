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

    nx::Int64 = length(X)
    np::Int64 = length(P)
    szero::SVector{np,T} = @SVector zeros(np)
    sone::SVector{np,T} = @SVector ones(np)
    exp_opt::Vector{Any} = Any[nx,np,mc_opts.lambda]

    x_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,@interval(X[i].lo,X[i].hi),true,[∅],[1.0]) for i=1:nx]
    xa_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,@interval(X[i].lo,X[i].lo),true,[∅],[1.0]) for i=1:nx]
    xA_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,@interval(X[i].hi,X[i].hi),true,[∅],[1.0]) for i=1:nx]
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
      if mc_opts.z_rnd_all == true
        z_mc = Rnd_Out_Z_All(z_mc,mc_opts.z_rnd_all_eps)
      elseif mc_opts.z_rnd == true
        z_mc = Rnd_Out_Z_Intv(z_mc,mc_opts.z_rnd_eps)
      end

      # sets up affine relaxations and preconditioning
      aff_mc = [SMCg(xA_mc[i].cc,xa_mc[i].cv,xA_mc[i].cc_grad,xa_mc[i].cv_grad,
                     @interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
      if mc_opts.aff_rnd_all == true
         aff_mc = Rnd_Out_Z_All(aff_mc,mc_opts.aff_rnd_all_eps)
      elseif mc_opts.aff_rnd == true
         aff_mc = Rnd_Out_Z_Intv(aff_mc,mc_opts.aff_rnd_eps)
      end

      if mc_opts.hhj_rnd_all == true
        H_mc,dH_mc = Rnd_Out_H_All(h(z_mc,p_mc),hj(aff_mc,p_mc),mc_opts.hhj_rnd_all_eps)
      elseif mc_opts.hhj_rnd == true
        H_mc,dH_mc = Rnd_Out_H_Intv(h(z_mc,p_mc),hj(aff_mc,p_mc),mc_opts.hhj_rnd_eps)
      else
        H_mc,dH_mc = h(z_mc,p_mc),hj(aff_mc,p_mc)
      end
      Y = (nx == 1) ? [one(T)/mid(dH_mc[1].Intv)] : mid.(Intv.(dH_mc))

      # applies preconditioner
      Precondition!(H_mc,dH_mc,Y,nx)

      if mc_opts.hhj_rnd_all == true
        H_mc,dH_mc = Rnd_Out_H_All(H_mc,dH_mc,mc_opts.hhj_rnd_all_eps)
      elseif mc_opts.hhj_rnd == true
        H_mc,dH_mc = Rnd_Out_H_Intv(H_mc,dH_mc,mc_opts.hhj_rnd_eps)
      end

      # applies parametric iteration
      if (mc_opts.style == "NewtonGS")
        MC_NewtonGS!(z_mc,x_mc,dH_mc,H_mc,nx)
      elseif (mc_opts.style == "KrawczykCW")
        MC_KrawczykCW!(z_mc,x_mc,dH_mc,H_mc,nx)
      else
        error("Unsupported Style of Implicit Relaxation")
      end
    end
    return x_mc
end

"""
    MC_NimpRelax(h::Function, hj::Function, p::Vector{SMCg{N,T}}, pmid::Vector{T},
                X::Vector{Interval{T}}, P::Vector{Interval{T}},mc_opts::mc_opts{T},
                param::Vector{Vector{SMCg{N,T}}})

Relaxes the implicit function determined by `h(x,p)` with `x` in `X` and `p` in
`P`. The reference point for the affine relaxations is `pmid`. The parameters
generated from the relaxation at `pmid` are `param` and the basic parameters of the
fixed point method are `mc_opt`. Preconditioning is using a sparse LDU with full pivoting approach,
the jacobian with respect to X and the preconditioner are stored in a sparse format,
sparse contractors are applied, and h & hj are evaluated in place.
"""
function MC_NimpRelax(h!::Function, hj!::Function, p::Vector{SMCg{N,V,T}}, pmid::Vector{T},
                     X::Vector{V}, P::Vector{V},
                     mc_opts::mc_opts{T},param::Vector{Vector{SMCg{N,V,T}}}) where {N,V,T<:AbstractFloat}

    nx::Int64 = length(X)
    np::Int64 = length(P)
    szero::SVector{np,T} = @SVector zeros(np)
    sone::SVector{np,T} = @SVector ones(np)
    exp_opt::Vector{Any} = Any[nx,np,mc_opts.lambda]

    x_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,@interval(X[i].lo,X[i].hi),true,[∅],[1.0]) for i=1:nx]
    xa_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,@interval(X[i].lo,X[i].lo),true,[∅],[1.0]) for i=1:nx]
    xA_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,@interval(X[i].hi,X[i].hi),true,[∅],[1.0]) for i=1:nx]
    if mc_opts.z_rnd_all == true
      z_mc::Vector{SMCg{np,V,T}} = Rnd_Out_Z_All(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_all_eps)
    elseif mc_opts.z_rnd == true
      z_mc = Rnd_Out_Z_Intv(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_eps)
    else
      z_mc = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc
    end

    p_mc::Vector{SMCg{np,V,T}} = copy(p)
    pref_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(pmid[i],pmid[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
    aff_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,Float64}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,@interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
    H_mc::Vector{SMCg{np,V,T}} = zeros(SMCg{np,T},nx)
    dH_mc::SparseMatrixCSC{SMCg{np,V,T},Int64} = spzeros(SMCg{np,T},nx,nx)
    Y::SparseMatrixCSC{T,Int64} = spzeros(T,nx,nx)

    # initalizes storage object for in-place sparse calculations
    SSto = SparseInSto()
    SSto.Xh = copy(H_mc)
    SSto.Yh = copy(H_mc)
    SSto.Zh = copy(H_mc)
    SSto.Xj = spzeros(SMCg{np,T},nx,nx)
    SSto.Yj = spzeros(SMCg{np,T},nx,nx)
    SSto.Zj = spzeros(SMCg{np,T},nx,nx)
    SSto.nx = copy(nx)

    # stores some things for repeated us
    optc = Any[szero,sone]

    # Begins loop to generate parameters
    for k=1:mc_opts.kmax
      Affine_Exp!(param[k],p_mc,pref_mc,xa_mc,xA_mc,z_mc,exp_opt)
      if mc_opts.z_rnd_all == true
        z_mc = Rnd_Out_Z_All(z_mc,mc_opts.z_rnd_all_eps)
      elseif mc_opts.z_rnd == true
        z_mc = Rnd_Out_Z_Intv(z_mc,mc_opts.z_rnd_eps)
      end

      # sets up affine relaxations and preconditioning
      aff_mc = [SMCg(xA_mc[i].cc,xa_mc[i].cv,xA_mc[i].cc_grad,xa_mc[i].cv_grad,
                     @interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
      if mc_opts.aff_rnd_all == true
         aff_mc = Rnd_Out_Z_All(aff_mc,mc_opts.aff_rnd_all_eps)
      elseif mc_opts.aff_rnd == true
         aff_mc = Rnd_Out_Z_Intv(aff_mc,mc_opts.aff_rnd_eps)
      end

      h!(H_mc,z_mc,p_mc),
      hj!(dH_mc,aff_mc,p_mc)
      Sparse_Precondition!(H_mc,dH_mc,mid.(Intv.(dH_mc)),SSto)
      dH_mc = transpose(dH_mc)

      if mc_opts.hhj_rnd_all == true
        H_mc,dH_mc = Rnd_Out_H_All(H_mc,dH_mc,mc_opts.hhj_rnd_all_eps)
      elseif mc_opts.hhj_rnd == true
        H_mc,dH_mc = Rnd_Out_H_Intv(H_mc,dH_mc,mc_opts.hhj_rnd_eps)
      end

      # applies parametric iteration
      if (mc_opts.style == "NewtonGS")
        MCn_NewtonGS!(z_mc,x_mc,dH_mc,H_mc,nx)
      elseif (mc_opts.style == "KrawczykCW")
        MCn_KrawczykCW!(z_mc,x_mc,dH_mc,H_mc,nx)
      else
        error("Unsupported Style of Implicit Relaxation")
      end
      dH_mc = transpose(dH_mc)
    end
    return x_mc
end

"""
    MC_NdimpRelax(h::Function, hj::Function, p::Vector{SMCg{N,T}}, pmid::Vector{T},
                X::Vector{Interval{T}}, P::Vector{Interval{T}},mc_opts::mc_opts{T},
                param::Vector{Vector{SMCg{N,T}}})

Relaxes the implicit function determined by `h(x,p)` with `x` in `X` and `p` in
`P`. The reference point for the affine relaxations is `pmid`. The parameters
generated from the relaxation at `pmid` are `param` and the basic parameters of the
fixed point method are `mc_opt`. Preconditioning is using a sparse LU with full pivoting approach,
the jacobian with respect to X, and h & hj are evaluated in place.
"""
function MC_NdimpRelax(h!::Function, hj!::Function, p::Vector{SMCg{N,V,T}}, pmid::Vector{T},
                     X::Vector{V}, P::Vector{V},
                     mc_opts::mc_opts{T},param::Vector{Vector{SMCg{N,V,T}}}) where {N,V,T<:AbstractFloat}

    nx::Int64 = length(X)
    np::Int64 = length(P)
    szero::SVector{np,T} = @SVector zeros(np)
    sone::SVector{np,T} = @SVector ones(np)
    exp_opt::Vector{Any} = Any[nx,np,mc_opts.lambda]

    x_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),true,[∅],[1.0]) for i=1:nx]
    xa_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),true,[∅],[1.0]) for i=1:nx]
    xA_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),true,[∅],[1.0]) for i=1:nx]
    if mc_opts.z_rnd_all == true
      z_mc::Vector{SMCg{np,V,T}} = Rnd_Out_Z_All(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_all_eps)
    elseif mc_opts.z_rnd == true
      z_mc = Rnd_Out_Z_Intv(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_eps)
    else
      z_mc = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc
    end

    p_mc::Vector{SMCg{np,V,T}} = copy(p)
    pref_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(pmid[i],pmid[i],sone,sone,V(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
    aff_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,Float64}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,V(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
    H_mc::Vector{SMCg{np,V,T}} = zeros(SMCg{np,T},nx)
    dH_mc::Array{SMCg{np,V,T},2} = zeros(SMCg{np,T},nx,nx)
    Y::Array{T,2} = zeros(T,nx,nx)

    # stores some things for repeated us
    optc = Any[szero,sone]

    # Begins loop to generate parameters
    for k=1:mc_opts.kmax
      Affine_Exp!(param[k],p_mc,pref_mc,xa_mc,xA_mc,z_mc,exp_opt)
      if mc_opts.z_rnd_all == true
        z_mc = Rnd_Out_Z_All(z_mc,mc_opts.z_rnd_all_eps)
      elseif mc_opts.z_rnd == true
        z_mc = Rnd_Out_Z_Intv(z_mc,mc_opts.z_rnd_eps)
      end

      # sets up affine relaxations and preconditioning
      aff_mc = [SMCg(xA_mc[i].cc,xa_mc[i].cv,xA_mc[i].cc_grad,xa_mc[i].cv_grad,
                     V(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
      if mc_opts.aff_rnd_all == true
         aff_mc = Rnd_Out_Z_All(aff_mc,mc_opts.aff_rnd_all_eps)
      elseif mc_opts.aff_rnd == true
         aff_mc = Rnd_Out_Z_Intv(aff_mc,mc_opts.aff_rnd_eps)
      end

      h!(H_mc,z_mc,p_mc),
      hj!(dH_mc,aff_mc,p_mc)
      Dense_Precondition!(H_mc,dH_mc,mid.(Intv.(dH_mc)),nx)

      if mc_opts.hhj_rnd_all == true
        H_mc,dH_mc = Rnd_Out_H_All(H_mc,dH_mc,mc_opts.hhj_rnd_all_eps)
      elseif mc_opts.hhj_rnd == true
        H_mc,dH_mc = Rnd_Out_H_Intv(H_mc,dH_mc,mc_opts.hhj_rnd_eps)
      end

      # applies parametric iteration
      if (mc_opts.style == "NewtonGS")
        MC_NewtonGS!(z_mc,x_mc,dH_mc,H_mc,nx)
      elseif (mc_opts.style == "KrawczykCW")
        MC_KrawczykCW!(z_mc,x_mc,dH_mc,H_mc,nx)
      else
        error("Unsupported Style of Implicit Relaxation")
      end
    end
    return x_mc
end
