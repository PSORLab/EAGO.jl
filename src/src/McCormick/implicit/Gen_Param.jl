#=

"""
    GenExpansionParams(h::Function, hj::Function, X::Vector{Interval{T}},
                       P::Vector{Interval{T}},pmid::Vector{T},mc_opts::mc_opts)

Generates relaxation of state variable at pmid in P and outputs results to use
as the parameters in subsequent calculations. Use a direct calculation of the
inverse function for preconditioning. Inputs are:
* `h::Function`: h(z,p) which implicit defines state function
* `hj::Function`: hj(z,p) w.r.t z which implicit defines state function
* `P::Vector{Interval}`: State variable bounds
* `X::Vector{Interval{T}}`: Decision variable bounds
* `pmid::Vector{T}`: Point at which to generate relaxation
* `mc_opts`: Options for generating implicit function relaxation
Outs the following the tuple `sto_out`:
* `sto_out::Vector{SMCg{np,T}}` - McCormick parameter object stack
--------------------------------------------------------------------------------
"""
function GenExpansionParams(h::Function, hj::Function,
                            X::Vector{V},
                            P::Vector{V},
                            pmid::Vector{T},mc_opts::mc_opts{T}) where {T<:AbstractFloat,V}

  nxi::Int64 = length(X)
  np::Int64 = length(P)
  szero::SVector{np,T} = @SVector zeros(np)
  sone::SVector{np,T} = @SVector ones(np)
  exp_opt::Vector{Any} = Any[nxi,np,mc_opts.lambda]
  sto_out::Vector{Vector{SMCg{np,V,T}}} = [[zero(SMCg{np,V,T}) for i=1:nxi] for j=1:(mc_opts.kmax+1)]

  x_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),false,[∅],[1.0]) for i=1:nxi]
  xa_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),false,[∅],[1.0]) for i=1:nxi]
  xA_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),false,[∅],[1.0]) for i=1:nxi]
  if mc_opts.z_rnd_all == true
    z_mc::Vector{SMCg{np,T}} = Rnd_Out_Z_All(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_all_eps)
  elseif mc_opts.z_rnd == true
    z_mc = Rnd_Out_Z_Intv(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_eps)
  else
    z_mc = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc
  end

  p_mc::Vector{SMCg{np,V,T}} = [SMCg(pmid[i],pmid[i],sone,sone,V(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  pref_mc::Vector{SMCg{np,V,T}} = copy(p_mc)
  aff_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,V(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nxi]

  H_mc::Vector{SMCg{np,V,T}} = h(z_mc,p_mc)
  dH_mc::Union{Vector{SMCg{np,V,T}},Array{SMCg{np,V,T},2}} = hj(aff_mc,p_mc)

  Y = (nxi == 1) ? [one(T)/mid(dH_mc[1].Intv)] : mid.(Intv.(dH_mc))
  Precondition!(H_mc,dH_mc,Y,nxi)
  # stores some things for repeated us
  optc = Any[szero,sone]

  # Begins loop to generate parameters
  sto_out[1] = copy(x_mc)
  for k=1:mc_opts.kmax

    aff_mc = [SMCg{np,Float64}(x_mc[i].cc,x_mc[i].cv,szero,szero,
                   @interval(x_mc[i].cv,x_mc[i].cc),false,[∅],[1.0]) for i=1:nxi]
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
    #println("ran to me 1")
    Y = (nxi == 1) ? [one(T)/mid(dH_mc[1].Intv)] : mid.(Intv.(dH_mc))
    #println("ran to me 2")
    # applies preconditioner
    Precondition!(H_mc,dH_mc,Y,nxi)
    #println("ran to me 3")

    for i=1:nxi
          println("x_mc[$i]: $(x_mc[i])")
          #println("xa_mc[$i]: $(xa_mc[i])")
          #println("xA_mc[$i]: $(xA_mc[i])")
          #println("z_mc[$i]: $(z_mc[i])")
    end

    # applies parametric iteration
    if (mc_opts.style == "NewtonGS")
      MC_NewtonGS!(z_mc,x_mc,dH_mc,H_mc,nxi)
    elseif (mc_opts.style == "KrawczykCW")
      MC_KrawczykCW!(z_mc,x_mc,dH_mc,H_mc,nxi)
    else
      error("Unsupported Style of Implicit Relaxation")
    end
    #println("ran to me 4")
    # update affine relaxations & correct

    for i=1:nxi
          println("x_mc[$i]: $(x_mc[i])")
          #println("xa_mc[$i]: $(xa_mc[i])")
          #println("xA_mc[$i]: $(xA_mc[i])")
          #println("z_mc[$i]: $(z_mc[i])")
    end
    for i=1:np
      println("p_mc[$i]: $(p_mc[i])")
    end

    Affine_Exp!(x_mc,p_mc,p_mc,xa_mc,xA_mc,z_mc,exp_opt)
    #println("ran to me 5")
    Correct_Exp!(z_mc,x_mc,X,nxi,np,mc_opts.aff_correct_eps)
    #println("ran to me 6")
    if mc_opts.z_rnd_all == true
      z_mc = Rnd_Out_Z_All(z_mc,mc_opts.z_rnd_all_eps)
    elseif mc_opts.z_rnd == true
      z_mc = Rnd_Out_Z_Intv(z_mc,mc_opts.z_rnd_eps)
    end

    # store relaxation
    sto_out[k+1] = copy(x_mc)
    end
  #return sto_out,sto_z,sto_x
  return sto_out
end

"""
    InGenExpansionParams(h::Function, hj::Function, X::Vector{Interval{T}},
                       P::Vector{Interval{T}},pmid::Vector{T},mc_opts::mc_opts)

Generates relaxation of state variable at pmid in P and outputs results to use
as the parameters in subsequent calculations. Use a sparse LDU factorization with
full pivoting to precondition and performs function evaluations in-place. Inputs are:
* `h::Function`: h!(H,z,p) which implicit defines state function
* `hj::Function`: hj!(H,z,p) w.r.t z which implicit defines state function
* `P::Vector{Interval}`: State variable bounds
* `X::Vector{Interval{T}}`: Decision variable bounds
* `pmid::Vector{T}`: Point at which to generate relaxation
* `mc_opts`: Options for generating implicit function relaxation
Outs the following the tuple `sto_out`:
* `sto_out::Vector{SMCg{np,T}}` - McCormick parameter object stack
--------------------------------------------------------------------------------
"""
function InGenExpansionParams(h!::Function, hj!::Function,
                            X::Vector{V},
                            P::Vector{V},
                            pmid::Vector{T},mc_opts::mc_opts{T}) where {T<:AbstractFloat,V}

  nx::Int64 = length(X)
  np::Int64 = length(P)
  szero::SVector{np,T} = @SVector zeros(np)
  sone::SVector{np,T} = @SVector ones(np)
  exp_opt::Vector{Any} = Any[nx,np,mc_opts.lambda]
  sto_out::Vector{Vector{SMCg{np,V,T}}} = [[zero(SMCg{np,V,T}) for i=1:nx] for j=1:(mc_opts.kmax+1)]

  x_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),false,[∅],[1.0]) for i=1:nx]
  xa_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),false,[∅],[1.0]) for i=1:nx]
  xA_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),false,[∅],[1.0]) for i=1:nx]

  if mc_opts.z_rnd_all == true
    z_mc::Vector{SMCg{np,V,T}} = Rnd_Out_Z_All(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_all_eps)
  elseif mc_opts.z_rnd == true
    z_mc = Rnd_Out_Z_Intv(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_eps)
  else
    z_mc = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc
  end

  p_mc::Vector{SMCg{np,V,T}} = [SMCg(pmid[i],pmid[i],sone,sone,V(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  pref_mc::Vector{SMCg{np,V,T}} = copy(p_mc)
  aff_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,V(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
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
  sto_out[1] = copy(x_mc)
  for k=1:mc_opts.kmax
    aff_mc = [SMCg{np,V,Float64}(x_mc[i].cc,x_mc[i].cv,szero,szero,
                   V(x_mc[i].cv,x_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
    if mc_opts.aff_rnd_all == true
      aff_mc = Rnd_Out_Z_All(aff_mc,mc_opts.aff_rnd_all_eps)
    elseif mc_opts.aff_rnd == true
      aff_mc = Rnd_Out_Z_Intv(aff_mc,mc_opts.aff_rnd_eps)
    end
    h!(H_mc,z_mc,p_mc)
    hj!(dH_mc,aff_mc,p_mc)
    Sparse_Precondition!(H_mc,dH_mc,mid.(Intv.(dH_mc)),SSto)
    dH_mc = transpose(dH_mc)

    # applies parametric iteration
    if (mc_opts.style == "NewtonGS")
      MCn_NewtonGS!(z_mc,x_mc,dH_mc,H_mc,nx)
    elseif (mc_opts.style == "KrawczykCW")
      MCn_KrawczykCW!(z_mc,x_mc,dH_mc,H_mc,nx)
    else
      error("Unsupported Style of Implicit Relaxation")
    end
    dH_mc = transpose(dH_mc)

    # update affine relaxations & correct
    Affine_Exp!(x_mc,p_mc,p_mc,xa_mc,xA_mc,z_mc,exp_opt)
    Correct_Exp!(z_mc,x_mc,X,nx,np,mc_opts.aff_correct_eps)
    if mc_opts.z_rnd_all == true
      z_mc = Rnd_Out_Z_All(z_mc,mc_opts.z_rnd_all_eps)
    elseif mc_opts.z_rnd == true
      z_mc = Rnd_Out_Z_Intv(z_mc,mc_opts.z_rnd_eps)
    end
    # store relaxation
    sto_out[k+1] = copy(x_mc)
    end
  return sto_out
end

"""
    InGenExpansionParams(h::Function, hj::Function, X::Vector{Interval{T}},
                       P::Vector{Interval{T}},pmid::Vector{T},mc_opts::mc_opts)

Generates relaxation of state variable at pmid in P and outputs results to use
as the parameters in subsequent calculations. Use a sparse LDU factorization with
full pivoting to precondition and performs function evaluations in-place. Inputs are:
* `h::Function`: h!(H,z,p) which implicit defines state function
* `hj::Function`: hj!(H,z,p) w.r.t z which implicit defines state function
* `P::Vector{Interval}`: State variable bounds
* `X::Vector{Interval{T}}`: Decision variable bounds
* `pmid::Vector{T}`: Point at which to generate relaxation
* `mc_opts`: Options for generating implicit function relaxation
Outs the following the tuple `sto_out`:
* `sto_out::Vector{SMCg{np,T}}` - McCormick parameter object stack
--------------------------------------------------------------------------------
"""
function IndGenExpansionParams(h!::Function, hj!::Function,
                            X::Vector{V},
                            P::Vector{V},
                            pmid::Vector{T},mc_opts::mc_opts{T}) where {V,T<:AbstractFloat}

  nx::Int64 = length(X)
  np::Int64 = length(P)
  szero::SVector{np,T} = @SVector zeros(np)
  sone::SVector{np,T} = @SVector ones(np)
  exp_opt::Vector{Any} = Any[nx,np,mc_opts.lambda]
  sto_out::Vector{Vector{SMCg{np,V,T}}} = [[zero(SMCg{np,V,T}) for i=1:nx] for j=1:(mc_opts.kmax+1)]

  x_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),false,[∅],[1.0]) for i=1:nx]
  xa_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),false,[∅],[1.0]) for i=1:nx]
  xA_mc::Vector{SMCg{np,V,T}} = [SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),false,[∅],[1.0]) for i=1:nx]

  if mc_opts.z_rnd_all == true
    z_mc::Vector{SMCg{np,T}} = Rnd_Out_Z_All(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_all_eps)
  elseif mc_opts.z_rnd == true
    z_mc = Rnd_Out_Z_Intv(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_eps)
  else
    z_mc = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc
  end

  p_mc::Vector{SMCg{np,T}} = [SMCg(pmid[i],pmid[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  pref_mc::Vector{SMCg{np,T}} = copy(p_mc)
  aff_mc::Vector{SMCg{np,T}} = [SMCg{np,T}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,@interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
  H_mc::Vector{SMCg{np,T}} = zeros(SMCg{np,T},nx)
  dH_mc::Array{SMCg{np,T},2} = zeros(SMCg{np,T},nx,nx)
  Y::Array{T,2} = zeros(T,nx,nx)

  # stores some things for repeated us
  optc = Any[szero,sone]

  # Begins loop to generate parameters
  sto_out[1] = copy(x_mc)
  for k=1:mc_opts.kmax
    aff_mc = [SMCg{np,V,Float64}(x_mc[i].cc,x_mc[i].cv,szero,szero,
                   V(x_mc[i].cv,x_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
    if mc_opts.aff_rnd_all == true
      aff_mc = Rnd_Out_Z_All(aff_mc,mc_opts.aff_rnd_all_eps)
    elseif mc_opts.aff_rnd == true
      aff_mc = Rnd_Out_Z_Intv(aff_mc,mc_opts.aff_rnd_eps)
    end

    h!(H_mc,z_mc,p_mc)
    hj!(dH_mc,aff_mc,p_mc)
    Dense_Precondition!(H_mc,dH_mc,mid.(Intv.(dH_mc)),nx)

    # applies parametric iteration
    if (mc_opts.style == "NewtonGS")
      MC_NewtonGS!(z_mc,x_mc,dH_mc,H_mc,nx)
    elseif (mc_opts.style == "KrawczykCW")
      MC_KrawczykCW!(z_mc,x_mc,dH_mc,H_mc,nx)
    else
      error("Unsupported Style of Implicit Relaxation")
    end

    # update affine relaxations & correct
    Affine_Exp!(x_mc,p_mc,p_mc,xa_mc,xA_mc,z_mc,exp_opt)
    Correct_Exp!(z_mc,x_mc,X,nx,np,mc_opts.aff_correct_eps)
    if mc_opts.z_rnd_all == true
      z_mc = Rnd_Out_Z_All(z_mc,mc_opts.z_rnd_all_eps)
    elseif mc_opts.z_rnd == true
      z_mc = Rnd_Out_Z_Intv(z_mc,mc_opts.z_rnd_eps)
    end
    # store relaxation
    sto_out[k+1] = copy(x_mc)
    end
  #return sto_out,sto_z,sto_x
  return sto_out
end
=#

function PreconditionSMCg(h::Function,
                          hj::Function,
                          z_mc::Vector{SMCg{N,V,T}},
                          aff_mc::Vector{SMCg{N,V,T}},
                          p_mc::Vector{SMCg{N,V,T}},
                          opt::mc_opts{T}) where {N,V,T<:AbstractFloat}
    if (opt.LAlg == :DenseBand)
        H,J = SMCg_DenseBand_Precondition!(h,hj,z_mc,aff_mc,p_mc,opt)
    elseif (opt.LAlg == :DenseBlockDiag)
        H,J = SMCg_DenseBlockDiag_Precondition!(h,hj,z_mc,aff_mc,p_mc,opt)
    elseif (opt.LAlg == :Dense)
        H,J = SMCg_Dense_Precondition!(h,hj,z_mc,aff_mc,p_mc,opt)
    else
        error("Unsupported Linear Algebra Style")
    end
    return H,J
end

function SMCg_Dense_Precondition!(h::Function,
                             hj::Function,
                             z_mc::Vector{SMCg{N,V,T}},
                             aff_mc::Vector{SMCg{N,V,T}},
                             p_mc::Vector{SMCg{N,V,T}},
                             opt::mc_opts{T}) where {N,V,T<:AbstractFloat}
    #println("start dense precond")
    H::Vector{SMCg{N,V,T}} = h(z_mc,p_mc)
    #println("finish H eval")
    J::Union{Array{SMCg{N,V,T},2},Vector{SMCg{N,V,T}}} = hj(aff_mc,p_mc)
    #println("finish J eval")
    #println("typeof J: $(J)")
    #println("typeof J: $(typeof(J))")
    #println("J.Intv: $(J[1].Intv)")
    #println("J.Intv mid: $(mid(J[1].Intv))")
  #  println("mid.(J): $(mid.(J))")
    Y = [mid(J[i,j].Intv) for i=1:opt.nx, j=1:opt.nx]
    #println("finish Y eval")
    if (opt.nx == 1)
        #println("start 1D")
        YH::Vector{SMCg{N,V,T}} = H/Y[1]
        #println("finish YH eval")
        YJ::Union{Array{SMCg{N,V,T},2},Vector{SMCg{N,V,T}}} = J/Y[1]
        #println("finish YJ eval")
    else
        #println("start multD")
        F = lufact(Y)
        #println("finish F eval")
        YH = F\H
        #println("finish YH eval")
        YJ = F\J
        #println("finish YJ eval")
    end
    return YH,YJ
end

function SMCg_Dense_Newton_GS!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                        YdH_mc::Vector{SMCg{N,V,T}},
                        YH_mc::Vector{SMCg{N,V,T}},mc_opts::mc_opts{T}) where {N,V,T<:AbstractFloat}
    S1::SMCg{N,V,T} = zero(x_mc[1])
    x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
    for i=1:mc_opts.nx
      S1 = zero(x_mc[1])
      for j=1:mc_opts.nx
        if (i<j)
          S1 = S1 + YdH_mc[i,j]*(x_mc[j]-z_mc[j])
        elseif (j<i)
          S1 = S1 + YdH_mc[i,j]*(x_mc[j]-z_mc[j])
        end
      end
      println("YH_mc[i]: $(YH_mc[i])")
      println("YdH_mc[i,i]: $(YdH_mc[i,i])")
      println("dividy: $(YH_mc[i]/YdH_mc[i,i])")
      x_mc[i] = z_mc[i] - (YH_mc[i]+S1)/YdH_mc[i,i]
      println("x_mc[i]: $(x_mc[i])")
      x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
      println("x_mc[i]: $(x_mc[i])")
    end
end

function SMCg_Dense_Newton_GS!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                        YdH_mc::Array{SMCg{N,V,T},2},
                        YH_mc::Vector{SMCg{N,V,T}},mc_opts::mc_opts{T}) where {N,V,T<:AbstractFloat}
    S1::SMCg{N,V,T} = zero(x_mc[1])
    x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
    for i=1:mc_opts.nx
      S1 = zero(x_mc[1])
      for j=1:mc_opts.nx
        if (i<j)
          S1 = S1 + YdH_mc[i,j]*(x_mc[j]-z_mc[j])
        elseif (j<i)
          S1 = S1 + YdH_mc[i,j]*(x_mc[j]-z_mc[j])
        end
      end
      println("YH_mc[i]: $(YH_mc[i])")
      println("YdH_mc[i,i]: $(YdH_mc[i,i])")
      println("dividy: $(YH_mc[i]/YdH_mc[i,i])")
      x_mc[i] = z_mc[i] - (YH_mc[i]+S1)/YdH_mc[i,i]
      x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
    end
end

function SMCg_Dense_Krawczyk_CW!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                                 YdH_mc::Vector{SMCg{N,V,T}},
                                 YH_mc::Vector{SMCg{N,V,T}},mc_opts::mc_opts{T}) where{N,V,T<:AbstractFloat}
  println("SMCg1 dense Kraw Start")
  S1::SMCg{N,V,T} = zero(x_mc[1])
  println("SMCg1 dense Kraw S1")
  x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
  println("SMCg1 dense Kraw int")
  for i=1:mc_opts.nx
    S1 = zero(x_mc[1])
    println("SMCg1 dense Kraw S1 $i")
    for j=1:mc_opts.nx
      if (i<j)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      elseif (j<i)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      else
        S1 = S1 + (one(T)-YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      end
    end
    println("SMCg1 dense Kraw S1 out $i")
    x_mc[i] =  z_mc[i] - YH_mc[i] + S1
    println("start final cut")
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
    println("end final cut")
  end
end

function SMCg_Dense_Krawczyk_CW!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                                 YdH_mc::Array{SMCg{N,V,T},2},
                                 YH_mc::Vector{SMCg{N,V,T}},mc_opts::mc_opts{T}) where{N,V,T<:AbstractFloat}
  println("SMCg1 dense Kraw Start")
  S1::SMCg{N,V,T} = zero(x_mc[1])
  println("SMCg1 dense Kraw S1")
  x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
  println("SMCg1 dense Kraw int")
  for i=1:mc_opts.nx
    S1 = zero(x_mc[1])
    println("SMCg1 dense Kraw S1 $i")
    for j=1:mc_opts.nx
      if (i<j)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      elseif (j<i)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      else
        S1 = S1 + (one(T)-YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      end
    end
    println("SMCg1 dense Kraw S1 out $i")
    x_mc[i] =  z_mc[i] - YH_mc[i] + S1
    println("start final cut")
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
    println("end final cut")
  end
end

function PSMCg_Kernel!(h,hj,z_mc,aff_mc,p_mc,x_mc,opt::mc_opts{T}) where {T}
  #println("start Kernel")
  H,J = PreconditionSMCg(h,hj,z_mc,aff_mc,p_mc,opt)
  #println("ran preconditioner")
  if (opt.CTyp == :Newton)
    if (opt.LAlg == :Dense)
      SMCg_Dense_Newton_GS!(z_mc,x_mc,J,H,opt)
    elseif (opt.LAlg == :DenseBand)
      SMCg_DenseBand_Newton_GS!(z_mc,x_mc,J,H,opt)
    else
      error("The linear algebra type $(LAlg) is not currently supported. The
             linear algebra styles currently supported are :Dense and :DenseBanded.")
    end
  elseif (opt.CTyp == :Krawczyk)
    if (opt.LAlg == :Dense)
      #println("start dense kraw")
      SMCg_Dense_Krawczyk_CW!(z_mc,x_mc,J,H,opt)
      #println("finish dense kraw")
    elseif (opt.LAlg == :DenseBand)
      SMCg_DenseBand_Krawczyk_CW!(z_mc,x_mc,J,H,opt)
    else
      error("The linear algebra type $(LAlg) is not currently supported. The
             linear algebra styles currently supported are :Dense and
             :DenseBanded.")
    end
  else
      error("The contractor type $(CTyp) is not currently supported. The
             contractors :Newton and :Krawczyk are currently supported.")
  end
end

function GenExpansionParams(h::Function, hj::Function,
                      X::Vector{V},
                      P::Vector{V},
                      pmid::Vector{T},mc_opts::mc_opts{T}) where {T<:AbstractFloat,V}

  println("start function")
  nxi::Int64 = length(X)
  np::Int64 = length(P)
  szero::SVector{np,T} = zeros(SVector{np,T})
  println("szero")
  sone::SVector{np,T} = ones(SVector{np,T})
  println("one")
  SP = SVector{np,V}(P)
  println("SP")
  println("pmid: $pmid")
  println("size(pmid): $(size(pmid))")
  SP0 = SVector{np,T}(pmid)
  println("ran declaration block 1")
  exp_opt::Vector{Any} = Any[nxi,np,mc_opts.lambda]
  #sto_out::Vector{Vector{SMCg{np,V,T}}} = [SMCg{np,V,T}[zero(SMCg{np,V,T}) for i=1:nxi] for j=1:(mc_opts.kmax+1)]

  println("ran declaration block 2")
  #x_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),false,SP,SP0) for i=1:nxi]
  x_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(X[i].hi,X[i].lo,sone,sone,V(X[i].lo,X[i].hi),false,SP,SP0) for i=1:nxi]
  xa_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),false,SP,SP0) for i=1:nxi]
  xA_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),false,SP,SP0) for i=1:nxi]
  z_mc::Vector{SMCg{np,V,T}} = mc_opts.lambda*xa_mc+(one(T)-mc_opts.lambda)*xA_mc

  println("ran declaration block 3")
  p_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(pmid[i],pmid[i],sone,sone,V(P[i].lo,P[i].hi),false,SP,SP0) for i=1:np]
  pref_mc::Vector{SMCg{np,V,T}} = copy(p_mc)
  aff_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,V(xa_mc[i].cv,xA_mc[i].cc),false,SP,SP0) for i=1:nxi]
  sto_out::Vector{Vector{SMCg{np,V,T}}} = Vector{SMCg{np,V,T}}[x_mc for j=1:(mc_opts.kmax+1)]
  sto_out[1] = copy(x_mc)
  optc = Any[szero,sone]

  println("xa_mc: $(xa_mc)")
  println("xA_mc: $(xA_mc)")
  println("z_mc: $(z_mc)")
  println("aff_mc: $(aff_mc)")
  println("ran declaration")
  for k=1:mc_opts.kmax
    println("iteration number: $k")
    #aff_mc = SMCg{np,V,Float64}[SMCg{np,V,Float64}(x_mc[i].cc,x_mc[i].cv,szero,szero,
    #               V(x_mc[i].cv,x_mc[i].cc),false,SP,SP0) for i=1:nxi]
    println("aff_mc: $(aff_mc)")
    println("z_mc: $(z_mc)")
    println("x_mc: $(x_mc)")
    println("ran affine")
    PSMCg_Kernel!(h,hj,z_mc,aff_mc,p_mc,x_mc,mc_opts)
    println("ran kernel")
    println("z_mc: $(z_mc)")
    println("x_mc: $(x_mc)")
    Affine_Exp!(x_mc,p_mc,pmid,xa_mc,xA_mc,z_mc,mc_opts)
    println("ran expansion")
    println("z_mc: $(z_mc)")
    println("x_mc: $(x_mc)")
    Correct_Exp!(z_mc,x_mc,X,nxi,np,mc_opts.aff_correct_eps)
    aff_mc = SMCg{np,V,Float64}[SMCg{np,V,Float64}(xA_mc[i].cc,xa_mc[i].cv,xA_mc[i].cc_grad,xa_mc[i].cv_grad,
                   V(x_mc[i].cv,x_mc[i].cc),false,SP,SP0) for i=1:nxi]
    println("ran correct")
    println("x_mc: $(x_mc)")
    # store relaxation
    sto_out[k+1] = copy(x_mc)
    println("ran sto")
    end
  println("finished loop")
  #return sto_out,sto_z,sto_x
  return sto_out
end
