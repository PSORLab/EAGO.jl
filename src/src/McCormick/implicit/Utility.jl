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
    #println("opt.nx: $(opt.nx)")
    H::Vector{SMCg{N,V,T}} = h(z_mc,p_mc)
    #println("size H: $(size(H))")
    J::VecOrMat{SMCg{N,V,T}} = hj(aff_mc,p_mc)
    #println("size J: $(size(J))")
    Y = [mid(J[i,j].Intv) for i=1:opt.nx, j=1:opt.nx]
    #println("size Y: $(size(Y))")
    if (opt.nx == 1)
        YH::Vector{SMCg{N,V,T}} = H/Y[1]
        YJ::VecOrMat{SMCg{N,V,T}} = J/Y[1]
    else
        F = lufact(Y)
        YH = F\H
        YJ = F\J
    end
    return YH,YJ
end


"""
    Smooth_Cut(x_mc::SMCg{N,T},x_mc_int::SMCg{N,T})

An operator that cuts the `x_mc` object using the `x_mc_int bounds` in a
differentiable fashion.
"""
function Smooth_Cut(x_mc::SMCg{N,V,T},x_mc_int::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  t_cv::SMCg{N,V,T} = x_mc + max(zero(T),x_mc_int-x_mc)
  t_cc::SMCg{N,V,T} = x_mc + min(zero(T),x_mc-x_mc_int)
  return SMCg{N,V,T}(t_cc.cc,t_cv.cv,t_cc.cc_grad,t_cv.cv_grad,
                   (x_mc.Intv ∩ x_mc_int.Intv),(t_cv.cnst && t_cc.cnst))
end

"""
    Final_Cut(x_mc::SMCg{N,T},x_mc_int::SMCg{N,T})

An operator that cuts the `x_mc` object using the `x_mc_int bounds` in a
differentiable or nonsmooth fashion as specified by the `MC_param.mu flag`.
"""
function Final_Cut(x_mc::SMCg{N,V,T},x_mc_int::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  if (MC_param.mu < 1)
    Intv::V = x_mc.Intv ∩ x_mc_int.Intv
    if (x_mc.cc <= x_mc_int.cc)
      cc::T = x_mc.cc
      cc_grad::SVector{N,T} = x_mc.cc_grad
    else
      cc = x_mc_int.cc
      cc_grad = x_mc_int.cc_grad
    end
    if (x_mc.cv >= x_mc_int.cv)
      cv::T = x_mc.cv
      cv_grad::SVector{N,T} = x_mc.cv_grad
    else
      cv = x_mc_int.cv
      cv_grad = x_mc_int.cv_grad
    end
    x_mc::SMCg{N,V,T} = SMCg{N,V,T}(cc,cv,cc_grad,cv_grad,(x_mc.Intv ∩ x_mc_int.Intv),x_mc.cnst)
  else
    x_mc = Smooth_Cut(x_mc,x_mc_int)
  end
  return x_mc
end

"""
    Rnd_Out_Z_Intv(z_mct::SMCg{N,T},epsvi::Float64)

Rounds the interval of the `z_mct` vector elements out by `epsvi`.
"""
function Rnd_Out_Z_Intv(z_mct::Vector{SMCg{N,V,T}},epsvi::S) where {N,V,S<:AbstractFloat,T<:AbstractFloat}
  epsv::T = convert(T,epsvi)
  return [SMCg{N,V,T}(z_mct[i].cc,z_mct[i].cv,
             z_mct[i].cc_grad, z_mct[i].cv_grad,
             V(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),z_mct[i].cnst) for i=1:length(z_mct)]
end

"""
    Rnd_Out_Z_All(z_mct::Vector{SMCg{N,T}},epsvi::S)

Rounds the interval and relaxation bounds of the `z_mct` vector elements out by `epsvi`.
"""
function Rnd_Out_Z_All(z_mct::Vector{SMCg{N,V,T}},epsvi::S) where {N,V,S<:AbstractFloat,T<:AbstractFloat}
  epsv::T = convert(T,epsvi)
  return [SMCg{N,V,T}(z_mct[i].cc+epsv,z_mct[i].cv-epsv,
             z_mct[i].cc_grad, z_mct[i].cv_grad,
             V(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),
             z_mct[i].cnst) for i=1:length(z_mct)]
end

"""
    Rnd_Out_H_Intv(z_mct::Vector{SMCg{N,T}},Y_mct::Array{SMCg{N,T},2},epsvi::S)

Rounds the interval bounds of the `z_mct` and `Y_mct` elements out by `epsvi`.
"""
function Rnd_Out_H_All(z_mct::Vector{SMCg{N,V,T}},Y_mct::Array{SMCg{N,V,T},2},epsvi::S) where {N,V,S<:AbstractFloat,T<:AbstractFloat}
  epsv::T = convert(T,epsvi)
  temp1::Vector{SMCg{N,V,T}} = [SMCg{N,V,T}(z_mct[i].cc+epsv,z_mct[i].cv-epsv,
                                        z_mct[i].cc_grad, z_mct[i].cv_grad,
                                        V(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),
                                        z_mct[i].cnst) for i=1:length(z_mct)]
  temp2::Array{SMCg{N,V,T},2} = [SMCg{N,V,T}(Y_mct[i,j].cc+epsv,Y_mct[i,j].cv-epsv,
                                        Y_mct[i,j].cc_grad, Y_mct[i,j].cv_grad,
                                        V(Y_mct[i,j].Intv.lo-epsv, Y_mct[i,j].Intv.hi+epsv),
                                        Y_mct[i,j].cnst) for i=1:length(z_mct), j=1:length(z_mct)]
  return temp1,temp2
end

"""
    Rnd_Out_H_All(z_mct::Vector{SMCg{N,T}},Y_mct::Array{SMCg{N,T},2},epsvi::S)

Rounds the interval and relaxation bounds of the `z_mct` and `Y_mct` elements out by `epsvi`.
"""
function Rnd_Out_H_Intv(z_mct::Vector{SMCg{N,V,T}},Y_mct::Array{SMCg{N,V,T},2},epsvi::S) where {N,V,S<:AbstractFloat,T<:AbstractFloat}
  epsv::T = convert(T,epsvi)
  temp1::Vector{SMCg{N,V,T}} = [SMCg{N,V,T}(z_mct[i].cc,z_mct[i].cv,
             z_mct[i].cc_grad, z_mct[i].cv_grad,
             V(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),
             z_mct[i].cnst) for i=1:length(z_mct)]
  temp2::Array{SMCg{N,V,T},2} = [SMCg{N,V,T}(Y_mct[i,j].cc,Y_mct[i,j].cv,
             Y_mct[i,j].cc_grad, Y_mct[i,j].cv_grad,
             V(Y_mct[i,j].Intv.lo-epsv, Y_mct[i,j].Intv.hi+epsv),
             Y_mct[i,j].cnst) for i=1:length(z_mct), j=1:length(z_mct)]
  return temp1,temp2
end



#=
"""
    Precondition(hm::Vector{SMCg{N,T}},hJm::Union{Vector{SMCg{N,T}},Array{SMCg{N,T},2}},
                 Y::Union{Vector{T},Array{T,2}},nx::Int64)

Preconditions `hm` and `hJm` by `Y` in place where all dimensions are `nx`.
"""
function Precondition!(hm::Vector{SMCg{N,V,T}},hJm::Vector{SMCg{N,V,T}},
                      Y::Vector{T},nx::Int64) where {N,V,T<:AbstractFloat}
  S1::SMCg{N,V,T},S2::SMCg{N,V,T} = zero(SMCg{N,V,T}),zero(SMCg{N,V,T})
  for i=1:nx
    S2 = zero(SMCg{N,V,T})
    for j=1:nx
      S1 = zero(SMCg{N,V,T})
      for k=1:nx
        S1 = S1 + Y[i,k]*hJm[k,j]
      end
      hJm[i,j] = S1
      S2 += Y[i,j]*hm[j]
    end
    hm[i] = S2
  end
  return hm,hJm
end

function Precondition!(hm::Vector{SMCg{N,V,T}},hJm::Vector{SMCg{N,V,T}},
                      Y::Array{T,2},nx::Int64) where {N,V,T<:AbstractFloat}
  S1::SMCg{N,V,T},S2::SMCg{N,V,T} = zero(SMCg{N,V,T}),zero(SMCg{N,V,T})
  for i=1:nx
    S2 = zero(SMCg{N,V,T})
    for j=1:nx
      S1 = zero(SMCg{N,V,T})
      for k=1:nx
        S1 = S1 + Y[i,k]*hJm[k,j]
      end
      hJm[i,j] = S1
      S2 += Y[i,j]*hm[j]
    end
    hm[i] = S2
  end
  return hm,hJm
end

function Precondition!(hm::Vector{SMCg{N,V,T}},hJm::Array{SMCg{N,V,T},2},
                      Y::Array{T,2},nx::Int64) where {N,V,T<:AbstractFloat}
  S1::SMCg{N,V,T},S2::SMCg{N,V,T} = zero(SMCg{N,V,T}),zero(SMCg{N,V,T})
  for i=1:nx
    S2 = zero(SMCg{N,V,T})
    for j=1:nx
      S1 = zero(SMCg{N,V,T})
      for k=1:nx
        S1 = S1 + Y[i,k]*hJm[k,j]
      end
      hJm[i,j] = S1
      S2 += Y[i,j]*hm[j]
    end
    hm[i] = S2
  end
  return hm,hJm
end

function Precondition!(hm::Vector{SMCg{N,V,T}},hJm::Array{SMCg{N,V,T},2},
                      Y::Vector{T},nx::Int64) where {N,V,T<:AbstractFloat}
  S1::SMCg{N,V,T},S2::SMCg{N,V,T} = zero(SMCg{N,V,T}),zero(SMCg{N,V,T})
  for i=1:nx
    S2 = zero(SMCg{N,V,T})
    for j=1:nx
      S1 = zero(SMCg{N,V,T})
      for k=1:nx
        S1 = S1 + Y[i,k]*hJm[k,j]
      end
      hJm[i,j] = S1
      S2 += Y[i,j]*hm[j]
    end
    hm[i] = S2
  end
  return hm,hJm
end
=#
