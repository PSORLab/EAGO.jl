struct poorLPOpts
    repts::Int
    mrow::Int
end


# Linear Program Contractor
function poorLP_contractor(A::SparseMatrixCSC,b::Vector{Float64},
                           X::Array{Interval{Float64}},opt::poorLPOpts)

    # sets variable bounds & gets typical bounds
    nx::Int = length(X)
    na::Int = length(A)
    l::Vector{Float64} = Float64[X[i].lo for i=1:nx]
    u::Vector{Float64} = Float64[X[i].hi for i=1:nx]

    # storage array for row
    for r = 1:opts.repts
        for i=1:nx
            for j=(A.rowptr[i]):(A.rowptr[i+1]-1)
                # Loads the row
                # Apply inequality for each nonzero entry
            end
        end
    end

    Xout::Array{Interval{Float64}} = [Interval(l_out[i],u_out[i]) for i=1:nx]
    return temp_arr,feas
end

# Banded Matrix GaussSiedel
function BandedImplicitGS(p::Vector{SMCg{N,V,T}, pref::Vector{SMCg{N,V,T},
                          X::Vector{V}, A::Function, b::Function, opt::mc_opts{T}) where {N,T<:AbstractFloat,V<:AbstractInterval}
    l::Vector{T} = T[X[i].lo for i=1:nx]
    u::Vector{T} = T[X[i].hi for i=1:nx]

    x_mc::Vector{SMCg{N,V,T}} = SMCg{N,V,T}[SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),false) for i=1:N]
    xa_mc::Vector{SMCg{N,V,T}} = SMCg{N,V,T}[SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),false) for i=1:N]
    xA_mc::Vector{SMCg{N,V,T}} = SMCg{N,V,T}[SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),false) for i=1:N]
    z_mc::Vector{SMCg{N,V,T}} = mc_opts.lambda*xa_mc+(one(T)-mc_opts.lambda)*xA_mc
    aff_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,V(xa_mc[i].Intv.lo,xA_mc[i].Intv.hi),false) for i=1:N]

    H,J = PLP_Precondition(b,A,z_mc,aff_mc,p_mc,opt)
    for k=1:mc_opts.kmax
      PLP_Kernel!(h,hj,z_mc,aff_mc,p_mc,x_mc,mc_opts)
      Affine_Exp!(x_mc,p_mc,p_mc,xa_mc,xA_mc,z_mc,mc_opts)
      Correct_Exp!(z_mc,x_mc,X,nxi,np,mc_opts.aff_correct_eps)
      aff_mc = SMCg{np,V,Float64}[SMCg{np,V,Float64}(xA_mc[i].cc,xa_mc[i].cv,xA_mc[i].cc_grad,xa_mc[i].cv_grad,
                     V(xA_mc[i].Intv.lo,xA_mc[i].Intv.hi),false) for i=1:nxi]
      # store relaxation
      sto_out[k+1] = copy(x_mc)
      end
    return sto_out
end

function SMCg_Dense_Newton_GS!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                        YdH_mc::VecOrMat{SMCg{N,V,T}},
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
      x_mc[i] = z_mc[i] - (YH_mc[i]+S1)/YdH_mc[i,i]
      x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
    end
end

function PLP_Kernel!(J,H,z_mc,aff_mc,p_mc,x_mc,opt::mc_opts{T}) where {T}
     if (opt.LAlg == :Dense)
         SMCg_Dense_Newton_GS!(z_mc,x_mc,J,H,opt)
     elseif (opt.LAlg == :DenseBand)
         SMCg_DenseBand_Newton_GS!(z_mc,x_mc,J,H,opt)
     else
         error("The contractor type $(CTyp) is not currently supported. The
                contractors :Newton and :Krawczyk are currently supported.")
     end
end
