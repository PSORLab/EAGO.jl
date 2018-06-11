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

struct poorSettings
	nrow

end

#=
Best script for Poor Man's LP

'''
	poorLP_contractor

Takes a number of reptitions, a interval vector `X`, a sparse matrix `A`, and

'''

for T in (:Interval,:MCInterval)
	@eval function poorLP_contractor(X::Vector{$T{Float64}},A::SparseMatrixCSC{Float64,Int},b::Vector{Float64},rept::Int)
		xL::Vector{Float64} = lo.(X)
		xU::Vector{Float64} = hi.(X)
		sub_term::Float64 = 0.0
		temp_sum::Float64 = 0.0
		for r = 1:rept
			for i = 1:s.n
				# Sum the rows times the vector
				temp_sum = b[i]
				for j=A.colptr[i]:(A.colptr[i+1])
					if (A.nzval[j] > 0.0 )
						temp_sum += A.nzval[j]*min(xL[A.rowval[j]], xU[A.rowval[j]])
					else
						temp_sum += A.nzval[j]*max(xL[A.rowval[j]], xU[A.rowval[j]])
					end
				end
				# For each element in the row apply the Poor Man’s LP
				for j=A.colptr[i]:(A.colptr[i+1])
					if (A.nzval[j] > 0.0 )
						sub_term = A.val[j]*min(xL[A.rowval[j]], xU[A.rowval[j]])
						temp_sum -= sub_term                                                     # Subtracts term formed from current nonzero values
						xU[A.rowval[j]] = min(xU[A.rowval[j]],temp_sum/A.nzval[j])               # Applies contractor
						temp_sum += sub_term                                                     # Recovers term
					else
						sub_term = A.nzval[j]*min(xL[A.rowval[j]], xU[A.rowval[j]])
						temp_sum -= sub_term
						xL[A.rowval[j]] = min(xL[A.rowval[j]],temp_sum/A.nzval[j])
						temp_sum += sub_term
					end
					(xU[A.rowval[j]] < xL[A.rowval[j]]) && (return X, false)                      # Returns contractor fathoms box
				end
			end
		end
		$T{Float64}.(xL,xU), true
	end
end

X = [Interval(-2.0,4.0),Interval(0.0,4.0),Interval(-1.0,1.0)]
S = sparse([-1.0 -1.0 0.0; 0.0 1.0 1.0])'
b = [-4.0; 1.0]

out = poorLP_contractor(X,S,b,2)
=#

function poorLP_contractor(s::poorSettings,X::Vector{Interval{Float64}})
	xL::Vector{Float64} = lo.(X)
	xU::Vector{Float64} = hi.(X)
	sub_term::Float64 = 0.0
	temp_sum::Float64 = 0.0
	for r = 1:s.rept
		for i = 1:s.nrow
			# Sum the rows times the vector
			temp_sum = b[i]
			for j=A.row_ptr[i]:(A.row_ptr[i+1])
				if (A.val[j] > 0.0 )
					temp_sum += A.val[j]*min(xL[A.col_ind], A.val[j]*xU[A.col_ind])
				else
					temp_sum += A.val[j]*max(xL[A.col_ind], A.val[j]*xU[A.col_ind])
				end
			end
			# For each element in the row apply the Poor Man’s LP
			for j=A.row_ptr[i]:(A.row_ptr[i+1])
				if (A.val[j] > 0.0 )
					sub_term = A.val[j]*min(xL[A.col_ind[j]], A.val[j]*xU[A.col_ind[j]])
					temp_sum -= sub_term                                                   # Subtracts term formed from current nonzero values
					xU[col_ind[j]] = min(xU[col_ind[j]],temp_sum/A.val[j])                 # Applies contractor
					temp_sum += sub_term                                                   # Recovers term
				else
					sub_term = A.val[j]*min(xL[A.col_ind[j]], A.val[j]*xU[A.col_ind[j]])
					temp_sum -= sub_term
					xL[col_ind[j]] = min(xL[col_ind[j]],temp_sum/A.val[j])
					temp_sum += sub_term
				end
			end
		end
	end
	Interval{Float64}.(xL,xU)
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
