"""
    MC_NewtonGS(z_mc::Vector{SMCg{N,T}},x_mc::Vector{SMCg{N,T}},
                         YdH_mc::Array{SMCg{N,T},2},YH_mc::Array{SMCg{N,T},2},
                         nx::Int64,np::Int64,optc::Vector{Any})

Performs a Newton Gauss-Seidel bounding iteration using McCormick relaxations
where `x_mc` is the initial state variable, `z_mc` are relaxations of the affine
function in `x_mc`, `YdH_mc` is the preconditioned jacobian, `YH_mc` is the
equality constraints jacobian, `nx` size of the state variables, `np` is the
size of the decision space, and `optc` is a vector with stored values.
"""
function MC_NewtonGS!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                     YdH_mc::Array{SMCg{N,V,T},2},
                     YH_mc::Vector{SMCg{N,V,T}},nx::Int64) where {N,V,T<:AbstractFloat}

  S1::SMCg{N,V,T} = zero(SMCg{N,V,T})
  x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
  for i=1:nx
    S1 = zero(SMCg{N,V,T})
    for j=1:nx
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

function MC_NewtonGS!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                     YdH_mc::Vector{SMCg{N,V,T}},
                     YH_mc::Vector{SMCg{N,V,T}},nx::Int64) where {N,V,T<:AbstractFloat}

  S1::SMCg{N,V,T} = zero(SMCg{N,V,T})
  x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
  for i=1:nx
    S1 = zero(SMCg{N,V,T})
    for j=1:nx
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

"""
    MC_KrawczykCW(z_mc::Vector{SMCg{N,T}},x_mc::Vector{SMCg{N,T}},
                           YdH_mc::Array{SMCg{N,T},2},YH_mc::Array{SMCg{N,T},2},
                           nx::Int64,np::Int64,optc::Vector{Any})

Performs a componentwise Krawczyk bounding iteration using McCormick relaxations
where `x_mc` is the initial state variable, `z_mc` are relaxations of the affine
function in `x_mc`, `YdH_mc` is the preconditioned jacobian, `YH_mc` is the
equality constraints jacobian, `nx` size of the state variables, `np` is the
size of the decision space, and `optc` is a vector with stored values.
"""
function MC_KrawczykCW!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                        YdH_mc::Array{SMCg{N,V,T},2},
                        YH_mc::Vector{SMCg{N,V,T}},nx::Int64) where{N,V,T<:AbstractFloat}
  #println("z_mc: $z_mc")
  #println("x_mc: $x_mc")
  #println("YdH_mc: $YdH_mc")
  #println("YH_mc: $YH_mc")
  #println("nx: $nx")
  S1::SMCg{N,V,T} = zero(SMCg{N,V,T})
  x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
  #println("YH_mc: $S1")
  #println("x_mc_int: $x_mc_int")
  for i=1:nx
    S1 = zero(SMCg{N,V,T})
    for j=1:nx
      #println("j: $j")
      if (i<j)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      elseif (j<i)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      else
        S1 = S1 + (one(T)-YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      end
    end
    x_mc[i] =  z_mc[i] - YH_mc[i] + S1
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
  end
end

function MC_KrawczykCW!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                        YdH_mc::Vector{SMCg{N,V,T}},
                        YH_mc::Vector{SMCg{N,V,T}},nx::Int64) where{N,V,T<:AbstractFloat}
  #println("z_mc: $z_mc")
  #println("x_mc: $x_mc")
  #println("YdH_mc: $YdH_mc")
  #println("YH_mc: $YH_mc")
  #println("nx: $nx")
  S1::SMCg{N,V,T} = zero(SMCg{N,V,T})
  x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
  #println("YH_mc: $S1")
  #println("x_mc_int: $x_mc_int")
  for i=1:nx
    S1 = zero(SMCg{N,V,T})
    for j=1:nx
      #println("j: $j")
      if (i<j)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      elseif (j<i)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      else
        S1 = S1 + (one(T)-YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      end
    end
    x_mc[i] =  z_mc[i] - YH_mc[i] + S1
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
  end
end

"""
    MCn_NewtonGS(z_mc::Vector{SMCg{N,T}},x_mc::Vector{SMCg{N,T}},
                 YdH_mc::Union{Vector{SMCg{N,T}},SparseMatrixCSC{SMCg{N,T},Int64}},
                 YH_mc::Array{SMCg{N,T},2},nx::Int64,np::Int64,optc::Vector{Any})

Performs a Newton Gauss-Seidel bounding iteration using McCormick relaxations
where `x_mc` is the initial state variable, `z_mc` are relaxations of the affine
function in `x_mc`, `YdH_mc` is the preconditioned jacobian, `YH_mc` is the
equality constraints jacobian, `nx` size of the state variables, `np` is the
size of the decision space, and `optc` is a vector with stored values.
"""
function MCn_NewtonGS!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                     YdH_mc::Vector{SMCg{N,V,T}},
                     YH_mc::Vector{SMCg{N,V,T}},nx::Int64) where {N,V,T<:AbstractFloat}

  S1::SMCg{N,V,T} = zero(SMCg{N,V,T})
  S2::SMCg{N,V,T} = zero(SMCg{N,V,T})
  S3::SMCg{N,V,T} = zero(SMCg{N,V,T})
  x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
  for i=1:nx
    S1 = zero(SMCg{N,V,T})
    S2 = zero(SMCg{N,V,T})
    S3 = zero(SMCg{N,V,T})
    for q=(YdH_mc.colptr[i]):(YdH_mc.colptr[i+1]-1)
      if (i < YdH_mc.rowval[q])
        S1 += YdH_mc.nzval[q]*(x_mc[YdH_mc.rowval[q]]-z_mc[YdH_mc.rowval[q]])
      elseif (i > YdH_mc.rowval[q])
        S2 += YdH_mc.nzval[q]*(x_mc[YdH_mc.rowval[q]]-z_mc[YdH_mc.rowval[q]])
      else
        S3 = YdH_mc.nzval[q]
      end
    end
    x_mc[i] = z_mc[i] - (YH_mc[i]+S1+S2)/S3
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
  end
end

function MCn_NewtonGS!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                     YdH_mc::SparseMatrixCSC{SMCg{N,V,T},Int64},
                     YH_mc::Vector{SMCg{N,V,T}},nx::Int64) where {N,V,T<:AbstractFloat}

  S1::SMCg{N,V,T} = zero(SMCg{N,V,T})
  S2::SMCg{N,V,T} = zero(SMCg{N,V,T})
  S3::SMCg{N,V,T} = zero(SMCg{N,V,T})
  x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
  for i=1:nx
    S1 = zero(SMCg{N,V,T})
    S2 = zero(SMCg{N,V,T})
    S3 = zero(SMCg{N,V,T})
    for q=(YdH_mc.colptr[i]):(YdH_mc.colptr[i+1]-1)
      if (i < YdH_mc.rowval[q])
        S1 += YdH_mc.nzval[q]*(x_mc[YdH_mc.rowval[q]]-z_mc[YdH_mc.rowval[q]])
      elseif (i > YdH_mc.rowval[q])
        S2 += YdH_mc.nzval[q]*(x_mc[YdH_mc.rowval[q]]-z_mc[YdH_mc.rowval[q]])
      else
        S3 = YdH_mc.nzval[q]
      end
    end
    x_mc[i] = z_mc[i] - (YH_mc[i]+S1+S2)/S3
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
  end
end

"""
    MCn_KrawczykCW(z_mc::Vector{SMCg{N,T}},x_mc::Vector{SMCg{N,T}},
                           YdH_mc::Array{SMCg{N,T},2},YH_mc::Array{SMCg{N,T},2},
                           nx::Int64,np::Int64,optc::Vector{Any})

Performs a componentwise Krawczyk bounding iteration using McCormick relaxations
where `x_mc` is the initial state variable, `z_mc` are relaxations of the affine
function in `x_mc`, `YdH_mc` is the preconditioned jacobian, `YH_mc` is the
equality constraints jacobian, `nx` size of the state variables, `np` is the
size of the decision space, and `optc` is a vector with stored values.
"""
function MCn_KrawczykCW!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                        YdH_mc::Vector{SMCg{N,V,T}},
                        YH_mc::Vector{SMCg{N,V,T}},nx::Int64) where{N,V,T<:AbstractFloat}
  S1::SMCg{N,V,T} = zero(SMCg{N,V,T})
  x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
  for i=1:nx
    S1 = zero(SMCg{N,T})
    for q=(YdH_mc.colptr[i]):(YdH_mc.colptr[i+1]-1)
      if (i == YdH_mc.rowval[q])
        S1 += (one(T)-YdH_mc.nzval[q])*(x_mc[YdH_mc.rowval[q]]-z_mc[YdH_mc.rowval[q]])
      else
        S1 += -YdH_mc.nzval[q]*(x_mc[YdH_mc.rowval[q]]-z_mc[YdH_mc.rowval[q]])
      end
    end
    x_mc[i] =  z_mc[i] - YH_mc[i] + S1
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
  end
end

function MCn_KrawczykCW!(z_mc::Vector{SMCg{N,V,T}},x_mc::Vector{SMCg{N,V,T}},
                        YdH_mc::SparseMatrixCSC{SMCg{N,V,T},Int64},
                        YH_mc::Vector{SMCg{N,V,T}},nx::Int64) where{N,V,T<:AbstractFloat}
  S1::SMCg{N,V,T} = zero(SMCg{N,V,T})
  x_mc_int::Vector{SMCg{N,V,T}} = copy(x_mc)
  for i=1:nx
    S1 = zero(SMCg{N,T})
    for q=(YdH_mc.colptr[i]):(YdH_mc.colptr[i+1]-1)
      if (i == YdH_mc.rowval[q])
        S1 += (one(T)-YdH_mc.nzval[q])*(x_mc[YdH_mc.rowval[q]]-z_mc[YdH_mc.rowval[q]])
      else
        S1 += -YdH_mc.nzval[q]*(x_mc[YdH_mc.rowval[q]]-z_mc[YdH_mc.rowval[q]])
      end
    end
    x_mc[i] =  z_mc[i] - YH_mc[i] + S1
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
  end
end
