function Dense_Precondition!(H::Vector{S},J::Array{S,2},Y::Array{U,2},nx::Int64) where {S,U<:AbstractFloat}
         lu = lufact(Y)
         HJ::Array{S,2} = [H, J][lu[:p],:]
         LowSolveult!(lu[:L],HJ,nx+1,nx)
         UppSolveMult!(lu[:U],HJ,nx+1,nx)
         H[:] = HJ[:,1]
         J[:] = HJ[:,2:(nx+1)]
end

"""
     Preconditioner(h,X,P;jac="User")

Directly applies inverse preconditioning matrix.
"""
function Preconditioner(h,X,P;jac="User")
  J = h(X,P)
  #println("J:   ",J)
  if (length(X)>1)
    #println("mid.(J)", mid.(J))
    Y = inv(mid.(J))
  else
    Y = [1.0/(mid(J[1]))]
  end
  return Y
end
