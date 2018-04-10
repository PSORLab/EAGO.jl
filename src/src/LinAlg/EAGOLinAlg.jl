
include("src/BandMatrix/TriangularSolve.jl")
include("src/BandMatrix/LU_Precond.jl")

include("src/BlockDiagMatrix/TriangularSolve.jl")
include("src/BlockDiagMatrix/LU_Precond.jl")

include("src/Dense/TriangularSolve.jl")
include("src/Dense/LU_Precond.jl")

include("src/Sparse/TriangularSolve.jl")
include("src/Sparse/LU_Precond.jl")

#=
"""
    Defines in-place preconditioner for use in ParametricInterval and McCormick schemes.
"""
function Precondition!(sym)
    if (sym == :Band)
    elseif (sym == :BlockDiag)
    elseif (sym == :Sparse)
    elseif (sym == :Dense)
    elseif (sym == :Inv)
    end
end
=#
