"""
    PIntvParams

Option structure for parametric interval contractor with the following fields:
- LAlg::Symbol       # Indicates type of linear algebra for preconditioner (:None, :Dense)
- CTyp::Symbol       # Contractor type to be used (:None, :Dense)
- etol::T            # Equality tolerance for interval checks
- rtol::T            # Amount added to interval during extended division
- nx::Int64          # Dimension of state space
- tband::Integer     # Total bandwidth u+l+1
- kmax::Int64        # Number of iterations to take for contractor
"""
mutable struct PIntvParams{T}
    LAlg::Symbol       # Indicates type of linear algebra for preconditioner (:None, :Dense)
    CTyp::Symbol       # Contractor type to be used (:None, :Dense)
    etol::T            # Equality tolerance for interval checks
    rtol::T            # Amount added to interval during extended division
    nx::Integer          # Dimension of state space
    tband::Integer     # Total bandwidth u+l+1
    kmax::Integer        # Number of iterations to take for contractor
end
PIntvParams(nx::Integer,band::Integer) = PIntvParams(:Dense,:Newton,1E-12,1E-12,nx,band,5)
