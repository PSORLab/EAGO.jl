"""
    parametric_interval_params

Option structure for parametric interval contractor with the following fields:
- `LAlg::Symbol`: Indicates type of linear algebra for preconditioner (:None, :Dense)
- `CTyp::Symbol`: Contractor type to be used (:None, :Dense)
- `etol::Float64`: Equality tolerance for interval checks
- `rtol::Float64:`: Amount added to interval during extended division
- `nx::Int`: Dimension of state space
- `tband::Int`: Total bandwidth u+l+1
- `kmax::Int`: Number of iterations to take for contractor
"""
mutable struct parametric_interval_params
    LAlg::Symbol       # Indicates type of linear algebra for preconditioner (:None, :Dense)
    CTyp::Symbol       # Contractor type to be used (:None, :Dense)
    etol::Float64            # Equality tolerance for interval checks
    rtol::Float64            # Amount added to interval during extended division
    nx::Int          # Dimension of state space
    tband::Int     # Total bandwidth u+l+1
    kmax::Int        # Number of iterations to take for contractor
end
parametric_interval_params(nx::Int,band::Int) = parametric_interval_params(:Dense,:Newton,1E-12,1E-12,nx,band,5)
