export EAGO_Inner_NLP, EAGO_NLP_Model, EAGO_NLPSolver, ImplicitSolver, ImplicitOptions,
       Solve_Implicit, registerEAGO, Optimize_Script

include("src/Utilities/Composite_DR.jl")
include("src/Utilities/General.jl")
include("src/Utilities/Vector_Hessian.jl")
include("src/Utilities/BnBRoutines.jl")

include("src/MathProgBase/Solver.jl")
include("src/MathProgBase/Checks.jl")
include("src/MathProgBase/Model.jl")
include("src/MathProgBase/Load.jl")
include("src/MathProgBase/Optimize.jl")
include("src/MathProgBase/Access.jl")

include("src/Lower_Explicit/AlphaBB_LBD.jl")
include("src/Lower_Explicit/Interval_LBD.jl")
include("src/Lower_Explicit/LP_Relax_LBD.jl")
include("src/Lower_Explicit/Ipopt_LBD.jl")
include("src/Lower_Explicit/SNOPT_LBD.jl")
#include("src/Lower_Explicit/Quadratic_LBD.jl")

#include("src/Lower_Implicit/Imp_Interval_LBD.jl")
include("src/Lower_Implicit/LP_Relax_LBD_Imp.jl")
#include("src/Lower_Implicit/Imp_Ipopt_LBD.jl")
include("src/Lower_Implicit/SNOPT_LBD_Imp.jl")

include("src/Upper_Explicit/Imp_Interval_UBD.jl")
include("src/Upper_Explicit/Interval_UBD.jl")
include("src/Upper_Explicit/Ipopt_UBD.jl")
include("src/Upper_Explicit/SNOPT_UBD.jl")
include("src/Upper_Explicit/NLPSolve_UBD.jl")
