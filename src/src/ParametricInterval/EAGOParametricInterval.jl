#=
type Param_Bisect_Opts
  DAGflag::Bool
  LPflag::Bool
  kmax_main
  kmax_cntr
  style
  display
  ptol
  etol
  rtol
  DAGpass
  p_rel_bisect
  DAGh
  DAGg
  DAGsym
end
Param_Bisect_Opts() = Param_Bisect_Opts(false, #DAGflag::Bool
                                        false, #LPflag::Bool
                                        1E2, #kmax_main
                                        1E2, #kmax_cntr
                                        "KrawczykCW", #style
                                        "Full", #display
                                        1E-3, #ptol
                                        1E-3, #etol
                                        1E-6, #rtol
                                        5, #DAGpass
                                        false, # p_rel_bisect
                                        [],
                                        [],
                                        [])

include("lib/Parametric_Utility.jl")
include("lib/Parametric_Contractor.jl")
#include("lib/Parametric_Test.jl")
#include("lib/Parametric_Bisection.jl") # Used for generalized bisection
#include("lib/Parametric_Main.jl")
=#
#=
export Param_Bisect_Opts,setprec, Miranda, MirandaExc, partialIncTop,
       partialIncBot, Strict_XinY, isEqual, extDivide, extProcess,
=#

#=
export  PI_NewtonGS, PI_KrawczykCW,
        PIn_NewtonGS, PIn_KrawczykCW,
        PId_NewtonGS, PId_KrawczykCW,
        Sparse_Precondition!, SparseInSto
=#
include("src/contractor/PIntvParams.jl")

include("src/utils/general.jl")
include("src/utils/inclusion.jl")
include("src/utils/extdivision.jl")
include("src/utils/bandwidth.jl")

include("src/preconditioner/densebanded.jl")
include("src/preconditioner/denseblockdiag.jl")
include("src/preconditioner/precondition.jl")

include("src/iteration/Krawczyk.jl")
include("src/iteration/Newton.jl")

include("src/contractor/Contractor.jl")

export PIntvParams, Param_Intv_Contractor
