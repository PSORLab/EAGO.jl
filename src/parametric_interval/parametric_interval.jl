include("src/contractor/parametric_interval_params.jl")

include("src/utils/general.jl")
include("src/utils/inclusion.jl")
include("src/utils/extdivision.jl")
include("src/utils/bandwidth.jl")

include("src/preconditioner/densebanded.jl")
include("src/preconditioner/denseblockdiag.jl")
include("src/preconditioner/precondition.jl")

include("src/iteration/krawczyk.jl")
include("src/iteration/newton.jl")

include("src/contractor/contractor.jl")
