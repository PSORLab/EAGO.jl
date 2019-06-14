module Check_McCormick

    using Compat
    using Compat.Test
    using EAGO, JuMP, StaticArrays, IntervalArithmetic

    include("operator_library.jl")
    include("reverse_operators.jl")
    include("mccormick_utilities.jl")
    include("implicit_routines.jl")

end
