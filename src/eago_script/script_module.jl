# Case 1 (Pass): sin(x[1]) + y * 3.0*2.0*x[2]
# Next Case: Assignment in Script (To Test)
# Next Case: For loop (To Test)
# Next Case: Through right assertion (To Test)
# Next Case: While loop
# Next Case: User-defined...
# Next Case: Multiple output

# Next Case: Through left assertion
# Next Case: If-else (NEXT BASED ON CONCOLIC FUZZER)

include("functional_control_flow.jl")

@reexport module Tracer

    import Base: abs, sin, cos, tan, sec, csc, cot, asin, acos, atan, asec, acsc,
              acot, sinh, cosh, tanh, asinh, acosh, atanh, sech, asech, csch,
              acsch, coth, acoth, sqrt, log, log2, log10, log1p, exp, exp2,
              exp10, expm1, +, -, step, sign, real, inv, *,
              +, -, / , min, max, >, <, ==, =>, <= ,^, getindex, afoldl

    using Cassette, JuMP, SparseArrays, MathOptInterface
    const MOI = MathOptInterface

    include("script_types.jl")
    include("script_utilities.jl")

    Cassette.@context TraceCtx

    include("execute_operations.jl")
    include("execute_utils.jl")

    # getindex(x,1) to getindex(x,n) correspond to 1..n variable
    function trace_script(f,n)
        tape = Tape(n)
        x = SetTraceSto([SetTrace(i) for i=1:n])
        Cassette.overdub(TraceCtx(metadata = tape), f, x)
        return tape
    end

    include("to_jump.jl")

    #export solve_script
    #include("solve_script.jl")

end # module
