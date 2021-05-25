function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    # Add precompile statements here
    # Base.precompile(Tuple{typeof(MOI.initialize),NLPEvaluator,Vector{Symbol}})
    return
end
