abstract type AbstractScheme end
#=
"""
Relaxation Scheme Linear is a function that takes a linear constraint and either
adds a relaxed constraint to an NLP evaluation block or to the optimizer directly.
"""
=#
mutable struct RelaxationScheme <: AbstractScheme
    optimizer_type::Symbol
    linear_relaxation::Function      # true -> NLP, false -> constraint
    linear_relaxation_check::Function
    quadratic_relaxation::Function
    quadratic_relaxation_check::Function
    expression_library::Dict{Symbol,Function}
    linear_relaxed_after_load::Dict{Int,Bool}
    quadratic_relaxed_after_load::Dict{Int,Bool}
    nonlinear_relaxed_after_load::Dict{Int,Bool}
end

macro DefineDefaultLinear(set)
    expr = quote
                default_linear(f::MOI.ScalarAffineFunction{Float64},s::$set,n::NodeBB) = (f,s)
                default_linear_check(f::MOI.ScalarAffineFunction{Float64},s::$set,n::NodeBB) = false
           end
    return esc(expr)
end

@DefineDefaultLinear(MOI.GreaterThan{Float64})
@DefineDefaultLinear(MOI.LessThan{Float64})
@DefineDefaultLinear(MOI.EqualTo{Float64})

default_relaxation_scheme() = RelaxationScheme(:LP,default_linear,default_linear_check,
                                             x->x, x->x,
                                             Dict{Symbol,Function}(),
                                             Dict{Int,Bool}(),
                                             Dict{Int,Bool}(),
                                             Dict{Int,Bool}())
