#=

Macro interface is inspired by InfiniteOpt.jl. Many of the code structures are
reused. Transitioning this over to a fully InfiniteOpt.jl in the future is desirable
but most of our use cases involve nonlinear functions so that InfiniteOpt.jl
nonlinear support is gating.

The general idea here is pretty simple. When we add a variable to the model
the user specifies whether or not it is a decision variable or an uncertain
variable we then classify constraints. See https://github.com/PSORLab/EAGO.jl/issues/61
for a discussion of the roadmap here.
=#

@enum(SIPCons, DECISION, UNCERTAIN, SEMIINFINITE, SIPNOTSET)

struct SIPModel <: JuMP.AbstractModel
    m::JuMP.Model
    x::Vector{JuMP.VariableRef}
    p::Vector{JuMP.VariableRef}
    constraint_type::Dict{JuMP.ConstraintRef, SIPCons}
    nl_constraint_type::Dict{Int, SIPCons}
    nl_expression_type::Dict{Int, SIPCons}
end
