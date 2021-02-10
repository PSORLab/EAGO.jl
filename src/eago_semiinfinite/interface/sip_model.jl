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

Base.@kwdef mutable struct SIPData
    p::Dict{JuMP.VariableRef,Bool}                  = Dict{JuMP.VariableRef,Bool}()
    constr_type::Dict{JuMP.ConstraintRef, SIPCons}  = Dict{JuMP.ConstraintRef, SIPCons}()
    nl_constraint_type::Dict{Int, SIPCons}          = Dict{Int, SIPCons}()
    nl_expression_type::Dict{Int, SIPCons}          = Dict{Int, SIPCons}()
end

function initialize_sip_data(m::JuMP.Model)
    m.ext[:sip] = SIPData()
end

function sip_optimizehook(m::JuMP.Model; kwargs...)
    data = _getsipdata(m)::SIPData
    if uncertain_variable_num(data) == 0.0
        ret = JuMP.optimize!(m::JuMP.Model, ignore_optimize_hook = true, kwargs...)
    else
        # ret  = model_sip_solve(m)
    end
    return ret
end

function enable_semiinfinite(m::JuMP.Model)
    haskey(m.ext, :sip) && error("Model already has semiinfinite programs enabled")
    initialize_sip_data(m)
    JuMP.set_optimize_hook(m, sip_optimizehook)
    return nothing
end

function ModelWithSIP(args...; kwargs...)
    m = JuMP.Model(args...; kwargs...)
    enable_semiinfinite(m)
    return m
end

"""
    @uncertain_variable

Add an *anonymous* variable to `m::SIPModel` described by the keyword
arguments `kw_args` and returns the variable.

    @uncertain_variable(m, expr, args..., kw_args...)

This macro simply denotes the variable as belonging to the set of uncertain
variables and then performs `@variable(m, expr, args..., kw_args...)` to
add the variable to the basic JuMP model for storage. All JuMP syntax is
supported.
"""
macro uncertain_variable(args...)
    esc(quote
        vi = @variable($(args...))
        $(args[1]).ext[:sip].p[vi] = true
        vi
        end)
end

"""
    @decision_variable

Add an *anonymous* variable to `m::SIPModel` described by the keyword
arguments `kw_args` and returns the variable.

    @decision_variable(m, expr, args..., kw_args...)

This macro simply denotes the variable as belonging to the set of decision
variables and then performs `@variable(m, expr, args..., kw_args...)` to
add the variable to the basic JuMP model for storage. All JuMP syntax is
supported.
"""
macro decision_variable(args...)
    esc(quote
        vi = @variable($(args...))
        $(args[1]).ext[:sip].p[vi] = false
        vi
        end)
end
