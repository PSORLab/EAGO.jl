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
const MOIU_AUTO_CACHE = MOIU.CachingOptimizerMode.AUTOMATIC

Base.@kwdef mutable struct SIPData
    p::Dict{JuMP.VariableRef,Bool}                  = Dict{JuMP.VariableRef,Bool}()
    constr_type::Dict{JuMP.ConstraintRef, SIPCons}  = Dict{JuMP.ConstraintRef, SIPCons}()
    nl_constraint_type::Dict{Int, SIPCons}          = Dict{Int, SIPCons}()
    nl_expression_type::Dict{Int, SIPCons}          = Dict{Int, SIPCons}()
    model_decision::Model                           = Model(caching_mode = MOIU_AUTO_CACHE)
    model_uncertain::Model                          = Model((aching_mode = MOIU_AUTO_CACHE)
    model_semiinfinite::Model                       = Model(caching_mode = MOIU_AUTO_CACHE)
end
uncertain_variable_num(d::SIPData) = any(x -> x, d.p)

function _copy_nlp_objective!(out::JuMP.Model, m::JuMP.Model, nlp_data::JuMP._NLPData)
    if nlp_data.nlobj !== nothing
        # TODO CHECK THAT THE FORMULATION IS PEACHY
        obj_expr = MOI.objective_expr(nlp_new.evaluator)
        # obj_expr_jref = substitute moi variables to jump refs...
        @NLobjective(m.ext[:sip].model_decision, Min, obj_expr_jref)
    end
    return nothing
end

function _copy_nlp_constraints!(m::JuMP.Model, nlp_data::JuMP._NLPData)
    for i = 1:length(nlp_data.nlconstr)
        con_expr = MOI.constraint_expr(nlp_data.evaluator, i)
    end
    return nothing
end

"""

Populates variables and constraints into one of three models based on classification:
1) store decision variables and constraints containing only decision variables in
`m.model_decision`, 2) store uncertain variables and constraints containing only
uncertain variables in `m.model_uncertain`, 3) store all variables and all
constraints containing both decision variables to `m.model_semiinfinite`.
"""
function initialize_pure_models!(m::JuMP.Model)
    if haskey(m.ext, :sip)

        # extract sip data and storage models
        sip_data = _get_sip_data(m)::SIPData
        m_dec = sip_data.model_decision
        m_unc = sip_data.model_uncertain
        m_sip = sip_data.model_semiinfinite

        # get index map
        indx_map_dec = MOI.copy_to(backend(m_dec), backend(m), copy_names = true)
        indx_map_unc = MOI.copy_to(backend(m_unc), backend(m), copy_names = true)
        indx_map_sip = MOI.copy_to(backend(m_sip), backend(m), copy_names = true)

        # copy extension data
        for (key, data) in m.ext
            if key != :sip
                m_dec.ext[key] = copy_extension_data(data, m_dec, m)
                m_unc.ext[key] = copy_extension_data(data, m_unc, m)
                m_sip.ext[key] = copy_extension_data(data, m_sip, m)
            end
        end

        # copy objects (TODO need to check on SIP, DEC, UNC)
        reference_map = ReferenceMap(new_model, index_map)
        for (name, value) in object_dictionary(model)
            new_model[name] = getindex.(reference_map, value)
        end

        if m.nlp_data !== nothing
            nlp_data = deepcopy(m.nlp_data)
            nlp_data.evaluator = NLPEvaluator(m)
            MOI.initialize(nlp_data.evaluator, Symbol[:ExprGraph])
            copy_nlp_objective!(m, nlp_data)
            copy_nlp_expressions!(m, nlp_data)
            copy_nlp_constraints!(m, nlp_data)
        end

        for (vi, is_uncertain) in m.ext.p
            model_semiinfinite
            model_decision_
    end
    return nothing
end


function initialize_sip_data(m::JuMP.Model)
    m.ext[:sip] = SIPData()
end

function sip_optimizehook(m::JuMP.Model; kwargs...)
    data = _getsipdata(m)::SIPData
    if uncertain_variable_num(data) == 0.0
        ret = JuMP.optimize!(m::JuMP.Model, ignore_optimize_hook = true, kwargs...)
    else
        initialize_pure_models!(m)
        ret = sip_solve(data)
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
    set_optimizer(SIPOptimizer)
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
