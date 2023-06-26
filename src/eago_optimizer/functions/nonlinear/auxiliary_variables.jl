
function _not_EAGO_error!(m::JuMP.Model)
    if JuMP.solver_name(m) !== "EAGO: Easy Advanced Global Optimization"
        error("Solver attached to model must be EAGO.Optimizer")
    end
end

# Reference for auxiliary variables
struct AuxiliaryVariableRef <: JuMP.AbstractVariableRef
    idx::Int
    model::JuMP.Model
end

Base.@kwdef mutable struct _AuxVarData
    aux::Dict{JuMP.VariableRef, AuxiliaryVariableRef} = Dict{JuMP.VariableRef, AuxiliaryVariableRef}()
    #mimo_expr::Vector{MIMOExpr} = MIMOExpr[]
    last_hook::Union{Nothing,Function} = nothing
end
function is_auxiliary_variable(m::_AuxVarData, i::Int)
    false
end
is_auxiliary_variable(::Nothing, ::Int) = false


#=
function aux_variable_optimizehook(model::JuMP.Model)
    initialize_auxiliary_variables!(model)
    model.optimize_hook = model.ext[:aux_var].last_hook
    optimize!(model)
    model.optimize_hook = aux_variable_optimizehook
    return
end
function _initialize_auxiliary_variable_data(model::Model)
    model.ext[:aux_var] = _AuxVarData()
    model.ext[:aux_var].last_hook = model.optimize_hook 
    model.optimize_hook = aux_variable_optimizehook
    return model
end
function enable_auxiliary_variables(model::Model)
    haskey(model.ext, :aux_var) && error("Model has auxiliary variables parameter enabled")
    return _initialize_auxiliary_variable_data(model)
end
EAGOModel(m::JuMP.Model) = enable_auxiliary_variables(m)
EAGOModel() = EAGOModel(Model(EAGO.Optimizer))


_getmodel(v::AuxiliaryVariableRef) = v.model
function _getauxdata(v::AuxiliaryVariableRef)::_AuxVarData
    return _getauxdata(_getmodel(v))::_AuxVarData
end
function _getauxdata(model::Model)::_AuxVarData
    auxvar = get(model.ext, :aux_var, nothing)
    if auxvar !== nothing
        return auxvar
    end
    return enable_auxiliary_variables(model)
end


JuMP.index(p::AuxiliaryVariableRef) = p.idx
Base.iszero(::AuxiliaryVariableRef) = false
Base.copy(p::AuxiliaryVariableRef) = ParameterRef(p.idx, p.model)

struct AuxiliaryVariableNotOwned <: Exception
    aux::AuxiliaryVariableRef
end
JuMP.owner_model(p::AuxiliaryVariableRef) = p.model
function JuMP.check_belongs_to_model(p::AuxiliaryVariableRef, model::AbstractModel)
    if owner_model(p) !== model
        throw(AuxiliaryVariableNotOwned(p))
    end
end
function JuMP.is_valid(model::Model, aux::AuxiliaryVariableRef)
    return model === owner_model(aux)
end

function Base.hash(p::AuxiliaryVariableRef, h::UInt)
    return hash(objectid(owner_model(p)), hash(p.ind, h))
end
function Base.isequal(p1::AuxiliaryVariableRef, p2::AuxiliaryVariableRef)
    return owner_model(p1) === owner_model(p2) && p1.ind == p2.ind
end

JuMP.name(p::AuxiliaryVariableRef) = get(_getauxdata(p).names, p, "")
JuMP.set_name(p::AuxiliaryVariableRef, s::String) = _getauxdata(p).names[p] = s


struct AuxVar end

_aux_msg(msg) = "Invalid initialization of auxiliary variable. " * msg * " not supported."
function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, ::AuxVar)
    info.has_lb    && _error(_aux_msg("Lower bound"))
    info.has_ub    && _error(_aux_msg("Upper bound"))
    info.binary    && _error(_aux_msg("Binary"))
    info.integer   && _error(_aux_msg("Integer"))
    info.has_start && _error(_aux_msg("Initial value"))
    info.has_fix   && _error(_aux_msg("Fixed value"))
    return AuxVar()
end
function JuMP.add_variable(m::JuMP.Model, v::AuxVar, name::String="")
    vref = _add_auxiliary_variable(m)
    if !isempty(name)
        JuMP.set_name(vref, name)
    end
    return vref
end

macro auxiliary_variable(m, args...)
    esc(quote 
            @show ($args)
            @show ($args...)
            @show $(args...)
            #@variable($m, ($args...), AuxVar())
        end
    )
end

macro mimo_expression(m, name, f!, y, x) end

macro mimo_expression(m, f!, y, x) 
    esc(quote
            name = Symbol($f!)
            @mimo_expression($m, name, $f!, $y, $x)
        end
    )
end
=#

#=
struct MIMOExpr
    m::Model
    f!::Function
    y::Vector{Union{AuxiliaryVariableRef, JuMP.VariableRef}}
    x::Vector{Union{AuxiliaryVariableRef, JuMP.VariableRef}}
end


function initialize_auxiliary_variables!(m::JuMP.Model)
    set_optimizer_attribute(m, "_auxiliary_variable_info", m.ext[:aux_var])
    return
end

function aux_variable_optimize!(m::JuMP.Model)
    initialize_auxiliary_variables!(m)
    m.optimize_hook = m.ext[:aux_var].last_hook
    optimize!(m)
    m.optimize_hook = aux_variable_optimize!
    return
end

function JuMP.add_variable(m::JuMP.Model, v::AuxVar, name::String="")
    x = JuMP.add_variable(m, JuMP.ScalarVariable(v.info), name)
    m.ext[:aux_var].var_to_aux[x] = v
    return x
end

JuMP.build_variable(_err::Function, info::JuMP.VariableInfo, ::Type{AuxVar}; kwargs...) = AuxVar(info)
function JuMP.add_variable(m::JuMP.Model, v::AuxVar, name::String)
    x = JuMP.add_variable(m, JuMP.ScalarVariable(v.info), name)
    m.ext[:aux_var].var_to_aux[x] = v
    return x
end

"
m = EAGOModel()
@variable(m, -2 <= x[i=1:2] <= 2)
"
function add_auxiliary_variable(m::JuMP.Model, x, l, u)
end

"""
"""
function add_mimo_expression(m::JuMP.Model, name::String, f!::Function, y, x)
    MIMOExpr
    m = mimo_expr
end
=#