const CR = JuMP.ConstraintRef

function classify_sip!(mSIP::SIPModel, c::AbstractConstraint)
    error("Function classify_sip! not defined for a constraint of type = $(typeof(c))")
end

function update_status(status::SIPCons, term_status::SIPCons)
    if status == SIPNOTSET
        return term_status
    elseif (status == DECISION) && (term_status == DECISION)
        return DECISION
    elseif (status == UNCERTAIN) && (term_status == UNCERTAIN)
        return UNCERTAIN
    end
    return SEMIINFINITE
end

classify_sip(mSIP::SIPModel, v::VariableRef) = v ∈ mSIP.p ? UNCERTAIN : DECISION
function classify_sip(mSIP::SIPModel, expr::GenericAffExpr{T,VariableRef}) where T
    status = SIPNOTSET
    for term in linear_terms(expr)
        status = update_status(status, classify_sip(mSIP, term[1]))
        if status == SEMIINFINITE
            break
        end
    end
    return status
end
function classify_sip(mSIP::SIPModel, expr::GenericQuadExpr{T,VariableRef}) where T
    status = SIPNOTSET
    for term in quad_terms(expr)
        status = update_status(status, classify_sip(mSIP, term[2]))
        status = update_status(status, classify_sip(mSIP, term[3]))
        if status == SEMIINFINITE
            break
        end
    end
    for term in linear_terms(expr)
        status = update_status(status, classify_sip(mSIP, term[1]))
        if status == SEMIINFINITE
            break
        end
    end
    return status
end

for typ in (:(Vector{GenericQuadExpr{T,VariableRef}}),
            :(Vector{GenericAffExpr{T,VariableRef}}),
            :(Vector{VariableRef}))
    @eval function classify_sip(mSIP::SIPModel, expr::Vector{$typ{T,VariableRef}}) where T
        status = SIPNOTSET
        for ex in expr
            status = update_status(status, classify_sip(mSIP, ex))
            if status == SEMIINFINITE
                break
            end
        end
        return status
    end
end

function classify_sip!(mSIP::SIPModel, cr::ScalarConstraint{F,S}) where {F<:JuMP.AbstractJuMPScalar, S<:JuMP.AbstractScalarSet}
    mSIP.constraint_type[cr] = classify_sip(mSIP, constraint_object(cr).func)
    nothing
end

# TODO: Define for constraint type...
function classify_sip!(mSIP::SIPModel, cr::CR{M,C,S}) where {M<:JuMP.AbstractModel, C ,S<:JuMP.AbstractShape}
    mSIP.constraint_type[cr] = classify_sip(constraint_object(cr))
    nothing
end

function classify_sip!(mSIP::SIPModel)
    constraint_types = list_of_constraint_types(mSIP.m)
    for (func_typ, set_typ) in constraint_types
        constraint_refs = all_constraints(mSIP.m, func_typ, set_typ)
        for cr in constraint_refs
            classify_sip!(mSIP, cr)
        end
    end
    return nothing
end
