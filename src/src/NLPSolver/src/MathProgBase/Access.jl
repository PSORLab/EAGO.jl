MathProgBase.status(m::EAGO_NLP_Model) = m.status
MathProgBase.getsolution(m::EAGO_NLP_Model) = m.BnBModel.soln
MathProgBase.getobjval(m::EAGO_NLP_Model) = m.BnBModel.UBDg
MathProgBase.getobjbound(m::EAGO_NLP_Model) = m.BnBModel.UBDg
getfeasibility(m::EAGO_NLP_Model) = m.BnBModel.feas_fnd
MathProgBase.getsolvetime(m::EAGO_NLP_Model) = LBDtime(m.BnBModel) + UBDtime(m.BnBModel)
MathProgBase.setwarmstart!(m::EAGO_NLP_Model,x) = []
MathProgBase.setsense!(m::EAGO_NLP_Model,sense) = (m.Opts.sense = sense)
MathProgBase.getsense(m::EAGO_NLP_Model) = m.Opts.sense

"""
    set_Implicit_Model!(m::EAGO.EAGO_NLP_Model,...)

Sets the functions and dimensionality of the implicit model corresponding to
the EAGO_NLP_model.
"""
function set_Implicit_Model!(m::EAGO.EAGO_NLP_Model,f::Function,h::Function,
                             hj::Function,g::Function,nx::Int64)
    m.Imp_f = f
    m.Imp_g = g
    m.Imp_h = h
    m.Imp_hj = hj
    m.Imp_nx = nx
end

"""
    Solve_Implicit(jm::JuMP.Model)

Sets the functions and dimensionality of the implicit model corresponding to
the JuMP model containing the EAGO_NLP_model.
"""
function Solve_Implicit(jm::JuMP.Model,f::Function,h::Function,
                             hj::Function,g::Function,nx::Int64)

    # Builds the JuMP model and adds the appropriate equations
    JuMP.build(jm)
    m = internalmodel(jm)
    m.Opts.Imp_f = f
    m.Opts.Imp_g = g
    m.Opts.Imp_h = h
    m.Opts.Imp_hj = hj
    m.Opts.Imp_nx = nx

    # Unpacks inner model and solves
    MathProgBase.optimize!(m)

    # Repacks the JuMP model with appropriate info
    return MathProgBase.status(m)
end

#=
"""
    Save_History(m::EAGO_NLP_Model)
"""
function Save_History(f::String,m::EAGO_NLP_Model)
     if (EAGO_NLP_Model.status == :Optimal)
     elseif (EAGO_NLP_Model.status == :Infeasible)
     end
end

=#
