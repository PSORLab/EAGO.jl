MathProgBase.status(m::EAGO_NLP_Model) = m.status
MathProgBase.getsolution(m::EAGO_NLP_Model) = m.BnBModel.soln
MathProgBase.getobjval(m::EAGO_NLP_Model) = m.BnBModel.UBDg
MathProgBase.getobjbound(m::EAGO_NLP_Model) = m.BnBModel.UBDg
getfeasibility(m::EAGO_NLP_Model) = m.BnBModel.feas_fnd
MathProgBase.getsolvetime(m::EAGO_NLP_Model) = LBDtime(m.BnBModel) + UBDtime(m.BnBModel)
MathProgBase.setwarmstart!(m::EAGO_NLP_Model,x) = []
MathProgBase.setsense!(m::EAGO_NLP_Model,sense) = (m.Opts.sense = sense)
MathProgBase.getsense(m::EAGO_NLP_Model) = m.Opts.sense

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
