"""
    getsolution(x::BnBModel)

Returns the solution stored in the BnBModel.
"""
getsolution(x::BnBModel) = x.soln

"""
    getobjval(x::BnBModel)

Returns the objective value stored in BnBModel (global upper bound).
"""
getobjval(x::BnBModel) = x.UBDg

"""
    getobjbound(x::BnBModel)

Returns the objective value stored in BnBModel (global upper bound).
"""
getobjbound(x::BnBModel) = x.UBDg

"""
    getfeasibility(x::BnBModel)

Returns feasibility of problem (feasible point found?).
"""
getfeasibility(x::BnBModel) = x.feas_fnd

"""
    LBDtime(x::BnBModel)

Returns time spent solving lower bounding problem.
"""
LBDtime(x::BnBModel) = x.LBDgtime[end]

"""
    UBDtime(x::BnBModel)

Returns time spent solving upper bounding problem.
"""
UBDtime(x::BnBModel) = x.UBDgtime[end]
