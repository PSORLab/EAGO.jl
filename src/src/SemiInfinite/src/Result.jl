"""
--------------------------------------------------------------------------------
Type: SIP_result
--------------------------------------------------------------------------------
Description:
Composite type for storing the resulting form SIP solution routine.
--------------------------------------------------------------------------------
Fields:
k             Int64 - Number of iterations run
UBD           Float64 - Upper bound of SIP solution (optimal value)
LBD           Float64 - Lower bound of SIP solution
feas          Bool - Feasibility of SIP
LBP_time      Float64 - Time spent solving the lower bounding problem (sec)
LLP_time      Float64 - Time spent solving the lower level problem (sec)
UBP_time      Float64 - Time spent solving the upper bounding problem (sec)
xbar          Array{Float64} - Solution point
--------------------------------------------------------------------------------
"""
type SIP_result
    k::Int64
    UBD::Float64
    LBD::Float64
    feas::Bool
    LBP_time::Float64
    LLP_time::Float64
    UBP_time::Float64
    xbar
end
SIP_result() = SIP_result(1,Inf,-Inf,true,0.0,0.0,0.0,[])
