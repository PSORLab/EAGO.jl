## Functions for generating console displayed
```@docs
print_int!(B::BnBSolver,k_int::Int64,k_nod::Int64,nid::Int64,lbdp::Float64,lbd::Float64,ubd::Float64,feasL::Bool,feasU::Bool)
print_results!(B::BnBSolver,sol::Float64,pnt,feas::Bool,lbd_bool::Bool)
print_node!(x::BnBSolver,id::Int64,lbd::Float64,box)
print_sol!(x::BnBSolver,y::BnBModel,ubdcnt::Int64,lbdcnt::Int64,ubdtime::Float64,lbdtime::Float64)
```
