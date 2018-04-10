"""
    BM_breadth!(S::BnBSolver,B::BnBModel,tL::Float64,tU::Float64,
                X1::Vector{Interval{Float64}},X2::Vector{Interval{Float64}},
                pos::Int64)

Stores two interval boxes `X1,X2` to the bottom of the stack along with their
respective lower, `tL` and upper bounds, `tU` and their position number in the
BnB tree. Also, assigns node numbers.
"""
function BM_breadth!(S::BnBSolver,B::BnBModel{T},tL::Float64,tU::Float64,
                             X1::Vector{T},
                             X2::Vector{T},pos::Int64) where {T}
  unshift!(B.box,X1,X2)
  unshift!(B.LBD,tL,tL)
  unshift!(B.UBD,tU,tU)
  unshift!(B.id,B.max_id+1,B.max_id+2)
  unshift!(B.pos,pos+1,pos+1)
  B.max_id += 2
end

"""
    BM_depth_best!(S::BnBSolver,B::BnBModel,tL::Float64,tU::Float64,
                X1::Vector{Interval{Float64}},X2::Vector{Interval{Float64}},
                pos::Int64)

Stores two interval boxes `X1,X2` to the top of the stack along with their
respective lower, `tL` and upper bounds, `tU` and their position number in the
BnB tree. Also, assigns node numbers.
"""
function BM_depth_best!(S::BnBSolver,B::BnBModel{T},tL::Float64,tU::Float64,
                             X1::Vector{T},
                             X2::Vector{T},pos::Int64) where {T}
  push!(B.box,X1,X2)
  push!(B.LBD,tL,tL)
  push!(B.UBD,tU,tU)
  push!(B.id,B.max_id+1,B.max_id+2)
  push!(B.pos,pos+1,pos+1)
  B.max_id += 2
end

"""
    BM_Single!(S::BnBSolver,B::BnBModel,tL::Float64,tU::Float64,
                             X::Vector{Interval{Float64}},pos::Int64)

Stores interval box `X` to the top of the stack along with their
respective lower, `tL` and upper bounds, `tU` and their position number in the
BnB tree.
"""
function BM_Single!(S::BnBSolver,B::BnBModel{T},tL::Float64,tU::Float64,
                             X::Vector{T},pos::Int64) where {T}
  push!(B.box,X)
  push!(B.LBD,tL)
  push!(B.UBD,tU)
  push!(B.id,B.max_id)
  push!(B.pos,pos)
end
