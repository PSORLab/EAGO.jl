function ExtendDivStorage!(S::BnBSolver,B::BnBModel{T},tL::Float64,tU::Float64,
                           X1::Vector{V},
                           X2::Vector{V},pos::Int64,
                           LBDorig::Float64,UBDorig::Float64) where {T<:AbstractFloat,V<:AbstractInterval}
    if (S.opt[3].ed_flag)
        push!(B.box,X1,S.opt[3].box2)
        push!(B.LBD,tL,LBDorig)
        push!(B.UBD,tU,UBDorig)
        push!(B.id,B.max_id+1,B.max_id+2)
        push!(B.pos,pos+1,pos+1)
        B.max_id += 2
    else
        push!(B.box,X1,X2)
        push!(B.LBD,tL,tL)
        push!(B.UBD,tU,tU)
        push!(B.id,B.max_id+1,B.max_id+2)
        push!(B.pos,pos+1,pos+1)
        B.max_id += 2
    end
end
