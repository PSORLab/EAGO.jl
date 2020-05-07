# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/evaluator/load.jl
# Create FunctionSetStorage and SubexpressionSetStorage from JuMP nonlinear
# backend. We'll evenetually want to build this from an expression drag instead.
#############################################################################


function copy_to_function!(d::Evaluator{N,T}, i::Int64, y::JuMP._FunctionStorage) where {N,T<:RelaxTag}
    x = deepcopy(y)
    lenx = length(x.nd)
    temp_set = fill(MC{N,T}(Interval(-Inf, Inf)), (lenx,))
    temp_flt = Array{Float64}(undef, lenx)
    temp_bool = Array{Bool}(undef, lenx)

    tpdict = Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}()
    tp1_count = 0
    tp2_count = 0
    tp3_count = 0
    tp4_count = 0
    for i = 1:lenx
        op = x.nd[i].index
        if double_tp(op)
            tp1_count += 1
            tp2_count += 1
            tp3_count += 1
            tp4_count += 1
            tpdict[i] = (tp1_count, tp2_count, tp3_count, tp4_count)
        elseif single_tp(op)
            tp1_count += 1
            tp2_count += 1
            tpdict[i] = (tp1_count, tp2_count, -1, -1)
        end
    end
    tp1storage = zeros(tp1_count)
    tp2storage = zeros(tp2_count)
    tp3storage = zeros(tp3_count)
    tp4storage = zeros(tp4_count)

    sto = FunctionSetStorage{N,T}(x.nd, x.adj, x.const_values, temp_set, temp_flt,
                                temp_bool, tp1storage, tp2storage, tp3storage, tp4storage, tpdict,
                                x.grad_sparsity, x.hess_I, x.hess_J, x.dependent_subexpressions)
    if i == 1
        d.objective = sto
    else
        push!(d.constraints, sto)
    end
    return
end
function copy_to_subexpr!(d::Evaluator{N,T}, y::JuMP._SubexpressionStorage) where {N,T<:RelaxTag}
    x = deepcopy(y)
    lenx = length(x.nd)
    temp_set = fill(MC{N,T}(Interval(-Inf, Inf)), (lenx,))
    temp_flt = Array{Float64}(undef, lenx)
    temp_bool = Array{Bool}(undef, lenx)

    tpdict = Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}()
    tp1_count = 0
    tp2_count = 0
    tp3_count = 0
    tp4_count = 0
    for i = 1:lenx
        op = x.nd[i].index
        if double_tp(op)
            tp1_count += 1
            tp2_count += 1
            tp3_count += 1
            tp4_count += 1
            tpdict[i] = (tp1_count, tp2_count, tp3_count, tp4_count)
        elseif single_tp(op)
            tp1_count += 1
            tp2_count += 1
            tpdict[i] = (tp1_count, tp2_count, -1, -1)
        end
    end
    tp1storage = zeros(tp1_count)
    tp2storage = zeros(tp2_count)
    tp3storage = zeros(tp3_count)
    tp4storage = zeros(tp4_count)

    sto = SubexpressionSetStorage{N,T}(x.nd, x.adj, x.const_values, temp_set, temp_flt,
                                     temp_bool, tp1storage, tp2storage, tp3storage, tp4storage, tpdict,
                                     x.linearity)
    push!(d.subexpressions, sto)
end
get_node(d::Evaluator) = d.current_node
