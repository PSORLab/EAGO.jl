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
    for i in 1:lenx
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
                                x.grad_sparsity, x.dependent_subexpressions)
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
    for i in 1:lenx
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
