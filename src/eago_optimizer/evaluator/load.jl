function copy_to_function(N::Int, x::JuMP._FunctionStorage)
    lenx = length(x.nd)
    temp_set = fill(MC{N}(Interval(-Inf, Inf)), (lenx,))
    temp_flt = Array{Float64}(undef, lenx)
    temp_bool = Array{Bool}(undef, lenx)

    tpdict = Dict{Int,Tuple{Int,Int,Int,Int}}()
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

    FunctionSetStorage{N}(x.nd, x.adj, x.const_values, temp_set, temp_flt,
                              temp_bool, tp1storage, tp2storage, tp3storage, tp4storage, tpdict,
                              x.grad_sparsity, x.hess_I, x.hess_J, x.dependent_subexpressions)
end
function copy_to_subexpr(N::Int, x::JuMP._SubexpressionStorage)
    lenx = length(x.nd)
    temp_set = fill(MC{N}(Interval(-Inf, Inf)), (lenx,))
    temp_flt = Array{Float64}(undef, lenx)
    temp_bool = Array{Bool}(undef, lenx)

    tpdict = Dict{Int,Tuple{Int,Int,Int,Int}}()
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

    SubexpressionSetStorage{N}(x.nd, x.adj, x.const_values, temp_set, temp_flt,
                               temp_bool, tp1storage, tp2storage, tp3storage, tp4storage, tpdict,
                               x.linearity)
end
get_node(d::Evaluator) = d.current_node
