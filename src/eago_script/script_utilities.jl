function add_constant(x::Tape,y)
    x.set_trace_count += 1; x.const_count += 1
    node = NodeData(JuMP.VALUE,x.const_count,[-2]); push!(x.nd,node)
    x.num_valued[x.set_trace_count] = true
    push!(x.const_values,Float64(y))
    x.set_trace_count
end

function add_set_node!(x::Tape,node::NodeData)
    push!(x.nd,node)
    x.set_trace_count += 1
    x.num_valued[x.set_trace_count] = false
end
