# states used in depth first search
@enum VISIT_STATUS VISITED UNVISITED INPROGRESS

# performs a depth-first search label nodes with distance from expression
function dfs_variable(x::Tape)
    status = Dict{Int,VISIT_STATUS}()
    parent = Dict{Int,Int}(1 => -1)
    for i in 1:length(x.nd)
        status[i] = UNVISITED
    end
    dfs_variable_kernel!(x,x.nd[1],1,status,parent)
    return parent
end

# kernel for depth first search
function dfs_variable_kernel!(x::Tape, n::Tracer.NodeData, indx, status, parent)
    first_child = children(n)[1]
    if (first_child == -1) || (first_child == -2)
        return
    else
        status[indx] = INPROGRESS
        for c in children(n)
            if status[c] == UNVISITED
                parent[c] = indx
                child_node = x.nd[c]
                dfs_variable_kernel!(x, child_node, c, status, parent)
            end
        end
        status[indx] = VISITED
    end
end

# takes a tracer tape, rearranges it to have the terminal node at the first
# index, replacing the array of child with an one elemnt array containing
# parents (PASSING)
function child_to_parent!(x::Tape)
    println("input x: $x")

    # reverse array order
    idx_arr = [i for i in x.set_trace_count:-1:1]
    nd_rev_temp = reverse(x.nd)
    nd_rev = Tracer.NodeData[]
    for (i,node) in enumerate(nd_rev_temp)
        if node.children[1] == -2
            push!(nd_rev,node)
        elseif node.children[1] == -1
            push!(nd_rev,node)
        else
            child_temp = [idx_arr[node.children[j]] for j in 1:length(node.children)]
            nd_temp = NodeData(node.nodetype, node.index, child_temp)
            push!(nd_rev,nd_temp)
        end
    end
    x.nd = nd_rev

    println("reversed x: $x")
    num_valued = Dict{Int,Bool}()
    for i in 1:x.set_trace_count
        num_valued[i] = x.num_valued[x.set_trace_count-i+1]
    end
    x.num_valued = num_valued

    # rewrite node list to have only parents not children
    parent_dict = dfs_variable(x)
    println("parent_dict: $parent_dict")

    # delete disconnected nodes and shift overs as appropriate
    # labels nodes to be deleted
    to_delete = fill(true,(length(x.nd)))
    keyvals = [k for k in keys(parent_dict)]
    value = [v for v in values(parent_dict)]
    for k in keyvals
        to_delete[k] = false
    end
    println("to_delete: $to_delete")

    println("GOOD TO HERE")

    # just redo the below section
    #=
    count_shift = 0
    to_shift = zeros(Int,length(x.nd))
    for i in 1:length(to_shift)
        if to_delete[i]
            count_shift += 1
        end
        to_shift[i] = count_shift
    end
    println("to_shift: $to_shift")
    println("LOOKS VERY WRONG")

    keyshift = Dict{Int,Int}()
    valshift = Dict{Int,Int}()
    valshift[-1] = -1

    for k in keys(parent_dict)
        v = parent_dict[k]
        keyshift[k] = k - to_shift[k]
        if v > 0
            valshift[v] = v - to_shift[v]
        end
    end

    println("keyshift: $(keyshift)")
    println("valshift: $(valshift)")

    cut_dictionary = Dict{Int,Int}()
    for k in keys(parent_dict)
        cut_dictionary[keyshift[k]] = valshift[parent_dict[k]]
    end

    println("cut_dictionary: $(cut_dictionary)")
    =#

    # Recreates the node list
    nd_list = []
    count = 1
    for i in 1:length(parent_dict)
        if ~to_delete[i]
            index = x.nd[i].index
            temp_node = Tracer.NodeData(x.nd[i].nodetype, index, [cut_dictionary[count]])
            push!(nd_list, temp_node)
            count += 1
        end
    end
    x.nd = nd_list
end


"""
    get_component_tape

Take a list of tapes for the i to n components.
"""
function get_component_tapes(g::Function, ng::Int, nx::Int)

    Nodes = Vector{JuMP.NodeData}[]
    Values = Vector{Float64}[]

    g_plus = x -> g(x) .+ [i for i in 1:ng]
    g_tape = Tracer.trace_script(g, nx);
    gp_tape = Tracer.trace_script(g_plus, nx);

    # gets list of terminal nodes
    terminal_nodes = zeros(Int,ng)
    position = ng
    for i in length(gp_tape.nd):-2:(length(g_tape.nd)+1)
        terminal_nodes[position] = gp_tape.nd[i].children[1]
        position -= 1
    end

    temp_tape_1 = []
    temp_tape_2 = []

    for i in 1:ng
        tape_i = Tape(g_tape.nd[1:terminal_nodes[i]], g_tape.const_values,
                      g_tape.num_valued, terminal_nodes[i], g_tape.const_count)
        (i == 1) && (temp_tape_1 = deepcopy(tape_i))
        Tracer.child_to_parent!(tape_i)
        (i == 1) && (temp_tape_2 = deepcopy(tape_i))
        push!(Nodes, convert.(JuMP.NodeData,tape_i.nd))
        push!(Values, tape_i.const_values)
    end
    return Nodes, Values, temp_tape_1, temp_tape_2
end

function MOI.initialize(d::JuMP.NLPEvaluator, f, g, ng, nx, requested_features::Vector{Symbol}, def_nlexpr)

    # get expression graph representation of the objective functions
    f_tape = Tracer.trace_script(f,nx)
    Tracer.child_to_parent!(f_tape)
    f_jnd = convert.(JuMP.NodeData, f_tape.nd)
    nlobj = JuMP._NonlinearExprData(f_jnd, f_tape.const_values); d.has_nlobj = true

    # get expression graph representation for constraints
    nlconstr_list = JuMP._NonlinearConstraint[]
    if isa(g, Function)
        gi_nd, gi_values = get_component_tapes(g, ng, nx);
        for i in 1:ng
            push!(nlconstr_list, JuMP._NonlinearConstraint(JuMP._NonlinearExprData(gi_nd[i], gi_values[i]), -Inf, 0.0))
        end
    end

    # get expression graph representation of nlexpr objects

    # checks for user_mv_operators
    has_user_mv_operator = false
    has_user_mv_operator |= JuMP._Derivatives.has_user_multivariate_operators(f_jnd)
    for nlexpr in def_nlexpr
        has_user_mv_operator |= JuMP._Derivatives.has_user_multivariate_operators(nlexpr.nd)
    end
    for nlconstr in nlconstr_list
        has_user_mv_operator |= JuMP._Derivatives.has_user_multivariate_operators(nlconstr.terms.nd)
    end
    d.disable_2ndorder = has_user_mv_operator

    for feat in requested_features
        if !(feat in MOI.features_available(d))
            error("Unsupported feature $feat")
        end
    end
    if d.eval_objective_timer != 0
        return
    end

    # assumes MOI index and consecutive index are the same as x[1:n] are the only variables used
    moi_index_to_consecutive_index = Dict(MOI.VariableIndex(i) => i for i in 1:nx)
    largest_user_input_dimension = nx

    d.user_output_buffer = Array{Float64}(undef,largest_user_input_dimension)
    d.jac_storage = Array{Float64}(undef,max(nx,largest_user_input_dimension))

    d.constraints = JuMP._FunctionStorage[]
    d.last_x = fill(NaN, nx)

    d.parameter_values = Float64[]

    d.want_hess = (:Hess in requested_features)
    want_hess_storage = (:HessVec in requested_features) || d.want_hess
    coloring_storage = JuMP._Derivatives.Coloring.IndexedSet(nx)

    max_expr_length = 0
    main_expressions = Array{Vector{JuMP.NodeData}}(undef,0)
    subexpr = Array{Vector{JuMP.NodeData}}(undef,0)

    for nlexpr in def_nlexpr
        push!(subexpr, nlexpr.nd)
    end
    if d.has_nlobj
        push!(main_expressions, f_jnd)
    end
    for nlconstr in nlconstr_list
        push!(main_expressions, nlconstr.terms.nd)
    end
    d.subexpression_order, individual_order = JuMP.order_subexpressions(main_expressions,subexpr)

    d.subexpression_linearity = Array{JuMP.Linearity}(undef,length(def_nlexpr))
    subexpression_variables = Array{Vector{Int}}(undef,length(def_nlexpr))
    subexpression_edgelist = Array{Set{Tuple{Int,Int}}}(undef,length(def_nlexpr))
    d.subexpressions = Array{JuMP._SubexpressionStorage}(undef,length(def_nlexpr))
    d.subexpression_forward_values = Array{Float64}(undef,length(d.subexpressions))
    d.subexpression_reverse_values = Array{Float64}(undef,length(d.subexpressions))

    # Load all subexpressions that are used
    empty_edgelist = Set{Tuple{Int,Int}}()
    for k in d.subexpression_order
        d.subexpression_forward_values[k] = NaN
        d.subexpressions[k] = JuMP._SubexpressionStorage(def_nlexpr[k].nd, def_nlexpr[k].const_values, nx, d.subexpression_linearity, moi_index_to_consecutive_index)
        subex = d.subexpressions[k]
        d.subexpression_linearity[k] = subex.linearity
        if d.want_hess
            empty!(coloring_storage)
            JuMP.compute_gradient_sparsity!(coloring_storage,subex.nd)
            for idx in JuMP.list_subexpressions(subex.nd)
                union!(coloring_storage, subexpression_variables[idx])
            end
            subexpression_variables[k] = collect(coloring_storage)
            empty!(coloring_storage)
            linearity = JuMP.classify_linearity(subex.nd, subex.adj, d.subexpression_linearity)
            edgelist = JuMP.compute_hessian_sparsity(subex.nd, subex.adj, linearity, coloring_storage, subexpression_edgelist, subexpression_variables)
            subexpression_edgelist[k] = edgelist
        end
    end

    if :ExprGraph in requested_features
        error("Script optimizer currently doesn't support ")
    end

    max_chunk = 1

    nd = main_expressions[1]
    d.objective = JuMP._FunctionStorage(nd, f_tape.const_values, nx, coloring_storage, d.want_hess, d.subexpressions, individual_order[1], d.subexpression_linearity, subexpression_edgelist, subexpression_variables, moi_index_to_consecutive_index)
    max_expr_length = max(max_expr_length, length(d.objective.nd))
    max_chunk = max(max_chunk, size(d.objective.seed_matrix,2))

    for k in 1:length(nlconstr_list)
        nlconstr = nlconstr_list[k]
        idx = (d.has_nlobj) ? k+1 : k
        nd = main_expressions[idx]
        push!(d.constraints, JuMP._FunctionStorage(nd, nlconstr.terms.const_values, nx, coloring_storage, d.want_hess, d.subexpressions, individual_order[idx], d.subexpression_linearity, subexpression_edgelist, subexpression_variables, moi_index_to_consecutive_index))
        max_expr_length = max(max_expr_length, length(d.constraints[end].nd))
        max_chunk = max(max_chunk, size(d.constraints[end].seed_matrix,2))
    end

    # hard coded chunk upper bound
    max_chunk = min(max_chunk, 10)

    if d.want_hess || want_hess_storage # storage for Hess or HessVec
        d.input_ϵ = Array{Float64}(undef,max_chunk*num_variables_)
        d.output_ϵ = Array{Float64}(undef,max_chunk*num_variables_)
        d.forward_storage_ϵ = Array{Float64}(undef,max_chunk*max_expr_length)
        d.partials_storage_ϵ = Array{Float64}(undef,max_chunk*max_expr_length)
        d.reverse_storage_ϵ = Array{Float64}(undef,max_chunk*max_expr_length)
        d.subexpression_forward_values_ϵ = Array{Float64}(undef,max_chunk*length(d.subexpressions))
        d.subexpression_reverse_values_ϵ = Array{Float64}(undef,max_chunk*length(d.subexpressions))
        for k in d.subexpression_order
            subex = d.subexpressions[k]
            subex.forward_storage_ϵ = zeros(Float64,max_chunk*length(subex.nd))
            subex.partials_storage_ϵ = zeros(Float64,max_chunk*length(subex.nd))
            subex.reverse_storage_ϵ = zeros(Float64,max_chunk*length(subex.nd))
        end
        d.max_chunk = max_chunk
        if d.want_hess
            d.hessian_sparsity = JuMP._hessian_lagrangian_structure(d)
        end
    end

    # reset timers
    d.eval_objective_timer = 0
    d.eval_objective_gradient_timer = 0
    d.eval_constraint_timer = 0
    d.eval_constraint_jacobian_timer = 0
    d.eval_hessian_lagrangian_timer = 0

    return nlobj, nlconstr_list, def_nlexpr
end
