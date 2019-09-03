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
function neg_objective!(d::Evaluator)
    if (d.has_nlobj)
        # shifts the adjacency matrix to introduce -f(x) as first element of nd array
        rowval = rowvals(d.objective.adj) .+ 1; pushfirst!(rowval, 2)
        colptr = copy(d.objective.adj.colptr) .+ 1; pushfirst!(colptr, 1)
        nzval = nonzeros(d.objective.adj); pushfirst!(nzval, true)
        m, n = size(d.objective.adj)
        d.objective.adj = SparseMatrixCSC{Bool,Int}(m+1, n+1, colptr, rowval, nzval)

        # shifts the node list (and parents)
        shift_nd = [JuMP.NodeData(JuMP.CALLUNIVAR, 2, 2)]
        for nd in d.objective.nd
            push!(shift_nd, JuMP.NodeData(nd.nodetype, nd.index, nd.parent+1))
        end
        d.objective.nd = shift_nd

        new_tpdict = Dict{Int,Tuple{Int,Int,Int,Int}}()
        for key in keys(d.tpdict)
            new_tpdict[key+2] = d.tpdict[key]
        end
        d.objective.tpdict = new_tpdict

        nvflag = length(d.objective.numvalued) > 0 ? d.objective.numvalued[1] : false
        pushfirst!(d.objective.numvalued, nvflag)
        pushfirst!(d.objective.numberstorage ,0.0)
        pushfirst!(d.objective.setstorage, zero(MC{eltype(d)}))
    end
end

"""
    build_nlp_evaluator

Builds the evaluator used to generate relaxations of the nonlinear equations
and constraints from a source model.
"""
function build_nlp_evaluator(N::Int, src::T, x::Optimizer, bool_flag::Bool) where {T<:MOI.AbstractNLPEvaluator}

    # Checks to see if nlp data block evaluator is present
    if ~isa(src, EAGO.EmptyNLPEvaluator)

        # Creates the empty evaluator
        d = Evaluator{N}(src.m)

        num_variables_ = JuMP.num_variables(d.m)
        d.variable_number = num_variables_
        nldata::JuMP._NLPData = deepcopy(d.m.nlp_data)

        # Copy state of user-defined multivariable functions
        d.has_user_mv_operator = src.disable_2ndorder
        d.parameter_values = nldata.nlparamvalues
        d.last_x = fill(NaN, d.variable_number)
        d.last_node = NodeBB()

        # Set valued operators must contain a (sub)gradient and/or (sub)hessian field to support 1st/2nd order eval
        d.disable_1storder = false
        d.disable_2ndorder = true

        # Add objective functions, constraints, subexpressions
        d.has_nlobj = isa(nldata.nlobj, JuMP._NonlinearExprData)
        if (src.has_nlobj)
            d.objective = copy_to_function(N, src.objective)
        end

        for i in 1:length(src.constraints)
            push!(d.constraints, copy_to_function(N, src.constraints[i]))
        end

        d.subexpression_order = src.subexpression_order
        d.subexpression_linearity = src.subexpression_linearity
        d.subexpressions_as_julia_expressions = Any[]
        if isdefined(src, :subexpressions_as_julia_expressions)
            d.subexpressions_as_julia_expressions = src.subexpressions_as_julia_expressions
        end

        d.subexpression_values_set = MC{N}[]
        d.subexpression_values_flt = Float64[]
        d.subexpressions = SubexpressionSetStorage{N}[]
        for i in 1:length(src.subexpressions)
            temp = copy_to_subexpr(N, src.subexpressions[i])
            push!(d.subexpressions,temp)
        end
        d.subexpression_values_set = fill(NaN, length(d.subexpressions))
        d.subexpression_values_flt = fill(NaN, length(d.subexpressions))

        # Add bounds to evaluator
        for bnds in x.nlp_data.constraint_bounds
            push!(d.constraints_lbd, bnds.lower)
            push!(d.constraints_ubd, bnds.upper)
        end

        # USER OUTPUT BUFFERS??????
        d.cp_tolerance = x.cp_interval_tolerance
        d.cp_reptitions = x.cp_interval_reptitions
        d.has_reverse = x.evaluation_reverse
        d.subgrad_tighten = x.subgrad_tighten
        d.subgrad_tighten_reverse = x.subgrad_tighten_reverse
        d.jac_storage = Array{Float64}(undef, max(num_variables_, d.m.nlp_data.largest_user_input_dimension))

        d.constraint_number = length(d.constraints)
        d.subexpression_number = length(d.subexpressions)

        # calculate an index for each variable via search on
        d.index_to_variable = fill((-1,-1,-1), (d.variable_number,))
        for (oindx,node) in enumerate(d.objective.nd)
            if (node.nodetype == JuMP._Derivatives.VARIABLE)
                current_value = d.index_to_variable[node.index]
                if (current_value[1] == current_value[2] == current_value[3] == -1)
                    d.index_to_variable[node.index] = (oindx, 1, 1)
                end
            end
        end
        for (cindx,constraint) in enumerate(d.constraints)
            for (indx,node) in enumerate(constraint.nd)
                if (node.nodetype == JuMP._Derivatives.VARIABLE)
                    current_value = d.index_to_variable[node.index]
                    if (current_value[1] == current_value[2] == current_value[3] == -1)
                        d.index_to_variable[node.index] = (indx, cindx, 2)
                    end
                end
            end
        end
        for (cindx,subexpress) in enumerate(d.subexpressions)
            for (indx,node) in enumerate(subexpress.nd)
                if (node.nodetype == JuMP._Derivatives.VARIABLE)
                    current_value = d.index_to_variable[node.index]
                    if (current_value[1] == current_value[2] == current_value[3] == -1)
                        d.index_to_variable[node.index] = (indx, cindx, 3)
                    end
                end
            end
        end

        d.subexpression_isnum = fill(true, (d.subexpression_number,))

        return d #deepcopy(d)
    end
end
get_node(d::Evaluator) = d.current_node
