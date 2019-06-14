function copy_to_function(T::S,x::JuMP._FunctionStorage) where {S<:DataType}
    lenx = length(x.nd)
    mc_inf = T(IntervalType(-Inf,Inf))
    temp_set = T[mc_inf for i=1:lenx]
    temp_flt = Array{Float64}(undef,lenx)
    temp_bool = Array{Bool}(undef,lenx)
    FunctionSetStorage{T}(x.nd, x.adj, x.const_values, temp_set, temp_flt, temp_bool,
                          x.grad_sparsity, x.hess_I, x.hess_J, x.dependent_subexpressions)
end

function copy_to_subexpr(T::S,x::JuMP._SubexpressionStorage) where {S<:DataType}
    lenx = length(x.nd)
    mc_inf = T(IntervalType(-Inf,Inf))
    temp_set = T[mc_inf for i=1:lenx]
    temp_flt = Array{Float64}(undef,lenx)
    temp_bool = Array{Bool}(undef,lenx)
    SubexpressionSetStorage{T}(x.nd, x.adj, x.const_values, temp_set, temp_flt,
                               temp_bool, x.linearity)
end

function neg_objective!(d::Evaluator{T}) where T<:Real
    if (d.has_nlobj)
        # shifts the adjacency matrix to introduce -f(x) as first element of nd array
        rowval = rowvals(d.objective.adj) .+ 1; pushfirst!(rowval, 2)
        colptr = copy(d.objective.adj.colptr) .+ 1; pushfirst!(colptr, 1)
        nzval = nonzeros(d.objective.adj); pushfirst!(nzval, true)
        m, n = size(d.objective.adj)
        d.objective.adj = SparseMatrixCSC{Bool,Int}(m+1,n+1,colptr,rowval,nzval)

        # shifts the node list (and parents)
        shift_nd = [JuMP.NodeData(JuMP.CALLUNIVAR,2,2)]
        for nd in d.objective.nd
            push!(shift_nd,JuMP.NodeData(nd.nodetype,nd.index,nd.parent+1))
        end
        d.objective.nd = shift_nd

        nvflag = length(d.objective.numvalued) > 0 ? d.objective.numvalued[1] : false
        pushfirst!(d.objective.numvalued,nvflag)
        pushfirst!(d.objective.numberstorage,0.0)
        pushfirst!(d.objective.setstorage,zero(T))
    end
end

"""
    build_nlp_evaluator

Builds the evaluator used to generate relaxations of the nonlinear equations
and constraints from a source model.
"""
function build_nlp_evaluator(S::R, src::T, x::Optimizer, bool_flag::Bool) where {R<:Type, T<:MOI.AbstractNLPEvaluator}

    # Checks to see if nlp data block evaluator is present
    if ~isa(src,EAGO.EmptyNLPEvaluator)

        # Creates the empty evaluator
        d = Evaluator{S}(src.m)

        num_variables_ = JuMP.num_variables(d.m)
        d.variable_number = num_variables_
        nldata::JuMP._NLPData = deepcopy(d.m.nlp_data)

        # Copy state of user-defined multivariable functions
        d.has_user_mv_operator = src.disable_2ndorder
        d.parameter_values = nldata.nlparamvalues
        d.last_x = fill(NaN, d.variable_number)
        d.last_node = NodeBB()

        # Set valued operators must contain a (sub)gradient and/or (sub)hessian field to support 1st/2nd order eval
        d.disable_1storder = !in(:cv_grad,fieldnames(eltype(d)))
        d.disable_2ndorder = !in(:cv_hess,fieldnames(eltype(d)))

        # Add objective functions, constraints, subexpressions
        d.has_nlobj = isa(nldata.nlobj, JuMP._NonlinearExprData)
        if (src.has_nlobj)
            d.objective = copy_to_function(S,src.objective)
        end

        for i in 1:length(src.constraints)
            push!(d.constraints, copy_to_function(S,src.constraints[i]))
        end

        d.subexpression_order = src.subexpression_order
        d.subexpression_linearity = src.subexpression_linearity
        d.subexpressions_as_julia_expressions = Any[]
        if isdefined(src, :subexpressions_as_julia_expressions)
            d.subexpressions_as_julia_expressions = src.subexpressions_as_julia_expressions
        end

        d.subexpression_values_set = S[]
        d.subexpression_values_flt = Float64[]
        d.subexpressions = SubexpressionSetStorage{S}[]
        for i in 1:length(src.subexpressions)
            temp = copy_to_subexpr(S,src.subexpressions[i])
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
        d.jac_storage = Array{Float64}(undef,max(num_variables_, d.m.nlp_data.largest_user_input_dimension)) # DO I NEED THIS

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

        d.subexpression_isnum = fill(true,(d.subexpression_number,))

        return deepcopy(d)
    end
end

function get_node(d::Evaluator)
    #=
    n_var = d.variable_number
    n = NodeBB(fill(-Inf, (n_var,)), fill(Inf, (n_var,)),
                    d.current_node.lower_bound, d.current_node.upper_bound,
                    d.current_node.depth, d.current_node.last_branch,
                    d.current_node.branch_direction)

    for i in 1:d.variable_number
        indx_num, eqn_num, sto = d.index_to_variable[i]
        if ~(eqn_num == -1)
            if (sto == 1)
                set = d.objective.setstorage[indx_num]
            elseif (sto == 2)
                set = d.constraints[eqn_num].setstorage[indx_num]
            else
                set = d.subexpressions[eqn_num].setstorage[indx_num]
            end
            n.lower_variable_bounds[i] = set.Intv.lo
            n.upper_variable_bounds[i] = set.Intv.hi
        end
    end
    =#
    n = d.current_node
    return n
end
