"""
    update_lower_variable_bounds!

Updates bounds on optimizer `z` by reset MOI set values for each constraint.
"""
function update_lower_variable_bounds!(x::Optimizer,y::NodeBB,z::T) where {T<:MOI.AbstractOptimizer}
    # Updates variables bounds
    for i=1:x.variable_number
        var = x.variable_info[i]

        if (~var.is_integer)
            ci1,ci2,num = x.lower_variable_index[i] # x.VariableIndexLow[i]
            if var.is_fixed
                MOI.set(z, MOI.ConstraintSet(), ci1, MOI.EqualTo{Float64}(y.lower_variable_bounds[i]))
            elseif var.has_lower_bound
                if var.has_upper_bound
                    MOI.set(z, MOI.ConstraintSet(), ci1, MOI.LessThan{Float64}(y.upper_variable_bounds[i]))
                    MOI.set(z, MOI.ConstraintSet(), ci2, MOI.GreaterThan{Float64}(y.lower_variable_bounds[i]))
                else
                    MOI.set(z, MOI.ConstraintSet(), ci1, MOI.GreaterThan{Float64}(y.lower_variable_bounds[i]))
                end
            elseif var.has_upper_bound
                MOI.set(z, MOI.ConstraintSet(), ci1, MOI.LessThan{Float64}(y.upper_variable_bounds[i]))
            end
        else
            #=
            MOI.set(z, MOI.ConstraintSet(), x.VariableIndex[i], MOI.Interval(x.CurrentLowerInfo.Solution[i],
                                                                             x.CurrentLowerInfo.Solution[i]))
            if var.is_fixed
                MOI.set(z, MOI.ConstraintSet(), MOI.VariableIndex(x.VariableNumber), MOI.EqualTo(var.upper_bound))
            elseif var.has_lower_bound
                if var.has_upper_bound
                    MOI.set(z, MOI.ConstraintSet(), MOI.VariableIndex(x.VariableNumber), MOI.LessThan(var.upper_bound))
                    MOI.set(z, MOI.ConstraintSet(), MOI.VariableIndex(x.VariableNumber), MOI.GreaterThan(var.lower_bound))
                else
                    MOI.set(z, MOI.ConstraintSet(), MOI.VariableIndex(x.VariableNumber), MOI.GreaterThan(var.lower_bound))
                end
            elseif var.has_upper_bound
                MOI.set(z, MOI.ConstraintSet(), MOI.VariableIndex(x.VariableNumber), MOI.LessThan(var.upper_bound)
            end
            =#
        end
    end
end

"""
    update_lower_variable_bounds1!

Updates bounds on optimizer `z` by adding new box valued constraints.
"""
function update_lower_variable_bounds1!(m::Optimizer,y::NodeBB,z::T) where {T<:MOI.AbstractOptimizer}
    # Updates variables bounds
    typevar = m.lower_variables
    for (i,var) in enumerate(m.variable_info)
        var_xi = MOI.SingleVariable(typevar[i])
        if var.is_integer
        else
            if var.is_fixed
                MOI.add_constraint(z, var_xi, MOI.EqualTo(y.lower_variable_bounds[i]))
            elseif var.has_lower_bound
                if var.has_upper_bound
                    MOI.add_constraint(z, var_xi, MOI.LessThan(y.upper_variable_bounds[i]))
                    MOI.add_constraint(z, var_xi, MOI.GreaterThan(y.lower_variable_bounds[i]))
                else
                    MOI.add_constraint(z, var_xi, MOI.GreaterThan(y.lower_variable_bounds[i]))
                end
            elseif var.has_upper_bound
                MOI.add_constraint(z, var_xi, MOI.LessThan(y.upper_variable_bounds[i]))
            end
        end
    end
end


"""
    update_upper_variable_bounds!

Updates bounds on optimizer `z` by adding new box valued constraints.
"""
function update_upper_variable_bounds!(m::Optimizer, y::NodeBB, z::T) where {T<:MOI.AbstractOptimizer}
    # Updates variables bounds
    midxi = 0.0
    typevar = m.upper_variables
    for (i,var) in enumerate(m.variable_info)
        var_ind = typevar[i]
        var_xi = MOI.SingleVariable(typevar[i])
        if var.is_integer
        else
            if var.is_fixed
                MOI.add_constraint(z, var_xi, MOI.EqualTo(y.lower_variable_bounds[i]))
            elseif var.has_lower_bound
                if var.has_upper_bound
                    MOI.add_constraint(z, var_xi, MOI.LessThan(y.upper_variable_bounds[i]))
                    MOI.add_constraint(z, var_xi, MOI.GreaterThan(y.lower_variable_bounds[i]))
                else
                    MOI.add_constraint(z, var_xi, MOI.GreaterThan(y.lower_variable_bounds[i]))
                end
            elseif var.has_upper_bound
                MOI.add_constraint(z, var_xi, MOI.LessThan(y.upper_variable_bounds[i]))
            end
            midxi = (y.upper_variable_bounds[i]+y.lower_variable_bounds[i])/2.0
            MOI.set(z, MOI.VariablePrimalStart(), var_ind, midxi)
        end
        #push!(m.VariableIndexUpp,VarTupleUpp)
    end
end

"""
    minus_objective!

Takes `d::JuMP.NLPEvaluator` which evaluates to the objective f(x) and transforms
it to evaluate -f(x)
"""
function minus_objective!(d::JuMP.NLPEvaluator,want_hess_storage)
    if (d.has_nlobj)
        # shifts the adjacency matrix to introduce -f(x) as first element of nd array
        rowval = rowvals(d.objective.adj) .+ 1; pushfirst!(rowval, 2)
        colptr = Vector{Int}(d.objective.adj.colptr) .+ 1; pushfirst!(colptr, 1)
        nzval = nonzeros(d.objective.adj); pushfirst!(nzval, true)
        m, n = size(d.objective.adj)
        d.objective.adj = SparseMatrixCSC{Bool,Int}(m+1,n+1,colptr,rowval,nzval)

        # shifts the node list (and parents)
        shift_nd = [JuMP.NodeData(JuMP.CALLUNIVAR,2,2)]
        for nd in d.objective.nd
            push!(shift_nd,JuMP.NodeData(nd.nodetype,nd.index,nd.parent+1))
        end
        d.objective.nd = shift_nd
        pushfirst!(d.objective.forward_storage,0.0)
        pushfirst!(d.objective.partials_storage,0.0)
        pushfirst!(d.objective.reverse_storage,0.0)

        # fix object seed matrix TODO: ADD CODE
        # d.objective.seed_matrix =

        # recompute max_chuck, expr_length, num_variables_
        max_expr_length = length(d.objective.nd)
        max_chunk = size(d.objective.seed_matrix,2)
        for k in 1:length(d.constraints)
            max_expr_length = max(max_expr_length, length(d.constraints[end].nd))
            max_chunk = max(max_chunk, size(d.constraints[end].seed_matrix,2))
        end
        max_chunk = min(max_chunk, 10)
        num_variables_ = JuMP.num_variables(d.m)

        if d.want_hess || want_hess_storage # storage for Hess or HessVec
        #if length(d.input_ϵ) > 0  # CURRENT ERROR (SOMETHING IS undefined here)
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
    end
end

"""
    set_local_nlp!

Copies constraints, objective, and variable into a structure reference by the
local nlp solver `initial_upper_optimizer`.
"""
function set_local_nlp!(m::Optimizer)

    # Add linear and quadratic constraints to model
    for (func, set) in m.linear_leq_constraints
         MOI.add_constraint(m.initial_upper_optimizer,func,set)
    end
    for (func, set) in m.linear_geq_constraints
        MOI.add_constraint(m.initial_upper_optimizer,func,set)
    end
    for (func, set) in m.linear_eq_constraints
        MOI.add_constraint(m.initial_upper_optimizer,func,set)
    end
    for (func,set,ind) in m.linear_interval_constraints
        MOI.add_constraint(m.initial_upper_optimizer, func, MOI.GreaterThan{Float64}(set.lower))
        MOI.add_constraint(m.initial_upper_optimizer, func, MOI.LessThan{Float64}(set.upper))
    end

    for (func, set) in m.quadratic_leq_constraints
        MOI.add_constraint(m.initial_upper_optimizer,func,set)
    end
    for (func, set) in m.quadratic_geq_constraints
        MOI.add_constraint(m.initial_upper_optimizer,func,set)
    end
    for (func, set) in m.quadratic_eq_constraints
        MOI.add_constraint(m.initial_upper_optimizer,func,set)
    end
    for (func,set,ind) in m.quadratic_interval_constraints
        MOI.add_constraint(m.initial_upper_optimizer, func, MOI.GreaterThan{Float64}(set.lower))
        MOI.add_constraint(m.initial_upper_optimizer, func, MOI.LessThan{Float64}(set.upper))
    end

    # Add nonlinear evaluation block
    MOI.set(m.initial_upper_optimizer, MOI.NLPBlock(), m.nlp_data)

    # Add objective sense
    #MOI.set(m.InitialUpperOptimizer, MOI.ObjectiveSense(), m.OptimizationSense)
    #MOI.set(m.initial_upper_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    # Specifies variables in upper problem
    m.upper_variables = MOI.VariableIndex.(1:m.variable_number)

    objmult = (m.optimization_sense == MOI.MIN_SENSE) ? 1.0 : -1.0

    # Add objective function (if any)
    if (m.objective != nothing)
        MOI.set(m.initial_upper_optimizer, MOI.ObjectiveSense(), m.optimization_sense)
        if isa(m.objective, MOI.SingleVariable)
            if (m.optimization_sense == MOI.MIN_SENSE)
                MOI.set(m.initial_upper_optimizer, MOI.ObjectiveFunction{MOI.SingleVariable}(), m.objective)
            elseif (m.optimization_sense == MOI.MAX_SENSE)
                neg_obj_var = MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm{Float64}(-1.0, m.objective.variable)], 0.0)
                MOI.set(m.initial_upper_optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), neg_obj_var)
            else
                error("Objective sense must be MOI.MIN_SENSE or MOI.MAX_SENSE")
            end
        elseif isa(m.objective, MOI.ScalarAffineFunction{Float64})
            if (m.optimization_sense == MOI.MIN_SENSE)
                MOI.set(m.initial_upper_optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), m.objective)
            elseif (m.optimization_sense == MOI.MAX_SENSE)
                neg_obj_aff_terms = []
                for term in m.objective.terms
                    push!(neg_obj_aff_terms,MOI.ScalarAffineTerm{Float64}(-term.coefficient,term.variable_index))
                end
                neg_obj_aff = MOI.ScalarAffineFunction{Float64}(neg_obj_aff_terms, -m.objective.constant)
                MOI.set(m.initial_upper_optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), neg_obj_aff)
            else
                error("Objective sense must be MOI.MinSense or MOI.MaxSense")
            end
        elseif isa(m.objective, MOI.ScalarQuadraticFunction{Float64})
            if (m.optimization_sense == MOI.MIN_SENSE) || (m.optimization_sense == MOI.MAX_SENSE)
                MOI.set(m.initial_upper_optimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), m.objective)
            #=
            elseif (m.optimization_sense == MOI.MAX_SENSE)
                neg_obj_qda_terms = []
                for term in m.objective.affine_terms
                    push!(neg_obj_aff_terms,MOI.ScalarAffineTerm{Float64}(-term.coefficient,term.variable_index))
                end
                neg_obj_qdq_terms = []
                for term in m.objective.quadratic_terms
                    push!(neg_obj_aff_terms,MOI.ScalarQuadraticTerm{Float64}(-term.coefficient,term.variable_index_1,term.variable_index_2))
                end
                neg_obj_qd = ScalarQuadraticFunction{Float64}(neg_obj_qda_terms,neg_obj_qdq_terms,-m.objective.constant)
                MOI.set(m.initial_upper_optimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), neg_obj_qd)
                =#
            else
                error("Objective sense must be MOI.MIN_SENSE or MOI.MAX_SENSE")
            end
        end
    else
        @assert m.nlp_data != empty_nlp_data()
        MOI.set(m.initial_upper_optimizer, MOI.ObjectiveSense(), m.optimization_sense)
        #=
        if (m.optimization_sense == MOI.MAX_SENSE)
            #println("ran minus objective")
            #minus_objective!(m.nlp_data.evaluator, false)
        elseif (m.optimization_sense != MOI.MIN_SENSE)
            #error("Objective sense must be MOI.MIN_SENSE or MOI.MAX_SENSE")
        end
        =#
    end
end

"""
    set_to_default!

Checks to see if the user set any of EAGO's subroutines. If not, the subroutines
are set to default values.
"""
function set_to_default!(m::Optimizer)
    (m.lower_problem! == dummy_function)       &&   (m.lower_problem! = default_lower_bounding!)
    (m.upper_problem! == dummy_function)       &&   (m.upper_problem! = default_upper_bounding!)
    (m.preprocess! == dummy_function)         &&   (m.preprocess! = default_preprocess!)
    (m.postprocess! == dummy_function)        &&   (m.postprocess! = default_postprocess!)
    (m.single_check == dummy_function)        &&   (m.single_check = default_repeat_check)
    (m.convergence_check == dummy_function)   &&   (m.convergence_check = default_convergence_check)
    (m.termination_check == dummy_function)   &&   (m.termination_check = default_termination_check)
    (m.node_storage! == dummy_function)        &&   (m.node_storage! = default_storage!)
    (m.node_selection == dummy_function)      &&   (m.node_selection = node_select_best!)
    (m.bisection_function == dummy_function)  &&   (m.bisection_function = continuous_relative_bisect)
    (m.cut_condition == dummy_function)       &&   (m.cut_condition = default_cut_condition)
    (m.add_cut! == dummy_function)             &&   (m.add_cut! = default_add_cut!)
    (m.relax_function! == dummy_function)      &&   (m.relax_function! = relax_model!)
    isa(m.initial_relaxed_optimizer, DummyOptimizer) && (m.initial_relaxed_optimizer = GLPK.Optimizer())
    isa(m.working_relaxed_optimizer, DummyOptimizer) && (m.working_relaxed_optimizer = GLPK.Optimizer())
    #isa(m.initial_relaxed_optimizer, DummyOptimizer) && (m.initial_relaxed_optimizer = CPLEX.Optimizer())
    #isa(m.working_relaxed_optimizer, DummyOptimizer) && (m.working_relaxed_optimizer = CPLEX.Optimizer())
    isa(m.initial_upper_optimizer, DummyOptimizer) && (m.initial_upper_optimizer = Ipopt.Optimizer())
    if isa(m.lower_factory.constructor, DummyOptimizer)
        m.lower_factory = JuMP.with_optimizer(GLPK.Optimizer; m.lower_optimizer_options...)
        #m.lower_factory = JuMP.with_optimizer(Clp.Optimizer)
        m.use_lower_factory = true
    end
    if isa(m.upper_factory.constructor, DummyOptimizer)
        m.upper_factory = JuMP.with_optimizer(Ipopt.Optimizer; m.upper_optimizer_options...)
        m.use_upper_factory = true
    end
    isa(m.working_upper_optimizer, DummyOptimizer) && (m.working_upper_optimizer = Ipopt.Optimizer())
    isa(m.linear_optimizer, DummyOptimizer) && (m.linear_optimizer = GLPK.Optimizer())
    #isa(m.linear_optimizer, DummyOptimizer) && (m.linear_optimizer = CPLEX.Optimizer())
    isa(m.nlp_optimizer, DummyOptimizer) && (m.nlp_optimizer = Ipopt.Optimizer())
end

function push_variable_bounds!(var::VariableInfo,var_xi,m)
    if var.is_integer
    else
        if var.is_fixed
            ci1 = MOI.add_constraint(m, var_xi, MOI.EqualTo(var.upper_bound))
            return ci1,ci1,1
        elseif var.has_lower_bound
            if var.has_upper_bound
                ci1 = MOI.add_constraint(m, var_xi, MOI.LessThan(var.upper_bound))
                ci2 = MOI.add_constraint(m, var_xi, MOI.GreaterThan(var.lower_bound))
                return ci1,ci2,2
            else
                ci1 = MOI.add_constraint(m, var_xi, MOI.GreaterThan(var.lower_bound))
                return ci1,ci1,1
            end
        elseif var.has_upper_bound
            ci1 = MOI.add_constraint(m, var_xi, MOI.LessThan(var.upper_bound))
            return ci1,ci1,1
        end
    end
end

function push_lower_variables!(m::Optimizer)
    # Copies the same variables to every submodel
    MOI.add_variables(m.initial_relaxed_optimizer, m.variable_number)
    x = MOI.add_variables(m.working_relaxed_optimizer, m.variable_number)
    for (i,var) in enumerate(m.variable_info)
        var_xi = MOI.SingleVariable(x[i])
        push_variable_bounds!(var, var_xi, m.initial_relaxed_optimizer)
        VarTupleLow = push_variable_bounds!(var, var_xi, m.working_relaxed_optimizer)
        push!(m.lower_variable_index,VarTupleLow)
    end
end

"""
    is_globally_optimal

Takes an `MOI.TerminationStatusCode` and a `MOI.ResultStatusCode` and returns
the tuple `(valid_result::Bool, feasible::Bool)`. The value `valid_result` is
`true` if the pair of codes prove that either the subproblem solution was solved
to global optimality or the subproblem solution is infeasible. The value of
`feasible` is true if the problem is feasible and false if the problem is infeasible.
"""
function is_globally_optimal(t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)

    feasible = false; valid_result = false

    if (t == MOI.INFEASIBLE && r == MOI.INFEASIBILITY_CERTIFICATE)
        valid_result = true
    elseif (t == MOI.INFEASIBLE && r == MOI.NO_SOLUTION)
        valid_result = true
    elseif (t == MOI.INFEASIBLE && r == MOI.UNKNOWN_RESULT_STATUS)
        valid_result = true
    elseif (t == MOI.OPTIMAL && r == MOI.FEASIBLE_POINT)
        valid_result = true
        feasible = true
    elseif (t == MOI.INFEASIBLE_OR_UNBOUNDED && r == MOI.NO_SOLUTION)
        valid_result = true
        feasible = false
    end

    return valid_result, feasible
end

"""
    is_feasible_solution

Takes an `MOI.TerminationStatusCode` and a `MOI.ResultStatusCode` and returns `true`
if this corresponds to a solution that is proven feasible. Returns `false` otherwise.
"""
function is_feasible_solution(t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)

    termination_flag = false; result_flag = false

    (t == MOI.OPTIMAL) && (termination_flag = true)
    (t == MOI.LOCALLY_SOLVED) && (termination_flag = true)

    (r == MOI.FEASIBLE_POINT) && (result_flag = true)

    return (termination_flag && result_flag)
end

"""
    set_dual!

Retrieves the lower and upper duals for variable bounds from the
`working_relaxed_optimizer` and sets the appropriate values in the
`current_lower_info` field.
"""
function set_dual!(x::Optimizer)
    for (vi,VarIndxTuple) in enumerate(x.lower_variable_index)
        (ci1,ci2,n) = VarIndxTuple
        if (n == 2)
            if isa(ci1,MOI.ConstraintIndex{MOI.SingleVariable,MOI.GreaterThan{Float64}})
                x.current_lower_info.lower_variable_dual[vi] = MOI.get(x.working_relaxed_optimizer, MOI.ConstraintDual(), ci1)
                x.current_lower_info.upper_variable_dual[vi] = MOI.get(x.working_relaxed_optimizer, MOI.ConstraintDual(), ci2)
            else
                x.current_lower_info.lower_variable_dual[vi] = MOI.get(x.working_relaxed_optimizer, MOI.ConstraintDual(), ci2)
                x.current_lower_info.upper_variable_dual[vi] = MOI.get(x.working_relaxed_optimizer, MOI.ConstraintDual(), ci1)
            end
        else
            if isa(ci1,MOI.ConstraintIndex{MOI.SingleVariable,MOI.GreaterThan{Float64}})
                x.current_lower_info.lower_variable_dual[vi] = MOI.get(x.working_relaxed_optimizer, MOI.ConstraintDual(), ci1)
            elseif isa(ci1,MOI.ConstraintIndex{MOI.SingleVariable,MOI.LessThan{Float64}})
                x.current_lower_info.upper_variable_dual[vi] = MOI.get(x.working_relaxed_optimizer, MOI.ConstraintDual(), ci1)
            else
                error("Not a variable constraint index.")
            end
        end
    end
end
