#=
"""
    dag_cse_simplify!

"""
function dag_cse_simplify!(d::NLPEvaluator) where T<:Real
end
=#

#=
function post_initialize_modification!(d::NLPEvaluator)

    nldata::_NLPData = d.m.nlp_data
    num_variables_ = num_variables(d.m) + 1

    moi_index_to_consecutive_index[last_indx] = last_indx
    #Dict(moi_index => consecutive_index for (consecutive_index, moi_index) in enumerate(MOI.get(d.m, MOI.ListOfVariableIndices())))

    d.user_output_buffer = Array{Float64}(undef,d.m.nlp_data.largest_user_input_dimension)
    d.jac_storage = Array{Float64}(undef,max(num_variables_, d.m.nlp_data.largest_user_input_dimension))

    d.last_x = fill(NaN, num_variables_)

end
=#

#=
"""
    reform_epigraph!

Reforms the optimization problem in EAGO. Modifies the optimizer and evaluators but doesn't
modify the associated JuMP model included in either
"""
function reform_epigraph!(m::Optimizer)
    println("start reform epigraph")
    orig_num_variables = num_variables(m.nlp_data.evaluator.m)
    num_variables_ = orig_num_variables + 1
    if (m.nlp_data.has_objective || isa(m.objective, MOI.ScalarQuadraticFunction))

        # calculate objective bounds
        d =  m.working_evaluator_block.evaluator
        x = ones(d.variable_number)
        d.current_node = NodeBB()
        d.current_node.lower_variable_bounds = lower_bound.(m.variable_info)
        d.current_node.upper_variable_bounds = upper_bound.(m.variable_info)
        forward_eval_obj(d,x)
        lower = get_node_lower(d.objective,1)
        upper = get_node_upper(d.objective,1)
        gt = MOI.GreaterThan{Float64}(lower)
        lt = MOI.LessThan{Float64}(upper)

        # Add objective variable and bounds
        vi = MOI.add_variable(m)
        var = MOI.SingleVariable(vi)
        cli = MOI.add_constraint(m, var, gt)
        cui = MOI.add_constraint(m, var, lt)

        # Update root node bounds & problem storage
        push!(m.stack[1].lower_variable_bounds, lower)
        push!(m.stack[1].upper_variable_bounds, upper)
        push!(m.current_lower_info.solution, 0.0)
        push!(m.current_lower_info.lower_variable_dual, 0.0)
        push!(m.current_lower_info.upper_variable_dual, 0.0)
        push!(m.current_upper_info.solution, 0.0)

        # Add replace nlp evaluator (NLconstraint mod) or replace quadratic constraint
        if m.nlp_data.has_objective
            if (m.optimization_sense == MOI.MIN_SENSE)
                lower = -Inf
                upper = 0.0
                op_index = 2
            elseif (m.optimization_sense == MOI.MAX_SENSE)
                lower = 0.0
                upper = Inf
                op_index = 1
            end

            # increase variable number
            orig_num_variables = num_variables(m.nlp_data.evaluator.m)
            num_variables_ = orig_num_variables + 1
            # assumes that only evaluator needs to be updated...
            # not attempting to modify the remaining model object

            # Build and add the constraint
            new_constraint_bounds = m.nlp_data.constraint_bounds
            push!(new_constraint_bounds, MOI.NLPBoundsPair(lower, upper))

            new_evaluator = m.nlp_data.evaluator
            t1 = typeof(new_evaluator.m)
            println("t1: $t1")

            new_constraint = new_evaluator.objective
            new_constraint.nd[1] = JuMP._Derivatives.NodeData(new_constraint.nd[1].nodetype, new_constraint.nd[1].index, 1)
            for i in 2:length(new_constraint.nd)
                new_constraint.nd[i] = JuMP._Derivatives.NodeData(new_constraint.nd[i].nodetype,
                                                                  new_constraint.nd[i].index,
                                                                  new_constraint.nd[i].parent+2)
            end
            pushfirst!(new_constraint.nd, JuMP._Derivatives.NodeData(JuMP._Derivatives.CALL, op_index, -1))
            insert!(new_constraint.nd, 3, JuMP._Derivatives.NodeData(JuMP._Derivatives.VARIABLE, num_variables_, 1))

            new_constraint.adj = JuMP.adjmat(new_constraint.nd)
            push!(new_constraint.forward_storage, 0.0)
            push!(new_constraint.partials_storage, 0.0)
            push!(new_constraint.grad_sparsity, num_variables_)

            new_evaluator.want_hess = false
            new_constraint.hess_I = Int[]
            new_constraint.hess_J = Int[]
            new_constraint.rinfo = JuMP.Coloring.RecoveryInfo()

            push!(new_evaluator.constraints, new_constraint)

            # Remove the objective
            # new_evaluator.objective = nothing       # don't bother rewritting objective yet
            new_evaluator.has_nlobj = false

            # Finish evaluator update
            push!(new_evaluator.last_x, 0.0)
            new_evaluator.jac_storage = Array{Float64}(undef,max(num_variables_, d.m.nlp_data.largest_user_input_dimension))

            # UNDEFINED ERORR IN BELOW BLOCK
            #chunk_calc = length(new_evaluator.input_ϵ)/orig_num_variables
            #for i in 1:chunk_calc
            #    push!(new_evaluator.input_ϵ, undef)
            #    push!(new_evaluator.output_ϵ, undef)
            #end
            #
            if new_evaluator.want_hess
                new_evaluator.hessian_sparsity = JuMP._hessian_lagrangian_structure(new_evaluator)
            end

            # Replace block
            new_block = MOI.NLPBlockData(new_constraint_bounds, new_evaluator, false)
            m.nlp_data = new_block
        elseif isa(m.objective, MOI.ScalarQuadraticFunction)
            affine_terms = m.objective.affine_terms
            quadratic_terms = m.objective.quadratic_terms
            constantv = m.objective.constant
            if (m.optimization_sense == MOI.MIN_SENSE)
                push!(affine_terms, MOI.ScalarAffineTerm{Float64}(-1.0,vi))
                set = MOI.LessThan{Float64}(0.0)
            else
                push!(affine_terms, MOI.ScalarAffineTerm{Float64}(1.0,vi))
                set = MOI.GreaterThan{Float64}(0.0)
            end
            func = MOI.ScalarQuadraticFunction{Float64}(affine_terms, quadratic_terms, constantv)
            MOI.add_constraint(m, func, set)
        end

        # Set Objective function to single variable
        MOI.set(m, MOI.ObjectiveFunction{MOI.SingleVariable}(), var)
    end
    println("end reform epigraph1")

    # Sets up relaxations terms that don't vary during iterations (mainly linear)
    push!(m.lower_variables, MOI.VariableIndex.(m.variable_number))
    println("end reform epigraph2")

    #=
    # Build the JuMP NLP (upper) evaluator
    evaluator = m.nlp_data.evaluator
    features = MOI.features_available(evaluator)
    has_hessian = (:Hess in features)
    init_feat = [:Grad]
    #has_hessian && push!(init_feat, :Hess)
    num_nlp_constraints = length(m.nlp_data.constraint_bounds)
    num_nlp_constraints > 0 && push!(init_feat, :Jac)
    MOI.initialize(evaluator,init_feat)
    post_initialize_modification!(evaluator)
    push!(m.nlp_data.evaluator.last_x, 0.0)
    println("end reform epigraph3")

    println("m.nlp_data.evaluator.last_x: $(m.nlp_data.evaluator.last_x)")
    num_variables_ = num_variables(m.nlp_data.evaluator.m)
    type_inner_m = typeof(m.nlp_data.evaluator.m)
    println("type_inner_m: $(type_inner_m)")
    type_inner_m_back = typeof(backend(m.nlp_data.evaluator.m))
    println("type_inner_m_back: $(type_inner_m_back)")
    println("num_variables_: $(num_variables_)")
    =#

    # Rebuild lower nlp evaluator
    m.working_evaluator_block = m.nlp_data
    if ~isa(m.nlp_data.evaluator, EAGO.EmptyNLPEvaluator)
        built_evaluator = build_nlp_evaluator(MC{m.variable_number}, m.nlp_data.evaluator, m, m.reform_epigraph_flag)
        #(m.optimization_sense == MOI.MAX_SENSE) && neg_objective!(built_evaluator)
        m.working_evaluator_block = MOI.NLPBlockData(m.nlp_data.constraint_bounds, built_evaluator, m.nlp_data.has_objective)
    end
    println("end reform epigraph4")
end
=#

function reform_epigraph!(m::Optimizer)

    # calculate objective bounds
    d =  m.working_evaluator_block.evaluator
    x = ones(d.variable_number)
    d.current_node = NodeBB()
    d.current_node.lower_variable_bounds = lower_bound.(m.variable_info)
    d.current_node.upper_variable_bounds = upper_bound.(m.variable_info)
    forward_eval_obj(d,x)
    lower = get_node_lower(d.objective,1)
    upper = get_node_upper(d.objective,1)

    # get model
    inner_model = m.nlp_data.evaluator.m

    # add variable
    @variable(inner_model, lower <= nu <= upper)

    # get model info
    nvar = JuMP.num_variables(inner_model)
    sense = JuMP.objective_sense(inner_model)
    expr = inner_model.nlp_data.nlobj

    # add objective
    @objective(inner_model, sense, nu)

    # what is largest user input dimensions???
    # add i
    if sense == MOI.MIN_SENSE
        pushfirst!(expr.nd, JuMP.NodeData(JuMP._Derivatives.CALL, 2, -1),
                            JuMP.NodeData(JuMP._Derivatives.VARIABLE, nvar+1, 1))
        for i in 3:length(expr.nd)
            expr.nd[i] = JuMP.NodeData(expr.nd[i].nodetype, expr.nd[i].index, expr.nd[i].parent+1)
        end
    elseif sense == MOI.MAX_SENSE
        pushfirst!(expr.nd, JuMP.NodeData(JuMP._Derivatives.CALL, 1, -1),
                            JuMP.NodeData(JuMP._Derivatives.VARIABLE, nvar+1, 1))
        for i in 3:length(expr.nd)
            expr.nd[i] = JuMP.NodeData(expr.nd[i].nodetype, expr.nd[i].index, expr.nd[i].parent+1)
        end
    end
    c =  JuMP._NonlinearConstraint(expr, 0.0, Inf)
    #inner_model.nlp_data.constraint_bounds
    push!(inner_model.nlp_data.nlconstr, c)
end
