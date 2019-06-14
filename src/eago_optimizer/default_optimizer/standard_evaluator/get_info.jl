function MOI.eval_objective(d::Evaluator, x)
    d.eval_objective_timer += @elapsed begin
        forward_reverse_pass(d,x)
        val = zero(eltype(x))
        if d.has_nlobj
            if d.objective.numvalued[1]
                val = d.objective.numberstorage[1]
            else
                val = d.objective.setstorage[1].cv
            end
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

get_node_lower(d::FunctionSetStorage, i::Int) = d.setstorage[i].Intv.lo
get_node_upper(d::FunctionSetStorage, i::Int) = d.setstorage[i].Intv.hi

function MOI.eval_constraint(d::Evaluator, g, x)
    d.eval_constraint_timer += @elapsed begin
        forward_reverse_pass(d,x)
        for i in 1:length(d.constraints)
            if d.constraints[i].numvalued[1]
                g[i] = d.constraints[i].numberstorage[1]
            else
                g[i] = d.constraints[i].setstorage[1].cv
            end
        end
    end
    return
end

function MOI.eval_objective_gradient(d::Evaluator, df, x)
    d.eval_objective_timer += @elapsed begin
        forward_reverse_pass(d,x)
        if d.has_nlobj
            if ~d.objective.numvalued[1]
                for j in 1:length(d.objective.setstorage[1].cv_grad)
                    df[j] = d.objective.setstorage[1].cv_grad[j]
                end
            end
        else
            error("No nonlinar objective.")
        end
    end
    return
end

function MOI.jacobian_structure(d::Evaluator)
    jacobian_sparsity = Tuple{Int64,Int64}[]
    for row in 1:length(d.constraints)
        row_sparsity = d.constraints[row].grad_sparsity
        for idx in row_sparsity
            push!(jacobian_sparsity, (row, idx))
        end
    end
    return jacobian_sparsity
end

function MOI.hessian_lagrangian_structure(d::Evaluator)
    !d.disable_2ndorder || error("Hessian computations were not requested on the call to initialize!.")
    return d.hessian_sparsity
end

function _hessian_lagrangian_structure(d::Evaluator)
    hessian_sparsity = Tuple{Int64,Int64}[]
    if d.has_nlobj
        for idx in 1:length(d.objective.hess_I)
            push!(hessian_sparsity, (d.objective.hess_I[idx], d.objective.hess_J[idx]))
        end
    end
    for ex in d.constraints
        for idx in 1:length(ex.hess_I)
            push!(hessian_sparsity, (ex.hess_I[idx], ex.hess_J[idx]))
        end
    end
    return hessian_sparsity
end

function MOI.eval_constraint_jacobian(d::Evaluator,g,x)
    #d.eval_constraint_jacobian_timer += @elapsed begin
        forward_reverse_pass(d,x)
        #t = typeof(d.constraints[1].setstorage[1])
        #g = zero.(g)
        for i in 1:length(d.constraints)
            if ~d.constraints[i].numvalued[1]
                for j in 1:d.variable_number
                    g[i,j] = d.constraints[i].setstorage[1].cv_grad[j]
                end
            else
                for j in 1:length(d.constraints[i].setstorage[1].cv_grad)
                    g[i,j] = 0.0
                end
            end
        end
    #end
    return
end

function eval_constraint_cc_grad(d::Evaluator, g, y)
    #d.eval_constraint_jacobian_timer += @elapsed begin
        forward_reverse_pass(d,y)
        #t = typeof(d.constraints[1].setstorage[1])
        #g = zero.(g)
        for i in 1:length(d.constraints)
            if ~d.constraints[i].numvalued[1]
                for j in 1:d.variable_number
                    g[i,j] = d.constraints[i].setstorage[1].cc_grad[j]
                end
            else
                for j in 1:d.variable_number
                    g[i,j] = 0.0
                end
            end
        end
    #end
    return
end

# TO DO (CHECK GRADIENT DIMS)
function MOI.eval_constraint_jacobian_product(d::Evaluator, y, x, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            forward_reverse_pass(d,x)
            t = typeof(d.constraints[1].setstorage[1])
            y = zeros(t,length(d.constraints[1].setstorage[1].cv_grad),length(d.constraints))
            for i in 1:length(d.constraints)
                if ~d.constraints[i].numvalued[1]
                    for j in 1:d.variable_number
                        y[i] += d.constraints[i].setstorage[1].cv_grad[j]*w[j]
                    end
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# TO DO
function MOI.eval_constraint_jacobian_transpose_product(d::Evaluator, y, x, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            forward_reverse_pass(d,x)
            #t = typeof(d.constraints[1].setstorage[1])
            y = zeros(Float64,length(d.constraints[1].setstorage[1].cv_grad),length(d.constraints))
            for i in 1:length(d.constraints)
                if ~d.constraints[i].numvalued[1]
                    for j in 1:d.variable_number
                        y[i] += d.constraints[i].setstorage[1].cv_grad[j]*w[j]
                    end
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# TO DO
function MOI.eval_hessian_lagrangian_product(d::Evaluator, h, x, v, σ, μ)

    return h
end

# TO DO
function MOI.eval_hessian_lagrangian(d::Evaluator, H, x, σ, μ)
    if (!d.disable_2ndorder)
        d.eval_hessian_lagrangian_timer += @elapsed begin
            forward_reverse_pass(d,x)
            t = typeof(d.constraints[i].setstorage[1])
            H[:,:] =  σ*d.objective.setstorage[1].cv_hess
            for i in 1:length(d.constraints)
                if ~d.constraints[i].numvalued[1]
                    H[:,:] += μ[i]*d.constraints[i].setstorage[1].cv_hess
                end
            end
        end
    end
end

# looks good
function MOI.features_available(d::Evaluator)
    features = Symbol[]
    if !d.disable_1storder
        push!(features,:Grad)
        push!(features,:Jac)
    end
    if !d.disable_2ndorder
        push!(features,:Hess)
        push!(features,:HessVec)
    end
    return features
end

# looks good, doesn't do anything, EAGO builds the evaluator and attaches it to lower problems
function MOI.initialize(d::Evaluator, requested_features::Vector{Symbol})
end

# looks good
MOI.objective_expr(d::Evaluator) = error("EAGO.Evaluator doesn't provide expression graphs of constraint functions.")
# looks good
MOI.constraint_expr(d::Evaluator) = error("EAGO.Evaluator doesn't provide expression graphs of constraint functions.")

function eval_constraint_cc(d::Evaluator, g, y)
    d.eval_constraint_timer += @elapsed begin
        forward_reverse_pass(d,y)
        for i in 1:length(d.constraints)
            if d.constraints[i].numvalued[1]
                g[i] = d.constraints[i].numberstorage[1]
            else
                g[i] = d.constraints[i].setstorage[1].cc
            end
        end
    end
    return
end
