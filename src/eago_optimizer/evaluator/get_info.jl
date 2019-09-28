function MOI.eval_objective(d::Evaluator, x::Vector{Float64})
    forward_reverse_pass(d,x)
    val = 0.0
    if d.has_nlobj
        if d.objective.numvalued[1]
            val = d.objective.numberstorage[1]
        else
            val = d.objective.setstorage[1].cv
        end
    else
        error("No nonlinear objective.")
    end
    return val
end

function eval_objective_lo(d::Evaluator)
    val = 0.0
    if d.has_nlobj
        if d.objective.numvalued[1]
            val = d.objective.numberstorage[1]
        else
            val = d.objective.setstorage[1].cv
        end
    else
        error("No nonlinear objective.")
    end
    return val
end

get_node_lower(d::FunctionSetStorage, i::Int64) = d.setstorage[i].Intv.lo
get_node_upper(d::FunctionSetStorage, i::Int64) = d.setstorage[i].Intv.hi

function MOI.eval_constraint(d::Evaluator, g::Vector{Float64}, x::Vector{Float64})
    forward_reverse_pass(d,x)
    for i in 1:length(d.constraints)
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].cv
        end
    end
    return
end

function eval_constraint_cc(d::Evaluator, g::Vector{Float64}, y::Vector{Float64})
    forward_reverse_pass(d,y)
    for i in 1:length(d.constraints)
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].cc
        end
    end
    return
end

function eval_constraint_lo!(d::Evaluator, g::Vector{Float64})
    for i in 1:length(d.constraints)
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].lo
        end
    end
    return
end

function eval_constraint_hi!(d::Evaluator, g::Vector{Float64})
    for i in 1:length(d.constraints)
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].hi
        end
    end
    return
end

function MOI.eval_objective_gradient(d::Evaluator, df::Vector{Float64}, x::Vector{Float64})
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
    else
        error("First order information unavailable.")
    end
    return
end

# TO DO
function MOI.eval_constraint_jacobian_transpose_product(d::Evaluator, y, x, w)
    if (!d.disable_1storder)
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
    else
        error("First order information unavailable.")
    end
    return
end

# looks good
function MOI.features_available(d::Evaluator)
    features = Symbol[]
    if !d.disable_1storder
        push!(features,:Grad)
        push!(features,:Jac)
    end
    return features
end

# looks good, doesn't do anything, EAGO builds the evaluator and attaches it to lower problems
function MOI.initialize(d::Evaluator, requested_features::Vector{Symbol})
end

function grad_sparsity(d::Evaluator, j::Int64)
    sparsity = Int64[]
    if j == 1
        sparsity = d.objective.grad_sparsity
    else
        sparsity = d.constraints[j-1].grad_sparsity
    end
    return sparsity
end
