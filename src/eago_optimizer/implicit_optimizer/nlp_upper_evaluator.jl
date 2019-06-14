const RD_COMPILE_SWITCH = 3000

mutable struct ImplicitUpperEvaluator <: MOI.AbstractNLPEvaluator
    current_node::NodeBB

    disable_1storder::Bool
    disable_2ndorder::Bool
    has_nlobj::Bool

    fg::Function
    func_eval::Bool

    np::Int
    ny::Int
    nx::Int
    ng::Int
    last_y::Vector{Float64}
    diff_result
    diff_tape

    value_storage::Vector{Float64}
    jacobian_storage::VecOrMat{Float64}
    jacobian_sparsity::Vector{Tuple{Int64,Int64}}
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64

    function ImplicitUpperEvaluator()

        d = new()

        d.disable_1storder = false
        d.disable_2ndorder = true
        d.has_nlobj = false
        d.func_eval = false
        d.jacobian_sparsity = Tuple{Int64,Int64}[]

        d.eval_objective_timer = 0.0
        d.eval_constraint_timer = 0.0
        d.eval_objective_gradient_timer = 0.0
        d.eval_constraint_jacobian_timer = 0.0
        d.eval_hessian_lagrangian_timer = 0.0

        return d
    end
end

function implicit_reform_1!(out,y,f,g,h,np,ng,nx)
    out[1] = f(y[1:nx], y[(nx+1):(nx+np)])
    out[2:(ng+1)] = g(y[1:nx], y[(nx+1):(nx+np)])
    h(view(out,(ng+2):(nx+ng+1)), y[1:nx], y[(nx+1):(nx+np)])
    view(out,(2+ng+nx):(1+ng+2*nx))[1:nx] = -view(out,(ng+2):(nx+ng+1))
end

function implicit_reform_2!(out,y,g,h,np,ng,nx)
    out[1:ng] = g(y[1:nx], y[(nx+1):(nx+np)])
    h(view(out,(ng+1):(nx+ng)), y[1:nx], y[(nx+1):(nx+np)])
    view(out,(1+ng+nx):(ng+2*nx))[1:nx]  = -view(out,(ng+1):(nx+ng))
end

function implicit_reform_3!(out,y,f,h,np,nx)
    out[1] = f(y[1:nx], y[(nx+1):(nx+np)])
    h(view(out,2:(nx+1)), y[1:nx], y[(nx+1):(nx+np)])
    view(out,(nx+2):(2*nx+1))[1:nx]  = -view(out,2:(nx+1))
end

function implicit_reform_4!(out,y,h,np,nx)
    h(view(out,1:nx),y[1:nx],y[(nx+1):(nx+np)])
    view(out,(nx+1):(2*nx))[1:nx] = -view(out,1:nx)
end

function set_current_node!(x::ImplicitUpperEvaluator,n::NodeBB)
    x.current_node = n
end

# LOOKS GREAT!
function build_evaluator!(d::ImplicitUpperEvaluator, f, h, np, nx; g = dummy_function,
                          ng::Int = 0, user_sparse = nothing, hj = dummy_function)
    d.nx = nx
    d.ny = np + nx
    d.ng = ng
    d.np = np

    if (d.ng > 0 && f != dummy_function)
        d.fg = (out,x) -> implicit_reform_1!(out,x,f,g,h,np,ng,nx)
        d.has_nlobj = true
    else
        d.fg = (out,x) -> implicit_reform_3!(out,x,f,h,np,nx)
        d.has_nlobj = true
    end

    d.last_y = zeros(d.ny)
    d.value_storage = zeros(1+ng+2*nx)
    d.diff_result = zeros(1+ng+2*nx, nx+np)

    d.diff_tape = ReverseDiff.JacobianTape(d.fg, d.value_storage, d.last_y)
    if length(d.diff_tape) < RD_COMPILE_SWITCH
        d.diff_tape = ReverseDiff.compile(d.diff_tape)
    end

end

# LOOKS GREAT!
function calc_functions!(d::ImplicitUpperEvaluator,y)
    cond = (d.last_y != y)
    if (d.last_y != y)
        if ~d.disable_1storder
            ReverseDiff.jacobian!(d.diff_result, d.diff_tape, y)
        end
        d.fg(d.value_storage,y)
        d.func_eval = true
        d.last_y[:] = y
    end
end

# LOOKS GREAT!
function MOI.eval_objective(d::ImplicitUpperEvaluator, y)
    d.eval_objective_timer += @elapsed begin
        val = 0.0
        if (d.has_nlobj)
            calc_functions!(d,y)
            val = d.value_storage[1]
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

# LOOKS GREAT!
function MOI.eval_constraint(d::ImplicitUpperEvaluator, g, y)
    d.eval_constraint_timer += @elapsed begin
        if (d.ng+d.nx) > 0
            calc_functions!(d,y)
            g[:] = d.value_storage[2:end]
        end
    end
    return
end

# LOOKS GREAT!
function MOI.eval_objective_gradient(d::ImplicitUpperEvaluator, df, y)
    d.eval_objective_timer += @elapsed begin
        if d.has_nlobj
            calc_functions!(d,y)
            df[:] = d.diff_result[1,1:d.ny]
        else
            error("No nonlinear objective.")
        end
    end
    return
end

# LOOKS GREAT!
function MOI.jacobian_structure(d::ImplicitUpperEvaluator)
    # use user-defined sparsity pattern if possible
    if length(d.jacobian_sparsity) > 0
        return d.jacobian_sparsity
    else # else assume dense pattern
        d.jacobian_sparsity = Tuple{Int64,Int64}[(row, idx) for row in 1:(d.ng+2*d.nx) for idx in 1:(d.nx+d.np)]
        return d.jacobian_sparsity
    end
end

# LOOKS GREAT!
function MOI.hessian_lagrangian_structure(d::ImplicitUpperEvaluator)
    error("Hessian computations not currently supported by Implicit optimizer.")
end

# LOOKS GREAT!
function _hessian_lagrangian_structure(d::ImplicitUpperEvaluator)
    error("Hessian lagrangian structure not supported by Implicit optimizer.")
end

# LOOKS GREAT!
function MOI.eval_constraint_jacobian(d::ImplicitUpperEvaluator, dg, y)
    #d.eval_constraint_jacobian_timer += @elapsed begin
        if d.ng+d.nx > 0
            calc_functions!(d,y)
            for (indx,(row,col)) in enumerate(d.jacobian_sparsity)
                dg[indx] = d.diff_result[row+1,col]
            end
        end
    #end
    return
end
#=
# FIX ME
function MOI.eval_constraint_jacobian_product(d::ImplicitUpperEvaluator, out, y, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            if d.ng+d.nx > 0
                calc_functions!(d,y)
                fill!(out,0.0)
                for j in 1:d.ny
                    for i in 1:(d.ng+2*d.nx)
                        y[i] += d.diff_result[i+1,j]*w[j]
                    end
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# FIX ME
function MOI.eval_constraint_jacobian_transpose_product(d::ImplicitUpperEvaluator, y, p, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            if d.ng+d.nx > 0
                calc_functions!(d,y)
                fill!(out,0.0)
                for i in 1:(d.ng+2*d.nx)
                    y[i] += d.diff_result[i+1,:]*w
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end
=#
# LOOKS GREAT
function MOI.features_available(d::ImplicitUpperEvaluator)
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

# LOOKS GREAT
function MOI.initialize(d::ImplicitUpperEvaluator, requested_features::Vector{Symbol}) end
# LOOKS GREAT
MOI.objective_expr(d::ImplicitUpperEvaluator) = error("EAGO.ImplicitUpperEvaluator doesn't provide expression graphs of constraint functions.")
#LOOKS GREAT
MOI.constraint_expr(d::ImplicitUpperEvaluator) = error("EAGO.ImplicitUpperEvaluator doesn't provide expression graphs of constraint functions.")
