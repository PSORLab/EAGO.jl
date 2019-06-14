const DEBUG_IMPLICIT_LOWER = false

mutable struct ImplicitLowerEvaluator{N} <: MOI.AbstractNLPEvaluator
    current_node::NodeBB
    last_node::NodeBB

    disable_1storder::Bool
    disable_2ndorder::Bool
    has_nlobj::Bool
    has_ineq::Bool
    objective_fun
    constraints_fun
    state_fun
    state_jac_fun

    np::Int
    ng::Int
    nx::Int

    objective_ubd::Float64
    constraints_lbd::Vector{Float64}
    constraints_ubd::Vector{Float64}

    jacobian_sparsity::Vector{Tuple{Int64,Int64}}
    obj_eval::Bool
    cnstr_eval::Bool
    init_relax_run::Bool

    P::Vector{IntervalType}
    X::Vector{IntervalType}

    H::Vector{MC{N}}
    J::Array{MC{N},2}
    Y::Array{Float64,2}

    ref_pMC::Vector{MC{N}}
    xa_mc::Vector{MC{N}}
    xA_mc::Vector{MC{N}}
    z_mc::Vector{MC{N}}
    aff_mc::Vector{MC{N}}

    last_p::Vector{Float64}
    ref_p::Vector{Float64}
    var_relax::Vector{MC{N}}

    state_relax_1::Vector{MC{N}}
    state_ref_relaxation_1::Array{MC{N},3}

    state_relax_n::Array{MC{N},2}
    state_ref_relaxation_n::Array{MC{N},3}

    obj_relax::MC{N}
    cnstr_relax::Vector{MC{N}}
    imp_opts::mc_opts

    # timer
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64

    function ImplicitLowerEvaluator{N}() where N
        d = new()

        d.obj_eval = false
        d.cnstr_eval = false
        d.init_relax_run = false
        d.jacobian_sparsity = Tuple{Int64,Int64}[]

        d.np = N
        d.nx = 0
        d.ng = 0

        d.objective_ubd = Inf
        d.constraints_lbd = Float64[]
        d.constraints_ubd = Float64[]

        d.eval_objective_timer = 0.0
        d.eval_constraint_timer = 0.0
        d.eval_objective_gradient_timer = 0.0
        d.eval_constraint_jacobian_timer = 0.0
        d.eval_hessian_lagrangian_timer = 0.0

        d.var_relax = fill(zero(MC{N}),(N,))

        d.state_relax_1 = [zero(MC{N})]
        d.state_ref_relaxation_1 = fill(zero(MC{N}), (2,2,2,))
        d.state_relax_n = fill(zero(MC{N}), (2,2,))
        d.state_ref_relaxation_n = fill(zero(MC{N}), (2,2,2,))

        d.obj_relax = zero(MC{N})
        d.cnstr_relax = [zero(MC{N})]
        d.P = fill(IntervalType(0.0),(1,))
        d.X = fill(IntervalType(0.0),(1,))

        d.H = fill(zero(MC{N}), (1,))
        d.J = fill(zero(MC{N}), (2,2))
        d.Y = fill(0.0, (1,1))

        d.ref_pMC = fill(zero(MC{N}), (1,))
        d.xa_mc = fill(zero(MC{N}), (1,))
        d.xA_mc = fill(zero(MC{N}), (1,))
        d.z_mc = fill(zero(MC{N}), (1,))
        d.aff_mc = fill(zero(MC{N}), (1,))

        d.last_p = Float64[0.0]
        d.ref_p = Float64[0.0]

        d.disable_1storder = false
        d.disable_2ndorder = true
        d.has_nlobj = false

        d.current_node = NodeBB()
        d.last_node = NodeBB()

        d.imp_opts = mc_opts()

        return d
    end
end

num_state_variables(eval::ImplicitLowerEvaluator) = eval.nx
num_decision_variables(eval::ImplicitLowerEvaluator) = eval.np
"""
    build_evaluator!

A function that builds the ImplicitLowerEvaluator form user-supplied functions.
"""
function build_evaluator!(d::ImplicitLowerEvaluator, impfun::Function, np::Int, nx::Int;
                                obj = nothing, constr = nothing, ng::Int = 0, user_sparse = nothing,
                                state_jac = nothing)

    # setup objective and constraint functions
    d.has_nlobj = obj != nothing
    d.has_ineq = constr != nothing
    d.objective_fun = obj
    d.constraints_fun = constr
    d.state_fun = impfun
    yimp = (out, y) -> d.state_fun(out, y[1:nx], y[(nx+1):(nx+np)])

    if (state_jac == nothing)
        d.state_jac_fun = (out, x, p) -> ForwardDiff.jacobian(yimp, out, vcat(x, p))[:,1:nx]
    else
        d.state_jac_fun = state_jac
    end

    # has a nonlinear objective?
    (obj == nothing) || (d.has_nlobj = true)

    # get dimension sizes
    d.np = np; d.nx = nx;
    if (constr == nothing)
        d.ng = 0
    else
        if (ng > 0 )
            d.ng = ng
        else
            d.ng = length(constr(ones(nx),ones(np)))
        end
    end

    # set implicit routine information
    d.imp_opts.nx = nx
    d.imp_opts.np = np

    # preallocates the storage variables
    temp = zero(MC{np})
    d.cnstr_relax = fill(temp,(ng,))
    if (nx == 1)
        d.state_relax_1 = fill(temp, (1,))
        d.state_ref_relaxation_1 = fill(temp, (1, d.imp_opts.kmax, 1))
    else
        d.state_relax_n = fill(temp, (nx, 1))
        d.state_ref_relaxation_n = fill(temp, (nx, d.imp_opts.kmax, 1))
    end

    d.P = fill(IntervalType(0.0),(np,))
    d.X = fill(IntervalType(0.0),(nx,))

    d.H = fill(zero(MC{np}), (nx,))
    d.J = fill(zero(MC{np}), (nx,nx))
    d.Y = fill(0.0, (nx,nx))

    d.ref_pMC = fill(zero(MC{np}), (np,))
    d.xa_mc = fill(zero(MC{np}), (nx,))
    d.xA_mc = fill(zero(MC{np}), (nx,))
    d.z_mc = fill(zero(MC{np}), (nx,))
    d.aff_mc = fill(zero(MC{np}), (nx,))

    # allocates the reference points
    d.last_p = zeros(Float64,np); fill!(d.last_p,NaN)
    d.ref_p = zeros(Float64,np)

    if (user_sparse == nothing)
        sparse_pattern = Tuple{Int64,Int64}[]
        for i in 1:nx
            for j in 1:np
                push!(sparse_pattern,(i,j))
            end
        end
        d.jacobian_sparsity = sparse_pattern
    else
        d.jacobian_sparsity = user_sparse
    end
end

function set_current_node!(x::ImplicitLowerEvaluator,n::NodeBB)
    x.current_node = n
end

function set_last_node!(x::ImplicitLowerEvaluator,n::NodeBB)
    x.current_node = n
end

# Is p_ref value or mc? implicit_relax_h! ?
function relax_implicit!(d::ImplicitLowerEvaluator,y)
    nx = d.nx
    np = d.np
    is_same_box = same_box(d.current_node, d.last_node, 0.0)
    # Generate new parameters for implicit relaxation if necessary
    if ~is_same_box
        d.init_relax_run = true
        d.obj_eval = false
        d.cnstr_eval = false
        d.last_node = d.current_node
        for i in 1:nx
            d.X[i] = IntervalType(d.current_node.lower_variable_bounds[i], d.current_node.upper_variable_bounds[i])
        end
        for j in 1:np
            shiftj = j + nx
            d.ref_p[j] = (d.current_node.lower_variable_bounds[shiftj] + d.current_node.upper_variable_bounds[shiftj])/2.0
            d.P[j] = IntervalType(d.current_node.lower_variable_bounds[shiftj], d.current_node.upper_variable_bounds[shiftj])
            d.ref_pMC[j] = MC{np}(d.ref_p[j], d.P[j], j)
            d.var_relax[j] = MC{np}(d.ref_p[j], d.P[j], j)
        end
        if (nx > 0)
            if (nx == 1)
                x0 = d.state_relax_1
                state_relax_1 = d.state_relax_1
                println("x0: $x0")
                println("state_relax: $state_relax_1")
                # view(d.state_relax_1, i:i)             # state_relax_1::Vector{MC{N}}
                # view(d.state_relax_n, 1:nx, 1:1)       # state_relax_n::Array{MC{N},2}
                gen_expansion_params!(d.state_fun, d.state_jac_fun, d.ref_pMC,
                                      x0, state_relax_1, d.xa_mc, d.xA_mc, d.z_mc,
                                      d.aff_mc, d.X, d.P, d.imp_opts,
                                      d.state_ref_relaxation_1,
                                      d.H, d.J, d.Y, true, Float64[], true)
            else
                x0 = d.state_relax_n
                state_relax_n = d.state_relax_n
                println("x0: $x0")
                println("state_relax: $state_relax_n")
                gen_expansion_params!(d.state_fun, d.state_jac_fun, d.ref_pMC,
                                      x0, state_relax_n, d.xa_mc, d.xA_mc, d.z_mc,
                                      d.aff_mc, d.X, d.P, d.imp_opts,
                                      d.state_ref_relaxation_n,
                                      d.H, d.J, d.Y, true, Float64[], true)
            end
        end
    end
    # Generate new value of implicit relaxation
    new_point_flag = false
    for i in 1:nx
        ref_vs_y = (d.ref_p[i] != y[i+np])
        if ref_vs_y
            new_point_flag = true; break
        end
    end

    #if new_point_flag && (~d.init_relax_run)
        d.obj_eval = false
        d.cnstr_eval = false
        pMC = fill(zero(MC{np}), (np,))
        for j in 1:np
            shiftj = j + nx
            sg_pi = seed_gradient(Float64, j, np)
            pMC[j] = MC{np}(y[shiftj], y[shiftj], d.P[j], sg_pi, sg_pi, false)
            d.var_relax[j] = pMC[j]
        end
        if (nx > 0)
            if (nx == 1)
                implicit_relax_h!(d.state_fun, d.state_jac_fun, pMC, d.ref_pMC,
                                  d.state_relax_1, d.state_relax_1,
                                  d.xa_mc, d.xA_mc, d.z_mc, d.aff_mc, d.X,
                                  d.P, d.imp_opts, d.state_ref_relaxation_1,
                                  d.H, d.J, d.Y, true, Float64[], true)
            else
                implicit_relax_h!(d.state_fun, d.state_jac_fun, pMC, d.ref_pMC,
                                  d.state_relax_n, d.state_relax_n,
                                  d.xa_mc, d.xA_mc, d.z_mc, d.aff_mc, d.X,
                                  d.P, d.imp_opts, d.state_ref_relaxation_n,
                                  d.H, d.J, d.Y, true, Float64[], true)
            end
        end
#    else
#        d.state_relax[:] = d.state_ref_relaxation[end][:]
#    end
end

# LOOKS GREAT!
function relax_objective!(d::ImplicitLowerEvaluator)
    if ~d.obj_eval
        if d.nx == 1
            d.obj_relax = d.objective_fun(d.state_relax_1, d.var_relax)
        else
            d.obj_relax = d.objective_fun(d.state_relax_n, d.var_relax)
        end
        d.obj_eval = true
        DEBUG_IMPLICIT_LOWER && println("mc objective: $(d.obj_relax)")
    end
end

# LOOKS GREAT!
function relax_constraints!(d::ImplicitLowerEvaluator)
    if ~d.cnstr_eval
        if d.nx == 1
            eval_val  = d.constraints_fun(d.state_relax_1, d.var_relax)
            d.cnstr_relax[:] = eval_val
        else
            eval_val  = d.constraints_fun(d.state_relax_n, d.var_relax)
            d.cnstr_relax[:] = eval_val
        end
        d.cnstr_eval = true
        DEBUG_IMPLICIT_LOWER && println("mc constraint: $(d.cnstr_relax)")
    end
end

# LOOKS GREAT!
function MOI.eval_objective(d::ImplicitLowerEvaluator, y)
    d.eval_objective_timer += @elapsed begin
        val = zero(eltype(y))
        if d.has_nlobj
            relax_implicit!(d,y)
            relax_objective!(d)
            val = d.obj_relax.cv
            DEBUG_IMPLICIT_LOWER && println("float objective: $(val)")
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

# LOOKS GREAT!
function MOI.eval_constraint(d::ImplicitLowerEvaluator, g, y)
    d.eval_constraint_timer += @elapsed begin
        if d.ng > 0
            relax_implicit!(d,y)
            relax_constraints!(d)
            for i in 1:d.ng
                g[i] = d.cnstr_relax[i].cv
            end
            DEBUG_IMPLICIT_LOWER && println("float constraints: $(g)")
        end
    end
    return
end

# LOOKS GREAT!
function MOI.eval_objective_gradient(d::ImplicitLowerEvaluator, df, y)
    d.eval_objective_timer += @elapsed begin
        if d.has_nlobj
            relax_implicit!(d,y)
            relax_objective!(d)
            for j in 1:d.np
                df[j] = d.obj_relax.cv_grad[j]
            end
            DEBUG_IMPLICIT_LOWER && println("float objective gradient: $(df)")
        else
            error("No nonlinear objective.")
        end
    end
    return
end

# LOOKS GREAT!
function MOI.jacobian_structure(d::ImplicitLowerEvaluator)
    # use user-defined sparsity pattern if possible
    if length(d.jacobian_sparsity) > 0
        return d.jacobian_sparsity
    else # else assume dense pattern
        d.jacobian_sparsity = Tuple{Int64,Int64}[(row, idx) for row in 1:d.ng for idx in 1:d.np]
        return d.jacobian_sparsity
    end
end

# LOOKS GREAT!
function MOI.hessian_lagrangian_structure(d::ImplicitLowerEvaluator)
    error("Hessian computations not currently supported by Implicit optimizer.")
end

# LOOKS GREAT!
function _hessian_lagrangian_structure(d::ImplicitLowerEvaluator)
    error("Hessian lagrangian structure not supported by Implicit optimizer.")
end

# LOOKS GREAT!
function MOI.eval_constraint_jacobian(d::ImplicitLowerEvaluator,g,y)
    d.eval_constraint_jacobian_timer += @elapsed begin
        if d.ng > 0
            relax_implicit!(d,y)
            relax_constraints!(d)
            fill!(g, 0.0)
            if (d.ng == 1) && (d.np == 1)
                g[1] = d.cnstr_relax[1].cv_grad[1]
            else
                for (i,j) in d.jacobian_sparsity
                    temp = d.cnstr_relax[i].cv_grad[j]
                    g[i,j] = d.cnstr_relax[i].cv_grad[j]
                end
            end
            DEBUG_IMPLICIT_LOWER && println("float constraint jacobian: $(g)")
        end
    end
    return g
end

# LOOKS GREAT!
function MOI.eval_constraint_jacobian_product(d::ImplicitLowerEvaluator, out, y, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            if d.ng > 0
                relax_implicit!(d,y)
                relax_constraints!(d)
                fill!(out, 0.0)
                for (i,j) in d.jacobian_sparsity
                    out[i] += d.cnstr_relax[i].cv_grad[j]*w[j]
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# PROBABLY DONE
function MOI.eval_constraint_jacobian_transpose_product(d::ImplicitLowerEvaluator, out, y, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            if d.ng > 0
                relax_implicit!(d,y)
                relax_constraints!(d)
                fill!(out, 0.0)
                for (i,j) in d.jacobian_sparsity
                    out[i] += d.cnstr_relax[i].cv_grad[j]*w[j]
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# LOOKS GREAT
function MOI.features_available(d::ImplicitLowerEvaluator)
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
function MOI.initialize(d::ImplicitLowerEvaluator, requested_features::Vector{Symbol}) end

# LOOKS GREAT
MOI.objective_expr(d::ImplicitLowerEvaluator) = error("EAGO.ImplicitLowerEvaluator doesn't provide expression graphs of constraint functions.")
#LOOKS GREAT
MOI.constraint_expr(d::ImplicitLowerEvaluator) = error("EAGO.ImplicitLowerEvaluator doesn't provide expression graphs of constraint functions.")

function eval_constraint_cc(d::ImplicitLowerEvaluator, g, y)
    d.eval_constraint_timer += @elapsed begin
        if d.ng > 0
            relax_implicit!(d,y)
            relax_constraints!(d)
            for i in 1:d.ng
                g[i] = d.cnstr_relax[i].cc
            end
        end
    end
    return
end

function eval_constraint_cc_grad(d::ImplicitLowerEvaluator, g, y)
    d.eval_constraint_jacobian_timer += @elapsed begin
        if d.ng > 0
            relax_implicit!(d,y)
            relax_constraints!(d)
            fill!(g, 0.0)
            if (d.ng == 1) && (d.np == 1)
                g[1] = d.cnstr_relax[1].cc_grad[1]
            else
                for (i,j) in d.jacobian_sparsity
                    temp = d.cnstr_relax[i].cc_grad[j]
                    g[i,j] = d.cnstr_relax[i].cc_grad[j]
                end
            end
        end
    end
    return g
end
