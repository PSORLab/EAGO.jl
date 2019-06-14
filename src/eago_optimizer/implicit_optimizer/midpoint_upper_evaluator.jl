"""
    MidPointUpperEvaluator

A mutable struct that holds all information required when calculating upper
bounds via a interval-newton method taken at the midpoint of the decisions space
P.
"""
mutable struct MidPointUpperEvaluator <: MOI.AbstractNLPEvaluator
    current_node::NodeBB

    disable_1storder::Bool
    disable_2ndorder::Bool
    has_nlobj::Bool

    fg::Function
    h::Function
    hj::Function
    func_eval::Bool

    interval_opts::interval_newton_opt
    param_interval_opts::parametric_interval_params

    np::Int
    ny::Int
    nx::Int
    ng::Int
    last_box::Vector{IntervalType}
    interval_P::Vector{IntervalType}
    interval_X::Vector{IntervalType}

    H::Vector{IntervalType}
    J::Array{IntervalType,2}
    Y::Array{Float64,2}
    B::Vector{IntervalType}
    M::Array{IntervalType,2}

    x_lower::Vector{Float64}
    x_upper::Vector{Float64}
    p_lower::Vector{Float64}
    p_upper::Vector{Float64}

    kmax::Int
    etol::Float64
    rtol::Float64

    N::Vector{EAGO.IntervalType}
    Ntemp::Vector{EAGO.IntervalType}
    Xi::Vector{EAGO.IntervalType}
    X1::Vector{EAGO.IntervalType}

    inc::Vector{Bool}
    incLow::Vector{Bool}
    incHigh::Vector{Bool}

    value_storage::Vector{IntervalType}
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64

    m::JuMP.Model

    function MidPointUpperEvaluator()

        d = new()

        d.disable_1storder = false
        d.disable_2ndorder = true
        d.has_nlobj = false
        d.func_eval = false
        d.fg = x->x

        d.np = 0
        d.ny = 0
        d.nx = 0
        d.ng = 0

        d.interval_opts = interval_newton_opt()
        d.param_interval_opts = parametric_interval_params(0,0)

        d.x_lower = Float64[]
        d.x_upper = Float64[]
        d.p_lower = Float64[]
        d.p_upper = Float64[]

        d.last_box = IntervalType[]
        d.interval_P = IntervalType[]
        d.interval_X = IntervalType[]
        d.value_storage = Float64[]

        d.H = IntervalType[]
        d.J = fill(zero(IntervalType), (2,2))
        d.Y = fill(0.0, (2,2))
        d.B = IntervalType[]
        d.M = fill(zero(IntervalType), (2,2))

        d.eval_objective_timer = 0.0
        d.eval_constraint_timer = 0.0
        d.eval_objective_gradient_timer = 0.0
        d.eval_constraint_jacobian_timer = 0.0
        d.eval_hessian_lagrangian_timer = 0.0

        d.m = JuMP.Model()

        return d
    end
end

function implicit_reform_5!(out,x,p,f,g,ng::Int)
    out[1] = f(x,p)
    out[2:(ng+1)] = g(x,p)
end

function implicit_reform_6!(out,x,p,f)
    out[1] = f(x,p)
end

function set_current_node!(x::MidPointUpperEvaluator,n::NodeBB)
    x.current_node = n
end

function set_intervals!(x::MidPointUpperEvaluator)
    x.x_lower[:] = lower_variable_bounds(x.current_node,1,x.nx)
    x.x_upper[:] = upper_variable_bounds(x.current_node,1,x.nx)
    x.p_lower[:] = lower_variable_bounds(x.current_node,x.nx+1,x.ny)
    x.p_upper[:] = upper_variable_bounds(x.current_node,x.nx+1,x.ny)
    x.interval_X[:] = IntervalType.(x.x_lower,x.x_upper)
    x.interval_P[:] = IntervalType.((x.p_lower+x.p_upper)/2.0)
end

"""
    build_midpointupperevaluator!

A function that builds the MidPointUpperEvaluator form user-supplied functions.
"""
function build_evaluator!(d::MidPointUpperEvaluator, f, h, np, nx; hj = dummy_function,
                          g = nothing, ng = 0, user_sparse = nothing)

    if (g != nothing && ng == 0)
        calc_ng = length(g(zeros(Float64,nx),zeros(Float64,np)))
    else
        calc_ng = ng
    end

    d.nx = nx
    d.ny = np + nx
    d.ng = calc_ng
    d.np = np

    d.param_interval_opts = parametric_interval_params(nx,0)
    d.kmax = d.param_interval_opts.kmax
    d.etol = d.param_interval_opts.etol
    d.rtol = d.param_interval_opts.rtol

    d.N = fill(zero(IntervalType), (nx,))
    d.Ntemp = fill(zero(IntervalType), (nx,))
    d.Xi = fill(zero(IntervalType), (nx,))
    d.X1 = fill(zero(IntervalType), (nx,))

    d.inc = fill(false, (nx,))
    d.incLow = fill(false, (nx,))
    d.incHigh = fill(false, (nx,))

    d.H = fill(zero(IntervalType), (nx,))
    d.J = fill(zero(IntervalType), (nx,nx))
    d.Y = fill(0.0, (nx,nx))
    d.B = fill(zero(IntervalType), (nx,))
    d.M = fill(zero(IntervalType), (nx,nx))

    if (d.ng > 0 && f != dummy_function)
        d.fg = (out,x,p) -> implicit_reform_5!(out,x,p,f,g,d.ng)
        d.has_nlobj = true
    else
        d.fg = (out,x,p) -> implicit_reform_6!(out,x,p,f)
        d.has_nlobj = true
    end

    d.param_interval_opts= parametric_interval_params(nx,nx)
    d.x_lower = zeros(Float64,nx)
    d.x_upper = zeros(Float64,nx)
    d.p_lower = zeros(Float64,np)
    d.p_upper = zeros(Float64,np)

    if hj == dummy_function
        yimp = (out, y) -> h(out, y[1:nx], y[(nx+1):(nx+np)])
        hj_ad = (out, x, p) -> ForwardDiff.jacobian(yimp, out, vcat(x, p))[:,1:nx]
        d.interval_opts = interval_newton_opt(d.nx,h,hj_ad)
    else
        d.interval_opts = interval_newton_opt(d.nx,h,hj)
    end

    d.interval_X = zeros(IntervalType,nx)
    d.interval_P = zeros(IntervalType,np)

    if (f != dummy_function)
        d.value_storage = zeros(1+d.ng)
    else
        d.value_storage = zeros(d.ng)
    end
end

"""
    calc_functions!(d::MidPointUpperEvaluator,y)

Main calculation kernel used for midpoint upper bounding calculation.
"""
function calc_functions!(d::MidPointUpperEvaluator,y)
    set_intervals!(d)
    interval_newton_nb!(d.H,d.J,d.Y,d.B,d.M,d.interval_X,d.interval_P,d.interval_opts)
    d.fg(d.value_storage,d.interval_X,d.interval_P)
    d.func_eval = true
end

function MOI.eval_objective(d::MidPointUpperEvaluator, y)
    d.eval_objective_timer += @elapsed begin
        val = 0.0
        if (d.has_nlobj)
            calc_functions!(d,y)
            val = d.value_storage[1].hi
            #println("interval val: $val")
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

function MOI.eval_constraint(d::MidPointUpperEvaluator, g, y)
    d.eval_constraint_timer += @elapsed begin
        if (d.ng) > 0
            calc_functions!(d,y)
            g[:] = lo.(d.value_storage[2:end])
            #println("interval g: $g")
        end
    end
    return
end

MOI.features_available(d::MidPointUpperEvaluator) = Symbol[]
function MOI.initialize(d::MidPointUpperEvaluator, requested_features::Vector{Symbol}) end
