function set_to_mid!(x::Vector{Float64},n::NodeBB)
    x[:] = (n.lower_variable_bounds + n.upper_variable_bounds)/2.0
end

function interval_preprocess(x::Optimizer,y::NodeBB)

    Eflag = false
    Iflag = false
    eDflag = false

    set_current_node!(x.nlp_data.evaluator,y)
    nx = x.state_variables
    np = x.variable_number - x.state_variables
    current_node = x.nlp_data.evaluator.current_node

    x.nlp_data.evaluator.x_lower[:] = lower_variable_bounds(current_node, 1, nx)
    x.nlp_data.evaluator.x_upper[:] = upper_variable_bounds(current_node, 1, nx)
    x.nlp_data.evaluator.p_lower[:] = lower_variable_bounds(current_node, nx+1, nx+np)
    x.nlp_data.evaluator.p_upper[:] = upper_variable_bounds(current_node, nx+1, nx+np)
    x.nlp_data.evaluator.interval_X[:] = IntervalType.(x.nlp_data.evaluator.x_lower, x.nlp_data.evaluator.x_upper)
    x.nlp_data.evaluator.interval_P[:] = IntervalType.(x.nlp_data.evaluator.p_lower, x.nlp_data.evaluator.p_upper)

    Eflag, Iflag, eDflag = param_intv_contractor(x.nlp_data.evaluator.interval_opts.h,
                                                 x.nlp_data.evaluator.interval_opts.hj,
                                                 x.nlp_data.evaluator.interval_X,
                                                 x.nlp_data.evaluator.Ntemp,
                                                 x.nlp_data.evaluator.N,
                                                 x.nlp_data.evaluator.Xi,
                                                 x.nlp_data.evaluator.X1,
                                                 EAGO.IntervalType[], Float64[],
                                                 x.nlp_data.evaluator.Y,
                                                 x.nlp_data.evaluator.J,
                                                 x.nlp_data.evaluator.H,
                                                 x.nlp_data.evaluator.interval_P,
                                                 x.nlp_data.evaluator.inc,
                                                 x.nlp_data.evaluator.incLow,
                                                 x.nlp_data.evaluator.incHigh,
                                                 x.nlp_data.evaluator.nx,
                                                 x.nlp_data.evaluator.kmax,
                                                 x.nlp_data.evaluator.etol,
                                                 x.nlp_data.evaluator.rtol)

    x.current_preprocess_info.feasibility = ~Eflag
    if ~Eflag
        y.lower_variable_bounds[1:nx] = lo.(x.nlp_data.evaluator.X1)
        y.upper_variable_bounds[1:nx] = hi.(x.nlp_data.evaluator.X1)
    end
end

function midpoint_upper_bnd(x::Optimizer,y::NodeBB)
    if is_integer_feasible(x) #&& mod(x.CurrentIterationCount,x.UpperBoundingInterval) == 1

        set_current_node!(x.nlp_data.evaluator,y)
        set_to_mid!(x.current_upper_info.solution,y)
        if x.nlp_data.evaluator.ng > 0
            g = zeros(x.nlp_data.evaluator.ng)
            MOI.eval_constraint(x.nlp_data.evaluator, g, x.current_upper_info.solution)
            result_status = any(i-> (i > 0), g) ? MOI.INFEASIBLE_POINT : MOI.FEASIBLE_POINT
        else
            result_status = MOI.FEASIBLE_POINT
        end
        if (result_status == MOI.FEASIBLE_POINT)
            x.current_upper_info.feasibility = true
            val = MOI.eval_objective(x.nlp_data.evaluator,x.current_upper_info.solution)
            x.current_upper_info.value = val
        else
            x.current_upper_info.feasibility = false
            x.current_upper_info.value = Inf
        end
    else
        x.current_upper_info.feasibility = false
        x.current_upper_info.value = Inf
    end
end

# Modifies functions post initial relaxation to use appropriate nlp evaluator
function implicit_mod!(opt::Optimizer,args)

    ImpLowerEval = args[1]; ImpUpperEval = args[2]; has_obj = args[3];
    lower_bnds = args[4]; upper_bnds = args[5]; alt_upper_flag = args[6];
    alt_upper = args[7]; nx = args[8]

    opt.preprocess! = interval_preprocess
    opt.relax_function! = implicit_relax_model!

    # load lower nlp data block
    lower_eval_block = MOI.NLPBlockData(lower_bnds, ImpLowerEval, has_obj)
    opt.working_evaluator_block = deepcopy(lower_eval_block)
    if MOI.supports(opt.initial_relaxed_optimizer, MOI.NLPBlock())
        opt.initial_relaxed_optimizer.nlp_data = deepcopy(lower_eval_block)
        opt.working_relaxed_optimizer.nlp_data = deepcopy(lower_eval_block)
    end

    # load upper nlp data block
    upper_eval_block = MOI.NLPBlockData(upper_bnds, ImpUpperEval, has_obj)
    # if using the midpoint evaluator don't setup upper optimizers &
    if alt_upper_flag
        opt.upper_problem! = alt_upper
    else
        if MOI.supports(opt.initial_relaxed_optimizer, MOI.NLPBlock())
            opt.initial_upper_optimizer.nlp_data = deepcopy(upper_eval_block)
            opt.working_upper_optimizer.nlp_data = deepcopy(upper_eval_block)
        end
    end

    opt.nlp_data = upper_eval_block
    opt.nlp_data.evaluator.param_interval_opts = parametric_interval_params(:Dense,:Newton,1E-12,1E-12,nx,nx,8)
end

"""
    solve_implicit



Inputs:
* `f::Function`: Objective in the decision variable. Takes a single argument
                 vector that must be untyped.
* `h!::Function`: The semi-infinite constraint. Takes two arguments: the first
                    being a vector containing the decision variable and the
                    second being a vector containing the uncertainity
                    variables. The function must be untyped.
* `x_l::Vector{Float64}`: Lower bounds on the state variables.
* `x_u::Vector{Float64}`: Upper bounds on the state variables.
* `p_l::Vector{Float64}`: Lower bounds on the decision variables.
* `p_u::Vector{Float64}`: Upper bounds on the decision variables.
* `m::JuMP.Model`: JuMP Model
* `hj!::Function`: Equality constraints to be relaxed via an implicit approach
* `g!::Function`: Jacobian of equality constraints to be relaxed

Returns: Variable vector and models.
"""
function solve_implicit(f, h!, xl, xu, pl, pu, m, hj!, g!; upper_sym = :MidPointUpperEvaluator,
                        lower = ImplicitLowerEvaluator{1}())

    opt = backend(m).optimizer.model.optimizer
    h_xold!(Hout, X, Xold, P, t) = h!(Hout, X, P)
    hj_xold!(Jout, X, Xold, P, t) = hj!(Jout, X, P)

    if upper_sym == :MidPointUpperEvaluator
        upper = MidPointUpperEvaluator()
    elseif upper_sym == :ImplicitUpperEvaluator
        upper = ImplicitUpperEvaluator()
    else
        error("unsupported evaluator type")
    end

    # get dimensions
    @assert length(pl) == length(pu)
    @assert length(xl) == length(xu)
    np = length(pl); nx = length(xl);

    if (g! == nothing)
        ng = 0
        g_xold! = g!
    else
        ng = length(g_func(ones(nx),ones(np)))
        g_xold!(Gout, X, Xold, P, t) = g!(Gout, X, P)
    end

    # sets most routines to default (expect bisection)
    set_to_default!(opt)
    opt.bisection_function = implicit_bisection

    # add variables to lower, upper, and EAGO models
    var_EAGO = MOI.add_variables(opt, np+nx)

    for j in 1:nx
        MOI.add_constraint(opt, var_EAGO[j], MOI.GreaterThan(xl[j]))
        MOI.add_constraint(opt, var_EAGO[j], MOI.LessThan(xu[j]))
    end

    for i in 1:np
        MOI.add_constraint(opt, var_EAGO[i+nx], MOI.GreaterThan(pl[i]))
        MOI.add_constraint(opt, var_EAGO[i+nx], MOI.LessThan(pu[i]))
    end

    # Build the lower implicit evaluator
    if isa(lower, ImplicitLowerEvaluator{1})
        lower = ImplicitLowerEvaluator{np}()
    end
    build_evaluator!(lower, h_xold!, np, nx, obj = f, constr = g_xold!, state_jac = hj_xold!, ng = ng)
                           #, user_sparse = user_sparsity)

    # Build the upper evaluator
    if g! == nothing
        g_xold! = nothing
    end
    build_evaluator!(upper, f, h_xold!, np, nx, g = g_xold!, ng = ng, hj = hj_xold!) #, user_sparse = user_sparsity)

    # Specify number of state variables
    opt.state_variables = nx
    opt.treat_as_nonlinear = [i for i in (nx+1):(nx+np)]
    opt.upper_has_node = true

    # Add nlp data blocks ("SHOULD" BE THE LAST THING TO DO)
    has_obj = (f != nothing)
    bnd_pair = MOI.NLPBoundsPair(-Inf,0.0)
    lower_bnds = [bnd_pair for i=1:ng]
    upper_bnds = [bnd_pair for i=1:(ng+2*nx)]

    # Set the objective sense
    MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    # Optimizes the model with load function
    alt_upper = ~isa(upper,ImplicitUpperEvaluator)
    custom_mod_args = (lower, upper, has_obj, lower_bnds, upper_bnds, alt_upper, midpoint_upper_bnd, nx)
    MOI.optimize!(opt, custom_mod! = implicit_mod!, custom_mod_args = custom_mod_args)

    return var_EAGO,opt
end
