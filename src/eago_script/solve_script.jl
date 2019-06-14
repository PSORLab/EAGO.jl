
# Modifies functions post initial relaxation to use appropriate nlp evaluator
function script_mod(opt::Optimizer,args)
    # loads relaxation, builds relaxed evaluator if necessary
    Eval = args[1]; has_obj = args[2]; c_bnds = args[3]; nx = args[4];
    ng = args[5]; nlobj = args[6]; nlconstr = args[7]; nlexpr = args[8];
    user_operators = args[9]

    println("Eval: $Eval")

    lower_eval_block = MOI.NLPBlockData(c_bnds, Eval, has_obj)
    opt.relax_function! = EAGO.relax_model!
    opt.optimization_sense = MOI.MIN_SENSE
    Eval.has_nlobj = true
    for i in 1:nx
        opt.nonlinear_variable[i] = true
    end

    nlconstr_duals = zeros(Float64, ng)
    nlparamvalues = Float64[]

    Eval.m.nlp_data = JuMP._NLPData(nlobj, nlconstr, nlexpr, nlconstr_duals, nlparamvalues,
                                   user_operators, nx, Eval)
    #opt.nlp_data = MOI.NLPBlockData(opt.nlp_data.constraint_bounds, Eval, has_obj)
    opt.nlp_data = MOI.NLPBlockData(c_bnds, Eval, has_obj)
    println("opt.nlp_data: $(opt.nlp_data)")
    opt.working_evaluator_block = deepcopy(lower_eval_block)
    if (typeof(opt.nlp_data.evaluator) != EAGO.EmptyNLPEvaluator)
        built_evaluator = EAGO.build_nlp_evaluator(MC{nx}, opt.nlp_data.evaluator, opt)         # failing due to type nothing
        opt.working_evaluator_block = MOI.NLPBlockData(opt.nlp_data.constraint_bounds, built_evaluator, has_obj)
    end

    # load lower nlp data block
    if MOI.supports(opt.initial_relaxed_optimizer, MOI.NLPBlock())
        if ~isempty(opt.nonlinear_variable)
            opt.initial_relaxed_optimizer.nlp_data = deepcopy(lower_eval_block)
            opt.working_relaxed_optimizer.nlp_data = deepcopy(lower_eval_block)
        end
    end

    # load upper nlp data block
    upper_eval_block = MOI.NLPBlockData(c_bnds, Eval, has_obj)
    if MOI.supports(opt.initial_relaxed_optimizer, MOI.NLPBlock())
        if ~isempty(opt.nonlinear_variable)
            opt.initial_relaxed_optimizer.nlp_data = deepcopy(upper_eval_block)
            opt.working_relaxed_optimizer.nlp_data = deepcopy(upper_eval_block)
        end
    end
end

#=
SolveScript(f, g, h, xl, xu, opt)
=#
function solve_script(f, xl, xu, opt; g = nothing, gL = nothing, gU = nothing)

    # get dimensions
    @assert length(xl) == length(xu)
    nx = length(xl);
    ng = (g == nothing) ? 0 : length(g(xl))

    # sets most routines to default
    set_to_default!(opt)

    # add variables to lower, upper, and EAGO models
    var_EAGO = MOI.add_variables(opt, nx)
    for j in 1:nx
        MOI.add_constraint(opt, var_EAGO[j], MOI.GreaterThan(xl[j]))
        MOI.add_constraint(opt, var_EAGO[j], MOI.LessThan(xu[j]))
    end

    # Build the JuMP evaluator
    # NEED TO ASSIGN JuMP.num_variables(d.m)
    # NEED TO ASSIGN .m.nlp_data to d_script
    m = JuMP.Model()
    JuMP.@variable(m, xl[i] <= x[i=1:nx] <= xu[i])
    d_script = JuMP.NLPEvaluator(m)
    requested_features = Symbol[]
    def_nlexpr = []
    nlobj, nlconstr_list, def_nlexpr = MOI.initialize(d_script, f, g, ng, nx, requested_features, def_nlexpr)

    has_obj = (f != nothing)
    bnd_pair = MOI.NLPBoundsPair(-Inf, 0.0)
    nlp_bnds = [bnd_pair for i=1:ng]

    println("d_script: $d_script")

    # Optimizes the model with load function
    user_ops = JuMP._Derivatives.UserOperatorRegistry()
    custom_mod_input = (d_script, has_obj, nlp_bnds, nx, ng,
                        nlobj, nlconstr_list, def_nlexpr, user_ops)
    MOI.optimize!(opt, custom_mod! = script_mod, custom_mod_args = custom_mod_input)

    return var_EAGO,opt
end
