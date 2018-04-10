mutable struct callback_storage
    fg_SNOPT_LBD::Function
    col_temp_Ipopt_LBD::Vector{Int32}
    row_temp_Ipopt_LBD::Vector{Int32}
    IPOPT_LBD_eval_f::Function
    IPOPT_LBD_eval_grad_f!::Function

    IPOPT_LBD_eval_g!::Function
    IPOPT_LBD_eval_jac_g!::Function
    IPOPT_LBD_eval_h::Function
    IPOPT_LBD_eval_g::Function
    IPOPT_UBD_eval_g!::Function

    IPOPT_UBD_eval_grad_f!::Function
    IPOPT_UBD_eval_jac_g!::Function
    IPOPT_UBD_eval_h::Function
    fg_SNOPT_UBD::Function
    fg_SNOPT_LBD_Imp

    alphaBB_Tape
    alphaBB_Result
    aBB_IPOPT_LBD_eval_h::Function

    Quadratic_ObjCV_Callback::Function
end
callback_storage() = callback_storage(x->x,[1],[1],x->x,x->x,
                                      x->x,x->x,x->x,x->x,x->x,
                                      x->x,x->x,x->x,x->x,x->x,
                                      [],[],x->x,
                                      x->x)

"""
    MathProgBase.optimize!(s::EAGO_NLP_Model)

Optimizes the `s::EAGO_NLP_Model`. May print console outputs depending on settings
and other solution information is accessible via start solver interface functions.
"""
function MathProgBase.optimize!(s::EAGO_NLP_Model)

    println("Start Optimization Data Structure Setup")
    call_sto = callback_storage()

    # sets the McCormick library options as necessary
    if (s.Opts.solver.LBD_func_relax == "NS-STD-OFF")
        EAGO.set_diff_relax(0)
        EAGO.set_multivar_refine(false,1E-15)
        EAGO.set_subgrad_refine(false)
    elseif (s.Opts.solver.LBD_func_relax == "NS-STD-ON")
        EAGO.set_diff_relax(0)
        EAGO.set_multivar_refine(true,1E-15)
        EAGO.set_subgrad_refine(true)
    elseif (s.Opts.solver.LBD_func_relax == "NS-MV-OFF")
        EAGO.set_diff_relax(0)
        EAGO.set_multivar_refine(false,1E-15)
        EAGO.set_subgrad_refine(false)
    elseif (s.Opts.solver.LBD_func_relax == "NS-MV-ON")
        EAGO.set_diff_relax(0)
        EAGO.set_multivar_refine(true,1E-15)
        EAGO.set_subgrad_refine(true)
    elseif (s.Opts.solver.LBD_func_relax == "Diff1-MV-ON")
        EAGO.set_diff_relax(1)
        EAGO.set_multivar_refine(true,1E-15)
        EAGO.set_subgrad_refine(true)
    elseif (s.Opts.solver.LBD_func_relax == "Diff1-MV-OFF")
        EAGO.set_diff_relax(1)
        EAGO.set_multivar_refine(false,1E-15)
        EAGO.set_subgrad_refine(false)
    elseif (s.Opts.solver.LBD_func_relax == "Diff2-MV-ON")
        EAGO.set_diff_relax(2)
        EAGO.set_multivar_refine(true,1E-15)
        EAGO.set_subgrad_refine(true)
    elseif (s.Opts.solver.LBD_func_relax == "Diff2-MV-OFF")
        EAGO.set_diff_relax(2)
        EAGO.set_multivar_refine(false,1E-15)
        EAGO.set_subgrad_refine(false)
    end

    cols_temp = zeros(Int32,s.Opts.numVar*s.Opts.numConstr)
    rows_temp = zeros(Int32,s.Opts.numVar*s.Opts.numConstr)
    for i = 1:s.Opts.numConstr
      cols_temp[(s.Opts.numVar*(i-1)+1):(s.Opts.numVar*i)] = 1:s.Opts.numVar
      rows_temp[(s.Opts.numVar*(i-1)+1):(s.Opts.numVar*i)] = ones(s.Opts.numVar)*i
    end
    call_sto.col_temp_Ipopt_LBD = cols_temp
    call_sto.row_temp_Ipopt_LBD = rows_temp

    call_sto.IPOPT_LBD_eval_f = (x::Vector{Float64}, X) -> IPOPT_LBD_eval_f(x, X, s.Opts)
    call_sto.IPOPT_LBD_eval_grad_f! = (x::Vector{Float64}, X, f_grad::Vector{Float64}) -> IPOPT_LBD_eval_grad_f!(x, X, f_grad, s.Opts)
    call_sto.IPOPT_LBD_eval_g! = (x::Vector{Float64}, X, g::Vector{Float64}) -> IPOPT_LBD_eval_g!(x, X, g, s.Opts)
    call_sto.IPOPT_LBD_eval_jac_g! = (x::Vector{Float64}, X, mode::Symbol, rows::Vector{Int32}, cols::Vector{Int32}, values::Array{Float64,1}) -> IPOPT_LBD_eval_jac_g!(x, X, mode, rows, cols, values, s.Opts, call_sto)
    call_sto.IPOPT_LBD_eval_h = (x::Vector{Float64}, X, mode::Symbol,
                                 rows::Vector{Int32}, cols::Vector{Int32}, obj_factor::Float64, lambda::Vector{Float64}, values::Array{Float64,1}) -> IPOPT_LBD_eval_h(x, X, mode, rows, cols, obj_factor, lambda, values, s.Opts, call_sto)
    call_sto.IPOPT_LBD_eval_g = (x::Vector{Float64}, X::Vector{Interval{Float64}}) -> IPOPT_LBD_eval_g(x, X, s.Opts)
    call_sto.fg_SNOPT_LBD  = (y::Vector{Float64},X) -> snopt_callback_LBD(y,X,s.Opts)
    call_sto.fg_SNOPT_LBD_Imp  = (y::Vector{Float64},X,param,pmid) -> snopt_callback_LBD_Imp(y,X,s.Opts,param,pmid)
    call_sto.fg_SNOPT_UBD  = (y::Vector{Float64}) -> snopt_callback_UBD(y,s.Opts)
    call_sto.IPOPT_UBD_eval_grad_f! = (x::Vector{Float64}, f_grad::Vector{Float64}) -> IPOPT_UBD_eval_grad_f!(x, f_grad, s.Opts)
    call_sto.IPOPT_UBD_eval_g! = (x::Vector{Float64}, g::Vector{Float64}) -> IPOPT_UBD_eval_g!(x, g, s.Opts)
    call_sto.IPOPT_UBD_eval_jac_g! = (x::Vector{Float64}, mode::Symbol, rows::Vector{Int32}, cols::Vector{Int32}, values::Array{Float64,1}) -> IPOPT_UBD_eval_jac_g!(x, mode, rows, cols, values, s.Opts, call_sto)
    call_sto.IPOPT_UBD_eval_h = (x::Vector{Float64}, mode::Symbol,
                                 rows::Vector{Int32}, cols::Vector{Int32}, obj_factor::Float64, lambda::Vector{Float64}, values::Array{Float64,1}) -> IPOPT_UBD_eval_h(x, mode, rows, cols, obj_factor, lambda, values, s.Opts)

    #=
    X1 = [Interval(0.0,400.0),Interval(0.0,200.0)]
    x1 = mid.(X1)
    out1 = call_sto.IPOPT_LBD_eval_f(x1,X1)
    println("IPOPT_LBD_eval_f: ", out1)
    println("IPOPT_LBD_eval_f type: ", typeof(out1))

    X2 = [Interval(0.0,400.0),Interval(0.0,200.0)]
    x2 = [400.0,200.0]
    out2 = zeros(x1)
    call_sto.IPOPT_LBD_eval_grad_f!(x2,X2,out2)
    println("IPOPT_LBD_eval_grad_f!: ", out2)
    println("IPOPT_LBD_eval_f type: ", typeof(out2))
    =#
 #=
    X3 = [Interval(0.0,400.0),Interval(0.0,200.0)]
    x3 = [400.0,200.0]
    out3 = zeros(2)
    call_sto.IPOPT_LBD_eval_g!(x3,X3,out3)
    println("IPOPT_LBD_eval_g!")
    call_sto.IPOPT_LBD_eval_jac_g!()
    println("IPOPT_LBD_eval_jac_g!")
    call_sto.IPOPT_LBD_eval_h()
    println("IPOPT_LBD_eval_h")
    call_sto.IPOPT_LBD_eval_g()
    println("IPOPT_LBD_eval_g")
    =#

    UBD_error_flag = false

    # checks to see whether an implicit solver should be used
    if (s.Opts.solver.Implicit_Options.flag == true)
        # checks that implicit solver options are valid
        set_Bisect_Func!(s.Opts.solver.BnBSolver,"relative midpoint",s.Opts.solver.Implicit_Options.nx)

        # loads lower problem
        if s.Opts.solver.LBDsolvertype == "LP"
            #temp_LBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
            #            opt,UBD::Float64) -> LP_Relax_LBD_Imp(X,k,pos,opt,UBD)
            #            println("LP Assigned")
            temp_LBP = LP_Relax_LBD_Imp
        elseif s.Opts.solver.LBDsolvertype == "SNOPT"
            temp_LBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                        opt,UBD::Float64) -> SNOPT_LBD_Imp(X,k,pos,opt,UBD)
        elseif s.Opts.solver.LBDsolvertype == "Ipopt"
            temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                        opt,UBD) -> Ipopt_LBD_Imp(X,k,pos,opt,UBD)
        else
            error("Unsupported problem relaxation (Lower Problem).")
        end

        # loads upper problem
        if s.Opts.solver.UBDsolvertype == "LP"
            temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                        opt,UBD) -> LP_Relax_UBD_Imp(X,k,pos,opt,UBD)
                        println("LP Assigned 2")
        else
            UBD_error_flag = true
        end
    else
        # loads lower problem
        if s.Opts.solver.LBDsolvertype == "LP"
            temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> LP_Relax_LBD(X,k,pos,opt,UBD)
        elseif s.Opts.solver.LBDsolvertype == "SNOPT"
            #EAGOSmoothMcCormickGrad.set_outer_rnd(true,1E-9)
            temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> SNOPT_LBD(X,k,pos,opt,UBD)
        elseif s.Opts.solver.LBDsolvertype == "Ipopt"
                        temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> Ipopt_LBD(X,k,pos,opt,UBD)
        elseif s.Opts.solver.LBDsolvertype == "Interval"
            temp_LBP = (X,k::Int64,pos::Int64,opt,UBD) -> Interval_LBD(X,k,pos,opt,UBD)
        elseif s.Opts.solver.LBDsolvertype == "AlphaBB"
            call_sto.alphaBB_Tape = ReverseDiff.HessianTape(s.Opts.f,s.BnBModel.Init_Box)
            call_sto.alphaBB_Result = [Interval(0.0) for i=1:s.Opts.numVar,j=1:s.Opts.numVar]
            call_sto.aBB_IPOPT_LBD_eval_h = (x::Vector{Float64},X::Vector{Interval{Float64}},mode::Symbol,rows::Vector{Int32},cols::Vector{Int32},obj_factor::Float64,
                                             lambda::Vector{Float64},values::Array{Float64,1},opts::EAGO_Inner_NLP,cb::callback_storage,finput::Function) ->
                                             aBB_IPOPT_LBD_eval_h(x,X,mode,rows,cols,obj_factor,lambda,values,opts,cb,finput)
            temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> AlphaBB_LBD(X,k,pos,opt,UBD)
        elseif s.Opts.solver.LBDsolvertype == "Quadratic"
            call_sto.Quadratic_ObjCV_Callback = (x::Vector{Float64},X::Vector{Interval{Float64}},a::Vector{Float64})->Quadratic_ObjCV_Callback(x,X,s.Opts,a)
            temp_LBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                        opt,UBD::Float64) -> Quadratic_LBD(X,k,pos,opt,UBD)
        else
            error("Unsupported problem relaxation.")
        end

        # loads upper problem
        if s.Opts.solver.UBDsolvertype == "LP"
            temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                        opt,UBD) -> LP_Relax_UBD(X,k,pos,opt,UBD)
        elseif s.Opts.solver.UBDsolvertype == "Interval"
            temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                        opt::Any,temp) -> Interval_UBD(X,k,pos,opt,temp)
        else
            UBD_error_flag = true
        end
    end

    if s.Opts.solver.UBDsolvertype == "SNOPT"
        temp_UBP = (X,k::Int64,pos::Int64,opt,UBD) -> SNOPT_UBD(X,k,pos,opt,UBD)
    elseif s.Opts.solver.UBDsolvertype == "Ipopt"
        temp_UBP = (X,k::Int64,pos::Int64,opt,UBD) -> Ipopt_UBD(X,k,pos,opt,UBD)
    elseif s.Opts.solver.UBDsolvertype == "MPBNonlinear"
        #temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,opt,UBD) -> NLPSolve_UBD(X,k,pos,opt,UBD)
        temp_UBP = NLPSolve_UBD
    elseif UBD_error_flag
        error("Unsupported problem relaxation (Upper Problem).")
    end
    s.Opts.solver.BnBSolver.Lower_Prob = deepcopy(temp_LBP)
    s.Opts.solver.BnBSolver.Upper_Prob = deepcopy(temp_UBP)

    s.Opts.solver.BnBSolver.Preprocess = composite_DR_pre
    s.Opts.solver.BnBSolver.Postprocess = composite_DR_post

    s.Opts.solver.BnBSolver.opt = deepcopy([s.Opts,call_sto])
    s.Opts.solver.BnBSolver.Verbosity = s.Opts.solver.verbosity
    println("End Optimization Data Structure Setup")

    println("ran to solve BnB")
    solveBnB!(s.Opts.solver.BnBSolver,s.BnBModel)
    println("ran post solve BnB")
    if (s.Opts.sense == :Max)
        s.BnBModel.soln_val = -s.BnBModel.soln_val
    end

    if (s.BnBModel.feas_fnd == true)
        s.status = :Optimal
    else
        s.status = :Infeasible
    end
end
