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
callback_storage() = callback_storage(x->x,[1],[1],x->x,x->x,x->x,x->x,x->x,x->x,x->x,x->x,x->x,x->x,x->x,x->x,[],[],x->x,x->x)

"""
    loadproblem!(m::EAGO_NLP_Model, nvar::Int64, ncon::Int64,xL::Vector{Float64},
                 xU::Vector{Float64}, gL::Vector{Float64}, gU::Vector{Float64},
                 sense::Symbol, d::MathProgBase.AbstractNLPEvaluator)

This variant is primarily used by the JuMP interface. Loads the model `m` with
the optimization parameters given by:
* `m::EAGO_NLP_Model`: The model type supported by EAGO_NLP_Solver
* `nvar::Int64`: Number of variables in the model
* `ncon::Int64`: Number of constraints in the model
* `xL::Vector{Float64}`: Variable lower bound
* `xU::Vector{Float64}`: Variable upper bound
* `gL::Vector{Float64}`: Constraint lower bound
* `gU::Vector{Float64}`: Constraint upper bound
* `sense::Symbol`: Min or Max
* `d::MathProgBase.AbstractNLPEvaluator`: Use for JuMP interface
"""
function MathProgBase.loadproblem!(m::EAGO_NLP_Model, nvar::Int64, ncon::Int64,
                                   xL::Vector{Float64},
                                   xU::Vector{Float64},
                                   gL::Vector{Float64},
                                   gU::Vector{Float64},
                                   sense::Symbol, d::MathProgBase.AbstractNLPEvaluator)

         #println("Began Loading Problem")
         @assert nvar == length(xL) == length(xU)
         @assert ncon == length(gL) == length(gU)

         if m.Opts.solver.UBDsolvertype == "MPBNonlinear"
             loadproblem!(m.Opts.UBDmodel,nvar,ncon,xL,xU,gL,gU,sense,d)
         end

         # counts number of inequality constraints (NonInf)
         m.Opts.gexp = count(x->x<Inf,gU)+count(x->x>-Inf,gL)
         m.Opts.gL_loc = filter!(x->x!=-1,[ gL[i]>-Inf ? i : -1 for i=1:ncon])
         m.Opts.gU_loc = filter!(x->x!=-1,[ gU[i]<Inf ? i : -1 for i=1:ncon])

         m.Opts.gL, m.Opts.gU = gL, gU
         if m.Opts.validated
                  m.BnBModel.Init_Box = [Interval(xL[i],xU[i]) for i=1:nvar]
                  m.BnBModel.box = [[Interval(xL[i],xU[i]) for i=1:nvar]]
         else
                  m.BnBModel.Init_Box = [MCInterval(xL[i],xU[i]) for i=1:nvar]
                  m.BnBModel.box = [[MCInterval(xL[i],xU[i]) for i=1:nvar]]
         end
         m.Opts.sense = sense

         m.Opts.numVar, m.Opts.numConstr = nvar, ncon
         m.d = d
         m.Opts.d = d

         MathProgBase.initialize(d,[:ExprGraph])
         m.Opts.obj = verify_support(MathProgBase.obj_expr(d))
         m.Opts.vartypes = fill(:Cont,nvar)

         m.Opts.constrs = map(1:ncon) do c
                  verify_support(MathProgBase.constr_expr(d,c))
         end

         # sets up DAG contractor
         expr_Array::Vector{Expr} = Expr[]
         gLt::Vector{Float64} = [-Inf for i=1:ncon]
         gUt::Vector{Float64} = [Inf for i=1:ncon]
         for i=1:ncon
             # sets up double-sided constraint
             if length(m.Opts.constrs[i].args) == 3
                   if (m.Opts.constrs[i].args[2] == :<)
                            gUt[i] = m.Opts.constrs[i].args[3]
                   elseif (m.Opts.constrs[i].args[2] == :>)
                            gLt[i] = m.Opts.constrs[i].args[3]
                   else
                            gUt[i] = m.Opts.constrs[i].args[3]
                            gLt[i] = m.Opts.constrs[i].args[3]
                   end
                   push!(expr_Array,m.Opts.constrs[i].args[2])
             # sets up single-sided constraints
             elseif length(m.Opts.constrs[i].args) == 5
                      gLt[i] = m.Opts.constrs[i].args[1]
                      gUt[i] = m.Opts.constrs[i].args[5]
                      push!(expr_Array,m.Opts.constrs[i].args[3])
             end
         end
         isempty(expr_Array) && (m.Opts.solver.DAG_depth = -1)

         if (m.Opts.solver.DAG_depth>0)
             if (m.Opts.solver.validated)
                 Generate_TapeList(expr_Array,nvar,gL,gU,Interval{Float64})
                 m.Opts.DAG_tlist = Generate_TapeList(expr_Array,nvar,gL,gU,Interval{Float64})
             else
                 Generate_TapeList(expr_Array,nvar,gL,gU,MCInterval{Float64})
                 m.Opts.DAG_tlist = Generate_TapeList(expr_Array,nvar,gL,gU,MCInterval{Float64})
             end
         end

         # sets implicit function values as appropriate
         m.Opts.Imp_np = m.Opts.numVar - m.Opts.Imp_nx
         m.Opts.solver.PSmcOpt.np = m.Opts.Imp_np
         m.Opts.solver.PSmcOpt.nx = m.Opts.Imp_nx
         m.Opts.solver.PIntOpt.nx = m.Opts.Imp_nx

         if sense == :Min
                  @eval f(x) = $(m.Opts.obj)
         else
                  @eval f(x) = -$(m.Opts.obj)
         end
         g_arr = Expr(:vect)
         g_arr.args = expr_Array
         @eval g(x) = $g_arr
         m.Opts.f = x -> Base.invokelatest(f,x)
         m.Opts.g = x -> Base.invokelatest(g,x)
         #println("End Load!")

         #println("Start Optimization Data Structure Setup")
         call_sto = callback_storage()

             cols_temp::VecOrMat{Int32} = zeros(Int32,m.Opts.numVar*m.Opts.numConstr)
             rows_temp::VecOrMat{Int32} = zeros(Int32,m.Opts.numVar*m.Opts.numConstr)
             for i = 1:m.Opts.numConstr
               cols_temp[(m.Opts.numVar*(i-1)+1):(m.Opts.numVar*i)] = 1:m.Opts.numVar
               rows_temp[(m.Opts.numVar*(i-1)+1):(m.Opts.numVar*i)] = ones(m.Opts.numVar)*i
             end
             call_sto.col_temp_Ipopt_LBD = cols_temp
             call_sto.row_temp_Ipopt_LBD = rows_temp

             call_sto.IPOPT_LBD_eval_f = (x::Vector{Float64}, X) -> IPOPT_LBD_eval_f(x, X, m.Opts)
             call_sto.IPOPT_LBD_eval_grad_f! = (x::Vector{Float64}, X, f_grad::Vector{Float64}) -> IPOPT_LBD_eval_grad_f!(x, X, f_grad, m.Opts)
             call_sto.IPOPT_LBD_eval_g! = (x::Vector{Float64}, X, g::Vector{Float64}) -> IPOPT_LBD_eval_g!(x, X, g, m.Opts)
             call_sto.IPOPT_LBD_eval_jac_g! = (x::Vector{Float64}, X, mode::Symbol, rows::Vector{Int32}, cols::Vector{Int32}, values::Array{Float64,1}) -> IPOPT_LBD_eval_jac_g!(x, X, mode, rows, cols, values, m.Opts, call_sto)
             call_sto.IPOPT_LBD_eval_h = (x::Vector{Float64}, X, mode::Symbol,rows::Vector{Int32}, cols::Vector{Int32}, obj_factor::Float64, lambda::Vector{Float64}, values::Array{Float64,1}) -> IPOPT_LBD_eval_h(x, X, mode, rows, cols, obj_factor, lambda, values, m.Opts, call_sto)
             call_sto.IPOPT_LBD_eval_g = (x::Vector{Float64}, X::Vector{Interval{Float64}}) -> IPOPT_LBD_eval_g(x, X, m.Opts)
             call_sto.fg_SNOPT_LBD  = (y::Vector{Float64},X) -> snopt_callback_LBD(y,X,m.Opts)
             call_sto.fg_SNOPT_LBD_Imp  = (y::Vector{Float64},X,param,pmid) -> snopt_callback_LBD_Imp(y,X,m.Opts,param,pmid)
             call_sto.fg_SNOPT_UBD  = (y::Vector{Float64}) -> snopt_callback_UBD(y,m.Opts)
             call_sto.IPOPT_UBD_eval_grad_f! = (x::Vector{Float64}, f_grad::Vector{Float64}) -> IPOPT_UBD_eval_grad_f!(x, f_grad, m.Opts)
             call_sto.IPOPT_UBD_eval_g! = (x::Vector{Float64}, g::Vector{Float64}) -> IPOPT_UBD_eval_g!(x, g, m.Opts)
             call_sto.IPOPT_UBD_eval_jac_g! = (x::Vector{Float64}, mode::Symbol, rows::Vector{Int32}, cols::Vector{Int32}, values::Array{Float64,1}) -> IPOPT_UBD_eval_jac_g!(x, mode, rows, cols, values, m.Opts, call_sto)
             call_sto.IPOPT_UBD_eval_h = (x::Vector{Float64}, mode::Symbol,rows::Vector{Int32}, cols::Vector{Int32}, obj_factor::Float64, lambda::Vector{Float64}, values::Array{Float64,1}) -> IPOPT_UBD_eval_h(x, mode, rows, cols, obj_factor, lambda, values, m.Opts)

             UBD_error_flag = false

             # checks to see whether an implicit solver should be used
             if (m.Opts.solver.ImplicitFlag == true)
                 # checks that implicit solver options are valid
                 set_Bisect_Func!(m.Opts.solver.BnBSolver,"relative midpoint",m.Opts.Imp_nx)

                 # loads lower problem
                 if m.Opts.solver.LBDsolvertype == "LP"
                     #temp_LBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                     #            opt,UBD::Float64) -> LP_Relax_LBD_Imp(X,k,pos,opt,UBD)
                     #            println("LP Assigned")
                     temp_LBP = LP_Relax_LBD_Imp
                 elseif m.Opts.solver.LBDsolvertype == "SNOPT"
                     temp_LBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt,UBD::Float64) -> SNOPT_LBD_Imp(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.LBDsolvertype == "Ipopt"
                     temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt,UBD) -> Ipopt_LBD_Imp(X,k,pos,opt,UBD)
                 else
                     error("Unsupported problem relaxation (Lower Problem).")
                 end

                 # loads upper problem
                 if m.Opts.solver.UBDsolvertype == "Interval"
                     temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt,UBD) -> Imp_Interval_UBD(X,k,pos,opt,UBD)
                                 println("LP Assigned 2")
                 else
                     UBD_error_flag = true
                 end
             else
                 # loads lower problem
                 if m.Opts.solver.LBDsolvertype == "LP"
                     temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> LP_Relax_LBD(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.LBDsolvertype == "SNOPT"
                     #EAGOSmoothMcCormickGrad.set_outer_rnd(true,1E-9)
                     temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> SNOPT_LBD(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.LBDsolvertype == "Ipopt"
                     #temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> Ipopt_LBD(X,k,pos,opt,UBD)
                     temp_LBP = Ipopt_LBD
                 elseif m.Opts.solver.LBDsolvertype == "Interval"
                     temp_LBP = (X,k::Int64,pos::Int64,opt,UBD) -> Interval_LBD(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.LBDsolvertype == "AlphaBB"
                     call_sto.alphaBB_Tape = ReverseDiff.HessianTape(m.Opts.f,m.BnBModel.Init_Box)
                     call_sto.alphaBB_Result = [Interval(0.0) for i=1:m.Opts.numVar,j=1:m.Opts.numVar]
                     call_sto.aBB_IPOPT_LBD_eval_h = (x::Vector{Float64},X::Vector{Interval{Float64}},mode::Symbol,rows::Vector{Int32},cols::Vector{Int32},obj_factor::Float64,
                                                      lambda::Vector{Float64},values::Array{Float64,1},opts::EAGO_Inner_NLP,cb::callback_storage,finput::Function) ->
                                                      aBB_IPOPT_LBD_eval_h(x,X,mode,rows,cols,obj_factor,lambda,values,opts,cb,finput)
                     temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> AlphaBB_LBD(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.LBDsolvertype == "Quadratic"
                     call_sto.Quadratic_ObjCV_Callback = (x::Vector{Float64},X::Vector{Interval{Float64}},a::Vector{Float64})->Quadratic_ObjCV_Callback(x,X,s.Opts,a)
                     temp_LBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt,UBD::Float64) -> Quadratic_LBD(X,k,pos,opt,UBD)
                 else
                     error("Unsupported problem relaxation.")
                 end

                 # loads upper problem
                 if m.Opts.solver.UBDsolvertype == "LP"
                     temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt,UBD) -> LP_Relax_UBD(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.UBDsolvertype == "Interval"
                     temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt::Any,temp) -> Interval_UBD(X,k,pos,opt,temp)
                 else
                     UBD_error_flag = true
                 end
             end

             if m.Opts.solver.UBDsolvertype == "SNOPT"
                 temp_UBP = (X,k::Int64,pos::Int64,opt,UBD) -> SNOPT_UBD(X,k,pos,opt,UBD)
             elseif m.Opts.solver.UBDsolvertype == "Ipopt"
                 temp_UBP = (X,k::Int64,pos::Int64,opt,UBD) -> Ipopt_UBD(X,k,pos,opt,UBD)
             elseif m.Opts.solver.UBDsolvertype == "MPBNonlinear"
                 #temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,opt,UBD) -> NLPSolve_UBD(X,k,pos,opt,UBD)
                 temp_UBP = NLPSolve_UBD
             elseif UBD_error_flag
                 error("Unsupported problem relaxation (Upper Problem).")
             end
             m.Opts.solver.BnBSolver.Lower_Prob = deepcopy(temp_LBP)
             m.Opts.solver.BnBSolver.Upper_Prob = deepcopy(temp_UBP)
             m.Opts.solver.BnBSolver.Preprocess = composite_DR_pre
             m.Opts.solver.BnBSolver.Postprocess = composite_DR_post

             m.Opts.solver.BnBSolver.opt = deepcopy(Any[call_sto])
             m.Opts.solver.BnBSolver.Verbosity = m.Opts.solver.verbosity

         return m
end

"""
    loadproblem!(m::EAGO_NLP_Model, nvar::Int64, ncon::Int64,xL::Vector{Float64},
                 xU::Vector{Float64}, gL::Vector{Float64}, gU::Vector{Float64},
                 f,g)

This variant is primarily used by the JuMP interface. Loads the model `m` with
the optimization parameters given by:
* `m::EAGO_NLP_Model`: The model type supported by EAGO_NLP_Solver
* `nvar::Int64`: Number of variables in the model
* `ncon::Int64`: Number of constraints in the model
* `xL::Vector{Float64}`: Variable lower bound
* `xU::Vector{Float64}`: Variable upper bound
* `gL::Vector{Float64}`: Constraint lower bound
* `gU::Vector{Float64}`: Constraint upper bound
* `sense::Symbol`: Min or Max
* `f`: Objective function
* `g`: Constraint function
"""
function MathProgBase.loadproblem!(m::EAGO_NLP_Model, nvar::Int64, ncon::Int64,
                                   xL::Vector{Float64},
                                   xU::Vector{Float64},
                                   gL,
                                   gU,
                                   sense::Symbol, f, g)

         @assert nvar == length(xL) == length(xU)
         @assert ncon == length(gL) == length(gU)

         m.Opts.gexp = count(x->x<Inf,gU)+count(x->x>-Inf,gL)
         m.Opts.gL_loc = filter!(x->x!=-1,[ gL[i]>-Inf ? i : -1 for i=1:ncon])
         m.Opts.gU_loc = filter!(x->x!=-1,[ gU[i]<Inf ? i : -1 for i=1:ncon])

         m.Opts.gL, m.Opts.gU = gL, gU
         if m.Opts.validated
                  m.BnBModel.Init_Box = [Interval(xL[i],xU[i]) for i=1:nvar]
                  m.BnBModel.box = [[Interval(xL[i],xU[i]) for i=1:nvar]]
         else
                  m.BnBModel.Init_Box = [MCInterval(xL[i],xU[i]) for i=1:nvar]
                  m.BnBModel.box = [[MCInterval(xL[i],xU[i]) for i=1:nvar]]
         end
         m.Opts.sense = sense

         m.Opts.numVar, m.Opts.numConstr = nvar, ncon

         m.Opts.vartypes = fill(:Cont,nvar)

         # sets implicit function values as appropriate

         if sense == :Min
                  m.Opts.f = x -> f(x)
         elseif sense == :Max
                  m.Opts.f = x -> -f(x)
         else
                  error("Sense must be :Max or :Min")
         end
         m.Opts.g = x -> g(x)

         call_sto = callback_storage()

             cols_temp::VecOrMat{Int32} = zeros(Int32,m.Opts.numVar*m.Opts.numConstr)
             rows_temp::VecOrMat{Int32} = zeros(Int32,m.Opts.numVar*m.Opts.numConstr)
             for i = 1:m.Opts.numConstr
               cols_temp[(m.Opts.numVar*(i-1)+1):(m.Opts.numVar*i)] = 1:m.Opts.numVar
               rows_temp[(m.Opts.numVar*(i-1)+1):(m.Opts.numVar*i)] = ones(m.Opts.numVar)*i
             end
             call_sto.col_temp_Ipopt_LBD = cols_temp
             call_sto.row_temp_Ipopt_LBD = rows_temp

             call_sto.IPOPT_LBD_eval_f = (x::Vector{Float64}, X) -> IPOPT_LBD_eval_f(x, X, s.Opts)
             call_sto.IPOPT_LBD_eval_grad_f! = (x::Vector{Float64}, X, f_grad::Vector{Float64}) -> IPOPT_LBD_eval_grad_f!(x, X, f_grad, m.Opts)
             call_sto.IPOPT_LBD_eval_g! = (x::Vector{Float64}, X, g::Vector{Float64}) -> IPOPT_LBD_eval_g!(x, X, g, m.Opts)
             call_sto.IPOPT_LBD_eval_jac_g! = (x::Vector{Float64}, X, mode::Symbol, rows::Vector{Int32}, cols::Vector{Int32}, values::Array{Float64,1}) -> IPOPT_LBD_eval_jac_g!(x, X, mode, rows, cols, values, m.Opts, call_sto)
             call_sto.IPOPT_LBD_eval_h = (x::Vector{Float64}, X, mode::Symbol,rows::Vector{Int32}, cols::Vector{Int32}, obj_factor::Float64, lambda::Vector{Float64}, values::Array{Float64,1}) -> IPOPT_LBD_eval_h(x, X, mode, rows, cols, obj_factor, lambda, values, m.Opts, call_sto)
             call_sto.IPOPT_LBD_eval_g = (x::Vector{Float64}, X::Vector{Interval{Float64}}) -> IPOPT_LBD_eval_g(x, X, m.Opts)
             call_sto.fg_SNOPT_LBD  = (y::Vector{Float64},X) -> snopt_callback_LBD(y,X,s.Opts)
             call_sto.fg_SNOPT_LBD_Imp  = (y::Vector{Float64},X,param,pmid) -> snopt_callback_LBD_Imp(y,X,m.Opts,param,pmid)
             call_sto.fg_SNOPT_UBD  = (y::Vector{Float64}) -> snopt_callback_UBD(y,s.Opts)
             call_sto.IPOPT_UBD_eval_grad_f! = (x::Vector{Float64}, f_grad::Vector{Float64}) -> IPOPT_UBD_eval_grad_f!(x, f_grad, m.Opts)
             call_sto.IPOPT_UBD_eval_g! = (x::Vector{Float64}, g::Vector{Float64}) -> IPOPT_UBD_eval_g!(x, g, m.Opts)
             call_sto.IPOPT_UBD_eval_jac_g! = (x::Vector{Float64}, mode::Symbol, rows::Vector{Int32}, cols::Vector{Int32}, values::Array{Float64,1}) -> IPOPT_UBD_eval_jac_g!(x, mode, rows, cols, values, m.Opts, call_sto)
             call_sto.IPOPT_UBD_eval_h = (x::Vector{Float64}, mode::Symbol,rows::Vector{Int32}, cols::Vector{Int32}, obj_factor::Float64, lambda::Vector{Float64}, values::Array{Float64,1}) -> IPOPT_UBD_eval_h(x, mode, rows, cols, obj_factor, lambda, values, m.Opts)

             UBD_error_flag = false

             # checks to see whether an implicit solver should be used
             if (m.Opts.solver.ImplicitFlag == true)
                 # checks that implicit solver options are valid
                 set_Bisect_Func!(m.Opts.solver.BnBSolver,"relative midpoint",m.Opts.Imp_nx)

                 # loads lower problem
                 if m.Opts.solver.LBDsolvertype == "LP"
                     #temp_LBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                     #            opt,UBD::Float64) -> LP_Relax_LBD_Imp(X,k,pos,opt,UBD)
                     #            println("LP Assigned")
                     temp_LBP = LP_Relax_LBD_Imp
                 elseif m.Opts.solver.LBDsolvertype == "SNOPT"
                     temp_LBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt,UBD::Float64) -> SNOPT_LBD_Imp(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.LBDsolvertype == "Ipopt"
                     temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt,UBD) -> Ipopt_LBD_Imp(X,k,pos,opt,UBD)
                 else
                     error("Unsupported problem relaxation (Lower Problem).")
                 end

                 # loads upper problem
                 if m.Opts.solver.UBDsolvertype == "Interval"
                     temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt,UBD) -> Imp_Interval_UBD(X,k,pos,opt,UBD)
                                 println("LP Assigned 2")
                 else
                     UBD_error_flag = true
                 end
             else
                 # loads lower problem
                 if m.Opts.solver.LBDsolvertype == "LP"
                     temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> LP_Relax_LBD(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.LBDsolvertype == "SNOPT"
                     #EAGOSmoothMcCormickGrad.set_outer_rnd(true,1E-9)
                     temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> SNOPT_LBD(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.LBDsolvertype == "Ipopt"
                     temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> Ipopt_LBD(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.LBDsolvertype == "Interval"
                     temp_LBP = (X,k::Int64,pos::Int64,opt,UBD) -> Interval_LBD(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.LBDsolvertype == "AlphaBB"
                     call_sto.alphaBB_Tape = ReverseDiff.HessianTape(m.Opts.f,m.BnBModel.Init_Box)
                     call_sto.alphaBB_Result = [Interval(0.0) for i=1:m.Opts.numVar,j=1:m.Opts.numVar]
                     call_sto.aBB_IPOPT_LBD_eval_h = (x::Vector{Float64},X::Vector{Interval{Float64}},mode::Symbol,rows::Vector{Int32},cols::Vector{Int32},obj_factor::Float64,
                                                      lambda::Vector{Float64},values::Array{Float64,1},opts::EAGO_Inner_NLP,cb::callback_storage,finput::Function) ->
                                                      aBB_IPOPT_LBD_eval_h(x,X,mode,rows,cols,obj_factor,lambda,values,opts,cb,finput)
                     temp_LBP = (X,k::Int64,pos::Int64,opt,UBD::Float64) -> AlphaBB_LBD(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.LBDsolvertype == "Quadratic"
                     call_sto.Quadratic_ObjCV_Callback = (x::Vector{Float64},X::Vector{Interval{Float64}},a::Vector{Float64})->Quadratic_ObjCV_Callback(x,X,s.Opts,a)
                     temp_LBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt,UBD::Float64) -> Quadratic_LBD(X,k,pos,opt,UBD)
                 else
                     error("Unsupported problem relaxation.")
                 end

                 # loads upper problem
                 if m.Opts.solver.UBDsolvertype == "LP"
                     temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt,UBD) -> LP_Relax_UBD(X,k,pos,opt,UBD)
                 elseif m.Opts.solver.UBDsolvertype == "Interval"
                     temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                                 opt::Any,temp) -> Interval_UBD(X,k,pos,opt,temp)
                 else
                     UBD_error_flag = true
                 end
             end

             if m.Opts.solver.UBDsolvertype == "SNOPT"
                 temp_UBP = (X,k::Int64,pos::Int64,opt,UBD) -> SNOPT_UBD(X,k,pos,opt,UBD)
             elseif m.Opts.solver.UBDsolvertype == "Ipopt"
                 temp_UBP = (X,k::Int64,pos::Int64,opt,UBD) -> Ipopt_UBD(X,k,pos,opt,UBD)
             elseif m.Opts.solver.UBDsolvertype == "MPBNonlinear"
                 #temp_UBP = (X::Vector{Interval{Float64}},k::Int64,pos::Int64,opt,UBD) -> NLPSolve_UBD(X,k,pos,opt,UBD)
                 temp_UBP = NLPSolve_UBD
             elseif UBD_error_flag
                 error("Unsupported problem relaxation (Upper Problem).")
             end
             m.Opts.solver.BnBSolver.Lower_Prob = deepcopy(temp_LBP)
             m.Opts.solver.BnBSolver.Upper_Prob = deepcopy(temp_UBP)
             m.Opts.solver.BnBSolver.Preprocess = composite_DR_pre
             m.Opts.solver.BnBSolver.Postprocess = composite_DR_post

             m.Opts.solver.BnBSolver.opt = deepcopy(Any[call_sto])
             m.Opts.solver.BnBSolver.Verbosity = m.Opts.solver.verbosity

         return m
end

function addEqnConstrDAG!(s::EAGO_NLP_Model,expr)
         s.Opts.eqn_constrs = expr
end
