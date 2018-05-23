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
         #println("Start Load!")

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
         expr_Array = Expr[]
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
         #println("Finished Loading Problem")
         #println("End Load!")
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

         return m
end

function addEqnConstrDAG!(s::EAGO_NLP_Model,expr)
         s.Opts.eqn_constrs = expr
end
