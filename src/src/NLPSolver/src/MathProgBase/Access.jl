MathProgBase.status(m::EAGO_NLP_Model) = m.status
MathProgBase.getsolution(m::EAGO_NLP_Model) = m.BnBModel.soln
MathProgBase.getobjval(m::EAGO_NLP_Model) = m.BnBModel.UBDg
MathProgBase.getobjbound(m::EAGO_NLP_Model) = m.BnBModel.UBDg
getfeasibility(m::EAGO_NLP_Model) = m.BnBModel.feas_fnd
MathProgBase.getsolvetime(m::EAGO_NLP_Model) = LBDtime(m.BnBModel) + UBDtime(m.BnBModel)
MathProgBase.setwarmstart!(m::EAGO_NLP_Model,x) = []
MathProgBase.setsense!(m::EAGO_NLP_Model,sense) = (m.Opts.sense = sense)
MathProgBase.getsense(m::EAGO_NLP_Model) = m.Opts.sense

"""
    set_Implicit_Model!(m::EAGO.EAGO_NLP_Model,...)

Sets the functions and dimensionality of the implicit model corresponding to
the EAGO_NLP_model.
"""
function set_Implicit_Model!(m::EAGO.EAGO_NLP_Model,f::Function,h::Function,
                             hj::Function,g::Union,nx::Int64)
    m.Imp_f = f
    m.Imp_g = g
    m.Imp_h = h
    m.Imp_hj = hj
    m.Imp_nx = nx
end

"""
    Solve_Implicit(jm::JuMP.Model)

Sets the functions and dimensionality of the implicit model corresponding to
the JuMP model containing the EAGO_NLP_model.
"""
function Solve_Implicit(jm::JuMP.Model,f::Function,h::Function,
                             hj::Function,g::Function,nx::Int64;
                             Imp_gL_Loc = [Int64(0)],
                             Imp_gU_Loc = [Int64(0)],
                             Imp_gL = [Float64(0)],
                             Imp_gU = [Float64(0)],
                             Imp_nCons = Int64(0))

    # Builds the JuMP model and adds the appropriate equations
    JuMP.build(jm)
    m = internalmodel(jm)
    m.Opts.Imp_f = f
    m.Opts.Imp_g = g
    m.Opts.Imp_h = h
    m.Opts.Imp_hj = hj
    m.Opts.Imp_nx = nx
    m.Opts.Imp_gL_Loc = Imp_gL_Loc
    m.Opts.Imp_gU_Loc = Imp_gU_Loc
    m.Opts.Imp_gL = Imp_gL
    m.Opts.Imp_gU = Imp_gU
    m.Opts.Imp_nCons = Imp_nCons

    # Unpacks inner model and solves
    MathProgBase.optimize!(m)

    # Repacks the JuMP model with appropriate info
    return MathProgBase.status(m)
end

"""
    Optimize_Script(f,xL,xU;g=[],h=[],p=[],mode=:Min,solver=EAGO_NLPSolver())

Sets the functions and dimensionality of the implicit model corresponding to
the JuMP model containing the EAGO_NLP_model.
"""
function Optimize_Script(f::Function,xL::Vector{Float64},xU::Vector{Float64};
                         g = [], h = [], p = [], mode::Symbol = :Min, solver = EAGO_NLPSolver())

    @assert length(xL) == length(xU)
    nx = length(xL)

    # creates functions with parameterization if necessary
    if (p != [])
        fp = x -> f(x,p)
        (g != []) && (gp = x -> g(x,p))
        (h != []) && (hp = x -> h(x,p))
    else
        fp = x -> f(x)
        (g != []) && (gp = x -> g(x))
        (h != []) && (hp = x -> h(x))
    end

    # converts equality constrains to inequality constraints and sets bounds
    if ((g == []) && (h == []))
        gh = []; gL = Float64[]; gU = Float64[]; glen = 0
    elseif (h == [])
        gh = x -> gp(x)
        glen = length(gh(xL))
        gL = Float64[-Inf for i=1:glen]
        gU = Float64[0.0 for i=1:glen]
    elseif (g == [])
        gh = x -> vcat(-hp(x),hp(x))
        glen = length(gh(xL))
        gL = Float64[-Inf for i=1:glen]
        gU = Float64[0.0 for i=1:glen]
    else
        gh = x -> vcat(-hp(x),hp(x),gp(x))
        glen = length(gh(xL))
        gL = Float64[-Inf for i=1:glen]
        gU = Float64[0.0 for i=1:glen]
    end

    # creates and loads model
    solver.UBDsolvertype = "Ipopt"
    m = MathProgBase.NonlinearModel(solver)
    MathProgBase.loadproblem!(m, nx, glen, xL, xU, gL, gU, mode, fp, gh)
    MathProgBase.optimize!(m)

    # extracts objective value, solution, feasibility
    stat = MathProgBase.status(m)
    soln = MathProgBase.getsolution(m)
    objv = MathProgBase.getobjval(m)
    return objv,soln,stat
end

#=
"""
    Save_History(m::EAGO_NLP_Model)
"""
function Save_History(f::String,m::EAGO_NLP_Model)
     if (EAGO_NLP_Model.status == :Optimal)
     elseif (EAGO_NLP_Model.status == :Infeasible)
     end
end

=#
