type ImplicitModel
    f::Function
    g::Function
    h::Function
    hj::Function
    nx::Int64
    Imp_RR_depth::Int64
    Imp_probe_depth::Int64
end
ImplicitModel() = ImplicitModel(x->x,x->x,x->x,x->x,0,-1,-1)


"""
    EAGO_Inner_NLP

This is the inner storage object for the model interface to JuMP and MathProgBase.
Effectively empty when initialized. It has the following fields:
* `gL::Vector{Float64}`: Lower constraint bounds
* `gU::Vector{Float64}`: Upper constraint bounds
* `gexp::Int64`: Expression for constraints
* `gL_loc::Vector{Int64}`: Index at which constraints have finite lower bounds
* `gU_loc::Vector{Int64}`: Index at which constraints have finite upper bounds
* `numVar::Int64`: Number of variables in full problem.
* `numConstr::Int64`: Number of constraints in full problem.
* `vartypes::Vector{Symbol}`: Type of variable array (currently not used)
* `sense::Symbol`: Minimization (:Min) or maximization (:Max)
* `obj::Expr`: Expression for objective provide by AbstractNLPEvaluator
* `constrs::Vector{Expr}`: Expressions for constraints provide by AbstractNLPEvaluator
* `DAG_tlist::TapeList`: Tape list for DAG propagation on constraints
* `f::Function`: Objective function
* `g::Function`: Constraint function
* `solver::EAGO_NLPSolver`: Solve storage object
* `impmod::ImplicitModel`: Implicit model storage object (not currently used)
"""
type EAGO_Inner_NLP
    gL::Vector{Float64}
    gU::Vector{Float64}
    gexp::Int64
    gL_loc::Vector{Int64}
    gU_loc::Vector{Int64}
    numVar::Int64
    numConstr::Int64
    vartypes::Vector{Symbol}
    sense::Symbol
    obj::Expr
    constrs::Vector{Expr}
    DAG_tlist::TapeList
    f::Function
    g::Function
    solver::EAGO_NLPSolver
    impmod::ImplicitModel
    UBDmodel
    validated
    d
end
EAGO_Inner_NLP(s::EAGO_NLPSolver) = EAGO_Inner_NLP([0.0],
                                                    [0.0],
                                                    0,
                                                    [0],
                                                    [0],
                                                    0,
                                                    0,
                                                    Symbol[],
                                                    Symbol(),
                                                    Expr(:call),
                                                    Expr[],
                                                    TapeList(),
                                                    x -> x,
                                                    x -> x,
                                                    s,
                                                    ImplicitModel(),
                                                    [],
                                                    false,
                                                    nothing)

"""
    EAGO_NLP_Model

The model type to interface with JuMP and MathProgBase. This has fields:
* `BnBModel::BnBModel`: which stores information from the Branch and Bound Problem.
* `Opts::EAGO_Inner_NLP`: storage for various problem descriptors
* `d`: The storage for an AbstractNLPEvaluator
* `status::Symbol`: descriptor for model status (e.g. :optimal)
"""
type EAGO_NLP_Model <: MathProgBase.AbstractNonlinearModel
    BnBModel::BnBModel
    Opts::EAGO_Inner_NLP
    d#::AbstractNLPEvaluator
    status::Symbol
end

function MathProgBase.NonlinearModel(s::EAGO_NLPSolver)
    eiNLP = EAGO_Inner_NLP(s)
    eiNLP.UBDmodel = NonlinearModel(s.UBDsolver)
    if s.validated == true
        bnb_mod = BnBModel([Interval{Float64}(0.0,1.0)])
    else
        bnb_mod = BnBModel([MCInterval{Float64}(0.0,1.0)])
    end
    EAGOmodel = EAGO_NLP_Model(bnb_mod,eiNLP,nothing,:Uninitialized)
    EAGOmodel.Opts.validated = s.validated
    return EAGOmodel
end

LinearQuadraticModel(s::EAGO_NLPSolver) = NonlinearToLPQPBridge(NonlinearModel(s))
