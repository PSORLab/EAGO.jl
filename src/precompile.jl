function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(MathOptInterface.optimize!),Optimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt}})   # time: 11.415445
    Base.precompile(Tuple{typeof(reform_epigraph_min!),GlobalOptimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt},ParsedProblem,BufferedNonlinearFunction{1, NS}})   # time: 1.4057246
    Base.precompile(Tuple{typeof(forward_pass!),Evaluator,BufferedNonlinearFunction{2, NS}})   # time: 0.881403
    Base.precompile(Tuple{Type{Optimizer}})   # time: 0.5019401
    Base.precompile(Tuple{typeof(MathOptInterface.empty!),Optimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt}})   # time: 0.0625456
    Base.precompile(Tuple{typeof(initialize!),RelaxCache{1, NS},DirectedTree})   # time: 0.0522575
    Base.precompile(Tuple{Type{BufferedNonlinearFunction},JuMP._FunctionStorage,MathOptInterface.NLPBoundsPair,Dict{Int64, Vector{Int64}},Vector{JuMP._Derivatives.Linearity},OperatorRegistry,Vector{Float64},NS})   # time: 0.030106
    Base.precompile(Tuple{typeof(MathOptInterface.is_empty),Optimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt}})   # time: 0.0280181
    Base.precompile(Tuple{typeof(initialize!),RelaxCache{2, NS},DirectedTree})   # time: 0.0236066
    Base.precompile(Tuple{typeof(_propagate_constraint!),Evaluator,BufferedNonlinearFunction{1, NS}})   # time: 0.0191075
    Base.precompile(Tuple{typeof(rprop!),Relax,Evaluator,BufferedNonlinearFunction{2, NS}})   # time: 0.0188975
    isdefined(EAGO, Symbol("#157#160")) && Base.precompile(Tuple{getfield(EAGO, Symbol("#157#160")),BufferedNonlinearFunction{2, NS}})   # time: 0.0069586
    Base.precompile(Tuple{typeof(relax!),GlobalOptimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt},BufferedNonlinearFunction{1, NS},Int64,Bool})   # time: 0.005456
    Base.precompile(Tuple{typeof(MathOptInterface.add_variable),Optimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt}})   # time: 0.0048171
    Base.precompile(Tuple{typeof(MathOptInterface.add_constraint),Optimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt},MathOptInterface.SingleVariable,MathOptInterface.LessThan{Float64}})   # time: 0.0040936
    Base.precompile(Tuple{typeof(MathOptInterface.add_constraint),Optimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt},MathOptInterface.SingleVariable,MathOptInterface.GreaterThan{Float64}})   # time: 0.0038925
    Base.precompile(Tuple{typeof(MathOptInterface.set),Optimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt},MathOptInterface.RawParameter,Int64})   # time: 0.0033517
    Base.precompile(Tuple{typeof(MathOptInterface.set),Optimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt},MathOptInterface.NLPBlock,MathOptInterface.NLPBlockData})   # time: 0.0032988
    Base.precompile(Tuple{typeof(MathOptInterface.set),Optimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt},MathOptInterface.ObjectiveFunction{MathOptInterface.SingleVariable},MathOptInterface.SingleVariable})   # time: 0.0022907
    Base.precompile(Tuple{Type{RelaxCache{1, NS}}})   # time: 0.0020124
    Base.precompile(Tuple{Type{RelaxCache{2, NS}}})   # time: 0.0017719
    isdefined(EAGO, Symbol("#133#135")) && Base.precompile(Tuple{getfield(EAGO, Symbol("#133#135")),BufferedNonlinearFunction{2, NS}})   # time: 0.0014442
    Base.precompile(Tuple{typeof(MathOptInterface.set),Optimizer{Incremental{GLPK.Optimizer}, Incremental{Ipopt.Optimizer}, DefaultExt},MathOptInterface.ObjectiveSense,MathOptInterface.OptimizationSense})   # time: 0.0012987
    return
end
