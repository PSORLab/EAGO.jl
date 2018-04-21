
# Sets up reverse dictionary
R_Dict = Dict{Symbol,Symbol}()
R_Dict[:+] = :plus_rev
R_Dict[:-] = :minus_rev
R_Dict[:*] = :mul_revDR
R_Dict[:/] = :div_rev
R_Dict[:^] = :power_rev
R_Dict[:pow] = :power_rev
R_Dict[:exp] = :exp_revDR
for f in (  :sqrt, :sqr, :abs,
            :sin, :cos, :tan,
            :asin, :log,
            :acos, :atan, :sinh,
            :cosh, :tanh, :asinh,
            :acosh, :atanh, :min,
            :max, :inv)
    R_Dict[f] = Symbol(f, "_rev")
end

"""
    Tape

A storage object to record expressions in a tape. The fields are:
* `Edge_List::Vector`: Contains a list of edges corresponding the directed graph
* `Head_List::Vector{Symbol}`: Contains a list of the expressions used
* `Last_Node::Int64`: The last node in directed graph
* `Intv_Storage::Vector{Interval}`: A tape of interval values corresponding to each node.
* `FW_Arg::Vector{Array{Int64}}`: Vector point to argument node of expressions
* `FW_Expr::Vector`: List of expressions in forward propagation
* `nvar::Float64`: Number of variables in expression
* `RV_Arg::Vector{Vector{Int64}}`: Arguement list for reverse propagation
* `RV_Expr::Vector`: List of expressions in reverse propagation
* `gL::Float64`: Lower bound on the expression gL<g(x)
* `gU::Float64`: Upper bound on the expression gU>g(x)
* `cnsts::Vector`: List of constant values and their location
"""
immutable Tape{V}
  Edge_List::Vector
  Head_List::Vector{Symbol}
  Last_Node::Int64
  Intv_Storage::Vector{V}
  FW_Arg::Vector{Array{Int64}}
  FW_Expr::Vector
  nvar::Float64
  RV_Arg::Vector{Vector{Int64}}
  RV_Expr::Vector
  gL::Float64
  gU::Float64
  cnsts::Vector
  fix_val
end

"""
    Tape()

Initializes an empty tape for storage.
"""
Tape(V) =  ([],
           Symbol[],
           0,
           V[],
           [[1]],
           [],
           0.0,
           [[1]],
           [],
           0.0,
           0.0,
           [],
           [])
Tape() = Tape(MCInterval{Float64})
"""
    TapeList

A storage object for a list of tapes. Has a single field sto which is an array
of Tape objects.
"""
immutable TapeList{V}
  sto::Any
end

"""
    TapeList

Initializes an empty tapelist object
"""
TapeList(::Type{V}) where {V} = TapeList{V}([Tape(V)])
TapeList() = TapeList{MCInterval{Float64}}([Tape()])

"""
    Generate_Tape(exp::Expr,nx::Int64,gL,gU)

Generates the tape of the provided expression `exp::Expr` assuming the expression
depends on `nx::Int64` variables and has lower bounds `gL` and upper bounds `gU`.
The variables in the expression must be off the form `x[1],...,x[nx]`.
"""
function Generate_Tape(exp::Expr,nx::Q,gL,gU,V) where {Q<:Integer}

  gL,gU = Float64(gL),Float64(gU)

  # Sets up node finding object and generates graph to global
  X_NodeF::Vector{NodeFinder} = NodeFinder[NodeFinder(i) for i=1:nx]
  global NodeCounter = deepcopy(nx)
  fX! = @eval x::Vector{NodeFinder} -> $exp
  Base.invokelatest(fX!,X_NodeF)

  # generates forward tape
  FW_Arg = [[EdgeList[i][j][1] for j=1:length(EdgeList[i])] for i=1:length(EdgeList)]
  RW_Arg = [vcat([EdgeList[i][1][2]],[EdgeList[i][j][1] for j=1:length(EdgeList[i])]) for i=1:length(EdgeList)]
  tape = Tape{V}(deepcopy(EdgeList),deepcopy(HeaderList),
                      NodeCounter,V[V(-Inf,Inf) for i=1:NodeCounter],
                      FW_Arg,[Expr(:call) for i=1:NodeCounter],nx,
                      RW_Arg,[Expr(:call) for i=1:NodeCounter],gL,gU,
                      deepcopy(ConstList),[])
  global EdgeList = []
  global HeaderList = Symbol[]
  global ConstList = []
  return tape
end

"""
    Generate_Fixed_Tape(exp::Expr,nx::Int64,gL,gU,vals)

Generates the tape of the provided expression `exp::Expr` assuming the expression
depends on `nx::Int64` variables and has lower bounds `gL` and upper bounds `gU`.
The variables in the expression must be off the form `x[1],...,x[nx]`. Variables `x[nx+1]`
to `x[end]` are fixed to the values in the `vals` array.
"""
function Generate_Fixed_Tape(exp::Expr,nx::Q,gL,gU,vals,V) where {Q<:Integer}

  gL,gU = Float64(gL),Float64(gU)

  # Sets up node finding object and generates graph to global
  X_NodeF::Vector{Any} = vcat(Any[NodeFinder(i) for i=1:nx],vals)
  global NodeCounter = deepcopy(nx)
  fX! = @eval x::Vector{Any} -> $exp
  Base.invokelatest(fX!,X_NodeF)

  # generates forward tape
  FW_Arg = [[EdgeList[i][j][1] for j=1:length(EdgeList[i])] for i=1:length(EdgeList)]
  RW_Arg = [vcat([EdgeList[i][1][2]],[EdgeList[i][j][1] for j=1:length(EdgeList[i])]) for i=1:length(EdgeList)]
  tape = Tape{V}(deepcopy(EdgeList),deepcopy(HeaderList),
                      NodeCounter,V[V(-Inf,Inf) for i=1:NodeCounter],
                      FW_Arg,[Expr(:call) for i=1:NodeCounter],nx,
                      RW_Arg,[Expr(:call) for i=1:NodeCounter],gL,gU,
                      deepcopy(ConstList),vals)
  global EdgeList = []
  global HeaderList = Symbol[]
  global ConstList = []
  return tape
end

"""
    Generate_TapeList(exprs::Vector{Expr},nx::Int64,gL::Vector{Float64},gU::Vector{Float64})

Generates the tape list for each provided expression `exprs[i]` in
`exprs::Vector{Expr}` assuming all expressions depend on `nx::Int64` variables
and has lower bounds `gL[i]` and upper bounds `gU[i]`. The variables in the
expression must be off the form `x[1],...,x[nx]`.
"""
function Generate_TapeList(exprs::Vector{Expr},nx::Q,gL::Vector{Float64},gU::Vector{Float64},V) where {Q<:Integer}
  @assert length(exprs) == length(gL) == length(gU)
  tapelist = []
  for i=1:length(exprs)
    push!(tapelist,Generate_Tape(exprs[i],nx,gL[i],gU[i]),V)
  end
  return TapeList(tapelist)
end

"""
    Generate_Fixed_TapeList(exprs::Vector{Expr},nx::Int64,gL::Vector{Float64},gU::Vector{Float64},val_arr)

Generates the tape list for each provided expression `exprs[i]` in
`exprs::Vector{Expr}` assuming all expressions depend on `nx::Int64` variables
and has lower bounds `gL[i]` and upper bounds `gU[i]`. The variables in the
expression must be off the form `x[1],...,x[nx]`. Variables `x[nx+1]`
to `x[end]` are fixed to the values in the `vals` array.
"""
function Generate_Fixed_TapeList(exprs::Vector{Expr},nx::Q,gL::Vector{Float64},
                                 gU::Vector{Float64},val_arr,V) where {Q<:Integer}
  @assert length(exprs) == length(gL) == length(gU)
  tapelist = []
  for i=1:length(exprs)
      for j=1:length(val_arr)
          push!(tapelist,Generate_Fixed_Tape(exprs[i],nx,gL[i],gU[i],val_arr[j],V))
      end
  end
  return TapeList{V}(tapelist)
end

"""
    SetConstraintNode!(x::Tape)

Sets the terminal interval to `[x.gL,x.gU]` using the bounds in the `x::Tape`.
"""
function SetConstraintNode!(x::Tape{V}) where {V}
  x.Intv_Storage[x.Last_Node] = V(x.gL,x.gU)
end
"""
    SetConstraintNode!(x::TapeList)

Sets the terminal interval to `[x.gL[i],x.gU[i]]`  in each tape using the bounds
in the `x::TapeList[i]`.
"""
function SetConstraintNode!(x::TapeList{V}) where {V}
  for i=1:length(x.sto)
    SetConstraintNode!(x.sto[i])
  end
end

"""
    SetConstantNode!(x::Tape)

Sets nodes in constaint value list of tape according to their value.
"""
function SetConstantNode!(x::Tape{V}) where {V}
  for i=1:length(x.cnsts)
    x.Intv_Storage[Int64(x.cnsts[i][1])] = V(x.cnsts[i][2])
  end
end
"""
    SetConstantNode!(x::TapeList)

Sets nodes in constaint value list of tape according to their value for
each tape in the tape list.
"""
function SetConstantNode!(x::TapeList{V}) where {V}
  for i=1:length(x.sto)
    SetConstantNode!(x.sto[i])
  end
end

"""
    SetVarBounds!(x::Tape,X::Vector{Interval})

Sets variable nodes to the interval value in `X::Vector{Interval}` in tape.
"""
function SetVarBounds!(x::Tape{V},X::Vector{V}) where {V}
  x.Intv_Storage[1:length(X)] = X
end

"""
    SetVarBounds!(x::TapeList,X::Vector{Interval})

Sets variable nodes to the interval value in `X::Vector{Interval}` for all
elements of the TapeList.
"""
function SetVarBounds!(x::TapeList{V},X::Vector{V}) where {V}
  for i=1:length(x)
    SetVarBounds!(x.sto[i],X)
  end
end

"""
    GetVarBounds(x::Tape)

Unpacks the variable node values and returns an array of intervals.
"""
function GetVarBounds(x::Tape)
  return x.Intv_Storage[1:x.Var_Num]
end

"""
    ForwardPass!(x::Tape)

Performs a forward contractor pass using tape `x::Tape`.
"""
function ForwardPass!(x::Tape)
  for i=1:length(x.Edge_List)
    x.FW_Expr[i].args = vcat([x.Head_List[i]],x.Intv_Storage[x.FW_Arg[i]])
    x.Intv_Storage[x.Edge_List[i][1][2]] = eval(x.FW_Expr[i]) ∩ x.Intv_Storage[x.Edge_List[i][1][2]]
  end
end

"""
    ForwardPassOne!(x::Tape)

Performs a forward contractor pass using tape `x::Tape` with some uninitialized nodes.
"""
function ForwardPassOne!(x::Tape)
  for i=1:length(x.Edge_List)
    x.FW_Expr[i].args = vcat([x.Head_List[i]],x.Intv_Storage[x.FW_Arg[i]])
    x.Intv_Storage[x.Edge_List[i][1][2]] = eval(x.FW_Expr[i])
  end
end

"""
    ReversePass!(x::Tape)

Performs a reverse contractor pass using tape `x::Tape`.
"""
function ReversePass!(x::Tape)
  for i=length(x.Edge_List):-1:1
    x.RV_Expr[i].args = vcat([R_Dict[x.Head_List[i]]],x.Intv_Storage[x.RV_Arg[i]])
    x.Intv_Storage[x.RV_Arg[i]] = [eval(x.RV_Expr[i])...]
  end
end

"""
    DAGContractor!(X::Vector{Interval{T}},x::TapeList,r)

Performs a forward-interval contactor propagation `r` times using the `x::TapeList`
and the initial interval bounds `X::Vector{Interval{T}}`.
"""
function DAGContractor!(X::Vector{V},x::TapeList{V},r) where {V}
  Xprev::Vector{V} = copy(X) # sets variable bounds on first Array to Box Bounds
  SetConstraintNode!(x)
  for i=1:r
    for j=1:length(x.sto)
      SetVarBounds!(x.sto[j],Xprev) # take refined bounds from previous graph
      SetConstantNode!(x)
      if (i == 1)
        ForwardPassOne!(x.sto[j])
      else
        ForwardPass!(x.sto[j]) # run forward and reverse pass
      end
      ReversePass!(x.sto[j])
      Xprev = x.sto[j].Intv_Storage[1:length(X)]
    end
  end
  X[:] = Xprev
end

"""
    getDAG(expr::Expr,root_list::Vector{Int64})

Generates the directed graph of expression `expr` using the root_list corresponding
to the indices of `x[i]` that appear in the expression.
"""
function getDAG(expr::Expr,root_list::Vector{Q},V) where {Q<:Integer}
  # generates tape
  tape = Generate_Tape(expr,maximum(root_list),0.0,0.0,V)

  # generates vertix list
  not_root = setdiff([i for i=1:Int64(tape.nvar)],root_list)
  vertix_list = filter!(x->~in(x,not_root),[i for i=1:(tape.Last_Node)])
  vertix_name = Array{Any}(tape.Last_Node)

  # generates edge list
  # to do: record edge number while taping and preallocate
  edge_list = []
  edge_temp = tape.Edge_List
  for i=1:length(edge_temp)
    for j=1:length(edge_temp[i])
      push!(edge_list,edge_temp[i][j])
      if (j == 1)
        vertix_name[edge_temp[i][j][2]] = tape.Head_List[i]
      end
    end
  end

  for i=1:length(tape.cnsts)
    vertix_name[Int64(tape.cnsts[i][1])] = tape.cnsts[i][2]
  end
  for i=1:Int64(tape.nvar)
    vertix_name[i] = :(x[$i])
  end
  deleteat!(vertix_name,not_root)

  return vertix_list,vertix_name,edge_list
end
