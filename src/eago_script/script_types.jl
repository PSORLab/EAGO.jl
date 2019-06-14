  struct NodeData
      nodetype::JuMP._Derivatives.NodeType
      index::Int64
      children::Vector{Int64}
  end
  children(x::NodeData) = x.children
  # convert method assumes that the children connectivity using in Tracer has
  # has been converted to an area hold a single parent value as is used in JuMP
  function convert(::Type{JuMP.NodeData}, x::Tracer.NodeData)
      JuMP.NodeData(x.nodetype,x.index,x.children[1])
  end
  # val = 1..n corresponds to 1..n variable
  struct SetTrace <: Real
      val::Int
  end
  SetTrace() = SetTrace(0)

  val(x::SetTrace) = x.val

  struct SetTraceSto
      storage::Vector{SetTrace}
  end

  getindex(A::SetTraceSto,i::Int64) = getindex(A.storage,i)



  export Tape

  # JuMP convention is to store child from function... function call = -1,
  # next highest call is has parent 1 and so on... This is a forward tape recording
  # architecture, so we'll store children then invert...
  mutable struct Tape
      nd::Vector{Tracer.NodeData}
      const_values::Vector{Float64}
      num_valued::Dict{Int,Bool}
      set_trace_count::Int
      const_count::Int
  end

  Tape() = Tape(Tracer.NodeData[],Float64[],Dict{Int,Bool}(),0,0)

  function Tape(n::Int)
      node_list = [Tracer.NodeData(JuMP.VARIABLE,i,[-1]) for i in 1:n]
      num_valued = Dict{Int,Bool}()
      for i in 1:n
          num_valued[i] = false
      end
      Tape(node_list,Float64[],num_valued,n,0)
  end
