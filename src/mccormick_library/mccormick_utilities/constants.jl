const IntervalType = Interval{Float64}

mutable struct McCormickParamters
  env_max_int::Int
  env_tol::Float64
  mu::Int
  mu_flt::Float64
  valid_check::Bool
  valid_tol::Float64
  multivar_refine::Bool
  mv_tol::Float64
  outer_rnding::Bool
  outer_param::Float64
  reference_point::Vector{Float64}
  reference_domain::Vector{IntervalType}
  use_reference::Bool
  McCormickParamters() = new(100,
                             1E-10,
                             0,
                             0.0,
                             false,
                             1E-8,
                             false,
                             1E-15,
                             false,
                             0.0,
                             Float64[],
                             IntervalType[],
                             false)
end

const MC_param = McCormickParamters()
const IntervalConstr = interval
const Half64 = Float64(0.5)
const Two64 = Float64(2.0)
const Three64 = Float64(3.0)
const EqualityTolerance = Float64(1E-12)
const DegToRadIntv = pi_interval(Float64)/Interval(180.0)


########### number list
int_list = [Int8,UInt8,Int16,UInt16,
            Int32,UInt32,Int64,UInt64,Int128,UInt128]
float_list = [Float16,Float32,Float64]

########### differentiable functions unitary functions
CVList = [:cosh,:exp,:exp2,:exp10] ### function is convex
CCList = [:acosh,:log,:log2,:log10,:sqrt] ### function is concave
CCtoCVList = [:asin,:sinh,:atanh,:tan] ### function is concave then convex
CVtoCCList = [:atan,:acos,:tanh,:asinh] ### function is convex then concave
Template_List = union(CVList,CCList,CCtoCVList,CVtoCCList)

########### non differentiable and non-unitary functions
OtherList = [:sin,:cos,:min,:max,:abs,:step, :sign, :inv, :*, :+, :-, :/,
:promote_rule, :convert, :one, :zero, :real, :dist, :eps, :fma, :^]
