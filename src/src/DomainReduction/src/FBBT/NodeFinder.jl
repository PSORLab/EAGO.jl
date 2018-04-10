global NodeCounter = 0
global EdgeList = []
global ExprList = Expr[]
global HeaderList = Symbol[]
global ConstList = []

immutable NodeFinder
  ind::Int64
end

step(x::T) where T = x>0 ? one(T) : zero(T)

# Defines conversion for NodeFinder Type
for flt in (Float16,Float32,Float64,Int8,Int16,Int32,Int64,Int128,
            UInt8,UInt16,UInt32,UInt64,UInt128)
     @eval function convert(::Type{NodeFinder}, x::$flt)
              global NodeCounter += 1
              push!(ConstList,[NodeCounter,x])
              return NodeFinder(NodeCounter)
           end
     @eval promote_rule(::Type{NodeFinder}, ::Type{$flt}) = NodeFinder
end


for flt in (Int8,Int16,Int32,Int64,Int128,
            UInt8,UInt16,UInt32,UInt64,UInt128)
    # Utility function for updating
    @eval function Update_Tracker_Unary(hd::Symbol,n::$flt)
              push!(EdgeList,[[n,NodeCounter+1]])
              push!(HeaderList,hd)
              global NodeCounter += 1
              return NodeFinder(NodeCounter)
          end

    for flt1 in (Int8,Int16,Int32,Int64,Int128,
                UInt8,UInt16,UInt32,UInt64,UInt128)

        @eval function Update_Tracker_Binary(hd::Symbol,n1::$flt,n2::$flt1)
                  push!(EdgeList,[[n1,NodeCounter+1],[n2,NodeCounter+1]])
                  push!(HeaderList,hd)
                  global NodeCounter += 1
                  return NodeFinder(NodeCounter)
              end
    end
end

# Unary NodeFinder Operation Definitions
exp(x::NodeFinder) = Update_Tracker_Unary(:exp,x.ind)
exp2(x::NodeFinder) = Update_Tracker_Unary(:exp2,x.ind)
exp10(x::NodeFinder) = Update_Tracker_Unary(:exp10,x.ind)
log(x::NodeFinder) = Update_Tracker_Unary(:log,x.ind)
log2(x::NodeFinder) = Update_Tracker_Unary(:log2,x.ind)
log10(x::NodeFinder) = Update_Tracker_Unary(:log10,x.ind)
acosh(x::NodeFinder) = Update_Tracker_Unary(:acosh,x.ind)
cosh(x::NodeFinder) = Update_Tracker_Unary(:cosh,x.ind)
sqrt(x::NodeFinder) = Update_Tracker_Unary(:sqrt,x.ind)
asin(x::NodeFinder) = Update_Tracker_Unary(:asin,x.ind)
sinh(x::NodeFinder) = Update_Tracker_Unary(:sinh,x.ind)
atanh(x::NodeFinder) = Update_Tracker_Unary(:atanh,x.ind)
tan(x::NodeFinder) = Update_Tracker_Unary(:tan,x.ind)
atan(x::NodeFinder) = Update_Tracker_Unary(:atan,x.ind)
acos(x::NodeFinder) = Update_Tracker_Unary(:acos,x.ind)
tanh(x::NodeFinder) = Update_Tracker_Unary(:tanh,x.ind)
asinh(x::NodeFinder) = Update_Tracker_Unary(:asinh,x.ind)
sin(x::NodeFinder) = Update_Tracker_Unary(:sin,x.ind)
cos(x::NodeFinder) = Update_Tracker_Unary(:cos,x.ind)
abs(x::NodeFinder) = Update_Tracker_Unary(:abs,x.ind)
step(x::NodeFinder) = Update_Tracker_Unary(:step,x.ind)
sign(x::NodeFinder) = Update_Tracker_Unary(:sign,x.ind)
inv(x::NodeFinder) = Update_Tracker_Unary(:inv,x.ind)
-(x::NodeFinder) = Update_Tracker_Unary(:-,x.ind)
one(x::NodeFinder) = Update_Tracker_Unary(:one,x.ind)
zero(x::NodeFinder) = Update_Tracker_Unary(:zero,x.ind)

# Binary Nodefinder Operator Definitions
/(x::NodeFinder,y::NodeFinder) = Update_Tracker_Binary(:/,x.ind,y.ind)
-(x::NodeFinder,y::NodeFinder) = Update_Tracker_Binary(:-,x.ind,y.ind)
+(x::NodeFinder,y::NodeFinder) = Update_Tracker_Binary(:+,x.ind,y.ind)
*(x::NodeFinder,y::NodeFinder) = Update_Tracker_Binary(:*,x.ind,y.ind)
min(x::NodeFinder,y::NodeFinder) = Update_Tracker_Binary(:min,x.ind,y.ind)
max(x::NodeFinder,y::NodeFinder) = Update_Tracker_Binary(:max,x.ind,y.ind)
^(x::NodeFinder,y::NodeFinder) = Update_Tracker_Binary(:^,x.ind,y.ind)

# Defines binary interactions with Float64 & Int64
for g in (:+,:-,:+,:*,:min,:max,:^)
    for t in (Float16, Float32, Float64,
              Int8,Int16,Int32,Int64,Int128,
              UInt8,UInt16,UInt32,UInt64,UInt128)
        @eval ($g)(x::NodeFinder,y::$t) = ($g)(x::NodeFinder,convert(NodeFinder,y))
        @eval ($g)(x::$t,y::NodeFinder) = ($g)(convert(NodeFinder,x),y::NodeFinder)
    end
end
