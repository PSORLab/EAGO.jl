# Temporary copy of resources from experimental code transformation package
# https://github.com/perrutquist/CodeTransformation.jl until it becomes tagged.


import Core: SimpleVector, svec, CodeInfo
import Base: uncompressed_ast, unwrap_unionall

jl_method_def(argdata::SimpleVector, ci::CodeInfo, mod::Module) =
    ccall(:jl_method_def, Cvoid, (SimpleVector, Any, Ptr{Module}), argdata, ci, pointer_from_objref(mod))

typevars(T::UnionAll) = (T.var, typevars(T.body)...)
typevars(T::DataType) = ()

@nospecialize # The functions below need not specialize on arguments

getmodule(F::Type{<:Function}) = F.name.mt.module
getmodule(f::Function) = getmodule(typeof(f))

makesig(f::Function, args) = Tuple{typeof(f), args...}

argdata(sig) = svec(unwrap_unionall(sig).parameters::SimpleVector, svec(typevars(sig)...))
argdata(sig, f::Function) = svec(svec(typeof(f), unwrap_unionall(sig).parameters[2:end]...), svec(typevars(sig)...))

addmethod!(f::Function, argtypes::Tuple, ci::CodeInfo) = addmethod!(makesig(f, argtypes), ci)

function addmethod!(sig::Type{<:Tuple{F, Vararg}}, ci::CodeInfo) where {F<:Function}
    jl_method_def(argdata(sig), ci, getmodule(F))
end
