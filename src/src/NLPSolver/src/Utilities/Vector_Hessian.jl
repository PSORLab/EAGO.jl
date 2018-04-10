"""
    vector_hessian(f,x,b)

Computes the hessian of the lagrangian of `f` at `x` with factors `b` via FowardDiff.jl.
"""
function vector_hessian(f, x, b)
  n1 = length(x)
  out = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(f, x), x)
  out1 =  reshape(out, n1, n1, length(b))
  c = zeros(n1,n1)
  for i = 1:length(b)
    c += b[i]*out1[:,:,i]
  end
  return c
end
