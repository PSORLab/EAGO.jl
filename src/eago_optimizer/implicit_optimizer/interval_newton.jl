"""
     interval_newton_opt
A structure which holds options for interval newton method.

- interval_newton_opt() creates the structure with dummy values.
- interval_newton_opt(nx,h,hj) creates the structure with default values:
etol is the equality tolerance, kmax is the number of iterations, nx is the state
of the state space x, h(x,p) are the equality constraints that define x(p), is
the jacobian of h w.r.t. x only
"""
struct interval_newton_opt
  etol::Float64
  kmax::Int
  nx::Int
  h::Function
  hj::Function
end
interval_newton_opt() = interval_newton_opt(0.0, 0, 0, dummy_function, dummy_function) # TO ADD
interval_newton_opt(nx,h,hj) = interval_newton_opt(1E-8,50,nx,h,hj)

"""
     preconditioner_1!(h,X,P)
Directly applies inverse preconditioning matrix for 1-by-matrix.
"""
function preconditioner_1!(J,Y)
  Y[1] = 1.0/(mid(J[1]))
end

"""
    preconditioner_n(h,X,P)

Directly applies inverse preconditioning matrix for n-by-n matrices.
"""
function preconditioner_n!(J,Y)
  Y[:,:] = inv(mid.(J))
end


"""
    strict_x_in_y(X::Vector{IntervalType},Y::Vector{IntervalType},nx::Int)

Returns true if X is strictly in Y (X.lo>Y.lo && X.hi<Y.hi).
"""
function strict_x_in_y(X::Vector{IntervalType},Y::Vector{IntervalType},nx::Int)
  k::Bool = true
  for i=1:nx
    if ((lo(X[i])<=lo(Y[i])) || (hi(X[i])>=hi(Y[i])))
      k = false
    end
  end
  return k
end

"""
    strict_x_in_y(X::IntervalType,Y::IntervalType)

Returns true if X is strictly in Y (X.lo>Y.lo && X.hi<Y.hi).
"""
function strict_x_in_y(X::IntervalType,Y::IntervalType)
  (lo(X)<=lo(Y)) && return false
  (hi(X)>=hi(Y)) && return false
  return true
end

"""
    is_equal(X1::Vector{Interval{T}},X2::Vector{Interval{T}},atol::Float64)

Returns true if X1 and X2 are equal to within tolerance atol in all dimensions.
"""
function is_equal(X1::Vector{IntervalType},X2::Vector{IntervalType},atol::Float64,nx::Int)
  out::Bool = true
  for i=1:nx
    if (abs(lo(X1[i])-lo(X2[i]))>=atol || abs(hi(X1[i])-hi(X2[i]))>=atol)
        out = false
        break
    end
  end
  return out
end

"""
    interval_newton_nb!

Implements a non-branching form of interval newton's method for use in solving
the upper bounding problem. The method will thow an error if it divides by an
interval containing zero as these cases are assumed to be removed during preprocessing.
"""
function interval_newton_nb!(H,J,Y,B,M,x::Vector{IntervalType},p::Vector{IntervalType}, opts::interval_newton_opt)
  if opts.nx == 1
    x[:] = interval_newton_1_nb!(H,J,Y,B,M,x,p,opts)
  elseif opts.nx > 1
    x[:] = interval_newton_n_nb!(H,J,Y,B,M,x,p,opts)
  else
    #error("State variable size must be greater than zero.")
  end
end

function precondition_block_1!(H,J,Y,B,M,x_mid,X,P,opts::interval_newton_opt)
  opts.h(H,x_mid,x_mid,P,Float64[])
  opts.hj(J,X,X,P,Float64[])
  preconditioner_1!(J,Y)
  B[1] = Y[1]*H[1]
  M[1] = Y[1]*J[1]
end

function interval_newton_1_kernel_nb!(M,B,x_mid,X,N)
  if lo(M[1])*hi(M[1]) > 0.0
      N[1] = x_mid[1] - (B[1])/M[1]
  else
      error("Extended division occured in upper bounding problem.")
  end
  if ~isdisjoint(N[1],X[1])
      X[1] = N[1] ∩ X[1]
  else
    error("Interval Newton shows non-existence of x(p) for current value of p.")
  end
end

"""
    interval_newton_1_nb!

Interval newton calculation for single-valued h.
"""
function interval_newton_1_nb!(H,J,Y,B,M,Xin::Vector{IntervalType},P::Vector{IntervalType}, opts::interval_newton_opt)

  # Initializes variables
  X = IntervalType[Xin[1]]
  Xi = IntervalType[Xin[1]]
  N = IntervalType[Xin[1]]
  x_mid = mid.(Xin[1])
  k::Int = 1

  continuation_condition::Bool = true
  while continuation_condition
    Xi = IntervalType.(lo.(X),hi.(X))
    x_mid = mid.(X)
    precondition_block_1!(H,J,Y,B,M,x_mid,X,P,opts)
    interval_newton_1_kernel_nb!(M,B,x_mid,X,N)
    continuation_condition = (k<opts.kmax) && (is_equal(X,Xi,opts.etol,opts.nx) == false)
    if (~continuation_condition)
      return X
    end
    k += 1
  end
end

"""
    precondition_block_n!

Preconditions h(x_mid,P) and hj(x_mid,P) using an inverse midpoint Jacobian
for a nx-by-nx matrix.
"""
function precondition_block_n!(H,J,Y,B,M,x_mid,X,P,opts::interval_newton_opt)
  opts.h(H,x_mid,x_mid,P,Float64[])
  opts.hj(J,X,X,P,Float64[])
  preconditioner_n!(J,Y)
  B[:] = Y*H
  M[:,:] = Y*J
end

function interval_newton_n_kernel_nb!(M,B,x_mid,X,N,nx)
  for i=1:nx
    S1 = IntervalType(0.0)
    S2 = IntervalType(0.0)
    for j=1:nx
      if (j<i)
        S1 += M[i,j]*(X[j]-x_mid[j])
      elseif (j>i)
        S2 += M[i,j]*(X[j]-x_mid[j])
      end
    end
    if lo(M[i,i])*hi(M[i,i]) > 0.0
      N[i] = x_mid[i] - (B[i]+S1+S2)/M[i,i]
    else
      error("Extended division occured in upper bounding problem.")
    end
    if ~isdisjoint(N[i],X[i])
        X[i] = N[i] ∩ X[i]
    else
      error("Interval Newton shows non-existence of x(p) for current value of p.")
    end
  end
end

"""
    interval_newton_n!

Interval newton calculation for vector-valued h.
"""
function interval_newton_n_nb!(H,J,Y,B,M,x::Vector{IntervalType},P::Vector{IntervalType},opts::interval_newton_opt)

  # Initializes variables
  nx::Int = opts.nx
  S1::IntervalType = IntervalType(0.0)
  S2::IntervalType = IntervalType(0.0)

  X::Vector{IntervalType} = copy(x)
  Xi::Vector{IntervalType} = copy(x)
  N::Vector{IntervalType} = copy(x)
  x_mid::Vector{Float64} = mid.(x)
  k::Int = 1

  continuation_condition::Bool = true
  while ((k<opts.kmax) && continuation_condition)

    Xi = copy(X)
    x_mid = mid.(X)

    precondition_block_n!(H,J,Y,B,M,x_mid,X,P,opts)
    interval_newton_n_kernel_nb!(M,B,x_mid,X,N,nx)
    k += 1
    continuation_condition = (k<opts.kmax) && (is_equal(X,Xi,opts.etol,opts.nx) == false)
    if (~continuation_condition)
      return X
    end
  end
end
