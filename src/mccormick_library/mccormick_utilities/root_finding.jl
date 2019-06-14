"""
    newton(x0::T,xL::T,xU::T,f::Function,df::Function,envp1::T,envp2::T)

Defines a local 1D newton method to solve for the root of `f` between the bounds
`xL` and `xU` using `x0` as a starting point. The derivative of `f` is `df`. The
inputs `envp1` and `envp2` are the envelope calculation parameters.
"""
function newton(x0::Float64, xL::Float64, xU::Float64, f::Function, df::Function, envp1::Float64, envp2::Float64)
  dfk::Float64 = 0.0

  xk::Float64 = max(xL,min(x0,xU))
  fk::Float64 = f(xk,envp1,envp2)

  for i=1:MC_param.env_max_int
    dfk = df(xk, envp1, envp2)
    if (abs(fk) < MC_param.env_tol)
      return (xk,false)
    end
    (dfk == 0.0) && return (0.0,true)
    if (xk == xL && fk/dfk > 0.0)
      return (xk,false)
    elseif (xk == xU && fk/dfk < 0.0)
      return (xk,false)
    end
    xk = max(xL,min(xU,xk-fk/dfk))
    fk = f(xk,envp1,envp2)
  end
  (0.0,true)
end


"""
    secant(x0::T,x1::T,xL::T,xU::T,f::Function,envp1::T,envp2::T)  where {T<:Real}

Defines a local 1D secant method to solve for the root of `f` between
the bounds `xL` and `xU` using `x0` and `x1` as a starting points. The inputs
`envp1` and `envp2` are the envelope calculation parameters.
"""
function secant(x0::Float64, x1::Float64, xL::Float64, xU::Float64, f::Function, envp1::Float64, envp2::Float64)
  xkm::Float64 = max(xL,min(xU,x0))
  xk::Float64 = max(xL,min(xU,x1))
  fkm::Float64 = f(xkm,envp1,envp2)

  for i=1:MC_param.env_max_int
    fk = f(xk,envp1,envp2)
    Bk::Float64 = (fk-fkm)/(xk-xkm)
    if (abs(fk)<MC_param.env_tol)
      return (xk,false)
    end
    (Bk == 0.0) && return (0.0, true)
    if ((xk == xL) && (fk/Bk > 0.0))
      return (xk,false)
    elseif ((xk == xU) && (fk/Bk < 0.0))
      return (xk,false)
    end
    xkm = xk
    fkm = fk
    xk = max(xL,min(xU,xk-fk/Bk))
  end
  (0.0,true)
end


"""
    golden_section(xL::T,xU::T,f::Function,envp1::T,envp2::T) where {T<:Real}

Defines a local 1D golden section method to solve for the root of `f` between
the bounds `xL` and `xU` using `x0` as a starting point. Define iteration used
in golden section method. The inputs `envp1` and `envp2` are the envelope
calculation parameters.
"""
function golden_section(xL::Float64,xU::Float64,f::Function,envp1::Float64,envp2::Float64)
  fL::Float64 = f(xL,envp1,envp2)
  fU::Float64 = f(xU,envp1,envp2)

  (fL*fU > 0.0) && error("GOLDEN EXCEPTION")
  xm::Float64 = xU-(2.0-golden)*(xU-xL)
  fm::Float64 = f(xm,envp1,envp2)
  return golden_section_it(1,xL,fL,xm,fm,xU,fU,f,envp1,envp2)
end
"""
    golden_section_it(init::Int64,a::T,fa::T,b::T,fb::T,c::T,
                      fc::T,f::Function,envp1::T,envp2::T) where {T<:Real}

Define iteration used in golden section method. The inputs `fa`,`fb`, and `fc`,
are the function `f` evaluated at `a`,`b`, and `c` respectively. The inputs
`envp1` and `envp2` are the envelope calculation parameters. The value `init` is
the iteration number of the golden section method.
"""
function golden_section_it(init::Int,a::Float64,fa::Float64,b::Float64,fb::Float64,c::Float64,
                                   fc::Float64,f::Function,envp1::Float64,envp2::Float64)
  b_t_x::Bool = (c-b > b-a)
  if (b_t_x)
    x::Float64 = b + (2.0-golden)*(c-b)
  else
    x = b - (2.0-golden)*(b-a)
  end
  itr::Int = init
  if (abs(c-a)<MC_param.env_tol*(abs(b)+abs(x))||(itr>MC_param.env_max_int))
    return (c+a)/2.0
  end
  itr += 1
  fx::Float64 = f(x,envp1,envp2)
  if (b_t_x)
    if (fa*fx < 0.0)
      golden_section_it(itr,a,fa,b,fb,x,fx,f,envp1,envp2)
    else
      golden_section_it(itr,b,fb,x,fx,c,fc,f,envp1,envp2)
    end
  else
    if (fa*fb<(fa))
      golden_section_it(itr,a,fa,x,fx,b,fb,f,envp1,envp2)
    else
      golden_section_it(itr,x,fx,b,fb,c,fc,f,envp1,envp2)
    end
  end
end
