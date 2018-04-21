"""
    newton(x0::T,xL::T,xU::T,f::Function,df::Function,envp1::T,envp2::T)

Defines a local 1D newton method to solve for the root of `f` between the bounds
`xL` and `xU` using `x0` as a starting point. The derivative of `f` is `df`. The
inputs `envp1` and `envp2` are the envelope calculation parameters.
"""
@inline function newton(x0::T,xL::T,xU::T,f::Function,df::Function,envp1::T,envp2::T) where {T<:Real}
  dfk::T = zero(x0)

  xk::T = max(xL,min(x0,xU))
  fk::T = f(xk,envp1,envp2)

  for i=1:MC_param.env_max_int
    dfk = df(xk,envp1,envp2)
    if (abs(fk)<MC_param.env_tol)
      return xk
    end
    if (dfk == zero(x0))
      error("NEWTON EXCEPTION")
    elseif (xk==xL && fk/dfk>zero(x0))
      return xk
    elseif (xk==xU && fk/dfk<zero(x0))
      return xk
    end
    xk = max(xL,min(xU,xk-fk/dfk))
    fk = f(xk,envp1,envp2)
  end
  error("NEWTON EXCEPTION")
end


"""
    secant(x0::T,x1::T,xL::T,xU::T,f::Function,envp1::T,envp2::T)  where {T<:Real}

Defines a local 1D secant method to solve for the root of `f` between
the bounds `xL` and `xU` using `x0` and `x1` as a starting points. The inputs
`envp1` and `envp2` are the envelope calculation parameters.
"""
@inline function secant(x0::T,x1::T,xL::T,xU::T,f::Function,envp1::T,envp2::T) where {T<:Real}
  xkm::T = max(xL,min(xU,x0))
  xk::T = max(xL,min(xU,x1))
  fkm::T = f(xkm,envp1,envp2)

  for i=1:MC_param.env_max_int
    fk = f(xk,envp1,envp2)
    Bk::T = (fk-fkm)/(xk-xkm)
    if (abs(fk)<MC_param.env_tol)
      return xk
    end
    if (Bk == zero(x0))
      error("SECANT EXCEPTION")
    elseif ((xk==xL)&(fk/Bk>zero(x0)))
      return xk
    elseif ((xk==xU)&(fk/Bkzero(x0)))
      return xk
    end
    xkm = xk
    fkm = fk
    xk = max(xL,min(xU,xk-fk/Bk))
  end
  error("SECANT EXCEPTION")
end


"""
    golden_section(xL::T,xU::T,f::Function,envp1::T,envp2::T) where {T<:Real}

Defines a local 1D golden section method to solve for the root of `f` between
the bounds `xL` and `xU` using `x0` as a starting point. Define iteration used
in golden section method. The inputs `envp1` and `envp2` are the envelope
calculation parameters.
"""
@inline function golden_section(xL::T,xU::T,f::Function,envp1::T,envp2::T) where {T<:Real}
  fL::T = f(xL,envp1,envp2)
  fU::T = f(xU,envp1,envp2)

  if (fL*fU > zero(fL))
    error("GOLDEN EXCEPTION")
  end
  xm::T = xU-(2.0-golden)*(xU-xL)
  fm::T = f(xm,envp1,envp2)
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
@inline function golden_section_it(init::Q,a::T,fa::T,b::T,fb::T,c::T,
                                   fc::T,f::Function,envp1::T,envp2::T) where {Q<:Integer,T<:Real}
  b_t_x::Bool = (c-b > b-a)
  if (b_t_x)
    x::T = b + (2.0-golden)*(c-b)
  else
    x = b - (2.0-golden)*(b-a)
  end
  itr::Int64 = init
  if (abs(c-a)<MC_param.env_tol*(abs(b)+abs(x))||(itr>MC_param.env_max_int))
    return (c+a)/2.0
  end
  itr += 1
  fx::T = f(x,envp1,envp2)
  if (b_t_x)
    if (fa*fx<zero(fa))
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
