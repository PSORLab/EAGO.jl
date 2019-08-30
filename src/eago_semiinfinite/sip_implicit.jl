function set_xpbar_implicit(problem_storage::SIP_Problem_Storage)
  xbar = (problem_storage.x_u + problem_storage.x_l)/2.0
  pbar = (problem_storage.p_u + problem_storage.p_l)/2.0
  return xbar, pbar, problem_storage.nx, problem_storage.np
end

function implicit_llp(xbar::Vector{Float64}, sip_storage::SIP_Result, problem_storage::SIP_Problem_Storage)       # SHOULD BE GOOD

  println("implicit_llp start:")
  np = problem_storage.np
  p_L = problem_storage.p_l
  p_U = problem_storage.p_u
  y_L = problem_storage.y_l
  y_U = problem_storage.y_u

  println("xbar: $xbar")
  println("problem_storage.opts.model: $(problem_storage.opts.model)")

  gLLP(y, p) = problem_storage.gSIP(xbar, y, p)
  hLLP(H, y, p) = problem_storage.h(H, xbar, y, p)
  hjLLP(J, y, p) = problem_storage.hj(J, xbar, y, p)
  var, model_llp = solve_implicit(gLLP, hLLP, y_L, y_U, p_L, p_U, problem_storage.opts.model, hjLLP, nothing)

  termination_status = MOI.get(model_llp, MOI.TerminationStatus())
  result_status = MOI.get(model_llp, MOI.PrimalStatus())
  valid_result, is_feasible = is_globally_optimal(termination_status, result_status)
  sip_storage.lower_level_time += MOI.get(model_llp, MOI.SolveTime())

  if valid_result
      if is_feasible
        INNg2 = -MOI.get(model_llp, MOI.ObjectiveValue())
        sip_storage.p_bar[:] = MOI.get(optimizer, MOI.VariablePrimal(), var)
      end
  else
    error("Lower level problem.")
  end
  println("implicit_llp end:")

  return INNg2, is_feasible
end

# should be done
function implicit_bnd(disc_set::Vector{Vector{Float64}}, eps_g::Float64, sip_storage::SIP_Result, problem_storage::SIP_Problem_Storage, flag::Bool)

  println("implicit_bnd start:")

  ng = length(disc_set)
  nx = problem_storage.nx
  np = problem_storage.np
  ny = problem_storage.ny

  x_L = problem_storage.x_l
  x_U = problem_storage.x_u
  p_L = problem_storage.p_l
  p_U = problem_storage.p_u
  y_L = problem_storage.y_l
  y_U = problem_storage.y_u

  gBND(y, x) = [problem_storage.gSIP(x, y[(i*(ny-1)+1):(i*ny)], pdisc[i]) + eps_g for i=1:ng]
  function hBND!(H, y, x)
    fill!(H, 0.0)
    for i in 1:ng
      hview = view(H, (1+ny*(i-1)):(ny*i))
      problem_storage.h(hview, x, y[(i*(ny-1)+1):(i*ny)], pdisc[i])
    end
  end
  function hjBND!(J, y, x)
    fill!(J, 0.0)
    for i in 1:ng
      hjview = view(J, (1+ny*(i-1)):(ny*i), (1+ny*(i-1)):(ny*i))
      problem_storage.hj(hjview, x, y[(i*(ny-1)+1):(i*ny)], pdisc[i])
    end
  end

  yL_BND = Float64[]; yU_BND = Float64[]
  for i in 1:ng
    push!(yL_BND, y_L); push!(yU_BND, y_U)
  end

  println("problem_storage.opts.model: $(problem_storage.opts.model)")
  if ng > 0
    var, model_bnd = solve_implicit(problem_storage.f, hBND!, yL_BND, yU_BND, p_L, p_U, problem_storage.opts.model, hjBND!, gBND)
  else
    var, model_bnd = solve_implicit(problem_storage.f, hBND!, yL_BND, yU_BND, p_L, p_U, problem_storage.opts.model, hjBND!, nothing)
  end

  termination_status = MOI.get(model_bnd, MOI.TerminationStatus())
  result_status = MOI.get(model_bnd, MOI.PrimalStatus())
  valid_result, is_feasible = is_globally_optimal(termination_status, result_status)
  sip_storage.lower_bounding_time += MOI.get(model_bnd, MOI.SolveTime())

  if valid_result
      if is_feasible
        if flag
          sip_storage.lower_bound = MOI.get(model_bnd, MOI.ObjectiveValue())
        else
          sip_storage.upper_bound = MOI.get(model_bnd, MOI.ObjectiveValue())
        end
        sip_storage.x_bar[:] = MOI.get(model_bnd, MOI.VariablePrimal(), var[(ng*ny+1):(ng*ny+np)]) # NEEDS TO BE JUST STATE VARIABLES (SHOULD BE THIS...)
      end
  else
    error("Lower problem did not solve to global optimality.")
  end
  println("implicit_bnd end:")

  return is_feasible
end

"""
    implicit_sip_solve

Solves a semi-infinite program via the algorithm presented in Mitsos2011 using
the EAGO solver to solve the lower bounding problem, lower level problem,
and the upper bounding problem. The options for the algorithm and the global
solvers utilized are set by manipulating a sip_options containing the options info.
Inputs:
* `f::Function`: Objective in the decision variable. Takes a single argument
                 vector that must be untyped.
* `gSIP::Function`: The semi-infinite constraint. Takes two arguments: the first
                    being a vector containing the decision variable and the
                    second being a vector containing the uncertainity
                    variables. The function must be untyped.
* `h::Function`: Equality constraints to be relaxed via an implicit approach
* `hj::Function`: Jacobian of equality constraints to be relaxed
* `x_l::Vector{Float64}`: Lower bounds on the decision variables.
* `x_u::Vector{Float64}`: Upper bounds on the decision variables.
* `p_l::Vector{Float64}`: Lower bounds on the uncertain variables.
* `p_u::Vector{Float64}`: Upper bounds on the uncertain variables.
* `y_l::Vector{Float64}`: Lower bounds on the state variables.
* `y_u::Vector{Float64}`: Upper bounds on the state variables.
* `opts::SIP_opts`: Option type containing problem information
Returns: A SIP_result composite type containing solution information. (SEEMS GOOD )
"""
function implicit_sip_solve(f::Function, gSIP::Function, h::Function, hj::Function,
                            x_l::Vector{Float64}, x_u::Vector{Float64},
                            p_l::Vector{Float64}, p_u::Vector{Float64},
                            y_l::Vector{Float64}, y_u::Vector{Float64},
                            m::EAGO.Optimizer; opts = nothing)

  @assert length(p_l) == length(p_u)
  @assert length(x_l) == length(x_u)
  @assert length(y_l) == length(y_u)

  n_p = length(p_l)
  n_x = length(x_l)
  n_y = length(y_l)

  if opts == nothing
      opts = SIP_Options()
      opts.model = m
  end

  f2(y,x) = f(x)

  problem_storage = SIP_Problem_Storage(f2, gSIP, x_l, x_u, p_l, p_u, n_p, n_x, opts,
                                        h, hj, y_l, y_u, n_y, :nothing, Float64[], Float64[], 0)

  sip_sto = core_sip_routine(implicit_llp, implicit_bnd, set_xpbar_implicit, problem_storage)
  return sip_sto
end
