function set_xpbar(problem_storage::SIP_Problem_Storage)
  xbar = (problem_storage.x_u + problem_storage.x_l)/2.0
  pbar = (problem_storage.p_u + problem_storage.p_l)/2.0
  return xbar, pbar, problem_storage.nx, problem_storage.np
end

function explicit_llp(xbar::Vector{Float64}, sip_storage::SIP_Result, problem_storage::SIP_Problem_Storage)

  np = problem_storage.np
  pL = problem_storage.p_l
  pU = problem_storage.p_u

  model_llp = deepcopy(problem_storage.opts.model)
  if np == 1
    g(p) = problem_storage.gSIP(xbar, p)
    register(model_llp, :g, np, g, autodiff=true)
    @variable(model_llp, pL[i] <= p[i=1:np] <= pU[i])
    @NLobjective(model_llp, Min, -g(p[1]))
  else
    gmulti(p...) = problem_storage.gSIP(xbar, p)
    register(model_llp, :gmulti, np, gmulti, autodiff=true)
    @variable(model_llp, pL[i] <= p[i=1:np] <= pU[i])
    @NLobjective(model_llp, Min, -gmulti(p...))
  end


  optimize!(model_llp)

  termination_status = JuMP.termination_status(model_llp)
  result_status = JuMP.primal_status(model_llp)
  valid_result, is_feasible = is_globally_optimal(termination_status, result_status)

  sip_storage.lower_level_time += MOI.get(model_llp, MOI.SolveTime())

  if valid_result
      if is_feasible
        INNg2 = -JuMP.objective_value(model_llp)
        sip_storage.p_bar[:] = JuMP.value.(p)
      end
  else
    error("Lower level problem.")
  end

  return INNg2, is_feasible
end

# should be done
function explicit_bnd(disc_set::Vector{Vector{Float64}}, eps_g::Float64, sip_storage::SIP_Result, problem_storage::SIP_Problem_Storage, flag::Bool)
  ng = length(disc_set)
  nx = problem_storage.nx
  xL = problem_storage.x_l
  xU = problem_storage.x_u

  # create JuMP model
  if problem_storage.opts.initialize_bnd_prob == nothing
    model_bnd = deepcopy(problem_storage.opts.model)
    @variable(model_bnd, xL[i] <= x[i=1:nx] <= xU[i])
  else
    model_bnd, x = problem_storage.opts.initialize_bnd_prob()
  end

  for i in 1:ng
      gi = Symbol("g$i")
      gtemp = x -> problem_storage.gSIP(x, disc_set[i])
      g(x...) = gtemp(x)
      register(model_bnd, gi, nx, g, autodiff=true)
      func_call = Expr(:call)
      args = []
      push!(args, gi)
      for i in 1:nx
        push!(args, JuMP.VariableRef(model_bnd, MathOptInterface.VariableIndex(i)))
      end
      func_call.args = args
      ineq_call = Expr(:call)
      ineq_call.args = [:(<=); func_call; -eps_g]
      JuMP.add_NL_constraint(model_bnd, ineq_call)
  end

  f(x...) =  problem_storage.f(x)
  register(model_bnd, :f, nx, f, autodiff=true)
  if nx == 1
    if problem_storage.sense == :min
      @NLobjective(model_bnd, Min,  f(x[1]))
    else
      @NLobjective(model_bnd, Max,  f(x[1]))
    end
  else
    if problem_storage.sense == :min
      @NLobjective(model_bnd, Min, f(x...))
    else
      @NLobjective(model_bnd, Max, f(x...))
    end
  end

  optimize!(model_bnd)

  termination_status = JuMP.termination_status(model_bnd)
  result_status = JuMP.primal_status(model_bnd)
  valid_result, is_feasible = is_globally_optimal(termination_status, result_status)

  sip_storage.lower_bounding_time += MOI.get(model_bnd, MOI.SolveTime())

  if valid_result
      if is_feasible
        if flag
          sip_storage.lower_bound = JuMP.objective_value(model_bnd)
        else
          sip_storage.upper_bound = JuMP.objective_value(model_bnd)
        end
        sip_storage.x_bar[:] = JuMP.value.(x)
      end
  else
    error("Lower problem did not solve to global optimality.")
  end

  return is_feasible
end

"""
    explicit_sip_solve

Solves a semi-infinite program via the algorithm presented in Mitsos2011 using
the EAGOGlobalSolver to solve the lower bounding problem, lower level problem,
and the upper bounding problem. The options for the algorithm and the global
solvers utilized are set by manipulating a sip_options containing the options info.
Inputs:
* `f::Function`: Objective in the decision variable. Takes a single argument
                 vector that must be untyped.
* `gSIP::Function`: The semi-infinite constraint. Takes two arguments: the first
                    being a vector containing the decision variable and the
                    second being a vector containing the uncertainity
                    variables. The function must be untyped.
* `x_l::Vector{Float64}`: Lower bounds on the decision variables.
* `x_u::Vector{Float64}`: Upper bounds on the decision variables.
* `p_l::Vector{Float64}`: Lower bounds on the uncertain variables.
* `p_u::Vector{Float64}`: Upper bounds on the uncertain variables.
* `m::JuMP.Model`: A JuMP model containing the optimizer to be used.
* `opts`: Option type a keyword argument containing problem information
Returns: A SIP_result composite type containing solution information.
"""
function explicit_sip_solve(f::Function, gSIP::Function, x_l::Vector{Float64},
                            x_u::Vector{Float64}, p_l::Vector{Float64},
                            p_u::Vector{Float64}, m::JuMP.Model; opts = nothing,
                            initialize_bnd_prob = nothing, sense = :min)

  @assert length(p_l) == length(p_u)
  @assert length(x_l) == length(x_u)
  n_p = length(p_l)
  n_x = length(x_l)

  println("ran update 1")

  if opts == nothing
      opts = SIP_Options()
      opts.model = m
  end
  if initialize_bnd_prob != nothing
    opts.initialize_bnd_prob = initialize_bnd_prob
  end

  problem_storage = SIP_Problem_Storage(f, gSIP, x_l, x_u, p_l, p_u, n_p, n_x, opts,
                                     nothing, nothing, Float64[], Float64[], 0,
                                     :nothing, sense, Float64[], Float64[], 0)

  sip_sto = core_sip_routine(explicit_llp, explicit_bnd, set_xpbar, problem_storage)
  return sip_sto
end
