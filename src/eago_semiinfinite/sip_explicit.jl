function set_xpbar(problem_storage::SIPProblem)
  xbar = (problem_storage.x_u + problem_storage.x_l)/2.0
  pbar = (problem_storage.p_u + problem_storage.p_l)/2.0
  return xbar, pbar, problem_storage.nx, problem_storage.np
end

function sipRes_llp(xbar::Vector{Float64}, model, np::Int, pL::Vector{Float64},
                    pU::Vector{Float64}, gSIP)

  model_llp = deepcopy(model)
  if np == 1
    g(p) = p -> gSIP(xbar, p)
    register(model_llp, :g, np, g, autodiff=true)
    @variable(model_llp, pL[i] <= p[i=1:np] <= pU[i])
    @NLobjective(model_llp, Min, -g(p[1]))
  else
    gmulti(p...) = p -> gSIP(xbar, p)
    register(model_llp, :gmulti, np, gmulti, autodiff=true)
    @variable(model_llp, pL[i] <= p[i=1:np] <= pU[i])
    @NLobjective(model_llp, Min, -gmulti(p...))
  end

  optimize!(model_llp)
  termination_status = JuMP.termination_status(model_llp)
  result_status = JuMP.primal_status(model_llp)
  valid_result, is_feasible = is_globally_optimal(termination_status, result_status)
  (~valid_result || ~is_feasible) && error("Error encountered in lower level problem.")

  obj_val = -JuMP.objective_value(model_llp)
  p_bar = JuMP.value.(p)
  time += MOI.get(model_llp, MOI.SolveTime())

  return obj_val, is_feasible, p_bar, time
end

# should be done
function sipRes_bnd(disc_set::Vector{Vector{Float64}}, eps_g::Float64, sip_storage::SIPResult, problem_storage::SIPProblem, flag::Bool, f, gSIP)
  ng = length(disc_set)
  nx = problem_storage.nx
  xL = problem_storage.x_l
  xU = problem_storage.x_u

  # create JuMP model
  if problem_storage.opts.initialize_bnd_prob === nothing
    model_bnd = deepcopy(problem_storage.opts.model)
    @variable(model_bnd, xL[i] <= x[i=1:nx] <= xU[i])
  else
    model_bnd, x = problem_storage.opts.initialize_bnd_prob()
  end

  if problem_storage.opts.initialize_extras !== nothing
    problem_storage.opts.initialize_extras(model_bnd, x)
  end

  for i in 1:ng
      gi = Symbol("g$i")
      gtemp = x -> gSIP(x, disc_set[i])
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

  obj(x...) =  f(x)
  register(model_bnd, :obj, nx, obj, autodiff=true)
  if nx == 1
    if problem_storage.sense == :min
      @NLobjective(model_bnd, Min, obj(x[1]))
    else
      @NLobjective(model_bnd, Max, obj(x[1]))
    end
  else
    if problem_storage.sense == :min
      @NLobjective(model_bnd, Min, obj(x...))
    else
      @NLobjective(model_bnd, Max, obj(x...))
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

const DEFINED_KWARGS_SIPRES = Symbol[:initialize_extras, :initialize_bnd_prob,
                                     :sense, :algo, :m]
function allowed_options_check!(kwarg::Dict{Symbol,Any})
  flag = true
  for (k,v) in kwargs
    ~in(k, DEFINED_KWARGS_SIPRES)
     if (String(x[1])[1:4] === "gSIP")
       flag = false
       println("Allowed keyword arguements are $DEFINED_KWARGS_SIPRES and any additional
                SIP constraints of the form gSIPX...X.")
       break
     end
  end
  flag
end

function explicit_sip_solve(x_l::Vector{Float64}, x_u::Vector{Float64},
                            p_l::Vector{Float64}, p_u::Vector{Float64},
                            gSIP::Function, f::Function; kwargs...)

  @assert length(p_l) == length(p_u)
  @assert length(x_l) == length(x_u)
  allowed_options_check!(kwargs)

  # collects all keyword arguments of the form :gSIP1, gSIP24234, :gSIPHello
  kv_gSIP = filter(x -> (String(x[1])[1:4] === "gSIP"), kwargs)
  gSIP = values(kv_gSIP)

  haskey(kwargs, :init_bnd) && (initialize_bnd_prob = kwargs[:init_bnd])
  haskey(kwargs, :init_llp) && (initialize_extras = kwargs[:init_llp])

  prob = SIPProblem(x_l, x_u, p_l, p_u, kwargs)
  result = SIPResult()

  if (algo === :sipRes)
    sip_sto = SIPres(prob, result, f, gSIP)
  else
    error("Desired algorithm unsupported.")
  end
  return sip_sto
end


function sipRes(prob::SIPProblem, result, f, gSIP)


  verbosity = problem_storage.verbosity
  header_interval = problem_storage.header_interval
  print_interval = problem_storage.print_interval

  # initializes solution
  UBDg = Inf;
  LBDg = -Inf;
  k = 0

  result.xbar = (prob.x_u + prob.x_l)/2.0
  result.pbar = (prob.p_u + prob.p_l)/2.0
  nx = prob.nx
  np = prob.np

  lower_disc = prob.init_lower_disc
  upper_disc = prob.init_upper_disc
  llp1_obj_val = Inf;
  llp2_obj_val = Inf;
  feas = true
  xstar = fill(NaN,(nx,))

  # checks for convergence
  for k = 1:problem_storage.opts.max_iterations

    # check for termination
    check_convergence(sto.lower_bound, sto.upper_bound, problem_storage.opts.tolerance) && (break)

    # solve lower bounding problem and check feasibility
    feas = sipRes_bnd(lower_disc_set, 0.0, sto, problem_storage, true, f, gSIP)::Bool
    print_summary!(verbosity, sto.lower_bound, sto.x_bar, feas, "LBD")
    sto.feasibility = check_lbp_feasible(feas)
    ~sto.feasibility && (break)

    # solve inner program  and update lower discretization set
    for (i, gsip) in enumerate(gSIP)
      llp1_objval, feas = sipRes_llp(sto.x_bar, sto, problem_storage, gsip)
      print_summary!(verbosity, llp1_obj_val, sto.p_bar, feas, "LLP$i")
      xstar, UBDg, return_flag = llp1_update(llp1_objval,
                                            sto.x_bar, sto.p_bar, sto.lower_bound,
                                            sto.upper_bound, lower_disc)
      return_flag && (return sto.lower_bound, sto.upper_bound, sto)
    end

    # solve upper bounding problem, if feasible solve lower level problem,
    # and potentially update upper discretization set
    feas = sipRes_bnd(upper_disc_set, eps_g, sto, problem_storage, false, f, gSIP)
    print_summary!(verbosity, sto.upper_bound, sto.x_bar, feas, "UBD")
    if feas
      llp2_obj_val, feas = lower_level_problem(sto.x_bar,  sto, problem_storage, gSIP)::Tuple{Float64,Bool}
      print_llp2!(problem_storage.opts, llp2_obj_val, sto.p_bar, feas)
      sto.upper_bound, xstar, eps_g = llp2_update(llp2_obj_val, r, eps_g,
                                                  sto.x_bar, sto.p_bar, sto.upper_bound,
                                                  UBDg, upper_disc)
    else
      eps_g = eps_g/r
    end

    # print iteration information and advance
    print_int!(verbosity, header_interval, print_interval, k, sto.lower_bound, sto.upper_bound, eps_g, r)
    sto.iteration_number = k
  end

  return sto
end
