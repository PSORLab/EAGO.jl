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
function sipRes_bnd(initialize_extras, disc_set::Vector{Vector{Vector{Float64}}}, eps_g::Float64, sip_storage::SIPResult,
  problem_storage::SIPProblem, flag::Bool, f, gSIP)
  nx = problem_storage.nx
  xL = problem_storage.x_l
  xU = problem_storage.x_u

  # create JuMP model
  model_bnd = deepcopy(problem_storage.opts.model)
  @variable(model_bnd, xL[i] <= x[i=1:nx] <= xU[i])
  initialize_extras(model_bnd, x)

  for i in 1:problem_storage.nSIP
    for j in 1:length(disc_set)
        gi = Symbol("g$i$j")
        gtemp = x -> gSIP[i](x, disc_set[j])
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

  if ~(valid_result && feasible)
    error("Lower problem did not solve to global optimality.")
  else
    objective_value = JuMP.objective_value(model_bnd)
    xsol = JuMP.value.(x)
  end

  return objective_value, xsol, is_feasible
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

function sipRes(init_bnd, prob::SIPProblem, result, f, gSIP)


  verbosity = prob.verbosity
  header_interval = prob.header_interval
  print_interval = prob.print_interval

  # initializes solution
  LBD = -Inf
  UBD = Inf
  k = 0

  xbar = (prob.x_u + prob.x_l)/2.0
  pbar = (prob.p_u + prob.p_l)/2.0
  nx = prob.nx
  np = prob.np
  nSIP = prob.nSIP

  lower_disc = prob.init_lower_disc
  upper_disc = prob.init_upper_disc
  llp1_obj_val = Inf
  llp2_obj_val = Inf
  feas = true
  xstar = fill(NaN,(nx,))
  new_disc_points = fill(zeros(Float64, np), )

  # checks for convergence
  for k = 1:prob.max_iterations

    # check for termination
    check_convergence(LBD, UBD, prob.abs_tolerance) && (break)

    # solve lower bounding problem and check feasibility
    val, xsol, feas = sipRes_bnd(init_bnd, lower_disc, 0.0, result, prob, true, f, gSIP)::Bool
    print_summary!(verbosity, LBD, xbar, result.feasibility, "LBD")
    if (~feas)
      result.feasibility = feas
      println("Lower Bounding Problem Not Feasible. Algorithm Terminated")
      break
    else
      LBD = val
    end

    # solve inner program  and update lower discretization set
    non_positive_flag = true
    for (i, gsip) in enumerate(gSIP)
      llp1_objval, pbar, feas = sipRes_llp(sto.x_bar, result, prob, gSIP)
      print_summary!(verbosity, llp1_obj_val, pbar, feas, "LLP$i")
      if (llp1_obj_val > 0.0)
        non_positive_flag = false
      end
      push!(lower_disc, pbar)
    end
    if non_positive_flag
      result.xsol[:] = xbar
      result.feasibility = true
      break
    end

    # solve upper bounding problem, if feasible solve lower level problem,
    # and potentially update upper discretization set
    feas = sipRes_bnd(init_bnd, upper_disc, eps_g, result, prob, false, f, gSIP)
    print_summary!(verbosity, sto.upper_bound, sto.x_bar, feas, "UBD")
    if feas
      llp2_obj_val, feas = lower_level_problem(sto.x_bar,  sto, problem_storage, gSIP)::Tuple{Float64,Bool}
      print_llp2!(problem_storage.opts, llp2_obj_val, sto.p_bar, feas)
      sto.upper_bound, xstar, eps_g = llp2_update(llp2_obj_val, r, eps_g,
                                                  xbar, pbar, UBD, UBD, upper_disc)
    else
      eps_g = eps_g/r
    end

    # print iteration information and advance
    print_int!(verbosity, header_interval, print_interval, k, sto.lower_bound, sto.upper_bound, eps_g, r)
    sto.iteration_number = k
  end

  return sto
end

no_init_bnd(m::JuMP.Model, x::Vector{JuMP.VariableRef}) = ()
function explicit_sip_solve(x_l::Vector{Float64}, x_u::Vector{Float64},
                            p_l::Vector{Float64}, p_u::Vector{Float64},
                            gSIP::Function, f::Function; kwargs...)

  @assert length(p_l) == length(p_u)
  @assert length(x_l) == length(x_u)
  ~isempty(kwargs) && allowed_options_check!(kwargs)

  # collects all keyword arguments of the form :gSIP1, gSIP24234, :gSIPHello
  kv_gSIP = filter(x -> (String(x[1])[1:4] === "gSIP"), kwargs)
  gSIP = values(kv_gSIP)

  init_bnd = haskey(kwargs, :init_bnd) ? kwargs[:init_bnd] : no_init_bnd
  algo = haskey(kwargs, :algo) ? kwargs[:algo] : :sipRes

  prob = SIPProblem(x_l, x_u, p_l, p_u, kwargs)
  result = SIPResult()

  if (algo === :sipRes)
    sip_sto = sipRes(init_bnd, prob, result, f, gSIP)
  else
    error("Desired algorithm unsupported.")
  end
  return sip_sto
end
