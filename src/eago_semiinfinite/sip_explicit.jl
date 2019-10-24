function sipRes_llp1(xbar::Vector{Float64}, result::SIPResult,
                    problem::SIPProblem, cb::SIPCallback, indx::Int64)

  println("xbar: $xbar")

  pL = problem.p_l
  pU = problem.p_u
  np = problem.np

  model_llp1 = Model(with_optimizer(EAGO.Optimizer,
                               lp_depth = 100000,
                               lp_reptitions = 0,
                               quad_uni_depth = -1,
                               obbt_depth = 0,
                               cp_depth = -1,
                               iteration_limit = 100000000,
                               verbosity = 0,
                               output_iterations = 1,
                               header_iterations = 400000,
                               relative_tolerance = 1E-4,
                               absolute_tolerance = 1E-4,
                               dbbt_depth = 100000000,
                               subgrad_tighten = true, #true,
                               objective_cut_on = true,
                               cut_max_iterations = 3,
                               upper_bounding_depth = 4,
                               time_limit = 1000.0,
                               presolve_epigraph_flag = false,
                               presolve_cse_flag = false,
                               presolve_flatten_flag = false,
                               relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))
  #model_llp1 = deepcopy(problem.model)
  @variable(model_llp1, pL[i] <= p1[i=1:np] <= pU[i])
  if np == 1
    g(p...) = cb.gSIP[indx](xbar, p)
    register(model_llp1, :g, np, g, autodiff=true)
    @NLobjective(model_llp1, Min, -g(p1[1]))
  else
    gmulti(p...) = cb.gSIP[indx](xbar, p)
    register(model_llp1, :gmulti, np, gmulti, autodiff=true)
    @NLobjective(model_llp1, Min, -gmulti(p1...))
  end

  optimize!(model_llp1)
  termination_status = JuMP.termination_status(model_llp1)
  result_status = JuMP.primal_status(model_llp1)
  valid_result, is_feasible = is_globally_optimal(termination_status, result_status)
  (~valid_result) && error("Error encountered in lower level problem.")

  obj_val = -JuMP.objective_value(model_llp1)
  if np == 1
    p_bar = JuMP.value.(p1)
  else
    p_bar = JuMP.value.(p1)
  end
  result.solution_time += MOI.get(model_llp1, MOI.SolveTime())

  return obj_val, p_bar, is_feasible
end

function sipRes_llp2(xbar::Vector{Float64}, result::SIPResult,
                    problem::SIPProblem, cb::SIPCallback, indx::Int64)

  pL = problem.p_l
  pU = problem.p_u
  np = problem.np

  println("xbar: $xbar")

  model_llp2 = Model(with_optimizer(EAGO.Optimizer,
                               lp_depth = 100000,
                               lp_reptitions = 0,
                               quad_uni_depth = -1,
                               obbt_depth = 0,
                               cp_depth = -1,
                               iteration_limit = 100000000,
                               verbosity = 0,
                               output_iterations = 1,
                               header_iterations = 400000,
                               relative_tolerance = 1E-4,
                               absolute_tolerance = 1E-4,
                               dbbt_depth = 100000000,
                               subgrad_tighten = true, #true,
                               objective_cut_on = true,
                               cut_max_iterations = 3,
                               upper_bounding_depth = 4,
                               time_limit = 1000.0,
                               presolve_epigraph_flag = false,
                               presolve_cse_flag = false,
                               presolve_flatten_flag = false,
                               relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))
  #model_llp2 = deepcopy(problem.model)
  @variable(model_llp2, pL[i] <= p2[i=1:np] <= pU[i])
  if np == 1
    g(p) = cb.gSIP[indx](xbar, p)
    register(model_llp2, :g, np, g, autodiff=true)
    @NLobjective(model_llp2, Min, -g(p2[1]))
  else
    gmulti(p...) = cb.gSIP[indx](xbar, p)
    register(model_llp2, :gmulti, np, gmulti, autodiff=true)
    @NLobjective(model_llp2, Min, -gmulti(p2...))
  end

  optimize!(model_llp2)
  termination_status = JuMP.termination_status(model_llp2)
  result_status = JuMP.primal_status(model_llp2)
  valid_result, is_feasible = is_globally_optimal(termination_status, result_status)
  (~valid_result) && error("Error encountered in lower level problem.")

  obj_val = -JuMP.objective_value(model_llp2)
  if np == 1
    p_bar = JuMP.value.(p2)
  else
    p_bar = JuMP.value.(p2)
  end
  result.solution_time += MOI.get(model_llp2, MOI.SolveTime())

  return obj_val, p_bar, is_feasible
end

# should be done
function sipRes_bnd(initialize_extras, disc_set::Vector{Vector{Vector{Float64}}},
                    eps_g::Float64, sip_storage::SIPResult,
                    problem_storage::SIPProblem, flag::Bool, cb::SIPCallback)

  nx = problem_storage.nx
  xL = problem_storage.x_l
  xU = problem_storage.x_u

  # create JuMP model
  model_bnd = Model(with_optimizer(EAGO.Optimizer,
                               lp_depth = 100000,
                               lp_reptitions = 0,
                               quad_uni_depth = -1,
                               obbt_depth = 0,
                               cp_depth = -1,
                               iteration_limit = 100000000,
                               verbosity = 0,
                               output_iterations = 1,
                               header_iterations = 400000,
                               relative_tolerance = 1E-4,
                               absolute_tolerance = 1E-4,
                               dbbt_depth = 100000000,
                               subgrad_tighten = true, #true,
                               objective_cut_on = true,
                               cut_max_iterations = 3,
                               upper_bounding_depth = 4,
                               time_limit = 1000.0,
                               presolve_epigraph_flag = false,
                               presolve_cse_flag = false,
                               presolve_flatten_flag = false,
                               relaxed_optimizer = Gurobi.Optimizer(OutputFlag=0)))
  #model_bnd = deepcopy(problem_storage.model)
  @variable(model_bnd, xL[i] <= x[i=1:nx] <= xU[i])
  initialize_extras(model_bnd, x)
  println("length(disc_set): $(length(disc_set))")
  println("problem_storage.nSIP: $(problem_storage.nSIP)")
  for i in 1:problem_storage.nSIP
    for j in 1:length(disc_set)
        gi = Symbol("g$i$j")
        g(x...) = cb.gSIP[i](x, disc_set[j][i])
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

  obj(x...) =  cb.f(x)
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
  out_time = MOI.get(model_bnd, MOI.SolveTime())

  if ~(valid_result && is_feasible)
    error("Lower problem did not solve to global optimality.")
  else
    objective_value = JuMP.objective_value(model_bnd)
    xsol = JuMP.value.(x)
    println("xsol: $(xsol)")
  end

  return objective_value, xsol, is_feasible, out_time
end

const DEFINED_KWARGS_SIPRES = Symbol[:initialize_extras, :initialize_bnd_prob,
                                     :sense, :algo, :m]
function allowed_options_check!(kwargs::Base.Iterators.Pairs)
  flag = true
  for (k,v) in kwargs
    if ~in(k, DEFINED_KWARGS_SIPRES)
        flag = false
        println("Allowed keyword arguements are $DEFINED_KWARGS_SIPRES and any additional
                 SIP constraints of the form gSIPX...X.")
        break
    end
  end
  flag
end

function sipRes(init_bnd, prob::SIPProblem, result::SIPResult, cb::SIPCallback)


  verbosity = prob.verbosity
  header_interval = prob.header_interval
  print_interval = prob.print_interval
  abs_tolerance = prob.absolute_tolerance

  # initializes solution
  LBD = -Inf
  UBD = Inf
  k = 0

  xbar = (prob.x_u + prob.x_l)/2.0
  pbar = (prob.p_u + prob.p_l)/2.0
  nx = prob.nx
  np = prob.np
  nSIP = prob.nSIP
  tolerance = 1E-6

  eps_g =  prob.initial_eps_g
  r = prob.initial_r
  lower_disc = prob.init_lower_disc
  upper_disc = prob.init_upper_disc
  feas = true
  xstar = fill(NaN,(nx,))
  new_disc_points = fill(zeros(Float64, np), )

  # checks for convergence
  for k = 1:prob.iteration_limit

    println("lower_disc: $(lower_disc)")
    println("upper_disc: $(upper_disc)")

    # check for termination
    check_convergence(LBD, UBD, abs_tolerance, verbosity) && (break)

    # solve lower bounding problem and check feasibility
    println("---- LOWER BOUND START----")
    val, xsol, feas = sipRes_bnd(init_bnd, lower_disc, 0.0, result, prob, true, cb)
    if (~feas)
      result.feasibility = feas
      println("Lower Bounding Problem Not Feasible. Algorithm Terminated")
      break
    else
      result.lower_bound = val
      LBD = val
    end
    print_summary!(verbosity, LBD, xsol, result.feasibility, "LBD")
    println("---- LOWER BOUND END----")

    # solve inner program  and update lower discretization set
    println("---- LOWER LEVEL PROBLEM START 1 ----")
    non_positive_flag = true
    temp_lower_disc = Vector{Float64}[]
    for i in 1:nSIP
      llp_out = sipRes_llp1(xsol, result, prob, cb, i)
      print_summary!(verbosity, llp_out[1], llp_out[2], llp_out[3], "LLP$i")
      if (llp_out[1] + tolerance > 0.0)
        non_positive_flag = false
      end
      push!(temp_lower_disc, llp_out[2])
    end
    push!(lower_disc, temp_lower_disc)
    if non_positive_flag
      result.xsol[:] = xsol
      result.feasibility = true
      break
    end
    println("---- LOWER LEVEL PROBLEM END 1 ----")

    # solve upper bounding problem, if feasible solve lower level problem,
    # and potentially update upper discretization set
    println("---- UPPER BOUND START ----")
    val, xsol, feas = sipRes_bnd(init_bnd, upper_disc, eps_g, result, prob, false, cb)
    print_summary!(verbosity, val, xsol, feas, "UBD")
    println("---- UPPER BOUND END ----")
    if feas
      println("---- LOWER LEVEL PROBLEM START 2 ----")
      non_positive_flag = true
      temp_upper_disc = Vector{Float64}[]
      for i in 1:nSIP
        llp2_out = sipRes_llp2(xsol, result, prob, cb, i)
        print_summary!(verbosity, llp2_out[1], llp2_out[2], llp2_out[3], "LLP$i")
        if (llp2_out[1] + tolerance > 0.0)
          non_positive_flag = false
        end
        println("pbar: $pbar")
        push!(temp_upper_disc, llp2_out[2])
      end
      if non_positive_flag
          if (val <= UBD)
              UBD = val
              result.upper_bound = UBD
              xstar[:] = xsol[:]
          end
          eps_g = eps_g/r
          println("UPDATED eps_g 1: $eps_g")
      else
          push!(upper_disc, temp_upper_disc)
          println("UPDATED upper_disc")
      end
      println("---- LOWER LEVEL PROBLEM END 2 ----")
    else
      eps_g = eps_g/r
      println("UPDATED eps_g 2: $eps_g")
    end
    print_summary!(verbosity, result.upper_bound, result.xsol, feas, "UBD")

    # print iteration information and advance
    print_int!(verbosity, header_interval, print_interval, k, result.lower_bound, result.upper_bound, eps_g, r)
    result.iteration_number = k
  end

  return result
end

no_init_bnd(m::JuMP.Model, x::Vector{JuMP.VariableRef}) = ()

function explicit_sip_solve(x_l::Vector{Float64}, x_u::Vector{Float64},
                            p_l::Vector{Float64}, p_u::Vector{Float64},
                            f::Function, gSIP; kwargs...)

  @assert length(p_l) == length(p_u)
  @assert length(x_l) == length(x_u)
  ~isempty(kwargs) && allowed_options_check!(kwargs)

  # collects all keyword arguments of the form :gSIP1, gSIP24234, :gSIPHello
  init_bnd = haskey(kwargs, :init_bnd) ? kwargs[:init_bnd] : no_init_bnd
  algo = haskey(kwargs, :algo) ? kwargs[:algo] : :sipRes
  m = haskey(kwargs, :m) ? kwargs[:m] : EAGO.Optimizer()

  prob = SIPProblem(x_l, x_u, p_l, p_u, gSIP, m, kwargs)
  result = SIPResult()
  result.xsol = fill(0.0, (prob.nx,))
  result.psol = fill(0.0, (prob.np,))

  cb = SIPCallback(f, gSIP)

  if (algo === :sipRes)
    sip_sto = sipRes(init_bnd, prob, result, cb)
  else
    error("Desired algorithm unsupported.")
  end
  return sip_sto
end

function explicit_sip_solve(x_l::Vector{Float64}, x_u::Vector{Float64},
                            p_l::Vector{Float64}, p_u::Vector{Float64},
                            f::Function, gSIP::Function, m)
    explicit_sip_solve(x_l, x_u, p_l, p_u, f, [gSIP], m = m)
end
