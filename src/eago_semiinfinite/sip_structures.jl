"""SIP_opts an option type for the SIP routine. The fields it contains are as follows:
- tol: the absolute tolerance for convergence
- kmax: the maximum number of iterations
- eps_g0: initial eps_g restriction parameter, must be greater than zero
- r0: parameter used to reduce eps_g as each iteration (e.g. eps_g = eps_g/r0)
- verbosity: Sets verbosity of the solution routine. Either "None", "Normal", and "Full".
- return_hist: Sets the return type. True returns LBD/UBD and times at each iteration.
- LLP_Solve:
- LBP_Solve:
- UBP_Solve:
- LLP_Opt: LLP solver options.
- LBP_Opt: LBP solver options.
- UBP_Opt: UBP solver options.
- gSIPExp: Expression for semi-infinite inequality constraint
- hSIPExp: Expression for semi-infinite equality constraints
"""
mutable struct SIP_Options
  lower_disc_set::Vector{Vector{Float64}}
  upper_disc_set::Vector{Vector{Float64}}
  tolerance::Float64
  max_iterations::Int
  initial_eps_g::Float64
  initial_r::Float64
  return_hist::Bool
  header_interval::Int
  print_interval::Int
  verbosity::String
  inner_tolerance::Float64
  model
  initialize_bnd_prob
end
SIP_Options() = SIP_Options(Vector{Float64}[], Vector{Float64}[], 1E-3, 100000, 1.0,
                            2.0, false, 20, 1, "Normal", 1.0E-8, Model(), nothing)

"""
    SIP_Result
--------------------------------------------------------------------------------
Description:
Composite type for storing the resulting form SIP solution routine.
--------------------------------------------------------------------------------
Fields:
k             Int64 - Number of iterations run
UBD           Float64 - Upper bound of SIP solution (optimal value)
LBD           Float64 - Lower bound of SIP solution
feas          Bool - Feasibility of SIP
LBP_time      Float64 - Time spent solving the lower bounding problem (sec)
LLP_time      Float64 - Time spent solving the lower level problem (sec)
UBP_time      Float64 - Time spent solving the upper bounding problem (sec)
xbar          Array{Float64} - Solution point
--------------------------------------------------------------------------------
"""
mutable struct SIP_Result
    iteration_number::Int
    upper_bound::Float64
    lower_bound::Float64
    feasibility::Bool
    lower_bounding_time::Float64
    lower_level_time::Float64
    upper_bounding_time::Float64
    x_bar::Vector{Float64}
    p_bar::Vector{Float64}
end
SIP_Result() = SIP_Result(1, Inf, -Inf, true, 0.0, 0.0, 0.0, Float64[], Float64[])

"""
    SIP_Problem_Storage

Storage to pass problem information and solution routine options.
"""
mutable struct SIP_Problem_Storage
  f::Function
  gSIP::Function
  x_l::Vector{Float64}
  x_u::Vector{Float64}
  p_l::Vector{Float64}
  p_u::Vector{Float64}

  np::Int
  nx::Int
  opts::SIP_Options

  h
  hj
  y_l::Vector{Float64}
  y_u::Vector{Float64}
  ny::Int
  upper::Symbol

  xp_l::Vector{Float64}
  xp_u::Vector{Float64}
  n_xp::Int
end
