"""
    print_int!(sip_options,k_int,lbd,ubd,eps,r)
--------------------------------------------------------------------------------
Description:
Prints current iteration statistics.
--------------------------------------------------------------------------------
Inputs:
sip_options   SIP_opts - Options object
k_int    Int64   - Iteration Number
lbd      Float64 - Lower Bound
ubd      Float64 - Upper Bound
eps      Float64 - Restriction size
r        Float64 - Restriction Adjustment
--------------------------------------------------------------------------------
Returns:
No returned value but prints the iteration, lower bound, upper bound,
restriction value, r value, absolute ratio
and relative ratio if the verbosity is set to "Normal".
--------------------------------------------------------------------------------
"""
@inline function print_int!(sip_options::SIP_Options, k_int::Int, lbd::Float64, ubd::Float64, eps::Float64, r::Float64)
  if (sip_options.verbosity == "Normal")
    # prints header line every hdr_intv times
    if (mod(k_int,sip_options.header_interval)==0||k_int==1)
      println("Iteration   LBD    UBD     eps      r     Absolute_Gap    Absolute_Ratio")
    end
    # prints iteration summary every prnt_intv times
    if (mod(k_int,sip_options.print_interval)==0)
      ptr_arr_temp = [k_int lbd ubd eps r (ubd-lbd) (lbd/ubd)]
      ptr_arr1 = join([Printf.@sprintf("%6u",x) for x in ptr_arr_temp[1]], ",   ")
      ptr_arr2 = join([Printf.@sprintf("%3.7f",x) for x in ptr_arr_temp[2:5]], ",     ")
      ptr_arr3 = join([Printf.@sprintf("%6u",x) for x in ptr_arr_temp[6:7]], ",")
      println(string(ptr_arr1,",      ",ptr_arr2,",      ",ptr_arr3))
    end
  end
end

function print_llp1!(sip_options::SIP_Options, INNg1::Float64, pbar::Vector{Float64}, feas::Bool)
  if (sip_options.verbosity == "Full" || sip_options.verbosity == "Normal")
    println("solved INN #1: ",INNg1," ",pbar," ",feas)
  end
end

function print_llp2!(sip_options::SIP_Options, INNg2::Float64, pbar::Vector{Float64}, feas::Bool)
  if (sip_options.verbosity == "Full" || sip_options.verbosity == "Normal")
    println("solved INN #2: ",INNg2," ",pbar," ",feas)
  end
end

function print_lbp!(sip_options::SIP_Options, LBDg::Float64, xbar::Vector{Float64}, feas::Bool)
  if (sip_options.verbosity == "Full" || sip_options.verbosity == "Normal")
    println("solved LBD: ",LBDg," ",xbar," ",feas)
  end
end

function print_ubp!(sip_options::SIP_Options, UBD_temp::Float64, xbar::Vector{Float64}, feas::Bool)
  if (sip_options.verbosity == "Full" || sip_options.verbosity == "Normal")
    println("solved UBD: ",UBD_temp," ",xbar," ",feas)
  end
end
