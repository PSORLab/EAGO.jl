length(x::IntervalArithmetic.IntervalBox{N,T}) where {N,T} = N

"""
--------------------------------------------------------------------------------
Function: print_int!
--------------------------------------------------------------------------------
Description:
Prints current iteration statistics.
--------------------------------------------------------------------------------
Inputs:
SIPopt   SIP_opts - Options object
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
@inline function print_int!(SIPopt,k_int,lbd,ubd,eps,r)
  if (SIPopt.Verbosity == "Normal")
    # prints header line every hdr_intv times
    if (mod(k_int,SIPopt.hdr_intv)==0||k_int==1)
      println("Iteration   LBD    UBD     eps      r     Absolute_Gap    Absolute_Ratio")
    end
    # prints iteration summary every prnt_intv times
    if (mod(k_int,SIPopt.prnt_intv)==0)
      ptr_arr_temp = [k_int lbd ubd eps r (ubd-lbd) (lbd/ubd)]
      ptr_arr1 = join([@sprintf("%6u",x) for x in ptr_arr_temp[1]], ",   ")
      ptr_arr2 = join([@sprintf("%3.7f",x) for x in ptr_arr_temp[2:5]], ",     ")
      ptr_arr3 = join([@sprintf("%6u",x) for x in ptr_arr_temp[6:7]], ",")
      println(string(ptr_arr1,",      ",ptr_arr2,",      ",ptr_arr3))
    end
  end
end
