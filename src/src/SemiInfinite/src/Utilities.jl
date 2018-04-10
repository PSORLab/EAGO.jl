length(x::IntervalArithmetic.IntervalBox{N,T}) where {N,T} = N

function lo_hi(x::IntervalArithmetic.Interval)
  temp1 = x.lo
  temp2 = x.hi
  return temp1, temp2
end

function lo_hi(X::IntervalArithmetic.IntervalBox)
  if length(X)>1
    tuple_lo, tuple_hi = lo_hi.(X)
    arr_lo = collect(tuple_lo)g
    arr_hi = collect(tuple_hi)
  else
    arr_lo = [X[1].lo]
    arr_hi = [X[1].hi]
  end
  return arr_lo,arr_hi
end

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
