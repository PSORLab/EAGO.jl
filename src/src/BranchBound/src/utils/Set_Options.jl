"""
    set_Branch_Scheme!(x::BnBSolver,BM::String)

Sets the search scheme to "best", "breadth", or "depth" first schemes.
"""
function set_Branch_Scheme!(x::BnBSolver,BM::String)
  if (BM == "best")
    x.Node_Select = NS_best
    x.Branch_Sto = BM_depth_best!
  elseif (BM == "breadth")
    x.Node_Select = NS_depth_breadth
    x.Branch_Sto = BM_breadth!
  elseif (BM == "depth")
    x.Node_Select = NS_depth_breadth
    x.Branch_Sto = BM_depth_best!
  else
    error("Invalid branching method")
  end
end

"""
    set_Bisect_Func!(x::BnBSolver,BF::String,nx::Int64)

Sets the bisection function to BF = "relative midpoint" or
BF = "absolute midpoint" and disregards the first nx components of the interval
box storage.
"""
function set_Bisect_Func!(x::BnBSolver,BF::String,nx::Int64)
  if (nx<0)
    if (BF == "relative midpoint")
      x.Bisect_Func = Bisect_Rel
    elseif (BF == "absolute midpoint")
      x.Bisect_Func = Bisect_Abs
    else
      error("Invalid bisection method")
    end
  else
    if (BF == "relative midpoint")
      x.Bisect_Func = (S::BnBSolver,B::BnBModel,N) -> Bisect_Rel_Imp(S,B,N,nx)
    elseif (BF == "absolute midpoint")
      x.Bisect_Func = (S::BnBSolver,B::BnBModel,N) -> Bisect_Abs_Imp(S,B,N,nx)
    else
      error("Invalid bisection method")
    end
  end
end

"""
    set_Verbosity!(x::BnBSolver,VB::String)

Sets the verbosity (console output) to either "None", "Normal", or "Full".
"""
function set_Verbosity!(x::BnBSolver,VB::String)
  if (VB == "None"||VB == "Normal"||VB == "Full")
    x.Verbosity = VB
  else
    error("Invalid Verbosity Option")
  end
end

"""
    set_to_default!(x::BnBSolver)

Returns the B&B solver to the default settings (does not include problems or processing routines.
"""
function set_to_default!(x::BnBSolver)
  x.Term_Check = Term_Check
  x.Branch_Sto = BM_depth_best!
  x.Node_Select = NS_best
  x.Bisect_Func = Bisect_Abs
  x.Verbosity = "Normal"
  x.max_iter = Inf
  x.iter_lim = false
  x.max_nodes = 1E6
  x.BnB_atol = 1E-4
  x.BnB_rtol = 1E-4
  x.itr_intv = 1
  x.hdr_intv = 20
  x.converged = Conv_Check
  x.BnB_digits = 3
  x.hist_return = false
  x.opt = []
  x.exhaust = false
end
