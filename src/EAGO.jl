# Main module to import subroutines

module EAGO

# call mathprogbase inteface
using MathProgBase
importall MathProgBase.SolverInterface

# call ValidatedNumerics interval libraries
using IntervalArithmetic
using IntervalContractors
import IntervalArithmetic: mid, isdisjoint, isinterior, mig, mag, precedes,
                           strictprecedes, sup, inf, hull, diam, pow, dist,
                           isthin, isentire, radius, widen, sqr
import IntervalRootFinding.bisect
import IntervalContractors: plus_rev, minus_rev, mul_rev, power_rev, sqr_rev,
                            sqrt_rev, abs_rev, sin_rev, cos_rev, tan_rev,
                            asin_rev, log_rev, exp_rev, inv_rev

# call external autodiff and solver packages
using ForwardDiff
using ReverseDiff
using Ipopt
using Clp
using StaticArrays
using BandedMatrices

# imports base functions to extend
import Base:
    +, -, *, /, //, fma,
    <, >, ==, !=, ⊆, ^, <=, >=,
    in, zero, one, abs, abs2, real, min, max,
    sqrt, exp, log, sin, cos, tan, inv,
    exp2, exp10, log2, log10,
    asin, acos, atan, atan2, step,
    sinh, cosh, tanh, asinh, acosh, atanh,
    union, intersect, isempty,
    convert, promote_rule, eltype,
    BigFloat, float, widen, big,
    ∩, ∪, ⊆, eps, length, middle,
    floor, ceil, trunc, sign, round,
    expm1, log1p, isfinite, isnan, show, showall,
    isinteger, setdiff, parse

# attempts to call commerical solvers
try
       include("/home/mewilhel/.julia/v0.6/Snopt/src/Snopt.jl")
       using Snopt
end

include("src/Interval/EAGOIntervalArithmetic.jl")
include("src/ParametricInterval/EAGOParametricInterval.jl")
include("src/BranchBound/EAGOBranchBound.jl")
include("src/McCormick/EAGOSmoothMcCormickGrad.jl")
include("src/DomainReduction/EAGODomainReduction.jl")
include("src/NLPSolver/EAGONLPSolver.jl")
include("src/SemiInfinite/EAGOSemiInfinite.jl")

end # module
