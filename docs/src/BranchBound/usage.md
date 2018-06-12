## Example 1 - Setup and Solve a Basic Problem
In the below example, we solve for minima of $f(x)=x_1+x_2^2$ on the domain $x_1 \in [-1,1]$, $x_2 \in [2,9]$. Natural interval extensions are used to compute the upper and lower bounds. The natural interval extensions are provided by the Validated Numerics package.

First, we create a BnBModel object which contains all the relevant problem info and a BnBSolver object that contains all nodes and their associated values. We specify default conditions for the Branch and Bound problem. Default conditions are a best-first search, relative width bisection, normal verbosity, a maximum of 1E6 nodes, an absolute tolerance of 1E-6, and a relative tolerance of 1E-3.

```julia

using EAGO
using ValidatedNumerics
b = [Interval(-1,1),Interval(1,9)]
a = BnBModel(b)
c = BnBSolver()
EAGO.set_to_default!(c)
c.BnB_atol = 1E-4

```

Next, the lower and upper bounding problems are defined. These problems must return a tuple containing the upper/lower value, a point corresponding the upper/lower value, and the feasibility of the problem. We then set the lower/upper problem of the BnBModel object and solve the BnBModel & BnBSolver pair.

```julia

function ex_LBP(X,k,pos,opt,temp)
  ex_LBP_int = @interval X[1]+X[2]^2
  return ex_LBP_int.lo, mid.(X), true, []
end
function ex_UBP(X,k,pos,opt,temp)
  ex_UBP_int = @interval X[1]+X[2]^2
  return ex_UBP_int.hi, mid.(X), true, []
end

c.Lower_Prob = ex_LBP
c.Upper_Prob = ex_UBP

outy = solveBnB!(c,a)

```
The solution is then returned in b.soln and b.UBDg is it's value. The corresponding output displayed to the console is given below.

![BnBChart2](BnBChart2.png)

## Example 2 - Adjust Solver Tolerances
The absolute tolerance can be adjusted as shown below

```julia

julia> a.BnB_tol = 1E-4

```

The relative tolerance can be changed in a similar manner

```julia

julia> a.BnB_rtol = 1E-3

```

## Example 3 - Select Alternative Search Routines
In the above problem, the search routine could be set to a breadth-first or depth-first routine by using the set_Branch_Scheme command

```julia

julia> EAGO.set_Branch_Scheme!(a,"breadth")
julia> EAGO.set_Branch_Scheme!(a,"depth")

```
## Example 4 - Adjust Information Printed
In order to print, node information in addition to iteration information the verbosity of the BnB routine can be set to full as shown below

```julia

julia> EAGO.set_Verbosity!(a,"Full")

```
Similarly, if one wishes to suppress all command line outputs, the verbosity can be set to none.

```julia

julia> EAGO.set_Verbosity!(a,"None")

```
