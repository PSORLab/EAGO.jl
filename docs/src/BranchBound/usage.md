# Optimization with user-defined subroutines

[Matthew Wilhelm](https://psor.uconn.edu/person/matthew-wilhelm/)  
Department of Chemical and Biomolecular Engineering, University of Connecticut

### Overview  
In this section, we construct an optimizer that uses EAGO's basic nlp solution routine with user-defined lower and upper bounding problems. The **EAGO.Optimizer** structure supplies a number of parameters and stored structures that advanced users may find useful for constructing specialized solution routines. For a full review, of the EAGO.optimizer object the reader is directed to the **EAGO.Optimizer** docstring and documentation provided at [LINKTOWEBSITE]().

In this example, we'll forgo extensive integration into the EAGO.optimizer and simply replace the lower and upper-bounding problems to construct B&B routine that solves the following problem to global optimality using bounds obtained from interval arithmetic:

$
\begin{align}
&\min_{\mathbf x \in X} \;\; \sin(x_1)x_2^2 - \cos(x_3) / x_4 \label{eq1} \\
&X = [-10,10]\times[-1,1]\times[-10,10]\times[2,20].
\end{align}
$

We begin importing EAGO, ValidatedNumerics (for interval arithmetic), and JuMP.

```julia
using EAGO, IntervalArithmetic, JuMP
```

### Defining a custom lower bounding problem

A valid lower bound is obtained from the lower bound of the natural interval extension using the **ValidatedNumerics.jl** package. The LowerProblem used accepts the **EAGO.Optimizer** structure and a **EAGO.NodeData** structure, computes the bound by method overloading interval arithmetic, and stores the results to the appropriate field of the **EAGO.Optimizer's** which is a mutable structure of type **EAGO.LowerInfo**. Note that the problem is unconstrained on the domain so we can assume it is always feasible.

```julia
function LowerProblem!(opt::EAGO.Optimizer, n::EAGO.NodeBB)
    x = Interval.(n.lower_variable_bounds,n.upper_variable_bounds)
    FInterval = sin(x[1])x[2]^2-cos(x[3])/x[4]
    opt.current_lower_info.value = FInterval.lo
    opt.current_lower_info.solution = mid.(x)
    opt.current_lower_info.feasibility = true
end
```

### Defining a custom upper bounding problem

```julia
function UpperProblem!(opt::EAGO.Optimizer, n::EAGO.NodeBB)
    x_value = (n.upper_variable_bounds-n.lower_variable_bounds)/2.0
    x_interval = Interval.(x_value)
    FInterval = sin(x_interval[1])x_interval[2]^2-cos(x_interval[3])/x_interval[4]
    opt.current_upper_info.value = FInterval.hi
    opt.current_upper_info.solution = x_value
    opt.current_upper_info.feasibility = true
end
```

### Build the JuMP Model and optimize

We now add our optimizer to a JuMP model, provide variable bounds, user-defined functions, and optimize. Note that options can be provided to the EAGO optimizer using a series of keywords of a Dict{Symbol,Any} object. Both manners of providing options to EAGO are illustrated below.


```julia
# Creates a JuMP model with the the lower_problem, upper_problem, and absolute tolerance set by keyword arguments
m = JuMP.Model(with_optimizer(EAGO.Optimizer, lower_problem! = LowerProblem!, upper_problem! = UpperProblem!,
                              absolute_tolerance = 0.001, obbt_depth = 0, dbbt_depth = 0, cp_depth = 0,
                              treat_as_nonlinear = Int[1; 2; 3; 4]))

# Create the same model m using an options dictionary
opt_dict = Dict{Symbol, Any}()
opt_dict[:lower_problem!] = LowerProblem!
opt_dict[:upper_problem!] = UpperProblem!
opt_dict[:absolute_tolerance] = 0.001
opt_dict[:obbt_depth] = 0
opt_dict[:dbbt_depth] = 0
opt_dict[:cp_depth] = 0
opt_dict[:treat_as_nonlinear] = [1; 2; 3; 4]

m = JuMP.Model(with_optimizer(EAGO.Optimizer; opt_dict...))

# Adds variables and bounds
x_L = [-10, -1, -10, 2]
x_U = [10, 1, 10, 20]
@variable(m, x_L[i] <= x[i=1:4] <= x_U[i])

# Solves the problem
JuMP.optimize!(m)
```


### Get information from the JuMP Model object

The objective value, solution, termination status,  and primal status can then be accessed via the standard JuMP interface.


```julia
fval = JuMP.objective_value(m)
xsol = JuMP.value.(x)
status_term = JuMP.termination_status(m)
status_prim = JuMP.primal_status(m)

println("EAGO terminated with a status of $status_term and a result code of $status_prim")
println("The optimal value is: $fval, the solution found is $xsol.")
```

    EAGO terminated with a status of OPTIMAL and a result code of FEASIBLE_POINT
    The optimal value is: -1.4749615097161377, the solution found is [0.624924, 0.0605696, 0.634765, 0.545127].


### Advice for more advanced constructions

The *default_lower_problem* and *default_upper_problem* should be used templates for error handling and retreiving information from MOI models. Similarly, the other default routine are good starting points for building custom modifications.

Essentially all of EAGO's subroutines are stored to a field in the **EAGO.Optimizer** structure can be reset as user-defined functions.

![BnBChart2](BnBChart2.png)

# Example 2 - Adjust Solver Tolerances
The absolute tolerance can be adjusted as shown below

```julia
julia>  EAGO.Optimizer(abolute_tolerance = 0.001)
```

The relative tolerance can be changed in a similar manner

```julia
julia>  EAGO.Optimizer(relative_tolerance = 0.0001)
```

# Example 3 - Adjust Information Printed
In order to print, node information in addition to iteration information the verbosity
of the BnB routine can be set to full as shown below using keyword arguments

```julia
julia> EAGO.Optimizer(verbosity = 3)
```
Similarly, if one wishes to suppress all command line outputs:

```julia
julia> EAGO.Optimizer(verbosity = 0)
```

# Returning the solver to default settings.

```@docs
    EAGO.set_to_default!
```
