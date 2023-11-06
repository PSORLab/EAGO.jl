# Standard-Use Example 2

This example is also provided [here as a Jupyter Notebook](https://github.com/PSORLab/EAGO-notebooks/blob/master/notebooks/nlpopt_interval_bnb.ipynb).

### Using EAGO's Basic Optimizer With User-Defined Subroutines

In this section, we construct an optimizer that uses EAGO's basic NLP solution routine with user-defined lower and upper-bounding problems. The [`Optimizer`](@ref) structure supplies a number of parameters and stored structures that advanced users may find useful for constructing specialized solution routines.

In this example, we'll forgo extensive integration into the [`Optimizer`](@ref) and simply replace the lower and upper-bounding problems to construct a B&B routine that solves the following problem to global optimality using bounds obtained from interval arithmetic:

```math
\begin{aligned}
& \min_{\mathbf x \in X} \; \; \sin(x_{1}) x_{2}^{2} - \cos(x_{3}) / x_{4} \\
& X = [-10, 10] \times [-1, 1] \times [-10, 10] \times [2, 20].
\end{aligned}
```

We begin by importing EAGO, IntervalArithmetic [[1](#References)], and JuMP [[2](#References)].

```julia
using EAGO, IntervalArithmetic, JuMP
```

We now define the `IntervalExt` struct as a subtype of the [`ExtensionType`](@ref).

```julia
struct IntervalExt <: EAGO.ExtensionType end
```

## Define a Custom Lower-Bounding Problem

A valid lower bound is obtained from the lower bound of the natural interval extension using the [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) [[1](#References)] package. The `LowerProblem` is dispatched using the new `IntervalExt` structure and the [`GlobalOptimizer`](@ref) structure, computes the bound using interval arithmetic, and stores the results to the appropriate field of the [`GlobalOptimizer`](@ref). Note that the problem is unconstrained on the domain so we can assume it is always feasible. Further, since the interval bound is constrained along the entire domain associated with a node, no additional cuts will be beneficial and thus we've disabled them using the `_cut_add_flag` field.

```julia
import EAGO: lower_problem!
function lower_problem!(t::IntervalExt, x::EAGO.GlobalOptimizer)
    # Retrieve bounds at current node
    n = x._current_node
    lower = n.lower_variable_bounds
    upper = n.upper_variable_bounds
    
    # Define X for the node and compute the interval extension
    x_value = Interval.(lower, upper)
    F = sin(x_value[1])*x_value[2]^2 - cos(x_value[3])/x_value[4]
    x._lower_objective_value = F.lo
    x._lower_solution = IntervalArithmetic.mid.(x_value)
    x._lower_feasibility = true
    x._cut_add_flag = false
    
    return
end
```

## Define a Custom Upper-Bounding Problem

Since the problem is unconstrained, any feasible point represents a valid upper bound. Thus, if we arbitrarily evaluate the function at the midpoint, we obtain a valid upper bound. This function constructs an upper bound in this manner then stores the results to the appropriate field of the [`GlobalOptimizer`](@ref).

```julia
import EAGO.upper_problem!
function EAGO.upper_problem!(t::IntervalExt, x::EAGO.GlobalOptimizer)
    # Retrieve bounds at current node
    n = x._current_node
    lower = n.lower_variable_bounds
    upper = n.upper_variable_bounds
    
    # Compute midpoint value and evaluate at that point
    x_value = 0.5*(upper + lower)
    f_val = sin(x_value[1])*x_value[2]^2-cos(x_value[3])/x_value[4]
    x._upper_objective_value = f_val
    x._upper_solution = x_value
    x._upper_feasibility = true
    
    return
end
```

## Disable Unnecessary Routines

It is entirely possible to disable domain reduction by manipulating keyword arguments supplied to the optimizer. However, for simplicity's sake we'll simply overload the default preprocessing and postprocessing methods and indicate that there are no conditions under which EAGO should cut the node.

```julia
import EAGO: preprocess!, postprocess!, cut_condition
function EAGO.preprocess!(t::IntervalExt, x::EAGO.GlobalOptimizer)
    x._preprocess_feasibility = true
    return
end
function EAGO.postprocess!(t::IntervalExt, x::EAGO.GlobalOptimizer)
    x._postprocess_feasibility = true
    return
end
EAGO.cut_condition(t::IntervalExt, x::EAGO.GlobalOptimizer) = false
```

## Construct the JuMP Model and Optimize

We now add our optimizer to a JuMP [[2](#References)] model, provide variable bounds, and optimize.

```julia
# Create a factory that specifies the interval extension in EAGO's SubSolver
factory = () -> EAGO.Optimizer(SubSolvers(; t = IntervalExt()))

# Create a JuMP model using the factory, and with the absolute tolerance set by keyword argument
m = Model(optimizer_with_attributes(factory,
                    "absolute_tolerance" => 0.001
        ))

# Add variables, bounds, and the objective function
x_L = [-10.0, -1.0, -10.0, 2.0]
x_U = [10.0, 1.0, 10.0, 20.0]
@variable(m, x_L[i] <= x[i=1:4] <= x_U[i])
@NLobjective(m, Min, sin(x[1])*x[2]^2 - cos(x[3])/x[4])

# Perform the optimization
optimize!(m)
```

## Retrieve Results

The objective value, solution, termination status, and primal status can then be accessed via the standard JuMP interface.

```julia
fval = JuMP.objective_value(m)
xsol = JuMP.value.(x)
status_term = JuMP.termination_status(m)
status_prim = JuMP.primal_status(m)

println("EAGO terminated with a status of $status_term and a result code of $status_prim")
println("The optimal value is: $(round(fval, digits=3)), the solution found is $(round.(xsol, digits=4)).")
```

## Advice for More Advanced Constructions

The default [`lower_problem!`](@ref EAGO.lower_problem!) and [`upper_problem!`](@ref EAGO.upper_problem!) should be used as templates for error handling and retrieving information from MOI models.

Essentially all of EAGO's subroutines stored to a field in the [`Optimizer`](@ref) structure can be reset as user-defined functions.

## References

1. IntervalArithmetic.jl \[Computer software\] (2019). Retrieved from https://github.com/JuliaIntervals/IntervalArithmetic.jl
2. Iain Dunning and Joey Huchette and Miles Lubin. JuMP: A Modeling Language for Mathematical Optimization, SIAM Review, SIAM 59 (2017), pp. 295-320.
