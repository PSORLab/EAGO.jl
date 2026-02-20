# Customization Guidelines

This section contains general guidelines on how the functionality of EAGO can be extended for specific use cases. Many functions in EAGO are extensible, so examples are not provided for every possible case, but some important functions are covered. If there is a use case you do not see provided here, but would like to see, please post an issue using the GitHub [issue tracker](https://github.com/PSORLab/EAGO.jl/issues).

## 1) Creating an Extension

Extensibility in EAGO is based on the [`ExtensionType`](@ref). To create customized functions, the recommended method is to create a new structure, as follows:

```julia
using EAGO

struct MyNewStruct <: EAGO.ExtensionType end
```

To let EAGO know that you would like to use this extension (and any functions you overload), when you create the JuMP model, declare your new type in the `subsolver_block` field of the [`Optimizer`](@ref) as follows:

```julia
using JuMP

factory = () -> EAGO.Optimizer(SubSolvers(; t = MyNewStruct() ))
model = Model(factory)
```

The key point to note here is that the new structure is set to [`SubSolvers.t`](@ref EAGO.SubSolvers), which is the field that holds any user-defined extension type. Now, when EAGO calls any of its functions, it will check to see if a custom function for this new extension has been created. If so, it will use the new one; if not, it will use the default version of that function.

## 2) Preprocessing

In EAGO's branch-and-bound routine, preprocessing is performed prior to solving the lower and upper problems. By default, linear and quadratic contractor methods are performed, followed by interval constraint propagation, then optimization-based bounds tightening. An important outcome of EAGO's default preprocessing is that the `_preprocess_feasibility` field of the [`GlobalOptimizer`](@ref) is set to `true`, unless the subproblem at the current node has been proven infeasible. Therefore, to bypass preprocessing for your own problem, you must, at a minimum, set this field to `true`. For example:

```julia
import EAGO: preprocess!
function EAGO.preprocess!(t::MyNewStruct, x::EAGO.GlobalOptimizer)
    x._preprocess_feasibility = true
end
```

The user-defined preprocessing step can be as simple or complex as desired, but if `_preprocess_feasibility` is not set to `true`, EAGO will assume each node is infeasible.

## 3) Lower Problem

By default, EAGO applies Kelley's cutting-plane algorithm [[1](#References)] to solve the lower-problem. This can be overloaded using the same syntax as for the other functions. Necessary changes to the [`GlobalOptimizer`](@ref) that occur within the lower problem are changing the `_lower_objective_value` and `_lower_feasibility` fields. If the `_lower_objective_value` field is not changed, branch-and-bound will not update the lower bound. If `_lower_feasibility` is not set to `true`, the node will be discarded as infeasible. A minimum functional (though not useful) lower problem extension is as follows:

```julia
import EAGO: lower_problem!
function EAGO.lower_problem!(t::MyNewStruct, x::EAGO.GlobalOptimizer)
    x._lower_objective_value = -Inf
    x._lower_feasibility = true
end
```

Any arbitrarily complex lower problem can be substituted here in place of EAGO's default method, as long as these fields are updated. Note also that, although there is a separate upper problem function, if the lower problem is being replaced by an algorithm that also calculates an upper objective value, the necessary fields to update in [`upper_problem!`](@ref EAGO.upper_problem!) can simply be updated here, and the [`upper_problem!`](@ref EAGO.upper_problem!) can be overloaded by a function that does `nothing`.

## 4) Upper Problem

By default, the upper-bounding problem is run on every node up to depth `upper_bounding_depth`, and is triggered with a probability of `0.5^(depth - upper_bounding_depth)` afterwards for continuous problems. For integer problems, this approach is used in addition to running on every node up to depth `upper_bounding_depth + cont_depth`, with another trigger of probability `0.5^(depth - upper_bounding_depth - cont_depth)`. The `upper_bounding_depth` and `cont_depth` are fields of [`EAGOParameters`](@ref) and can be changed from their default values when the JuMP model is created, or at any time afterwards. If any of these trigger conditions are met, the default EAGO upper problem runs a local NLP solve. 

The important fields to update in the [`upper_problem!`](@ref EAGO.upper_problem!) are the `_upper_objective_value` and `_upper_feasibility` fields. If `_upper_objective_value` is not updated, the upper bound in the branch-and-bound algorithm will not update and the problem will not converge. If `_upper_feasibility` is not set to true, then any changes made to `_upper_objective_value` will not be updated for each node and the problem will not converge. A minimum functional (though not useful) upper problem extension is as follows:

```julia
import EAGO: upper_problem!
function upper_problem!(t::MyNewStruct, x::EAGO.GlobalOptimizer)
    x._upper_objective_value = Inf
    x._upper_feasibility = true
end
```

This example upper problem will set the upper bound on the objective value to `Inf` for each node. Given that no useful information is provided here, we could also have set `x_upper_feasibility` to `false` (or equivalently, remove the line where we set it to `true`), and the global upper bound would never be updated.

Note that if the `_upper_objective_value` is changed elsewhere, such as in the definition for the [`lower_problem!`](@ref EAGO.lower_problem!), the `_upper_feasibility` flag must be set to `true`. If this is not done, the change to the `_upper_objective_value` will be discarded.

## 5) Convergence Check

By default, EAGO checks to see if the lower and upper bounds have converged to within either the absolute or relative tolerance. This method of checking convergence may not be desired if, for example, only the absolute tolerance is relevant, and you want to ensure that the program does not end prematurely due to a relative tolerance limit being reached. The fields to check are `_lower_objective_value` and `_upper_objective_value` for the best-identified global lower and upper bounds, respectively. This function should return `true` if the lower and upper bounds have converged, and `false` otherwise. An example of how to specify that the convergence check only use the absolute tolerance is as follows:

```julia
import EAGO: convergence_check
function EAGO.convergence_check(t::MyNewStruct, x::EAGO.GlobalOptimizer)
    gap = (x._upper_objective_value - x._lower_objective_value)
    return (gap <= x._parameters.absolute_tolerance)
end
```

## 6) Postprocessing

Postprocessing is the final step before a node is branched on. By default, EAGO performs duality-based bounds tightening [[2](#References)] up to an iteration limit set by `dbbt_depth` in the [`EAGOParameters`](@ref) field. The important field to update in postprocessing is `_postprocess_feasibility`, which must be set to `true` for EAGO to branch on any given node. The minimum working postprocessing function is therefore:

```julia
import EAGO: postprocess!
function EAGO.postprocess!(t::MyNewStruct, x::EAGO.GlobalOptimizer)
    x._postprocess_feasibility = true
end
```

If `_postprocess_feasibility` is not set to `true`, no nodes will be branched on.

## 7) Termination Check

This is the check that occurs on each iteration of the branch-and-bound algorithm that determines whether the algorithm continues or not. By default, several conditions are checked for such as the satisfaction of absolute or relative tolerances, solution infeasibility, or other specified limits. This function returns `true` if any of the stopping conditions have been met, and  branch-and-bound should stop, and `false` otherwise. The important fields to update are `_end_state`, which takes values of `EAGO.GlobalEndState`, `_termination_status_code`, which takes values of [`MOI.TerminationStatusCode`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.TerminationStatusCode), and `_result_status_code`, which takes values of [`MOI.ResultStatusCode`](https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.ResultStatusCode). Combined, these fields provide information about the branch-and-bound completion status and result feasibility. 

As an example, suppose we have updated the [`convergence_check`](@ref EAGO.convergence_check) function to only check for absolute tolerance, and based on our knowledge of the problem, this is the only condition we care about, and we know the solution will be optimal. We could then overload the [`termination_check`](@ref EAGO.termination_check) function as follows:

```julia
import EAGO: termination_check
function EAGO.termination_check(t::MyNewStruct, x::EAGO.GlobalOptimizer)
    flag = EAGO.convergence_check(t, x)
    if flag
        x._end_state = EAGO.GS_OPTIMAL
        x._termination_status_code = MathOptInterface.OPTIMAL
        x._result_status_code = MathOptInterface.FEASIBLE_POINT
    end
    return flag
end
```

## References

1. Kelley, J. E. “The Cutting-Plane Method for Solving Convex Programs.” *Journal of the Society for Industrial and Applied Mathematics*, vol. 8, no. 4, pp. 703–12 (1960). 
2. Tawarmalani, M., Sahinidis, N. V. "Global optimization of mixed-integer nonlinear programs: A theoretical and computational study." *Math. Program., Ser. A*, 99, pp. 563-591 (2004).
