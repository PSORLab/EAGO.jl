# Quick Start

EAGO is a global optimizer primarily meant to be used with the JuMP algebraic modeling
language. Typical use will involve installing EAGO and JuMP, creating a problem using
JuMP syntax, and passing the problem to the EAGO optimizer. 

## Customization

EAGO is designed to be easily extensible. Some of the examples that follow include use
cases where the standard EAGO functionality is overloaded and readily incorporated into
the main optimization routine. Information on how to extend the main branch-and-bound
functions (including lower and upper bounding routines) can be found in the
[Customization Guidelines](https://psorlab.github.io/EAGO.jl/dev/quick_start/guidelines/) section.

## Examples

The following pages in this section include several representative examples of how EAGO
can be used. Additional (and in some cases, shortened) examples can be found in the
[EAGO-notebooks repository](https://github.com/PSORLab/EAGO-notebooks/blob/master/notebooks).
Examples and instructional pages in this section include:
- [Standard-Use Example 1](@ref): A base-case optimization problem solved using the EAGO optimizer. No extensions or function overloading required.
- [Standard-Use Example 2](@ref): user-defined functions (TODO)
- [Advanced-Use Example 1](@ref): A quasiconvex optimization problem solved by overloading some of EAGO's functionality to implement a bisection-based algorithm instead of typical branch-and-bound. (TODO, but see the [Jupyter Notebook version](https://github.com/PSORLab/EAGO-notebooks/blob/master/notebooks/custom_quasiconvex.ipynb))
- [Advanced-Use Example 2](@ref): Overloading the branch-and-bound algorithm with a custom extension type.
