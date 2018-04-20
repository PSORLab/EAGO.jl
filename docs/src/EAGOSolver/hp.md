## Configuring EAGO for High-Performance

### Solver Parameters

Parameter considerations for explicit optimization:
- **Validation**: The validated interval arithmetic option comes with a
                  significant performance decrease but can be useful for some
                  problems.
- **Range Reduction**: Recommend selecting an arbitrary high depth for range reduction for constrained problems.
- **Probing**: Recommend limiting this to the first few nodes due to the high computational cost.
- **DBBT**: Recommend selecting an arbitrary high depth for duality-based tightening.
- **Interval Constraint Propagation**: Recommended using for problems with highly nonlinear and complex constraints.

Parameter considerations specifically for implicit optimization:
- **Interval Contractor**: Run roughly 5-10 interval iterations per McCormick contractor iteration. Recommend starting with 10.
- **McCormick Contractor**: Limit iterations to three or fewer.

### Lower Bounding Problem

Currently, two (McCormick-based) lower bounding problems are available for use
in both the explicit and implicit formulations: an LP solve and a solve via
SNOPT. Using either problem type is recommended and the correct choice will
depend on the specific problem being solve.

### LP Solver Selection

By default, EAGO uses Clp for solving linear subproblems introduced. Using a
commercial linear solver is highly recommended such as Gurobi, CPLEX, or XPRESS
is highly recommended. Both Gurobi and CPLEX are free for academics and
installation information can be found through http://www.gurobi.com/academia/academia-center and
https://www.ibm.com/developerworks/community/blogs/jfp/entry/CPLEX_Is_Free_For_Students?lang=en, respectively.  

### Ipopt Build

Ipopt is the recommended solver for upper bounding problems and is supported for
algorithmically optimization problems in the MathProgBase interface as well as
problems defined through the JuMP interface. Ipopt's performance is highly
dependent on the linear algebra package used (up to 30x). By default MUMPS is used.
It's recommended that you either compile Ipopt with HSL MA57 or the Pardiso linear
algebra packages with a machine specific Blas library (for Intel users the JuliaPro
MKL version is recommended). For information on this, see the below links:

- Compiling Ipopt: https://www.coin-or.org/Ipopt/documentation/node13.html
- Julia Specifics:
   - Pointing Ipopt to a compiled version:
      - Ipopt Package Info: https://github.com/JuliaOpt/Ipopt.jl
      - Discourse discussion: https://discourse.julialang.org/t/use-ipopt-with-custom-version/9176
   - Issues using Pardiso:
      - Ubuntu: https://github.com/JuliaOpt/Ipopt.jl/issues/106
      - Windows: https://github.com/JuliaOpt/Ipopt.jl/issues/83
- HSL Website: http://www.hsl.rl.ac.uk/ipopt/
- Pardiso Website: https://pardiso-project.org/
