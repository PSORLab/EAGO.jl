## Current Activity:
* Upgrade the solver to use MathOptInteface.
* Add the ability to solve mixed-integer nonlinear programs via a
  Branch-and-Cut methodology.
* Adding support to support mixed relaxation types based on expressions
  and subexpressions.
* Add support for relaxing logical-statements.
* A GUI interface for flowsheeting problems.
* Support for propagating McCormick relaxations as expressions and taking
  higher order derivatives of the relaxations (port the previous
  EAGOSmoothMcCormick.jl package over the EAGO.jl package)

## Other things on the wishlist (Not actively being worked on but certainly desirable):
* Implement the interval constraint propagation scheme presented in Vu 2008. For
improved convergences.
* Implement dedicated single-precision linear algebra routines for intervals and
McCormick operators.
* Dedicated handling for linear and quadratic constraints.
* A parametric bisection routine will be updated that can divide the `(X,P)`
space into a a series of boxes that all contain unique branches of the implicit 
function `p->y(p)`.
* Additional support for structured linear algebra solution routines
  with non-standard types.
* Provide a better interface the nonconvex semi-infinite programs solvers
(JuMPeR extension?).
* Add additional McCormick relaxations.
* Provide means of relaxing Julia code that isn't flexibly typed (Cassette.jl?)
