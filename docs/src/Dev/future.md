# Future Work

## Current Activity:
* Update CI testing.
* Specialized algorithms for relaxing ODE constrained problems and solving
global and robust optimization problems.
* Development of dedicated linear algebra libraries for McCormick operators.

## Other things on the wishlist (Not actively being worked on but certainly desirable):
* Implement the interval constraint propagation scheme presented in Vu 2008. For
improved convergences.
* A parametric bisection routine will be updated that can divide the `(X,P)`
space into a series of boxes that all contain unique branches of the implicit
function `p->y(p)`.
* Provide a better interface the nonconvex semi-infinite programs solvers
(JuMPeR extension?).
* Add additional McCormick relaxations.
* Add handling for domain reduction of special expression forms.
