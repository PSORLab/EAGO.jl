# Future Work

## Current Activity

* Update CI testing.
* Specialized algorithms for relaxing ODE constrained problems and solving global and robust optimization problems.
* Extensions for nonconvex dynamic global & robust optimization.
* Provide support for mixed-integer problems.
* Update EAGO to support nonsmooth problems (requires: a nonsmooth local nlp optimizer or lexicographic AD, support for relaxations is already included).
* Evaluation and incorporation of implicit relaxation routines in basic solver.

## Other Things on the Wishlist (But Not Actively Being Worked On)

* Implement the interval constraint propagation scheme presented in Vu 2008. For improved convergences.
* A parametric bisection routine will be updated that can divide the `(X,P)` space into a series of boxes that all contain unique branches of the implicit function `p->y(p)`.
* Provide a better interface the nonconvex semi-infinite programs solvers (JuMPeR extension?).
* Add additional McCormick relaxations.
* Add handling for domain reduction of special expression forms.
