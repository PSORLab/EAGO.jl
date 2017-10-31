# EAGO
A development environment for robust and global optimization.

This project is intended to eventually take the form a gui driven analytic tool for advanced optimization problems. The preliminary work on this package consists of a series of global optimization routines based on McCormick relaxations and semi-infinite program solution routines.

The subpackages available are currently:
- EAGOBranchBound: A fully customizable branch and bound library.
- EAGOSmoothMcCormick: A library of smooth McCormick relaxations.
- [EAGODAGContractor](https://github.com/MatthewStuber/EAGODAGContractor): A forward-reverse interval contractor routine for constraint propagation on directed graphs with visualization tools.
- EAGOGlobalSolver: A series of global solvers based on McCormick relaxations.
- EAGOSemiInfinite: Implements a general nonconvex solver for semi-infinite programs via a restriction of the righthandside approach.
