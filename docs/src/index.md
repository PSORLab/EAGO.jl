
![full_Logo](full_Logo1.png)

# **EAGO - Easy Advanced Global Optimization in Julia**

A flexible framework for global and robust optimization in Julia.

## Authors
- [Matthew Wilhelm](https://psor.uconn.edu/person/matthew-wilhelm/), Department of Chemical and Biomolecular Engineering, University of Connecticut (UConn)
- [Robert Gottlieb](https://psor.uconn.edu/person/robert-gottlieb/), Department of Chemical and Biomolecular Engineering, University of Connecticut (UConn)
- [Dimitri Alston](https://psor.uconn.edu/person/dimitri-alston/), Department of Chemical and Biomelecular Engineering, University of Connecticut (UConn)

## Overview
EAGO is a global and robust optimization platform based on McCormick relaxations.
It contains the first widely accessible global optimization routine based on
generalized McCormick relaxations. With the exception of calls to local solvers
and linear algebra routines, EAGO is written entirely in native Julia.
The solver is flexibly arranged so the end user can easily customize low-level routines.

## Installing EAGO
EAGO is a registered Julia package and it can be installed using the Julia package manager.
From the Julia REPL, type `]` to enter the Package manager (Pkg) mode and run the following command:

```julia
pkg> add EAGO
```

Currently, EAGO is compatible with version 1.11+ of JuMP. This allows a replication of some of the internal features shared by EAGO and JuMP's automatic differentiation scheme, e.g., generation of Wergert Tapes, passing evaluators between JuMP and EAGO, etc.

```julia
pkg> add JuMP
```

EAGO v0.8.1 is the current tagged version and requires Julia 1.6+ for full functionality (however Julia 1.0+ versions support partial functionality). Use with version 1.8 is recommended as the majority of in-house testing has occurred using this version of Julia. The user is directed to the [**High-Performance Configuration**](https://psorlab.github.io/EAGO.jl/Optimizer/high_performance/) for instructions on how to install a high performance version of EAGO (rather than the basic entirely open-source version).
If any issues are encountered when loading EAGO (or when using it), please submit an issue using the Github [**issue tracker**](https://github.com/PSORLab/EAGO.jl/issues).

## Examples
Several examples are provided within this documentation, but additional examples
are provided in the form of Jupyter notebooks at [**EAGO-notebooks**](https://github.com/PSORLab/EAGO-notebooks),
which can be run using IJulia. To add IJulia, run the command:

```julia
pkg> add IJulia
```
Then launch the Jupyter notebook using the following command from the Julia terminal:

```julia
julia> using IJulia; notebook()
```

And then simply navigate to the example directory and run the example of most interest.
