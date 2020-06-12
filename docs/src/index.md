
![full_Logo](full_Logo1.png)

# **EAGO - Easy Advanced Global Optimization in Julia**

A flexible-framework for global and robust optimization in Julia

## Authors
- [Matthew Wilhelm](https://psor.uconn.edu/person/matthew-wilhelm/), Department of Chemical and Biomolecular Engineering,  University of Connecticut (UCONN)

## Overview
**EAGO** is a global and robust optimization platform based on McCormick relaxations.
It contains the first widely accessible global optimization routine based on
generalized McCormick relaxations. With the exception of calls to local solvers
and linear algebra routines, EAGO is written entirely in native Julia.
The solver is quite flexibly arranged so the end user can easily customize low-level routines.

## Installing EAGO
EAGO is registered Julia package. It can be installed using the Julia package manager.
From the Julia REPL, type ] to enter the Pkg REPL mode and run the following command

```julia
pkg> add EAGO
```

Currently, EAGO is tied to a 0.21.2+ or greater version of JuMP. This allows a replication
of some of the internal features shared by EAGO and JuMP's AD scheme aka
generation of Wergert Tapes pass evaluators between JuMP and EAGO etc.

```julia
pkg> add JuMP
```

EAGO v0.4 is the current version requires Julia 1.2+. Use with version 1.4 is recommended as the majority of in-house testing has occured using this version of Julia. The user is directed to the [**High-Performance Configuration**](https://psorlab.github.io/EAGO.jl/Optimizer/high_performance/) for instructions on how to install a high performance version of EAGO (rather than the basic entirely open-source version). If any issues are encountered when loading EAGO (or when using it), please submit an issue using the Github [**issue tracker**](https://github.com/PSORLab/EAGO.jl/issues).

## Examples
A few examples are provided in the documentation website. More involved
examples are provided at in the form of Jupyter Notebooks at [**EAGO-notebooks**](https://github.com/PSORLab/EAGO-notebooks) and can be run using
IJulia. To add IJulia

```julia
pkg> add IJulia
```
Then launch the Jupyter notebook using the following command from the Julia terminal,

```julia
julia> using IJulia; notebook()
```

Then simply navigate to the example directory and run the example of most interest.
