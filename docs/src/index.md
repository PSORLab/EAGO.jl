# EAGO - Easy Advanced Global Optimization in Julia

## Authors
- [Matthew Wilhelm](https://psor.uconn.edu/person/matthew-wilhelm/), Department of Chemical and Biomolecular Engineering,  University of Connecticut (UCONN)

## Overview
**EAGO** is a global and robust optimization platform based on McCormick relaxations.
It contains the first widely accessible global optimization routine based on
generalized McCormick relaxations. With the exception of calls to local solvers
and linear algebra routines, EAGO is written entirely in native Julia.
The solver is quite flexibly arranged so the end user can easily customize low-level routines.

## Installing EAGO
EAGO is registered Julia package and can be installed by running:

```julia
julia> Pkg.add("EAGO")
```
