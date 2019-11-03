
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

Currently, EAGO is tied to a 0.19 or greater version of JuMP. This allows a replication
of some of the internal features shared by EAGO and JuMP's AD scheme aka
generation of Wergert Tapes pass evaluators between JuMP and EAGO etc.

```julia
pkg> add JuMP
```

Once JuMP 0.19+ is registered EAGO will be updated to eliminate the need for this.

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
