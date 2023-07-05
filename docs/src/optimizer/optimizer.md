# EAGO Optimizer

The `Optimizer` object holds all algorithm solution information. A description of all user-facing options has been provided in the docstring.

## EAGO.Optimizer

```@docs
Optimizer
```

## EAGO Specific Functions and Operators

EAGO supports a number of functions and operators that for which specialized relaxation routines are available. These can be registered and added to a JuMP model using the function:

```@docs
EAGO.register_eago_operators!(m::JuMP.Model)
```

## Storage for Input Parameters

```@docs
EAGO.EAGOParameters
```

## Internal Storage Structures

```@docs
VariableInfo
ExtensionType
```

## Internal Problem Representations

```@docs
EAGO.InputProblem
EAGO.ParsedProblem
```

## Interval Optimizer Subroutines

```@docs
EAGO.initial_parse!(m::Optimizer{R,S,T}) where {R,S,T}
```

## Extending EAGO

Functionality has been included that allows for extensions to EAGO's [`Optimizer`](@ref) to be readily defined. This can be done in two ways first defining a new structure which is a subtype of [`ExtensionType`](@ref) and overloading methods associated with this new structure. An instance of this new structure is provided to the [`Optimizer`](@ref) using the `ext_type` keyword. This results in EAGO now dispatch to the new methods rather than the generally defined methods for the parent type. For a complete example, the reader is directed to the [interval bounding example](@ref "Standard-Use Example 2") and the [quasiconvex example](@ref "Advanced-Use Example 1"). Alternatively, the user can overload the `optimize_hook!` for this subtype which will entirely circumvent the default global solution routine. Additional information can be stored in the `ext` field of the [`Optimizer`](@ref). In order to allow for compatibility between packages the user is encouraged to append their extension name to the start of each variable name (e.g. `newext_newdata`).
