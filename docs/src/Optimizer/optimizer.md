# EAGO Optimizer

The `EAGO.Optimizer` object holds all algorithm solution information. A description
of all user-facing options has been provided in the docstring.

## EAGO.Optimizer
```@docs
Optimizer
```

## Storage for Input Parameters
```@docs
EAGO.EAGOParameters
```

## Internal Storage Structures
```@docs
VariableInfo
ExtensionType
Log
```

## Internal Problem Representations
```@docs
EAGO.InputProblem
EAGO.ParsedProblem
```

## Extending EAGO

Functionality has been included that allows for extension's to EAGO's optimizer
to readily be defined. This can be done in two ways first defining a new structure
which is a subtype of `EAGO.ExtensionType` and overloading methods associated with
this new structure. An instance of this new structure is provided to the `EAGO.Optimizer`
using the `ext_type` keyword. This results in EAGO now dispatch to the new
methods rather than the generally defined methods for the parent type. For a complete
example, the reader is directed to the [**interval bounding example**](https://github.com/PSORLab/EAGO-notebooks/blob/master/notebooks/nlpopt_interval_bnb.ipynb) and [**quasiconvex example**](https://github.com/PSORLab/EAGO-notebooks/blob/master/notebooks/custom_quasiconvex.ipynb). Alternatively, the user can overload the `optimize_hook!` for
this subtype which will entirely circumvent the default global solution routine. Additional
information can be stored in the `ext` field of EAGO. In order to allow for compatibility
between packages the user is encouraged to append their extension name to the start of each
variable name (e.g. `newext_newdata`).
