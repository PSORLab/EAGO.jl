# Relaxation of Implicit Functions

## High-Level Functions

```@docs
implicit_relax_h!
```

## Data Structures

```@docs
McCormick.AbstractContractorMC
NewtonGS
KrawczykCW
McCormick.AbstractPreconditionerMC
DenseMidInv
McCormick.AbstractMCCallback
MCCallback
```

## Subroutines

```@docs
preconditioner_storage
affine_exp!
correct_exp!
contract!
precond_and_contract!
final_cut
gen_expansion_params!
populate_affine!
```
