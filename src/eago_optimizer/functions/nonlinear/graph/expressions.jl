# Definitions borrow from https://github.com/FluxML/NNlib.jl (names used to
# standardize package). TODO: Decide how/if to incorporate NNlib depedency into
# McCormick.jl/EAGO.jl.
oftf(x, y) = oftype(float(x), y)
leakyrelu(x, a=oftf(x, 0.01)) = max(a * x, x)
swish(x) = x * sigmoid(x)

"""
    AtomType
"""
@enum(AtomType, VAR_ATOM, PARAM_ATOM, CONST_ATOM, SELECT_ATOM,
                ABS, ABS2, INV, RAD2DEG, DEG2RAD, ONE, ZERO,
                MULT, DIV, PLUS, MINUS, POW, MIN, MAX, STEP, SIGN,
                LOG, LOG2, LOG10, LOG1P, EXP, EXP2, EXP10, EXPM1,
                SIN, COS, TAN, CSC, SEC, COT,
                ASIN, ACOS, ATAN, ACSC, ASEC, ACOT,
                SINH, COSH, TANH, CSCH, SECH, COTH,
                ASINH, ACOSH, ATANH, ACSCH, ASECH, ACOTH,
                ERF, ERFC, ERFINV, ERFCINV, SQRT, CBRT,
                # ML specific functions
                RELU, LEAKYRELU, SOFTPLUS, SOFTSIGN, GELU, SILU, SIGMOID, SWISH,
                # EAGO modeling functions
                XLOGX, ARH,
                # EAGO bound enforcing functions
                POS, NEG, LOWER_BND, UPPER_BND, BND,
                # Atoms for nonexpression types
                USER, USERN, SUBEXPR)


# Define univariate functions that correspond to the AtomType
const UNIVARIATE_ATOM_DICT = Dict{AtomType, Symbol}()

UNIVARIATE_ATOM_DICT[INV]       = :inv
UNIVARIATE_ATOM_DICT[ONE]       = :one
UNIVARIATE_ATOM_DICT[ZERO]      = :zero
UNIVARIATE_ATOM_DICT[ABS]       = :abs
UNIVARIATE_ATOM_DICT[ABS2]      = :abs2
UNIVARIATE_ATOM_DICT[RAD2DEG]   = :rad2deg
UNIVARIATE_ATOM_DICT[DEG2RAD]   = :deg2rad
UNIVARIATE_ATOM_DICT[STEP]      = :step
UNIVARIATE_ATOM_DICT[SIGN]      = :sign

UNIVARIATE_ATOM_DICT[LOG]       = :log
UNIVARIATE_ATOM_DICT[LOG2]      = :log2
UNIVARIATE_ATOM_DICT[LOG10]     = :log10
UNIVARIATE_ATOM_DICT[LOG1P]     = :log1p
UNIVARIATE_ATOM_DICT[EXP]       = :exp
UNIVARIATE_ATOM_DICT[EXP2]      = :exp2
UNIVARIATE_ATOM_DICT[EXP10]     = :exp10
UNIVARIATE_ATOM_DICT[EXP10]     = :expm1

UNIVARIATE_ATOM_DICT[SIN]       = :sin
UNIVARIATE_ATOM_DICT[COS]       = :cos
UNIVARIATE_ATOM_DICT[TAN]       = :tan
UNIVARIATE_ATOM_DICT[CSC]       = :csc
UNIVARIATE_ATOM_DICT[SEC]       = :sec
UNIVARIATE_ATOM_DICT[COT]       = :cot
UNIVARIATE_ATOM_DICT[ASIN]      = :asin
UNIVARIATE_ATOM_DICT[ACOS]      = :acos
UNIVARIATE_ATOM_DICT[ATAN]      = :atan
UNIVARIATE_ATOM_DICT[ACSC]      = :acsc
UNIVARIATE_ATOM_DICT[ASEC]      = :asec
UNIVARIATE_ATOM_DICT[ACOT]      = :acot
UNIVARIATE_ATOM_DICT[SINH]      = :sinh
UNIVARIATE_ATOM_DICT[COSH]      = :cosh
UNIVARIATE_ATOM_DICT[TANH]      = :tanh
UNIVARIATE_ATOM_DICT[CSCH]      = :csch
UNIVARIATE_ATOM_DICT[SECH]      = :sech
UNIVARIATE_ATOM_DICT[COTH]      = :coth
UNIVARIATE_ATOM_DICT[ASINH]     = :asinh
UNIVARIATE_ATOM_DICT[ACOSH]     = :acosh
UNIVARIATE_ATOM_DICT[ATANH]     = :atanh
UNIVARIATE_ATOM_DICT[ACSCH]     = :acsch
UNIVARIATE_ATOM_DICT[ASECH]     = :asech
UNIVARIATE_ATOM_DICT[ACOTH]     = :acoth

UNIVARIATE_ATOM_DICT[ERF]       = :erf
UNIVARIATE_ATOM_DICT[ERFC]      = :erfc
UNIVARIATE_ATOM_DICT[ERFINV]    = :erfinv
UNIVARIATE_ATOM_DICT[ERFCINV]   = :erfcinv

UNIVARIATE_ATOM_DICT[SQRT]      = :sqrt
UNIVARIATE_ATOM_DICT[CBRT]      = :cbrt

UNIVARIATE_ATOM_DICT[POS]       = :positive
UNIVARIATE_ATOM_DICT[NEG]       = :negative
UNIVARIATE_ATOM_DICT[USER]      = :user

UNIVARIATE_ATOM_DICT[RELU]       = :relu
UNIVARIATE_ATOM_DICT[LEAKYRELU]  = :leakyrelu
UNIVARIATE_ATOM_DICT[SIGMOID]    = :sigmoid
UNIVARIATE_ATOM_DICT[GELU]       = :gelu
UNIVARIATE_ATOM_DICT[SOFTPLUS]   = :softplus
UNIVARIATE_ATOM_DICT[SOFTSIGN]   = :softsign
UNIVARIATE_ATOM_DICT[SWISH]      = :swish

UNIVARIATE_ATOM_DICT[XLOGX]      = :xlogx

# Define bivariate functions that correspond to the AtomType
const BIVARIATE_ATOM_DICT = Dict{AtomType, Symbol}()
BIVARIATE_ATOM_DICT[DIV]       = :/
BIVARIATE_ATOM_DICT[POW]       = :^
BIVARIATE_ATOM_DICT[LOWER_BND] = :lower_bnd
BIVARIATE_ATOM_DICT[UPPER_BND] = :upper_bnd
BIVARIATE_ATOM_DICT[ARH]       = :arh
#BIVARIATE_ATOM_DICT[EXPAX]     = :lower_bnd
#BIVARIATE_ATOM_DICT[EXPXY]     = :lower_bnd

# Define n-arity functions that correspond to the AtomType
const NARITY_ATOM_DICT = Dict{AtomType, Symbol}()
NARITY_ATOM_DICT[MIN]  = :min
NARITY_ATOM_DICT[MAX]  = :max
NARITY_ATOM_DICT[PLUS] = :+
NARITY_ATOM_DICT[MULT] = :*
NARITY_ATOM_DICT[USERN]   = :usern

const ALL_ATOM_DICT = Dict{AtomType, Symbol}()
foreach(x -> setindex!(ALL_ATOM_DICT, x[2], x[1]), UNIVARIATE_ATOM_DICT)
foreach(x -> setindex!(ALL_ATOM_DICT, x[2], x[1]), BIVARIATE_ATOM_DICT)
foreach(x -> setindex!(ALL_ATOM_DICT, x[2], x[1]), NARITY_ATOM_DICT)

# A functions that may be 1 to n-arity functions that correspond to the AtomType
ALL_ATOM_DICT[MINUS]     = :-
ALL_ATOM_DICT[BND]       = :bnd
#ATM_EVAL[QUAD1]     = :quad1
#ATM_EVAL[QUAD2]     = :quad2

# List of keys only
const UNIVARIATE_ATOM_TYPES = AtomType[k for k in keys(UNIVARIATE_ATOM_DICT)]
const BIVARIATE_ATOM_TYPES = AtomType[k for k in keys(BIVARIATE_ATOM_DICT)]
const NARITY_ATOM_TYPES = AtomType[k for k in keys(NARITY_ATOM_DICT)]
const ALL_ATOM_TYPES = AtomType[k for k in keys(ALL_ATOM_DICT)]

# Reverse lookup dicts
const REV_UNIVARIATE_ATOM_DICT = Dict(UNIVARIATE_ATOM_DICT[k] => k for k in keys(UNIVARIATE_ATOM_DICT))
const REV_BIVARIATE_ATOM_DICT = Dict(BIVARIATE_ATOM_DICT[k] => k for k in keys(BIVARIATE_ATOM_DICT))
const REV_NARITY_ATOM_DICT = Dict(NARITY_ATOM_DICT[k] => k for k in keys(NARITY_ATOM_DICT))
