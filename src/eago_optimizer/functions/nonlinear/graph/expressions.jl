const UNIVARIATE_EVAL = Dict{AtomType, Symbol}()

UNIVARIATE_EVAL[INV]       = :inv
UNIVARIATE_EVAL[ONE]       = :one
UNIVARIATE_EVAL[ZERO]      = :zero
UNIVARIATE_EVAL[REAL]      = :real
UNIVARIATE_EVAL[EPS]       = :eps
UNIVARIATE_EVAL[ABS]       = :abs
UNIVARIATE_EVAL[ABS2]      = :abs2
UNIVARIATE_EVAL[RAD2DEG]   = :rad2deg
UNIVARIATE_EVAL[DEG2RAD]   = :deg2rad
UNIVARIATE_EVAL[STEP]      = :step
UNIVARIATE_EVAL[SIGN]      = :sign

UNIVARIATE_EVAL[LOG]       = :log
UNIVARIATE_EVAL[LOG2]      = :log2
UNIVARIATE_EVAL[LOG10]     = :log10
UNIVARIATE_EVAL[LOG1P]     = :log1p
UNIVARIATE_EVAL[EXP]       = :exp
UNIVARIATE_EVAL[EXP2]      = :exp2
UNIVARIATE_EVAL[EXP10]     = :exp10
UNIVARIATE_EVAL[EXP10]     = :expm1

UNIVARIATE_EVAL[SIN]       = :sin
UNIVARIATE_EVAL[COS]       = :cos
UNIVARIATE_EVAL[TAN]       = :tan
UNIVARIATE_EVAL[CSC]       = :csc
UNIVARIATE_EVAL[SEC]       = :sec
UNIVARIATE_EVAL[COT]       = :cot
UNIVARIATE_EVAL[ASIN]      = :asin
UNIVARIATE_EVAL[ACOS]      = :acos
UNIVARIATE_EVAL[ATAN]      = :atan
UNIVARIATE_EVAL[ACSC]      = :acsc
UNIVARIATE_EVAL[ASEC]      = :asec
UNIVARIATE_EVAL[ACOT]      = :acot
UNIVARIATE_EVAL[SINH]      = :sinh
UNIVARIATE_EVAL[COSH]      = :cosh
UNIVARIATE_EVAL[TANH]      = :tanh
UNIVARIATE_EVAL[CSCH]      = :csch
UNIVARIATE_EVAL[SECH]      = :sech
UNIVARIATE_EVAL[COTH]      = :coth
UNIVARIATE_EVAL[ASINH]     = :asinh
UNIVARIATE_EVAL[ACOSH]     = :acosh
UNIVARIATE_EVAL[ATANH]     = :atanh
UNIVARIATE_EVAL[ACSCH]     = :acsch
UNIVARIATE_EVAL[ASECH]     = :asech
UNIVARIATE_EVAL[ACOTH]     = :acoth

UNIVARIATE_EVAL[ERF]       = :erf
UNIVARIATE_EVAL[ERFC]      = :erfc
UNIVARIATE_EVAL[ERFINV]    = :erfinv
UNIVARIATE_EVAL[ERFCINV]   = :erfcinv

UNIVARIATE_EVAL[SQRT]      = :sqrt
UNIVARIATE_EVAL[CBRT]      = :cbrt

UNIVARIATE_EVAL[POS]       = :positive
UNIVARIATE_EVAL[NEG]       = :negative
UNIVARIATE_EVAL[USER]      = :user

UNIVARIATE_EVAL[RELU]       = :relu
UNIVARIATE_EVAL[LEAKY_RELU] = :leaky_relu
UNIVARIATE_EVAL[MAXSIG]     = :maxsig
UNIVARIATE_EVAL[MAXTANH]    = :maxtanh
UNIVARIATE_EVAL[SOFTPLUS]   = :softplus
UNIVARIATE_EVAL[PENTANH]    = :pentanh
UNIVARIATE_EVAL[BISIGMOID]  = :bisigmoid
UNIVARIATE_EVAL[SOFTSIGN]   = :softsign
UNIVARIATE_EVAL[GELU]       = :gelu
UNIVARIATE_EVAL[SILU]       = :silu
UNIVARIATE_EVAL[XLOGX]      = :xlogx
UNIVARIATE_EVAL[XABSX]      = :xabsx

const UNIVARIATE_ATOM_INSTANCES = AtomType[k for k in keys(UNIVARIATE_EVAL)]

const BIVARIATE_EVAL = Dict{AtomType, Symbol}()
BIVARIATE_EVAL[DIV]       = :/
BIVARIATE_EVAL[POW]       = :^
BIVARIATE_EVAL[LOWER_BND] = :lower_bnd
BIVARIATE_EVAL[UPPER_BND] = :upper_bnd
BIVARIATE_EVAL[ARH]       = :arh
#BIVARIATE_EVAL[EXPAX]     = :lower_bnd
#BIVARIATE_EVAL[EXPXY]     = :lower_bnd

const NARITY_EVAL = Dict{AtomType, Symbol}()
NARITY_EVAL[MIN]  = :min
NARITY_EVAL[MAX]  = :max
NARITY_EVAL[PLUS] = :+
NARITY_EVAL[MULT] = :*
NARITY_EVAL[USERN]   = :usern

const ATM_EVAL = union(UNIVARIATE_EVAL, BIVARIATE_EVAL, NARITY_EVAL)
ATM_EVAL[MINUS]     = :-
ATM_EVAL[BND]       = :bnd
#ATM_EVAL[QUAD1]     = :quad1
#ATM_EVAL[QUAD2]     = :quad2
