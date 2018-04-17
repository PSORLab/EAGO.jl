## **McCormick Operator Capabilities**

EAGO provides a library of McCormick relaxations in native Julia code. It supports
relaxing functions using both **nonsmooth McCormick relaxations** (Mitsos2009), **smooth McCormick relaxations** (Khan2017), **multi-variant McCormick relaxations** (Tsoukalas2014), as well
as **subgradient-based interval refinement** (Najman2017). For functions with
arbitrarily differentiable relaxations, the differentiable can be modified by adjusting a constant value. Additionally, validated interval bounds are supported via **ValidatedNumerics.jl** and nonvalidated interval operators are available through use of the **MCInterval{T}** object included in EAGO (which is essentially just a copy of the **IntervalArithmetic.jl** and **IntervalContractor.jl** library with corrected rounding features removed).

For details on constructing relaxations of functions, please see the

## **Currently supported operators**

The operators currently supported are listed below. The operators with a check box
have been subject to a large degree of scrutiny and are near optimal implementations.

### **Univariate McCormick Operators**

Arbitrarily differentiable relaxations can be constructed for the following operators:

- [x] **Inverse** (inv)
- [ ] **Logarithms** (log, log2, log10)
- [ ] **Exponential Functions** (exp, exp2, exp10)
- [ ] **Square Root** (sqrt)
- [ ] **Absolute Value** (abs)

Both nonsmooth and Whitney-1 (once differentiable) relaxations are supported:

- [ ] **Step Functions** (step, sign)
- [ ] **Trignometric Functions** (sin, cos, tan)
- [ ] **Inverse Trignometric Functions** (asin, acos, atan)
- [ ] **Hyperbolic Functions** (sinh, cosh, tanh)
- [ ] **Inverse Hyperbolic Functions** (asinh, acosh, atanh)

### **Bivariate McCormick/McCormick Operators**

The following bivariant operators are supported for two **SMCg** objects. Both nonsmooth and Whitney-1 (once differentiable) relaxations are supported.

- [x] **multiplication** (\*)
- [x] **division** (/)

Arbitrarily differentiable relaxations can be constructed for the following operators:

- [x] **addition** (+)
- [x] **subtraction** (-)
- [ ] **minimization** (min)
- [ ] **maximization** (max)

### **Bivariate McCormick/(Integer/Float) Operators**

Arbitrarily differentiable relaxations can be constructed for the following operators:

- [x] **addition** (+)
- [x] **subtraction** (-)
- [x] **multiplication** (\*)
- [x] **division** (/)
- [x] **minimization** (min)
- [x] **maximization** (max)
- [x] **Exponentiation** (pow, ^)
