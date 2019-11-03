# **Currently supported operators**

The operators currently supported are listed below. The operators with a check box
have been subject to a large degree of scrutiny and have been implemented for
both forward and reverse McCormick relaxations.

### **Univariate McCormick Operators**

Arbitrarily differentiable relaxations can be constructed for the following operators:

- [x] **Inverse** (inv)
- [x] **Logarithms** (log, log2, log10)
- [x] **Exponential Functions** (exp, exp2, exp10)
- [x] **Square Root** (sqrt)
- [x] **Absolute Value** (abs)

Both nonsmooth and Whitney-1 (once differentiable) relaxations are supported:

- [x] **Step Functions** (step, sign)
- [x] **Trignometric Functions** (sin, cos, tan)
- [x] **Inverse Trignometric Functions** (asin, acos, atan)
- [x] **Hyperbolic Functions** (sinh, cosh, tanh)
- [x] **Inverse Hyperbolic Functions** (asinh, acosh, atanh)

### **Bivariate Operators: McCormick & McCormick**

The following bivariant operators are supported for two **MC** objects. Both nonsmooth and Whitney-1 (once differentiable) relaxations are supported.

- [x] **multiplication** (\*)
- [x] **division** (/)

Arbitrarily differentiable relaxations can be constructed for the following operators:

- [x] **addition** (+)
- [x] **subtraction** (-)
- [x] **minimization** (min)
- [x] **maximization** (max)

### **Bivariate Operators: McCormick & (Integer or Float)**

Arbitrarily differentiable relaxations can be constructed for the following operators:

- [x] **addition** (+)
- [x] **subtraction** (-)
- [x] **multiplication** (\*)
- [x] **division** (/)
- [x] **minimization** (min)
- [x] **maximization** (max)
- [x] **Exponentiation** (pow, ^)
