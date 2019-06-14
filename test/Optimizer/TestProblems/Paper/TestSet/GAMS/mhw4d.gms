$title Nonlinear Test Problem (MHW4D,SEQ=84)

$onText
Another popular test problem for NLP codes.


Wright, M H, Numerical Methods for Nonlinearly Constraint Optimization.
PhD thesis, Stanford University, 1976.

Keywords: nonlinear programming, mathematics
$offText

Variable m, x1, x2, x3, x4, x5;

Equation funct, eq1, eq2, eq3;

funct.. m =e= sqr(x1-1) + sqr(x1-x2) + power(x2-x3,3) + power(x3-x4,4) + power(x4-x5,4);

eq1.. x1 + sqr(x2) + power(x3,3) =e= 3*sqrt(2) + 2;

eq2.. x2 - sqr(x3) + x4          =e= 2*sqrt(2) - 2;

eq3.. x1*x5 =e= 2;

Model wright / all /;
wright.optca = 1E-3;
wright.optcr = 1E-3;

x1.l = -1;
x2.l =  2;
x3.l =  1;
x4.l = -2;
x5.l = -2;

solve wright using nlp minimizing m;
