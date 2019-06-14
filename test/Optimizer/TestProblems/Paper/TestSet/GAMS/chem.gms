$title Chemical Equilibrium Problem (CHEM,SEQ=21)

$onText
The problem is to find the equilibrium composition of a
mixture of different atoms.


Bracken, J, and McCormick, G P, Chapter 5. In Selected Applications
of Nonlinear Programming. John Wiley and Sons, New York, 1968, pp. 48-49.

Keywords: nonlinear programming, chemical equilibrium problem, chemical engineering
$offText

Set
   c 'compounds' / H, H2, H2O, N, N2, NH, NO, O, O2, OH   /
   i 'atoms'     / H 'hydrogen', N 'nitrogen', O 'oxygen' /;

Table a(i,c) 'atoms per compound'
       H  H2  H2O  N  N2  NH  NO  O  O2  OH
   H   1   2    2          1              1
   N               1   2   1   1
   O            1              1  1   2   1;

Parameter
   mix(i)   'number of elements in mixture' / h 2, n 1, o 1 /
   gibbs(c) 'gibbs free energy at 3500 k and 750 psi'
            / H  -10.021, H2  -21.096, H2O -37.986, N   -9.846, N2 -28.653
              NH -18.918, NO -28.032 , O   -14.640, o2 -30.594, OH -26.11  /
   gplus(c) 'gibbs energy plus pressure';

gplus(c) = gibbs(c) + log(750*.07031);

display gplus;

Variable
   x(c)    'number of mols in mixture'
   xb      'total number of mols in mixture'
   energy  'total free energy in mixture';

Positive Variable x, xb;

Equation
   cdef(i) 'component definition'
   edef    'energy definition'
   xdef    'total mol definition';

cdef(i)..  sum(c, a(i,c)*x(c)) =e= mix(i);

xdef..     xb     =e= sum(c, x(c));

edef..     energy =e= sum(c, x(c)*(gplus(c) + log(x(c)/xb)));

x.lo(c) = .001;
xb.lo   = .01;

Model mixer 'chemical mix for N2H4+O2' / all /;

solve mixer minimizing energy using nlp;
