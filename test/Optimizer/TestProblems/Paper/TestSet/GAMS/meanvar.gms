$title Financial Optimization: Risk Management (MEANVAR,SEQ=112)

$onText
The mean-variance model is often used to manage residual risk in
an elegant way. This class of models is due to H Markowitz who
introduced these concepts in 1952. Other models like [IMMUN] hedge
against different kinds of systematic risk.

The family of models presented in this file are:

VAR1   Minimum Variance Portfolio Model
       Global Asset Allocation - Efficient Frontier Model
VAR2   Riskless Asset Minimum Variance Model
SHARPE Minimum Variance Model using Sharpe Ratio
MEAN   Mean-Variance Model
MARG   Mean-Variance Model with Marginal Variance Definition


Dahl, H, Meeraus, A, and Zenios, S A, Some Financial Optimization
Models: Risk Management. In Zenios, S A, Ed, Financial Optimization.
Cambridge University Press, New York, NY, 1993.

keywords: nonlinear programming, risk management, portfolio optimization, finance
          financial optimization
$offText

Set
   s    'investment set'       / cn, fr, gr, jp, sw, uk, us, wr /
   i(s) 'analyzed investments' / cn, fr, gr, jp, sw, uk, us     /;

Alias (i,j);

* note: this instance represents investments of hedged foreign securities.
Parameter mu(s) 'expected return of security' / cn  0.1287
                                                fr  0.1096
                                                gr  0.0501
                                                jp  0.1524
                                                sw  0.0763
                                                uk  0.1854
                                                us  0.0620
                                                wr  0.0916 /;

Table q(i,j) 'covariance matrix'
           cn     fr     gr     jp     sw     uk     us
   cn   42.18
   fr   20.18  70.89
   gr   10.88  21.58  25.51
   jp    5.30  15.41   9.60  22.33
   sw   12.32  23.24  22.63  10.32  30.01
   uk   23.84  23.80  13.22  10.46  16.36  42.23
   us   17.41  12.62   4.70   1.00   7.20   9.90  16.42;

q(i,j)$(ord(j) > ord(i)) = q(j,i);

Scalar mup 'target expected return for the portfolio' / 0.115 /;

Variable
   v    'variance'
   m    'mean return'
   x(i) 'fraction of the portfolio that consists of security i';

Positive Variable x;
x.up(i) = 1;

Equation
   vbal   'variance definition'
   mbal   'mean balancing constraint'
   budget 'budget constraint';

vbal..   v =e= sum((i,j), x(i)*q(i,j)*x(j));

mbal..   m =e= sum(i, mu(i)*x(i));

budget.. sum(i, x(i)) =e= 1;

m.fx = mup;

Model var1 / vbal, mbal, budget /;

solve var1 minimizing v using nlp;

$sTitle Riskless Asset Minimum Variance Model
Scalar r 'return of riskless asset' / 0.0979 /;

Equation riskless 'mean balancing constraint using riskless rate';

riskless..  sum(i, (mu(i)-r)*x(i)) =e= mup - r;

Model var2 / vbal, riskless /;

solve var2 minimizing v using nlp;

display mup, r;

$sTitle Minimum Variance Model using Sharpe Ratio
Set      h(s)  'historical investments' / wr /;

Scalar   sf    'shrinkage factor';

Variable omega 'objective';

Equation
  obj2         'objective function'
  meanbal      'mean balancing constraint';

obj2..     omega  =e=  m/sqrt(v);

meanbal..  sum(i, x(i)*(sf*sum(h, mu(h)) + (1 - sf)*mu(i))) =e= m;

m.l  = 1;
v.lo = 0.01;
v.l  = 1;
sf   = .65;
m.lo = 0;
m.up = inf;

Model sharpe / obj2, meanbal, vbal, budget /;

solve sharpe maximizing omega using nlp;

$sTitle Global Asset Allocation - efficient frontier Model
Set
   p    'efficient frontier points' / minvar, p1*p4, maxvar /
   pp   'efficient frontier points' /         p1*p4         /;

Scalar
   vmin 'minimum variance'
   vmax 'maximum variance';

Parameter xres(*,p) 'portfolios at different points';

option limCol = 0, limRow = 0;
var1.solPrint=%solPrint.quiet%;

loop(p('maxvar'),
   solve var1 maximize m using nlp;
   xres(i,p)        = x.l(i);
   xres('mean',p)   = m.l;
   xres('var',p)    = v.l;
   xres('status',p) = var1.modelStat;
   vmax = v.l;
);

loop(p('minvar'),
   solve var1 minimize v using nlp;
   xres(i,p)        = x.l(i);
   xres('mean',p)   = m.l;
   xres('var',p)    = v.l;
   xres('status',p) = var1.modelStat;
   vmin = v.l;
);

loop(p(pp),
   v.fx = vmin + (vmax - vmin)/(card(pp) + 1)*ord(pp);
   solve var1 maximizing m using nlp;
   xres(i,p)        = x.l(i);
   xres('mean',p)   = m.l;
   xres('var',p)    = v.l;
   xres('status',p) = var1.modelStat;
);

display xres;
