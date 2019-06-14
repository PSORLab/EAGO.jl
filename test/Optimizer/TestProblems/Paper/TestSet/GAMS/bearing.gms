$title Hydrostatic Thrust Bearing Design for a Turbogenerator (BEARING,SEQ=202)

$onText
The objective of this problem is to minimize the power loss
during the operation of a hydrostatic thrust bearing subject
to a number of constraints.


Siddall, J N, Optimal Engineering Design. Marcel Dekker,
New York, 1982.

Deb, K, and Goyal, M, Optimizing Engineering Designs Using
a Combined Genetic Search. In Back, T, Ed, Proceedings of
the Seventh International Conference on Genetic Algorithms.
1997, pp. 521-528.

Coello Coello, C A, Treating Constraints as Objectives for
Single-Objective Evolutionary Optimization. Engineering
Optimization 32 (2000), 275-308.

Keywords: nonlinear programming, engineering design, hydrostatic thrust bearing
$offText

Variable
   R       'bearing step radius               (in.)'
   R0      'recess radius                     (in.)'
   mu      'oil viscosity in reynolds  (lb sec/in.)'
   Q       'flow rate of oil            (in.^3/sec)'
   PL      'power loss                 (in. lb/sec)'
   P0      'inlet pressure              (psi gauge)'
   Ef      'power loss due to friction (in. lb/sec)'
   Ep      'pumping energy             (in. lb/sec)'
   h       'film thickness                    (in.)'
   delta_t 'temperature rise               (degr F)'
   t       'temperature                    (degr F)'
   W       'weight                             (lb)'
   tmp1    'temp variable: log(R/R0)'
   tmp2    'temp variable: R^2-R0^2';

Equation
   power_loss      'objective function'
   pumping_energy  'included pump efficiency'
   friction        'friction loss'
   temp_rise       'temperature rise of the oil passing through the bearing'
   inlet_pressure  'relation between inlet pressure and flow rate'
   load_capacity   'load-carrying capacity'
   oil_viscosity   'oil viscosity is a function of temperature'
   temperature     'T is mean of ambient and oil outlet temperatures'
   radius          'inner radius should be the smallest'
   limit1          'limit exit loss to 0.1% of the pressure drop'
   limit2          'avoid surface damage in case of pressure loss'
   temp1           'calculation of tmp1'
   temp2           'calculation of tmp2';

Scalar
   gamma           'weight density of oil          (lb/in^3)' / 0.0307 /
   C               'specific heat of oil     (Btu/lb degr F)' / 0.5    /
   N               'angular speed of shaft             (rpm)' / 750    /
   Ws              'thrust load                         (lb)' / 101000 /
   Pmax            'maximum inlet pressure             (psi)' / 1000   /
   delta_t_max     'maximum temp change             (degr F)' / 50     /
   hmin            'minimum film thickness              (in)' / 0.001  /
   g               'gravitational constant          (in/s^2)' / 386.4  /
   taf             'ambient temperature             (degr F)' / 100.0  /
   P1              'outlet pressure, atmospheric (psi gauge)' / 0      /
   pump_efficiency '70% efficiency'                           / 0.7    /;

Table oil_constants(*,*) 'various oil grades'
                C1      n
   'SAE 5'   10.85  -3.91
   'SAE 10'  10.45  -3.72
   'SAE 20'  10.04  -3.55
   'SAE 30'   9.88  -3.48
   'SAE 40'   9.83  -3.46
   'SAE 50'   9.82  -3.44;

Parameter
   C1  'c factor for given oil grade'
   cn  'n factor for given oil grade'
   gr  'specific gravity of oil'
   tar 'ambient temperature (degr R)';

* calculate or extract parameters
C1  = oil_constants('SAE 20','C1');
cn  = oil_constants('SAE 20','n');
gr  = gamma/0.0361111;
tar = 459.7 + taf;

* model restrictions
*  P0 <= pmax
*  h  >= hmin
*  delta_t <= delta_t_max
*  w  >= Ws
P0.up      = pmax;
h.lo       = hmin;
delta_t.up = delta_t_max;
w.lo       = Ws;

* bounds to keep variables in reasonable ranges
Ep.lo =   1;
Ef.lo =   1;
P0.lo =   1;
R.lo  =   1;
R.up  =  16;
R0.lo =   1;
R0.up =  16;
Q.lo  =   1;
Q.up  =  16;
mu.lo =   1.0e-6;
mu.up =  16.0e-6;
t.lo  = 100;

* variable scaling
mu.scale = 1.0e-6;
h.scale  = hmin;
W.scale  = Ws;
PL.scale = 1.0e4;
Ep.scale = 1.0e4;
Ef.scale = 1.0e4;

* bounds to prevent multiplication by zero
tmp1.lo = 0.0001;
tmp2.lo = 0.01;

* Initial point. Some NLP solvers need a good starting point for this problem.
R.L       = 6;
R0.L      = 5;
mu.L      = 0.000006;
Q.L       = 3;
P0.L      = 1000.0000;
Ef.L      = 16000;
Ep.L      = 3000;
h.L       = 0.001;
delta_t.L = 50;
t.L       = 600;
W.L       = 101000;

power_loss..       PL   =e= Ep + Ef;

pumping_energy..   Ep   =e= Q*(P0 - P1)/pump_efficiency;

* friction..       Ef =e= sqr(2*pi*N/60)*[(2*pi*mu)/h]*(r**4 - r0**4)/4;
* is reformulated into:
friction..         Ef*h =e= sqr(2*pi*N/60)*[(2*pi*mu)]*(r**4 - r0**4)/4;

* temp_rise..      delta_t =e= Ef/(12*778*Q*gamma*C);
* is reformulated into:
temp_rise..        delta_t*(12*778*Q*gamma*C) =e= Ef;

* load_capacity..  W =e= [(pi*P0)/2]*tmp2/tmp1;
* is reformulated into:
load_capacity..    W*tmp1 =e= [(pi*P0)/2]*tmp2;

* inlet_pressure.. P0 =e= (6*mu*Q)/(pi*h**3)*tmp1;
* is reformulated into:
inlet_pressure..   P0*(pi*h**3) =e= (6*mu*Q)*tmp1;

oil_viscosity..    log10(8.112e6*mu + 0.8) =e= (T**cn)*(10**C1);

temperature..      T =e= 560 + delta_t/2;

radius..           R =g= R0;

limit1..           gamma*sqr(Q) =l= 0.001*g*sqr(2*pi*R*h)*P0;

limit2..           W    =l= 5000*pi*tmp2;

temp1..            tmp1 =e= log(R) - log(R0);

temp2..            tmp2 =e= sqr(R) - sqr(R0);

Model m / all /;

m.scaleOpt = 1;
m.optca = 1E-3;
m.optcr = 1E-3;

solve m minimizing pl using nlp;
