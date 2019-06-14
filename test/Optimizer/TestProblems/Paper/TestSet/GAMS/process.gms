$title Alkylation Process Optimization (PROCESS,SEQ=20)

$onText
Optimization of a alkylation process.


Bracken, J, and McCormick, G P, Chapter 4. In Selected Applications
of Nonlinear Programming. John Wiley and Sons, New York, 1968.

Keywords: nonlinear programming, alkylation process, chemical engineering
$offText

Positive Variable
   olefin   'olefin feed                   (bpd)'
   isor     'isobutane recycle             (bpd)'
   acid     'acid addition rate (1000lb per day)'
   alkylate 'alkylate yield                (bpd)'
   isom     'isobutane makeup              (bpd)'
   strength 'acid strength          (weight pct)'
   octane   'motor octane number'
   ratio    'isobutane makeup to olefin ratio'
   dilute   'acid dilution factor'
   f4       'f-4 performance number';

Variable
   profit
   rangey
   rangem
   ranged
   rangef;

Equation
   yield    'alkylate yield definition'
   rngyield 'ranged alkylate yield definition'
   makeup   'isobutane makeup definition'
   sdef     'acid strength definition'
   motor    'motor octane number'
   rngmotor 'ranged motor octane number'
   drat     'isobutane to olefin ratio'
   ddil     'dilution definition'
   rngddil  'ranged dilution definition'
   df4      'f-4 definition'
   rngdf4   'ranged f-4 definition'
   dprofit  'profit definition';

yield..    alkylate =e= olefin*(1.12+.13167*ratio-.00667*sqr(ratio));

makeup..   alkylate =e= olefin + isom -.22*alkylate;

sdef..     acid   =e= alkylate*dilute*strength/(98-strength)/1000;

motor..    octane =e= 86.35 + 1.098*ratio - .038*sqr(ratio) - .325*(89-strength);

drat..     ratio  =e= (isor+isom)/olefin;

ddil..     dilute =e= 35.82 - .222*f4;

df4..      f4     =e= -133 + 3*octane;

dprofit..  profit =e= .063*alkylate*octane - 5.04*olefin - .035*isor - 10*acid - 3.36*isom;

rngyield.. rangey*alkylate =e= olefin*(1.12+.13167*ratio-.00667*sqr(ratio));

rngmotor.. rangem*octane =e= 86.35 + 1.098*ratio - .038*sqr(ratio) - .325*(89-strength);

rngddil..  ranged*dilute =e= 35.82 - .222*f4;

rngdf4..   rangef*f4 =e= -133 + 3*octane;

Model
   process 'process model with equalities'
           / yield, makeup, sdef, motor, drat, ddil, df4, dprofit /
   rproc   'ranged process  model'
           / rngyield, makeup, sdef, rngmotor, drat, rngddil, rngdf4, dprofit /;

rangey.lo = .9; rangey.up = 1.1; rangey.l = 1;
rangem.lo = .9; rangem.up = 1.1; rangem.l = 1;
ranged.lo = .9; ranged.up = 1.1; ranged.l = 1;
rangef.lo = .9; rangef.up = 1.1; rangef.l = 1;

strength.lo = 85;    strength.up = 93;  octane.lo   = 90;   octane.up = 95;
ratio.lo    = 3;     ratio.up    = 12;  dilute.lo   = 1.2;  dilute.up = 4;
f4.lo       = 145;   f4.up       = 162; olefin.lo   = 10;   olefin.up = 2000;
isor.up     = 16000; acid.up     = 120; alkylate.up = 5000; isom.up   = 2000;

olefin.l = 1745; isor.l     = 12000; acid.l   = 110;  alkylate.l = 3048;
isom.l   = 1974; strength.l = 89.2;  octane.l = 92.8; ratio.l    = 8;
dilute.l = 3.6;  f4.l       = 145;   profit.l = 872;

solve process maximizing profit using nlp;

solve rproc   maximizing profit using nlp;
