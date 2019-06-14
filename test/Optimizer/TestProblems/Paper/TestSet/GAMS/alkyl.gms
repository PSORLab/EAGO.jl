$title Simplified Alkylation Process (ALKYL,SEQ=165)

$onText
This model describes a simplified alkylation process. Note the
modeling of error bounds on the estimated equations. This formulation
is very efficient in terms of problem comprehension and solution. The
additional nonlinearities are bounded in a narrow range and
introduce no additional computational burden.


Berna, T, Locke, M, and Westerberg, A, Simplified Alkylation Process.
AIChE Journal 26 (1980), 37.

Keywords: nonlinear programming, chemical engineering, alkylation process
$offText

Variable
   F          'objective variable'
   OlefinFeed 'Olefins feed'
   IsobutRec  'Isobutane recycle'
   AcidFeed   'Acid feed'
   AlkylYld   'Alkylate yield'
   IsobutMak  'Isobutane makeup'
   AcidStren  'Acid strength'
   Octane     'Octane number'
   Ratio      'iC4 Olefin ratio'
   AcidDilut  'Acid dilution factor'
   F4Perf     'F4 performance number'
   alkerr
   octerr
   aciderr
   F4err;

Equation
   Objective  'objective function'
   AlkylShrnk 'Alkylate volumetric shrinkage equation'
   AcidBal    'Acid material balance'
   IsobutBal  'Isobutane component balance'
   AlkylDef
   OctDef
   AcidDef
   F4Def;

Objective..  F =e= - 6.3*AlkylYld*Octane + 5.04*OlefinFeed + 0.35*IsobutRec + AcidFeed + 3.36*IsobutMak;

AlkylShrnk.. AlkylYld =e= (OlefinFeed+IsobutMak)/1.22;

AcidBal..    0.98*AcidFeed =e= AcidStren*((AlkylYld*AcidDilut)/100.0 + AcidFeed);

IsoButBal..  10.0*IsobutRec + IsobutMak =e= OlefinFeed*Ratio;

AlkylDef..   AlkylYld*AlkErr   =e= OlefinFeed*(1.12 + 0.13167*Ratio - 0.0067*Ratio*Ratio);

OctDef..     Octane*OctErr     =e= 0.8635+(1.098*Ratio - 0.038*Ratio*Ratio)/100 + 0.325*(AcidStren - 0.89);

AcidDef..    AcidDilut*AcidErr =e= 35.82 - 22.2*F4Perf;

F4Def..      F4Perf*F4Err      =e= -1.33 + 3*Octane;

alkerr.lo  = .99; alkerr.up  = 1/.99; alkerr.l  = 1;
octerr.lo  = .99; octerr.up  = 1/.99; octerr.l  = 1;
aciderr.lo = .90; aciderr.up = 1/.90; aciderr.l = 1;
F4err.lo   = .99; F4err.up   = 1/.99; F4err.l   = 1;

F.l = -0.90;

OlefinFeed.lo = 0;    OlefinFeed.up = 2.00; OlefinFeed.l = 1.745;
IsobutRec.lo  = 0;    IsobutRec.up  = 1.60; IsobutRec.l  = 1.2;
AcidFeed.lo   = 0;    AcidFeed.up   = 1.20; AcidFeed.l   = 1.10;
AlkylYld.lo   = 0;    AlkylYld.up   = 5.00; AlkylYld.l   = 3.048;
IsobutMak.lo  = 0;    IsobutMak.up  = 2.00; IsobutMak.l  = 1.974;
AcidStren.lo  = 0.85; AcidStren.up  = 0.93; AcidStren.l  = 0.893;
Octane.lo     = 0.90; Octane.up     = 0.95; Octane.l     = 0.928;
Ratio.lo      = 3;    Ratio.up      = 12;   Ratio.l      = 8;
AcidDilut.lo  = 1.2;  AcidDilut.up  = 4;    AcidDilut.l  = 3.6;
F4Perf.lo     = 1.45; F4Perf.up     = 1.62; F4Perf.l     = 1.45;

Model m / all /;
m.optca = 1E-3;
m.optcr = 1E-3;

$onText optimal solutions
OlefinFeed.fx =  1.70368;
IsobutRec.fx  =  1.58449;
AcidFeed.fx   =   .543165;
AlkylYld.fx   =  3.03581;
IsobutMak.fx  =  2.0;
AcidStren.fx  =   .90133;
Octane.fx     =   .950;
Ratio.fx      = 10.4743;
AcidDilut.fx  =  1.56164;
F4Perf.fx     =  1.53535;
m.holdFixed   =  1;
$offText

solve m using nlp minimizing f;
