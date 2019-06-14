*****************************************************************
*** Chapter 6
*** Test Problem 10
***
*** Ethanol - Benzene - Water -- Gibbs energy minimization(UNIFAC)
*****************************************************************

SETS
    i     components                /1*3/
    k     phases                    /1*2/
    m     groups                    /1*3/
    alias(i,j)
    alias(m,l);

*****************************************************************
* P = pressure (atm)
* T = temperature (K)
* liqphase(k) = is phase k a liquid phase?
* ntot(i) = number of moles of component i in feed
* z, zr(i), zrm, zb(i), za, psi(i), q(i), r(i), nu(i),
*    vhat(m,i), lambda(m,i) = pure component and calculated parameters
* v(m,i) = group-component matrix
* bigq(m), bigr(m) = group parameters

PARAMETERS
P, T, liqphase(k), ntot(i), gform(i,k),
minval,
z, bigq(m), bigr(m), v(m,i), zr(i), zrm, zb(i), za, psi(i), tau(m,l), bip(m,l),
q(i), r(i), vhat(m,i), temp(m,i), lambda(m,i), nu(i);

P = 1.0;
T = 298.0;

liqphase('1') = 1; liqphase('2') = 1;

gform('1','1') = -2.562164; gform('1','2') = -2.562164;
gform('2','1') = -2.084538; gform('2','2') = -2.084538;
gform('3','1') = -3.482138; gform('3','2') = -3.482138;

ntot('1') = 0.2;
ntot('2') = 0.4;
ntot('3') = 0.4;

z = 10.0;

bigq('1') = 1.972;
bigq('2') = 0.400;
bigq('3') = 1.400;

bigr('1') = 2.1055;
bigr('2') = 0.5313;
bigr('3') = 0.9200;

v('1','1') = 1; v('1','2') = 0; v('1','3') = 0;
v('2','1') = 0; v('2','2') = 6; v('2','3') = 0;
v('3','1') = 0; v('3','2') = 0; v('3','3') = 1;

tau('1','1') = 0.0; tau('1','2') = 89.6; tau('1','3') = 353.5;
tau('2','1') = 636.1; tau('2','2') = 0.0; tau('2','3') = 903.8;
tau('3','1') = -229.1; tau('3','2') = 362.3; tau('3','3') = 0.0;

bip(m,l) = EXP(-tau(m,l)/T);
q(i) = SUM(m, v(m,i)*bigq(m));
r(i) = SUM(m, v(m,i)*bigr(m));
zr(i) = (z*q(i)/2.0 - 1.0)/r(i);
minval = 100.0;
LOOP(i, IF( (zr(i) LT minval), minval = zr(i)));
zrm = minval;
za = zrm + SUM(i,zr(i)-zrm);
zb(i) = SUM(j$(ord(j) NE ord(i)), zr(j)-zrm);
psi(i) = q(i) + r(i)*za;

vhat(m,i) = SUM(l, bigq(l)*v(l,i)*bip(l,m));
temp(m,i) = 1-LOG(vhat(m,i)/q(i))-SUM(l,v(l,i)*bigq(l)*bip(m,l)/vhat(l,i));
lambda(m,i) = EXP(bigq(m)*temp(m,i));
*nu(i) = SUM(m, v(m,i)*lambda(m,i));
nu('1') = 0.0;
nu('2') = 0.0;
nu('3') = 0.0;


DISPLAY psi, zr, zb, za, bip, vhat, nu, temp, lambda;

***************************************************************************
* gfe = Gibbs free energy
* n(i,k) = number of moles of component i in phase k

VARIABLES gfe, n(i,k);

***************************************************************************
* obj = objective function
* molesum(i) = material balance on component i

EQUATIONS
    obj
    molesum(i);


obj.. gfe =e= SUM(k, SUM(i, n(i,k)*(gform(i,k)-zr(i)*r(i)*LOG(r(i))
				    + z*q(i)*LOG(q(i))/2.0 - nu(i)))
		    +SUM(i, za*r(i)*n(i,k))*LOG(SUM(j,r(j)*n(j,k)))
		    +SUM(i, zb(i)*r(i)*n(i,k)*LOG(n(i,k)/SUM(j, r(j)*n(j,k))))
		    +SUM(i, (z/2.0)*q(i)*n(i,k)*LOG(n(i,k)/SUM(j,q(j)*n(j,k))))
		    +SUM(i, q(i)*n(i,k))*LOG(SUM(j,q(j)*n(j,k)))
      +SUM(i, SUM(m, n(i,k)*v(m,i)*bigq(m)*LOG(n(i,k)/SUM(j,vhat(m,j)*n(j,k))))) )
	     +SUM(k,SUM(i, -psi(i)*n(i,k)*LOG(n(i,k))));

molesum(i).. SUM(k, n(i,k)) =e= ntot(i);

MODEL gmin / all /;

n.lo(i,k) = 0.0000001; n.up(i,k) = ntot(i);
n.l('1','1') = 0.19863; n.l('1','2') = 0.00137;
n.l('2','1') = 0.00428; n.l('2','2') = 0.39572;
n.l('3','1') = 0.39922; n.l('3','2') = 0.00078;

SOLVE gmin USING nlp MINIMIZING gfe;

