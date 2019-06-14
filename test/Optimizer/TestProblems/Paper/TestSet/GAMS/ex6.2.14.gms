*****************************************************************
*** Chapter 6
*** Test Problem 14
***
*** Methanol - Cyclohexane -- Gibbs energy minimization(Modified Wilson)
*****************************************************************

SETS
    i     components                /1*2/
    k     phases                    /1*2/
    alias(i,j)
    alias(i,l);

*****************************************************************
* P = pressure (atm)
* T = temperature (K)
* liqphase(k) = is phase k a liquid phase?
* ntot(i) = number of moles of component i in feed
* lambda(i,j), rho(i,j) = binary interaction parameters
* psi(i) = calculated parameter

PARAMETERS
gc, P, T, liqphase(k), ntot(i), lambda(i,j), rho(i,j), psi(i);

gc = 1.98721;
P = 1.0;
T = 298;

liqphase('1') = 1; liqphase('2') = 1;

ntot('1') = 0.5;
ntot('2') = 0.5;

lambda('1','1') = 1.0; lambda('1','2') = 0.30384;
lambda('2','1') = 0.095173; lambda('2','2') = 1.0;

rho('1','1') = 1.0; rho('1','2') = 0.374;
rho('2','1') = 2.6738; rho('2','2') = 1.0;

psi(i) = 1 + SUM(j$(ord(j) NE ord(i)), rho(j,i));
psi('1') = 3.6838;
psi('2') = 1.59549;

DISPLAY psi;

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


*obj.. gfe =e= SUM(k, SUM(i, n(i,k)*( 2.0*LOG(n(i,k)) - LOG(SUM(j,n(j,k)))
*				    -LOG(SUM(j,lambda(j,i)*n(j,k))) ) )
*                    +SUM(i, SUM(j,rho(j,i)*n(j,k))*LOG(SUM(j,rho(j,i)*n(j,k))) )
*		    +SUM(i, SUM(j$(ord(j) NE ord(i)), rho(j,i)*n(j,k)
*				     *LOG(n(j,k)/SUM(l,rho(l,i)*n(l,k))))) )
*	     -SUM(k, SUM(i, psi(i)*n(i,k)*LOG(n(i,k))));
obj.. gfe =e= SUM(k, SUM(i, n(i,k)*( LOG(n(i,k)/SUM(j,n(j,k)))
				    +LOG(n(i,k)/SUM(j,lambda(j,i)*n(j,k))) ) )
                    +SUM(i, SUM(j,rho(j,i)*n(j,k))*LOG(SUM(j,rho(j,i)*n(j,k))) )
		    +SUM(i, SUM(j$(ord(j) NE ord(i)), rho(j,i)*n(j,k)
				     *LOG(n(j,k)/SUM(l,rho(l,i)*n(l,k)))) ) )
	     +SUM(k, SUM(i, -psi(i)*n(i,k)*LOG(n(i,k))));

molesum(i).. SUM(k, n(i,k)) =e= ntot(i);

MODEL gmin / all /;

n.lo(i,k) = 0.0000001; n.up(i,k) = ntot(i);
n.l('1','1') = 0.0583; n.l('1','2') = 0.4417;
n.l('2','1') = 0.4080; n.l('2','2') = 0.0920;
*n.lo('1','1') = 0.0583; n.lo('1','2') = 0.4417;
*n.lo('2','1') = 0.4080; n.lo('2','2') = 0.0920;
*n.up('1','1') = 0.0583; n.up('1','2') = 0.4417;
*n.up('2','1') = 0.4080; n.up('2','2') = 0.0920;

SOLVE gmin USING nlp MINIMIZING gfe;
