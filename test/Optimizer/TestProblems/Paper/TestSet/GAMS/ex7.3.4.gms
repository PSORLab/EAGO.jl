* Generalized geometric programming
* Robust stability problem
* Vicino et al. (1990)

VARIABLES
	 q1  
	 q2  
	 q3  
	 q4  
	 q5  
	 a0
	 a1
	 a2
	 a3
	 a4
	 w     frequency
	 k     stability margin
     objval    objective function variable;

FREE VARIABLES    objval;


EQUATIONS
	 f Objective function
	 g1
	 g2
	 b1l  
	 b1u  
	 b2u  
	 b2l  
	 b3l  
	 b3u  
	 b4l  
	 b4u  
	 b5l  
	 b5u  
	 ga0  
	 ga1  
	 ga2  
	 ga3  
	 ga4  ;

f  .. objval =e=k;
g1  .. a4*POWER(w,4) - a2*POWER(w,2) + a0 =e= 0;
g2  .. a3*POWER(w,2) - a1 =e= 0;
b1l  .. -k - q1 =l= -10;
b1u  .. -k + q1 =l= 10;
b2l  .. -0.1*k - q2 =l= -1;
b2u  .. -0.1*k + q2 =l= 1;
b3l  .. -0.1*k - q3 =l= -1;
b3u  .. -0.1*k + q3 =l= 1;
b4l  .. -0.01*k - q4 =l= -0.2;
b4u  .. -0.01*k + q4 =l= 0.2;
b5l  .. -0.005*k - q5 =l= -0.05;
b5u  .. -0.005*k + q5 =l= 0.05;
ga0  .. a0 =e= 54.387*q3*q2;
ga1  .. a1 =e= 0.2*(-147.15*q4*q3*q2 + 1364.67*q3*q2 - 27.72*q5);
ga2  .. a2 =e= 3*(-9.81*q3*POWER(q2,2) - 9.81*q3*q1*q2 - 4.312*POWER(q3,2)*q2 + 264.896*q3*q2 + q4*q5 - 9.274*q5);
ga3  .. a3 =e= 7*q4*POWER(q3,2)*q2 - 64.918*POWER(q3,2)*q2 + 380.067*q3*q2 + 3*q5*q2 + 3*q5*q1;
ga4  .. a4 =e= POWER(q3,2)*q2*(4*q2 + 7*q1);

* Bounds
w.LO = 0;
w.UP = 10;

* Starting point (global solution)
* q1.L = 16.2746;
* q2.L = 1.6275;
* q3.L = 1.675;
* q4.L = 0.1373;
* q5.L = 0.0186;
* w.L = 0.9864;
* k.L = 6.2746;

MODEL test /ALL/;
SOLVE test USING NLP MINIMIZING objval;
