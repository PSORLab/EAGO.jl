*--------------------------------------------------------------*
* Test Problem 1 in Chapter 8, section 5.3                     *
* Linear Model                                                 *
*--------------------------------------------------------------*

Sets
    m  number of data sets      /1*10/
    n  number of variables      /1*2/
    p  number of parameters     /1*2/;

Parameters
   ze(m,n)  observed data values
   std(n)   standard deviations;

Variables
   z(m,n)    fitted data variables
   a(p)      model paramters
   c         cost function;

table ze(m,n)
      1      2
1    0.0    5.9
2    0.9    5.4
3    1.8    4.4
4    2.6    4.6
5    3.3    3.5
6    4.4    3.7
7    5.2    2.8
8    6.1    2.8
9    6.5    2.4
10   7.4    1.5;

std(n) = 1;


Equations
      obj        objective function
      con(m)     non-convex constraint equality;


obj.. c =e= sum(m,sum(n,sqr((z(m,n)-ze(m,n))/std(n))));

con(m) .. - z(m,'2') + a('1') + a('2')*z(m,'1') =e= 0;


model
     problem /obj,con/;


z.lo(m,n) = ze(m,n) - 0.5;
z.up(m,n) = ze(m,n) + 0.5;

a.lo('1') = 0;
a.lo('2') = -2;

a.up('1') = 10;
a.up('2') = 2;

z.l(m,n)  = uniform(z.lo(m,n),z.up(m,n));
a.l(p) = uniform(a.lo(p), a.up(p));


solve problem using nlp minimizing c;



