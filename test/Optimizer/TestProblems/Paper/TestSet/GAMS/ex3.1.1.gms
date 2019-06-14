*--------------------------------------------------------------*
* Quadratically Constrained - Test Problem 1                   *
*--------------------------------------------------------------*
File  
    res  / results / ;
Put res ;

Sets
    i    /1*8/;

Variables
    x(i)
    f;

x.lo('1') = 100;
x.up('1') = 10000;
x.lo('2') = 1000;
x.up('2') = 10000;
x.lo('3') = 1000;
x.up('3') = 10000;

Loop(i $ (ord(i) ge 4),
    x.lo(i) = 10;
    x.up(i) = 1000);

Equations
    Obj    objective function
    Con1   constraint function 1
    Con2   constraint function 2
    Con3   constraint function 3
    Con4   constraint function 4
    Con5   constraint function 5
    Con6   constraint function 6;

Obj .. 
    f =e= x('1') + x('2') + x('3');

Con1 ..
    -1 + 0.0025*(x('4') + x('6')) =l= 0;

Con2 ..
    -1 + 0.0025*(-x('4') + x('5') + x('7')) =l= 0;

Con3 ..
    -1 + 0.01*(-x('5') + x('8')) =l= 0;

Con4 ..
    100*x('1') - x('1')*x('6') + 833.33252*x('4') - 83333.333 =l= 0;

Con5 ..
    x('2')*x('4') - x('2')*x('7') - 1250*x('4') + 1250*x('5') =l= 0;

Con6 ..
    x('3')*x('5') - x('3')*x('8') - 2500*x('5') + 1250000 =l= 0;

Model
    problem /Obj, Con1, Con2, Con3, Con4, Con5, Con6/;

    x.l('1') = 579.19;
    x.l('2') = 1360.13;
    x.l('3') = 5109.92;
    x.l('4') = 182.01;
    x.l('5') = 295.60;
    x.l('6') = 217.99;
    x.l('7') = 286.40;
    x.l('8') = 395.60;

    solve problem using nlp minimizing f;
    PUT "Min f",f.l:16:10//;
    Loop(i, PUT "x   ",x.l(i):16:10//);
