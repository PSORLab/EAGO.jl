a = MCInterval(2.0,3.0)
b = MCInterval(-0.5,0.5)

t1 = sin(a)
t2 = cos(a)
t3 = tan(a)
#t4 = asin(b)
#t5 = acos(b)
#t6 = atan(a)
#t7 = atan2(b,a)

@test tan(@interval(0.5)) == Interval(0.54630248984379048, 0.5463024898437906)
@test tan(@interval(0.5, 1.67)) == entireinterval()
@test tan(@interval(1.67, 3.2)) == Interval(-10.047182299210307, 0.05847385445957865)

@test asin(@interval(1)) == @interval(pi/2)#pi_interval(Float64)/2
@test asin(@interval(0.9, 2)) == asin(@interval(0.9, 1))
@test asin(@interval(3, 4)) == ∅

@test acos(@interval(1)) == Interval(0., 0.)
@test acos(@interval(-2, -0.9)) == acos(@interval(-1, -0.9))
@test acos(@interval(3, 4)) == ∅

@test atan(@interval(-1,1)) ==
    Interval(-pi_interval(Float64).hi/4, pi_interval(Float64).hi/4)
@test atan(@interval(0)) == Interval(0.0, 0.0)
