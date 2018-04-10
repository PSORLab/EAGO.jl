a = MCInterval(2.0,3.0)
b = MCInterval(-0.5,0.5)

sinh(Interval(0.5))
@test sinh(Interval(0.5, 1.67))
@test sinh(Interval(-4.5, 0.1))

@test cosh(Interval(0.5)) == Interval(1.1276259652063807, 1.127625965206381)
@test cosh(Interval(0.5, 1.67)) == Interval(1.1276259652063807, 2.750207431409957)
@test cosh(Interval(-4.5, 0.1)) == Interval(1.0, 45.01412014853003)

@test tanh(emptyinterval()) == emptyinterval()
@test tanh(Interval(0.5)) == Interval(0.46211715726000974, 0.4621171572600098)
@test tanh(Interval(0.5, 1.67)) == Interval(0.46211715726000974, 0.9315516846152083)
@test tanh(Interval(-4.5, 0.1)) == Interval(-0.9997532108480276, 0.09966799462495583)

@test asinh(@biginterval(1)) ⊆ asinh(@interval(1))
@test asinh(@biginterval(0.9, 2)) ⊆ asinh(@interval(0.9, 2))
@test asinh(@biginterval(3, 4)) ⊆ asinh(@interval(3, 4))
