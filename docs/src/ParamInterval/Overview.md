## Parametric Interval methods

For those of you unfamiliar with basic interval arithmetic, please check out the
[ValidatedNumerics.jl](https://github.com/JuliaIntervals/ValidatedNumerics.jl) package.

Parametric interval contractor methods are quite similar to their standard interval
analogs. The key exception is that interval contractor methods, e.g. Newton/Krawczyk,
bound roots of a function, the parametric methods bound functions parametrically defined by these roots.

To elaborate further, the interval contractor method provides bounds, X, on $x \in X$ s.t. $f(x) = 0$.
The parameter interval contractor methods takes a function $y = x(p)$ defined implicitly by $h(y,p) = 0$ on (X,P)
and then subsequent computes a smaller box containing a solution branch of $h(y,p) = 0$.
