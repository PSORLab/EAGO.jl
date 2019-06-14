# Defines functions required for linear algebra packages
one(x::MC{N}) where N = MC{N}(1.0,1.0,IntervalConstr(1.0),zeros(SVector{N,Float64}),zeros(SVector{N,Float64}),x.cnst)
zero(x::MC{N}) where N = MC{N}(0.0,0.0,IntervalConstr(0.0),zeros(SVector{N,Float64}),zeros(SVector{N,Float64}),x.cnst)
real(x::MC) = x
dist(x1::MC, x2::MC) = max(abs(x1.cc-x2.cc), abs(x1.cv-x2.cv))
eps(x::MC) = max(eps(x.cc), eps(x.cv))
mid(x::MC) = mid(x.Intv)

# Addition and subtraction of McCormick objects
+(x::MC{N},y::MC{N}) where N = MC{N}(x.cv+y.cv, x.cc+y.cc, (x.Intv+y.Intv),
                                     x.cv_grad+y.cv_grad, x.cc_grad+y.cc_grad,
                                     (x.cnst && y.cnst))

-(x::MC{N}) where N = MC{N}(-x.cc, -x.cv, -x.Intv, -x.cc_grad, -x.cv_grad, x.cnst)

-(x::MC{N},y::MC{N}) where N = MC{N}(x.cv-y.cc, x.cc-y.cv, (x.Intv-y.Intv),
                                     x.cv_grad-y.cc_grad, x.cc_grad-y.cv_grad,
                                     (x.cnst && y.cnst))
