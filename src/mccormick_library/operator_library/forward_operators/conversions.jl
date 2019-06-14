###############################################################################
################## CONVERT THROUGH BINARY DEFINITIONS #########################
###############################################################################

############################ Addition ##########################################
# Addition (Matched Floats)
+(x::MC{N},y::Float64) where N = MC{N}(x.cv+y, x.cc+y, (x.Intv+y), x.cv_grad, x.cc_grad, x.cnst)
+(x::Float64,y::MC{N}) where N = MC{N}(x+y.cv, x+y.cc, (x+y.Intv), y.cv_grad, y.cc_grad, y.cnst)

# Addition (Mixed Floats)
function +(x::MC{N},y::C) where {N,C<:AbstractFloat}
	f::Float64 = convert(Float64,y)
	return MC{N}(x.cv+f, x.cc+f, (x.Intv+f), x.cv_grad, x.cc_grad, x.cnst)
end
function +(y::C,x::MC{N}) where {N,C<:AbstractFloat}
	f::Float64 = convert(Float64,y)
	return MC{N}(x.cv+f, x.cc+f, (x.Intv+f), x.cv_grad, x.cc_grad, x.cnst)
end

# Addition (Integers)
function +(x::MC{N},y::X) where {N,X<:Integer}
	f::Float64 = convert(Float64,y)
	return MC{N}(x.cv+f, x.cc+f, (x.Intv+f), x.cv_grad, x.cc_grad, x.cnst)
end
function +(y::X,x::MC{N}) where {N,X<:Integer}
	f::Float64 = convert(Float64,y)
	return MC{N}(x.cv+f, x.cc+f, (x.Intv+f), x.cv_grad, x.cc_grad, x.cnst)
end

############################ Subtraction #######################################
# Subtraction (Matched Floats)
-(x::MC{N},c::Float64) where N = MC{N}(x.cv-c, x.cc-c, (x.Intv-c), x.cv_grad, x.cc_grad, x.cnst)
-(c::Float64,x::MC{N}) where N = MC{N}(c-x.cc, c-x.cv, (c-x.Intv), x.cv_grad, x.cc_grad, x.cnst)

# Subtraction (Mixed Floats)
function -(x::MC{N},c::C) where {N,C<:AbstractFloat}
	f::Float64 = convert(Float64,c)
	return MC{N}(x.cv-f, x.cc-f, (x.Intv-f), x.cv_grad, x.cc_grad, x.cnst)
end
function -(c::C,x::MC{N}) where {N,C<:AbstractFloat}
	f::Float64 = convert(Float64,c)
	return MC{N}(f-x.cc, f-x.cv, (f-x.Intv), x.cv_grad, x.cc_grad, x.cnst)
end
# Subtraction (Integers)
function -(x::MC{N},c::X) where {N,X<:Integer}
	f::Float64 = convert(Float64,c)
	return MC{N}(x.cv-c, x.cc-c, (x.Intv-c), x.cv_grad, x.cc_grad, x.cnst)
end
function -(c::X,x::MC{N}) where {N,X<:Integer}
	f::Float64 = convert(Float64,c)
	return MC{N}(c-x.cc, c-x.cv, (c-x.Intv), x.cv_grad, x.cc_grad, x.cnst)
end

######################### Multiplication #######################################
# Multiplication (Matched Floats)
function *(x::MC{N},c::Float64) where N
	if (c>=zero(Float64))
		return MC{N}(c*x.cv,c*x.cc,c*x.Intv,c*x.cv_grad,c*x.cc_grad,x.cnst)
	else
		return MC{N}(c*x.cc,c*x.cv,c*x.Intv,c*x.cc_grad,c*x.cv_grad,x.cnst)
	end
end
function *(c::Float64,x::MC{N}) where N
	if (c>=zero(Float64))
		return MC{N}(c*x.cv,c*x.cc,c*x.Intv,c*x.cv_grad,c*x.cc_grad,x.cnst)
	else
		return MC{N}(c*x.cc,c*x.cv,c*x.Intv,c*x.cc_grad,c*x.cv_grad,x.cnst)
	end
end
# Multiplication (Mixed Floats)
function *(x::MC{N},c::C) where {N,C<:AbstractFloat}
	if (c>=zero(C))
		return MC{N}(convert(Float64,c*x.cv),convert(Float64,c*x.cc),c*x.Intv,c*x.cv_grad,c*x.cc_grad,x.cnst)
	else
		return MC{N}(convert(Float64,c*x.cc),convert(Float64,c*x.cv),c*x.Intv,c*x.cc_grad,c*x.cv_grad,x.cnst)
	end
end
function *(c::C,x::MC{N}) where {N,C<:AbstractFloat}
	if (c>=zero(C))
		return MC{N}(convert(Float64,c*x.cv),convert(Float64,c*x.cc),c*x.Intv,c*x.cv_grad,c*x.cc_grad,x.cnst)
	else
		return MC{N}(convert(Float64,c*x.cc),convert(Float64,c*x.cv),c*x.Intv,c*x.cc_grad,c*x.cv_grad,x.cnst)
	end
end
# Multiplication(Integers)
function *(x::MC{N},c::C) where {N,C<:Integer}
	if (c>=zero(C))
		return MC{N}(convert(Float64,c*x.cv),convert(Float64,c*x.cc),c*x.Intv,c*x.cv_grad,c*x.cc_grad,x.cnst)
	else
		return MC{N}(convert(Float64,c*x.cc),convert(Float64,c*x.cv),c*x.Intv,c*x.cc_grad,c*x.cv_grad,x.cnst)
	end
end
function *(c::C,x::MC{N}) where {N,C<:Integer}
	if (c>=zero(C))
		return MC{N}(convert(Float64,c*x.cv),convert(Float64,c*x.cc),c*x.Intv,c*x.cv_grad,c*x.cc_grad,x.cnst)
	else
		return MC{N}(convert(Float64,c*x.cc),convert(Float64,c*x.cv),c*x.Intv,c*x.cc_grad,c*x.cv_grad,x.cnst)
	end
end

##########################     Division  #######################################
/(x::MC,y::C) where {C<:AbstractFloat} = x*inv(y)
/(x::C,y::MC) where {C<:AbstractFloat} = x*inv(y)
/(x::MC,y::C) where {C<:Integer} = x*inv(y)
/(x::C,y::MC) where {C<:Integer} = x*inv(y)

########################## Maximization  #######################################
# Addition (Matched Floats)
max(c::Float64,x::MC) = max(x,c)

# Addition (Mixed Floats)
function max(x::MC,c::C) where {C<:AbstractFloat}
	f::Float64 = convert(Float64,c)
	return max(x,f)
end
function max(c::C,x::MC) where {C<:AbstractFloat}
	f::Float64 = convert(Float64,c)
	return max(x,f)
end
# Addition (Integers)
function max(x::MC,c::C) where {C<:Integer}
	f::Float64 = convert(Float64,c)
	return max(x,f)
end
function max(c::C,x::MC) where {C<:Integer}
	f::Float64 = convert(Float64,c)
	return max(x,f)
end

########################## Minimization  #######################################
# REPLACE DEFINITION WITH UNIVARIANT MIN CALL EVENTUALLY
# Minimization (Matched Floats)
min(c::Float64,x::MC) = -max(-x,-c)
min(x::MC,c::Float64) = -max(-x,-c)

# Addition (Mixed Floats)
function min(x::MC,c::C) where {C<:AbstractFloat}
	f::Float64 = convert(Float64,c)
	return -max(-x,-f)
end
function min(c::C,x::MC) where {C<:AbstractFloat}
	f::Float64 = convert(Float64,c)
	return -max(-x,-f)
end

# Addition (Integers)
function min(x::MC,c::C) where {C<:Integer}
	f::Float64 = convert(Float64,c)
	return -max(-x,-f)
end
function min(c::C,x::MC) where {C<:Integer}
	f::Float64 = convert(Float64,c)
	return -max(-x,-f)
end

###############################################################################
##################      CONVERT THROUGH PROMOTION     #########################
###############################################################################

promote_rule(::Type{MC{N}}, ::Type{S}) where {S<:Integer,N} = MC{N}
promote_rule(::Type{MC{N}}, ::Type{S}) where {S<:AbstractFloat,N} = MC{N}
promote_rule(::Type{MC{N}}, ::Type{S}) where {S<:Interval,N} = MC{N}
promote_rule(::Type{MC{N}}, ::Type{S}) where {S<:Real,N} = MC{N}

function convert(::Type{MC{N}},x::S) where {S<:Integer,N}
         MC{N}(IntervalConstr(x))
end
function convert(::Type{MC{N}},x::S) where {S<:AbstractFloat,N}
         MC{N}(IntervalConstr(x))
end
function convert(::Type{MC{N}},x::S) where {S<:Interval,N}
         MC{N}(IntervalType(x.lo,x.hi))
end
