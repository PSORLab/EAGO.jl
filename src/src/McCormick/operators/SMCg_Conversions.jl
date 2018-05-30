###############################################################################
################## CONVERT THROUGH BINARY DEFINITIONS #########################
###############################################################################

############################ Addition ##########################################
# Addition (Matched Floats)
function +(x::SMCg{N,V,T},y::T) where {N,V,T<:AbstractFloat}
	return SMCg{N,V,T}(x.cc+y, x.cv+y, x.cc_grad, x.cv_grad, (x.Intv+y),x.cnst)
end
function +(x::T,y::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	return SMCg{N,V,T}(x+y.cc, x+y.cv, y.cc_grad, y.cv_grad, (x+y.Intv),y.cnst)
end
# Addition (Mixed Floats)
function +(x::SMCg{N,V,T},y::C) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	f::T = convert(T,y)
	return SMCg{N,V,T}(x.cc+f, x.cv+f, x.cc_grad, x.cv_grad, (x.Intv+f), x.cnst)
end
function +(y::C,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	f::T = convert(T,y)
	return SMCg{N,V,T}(x.cc+f, x.cv+f, x.cc_grad, x.cv_grad, (x.Intv+f),x.cnst)
end
# Addition (Integers)
function +(x::SMCg{N,V,T},y::X) where {N,V,T<:AbstractFloat,X<:Integer}
	f::T = convert(T,y)
	return SMCg{N,V,T}(x.cc+f, x.cv+f, x.cc_grad, x.cv_grad, (x.Intv+f),x.cnst)
end
function +(y::X,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,X<:Integer}
	f::T = convert(T,y)
	return SMCg{N,V,T}(x.cc+f, x.cv+f, x.cc_grad, x.cv_grad, (x.Intv+f),x.cnst)
end

############################ Subtraction #######################################
# Subtraction (Matched Floats)
function -(x::SMCg{N,V,T},c::T) where {N,V,T<:AbstractFloat}
	return SMCg{N,V,T}(x.cc-c, x.cv-c, x.cc_grad, x.cv_grad, (x.Intv-c),x.cnst)
end
function -(c::T,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	return SMCg{N,V,T}(c-x.cv, c-x.cc, x.cc_grad, x.cv_grad, (c-x.Intv),x.cnst)
end
# Subtraction (Mixed Floats)
function -(x::SMCg{N,V,T},c::C) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	f::T = convert(T,c)
	return SMCg{N,V,T}(x.cc-f, x.cv-f, x.cc_grad, x.cv_grad, (x.Intv-f),x.cnst)
end
function -(c::C,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	f::T = convert(T,c)
	return SMCg{N,V,T}(f-x.cv, f-x.cc, x.cc_grad, x.cv_grad, (f-x.Intv),x.cnst)
end
# Subtraction (Integers)
function -(x::SMCg{N,V,T},c::X) where {N,V,T<:AbstractFloat,X<:Integer}
	f::T = convert(T,c)
	return SMCg{N,V,T}(x.cc-c, x.cv-c, x.cc_grad, x.cv_grad, (x.Intv-c),x.cnst)
end
function -(c::X,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,X<:Integer}
	f::T = convert(T,c)
	return SMCg{N,V,T}(c-x.cv, c-x.cc, x.cc_grad, x.cv_grad, (c-x.Intv),x.cnst)
end

######################### Multiplication #######################################
# Multiplication (Matched Floats)
function *(x::SMCg{N,V,T},c::T) where {N,V,T<:AbstractFloat}
	if (c>=zero(T))
		return SMCg{N,V,T}(c*x.cc,c*x.cv,c*x.cc_grad,c*x.cv_grad,c*x.Intv,x.cnst)
	else
		return SMCg{N,V,T}(c*x.cv,c*x.cc,c*x.cv_grad,c*x.cc_grad,c*x.Intv,x.cnst)
	end
end
function *(c::T,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
	if (c>=zero(T))
		return SMCg{N,V,T}(c*x.cc,c*x.cv,c*x.cc_grad,c*x.cv_grad,c*x.Intv,x.cnst)
	else
		return SMCg{N,V,T}(c*x.cv,c*x.cc,c*x.cv_grad,c*x.cc_grad,c*x.Intv,x.cnst)
	end
end
# Multiplication (Mixed Floats)
function *(x::SMCg{N,V,T},c::C) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	if (c>=zero(C))
		return SMCg{N,V,T}(convert(T,c*x.cc),convert(T,c*x.cv),c*x.cc_grad,c*x.cv_grad,c*x.Intv,x.cnst)
	else
		return SMCg{N,V,T}(convert(T,c*x.cv),convert(T,c*x.cc),c*x.cv_grad,c*x.cc_grad,c*x.Intv,x.cnst)
	end
end
function *(c::C,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	if (c>=zero(C))
		return SMCg{N,V,T}(convert(T,c*x.cc),convert(T,c*x.cv),c*x.cc_grad,c*x.cv_grad,c*x.Intv,x.cnst)
	else
		return SMCg{N,V,T}(convert(T,c*x.cv),convert(T,c*x.cc),c*x.cv_grad,c*x.cc_grad,c*x.Intv,x.cnst)
	end
end
# Multiplication(Integers)
function *(x::SMCg{N,V,T},c::C) where {N,V,T<:AbstractFloat,C<:Integer}
	if (c>=zero(C))
		return SMCg{N,V,T}(convert(T,c*x.cc),convert(T,c*x.cv),c*x.cc_grad,c*x.cv_grad,c*x.Intv,x.cnst)
	else
		return SMCg{N,V,T}(convert(T,c*x.cv),convert(T,c*x.cc),c*x.cv_grad,c*x.cc_grad,c*x.Intv,x.cnst)
	end
end
function *(c::C,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,C<:Integer}
	if (c>=zero(C))
		return SMCg{N,V,T}(convert(T,c*x.cc),convert(T,c*x.cv),c*x.cc_grad,c*x.cv_grad,c*x.Intv,x.cnst)
	else
		return SMCg{N,V,T}(convert(T,c*x.cv),convert(T,c*x.cc),c*x.cv_grad,c*x.cc_grad,c*x.Intv,x.cnst)
	end
end

##########################     Division  #######################################
function /(x::SMCg{N,V,T},y::C) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	x*inv(y)
end
function /(x::C,y::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	x*inv(y)
end
function /(x::SMCg{N,V,T},y::C) where {N,V,T<:AbstractFloat,C<:Integer}
	x*inv(y)
end
function /(x::C,y::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,C<:Integer}
	x*inv(y)
end

########################## Maximization  #######################################
# Addition (Matched Floats)
max(c::T,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = max(x,c)
# Addition (Mixed Floats)
function max(x::SMCg{N,V,T},c::C) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	f::T = convert(T,c)
	return max(x,f)
end
function max(c::C,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	f::T = convert(T,c)
	return max(x,f)
end
# Addition (Integers)
function max(x::SMCg{N,V,T},c::C) where {N,V,T<:AbstractFloat,C<:Integer}
	f::T = convert(T,c)
	return max(x,f)
end
function max(c::C,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,C<:Integer}
	f::T = convert(T,c)
	return max(x,f)
end

########################## Minimization  #######################################
# REPLACE DEFINITION WITH UNIVARIANT MIN CALL EVENTUALLY
# Minimization (Matched Floats)
min(c::T,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat} = -max(-x,-c)
min(x::SMCg{N,V,T},c::T) where {N,V,T<:AbstractFloat} = -max(-x,-c)
# Addition (Mixed Floats)
function min(x::SMCg{N,V,T},c::C) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	f::T = convert(T,c)
	return -max(-x,-f)
end
function min(c::C,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,C<:AbstractFloat}
	f::T = convert(T,c)
	return -max(-x,-f)
end
# Addition (Integers)
function min(x::SMCg{N,V,T},c::C) where {N,V,T<:AbstractFloat,C<:Integer}
	f::T = convert(T,c)
	return -max(-x,-f)
end
function min(c::C,x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat,C<:Integer}
	f::T = convert(T,c)
	return -max(-x,-f)
end

###############################################################################
##################      CONVERT THROUGH PROMOTION     #########################
###############################################################################

function convert(::Type{SMCg{N,V,T}},x::S) where {S<:Integer,N,V,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          SMCg{N,V,T}(convert(T,x),convert(T,x),seed,seed,V(convert(V,x)),false)
end
function convert(::Type{SMCg{N,V,T}},x::S) where {S<:AbstractFloat,N,V,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          SMCg{N,V,T}(convert(T,x),convert(T,x),seed,seed,V(convert(V,x)),false)
end
function convert(::Type{SMCg{N,V,T}},x::S) where {S<:Interval,N,V,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          SMCg{N,V,T}(convert(T,x.hi),convert(T,x.lo),seed,seed,convert(V,x),false)
end

promote_rule(::Type{SMCg{N,V,T}}, ::Type{S}) where {S<:Integer,N,V,T<:AbstractFloat} = SMCg{N,V,T}
promote_rule(::Type{SMCg{N,V,T}}, ::Type{S}) where {S<:AbstractFloat,N,V,T<:AbstractFloat} = SMCg{N,V,T}
promote_rule(::Type{SMCg{N,V,T}}, ::Type{S}) where {S<:Interval,N,V,T<:AbstractFloat} = SMCg{N,V,T}
promote_rule(::Type{SMCg{N,V,T}}, ::Type{S}) where {S<:Real,N,V,T<:AbstractFloat} = SMCg{N,V,T}
