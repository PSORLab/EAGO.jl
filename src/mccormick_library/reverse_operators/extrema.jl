# Currently no refinement on min/max... will add later

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `max`.
"""
max_rev(a::MC, b::MC, c::MC) = (a,b,c)
#max_rev(a,b,c) = max_rev(promote(a,b,c)...)

"""
$(FUNCTIONNAME)

Reverse McCormick operator for `min`.
"""
min_rev(a::MC, b::MC, c::MC) = (a,b,c)
#min_rev(a,b,c) = min_rev(promote(a,b,c)...)
