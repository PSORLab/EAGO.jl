# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_semiinfinite/algorithms/sip_hybrid.jl
# Defines the revised SIP-res algorithm which implements Algorithm #1 of Djelassi,
# Hatim, and Alexander Mitsos. "A hybrid discretization algorithm with guaranteed
# feasibility for the global solution of semi-infinite programs."
# Journal of Global Optimization 68.2 (2017): 227-253.
#############################################################################

struct SIPResRev <: AbstractSIPAlgo end
