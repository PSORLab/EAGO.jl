# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimize/optimize.jl
# Load subproblem type optimization routines and define general optimize!
# throw_optimize_hook!.
#############################################################################

include("optimize_lp.jl")
include("optimize_conic.jl")
include("optimize_convex.jl")
include("optimize_nonconvex.jl")

throw_optimize_hook!(m::Optimizer{Q,S,T}) where {Q,S,T}  = optimize_hook!(_ext_type(m), m)

function MOI.optimize!(m::Optimizer{Q,S,T}) where {Q,S,T}

    m._global_optimizer._start_time = time()

    # Runs the branch and bound routine
    if !m.enable_optimize_hook
        #println("input min sense 1 = $(m._input_problem._optimization_sense)")

        # Need to eliminate fixed variables prior to checking
        # for convex quadratic constraints
        initial_parse!(m)
        #println("input min sense initial parse = $(m._input_problem._optimization_sense)")

        # Determines if the problem is an LP, MILP, SOCP, MISCOP,
        # CONVEX, OF MINCVX PROBLEM TYPE
        parse_classify_problem!(m._global_optimizer)
        #println("input min sense post classify = $(m._input_problem._optimization_sense)")

        m._global_optimizer._parse_time = m._global_optimizer._start_time - time()

        # Throws the problem to the appropriate solution routine
        optimize!(m._working_problem._problem_type, m)
    else

        # throws to user-defined optimization hook
        throw_optimize_hook!(m)

    end

    return
end
