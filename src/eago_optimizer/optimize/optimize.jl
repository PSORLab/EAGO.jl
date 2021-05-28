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

include("optimize_lp_cone.jl")
include("optimize_convex.jl")
include("optimize_nonconvex.jl")

throw_optimize_hook!(m::Optimizer) = optimize_hook!(m.ext_type, m)

function MOI.optimize!(m::Optimizer)

    m._start_time = time()

    # Runs the branch and bound routine
    if !m.enable_optimize_hook

        # Need to eliminate fixed variables prior to checking
        # for convex quadratic constraints
        initial_parse!(m)

        # Determines if the problem is an LP, MILP, SOCP, MISCOP,
        # CONVEX, OF MINCVX PROBLEM TYPE
        parse_classify_problem!(m)

        m._parse_time = m._start_time - time()

        # Throws the problem to the appropriate solution routine
        optimize!(Val{m._working_problem._problem_type}(), m)
    else

        # throws to user-defined optimization hook
        throw_optimize_hook!(m)

    end

    return
end
