"""
    MathProgBase.optimize!(s::EAGO_NLP_Model)

Optimizes the `s::EAGO_NLP_Model`. May print console outputs depending on settings
and other solution information is accessible via start solver interface functions.
"""
function MathProgBase.optimize!(m::EAGO_NLP_Model)

    #println("Start Optimization")
    #  Implicit solver load
    unshift!(m.Opts.solver.BnBSolver.opt,m.Opts)

    m.Opts.Imp_np = m.Opts.numVar - m.Opts.Imp_nx
    m.Opts.solver.PSmcOpt.np = m.Opts.Imp_np
    m.Opts.solver.PSmcOpt.nx = m.Opts.Imp_nx
    m.Opts.solver.PIntOpt.nx = m.Opts.Imp_nx

    # sets the McCormick library options as necessary
    if (m.Opts.solver.LBD_func_relax == "NS-STD")
        EAGO.set_diff_relax(0)
        EAGO.set_multivar_refine(false,1E-15)
    elseif (m.Opts.solver.LBD_func_relax == "NS-MV")
        EAGO.set_diff_relax(0)
        EAGO.set_multivar_refine(true,1E-15)
    elseif (m.Opts.solver.LBD_func_relax == "Diff1-MV")
        EAGO.set_diff_relax(1)
        EAGO.set_multivar_refine(true,1E-15)
    elseif (m.Opts.solver.LBD_func_relax == "Diff2-MV")
        EAGO.set_diff_relax(2)
        EAGO.set_multivar_refine(true,1E-15)
    end

    solveBnB!(m.Opts.solver.BnBSolver,m.BnBModel)

    if (m.Opts.sense == :Max)
        m.BnBModel.soln_val = -m.BnBModel.soln_val
    end

    if (m.BnBModel.feas_fnd == true)
        m.status = :Optimal
    else
        m.status = :Infeasible
    end
    #println("Finished Optimization Routine")
end
