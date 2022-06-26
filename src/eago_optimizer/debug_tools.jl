"""
print_problem_summary

Internal script used to display all constraints, objectives in a linear program.
Added to functions for debug purposes while writing code.
"""

function print_problem_summary!(m, s, b = true)
    println(" ")
    println(" ----------------------------- ")
    println(s)
    println(" ----------------------------- ")
    print(m)
    b && println(" ----------------------------- ")
end

function print_problem_summary!(d, m, s::String)
    
    print_problem_summary!(m, s, false)

    t_status = MOI.get(m, MOI.TerminationStatus())
    p_status = MOI.get(m, MOI.PrimalStatus())
    d_status = MOI.get(m, MOI.DualStatus())

    obj_val = MOI.get(m, MOI.ObjectiveValue())
    obj_bnd = MOI.get(m, MOI.ObjectiveBound())

    println("Termination Status = $(t_status)")
    println("Primal Status      = $(p_status)")
    println("Dual Status        = $(d_status)")

    println("Objective Value = $(obj_val)")
    println("Objective Bound = $(obj_bnd)")

    println("Solution at ")
    for i = 1:_variable_num(FullVar(), d)
        x = MOI.get(m, MOI.VariablePrimal(), d._relaxed_variable_index[i])
        println("x[$i] = $x")
    end
    println(" ----------------------------- ")
end