"""
    default_preprocess!

Runs interval, linear, and quadratic contractor methods up to tolerances
specified in `EAGO.Optimizer` object.
"""
function preprocess!(x::Optimizer,y::NodeBB)

    # Sets initial feasibility
    feas = true; rept = 0

    x.initial_volume = prod(upper_variable_bounds(y) - lower_variable_bounds(y))

    # runs poor man's LP contractor
    if (x.poor_man_lp_depth >= x.current_iteration_count)
        for i in 1:x.poor_man_lp_reptitions
            feas = poor_man_lp(x,y)
            (~feas) && (break)
        end
    end

    # runs univariate quadratic contractor
    if ((x.univariate_quadratic_depth >= x.current_iteration_count) && feas)
        for i in 1:x.univariate_quadratic_reptitions
            feas = univariate_quadratic(x,y)
            (~feas) && (break)
        end
    end
    #println("pre-obbt: $feas")

    x.obbt_performed_flag = false
    if ((x.obbt_depth >= x.current_iteration_count) && feas)
        x.obbt_performed_flag = true
        for i in 1:x.obbt_reptitions
            feas = obbt(x,y)
            (~feas) && (break)
        end
    end
    #println("obbt: $feas")

    if ((x.cp_depth >= x.current_iteration_count) && feas)
        feas = cpwalk(x,y)
    end
    #println("cp walk: $feas")

    x.final_volume = prod(upper_variable_bounds(y) - lower_variable_bounds(y))

    x.current_preprocess_info.feasibility = feas
end
