# Performs constraint walking on nonlinear terms
function cpwalk(x::Optimizer, n::NodeBB)

    # runs at midpoint bound
    midx = (n.upper_variable_bounds + n.lower_variable_bounds)/2.0

    # set working node to n, copies pass parameters from EAGO optimizer
    x.working_evaluator_block.evaluator.current_node = n
    x.working_evaluator_block.evaluator.has_reverse = true
    prior_sg_tighten = x.working_evaluator_block.evaluator.subgrad_tighten

    x.working_evaluator_block.evaluator.subgrad_tighten = false
    x.working_evaluator_block.evaluator.cp_reptitions = x.cp_interval_reptitions
    x.working_evaluator_block.evaluator.cp_tolerance = x.cp_interval_tolerance

    # Run forward-reverse pass and retreive node for interval forward-reverse pass
    feas = forward_reverse_pass(x.working_evaluator_block.evaluator, midx)
    n.lower_variable_bounds[:] = x.working_evaluator_block.evaluator.current_node.lower_variable_bounds
    n.upper_variable_bounds[:] = x.working_evaluator_block.evaluator.current_node.upper_variable_bounds


    # Run forward-reverse pass and retreive node for mccormick forward-reverse pass
    if feas
        x.working_evaluator_block.evaluator.subgrad_tighten = true
        x.working_evaluator_block.evaluator.cp_reptitions = x.cp_mccormick_reptitions
        x.working_evaluator_block.evaluator.cp_tolerance = x.cp_mccormick_tolerance

        feas = forward_reverse_pass(x.working_evaluator_block.evaluator, midx)
        n.lower_variable_bounds[:] = x.working_evaluator_block.evaluator.current_node.lower_variable_bounds
        n.upper_variable_bounds[:] = x.working_evaluator_block.evaluator.current_node.upper_variable_bounds
    end

    # resets forward reverse scheme for lower bounding problem
    x.working_evaluator_block.evaluator.has_reverse = false
    x.working_evaluator_block.evaluator.subgrad_tighten = prior_sg_tighten

    return feas
end
