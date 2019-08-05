"""
    relax_model!

Takes an NLP data block structure, linear, quadratic constraints, and
variable bounds and subsequently builds the relaxed model
"""
function relax_model!(src::Optimizer, trg, n::NodeBB, r::RelaxationScheme, xpnt::Vector{Float64}; load::Bool = false)

    if load
        #println("during load")
        # add linear terms to model
        relax_linear!(src,trg)
        if ~isa(src.nlp_data.evaluator, EAGO.EmptyNLPEvaluator)
            # copy working evaluator into block if nonlinear block is needed
            if (r.optimizer_type == :NLP)
                if ~isempty(src.bisection_variable)
                    trg.nlp_data = src.working_evaluator_block
                end
            end
        end
    else
        set_current_node!(src.working_evaluator_block.evaluator, n)
        ~isinf(src.global_upper_bound) && objective_cut_linear!(src,trg)
        relax_quadratic!(trg,src,n,r)
        if ~isempty(src.bisection_variable)
            if MOI.supports(trg, MOI.NLPBlock())
                nlp_data = MOI.NLPBlockData(src.nlp_data.constraint_bounds,
                                            src.working_evaluator_block.evaluator,
                                            src.nlp_data.has_objective)
                MOI.set(trg, MOI.NLPBlock(), nlp_data)
            else
                midpoint_affine!(src,trg,n,r,xpnt)
            end
        end
    end
end
