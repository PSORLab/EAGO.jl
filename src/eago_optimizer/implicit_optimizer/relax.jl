"""
    implicit_relax_model!

Takes an NLP data block structure, linear, quadratic constraints, and
variable bounds and subsequently builds the relaxed implicit model.
"""
function implicit_relax_model!(src::Optimizer, trg, n::NodeBB, r::RelaxationScheme, x::Vector{Float64}; load::Bool = false)
    if load
        # add linear terms to model
        relax_linear!(src, trg)
        # build NLP evaluator and save to EAGO object
        if ~isa(src.nlp_data.evaluator, EAGO.EmptyNLPEvaluator)
            # copy working evaluator into block if nonlinear block is needed
            if (r.optimizer_type == :NLP)
                if ~isempty(src.NonlinearVariable)
                    trg.nlp_data = src.working_evaluator_block
                end
            end
        end
    else
        relax_quadratic!(trg,src,n,r)
        if ~isempty(src.nonlinear_variable)
            if MOI.supports(trg, MOI.NLPBlock())
                evaluator = src.working_evaluator_block.evaluator
                set_current_node!(evaluator, n)
                nlp_data = MOI.NLPBlockData(src.nlp_data.constraint_bounds,
                                            evaluator,
                                            src.nlp_data.has_objective)
                MOI.set(trg, MOI.NLPBlock(), nlp_data)
            else
                midpoint_affine!(src, trg, n, r, x)
            end
        end
    end
end
