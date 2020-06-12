push!(LOAD_PATH,"../src/")

using Documenter
using EAGO, IntervalArithmetic, MathOptInterface, McCormick

import McCormick: final_cut, mid3v, precond_and_contract!, AbstractMCCallback, populate_affine!,
                  dline_seg, correct_exp!, cut, mid_grad, preconditioner_storage, newton, secant,
                  MCCallback, contract!, affine_exp!

import EAGO: ExtensionType, Evaluator, FunctionSetStorage, SubexpressionSetStorage, variable_dbbt!, set_current_node!,
             relax_quadratic_gen_saf, eval_objective_lo, eval_constraint_cc, eval_constraint_lo!, eval_constraint_hi!,
             eval_constraint_cc_grad, forward_reverse_pass, VariableInfo, Log, classify_quadratics!, lp_bound_tighten,
             univariate_kernel, univariate_quadratic, cpwalk, aggressive_filtering!, aggressive_obbt_on_heurestic,
             bool_indx_diff,  obbt, trivial_filtering!, SIPResult, SIPProblem, relax_convex_kernel
import EAGO.Script: dag_flattening!, register_substitution!, Template_Graph,
                    Template_Node, scrub, scrub!, flatten_expression!

const MOI = MathOptInterface

makedocs(modules = [EAGO, McCormick],
         doctest = false,
         format = Documenter.HTML(),
         authors = "Matthew Wilhelm",
         sitename = "EAGO.jl: Easy Advanced Global Optimization",
         pages = Any["Introduction" => "index.md",
                     "McCormick Operator Library" => Any["McCormick/overview.md",
                                                         "McCormick/usage.md",
                                                         "McCormick/operators.md",
                                                         "McCormick/type.md",
                                                         "McCormick/implicit.md"
                                                        ],
                     "Global Optimizer" => Any["Optimizer/starting.md",
                                               "Optimizer/optimizer.md",
                                               "Optimizer/bnb_back.md",
                                               "Optimizer/relax_back.md",
                                               "Optimizer/domain_reduction.md",
                                               "Optimizer/high_performance.md",
                                               "Optimizer/udf_utilities.md"
                                               ],
                     "Semi-Infinite Programming" => "SemiInfinite/semiinfinite.md",
                     "Contributing to EAGO"      => Any["Dev/contributing.md",
                                                        "Dev/future.md"
                                                        ],
                     "References"  => "ref.md",
                     "Citing EAGO" => "cite.md"]
)

deploydocs(repo = "github.com/PSORLab/EAGO.jl.git")
