push!(LOAD_PATH,"../src/")

using Documenter, DocumenterTools
using EAGO, IntervalArithmetic, MathOptInterface, McCormick, JuMP

import McCormick: final_cut, mid3v, precond_and_contract!, AbstractMCCallback, populate_affine!,
                  dline_seg, correct_exp!, cut, mid_grad, preconditioner_storage, newton, secant,
                  MCCallback, contract!, affine_exp!

import EAGO: ExtensionType, Evaluator, variable_dbbt!,
             VariableInfo, Log, aggressive_filtering!,
             bool_indx_diff!, trivial_filtering!, SIPResult, SIPProblem, 
             GlobalOptimizer, InputProblem, ParsedProblem, is_integer_feasible_relaxed, 
             local_problem_status, default_upper_heuristic, label_branch_variables!,
             label_fixed_variables!, AbstractDirectedGraph, AbstractCache, 
             AbstractCacheAttribute, initialize!, f_init!, r_init!, fprop!, rprop!,
             Variable, Subexpression, Expression, Constant, Parameter
             
import EAGO.Script: dag_flattening!, register_substitution!, Template_Graph,
                    Template_Node, scrub, scrub!, flatten_expression!

@info "Making documentation..."
makedocs(modules = [EAGO, McCormick],
         doctest = false,
         warnonly = [:docs_block, :missing_docs],
         format = Documenter.HTML(
                  prettyurls = get(ENV, "CI", nothing) == "true",
                  canonical = "https://PSORLab.github.io/EAGO.jl/stable/",
                  collapselevel = 1,
                  assets = ["assets/favicon.ico"]
         ),
         authors = "Matthew Wilhelm, Robert Gottlieb, Dimitri Alston, and Matthew Stuber",
         sitename = "EAGO.jl",
         pages = Any["Introduction" => "index.md",
                     "Manual" => Any["Optimizer" => Any["optimizer/optimizer.md",
                                                        "optimizer/bnb_back.md",
                                                        "optimizer/relax_back.md",
                                                        "optimizer/domain_reduction.md",
                                                        "optimizer/high_performance.md",
                                                        "optimizer/udf_utilities.md"
                                                        ],
                                      "McCormick.jl" => Any["mccormick/overview.md",
                                                            "mccormick/usage.md",
                                                            "mccormick/operators.md",
                                                            "mccormick/type.md",
                                                            "mccormick/implicit.md"
                                                            ],
                                      "Semi-Infinite Programming" => "semiinfinite/semiinfinite.md",
                                     ],
                     "Customization" => "custom_guidelines.md",
                     "Examples" => Any["examples/explicit_ann.md",
                                       "examples/interval_bb.md",
                                       "examples/quasiconvex.md",
                                       "examples/alpha_bb.md",
                                       "examples/modelingtoolkit.md"
                                       ],
                     "API Reference" => Any["dev/api_types.md",
                                            "dev/api_functions.md"
                                            ],
                     "Contributing" => "dev/contributing.md",
                     "News" => "news.md",
                     "Citing EAGO" => "cite.md",
                     "References" => "ref.md"
                     ]
)

@info "Deploying documentation..."
deploydocs(repo = "github.com/PSORLab/EAGO.jl.git",
           push_preview  = true)
