push!(LOAD_PATH,"../src/")

using Documenter, DocumenterTools
using EAGO, IntervalArithmetic, MathOptInterface, McCormick, JuMP

import McCormick: final_cut, mid3v, precond_and_contract!, AbstractMCCallback, populate_affine!,
                  dline_seg, correct_exp!, cut, mid_grad, preconditioner_storage, newton, secant,
                  MCCallback, contract!, affine_exp!

import EAGO: ExtensionType, Evaluator, variable_dbbt!, set_current_node!,
             VariableInfo, Log, aggressive_filtering!,
             bool_indx_diff, trivial_filtering!, SIPResult, SIPProblem, 
             GlobalOptimizer, InputProblem, ParsedProblem, is_integer_feasible_relaxed, 
             local_problem_status, default_upper_heuristic, label_branch_variables!,
             label_fixed_variables!, AbstractDirectedGraph, AbstractCache, 
             AbstractCacheAttribute, initialize!, f_init!, r_init!, fprop!, rprop!,
             Variable, Subexpression, Expression, Constant, Parameter
             
import EAGO.Script: dag_flattening!, register_substitution!, Template_Graph,
                    Template_Node, scrub, scrub!, flatten_expression!

const MOI = MathOptInterface

@info "Making documentation..."
makedocs(modules = [EAGO, McCormick],
         doctest = false,
         format = Documenter.HTML(
                prettyurls = get(ENV, "CI", nothing) == "true",
                canonical = "https://PSORLab.github.io/EAGO.jl/stable/",
                collapselevel = 1,
         ),
         authors = "Matthew Wilhelm, Robert Gottlieb, Dimitri Alston, and Matthew Stuber",
         #repo = "https://github.com/PSORLab/EAGO.jl/blob/{commit}{path}#L{line}",
         sitename = "EAGO",
         pages = Any["Introduction" => "index.md",
                     "Quick Start" => Any["quick_start/qs_landing.md",
                                          "quick_start/guidelines.md",
                                          "quick_start/simple.md",
                                          "quick_start/medium.md",
                                          "quick_start/difficult.md"
                                          ],
                     "McCormick Operator Library" => Any["mccormick/overview.md",
                                                         "mccormick/usage.md",
                                                         "mccormick/operators.md",
                                                         "mccormick/type.md",
                                                         "mccormick/implicit.md"
                                                        ],
                     "Optimizer" => Any["optimizer/optimizer.md",
                                        "optimizer/bnb_back.md",
                                        "optimizer/relax_back.md",
                                        "optimizer/domain_reduction.md",
                                        "optimizer/high_performance.md",
                                        "optimizer/udf_utilities.md"
                                               ],
                     "Semi-Infinite Programming" => "semiinfinite/semiinfinite.md",
                     "Contributing to EAGO"      => Any["dev/contributing.md",
                                                        "dev/future.md"
                                                        ],
                     "References"  => "ref.md",
                     "Citing EAGO" => "cite.md",
                     "News"        => "news.md"]
)

@info "Deploying documentation..."
deploydocs(repo = "github.com/PSORLab/EAGO.jl.git",
           push_preview  = true)
