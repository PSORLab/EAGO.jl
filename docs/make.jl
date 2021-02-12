push!(LOAD_PATH,"../src/")

using Documenter
using EAGO, IntervalArithmetic, MathOptInterface, McCormick

import McCormick: final_cut, mid3v, precond_and_contract!, AbstractMCCallback, populate_affine!,
                  dline_seg, correct_exp!, cut, mid_grad, preconditioner_storage, newton, secant,
                  MCCallback, contract!, affine_exp!

import EAGO: ExtensionType, Evaluator, variable_dbbt!, set_current_node!,
             VariableInfo, Log, aggressive_filtering!,
             bool_indx_diff, trivial_filtering!, SIPResult, SIPProblem
import EAGO.Script: dag_flattening!, register_substitution!, Template_Graph,
                    Template_Node, scrub, scrub!, flatten_expression!

const MOI = MathOptInterface

makedocs(modules = [EAGO, McCormick],
         doctest = false,
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")),
         authors = "Matthew E. Wilhelm",
         #repo = "https://github.com/PSORLab/EAGO.jl/blob/{commit}{path}#L{line}",
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

deploydocs(repo = "github.com/PSORLab/EAGO.jl",
           target = "build",
           deps   = nothing,
           make   = nothing)
