push!(LOAD_PATH,"../src/")

using Documenter
using EAGO, IntervalArithmetic, MathOptInterface
using EAGO.McCormick
import EAGO.McCormick: affine_exp!, correct_exp!, final_cut, mc_dense_krawczyk_cw!,
                       mc_dense_newton_gs!, pmc_kernel!
import EAGO: ExtensionType, Evaluator, FunctionSetStorage,
             SubexpressionSetStorage, variable_dbbt!,
             set_current_node!, relax_quadratic_gen_saf, eval_objective_lo,
             eval_constraint_cc, eval_constraint_lo!, eval_constraint_hi!,
             eval_constraint_cc_grad, get_node_lower, get_node_upper,
             forward_reverse_pass, VariableInfo, Log, classify_quadratics!,
             lp_bound_tighten, univariate_kernel, univariate_quadratic, cpwalk,
             aggressive_filtering!, aggressive_obbt_on_heurestic, bool_indx_diff,
             obbt, trivial_filtering!, SIPResult, SIPProblem
import EAGO.Script: dag_flattening!, register_substitution!, Template_Graph,
                    Template_Node, scrub, scrub!, flatten_expression!

const MOI = MathOptInterface

makedocs(modules = [EAGO, McCormick],
         doctest = true,
         format = Documenter.HTML(),
         authors = "Matthew Wilhelm",
         sitename = "EAGO.jl: Easy Advanced Global Optimization",
         pages = ["Introduction" => "index.md",
                  "McCormick Operator Library" => ["McCormick/overview.md",
                                                   "McCormick/usage.md",
                                                   "McCormick/operators.md",
                                                   "McCormick/type.md",
                                                   "McCormick/implicit.md"],
                  "Global Optimizer" => ["Optimizer/starting.md",
                                         "Optimizer/optimizer.md",
                                         "Optimizer/bnb_back.md",
                                         "Optimizer/relax_back.md",
                                         "Optimizer/domain_reduction.md",
                                         "Optimizer/high_performance.md",
                                         "Optimizer/udf_utilities.md"],
                  "Semi-Infinite Programming" => "SemiInfinite/semiinfinite.md",
                  "Contributing to EAGO" => ["Dev/contributing.md",
                                             "Dev/future.md"],
                  "References" => "ref.md",
                  "Citing EAGO" => "cite.md"]
)


deploydocs(
    repo   = "github.com/PSORLab/EAGO.jl.git",
    target = "build",
    osname = "linux",
    julia = "1.1",
    deps = nothing,
    make = nothing
)
