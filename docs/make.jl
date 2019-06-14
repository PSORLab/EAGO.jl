push!(LOAD_PATH,"../src/")

using Documenter
using EAGO

makedocs(modules = [EAGO],
         doctest = true,
         format = Documenter.HTML(),
         authors = "Matthew Wilhelm",
         sitename = "EAGO.jl: Easy Advanced Global Optimization",
         pages = ["Introduction" => "index.md",
                  #"Why EAGO?" => "why.md",
                  "Branch & Bound" => ["BranchBound/overview.md",
                                       "BranchBound/usage.md",
                                       "BranchBound/back.md"],
                  "McCormick Operator Library" => ["McCormick/overview.md",
                                                   "McCormick/usage.md",
                                                   "McCormick/operators.md",
                                                   "McCormick/type.md",
                                                   "McCormick/implicit.md",
                                                   "McCormick/options.md"],
                  "Global Optimizer" => ["Optimizer/optimizers.md",
                                         "Optimizer/starting.md",
                                         "Optimizer/high_performance.md",
                                         "Optimizer/DomainReduction/feasibility.md",
                                         "Optimizer/DomainReduction/obbt.md"],
                  "Semi-Infinite Programming" => ["SemiInfinite/solving.md",
                                                  "SemiInfinite/api.md"],
                  "Contributing to EAGO" => ["Dev/contributing.md",
                                             "Dev/future.md"],
                  "References" => "ref.md",
                  "Citing EAGO" => "cite.md"]
)


deploydocs(
    deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo   = "github.com/PSORLab/EAGO.jl.git",
    target = "build",
    osname = "linux",
    julia = "1.1",
    deps = nothing,
    make = nothing)
