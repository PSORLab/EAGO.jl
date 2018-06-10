workspace()

using Documenter, EAGO

makedocs(modules = [EAGO],
         format = :html,
         sitename = "EAGO.jl",
         pages = ["Introduction" => "intro.md",
                  "EAGO Solver" => ["Getting Started" => "EAGOSolver/starting.md",
                                    "Setting Solver Options" => "EAGOSolver/SolverOpts.md",
                                    "MathProgBase Inferace" => "EAGOSolver/MPB.md",
                                    "SIP Solver Interface" => "EAGOSolver/SIP.md",
                                    "High Performance Builds" => "EAGOSolver/hp.md"],
                  "McCormick Relaxation Library" => ["Overview" => "McCormick/overview.md",
                                                     "Bounding Functions via McCormick Operators" => "McCormick/Usage.md",
                                                     "Supported Operators" => "McCormick/Operators.md",
                                                     "Options available for McCormick Operators" => "McCormick/Options.md"],
                  "Branch and Bound Library" => ["Overview" => "BranchBound/Overview.md",
                                                 "Usage" => "BranchBound/usage.md",
                                                 "Types" => "BranchBound/types.md",
                                                 "API" => "BranchBound/API.md"],
                  "Developer Documentation" => ["Future Plans" => "Dev/future.md",
                                                "How to Contribute" => "Dev/contributing.md"],
                  "References" => "ref.md",
                  "Citing EAGO" => "cite.md"]
)

deploydocs(
    repo   = "github.com/PSORLab/EAGO.jl.git",
    osname = "linux",
    julia = "0.6",
    target = "build",
    deps = nothing,
    make = nothing)
