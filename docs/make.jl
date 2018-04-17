workspace()

using Documenter, EAGO

makedocs(modules = [EAGO],
         format = :html,
         sitename = "EAGO.jl",
         pages = ["Introduction" => "intro.md",
                  "EAGO Solver" => ["Getting Started" => "EAGOSolver/starting.md",
                                    "JuMP Interface" => "EAGOSolver/JuMP.md",
                                    "MathProgBase Inferace" => "EAGOSolver/MPB.md",
                                    "Setting Solver Options" => "EAGOSolver/SolverOpts.md",
                                    "Implicit Function Handling" => "EAGOSolver/Implicit.md",
                                    "SIP Solver Interface" => "EAGOSolver/SIP.md",
                                    "High Performance Builds" = > "EAGOSolver/hp.md"],
                  "McCormick Relaxation Library" => ["Overview" => "McCormick/overview.md",
                                                     "Bounding Functions via McCormick Operators" => "McCormick/McCormick.md",
                                                     "Options available for McCormick Operators" => "McCormick/Options.md"],
                  "Parametric Interval Methods" => ["Parametric Contractor" => "ParamInterval/contractor.md",
                                                    "Parametric Tests" => "ParamInterval/test.md"],
                  "Domain Reduction Library" => ["Overview" => "DomainReduction/overview.md",
                                                 "Range Reduction" => "DomainReduction/RR.md",
                                                 "Duality Based Bound Tightening" => "DomainReduction/DBBT.md",
                                                 "Probing" => "DomainReduction/probing.md"],
                  "Branch and Bound Library" => ["Overview" => "BranchBound/overview.md",
                                                 "Usage" => "BranchBound/usage.md",
                                                 "Types" => "BranchBound/types.md",
                                                 "API" => "BranchBound/API.md"],
                  "Developer Documentation" => ["Future Plans" => "Dev/future.md",
                                                "How to Contribute" => "Dev/contributing.md"],
                  "References" => "ref.md",
                  "Citing EAGO" => "cite.md"]
)
