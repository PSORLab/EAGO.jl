using MINLPTests

solver = JuMP.optimizer_with_attributes(EAGO.Optimizer, "relative_tolerance" => 1E-9)
 
minlp_nlp_exclude = String[
    "001_010", # Unbounded box, check solution bad if not gradient-based....
    "002_010", # Unbounded box   
    "004_010", # Unbounded box
    "004_011", # Unbounded box
    "005_010", # Unbounded box
    "006_010", # Unbounded box
    "007_010", # Unbounded box
    "008_010", # Unbounded box
    "008_011", # Unbounded box
    "005_011"  # \ operator not in JuMP
]
MINLPTests.test_nlp(solver, exclude = minlp_nlp_exclude, 
                            objective_tol = 1E-3,
                            #primal_tol = PRIMAL_TOL,
                            dual_tol = NaN,
                            termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
                            primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL)

minlp_nlp_cvx_exclude = String[

    "001_011", # convex quadratic objective... (linear unbounded...)
    "002_011", # unbounded linear problem & convex quadratic objective
    "101_010", # convex quadratic nl constraints...
    "101_011", 
    "101_012",
    "102_010",
    "102_011",
    "102_012",
    "102_013",
    "102_014",
    "103_010",
    "103_011",
    "103_012",
    "103_013",
    "103_014",
    "104_010",
    "105_010",
    "105_011",
    "105_012",
    "105_013",
    "106_010", # simple bounded domain
    "106_011", #
    "107_010",
    "107_011",
    "107_012",
    "108_010",
    "108_011",
    "108_012",
    "108_013",
    "109_010",
    "109_011",
    "109_012",
    "110_010",
    "110_011",
    "110_012",
    "201_010",
    "201_011",
    "202_010",
    "202_011",
    "202_012",
    "202_013",
    "202_014",
    "203_010",
    "204_010",
    "205_010",
    "206_010",
    "210_010",
    "210_011",
    "210_012",
    "501_010",
    "501_011"
]
MINLPTests.test_nlp_cvx(solver, exclude = minlp_nlp_cvx_exclude,
                                objective_tol = 1E-3,
                                #primal_tol = PRIMAL_TOL,
                                dual_tol = NaN,
                                termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
                                primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL)

minlp_nlp_mi_exclude = String[
    "001_010",  # no box constraints
    "002_010",

    #"003_010",   # TODO: Fix 003_010 - 003_016
    "003_011",  # FAIL
    "003_012",  # FAIL
    "003_013",  # FAIL
    "003_014",  # FAIL Never converges...
    "003_015",  #FAIL
    "003_016",

    "004_010",
    "004_011",
    "004_012",
    "005_010",
    "005_011",   # \ operator not in JuMP
    "006_010",
    "007_010",
    #"007_020"    # no way of specifying 
]
MINLPTests.test_nlp_mi(solver, exclude = minlp_nlp_mi_exclude, 
objective_tol = 1E-3,
termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL)