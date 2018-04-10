# List of valid function relaxations
relax_list = ["NS-MV-ON", "NS-MV-OFF", "NS-STD-ON", "NS-STD-OFF",
              "Diff1-MV-ON", "Diff1-MV-OFF", "Diff2-MV-ON", "Diff2-MV-OFF",
              "Interval","AlphaBB"]
relax_list_UBD = ["NS-MV-ON", "NS-MV-OFF", "NS-STD-ON", "NS-STD-OFF",
                  "Diff1-MV-ON", "Diff1-MV-OFF", "Diff2-MV-ON",
                  "Diff2-MV-OFF", "Original","Interval","AlphaBB"]

# Defines supported problem relaxations for given function relaxation
const LBD_f_relax = Dict{Any,Any}()
    LBD_f_relax["NS-MV-ON"] = ["LP","BNLP"]
    LBD_f_relax["NS-MV-OFF"] = ["LP","BNLP"]
    LBD_f_relax["NS-STD-ON"] = ["LP","BNLP"]
    LBD_f_relax["NS-STD-OFF"] = ["LP","BNLP"]
    LBD_f_relax["Diff1-MV-ON"] = ["LP","NLP1","BNLP"]
    LBD_f_relax["Diff1-MV-OFF"] = ["LP","NLP1","BNLP"]
    LBD_f_relax["Diff2-MV-ON"] = ["LP","QP","QCQP","NLP1","NLP2","BNLP"]
    LBD_f_relax["Diff2-MV-OFF"] = ["LP","QP","QCQP","NLP1","NLP2","BNLP"]
    LBD_f_relax["Interval"] = ["Interval"]
    LBD_f_relax["AlphaBB"] = ["AlphaBB"]
const UBD_f_relax = Dict{Any,Any}()
    UBD_f_relax["NS-MV-ON"] = ["LP","BNLP"]
    UBD_f_relax["NS-MV-OFF"] = ["LP","BNLP"]
    UBD_f_relax["NS-STD-ON"] = ["LP","BNLP"]
    UBD_f_relax["NS-STD-OFF"] = ["LP","BNLP"]
    UBD_f_relax["Diff1-MV-ON"] = ["LP","NLP1","BNLP"]
    UBD_f_relax["Diff1-MV-OFF"] = ["LP","NLP1","BNLP"]
    UBD_f_relax["Diff2-MV-ON"] = ["LP","QP","QCQP","NLP1","NLP2","BNLP"]
    UBD_f_relax["Diff2-MV-OFF"] = ["LP","QP","QCQP","NLP1","NLP2","BNLP"]
    UBD_f_relax["Original"] = ["NLP1","NLP2","BNLP"]
    UBD_f_relax["Interval"] = ["Interval"]

# Defines supported solvers for given problem relaxation
const LBD_p_relax = Dict{Any,Any}()
    LBD_p_relax["LP"] = ["CPLEX","Gurobi","Clp"]
    LBD_p_relax["QP"] = ["CPLEX","Gurobi"]
    LBD_p_relax["QCQP"] = ["CPLEX","Gurobi"]
    LBD_p_relax["NLP1"] = ["SNOPT"]
    LBD_p_relax["NLP2"] = ["Ipopt","SNOPT"]
    LBD_p_relax["BNLP"] = []
    LBD_p_relax["Interval"] = ["Interval"]
    LBD_p_relax["AlphaBB"] = ["AlphaBB"]
const UBD_p_relax = Dict{Any,Any}()
    UBD_p_relax["LP"] = ["CPLEX","Gurobi","Clp"]
    UBD_p_relax["QP"] = ["CPLEX","Gurobi"]
    UBD_p_relax["QCQP"] = ["CPLEX","Gurobi"]
    UBD_p_relax["NLP1"] = ["SNOPT"]
    UBD_p_relax["NLP2"] = ["Ipopt","SNOPT"]
    UBD_p_relax["BNLP"] = []
    UBD_p_relax["Interval"] = ["Interval"]

# checks to see if options are valid for problem, function, & solver settings
"""
    Solver_Relax_Valid_LBD!(opt)

Checks that relaxations used in lower bounding are appropriate for the solvers used.
"""
function Solver_Relax_Valid_LBD!(opt)

    # checks to see if function relaxation is supported
    if (~in(opt.LBD_func_relax,relax_list))
        error("Unsupported function relaxation.")
    end

    # checks that problem relaxation is valid for function relaxation
    if (~in(opt.LBD_problem_relax,LBD_f_relax[opt.LBD_func_relax]))
        error("Invalid combintation of function relaxation and problem relaxation.
                 Nondifferentiable McCormick relaxations support LP and bundle NLP relaxations (BNLP).
                 Differentiable relaxations support LP and NLP that require only first derivative info (BNLP, NLP1).
                 Twice differentiable relaxations support LP, QP, QPQC, and all NLP relaxations.")
    end

    # checks that solver is valid for problem relaxation
    if (~in(opt.LBD_problem_solver,LBD_p_relax[opt.LBD_problem_relax]))
        error("This solver does not support this type of problem relaxation.")
    end
end

"""
    Solver_Relax_Valid_UBD!(opt)

Checks that relaxations used in upper bounding are appropriate for the solvers used.
"""
function Solver_Relax_Valid_UBD!(opt)

    # checks to see if function relaxation is supported
    if (~in(opt.UBD_func_relax,relax_list_UBD))
        error("Unsupported function relaxation.")
    end

    # checks that problem relaxation is valid for function relaxation
    if (~in(opt.UBD_problem_relax,UBD_f_relax[opt.UBD_func_relax]))
        println("opt.UBD_problem_relax: ", opt.UBD_problem_relax)
        println("opt.UBD_func_relax: ", opt.UBD_func_relax)
        println("UBD_f_relax: ", UBD_f_relax)
        println("UBD_f_relax[opt.UBD_func_relax]: ", UBD_f_relax[opt.UBD_func_relax])
        error("Invalid combintation of function relaxation and problem relaxation.
         Nondifferentiable McCormick relaxations support LP and bundle NLP relaxations (BNLP).
         Differentiable relaxations support LP and NLP that require only first derivative info (BNLP, NLP1).
         Twice differentiable relaxations support LP, QP, QPQC, and all NLP relaxations.")
    end

    # checks that solver is valid for problem relaxation
    if (~in(opt.UBD_problem_solver,UBD_p_relax[opt.UBD_problem_relax]))
        println("opt.UBD_problem_relax: ", opt.UBD_problem_relax)
        println("opt.UBD_func_relax: ", opt.UBD_func_relax)
        println("UBD_p_relax: ", UBD_p_relax)
        println("UBD_f_relax[opt.UBD_func_relax]: ", UBD_p_relax[opt.UBD_func_relax])
        error("This solver does not support this type of problem relaxation.")
    end
end

verify_support(c) = c

"""
    verify_support(c::Expr)

Verifies that expression contains only supported operators.
"""
function verify_support(c::Expr)
    if c.head == :comparison
        map(verify_support, c.args)
        return c
    end
    if c.head == :call
        if c.args[1] in (:+, :-, :*, :/, :exp, :log, :sin, :cos, :tan, :sinh,
                         :cosh, :tanh, :exp2, :log2, :exp10, :log10, :min, :max,
                         :abs, :sqrt, :sqr, :sign, :asin, :acos, :atan, :asinh,
                         :acosh, :atanh)
            return c
        elseif c.args[1] in (:<=, :>=, :(==))
            map(verify_support, c.args[2:end])
            return c
        elseif c.args[1] == :^
            @assert isa(c.args[2], Real) || isa(c.args[3], Real)
            return c
        else
            error("Unsupported expression $c")
        end
    end
    return c
end
