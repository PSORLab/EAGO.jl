# List of valid function relaxations
relax_list = ["NS-MV-ON", "NS-MV-OFF", "NS-STD-ON", "NS-STD-OFF",
              "Diff1-MV-ON", "Diff1-MV-OFF", "Diff2-MV-ON", "Diff2-MV-OFF",
              "Interval","AlphaBB"]

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
end


const SupportedExpr = [:+, :-, :*, :/, :exp, :log, :sin, :cos, :tan, :sinh,
                 :cosh, :tanh, :exp2, :log2, :exp10, :log10, :min, :max,
                 :abs, :sqrt, :sqr, :sign, :asin, :acos, :atan, :asinh,
                 :acosh, :atanh]

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
        if c.args[1] in SupportedExpr
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

function registerEAGO(m::JuMP.Model, s::Symbol, dimension::Integer, f::Function; ad::Bool=false)
    push!(SupportedExpr,s)
    JuMP.register(m, s, dimension, f, autodiff=ad)
end
