# OBBT
# Integer Variables: Fix one value and solve relaxation, fix other value and solve relaxtion, can dual simplex warm start potentially
# Apply to all nonlinear nonbinary variables
# Strong branching MINLP...

"""
    Excludes OBBT on variable indices that are tight for the solution of the relaxation (Should be done (for continuous)).
"""
function trival_filtering!(x::Optimizer, y::NodeBB)

    termination_status = MOI.get(x.working_relaxed_optimizer, MOI.TerminationStatus())
    result_status_code = MOI.get(x.working_relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(termination_status, result_status_code)

    if valid_flag
        if feasible_flag
            low_delete = Int[]
            upp_delete = Int[]

            for i in 1:length(x.obbt_working_lower_index)
                var_ref = x.obbt_working_lower_index[i]
                var_indx = var_ref.value
                variable_primal = MOI.get(x.working_relaxed_optimizer, MOI.VariablePrimal(), var_ref)

                val1 = abs(variable_primal - y.lower_variable_bounds[var_indx])
                val2 = x.obbt_tolerance

                (abs(variable_primal - y.lower_variable_bounds[var_indx]) <= x.obbt_tolerance) && push!(low_delete, i)
            end

            for i in 1:length(x.obbt_working_upper_index)
                var_ref = x.obbt_working_upper_index[i]
                var_indx = var_ref.value
                variable_primal = MOI.get(x.working_relaxed_optimizer, MOI.VariablePrimal(), var_ref)

                val1 = abs(y.upper_variable_bounds[var_indx] - variable_primal)
                val2 = x.obbt_tolerance

                (abs(y.upper_variable_bounds[var_indx] - variable_primal) <= x.obbt_tolerance) && push!(upp_delete, i)
            end
            deleteat!(x.obbt_working_lower_index, low_delete)
            deleteat!(x.obbt_working_upper_index, upp_delete)
        end
    end
end

#(Should be done (for continuous))
function aggressive_filtering!(x::Optimizer, y::NodeBB)

    # Initial filtering vector (negative one direction per remark in Gleixner2017)
    variable_number = x.variable_number
    v = -ones(Float64, variable_number)

    # Copy prior index set (ignores linear and binary terms)
    old_low_index = x.obbt_working_lower_index
    old_upp_index = x.obbt_working_upper_index
    new_low_index = x.obbt_working_lower_index
    new_upp_index = x.obbt_working_upper_index

    # Exclude unbounded directions
    low_delete = Int[]
    upp_delete = Int[]
    for i in new_low_index
        (y.lower_variable_bounds[i.value] != -Inf) && push!(low_delete, i.value)
    end
    for i in new_upp_index
        (y.upper_variable_bounds[i.value] != Inf) && push!(upp_delete, i.value)
    end
    deleteat!(new_low_index, low_delete); low_delete = Int[]
    deleteat!(new_upp_index, upp_delete); upp_delete = Int[]

    # Begin the main algorithm
    for k in 1:x.obbt_aggressive_max_iteration

        # Set index differences and vector for filtering direction
        lower_indx_diff = setdiff(old_low_index, new_low_index)
        upper_indx_diff = setdiff(old_upp_index, new_upp_index)
        low_delete = Int[]
        upp_delete = Int[]
        for i in lower_indx_diff
            (v[i.value] < 0.0) && (v[i.value] = 0.0)
        end
        for i in upper_indx_diff
            (v[i.value] > 0.0) && (v[i.value] = 0.0)
        end

        # Termination Condition
        (isempty(union(new_low_index,new_upp_index)) || (v == zeros(Float64,variable_number))) && break
        ((k >= 2) && (length(lower_indx_diff) + length(upper_indx_diff) < x.obbt_aggressive_min_dimension)) && break

        # Set objective in OBBT problem to filtering vector

        MOI.set(x.working_relaxed_optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(v,x.obbt_variables),0.0)
        MOI.set(x.working_relaxed_optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), saf)

        # Optimizes the problem and if successful filter additional bounds
        MOI.optimize!(x.working_relaxed_optimizer)

        termination_status = MOI.get(x.working_relaxed_optimizer, MOI.TerminationStatus())
        result_status_code = MOI.get(x.working_relaxed_optimizer, MOI.PrimalStatus())
        valid_flag, feasible_flag = is_globally_optimal(termination_status, result_status_code)

        if valid_flag
            if feasible_flag
                variable_primal = MOI.get(x.working_relaxed_optimizer, MOI.VariablePrimal(), x.obbt_variables)
                new_low_index = copy(old_low_index)
                new_upp_index = copy(old_upp_index)
                for i in old_low_index
                    (variable_primal[i.value] == y.lower_variable_bounds[i.value]) && push!(low_delete, i.value)
                end
                for i in old_upp_index
                    (variable_primal[i.value] == y.upper_variable_bounds[i.value]) && push!(upp_delete, i.value)
                end
                deleteat!(new_low_index, low_delete)
                deleteat!(new_upp_index, upp_delete)
            end
        else
            return false
        end
    end
    x.obbt_working_lower_index = new_low_index
    x.obbt_working_upper_index = new_upp_index
    return true
end

function calculate_initial_relaxation!(x::Optimizer, y::NodeBB)


    # build optimizer
    if x.use_lower_factory
        factory = x.lower_factory(;x.lower_optimizer_options...)
        x.working_relaxed_optimizer = factory
        MOI.add_variables(x.working_relaxed_optimizer, x.variable_number)
    else
        if x.initial_relaxed_optimizer != DummyOptimizer()
            x.working_relaxed_optimizer = deepcopy(x.initial_relaxed_optimizer)
        end
    end

    # relax functions
    xpnt = 0.5*(lower_variable_bounds(y) + upper_variable_bounds(y))
    update_lower_variable_bounds1!(x, y, x.working_relaxed_optimizer)
    x.relax_function!(x, x.working_relaxed_optimizer, y, x.relaxation, xpnt, load = true)
    x.relax_function!(x, x.working_relaxed_optimizer, y, x.relaxation, xpnt, load = false)

    MOI.set(x.working_relaxed_optimizer, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    MOI.optimize!(x.working_relaxed_optimizer)
end

function calculate_single_variable_relaxation!(opt::Optimizer, y::NodeBB, vari::MOI.VariableIndex, sense)

    # build optimizer
    if opt.use_lower_factory
        factory = opt.lower_factory(;opt.lower_optimizer_options...)
        opt.working_relaxed_optimizer = factory
        MOI.add_variables(opt.working_relaxed_optimizer, opt.variable_number)
    else
        if opt.initial_relaxed_optimizer != DummyOptimizer()
            opt.working_relaxed_optimizer = deepcopy(opt.initial_relaxed_optimizer)
        end
    end

    # relax functions
    xmid = 0.5*(lower_variable_bounds(y) + upper_variable_bounds(y))
    update_lower_variable_bounds1!(opt, y, opt.working_relaxed_optimizer)
    opt.relax_function!(opt, opt.working_relaxed_optimizer, y, opt.relaxation, xmid, load = true)
    opt.relax_function!(opt, opt.working_relaxed_optimizer, y, opt.relaxation, xmid, load = false)

    # set objectives
    var = MOI.SingleVariable(vari)
    MOI.set(opt.working_relaxed_optimizer, MOI.ObjectiveSense(), sense)
    MOI.set(opt.working_relaxed_optimizer, MOI.ObjectiveFunction{MOI.SingleVariable}(), var)

    modelsense = opt.working_relaxed_optimizer.obj_sense
    modeltype = opt.working_relaxed_optimizer.obj_type
    modelsingobjvar = opt.working_relaxed_optimizer.single_obj_var
    modelobj_const = opt.working_relaxed_optimizer.objective_constant

    MOI.optimize!(opt.working_relaxed_optimizer)
end

function aggressive_obbt_on_heurestic(x::Optimizer)
    x.obbt_aggressive_on
end

"""
    obbt

Performs OBBT with Filtering and Greedy Ordering
"""
function obbt(x::Optimizer, y::NodeBB)

    feasibility = true

    calculate_initial_relaxation!(x, y)

    # Sets indices to attempt OBBT on
    append!(x.obbt_working_lower_index, x.obbt_variables)
    append!(x.obbt_working_upper_index, x.obbt_variables)

    # Prefiltering steps && and sets initial LP values
    trival_filtering!(x, y)

    if aggressive_obbt_on_heurestic(x)
        feasibility = aggressive_filtering!(x, y)
    end
    xLP = MOI.get(x.working_relaxed_optimizer, MOI.VariablePrimal(), [MOI.VariableIndex(i) for i in 1:x.variable_number])

    while ~(isempty(x.obbt_working_lower_index) && isempty(x.obbt_working_upper_index)) && ~isempty(y)

        # Get lower value
        lower_indx = 0
        upper_indx = 0
        if (isempty(x.obbt_working_lower_index))
            lower_value = Inf
        else
            lower_value = xLP[x.obbt_working_lower_index[1].value] - y.lower_variable_bounds[x.obbt_working_lower_index[1].value]
            lower_indx = 1
            for i in 2:length(x.obbt_working_lower_index)
                temp_value = xLP[x.obbt_working_lower_index[i].value] - y.lower_variable_bounds[x.obbt_working_lower_index[i].value]
                if temp_value < lower_value
                    lower_value = temp_value
                    lower_indx = i
                end
            end
        end
        if (isempty(x.obbt_working_upper_index))
            upper_value = Inf
        else
            upper_value = y.upper_variable_bounds[x.obbt_working_upper_index[1].value] - xLP[x.obbt_working_upper_index[1].value]
            upper_indx = 1
            for i in 2:length(x.obbt_working_upper_index)
                temp_value = y.upper_variable_bounds[x.obbt_working_upper_index[i].value] - xLP[x.obbt_working_upper_index[i].value]
                if temp_value < upper_value
                    upper_value = temp_value
                    upper_indx = 1
                end
            end
        end

        if (lower_value <= upper_value)

            var_ref = x.obbt_working_lower_index[lower_indx]
            deleteat!(x.obbt_working_lower_index, lower_indx)

            calculate_single_variable_relaxation!(x, y, var_ref, MOI.MIN_SENSE)
            termination_status = MOI.get(x.working_relaxed_optimizer, MOI.TerminationStatus())
            result_status_code = MOI.get(x.working_relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(termination_status, result_status_code)

            if valid_flag
                if feasible_flag
                    xLP[:] = MOI.get(x.working_relaxed_optimizer, MOI.VariablePrimal(), [MOI.VariableIndex(i) for i in 1:x.variable_number])
                    y.lower_variable_bounds[var_ref.value] = is_integer_variable(x,var_ref.value) ? ceil(xLP[var_ref.value]) : xLP[var_ref.value]
                    if isempty(y)
                        feasibility = false
                        break
                    end
                else
                    feasibility = false
                    x.obbt_working_lower_index = []
                    x.obbt_working_upper_index = []
                end
            else
                break
            end

        else

            # Get index
            var_ref = x.obbt_working_upper_index[upper_indx]
            deleteat!(x.obbt_working_upper_index, upper_indx)

            # Solve optimization model
            calculate_single_variable_relaxation!(x, y, var_ref, MOI.MAX_SENSE)
            termination_status = MOI.get(x.working_relaxed_optimizer, MOI.TerminationStatus())
            result_status_code = MOI.get(x.working_relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(termination_status, result_status_code)

            if valid_flag
                if feasible_flag
                    xLP[:] = MOI.get(x.working_relaxed_optimizer, MOI.VariablePrimal(), [MOI.VariableIndex(i) for i in 1:x.variable_number])
                    y.upper_variable_bounds[var_ref.value] = is_integer_variable(x,var_ref.value) ? ceil(xLP[var_ref.value]) : xLP[var_ref.value]
                    if isempty(y)
                        feasibility = false
                        break
                    end
                else
                    feasibility = false
                    x.obbt_working_lower_index = []
                    x.obbt_working_upper_index = []
                end
            else
                break
            end
        end
        #GenerateLVB!(x)
        trival_filtering!(x,y)
    end

    return feasibility
end
