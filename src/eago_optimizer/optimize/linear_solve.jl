#=
function add_variables(m::Optimizer, optimizer::T, variable_number::Int) where T

    variable_index = fill(VI(1), variable_number)
    for i = 1:variable_number
        @inbounds variable_index[i] = MOI.add_variable(optimizer)
        relaxed_variable = SV(@inbounds variable_index[i])
        v_info = @inbounds m._working_problem._variable_info[i]
        if v_info.is_integer && v_info.is_fixed
            MOI.add_constraint(optimizer, relaxed_variable, ET(v_info.lower_bound))
        elseif v_info.is_integer
            MOI.add_constraint(optimizer, relaxed_variable, ZO())
        elseif v_info.is_fixed
            MOI.add_constraint(optimizer, relaxed_variable, ET(v_info.lower_bound))
        elseif v_info.has_lower_bound && v_info.has_upper_bound
            MOI.add_constraint(optimizer, relaxed_variable, IT(v_info.lower_bound, v_info.upper_bound))
        elseif v_info.has_lower_bound
            MOI.add_constraint(optimizer, relaxed_variable, GT(v_info.lower_bound))
        elseif v_info.has_upper_bound
            MOI.add_constraint(optimizer, relaxed_variable, LT(v_info.upper_bound))
        end
    end

    return variable_index
end

### LP and MILP routines
function _add_linear_constraints!(ip::InputProblem, opt::T) where T
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _linear_leq_constraints(ip))
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _linear_geq_constraints(ip))
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _linear_eq_constraints(ip))
    return nothing
end


function add_sv_or_aff_obj!(m::Optimizer, opt::T) where T
    if _objective_type(_input_problem(m)) == SINGLE_VARIABLE
        MOI.set(opt, MOI.ObjectiveFunction{SV}(), m._input_problem._objective_sv)
    elseif _objective_type(_input_problem(m)) == SCALAR_AFFINE
        MOI.set(opt, MOI.ObjectiveFunction{SAF}(), m._input_problem._objective_saf)
    end
    return nothing
end

function _set_sense!(ip::InputProblem, opt::T) where T
    MOI.set(opt, MOI.ObjectiveSense(), _optimization_sense(ip))
end

function _unpack_local_solve!(m::Optimizer, opt::T) where T

    m._maximum_node_id = 0
    m._termination_status_code = MOI.get(opt, MOI.TerminationStatus())
    m._result_status_code = MOI.get(opt, MOI.PrimalStatus())
    m._feasible_solution_found = m._result_status_code == MOI.FEASIBLE_POINT

    if MOI.get(opt, MOI.ResultCount()) > 0

        objective_value = MOI.get(opt, MOI.ObjectiveValue())

        # corrects for standard printing multiplier
        if _optimization_sense(m._input_problem) == MOI.MAX_SENSE
            objective_value *= -1.0
        end

        m._global_lower_bound = objective_value
        m._global_upper_bound = objective_value
        m._objective_value = objective_value
        m._best_upper_value = objective_value
        m._solution_value = objective_value
    end

    map!(x -> MOI.get(opt, MOI.VariablePrimal(), x), m._continuous_solution, m._relaxed_variable_index)

    m._run_time = time() - m._start_time

    return nothing
end

function optimize!(::Val{LP}, m::Optimizer)

    lp_optimizer = m.lp_optimizer
    set_default_config!(lp_optimizer)
    MOI.empty!(lp_optimizer)

    ip = m._input_problem
    m._relaxed_variable_index = add_variables(m, lp_optimizer, _variable_count(ip))
    _add_linear_constraints!(ip, lp_optimizer)
    _add_sv_or_aff_obj!(ip, lp_optimizer)
    _set_sense!(ip, lp_optimizer)

    if m._parameters.verbosity < 5
        MOI.set(lp_optimizer, MOI.Silent(), true)
    end
    m._parse_time = time() - m._start_time

    MOI.optimize!(lp_optimizer)
    _unpack_local_solve!(m, lp_optimizer)

    return nothing
end

optimize!(::Val{MILP}, m::Optimizer) = optimize!(Val{LP}(), m)


function _add_quadratic_constraints!(ip::InputProblem, opt::T) where T
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _quadratic_leq_constraints(ip))
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _quadratic_geq_constraints(ip))
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _quadratic_eq_constraints(ip))
    return nothing
end
=#
