# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eeago_optimizer/optimize/optimize_lp_cone.jl
# Contains the optimize! routines for LP, SOCP, (and in the future MILP and
# MISOCP) type problems. This also includes functions to add variables,
# linear constraints, soc constraints, and unpack solutions.
#############################################################################

#=
LP          -> COPY TO RELAXED SOLVER AND SOLVE
MILP        -> COPY TO RELAXED SOLVER AND SOLVE
SOCP        -> COPY TO RELAXED SOLVER AND SOLVE
MISOCP      -> COPY TO RELAXED SOLVER AND SOLVE
DIFF_CVX    -> COPY TO NLP SOLVER AND SOLVE (POTENTIAL MULTISTART)
NS_CVX      -> COPY TO NLP SOLVER AND SOLVE (POTENTIAL MULTISTART)
DIFF_NCVX   -> APPLY GLOBAL SOLVER (UNLESS USER REQUEST LOCAL SOLVE THEN NLP)
NS_NCVX     -> APPLY GLOBAL SOLVER (UNLESS USER REQUEST LOCAL SOLVE THEN NLP)
MINCVX      -> APPLY GLOBAL SOLVER (LOCAL SOLVE OPTION FUTURE FEATURE)
=#

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
function _add_linear_constraints!(m::Optimizer, opt::T) where T
    ip = m._input_problem
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _linear_leq_constraints(ip))
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _linear_geq_constraints(ip))
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _linear_eq_constraints(ip))
    return nothing
end

function _add_quadratic_constraints!(m::Optimizer, opt::T) where T
    ip = m._input_problem
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _quadratic_leq_constraints(ip))
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _quadratic_geq_constraints(ip))
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _quadratic_eq_constraints(ip))
    return nothing
end

### LP and MILP routines
function _add_soc_constraints!(m::Optimizer, opt::T) where T
    if MOI.supports_constraint(opt, VECOFVAR, SECOND_ORDER_CONE)
        ip = m._input_problem
        foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _second_order_cone_constraints(ip))
    else
        # TODO: FIX
        for (func, set) in m._input_problem._conic_second_order
            quad_terms = SQT[SQT((), func.variables[i], func.variables[i]) for i = 1:length(func.variables)]
            sqf = SQF(SQT[], SAF[], 0.0)
            MOI.add_constraint(opt, sqf, LT_ZERO)
        end
    end
    return nothing
end

function add_sv_or_aff_obj!(m::Optimizer, opt::T) where T
    if m._input_problem._objective_type == SINGLE_VARIABLE
        MOI.set(opt, MOI.ObjectiveFunction{SV}(), m._input_problem._objective_sv)
    elseif m._input_problem._objective_type == SCALAR_AFFINE
        MOI.set(opt, MOI.ObjectiveFunction{SAF}(), m._input_problem._objective_saf)
    end
    return nothing
end

function unpack_local_solve!(m::Optimizer, opt::T) where T

    m._maximum_node_id = 0

    m._termination_status_code = MOI.get(opt, MOI.TerminationStatus())
    m._result_status_code = MOI.get(opt, MOI.PrimalStatus())
    m._feasible_solution_found = m._result_status_code === MOI.FEASIBLE_POINT

    if MOI.get(opt, MOI.ResultCount()) > 0

        objective_value = MOI.get(opt, MOI.ObjectiveValue())

        # corrects for standard printing multiplier
        if m._input_problem._optimization_sense === MOI.MAX_SENSE
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

    m._relaxed_variable_index = add_variables(m, lp_optimizer, m._input_problem._variable_count)
    _add_linear_constraints!(m, lp_optimizer)
    _add_sv_or_aff_obj!(m, lp_optimizer)
    MOI.set(lp_optimizer, MOI.ObjectiveSense(), m._input_problem._optimization_sense)

    if m._parameters.verbosity < 5
        MOI.set(lp_optimizer, MOI.Silent(), true)
    end
    m._parse_time = time() - m._start_time

    MOI.optimize!(lp_optimizer)

    unpack_local_solve!(m, lp_optimizer)

    return nothing
end

optimize!(::Val{MILP}, m::Optimizer) = optimize!(Val{LP}(), m)

function optimize!(::Val{SOCP}, m::Optimizer)

    socp_optimizer = m.second_order_cone_optimizer
    set_default_config!(socp_optimizer)
    MOI.empty!(socp_optimizer)

    m._relaxed_variable_index = add_variables(m, socp_optimizer, m._input_problem._variable_count)
    _add_linear_constraints!(m, socp_optimizer)
    _add_soc_constraints!(m, socp_optimizer)
    _add_sv_or_aff_obj!(m, socp_optimizer)
    MOI.set(socp_optimizer, MOI.ObjectiveSense(), m._input_problem._optimization_sense)

    if m._parameters.verbosity < 5
        MOI.set(socp_optimizer, MOI.Silent(), true)
    end
    m._parse_time = time() - m._start_time

    MOI.optimize!(socp_optimizer)

    unpack_local_solve!(m, socp_optimizer)

    return nothing
end

optimize!(::Val{MISOCP}, m::Optimizer) = optimize!(Val{SOCP}(), m)
