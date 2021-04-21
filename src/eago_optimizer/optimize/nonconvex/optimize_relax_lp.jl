"""
$(SIGNATURES)

Takes an `MOI.TerminationStatusCode` and a `MOI.ResultStatusCode` and returns
the tuple `(valid_result::Bool, feasible::Bool)`. The value `valid_result` is
`true` if the pair of codes prove that either the subproblem solution was solved
to global optimality or the subproblem solution is infeasible. The value of
`feasible` is true if the problem is feasible and false if the problem is infeasible.
"""
function is_globally_optimal(t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)

    if t == MOI.INFEASIBLE && r == MOI.INFEASIBILITY_CERTIFICATE
        valid_result = true
        feasible = false

    elseif t == MOI.INFEASIBLE && r == MOI.NO_SOLUTION
        valid_result = true
        feasible = false

    elseif t == MOI.INFEASIBLE && r == MOI.UNKNOWN_RESULT_STATUS
        valid_result = true
        feasible = false

    elseif t == MOI.OPTIMAL && r == MOI.FEASIBLE_POINT
        valid_result = true
        feasible = true

    elseif t == MOI.INFEASIBLE_OR_UNBOUNDED && r == MOI.NO_SOLUTION
        valid_result = true
        feasible = false

    else
        valid_result = false
        feasible = false
    end

    return valid_result, feasible
end

"""
$(SIGNATURES)

Retrieves the lower and upper duals for variable bounds from the
`relaxed_optimizer` and sets the appropriate values in the
`_lower_lvd` and `_lower_uvd` storage fields.
"""
function _set_dual!(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
    opt = m.relaxed_optimizer
    foreach(t -> (m._lower_uvd[t[2]] = MOI.get(opt, MOI.ConstraintDual(), t[1]);), m._relaxed_variable_lt)
    foreach(t -> (m._lower_lvd[t[2]] = MOI.get(opt, MOI.ConstraintDual(), t[1]);), m._relaxed_variable_gt)
    return nothing
end

function interval_objective_bound(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
    interval_objective_bound = bound_objective(m)
    if interval_objective_bound > m._lower_objective_value
        m._lower_objective_value = interval_objective_bound
        fill!(m._lower_lvd, zero(T))
        fill!(m._lower_uvd, zero(T))
        return true
    end
    return false
end

"""
$(SIGNATURES)

A fallback lower bounding problem that consists of an natural interval extension
calculation. This is called when the optimizer used to compute the lower bound
does not return a termination and primal status code indicating that it
successfully solved the relaxation to a globally optimal point.
"""
function fallback_interval_lower_bound!(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}

    feas = true
    wp = _working_problem(m)
    if !cp_condition(m)
        feas = foreach_until(x -> _is_feas(Val{:leq}, m, x), wp._saf_leq, feas)
        feas = foreach_until(x -> _is_feas(Val{:eq}, m, x), wp._saf_eq, feas)
        feas = foreach_until(x -> _is_feas(m, x), wp._sqf_leq, feas)
        feas = foreach_until(x -> _is_feas(m, x), wp._sqf_eq,  feas)
        feas = foreach_until(x -> _is_feas(m, x), wp._nonlinear_constr, feas)
    end

    if feas
        interval_objective_used = interval_objective_bound(m)
        @__dot__ m._current_xref = _mid(BranchVar, m, 1:N)
        unsafe_check_fill!(isnan, m._current_xref, zero(T), N)
    else
        m._lower_objective_value = -Inf
    end
    m._lower_feasibility = feas

    return
end

function _update_box!(m::GlobalOptimizer{N,T,S}, v::Val{:cont}) where {N,T<:AbstractFloat,S}
    wp = _working_problem(m)
    lb = _lower_bound(FullVar, m)
    ub = _upper_bound(FullVar, m)
    opt = m.mip_optimizer
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], ET(lb[x[2]])),           m._relaxed_variable_eq)
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], LT(ub[x[2]])),           m._relaxed_variable_lt)
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], GT(ub[x[2]])),           m._relaxed_variable_gt)
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], IT(lb[x[2]], ub[x[2]])), m._relaxed_variable_it)
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], IT(zero(T), one(T))),    m._relaxed_variable_zo)
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], IT(lb[x[2]], ub[x[2]])), m._relaxed_variable_int)
    return nothing
end

function _update_box!(m::GlobalOptimizer{N,T,S}, v::Val{:mixed}) where {N,T<:AbstractFloat,S}
    wp = _working_problem(m)
    for x in m._relaxed_variable_zo
        l = _lower_bound(FullVar, m, x[2])
        u = _upper_bound(FullVar, m, x[2])
        if u < one(T)
            MOI.set(wp.mip_optimizer, MOI.ConstraintSet(), x[1], ET(zero(T)))
        elseif l > zero(T)
            MOI.set(wp.mip_optimizer, MOI.ConstraintSet(), x[1], ET(one(T)))
        else
            MOI.set(wp.mip_optimizer, MOI.ConstraintSet(), x[1], ZO())
        end
    end
    for x in m._relaxed_variable_int
        l = _lower_bound(FullVar, m, x[2])
        u = _upper_bound(FullVar, m, x[2])
        if iszero(l - floor(u))
            MOI.set(wp.mip_optimizer, MOI.ConstraintSet(), x[1], ET(l))
        elseif iszero(u - ceil(l))
            MOI.set(wp.mip_optimizer, MOI.ConstraintSet(), x[1], ET(u))
        else
            MOI.set(wp.mip_optimizer, MOI.ConstraintSet(), x[1], MOI.Semiinteger(l,u))
        end
    end
    return nothing
end

function _reset_objective!(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
    wp = _working_problem(m)
    MOI.set(wp.mip_optimizer, MOI.ObjectiveFunction{SAF}(), wp._objective)
    return nothing
end

"""
$(SIGNATURES)

Constructs and solves the relaxation using the default EAGO relaxation scheme
and optimizer on node `y`.
"""
function lower_problem!(t::ExtensionType, m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}

    if !_relaxation_available(m)
        _update_box!(m, Val(:cont))
        _relax_problem!(m, Val(:cont))
    else
        _reset_objective!(m)
    end

    # solve LPs, adding linear cuts
    MOI.set(m.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(m.relaxed_optimizer)
    valid_flag, feasible_flag = is_globally_optimal(m._lower_termination_status,
                                                    m._lower_result_status)

    # add continuous linear cuts
    linear_cut_condition = true
    while valid_flag && feasible_flag && linear_cut_condition
        _relax_problem!(m, Val(:cont))
        MOI.optimize!(m.relaxed_optimizer)
        valid_flag, feasible_flag = is_globally_optimal(m._lower_termination_status,
                                                        m._lower_result_status)
        linear_cut_condition = false # TODO: ADD CUT CHECK...
    end

    # activate integer variables and solve mixed integer formulation
    is_integer_prob = _integer_variable_num(m) > 0
    if valid_flag && feasible_flag && is_integer_prob
        _update_box!(m, Val(:mixed))      # need to update qp and conic relaxations as well
        _relax_problem!(m, Val(:mixed))
        MOI.optimize!(m.relaxed_optimizer)
    end

    if valid_flag && feasible_flag
        _set_dual!(m)
        m._lower_feasibility = true
        m._lower_objective_value = MOI.get(d, MOI.ObjectiveValue())
        result_count = MOI.get(d, MOI.ResultCount())
        for i = m._lower_result_count_max:(result_count-1)
            push!(m._lower_solution, zeros(m._working_problem._variable_num))
        end
        for j = 1:result_count
            for i = 1:_variable_num(_working_problem(m))
                 m._lower_solution[j][i] = MOI.get(d, MOI.VariablePrimal(j), m._relaxed_variable_index[i])
             end
         end
         m._lower_result_count = result_count
    elseif valid_flag && !feasible_flag
        m._lower_feasibility  = false
        m._lower_objective_value = -Inf
    else
        fallback_interval_lower_bound!(m)
    end

    return nothing
end

"""
$(SIGNATURES)

Checks if additional cuts should be added by the `add_cut!` routine.
"""
cut_condition(t::ExtensionType, m::GlobalOptimizer) = false

"""
$(SIGNATURES)

Add additionals cuts to subproblem while `cut_condition` is satisfied subproblem.
"""
add_cut!(t::ExtensionType, m::GlobalOptimizer) = nothing


#=
m._working_problem._relaxed_evaluator.is_post = m.subgrad_tighten
if !m._obbt_performed_flag
    if m._nonlinear_evaluator_created
        set_node!(m._working_problem._relaxed_evaluator, n)
        set_node_flag!(m)
        set_reference_point!(m)
        fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)
    end
    update_relaxed_problem_box!(m)
end
m._working_problem._relaxed_evaluator.interval_intersect = false

if !m._obbt_performed_flag
    relax_constraints!(m, 1)
end

# Optimizes the object
MOI.set(m.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
MOI.optimize!(m.relaxed_optimizer)
if !intrepret_relaxed_solution(m, m.relaxed_optimizer)
    fallback_interval_lower_bound!(m, n)
end
=#
