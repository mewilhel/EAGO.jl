function reset_relaxation!(m::GlobalOptimizer{N,T})

    m._working_problem._relaxed_evaluator.is_first_eval = true
    fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)

    m._new_eval_constraint = true

    delete_nl_constraints!(m)
    delete_objective_cuts!(m)

    return nothing
end

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
function set_dual!(m::GlobalOptimizer{N,T})
    opt = m.relaxed_optimizer
    foreach(t -> (m._lower_uvd[t[2]] = MOI.get(opt, MOI.ConstraintDual(), t[1]);), m._relaxed_variable_lt)
    foreach(t -> (m._lower_lvd[t[2]] = MOI.get(opt, MOI.ConstraintDual(), t[1]);), m._relaxed_variable_gt)
    return nothing
end

"""
$(SIGNATURES)

Updates the relaxed constraint by setting the constraint set of `v == x*`` ,
`xL_i <= x_i`, and `x_i <= xU_i` for each such constraint added to the relaxed
optimizer.
"""
function update_relaxed_problem_box!(m::GlobalOptimizer{N,T})
    opt = m.relaxed_optimizer
    lb = _lower_variable_bound(m)
    ub = _upper_variable_bound(m)
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], ET(lb[x[2]])), m._relaxed_variable_eq)
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], LT(ub[x[2]])), m._relaxed_variable_lt)
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], GT(ub[x[2]])), m._relaxed_variable_gt)
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], IT(lb[x[2]], ub[x[2]])), m._relaxed_variable_it)
    # Interval... & ZO interaction need to fix
    return nothing
end

function interval_objective_bound(m::GlobalOptimizer{N,T}, n::NodeBB)

    interval_objective_bound = bound_objective(m)

    if interval_objective_bound > m._lower_objective_value
        m._lower_objective_value = interval_objective_bound
        fill!(m._lower_lvd, 0.0)
        fill!(m._lower_uvd, 0.0)
        m._cut_add_flag = false
        return true
    end

    return false
end

_is_feas(m::GlobalOptimizer{N,T}, x::AffineFunctionIneq, n) = lower_interval_bound(m, x, n) <= 0.0
_is_feas(m::GlobalOptimizer{N,T}, x::BufferedQuadraticIneq, n) = lower_interval_bound(m, x, n) <= 0.0
function _is_feas(m::GlobalOptimizer{N,T}, x::AffineFunctionEq, n)
    lower_value, upper_value = interval_bound(m, x, n)
    return lower_value <= 0.0 <= upper_value
end
function _is_feas(m::GlobalOptimizer{N,T}, x::BufferedQuadraticEq, n)
    lower_value, upper_value = interval_bound(m, x, n)
    return lower_value <= 0.0 <= upper_value
end
function _is_feas(m::GlobalOptimizer{N,T}, x::NonlinearExpression, n)
    lower_value, upper_value = interval_bound(m, x, n)
    feasible_flag &= upper_value < x.lower_bound
    feasible_flag &= lower_value > x.upper_bound
    return !feasible_flag
end

"""
$(SIGNATURES)

A fallback lower bounding problem that consists of an natural interval extension
calculation. This is called when the optimizer used to compute the lower bound
does not return a termination and primal status code indicating that it
successfully solved the relaxation to a globally optimal point.
"""
function fallback_interval_lower_bound!(m::GlobalOptimizer{N,T}, n::NodeBB)

    feas = true
    wp = _working_problem(m)
    if !cp_condition(m)
        feas = foreach_until(x -> _is_feas(m, x, n), wp._saf_leq, feas)
        feas = foreach_until(x -> _is_feas(m, x, n), wp._saf_eq,  feas)
        feas = foreach_until(x -> _is_feas(m, x, n), wp._sqf_leq, feas)
        feas = foreach_until(x -> _is_feas(m, x, n), wp._sqf_eq,  feas)
        feas = foreach_until(x -> _is_feas(m, x, n), wp._nonlinear_constr, feas)
    end

    if feas
        interval_objective_used = interval_objective_bound(m, n)
        @__dot__ m._current_xref = 0.5*(n.upper_variable_bound + n.lower_variable_bound)
        unsafe_check_fill!(isnan, m._current_xref, 0.0, length(m._current_xref))
    else
        m._lower_objective_value = -Inf
    end
    m._lower_feasibility = feas

    return
end

function intrepret_relaxed_solution(m::GlobalOptimizer, d::T) where T
    valid_flag, feasible_flag = is_globally_optimal(m._lower_termination_status,
                                                    m._lower_result_status)
    if valid_flag && feasible_flag
        set_dual!(m)
        m._cut_add_flag = true
        m._lower_feasibility = true
        m._lower_objective_value = MOI.get(d, MOI.ObjectiveValue())
        for i = 1:_variable_num(_working_problem(m))
             m._lower_solution[i] = MOI.get(d, MOI.VariablePrimal(),
                                               m._relaxed_variable_index[i])
        end
    elseif valid_flag
        m._cut_add_flag = false
        m._lower_feasibility  = false
        m._lower_objective_value = -Inf
    end
    return valid_flag
end

"""
$(SIGNATURES)

Constructs and solves the relaxation using the default EAGO relaxation scheme
and optimizer on node `y`.
"""
function lower_problem!(t::ExtensionType, m::GlobalOptimizer{N,T})

    n = m._current_node

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

    return nothing
end

"""
$(SIGNATURES)

Updates the internal storage in the optimizer after a valid feasible cut is added.
"""
function cut_update!(m::GlobalOptimizer{N,T})

    m._cut_feasibility = true

    relaxed_optimizer = m.relaxed_optimizer
    obj_val = MOI.get(relaxed_optimizer, MOI.ObjectiveValue())
    prior_obj_val = (m._cut_iterations == 2) ? m._lower_objective_value : m._cut_objective_value

    m._cut_add_flag = true
    m._lower_termination_status = m._cut_termination_status
    m._lower_result_status = m._cut_result_status
    m._cut_solution[:] = MOI.get(relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index)

    if prior_obj_val < obj_val
        m._cut_objective_value = obj_val
        m._lower_objective_value = obj_val
        set_dual!(m)
        copyto!(m._lower_solution, m._cut_solution)

    else
        m._cut_objective_value = prior_obj_val
        m._lower_objective_value = prior_obj_val
        m._cut_add_flag = false
    end

    return nothing
end


"""
$(SIGNATURES)

Checks if a cut should be added and computes a new reference point to add the
cut at. If no cut should be added the constraints not modified in place are
deleted from the relaxed optimizer and the solution is compared with the
interval lower bound. The best lower bound is then used.
"""
function cut_condition(t::ExtensionType, m::GlobalOptimizer{N,T})

    # always add cut if below the minimum iteration limit, otherwise add cut
    # the number of cuts is less than the maximum and the distance between
    # prior solutions exceeded a tolerance.
    continue_cut_flag = m._cut_add_flag
    continue_cut_flag &= (m._cut_iterations < m.cut_max_iterations)

    # compute distance between prior solutions and compare to tolerances
    n = m._current_node
    ns_indx = m._branch_to_sol_map

    cvx_factor =  m.cut_cvx
    xsol = (m._cut_iterations > 1) ? m._cut_solution[ns_indx] : m._lower_solution[ns_indx]
    xnew = (1.0 - cvx_factor)*mid(n) + cvx_factor*xsol

    continue_cut_flag &= (norm((xsol - xnew)/diam(n), 1) > m.cut_tolerance)
    continue_cut_flag |= (m._cut_iterations < m.cut_min_iterations)

    # update reference point for new cut
    if continue_cut_flag
        copyto!(m._current_xref, xnew)
        if m._nonlinear_evaluator_created
            set_reference_point!(m)
            fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)
        end
    end

    # check to see if interval bound is preferable and replaces the objective
    # value with the interval value if so. Any available dual values are then
    # set to zero since the interval bounds are by definition constant
    if m._lower_feasibility && !continue_cut_flag
        objective_lo = lower_interval_bound(m, m._working_problem._objective_parsed, n)
        if objective_lo > m._lower_objective_value
            m._lower_objective_value = objective_lo
            fill!(m._lower_lvd, 0.0)
            fill!(m._lower_uvd, 0.0)
        end
    end

    m._cut_iterations += 1

    return continue_cut_flag
end

"""
$(SIGNATURES)

Adds a cut for each constraint and the objective function to the subproblem.
"""
function add_cut!(t::ExtensionType, m::GlobalOptimizer{N,T})

    fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)
    m._working_problem._relaxed_evaluator.is_first_eval = true
    m._working_problem._relaxed_evaluator.is_intersect = false
    m._new_eval_constraint = true

    relax_constraints!(m, m._cut_iterations)

    # Optimizes the object
    relaxed_optimizer = m.relaxed_optimizer
    MOI.optimize!(relaxed_optimizer)

    m._cut_termination_status = MOI.get(relaxed_optimizer, MOI.TerminationStatus())
    m._cut_result_status = MOI.get(relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(m._cut_termination_status, m._cut_result_status)

    if valid_flag && feasible_flag
        cut_update!(m)

    elseif valid_flag
        m._cut_add_flag = false
        m._lower_feasibility  = false
        m._lower_objective_value = -Inf

    else
        m._cut_add_flag = false
    end

    return nothing
end
