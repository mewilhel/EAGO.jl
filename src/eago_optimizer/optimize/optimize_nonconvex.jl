# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimize/optimize_nonconvex.jl
# Contains the optimize! routine and subroutines needed in the branch and
# bound routine called by EAGO.
#############################################################################

function _add_linear_constraints!(ip::InputProblem, opt::T) where T
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _linear_leq_constraint(ip))
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _linear_geq_constraint(ip))
    foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _linear_eq_constraint(ip))
    return nothing
end


function set_evaluator_flags!(d, is_post, is_intersect, is_first_eval, interval_intersect)

    d.is_post = is_post
    d.is_intersect = is_intersect
    d.is_first_eval = is_first_eval
    d.interval_intersect = interval_intersect

    return nothing
end

function reset_relaxation!(m::Optimizer)

    m._working_problem._relaxed_evaluator.is_first_eval = true
    fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)

    m._new_eval_objective = true
    m._new_eval_constraint = true

    delete_nl_constraints!(m)
    delete_objective_cuts!(m)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Creates an initial node with initial box constraints and adds it to the stack.
"""
function create_initial_node!(m::Optimizer)

    branch_variable_num = m._branch_variable_num
    lower_bound = zeros(Float64, branch_variable_num)
    upper_bound = zeros(Float64, branch_variable_num)
    branch_count = 1

    for i = 1:m._working_problem._variable_num
        vi = m._working_problem._variable_info[i]
        if vi.branch_on === BRANCH
            lower_bound[branch_count] = vi.lower_bound
            upper_bound[branch_count] = vi.upper_bound
            branch_count += 1
        end
    end

    n = NodeBB(lower_bound, upper_bound, -Inf, Inf, 1, 1)
    push!(m._stack, n)
    m._node_count = 1
    m._maximum_node_id += 1

    return nothing
end

"""
$(TYPEDSIGNATURES)

Loads variables, linear constraints, and empty storage for first nlp and
quadratic cut.
"""
function load_relaxed_problem!(m::Optimizer)

    relaxed_optimizer = m.relaxed_optimizer

    # add variables and indices and constraints
    wp = m._working_problem
    branch_variable_num = 0

    variable_count = wp._variable_num
    for i = 1:variable_count

        relaxed_variable_indx = MOI.add_variable(relaxed_optimizer)
        relaxed_variable = SV(relaxed_variable_indx)
        push!(m._relaxed_variable_index, relaxed_variable_indx)

        vinfo =  wp._variable_info[i]

        is_branch_variable =  m._branch_variables[i]
        vinfo.branch_on = is_branch_variable ? BRANCH : NO_BRANCH
        is_branch_variable && (branch_variable_num += 1)

        if vinfo.is_integer

        elseif vinfo.is_fixed
            ci_sv_et = MOI.add_constraint(relaxed_optimizer, relaxed_variable, ET(vinfo.lower_bound))
            if is_branch_variable
                push!(m._relaxed_variable_eq, (ci_sv_et, branch_variable_num))
                wp._var_eq_count += 1
            end

        else
            if vinfo.has_lower_bound
                ci_sv_gt = MOI.add_constraint(relaxed_optimizer, relaxed_variable, GT(vinfo.lower_bound))
                if is_branch_variable
                    push!(m._relaxed_variable_gt, (ci_sv_gt, branch_variable_num))
                    wp._var_geq_count += 1
                end
            end

            if vinfo.has_upper_bound
                ci_sv_lt = MOI.add_constraint(relaxed_optimizer, relaxed_variable, LT(vinfo.upper_bound))
                if is_branch_variable
                    push!(m._relaxed_variable_lt, (ci_sv_lt, branch_variable_num))
                    wp._var_leq_count += 1
                end
            end
        end
    end

    # set node index to single variable constraint index maps
    m._node_to_sv_leq_ci = fill(CI{SV,LT}(-1), branch_variable_num)
    m._node_to_sv_geq_ci = fill(CI{SV,GT}(-1), branch_variable_num)
    for i = 1:wp._var_leq_count
        ci_sv_lt, branch_index = m._relaxed_variable_lt[i]
        m._node_to_sv_leq_ci[branch_index] = ci_sv_lt
    end
    for i = 1:wp._var_geq_count
        ci_sv_gt, branch_index = m._relaxed_variable_gt[i]
        m._node_to_sv_geq_ci[branch_index] = ci_sv_gt
    end

    # set number of variables to branch on
    m._branch_variable_num = branch_variable_num

    # add linear constraints
    _add_linear_constraints!(m._input_problem, relaxed_optimizer)

    # sets relaxed problem objective sense to Min as all problems
    # are internally converted in Min problems in EAGO
    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    return nothing
end

function presolve_global!(t::ExtensionType, m::Optimizer)

    load_relaxed_problem!(m)
    create_initial_node!(m)

    branch_variable_num = m._branch_variable_num

    m._current_xref             = fill(0.0, branch_variable_num)
    m._candidate_xref           = fill(0.0, branch_variable_num)
    m._current_objective_xref   = fill(0.0, branch_variable_num)
    m._prior_objective_xref     = fill(0.0, branch_variable_num)
    m._lower_lvd                = fill(0.0, branch_variable_num)
    m._lower_uvd                = fill(0.0, branch_variable_num)

    # populate in full space until local MOI nlp solves support constraint deletion
    # uses input model for local nlp solves... may adjust this if a convincing reason
    # to use a reformulated upper problem presents itself
    m._lower_solution      = zeros(Float64, m._working_problem._variable_num)
    m._cut_solution        = zeros(Float64, m._working_problem._variable_num)
    m._solution = zeros(Float64, m._working_problem._variable_num)
    m._upper_solution      = zeros(Float64, m._working_problem._variable_num)
    m._upper_variables     = fill(VI(-1), m._working_problem._variable_num)

    # add storage for fbbt
    m._lower_fbbt_buffer   = zeros(Float64, m._working_problem._variable_num)
    m._upper_fbbt_buffer   = zeros(Float64, m._working_problem._variable_num)

    # add storage for obbt ( perform obbt on all relaxed variables, potentially)
    m._obbt_working_lower_index = fill(false, branch_variable_num)
    m._obbt_working_upper_index = fill(false, branch_variable_num)
    m._old_low_index            = fill(false, branch_variable_num)
    m._old_upp_index            = fill(false, branch_variable_num)
    m._new_low_index            = fill(false, branch_variable_num)
    m._new_upp_index            = fill(false, branch_variable_num)
    m._lower_indx_diff          = fill(false, branch_variable_num)
    m._upper_indx_diff          = fill(false, branch_variable_num)
    m._obbt_variable_num      = branch_variable_num

    # add storage for objective cut if quadratic or nonlinear
    wp = m._working_problem
    obj_type = wp._objective_type
    if obj_type === SCALAR_QUADRATIC
        wp._objective_saf.terms = copy(wp._objective_sqf.saf.terms)
    elseif obj_type === NONLINEAR
        wp._objective_saf.terms = copy(wp._objective_nl.saf.terms)
    end

    # set subgradient refinement flag
    wp._relaxed_evaluator.is_post = m.subgrad_tighten
    wp._relaxed_evaluator.subgrad_tighten = m.subgrad_tighten
    wp._relaxed_evaluator.reverse_subgrad_tighten =  m.reverse_subgrad_tighten

    m._presolve_time = time() - m._parse_time

    return nothing
end

"""
$(SIGNATURES)

Selects node with the lowest lower bound in stack.
"""
function node_selection!(t::ExtensionType, m::Optimizer)
    m._node_count -= 1
    m._current_node = popmin!(m._stack)
    return nothing
end

"""
$(SIGNATURES)

Creates two nodes from `current_node` using information available the `x`
and stores them to the stack. By default, relative width bisection is perfomed
at a point `branch_pnt` which is a convex combination
(parameter: `branch_cvx_factor`) of the solution to the relaxation and
the midpoint of the node. If this solution lies within `branch_offset/width` of
a bound then the branch point is moved to a distance of `branch_offset/width`
from the bound.
"""
function branch_node!(t::ExtensionType, m::Optimizer)

    n = m._current_node

    lvbs = n.lower_variable_bounds
    uvbs = n.upper_variable_bounds

    max_pos = 0
    max_val = -Inf
    temp_max = 0.0

    flag = true
    for i = 1:m._branch_variable_num
        si = m._branch_to_sol_map[i]
        vi = m._working_problem._variable_info[si]
        if vi.branch_on === BRANCH
            temp_max =  uvbs[i] - lvbs[i]
            temp_max /= vi.upper_bound - vi.lower_bound
            if temp_max > max_val
                max_pos = i
                max_val = temp_max
            end
        end
    end

    lvb  = lvbs[max_pos]
    uvb  = uvbs[max_pos]
    si   = m._branch_to_sol_map[max_pos]
    lsol = m._lower_solution[si]

    cvx_f = m.branch_cvx_factor
    cvx_g = m.branch_offset

    branch_pnt = cvx_f*lsol + (1.0 - cvx_f)*(lvb + uvb)/2.0
    if branch_pnt < lvb*(1.0 - cvx_g) + cvx_g*uvb
        branch_pnt = (1.0 - cvx_g)*lvb + cvx_g*uvb
    elseif branch_pnt > cvx_g*lvb + (1.0 - cvx_g)*uvb
        branch_pnt = cvx_g*lvb + (1.0 - cvx_g)*uvb
    end

    # rounds into branch points, which in turn prevents the
    # solution at the branch point from being discarded
    N1::Interval{Float64} = Interval{Float64}(lvb, branch_pnt)
    N2::Interval{Float64} = Interval{Float64}(branch_pnt, uvb)
    lvb_1 = copy(lvbs)
    uvb_1 = copy(uvbs)
    lvb_2 = copy(lvbs)
    uvb_2 = copy(uvbs)
    lvb_1[max_pos] = N1.lo
    uvb_1[max_pos] = N1.hi
    lvb_2[max_pos] = N2.lo
    uvb_2[max_pos] = N2.hi

    lower_bound = max(n.lower_bound, m._lower_objective_value)
    upper_bound = min(n.upper_bound, m._upper_objective_value)
    new_depth = n.depth + 1

    m._maximum_node_id += 1
    X1 = NodeBB(lvb_1, uvb_1, lower_bound, upper_bound, new_depth, m._maximum_node_id)
    m._maximum_node_id += 1
    X2 = NodeBB(lvb_2, uvb_2, lower_bound, upper_bound, new_depth, m._maximum_node_id)

    push!(m._stack, X1)
    push!(m._stack, X2)

    m._node_repetitions = 1
    m._node_count += 2

    return nothing
end

"""
$(SIGNATURES)

Stores the current node to the stack after updating lower/upper bounds.
"""
function single_storage!(t::ExtensionType, m::Optimizer)
    y = m._current_node
    m._node_repetitions += 1
    m._node_count += 1
    lower_bound = max(y.lower_bound, m._lower_objective_value)
    upper_bound = min(y.upper_bound, m._upper_objective_value)
    push!(m._stack, NodeBB(y.lower_variable_bounds, y.upper_variable_bounds,
                           lower_bound, upper_bound, y.depth, y.id))

    return nothing
end

"""
$(SIGNATURES)

Selects and deletes nodes from stack with lower bounds greater than global
upper bound.
"""
function fathom!(t::ExtensionType, m::Optimizer)

    upper = m._global_upper_bound
    continue_flag = !isempty(m._stack)

    while continue_flag
        max_node = maximum(m._stack)
        max_check = (max_node.lower_bound > upper)

        if max_check
            popmax!(m._stack)
            m._node_count -= 1
            if isempty(m._stack)
                continue_flag = false
            end

        else
            if !max_check
                continue_flag = false
            elseif isempty(m._stack)
                continue_flag = false
            end

        end
    end

    return nothing
end

"""
$(SIGNATURES)

Checks to see if current node should be reprocessed.
"""
repeat_check(t::ExtensionType, m::Optimizer) = false

relative_gap(L::Float64, U::Float64) = ((L > -Inf) && (U < Inf)) ?  abs(U - L)/(max(abs(L), abs(U))) : Inf
relative_tolerance(L::Float64, U::Float64, tol::Float64) = relative_gap(L, U)  > tol || ~(L > -Inf)

"""
$(SIGNATURES)

Checks for termination of algorithm due to satisfying absolute or relative
tolerance, infeasibility, or a specified limit, returns a boolean valued true
if algorithm should continue.
"""
function termination_check(t::ExtensionType, m::Optimizer)

    node_in_stack = length(m._stack)
    L = m._global_lower_bound
    U = m._global_upper_bound

    if node_in_stack === 0

        if m._first_solution_node > 0
            m._termination_status_code = MOI.OPTIMAL
            m._result_status_code = MOI.FEASIBLE_POINT
            (m.verbosity >= 3) && println("Empty Stack: Exhaustive Search Finished")

        else
            m._termination_status_code = MOI.INFEASIBLE
            m._result_status_code = MOI.INFEASIBILITY_CERTIFICATE
            (m.verbosity >= 3) && println("Empty Stack: Infeasible")
        end

    elseif node_in_stack >= m.node_limit

        m._termination_status_code = MOI.NODE_LIMIT
        m._result_status_code = MOI.UNKNOWN_RESULT_STATUS
        (m.verbosity >= 3) && println("Node Limit Exceeded")

    elseif m._iteration_count >= m.iteration_limit

        m._termination_status_code = MOI.ITERATION_LIMIT
        m._result_status_code = MOI.UNKNOWN_RESULT_STATUS
        (m.verbosity >= 3) && println("Maximum Iteration Exceeded")

    elseif ~relative_tolerance(L, U, m.relative_tolerance)

        m._termination_status_code = MOI.OPTIMAL
        m._result_status_code = MOI.FEASIBLE_POINT
        (m.verbosity >= 3) && println("Relative Tolerance Achieved")

    elseif (U - L) < m.absolute_tolerance

        m._termination_status_code = MOI.OPTIMAL
        m._result_status_code = MOI.FEASIBLE_POINT
        (m.verbosity >= 3) && println("Absolute Tolerance Achieved")

    elseif m._run_time > m.time_limit

        m._termination_status_code = MOI.TIME_LIMIT
        m._result_status_code = MOI.UNKNOWN_RESULT_STATUS
        (m.verbosity >= 3) && println("Time Limit Exceeded")

    else

        return false
    end

    return true
end

"""
$(SIGNATURES)

Checks for convergence of algorithm with respect to absolute and/or relative
tolerances.
"""
function convergence_check(t::ExtensionType, m::Optimizer)

  L = m._lower_objective_value
  U = m._global_upper_bound
  t = (U - L) <= m.absolute_tolerance
  if (U < Inf) && (L > Inf)
      t |= (abs(U - L)/(max(abs(L), abs(U))) <= m.relative_tolerance)
  end

  if t && m._min_converged_value < Inf
      m._min_converged_value = min(m._min_converged_value, L)
  else
      m._min_converged_value = L
  end

  return t
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

Takes an `MOI.TerminationStatusCode` and a `MOI.ResultStatusCode` and returns `true`
if this corresponds to a solution that is proven to be feasible.
Returns `false` otherwise.
"""
function is_feasible_solution(t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)

    termination_flag = false
    result_flag = false

    (t === MOI.OPTIMAL) && (termination_flag = true)
    (t === MOI.LOCALLY_SOLVED) && (termination_flag = true)

    # This is default solver specific... the acceptable constraint tolerances
    # are set to the same values as the basic tolerance. As a result, an
    # acceptably solved solution is feasible but non necessarily optimal
    # so it should be treated as a feasible point
    if (t === MOI.ALMOST_LOCALLY_SOLVED) && (r === MOI.NEARLY_FEASIBLE_POINT)
        termination_flag = true
        result_flag = true
    end

    (r === MOI.FEASIBLE_POINT) && (result_flag = true)

    return (termination_flag && result_flag)
end

"""
$(SIGNATURES)

Retrieves the lower and upper duals for variable bounds from the
`relaxed_optimizer` and sets the appropriate values in the
`_lower_lvd` and `_lower_uvd` storage fields.
"""
function set_dual!(m::Optimizer)

    relaxed_optimizer = m.relaxed_optimizer
    relaxed_variable_lt = m._relaxed_variable_lt
    relaxed_variable_gt = m._relaxed_variable_gt

    for i = 1:m._working_problem._var_leq_count
        ci_lt, i_lt = @inbounds relaxed_variable_lt[i]
        @inbounds m._lower_uvd[i_lt] = MOI.get(relaxed_optimizer, MOI.ConstraintDual(), ci_lt)
    end
    for i = 1:m._working_problem._var_geq_count
        ci_gt, i_gt = @inbounds relaxed_variable_gt[i]
        @inbounds m._lower_lvd[i_gt] = MOI.get(relaxed_optimizer, MOI.ConstraintDual(), ci_gt)
    end

    return nothing
end

"""
$(SIGNATURES)

Runs interval, linear, quadratic contractor methods followed by obbt and a
constraint programming walk up to tolerances specified in
`EAGO.Optimizer` object.
"""
function preprocess!(t::ExtensionType, m::Optimizer)

    reset_relaxation!(m)

    wp = m._working_problem
    params = m._parameters

    # Sets initial feasibility
    feasible_flag = true
    m._obbt_performed_flag = false

    # compute initial volume
    m._initial_volume = prod(upper_variable_bounds(m._current_node) -
                             lower_variable_bounds(m._current_node))

    if params.fbbt_lp_depth >= m._iteration_count
        load_fbbt_buffer!(m)
        for i = 1:m.fbbt_lp_repetitions
            if feasible_flag
                for j = 1:wp._saf_leq_count
                    !feasible_flag && break
                    saf_leq =  wp._saf_leq[j]
                    feasible_flag &= fbbt!(m, saf_leq)
                end
                !feasible_flag && break

                for j = 1:wp._saf_eq_count
                    !feasible_flag && break
                    saf_eq = wp._saf_eq[j]
                    feasible_flag &= fbbt!(m, saf_eq)
                end
                !feasible_flag && break
            end
        end
        unpack_fbbt_buffer!(m)
    end

    # done after cp to prevent using cp specific flags in cut generation
    set_first_relax_point!(m)

    cp_walk_count = 0
    perform_cp_walk_flag = feasible_flag
    perform_cp_walk_flag &= (params.cp_depth >= m._iteration_count)
    perform_cp_walk_flag &= (cp_walk_count < m.cp_repetitions)
    while perform_cp_walk_flag
        feasible_flag &= set_constraint_propagation_fbbt!(m)
        !feasible_flag && break
        cp_walk_count += 1
        perform_cp_walk_flag = (cp_walk_count < m.cp_repetitions)
    end

    obbt_count = 0
    perform_obbt_flag = feasible_flag
    perform_obbt_flag &= (params.obbt_depth >= m._iteration_count)
    perform_obbt_flag &= (obbt_count < m.obbt_repetitions)

    while perform_obbt_flag
        feasible_flag &= obbt!(m)
        m._obbt_performed_flag = true
        !feasible_flag && break
        obbt_count += 1
        perform_obbt_flag     = (obbt_count < m.obbt_repetitions)
    end

    m._final_volume = prod(upper_variable_bounds(m._current_node) -
                           lower_variable_bounds(m._current_node))

    m._preprocess_feasibility = feasible_flag

    return nothing
end

"""
$(SIGNATURES)

Updates the relaxed constraint by setting the constraint set of `v == x*`` ,
`xL_i <= x_i`, and `x_i <= xU_i` for each such constraint added to the relaxed
optimizer.
"""
function update_relaxed_problem_box!(m::Optimizer)
    opt = m.relaxed_optimizer
    lb = m._current_node.lower_variable_bounds
    ub = m._current_node.upper_variable_bounds
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], ET(lb[x[2]])), m._relaxed_variable_eq)
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], LT(ub[x[2]])), m._relaxed_variable_lt)
    foreach(x -> MOI.set(opt, MOI.ConstraintSet(), x[1], GT(ub[x[2]])), m._relaxed_variable_gt)
    return nothing
end

function interval_objective_bound(m::Optimizer, n::NodeBB)

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

"""
$(SIGNATURES)

A fallback lower bounding problem that consists of an natural interval extension
calculation. This is called when the optimizer used to compute the lower bound
does not return a termination and primal status code indicating that it
successfully solved the relaxation to a globally optimal point.
"""
function fallback_interval_lower_bound!(m::Optimizer, n::NodeBB)

    feasible_flag = true

    if !cp_condition(m)
        for i = 1:m._working_problem._saf_leq_count
            saf_leq =  m._working_problem._saf_leq[i]
            feasible_flag &= (lower_interval_bound(m, saf_leq, n) <= 0.0)
            !feasible_flag && break
        end

        if feasible_flag
            for i = 1:m._working_problem._saf_eq_count
                saf_eq =  m._working_problem._saf_eq[i]
                lower_value, upper_value = interval_bound(m, saf_eq, n)
                feasible_flag &= (lower_value <= 0.0 <= upper_value)
                !feasible_flag && break
            end
        end

        if feasible_flag
            for i = 1:m._working_problem._sqf_leq_count
                sqf_leq =  m._working_problem._sqf_leq[i]
                feasible_flag &= (lower_interval_bound(m, sqf_leq, n) <= 0.0)
                !feasible_flag && break
            end
        end

        if feasible_flag
            for i = 1:m._working_problem._sqf_eq_count
                sqf_eq =  m._working_problem._sqf_eq[i]
                lower_value, upper_value = interval_bound(m, sqf_eq, n)
                feasible_flag &= (lower_value <= 0.0 <= upper_value)
                !feasible_flag && break
            end
        end

        if feasible_flag
            for i = 1:m._working_problem._nonlinear_count
                nl_constr =  m._working_problem._nonlinear_constr[i]
                lower_value, upper_value = interval_bound(m, nl_constr, n)
                feasible_flag &= upper_value < nl_constr.lower_bound
                feasible_flag &= lower_value > nl_constr.upper_bound
                !feasible_flag && break
            end
        end
    end

    if feasible_flag
        interval_objective_used = interval_objective_bound(m, n)
        @__dot__ m._current_xref = 0.5*(n.upper_variable_bounds + n.lower_variable_bounds)
        unsafe_check_fill!(isnan, m._current_xref, 0.0, length(m._current_xref))
    else
        m._lower_objective_value = -Inf
    end
    m._lower_feasibility = feasible_flag

    return nothing
end

function intrepret_relaxed_solution!(m::Optimizer, d::T) where T
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
function lower_problem!(t::ExtensionType, m::Optimizer)

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
    m._working_problem._objective_nl.has_value = false
    m._working_problem._relaxed_evaluator.interval_intersect = false

    if !m._obbt_performed_flag
        relax_constraints!(m, 1)
    end
    relax_objective!(m, 1)

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
function cut_update!(m::Optimizer)

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
function cut_condition(t::ExtensionType, m::Optimizer)

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
        objective_lo = -Inf
        obj_type = m._working_problem._objective_type
        if obj_type === SINGLE_VARIABLE
            var_index = m._working_problem._objective_sv.variable.value
            if m._branch_variables[var_index]
                obj_indx = m._sol_to_branch_map[var_index]
                lower_variable_bnd = n.lower_variable_bounds[obj_indx]
                if !isinf(lower_variable_bnd)
                    objective_lo = lower_variable_bnd
                end
            end

        elseif obj_type === SCALAR_AFFINE
            objective_lo = lower_interval_bound(m, m._working_problem._objective_saf_parsed, n)

        elseif obj_type === SCALAR_QUADRATIC
            objective_lo = lower_interval_bound(m, m._working_problem._objective_sqf, n)

        elseif obj_type === NONLINEAR
            objective_lo = lower_interval_bound(m, m._working_problem._objective_nl, n)

        end

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
function add_cut!(t::ExtensionType, m::Optimizer)

    fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)
    m._working_problem._relaxed_evaluator.is_first_eval = true
    m._working_problem._relaxed_evaluator.is_intersect = false
    m._new_eval_objective = true
    m._new_eval_constraint = true

    relax_constraints!(m, m._cut_iterations)
    relax_objective!(m, m._cut_iterations)

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

"""
$(SIGNATURES)

Default check to see if the upper bounding problem should be run. By default,
The upper bounding problem is run on every node up to depth `upper_bounding_depth`
and is triggered with a probability of `0.5^(depth - upper_bounding_depth)`
afterwards.
"""
function default_nlp_heurestic(m::Optimizer)
    bool = false
    ubd_limit = m.upper_bounding_depth
    depth = m._current_node.depth
    bool |= (depth <= ubd_limit)
    bool |= (rand() < 0.5^(depth - m.upper_bounding_depth))
    return bool
end

"""
$(SIGNATURES)

Default upper bounding problem which simply calls `solve_local_nlp!` to solve
the nlp locally.
"""
function upper_problem!(t::ExtensionType, m::Optimizer)

    if !default_nlp_heurestic(m)
        m._upper_feasibility = false
        m._upper_objective_value = Inf

    else
        solve_local_nlp!(m)

    end

    return nothing
end


"""
$(SIGNATURES)

Default postprocess perfoms duality-based bound tightening on the `y`.
"""
function postprocess!(t::ExtensionType, m::Optimizer)

    if m.dbbt_depth > m._iteration_count
        variable_dbbt!(m._current_node, m._lower_lvd, m._lower_uvd,
                       m._lower_objective_value, m._global_upper_bound,
                       m._branch_variable_num)
    end

    return nothing
end

"""
$(SIGNATURES)

Provides a hook for extensions to EAGO as opposed to standard global, local,
or linear solvers.
"""
optimize_hook!(t::ExtensionType, m::Optimizer) = nothing

function store_candidate_solution!(m::Optimizer)

    if m._upper_feasibility && (m._upper_objective_value < m._global_upper_bound)

        m._feasible_solution_found = true
        m._first_solution_node = m._maximum_node_id
        m._solution_value = m._upper_objective_value
        m._global_upper_bound = m._upper_objective_value
        @__dot__ m._solution = m._upper_solution

    end
    return nothing
end

function set_global_lower_bound!(m::Optimizer)

    if !isempty(m._stack)

        min_node = minimum(m._stack)
        lower_bound = min_node.lower_bound
        if m._global_lower_bound < lower_bound
            m._global_lower_bound = lower_bound
        end

    end

    return nothing
end

# wraps subroutine call to isolate ExtensionType
parse_global!(m::Optimizer) = parse_global!(m.ext_type, m)
presolve_global!(m::Optimizer) = presolve_global!(m.ext_type, m)
termination_check(m::Optimizer) = termination_check(m.ext_type, m)
cut_condition(m::Optimizer) = cut_condition(m.ext_type, m)
convergence_check(m::Optimizer) = convergence_check(m.ext_type, m)
repeat_check(m::Optimizer) = repeat_check(m.ext_type, m)
node_selection!(m::Optimizer) = node_selection!(m.ext_type, m)
preprocess!(m::Optimizer) = preprocess!(m.ext_type, m)
lower_problem!(m::Optimizer) = lower_problem!(m.ext_type, m)
add_cut!(m::Optimizer) = add_cut!(m.ext_type, m)
upper_problem!(m::Optimizer) = upper_problem!(m.ext_type, m)
postprocess!(m::Optimizer) = postprocess!(m.ext_type, m)
single_storage!(m::Optimizer) = single_storage!(m.ext_type, m)
branch_node!(m::Optimizer) = branch_node!(m.ext_type, m)
fathom!(m::Optimizer) = fathom!(m.ext_type, m)
revert_adjusted_upper_bound!(m::Optimizer) = revert_adjusted_upper_bound!(m.ext_type, m)

"""
$(TYPEDSIGNATURES)

Solves the branch and bound problem with the input EAGO optimizer object.
"""
function global_solve!(m::Optimizer)

    m._iteration_count = 1
    m._node_count = 1

    parse_global!(m)
    presolve_global!(m)

    logging_on = m.log_on
    verbosity = m.verbosity

    # terminates when max nodes or iteration is reach, or when node stack is empty
    while !termination_check(m)

        # Selects node, deletes it from stack, prints based on verbosity
        node_selection!(m)
        (verbosity >= 3) && print_node!(m)

        # Performs prepocessing and times
        logging_on && (start_time = time())
        preprocess!(m)
        if logging_on
            m._last_preprocess_time = time() - start_time
        end

        if m._preprocess_feasibility

            # solves & times lower bounding problem
            logging_on && (start_time = time())
            m._cut_iterations = 1
            lower_problem!(m)
            while cut_condition(m)
                add_cut!(m)
            end
            if logging_on
                m._last_lower_problem_time = time() - start_time
            end
            print_results!(m, true)
            print_results_post_cut!(m)

            # checks for infeasibility stores solution
            if m._lower_feasibility
                if !convergence_check(m)

                    logging_on && (start_time = time())
                    upper_problem!(m)
                    if logging_on
                        m._last_upper_problem_time = time() - start_time
                    end
                    print_results!(m, false)
                    store_candidate_solution!(m)
                    if m._input_problem._optimization_sense === MOI.FEASIBILITY_SENSE
                        if !m.feasible_local_continue || m.local_solve_only
                            break
                        end
                    end

                    # Performs and times post processing
                    logging_on && (start_time = time())
                    postprocess!(m)
                    if logging_on
                        m._last_postprocessing_time = time() - start_time
                    end

                    # Checks to see if the node
                    if m._postprocess_feasibility
                        if repeat_check(m)
                            single_storage!(m)
                        else
                            branch_node!(m)
                        end
                    end
                end
            end
            fathom!(m)
        else
            m._lower_objective_value = -Inf
            m._lower_feasibility = false
            m._upper_feasibility = false
        end
        set_global_lower_bound!(m)
        m._run_time = time() - m._start_time
        m._time_left = m.time_limit - m._run_time
        log_iteration!(m)
        print_iteration!(m)
        m._iteration_count += 1
    end

    revert_adjusted_upper_bound!(m)
    m._objective_value = m._global_upper_bound

    # Prints the solution
    print_solution!(m)

    return nothing
end

optimize!(::Val{MINCVX}, m::Optimizer) = global_solve!(m)


#=
TODO: SET SOLVE AT END OF NONCONVEX
if MOI.get(opt, MOI.ObjectiveSense()) == MOI.MIN_SENSE
    m._global_lower_bound = MOI.get(opt, MOI.ObjectiveBound())
    m._global_upper_bound = MOI.get(opt, MOI.ObjectiveValue())
else
    m._global_lower_bound = MOI.get(opt, MOI.ObjectiveValue())
    m._global_upper_bound = MOI.get(opt, MOI.ObjectiveBound())
end
=#
