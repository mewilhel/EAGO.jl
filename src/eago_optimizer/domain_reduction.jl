# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/domain_reduction.jl
# Contains subroutines used for domain reduction.
#############################################################################
"""
    trivial_filtering!

Excludes OBBT on variable indices that are tight for the solution of the relaxation.
"""
function trivial_filtering!(m::GlobalOptimizer)

    m._preprocess_termination_status = MOI.get(m.relaxed_optimizer, MOI.TerminationStatus())
    m._preprocess_result_status = MOI.get(m.relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(m._preprocess_termination_status,
                                                    m._preprocess_result_status)

    if valid_flag && feasible_flag
        @inbounds for j = 1:length(m._obbt_working_lower_index)
            if m._obbt_working_lower_index[j]
                vi = m._relaxed_variable_index[j]
                diff = MOI.get(m.relaxed_optimizer, MOI.VariablePrimal(), vi)
                diff -= _lower_bound(FullVar, m, j)
                if abs(diff) <= _obbt_tolerance(m)
                    m._obbt_working_lower_index[j] = false
                end
            end
        end
        @inbounds for j = 1:length(m._obbt_working_upper_index)
            if m._obbt_working_upper_index[j]
                vi = m._relaxed_variable_index[j]
                diff = -MOI.get(m.relaxed_optimizer, MOI.VariablePrimal(), vi)
                diff += _upper_bound(FullVar, m, j)
                if abs(diff) <= _obbt_tolerance(m)
                    m._obbt_working_upper_index[j] = false
                end
            end
        end
    end

    return
end

"""
    bool_indx_diff

Utility function used to set vector of booleans z to x & ~y. Avoids the
generation of conversion of the BitArray created by broadcasting logical operators.
"""
function bool_indx_diff(z::Vector{Bool}, x::Vector{Bool}, y::Vector{Bool})
    for i = 1:length(z)
        @inbounds z[i] = (x[i] & ~y[i])
    end
    return
end

"""
    aggressive_filtering!

Excludes OBBT on variable indices after a search in a filtering direction.
"""
function aggressive_filtering!(m::GlobalOptimizer)

    # Initial filtering vector (negative one direction per remark in Gleixner2017)
    variable_number = m._working_problem._variable_num
    v = -ones(variable_number)

    # Copy prior index set (ignores linear and binary terms)
    obbt_variable_num = m._obbt_variable_num
    copyto!(m._old_low_index, m._obbt_working_lower_index)
    copyto!(m._old_upp_index, m._obbt_working_upper_index)
    copyto!(m._new_low_index, m._obbt_working_lower_index)
    copyto!(m._new_upp_index, m._obbt_working_upper_index)

    # Exclude unbounded directions
    for i = 1:obbt_variable_num
        m._new_low_index[i] &= !(_lower_bound(FullVar, m, i) == -Inf)
        m._new_upp_index[i] &= !(_upper_bound(FullVar, m, i) == Inf)
    end

    # Begin the main algorithm
    for k = 1:m.obbt_aggressive_max_iteration

        # Set index differences and vector for filtering direction
        bool_indx_diff(m._lower_indx_diff, m._old_low_index, m._new_low_index)
        bool_indx_diff(m._upper_indx_diff, m._old_upp_index, m._new_upp_index)

        for i = 1:obbt_variable_num
            if m._lower_indx_diff[i] && (v[i] < 0.0)
                v[i] = 0.0
            elseif m._upper_indx_diff[i] && (v[i] > 0.0)
                v[i] = 0.0
            end
        end

        # Termination Condition
        ((~any(m._new_low_index) & ~any(m._new_upp_index)) || (iszero(v))) && break
        if k >= 2
            if (count(m._lower_indx_diff) + count(m._upper_indx_diff)) < m.obbt_aggressive_min_dimension
                break
            end
        end

        # Set objective in OBBT problem to filtering vector
        MOI.set(m.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        saf = SAF(SAT.(v, m._relaxed_variable_index), 0.0)
        MOI.set(m.relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), saf)

        # Optimizes the problem and if successful filter additional bounds
        MOI.optimize!(m.relaxed_optimizer)

        m._preprocess_termination_status = MOI.get(m.relaxed_optimizer, MOI.TerminationStatus())
        m._preprocess_result_status = MOI.get(m.relaxed_optimizer, MOI.PrimalStatus())
        valid_flag, feasible_flag = is_globally_optimal(m._preprocess_termination_status,
                                                        m._preprocess_result_status)

        if valid_flag && feasible_flag
            variable_primal = MOI.get(m.relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index)
            copyto!(m._new_low_index, m._old_low_index)
            copyto!(m._new_upp_index, m._old_upp_index)
            for i = 1:obbt_variable_num
                vp_value = variable_primal[i]
                if m._old_low_index[i] && vp_value == m.lower_variable_bound[i]
                    m._new_low_index[i] = false
                end
                if m._old_upp_index[i] && vp_value == m.upper_variable_bound[i]
                    m._new_upp_index[i] = false
                end
            end
        else
            return false
        end
    end
    copyto!(m._obbt_working_lower_index, m._new_low_index)
    copyto!(m._obbt_working_upper_index, m._new_upp_index)
    return true
end

"""
    set_node_flag!
"""
function set_node_flag!(m::GlobalOptimizer)
    foreach(set_node_flag!, m._working_problem._nonlinear_constr)
    return nothing
end

"""
    set_reference_point!
"""
function set_reference_point!(m::GlobalOptimizer)

    evaluator = m._working_problem._relaxed_evaluator
    evaluator_x = evaluator.x
    current_xref = m._current_xref

    new_reference_point = false
    for node_i = 1:m._branch_variable_num
        solution_i = m._branch_to_sol_map[node_i]

        node_x = current_xref[node_i]
        solution_x = evaluator_x[solution_i]

        if node_x !== solution_x
            evaluator_x[solution_i] = node_x
            new_reference_point = true
        end
    end

    if new_reference_point
        foreach(x -> (x.has_value = false;), m._working_problem._nonlinear_constr)
    end
    fill!(evaluator.subexpressions_eval, false)

    return nothing
end

"""
    obbt!

Performs OBBT with filtering and greedy ordering as detailed in:
Gleixner, A.M., Berthold, T., Müller, B. et al. J Glob Optim (2017) 67: 731.
https://doi.org/10.1007/s10898-016-0450-4
"""
function obbt!(m::GlobalOptimizer)

    feasibility = true

    # set node and reference point if necessary then solve initial problem to
    # feasibility. This is repeated `obbt_repetitions` number of times in the
    # following fashion. Relax the problem, populate affine constraints, run
    # obbt which contracts variable bounds, delete affine constraints...
    # update variable bounds and repeat. TODO: Keep track of which variables
    # participate in which functions and only delete a constraint if a
    # variable participating in a nonlinear term changes it bounds.
    if m._obbt_performed_flag
        reset_relaxation!(m)
        set_first_relax_point!(m)
    end

    update_relaxed_problem_box!(m)
    if m._nonlinear_evaluator_created
        set_node!(m._working_problem._relaxed_evaluator, n)
        set_node_flag!(m)
        set_reference_point!(m)
    end
    relax_constraints!(m, 1)
    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(relaxed_optimizer)

    # Sets indices to attempt OBBT on
    obbt_variable_num = m._obbt_variable_num
    fill!(m._obbt_working_lower_index, true)
    fill!(m._obbt_working_upper_index, true)

    # Filters out any indicies with active bounds on variables
    # determined by solving the feasibility problem
    trivial_filtering!(m, n)

    # Applies an aggressive filter to eliminate indices that
    # cannot be tightened by obbt
    if m._obbt_aggressive_on
        feasibility = aggressive_filtering!(m, n)
    end

    # extracts info from relaxed problem (redundant if aggressive obbt is called)
    m._preprocess_termination_status = MOI.get(relaxed_optimizer, MOI.TerminationStatus())
    m._preprocess_result_status = MOI.get(relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(m._preprocess_termination_status,
                                                    m._preprocess_result_status)

    if valid_flag && feasible_flag
        xLP = MOI.get(relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index)
    else
        return false
    end

    # continue tightening bounds by optimization until all indices have been checked
    # or the node is empty and the problem is thus proven infeasible
    while (any(m._obbt_working_lower_index) || any(m._obbt_working_upper_index)) && !isempty(n)

        # Get lower value
        lower_indx = -1;     upper_indx = -1
        lower_value = Inf;   upper_value = Inf

        # min of xLP - yL on active
        if any(m._obbt_working_lower_index)
            for i = 1:obbt_variable_num
                if m._obbt_working_lower_index[i]
                    temp_value = _lower_solution(FullVar, m, i) - _lower_bound(FullVar, m, i)
                    # Need less than or equal to handle unbounded cases
                    if temp_value <= lower_value
                        lower_value = temp_value
                        lower_indx = i
                    end
                end
            end
        end

        # min of yU - xLP on active
        if any(m._obbt_working_upper_index)
            for i = 1:obbt_variable_num
                if m._obbt_working_upper_index[i]
                    temp_value = _upper_bound(FullVar, m, i) - _lower_solution(FullVar, m, i)
                    if temp_value <= upper_value
                        upper_value = temp_value
                        upper_indx = i
                    end
                end
            end
        end

        # default to upper bound if no lower bound is found, use maximum distance otherwise
        if lower_value <= upper_value && lower_indx > 0

            m._obbt_working_lower_index[lower_indx] = false
            var = SV(m._relaxed_variable_index[lower_indx])

            MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(relaxed_optimizer, MOI.ObjectiveFunction{SV}(), var)

            MOI.optimize!(m.relaxed_optimizer)
            m._preprocess_termination_status = MOI.get(relaxed_optimizer, MOI.TerminationStatus())
            m._preprocess_result_status = MOI.get(relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(m._preprocess_termination_status,
                                                            m._preprocess_result_status)

            if valid_flag && feasible_flag
                xLP .= MOI.get(relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index)

                node_index = branch_to_sol_map[lower_indx]
                updated_value = xLP[node_index]
                previous_value = n.lower_variable_bound[lower_indx]

                # if bound is improved update node and corresponding constraint update
                # the node bounds and the single variable bound in the relaxation
                # we assume branching does not occur on fixed variables and interval
                # constraints are internally bridged by EAGO. So the only L <= x
                # constraint in the model is a GreaterThan.
                if updated_value > previous_value && (updated_value - previous_value) > 1E-6
                    ci_list = MOI.get(relaxed_optimizer, MOI.ListOfConstraintIndices{SAF,LT}())
                    sv_geq_ci = m._node_to_sv_geq_ci[lower_indx]
                    MOI.set(relaxed_optimizer, MOI.ConstraintSet(), sv_geq_ci, GT(updated_value))
                    n.lower_variable_bound[lower_indx] = updated_value
                end

                if isempty(n)
                    feasibility = false
                    break
                end

            elseif valid_flag && !feasible_flag
                feasibility = false
                break

            else
                break
            end

        elseif upper_indx > 0

            m._obbt_working_upper_index[upper_indx] = false
            var = SV(m._relaxed_variable_index[upper_indx])
            MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
            MOI.set(relaxed_optimizer, MOI.ObjectiveFunction{SV}(), var)
            MOI.optimize!(relaxed_optimizer)
            m._preprocess_termination_status = MOI.get(m.relaxed_optimizer, MOI.TerminationStatus())
            m._preprocess_result_status = MOI.get(m.relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(m._preprocess_termination_status,
                                                            m._preprocess_result_status)

            if valid_flag && feasible_flag
                xLP .= MOI.get(relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index)
                node_index = branch_to_sol_map[upper_indx]
                updated_value = xLP[node_index]
                previous_value = n.upper_variable_bound[upper_indx]

                # if bound is improved update node and corresponding constraint update
                # the node bounds and the single variable bound in the relaxation
                # we assume branching does not occur on fixed variables and interval
                # constraints are internally bridged by EAGO. So the only U => x
                # constraint in the model is a LessThan.
                if updated_value < previous_value && (previous_value - updated_value) > 1E-6
                    sv_leq_ci = m._node_to_sv_leq_ci[upper_indx]
                    MOI.set(relaxed_optimizer, MOI.ConstraintSet(), sv_leq_ci, LT(updated_value))
                    n.upper_variable_bound[upper_indx] = updated_value
                end

                if isempty(n)
                    feasibility = false
                    break
                end

            elseif valid_flag && !feasible_flag
                feasibility = false
                break

            else
                break
            end

        else
            break
        end
        trivial_filtering!(m)
    end

    return feasibility
end

"""
    fbbt!

Performs feasibility-based bound tightening on a back-end constraint and returns `true` if it is feasible or
`false` if it is infeasible.
"""
function fbbt! end

function fbbt!(::Type{Val{:leq}}, m::Optimizer, f::AffineFunction)

    # compute full sum
    temp_sum = -f.constant
    @inbounds for k = 1:f.len
        aik, indx_k = f.terms[k]
        if aik !== 0.0
            aik_xL = aik*_lower_bound(FullVar, m, indx_k)
            aik_xU = aik*_upper_bound(FullVar, m, indx_k)
            temp_sum -= min(aik_xL, aik_xU)
        end
    end

    # subtract extra term, check to see if implied bound is better, if so update the node and
    # the working sum if the node is now empty then break
    @inbounds for k = 1:f.len
        aik, indx_k = f.terms[k]
        if aik !== 0.0
            xL = _lower_bound(FullVar, m, indx_k)
            xU = _upper_bound(FullVar, m, indx_k)
            aik_xL = aik*xL
            aik_xU = aik*xU
            temp_sum += min(aik_xL, aik_xU)
            xh = temp_sum/aik

            if aik > 0.0
                (xh < xL) && return false
                if xh > xL
                    _set_upper_bound(FullVar, m, indx_k, xh)
                end
            elseif aik < 0.0
                (xh > xU) && return false
                if xh < xU
                    _set_lower_bound(FullVar, m, indx_k, xh)
                end
            else
                temp_sum -= min(aik_xL, aik_xU)
                continue
            end

            aik_xL = aik*_lower_bound(FullVar, m, indx_k)
            aik_xU = aik*_upper_bound(FullVar, m, indx_k)
            temp_sum -= min(aik_xL, aik_xU)
        end
    end

    return true
end

function fbbt!(::Type{Val{:eq}}, m::Optimizer, f::AffineFunction)

    # compute full sum
    terms = f.terms
    temp_sum_leq = -f.constant
    temp_sum_geq = -f.constant

    @inbounds for k = 1:f.len
        aik, indx_k = terms[k]
        if aik !== 0.0
            aik_xL = aik*_lower_bound(FullVar, m, indx_k)
            aik_xU = aik*_upper_bound(FullVar, m, indx_k)
            temp_sum_leq -= min(aik_xL, aik_xU)
            temp_sum_geq -= max(aik_xL, aik_xU)
        end
    end

    # subtract extra term, check to see if implied bound is better, if so update the node and
    # the working sum if the node is now empty then break
    @inbounds for k = 1:f.len

        aik, indx_k = terms[k]
        if aik !== 0.0

            xL = _lower_bound(FullVar, m, indx_k)
            xU = _upper_bound(FullVar, m, indx_k)

            aik_xL = aik*xL
            aik_xU = aik*xU

            temp_sum_leq += min(aik_xL, aik_xU)
            temp_sum_geq += max(aik_xL, aik_xU)

            xh_leq = temp_sum_leq/aik
            xh_geq = temp_sum_geq/aik

            if aik > 0.0
                (xh_leq < xL) && return false
                if xh_leq > xL
                    upper_bounds[indx_k] = xh_leq
                end
                (xh_geq > xU) && return false
                if (xh_geq < xU)
                    lower_bounds[indx_k] = xh_geq
                end

            elseif aik < 0.0
                (xh_leq > xU) && return false
                if xh_leq < xU
                    lower_bounds[indx_k] = xh_leq
                end
                (xh_geq < xL) && return false
                if (xh_geq > xL)
                    upper_bounds[indx_k] = xh_geq
                end

            else
                temp_sum_leq -= min(aik_xL, aik_xU)
                temp_sum_geq -= max(aik_xL, aik_xU)
                continue

            end
            aik_xL = aik*_lower_bound(FullVar, m, indx_k)
            aik_xU = aik*_upper_bound(FullVar, m, indx_k)

            temp_sum_leq -= min(aik_xL, aik_xU)
            temp_sum_geq -= max(aik_xL, aik_xU)
        end
    end

    return true
end

cp_condition(m::Optimizer) = false

"""
    set_constraint_propagation_fbbt!

Performs bound tightening based on forward/reverse interval and/or McCormick passes. This routine
resets the current node with new interval bounds.
"""
function set_constraint_propagation_fbbt!(m::Optimizer)
    feasible_flag = true

    if m._nonlinear_evaluator_created
        evaluator = m._working_problem._relaxed_evaluator
        set_node!(m._working_problem._relaxed_evaluator, m._current_node)
        set_node_flag!(m)
        set_reference_point!(m)

        m._working_problem._relaxed_evaluator.is_first_eval = m._new_eval_constraint
        for constr in m._working_problem._nonlinear_constr
            if feasible_flag
                forward_pass!(evaluator, constr)
                feasible_flag &= reverse_pass!(evaluator, constr)
                evaluator.interval_intersect = true
            end
        end

        m._working_problem._relaxed_evaluator.is_first_eval = m._new_eval_constraint
        for constr in m._working_problem._nonlinear_constr
            if feasible_flag
                forward_pass!(evaluator, constr)
            end
        end

        evaluator.is_post = m.subgrad_tighten

        m._new_eval_constraint = false
        m._new_eval_objective = false

        retrieve_x!(m._current_xref, evaluator)
        m._current_node = retrieve_node(evaluator)
    end

    return feasible_flag
end
