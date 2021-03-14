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
    variable_dbbt!

Tightens the bounds of the `_current_node` using the current global upper bound
and the duality information obtained from the relaxation.
"""
function variable_dbbt!(n::NodeBB, mult_lo::Vector{Float64}, mult_hi::Vector{Float64},
                        LBD::Float64, UBD::Float64, nx::Int64)

    delta = UBD - LBD
    lvbs = n.lower_variable_bounds
    uvbs = n.upper_variable_bounds
    if LBD <= UBD
        for i = 1:nx
            ml = @inbounds mult_lo[i]
            if ml > 0.0
                cut = @inbounds lvbs[i] + delta/ml
                if cut < @inbounds uvbs[i]
                    @inbounds uvbs[i] = cut
                end
            else
                mh = @inbounds mult_hi[i]
                if mh > 0.0
                    cut = @inbounds uvbs[i] - delta/mh
                    if cut > @inbounds lvbs[i]
                        @inbounds lvbs[i] = cut
                    end
                end
            end
         end
    end

    return nothing
end

"""
    trivial_filtering!

Excludes OBBT on variable indices that are tight for the solution of the relaxation.
"""
function trivial_filtering!(m::Optimizer, n::NodeBB)

    obbt_tolerance = m.obbt_tolerance
    m._preprocess_termination_status = MOI.get(m.relaxed_optimizer, MOI.TerminationStatus())
    m._preprocess_result_status = MOI.get(m.relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(m._preprocess_termination_status,
                                                    m._preprocess_result_status)

    if valid_flag && feasible_flag
        for j = 1:length(m._obbt_working_lower_index)
            if @inbounds m._obbt_working_lower_index[j]
                vi = @inbounds m._relaxed_variable_index[j]
                diff = MOI.get(m.relaxed_optimizer, MOI.VariablePrimal(), vi)
                diff -= @inbounds n.lower_variable_bounds[j]
                if abs(diff) <= obbt_tolerance
                    @inbounds m._obbt_working_lower_index[j] = false
                end
            end
        end
        for j = 1:length(m._obbt_working_upper_index)
            if @inbounds m._obbt_working_upper_index[j]
                vi = @inbounds m._relaxed_variable_index[j]
                diff = -MOI.get(m.relaxed_optimizer, MOI.VariablePrimal(), vi)
                diff += @inbounds n.upper_variable_bounds[j]
                if abs(diff) <= obbt_tolerance
                    @inbounds m._obbt_working_upper_index[j] = false
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
function aggressive_filtering!(m::Optimizer, n::NodeBB)

    # Initial filtering vector (negative one direction per remark in Gleixner2017)
    variable_number = m._working_problem._variable_count
    v = -ones(variable_number)

    # Copy prior index set (ignores linear and binary terms)
    obbt_variable_count = m._obbt_variable_count
    copyto!(m._old_low_index, m._obbt_working_lower_index)
    copyto!(m._old_upp_index, m._obbt_working_upper_index)
    copyto!(m._new_low_index, m._obbt_working_lower_index)
    copyto!(m._new_upp_index, m._obbt_working_upper_index)

    # Exclude unbounded directions
    for i = 1:obbt_variable_count
        if @inbounds m._new_low_index[i] && @inbounds n.lower_variable_bounds[i] === -Inf
            @inbounds m._new_low_index[i] = false
        end
        if @inbounds m._new_low_index[i] && @inbounds n.upper_variable_bounds[i] === Inf
            @inbounds m._new_low_index[i] = false
        end
    end

    # Begin the main algorithm
    for k = 1:m.obbt_aggressive_max_iteration

        # Set index differences and vector for filtering direction
        bool_indx_diff(m._lower_indx_diff, m._old_low_index, m._new_low_index)
        bool_indx_diff(m._upper_indx_diff, m._old_upp_index, m._new_upp_index)

        for i = 1:obbt_variable_count
            vi = @inbounds v[i]
            if @inbounds m._lower_indx_diff[i] && vi < 0.0
                @inbounds v[i] = 0.0
            end
            if @inbounds m._upper_indx_diff[i] && vi > 0.0
                @inbounds v[i] = 0.0
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
            for i = 1:obbt_variable_count
                vp_value =  @inbounds variable_primal[i]
                if @inbounds m._old_low_index[i] && vp_value == @inbounds n.lower_variable_bounds[i]
                    @inbounds m._new_low_index[i] = false
                end
                if @inbounds m._old_upp_index[i] && vp_value == @inbounds n.upper_variable_bounds[i]
                    @inbounds m._new_upp_index[i] = false
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
function set_node_flag!(m::Optimizer)
    for constr in m._working_problem._nonlinear_constr
        set_node_flag!(constr)
    end

    if m._working_problem._objective_type === NONLINEAR
        set_node_flag!(m._working_problem._objective_nl)
    end

    return nothing
end

"""
    set_reference_point!
"""
function set_reference_point!(m::Optimizer)

    evaluator = m._working_problem._relaxed_evaluator
    evaluator_x = evaluator.x
    current_xref = m._current_xref

    new_reference_point = false
    for node_i = 1:m._branch_variable_count
        solution_i = m._branch_to_sol_map[node_i]

        node_x = current_xref[node_i]
        solution_x = evaluator_x[solution_i]

        if node_x !== solution_x
            evaluator_x[solution_i] = node_x
            new_reference_point = true
        end
    end

    if new_reference_point

        for constr in m._working_problem._nonlinear_constr
            constr.has_value = false
        end

        if m._working_problem._objective_type === NONLINEAR
            m._working_problem._objective_nl.has_value = false
        end
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
function obbt!(m::Optimizer)

    feasibility = true

    n = m._current_node
    branch_to_sol_map = m._branch_to_sol_map
    relaxed_optimizer = m.relaxed_optimizer

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
    relax_objective!(m, 1)
    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(relaxed_optimizer)

    # Sets indices to attempt OBBT on
    obbt_variable_count = m._obbt_variable_count
    fill!(m._obbt_working_lower_index, true)
    fill!(m._obbt_working_upper_index, true)

    # Filters out any indicies with active bounds on variables
    # determined by solving the feasibility problem
    trivial_filtering!(m, n)

    # Applies an aggressive filter to eliminate indices that
    # cannot be tightened by obbt
    if m.obbt_aggressive_on
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
            for i = 1:obbt_variable_count
                if @inbounds m._obbt_working_lower_index[i]
                    sol_indx = branch_to_sol_map[i]
                    temp_value = @inbounds xLP[sol_indx] - n.lower_variable_bounds[i]
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
            for i = 1:obbt_variable_count
                if @inbounds m._obbt_working_upper_index[i]
                    sol_indx = branch_to_sol_map[i]
                    temp_value = @inbounds n.upper_variable_bounds[i] - xLP[sol_indx]
                    if temp_value <= upper_value
                        upper_value = temp_value
                        upper_indx = i
                    end
                end
            end
        end

        # default to upper bound if no lower bound is found, use maximum distance otherwise
        if lower_value <= upper_value && lower_indx > 0

            @inbounds m._obbt_working_lower_index[lower_indx] = false
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
                previous_value = n.lower_variable_bounds[lower_indx]

                # if bound is improved update node and corresponding constraint update
                # the node bounds and the single variable bound in the relaxation
                # we assume branching does not occur on fixed variables and interval
                # constraints are internally bridged by EAGO. So the only L <= x
                # constraint in the model is a GreaterThan.
                if updated_value > previous_value && (updated_value - previous_value) > 1E-6
                    ci_list = MOI.get(relaxed_optimizer, MOI.ListOfConstraintIndices{SAF,LT}())
                    sv_geq_ci = m._node_to_sv_geq_ci[lower_indx]
                    MOI.set(relaxed_optimizer, MOI.ConstraintSet(), sv_geq_ci, GT(updated_value))
                    @inbounds n.lower_variable_bounds[lower_indx] = updated_value
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
                previous_value = n.upper_variable_bounds[upper_indx]

                # if bound is improved update node and corresponding constraint update
                # the node bounds and the single variable bound in the relaxation
                # we assume branching does not occur on fixed variables and interval
                # constraints are internally bridged by EAGO. So the only U => x
                # constraint in the model is a LessThan.
                if updated_value < previous_value && (previous_value - updated_value) > 1E-6
                    sv_leq_ci = m._node_to_sv_leq_ci[upper_indx]
                    MOI.set(relaxed_optimizer, MOI.ConstraintSet(), sv_leq_ci, LT(updated_value))
                    @inbounds n.upper_variable_bounds[upper_indx] = updated_value
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
        trivial_filtering!(m, n)
    end

    return feasibility
end

"""
    load_fbbt_buffer!
"""
function load_fbbt_buffer!(m::Optimizer)

    n = m._current_node
    sol_to_branch = m._sol_to_branch_map
    lower_variable_bounds = n.lower_variable_bounds
    upper_variable_bounds = n.upper_variable_bounds

    for i = 1:m._working_problem._variable_count
        if @inbounds m._branch_variables[i]
            indx = @inbounds sol_to_branch[i]
            @inbounds m._lower_fbbt_buffer[i] = lower_variable_bounds[indx]
            @inbounds m._upper_fbbt_buffer[i] = upper_variable_bounds[indx]

        else
            @inbounds m._lower_fbbt_buffer[i] = m._working_problem._variable_info[i].lower_bound
            @inbounds m._upper_fbbt_buffer[i] = m._working_problem._variable_info[i].upper_bound

        end
    end

    return nothing
end

"""
    unpack_fbbt_buffer!
"""
function unpack_fbbt_buffer!(m::Optimizer)

    n = m._current_node
    sol_to_branch = m._sol_to_branch_map
    lower_variable_bounds = n.lower_variable_bounds
    upper_variable_bounds = n.upper_variable_bounds

    for i = 1:m._working_problem._variable_count
        if m._branch_variables[i]
            indx = sol_to_branch[i]
            if m._lower_fbbt_buffer[i] > lower_variable_bounds[indx]
                lower_variable_bounds[indx] = m._lower_fbbt_buffer[i]
            end
            if upper_variable_bounds[indx] > m._upper_fbbt_buffer[i]
                upper_variable_bounds[indx] = m._upper_fbbt_buffer[i]
            end
        end
    end

    return nothing
end

"""
    fbbt!

Performs feasibility-based bound tightening on a back-end constraint and returns `true` if it is feasible or
`false` if it is infeasible.
"""
function fbbt! end

function fbbt!(m::Optimizer, f::AffineFunctionIneq)

    # compute full sum
    lower_bounds = m._lower_fbbt_buffer
    upper_bounds = m._upper_fbbt_buffer

    terms = f.terms
    temp_sum = -f.constant

    for k = 1:f.len

        aik, indx_k = @inbounds terms[k]
        if aik !== 0.0
            aik_xL = aik*(@inbounds lower_bounds[indx_k])
            aik_xU = aik*(@inbounds upper_bounds[indx_k])
            temp_sum -= min(aik_xL, aik_xU)
        end
    end

    # subtract extra term, check to see if implied bound is better, if so update the node and
    # the working sum if the node is now empty then break

    for k = 1:f.len

        aik, indx_k = @inbounds terms[k]
        if aik !== 0.0

            xL = @inbounds lower_bounds[indx_k]
            xU = @inbounds upper_bounds[indx_k]

            aik_xL = aik*xL
            aik_xU = aik*xU

            temp_sum += min(aik_xL, aik_xU)
            xh = temp_sum/aik

            if aik > 0.0
                (xh < xL) && return false
                if xh > xL
                    @inbounds upper_bounds[indx_k] = xh
                end

            elseif aik < 0.0
                (xh > xU) && return false
                if xh < xU
                    @inbounds lower_bounds[indx_k] = xh
                end

            else
                temp_sum -= min(aik_xL, aik_xU)
                continue

            end

            aik_xL = aik*(@inbounds lower_bounds[indx_k])
            aik_xU = aik*(@inbounds upper_bounds[indx_k])
            temp_sum -= min(aik_xL, aik_xU)

        end
    end

    return true
end

function fbbt!(m::Optimizer, f::AffineFunctionEq)

    # compute full sum
    lower_bounds = m._lower_fbbt_buffer
    upper_bounds = m._upper_fbbt_buffer

    terms = f.terms
    temp_sum_leq = -f.constant
    temp_sum_geq = -f.constant

    for k = 1:f.len
        aik, indx_k = @inbounds terms[k]

        if aik !== 0.0
            aik_xL = aik*(@inbounds lower_bounds[indx_k])
            aik_xU = aik*(@inbounds upper_bounds[indx_k])
            temp_sum_leq -= min(aik_xL, aik_xU)
            temp_sum_geq -= max(aik_xL, aik_xU)

        end
    end

    # subtract extra term, check to see if implied bound is better, if so update the node and
    # the working sum if the node is now empty then break
    for k = 1:f.len

        aik, indx_k = @inbounds terms[k]
        if aik !== 0.0

            xL = @inbounds lower_bounds[indx_k]
            xU = @inbounds upper_bounds[indx_k]

            aik_xL = aik*xL
            aik_xU = aik*xU

            temp_sum_leq += min(aik_xL, aik_xU)
            temp_sum_geq += max(aik_xL, aik_xU)

            xh_leq = temp_sum_leq/aik
            xh_geq = temp_sum_geq/aik

            if aik > 0.0
                (xh_leq < xL) && return false
                if xh_leq > xL
                    @inbounds upper_bounds[indx_k] = xh_leq
                end
                (xh_geq > xU) && return false
                if (xh_geq < xU)
                    @inbounds lower_bounds[indx_k] = xh_geq
                end

            elseif aik < 0.0
                (xh_leq > xU) && return false
                if xh_leq < xU
                    @inbounds lower_bounds[indx_k] = xh_leq
                end
                (xh_geq < xL) && return false
                if (xh_geq > xL)
                    @inbounds upper_bounds[indx_k] = xh_geq
                end

            else
                temp_sum_leq -= min(aik_xL, aik_xU)
                temp_sum_geq -= max(aik_xL, aik_xU)
                continue

            end
            aik_xL = aik*(@inbounds lower_bounds[indx_k])
            aik_xU = aik*(@inbounds upper_bounds[indx_k])

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

        m._working_problem._relaxed_evaluator.is_first_eval = m._new_eval_objective
        if feasible_flag && (m._working_problem._objective_type === NONLINEAR)
            obj_nonlinear = m._working_problem._objective_nl
            set_node_flag!(obj_nonlinear)
            forward_pass!(evaluator, obj_nonlinear)
            feasible_flag &= reverse_pass!(evaluator, obj_nonlinear)
            evaluator.interval_intersect = true
            if feasible_flag
                set_node_flag!(obj_nonlinear)
                forward_pass!(evaluator, obj_nonlinear)
            end
        end

        m._new_eval_constraint = false
        m._new_eval_objective = false

        retrieve_x!(m._current_xref, evaluator)
        m._current_node = retrieve_node(evaluator)
    end

    return feasible_flag
end
