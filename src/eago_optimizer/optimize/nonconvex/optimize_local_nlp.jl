# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimize/optimize_convex.jl
# Contains the solve_local_nlp! routine which computes the optimal value
# of a convex function. This is used to compute the upper bound in the
# branch and bound routine. A number of utility function required for
# solve_local_nlp! are also included.
#############################################################################


function _update_branch_variables!(nlp_opt, m)
    n = m._current_node
    for i = 1:m._working_problem._variable_num
        vinfo = @inbounds m._working_problem._variable_info[i]
        if m._branch_variables[i]
            if vinfo.is_integer
            else
                indx = @inbounds m._sol_to_branch_map[i]
                lvb  = @inbounds n.lower_variable_bound[indx]
                uvb  = @inbounds n.upper_variable_bound[indx]
                if vinfo.is_fixed
                    MOI.add_constraint(nlp_opt, SV(m.upper_variables[i]), ET(lvb))
                elseif vinfo.has_lower_bound
                    if vinfo.has_upper_bound
                        MOI.add_constraint(nlp_opt, SV(m.upper_variables[i]), LT(uvb))
                        MOI.add_constraint(nlp_opt, SV(m.upper_variables[i]), GT(lvb))
                    else
                        MOI.add_constraint(nlp_opt, SV(m.upper_variables[i]), GT(lvb))
                    end
                elseif vinfo.has_upper_bound
                    MOI.add_constraint(nlp_opt, SV(m.upper_variables[i]), LT(uvb))
                end
            end
        end
    end
    return nothing
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

    (t == MOI.OPTIMAL) && (termination_flag = true)
    (t == MOI.LOCALLY_SOLVED) && (termination_flag = true)

    # This is default solver specific... the acceptable constraint tolerances
    # are set to the same values as the basic tolerance. As a result, an
    # acceptably solved solution is feasible but non necessarily optimal
    # so it should be treated as a feasible point
    if (t == MOI.ALMOST_LOCALLY_SOLVED) && (r == MOI.NEARLY_FEASIBLE_POINT)
        termination_flag = true
        result_flag = true
    end

    (r == MOI.FEASIBLE_POINT) && (result_flag = true)

    return (termination_flag && result_flag)
end

# TODO: Select between infeasible starting point reduction strategies
function _set_starting_point!(opt, m)
    for i = 1:m._working_problem._variable_num
        if m._branch_variables[i]
            vp = mid(m, m._sol_to_branch_map[i])
        else
            vp = mid(m._working_problem.variable_info[i])
            if isinf(vp)
                if isinf(_lower_bound(m._working_problem.variable_info[i]))
                    vp = _lower_bound(m._working_problem.variable_info[i])
                else
                    vp = _upper_bound(m._working_problem.variable_info[i])
                end
            end
        end
        MOI.set(opt, MOI.VariablePrimalStart(m._upper_variables[i]), vp)
    end
    return nothing
end

function _unpack_local_nlp_solve!(m::Optimizer, opt::T, idx_map; adjust_bnd::Bool = true) where T
    # Process output info and save to CurrentUpperInfo object
    m._upper_termination_status = MOI.get(nlp_optimizer, MOI.TerminationStatus())
    m._upper_result_status = MOI.get(nlp_optimizer, MOI.PrimalStatus())
    if is_feasible_solution(m._upper_termination_status, m._upper_result_status)
        m._upper_feasibility = true
        value = MOI.get(nlp_optimizer, MOI.ObjectiveValue())
        m._best_upper_value = min(value, m._best_upper_value)
        m._upper_solution .= MOI.get(nlp_optimizer, MOI.VariablePrimal(), upper_variables)
    else
        m._upper_feasibility = false
        m._upper_objective_value = Inf
    end
    return nothing
end

"""

Constructs and solves the problem locally on on node `y` updated the upper
solution informaton in the optimizer.
"""
function solve_local_nlp!(m::Optimizer)

    opt = m.nlp_optimizer
    #set_config!(m, opt)
    bridged_opt = _bridge_optimizer(Val{DIFF_CVX}(), opt)
    idx_map = MOIU.default_copy_to(bridged_opt, m._input_problem, false)
    (m.verbosity < 5) && MOI.set(bridged_opt, MOI.Silent(), true)

    _update_branch_variables!(bridged_opt, m)
    _set_starting_point!(bridged_opt, m)
    MOI.optimize!(bridged_opt)
    _unpack_local_nlp_solve!(m, bridged_opt, idx_map)

    return nothing
end

"""
$(SIGNATURES)

Default check to see if the upper bounding problem should be run. By default,
The upper bounding problem is run on every node up to depth `upper_bounding_depth`
and is triggered with a probability of `0.5^(depth - upper_bounding_depth)`
afterwards.
"""
function default_nlp_heurestic(m::GlobalOptimizer{N,T})
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
function upper_problem!(t::ExtensionType, m::GlobalOptimizer{N,T,S}) where {N,T<:Real,S}

    if !default_nlp_heurestic(m)
        m._upper_feasibility = false
        m._upper_objective_value = Inf
    else
        solve_local_nlp!(m)
    end
    return nothing
end
