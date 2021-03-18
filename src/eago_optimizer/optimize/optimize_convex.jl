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

"""

Shifts the resulting local nlp objective value `f*` by `(1.0 + relative_tolerance/100.0)*f* + absolute_tolerance/100.0`.
This assumes that the local solvers relative tolerance and absolute tolerance is significantly lower than the global
tolerance (local problem is minimum).
"""
function stored_adjusted_upper_bound!(d::Optimizer, v::Float64)
    adj_atol = d.absolute_tolerance/100.0
    adj_rtol = d.relative_tolerance/100.0
    if v > 0.0
        d._upper_objective_value = v*(1.0 + adj_rtol) + adj_atol
    else
        d._upper_objective_value = v*(1.0 - adj_rtol) + adj_atol
    end

    return nothing
end

revert_adjusted_upper_bound!(t::ExtensionType, d::Optimizer) = nothing

function revert_adjusted_upper_bound!(t::DefaultExt, d::Optimizer)

    adj_atol = d.absolute_tolerance/100.0
    adj_rtol = d.relative_tolerance/100.0

    adj_objective_value = d._global_upper_bound
    adj_objective_value -= adj_atol
    if adj_objective_value > 0.0
        adj_objective_value /= (1.0 + adj_rtol)
    else
        adj_objective_value /= (1.0 - adj_rtol)
    end
    d._global_upper_bound = adj_objective_value

    return nothing
end

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
        stored_adjusted_upper_bound!(m, value)
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
