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
    adj_atol = d._parameters.absolute_tolerance/100.0
    adj_rtol = d._parameters.relative_tolerance/100.0
    if v > 0.0
        d._upper_objective_value = v*(1.0 + adj_rtol) + adj_atol
    else
        d._upper_objective_value = v*(1.0 - adj_rtol) + adj_atol
    end

    return nothing
end


revert_adjusted_upper_bound!(t::ExtensionType, d::Optimizer) = nothing

function revert_adjusted_upper_bound!(t::DefaultExt, d::Optimizer)

    adj_atol = d._parameters.absolute_tolerance/100.0
    adj_rtol = d._parameters.relative_tolerance/100.0

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

"""

Constructs and solves the problem locally on on node `y` updated the upper
solution informaton in the optimizer.
"""
function solve_local_nlp!(m::Optimizer)

    nlp_optimizer = m.nlp_optimizer
    MOI.empty!(nlp_optimizer)
    set_default_config!(nlp_optimizer)

    upper_variables = m._upper_variables
    for i = 1:m._working_problem._variable_count
        @inbounds upper_variables[i] = MOI.add_variable(nlp_optimizer)
    end

    n = m._current_node
    sol_to_branch_map = m._sol_to_branch_map
    lower_variable_bounds = n.lower_variable_bounds
    upper_variable_bounds = n.upper_variable_bounds
    variable_info = m._input_problem._variable_info

    lvb = 0.0
    uvb = 0.0
    x0 = 0.0

    for i = 1:m._input_problem._variable_count
        vinfo = @inbounds variable_info[i]
        single_variable = MOI.SingleVariable(@inbounds upper_variables[i])

        if vinfo.branch_on === BRANCH
            if vinfo.is_integer
            else
                indx = @inbounds sol_to_branch_map[i]
                lvb  = @inbounds lower_variable_bounds[indx]
                uvb  = @inbounds upper_variable_bounds[indx]
                if vinfo.is_fixed
                    MOI.add_constraint(nlp_optimizer, single_variable, ET(lvb))

                elseif vinfo.has_lower_bound
                    if vinfo.has_upper_bound
                        MOI.add_constraint(nlp_optimizer, single_variable, LT(uvb))
                        MOI.add_constraint(nlp_optimizer, single_variable, GT(lvb))

                    else
                        MOI.add_constraint(upper_optimizer, single_variable, GT(lvb))

                    end
                elseif vinfo.has_upper_bound
                    MOI.add_constraint(nlp_optimizer, single_variable, LT(uvb))

                end
            end
            x0 = 0.5*(lvb + uvb)
            upper_variable_index = @inbounds upper_variables[i]
            MOI.set(nlp_optimizer, MOI.VariablePrimalStart(), upper_variable_index, x0)

        else
            # not branch variable
            if vinfo.is_integer
            else
                lvb  = vinfo.lower_bound
                uvb  = vinfo.upper_bound
                if vinfo.is_fixed
                    MOI.add_constraint(nlp_optimizer, single_variable, ET(lvb))

                elseif vinfo.has_lower_bound
                    if vinfo.has_upper_bound
                        MOI.add_constraint(nlp_optimizer, single_variable, LT(uvb))
                        MOI.add_constraint(nlp_optimizer, single_variable, GT(lvb))

                    else
                        MOI.add_constraint(nlp_optimizer, single_variable, GT(lvb))

                    end
                elseif vinfo.has_upper_bound
                    MOI.add_constraint(nlp_optimizer, single_variable, LT(uvb))
                end
                x0 = 0.5*(lvb + uvb)
                upper_variable_index = @inbounds upper_variables[i]
                MOI.set(nlp_optimizer, MOI.VariablePrimalStart(), upper_variable_index, x0)
            end
        end
    end

    # Add linear, quadratic, and conic constraints constraints to model
    _add_linear_constraints!(m, nlp_optimizer)
    _add_quadratic_constraints!(m, nlp_optimizer)
    _add_soc_constraints!(m, nlp_optimizer)

    # Add nonlinear evaluation block
    MOI.set(nlp_optimizer, MOI.NLPBlock(), m._working_problem._nlp_data)
    MOI.set(nlp_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    # set objective as NECESSARY
    _add_sv_or_aff_obj!(m, nlp_optimizer)
    if m._input_problem._objective_type === SCALAR_QUADRATIC
        MOI.set(nlp_optimizer, MOI.ObjectiveFunction{SQF}(), m._input_problem._objective_sqf)
    end

    # Optimizes the object
    MOI.optimize!(nlp_optimizer)

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
