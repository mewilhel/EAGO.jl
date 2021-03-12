


if !m._obbt_performed_flag
    if m._nonlinear_evaluator_created
        set_node!(m._working_problem._relaxed_evaluator, m._current_node)
        set_node_flag!(m)
        set_reference_point!(m)
        fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)
    end
    update_relaxed_problem_box!(m)
    relax_constraints!(m)
end
relax_objective!(m)

# Optimizes the object
MOI.set(m.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

"""

"""
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

function _modify_existing_relax(m)
    flag = MOIU.supports_default_copy_to(m)
    flag &=
    return flag
end

function MOIU.supports_default_copy_to(opt::T)
end

function lower_problem!(t::ExtensionType, m::Optimizer)
    if _modify_existing_relax(m)
        _modify_relax!(m)
    elseif _add_delete_existing_relax(m)
        _add_delete_relax!(m)
    else
        MOI.empty!(opt)
        if MOIU.supports_default_copy_to(m.relaxed_optimizer)
            _add_relax!(m, m.relaxed_optimizer)
        else
            MOI.empty!(m.cached_opt)
            _add_relax!(m, cached_opt)
            MOIU.copyto!(m.relaxed_optimizer, m.cached_opt)
        end
    end

    MOI.optimize!(m.relaxed_optimizer)
    if !intrepret_relaxed_solution(m, m.relaxed_optimizer)
        fallback_interval_lower_bound!(m, n)
    end
    return
end




if supports_constraint_mod(opt)
    # modify
else
    MOI.empty!(opt)
    if MOIU.supports_default_copy_to(opt::T)
        _relax!(opt)
    else
        _relax!(cached_opt)
        MOIU.copyto!(opt, cached_opt)
    end
end
