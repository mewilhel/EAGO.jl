function _reset_relaxation!(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}

    m._working_problem._relaxed_evaluator.is_first_eval = true
    fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)

    m._new_eval_constraint = true

    delete_nl_constraints!(m)
    delete_objective_cuts!(m)

    return nothing
end

"""
$(SIGNATURES)

Runs interval, linear, quadratic contractor methods followed by obbt and a
constraint programming walk up to tolerances specified in
`EAGO.GlobalOptimizer` object.
"""
function preprocess!(t::ExtensionType, m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}

    _reset_relaxation!(m)

    wp = m._working_problem

    # Sets initial feasibility
    feasible_flag = true
    m._obbt_performed_flag = false

    # compute initial volume
    m._initial_volume = prod(i -> _mid(BranchVar, m, i), 1:N)

    if m.fbbt_lp_depth >= m._iteration_count
        for i = 1:m.fbbt_lp_repetitions
            if feasible_flag
                for saf_leq in wp._saf_leq
                    !feasible_flag && break
                    feasible_flag &= fbbt!(m, saf_leq)
                end
                !feasible_flag && break

                for saf_eq in  wp._saf_eq
                    !feasible_flag && break
                    feasible_flag &= fbbt!(m, saf_eq)
                end
                !feasible_flag && break
            end
        end
    end

    # done after cp to prevent using cp specific flags in cut generation
    set_first_relax_point!(m)

    cp_walk_count = 0
    perform_cp_walk_flag = feasible_flag
    perform_cp_walk_flag &= (m.cp_depth >= m._iteration_count)
    perform_cp_walk_flag &= (cp_walk_count < m.cp_repetitions)
    while perform_cp_walk_flag
        feasible_flag &= set_constraint_propagation_fbbt!(m)
        !feasible_flag && break
        cp_walk_count += 1
        perform_cp_walk_flag = (cp_walk_count < m.cp_repetitions)
    end

    obbt_count = 0
    perform_obbt_flag = feasible_flag
    perform_obbt_flag &= (m.obbt_depth >= m._iteration_count)
    perform_obbt_flag &= (obbt_count < m.obbt_repetitions)

    while perform_obbt_flag
        feasible_flag &= obbt!(m)
        m._obbt_performed_flag = true
        !feasible_flag && break
        obbt_count += 1
        perform_obbt_flag = (obbt_count < m.obbt_repetitions)
    end

    m._final_volume = prod(i -> _mid(BranchVar, m, i), 1:N)
    m._preprocess_feasibility = feasible_flag

    return nothing
end
