"""
$(SIGNATURES)

Default postprocess perfoms duality-based bound tightening on the `y`.
"""
function postprocess!(t::ExtensionType, m::GlobalOptimizer{N,T}) where {N,T<:Real}

    if m.dbbt_depth > m._iteration_count
        variable_dbbt!(m._current_node, m._lower_lvd, m._lower_uvd,
                       m._lower_objective_value, m._global_upper_bound,
                       m._branch_variable_num)
    end

    return nothing
end
