function _conic_branch!(t::ExtensionType, m::GlobalOptimizer{N,T}) where {N,T}
    return false
end

_branchvar_lower_solution(m,i) = m._lower_solution[m._branch_to_sol_map[i]]
function _fractional_integer_branch!(t::ExtensionType, m::GlobalOptimizer{N,T}) where {N,T}
    if !m._current_node.cont
        integer_count = 0
        for i in 1:m._branch_variable_num
            if _is_branchvar_unfixed_integer(m, i)
                isinteger(_branchvar_lower_solution(m, i))
                integer_count += 1
            end
        end
    end
    return false
end

function _select_branch_variable!(m)
    max_pos = 0
    max_val = -Inf
    temp_max = 0.0
    flag = true
    for i = 1:m._branch_variable_num
        si = m._branch_to_sol_map[i]
        vi = m._working_problem._variable_info[si]
        if m._branch_variables[i]
            temp_max =  _current_width(m, i)    # uvbs[i] - lvbs[i]
            temp_max /= _initial_width(m, i)    # vi.upper_bound - vi.lower_bound
            if temp_max > max_val
                max_pos = i
                max_val = temp_max
            end
        end
    end
    return max_pos
end

function _select_branch_point!(m)
    n = m._current_node
    lvbs = m._lower_variable_bound
    uvbs = m._upper_variable_bound
    lvb  = lvbs[bi]
    uvb  = uvbs[bi]
    si   = m._branch_to_sol_map[bi]
    lsol = m._lower_solution[si]
    cvx_f = m.branch_cvx_factor
    cvx_g = m.branch_offset
    bp = cvx_f*lsol + (1.0 - cvx_f)*(lvb + uvb)/2.0
    if bp < lvb*(1.0 - cvx_g) + cvx_g*uvb
        bp = (1.0 - cvx_g)*lvb + cvx_g*uvb
    elseif bp > cvx_g*lvb + (1.0 - cvx_g)*uvb
        bp = cvx_g*lvb + (1.0 - cvx_g)*uvb
    end
    return bp
end

function _continuous_branch!(t::ExtensionType, m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}

    bi = _select_branch_variable!(m)
    bp = _select_branch_point!(m)

    # rounds into branch points, which in turn prevents the
    # solution at the branch point from being discarded
    lower_bound = max(n.lower_bound, m._lower_objective_value)
    upper_bound = min(n.upper_bound, m._upper_objective_value)
    new_depth = n.depth + 1

    m._maximum_node_id += 1
    push!(m._stack, NodeBB(_branch_lo(m),
                           ntuple(i -> (i == bi ? nextfloat(bp) : _branch_hi(m, i)), N),
                           lower_bound, upper_bound, new_depth, m._maximum_node_id))

    m._maximum_node_id += 1
    push!(m._stack, NodeBB(ntuple(i -> (i == bi ? prevfloat(bp) : _branch_lo(m, i)), N),
                           _branch_hi(m),
                           lower_bound, upper_bound, new_depth, m._maximum_node_id))

    m._node_repetitions = 1
    m._node_count += 2
    return true
end
