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

function _cvx_relative_width_branch!(t::ExtensionType, m::GlobalOptimizer{N,T}) where {N,T}
    n = m._current_node

    lvbs = m._lower_variable_bound
    uvbs = m._upper_variable_bound

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
    lower_bound = max(n.lower_bound, m._lower_objective_value)
    upper_bound = min(n.upper_bound, m._upper_objective_value)
    new_depth = n.depth + 1

    m._maximum_node_id += 1
    uvb_1 = ntuple(i -> (i == max_pos ? nextfloat(branch_pnt) : _lower_branch_bound(m, i)), N)
    X1 = NodeBB(_lower_branch_bound(m), uvb_1, lower_bound, upper_bound, new_depth, m._maximum_node_id)
    push!(m._stack, X1)

    m._maximum_node_id += 1
    lvb_2 = ntuple(i -> (i == max_pos ? prevfloat(branch_pnt) : _lower_branch_bound(m, i)), N)
    X2 = NodeBB(lvb_2, _upper_branch_bound(m), lower_bound, upper_bound, new_depth, m._maximum_node_id)
    push!(m._stack, X2)

    m._node_repetitions = 1
    m._node_count += 2
    return true
end
