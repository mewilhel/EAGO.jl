include(joinpath(@__DIR__, "branch_variable.jl"))
include(joinpath(@__DIR__, "branch_special.jl"))

_branchvar_lower_solution(m,i) = m._lower_solution[m._branch_to_sol_map[i]]
function _fractional_integer_branch!(t::ExtensionType, m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
    if !m._current_node.cont
        integer_count = 0
        for i in m.branch_variables
            if _is_unfixed_integer(m, i)
                isinteger(_lower_solution(BranchVar, m, i))
                integer_count += 1
            end
        end
    end
    return false
end


#=
5.3.1 of Belotti, Pietro, et al. "Branching and bounds tighteningtechniques for
non-convex MINLP." Optimization Methods & Software 24.4-5 (2009): 597-634.
=#
function _select_branch_point!(m::GlobalOptimizer{N,T,S}, i) where {N,T<:Real,S}
    l = _lower_bound(BranchVar, m, i)
    u = _upper_bound(BranchVar, m, i)
    s = _lower_solution(BranchVar, m, i)
    α = _branch_cvx_α(m)
    b = _branch_offset_β(m)*(u - l)
    return max(l + b, min(u - b, α*s + (one(T) - α)*_mid(BranchVar, m, i)))
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
    push!(m._stack, NodeBB(_upper_bound(:branch, m),
                           ntuple(i -> (i == bi ? nextfloat(bp) : _upper_bound(BranchVar, m, i)), N),
                           lower_bound, upper_bound, new_depth, m._maximum_node_id,
                           bi, BD_NEG))

    m._maximum_node_id += 1
    push!(m._stack, NodeBB(ntuple(i -> (i == bi ? prevfloat(bp) : _lower_bound(BranchVar, m, i)), N),
                           _upper_bound(:branch, m),
                           lower_bound, upper_bound, new_depth, m._maximum_node_id,
                           bi, BD_POS))

    m._node_repetitions = 1
    m._node_count += 2
    return true
end

"""
$(SIGNATURES)

Creates two nodes from `current_node` using information available the `x`
and stores them to the stack. By default, relative width bisection is perfomed
at a point `branch_pnt` which is a convex combination
(parameter: `branch_cvx_factor`) of the solution to the relaxation and
the midpoint of the node. If this solution lies within `branch_offset/width` of
a bound then the branch point is moved to a distance of `branch_offset/width`
from the bound.
"""
function branch_node!(t::ExtensionType, m::GlobalOptimizer)
    was_branched = _conic_branch!(t, m)
    !was_branched && (was_branched &= _fractional_integer_branch!(t, m))
    !was_branched && (was_branched &= _continuous_branch!(t, m))
    return nothing
end
