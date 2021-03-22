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
function branch_node!(t::ExtensionType, m::GlobalOptimizer{N,T})
    was_branched = _conic_branch!(t, m)
    !was_branched && (was_branched &= _fractional_integer_branch!(t, m))
    !was_branched && (was_branched &= _cvx_relative_width_branch!(t, m))
    return nothing
end
