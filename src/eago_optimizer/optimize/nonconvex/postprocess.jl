"""
    _variable_dbbt!

Tightens the bounds of the `_current_node` using the current global upper bound
and the duality information obtained from the relaxation.
"""
function _variable_dbbt!(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}

    LBD = m._lower_objective_value
    UBD = m._global_upper_bound
    delta = UBD - LBD
    if LBD <= UBD
        @inbounds for i = 1:m._branch_variable_num
            ml = m._lower_lvd[i]
            if ml > zero(T)
                cut = _lower_bound(BranchVar,m,i) + delta/ml
                if cut < _upper_bound(BranchVar,m,i)
                    _set_upper_bound!(BranchVar,m,i,cut)
                end
            else
                mh = m._lower_uvd[i]
                if mh > zero(T)
                    cut = _upper_bound(BranchVar,m,i) - delta/mh
                    if cut > _lower_bound(BranchVar,m,i)
                        _set_lower_bound!(BranchVar,m,i,cut)
                    end
                end
            end
        end
    end
    return nothing
end

"""
$(SIGNATURES)

Default postprocess perfoms duality-based bound tightening only.
"""
function postprocess!(t::ExtensionType, m::GlobalOptimizer)
    (m.dbbt_depth > m._iteration_count) && _variable_dbbt!(m)
    return nothing
end
