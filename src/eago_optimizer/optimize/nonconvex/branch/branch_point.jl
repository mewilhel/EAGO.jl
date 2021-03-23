function _select_branch_point!(m::GlobalOptimizer{N,T,S}, i) where {N,T<:Real,S}
    l = _branch_lo(m, i)
    u = _branch_hi(m, i)
    s = _branch_lower_solution(m, i)
    α = _branch_cvx_α(m)
    b = _branch_offset_β(m)*(u - l)
    return max(l + b, min(u - b, α*s + (one(T) - α)*_branch_mid(m, i))
end
