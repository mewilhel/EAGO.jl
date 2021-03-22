abstract type AbstractPseudoCostStrategy end
struct PCInfeasible <: AbstractPseudoCostStrategy end
struct PCIntervalBranch <: AbstractPseudoCostStrategy end
struct PCIntervalBranchRev <: AbstractPseudoCostStrategy end
struct PCIntervalLP <: AbstractPseudoCostStrategy end
struct PCIntervalLPRev <: AbstractPseudoCostStrategy end

mutable struct PseudoCost{S<:PseudoCostStrategy,T<:Real}
    strategy::S         = S()
    𝛹_minus::Vector{T}  = T[]
    𝛹_plus::Vector{T}   = T[]
    δ_minus::Vector{T}  = T[]
    δ_plus::Vector{T}   = T[]
    μ1::T               = 0.1
    μ2::T               = 1.3
    μ3::T               = 0.8
    β::T                = 0.05
    μ_score::T          = 0.15
end

_score(x, y, μ) = (1 - μ)*min(x, y) + max(x, y)
function _pseudocost(d::PseudoCost, i)
    _score(𝛹_minus[i]*δ_minus[i], 𝛹_plus[i]*δ_plus[i], d.μ_score)
end

function _variable_infeasibility(m, d::PseudoCost{S,T}, i) where {S,T}
    tsum = zero(T); tmin = Inf; tmax = -Inf
    for j in m._constr_dep[i]
        constraint_infeasiblity = m._constr_inf[j]
        tsum += constraint_infeasiblity
        if constraint_infeasiblity > tmax
            tmax = constraint_infeasiblity
        end
        if constraint_infeasiblity < tmin
            tmin = constraint_infeasiblity
        end
    end
    return d.μ1*tsum + d.μ2*tmin + d.μ3*tmax
end

function _set_δ!(m, d::PseudoCost{PCInfeasible,T}, i) where {T<:Real}
    _set_constraint_infeasibility!(m, i) # sets m._constr_inf
    v = _variable_infeasibility(m, d, i)
    m._variable_inf[i] = v
    d.δ_minus[i] = v
    d.δ_plus[i]  = v
    return
end

function _set_δ!(d::PseudoCost{PCBranchInterval,T}, xb, i) where {T<:Real}
    l = _lower_branch_bound(m, i)
    u = _upper_branch_bound(m, i)
    d.δ_minus[i] = isfinite(l) ? (xb - l) : _variable_infeasibility(m, d, i)
    d.δ_plus[i]  = isfinite(u) ? (u - xb) : _variable_infeasibility(m, d, i)
    return
end

function _set_δ!(d::PseudoCost{PCBranchIntervalRev,T}, xb, i) where {T<:Real}
    l = _lower_branch_bound(m, i)
    u = _upper_branch_bound(m, i)
    d.δ_minus[i] = isfinite(l) ? (u - xb) : _variable_infeasibility(m, d, i)
    d.δ_plus[i]  = isfinite(u) ? (xb - l) : _variable_infeasibility(m, d, i)
    return
end

function _set_δ!(d::PseudoCost{PCIntervalLP,T}, xlp, i) where {T<:Real}
    l = _lower_branch_bound(m, i)
    u = _upper_branch_bound(m, i)
    ρ = d.β*(u - l)
    xlp_adj = max(min(xlp, u - ρ), l + ρ)
    d.δ_minus[i] = isfinite(l) ? (xlp_adj - l) : _variable_infeasibility(m, d, i)
    d.δ_plus[i]  = isfinite(u) ? (u - xlp_adj) : _variable_infeasibility(m, d, i)
    return
end

function _set_δ!(d::PseudoCost{PCIntervalLPRev,T}) where {T<:Real}
    l = _lower_branch_bound(m, i)
    u = _upper_branch_bound(m, i)
    ρ = d.β*(u - l)
    xlp_adj = max(min(xlp, u - ρ), l + ρ)
    d.δ_minus[i] = isfinite(l) ? (u - xlp_adj) : _variable_infeasibility(m, d, i)
    d.δ_plus[i]  = isfinite(u) ? (xlp_adj - l) : _variable_infeasibility(m, d, i)
end

function _pseudocost_branch!()
end
