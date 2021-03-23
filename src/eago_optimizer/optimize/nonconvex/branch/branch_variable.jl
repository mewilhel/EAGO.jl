@enum(BRANCH_WEIGHT_STRATEGY, BW_INFEASIBLE, BW_INTERVAL_BRANCH,
                              BW_INTERVAL_BRANCH_REV, BW_INTERVAL_LP,
                              BW_INTERVAL_LP_REV)

Base.@kwdef mutable struct BranchOracle{T<:Real}
    strategy::BRANCH_WEIGHT_STRATEGY = BW_INTERVAL_LP
    𝛹m::Vector{T}     = T[]
    𝛹p::Vector{T}      = T[]
    δm::Vector{T}     = T[]
    δp::Vector{T}      = T[]
    ηm::Vector{T}     = T[]
    ηp::Vector{T}        = T[]
    score::Vector{T}         = T[]
    score_indx::Vector{Int}  = Int[]
    μ1::T                  = 0.1
    μ2::T                  = 1.3
    μ3::T                  = 0.8
    β::T                   = 0.05
    μ_score::T             = 0.15
    η_rel::T               = 3
end
function BranchOracle(T::Type{<:AbstractFloat}, n)
    BranchOracle{T}(𝛹m = ones(T,n),
                    𝛹p = ones(T,n),
                    δm = zeros(T,n),
                    δp = zeros(T,n),
                    ηm = zeros(T,n),
                    ηp = zeros(T,n)
                    )
end

function _variable_infeasibility(m, d::BranchOracle{T}, i) where T<:Real
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

function _set_δ!(::Val{BW_INFEASIBLE}, m::GlobalOptimizer, d::BranchOracle, i) where T<:Real
    _set_constraint_infeasibility!(m, i) # sets m._constr_inf
    v = _variable_infeasibility(m, d, i)
    m._variable_inf[i] = v
    d.δm[i] = v
    d.δp[i] = v
    return
end

function _set_δ!(::Val{BW_INTERVAL_BRANCH}, m::GlobalOptimizer, d::BranchOracle, i) where {T<:Real}
    l = _lower_bound(BranchVar, m, i)
    u = _upper_bound(BranchVar, m, i)
    d.δm[i] = isfinite(l) ? (xb - l) : _variable_infeasibility(m, d, i)
    d.δp[i]  = isfinite(u) ? (u - xb) : _variable_infeasibility(m, d, i)
    return
end

function _set_δ!(::Val{BW_INTERVAL_BRANCH_REV}, m::GlobalOptimizer, d::BranchOracle, i) where {T<:Real}
    l = _lower_bound(BranchVar, m, i)
    u = _upper_bound(BranchVar, m, i)
    d.δm[i] = isfinite(l) ? (u - xb) : _variable_infeasibility(m, d, i)
    d.δp[i] = isfinite(u) ? (xb - l) : _variable_infeasibility(m, d, i)
    return
end

function _set_δ!(::Val{BW_INTERVAL_LP}, m::GlobalOptimizer, d::BranchOracle, i) where {T<:Real}
    l = _lower_bound(BranchVar, m, i)
    u = _upper_bound(BranchVar, m, i)
    ρ = d.β*(u - l)
    xlp_adj = max(min(xlp, u - ρ), l + ρ)
    d.δm[i] = isfinite(l) ? (xlp_adj - l) : _variable_infeasibility(m, d, i)
    d.δp[i] = isfinite(u) ? (u - xlp_adj) : _variable_infeasibility(m, d, i)
    return
end

function _set_δ!(::Val{BW_INTERVAL_LP_REV}, m::GlobalOptimizer, d::BranchOracle, i) where {T<:Real}
    l = _lower_bound(BranchVar, m, i)
    u = _upper_bound(BranchVar, m, i)
    ρ = d.β*(u - l)
    xlp_adj = max(min(xlp, u - ρ), l + ρ)
    d.δm[i] = isfinite(l) ? (u - xlp_adj) : _variable_infeasibility(m, d, i)
    d.δp[i] = isfinite(u) ? (xlp_adj - l) : _variable_infeasibility(m, d, i)
    return
end

function _set_δ!(m::GlobalOptimizer, i)
    s = m._branch_oracle.strategy
    if s == BW_INTERVAL_LP
        _set_δ!(Val(BW_INTERVAL_LP), m, m._branch_oracle, i)
    elseif s == BW_INTERVAL_LP_REV
        _set_δ!(Val(BW_INTERVAL_LP_REV), m, m._branch_oracle, i)
    elseif s == BW_INTERVAL_BRANCH
        _set_δ!(Val(BW_INTERVAL_BRANCH), m, m._branch_oracle, i)
    elseif s == BW_INTERVAL_BRANCH_REV
        _set_δ!(Val(BW_INTERVAL_BRANCH_REV), m, m._branch_oracle, i)
    elseif s == BW_INFEASIBLE
        _set_δ!(Val(BW_INFEASIBLE), m, m._branch_oracle, i)
    end
    return nothing
end

function _pseudocost_branch!()
end

#=
Algorithm 3, Section 2.6 of Achterberg, Tobias, Thorsten Koch, and Alexander
Martin. "Branching rules revisited." Operations Research Letters 33.1 (2005):
42-54.
=#
function _reliability_branch!()
    MOI.set(m, MOI.InterationNumber(), k)
    _score!(m)
    sort_perm!(score_indx, score, rev = true)

    max_score = maximum(m._branch_score)
    no_improve_count = 0
    for i in unreliable
        if no_improve_count > no_improve_limit
            break
        end

        # solve subproblem lower
        obj_minus = MOI.get(m, MOI.ObjectiveValue()) - obj_val
        # solve subproblem upper
        obj_plus = MOI.get(m, MOI.ObjectiveValue()) - obj_val

        _set_pseudo_cost!(obj_plus, obj_minus, i)

        new_score = _score(obj_plus, obj_minus, d.μ_score)
        if new_score >= max_score
            no_improve_count += 1
        else
            m._branch_score[i] = new_score
            no_improve_count = 0
        end
    end
    MOI.set(m, MOI.InterationNumber(), k_orig)
end


_score(x::T, y, μ) where T<:Real = (one(T) - μ)*min(x, y) + max(x, y)
_score(d::BranchOracle{T}, i) where T<:Real = _score(𝛹m[i]*δm[i], 𝛹p[i]*δp[i], d.μ_score)
function _score!(m)
    map(i -> _score(m._branch_oracle, i), m._branch_score, m._branch_variable)
    return
end

function _select_branch_variable!(m)
    _score!(m)
    return argmax(m._branch_score)
end
