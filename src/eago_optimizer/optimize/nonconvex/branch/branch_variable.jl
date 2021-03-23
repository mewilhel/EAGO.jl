@enum(BRANCH_WEIGHT_STRATEGY, BW_INFEASIBLE, BW_INTERVAL_BRANCH,
                              BW_INTERVAL_BRANCH_REV, BW_INTERVAL_LP,
                              BW_INTERVAL_LP_REV)

Base.@kwdef mutable struct BranchOracle{T<:Real}
    strategy::BRANCH_WEIGHT_STRATEGY = BW_INTERVAL_LP
    𝛹_minus::Vector{T} = T[]
    𝛹_plus::Vector{T}  = T[]
    δ_minus::Vector{T} = T[]
    δ_plus::Vector{T}  = T[]
    η_minus::Vector{T} = T[]
    η_plus::Vector{T}  = T[]
    μ1::T               = 0.1
    μ2::T               = 1.3
    μ3::T               = 0.8
    β::T                = 0.05
    μ_score::T          = 0.15
    η_rel::T            = 0.01 # TODO: FIND BASE
end
function BranchOracle(T::Type{<:AbstractFloat}, n)
    BranchOracle{T}(𝛹_minus = ones(T,n),
                    𝛹_plus = ones(T,n),
                    δ_minus = zeros(T,n),
                    δ_plus = zeros(T,n),
                    η_minus = zeros(T,n),
                    η_plus = zeros(T,n)
                    )
end

_score(x::T, y, μ) where T<:Real = (one(T) - μ)*min(x, y) + max(x, y)
function _score(d::BranchOracle{T}, i) where T<:Real
    _score(𝛹_minus[i]*δ_minus[i], 𝛹_plus[i]*δ_plus[i], d.μ_score)
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

function _set_δ!(m, d::BranchOracle{T}, i) where T<:Real
    _set_constraint_infeasibility!(m, i) # sets m._constr_inf
    v = _variable_infeasibility(m, d, i)
    m._variable_inf[i] = v
    d.δ_minus[i] = v
    d.δ_plus[i]  = v
    return
end

function _set_δ!(d::BranchOracle{T}, xb, i) where {T<:Real}
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

function _reliability_branch!()
    MOI.set(m, MOI.InterationNumber(), k)
    _pseudocost_branch!()
    MOI.set(m, MOI.InterationNumber(), k_orig)
end


function _select_branch_variable!(m)
    max_pos = 0
    max_val = -Inf
    temp_max = 0.0
    for i = 1:m._branch_variable_num
        si = m._branch_to_sol_map[i]
        vi = m._working_problem._variable_info[si]
        if m._branch_variables[i]
            temp_max = _score(m, i)
            #temp_max =  _current_width(m, i)    # uvbs[i] - lvbs[i]
            #temp_max /= _initial_width(m, i)    # vi.upper_bound - vi.lower_bound
            if temp_max > max_val
                max_pos = i
                max_val = temp_max
            end
        end
    end
    return max_pos
end
