# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/bound.jl
# Computes interval bounds of various functions.
#############################################################################

@inline function _get_lo(m::GlobalOptimizer{N,T,S}, t::SAT) where {N,T<:AbstractFloat,S}
    c = t.coefficient
    i = t.variable_index.value
    xL = _lower_bound(FullVar, m, i)
    xU = _upper_bound(FullVar, m, i)
    return (c > 0.0) ? c*xL : c*xU
end

@inline function _get_bnds(m::GlobalOptimizer{N,T,S}, t::SAT) where {N,T<:AbstractFloat,S}
    c = t.coefficient
    i = t.variable_index.value
    xL = _lower_bound(FullVar, m, i)
    xU = _upper_bound(FullVar, m, i)
    l = (c > 0.0) ? c*xL : c*xU
    u = (c > 0.0) ? c*xU : c*xL
    return l, u
end

@inline function _get_bnds(m::GlobalOptimizer{N,T,S}, t::SQT) where {N,T<:AbstractFloat,S}
    c = t.coefficient
    v1 = t.variable_index_1.value
    v2 = t.variable_index_2.value
    xL1 = _lower_bound(FullVar, m, v1)
    xU1 = _upper_bound(FullVar, m, v1)

    if v1 == v2
        if c > zero(T)
            bnds = (zero(T) < xL1) ? 0.5*c*xL1*xL1 : ((xU1 <= zero(T)) ? 0.5*c*xU1*xU1 : zero(T))
        else
            bnds = (xL1 < xU1) ? 0.5*c*xU1*xU1 : 0.5*c*xL1*xL1
        end
    else
        xL2 = _lower_bound(FullVar, m, v2)
        xU2 = _upper_bound(FullVar, m, v2)
        bnds = c*Interval{Float64}(xL, xU)*Interval{Float64}(xL2, xU2)
    end
    return bnds
end

function lower_interval_bound(m::GlobalOptimizer{N,T,S}, t::Vector{SAT}, c) where {N,T<:AbstractFloat,S}
    mapreduce(i -> _get_lo(m, i), +, t; init = c)
end
function interval_bound(m::GlobalOptimizer{N,T,S}, t::Vector{SQT}, c) where {N,T<:AbstractFloat,S}
    mapreduce(i -> _get_bnds(m, i), +, t; init = c)
end

function interval_bound(m::GlobalOptimizer{N,T,S}, f::AffineFunctionEq) where {N,T<:AbstractFloat,S}
    mapreduce(i -> _get_bnds(m, i), +, f.func.terms; init = (f.func.constant,
                                                             f.func.constant))
end

function lower_interval_bound(m::GlobalOptimizer{N,T,S}, f::AffineFunctionIneq) where {N,T<:AbstractFloat,S}
    lower_interval_bound(m, f.func.terms, f.func.constant)
end

function interval_bound(m::GlobalOptimizer{N,T,S}, f::BufferedQuadraticEq) where {N,T<:AbstractFloat,S}
    bnds = lower_interval_bound(m, f.func.affine_terms, Interval{Float64}(f.func.constant))
    interval_bound(m, f.func.quadratic_terms, lower)
end
function lower_interval_bound(m::GlobalOptimizer{N,T,S}, f::BufferedQuadraticIneq) where {N,T<:AbstractFloat,S}
    bnds = interval_bound(m, f.func.quadratic_terms, Interval{Float64}(f.func.constant))
    return bnds.lo
end

###
### SECOND-ORDER CONE
###
function lower_interval_bound(m::GlobalOptimizer{N,T,S}, d::BufferedSOC) where {N,T<:AbstractFloat,S}
    vec_of_vi = d.variables.variables
    norm_bound = zero(Interval{T})
    for i = 2:length(vec_of_vi)
        vi = vec_of_vi[i].value
        x = Interval{T}(_lower_bound(FullVar, m, vi), _upper_bound(FullVar, m, vi))
        norm_bound += pow(x, 2)
    end
    norm_bound = sqrt(norm_bound)
    return norm_bound.lo - _upper_bound(FullVar, m, vec_of_vi[1].value)
end


###
### NONLINEAR FUNCTIONS
###
function lower_interval_bound(m::GlobalOptimizer{N,T,S}, d::BufferedNonlinearFunction{V}) where {N,T<:AbstractFloat,S,V}
    !d.has_value && forward_pass!(m._working_problem._relaxed_evaluator, d)
    ex = d.expr
    v = ex.isnumber[1] ? ex.numberstorage[1] : ex.setstorage[1].Intv.lo
    return v
end
function interval_bound(m::GlobalOptimizer{N,T,S}, d::BufferedNonlinearFunction{V}) where {N,T<:AbstractFloat,S,V}
    !d.has_value && forward_pass!(d.evaluator, d)
    ex = d.expr
    v = ex.isnumber[1] ? Interval(ex.numberstorage[1]) : ex.setstorage[1].Intv
    return v.lo, v.hi
end

_is_feas(m::GlobalOptimizer{N,T,S}, x::AffineFunctionIneq) where {N,T<:AbstractFloat,S} = lower_interval_bound(m, x) <= 0.0
_is_feas(m::GlobalOptimizer{N,T,S}, x::BufferedQuadraticIneq) where {N,T<:AbstractFloat,S} = lower_interval_bound(m, x) <= 0.0
function _is_feas(m::GlobalOptimizer{N,T,S}, x::AffineFunctionEq) where {N,T<:AbstractFloat,S}
    lower_value, upper_value = interval_bound(m, x)
    return lower_value <= 0.0 <= upper_value
end
function _is_feas(m::GlobalOptimizer{N,T,S}, x::BufferedQuadraticEq) where {N,T<:AbstractFloat,S}
    lower_value, upper_value = interval_bound(m, x)
    return lower_value <= 0.0 <= upper_value
end
function _is_feas(m::GlobalOptimizer{N,T,S}, x::NonlinearExpression) where {N,T<:AbstractFloat,S}
    lower_value, upper_value = interval_bound(m, x)
    feasible_flag &= upper_value < x.lower_bound
    feasible_flag &= lower_value > x.upper_bound
    return !feasible_flag
end
