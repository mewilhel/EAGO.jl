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

@inline function _get_lo(m, n, t::SAT)
    coeff = t.coefficient
    i = t.variable_index.value
    if m._branch_variables[i]
        mapped_vi = m._sol_to_branch_map[i]
        xL = n.lower_variable_bound[mapped_vi]
        xU = n.upper_variable_bound[mapped_vi]
    else
        xL = m._working_problem._variable_info[i].lower_bound
        xU = m._working_problem._variable_info[i].upper_bound
    end
    return (coeff > 0.0) ? coeff*xL : coeff*xU
end

@inline function _get_bnds(m, n, t::SAT)
    coeff = t.coefficient
    i = t.variable_index.value
    if m._branch_variables[i]
        mapped_vi = m._sol_to_branch_map[i]
        xL = n.lower_variable_bound[mapped_vi]
        xU = n.upper_variable_bound[mapped_vi]
    else
        xL = m._working_problem._variable_info[i].lower_bound
        xU = m._working_problem._variable_info[i].upper_bound
    end
    lower = (coeff > 0.0) ? coeff*xL : coeff*xU
    upper = (coeff > 0.0) ? coeff*xU : coeff*xL
    return lower, upper
end

@inline function _get_bnds(m, n, t::SQT)
    coeff = t.coefficient
    vi1 = t.variable_index_1.value
    vi2 = t.variable_index_2.value

    mapped_vi1 = m._sol_to_branch_map[vi1]
    xL = n.lower_variable_bound[mapped_vi1]
    xU = n.upper_variable_bound[mapped_vi1]

    if vi1 == vi2
        if coeff > 0.0
            bnds = (0.0 < xL) ? 0.5*coeff*xL*xL : ((xU <= 0.0) ? 0.5*coeff*xU*xU : 0.0)
        else
            bnds = (xL < xU) ? 0.5*coeff*xU*xU : 0.5*coeff*xL*xL
        end
    else
        mapped_vi2 = m._sol_to_branch_map[vi2]
        il2b = n.lower_variable_bound[mapped_vi2]
        iu2b = n.upper_variable_bound[mapped_vi2]
        bnds = coeff*Interval{Float64}(xL, xU)*Interval{Float64}(il2b, iu2b)
    end
    return bnds
end

function lower_interval_bound(m::Optimizer, n::NodeBB, t::Vector{SAT}, c)
    mapreduce(i -> _get_lo(m, n, i), +, t; init = c)
end
function interval_bound(m::Optimizer, n::NodeBB, t::Vector{SQT}, c)
    mapreduce(i -> _get_bnds(m, n, i), +, t; init = c)
end

function interval_bound(m::Optimizer, f::AffineFunctionEq, y::NodeBB)
    mapreduce(i -> _get_bnds(m, n, i), +, f.func.terms; init = (f.func.constant,
                                                                f.func.constant))
end

function lower_interval_bound(m::Optimizer, f::AffineFunctionIneq, n::NodeBB)
    lower_interval_bound(m, n, f.func.terms, f.func.constant)
end

function interval_bound(m::Optimizer, f::BufferedQuadraticEq, n::NodeBB)
    bnds = lower_interval_bound(m, n, f.func.affine_terms, Interval{Float64}(f.func.constant))
    interval_bound(m, n, f.func.quadratic_terms, lower)
end
function lower_interval_bound(m::Optimizer, f::BufferedQuadraticIneq, n::NodeBB)
    bnds = interval_bound(m, n, f.func.quadratic_terms, Interval{Float64}(f.func.constant))
    bnds.lo
end

###
### SECOND-ORDER CONE
###
function lower_interval_bound(m::Optimizer, d::BufferedSOC, n::NodeBB)

    sol_branch_map = m._sol_to_branch_map
    lo_bnds = n.lower_variable_bound
    up_bnds = n.upper_variable_bound
    vec_of_vi = d.variables.variables

    norm_bound = Interval(0.0)
    for i = 2:length(vec_of_vi)
        mapped_vi = @inbounds sol_branch_map[vec_of_vi[i].value]
        x = Interval{Float64}(lo_bnds[mapped_vi], up_bnds[mapped_vi])
        norm_bound += pow(x, 2)
    end
    norm_bound = sqrt(norm_bound)

    mapped_vi = @inbounds sol_branch_map[vec_of_vi[1].value]
    lower_bound = norm_bound.lo - (@inbounds up_bnds[mapped_vi])

    return lower_bound
end


###
### NONLINEAR FUNCTIONS
###
function lower_interval_bound(m::Optimizer, d::BufferedNonlinearFunction{V}, n::NodeBB) where V
    !d.has_value && forward_pass!(m._working_problem._relaxed_evaluator, d)
    ex = d.expr
    v = ex.isnumber[1] ? ex.numberstorage[1] : ex.setstorage[1].Intv.lo
    return v
end
function interval_bound(m::Optimizer, d::BufferedNonlinearFunction{V}, n::NodeBB) where V
    !d.has_value && forward_pass!(d.evaluator, d)
    ex = d.expr
    v = ex.isnumber[1] ? Interval(ex.numberstorage[1]) : ex.setstorage[1].Intv
    return v.lo, v.hi
end

_is_feas(m::GlobalOptimizer{N,T,S}, x::AffineFunctionIneq, n) where {N,T<:AbstractFloat,S} = lower_interval_bound(m, x, n) <= 0.0
_is_feas(m::GlobalOptimizer{N,T,S}, x::BufferedQuadraticIneq, n) where {N,T<:AbstractFloat,S} = lower_interval_bound(m, x, n) <= 0.0
function _is_feas(m::GlobalOptimizer{N,T,S}, x::AffineFunctionEq, n) where {N,T<:AbstractFloat,S}
    lower_value, upper_value = interval_bound(m, x, n)
    return lower_value <= 0.0 <= upper_value
end
function _is_feas(m::GlobalOptimizer{N,T,S}, x::BufferedQuadraticEq, n) where {N,T<:AbstractFloat,S}
    lower_value, upper_value = interval_bound(m, x, n)
    return lower_value <= 0.0 <= upper_value
end
function _is_feas(m::GlobalOptimizer{N,T,S}, x::NonlinearExpression, n) where {N,T<:AbstractFloat,S}
    lower_value, upper_value = interval_bound(m, x, n)
    feasible_flag &= upper_value < x.lower_bound
    feasible_flag &= lower_value > x.upper_bound
    return !feasible_flag
end
