# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/functions/nonlinear/quadratic.jl
# Defines buffered structures to store quadratic functions:
# BufferedQuadraticIneq, BufferedQuadraticEq, as well as the
# lower_interval_bound, interval_bound, and eliminate_fixed_variables!
# functions associated with each structure.
#############################################################################

###
### Structure definitions
###

"""
$(TYPEDEF)

Stores a general quadratic inequality constraint with a buffer.
"""
mutable struct BufferedQuadraticIneq{T<:AbstractFloat} <: AbstractEAGOConstraint
    func::SQF{T}
    buffer::Dict{Int,T}
    saf::SAF{T}
    len::Int
end

"""
$(TYPEDEF)

Stores a general quadratic equality constraint with a buffer.
"""
mutable struct BufferedQuadraticEq{T<:AbstractFloat} <: AbstractEAGOConstraint
    func::SQF{T}
    minus_func::SQF{T}
    buffer::Dict{Int,T}
    saf::SAF{T}
    len::Int
end

#=
mutable struct BufferedConvexQuadratic <: AbstractEAGOConstraint
    func::SQF
    buffer::Dict{Int, Float64}
    saf::SAF
    len::Int
end
=#

###
### Constructor definitions
###

function create_buffer_dict(func::SQF{T}) where {T<:AbstractFloat}

    d = Dict{Int,T}()

    for term in func.quadratic_terms
        d[term.variable_index_1.value] = zero(T)
        d[term.variable_index_2.value] = zero(T)
    end

    for term in func.affine_terms
        d[term.variable_index.value] = zero(T)
    end

    return d
end

function BufferedQuadraticIneq{T}() where T<:AbstractFloat
    BufferedQuadraticIneq{T}(SQF{T}(SQT{T}[], SAT{T}[], zero(T)),
                             Dict{Int, T}(), SAF{T}(SAT{T}[], zero(T)), zero(T))
end

function BufferedQuadraticIneq(func::SQF{T}, set::LT{T}) where T<:AbstractFloat
    b = create_buffer_dict(func)
    saf = SAF{T}([SAT{T}(zero(T), VI(k)) for k in keys(b)], zero(T))
    cfunc = copy(func)
    cfunc.constant -= set.upper
    return BufferedQuadraticIneq{T}(cfunc, b, saf, length(b))
end

function BufferedQuadraticIneq(func::SQF{T}, set::GT{T}) where T<:AbstractFloat
    b = create_buffer_dict(func)
    saf = SAF{T}([SAT{T}(zero(T), VI(k)) for k in keys(b)], zero(T))
    cfunc = MOIU.operate(-, T, func)
    cfunc.constant += set.lower
    return BufferedQuadraticIneq(cfunc, b, saf, length(b))
end

function BufferedQuadraticEq{T}() where T<:AbstractFloat
    BufferedQuadraticEq{T}(SQF{T}(SQT{T}[], SAT{T}[], zero(T)),
                           SQF{T}(SQT{T}[], SAT{T}[], zero(T)),
                           Dict{Int,T}(), SAF{T}(SAT{T}[], zero(T)), zero(T))
end

function BufferedQuadraticEq(func::SQF{T}, set::ET{T}) where T<:AbstractFloat
    b = create_buffer_dict(func)
    saf = SAF{T}([SAT{T}(zero(T), VI(k)) for k in keys(b)], zero(T))
    cfunc1 = copy(func)
    cfunc1.constant -= set.value
    cfunc2 = MOIU.operate(-, T, func)
    cfunc2.constant += set.value
    return BufferedQuadraticEq(cfunc1, cfunc2, buffer, saf, length(b))
end

#=
function BufferedConvexQuadratic(f::BufferedQuadraticIneq)
    BufferedConvexQuadratic(copy(f.func), copy(f.buffer), copy(f.saf), f.len)
end
=#

###
### Parsing definitions
###

function eliminate_fixed_variables!(f::T, v::Vector{VariableInfo}) where T <: Union{BufferedQuadraticIneq,
                                                                                    BufferedQuadraticIneq}
    deleted_count = 0
    index = 1
    while i + deleted_count <= f.len
        term = @inbounds f.sqf.terms[i]
        variable_info_1 = @inbounds v[term.variable_index_1.value]
        variable_info_2 = @inbounds v[term.variable_index_2.value]
        if variable_info_1.is_fixed && variable_index_2.is_fixed
            f.sqf.constant += coeff*variable_info_1.lower_bound*variable_index_2.lower_bound
            deleteat!(f.sqf.terms, i)
            deleted_count += 1
        else
            i += 1
        end
    end
    f.len -= deleted_count

    return nothing
end
