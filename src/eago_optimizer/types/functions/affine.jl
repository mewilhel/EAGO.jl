# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/functions/nonlinear/affine.jl
# Defines buffered structures to store quadratic functions:
# AffineFunctionIneq, AffineFunctionEq, as well as the
# lower_interval_bound, interval_bound, and eliminate_fixed_variables!
# functions associated with each structure.
#############################################################################

###
### Structure definitions
###

"""
$(TYPEDEF)

Current only used for bound tightening. Stores a representation
of an affine inequality.
"""
mutable struct AffineFunctionIneq{T<:AbstractFloat} <: AbstractEAGOConstraint
    func::SAF{T}
    len::Int
end
AffineFunctionIneq{T}() where {T<:AbstractFloat} = AffineFunctionIneq{T}(SAF{T}(SAT{T}[], zero(T)), zero(T))

function AffineFunctionIneq(func::SAF{T}, set::LT{T}) where {T<:AbstractFloat}
    func.constant -= set.upper
    return AffineFunctionIneq{T}(func, length(func.terms))
end
function AffineFunctionIneq(func::SAF{T}, set::GT{T}) where {T<:AbstractFloat}
    func.constant = set.lower - func.constant
    return AffineFunctionIneq{T}(func, length(func.terms))
end

"""
$(TYPEDEF)

Current only used for bound tightening. Stores a representation
of an affine equality.
"""
mutable struct AffineFunctionEq{T<:AbstractFloat} <: AbstractEAGOConstraint
    func::SAF{T}
    len::Int
end
AffineFunctionEq{T}() where T<:AbstractFloat = AffineFunctionEq{T}(Tuple{T,Int}[], zero(T), 0)
function AffineFunctionEq(func::SAF{T}, set::ET{T}) where {T<:AbstractFloat}
    func.constant -= set.value
    return AffineFunctionEq{T}(func, length(func.terms))
end

###
### Parsing definitions
###

function eliminate_fixed_variables!(f::T, v::Vector{VariableInfo}) where T <: Union{AffineFunctionIneq,
                                                                                    AffineFunctionEq}
    deleted_count = 0
    for i = 1:f.len
        sat = @inbounds f.func.terms[i - deleted_count]
        vi = @inbounds v[sat.value]
        if vi.is_fixed
            f.constant += (sat.coefficient)*vi.lower_bound
            deleteat!(f.terms, i - deleted_count)
            deleted_count += 1
        end
    end
    f.len -= deleted_count
    return nothing
end
