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
struct AffineFunction{T<:Real} <: AbstractEAGOConstraint
    func::SAF{T}
    len::Int
end
function AffineFunction(func::SAF{T}, c::T) where {T<:Real}
    fc = copy(func)
    fc.constant = c
    return AffineFunction(fc, length(fc.terms))
end
AffineFunction{T}() where {T<:Real} = AffineFunction{T}(SAF{T}(SAT{T}[], zero(T)), zero(T))
AffineFunction(f::SAF{T}, s::LT{T}) where {T<:Real} = AffineFunction(f, f.constant - s.upper)
AffineFunction(f::SAF{T}, s::GT{T}) where {T<:Real} = AffineFunction(f, s.lower - f.constant)
AffineFunction(f::SAF{T}, s::ET{T}) where {T<:Real} = AffineFunction(f, f.constant - s.value)

###
### Parsing definitions
###

function eliminate_fixed_variables(f::AffineFunction{T}, v::Vector{VariableInfo{T}}) where T
    func = f.func
    len = f.len
    deleted_count = 0
    for i = 1:len
        sat = @inbounds func.terms[i - deleted_count]
        vi = @inbounds v[sat.value]
        if vi.is_fixed
            f.constant += (sat.coefficient)*vi.lower_bound
            deleteat!(func.terms, i - deleted_count)
            deleted_count += 1
        end
    end
    len -= deleted_count
    return AffineFunction(func, len)
end
