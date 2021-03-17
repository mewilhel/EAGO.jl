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
mutable struct AffineFunctionIneq <: AbstractEAGOConstraint
    func::SAF
    len::Int
end
AffineFunctionIneq() = AffineFunctionIneq(SAF(SAT[],0.0), 0)
function AffineFunctionIneq(func::SAF, set::LT)
    func.constant -= set.upper
    return AffineFunctionIneq(func, length(func.terms))
end

function AffineFunctionIneq(func::SAF, set::GT)
    func.constant = set.lower - func.constant
    return AffineFunctionIneq(func, length(func.terms))
end

"""
$(TYPEDEF)

Current only used for bound tightening. Stores a representation
of an affine equality.
"""
mutable struct AffineFunctionEq <: AbstractEAGOConstraint
    func::SAF
    len::Int
end
AffineFunctionEq() = AffineFunctionEq(Tuple{Float64,Int}[], 0.0, 0)
function AffineFunctionEq(func::SAF, set::ET)
    func.constant -= set.value
    return AffineFunctionEq(func, length(func.terms))
end

###
### Parsing definitions
###

function eliminate_fixed_variables!(f::T, v::Vector{VariableInfo}) where T <: Union{AffineFunctionIneq,
                                                                                    AffineFunctionEq}
    for i = 1:f.len
        sat = @inbounds f.func.terms[i - deleted_count]
        coeff = sat.coefficient
        indx = sat.value
        variable_info = @inbounds v[indx]
        if variable_info.is_fixed
            f.constant += coeff*variable_info.lower_bound
            deleteat!(f.terms, i - deleted_count)
            deleted_count += 1
        end
    end
    f.len -= deleted_count
    return nothing
end
