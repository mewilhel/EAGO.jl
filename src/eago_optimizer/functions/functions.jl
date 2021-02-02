# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/functions/nonlinear/functions.jl
# Defines variable info and function types.
#############################################################################

"""
$(TYPEDEF)

An abstract super-type used for representing constraints built by EAGO's backend.
"""
abstract type AbstractEAGOConstraint end

"""

Computes the lower interval bound for `AbstractEAGOConstraint` representing an
inequality constraint.
"""
function lower_interval_bound end

"""

Computes a tuple representing the lower and upper interval bounds for a
`AbstractEAGOConstraint` representing an equality constraint.
"""
function interval_bound end

"""

Eliminate fixed variables by rearrangment or restructuring of `AbstractEAGOConstraint`.
"""
function eliminate_fixed_variables! end

@enum(BRANCH_VARIABLE, UNSPEC_BRANCH, NO_BRANCH, BRANCH)

"""
$(TYPEDEF)

A structure used to store information related to the bounds assigned to each
variable.

$(TYPEDFIELDS)
"""
mutable struct VariableInfo
    "Is the variable integer valued?"
    is_integer::Bool
    "Lower bounds. May be -Inf."
    lower_bound::Float64
    "Boolean indicating whether finite lower bound exists."
    has_lower_bound::Bool
    "Upper bounds. May be Inf."
    upper_bound::Float64
    "Boolean indicating whether finite upper bound exists."
    has_upper_bound::Bool
    "Boolean indicating variable is fixed to a finite value."
    is_fixed::Bool
    "Is variable used for branching (unset in input model)"
    branch_on::BRANCH_VARIABLE
end
VariableInfo() = VariableInfo(false,-Inf, false, Inf, false, false, UNSPEC_BRANCH)
lower_bound(x::VariableInfo) = x.lower_bound
upper_bound(x::VariableInfo) = x.upper_bound
is_integer(x::VariableInfo) = x.is_integer

include("affine.jl")
include("quadratic.jl")
include("cone.jl")
include("nonlinear/nonlinear.jl")
