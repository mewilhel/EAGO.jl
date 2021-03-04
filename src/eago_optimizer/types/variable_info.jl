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
