"""
$(TYPEDEF)

A structure used to store information related to the bounds assigned to each
variable.

$(TYPEDFIELDS)
"""
Base.@kwdef struct VariableInfo
    "Is the variable integer valued?"
    is_integer::Bool                                            = false
    "Boolean indicating whether finite lower bound exists."
    has_lower_bound::Bool                                       = false
    "Boolean indicating whether finite upper bound exists."
    has_upper_bound::Bool                                       = false
    "Boolean indicating variable is fixed to a finite value."
    is_fixed::Bool                                              = false
    "Indicates that no constraints have been set"
    has_constraints::Bool                                       = false
    "Lower bounds. May be -Inf."
    lower_bound::Float64                                        = -Inf
    "Upper bounds. May be Inf."
    upper_bound::Float64                                        = Inf
    ""
end
_is_integer(x::VariableInfo) = x.is_integer
_has_lower_bound(x::VariableInfo) = x.has_lower_bound
_has_upper_bound(x::VariableInfo) = x.has_upper_bound
_is_fixed(x::VariableInfo) = x.is_fixed
_lower_bound(x::VariableInfo) = x.lower_bound
_upper_bound(x::VariableInfo) = x.upper_bound

function VariableInfo(::ZO)
    return VariableInfo(is_integer = true,
                        has_lower_bound = true,
                        has_upper_bound = true,
                        has_constraints = true,
                        lower_bound = 0.0,
                        upper_bound = 1.0)
end

VariableInfo(lt::LT) = VariableInfo(has_upper_bound = true,
                                    has_constraints = true,
                                    upper_bound = lt.upper)

VariableInfo(gt::GT) = VariableInfo(has_lower_bound = true,
                                    has_constraints = true,
                                    lower_bound = gt.lower)

VariableInfo(et::ET) = VariableInfo(has_lower_bound = true,
                                    has_upper_bound = true,
                                    has_constraints = true,
                                    is_fixed = true
                                    lower_bound = et.value,
                                    upper_bound = et.value)
