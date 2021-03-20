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
_lower_bound(x::VariableInfo) = x.lower_bound
_upper_bound(x::VariableInfo) = x.upper_bound

_is_fixed(x::VariableInfo) = x.is_fixed
function _is_less_than(x::VariableInfo)
    flag = x.has_upper_bound
    flag &= !x.has_lower_bound
    flag &= !x.is_integer
    return flag
end
function _is_greater_than(x::VariableInfo)
    flag = x.has_lower_bound
    flag &= !x.has_upper_bound
    flag &= !x.is_integer
    return flag
end
function _is_zero_one(x::VariableInfo)
    flag = iszero(x.lower_bound)
    flag &= isone(x.upper_bound)
    flag &= x.is_integer
    return flag
end
function _is_int_interval(x::VariableInfo)
    flag = x.has_lower_bound)
    flag &= x.has_upper_bound
    flag &= x.is_integer
    return flag
end
function _is_real_interval(x::VariableInfo)
    flag = x.has_lower_bound)
    flag &= x.has_upper_bound
    flag &= !x.is_integer
    return flag
end

mid(x::VariableInfo) = 0.5*(_upper_bound(x) - _lower_bound(x))

empty_variable_info() = VariableInfo(lower_bound = Inf,
                                     upper_bound = -Inf)

Base.isempty(v::VariableInfo) = _lower_bound(v) > _upper_bound(v)
function check_isempty(l, u, b)
    flag = l < u
    if b
        flag &= (l <= 0.0) | (u >= 1.0)
    end
    return !flag
end

function VariableInfo(::ZO)
    return VariableInfo(is_integer = true,
                        has_lower_bound = true,
                        has_upper_bound = true,
                        has_constraints = true,
                        lower_bound = 0.0,
                        upper_bound = 1.0)
end

function VariableInfo(it::IT)
    VariableInfo(has_lower_bound = !isinf(it.lower),
                 has_upper_bound = !isinf(it.upper),
                 has_constraints = !isinf(it.lower) | !isinf(it.upper),
                 is_fixed = it.lower == it.upper,
                 lower_bound = it.lower,
                 upper_bound = it.upper)
end

function VariableInfo(v::VariableInfo, ::ZO)
    isempty(v) && (return v)
    l = max(0.0, _lower_bound(v))
    u = min(1.0, _upper_bound(v))
    if check_isempty(l, u, _is_integer(v))
        return empty_variable_info()
    end
    return VariableInfo(is_integer = true,
                        has_lower_bound = true,
                        has_upper_bound = true,
                        has_constraints = true,
                        is_fixed = l == u,
                        lower_bound = l,
                        upper_bound = u)
end

function VariableInfo(v::VariableInfo, it::IT)
    isempty(v) && (return v)
    l = max(it.lower, _lower_bound(v))
    u = max(it.upper, _upper_bound(v))
    if check_isempty(l, u, _is_integer(v))
        return empty_variable_info()
    end
    return VariableInfo(is_integer = _is_integer(v),
                        has_lower_bound = !isinf(l),
                        has_upper_bound = !isinf(u),
                        has_constraints = !isinf(l) | !isinf(u),
                        is_fixed = l == u,
                        lower_bound = l,
                        upper_bound = u)
end

_ET(v::VariableInfo)  = ET(v.lower_bound)
_IT(v::VariableInfo)  = IT(v.lower_bound)
_GT(v::VariableInfo)  = GT(v.lower_bound)
_LT(v::VariableInfo)  = LT(v.upper_bound)
_ZO(v::VariableInfo)  = ZO()
_INT(v::VariableInfo) = MOI.Semiinteger(v.lower_bound, v.upper_bound)
