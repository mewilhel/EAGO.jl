
function is_linear(::Val{scalar}, ex::Expr) # standard scalar linear
    return false
end
function is_linear(::Val{L1}, ex::Expr)     # L1 norm
    return false
end
function is_linear(::Val{Linf}, ex::Expr)   # Linf norm
    return false
end
function is_linear(::Val{min}, ex::Expr)    # Linf min
    return false
end
function is_linear(::Val{max}, ex::Expr)    # Linf max
    return
end

function is_linear(ex::Expr)
    is_linear(Val(:scalar), ex) && (return true)
    is_linear(Val(:L1), ex)
    is_linear(Val(:Linf), ex)
    is_linear(Val(:min), ex)
    is_linear(Val(:max), ex)
end
