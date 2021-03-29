# check ex1 is norm, ex2 is prod of exponents, exponents = 1.0
function _inner_power_cone_chk2(ex1::Expr, ex2::Expr)
end

# check sqrt(...) - prod(...) <= 0.0
function _inner_power_cone_chk(ex::Expr)
    if (ex.head == :call) && (ex.args[1] == :-) && (length(ex.args) == 3)
        one_ex = ex.args[2]
        two_ex = ex.args[3]
        one_is_ex = one_ex isa Expr
        two_is_ex = two_ex isa Expr
        (!one_is_ex || !two_is_ex) && (return false)
        if (one_is_ex.args[1] == :sqrt) && (two_ex.args[1] == :*)
            return _inner_power_cone_chk2(one_is_ex, two_ex)
        elseif (one_is_ex.args[1] == :*) && (two_ex.args[1] == :sqrt)
            return _inner_power_cone_chk2(two_expr, one_is_ex)
        else
            return false
        end
    end
    return false
end

function _is_power_cone(ex::Expr)
   # flatten if necessary && make leq
   # is sqrt(x_3^2 + x_4^2 + .... x_n^2) - x_1^\alpha*x_2^(1-\alpha)...x^m leq 0
   _inner_power_cone_chk(ex)
end
