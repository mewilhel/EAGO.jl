
function _is_exp_form(::Val{:exp1}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:exp2}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:exp3}, ex::Expr)
   return false
end

function _is_exp_form(::Val{:log1}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:log2}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:log3}, ex::Expr)
   return false
end

function _is_exp_form(::Val{:arbexp}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:arblog}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:logsuminv}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:softplus}, ex::Expr)
   return false
end

function _is_exp_form(::Val{:alg1}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:alg2}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:alg3}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:alg4}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:alg5}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:alg6}, ex::Expr)
   return false
end
function _is_exp_form(::Val{:alg7}, ex::Expr)
   return false
end

function _is_exp(ex::Expr)

   # basic exp cones
   _is_exp_form(Val(:exp1), ex) && (return true) # is x_2*exp(x_3/x_2) - x_1 <= 0.0
   _is_exp_form(Val(:exp2), ex) && (return true) # is exp(x_3/x_2) - x_1/x_2 <= 0.0
   _is_exp_form(Val(:exp3), ex) && (return true) # is (x_2/x_1)*exp(x_3/x_2) - 1 <= 0.0

   # log exp cones
   _is_exp_form(Val(:log1), ex) && (return true) # is x_3 - x_2*log(x_1/x_2) <= 0.0
   _is_exp_form(Val(:log2), ex) && (return true) # is (x_3/x_2) - log(x_1/x_2) <= 0.0
   _is_exp_form(Val(:log3), ex) && (return true) # is 1 - (x_2/x_3)*log(x_1/x_2) <= 0.0

   _is_exp_form(Val(:arbexp), ex) && (return true)    # is arbitrary exponential
   _is_exp_form(Val(:arblog), ex) && (return true)    # is arbitrary logarithm
   _is_exp_form(Val(:logsuminv), ex) && (return true) # is log sum inv

   _is_exp_form(Val(:softplus), ex) && (return true)  # is softplus

   _is_exp_form(Val(:alg1), ex) && (return true)     # log(x)^2 - t,    0 < x < 1
   _is_exp_form(Val(:alg2), ex) && (return true)     # t - log(log(x)), x > 1
   _is_exp_form(Val(:alg3), ex) && (return true)     # log(x)^-1 - t,   x > 1
   _is_exp_form(Val(:alg4), ex) && (return true)     # log(x)^2 - t, 0 < x < 1
   _is_exp_form(Val(:alg5), ex) && (return true)     # log(x)^2 - t, 0 < x < 1
   _is_exp_form(Val(:alg6), ex) && (return true)     # log(x)^2 - t, 0 < x < 1
   _is_exp_form(Val(:alg7), ex) && (return true)     # log(x)^2 - t, 0 < x < 1

    return false
end
