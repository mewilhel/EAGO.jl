
function _label_branch!(d, m, ci::CI{SQF,S}) where S
    sqf = MOI.get(m, MOI.ConstraintFunction(), ci)
    for term in sqf.quadratic_terms
        d[term.variable_index_1.value] = true
        d[term.variable_index_2.value] = true
    end
    return
end
function _label_branch!(d, m, ci::CI{VECVAR,SOC_CONE}, start::Int = 1)
    vv = MOI.get(m, MOI.ConstraintFunction(), ci)
    for i = start:length(vv)
        d[vv.variables[i].value] = true
    end
    return
end

# Iterative deepening depth-first search for variable references
_label_branch_nl!(d, expr::T) where T<:Number = nothing
function _label_branch_nl!(d, ex::Expr)
    if ex.head == :ref
        d[ex.args[2]] = true
        return
    end
    for i = 2:length(ex.args)
        _label_branch_nl!(d, ex_arg)
    end
    return
end

_label_nl!(d::Dict{Int,Bool}, ::Nothing) = nothing
function _label_nl!(d::Dict{Int,Bool}, nl_data)
    nl_eval = nl_data.evaluator
    nl_bnds = nl_data.constraint_bounds
    MOI.initialize(nl_eval, Symbol[:ExprGraph])

    if nl_data.has_objective
        obj_expr = MOI.objective_expr(nl_eval)
        _label_branch_nl!(d, obj_expr)
    end
    for i = 1:length(nl_bnds)
        _label_branch_nl!(d, MOI.constraint_expr(nl_eval, i))
    end
    return nothing
end
