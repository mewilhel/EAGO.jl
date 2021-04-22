
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

is_scaled_ref(d::T) where T <: Real = true
function is_scaled_ref(d::Expr)
    is_affine = true
    h = ex.head
    if h == :*
        is_affine = is_scaled_ref(d.args[1]) && is_scaled_ref(d.args[1])
    elseif h != :ref
        is_affine = false
    end
    return is_affine
end
is_affine_expr(d::T) where T <: Real = true
function is_affine_expr(d::Expr)
    is_affine = true
    h = ex.head
    if h == :+
        is_affine = all(is_scaled_ref, d.args)
    elseif h == :*
        is_affine = is_scaled_ref(d)
    elseif h != :ref
        is_affine = false
    end
    return is_affine
end

_label_nl(d::Dict{Int,Bool}, ::Nothing) = false
function _label_nl(d::Dict{Int,Bool}, nl_data)
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
    return is_affine_expr(obj_expr)
end
