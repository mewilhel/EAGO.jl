# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/relax.jl
# Defines routines used construct the relaxed subproblem.
#############################################################################

"""
$(FUNCTIONNAME)

Applies the safe cut checks detailed in Khajavirad, 2018 [Khajavirad, Aida,
and Nikolaos V. Sahinidis. "A hybrid LP/NLP paradigm for global optimization
relaxations." Mathematical Programming Computation 10.3 (2018): 383-421] to
ensure that only numerically safe affine relaxations are added. Checks that
i) ``|b| <= safe b`, ii) `safe_l <= abs(ai) <= safe u`, and iii) violates
`safe_l <= abs(ai/aj) <= safe_u`.
"""
function is_safe_cut!(m::GlobalOptimizer{N,T,S}, f::SAF) where {N,T<:AbstractFloat,S}

    safe_l = m.cut_safe_l
    safe_u = m.cut_safe_u

    # violates |b| <= safe_b
    (abs(f.constant) > m.cut_safe_b) && return false

    term_count = length(f.terms)
    for i = 1:term_count
        ai = (@inbounds f.terms[i]).coefficient
        if !iszero(ai)
            ai_abs = abs(ai)
            # violates safe_l <= abs(ai) <= safe_u
            !(safe_l <= abs(ai) <= safe_u) && return false
            for j = i:term_count
                aj = (@inbounds f.terms[j]).coefficient
                if !iszero(aj)
                    # violates safe_l <= abs(ai/aj) <= safe_u
                    !(safe_l <= abs(ai/aj) <= safe_u) && return false
                end
            end
        end
    end

    return true
end

"""
$(FUNCTIONNAME)

Relaxs the constraint by adding an affine constraint to the model.
"""
function relax! end

"""
$(FUNCTIONNAME)

Default routine for relaxing quadratic constraint `func < 0.0` on node `n`.
Takes affine bounds of convex part at point `x0` and secant line bounds on
concave parts.
"""
function affine_relax_quadratic!(func::SQF, d::Dict{Int,Float64}, saf::SAF,
                                 n::NodeBB, sol_to_branch_map::Vector{Int},
                                 x::Vector{Float64})

    quadratic_constant = func.constant

    # Affine terms only contribute coefficients, so the respective
    # values do not contribute to the cut. Since all quadratic terms
    # are considered to be branch variables we exclude any potential
    # need to retrieve variable bounds from locations other than
    # the node.
    for term in func.quadratic_terms

        a = term.coefficient
        idx1 = term.variable_index_1.value
        idx2 = term.variable_index_2.value
        x1 = x[sol_idx1]                    # TODO FIX THIS
        xL1 = _lower_bound(FullVar, m, idx1)
        xU1 = _upper_bound(FullVar, m, idx1)

        if idx1 == idx2
            if a > zero(T)
                d[idx1] += a*x1
                quadratic_constant -= 0.5*a*x1*x1
            else
                if !isinf(xL1) && !isinf(xU1)
                    d[idx1] += 0.5*a*(xL1 + xU1)
                    quadratic_constant -= 0.5*a*xL1*xU1
                else
                    return false
                end
            end
        else
            x2 = x[sol_idx2]
            xL2 = _lower_bound(FullVar, m, idx2)
            xU2 = _upper_bound(FullVar, m, idx2)
            if a > 0.0
                if (!isinf(xL1) && !isinf(xL2)) &&
                   ((xU1 - xL1)*x2 + (xU2 - xL2)*x1 <= xU1*xU2 - xL1*xL2)
                    d[idx1] += a*xL2
                    d[idx2] += a*xL1
                    quadratic_constant -= a*xL1*xL2
                elseif !isinf(xU1) && !isinf(xU2)
                    d[idx1] += a*xU2
                    d[idx2] += a*xU1
                    quadratic_constant -= a*xU1*xU2
                else
                    return false
                end
            else
                if (!isinf(xU1) && !isinf(xL2)) &&
                   ((xU1 - xL1)*x2 - (xU2 - xL2)*x1 <= xU1*xL2 - xL1*xU2)
                    d[idx1] += a*xL2
                    d[idx2] += a*xU1
                    quadratic_constant -= a*xU1*xL2
                elseif !isinf(xL1) && !isinf(xU2)
                    d[idx1] += a*xU2
                    d[idx2] += a*xL1
                    quadratic_constant -= a*xL1*xU2
                else
                    return false
                end
            end
        end
    end
    foreach(x -> (d[x.variable_index.value] = x.coefficient;), func.affine_terms)

    i = 1
    for (k, v) in d
        saf.terms[i] = SAT(v, VI(k))
        d[k] = 0.0
        i += 1
    end
    saf.constant = quadratic_constant

    return true
end

"""
$(TYPEDSIGNATURES)
"""
function relax!(m::Optimizer, f::BufferedQuadraticIneq, check_safe::Bool)

    constraint_tol = m.absolute_constraint_feas_tolerance
    finite_cut_generated = affine_relax_quadratic!(f.func, f.buffer, f.saf, m._current_node, m._sol_to_branch_map, m._current_xref)
    if finite_cut_generated
        if !check_safe || is_safe_cut!(m, f.saf)
            lt = LT(-f.saf.constant + constraint_tol)
            f.saf.constant = 0.0
            ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, lt)
            push!(m._buffered_quadratic_ineq_ci, ci)
        end
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function relax!(m::Optimizer, f::BufferedQuadraticEq, check_safe::Bool)

    constraint_tol = m.absolute_constraint_feas_tolerance
    finite_cut_generated = affine_relax_quadratic!(f.func, f.buffer, f.saf, m._current_node, m._sol_to_branch_map, m._current_xref)
    if finite_cut_generated
        if !check_safe || is_safe_cut!(m, f.saf)
            lt = LT(-f.saf.constant + constraint_tol)
            f.saf.constant = 0.0
            ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, lt)
            push!(m._buffered_quadratic_eq_ci, ci)
        end
    end

    finite_cut_generated = affine_relax_quadratic!(f.minus_func, f.buffer, f.saf, m._current_node, m._sol_to_branch_map, m._current_xref)
    if finite_cut_generated
        if !check_safe || is_safe_cut!(m, f.saf)
            lt = LT(-f.saf.constant + constraint_tol)
            f.saf.constant = 0.0
            ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, lt)
            push!(m._buffered_quadratic_eq_ci, ci)
        end
    end

    return nothing
end

"""
$(FUNCTIONNAME)
"""
function affine_relax_nonlinear!(f::BufferedNonlinearFunction{MC{N,T}}, evaluator::Evaluator,
                                 use_cvx::Bool, new_pass::Bool) where {N,T<:RelaxTag}

    new_pass && forward_pass!(evaluator, f)
    x = evaluator.x
    finite_cut = true

    if f.expr.isnumber[1]
        f.saf.constant = f.expr.numberstorage[1]
        map!(i -> SAT(0.0, VI(i)), f.saf.terms, f.expr.grad_sparsity)
    else
        setvalue = f.expr.setstorage[1]
        finite_cut &= !(isempty(setvalue) || isnan(setvalue))
        if finite_cut
            value = f.expr.setstorage[1]
            f.saf.constant = use_cvx ? value.cv : -value.cc
            for i = 1:N
                vi = f.expr.grad_sparsity[i]
                c = use_cvx ? value.cv_grad[i] : -value.cc_grad[i]
                f.saf.terms[i] = SAT(c, VI(vi))
                f.saf.constant = sub_round(f.saf.constant , mul_round(coef, x[vi], RoundUp), RoundDown)
            end
            bnd_used =  use_cvx ? -f.upper_bound : f.lower_bound
            f.saf.constant = add_round(f.saf.constant, bnd_used, RoundDown)
        end
    end
    return finite_cut
end

"""
$(TYPEDSIGNATURES)
"""
function check_set_affine_nl!(m::Optimizer, f::BufferedNonlinearFunction{MC{N,T}}, finite_cut_generated::Bool, check_safe::Bool) where {N,T<:RelaxTag}
    if finite_cut_generated && (!check_safe || is_safe_cut!(m, f.saf))
        lt = LT(-f.saf.constant + m.absolute_constraint_feas_tolerance)
        f.saf.constant = 0.0
        ci = MOI.add_constraint(m.relaxed_optimizer, f.saf, lt)
        push!(m._buffered_nonlinear_ci, ci)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function relax!(m::Optimizer, f::BufferedNonlinearFunction{MC{N,T}}, check_safe::Bool) where {N,T<:RelaxTag}
    evaluator = m._working_problem._relaxed_evaluator

    finite_cut_generated = affine_relax_nonlinear!(f, evaluator, true, true)
    check_set_affine_nl!(m, f, finite_cut_generated, check_safe)

    finite_cut_generated = affine_relax_nonlinear!(f, evaluator, false, false)
    check_set_affine_nl!(m, f, finite_cut_generated, check_safe)

    return nothing
end

"""
$(FUNCTIONNAME)

Adds linear objective cut constraint to the `x.relaxed_optimizer`.
"""
function objective_cut!(m::Optimizer, check_safe::Bool)
    u = m._global_upper_bound
    if u < Inf
        wp = m._working_problem
        if check_safe && is_safe_cut!(m, wp._objective)
            val = u - wp._objective.constant + m.abs_constraint_feas_tol
            MOI.set(m.mip_solver, MOI.ConstraintSet(), m._obj_ci, LT(val))
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

A routine that adds relaxations for all nonlinear constraints and quadratic constraints
corresponding to the current node to the relaxed problem. This adds an objective cut
(if specified by `objective_cut_on`) and then sets the `_new_eval_constraint` flag
to false indicating that an initial evaluation of the constraints has occurred. If
the `objective_cut_on` flag is `true` then the `_new_eval_objective` flag is also
set to `false` indicating that the objective expression was evaluated.
"""
function relax_all_constraints!(t::ExtensionType, m::Optimizer, q::Int64)

    wp = _working_problem(m)
    check_safe = (q === 1) ? false : m.cut_safe_on
    wp._relaxed_evaluator.is_first_eval = m._new_eval_constraint

    foreach(x -> relax!(m, x, check_safe), wp._sqf_leq)
    foreach(x -> relax!(m, x, check_safe), wp._sqf_eq)
    foreach(x -> relax!(m, x, check_safe), wp._nonlinear_constr)

    m._new_eval_constraint = false
    objective_cut!(m, check_safe)

    return nothing
end
relax_constraints!(t::ExtensionType, m::Optimizer, q::Int64) = relax_all_constraints!(t, m, q)
relax_constraints!(m::Optimizer, q::Int64) = relax_constraints!(m.ext_type, m, q)

"""
$(FUNCTIONNAME)

Deletes all nonlinear constraints added to the relaxed optimizer.
"""
function delete_nl_constraints!(m::Optimizer)
    MOI.delete.(m.relaxed_optimizer, m._buffered_quadratic_ineq_ci)
    MOI.delete.(m.relaxed_optimizer, m._buffered_quadratic_eq_ci)
    MOI.delete.(m.relaxed_optimizer, m._buffered_nonlinear_ci)
    empty!(m._buffered_quadratic_ineq_ci)
    empty!(m._buffered_quadratic_eq_ci)
    empty!(m._buffered_nonlinear_ci)
    return
end

"""
$(FUNCTIONNAME)

"""
function set_first_relax_point!(m::Optimizer)

    m._working_problem._relaxed_evaluator.is_first_eval = true
    m._new_eval_constraint = true
    m._new_eval_objective = true

    n = m._current_node
    @__dot__ m._current_xref = 0.5*(n.upper_variable_bounds + n.lower_variable_bounds)
    unsafe_check_fill!(isnan, m._current_xref, 0.0, length(m._current_xref))

    return nothing
end
