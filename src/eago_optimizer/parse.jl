# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# Defines functions used to parse the input optimization problem into
# a solvable form including routines used to classify input problems as a
# LP, SOCP, MILP, MISOCP, and convex problem types.
#############################################################################


function _max_to_min!(wp::ParsedProblem, obj::SV)
    wp._objective_type = SCALAR_AFFINE
    wp._objective_saf = MOIU.operate(-, Float64, obj)
    wp._objective_saf_parsed = AffineFunctionIneq(wp._objective_saf, LT_ZERO)
    return
end
function _max_to_min!(wp::ParsedProblem, obj::SAF)
    wp._objective_saf = MOIU.operate(-, Float64, wp._objective_saf)
    wp._objective_saf_parsed = AffineFunctionIneq(wp._objective_saf, LT_ZERO)
    return
end
function _max_to_min!(wp::ParsedProblem, obj::SQF)
    obj.sqf = MOIU.operate(-, Float64, obj.sqf)
    return
end
function _max_to_min!(wp::ParsedProblem)
    # updates tape for nlp_data block (used by local optimizer)
    nd_nlp = wp._nlp_data.evaluator.m.nlp_data.nlobj.nd
    pushfirst!(nd_nlp, NodeData(JuMP._Derivatives.CALLUNIVAR, 2, -1))
    nd_nlp[2] = NodeData(nd_nlp[2].nodetype, nd_nlp[2].index, 1)
    for i = 3:length(nd_nlp)
        @inbounds nd_nlp[i] = NodeData(nd_nlp[i].nodetype, nd_nlp[i].index, nd_nlp[i].parent + 1)
    end

    # updates tape used by evaluator for the nonlinear objective (used by the relaxed optimizer)
    nd_obj = wp._objective_nl.expr.nd
    pushfirst!(nd_obj, NodeData(JuMP._Derivatives.CALLUNIVAR, 2, -1))
    nd_obj[2] = NodeData(nd_obj[2].nodetype, nd_obj[2].index, 1)
    for i = 3:length(nd_obj)
        @inbounds nd_obj[i] = NodeData(nd_obj[i].nodetype, nd_obj[i].index, nd_obj[i].parent + 1)
    end
    I, J, V = findnz(wp._objective_nl.expr.adj)
    I .+= 1
    J .+= 1
    pushfirst!(I, 2)
    pushfirst!(J, 1)
    pushfirst!(V, true)
    wp._objective_nl.expr.adj = sparse(I, J, V)

    set_val = copy(wp._objective_nl.expr.setstorage[1])
    pushfirst!(wp._objective_nl.expr.setstorage, set_val)
    pushfirst!(wp._objective_nl.expr.numberstorage, 0.0)
    pushfirst!(wp._objective_nl.expr.isnumber, false)
    return
end

"""
$(TYPEDSIGNATURES)

Converts `MOI.MAX_SENSE` objective to equivalent `MOI.MIN_SENSE` objective
`max(f) = -min(-f)`.
"""
function _max_to_min!(m::Optimizer)
    wp = m._working_problem
    wp._optimization_sense = MOI.MIN_SENSE
    if _optimization_sense(m._input_problem) == MOI.MAX_SENSE
        obj_type = _objective_type(m._input_problem)
        (obj_type == SINGLE_VARIABLE)   && _max_to_min!(wp, wp._objective_sv)
        (obj_type == SCALAR_AFFINE)     && _max_to_min!(wp, wp._objective_saf)
        (obj_type == SCALAR_QUADRATIC)  && _max_to_min!(wp, wp._objective_sqf)
        (obj_type == NONLINEAR)         && _max_to_min!(wp)
    end
    return nothing
end

"""

Performs an epigraph reformulation assuming the working_problem is a minimization problem.
"""
function reform_epigraph!(m::Optimizer)

    if m._parameters.presolve_epigraph_flag
        #=
        # add epigraph variable
        obj_variable_index = MOI.add_variable(m)

        # converts ax + b objective to ax - y <= -b constraint with y objective
        obj_type = m._working_problem._objective_type
        if obj_type === SCALAR_AFFINE

            # update unparsed expression
            objective_saf = m._working_problem._objective_saf
            push!(objective_saf, SAT(-1.0, obj_variable_index))
            obj_ci = MOI.add_constraint(m, saf, LT(-objective_saf.constant))

            # update parsed expression (needed for interval bounds)

        # converts ax + b objective to ax - y <= -b constraint with y objective
        elseif obj_type === SCALAR_QUADRATIC

            # update parsed expression
            objective_sqf = m._working_problem._objective_sqf
            obj_ci = MOI.add_constraint(m, saf, LT())

        elseif obj_type === NONLINEAR

            # updated parsed expressions
            objective_nl = m._working_problem._objective_nl

        end

        MOI.set(m, MOI.ObjectiveFunction{SV}(), SV(obj_variable_index))
        =#
    end

    return nothing
end

function _check_set_is_fixed(v::VariableInfo)
    v.is_fixed && return true
    v.is_fixed = x.lower_bound === x.upper_bound
    return v.is_fixed
end

"""
$(TYPEDSIGNATURES)

Detects any variables set to a fixed value by equality or inequality constraints
and populates the `_fixed_variable` storage array.
"""
function label_fixed_variables!(m::Optimizer)
    map!(x -> _check_set_is_fixed(x), m._fixed_variable, m._working_problem._variable_info)
end

"""
$(TYPEDSIGNATURES)

Detects any variables participating in nonconvex terms and populates the
`_branch_variables` storage array.
"""
function label_branch_variables!(m::Optimizer)

    m._user_branch_variables = !isempty(m._parameters.branch_variable)
    if m._user_branch_variables
        append!(m._branch_variables, m._parameters.branch_variable)
    else

        append!(m._branch_variables, fill(false, m._working_problem._variable_count))

        # adds nonlinear terms in quadratic constraints
        sqf_leq = m._working_problem._sqf_leq
        for i = 1:m._working_problem._sqf_leq_count
            quad_ineq = @inbounds sqf_leq[i]
            for term in quad_ineq.func.quadratic_terms
                variable_index_1 = term.variable_index_1.value
                variable_index_2 = term.variable_index_2.value
                @inbounds m._branch_variables[variable_index_1] = true
                @inbounds m._branch_variables[variable_index_2] = true
            end
        end

        sqf_eq = m._working_problem._sqf_eq
        for i = 1:m._working_problem._sqf_eq_count
            quad_eq = @inbounds sqf_eq[i]
            for term in quad_eq.func.quadratic_terms
                variable_index_1 = term.variable_index_1.value
                variable_index_2 = term.variable_index_2.value
                @inbounds m._branch_variables[variable_index_1] = true
                @inbounds m._branch_variables[variable_index_2] = true
            end
        end

        obj_type = m._working_problem._objective_type
        if obj_type === SCALAR_QUADRATIC
            for term in m._working_problem._objective_sqf.func.quadratic_terms
                variable_index_1 = term.variable_index_1.value
                variable_index_2 = term.variable_index_2.value
                @inbounds m._branch_variables[variable_index_1] = true
                @inbounds m._branch_variables[variable_index_2] = true
            end
        end

        # label nonlinear branch variables (assumes affine terms have been extracted)
        nl_constr = m._working_problem._nonlinear_constr
        for i = 1:m._working_problem._nonlinear_count
            nl_constr_eq = @inbounds nl_constr[i]
            grad_sparsity = nl_constr_eq.expr.grad_sparsity
            for indx in grad_sparsity
                @inbounds m._branch_variables[indx] = true
            end
        end

        if obj_type === NONLINEAR
            grad_sparsity = m._working_problem._objective_nl.expr.grad_sparsity
            for indx in grad_sparsity
                @inbounds m._branch_variables[indx] = true
            end
        end
    end

    # add a map of branch/node index to variables in the continuous solution
    for i = 1:m._working_problem._variable_count
        if m._working_problem._variable_info[i].is_fixed
            m._branch_variables[i] = false
            continue
        end
        if m._branch_variables[i]
            push!(m._branch_to_sol_map, i)
        end
    end

    # creates reverse map
    m._sol_to_branch_map = zeros(m._working_problem._variable_count)
    for i = 1:length(m._branch_to_sol_map)
        j = m._branch_to_sol_map[i]
        m._sol_to_branch_map[j] = i
    end

    # adds branch solution to branch map to evaluator
    m._working_problem._relaxed_evaluator.node_to_variable_map = m._branch_to_sol_map
    m._working_problem._relaxed_evaluator.variable_to_node_map = m._sol_to_branch_map
    m._working_problem._relaxed_evaluator.node_count = length(m._branch_to_sol_map)

    return nothing
end


add_nonlinear_functions!(m::Optimizer) = add_nonlinear_functions!(m, m._input_problem._nlp_data.evaluator)

add_nonlinear_functions!(m::Optimizer, evaluator::Nothing) = nothing
add_nonlinear_functions!(m::Optimizer, evaluator::EmptyNLPEvaluator) = nothing
function add_nonlinear_functions!(m::Optimizer, evaluator::JuMP.NLPEvaluator)

    nlp_data = m._input_problem._nlp_data
    MOI.initialize(evaluator, Symbol[:Grad, :ExprGraph])

    # set nlp data structure
    m._working_problem._nlp_data = nlp_data

    # add subexpressions (assumes they are already ordered by JuMP)
    # creates a dictionary that lists the subexpression sparsity
    # by search each node for variables dict[2] = [2,3] indicates
    # that subexpression 2 depends on variables 2 and 3
    # this is referenced when subexpressions are called by other
    # subexpressions or functions to determine overall sparsity
    # the sparsity of a function is the collection of indices
    # in all participating subexpressions and the function itself
    # it is necessary to define this as such to enable reverse
    # McCormick constraint propagation
    relax_evaluator = m._working_problem._relaxed_evaluator
    has_subexpressions = length(evaluator.m.nlp_data.nlexpr) > 0
    dict_sparsity = Dict{Int64,Vector{Int64}}()
    if has_subexpressions
        for i = 1:length(evaluator.subexpressions)
            subexpr = evaluator.subexpressions[i]
            push!(relax_evaluator.subexpressions, NonlinearExpression!(subexpr, dict_sparsity, i,
                                                                      evaluator.subexpression_linearity,
                                                                      m._parameters.relax_tag))
        end
    end

    # scrubs udf functions using Cassette to remove odd data structures...
    # alternatively convert udfs to JuMP scripts...
    m._parameters.presolve_scrubber_flag && Script.scrub!(m._working_problem._nlp_data)
    if m._parameters.presolve_to_JuMP_flag
        Script.udf_loader!(m)
    end

    parameter_values = copy(evaluator.parameter_values)

    # add nonlinear objective
    if evaluator.has_nlobj
        m._working_problem._objective_nl = BufferedNonlinearFunction(evaluator.objective, MOI.NLPBoundsPair(-Inf, Inf),
                                                                     dict_sparsity, evaluator.subexpression_linearity,
                                                                     m._parameters.relax_tag)
    end

    # add nonlinear constraints
    constraint_bounds = m._working_problem._nlp_data.constraint_bounds
    for i = 1:length(evaluator.constraints)
        constraint = evaluator.constraints[i]
        bnds = constraint_bounds[i]
        push!(m._working_problem._nonlinear_constr, BufferedNonlinearFunction(constraint, bnds, dict_sparsity,
                                                                              evaluator.subexpression_linearity,
                                                                              m._parameters.relax_tag))
    end

    m._input_problem._nonlinear_count = length(m._working_problem._nonlinear_constr)
    m._working_problem._nonlinear_count = length(m._working_problem._nonlinear_constr)

    return nothing
end

function add_nonlinear_evaluator!(m::Optimizer)
    evaluator = m._input_problem._nlp_data.evaluator
    add_nonlinear_evaluator!(m, evaluator)
    return nothing
end

add_nonlinear_evaluator!(m::Optimizer, evaluator::Nothing) = nothing
add_nonlinear_evaluator!(m::Optimizer, evaluator::EmptyNLPEvaluator) = nothing
function add_nonlinear_evaluator!(m::Optimizer, evaluator::JuMP.NLPEvaluator)
    m._working_problem._relaxed_evaluator = Evaluator()

    relax_evaluator = m._working_problem._relaxed_evaluator
    relax_evaluator.variable_count = length(m._working_problem._variable_info)
    relax_evaluator.user_operators = evaluator.m.nlp_data.user_operators

    relax_evaluator.lower_variable_bounds = zeros(relax_evaluator.variable_count)
    relax_evaluator.upper_variable_bounds = zeros(relax_evaluator.variable_count)
    relax_evaluator.x                     = zeros(relax_evaluator.variable_count)
    relax_evaluator.num_mv_buffer         = zeros(relax_evaluator.variable_count)
    relax_evaluator.treat_x_as_number     = fill(false, relax_evaluator.variable_count)
    relax_evaluator.ctx                   = GuardCtx(metadata = GuardTracker(m._parameters.domain_violation_ϵ,
                                                                             m._parameters.domain_violation_guard_on))
    relax_evaluator.subgrad_tol           = m._parameters.subgrad_tol

    m._nonlinear_evaluator_created = true

    return nothing
end

function add_subexpression_buffers!(m::Optimizer)
    relax_evaluator = m._working_problem._relaxed_evaluator
    relax_evaluator.subexpressions_eval = fill(false, length(relax_evaluator.subexpressions))

    return nothing
end

"""
Translates input problem to working problem. Routines and checks and optional manipulation is left to the presolve stage.
"""
function initial_parse!(m::Optimizer)

    # reset initial time and solution statistics
    m._time_left = m._parameters.time_limit

    # add variables to working model
    ip = m._input_problem
    append!(m._working_problem._variable_info, ip._variable_info)
    m._working_problem._variable_count = ip._variable_count

    # add linear constraints to the working problem
    linear_leq = ip._linear_leq_constraints
    for i = 1:ip._linear_leq_count
        linear_func, leq_set = @inbounds linear_leq[i]
        push!(m._working_problem._saf_leq, AffineFunctionIneq(linear_func, leq_set))
        m._working_problem._saf_leq_count += 1
    end

    linear_geq = ip._linear_geq_constraints
    for i = 1:ip._linear_geq_count
        linear_func, geq_set = @inbounds linear_geq[i]
        push!(m._working_problem._saf_leq, AffineFunctionIneq(linear_func, geq_set))
        m._working_problem._saf_leq_count += 1
    end

    linear_eq = ip._linear_eq_constraints
    for i = 1:ip._linear_eq_count
        linear_func, eq_set = @inbounds linear_eq[i]
        push!(m._working_problem._saf_eq, AffineFunctionEq(linear_func, eq_set))
        m._working_problem._saf_eq_count += 1
    end

    # add quadratic constraints to the working problem
    quad_leq = ip._quadratic_leq_constraints
    for i = 1:ip._quadratic_leq_count
        quad_func, leq_set = @inbounds quad_leq[i]
        push!(m._working_problem._sqf_leq, BufferedQuadraticIneq(quad_func, leq_set))
        m._working_problem._sqf_leq_count += 1
    end

    quad_geq = ip._quadratic_geq_constraints
    for i = 1:ip._quadratic_geq_count
        quad_func, geq_set = @inbounds quad_geq[i]
        push!(m._working_problem._sqf_leq, BufferedQuadraticIneq(quad_func, geq_set))
        m._working_problem._sqf_leq_count += 1
    end

    quad_eq = ip._quadratic_eq_constraints
    for i = 1:ip._quadratic_eq_count
        quad_func, eq_set = @inbounds quad_eq[i]
        push!(m._working_problem._sqf_eq, BufferedQuadraticEq(quad_func, eq_set))
        m._working_problem._sqf_eq_count += 1
    end

    # add conic constraints to the working problem
    soc_vec = m._input_problem._conic_second_order
    for i = 1:ip._conic_second_order_count
        soc_func, soc_set = @inbounds soc_vec[i]
        first_variable_loc = soc_func.variables[1].value
        prior_lbnd = m._working_problem._variable_info[first_variable_loc].lower_bound
        m._working_problem._variable_info[first_variable_loc].lower_bound = max(prior_lbnd, 0.0)
        push!(m._working_problem._conic_second_order, BufferedSOC(soc_func, soc_set))
        m._working_problem._conic_second_order_count += 1
    end

    # set objective function
    m._working_problem._objective_type = ip._objective_type
    m._working_problem._objective_sv = ip._objective_sv
    m._working_problem._objective_saf = ip._objective_saf
    m._working_problem._objective_saf_parsed = AffineFunctionIneq(ip._objective_saf, LT_ZERO)
    m._working_problem._objective_sqf = BufferedQuadraticIneq(ip._objective_sqf, LT_ZERO)

    # add nonlinear constraints
    # the nonlinear evaluator loads with populated subexpressions which are then used
    # to asssess the linearity of subexpressions
    add_nonlinear_evaluator!(m)
    add_nonlinear_functions!(m)
    add_subexpression_buffers!(m)

    # converts a maximum problem to a minimum problem (internally) if necessary
    # this is placed after adding nonlinear functions as this prior routine
    # copies the nlp_block from the input_problem to the working problem
    _max_to_min!(m)
    reform_epigraph!(m)

    # labels the variable info and the _fixed_variable vector for each fixed variable
    label_fixed_variables!(m)

    # labels variables to branch on
    label_branch_variables!(m)

    # updates run and parse times
    new_time = time() - m._start_time
    m._parse_time = new_time
    m._run_time = new_time

    return nothing
end

function parse_classify_problem(::Val{LP}, m::Optimizer)

    input_problem = m._input_problem
    obj_typ = input_problem._objective_type

    is_lp = (obj_typ == SINGLE_VARIABLE) || (obj_typ == SCALAR_AFFINE)
    is_lp &= _second_order_cone_num(input_problem)       == 0
    is_lp &= _quadratic_num(input_problem)               == 0
    is_lp &= _nl_expr_num(input_problem)                 == 0
    is_lp &= _integer_variable_num(input_problem)        == 0

    if is_lp
        m._working_problem._problem_type = LP
        return true
    end

    return false
end

function parse_classify_problem(::Val{MILP}, m::Optimizer)
    input_problem = m._input_problem
    obj_typ = input_problem._objective_type

    is_milp = (obj_typ == SINGLE_VARIABLE) || (obj_typ == SCALAR_AFFINE)
    is_milp &= _second_order_cone_num(input_problem)       == 0
    is_milp &= _quadratic_num(input_problem)               == 0
    is_milp &= _nl_expr_num(input_problem)                 == 0
    is_milp &= _integer_variable_num(input_problem)        > 0

    if is_milp
        m._working_problem._problem_type = MILP
        return true
    end

    return false
end

function parse_classify_problem(::Val{SOCP}, m::Optimizer)
    input_problem = m._input_problem
    obj_typ = input_problem._objective_type

    is_socp = (obj_typ == SINGLE_VARIABLE) || (obj_typ == SCALAR_AFFINE)
    is_socp &= _second_order_cone_num(input_problem)       > 0
    is_socp &= _quadratic_num(input_problem)               == 0
    is_socp &= _nl_expr_num(input_problem)                 == 0
    is_socp &= _integer_variable_num(input_problem)        == 0

    if is_socp
        m._working_problem._problem_type = SOCP
        return true
    end

    return false
end

function parse_classify_problem(::Val{MISOCP}, m::Optimizer)
    input_problem = m._input_problem
    obj_typ = input_problem._objective_type

    is_misocp = (obj_typ == SINGLE_VARIABLE) || (obj_typ == SCALAR_AFFINE)
    is_misocp &= _second_order_cone_num(input_problem)       > 0
    is_misocp &= _quadratic_num(input_problem)               == 0
    is_misocp &= _nl_expr_num(input_problem)                 == 0
    is_misocp &= _integer_variable_num(input_problem)        > 0

    if is_misocp
        m._working_problem._problem_type = MISOCP
        return true
    end

    return false
end

#function parse_classify_problem(::typeof(SDP), m::Optimizer)
#    m._working_problem._problem_type = MINCVX
#    return false
#end

# is most general type support aka always true...
function parse_classify_problem(::Val{MINCVX}, m::Optimizer)
    m._working_problem._problem_type = MINCVX
    return true
end


"""
Classifies the problem type
"""
function parse_classify_problem!(m::Optimizer)

    # TODO: Transform Nonlinear Expressions into Linear or SOCP where appropriate
    parse_classify_problem(Val{LP}(), m)       && return
    parse_classify_problem(Val{MILP}(), m)     && return
    parse_classify_problem(Val{SOCP}(), m)     && return
    #parse_classify_problem(Val{SDP}(), m)      && return
    parse_classify_problem(Val{MISOCP}(), m)   && return
    parse_classify_problem(Val{MINCVX}(), m)   && return

    return
end

"""

Basic parsing for global solutions (no extensive manipulation)
"""
function parse_global!(t::ExtensionType, m::Optimizer)
    return nothing
end
