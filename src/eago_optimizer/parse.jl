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
#=
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

function _check_set_is_fixed(v::VariableInfo)
    v.is_fixed && return true
    v.is_fixed = x.lower_bound === x.upper_bound
    return v.is_fixed
end
=#

"""
$(TYPEDSIGNATURES)

Detects any variables set to a fixed value by equality or inequality constraints
and populates the `_fixed_variable` storage array.
"""
function label_fixed_variables!(m::Optimizer)
    map!(x -> _check_set_is_fixed(x), m._fixed_variable, m._working_problem._variable_info)
end

function _label_branch_quad!(d, sqf)
    for term in sqf.func.quadratic_terms
        d[term.variable_index_1.value] = true
        d[term.variable_index_2.value] = true
    end
    return
end

function _label_branch_nl!(d, f)
    for i in f.expr.grad_sparsity
        m._branch_variables[i] = true
    end
    return
end

"""
$(TYPEDSIGNATURES)

Detects any variables participating in nonconvex terms and populates the
`_branch_variables` storage array.
"""
function label_branch_variables!(m::Optimizer)

    wp = _working_problem(m)
    m._user_branch_variables = !isempty(m.branch_variable)

    if m._user_branch_variables
        append!(m._branch_variables, m.branch_variable)
    else
        append!(m._branch_variables, fill(false, wp._variable_num))

        # adds nonlinear terms in quadratic constraints
        foreach(x -> _label_branch_quad!(m._branch_variables, x), wp._sqf_leq)
        foreach(x -> _label_branch_quad!(m._branch_variables, x), wp._sqf_eq)

        # TODO: label branch variable in conic forms

        # label nonlinear branch variables (assumes affine terms have been extracted)
        foreach(x -> _label_branch_quad!(m._branch_variables, x), wp._sqf_leq)
        foreach(x -> _label_branch_nl!(m._branch_variables, x), wp._nonlinear_constr)
    end

    # add a map of branch/node index to variables in the continuous solution
    for i = 1:wp._variable_num
        if wp._variable_info[i].is_fixed
            m._branch_variables[i] = false
            continue
        end
        if m._branch_variables[i]
            push!(m._branch_to_sol_map, i)
        end
    end

    # creates reverse map
    m._sol_to_branch_map = zeros(wp._variable_num)
    for i = 1:length(m._branch_to_sol_map)
        j = m._branch_to_sol_map[i]
        m._sol_to_branch_map[j] = i
    end

    # adds branch solution to branch map to evaluator
    wp._relaxed_evaluator.node_to_variable_map = m._branch_to_sol_map
    wp._relaxed_evaluator.variable_to_node_map = m._sol_to_branch_map
    wp._relaxed_evaluator.node_count = length(m._branch_to_sol_map)

    return nothing
end


add_nonlinear_functions!(m::Optimizer) = add_nonlinear_functions!(m, m._input_problem._nlp_data)
add_nonlinear_functions!(m::Optimizer, b::Nothing) = nothing
add_nonlinear_functions!(m::Optimizer, b::MOI.NLPBlockData) = add_nonlinear_functions!(m, b.evaluator)
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
                                                                      m.relax_tag))
        end
    end

    # scrubs udf functions using Cassette to remove odd data structures...
    # alternatively convert udfs to JuMP scripts...
    m.presolve_scrubber_flag && Script.scrub!(m._working_problem._nlp_data)
    if m.presolve_to_JuMP_flag
        Script.udf_loader!(m)
    end

    parameter_values = copy(evaluator.parameter_values)

    # add nonlinear objective
    if evaluator.has_nlobj
        # TODO: Should add a constraint to nlp_data... or evaluator not objective
        m._working_problem._objective_nl = BufferedNonlinearFunction(evaluator.objective, MOI.NLPBoundsPair(-Inf, Inf),
                                                                     dict_sparsity, evaluator.subexpression_linearity,
                                                                     m.relax_tag)
    end

    # add nonlinear constraints
    constraint_bounds = m._working_problem._nlp_data.constraint_bounds
    for i = 1:length(evaluator.constraints)
        constraint = evaluator.constraints[i]
        bnds = constraint_bounds[i]
        push!(m._working_problem._nonlinear_constr, BufferedNonlinearFunction(constraint, bnds, dict_sparsity,
                                                                              evaluator.subexpression_linearity,
                                                                              m.relax_tag))
    end

    m._input_problem._nonlinear_count = length(m._working_problem._nonlinear_constr)
    m._working_problem._nonlinear_count = length(m._working_problem._nonlinear_constr)

    return nothing
end

add_nonlinear_evaluator!(m::Optimizer) = add_nonlinear_evaluator!(m, m._input_problem._nlp_data)
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
    relax_evaluator.ctx                   = GuardCtx(metadata = GuardTracker(m.domain_violation_ϵ,
                                                                             m.domain_violation_guard_on))
    relax_evaluator.subgrad_tol           = m.subgrad_tol

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

    # updates run and parse times
    new_time = time() - m._start_time
    m._parse_time = new_time
    m._run_time = new_time

    return nothing
end

const INTEGER_FUNC_SET{T <: Real} = Union{Tuple{MOI.SingleVariable, MOI.ZeroOne},
                                          Tuple{MOI.SingleVariable, MOI.Integer},
                                          Tuple{MOI.SingleVariable, MOI.Semicontinuous},
                                          Tuple{MOI.SingleVariable, MOI.Semiinteger},
                                          Tuple{MOI.VectorOfVariables, MOI.SOS1{T}},
                                          Tuple{MOI.VectorOfVariables, MOI.SOS2{T}}}

const LINEAR_FUNC_SET{T <: Real} = Union{Tuple{MOI.SingleVariable, MOI.LessThan{T}},
                                         Tuple{MOI.SingleVariable, MOI.GreaterThan{T}},
                                         Tuple{MOI.SingleVariable, MOI.EqualTo{T}},
                                         Tuple{MOI.SingleVariable, MOI.Interval{T}},
                                         Tuple{MOI.ScalarAffineFunction, MOI.LessThan{T}},
                                         Tuple{MOI.ScalarAffineFunction, MOI.GreaterThan{T}},
                                         Tuple{MOI.ScalarAffineFunction, MOI.EqualTo{T}},
                                         Tuple{MOI.ScalarAffineFunction, MOI.Interval{T}},
                                         Tuple{MOI.VectorOfVariables,    MOI.Nonnegatives},
                                         Tuple{MOI.VectorOfVariables,    MOI.Nonpositives},
                                         Tuple{MOI.VectorOfVariables,    MOI.Zeros},
                                         Tuple{MOI.VectorAffineFunction, MOI.Nonnegatives},
                                         Tuple{MOI.VectorAffineFunction, MOI.Nonpositives},
                                         Tuple{MOI.VectorAffineFunction, MOI.Zeros}}

const CONE_FUNC_SET =  Union{Tuple{MOI.VectorOfVariables, MOI.NormInfinityCone},
                             Tuple{MOI.VectorOfVariables, MOI.NormOneCone},
                             Tuple{MOI.VectorOfVariables, MOI.SecondOrderCone},
                             Tuple{MOI.VectorOfVariables, MOI.RotatedSecondOrderCone},
                             Tuple{MOI.VectorOfVariables, MOI.GeometricMeanCone},
                             Tuple{MOI.VectorOfVariables, MOI.ExponentialCone},
                             Tuple{MOI.VectorOfVariables, MOI.DualExponentialCone},
                             Tuple{MOI.VectorOfVariables, MOI.RelativeEntropyCone},
                             Tuple{MOI.VectorOfVariables, MOI.NormSpectralCone},
                             Tuple{MOI.VectorOfVariables, MOI.NormNuclearCone}}

const SDP_FUNC_SET = Union{Tuple{MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle},
                           Tuple{MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeSquare},
                           Tuple{MOI.VectorOfVariables, MOI.RootDetConeTriangle},
                           Tuple{MOI.VectorOfVariables, MOI.RootDetConeSquare},
                           Tuple{MOI.VectorOfVariables, MOI.LogDetConeTriangle},
                           Tuple{MOI.VectorOfVariables, MOI.LogDetConeSquare}}

_in_linear(S, T) = Tuple{S, T} <: LINEAR_FUNC_SET{Float64}
_in_milp(S, T) = _in_linear(S, T) | (Tuple{S, T} <: INTEGER_FUNC_SET{Float64})
_in_socp(S, T) = _in_linear(S, T) | (Tuple{S, T} <: CONE_FUNC_SET)
_in_sdp(S, T) = _in_socp(S, T) | (Tuple{S, T} <: SDP_FUNC_SET)

function _parse_classify_problem(::Val{LP}, m::Optimizer)
    ip = _input_model(m)
    if MOI.get(ip, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE
        is_lp = true
    else
        is_lp = MOI.get(ip, MOI.ObjectiveFunctionType()) == MOI.SingleVariable
        is_lp |= MOI.get(ip, MOI.ObjectiveFunctionType()) == MOI.ScalarAffineFunction{Float64}
    end
    is_lp &= isempty(filter(x -> !_in_linear(x[1], x[2]), MOI.get(ip, MOI.ListOfConstraints())))
    is_lp &= _input_nlp_data(m) === nothing
    if is_lp
        m._problem_type = LP
    end
    return is_lp
end

function _parse_classify_problem(::Val{MILP}, m::Optimizer)
    ip = _input_model(m)
    if MOI.get(ip, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE
        is_milp = true
    else
        is_milp = MOI.get(ip, MOI.ObjectiveFunctionType()) == MOI.SingleVariable
        is_milp |= MOI.get(ip, MOI.ObjectiveFunctionType()) == MOI.ScalarAffineFunction{Float64}
    end
    is_milp &= isempty(filter(x -> !_in_milp(x[1], x[2]), MOI.get(ip, MOI.ListOfConstraints())))
    is_milp &= _input_nlp_data(m) === nothing
    if is_milp
        m._problem_type = MILP
    end
    return is_milp
end

function _parse_classify_problem(::Val{SOCP}, m::Optimizer)
    ip = _input_model(m)
    if MOI.get(ip, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE
        is_socp = true
    else
        is_socp = MOI.get(ip, MOI.ObjectiveFunctionType()) == MOI.SingleVariable
        is_socp |= MOI.get(ip, MOI.ObjectiveFunctionType()) == MOI.ScalarAffineFunction{Float64}
    end
    is_socp &= isempty(filter(x -> !_in_socp(x[1], x[2]), MOI.get(ip, MOI.ListOfConstraints())))
    is_socp &= _input_nlp_data(m) === nothing
    if is_socp
        m._problem_type = SOCP
    end
    return is_socp
end

function _parse_classify_problem(::Val{SDP}, m::Optimizer)
    ip = _input_model(m)
    if MOI.get(ip, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE
        is_sdp = true
    else
        is_sdp = MOI.get(ip, MOI.ObjectiveFunctionType()) == MOI.SingleVariable
        is_sdp |= MOI.get(ip, MOI.ObjectiveFunctionType()) == MOI.ScalarAffineFunction{Float64}
    end
    is_sdp &= isempty(filter(x -> !_in_sdp(x[1], x[2]), MOI.get(ip, MOI.ListOfConstraints())))
    is_sdp &= _input_nlp_data(m) === nothing
    if is_sdp
        m._problem_type = SDP
    end
    return is_sdp
end

#=
function _parse_classify_problem(::Val{MISOCP}, ip::InputProblem, wp::ParsedProblem)
    is_misocp =  _objective_type(ip) == SINGLE_VARIABLE
    is_misocp |= _objective_type(ip) == SCALAR_AFFINE
    is_misocp &= _second_order_cone_num(ip)       > 0
    is_misocp &= _quadratic_num(ip)               == 0
    is_misocp &= _nl_expr_num(ip)                 == 0
    is_misocp &= _integer_variable_num(ip)        > 0
    if is_misocp
        wp._problem_type = MISOCP
    end
    return is_misocp
end
=#

function _parse_classify_problem(::Val{MINCVX}, m::Optimizer)
    m._problem_type = MINCVX
    return true
end


"""
Classifies the problem type
"""
function _parse_classify_problem!(m::Optimizer)
    _parse_classify_problem(Val{LP}(), m)       && return
    _parse_classify_problem(Val{MILP}(), m)     && return
    _parse_classify_problem(Val{SOCP}(), m)     && return
    _parse_classify_problem(Val{SDP}(), m)      && return
    #_parse_classify_problem(Val{MISOCP}(), m)   && return
    _parse_classify_problem(Val{MINCVX}(), m)   && return
    return
end

"""

Basic parsing for global solutions (no extensive manipulation)
"""
function parse_global!(t::ExtensionType, m::Optimizer)
    return nothing
end
