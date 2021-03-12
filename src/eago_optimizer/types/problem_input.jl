"""
$(TYPEDEF)

A structure used to hold objectives and constraints added to EAGO model.
The constraints generally aren't used for relaxations.
"""
Base.@kwdef mutable struct InputProblem

    # variables (set by MOI.add_variable in variables.jl)
    _variable_info::Vector{VariableInfo} = VariableInfo[]
    _variable_num::Int64 = 0

    # last constraint index added
    _constraint_index_num::Int = 0

    # linear constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _linear_leq_constraint::Vector{Tuple{SAF, LT}} = Tuple{SAF, LT}[]
    _linear_geq_constraint::Vector{Tuple{SAF, GT}} = Tuple{SAF, GT}[]
    _linear_eq_constraint::Vector{Tuple{SAF, ET}} = Tuple{SAF, ET}[]

    # quadratic constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _quadratic_leq_constraint::Vector{Tuple{SQF, LT}} = Tuple{SQF, LT}[]
    _quadratic_geq_constraint::Vector{Tuple{SQF, GT}} = Tuple{SQF, GT}[]
    _quadratic_eq_constraint::Vector{Tuple{SQF, ET}} = Tuple{SQF, ET}[]

    # conic constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _conic_second_order_constraint::Vector{Tuple{VECOFVAR, SECOND_ORDER_CONE}} = Tuple{VECOFVAR, SECOND_ORDER_CONE}[]

    # nonlinear constraint storage
    _nonlinear_count::Int = 0

    # objective information (set by MOI.set(m, ::ObjectiveFunction...) in optimizer.jl)
    _objective_sv::SV = SV(VI(-1))
    _objective_saf::SAF = SAF(SAT[], 0.0)
    _objective_sqf::SQF = SQF(SAT[], SQT[], 0.0)
    _objective_type::ObjectiveType = UNSET

    # nlp constraints (set by MOI.set(m, ::NLPBlockData...) in optimizer.jl)
    _nlp_data::MOI.NLPBlockData = empty_nlp_data()

    # objective sense information (set by MOI.set(m, ::ObjectiveSense...) in optimizer.jl)
    _optimization_sense::MOI.OptimizationSense = MOI.MIN_SENSE
end

@inline _integer_variable_num(d::InputProblem) = count(is_integer.(d._variable_info))
@inline function _second_order_cone_num(d::InputProblem)
    d._conic_second_order_count
end
@inline function _quadratic_num(d::InputProblem)
    d._quadratic_leq_count + d._quadratic_geq_count + d._quadratic_eq_count
end
@inline function _nl_expr_num(d::InputProblem)
    nl_expr_number = d._objective_type === NONLINEAR ? 1 : 0
    nl_expr_number += d._nonlinear_count
    return nl_expr_number
end
@inline _optimization_sense(d::InputProblem) = d._optimization_sense
@inline _objective_type(d::InputProblem)     = d._objective_type

function _add_var(d::InputProblem)
    d._variable_num += 1
    push!(d._variable_info, VariableInfo())
    return VI(d._variable_num)
end

function _add_var_constraint!(d::InputProblem, f::SV, s::LT)
    vi = v.variable
    d._variable_info[vi.value].upper_bound = lt.upper
    d._variable_info[vi.value].has_upper_bound = true
    return CI{SV, LT}(vi.value)
end

function add_var_constraint!(d::InputProblem, v::SV, gt::GT)
    vi = v.variable
    d._variable_info[vi.value].lower_bound = gt.lower
    d._variable_info[vi.value].has_lower_bound = true
    return CI{SV, GT}(vi.value)
end

function _add_var_constraint!(d::InputProblem, v::SV, eq::ET)
    vi = v.variable
    d._variable_info[vi.value].lower_bound = eq.value
    d._variable_info[vi.value].upper_bound = eq.value
    d._variable_info[vi.value].has_lower_bound = true
    d._variable_info[vi.value].has_upper_bound = true
    d._variable_info[vi.value].is_fixed = true
    return CI{SV, ET}(vi.value)
end

macro define_add_variable_constraint(F, S, array_name)
    quote
        function _add_variable_constraint!(d::InputProblem, f::$F, s::$S)
            push!(d.$(array_name), (copy(f), s))
            return CI{F, S}(d._constraint_index_num)
        end
    end
end

macro define_addconstraint(F, S, array_name)
    quote
        function _add_constraint!(d::InputProblem, f::$F, s::$S)
            push!(d.$(array_name), (copy(f), s))
            return CI{F, S}(d._constraint_index_num)
        end
    end
end

macro define_addobjective(F, field_name)
    quote
        function _setobjective!(d::InputProblem, func::F)
            d.$field_name = func
            d._objective_type = _moi_to_obj_type(func)
            return nothing
        end
    end
end

@define_addconstraint SAF LT _linear_leq_constraint
@define_addconstraint SAF GT _linear_geq_constraint
@define_addconstraint SAF ET _linear_eq_constraint

@define_addconstraint SQF LT _quadratic_leq_constraint
@define_addconstraint SQF GT _quadratic_geq_constraint
@define_addconstraint SQF ET _quadratic_eq_constraint
@define_addconstraint VECOFVAR SECOND_ORDER_CONE _conic_second_order_constraint

_moi_to_obj_type(d::SV) = SINGLE_VARIABLE
_moi_to_obj_type(d::SAF) = SCALAR_AFFINE
_moi_to_obj_type(d::SQF) = SCALAR_QUADRATIC

@define_addobjective SV  _objective_sv
@define_addobjective SAF _objective_saf
@define_addobjective SQF _objective_sqf

function _set!(d::InputProblem, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    if nlp_data.has_objective
        d._objective_type = NONLINEAR
    end
    d._nlp_data = nlp_data
    return
end

function _set!(d::InputProblem, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    d._optimization_sense = sense
    return
end

function Base.isempty(x::InputProblem)

    is_empty_flag = true
    new_input_problem = InputProblem()

    for field in fieldnames(InputProblem)

        field_value = getfield(x, field)

        if field_value isa Array
            if !isempty(field_value)
                is_empty_flag = false
                break
            end

        elseif field_value isa Number
            if getfield(new_input_problem, field) !== field_value
                is_empty_flag = false
                break
            end

        end
    end

    is_empty_flag &= x._nlp_data.evaluator isa EmptyNLPEvaluator
    is_empty_flag &= !x._nlp_data.has_objective
    is_empty_flag &= isempty(x._nlp_data.constraint_bounds)

    is_empty_flag &= isempty(x._objective_saf.terms)
    is_empty_flag &= x._objective_saf.constant === 0.0

    is_empty_flag &= isempty(x._objective_sqf.quadratic_terms)
    is_empty_flag &= isempty(x._objective_sqf.affine_terms)
    is_empty_flag &= x._objective_sqf.constant === 0.0

    return is_empty_flag
end
