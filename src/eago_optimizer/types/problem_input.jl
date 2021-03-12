"""
$(TYPEDEF)

A structure used to hold objectives and constraints added to EAGO model.
The constraints generally aren't used for relaxations.
"""
Base.@kwdef mutable struct InputProblem <: MOI.ModelLike

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
    _nlp_data::Union{MOI.NLPBlockData,Nothing} = nothing

    # objective sense information (set by MOI.set(m, ::ObjectiveSense...) in optimizer.jl)
    _optimization_sense::MOI.OptimizationSense = MOI.MIN_SENSE
end

@inline _variable_num(d::InputProblem) = d._variable_num
@inline _integer_variable_num(d::InputProblem) = count(is_integer.(d._variable_info))
@inline function _second_order_cone_num(d::InputProblem)
    length(d._conic_second_order_constraint)
end
@inline function _quadratic_num(d::InputProblem)
    num = length(d._quadratic_leq_constraint)
    num += length(d._quadratic_geq_constraint)
    num += length(d._quadratic_eq_constraint)
    return num
end
@inline function _nl_expr_num(d::InputProblem)
    num = d._objective_type === NONLINEAR ? 1 : 0
    num += d._nonlinear_count
    return num
end
@inline _optimization_sense(d::InputProblem) = d._optimization_sense
@inline _objective_type(d::InputProblem)     = d._objective_type

function _check_inbounds!(d::InputProblem, vi::VI)
    if !(1 <= vi.value <= d._variable_num)
        error("Invalid variable index $vi. ($(m._input_problem._variable_num) variables in the model.)")
    end
    return nothing
end
_check_inbounds!(d::InputProblem, var::SV) = _check_inbounds!(d, var.variable)
_check_inbounds!(d::InputProblem, aff::SAF) = foreach(x -> _check_inbounds!(d, x.variable_index), aff.terms)
function _check_inbounds!(d::InputProblem, quad::SQF)
    foreach(x -> _check_inbounds!(d, x.variable_index), quad.affine_terms)
    for term in quad.quadratic_terms
        _check_inbounds!(d, term.variable_index_1)
        _check_inbounds!(d, term.variable_index_2)
    end
    return nothing
end
_check_inbounds!(d::InputProblem, v::VECOFVAR) = foreach(x -> check_inbounds!(d, x), v.variables)

@inline _has_upper_bound(d::InputProblem, vi::MOI.VariableIndex) = d._variable_info[vi.value].has_upper_bound
@inline _has_lower_bound(d::InputProblem, vi::MOI.VariableIndex) = d._variable_info[vi.value].has_lower_bound
@inline _is_fixed(d::InputProblem, vi::MOI.VariableIndex) = d._variable_info[vi.value].is_fixed
@inline _is_integer(d::InputProblem, i::Int) = is_integer(d._variable_info[i])

function MOI.add_variable(d::InputProblem)
    d._variable_num += 1
    push!(d._variable_info, VariableInfo())
    return VI(d._variable_num)
end

function MOI.add_constraint(d::InputProblem, v::SV, lt::LT)
    vi = v.variable
    _check_inbounds!(d, vi)
    if isnan(lt.upper)
        error("Invalid upper bound value $(lt.upper).")
    elseif _has_upper_bound(d, vi)
        error("Upper bound on variable $vi already exists.")
    elseif _is_fixed(d, vi)
        error("Variable $vi is fixed. Cannot also set upper bound.")
    end
    d._variable_info[vi.value].upper_bound = lt.upper
    d._variable_info[vi.value].has_upper_bound = true
    return CI{SV, LT}(vi.value)
end

function MOI.add_constraint(d::InputProblem, v::SV, gt::GT)
    vi = v.variable
    _check_inbounds!(d, vi)
    if isnan(gt.lower)
        error("Invalid lower bound value $(gt.lower).")
    elseif _has_lower_bound(d, vi)
        error("Lower bound on variable $vi already exists.")
    elseif _is_fixed(d, vi)
        error("Variable $vi is fixed. Cannot also set lower bound.")
    end
    d._variable_info[vi.value].lower_bound = gt.lower
    d._variable_info[vi.value].has_lower_bound = true
    return CI{SV, GT}(vi.value)
end

function MOI.add_constraint(d::InputProblem, v::SV, eq::ET)
    vi = v.variable
    _check_inbounds!(d, vi)
    if isnan(eq.value)
        error("Invalid fixed value $(gt.lower).")
    elseif _has_lower_bound(d, vi)
        error("Variable $vi has a lower bound. Cannot be fixed.")
    elseif _has_upper_bound(d, vi)
        error("Variable $vi has an upper bound. Cannot be fixed.")
    elseif _is_fixed(d, vi)
        error("Variable $vi is already fixed.")
    end
    d._variable_info[vi.value].lower_bound = eq.value
    d._variable_info[vi.value].upper_bound = eq.value
    d._variable_info[vi.value].has_lower_bound = true
    d._variable_info[vi.value].has_upper_bound = true
    d._variable_info[vi.value].is_fixed = true
    return CI{SV, ET}(vi.value)
end

macro define_addconstraint(F, S, array_name)
    esc(quote
            function MOI.add_constraint(d::InputProblem, f::$F, s::$S)
                _check_inbounds!(d, f)
                d._constraint_index_num += 1
                push!(d.$(array_name), (copy(f), s))
                return CI{$F, $S}(d._constraint_index_num)
            end
            function ($array_name)(d::InputProblem)
                return d.$(array_name)
            end
        end)
end

macro define_addobjective(F, field_name)
    esc(quote
            function MOI.set(d::InputProblem, ::MOI.ObjectiveFunction{$F}, f::$F)
                _check_inbounds!(d, f)
                d.$field_name = f
                d._objective_type = _moi_to_obj_type(f)
                return nothing
            end
    end)
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

function MOI.set(d::InputProblem, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    if nlp_data.has_objective
        d._objective_type = NONLINEAR
    end
    d._nlp_data = nlp_data
    return
end
MOI.set(d::InputProblem, ::MOI.NLPBlock, ::Nothing) = nothing

function MOI.set(d::InputProblem, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    d._optimization_sense = sense
    return
end

function MOI.isempty(x::InputProblem)

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

    is_empty_flag &= x._nlp_data === nothing
    is_empty_flag &= isempty(x._objective_saf.terms)
    is_empty_flag &= x._objective_saf.constant === 0.0

    is_empty_flag &= isempty(x._objective_sqf.quadratic_terms)
    is_empty_flag &= isempty(x._objective_sqf.affine_terms)
    is_empty_flag &= x._objective_sqf.constant === 0.0

    return is_empty_flag
end
