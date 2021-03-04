"""
$(TYPEDEF)

A structure used to hold objectives and constraints added to EAGO model.
The constraints generally aren't used for relaxations.
"""
Base.@kwdef mutable struct InputProblem

    # variables (set by MOI.add_variable in variables.jl)
    _variable_info::Vector{VariableInfo} = VariableInfo[]
    _variable_count::Int64 = 0

    # last constraint index added
    _last_constraint_index::Int = 0

    # linear constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _linear_leq_constraints::Vector{Tuple{SAF, LT}} = Tuple{SAF, LT}[]
    _linear_geq_constraints::Vector{Tuple{SAF, GT}} = Tuple{SAF, GT}[]
    _linear_eq_constraints::Vector{Tuple{SAF, ET}} = Tuple{SAF, ET}[]

    _linear_leq_count::Int = 0
    _linear_geq_count::Int = 0
    _linear_eq_count::Int = 0

    # quadratic constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _quadratic_leq_constraints::Vector{Tuple{SQF, LT}} = Tuple{SQF, LT}[]
    _quadratic_geq_constraints::Vector{Tuple{SQF, GT}} = Tuple{SQF, GT}[]
    _quadratic_eq_constraints::Vector{Tuple{SQF, ET}} = Tuple{SQF, ET}[]

    _quadratic_leq_count::Int = 0
    _quadratic_geq_count::Int = 0
    _quadratic_eq_count::Int = 0

    # conic constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _conic_second_order::Vector{Tuple{VECOFVAR, SECOND_ORDER_CONE}} = Tuple{VECOFVAR, SECOND_ORDER_CONE}[]
    _conic_power_corder::Vector{Tuple{VECOFVAR, POWER_CONE}} = Tuple{VECOFVAR, POWER_CONE}[]
    _conic_exp_order::Vector{Tuple{VECOFVAR, EXP_CONE}} = Tuple{VECOFVAR, EXP_CONE}[]

    _conic_second_order_count::Int = 0
    _conic_power_count::Int = 0
    _conic_exp_count::Int = 0

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

integer_variable_num(d::InputProblem) = count(is_integer.(d._variable_info))
function second_order_cone_num(d::InputProblem)
    d._conic_second_order_count
end
function quadratic_num(d::InputProblem)
    d._quadratic_leq_count + d._quadratic_geq_count + d._quadratic_eq_count
end
function nl_expr_num(d::InputProblem)
    nl_expr_number = d._objective_type === NONLINEAR ? 1 : 0
    nl_expr_number += d._nonlinear_count
    return nl_expr_number
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
