"""
$(TYPEDEF)

A structure used to expressions and problem descriptions EAGO uses to formulate
relaxed problems.
"""
Base.@kwdef mutable struct ParsedProblem

    # Problem classification (set in parse_classify_problem!)
    _problem_type::ProblemType = UNCLASSIFIED

    # objectives (set in initial_parse)
    _objective_sv::SV = SV(VI(-1))
    "_objective_saf stores the objective and is used for constructing linear affine cuts
     of any ObjectiveType"
    _objective_saf::SAF = SAF(SAT[], 0.0)
    _objective_saf_parsed::AffineFunctionIneq = AffineFunctionIneq()
    _objective_sqf::BufferedQuadraticIneq = BufferedQuadraticIneq()
    _objective_nl::BufferedNonlinearFunction = BufferedNonlinearFunction()
    _objective_type::ObjectiveType = UNSET

    # objective sense information (set by convert_to_min in parse.jl)
    _optimization_sense::MOI.OptimizationSense = MOI.MIN_SENSE

    # non-single variable constraints (set in initial_parse)
    _saf_leq::Vector{AffineFunctionIneq} = AffineFunctionIneq[]
    _saf_eq::Vector{AffineFunctionEq} = AffineFunctionEq[]
    _sqf_leq::Vector{BufferedQuadraticIneq} = BufferedQuadraticIneq[]
    _sqf_eq::Vector{BufferedQuadraticEq} = BufferedQuadraticEq[]
    _conic_second_order::Vector{BufferedSOC} = BufferedSOC[]

    # count of non-single variable constraints (set in initial_parse)
    _saf_leq_count::Int = 0
    _saf_eq_count::Int = 0
    _sqf_leq_count::Int = 0
    _sqf_eq_count::Int = 0
    _conic_second_order_count::Int = 0

    # nlp constraints (set in initial_parse)
    _nlp_data::MOI.NLPBlockData = empty_nlp_data()

    # storage for nonlinear functions
    _nonlinear_constr::Vector{BufferedNonlinearFunction} = BufferedNonlinearFunction[]

    # nonlinear constraint storage
    _nonlinear_count::Int = 0

    # nonlinear evaluator
    _relaxed_evaluator = Evaluator()

    # variables (set in initial_parse)
    _variable_info::Vector{VariableInfo} = VariableInfo[]
    _variable_count::Int = 0

    # count of single variable constraint types (set in load_relaxed_problem!)
    _var_leq_count::Int = 0
    _var_geq_count::Int = 0
    _var_eq_count::Int = 0
end

function Base.isempty(x::ParsedProblem)

    is_empty_flag = true

    new_input_problem = ParsedProblem()
    for field in fieldnames(ParsedProblem)

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

    is_empty_flag &= isempty(x._objective_saf.terms)
    is_empty_flag &= x._objective_saf.constant === 0.0
    is_empty_flag &= isempty(x._objective_saf_parsed.terms)
    is_empty_flag &= x._objective_saf_parsed.constant === 0.0
    is_empty_flag &= x._objective_saf_parsed.len === 0

    is_empty_flag &= isempty(x._objective_sqf.func.quadratic_terms)
    is_empty_flag &= isempty(x._objective_sqf.func.affine_terms)
    is_empty_flag &= x._objective_sqf.func.constant === 0.0

    is_empty_flag &= isempty(x._objective_sqf.buffer)
    is_empty_flag &= isempty(x._objective_sqf.saf.terms)
    is_empty_flag &= x._objective_sqf.saf.constant === 0.0
    is_empty_flag &= x._objective_sqf.len === 0

    return is_empty_flag
end