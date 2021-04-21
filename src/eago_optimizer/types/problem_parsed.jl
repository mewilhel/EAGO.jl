


# Holds a copy of the problem passed to EAGO
# Needs separate storage for the MOI compilant model and the nlp_data
Base.@kwdef mutable struct InputModel{T<:Real}
    _input_model::MOIU.Model{T} = MOIU.Model{T}()
    _nlp_data::Union{MOI.NLPBlockData, Nothing} = nothing
end

MOI.is_empty(m::InputModel) = (m._nlp_data == nothing) && MOI.is_empty(m._input_model)

@enum(ObjectiveType, UNSET, SCALAR_AFFINE, NONLINEAR)

"""
$(TYPEDEF)

A structure used to expressions and problem descriptions EAGO uses to formulate
relaxed problems.
"""
Base.@kwdef mutable struct ParsedProblem{N,T} <: MOI.ModelLike

    # Problem classification (set in parse_classify_problem!)
    _problem_type::ProblemType = UNCLASSIFIED

    # objectives (set in initial_parse)
    _objective::SAF = SAF(SAT{T}[], zero(T))
    _objective_parsed::AffineFunction{T} = AffineFunction{T}()
    _optimization_sense::MOI.OptimizationSense = MOI.MIN_SENSE

    # non-single variable constraints (set in initial_parse)
    _saf_leq::Vector{AffineFunction{T}} = AffineFunction{T}[]
    _saf_eq::Vector{AffineFunction{T}} = AffineFunction{T}[]
    _sqf_leq::Vector{BufferedQuadraticIneq{T}} = BufferedQuadraticIneq{T}[]
    _sqf_eq::Vector{BufferedQuadraticEq{T}} = BufferedQuadraticEq{T}[]
    _conic_second_order::Vector{BufferedSOC{T}} = BufferedSOC{T}[]

    # nlp constraints (set in initial_parse)
    _nlp_data::Union{MOI.NLPBlockData,Nothing} = nothing

    # storage for nonlinear functions
    _nonlinear_constr::Vector{BufferedNonlinearFunction{T}} = BufferedNonlinearFunction{T}[]

    # nonlinear evaluator
    _relaxed_evaluator = Evaluator{N,T}()

    # variables (set in initial_parse)
    _variable_info::Vector{VariableInfo} = VariableInfo[]
    _variable_num::Int = 0
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

    is_empty_flag &= isempty(x._objective.terms)
    is_empty_flag &= x._objective.constant == 0.0

    is_empty_flag &= isempty(x._objective.terms)
    is_empty_flag &= x._objective.constant == 0.0
    is_empty_flag &= isempty(x._objective_parsed.func.terms)
    is_empty_flag &= x._objective_parsed.func.constant == 0.0
    is_empty_flag &= x._objective_parsed.len == 0

    return is_empty_flag
end

#
# `GeometricMeanCone` is `SecondOrderCone` representable
# `RootDetConeTriangle` is representable by a `PositiveSemidefiniteConeTriangle` and an `GeometricMeanCone`
# `LogDetConeTriangle` is representable by a `PositiveSemidefiniteConeTriangle` and `ExponentialCone`
# Transforms a `G`-in-`GreaterThan{T}` constraint into an `F`-in-`LessThan{T}` drop support for GT
# IndicatorSOS1Bridge not need for
# The `NormSpectralCone` is representable with a PSD constraint, since
# ``t \\ge \\sigma_1(X)`` if and only if ``[tI X^\\top; X tI] \\succ 0``.
# The `RelativeEntropyCone` is representable with exponential cone and LP constraints
# The `RotatedSecondOrderCone` is `SecondOrderCone` representable
#=
The `geometric mean cone` is representable with a relative entropy constraint and a
nonnegative auxiliary variable
=#

#=
The `SemiToBinaryBridge` replaces an Semicontinuous constraint:
``x \\in \\mathsf{Semicontinuous}(l, u)``
is replaced by:
``z \\in \\{0, 1\\}``,
``x \\leq z \\cdot u ``,
``x \\geq z \\cdot l ``.
The `SemiToBinaryBridge` replaces an Semiinteger constraint:
``x \\in Semiinteger(l, u)``
is replaced by:
``z \\in \\{0, 1\\}``,
``x \\in \\mathbb{Z}``,
``x \\leq z \\cdot u ``,
``x \\geq z \\cdot l ``.
=#

#=
No need for Semicontinuous or Semiinteger wp representation... but may speed
up mixed integer solves.
=#

# LP relax bridges ->
# The `NormInfinityCone` is representable with LP constraints
# The `NormOneCone` is representable with LP constraint

# SlackBridge{T, F, G} to perform epigraph rearrangement -> FunctionizeBridge{T}
# only support scalar_affine_objectives
