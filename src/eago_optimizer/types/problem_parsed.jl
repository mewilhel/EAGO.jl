@enum(ObjectiveType, UNSET, SCALAR_AFFINE, NONLINEAR)

"""
$(TYPEDEF)

A structure used to expressions and problem descriptions EAGO uses to formulate
relaxed problems.
"""
Base.@kwdef mutable struct ParsedProblem <: MOI.ModelLike

    # Problem classification (set in parse_classify_problem!)
    _problem_type::ProblemType = UNCLASSIFIED

    # objectives (set in initial_parse)
    _objective_saf::SAF = SAF(0.0, SAT(0.0, VI(-1), 0.0)
    _optimization_sense::MOI.OptimizationSense = MOI.MIN_SENSE

    # non-single variable constraints (set in initial_parse)
    _saf_leq::Vector{AffineFunctionIneq} = AffineFunctionIneq[]
    _saf_eq::Vector{AffineFunctionEq} = AffineFunctionEq[]
    _sqf_leq::Vector{BufferedQuadraticIneq} = BufferedQuadraticIneq[]
    _sqf_eq::Vector{BufferedQuadraticEq} = BufferedQuadraticEq[]
    _conic_second_order::Vector{BufferedSOC} = BufferedSOC[]

    # nlp constraints (set in initial_parse)
    _nlp_data::MOI.NLPBlockData = empty_nlp_data()

    # storage for nonlinear functions
    _nonlinear_constr::Vector{BufferedNonlinearFunction} = BufferedNonlinearFunction[]

    # nonlinear evaluator
    _relaxed_evaluator = Evaluator()

    # variables (set in initial_parse)
    _variable_info::Vector{VariableInfo} = VariableInfo[]
    _variable_num::Int = 0
end

function _add_constraint(wp::ParsedProblem, fs::Tuple{SAF,LT})
    push!(wp._saf_leq, AffineFunctionIneq(fs[1], fs[2])); return
end
function _add_constraint(wp::ParsedProblem, fs::Tuple{SAF,GT})
    push!(wp._saf_leq, AffineFunctionIneq(fs[1], fs[2])); return
end
function _add_constraint(wp::ParsedProblem, fs::Tuple{SAF,ET})
    push!(wp._saf_eq, AffineFunctionEq(fs[1], fs[2])); return
end

function _add_constraint(wp::ParsedProblem, fs::Tuple{SQF,LT})
    push!(wp._sqf_leq, BufferedQuadraticIneq(fs[1], fs[2])); return
end
function _add_constraint(wp::ParsedProblem, fs::Tuple{SQF,GT})
    push!(wp._sqf_leq, BufferedQuadraticIneq(fs[1], fs[2])); return
end
function _add_constraint(wp::ParsedProblem, fs::Tuple{SQF,ET})
    push!(wp._sqf_eq, BufferedQuadraticEq(fs[1], fs[2])); return
end

function _add_constraint(wp::ParsedProblem, fs::Tuple{VECVAR,SOC_CONE})
    f, s = fs
    first_variable_loc = f.variables[1].value
    prior_lbnd = wp._variable_info[first_variable_loc].lower_bound
    wp._variable_info[first_variable_loc].lower_bound = max(prior_lbnd, 0.0)
    push!(wp._conic_second_order, BufferedSOC(f, s))
    return
end

function _mod_decision_variables!(wp::ParsedProblem, d::Dict{S,T}) where {S,T}
    for (k,v) in d
        i = v[1].variable.value
        wp._variable_info[i] = VariableInfo(wp._variable_info[i], v[2])
    end
    return
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
    is_empty_flag &= x._objective_saf.constant == 0.0

    is_empty_flag &= isempty(x._objective_saf.terms)
    is_empty_flag &= x._objective_saf.constant == 0.0
    is_empty_flag &= isempty(x._objective_saf_parsed.func.terms)
    is_empty_flag &= x._objective_saf_parsed.func.constant == 0.0
    is_empty_flag &= x._objective_saf_parsed.len == 0

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
