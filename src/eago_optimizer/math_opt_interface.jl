
#=
Directly bridging SV in GT -> LT triggers
=#
MOI.supports_constraint(::Optimizer,
                        ::Type{<:Union{SV, SAF, SQF}},
                        ::Type{<:Union{LT,GT,ET,IT}},
                        ) = true

# No automatic bridge from ZO to Integers exist. Closest is semi-integer to
# ZO + 2x SAF.
MOI.supports_constraint(::Optimizer,
                        ::Type{SV},
                        ::Type{<:Union{ZO, MOI.Integer}}) = true

MOI.supports_constraint(::Optimizer,
                        ::Type{<:Union{VECVAR}},
                        ::Type{<:Union{SOC_CONE, PSD_CONE}},
                        ) = true

# Using SingleVariable or ScalarAffineFunction as the objective greatly
# simplifies objective bound cuts and access to objective values. Using
# SAF prevents the introduce of a new auxilliary variable if the user
# defines a SAF objective type.
MOI.supports(::Optimizer,
             ::Union{MOI.ObjectiveSense,
                     MOI.ObjectiveFunction{SAF}}
                     ) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; copy_names = false)
    return MOIU.default_copy_to(model, src, copy_names)
end

MOI.get(m::Optimizer, ::MOI.SolverName) = "EAGO: Easy Advanced Global Optimization"
MOI.get(m::Optimizer, ::MOI.TerminationStatus) = m._termination_status_code
function MOI.get(m::Optimizer, v::MOI.PrimalStatus)
    if v.N > 1
        MOI.check_result_index_bounds(m, v)
    end
    m._primal_status_code
end
MOI.get(m::Optimizer, ::MOI.SolveTime) = m._run_time
MOI.get(m::Optimizer, ::MOI.NodeCount) = m._node_count
MOI.get(m::Optimizer, ::MOI.ResultCount) = (m._result_status_code === MOI.FEASIBLE_POINT) ? 1 : 0

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(m::Optimizer, ::MOI.Silent, b::Bool)
    m._is_silent = b
    return
end
MOI.get(m::Optimizer, ::MOI.Silent) = m._is_silent

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
MOI.get(m::Optimizer, ::MOI.TimeLimitSec) = isfinite(m._time_limit) ? m._time_limit : nothing
MOI.set(m::Optimizer, ::MOI.TimeLimitSec, ::Nothing)  = (m._time_limit = Inf;   return)
MOI.set(m::Optimizer, ::MOI.TimeLimitSec, v::Float64) = (m._time_limit = v; return)

function MOI.set(m::Optimizer, p::MOI.RawParameter, value)
    if !(p.name isa String) && !(p.name isa Symbol)
        error("EAGO only supports raw parameters with Symbol or String names.")
    end
    setfield!(m, Symbol(p.name), value)

    return nothing
end

for attr in (MOI.ListOfVariableIndices, MOI.ListOfConstraints)
    @eval function MOI.get(m::Optimizer, ::$attr)
        MOI.get(m._model, $attr())
    end
end

for attr in (MOI.ListOfConstraintAttributesSet, MOI.ListOfConstraintIndices)
    @eval function MOI.get(d::Optimizer, ::$attr{F,S}) where {F,S}
        return MOI.get(d._model, $attr{F,S}())
    end
end

for attr in (MOI.ConstraintFunction, MOI.ConstraintSet)
    @eval function MOI.get(d::Optimizer, ::$attr, ci::CI{SV,ZO})
        return MOI.get(d._model, $attr(), ci)
    end
    @eval function MOI.get(d::Optimizer, ::$attr, ci::CI{F,S}) where {F <: Union{SV, SAF, SQF}, S <: Union{LT,GT,ET,IT}}
        return MOI.get(d._model, $attr(), ci)
    end
    @eval function MOI.get(d::Optimizer, ::$attr, ci::CI{F,S}) where {F <: Union{VECVAR}, S <: Union{SOC_CONE, PSD_CONE}}
        return MOI.get(d._model, $attr(), ci)
    end
    @eval function MOI.set(d::Optimizer, ::$attr, ci::CI{SV,S}, v) where  S <: Union{ZO, MOI.Integer}
        return MOI.get(d._model, $attr(), ci, v)
    end
    @eval function MOI.set(d::Optimizer, ::$attr, ci::CI{F,S}) where {F <: Union{SV, SAF, SQF}, S <: Union{LT,GT,ET,IT}}
        return MOI.get(d._model, $attr(), ci, v)
    end
    @eval function MOI.set(d::Optimizer, ::$attr, ci::CI{F,S}) where {F <: Union{VECVAR}, S <: Union{SOC_CONE, PSD_CONE}}
        return MOI.get(d._model, $attr(), ci, v)
    end
end

function MOI.get(m::Optimizer, v::MOI.VariablePrimal, vi::MOI.VariableIndex)
    MOI.check_result_index_bounds(m, v)
    m._solution[vi.value]
end
MOI.get(m::Optimizer, p::MOI.VariablePrimal, vi::Vector{MOI.VariableIndex}) = MOI.get.(m, p, vi)
function MOI.get(m::Optimizer, v::MOI.ConstraintPrimal, ci::MOI.ConstraintIndex{SV, S}) where {S <: Union{IT, ZO, MOI.Integer}}
    MOI.check_result_index_bounds(m, v)
    return m._solution[ci.value]
end

function MOI.get(m::Optimizer, v::MOI.ConstraintPrimal, ci::MOI.ConstraintIndex{F, S}) where {F <: Union{SAF, SQF}, S <: Union{LT,GT,ET,IT}}
    MOI.check_result_index_bounds(m, v)
    return m._constraint_primal[ci]
end
function MOI.get(m::Optimizer, v::MOI.ConstraintPrimal, ci::MOI.ConstraintIndex{F, S}) where {F <: Union{VECVAR}, S <: Union{SOC_CONE, PSD_CONE}}
    MOI.check_result_index_bounds(m, v)
    return m._constraint_primal[ci]
end
MOI.get(opt::Optimizer, a::MOI.ConstraintPrimal, ci::Vector{MOI.ConstraintIndex}) = MOI.get.(opt, a, ci)

const EAGO_OPTIMIZER_ATTRIBUTES = Symbol[:relaxed_optimizer, :lp_optimizer, :mip_optimizer,
                                         :socp_optimizer, :semidefinite_optimizer,
                                         :nlp_optimizer, :minlp_optimizer,
                                         :relaxed_optimizer_kwargs, :upper_optimizer,
                                         :enable_optimize_hook, :ext, :ext_type, :_parameters]
const EAGO_MODEL_STRUCT_ATTRIBUTES = Symbol[:_stack, :_log, :_current_node, :_working_problem, :_model]
const EAGO_MODEL_NOT_STRUCT_ATTRIBUTES = setdiff(fieldnames(Optimizer), union(EAGO_OPTIMIZER_ATTRIBUTES,
                                                                              EAGO_MODEL_STRUCT_ATTRIBUTES))
const EAGO_MODEL_EITHER_ATTRIBUTE = union(EAGO_MODEL_STRUCT_ATTRIBUTES, EAGO_MODEL_NOT_STRUCT_ATTRIBUTES)

function MOI.empty!(m::Optimizer{T}) where T
    m._model = InputModel{T}()
    m._solver = GlobalOptimizer()
    m._problem_type = UNCLASSIFIED

    m._termination_status_code = MOI.OPTIMIZE_NOT_CALLED
    m._primal_status_code      = MOI.OTHER_RESULT_STATUS

    m._node_count      = 0
    m._objective_value = Inf
    m._objective_bound = -Inf
    m._relative_gap = Inf

    m._start_time   = 0.0
    m._time_limit   = Inf
    m._time_left    = 1000.0
    m._run_time     = 0.0
    m._parse_time   = 0.0
    return
end

function MOI.is_empty(m::Optimizer)

    is_empty_flag = MOI.is_empty(m._model)
    is_empty_flag &= MOI.is_empty(m._solver)
    is_empty_flag &= (m._problem_type == UNCLASSIFIED)

    is_empty_flag &= (m._termination_status_code == MOI.OPTIMIZE_NOT_CALLED)
    is_empty_flag &= (m._primal_status_code      == MOI.OTHER_RESULT_STATUS)

    is_empty_flag &= (m._node_count      == 0)
    is_empty_flag &= (m._objective_value == Inf)
    is_empty_flag &= (m._objective_bound == -Inf)
    is_empty_flag &= (m._relative_gap    == Inf)

    is_empty_flag &= (m._start_time == 0.0)
    is_empty_flag &= (m._time_limit == Inf)
    is_empty_flag &= (m._time_left  == 1000.0)
    is_empty_flag &= (m._run_time   == 0.0)
    is_empty_flag &= (m._parse_time == 0.0)

    return is_empty_flag
end

function MOI.get(m::Optimizer, v::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(m, v)
    m._objective_value
end
MOI.get(m::Optimizer, ::MOI.ObjectiveBound) = m._objective_bound
MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = MOI.get(m._model._input_model, MOI.NumberOfVariables())

function MOI.get(m::Optimizer, ::MOI.RelativeGap)
    LBD = m._objective_value
    UBD = m._objective_bound
    if m._model._optimization_sense === MOI.MAX_SENSE
        gap = abs(UBD - LBD)/abs(LBD)
    else
        gap = abs(UBD - LBD)/abs(UBD)
    end
    return gap
end

function MOI.get(m::Optimizer, p::MOI.RawParameter)
    if p.name isa String
        psym = Symbol(p.name)
    elseif p.name isa Symbol
        psym = p.name
    else
        error("EAGO only supports raw parameters with Symbol or String names.")
    end

    if psym in fieldnames(EAGOParameters)
        return getfield(m._parameters, psym)
    else
        return getfield(m, psym)
    end
end

MOI.add_variable(d::Optimizer) = MOI.add_variable(d._model._input_model)

function MOI.add_constraint(d::Optimizer, f::F, s::S) where {F <: Union{SV, SAF, SQF}, S <: Union{LT,GT,ET,IT}}
    MOI.add_constraint(d._model._input_model, f, s)
end
function MOI.add_constraint(d::Optimizer, f::SV, s::S) where {S<:Union{ZO, MOI.Integer}}
    MOI.add_constraint(d._model._input_model, f, s)
end
function MOI.add_constraint(d::Optimizer, f::F, s::S) where {F<:Union{VECVAR},
                                                             S<:Union{SOC_CONE, PSD_CONE}}
    append!(d._constraint_primal, zeros(MOI.dimension(s)))
    MOI.add_constraint(d._model._input_model, f, s)
end

function MOI.set(d::Optimizer, ::MOI.ObjectiveFunction{SAF}, func::SAF)
    MOI.set(d._model._input_model, MOI.ObjectiveFunction{SAF}(), func)
end

MOI.supports(d::Optimizer, ::MOI.NLPBlock) = true
function MOI.set(d::Optimizer, ::MOI.NLPBlock, nlp_data)
    d._model._nlp_data = nlp_data
    return
end
function MOI.set(d::Optimizer, ::MOI.NLPBlock, ::Nothing)
    d._model._nlp_data = nothing
    return
end

function MOI.set(d::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    MOI.set(d._model._input_model, MOI.ObjectiveSense(), sense)
end
