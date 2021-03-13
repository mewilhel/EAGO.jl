
MOI.supports_constraint(::Optimizer,
                        ::Type{<:Union{SV, SAF, SQF}},
                        ::Type{<:Union{ET, GT, LT}},
                        ) = true
MOI.supports_constraint(::Optimizer, ::Type{SV}, ::Type{ZO}) = true

MOI.supports_constraint(::Optimizer,
                        ::Type{<:Union{VECOFVAR}},
                        ::Type{<:Union{SECOND_ORDER_CONE}},
                        ) = true

MOI.supports(::Optimizer,
             ::Union{MOI.ObjectiveSense,
                     MOI.ObjectiveFunction{SV},
                     MOI.ObjectiveFunction{SAF},
                     MOI.ObjectiveFunction{SQF}},
                     ) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; copy_names = false)
    return MOIU.default_copy_to(model, src, copy_names)
end

MOI.get(m::Optimizer, ::MOI.SolverName) = "EAGO: Easy Advanced Global Optimization"
MOI.get(m::Optimizer, ::MOI.TerminationStatus) = m._termination_status_code
MOI.get(m::Optimizer, ::MOI.PrimalStatus) = m._result_status_code
MOI.get(m::Optimizer, ::MOI.SolveTime) = m._run_time
MOI.get(m::Optimizer, ::MOI.NodeCount) = m._maximum_node_id
MOI.get(m::Optimizer, ::MOI.ResultCount) = (m._result_status_code === MOI.FEASIBLE_POINT) ? 1 : 0

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(m::Optimizer, ::MOI.Silent, value::Bool)
    if value
        m._parameters.verbosity = 0
        m._parameters.log_on = false
    end
    return
end
MOI.get(opt::Optimizer, ::MOI.Silent) = iszero(m._parameters.verbosity) && !m._parameters.log_on

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
MOI.get(m::Optimizer, ::MOI.TimeLimitSec) = isfinite(m._parameters.time_limit) ? m._parameters.time_limit : nothing
MOI.set(m::Optimizer, ::MOI.TimeLimitSec, ::Nothing)      = (m._parameters.time_limit = Inf;   return)
MOI.set(m::Optimizer, ::MOI.TimeLimitSec, value::Float64) = (m._parameters.time_limit = value; return)

function MOI.set(m::Optimizer, p::MOI.RawParameter, value)
    if p.name isa String
        psym = Symbol(p.name)
    elseif p.name isa Symbol
        psym = p.name
    else
        error("EAGO only supports raw parameters with Symbol or String names.")
    end

    if psym in fieldnames(EAGOParameters)
        setfield!(m._parameters, psym, value)
    else
        setfield!(m, psym, value)
    end

    return nothing
end

function MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices)
    MOI.get(m._input_problem, MOI.ListOfVariableIndices())
end

MOI.get(m::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex) = m._continuous_solution[vi.value]
MOI.get(m::Optimizer, p::MOI.VariablePrimal, vi::Vector{MOI.VariableIndex}) = MOI.get.(m, p, vi)

const EAGO_OPTIMIZER_ATTRIBUTES = Symbol[:relaxed_optimizer, :lp_optimizer, :mip_optimizer,
                                         :socp_optimizer, :semidefinite_optimizer,
                                         :nlp_optimizer, :minlp_optimizer,
                                         :relaxed_optimizer_kwargs, :upper_optimizer,
                                         :enable_optimize_hook, :ext, :ext_type, :_parameters]
const EAGO_MODEL_STRUCT_ATTRIBUTES = Symbol[:_stack, :_log, :_current_node, :_working_problem, :_input_problem]
const EAGO_MODEL_NOT_STRUCT_ATTRIBUTES = setdiff(fieldnames(Optimizer), union(EAGO_OPTIMIZER_ATTRIBUTES,
                                                                              EAGO_MODEL_STRUCT_ATTRIBUTES))
const EAGO_MODEL_EITHER_ATTRIBUTE = union(EAGO_MODEL_STRUCT_ATTRIBUTES, EAGO_MODEL_NOT_STRUCT_ATTRIBUTES)

function MOI.empty!(m::Optimizer)
    new_optimizer = Optimizer()
    for field in EAGO_MODEL_EITHER_ATTRIBUTE
        setfield!(m, field, getfield(new_optimizer, field))
    end
    return
end

function MOI.is_empty(m::Optimizer)

    is_empty_flag = uninitialized(m._current_node)
    is_empty_flag &= isempty(m._stack)
    is_empty_flag &= isempty(m._log)
    is_empty_flag &= isempty(m._input_problem)
    is_empty_flag &= isempty(m._working_problem)

    new_optimizer = Optimizer()
    for field in EAGO_MODEL_NOT_STRUCT_ATTRIBUTES
        if getfield(m, field) != getfield(new_optimizer, field)
            is_empty_flag = false
            break
        end
    end
    return is_empty_flag
end

MOI.get(m::Optimizer, ::MOI.ObjectiveValue) = m._objective_value
MOI.get(m::Optimizer, ::MOI.ObjectiveBound) = m._objective_bound
MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = m._input_problem._variable_num

function MOI.get(m::Optimizer, ::MOI.RelativeGap)
    LBD = m._objective_value
    UBD = m._objective_bound
    if m._input_problem._optimization_sense === MOI.MAX_SENSE
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

MOI.add_variable(d::Optimizer) = MOI.add_variable(d._input_problem)

function MOI.add_constraint(d::Optimizer, f::F, s::S) where {F<:Union{SV, SAF, SQF},
                                                             S<:Union{ET, GT, LT}}
    MOI.add_constraint(d._input_problem, f, s)
end
MOI.add_constraint(d::Optimizer, f::SV, s::ZO) = MOI.add_constraint(d._input_problem, f, s)

function MOI.add_constraint(d::Optimizer, f::F, s::S) where {F<:Union{VECOFVAR},
                                                             S<:Union{SECOND_ORDER_CONE}}
    MOI.add_constraint(d._input_problem, f, s)
end

function MOI.set(d::Optimizer, ::MOI.ObjectiveFunction{F}, func::F) where F<:Union{SV, SAF, SQF}
    MOI.set(d._input_problem, MOI.ObjectiveFunction{F}(), func)
end

MOI.supports(d::Optimizer, ::MOI.NLPBlock) = true
function MOI.set(d::Optimizer, ::MOI.NLPBlock, nlp_data)
    MOI.set(d._input_problem, MOI.NLPBlock(), nlp_data)
end
function MOI.set(d::Optimizer, ::MOI.NLPBlock, ::Nothing)
    MOI.set(d._input_problem, MOI.NLPBlock(), nothing)
end

function MOI.set(d::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    MOI.set(d._input_problem, MOI.ObjectiveSense(), sense)
end
