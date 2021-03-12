
MOI.supports_constraint(::Optimizer,
                        ::Type{<:Union{SV, SAF, SQF}},
                        ::Type{<:Union{ET, GT, LT}},
                        ) = true

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

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    @assert !copy_names
    input_prob = dest._input_problem
    idx_map = MOI.Utilities.IndexMap()

    # add variables to index map & input problem
    list_of_vi = MOI.ListOfVariableIndices()
    foreach(i -> setindex!(idx_map, i, _add_var(input_prob)), list_of_vi)

    get_ci_list(F, S) = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    get_con_fun(ci) = MOI.get(src, MOI.ConstraintFunction(), ci)
    get_con_set(ci) = MOI.get(src, MOI.ConstraintSet(), ci)
    function _populate!(F, S, imap)
        for ci in get_ci_list(F, S)
            input_prob._constraint_index_num += 1
            imap[ci] = CI{F, S}(input_prob._constraint_index_num)
            _add_constraint!(input_prob, get_con_fun(ci), get_con_set(ci))
        end
    end
    function _populate_var!(F, S, imap)
        for ci in get_ci_list(F, S)
            imap[ci] = _add_var_constraint!(input_prob, get_con_fun(ci), get_con_set(ci))
        end
    end

    _populate_var!(SV, LT, idx_map)
    _populate_var!(SV, GT, idx_map)
    _populate_var!(SV, ET, idx_map)

    _populate!(SAF, LT, idx_map)
    _populate!(SAF, GT, idx_map)
    _populate!(SAF, ET, idx_map)

    _populate!(SQF, LT, idx_map)
    _populate!(SQF, GT, idx_map)
    _populate!(SQF, ET, idx_map)

    _populate!(VECOFVAR, SECOND_ORDER_CONE, idx_map)

    _set!(input_prob, MOI.ObjectiveSense(), MOI.get(src, MOI.ObjectiveSense()))

    # TODO NLP...

    return idx_map
end
MOIU.supports_allocate_load(::Optimizer, copy_names::Bool) = !copy_names

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

MOI.get(m::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex) = m._continuous_solution[vi.value]
MOI.get(m::Optimizer, p::MOI.VariablePrimal, vi::Vector{MOI.VariableIndex}) = MOI.get.(m, p, vi)

const EAGO_OPTIMIZER_ATTRIBUTES = Symbol[:relaxed_optimizer, :relaxed_optimizer_kwargs, :upper_optimizer,
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

function MOI.get(m::Optimizer, ::MOI.ObjectiveValue)
    mult = 1.0
    if m._input_problem._optimization_sense == MOI.MAX_SENSE
        mult *= -1.0
    end
    return mult*m._objective_value
end

MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = m._input_problem._variable_num

function MOI.get(m::Optimizer, ::MOI.ObjectiveBound)
    if m._input_problem._optimization_sense === MOI.MAX_SENSE
        bound = -m._global_lower_bound
    else
        bound = m._global_upper_bound
    end
    return bound
end

function MOI.get(m::Optimizer, ::MOI.RelativeGap)
    LBD = m._global_lower_bound
    UBD = m._global_upper_bound
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
