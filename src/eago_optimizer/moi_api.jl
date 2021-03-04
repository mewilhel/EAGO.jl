
#####
#####
##### General MOI utilities required
#####
#####

const INEQ_SETS = Union{LT, GT, ET}
const EAGO_OPTIMIZER_ATTRIBUTES = Symbol[:relaxed_optimizer, :relaxed_optimizer_kwargs, :upper_optimizer,
                                         :enable_optimize_hook, :ext, :ext_type, :_parameters]
const EAGO_MODEL_STRUCT_ATTRIBUTES = Symbol[:_stack, :_log, :_current_node, :_working_problem, :_input_problem]
const EAGO_MODEL_NOT_STRUCT_ATTRIBUTES = setdiff(fieldnames(Optimizer), union(EAGO_OPTIMIZER_ATTRIBUTES,
                                                                              EAGO_MODEL_STRUCT_ATTRIBUTES))

function MOI.empty!(m::Optimizer)

    # create a new empty optimizer and copy fields to m
    new_optimizer = Optimizer()
    for field in union(EAGO_MODEL_STRUCT_ATTRIBUTES, EAGO_MODEL_NOT_STRUCT_ATTRIBUTES)
        setfield!(m, field, getfield(new_optimizer, field))
    end

    return nothing
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

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; copy_names = false)
    return MOI.Utilities.default_copy_to(model, src, copy_names)
end

#####
#####
##### Utilities for checking that JuMP model contains variables used in expression
#####
#####

function check_inbounds!(m::Optimizer, vi::VI)
    if !(1 <= vi.value <= m._input_problem._variable_count)
        error("Invalid variable index $vi. ($(m._input_problem._variable_count) variables in the model.)")
    end
    return nothing
end

check_inbounds!(m::Optimizer, var::SV) = check_inbounds!(m, var.variable)

function check_inbounds!(m::Optimizer, aff::SAF)
    for term in aff.terms
        check_inbounds!(m, term.variable_index)
    end
    return nothing
end

function check_inbounds!(m::Optimizer, quad::SQF)
    for term in quad.affine_terms
        check_inbounds!(m, term.variable_index)
    end
    for term in quad.quadratic_terms
        check_inbounds!(m, term.variable_index_1)
        check_inbounds!(m, term.variable_index_2)
    end
    return nothing
end

function check_inbounds!(m::Optimizer, vov::VECOFVAR)
    for vi in vov.variables
        check_inbounds!(m, vi)
    end
    return
end

#####
#####
##### Set & get attributes of model
#####
#####

function MOI.set(m::Optimizer, ::MOI.Silent, value)

     m._parameters.verbosity = 0
     m._parameters.log_on = false
     return nothing

end

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

function MOI.set(m::Optimizer, ::MOI.TimeLimitSec, value::Nothing)
    m._parameters.time_limit = Inf
    return nothing
end

function MOI.set(m::Optimizer, ::MOI.TimeLimitSec, value::Float64)
    m._parameters.time_limit = value
    return nothing
end

function MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices)
    return [MOI.VariableIndex(i) for i = 1:length(m._input_problem._variable_info)]
end

function MOI.get(m::Optimizer, ::MOI.ObjectiveValue)
    mult = 1.0
    if m._input_problem._optimization_sense === MOI.MAX_SENSE
        mult *= -1.0
    end
    return mult*m._objective_value
end

MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = m._input_problem._variable_count

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

MOI.get(m::Optimizer, ::MOI.SolverName) = "EAGO: Easy Advanced Global Optimization"
MOI.get(m::Optimizer, ::MOI.TerminationStatus) = m._termination_status_code
MOI.get(m::Optimizer, ::MOI.PrimalStatus) = m._result_status_code
MOI.get(m::Optimizer, ::MOI.SolveTime) = m._run_time
MOI.get(m::Optimizer, ::MOI.NodeCount) = m._maximum_node_id
MOI.get(m::Optimizer, ::MOI.ResultCount) = (m._result_status_code === MOI.FEASIBLE_POINT) ? 1 : 0
MOI.get(m::Optimizer, ::MOI.TimeLimitSec) = m.time_limit

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    check_inbounds!(model, vi)
    return model._continuous_solution[vi.value]
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

#####
#####
##### Support, set, and evaluate objective functions
#####
#####
MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{F}) where {F <: Union{SV, SAF, SQF}} = true

function MOI.set(m::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    if nlp_data.has_objective
        m._input_problem._objective_type = NONLINEAR
    end
    m._input_problem._nlp_data = nlp_data
    return nothing
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction{SV}, func::SV)
    check_inbounds!(m, func)
    m._input_problem._objective_sv = func
    m._input_problem._objective_type = SINGLE_VARIABLE
    return nothing
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction{SAF}, func::SAF)
    check_inbounds!(m, func)
    m._input_problem._objective_saf = func
    m._input_problem._objective_type = SCALAR_AFFINE
    return nothing
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction{SQF}, func::SQF)
    check_inbounds!(m, func)
    m._input_problem._objective_sqf = func
    m._input_problem._objective_type = SCALAR_QUADRATIC
    return nothing
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    m._input_problem._optimization_sense = sense
    return nothing
end

##### Supports function and add_constraint for scalar affine functions
MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{S}) where {S <: INEQ_SETS} = true

macro define_addconstraint_linear(function_type, set_type, array_name, count_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            push!(m._input_problem.$(array_name), (func, set))
            m._input_problem._last_constraint_index += 1
            m._input_problem.$(count_name) += 1
            indx = CI{$function_type, $set_type}(m._input_problem._last_constraint_index)
            return indx
        end
    end
end

@define_addconstraint_linear SAF LT _linear_leq_constraints _linear_leq_count
@define_addconstraint_linear SAF GT _linear_geq_constraints _linear_geq_count
@define_addconstraint_linear SAF ET _linear_eq_constraints _linear_eq_count

##### Supports function and add_constraint for scalar quadratic functions
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{S}) where {S <: INEQ_SETS} = true

macro define_addconstraint_quadratic(function_type, set_type, array_name, count_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            push!(m._input_problem.$(array_name), (func, set))
            m._input_problem._last_constraint_index += 1
            m._input_problem.$(count_name) += 1
            indx = CI{$function_type, $set_type}(m._input_problem._last_constraint_index)
            return indx
        end
    end
end

@define_addconstraint_quadratic SQF LT _quadratic_leq_constraints _quadratic_leq_count
@define_addconstraint_quadratic SQF GT _quadratic_geq_constraints _quadratic_geq_count
@define_addconstraint_quadratic SQF ET _quadratic_eq_constraints _quadratic_eq_count

##### Supports function and add_constraint for conic functions
const CONE_SETS = Union{SECOND_ORDER_CONE}
MOI.supports_constraint(::Optimizer, ::Type{VECOFVAR}, ::Type{S}) where {S <: CONE_SETS} = true

macro define_addconstraint_conic(function_type, set_type, array_name, count_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            if length(func.variables) !== set.dimension
                error("Dimension of $(s) does not match number of terms in $(f)")
            end
            check_inbounds!(m, func)
            push!(m._input_problem.$(array_name), (func, set))
            m._input_problem._last_constraint_index += 1
            m._input_problem.$(count_name) += 1
            indx = CI{$function_type, $set_type}(m._input_problem._last_constraint_index)
            return indx
        end
    end
end

@define_addconstraint_quadratic VECOFVAR SECOND_ORDER_CONE _conic_second_order _conic_second_order_count

##### Access variable information from MOI variable index
has_upper_bound(m::Optimizer, vi::MOI.VariableIndex) = m._input_problem._variable_info[vi.value].has_upper_bound
has_lower_bound(m::Optimizer, vi::MOI.VariableIndex) = m._input_problem._variable_info[vi.value].has_lower_bound
is_fixed(m::Optimizer, vi::MOI.VariableIndex) = m._input_problem._variable_info[vi.value].is_fixed
is_integer(m::Optimizer, i::Int64) = is_integer(m._input_problem._variable_info[i])

##### Add unconstrained variables
function MOI.add_variable(m::Optimizer)
    m._input_problem._variable_count += 1
    push!(m._input_problem._variable_info, VariableInfo())
    return VI(m._input_problem._variable_count)
end
MOI.add_variables(m::Optimizer, n::Int) = [MOI.add_variable(m) for i in 1:n]

##### Supports function and add_constraint for single variable functions
MOI.supports_constraint(::Optimizer, ::Type{SV}, ::Type{S}) where {S <: INEQ_SETS} = true

#=
function MOI.add_constraint(m::Optimizer, v::SV, zo::ZO)
    vi = v.variable
    check_inbounds!(m, vi)
    has_upper_bound(m, vi) && error("Upper bound on variable $vi already exists.")
    has_lower_bound(m, vi) && error("Lower bound on variable $vi already exists.")
    is_fixed(m, vi) && error("Variable $vi is fixed. Cannot also set upper bound.")
    m._input_problem._variable_info[vi.value].lower_bound = 0.0
    m._input_problem._variable_info[vi.value].upper_bound = 1.0
    m._input_problem._variable_info[vi.value].has_lower_bound = true
    m._input_problem._variable_info[vi.value].has_upper_bound = true
    m._input_problem._variable_info[vi.value].is_integer = true
    return CI{SV, ZO}(vi.value)
end
=#

function MOI.add_constraint(m::Optimizer, v::SV, lt::LT)
    vi = v.variable
    check_inbounds!(m, vi)
    if isnan(lt.upper)
        error("Invalid upper bound value $(lt.upper).")
    end
    if has_upper_bound(m, vi)
        error("Upper bound on variable $vi already exists.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is fixed. Cannot also set upper bound.")
    end
    m._input_problem._variable_info[vi.value].upper_bound = lt.upper
    m._input_problem._variable_info[vi.value].has_upper_bound = true
    return CI{SV, LT}(vi.value)
end

function MOI.add_constraint(m::Optimizer, v::SV, gt::GT)
    vi = v.variable
    check_inbounds!(m, vi)
    if isnan(gt.lower)
        error("Invalid lower bound value $(gt.lower).")
    end
    if has_lower_bound(m, vi)
        error("Lower bound on variable $vi already exists.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is fixed. Cannot also set lower bound.")
    end
    m._input_problem._variable_info[vi.value].lower_bound = gt.lower
    m._input_problem._variable_info[vi.value].has_lower_bound = true
    return CI{SV, GT}(vi.value)
end

function MOI.add_constraint(m::Optimizer, v::SV, eq::ET)
    vi = v.variable
    check_inbounds!(m, vi)
    if isnan(eq.value)
        error("Invalid fixed value $(gt.lower).")
    end
    if has_lower_bound(m, vi)
        error("Variable $vi has a lower bound. Cannot be fixed.")
    end
    if has_upper_bound(m, vi)
        error("Variable $vi has an upper bound. Cannot be fixed.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is already fixed.")
    end
    m._input_problem._variable_info[vi.value].lower_bound = eq.value
    m._input_problem._variable_info[vi.value].upper_bound = eq.value
    m._input_problem._variable_info[vi.value].has_lower_bound = true
    m._input_problem._variable_info[vi.value].has_upper_bound = true
    m._input_problem._variable_info[vi.value].is_fixed = true
    return CI{SV, ET}(vi.value)
end
